#!/usr/bin/env python3

# Crashes of sicer_df, recognicer, recognicer_df and single-end --significant_reads.
# Run from anywhere: python test/test_scripts/test_df_recognicer_significant_reads.py
# It imports the sicer package of this checkout, not the installed one.
# Each part fails on upstream v2.1.1 and passes with the fix; the script exits 1 if any part fails.

import contextlib
import inspect
import io
import math
import multiprocessing as mp
import os
import sys
import tempfile
import traceback
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from sicer.src import compare_two_libraries_on_islands
from sicer.src import filter_raw_tags_by_islands
from sicer.src import remove_redundant_reads
from sicer.src import run_make_graph_file_by_chrom

CHROM = 'chr1'
CHROM_LENGTH = 100000
FRAGMENT_SIZE = 150
READ_DTYPE = np.dtype([('chrom', 'U6'), ('start', np.int32), ('end', np.int32), ('name', 'U20'),
                       ('score', np.int32), ('strand', 'U1')])


def in_tmp(function):
    # SICER writes its .npy files in the working directory.
    def wrapped():
        old_dir = os.getcwd()
        with tempfile.TemporaryDirectory() as tmp:
            os.chdir(tmp)
            try:
                return function(tmp)
            finally:
                os.chdir(old_dir)
    wrapped.__name__ = function.__name__
    return wrapped


@in_tmp
def test_no_input_type(tmp):
    # Fix A: sicer_df, recognicer and recognicer_df never set args.input_type.
    bed = os.path.join(tmp, 'treatment.bed')
    with open(bed, 'w') as f:
        f.write('chr1\t100\t136\tr1\t0\t+\n')
        f.write('chr1\t100\t136\tr2\t0\t+\n')  # redundant with r1, dropped at threshold 1
        f.write('chr1\t500\t536\tr3\t0\t-\n')
        f.write('chr1\t900\t936\tr4\t0\t+\n')
        f.write('chr2\t100\t136\tr5\t0\t-\n')
    args = SimpleNamespace(species='mm10', redundancy_threshold=1)
    with mp.Pool(2) as pool, contextlib.redirect_stdout(io.StringIO()):  # its per-chromosome table
        total = remove_redundant_reads.main(args, bed, pool)
    assert total == 4, 'expected 4 retained reads, got ' + str(total)


def filter_reads(file_name, chrom):
    # The single-end call, with or without the input_type argument of the PE island-count fix.
    if 'input_type' in inspect.signature(filter_raw_tags_by_islands.filter_tags_by_islands).parameters:
        return filter_raw_tags_by_islands.filter_tags_by_islands(file_name, FRAGMENT_SIZE, 'SE', chrom)
    return filter_raw_tags_by_islands.filter_tags_by_islands(file_name, FRAGMENT_SIZE, chrom)


def check_filtered(tmp, islands, expected_names, expected_tags):
    # Fix B: the kept reads keep their record type, so get_bed_coords' SE guard can read their fields.
    reads = np.array([(CHROM, 1000, 1036, 'r1', 0, '+'),    # tag 1075, in the island
                      (CHROM, 1500, 1536, 'r2', 0, '-'),    # tag 1460, in the island
                      (CHROM, 3000, 3036, 'r3', 0, '+')],   # tag 3075, outside
                     dtype=READ_DTYPE)
    np.save(os.path.join(tmp, 'treatment_' + CHROM + '.npy'), reads)
    np.save(os.path.join(tmp, 'treatment_' + CHROM + '_island_summary.npy'), np.array(islands, dtype=object))
    filter_reads(os.path.join(tmp, 'treatment'), CHROM)
    filtered = np.load(os.path.join(tmp, 'treatment_' + CHROM + '_filtered.npy'), allow_pickle=True)
    assert filtered.dtype == READ_DTYPE, 'filtered reads have dtype ' + str(filtered.dtype)
    assert list(filtered['name']) == expected_names, 'kept ' + str(list(filtered['name']))
    tags, _ = run_make_graph_file_by_chrom.get_bed_coords(filtered, CHROM_LENGTH, FRAGMENT_SIZE, CHROM, False, 'SE')
    assert list(tags) == expected_tags, 'tags ' + str(list(tags))


@in_tmp
def test_filtered_reads_with_islands(tmp):
    check_filtered(tmp, [(CHROM, 1000, 1599, 5, 0, 1e-5, 10.0, 1e-4)], ['r1', 'r2'], [1075, 1460])


@in_tmp
def test_filtered_reads_without_islands(tmp):
    check_filtered(tmp, [], [], [])


def test_fdr():
    # Fix C: scipy.array is gone from recent SciPy.
    observed = compare_two_libraries_on_islands.fdr([0.01, 0.04, 0.03])
    expected = [0.03, 0.04, 0.045]
    assert all(math.isclose(o, e) for o, e in zip(observed, expected)) and len(observed) == 3, str(observed)


if __name__ == '__main__':
    failed = 0
    for test in (test_no_input_type, test_filtered_reads_with_islands, test_filtered_reads_without_islands, test_fdr):
        try:
            test()
            print(test.__name__ + ': OK')
        except Exception:
            failed += 1
            print(test.__name__ + ': FAILED')
            traceback.print_exc(limit=1)
    sys.exit(1 if failed else 0)
