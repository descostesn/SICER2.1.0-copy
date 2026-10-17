#!/usr/bin/env python3

# Paired-end island read counts at the fragment midpoint.
# Run from anywhere: python test/test_scripts/test_pe_island_count.py
# It imports the sicer package of this checkout, not the installed one.

import os
import sys
import tempfile
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from sicer.lib import associate_tags_with_regions
from sicer.src import associate_tags_with_chip_and_control_w_fc_q
from sicer.src import filter_raw_tags_by_islands

FRAGMENT_SIZE = 150
CHROM = 'chr1'
READ_DTYPE = np.dtype([('chrom', 'U6'), ('start', np.int32), ('end', np.int32), ('name', 'U20'),
                       ('score', np.int32), ('strand', 'U1')])

# Two islands; both ends are inclusive (find_readcount_on_islands).
ISLANDS = [(CHROM, 1000, 1599, 10.0), (CHROM, 2000, 2399, 10.0)]

# PE fragments, one BED line each, all on '+' (as convert_bam_to_bed writes them).
# (start, end): midpoint -> island | old position start + 75 -> island
FRAGMENTS = [
    (950, 1100),   # 150 bp: 1025 -> 1    | 1025 -> 1
    (1449, 1749),  # 300 bp: 1599 -> 1    | 1524 -> 1   (midpoint on the island end)
    (1450, 1750),  # 300 bp: 1600 -> none | 1525 -> 1
    (1900, 2200),  # 300 bp: 2050 -> 2    | 1975 -> none
    (700, 1300),   # 600 bp: 1000 -> 1    | 775 -> none (midpoint on the island start)
    (1500, 2100),  # 600 bp: 1800 -> none | 1575 -> 1
    (1900, 2500),  # 600 bp: 2200 -> 2    | 1975 -> none
]
CONTROL = [FRAGMENTS[2], FRAGMENTS[3], FRAGMENTS[4]]

# Expected (ChIP counts per island, control counts per island, fragments kept by the filter)
EXPECTED = {
    'PE': ([3, 2], [1, 1], [0, 1, 3, 4, 6]),
    'SE': ([4, 0], [1, 0], [0, 1, 2, 5]),
}


def save_reads(path, fragments):
    reads = np.array([(CHROM, s, e, 'f' + str(i), 0, '+') for i, (s, e) in enumerate(fragments)], dtype=READ_DTYPE)
    np.save(path, reads)


def run(input_type):
    old_dir = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        # associate_tag_count_to_regions writes <chrom>_pvalue.npy in the working directory.
        os.chdir(tmp)
        try:
            treatment = os.path.join(tmp, 'treatment.bed')
            control = os.path.join(tmp, 'control.bed')
            np.save(os.path.join(tmp, 'treatment_' + CHROM + '_graph.npy'), np.array(ISLANDS, dtype=object))
            save_reads(os.path.join(tmp, 'treatment_' + CHROM + '.npy'), FRAGMENTS)
            save_reads(os.path.join(tmp, 'control_' + CHROM + '.npy'), CONTROL)

            args = SimpleNamespace(treatment_file=treatment, control_file=control,
                                   fragment_size=FRAGMENT_SIZE, input_type=input_type)
            associate_tags_with_chip_and_control_w_fc_q.associate_tag_count_to_regions(
                args, 1.0, len(CONTROL), 1e6, CHROM)
            summary = np.load(os.path.join(tmp, 'treatment_' + CHROM + '_island_summary.npy'), allow_pickle=True)
            chip_counts = [int(island[3]) for island in summary]
            control_counts = [int(island[4]) for island in summary]

            filter_raw_tags_by_islands.filter_tags_by_islands(
                os.path.join(tmp, 'treatment'), FRAGMENT_SIZE, input_type, CHROM)
            filtered = np.load(os.path.join(tmp, 'treatment_' + CHROM + '_filtered.npy'), allow_pickle=True)
            kept = sorted(int(read[3][1:]) for read in filtered)
        finally:
            os.chdir(old_dir)
    return chip_counts, control_counts, kept


def test_tag_position():
    plus = ('chr1', 100, 400, 'f', 0, '+')
    minus = ('chr1', 100, 400, 'f', 0, '-')
    for tag_position in (associate_tags_with_regions.tag_position, filter_raw_tags_by_islands.tag_position):
        assert tag_position(plus, FRAGMENT_SIZE, 'PE') == 250
        assert tag_position(minus, FRAGMENT_SIZE, 'PE') == 250
        assert tag_position(plus, FRAGMENT_SIZE) == 175
        assert tag_position(minus, FRAGMENT_SIZE) == 324
        assert tag_position(plus, FRAGMENT_SIZE, 'SE') == 175


def test_island_counts():
    for input_type, expected in EXPECTED.items():
        observed = run(input_type)
        assert observed == expected, input_type + ': expected ' + str(expected) + ', got ' + str(observed)


if __name__ == '__main__':
    test_tag_position()
    print('tag_position: OK')
    test_island_counts()
    print('island counts and filter, PE at the midpoint and SE unchanged: OK')
