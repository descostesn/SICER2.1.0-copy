# Changelog
All notable changes to this project will be documented in this file.

## [2.1.1] - 2026-10-19
### Fixed
- Paired-end island read counts now count each fragment at its midpoint, the position the islands are built from, instead of start + fragment_size/2.
- sicer_df, recognicer and recognicer_df no longer stop with "'Namespace' object has no attribute 'input_type'": a run without input_type is single-end.
- sicer --significant_reads on single-end data no longer stops with an IndexError in the filtered graph step: the reads kept in islands keep their record type.
- sicer_df and recognicer_df no longer stop on "module 'scipy' has no attribute 'array'" with recent SciPy: the FDR step uses numpy.array.

## [2.1] - 2026-02-05
### Fixed
- Fixed generating bed file from bam file duplicate issue.

## [2.1] - 2026-01-03
### Added
- SICER2.1 now accepts BAM file as input data.

## [2.1] - 2025-12-08
### Added
- SICER2.1 now accepts paired-end data.

