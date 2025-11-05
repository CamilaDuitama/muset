# Changelog

All notable changes to muset will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.6.0] - 2025-11-05

### Added
- `--abundance-metric` option to choose between mean (default) or median for unitig abundance calculation
- `--output-format` option to output matrices as txt (default) or compressed tsv format
- `-e/--generate-maximal-unitigs-links` flag to generate GGCAT CDBG with edges in BCALM2 format
- `generate_fof` script to automatically create file-of-files from a directory with paired-end support
- Unit tests for aggregator interfaces and statistical metrics
- Separated interfaces for statistical metrics computation and matrix writing for better extensibility

### Changed
- Skip k-mer matrix filtering when all filter parameters are zero (f=0, F=0, n=0, N=0) for improved performance
- Refactored unitig metric computation and matrix writing with cleaner separation of concerns

### Performance
- Filtering step now directly aggregates matrices when no actual filtering is needed, reducing I/O overhead

## [0.5.1] - 2024-XX-XX
- Previous release

[Unreleased]: https://github.com/CamilaDuitama/muset/compare/v0.6.0...HEAD
[0.6.0]: https://github.com/CamilaDuitama/muset/compare/v0.5.1...v0.6.0
[0.5.1]: https://github.com/CamilaDuitama/muset/releases/tag/v0.5.1
