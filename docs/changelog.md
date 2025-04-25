# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Support for advanced matrix inversion methods in `MatrixOperations.cpp`.
- New Python script `advanced_plot.py` for enhanced visualization.

### Fixed
- Memory leak issue in `LyapunovAnalysis.cpp` when handling large datasets.

## [1.1.0] - 2024-04-27

### Added
- Implemented `CalculateDeterminant` function in `MatrixOperations.cpp`.
- Added plotting capability for eigenvalues in `plot_results.py`.

### Changed
- Updated `plot_results.py` to use Matplotlib 3.5 for better performance.

### Deprecated
- Marked `OldLyapunovFunction` as deprecated; use `NewLyapunovFunction` instead.

## [1.0.0] - 2024-03-15

### Added
- Initial release with core functionalities:
  - `MatrixOperations.cpp` and `MatrixOperations.h` for basic matrix computations.
  - `LyapunovAnalysis.cpp` and `LyapunovAnalysis.h` for stability analysis.
  - Python script `plot_results.py` for plotting CSV data.
- Documentation in the `docs/` directory.
