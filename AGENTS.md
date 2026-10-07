# Agent Guidelines for dist-m4ri

## Interaction Rules
- **No unsolicited code edits on questions**: When the user asks a question, requests an explanation, discusses
  architecture, or asks for advice, provide a clear, concise answer first. Do NOT edit code or modify repository files
  without explicit instructions from the user to proceed.

## Code and Formatting Standards
- **Line Length**: Wrap all lines in source code (C, Python), documentation (Markdown), workflow files (YAML), and test
  scripts to a maximum of 120 columns.
- **Documentation Integrity**: Maintain documentation integrity. Keep comments, docstrings, and help messages
  up-to-date and consistent between the C binary and the Python library.

## Version Management
- **Synchronized Versioning**: When the version is incremented, ensure it is updated consistently across all four
  project metadata files:
  1. `CITATION.cff` (`version: vX.Y.Z`)
  2. `Doxyfile` (`PROJECT_NUMBER = X.Y.Z`)
  3. `src/util_io.h` (`#define DIST_M4RI_VERSION "X.Y.Z"`)
  4. `dist_m4ri.py` (`__version__ = "X.Y.Z"`)
- **Compatibility**: Ensure binary CLI arguments and Python CLI arguments remain backward-compatible, and cache file
  versioning (`"__version__"`) is preserved.
