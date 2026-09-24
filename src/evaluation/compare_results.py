#!/usr/bin/env python
"""Backward-compatible entry point for result figures and comparisons.

The implementation lives in generate_results.py so the paper comparison and
all figures/tables use exactly the same best-weighted-F1 selection rule.
"""
from src.evaluation.generate_results import main


if __name__ == "__main__":
    main()
