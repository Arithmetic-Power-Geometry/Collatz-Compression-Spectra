#!/usr/bin/env bash
set -euo pipefail
python -m pip install -r requirements.txt
python -m pytest -q
python run_analysis.py --out analysis_output --n-raw 200000 --n-spectrum 200000 --max-m 128 --n-validation 10000
