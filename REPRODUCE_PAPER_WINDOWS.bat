@echo off
python -m pip install -r requirements.txt || exit /b 1
python -m pytest -q || exit /b 1
python run_analysis.py --out analysis_output --n-raw 200000 --n-spectrum 200000 --max-m 128 --n-validation 10000 || exit /b 1
