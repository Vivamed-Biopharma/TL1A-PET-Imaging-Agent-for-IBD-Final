#!/usr/bin/env bash
set -euo pipefail
echo "[1/1] Regenerating report and artifacts..."
python report.py
echo "Done. See REPORT.md and root CSV/PNG outputs."
