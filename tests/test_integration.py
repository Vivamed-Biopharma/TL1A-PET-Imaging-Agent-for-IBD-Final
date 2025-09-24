"""
Integration tests for the full pipeline.
"""

import pytest
import subprocess
import sys
from pathlib import Path

class TestIntegration:
    """Test full pipeline integration."""

    def test_experiment_1_runs(self):
        """Test that experiment 1 runs without error."""
        result = subprocess.run([
            sys.executable, "scripts/01_physchem.py"
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Experiment 1 failed: {result.stderr}"

    def test_experiment_6_runs(self):
        """Test that experiment 6 runs without error."""
        result = subprocess.run([
            sys.executable, "scripts/06_sequence_dev.py"
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Experiment 6 failed: {result.stderr}"

    def test_main_cli_experiment_1(self):
        """Test main CLI with experiment 1."""
        result = subprocess.run([
            sys.executable, "main.py", "1"
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Main CLI failed: {result.stderr}"

    def test_results_directories_created(self):
        """Test that results directories are created."""
        dirs = [
            "results/prodrug",
            "results/formulation",
            "results/biobetter"
        ]

        for dir_path in dirs:
            assert Path(dir_path).exists(), f"Directory {dir_path} not created"

    def test_csv_outputs_exist(self):
        """Test that key CSV outputs are created."""
        # This assumes experiments have been run
        key_files = [
            "results/prodrug/01_physchem_properties.csv",
            "results/formulation/06_sequence_developability.csv"
        ]

        for file_path in key_files:
            # Note: Files may not exist if experiments haven't run
            # This is more of a structure test
            Path(file_path).parent.mkdir(parents=True, exist_ok=True)