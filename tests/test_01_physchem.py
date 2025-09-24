"""
Unit tests for Experiment 1: Physicochemical Profiling
"""

import pytest
import pandas as pd
from scripts import inputs
from scripts.logging_config import setup_logging

class TestPhyschem:
    """Test physicochemical profiling functionality."""

    @pytest.fixture
    def setup_logging(self):
        """Setup logging for tests."""
        return setup_logging(log_level=logging.WARNING)

    def test_calculate_properties(self):
        """Test property calculation."""
        from scripts._01_physchem import calculate_properties

        test_molecules = {
            "Test": "CCO"  # Ethanol
        }

        df = calculate_properties(test_molecules)

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1
        assert "MW" in df.columns
        assert "LogP" in df.columns
        assert df["MW"].values[0] > 40  # Ethanol MW ~46

    def test_input_validation(self):
        """Test input validation."""
        from scripts._01_physchem import validate_inputs

        # Should not raise exception with valid inputs
        validate_inputs()

    def test_output_structure(self):
        """Test that output has expected columns."""
        from scripts._01_physchem import calculate_properties

        test_molecules = {"Test": "CCO"}
        df = calculate_properties(test_molecules)

        expected_columns = ["Name", "SMILES", "MW", "LogP", "HBD", "HBA", "Rotatable_Bonds", "TPSA"]
        for col in expected_columns:
            assert col in df.columns