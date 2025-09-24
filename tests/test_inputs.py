"""
Unit tests for inputs.py
"""

import pytest
import scripts.inputs as inputs
import scripts.error_handling as eh

class TestInputs:
    """Test input data validation."""

    def test_smiles_validity(self):
        """Test that SMILES strings are valid."""
        from rdkit import Chem

        smiles_list = [
            inputs.NOTA_CHELATOR_SMILES,
            inputs.LINKER_CHELATOR_SMILES
        ]

        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            assert mol is not None, f"Invalid SMILES: {smiles}"
            assert mol.GetNumAtoms() > 0, f"Empty molecule: {smiles}"

    def test_fab_sequences_present(self):
        """Test that Fab sequences are loaded."""
        assert hasattr(inputs, 'fab_sequences')
        assert len(inputs.fab_sequences) == 4  # Two Fabs with VH/VL each
        assert 'Fab06_VH' in inputs.fab_sequences
        assert 'Fab06_VL' in inputs.fab_sequences
        assert 'Fab11_VH' in inputs.fab_sequences
        assert 'Fab11_VL' in inputs.fab_sequences

    def test_sequence_length(self):
        """Test that sequences have reasonable lengths."""
        for name, seq in inputs.fab_sequences.items():
            assert len(seq) > 100, f"Sequence {name} too short: {len(seq)}"
            assert len(seq) < 150, f"Sequence {name} too long: {len(seq)}"
            assert len(seq) % 3 == 0, f"Sequence {name} length not divisible by 3: {len(seq)}"

    def test_amino_acids_valid(self):
        """Test that sequences contain only valid amino acids."""
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        for name, seq in inputs.fab_sequences.items():
            seq_aa = set(seq.upper())
            invalid = seq_aa - valid_aa
            assert len(invalid) == 0, f"Invalid amino acids in {name}: {invalid}"

    def test_validation_functions(self):
        """Test validation functions work."""
        # Valid SMILES
        eh.validate_smiles("CCO", "ethanol")

        # Valid sequence
        eh.validate_sequence("ACDEFGHIKLMNPQRSTVWY", "test_seq")

        # Invalid SMILES
        with pytest.raises(eh.ValidationError):
            eh.validate_smiles("invalid", "test")

        # Invalid sequence
        with pytest.raises(eh.ValidationError):
            eh.validate_sequence("XYZ", "test")

    def test_biotin_sequences(self):
        """Test Fc and FcRn sequences."""
        assert hasattr(inputs, 'FC_SEQUENCE')
        assert hasattr(inputs, 'FCRN_ALPHA_CHAIN_SEQUENCE')
        assert hasattr(inputs, 'FCRN_BETA2M_SEQUENCE')

        eh.validate_sequence(inputs.FC_SEQUENCE, "FC")
        eh.validate_sequence(inputs.FCRN_ALPHA_CHAIN_SEQUENCE, "FcRn_Alpha")
        eh.validate_sequence(inputs.FCRN_BETA2M_SEQUENCE, "FcRn_Beta2M")

    def test_safe_conversions(self):
        """Test safe conversion functions."""
        assert eh.safe_float_convert("1.5") == 1.5
        assert eh.safe_float_convert("invalid") == 0.0
        assert eh.safe_int_convert("42") == 42
        assert eh.safe_int_convert("invalid") == 0