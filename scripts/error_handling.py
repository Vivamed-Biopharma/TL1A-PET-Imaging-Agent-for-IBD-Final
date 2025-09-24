#!/usr/bin/env python3
"""
Common error handling utilities for the TL1A platform.
"""

import logging
import functools
import sys
from typing import Callable, Any

logger = logging.getLogger(__name__)

class TL1AError(Exception):
    """Base exception for TL1A platform errors."""
    pass

class ValidationError(TL1AError):
    """Raised when input validation fails."""
    pass

class APIError(TL1AError):
    """Raised when API calls fail."""
    pass

class FileNotFoundError(TL1AError):
    """Raised when required files are missing."""
    pass

def handle_errors(func: Callable) -> Callable:
    """
    Decorator to handle common errors in experiment functions.

    Args:
        func: Function to decorate

    Returns:
        Decorated function
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except ValidationError as e:
            logger.error(f"Validation error in {func.__name__}: {str(e)}")
            raise
        except APIError as e:
            logger.error(f"API error in {func.__name__}: {str(e)}")
            raise
        except FileNotFoundError as e:
            logger.error(f"File not found in {func.__name__}: {str(e)}")
            raise
        except Exception as e:
            logger.error(f"Unexpected error in {func.__name__}: {str(e)}")
            logger.error(f"Error type: {type(e).__name__}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            raise TL1AError(f"Unexpected error: {str(e)}") from e

    return wrapper

def validate_smiles(smiles: str, name: str = "molecule") -> None:
    """
    Validate SMILES string.

    Args:
        smiles: SMILES string to validate
        name: Name for error messages

    Raises:
        ValidationError: If SMILES is invalid
    """
    if not smiles or not isinstance(smiles, str):
        raise ValidationError(f"Invalid SMILES for {name}: must be non-empty string")

    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValidationError(f"Invalid SMILES syntax for {name}: {smiles}")

def validate_sequence(sequence: str, name: str = "sequence", allow_gaps: bool = False) -> None:
    """
    Validate amino acid sequence.

    Args:
        sequence: Sequence to validate
        name: Name for error messages
        allow_gaps: Whether to allow gaps (-)

    Raises:
        ValidationError: If sequence is invalid
    """
    if not sequence or not isinstance(sequence, str):
        raise ValidationError(f"Invalid sequence for {name}: must be non-empty string")

    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    if allow_gaps:
        valid_aa.add('-')

    seq_upper = sequence.upper()
    invalid_chars = set(seq_upper) - valid_aa

    if invalid_chars:
        raise ValidationError(f"Invalid amino acids in {name}: {invalid_chars}")

def validate_file_exists(filepath: str) -> None:
    """
    Validate that a file exists.

    Args:
        filepath: Path to file

    Raises:
        FileNotFoundError: If file doesn't exist
    """
    import os
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Required file not found: {filepath}")

def safe_float_convert(value: Any, default: float = 0.0) -> float:
    """
    Safely convert value to float.

    Args:
        value: Value to convert
        default: Default value if conversion fails

    Returns:
        Float value or default
    """
    try:
        return float(value)
    except (ValueError, TypeError):
        logger.warning(f"Could not convert {value} to float, using default {default}")
        return default

def safe_int_convert(value: Any, default: int = 0) -> int:
    """
    Safely convert value to int.

    Args:
        value: Value to convert
        default: Default value if conversion fails

    Returns:
        Int value or default
    """
    try:
        return int(value)
    except (ValueError, TypeError):
        logger.warning(f"Could not convert {value} to int, using default {default}")
        return default