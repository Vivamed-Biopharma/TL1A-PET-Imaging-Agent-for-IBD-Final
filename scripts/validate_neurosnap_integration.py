#!/usr/bin/env python3
"""
Validation script for NeuroSnap API integration.

This script tests all NeuroSnap integrations with real API calls to ensure
they work correctly with actual data.
"""

import logging
import sys
from pathlib import Path
import scripts.inputs as inputs
from scripts.neurosnap_wrappers import NeuroSnapWrapper

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_admet_integration():
    """Test ADMET prediction integration."""
    logger.info("Testing ADMET prediction integration...")

    wrapper = NeuroSnapWrapper()
    try:
        results = wrapper.predict_admet(inputs.LINKER_CHELATOR_SMILES, max_wait_time=600)
        logger.info(f"ADMET prediction successful: {len(results)} properties predicted")
        return True, results
    except Exception as e:
        logger.error(f"ADMET prediction failed: {str(e)}")
        return False, str(e)

def test_toxicity_integration():
    """Test toxicity prediction integration."""
    logger.info("Testing toxicity prediction integration...")

    wrapper = NeuroSnapWrapper()
    try:
        results = wrapper.predict_toxicity(inputs.LINKER_CHELATOR_SMILES, max_wait_time=600)
        logger.info(f"Toxicity prediction successful: {results.get('toxic', 'unknown')} toxicity predicted")
        return True, results
    except Exception as e:
        logger.error(f"Toxicity prediction failed: {str(e)}")
        return False, str(e)

def test_aggregation_integration():
    """Test aggregation prediction integration."""
    logger.info("Testing aggregation prediction integration...")

    wrapper = NeuroSnapWrapper()
    sequence = inputs.fab_sequences["Fab06_VH"]

    try:
        results = wrapper.predict_aggregation(sequence, max_wait_time=600)
        score = results.get('aggregation_score', 'unknown')
        logger.info(f"Aggregation prediction successful: score {score}")
        return True, results
    except Exception as e:
        logger.error(f"Aggregation prediction failed: {str(e)}")
        return False, str(e)

def test_thermostability_integration():
    """Test thermostability prediction integration."""
    logger.info("Testing thermostability prediction integration...")

    wrapper = NeuroSnapWrapper()
    sequence = inputs.fab_sequences["Fab06_VH"]

    try:
        results = wrapper.predict_thermostability(sequence, max_wait_time=600)
        tm = results.get('melting_temperature', 'unknown')
        logger.info(f"Thermostability prediction successful: Tm {tm}°C")
        return True, results
    except Exception as e:
        logger.error(f"Thermostability prediction failed: {str(e)}")
        return False, str(e)

def test_structure_integration():
    """Test structure prediction integration."""
    logger.info("Testing structure prediction integration...")

    wrapper = NeuroSnapWrapper()
    sequences = [inputs.fab_sequences["Fab06_VH"], inputs.fab_sequences["Fab06_VL"]]

    try:
        results = wrapper.predict_structure(sequences, max_wait_time=1200)  # Longer timeout for structure prediction
        logger.info(f"Structure prediction successful: {len(results.get('structures', []))} structures predicted")
        return True, results
    except Exception as e:
        logger.error(f"Structure prediction failed: {str(e)}")
        return False, str(e)

def test_stability_change_integration():
    """Test stability change prediction integration."""
    logger.info("Testing stability change prediction integration...")

    wrapper = NeuroSnapWrapper()
    sequence = inputs.fab_sequences["Fab06_VH"]
    mutations = ["A1V", "K2M"]  # Example mutations

    try:
        results = wrapper.predict_stability_change(sequence, mutations, max_wait_time=600)
        logger.info(f"Stability change prediction successful: {len(mutations)} mutations analyzed")
        return True, results
    except Exception as e:
        logger.error(f"Stability change prediction failed: {str(e)}")
        return False, str(e)

def test_immunogenicity_integration():
    """Test immunogenicity prediction integration."""
    logger.info("Testing immunogenicity prediction integration...")

    wrapper = NeuroSnapWrapper()
    sequence = inputs.fab_sequences["Fab06_VH"]

    try:
        results = wrapper.predict_immunogenicity(sequence, max_wait_time=600)
        score = results.get('immunogenicity_score', 'unknown')
        logger.info(f"Immunogenicity prediction successful: score {score}")
        return True, results
    except Exception as e:
        logger.error(f"Immunogenicity prediction failed: {str(e)}")
        return False, str(e)

def main():
    """Run all NeuroSnap integration validations."""
    logger.info("Starting NeuroSnap API integration validation...")

    # Test data validation
    if not hasattr(inputs, 'LINKER_CHELATOR_SMILES'):
        logger.error("LINKER_CHELATOR_SMILES not found in inputs")
        sys.exit(1)

    if not hasattr(inputs, 'fab_sequences'):
        logger.error("fab_sequences not found in inputs")
        sys.exit(1)

    # Run all integration tests
    tests = [
        ("ADMET Prediction", test_admet_integration),
        ("Toxicity Prediction", test_toxicity_integration),
        ("Aggregation Prediction", test_aggregation_integration),
        ("Thermostability Prediction", test_thermostability_integration),
        ("Structure Prediction", test_structure_integration),
        ("Stability Change Prediction", test_stability_change_integration),
        ("Immunogenicity Prediction", test_immunogenicity_integration),
    ]

    results = {}
    passed = 0
    failed = 0

    for test_name, test_func in tests:
        logger.info(f"\n{'='*50}")
        logger.info(f"Running {test_name}")
        logger.info(f"{'='*50}")

        try:
            success, result = test_func()
            results[test_name] = {"success": success, "result": result}

            if success:
                logger.info(f"✅ {test_name} PASSED")
                passed += 1
            else:
                logger.error(f"❌ {test_name} FAILED: {result}")
                failed += 1

        except Exception as e:
            logger.error(f"❌ {test_name} CRASHED: {str(e)}")
            results[test_name] = {"success": False, "result": str(e)}
            failed += 1

    # Summary
    logger.info(f"\n{'='*60}")
    logger.info("VALIDATION SUMMARY")
    logger.info(f"{'='*60}")
    logger.info(f"Total tests: {len(tests)}")
    logger.info(f"Passed: {passed}")
    logger.info(f"Failed: {failed}")

    if failed == 0:
        logger.info("🎉 All NeuroSnap integrations validated successfully!")
        logger.info("The platform is ready for production use with real API calls.")
        sys.exit(0)
    else:
        logger.error(f"❌ {failed} integration(s) failed validation.")
        logger.error("Please check the API key, network connectivity, and API service status.")
        sys.exit(1)

if __name__ == "__main__":
    main()