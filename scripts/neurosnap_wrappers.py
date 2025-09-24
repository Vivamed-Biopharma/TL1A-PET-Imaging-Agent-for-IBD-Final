#!/usr/bin/env python3
"""
NeuroSnap Model Wrappers

This module provides high-level wrapper functions for specific NeuroSnap models,
integrating them into the TL1A platform workflow.
"""

import logging
from typing import Dict, Any, List, Optional
from .neurosnap_client import NeuroSnapClient

logger = logging.getLogger(__name__)

class NeuroSnapWrapper:
    """
    Wrapper class for NeuroSnap model integrations.
    """

    def __init__(self, client: Optional[NeuroSnapClient] = None):
        """
        Initialize wrapper with NeuroSnap client.

        Args:
            client: NeuroSnap client instance (optional)
        """
        self.client = client or NeuroSnapClient()

    def predict_admet(self, smiles: str, properties: List[str] = None,
                     max_wait_time: int = 1800) -> Dict[str, Any]:
        """
        Predict ADMET properties using ADMET-AI.

        Args:
            smiles: SMILES string of the molecule
            properties: Properties to predict
            max_wait_time: Maximum wait time in seconds

        Returns:
            Dictionary with ADMET predictions
        """
        logger.info(f"Predicting ADMET for molecule: {smiles[:50]}...")

        job = self.client.submit_admet_job(smiles, properties)
        results = self.client.submit_and_wait("admet_ai",
                                             {"smiles": smiles, "properties": properties},
                                             max_wait_time)

        logger.info("ADMET prediction completed")
        return results

    def predict_toxicity(self, smiles: str, max_wait_time: int = 1800) -> Dict[str, Any]:
        """
        Predict toxicity using eTox.

        Args:
            smiles: SMILES string of the molecule
            max_wait_time: Maximum wait time in seconds

        Returns:
            Dictionary with toxicity predictions
        """
        logger.info(f"Predicting toxicity for molecule: {smiles[:50]}...")

        job = self.client.submit_etox_job(smiles)
        results = self.client.submit_and_wait("etox", {"smiles": smiles}, max_wait_time)

        logger.info("Toxicity prediction completed")
        return results

    def predict_aggregation(self, sequence: str, max_wait_time: int = 1800) -> Dict[str, Any]:
        """
        Predict aggregation hotspots using Aggrescan3D.

        Args:
            sequence: Protein sequence
            max_wait_time: Maximum wait time in seconds

        Returns:
            Dictionary with aggregation predictions
        """
        logger.info(f"Predicting aggregation for sequence: {sequence[:20]}...")

        job = self.client.submit_aggrescan_job(sequence)
        results = self.client.submit_and_wait("aggrescan3d", {"sequence": sequence}, max_wait_time)

        logger.info("Aggregation prediction completed")
        return results

    def predict_thermostability(self, sequence: str, temperature: float = 25.0,
                              max_wait_time: int = 1800) -> Dict[str, Any]:
        """
        Predict thermostability using ThermoMPNN.

        Args:
            sequence: Protein sequence
            temperature: Temperature in Celsius
            max_wait_time: Maximum wait time in seconds

        Returns:
            Dictionary with thermostability predictions
        """
        logger.info(f"Predicting thermostability for sequence: {sequence[:20]}...")

        job = self.client.submit_thermompn_job(sequence, temperature)
        results = self.client.submit_and_wait("thermompn",
                                             {"sequence": sequence, "temperature": temperature},
                                             max_wait_time)

        logger.info("Thermostability prediction completed")
        return results

    def predict_structure(self, sequences: List[str], max_wait_time: int = 3600) -> Dict[str, Any]:
        """
        Predict protein structure using Boltz-2.

        Args:
            sequences: List of protein sequences
            max_wait_time: Maximum wait time in seconds

        Returns:
            Dictionary with structure predictions
        """
        logger.info(f"Predicting structure for {len(sequences)} sequence(s)")

        job = self.client.submit_boltz_job(sequences)
        results = self.client.submit_and_wait("boltz2", {"sequences": sequences}, max_wait_time)

        logger.info("Structure prediction completed")
        return results

    def predict_stability_change(self, sequence: str, mutations: List[str],
                               max_wait_time: int = 1800) -> Dict[str, Any]:
        """
        Predict stability changes using StaB-ddG.

        Args:
            sequence: Wild-type protein sequence
            mutations: List of mutations (e.g., ["A1V", "K2M"])
            max_wait_time: Maximum wait time in seconds

        Returns:
            Dictionary with stability predictions
        """
        logger.info(f"Predicting stability changes for {len(mutations)} mutation(s)")

        job = self.client.submit_stab_ddg_job(sequence, mutations)
        results = self.client.submit_and_wait("stab_ddg",
                                             {"sequence": sequence, "mutations": mutations},
                                             max_wait_time)

        logger.info("Stability prediction completed")
        return results

    def predict_immunogenicity(self, sequence: str, max_wait_time: int = 1800) -> Dict[str, Any]:
        """
        Predict immunogenicity using DeepImmuno.

        Args:
            sequence: Protein sequence
            max_wait_time: Maximum wait time in seconds

        Returns:
            Dictionary with immunogenicity predictions
        """
        logger.info(f"Predicting immunogenicity for sequence: {sequence[:20]}...")

        job = self.client.submit_deepimmuno_job(sequence)
        results = self.client.submit_and_wait("deepimmuno", {"sequence": sequence}, max_wait_time)

        logger.info("Immunogenicity prediction completed")
        return results

# Convenience functions for direct use
def predict_admet(smiles: str, properties: List[str] = None) -> Dict[str, Any]:
    """Convenience function for ADMET prediction."""
    wrapper = NeuroSnapWrapper()
    return wrapper.predict_admet(smiles, properties)

def predict_toxicity(smiles: str) -> Dict[str, Any]:
    """Convenience function for toxicity prediction."""
    wrapper = NeuroSnapWrapper()
    return wrapper.predict_toxicity(smiles)

def predict_aggregation(sequence: str) -> Dict[str, Any]:
    """Convenience function for aggregation prediction."""
    wrapper = NeuroSnapWrapper()
    return wrapper.predict_aggregation(sequence)

def predict_thermostability(sequence: str, temperature: float = 25.0) -> Dict[str, Any]:
    """Convenience function for thermostability prediction."""
    wrapper = NeuroSnapWrapper()
    return wrapper.predict_thermostability(sequence, temperature)

def predict_structure(sequences: List[str]) -> Dict[str, Any]:
    """Convenience function for structure prediction."""
    wrapper = NeuroSnapWrapper()
    return wrapper.predict_structure(sequences)

def predict_stability_change(sequence: str, mutations: List[str]) -> Dict[str, Any]:
    """Convenience function for stability change prediction."""
    wrapper = NeuroSnapWrapper()
    return wrapper.predict_stability_change(sequence, mutations)

def predict_immunogenicity(sequence: str) -> Dict[str, Any]:
    """Convenience function for immunogenicity prediction."""
    wrapper = NeuroSnapWrapper()
    return wrapper.predict_immunogenicity(sequence)