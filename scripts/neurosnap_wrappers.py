#!/usr/bin/env python3
"""
NeuroSnap Model Wrappers

High-level wrappers that prepare inputs for NeuroSnap services, submit jobs via
NeuroSnapClient, reuse existing jobs by note hash, and parse downloaded results
into convenient Python dictionaries for downstream scripts.

These wrappers preserve the existing public API used by experiment scripts
(returning dictionaries), while under the hood using the real NeuroSnap API.
"""

import json
import os
import hashlib
import logging
from typing import Dict, Any, List, Optional

try:
    from .neurosnap_client import NeuroSnapClient
except Exception:
    # Fallback when executed as a script without package context
    from neurosnap_client import NeuroSnapClient

logger = logging.getLogger(__name__)


def _get_input_hash(data: Any) -> str:
    """Create a SHA256 hash of the note data for idempotent job reuse."""
    try:
        serialized = json.dumps(data, sort_keys=True, default=str)
    except TypeError:
        serialized = str(data)
    return hashlib.sha256(serialized.encode()).hexdigest()


def _parse_downloaded_files(file_paths: List[str]) -> Dict[str, Any]:
    """Best-effort parse of downloaded result files into a dictionary.
    Prefers JSON; otherwise returns minimal metadata.
    """
    results: Dict[str, Any] = {}
    for path in file_paths:
        lower = path.lower()
        if lower.endswith(".json"):
            try:
                with open(path, "r") as f:
                    payload = json.load(f)
                if isinstance(payload, dict):
                    results.update(payload)
                else:
                    results.setdefault("files", []).append({"path": path, "content": payload})
            except Exception as exc:
                logger.warning(f"Failed to parse JSON file {path}: {exc}")
        else:
            results.setdefault("files", []).append({"path": path})
    results["downloaded_files"] = file_paths
    return results


class NeuroSnapWrapper:
    """Wrapper orchestrating NeuroSnap service calls and result parsing."""

    def __init__(self, client: Optional[NeuroSnapClient] = None):
        self.client = client or NeuroSnapClient()

    def _run_service(self, service_name: str, fields: Dict[str, Any], note_data: Any,
                     output_dir: Optional[str] = None, max_wait_time: int = 3600) -> Dict[str, Any]:
        note = _get_input_hash(note_data)
        job_id = self.client.find_existing_job(service_name, note)
        if not job_id:
            job_id = self.client.submit_job(service_name, fields, note)
        completed = self.client.wait_for_job_completion(job_id, max_wait_time=max_wait_time)
        if not completed:
            raise RuntimeError(f"NeuroSnap job {job_id} for {service_name} did not complete successfully")
        out_dir = output_dir or os.path.join("results", "neurosnap", service_name.replace(" ", "_"), note)
        files = self.client.download_job_files(job_id, out_dir)
        return _parse_downloaded_files(files)

    # ----------------------------
    # Service-specific wrappers
    # ----------------------------
    def predict_admet(self, smiles: str, properties: Optional[List[str]] = None,
                      max_wait_time: int = 1800) -> Dict[str, Any]:
        service = "ADMET-AI"
        payload = {"smiles": smiles, "properties": properties or []}
        # Many services expect multipart even for JSON fields
        fields = {"Input Molecule": json.dumps([{"data": smiles, "type": "smiles"}])}
        return self._run_service(service, fields, {"service": service, **payload}, max_wait_time=max_wait_time)

    def predict_toxicity(self, smiles: str, max_wait_time: int = 1800) -> Dict[str, Any]:
        service = "eTox"
        fields = {"Input Molecule": json.dumps([{"data": smiles, "type": "smiles"}])}
        return self._run_service(service, fields, {"service": service, "smiles": smiles}, max_wait_time=max_wait_time)

    def predict_aggregation(self, sequence: str, max_wait_time: int = 1800) -> Dict[str, Any]:
        service = "Aggrescan3D"  # keep as-is if present; fallback by services lookup not implemented here
        fasta = f">protein\n{sequence}"
        fields = {"Input Molecule": fasta}
        return self._run_service(service, fields, {"service": service, "sequence": sequence}, max_wait_time=max_wait_time)

    def predict_thermostability(self, sequence: str, temperature: float = 25.0,
                                max_wait_time: int = 1800) -> Dict[str, Any]:
        # Sequence-based alternative per instructions
        service = "TemStaPro"
        fasta = f">protein\n{sequence}"
        fields = {"Input Molecule": fasta}
        return self._run_service(service, fields, {"service": service, "sequence": sequence, "temp": temperature}, max_wait_time=max_wait_time)

    def predict_structure(self, sequences: List[str], max_wait_time: int = 7200) -> Dict[str, Any]:
        service = "Boltz-2 (AlphaFold3)"
        # Convert list to mapping expected by API: {"aa": {"seq1": s1, ...}}
        seq_map = {f"seq{i+1}": seq for i, seq in enumerate(sequences)}
        input_data = {"aa": seq_map}
        fields = {
            "Input Sequences": json.dumps(input_data),
            "Number Recycles": "3",
            "Diffusion Samples": "1",
            "Diffusion Samples Affinity": "1",
        }
        return self._run_service(service, fields, {"service": service, "sequences": seq_map}, max_wait_time=max_wait_time)

    def predict_stability_change(self, pdb_path: str, mutations: List[str], max_wait_time: int = 1800) -> Dict[str, Any]:
        service = "StaB-ddG"
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"PDB file not found: {pdb_path}")
        with open(pdb_path, "rb") as f:
            pdb_bytes = f.read()
        fields = {
            "Input Molecule": (os.path.basename(pdb_path), pdb_bytes, "application/octet-stream"),
            "Mutations": "\n".join(mutations),
        }
        return self._run_service(service, fields, {"service": service, "mutations": sorted(mutations)}, max_wait_time=max_wait_time)

    def predict_immunogenicity(self, sequence: str, max_wait_time: int = 1800) -> Dict[str, Any]:
        service = "DeepImmuno"
        fasta = f">protein\n{sequence}"
        fields = {"Input Molecule": fasta}
        return self._run_service(service, fields, {"service": service, "sequence": sequence}, max_wait_time=max_wait_time)


# ------------------------------
# Convenience functions (legacy)
# ------------------------------

def predict_admet(smiles: str, properties: Optional[List[str]] = None, max_wait_time: int = 1800) -> Dict[str, Any]:
    return NeuroSnapWrapper().predict_admet(smiles, properties, max_wait_time=max_wait_time)


def predict_toxicity(smiles: str, max_wait_time: int = 1800) -> Dict[str, Any]:
    return NeuroSnapWrapper().predict_toxicity(smiles, max_wait_time=max_wait_time)


def predict_aggregation(sequence: str, max_wait_time: int = 1800) -> Dict[str, Any]:
    return NeuroSnapWrapper().predict_aggregation(sequence, max_wait_time=max_wait_time)


def predict_thermostability(sequence: str, temperature: float = 25.0, max_wait_time: int = 1800) -> Dict[str, Any]:
    return NeuroSnapWrapper().predict_thermostability(sequence, temperature=temperature, max_wait_time=max_wait_time)


def predict_structure(sequences: List[str], max_wait_time: int = 7200) -> Dict[str, Any]:
    return NeuroSnapWrapper().predict_structure(sequences, max_wait_time=max_wait_time)


def predict_stability_change(pdb_path: str, mutations: List[str], max_wait_time: int = 1800) -> Dict[str, Any]:
    return NeuroSnapWrapper().predict_stability_change(pdb_path, mutations, max_wait_time=max_wait_time)


def predict_immunogenicity(sequence: str, max_wait_time: int = 1800) -> Dict[str, Any]:
    return NeuroSnapWrapper().predict_immunogenicity(sequence, max_wait_time=max_wait_time)
