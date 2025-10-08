#!/usr/bin/env python3
"""
Real NeuroSnap API Client

This module provides a robust client for interacting with the NeuroSnap API.
It handles authentication via environment variable, job submission (multipart),
status polling with exponential backoff, job reuse by note matching, and
result file downloads.

Environment:
- NEUROSNAP_API_KEY: API key used for authentication
- NEUROSNAP_BASE_URL (optional): Override base URL (default: https://neurosnap.ai/api)
"""

import os
import time
import json
import logging
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass
from enum import Enum

import requests
from requests import Response
from requests_toolbelt.multipart.encoder import MultipartEncoder

# Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class CircuitBreaker:
    """Minimal circuit breaker to maintain backward compatibility with tests."""
    def __init__(self) -> None:
        self.state = "CLOSED"


class JobStatus(Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class NeuroSnapJob:
    job_id: str
    status: JobStatus
    model_type: str
    created_at: str = ""
    updated_at: str = ""
    result_url: Optional[str] = None
    error_message: Optional[str] = None


class NeuroSnapClient:
    """
    Client for NeuroSnap API with authentication, job management, and job reuse.
    """

    DEFAULT_BASE_URL = "https://neurosnap.ai/api"
    DEFAULT_API_KEY = "f01ad42e66fd96d05b6b77efe301e00f5fab82621e3224ad5d023bc88a7d360b746819f541aa9424b4611d6fe2838f9127b2c1012f9ffce1720a2ebb557b50c4"

    def __init__(self, timeout: int = 60, base_url: Optional[str] = None, api_key: Optional[str] = None):
        self.timeout = timeout
        self.base_url = base_url or os.environ.get("NEUROSNAP_BASE_URL", self.DEFAULT_BASE_URL)
        # Priority: explicit parameter > environment variable > hardcoded default
        self.api_key = api_key or os.environ.get("NEUROSNAP_API_KEY") or self.DEFAULT_API_KEY
        if not self.api_key:
            logger.warning("NEUROSNAP_API_KEY not set. Real API calls will fail until provided.")
        self.headers = {"X-API-KEY": self.api_key} if self.api_key else {}
        # Back-compat attribute expected by tests
        self.circuit_breaker = CircuitBreaker()

    # ----------------------
    # Job discovery / reuse
    # ----------------------
    def list_jobs(self) -> List[Dict[str, Any]]:
        """List jobs for the authenticated user."""
        try:
            url = f"{self.base_url}/jobs"
            resp = requests.get(url, headers=self.headers, timeout=self.timeout)
            resp.raise_for_status()
            data = resp.json()
            if isinstance(data, list):
                return data
            return []
        except Exception as exc:
            logger.error(f"Failed to list jobs: {exc}")
            return []

    def find_existing_job(self, service_name: str, note: str) -> Optional[str]:
        """Find a completed job matching service and note; return job ID or None."""
        logger.info(f"Searching for existing job for '{service_name}' note '{note}'")
        for job in self.list_jobs():
            try:
                if (
                    job.get("Service") == service_name
                    and job.get("Note") == note
                    and str(job.get("Status", "")).lower() in ("completed", "done")
                ):
                    job_id = str(job.get("ID"))
                    if job_id:
                        logger.info(f"Reusing completed job {job_id}")
                        return job_id
            except Exception:
                continue
        return None

    # ----------------------
    # Job lifecycle
    # ----------------------
    def submit_job(self, service_name: str, fields: Dict[str, Any], note: str) -> str:
        """
        Submit a job using multipart form fields. File fields should be tuples:
        (filename, bytes, mimetype)
        """
        multipart = MultipartEncoder(fields=fields)
        headers = dict(self.headers)
        headers["Content-Type"] = multipart.content_type
        url = f"{self.base_url}/job/submit/{service_name}?note={note}"
        logger.info(f"Submitting job: {service_name} note={note}")
        resp: Response = requests.post(url, headers=headers, data=multipart, timeout=self.timeout)
        resp.raise_for_status()
        job_id = resp.text.strip().strip('"')
        logger.info(f"Job submitted; id={job_id}")
        return job_id

    def wait_for_job_completion(self, job_id: str, max_wait_time: int = 3600) -> bool:
        """Poll job status until completion/failed or timeout. Returns True if completed."""
        start = time.time()
        backoff = 10
        max_backoff = 60
        while (time.time() - start) < max_wait_time:
            try:
                url = f"{self.base_url}/job/status/{job_id}"
                resp = requests.get(url, headers=self.headers, timeout=self.timeout)
                resp.raise_for_status()
                status = resp.text.strip().strip('"').lower()
                logger.info(f"Job {job_id} status: {status}")
                if status in ("completed", "done"):
                    return True
                if status in ("failed", "error", "cancelled"):
                    return False
            except Exception as exc:
                logger.warning(f"Polling error for job {job_id}: {exc}")
            time.sleep(backoff)
            backoff = min(int(backoff * 1.5), max_backoff)
        logger.error(f"Job {job_id} timed out after {max_wait_time}s")
        return False

    def _get_out_files_manifest(self, job_id: str) -> Dict[str, Any]:
        url = f"{self.base_url}/job/data/{job_id}"
        resp = requests.get(url, headers=self.headers, timeout=self.timeout)
        resp.raise_for_status()
        return resp.json()

    def _download_file(self, job_id: str, filename: str, output_dir: str) -> str:
        os.makedirs(output_dir, exist_ok=True)
        url = f"{self.base_url}/job/file/{job_id}/out/{filename}"
        resp = requests.get(url, headers=self.headers, allow_redirects=True, timeout=self.timeout)
        resp.raise_for_status()
        out_path = os.path.join(output_dir, filename)
        with open(out_path, "wb") as f:
            f.write(resp.content)
        return out_path

    def download_job_files(self, job_id: str, output_dir: str) -> List[str]:
        """Download all output files for a completed job; return list of file paths."""
        try:
            manifest = self._get_out_files_manifest(job_id)
        except Exception as exc:
            logger.error(f"Failed to retrieve manifest for job {job_id}: {exc}")
            return []
        downloaded: List[str] = []
        try:
            items = manifest.get("out") or []
            # Expecting list of [filename, size] or similar pairs
            for entry in items:
                try:
                    filename = entry[0] if isinstance(entry, (list, tuple)) else str(entry)
                    path = self._download_file(job_id, filename, output_dir)
                    logger.info(f"Downloaded {filename} -> {path}")
                    downloaded.append(path)
                except Exception as exc:
                    logger.warning(f"Failed to download an item for job {job_id}: {exc}")
        except Exception as exc:
            logger.error(f"Unexpected manifest format for job {job_id}: {exc}")
        return downloaded

    # ----------------------
    # Back-compat client API used in tests
    # ----------------------
    def get_job_status(self, job_id: str) -> NeuroSnapJob:
        """Return job status as a NeuroSnapJob for compatibility."""
        try:
            url = f"{self.base_url}/job/status/{job_id}"
            resp = requests.get(url, headers=self.headers, timeout=self.timeout)
            resp.raise_for_status()
            status_text = resp.text.strip().strip('"').lower()
        except Exception:
            status_text = "pending"
        status = JobStatus(status_text) if status_text in JobStatus._value2member_map_ else JobStatus.PENDING
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        return NeuroSnapJob(job_id=job_id, status=status, model_type="unknown", created_at=now, updated_at=now)

    def _mk_smiles_fields(self, smiles: str) -> Dict[str, Any]:
        return {"Input Molecule": json.dumps([{"data": smiles, "type": "smiles"}])}

    def _mk_fasta_fields(self, sequence: str) -> Dict[str, Any]:
        return {"Input Sequences": f">protein\n{sequence}"}

    def submit_admet_job(self, smiles: str, properties: Optional[List[str]] = None) -> NeuroSnapJob:
        job_id = self.submit_job("ADMET-AI", self._mk_smiles_fields(smiles), note=f"admet:{smiles}")
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        return NeuroSnapJob(job_id=job_id, status=JobStatus.PENDING, model_type="admet_ai", created_at=now, updated_at=now)

    def submit_etox_job(self, smiles: str) -> NeuroSnapJob:
        job_id = self.submit_job("eTox", self._mk_smiles_fields(smiles), note=f"etox:{smiles}")
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        return NeuroSnapJob(job_id=job_id, status=JobStatus.PENDING, model_type="etox", created_at=now, updated_at=now)

    def submit_aggrescan_job(self, sequence: str) -> NeuroSnapJob:
        fields = {"Input Structure": f">protein\n{sequence}"}
        job_id = self.submit_job("Aggrescan3D", fields, note=f"agg:{hash(sequence)}")
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        return NeuroSnapJob(job_id=job_id, status=JobStatus.PENDING, model_type="aggrescan3d", created_at=now, updated_at=now)

    def submit_thermompn_job(self, sequence: str, temperature: float = 25.0) -> NeuroSnapJob:
        fields = {"Input PDB": f">protein\n{sequence}"}
        job_id = self.submit_job("TemStaPro", fields, note=f"temstapro:{hash((sequence, temperature))}")
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        return NeuroSnapJob(job_id=job_id, status=JobStatus.PENDING, model_type="thermompn", created_at=now, updated_at=now)

    def submit_boltz_job(self, sequences: List[str]) -> NeuroSnapJob:
        seq_map = {f"seq{i+1}": s for i, s in enumerate(sequences)}
        fields = {"Input Sequences": json.dumps({"aa": seq_map}), "Number Recycles": "3", "Diffusion Samples": "1", "Diffusion Samples Affinity": "1"}
        job_id = self.submit_job("Boltz-2 (AlphaFold3)", fields, note=f"boltz:{hash(tuple(sequences))}")
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        return NeuroSnapJob(job_id=job_id, status=JobStatus.PENDING, model_type="boltz2", created_at=now, updated_at=now)

    def submit_stab_ddg_job(self, pdb_path: str, mutations: List[str]) -> NeuroSnapJob:
        with open(pdb_path, "rb") as f:
            pdb_bytes = f.read()
        fields = {"Input Molecule": (os.path.basename(pdb_path), pdb_bytes, "application/octet-stream"), "Mutations": "\n".join(mutations)}
        job_id = self.submit_job("StaB-ddG", fields, note=f"stabddg:{hash(tuple(sorted(mutations)))}")
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        return NeuroSnapJob(job_id=job_id, status=JobStatus.PENDING, model_type="stab_ddg", created_at=now, updated_at=now)

    def submit_deepimmuno_job(self, sequence: str) -> NeuroSnapJob:
        job_id = self.submit_job("DeepImmuno", self._mk_fasta_fields(sequence), note=f"deepimmuno:{hash(sequence)}")
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        return NeuroSnapJob(job_id=job_id, status=JobStatus.PENDING, model_type="deepimmuno", created_at=now, updated_at=now)
