#!/usr/bin/env python3
"""
NeuroSnap API Client

This module provides a client for interacting with the NeuroSnap API,
including job submission, polling, and result retrieval.
"""

import requests
import time
import json
import logging
import os
from typing import Dict, Any, Optional, List
from dataclasses import dataclass
from enum import Enum
import random

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class CircuitBreaker:
    """
    Circuit breaker pattern for handling API failures.
    """

    def __init__(self, failure_threshold: int = 5, recovery_timeout: int = 60):
        """
        Initialize circuit breaker.

        Args:
            failure_threshold: Number of failures before opening circuit
            recovery_timeout: Time to wait before trying again
        """
        self.failure_threshold = failure_threshold
        self.recovery_timeout = recovery_timeout
        self.failure_count = 0
        self.last_failure_time = None
        self.state = "CLOSED"  # CLOSED, OPEN, HALF_OPEN

    def call(self, func, *args, **kwargs):
        """
        Execute function with circuit breaker protection.

        Args:
            func: Function to call
            *args: Positional arguments
            **kwargs: Keyword arguments

        Returns:
            Function result

        Raises:
            Exception: If circuit is open or function fails
        """
        if self.state == "OPEN":
            if time.time() - self.last_failure_time > self.recovery_timeout:
                self.state = "HALF_OPEN"
                logger.info("Circuit breaker entering HALF_OPEN state")
            else:
                raise Exception("Circuit breaker is OPEN")

        try:
            result = func(*args, **kwargs)
            if self.state == "HALF_OPEN":
                self.state = "CLOSED"
                self.failure_count = 0
                logger.info("Circuit breaker reset to CLOSED")
            return result
        except Exception as e:
            self.failure_count += 1
            self.last_failure_time = time.time()
            if self.failure_count >= self.failure_threshold:
                self.state = "OPEN"
                logger.error(f"Circuit breaker opened after {self.failure_count} failures")
            raise e

class JobStatus(Enum):
    """Job status enumeration."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"

@dataclass
class NeuroSnapJob:
    """Represents a NeuroSnap job."""
    job_id: str
    status: JobStatus
    model_type: str
    created_at: str
    updated_at: str
    result_url: Optional[str] = None
    error_message: Optional[str] = None

class NeuroSnapClient:
    """
    Client for NeuroSnap API with authentication and job management.
    """

    BASE_URL = "https://api.neurosnap.ai"  # Placeholder - adjust as needed
    API_KEY = os.environ.get('NEUROSNAP_API_KEY', 'bd0e1ed66ab2b0e73dfa1d2eba2ddf5d5aaa39d90d7d751547f81d616ddcdc30565d092bc4021535b504e6c30df5eb09a84bffc02ddafc6fb1c7abb53f123c1b')

    def __init__(self, base_url: Optional[str] = None, timeout: int = 30):
        """
        Initialize the NeuroSnap client.

        Args:
            base_url: Base URL for the API (optional)
            timeout: Request timeout in seconds
        """
        self.base_url = base_url or self.BASE_URL
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({
            "Authorization": f"Bearer {self.API_KEY}",
            "Content-Type": "application/json",
            "User-Agent": "TL1A-Platform/1.0"
        })
        self.circuit_breaker = CircuitBreaker()

    def _make_request(self, method: str, endpoint: str, data: Optional[Dict] = None, **kwargs) -> Dict[str, Any]:
        """
        Make an authenticated request to the API with circuit breaker protection.

        Args:
            method: HTTP method
            endpoint: API endpoint
            data: Request data
            **kwargs: Additional request parameters

        Returns:
            Response data as dictionary

        Raises:
            requests.RequestException: On request failure
            ValueError: On API error response
        """
        def _request():
            url = f"{self.base_url}{endpoint}"

            logger.debug(f"Making request: {method} {url}")

            if data:
                kwargs['json'] = data

            response = self.session.request(method, url, timeout=self.timeout, **kwargs)
            response.raise_for_status()

            # Validate response
            if not response.content:
                raise ValueError("Empty response from API")

            try:
                result = response.json()
                self._validate_response(result)
                return result
            except json.JSONDecodeError as e:
                logger.error(f"Invalid JSON response: {str(e)}")
                raise ValueError(f"Invalid API response: {str(e)}")

        return self.circuit_breaker.call(_request)

    def _validate_response(self, response: Dict[str, Any]):
        """
        Validate API response structure.

        Args:
            response: Response dictionary

        Raises:
            ValueError: If response is invalid
        """
        if not isinstance(response, dict):
            raise ValueError("Response must be a dictionary")

        # Add specific validation based on expected response structure
        # This can be extended based on API documentation

    def submit_job(self, model_type: str, input_data: Dict[str, Any]) -> NeuroSnapJob:
        """
        Submit a job to NeuroSnap.

        Args:
            model_type: Type of model to run
            input_data: Input data for the model

        Returns:
            NeuroSnapJob: Submitted job information

        Raises:
            ValueError: On invalid input or API error
        """
        logger.info(f"Submitting job for model: {model_type}")

        request_data = {
            "model_type": model_type,
            "input": input_data
        }

        try:
            response = self._make_request("POST", "/jobs", request_data)

            job = NeuroSnapJob(
                job_id=response["job_id"],
                status=JobStatus(response["status"]),
                model_type=response["model_type"],
                created_at=response["created_at"],
                updated_at=response["updated_at"]
            )

            logger.info(f"Job submitted successfully: {job.job_id}")
            return job

        except KeyError as e:
            raise ValueError(f"Invalid API response structure: missing {str(e)}")

    def get_job_status(self, job_id: str) -> NeuroSnapJob:
        """
        Get the status of a job.

        Args:
            job_id: Job ID to check

        Returns:
            NeuroSnapJob: Updated job information

        Raises:
            ValueError: On invalid job ID or API error
        """
        try:
            response = self._make_request("GET", f"/jobs/{job_id}")

            job = NeuroSnapJob(
                job_id=response["job_id"],
                status=JobStatus(response["status"]),
                model_type=response["model_type"],
                created_at=response["created_at"],
                updated_at=response["updated_at"],
                result_url=response.get("result_url"),
                error_message=response.get("error_message")
            )

            return job

        except KeyError as e:
            raise ValueError(f"Invalid API response structure: missing {str(e)}")

    def wait_for_job_completion(self, job_id: str, max_wait_time: int = 3600,
                               initial_poll_interval: int = 5, max_poll_interval: int = 60) -> NeuroSnapJob:
        """
        Wait for a job to complete with exponential backoff and jitter.

        Args:
            job_id: Job ID to wait for
            max_wait_time: Maximum time to wait in seconds
            initial_poll_interval: Initial polling interval in seconds
            max_poll_interval: Maximum polling interval in seconds

        Returns:
            NeuroSnapJob: Completed job information

        Raises:
            TimeoutError: If job doesn't complete within max_wait_time
            RuntimeError: If job fails
        """
        start_time = time.time()
        current_interval = initial_poll_interval
        attempt = 0

        while time.time() - start_time < max_wait_time:
            try:
                job = self.get_job_status(job_id)
                attempt += 1

                logger.info(f"Job {job_id} status: {job.status.value} (attempt {attempt})")

                if job.status == JobStatus.COMPLETED:
                    logger.info(f"Job {job_id} completed successfully")
                    return job
                elif job.status == JobStatus.FAILED:
                    error_msg = job.error_message or "Unknown error"
                    logger.error(f"Job {job_id} failed: {error_msg}")
                    raise RuntimeError(f"Job failed: {error_msg}")
                elif job.status == JobStatus.CANCELLED:
                    logger.error(f"Job {job_id} was cancelled")
                    raise RuntimeError("Job was cancelled")

                # Exponential backoff with jitter
                sleep_time = min(current_interval, max_poll_interval)
                sleep_time += random.uniform(0, sleep_time * 0.1)  # Add 10% jitter

                logger.debug(f"Sleeping for {sleep_time:.2f} seconds")
                time.sleep(sleep_time)
                current_interval = min(current_interval * 1.5, max_poll_interval)

            except Exception as e:
                logger.warning(f"Error polling job status (attempt {attempt}): {str(e)}")
                # Continue polling despite errors, but with increased interval
                time.sleep(current_interval)
                current_interval = min(current_interval * 1.2, max_poll_interval)

        raise TimeoutError(f"Job {job_id} did not complete within {max_wait_time} seconds")

    def get_job_results(self, job: NeuroSnapJob) -> Dict[str, Any]:
        """
        Retrieve results for a completed job.

        Args:
            job: Completed job information

        Returns:
            Dictionary containing job results

        Raises:
            ValueError: If job is not completed or result_url is missing
            RuntimeError: On result retrieval failure
        """
        if job.status != JobStatus.COMPLETED:
            raise ValueError(f"Job {job.job_id} is not completed (status: {job.status.value})")

        if not job.result_url:
            raise ValueError(f"No result URL available for job {job.job_id}")

        try:
            logger.info(f"Retrieving results for job {job.job_id}")
            response = self._make_request("GET", job.result_url)

            return response

        except Exception as e:
            logger.error(f"Failed to retrieve results for job {job.job_id}: {str(e)}")
            raise RuntimeError(f"Result retrieval failed: {str(e)}")

    def submit_admet_job(self, smiles: str, properties: List[str] = None) -> NeuroSnapJob:
        """
        Submit ADMET prediction job.

        Args:
            smiles: SMILES string of the molecule
            properties: List of properties to predict (optional)

        Returns:
            NeuroSnapJob: Submitted job
        """
        input_data = {
            "smiles": smiles,
            "properties": properties or ["solubility", "permeability", "toxicity", "metabolism"]
        }
        return self.submit_job("admet_ai", input_data)

    def submit_etox_job(self, smiles: str) -> NeuroSnapJob:
        """
        Submit eTox toxicity prediction job.

        Args:
            smiles: SMILES string of the molecule

        Returns:
            NeuroSnapJob: Submitted job
        """
        input_data = {"smiles": smiles}
        return self.submit_job("etox", input_data)

    def submit_aggrescan_job(self, sequence: str) -> NeuroSnapJob:
        """
        Submit Aggrescan3D aggregation prediction job.

        Args:
            sequence: Protein sequence

        Returns:
            NeuroSnapJob: Submitted job
        """
        input_data = {"sequence": sequence}
        return self.submit_job("aggrescan3d", input_data)

    def submit_thermompn_job(self, sequence: str, temperature: float = 25.0) -> NeuroSnapJob:
        """
        Submit ThermoMPNN thermostability prediction job.

        Args:
            sequence: Protein sequence
            temperature: Temperature in Celsius

        Returns:
            NeuroSnapJob: Submitted job
        """
        input_data = {
            "sequence": sequence,
            "temperature": temperature
        }
        return self.submit_job("thermompn", input_data)

    def submit_boltz_job(self, sequences: List[str]) -> NeuroSnapJob:
        """
        Submit Boltz-2 structure prediction job.

        Args:
            sequences: List of protein sequences

        Returns:
            NeuroSnapJob: Submitted job
        """
        input_data = {"sequences": sequences}
        return self.submit_job("boltz2", input_data)

    def submit_stab_ddg_job(self, sequence: str, mutations: List[str]) -> NeuroSnapJob:
        """
        Submit StaB-ddG stability prediction job.

        Args:
            sequence: Wild-type protein sequence
            mutations: List of mutations (e.g., ["A1V", "K2M"])

        Returns:
            NeuroSnapJob: Submitted job
        """
        input_data = {
            "sequence": sequence,
            "mutations": mutations
        }
        return self.submit_job("stab_ddg", input_data)

    def submit_deepimmuno_job(self, sequence: str) -> NeuroSnapJob:
        """
        Submit DeepImmuno immunogenicity prediction job.

        Args:
            sequence: Protein sequence

        Returns:
            NeuroSnapJob: Submitted job
        """
        input_data = {"sequence": sequence}
        return self.submit_job("deepimmuno", input_data)

    def submit_and_wait(self, model_type: str, input_data: Dict[str, Any],
                       max_wait_time: int = 3600) -> Dict[str, Any]:
        """
        Submit a job and wait for completion.

        Args:
            model_type: Type of model to run
            input_data: Input data for the model
            max_wait_time: Maximum time to wait for completion

        Returns:
            Job results as dictionary
        """
        job = self.submit_job(model_type, input_data)
        completed_job = self.wait_for_job_completion(job.job_id, max_wait_time)
        results = self.get_job_results(completed_job)

        return results