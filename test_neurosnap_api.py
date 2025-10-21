#!/usr/bin/env python3
"""
Quick test to verify NeuroSnap API connectivity and authentication
"""

import sys
from scripts.neurosnap_client import NeuroSnapClient

def test_api_connection():
    """Test basic API connectivity"""
    print("Testing NeuroSnap API connection...")

    client = NeuroSnapClient()

    # Check API key is set
    if not client.api_key:
        print("❌ ERROR: No API key configured")
        return False

    print(f"✓ API key configured (ends with: ...{client.api_key[-8:]})")
    print(f"✓ Base URL: {client.base_url}")

    # Try to list jobs
    try:
        jobs = client.list_jobs()
        print(f"✓ API connection successful!")
        print(f"✓ Found {len(jobs)} existing jobs")

        if jobs:
            print("\nRecent jobs:")
            for job in jobs[:5]:  # Show first 5
                print(f"  - Job {job.get('ID')}: {job.get('Service')} ({job.get('Status')})")

        return True

    except Exception as e:
        print(f"❌ API connection failed: {e}")
        return False

if __name__ == "__main__":
    success = test_api_connection()
    sys.exit(0 if success else 1)
