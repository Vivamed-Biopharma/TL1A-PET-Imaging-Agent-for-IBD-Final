#!/usr/bin/env python3
"""
Main CLI script for TL1A PET Imaging Agent Computational Platform
"""

import argparse
import subprocess
import sys
import os
from pathlib import Path

def run_experiment(exp_num):
    """Run a specific experiment script"""
    script_path = f"scripts/{exp_num:02d}_*.py"
    # Find the script
    scripts_dir = Path("scripts")
    matching_scripts = list(scripts_dir.glob(f"{exp_num:02d}_*.py"))
    if not matching_scripts:
        print(f"No script found for experiment {exp_num}")
        return
    script = matching_scripts[0]
    print(f"Running {script}")
    result = subprocess.run([sys.executable, str(script)], cwd=".")
    return result.returncode

def main():
    parser = argparse.ArgumentParser(description="TL1A PET Imaging Agent Platform")
    parser.add_argument("experiment", type=int, choices=range(1,16), help="Experiment number to run (1-15)")
    parser.add_argument("--all", action="store_true", help="Run all experiments")

    args = parser.parse_args()

    if args.all:
        for i in range(1,16):
            code = run_experiment(i)
            if code != 0:
                print(f"Experiment {i} failed")
                sys.exit(1)
    else:
        code = run_experiment(args.experiment)
        sys.exit(code)

if __name__ == "__main__":
    main()