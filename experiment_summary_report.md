# TL1A PET Imaging Agent - Experiment Run Summary

**Total Experiments:** 15
**Successful:** 12
**Failed:** 3

## Results

| Experiment | Status | Details |
|------------|--------|---------|
| 01_physchem | ✅ Pass | Completed successfully |
| 02_metabolism | ❌ Fail | Script failed with return code 1<br>Stderr: Traceback (most recent call last):<br>  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/02_metabolism.py", line 11, in <module><br>    import scripts.inputs as inputs<br>ModuleNotFoundError: No module named 'scripts'<br> |
| 03_toxicity_flags | ✅ Pass | Completed successfully |
| 04_linker_flex | ✅ Pass | Completed successfully |
| 05_mmp_analysis | ✅ Pass | Completed successfully |
| 06_sequence_dev | ✅ Pass | Completed successfully |
| 07_agg_hotspots | ✅ Pass | Completed successfully |
| 08_charge_pI | ✅ Pass | Completed successfully |
| 09_flexibility_anm | ✅ Pass | Completed successfully |
| 10_thermostability | ✅ Pass | Completed successfully |
| 11_complex_model | ✅ Pass | Completed successfully |
| 12_interface_fingerprint | ✅ Pass | Completed successfully |
| 13_ala_scanning | ❌ Fail | Script failed with return code 1<br>Stderr: 2025-10-07 20:23:08,568 - INFO - Performing alanine scanning for Fab06_VH<br>2025-10-07 20:23:08,568 - INFO - Searching for existing job for 'StaB-ddG' note 'f2a3f3a5580c825aeb75d915fc0946395d74ae640c5385efc2a4642e2ea2b257'<br>2025-10-07 20:23:08,979 - INFO - Submitting job: StaB-ddG note=f2a3f3a5580c825aeb75d915fc0946395d74ae640c5385efc2a4642e2ea2b257<br>2025-10-07 20:23:10,150 - ERROR - StaB-ddG prediction failed for Fab06_VH: 400 Client Error: Bad Request for url: https://neurosnap.ai/api/job/submit/StaB-ddG?note=f2a3f3a5580c825aeb75d915fc0946395d74ae640c5385efc2a4642e2ea2b257<br>2025-10-07 20:23:10,150 - ERROR - Error in Experiment 13: Alanine scanning failed for Fab06_VH: 400 Client Error: Bad Request for url: https://neurosnap.ai/api/job/submit/StaB-ddG?note=f2a3f3a5580c825aeb75d915fc0946395d74ae640c5385efc2a4642e2ea2b257<br>Traceback (most recent call last):<br>  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/13_ala_scanning.py", line 67, in perform_alanine_scanning<br>    stab_results = ns_predict_stability_change(pdb_path, all_mutations, max_wait_time=1800)<br>  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/neurosnap_wrappers.py", line 166, in predict_stability_change<br>    return NeuroSnapWrapper().predict_stability_change(pdb_path, mutations, max_wait_time=max_wait_time)<br>           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^<br>  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/neurosnap_wrappers.py", line 132, in predict_stability_change<br>    return self._run_service(service, fields, {"service": service, "mutations": sorted(mutations)}, max_wait_time=max_wait_time)<br>           ~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^<br>  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/neurosnap_wrappers.py", line 71, in _run_service<br>    job_id = self.client.submit_job(service_name, fields, note)<br>  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/neurosnap_client.py", line 125, in submit_job<br>    resp.raise_for_status()<br>    ~~~~~~~~~~~~~~~~~~~~~^^<br>  File "/opt/homebrew/Caskroom/miniconda/base/lib/python3.13/site-packages/requests/models.py", line 1026, in raise_for_status<br>    raise HTTPError(http_error_msg, response=self)<br>requests.exceptions.HTTPError: 400 Client Error: Bad Request for url: https://neurosnap.ai/api/job/submit/StaB-ddG?note=f2a3f3a5580c825aeb75d915fc0946395d74ae640c5385efc2a4642e2ea2b257<br><br>During handling of the above exception, another exception occurred:<br><br>Traceback (most recent call last):<br>  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/13_ala_scanning.py", line 154, in <module><br>    main()<br>    ~~~~^^<br>  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/13_ala_scanning.py", line 117, in main<br>    results = perform_alanine_scanning(name, seq)<br>  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/13_ala_scanning.py", line 99, in perform_alanine_scanning<br>    raise RuntimeError(f"Alanine scanning failed for {fab_name}: {str(e)}")<br>RuntimeError: Alanine scanning failed for Fab06_VH: 400 Client Error: Bad Request for url: https://neurosnap.ai/api/job/submit/StaB-ddG?note=f2a3f3a5580c825aeb75d915fc0946395d74ae640c5385efc2a4642e2ea2b257<br> |
| 14_immunogenicity | ✅ Pass | Completed successfully |
| 15_decision_scorecard | ❌ Fail | Script failed with return code 1<br>Stderr: 2025-10-07 20:23:13,954 - INFO - Generating integrated decision scorecard<br>2025-10-07 20:23:13,957 - INFO - Loaded physchem results: 2 rows<br>2025-10-07 20:23:13,957 - WARNING - Missing results file: results/prodrug/02_biotransformer_metabolites.csv<br>2025-10-07 20:23:13,958 - INFO - Loaded toxicity results: 2 rows<br>2025-10-07 20:23:13,959 - INFO - Loaded flexibility results: 2 rows<br>2025-10-07 20:23:13,959 - INFO - Loaded mmp results: 1 rows<br>2025-10-07 20:23:13,960 - INFO - Loaded developability results: 4 rows<br>2025-10-07 20:23:13,960 - INFO - Loaded aggregation results: 4 rows<br>2025-10-07 20:23:13,960 - INFO - Loaded charge results: 4 rows<br>2025-10-07 20:23:13,961 - INFO - Loaded anm_flex results: 1 rows<br>2025-10-07 20:23:13,961 - INFO - Loaded thermostability results: 4 rows<br>2025-10-07 20:23:13,961 - INFO - Loaded complex results: 2 rows<br>2025-10-07 20:23:13,961 - INFO - Loaded interface results: 1 rows<br>2025-10-07 20:23:13,961 - WARNING - Missing results file: results/biobetter/13_scanning_summary.csv<br>2025-10-07 20:23:13,962 - INFO - Loaded immunogenicity results: 4 rows<br>2025-10-07 20:23:13,963 - ERROR - Error in Experiment 15: 'Thermo_Score'<br>Traceback (most recent call last):<br>  File "/opt/homebrew/Caskroom/miniconda/base/lib/python3.13/site-packages/pandas/core/indexes/base.py", line 3805, in get_loc<br>    return self._engine.get_loc(casted_key)<br>           ~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^<br>  File "index.pyx", line 167, in pandas._libs.index.IndexEngine.get_loc<br>  File "index.pyx", line 196, in pandas._libs.index.IndexEngine.get_loc<br>  File "pandas/_libs/hashtable_class_helper.pxi", line 7081, in pandas._libs.hashtable.PyObjectHashTable.get_item<br>  File "pandas/_libs/hashtable_class_helper.pxi", line 7089, in pandas._libs.hashtable.PyObjectHashTable.get_item<br>KeyError: 'Thermo_Score'<br><br>The above exception was the direct cause of the following exception:<br><br>Traceback (most recent call last):<br>  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/15_decision_scorecard.py", line 217, in <module><br>    main()<br>    ~~~~^^<br>  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/15_decision_scorecard.py", line 193, in main<br>    scorecard_df = calculate_decision_scores(experiment_results)<br>  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/15_decision_scorecard.py", line 116, in calculate_decision_scores<br>    scores["Thermo_Score"] = fab_row["Thermo_Score"].values[0]<br>                             ~~~~~~~^^^^^^^^^^^^^^^^<br>  File "/opt/homebrew/Caskroom/miniconda/base/lib/python3.13/site-packages/pandas/core/frame.py", line 4102, in __getitem__<br>    indexer = self.columns.get_loc(key)<br>  File "/opt/homebrew/Caskroom/miniconda/base/lib/python3.13/site-packages/pandas/core/indexes/base.py", line 3812, in get_loc<br>    raise KeyError(key) from err<br>KeyError: 'Thermo_Score'<br> |

## Failed Experiments

### 02_metabolism: Activation & Metabolism Prediction
**Error:** Script failed with return code 1
Stderr: Traceback (most recent call last):
  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/02_metabolism.py", line 11, in <module>
    import scripts.inputs as inputs
ModuleNotFoundError: No module named 'scripts'


### 13_ala_scanning: Alanine Scanning
**Error:** Script failed with return code 1
Stderr: 2025-10-07 20:23:08,568 - INFO - Performing alanine scanning for Fab06_VH
2025-10-07 20:23:08,568 - INFO - Searching for existing job for 'StaB-ddG' note 'f2a3f3a5580c825aeb75d915fc0946395d74ae640c5385efc2a4642e2ea2b257'
2025-10-07 20:23:08,979 - INFO - Submitting job: StaB-ddG note=f2a3f3a5580c825aeb75d915fc0946395d74ae640c5385efc2a4642e2ea2b257
2025-10-07 20:23:10,150 - ERROR - StaB-ddG prediction failed for Fab06_VH: 400 Client Error: Bad Request for url: https://neurosnap.ai/api/job/submit/StaB-ddG?note=f2a3f3a5580c825aeb75d915fc0946395d74ae640c5385efc2a4642e2ea2b257
2025-10-07 20:23:10,150 - ERROR - Error in Experiment 13: Alanine scanning failed for Fab06_VH: 400 Client Error: Bad Request for url: https://neurosnap.ai/api/job/submit/StaB-ddG?note=f2a3f3a5580c825aeb75d915fc0946395d74ae640c5385efc2a4642e2ea2b257
Traceback (most recent call last):
  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/13_ala_scanning.py", line 67, in perform_alanine_scanning
    stab_results = ns_predict_stability_change(pdb_path, all_mutations, max_wait_time=1800)
  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/neurosnap_wrappers.py", line 166, in predict_stability_change
    return NeuroSnapWrapper().predict_stability_change(pdb_path, mutations, max_wait_time=max_wait_time)
           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/neurosnap_wrappers.py", line 132, in predict_stability_change
    return self._run_service(service, fields, {"service": service, "mutations": sorted(mutations)}, max_wait_time=max_wait_time)
           ~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/neurosnap_wrappers.py", line 71, in _run_service
    job_id = self.client.submit_job(service_name, fields, note)
  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/neurosnap_client.py", line 125, in submit_job
    resp.raise_for_status()
    ~~~~~~~~~~~~~~~~~~~~~^^
  File "/opt/homebrew/Caskroom/miniconda/base/lib/python3.13/site-packages/requests/models.py", line 1026, in raise_for_status
    raise HTTPError(http_error_msg, response=self)
requests.exceptions.HTTPError: 400 Client Error: Bad Request for url: https://neurosnap.ai/api/job/submit/StaB-ddG?note=f2a3f3a5580c825aeb75d915fc0946395d74ae640c5385efc2a4642e2ea2b257

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/13_ala_scanning.py", line 154, in <module>
    main()
    ~~~~^^
  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/13_ala_scanning.py", line 117, in main
    results = perform_alanine_scanning(name, seq)
  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/13_ala_scanning.py", line 99, in perform_alanine_scanning
    raise RuntimeError(f"Alanine scanning failed for {fab_name}: {str(e)}")
RuntimeError: Alanine scanning failed for Fab06_VH: 400 Client Error: Bad Request for url: https://neurosnap.ai/api/job/submit/StaB-ddG?note=f2a3f3a5580c825aeb75d915fc0946395d74ae640c5385efc2a4642e2ea2b257


### 15_decision_scorecard: Integrated Decision Scorecard
**Error:** Script failed with return code 1
Stderr: 2025-10-07 20:23:13,954 - INFO - Generating integrated decision scorecard
2025-10-07 20:23:13,957 - INFO - Loaded physchem results: 2 rows
2025-10-07 20:23:13,957 - WARNING - Missing results file: results/prodrug/02_biotransformer_metabolites.csv
2025-10-07 20:23:13,958 - INFO - Loaded toxicity results: 2 rows
2025-10-07 20:23:13,959 - INFO - Loaded flexibility results: 2 rows
2025-10-07 20:23:13,959 - INFO - Loaded mmp results: 1 rows
2025-10-07 20:23:13,960 - INFO - Loaded developability results: 4 rows
2025-10-07 20:23:13,960 - INFO - Loaded aggregation results: 4 rows
2025-10-07 20:23:13,960 - INFO - Loaded charge results: 4 rows
2025-10-07 20:23:13,961 - INFO - Loaded anm_flex results: 1 rows
2025-10-07 20:23:13,961 - INFO - Loaded thermostability results: 4 rows
2025-10-07 20:23:13,961 - INFO - Loaded complex results: 2 rows
2025-10-07 20:23:13,961 - INFO - Loaded interface results: 1 rows
2025-10-07 20:23:13,961 - WARNING - Missing results file: results/biobetter/13_scanning_summary.csv
2025-10-07 20:23:13,962 - INFO - Loaded immunogenicity results: 4 rows
2025-10-07 20:23:13,963 - ERROR - Error in Experiment 15: 'Thermo_Score'
Traceback (most recent call last):
  File "/opt/homebrew/Caskroom/miniconda/base/lib/python3.13/site-packages/pandas/core/indexes/base.py", line 3805, in get_loc
    return self._engine.get_loc(casted_key)
           ~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^
  File "index.pyx", line 167, in pandas._libs.index.IndexEngine.get_loc
  File "index.pyx", line 196, in pandas._libs.index.IndexEngine.get_loc
  File "pandas/_libs/hashtable_class_helper.pxi", line 7081, in pandas._libs.hashtable.PyObjectHashTable.get_item
  File "pandas/_libs/hashtable_class_helper.pxi", line 7089, in pandas._libs.hashtable.PyObjectHashTable.get_item
KeyError: 'Thermo_Score'

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/15_decision_scorecard.py", line 217, in <module>
    main()
    ~~~~^^
  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/15_decision_scorecard.py", line 193, in main
    scorecard_df = calculate_decision_scores(experiment_results)
  File "/Users/nicholasharris/Documents/Vivamed in silico/v25_clones/TL1A-PET-Imaging-Agent-for-IBD-1.455-v4/scripts/15_decision_scorecard.py", line 116, in calculate_decision_scores
    scores["Thermo_Score"] = fab_row["Thermo_Score"].values[0]
                             ~~~~~~~^^^^^^^^^^^^^^^^
  File "/opt/homebrew/Caskroom/miniconda/base/lib/python3.13/site-packages/pandas/core/frame.py", line 4102, in __getitem__
    indexer = self.columns.get_loc(key)
  File "/opt/homebrew/Caskroom/miniconda/base/lib/python3.13/site-packages/pandas/core/indexes/base.py", line 3812, in get_loc
    raise KeyError(key) from err
KeyError: 'Thermo_Score'


