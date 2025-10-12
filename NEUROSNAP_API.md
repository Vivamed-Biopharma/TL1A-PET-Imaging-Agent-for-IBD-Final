

# NeuroSnap API Documentation (VERIFIED & TESTED)

## ⚠️ CRITICAL CORRECTIONS (Verified 2025-10-09)

**The original NeuroSnap documentation has ERRORS. Use these corrections:**

1. **DiffDock-L Receptor:** Must be file-like tuple, NOT JSON
2. **Job Data Endpoint:** Use `/job/data/`, NOT `/job/files/`
3. **Number Samples:** Must be STRING ("100"), not integer
4. **ADMET-AI Output:** Returns `results.csv`, not `output.json`

## Authentication & Base URL
- **Base URL**: `https://neurosnap.ai/api/`
- **Authentication**: Use `X-API-KEY` header with your API key
- **Content-Type**: `multipart/form-data` for job submissions

## ✅ VERIFIED WORKING EXAMPLES

### DiffDock-L (CORRECT FORMAT)
```python
from io import BytesIO
from requests_toolbelt.multipart.encoder import MultipartEncoder
import requests, json

API_KEY = "your_key_here"

# Read receptor
with open("receptor.pdb", "r") as f:
    receptor_pdb = f.read()

# CRITICAL: Use file-like tuple, NOT json.dumps()!
receptor_file = ("receptor.pdb", BytesIO(receptor_pdb.encode('utf-8')), 'chemical/x-pdb')

multipart_data = MultipartEncoder(
    fields={
        "Input Receptor": receptor_file,
        "Input Ligand": json.dumps([{"type": "smiles", "data": "CCO"}]),
        "Number Samples": "100",  # STRING, not int!
    }
)

r = requests.post(
    "https://neurosnap.ai/api/job/submit/DiffDock-L?note=Test",
    headers={"X-API-KEY": API_KEY, "Content-Type": multipart_data.content_type},
    data=multipart_data,
)
job_id = r.json()  # Returns job ID string
```

### ADMET-AI (CORRECT FORMAT)
```python
import json, requests
from requests_toolbelt.multipart.encoder import MultipartEncoder

API_KEY = "your_key_here"

multipart_data = MultipartEncoder(
    fields={
        "Input Molecules": json.dumps([{"type": "smiles", "data": "CCO"}])
    }
)

r = requests.post(
    "https://neurosnap.ai/api/job/submit/ADMET-AI?note=Test",
    headers={"X-API-KEY": API_KEY, "Content-Type": multipart_data.content_type},
    data=multipart_data,
)
job_id = r.json()

# Wait for completion, then download results.csv (NOT output.json!)
r = requests.get(f"https://neurosnap.ai/api/job/file/{job_id}/out/results.csv",
                 headers={"X-API-KEY": API_KEY})
with open("results.csv", "wb") as f:
    f.write(r.content)
```

## Available Services (Choose Based on Your Project Needs)
- **Boltz-2 (AlphaFold3)** - Protein structure prediction
- **Conformer Generator** - Molecular conformer generation  
- **Aggrescan3D** - Protein aggregation analysis
- **ChemBERTa** - Chemical language model
- **ESM-2** - Evolutionary scale modeling
- **ProtBERT** - Protein language model
- **Molecular Transformer** - Chemical transformations

**Note**: Select the appropriate service(s) based on your specific project requirements. The examples below show tested implementations for common use cases.
### Job Endpoints

Job endpoints are HTTP endpoints related to neurosnap services and submitted jobs for those services.

#### `/api/services`

Lists all currently available Neurosnap services.

```py
r = requests.get("https://neurosnap.ai/api/services", headers={"X-API-KEY": YOUR_API_KEY})
r.json()
# [
#   {
#     "beta": false,
#     "title": "NeuroFold",
#     "desc": "The most accurate approach for in silico enzyme optimization. NeuroFold makes it trivial to optimize multiple enzyme properties, simultaneously.",
#     "desc_short": "Optimize enzyme thermostability, pH stability, solubility, and reaction rate with high accuracy.",
#     "paper": "https://doi.org/10.1101/2024.03.12.584504",
#     "video": "",
#     "tags": [
#       "Protein Design",
#       "Protein Stability",
#       "Protein Solubility"
#     ],
#     "features": [
#       "Utilizes our proprietary NeuroFold model to optimize enzymes.",
#       "40-fold increase in Spearman rank correlation when compared to Meta's ESM-1v model.",
#       "Can optimize the thermostability, pH stability, and reaction rate of a target enzyme.",
#       "Output sequences tend to be evolutionarily distinct and diverse.",
#       "Typical sequence identity to wild type is around 49-56%.",
#       "When tested on \u03b2-lactamase, NeuroFold was found to have a 100% experimental success rate.",
#       "Currently only supports monomeric enzymes that do not form complex interactions with nucleotides."
#     ],
#     "citations": [
#       "Keaun Amani, Michael Fish, Matthew D. Smith, Christian Danve M. Castroverde\nbioRxiv 2024.03.12.584504; doi: https://doi.org/10.1101/2024.03.12.584504"
#     ]
#   },
#   ...
# ]
```

#### `/api/jobs`

Lists your accounts submitted jobs.

```py
r = requests.get("https://neurosnap.ai/api/jobs", headers={"X-API-KEY": YOUR_API_KEY})
r.json()
# [
#   {
#     "Status": "completed",
#     "Job ID": "66aea7b235c83718af6161d5",
#     "Service": "PDBFixer",
#     "Submitted": 1722722226000,
#     "Runtime": "31s",
#     "Note": "test job binder + membrane",
#     "Shared": false
#   },
#   {
#     "Status": "completed",
#     "Job ID": "66ae52068d3b43e471e39c13",
#     "Service": "AntiFold",
#     "Submitted": 1722700294000,
#     "Runtime": "35s",
#     "Note": "test job using 8DXT",
#     "Shared": false
#   },
#   ...
# ]
```

#### `/api/job/submit/SERVICE_NAME?note=NOTE`

Submit a Neurosnap job. SERVICE\_NAME needs to be replaced with the exact name of the service. For example consider the following DiffDock-L job example. Note that this endpoint receives information in the format of `multipart/form-data`. Every field input field needs to be present in order for the endpoint to accept the request. Optional inputs do not need to be supplied and will use their default values. If you're unsure of what fields to provide as input, [always reference the service's submission panel page](https://neurosnap.ai/service/DiffDock-L). Job names or descriptions can be set via the optional note query parameter. If note is not provided then no job note will be set. If the job is successfully submitted this endpoint will return a job ID in JSON format.

```py
import json
import requests
from requests_toolbelt.multipart.encoder import MultipartEncoder
from io import BytesIO

# ⚠️ CRITICAL: Receptor must be sent as a file-like object, NOT JSON!
with open("receptor.pdb", "r") as f:
    receptor_pdb = f.read()

receptor_file = ("receptor.pdb", BytesIO(receptor_pdb.encode('utf-8')), 'chemical/x-pdb')

multipart_data = MultipartEncoder(
  fields={
    "Input Receptor": receptor_file,  # File-like object, NOT json.dumps()!
    "Input Ligand": json.dumps([{"data": "NC1CC(C(O)Cl)OC(N)P1", "type": "smiles"}]),
    "Number Samples": "100",  # since this field is a multiple choice field the values must be provided as a string, even if they are numbers.
  }
)

r = requests.post(
  "https://neurosnap.ai/api/job/submit/DiffDock-L",
  headers={
    "X-API-KEY": YOUR_API_KEY,
    "Content-Type": multipart_data.content_type,  # Set the correct content type
  },
  data=multipart_data,
)
r.json()
# '66affe9ca8c3d6ab6d8e377e'
```

#### `/api/job/status/JOB_ID`;

Fetches the status of a specified job. Job statuses can include pending, running, failed, and completed.

```py
r = requests.get(
  "https://neurosnap.ai/api/job/status/66affe9ca8c3d6ab6d8e377e",
  headers={"X-API-KEY": YOUR_API_KEY},
)
r.json()
# 'completed'
```

#### `/api/job/data/JOB_ID/?share=SHARE_ID` ⚠️ CORRECTED

**IMPORTANT:** This is the ONLY endpoint that lists job files. The old `/api/job/files/` endpoint does NOT exist!

Fetches all the configuration and files associated with a Neurosnap job. For shared / public jobs you can optionally provide the share ID as a query parameter as well. This is not required and job sharing is disabled by default.

```py
r = requests.get(
  "https://neurosnap.ai/api/job/data/68dd7c98c435f1b0ee471c17", # fetch output files of a job
  headers={"X-API-KEY": YOUR_API_KEY},
)
r.json()
# {
#   "config": {
#     "Density Map": "Density_Map.map"
#   },
#   "in": [
#     [
#       "Density_Map.map",
#       "64.00 MB"
#     ]
#   ],
#   "out": [
#     [
#       "output.mrc",
#       "42.30 MB"
#     ]
#   ]
# }
```

#### `/api/job/file/JOB_ID/[in/out]/FILE_NAME?share=SHARE_ID`;

Fetches a specific file to download from a completed Neurosnap job. This endpoint requires you to specify whether you want to fetch an input file (file you uploaded to run the job) or fetch an output file (file produced during the Neurosnap job by the service). For shared / public jobs you can optionally provide the share ID as a query parameter as well. This is not required and job sharing is disabled by default.

```py
r = requests.get(
  "https://neurosnap.ai/api/job/file/66affe9ca8c3d6ab6d8e377e/out/output.json", # fetch a file called output.json
  headers={"X-API-KEY": YOUR_API_KEY},
)
# save the file to disk
with open("output.json", "wb") as f:
  f.write(r.content)
```

#### `/api/job/note/set`;

Set a note for a submitted job. Notes can be thought of as mini descriptions that you can assign to your jobs for convenience and organizational purposes.

```py
r = requests.post(
  "https://neurosnap.ai/api/job/status/66affe9ca8c3d6ab6d8e377e",
  headers={"X-API-KEY": YOUR_API_KEY},
  json={"job_id": "66affe9ca8c3d6ab6d8e377e", "note": "Binder design against target PDCD1."}
)
r.raise_for_status() # no content in response but can raise if status code is not 200
```

#### `/api/job/share/set/JOB_ID`;

Enables the sharing feature of a job and makes it public to anybody with the job ID and share ID. Note that sharing is always disabled by default and that enabling it means that **anybody** with the job ID and share ID will be able to access the job results and inputs (only for jobs with sharing enabled).

```py
r = requests.get(
  "https://neurosnap.ai/api/job/share/set/66affe9ca8c3d6ab6d8e377e",
  headers={"X-API-KEY": YOUR_API_KEY},
)
r.json()
# '66affe9ca8c3d6ab6d8e377e' #the new share ID of the job
```

#### `/api/job/share/delete/JOB_ID`;

Disables the sharing feature of a job and makes the job private. Note that sharing is always disabled by default.

```py
r = requests.get(
  "https://neurosnap.ai/api/job/share/delete/66affe9ca8c3d6ab6d8e377e",
  headers={"X-API-KEY": YOUR_API_KEY},
)
r.raise_for_status() # no content in response but can raise if status code is not 200
```

### Team Endpoints

Team endpoints are HTTP endpoints related to viewing and managing [Neurosnap Teams](https://neurosnap.ai/overview?view=teams).

#### `/api/teams/info`;

Fetches your team's information if you are part of a Neurosnap Team.

```py
r = requests.get(
  "https://neurosnap.ai/api/teams/info",
  headers={"X-API-KEY": YOUR_API_KEY},
)
r.json()
# {
#   "name": "Mario Brothers Plumbing Inc.",
#   "is_leader": true,
#   "seatsMax": 2
#   "computeCreditsRemaining": 250,
#   "members": [
#     {
#       "ID": "66affe9ca8c3d6ab6d8e377e",
#       "Name": "Mario",
#       "Email": "mario@neurosnap.ai",
#       "LastLogin": 1722810058,
#       "EmailVerified": true,
#       "Leader": true
#     },
#     {
#       "ID": "64da8b62ad82b9c2bc36b143",
#       "Name": "Luigi",
#       "Email": "luigi@neurosnap.ai",
#       "LastLogin": 1722616210,
#       "EmailVerified": true,
#       "Leader": false
#     }
#   ],
# }
```

#### `/api/teams/jobs`;

Fetches all the jobs submitted by all members of your Neurosnap Team.

```py
r = requests.get(
  "https://neurosnap.ai/api/teams/jobs",
  headers={"X-API-KEY": YOUR_API_KEY},
)
r.json()
# [
#   {
#     "Status": "completed",
#     "Job ID": "66aea7b235c83718af6161d5",
#     "Service": "PDBFixer",
#     "Submitted": 1722722226000,
#     "Runtime": "31s",
#     "Note": "test job binder + membrane",
#     "Shared": false
#   },
#   {
#     "Status": "completed",
#     "Job ID": "66ae52068d3b43e471e39c13",
#     "Service": "AntiFold",
#     "Submitted": 1722700294000,
#     "Runtime": "35s",
#     "Note": "test job using 8DXT",
#     "Shared": false
#   },
#   ...
# ]
```

## Conclusion

Congrats on making it this far! You should now be equipped with everything you need to know in order to successfully and effectively utilize the Neurosnap API. Just keep in mind that we have plans to greatly expand the capabilities of our API and that changes to the existing API structure are certainly possible. Make sure to keep referencing this page for updates.

API KEY: bd0e1ed66ab2b0e73dfa1d2eba2ddf5d5aaa39d90d7d751547f81d616ddcdc30565d092bc4021535b504e6c30df5eb09a84bffc02ddafc6fb1c7abb53f123c1b


---
### ADMET-AI ⚠️ CORRECTED

**Description:** Predicts absorption, distribution, metabolism, excretion, and toxicity (ADMET) endpoints for small molecules using machine learning models. Ideal for quickly triaging prodrug ideas.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/ADMET-AI`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Molecules` | JSON string | Yes | A JSON array of molecule objects. Each object must specify a `type` (`"sdf"` or `"smiles"`) and `data` (the file content or SMILES string). | `json.dumps([{"type": "smiles", "data": "CCO"}])` |

**Output:** Returns `results.csv` (NOT `output.json`!) with columns:
- `molecule` (SMILES)
- `hERG`, `DILI`, `AMES` (toxicity probabilities)
- `CYP3A4_Veith`, `Solubility_AqSolDB`, `Caco2_Wang`, `ClinTox`

---
### DiffDock-L ⚠️ CORRECTED

**Description:** A diffusion-model based docking tool that generates ligand poses, effectively capturing induced-fit effects. Used to place prodrugs in activation enzymes or probe linker exposure.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/DiffDock-L`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Receptor` | **File-like tuple** | Yes | **CRITICAL:** Must be `(filename, BytesIO(pdb_content.encode()), 'chemical/x-pdb')`, NOT JSON! | `("receptor.pdb", BytesIO(pdb_content.encode('utf-8')), 'chemical/x-pdb')` |
| `Input Ligand` | JSON string | Yes | **CRITICAL:** Array must contain EXACTLY 1 ligand! For multiple ligands, submit separate jobs. | `json.dumps([{"type": "smiles", "data": "CCO"}])` |
| `Number Samples` | String | Yes | **Must be string**, not int! Valid values: "100", "200", "500". | `"100"` |

**⚠️ IMPORTANT LIMITATION:** DiffDock-L accepts **only 1 ligand per job**. To dock multiple ligands, you must submit separate jobs for each ligand.

---
### GNINA

**Description:** A classic docking tool that uses convolutional neural networks (CNNs) to rescore poses. Excellent for consensus rescoring after DiffDock-L to improve enrichment.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/GNINA`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Receptor` | File (.pdb) | Yes | The protein receptor structure file. | `("receptor.pdb", open("receptor.pdb", "rb"))` |
| `Input Ligand` | JSON string | Yes | A JSON array of ligand objects to dock. Each object specifies a `type` (`"sdf"` or `"smiles"`) and `data`. | `[{"type": "smiles", "data": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"}]` |

---
### FlowDock

**Description:** A geometric flow-matching framework that co-predicts protein-ligand poses and relative binding affinity. Useful as a secondary ranker for enzyme-activation hypotheses.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/FlowDock`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Receptor` | JSON string | Yes | A JSON object containing the protein sequence. | `{"aa": {"prot1": "M..."}}` |
| `Input Molecules` | JSON string | Yes | A JSON array of ligand objects to dock. Each object specifies a `type` (`"sdf"` or `"smiles"`) and `data`. | `[{"type": "smiles", "data": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"}]` |
| `Input Template` | File (.pdb) | No | An optional template structure to guide docking. | `("template.pdb", open("template.pdb", "rb"))` |
| `Number of Samples` | Integer | No | Number of docking poses to generate. | `5` |
| `Number of Steps` | Integer | No | Number of steps in the diffusion process. | `40` |
| `Sampler` | String | No | The sampling algorithm to use. | `"VDODE"` |
| `Sampler ETA` | Float | No | Sampler-specific parameter. | `1.0` |
| `Prior Type` | String | No | Method for generating the initial protein structure. | `"esmfold"` |
| `Exact Prior` | Boolean | No | If true, use the exact prior. | `false` |
| `Detect Covalent` | Boolean | No | If true, attempt to detect and model covalent bonds. | `false` |

---
### Mordred

**Description:** A cheminformatics tool that generates approximately 1800 molecular descriptors and fingerprints from a molecule's structure. Serves as a feature engine for building custom QSAR models.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/Mordred Molecular Descriptor Calculator`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Molecules` | JSON string | Yes | A JSON array of molecule objects. Each object must specify a `type` (`"sdf"` or `"smiles"`) and `data`. | `[{"type": "smiles", "data": "c1ccccc1"}]` |

---
### eTox

**Description:** Predicts toxicity and synthetic accessibility from SMILES using machine learning models. A quick sanity check to complement ADMET-AI.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/eTox Drug Toxicity Prediction`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Molecules` | JSON string | Yes | A JSON array of molecule objects. Each object must specify a `type` (`"sdf"` or `"smiles"`) and `data`. | `[{"type": "smiles", "data": "CC(=O)NC1=CC=C(O)C=C1"}]` |

---
### PDBFixer

**Description:** A tool for preparing protein structures by fixing missing atoms, residues, and formatting issues. A crucial first step for reliable downstream modeling.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/PDBFixer`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Structure` | File (.pdb) | Yes | The PDB structure file to be fixed. | `("input.pdb", open("input.pdb", "rb"))` |
| `Add Missing Residues` | Boolean | No | If true, attempts to model missing residues. | `true` |
| `Add Missing Heavy Atoms` | Boolean | No | If true, adds missing heavy atoms to existing residues. | `true` |
| `Replace Nonstandard Residues` | Boolean | No | If true, replaces non-standard residues with standard equivalents. | `true` |
| `Remove Heterogens` | String | No | Action for heterogens. Options: `"Keep Heterogens"`, `"Remove Water"`, `"Remove All"`. | `"Keep Heterogens"` |
| `Add Missing Hydrogens` | Boolean | No | If true, adds hydrogen atoms at a specified pH. | `true` |
| `pH` | Float | No | The pH to use for protonation when adding hydrogens. | `7.4` |

---
### PDB2PQR

**Description:** Assigns atomic charges and radii to a PDB structure according to a specified force field, producing a PQR file suitable for electrostatic calculations.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/PDB2PQR`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Structure` | File (.pdb) | Yes | The PDB structure file to convert. | `("input.pdb", open("input.pdb", "rb"))` |
| `Force Field` | String | No | The force field to use for parameter assignment. Options: `"PARSE"`, `"AMBER"`, `"CHARMM"`, etc. | `"AMBER"` |
| `Drop Water` | Boolean | No | If true, removes water molecules from the structure. | `true` |
| `Optimize HBN` | Boolean | No | If true, optimizes the hydrogen bond network. | `true` |
| `pH` | Float | No | The pH to use for protonation state calculations (if enabled). | `7.0` |

---
### AMBER Relaxation

**Description:** Performs a short energy minimization using the AMBER force field to validate pose plausibility and relieve steric clashes in a structure.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/AMBER Relaxation`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Structure` | File (.pdb) | Yes | The PDB structure file to relax. | `("complex.pdb", open("complex.pdb", "rb"))` |
| `Max Iterations` | Integer | No | The maximum number of minimization steps. | `2500` |
| `Tolerance` | Float | No | The energy convergence tolerance for minimization. | `1.0` |

---
### GROMACS Molecular Dynamics

**Description:** Runs a short molecular dynamics (MD) simulation to validate docked complexes, probe linker strain, and assess solvent exposure.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/GROMACS Molecular Dynamics`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Structure` | File (.pdb) | Yes | The starting structure for the MD simulation. | `("complex.pdb", open("complex.pdb", "rb"))` |
| `Simulation Duration` | Float | No | The length of the simulation in nanoseconds (ns). | `1.0` |
| `Simulation Temperature` | Integer | No | The simulation temperature in Kelvin (K). | `300` |
| `Forcefield` | String | No | The force field to use for the simulation. | `"AMBER99SB-ILDN"` |
| `Binding Energetics` | String | No | Method for calculating binding energy. | `"MMGBSA"` |
| `Solvent Box` | String | No | The shape of the solvent box. | `"Cubic"` |

---
### NetSolP-1.0

**Description:** A sequence-only predictor of protein solubility and expression "usability." Used to quickly screen variants for high-viscosity or aggregation risks.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/NetSolP-1.0`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Sequences` | JSON string | Yes | A JSON object containing one or more named protein sequences. | `{"aa": {"prot1": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"}}` |

---
### ThermoMPNN

**Description:** A deep learning model that predicts the change in stability (ΔΔG) for point mutations from a protein structure. Used to identify stabilizing mutations that improve formulation robustness.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/ThermoMPNN`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Structure` | File (.pdb) | Yes | The PDB structure file for stability analysis. | `("protein.pdb", open("protein.pdb", "rb"))` |
| `Input Chain(s)` | String | No | Comma-separated list of chains to analyze. | `"A,B"` |
| `Mode` | String | No | Analysis mode. | `"Single"` |
| `Threshold` | Float | No | ΔΔG threshold for reporting mutations. | `-0.5` |
| `Distance` | Integer | No | Distance cutoff for considering residue interactions. | `5` |
| `Include Cysteine` | Boolean | No | If true, allows mutations to cysteine. | `false` |
| `Penalize` | Boolean | No | If true, applies penalties for certain mutations. | `true` |

---
### Aggrescan3D

**Description:** Maps intrinsic aggregation propensities onto a 3D structure to identify solvent-exposed aggregation hotspots, guiding mitigation strategies.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/Aggrescan3D`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Structure` | File (.pdb) | Yes | The PDB structure file to analyze for aggregation hotspots. | `("protein.pdb", open("protein.pdb", "rb"))` |
| `Relax` | Boolean | No | If true, performs a quick energy relaxation before analysis. | `false` |

---
### TemStaPro

**Description:** A protein language model that predicts thermostability (in the 40–65 °C range) from sequence alone. A solid fallback when a reliable structure is unavailable.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/TemStaPro Protein Thermostability Prediction`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Sequences` | JSON string | Yes | A JSON object containing one or more named protein sequences. | `{"aa": {"variant1": "M..."}}` |

---
### TIsigner

**Description:** Optimizes the translation initiation site (TIS) of an mRNA sequence to improve protein expression without altering the amino acid sequence.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/TIsigner Expression Optimization`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Sequence` | String | Yes | The DNA coding sequence to be optimized. | `"ATGGCC...TAA"` |
| `Target Expression Score` | Integer | No | The desired expression score to aim for. | `70` |
| `Host Species` | String | No | The expression host. | `"Escherichia coli"` |
| `Promoter` | String | No | The promoter used in the expression system. | `"T7lac promoter (pET21)"` |
| `Restriction Sites to Avoid` | String | No | Comma-separated list of restriction sites to avoid (e.g., EcoRI, BamHI). | `"GAATTC,GGATTC"` |

---
### DockQ

**Description:** Assesses the quality of a docked protein complex by comparing it to a native (reference) structure. A lightweight QA step to catch modeling errors.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/DockQ`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Native Structure` | File (.pdb) | Yes | The reference (ground truth) complex structure. | `("native.pdb", open("native.pdb", "rb"))` |
| `Model Structure` | File (.pdb) | Yes | The predicted (docked) complex structure to evaluate. | `("model.pdb", open("model.pdb", "rb"))` |
| `Use Alpha-Carbon` | Boolean | No | If true, performs comparison using only C-alpha atoms. | `false` |

---
### USalign

**Description:** A tool for structural alignment of proteins and nucleic acids. Used to compare structures to ensure mutations did not distort the overall fold.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/USalign Structural Alignment`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `PDB File 1` | File (.pdb) | Yes | The first structure for alignment (e.g., wild-type). | `("wt.pdb", open("wt.pdb", "rb"))` |
| `PDB File 2` | File (.pdb) | Yes | The second structure for alignment (e.g., mutant). | `("mutant.pdb", open("mutant.pdb", "rb"))` |
| `TM-Score Superposition` | String | No | The algorithm for superposition. | `"Sequence independent structure alignment"` |

---
### RFantibody

**Description:** A structure-aware antibody design pipeline that proposes backbone and CDR-support edits targeted to a specific epitope, preserving antigen binding during engineering.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/RFantibody`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Framework` | File (.pdb) | Yes | The starting antibody framework structure. | `("framework.pdb", open("framework.pdb", "rb"))` |
| `Input Target` | File (.pdb) | Yes | The antigen (target) structure. | `("antigen.pdb", open("antigen.pdb", "rb"))` |
| `Hotspots` | String | Yes | Comma-separated list of target residues defining the epitope. | `"A21,A25-30"` |
| `Mode` | String | No | The type of binder to design. | `"Antibody"` |
| `Number of RFdiffusion Backbone Designs` | Integer | No | Number of backbone structures to generate. | `1` |
| `Number of ProteinMPNN Designs` | Integer | No | Number of sequences to design per backbone. | `1` |
| `Design Loops` | String | No | Specifies which CDR loops to redesign and their lengths. | `"H3:5-13,L3:10-15"` |

---
### StaB-ddG

**Description:** A deep learning model for predicting the change in binding affinity (ΔΔG) for mutations at a protein-protein interface. Ideal for scanning Fc:FcRn interface mutations.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/StaB-ddG`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Structure` | File (.pdb) | Yes | The PDB structure of the protein complex. | `("fc_fcrn.pdb", open("fc_fcrn.pdb", "rb"))` |
| `Interface Chain partner(s) 1` | String | Yes | Chains belonging to the first binding partner. | `"A"` |
| `Interface Chain partner(s) 2` | String | Yes | Chains belonging to the second binding partner. | `"B"` |
| `Mutants` | String | Yes | Comma-separated list of mutations to evaluate (e.g., WT_AA + Chain + ResID + MUT_AA). | `"YA434S,TA252Y"` |
| `Number of Monte Carlo samples` | Integer | No | Number of samples for robust ranking. | `20` |

---
### Boltz-2 (AlphaFold3-class)

**Description:** An AlphaFold3-level model for predicting complex structures (protein, nucleic acid, ligand) with an added affinity prediction head. Used to build Fc:FcRn and Fc:FcγR complexes.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/Boltz-2 (AlphaFold3)`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Sequences` | JSON string | Yes | A JSON object containing named sequences for proteins, DNA, and/or RNA. | `{"aa": {"Fc": "...", "FcRn": "..."}}` |
| `Input Molecules` | JSON string | No | A JSON array of small molecule objects (`"sdf"` or `"smiles"`). | `[{"type": "smiles", "data": "C1=CC=CC=C1"}]` |
| `MSA Mode` | String | No | The method for generating Multiple Sequence Alignments. | `"mmseqs2_uniref_env"` |
| `Number Recycles` | Integer | No | Number of refinement cycles in the structure prediction. | `6` |
| `Diffusion Samples` | Integer | No | Number of samples for the diffusion model. | `5` |
| `Diffusion Samples Affinity` | Integer | No | Number of samples for the affinity prediction head. | `5` |

---
### Prodigy Binding Affinity

**Description:** A contact-based model that predicts absolute binding affinity (ΔG and Kd) from a complex structure. Complements ΔΔG tools by providing absolute estimates.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/Prodigy Binding Affinity Prediction`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Structure` | File (.pdb) | Yes | The PDB structure of the protein-protein complex. | `("complex.pdb", open("complex.pdb", "rb"))` |
| `Temprature` | Float | No | The temperature in Celsius for Kd calculation. | `25.0` |
| `Interactor 1` | String | No | Chain ID(s) for the first interactor. | `"A"` |
| `Interactor 2` | String | No | Chain ID(s) for the second interactor. | `"B"` |

---
### Humatch

**Description:** A CNN-based tool for antibody humanization, classification, and germline pairing recommendations. Used to reduce immunogenicity risk when engineering antibodies.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/Humatch`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Sequences` | JSON string | Yes | A JSON object containing the antibody heavy and light chain sequences. | `{"aa": {"heavy": "EVQL...", "light": "DIQ..."}}` |
| `Mode` | String | No | The analysis mode. | `"Classification"` |

---
### DeepImmuno

**Description:** A CNN model that predicts peptide-MHC immunogenicity across a broad range of HLA types. Used to screen engineered protein patches for potential T-cell epitopes.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/DeepImmuno Immunogenicity Prediction`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Input Sequences` | JSON string | Yes | A JSON object containing one or more named protein sequences to screen. | `{"aa": {"Fc_variant_1": "..."}}` |

---
### NeuroBind

**Description:** An end-to-end platform for designing binders (antibodies, nanobodies, etc.) with a focus on both affinity and low immunogenicity.

**API Endpoint:** `POST https://neurosnap.ai/api/job/submit/NeuroBind`

**Headers:**
- `X-API-KEY`: Your unique API key.
- `Content-Type`: `multipart/form-data`

**Parameters:**

| Parameter Name | Type | Required? | Description | Example |
|---|---|---|---|---|
| `Target Sequences` | JSON string | Yes | A JSON object containing the target protein sequence. | `{"aa": {"prot1": "LCLYTHIGRNIYYGSYLYSETWN"}}` |
| `Design Mode` | String | No | The type of binder to design. | `"Nanobody"` |
| `Binder Length` | Integer | No | The desired length of the binder sequence. | `30` |
| `Hotspots` | String | No | Residues on the target to focus the binding interaction. | `"Protein_1:6"` |
| `Custom Template` | File (.pdb) | No | An optional structural template for the target. | `("structure.pdb", open("structure.pdb", "rb"))` |
| `Paratopes` | String | No | Residues on the binder to involve in the interaction. | `"A:6"` |