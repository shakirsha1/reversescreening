# Complete File Structure & Deployment Package

## 📦 Package Contents

This package contains everything you need to set up PharmacoNet reverse screening on GitHub Actions.

### Core Files

| File | Purpose | Location in Your Repo |
|------|---------|----------------------|
| **fixed_workflow.yml** | GitHub Actions workflow | `.github/workflows/reverse_screening.yml` |
| **batch_modeling.py** | Build pharmacophore database | Root directory |
| **reverse_screening.py** | Target fishing script | Root directory |
| **analyze_results.py** | Results analysis | `scripts/analyze_results.py` |
| **environment.yml** | Conda environment | Root directory |
| **setup.py** | Python package setup | Root directory |
| **requirements.txt** | Pip dependencies | Root directory |
| **.gitignore** | Git ignore patterns | Root directory |
| **quickstart.sh** | Quick start script | Root directory |

### Input Templates

| File | Purpose | Action Required |
|------|---------|----------------|
| **input_pdb_database.csv** | Protein structures | ✏️ Edit with your PDB codes |
| **input_query_molecules.csv** | Query molecules | ✏️ Edit with your SMILES |

### Documentation

| File | Contents |
|------|----------|
| **README.md** | Overview and quick reference |
| **INSTALLATION.md** | Detailed installation guide |
| **DEPLOYMENT.md** | Step-by-step deployment |
| **FILE_STRUCTURE.md** | This file |

---

## 🗂️ Target Repository Structure

After deployment, your repository should look like this:

```
your-pharmaconet-repo/
│
├── .github/
│   └── workflows/
│       └── reverse_screening.yml        ← Rename from fixed_workflow.yml
│
├── input/
│   ├── pdb_database.csv                 ← From input_pdb_database.csv
│   └── query_molecules.csv              ← From input_query_molecules.csv
│
├── scripts/
│   └── analyze_results.py               ← Copy from package
│
├── batch_modeling.py                     ← Copy from package
├── reverse_screening.py                  ← Copy from package (or use existing)
├── environment.yml                       ← Copy from package
├── setup.py                              ← Copy from package
├── requirements.txt                      ← Copy from package
├── .gitignore                           ← Copy from package
├── quickstart.sh                        ← Copy from package
│
├── README.md                            ← Copy from package
├── INSTALLATION.md                      ← Copy from package
└── DEPLOYMENT.md                        ← Copy from package
```

---

## 📋 Quick Deployment Checklist

### Phase 1: Setup (5 minutes)

- [ ] Create GitHub repository
- [ ] Clone repository locally
- [ ] Create directory structure:
  ```bash
  mkdir -p .github/workflows input scripts
  ```

### Phase 2: Add Files (10 minutes)

- [ ] Copy `fixed_workflow.yml` → `.github/workflows/reverse_screening.yml`
- [ ] Copy `batch_modeling.py` → root
- [ ] Copy `reverse_screening.py` → root (if not already present)
- [ ] Copy `analyze_results.py` → `scripts/`
- [ ] Copy `environment.yml` → root
- [ ] Copy `setup.py` → root
- [ ] Copy `requirements.txt` → root
- [ ] Copy `.gitignore` → root
- [ ] Copy `quickstart.sh` → root
- [ ] Copy all `.md` files → root

### Phase 3: Customize Inputs (15 minutes)

- [ ] Edit `input/pdb_database.csv`:
  - Add your PDB codes
  - Verify format: `PDB_code,Ligand_ID,Chain`
  - Minimum 1 protein, recommended 10+

- [ ] Edit `input/query_molecules.csv`:
  - Add your query molecules
  - Verify format: `Name,SMILES`
  - Validate SMILES using RDKit

### Phase 4: Test Locally (Optional, 30 minutes)

- [ ] Install conda/mamba
- [ ] Run `./quickstart.sh`
- [ ] Verify results generated
- [ ] Fix any errors before pushing

### Phase 5: Deploy to GitHub (5 minutes)

- [ ] Commit all files:
  ```bash
  git add .
  git commit -m "Setup PharmacoNet reverse screening"
  git push
  ```

- [ ] Enable GitHub Actions
- [ ] Verify workflow appears in Actions tab

### Phase 6: Run Workflow (Variable time)

- [ ] Trigger workflow manually
- [ ] Monitor progress
- [ ] Wait for completion
  - 10 proteins: ~10 minutes
  - 100 proteins: ~1-2 hours
  - 1000 proteins: ~8-12 hours

### Phase 7: Download & Analyze (10 minutes)

- [ ] Download artifacts from GitHub Actions
- [ ] Extract ZIP file
- [ ] Review `screening_results.csv`
- [ ] Check `analysis/` plots
- [ ] Read `SUMMARY.md`

---

## 🔧 File-by-File Guide

### 1. fixed_workflow.yml
**What it does:** GitHub Actions workflow definition
**Rename to:** `.github/workflows/reverse_screening.yml`
**Edit needed:** ❌ No (works as-is)
**Key features:**
- Auto-triggers on CSV changes
- Manual trigger with parameters
- Parallel processing
- Artifact uploads

### 2. batch_modeling.py
**What it does:** Builds pharmacophore models from PDB structures
**Location:** Root directory
**Edit needed:** ❌ No
**Usage:**
```bash
python batch_modeling.py \
  --input_csv input/pdb_database.csv \
  --output_dir output \
  --cuda  # Optional GPU flag
```

### 3. reverse_screening.py
**What it does:** Screens query molecules against protein database
**Location:** Root directory
**Edit needed:** ❌ No (unless customizing weights)
**Usage:**
```bash
python reverse_screening.py \
  --query_csv input/query_molecules.csv \
  --model_database_dir output \
  --out results.csv \
  --cpus 8
```

### 4. analyze_results.py
**What it does:** Generates plots and statistics
**Location:** `scripts/analyze_results.py`
**Edit needed:** ❌ No
**Usage:**
```bash
python scripts/analyze_results.py \
  results.csv \
  --output_dir analysis/
```

### 5. environment.yml
**What it does:** Defines conda environment
**Location:** Root directory
**Edit needed:** ❌ No (unless adding packages)
**Contents:**
- Python 3.11
- PyTorch
- RDKit
- Scientific packages

### 6. input_pdb_database.csv
**What it does:** Template for protein database
**Rename to:** `input/pdb_database.csv`
**Edit needed:** ✅ YES - Add your PDB codes
**Format:**
```csv
PDB_code,Ligand_ID,Chain
5XRA,8D3,A
6LU7,N3J,A
```

**How to populate:**
1. Visit https://www.rcsb.org/
2. Search for proteins (e.g., "kinase")
3. Filter: "Has Ligand" = Yes
4. Download PDB IDs
5. Add to CSV

### 7. input_query_molecules.csv
**What it does:** Template for query molecules
**Rename to:** `input/query_molecules.csv`
**Edit needed:** ✅ YES - Add your molecules
**Format:**
```csv
Name,SMILES
aspirin,CC(=O)Oc1ccccc1C(=O)O
```

**How to get SMILES:**
1. Go to PubChem: https://pubchem.ncbi.nlm.nih.gov/
2. Search compound name
3. Copy "Canonical SMILES"
4. Add to CSV

### 8. quickstart.sh
**What it does:** Interactive setup and run script
**Location:** Root directory
**Edit needed:** ❌ No
**Usage:**
```bash
chmod +x quickstart.sh
./quickstart.sh
```

**What it does:**
1. Checks inputs
2. Creates environment
3. Builds database
4. Runs screening
5. Analyzes results

---

## 🎯 Common File Operations

### Copying Files to Your Repository

**Option 1: Manual Copy**
```bash
# From download location to your repo
cp ~/Downloads/fixed_workflow.yml your-repo/.github/workflows/reverse_screening.yml
cp ~/Downloads/batch_modeling.py your-repo/
cp ~/Downloads/analyze_results.py your-repo/scripts/
# ... etc
```

**Option 2: Using Script**
```bash
# Create a deployment script
cat > deploy_files.sh << 'EOF'
#!/bin/bash
PACKAGE_DIR="~/Downloads/pharmaconet-package"
REPO_DIR="./your-repo"

# Core files
cp $PACKAGE_DIR/fixed_workflow.yml $REPO_DIR/.github/workflows/reverse_screening.yml
cp $PACKAGE_DIR/batch_modeling.py $REPO_DIR/
cp $PACKAGE_DIR/analyze_results.py $REPO_DIR/scripts/
cp $PACKAGE_DIR/environment.yml $REPO_DIR/
cp $PACKAGE_DIR/setup.py $REPO_DIR/
cp $PACKAGE_DIR/requirements.txt $REPO_DIR/
cp $PACKAGE_DIR/.gitignore $REPO_DIR/
cp $PACKAGE_DIR/quickstart.sh $REPO_DIR/

# Documentation
cp $PACKAGE_DIR/*.md $REPO_DIR/

# Input templates
cp $PACKAGE_DIR/input_pdb_database.csv $REPO_DIR/input/pdb_database.csv
cp $PACKAGE_DIR/input_query_molecules.csv $REPO_DIR/input/query_molecules.csv

echo "✓ All files copied!"
EOF

chmod +x deploy_files.sh
./deploy_files.sh
```

### Verifying File Structure

```bash
# Check all required files are present
cd your-repo

echo "Checking file structure..."

[ -f .github/workflows/reverse_screening.yml ] && echo "✓ Workflow" || echo "✗ Missing workflow"
[ -f batch_modeling.py ] && echo "✓ batch_modeling.py" || echo "✗ Missing batch_modeling.py"
[ -f scripts/analyze_results.py ] && echo "✓ analyze_results.py" || echo "✗ Missing analyze_results.py"
[ -f input/pdb_database.csv ] && echo "✓ pdb_database.csv" || echo "✗ Missing pdb_database.csv"
[ -f input/query_molecules.csv ] && echo "✓ query_molecules.csv" || echo "✗ Missing query_molecules.csv"
[ -f environment.yml ] && echo "✓ environment.yml" || echo "✗ Missing environment.yml"

echo "Done!"
```

---

## 📊 File Size Reference

| File | Typical Size |
|------|--------------|
| fixed_workflow.yml | ~10 KB |
| batch_modeling.py | ~6 KB |
| reverse_screening.py | ~6 KB (yours may vary) |
| analyze_results.py | ~6 KB |
| environment.yml | ~0.3 KB |
| setup.py | ~1 KB |
| requirements.txt | ~0.2 KB |
| .gitignore | ~0.3 KB |
| quickstart.sh | ~5 KB |
| README.md | ~4 KB |
| INSTALLATION.md | ~7 KB |
| DEPLOYMENT.md | ~11 KB |
| input_pdb_database.csv | Varies (your data) |
| input_query_molecules.csv | Varies (your data) |

**Total package size:** ~50 KB (excluding your input data)

---

## ✅ Pre-deployment Verification

Run this checklist before pushing to GitHub:

```bash
# 1. Check all files exist
find . -name "*.py" -o -name "*.yml" -o -name "*.csv" -o -name "*.md"

# 2. Verify Python syntax
python -m py_compile batch_modeling.py
python -m py_compile scripts/analyze_results.py

# 3. Validate YAML
# Install: pip install yamllint
yamllint .github/workflows/reverse_screening.yml

# 4. Check CSV format
head -3 input/pdb_database.csv
head -3 input/query_molecules.csv

# 5. Test imports (if pmnet available)
python -c "from rdkit import Chem; import pandas; import torch"

# If all pass → Ready to deploy!
```

---

## 🆘 Troubleshooting File Issues

### Issue: "File not found" errors

**Check:**
```bash
# List all files recursively
find . -type f -not -path '*/\.*'

# Compare with expected structure above
```

### Issue: "Permission denied" when running quickstart.sh

**Fix:**
```bash
chmod +x quickstart.sh
```

### Issue: YAML syntax errors

**Fix:**
```bash
# Use online validator
# Copy .github/workflows/reverse_screening.yml contents
# Paste into: http://www.yamllint.com/
# Fix any reported errors
```

### Issue: CSV not recognized by workflow

**Check format:**
```bash
# Must have Unix line endings (LF not CRLF)
file input/pdb_database.csv

# Should show: ASCII text
# If shows: ASCII text, with CRLF line terminators

# Fix on Linux/Mac:
dos2unix input/pdb_database.csv

# Fix on Windows (Git Bash):
sed -i 's/\r$//' input/pdb_database.csv
```

---

## 📚 Additional Information

### File Dependencies

```
reverse_screening.yml (workflow)
├─ Requires: batch_modeling.py
├─ Requires: reverse_screening.py
├─ Requires: scripts/analyze_results.py
├─ Requires: input/pdb_database.csv
├─ Requires: input/query_molecules.csv
└─ Optional: environment.yml, setup.py

batch_modeling.py
└─ Uses: pmnet module (from setup.py)

reverse_screening.py
└─ Uses: pmnet module (from setup.py)

analyze_results.py
└─ Requires: pandas, matplotlib, seaborn
```

### Output Files (Generated by Workflow)

```
output/                          ← Pharmacophore models (.pm files)
├── 5XRA/
│   ├── 5XRA.pdb
│   ├── 5XRA_A_8D3.sdf
│   └── 5XRA_A_8D3_model.pm
├── 6LU7/
│   └── ...
└── ...

results/                         ← Screening results
├── screening_results.csv
├── SUMMARY.md
└── analysis/
    ├── analysis_report.txt
    ├── score_distribution.png
    ├── score_by_query_boxplot.png
    └── query_target_heatmap.png

*.log                           ← Log files
├── build_database.log
├── screening.log
└── analysis.log
```

---

## 🎓 Learning Resources

### Understanding Each File

1. **Workflow files (.yml)**
   - Tutorial: https://docs.github.com/en/actions/learn-github-actions
   - YAML syntax: https://yaml.org/

2. **Python scripts (.py)**
   - PharmacoNet docs: [Add link]
   - RDKit tutorial: https://www.rdkit.org/docs/GettingStartedInPython.html

3. **CSV files (.csv)**
   - RFC 4180 standard: https://tools.ietf.org/html/rfc4180
   - Pandas CSV guide: https://pandas.pydata.org/docs/user_guide/io.html

---

**You now have all the files and knowledge to deploy PharmacoNet reverse screening! 🚀**

For detailed instructions, see:
- Quick start: README.md
- Local setup: INSTALLATION.md
- GitHub deployment: DEPLOYMENT.md
