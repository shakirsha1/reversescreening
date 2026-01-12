# Complete Deployment Guide

This guide walks you through deploying PharmacoNet reverse screening from scratch.

---

## 📋 Checklist

Before starting, ensure you have:

- [ ] GitHub account
- [ ] Git installed locally
- [ ] Python 3.11+ installed (for local testing)
- [ ] Text editor (VS Code, Sublime, etc.)
- [ ] Input data ready (PDB codes and query molecules)

---

## 🚀 Deployment Steps

### Step 1: Create GitHub Repository

1. Go to https://github.com
2. Click **New repository**
3. Repository name: `pharmaconet-screening` (or your choice)
4. Privacy: Public or Private
5. ✅ Check "Add a README file"
6. Click **Create repository**

### Step 2: Clone Repository Locally

```bash
# Clone your new repository
git clone https://github.com/YOUR_USERNAME/pharmaconet-screening.git
cd pharmaconet-screening
```

### Step 3: Add All Files

**Copy all files from this package to your repository:**

```bash
# Create directory structure
mkdir -p .github/workflows input scripts

# Copy files to their locations:
# From this package → To your repository

# Workflow file
cp fixed_workflow.yml .github/workflows/reverse_screening.yml

# Python scripts
cp batch_modeling.py .
cp reverse_screening.py .
cp analyze_results.py scripts/

# Configuration files
cp environment.yml .
cp setup.py .
cp requirements.txt .
cp .gitignore .

# Documentation
cp README.md .
cp INSTALLATION.md .

# Input templates (edit these!)
cp input_pdb_database.csv input/pdb_database.csv
cp input_query_molecules.csv input/query_molecules.csv

# Quick start script
cp quickstart.sh .
chmod +x quickstart.sh
```

**Your structure should look like:**
```
pharmaconet-screening/
├── .github/
│   └── workflows/
│       └── reverse_screening.yml
├── input/
│   ├── pdb_database.csv
│   └── query_molecules.csv
├── scripts/
│   └── analyze_results.py
├── batch_modeling.py
├── reverse_screening.py
├── environment.yml
├── setup.py
├── requirements.txt
├── .gitignore
├── quickstart.sh
├── README.md
└── INSTALLATION.md
```

### Step 4: Customize Input Files

**Edit `input/pdb_database.csv`** with your proteins:

```csv
PDB_code,Ligand_ID,Chain
5XRA,8D3,A
6LU7,N3J,A
1ATP,ATP,
2BRC,BRC,B
3COX,ASA,A
```

**Tips for finding PDB codes:**
- Visit https://www.rcsb.org/
- Search by protein name (e.g., "kinase", "protease")
- Filter by: "Has Ligand" = Yes, Resolution < 3.0 Å
- Download search results as CSV
- Extract PDB codes

**Edit `input/query_molecules.csv`** with your molecules:

```csv
Name,SMILES
acetazolamide,CC(=O)NC1=NN=C(S1)S(=O)(=O)N
aspirin,CC(=O)Oc1ccccc1C(=O)O
ibuprofen,CC(C)Cc1ccc(cc1)C(C)C(=O)O
```

**Tips for getting SMILES:**
- PubChem: https://pubchem.ncbi.nlm.nih.gov/
- Search by compound name → Copy "Canonical SMILES"
- DrugBank: https://go.drugbank.com/
- ChEMBL: https://www.ebi.ac.uk/chembl/

### Step 5: Test Locally (Optional but Recommended)

```bash
# Install conda/mamba if needed
# See INSTALLATION.md for details

# Run quick start
./quickstart.sh

# This will:
# 1. Create environment
# 2. Build database (small test)
# 3. Run screening
# 4. Generate results
```

**If local test succeeds → proceed to GitHub Actions**
**If local test fails → fix issues before deploying**

### Step 6: Commit and Push to GitHub

```bash
# Check what will be committed
git status

# Add all files
git add .

# Commit
git commit -m "Setup PharmacoNet reverse screening workflow"

# Push to GitHub
git push origin main
```

**Verify on GitHub:**
1. Go to your repository URL
2. Check all files are uploaded
3. Verify `.github/workflows/reverse_screening.yml` exists

### Step 7: Enable GitHub Actions

1. Go to your repository on GitHub
2. Click **Actions** tab
3. If prompted, click **I understand my workflows, go ahead and enable them**

### Step 8: Run Workflow

**Option A: Manual Trigger (Recommended for first run)**

1. Click **Actions** tab
2. Select **PharmacoNet Reverse Screening** from left sidebar
3. Click **Run workflow** button (right side)
4. Configure parameters (or use defaults):
   - Number of conformers: `50`
   - Minimum score: `5.0`
   - Top N results: `50`
5. Click green **Run workflow** button
6. Wait for workflow to start

**Option B: Auto-trigger (Push to main)**

Workflow automatically runs when you push changes to:
- `input/pdb_database.csv`
- `input/query_molecules.csv`
- `.github/workflows/reverse_screening.yml`

### Step 9: Monitor Progress

1. Click on the running workflow (blue dot = running)
2. Click on **reverse-screening** job
3. Expand steps to see progress:
   - **Build protein pharmacophore database** ← Takes longest
   - **Run reverse screening** ← Quick if database exists
   - **Analyze results**

**Estimated times:**
- 10 proteins: ~5-10 minutes
- 100 proteins: ~1-2 hours
- 1000 proteins: ~8-12 hours

**Logs you'll see:**
```
Processing 1/10: 5XRA
  → Success: 5XRA_A_8D3_model.pm
Processing 2/10: 6LU7
  → Success: 6LU7_A_N3J_model.pm
...
```

### Step 10: Download Results

Once completed (✅ green check):

1. Scroll to bottom of workflow page
2. **Artifacts** section shows:
   - 📦 `screening-results-XXXXXX` ← **Download this**
   - 📦 `pharmacophore-models-XXXXXX` (optional)
   - 📦 `logs-XXXXXX` (for debugging)

3. Click artifact name to download ZIP
4. Extract ZIP file
5. Open `results/screening_results.csv`

**Results structure:**
```
screening-results-XXXXXX/
├── screening_results.csv          ← Main results
├── SUMMARY.md                      ← Quick overview
└── analysis/
    ├── analysis_report.txt         ← Statistics
    ├── score_distribution.png      ← Histogram
    ├── score_by_query_boxplot.png ← Per-query
    └── query_target_heatmap.png   ← Heatmap
```

---

## 🔧 Troubleshooting Deployment

### Issue: Workflow doesn't appear in Actions tab

**Solution:**
- Verify `.github/workflows/reverse_screening.yml` exists
- Check YAML syntax (use online YAML validator)
- Enable Actions if disabled (Actions → Enable)

### Issue: Workflow fails immediately

**Most common causes:**

1. **Invalid YAML syntax**
   - Use YAML validator: http://www.yamllint.com/
   - Check indentation (use spaces, not tabs)

2. **Empty input files**
   ```bash
   # Check locally
   wc -l input/pdb_database.csv
   # Should show > 1 (header + data)
   ```

3. **Missing files**
   - Ensure `batch_modeling.py` exists
   - Ensure `reverse_screening.py` exists
   - Ensure `scripts/analyze_results.py` exists

### Issue: Database building step fails

**Check workflow logs:**
```
Actions → Your workflow → reverse-screening job → Build protein pharmacophore database
```

**Common causes:**

1. **Invalid PDB codes**
   ```csv
   # Bad:
   PDB_code,Ligand_ID,Chain
   INVALID,XXX,A

   # Good:
   PDB_code,Ligand_ID,Chain
   5XRA,8D3,A
   ```

2. **Network issues** (GitHub runners can't reach RCSB PDB)
   - Retry workflow
   - Check RCSB PDB status: https://www.rcsb.org/

3. **Python package issues**
   - Check "Install Python dependencies" step for errors
   - Verify all packages installed

### Issue: "No .pm files found in database directory"

**Fix:**
1. Review `build_database.log` in artifacts
2. Check for specific PDB download errors
3. Remove problematic PDB codes from CSV
4. Re-run workflow

### Issue: No results / all scores below threshold

**Solutions:**

1. **Lower threshold:**
   - Re-run with `min_score: 0.0` to see all scores

2. **Check SMILES validity:**
   ```python
   from rdkit import Chem
   mol = Chem.MolFromSmiles("YOUR_SMILES")
   print(mol)  # Should not be None
   ```

3. **Increase conformers:**
   - Re-run with `num_conformers: 100`

### Issue: Workflow timeout (> 12 hours)

**Solutions:**

1. **Split database:**
   ```bash
   # Create pdb_database_part1.csv (first 500)
   # Create pdb_database_part2.csv (next 500)
   # Run workflow twice
   ```

2. **Increase timeout in workflow:**
   ```yaml
   # In .github/workflows/reverse_screening.yml
   timeout-minutes: 1440  # 24 hours
   ```

---

## 📊 Interpreting Results

### Quick Results Check

**Look at `SUMMARY.md`:**
```markdown
## Input Summary
- Proteins in database: 100
- Pharmacophore models generated: 98  ← Should be close to input
- Query molecules: 5

## Top 5 Matches
query_name         pharmacophore_model              score
acetazolamide      output/5XRA/5XRA_A_8D3_model.pm  42.1
acetazolamide      output/6LU7/6LU7_A_N3J_model.pm  35.4
```

### Full Results Analysis

**Open `results/analysis/analysis_report.txt`:**
```
SUMMARY STATISTICS
Total matches: 250
Unique queries: 5
Unique targets: 98

Score Statistics:
  Mean:   15.2341
  Median: 12.4567
  Max:    42.1000

TOP 10 QUERY-TARGET PAIRS
acetazolamide → 5XRA_A_8D3_model    Score: 42.1000
acetazolamide → 6LU7_A_N3J_model    Score: 35.4000
...
```

### Score Meaning

| Score | Interpretation | Next Steps |
|-------|---------------|------------|
| > 40 | **Strong match** | High priority for experimental validation |
| 30-40 | **Good match** | Worth investigating further |
| 20-30 | **Moderate** | Consider molecular docking first |
| 10-20 | **Weak** | Low confidence, secondary priority |
| < 10 | **Unlikely** | Probably not a real interaction |

---

## 🎯 Next Steps After Results

### 1. Validate Top Hits

**For scores > 30:**

```bash
# Run molecular docking (e.g., AutoDock Vina)
# This gives binding pose + affinity estimate

# Example:
autodock_vina \
  --receptor protein.pdbqt \
  --ligand molecule.pdbqt \
  --out docking_result.pdbqt
```

### 2. Literature Search

```bash
# Check if interactions are known
# Search PubMed, Google Scholar for:
# "[Query molecule] [Target protein] binding"
```

### 3. Experimental Validation

**In vitro binding assays:**
- Surface Plasmon Resonance (SPR)
- Isothermal Titration Calorimetry (ITC)
- Microscale Thermophoresis (MST)

**Functional assays:**
- Enzyme inhibition (IC50)
- Cell-based assays
- Phenotypic screens

### 4. Structural Confirmation

**Get binding mode:**
- X-ray crystallography
- Cryo-EM
- NMR spectroscopy

---

## 📚 Additional Resources

- **PharmacoNet Documentation**: [Add link]
- **RCSB PDB Help**: https://www.rcsb.org/docs/
- **RDKit Documentation**: https://www.rdkit.org/docs/
- **GitHub Actions Docs**: https://docs.github.com/en/actions

---

## ✅ Deployment Checklist

- [ ] GitHub repository created
- [ ] All files added and pushed
- [ ] Input CSVs customized with your data
- [ ] Tested locally (optional)
- [ ] Workflow runs successfully
- [ ] Results downloaded
- [ ] Top hits identified
- [ ] Ready for validation

**Congratulations! You've successfully deployed PharmacoNet reverse screening! 🎉**

---

## 🆘 Need Help?

If you encounter issues not covered here:

1. Check workflow logs carefully
2. Review INSTALLATION.md for environment issues
3. Test problematic steps locally
4. Check GitHub Actions documentation
5. Verify input file formats

**Happy target fishing! 🎯🧬**
