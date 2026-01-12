# Installation Guide

## Prerequisites

- **Operating System**: Linux, macOS, or Windows
- **Python**: 3.11 or higher
- **Git**: For cloning repository
- **GitHub Account**: For running GitHub Actions

## Method 1: GitHub Actions (Recommended for Large-Scale)

### Step 1: Create Repository Structure

```bash
# Create new repository or use existing
mkdir pharmaconet-screening
cd pharmaconet-screening
git init

# Create directory structure
mkdir -p .github/workflows input scripts
```

### Step 2: Add All Files

Copy the following files to your repository:

```
pharmaconet-screening/
├── .github/workflows/reverse_screening.yml
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
└── README.md
```

### Step 3: Configure Input Files

**Edit `input/pdb_database.csv`:**
```csv
PDB_code,Ligand_ID,Chain
5XRA,8D3,A
6LU7,N3J,A
1ATP,ATP,
```

**Edit `input/query_molecules.csv`:**
```csv
Name,SMILES
acetazolamide,CC(=O)NC1=NN=C(S1)S(=O)(=O)N
THC,CCCCCC1=CC(=C2[C@@H]3C=C(CC[C@H]3C(OC2=C1)(C)C)C)O
```

### Step 4: Push to GitHub

```bash
git add .
git commit -m "Initial setup for PharmacoNet reverse screening"
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPO.git
git branch -M main
git push -u origin main
```

### Step 5: Run Workflow

1. Go to your repository on GitHub
2. Click **Actions** tab
3. Select **PharmacoNet Reverse Screening**
4. Click **Run workflow**
5. Optionally configure parameters:
   - Number of conformers (default: 50)
   - Minimum score (default: 5.0)
   - Top N results (default: 50)
6. Click **Run workflow** button

### Step 6: Download Results

Once completed:
1. Go to Actions → Your workflow run
2. Scroll to **Artifacts** section
3. Download:
   - `screening-results-XXXXXX` (main results)
   - `pharmacophore-models-XXXXXX` (optional, for reuse)
   - `logs-XXXXXX` (for debugging)

---

## Method 2: Local Installation

### Step 1: Install Conda/Mamba (if not already installed)

```bash
# Download Miniforge (includes mamba)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh

# Or use existing conda
conda install mamba -c conda-forge
```

### Step 2: Create Environment

```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/YOUR_REPO.git
cd YOUR_REPO

# Create conda environment
conda env create -f environment.yml

# OR use pip
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt

# Install PharmacoNet package
pip install -e .
```

### Step 3: Verify Installation

```bash
# Activate environment
conda activate pmnet  # or: source venv/bin/activate

# Test imports
python -c "import torch; print(f'PyTorch: {torch.__version__}')"
python -c "import rdkit; print(f'RDKit: {rdkit.__version__}')"
python -c "import pmnet; print(f'PharmacoNet: {pmnet.__version__}')"
```

### Step 4: Build Protein Database

```bash
# Prepare input CSV
# Edit input/pdb_database.csv with your PDB codes

# Run batch modeling
python batch_modeling.py \
  --input_csv input/pdb_database.csv \
  --output_dir output \
  --cuda  # Remove if no GPU available

# Check output
find output -name "*.pm" | wc -l
```

Expected output:
```
Processing 1/10: 5XRA
  → Success: 5XRA_A_8D3_model.pm
Processing 2/10: 6LU7
  → Success: 6LU7_A_N3J_model.pm
...
Generated 10 pharmacophore models
```

### Step 5: Run Reverse Screening

```bash
# Edit input/query_molecules.csv with your molecules

# Run screening
python reverse_screening.py \
  --query_csv input/query_molecules.csv \
  --model_database_dir output \
  --out results/screening_results.csv \
  --num_conformers 50 \
  --min_score 5.0 \
  --top_n 50 \
  --cpus 8  # Adjust based on your CPU cores
```

Expected output:
```
Loading queries from CSV: input/query_molecules.csv
Found 2 query molecules
Found 10 pharmacophore models in database

Processing query 1/2: acetazolamide
Screening against 10 protein models...

TOP 10 MATCHING PROTEINS for acetazolamide:
  1. 5XRA_A_8D3_model                       Score: 35.4000
  2. 6LU7_A_N3J_model                       Score: 28.7000
  ...

Reverse screening complete! ✓
Results saved to: results/screening_results.csv
```

### Step 6: Analyze Results

```bash
# Run analysis script
python scripts/analyze_results.py \
  results/screening_results.csv \
  --output_dir results/analysis

# View results
cat results/analysis/analysis_report.txt
```

---

## Troubleshooting

### Issue: "ModuleNotFoundError: No module named 'pmnet'"

**Solution:**
```bash
# Install in development mode
pip install -e .

# OR skip if using standalone scripts
# (batch_modeling.py and reverse_screening.py can work without pmnet package)
```

### Issue: "CUDA out of memory"

**Solution:**
```bash
# Use CPU mode (slower but works)
python batch_modeling.py --input_csv input/pdb_database.csv --output_dir output
# (Remove --cuda flag)
```

### Issue: "No .pm files generated"

**Possible causes:**
1. Invalid PDB codes in CSV
2. Network issues (cannot download from RCSB PDB)
3. Python environment issues

**Debug:**
```bash
# Test single protein manually
python -c "
from pmnet import modeling
modeling.main(['-p', '5XRA', '-l', '8D3', '-c', 'A', '-o', 'test_output'])
"

# Check if file was created
ls test_output/5XRA/
```

### Issue: "Low/no scores for all proteins"

**Solution:**
1. Check SMILES validity:
```python
from rdkit import Chem
mol = Chem.MolFromSmiles("YOUR_SMILES_HERE")
print(mol)  # Should not be None
```

2. Increase conformers:
```bash
python reverse_screening.py --num_conformers 100 ...
```

3. Lower score threshold:
```bash
python reverse_screening.py --min_score 0.0 ...
```

---

## Performance Optimization

### For Large Databases (>1000 proteins)

**1. Use GPU for modeling:**
```bash
python batch_modeling.py --cuda ...
```

**2. Parallelize screening:**
```bash
python reverse_screening.py --cpus $(nproc) ...
```

**3. Split database:**
```bash
# Split CSV into chunks
split -l 500 input/pdb_database.csv pdb_chunk_

# Process each chunk
for chunk in pdb_chunk_*; do
  python batch_modeling.py --input_csv $chunk --output_dir output_$chunk
done

# Merge outputs
mkdir -p output
cp -r output_*/* output/
```

### Memory Management

**Reduce memory usage:**
```bash
# Fewer parallel processes
python reverse_screening.py --cpus 4 ...

# Fewer conformers
python reverse_screening.py --num_conformers 20 ...
```

---

## Next Steps

1. ✅ Install environment
2. ✅ Build protein database
3. ✅ Run reverse screening
4. ✅ Analyze results
5. 🔬 Validate top hits with molecular docking
6. 🧪 Perform experimental validation

**Good luck! 🎯**
