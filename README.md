# PharmacoNet Reverse Screening - Complete Setup Guide

## 📁 Repository Structure

```
your-repo/
├── .github/
│   └── workflows/
│       └── reverse_screening.yml          # GitHub Actions workflow
├── input/
│   ├── pdb_database.csv                   # Protein database (PDB codes)
│   └── query_molecules.csv                # Query molecules (SMILES)
├── scripts/
│   └── analyze_results.py                 # Results analysis script
├── batch_modeling.py                      # Batch pharmacophore modeling
├── reverse_screening.py                   # Target fishing script
├── environment.yml                        # Conda environment
├── setup.py                               # PharmacoNet installation (if available)
└── README.md                              # This file
```

## 🚀 Quick Start

### 1. Setup Repository

```bash
# Create directory structure
mkdir -p .github/workflows input scripts

# Add all files from this package
# (Copy all provided files to their respective locations)
```

### 2. Prepare Input Files

**input/pdb_database.csv:**
```csv
PDB_code,Ligand_ID,Chain
5XRA,8D3,A
6LU7,N3J,A
1ATP,ATP,
```

**input/query_molecules.csv:**
```csv
Name,SMILES
acetazolamide,CC(=O)NC1=NN=C(S1)S(=O)(=O)N
THC,CCCCCC1=CC(=C2[C@@H]3C=C(CC[C@H]3C(OC2=C1)(C)C)C)O
aspirin,CC(=O)Oc1ccccc1C(=O)O
```

### 3. Test Locally (Recommended)

```bash
# Create conda environment
conda env create -f environment.yml
conda activate pmnet

# Test batch modeling (single protein)
python batch_modeling.py \
  --input_csv input/pdb_database.csv \
  --output_dir output \
  --verbose

# Test reverse screening
python reverse_screening.py \
  --query_csv input/query_molecules.csv \
  --model_database_dir output \
  --out results.csv \
  --num_conformers 50 \
  --cpus 4
```

### 4. Deploy to GitHub Actions

```bash
# Commit and push
git add .
git commit -m "Setup PharmacoNet reverse screening workflow"
git push

# Trigger workflow manually
# Go to: Actions → PharmacoNet Reverse Screening → Run workflow
```

## 📊 Understanding Results

### Output Files

After workflow completes:

```
results/
├── screening_results.csv              # Raw scores
├── SUMMARY.md                         # Quick overview
└── analysis/
    ├── analysis_report.txt            # Statistics
    ├── score_distribution.png         # Histogram
    ├── score_by_query_boxplot.png    # Per-query comparison
    └── query_target_heatmap.png      # Heatmap
```

### Score Interpretation

| Score | Meaning | Action |
|-------|---------|--------|
| > 40  | Strong match | High priority for validation |
| 30-40 | Good match | Worth investigating |
| 20-30 | Moderate | Consider secondary validation |
| < 20  | Weak/No match | Low confidence |

## ⚙️ Configuration

### Workflow Parameters

Edit `.github/workflows/reverse_screening.yml` or use manual trigger:

- **num_conformers**: 50-100 (more = slower but better coverage)
- **min_score**: 5.0 (filter threshold)
- **top_n**: 50 (max results per query)

### Feature Weights

Edit `reverse_screening.py` defaults:

```python
--cation 8.0      # Electrostatic (strongest)
--anion 8.0       # Electrostatic (strongest)
--aromatic 4.0    # π-π stacking
--hba 4.0         # H-bond acceptor
--hbd 4.0         # H-bond donor
--halogen 4.0     # Halogen bonding
--hydrophobic 1.0 # Van der Waals (weakest)
```

## 🔧 Troubleshooting

### "No .pm files found"
- Check `input/pdb_database.csv` has valid PDB codes
- Verify network connectivity (downloads from RCSB PDB)
- Check workflow logs: Actions → Run → build_database.log

### "No screening results"
- All scores below threshold → Lower `min_score`
- Check query SMILES are valid
- Verify .pm models exist: `find output -name "*.pm"`

### "Out of memory"
- Reduce conformers: `--num_conformers 20`
- Reduce parallel processes: `--cpus 2`
- Split large databases

### Workflow timeout
- Default: 12 hours
- For >1000 proteins, consider splitting database
- Or increase timeout in workflow file

## 📚 Additional Resources

- **PharmacoNet Paper**: [Add citation]
- **PDB Database**: https://www.rcsb.org/
- **SMILES Reference**: https://pubchem.ncbi.nlm.nih.gov/

## 🆘 Support

If you encounter issues:

1. Check workflow logs (Actions tab)
2. Review troubleshooting section
3. Test locally before pushing
4. Verify input CSV format

## 📝 License

MIT License

---

**Happy Target Fishing! 🎯🧬**
