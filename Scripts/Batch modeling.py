#!/usr/bin/env python3
"""
Batch Pharmacophore Modeling for PharmacoNet
============================================

Process multiple PDB structures to build a pharmacophore model database.

Author: PharmacoNet Team
License: MIT
"""

import argparse
import csv
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Batch build pharmacophore models from PDB database",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input_csv",
        type=str,
        required=True,
        help="CSV file with columns: PDB_code,Ligand_ID,Chain",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="output directory for pharmacophore models",
    )
    parser.add_argument(
        "--cuda", action="store_true", help="use GPU acceleration"
    )
    parser.add_argument(
        "--force", action="store_true", help="overwrite existing models"
    )
    parser.add_argument(
        "--verbose", action="store_true", help="verbose output"
    )

    args = parser.parse_args()

    # Read input CSV
    csv_path = Path(args.input_csv)
    if not csv_path.exists():
        print(f"ERROR: Input CSV not found: {args.input_csv}")
        sys.exit(1)

    print(f"Reading input from: {args.input_csv}")
    
    entries = []
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pdb_code = row.get("PDB_code", "").strip()
            ligand_id = row.get("Ligand_ID", "").strip()
            chain = row.get("Chain", "").strip()
            
            if not pdb_code:
                continue
                
            entries.append({
                "pdb_code": pdb_code,
                "ligand_id": ligand_id if ligand_id else None,
                "chain": chain if chain else None,
            })

    if not entries:
        print("ERROR: No valid entries found in CSV!")
        print("Expected format:")
        print("PDB_code,Ligand_ID,Chain")
        print("5XRA,8D3,A")
        sys.exit(1)

    print(f"Found {len(entries)} protein structures to process")

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Import pmnet (after argument parsing to fail fast if CSV is wrong)
    try:
        from pmnet import modeling
    except ImportError as e:
        print(f"ERROR: Cannot import pmnet. Is it installed?")
        print(f"Run: pip install -e .")
        print(f"Error: {e}")
        sys.exit(1)

    # Process each entry
    success_count = 0
    error_count = 0
    skipped_count = 0

    for idx, entry in enumerate(entries, start=1):
        pdb_code = entry["pdb_code"]
        ligand_id = entry["ligand_id"]
        chain = entry["chain"]

        print(f"\n[{idx}/{len(entries)}] Processing {pdb_code}...")

        # Check if model already exists
        pdb_output_dir = output_dir / pdb_code
        existing_models = list(pdb_output_dir.glob("*.pm")) if pdb_output_dir.exists() else []
        
        if existing_models and not args.force:
            print(f"  → Skipped (model already exists: {existing_models[0].name})")
            skipped_count += 1
            continue

        # Build command
        cmd_args = [
            "modeling.py" if args.verbose else None,  # Script name for display
            "-p", pdb_code,
            "-o", str(output_dir),
        ]
        
        if ligand_id:
            cmd_args.extend(["-l", ligand_id])
        if chain:
            cmd_args.extend(["-c", chain])
        if args.cuda:
            cmd_args.append("--cuda")

        # Run modeling
        try:
            # Call the modeling module directly
            import sys
            from io import StringIO
            
            # Capture output
            old_stdout = sys.stdout
            old_stderr = sys.stderr
            
            if not args.verbose:
                sys.stdout = StringIO()
                sys.stderr = StringIO()
            
            # Build argument list for modeling
            modeling_args = [
                "-p", pdb_code,
                "-o", str(output_dir),
            ]
            if ligand_id:
                modeling_args.extend(["-l", ligand_id])
            if chain:
                modeling_args.extend(["-c", chain])
            if args.cuda:
                modeling_args.append("--cuda")
            
            # Call modeling
            modeling.main(modeling_args)
            
            # Restore output
            sys.stdout = old_stdout
            sys.stderr = old_stderr
            
            # Check if model was created
            pdb_output_dir = output_dir / pdb_code
            models = list(pdb_output_dir.glob("*.pm")) if pdb_output_dir.exists() else []
            
            if models:
                print(f"  → Success: {models[0].name}")
                success_count += 1
            else:
                print(f"  → Failed: No model file created")
                error_count += 1
                
        except Exception as e:
            sys.stdout = old_stdout
            sys.stderr = old_stderr
            print(f"  → Error: {str(e)}")
            error_count += 1

    # Summary
    print(f"\n{'='*80}")
    print(f"BATCH MODELING COMPLETE")
    print(f"{'='*80}")
    print(f"Total entries: {len(entries)}")
    print(f"Success: {success_count}")
    print(f"Skipped (already exist): {skipped_count}")
    print(f"Errors: {error_count}")
    print(f"{'='*80}")

    if error_count > 0:
        print(f"\nWARNING: {error_count} structures failed to process")
        print("Check error messages above for details")
    
    if success_count == 0 and skipped_count == 0:
        print("\nERROR: No models were generated!")
        sys.exit(1)

    print(f"\nOutput directory: {output_dir}")
    print("✓ Database building complete!")


if __name__ == "__main__":
    main()
