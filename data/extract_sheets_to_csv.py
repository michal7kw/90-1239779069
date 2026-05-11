#!/usr/bin/env python3
"""
Extract each sheet from an Excel file as a separate CSV file.

Usage:
    python extract_sheets_to_csv.py <excel_file_path> [output_directory]

Arguments:
    excel_file_path: Path to the input Excel (.xlsx) file
    output_directory: (Optional) Directory to save CSV files. Defaults to same directory as Excel file.
"""

import pandas as pd
import sys
import os
from pathlib import Path


def extract_sheets_to_csv(excel_path, output_dir=None):
    """
    Extract each sheet from an Excel file and save as separate CSV files.
    
    Parameters:
    -----------
    excel_path : str or Path
        Path to the input Excel file
    output_dir : str or Path, optional
        Directory to save CSV files. If None, saves in same directory as Excel file.
    """
    excel_path = Path(excel_path)
    
    # Validate input file exists
    if not excel_path.exists():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    
    # Set output directory
    if output_dir is None:
        output_dir = excel_path.parent
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read Excel file with all sheets
    print(f"Reading Excel file: {excel_path}")
    excel_file = pd.ExcelFile(excel_path)
    
    # Process each sheet
    print(f"Found {len(excel_file.sheet_names)} sheet(s)")
    
    for sheet_name in excel_file.sheet_names:
        # Read the sheet
        df = pd.read_excel(excel_path, sheet_name=sheet_name)
        
        # Create safe filename from sheet name
        safe_sheet_name = sheet_name.replace('/', '-').replace('\\', '-').replace(':', '-')
        csv_filename = f"{safe_sheet_name}.csv"
        csv_path = output_dir / csv_filename
        
        # Save to CSV
        df.to_csv(csv_path, index=False, encoding='utf-8')
        print(f"  [OK] Saved sheet '{sheet_name}' -> {csv_path} ({len(df)} rows, {len(df.columns)} columns)")
    
    print(f"\n[OK] Successfully extracted {len(excel_file.sheet_names)} sheet(s) to: {output_dir}")


def main():
    """Main entry point for command-line usage."""
    if len(sys.argv) < 2:
        print("Error: Excel file path is required")
        print("\nUsage:")
        print("  python extract_sheets_to_csv.py <excel_file_path> [output_directory]")
        print("\nExample:")
        print("  python extract_sheets_to_csv.py data.xlsx")
        print("  python extract_sheets_to_csv.py data.xlsx output_csv/")
        sys.exit(1)
    
    excel_path = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    
    try:
        extract_sheets_to_csv(excel_path, output_dir)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
