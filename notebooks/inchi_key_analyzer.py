import os
import pandas as pd
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent

def extract_file_info(file_path):
    """Extract chromatography, polarity, and atlas type from file name (case-insensitive, more tolerant)."""
    filename = file_path.name
    fname = filename.lower()
    
    # Extract chromatography
    if fname.startswith('c18_'):
        chrom = 'C18'
    elif fname.startswith('hilic_'):
        chrom = 'HILIC'
    else:
        return None, None, None
    
    # Extract polarity
    if '_negative' in fname:
        polarity = 'negative'
    elif '_positive' in fname:
        polarity = 'positive'
    else:
        return None, None, None
    
    # Extract atlas type
    if 'ema-standards' in fname:
        atlas_type = 'EMA'
    elif 'istd' in fname:
        atlas_type = 'ISTD'
    elif 'qc' in fname:
        # Only use QCv7 for QC category; accept variants like 'qcv7' case-insensitively
        if 'qcv7' in fname:
            atlas_type = 'QC'
        else:
            return None, None, None  # Skip non-QCv7 QC files
    else:
        return None, None, None
    
    return chrom, polarity, atlas_type

def process_files_in_subdirectory(subdir):
    """Process all files in a specific subdirectory and categorize inchi_key values."""
    data = defaultdict(set)
    file_counts = defaultdict(int)
    
    # Process both C18 and HILIC directories
    for base_dir in ['C18', 'HILIC']:
        target_dir = ROOT / base_dir / subdir
        if not target_dir.exists():
            continue
            
        for file_path in target_dir.glob('*.tsv'):
            if not file_path.is_file():
                continue
                
            chrom, polarity, atlas_type = extract_file_info(file_path)
            if not all([chrom, polarity, atlas_type]):
                print(f"Skipping {file_path.name}: filename did not match expected pattern")
                continue
                
            try:
                df = pd.read_csv(file_path, sep='\t', low_memory=False)
                # Normalize column names and find inchi_key case-insensitively
                df.columns = df.columns.astype(str)
                df.columns = df.columns.str.strip()
                ik_col = next((c for c in df.columns if c.strip().lower() == 'inchi_key'), None)
                
                if ik_col is None:
                    print(f"Skipping {file_path.name}: no InChI key column found (columns: {list(df.columns)})")
                    continue
                
                inchi_keys = set(df[ik_col].dropna().astype(str))
                # Remove empty strings and 'nan' (case-insensitive)
                inchi_keys = {k for k in inchi_keys if k and k.strip() and k.strip().lower() != 'nan'}
                
                category = (chrom, polarity, atlas_type)
                data[category].update(inchi_keys)
                file_counts[category] += 1
                
                print(f"Processed {file_path.name}: found {len(inchi_keys)} unique InChI keys")
                    
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
    
    return data, file_counts

def generate_report_for_subdirectory(subdir):
    """Generate report data for a specific subdirectory."""
    data, file_counts = process_files_in_subdirectory(subdir)
    
    # Create detailed report data
    report_rows = []
    
    # Individual category rows
    category_totals = {}
    for (chrom, polarity, atlas_type), inchi_keys in data.items():
        count = len(inchi_keys)
        files = file_counts[(chrom, polarity, atlas_type)]
        
        report_rows.append({
            'Category': f"{chrom}_{polarity}_{atlas_type}",
            'Chromatography': chrom,
            'Polarity': polarity,
            'Atlas_Type': atlas_type,
            'Unique_InChI_Keys': count,
            'Files_Processed': files
        })
        
        category_totals[(chrom, polarity, atlas_type)] = inchi_keys
    
    # Add empty row for separation
    report_rows.append({'Category': '', 'Chromatography': '', 'Polarity': '', 'Atlas_Type': '', 'Unique_InChI_Keys': '', 'Files_Processed': ''})
    report_rows.append({'Category': 'SUMMARY TOTALS', 'Chromatography': '', 'Polarity': '', 'Atlas_Type': '', 'Unique_InChI_Keys': '', 'Files_Processed': ''})
    
    # Chromatography totals
    for chrom in ['C18', 'HILIC']:
        chrom_keys = set()
        chrom_files = 0
        for (c, p, a), keys in category_totals.items():
            if c == chrom:
                chrom_keys.update(keys)
                chrom_files += file_counts[(c, p, a)]
        
        if chrom_files > 0:  # Only add if we have files for this chromatography
            report_rows.append({
                'Category': f"{chrom}_Total",
                'Chromatography': chrom,
                'Polarity': 'All',
                'Atlas_Type': 'All',
                'Unique_InChI_Keys': len(chrom_keys),
                'Files_Processed': chrom_files
            })
    
    # Polarity totals
    for polarity in ['negative', 'positive']:
        pol_keys = set()
        pol_files = 0
        for (c, p, a), keys in category_totals.items():
            if p == polarity:
                pol_keys.update(keys)
                pol_files += file_counts[(c, p, a)]
        
        if pol_files > 0:  # Only add if we have files for this polarity
            report_rows.append({
                'Category': f"{polarity}_Total",
                'Chromatography': 'All',
                'Polarity': polarity,
                'Atlas_Type': 'All',
                'Unique_InChI_Keys': len(pol_keys),
                'Files_Processed': pol_files
            })
    
    # Atlas type totals
    for atlas_type in ['EMA', 'ISTD', 'QC']:
        atlas_keys = set()
        atlas_files = 0
        for (c, p, a), keys in category_totals.items():
            if a == atlas_type:
                atlas_keys.update(keys)
                atlas_files += file_counts[(c, p, a)]
        
        if atlas_files > 0:  # Only add if we have files for this atlas type
            report_rows.append({
                'Category': f"{atlas_type}_Total",
                'Chromatography': 'All',
                'Polarity': 'All',
                'Atlas_Type': atlas_type,
                'Unique_InChI_Keys': len(atlas_keys),
                'Files_Processed': atlas_files
            })
    
    # Grand total
    all_keys = set()
    all_files = 0
    for keys in category_totals.values():
        all_keys.update(keys)
    for files in file_counts.values():
        all_files += files
    
    if all_files > 0:  # Only add grand total if we have files
        report_rows.append({
            'Category': 'GRAND_TOTAL',
            'Chromatography': 'All',
            'Polarity': 'All',
            'Atlas_Type': 'All',
            'Unique_InChI_Keys': len(all_keys),
            'Files_Processed': all_files
        })
    
    return report_rows, len(all_keys), all_files

def generate_excel_report():
    """Generate Excel report with separate tabs for everything and production."""
    # Create Excel writer
    with pd.ExcelWriter('inchi_key_analysis_report.xlsx', engine='openpyxl') as writer:
        
        for subdir in ['everything', 'production']:
            print(f"\nProcessing {subdir} subdirectory...")
            report_rows, total_keys, total_files = generate_report_for_subdirectory(subdir)
            
            if total_files > 0:
                # Create DataFrame and write to Excel tab
                df_report = pd.DataFrame(report_rows)
                df_report.to_excel(writer, sheet_name=subdir, index=False)
                
                print(f"  Total unique InChI keys in {subdir}: {total_keys}")
                print(f"  Total files processed in {subdir}: {total_files}")
                
                # Print quick summary for this subdirectory
                print(f"\n  QUICK SUMMARY for {subdir}:")
                print(f"  {'Category':<25} | {'Keys':<10} | {'Files':<5}")
                print(f"  {'-'*25} | {'-'*10} | {'-'*5}")
                for row in report_rows:
                    if row['Category'] and row['Unique_InChI_Keys'] != '':
                        print(f"  {row['Category']:<25} | {str(row['Unique_InChI_Keys']):<10} | {str(row['Files_Processed']):<5}")
            else:
                print(f"  No files found in {subdir} subdirectory")
    
    print(f"\nExcel report generated: inchi_key_analysis_report.xlsx")

if __name__ == "__main__":
    generate_excel_report()