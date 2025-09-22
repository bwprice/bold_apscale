#!/usr/bin/env python3
import pandas as pd
import argparse
import os

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Compare taxonomic identifications from MIDORI and BOLD')
    parser.add_argument('midori_file', help='Path to the BLASTN MIDORI results CSV')
    parser.add_argument('bold_file', help='Path to the BOLD results CSV')
    parser.add_argument('output_file', help='Path to output CSV file')
    return parser.parse_args()

def get_taxonomy_level(row):
    """
    Determine the taxonomic resolution level for a row
    Returns a numeric value (higher means more specific)
    """
    # Taxonomic levels in order of increasing specificity
    levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    
    # Start from most specific
    for i, level in enumerate(reversed(levels)):
        if level in row and row[level] and pd.notna(row[level]) and row[level] != '' and row[level] != 'no-match':
            return 7 - i  # 7 for Species, 6 for Genus, etc.
    return 0  # No taxonomic info

def find_common_level(midori_row, bold_row):
    """
    Find the lowest common taxonomic level between two identifications
    Returns the level name and both values at that level
    """
    # Taxonomic levels in order of increasing specificity
    levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    
    common_level = None
    midori_value = None
    bold_value = None
    
    # Check each level from most specific to least specific
    for level in reversed(levels):
        # Skip Kingdom as it might be missing in BOLD
        if level == 'Kingdom' and level not in bold_row:
            continue
            
        if level in midori_row and level in bold_row:
            midori_value = midori_row[level]
            bold_value = bold_row[level]
            
            # If both have values at this level (and not 'no-match')
            if (pd.notna(midori_value) and midori_value != '' and midori_value != 'no-match' and 
                pd.notna(bold_value) and bold_value != '' and bold_value != 'no-match'):
                
                # If they match, continue to more specific level
                if midori_value.lower() == bold_value.lower():
                    common_level = level
                    continue
                # If they don't match, return the previous level
                else:
                    # If we're at Phylum and it doesn't match, return Phylum
                    if level == 'Phylum':
                        return 'Phylum', midori_value, bold_value
                    # For other levels, go up one level
                    else:
                        previous_idx = levels.index(level) - 1
                        if previous_idx >= 0:
                            previous_level = levels[previous_idx]
                            if previous_level in midori_row and previous_level in bold_row:
                                return previous_level, midori_row[previous_level], bold_row[previous_level]
                        return None, None, None
    
    # If we got through all levels without finding a conflict, return the most specific common level
    return common_level, midori_value, bold_value

def get_best_identification(midori_row, bold_row):
    """
    Determine the best identification between MIDORI and BOLD results
    Based on:
    1. If IDs differ, use lowest common level
    2. If resolutions differ, use the most specific
    """
    # Check if taxonomic identifications differ
    midori_level = get_taxonomy_level(midori_row)
    bold_level = get_taxonomy_level(bold_row)
    
    # Special case: if one source has all 'no-match' entries, use the other source
    bold_no_match = all(
        (level in bold_row and (bold_row[level] == 'no-match' or bold_row[level] == '' or pd.isna(bold_row[level])))
        for level in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        if level in bold_row
    )
    
    midori_no_match = all(
        (level in midori_row and (midori_row[level] == 'no-match' or midori_row[level] == '' or pd.isna(midori_row[level])))
        for level in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        if level in midori_row
    )
    
    if bold_no_match and not midori_no_match:
        # Use MIDORI identification
        result = {}
        source = "MIDORI"
        for level in ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:
            if level in midori_row:
                result[level] = midori_row[level]
        result['ID_source'] = source
        
        # Add MIDORI-specific columns
        result['MIDORI_Similarity'] = midori_row.get('Similarity', '')
        result['MIDORI_evalue'] = midori_row.get('evalue', '')
        result['MIDORI_Flag'] = midori_row.get('Flag', '')
        result['MIDORI_Ambigous_taxa'] = midori_row.get('Ambigous taxa', '')
        result['MIDORI_Status'] = midori_row.get('Status', '')
        
        # Add empty BOLD columns
        result['BOLD_pct_identity'] = ''
        result['BOLD_status'] = ''
        result['BOLD_records'] = ''
        result['BOLD_selected_level'] = ''
        result['BOLD_BIN'] = ''
        result['BOLD_flags'] = ''
        
        return result
    
    if midori_no_match and not bold_no_match:
        # Use BOLD identification
        result = {}
        source = "BOLD"
        # Add Kingdom from MIDORI if available (since BOLD might not have it)
        if ('Kingdom' in midori_row and pd.notna(midori_row['Kingdom']) and 
            midori_row['Kingdom'] != '' and midori_row['Kingdom'] != 'no-match'):
            result['Kingdom'] = midori_row['Kingdom']
        
        for level in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:
            if level in bold_row:
                result[level] = bold_row[level]
        result['ID_source'] = source
        
        # Add empty MIDORI columns
        result['MIDORI_Similarity'] = ''
        result['MIDORI_evalue'] = ''
        result['MIDORI_Flag'] = ''
        result['MIDORI_Ambigous_taxa'] = ''
        result['MIDORI_Status'] = ''
        
        # Add BOLD-specific columns
        result['BOLD_pct_identity'] = bold_row.get('pct_identity', '')
        result['BOLD_status'] = bold_row.get('status', '')
        result['BOLD_records'] = bold_row.get('records', '')
        result['BOLD_selected_level'] = bold_row.get('selected_level', '')
        result['BOLD_BIN'] = bold_row.get('BIN', '')
        result['BOLD_flags'] = bold_row.get('flags', '')
        
        return result
    
    # Find the common taxonomic level
    common_level, midori_value, bold_value = find_common_level(midori_row, bold_row)
    
    # If no common level found, use highest common taxonomy (e.g., at Kingdom or Phylum level)
    if common_level is None:
        # Default to the first common taxonomy level
        for level in ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:
            if level in midori_row and level in bold_row:
                if (pd.notna(midori_row[level]) and pd.notna(bold_row[level]) and 
                    midori_row[level] != '' and bold_row[level] != '' and
                    midori_row[level] != 'no-match' and bold_row[level] != 'no-match'):
                    if midori_row[level].lower() == bold_row[level].lower():
                        common_level = level
                        midori_value = midori_row[level]
                        bold_value = bold_row[level]
                        break
    
    # Create a new row with the best identification
    result = {}
    
    # Copy all taxonomic levels from the identification with better resolution
    if midori_level > bold_level:
        source = "MIDORI"
        for level in ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:
            if level in midori_row:
                result[level] = midori_row[level]
    else:
        source = "BOLD"
        # Add Kingdom from MIDORI if available (since BOLD might not have it)
        if ('Kingdom' in midori_row and pd.notna(midori_row['Kingdom']) and 
            midori_row['Kingdom'] != '' and midori_row['Kingdom'] != 'no-match'):
            result['Kingdom'] = midori_row['Kingdom']
        
        for level in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:
            if level in bold_row:
                result[level] = bold_row[level]
    
    # If identifications differ but are at the same level, use the lowest common taxonomy
    if common_level and midori_level == bold_level and midori_value != bold_value:
        # Reset taxonomic levels after the common level
        levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        common_index = levels.index(common_level)
        
        # Clear values for more specific levels
        for i in range(common_index + 1, len(levels)):
            if levels[i] in result:
                result[levels[i]] = ''
    
    # Add source information
    result['ID_source'] = source
    
    # Add additional columns that differ between the two sources
    result['MIDORI_Similarity'] = midori_row.get('Similarity', '')
    result['MIDORI_evalue'] = midori_row.get('evalue', '')
    result['MIDORI_Flag'] = midori_row.get('Flag', '')
    result['MIDORI_Ambigous_taxa'] = midori_row.get('Ambigous taxa', '')
    result['MIDORI_Status'] = midori_row.get('Status', '')
    
    result['BOLD_pct_identity'] = bold_row.get('pct_identity', '')
    result['BOLD_status'] = bold_row.get('status', '')
    result['BOLD_records'] = bold_row.get('records', '')
    result['BOLD_selected_level'] = bold_row.get('selected_level', '')
    result['BOLD_BIN'] = bold_row.get('BIN', '')
    result['BOLD_flags'] = bold_row.get('flags', '')
    
    return result

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Check if input files exist
    if not os.path.isfile(args.midori_file):
        print(f"Error: MIDORI file '{args.midori_file}' does not exist")
        return
    if not os.path.isfile(args.bold_file):
        print(f"Error: BOLD file '{args.bold_file}' does not exist")
        return
    
    # Read the input CSV files
    print(f"Reading MIDORI file: {args.midori_file}")
    midori_df = pd.read_csv(args.midori_file)
    print(f"Reading BOLD file: {args.bold_file}")
    bold_df = pd.read_csv(args.bold_file)
    
    # Ensure the unique ID columns match
    if 'unique ID' in midori_df.columns and 'id' in bold_df.columns:
        # Rename 'id' to 'unique ID' in the BOLD dataframe for consistency
        bold_df.rename(columns={'id': 'unique ID'}, inplace=True)
    
    # Check if 'unique ID' exists in both dataframes
    if 'unique ID' not in midori_df.columns or 'unique ID' not in bold_df.columns:
        print("Error: 'unique ID' column not found in one or both files")
        return
    
    # Create a new dataframe for the combined results
    combined_results = []
    
    # Count for progress updates
    total_records = len(midori_df)
    processed = 0
    
    # Process each record in the MIDORI dataframe
    for midori_idx, midori_row in midori_df.iterrows():
        unique_id = midori_row['unique ID']
        
        # Find corresponding row in BOLD dataframe
        bold_match = bold_df[bold_df['unique ID'] == unique_id]
        
        if len(bold_match) > 0:
            # Get the best identification between the two sources
            bold_row = bold_match.iloc[0].to_dict()
            best_id = get_best_identification(midori_row.to_dict(), bold_row)
            
            # Add the unique ID
            best_id['unique ID'] = unique_id
            
            # Add to results
            combined_results.append(best_id)
        else:
            # If no match in BOLD, use MIDORI identification
            result = midori_row.to_dict()
            result['ID_source'] = 'MIDORI'
            combined_results.append(result)
        
        # Update progress
        processed += 1
        if processed % 100 == 0 or processed == total_records:
            print(f"Processed {processed}/{total_records} records")
    
    # Convert to DataFrame and save to CSV
    result_df = pd.DataFrame(combined_results)
    
    # Reorder columns for better readability
    column_order = ['unique ID', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'ID_source']
    # Add all remaining columns
    for col in result_df.columns:
        if col not in column_order:
            column_order.append(col)
    
    result_df = result_df[column_order]
    
    # Save to output file
    print(f"Saving results to: {args.output_file}")
    result_df.to_csv(args.output_file, index=False)
    print("Done!")

if __name__ == "__main__":
    main()
