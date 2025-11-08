#!/usr/bin/env python3
"""
BOLDistilled to APSCALE Converter
========================

Converts BOLD FASTA files with ProcessID|BIN headers to:
1. Clean BLAST database (ProcessID only headers)
2. APSCALE-format taxonomy parquet file

Author: Claude (Anthropic)
Version: 1.0
"""

import argparse
import os
import sys
import pandas as pd
import subprocess
from pathlib import Path
import logging
from typing import Dict, Optional, Tuple

def setup_logging(verbose: bool = False) -> None:
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Convert BOLD FASTA to BLAST database with APSCALE-format taxonomy',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i sequences.fasta -t taxonomy.tsv
  %(prog)s -i sequences.fasta -t taxonomy.tsv -o my_database
  %(prog)s -i sequences.fasta -t taxonomy.tsv -o my_db --output-dir /path/to/output
        """
    )
    
    # Required arguments
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input FASTA file with ProcessID|BIN headers (required)'
    )
    
    parser.add_argument(
        '-t', '--taxonomy',
        required=True,
        help='Input taxonomy TSV file with BIN mappings (required)'
    )
    
    # Optional arguments
    parser.add_argument(
        '-o', '--output',
        default='db',
        help='Output database name (default: db)'
    )
    
    parser.add_argument(
        '--output-dir',
        default='.',
        help='Output directory (default: current directory)'
    )
    
    parser.add_argument(
        '--title',
        default='BOLD Database',
        help='Database title for BLAST (default: "BOLD Database")'
    )
    
    parser.add_argument(
        '--temp-dir',
        help='Temporary directory for intermediate files (default: system temp)'
    )
    
    parser.add_argument(
        '--keep-temp',
        action='store_true',
        help='Keep temporary files after completion'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be done without executing'
    )
    
    return parser.parse_args()

def validate_inputs(args: argparse.Namespace) -> None:
    """Validate input files and directories."""
    # Check input FASTA file
    if not os.path.isfile(args.input):
        raise FileNotFoundError(f"Input FASTA file not found: {args.input}")
    
    # Check taxonomy file
    if not os.path.isfile(args.taxonomy):
        raise FileNotFoundError(f"Taxonomy file not found: {args.taxonomy}")
    
    # Check/create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Check if makeblastdb is available
    try:
        result = subprocess.run(['makeblastdb', '-help'], 
                              capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError("makeblastdb not found or not working")
    except FileNotFoundError:
        raise RuntimeError("makeblastdb not found in PATH. Please install BLAST+")
    
    logging.info("Input validation completed successfully")

def clean_fasta_headers(input_fasta: str, output_fasta: str) -> Tuple[int, Dict[str, str]]:
    """
    Clean FASTA headers by removing BIN IDs and create mapping.
    Handles duplicate ProcessIDs by adding unique suffixes (_1, _2, etc.).
    
    Args:
        input_fasta: Path to input FASTA with ProcessID|BIN headers
        output_fasta: Path to output FASTA with clean ProcessID headers
        
    Returns:
        Tuple of (sequence_count, process_to_bin_mapping)
    """
    logging.info("Cleaning FASTA headers...")
    
    process_to_bin = {}
    sequence_count = 0
    seen_process_ids = {}  # Track occurrences of each process ID
    duplicate_count = 0
    
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                header = line.strip()[1:]  # Remove '>'
                if '|' in header:
                    process_id, bin_id = header.split('|', 1)
                    
                    # Handle duplicates by adding suffix
                    if process_id in seen_process_ids:
                        seen_process_ids[process_id] += 1
                        unique_process_id = f"{process_id}_{seen_process_ids[process_id]}"
                        duplicate_count += 1
                        
                        if duplicate_count <= 10:  # Log first few duplicates
                            logging.warning(f"  Duplicate ProcessID '{process_id}' renamed to '{unique_process_id}'")
                        elif duplicate_count == 11:
                            logging.warning(f"  (Further duplicate warnings suppressed...)")
                    else:
                        seen_process_ids[process_id] = 0
                        unique_process_id = process_id
                    
                    process_to_bin[unique_process_id] = bin_id
                    outfile.write(f">{unique_process_id}\n")
                    sequence_count += 1
                    
                    if sequence_count % 100000 == 0:
                        logging.info(f"  Processed {sequence_count:,} sequences...")
                else:
                    # Header without '|' - use as-is but still check for duplicates
                    if header in seen_process_ids:
                        seen_process_ids[header] += 1
                        unique_header = f"{header}_{seen_process_ids[header]}"
                        duplicate_count += 1
                        
                        if duplicate_count <= 10:
                            logging.warning(f"  Duplicate ProcessID '{header}' renamed to '{unique_header}'")
                        elif duplicate_count == 11:
                            logging.warning(f"  (Further duplicate warnings suppressed...)")
                        
                        process_to_bin[unique_header] = header
                        outfile.write(f">{unique_header}\n")
                    else:
                        seen_process_ids[header] = 0
                        process_to_bin[header] = header
                        outfile.write(line)
                    
                    sequence_count += 1
            else:
                outfile.write(line)
    
    if duplicate_count > 0:
        logging.warning(f"Found and renamed {duplicate_count:,} duplicate ProcessIDs")
    
    logging.info(f"Cleaned {sequence_count:,} sequences")
    return sequence_count, process_to_bin

def load_taxonomy_data(taxonomy_file: str) -> pd.DataFrame:
    """Load and validate taxonomy data."""
    logging.info("Loading taxonomy data...")
    
    try:
        taxonomy_df = pd.read_csv(taxonomy_file, sep='\t')
    except Exception as e:
        raise RuntimeError(f"Failed to load taxonomy file: {e}")
    
    # Validate required columns
    required_columns = ['bin', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    missing_columns = [col for col in required_columns if col not in taxonomy_df.columns]
    
    if missing_columns:
        raise ValueError(f"Missing required columns in taxonomy file: {missing_columns}")
    
    logging.info(f"Loaded {len(taxonomy_df):,} taxonomy records")
    return taxonomy_df

def create_apscale_taxonomy(process_to_bin: Dict[str, str], 
                          taxonomy_df: pd.DataFrame, 
                          output_file: str) -> int:
    """
    Create APSCALE-format taxonomy parquet file.
    
    Args:
        process_to_bin: Mapping from ProcessID to BIN
        taxonomy_df: Taxonomy data DataFrame
        output_file: Output parquet file path
        
    Returns:
        Number of taxonomy records created
    """
    logging.info("Creating APSCALE-format taxonomy file...")
    
    # Create BIN to taxonomy lookup
    bin_to_taxonomy = {}
    for _, row in taxonomy_df.iterrows():
        bin_id = row['bin']
        bin_to_taxonomy[bin_id] = {
            'superkingdom': row['kingdom'],
            'phylum': row['phylum'],
            'class': row['class'],
            'order': row['order'],
            'family': row['family'],
            'genus': row['genus'],
            'species': row['species']
        }
    
    logging.info(f"Created taxonomy lookup for {len(bin_to_taxonomy):,} BIN IDs")
    
    # Create final taxonomy records
    final_records = []
    missing_bins = set()
    
    for process_id, bin_id in process_to_bin.items():
        if bin_id in bin_to_taxonomy:
            tax_data = bin_to_taxonomy[bin_id]
            record = {
                'Accession': process_id,  # Process ID as accession (APSCALE style)
                'superkingdom': tax_data['superkingdom'],
                'phylum': tax_data['phylum'],
                'class': tax_data['class'],
                'order': tax_data['order'],
                'family': tax_data['family'],
                'genus': tax_data['genus'],
                'species': tax_data['species']
            }
            final_records.append(record)
        else:
            missing_bins.add(bin_id)
    
    if missing_bins:
        logging.warning(f"Missing taxonomy for {len(missing_bins)} BINs")
        if len(missing_bins) <= 10:
            logging.warning(f"Missing BINs: {list(missing_bins)}")
    
    # Create DataFrame and save
    final_df = pd.DataFrame(final_records)
    final_df = final_df.fillna('')  # Clean up NaN values
    
    # Save as parquet with Snappy compression
    final_df.to_parquet(output_file, compression='snappy', index=False)
    
    logging.info(f"Created taxonomy file with {len(final_df):,} records")
    return len(final_df)

def create_blast_database(clean_fasta: str, output_name: str, title: str) -> None:
    """Create BLAST database from clean FASTA file."""
    logging.info("Creating BLAST database...")
    
    cmd = [
        'makeblastdb',
        '-in', clean_fasta,
        '-dbtype', 'nucl',
        '-out', output_name,
        '-title', title,
        '-parse_seqids'
    ]
    
    logging.debug(f"Running command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logging.debug(f"makeblastdb output: {result.stdout}")
        logging.info("BLAST database created successfully")
    except subprocess.CalledProcessError as e:
        logging.error(f"makeblastdb failed: {e.stderr}")
        raise RuntimeError(f"Failed to create BLAST database: {e}")

def verify_outputs(db_name: str, taxonomy_file: str) -> None:
    """Verify that outputs were created correctly."""
    logging.info("Verifying outputs...")
    
    # Check BLAST database files
    db_extensions = ['.ndb', '.nhr', '.nin', '.nsq', '.ntf', '.nto']
    for ext in db_extensions:
        db_file = f"{db_name}{ext}"
        if not os.path.isfile(db_file):
            raise RuntimeError(f"BLAST database file missing: {db_file}")
    
    # Check taxonomy file
    if not os.path.isfile(taxonomy_file):
        raise RuntimeError(f"Taxonomy file missing: {taxonomy_file}")
    
    # Verify taxonomy file can be loaded
    try:
        tax_df = pd.read_parquet(taxonomy_file)
        expected_columns = ['Accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        if list(tax_df.columns) != expected_columns:
            raise RuntimeError(f"Taxonomy file has incorrect columns: {list(tax_df.columns)}")
        logging.info(f"Taxonomy file verification: {len(tax_df):,} records with correct structure")
    except Exception as e:
        raise RuntimeError(f"Failed to verify taxonomy file: {e}")
    
    logging.info("Output verification completed successfully")

def cleanup_temp_files(temp_files: list, keep_temp: bool) -> None:
    """Clean up temporary files if requested."""
    if keep_temp:
        logging.info(f"Keeping temporary files: {temp_files}")
        return
    
    for temp_file in temp_files:
        try:
            if os.path.isfile(temp_file):
                os.remove(temp_file)
                logging.debug(f"Removed temporary file: {temp_file}")
        except Exception as e:
            logging.warning(f"Failed to remove temporary file {temp_file}: {e}")

def main():
    """Main function."""
    args = parse_arguments()
    setup_logging(args.verbose)
    
    try:
        # Validate inputs
        validate_inputs(args)
        
        if args.dry_run:
            logging.info("DRY RUN MODE - showing what would be done:")
            logging.info(f"  Input FASTA: {args.input}")
            logging.info(f"  Taxonomy file: {args.taxonomy}")
            logging.info(f"  Output database: {args.output}")
            logging.info(f"  Output directory: {args.output_dir}")
            logging.info(f"  Database title: {args.title}")
            return 0
        
        # Setup paths
        output_dir = Path(args.output_dir)
        db_name = output_dir / args.output
        taxonomy_file = output_dir / f"{args.output}_taxonomy.parquet.snappy"
        
        # Setup temporary directory
        if args.temp_dir:
            temp_dir = Path(args.temp_dir)
            temp_dir.mkdir(exist_ok=True)
        else:
            temp_dir = output_dir
        
        clean_fasta = temp_dir / f"{args.output}_clean.fasta"
        temp_files = [str(clean_fasta)]
        
        logging.info("Starting BOLD to APSCALE conversion...")
        logging.info(f"Input FASTA: {args.input}")
        logging.info(f"Taxonomy file: {args.taxonomy}")
        logging.info(f"Output database: {db_name}")
        logging.info(f"Output taxonomy: {taxonomy_file}")
        
        # Step 1: Clean FASTA headers
        sequence_count, process_to_bin = clean_fasta_headers(args.input, str(clean_fasta))
        
        # Step 2: Load taxonomy data
        taxonomy_df = load_taxonomy_data(args.taxonomy)
        
        # Step 3: Create APSCALE-format taxonomy file
        taxonomy_count = create_apscale_taxonomy(process_to_bin, taxonomy_df, str(taxonomy_file))
        
        # Step 4: Create BLAST database
        create_blast_database(str(clean_fasta), str(db_name), args.title)
        
        # Step 5: Verify outputs
        verify_outputs(str(db_name), str(taxonomy_file))
        
        # Step 6: Cleanup
        cleanup_temp_files(temp_files, args.keep_temp)
        
        # Success summary
        logging.info("=" * 50)
        logging.info("CONVERSION COMPLETED SUCCESSFULLY!")
        logging.info("=" * 50)
        logging.info(f"Sequences processed: {sequence_count:,}")
        logging.info(f"Taxonomy records: {taxonomy_count:,}")
        logging.info(f"BLAST database: {db_name}")
        logging.info(f"Taxonomy file: {taxonomy_file}")
        logging.info("Ready for use with apscale or other BLAST-based software!")
        
        return 0
        
    except Exception as e:
        logging.error(f"Conversion failed: {e}")
        if args.verbose:
            import traceback
            logging.error(traceback.format_exc())
        return 1

if __name__ == '__main__':
    sys.exit(main())
