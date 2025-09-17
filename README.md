# BOLDistilled to APSCALE Converter

This script converts [BOLDistilled](https://boldsystems.org/data/boldistilled/) FASTA files with `ProcessID|BIN` headers and the associated taxonomy TSV to [APSCALE](https://github.com/DominikBuchner/apscale) formatted blast+ database, specifically:
1. Clean BLAST database (ProcessID only headers)  
2. APSCALE-format taxonomy parquet file

## Requirements

- Python 3.6+
- pandas
- pyarrow (for parquet support)
- BLAST+ tools (makeblastdb must be in PATH)

Install dependencies:
```bash
pip install pandas pyarrow
```

## Usage

### Basic Usage
```bash
python bold_to_apscale.py -i sequences.fasta -t taxonomy.tsv
```

### Custom Database Name
```bash
python bold_to_apscale.py -i sequences.fasta -t taxonomy.tsv -o my_database
```

### Specify Output Directory
```bash
python bold_to_apscale.py -i sequences.fasta -t taxonomy.tsv -o db --output-dir /path/to/output
```

### Verbose Mode
```bash
python bold_to_apscale.py -i sequences.fasta -t taxonomy.tsv -v
```

### Dry Run (Show What Would Be Done)
```bash
python bold_to_apscale.py -i sequences.fasta -t taxonomy.tsv --dry-run
```

## Arguments

### Required
- `-i, --input`: Input FASTA file with ProcessID|BIN headers
- `-t, --taxonomy`: Input taxonomy TSV file with BIN mappings

### Optional
- `-o, --output`: Output database name (default: db)
- `--output-dir`: Output directory (default: current directory)
- `--title`: Database title for BLAST (default: "BOLD Database")
- `--temp-dir`: Temporary directory for intermediate files
- `--keep-temp`: Keep temporary files after completion
- `-v, --verbose`: Enable verbose logging
- `--dry-run`: Show what would be done without executing

## Input File Formats

### FASTA File
Headers must be in format: `>ProcessID|BIN`
```
>BLPAA17045-20|BOLD:AAA0017
TACACTATACTTCATTTTTGGAATTTGAGCTGGAATAGTTGGAACTTCTCTTAGTTT...
>MHMYD2007-09|BOLD:AAA0037
TACTTTATATTTTATTTTTGGAATTTGAGCAGGTATATTAGGTACTTCTCTTAGTTT...
```

### Taxonomy TSV File
Must contain columns: `bin`, `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`
```
bin	kingdom	phylum	class	order	family	genus	species
BOLD:AAA0017	Animalia	Arthropoda	Insecta	Lepidoptera	Depressariidae	Cerconota	
BOLD:AAA0037	Animalia	Arthropoda	Insecta	Lepidoptera	Gelechiidae	Anacampsis	
```

## Output Files

### BLAST Database
Multiple files with specified name:
- `db.ndb`, `db.nhr`, `db.nin`, `db.nsq`, etc.

### Taxonomy File
APSCALE-format parquet file:
- `db_taxonomy.parquet.snappy`

## Example
```bash
# Convert BOL data to BLAST database
python bold_to_apscale.py \
  -i BOLDistilled_COI_Jul2025_SEQUENCES.fasta \
  -t BOLDistilled_COI_Jul2025_TAXONOMY.tsv \
  -o db \
  --output-dir ./databases \
  -v

# Output:
# ./databases/db.* (BLAST database files)
# ./databases/db_taxonomy.parquet.snappy
```

