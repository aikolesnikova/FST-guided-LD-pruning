# FST-guided LD pruning

A Python command-line tool for prioritising population-differentiated genetic variants using an iterative, FST-guided linkage disequilibrium (LD) pruning algorithm. Candidate variants are ranked by **MaxFst**, and pairwise LD is calculated with **PLINK 1.9**. By default, candidates with **r² > 0.05** within a **10 Mb window** of a selected lead variant are excluded.

## Algorithm

1. Rank candidate variants by MaxFst in descending order.
2. Select the highest-ranked remaining variant as a lead SNP.
3. Calculate pairwise LD between the lead SNP and nearby variants using PLINK 1.9.
4. Exclude remaining candidates with **r² > 0.05** relative to the lead SNP, using the default 10 Mb window.
5. Repeat until all candidates have been selected as leads or excluded by a lead SNP.

Ties in MaxFst are resolved by chromosome, genomic position and SNP ID. Only variants provided in the candidate score file can be selected or excluded; the PLINK genotype dataset is used to calculate LD and is not modified. The LD threshold is strict: a variant with r² exactly equal to 0.05 is not excluded.

## Requirements

- Python 3.9+ (standard library only)
- PLINK 1.9, available on `PATH` or specified with `--plink`
- A PLINK binary genotype dataset (`.bed`, `.bim`, `.fam`) with unique variant IDs
- A whitespace-delimited candidate score file containing seven columns:

  ```text
  CHR  SNP  CM  POS  A1  A2  MaxFst
  ```

The first six fields follow the PLINK BIM format, and the seventh provides the numeric MaxFst score used for prioritisation. Candidate SNP IDs must be present in the genotype dataset's `.bim` file.

An illustrative score file is included at [`examples/candidate_scores.tsv`](examples/candidate_scores.tsv). Genotype data are not included.

## Usage

Run the script with a PLINK dataset prefix, a candidate score file and an output directory:

```bash
python scripts/fst_ld_pruning.py \
  --bfile data/chr1 \
  --scores data/chr1_with_MaxFst_0.bim \
  --outdir results/chr1
```

The `--bfile` argument takes the shared prefix of the `.bed`, `.bim` and `.fam` files, without the extension. The input paths above are examples and should be replaced with paths to your data.

### Options

| Option | Default | Description |
| --- | --- | --- |
| `--bfile` | Required | PLINK BED/BIM/FAM file prefix |
| `--scores` | Required | Seven-column candidate file with MaxFst scores |
| `--outdir` | Required | Directory for output tables |
| `--r2-threshold` | `0.05` | Exclude candidates with r² strictly greater than this value |
| `--window-kb` | `10000` | LD window around each lead SNP, in kb |
| `--ld-window` | `99999` | PLINK count-based LD window limit |
| `--threads` | `1` | Threads passed to PLINK |
| `--memory-mb` | `8000` | Memory limit passed to PLINK, in MB |
| `--plink` | `plink` | PLINK 1.9 executable or path |
| `--overwrite` | Off | Allow output to an existing non-empty directory |

## Outputs

The script writes three tab-separated files to `--outdir`:

- **`lead_variants.tsv`**: retained lead SNPs in selection order, with iteration number, chromosome, position, SNP ID and MaxFst.
- **`excluded_variants.tsv`**: excluded candidate SNPs, the lead SNP responsible for exclusion, pairwise r² and lead iteration.
- **`summary.tsv`**: numbers of input candidates, lead variants and excluded variants, together with the r² threshold and LD window.

Every input candidate appears exactly once in either the retained or excluded set.

## Testing

Run the included unit tests from the repository root:

```bash
python -m unittest discover -s tests -v
```

The tests cover candidate ranking, tie-breaking, LD output parsing, the strict r² threshold and an end-to-end workflow using a mock PLINK executable. They do not require genotype data or an installed PLINK binary.

## SLURM execution

An example batch submission script is provided at [`profiles/slurm.sh`](profiles/slurm.sh). Adapt the working directory, environment and resource settings for your HPC system, and set `BFILE`, `SCORES` and `OUTDIR` before submission.

## Data availability

The repository includes illustrative candidate scores and automated tests. No individual-level genotype data are distributed.
