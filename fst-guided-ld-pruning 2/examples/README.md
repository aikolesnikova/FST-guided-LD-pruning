# Example candidate scores

`candidate_scores.tsv` contains four **illustrative** chromosome 1 entries in the seven-column BIM-style format expected by `--scores`. The identifiers, positions and scores are fabricated and are provided for input-format inspection only. Running the full pipeline additionally requires PLINK BED/BIM/FAM genotype data with matching SNP IDs; no genotype data are distributed here.

The automated integration test (`python -m unittest discover -s tests -v`) uses a mock PLINK executable, and is self-contained. It verifies a four-candidate example, producing two leads and two exclusions, including strict handling of the r² = 0.05 boundary.
