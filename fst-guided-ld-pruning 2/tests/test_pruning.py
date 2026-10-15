import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("fst_ld_pruning", ROOT / "scripts" / "fst_ld_pruning.py")
mod = importlib.util.module_from_spec(spec)
import sys
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


class PruningTests(unittest.TestCase):
    def test_ranking_and_tie_breaking(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "score.bim"
            p.write_text("1 B 0 200 A G 0.8\n1 C 0 100 A G 0.8\n1 A 0 50 A G 0.9\n")
            self.assertEqual([x.snp for x in mod.load_candidates(p)], ["A", "C", "B"])

    def test_duplicate_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "score.bim"
            p.write_text("1 A 0 200 A G 0.8\n1 A 0 100 A G 0.7\n")
            with self.assertRaisesRegex(ValueError, "duplicate"):
                mod.load_candidates(p)

    def test_ld_parsing_and_self_pair(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "result.ld"
            p.write_text("CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\n1 100 A 1 100 A 1\n1 100 A 1 200 B 0.05\n1 100 A 1 300 C 0.15\n")
            self.assertEqual(mod.parse_ld(p, "A"), {"A": 1.0, "B": 0.05, "C": 0.15})

    def test_full_workflow_with_mock_plink(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            bfile = root / "cohort"
            (root / "cohort.bed").write_bytes(b"mockbed")
            (root / "cohort.fam").write_text("F I 0 0 0 -9\n")
            (root / "cohort.bim").write_text("1 A 0 100 A G\n1 B 0 200 A G\n1 C 0 300 A G\n1 D 0 400 A G\n")
            score = root / "scores.bim"
            score.write_text("1 A 0 100 A G 0.9\n1 B 0 200 A G 0.8\n1 C 0 300 A G 0.7\n1 D 0 400 A G 0.6\n")
            executable = root / "mock_plink"
            executable.write_text("""#!/usr/bin/env python3
import sys
from pathlib import Path
args=sys.argv
lead=args[args.index('--ld-snp')+1]
out=Path(args[args.index('--out')+1]+'.ld')
rows={'A':['1 100 A 1 100 A 1','1 100 A 1 200 B 0.2','1 100 A 1 300 C 0.05'],
      'C':['1 300 C 1 300 C 1','1 300 C 1 400 D 0.15']}
out.write_text('CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\\n'+'\\n'.join(rows.get(lead, []))+'\\n')
""")
            executable.chmod(0o755)
            out = root / "out"
            mod.run(SimpleNamespace(bfile=str(bfile), scores=str(score), outdir=str(out),
                r2_threshold=0.05, window_kb=10000, ld_window=99999,
                threads=1, memory_mb=8000, plink=str(executable), overwrite=False))
            rows = (out / "lead_variants.tsv").read_text().splitlines()
            self.assertEqual([r.split("\t")[3] for r in rows[1:]], ["A", "C"])
            removed = (out / "excluded_variants.tsv").read_text()
            self.assertIn("B\tA\t0.2", removed)
            self.assertIn("D\tC\t0.15", removed)
            self.assertNotIn("\t0.05\t", removed)
            self.assertIn("4\t2\t2\t0.05\t10000", (out / "summary.tsv").read_text())


if __name__ == "__main__":
    unittest.main()
