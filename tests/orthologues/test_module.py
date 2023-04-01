from pathlib import Path
import pandas as pd
import unittest
import sys

sys.path.append("../orthologues")
from orthologues import Orthologues_identifier

class TestModule(unittest.TestCase):
    def test_create_directories(self):
        fasta_dir = Path("../data/Protein_fasta_files/")
        out_dir = Path("../data/")
        blast_threads = 4
        blast_evalue = 1e5
        blast_method = "diamond"
        filter = True
        refs = ["GCF_000009045.1", "GCF_019704415.1"]
        for ref in refs:
            instance = Orthologues_identifier(fasta_dir, out_dir, ref, blast_threads, blast_evalue, blast_method, filter)
            instance._set_directories()
            directories = [
                    out_dir / "Orthologues"/ ref / "Blast_results",
                    out_dir / "Orthologues"/ ref / "Best_reciprocal_hits",
                    out_dir / "Orthologues"/ ref / "Blast_DB",
                    ]
            for directory in directories:
                self.assertIs(directory.exists(), True)

    def test_get_blast_binaries(self):
        fasta_dir = Path("data/Protein_fasta_files/")
        out_dir = Path("data/")
        blast_threads = 4
        blast_evalue = 1e5
        blast_method = "diamond"
        filter = True
        ref = "GCF_000009045.1"
        instance = Orthologues_identifier(fasta_dir, out_dir, ref, blast_threads, blast_evalue, blast_method, filter)
        instance._get_blast_binaries()
        self.assertIs(instance.blast_bin.exists(), True)
        cmd = f"{instance.blast_bin} help"
        status = instance._execute_cmd(cmd)
        self.assertEqual(status, 0)

    def test_create_blast_db(self):
        fasta_dir = Path("data/Protein_fasta_files/")
        out_dir = Path("data/")
        blast_threads = 4
        blast_evalue = 1e5
        blast_method = "diamond"
        filter = True
        ref = "GCF_000009045.1"
        instance = Orthologues_identifier(fasta_dir, out_dir, ref, blast_threads, blast_evalue, blast_method, filter)
        fasta_file = fasta_dir / (ref + ".faa")
        instance.setup()
        return_code, database_f = instance._create_blast_db(fasta_file)
        self.assertEqual(return_code, 0)
        self.assertEqual(database_f.exists(), True)
    
    def test_reciprocal_blast(self):
        fasta_dir = Path("data/Protein_fasta_files/")
        out_dir = Path("data/")
        blast_threads = 4
        blast_evalue = 1e5
        blast_method = "diamond"
        filter = True
        ref = "GCF_000009045.1"
        instance = Orthologues_identifier(fasta_dir, out_dir, ref, blast_threads, blast_evalue, blast_method, filter)
        instance.setup()
        instance.reciprocal_blast()

    def test_parse_blast_results(self):
        fasta_dir = Path("data/Protein_fasta_files/")
        out_dir = Path("data/")
        blast_threads = 4
        blast_evalue = 1e5
        blast_method = "diamond"
        filter = True
        ref = "GCF_000009045.1"
        instance = Orthologues_identifier(fasta_dir, out_dir, ref, blast_threads, blast_evalue, blast_method, filter)
        instance.setup()
        instance.parse_blast_results()

    def test_create_orthology_matrix(self):
        fasta_dir = Path("data/Protein_fasta_files/")
        out_dir = Path("data/")
        blast_threads = 4
        blast_evalue = 1e5
        blast_method = "diamond"
        filter = True
        ref = "GCF_000009045.1"
        instance = Orthologues_identifier(fasta_dir, out_dir, ref, blast_threads, blast_evalue, blast_method, filter)
        instance.setup()
        instance.create_orthology_matrix()
        orthology_matrix_f = out_dir / "Orthologues" / ref / "OGmatrix.csv"
        orthology_matrix_df = pd.read_csv(orthology_matrix_f, sep="\t", index_col=0)
        self.assertIs(orthology_matrix_f.exists(), True)
        self.assertEqual(orthology_matrix_df.shape[0], 4237)
        self.assertEqual(orthology_matrix_df.shape[1], 2)

