import unittest
import sys
import pandas as pd
sys.path.append("../eggNOG")
from eggNOG import (
    eggNOGInstaller,
    eggNOGRunner,
    eggNOGParser
)
from pathlib import Path

class TestModule(unittest.TestCase):

    def test_installation(self):
        resources_dir = Path("data/db")
        resources_dir.mkdir(exist_ok=True)
        installer = eggNOGInstaller(resources_dir, debug=True)
        installer.install_emapper()
        import time;
        time.sleep(10)
        installer.download_databases()

    def test_get_reference_fasta(self):
        protein_fasta_dir = Path("data/Protein_fasta_files")
        core_protein_table_f = Path("data/GCF_000009045.1_species_core.xlsx")
        out_dir = Path("data")
        cores = 4
        pident = 40
        qcov = 20
        scov = 20
        resources_dir = Path("../../resources")
        debug = True
        runner = eggNOGRunner(resources_dir, protein_fasta_dir, core_protein_table_f, out_dir,
                              cores, pident, qcov, scov, debug)
        runner.get_reference_fasta_file(self)
        should_be = Path("data/Protein_fasta_files/GCF_000009045.1.faa")
        self.assertEqual(runner.fasta_file, should_be)

    def test_get_reference_fasta_fail(self):
        protein_fasta_dir = Path("data/Protein_fasta_files")
        core_protein_table_f = Path("data/GCF_000009045.1_species_core.xlsx")
        out_dir = Path("data")
        cores = 4
        pident = 40
        qcov = 20
        scov = 20
        resources_dir = Path("../../resources")
        debug = True
        runner = eggNOGRunner(resources_dir, protein_fasta_dir, core_protein_table_f, out_dir,
                              cores, pident, qcov, scov, debug)
        runner.ref = "GCF009"
        with self.assertRaises(FileNotFoundError):
            runner.get_reference_fasta_file()

    def test_create_eggnog_cmd(self):
        protein_fasta_dir = Path("data/Protein_fasta_files")
        core_protein_table_f = Path("data/GCF_000009045.1_species_core.xlsx")
        out_dir = Path("data")
        cores = 4
        pident = 40
        qcov = 20
        scov = 20
        debug = True
        resources_dir = Path("../../resources")
        runner = eggNOGRunner(resources_dir, protein_fasta_dir, core_protein_table_f, out_dir,
                              cores, pident, qcov, scov, debug)
        runner.get_reference_fasta_file()
        expected = "emapper.py --data_dir ../../resources --cpu 4 --pident 40 --query_cov 20 --subject_cov 20 -o data/eggNOG/eggNOG_results.csv -i data/Protein_fasta_files/GCF_000009045.1.faa"
        runner.create_eggnog_cmd()
        self.assertEqual(runner.eggnog_cmd, expected)

    def test_execute_eggnog_mapper(self):
        protein_fasta_dir = Path("data/Protein_fasta_files")
        core_protein_table_f = Path("data/GCF_000009045.1_species_core.xlsx")
        out_dir = Path("data")
        cores = 4
        pident = 40
        qcov = 20
        scov = 20
        debug = True
        resources_dir = Path("../../resources")
        runner = eggNOGRunner(resources_dir, protein_fasta_dir, core_protein_table_f, out_dir,
                              cores, pident, qcov, scov, debug)

        runner.execute_eggnog_mapper() 
        self.assertEqual(runner.execute_status, 0)
        self.assertIs(runner.eggnog_raw_results_file.exists(), True)
        self.assertIs(runner.eggnog_raw_results_file.with_suffix(".csv.emmaper.hits").exists(), False)
        
    def test_parse_eggnog_results(self):
        protein_fasta_dir = Path("data/Protein_fasta_files")
        core_protein_table_f = Path("data/GCF_000009045.1_species_core.xlsx")
        out_dir = Path("data")
        parser = eggNOGParser(protein_fasta_dir, core_protein_table_f, out_dir)
        parser.parse_eggnog_results()
        self.assertEqual(len(parser.proteome_cog.keys()), 39)
        proteins = [
            "NP_387908.1", "NP_387890.1",
            "NP_387895.1", "NP_387896.1",
            "NP_387909.2", "NP_387899.1",
        ]
        values = [
            ["E"], ["F"],
            ["F"], ["F"],
            ["F"], ["F", "J"]
        ]
        for prot, val in zip(proteins, values):
            self.assertEqual(parser.proteome_cog[prot], val)

    def test_get_protein_subsets(self):
        protein_fasta_dir = Path("data/Protein_fasta_files")
        core_protein_table_f = Path("data/GCF_000009045.1_species_core.xlsx")
        out_dir = Path("data")
        parser = eggNOGParser(protein_fasta_dir, core_protein_table_f, out_dir)
        parser.get_protein_subsets()
        self.assertIn("fingerprints", parser.protein_subsets)
        self.assertEqual(len(parser.protein_subsets["fingerprints"]), 1936)
        self.assertIn("Core_100%", parser.protein_subsets)
        self.assertEqual(len(parser.protein_subsets["Core_100%"]), 4073)

    def test_gather_eggnog_results(self):
        protein_fasta_dir = Path("data/Protein_fasta_files")
        core_protein_table_f = Path("data/GCF_000009045.1_species_core.xlsx")
        out_dir = Path("data")
        parser = eggNOGParser(protein_fasta_dir, core_protein_table_f, out_dir)
        parser.gather_eggnog_results()
        categories = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "S", "T", "U", "V", "W", "Z"]
        proteome_values = [0, 0, 0, 0, 1, 5, 0, 3, 0, 3, 1, 11, 2, 0, 0, 1, 1, 12, 0, 0, 0, 0, 0]
        core_values = proteome_values[:]
        fingerprint_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 6, 0, 0, 0, 0, 0]
        sets = ["proteome", "Core_100%", "fingerprints"]
        set_values = [proteome_values, core_values, fingerprint_values]
        for protein_set, values in zip(sets, set_values):
            for category, val in zip(categories, values):
                self.assertEqual(parser.cog_categories_data[category][protein_set], val)

    
    def test_compare_sets_with_hypergeometric_test(self):
        protein_fasta_dir = Path("data/Protein_fasta_files")
        core_protein_table_f = Path("data/GCF_000009045.1_species_core.xlsx")
        out_dir = Path("data")
        parser = eggNOGParser(protein_fasta_dir, core_protein_table_f, out_dir)
        parser.gather_eggnog_results()
        parser.compare_sets_with_hypergeometric_test()
        print(parser.hypergeometric_dfs)
        # TODO: Need to manually test some values

    def test_write_hypergeometric_dfs(self):
        protein_fasta_dir = Path("data/Protein_fasta_files")
        core_protein_table_f = Path("data/GCF_000009045.1_species_core.xlsx")
        out_dir = Path("data")
        parser = eggNOGParser(protein_fasta_dir, core_protein_table_f, out_dir)
        parser.gather_eggnog_results()
        parser.compare_sets_with_hypergeometric_test()
        parser.write_hypergeometric_dfs()
        xl = pd.read_excel("data/eggNOG/eggNOG_hypergeometric.xlsx", index_col=0)
        print(xl.style.)
        #self.ass

