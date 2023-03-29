import unittest
import sys

sys.path.append("../smbgc")

from smbgc import smBGCInstaller, smBGCLocalRunner, smBGCParser
from pathlib import Path

class TestModule(unittest.TestCase):
    pass
    def test_installation(self):
        debug = True
        installer = smBGCInstaller(debug)
        installer.install_antismash()
        installer.install_databases()

    def test_local_runner(self):
        genome_fasta_dir = Path("data/Genome_fasta_files")
        out_dir = Path("data")
        cores = 4
        strictness = "strict"
        genefinding_tool = "prodigal"
        runner = smBGCLocalRunner(genome_fasta_dir, out_dir, cores, strictness, genefinding_tool)
        runner.analyze_genomes()


    def test_parser_gather_results(self):
        out_dir = Path("data/")
        cores = 4
        parser = smBGCParser(out_dir, cores)
        parser.gather_results()
        # Check if dataframe has correct shape
        self.assertEqual(parser.final_df.shape[0], 52)
        # Check if file was created
        xl_file = out_dir / "smBGC" / "antiSMASH_results.xlsx"
        self.assertEqual(xl_file.exists(), True)

