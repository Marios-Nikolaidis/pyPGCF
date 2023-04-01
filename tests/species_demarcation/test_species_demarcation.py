import unittest
import sys
from dataclasses import dataclass

sys.path.append("../species_demarcation")

from species_demarcation import (
    create_input_for_fastani, perform_fastani,
    prepare_input_for_mcl, run_mcl, parse_mcx_output
)
from pathlib import Path

@dataclass
class Data:
    def __init__(self):
        self.cores = 4
        self.kmer = 16
        self.fraglen = 3000
        self.minfraction = 0.2
        self.mcl_inflation = 2
        self.data_dir = Path(__file__).parent / "../data"


class TestModule(unittest.TestCase):
    def test_create_input_for_fastani(self):
        data = Data()
        files_for_fastani = list(data.data_dir.glob("*.fasta"))
        org_list = data.data_dir/"org_list.txt"
        create_input_for_fastani(files_for_fastani, org_list)
        self.assertIs(org_list.exists(),True)


#fastani_fout = data_dir/"FastANI.tsv"
#perform_fastani(org_list, cores, fraglen, minfraction, kmer, fastani_fout)
## Need to read the output and check that all lines are correct
#
#fastani_for_mcl = prepare_input_for_mcl(fastani_fout)
## Need to read the file and check that all lines are correct
#
#mcl_out = run_mcl(fastani_for_mcl, mcl_inflation, cores)
#
#parse_mcx_output(mcl_out)
#
