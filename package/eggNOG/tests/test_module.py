import unittest
import sys
from dataclasses import dataclass

sys.path.append("../eggNOG")

from eggNOG import (
        eggNOGInstaller
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

#prepare_input_for_mcl

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
