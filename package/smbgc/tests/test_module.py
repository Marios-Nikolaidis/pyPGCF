import unittest
import sys
from pathlib import Path
from dataclasses import dataclass

sys.path.append("../smbgc")

from smbgc import (
        smBGCWeb
)
from pathlib import Path

class TestModule(unittest.TestCase):
    def test_create_input_for_fastani(self):
        web_mod = smBGCWeb()
        input_file = Path("data/GCF_000009045.1.fna")
        tracking_dict = web_mod.post(input_file)
        print(tracking_dict)
        web_mod.keep_track_of_pending(tracking_dict)









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
