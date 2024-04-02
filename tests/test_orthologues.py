from dataclasses import dataclass
from pathlib import Path
import unittest
from pypgcf import orthologues


@dataclass
class Data:
    fasta_dir = Path(__file__).parent / "../data/protein/"
    out_dir = Path(__file__).parent / "../data/s/orthologues/"
    ref = "GCF_000009045"
    evalue = 1e-5
    dmnd_sensitivity = "sensitive"
    input_type = "prot"


class TestModule(unittest.TestCase):
    def test_create_blast_db(self):
        pass

    # def addCleanup(
    #     self, __function: Callable[_P, Any], *args: _P.args, **kwargs: _P.kwargs
    # ) -> None:
    #     return super().addCleanup(__function, *args, **kwargs)
