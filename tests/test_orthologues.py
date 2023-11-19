from collections.abc import Callable
from typing import Any
import unittest
from pypgcf import orthologues


class TestModule(unittest.TestCase):
    def setup(self):
        protein_names = []
        protein_sequences = []
        nucleotide_names = []
        nucleotide_sequences = []
        wf = open("test_prot.fa", "w")
        for protein, sequence in zip(protein_names, protein_sequences):
            to_write = f">{protein}\n{sequence}"
            wf.write(to_write)

    def test_create_blast_db(self):
        pass

    # def addCleanup(
    #     self, __function: Callable[_P, Any], *args: _P.args, **kwargs: _P.kwargs
    # ) -> None:
    #     return super().addCleanup(__function, *args, **kwargs)
