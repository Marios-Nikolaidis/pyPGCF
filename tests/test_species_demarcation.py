import unittest
import sys
from dataclasses import dataclass


from pypgcf import species_demarcation
from pathlib import Path

@dataclass
class Data:
    def __init__(self):
        self.cores = 4
        self.kmer = 16
        self.fraglen = 3000
        self.minfraction = 0.2
        self.mcl_inflation = 2
        self.data_dir = Path(__file__).parent / "../data/genomes/"
        self.outdir = Path(__file__).parent / "test_output"
        self.outdir.mkdir(exist_ok=True)


class TestModule(unittest.TestCase):
    def test_create_input_for_fastani(self):
        data = Data()
        demarcator =  species_demarcation.SpeciesDemarcator(
            in_dir = data.data_dir,
            out_dir = data.outdir,
            fastani_cores = data.cores,
            kmer = data.kmer,
            fraglen = data.fraglen,
            minfrac = data.minfraction,
            inflation = data.mcl_inflation,
            mcl_cores = data.cores
        )
        files_for_fastani = list(data.data_dir.glob("*.fna"))
        org_list = data.outdir / "org_list.txt"
        demarcator.create_input_for_fastani(files_for_fastani, org_list)
        self.assertIs(org_list.exists(),True)

    def test_perform_fastani(self):
        data = Data()
        demarcator =  species_demarcation.SpeciesDemarcator(
            in_dir = data.data_dir,
            out_dir = data.outdir,
            fastani_cores = data.cores,
            kmer = data.kmer,
            fraglen = data.fraglen,
            minfrac = data.minfraction,
            inflation = data.mcl_inflation,
            mcl_cores = data.cores
        )
        org_list = data.outdir / "org_list.txt"
        files_for_fastani = list(data.data_dir.glob("*.fna"))
        fastani_fout = data.outdir / "FastANI.tsv"
        demarcator.create_input_for_fastani(files_for_fastani, org_list)
        demarcator.perform_fastani(org_list, fastani_fout)
        self.assertIs(fastani_fout.exists(),True)

    def prepare_input_for_mcl(self):
        data = Data()
        demarcator =  species_demarcation.SpeciesDemarcator(
            in_dir = data.data_dir,
            out_dir = data.outdir,
            fastani_cores = data.cores,
            kmer = data.kmer,
            fraglen = data.fraglen,
            minfrac = data.minfraction,
            inflation = data.mcl_inflation,
            mcl_cores = data.cores
        )
        fastani_fout = data.outdir / "FastANI.tsv"
        mcl_input = demarcator.prepare_input_for_mcl(fastani_fout)
        self.assertIs(mcl_input.exists(),True)

    def test_run_mcl(self):
        data = Data()
        demarcator =  species_demarcation.SpeciesDemarcator(
            in_dir = data.data_dir,
            out_dir = data.outdir,
            fastani_cores = data.cores,
            kmer = data.kmer,
            fraglen = data.fraglen,
            minfrac = data.minfraction,
            inflation = data.mcl_inflation,
            mcl_cores = data.cores
        )
        fastani_fout = data.outdir / "FastANI.tsv"
        mcl_input = demarcator.prepare_input_for_mcl(fastani_fout)
        mcl_output = demarcator.run_mcl(mcl_input)
        self.assertIs(mcl_output.exists(),True)

    def test_parse_mcx_output(self):
    # , fastani_from_mcl: Path) -> None:
        data = Data()
        demarcator =  species_demarcation.SpeciesDemarcator(
            in_dir = data.data_dir,
            out_dir = data.outdir,
            fastani_cores = data.cores,
            kmer = data.kmer,
            fraglen = data.fraglen,
            minfrac = data.minfraction,
            inflation = data.mcl_inflation,
            mcl_cores = data.cores
        )
        fastani_fout = data.outdir / "FastANI.tsv"
        mcl_input = demarcator.prepare_input_for_mcl(fastani_fout)
        mcl_output = demarcator.run_mcl(mcl_input)
        demarcator.parse_mcx_output(mcl_output)
        fastani_from_mcl = data.outdir / "FastANI_species_clusters.xlsx"
        self.assertIs(fastani_from_mcl.exists(),True)


