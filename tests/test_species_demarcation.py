from dataclasses import dataclass
from pathlib import Path
import unittest
from pypgcf import species_demarcation


@dataclass
class Data:
     cores = 4
     kmer = 16
     fraglen = 3000
     minfraction = 0.2
     mcl_inflation = 2
     data_dir = Path(__file__).parent / "../data/genomes/"
     outdir = Path(__file__).parent / "../data/species_demarcation/"
     outdir.mkdir(exist_ok=True)


class TestModule(unittest.TestCase):

    def test_assign_species(self):
        data = Data()
        demarcator = species_demarcation.SpeciesDemarcator(
            in_dir=data.data_dir,
            out_dir=data.outdir,
            fastani_cores=data.cores,
            kmer=data.kmer,
            fraglen=data.fraglen,
            minfrac=data.minfraction,
            inflation=data.mcl_inflation,
            mcl_cores=data.cores,
            debug=True,
        )
        demarcator.assign_species()

    def test_create_input_for_fastani(self):
        data = Data()
        demarcator = species_demarcation.SpeciesDemarcator(
            in_dir=data.data_dir,
            out_dir=data.outdir,
            fastani_cores=data.cores,
            kmer=data.kmer,
            fraglen=data.fraglen,
            minfrac=data.minfraction,
            inflation=data.mcl_inflation,
            mcl_cores=data.cores,
        )
        files_for_fastani = list(data.data_dir.glob("*.fna"))
        org_list = data.outdir / "org_list.txt"
        demarcator.create_input_for_fastani(files_for_fastani, org_list)
        self.assertIs(org_list.exists(), True)

    def test_perform_fastani(self):
        data = Data()
        demarcator = species_demarcation.SpeciesDemarcator(
            in_dir=data.data_dir,
            out_dir=data.outdir,
            fastani_cores=data.cores,
            kmer=data.kmer,
            fraglen=data.fraglen,
            minfrac=data.minfraction,
            inflation=data.mcl_inflation,
            mcl_cores=data.cores,
        )
        org_list = data.outdir / "org_list.txt"
        files_for_fastani = list(data.data_dir.glob("*.fna"))
        fastani_fout = data.outdir / "FastANI.tsv"
        demarcator.create_input_for_fastani(files_for_fastani, org_list)
        demarcator.perform_fastani(org_list, fastani_fout)
        self.assertIs(fastani_fout.exists(), True)

    def test_prepare_input_for_mcl(self):
        data = Data()
        demarcator = species_demarcation.SpeciesDemarcator(
            in_dir=data.data_dir,
            out_dir=data.outdir,
            fastani_cores=data.cores,
            kmer=data.kmer,
            fraglen=data.fraglen,
            minfrac=data.minfraction,
            inflation=data.mcl_inflation,
            mcl_cores=data.cores,
        )
        fastani_fout = data.outdir / "FastANI.tsv"
        mcl_input = demarcator.prepare_input_for_mcl(fastani_fout)
        self.assertIs(mcl_input.exists(), True)

    def test_run_mcl(self):
        data = Data()
        demarcator = species_demarcation.SpeciesDemarcator(
            in_dir=data.data_dir,
            out_dir=data.outdir,
            fastani_cores=data.cores,
            kmer=data.kmer,
            fraglen=data.fraglen,
            minfrac=data.minfraction,
            inflation=data.mcl_inflation,
            mcl_cores=data.cores,
            debug=True
        )
        fastani_fout = data.outdir / "FastANI.tsv"
        mcl_input = demarcator.prepare_input_for_mcl(fastani_fout)
        mcl_output = demarcator.run_mcl(mcl_input)
        self.assertIs(mcl_output.exists(), True)

    def test_parse_mcx_output(self):
        data = Data()
        demarcator = species_demarcation.SpeciesDemarcator(
            in_dir=data.data_dir,
            out_dir=data.outdir,
            fastani_cores=data.cores,
            kmer=data.kmer,
            fraglen=data.fraglen,
            minfrac=data.minfraction,
            inflation=data.mcl_inflation,
            mcl_cores=data.cores,
        )
        fastani_fout = data.outdir / "FastANI.tsv"
        mcl_input = demarcator.prepare_input_for_mcl(fastani_fout)
        mcl_output = demarcator.run_mcl(mcl_input)
        demarcator.parse_mcx_output(mcl_output)
        fastani_from_mcl = data.outdir / "FastANI_species_clusters.xlsx"
        self.assertIs(fastani_from_mcl.exists(), True)

    def test_read_fastani_output(self):
        ...
        # data = Data()
        # demarcator = species_demarcation.SpeciesDemarcator(
        #     in_dir=data.data_dir,
        #     out_dir=data.outdir,
        #     fastani_cores=data.cores,
        #     kmer=data.kmer,
        #     fraglen=data.fraglen,
        #     minfrac=data.minfraction,
        #     inflation=data.mcl_inflation,
        #     mcl_cores=data.cores,
        # )
        # fastani_fout = data.outdir / "FastANI.tsv"
        # mcl_input = demarcator.prepare_input_for_mcl(fastani_fout)
        # mcl_output = demarcator.run_mcl(mcl_input)
        # demarcator.parse_mcx_output(mcl_output)
        # fastani_from_mcl = data.outdir / "FastANI_species_clusters.xlsx"
        # df = pd.read_excel(fastani_from_mcl, index_col=0)
        # test_df = pd.DataFrame.from_dict(
        #     {
        #         "GCA_000769555": {"FastANI_species": "C0"},
        #         "GCA_002220285": {"FastANI_species": "C1"},
        #         "GCA_000009045": {"FastANI_species": "C2"},
        #     },
        #     orient="index",
        # )
        # print(test_df)
        # print(df)
        # self.assertTrue(df.equals(test_df))
        # return None
