from pathlib import Path

import pandas as pd
import unittest
import sys

sys.path.append("../core")
from core import Core_identifier

class TestModule(unittest.TestCase):
    def test_split_genomes_into_groups(self):
        indir = Path("data")
        og_matrix = indir / "OGmatrix.csv"
        species_file = indir / "FastANI_species_clusters.xlsx"
        out_dir = Path("data")
        core_perc = 100
        instance = Core_identifier(og_matrix, species_file, core_perc, out_dir)
        instance.split_genomes_into_groups()
        self.assertEqual(len(instance.group_orgs), 1)
        self.assertNotIn("GCF_000009045.1", instance.group_orgs)
        self.assertIn("GCF_019704415.1", instance.group_orgs)
        self.assertEqual("GCF_000009045.1", instance.ref)
        self.assertEqual(len(instance.non_group_orgs), 1)

    def test_calculate_protein_presence(self):
        indir = Path("data")
        og_matrix = indir / "OGmatrix.csv"
        species_file = indir / "FastANI_species_clusters.xlsx"
        out_dir = Path("data")
        core_perc = 100
        instance = Core_identifier(og_matrix, species_file, core_perc, out_dir)
        instance.split_genomes_into_groups()
        presence_df = instance.calculate_protein_presence()
        num_core = len(presence_df[presence_df["Group orthologues %"] == 100].index.tolist())
        self.assertEqual(num_core, 4073)
        specific_prot = "NP_387886.2"
        print(instance.orthology_df.loc[specific_prot])
        print(instance.group_orgs)
        specific_prot_presence = presence_df.loc[specific_prot, "Group orthologues %"]
        specific_prot_presence_outside = presence_df.loc[specific_prot, "Non group orthologues %"]
        print(presence_df.loc[specific_prot])
        self.assertEqual(specific_prot_presence, 100)
        self.assertEqual(specific_prot_presence_outside, 0)

    # def test_identify_core_proteins(self):
    #     indir = Path("data")
    #     og_matrix = indir / "OGmatrix.csv"
    #     species_file = indir / "FastANI_species_clusters.xlsx"
    #     out_dir = Path("data")
    #     core_perc = 100
    #     instance = Core_identifier(og_matrix, species_file, core_perc, out_dir)
    #     instance.split_genomes_into_groups()
    #     presence_df = instance.calculate_protein_presence()
    #     core_protein_data = instance.identify_core_proteins(presence_df)
    #

    # def test_split_genomes_into_groups_various_perc(self):
    #     indir = Path("data")
    #     og_matrix = indir / "OGmatrix.csv"
    #     species_file = indir / "FastANI_species_clusters.xlsx"
    #     out_dir = Path("data")
    #     core_percents = [100,98,95,90]
    #     for core_perc in core_percents:
    #         instance = Core_identifier(og_matrix, species_file, core_perc, out_dir)
    #         instance.calculate_core()

if __name__ == "__main__":
    test = TestModule()
    test.test_identify_core_proteins()
