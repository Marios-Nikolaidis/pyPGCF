"""Add functional annotation to the proteins in the core genes excels"""
import pandas as pd
from collections import Counter
from Bio import SeqIO
import pathlib
from tqdm import tqdm


def initialize_ref_data(repr_file: str, sep: str = "\t") -> dict:
    """"""
    data = {}
    with open(repr_file, "r", encoding="UTF-8") as rf:
        for line in rf:
            clust, ref_org = line.rstrip().split(sep)
            data[ref_org] = clust
    return data


def initialize_categories(annot_df: pd.DataFrame, ref_data: dict) -> dict:
    """"""
    data_dict = {}
    for clust in ref_data.values():
        data_dict[clust] = {}
        for cog_category in annot_df.index:
            data_dict[clust][cog_category] = 0
    return data_dict


def get_assembly(x: str) -> str:
    return x.split("-")[1]


def get_category(x: str) -> str:
    if x == "-":
        return ["S"]
    return list(x)


print("Loading all cog categories ")
cog_annot_f = "all_COG_categories_in_dataset.tsv"
cog_annot_df = pd.read_csv(cog_annot_f, sep="\t", index_col=0)
ref_file = "Reference_organisms.txt"
ref_data = initialize_ref_data(ref_file)

print("Loading eggNOG annotation data")
annotation_f = "local_eggnog_output/rep_proteomes_eggnog.emapper.annotations"
annotation_df = pd.read_csv(annotation_f, sep="\t", index_col=0)

print("Performing transformations")
annotation_df["Protein"] = annotation_df.index
annotation_df["Assembly"] = annotation_df["Protein"].apply(lambda x: get_assembly(x))
annotation_df["COG_category"] = annotation_df["COG_category"].apply(
    lambda x: get_category(x)
)
print("Counting categories")
data_dictionary = initialize_categories(cog_annot_df, ref_data)
for assembly, assembly_df in tqdm(annotation_df.groupby("Assembly")):
    clust = ref_data[assembly]

    data_dictionary[clust] = {idx: 0 for idx in cog_annot_df.index}
    for _, row in assembly_df.iterrows():
        for category in row["COG_category"]:
            data_dictionary[clust][category] += 1


fasta_file = pathlib.Path("Representatives.faa")
record = SeqIO.parse(str(fasta_file), "fasta")
print("Adding non identified proteins as unknown")
for rec in tqdm(record):
    if rec.id not in annotation_df.index:
        _, assembly = rec.id.split("-")
        clust = ref_data[assembly]
        data_dictionary[clust]["S"] += 1

fout = "Genome_functional_categories_20221117.xlsx"
data_df = pd.DataFrame.from_dict(data_dictionary, orient="index")
data_df.to_excel(fout)


Category	Annotation
A	RNA processing and modification
B	Chromatin structure and dynamics
C	Energy production and conversion
D	Cell cycle control, cell division, chromosome partitioning
E	Amino acid transport and metabolism
F	Nucleotide transport and metabolism
G	Carbohydrate transport and metabolism
H	Coenzyme transport and metabolism
I	Lipid transport and metabolism
J	Translation, ribosomal structure and biogenesis
K	Transcription
L	Replication, recombination and repair
M	Cell wall/membrane/envelope biogenesis
N	Cell motility
O	Post-translational modification, protein turnover, and chaperones
P	Inorganic ion transport and metabolism
Q	Secondary metabolites biosynthesis, transport, and catabolism
S	Function unknown
T	Signal transduction mechanisms
U	Intracellular trafficking, secretion, and vesicular transport
V	Defense mechanisms
W	Extracellular structures
Z	Cytoskeleton
