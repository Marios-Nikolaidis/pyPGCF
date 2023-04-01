import os
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import logging
from pathlib import Path
import numpy as np
from typing import Tuple, Union

class Orthologues_identifier():
    def __init__(self, fasta_dir: Path, out_dir: Path, ref: Union[str, None], 
                 ref_list: Union[str, None], input_type: str,
                 cores: int, evalue: float, dmnd_sensitivity: str,
                 no_filter: bool):
        self.ref = ref
        self.ref_list = ref_list
        self.in_dir = fasta_dir.parent
        self.out_dir = out_dir / "Orthologues"
        self.blast_cores = cores
        self.blast_evalue = evalue
        self.dmnd_sensitivity = dmnd_sensitivity
        self.fasta_files = list(fasta_dir.glob("*"))
        self.input_type = input_type
        #self.threads = cores
        self.no_filter = no_filter

    def _set_directories(self, ref) -> None:
        directories = [
                self.out_dir / ref / "Blast_results",
                self.out_dir / ref / "Best_reciprocal_hits",
                self.out_dir / ref / "Blast_DB",
                ]
        for directory in directories:
            directory.mkdir(exist_ok=True, parents=True)
        

    def _get_blast_binaries(self) -> None:
        #resources_dir = Path('../../resources/')
        #makeblastdb_bin = Path(__file__).parent/resources_dir/"makeblastdb" #/ .replace(os.path.basename(__file__) + "/", "")
        #blastn_bin = Path(__file__).parent/resources_dir/"blastn" #/ .replace(os.path.basename(__file__) + "/", "")
        #diamond_bin = Path(__file__).parent/resources_dir/"diamond" #/ .replace(os.path.basename(__file__) + "/", "")
        diamond_bin = "diamond"
        blastn_bin = "blastn"
        makeblastdb_bin = "makeblastdb"

        self.blast_bin = diamond_bin
        if self.input_type == "nucl":
            self.blast_bin = blastn_bin
            self.blast_db_bin = makeblastdb_bin

    def _execute_cmd(self, cmd: str) -> int:
        """Generic function to execute a command"""
        return os.system(cmd)

    def _create_blast_db(self, ref: str, fasta_file: Path) -> Tuple[int, Path] :
        """
        Create blast database for a fasta file
        """
        database_f = self.out_dir / ref / "Blast_DB" / fasta_file.name
        cmd = f"{self.blast_bin} makedb --in {fasta_file} --quiet --db {database_f} --threads {self.blast_cores}" # DIAMOND
        if self.input_type == "nucl":
            cmd = f"{self.blast_db_bin} -in {fasta_file} -dbtype nucl -out {database_f}"
        database_f = database_f.with_suffix(".faa.dmnd")
        return self._execute_cmd(cmd), database_f

    def _create_blast_cmd(self, fasta_file: Path, database_f: Path, out_file: Path) -> str:
        # DIAMOND case first
        added_sensitivity = " --very-sensitive"
        if self.dmnd_sensitivity == "sensitive":
            added_sensitivity = " --sensitive"
        if self.dmnd_sensitivity == "mid_sensitive":
            added_sensitivity = " --mid-sensitive"
        if self.dmnd_sensitivity == "more_sensitive":
            added_sensitivity = " --more-sensitive"
        if self.dmnd_sensitivity == "ultra_sensitive":
            added_sensitivity = " --ultra-sensitive"
        cmd = f"{self.blast_bin} blastp --query {fasta_file} --quiet --db {database_f} --outfmt 6 --out {out_file} --evalue {self.blast_evalue} --threads {self.blast_cores}"
        cmd += added_sensitivity
        if self.input_type == "nucl":
            cmd = f"{self.blast_bin} -query {fasta_file} -db {database_f} -outfmt 6 -out {out_file} -evalue {self.blast_evalue} -num_threads {self.blast_cores}"
        return cmd

    def setup(self):
        if self.ref != None:
            self._set_directories(self.ref)
        if self.ref_list != None:
            list_in = open(self.ref_list, "r")
            for ref in list_in:
                ref = ref.rstrip()
                self._set_directories(ref)
            list_in.close()
        self._get_blast_binaries()

    def reciprocal_blast(self, ref:str):
        """ """
        # Create the refseq blast database
        ref_fasta = None
        for fasta_file in self.fasta_files:
            if ref in fasta_file.name:
                ref_fasta = fasta_file
        if ref_fasta == None:
            raise FileNotFoundError(f"{ref} is not in input fasta file names")

        _, ref_db = self._create_blast_db(ref, ref_fasta)

        for fasta_file in tqdm(self.fasta_files, ascii=True, leave=True, desc="Performing reciprocal BLAST"):
            if ref == fasta_file.stem:
                continue
            # Create the blast database for the fasta file
            _, database_f = self._create_blast_db(ref, fasta_file)
            # Perform the blast one way
            fout_f = self.out_dir / ref / "Blast_results" / (fasta_file.name + "_vs_reference.txt")
            fout_r = self.out_dir / ref / "Blast_results" / (fasta_file.name + "_vs_reference_reverse.txt")
            forward_blast = self._create_blast_cmd(ref_fasta, database_f, fout_f)
            reverse_blast = self._create_blast_cmd(fasta_file, ref_db, fout_r)
            self._execute_cmd(forward_blast)
            self._execute_cmd(reverse_blast)
        return ref_fasta

    def parse_blast_results(self, ref:str):
        """
        Reads the txt blast output and create the reciprocal table
        First read the files and filter the redundant information for each protein hit, keep the best hit
        Create the reciprocal blast matrices and write into files
        Input: The ref and fasta files
        Output: Best reciprocal hits per strain
        Return: None
        """
        def _get_best_subject(df: pd.DataFrame) -> pd.DataFrame:
            """
            Helper function to keep the best subject hit per query
            Input: Pandas Data frame
            Return: Filtered Pandas data frame
            """
            indeces_to_keep = []
            first_col = df.columns[0]
            for grp in df.groupby(first_col):
                indeces_to_keep.append(grp[1].index[0])
            return df.loc[indeces_to_keep].reset_index(drop=True)
    
        def _orthologue_filter(df: pd.DataFrame) -> pd.DataFrame:
            """
            Helper function to create the distribution of percent identities for all the genes in the data frame
            Filter the hits that have perc. identity lower than -2 SD
            """
            pid_std  = df.Pident.std( axis=0, numeric_only=True)
            pid_mean = df.Pident.mean(axis=0, numeric_only=True)
            drop_condition = pid_mean - (2 * pid_std)
            return df.drop(df[df.Pident < drop_condition].index)
    
        def _create_reciprocal_matrix(ref_vs_query_df, query_vs_ref_df):
            """
            Helper function which compare the two dataframes given line by line 
            to identify reciprocal hits
            Input: The one way and the reverse way dataframe
            Return a new data frame which contains only the reciprocal hits
            """
            tmp_df = pd.merge(ref_vs_query_df, query_vs_ref_df, on="RefSeq")
            df = tmp_df.drop(tmp_df[tmp_df.QuerySeq_x != tmp_df.QuerySeq_y].index)
            case1 = df[df['Pident_x'] >= df['Pident_y']]
            case1 = case1[['RefSeq', 'QuerySeq_x', 'Pident_x', 'Evalue_x']]
            case1 = case1.rename(columns={'QuerySeq_x':'QuerySeq', 'Pident_x':'Pident', 'Evalue_x':'Evalue'})
            case2 = df[df['Pident_x'] < df['Pident_y']]
            case2 = case2.rename(columns={'QuerySeq_y':'QuerySeq', 'Pident_y':'Pident', 'Evalue_y':'Evalue'})
            df = pd.concat([case1,case2])
            return df[['RefSeq', 'QuerySeq', 'Pident', 'Evalue']]
    
        logging.debug("Parsing BLAST results")
        fasta_files = self.fasta_files
        for fasta_file in tqdm(fasta_files, ascii=True, leave=True, desc="Parsing BLAST output"):
            if ref in fasta_file.stem:
                continue
            ref_vs_query_file = self.out_dir / ref / "Blast_results" / (fasta_file.name + "_vs_reference.txt")
            query_vs_ref_file = self.out_dir / ref / "Blast_results" / (fasta_file.name + "_vs_reference_reverse.txt")
    
            c1 = ["RefSeq", "QuerySeq", "Pident", "Evalue"]
            c2 = ["QuerySeq", "RefSeq", "Pident", "Evalue"]
            
            ref_vs_query_df = pd.read_csv(ref_vs_query_file,  sep="\t", names=c1, usecols=[0, 1, 2, 10] )
            query_vs_ref_df = pd.read_csv(query_vs_ref_file,  sep="\t", names=c2, usecols=[0, 1, 2, 10] )
            
            ref_vs_query_df.iloc[:, 0] = ref_vs_query_df.iloc[:, 0].astype(str)
            ref_vs_query_df.iloc[:, 1] = ref_vs_query_df.iloc[:, 1].astype(str)
    
            query_vs_ref_df.iloc[:, 0] = query_vs_ref_df.iloc[:, 0].astype(str)
            query_vs_ref_df.iloc[:, 1] = query_vs_ref_df.iloc[:, 1].astype(str)
    
            # Blast output
            # qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
            ref_vs_query_df = _get_best_subject(ref_vs_query_df)
            query_vs_ref_df = _get_best_subject(query_vs_ref_df)
            
            reciprocal_df = _create_reciprocal_matrix(ref_vs_query_df, query_vs_ref_df)
            if not self.no_filter: # No_filtering is false
                reciprocal_df = _orthologue_filter(reciprocal_df)
            # Create reciprocal files with headers
            rec_fout = self.out_dir / ref / "Best_reciprocal_hits" / (fasta_file.stem + "_BRH.txt")
            reciprocal_df.to_csv(rec_fout, sep='\t', index=False)
        logging.debug("Finished parsing results")
    
    def create_orthology_matrix(self, ref, ref_fasta):
        """
        Create the initial cog matrix with empty (NaN) vectors and replace with the true values
        :return: None
        """
        # Pass all genes of the Refseq into a list
        refgenes = SeqIO.parse(str(ref_fasta), "fasta")
        # TODO: Writing the files and then reopening them seems kind of waste of resources
        reciprocal_f_dir = self.out_dir / ref / "Best_reciprocal_hits"
        reciprocal_files = list(reciprocal_f_dir.glob("*"))
    
        def _init_orthology_matrix(refgenes, reciprocal_files):
            orthology_matrix_dict = {gene.id: {} for gene in refgenes}
            for reciprocal_file in reciprocal_files:
                query_genome = reciprocal_file.stem.replace("_BRH", "")
                for gene in orthology_matrix_dict:
                    orthology_matrix_dict[gene][query_genome] = np.nan
            return orthology_matrix_dict
    
        def _expand_orthology_matrix(init_orthology_matrix_dict, reciprocal_files):
            """
            Replace the NaN values in the Pandas dataframe with the corresponding values
            :return: expanded cogmatrix as pandas dataframe
            """
            for reciprocal_file in reciprocal_files:
                query_genome = reciprocal_file.stem.replace("_BRH", "")
                header = True
                for lines in open(reciprocal_file, mode="r"):
                    if header:
                        header=False
                        continue
                    refseq, queryseq, *_ = lines.rstrip().split("\t")
                    init_orthology_matrix_dict[refseq][query_genome] = queryseq
            init_orthology_matrix = pd.DataFrame.from_dict(init_orthology_matrix_dict, orient="index")
            init_orthology_matrix.index.name = ref
            return init_orthology_matrix
    
        init_orthology_matrix = _init_orthology_matrix(refgenes, reciprocal_files)
        orthology_matrix_df = _expand_orthology_matrix(init_orthology_matrix, reciprocal_files)
        orthology_matrix_f = self.out_dir / ref / "OGmatrix.csv"
        orthology_matrix_df.to_csv(orthology_matrix_f, sep = '\t', na_rep="X") # Contains all the COGs

    def calculate_orthologues(self):
        refs = []
        self.setup()
        if self.ref != None:
            ref = self.ref
            refs.append(ref)
        if self.ref_list != None:
            list_in = open(self.ref_list, "r")
            for ref in list_in:
                ref = ref.rstrip()
                refs.append(ref)
            list_in.close()
        for idx, ref in enumerate(refs):
            if idx > 0:
                print("-" * 200)
            print(f"Reference strain: {ref}")
            ref_fasta = self.reciprocal_blast(ref)
            self.parse_blast_results(ref)
            print("Creating orthology matrix")
            self.create_orthology_matrix(ref, ref_fasta)
            print("Done")

