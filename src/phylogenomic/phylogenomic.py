import pandas as pd
import pathlib
import os
import re
import logging
import multiprocessing
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


def _alignment(filelist:list):
    for file in filelist:
        cmd = " ".join([muscle_bin_path,
                        "-align -align_algo 1 -o",
                        str(outdir / "COGs_aln" / file.name),
                        str(file),
                        "> /dev/null 2>&1"
                    ])
        os.system(cmd)


class Phylogenomics_instance():

    def __init__(self, params, ref, name):
        self.system = params['sys']
        self.in_dir = params['in']
        self.out_dir = pathlib.Path(params['out']) / name
        self.ref = ref
        #self.bs_replicates = params['bootstrap']
        #self.distance = params['distance']
        self.threads = params['phylogenomic_threads']

    def createCOGFasta(self,**kwargs):
        """
        Create a fasta file for each cluster of orthologous genes
        :return: None
        """
        cog_matrix_f = self.out_dir / 'completeCOGs.csv'
        cog_matrix = pd.read_csv(cog_matrix_f, sep='\t', header=0, index_col=0)
        num_of_rows = len(cog_matrix.index)
        organisms = [self.ref[:].replace(".faa", "")]
        organisms.extend(list(cog_matrix.columns))
        for org_count, organism in enumerate(organisms):
            filename = organism
            if ".faa" not in organism:
                filename += ".faa"
            fasta_file = self.in_dir / "Fasta_files" / filename
            protein_records = SeqIO.to_dict(SeqIO.parse(open(fasta_file,"r"), format="fasta"))
            if org_count == 0:
                print(fasta_file)
                genes = list(cog_matrix.index)
            else:
                genes = list(cog_matrix[organism].values)
            for x, gene in enumerate(genes):
                info = protein_records[gene]
                info.description = ""
                info.name = ""
                info.id = info.id + "-" + organism
                fout = self.out_dir / "COGs"/ ("COG" + str(x) + ".fa")
                fout_handle = open(fout, "a")
                SeqIO.write(info, fout_handle, "fasta")
                fout_handle.close()
        logging.debug(" Created a file for each orthologous group. Dir %s"  %(self.out_dir))
    # Write presequity for dash in gene name. Need to find a way to fix this

    def alignCOGs(self):
        """
        Use muscle to align the COGs
        for muscle
        r = p / ".." / "resources" / params['sys']
        muscle_bin_path = r / "muscle3.8.31"
        """

        resources_path = pathlib.Path(__file__).parents[1] / "resources" / platform
        # THESE ARE USED BECAUSE THERE IS A pickling ERROR
        # if i put the _alignemtn inside this function
        # TODO: Try to put it inside the class?
        global muscle_bin_path 
        global outdir
        outdir = self.out_dir
        muscle_bin_path = str(r / "seaview")
        cog_files = list((self.out_dir / "COGs").glob("*"))
        logging.debug(" Initiating alignment of individual orthologous groups ")
        # pbar = tqdm(total=len(cog_files))
        # for file in cog_files:
        #   in_file = self.out_dir + "/COGs/" + file # p / 'COGs' / file
        #   # f = file + ".aln"
        #   out_file = self.out_dir + "/COGs_aln/" + file + ".aln" # p / 'COGs_aln' / f
        #   muscle_cline = muscle_bin_path + " -align -align_algo 1 -o " + out_file + " " + in_file
        #   os.system(muscle_cline)
        #   pbar.update(1)
    
        step = int(len(cog_files)/self.threads)
        if step == 0: step = 4
        files_lists = [cog_files[i:i + step] for i in range(0, len(cog_files), step)]
        with multiprocessing.Pool(processes=self.threads) as pool:
            res = pool.map(_alignment,files_lists)
        logging.debug(" Finished COG alignments in dir %s" %(self.out_dir))
    
    def createSupersequence(self):
        """
        Join the aligned COGs into a supersequence
        :return: None
        """
        def _sortAlphanum(iteratable):
            # Helper func Sort the given motif list alphanumerically :return: sorted list 
            int_convert = lambda text: int(text) if text.isdigit() else text 
            sorting_key = lambda key: [ int_convert(c) for c in re.split('([0-9]+)', key) ] 
            return sorted(iteratable, key = sorting_key)
        
        def _initSuperseqFile(superseq_file):
            """
            Initialise the supersequence file with empty sequences for each organism
            """
            cog_matrix_f = self.out_dir / 'completeCOGs.csv'
            with open(cog_matrix_f, "r") as rf:
                first_line = rf.readline().rstrip()

            organisms = [self.ref.replace(".faa", "")]
            organisms.extend(first_line.split("\t")[1:])
            seqrecord_list = []
            for org in organisms:
                record = SeqRecord(Seq(""), id = org, name = org, description="")
                seqrecord_list.append(record)
            superseq_file_handle_out = open(superseq_file, "w")
            SeqIO.write(seqrecord_list,superseq_file_handle_out, "fasta")
    

        aln_files_dir = self.out_dir / "COGs_aln"
        aln_files = aln_files_dir.glob("*")
        # aln_file = [str(f) for f in aln_files]
        # aln_files = _sortAlphanum(aln_files) 
        # To know that this is the correct order of genes
        superseq_file = self.out_dir / "supersequence.fa"
    
        _initSuperseqFile(superseq_file)
        superseq_file_handle_in = open(superseq_file, "r")
        superseq_records = SeqIO.to_dict(AlignIO.read(superseq_file_handle_in, "fasta"))
        for aln_file in aln_files:
            # fin = aln_files_dir + files
            fin_handle = open(str(aln_file), "r")
            aln_records = SeqIO.to_dict(AlignIO.read(fin_handle, "fasta"))
            aln_records = {records.split("-")[1] : aln_records[records] for records in aln_records} # Rename the keys, need to use the org name
            aln_organisms = list(aln_records.keys())
            for organism in superseq_records:
                if organism in aln_organisms:
                    superseq_records[organism].seq = superseq_records[organism].seq + aln_records[organism].seq
        superseq_seqrecords = list(superseq_records.values())
        superseq_file_handle_out = open(superseq_file, "w")
        SeqIO.write(superseq_seqrecords, superseq_file_handle_out, "fasta")
        logging.debug(" Created the supersequence file")
        return None
    
    
    def filterAln(self):
        # Filter the supersequence alignment using Gblocks with default parameters
        gblock_bin_path = os.path.join(str(pathlib.Path(__file__)),str(pathlib.Path("../resources/Gblocks_0.91b/Gblocks"))).replace(os.path.basename(__file__) + "/", "")
        cmd = " ".join([gblock_bin_path, 
                        str(self.out_dir / "supersequence.fa"),
                        "-s=y -e=-gb -p=y" 
                        # Need to check if -p=y is needed
                        # Maybe -p=n is better
                        # Allowing for gaps in the alignment?
                        ])
        os.system(cmd)
        logging.debug(" Filtered the alignment with GBlocks")
        return None
    
    def computeTree(self):
        """
        Compute the phylogenomic tree using either IQtree
        TODO: Add support for other tree building methods
        """
        supersequence_file = self.out_dir / "supersequence.fa-gb"
        resources_path = pathlib.Path(__file__).parents[1] / "resources" / platform
        iqtree = str(resources_path / "iqtree2")
        cmd = "".join([iqtree,
                        " -m TEST -merit AIC -alrt 1000 -T ",
                        str(self.threads), " -s ", str(supersequence_file)])
        method = "Maximum Likelihood"
        logging.debug(f" Initiating computation of {method} phylogenomic tree")
        os.system(cmd)
        logging.debug(" Computed Phylogenomic Tree ")
