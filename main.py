#from Pipeline import Core_genome_instance
#from Phylogenomics import Phylogenomics_instance
import argparse
import logging
import os
from src import args


def message(msg: str, level: int = 0):
    """
    :param msg: Message to be printed
    :param level: Level of importance of the message. 0: Info, 1: Warning, 2: Error
    :return: None
    """
    if level == 0:
        logging.info(msg)
    elif level == 1:
        logging.warning(msg)
    elif level == 2:
        logging.error(msg)
    else:
        logging.debug(msg)

def directory_setup():
        """
        Perform the necessary actions to set up the script working environment
        """
        directories = ["Blastp_output","Reciprocal_files",
                        "Blast_DB_dir","COGs","COGs_aln"]
        for directory in directories:
            tmp = self.out_dir / directory
            tmp.mkdir(parents=True,exist_ok=False)
        logging.debug(" Environment set successfully in dir %s" %(self.out_dir))

if __name__ == "__main__":
    # Setup logging
    #logfile = out_dir / "LOGFILE.log"
    logging.basicConfig(filename=logfile, format='%(asctime)s: %(message)s', level=logging.INFO)
    directory_setup()
    # 
    instance_n = Core_genome_instance(execution_params, ref)
    instance_n.reciprocal_blast()
    instance_n.parse_blast_results()
    instance_n.create_cog_matrix()
    instance_n_phylo = Phylogenomics_instance(execution_params, ref, instance_name)
    instance_n_phylo.createCOGFasta()
    instance_n_phylo.alignCOGs()
    instance_n_phylo.createSupersequence()
    instance_n_phylo.filterAln() 
    instance_n_phylo.computeTree(iqtree=True)
   
