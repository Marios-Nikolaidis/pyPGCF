# Default configurations for the modules
# By M. Nikolaidis

from sys import platform as platform
from multiprocessing import cpu_count

system = platform
system_cores = cpu_count()

## 1. Species demarcation configs
species_demarcation_cores = system_cores - 2
species_demarcation_kmer = 16
species_demarcation_fraglen = 3000
species_demarcation_minfraction = 0.2
species_demarcation_mcl_inflation = 2

## 2. Orthologues configs
orthologues_cores = 6
orthologues_evalue = 1e-5
orthologues_dmnd_sensitivity = "very-sensitive"

## 3. Core configs
core_core_perc = 100

## 4. Phylogenomic configs
phylogenomic_cores = system_cores - 2

## 5. eggNOG configs
eggnog_cores = system_cores - 2
eggnog_pident = 40
eggnog_scov = 20
eggnog_qcov = 20

## 6. smBGC configs
smbgc_cores = 6
smbgc_strictness = "strict"
smbgc_genefinding_tool = "prodigal"
smbgc_valid_strictness = ["loose", "relaxed", "strict"]
smbgc_genefinding_tools = ["prodigal"]
