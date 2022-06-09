# EEfinder
## Endogenous Element Finder, Alpha dev version

Run steps
===

Install dependencies with conda enviroment:

      conda env create -f env.yml
      conda activate env
      
Install diamond: https://github.com/bbuchfink/diamond/releases/tag/v2.0.15

Running:

      python EEfinder.py -in <fasta_genome> -mt <virus_proteins_table> -db <virus_db> -db2 <host_protein_db> -db3 <TE_protein_db> -p <threads>

where:

      in test_data_dev dir:
      fasta_genome -> Ae_aeg_Aag2_ctg_1913.fasta
      virus_proteins_table -> virus_subset2.csv
      virus_db -> virus_subset.fa
      host_protein_db -> filter_subset.fa
      TE_protein_db-> TEs_subset.fa
