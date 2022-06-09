# EEfinder
## Endogenous Element Finder

Run steps
===

To use conda enviroment:
> conda env create -f env.yml

> conda activate env

Running:
> python EEfinder.py -in <fasta_genome> -mt <virus_proteins_table> -db <virus_db> -db2 <heatshock_db> -db3 <TE_db> -p <threads>
    
      in test_files_dev
      fasta_genome -> Ae_aeg_Aag2_ctg_1913.fasta
      virus_proteins_table -> virus_subset2.csv
      virus_db -> virus_subset.fa
      TE_db -> TEs_subset.fa
