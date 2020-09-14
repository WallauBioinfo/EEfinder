#!/usr/bin/python3
# -*- coding: utf-8 -*-
###############################>GENERAL-INFORMATIONS<###############################
"""
Build in Python 3.6.5+

Author:
Filipe Dezordi
zimmer.filipe@gmail.com
https://github.com/dezordi

"""
###############################>LIBRARIES<###############################
import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.style as style
###############################>ARGUMENTS<###############################
parser = argparse.ArgumentParser(description = 'Create beeswarm plots of .tax files')
parser.add_argument("-in", "--input", help=".Tax file generated by PEVEI.py", required=True)
parser.add_argument("-md","--mode",help="Create plots with family or order taxonomy? default=Family",default='Family', choices=['Family','Order'])
parser.add_argument("-st","--specifictax",help="Pass specific taxon to generate plots, e.g. Rhabdoviridae Flaviviridae",nargs='+')
args = parser.parse_args()
input_file = args.input
plot_mode = args.mode
specific_tax = args.specifictax
###############################>SNS-STYLE<###############################
sns.set(style="ticks")
style.use('seaborn-poster')
###############################>DATAFRAME<###############################
df = pd.read_csv(input_file, sep=',',header = 0)
if specific_tax != None:
    df = df.loc[df['Family'].isin(specific_tax)]
#grouping families with less than 5 elements to 'Others'
family_count = df['Family'].value_counts()
family_low = family_count[family_count <= 5]
family_low_list = family_low.index
for i in family_low_list:
    df['Family'] = df['Family'].replace({i:'Others'})
#grouping proteins terms with less than 5 elements to 'Other proteins'
protein_count = df['Protein-Product'].value_counts()
protein_low = protein_count[protein_count <= 20]
protein_low_list = protein_low.index
for i in protein_low_list:
    df['Protein-Product'] = df['Protein-Product'].replace({i:'Other proteins'})
#get number of elements by strand setnse
sense_count = df['Sense'].value_counts()
pos_count = str(sense_count['pos'])
neg_count = str(sense_count['neg'])
###############################>GENERAL-BEESWARNPLOT<###############################
beeswarn_plot =sns.swarmplot(x='EVE-length',y='Super-Kingdom', hue=plot_mode, data=df)
sns.despine(fig=None, top=True, right=True, left=True, bottom=False, offset=None, trim=False)
beeswarn_plot.set_ylabel('')
beeswarn_plot.set_yticks([])
beeswarn_plot.set_xlabel('EVE length (pb)')
beeswarn_plot.legend(loc='upper center', bbox_to_anchor=(0.7, 1.05), ncol=3, fancybox=True, prop={'size': 14},title='EVEs '+plot_mode)
plt.savefig(input_file+'.beeswarn_plot.pdf',dpi=300)
plt.clf()
###############################>SENSE-BEESWARNPLOT<###############################
beeswarn_plot =sns.swarmplot(x='Sense',y='EVE-length', hue=plot_mode, data=df)
sns.despine(fig=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
beeswarn_plot.set_xlabel('Strand Sense')
beeswarn_plot.set_ylabel('EVE length (pb)')
beeswarn_plot.set_xticklabels(['Negative (n = '+neg_count+')','Positive (n = '+pos_count+')'])
beeswarn_plot.legend(loc='upper center', ncol=2, fancybox=True, prop={'size': 12},title='EVEs '+plot_mode)
plt.savefig(input_file+'.beeswarn_plot_strand.pdf',dpi=300,bbox_inches='tight')
plt.clf()
###############################>PROTEIN-BEESWARNPLOT<###############################
beeswarn_plot =sns.swarmplot(x='EVE-length',y='Protein-Product', hue=plot_mode, data=df)
sns.despine(fig=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
beeswarn_plot.set_ylabel('Protein term')
beeswarn_plot.set_xlabel('EVE length (pb)')
beeswarn_plot.legend(loc='upper center', bbox_to_anchor=(0.7, 1.25), ncol=3, fancybox=True, prop={'size': 14},title='EVEs '+plot_mode)
#beeswarn_plot.text('center','top'   ,'Text Here', fontsize=12)
plt.savefig(input_file+'.beeswarn_plot_product.pdf',dpi=300,bbox_inches='tight')

## verificar os termos, corrigir os incorretos
## criar parametros para opções de save, e de dpi