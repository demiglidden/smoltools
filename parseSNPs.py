#!/usr/bin/env python

##################################################################
## This program takes a CVS of SNPS across multiple samples and ##
## adds empty space so that each column contains only 1 SNP     ##
##################################################################


#############
## IMPORTS ##
#############

import argparse
import pandas as pd
from collections import OrderedDict
import numpy as np
import datetime

#############
##  ARGS   ##
#############

def get_args():
    parser = argparse.ArgumentParser("a program to space out SNPS in a nextclade CSV, can take up to 3 CSVs using -f --file; -f2 --file2; -f3 --file3, can optionally name output using -o --output")
    parser.add_argument("-f", "--file", type=str, help="csv file to pass through program", required=True)
    parser.add_argument("-f2", "--file2", type=str, help="optional second csv file to pass through program", required=False)
    parser.add_argument("-f3", "--file3", type=str, help="optional third csv file to pass through program", required=False)
    parser.add_argument("-o", "--output", type=str, help="name of output", required=False)

    return parser.parse_args()
args = get_args()

# assign arguments to variables inside of program
CSV_IN = args.file
CSV_OUT = args.output



#######################
## READ IN, CLEAN UP ##
#######################

#columns to keep
keep_cols = ['sample','substitutions', 'deletions', 'insertions']

# read nextclade output in
nc_df1 = pd.read_csv(CSV_IN, index_col = 0)
#this renames vlaad header to match dragen, so that this script works for both dragen and VLAAD nextclade output
#rename_cols = {'seqName' : 'sample'}
#add samps as column
nc_df1.index.name = 'sample'
nc_df1.reset_index(inplace=True)
#nc_df1.rename(columns=rename_cols, inplace=True)
#drop uneccesary columns 
nc_df1 = nc_df1[nc_df1.columns.intersection(keep_cols)]

# include additional files if neccessary


if args.file2:
    if args.file3:
        CSV_IN2 = args.file2
        nc_df2 = pd.read_csv(CSV_IN2, index_col = 0)
        #nc_df2.rename(columns=rename_cols, inplace=True)
        #add samps as column
        nc_df2.index.name = 'sample'
        nc_df2.reset_index(inplace=True)
        nc_df2 = nc_df2[nc_df2.columns.intersection(keep_cols)]
        CSV_IN3 = args.file3
        nc_df3 = pd.read_csv(CSV_IN3, index_col = 0)
        #nc_df3.rename(columns=rename_cols, inplace=True)
        #add samps as column
        nc_df3.index.name = 'sample'
        nc_df3.reset_index(inplace=True)
        nc_df3 = nc_df3[nc_df3.columns.intersection(keep_cols)]
        nc_df = pd.concat([nc_df1, nc_df2, nc_df3], axis=0)
    else:
        CSV_IN2 = args.file2
        nc_df2 = pd.read_csv(CSV_IN2, index_col = 0)
        #nc_df2.rename(columns=rename_cols, inplace=True)
        #add samps as column
        nc_df2.index.name = 'sample'
        nc_df2.reset_index(inplace=True)
        nc_df2 = nc_df2[nc_df2.columns.intersection(keep_cols)]
        nc_df = pd.concat([nc_df1, nc_df2], axis=0, ignore_index=True)
else:
    nc_df = nc_df1

# change snps, and indels to strings
nc_df['substitutions'] = nc_df['substitutions'].str.strip('()').str.split(',')

# get snps into df
snps_df = nc_df.copy()[["sample","substitutions","deletions","insertions"]]
# Replace all NaN's with empty string
snps_df = snps_df.replace(np.nan, '', regex=True)


###############
## THE MAGIC ##
###############

# get snps
snps_from_df = list(snps_df['substitutions'])
# flatten list of lists
SNPs_big = []
for sublist in snps_from_df:
    for item in sublist:
        SNPs_big.append(item)
# remove duplicate snps from list
SNPS = []
for i in SNPs_big:
    if i not in SNPS:
        SNPS.append(i)

# split snps into their own cells
split_df = pd.DataFrame(snps_df['substitutions'].tolist())
# expand snps so they are each in their own column
expanded = pd.get_dummies(split_df.stack()).groupby(level=0).sum().astype(int)
expanded_with_snp_names = expanded * expanded.columns.values

# add map names back in
concat_df = pd.concat([nc_df['sample'], expanded_with_snp_names], axis=1)
snp_list = list(concat_df.columns)
# drop the sample column name
snp_list.remove('sample')
# sort snps
snp_list_sorted = sorted(snp_list, key=lambda x: int("".join([i for i in x if i.isdigit()])))
#add sample names back in
snp_list_sorted.insert(0,'sample')
#rearrange columns sorted by snp position
final_df=concat_df.reindex(columns=snp_list_sorted)


###############
## WRITE OUT ##
###############

NOW_STAMP = datetime.datetime.now().strftime("%Y%m%d")

if not args.output:
    SPACED_SNPS_FINAL = f'{NOW_STAMP}_spaced_snps.csv'
    final_df.to_csv(SPACED_SNPS_FINAL, index=False)
else:
    SPACED_SNPS_FINAL = f'{CSV_OUT}.csv'
    final_df.to_csv(SPACED_SNPS_FINAL, index=False)
