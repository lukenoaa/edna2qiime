#!/usr/bin/env python

# Code to create a metadata file and manifest files from the NOAA_eDNA metadata template
# v2.0 (v7_20241004)
# Author: Katherine Silliman

# import libraries
import pandas as pd
import sys
import os
import argparse


# usage
usage = '''
manifest_from_noaa_template.py -s sampleMetadata_file -l libraryMetadata_file -p paired -o output_prefix

    -s, --sample_file - path name of the sampleMetadata Excel file (required)
    -l, --library_file - path name of the libraryMetadata Excel file (required)
    -p, --paired - Make a manifest file for [paired, single, both], default = paired
    -o, --output_prefix - Prefix for output files, default is ""

'''


argParser = argparse.ArgumentParser()
argParser.add_argument("-s", "--sample_file",help="path name of the sampleMetadata Excel file (required)")
argParser.add_argument("-l", "--library_file", help="path name of the libraryMetadata Excel file (required)")
argParser.add_argument("-p", "--paired", default="paired",help="Make a manifest file for [paired, single, both], default = paired")
argParser.add_argument("-o", "--output_prefix", default="",help="Prefix for output files, default is "" ")


args = argParser.parse_args()

# Load the first sheet from each Excel file
sampleMetadata = pd.read_excel(args.sample_file, sheet_name=0, comment='#',na_values=[""],index_col=None)
libraryMetadata = pd.read_excel(args.library_file, sheet_name=0, comment='#',na_values=[""],index_col=None)

# From libraryMetadata get the unique values of the 'seq_run_id' and 'assay_name' columns
seq_run_ids = libraryMetadata['seq_run_id'].unique()
assay_names = libraryMetadata['assay_name'].unique()

for runID in seq_run_ids:
    for assay in assay_names:
        # Get the subset of libraryMetadata for the current 'seq_run_id' and 'assay_name'
        libraryMetadata_subset = libraryMetadata[(libraryMetadata['seq_run_id'] == runID) & (libraryMetadata['assay_name'] == assay)]
        # Merge the subset of libraryMetadata with sampleMetadata using 'samp_name' as the key and add to list of merged DataFrames
        merged = pd.merge(libraryMetadata_subset, sampleMetadata, left_on='samp_name', right_on='samp_name')
        # If the merged DataFrame is not empty, continue
        if merged.shape[0] > 0:
            # Drop columns that contain only NaN values
            merged = merged.dropna(axis=1, how='all')
            #format metadata file and save
            meta = merged.drop(columns=[col for col in merged.columns if 'filename' in col])
            # TODO: add code to split columns in metadata file to extract numeric values
            meta.to_csv(args.output_prefix+f'-{runID}-{assay}_metadata.tsv',sep='\t',index=False)
            # user inputs fastq file path interactively
            file_path = input(f"Please provide an absolute path to fastq files for runID '{runID}' and assay '{assay}': ")
            #save manifest
            # need to split up filenames into 2 rows
            if args.paired in ["paired","both"]:
                man_pe = merged.loc[:,['samp_name','filename','filename2']]
                man_pe = pd.melt(man_pe, id_vars=['samp_name'], value_vars=['filename', 'filename2'], 
                    var_name='direction',value_name='absolute-filepath')
                man_pe['direction'] = man_pe['direction'].replace({'filename': 'forward', 'filename2': 'reverse'})
                man_pe['absolute-filepath'] = file_path + '/'+man_pe['absolute-filepath'].astype(str)
                man_pe = man_pe.rename(columns={'samp_name': 'sample-id'})
                # write paired manifest file
                man_pe = man_pe.iloc[:,[0,2,1]]
                man_pe.to_csv(args.output_prefix+f'-{runID}-{assay}_manifest_pe.csv',index=False)
                if args.paired == "both":
                    #print SE manifest
                    man_se = man_pe[man_pe['direction'] != 'reverse']
                    man_se.to_csv(args.output_prefix+f'-{runID}-{assay}_manifest_se.csv',index=False)
            #if not paired-end sequencing, only make SE manifest
            else:
                man_se = merged.loc[:,['samp_name','filename']]
                man_se = man_se.rename(columns={'samp_name': 'sample-id', 'filename': 'absolute-filepath'})
                man_se['direction'] = 'forward'
                man_se['absolute-filepath'] = file_path + '/'+man_se['absolute-filepath'].astype(str)

                man_se = man_se.iloc[:,[0,2,1]]

                man_se.to_csv(args.output_prefix+f'-{runID}-{assay}_manifest_se.csv',index=False)

print('Done!')
