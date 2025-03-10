####################################
# Subsetting the full imputed SNPs #
####################################

# load libraries
import h5py
import pandas as pd
import numpy as np
import csv

# load data
# Please download the full_imputed_SNP from https://aragwas.1001genomes.org/#/download-center
input_file = "./geno/GENOTYPES/4.hdf5"
h5file = h5py.File(input_file,"r")

snps = h5file["snps"]
accs = h5file["accessions"]
pos = h5file["positions"]

acID = []
for i in accs:
    acID.append(int(i))

acID = np.array(acID)
gwasID = pd.read_csv("gwasIDlist.csv") # load gwasIDs

acc_list = np.array([])
for i in gwasID["GWASid"]:
    place = np.where(acID == i) # search a focal accession
    acc_list = np.append(acc_list,place[0])

acc_list = acc_list.astype(np.int64)
    
# export the subset data
snps = np.array(snps)
sub_snps = snps[:,acc_list] # note: needs to be an increasing order if hdf5
del snps
sub_snps = pd.DataFrame(sub_snps)
sub_snps.to_csv("sub_snps.csv") # export SNPs

pos = pd.DataFrame(pos)
pos.to_csv("positions.csv") # export positions

accs = pd.DataFrame(accs)
accs.to_csv("accs.csv")

