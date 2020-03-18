"""Import Modules"""
import os, sys, subprocess

"""File Handling"""
user_profile = os.environ ['USERPROFILE']
File_Path = input("Trimmed FASTQ folder (run samples in here): ")
parentFolder = os.path.split(File_Path)[0]
FedSPCR_Path = '%s/extractedSpacers' % (parentFolder)
#Create Results folder
if not os.path.exists(FedSPCR_Path): os.makedirs(FedSPCR_Path)

"""Globals"""

"""Run"""
for file in os.listdir(File_Path):
    #loop through File_Path folder, running Fed_SPCRs_v2.py on all files
    #inside.
    print()
    print("Extracting spacers from: " + file)
    subprocess.call("py Fed_SPCRs_v2.py " + File_Path + " " + file + " " +
    FedSPCR_Path)
