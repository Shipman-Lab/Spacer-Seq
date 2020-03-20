"""Import Modules"""
import os, sys, subprocess, xlsxwriter

"""File Handling"""
user_profile = os.environ ['USERPROFILE']
# File_Path = input("Trimmed FASTQ folder (run samples in here): ")
# parentFolder = os.path.split(File_Path)[0]
# FedSPCR_Path = '%s/Results' % (parentFolder)

"""Globals"""

"""Run"""
def run_FedSPCRs_v2_all():
    print()
    File_Path = input("Trimmed FASTQ folder (run samples in here): ")
    parentFolder = os.path.split(File_Path)[0]
    FedSPCR_Path = '%s/Results' % (parentFolder)
    #Create Results folder
    if not os.path.exists(FedSPCR_Path): os.makedirs(FedSPCR_Path)
    for file in os.listdir(File_Path):
        #loop through File_Path folder, running Fed_SPCRs_v2.py on all files
        #inside.
        print()
        print("Extracting spacers from: " + file)
        subprocess.call("py Fed_SPCRs_v2.py " + File_Path + " " + file + " " +
        FedSPCR_Path)

def run_oComp_Ordering():
    print()
    FedSPCR_Path = input("Results folder path: ")
    print()
    configFile = input("Config file (full path with filename): ")
    ResultsFolderList = os.listdir(FedSPCR_Path)
    for ResultsFolder in ResultsFolderList:
        print()
        print("Determining spacer order: " + ResultsFolder)
        subprocess.call("py oComp_Ordering_SBK.py " + FedSPCR_Path + "\\" +
        ResultsFolder + " " + configFile)


def inquire_about_spacers():
    print()
    reply = input("Extract spacers? (y = yes, n = no) -- ")
    if reply == "y":
        run_FedSPCRs_v2_all()
        inquire_about_ordering()
    elif reply == "n":
        inquire_about_ordering()
    else:
        inquire_about_spacers()

def inquire_about_ordering():
    print()
    reply = input("Determine spacer ordering? (y = yes, n = no) -- ")
    if reply == "y":
        print()
        reply = input("Have you completed the configuration sheet? (y = yes, "\
        "n = no) -- ")
        if reply == "y":
            run_oComp_Ordering()
        elif reply == "n":
            print()
            print("Quitting")
            quit()
        else:
            inquire_about_ordering()
    elif reply == "n":
        print()
        print("Quitting")
        quit()
    else:
        inquire_about_ordering()


inquire_about_spacers()
