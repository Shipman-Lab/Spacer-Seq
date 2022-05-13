"""Import Modules"""
import os, sys, subprocess, xlsxwriter, platform

"""File Handling"""
File_Path = sys.argv[1] #input("Trimmed FASTQ folder (run samples in here): ")
parentFolder = os.path.split(File_Path)[0]
FedSPCR_Path = '%s/Results' % (parentFolder)
# FedSPCR_Path = sys.argv[2] #input("Results folder path: ")
logic = sys.argv[2]
OconfigFile = sys.argv[3] #input("Config file (full path with filename): ")
VarConfigFile = sys.argv[4] #input("Config file for variable region analysis")

"""Globals"""
sysPlatform = platform.system()

"""Run"""
def run_FedSPCRs_v2_all(File_Path):
    print()
    # parentFolder = os.path.split(File_Path)[0]
    # FedSPCR_Path = '%s/Results' % (parentFolder)
    #Create Results folder
    if not os.path.exists(FedSPCR_Path): os.makedirs(FedSPCR_Path)
    for file in os.listdir(File_Path):
        #loop through File_Path folder, running Fed_SPCRs_v2.py on all files
        #inside.
        print()
        print("Extracting spacers from: " + file)
        if sysPlatform == 'Windows':
            subprocess.call(['py', 'Fed_SPCRs_v2.py', File_Path, file, \
            FedSPCR_Path])
        else:
            subprocess.call(['python3', 'Fed_SPCRs_v2.py', File_Path, file, \
            FedSPCR_Path])

def run_oComp_Ordering(FedSPCR_Path):
    print()
    print()
    ResultsFolderList = os.listdir(FedSPCR_Path)
    for ResultsFolder in ResultsFolderList:
        print()
        print("Determining spacer order: " + ResultsFolder)
        if sysPlatform == 'Windows':
            subprocess.call(['py', 'oComp_Ordering_SBK.py', FedSPCR_Path + "/" +
            ResultsFolder, OconfigFile])
        else:
            subprocess.call(['python3', 'oComp_Ordering_SBK.py', FedSPCR_Path + "/" +
            ResultsFolder, OconfigFile])

def run_oComp_Ordering_11baseBC(FedSPCR_Path):
    print()
    print()
    ResultsFolderList = os.listdir(FedSPCR_Path)
    for ResultsFolder in ResultsFolderList:
        print()
        print("Determining spacer order (11 base barcode): " + ResultsFolder)
        if sysPlatform == 'Windows':
            subprocess.call(['py', 'oComp_Ordering_11baseBC.py', FedSPCR_Path + "/" +
            ResultsFolder, OconfigFile])
        else:
            subprocess.call(['python3', 'oComp_Ordering_11baseBC.py', FedSPCR_Path + "/" +
            ResultsFolder, OconfigFile])

def run_VariableRegion(FedSPCR_Path):
    print()
    print()
    ResultsFolderList = os.listdir(FedSPCR_Path)
    for ResultsFolder in ResultsFolderList:
        print()
        print("Analyzing variable regions: " + ResultsFolder)
        if sysPlatform == 'Windows':
            subprocess.call(['py', 'VariableRegs.py', FedSPCR_Path + "/" +
            ResultsFolder, VarConfigFile])
        else:
            subprocess.call(['python3', 'VariableRegs.py', FedSPCR_Path + "/" +
            ResultsFolder, VarConfigFile])

if logic == 'e':
    run_FedSPCRs_v2_all(File_Path)
if logic == 'o':
    run_oComp_Ordering(FedSPCR_Path)
if logic == 'o11':
    run_oComp_Ordering_11baseBC(FedSPCR_Path)
if logic == 'v':
    run_VariableRegion(FedSPCR_Path)
else:
    print('Quitting')
# def inquire_about_spacers():
#     print()
#     reply = input("Extract spacers? (y = yes, n = no) -- ")
#     if reply == "y":
#         run_FedSPCRs_v2_all()
#         inquire_about_ordering()
#     elif reply == "n":
#         inquire_about_ordering()
#     else:
#         inquire_about_spacers()
#
# def inquire_about_ordering():
#     print()
#     reply = input("Determine spacer ordering? (y = yes, n = no) -- ")
#     if reply == "y":
#         print()
#         reply = input("Have you completed the configuration sheet? (y = yes, "\
#         "n = no) -- ")
#         if reply == "y":
#             run_oComp_Ordering()
#         elif reply == "n":
#             print()
#             print("Quitting")
#             quit()
#         else:
#             inquire_about_ordering()
#     elif reply == "n":
#         print()
#         print("Quitting")
#         quit()
#     else:
#         inquire_about_ordering()
#
#
# inquire_about_spacers()
