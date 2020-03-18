import os, sys, subprocess

user_profile = os.environ ['USERPROFILE']
inFolderPath = input("Paste file path here (folder containing raw FASTQs): ")
parentFolder = os.path.split(inFolderPath)[0]
outFolderPath = parentFolder + '\\trimmedFASTQs'

def trim():
    for file in os.listdir(inFolderPath):
        print("Trimming " + file)
        inFilePath = inFolderPath + '\\' + file
        sampleID = os.path.splitext(file)[0]
        outFilePath = outFolderPath + '\\' + sampleID + '_trimmed.fastq'
        subprocess.run("java -jar " + os.getcwd() + '\\trimmomatic-0.39.jar SE ' + inFilePath + ' ' + outFilePath + " LEADING:3 TRAILING:3" + " SLIDINGWINDOW:6:25 MINLEN:60")
        print()
if not os.path.exists(outFolderPath):
    os.makedirs(outFolderPath)
    trim()
else:
    print("You have already trimmed these files.")
    response = input("Do you wish to proceed? (y = yes, n = no) -- ")
    if response == "y":
        trim()
