import gzip
import shutil
import os, sys

user_profile = os.environ ['USERPROFILE']
File_Path = sys.argv[1]  #input("File Path (MiSeq Run \"FASTQ Generation\" Folder):")
parentFolder = os.path.split(File_Path)[0]
Data_Path = '%s/rawIndexedFASTQs' % (parentFolder)
newpath = ((r'%s') % (Data_Path))

def unzip():
# Sort through MiSeq run output folder and unzip each file inside.
# Unzipped files written to global "Data_Path"


	folder_list = os.listdir('%s' % (File_Path))
	for folder in folder_list:
		junkID = folder.split('_')[-1]
		sampleID = folder[:len(folder) - len(junkID) - 1]
		print('Unzipping ' + sampleID)
		file_name = os.listdir('%s/%s' % (File_Path,folder))[0]
		with gzip.open('%s/%s/%s' % (File_Path,folder,file_name), 'rb') as f_in:
			with open('%s/%s.fastq' % (Data_Path,sampleID), 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)


if not os.path.exists(newpath):
	os.makedirs(newpath)
	unzip()
else:
	# print("You have already unzipped these files.")
	# response = input("Do you wish to repeat? (y = yes, n = no) -- ")
	# if response == "y":
	unzip()
