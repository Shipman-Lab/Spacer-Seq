import gzip
import shutil
import os, sys

user_profile = os.environ ['USERPROFILE']
File_Path = input("File Path (MiSeq Run Folder):")
runName = input("Run Name (No Spaces):")
Data_Path = '%s/%s_rawIndexedFASTQs' % (File_Path, runName)

# Create Data folder
newpath = ((r'%s') % (Data_Path))
if not os.path.exists(newpath): os.makedirs(newpath)

interior_folder = os.listdir(File_Path)[0]

folder_list = os.listdir('%s/%s' % (File_Path,interior_folder))
for folder in folder_list:
	junkID = folder.split('_')[-1]
	sampleID = folder[:len(folder) - len(junkID) - 1]
	print('Unzipping ' + sampleID)
	file_name = os.listdir('%s/%s/%s' % (File_Path,interior_folder,folder))[0]
	with gzip.open('%s/%s/%s/%s' % (File_Path,interior_folder,folder,file_name), 'rb') as f_in:
		with open('%s/%s.fastq' % (Data_Path,sampleID), 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)











# with gzip.open('file.txt.gz', 'rb') as f_in:
#     with open('file.txt', 'wb') as f_out:
#         shutil.copyfileobj(f_in, f_out)
