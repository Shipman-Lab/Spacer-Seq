import gzip
import shutil
import os, sys

user_profile = os.environ ['USERPROFILE']
File_Path = '%s/Dropbox (Gladstone)/eVOLVER_and_Retro_Record/MiSeq_Data/msSBK_2-141989848'  % user_profile
Data_Path = '%s/Dropbox (Gladstone)/eVOLVER_and_Retro_Record/MiSeq_Data/msSBK_2'  % user_profile

#Create Data folder
newpath = ((r'%s') % (Data_Path))
if not os.path.exists(newpath): os.makedirs(newpath)

interior_folder = os.listdir(File_Path)[0]

folder_list = os.listdir('%s/%s' % (File_Path,interior_folder))
for folder in folder_list:
	number = folder.split('_')[2]
	file_name = os.listdir('%s/%s/%s' % (File_Path,interior_folder,folder))[0]
	with gzip.open('%s/%s/%s/%s' % (File_Path,interior_folder,folder,file_name), 'rb') as f_in:
		with open('%s/msSBK_2_%s.fastq' % (Data_Path,number), 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)











# with gzip.open('file.txt.gz', 'rb') as f_in:
#     with open('file.txt', 'wb') as f_out:
#         shutil.copyfileobj(f_in, f_out)
