"""Import Modules"""
from Bio import SeqIO
import os, sys
import fuzzysearch
import itertools
import xlsxwriter

"""Line Arguments"""
run_number = sys.argv[1] #MiSeq run number for Data_Path
condition = sys.argv[2] #file name (minus .fastq) shoud be the only argument following SPCR_blast.py when running script
running_location = sys.argv[3] #either local or orchestra

"""Globals"""
if running_location == 'local':
    user_profile = os.environ ['USERPROFILE']
    Data_Path = ('%s/Dropbox (Gladstone)/eVOLVER_and_Retro_Record/MiSeq_Data/'
        + 'msSBK_%s/msSBK_%s_%s_trimmed_Results') % (user_profile,run_number,run_number,condition) #for running locally
    # Data_Path = '%s/Dropbox/Data Analysis/MS%s/%s_Results' % (user_profile,run_number,condition) #for running locally
    Blast_Data_Path = '%s/Dropbox/Blast_Databases' % user_profile  #for running locally
elif running_location == 'orchestra_group':
    Data_Path = '/n/groups/church/Seth/MS%s_Data_Analysis/%s_Results' % (run_number,condition)  #for running on orchestra
    Blast_Data_Path = '/home/ss695/Blast_Databases' #for running on orchestra
#b3v2 is A, b3v35 is B,  psAA33 is C
Target_dict = {'A': 'GCTGTTTGTCGCTCACTGAGTCAGACTCAGTGA',
    'B': 'GAAAATGGAGAGGTTGCTGCAACCTCTCCATTT',
    'C': 'GCCCAATTTACTACTCGTTCTGGTGTTTCTCGT'}
count_dict = {'A': 0, 'B': 0, 'C': 0}
percent_dict = {}
Repeat = 'GTGTTCCCCGCGCCAGCGGGGATAAACC'
dist_repeat = 4
dist_SPCRs = 5
total_spcrs = 0

double_list = []
triple_list = []
double_options = []
triple_options = []
double_dict = {}
triple_dict = {}
for prod in itertools.product('ABCN', repeat=2):
	double_options.append(prod)
for prod in itertools.product('ABCN', repeat=3):
	triple_options.append(prod)
for double in double_options:
	double_dict[double] = 0
for triple in triple_options:
	triple_dict[triple] = 0

"""Defs"""
def get_spcrs(sequence):
    last_rep = fuzzysearch.find_near_matches(Repeat[0:15], sequence.seq, max_l_dist=3)
    results = fuzzysearch.find_near_matches(Repeat, sequence.seq, max_l_dist=6)
    if len(results) == 3 and len(last_rep) >= 1 and last_rep[len(last_rep)-1].start > results[len(results)-1].start:
        spacer_list = [sequence[results[0].end:results[1].start]]
        spacer_list.append (sequence[results[1].end:results[2].start])
        spacer_list.append (sequence[results[2].end:last_rep[len(last_rep)-1].start])
    elif len(results) == 2 and len(last_rep) >= 1 and last_rep[len(last_rep)-1].start > results[len(results)-1].start:
        spacer_list = [sequence[results[0].end:results[1].start]]
        spacer_list.append (sequence[results[1].end:last_rep[len(last_rep)-1].start])
    elif len(results) == 1 and len(last_rep) >= 1 and last_rep[len(last_rep)-1].start > results[len(results)-1].start:
        spacer_list = [sequence[results[0].end:last_rep[len(last_rep)-1].start]]

    elif len(results) == 4:
        spacer_list = [sequence[results[0].end:results[1].start]]
        spacer_list.append (sequence[results[1].end:results[2].start])
        spacer_list.append (sequence[results[2].end:results[3].start])
    elif len(results) == 3:
        spacer_list = [sequence[results[0].end:results[1].start]]
        spacer_list.append (sequence[results[1].end:results[2].start])
    elif len(results) == 2:
        spacer_list = [sequence[results[0].end:results[1].start]]
    else:
        spacer_list = []
    return spacer_list


def double_order(double):
	"""input is two spacers from a double expansion
		returns tuple of coded spacers, e.g. ('A','C') or ('C','N')"""
	First = 'N'
	Second = 'N'
	for target in ['A','B','C']:
		if len(fuzzysearch.find_near_matches(Target_dict[target],str(double[0].seq), max_l_dist=6)):
			First = target
		if len(fuzzysearch.find_near_matches(Target_dict[target],str(double[1].seq), max_l_dist=6)):
			Second = target
	order = (First,Second)
	return order

def triple_order(triplet):
	"""input is three spacers from a triple expansion
		returns tuple of coded spacers"""
	First = 'N'
	Second = 'N'
	Third = 'N'
	for target in ['A','B','C']:
		if len(fuzzysearch.find_near_matches(Target_dict[target],str(double[0].seq), max_l_dist=6)):
			First = target
		if len(fuzzysearch.find_near_matches(Target_dict[target],str(double[1].seq), max_l_dist=6)):
			Second = target
		if len(fuzzysearch.find_near_matches(Target_dict[target],str(double[2].seq), max_l_dist=6)):
			Third = target
	order = (First,Second,Third)
	return order

"""Run"""
#Pull out spacer pairs and triplets
for seq_record in SeqIO.parse("%s/double_expansion_sequences_two_read_seqs.fastq" % Data_Path, "fastq"):
	double_list.append(get_spcrs(seq_record))
for seq_record in SeqIO.parse("%s/double_expansion_sequences_three_read_seqs.fastq" % Data_Path, "fastq"):
	double_list.append(get_spcrs(seq_record))
for seq_record in SeqIO.parse("%s/triple_expansion_sequences_three_read_seqs.fastq" % Data_Path, "fastq"):
	triple_list.append(get_spcrs(seq_record))

#get order for doubles and triplets
for double in double_list:
	double_dict[double_order(double)] += 1
for triplet in triple_list:
	triple_dict[triple_order(triplet)] += 1

#get spacer counts for each target
for seq_record in SeqIO.parse("%s/new_SPCRs_seqs.fasta" % Data_Path, "fasta"):
	total_spcrs += 1
	for key in Target_dict:
		if len(fuzzysearch.find_near_matches(Target_dict[key], seq_record.seq, max_l_dist=6)):
			count_dict[key] += 1
count_dict['N'] = total_spcrs-sum(count_dict.values())
for key in count_dict:
    if total_spcrs == 0:
        percent_dict[key] = 0
    else:
        percent_dict[key] = float(count_dict[key])/total_spcrs


"""Output"""
#excel file with relevant data
workbook = xlsxwriter.Workbook('%s/oComp_Orders.xlsx' % Data_Path)
worksheet = workbook.add_worksheet()
bold = workbook.add_format({'bold': True})
#Add titles (row, col: zero referenced)
worksheet.write(0,0,'Doubles',bold)
worksheet.write(0,2,'expected')
worksheet.write(0,3,'Triples',bold)
worksheet.write(0,6,'From Doubles')
worksheet.write(0,7,'First',bold)
worksheet.write(0,8,'Second',bold)
worksheet.write(5,6,'Total Spacers')
worksheet.write(1,6,'A')
worksheet.write(2,6,'B')
worksheet.write(3,6,'C')

#Add data
row = 1
col = 0
for double in double_options:
	worksheet.write(row,col,''.join(double))
	worksheet.write(row,col+1,double_dict[double])
	worksheet.write(row,col+2,sum(double_dict.values())*percent_dict[double[0]]*percent_dict[double[1]])
	row += 1
row = 1
col = 3
for triple in triple_options:
	worksheet.write(row,col,''.join(triple))
	worksheet.write(row,col+1,triple_dict[triple])
	row += 1
worksheet.write(1,7,double_dict[('A','B')]+double_dict[('A','C')]+double_dict[('A','N')])
worksheet.write(1,8,double_dict[('B','A')]+double_dict[('C','A')]+double_dict[('N','A')])
worksheet.write(2,7,double_dict[('B','A')]+double_dict[('B','C')]+double_dict[('B','N')])
worksheet.write(2,8,double_dict[('A','B')]+double_dict[('C','B')]+double_dict[('N','B')])
worksheet.write(3,7,double_dict[('C','A')]+double_dict[('C','B')]+double_dict[('C','N')])
worksheet.write(3,8,double_dict[('A','C')]+double_dict[('B','C')]+double_dict[('N','C')])
row = 6
col = 6
for key in count_dict:
	worksheet.write(row,col,key,bold)
	worksheet.write(row,col+1,count_dict[key])
	row += 1

workbook.close()