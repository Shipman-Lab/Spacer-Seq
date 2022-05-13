"""Import Modules"""
from Bio import SeqIO
import os, sys
import fuzzysearch
import itertools
import xlsxwriter, xlrd
from xlrd import open_workbook

"""Line Arguments"""
Data_Path = sys.argv[1] # path to trimmed_Results FASTQ folder
#example: C:\Users\santi.bhattaraikline\Shipman_Lab_Dev\MiSeq_Data_Dev\
#msSBK_3-145984852\Results\msSBK_3_37_trimmed_Results
sampleID = os.path.split(Data_Path)[1][0:-16]
#example: C:\Users\santi.bhattaraikline\Shipman_Lab_Dev\MiSeq_Data_Dev\
#msSBK_3-145984852\Results\msSBK_3_37_trimmed_Results
# ->msSBK_3_37
configFile = sys.argv[2]

"""Globals"""

def Target_dict(sampleID):
	from xlrd import open_workbook
	book = open_workbook(str(configFile))
	sheet = book.sheet_by_index(0)
	column0 = sheet.col_values(0)
	if column0.index(sampleID):
		sampleRow = column0.index(sampleID)
		Target_dict = {'A': [sheet.cell(sampleRow, 1).value,
			sheet.cell(sampleRow, 2).value],
			'B': [sheet.cell(sampleRow, 3).value,
			sheet.cell(sampleRow, 4).value],
			'C': [sheet.cell(sampleRow, 5).value,
			sheet.cell(sampleRow, 6).value]}
	print(Target_dict)
	for target in Target_dict:
		for target_seq in Target_dict[target]:
			if len(target_seq) != 11:
				raise Exception("Only 11 base barcodes allowed. Please check ordering excel sheet")
	return Target_dict

Target_dict = Target_dict(sampleID)
mismatch_tol = 1
count_dict = {'A': 0, 'B': 0, 'C': 0, 'N': 0}
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
		spacer_list = [sequence[(results[0].end + 5) : (results[1].start - 5)]]
		# print(last_rep)
		# print(sequence.seq)
		# print("results[0] = ")
		# print(results[0])
		# print()
		# print("results[0].end = ")
		# print(results[0].end)
		# print()
		# print("results[1] = ")
		# print(results[1])
		# print()
		# print("results[1].start = ")
		# print(results[1].start)
		# print()
		# print("spacer_list = ")
		# print(spacer_list)
		# print("results[0].end + 5 = ")
		# print(results[0].end + 5)
		# print("results[1].start - 5 = ")
		# print(results[1].start - 5)
		# print("sequence[(results[0].end + 5):(results[1].start - 5)] = ")
		# print([sequence[(results[0].end + 5):(results[1].start - 5)]])
		spacer_list.append (sequence[(results[1].end + 5) :
		(results[2].start - 5)])
		spacer_list.append (sequence[(results[2].end + 5) :
		(last_rep[len(last_rep)-1].start - 5)])
		# print(spacer_list)

	elif len(results) == 2 and len(last_rep) >= 1 and last_rep[len(last_rep)-1].start > results[len(results)-1].start:
		spacer_list = [sequence[(results[0].end + 5) : (results[1].start - 5)]]
		spacer_list.append (sequence[(results[1].end + 5) :
		(last_rep[len(last_rep)-1].start - 5)])

	elif len(results) == 1 and len(last_rep) >= 1 and last_rep[len(last_rep)-1].start > results[len(results)-1].start:
		spacer_list = [sequence[(results[0].end + 5) :
		(last_rep[len(last_rep)-1].start - 5)]]

	elif len(results) == 4:
		spacer_list = [sequence[(results[0].end + 5) : (results[1].start - 5)]]
		spacer_list.append (sequence[(results[1].end + 5) :
		(results[2].start - 5)])
		spacer_list.append (sequence[(results[2].end + 5) :
		(results[3].start - 5)])
	elif len(results) == 3:
		spacer_list = [sequence[(results[0].end + 5) : (results[1].start - 5)]]
		spacer_list.append (sequence[(results[1].end + 5) :
		(results[2].start - 5)])
	elif len(results) == 2:
		spacer_list = [sequence[(results[0].end + 5) : (results[1].start - 5)]]
	else:
		spacer_list = []
	return spacer_list

def matchesTarget(target,seq):
	if bool(fuzzysearch.find_near_matches(Target_dict[target][0],
	str(seq), max_l_dist=mismatch_tol)) \
	or bool(fuzzysearch.find_near_matches(Target_dict[target][1],
	str(seq), max_l_dist=mismatch_tol)):
		return True
	else:
		return False

def double_order(double):
	"""input is two spacers from a double expansion
		returns tuple of coded spacers, e.g. ('A','C') or ('C','N')"""
	First = 'N'
	Second = 'N'
	for target in ['A','B','C']:
		# print(target)
		# print(double[0].seq)
		if matchesTarget(target,double[0].seq):
			First = target
		if matchesTarget(target,double[1].seq):
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
		if matchesTarget(target,triplet[0].seq):
			First = target
		if matchesTarget(target,triplet[1].seq):
			Second = target
		if matchesTarget(target,triplet[2].seq):
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
double_counts = []
A_spacers = []
B_spacers = []
C_spacers = []
for seq_record in SeqIO.parse("%s/new_SPCRs_seqs.fasta" % Data_Path, "fasta"):
	total_spcrs += 1
	identity = []
	for target in Target_dict:
		if matchesTarget(target, seq_record.seq):
			identity.append(target)
	if identity == []:
		identity.append('N')
	if identity == ['A']:
		A_spacers.append(seq_record)
	if identity == ['B']:
		B_spacers.append(seq_record)
	if identity == ['C']:
		C_spacers.append(seq_record)
	# if identity != ['N']:
	# 	print(seq_record.seq,identity)
	if len(identity) > 1:
		double_counts.append(seq_record)
	if len(identity) == 1:
		count_dict[identity[0]] += 1
print('#A spacers = ',len(A_spacers))
print('#B Spacers = ',len(B_spacers))
print('#C Spacers = ',len(C_spacers))
SeqIO.write(A_spacers, "%s/A_spacers_11baseBC.fasta" % Data_Path,"fasta")
SeqIO.write(B_spacers, "%s/B_spacers_11baseBC.fasta" % Data_Path,"fasta")
SeqIO.write(C_spacers, "%s/C_spacers_11baseBC.fasta" % Data_Path,"fasta")
SeqIO.write(double_counts, "%s/double_counted_spacers_11baseBC.fasta" % Data_Path,"fasta")

for target in count_dict:
	if total_spcrs == 0:
		percent_dict[target] = 0
	else:
		percent_dict[target] = float(count_dict[target])/total_spcrs


"""Output"""
#excel file with relevant data
workbook = xlsxwriter.Workbook('%s/%s_oComp_Orders_11baseBC.xlsx' % (Data_Path,
sampleID))
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
for target in count_dict:
	worksheet.write(row,col,target,bold)
	worksheet.write(row,col+1,count_dict[target])
	row += 1

workbook.close()
