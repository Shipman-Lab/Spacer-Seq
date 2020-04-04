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

def build_dict(sampleID):
    from xlrd import open_workbook
    book = open_workbook(str(configFile))
    sheet = book.sheet_by_index(0)
    column0 = sheet.col_values(0)
    if column0.index(sampleID):
        sampleRow = column0.index(sampleID)
        dictionary = {'LeadProx_A': sheet.cell(sampleRow, 1).value,
        'LeadProx_B': sheet.cell(sampleRow, 2).value,
        'LeadProxCore': sheet.cell(sampleRow, 3).value,
        'Mismatch_A': sheet.cell(sampleRow, 4).value,
        'Mismatch_B': sheet.cell(sampleRow, 5).value,
        'LeadDistCore': sheet.cell(sampleRow, 6).value,
        'LeadDist_A': sheet.cell(sampleRow, 7).value,
        'LeadDist_B': sheet.cell(sampleRow, 8).value,}
        return dictionary

def runVariableRegionAnalysis():
    Target_dict = build_dict(sampleID)
    ID_options = []
    ID_dict = {}
    for prod in itertools.product('ABN', repeat=3):
        ID_options.append(prod)
    for ID in ID_options:
        ID = ''.join(ID)
        ID_dict[ID] = 0
        coreSeq = Target_dict['LeadProxCore'] + 'N' + Target_dict['LeadDistCore']

    total_spcrs = 0
    target_derived = 0
    for seq_record in SeqIO.parse("%s/new_SPCRs_seqs.fasta" % Data_Path, "fasta"):
        total_spcrs += 1
        spacer = seq_record.seq
        spacerKey = []
        match = fuzzysearch.find_near_matches(coreSeq, spacer,
        max_l_dist=3)
        if match:
            target_derived += 1
            fivePrime = spacer[0:5]
            mismatch = spacer[5 + len(Target_dict['LeadProxCore'])]
            threePrime = spacer[-5:len(spacer)]
            if fuzzysearch.find_near_matches(Target_dict['LeadProx_A'], fivePrime,
            max_l_dist=1):
                spacerKey.append('A')
            elif fuzzysearch.find_near_matches(Target_dict['LeadProx_B'], fivePrime,
            max_l_dist=1):
                spacerKey.append('B')
            else:
                spacerKey.append('N')
            if mismatch == Target_dict['Mismatch_A']:
                spacerKey.append('A')
            elif mismatch == Target_dict['Mismatch_B']:
                spacerKey.append('B')
            else:
                spacerKey.append('N')
            if fuzzysearch.find_near_matches(Target_dict['LeadDist_A'], threePrime,
            max_l_dist=1):
                spacerKey.append('A')
            elif fuzzysearch.find_near_matches(Target_dict['LeadDist_B'], threePrime,
            max_l_dist=1):
                spacerKey.append('B')
            else:
                spacerKey.append('N')
            keyString = ''.join(spacerKey)
            ID_dict[keyString] += 1

    workbook = xlsxwriter.Workbook('%s/%s_VariableRegs.xlsx' % (Data_Path,
    sampleID))
    worksheet = workbook.add_worksheet()
    #Add titles (row, col: zero referenced)
    worksheetHeaders = ['Spacer Identity', 'Number', 'Sequence (if defined)']
    worksheetKey = ['Key:', 'First Letter: 5\' end source strand',
    'Second Letter: Interior mismatch source strand',
    'Third Letter: 3\' end source strand']
    spacerIDList = [key for key in ID_dict]
    spacerCountList = [ID_dict[key] for key in ID_dict]
    spacerSeqDict = {}
    for key in ID_dict:
        if 'N' in key:
            spacerSeqDict[key] = 'Undefined'
        else:
            if key[0] == 'A':
                spacerSeqDict[key] = Target_dict['LeadProx_A']
            else:
                spacerSeqDict[key] = Target_dict['LeadProx_B']
            spacerSeqDict[key] += Target_dict['LeadProxCore']
            if key[1] == 'A':
                spacerSeqDict[key] += Target_dict['Mismatch_A']
            else: spacerSeqDict[key] += Target_dict['Mismatch_B']
            spacerSeqDict[key] += Target_dict['LeadDistCore']
            if key[2] == 'A':
                spacerSeqDict[key] += Target_dict['LeadDist_A']
            else:
                spacerSeqDict[key] += Target_dict['LeadDist_B']

    worksheet.write(0, 0, sampleID)
    worksheet.write_row(1, 0, worksheetHeaders)
    worksheet.write_column(1, 5, worksheetKey)
    worksheet.write_column(2, 0, spacerIDList)
    worksheet.write_column(2, 1, spacerCountList)
    worksheet.write_column(2, 2, [spacerSeqDict[key] for key in spacerSeqDict])
    worksheet.write(31, 0, 'Total Retron-Derived:')
    worksheet.write(31, 1, target_derived)
    worksheet.write(32, 0, 'Total Nonretron-Derived:')
    worksheet.write(32, 1, total_spcrs - target_derived)
    workbook.close()

runVariableRegionAnalysis()
