
"""Simulation of CRISPR spacer acquisition by cells
Santi Bhattarai-Kline 6/23/21"""

import numpy as np
import matplotlib.pyplot as plt
import random
import itertools
import xlsxwriter


"""Globals"""
# #method of simulating integrations. specify 'full_interval' or 'subinterval'
# method = 'full_interval'
# #number of cells in the simulated experiment
# population_size = 1000000
#
# #rates of acquisition of spacer types A, B, and N with and without inducers
# #present. rates should be given in units of integrations per hour per cell
# rate_A_on = 6.46E-5
# rate_A_off = 1.21E-5
# rate_B_on = 6.47E-5
# rate_B_off = 4.37E-5
# rate_N = 1.41E-3
# #times of expression of spacer should be given in units of hours
# time_of_expression_A = 24
# time_of_expression_B = 24
#
# #unit of time (in hours) over which integration probabilities are assessed.
# #the user assumes that over the time specified, a single cell may receive no
# #more than a single integration
# subinterval_length = 0.1

"""Classes"""

class Cell:

    def __init__(self, name):
        self.name = name
        #arrays are arranged in temporal order with oldest acquisition in
        #position array[0]
        self.array = ''

    def integrate_spacer(self, spacer):
        self.array += spacer

"""Functions"""

def create_cells(number):
    list_of_cells = list()
    for i in range(number):
        list_of_cells.append(Cell(str(i)))
    return list_of_cells


def run_single_experiment(method='full_interval', population_size=1000000,
rate_A_on=6.46E-5, rate_A_off=1.21E-5, rate_B_on=6.47E-5, rate_B_off=4.37E-5,
rate_N=1.41E-3, time_of_expression_A=24, time_of_expression_B=24,
order_of_expression='A_before_B', subinterval_length = 0.1):
#'method' defines method of simulating integrations. specify 'full_interval' or 'subinterval'
#'population' defines number of cells in the simulated experiment
#'rates' of acquisition of spacer types A, B, and N with and without inducers
#present. rates should be given in units of integrations per hour per cell
#times of expression of spacer should be given in units of hours.
#'order_of_expression': specify 'A_before_B' or 'B_before_A'
#'subinterval_length' defines unit of time (in hours) over which integration
#probabilities are assessed. the user assumes that over the time specified, a
#single cell may receive no more than a single integration
    if order_of_expression == 'A_before_B':
        retron_1 = 'A'
        retron_2 = 'B'
        rate_1_on = rate_A_on
        rate_1_off = rate_A_off
        rate_2_on = rate_B_on
        rate_2_off = rate_B_off
        time_of_expression_1 = time_of_expression_A
        time_of_expression_2 = time_of_expression_B
    elif order_of_expression == 'B_before_A':
        retron_1 = 'B'
        retron_2 = 'A'
        rate_1_on = rate_B_on
        rate_1_off = rate_B_off
        rate_2_on = rate_A_on
        rate_2_off = rate_A_off
        time_of_expression_1 = time_of_expression_B
        time_of_expression_2 = time_of_expression_A

    cell_list = create_cells(population_size)
    if method == 'subinterval':
        #first epoch
        for i in range(int(time_of_expression_1//subinterval_length)):
            for cell in cell_list:
                integrations = ''
                integrations += retron_1 * int(np.random.poisson((rate_1_on *
                subinterval_length), 1))
                integrations += retron_2 * int(np.random.poisson((rate_2_off *
                subinterval_length), 1))
                integrations += 'N' * int(np.random.poisson((rate_N *
                subinterval_length), 1))
                if len(integrations) > 1:
                    cell.integrate_spacer(random.choice(integrations))
                else:
                    cell.integrate_spacer(integrations)
        #second epoch
        for i in range(int(time_of_expression_2//subinterval_length)):
            for cell in cell_list:
                integrations = ''
                integrations += retron_1 * int(np.random.poisson((rate_1_off *
                subinterval_length), 1))
                integrations += retron_2 * int(np.random.poisson((rate_2_on *
                subinterval_length), 1))
                integrations += 'N' * int(np.random.poisson((rate_N *
                subinterval_length), 1))
                if len(integrations) > 1:
                    cell.integrate_spacer(random.choice(integrations))
                else:
                    cell.integrate_spacer(integrations)

    if method == 'full_interval':
        for cell in cell_list:
            integrations = ''
            integrations += retron_1 * int(np.random.poisson((rate_1_on *
            time_of_expression_1), 1))
            integrations += retron_2 * int(np.random.poisson((rate_2_off *
            time_of_expression_1), 1))
            integrations += 'N' * int(np.random.poisson((rate_N *
            time_of_expression_1), 1))
            shuffled_integrations = ''.join(random.sample(integrations,
            len(integrations)))
            cell.integrate_spacer(shuffled_integrations)

        for cell in cell_list:
            integrations = ''
            integrations += retron_1 * int(np.random.poisson((rate_1_off *
            time_of_expression_2), 1))
            integrations += retron_2 * int(np.random.poisson((rate_2_on *
            time_of_expression_2), 1))
            integrations += 'N' * int(np.random.poisson((rate_N *
            time_of_expression_2), 1))
            shuffled_integrations = ''.join(random.sample(integrations,
            len(integrations)))
            cell.integrate_spacer(shuffled_integrations)

    double_dict = {}
    for prod in itertools.product('ABN', repeat=2):
    	double_dict[''.join(prod)] = 0

    for cell in cell_list:
        #trim array to 3 spacers to mimic sequence length limit=
        if len(cell.array) > 3:
            cell.array = cell.array[-3:]
        #sort triple acquisition into doubles and tally
        if len(cell.array) == 3:
            for i in range(3):
                extracted_double = cell.array[:i] + cell.array[(i + 1):]
                double_dict[extracted_double] += 1
        #tally double acquisition
        elif len(cell.array) == 2:
            double_dict[cell.array] += 1

    # print()
    print(double_dict)

    max_composite_score = 0
    non_normalized_composite = 0

    if (double_dict['AB'] + double_dict['BA']) > 0:
        AB_score = ((double_dict['AB'] - double_dict['BA'])
        / (double_dict['AB'] + double_dict['BA']))
        max_composite_score += 1
        non_normalized_composite += AB_score
    else:
        AB_score = None

    if (double_dict['AN'] + double_dict['NA']) > 0:
        AN_score = ((double_dict['AN'] - double_dict['NA'])
        / (double_dict['AN'] + double_dict['NA']))
        max_composite_score += 0.5
        non_normalized_composite += AN_score
    else:
        AN_score = None

    if (double_dict['NB'] + double_dict['BN']) > 0:
        BN_score = ((double_dict['NB'] - double_dict['BN'])
        / (double_dict['NB'] + double_dict['BN']))
        max_composite_score += 0.5
        non_normalized_composite += BN_score
    else:
        BN_score = None

    if max_composite_score > 0:
        composite_score = non_normalized_composite / max_composite_score
    else:
        composite_score = None

    score_dict = {'AN_score': AN_score, 'BN_score': BN_score, 'AB_score': AB_score,
    'composite_score': composite_score}
    # print(score_dict)
    # print(non_normalized_composite, max_composite_score)
    return(score_dict)


"""Execute"""

data_AB = []
data_BA = []
for i in range(24):
	data_AB.append([run_single_experiment(population_size = (1000000),
	time_of_expression_A=(i + 1)*(72/24), time_of_expression_B=(i + 1)*(72/24),
	order_of_expression='A_before_B')['composite_score'] for j in range(100)])
	data_BA.append([run_single_experiment(population_size = (1000000),
	time_of_expression_A=(i + 1)*(72/24), time_of_expression_B=(i + 1)*(72/24),
	order_of_expression='B_before_A')['composite_score'] for j in range(100)])

print('data_AB', data_AB)
print('data_BA', data_BA)

workbook_134 = xlsxwriter.Workbook('134_sim_time_vary_comp2.xlsx')
sheet1_134_AB = workbook_134.add_worksheet()
sheet2_134_BA = workbook_134.add_worksheet()

sheet1_134_AB.write(0, 0, 'Hours per epoch')
sheet2_134_BA.write(0, 0, 'Hours per epoch')
row, col = 1, 0
for i in range(len(data_AB)):
	sheet1_134_AB.write(row, col, (i + 1)*(72/24))
	row += 1

row, col = 1, 0
for i in range(len(data_BA)):
	sheet2_134_BA.write(row, col, (i + 1)*(72/24))
	row += 1

row, col = 1, 1
for read_number in data_AB:
	for score in read_number:
		sheet1_134_AB.write(row, col, score)
		col += 1
	row += 1
	col = 1

row, col = 1, 1
for read_number in data_BA:
	for score in read_number:
		sheet2_134_BA.write(row, col, score)
		col += 1
	row += 1
	col = 1

workbook_134.close()
