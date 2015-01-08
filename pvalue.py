#!/usr/bin/python
import sys
import random
import argparse
from chemoUtils import Drugs
from chemoUtils import Targets

#GET_ARGS_____________________________________________________________________________
#USE: parser = get_args()
#FUNCTION: Parses command line arguments and provides error if incorrect/not enough are given
#Output: parser object containing all supplied arguments
#NOTE: this function will quit the program if not enough input files are given
def get_args():
	parser = argparse.ArgumentParser(description="Generates a bootstrap p-value for the comparison of two proteins.")
	parser.add_argument("-r", type=int, help="pseudo-random number generator seed.")
	parser.add_argument("-n", type=int, default = 100, help="Number of iterations")
	parser.add_argument("drugs_fn", type=str, help="Name of the file containing all the drug information")
	parser.add_argument("targets_fn", type=str, help="Name of the file containing all the target information")
	parser.add_argument("proteinA", type=str, help="UniProt accession number for the first protein")
	parser.add_argument("proteinB", type=str, help="UniProt accession number for the second protein")
	return parser

#RANDOM_LIST_____________________________________________________________________________
#USE: list1 = random_list(1000,50)
#FUNCTION: Creates random list of n ints in the range [0,upper)
#INPUT: size of output list and the upper bound on the element value (int,int)
#OUTPUT: random list of integers (list[<int>])
def random_list(n,upper):
	result = []
	for i in range(n):
		result.append(int(upper*random.random()))
	return result

#CALCULATE_P_BOOTSTRAP____________________________________________________________________
#USE: p = calculate_p_bootstrap(n,drugs,targets,proteinA,proteinB)
#FUNCTION: Computes the p-value of the two proteins
#INPUT: number of iterations (int), Drug object, Target object, uniprot code of first protein(str),
#uniprot code of second protein(str)
#OUTPUT: resulting p-value(float)
def calculate_p_bootstrap(n,drugs,targets,proteinA,proteinB):
	#get T-summary of the given proteins
		ligandsA = [drugs.get_index(x) for x in targets.get_ligands(proteinA)]
		ligandsA_size = len(ligandsA)
		ligandsB = [drugs.get_index(x) for x in targets.get_ligands(proteinB)]
		ligandsB_size = len(ligandsB)
		tSummaryAB = drugs.calculate_t_summary(ligandsA,ligandsB)
	#get all threshold T-summary scores
		sum_tSummary = 0.0
		for i in range(n):
			tSummary = drugs.calculate_t_summary(random_list(ligandsA_size,drugs.size),random_list(ligandsB_size,drugs.size))
			if tSummary >= tSummaryAB:
				sum_tSummary += 1
		return sum_tSummary/n

#MAIN_____________________________________________________________________________
#USE: main()
#FUNCTION: Reads in the (specified) drugs and targets files, and calculates and prints the
#p-score of the two given protein arguments
#INPUT/OUTPUT: none
def main():
	args = get_args().parse_args()
	if args.r:
		random.seed(args.r)
	drugs = Drugs(args.drugs_fn)
	targets = Targets(args.targets_fn)
	print calculate_p_bootstrap(args.n,drugs,targets,args.proteinA,args.proteinB)

#CODE_______________________________________________________________________________
if __name__ == '__main__':
	main()