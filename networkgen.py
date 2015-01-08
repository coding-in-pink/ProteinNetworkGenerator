#!/usr/bin/python
import sys
import random
import argparse
from chemoUtils import Drugs
from chemoUtils import Targets
from chemoUtils import Proteins
import pvalue

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
	parser.add_argument("proteins_fn", type=str, help="Name of file containing all protein targets")
	return parser

#WRITE_OUT_SIF_____________________________________________________________________________
#USE: write_out_sif(result)
#FUNCTION: Creates and writes to the sif out file
#INPUT: List of lines to write to out file (list[<str>])
def write_out_sif(result):
	out_sif = open('network.sif', 'w+')
	result.sort()
	for line in result:
		out_sif.write(line + '\n')

#WRITE_OUT_REST_____________________________________________________________________________
#USE: write_out_rest(proteins,indices)
#FUNCTION: Creates and writes to the name and indicator out files
#INPUT: Set of indices of proteins to write the information of (set(<int>))
def write_out_rest(proteins,indices):
	out_name = open('name.nodeAttr', 'w+')
	out_name.write("name\n")
	out_indication = open('indication.nodeAttr', 'w+')
	out_indication.write("indication\n")
	for index in indices:
		out_name.write(proteins.get_code(index) + " = " + proteins.get_name(index) + '\n')
		out_indication.write(proteins.get_code(index) + " = " + proteins.get_indications(index) + '\n')

#MAIN_____________________________________________________________________________
#USE: main()
#FUNCTION: Reads in the (specified) drug, protein, and target files, and from the data given in this, 
#1)computes which share edges using the pvalue module
#2)outputs all this information to the three specified out files
#INPUT/OUTPUT: none
def main():
	args = get_args().parse_args()
	if args.r:
		random.seed(args.r)
		drugs = Drugs(args.drugs_fn)
		targets = Targets(args.targets_fn)
		proteins = Proteins(args.proteins_fn)
		result = []
		protein_indices = []
	#iterates through each protein pair and computes the p-value
	for i in range(proteins.size):
		for j in range(proteins.size):
			if j > i:
				p = pvalue.calculate_p_bootstrap(args.n, drugs, targets, proteins.get_code(i), proteins.get_code(j))
				if p <= 0.05:
					result.append(proteins.get_code(i) + " edge " + proteins.get_code(j))
					protein_indices.append(i)
					protein_indices.append(j)
	#write to output files
	write_out_sif(result)
	write_out_rest(proteins,set(protein_indices))

#CODE_________________________________________________________________________________
if __name__ == '__main__':
	main()