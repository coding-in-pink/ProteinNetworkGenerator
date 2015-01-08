#!/usr/bin/python
import sys
from chemoUtils import Drugs
from chemoUtils import Targets

#MAIN_____________________________________________________________________________
#USE: main()
#FUNCTION: Reads in the (specified) drugs and targets files, and from the data given in this, 
#1)computes the tanimoto score of each drug pair and determines whether they share a target and 
#2)outputs all this information to the specified out file
#INPUT/OUTPUT: none
#NOTE: this function will quit the program if not enough input files are given
def main():
	if len(sys.argv) < 4:
		print("Missing input files. Program exited.")
		exit()
		drugs = Drugs(sys.argv[1])
		targets = Targets(sys.argv[2])
		out = open(sys.argv[3], 'w+')

		for i in range(drugs.size):
			for j in range(drugs.size):
				if j > i:
					score = drugs.calculate_tanimoto(i,j)
					shared = targets.has_shared_target(drugs.get_code(i),drugs.get_code(j))
					out.write(drugs.get_code(i) + ',' + drugs.get_code(j) + ',' + str(round(score,6)) + ',' + str(shared) + '\n')

#CODE___________________________________________________________________________
if __name__ == '__main__':
	main()