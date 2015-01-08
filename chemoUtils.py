#Drug object
class Drugs():
	#INITIALIZE: DRUGS()__________________________________________________________
		#USE: object = Drugs(filename) 
		#file containing all the drug information
		#FUNCTION: reads in the input files and inputs data in to an internal data structure
		#INPUT: filename of the data file, which should be in the specified format (str)
		def __init__(self, filename):
			self.all = []
			self.size = 0
			for line in file(filename).readlines()[1:]:
				line = line.rstrip('\r\n')
				line = line.split(',')
				self.all.append(line[:2] + [float(x) for x in line[2].split()])
				self.size += 1
				self.all_t = zip(*self.all)

	#GET_CODE_____________________________________________________________________
		#USE: drug_code = self.get_code(i)
		#FUNCTION: returns the accession number of the drug, given its index in the list
		#INPUT: index in file (int)
		#OUTPUT: drug code, i.e. accession number (str)
		def get_code(self,index):
			return self.all[index][0]

    #GET_INDEX_____________________________________________________________________
		#USE: drug_index = self.get_index("DB00714")
		#FUNCTION: returns the index of the drug in the file, given its accession number
		#INPUT: accession number (str)
		#OUTPUT: index (int)
		def get_index(self,name):
			return self.all_t[0].index(name)

    #GET_ATTRIBUTES_____________________________________________________________________
		#USE: attributes = self.get_attributes(i)
		#FUNCTION: returns the list of the drug attributes, given the index of the drug
		#INPUT: index (int)
		#OUTPUT: list of attributes of the drug at that index (list[<int>])
		def get_attributes(self,index):
			return self.all[index][2:]

    #CALCULATE_TANIMOTO____________________________________________________
		#USE: tanimoto_score = self.calculate_tanimoto(index1,index2)
		#FUNCTION: calculates the fraction of the intersecting attributes over the union
		#of the union of the attributes of two drugs
		#INPUT: the indices of the 2 drugs you wish to compare (int,int)
		#OUTPUT: the tanimoto score of the attributes of those 2 drugs (float)
		def calculate_tanimoto(self,index1,index2):
			arr1 = self.get_attributes(index1)
			arr2 = self.get_attributes(index2)
			intersect = len(set.intersection(set(arr1),set(arr2)))
			return float(intersect)/(len(arr1) + len(arr2) - intersect)

    #CALCULATE_T_SUMMARY____________________________________________________
    	#USE: t_summary = self.calculate_t_summary(list1,list2)
    	#FUNCTION: calculates the t_summary (sum of all pair wise taninomoto scores between
    	#two lists of ligands that are above 0.5)
		#INPUT: two lists of ligand accession numbers (list[<str>],list[<str>])
		#OUTPUT: the resutling t_summary of those two lists of ligands (float)
		def calculate_t_summary(self,list1,list2):
			t_summary = 0.0
			for l1 in list1:
				for l2 in list2:
					score = self.calculate_tanimoto(l1,l2)
					if score > 0.5:
						t_summary += score
			return t_summary

#Targets object
class Targets():
	#Initialize: Targets__________________________________________________________
		#USE: object = Targets(filename) 
		#FUNCTION: reads in the input files and inputs data in to an internal data structure
		#INPUT: name of file containing all the target data, which should be in the specified format (str)
		def __init__(self, filename):
			self.all = []
			for line in file(filename).readlines()[1:]:
				line = line.rstrip('\r\n')
				self.all.append(line.split(','))
				self.all_t = zip(*self.all)

	#GET_TARGETS_____________________________________________________________________
		#USE: targets = self.get_targets("DB00714")
		#FUNCTION: returns the list of the drug targets, given the accession number of the drug
		#INPUT: accession number of drug (str)
		#OUTPUT: list of uniprot accession numbers of proteins that drug targets (list[<str>])
		def get_targets(self,drug_code):
			return [self.all_t[1][i] for i, x in enumerate(self.all_t[0]) if x == drug_code]

	#GET_LIGANDS_____________________________________________________________________
		#USE: ligands = self.get_targets("P30968")
		#FUNCTION: returns the list of the drugs that bind to a specified protein
		#INPUT: uniprot accession number of protein (str)
		#OUTPUT: list of accession numbers of drugs that target that protein (list[<str>])
		def get_ligands(self,protein_code):
			return [self.all_t[0][i] for i, x in enumerate(self.all_t[1]) if x == protein_code]

	#SHARED_TARGET________________________________________________________________
		#USE: indicator = self.shared_target("DB00714", "DB00718")
		#FUNCTION: determines whether the two drugs share a protein target by comparing the set
		#of protein targets of drugs 1 and 2 for overlaps
		#INPUT: the codes of the two drugs of interest (str,str)
		#OUTPUT: a binary indicator (1=yes,0=no) of whether the two drugs share a target (int)
		def has_shared_target(self,drug_code1,drug_code2):
			if len(set.intersection(set(self.get_targets(drug_code1)),set(self.get_targets(drug_code2)))) > 0:
				return 1
			else:
				return 0

#Proteins object
class Proteins():
	#INITIALIZE: DRUGS()__________________________________________________________
		#USE: object = Drugs(filename) 
		#file containing all the drug information
		#FUNCTION: reads in the input files and inputs data in to an internal data structure
		#INPUT: filename of the data file, which should be in the specified format (str)
		def __init__(self, filename):
			self.all = []
			self.size = 0
			for line in file(filename).readlines()[1:]:
				line = line.rstrip('\r\n')
				self.all.append(line.split(','))
				self.size += 1
			self.all_t = zip(*self.all)

    #GET_code_____________________________________________________________________
		#USE: protein_code = self.get_code(i)
		#FUNCTION: returns the uniprot accession number, given the index of the protein in the file
		#INPUT: index (int)
		#OUTPUT: uniprot accession number (str)
		def get_code(self,index):
			return self.all[index][0]

	#GET_NAME_____________________________________________________________________
		#USE: name = self.get_name("P30968")
		#FUNCTION: returns the name of the specified protein
		#INPUT: protein index in file (int)
		#OUTPUT: protein name (str)
		def get_name(self,index):
			return self.all[index][1]

	#GET_INDICATIONS_____________________________________________________________________
		#USE: indications = self.get_indications()
		#FUNCTION: returns a list of all the indications of a specified protein
		#INPUT: protein index in file (int)
		#OUTPUT: list of associated indications (list[<str>])
		def get_indications(self,index):
			return self.all[index][2]
