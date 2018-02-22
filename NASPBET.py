import re
import sys

args = sys.argv

## Relative path
defaultdb = "resources/default.ndb"

if "-h" in args or len(args) < 2:
	print "NASPBET_v1.0"
	print "\nThis script will isolate binding sites related to a \ngiven protein in all sequences of a given multi-fasta file."
	print "\nHelp Parameters:\nThese parameters can be used to print help menus.\nOnly one of these parameters should be used at a time."
	print "		-h"
	print "			prints this help menu listing script parameters."
	print "		-hd"
	print "			prints an explanation and example of the way in\n			which user databases need to be set up for the script."
	print "		-hp"
	print "			prints the list of binding proteins that are\n			currently available in the default database."
	print "		-hdp database.eg"
	print "			prints the list of binding proteins that are\n			currently available in the given database."
	print "\n\nMandatory Parameters:\nThese parameters are mandatory for the script to run."
	print "		-i \"fasta-file.example\""
	print "			points to the relative path of the target multi-fasta\n			file in which the script should look for binding sequences."
	print "			Input sequences must be at least 5 letters long to be considered."
	print "		-o \"output.tsv\""
	print "			points to the relative path of the file to create to write the results TSV file.\n"
	print "			The file will be overwritten if it already exists."
	print "\n\nOptional Parameters:\nThese parameters are optional, with respective default values if not specified."
	print "		-d \"database.eg\""
	print "			points to the relative path of the database to be used for the search."
	print "			If not specified the default database will be used.\n			For more information on database format use the -hd parameter."
	print "		-p \"prot1,prot2,prot3\""
	print "			represents the proteins for which binding sites will\n			be searched for in the target file. These must\n			be seperated by commas, with no spaces between them.\n			If no protein is indicated, all proteins will be considered"

elif "-hd" in args:
	print "The database mimics a fasta-file, announcing each new binding\nprotein with a >Prot-Name line."
	print "Following the name of the protein, all binding sequences or\npatterns follow successively, listing all matches for that protein."
	print "\nFor a list of pattern translation for nucleic acid sequences, please visit:\n"
	print "https://www.bioinformatics.org/sms/iupac.html"
	print "\n\nAn example database file is shown below:"
	print "\n>HSF2\nAGAANNTTCG\nAGAANNTTCT\n\n>GST-B6ZF\nTNCTTTCNAGGAAT\nTNCTTTCNAGGGAT"

elif "-hp" in args:
	print "The following is a list of proteins available in the default database."
	print "This database combines information from HTPSelex and SELEX_DB."
	print "\nDatabase\nProtein\nInfo Link"
	for line in open(defaultdb,'r'):		
		match = re.search(r'^>([^\n]+)', line)
		if match:
			gene = match.group(1).split(" ")
			db = gene[0].replace(":", "")
			name = "\t" + gene[1]
			link = "\thttp://www.uniprot.org/uniprot/?query=" + name.strip() + "&sort=score"
			print db + name + link

elif "-hdp" in args:
	print "The following is a list of proteins available in the given database."
	print "\nDatabase\nProtein\nInfo Link"
	d = args.index("-hdp") + 1
	for line in open(args[d],'r'):		
		match = re.search(r'^>([^\n]+)', line)
		if match:
			gene = match.group(1).split(" ")
			db = gene[0].replace(":", "")
			name = "\t" + gene[1]
			link = "\thttp://www.uniprot.org/uniprot/?query=" + name.strip() + "&sort=score"
			print db + name + link

else:

	## Global variables
	inputNames = []
	inputSequences = []

	dbGenes = []
	dbLinks = []
	dbSequences = []
	curList = []

	givenGenes = []
	currentGene = ""

	currentSeq = ""

	## Starts with 1 to account for name which is automatically in args
	## This will be used to make sure all arguments were used at the end of processing and output a warning if needed
	count = 1


	## An input file is mandatory
	if "-i" not in args:
		print "ERROR: Input file not given. Use -h parameter for full usage information."
	else:
		count += 2
		i = args.index("-i") + 1
		## Read the input file and extract sequences
		openInput = open(args[i], 'r')

		for eachLine in openInput:
			line = eachLine.strip()
			match = re.search(r'^>([^\n]+)', line)
			if match:
				## Input old sequence when new label is found
				if(currentSeq != ""):
					inputSequences.append(currentSeq)
					currentSeq = ""
				inputNames.append(match.group(1))
			## This is to also accept FASTA files that have line
			## breaks in long sequences
			elif len(line) > 0:
				currentSeq += line

		## Also save the last one
		if(currentSeq != ""):
			inputSequences.append(currentSeq)
			currentSeq = ""


		## Output file is also necessary
		if "-o" not in args:
			print "ERROR: Output file not given. Use -h parameter for full usage information."
		else:
			count += 2

			## Will save the list of proteins when they are given
			if "-p" in args:
				count += 2
				p = args.index("-p") + 1
				givenProts = [pe.strip() for pe in args[p].split(',')]
				if(len(givenProts) <= 10):
					print "Looking for the binding sequences associated with these genes:"
					for eachProt in givenProts:
						print "\t\t" + eachProt
				else:
					print "Looking for the binding sequences associated with " + len(givenProts) + " given genes."

			## Announce either names of proteins (if <= 10), number (if > 10),
			## or that it is looking through all proteins (if not given)
			else:
				print "Looking for all binding proteins in the database."


			## Change database
			if "-d" in args:
				count += 2
				d = args.index("-d") + 1
				## Read the database file
				openDatabase = open(args[d], 'r')


			else:
				openDatabase = open(defaultdb,'r')	
			


			## Either look for all or just the proteins given
			## These loops put the database information in lists
			## to iterate through them afterwards
			if len(givenGenes) == 0:
				for line in openDatabase:		
					match = re.search(r'^>([^\n]+)', line)
					if match:
						if len(curList) > 0:
							dbSequences.append(curList)
						curList = []
						currentGene = match.group(1)
						dbGenes.append(currentGene.replace("\r",""));
						link = "http://www.uniprot.org/uniprot/?query=" + currentGene.split(" ")[1].strip() + "&sort=score"
						dbLinks.append(link)

					elif len(line.strip()) > 0 and currentGene != "":
						curSeq = line.strip()
						curList.append(curSeq)
				if len(curList) > 0:
					dbSequences.append(curList)

			else:
				skip = True
				for line in openDatabase:
					match = re.search(r'^>([^\n]+)', line)
					if match:
						if len(curList) > 0:
							dbSequences.append(curList)
						curList = []
						skip = True
						currentGene = match.group(1)
						if currentGene.split(" ")[0] in givenGenes:
							dbGenes.append(currentGene.replace("\r",""));
							skip = False
							link = "http://www.uniprot.org/uniprot/?query=" + currentGene.split(" ")[1].strip() + "&sort=score"
							dbLinks.append(link)

					elif not skip and len(line.strip()) > 0 and currentGene != "":
						curSeq = line.strip()
						curList.append(curSeq)
				if len(curList) > 0:
					dbSequences.append(curList)


			## Should never happen with default, but can happen with given database
			if(len(dbSequences) != len(dbGenes)):
				print "ERROR READING DATABASE"
				print "Number of given names: " + str(len(dbGenes))
				print "Number of given sequence lists: " + str(len(dbSequences))
				print "Make sure there are no lines labeled '>' with no sequences following them"



			## Global variables to be accessible anywhere
			outputGenes = []
			outputGeneCounts = []
			outputGeneAvg = []
			# outputSequences = []
			# outputSeqCounts = []


			## Will create the regular expressions
			fromChange = ["N","R","Y","S","W","K","M","B","D","H","V","U"]
			toChange = ["(A|C|G|T)","(A|G)","(C|T)","(C|G)","(A|T)","(G|T)","(A|C)","(C|G|T)","(A|G|T)","(A|C|T)","(A|C|G)","T"]


			for w in range(len(inputSequences)):
				eachSeq = inputSequences[w]
				name = inputNames[w]
				tempOutGenes = []
				tempOutGeneCounts = []
				tempOutGeneAvg = []
				# tempOutGeneSeqs = []
				# tempOutGeneSeqCounts = []

				bindingSize = 0

				for i in range(len(dbGenes)):
					curGeneCount = 0
					currentGene = dbGenes[i]

					pos = []
					# tempOutSequences = []
					# tempOutSeqCount = []

					curList = dbSequences[i]
					for curSeq in curList:
						curSeqCount = 0

						regEx = curSeq

						for r in range(len(fromChange)):
							regEx.replace(fromChange[r],toChange[r])

						## This allows for the regular exp to be made dynamically
						p = re.compile(regEx)
						for match in p.finditer(eachSeq):
							start = match.start()
							old = False
							for posit in pos:
								if start => posit[0] and start <= posit[1]:
									old = True
									pos.append([start,len(match.group())])

							if not old:
								# curSeqCount += 1
								curGeneCount += 1
								bindingSize += len(curSeq)
								pos.append([start,len(match.group())])


						# if(curSeqCount > 0):

						# 	tempOutSequences.append(curSeq)
						# 	tempOutSeqCount.append(curSeqCount)

					if(curGeneCount > 0):
						curGeneAvg = bindingSize*1.0/curGeneCount
						tempOutGenes.append(currentGene)
						tempOutGeneCounts.append(curGeneCount)
						tempOutGeneAvg.append(curGeneAvg)
						# tempOutGeneSeqs.append(tempOutSequences)
						# tempOutGeneSeqCounts.append(tempOutSeqCount)

				outputGenes.append(tempOutGenes)
				outputGeneCounts.append(tempOutGeneCounts)
				outputGeneAvg.append(tempOutGeneAvg)
				# outputSequences.append(tempOutGeneSeqs)
				# outputSeqCounts.append(tempOutGeneSeqCounts)

				

				i = args.index("-o") + 1
				openOutput = open(args[i], 'w')
				openOutput.write("Input Sequence\tProtein\tHits\tAvg Seq Size\tInfo URL\n")
				for j in range(len(tempOutGenes)):
					gene = tempOutGenes[j]
					openOutput.write("\n"+ name + "\t" + str(gene) + "\t" + str(tempOutGeneCounts[j]) + "\t" + str(tempOutGeneAvg[j]) + "\t" + dbLinks[dbGenes.index(gene)])



		#	for i in range(len(outputGenes)):
		#		name = inputNames[i]
		#		geneList = outputGenes[i]
		#		geneCount = outputGeneCounts[i]
		#		for j in range(len(geneList)):
		#			print "{:10}{:10}{:10}".format(name,geneList[j],geneCount[j])




			if(count < len(args)):
				print "WARNING: Some arguments might have been ommitted. Make sure you enter your arguments are properly written."
				print "This is often caused by spaces in the paths or gene list. Make sure to use quotation marks when paths contain spaces."
				print "For full usage information, use the -h parameter."
