import sequence
import os
import argparse
import csv

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", help="FASTA file to query from", required=True)
	parser.add_argument("-q", "--query", help="Query FASTA file", required=True)
	parser.add_argument("-db", "--database", help="Database output file name", required=True)
	parser.add_argument("-r", "--reference", help="Reference database ", default="uniprotkb")
	parser.add_argument("-o", "--output", help="Output path", default="matchmyseqs")

	args = parser.parse_args()

	seqDict = {}
	tier1seq = ''
	representative = ''
	fasta = {}
	seqsforCSV = {}
	progress = 0
	tier1 = {}
	tier1_annots = {} # annotations that we want to include in the final dataset

	os.system('makeblastdb -dbtype prot -in ' + args.input + ' -out ' + args.database)

	db = sequence.readFastaFile(args.input, sequence.Protein_Alphabet, ignore=True, parse_defline=False)
	db_map = {} # map from "long" name to actual entry
	db_map_short = {} # map from "short" name to entry
	for s in db:
	    db_map[s.name] = s
	    db_map_short[sequence.parseDefline(s.name)[0]] = s
	print("Database size is " + str(len(db_map)))

	print("Blast started, this might take a bit depending on your dataset size")
	os.system("blastp -db " + args.database + " -outfmt 3 -num_descriptions 1 -num_alignments 0 -query " + args.query + " -out query.txt")
	

	if args.reference == 'uniprotkb':
		os.system("grep -e \"^[st][pr]|\" query.txt | cut -d\' \' -f1 > UniProt_query.tab")

		# Extract the resulting sequence identifiers
		repSeqNames = set([])
		f = open('UniProt_query.tab', 'rt')
		for row in f:
		    repSeqNames.add(sequence.parseDefline(row.strip())[0])
		f.close()
		print(str(len(repSeqNames)), " representative sequences have been found")

		#Annot the representative sequences 
		notfound = []
		for name in repSeqNames:
			if name in db_map_short:
				s = db_map_short[name]
				seqsforCSV[s.name] = "".join(s)
			else:
				notfound.append(name)
		print('Matched', len(repSeqNames)-len(notfound), 'of', len(repSeqNames))

		with open("query.txt", newline='') as f:
			reader = csv.reader(f)
			for row in reader:
				if len(row) > 0 and row[0].startswith('Query'):
					querySeq = (str(row).split("=")[1][:-2].strip())
				elif len(row) > 0 and (row[0].startswith('tr|') or row[0].startswith('sp|')):
					representative = (str(row).split(" ")[0][2:].strip())
					seqDict[querySeq] = representative

	elif args.reference == 'refseq':
		grab = False
		repSeqNames = set([])

		with open("query.txt", newline='') as f:
			reader = csv.reader(f)
			for row in reader:
				if len(row) > 0 and row[0].startswith('Query'):
					querySeq = (str(row[0]).split("=")[1][:-2].strip().split(" ")[0])
				elif len(row) > 0 and row[0].startswith('Sequences'):
					grab = True
					continue
				elif grab == True:
					if len(row) > 0 and not row[0].strip() == "":
						representative = (row[0].split('.')[0]+"."+row[0].split('.')[1].split(" ")[0])
						repSeqNames.add(representative)
						seqDict[querySeq] = representative
						grab = False
			#print(len(repSeqNames))

			notfound = []

			for name in repSeqNames:
				if name in db_map_short:
					s = db_map_short[name]
					seqsforCSV[s.name] = "".join(s)
				else:
					notfound.append(name)
			print('Matched', len(repSeqNames)-len(notfound), 'of', len(repSeqNames))
			    
			print(len(repSeqNames) , " representative sequences found for " + args.query)






	# done25 = False
	# done50 = False
	# done75 = False
	# for s,rep in seqDict.items():
	# 	total = (len(seqDict))
	# 	seq = (sequence.getSequence(rep,'uniprotkb'))
	# 	seqsforCSV[rep] = str(seq).split(":")[1].strip()
	# 	elem = rep + str(seq)
	# 	progress+=1
	# 	if (progress/total)*100 > 25 and not done25:
	# 		print("25% done")
	# 		done25 = True
	# 	elif (progress/total)*100 > 50 and not done50:
	# 		print("50% done")
	# 		done50 = True
	# 	elif (progress/total)*100 > 75 and not done75:
	# 		print("75% done")
	# 		done75 = True

	faOut = args.output + '.fa'

	seq_list = [sequence.Sequence(sequence=seq, name=seqname) for seqname, seq in seqsforCSV.items()]

	sequence.writeFastaFile(faOut, seq_list)
	
	csvOut = args.output + '.csv'

	with open(csvOut, 'w', newline='') as f:
	    fieldnames = ['Name', 'Representative', 'Sequence']
	    thewriter = csv.DictWriter(f, fieldnames=fieldnames)
	    
	    thewriter.writeheader()
	    for given,rep in seqDict.items():
	    	thewriter.writerow({'Name' : given, 'Representative' : rep, 'Sequence' : seqsforCSV[rep]})


if __name__ == "__main__":
    main()
