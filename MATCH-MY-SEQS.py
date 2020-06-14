import sequence
import os
import argparse
import csv

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", help="FASTA file to query from", required=True)
	parser.add_argument("-q", "--query", help="Query FASTA file", required=True)
	parser.add_argument("-db", "--database", help="Database output file name", required=True)
	parser.add_argument("-o", "--output", help="Output path", default="matchmyseqs")

	args = parser.parse_args()

	seqDict = {}
	tier1seq = ''
	representative = ''
	fasta = {}
	seqsforCSV = {}
	progress = 0

	os.system('makeblastdb -dbtype prot -in ' + args.input + ' -out ' + args.database)

	os.system("blastp -db " + args.database + " -outfmt 3 -num_descriptions 1 -num_alignments 0 -query " + args.query + " -out UniProt_query.txt")
	os.system("grep -e \"^[st][pr]|\" UniProt_query.txt | cut -d\' \' -f1 > UniProt_query.tab")

	with open("UniProt_query.txt", newline='') as f:
		reader = csv.reader(f)
		for row in reader:
			if len(row) > 0 and row[0].startswith('Query'):
				tier1seq = (str(row).split("=")[1][:-2].strip())
			elif len(row) > 0 and (row[0].startswith('tr|') or row[0].startswith('sp|')):
				representative = (str(row).split(" ")[0][2:].strip())
				seqDict[tier1seq] = representative

	done25 = False
	done50 = False
	done75 = False
	for s,rep in seqDict.items():
		total = (len(seqDict))
		seq = (sequence.getSequence(rep,'uniprotkb'))
		seqsforCSV[rep] = str(seq).split(":")[1].strip()
		elem = rep + str(seq)
		progress+=1
		if (progress/total)*100 > 25 and not done25:
			print("25% done")
			done25 = True
		elif (progress/total)*100 > 50 and not done50:
			print("50% done")
			done50 = True
		elif (progress/total)*100 > 75 and not done75:
			print("75% done")
			done75 = True

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
