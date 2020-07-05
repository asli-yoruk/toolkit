import sequence
import os
import argparse
import csv
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", help="Input FASTA file", required=True)
	parser.add_argument("-db", "--database", help="Database output file name", required=True)
	parser.add_argument("-r", "--redundancy", nargs='*', help="List of redundancy levels", default=[99,98,95,90,80,70])
	#if user provided the CSV file after MATCHMYSEQS -af, no blastp needed
	parser.add_argument("-af1", "--annotateFile1", help="Tier1 representative sequences")
	#if user provided CSV file with their own tier1 seqs, we blastp
	parser.add_argument("-t1", "--tier1", help="User's Tier1 sequences")
	parser.add_argument("-af2", "--annotateFile2", help="Tier2 representative sequences")
	parser.add_argument("-t2", "--tier2", help="User's Tier2 sequences")
	parser.add_argument("-ml", "--maxlength", help="Max length that the sequence can be", default=800)
	parser.add_argument("-e", "--eval", nargs='*', help="List of evalues", default=[1e-100, 1e-75, 1e-50, 1e-20, 1e-10, 1e-5])
	args = parser.parse_args()

	header = True
	tier1 = {}
	seqNameCol = 0
	#As long as the csv has Name and Sequence, it is fine
	if args.annotateFile1:
		print("these are representative")
		with open(args.annotateFile1, newline='') as f:
				reader = csv.reader(f)
				for row in reader:
					if header:
						for i in range (len(row)):
							if row[i] == 'Representative':
								seqNameCol = i
							elif row[i] == 'Sequence':
								seqCol = i
						header = False
						continue
					tier1[row[seqNameCol]] = row[seqCol]
		print(len(tier1))

	elif args.tier1:

		print("these are special seqs")
		if '.fa' in args.tier1:		
			print("tier 1 sequences are FASTA file")
			tier1 = {}
			tier1_annots = {} # annotations that we want to include in the final dataset
			seq_list = sequence.readFastaFile(args.tier1, sequence.Protein_Alphabet, ignore=True, parse_defline=False)
			print(seq_list)
			sequence.writeFastaFile('tier1_query.fa', seq_list)

			#print(t)
			# tier1 = {} # map from "long" name to actual entry
			# for s in t:
			# # #     #print(s.name)
			# #     tier1[s.name] = s
			# # # #     print(s)
			# 	db100_map_short[sequence.parseDefline(s.name)[0]] = s
			# # for i in args.tier1:
			# # 	print(i)
			# print(tier1)
			
		elif '.csv' in args.tier1:
			tier1 = {}
			tier1_annots = {} # annotations that we want to include in the final dataset
			print("tier 1 sequences are CSV file")
			with open(args.tier1, newline='') as f:
					reader = csv.reader(f)
					headerInput = next(reader)  # gets the first line
					print(headerInput)
					for row in reader:
						#print(row)
						tier1[row[0]] = row[1]

			seq_list = [sequence.Sequence(sequence=seq, name=seqname) for seqname, seq in tier1.items()]
			sequence.writeFastaFile('tier1_query.fa', seq_list)

	tier2 = {}
	tier2_annots = {} # annotations that we want to include in the final dataset

	if args.tier2:
		print("tier2 special sequences have been provided")
		if '.fa' in args.tier2:
			print("tier 2 sequences are FASTA file")
			# for name in copunames:
		 #    	if name in db100_map_short:
			#         s = db100_map_short[name]
			#         tier2[s.name] = s
			#         tier2_annots[s.name] = 'T2=yes'
	  #       print(tier2)
		elif '.csv' in args.tier2:
			print("tier 2 sequences are CSV file")

	elif args.annotateFile2:
		print("tier2 representative sequences have been provided")
		if '.fa' in args.tier2:
			print("tier 2 sequences are FASTA file")
		elif '.csv' in args.tier2:
			print("tier 2 sequences are CSV file")


	db100 = sequence.readFastaFile(args.input, sequence.Protein_Alphabet, ignore=True, parse_defline=False)
	db100_map = {} # map from "long" name to actual entry
	db100_map_short = {} # map from "short" name to entry
	for s in db100:
	    #print(s.name)
	    db100_map[s.name] = s
	#     print(s)
	    db100_map_short[sequence.parseDefline(s.name)[0]] = s
	    #print(sequence.parseDefline(s.name)[0])
	print(len(db100_map))
	print(len(db100_map_short))

	os.system('makeblastdb -dbtype prot -in ' + args.input + ' -out ' + args.database)


	redundancy = []
	for i in range (len(args.redundancy)):
		redundancy.append(int(args.redundancy[i]))


	if args.tier1:
		os.system("blastp -db " + args.database + " -outfmt 3 -num_descriptions 1 -num_alignments 0 -query tier1_query.fa -out UniProt_rep.txt")
		command = "grep -e \"^[st][pr]|\" UniProt_rep.txt | cut -d\' \' -f1 > UniProt_copu.tab"
		print (command)
		os.system(command)

		# Extract the resulting sequence identifiers
		copunames = set([])
		c = 0
		f = open('UniProt_copu.tab', 'rt')
		for row in f:
		    copunames.add(sequence.parseDefline(row.strip())[0])
		    c +=1
		f.close()
		print(c)
		print(len(copunames))

		notfound = []
		for name in copunames:
		    if name in db100_map_short:
		        s = db100_map_short[name]
		        tier1[s.name] = s
		        tier1_annots[s.name] = 'CP=yes'
		    else:
		        notfound.append(name)
		print('Matched', len(copunames)-len(notfound), 'of', len(copunames))
		print(len(tier1))
		print(len(tier1_annots))
		print(len(notfound))
		print(tier1)

		tier1_seqs = [tier1[name] for name in tier1]
		sequence.writeFastaFile('tier1.fa', tier1_seqs)

	print('Ultimately', len(tier1), 'Tier-1 entries')


	if args.annotateFile1:
		seq_list = [sequence.Sequence(sequence=seq, name=seqname) for seqname, seq in tier1.items()]
		print(seq_list)
		sequence.writeFastaFile('tier1.fa', seq_list)

	for rr in redundancy:
	    rs = str(rr)
	    os.system('cd-hit -i ' + args.input + ' -c 0.'+rs+' -T 5 -o db'+rs+'_tps -d 0')

	print(len(db100_map))
	print(len(tier1))
	selected = {}
	for rr in redundancy:
	    selected[rr] = []
	    filename = 'db' + str(rr) + '_tps.clstr'
	    clusters = readCDHIT(filename)
	    for c in clusters: 
	        picked_one = False
	        shortest = None
	        reviewed = None
	        for name in clusters[c]:
	            if name in db100_map:
	                seq = db100_map[name]
	                if shortest:
	                    if len(seq) < len(shortest) and not disqualified(args,seq):
	                        shortest = seq
	                elif not disqualified(args,seq):
	                    shortest = seq
	                if seq.name.startswith('sp|') and not disqualified(args,seq):
	                    reviewed = seq
	                if name in tier1:
	                    selected[rr].append(seq)
	                    picked_one = True
	                    #print("found one!")
	                elif name in tier2:
	                	selected[rr].append(seq)
	                	picked_one = True
	            else:
	            	pass
	                #print('Did not find', name)
	        # If no Tier-1, prefer "reviewed", then shortest length
	        if not picked_one and reviewed:
	            selected[rr].append(reviewed)
	        elif not picked_one and shortest:
	            selected[rr].append(shortest)
	print("FINISHED with redundancy levels")
	print(len(selected[99]))
	print(len(selected[98]))
	print(len(selected[95]))
	print(len(selected[90]))
	print(len(selected[80]))
	print(len(selected[70]))

	for rr in redundancy:
	    print(len(selected[rr]))
	    filename = 'db' + str(rr) + '_tps.fa'
	    sequence.writeFastaFile(filename, selected[rr])

	for rr in redundancy:
		os.system('makeblastdb -dbtype prot -in db' + str(rr) + '_tps.fa -out db-' + str(rr))

	for rr in redundancy:
	    for evalue in [1e-100, 1e-75, 1e-50, 1e-20, 1e-10, 1e-5]:
	        result_file = "dataset-" + str(rr) + '-'+ str(evalue)
	        cmd1 = "blastp -db db-" + str(rr) + " -outfmt 3 -num_descriptions 20000 -num_alignments 0 -num_threads 5 -query tier1.fa -out " + result_file + ".txt -evalue " + str(evalue)
	        print(cmd1)
	        os.system(cmd1)
	        cmd2 = "grep -e \"^[st][pr]|\" " + result_file + ".txt | cut -d\' \' -f1 > " + result_file + ".tab"
	        print(cmd2)
	        os.system(cmd2)
	print('Done')



	totalSeqCount = []
	c = 0
	print(args.eval)
	for evalue in args.eval:
	    for rr in redundancy: 
	        output = []
	        ev = str(evalue)
	        ev = ev[1:]
	        red = str(rr)
	        result_file = "dataset-" + str(rr) + '-'+ str(evalue)
	        # Extract the resulting sequence identifiers
	        names = set([]) # there are duplicate names 
	        f = open(result_file + '.tab', 'rt')
	        for row in f:
	            names.add(row.strip())

	        f.close()
	        tier1_cnt = 0
	        seqs = []

	        #add the user sequences to all of our datasets
	        for seq in seq_list:
	        	c+=1
	        	seqs.append(seq)

	        for name in names:
	            try:
	                seq = db100_map[name]
	                info = ''
	                if name in tier1:
	                    tier1_cnt += 1
	                    info = seq.info + ' ' + tier1_annots[name]
	                elif name in tier2:
	                	tier2_cnt += 1
	                	info = seq.info + ' ' + tier2_annots[name]
	                seqs.append(sequence.Sequence(seq.sequence, seq.alphabet, seq.name, info))
	            except:
	            	pass
	                #print('Did not find', name)
	        sequence.writeFastaFile(result_file + ".fa", seqs)
	        print('Processed', len(seqs), 'for', result_file, ' Tier-1:', tier1_cnt)
	        output = [ev,red,len(seqs)]
	        totalSeqCount.append(output)


	df = pd.DataFrame(totalSeqCount, columns=["E-value", "RedundancyLevel", "Quantity"])

	g = sns.catplot(x="RedundancyLevel", y="Quantity", hue="E-value", data=df,
	                height=20,kind="bar", palette="muted")

	g.despine(left=True)
	g.set_ylabels("Sequence Count")
	g.savefig("sequenceCount.png")

def readCDHIT(filename):
	m = {}
	with open(filename) as fd:
	    rd = csv.reader(fd, delimiter="\t")
	    current = []
	    clustername = None
	    for row in rd:
	        if len(row) > 0: # there is content
	            if row[0].startswith('>'): # new cluster
	                if len(current) > 0 and clustername: # old cluster needs taking care of
	                    m[clustername] = current
	                    current = []
	                clustername = row[0][1:].strip()
	            elif len(row) > 1:
	                field = row[1].split('>')
	                if len(field) > 1:
	                    name = field[1].split('.')[0]
	                    current.append(name)
	print(len(m))
	return m

def disqualified(args,seq):
    # Potential parameters to restrict inclusion
    # Must be set manually, by someone informed
    maxlength = args.maxlength
    if len(seq) > maxlength:
        return True
    return False



if __name__ == "__main__":
	main()