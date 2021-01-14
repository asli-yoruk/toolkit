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
	parser.add_argument("-r", "--redundancy", nargs='*', help="List of redundancy levels", default=[90,80,70])
	parser.add_argument("-t1", "--tier1", help="User's Tier1 sequences")
	parser.add_argument("-t2", "--tier2", help="User's Tier2 sequences")
	parser.add_argument("-ml", "--maxlength", help="Max length that the sequence can be", default=800)
	parser.add_argument("-e", "--eval", nargs='*', help="List of evalues", default=[1e-100, 1e-75, 1e-50, 1e-20, 1e-10, 1e-5])
	args = parser.parse_args()

	tier2 = {}
	tier2_short = {}
	tier2_annots = {} # annotations that we want to include in the final dataset


	if args.tier2:
		print("tier2 sequences have been provided")

		if '.fa' in args.tier2 or '.fasta' in args.tier2:
			print("tier2 sequences are FASTA file")
			tier2db = sequence.readFastaFile(args.tier2, sequence.Protein_Alphabet, ignore=True, parse_defline=False)
			print(str(len(tier2_list)) + " sequences in tier2")
			tier2_list = {} # map from "long" name to actual entry
			tier2_map_short = {} # map from "short" name to entry
			for s in tier2db:
			    tier2_list[s.name] = s
			    tier2_map_short[sequence.parseDefline(s.name)[0]] = s
		else:
			print("Please provide FASTA file for tier-2")

	if args.tier1:
		tier1 = {}
		tier1_annots = {} # annotations that we want to include in the final dataset
		print("Tier-1 sequences have been provided")
		if '.fa' in args.tier1 or '.fasta' in args.tier1:

			print("Tier-1 sequences are provided as a FASTA file")
			tier1db = sequence.readFastaFile(args.tier1, sequence.Protein_Alphabet, ignore=True, parse_defline=False)
			tier1_list = {}
			for s in tier1db:
				tier1_list[s.name] = "".join(s.sequence)
			print("Tier-1 has " + str(len(tier1_list)) + " sequences") 


		else:
			print("Please provide FASTA file for tier-1")


	db100 = sequence.readFastaFile(args.input, sequence.Protein_Alphabet, ignore=True, parse_defline=False)
	db100_map = {} # map from "long" name to actual entry
	db100_map_short = {} # map from "short" name to entry
	for s in db100:
	    db100_map[s.name] = s
	    db100_map_short[sequence.parseDefline(s.name)[0]] = s
	print("Database has " + str(len(db100_map)) + " sequences")

	for rr in args.redundancy:
	    rs = str(rr)
	    os.system('cd-hit -i ' + args.input + ' -c 0.'+rs+' -T 5 -o db'+rs+' -d 0')


	selected = {}
	for rr in args.redundancy:
	    selected[rr] = []
	    filename = 'db' + str(rr) + '.clstr'
	    clusters = readCDHIT(filename)
	    for c in clusters: 
	        picked_one = False
	        shortest = None
	        reviewed = None
	        for name in clusters[c]:
	            if name in db100_map:
	                seq = db100_map[name]
	                if shortest:
	                    if len(seq) < len(shortest) and not disqualified(seq,args):
	                        shortest = seq
	                elif not disqualified(seq,args):
	                    shortest = seq
	                if seq.name.startswith('sp|') and not disqualified(seq,args):
	                    reviewed = seq
	                if name in tier1_list:
	                    #print("this one orig" + str(seq))
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

	for rr in args.redundancy:
	    filename = 'db' + str(rr) + '.fa'
	    sequence.writeFastaFile(filename, selected[rr])

	for rr in args.redundancy:
		os.system('makeblastdb -dbtype prot -in db' + str(rr) + '.fa -out db-' + str(rr))

	# for rr in args.redundancy:
	#     for evalue in args.evalue:
	#         result_file = "dataset-" + str(rr) + '-'+ str(evalue)
	#         cmd1 = "blastp -db db-" + str(rr) + " -outfmt 3 -num_descriptions 20000 -num_alignments 0 -num_threads 5 -query " + args.tier1 + " -out " + result_file + ".txt -evalue " + str(evalue)
	#         print(cmd1)
	#         os.system(cmd1)

	grab = False

	for rr in args.redundancy:
	    for evalue in args.eval:
	        c=0
	        tpsIdentifier = set([])
	        seqs = []
	        result_file = "dataset-" + str(rr) + '-'+ str(evalue)
	        f = open(result_file + '.txt', 'rt')
	        for row in f:
	            if row.startswith('Sequences'):
	                grab = True
	                continue
	            if grab == True:
	                if row.startswith('Lambda'):
	                    grab = False
	                if not row.strip() == "":
	                    identifier = row.split(' ')[0]
	                    if identifier != "Lambda":
	                        tpsIdentifier.add(identifier)
	                    

	        
	        for name in tpsIdentifier:
	            try:
	                seq = db100_map[name]
	                info = ''
	                seqs.append(sequence.Sequence(seq.sequence, seq.alphabet, seq.name, info))
	            except:
	                pass
	        sequence.writeFastaFile(result_file + ".fa", seqs)
	        print(result_file + " has " + str(len(seqs)) + "sequences")

	print('Done')


	totalSeqCount = []
	c = 0
	for evalue in args.eval:
	    for rr in args.redundancy: 
	        output = []
	        ev = str(evalue)
	        ev = ev[1:]
	        red = str(rr)
	        result_file = "dataset-" + str(rr) + '-'+ str(evalue)
	        a = sequence.readFastaFile(result_file + '.fa', sequence.Protein_Alphabet, ignore=True, parse_defline=False)

	        names = set([]) 
	        for s in a:
	            names.add(s.name)
	        tier1_cnt = 0
	        tier2_cnt = 0
	        seqs = []
	        for name in names:
	            try:
	                seq = db100_map[name]
	                info = ''
	                if name in tier1_list:
	                    tier1_cnt += 1
	                    #info = seq.info + ' ' + tier1_annots[name]
	                elif name in tier2:
	                    tier2_cnt += 1
	                    #info = seq.info + ' ' + tier2_annots[name]
	                seqs.append(sequence.Sequence(seq.sequence, seq.alphabet, seq.name, info))
	            except:
	            	pass
	                #print('Did not find', name)
	        print('Processed', len(seqs), 'for', result_file, ' Tier-1:', tier1_cnt, ' Tier-2:', tier2_cnt)
	        output = [ev,red,len(seqs)]
	        totalSeqCount.append(output)


	plotSeqs(totalSeqCount)



def plotSeqs(totalSeqCount):
	plt.figure()
	df = pd.DataFrame(totalSeqCount, columns=["E-value", "RedundancyLevel", "Quantity"])

	g = sns.catplot(x="RedundancyLevel", y="Quantity", hue="E-value", data=df,
	                height=20,kind="bar", palette="muted")
	#plt.legend(loc="lower right", prop={'size': 34})
	sns.set(font_scale=30)

	g.despine(left=True)
	g.set_ylabels("Sequence Count")
	g.savefig("sequenceCount.png")




def readCDHIT(filename):
    ''' Function to place all clusters from CD-HIT in a map.
        Each cluster is a list of names. '''
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
        m[clustername] = current
    return m



def disqualified(seq,args):
    # Potential parameters to restrict inclusion
    # Must be set manually, by someone informed
    maxlength = args.maxlength
    if len(seq) > maxlength:
        return True
    return False




if __name__ == "__main__":
	main()