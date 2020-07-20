import argparse
import sequence
import csv
#TODO havent added graphs yet

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-t", "--type", help="Aligned or not aligned file", default="a")
	parser.add_argument("-a", "--add", help="File to be merged with the matches")
	parser.add_argument("-i", "--input", nargs="+", help="Input FASTA file", required=True)
	parser.add_argument("-m", "--motif",nargs="*", help="Motif regex")

	args = parser.parse_args()
	files = args.input

	motifs = getMotifs(args)

	matches, d = motifSearch(args,motifs)

	writeMatches(args, matches, d)
	d = {}

def motifSearch(args, motifs):
	names = dict()
    # Create a dictionary for each motif, 
    # which associates columns with a match to the number of sequences with that match
	for s in args.input:
		if (args.type == 'a'):
		    aln = sequence.readFastaFile(s, sequence.Protein_Alphabet, gappy=True, ignore=True, parse_defline=False)
		    ali = sequence.Alignment(aln)
		else:
			ali = sequence.readFastaFile(s, sequence.Protein_Alphabet, ignore=True, parse_defline=False)

		dictionary = dict()

		for a in ali:
			T1 = False
			if 'T1' in a.info:
				seqName = a.name+"+T1=yes"
				T1 = True
			if 'T2' in a.info:
				seqName = a.name+"+T2=yes"
			elif T1 == False:
				seqName = a.name
			    #print(seqName)
			#seqName = str(a).split(":")[0]
			seqSequence = str(a).split(":")[1].strip()
			thisset = set()

			for m in motifs:
				number = 0
				for i in range(len(args.motif)):
					if str(m) == args.motif[i]:
						number = i
				if (args.type == 'a'):
					result = m.search(a, gappy=True)
				else:
					result = m.search(a)

				if (len(result) > 1):
					for r in result: #position, matched string, score
						motifStart, foundMotif, n = r
						addThis = ('motif' + str(number+1) + ',' + str(motifStart) + ',' + str(foundMotif) + ',' + str(len(foundMotif)+motifStart))
						thisset.add(addThis)
				elif (len(result) == 1):
					motifStart, foundMotif, n = str(result).split(",")
					addThis = ('motif' + str(number+1) + ',' + str(motifStart[2:]) + ',' + str(foundMotif) + ',' + str(len(foundMotif)+int(motifStart[2:])))
					thisset.add(addThis)
				else:
					pass
			thisset.add('Sequence,' + seqSequence+', '+', ')
			dictionary[seqName] = thisset
		names[s] = dictionary
	#print(names)
	return(names,dictionary)

def writeMatches(args, matches,d):
	count = 0
	for keyBigD,valueBigD in matches.items():
		#will create same number of input file has been provided
		filename = str(args.input[count]).split(".")[0]
		filename += '.csv'

		with open(filename, 'w', newline='') as f:
			fieldnames = ['Name']
			for i in range (len(args.motif)):
				fieldnames.append('motif' + str(i+1))
			fieldnames.append('Sequence')

			thewriter = csv.DictWriter(f, fieldnames=fieldnames)
			thewriter.writeheader()
			for keySmallD,valueSmallD in valueBigD.items():
				dictionary = dict()
				for a in valueSmallD:
					whichMotif, start, found, end = str(a).split(",")
					for x in range (len(fieldnames)):
						if whichMotif == fieldnames[x]:
							if dictionary.get(whichMotif) != None:
								dictionary[whichMotif] += ',' + start
							else:
								dictionary[whichMotif] = start

					dictionary['Name'] = keySmallD
				if bool(dictionary):
					thewriter.writerow(dictionary)
		count += 1
		print("Done")



def getMotifs(args):
    motifObj = []

    for m in args.motif:
    	motifsRegex = sequence.Regexp(m)
    	motifObj.append(motifsRegex)	

    return motifObj


if __name__ == '__main__':
	main()