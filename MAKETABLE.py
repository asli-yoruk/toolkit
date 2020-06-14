import argparse
import sequence
import csv
import pandas as pd

def main():
	args = getInput()

	#if annotateFile and annotateKeyword exists in the input
	#input file extension MUST be a csv
	#TODO not handling every case
	if args.input and args.annotateFile and args.annotateKeyword:
		annotateMerge(args)
	elif args.annotateKeyword and args.annotateFile:
		annotateThis(args)
	else:
		read(args)



def getInput():
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", help="Input FASTA file")
	parser.add_argument("-ak", "--annotateKeyword", help="Annotation keyword to add")
	parser.add_argument("-af", "--annotateFile", help="Annotation file")
	parser.add_argument("-o", "--output", help="Output path", default="table_output.csv")

	args = parser.parse_args()
	return args


def read(args):

	outputfile = output(args)

	db100 = sequence.readFastaFile(args.input, sequence.Protein_Alphabet, ignore=True, parse_defline=False)

	with open(outputfile, 'w', newline='') as f:
	    fieldnames = ['Name', 'Sequence']
	    thewriter = csv.DictWriter(f, fieldnames=fieldnames)
	    
	    thewriter.writeheader()
	    for seq in db100:
	    	s = ''.join(seq.sequence)
	    	thewriter.writerow({'Name' : seq.name, 'Sequence' : s})


#IF GIVEN ANNOTATION FILE IS A CSV, go through the header names 
#and write a csv file accordingly
def annotateMerge(args):
	outputfile = output(args)

	commonCols = []

	if '.fa' in args.input:
		db100 = sequence.readFastaFile(args.input, sequence.Protein_Alphabet, ignore=True, parse_defline=False)

		with open('input.csv', 'w', newline='') as f:
		    fieldnames = ['Name', 'Sequence']
		    thewriter = csv.DictWriter(f, fieldnames=fieldnames)
		    
		    thewriter.writeheader()
		    for seq in db100:
		    	s = ''.join(seq.sequence)
		    	thewriter.writerow({'Name' : seq.name, 'Sequence' : s})
		args.input = 'input.csv'


	if '.fa' in args.annotateFile:
		annot = sequence.readFastaFile(args.annotateFile, sequence.Protein_Alphabet, ignore=True, parse_defline=False)
		with open('annot.csv', 'w', newline='') as f:
		    fieldnames = ['Name', 'Sequence', args.annotateKeyword]
		    thewriter = csv.DictWriter(f, fieldnames=fieldnames)
		    
		    thewriter.writeheader()
		    for seq in annot:
			    	s = ''.join(seq.sequence)
			    	thewriter.writerow({'Name' : seq.name, 'Sequence' : s, args.annotateKeyword : args.annotateKeyword})
		headerAnnot = fieldnames

		with open(args.input, newline='') as f:
			reader = csv.reader(f)
			headerInput = next(reader)  # gets the first line

		for h1 in headerInput:
			for h2 in headerAnnot:
				if h1 == h2:
					commonCols.append(h1)

		a = pd.read_csv(args.input)
		b = pd.read_csv('annot.csv')
		merged = pd.merge(a, b,  how='outer', left_on=commonCols, right_on = commonCols)
		merged.to_csv(outputfile, index=False)

	elif '.csv' in args.annotateFile:
		print("this is csv")

		with open(args.annotateFile, newline='') as f:
			reader = csv.reader(f)
			headerAnnot = next(reader)  

		with open(args.input, newline='') as f:
			reader = csv.reader(f)
			headerInput = next(reader)  

		for h1 in headerInput:
			for h2 in headerAnnot:
				if h1 == h2:
					commonCols.append(h1)

		a = pd.read_csv(args.input)
		b = pd.read_csv(args.annotateFile)
		merged = pd.merge(a, b,  how='outer', left_on=commonCols, right_on = commonCols)
		merged.to_csv(outputfile, index=False)


def annotateThis(args):
	outputfile = output(args)

	annot = sequence.readFastaFile(args.annotateFile, sequence.Protein_Alphabet, ignore=True, parse_defline=False)

	with open(outputfile, 'w', newline='') as f:
	    fieldnames = ['Name', 'Sequence', args.annotateKeyword]
	    thewriter = csv.DictWriter(f, fieldnames=fieldnames)
	    
	    thewriter.writeheader()
	    for seq in annot:
		    	s = ''.join(seq.sequence)
		    	thewriter.writerow({'Name' : seq.name, 'Sequence' : s, args.annotateKeyword : args.annotateKeyword})



def output(args):
	if args.output != 'table_output.csv':
		output = args.output
	else:
		output = "table_output.csv"
	return output


if __name__ == "__main__":
    main()

