import argparse
import sequence
import csv
import pandas as pd

def main():
	args = getInput()

	if args.input and args.annotateFile:
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


#Force the user to give CSV for merging tables
def annotateMerge(args):
	outputfile = output(args)

	# commonCols = []

	# if '.fa' in args.input:
	# 	args.input = 'table1.csv'
	# 	db100 = sequence.readFastaFile(args.input, sequence.Protein_Alphabet, ignore=True, parse_defline=False)

	# 	with open('table1.csv', 'w', newline='') as f:
	# 	    fieldnames = ['Name', 'Sequence']
	# 	    thewriter = csv.DictWriter(f, fieldnames=fieldnames)
		    
	# 	    thewriter.writeheader()
	# 	    for seq in db100:
	# 	    	s = ''.join(seq.sequence)
	# 	    	thewriter.writerow({'Name' : seq.name, 'Sequence' : s})



	# if '.fa' in args.annotateFile:
	# 	print("this is a FASTA file")
	# 	annot = sequence.readFastaFile(args.annotateFile, sequence.Protein_Alphabet, ignore=True, parse_defline=False)
	# 	args.annotateFile = "table2.csv"
	# 	with open("table2.csv", 'w', newline='') as f:
	# 	    fieldnames = ['Name', 'Sequence', args.annotateKeyword]
	# 	    thewriter = csv.DictWriter(f, fieldnames=fieldnames)
		    
	# 	    thewriter.writeheader()
	# 	    for seq in annot:
	# 		    	s = ''.join(seq.sequence)
	# 		    	thewriter.writerow({'Name' : seq.name, 'Sequence' : s, args.annotateKeyword : args.annotateKeyword})


	if '.csv' in args.annotateFile:
		print("this is csv")


	a = pd.read_csv(args.input)
	f = pd.read_csv(args.annotateFile)
	cols = []
	print("First table")
	print(a)
	print("Second table")
	print(f)

	for col in a.columns:
	    if (col != 'Name') and (col != 'Sequence'):
	        cols.append(col)

	for col in f.columns:
	    if (col != 'Name') and (col != 'Sequence'):
	        cols.append(col)

	print (cols)

	merged = pd.merge(a, f,  how='outer', on = ['Name','Sequence'])
	print(merged)

	if len(cols) == 1:
		merged.to_csv("annot_file.csv", index=False)
	else:
		cols = []


		for col in merged.columns:
		    if (col != 'Name') and (col != 'Sequence'):
		        cols.append(col)


		print(cols)
		merged['combined'] = merged[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)


		merged = merged.drop(cols,1)
		print(merged)

		setAnnots = []

		for v in merged['combined']:
			mergedannots = ""
			v = v.split("_")
			annots = []
			for i in range (len(v)):
			    if v[i] not in annots:
			        if v[i] != 'nan':                
			            annots.append(v[i])
			count = 0
			mergedannots = annots[0]
			for x in annots:
			    if count != 0:
			        mergedannots += '_'+x
			    count += 1
			setAnnots.append(mergedannots)
		    # v = v.split("_")
		    # annots = []
		    # for i in range (len(v)):
		    #     if v[i] not in annots:
		    #         if v[i] != 'nan':                
		    #             annots.append(v[i])
		    # setAnnots.append(annots)

		    
		merged = merged.drop("combined",1)
		merged['Annots'] = setAnnots
		merged
		merged.to_csv("new_annot_file.csv", index=False)



def annotateThis(args):

	outputfile = output(args)

	if '.fa' in args.annotateFile:
		print("this is a FASTA file")
		annot = sequence.readFastaFile(args.annotateFile, sequence.Protein_Alphabet, ignore=True, parse_defline=False)

		with open(outputfile, 'w', newline='') as f:
		    fieldnames = ['Name', 'Sequence', args.annotateKeyword]
		    thewriter = csv.DictWriter(f, fieldnames=fieldnames)
		    
		    thewriter.writeheader()
		    for seq in annot:
			    	s = ''.join(seq.sequence)
			    	thewriter.writerow({'Name' : seq.name, 'Sequence' : s, args.annotateKeyword : args.annotateKeyword})

	elif '.csv' in args.annotateFile:
		print("this is a CSV file")
		with open(args.annotateFile, newline='') as f:
			reader = csv.reader(f)
			header = next(reader)  
			nameCol = 0
			seqCol = 0
			dictionary = {}

			for h in range (len(header)):
				if header[h] == 'Name':
					nameCol = h
				elif header[h] == 'Sequence':
					seqCol = h

			for row in reader:
				dictionary[row[nameCol]] = row[seqCol]

		
		with open(outputfile, 'w', newline='') as f:
			fieldnames = ['Name', 'Sequence', 'Annots']
			thewriter = csv.DictWriter(f, fieldnames=fieldnames)
			thewriter.writeheader()
			for name,seq in dictionary.items():
				thewriter.writerow({'Name' : name, 'Sequence' : seq, 'Annots' : args.annotateKeyword})




def output(args):
	if args.output != 'table_output.csv':
		output = args.output
	else:
		output = "table_output.csv"
	return output


if __name__ == "__main__":
    main()

