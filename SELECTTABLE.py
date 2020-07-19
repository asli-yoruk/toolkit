import argparse
import sequence
import csv

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input",nargs="*", help="Input FASTA file", required=True)
	parser.add_argument("-req", "--reqr", nargs="*", help="Required Motifs/info")
	#parser.add_argument("-o", "--output", help="Output path", default="selecttable.csv")

	args = parser.parse_args()
	equals = []
	notEquals = []
	betweens = []

	#outputfile = output(args)

	for req in args.reqr:
		if "==" in req:
			equals.append(req.split("=="))
			posEqual_csv(args, equals)
			posEqual_fa(args)

		elif "!=" in req:
			notEquals.append(req.split("!="))

		elif "~" in req:
			betweens.append(req.split("~"))

		elif "all":
			print("All motifs must be present")
			allMotifs_fa(args)
			allMotifs_csv(args,equals)

		else:
			print("Your requirement input is wrong")


def posEqual_csv(args, equals):

	for i in range (len(args.input)):
		name = (args.input[i])
		name = name.split('.')[0]
		name = name + '_reduced.csv'

		colsToSearch = []
		lst = []
		take = True

		with open(args.input[i], newline='') as f:
			reader = csv.reader(f)
			header = next(reader)  
			lst.append(header)

			for h in range (len(header)):
				for e in equals:
					if e[0] == header[h]:
						colsToSearch.append(str(h)+','+e[1])

			for row in reader:
				for c in colsToSearch:

					col, pos = c.split(',')
					elems = (row[(int(col))]).split(",")

					if len(elems) > 1:
						for i in range (len(elems)):
							elems[i] = int(elems[i])

					elif len(elems) == 1 and elems[0] != '':
						elems[0] = int(elems[0])

					
					if int(pos) in elems:
						take = True
					else:
						take = False
						break
				if take:
					lst.append(row)

		with open(name, 'w', newline='') as f:
			c=0
			fieldnames = lst[0]
			thewriter = csv.DictWriter(f, fieldnames=fieldnames)
			thewriter.writeheader()
			first = True
			for ab in lst:	
				motif_pos = dict()
				count = 0
				if first == True:
					first = False	
				else:	
					for a in ab:
						motif_pos[lst[0][count]] = a
						count += 1
					c+=1
					thewriter.writerow(motif_pos)
		print(str(c) + " sequences kept after applying the requirements for " + name.split('.')[0])

		


def allMotifs_csv(args,equals):
	for i in range (len(args.input)):
		name = (args.input[i])
		name = name.split('.')[0]
		name = name + '_reduced.csv'

		colsToSearch = []
		lst = []
		take = True

		with open(args.input[i], newline='') as f:
			reader = csv.reader(f)
			header = next(reader)  
			lst.append(header)

			for h in range (len(header)):
				for e in equals:
					if e[0] == header[h]:
						colsToSearch.append(str(h)+','+e[1])

			for row in reader:
				for c in colsToSearch:

					col, pos = c.split(',')
					elems = (row[(int(col))]).split(",")

					if len(elems) > 1:
						for i in range (len(elems)):
							elems[i] = int(elems[i])

					elif len(elems) == 1 and elems[0] != '':
						elems[0] = int(elems[0])

					
					if int(pos) in elems:
						take = True
					else:
						take = False
						break
				if take:
					lst.append(row)



		with open(name, 'w', newline='') as f:
			c = 0
			fieldnames = lst[0]
			thewriter = csv.DictWriter(f, fieldnames=fieldnames)
			thewriter.writeheader()
			first = True
			for ab in lst:	
				motif_pos = dict()
				count = 0
				if first == True:
					first = False	
				else:
					empty = False
					for a in ab:
						motif_pos[lst[0][count]] = a
						count += 1
						if a == "":
							empty = True
					if empty == False:
						c+=1
						thewriter.writerow(motif_pos)
			print(str(c) + " sequences kept after applying the requirements for " + name)
			#print(str(c) + " sequences kept after applying the requirements for " + name.split('.')[0])



def allMotifs_fa(args):
	#check hoow many cols, check if all of them has a value
	#make a FASTA and make a CSV

	
	for i in range (len(args.input)):
		c=0
		fasta = {}
		name = (args.input[i])
		name = name.split('.')[0]
		name = name + '_reduced.fa'

		with open(args.input[i], newline='') as f:
			reader = csv.reader(f)
			header = next(reader)

			for row in reader:
				isEmpty = False
				for i in range(1, len(header)-1):
					if row[i] == "":
						isEmpty = True
						break
				if isEmpty == False:
					fasta[row[0]] = row[len(header)-1]

		seq_list = [sequence.Sequence(sequence=seq, name=seqname) for seqname, seq in fasta.items()]
		sequence.writeFastaFile(name, seq_list)
		c+=1
		print(str(len(seq_list)) + " sequences kept after applying the requirements for " + name)


def posEqual_fa(args):
	fasta = {}
	
	for i in range (len(args.input)):
		name = (args.input[i])
		name = name.split('.')[0]
		name = name + '_reduced.fa'

		with open(args.input[i], newline='') as f:
			reader = csv.reader(f)
			header = next(reader)

			for row in reader:
				isEmpty = False
				for i in range(1, len(header)-1):
					if row[i] == "":
						isEmpty = True
						break
				if isEmpty == False:
					fasta[row[0]] = row[len(header)-1]

		seq_list = [sequence.Sequence(sequence=seq, name=seqname) for seqname, seq in fasta.items()]
		sequence.writeFastaFile(name, seq_list)


if __name__ == '__main__':
	main()