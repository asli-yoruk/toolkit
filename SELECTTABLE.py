import argparse
import sequence
import csv

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", help="Input FASTA file", required=True)
	parser.add_argument("-req", "--reqr", nargs="*", help="Required Motifs/info")
	parser.add_argument("-o", "--output", help="Output path", default="selecttable.csv")

	args = parser.parse_args()
	equals = []
	notEquals = []
	betweens = []

	outputfile = output(args)

	for req in args.reqr:
		if "==" in req:
			equals.append(req.split("=="))
			#print(equals)
		elif "!=" in req:
			notEquals.append(req.split("!="))
			#print(notEquals)
		elif "~" in req:
			betweens.append(req.split("~"))
			#print(betweens)

	colsToSearch = []
	lst = []

	with open(args.input, newline='') as f:
		reader = csv.reader(f)
		header = next(reader)  
		lst.append(header)

		for h in range (len(header)):
			for e in equals:
				if e[0] == header[h]:
					colsToSearch.append(str(h)+','+e[1])

		for row in reader:
			#print("this is the row " + str(row))
			for c in colsToSearch:
				col, pos = c.split(',')
				elems = (row[(int(col))]).split(",")

				if len(elems) > 1:
					for i in range (len(elems)):
						elems[i] = int(elems[i])

				elif len(elems) == 1 and elems[0] != '':
					elems[0] = int(elems[0])

				#print(elems)
				if int(pos) in elems:
					take = True
				else:
					take = False
					break
			if take:
				lst.append(row)


	with open(outputfile, 'w', newline='') as f:
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
	    		thewriter.writerow(motif_pos)


def output(args):
	if args.output != "motifsearch.csv":
		output = args.output
	else:
		output = "motifsearch.csv"
	return output

if __name__ == '__main__':
	main()