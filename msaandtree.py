import os

redundancy = [99,98,95,90,80,70]
def main():
        for rr in redundancy:
                for evalue in [1e-100, 1e-75, 1e-50, 1e-20, 1e-10, 1e-5]:
                        filename =  str(rr) + '-'+ str(evalue)
                        os.system('clustalo -i dataset-' + filename + '_reduced.fa -o dataset-' + filename + '.aln --output-order=tree-order --threads=5 --force')
        for rr in redundancy:
                for evalue in [1e-100, 1e-75, 1e-50, 1e-20, 1e-10, 1e-5]:
                        filename =  str(rr) + '-'+ str(evalue)
                        os.system("FastTree -out dataset-" + filename + ".nwk dataset-" + filename + ".aln")

if __name__ == "__main__":
    main()
