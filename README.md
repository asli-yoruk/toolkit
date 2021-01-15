# toolkit

	⁃	Make sure that the desired applications are installed on your machine.
		- BLAST+ (Please don't skip the configuration step which has been mentioned
		in the https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download .
		Find BLAST Help manual there and click it. If your machine is Windows then go to Standalone BLAST
		Setup for Windows PC, https://www.ncbi.nlm.nih.gov/books/NBK52637/ ,if then go to Standalone BLAST
		Setup for Unix https://www.ncbi.nlm.nih.gov/books/NBK52640/?fbclid=IwAR0O3X6ZfMdlojxPcxLFDn02QrpHau4QbICznOsFwNh2weBajAkklRVRJtU. 
		And make sure you do the Configuration which has been mentioned there.)
		- CD-HIT (commandline tool; https://github.com/weizhongli/cdhit/wiki); install is easy if Anaconda has been installed earlier -> conda install -c bioconda cd-hit
		- Clustal Omega
		- Python version 3 or higher (you can check your version by simply typing python --version on your command line)
	⁃	Please keep all of the files in the same directory. (script.sh, all of the python files .py, all of the fasta files .fa)
	
	Other programs that are needed to do further work:
		Sequence alignment program, recommend Clustal Omega (commandline tool; http://www.clustal.org/omega/) or MAFFT (commandline https://mafft.cbrc.jp/alignment/software/)

		Alignment viewer, recommend AliView (http://ormbunkar.se/aliview/)

		Phylogenetic tree inference program, recommend FastTree (commandline tool; http://www.microbesonline.org/fasttree/)

		Phylogenetic tree viewer, recommend FigTree (well-designed; http://tree.bio.ed.ac.uk/software/figtree/) or Archaeopteryx (for large trees; https://sites.google.com/site/cmzmasek/home/software/archaeopteryx)

		Below, instructions assume a UNIX like environment; MacOS is one and so is Linux, but Windows is not. Cygwin is a suite of tools that installs under Windows to provide a lot of the same command-like functionality.
  - Happy coding.
