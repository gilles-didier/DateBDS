Divergence time distribution and sampling

Two software are provided
 - 'dist'
	compute the divergence time distributions of a phylogenetic tree under the birth-death-sampling model.
 - 'samp'
	sample the divergence times of a phylogenetic tree under the birth-death-sampling model.

type
	> make all
	in a console opened on the src directory for compiling the software.

Directory "src" contains the C sources of the software and Directory "data" contains 3 examples of trees.

In a console opened in the directory "src", you could write

./dist -f 1 -o 0 10 -p 0 0.1 0.01 0.1 -p 5 0.9 0.1 0.7 -e -d 500 ../data/A.phy

./samp -o 0 10 -p 0 0.1 0.01 0.1 -p 5 0.9 0.1 0.7 ../data/A.phy

./asse -a 0.6 -b 1. -s 5  -t 10 -i 25000 -m 10 -M 50000 -x 8 Test


A complete description of the options allowed is given below.

--------
| dist |
--------

--------------------------
REQUIREMENT

	The software requires the GSL and Cairo libraries.

--------------------------
COMPILING

	Just type
	> make dist
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'dist' compute the divergence time distributions of a phylogenetic tree under the birth-death-sampling model.


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	dist - divergence time distributions of a phylogenetic tree under the birth-death-sampling model.
	
SYNOPSIS
	dist [OPTIONS] <input Tree File> <output Ident>

DESCRIPTION
	Compute the divergence time distribution of the phylogenetic tree contained in <input Tree File> (which must be in Newick format) with regards  to the parameters provided as options and output the results as array/text files '.csv' and as a figure if the option '-f' is set. Result files are save in the working directory.

	Options are
	-o <origin time> <end time>
		set the origin and end time of the diversification process resulting to the phylogenetic tree. 
	-p <start time> <speciation rate> <extinction rate> <sampling probability>
		set the parameters of the piece starting at <starting time> of a piecewise-constant-birth-death-sampling model. 
	-e
		display probability densities (distributions otherwise).
	-z <input Tree File>
		output the tree in 'text' format in the console and save it in pdf format as 'tree.pdf' with the internal idents of the nodes (for debug purposes) next exit
	-n <node number>
		compute distribution of only one node (use -z option to see the node numbers).
	-d <number>
		set the number of points of the distributions computed.
	-f <number>
		set the graphic format of the output (option is required if one wants a graphic output)
			-f 1 -> pdf
			-f 2 -> postscript
			-f 3 -> png
			-f 4 -> svg
			-f 5 -> LaTeX (psTricks)
			-f 6 -> LaTeX (TikZ)
	-c <r1> <g1> <b1> <r2> <g2> <b2>
		set the color scale to go from the color (r1,g1,b1) to the color (r2,g2,b2) (rgb codes of colors have their components in [0,1])
	-h
		display help

--------------------------


--------
| samp |
--------

--------------------------
REQUIREMENT

	The software requires the GSL library.

--------------------------
COMPILING

	Just type
	> make samp
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'samp' fix data to be used by ParSplit.


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	samp - divergence times sampling under the birth-death-sampling model
	
SYNOPSIS
	samp [OPTIONS] <input Tree File> <output Ident>

DESCRIPTION
	Sample the divergence time distribution of the phylogenetic tree contained in <input Tree File> (which must be in Newick format) with regards  to the parameters provided as options and output a dated tree in Newick format (with extension ".phy"). Result files are saved in the directory of the input file.

	Options are
	-o <origin time> <end time>
		set the origin and end time of the diversification process resulting to the phylogenetic tree. 
	-p <start time> <speciation rate> <extinction rate> <sampling probability>
		set the parameters of the piece starting at <starting time> of a piecewise-constant-birth-death-sampling model. 
	-z <input Tree File>
		output the tree in 'text' format in the console and save it in pdf format as 'tree.pdf' with the internal idents of the nodes (for debug purposes) next exit
	-e <value>
		set the precision used for sampling.
	-h
		display help

--------------------------


--------
| shif |
--------

--------------------------
REQUIREMENT

	The software requires the GSL and the NLOpt libraries.

--------------------------
COMPILING

	Just type
	> make shif
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'shif' computes ROC plots for 4 tests of diversification shifts from Yule tree simulations


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	shif - ROC plots of diversification shift tests
	
SYNOPSIS
	shif [OPTIONS] <output Ident>

DESCRIPTION
	Compute ROC plots for 4 tests of diversification shifts from Yule tree simulations. Results are returned as 4 table files with names "ROC_?_<output Ident>.csv" where columns 1 and 2 is the true and the false positive rates respectively. 

	Options are
	-o <NLopt option file>
		read the numerical optimizers parameters from  <NLopt option file>
	-a <birth rate>
		set the general birth rate
	-b <birth rate>
		set the shifted birth rate
	-s <time>
		set the shift time
	-t <time>
		set the total diversification time (it starts from time 0)
	-i <number>
		set the number of simulations
	-x <number>
		set the number of simultaneous threads in parallel
	-m <number>
		set the minimum number of tips required
	-M <number>
		set the maximum number of nodes allowed
	-h
		display help

--------------------------

