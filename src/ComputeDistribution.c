#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include "Utils.h"
#include "Tree.h"
#include "TreeBounds.h"
#include "SimulTree.h"
#include "Model.h"
#include "Uncertainty.h"
#include "Distribution.h"
#include "DrawTreeCairo.h"
#include "DrawTreePSTricks.h"
#include "DrawTreeTikz.h"
#include "DrawTreeGeneric.h"
#include "DrawDensity.h"

#ifdef DO_PS
#endif

#define STRING_SIZE 300
#define INC_CUT 5
#define SEQ_SIZE 30
#define EXT_OUTPUT "_added.phy"
#define MAX_PRED 7.5E+07

#define SIZE_BUFFER_CHAR 300
#define INFTY 1E99
#define RINFTY 1E99
#define DEF 10
#define MIN_VAL 0.000001
#define DELTA 0.000001

#define MINVAL 0.01
#define TRIAL 10
#define FIX_VAL(x) (((x)<=0)?MINVAL:(x))

#define MAX_ITER 1000

#define M_MAX 6
#define M_MAX_F 4
#define MIN_TREE 20
#define MAX_TRIALS 1000
#define MAX_PIECES 20

#define PREFIX "table"

#define HELPMESSAGE "NAME\n	dist - divergence time distributions of a phylogenetic tree under the birth-death-sampling model.\n	\nSYNOPSIS\n	dist [OPTIONS] <input Tree File> <output Ident>\n\nDESCRIPTION\n	Compute the divergence time distribution of the phylogenetic tree contained in <input Tree File> (which must be in Newick format) with regards  to the parameters provided as options and output the results as array/text files '.csv' and as a figure if the option '-f' is set. Result files are save in the working directory.\n\n	Options are\n	-o <origin time> <end time>\n		set the origin and end time of the diversification process resulting to the phylogenetic tree. \n	-p <start time> <speciation rate> <extinction rate> <sampling probability>\n		set the parameters of the piece starting at <starting time> of a piecewise-constant-birth-death-sampling model. \n	-e\n		display probability densities (distributions otherwise).\n	-z <input Tree File>\n		output the tree in 'text' format in the console and save it in pdf format as 'tree.pdf' with the internal idents of the nodes (for debug purposes) next exit\n	-n <node number>\n		compute distribution of only one node (use -z option to see the node numbers).\n	-d <number>\n		set the number of points of the distributions computed.\n	-f <number>\n		set the graphic format of the output (option is required if one wants a graphic output)\n			-f 1 -> pdf\n			-f 2 -> postscript\n			-f 3 -> png\n			-f 4 -> svg\n			-f 5 -> LaTeX (psTricks)\n			-f 6 -> LaTeX (TikZ)\n	-c <r1> <g1> <b1> <r2> <g2> <b2>\n		set the color scale to go from the color (r1,g1,b1) to the color (r2,g2,b2) (rgb codes of colors have their components in [0,1])\n	-h\n		display help\n"
//./dist -p 0.1 0.02 0.5 -e -o 0 10 -d 1000 ../figures/trees/treeA.newick 

//./dist -p 0.1 0.02 0.5 TEST4.txt 

int main(int argc, char **argv) {	
	char inputFileNameTree[STRING_SIZE], *outputPrefix = PREFIX, outputDistribution[STRING_SIZE], option[256], format = '1';
	FILE *fi, *fo;
	int i, j, node = NOSUCH, def=10, outDens = 0, max_size = INT_MAX;
	double val = NEG_INFTY,  contemp = 0.,  origin = -100., *startTmp;
	TypeModelParam *paramTmp;
	TypePiecewiseModelParam piece;
	TypeTree *tree;
	size_t *pieceIndex;
	
	piece.size = 0;
	piece.startTime = (double*) malloc(MAX_PIECES*sizeof(double));
	piece.param = (TypeModelParam*) malloc(MAX_PIECES*sizeof(TypeModelParam));
	
	for(i=0; i<256; i++)
		option[i] = 0;
	   
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &format) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a character is required after option -f");
		}
		if(option['z']) {
			TypeInfoDrawTreeGeneric info;
			option['z'] = 0;
			if((i+1)<argc && (fi = fopen(argv[++i], "r"))) {
				TypeTree *tree;
				tree = readTree(fi);
				printTreeDebug(stdout, tree->root, tree, tree->name);
//				bltoabsTime(tree);
				reorderTree(tree->name, tree);
				if(tree->minTime == NO_TIME || tree->minTime == 0.)
					tree->minTime = tree->time[tree->root]*0.9;
				if(tree->maxTime == NO_TIME) {
					int n;
					tree->maxTime = 0.;
					for(n=0; n<tree->size; n++)
						if(tree->time[n]>tree->maxTime)
							tree->maxTime = tree->time[n];
				}
				setFunctPDF(&(info.funct));
				drawTreeFileGenericDebug("tree.pdf", tree, &info, NULL);
			} else {
				fprintf(stderr, "Error while reading tree.\n");
				exit(1);
			}
			exit(0);
		}
		if(option['o']) {
			option['o'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &(origin)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "two times are expected after -o");
			if((i+1)<argc && sscanf(argv[i+1], "%le", &(contemp)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "two times are expected after -o");
		}
		if(option['p']) {
			option['p'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &(piece.startTime[piece.size])) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -p");
			if((i+1)<argc && sscanf(argv[i+1], "%le", &(piece.param[piece.size].birth)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -p");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(piece.param[piece.size].death)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -p");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(piece.param[piece.size].sampl)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -p");
			piece.size++;
		}
		if(option['e']) {
			option['e'] = 0;
			outDens = 1;
		}
		if(option['n']) {
			option['n'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &(node)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "1 node index are expected after -n");
		}
		if(option['m']) {
			option['m'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &(max_size)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "1 node index are expected after -n");
		}
		if(option['d']) {
			option['d'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &(def)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "1 number are expected after -d");
		}
		if(option['v']) {
			option['v'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &(val)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "1 number are expected after -v");
		}
		if(option['h']) {
			option['h'] = 0;
			printf("%s\n", HELPMESSAGE);
			exit(0);
		}
	}
	if(!(i<argc && sscanf(argv[i++], "%s", inputFileNameTree) == 1)) {
		fprintf(stderr, "%s\n",  HELPMESSAGE);
		exit(1);
	}
	if(i<argc)
		outputPrefix = argv[i++];
	pieceIndex = qsortTable(piece.startTime, piece.size, sizeof(double), compareDouble);
	startTmp = (double*) malloc((piece.size+1)*sizeof(double));
	paramTmp = (TypeModelParam*) malloc(piece.size*sizeof(TypeModelParam));
	for(i=0; i<piece.size; i++) {
		startTmp[pieceIndex[i]] = piece.startTime[i];
		paramTmp[pieceIndex[i]] = piece.param[i];
	}
	startTmp[0] = origin;
	startTmp[piece.size] = contemp;
	free((void*)pieceIndex);
	free((void*)piece.startTime);
	free((void*)piece.param);
	piece.startTime = startTmp;
	piece.param = paramTmp;
    if((fi = fopen(inputFileNameTree, "r"))) {
		TypeDistribution *d;
		TypeInfoDrawTreeGeneric info;
		TypeAdditionalDrawTreeGeneric add;
		TypeDataDrawDensity dataD;		
		double *timeSave, *timeSaveX, *tmin, *tmax;
		int n, *todo;
        tree = readTree(fi);
        fclose(fi);
        if(tree == NULL)
            return 1;
        toBinary(tree);
//        reorderTree(tree->name, tree);
        reorderTree(NULL, tree);
		tree->minTime = origin;
		tree->maxTime = contemp;
		timeSave = tree->time;
		tree->time = (double*) malloc(tree->size*sizeof(double));
		for(n=0; n<tree->size; n++) {
			if(tree->node[n].child == NOSUCH)
				tree->time[n] = contemp;
			else
				tree->time[n] = NO_TIME;
		}
		tmin = (double*) malloc(tree->size*sizeof(double));
		tmax = (double*) malloc(tree->size*sizeof(double));
		fillBoundsFromComments(tmin, tmax, tree);
		//for(n=0; n<tree->size; n++)
			//printf("node %d min %.2lf max %.2lf\n", n, (tmin[n]!=NO_TIME)?tmin[n]:sqrt(-1), (tmax[n]!=NO_TIME)?tmax[n]:sqrt(-1));
		todo = (int*) malloc(tree->size*sizeof(int));
		if(node != NOSUCH && node>=tree->size)
			node = NOSUCH;
		for(n=0; n<tree->size; n++)
			todo[n] = 0;
		if(node == NOSUCH) {
			for(n=0; n<tree->size; n++)
				if(tree->node[n].child != NOSUCH)
					todo[n] = 1;
		} else
			if(tree->node[node].child != NOSUCH)
				todo[node] = 1;
		d = (TypeDistribution*) malloc(tree->size*sizeof(TypeDistribution));
		fillDistributionAllGeneral(def, tmin, tmax, tree, &piece, todo, d);
//		fillDistributionAll(def, tree, &param, todo, d);
//printf("\nOK\n");
//printf("OK dist %d\n", d[5].size);
//deriveDistribution(&(d[5]));
//if((fo = fopen("out.csv", "w"))) {
//fprintDistribution(fo, d[5]);
//fclose(fo);
//}
//exit(0);
		for(n=0; n<tree->size; n++) {
			if(todo[n]) {
				sprintf(outputDistribution, "%s%d.csv", outputPrefix, n);
				if((fo = fopen(outputDistribution, "w"))) {
					if(outDens)
						fprintDerive(fo, d[n]);
					else
						fprintDistribution(fo, d[n]);
					fclose(fo);
				}
			} else {
				d[n].size = 0;
				d[n].item = NULL;
			}
		}
		timeSaveX = tree->time;
		tree->time = (double*) malloc(tree->size*sizeof(double));
		for(n=0; n<tree->size; n++) {
			tree->time[n] = timeSaveX[n];
			if(tree->time[n] == NO_TIME && d[n].size != 0)
				tree->time[n] = getMedian(d[n]);
		}
		for(n=0; n<tree->size; n++)
			deriveDistribution(&(d[n]));
//		fprintTreeNewick(stdout, tree);
		if(format != '0') {
			char outputFileNameG[SIZE_BUFFER_CHAR], *outputFileName;
//			outputFileName = inputFileNameTree;
			outputFileName = outputPrefix;
			//if((tmp = strrchr(outputFileName, '.')) != NULL)
				//tmp[0] = '\0';
			info.param.tmin = tree->minTime;
			info.param.tmax = tree->maxTime;
//			dataD.color = (TypeRGB) {.red = 1., .green = 0., .blue = 0.};
			dataD.color = (TypeRGB) {.red = 0.5, .green = 0.5, .blue = 0.5};
			dataD.alpha = 0.5;
			dataD.dens = d;
			add.data = &dataD;
			add.draw = drawDensity;
			switch(format) {
				case '1':
					sprintf(outputFileNameG, "%s_tree.pdf", outputFileName);
					setFunctPDF(&(info.funct));
					dataD.fillPolygon = fillPolygonCairo;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				case '2':
					sprintf(outputFileNameG, "%s_tree.ps", outputFileName);
					setFunctPS(&(info.funct));
					dataD.fillPolygon = fillPolygonCairo;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				case '3':
					sprintf(outputFileNameG, "%s_tree.png", outputFileName);
					setFunctPNG(&(info.funct));
					dataD.fillPolygon = fillPolygonCairo;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				case '4':
					sprintf(outputFileNameG, "%s_tree.svg", outputFileName);
					setFunctSVG(&(info.funct));
					dataD.fillPolygon = fillPolygonCairo;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				case '5':
					sprintf(outputFileNameG, "%s_tree_pst.tex", outputFileName);
					setFunctPSTricks(&(info.funct));
					dataD.fillPolygon = fillPolygonPSTricks;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				case '6':
					sprintf(outputFileNameG, "%s_tree_tikz.tex", outputFileName);
					setFunctTikz(&(info.funct));
					dataD.fillPolygon = fillPolygonTikz;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				default:
					;
			}
		}
		free((void*) tree->time);
		tree->time = timeSaveX;
		fprintf(stdout,"\nplot");
		for(n=0; n<tree->size; n++) {
			if(todo[n]) {
				if(tree->name != NULL && tree->name[n] != NULL)
					fprintf(stdout," '%s%d.csv' using 1:2 with lines title 'node %s',", outputPrefix, n, tree->name[n]);
				else
					fprintf(stdout," '%s%d.csv' using 1:2 with lines title 'node %d',", outputPrefix, n, n);
			}
		}
		fprintf(stdout,"\n\n");
		if(d != NULL) {
			for(n=0; n<tree->size; n++)
				if(d[n].item != NULL)
					free((void*)d[n].item);
			free((void*)d);
		}
		free((void*)todo);
		free((void*)tree->time);
		tree->time = timeSave;
		freeTree(tree);
	} else {
		fprintf(stderr, "Cannot read %s\n", inputFileNameTree);
		exit(1);
	}
	return 0;
}
