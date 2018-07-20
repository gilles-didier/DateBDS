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

#define PREFIX "table"

#define HELPMESSAGE "NAME\n\tdist - divergence time distributions of a phylogenetic tree under the birth-death-sampling model.\n\t\nSYNOPSIS\n\tdist [OPTIONS] <input Tree File> <output Ident>\n\nDESCRIPTION\n\tCompute the divergence time distribution of the phylogenetic tree contained in <input Tree File> (which must be in Newick format) with regards  to the parameters provided as options and output the results as array/text files '.csv' and as a figure if the option '-f' is set.\n\n\tOptions are\n\t-o <origin time> <end time>\n\t\tset the origin and end time of the diversification process resulting to the phylogenetic tree. \n\t-p <speciation rate> <extinction rate> <sampling probability>\n\t\tset the parameters of the birth-death-sampling model. \n\t-e\n\t\tdisplay probability densities (distributions otherwise).\n\t-z <input Tree File>\n\t\toutput the tree in 'text' format in the console and save it in pdf format as 'tree.pdf' with the internal idents of the nodes (for debug purposes) next exit\n\t-n <node number>\n\t\tcompute distribution of only one node (use -z option to see the node numbers).\n\t-d <number>\n\t\tset the number of points of the distributions computed.\n\t-f <number>\n\t\tset the graphic format of the output (option is required if one wants a graphic output)\n\t\t\t-f 1 -> pdf\n\t\t\t-f 2 -> postscript\n\t\t\t-f 3 -> png\n\t\t\t-f 4 -> svg\n\t\t\t-f 5 -> LaTeX (psTricks)\n\t\t\t-f 6 -> LaTeX (TikZ)\n\t-c <r1> <g1> <b1> <r2> <g2> <b2>\n\t\tset the color scale to go from the color (r1,g1,b1) to the color (r2,g2,b2) (rgb codes of colors have their components in [0,1])\n\t-h\n\t\tdisplay help\n"

//./dist -p 0.1 0.02 0.5 -e -o 0 10 -d 1000 ../figures/trees/treeA.newick 

//./dist -p 0.1 0.02 0.5 TEST4.txt 

int main(int argc, char **argv) {	
	char inputFileNameTree[STRING_SIZE], *outputPrefix = PREFIX, outputDistribution[STRING_SIZE], option[256], format = '1';
	FILE *fi, *fo;
	int i, j, node = NOSUCH, def=10, outDens = 0, max_size = INT_MAX;
	double contemp = 0.,  origin = -100.;
	TypeModelParam param = {.birth=0.3, .death = 0.1, .sampl=1.};
	TypeTree *tree;
	
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
			if((i+1)<argc && sscanf(argv[i+1], "%le", &(param.birth)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -p");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(param.death)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -p");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(param.sampl)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -p");
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
    if((fi = fopen(inputFileNameTree, "r"))) {
		TypeDistribution *d;
		TypeInfoDrawTreeGeneric info;
		TypeAdditionalDrawTreeGeneric add;
		TypeDataDrawDensity dataD;		
		double *timeSave, *timeSaveX;
		int n, *todo;
        tree = readTree(fi);
        fclose(fi);
        if(tree == NULL)
            return 1;
        toBinary(tree);
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
		fillDistributionAll(def, tree, &param, todo, d);
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
if((fo = fopen("out.csv", "w"))) {
fprintDistribution(fo, d[1]);
fclose(fo);
}
		if(format != '0') {
			char outputFileNameG[SIZE_BUFFER_CHAR], *outputFileName;
			outputFileName = outputPrefix;
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
			if(todo[n])
				fprintf(stdout," '%s%d.csv' using 1:2 with lines title 'node %d',", outputPrefix, n, n);
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
