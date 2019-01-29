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
#include "Sampling.h"

#include "Distribution.h"


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
#define MAX_PIECES 20

#define HELPMESSAGE "NAME\n	samp - divergence times sampling under the birth-death-sampling model\n	\nSYNOPSIS\n	samp [OPTIONS] <input Tree File> <output Ident>\n\nDESCRIPTION\n	Sample the divergence time distribution of the phylogenetic tree contained in <input Tree File> (which must be in Newick format) with regards  to the parameters provided as options and output a dated tree in Newick format (with extension \".phy\"). Result files are saved in the directory of the input file.\n\n	Options are\n	-o <origin time> <end time>\n		set the origin and end time of the diversification process resulting to the phylogenetic tree. \n	-p <start time> <speciation rate> <extinction rate> <sampling probability>\n		set the parameters of the piece starting at <starting time> of a piecewise-constant-birth-death-sampling model. \n	-z <input Tree File>\n		output the tree in 'text' format in the console and save it in pdf format as 'tree.pdf' with the internal idents of the nodes (for debug purposes) next exit\n	-e <value>\n		set the precision used for sampling.\n	-h\n		display help\n"

//./dist -p 0.1 0.02 0.5 TEST4.txt 

int main(int argc, char **argv) {	
	char inputFileNameTree[STRING_SIZE], option[256];
	FILE *fi, *fo;
	int i, j;
	double eps = 0.0001,  contemp = 0.,  origin = -100., *startTmp;
	TypePiecewiseModelParam piece;
	TypeModelParam *paramTmp;
	TypeTree *tree;
	gsl_rng *rg = gsl_rng_alloc(gsl_rng_random_glibc2);
	size_t *pieceIndex;
	
	piece.size = 0;
	piece.startTime = (double*) malloc(MAX_PIECES*sizeof(double));
	piece.param = (TypeModelParam*) malloc(MAX_PIECES*sizeof(TypeModelParam));

	for(i=0; i<256; i++)
		option[i] = 0;
	   
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['z']) {
			option['z'] = 0;
			if((i+1)<argc && (fi = fopen(argv[++i], "r"))) {
				TypeTree *tree;
				tree = readTree(fi);
				printTreeDebug(stdout, tree->root, tree, tree->name);
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
				exitProg(ErrorArgument, "4 values are expected after -p");
			if((i+1)<argc && sscanf(argv[i+1], "%le", &(piece.param[piece.size].birth)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "4 values are expected after -p");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(piece.param[piece.size].death)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "4 values are expected after -p");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(piece.param[piece.size].sampl)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "4 values are expected after -p");
			piece.size++;
		}
		if(option['e']) {
			option['e'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &(eps)) == 1)
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
	pieceIndex = qsortTable(piece.startTime, piece.size, sizeof(double), compareDouble);
	for(i=0; i<piece.size; i++)
		printf("start %d %.1lf (%.1lf %.1lf %.1lf)\n", i, piece.startTime[i], piece.param[i].birth, piece.param[i].death, piece.param[i].sampl);
	printf("\n");
	for(i=0; i<piece.size; i++)
		printf("index %d %d\n", i, (int) pieceIndex[i]);
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
		char *tmp, outputFileNameG[SIZE_BUFFER_CHAR], *outputFileName;
		int n;
		double *timeSave;
		tree = readTree(fi);
        fclose(fi);
        if(tree == NULL)
            return 1;
        toBinary(tree);
        reorderTree(tree->name, tree);
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
		sampleTimes(eps, tree, &piece, rg);
		outputFileName = inputFileNameTree;
		if((tmp = strrchr(outputFileName, '.')) != NULL)
			tmp[0] = '\0';
		sprintf(outputFileNameG, "%s_sample.phy", outputFileName);
		if((fo = fopen(outputFileNameG, "w"))) {
			fprintTreeNewick(fo, tree);
			fclose(fo);
		} else {
			fprintf(stderr, "Cannot open %s\n", outputFileNameG);
			exit(1);
		}
		free((void*)tree->time);
		tree->time = timeSave;
		freeTree(tree);
	} else {
		fprintf(stderr, "Cannot read %s\n", inputFileNameTree);
		exit(1);
	}
	return 0;
}
