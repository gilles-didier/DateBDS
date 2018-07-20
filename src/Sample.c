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

#define HELPMESSAGE "NAME\n\tsamp - divergence times sampling under the birth-death-sampling model\n\t\nSYNOPSIS\n\tsamp [OPTIONS] <input Tree File> <output Ident>\n\nDESCRIPTION\n\tSample the divergence time distribution of the phylogenetic tree contained in <input Tree File> (which must be in Newick format) with regards  to the parameters provided as options and output a dated tree in Newick format (with extension '.phi').\n\n\tOptions are\n\t-o <origin time> <end time>\n\t\tset the origin and end time of the diversification process resulting to the phylogenetic tree. \n\t-p <speciation rate> <extinction rate> <sampling probability>\n\t\tset the parameters of the birth-death-sampling model. \n\t-z <input Tree File>\n\t\toutput the tree in 'text' format in the console and save it in pdf format as 'tree.pdf' with the internal idents of the nodes (for debug purposes) next exit\n\t-e <value>\n\t\tset the precision used for sampling.\n\t-h\n\t\tdisplay help\n"

int main(int argc, char **argv) {	
	char inputFileNameTree[STRING_SIZE], option[256];
	FILE *fi, *fo;
	int i, j;
	double eps = 0.0001,  contemp = 0.,  origin = -100.;
	TypeModelParam param = {.birth=0.3, .death = 0.1, .sampl=1.};
	TypeTree *tree;
	gsl_rng *rg = gsl_rng_alloc(gsl_rng_random_glibc2);
	
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
				exitProg(ErrorArgument, "a value is expected after -o");
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
		sampleTimes(eps, tree, &param, rg);
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
