#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <signal.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "Utils.h"
#include "Tree.h"
#include "SimulTree.h"
#include "Model.h"
#include "Uncertainty.h"
#include "NLOpt.h"

#define END_DOUBLE NO_TIME
#define STRING_SIZE 300
#define HELP_MESSAGE "NAME\n	shif - ROC plots of diversification shift tests\n	\nSYNOPSIS\n	shif [OPTIONS] <output Ident>\n\nDESCRIPTION\n	Compute ROC plots for 4 tests of diversification shifts from Yule tree simulations. Results are returned as 4 table files with names \"ROC_?_<output Ident>.csv\" where columns 1 and 2 is the true and the false positive rates respectively. \n\n	Options are\n	-o <NLopt option file>\n		read the numerical optimizers parameters from  <NLopt option file>\n	-a <birth rate>\n		set the general birth rate\n	-b <birth rate>\n		set the shifted birth rate\n	-s <time>\n		set the shift time\n	-t <time>\n		set the total diversification time (it starts from time 0)\n	-i <number>\n		set the number of simulations\n	-x <number>\n		set the number of simultaneous threads in parallel\n	-m <number>\n		set the minimum number of tips required\n	-M <number>\n		set the maximum number of nodes allowed\n	-h\n		display help\n"
//./asse -a 0.5 0 1 -b 0.75 0 1 -s 5 -i 100 -m 10 -M 5000 -t 10 -x 40

typedef struct THREAD_PARAMETER {
	int *number, ind, new, old;
	pthread_mutex_t *mutex_result;
	pthread_mutex_t *mutex_number;
	pthread_cond_t *cond_number;
	double **result, shiftTime;
	TypeTree *tree;
	TypeNLOptOption *nloptOption;
} TypeThreadParameter;


typedef struct MAXIMIZATION_DATA_4 {
    TypeTree *tree;
    int node;
    double deathA, deathB, samplA, samplB, shiftTime;
} TypeMaximizationData4;

typedef struct MAXIMIZATION_DATA_5 {
    TypeTree *tree;
    int node;
    double death, sampl, shiftTime;
} TypeMaximizationData5;

typedef struct MAXIMIZATION_DATA_6 {
    TypeTree *tree;
    int node;
    double deathA, samplA, birthB, deathB, samplB, shiftTime;
} TypeMaximizationData6;

#define NTEST 10

//static double toMaximizeData4(unsigned n, const double *x, double *grad, void *data);
//static double toMaximizeData5(unsigned n, const double *x, double *grad, void *data);
static double toMaximizeData6(unsigned n, const double *x, double *grad, void *data);
static double h(double x);
static void *threadCompute(void *data);
static double getTotalTimeSubtree(int n, TypeTree *tree);
static double *getROCPlot(double *list0, int size0, double *list1, int size1);

int main(int argc, char **argv) {	
	char outputName[STRING_SIZE], *outputPrefix, option[256];
	int i, niter = 5, minContemp = 5, maxSizeTree = 1000, minSizeTree = 100, nT, cont, ind, maxT = 8;
	TypeModelParam paramA = {.birth=0.75, .death=0.1, .sampl=0.5}, paramB = {.birth=0.75, .death = 0.1, .sampl=1.};
	double maxTime = 10., shiftTime = NO_TIME;
	gsl_rng *rg = gsl_rng_alloc(gsl_rng_random_glibc2);	
	TypeNLOptOption nloptOption;
	pthread_mutex_t mutexR = PTHREAD_MUTEX_INITIALIZER, mutexN = PTHREAD_MUTEX_INITIALIZER;
	pthread_cond_t condN = PTHREAD_COND_INITIALIZER;
	time_t t0, t1;
	double **result, *table;
	FILE *fo;
		
    nloptOption.trials = 5;
    nloptOption.tolOptim = 0.00001;
    nloptOption.maxIter = 500;
	paramA.death = 0.;
	paramA.sampl = 1.;
	paramB.death = 0.;
	paramB.sampl = 1.;

	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['o']) {
			FILE *fopt;
			option['o'] = 0;
			if((i+1)<argc) {
				if((fopt = fopen(argv[++i], "r"))) {
					fscanNLoptOptionTag(fopt, &nloptOption);
					fclose(fopt);
				} else {
					fprintf(stderr, "Can't open option file %s\n", argv[++i]);
					exit(EXIT_FAILURE);
				}
			} else {
				fprintf(stderr, "File name missing after option -o\n");
				exit(EXIT_FAILURE);
			}
		}
		if(option['a']) {
			option['a'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &(paramA.birth)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "1 value is expected after -p");
		}
		if(option['b']) {
			option['b'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &(paramB.birth)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "1 value is expected after -p");
		}
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &shiftTime) == 1)
				i++;
		}
		if(option['i']) {
			option['i'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &niter) == 1)
				i++;
		}
		if(option['x']) {
			option['x'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &maxT) == 1)
				i++;
		}
		if(option['m']) {
			option['m'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &minContemp) == 1)
				i++;
		}
		if(option['M']) {
			option['M'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &maxSizeTree) == 1)
				i++;
		}
		if(option['t']) {
			option['t'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &maxTime) == 1)
				i++;
		}
		if(option['h']) {
			printf("%s\n", HELP_MESSAGE);
			exit(EXIT_SUCCESS);
		}
	}
	if(i<argc) {
		if((outputPrefix = strrchr(argv[i], '.')) != NULL)
			outputPrefix[0] = '\0';
		if((outputPrefix=strrchr(argv[i], '/')) == NULL)
			outputPrefix = argv[i];
		else
			outputPrefix++;
	} else
		outputPrefix = "";
	time(&t0);
	nT = 0;
	cont = 1;
	ind = 0;
	result = (double**) malloc(NTEST*sizeof(double*));
	for(i=0; i<NTEST; i++)
		result[i] = (double*) malloc(niter*sizeof(double));
	while(cont) {
		int new, old;
		pthread_mutex_lock(&mutexN);
		while(ind < niter && nT < maxT) {
			pthread_t thread;	
			int ret = 0;
			TypeThreadParameter *param;
			TypeTree *tree = NULL;
//printf("paramA %.2lf %.2lf %.2lf\nparamB %.2lf %.2lf %.2lf\n", paramA.birth, paramA.death, paramA.sampl, paramB.birth, paramB.death, paramB.sampl);
			do {
				if(tree != NULL)
					freeTree(tree);
				tree = simulTree(rg, 0., maxTime, &paramA);
				addTree(rg, shiftTime, &old, &new, tree, &paramB);
			} while(tree == NULL || new == tree->root || tree->size<=minSizeTree || tree->size>=maxSizeTree);
			param = (TypeThreadParameter*) malloc(sizeof(TypeThreadParameter));
			param->tree = tree;
			param->old = old;
			param->new = new;
			param->ind = ind++;
			param->result = result;
			param->shiftTime = shiftTime;
			param->mutex_result = &mutexR;
			param->mutex_number = &mutexN;
			param->cond_number = &condN;
			param->number = &nT;
			param->nloptOption = &nloptOption;
fprintf(stdout, "simulation %d/%d\r", ind+1, niter); fflush(stdout);
			if((ret = pthread_create(&thread, NULL, threadCompute, (void*) param)) == 0) {
				int err;
				if((err = pthread_detach(thread)) == 0) {
					nT++;
				} else {
					fprintf (stderr, "Error %d while detaching thread: %s\n", err, (char*) strerror(err));
//					pthread_kill(thread, 0);
				}
			} else
				fprintf (stderr, "Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
		}
		cont = (nT > 0);
		if(cont)
			pthread_cond_wait(&condN, &mutexN);
		pthread_mutex_unlock(&mutexN);
	}
	gsl_rng_free(rg);
	for(i=0; i<NTEST; i++)
		qsort((void*)result[i], niter, sizeof(double), compareDouble);
	table = getROCPlot(result[5], niter, result[0], niter);
	sprintf(outputName, "ROC_N_%s.csv", outputPrefix);
	if((fo = fopen(outputName, "w"))) {
		for(i=0; table[i] != END_DOUBLE; i+=3)
			fprintf(fo, "%le\t%le\t%le\n", table[i], table[i+1], table[i+2]);
		free((void*)table);
		fclose(fo);
	}
	sprintf(outputName, "ROC_X_%s.csv", outputPrefix);
	if((fo = fopen(outputName, "w"))) {
		table = getROCPlot(result[6], niter, result[1], niter);
		for(i=0; table[i] != END_DOUBLE; i+=3)
			fprintf(fo, "%le\t%le\t%le\n", table[i], table[i+1], table[i+2]);
		free((void*)table);
		fclose(fo);
	}
	sprintf(outputName, "ROC_O_%s.csv", outputPrefix);
	if((fo = fopen(outputName, "w"))) {
		table = getROCPlot(result[7], niter, result[2], niter);
		for(i=0; table[i] != END_DOUBLE; i+=3)
			fprintf(fo, "%le\t%le\t%le\n", table[i], table[i+1], table[i+2]);
		free((void*)table);
		fclose(fo);
	}
	sprintf(outputName, "ROC_A_%s.csv", outputPrefix);
	if((fo = fopen(outputName, "w"))) {
		table = getROCPlot(result[8], niter, result[3], niter);
		for(i=0; table[i] != END_DOUBLE; i+=3)
			fprintf(fo, "%le\t%le\t%le\n", table[i], table[i+1], table[i+2]);
		free((void*)table);
		fclose(fo);
	}
	sprintf(outputName, "ROC_P_%s.csv", outputPrefix);
	if((fo = fopen(outputName, "w"))) {
		table = getROCPlot(result[9], niter, result[4], niter);
		for(i=0; table[i] != END_DOUBLE; i+=3)
			fprintf(fo, "%le\t%le\t%le\n", table[i], table[i+1], table[i+2]);
		free((void*)table);
		fclose(fo);
	}
	time(&t1);
	fprintf(stderr, "Computation time: %.0lf hours %.0lf mins %.0lf secs\n", floor(difftime(t1, t0)/(3600.)), floor(difftime(t1, t0)/(60.))-60.*floor(difftime(t1, t0)/(3600.)), difftime(t1, t0)-60.*floor(difftime(t1, t0)/(60.)));
	return 0;
}


//double toMaximizeData4(unsigned n, const double *x, double *grad, void *data) {
	//TypePiecewiseModelParam modelA;
	//double startA[2];
	//TypeModelParam paramA;
	//paramA.birth = x[0];
	//paramA.death = ((TypeMaximizationData4*)data)->deathA;
	//paramA.sampl = ((TypeMaximizationData4*)data)->samplA;
	//modelA.size = 1;
	//startA[0] = ((TypeMaximizationData4*)data)->tree->minTime;
	//startA[1] = ((TypeMaximizationData4*)data)->tree->maxTime;
	//modelA.startTime = startA;
	//modelA.param = &paramA;
	//TypePiecewiseModelParam modelB;
	//double startB[2];
	//TypeModelParam paramB;
	//paramB.birth = x[1];
	//paramB.death = ((TypeMaximizationData4*)data)->deathB;
	//paramB.sampl = ((TypeMaximizationData4*)data)->samplB;
	//modelB.size = 1;
	//startB[0] = ((TypeMaximizationData4*)data)->shiftTime;
	//startB[1] = ((TypeMaximizationData4*)data)->tree->maxTime;
	//modelB.startTime = startB;
	//modelB.param = &paramB;
	//return getLogProbabilityDensityShiftX(((TypeMaximizationData4*)data)->node, ((TypeMaximizationData4*)data)->shiftTime, ((TypeMaximizationData4*)data)->tree, &modelA, &modelB);
//}

//double toMaximizeData5(unsigned n, const double *x, double *grad, void *data) {
	//TypePiecewiseModelParam modelA, modelB;
	//double startA[2], startB[2];
	//TypeModelParam param;
	//param.birth = x[0];
	//param.death = ((TypeMaximizationData5*)data)->death;
	//param.sampl = ((TypeMaximizationData5*)data)->sampl;
	//modelA.size = 1;
	//startA[0] = ((TypeMaximizationData4*)data)->tree->minTime;
	//startA[1] = ((TypeMaximizationData4*)data)->tree->maxTime;
	//modelA.startTime = startA;
	//modelA.param = &param;
	//modelB.size = 1;
	//startB[0] = ((TypeMaximizationData4*)data)->shiftTime;
	//startB[1] = ((TypeMaximizationData4*)data)->tree->maxTime;
	//modelB.startTime = startB;
	//modelB.param = &param;
	//return getLogProbabilityDensityShiftX(((TypeMaximizationData5*)data)->node, ((TypeMaximizationData5*)data)->shiftTime, ((TypeMaximizationData5*)data)->tree, &modelA, &modelB);
//}

double toMaximizeData6(unsigned n, const double *x, double *grad, void *data) {
	TypePiecewiseModelParam modelA;
	double startA[2];
    TypeModelParam paramA;
    paramA.birth = x[0];
    paramA.death = ((TypeMaximizationData6*)data)->deathA;
    paramA.sampl = ((TypeMaximizationData6*)data)->samplA;
	modelA.size = 1;
	startA[0] = ((TypeMaximizationData4*)data)->tree->minTime;
	startA[1] = ((TypeMaximizationData4*)data)->tree->maxTime;
	modelA.startTime = startA;
	modelA.param = &paramA;
	TypePiecewiseModelParam modelB;
	double startB[2];
    TypeModelParam paramB;
    paramB.birth = ((TypeMaximizationData6*)data)->birthB;
    paramB.death = ((TypeMaximizationData6*)data)->deathB;
    paramB.sampl = ((TypeMaximizationData6*)data)->samplB;
	modelB.size = 1;
	startB[0] = ((TypeMaximizationData4*)data)->shiftTime;
	startB[1] = ((TypeMaximizationData4*)data)->tree->maxTime;
	modelB.startTime = startB;
	modelB.param = &paramB;
    return getLogProbabilityDensityShiftX(((TypeMaximizationData6*)data)->node, ((TypeMaximizationData6*)data)->shiftTime, ((TypeMaximizationData6*)data)->tree, &modelA, &modelB);
}

double h(double x) {
	if(x>0.)
		return x*log(x);
	return 0.;
}

void *threadCompute(void *data) {
	//TypeMaximizationData4 data4 = {.samplA=1., .deathA=0., .samplB=1., .deathB=0., .shiftTime=((TypeThreadParameter*)data)->shiftTime, .node=((TypeThreadParameter*)data)->new, .tree=((TypeThreadParameter*)data)->tree};
	//TypeMaximizationData5 data5 = {.sampl=1., .death=0., .shiftTime=((TypeThreadParameter*)data)->shiftTime, .node=((TypeThreadParameter*)data)->new, .tree=((TypeThreadParameter*)data)->tree};
	TypeMaximizationData6 data6 = {.samplA=1., .deathA=0., .samplB=1., .deathB=0., .shiftTime=((TypeThreadParameter*)data)->shiftTime, .node=((TypeThreadParameter*)data)->new, .tree=((TypeThreadParameter*)data)->tree};	int *parent, *nleaves, N1, N2, indRes, m;
	TypeModelParam paramA = {.death=0., .sampl=1.}, paramB = {.death=0., .sampl=1.};
	TypePiecewiseModelParam modelA, modelB;
	double startA[2], startB[2];
	double lowerBounds[2]= {0.,0.}, upperBounds[2]= {5.,5.}, maxN, t1, t2, lm1, lm2, lm, T;
	indRes = 0;
	T = (((TypeThreadParameter*)data)->tree->maxTime-((TypeThreadParameter*)data)->shiftTime),
	parent = getParent(((TypeThreadParameter*)data)->tree);
	nleaves = getNLeaves(((TypeThreadParameter*)data)->tree);
	if(((TypeThreadParameter*)data)->tree->node[parent[((TypeThreadParameter*)data)->new]].child == ((TypeThreadParameter*)data)->new)
		m = ((TypeThreadParameter*)data)->tree->node[((TypeThreadParameter*)data)->new].sibling;
	else
		m = ((TypeThreadParameter*)data)->tree->node[parent[((TypeThreadParameter*)data)->new]].child;
	if(nleaves[((TypeThreadParameter*)data)->new] > nleaves[m]) {
		N1 = nleaves[((TypeThreadParameter*)data)->new];
		N2 = nleaves[m];
	} else {
		N1 = nleaves[m];
		N2 = nleaves[((TypeThreadParameter*)data)->new];
	}
	lm1 = log((double)N1)/T;
	lm2 = log((double)N2)/T;
	lm = (log((double) (N1+N2))-log(2.))/T;
	paramB.birth = lm1;
	paramA.birth = lm;
	modelA.size = 1;
	startA[0] = ((TypeThreadParameter*)data)->tree->minTime;
	startA[1] = ((TypeThreadParameter*)data)->tree->maxTime;
	modelA.startTime = startA;
	modelA.param = &paramA;
	modelB.size = 1;
	startB[0] = ((TypeThreadParameter*)data)->shiftTime;
	startB[1] = ((TypeThreadParameter*)data)->tree->maxTime;
	modelB.startTime = startB;
	modelB.param = &paramB;
	data6.birthB = lm1;
	lm = (lm1+lm2)/2.;
	maximizeNLOPT(1, &(paramA.birth), lowerBounds, upperBounds, &maxN, toMaximizeData6, (void*) &data6, ((TypeThreadParameter*)data)->nloptOption);
	pthread_mutex_lock(((TypeThreadParameter*)data)->mutex_result);
		((TypeThreadParameter*)data)->result[indRes++][((TypeThreadParameter*)data)->ind] = 2.*(getLogProbabilityDensityShiftX(((TypeThreadParameter*)data)->new, ((TypeThreadParameter*)data)->shiftTime, ((TypeThreadParameter*)data)->tree, &modelA, &modelB)-getLogProbabilityDensityShiftX(((TypeThreadParameter*)data)->new, ((TypeThreadParameter*)data)->shiftTime, ((TypeThreadParameter*)data)->tree, &modelA, &modelA));
		((TypeThreadParameter*)data)->result[indRes++][((TypeThreadParameter*)data)->ind] = -(2.*lm-lm1-lm2)*T+(N1>1?((double)N1-1.)*log(1-exp(-lm1*T)):0.)+(N2>1?((double)N2-1)*log(1-exp(-lm2*T)):0.)-((N1+N2)>2?((double)N1+N2-2.)*log(1-exp(-lm*T)):0.);
		((TypeThreadParameter*)data)->result[indRes++][((TypeThreadParameter*)data)->ind] = -2*((double)N2)/((double)N1+N2-1.);
		((TypeThreadParameter*)data)->result[indRes++][((TypeThreadParameter*)data)->ind] = -1.629*(h((double)N1-1.)-h((double)N1)+h((double)N2-1.)-h((double)N2)-h(2.)-h((double)N1+N2-2.)+h((double)N1+N2));
	pthread_mutex_unlock(((TypeThreadParameter*)data)->mutex_result);
	t1 = getTotalTimeSubtree(((TypeThreadParameter*)data)->new, ((TypeThreadParameter*)data)->tree);
	t2 = getTotalTimeSubtree(m, ((TypeThreadParameter*)data)->tree);
	//l1 = ((double) N1-1.)/t1;
	//l2 = ((double) N2-1.)/t2;
	//l = ((double) N1+N2-2.)/(t1+t2);
	pthread_mutex_lock(((TypeThreadParameter*)data)->mutex_result);
//		((TypeThreadParameter*)data)->result[indRes++][((TypeThreadParameter*)data)->ind] = -2.*pow(l1/l, (double)N1-1.)*pow(l2/l, (double)N2-1.);
		((TypeThreadParameter*)data)->result[indRes++][((TypeThreadParameter*)data)->ind] = -2.*(h((double) N1-1.)-((double) N1-1.)*(log(t1)-1.)+h((double) N2-1.)-((double) N2-1.)*(log(t2)-1.)-h((double) N1+N2-2.)+((double) N1+N2-2.)*(t1+t2));
	pthread_mutex_unlock(((TypeThreadParameter*)data)->mutex_result);	
	switchNodes(((TypeThreadParameter*)data)->new, ((TypeThreadParameter*)data)->old, parent[((TypeThreadParameter*)data)->new], ((TypeThreadParameter*)data)->tree);
	free((void*)parent);
	free((void*)nleaves);
	parent = getParent(((TypeThreadParameter*)data)->tree);
	nleaves = getNLeaves(((TypeThreadParameter*)data)->tree);
	if(((TypeThreadParameter*)data)->tree->node[parent[((TypeThreadParameter*)data)->old]].child == ((TypeThreadParameter*)data)->old)
		m = ((TypeThreadParameter*)data)->tree->node[((TypeThreadParameter*)data)->old].sibling;
	else
		m = ((TypeThreadParameter*)data)->tree->node[parent[((TypeThreadParameter*)data)->old]].child;
	if(nleaves[((TypeThreadParameter*)data)->old] > nleaves[m]) {
		N1 = nleaves[((TypeThreadParameter*)data)->old];
		N2 = nleaves[m];
	} else {
		N1 = nleaves[m];
		N2 = nleaves[((TypeThreadParameter*)data)->old];
	}
	lm1 = log((double)N1)/(((TypeThreadParameter*)data)->tree->maxTime-((TypeThreadParameter*)data)->shiftTime);
	lm2 = log((double)N2)/(((TypeThreadParameter*)data)->tree->maxTime-((TypeThreadParameter*)data)->tree->time[parent[((TypeThreadParameter*)data)->old]]);
	lm = log((double)nleaves[((TypeThreadParameter*)data)->tree->root])/(((TypeThreadParameter*)data)->tree->maxTime-((TypeThreadParameter*)data)->tree->minTime);
	paramB.birth = lm1;
	paramA.birth = lm;
	lm = (lm1+lm2)/2.;
	maximizeNLOPT(1, &(paramA.birth), lowerBounds, upperBounds, &maxN, toMaximizeData6, (void*) &data6, ((TypeThreadParameter*)data)->nloptOption);
	pthread_mutex_lock(((TypeThreadParameter*)data)->mutex_result);
		((TypeThreadParameter*)data)->result[indRes++][((TypeThreadParameter*)data)->ind] = 2.*(getLogProbabilityDensityShiftX(((TypeThreadParameter*)data)->old, ((TypeThreadParameter*)data)->shiftTime, ((TypeThreadParameter*)data)->tree, &modelA, &modelB)-getLogProbabilityDensityShiftX(((TypeThreadParameter*)data)->old, ((TypeThreadParameter*)data)->shiftTime, ((TypeThreadParameter*)data)->tree, &modelA, &modelA));
		((TypeThreadParameter*)data)->result[indRes++][((TypeThreadParameter*)data)->ind] = -(2.*lm-lm1-lm2)*T+(N1>1?((double)N1-1.)*log(1-exp(-lm1*T)):0.)+(N2>1?((double)N2-1)*log(1-exp(-lm2*T)):0.)-((N1+N2)>2?((double)N1+N2-2.)*log(1-exp(-lm*T)):0.);
		((TypeThreadParameter*)data)->result[indRes++][((TypeThreadParameter*)data)->ind] = -2*((double)N2)/((double)N1+N2-1.);
		((TypeThreadParameter*)data)->result[indRes++][((TypeThreadParameter*)data)->ind] = -1.629*(h((double)N1-1.)-h((double)N1)+h((double)N2-1.)-h((double)N2)-h(2.)-h((double)N1+N2-2.)+h((double)N1+N2));
	pthread_mutex_unlock(((TypeThreadParameter*)data)->mutex_result);
	t1 = getTotalTimeSubtree(((TypeThreadParameter*)data)->old, ((TypeThreadParameter*)data)->tree);
	t2 = getTotalTimeSubtree(m, ((TypeThreadParameter*)data)->tree);
	//l1 = ((double) N1-1.)/t1;
	//l2 = ((double) N2-1.)/t2;
	//l = ((double) N1+N2-2.)/(t1+t2);
	pthread_mutex_lock(((TypeThreadParameter*)data)->mutex_result);
//		((TypeThreadParameter*)data)->result[indRes++][((TypeThreadParameter*)data)->ind] = -2.*pow(l1/l, (double)N1-1.)*pow(l2/l, (double)N2-1.);
		((TypeThreadParameter*)data)->result[indRes++][((TypeThreadParameter*)data)->ind] = -2.*(h((double) N1-1.)-((double) N1-1.)*(log(t1)-1.)+h((double) N2-1.)-((double) N2-1.)*(log(t2)-1.)-h((double) N1+N2-2.)+((double) N1+N2-2.)*(t1+t2));
	pthread_mutex_unlock(((TypeThreadParameter*)data)->mutex_result);	
	free((void*)parent);
	free((void*)nleaves);
	freeTree(((TypeThreadParameter*)data)->tree);
	pthread_mutex_lock(((TypeThreadParameter*)data)->mutex_number);
		(*((TypeThreadParameter*)data)->number)--;
		pthread_cond_signal(((TypeThreadParameter*)data)->cond_number);
	pthread_mutex_unlock(((TypeThreadParameter*)data)->mutex_number);
	free(data);
	return NULL;	
}

double getTotalTimeSubtree(int n, TypeTree *tree) {
	int c;
	double res = tree->time[n];
	for(c=tree->node[n].child; c!=NOSUCH; c=tree->node[c].sibling)
		res += getTotalTimeSubtree(c, tree);
	return res;
}

double *getROCPlot(double *list0, int size0, double *list1, int size1) {
    double *res, cur;
    int ind0, ind1, ind;
    res = (double*) malloc((3*(size0+size1+1)+1)*sizeof(double));	
    cur = utils_MIN(list0[0],list1[0]);
    for(ind=0, ind0=0, ind1=0; ind0<size0 || ind1<size1; ind+=3) {
//printf("ind %d ind0 %d ind1 %d %.3le %.3le\n", ind, ind0, ind1, list0[ind0], list1[ind1]);
if(isnan(list0[ind0]) || isnan(list1[ind1])) {
	fprintf(stderr, "Error: Nan in an array\n");
	exit(1);
}
        res[ind] = ((double)size1-ind1)/((double)size1);
        res[ind+1] = ((double)size0-ind0)/((double)size0);
        res[ind+2] = cur;
        while(ind0<size0 && list0[ind0]<=cur)
            ind0++;
        while(ind1<size1 && list1[ind1]<=cur)
            ind1++;
 		if(ind0<size0) {
			if(ind1<size1)
				cur = utils_MIN(list0[ind0],list1[ind1]);
			else
				cur = list0[ind0];
		} else
			if(ind1<size1)
				cur = list1[ind1];
    }
    res[ind] = END_DOUBLE;
    res = (double*) realloc((void*)res, (3*ind+1)*sizeof(double));
    return res;
}
