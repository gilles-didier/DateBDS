#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "Utils.h"
#include "Uncertainty.h"
#include "Sampling.h"

static double inverseFuncDistribution(double x, double inf, double sup, double eps, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeCoeff *param);
static void sampleRootRec(double eps, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeCoeff *param, gsl_rng *rg);

double inverseFuncDistribution(double x, double inf, double sup, double eps, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeCoeff *param) {
	double a=inf, b=sup;
	while((b-a) > eps) {
		double c = (b+a)/2., y = getDistributionRoot(getPieceIndex(c, &(param->model)), c, tree, lr, param);
		if(y>x)
			b = c;
		else {
			if(y<x)
				a = c;
			else
				return c;
		}
	}
	return (a+b)/2.;	
}

void test(double eps, TypeTree *tree, TypePiecewiseModelParam *param) {
	if(tree != NULL && tree->size > 0) {
		TypeUncertaintyLikeCoeff pu = getUncertaintyLikeCoeff(param);
		TypeUncertaintyLogRankData *lr = (TypeUncertaintyLogRankData*) malloc(tree->size*sizeof(TypeUncertaintyLogRankData));
		fillUncertaintyLogRankData(tree->root, tree, lr);
		double x;
		for(x=0.01; x<1; x+=0.01)
			printf("inv node %d x %lf time %lf\n", tree->root, x, inverseFuncDistribution(x, tree->minTime, tree->maxTime, eps, tree, lr, &pu));
		printf("\n\n");
		for(x=tree->minTime; x<tree->maxTime; x+=0.01)
			printf("dist node %d x %lf %d time %lf\n", tree->root, x, getPieceIndex(x, &(pu.model)), getDistributionRoot(getPieceIndex(x, &(pu.model)), x, tree, lr, &pu));
		free((void*)lr);
	}
}

void sampleRootRec(double eps, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeCoeff *param, gsl_rng *rg) {
	if(tree->node[tree->root].child != NOSUCH) {
		int c;
		double x = gsl_rng_uniform_pos(rg);
		tree->time[tree->root] = inverseFuncDistribution(x, tree->minTime, tree->maxTime, eps, tree, lr, param);
		for(c=tree->node[tree->root].child; c!=NOSUCH; c=tree->node[c].sibling) {
			int rootSave;
			double minTimeSave;
			rootSave = tree->root;
			minTimeSave = tree->minTime;
			tree->root = c;
			tree->minTime = tree->time[rootSave];
			sampleRootRec(eps, tree, lr, param, rg);
			tree->root = rootSave;
			tree->minTime = minTimeSave;
		}
	}
}

void sampleTimes(double eps, TypeTree *tree, TypePiecewiseModelParam *param, gsl_rng *rg) {
	if(tree != NULL && tree->size > 0) {
		TypeUncertaintyLikeCoeff pu = getUncertaintyLikeCoeff(param);
		TypeUncertaintyLogRankData *lr = (TypeUncertaintyLogRankData*) malloc(tree->size*sizeof(TypeUncertaintyLogRankData));
		fillUncertaintyLogRankData(tree->root, tree, lr);
		sampleRootRec(eps, tree, lr, &pu, rg);
		free((void*)lr);
	}
}
