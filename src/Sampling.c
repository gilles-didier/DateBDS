#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "Utils.h"
#include "Uncertainty.h"
#include "Sampling.h"

static double inverseFuncDistribution(double x, double inf, double sup, double eps, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeParameters *param);
static void sampleRootRec(double eps, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeParameters *param, gsl_rng *rg);

double inverseFuncDistribution(double x, double inf, double sup, double eps, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeParameters *param) {
	double a=inf, b=sup;
	while((b-a) > eps) {
		double c = (b+a)/2., y = getDistributionRoot(c, tree, lr, param);
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

void sampleRootRec(double eps, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeParameters *param, gsl_rng *rg) {
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

void sampleTimes(double eps, TypeTree *tree, TypeModelParam *param, gsl_rng *rg) {
	if(tree != NULL && tree->size > 0) {
		TypeUncertaintyLikeParameters pu = getUncertaintyLikeParameters(param);
		TypeUncertaintyLogRankData *lr = (TypeUncertaintyLogRankData*) malloc(tree->size*sizeof(TypeUncertaintyLogRankData));
		fillUncertaintyLogRankData(tree->root, tree, lr);
		sampleRootRec(eps, tree, lr, &pu, rg);
		free((void*)lr);
	}
}
