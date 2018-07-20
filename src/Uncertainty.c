#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <gsl/gsl_sf_gamma.h>

#include "Utils.h"
#include "TreeExtras.h"
#include "Uncertainty.h"

#define NO_CONF 0
#define SP_CONF INT_MAX


static void fillIndexMin(int n, TypeTree *tree, int *indexMin);
static double logBinomial(unsigned int k, unsigned int n);
static double getTreeShapeLogLikelihood(double rank, unsigned int nleaf);
static double getLogProbStart(double startTime, double endTime, double maxTime, int k, TypeUncertaintyLikeParameters *param);
static double getLogProbTypeA(double startTime, double maxTime, int k,  TypeUncertaintyLikeParameters *param);
static double *fillUncertaintyInfo(int n, double val, int node, TypeTree *tree, TypeUncertaintyLogRankData *lr, int *indexMin, TypeUncertaintyLikeParameters *param);
static double getLogLikelihoodConstraint(double val, int node, TypeTree *tree, TypeUncertaintyLogRankData *lr, int *indexMin, TypeUncertaintyLikeParameters *param);

double logBinomial(unsigned int k, unsigned int n) {
    return gsl_sf_lnfact(n)-(gsl_sf_lnfact(k)+gsl_sf_lnfact(n-k));
}

void fillUncertaintyLogRankData(int node, TypeTree *tree, TypeUncertaintyLogRankData *lr) {
	if(tree->node[node].child == NOSUCH) {
		lr[node].nLeaves = 1;
		lr[node].logRank = 0.;
	} else {
		fillUncertaintyLogRankData(tree->node[node].child, tree, lr);
		fillUncertaintyLogRankData(tree->node[tree->node[node].child].sibling, tree, lr);		
		lr[node].nLeaves = lr[tree->node[node].child].nLeaves+lr[tree->node[tree->node[node].child].sibling].nLeaves;
		lr[node].logRank = lr[tree->node[node].child].logRank+lr[tree->node[tree->node[node].child].sibling].logRank+log(2)+logBinomial((unsigned int) (lr[tree->node[node].child].nLeaves-1), (unsigned int) (lr[node].nLeaves-2));
	}
}

TypeUncertaintyLogRankData getLogRankDataShape(int index, int node, int *sizeShape, TypeTree *tree) {
	TypeUncertaintyLogRankData res;
	if(tree->node[node].child == NOSUCH || index == 0) {
		res.nLeaves = 1;
		res.logRank = 0.;
	} else {
		TypeUncertaintyLogRankData d0, d1;
		d0 = getLogRankDataShape((index-1)%sizeShape[tree->node[node].child], tree->node[node].child, sizeShape, tree);
		d1 = getLogRankDataShape((index-1)/sizeShape[tree->node[node].child], tree->node[tree->node[node].child].sibling, sizeShape, tree);
		res.nLeaves = d0.nLeaves+d1.nLeaves;
		res.logRank = d0.logRank+d1.logRank+log(2)+logBinomial((unsigned int) (d0.nLeaves-1), (unsigned int) (res.nLeaves-2));
	}
	return res;
}

double *fillUncertaintyInfo(int n, double val, int node, TypeTree *tree, TypeUncertaintyLogRankData *lr, int *indexMin, TypeUncertaintyLikeParameters *param) {
	double *res;
	if(tree->node[n].child == NOSUCH) {
		res = (double*) malloc(sizeof(double));
		res[0] = getLogProbTypeA(val, tree->maxTime, 1,  param);
	} else {
		int i, j;
		double *logLike0, *logLike1, *like, *offset;
		logLike0 = fillUncertaintyInfo(tree->node[n].child, val, node, tree, lr, indexMin, param);
		logLike1 = fillUncertaintyInfo(tree->node[tree->node[n].child].sibling, val, node, tree, lr, indexMin, param);		
		res = (double*) malloc(lr[n].nLeaves*sizeof(double));
		like = (double*) malloc(lr[n].nLeaves*sizeof(double));
		offset = (double*) malloc(lr[n].nLeaves*sizeof(double));
		for(i=0; i<lr[n].nLeaves; i++) {
			offset[i] = NEG_INFTY;
			like[i] = 0.;
		}
		if(tree->time[indexMin[n]] == tree->maxTime)
			res[0] = getLogProbTypeA(val, tree->maxTime, lr[n].nLeaves,  param)+getTreeShapeLogLikelihood(lr[n].logRank, lr[n].nLeaves);
		else
			res[0] = NEG_INFTY;
		for(i=0; i<lr[tree->node[n].child].nLeaves; i++) 
			if(logLike0[i] != NEG_INFTY) {
				for(j=0; j<lr[tree->node[tree->node[n].child].sibling].nLeaves; j++)
					if(logLike1[j] != NEG_INFTY) {
						double logLike = logLike0[i]+logLike1[j]+log(2.)+logBinomial((unsigned int) i, (unsigned int) (i+j));
						if(offset[i+j+1] == NEG_INFTY) {
							offset[i+j+1] = logLike;
							like[i+j+1] = 1.;
						} else {
							if(logLike>offset[i+j+1]) { /*compute max in offset just to avoid numerical precision issues*/
								like[i+j+1] *= exp(offset[i+j+1]-logLike);
								offset[i+j+1] = logLike;
								like[i+j+1]++;
							} else
								like[i+j+1] += exp(logLike-offset[i+j+1]);
						}
					}
			}
		for(i=1; i<lr[n].nLeaves; i++)
			if(like[i]>0.)
				res[i] = log(like[i])+offset[i];
			else
				res[i] = NEG_INFTY;
		free((void*)like);
		free((void*)offset);
		free((void*)logLike0);
		free((void*)logLike1);
	}
	return res;
}

void fillDistributionAll(int def, TypeTree *tree, TypeModelParam *param, int *todo, TypeDistribution *d) {
	if(tree != NULL && tree->size > 0) {
		TypeUncertaintyLikeParameters pu = getUncertaintyLikeParameters(param);
		int *indexMin, i, n;
		double ll0, ll1, step;
		TypeUncertaintyLogRankData *lr;
		step = (tree->maxTime-tree->minTime)/((double)def-1.);
		if(tree->parent == NULL)
			setParent(tree);
		lr = (TypeUncertaintyLogRankData*) malloc(tree->size*sizeof(TypeUncertaintyLogRankData));
		fillUncertaintyLogRankData(tree->root, tree, lr);
		ll0 = getLogProbTypeA(tree->minTime, tree->maxTime, lr[tree->root].nLeaves,  &pu)+getTreeShapeLogLikelihood(lr[tree->root].logRank, lr[tree->root].nLeaves);
		indexMin = (int*) malloc(tree->size*sizeof(int));
		for(n=0; n<tree->size; n++) {
			if(todo[n]) {
printf("\rcomputing node %d/%d      ", n, tree->size-1); fflush(stdout);
				tree->time[n] = (tree->maxTime+tree->minTime)/2.;
				fillIndexMin(tree->root, tree, indexMin);
				d[n].size = def;
				d[n].item = (TypeDistributionItem*) malloc(def*sizeof(TypeDistributionItem));
				for(i=0; i<def; i++) {
						d[n].item[i].val = tree->minTime+((double)i)*step;
						ll1 = getLogLikelihoodConstraint(d[n].item[i].val, n, tree, lr, indexMin, &pu);
						d[n].item[i].dens = exp(ll1-ll0);
				}
				tree->time[n] = NO_TIME;
			} else {
				d[n].size = 0;
				d[n].item = NULL;
			}
		}
		free((void*)lr);
		free((void*)indexMin);
	}
}

TypeUncertaintyLikeParameters getUncertaintyLikeParameters(TypeModelParam *param) {
	TypeUncertaintyLikeParameters res;
	res.birth = param->birth;
	res.death = param->death;
	res.sampl = param->sampl;
	res.bmd = res.birth-res.death;
	res.bs = res.birth*res.sampl;
	res.bsd = res.birth*(1.-res.sampl)-res.death;
	return res;
}

double getLogLikelihoodConstraint(double val, int node, TypeTree *tree, TypeUncertaintyLogRankData *lr, int *indexMin, TypeUncertaintyLikeParameters *param) {
	int i;
	double like = 0., offset = NEG_INFTY, *ll = fillUncertaintyInfo(tree->root, val, node, tree, lr, indexMin, param);
	for(i=0; i<lr[tree->root].nLeaves; i++) {
		double logLike;
		if(ll[i] != NEG_INFTY)
			logLike = getLogProbStart(tree->minTime, val, tree->maxTime, i+1, param) + ll[i] - gsl_sf_lnfact(i);
		else
			logLike = NEG_INFTY;
		if(logLike != NEG_INFTY && !isinf(logLike) && !isnan(logLike)) {
			if(offset == NEG_INFTY) {
				offset = logLike;
				like = 1.;
			} else {
				if(logLike>offset) { /*compute max in offset just to avoid numerical precision issues*/
					like *= exp(offset-logLike);
					offset = logLike;
					like++;
				} else
					like += exp(logLike-offset);
			}
		}
	}
	free((void*)ll);
	if(like>0.)
		return log(like)+offset;
	else
		return NEG_INFTY;
}

double getTreeShapeLogLikelihood(double rank, unsigned int nleaf) {
    return rank-gsl_sf_lnfact(nleaf-1);
}

/*return p(k, startTime)*/
double getLogProbTypeA(double startTime, double maxTime, int k,  TypeUncertaintyLikeParameters *param) {
    double time = maxTime-startTime;
    if(k == 1)
		return
			log(param->sampl)
			+2.*(log(param->bmd)-log(param->bs+param->bsd*exp(-param->bmd*time)))
			-param->bmd*time;
	else
		return
			((double)k)*log(param->sampl)
			+2.*log(param->bmd)
			-param->bmd*time
			+((double)k-1.)*(log(param->birth)+log(1.-exp(-param->bmd*time)))
			-((double)k+1.)*log(param->bs+param->bsd*exp(-param->bmd*time));
}

/* return P_y(k, startTime, endTime)/(P^\star_o)^(k-l)*/
double getLogProbStart(double startTime, double endTime, double maxTime, int k, TypeUncertaintyLikeParameters *param) {
    if(k == 1)
		return
			-param->bmd*(endTime-startTime)
			+2.*(log(param->bs+param->bsd*exp(-param->bmd*(maxTime-endTime)))-log(param->bs+param->bsd*exp(-param->bmd*(maxTime-startTime))));
	else
		return
			-param->bmd*(endTime-startTime)
			+((double)k-1.)*(log(param->birth)+log(1.-exp(-param->bmd*(endTime-startTime)))-log(param->bmd))
			+((double)k+1)*(log(param->bs+param->bsd*exp(-param->bmd*(maxTime-endTime)))-log(param->bs+param->bsd*exp(-param->bmd*(maxTime-startTime))));
}

/*For all nodes n of tree, indexMin[n] = leaf with the smallest time in subtree n or descendent diff from n with the smallest time*/
void fillIndexMin(int n, TypeTree *tree, int *indexMin) {
	if(tree->node[n].child == NOSUCH) {
		indexMin[n] = n;
	} else {
		int c;
		for(c=tree->node[n].child; c != NOSUCH; c=tree->node[c].sibling)
			fillIndexMin(c, tree, indexMin);
		if(tree->time[n] != NO_TIME)
			indexMin[n] = n;
		else {
			indexMin[n] = indexMin[tree->node[n].child];
			for(c=tree->node[tree->node[n].child].sibling; c != NOSUCH; c=tree->node[c].sibling)
				if(tree->time[indexMin[c]]<tree->time[indexMin[n]])
					indexMin[n] = indexMin[c];
		}
	}
}

double getDistributionRoot(double val, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeParameters *param) {
//printf("%lf (%lf, %lf)\t%le\t%le\t%le\t%le\n", val, tree->minTime, tree->maxTime, getLogProbStart(tree->minTime, val, tree->maxTime, 1, param), getLogProbTypeA(val, tree->maxTime, lr[tree->root].nLeaves, param), getLogProbTypeA(tree->minTime, tree->maxTime, lr[tree->root].nLeaves, param), getLogProbStart(tree->minTime, val, tree->maxTime, 1, param) - getLogProbTypeA(val, tree->maxTime, lr[tree->root].nLeaves, param)+getLogProbTypeA(tree->minTime, tree->maxTime, lr[tree->root].nLeaves, param));
	return 1.-exp(getLogProbStart(tree->minTime, val, tree->maxTime, 1, param) 
	+getLogProbTypeA(val, tree->maxTime, lr[tree->root].nLeaves, param)-getLogProbTypeA(tree->minTime, tree->maxTime, lr[tree->root].nLeaves, param));
}
