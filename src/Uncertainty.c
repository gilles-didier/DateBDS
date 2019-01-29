#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>

#include "Utils.h"
#include "TreeExtras.h"
#include "Uncertainty.h"

#define NO_CONF 0
#define SP_CONF INT_MAX


typedef struct UNCERTAINTY_NUMBER_RANK_DATA {
	int nLeaves;
	double lnRank;
} TypeUncertaintyLogNumberRankData;

static void fillIndexMin(int n, TypeTree *tree, int *indexMin);
static double logBinomial(unsigned int k, unsigned int n);
static double getTreeShapeLogProbabilityDensity(double rank, unsigned int nleaf);
static double getLogProbStart(double startTime, double endTime, double maxTime, int k, TypeUncertaintyLikeCoeff *param);
static double getLogProbTypeA(double startTime, double maxTime, int k,  TypeUncertaintyLikeCoeff *param);
static double *fillUncertaintyInfo(int n, double val, TypeTree *tree, TypeUncertaintyLogRankData *lr, int *indexMin, TypeUncertaintyLikeCoeff *param);
static double getLogProbabilityDensityConstraint(double val, TypeTree *tree, TypeUncertaintyLogRankData *lr, int *indexMin, TypeUncertaintyLikeCoeff *param);
static double *fillUncertaintyInfoShift(int n, double val, int *shift, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeCoeff *paramA, TypeUncertaintyLikeCoeff *paramB);
static double getLogProbabilityDensityShift(double val, int *shift, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeCoeff *coeffA, TypeUncertaintyLikeCoeff *coeffB);
static double lfA(double s, double t, int k, int n, TypeUncertaintyLikeCoeff *coeff);
static double lf(double s, double t, TypeUncertaintyLikeCoeff *coeff);
static double lF(double s, double t, TypeUncertaintyLikeCoeff *coeff);
static double *fillUncertaintyInfoGeneral(int n, int icur, int imod, double *tstop, double *tinf, double *tsup, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeCoeff *coeff);
static double getLogProbabilityDensityConstraintGeneral(int n, int icur, int imod, double *tstop, double *tinf, double *tsup, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeCoeff *coeff);
static void fillInfSup(int n, double tcur, double *tinf, double *tsup, double *tmin, double *tmax, TypeTree *tree);

double getLogProbX(double s, double e, int i, int k, TypeUncertaintyLikeCoeff *coeff);
double getLogProbY(double s, double e, int i, int k, TypeUncertaintyLikeCoeff *coeff);
double getLogProbObs(double t, int i, TypeUncertaintyLikeCoeff *coeff);
double getProbObs(double t, int i, TypeUncertaintyLikeCoeff *coeff);
double getLogProbNotObs(double t, int i, TypeUncertaintyLikeCoeff *coeff);
double getProbNotObs(double t, int i, TypeUncertaintyLikeCoeff *coeff);

void fillDistributionAll(int def, TypeTree *tree, TypePiecewiseModelParam *param, int *todo, TypeDistribution *d) {
	if(tree != NULL && tree->size > 0) {
		TypeUncertaintyLikeCoeff pu = getUncertaintyLikeCoeff(param);
		int *indexMin, i, n;
		double ll0, ll1, step;
		TypeUncertaintyLogRankData *lr;
		step = (tree->maxTime-tree->minTime)/((double)def-1.);
		if(tree->parent == NULL)
			setParent(tree);
		lr = (TypeUncertaintyLogRankData*) malloc(tree->size*sizeof(TypeUncertaintyLogRankData));
		fillUncertaintyLogRankData(tree->root, tree, lr);
		ll0 = getLogProbTypeA(tree->minTime, tree->maxTime, lr[tree->root].nLeaves,  &pu)+getTreeShapeLogProbabilityDensity(lr[tree->root].logRank, lr[tree->root].nLeaves);
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
						ll1 = getLogProbabilityDensityConstraint(d[n].item[i].val, tree, lr, indexMin, &pu);
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

double getLogProbabilityDensity(TypeTree *tree, TypePiecewiseModelParam *param) {
	TypeUncertaintyLogRankData *lr;
	TypeUncertaintyLikeCoeff coeff=getUncertaintyLikeCoeff(param);
	double res;
	lr = (TypeUncertaintyLogRankData*) malloc(tree->size*sizeof(TypeUncertaintyLogRankData));
	fillUncertaintyLogRankData(tree->root, tree, lr);
	res = getLogProbTypeA(tree->minTime, tree->maxTime, lr[tree->root].nLeaves,  &coeff)+getTreeShapeLogProbabilityDensity(lr[tree->root].logRank, lr[tree->root].nLeaves);
	free((void*)lr);
	return res;
}

double getLogProbabilityDensityShiftX(int node, double val, TypeTree *tree, TypePiecewiseModelParam *paramA, TypePiecewiseModelParam *paramB) {
	TypeUncertaintyLogRankData *lr;
	TypeUncertaintyLikeCoeff coeffA=getUncertaintyLikeCoeff(paramA), coeffB=getUncertaintyLikeCoeff(paramB);
	int *shift, n;
	double res;
	shift = (int*) malloc(tree->size*sizeof(int));
	for(n=0; n<tree->size; n++)
		shift[n] = 0;
	shift[node] = 1;
	lr = (TypeUncertaintyLogRankData*) malloc(tree->size*sizeof(TypeUncertaintyLogRankData));
	fillUncertaintyLogRankData(tree->root, tree, lr);
	res = getLogProbabilityDensityShift(val, shift, tree, lr, &coeffA, &coeffB);
	free((void*)lr);
	free((void*)shift);
	return res;
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

double *fillUncertaintyInfo(int n, double val, TypeTree *tree, TypeUncertaintyLogRankData *lr, int *indexMin, TypeUncertaintyLikeCoeff *coeff) {
	double *res;
	if(tree->node[n].child == NOSUCH) {
		res = (double*) malloc(sizeof(double));
		res[0] = getLogProbTypeA(val, tree->maxTime, 1,  coeff);
	} else {
		int i, j;
		double *logLike0, *logLike1, *like, *offset;
		logLike0 = fillUncertaintyInfo(tree->node[n].child, val,tree, lr, indexMin, coeff);
		logLike1 = fillUncertaintyInfo(tree->node[tree->node[n].child].sibling, val, tree, lr, indexMin, coeff);		
		res = (double*) malloc(lr[n].nLeaves*sizeof(double));
		like = (double*) malloc(lr[n].nLeaves*sizeof(double));
		offset = (double*) malloc(lr[n].nLeaves*sizeof(double));
		for(i=0; i<lr[n].nLeaves; i++) {
			offset[i] = NEG_INFTY;
			like[i] = 0.;
		}
		if(tree->time[indexMin[n]] == tree->maxTime)
			res[0] = getLogProbTypeA(val, tree->maxTime, lr[n].nLeaves,  coeff)+getTreeShapeLogProbabilityDensity(lr[n].logRank, lr[n].nLeaves);
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

double getLogProbabilityDensityConstraint(double val, TypeTree *tree, TypeUncertaintyLogRankData *lr, int *indexMin, TypeUncertaintyLikeCoeff *coeff) {
	int i;
	double like = 0., offset = NEG_INFTY, *ll = fillUncertaintyInfo(tree->root, val, tree, lr, indexMin, coeff);
	for(i=0; i<lr[tree->root].nLeaves; i++) {
		double logLike;
		if(ll[i] != NEG_INFTY)
			logLike = getLogProbStart(tree->minTime, val, tree->maxTime, i+1, coeff) + ll[i] - gsl_sf_lnfact(i);
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


double logProbabilityDensityConstraintGeneral(double *tmin, double *tmax, TypeTree *tree, TypePiecewiseModelParam *param) {
	int i, j, n;
	double ll0, *tinf, *tsup, *tinfX, *tsupX, *tstop;
	TypeUncertaintyLogRankData *lr;
	TypeUncertaintyLikeCoeff coeff = getUncertaintyLikeCoeff(param);
	tinf = (double*) malloc(tree->size*sizeof(double));
	tsup = (double*) malloc(tree->size*sizeof(double));
	tinfX = (double*) malloc(tree->size*sizeof(double));
	tsupX = (double*) malloc(tree->size*sizeof(double));
	fillInfSup(tree->root, tree->minTime, tinfX, tsupX, tmin, tmax, tree);
	tstop = (double*) malloc((2*(tree->size+1)+param->size)*sizeof(double));
	i=0;
	tstop[i++] = tree->minTime;
	tstop[i++] = tree->maxTime;
	for(n=0; n<param->size; n++)
		tstop[i++] = param->startTime[n];
	for(n=0; n<tree->size; n++) {
		if(tinfX[n] != NO_TIME)
			tstop[i++] = tinfX[n];			
		if(tsupX[n] != NO_TIME)
			tstop[i++] = tsupX[n];
	}
	qsort(tstop, i, sizeof(double), compareDouble);
	j=1;
	for(n=1; n<i; n++)
		if(tstop[n] != tstop[j-1])
			tstop[j++] = tstop[n];
	tstop = (double*) realloc((void*) tstop, (j+2)*sizeof(double));
	if(tree->parent == NULL)
		setParent(tree);
	lr = (TypeUncertaintyLogRankData*) malloc(tree->size*sizeof(TypeUncertaintyLogRankData));
	fillUncertaintyLogRankData(tree->root, tree, lr);
	ll0 = getLogProbabilityDensityConstraintGeneral(tree->root, 0, 0, tstop, tinfX, tsupX, tree, lr, &coeff);
	free((void*)tstop);
	free((void*)tinfX);
	free((void*)tsupX);
	free((void*)tinf);
	free((void*)tsup);
	free((void*)lr);
	return ll0;
}
void fillDistributionAllGeneral(int def, double *tmin, double *tmax, TypeTree *tree, TypePiecewiseModelParam *param, int *todo, TypeDistribution *d) {
	if(tree != NULL && tree->size > 0) {
		TypeUncertaintyLikeCoeff coeff = getUncertaintyLikeCoeff(param);
		int i, j, k, n;
		double ll0, ll1, step, *tinf, *tsup, *tinfX, *tsupX, *tstop;
		TypeUncertaintyLogRankData *lr;
		tinf = (double*) malloc(tree->size*sizeof(double));
		tsup = (double*) malloc(tree->size*sizeof(double));
		tinfX = (double*) malloc(tree->size*sizeof(double));
		tsupX = (double*) malloc(tree->size*sizeof(double));
		fillInfSup(tree->root, tree->minTime, tinfX, tsupX, tmin, tmax, tree);
		tstop = (double*) malloc((2*(tree->size+1)+param->size)*sizeof(double));
		i=0;
		tstop[i++] = tree->minTime;
		tstop[i++] = tree->maxTime;
		for(n=0; n<param->size; n++)
			tstop[i++] = param->startTime[n];
		for(n=0; n<tree->size; n++) {
			if(tinfX[n] != NO_TIME)
				tstop[i++] = tinfX[n];			
			if(tsupX[n] != NO_TIME)
				tstop[i++] = tsupX[n];
		}
		qsort(tstop, i, sizeof(double), compareDouble);
		j=1;
		for(n=1; n<i; n++)
			if(tstop[n] != tstop[j-1])
				tstop[j++] = tstop[n];
		tstop = (double*) realloc((void*) tstop, (j+2)*sizeof(double));
		step = (tree->maxTime-tree->minTime)/((double)def-1.);
		if(tree->parent == NULL)
			setParent(tree);
		lr = (TypeUncertaintyLogRankData*) malloc(tree->size*sizeof(TypeUncertaintyLogRankData));
		fillUncertaintyLogRankData(tree->root, tree, lr);
		ll0 = getLogProbabilityDensityConstraintGeneral(tree->root, 0, 0, tstop, tinfX, tsupX, tree, lr, &coeff);
		for(n=0; n<tree->size; n++) {
			if(todo[n]) {
printf("\n\n\rcomputing node %d/%d      \n", n, tree->size-1); fflush(stdout);
				double maxTmp = tmax[n];
				d[n].size = def;
				d[n].item = (TypeDistributionItem*) malloc((def+1)*sizeof(TypeDistributionItem));
				d[n].item[0].val = tree->minTime;
				d[n].item[0].dens = 0.;
				for(i=1; i<def && tree->minTime+((double)i)*step<tinfX[n]; i++) {
					d[n].item[i].val = tree->minTime+((double)i)*step;
					d[n].item[i].dens = 0.;
				}		
				for(; i<def && tree->minTime+((double)i)*step<tsupX[n]; i++) {
					int v;
					d[n].item[i].val = tree->minTime+((double)i)*step;
					tmax[n] = d[n].item[i].val;
					fillInfSup(tree->root, tree->minTime, tinf, tsup, tmin, tmax, tree);
					free((void*)tstop);
					tstop = (double*) malloc((2*(tree->size+1)+param->size)*sizeof(double));
					k=0;
					tstop[k++]= tree->minTime;
					tstop[k++]= tree->maxTime;
					for(v=0; v<param->size; v++)
						tstop[k++] = param->startTime[v];
					for(v=0; v<tree->size; v++) {
						if(tmin[v] != NO_TIME)
							tstop[k++] = tmin[v];			
						if(tmax[v] != NO_TIME)
							tstop[k++] = tmax[v];
					}
					qsort(tstop, k, sizeof(double), compareDouble);
					j=1;
					for(v=1; v<k; v++)
						if(tstop[v] != tstop[j-1])
							tstop[j++] = tstop[v];
					ll1 = getLogProbabilityDensityConstraintGeneral(tree->root, 0, 0, tstop, tinf, tsup, tree, lr, &coeff);
					d[n].item[i].dens = exp(ll1-ll0);
				}
				for(; i<def; i++) {
					d[n].item[i].val = tree->minTime+((double)i)*step;
					d[n].item[i].dens = 1.;
				}					
				tmax[n] = maxTmp;
			} else {
				d[n].size = 0;
				d[n].item = NULL;
			}
		}
		free((void*)tstop);
		free((void*)tinfX);
		free((void*)tsupX);
		free((void*)tinf);
		free((void*)tsup);
		free((void*)lr);
	}
}

double *fillUncertaintyInfoGeneral(int n, int icur, int imod, double *tstop, double *tinf, double *tsup, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeCoeff *coeff) {
	double *res;
	if(tree->node[n].child == NOSUCH) {
		res = (double*) malloc(sizeof(double));
		res[0] = getLogProbabilityDensityConstraintGeneral(n, icur, imod, tstop, tinf, tsup, tree, lr, coeff);
	} else {
		int i, j;
		double *logLike0, *logLike1, *like, *offset;
		logLike0 = fillUncertaintyInfoGeneral(tree->node[n].child, icur, imod, tstop, tinf, tsup, tree, lr, coeff);
		logLike1 = fillUncertaintyInfoGeneral(tree->node[tree->node[n].child].sibling, icur, imod, tstop, tinf, tsup, tree, lr, coeff);		
		res = (double*) malloc(lr[n].nLeaves*sizeof(double));
		like = (double*) malloc(lr[n].nLeaves*sizeof(double));
		offset = (double*) malloc(lr[n].nLeaves*sizeof(double));
		for(i=0; i<lr[n].nLeaves; i++) {
			offset[i] = NEG_INFTY;
			like[i] = 0.;
		}
		if(tsup[n] > tstop[icur])
			res[0] = getLogProbabilityDensityConstraintGeneral(n, icur, imod, tstop, tinf, tsup, tree, lr, coeff);
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

double getLogProbabilityDensityConstraintGeneral(int n, int icur, int imod, double *tstop, double *tinf, double *tsup, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeCoeff *coeff) {
	if(tstop[icur+1] < tree->maxTime) {
		if(tinf[n] > tstop[icur+1]) {
			if(tstop[icur+1] >= coeff->model.startTime[imod+1])		
				return getLogProbX(tstop[icur], tstop[icur+1], imod, 1, coeff)+getLogProbabilityDensityConstraintGeneral(n, icur+1, imod+1, tstop, tinf, tsup, tree, lr, coeff);
			else
				return getLogProbX(tstop[icur], tstop[icur+1], imod, 1, coeff)+getLogProbabilityDensityConstraintGeneral(n, icur+1, imod, tstop, tinf, tsup, tree, lr, coeff);
		} else {
			int i;
			double like = 0., offset = NEG_INFTY, *ll;
			if(tstop[icur+1] >= coeff->model.startTime[imod+1])
				ll = fillUncertaintyInfoGeneral(n, icur+1, imod+1, tstop, tinf, tsup, tree, lr,  coeff);
			else
				ll = fillUncertaintyInfoGeneral(n, icur+1, imod, tstop, tinf, tsup, tree, lr,  coeff);
			for(i=0; i<lr[n].nLeaves; i++) {
				double logLike;
				if(ll[i] != NEG_INFTY)
					logLike = getLogProbX(tstop[icur], tstop[icur+1], imod, i+1, coeff) + ll[i] - gsl_sf_lnfact(i);
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
	} else {
		return getLogProbX(tstop[icur], tree->maxTime, imod, lr[n].nLeaves,  coeff)+((double)lr[n].nLeaves)*getLogProbObs(tree->maxTime, imod, coeff)+getTreeShapeLogProbabilityDensity(lr[n].logRank, lr[n].nLeaves);
	}
}

void fillInfSup(int n, double tcur, double *tinf, double *tsup, double *tmin, double *tmax, TypeTree *tree) {
	if(tree->node[n].child == NOSUCH) {
		tinf[n] = tree->maxTime;
		tsup[n] = tree->maxTime;
	} else {
		int c;
		tinf[n] = (tmin[n]!=NO_TIME)?utils_MAX(tmin[n], tcur):tcur;
		tsup[n] = (tmax[n]!=NO_TIME)?tmax[n]:tree->maxTime;
		for(c=tree->node[n].child; c!=NOSUCH; c=tree->node[c].sibling) {
			fillInfSup(c, tinf[n], tinf, tsup, tmin, tmax, tree);
			tsup[n] = utils_MIN(tsup[n], tsup[c]);
		}
	}
}


double *fillUncertaintyInfoShift(int n, double val, int *shift, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeCoeff *coeffA, TypeUncertaintyLikeCoeff *coeffB) {
	double *res;
	if(shift[n]) {
		int i;
		res = (double*) malloc(lr[n].nLeaves*sizeof(double));
		res[0] = getLogProbTypeA(val, tree->maxTime, lr[n].nLeaves,  coeffB)+getTreeShapeLogProbabilityDensity(lr[n].logRank, lr[n].nLeaves);
		for(i=1; i<lr[n].nLeaves; i++) 
			res[i] = NEG_INFTY;
	} else {
		if(tree->node[n].child == NOSUCH) {
			res = (double*) malloc(sizeof(double));
			res[0] = getLogProbTypeA(val, tree->maxTime, 1,  coeffA);
		} else {
			int i, j;
			double *logLike0, *logLike1, *like, *offset;
			logLike0 = fillUncertaintyInfoShift(tree->node[n].child, val, shift, tree, lr, coeffA, coeffB);
			logLike1 = fillUncertaintyInfoShift(tree->node[tree->node[n].child].sibling, val, shift, tree, lr, coeffA, coeffB);		
			res = (double*) malloc(lr[n].nLeaves*sizeof(double));
			like = (double*) malloc(lr[n].nLeaves*sizeof(double));
			offset = (double*) malloc(lr[n].nLeaves*sizeof(double));
			for(i=0; i<lr[n].nLeaves; i++) {
				offset[i] = NEG_INFTY;
				like[i] = 0.;
			}
			res[0] = getLogProbTypeA(val, tree->maxTime, lr[n].nLeaves,  coeffA)+getTreeShapeLogProbabilityDensity(lr[n].logRank, lr[n].nLeaves);
			for(i=0; i<lr[tree->node[n].child].nLeaves; i++) 
				if(logLike0[i] != NEG_INFTY)
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
	}
	return res;
}

double getLogProbabilityDensityShift(double val, int *shift, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeCoeff *coeffA, TypeUncertaintyLikeCoeff *coeffB) {
	int i;
	double like = 0., offset = NEG_INFTY, *ll = fillUncertaintyInfoShift(tree->root, val, shift, tree, lr, coeffA, coeffB);
	for(i=0; i<lr[tree->root].nLeaves; i++) {
		double logLike;
		if(ll[i] != NEG_INFTY)
			logLike = getLogProbStart(tree->minTime, val, tree->maxTime, i+1, coeffA) + ll[i] - gsl_sf_lnfact(i);
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

/*return p(k, startTime)*/
double getLogProbTypeA(double startTime, double maxTime, int k,  TypeUncertaintyLikeCoeff *coeff) {
    double time = maxTime-startTime;
    if(k == 1)
		return
			log(coeff->model.param[0].sampl)
			+2.*(log(coeff->bmd[0])-log(coeff->bs[0]+coeff->bsd[0]*exp(-coeff->bmd[0]*time)))
			-coeff->bmd[0]*time;
	else
		return
			((double)k)*log(coeff->model.param[0].sampl)
			+2.*log(coeff->bmd[0])
			-coeff->bmd[0]*time
			+((double)k-1.)*(log(coeff->model.param[0].birth)+log(1.-exp(-coeff->bmd[0]*time)))
			-((double)k+1.)*log(coeff->bs[0]+coeff->bsd[0]*exp(-coeff->bmd[0]*time));
}

double getLogProbX(double s, double e, int i, int k, TypeUncertaintyLikeCoeff *coeff) {
	if(k == 1)
		return 2.*log(coeff->bmd[i])
		-coeff->bmd[i]*(e-s)
		-2.*log(coeff->model.param[i].birth*getProbObs(e, i, coeff)+(coeff->model.param[i].birth*getProbNotObs(e, i, coeff)-coeff->model.param[i].death)*exp(-coeff->bmd[i]*(e-s)));
	else
		return
		2.*log(coeff->bmd[i])
		-coeff->bmd[i]*(e-s)
		+((double)k-1.)*(log(coeff->model.param[i].birth)+log(1.-exp(-coeff->bmd[i]*(e-s))))
		-((double)k+1.)*log(coeff->model.param[i].birth*getProbObs(e, i, coeff)+(coeff->model.param[i].birth*getProbNotObs(e, i, coeff)-coeff->model.param[i].death)*exp(-coeff->bmd[i]*(e-s)));
}

double getLogProbY(double s, double e, int i, int k, TypeUncertaintyLikeCoeff *coeff) {
	if(k == 1)
		return 2.*log(coeff->bmd[i])
		-coeff->bmd[i]*(e-s)
		-((double)k+1.)*(log(1-coeff->model.param[i].birth*exp(-coeff->bmd[i]*(e-s)))+getLogProbNotObs(e, i, coeff));
	else
		return 2.*log(coeff->bmd[i])
		-coeff->bmd[i]*(e-s)
		+((double)k-1.)*(log(coeff->model.param[i].birth)+log(1.-exp(-coeff->bmd[i]*(e-s))))
		-((double)k+1.)*(log(1-coeff->model.param[i].birth*exp(-coeff->bmd[i]*(e-s)))+getLogProbNotObs(e, i, coeff));
}

double getLogProbObs(double t, int i, TypeUncertaintyLikeCoeff *coeff) {
	return 
		log(coeff->bmd[i])
		+log(1.-coeff->c[i])
		-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-coeff->model.startTime[i+1])
		-log(coeff->model.param[i].birth*(1.-coeff->c[i])*exp(-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-coeff->model.startTime[i+1]))+(coeff->c[i]*coeff->model.param[i].birth-coeff->model.param[i].death)*exp(-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-t)));
}

double getProbObs(double t, int i, TypeUncertaintyLikeCoeff *coeff) {
	return (coeff->bmd[i]*(1.-coeff->c[i])*exp(-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-coeff->model.startTime[i+1])))/(coeff->model.param[i].birth*(1.-coeff->c[i])*exp(-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-coeff->model.startTime[i+1]))+(coeff->c[i]*coeff->model.param[i].birth-coeff->model.param[i].death)*exp(-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-t)));
}

double getLogProbNotObs(double t, int i, TypeUncertaintyLikeCoeff *coeff) {
	return log(-(coeff->model.param[i].death*(coeff->c[i]-1.)*exp(-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-coeff->model.startTime[i+1]))+(coeff->model.param[i].death-coeff->c[i]*coeff->model.param[i].birth)*exp(-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-t))))
	-log(-(coeff->model.param[i].birth*(coeff->c[i]-1.)*exp(-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-coeff->model.startTime[i+1]))+(coeff->model.param[i].death-coeff->c[i]*coeff->model.param[i].birth)*exp(-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-t))));
}

double getProbNotObs(double t, int i, TypeUncertaintyLikeCoeff *coeff) {
	return (coeff->model.param[i].death*(coeff->c[i]-1.)*exp(-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-coeff->model.startTime[i+1]))+(coeff->model.param[i].death-coeff->c[i]*coeff->model.param[i].birth)*exp(-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-t)))/(coeff->model.param[i].birth*(coeff->c[i]-1.)*exp(-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-coeff->model.startTime[i+1]))+(coeff->model.param[i].death-coeff->c[i]*coeff->model.param[i].birth)*exp(-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-t)));
}

TypeUncertaintyLikeCoeff getUncertaintyLikeCoeff(TypePiecewiseModelParam *param) {
	TypeUncertaintyLikeCoeff res;
	int i;
	res.model = *param;
	res.c = (double*) malloc(res.model.size*sizeof(double));
	res.bmd = (double*) malloc(res.model.size*sizeof(double));
	res.bs = (double*) malloc(res.model.size*sizeof(double));
	res.bsd = (double*) malloc(res.model.size*sizeof(double));
	for(i=0; i<res.model.size; i++) {
		res.bmd[i] = res.model.param[i].birth-res.model.param[i].death;
		res.bs[i] = res.model.param[i].birth*res.model.param[i].sampl;
		res.bsd[i] = res.model.param[i].birth*(1.-res.model.param[i].sampl)-res.model.param[i].death;
	}
	res.c[res.model.size-1] = 1.-res.model.param[res.model.size-1].sampl;
	for(i=res.model.size-2; i>=0; i--)
		res.c[i] = 1.-res.model.param[i].sampl*getProbObs(res.model.startTime[i+1],i+1, &res);
	return res;
}

double logBinomial(unsigned int k, unsigned int n) {
    return gsl_sf_lnfact(n)-(gsl_sf_lnfact(k)+gsl_sf_lnfact(n-k));
}

double getTreeShapeLogProbabilityDensity(double rank, unsigned int nleaf) {
    return rank-gsl_sf_lnfact(nleaf-1);
}

/* return P_y(k, startTime, endTime)/(P^\star_o)^(k-l)*/
double getLogProbStart(double startTime, double endTime, double maxTime, int k, TypeUncertaintyLikeCoeff *coeff) {
    if(k == 1)
		return
			-coeff->bmd[0]*(endTime-startTime)
			+2.*(log(coeff->bs[0]+coeff->bsd[0]*exp(-coeff->bmd[0]*(maxTime-endTime)))-log(coeff->bs[0]+coeff->bsd[0]*exp(-coeff->bmd[0]*(maxTime-startTime))));
	else
		return
			-coeff->bmd[0]*(endTime-startTime)
			+((double)k-1.)*(log(coeff->model.param[0].birth)+log(1.-exp(-coeff->bmd[0]*(endTime-startTime)))-log(coeff->bmd[0]))
			+((double)k+1)*(log(coeff->bs[0]+coeff->bsd[0]*exp(-coeff->bmd[0]*(maxTime-endTime)))-log(coeff->bs[0]+coeff->bsd[0]*exp(-coeff->bmd[0]*(maxTime-startTime))));
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

double ghat(int i, double t, TypeUncertaintyLikeCoeff *coeff) {
	return exp(-coeff->bsd[i]*(2*coeff->model.startTime[coeff->model.size]-t-coeff->model.startTime[i+1]))/pow((coeff->model.param[i].birth*(coeff->c[i]-1.)*exp(-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-coeff->model.startTime[i+1]))+(coeff->model.param[i].death-coeff->c[i]*coeff->model.param[i].birth)*exp(-coeff->bmd[i]*(coeff->model.startTime[coeff->model.size]-t))), 2.);
}

double probSingle(int i, double tstart, int j, double tend, TypeUncertaintyLikeCoeff *coeff) {
	int k;
	double res;
	res = ghat(i, tstart, coeff);
	for(k=i; k<j; k++)
		res *= coeff->model.param[k].sampl*pow(coeff->bsd[k], 2.)*ghat(k+1, coeff->model.startTime[k+1], coeff);
	res /= ghat(j, tend, coeff);
	if(tend == coeff->model.startTime[j+1])
		res *= coeff->model.param[j].sampl;
	return res;
}

double getDistributionRoot(int j, double val, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeCoeff *coeff) {
	TypePiecewiseModelParam piece;
	TypeUncertaintyLikeCoeff coeffBis;
	double res, *tinf, *tsup;
	int k;
	double ll, ll0;

	piece.size = coeff->model.size-j;
	piece.startTime = (double*) malloc((piece.size+1)*sizeof(double));
	piece.param = &(coeff->model.param[j]);
	piece.startTime[0] = val;
	for(k=1; k<=piece.size; k++)
		piece.startTime[k] = coeff->model.startTime[k+j];
	coeffBis = getUncertaintyLikeCoeff(&piece);
	if(tree->parent == NULL)
		setParent(tree);
	tinf = (double*) malloc(tree->size*sizeof(double));
	tsup = (double*) malloc(tree->size*sizeof(double));
	for(k=0; k<tree->size; k++)
		if(tree->node[k].child == NOSUCH) {
			tinf[k] = tree->maxTime;
			tsup[k] = tree->maxTime;
		} else {
			tinf[k] = tree->minTime;
			tsup[k] = tree->maxTime;
		}
	ll0 = getLogProbabilityDensityConstraintGeneral(tree->root, 0, 0, coeff->model.startTime, tinf, tsup, tree, lr, coeff);
	ll = getLogProbabilityDensityConstraintGeneral(tree->root, 0, 0, piece.startTime, tinf, tsup, tree, lr, &coeffBis);
	res = 1.-(probSingle(0, coeff->model.startTime[0], j, val, coeff)*exp(ll))/exp(ll0);
	free((void*)piece.startTime);
	free((void*)tinf);
	free((void*)tsup);
	return res;
}

double logNLeaves(int n, int *nleaves, TypeTree *tree) {
	int c;
	double res = log((double) nleaves[n]-1.);
	for(c=tree->node[n].child; c!=NOSUCH; c=tree->node[c].sibling)
		if(tree->node[c].child != NOSUCH)
			res += logNLeaves(c, nleaves, tree);
	return res;
}

double logRank(int n, TypeTree *tree) {
	int *nleaves;
	double res;
	if(tree->node[n].child == NOSUCH)
		return 0.;
	nleaves = getNLeaves(tree);
	res = gsl_sf_lnfact((unsigned int)(nleaves[n]-1.))-logNLeaves(n, nleaves, tree);
	free((void*)nleaves);
	return res;
}

double logRankTot(TypeTree *tree) {
	int *nleaves, n;
	double res;
	if(tree->node[tree->root].child == NOSUCH)
		return 0.;
	nleaves = getNLeaves(tree);
	res = gsl_sf_lnfact((unsigned int)(nleaves[tree->root]-1));
	for(n=0; n<tree->size; n++)
		if(tree->node[n].child != NOSUCH)
			res -= log((double) nleaves[n]-1.);
	free((void*)nleaves);
	return res;
}

void fillUncertaintyLogNumberRankData(int n, TypeTree *tree, TypeUncertaintyLogNumberRankData *nr) {
	if(tree->node[n].child == NOSUCH) {
		nr[n].nLeaves = 1;
		nr[n].lnRank = 0.;
	} else {
		fillUncertaintyLogNumberRankData(tree->node[n].child, tree, nr);
		fillUncertaintyLogNumberRankData(tree->node[tree->node[n].child].sibling, tree, nr);		
		nr[n].nLeaves = nr[tree->node[n].child].nLeaves+nr[tree->node[tree->node[n].child].sibling].nLeaves;
		nr[n].lnRank = nr[tree->node[n].child].lnRank+nr[tree->node[tree->node[n].child].sibling].lnRank+logBinomial((unsigned int) (nr[tree->node[n].child].nLeaves-1), (unsigned int) (nr[n].nLeaves-2));
	}
}

double *probRank(int v, TypeTree *tree) {
	int *parent, m, n, u, *x;
	TypeUncertaintyLogNumberRankData *nr;
	nr = (TypeUncertaintyLogNumberRankData*) malloc(tree->size*sizeof(TypeUncertaintyLogNumberRankData));
	fillUncertaintyLogNumberRankData(tree->root, tree, nr);
	double **a, *res;
	int *rmin, *rmax;
	parent = getParent(tree);
	x = (int*) malloc(tree->size*sizeof(int));
	n = 0;
	u = v;
	do {
		x[n++] = u;
		u = parent[u];
	} while(u != NOSUCH);
	a = (double**) malloc(n*sizeof(double*));
	rmin = (int*) malloc(n*sizeof(int));
	rmax = (int*) malloc(n*sizeof(int));
	rmin[0] = 1;
	rmax[0] = 1;
	a[0] = (double*) malloc((rmax[0]-rmin[0]+1)*sizeof(double));
	a[0][0] = exp(nr[x[0]].lnRank);
	for(m=1; m<n; m++) {
		int s, nintm, nints, i, j;
		if(tree->node[x[m]].child == x[m-1])
			s = tree->node[x[m-1]].sibling;
		else
			s = tree->node[x[m]].child;
		nints = nr[s].nLeaves-1;
		nintm = nr[x[m-1]].nLeaves-1;
		rmin[m] = rmin[m-1]+1;
		rmax[m] = rmax[m-1]+nints+1;
		a[m] = (double*) malloc((rmax[m]-rmin[m]+1)*sizeof(double));
		for(i=0; i<=(rmax[m]-rmin[m]); i++)
			a[m][i] = 0.;
		for(i=0; i<=(rmax[m-1]-rmin[m-1]); i++)
			for(j=0; j<=nints; j++)
				a[m][i+j] += a[m-1][i]*exp(logBinomial((unsigned int) rmin[m-1]+i-1, (unsigned int) (rmin[m-1]+i-1+j))+logBinomial((unsigned int) (nintm-rmin[m-1]-i), (unsigned int) (nintm-i-rmin[m-1]+nints-j)));
		for(i=0; i<=(rmax[m]-rmin[m]); i++)
			a[m][i] *= exp(nr[s].lnRank);
	}
	res = (double*) malloc((nr[tree->root].nLeaves-1)*sizeof(double));
	for(m=1; m<rmin[n-1]; m++)
		res[m-1] = 0.;
	for(m=rmin[n-1]; m<=rmax[n-1]; m++)
		res[m-1] = a[n-1][m-rmin[n-1]]/exp(nr[tree->root].lnRank);
	for(m=rmax[n-1]+1; m<nr[tree->root].nLeaves; m++)
		res[m-1] = 0.;
	for(m=0; m<n; m++)
		free((void*)a[m]);
	free((void*)a);
	free((void*)x);
	free((void*)rmin);
	free((void*)rmax);
	free((void*)nr);
	return res;
}

double lfA(double s, double t, int k, int n, TypeUncertaintyLikeCoeff *coeff) {
	if(n<=k)
		return NEG_INFTY;
	return log((double)(n-k))+logBinomial((unsigned int) (n-k), (unsigned int) (n-1))+((double)n-k-1)*lF(s, t, coeff)+((double)k-1)*log(1-exp(lF(s, t, coeff)))+lf(s, t, coeff);
}

double lf(double s, double t, TypeUncertaintyLikeCoeff *coeff) {
	if(s>t)
		return NEG_INFTY;
	return 2.*log(coeff->bmd[0])-coeff->bmd[0]*s-2.*log(coeff->model.param[0].birth-coeff->model.param[0].death*exp(-coeff->bmd[0]*s))+log(coeff->model.param[0].birth-coeff->model.param[0].death*exp(-coeff->bmd[0]*t))-log(1.-exp(-coeff->bmd[0]*t));
}

double lF(double s, double t, TypeUncertaintyLikeCoeff *coeff) {
	if(s>t)
		return 0.;
	return log(1.-exp(-coeff->bmd[0]*s))-log(coeff->model.param[0].birth-coeff->model.param[0].death*exp(-coeff->bmd[0]*s))+log(coeff->model.param[0].birth-coeff->model.param[0].death*exp(-coeff->bmd[0]*t))-log(1.-exp(-coeff->bmd[0]*t));
}

double getDivergenceTimeDensityStadler(int v, double s, TypeTree *tree, TypePiecewiseModelParam *param) {
	TypeUncertaintyLikeCoeff coeff=getUncertaintyLikeCoeff(param);
	double *pr, res;
	int i, n;
	n = tree->size/2;
	pr = probRank(v, tree);
	res = 0.;
	for(i=1; i<=n; i++)
		res += pr[i-1]*exp(lfA(tree->maxTime-s, tree->maxTime, i, n+1, &coeff));
	free((void*)pr);
	return res;
}
