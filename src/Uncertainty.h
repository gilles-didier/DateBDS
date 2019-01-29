#ifndef UncertaintyF
#define UncertaintyF

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include "Utils.h"
#include "Tree.h"
#include "Model.h"
#include "Distribution.h"


typedef struct UNCERTAINTY_LR_DATA {
	int nLeaves;
	double logRank;
} TypeUncertaintyLogRankData;


typedef struct UNCERTAINTY_LIKE_PARAMETERS {
    TypePiecewiseModelParam model;
    double *c, *bmd, *bs, *bsd;
} TypeUncertaintyLikeCoeff;


#ifdef __cplusplus
extern "C" {
#endif
void fillUncertaintyLogRankData(int node, TypeTree *tree, TypeUncertaintyLogRankData *lr);
TypeUncertaintyLikeCoeff getUncertaintyLikeCoeff(TypePiecewiseModelParam *param);
TypeDistribution getDistributionUn(int n, double val, TypeTree *tree, TypePiecewiseModelParam *param);
TypeDistribution getDistribution(int n, int def, TypeTree *tree, TypePiecewiseModelParam *param);
void fillDistributionAll(int def, TypeTree *tree, TypePiecewiseModelParam *param, int *todo, TypeDistribution *d);
double getDistributionRoot(int j, double val, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeCoeff *coeff);
double getLogSumProbabilityDensityShift(int node, TypeTree *tree, TypePiecewiseModelParam *paramA, TypePiecewiseModelParam *paramB);
double getLogProbabilityDensityShiftX(int node, double val, TypeTree *tree, TypePiecewiseModelParam *paramA, TypePiecewiseModelParam *paramB);
double getLogProbabilityDensity(TypeTree *tree, TypePiecewiseModelParam *param);
double *probRank(int v, TypeTree *tree);
double getDivergenceTimeDensityStadler(int v, double s, TypeTree *tree, TypePiecewiseModelParam *param);
void fillDistributionAllGeneral(int def, double *tmin, double *tmax, TypeTree *tree, TypePiecewiseModelParam *param, int *todo, TypeDistribution *d);
double logProbabilityDensityConstraintGeneral(double *tmin, double *tmax, TypeTree *tree, TypePiecewiseModelParam *param);

#ifdef __cplusplus
}
#endif


#endif
