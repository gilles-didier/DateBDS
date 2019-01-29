#ifndef SimulTreeF
#define SimulTreeF

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Tree.h"
#include "Model.h"



/*simulate a random tree with specified birth and death rates and fossil finds on this tree with rate "fossil"*/
TypeTree *simulTree(gsl_rng *rgen, double startTime, double endTime, TypeModelParam *param);
TypeTree *simulTreeShift(gsl_rng *rgen, int *shift, double shiftTime, double startTime, double endTime, TypeModelParam *paramA, TypeModelParam *paramB);
void addTree(gsl_rng *rgen, double shiftTime, int *old, int *new, TypeTree *tree, TypeModelParam *param);
void switchNodes(int a, int b, int anc, TypeTree *tree);

#endif
