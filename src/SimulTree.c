#include "SimulTree.h"
#include <math.h>
#include "Utils.h"


#define INC_FOSSIL_ITEM 50
#define INFTY 9.9E99;

static int getSampled(gsl_rng *r, double p);
/* return  a random type of event wrt the rates*/
static char getType(gsl_rng *rgen, double birth, double death);
/* return a lineage uniformly drawn among 1,.., n*/
static int getWhich(gsl_rng *rgen, int n);


/* return 1 under p, 1-p*/
int getSampled(gsl_rng *rgen, double p) {
	if(gsl_rng_uniform(rgen) <= p)
		return 1;
	else
		return 0;
}

/* return  a random type of event wrt the rates*/
char getType(gsl_rng *rgen, double birth, double death) {
	double uni = gsl_rng_uniform(rgen);
	if(uni<birth/(birth+death))
		return 'b';
	else
		return 'd';
}

/* return a lineage uniformly drawn among 1,.., n*/
int getWhich(gsl_rng *rgen, int n) {
	return gsl_rng_uniform_int(rgen, n);
}


/*simulate a random tree with specified birth and death rates and fossil finds on this tree with rate "fossil"*/
TypeTree *simulTree(gsl_rng *rgen, double birth, double death, double maxTime) {
	TypeTree *tree;
	int *cur, ncur, i;
	double time = 0.;
	if((cur = (int*) malloc((MAX_CURRENT+1)*sizeof(int))) == NULL)
		return NULL;
	tree = newTree(INC_SIZE);
	tree->maxTime = maxTime;
	tree->minTime = 0.;
	tree->size = 1;
	tree->root = 0;
	tree->time[0] = INFTY;
	tree->node[0].child = NOSUCH;
	tree->node[0].sibling = NOSUCH;
	cur[0] = 0;
	ncur = 1;
	while(time<maxTime && ncur>0) {
		int type = getType(rgen, birth, death);
		int which = gsl_rng_uniform_int(rgen, ncur);
		double wait = gsl_ran_exponential(rgen, 1./(ncur*(birth+death)));
		time += wait;
		if(time < maxTime) {
			switch(type) {
				case 'b':
					if(ncur > MAX_CURRENT) {
						printf("too much lineages generated during simulations (%d - max %d)\n", ncur, MAX_CURRENT);
						freeTree(tree);
						return NULL;
					}
					tree->time[cur[which]] = time;
					if((tree->size+1)>=tree->sizeBuf) {
						tree->sizeBuf += INC_SIZE;
						tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
						tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
					}
					tree->node[cur[which]].child = tree->size;
					tree->time[tree->size] = INFTY;
					tree->node[tree->size].child = NOSUCH;
					tree->node[tree->size].sibling = tree->size+1;
					tree->time[tree->size+1] = INFTY;
					tree->node[tree->size+1].child = NOSUCH;
					tree->node[tree->size+1].sibling = NOSUCH;
					cur[which] = tree->size;
					cur[ncur] = tree->size+1;
					ncur++;
					tree->size += 2;
					break;
				case 'd':
					tree->time[cur[which]] = time;
					for(i=which+1; i<ncur; i++)
						cur[i-1] = cur[i];
					ncur--;
					break;
				default:
					break;
			}
		}
	}
	for(i=0; i<ncur; i++)
		tree->time[cur[i]] = maxTime;
	free((void*)cur);
	tree->sizeBuf = tree->size;
	if(tree->size) {
		tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
		tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
	} else {
		free((void*)tree->node);
		free((void*)tree->time);
		tree->node = NULL;
		tree->time = NULL;
	}
	return tree;
}
