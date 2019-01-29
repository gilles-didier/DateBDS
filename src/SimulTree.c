#include "SimulTree.h"
#include <math.h>
#include "Utils.h"


#define INC_FOSSIL_ITEM 50
#define INFTY 9.9E99;
#define MAX_SIZE_TREE 100000

static int getSampled(gsl_rng *r, double p);
/* return  a random type of event wrt the rates*/
static char getType(gsl_rng *rgen, double birth, double death);
static void simulSubTreeRec(gsl_rng *rgen, int n, TypeTree *tree, int **sampled, double startTime, double endTime, TypeModelParam *param);
static void simulSubTreeShiftRec(gsl_rng *rgen, int n, TypeTree *tree, int **sampled, int *shift, double shiftTime, double startTime, double endTime, TypeModelParam *paramA, TypeModelParam *paramB);


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


/*simulate a random tree with specified param->birth and param->death rates and fossil finds on this tree with rate "fossil"*/
TypeTree *simulTree(gsl_rng *rgen, double startTime, double endTime, TypeModelParam *param) {
	TypeTree *tree;
	int *sampled;
	tree = newTree(INC_SIZE);
	tree->maxTime = endTime;
	tree->minTime = startTime;
	tree->size = 1;
	tree->root = 0;
	tree->time[0] = INFTY;
	tree->node[0].child = NOSUCH;
	tree->node[0].sibling = NOSUCH;
	sampled = (int*) malloc(tree->sizeBuf*sizeof(int));
	simulSubTreeRec(rgen, tree->root, tree, &sampled, startTime, endTime, param);
	tree->sizeBuf = tree->size;
	if(tree->size) {
		TypeTree *tree1, *tree2;
		tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
		tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
		tree1 = pruneLeavesFromTable(tree, sampled);
		freeTree(tree);
		tree2 = fixBinary(tree1);
		freeTree(tree1);
		tree = tree2;
	} else {
		free((void*)tree->node);
		free((void*)tree->time);
		tree->node = NULL;
		tree->time = NULL;
	}
	free((void*)sampled);
	return tree;
}

/*simulate a random tree with specified param->birth and param->death rates and fossil finds on this tree with rate "fossil"*/
TypeTree *simulTreeShift(gsl_rng *rgen, int *shift, double shiftTime, double startTime, double endTime, TypeModelParam *paramA, TypeModelParam *paramB) {
	TypeTree *tree;
	int *sampled;
	tree = newTree(INC_SIZE);
	tree->maxTime = endTime;
	tree->minTime = startTime;
	tree->size = 1;
	tree->root = 0;
	tree->time[0] = INFTY;
	tree->node[0].child = NOSUCH;
	tree->node[0].sibling = NOSUCH;
	sampled = (int*) malloc(tree->sizeBuf*sizeof(int));
	*shift = NOSUCH;
	simulSubTreeShiftRec(rgen, tree->root, tree, &sampled, shift, shiftTime, startTime, endTime, paramA, paramB);
	tree->sizeBuf = tree->size;
	if(tree->size>0) {
		TypeTree *tree1, *tree2;
		int n;
		tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
		tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
		tree->name = (char**) malloc(tree->size*sizeof(char*));
		for(n=0; n<tree->size; n++)
			if(tree->node[n].child == NOSUCH)
				tree->name[n] = "o";
			else
				tree->name[n] = NULL;
		//for(n=0; n<tree->size; n++)
			//tree->name[n] = "o";
		if(*shift != NOSUCH)
			tree->name[*shift] = "*";
//printf("\nTree\n");
//printTreeDebug(stdout, tree->root, tree, tree->name);
		tree1 = pruneLeavesFromTable(tree, sampled);
		free((void*)tree->name);
		tree->name = NULL;
		freeTree(tree);
//printf("\nTree1\n");
//printTreeDebug(stdout, tree1->root, tree1, tree1->name);
		tree2 = fixBinary(tree1);
		freeTree(tree1);
//printf("\nTree2\n");
//printTreeDebug(stdout, tree2->root, tree2, tree2->name);
		tree = tree2;
	} else {
		free((void*)tree->node);
		free((void*)tree->time);
		tree->node = NULL;
		tree->time = NULL;
	}
	free((void*)sampled);
	return tree;
}

/*simulate a random tree with specified param->birth and param->death rates and fossil finds on this tree with rate "fossil"*/
void simulSubTreeRec(gsl_rng *rgen, int n, TypeTree *tree, int **sampled, double startTime, double endTime, TypeModelParam *param) {
	double time = startTime+gsl_ran_exponential(rgen, 1./(param->birth+param->death));
	if(time<endTime) {
		int type = getType(rgen, param->birth, param->death);
		switch(type) {
			case 'b':
				tree->time[n] = time;
				if(tree->size>=tree->sizeBuf) {
					tree->sizeBuf += INC_SIZE;
					tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
					tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
					*sampled = (int*) realloc((void*)*sampled, tree->sizeBuf*sizeof(int));
					if(tree->size>MAX_SIZE_TREE)
						printf("Warning size tree of %d\n", tree->size);
				}
				tree->node[n].child = tree->size;
				tree->time[tree->size] = INFTY;
				tree->node[tree->size].child = NOSUCH;
				tree->size++;
				simulSubTreeRec(rgen, tree->node[n].child, tree, sampled, time, endTime, param);
				if(tree->size>=tree->sizeBuf) {
					tree->sizeBuf += INC_SIZE;
					tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
					tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
					*sampled = (int*) realloc((void*)*sampled, tree->sizeBuf*sizeof(int));
				}
				tree->node[tree->node[n].child].sibling = tree->size;
				tree->time[tree->size] = INFTY;
				tree->node[tree->size].child = NOSUCH;
				tree->node[tree->size].sibling = NOSUCH;
				tree->size++;
				simulSubTreeRec(rgen, tree->node[tree->node[n].child].sibling, tree, sampled, time, endTime, param);
				break;
			case 'd':
				tree->time[n] = time;
				(*sampled)[n] = 0;
				break;
			default:
				break;
		}
	} else {
		tree->time[n] = endTime;
		(*sampled)[n] = getSampled(rgen, param->sampl);
	}	
}

/*simulate a random tree with specified param->birth and param->death rates and fossil finds on this tree with rate "fossil"*/
void simulSubTreeShiftRec(gsl_rng *rgen, int n, TypeTree *tree, int **sampled, int *shift, double shiftTime, double startTime, double endTime, TypeModelParam *paramA, TypeModelParam *paramB) {
	double time = startTime+gsl_ran_exponential(rgen, 1./(paramA->birth+paramA->death));
	if(*shift == NOSUCH && shiftTime>startTime && shiftTime<time) {
printf("node %d (root %d size %d) start %.2lf end %.2lf shift %.2lf\n", n, tree->root, tree->size, startTime, time, shiftTime);
		simulSubTreeRec(rgen, n, tree, sampled, shiftTime, endTime, paramB);
		*shift = n;
	} else { 
		if(time<endTime) {
			int type = getType(rgen, paramA->birth, paramA->death);
			switch(type) {
				case 'b':
					tree->time[n] = time;
					if(tree->size>=tree->sizeBuf) {
						tree->sizeBuf += INC_SIZE;
						tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
						tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
						*sampled = (int*) realloc((void*)*sampled, tree->sizeBuf*sizeof(int));
						if(tree->size>MAX_SIZE_TREE)
							printf("Warning size tree of %d\n", tree->size);
					}
					tree->node[n].child = tree->size;
					tree->time[tree->size] = INFTY;
					tree->node[tree->size].child = NOSUCH;
					tree->size++;
					if(*shift == NOSUCH)
						simulSubTreeShiftRec(rgen, tree->node[n].child, tree, sampled, shift, shiftTime, time, endTime, paramA, paramB);
					else
						simulSubTreeRec(rgen, tree->node[n].child, tree, sampled, time, endTime, paramA);
					if(tree->size>=tree->sizeBuf) {
						tree->sizeBuf += INC_SIZE;
						tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
						tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
						*sampled = (int*) realloc((void*)*sampled, tree->sizeBuf*sizeof(int));
					}
					tree->node[tree->node[n].child].sibling = tree->size;
					tree->time[tree->size] = INFTY;
					tree->node[tree->size].child = NOSUCH;
					tree->node[tree->size].sibling = NOSUCH;
					tree->size ++;
					if(*shift == NOSUCH)
						simulSubTreeShiftRec(rgen, tree->node[tree->node[n].child].sibling, tree, sampled, shift, shiftTime, time, endTime, paramA, paramB);
					else
						simulSubTreeRec(rgen, tree->node[tree->node[n].child].sibling, tree, sampled, time, endTime, paramA);
					break;
				case 'd':
					tree->time[n] = time;
					(*sampled)[n] = 0;
					break;
				default:
					break;
			}
		} else {
			tree->time[n] = endTime;
			(*sampled)[n] = getSampled(rgen, paramA->sampl);			
		}
	}
}	

void switchNodes(int a, int b, int anc, TypeTree *tree) {
	if(tree->node[anc].child == a) {
		tree->node[anc].child = b;
		tree->node[b].sibling = tree->node[a].sibling;
	} else
		tree->node[tree->node[anc].child].sibling = b;
}

int getNL(int n, TypeTree *tree) {
	if(tree->node[n].child == NOSUCH)
		return 1;
	else {
		int c, r=0;
		for(c=tree->node[n].child; c != NOSUCH; c=tree->node[c].sibling)
			r += getNL(c, tree);
		return r;
	}
}

void addTree(gsl_rng *rgen, double shiftTime, int *old, int *new, TypeTree *tree, TypeModelParam *param) {
	int *cur, ncur, *sampled, *parent;
	double time;
	if(tree == NULL || tree->size == 0)
		return;
	cur = (int*) malloc(tree->size*sizeof(int));
	cur[0] = tree->root;
	time = tree->time[cur[0]];
//printTreeDebug(stdout, tree->root, tree, NULL);
//printf("size %d %.3lf/%.3lf\n", tree->size, tree->time[cur[0]], shiftTime);
	ncur = 1;
	while(time<shiftTime && ncur>0) {
		int i, min = 0;
		for(i=1; i<ncur; i++)
			if(tree->time[cur[i]]<tree->time[cur[min]])
				min = i;
//printf("%d cur:", ncur);
//for(i=0; i<ncur; i++)
	//printf(" %d", cur[i]);
//printf(" %.3lf/%.3lf\n", tree->time[cur[min]], shiftTime);

		if(tree->node[cur[min]].child == NOSUCH) {
			for(i=min; i<ncur-1; i++)
				cur[i] = cur[i+1];
			ncur--;
		} else {
			int n = cur[min];
			cur[min] = tree->node[n].child;
			cur[ncur++] = tree->node[tree->node[n].child].sibling;
		}
		time = tree->time[0];
		for(i=1; i<ncur; i++)
			if(tree->time[cur[i]]<time)
				time = tree->time[cur[i]];		
	}
	
	//printf("nodes (%d, %d):", ncur, tree->size);
	//for(i=0; i<ncur; i++)
		//printf(" %d-%d", cur[i], getNL(cur[i],tree));
	//printf("\n");
	if(cur[0] != tree->root) {
		parent = getParent(tree);
		if(tree->size>=tree->sizeBuf) {
			tree->sizeBuf += INC_SIZE;
			tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
			tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
		}
		sampled = (int*) malloc(tree->sizeBuf*sizeof(int));
		*old = cur[gsl_rng_uniform_int(rgen, ncur)];
		*new = tree->size++;
		tree->time[*new] = INFTY;
		tree->node[*new].child = NOSUCH;
		tree->node[*new].sibling = NOSUCH;
	//printf("Size %d/%d\n", tree->size, tree->sizeBuf);
		simulSubTreeRec(rgen, *new, tree, &sampled, shiftTime, tree->maxTime, param);
		switchNodes(*old, *new, parent[*old], tree);
		free((void*)parent);
		free((void*)sampled);
	} else {
		*old = tree->root;
		*new = tree->root;
		printf("Warning cut at root\n");
	}
	//printf(" %d-%d -- %d-%d\n", *old, getNL(*old,tree), *new, getNL(*new,tree));
	free((void*)cur);
}
				
		
	
///*simulate a random tree with specified param->birth and param->death rates and fossil finds on this tree with rate "fossil"*/
//TypeTree *simulTree(gsl_rng *rgen, double startTime, double endTime, TypeModelParam *param) {
	//TypeTree *tree;
	//int *cur, ncur, i;
	//double time = startTime;
	//if((cur = (int*) malloc((MAX_CURRENT+1)*sizeof(int))) == NULL)
		//return NULL;
	//tree = newTree(INC_SIZE);
	//tree->maxTime = endTime;
	//tree->minTime = startTime;
	//tree->size = 1;
	//tree->root = 0;
	//tree->time[0] = INFTY;
	//tree->node[0].child = NOSUCH;
	//tree->node[0].sibling = NOSUCH;
	//cur[0] = 0;
	//ncur = 1;
	//while(time<endTime && ncur>0) {
		//int type = getType(rgen, param->birth, param->death);
		//int which = gsl_rng_uniform_int(rgen, ncur);
		//double wait = gsl_ran_exponential(rgen, 1./(ncur*(param->birth+param->death)));
		//time += wait;
		//if(time < endTime) {
			//switch(type) {
				//case 'b':
					//if(ncur > MAX_CURRENT) {
						//printf("too much lineages generated during simulations (%d - max %d)\n", ncur, MAX_CURRENT);
						//freeTree(tree);
						//return NULL;
					//}
					//tree->time[cur[which]] = time;
					//if((tree->size+1)>=tree->sizeBuf) {
						//tree->sizeBuf += INC_SIZE;
						//tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
						//tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
					//}
					//tree->node[cur[which]].child = tree->size;
					//tree->time[tree->size] = INFTY;
					//tree->node[tree->size].child = NOSUCH;
					//tree->node[tree->size].sibling = tree->size+1;
					//tree->time[tree->size+1] = INFTY;
					//tree->node[tree->size+1].child = NOSUCH;
					//tree->node[tree->size+1].sibling = NOSUCH;
					//cur[which] = tree->size;
					//cur[ncur] = tree->size+1;
					//ncur++;
					//tree->size += 2;
					//break;
				//case 'd':
					//tree->time[cur[which]] = time;
					//for(i=which+1; i<ncur; i++)
						//cur[i-1] = cur[i];
					//ncur--;
					//break;
				//default:
					//break;
			//}
		//}
	//}
	//for(i=0; i<ncur; i++)
		//tree->time[cur[i]] = endTime;
	//free((void*)cur);
	//tree->sizeBuf = tree->size;
	//if(tree->size) {
		//int *keep, n;
		//TypeTree *tree1, *tree2;
		//tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
		//tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
		//keep = malloc(tree->size*sizeof(int));
		//i = 0;
		//for(n=0; n<tree->size; n++)
			//if(tree->node[n].child == NOSUCH)
				//keep[n] = getSampled(rgen, param->sampl);
		//keep[i] = NOSUCH;
		//tree1 = pruneLeavesFromTable(tree, keep);
		//free((void*)keep);
		//freeTree(tree);
		//tree2 = fixBinary(tree1);
		//freeTree(tree1);
		//tree = tree2;
	//} else {
		//free((void*)tree->node);
		//free((void*)tree->time);
		//tree->node = NULL;
		//tree->time = NULL;
	//}
	//return tree;
//}
