#include <stdlib.h>
#include <math.h>

#include "Model.h"
#include "TreeExtras.h"


//static double factorial(int n);

//static double real_fossil = 0.1, real_sample = 1.;

//static int minContemp = 20;

//static double getRate(double birth, double death, double fossil);
//static double getWaitSamFosProbLog(int n, double t, double x, double birth, double death, double fossil, double sample, double rate);
//static double getEndSamFosProbLog(int n, double t, double birth, double death, double fossil, double sample, double rate);
//static double getObsSamFosProb(double t, double birth, double death, double fossil, double sample, double rate);
//static double getObsSamFosProbLog(double t, double birth, double death, double fossil, double sample, double rate);
//static double getBirthSamFosProbLog(double t, double birth, double death, double fossil, double sample, double rate);
//static double getDeathSamFosProbLog(double t, double birth, double death, double fossil, double sample, double rate);
//static double getFossilSamFosProbLog(double t, double birth, double death, double fossil, double sample, double rate);

//static double getWaitFossilProbLog(int n, double t, double x, double birth, double death, double fossil, double alpha, double beta);
//static double getEndFossilProbLog(int n, double t, double birth, double death, double fossil, double alpha, double beta);

//static double getWaitCompleteProbLog(int n, double t, double x, double birth, double death, double fossil);
//static double getBirthCompleteProbLog(double birth, double death, double fossil);
//static double getDeathCompleteProbLog(double birth, double death, double fossil);
//static double getFossilCompleteProbLog(double birth, double death, double fossil);

//static double getWaitPhyloProbLog(int n, double t, double x, double birth, double death);
//static double getEndPhyloProbLog(int n, double t,double birth, double death);
//static double getfuncULog(double t, double birth, double death);
//static double getfuncPLog(double t, double birth, double death);

//static double getWaitPureProbLog(int n, double t, double x, double birth, double death);
//static double getBirthPureProbLog(double birth, double death);
//static double getDeathPureProbLog(double birth, double death);

//static double factorial(int n);

//double factorial(int n) {
    //int i;
    //double res = 1.;
    //for(i=1.; i<=n; i++)
        //res *= (double) i;
    //return res;
//}


///*get the total time*/
//double getTotalTimeEvent(TypeListEvent *event) {
    //double tot = 0., tprec = 0.;
    //int n;
    //for(n=0; n<event->size; n++) {
        //tot += (event->list[n].time-tprec)*event->list[n].n;
        //tprec = event->list[n].time;
    //}
    //if(event->size>0) {
        //n = event->list[event->size-1].n;
        //if(event->list[event->size-1].type == 'b')
            //n++;
        //if(event->list[event->size-1].type == 'd')
            //n--;
    //} else
        //n = 1;
    //return tot+n*(event->maxTime-tprec);
//}

///*returns the sequence of speciation/extinction/fossil find events occurring in "tree", chronologically ordered in ascending order*/
///*event.list[i].n refers to the number of lineages alive just before the ith event occurs*/
///*fossilList has to be sorted in ascending order*/
//TypeListEvent *getEventSequence(TypeTree *tree, TypeFossilFeature *fos) {
    //double time;
    //int cur[MAX_CURRENT], ncur, i, ind = 0, nbuf = INC_SIZE;
    //TypeListEvent *res;

    //res = (TypeListEvent*) malloc(sizeof(TypeListEvent));
    //res->size = 0;
    //if(tree->size == 0) {
        //res->list = NULL;
        //return res;
    //}
    //res->list = (TypeEvent*) malloc(nbuf*sizeof(TypeEvent));
    //if(tree->time[tree->root] == 0.) {
        //ncur = 2;
        //cur[0] = tree->node[tree->root].child;
        //cur[1] = tree->node[tree->node[tree->root].child].sibling;
    //} else {
        //ncur = 1;
        //cur[0] = tree->root;
    //}
    //time = 0.;
    //while(time<tree->maxTime) {
        //double tf, tbd;
        //int which;
        //tf = tree->maxTime;
        //if(fos != NULL && ind<fos->size)
            //tf = fos->fossilList[ind].time;
        //tbd = tree->maxTime;
        //for(i=0; i<ncur; i++) {
            //if(tree->time[cur[i]] < tbd) {
                //tbd = tree->time[cur[i]];
                //which = i;
            //}
        //}
        //time = utils_MIN(tf,tbd);
        //if(time<tree->maxTime) {
            //if(res->size>=nbuf) {
                //nbuf += INC_SIZE;
                //res->list = (TypeEvent*) realloc((void*)res->list, nbuf*sizeof(TypeEvent));
            //}
            //if(tf<tbd) {
                //res->list[res->size].time = tf;
                //res->list[res->size].n = ncur;
                //res->list[res->size].type = 'f';
                //(res->size)++;
                //ind++;
            //} else {
                //if(tf==tbd)
                    //ind++;
                //res->list[res->size].time = tbd;
                //res->list[res->size].n = ncur;
                //if(tree->node[cur[which]].child>=0 && tree->node[tree->node[cur[which]].child].sibling>=0) {
                    //if(ncur >= MAX_CURRENT)
                        //break;
                    //res->list[res->size].type = 'b';
                    //(res->size)++;
                    //cur[ncur] = tree->node[tree->node[cur[which]].child].sibling;
                    //cur[which] = tree->node[cur[which]].child;
                    //ncur++;
                //} else {
                    //res->list[res->size].type = 'd';
                    //(res->size)++;
                    //for(i=which+1; i<ncur; i++)
                        //cur[i-1] = cur[i];
                    //ncur--;
                //}
            //}
        //}
    //}
    //if(res->size)
        //res->list = (TypeEvent*) realloc((void*)res->list, (res->size)*sizeof(TypeEvent));
    //else {
        //free((void*)res->list);
        //res->list = NULL;
    //}
    //res->maxTime = tree->maxTime;
    //return res;
//}

///*shift all the times to set first event times to 0*/
//void offsetTimeEvents(TypeListEvent *event) {
    //int i;
    //double zero;

    //if(event->size ==0)
        //return;
    //zero = event->list[0].time;
    //for(i=0; i<event->size; i++)
        //event->list[i].time -= zero;
    //event->maxTime -= zero;
//}

///*set the number of events variables b, d, f for each type*/
//void getStatEvent(TypeListEvent *event, int *b,  int *d,  int *f) {
    //int n;

    //*b = 0;
    //*d = 0;
    //*f = 0;
    //for(n=0; n<event->size; n++) {
        //switch(event->list[n].type) {
            //case 'b':
                //(*b)++;
                //break;
            //case 'd':
                //(*d)++;
                //break;
            //case 'f':
                //(*f)++;
                //break;
            //default:
            //;
        //}
    //}
//}

///*free list of events*/
//void freeListEvent(TypeListEvent *event) {
    //if(event->size)
        //free((void*)event->list);
    //free((void*)event);
//}

///*print event*/
//void fprintEvent(FILE *f, TypeEvent *event) {
    //fprintf(f, "%c\t%d\t%.2lE", event->type, event->n, event->time);
//}

///*print list of events*/
//void fprintListEvent(FILE *f, TypeListEvent *event) {
    //int i;
    //for(i=0; i<event->size; i++) {
        //fprintEvent(f, &(event->list[i]));
        //fprintf(f, "\n");
    //}
//}


///*fossil reconstruction likelihood*/
//double getAlpha(double birth, double death, double fossil) {
    //return (birth+death+fossil-sqrt(pow(birth+death+fossil, 2.)-4.*birth*death))/(2*birth);
//}

//double getBeta(double birth, double death, double fossil) {
    //return (birth+death+fossil+sqrt(pow(birth+death+fossil, 2.)-4.*birth*death))/(2*birth);
//}

//double getWaitFossilProbLog(int n, double t, double x, double birth, double death, double fossil, double alpha, double beta) {
    //return log(n)-((double)n)*(birth*(1.-alpha)+fossil)*x+((double)n-1.)*log(beta-alpha*exp(-birth*(beta-alpha)*(t-x)))-((double)n)*log(beta-alpha*exp(-birth*(beta-alpha)*t))+log(beta*(birth*(1.-alpha)+fossil)-alpha*(birth*(1.-beta)+fossil)*exp(-birth*(beta-alpha)*(t-x)));
//}

//double getEndFossilProbLog(int n, double t, double birth, double death, double fossil, double alpha, double beta) {
    //return ((double)n)*(-(birth*(1.-alpha)+fossil)*t+log(beta-alpha)-log(beta-alpha*exp(-birth*(beta-alpha)*t)));
//}

//double getObsProb(double t, double birth, double death, double fossil, double alpha, double beta) {
    //return 1.-(alpha*beta*(1.-exp(-birth*(beta-alpha)*t)))/(beta-alpha*exp(-birth*(beta-alpha)*t));
//}

//double getLogCompProbObs(double t, double birth, double alpha, double beta) {
    //return log(alpha)+log(beta)+log((1.-exp(-birth*(beta-alpha)*t)))-log(beta-alpha*exp(-birth*(beta-alpha)*t));
//}

//double getLogProbObs(double t, double birth, double alpha, double beta) {
    //if(t==0.)
        //return 0.;
    //return log(1.-exp(getLogCompProbObs(t, birth, alpha, beta)));
//}

//double getObsProbLog(double t, double birth, double death, double fossil, double alpha, double beta) {
    //return log(1.-(alpha*beta*(1.-exp(-birth*(beta-alpha)*t)))/(beta-alpha*exp(-birth*(beta-alpha)*t)));
//}

//double getBirthProbLog(double t, double birth, double death, double fossil, double alpha, double beta) {
    //return log(birth)+ getObsProbLog(t, birth, death, fossil, alpha, beta)-log(birth*getObsProb(t, birth, death, fossil, alpha, beta)+fossil);
//}

//double getDeathProbLog(double t, double birth, double death, double fossil, double alpha, double beta) {
    //return log(fossil)+log(1.-getObsProb(t, birth, death, fossil, alpha, beta))-log(birth*getObsProb(t, birth, death, fossil, alpha, beta)+fossil);
//}

//double getFossilProbLog(double t, double birth, double death, double fossil, double alpha, double beta) {
    //return log(fossil) + getObsProbLog(t, birth, death, fossil, alpha, beta)-log(birth*getObsProb(t, birth, death, fossil, alpha, beta)+fossil);
//}

//double getLogLikelihoodWithoutFossil(double startTime, double endTime, double maxTime, int k,  double birth, double death, double fossil) {
    //double alpha, beta, D, E;
    //alpha = getAlpha(birth, death,fossil);
    //beta = getBeta(birth, death,fossil);
    //E = exp(-birth*(beta-alpha)*(endTime-startTime));
    //D = exp(getLogProbObs(maxTime-endTime, birth, alpha, beta))*(1.-E);
////	printf("E %lf\t D %lf\n", E, D);
    //return 2*log(beta-alpha)
    //-birth*(beta-alpha)*(endTime-startTime)
    //-log(1.-E)
    //-log(D)
    //+ (k+1)*(log(D)-log(beta-alpha*E-D));
///*	return 2*log(beta-alpha)
    //-birth*(beta-alpha)*(endTime-startTime)
    //+((double)k-1.)*log(1.-exp(-birth*(beta-alpha)*(endTime-startTime)))
    //+((double)k)*getLogProbObs(maxTime-endTime, birth, alpha, beta)
    //-((double)k+1.)*log(beta-alpha*exp(-birth*(beta-alpha)*(endTime-startTime))-(1.-exp(-birth*(beta-alpha)*(endTime-startTime)))*exp(getLogCompProbObs(maxTime-endTime, birth, alpha, beta)));
//*/
///*	if(isinf(((double)k)*getLogProbObs(maxTime-endTime, birth, alpha, beta)-((double)k+1.)*getLogCompProbObs(maxTime-endTime, birth, alpha, beta)-log(1.-exp(-birth*(beta-alpha)*(endTime-startTime)))))
    //printf("%.2lE %.2lE %.2lE %d %.2lE %.2lE %.2lE\n", startTime, endTime, maxTime, k,  birth, death, fossil);
    //return ((double)k)*getLogProbObs(maxTime-endTime, birth, alpha, beta)-((double)k+1.)*getLogCompProbObs(maxTime-endTime, birth, alpha, beta)-log(1.-exp(-birth*(beta-alpha)*(endTime-startTime)));
//*/}

//double getSpecialLogLikelihoodWithoutFossil(double startTime, double endTime, double maxTime, int k,  double birth, double death, double fossil) {
    //double alpha, beta;
    //alpha = getAlpha(birth, death,fossil);
    //beta = getBeta(birth, death,fossil);
    //return 2*log(beta-alpha)
    //-birth*(beta-alpha)*(endTime-startTime)
    //+((double)k-1.)*log(1.-exp(-birth*(beta-alpha)*(endTime-startTime)))
    //-((double)k+1.)*log(beta-alpha*exp(-birth*(beta-alpha)*(endTime-startTime)));
///*	if(isinf(((double)k)*getLogProbObs(maxTime-endTime, birth, alpha, beta)-((double)k+1.)*getLogCompProbObs(maxTime-endTime, birth, alpha, beta)-log(1.-exp(-birth*(beta-alpha)*(endTime-startTime)))))
    //printf("%.2lE %.2lE %.2lE %d %.2lE %.2lE %.2lE\n", startTime, endTime, maxTime, k,  birth, death, fossil);
    //return ((double)k)*getLogProbObs(maxTime-endTime, birth, alpha, beta)-((double)k+1.)*getLogCompProbObs(maxTime-endTime, birth, alpha, beta)-log(1.-exp(-birth*(beta-alpha)*(endTime-startTime)));
//*/}


///*return P_o(k, startTime, endTime)*(P^\star_o)^k*/
//double getLogProbWaitFossilExact(double startTime, double endTime, double maxTime, int k,  double birth, double death, double fossil) {
    //double alpha, beta, delta;
    //alpha = getAlpha(birth, death,fossil);
    //beta = getBeta(birth, death,fossil);
    //delta = birth*(beta-alpha);
    //return
        //log((double)k+1.)
        //+ log(fossil)
        //-delta*(endTime-startTime)
        //+((double)k)*(log(1-exp(-delta*(endTime-startTime))))
        //+((double)k+2.)*log(beta-alpha*exp(-delta*(maxTime-endTime)))
        //-((double)k+2.)*log(beta-alpha*exp(-delta*(maxTime-startTime)))
        //-((double)k)*log(beta-alpha);
//}

///*return p(k, startTime)*/
//double getLogProbWaitFossilMore(double startTime, double maxTime, int k,  double birth, double death, double fossil) {
    //double alpha, beta, delta, time;
    //alpha = getAlpha(birth, death,fossil);
    //beta = getBeta(birth, death,fossil);
    //delta = birth*(beta-alpha);
    //time = maxTime-startTime;
    //return
        //2.*log(beta-alpha)
        //-delta*time
        //+((double)k-1.)*log(1.-exp(-delta*time))
        //-((double)k+1.)*log(beta-alpha*exp(-delta*time));
//}

//double getEventFossilProbLogX(TypeListEvent *event, TypeModelParam *pa) {
    //return getEventFossilProbLog(event, pa->birth,  pa->death, pa->fossil);
//}

//double getEventFossilProbLogY(TypeListEvent *event, TypeModelParam *pa) {
    //return getEventFossilProbLog(event, pa->birth,  pa->death, real_fossil);
//}

//double getEventFossilProbLog(TypeListEvent *event, double birth, double death, double fossil) {
    //double res = 0., tprec = event->minTime;
    //int i;
    //double alpha, beta;
    //alpha = getAlpha(birth, death,fossil);
    //beta = getBeta(birth, death,fossil);
    //for(i=0; i<event->size; i++) {
        //res += getWaitFossilProbLog(event->list[i].n, event->maxTime-tprec, event->list[i].time-tprec,  birth,  death,  fossil, alpha, beta);
        //tprec = event->list[i].time;
        //switch(event->list[i].type) {
            //case 'b':
                //res += getBirthProbLog(event->maxTime-event->list[i].time, birth, death, fossil, alpha, beta);
                //break;
            //case 'd':
                //res += getDeathProbLog(event->maxTime-event->list[i].time, birth, death, fossil, alpha, beta);
                //break;
            //case 'f':
                //res += getFossilProbLog(event->maxTime-event->list[i].time, birth, death, fossil, alpha, beta);
                //break;
        //}
    //}
    //i = event->list[event->size-1].n;
    //if(event->list[event->size-1].type == 'b')
        //i++;
    //if(event->list[event->size-1].type == 'd')
        //i--;
    //return res+getEndFossilProbLog(i, event->maxTime-event->list[event->size-1].time,  birth,  death,  fossil, alpha, beta)+getObsProbLog(event->maxTime,  birth, death, fossil, alpha, beta);
//}

///*complete reconstruction likelihood*/

//double getWaitCompleteProbLog(int n, double t, double x, double birth, double death, double fossil) {
    //return -n*(birth+death+fossil)*x+log(n)+log(birth+death+fossil);
//}

//double getBirthCompleteProbLog(double birth, double death, double fossil) {
    //return log(birth)-log(birth+death+fossil);
//}

//double getDeathCompleteProbLog(double birth, double death, double fossil) {
    //return log(death)-log(birth+death+fossil);
//}

//double getFossilCompleteProbLog(double birth, double death, double fossil) {
    //return log(fossil)-log(birth+death+fossil);
//}

//double getEventCompleteProbLog(TypeListEvent *event, double birth, double death, double fossil) {
    //double res = 0., tprec = 0;
    //int i;
    //for(i=0; i<event->size; i++) {
        //res += getWaitCompleteProbLog(event->list[i].n, event->maxTime-tprec, event->list[i].time-tprec,  birth,  death,  fossil);
        //tprec = event->list[i].time;
        //switch(event->list[i].type) {
            //case 'b':
                //res += getBirthCompleteProbLog(birth, death, fossil);
                //break;
            //case 'd':
                //res += getDeathCompleteProbLog(birth, death, fossil);
                //break;
            //case 'f':
                //res += getFossilCompleteProbLog(birth, death, fossil);
                //break;
        //}
    //}
    //i = event->list[event->size-1].n;
    //if(event->list[event->size-1].type == 'b')
        //i++;
    //if(event->list[event->size-1].type == 'd')
        //i--;
    //return res-i*(birth+death+fossil)*(event->maxTime-event->list[event->size-1].time);
//}

///*contemporary reconstruction likelihood*/

//double getfuncULog(double t, double birth, double death) {
    //return log(birth)+log(1.-exp(-(birth-death)*t))-log(birth-death*exp(-(birth-death)*t));
//}

//double getfuncPLog(double t, double birth, double death) {
    //return log(birth-death)-log(birth-death*exp(-(birth-death)*t));
//}

//double getWaitPhyloProbLog(int n, double t, double x, double birth, double death) {
    //return (log(n)+log(birth-death)-n*(birth-death)*x+(n-1)*log(1.-(death*exp(-(birth-death)*(t-x)))/birth)-n*log(1.-(death*exp(-(birth-death)*t))/birth));
//}

//double getEndPhyloProbLog(int n, double t, double birth, double death) {
    //return n*(-(birth-death)*t+log(1.-death/birth)-log(1.-(death*exp(-(birth-death)*t))/birth));
//}

//double getEventPhyloProbLog(TypeListEvent *event, double birth, double death) {
    //double res = 0., tprec = event->minTime;
    //int i;
    //if(event->size == 0)
        //return log(death)+log(1-exp(-(birth-death)*(event->maxTime-event->minTime)))-log(birth-death*exp(-(birth-death)*(event->maxTime-event->minTime)));
    //for(i=0; i<event->size; i++) {
        //res += getWaitPhyloProbLog(event->list[i].n, event->maxTime-tprec, event->list[i].time-tprec,  birth,  death);
        //tprec = event->list[i].time;
    //}
    //return res+getEndPhyloProbLog(event->list[event->size-1].n+1, event->maxTime-event->list[event->size-1].time,  birth,  death)+getfuncPLog(event->maxTime, birth, death);
//}

///*Complete without fossils reconstruction*/

//double getWaitPureProbLog(int n, double t, double x, double birth, double death) {
    //return -n*(birth+death)*x+log(n)+log(birth+death);
//}

//double getBirthPureProbLog(double birth, double death) {
    //return log(birth)-log(birth+death);
//}

//double getDeathPureProbLog(double birth, double death) {
    //return log(death)-log(birth+death);
//}

//double getEventPureProbLog(TypeListEvent *event, double birth, double death) {
    //double res = 0., tprec = 0;
    //int i;

    //if(event->size == 0)
        //return getWaitPureProbLog(1, 0., event->maxTime,  birth,  death);
    //for(i=0; i<event->size; i++) {
        //if(event->list[i].type != 'f') {
            //res += getWaitPureProbLog(event->list[i].n, event->maxTime-tprec, event->list[i].time-tprec,  birth,  death);
            //tprec = event->list[i].time;
            //switch(event->list[i].type) {
                //case 'b':
                    //res += getBirthPureProbLog(birth, death);
                    //break;
                //case 'd':
                    //res += getDeathPureProbLog(birth, death);
                    //break;
            //}
        //}
    //}
    //i = event->list[event->size-1].n;
    //if(event->list[event->size-1].type == 'b')
        //i++;
    //if(event->list[event->size-1].type == 'd')
        //i--;
    //return res-i*(birth+death)*(event->maxTime-tprec);
//}



///*************/
///*Estimations*/
///*************/

//int estimateML(TypeListEvent *event, double *p) {
    //int b, d, f;
    //double totTime;

    //getStatEvent(event, &b,  &d,  &f);
    //totTime = getTotalTimeEvent(event);
    //if(totTime == 0.) {
        //printf("Problem\n");
        //return 0;
    //}
    //p[0] = ((double)b)/totTime;
    //p[1] = ((double)d)/totTime;
    //p[2] = ((double)f)/totTime;
    //return 1;
//}

//void setFossilRate(double fos) {
    //real_fossil = fos;
//}

//void setSample(double sam) {
    //real_sample = sam;
//}

//double funcLogU(double t, double birth, double death) {
    //return log(birth)+log(1.-exp(-(birth-death)*t))-log(birth-death*exp(-(birth-death)*t));
//}

//double funcLogP(double t, double birth, double death) {
    //return log(birth-death)-log(birth-death*exp(-(birth-death)*t));
//}

//double funcLogCondition(double t, double birth, double death) {
    //return funcLogP(t,birth,death)+((double)minContemp-1.)*funcLogU(t,birth,death);
//}

//double funcFossilGSL(const gsl_vector *x, void *params) {
    //double birth, death, fossil;

    //birth = gsl_vector_get(x, 0);
    //death = gsl_vector_get(x, 1);
    //fossil = gsl_vector_get(x,2);
    //if(death>=birth || birth<0 || death<0 || fossil<0)
        //return RINFTY;
    //return -getEventFossilProbLog((TypeListEvent*)params, birth, death, fossil);
//}

//double funcSamFosFixedGSL(const gsl_vector *x, void *params) {
    //double birth, death;

    //birth = gsl_vector_get(x, 0);
    //death = gsl_vector_get(x, 1);
    //if(death>=birth || birth<0 || death<0)
        //return RINFTY;
    //return -getEventSamFosProbLog((TypeListEvent*)params, birth, death, real_fossil, real_sample);
//}



//double funcFossilFixedGSL(const gsl_vector *x, void *params) {
    //double birth, death;

    //birth = gsl_vector_get(x, 0);
    //death = gsl_vector_get(x, 1);
    //if(death>=birth || birth<0 || death<0)
        //return RINFTY;
    //return -getEventFossilProbLog((TypeListEvent*)params, birth, death, real_fossil);
//}

//double funcFossilFixedNeeGSL(const gsl_vector *x, void *params) {
    //double birth, death;

    //birth = gsl_vector_get(x, 0);
    //death = gsl_vector_get(x, 1);
    //if(death>=birth || birth<0 || death<0)
        //return RINFTY;
    //return -getEventFossilProbLogNee((TypeListEvent*)params, birth, death, real_fossil);
//}

//double funcPhyloGSL(const gsl_vector *x, void *params) {
    //double birth, death;
    //birth = gsl_vector_get(x, 0);
    //death = gsl_vector_get(x, 1);
    //if(death>=birth || birth<0 || death<0)
        //return RINFTY;
    //return -getEventPhyloProbLog((TypeListEvent*)params, birth, death);
//}

//double getEventPhyloProbLogX(TypeListEvent *event, TypeModelParam *pa) {
    //return getEventPhyloProbLog(event, pa->birth, pa->death);
//}

//double funcPhyloNeeGSL(const gsl_vector *x, void *params) {
    //double birth, death;
    //birth = gsl_vector_get(x, 0);
    //death = gsl_vector_get(x, 1);
    //if(death>=birth || birth<0 || death<0)
        //return RINFTY;
    //return -getEventPhyloProbLogNee((TypeListEvent*)params, birth, death);
//}

//double funcCompleteGSL(const gsl_vector *x, void *params) {
    //double birth, death, fossil;
    //birth = gsl_vector_get(x,0);
    //death = gsl_vector_get(x,1);
    //fossil = gsl_vector_get(x,2);
    //if(death>=birth || birth<0 || death<0 || fossil<0)
        //return RINFTY;
    //return -getEventCompleteProbLog((TypeListEvent*)params, birth, death, fossil);
//}

//double funcPureGSL(const gsl_vector *x, void *params) {
    //double birth, death;
    //birth = gsl_vector_get(x,0);
    //death = gsl_vector_get(x,1);
    //if(death>=birth || birth<0 || death<0)
        //return RINFTY;
    //return -getEventPureProbLog((TypeListEvent*)params, birth, death);
//}

///*fossil + sampled reconstruction likelihood*/
//double getRate(double birth, double death, double fossil) {
    //return sqrt(pow(birth+death+fossil, 2.)-4.*birth*death);
//}


//double getWaitSamFosProbLog(int n, double t, double x, double birth, double death, double fossil, double sample, double rate) {
    //return log(n)-log(2)-(((double)n)/2.)*(birth+fossil-death+rate)*x+
    //((double)n-1.)*log(((2*sample-1)*birth+death+fossil+rate)-((2*sample-1)*birth+death+fossil-rate)*exp(-rate*(t-x)))-
    //((double)n)*log(((2*sample-1)*birth+death+fossil+rate)-((2*sample-1)*birth+death+fossil-rate)*exp(-rate*(t)))+
    //log(((2*sample-1.)*birth+death+fossil+rate)*(birth+fossil-death+rate)-((2*sample-1.)*birth+death+fossil-rate)*(birth+fossil-death-rate)*exp(-rate*(t-x)));
//}

//double getEndSamFosProbLog(int n, double t, double birth, double death, double fossil, double sample, double rate) {
    //return -(((double)n)/2.)*(birth+fossil-death+rate)*t+((double)n)*(log(-(((2*sample-1)*birth+death+fossil-rate)-((2*sample-1)*birth+death+fossil+rate)))-log(-(((2*sample-1)*birth+death+fossil-rate)*exp(-rate*t)-((2*sample-1)*birth+death+fossil+rate))));
//}

//double getObsSamFosProb(double t, double birth, double death, double fossil, double sample, double rate) {
    //return exp(getObsSamFosProbLog( t,  birth,  death,  fossil,  sample,  rate));
//}

//double getObsSamFosProbLog(double t, double birth, double death, double fossil, double sample, double rate) {
    //return -log(2)-log(birth)+log(-(((2*sample-1)*birth+death+fossil-rate)*(birth-death-fossil-rate)*exp(-rate*t)-((2*sample-1)*birth+death+fossil+rate)*(birth-death-fossil+rate)))-log(-(((2*sample-1)*birth+death+fossil-rate)*exp(-rate*t)-((2*sample-1)*birth+death+fossil+rate)));
//}

//double getBirthSamFosProbLog(double t, double birth, double death, double fossil, double sample, double rate) {
    //return log(birth)+ getObsSamFosProb(t, birth, death, fossil, sample, rate)-log(birth*getObsSamFosProb(t, birth, death, fossil, sample, rate)+fossil);
//}
//double getDeathSamFosProbLog(double t, double birth, double death, double fossil, double sample, double rate) {
    //return log(fossil)+log(1.-getObsSamFosProb(t, birth, death, fossil, sample, rate))-log(birth*getObsSamFosProb(t, birth, death, fossil, sample, rate)+fossil);
//}
//double getFossilSamFosProbLog(double t, double birth, double death, double fossil, double sample, double rate) {
    //return log(fossil) + getObsSamFosProb(t, birth, death, fossil, sample, rate)-log(birth*getObsSamFosProb(t, birth, death, fossil, sample, rate)+fossil);
//}

//double getEventSamFosProbLog(TypeListEvent *event, double birth, double death, double fossil, double sample) {
    //double res = 0., tprec = 0;
    //int i;
    //double rate;
    //rate = getRate(birth, death, fossil);
////printf("birth %lf death %lf fossil %lf sample %lf rate %lf\n", birth, death, fossil, sample, rate);
    //if(event->size == 0)
        //return 1.-getObsSamFosProbLog(event->maxTime,  birth,  death, fossil, sample, rate);
    //for(i=0; i<event->size; i++) {
////printf("%d\twait %lf\twait %lf\n",i, getWaitSamFosProbLog(event->list[i].n, event->maxTime-tprec, event->list[i].time-tprec,  birth,  death,  fossil, sample, rate), getWaitFossilProbLog(event->list[i].n, event->maxTime-tprec, event->list[i].time-tprec,  birth,  death,  fossil, getAlpha(birth, death, fossil), getBeta(birth, death, fossil)));
////printf("%d\tend %lf\tend %lf\n",i, getEndSamFosProbLog(event->list[i].n, event->maxTime-tprec,  birth,  death,  fossil, sample, rate), getEndFossilProbLog(event->list[i].n, event->maxTime-tprec,  birth,  death,  fossil, getAlpha(birth, death, fossil), getBeta(birth, death, fossil)));
////printf("%d\tt %lf\tobss %lf\tobsf %lf\n",i, event->maxTime-event->list[i].time, getObsSamFosProbLog(event->maxTime-event->list[i].time, birth, death, fossil, sample, rate), getObsProbLog(event->maxTime-event->list[i].time, birth, death, fossil, getAlpha(birth, death, fossil), getBeta(birth, death, fossil)));
        //res += getWaitSamFosProbLog(event->list[i].n, event->maxTime-tprec, event->list[i].time-tprec,  birth,  death,  fossil, sample, rate);
        //tprec = event->list[i].time;
        //switch(event->list[i].type) {
            //case 'b':
                //res += getBirthSamFosProbLog(event->maxTime-event->list[i].time, birth, death, fossil, sample, rate);
                //break;
            //case 'd':
                //res += getDeathSamFosProbLog(event->maxTime-event->list[i].time, birth, death, fossil, sample, rate);
                //break;
            //case 'f':
                //res += getFossilSamFosProbLog(event->maxTime-event->list[i].time, birth, death, fossil, sample, rate);
                //break;
        //}
    //}
    //i = event->list[event->size-1].n;
    //if(event->list[event->size-1].type == 'b')
        //i++;
    //if(event->list[event->size-1].type == 'd')
        //i--;
    //return res+getObsSamFosProbLog(event->maxTime,  birth,  death, fossil, sample, rate)+getEndSamFosProbLog(i, event->maxTime-tprec,  birth,  death, fossil, sample, rate);
//}


//double getEventFossilProbLogNee(TypeListEvent *event, double birth, double death, double fossil) {
    //double res = 0., tprec = 0;
    //int i;
    //double alpha, beta;
    //alpha = getAlpha(birth, death, fossil);
    //beta = getBeta(birth, death, fossil);
    //if(event->size == 0)
        //return 1.-getObsProbLog(event->maxTime,  birth,  death, fossil, alpha, beta);
    //tprec = event->list[0].time;
    //switch(event->list[0].type) {
        //case 'b':
            //res = getBirthProbLog(event->maxTime-event->list[0].time, birth, death, fossil, alpha, beta);
            //break;
        //case 'd':
            //res = getDeathProbLog(event->maxTime-event->list[0].time, birth, death, fossil, alpha, beta);
            //break;
        //case 'f':
            //res = getFossilProbLog(event->maxTime-event->list[0].time, birth, death, fossil, alpha, beta);
            //break;
    //}
    //for(i=1; i<event->size; i++) {
        //res += getWaitFossilProbLog(event->list[i].n, event->maxTime-tprec, event->list[i].time-tprec,  birth,  death,  fossil, alpha, beta);
        //tprec = event->list[i].time;
        //switch(event->list[i].type) {
            //case 'b':
                //res += getBirthProbLog(event->maxTime-event->list[i].time, birth, death, fossil, alpha, beta);
                //break;
            //case 'd':
                //res += getDeathProbLog(event->maxTime-event->list[i].time, birth, death, fossil, alpha, beta);
                //break;
            //case 'f':
                //res += getFossilProbLog(event->maxTime-event->list[i].time, birth, death, fossil, alpha, beta);
                //break;
        //}
    //}
    //i = event->list[event->size-1].n;
    //if(event->list[event->size-1].type == 'b')
        //i++;
    //if(event->list[event->size-1].type == 'd')
        //i--;
    //return res+getEndFossilProbLog(i, event->maxTime-tprec,  birth,  death,  fossil, alpha, beta);
//}

//double getEventPhyloProbLogNee(TypeListEvent *event, double birth, double death) {
    //double res = 0., tprec;
    //int i;
    //if(event->size == 0)
        //return 0.;
    //tprec = event->list[0].time;
    //for(i=1; i<event->size; i++) {
        //res += getWaitPhyloProbLog(event->list[i].n, event->maxTime-tprec, event->list[i].time-tprec,  birth,  death);
////printf("Nee%d/%d %lf %lf\t%lf\n", i, event->size, event->maxTime-tprec, event->list[i].time-tprec, res);
        //tprec = event->list[i].time;
    //}
////printf("Nee %d %lf\t%lf\t%lf\n", event->list[event->size-1].n+1, event->maxTime-tprec, getEndPhyloProbLog(event->list[event->size-1].n+1,  event->maxTime-tprec,  birth,  death), res+getEndPhyloProbLog(event->list[event->size-1].n+1,  event->maxTime-tprec,  birth,  death));
    //return res+getEndPhyloProbLog(event->list[event->size-1].n+1,  event->maxTime-tprec,  birth,  death);
//}
