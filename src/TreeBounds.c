#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "Utils.h"
#include "TreeBounds.h"

#define MAX_NAME_SIZE 3000
#define TAG_ORIGIN "ORI"
#define TAG_END "END"

void getBounds(double *min, double *max, char *s) {
    int i, ind;
    char tmp[MAX_NAME_SIZE];
    for(i=0; s[i] != '\0'; i++)
		if(s[i] == ',')
			s[i] = '.';
    for(i=0; s[i] != '\0' && isspace(s[i]); i++);
    if(s[i] == '(' || s[i] == '[') {
        i++;
        for(; s[i] != '\0' && isspace(s[i]); i++);
        ind = 0;
        for(; s[i] != '\0' && s[i] != ':' && s[i] != ';' && !isspace(s[i]); i++)
            tmp[ind++] = s[i];
        tmp[ind] = '\0';
		if(ind>0)
			*min = atof(tmp);
		else
			*min = NO_TIME;
        for(; s[i] != '\0' && isspace(s[i]); i++);
        if(s[i] != ':' && s[i] != ';')
            exitProg(ErrorReading, "Missing ':' or ';' while reading a time interval.");
        i++;
        for(; s[i] != '\0' && isspace(s[i]); i++);
        ind = 0;
        for(; s[i] != '\0' && s[i] != ')' && s[i] != ']' && !isspace(s[i]); i++)
            tmp[ind++] = s[i];
        tmp[ind] = '\0';
		if(ind>0)
			*max =	atof(tmp);
		else
			*max = NO_TIME;
        for(; s[i] != '\0' && isspace(s[i]); i++);
        if(s[i] != ')' && s[i] != ']')
            exitProg(ErrorReading, "Missing ')' or ']' while reading a time interval.");
        if(*min != NO_TIME && *max != NO_TIME && *min > *max) {
            double tmp = *max;
            *max = *min;
            *min = tmp;
        }
    } else {
		*min = NO_TIME;
		*max = NO_TIME;
    }
}


void fillBoundsFromComments(double *tmin, double *tmax, TypeTree *tree) {
	int n;
	for(n=0; n<tree->size; n++) {
		tmin[n] = NO_TIME;
		tmax[n] = NO_TIME;
	}
	if(tree->comment != NULL) {
		char *s;
		for(n=0; n<tree->size; n++) {
			if(tree->node[n].child != NOSUCH && tree->comment[n] != NULL) {
				s = tree->comment[n];
				while(s!=NULL && s[0]!='\0') {
					char *tag, *val;
					s = nextTag(s, &tag, &val);
					if(tag != NULL && strcmp(tag, "BND")==0 && val != NULL)
						getBounds(&(tmin[n]), &(tmax[n]), val);
				}
            }
        }
        if(tree->comment[tree->root] != NULL) {
            s = tree->comment[tree->root];
            while(s!=NULL && s[0]!='\0') {
                char *tag, *val;
                s = nextTag(s, &tag, &val);
                if(tag != NULL && strcmp(tag, TAG_ORIGIN)==0 && val != NULL)
                    tree->minTime = atof(val);
                if(tag != NULL && strcmp(tag, TAG_END)==0 && val != NULL)
                    tree->maxTime = atof(val);
            }
        }
    }
}
