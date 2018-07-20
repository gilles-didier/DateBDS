#ifndef DrawTreeGenericF
#define DrawTreeGenericF
#include <stdio.h>
#include "Tree.h"

typedef struct POINT {
	double x, y;
} TypePoint;

typedef struct RGB_COLOR {
	double red, green, blue;
} TypeRGB;

typedef struct PARAM_DRAW_TREE_GENERIC {
	double scale, leafSep, leafCur, labelSep, ycenter, ydec, radius, roffset, xoffset, yoffset, height, width, labelWidth, tickLength;
	double tmin, tmax;
	TypeRGB curgb;
	void *info;
} TypeParamDrawTreeGeneric;

typedef struct FUNCT_DRAW_TREE_GENERIC {
	void (*start)(char *, TypeTree *, TypeParamDrawTreeGeneric *), (*end)(TypeParamDrawTreeGeneric *), (*drawText)(double, double, char*, char*, TypeParamDrawTreeGeneric*), 
	(*drawLine)(double, double, double, double, TypeParamDrawTreeGeneric*);
} TypeFunctDrawTreeGeneric;


typedef struct INFO_DRAW_TREE_GENERIC {
	TypeParamDrawTreeGeneric param;
	TypeFunctDrawTreeGeneric funct;
} TypeInfoDrawTreeGeneric;

typedef struct ADDITIONAL_DRAW_TREE_GENERIC {
	void (*draw)(int, double, double, TypeInfoDrawTreeGeneric*, void *), *data;
} TypeAdditionalDrawTreeGeneric;

void fillTime(int n, double tanc, TypeTree *tree, double *min, double *max, int *dmax);
void fillBounds(int n, double tmin, double tmax, TypeTree *tree, double *min, double *max, int *dmax);
void fillUnknownTimes(double tmin, double tmax, TypeTree *tree);

void drawScaleGeneric(double x, double y, TypeInfoDrawTreeGeneric *info);
void drawTreeFileGeneric(char *filename, TypeTree *tree, TypeInfoDrawTreeGeneric *info, TypeAdditionalDrawTreeGeneric *add);
void drawTreeFileGenericDebug(char *filename, TypeTree *tree, TypeInfoDrawTreeGeneric *info, TypeAdditionalDrawTreeGeneric *add);
#endif
