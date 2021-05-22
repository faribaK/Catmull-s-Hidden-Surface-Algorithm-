
/*

Dummy routines for matrix transformations.

These are for you to write!

*/


#include <stdio.h>
#include <math.h>
#include <fstream>
#include <Windows.h>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <iostream>
#include "icVector.H"

//#include "osuGraphics.h"
#include "lines.h"
#include "matlib.h"
using namespace std;
//-------------------------------------------------
#define stackSize 20
int X = 1, Y = 2, Z = 3, W = 0;

int drawPolygon = 0, drawLine = 0, draw2D = 0, draw3D = 0, projectionPers = 2, projectionOrtho = 2, orhtoPers = 1;
int screenPixelBuffer[500][500] = { 0 };
typedef double matrix[4][4];
typedef double vector4[4];
static matrix currentTransMatrix, currentProjectMatrix, currentViewMAtrix;
int stackCounter = -1;
static double matrixStack[stackSize][4][4];
double **zBuffer = NULL;
int depthTest = 0;
int shadeModel = 2;
float ambientIntensity, rDiffuse, gDiffuse, bDiffuse, rSpecular, gSpecular, bSpecular, phongCoefficient;
double globalNormal[3] = { 0 };
int diffuse = 0, ambient = 0;
icVector3 eyePos(0, 0, -1);
struct osuVertex2 {
	double x, y, w = 1;
	double rV, gV, bV;
};

//struct osuVertex3 {
//	double x, y, z, w = 1;
//	double rV, gV, bV;
//};

struct osuColor3 {
	double rC, gC, bC;
};

struct osuColor4 {
	double rC, gC, bC, aC;
};

typedef struct osuVertex {
	float x, y, z;
	float r, g, b;
	double nX, nY, nZ;
} osuVertex;

typedef struct pointLight {
	float pl0, pl1, pl2;
	float pI;
} pointLight;

typedef struct directLight {
	float dl0, dl1, dl2;
	float dI;
} directLight;

typedef struct Light {
	float l0, l1, l2;
	float Intensity;
} Light;

typedef struct Poly {
	vector <osuVertex> vertexList;
	//int vertexCount;
}Poly;

typedef struct pixelBucket {
	Poly piecePoly;
	int isSolid;
}pixelBucket;

static vector <Poly> polyList;
static vector <int> highest_y;
vector <Poly> activePolyList;
vector <pixelBucket> xBucket[360];

static vector <osuVertex2> vertexList;
static vector <osuVertex> vertexList3D;
static osuColor3 currentColorState;

static vector <pointLight> pointLightList;
static vector <directLight> directLightList;

static 

/******************************************************************************
Copy the contents of a source vertex to a destination vertex.
******************************************************************************/

void vertex_copy(osuVertex *dest, osuVertex *source)
{
	dest->x = source->x;
	dest->y = source->y;
	dest->z = source->z;

	dest->r = source->r;
	dest->g = source->g;
	dest->b = source->b;
}


/******************************************************************************
Create a new vertex that is the intersection between a plane and a line
segment between two given vertices.

Entry:
v0,v1 - two vertex endpoints of the line segment
a,b,c,d - coefficients for the plane ax + by + cz + d = 0

Exit:
vnew - the new vertex at the intersection between the line and plane
******************************************************************************/

void create_vertex(
	osuVertex *v0,
	osuVertex *v1,
	osuVertex *newv,
	float a,
	float b,
	float c,
	float d
)
{
	float t;
	float x0, y0, z0;
	float x1, y1, z1;
	float r0, g0, b0;
	float r1, g1, b1;
	float dx, dy, dz;

	/* shorthands */

	x0 = v0->x;
	y0 = v0->y;
	z0 = v0->z;

	x1 = v1->x;
	y1 = v1->y;
	z1 = v1->z;

	dx = x1 - x0;
	dy = y1 - y0;
	dz = z1 - z0;

	/* find parameter t saying how far between v0 and v1 the intersection is */

	t = -1.0 * (a*x0 + b*y0 + c*z0 + d) / (a*dx + b*dy + c*dz);

	/* interpolate between values in v0 and v1 for location and color */

	newv->x = x0 + t * (x1 - x0);
	newv->y = y0 + t * (y1 - y0);
	newv->z = z0 + t * (z1 - z0);

	newv->r = v0->r + t * (v1->r - v0->r);
	newv->g = v0->g + t * (v1->g - v0->g);
	newv->b = v0->b + t * (v1->b - v0->b);
}


/******************************************************************************
Clip a polygon to a plane.

Entry:
verts   - vertices of polygon to clip
count   - number of vertices in polygon
a,b,c,d - coefficients of plane equation against which to clip:
positive side described by ax + by + cz + d > 0 are kept

Exit:
out_verts - vertices of clipped polygon
out_count - number of vertices in the clipped polygon, or 0 if the entire
polygon is on the wrong side of the clipping plane
******************************************************************************/

void poly_clip(
	osuVertex *verts,
	int count,
	osuVertex *out_verts,
	int *out_count,
	float a,
	float b,
	float c,
	float d
)
{
	int i, ii;
	int new_count = 0;
	osuVertex *v0, *v1;
	int in0, in1;  /* are v0 or v1 in the proper half-space */

	v0 = &verts[0];
	in0 = (a * v0->x + b * v0->y + c * v0->z + d > 0);

	for (i = 0; i < count; i++) {

		v0 = &verts[i];
		v1 = &verts[(i + 1) % count];
		in1 = (a * v1->x + b * v1->y + c * v1->z + d > 0);

		if (in0 && in1) {
			vertex_copy(&out_verts[new_count++], v1);
		}
		else if (!in0 && in1) {
			create_vertex(v0, v1, &out_verts[new_count++], a, b, c, d);
			vertex_copy(&out_verts[new_count++], v1);
		}
		else if (in0 && !in1) {
			create_vertex(v0, v1, &out_verts[new_count++], a, b, c, d);
		}
		else {
			/* both are not in, so we add no vertices to the clipped polygon */
		}

		in0 = in1;
	}

	*out_count = new_count;
}


void poly_clip_negative(
	osuVertex *verts,
	int count,
	osuVertex *out_verts,
	int *out_count,
	float a,
	float b,
	float c,
	float d
)
{
	int i, ii;
	int new_count = 0;
	osuVertex *v0, *v1;
	int in0, in1;  /* are v0 or v1 in the proper half-space */
	v0 = &verts[0];
	in0 = (a * v0->x + b * v0->y + c * v0->z + d < 0);
	for (i = 0; i < count; i++) {
		v0 = &verts[i];
		v1 = &verts[(i + 1) % count];
		in1 = (a * v1->x + b * v1->y + c * v1->z + d < 0);

		if (in0 && in1) {
			vertex_copy(&out_verts[new_count++], v1);
		}
		else if (!in0 && in1) {
			create_vertex(v0, v1, &out_verts[new_count++], a, b, c, d);
			vertex_copy(&out_verts[new_count++], v1);
		}
		else if (in0 && !in1) {
			create_vertex(v0, v1, &out_verts[new_count++], a, b, c, d);
		}
		else {
			/* both are not in, so we add no vertices to the clipped polygon */
		}
		in0 = in1;
	}
	*out_count = new_count;
}

void loadIdentity(matrix &identity);
void rightmultMatrix(matrix &mult, matrix A, matrix B);
void matvectMult(matrix A, vector4 &V);
osuColor3 getColorLinearInterpolation(double x, double y, osuVertex2 V1, osuVertex2 V2, osuVertex2 V3);
//int checkinside1(double x, double y, osuVertex2 V1, osuVertex2 V2, osuVertex2 V3);
int checkinside2(double x, double y, osuVertex2 Va, osuVertex2 Vb, osuVertex2 Vc, int excludeCorners);
void getBaryCoord(double x, double y, osuVertex V1, osuVertex V2, osuVertex V3, double &bary_a, double &bary_b);


void rightmultMatrix(matrix &mult, matrix a, matrix b) {
	int i, j;
	matrix m;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			//m[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j] + A[i][3] * B[3][j];
			m[i][j] = a[0][j] * b[i][0] + a[1][j] * b[i][1] + a[2][j] * b[i][2] + a[3][j] * b[i][3];
		}
	}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			mult[i][j] = m[i][j];
		}
	}
}

void loadIdentity(matrix &identity) {
	int i, j;
	//matrix id;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			if (i == j)
				identity[i][j] = 1;
			else
				identity[i][j] = 0;
		}
	}
}
//osuColor3 getColorLinearInterpolation(double x, double y, osuVertex V1, osuVertex V2, osuVertex V3);
//int checkinside2(double x, double y, osuVertex Va, osuVertex Vb, osuVertex Vc);

void osuOrtho(double left, double right, double bottom, double top, double nearp, double farp)
{
	matrix scale = {
		2 / (right - left),  0,                   0,                   0,
		0,               2 / (top - bottom),  0,                   0,
		0,               0,                   2 / (nearp - farp),  0,
		0,               0,                   0,                   1
	};

	/*matrix trans = {
	1,               0,                  0,                   0,
	0,               1,                  0,                   0,
	0,               0,                  1,                   0,
	0,               0,                  0,                   1
	};*/

	matrix trans = {
		1,               0,                  0,                  0,
		0,               1,                  0,                  0,
		0,               0,                  1,                  0,
		-1 * (left + right) / 2,               -1 * (bottom + top) / 2,                  -1 * (nearp + farp) / 2,                    1
	};

	matrix scalemultTrans;  rightmultMatrix(scalemultTrans, scale, trans);
	//copyMatrix(*scalemultTrans, currentProjectMatrix);
	int i, j;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			currentProjectMatrix[i][j] = scalemultTrans[i][j];
		}
	}
	currentProjectMatrix[3][2] = currentProjectMatrix[3][2] * (-1);
	projectionPers = 0;
	projectionOrtho = 1;
	orhtoPers = 1;
	

}

void osuBegin(OSUDrawable mode)
{
	if (mode == OSU_LINES) {
		drawLine = 1;
	}
	else if (mode == OSU_TRIANGLE || mode == OSU_POLYGON)
		drawPolygon = 1;
	vertexList.clear();
	vertexList3D.clear();
	draw2D = 0; draw3D = 0;

}

//void osuEnd1()
//{
//	int w, h;
//	osuVertex2 tempV[3];
//	osuVertex tempV3D[3];
//	osuGetFramebufferSize(&w, &h);
//	double *zValues;
//
//	if (drawLine == 1)
//		drawLine = 0;
//	else if (drawPolygon == 1) {
//		if (draw2D == 1) {
//			int vCount = vertexList.size();
//			if (vCount > 2) {
//				for (int i = 0; i < vCount - 2; i++) {
//					tempV[0] = vertexList[0];
//					tempV[1] = vertexList[i + 1];
//					tempV[2] = vertexList[i + 2];
//
//					tempV[0].x = tempV[0].x * w;
//					tempV[0].y = tempV[0].y * (w / h) * h;
//
//					tempV[2].x = tempV[2].x * w;
//					tempV[2].y = tempV[2].y * (w / h) * h;
//
//					tempV[1].x = tempV[1].x * w;
//					tempV[1].y = tempV[1].y * (w / h) * h;
//					int j;
//					double x_min = 30000, x_max = -1, y_min = 30000, y_max = -1;
//					// finding bounding box of triangle
//					for (j = 0; j < 3; j++) {
//						if (tempV[j].x < x_min) {
//							x_min = tempV[j].x;
//						}
//
//						if (tempV[j].x > x_max) {
//							x_max = tempV[j].x;
//						}
//
//						if (tempV[j].y < y_min) {
//							y_min = tempV[j].y;
//						}
//
//						if (tempV[j].y > y_max) {
//							y_max = tempV[j].y;
//						}
//					}
//					/*x_min = floor(x_min)-1; y_min = floor(y_min);
//					x_max = ceil(x_max) ; y_max = ceil(y_max);*/
//					double posx, posy;
//					//	find intersection points with each pixel row and turn pixels on betn those points 
//
//					for (j = y_min; j <= y_max; j++) { // j = pixel row
//						for (int k = x_min; k <= x_max; k++) {
//							int check = checkinside2(k, j, tempV[0], tempV[1], tempV[2],1);
//							if (check == 1) {
//								osuColor3 tempC;
//								//call getcolor
//								posx = double(k) / w;
//								posy = double(j) / w;
//								if (shadeModel == 1) {
//									tempC.rC = vertexList[0].rV * 255;
//									tempC.gC = vertexList[0].gV * 255;
//									tempC.bC = vertexList[0].bV * 255;
//								}
//								else if (shadeModel == 2) {
//									tempC = getColorLinearInterpolation(posx, posy, vertexList[0], vertexList[i + 1], vertexList[i + 2]);
//								}
//								osuWritePixel(k, j, tempC.rC, tempC.gC, tempC.bC);
//
//								//osuColor3 tempC;
//								////call getcolor
//								//posx = double(k) / w;
//								//posy = double(j) / w;
//								//tempC = getColorLinearInterpolation(posx, posy, vertexList[0], vertexList[i + 1], vertexList[i + 2]);
//								//osuWritePixel(k, j, tempC.rC, tempC.gC, tempC.bC);
//							}
//						} //loop to find 
//
//					}
//				}
//			}
//		}
//		else if (draw3D == 1) {
//			osuVertex *vertsInitial = new osuVertex[vertexList3D.size()];
//
//			osuVertex *verts = new osuVertex[vertexList3D.size()];
//
//			double *zValues = new double[vertexList3D.size()];
//
//			for (int i = 0; i < vertexList3D.size(); i++)
//			{
//				osuVertex tempVert = vertexList3D[i];
//
//				vertsInitial[i] = tempVert;
//
//				vector4 tempV4D = { tempVert.x,tempVert.y,tempVert.z,1.0 };
//
//				matvectMult(currentTransMatrix, tempV4D);
//
//				if (orhtoPers == 2)
//				{
//					zValues[i] = 1 / (tempV4D[2] / tempV4D[3]);
//				}
//
//				if (orhtoPers == 1)
//				{
//					zValues[i] = tempV4D[2] / tempV4D[3];
//				}
//				matvectMult(currentProjectMatrix, tempV4D);
//				double x1, y1, z1;
//				x1 = tempV4D[0] / tempV4D[3];
//				y1 = tempV4D[1] / tempV4D[3];
//				z1 = tempV4D[2] / tempV4D[3];
//				verts[i].x = x1;
//				verts[i].y = y1;
//				verts[i].z = z1;
//				verts[i].r = tempVert.r;
//				verts[i].g = tempVert.g;
//				verts[i].b = tempVert.b;
//			}
//
//			osuVertex *outputVertexNegativeX = new osuVertex[vertexList3D.size() * 2];
//			int* outputVertexCountNegativeX = new int;
//			poly_clip(verts, vertexList3D.size(), outputVertexNegativeX, outputVertexCountNegativeX, 0, 0, 1, 1);
//
//			osuVertex *outputVertexPositiveX = new osuVertex[*outputVertexCountNegativeX * 2];
//			int* outputVertexCountPositiveX = new int;
//			poly_clip_negative(outputVertexNegativeX, *outputVertexCountNegativeX, outputVertexPositiveX, outputVertexCountPositiveX, 1, 0, 0, -1);
//
//			osuVertex *outputVertexNegativeY = new osuVertex[*outputVertexCountPositiveX * 2];
//			int* outputVertexCountNegativeY = new int;
//			poly_clip(outputVertexPositiveX, *outputVertexCountPositiveX, outputVertexNegativeY, outputVertexCountNegativeY, 0, 1, 0, 1);
//
//			osuVertex *outputVertexPositiveY = new osuVertex[*outputVertexCountNegativeY * 2];
//			int* outputVertexCountPositiveY = new int;
//			poly_clip_negative(outputVertexNegativeY, *outputVertexCountNegativeY, outputVertexPositiveY, outputVertexCountPositiveY, 0, 1, 0, -1);
//
//			osuVertex *outputVertexNegativeZ = new osuVertex[*outputVertexCountPositiveY * 2];
//			int* outputVertexCountNegativeZ = new int;
//			poly_clip(outputVertexPositiveY, *outputVertexCountPositiveY, outputVertexNegativeZ, outputVertexCountNegativeZ, 0, 0, 1, 1);
//
//			osuVertex *outputVertexPositiveZ = new osuVertex[*outputVertexCountNegativeZ * 2];
//			int* outputVertexCountPositiveZ = new int;
//			poly_clip(outputVertexNegativeZ, *outputVertexCountNegativeZ, outputVertexPositiveZ, outputVertexCountPositiveZ, 0, 0, 1, -1);
//
//			int vertexCount = vertexList3D.size();
//			for (int j = 0; j < vertexCount - 2; j++)
//			{
//				osuVertex vertex1 = verts[0];
//				osuVertex vertex2 = verts[j + 1];
//				osuVertex vertex3 = verts[j + 2];
//
//				osuVertex vertex1ForColor = vertexList3D[0];
//				osuVertex vertex2ForColor = vertexList3D[j + 1];
//				osuVertex vertex3ForColor = vertexList3D[j + 2];
//
//				double z1 = zValues[0];
//
//				double z2 = zValues[j + 1];
//
//				double z3 = zValues[j + 2];
//
//				vertex1.x = vertex1.x * w / 2 + w / 2;
//				vertex1.y = vertex1.y *h / 2 + h / 2;
//
//				vertex2.x = vertex2.x *w / 2 + w / 2;
//				vertex2.y = vertex2.y *h / 2 + h / 2;
//
//				vertex3.x = vertex3.x *w / 2 + w / 2;
//				vertex3.y = vertex3.y *h / 2 + h / 2;
//
//				double xValues[3] = { vertex1.x,vertex2.x,vertex3.x };
//				double yValues[3] = { vertex1.y,vertex2.y,vertex3.y };
//
//				double x_min = 30000, x_max = -1, y_min = 30000, y_max = -1;
//				// finding bounding box of triangle
//				for (int j = 0; j < 3; j++) {
//					if (xValues[j] < x_min) {
//						x_min = xValues[j];
//					}
//
//					if (xValues[j] > x_max) {
//						x_max = xValues[j];
//					}
//
//					if (yValues[j] < y_min) {
//						y_min = yValues[j];
//					}
//
//					if (yValues[j] > y_max) {
//						y_max = yValues[j];
//					}
//				}
//
//				osuVertex2 vert1, vert2, vert3;
//
//				vert1.x = vertex1.x;
//				vert2.x = vertex2.x;
//				vert3.x = vertex3.x;
//				vert1.y = vertex1.y;
//				vert2.y = vertex2.y;
//				vert3.y = vertex3.y;
//
//				for (int pp = x_min; pp < x_max; pp++)
//				{
//					for (int qq = y_min; qq < y_max; qq++)
//					{
//						int check = checkinside2(pp, qq, vert1, vert2, vert3,1);
//						if (check == 1)
//						{
//							double bary1, bary2;
//							getBaryCoord(pp, qq, vertex1, vertex2, vertex3, bary1, bary2);
//							double zCurrent = 0;
//							if (depthTest == 1)
//							{
//								zCurrent = bary1*z1 + bary2*z2 + (1 - bary1 - bary2)*z3;
//								double 	zPrev = zBuffer[pp][qq];
//
//								if (zCurrent > zPrev)
//								{
//									osuColor3 colorCurrPoint;
//									colorCurrPoint.rC = bary1*vertex1.r + bary2*vertex2.r + (1 - bary1 - bary2)*vertex3.r;
//									colorCurrPoint.gC = bary1*vertex1.g + bary2*vertex2.g + (1 - bary1 - bary2)*vertex3.g;
//									colorCurrPoint.bC = bary1*vertex1.b + bary2*vertex2.b + (1 - bary1 - bary2)*vertex3.b;
//
//
//									if (shadeModel == 2)
//									{
//										colorCurrPoint.rC = colorCurrPoint.rC * 255;
//										colorCurrPoint.gC = colorCurrPoint.gC * 255;
//										colorCurrPoint.bC = colorCurrPoint.bC * 255;
//
//									}
//									else if (shadeModel == 1)
//									{
//										colorCurrPoint.rC = vertexList3D[0].r * 255;
//										colorCurrPoint.gC = vertexList3D[0].g * 255;
//										colorCurrPoint.bC = vertexList3D[0].b * 255;
//									}
//
//									if (pp < w && qq < h)
//										osuWritePixel(pp, qq, colorCurrPoint.rC, colorCurrPoint.gC, colorCurrPoint.bC);
//									if (1 == depthTest)
//									{
//										zBuffer[pp][qq] = zCurrent;
//									}
//								}
//							}
//							else {
//								osuColor3 colorCurrPoint;
//								colorCurrPoint.rC = bary1*vertex1.r + bary2*vertex2.r + (1 - bary1 - bary2)*vertex3.r;
//								colorCurrPoint.gC = bary1*vertex1.g + bary2*vertex2.g + (1 - bary1 - bary2)*vertex3.g;
//								colorCurrPoint.bC = bary1*vertex1.b + bary2*vertex2.b + (1 - bary1 - bary2)*vertex3.b;
//
//
//								if (shadeModel == 2)
//								{
//									colorCurrPoint.rC = colorCurrPoint.rC * 255;
//									colorCurrPoint.gC = colorCurrPoint.gC * 255;
//									colorCurrPoint.bC = colorCurrPoint.bC * 255;
//
//								}
//								else if (shadeModel == 1)
//								{
//									colorCurrPoint.rC = vertexList3D[0].r * 255;
//									colorCurrPoint.gC = vertexList3D[0].g * 255;
//									colorCurrPoint.bC = vertexList3D[0].b * 255;
//								}
//
//								if (pp < w && qq < h)
//									osuWritePixel(pp, qq, colorCurrPoint.rC, colorCurrPoint.gC, colorCurrPoint.bC);
//							}
//						}
//					}
//				}
//			}
//			vertexList3D.clear();
//			drawPolygon = 0;
//		}
//	}
//}

void osuEnd()
{
	int w, h;
	osuVertex2 tempV[3];
	osuVertex tempV3D[3];
	osuGetFramebufferSize(&w, &h);
	double *zValues;

	if (drawLine == 1)
		drawLine = 0;
	else if (drawPolygon == 1) {
		if (draw2D == 1) {
			int vCount = vertexList.size();
			if (vCount > 2) {
				for (int i = 0; i < vCount - 2; i++) {
					tempV[0] = vertexList[0];
					tempV[1] = vertexList[i + 1];
					tempV[2] = vertexList[i + 2];

					tempV[0].x = tempV[0].x * w;
					tempV[0].y = tempV[0].y * (w / h) * h;

					tempV[2].x = tempV[2].x * w;
					tempV[2].y = tempV[2].y * (w / h) * h;

					tempV[1].x = tempV[1].x * w;
					tempV[1].y = tempV[1].y * (w / h) * h;
					int j;
					double x_min = 30000, x_max = -1, y_min = 30000, y_max = -1;
					// finding bounding box of triangle
					for (j = 0; j < 3; j++) {
						if (tempV[j].x < x_min) {
							x_min = tempV[j].x;
						}

						if (tempV[j].x > x_max) {
							x_max = tempV[j].x;
						}

						if (tempV[j].y < y_min) {
							y_min = tempV[j].y;
						}

						if (tempV[j].y > y_max) {
							y_max = tempV[j].y;
						}
					}
					/*x_min = floor(x_min)-1; y_min = floor(y_min);
					x_max = ceil(x_max) ; y_max = ceil(y_max);*/
					double posx, posy;
					//	find intersection points with each pixel row and turn pixels on betn those points 

					for (j = y_min; j <= y_max; j++) { // j = pixel row
						for (int k = x_min; k <= x_max; k++) {
							int check = checkinside2(k, j, tempV[0], tempV[1], tempV[2],1);
							if (check == 1) {
								osuColor3 tempC;
								//call getcolor
								posx = double(k) / w;
								posy = double(j) / w;
								if (shadeModel == 1) {
									tempC.rC = vertexList[0].rV * 255;
									tempC.gC = vertexList[0].gV * 255;
									tempC.bC = vertexList[0].bV * 255;
								}
								else if (shadeModel == 2) {
									tempC = getColorLinearInterpolation(posx, posy, vertexList[0], vertexList[i + 1], vertexList[i + 2]);
								}
								osuWritePixel(k, j, tempC.rC, tempC.gC, tempC.bC);

								//osuColor3 tempC;
								////call getcolor
								//posx = double(k) / w;
								//posy = double(j) / w;
								//tempC = getColorLinearInterpolation(posx, posy, vertexList[0], vertexList[i + 1], vertexList[i + 2]);
								//osuWritePixel(k, j, tempC.rC, tempC.gC, tempC.bC);
							}
						} //loop to find 

					}
				}
			}
		}
		else if (draw3D == 1) {
			osuVertex *vertsInitial = new osuVertex[vertexList3D.size()];

			osuVertex *verts = new osuVertex[vertexList3D.size()];

			double *zValues = new double[vertexList3D.size()];

			Poly temp;
			//temp.vertexList = new osuVertex[vertexList3D.size()];
			//temp.vertexCount = vertexList3D.size();
			double y_max = -1000000;
			for (int i = 0; i < vertexList3D.size(); i++)
			{
				osuVertex tempVert = vertexList3D[i];

				vertsInitial[i] = tempVert;

				vector4 tempV4D = { tempVert.x,tempVert.y,tempVert.z,1.0 };

				matvectMult(currentTransMatrix, tempV4D);

				if (orhtoPers == 2)
				{
					zValues[i] = 1 / (tempV4D[2] / tempV4D[3]);
				}

				if (orhtoPers == 1)
				{
					zValues[i] = tempV4D[2] / tempV4D[3];
				}
				double z1 = zValues[i];
				matvectMult(currentProjectMatrix, tempV4D);
				double x1, y1;
				x1 = tempV4D[0] / tempV4D[3];
				y1 = tempV4D[1] / tempV4D[3];
				x1 = x1 * w / 2 + w / 2;
				y1 = y1 *h / 2 + h / 2;
				//z1 = tempV4D[2] / tempV4D[3];
				verts[i].x = x1;
				verts[i].y = y1;
				verts[i].z = z1;
				verts[i].r = tempVert.r;
				verts[i].g = tempVert.g;
				verts[i].b = tempVert.b;
				/*temp.vertexList[i].x = x1;
				temp.vertexList[i].y = y1;
				temp.vertexList[i].z = z1;
				temp.vertexList[i].r = tempVert.r;
				temp.vertexList[i].g = tempVert.g;
				temp.vertexList[i].b = tempVert.b;*/
				temp.vertexList.push_back(verts[i]);
				if (temp.vertexList[i].y > y_max)
					y_max = temp.vertexList[i].y;
			}
			polyList.push_back(temp);
			highest_y.push_back(y_max);

			osuVertex *outputVertexNegativeX = new osuVertex[vertexList3D.size() * 2];
			int* outputVertexCountNegativeX = new int;
			poly_clip(verts, vertexList3D.size(), outputVertexNegativeX, outputVertexCountNegativeX, 0, 0, 1, 1);

			osuVertex *outputVertexPositiveX = new osuVertex[*outputVertexCountNegativeX * 2];
			int* outputVertexCountPositiveX = new int;
			poly_clip_negative(outputVertexNegativeX, *outputVertexCountNegativeX, outputVertexPositiveX, outputVertexCountPositiveX, 1, 0, 0, -1);

			osuVertex *outputVertexNegativeY = new osuVertex[*outputVertexCountPositiveX * 2];
			int* outputVertexCountNegativeY = new int;
			poly_clip(outputVertexPositiveX, *outputVertexCountPositiveX, outputVertexNegativeY, outputVertexCountNegativeY, 0, 1, 0, 1);

			osuVertex *outputVertexPositiveY = new osuVertex[*outputVertexCountNegativeY * 2];
			int* outputVertexCountPositiveY = new int;
			poly_clip_negative(outputVertexNegativeY, *outputVertexCountNegativeY, outputVertexPositiveY, outputVertexCountPositiveY, 0, 1, 0, -1);

			osuVertex *outputVertexNegativeZ = new osuVertex[*outputVertexCountPositiveY * 2];
			int* outputVertexCountNegativeZ = new int;
			poly_clip(outputVertexPositiveY, *outputVertexCountPositiveY, outputVertexNegativeZ, outputVertexCountNegativeZ, 0, 0, 1, 1);

			osuVertex *outputVertexPositiveZ = new osuVertex[*outputVertexCountNegativeZ * 2];
			int* outputVertexCountPositiveZ = new int;
			poly_clip(outputVertexNegativeZ, *outputVertexCountNegativeZ, outputVertexPositiveZ, outputVertexCountPositiveZ, 0, 0, 1, -1);

			//int vertexCount = vertexList3D.size();

			//for (int j = 0; j < vertexCount - 2; j++)
			//{
			//	osuVertex vertex1 = verts[0];
			//	osuVertex vertex2 = verts[j + 1];
			//	osuVertex vertex3 = verts[j + 2];

			//	osuVertex vertex1ForColor = vertexList3D[0];
			//	osuVertex vertex2ForColor = vertexList3D[j + 1];
			//	osuVertex vertex3ForColor = vertexList3D[j + 2];

			//	double z1 = zValues[0];

			//	double z2 = zValues[j + 1];

			//	double z3 = zValues[j + 2];

			//	vertex1.x = vertex1.x * w / 2 + w / 2;
			//	vertex1.y = vertex1.y *h / 2 + h / 2;

			//	vertex2.x = vertex2.x *w / 2 + w / 2;
			//	vertex2.y = vertex2.y *h / 2 + h / 2;

			//	vertex3.x = vertex3.x *w / 2 + w / 2;
			//	vertex3.y = vertex3.y *h / 2 + h / 2;

			//	double xValues[3] = { vertex1.x,vertex2.x,vertex3.x };
			//	double yValues[3] = { vertex1.y,vertex2.y,vertex3.y };

			//	double x_min = 30000, x_max = -1, y_min = 30000, y_max = -1;
			//	// finding bounding box of triangle
			//	for (int j = 0; j < 3; j++) {
			//		if (xValues[j] < x_min) {
			//			x_min = xValues[j];
			//		}

			//		if (xValues[j] > x_max) {
			//			x_max = xValues[j];
			//		}

			//		if (yValues[j] < y_min) {
			//			y_min = yValues[j];
			//		}

			//		if (yValues[j] > y_max) {
			//			y_max = yValues[j];
			//		}
			//	}

			//	osuVertex2 vert1, vert2, vert3;

			//	vert1.x = vertex1.x;
			//	vert2.x = vertex2.x;
			//	vert3.x = vertex3.x;
			//	vert1.y = vertex1.y;
			//	vert2.y = vertex2.y;
			//	vert3.y = vertex3.y;

			//	/*for (int pp = x_min; pp < x_max; pp++)
			//	{
			//		for (int qq = y_min; qq < y_max; qq++)
			//		{
			//			int check = checkinside2(pp, qq, vert1, vert2, vert3);
			//			if (check == 1)
			//			{
			//				double bary1, bary2;
			//				getBaryCoord(pp, qq, vertex1, vertex2, vertex3, bary1, bary2);
			//				double zCurrent = 0;
			//				if (depthTest == 1)
			//				{
			//					zCurrent = bary1*z1 + bary2*z2 + (1 - bary1 - bary2)*z3;
			//					double 	zPrev = zBuffer[pp][qq];

			//					if (zCurrent > zPrev)
			//					{
			//						osuColor3 colorCurrPoint;
			//						colorCurrPoint.rC = bary1*vertex1.r + bary2*vertex2.r + (1 - bary1 - bary2)*vertex3.r;
			//						colorCurrPoint.gC = bary1*vertex1.g + bary2*vertex2.g + (1 - bary1 - bary2)*vertex3.g;
			//						colorCurrPoint.bC = bary1*vertex1.b + bary2*vertex2.b + (1 - bary1 - bary2)*vertex3.b;


			//						if (shadeModel == 2)
			//						{
			//							colorCurrPoint.rC = colorCurrPoint.rC * 255;
			//							colorCurrPoint.gC = colorCurrPoint.gC * 255;
			//							colorCurrPoint.bC = colorCurrPoint.bC * 255;

			//						}
			//						else if (shadeModel == 1)
			//						{
			//							colorCurrPoint.rC = vertexList3D[0].r * 255;
			//							colorCurrPoint.gC = vertexList3D[0].g * 255;
			//							colorCurrPoint.bC = vertexList3D[0].b * 255;
			//						}

			//						if (pp < w && qq < h)
			//							osuWritePixel(pp, qq, colorCurrPoint.rC, colorCurrPoint.gC, colorCurrPoint.bC);
			//						if (1 == depthTest)
			//						{
			//							zBuffer[pp][qq] = zCurrent;
			//						}
			//					}
			//				}
			//				else {
			//					osuColor3 colorCurrPoint;
			//					colorCurrPoint.rC = bary1*vertex1.r + bary2*vertex2.r + (1 - bary1 - bary2)*vertex3.r;
			//					colorCurrPoint.gC = bary1*vertex1.g + bary2*vertex2.g + (1 - bary1 - bary2)*vertex3.g;
			//					colorCurrPoint.bC = bary1*vertex1.b + bary2*vertex2.b + (1 - bary1 - bary2)*vertex3.b;


			//					if (shadeModel == 2)
			//					{
			//						colorCurrPoint.rC = colorCurrPoint.rC * 255;
			//						colorCurrPoint.gC = colorCurrPoint.gC * 255;
			//						colorCurrPoint.bC = colorCurrPoint.bC * 255;

			//					}
			//					else if (shadeModel == 1)
			//					{
			//						colorCurrPoint.rC = vertexList3D[0].r * 255;
			//						colorCurrPoint.gC = vertexList3D[0].g * 255;
			//						colorCurrPoint.bC = vertexList3D[0].b * 255;
			//					}

			//					if (pp < w && qq < h)
			//						osuWritePixel(pp, qq, colorCurrPoint.rC, colorCurrPoint.gC, colorCurrPoint.bC);
			//				}
			//			}
			//		}
			//	}*/
			//}
			vertexList3D.clear();
			drawPolygon = 0;
		}
	}
}

int checkOnScanLine(int s, Poly P) {

	int vertexCount = P.vertexList.size();

	for (int j = 0; j < vertexCount - 2; j++)
	{
		osuVertex vertex1 = P.vertexList[0];
		osuVertex vertex2 = P.vertexList[j + 1];
		osuVertex vertex3 = P.vertexList[j + 2];

		double xValues[3] = { vertex1.x,vertex2.x,vertex3.x };
		double yValues[3] = { vertex1.y,vertex2.y,vertex3.y };

		double x_min = 30000, x_max = -1, y_min = 30000, y_max = -1;
		// finding bounding box of triangle
		for (int j = 0; j < 3; j++) {
			if (xValues[j] < x_min) {
				x_min = xValues[j];
			}

			if (xValues[j] > x_max) {
				x_max = xValues[j];
			}

			if (yValues[j] < y_min) {
				y_min = yValues[j];
			}

			if (yValues[j] > y_max) {
				y_max = yValues[j];
			}
		}

		osuVertex2 vert1, vert2, vert3;

		vert1.x = vertex1.x;
		vert2.x = vertex2.x;
		vert3.x = vertex3.x;
		vert1.y = vertex1.y;
		vert2.y = vertex2.y;
		vert3.y = vertex3.y;
		int qq = s;
		for (int pp = x_min; pp <= x_max; pp++)
		{
			/*if (pp == 342)
				cout << "here;";*/
			int check = checkinside2(pp, qq, vert1, vert2, vert3, 0);
			if (check == 1)
			{
				return 1;
			}
		}
		qq = s+1;
		for (int pp = x_min; pp <= x_max; pp++)
		{
			int check = checkinside2(pp, qq, vert1, vert2, vert3, 1);
			if (check == 1)
			{
				return 1;
			}
		}
	}
	return 0;
}

void clipoff(Poly poly, float a, float b, float c, Poly &A, Poly &B) {
	//A left lower
	//B right upper
		//II.
		osuVertex prev;
		int n = poly.vertexList.size() - 1;
		prev = poly.vertexList[n];
		float dprev;
		dprev = a*prev.x + b*prev.y + c;
		
		//III.
		for (int i = 0; i <= n; i++) {
			osuVertex curr = poly.vertexList[i];
			float dcurr;
			dcurr = a*curr.x + b*curr.y + c;
			//A.
			if (dcurr < 0) {//below
				if (dprev < 0) {//below
					A.vertexList.push_back(curr);
				}
				else { //>=0 on or above
					//v2 prev
					float a1 = prev.y - curr.y;
					float b1 = curr.x - prev.x;
					float c1 = curr.y*prev.x - curr.x*prev.y;
					float cp = a*b1 - a1*b;
					if (cp != 0) {
						float ap = b*c1 - b1*c;
						float bp = a1*c - a*c1;
						osuVertex ip;
						ip.x = ap / cp;
						ip.y = bp / cp;
						double alpha = sqrt((ip.x - curr.x) *(ip.x - curr.x) + (ip.y - curr.y) *(ip.y - curr.y)) /
							sqrt((prev.x - curr.x) *(prev.x - curr.x) + (prev.y - curr.y) *(prev.y - curr.y));
						if (alpha >= 0 && alpha <= 1)
						{
							ip.z = alpha*curr.z + (1 - alpha)*prev.z;
							ip.r = alpha*curr.r + (1 - alpha)*prev.r;
							ip.g = alpha*curr.g + (1 - alpha)*prev.g;
							ip.b = alpha*curr.b + (1 - alpha)*prev.b;
							A.vertexList.push_back(ip);
							B.vertexList.push_back(ip);
							
						}
						A.vertexList.push_back(curr);
					}
				}
			}
			//B.
			else if (dcurr >= 0) {
				if (dprev >= 0) {//below
					B.vertexList.push_back(curr);
				}
				else { //>=0 on or above
					   //v2 prev
					float a1 = prev.y - curr.y;
					float b1 = curr.x - prev.x;
					float c1 = curr.y*prev.x - curr.x*prev.y;
					float cp = a*b1 - a1*b;
					if (cp != 0) {
						float ap = b*c1 - b1*c;
						float bp = a1*c - a*c1;
						osuVertex ip;
						ip.x = ap / cp;
						ip.y = bp / cp;
						double alpha = sqrt((ip.x - curr.x) *(ip.x - curr.x) + (ip.y - curr.y) *(ip.y - curr.y)) /
							sqrt((prev.x - curr.x) *(prev.x - curr.x) + (prev.y - curr.y) *(prev.y - curr.y));
						if (alpha >= 0 && alpha <= 1)
						{
							ip.z = alpha*curr.z + (1 - alpha)*prev.z;
							ip.r = alpha*curr.r + (1 - alpha)*prev.r;
							ip.g = alpha*curr.g + (1 - alpha)*prev.g;
							ip.b = alpha*curr.b + (1 - alpha)*prev.b;
							A.vertexList.push_back(ip);
							B.vertexList.push_back(ip);
							
						}
						B.vertexList.push_back(curr);
					}
				}
			}
			//C.
			prev = curr;
			dprev = dcurr;
		}
}

double findArea(Poly P) {
	double area = 0;
	int vertexCount = P.vertexList.size();
	for (int j = 0; j < vertexCount - 2; j++)
	{
		osuVertex v1 = P.vertexList[0];
		osuVertex v2 = P.vertexList[j + 1];
		osuVertex v3 = P.vertexList[j + 2];

		double a = sqrt((v1.x - v2.x) *(v1.x - v2.x) + (v1.y - v2.y) *(v1.y - v2.y));
		double b = sqrt((v3.x - v2.x) *(v3.x - v2.x) + (v3.y - v2.y) *(v3.y - v2.y));
		double c = sqrt((v1.x - v3.x) *(v1.x - v3.x) + (v1.y - v3.y) *(v1.y - v3.y));

		double p = (a + b + c) / 2;
		area += sqrt(p*(p-a)*(p-b)*(p-c));
	}
	return abs(area);
}

void drawScene() {

	int w, h;
	osuGetFramebufferSize(&w, &h);
	// I. sort polygonList in y
	for (int i = 0; i < highest_y.size(); i++) {
		for (int j = i + 1; j < highest_y.size(); j++) {
			if (highest_y[i] > highest_y[j])
			{
				int temp;
				Poly tempPoly;
				temp = highest_y[i];
				highest_y[i] = highest_y[j];
				highest_y[j] = temp;

				tempPoly = polyList[i];
				polyList[i] = polyList[j];
				polyList[j] = tempPoly;
			}
		}
	}
	
	Poly blackPoly;
	osuVertex blV;
	blV.x = 0; blV.y = 0; blV.z = 1000000;
	blV.r = 0; blV.g = 0; blV.b = 0;
	blackPoly.vertexList.push_back(blV);

	pixelBucket blackPixel;
	blackPixel.isSolid = 1;
	blackPixel.piecePoly = blackPoly;

	//III. for each scanline s
	for (int s = 0; s < h; s++) {
		// II. initialize activePolyList
		activePolyList.clear();
		/*if (s == 198)
			cout << "here\n";
		cout << "\n";*/
		// A. add active polys on s in activePolyList
		for (int i = 0; i < polyList.size(); i++) {
			if (highest_y[i] < s) //s is above i-th poly
				continue;
			int check = checkOnScanLine(s, polyList[i]);
			if (check == 1)
				activePolyList.push_back(polyList[i]);
		}

		// B. initialize x bucket to be empty and ??? scanline to background ???
		for (int i = 0; i < w; i++) {
			xBucket[i].clear();
			//xBucket[i].push_back(blackPixel);
		}
			
		//C.
		for (int i = 0; i < activePolyList.size(); i++) {

			Poly A, B, A1, B1;
			clipoff(activePolyList[i], 0, 1, -s, A, B);
			clipoff(B, 0, 1, -s - 1, A1, B1);
			activePolyList[i] = A1;
			float min_x = 100000, max_x = -100000;
			for (int j = 0; j < activePolyList[i].vertexList.size(); j++) {
				osuVertex v = activePolyList[i].vertexList[j];
				if (v.x < min_x)
					min_x = v.x;
				if (v.x > max_x)
					max_x = v.x;
			}
			int minfloor = floor(min_x);
			int maxfloor = floor(max_x);
			if (maxfloor == max_x)
				maxfloor--;
			int maxceil = ceil(max_x);
			pixelBucket tempPixel;
			if (maxceil - minfloor > 2) {
				/*for (int j = 0; j < activePolyList[i].vertexList.size(); j++) {*/
					A.vertexList.clear();
					B.vertexList.clear();
					clipoff(activePolyList[i], 1, 0, -(floor(min_x) + 1), A, B);
					tempPixel.piecePoly = A;
					double area = findArea(A);
					if(area==1)
						tempPixel.isSolid = 1;
					else
						tempPixel.isSolid = 0;
					if (A.vertexList.size() > 0)
						xBucket[minfloor].push_back(tempPixel);

					A1.vertexList.clear();
					B1.vertexList.clear();
					clipoff(B, 1, 0, -(ceil(max_x) - 1), A1, B1);
					tempPixel.piecePoly = B1;
					area = findArea(B1);
					if (area == 1)
						tempPixel.isSolid = 1;
					else
						tempPixel.isSolid = 0;
					if (B1.vertexList.size() > 0)
						xBucket[maxfloor].push_back(tempPixel);

					for (int k = minfloor + 1; k <= maxfloor - 1; k++) {
						tempPixel.piecePoly = activePolyList[i];
						tempPixel.isSolid = 1;
						xBucket[k].push_back(tempPixel);
					}
				//}
			}
			else if (maxceil - minfloor == 2) {
				A.vertexList.clear();
				B.vertexList.clear();
				clipoff(activePolyList[i], 1, 0, -ceil(min_x), A, B);
				/*xBucket[minfloor].piecePoly = A;
				xBucket[minfloor].isSolid = 0;*/
				tempPixel.piecePoly = A;
				double area = findArea(A);
				if (area == 1)
					tempPixel.isSolid = 1;
				else
					tempPixel.isSolid = 0;
				if (A.vertexList.size() > 0) {
					xBucket[minfloor].push_back(tempPixel);
				}
				
				A1.vertexList.clear();
				B1.vertexList.clear();
				clipoff(B, 1, 0, -floor(max_x), A1, B1);
				/*xBucket[maxfloor].piecePoly = B1;
				xBucket[maxfloor].isSolid = 0;*/
				tempPixel.piecePoly = B1;
				area = findArea(B1);
				if (area == 1)
					tempPixel.isSolid = 1;
				else
					tempPixel.isSolid = 0;
				if (B1.vertexList.size() > 0)
					xBucket[maxfloor].push_back(tempPixel);
			}
			else if (maxceil - minfloor > 0) {

				/*if (max_x == min_x)
					cout << "equal x\n";*/
				/*xBucket[minfloor].piecePoly = activePolyList[i];
				xBucket[minfloor].isSolid = 0;*/
				tempPixel.piecePoly = activePolyList[i];
				double area = findArea(activePolyList[i]);
				if (area == 1)
					tempPixel.isSolid = 1;
				else
					tempPixel.isSolid = 0;
				xBucket[minfloor].push_back(tempPixel);
			}
			else {

			}
		}

		//D.
		osuColor3 intensity;
		//E.
		for (int p = 0; p < w; p++) {
			int polyCount = xBucket[p].size();
			if (polyCount > 1) {
				pixelBucket tempPoly;// = new pixelBucket[polyCount];
				//1.
				for (int g = 0; g < polyCount - 1; g++) {
					for (int h1 = g + 1; h1 <= polyCount - 1; h1++) {
						if (xBucket[p][g].piecePoly.vertexList[0].z > xBucket[p][h1].piecePoly.vertexList[0].z) {
							tempPoly = xBucket[p][h1];
							xBucket[p][h1] = xBucket[p][g];
							xBucket[p][g] = tempPoly;
						}
					}
				}
				//2.a
				if (xBucket[p][0].isSolid == 1) {

					double area1 = findArea(xBucket[p][0].piecePoly);
					/*if (area1 < 1)
						cout << "solid but less\n";*/
					intensity.rC = xBucket[p][0].piecePoly.vertexList[0].r;
					intensity.gC = xBucket[p][0].piecePoly.vertexList[0].g;
					intensity.bC = xBucket[p][0].piecePoly.vertexList[0].b;
					osuWritePixel(p, s, intensity.rC * 255, intensity.gC * 255, intensity.bC * 255);
				}
				else {
					if (xBucket[p].size() > 1) {
						
						double area1 = findArea(xBucket[p][0].piecePoly);
						if (xBucket[p][1].isSolid == 1) {
							double area2 = findArea(xBucket[p][1].piecePoly);
							double w1 = area1;// / (area1 + area2);
							double w2 = (1 - w1);
							intensity.rC = xBucket[p][0].piecePoly.vertexList[0].r*w1 + xBucket[p][1].piecePoly.vertexList[0].r*w2;
							intensity.gC = xBucket[p][0].piecePoly.vertexList[0].g*w1 + xBucket[p][1].piecePoly.vertexList[0].g*w2;
							intensity.bC = xBucket[p][0].piecePoly.vertexList[0].b*w1 + xBucket[p][1].piecePoly.vertexList[0].b*w2;
						}
						else {
							double w1 = area1;// / (area1 + area2);
							//double w2 = (1 - w1);
							intensity.rC = xBucket[p][0].piecePoly.vertexList[0].r*w1;// +xBucket[p][1].piecePoly.vertexList[0].r*w2;
							intensity.gC = xBucket[p][0].piecePoly.vertexList[0].g*w1;// +xBucket[p][1].piecePoly.vertexList[0].g*w2;
							intensity.bC = xBucket[p][0].piecePoly.vertexList[0].b*w1;// +xBucket[p][1].piecePoly.vertexList[0].b*w2;
						}
						osuWritePixel(p, s, intensity.rC * 255, intensity.gC * 255, intensity.bC * 255);
					}
					else {
						/*double w1 = findArea(xBucket[p][0].piecePoly);
						if (w1 > 1 || w1<0) {
							cout << "abnormal case";
							w1 = 1;
						}
						intensity.rC = xBucket[p][0].piecePoly.vertexList[0].r*w1 ;
						intensity.gC = xBucket[p][0].piecePoly.vertexList[0].g*w1 ;
						intensity.bC = xBucket[p][0].piecePoly.vertexList[0].b*w1 ;
						osuWritePixel(p, s, intensity.rC * 255, intensity.gC * 255, intensity.bC * 255);*/
					}
				}

			}
			else if (polyCount == 1) {
				if (xBucket[p][0].isSolid == 0) {
					double w1 = findArea(xBucket[p][0].piecePoly);
					/*if(w1==1)
					cout << w1 << " ";*/
					if (w1 > 1) {
						//cout << "abnormal case\n";
						w1 = 1;
					}
					intensity.rC = xBucket[p][0].piecePoly.vertexList[0].r*w1;
					intensity.gC = xBucket[p][0].piecePoly.vertexList[0].g*w1;
					intensity.bC = xBucket[p][0].piecePoly.vertexList[0].b*w1;
				}
				else {
					double area1 = findArea(xBucket[p][0].piecePoly);
					/*if (area1 < 1)
						cout << "solid but less\n";*/
					intensity.rC = xBucket[p][0].piecePoly.vertexList[0].r;
					intensity.gC = xBucket[p][0].piecePoly.vertexList[0].g;
					intensity.bC = xBucket[p][0].piecePoly.vertexList[0].b;
				}
				osuWritePixel(p, s, intensity.rC*255, intensity.gC * 255, intensity.bC * 255);
			}
			else {
				osuWritePixel(p, s, 0, 0, 0);
			}
			
		}
	}
	
}

void osuNormal3f(double x, double y, double z) {
	globalNormal[0] = x;
	globalNormal[1] = y;
	globalNormal[2] = z;
}
//shading

void osuShadeModel(int model) {
	if (model == OSU_FLAT)
		shadeModel = 1;
	else if (model == OSU_SMOOTH)
		shadeModel = 2;
}

void osuPointLight(float pos[3], float i)
{
	pointLight tempPLight;
	tempPLight.pI = i;
	tempPLight.pl0 = pos[0];
	tempPLight.pl1 = pos[1];
	tempPLight.pl2 = pos[2];
	pointLightList.push_back(tempPLight);
}

void osuDirectionalLight(float dir[3], float i)
{
	directLight tempDLight;
	tempDLight.dI = i;
	tempDLight.dl0 = -dir[0];
	tempDLight.dl1 = -dir[1];
	tempDLight.dl2 = -dir[2];
	directLightList.push_back(tempDLight);
}

void osuAmbientLight(float i)
{
	ambientIntensity = i;
	ambient = 1;
}

void osuDiffuse(float r, float g, float b)
{
	rDiffuse = r;
	gDiffuse = g;
	bDiffuse = b;
	diffuse = 1;
}

void osuSpecular(float r, float g, float b, float s)
{
	rSpecular = r;
	gSpecular = g;
	bSpecular = b;
	phongCoefficient = s;
}

void osuColor3f(double red, double green, double blue)
{
	currentColorState.rC = red;
	currentColorState.gC = green;
	currentColorState.bC = blue;
}

void osuVertex2f(double x, double y) {
	osuVertex2 vertex;
	vertex.x = x; vertex.y = y;
	vertex.rV = currentColorState.rC;
	vertex.gV = currentColorState.gC;
	vertex.bV = currentColorState.bC;
	if (drawPolygon == 1) {
		draw2D = 1;
		draw3D = 0;
	}

	vertexList.push_back(vertex);
}

void osuEnable(int depthTestBit)
{

	int w, h;
	depthTest = 1;
	osuGetFramebufferSize(&w, &h);

	zBuffer = (double**)malloc(sizeof(double*)*h);
	for (int i = 0; i<h; i++)
		zBuffer[i] = (double*)malloc(sizeof(double)*w);

	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			zBuffer[i][j] = -10000;
		}
	}

}

void osuClearZ()
{
	int w, h;
	osuGetFramebufferSize(&w, &h);
	if (depthTest == 1) {
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				zBuffer[i][j] = -10000;
			}
		}
	}
}

void osuVertex3f(double x, double y, double z) {
	osuVertex vertex;
	vertex.x = x; vertex.y = y, vertex.z = z;
	vertex.r = currentColorState.rC;
	vertex.g = currentColorState.gC;
	vertex.b = currentColorState.bC;
	vertex.nX = globalNormal[0];
	vertex.nY = globalNormal[1];
	vertex.nZ = globalNormal[2];

	if (drawLine == 1) {
		vertexList3D.push_back(vertex);
		if (vertexList3D.size() == 2) {
			//do trans, projection, call drawlines
			vector4 v1 = { vertexList3D[0].x, vertexList3D[0].y, vertexList3D[0].z, 1.0 };
			vector4 v2 = { vertexList3D[1].x, vertexList3D[1].y, vertexList3D[1].z, 1.0 };

			matvectMult(currentTransMatrix, v1);
			matvectMult(currentViewMAtrix, v1);
			matvectMult(currentProjectMatrix, v1);

			matvectMult(currentTransMatrix, v2);
			matvectMult(currentViewMAtrix, v2);
			matvectMult(currentProjectMatrix, v2);

			double x1, y1, x2, y2, z1, z2;
			int widthWin, heightWin;
			x1 = v1[0] / v1[3]; y1 = v1[1] / v1[3]; z1 = v1[2] / v1[3];
			x2 = v2[0] / v2[3]; y2 = v2[1] / v2[3]; z2 = v2[2] / v2[3];
			getWindowSize(&widthWin, &heightWin);
			int clipped = near_far_clip(-1, 1, &x1, &y1, &z1, &x2, &y2, &z2);
			if (clipped) {
				x1 = widthWin / 2 + x1* (widthWin / 2);
				x2 = widthWin / 2 + x2* (widthWin / 2);
				y1 = heightWin / 2 + y1* (heightWin / 2);
				y2 = heightWin / 2 + y2* (heightWin / 2);

				draw_line(x1, y1, x2, y2);
			}

			vertexList3D.clear();
		}
	}
	else if (drawPolygon == 1) {
		draw2D = 0;
		draw3D = 1;
		float rL, gL, bL;

		icVector3 Normal(vertex.nX, vertex.nY, vertex.nZ);
		icVector3 pos(x, y, z);
		icVector3 eyeVec = eyePos - pos;
		normalize(Normal);
		normalize(eyeVec);

		rL = ambientIntensity*rDiffuse;
		gL = ambientIntensity*gDiffuse;
		bL = ambientIntensity*bDiffuse;

		for (int i = 0; i < pointLightList.size(); i++) {
			icVector3 lightPos(pointLightList[i].pl0, pointLightList[i].pl1, pointLightList[i].pl2);
			
			icVector3 lightVec = lightPos - pos;

			normalize(lightVec);

			rL += pointLightList[i].pI*rDiffuse*max(0.0, dot(Normal, lightVec));
			gL += pointLightList[i].pI*gDiffuse*max(0.0, dot(Normal, lightVec));
			bL += pointLightList[i].pI*bDiffuse*max(0.0, dot(Normal, lightVec));

			icVector3 r = 2 * Normal*(dot(Normal, lightVec)) - lightVec;
			float cp = pow(max(dot(eyeVec, r), 0.0), phongCoefficient);
			rL += pointLightList[i].pI*rSpecular*cp;
			gL += pointLightList[i].pI*gSpecular*cp;
			bL += pointLightList[i].pI*bSpecular*cp;
		}

		for (int i = 0; i < directLightList.size(); i++) {
			icVector3 lightVec(directLightList[i].dl0, directLightList[i].dl1, directLightList[i].dl2);
			//icVector3 pos(x, y, z);
			//icVector3 lightVec = lightPos - pos;

			normalize(lightVec);

			rL += directLightList[i].dI*rDiffuse*max(0.0, dot(Normal, lightVec));
			gL += directLightList[i].dI*gDiffuse*max(0.0, dot(Normal, lightVec));
			bL += directLightList[i].dI*bDiffuse*max(0.0, dot(Normal, lightVec));

			icVector3 r = 2 * Normal*(dot(Normal, lightVec)) - lightVec;
			float cp = pow(max(dot(eyeVec, r), 0.0), phongCoefficient);
			rL += directLightList[i].dI*rSpecular*cp;
			gL += directLightList[i].dI*gSpecular*cp;
			bL += directLightList[i].dI*bSpecular*cp;
		}
		if (pointLightList.size() != 0 || directLightList.size() != 0 || (ambient && diffuse)) {
			vertex.r = rL;
			vertex.g = gL;
			vertex.b = bL;
		}
		vertexList3D.push_back(vertex);
	}
}

void osuInitialize()
{
	stackCounter = 0;
	depthTest = 0;
	matrix identity;
	loadIdentity(identity);
	//copyMatrix(*identity, &currentTransMatrix);
	//osuPushMatrix();
	//copyMatrix(identity, matrixStack[stackCounter]);
	//currentColorState.rC = 0; currentColorState.gC = 0; currentColorState.bC = 0;
	int i, j;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			matrixStack[stackCounter][i][j] = identity[i][j];
		}
	}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			currentTransMatrix[i][j] = identity[i][j];
		}
	}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			currentProjectMatrix[i][j] = identity[i][j];
		}
	}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			currentViewMAtrix[i][j] = identity[i][j];
		}
	}
}

void osuPushMatrix()
{
	if (stackCounter >= stackSize)
		cout << "ERROR: STACK OVERFLOW";
	else {

		//copyMatrix(currentTransMatrix, matrixStack[stackCounter]);
		stackCounter++;
		int i, j;
		for (i = 0; i < 4; i++) {
			for (j = 0; j < 4; j++) {
				matrixStack[stackCounter][i][j] = currentTransMatrix[i][j];
			}
		}

	}

}

void osuPopMatrix()
{
	if (stackCounter <= 0) { //at least identity needs to be in stack
		cout << "ERROR: STACK UNDERFLOW";
	}
	else {

		int i, j;
		for (i = 0; i < 4; i++) {
			for (j = 0; j < 4; j++) {
				currentTransMatrix[i][j] = matrixStack[stackCounter][i][j];
			}
		}
		stackCounter--;
	}
}

void osuLoadIdentityMatrix()
{
	matrix identity;
	loadIdentity(identity);
	/*copyMatrix(identity, currentTransMatrix);
	copyMatrix(identity, currentProjectMatrix);*/

	int i, j;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			currentTransMatrix[i][j] = identity[i][j];
		}
	}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			currentProjectMatrix[i][j] = identity[i][j];
		}
	}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			currentViewMAtrix[i][j] = identity[i][j];
		}
	}
}

void osuTranslate(double tx, double ty, double tz)
{
	matrix id;
	loadIdentity(id);
	matrix translate, multM;
	//copyMatrix(id, translate);
	int i, j;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			translate[i][j] = id[i][j];
		}
	}
	translate[3][0] = tx;
	translate[3][1] = ty;
	translate[3][2] = tz;
	rightmultMatrix(multM, currentTransMatrix, translate);
	//copyMatrix(*multM, currentTransMatrix);

	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			currentTransMatrix[i][j] = multM[i][j];
		}
	}

	
}

void osuScale(double sx, double sy, double sz)
{
	matrix id;
	loadIdentity(id);
	matrix scale, multM;
	//copyMatrix(id, scale);
	int i, j;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			scale[i][j] = id[i][j];
		}
	}
	scale[0][0] = sx;
	scale[1][1] = sy;
	scale[2][2] = sz;
	rightmultMatrix(multM, currentTransMatrix, scale);
	//copyMatrix(*multM, currentTransMatrix);

	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			currentTransMatrix[i][j] = multM[i][j];
		}
	}

	

}

void osuRotate(double angle, double ax, double ay, double az)
{
	int i, j;
	
	icVector3 rotateAxis(ax, ay, az);
	normalize(rotateAxis);
	double xsquare, ysquare, zsquare, xy, yz, xz;
	xsquare = rotateAxis.x*rotateAxis.x;
	ysquare = rotateAxis.y*rotateAxis.y;
	zsquare = rotateAxis.z*rotateAxis.z;
	xy = rotateAxis.x*rotateAxis.y; 
	yz= rotateAxis.y*rotateAxis.z; 
	xz = rotateAxis.x*rotateAxis.z;
	

	double ang = angle*3.1416 / 180;
	double sn = sin(ang);
	double cs = cos(ang);
	

	matrix rotate = {
		xsquare * (1 - cs) + cs          ,  yz * (1 - cs) + rotateAxis.z * sn,  xz * (1 - cs) - rotateAxis.y * sn,  0,
		xy * (1 - cs) - rotateAxis.z * sn,  ysquare * (1 - cs) + cs          ,  yz * (1 - cs) + rotateAxis.x * sn,  0,
		xz * (1 - cs) + rotateAxis.y * sn,  yz * (1 - cs) - rotateAxis.x * sn,  zsquare * (1 - cs) + cs          ,  0,
		0                         ,  0                         , 0                          ,  1
	};
	matrix multMat;
	rightmultMatrix(multMat, currentTransMatrix, rotate);
	//copyMatrix(*multMat, currentTransMatrix);
	
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			currentTransMatrix[i][j] = multMat[i][j];
		}
	}

	
}

void matvectMult(matrix A, vector4 &V) {
	int i, j;
	vector4 t = { 0 };
	/*for (i = 0; i < 4; i++) {
	for (j = 0; j < 4; j++) {
	mult[i] += A[i][j] * V[j];
	}
	}*/
	for (j = 0; j <= 3; j++)
		t[j] = V[0] * A[0][j] + V[1] * A[1][j] + V[2] * A[2][j] + V[3] * A[3][j];

	for (j = 0; j < 4; j++) {
		V[j] = t[j];
	}
}

osuColor3 getColorLinearInterpolation(double x, double y, osuVertex2 V1, osuVertex2 V2, osuVertex2 V3) {
	osuColor3 xyColor;
	double bary1, bary2;
	bary1 = ((V2.y - V3.y)*(x - V3.x) + (V3.x - V2.x)*(y - V3.y))
		/ ((V2.y - V3.y)*(V1.x - V3.x) + (V3.x - V2.x)*(V1.y - V3.y));
	bary2 = ((V3.y - V1.y)*(x - V3.x) + (V1.x - V3.x)*(y - V3.y))
		/ ((V2.y - V3.y)*(V1.x - V3.x) + (V3.x - V2.x)*(V1.y - V3.y));
	xyColor.rC = bary1*V1.rV + bary2*V2.rV + (1 - bary1 - bary2)*V3.rV;
	xyColor.gC = bary1*V1.gV + bary2*V2.gV + (1 - bary1 - bary2)*V3.gV;
	xyColor.bC = bary1*V1.bV + bary2*V2.bV + (1 - bary1 - bary2)*V3.bV;

	xyColor.rC = xyColor.rC * 255;
	xyColor.gC = xyColor.gC * 255;
	xyColor.bC = xyColor.bC * 255;

	return xyColor;
}

double functionimplicit(double x, double y, osuVertex2 Va, osuVertex2 Vb) {
	double fab_xy = (Va.y - Vb.y)*x + (Vb.x - Va.x)*y + (Va.x*Vb.y - Vb.x*Va.y);
	return fab_xy;
}

int checkinside2(double x, double y, osuVertex2 Va, osuVertex2 Vb, osuVertex2 Vc, int excludeCorners) {

	////f(x, y) = (y0 - y1)*x + (x1 - x0)*y + (x0*y1 - x1*y0) abc 123
	////fac - V1V3 
	//double f13 = (V1.y - V3.y)*x + (V3.x - V1.x)*y + (V1.x*V3.y - V3.x*V1.y);
	////fbc - 23
	//double f23 = (V2.y - V3.y)*x + (V3.x - V2.x)*y + (V2.x*V3.y - V3.x*V2.y);
	////fab - 12
	//double f12 = (V1.y - V2.y)*x + (V2.x - V1.x)*y + (V1.x*V2.y - V2.x*V1.y);
	////fac(b)
	//double f13 = (V1.y - V3.y)*x + (V3.x - V1.x)*y + (V1.x*V3.y - V3.x*V1.y);

	double falpha = functionimplicit(Va.x, Va.y, Vb, Vc);
	double fbeta = functionimplicit(Vb.x, Vb.y, Vc, Va);
	double fgamma = functionimplicit(Vc.x, Vc.y, Va, Vb);

	double alpha = functionimplicit(x, y, Vb, Vc) / falpha;
	double beta = functionimplicit(x, y, Vc, Va) / fbeta;
	double gamma = functionimplicit(x, y, Va, Vb) / fgamma;

	int rx = -1, ry = -1;
	double frefbc = functionimplicit(rx, ry, Vb, Vc);
	double frefca = functionimplicit(rx, ry, Vc, Va);
	double frefab = functionimplicit(rx, ry, Va, Vb);

	while (frefbc == 0 || frefca == 0 || frefca == 0) {
		ry = ry + .01;
		frefbc = functionimplicit(rx, ry, Vb, Vc);
		frefca = functionimplicit(rx, ry, Vc, Va);
		frefab = functionimplicit(rx, ry, Va, Vb);
	}

	if (alpha >= 0 && beta >= 0 && gamma >= 0) {
		/*if( (alpha > 0 || falpha*functionimplicit(-1, -1, Vb, Vc) > 0) &&
		(beta  > 0 || fbeta *functionimplicit(-1, -1, Vc, Va) > 0) &&
		(gamma > 0 || fgamma*functionimplicit(-1, -1, Va, Vb) > 0) ) {

		return 1;
		}*/
		if ((alpha == 1 && beta == 0 && gamma == 0) ||
			(alpha == 0 && beta == 1 && gamma == 0) || 
			(alpha == 0 && beta == 0 && gamma == 1)) {
			if(!excludeCorners)
				return 1;
		}
		if ((alpha > 0 || falpha*frefbc > 0) &&
			(beta  > 0 || fbeta *frefca > 0) &&
			(gamma > 0 || fgamma*frefab > 0)) {

			return 1;
		}
	}
	else
		return 0;
}

//osuColor3 getColorLinearInterpolation(double x, double y, osuVertex V1, osuVertex V2, osuVertex V3)
//{
//	osuColor3 b;
//	return b;
//}
//
//int checkinside2(double x, double y, osuVertex Va, osuVertex Vb, osuVertex Vc)
//{
//	
//		return 0;
//}

void getBaryCoord(double x, double y, osuVertex V1, osuVertex V2, osuVertex V3, double &bary_a, double &bary_b) {

	bary_a = ((V2.y - V3.y)*(x - V3.x) + (V3.x - V2.x)*(y - V3.y))
		/ ((V2.y - V3.y)*(V1.x - V3.x) + (V3.x - V2.x)*(V1.y - V3.y));
	bary_b = ((V3.y - V1.y)*(x - V3.x) + (V1.x - V3.x)*(y - V3.y))
		/ ((V2.y - V3.y)*(V1.x - V3.x) + (V3.x - V2.x)*(V1.y - V3.y));
}

void osuPerspective(double fovy, double nearp, double farp)
{
	projectionOrtho = 1;
	orhtoPers = 2;
	double f = 1 / (tan((fovy / 2)*3.1416 / 180));
	currentProjectMatrix[0][0] = f;
	currentProjectMatrix[1][1] = f;
	currentProjectMatrix[2][2] = (nearp + farp) / (nearp - farp);
	currentProjectMatrix[2][3] = -1;
	currentProjectMatrix[3][2] = (2 * nearp*farp) / (nearp - farp);
	currentProjectMatrix[3][3] = 0;

}
void osuLookat(float from[3], float at[3], float up[3])
{
	icVector3 At(at[0], at[1], at[2]);
	icVector3 Eye(from[0], from[1], from[2]);
	icVector3 Up(up[0], up[1], up[2]);
	eyePos.x = Eye.x; eyePos.y = Eye.y; eyePos.z = Eye.z;
	icVector3 zaxis = At - Eye;
	normalize(zaxis);

	icVector3 xaxis = cross(Up, zaxis);
	normalize(xaxis);

	icVector3 yaxis = -cross(zaxis, xaxis);
	normalize(yaxis);
	
	matrix lMatrix = {
		xaxis.x ,  yaxis.x ,  zaxis.x ,  0,
		xaxis.y ,  yaxis.y ,  zaxis.y ,  0,
		xaxis.z  , yaxis.z ,  zaxis.z ,  0,
		-dot(xaxis, Eye)      ,  -dot(yaxis, Eye)       ,  -dot(zaxis, Eye)       ,   1
	};
	matrix currentmultV;
	rightmultMatrix(currentmultV, currentTransMatrix, lMatrix);
	//copyMatrix(*currentmultV, currentTransMatrix);
	int i, j;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			currentTransMatrix[i][j] = currentmultV[i][j];
		}
	}
}