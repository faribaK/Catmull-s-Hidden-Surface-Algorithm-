
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "osuGraphics.h"
#include "matlib.h"
#include "lines.h"
#include <Windows.h>
#include <iostream>
#include <vector>
#include "ObjLoader.h"
#include <gl\GL.h>
#include "glut.h"
using namespace std;


int xsize = 360;
int ysize = 360;
static int point1 = 0, ambient1 = 0, left1 = 0, right1 = 0, ftop = 0, fbottom = 0, yellowSpecular = 0, smooth1 = 1, flat1 = 0, dirL = 0;

void simpleTest()
{
	osuEnable(OSU_DEPTH_TEST);
	osuPerspective(90.0, 1.0, -1000.);

	float from[3] = { 3.0,3.0,3.0 };
	float at[3] = { 0.0,0.0,-8.0 };
	float up[3] = { 0.0, 1.0, 0.0 };

	osuLookat(from, at, up);
	osuClear(0.0, 0.0, 0.0);

	//osuSpecular(1.0, 1.0, 1.0, 1.0);

	float lpos[3]={0.0, 1.5, 5.0};
	osuPointLight(lpos,0.7);

	osuDiffuse(0.0, 0.0, 1.0);
	float dir[3] = { -2.0, -0.0, -3.0 };
	osuDirectionalLight(dir, 0.7);
	osuAmbientLight(0.2);
	osuNormal3f(0.0, 0.0, 1.0);

	//YOU MUST CONVERT THIS TO TWO TRIANGLES!!!

	osuBegin(OSU_POLYGON);
	osuVertex3f(-4.5, -1.75, -5.5);
	osuVertex3f(-4.5, 1.75, -5.5);
	osuVertex3f(4.5, 1.75, -5.5);
	osuVertex3f(4.5, -1.75, -5.5);
	osuEnd();

	osuNormal3f(1.0, 0.0, 0.0);
	osuDiffuse(1.0, 0.0, 0.0);
	
	osuColor3f(0.0, 0.0, 1.0);
	osuBegin(OSU_POLYGON);
	osuVertex3f(0.0, -1.75, -2.5);
	osuVertex3f(0.0, 1.75, -2.5);
	osuVertex3f(0.0, 1.75, -7.5);
	osuVertex3f(0.0, -1.75, -7.5);
	osuEnd();
}

void color_polygons()
{
	double x0 = 0.6;
	double y0 = 0.6;
	double x1 = 0.9;
	double y1 = 0.9;
	double x2, y2;
	double z1 = 0.2, z2 = .5, z3 = .7, z4 = .9;

	/* colorful triangle */

	osuBegin(OSU_TRIANGLE);
	osuColor3f(1.0, 0.0, 0.0);
	osuVertex3f(0.1, 0.1,z1);
	//osuColor3f(0.0, 1.0, 0.0);
	osuVertex3f(0.1, 0.5, z1);
	//osuColor3f(0.0, 0.0, 1.0);
	osuVertex3f(0.5, 0.3, z1);
	osuEnd();

	osuBegin(OSU_TRIANGLE);
	osuColor3f(1.0, 0.0, 0.0);
	osuVertex3f(x0-.5, y0, z2);
	//osuColor3f(0.0, 1.0, 0.0);
	osuVertex3f(x0-.5, y1, z2);
	//osuColor3f(0.0, 0.0, 1.0);
	osuVertex3f(x1-.5, y1, z2);
	osuEnd();


	/* colors for square */
	osuBegin(OSU_TRIANGLE);
	osuColor3f(0.0, 1.0, 0.0);
	osuVertex3f(x0, y0, z2);
	//osuColor3f(0.0, 1.0, 0.0);
	osuVertex3f(x0, y1, z2);
	//osuColor3f(0.0, 0.0, 1.0);
	osuVertex3f(x1, y1, z2);
	osuEnd();

	osuBegin(OSU_TRIANGLE);
	//osuColor3f(1.0, 0.0, 0.0);
	osuVertex3f(x0, y0, z2);
	//osuColor3f(0.0, 0.0, 1.0);
	osuVertex3f(x1, y1, z2);
	//osuColor3f(1.0, 1.0, 1.0);
	osuVertex3f(x1, y0, z2);
	osuEnd();


	x0 = 0.55;
	y0 = 0.15;
	x1 = 0.7;
	y1 = 0.3;
	x2 = 0.85;
	y2 = 0.45;

	osuBegin(OSU_TRIANGLE);
	osuColor3f(1.0, 0.0, 1.0);
	osuVertex3f(x0, y1, z3);
	//osuColor3f(1.0, 1.0, 1.0);
	osuVertex3f(x1, y0, z3);
	//osuColor3f(1.0, 0.0, 0.0);
	osuVertex3f(x2, y1, z3);
	osuEnd();

	osuBegin(OSU_TRIANGLE);
	//osuColor3f(1.0, 1.0, 0.0);
	osuVertex3f(x0, y1, z3);
	//osuColor3f(1.0, 0.0, 0.0);
	osuVertex3f(x2, y1, z3);
	//osuColor3f(1.0, 1.0, 1.0);
	osuVertex3f(x1, y2, z3);
	osuEnd();

	x0 = 0.15;
	y0 = 0.55;
	x1 = 0.3;
	y1 = 0.7;
	x2 = 0.45;
	y2 = 0.85;

	//osuBegin(OSU_TRIANGLE);
	//osuColor3f(0.0, 1.0, 1.0);
	//osuVertex3f(x0, y1, z4);
	////osuColor3f(1.0, 0.0, 0.0);
	//osuVertex3f(x1, y0, z4);
	////osuColor3f(1.0, 1.0, 1.0);
	//osuVertex3f(x2, y1, z4);
	//osuEnd();

	//osuBegin(OSU_TRIANGLE);
	////osuColor3f(1.0, 1.0, 1.0);
	//osuVertex3f(x0, y1, z4);
	////osuColor3f(1.0, 1.0, 1.0);
	//osuVertex3f(x2, y1, z4);
	////osuColor3f(1.0, 0.0, 0.0);
	//osuVertex3f(x1, y2, z4);
	//osuEnd();

	//glBegin(GL_TRIANGLES);
	//glColor3f(0.0, 1.0, 1.0);
	//glVertex3f(x0, y1, z4);
	////glColor3f(1.0, 0.0, 0.0);
	//glVertex3f(x1, y0, z4);
	////glColor3f(1.0, 1.0, 1.0);
	//glVertex3f(x2, y1, z4);
	//glEnd();

	//glBegin(GL_TRIANGLES);
	////glColor3f(1.0, 1.0, 1.0);
	//glVertex3f(x0, y1, z4);
	////glColor3f(1.0, 1.0, 1.0);
	//glVertex3f(x2, y1, z4);
	////glColor3f(1.0, 0.0, 0.0);
	//glVertex3f(x1, y2, z4);
	//glEnd();
}

void Cube()
{

	osuEnable(OSU_DEPTH_TEST);
	osuPerspective(40, 7.5, 100);

	float from[3] = { 5.0,5.0,3.0 };
	float at[3] = { 0.0,0.0,0.0 };
	float up[3] = { 0.0, 1.0, 0.0 };

	osuLookat(from, at, up);
	osuClear(0.0, 0.0, 0.0);
	osuClearZ();

	//osuDiffuse(0.0, 0.0, 1.0);
	//osuSpecular(1.0, 1.0, 1.0, 1.0);

	float lpos[3] = { 3.0, 3.5, 3.0 };
	float dir[3] = { -3.0, -1.0, -2.0 };
	float dir1[3] = { 2.0, -0.0, -3.0 };

//	osuDirectionalLight(dir1, 0.7);

	osuShadeModel(OSU_SMOOTH);
	//osuPointLight(lpos, 0.5);
	//osuDirectionalLight(dir, 0.7);
	//osuAmbientLight(0.2);
	osuNormal3f(0.0, 1.0, 0.0);

	//YOU MUST CONVERT THESE TO USE TRIANGLES!!!

	//back
	osuColor3f(1, 0, 0);
	osuNormal3f(1.0, 1.0, 0.0);
	osuBegin(OSU_POLYGON);
	osuVertex3f(-1, -1, -1);
	osuVertex3f(1, -1, -1);
	osuVertex3f(1, 1, -1);
	osuVertex3f(-1, 1, -1);
	osuEnd();
	//osuDiffuse(0.0, 0.0, 1.0);
	osuNormal3f(1.0, 0.0, 0.0);

	//right
	osuBegin(OSU_POLYGON);
	osuVertex3f(1, -1, -1);
	osuVertex3f(1, -1, 1);
	osuVertex3f(1, 1, 1);
	osuVertex3f(1, 1, -1);
	osuEnd();
	//osuDiffuse(0.5, 0.5, 0.0);
	osuNormal3f(0.0, 0.0, 1.0);

	//front
	osuBegin(OSU_POLYGON);
	osuVertex3f(-1, -1, 1);
	osuVertex3f(-1, 1, 1);
	osuVertex3f(1, 1, 1);
	osuVertex3f(1, -1, 1);
	osuEnd();
	//osuDiffuse(0.5, 0.5, 0.0);
	osuNormal3f(0.0, 1.0, 0.0);

	//top
	osuBegin(OSU_POLYGON);
	osuVertex3f(-1, 1, -1);
	osuVertex3f(1, 1, -1);
	osuVertex3f(1, 1, 1);
	osuVertex3f(-1, 1, 1);
	osuEnd();
	//osuDiffuse(0.5, 0.5, 0.0);
	osuNormal3f(0.0, 1.0, 0.0);

	//bottom
	osuBegin(OSU_POLYGON);
	osuVertex3f(-1, -1, -1);
	osuVertex3f(-1, -1, 1);
	osuVertex3f(1, -1, 1);
	osuVertex3f(1, -1, -1);
	osuEnd();
	//osuDiffuse(0.5, 0.5, 0.0);
	osuNormal3f(1.0, 0.0, 0.0);

	//left
	osuBegin(OSU_POLYGON);
	osuVertex3f(-1, -1, -1);
	osuVertex3f(-1, 1, -1);
	osuVertex3f(-1, 1, 1);
	osuVertex3f(-1, -1, 1);
	osuEnd();
}

void loadAndDrawObj(char *fname)
{
	ObjModel data;
	ObjLoader LoaderClass;

	LoaderClass.LoadObj(fname);
	data = LoaderClass.ReturnObj();

	//For flat shading, only make one call to osuNormal in the beginning, and use the following
	// to access the faceNormal
	//
	//
	

	for (int i = 0; i < data.NumTriangle; i++) {

		/*osuNormal3f(data.TriangleArray[i].faceNormal[0],
			data.TriangleArray[i].faceNormal[1],
				data.TriangleArray[i].faceNormal[2]);*/
		osuBegin(OSU_POLYGON);

		osuNormal3f(data.NormalArray[data.TriangleArray[i].Vertex[0]].X,
			data.NormalArray[data.TriangleArray[i].Vertex[0]].Y,
			data.NormalArray[data.TriangleArray[i].Vertex[0]].Z);

		osuVertex3f(data.VertexArray[data.TriangleArray[i].Vertex[0]].X,
			data.VertexArray[data.TriangleArray[i].Vertex[0]].Y,
			data.VertexArray[data.TriangleArray[i].Vertex[0]].Z);


		osuNormal3f(data.NormalArray[data.TriangleArray[i].Vertex[1]].X,
			data.NormalArray[data.TriangleArray[i].Vertex[1]].Y,
			data.NormalArray[data.TriangleArray[i].Vertex[1]].Z);
		osuVertex3f(data.VertexArray[data.TriangleArray[i].Vertex[1]].X,
			data.VertexArray[data.TriangleArray[i].Vertex[1]].Y,
			data.VertexArray[data.TriangleArray[i].Vertex[1]].Z);

		osuNormal3f(data.NormalArray[data.TriangleArray[i].Vertex[2]].X,
			data.NormalArray[data.TriangleArray[i].Vertex[2]].Y,
			data.NormalArray[data.TriangleArray[i].Vertex[2]].Z);
		osuVertex3f(data.VertexArray[data.TriangleArray[i].Vertex[2]].X,
			data.VertexArray[data.TriangleArray[i].Vertex[2]].Y,
			data.VertexArray[data.TriangleArray[i].Vertex[2]].Z);
		osuEnd();

	}

}


void objTest()
{
	osuEnable(OSU_DEPTH_TEST);
	osuPerspective(90.0, 1.0, 1000);
	osuClear(0.0, 0.0, 0.0);
	if(flat1==1)
		osuShadeModel(OSU_FLAT);
	else if(smooth1==1)
		osuShadeModel(OSU_SMOOTH);

	float from[3] = { 0.0,0,2 };
	float at[3] = { 0.0,0.0,0.0 };
	float up[3] = { 0.0, 1.0, 0.0 };


	osuLookat(from, at, up);
	//Diffuse blue color
	osuDiffuse(0.0, 0.0, 1.0);

	//Specular white color
	if(yellowSpecular==1)
		osuSpecular(1.0, 1.0, 0.0, 0.5);
	else
		osuSpecular(1.0, 1.0, 1.0, 0.5);

	float lpos1[3] = { 0.0, 3.0, 5.0 };
	if (point1 == 1)
	osuPointLight(lpos1, 0.7);

	float lpos[3] = { 204.0, 204.0, 204.0 };
	//osuPointLight(lpos, 0.5);
	if(ambient1==1)
		osuAmbientLight(0.2);
	else if(ambient1 == 3)
		osuAmbientLight(0.3);
	else 
		osuAmbientLight(0.0);
	if (dirL == 1) {
		if (left1 == 1) {
			float dir1[3] = { 2.0, -0.0, -3.0 };
			osuDirectionalLight(dir1, 0.7);
		}
			
		else if (right1 == 1) {
			float dir1[3] = { -2.0, -0.0, -3.0 };
			osuDirectionalLight(dir1, 0.7);
		}
			
		else  if (ftop == 1) {
			float dir1[3] = { 0.0, -2.0, -3.0 };
			osuDirectionalLight(dir1, 0.7);
		}
			
		else  if (fbottom == 1) {
			float dir1[3] = { 0.0, 2.0, -3.0 };
			osuDirectionalLight(dir1, 0.7);
		}
			
		else {
			float dir1[3] = { 0.0, 2.0, -3.0 };
			osuDirectionalLight(dir1, 0.7);
		}
		
	}
	
	

	// To test flath shading, use test.obj 
	// Look at it from behind... ie osuLookAt(0,1,-3, 0,0,10,0,1,0)
	//loadAndDrawObj("test.obj");

	// To test smooth shading, use face.ws.obj
	// Look at it from in front osuLookAt(0, 0, 2, 0, 0, 0, 0, 1, 0)
	loadAndDrawObj("face.ws.obj");
}

void overlapTest5()
{
	osuSetWriteMode(OSU_XOR);
	osuTranslate(-1, -1, 0);

	osuBegin(OSU_POLYGON);
	osuColor3f(1, 1, 1);
	osuVertex3f(2 * 0.5, 2 * 0.7, .5);
	osuVertex3f(2 * 0.3, 2 * 0.6, .5);
	osuVertex3f(2 * 0.7, 2 * 0.01, .5);
	osuEnd();

	osuBegin(OSU_POLYGON);
	osuColor3f(1, 1, 1);
	osuVertex3f(2 * 0.5, 2 * 0.2, .4);
	osuVertex3f(2 * 0.8, 2 * 0.2, .4);
	osuVertex3f(2 * 0.5, 2 * 0.4, .4);
	osuVertex3f(2 * 0.8, 2 * 0.4, .4);
	osuEnd();
	drawScene();
}

void overlapTest2()
{
	osuSetWriteMode(OSU_XOR);
	osuTranslate(-1, -1, 0);
	//case pink behind

	osuBegin(OSU_POLYGON);
	osuColor3f(1, 1, 1);
	osuVertex3f(2 * 0.5, 2 * 0.7, .5);
	osuVertex3f(2 * 0.3, 2 * 0.6, .5);
	osuVertex3f(2 * 0.7, 2 * 0.01, .5);
	osuEnd();

	osuBegin(OSU_POLYGON);
	osuColor3f(1, 0, 1);
	osuVertex3f(2 * 0.5, 2 * 0.2, .7);
	osuVertex3f(2 * 0.8, 2 * 0.2, .7);
	osuVertex3f(2 * 0.5, 2 * 0.4, .7);
	osuVertex3f(2 * 0.8, 2 * 0.4, .7);
	osuEnd();
	drawScene();


	//osuBegin(OSU_POLYGON);
	//osuColor3f(1, 1, 1);
	//osuVertex3f(2 * 0.4, 2 * 0.8, .5);
	//osuVertex3f(2 * 0.1, 2 * 0.8, .5);
	//osuVertex3f(2 * 0.85, 2 * 0.1, .5);
	////osuVertex3f(2 * 0.8, 2 * 0.1001, .5);
	//osuEnd();
	
}

void overlapTest6()
{
	osuSetWriteMode(OSU_XOR);
	osuTranslate(-1, -1, 0);
	//case pink in front

	osuBegin(OSU_POLYGON);
	osuColor3f(1, 1, 1);
	osuVertex3f(2 * 0.5, 2 * 0.7, .5);
	osuVertex3f(2 * 0.3, 2 * 0.6, .5);
	osuVertex3f(2 * 0.7, 2 * 0.01, .5);
	osuEnd();

	osuBegin(OSU_POLYGON);
	osuColor3f(1, 0, 1);
	osuVertex3f(2 * 0.5, 2 * 0.2, .4);
	osuVertex3f(2 * 0.8, 2 * 0.2, .4);
	osuVertex3f(2 * 0.5, 2 * 0.4, .4);
	osuVertex3f(2 * 0.8, 2 * 0.4, .4);
	osuEnd();
	drawScene();


	//osuBegin(OSU_POLYGON);
	//osuColor3f(1, 1, 1);
	//osuVertex3f(2 * 0.4, 2 * 0.8, .5);
	//osuVertex3f(2 * 0.1, 2 * 0.8, .5);
	//osuVertex3f(2 * 0.85, 2 * 0.1, .5);
	////osuVertex3f(2 * 0.8, 2 * 0.1001, .5);
	//osuEnd();

}

void overlapTest3()
{
	osuSetWriteMode(OSU_XOR);
	osuTranslate(-1, -1, 0);
	//case 3
	osuBegin(OSU_POLYGON);
	osuColor3f(1, 0, 0);
	osuVertex3f(2 * 0.35, 2 * 0.1, .5);
	osuVertex3f(2 * 0.1, 2 * 0.35, .5);
	osuVertex3f(2 * 0.35, 2 * 0.6, .5);
	osuVertex3f(2 * 0.6, 2 * 0.5, .5);
	osuEnd();

	osuColor3f(0, 1, 0);
	osuBegin(OSU_POLYGON);
	osuVertex3f(2 * 0.65, 2 * 0.1, .5);
	osuVertex3f(2 * 0.4, 2 * 0.35, .5);
	osuVertex3f(2 * 0.65, 2 * 0.6, .5);
	osuVertex3f(2 * 0.9, 2 * 0.5, .5);
	osuEnd(); 
	drawScene();
}

void overlapTest4()
{
	osuSetWriteMode(OSU_XOR);
	osuTranslate(-.5, -.5, 0);
	//case 2
	osuBegin(OSU_TRIANGLE);
	osuColor3f(1, 0, 0);
	osuVertex3f(0.1, 0.1, .5);
	osuVertex3f(0.6, 0.5, .5);
	osuVertex3f(0.1, 0.9, .5);
	//osuVertex3f(0.6, 0.5, .5);
//	osuVertex3f(0.6, 0.1, .5);
	osuEnd();

	osuColor3f(0, 1, 0);
	osuBegin(OSU_TRIANGLE);
	osuVertex3f(0.9, 0.1, .2);
	osuVertex3f(0.4, 0.5, .2);
	osuVertex3f(0.9, 0.9, .2);
	osuEnd();
	drawScene();
}

void overlapTest1()
{
	osuSetWriteMode(OSU_XOR);
	osuTranslate(-1, -1, 0);
	//case1:
	glBegin(OSU_TRIANGLE);
	glColor3f(1, 0, 0);
	glVertex3f(0.1, 0.1, .5);
	glVertex3f(0.6, 0.5, .5);
	glVertex3f(0.1, 0.9, .5);
	glEnd();

	osuColor3f(0, 1, 0);
	osuBegin(OSU_TRIANGLE);
	osuVertex3f(0.9, 0.1, .2);
	osuVertex3f(0.4, 0.5, .2);
	osuVertex3f(0.9, 0.9, .2);
	osuEnd();
	

	drawScene();
}

/******************************************************************************

Test out drawing routines.

******************************************************************************/
void main(int argc, char **argv)
{
	int num;
	cout << "Test Cases presented in Report:" << endl;
	cout << "1: Figure 1(b)" << endl;
	cout << "2: Figure 2(b)" << endl;
	cout << "3: Figure 3(b)" << endl;
	cout << "4: Figure 4(b)" << endl;
	cout << "5: Figure 4(c)" << endl;
	cout << "Enter the choice bewteen 1 to 3:" << endl;
	cin >> num;
	
	
	osuBeginGraphics(xsize, ysize);

	/* inialize the matrix stack*/
	osuInitialize();

	/* select routine to execute */
	switch (num) {
	case 1:
		overlapTest3();
		break;

	case 2:
		
		overlapTest4();
		break;

	case 3:
		overlapTest5();
		break;
	case 4:
		overlapTest2();
		break;
	case 5:
		overlapTest6();
		break;
	

	default:
		fprintf(stderr, "Please use a number from 1 to 5.\n");
		exit(-1);
	}

	osuFlush();

	printf("Press 'escape' to exit.\n");
	osuWaitOnEscape();
	osuEndGraphics();
}

//int main(int argc, char ** argv)
//{
//	
//
//	glutInit(&argc, argv);
//	glutInitWindowSize(400, 400);
//	glutInitWindowPosition(10, 10);
//	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
//	glutCreateWindow("opengL");
//
//	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//
//	glMatrixMode(GL_MODELVIEW);
//	glLoadIdentity();
//	glClearColor(0, 0, 0, 0);
//
//	glBegin(GL_POLYGON);
//	glColor3f(1, 0, 0);
//	glVertex3f(0.35, 0.1, .5);
//	glVertex3f(0.1, 0.5, .5);
//	glVertex3f(0.35, 0.9, .5);
//	glVertex3f(0.6, 0.5, .5);
//	glEnd();
//
//	glColor3f(0, 1, 0);
//	glBegin(GL_TRIANGLES);
//	glVertex3f(0.65, 0.1, .5);
//	glVertex3f(0.4, 0.5, .5);
//	glVertex3f(0.65, 0.9, .5);
//	glVertex3f(0.9, 0.5, .5);
//	glEnd();
//
//	glutSwapBuffers();
//	glFlush();
//	return 0;
//}