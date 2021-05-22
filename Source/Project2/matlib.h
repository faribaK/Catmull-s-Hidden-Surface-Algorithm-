/*

Line and matrix header.

*/
#include "osuGraphics.h"
void draw_line(double, double, double, double);
int near_far_clip(double,double,double*,double*,double*,double*,double*,double*);
void osuInitialize();
void osuPushMatrix();
void osuPopMatrix();
void osuLoadIdentityMatrix();
void osuTranslate(double, double, double);
void osuScale(double, double, double);
void osuRotate(double angle, double ax, double ay, double az);
void osuOrtho(double, double, double, double, double, double);
void osuPerspective(double, double, double);
void osuLookat(float from[3], float at[3], float up[3]);
void osuBegin(OSUDrawable);
void osuEnd();
void osuVertex3f(double, double, double);
void osuNormal3f(double x, double y, double z);
void osuShadeModel(int model);
void osuPointLight(float pos[3], float i);
void osuDirectionalLight(float dir[3], float i);
void osuAmbientLight(float i);
void osuDiffuse(float r, float g, float b);
void osuSpecular(float r, float g, float b, float s);
void osuColor3f(double red, double green, double blue);
void osuVertex2f(double x, double y);
void osuEnable(int depthTestBit);
void osuClearZ();
void drawScene();



