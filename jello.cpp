/*

USC/Viterbi/Computer Science
"Jello Cube" Assignment 1 starter code

Your name:
<write your name here>

*/

#include "jello.h"
#include "showCube.h"
#include "input.h"
#include "physics.h"
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <queue>

// camera parameters
double Theta = pi / 6;
double Phi = pi / 6;
double R = 6;

// mouse control
int g_iMenuId;
int g_vMousePos[2];
int g_iLeftMouseButton,g_iMiddleMouseButton,g_iRightMouseButton;

// number of images saved to disk so far
int sprite=0;

namespace Input {
	// these variables control what is displayed on screen
	int shear=0, bend=0, structural=1, pause=0, viewingMode=0, saveScreenToFile=0;
	int texture=0;
}

bool particleExists(int i, int j, int k) {
	return (0 <= i && i <= 7 && 0 <= j && j <= 7 && 0 <= k && k <= 7);
}

struct world jello;
PerformanceCounter counter;
std::deque<double> fpsBuffer;

int windowWidth, windowHeight;

void myinit()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(90.0,1.0,0.01,1000.0);

	// set background color to grey
	glClearColor(0.5, 0.5, 0.5, 0.0);

	glCullFace(GL_BACK);
	glEnable(GL_CULL_FACE);

	glShadeModel(GL_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
		
	GLuint myTexture;
	Pic* img = jello.textureImage;
	// allocate a texture name
	glGenTextures( 1, &myTexture );
	// select our current texture
	glBindTexture( GL_TEXTURE_2D, myTexture );
	
	// when texture area is small, bilinear filter the closest mipmap
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
	                 GL_LINEAR_MIPMAP_NEAREST );
	// when texture area is large, bilinear filter the original
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
	
	// the texture wraps over at the edges (repeat)
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
	
	// build our texture mipmaps
	gluBuild2DMipmaps( GL_TEXTURE_2D, 3, img->nx, img->ny,
	                   GL_RGB, GL_UNSIGNED_BYTE, img->pix );
	pic_free(img);
	jello.textureId = myTexture;
}

void reshape(int w, int h) 
{
// Prevent a divide by zero, when h is zero.
// You can't make a window of zero height.
	if(h == 0)
		h = 1;

	glViewport(0, 0, w, h);

// Reset the coordinate system before modifying
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
// Set the perspective
	double aspectRatio = 1.0 * w / h;
	gluPerspective(60.0f, aspectRatio, 0.01f, 1000.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity(); 

	windowWidth = w;
	windowHeight = h;

	glutPostRedisplay();
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

// camera parameters are Phi, Theta, R
	gluLookAt(R * cos(Phi) * cos (Theta), R * sin(Phi) * cos (Theta), R * sin (Theta),
		0.0,0.0,0.0, 0.0,0.0,1.0);


	/* Lighting */
	/* You are encouraged to change lighting parameters or make improvements/modifications
	to the lighting model . 
	This way, you will personalize your assignment and your assignment will stick out. 
*/

	// global ambient light
	GLfloat aGa[] = { 0.0, 0.0, 0.0, 0.0 };

	// light 's ambient, diffuse, specular
	GLfloat lKa0[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd0[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat lKs0[] = { 1.0, 1.0, 1.0, 1.0 };

	GLfloat lKa1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd1[] = { 1.0, 0.0, 0.0, 1.0 };
	GLfloat lKs1[] = { 1.0, 0.0, 0.0, 1.0 };

	GLfloat lKa2[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd2[] = { 1.0, 1.0, 0.0, 1.0 };
	GLfloat lKs2[] = { 1.0, 1.0, 0.0, 1.0 };

	GLfloat lKa3[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd3[] = { 0.0, 1.0, 1.0, 1.0 };
	GLfloat lKs3[] = { 0.0, 1.0, 1.0, 1.0 };

	GLfloat lKa4[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd4[] = { 0.0, 0.0, 1.0, 1.0 };
	GLfloat lKs4[] = { 0.0, 0.0, 1.0, 1.0 };

	GLfloat lKa5[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd5[] = { 1.0, 0.0, 1.0, 1.0 };
	GLfloat lKs5[] = { 1.0, 0.0, 1.0, 1.0 };

	GLfloat lKa6[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd6[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat lKs6[] = { 1.0, 1.0, 1.0, 1.0 };

	GLfloat lKa7[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd7[] = { 0.0, 1.0, 1.0, 1.0 };
	GLfloat lKs7[] = { 0.0, 1.0, 1.0, 1.0 };

	// light positions and directions
	GLfloat lP0[] = { -1.999, -1.999, -1.999, 1.0 };
	GLfloat lP1[] = { 1.999, -1.999, -1.999, 1.0 };
	GLfloat lP2[] = { 1.999, 1.999, -1.999, 1.0 };
	GLfloat lP3[] = { -1.999, 1.999, -1.999, 1.0 };
	GLfloat lP4[] = { -1.999, -1.999, 1.999, 1.0 };
	GLfloat lP5[] = { 1.999, -1.999, 1.999, 1.0 };
	GLfloat lP6[] = { 1.999, 1.999, 1.999, 1.0 };
	GLfloat lP7[] = { -1.999, 1.999, 1.999, 1.0 };

	// jelly material color

	GLfloat mKa[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat mKd[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat mKs[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mKe[] = { 0.0, 0.0, 0.0, 1.0 };

	/* set up lighting */
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, aGa);
	glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

	// set up cube color
	glMaterialfv(GL_FRONT, GL_AMBIENT, mKa);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mKd);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mKs);
	glMaterialfv(GL_FRONT, GL_EMISSION, mKe);
	glMaterialf(GL_FRONT, GL_SHININESS, 120);

	// macro to set up light i
	#define LIGHTSETUP(i)\
	glLightfv(GL_LIGHT##i, GL_POSITION, lP##i);\
		glLightfv(GL_LIGHT##i, GL_AMBIENT, lKa##i);\
		glLightfv(GL_LIGHT##i, GL_DIFFUSE, lKd##i);\
		glLightfv(GL_LIGHT##i, GL_SPECULAR, lKs##i);\
		glEnable(GL_LIGHT##i)

		LIGHTSETUP (0);
	LIGHTSETUP (1);
	LIGHTSETUP (2);
	LIGHTSETUP (3);
	LIGHTSETUP (4);
	LIGHTSETUP (5);
	LIGHTSETUP (6);
	LIGHTSETUP (7);

	// enable lighting
	glEnable(GL_LIGHTING);    
	glEnable(GL_DEPTH_TEST);

	// show the cube
	showCube(&jello);

	glDisable(GL_LIGHTING);

	// show the bounding box
	showBoundingBox();

	glutSwapBuffers();
}

double getAvgFps() {
	double sum = 0;
	for(int i = 0; i < fpsBuffer.size(); ++i)
		sum += fpsBuffer[i];
	return sum / fpsBuffer.size();
}

void doIdle()
{
	counter.StopCounter();
	fpsBuffer.push_front(1.0 / counter.GetElapsedTime());
	if(fpsBuffer.size() > 64) // limit buffer size
		fpsBuffer.pop_back();
	counter.StartCounter();
	double fps = getAvgFps();
	fps = truncf(fps * 10) / 10; // round to 1 decimal place
	std::ostringstream sfps;
	sfps << fps;
	std::string title = "Jello cube (FPS: " + sfps.str() + ")";
	glutSetWindowTitle(title.c_str());
	
	using Input::saveScreenToFile;
	using Input::pause;

	char s[20]="picxxxx.ppm";
	int i;

	// save screen to file
	s[3] = 48 + (sprite / 1000);
	s[4] = 48 + (sprite % 1000) / 100;
	s[5] = 48 + (sprite % 100 ) / 10;
	s[6] = 48 + sprite % 10;

	if (saveScreenToFile==1)
	{
		saveScreenshot(windowWidth, windowHeight, s);
		saveScreenToFile=0; // save only once, change this if you want continuos image generation (i.e. animation)
		sprite++;
	}

	if (sprite >= 300) // allow only 300 snapshots
	{
		exit(0);	
	}

	if (pause == 0) {
		// insert code which appropriately performs one step of the cube simulation:
		if(jello.integrator[0] == 'R')
			RK4(&jello);
		else if(jello.integrator[0] == 'E')
			Euler(&jello);
	}

	glutPostRedisplay();
}

// math
double norm(const Vec3d& a) {
	return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

double normalize(const Vec3d& src, Vec3d* dest) {
	double length = norm(src);
	dest->x = src.x / length;
	dest->y = src.y / length;
	dest->z = src.z / length;
}

void make(double x, double y, double z, Vec3d* dest) {
	dest->x = x; 
	dest->y = y;
	dest->z = z;
}

void copy(const Vec3d& src, Vec3d* dest) {
	dest->x = src.x;
	dest->y = src.y;
	dest->z = src.z;
}

double dot(const Vec3d& a, const Vec3d& b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

void plus(const Vec3d& src1, const Vec3d& src2, Vec3d* dest) {
	dest->x = src1.x + src2.x;
	dest->y = src1.y + src2.y;
	dest->z = src1.z + src2.z;
}

void minus(const Vec3d& src1, const Vec3d& src2, Vec3d* dest) {
	dest->x = src1.x - src2.x;
	dest->y = src1.y - src2.y;
	dest->z = src1.z - src2.z;
}

void times(const Vec3d& src, double scalar, Vec3d* dest) {
	dest->x = src.x * scalar;
	dest->y = src.y * scalar;
	dest->z = src.z * scalar;
}

void divide(const Vec3d& src, double scalar, Vec3d* dest) {
	dest->x = src.x / scalar;
	dest->y = src.y / scalar;
	dest->z = src.z / scalar;
}

int main (int argc, char ** argv)
{
	if (argc<2)
	{  
		printf ("Oops! You didn't say the jello world file!\n");
		printf ("Usage: %s [worldfile]\n", argv[0]);
		exit(0);
	}

	readWorld(argv[1],&jello);

	glutInit(&argc,argv);

/* double buffered window, use depth testing, 640x480 */
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

	windowWidth = 640;
	windowHeight = 480;
	glutInitWindowSize (windowWidth, windowHeight);
	glutInitWindowPosition (0,0);
	glutCreateWindow ("Jello cube");

/* tells glut to use a particular display function to redraw */
	glutDisplayFunc(display);

/* replace with any animate code */
	glutIdleFunc(doIdle);

/* callback for mouse drags */
	glutMotionFunc(mouseMotionDrag);

/* callback for window size changes */
	glutReshapeFunc(reshape);

/* callback for mouse movement */
	glutPassiveMotionFunc(mouseMotion);

/* callback for mouse button changes */
	glutMouseFunc(mouseButton);

/* register for keyboard events */
	glutKeyboardFunc(keyboardFunc);

/* do initialization */
	myinit();

/* forever sink in the black hole */
	glutMainLoop();

	return(0);
}

