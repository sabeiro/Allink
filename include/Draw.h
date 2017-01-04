#ifndef DISEGNA_H
#define DISEGNA_H

#include "../include/Matematica.h"
#ifdef USE_GL

#include <GL/glut.h>


#if !defined(GLUT_WHEEL_UP)
#  define GLUT_WHEEL_UP   3
#  define GLUT_WHEEL_DOWN 4
#endif
/** Draw provides the basic configuration of the openGL libraries used in every derived program */
class Draw{
 private:
  /// Switch the perspective view to orthogonal
  int pr;//prospettiva
  /// Boh
  int IfShift;
  /// Obsolete
  void Expand(void);
  /// Obsolete
  void spinDisplay(void);
  /// Definition of the first light
  void Illuminazione(void);
  /// Initializes the idle functions
  void init (void);
  /// Obsolete
  void Spot(void);
  /// Fog
  void Nebbia(void);
  /// Boh
  void processEvent(int value);
  /// Another reshape function
  void ChangeSize(GLsizei w,GLsizei h);
  /// Name of the last function called. For debugging
  void Shout(const char *s, ... );
  //Matematica *Mat;
  /// Choosing a GL_SRC value
  int BlendSource;
  /// Choosing a GL_DST value
  int BlendDest;
 public:
  //Draw(int argc,char** argv);
  Draw();
  /// A scene for debugging
  void Draw1(void);
  /// A empty scene (only the list)
  void DMinimal(void);
  /// Boh
  void DrTriangles(int NPoint);
  /// Displays the data stored in pixel
  void ShowImage();
  /// Initial definition of the window
  void Window(int argc,char **argv);
  /// Definition of the scene on which the objects will be drawn
  void DFigure(void);
  /// Not working
  void DTimer(int v);
  /// Principal reshape function
  void Dreshape(int w,int h);
  /// Apply the texture to a square
  int ApplyTexture();
  /// Call the texture in the scene
  int ShowTexture();
  /// Transform the system coordinates
  void Transform();
  /// Line in a cube
  void DrCube();
  /// Choose a different blending function
  void ChooseBlend(int Which);
  /// Put a string
  void PutString(double *Pos,char *String);
  /// Put a string
  void PutString(double Posx,double Posy,double Posz,char *String);
  /// Print the \param n number in the position \param Pos
  void Numera(double *Pos,int n);
  //------------Definition.cpp-------------
  /// Definition of the primitives
  void Lista(int NSquare);
  /// Point
  int DefPoint();
  /// Cube
  int DefCube(int NSquare);
  /// Quad
  int DefQuad(int NSquare);
  /// Cylinder
  int DefCylinder(double Rad,double Height);
  /// Metallic cylinder
  int DefMetalCylinder(double Rad,double Height);
  /// Hexagon
  int DefHexagon();
  /// Griglia
  int DefGriglia();
  /// Wall
  int DefWall();
  /// Arrow
  int DefArrow();
  /// Arrow
  int DefArrowThin();
  /// Define a simple texture
  int DefTexture();
  /// Initializes all the view constants
  void InitConstant();
  /// Calculate the normal
  double Normal(double *v,double *u,double *w,double *n);
  /// Data type for distance/field functions
  typedef void(Draw::*DEPTH_MAP)(double Val,GLfloat *Color);
  /// Pointer to a distance/field function
  DEPTH_MAP Depth_Map;
  /// Pointer to a generic function
  void DepthMap(double Val,GLfloat *Color){(*this.*Depth_Map)(Val,Color);}
  /// Depth map
  void DepthMap1(double Val,GLfloat *Color);
  /// Choose Depth map
  void ChooseDepthMap(int n);
  /// Depth map
  int DefLegend();
  //-------------Control.cpp
  /// To launch the menu
  void Dmouse(int button, int state,int x,int y);
  /// How the scene rotate (Camera view should be implemented)
  void DMouseMove(int x,int y);
  /// Boh
  void Dspecial(int k, int x, int y);
  /// Combines the key with the functions
  void keyboardDraw(unsigned char key);
  /// Quaternion camera implementation
  void CameraQuat();
  /// Movement up-down
  void ChangeSuGiu(GLfloat Movement);
  /// Movement right left
  void ChangeDxSx(GLfloat Movement);
  //----------File.cpp-------------
  /// Check
  void abort_(const char * s, ...);
  /// Write a tiff file of the data in pixel
  int Picture();
  /// Boh
  int WritePixel();
  /// Write a png file of the data in pixel (not working)
  int WritePng();
  /// Write a png file of the data in pixel (uses libpngwriter)
  int WritePngwriter();
  /// Open a image to be store in pixel
  int OpenImage(const char *FileName);
  /// Reads and draw a script file
  void ReadScript();
  /// Reads and applies external configurations
  void ReadConf();
  /// Width and height of the image
  int ImWidth,ImHeight;
  /// Obsolete
  GLfloat spin,angolo, dspin,
    /// Angles
    xa, ya, za,
    /// Orientation of the light
    xf, yf, zf,
    /// Translation, wheel
    xp, yp, zp, zw,
    /// Position of the info string
    xi, yi, zi,
    /// Position of the legend
    xLeg, yLeg, zLeg,
    /// Width of the legend
    dxLeg, dyLeg,
    /// Position of the light0
    xl0, yl0, zl0,
    /// Position of the light1
    xl1, yl1, zl1,
    /// Obsolete
    scale, dscale,tscale,
    /// Background color
    Rback, Gback,Bback,Aback;
  /// Increment visual DxSx, SuGiu
  GLfloat IncrVisDxSx, IncrVisSuGiu;
  /// Angle DxSx, SuGiu
  GLfloat AngleDxSx, AngleSuGiu;
  /// Rescale the three orthogonal directions
  double InvScaleUn;
  /// Finess of the grid
  double GridStep;//Griglia
  /// Refers to the list of a hexagon
  GLuint Hexagon;
  /// Refers to the list of the legend
  GLuint DrLegend;
  /// Refers to the list of the grid
  GLuint Griglia;
  /// Refers to the list of the square
  GLuint Quad;
  /// Refers to the list of the point
  GLuint Point;
  /// Refers to the list of the cylinder
  GLuint Cylinder;
  /// Refers to the list of another cylinder (obsolete)
  GLuint MetalCylinder;
  /** Refers to the list of the total position of the particles
      which will be generated in another program */
  GLuint Particles;
  /// Refers to the list of the objects called by the script file
  GLuint ScriptList;
  /// Refers to the list of a wall
  GLuint GlWall;
  /// Refers to the list of a arrow
  GLuint Arrow;
  /// Refers to the list of the texture
  GLuint Cube;
  /// Refers to the list of the texture
  GLuint Texture;
  /// Center of the frame
  GLuint XCenter;
  /// Center of the frame
  GLuint YCenter;
  /// Puts/removes the box edges
  int la;
  /// Puts/removes the grid
  int gr;
  /// Enables/disables illumination
  int lu;
  /// Enables/disables spot light
  int sp;
  /// Enables/disables fog
  int ne;
  /// Refers to optional different windows
  int MainWindow,SubWindow1,SubWindow2;
  /// Number of frames
  int Diap,tDiap,tDiapBase;
  /// Decides to draw points or spheres
  int IfPoint;//point or Sphere
  /// Removes the info line
  int IfInfo;
  /// Ignores the script file
  int IfScript;
  /// Boh
  int IfImage;
  /// Activate the blending
  int IfBlend;
  /// Activate the illumination for a specific material
  int IfMaterial;
  /// Number of values to divide the edge in squares
  int Values;
  /// Current step for the picture's name
  int Step;
  /// Width of the window
  int WinWidth;
  /// Height of the window
  int WinHeight;
  /// Old x position of the mouse
  int xRem; 
  /// Old y position of the mouse
  int yRem;
  /// Boh
  int ChangeMouse;
  /// Levels of the images data (usually 4=RGBA)
  int NLevel;
  /// Obsolete
  float Diameter,StepDiameter;
  /// Obsolete
  float NanoRad,ExtraDiam;
  /// Box size
  double Edge[3];
  /// Number of lines per edge
  int GridEdge[3];
  /// Cylinder radius
  double ExtRad;
  /// Cylinder height
  double ExtHeight;
  /// Principal image (always allocated)
  GLubyte *pixel;
  /// Characters for the grid
  char *Number;
  /// Boh
  char *frame;
  /// Info line
  char *info;
};

extern void keyboard(unsigned char key,int x, int y);
extern void MenuChoise(int option);
extern void MenuVisual(int option);
extern void ParticleList(void);
extern void ParticleRealTime();
extern void Figure();
extern void Figure1();
extern void Timer(int v);
extern void reshape(int w,int h);
extern void mouse(int button, int state,int x,int y);
extern void MouseMove(int x,int y);
extern void special(int k, int x, int y);
extern void DefCylinder(GLuint Cylinder,double CyEdge,double CyHeight,double InvScaleUn);
static const GLfloat DrAmbientWhite [] = {0.25, 0.25, 0.25, 1.00}; 
static const GLfloat DrAmbientRed   [] = {0.25, 0.00, 0.00, 1.00}; 
static const GLfloat DrAmbientGreen [] = {0.00, 0.25, 0.00, 1.00}; 
static const GLfloat DrAmbientBlue  [] = {0.00, 0.00, 0.25, 1.00}; 
static const GLfloat DrDiffuseWhite [] = {0.75, 0.75, 0.75, 1.00}; 
static const GLfloat DrDiffuseRed   [] = {0.75, 0.00, 0.00, 1.00}; 
static const GLfloat DrDiffuseGreen [] = {0.00, 0.75, 0.00, 1.00}; 
static const GLfloat DrDiffuseBlue  [] = {0.00, 0.00, 0.75, 1.00}; 
static const GLfloat DrSpecularWhite[] = {1.00, 1.00, 1.00, 1.00}; 
static const GLfloat DrSpecularRed  [] = {1.00, 0.25, 0.25, 1.00}; 
static const GLfloat DrSpecularGreen[] = {0.25, 1.00, 0.25, 1.00}; 
static const GLfloat DrSpecularBlue [] = {0.25, 0.25, 1.00, 1.00}; 
#endif//USE_GL
#endif//DISEGNA_H
static float ColorType[12][4] = 
  { {0.,.7,0.,1.},
    {0.,.0,.7,1.},
    {.7,.0,.0,1.},
    {.4,.1,.7,1.},
    {.4,.5,.0,1.},
    {.7,.3,.3,1.},
    {0.,.4,.3,1.},
    {0.,.0,.5,1.},
    {.0,.2,.5,1.},
    {.0,.4,.7,1.},
    {.7,.1,.4,1.},
    {.4,.7,.2,1.}};
