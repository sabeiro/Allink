#include "../../include/Draw.h"
#include <unistd.h>
#include <setjmp.h>

#ifdef __glut_h__
void Draw::ReadConf(){
  Shout("ReadConf");
  FILE *ConfFile;
  char cLine[256];
  double Val[9];
  if( (ConfFile = fopen("Draw.conf","r")) == 0 ){return;}
  for(int k=0;!(fgets(cLine,256,ConfFile)==NULL);k++){
    if(1 == sscanf(cLine,"Pre %lf",Val) )
      pr = (int) Val[0];
    else if(1 ==  sscanf(cLine,"Passo %lf",Val))
      GridStep = (int)Val[0];
    else if(1 == sscanf(cLine,"Box %lf",Val))
      la = (int) Val[0];
    else if(1 == sscanf(cLine,"Grid %lf",Val))
      gr = (int) Val[0];
    else if(3 == sscanf(cLine,"Offset %lf %lf %lf",Val,Val+1,Val+2)){
      xp = Val[0];yp = Val[1];zp = Val[2];
    }
    else if(3 == sscanf(cLine,"Rotate %lf %lf %lf",Val,Val+1,Val+2)){
      xa = Val[0];ya = Val[1];za = Val[2];
  }
    else if(3 == sscanf(cLine,"Sun %lf %lf %lf",Val,Val+1,Val+2)){
      xf = Val[0];yf = Val[1];zf = Val[2];
    }
    else if(3 == sscanf(cLine,"Info %lf %lf %lf",Val,Val+1,Val+2)){
      xi = Val[0];yi = Val[1];zi = Val[2];
    }
    else if(5 == sscanf(cLine,"Legend %lf %lf %lf %lf %lf",Val,Val+1,Val+2,Val+3,Val+4)){
      xLeg = Val[0];yLeg = Val[1];zLeg = Val[2];dxLeg = Val[3];dyLeg = Val[4];
    }
    else if(3 == sscanf(cLine,"Background %lf %lf %lf",Val,Val+1,Val+2)){
      Rback = Val[0];Gback = Val[1];Bback = Val[2];
    }
    else if(1 == sscanf(cLine,"Wheel %lf",Val)){
      zw = Val[0];
    }
    else if(3 == sscanf(cLine,"GridEdge %lf %lf %lf",Val,Val+1,Val+2)){
      GridEdge[0] = Val[0];GridEdge[1] = Val[1];GridEdge[2] = Val[2];
    }
  }
  fclose(ConfFile);
}
#ifdef USE_TIFF
#include <tiff.h>
#include <tiffio.h>
int Draw::OpenImage(const char *FileName){
  Shout("OpenImage");
//   pixel = (GLubyte *)realloc(pixel,4*width*height*sizeof(*pixel));
//   glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE,pixel);
//   return 0;
  TIFF *image;
  uint16 photo=0, bps=0, spp=0, fillorder=0,resUnit=0;
  uint32 ImWidth1=0,ImHeight1=0;
  tsize_t stripSize;
  int imageOffset=0, result=0;
  int stripMax=0, stripCount=0;
  float xRes=0,yRes=0;
  // char *buffer, tempbyte;
  uint32 *buffer=NULL;
  int bufferSize=0, count=0;
  int NPixel=0,Level=4;
  if((image = TIFFOpen(FileName, "r")) == NULL){
    fprintf(stderr, "Could not open incoming image\n");
    return 0;
  }
  TIFFGetField(image, TIFFTAG_IMAGELENGTH, &ImHeight1);
  TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &ImWidth1);
  TIFFGetField(image, TIFFTAG_PHOTOMETRIC, &photo);
  TIFFGetField(image, TIFFTAG_FILLORDER, &fillorder);
  TIFFGetField(image, TIFFTAG_XRESOLUTION, &xRes);
  TIFFGetField(image, TIFFTAG_YRESOLUTION, &yRes);
  TIFFGetField(image, TIFFTAG_RESOLUTIONUNIT, &resUnit);
  TIFFGetField(image, TIFFTAG_BITSPERSAMPLE, &bps);
  TIFFGetField(image, TIFFTAG_SAMPLESPERPIXEL, &spp);
  stripSize = TIFFStripSize (image);
  stripMax = TIFFNumberOfStrips (image);
  imageOffset = 0;
  bufferSize = TIFFNumberOfStrips (image) * stripSize;
  ImHeight = (GLuint)ImHeight1;
  ImWidth  = (GLuint) ImWidth1;
  NPixel = ImWidth*ImHeight*Level;
  //sprintf(info,"width %d height %d photo %d fillorder %d Res (%lf %lf) ResUnit %d bit/sample %d sample/pixel %d stripSize %d stripMax %d bufferSize %d NPixel %d\n",ImWidth,ImHeight,photo,fillorder,xRes,yRes,resUnit,bps,spp,stripSize,stripMax,bufferSize,NPixel);
  //buffer = (uint32 *) malloc(bufferSize);
  buffer = (uint32 *)_TIFFmalloc(ImHeight1 * ImWidth1 * sizeof (uint32));
  if(buffer == NULL){
    fprintf(stderr, "Could not allocate enough memory for the uncompressed image\n");
    return 0;
  }
  free(pixel);
  pixel = NULL;
  pixel = (GLubyte *)calloc(4*ImWidth*ImHeight,sizeof(*pixel));
  if(!TIFFReadRGBAImage(image,ImWidth1,ImHeight1,buffer,1)){
  //if(!TIFFReadRGBAImageOriented(image,ImWidth1,ImHeight1,buffer,ORIENTATION_TOPLEFT,0)){
    _TIFFfree(buffer);
    printf("The immage was not read\n");
    throw "Corrupted TIFF file!";
  }
  for(int h=0; h<ImHeight; h++) {
    for(int w=0; w<ImWidth; w++) {
      pixel[(h*ImWidth+w)*4+0] = buffer[(h*ImWidth+w)] & 0x000000ff;
      pixel[(h*ImWidth+w)*4+1] = ((buffer[(h*ImWidth+w)] >> 8  ) & 0x0000ff);
      pixel[(h*ImWidth+w)*4+2] = ((buffer[(h*ImWidth+w)] >> 16 ) & 0x00ff);
      pixel[(h*ImWidth+w)*4+3] = ((buffer[(h*ImWidth+w)] >> 24 ) & 0xff);
    }
  }
  TIFFClose(image);
  _TIFFfree(buffer);
   return 1;
}
int Draw::WritePixel(){
  Shout("WritePixel");
  char cImage[160];
  sprintf(cImage,"Image%05u.tif",Step);
  TIFF *tif = TIFFOpen( cImage, "w");
  if(tif) {
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, ImWidth);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, ImHeight);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 4);
  }
  TIFFWriteEncodedStrip(tif, 0, pixel, ImWidth * ImHeight *4);
  TIFFClose(tif);
  return 0;
}
int Draw::Picture(){
  Shout("Picture");
  //       mkOpenGLJPEGImage jpgScreen;
  //       jpgScreen.SetDebugMode(true);
  //       jpgScreen.GetOpenGLScreenImage(imageWidth, imageHeight);
  //       jpgScreen.SaveToFile("ScreenShot.jpg");
  //GLubyte pixel[Resx][Resy][3];
  //  GLubyte *pixel=new GLubyte [Resx * Resy];
  glutPostRedisplay();
  char cImage[160];
  sprintf(cImage,"Image%05u.tif",Step);
  TIFF *tif = TIFFOpen( cImage, "w");
  if(tif) {
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, WinWidth);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, WinHeight);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
    //TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 4);
  }
  char * data = (char *) calloc(WinWidth * WinHeight*4,sizeof(char));
  // 	pixel = (GLubyte *)realloc(pixel,3*width*height*sizeof(*pixel));
  // 	glPixelStorei(GL_UNPACK_ALIGNMENT,1);
  // 	glRasterPos3i(0,0,0);
  // 	//glutUseLayer(GLUT_NORMAL);
  // 	glReadPixels(-1,-1,width,height,GL_RGBA,GL_UNSIGNED_BYTE,pixel);
  //         glReadBuffer(GL_BACK);
  //         //glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0, WinWidth, WinHeight, GL_RGBA, GL_UNSIGNED_BYTE, data);
  char * MonaSoMare = (char *) calloc(WinWidth * WinHeight*4,sizeof(char));
  for(int l=0;l<4;l++)
    for(int w = 0;w<WinWidth;w++)
      for(int h=0;h<WinHeight;h++)
	MonaSoMare[(h*WinWidth+w)*4+l] = data[( (WinHeight-h)*WinWidth+w)*4 + l];
  TIFFWriteEncodedStrip(tif, 0,MonaSoMare, WinWidth * WinHeight *4);
  TIFFClose(tif);
  free(data);
  free(MonaSoMare);
  //       //UINT *pixels=new UINT[nWidth * nHeight];
  //glReadBuffer(GL_BACK_LEFT); 
  //       //glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,image);
  //GL_COLOR_INDEX, GL_STENCIL_INDEX, GL_DEPTH_COMPONENT, GL_RED, GL_GREEN, GL_BLUE, 
  //GL_ALPHA, GL_RGB, GL_RGBA, GL_LUMINANCE, and GL_LUMINANCE_ALPHA. 
  //GL_UNSIGNED_BYTE, GL_BYTE, GL_BITMAP, GL_UNSIGNED_SHORT, GL_SHORT, GL_UNSIGNED_INT, GL_INT, or GL_FLOAT.    
}
#else
int Draw::OpenImage(const char *FileName){
  printf("libtiff not supported\n");
}
int Draw::Picture(){
  printf("libtiff not supported\n");
}
int Draw::WritePixel(){
  printf("libtiff not supported\n");
}
#endif //USE_TIFF
#ifdef USE_PNG
#include <png.h>
#include <pngwriter.h>
int Draw::WritePngwriter(){
  char cImage[160];
  sprintf(cImage,"Image%05u.png",Step);
  pngwriter Image(ImWidth,ImHeight,1.0,cImage);
  double Norm = 1./256.;
  for(int w = 0;w<ImWidth;w++)
    for(int h=0;h<ImHeight;h++){
      double Color[3];
      for(int l=0;l<3;l++){
	Color[l] = Norm*pixel[(h*ImWidth+w)*4+l];
      }
      Image.plot(w,h,Color[0],Color[1],Color[2]);
    }
  Image.close();
}
int Draw::WritePng(){
  Shout("WritePng");
  int NLevel = 4;
  png_byte color_type;
  png_byte bit_depth;
  png_structp png_ptr;
  png_infop info_ptr;
  int number_of_passes;
  png_byte header[8];	// 8 is the maximum size that can be checked
  char cImage[160];
  sprintf(cImage,"Image%05u.png",Step);
  png_bytep *  row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * ImHeight);
  for (int y=0; y<ImHeight; y++)
    row_pointers[y] = (png_byte*) malloc(ImWidth);//info_ptr->rowbytes);
  for (int y=0; y<ImHeight; y++) {
    png_byte* row = row_pointers[y];
    for (int x=0; x<ImWidth; x++) {
      png_byte* ptr = &(row[x*4]);
      for(int l=0;l<NLevel;l++)
	ptr[l] = pixel[(y*ImWidth+x)*NLevel+l];
    }
  }
  /* initialize stuff */
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  SigErr(!png_ptr,"[write_png_file] png_create_write_struct failed");
  info_ptr = png_create_info_struct(png_ptr);
  SigErr(!info_ptr,"[write_png_file] png_create_info_struct failed");
  SigErr(setjmp(png_jmpbuf(png_ptr)),"[read_png_file] Error during init_io");
  png_set_sig_bytes(png_ptr, 8);
  //  png_read_info(png_ptr, info_ptr);
  color_type = 6;//info_ptr->color_type;
  bit_depth = 8;//info_ptr->bit_depth;
  //  number_of_passes = png_set_interlace_handling(png_ptr);
  //png_read_update_info(png_ptr, info_ptr);
  FILE *fp = fopen(cImage, "wb");
  SigErr(!fp,"[write_png_file] File %s could not be opened for writing",cImage);
  SigErr(setjmp(png_jmpbuf(png_ptr)),"[write_png_file] Error during init_io");
 /* setup libpng for using standard C fread() function
     with our FILE pointer */
  png_init_io(png_ptr, fp);
  /* convert index color images to RGB images */
   //  if (color_type == PNG_COLOR_TYPE_PALETTE)
   //png_set_palette_to_rgb (png_ptr);
 /* convert 1-2-4 bits grayscale images to 8 bits
     grayscale. */
//   if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
//     png_set_gray_1_2_4_to_8 (png_ptr);

//   if (png_get_valid (png_ptr, info_ptr, PNG_INFO_tRNS))
//     png_set_tRNS_to_alpha (png_ptr);
  /* write header */
  SigErr(setjmp(png_jmpbuf(png_ptr)),"[write_png_file] Error during writing header");
  png_set_IHDR(png_ptr, info_ptr, ImWidth, ImHeight,
	       bit_depth, color_type, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
  png_write_info(png_ptr, info_ptr);
  /* write bytes */
  SigErr(setjmp(png_jmpbuf(png_ptr)),"[write_png_file] Error during writing bytes");
  png_write_image(png_ptr, row_pointers);
  /* end write */
  SigErr(setjmp(png_jmpbuf(png_ptr)),"[write_png_file] Error during end of write");
  png_write_end(png_ptr, NULL);
  /* cleanup heap allocation */
  for (int y=0; y<ImHeight; y++) free(row_pointers[y]);
  free(row_pointers);
  fclose(fp);
  return 0;
}
#else //USE_PNG
int Draw::WritePngwriter(){
  printf("Pngwriter not installed\n");
}
int Draw::WritePng(){
  printf("Pngwriter not installed\n");
}
#endif

void Draw::ReadScript(){
  Shout("ReadScript");
  FILE *File2Read;
  if((File2Read = fopen("DrawScript.dr","r"))==0){
    return ;
  }
  glDeleteLists(ScriptList,1);
  ScriptList = glGenLists(1);
  glNewList(ScriptList,GL_COMPILE);
  double Expand[3] = {Edge[0]*InvScaleUn,Edge[1]*InvScaleUn,Edge[2]*InvScaleUn};
  double *From = (double *)calloc(3,sizeof(double));
  double *To = (double *)calloc(3,sizeof(double));
  double *Hue = (double *)calloc(5,sizeof(double));
  char *cLine = (char *)calloc(265,sizeof(char));
  for(int k=0;!(fgets(cLine,256,File2Read)==NULL);k++){
    if(strncmp(cLine,"Line",4)==0){
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf\n",Hue,Hue+1,Hue+2,Hue+3,Hue+4);
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf %lf %lf\n",From,From+1,From+2,To,To+1,To+2);
      glPushMatrix();
      glColor4f(Hue[0],Hue[1],Hue[2],Hue[3]);
      glLineWidth(Hue[4]);
      glBegin(GL_LINES);
      glNormal3f(0.,0.,1.);
      glVertex3d(From[0]*Expand[0],From[1]*Expand[1],From[2]*Expand[2]);
      glVertex3d(To[0]*Expand[0],To[1]*Expand[1],To[2]*Expand[2]);
      glEnd();
      glPopMatrix();
    }
    else if(strncmp(cLine,"Curve",4)==0){
      glPushMatrix();
      glBegin(GL_LINE_STRIP);
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf\n",Hue,Hue+1,Hue+2,Hue+3);
      glColor4f(Hue[0],Hue[1],Hue[2],Hue[3]);
      glNormal3f(0.,0.,1.);
      while(1==1) {
	fpos_t FPos;
	fgetpos(File2Read,&FPos);
	if(fgets(cLine,256,File2Read)==0) break;
	if(sscanf(cLine,"%lf %lf %lf\n",From,From+1,From+2)!=3){
	  fsetpos(File2Read,&FPos);
	  break;
	}
	glVertex3d(From[0]*Expand[0],From[1]*Expand[1],From[2]*Expand[2]);
      }
      glEnd();
      glPopMatrix();
    }
    else if(strncmp(cLine,"Polygon",4)==0){
      glPushMatrix();
      glBegin(GL_POLYGON);
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf\n",Hue,Hue+1,Hue+2,Hue+3);
      glColor4f(Hue[0],Hue[1],Hue[2],Hue[3]);
      glNormal3f(0.,0.,1.);
      while(1==1) {
	fpos_t FPos;
	fgetpos(File2Read,&FPos);
	if(fgets(cLine,256,File2Read)==0) break;
	if(sscanf(cLine,"%lf %lf %lf\n",From,From+1,From+2)!=3){
	  fsetpos(File2Read,&FPos);
	  break;
	}
	glVertex3d(From[0]*Expand[0],From[1]*Expand[1],From[2]*Expand[2]);
      }
      glEnd();
      glPopMatrix();
    }
    else if(strncmp(cLine,"Arrow",4)==0){
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf\n",Hue,Hue+1,Hue+2,Hue+3,Hue+4);
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf %lf %lf\n",From,From+1,From+2,To,To+1,To+2);
      glLineWidth(Hue[4]);
      glColor4f(Hue[0],Hue[1],Hue[2],Hue[3]);
      glPushMatrix();//Arrow
      double Dist[3];
      double Freccia[3] = {Expand[0],0.,0.};
      double Axis[3];
      for(int d=0;d<3;d++){
	Dist[d] = (To[d] - From[d])*Expand[d];
      }      
      double Length1 = sqrt( QUAD(Dist[0]) + QUAD(Dist[1]) + QUAD(Dist[2]) );
      double Length2 = sqrt( QUAD(Expand[0]) );
      for(int d=0;d<3;d++){
	Axis[d] = .5*(Dist[d]/Length1 + Freccia[d]/Length2);
      }
      glTranslated(From[0]*Expand[0],From[1]*Expand[1],From[2]*Expand[2]);
      glRotatef(180.,Axis[0],Axis[1],Axis[2]);
      glScaled(Length1/Length2,1.0,1.0);
      glNormal3f(Dist[0],Dist[1],Dist[2]);
      glCallList(Arrow);
      glPopMatrix();//Arrow
    }
    else if(strncmp(cLine,"Sphere",6)==0){
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf\n",Hue,Hue+1,Hue+2,Hue+3);
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf %lf %lf\n",From,From+1,From+2,To,To+1,To+2);
      glPushMatrix();
      glTranslated(From[0]*Expand[0],From[1]*Expand[1],From[2]*Expand[2]);
      glColor4f(Hue[0],Hue[1],Hue[2],Hue[3]);
      glutSolidSphere(To[0]*Expand[0],(int)To[1],(int)To[2]);
      glPopMatrix();
    }
     else if(strncmp(cLine,"Cylinder",7)==0){
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf\n",Hue,Hue+1,Hue+2,Hue+3);
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf %lf %lf\n",From,From+1,From+2,To,To+1,To+2);
      glPushMatrix();
      glTranslated(From[0]*Expand[0],From[1]*Expand[1],From[2]*Expand[2]);
      glColor4f(Hue[0],Hue[1],Hue[2],Hue[3]);
      glutSolidSphere(To[0]*Expand[0],(int)To[1],(int)To[2]);
      glPopMatrix();
    }
   else if(strncmp(cLine,"Text",4)==0){
      char *String = (char *)calloc(265,sizeof(char));
      fgets(cLine,256,File2Read);
      for(int c=0;c<strlen(cLine);c++)
	String[c] = cLine[c];
      //      sscanf(cLine,"%s\n",String);
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf\n",Hue,Hue+1,Hue+2,Hue+3);
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf %lf %lf\n",From,From+1,From+2,To,To+1,To+2);
      glPushMatrix();
      glNormal3f(0.,0.,-1.);
      glColor4f(Hue[0],Hue[1],Hue[2],Hue[3]);
      glRasterPos3f(From[0]*Expand[0],From[1]*Expand[1],From[2]*Expand[2]);
      for(int i=0;i<strlen(String);i++)
	glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24,String[i]);
      glCallList(Point);
      glPopMatrix();
      free(String);
    }
    else if(strncmp(cLine,"Fog",3)==0){
      GLuint filter;
      GLuint fogMode[]= { GL_EXP, GL_EXP2, GL_LINEAR };
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf\n",Hue,Hue+1,Hue+2,Hue+3);
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf %lf %lf\n",From,From+1,From+2,To,To+1,To+2);
      GLuint fogfilter = From[0] > 2. ? 0 : (int) From[0];
      GLfloat fogColor[4]= {Hue[0], Hue[1], Hue[2], Hue[3]};
      glFogi(GL_FOG_MODE, fogMode[fogfilter]);
      glFogfv(GL_FOG_COLOR, fogColor);
      glFogf(GL_FOG_DENSITY, From[1]);
      glHint(GL_FOG_HINT, GL_DONT_CARE);
      glFogf(GL_FOG_START, From[2]*Expand[2]);
      glFogf(GL_FOG_END, To[2]*Expand[2]);
      glEnable(GL_FOG);
    }
    else if(strncmp(cLine,"Background",10)==0){
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf\n",Hue,Hue+1,Hue+2,Hue+3);
      glClearColor(Hue[0],Hue[1],Hue[2],Hue[3]);
    }
  }
  glEndList();
  fclose(File2Read);
  free(From);
  free(To);
  free(Hue);
  free(cLine);
}
#endif
