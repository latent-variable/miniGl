/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstring>
#include <iostream>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h


MGLpoly_mode Mode;
MGLmatrix_mode MMode;
//Global color
vec3 color;
mat4 matrix;
//list of vertex
vector<vertex> vertices;


//list of triangles
vector<Triangle> Triangles;


/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}


/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
  memset(data, 0, sizeof(MGLpixel)*width * height);
  #if 0
  for(unsigned i = 0; i < width; i++){
    for(unsigned j = 0; j< height; j++){
      data[i+j*width] = Make_Pixel(color[0]*255,color[1]*255,color[2]*255);
    }
  }
  #endif

  for(unsigned i = 0; i < Triangles.size();i++){
    Rasterize_Triangle(Triangles[i],width, height,data);
  }
  Triangles.clear();
}

float triangle_area(vec2 a, vec2 b, vec2 c){
    return (a[0]*(b[1]-c[1])+ a[1]*(c[0] - b[0]) + (b[0]*c[1] - b[1]*c[0]));
}

void Rasterize_Triangle(const Triangle &tri,
                        int width,
                        int height,
                        MGLpixel* data)
{
  vec2 a, b, c;
  a[0] = (tri.a.pos[0] + 1)*width/2 -0.5;
  a[1] = (tri.a.pos[1] + 1)*height/2 - 0.5;

  b[0] = (tri.b.pos[0] + 1)*width/2 -0.5;
  b[1] = (tri.b.pos[1] + 1)*height/2 - 0.5;

  c[0] = (tri.c.pos[0] + 1)*width/2 -0.5;
  c[1] = (tri.c.pos[1] + 1)*height/2 - 0.5;

  float area = triangle_area(a,b,c);

  for(int i = 0; i < width; i++){
    for(int j = 0; j< height; j++){
      vec2 p;
      p[0] = i;
      p[1] = j;

      float alpha = triangle_area(p,b,c)/area;
      float beta = triangle_area(a,p,c)/area;
      float gamma = triangle_area(a,b,p)/area;

      if(alpha >= 0 && beta >= 0 && gamma >=0 ){
        data[i+j*width] = Make_Pixel(tri.a.color[0]*255,tri.a.color[1]*255,tri.a.color[2]*255);
      }

    }
  }


}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{

  Mode = mode;
}

/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{


  if(Mode == MGL_TRIANGLES){

      for(unsigned i = 0; i < vertices.size(); i+=3){
         Triangle triangle;

         triangle.a = vertices[i];
         triangle.b = vertices[i+1];
         triangle.c = vertices[i+2];

         Triangles.push_back(triangle);

         if( i + 6 > vertices.size()) break;
      }
  }
  if(Mode == MGL_QUADS){

      for(unsigned i = 0; i < vertices.size(); i+=4){
         Triangle triangle1;
         Triangle triangle2;

         triangle1.a = vertices[i];
         triangle2.a = vertices[i];
         triangle1.b = vertices[i+3];
         triangle2.b = vertices[i+1];
         triangle1.c = vertices[i+2];
         triangle2.c = vertices[i+2];

         Triangles.push_back(triangle2);
         Triangles.push_back(triangle1);


         if( i + 8 > vertices.size()) break;
      }
  }
  vertices.clear();

}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
  mglVertex3(x,y,0);

}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
  vec4 tempv;

  tempv[0] = x;
  tempv[1] = y;
  tempv[2] = z;
  tempv[3] = 1.0f;

  //cout <<"tempv: "<<tempv<<endl;
  tempv = matrix * tempv;
  //cout <<"Matrix "<< matrix << endl;
  //cout <<"tempv3: "<< tempv[3]<<endl;
  vertex temp(tempv/tempv[3],color);

 //cout << "temp.pos " << temp.pos<<endl;
  vertices.push_back(temp);

}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
  MMode = mode;


}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
    if(MMode == MGL_PROJECTION)
      matrix = {{1.0f,0.0f,0.0f,0.0f,
                 0.0f,1.0f,0.0f,0.0f,
                 0.0f,0.0f,1.0f,0.0f,
                 0.0f,0.0f,0.0f,1.0f}};

}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{

}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
    if(MMode == MGL_PROJECTION){
        float A = (right + left)/(right -left);
        float B = (top + bottom)/(top-bottom);
        float C = -(far + near)/(far - near);
        float D = - 2*(far*near)/(far - near);

        mat4 temp = {{2*near/(right-left),0,0,0,
                      0, 2*near/(top-bottom),0,0,
                      A,B,C,-1,
                      0,0,D,0}};

        matrix = temp*matrix;

    }
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{

  if (MMode == MGL_PROJECTION){
    float tx = -(right+left)/(right-left);
    float ty = -(top+bottom)/(top-bottom);
    float tz = -(far + near)/(far - near);

    /*mat4 temp = {{2.0f/(right-left),0.0f,0.0f,tx,
                 0.0f,2.0f/(top-bottom),0.0f, ty,
                 0.0f,0.0f, -2.0f/(far-near), tz,
                 0.0f,0.0f,0.0f,1.0f}};*/

      mat4 temp = {{2.0f/(right-left),0.0f,0.0f,0.0f,
                  0.0f,2.0f/(top-bottom),0.0f,0.0f,
                  0.0f,0.0f,-2.0f/(far-near),0.0f,
                  tx,ty,tz,1.0f}};

      matrix = temp*matrix;
  }


}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
  color[0] = red;
  color[1] = green;
  color[2] = blue;
}
