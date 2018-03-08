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

//******************************************************
//vertex structure
struct vertex
{

  vertex(){}
  vertex(vec4 p,vec3 c):  color(c),pos(p){}
  vec3 color;
  vec4 pos;
};

//structure that stores the 3 vertices of a triangle
struct Triangle
{
  vertex a,b,c;
};

//Helper fuction for mglReadPixels
void Rasterize_Triangle(const Triangle &tri,
                       int width,
                       int height,
                       MGLpixel* data);
//******************************************************
MGLpoly_mode Mode;
MGLmatrix_mode MMode;
//Global color
vec3 color;
//list of vertex
vector<vertex> vertices;
vector<MGLfloat> z_buffering;

//stack of matrix
mat4 I = {{1.0f,0.0f,0.0f,0.0f,
            0.0f,1.0f,0.0f,0.0f,
            0.0f,0.0f,1.0f,0.0f,
            0.0f,0.0f,0.0f,1.0f}};
vector<mat4> Matrix_stack_projection = {I};


vector<mat4> Matrix_stack_modelview = {I};

//list of triangles
vector<Triangle> Triangles;


/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}
//Helper function
/*
mat4& current_matrix(){
  if(MMode == MGL_PROJECTION){
    return projection_matrix;
  }else{
    return modelview_matrix;
  }
}*/

//helper fuction
mat4& top_of_active_matrix_stack(){
  if(MMode == MGL_PROJECTION){
    return Matrix_stack_projection.back();
  }else{
    return Matrix_stack_modelview.back();
  }
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
  z_buffering.assign(width*height,100000000);

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
  //z-buffering
  MGLfloat* zptr = z_buffering.data();

  MGLfloat Az, Bz, Cz, wa, wb, wc;
  wa = tri.a.pos[3];
  wb = tri.b.pos[3];
  wc = tri.c.pos[3];
  Az = tri.a.pos[2]/abs(wa);
  Bz = tri.b.pos[2]/abs(wb);
  Cz = tri.c.pos[2]/abs(wc);

  vec2 a, b, c;
  a[0] = (tri.a.pos[0] + 1)*width/2 -0.5;
  a[1] = (tri.a.pos[1] + 1)*height/2 - 0.5;

  b[0] = (tri.b.pos[0] + 1)*width/2 -0.5;
  b[1] = (tri.b.pos[1] + 1)*height/2 - 0.5;

  c[0] = (tri.c.pos[0] + 1)*width/2 -0.5;
  c[1] = (tri.c.pos[1] + 1)*height/2 - 0.5;

  // cout << "a[0] ="<< a[0]<< " a[1] ="<< a[1]<<endl;
  // cout << "b[0] ="<< b[0]<< " b[1] ="<< b[1]<<endl;
  // cout << "c[0] ="<< c[0]<< " c[1] ="<< c[1]<<endl;

  int mini, minj, maxi, maxj;

  mini = min(a[0], min(b[0],c[0]));
  mini = max(min(width,mini - 2), 0);

  minj = min(a[1], min(b[1],c[1]));
  minj = max(min(height,minj - 2), 0);

  maxi = max(a[0], max(b[0],c[0]));
  maxi = min(width , max(maxi + 2, 0));

  maxj = max(a[1], max(b[1],c[1]));
  maxj = max(0 , min(maxj + 2, height));

  // cout<< "mini " << mini << " maxi " << maxi << endl;
  //
  // cout<< "minj " << minj << " maxj " << maxj << endl;
  if(MMode == MGL_MODELVIEW){
    //cout << "Here"<<endl;
    mini = 0;
    maxi = width;
    minj = 0;
    maxj = height;
  }

  float area = triangle_area(a,b,c);
  for(int j = minj; j< maxj; ++j){
    for(int i = mini; i < maxi; ++i,++zptr ){
  // for(int j = 0; j< width; ++j){
  //   for(int i = 0; i < height; ++i,++zptr ){

      vec2 p;
      p[0] = i;
      p[1] = j;

      float alpha = triangle_area(p,b,c)/area;
      float beta = triangle_area(a,p,c)/area;
      float gamma = triangle_area(a,b,p)/area;



      float z_current = Az*alpha + Bz*beta+ Cz*gamma;

      //z_current = std::clamp(z_current,-1.0f, 1.0f);
      if(alpha >= 0 && beta >= 0 && gamma >=0 ){
        if(*zptr > z_current){
            if(z_current >= -1.0f && z_current <= 1.0f ){

                float K = 1 /(alpha *wa + beta * wb + gamma*wc);
                float alpha_p = alpha * wa * K;
                float beta_p = beta *wb *K;
                float gamma_p = gamma * wc * K;
                K =  1 / (alpha_p/ wa + beta_p/wb + gamma_p/wc);
                alpha = (alpha_p / wa)* K;
                beta = (beta_p / wb) * K ;
                gamma = (gamma_p/wc) * K;
                float percent_red = tri.a.color[0]*255*alpha + tri.b.color[0]*255*beta +tri.c.color[0]*255*gamma;
                float percent_green = tri.a.color[1]*255*alpha + tri.b.color[1]*255*beta +tri.c.color[1]*255*gamma;
                float percent_blue = tri.a.color[2]*255*alpha + tri.b.color[2]*255*beta +tri.c.color[2]*255*gamma;
                data[i+j*width] = Make_Pixel(percent_red,percent_green,percent_blue);
                *zptr = z_current;
          }
        }
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
  //tempv = projection_matrix* modelview_matrix* tempv;
  tempv = Matrix_stack_projection.back()* Matrix_stack_modelview.back()* tempv;
  //cout <<"Matrix "<< cmat << endl;
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
  if(MMode == MGL_PROJECTION){
    if(Matrix_stack_projection.empty() == false){
       //cout << " Here push projection \n";
        Matrix_stack_projection.push_back(Matrix_stack_projection.back());
    }
  }else{
    if(Matrix_stack_modelview.empty() == false){
        //cout << " Here push modelview \n";
        Matrix_stack_modelview.push_back(Matrix_stack_modelview.back());
    }
  }
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
  if(MMode == MGL_PROJECTION){
    if(Matrix_stack_projection.empty() == false){
      //cout << " Here pop projection \n";
        Matrix_stack_projection.pop_back();
    }
  }else{
    if(Matrix_stack_modelview.empty() == false){
        //cout << " Here pop modelview \n";
        Matrix_stack_modelview.pop_back();
    }
  }
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
   mat4 I = {{1.0f,0.0f,0.0f,0.0f,
             0.0f,1.0f,0.0f,0.0f,
             0.0f,0.0f,1.0f,0.0f,
             0.0f,0.0f,0.0f,1.0f}};

  mat4& Mstk = top_of_active_matrix_stack();
  Mstk = I;
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
void mglLoadMatrix(const MGLfloat *m)
{
  mat4& Mstk = top_of_active_matrix_stack();

  Mstk = {{m[0],m[1],m[2],m[3],
          m[4],m[5],m[6],m[7],
          m[8],m[9],m[10],m[11],
          m[12],m[13],m[14],m[15]}};

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
void mglMultMatrix(const MGLfloat *m)
{
  mat4& Mstk = top_of_active_matrix_stack();

  mat4 temp = {{m[0],m[1],m[2],m[3],
                m[4],m[5],m[6],m[7],
                m[8],m[9],m[10],m[11],
                m[12],m[13],m[14],m[15]}};



  Mstk = Mstk*temp;
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
    mat4& Mstk = top_of_active_matrix_stack();

    mat4 trans = {{1, 0, 0, 0,
                   0, 1, 0, 0,
                   0, 0, 1, 0,
                   x, y, z, 1}};
    Mstk = Mstk * trans;
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
    float A = angle *M_PI /180;
    mat4& Mstk = top_of_active_matrix_stack();

    float norm = sqrt(x*x + y*y + z*z);
    float c =  cos(A);
    float s = sin(A);
    x = x/norm;
    y = y/norm;
    z = z/norm;



    mat4 rotate = {{ x*x*(1-c) +c, y*x*(1-c)+z*s,x*z*(1-c)-y*s, 0,
                     x*y*(1-c)-z*s, y*y*(1-c)+c, y*z*(1-c)+x*s, 0,
                     x*z*(1-c)+y*z, y*z*(1-c)-x*s, z*z*(1-c)+c, 0,
                     0, 0, 0, 1}};

   Mstk = Mstk * rotate;

}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
  mat4& Mstk = top_of_active_matrix_stack();

  mat4 scale = {{x, 0, 0 ,0,
                0, y, 0, 0,
                0, 0, z, 0,
                0, 0, 0, 1}};

  Mstk = Mstk *scale;
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

        float A = (right + left)/(right -left);
        float B = (top + bottom)/(top-bottom);
        float C = -(far + near)/(far - near);
        float D = - 2*(far*near)/(far - near);
        mat4& Mstk = top_of_active_matrix_stack();

        mat4 temp = {{2*near/(right-left),0,0,0,
                      0, 2*near/(top-bottom),0,0,
                      A,B,C,-1,
                      0,0,D,0}};

        Mstk =  Mstk * temp;


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

    float tx = -(right+left)/(right-left);
    float ty = -(top+bottom)/(top-bottom);
    float tz = -(far + near)/(far - near);

    /*mat4 temp = {{2.0f/(right-left),0.0f,0.0f,tx,
                 0.0f,2.0f/(top-bottom),0.0f, ty,
                 0.0f,0.0f, -2.0f/(far-near), tz,
                 0.0f,0.0f,0.0f,1.0f}};*/

      mat4& Mstk = top_of_active_matrix_stack();

      mat4 temp = {{2.0f/(right-left),0.0f,0.0f,0.0f,
                  0.0f,2.0f/(top-bottom),0.0f,0.0f,
                  0.0f,0.0f,-2.0f/(far-near),0.0f,
                  tx,ty,tz,1.0f}};

     Mstk = Mstk*temp;

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
