//
// Torbert, 3 February 2014
//
// gcc -lm raytrace.c
// ./a.out
// display [filename].ppm



#include <stdio.h>
#include <math.h>
//Window Properties
#define WIDTH 640
#define HIGHT 480
//System Properties
int shading = 3;
int colorfulFlooring = 1;
int floorpatternsize = 7;
#define numObjects 3
int reflecting = 1;
double reflectence = .5;
//
//Object Properties
//Screen
double sz = 0;
//Floor
double fy = .333333;
//Color
typedef struct
{
  int r;
  int g;
  int b;
}color;
color rgb[HIGHT][WIDTH] = {0};
color colorFormat[numObjects + 2] = {255,0,0,0,255,0,0,0,255,0,0,0,255,255,255};
//Eye
typedef struct
{
  double x;
  double y;
  double z;
} eye;
eye e = {.5,.5,-1};
//Light
typedef struct
{
  double x;
  double y;
  double z;
} light;
light l = {0,1,-.5};
//Ray
typedef struct
{
  double x;
  double y;
  double z;
} Ray;
//Spheres
typedef struct
{
  double x;
  double y;
  double z;
  double r;
} sphere;
sphere spheres[3] = {0.333333,0.666667,0.666667,0.333333,0.833333,0.500000,0.500000,0.166667,0.500000,0.500000,0.166667,0.166667};
//
Ray r = {0,0,0};
double minT,a,b,c;
//Misc. Structs
typedef struct
{
  int object;
  double time;
} objectAndTime;
//
double findMin(double T[3])
{
  double min = T[0];
  int x;
  for(x = 0; x < numObjects;x++)
  {
    if(min > T[x]) min = T[x];
  }
  return min;
}
int getFloorPattern(double xpos, double zpos)
{
  int floorpattern = abs((((int)((xpos)*floorpatternsize))%2 + ((int)((zpos)*floorpatternsize))%2)%2) + numObjects;
  return floorpattern;
}
objectAndTime getMinTime(double temp[3],Ray ray, int currentObject)
{
  int objects;
  double a,b,c;
  a = ray.x*ray.x + ray.y*ray.y + ray.z*ray.z;
  double T[3];
  for(objects = 0; objects < 3; objects++)
  {
    b = 2*(ray.x*(temp[0] - spheres[objects].x) + ray.y*(temp[1]-spheres[objects].y) + ray.z*(temp[2] - spheres[objects].z));
    c = (temp[0]-spheres[objects].x)*(temp[0]-spheres[objects].x) + (temp[1] -spheres[objects].y)*(temp[1] -spheres[objects].y) + (temp[2] -spheres[objects].z)*(temp[2] -spheres[objects].z) - spheres[objects].r*spheres[objects].r;
    if((b*b - 4*a*c) < 0 || objects == currentObject)//miss
    {
      T[objects] = 9999999;
    }
    else
    {
      double timeValue = (-1*b - sqrt(b*b - 4*a*c))/(2*a);
      if((b*b - 4*a*c) < 0 || timeValue <= 0)
      {
        T[objects] = 9999999;
      }
      else
      {
        T[objects] = timeValue;
      }
    }
  }
  minT = findMin(T);
  if(minT == 9999999)
  {
   
    if(ray.y < 0)
    {
      minT = (fy - temp[1])/ray.y;
      int floorType = getFloorPattern(temp[0]+ray.x*minT,temp[2]+ray.z*minT);
      objectAndTime oT = {floorType,minT};
      return oT;
    }
    objectAndTime oT = {-1,minT};
    return oT;
  }
  for(objects = 0; objects < 3; objects++)
  {
    if(minT == T[objects])
    {
      objectAndTime oT = {objects,minT};
      return oT;
    }
  }
}
objectAndTime getObject(int y, int x)
{
  r.x = ((double)x/WIDTH)  - e.x;
  r.y = ((double)(HIGHT -y)/HIGHT)  - e.y;
  r.z = sz - e.z;
  //
  double eye[3] = {e.x,e.y,e.z};
  objectAndTime oT = getMinTime(eye,r,9999999);
  return oT;
  
}

double getIntensityFloor(void)
{
  double intensity = .5;
  double time = (fy - e.y)/r.y;
  if(shading > 0)
  {
    double T[numObjects];   
    double ix = (e.x + r.x*time);
    double iy = (e.y + r.y*time);
    double iz = (e.z + r.z*time);
    Ray rS = {ix-l.x, iy-l.y, iz-l.z};
    double light[3] = {l.x,l.y,l.z};
    objectAndTime oT = getMinTime(light,rS,9999999);
    if(oT.object == 4)
    {
      intensity = 1;
      if( shading > 1)
      {
        Ray normal = {0,1,0};
        double normalFloor = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
        double normalLight = sqrt(rS.x*rS.x + rS.y*rS.y + rS.z*rS.z);
        double dotProd = -1*(normal.x*(rS.x) + normal.y*(rS.y) + normal.z*(rS.z));     
        if(dotProd <= 0) 
        {
          dotProd = 0;
          intensity = .5;
        }
        if(shading > 2)
        {
          intensity = .5 + .5*dotProd/(normalFloor*normalLight);
        }
      }
    }
  }
  return intensity;
}
double getIntensitySphere(int color)
{
  double intensity = 1;
  if(shading > 0)
  {
    double T[numObjects];
    double ix = (e.x + r.x*minT);
    double iy = (e.y + r.y*minT);
    double iz = (e.z + r.z*minT);
    Ray rS = {ix-l.x, iy-l.y, iz-l.z};
    double light[3] = {l.x,l.y,l.z};
    objectAndTime oT = getMinTime(light,rS,9999999);
    if(oT.object != color)
    {
      intensity = .5; 
    }
    else if(shading > 1)
    {
      Ray normal = {(ix - spheres[oT.object].x)/spheres[oT.object].r,(iy - spheres[oT.object].y)/spheres[oT.object].r,(iz - spheres[oT.object].z)/spheres[oT.object].r};
      double normalCircle = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
      double normalLight = sqrt(rS.x*rS.x + rS.y*rS.y + rS.z*rS.z);
      double dotProd = -1*(normal.x*(rS.x) + normal.y*(rS.y) + normal.z*(rS.z));
      if(dotProd <= 0) 
      {
        dotProd = 0;
        intensity = .5;
      }
      if(shading > 2)
      {
        if(intensity != .5) intensity = .5 + .5*dotProd/(normalCircle*normalLight);
      }
    }
  }
  return intensity;
}
double getReflectedSpheres(int object,double time, int y,int x)
{
  double T[numObjects];
  double ix = (e.x + r.x*time);
  double iy = (e.y + r.y*time);
  double iz = (e.z + r.z*time);
  Ray normal = {(ix - spheres[object].x)/spheres[object].r,(iy - spheres[object].y)/spheres[object].r,(iz - spheres[object].z)/spheres[object].r};
  double normalCircle = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
  double normalEye = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
  double dotProd = (normal.x*(r.x) + normal.y*(r.y) + normal.z*(r.z));
  //
  double wx = r.x - 2*dotProd*normal.x;
  double wy = r.y - 2*dotProd*normal.y;
  double wz = r.z - 2*dotProd*normal.z;
  //
  Ray w = {wx,wy,wz};
  double i[3] = {ix,iy,iz};
  objectAndTime oT = getMinTime(i,w,object);
   //
  return oT.object;
}
int main(void)
{
  int y , x ;
  objectAndTime oT;
  for( y = 0 ; y < HIGHT ; y++ )
  {
      for( x = 0 ; x < WIDTH ; x++)
      {
        oT = getObject(y,x);
        if(oT.object == numObjects || oT.object == numObjects + 1)//Found the floor
        {
          double intensity = getIntensityFloor();
          //
          rgb[y][x].r = colorFormat[oT.object].r*intensity;
          rgb[y][x].g = colorFormat[oT.object].g*intensity;
          rgb[y][x].b = colorFormat[oT.object].b*intensity;
        }
        else if( oT.object >= 0)//Found a Circle
        {
          double intensity = getIntensitySphere(oT.object);
          double reflection = 0;
          if( reflecting > 0) reflection = reflectence;
          rgb[y][x].r = colorFormat[oT.object].r*(1-reflection);
          rgb[y][x].g = colorFormat[oT.object].g*(1-reflection);
          rgb[y][x].b = colorFormat[oT.object].b*(1-reflection);
          //
          if( reflecting > 0) 
          {
            int reflectedObject = getReflectedSpheres(oT.object,oT.time,y,x);
            //
            if(reflectedObject == -1)
            {
              if(oT.object == 0) rgb[y][x].r = 255*intensity;
              else if(oT.object == 1) rgb[y][x].g = 255*intensity;
              else if(oT.object == 2) rgb[y][x].b = 255*intensity;
            }
            else
            {
              rgb[y][x].r = (rgb[y][x].r + colorFormat[reflectedObject].r*reflectence);
              rgb[y][x].g = (rgb[y][x].g + colorFormat[reflectedObject].g*reflectence);
              rgb[y][x].b = (rgb[y][x].b + colorFormat[reflectedObject].b*reflectence);
            }            
          }
          rgb[y][x].r = rgb[y][x].r*intensity;
          rgb[y][x].g = rgb[y][x].g*intensity;
          rgb[y][x].b = rgb[y][x].b*intensity;
        }
      }
  }
  FILE* fout ;
  fout = fopen( "raytraceReflection.ppm" , "w" ) ;
  fprintf( fout , "P3\n" ) ;
  fprintf( fout , "%d %d\n" , WIDTH , HIGHT ) ;
  fprintf( fout , "255\n" ) ;
  for( y = 0 ; y < HIGHT ; y++ )
  {
    for( x = 0 ; x < WIDTH ; x++)
    {
        fprintf( fout , "%d %d %d\n" ,
        rgb[y][x].r , rgb[y][x].g , rgb[y][x].b ) ;
    }
  }
  close( fout ) ;
  return 0 ;
}
//
// end of file
//
