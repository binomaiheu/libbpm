/**
   @file
   @ingroup orbit
*/

#include <bpm/bpm_orbit.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void v_copy(struct v3 *v1, struct v3 *v2) {
  v1->x = v2->x;
  v1->y = v2->y;
  v1->z = v2->z;
}


double v_mag(struct v3 *v1) {
  return sqrt(v_dot(v1,v1));
}

void   v_scale(struct v3 *v1, double dscale) {
  v1->x = v1->x*dscale;
  v1->y = v1->y*dscale;
  v1->z = v1->z*dscale;
}

void   v_norm(struct v3 *v1) {
  v_scale(v1,1/v_mag(v1));
}

void   v_matmult(struct m33 *m1, struct v3 *v1) {
  double x = 0;
  double y = 0;
  double z = 0;
  x = m1->e[0][0]*v1->x + m1->e[0][1]*v1->y + m1->e[0][2]*v1->z;
  y = m1->e[1][0]*v1->x + m1->e[1][1]*v1->y + m1->e[1][2]*v1->z;
  z = m1->e[2][0]*v1->x + m1->e[2][1]*v1->y + m1->e[2][2]*v1->z;
  v1->x = x;
  v1->y = y;
  v1->z = z;
}

void   v_add(struct v3 *v1, struct v3 *v2) {
  v1->x = v1->x + v2->x;
  v1->y = v1->y + v2->y;
  v1->z = v1->z + v2->z;
}

void   v_sub(struct v3 *v1, struct v3 *v2) {
  v1->x = v1->x - v2->x;
  v1->y = v1->y - v2->y;
  v1->z = v1->z - v2->z;
}

double   v_dot(struct v3 *v1, struct v3 *v2) {
  return v1->x*v2->x + v1->y*v2->y + v1->z*v2->z;
}

void   v_cross(struct v3 *v1, struct v3 *v2) {
  double x = 0;
  double y = 0;
  double z = 0;

  x =  v1->y*v2->z - v1->z*v2->y;
  y = -v1->x*v2->z + v1->z*v2->x;
  z =  v1->x*v2->y - v1->y*v2->x;

  v1->x = x;
  v1->y = y;
  v1->z = z;
}

void v_print(struct v3 *v1) {
  printf("(%f,%f,%f)\n",v1->x,v1->y,v1->z);
}

void   m_rotmat(struct m33 *m1, double alpha, double beta, double gamma) {
  struct m33 *xrot;
  struct m33 *yrot; 
  struct m33 *zrot;
  struct m33 *mtem;

  /*  printf("m_rotmat> %f %f %f\n",alpha,beta,gamma); */

  xrot = calloc(1,sizeof(struct m33));
  yrot = calloc(1,sizeof(struct m33));
  zrot = calloc(1,sizeof(struct m33));
  mtem = calloc(1,sizeof(struct m33));

  xrot->e[0][0] = 1.0;
  xrot->e[1][1] = cos(alpha);
  xrot->e[1][2] = sin(alpha);
  xrot->e[2][1] = -sin(alpha);
  xrot->e[2][2] = cos(alpha);
  
  yrot->e[0][0] = cos(beta);
  yrot->e[0][2] = -sin(beta);
  yrot->e[1][1] = 1.0;
  yrot->e[2][0] = sin(beta);
  yrot->e[2][2] = cos(beta);

  zrot->e[0][0] = cos(gamma);
  zrot->e[0][1] = sin(gamma);
  zrot->e[1][0] = -sin(gamma);
  zrot->e[1][1] = cos(gamma);
  zrot->e[2][2] = 1.0;
  
  /*  m_print(xrot);
  m_print(yrot);
  m_print(zrot); */
  
  m_matmult(mtem,yrot,xrot);
  m_matmult(m1,zrot,mtem);  

  
  /* m_print(m1); */

  free(xrot);
  free(yrot);
  free(zrot);
  free(mtem);

}

void   m_matmult(struct m33 *m,struct m33 *m1, struct m33 *m2) {
  int i;
  int j;
  int k;
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      m->e[i][j] = 0.0;
      for(k=0;k<3;k++) {
	m->e[i][j] += m1->e[i][k]*m2->e[k][j];
      }
    }
  }
}

void   m_matadd(struct m33 *m1, struct m33 *m2) 
{
  int i;
  int j;
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      m1->e[i][j] = m1->e[i][j]+m2->e[i][j];
    }
  }
}

void m_print(struct m33 *m1) {
  printf("(%f %f %f\n",m1->e[0][0],m1->e[0][1],m1->e[0][2]);
  printf(" %f %f %f\n",m1->e[1][0],m1->e[1][1],m1->e[1][2]);
  printf(" %f %f %f)\n",m1->e[2][0],m1->e[2][1],m1->e[2][2]);
}
