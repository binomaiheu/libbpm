/**
   @file
   @ingroup orbit
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_orbit.h>
#include <bpm/bpm_interface.h>

int get_bpmhits( beamconf_t *beam, bpmconf_t *bpm ) {
  
  int i = 0;

  if ( ! beam ) {
    bpm_error( "Invalid pointer arguments in get_bpmhits(...)",
               __FILE__, __LINE__ ); 
    return BPM_FAILURE;
  }

  if ( ! bpm ) {
    bpm_error( "Invalid pointer arguments in get_bpmhits(...)",
               __FILE__, __LINE__ ); 
    return BPM_FAILURE;
  }

  for ( i = 0; i < beam->nbunches; i++ ) {
    if ( get_bpmhit( &(beam->bunch[i]), bpm ) == BPM_FAILURE ) return BPM_FAILURE;
  }
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------

int get_bpmhit( bunchconf_t *bunch, bpmconf_t *bpm ) {

  int i = 0;
  
  struct v3 nl;     /* direction of the line */
  struct v3 pl;     /* point on the line */
  struct v3 np;     /* normal of the plane */
  struct v3 pp;     /* point in the plane */
  struct v3 xp;     /* coordinate system x unit vector */
  struct v3 yp;     /* coordinate system y unit vector */
  struct v3 zp;     /* coordinate system z unit vector */
  struct m33 rotp;  /* rotation matrix in yaw,roll,pitch of BPM */
  double lambda;    /* parametric line / plane intersection solution */
  struct v3 vtemp;  /* temporary vector */
  struct v3 v;      /* solution */
  struct v3 lv;     /* local BPM coordinates solution */

  if ( ! bunch ) {
    bpm_error( "Invalid pointer arguments in get_bpmhit(...)",
               __FILE__, __LINE__ ); 
    return BPM_FAILURE;
  }

  if ( ! bpm ) {
    bpm_error( "Invalid pointer arguments in get_bpmhit(...)",
               __FILE__, __LINE__ ); 
    return BPM_FAILURE;
  }

  /* form BPM detector plane */
  pp.x = bpm->geom_pos[0]; 
  pp.y = bpm->geom_pos[1]; 
  pp.z = bpm->geom_pos[2]; 
  xp.x = 1.0; xp.y = 0;   xp.z = 0;
  yp.x = 0;   yp.y = 1.0; yp.z = 0;
  zp.x = 0;   zp.y = 0;   zp.z = 1.0;

  m_rotmat(&rotp, bpm->geom_tilt[0], bpm->geom_tilt[1], bpm->geom_tilt[2]);
  v_matmult(&rotp,&xp);
  v_matmult(&rotp,&yp);
  v_matmult(&rotp,&zp);
  
  v_copy(&np,&xp);
  v_cross(&np,&yp);
  
  /* find intersection of bunch with BPM */
  pl.x = bunch->position[0];
  pl.y = bunch->position[1];
  pl.z = bpm->geom_pos[2];
  
  nl.x = sin(bunch->slope[0]) * cos(bunch->slope[1]);
  nl.y = sin(bunch->slope[0]) * sin(bunch->slope[1]);
  nl.z = cos(bunch->slope[0]);
  
  /* intersection solution */    
  v_copy(&vtemp,&pp);
  v_sub(&vtemp,&pl);
  lambda = v_dot(&vtemp,&np)/v_dot(&np,&nl);
  
  v_copy(&v,&nl);
  v_scale(&v,lambda);
  v_add(&v,&pl);
  
  /* beam position in coordinate system of bpm */
  v_copy(&lv,&v);
  v_sub(&lv,&pp);
  
  bunch->bpmposition[0]   = v_dot(&lv,&xp);
  bunch->bpmposition[1]   = v_dot(&lv,&yp);
  bunch->bpmposition[2]   = v.z; // in global coo !!!!
  
  bunch->bpmslope[0] = bunch->slope[0] - bpm->geom_tilt[0];
  bunch->bpmslope[1] = bunch->slope[1] - bpm->geom_tilt[1];
  
  bunch->bpmtilt[0]  = 0.; // bunch tilt bunch->tilt[0], bunch->tilt[1]...
  bunch->bpmtilt[1]  = 0.; // bunch tilt...


  return BPM_SUCCESS;
}
