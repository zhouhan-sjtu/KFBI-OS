/*=============================================================================
*   
*   Filename : ReacDiffSolver1dN-D.h
*   Creator : Han Zhou
*   Date : 11/22/21
*   Description : 
*
=============================================================================*/
   
 
#ifndef _REACDIFFSOLVER1DN_D_H
#define _REACDIFFSOLVER1DN_D_H


#include "Variables.h"

extern void 
ReacDiffSolver1dN_D(const VectorXd &x, const VectorXb &interior, int n,
										const VectorXd &bdry_crd, int bdry_num, double kap_i,
										double kap_e, const VectorXd &F, 
										const VectorXd &ux_bdry, double ul, double ur,
										VectorXd &u, VectorX6d &u_bdry_data);

#endif

