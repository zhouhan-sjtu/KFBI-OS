/*=============================================================================
*   
*   Filename : ReacDiffSolver1d3-2.h
*   Creator : Han Zhou
*   Date : 11/06/21
*   Description : 
*
=============================================================================*/
 
#ifndef _REACDIFFSOLVER1D3_2_H
#define _REACDIFFSOLVER1D3_2_H

#include "Variables.h"

extern void 
ReacDiffSolver1d3_D(const VectorXd &x, const VectorXb &interior, int n,
										const VectorXd &bdry_crd, int bdry_num, double kap_i,
										double kap_e, const VectorXd &F, 
										const VectorXd &u_bdry, double ul, double ur,
										VectorXd &u, VectorX6d &u_bdry_data);

#endif

