/*=============================================================================
*   
*   Filename : DiffusionSolver1d-001-O4.h
*   Creator : Han Zhou
*   Date : 02/26/22
*   Description : 
*
=============================================================================*/
 
#ifndef _DIFFUSIONSOLVER1D_001_O4_H
#define _DIFFUSIONSOLVER1D_001_O4_H

#include "Variables.h"

extern void DiffusionSolver1d_001_O4(const VectorXd &x, 
																		 const VectorXb &interior, int n,
										 								 const VectorXd &bdry_crd, int bdry_num, 
																		 double kap,
										 								 const VectorXd &F, 
																		 const VectorXd &u_bdry, 
										 								 double ul, double ur,
										 								 VectorXd &u, 
																		 VectorX6d &u_bdry_data);
																		 //VectorX10d &u_bdry_data);

#endif

