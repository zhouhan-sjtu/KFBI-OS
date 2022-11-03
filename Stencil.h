/*=============================================================================
*   
*   Filename : Stencil.h
*   Creator : Han Zhou
*   Date : 02/10/22
*   Description : 
*
=============================================================================*/
 
#ifndef _STENCIL_H
#define _STENCIL_H

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// 5-point central difference scheme nearby stencil points

const int FD5P_2D_I[4][2] = {{1,0},{-1,0},{0,1},{0,-1}};

// 9-point compact scheme nearby stencil points

const int CDS_2D_I[8][2] = {{1,0},{-1,0},{0,1},{0,-1},
														{1,1},{1,-1},{-1,1},{-1,-1}};

const int CDS_2D_I1[4][2] = {{1,0},{-1,0},{0,1},{0,-1}};

const int CDS_2D_I2[4][2] = {{1,1},{1,-1},{-1,1},{-1,-1}};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// 7-point central difference scheme nearby stencil points

const int FD7P_3D_I[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};

// 27-point compact scheme nearby stencil points

const int CDS_3D_I[26][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
								  					 {1,1,0},{-1,1,0},{-1,-1,0},{1,-1,0},
								  					 {1,0,1},{-1,0,1},{-1,0,-1},{1,0,-1},
                             {0,1,1},{0,-1,1},{0,-1,-1},{0,1,-1},
								  					 {1,1,1},{-1,1,1},{-1,-1,1},{1,-1,1},
								  					 {1,1,-1},{-1,1,-1},{-1,-1,-1},{1,-1,-1}};

const int CDS_3D_I1[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};

const int CDS_3D_I2[12][3] = {{1,1,0},{-1,1,0},{-1,-1,0},{1,-1,0},
								   	 					{1,0,1},{-1,0,1},{-1,0,-1},{1,0,-1},
								   	 					{0,1,1},{0,-1,1},{0,-1,-1},{0,1,-1}};

const int CDS_3D_I3[8][3] = {{1,1,1},{-1,1,1},{-1,-1,1},{1,-1,1},
														 {1,1,-1},{-1,1,-1},{-1,-1,-1},{1,-1,-1}};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 1D interpolation stencil
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void get1DStencil_2(int d, int stc[2])
{
	stc[0] = 0;
	stc[1] = d;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void get1DStencil_3(int d, int stc[3])
{
	stc[0] = 0;
	stc[1] = d;
	stc[2] = -d;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void get1DStencil_4(int d, int stc[4])
{
	stc[0] = 0;
	stc[1] = d;
	stc[2] = -d;
	stc[3] = d+d;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void get1DStencil_5(int d, int stc[5])
{
	stc[0] = 0;
	stc[1] = d;
	stc[2] = -d;
	stc[3] = d+d;
	stc[4] = -(d+d);
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 2D interpolation stencil
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void get2DStencil_3(int d1, int d2, int stc[3][2])
{
	stc[0][0] = 0;					stc[0][1] = 0;

	stc[1][0] = d1;					stc[1][1] = 0;
	stc[2][0] = 0;					stc[2][1] = d2;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void get2DStencil_6(int d1, int d2, int stc[6][2])
{
	stc[0][0] = 0;					stc[0][1] = 0;

	stc[1][0] = d1;					stc[1][1] = 0;
	stc[2][0] = 0;					stc[2][1] = d2;

	stc[3][0] = d1;					stc[3][1] = d2;
	stc[4][0] = -d1;				stc[4][1] = 0;
	stc[5][0] = 0;					stc[5][1] = -d2;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void get2DStencil_10(int d1, int d2, int stc[10][2])
{
	stc[0][0] = 0;					stc[0][1] = 0;

	stc[1][0] = d1;					stc[1][1] = 0;
	stc[2][0] = 0;					stc[2][1] = d2;

	stc[3][0] = d1;					stc[3][1] = d2;
	stc[4][0] = -d1;				stc[4][1] = 0;
	stc[5][0] = 0;					stc[5][1] = -d2;

	stc[6][0] = d1;					stc[6][1] = -d2;
	stc[7][0] = d1+d1;			stc[7][1] = 0;
	stc[8][0] = 0;					stc[8][1] = d2+d2;
	stc[9][0] = -d1;				stc[9][1] = d2;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void get2DStencil_15(int d1, int d2, int stc[15][2])
{
	stc[0][0] = 0;					stc[0][1] = 0;

	stc[1][0] = d1;					stc[1][1] = 0;
	stc[2][0] = 0;					stc[2][1] = d2;

	stc[3][0] = d1;					stc[3][1] = d2;
	stc[4][0] = -d1;				stc[4][1] = 0;
	stc[5][0] = 0;					stc[5][1] = -d2;

	stc[6][0] = d1;					stc[6][1] = -d2;
	stc[7][0] = d1+d1;			stc[7][1] = 0;
	stc[8][0] = 0;					stc[8][1] = d2+d2;
	stc[9][0] = -d1;				stc[9][1] = d2;

	stc[10][0] = d1+d1;			stc[10][1] = d2;
	stc[11][0] = d1;				stc[11][1] = d2+d2;
	stc[12][0] = -(d1+d1);	stc[12][1] = 0;
	stc[13][0] = -d1;				stc[13][1] = -d2;
	stc[14][0] = 0;					stc[14][1] = -(d2+d2);
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 3D interpolation stencil
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void get3DStencil_4(int d1, int d2, int d3, int stc[4][3])
{
	stc[0][0] = 0;					stc[0][1] = 0;					stc[0][2] = 0;
                    			                  			     
	stc[1][0] = 0;					stc[1][1] = d2;					stc[1][2] = 0;
	stc[2][0] = 0;					stc[2][1] = 0;					stc[2][2] = d3;
	stc[3][0] = d1;					stc[3][1] = 0;					stc[3][2] = 0;
                    			                  			     
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void get3DStencil_10(int d1, int d2, int d3, int stc[10][3])
{
	stc[0][0] = 0;					stc[0][1] = 0;					stc[0][2] = 0;
                    			                  			     
	stc[1][0] = 0;					stc[1][1] = d2;					stc[1][2] = 0;
	stc[2][0] = 0;					stc[2][1] = 0;					stc[2][2] = d3;
	stc[3][0] = d1;					stc[3][1] = 0;					stc[3][2] = 0;
                    			                  			     
	stc[4][0] = 0;					stc[4][1] = d2;					stc[4][2] = d3;
	stc[5][0] = 0;					stc[5][1] = -d2;				stc[5][2] = 0;
	stc[6][0] = 0;					stc[6][1] = 0;					stc[6][2] = -d3;
	stc[7][0] = d1;					stc[7][1] = d2;					stc[7][2] = 0;
	stc[8][0] = d1;					stc[8][1] = 0;					stc[8][2] = d3;
	stc[9][0] = -d1;				stc[9][1] = 0;					stc[9][2] = 0;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void get3DStencil_20(int d1, int d2, int d3, int stc[20][3])
{
	stc[0][0] = 0;					stc[0][1] = 0;					stc[0][2] = 0;
                    			                  			     
	stc[1][0] = 0;					stc[1][1] = d2;					stc[1][2] = 0;
	stc[2][0] = 0;					stc[2][1] = 0;					stc[2][2] = d3;
	stc[3][0] = d1;					stc[3][1] = 0;					stc[3][2] = 0;
                    			                  			     
	stc[4][0] = 0;					stc[4][1] = d2;					stc[4][2] = d3;
	stc[5][0] = 0;					stc[5][1] = -d2;				stc[5][2] = 0;
	stc[6][0] = 0;					stc[6][1] = 0;					stc[6][2] = -d3;
	stc[7][0] = d1;					stc[7][1] = d2;					stc[7][2] = 0;
	stc[8][0] = d1;					stc[8][1] = 0;					stc[8][2] = d3;
	stc[9][0] = -d1;				stc[9][1] = 0;					stc[9][2] = 0;

	stc[10][0] = 0;					stc[10][1] = d2+d2;			stc[10][2] = 0;
	stc[11][0] = 0;					stc[11][1] = d2;				stc[11][2] = -d3;
	stc[12][0] = 0;					stc[12][1] = 0;					stc[12][2] = d3+d3;
	stc[13][0] = 0;					stc[13][1] = -d2;				stc[13][2] = d3;
	stc[14][0] = d1;				stc[14][1] = d2;				stc[14][2] = d3;
	stc[15][0] = d1;				stc[15][1] = -d2;				stc[15][2] = 0;
	stc[16][0] = d1;				stc[16][1] = 0;					stc[16][2] = -d3;
	stc[17][0] = -d1;				stc[17][1] = d2;				stc[17][2] = 0;
	stc[18][0] = -d1;				stc[18][1] = 0;					stc[18][2] = d3;
	stc[19][0] = d1+d1;			stc[19][1] = 0;					stc[19][2] = 0;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void get3DStencil_35(int d1, int d2, int d3, int stc[35][3])
{
	stc[0][0] = 0;					stc[0][1] = 0;					stc[0][2] = 0;
                    			                  			     
	stc[1][0] = 0;					stc[1][1] = d2;					stc[1][2] = 0;
	stc[2][0] = 0;					stc[2][1] = 0;					stc[2][2] = d3;
	stc[3][0] = d1;					stc[3][1] = 0;					stc[3][2] = 0;
                    			                  			     
	stc[4][0] = 0;					stc[4][1] = d2;					stc[4][2] = d3;
	stc[5][0] = 0;					stc[5][1] = -d2;				stc[5][2] = 0;
	stc[6][0] = 0;					stc[6][1] = 0;					stc[6][2] = -d3;
	stc[7][0] = d1;					stc[7][1] = d2;					stc[7][2] = 0;
	stc[8][0] = d1;					stc[8][1] = 0;					stc[8][2] = d3;
	stc[9][0] = -d1;				stc[9][1] = 0;					stc[9][2] = 0;

	stc[10][0] = 0;					stc[10][1] = d2+d2;			stc[10][2] = 0;
	stc[11][0] = 0;					stc[11][1] = d2;				stc[11][2] = -d3;
	stc[12][0] = 0;					stc[12][1] = 0;					stc[12][2] = d3+d3;
	stc[13][0] = 0;					stc[13][1] = -d2;				stc[13][2] = d3;
	stc[14][0] = d1;				stc[14][1] = d2;				stc[14][2] = d3;
	stc[15][0] = d1;				stc[15][1] = -d2;				stc[15][2] = 0;
	stc[16][0] = d1;				stc[16][1] = 0;					stc[16][2] = -d3;
	stc[17][0] = -d1;				stc[17][1] = d2;				stc[17][2] = 0;
	stc[18][0] = -d1;				stc[18][1] = 0;					stc[18][2] = d3;
	stc[19][0] = d1+d1;			stc[19][1] = 0;					stc[19][2] = 0;

	stc[20][0] = 0;					stc[20][1] = d2+d2;			stc[20][2] = d3;
	stc[21][0] = 0;					stc[21][1] = d2;				stc[21][2] = d3+d3;
	stc[22][0] = 0;					stc[22][1] = -(d2+d2);	stc[22][2] = 0;
	stc[23][0] = 0;					stc[23][1] = -d2;				stc[23][2] = -d3;
	stc[24][0] = 0;					stc[24][1] = 0;					stc[24][2] = -(d3+d3);
	stc[25][0] = d1;				stc[25][1] = d2+d2;			stc[25][2] = 0;
	stc[26][0] = d1;				stc[26][1] = d2;				stc[26][2] = -d3;
	stc[27][0] = d1;				stc[27][1] = 0;					stc[27][2] = d3+d3;
	stc[28][0] = d1;				stc[28][1] = -d2;				stc[28][2] = d3;
	stc[29][0] = -d1;				stc[29][1] = d2;				stc[29][2] = d3;
	stc[30][0] = -d1;				stc[30][1] = -d2;				stc[30][2] = 0;
	stc[31][0] = -d1;				stc[31][1] = 0;					stc[31][2] = -d3;
	stc[32][0] = d1+d1;			stc[32][1] = d2;				stc[32][2] = 0;
	stc[33][0] = d1+d1;			stc[33][1] = 0;					stc[33][2] = d3;
	stc[34][0] = -(d1+d1);	stc[34][1] = 0;					stc[34][2] = 0;

}

#endif
