/**********************[          TRIMBC2DK.cpp            ]********************/
#include"TRIMBC2DK.h"


/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: TRIMBC2DK.C                                                #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #        Routines for assembling of global stiffness matrix for             #
 #        Linear triangular elements cubic Bezier shape.                     #
 #        Element code = 8.                                                  #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  trimbc2dk()  : computation of triangular element stiffness matrix. #
 #                                                                           #
 #****************************************************************************/

/*******************************************[ USER LIBRARY : TRIM2DK.C  ]*****
 #                                                                           #
 #  FUNCTION :  trimbc2dk ()                                                 #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of element stiffness matrix in global coordinates        #
 #      and assembling into the global stiffness matrix.                     #
 #                                                                           #
 #      Linear triangular elements plane stress.                             #
 #      Element code = 8.                                                    #
 #                                                                           #
 #      Only the upper triangular part of the element matrix is computed     #
 #      and stored in vector form by rows. Then is assembled into the        #
 #      global stiffness matrix.                                             #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED LIBRARY FUNCTIONS USED ]--------------------------------#
 #                                                                           #
 #                                                                           #
 #  SKYLINE.C :                                                              #
 #    void  addban()  : assembles element matrix in compacted global matrix. #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  MATH.H  :                                                                #
 #    double sqrt() : calculates square root.                                #
 #                                                                           #
 #                                                                           #
 #***************************************************************************/


  void trimbc2dk(int feblk  , int numel     , int ndof      , int nudim,
	             int *id    , int *elm_sides, int *side_nodes, int xnsn,
	             double *xyz, int *mat      , double *prop   , int numpr,
	             double *sg , int *maxa     , int *lm        , double *sv,
	             double *Gij, double *Kij   , double *diag)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

/* int     feblk;    first element number in this block.                    */
/* int     numel;    number of elements in this block.                      */
/* int     ndof ;    maximum number of nodal degrees of freedom.            */
/* int     nudim;    number of spatial dimensions.                          */
/* int    *id   ;    identification vector [numnp * ndof].                  */
/* int *side_nodes;  nodes of each side       [sides * xnsn]                */
/* int   xnsn   ;    maximum number of side nodes.                          */
/* int *elm_sides ;  sides of each element    [numel * xnes]                */
/* double *xyz  ;    nodal coordinates [numnp * nudim].                     */
/* int    *mat  ;    element material type vector [numel].                  */
/* double *prop ;    material properties vector [numat * numpr].            */
/* int     numpr;    maximum number of material properties.                 */
/* double  *sg  ;    global stiffness matrix stored by skyline [nwk]        */
/* int     *maxa;    pointer of diagonal elements for skyline [neq+1]       */
/* int     *lm  ;    identification vector of element dof [nrw]             */
/* double  *sv  ;    upper triangular element matrix [nrw *(nrw +1)/2]      */
/* double *Gij  ;    element stiffness coefficients [10*numel]              */
/* double *Kij  ;    element stiffness coefficients [21*numel]              */
/* double *diag ;    stiffness diagonal coefficients                        */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int     nrw  ;   /* number of matrix rows = 6.                             */
 int     ie   ;   /* element counter for loop on elements.                  */
 int     L[4] ;   /* element side numbers.                                  */
 int     ni[11];   /* element node numbers.                                  */
 int     i, j ;   /* loop counters.                                         */
 int     i1   ;   /* auxiliar pointer for mesh arrays.                      */

/*____________________________________________________________________________

			  DECLARATION OF NODAL COORDINATES
____________________________________________________________________________*/

 double  x1   ;   /* X coordinate of node 1.                                */
 double  x2   ;   /* X coordinate of node 2.                                */
 double  x3   ;   /* X coordinate of node 3.                                */
 double  y1   ;   /* Y coordinate of node 1.                                */
 double  y2   ;   /* Y coordinate of node 2.                                */
 double  y3   ;   /* Y coordinate of node 3.                                */

 double x[11];   // vector of X coordinates
 double y[11];   // vector of Y coordinates

/*____________________________________________________________________________

		  DECLARATION OF MATERIAL PROPERTIES
____________________________________________________________________________*/

 double  E    ;   /* elasticity modulus                                     */
 double  nu   ;   /* Poisson's ratio                                        */
 double  h    ;   /* thickness                                              */

/*____________________________________________________________________________

		DECLARATION OF COMPUTED GEOMETRIC PROPERTIES
____________________________________________________________________________*/

 double  x13, y31, x32, y23;
 double  d    ;   /* two times the element area = x13*y23-x32*y31           */
 double  G11, G12, G13, G14, G22, G23, G24, G33, G34, G44;

 /*____________________________________________________________________________

         DECLARATION OF SHAPE FUNCTION DERIVATIVES
 ____________________________________________________________________________*/

 double  N1x, N1y, N2x, N2y;
/*____________________________________________________________________________

			DECLARATION OF SQUARE ELEMENT MATRIX
____________________________________________________________________________*/

 double sk[7][7];

 /*____________________________________________________________________________

            DECLARATION OF BEZIER ARRAYS
 ____________________________________________________________________________*/

 const double  BE[11] = {0., 1./27, 1./27, 1./27, 1./9, 1./9, 1./9, 1./9, 1./9, 1./9, 2./9};
 const double DBE[11][4] = { {0.,0.,0.,0.},\
                             {0.,1./3,0.,0.}, {0.,0.,1./3.,0.}, {0.,0.,0.,1./3}, \
                             {0.,2./3,1./3,0.}, {0.,1./3,2./3,0.}, {0.,0.,2./3,1./3}, \
                             {0.,0.,1./3,2./3}, {0.,1./3,0.,2./3}, {0.,2./3,0.,1./3}, \
                             {0.,2./3,2./3,2./3} };
 double DBW[4];
 double  DM[11][4];
 double  JX[4], JY[4];
 double  w[11], BW, JX21, JX13, JY31, JY12, JX32, JY23, J;

 /*____________________________________________________________________________

			DECLARATION OF MATRIX FACTORS AND COEFFICIENTS
____________________________________________________________________________*/

 double  f1;  /* 0.5*h/d*E/(1-nu*nu)                                        */
 double  f2;  /* 0.5*(1-nu)                                                 */

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

 void addban (double *sg, int *maxa, double *sv, int *lm, int nrw);

/*__________________________________________________________________________

			 SET ELEMENT PARAMETERS
___________________________________________________________________________*/

  nrw  = 6;
/*__________________________________________________________________________

				LOOP ON ELEMENTS
___________________________________________________________________________*/

 for (ie=feblk; ie<=numel ; ie++)
	 {
/*________________________________________________________________________

			NODAL INCIDENCES
 __________________________________________________________________________*/


		i1 = 3*(ie-1);
		L[1] = elm_sides [++i1];
		L[2] = elm_sides [++i1];
		L[3] = elm_sides [++i1];

		ni[1] = (L[1]>0) ? side_nodes[xnsn*(L[1]-1)+1] : side_nodes[xnsn*(-L[1]-1)+2];
		ni[2] = (L[2]>0) ? side_nodes[xnsn*(L[2]-1)+1] : side_nodes[xnsn*(-L[2]-1)+2];
		ni[3] = (L[3]>0) ? side_nodes[xnsn*(L[3]-1)+1] : side_nodes[xnsn*(-L[3]-1)+2];

		ni[4] = (L[1]>0) ? side_nodes[xnsn*(L[1] - 1) + 5] : side_nodes[xnsn*(-L[1] - 1) + 6];
		ni[5] = (L[1]>0) ? side_nodes[xnsn*(L[1] - 1) + 6] : side_nodes[xnsn*(-L[1] - 1) + 5];
		ni[6] = (L[2]>0) ? side_nodes[xnsn*(L[2] - 1) + 5] : side_nodes[xnsn*(-L[2] - 1) + 6];
		ni[7] = (L[2]>0) ? side_nodes[xnsn*(L[2] - 1) + 6] : side_nodes[xnsn*(-L[2] - 1) + 5];
		ni[8] = (L[3]>0) ? side_nodes[xnsn*(L[3] - 1) + 5] : side_nodes[xnsn*(-L[3] - 1) + 6];
		ni[9] = (L[3]>0) ? side_nodes[xnsn*(L[3] - 1) + 6] : side_nodes[xnsn*(-L[3] - 1) + 5];

/*____________________________________________________________________________

				 NODAL COORDINATES
____________________________________________________________________________*/

		i1 = nudim*(ni[1] - 1);
		x1 = xyz[++i1];
		y1 = xyz[++i1];

		i1 = nudim*(ni[2] - 1);
		x2 = xyz[++i1];
		y2 = xyz[++i1];

		i1 = nudim*(ni[3] - 1);
		x3 = xyz[++i1];
		y3 = xyz[++i1];

		for (i = 1; i <= 9; i++) {
			i1 = nudim * (ni[i] - 1);
			x[i] = xyz[++i1];
			y[i] = xyz[++i1];
			w[i] = xyz[++i1];
		}

		x[10] = (x[1] + x[2] + x[3]) / 3.0;
		y[10] = (y[1] + y[2] + y[3]) / 3.0;
		w[10] = 1;

/*__________________________________________________________________________

				MATERIAL PROPERTIES
___________________________________________________________________________*/


		i1 = numpr * (mat[ie] - 1);

		E  = prop[++i1];
		nu = prop[++i1];
		h  = prop[++i1];

/*__________________________________________________________________________

			 ELEMENT CONSTANTS
___________________________________________________________________________*/

		x13 = x1 - x3;
		y31 = y3 - y1;

		x32 = x3 - x2;
		y23 = y2 - y3;

		d = x13*y23 - x32*y31;

		f1 =  0.5*(h/d)*E/(1-nu*nu);
		f2 =  0.5*(1-nu);

/*__________________________________________________________________________

	 	           SHAPE FUNCTIONS DERIVATIVES
___________________________________________________________________________*/

		N1x = y23;
		N1y = x32;
		N2x = y31;
		N2y = x13;
/*__________________________________________________________________________

	 	       SHAPE FUNCTIONS DERIVATIVES FOR CUBIC BEZIER
___________________________________________________________________________*/

		BW = 0.0;

		for (i = 1; i <= 10; i++) BW += BE[i] * w[i];

		for (j = 1; j <= 3; j++) {
			DBW[j] = 0.0;
			for (i = 1; i <= 10; i++)  DBW[j] += DBE[i][j]*w[i];
			for (i = 1; i <= 10; i++)  DM[i][j] = (DBE[i][j]/BW - BE[i]*DBW[j]/BW/BW)*w[i];
			JX[j] = 0.0;
			JY[j] = 0.0;
			for (i = 1; i <= 10; i++) {
				JX[j] += DM[i][j] * x[i];
				JY[j] += DM[i][j] * y[i];
			}
		}

		JX21 = JX[2] - JX[1];
		JX13 = JX[1] - JX[3];
		JY31 = JY[3] - JY[1];
		JY12 = JY[1] - JY[2];

		d = (JX21*JY31-JY12*JX13);

		JX32 = JX[3] - JX[2];
		JY23 = JY[2] - JY[3];

		N1x = JY23;
		N1y = JX32;
		N2x = JY31;
		N2y = JX13;

		f1 = 0.5*(h / d)*E / (1 - nu * nu);
		f2 = 0.5*(1 - nu);



/*__________________________________________________________________________

          MATRIX G (ONLY LINEAR TRIANGLE) (HUGHES/BELYTSCHKO STYLE)
___________________________________________________________________________*/

//		G11 = f1 * (   y23*y23 + f2*x32*x32);   // K[1][1] =    N1x*N1x + f2*N1y*N1y
//		G12 = f1 * (   y23*y31 + f2*x32*x13);   // K[1][3] =    N1x*N2x + f2*N1y*N2y
//		G13 = f1 * (nu*y23*x32 + f2*x32*y23);   // K[1][2] = nu*N1x*N1y + f2*N1y*N1x
//		G14 = f1 * (nu*y23*x13 + f2*x32*y31);   // K[1][4] = nu*N1x*N2y + f2*N1y*N2x
//		G22 = f1 * (   y31*y31 + f2*x13*x13);   // K[3][3] =    N2x*N2x + f2*N2y*N2y
//		G23 = f1 * (nu*y31*x32 + f2*x13*y23);   // K[2][3] = nu*N2x*N1y + f2*N2y*N1x
//		G24 = f1 * (nu*y31*x13 + f2*x13*y31);   // K[3][4] = nu*N2x*N2y + f2*N2y*N2x
//		G33 = f1 * (   x32*x32 + f2*y23*y23);   // K[2][2] =    N1y*N1y + f2*N1x*N1x
//		G34 = f1 * (   x32*x13 + f2*y23*y31);   // K[2][4] =    N1y*N2y + f2*N1x*N2x
//		G44 = f1 * (   x13*x13 + f2*y31*y31);   // K[4][4] =    N2y*N2y + f2*N2x*N2x

		G11 = f1 * (   N1x*N1x + f2 * N1y*N1y);   // K[1][1] =f1 * (   N1x*N1x + f2*N1y*N1y)
		G12 = f1 * (   N1x*N2x + f2 * N1y*N2y);   // K[1][3] =f1 * (   N1x*N2x + f2*N1y*N2y)
		G13 = f1 * (nu*N1x*N1y + f2 * N1y*N1x);   // K[1][2] =f1 * (nu*N1x*N1y + f2*N1y*N1x)
		G14 = f1 * (nu*N1x*N2y + f2 * N1y*N2x);   // K[1][4] =f1 * (nu*N1x*N2y + f2*N1y*N2x)
		G22 = f1 * (   N2x*N2x + f2 * N2y*N2y);   // K[3][3] =f1 * (   N2x*N2x + f2*N2y*N2y)
		G23 = f1 * (nu*N2x*N1y + f2 * N2y*N1x);   // K[2][3] =f1 * (nu*N2x*N1y + f2*N2y*N1x)
		G24 = f1 * (nu*N2x*N2y + f2 * N2y*N2x);   // K[3][4] =f1 * (nu*N2x*N2y + f2*N2y*N2x)
		G33 = f1 * (   N1y*N1y + f2 * N1x*N1x);   // K[2][2] =f1 * (   N1y*N1y + f2*N1x*N1x)
		G34 = f1 * (   N1y*N2y + f2 * N1x*N2x);   // K[2][4] =f1 * (   N1y*N2y + f2*N1x*N2x)
		G44 = f1 * (   N2y*N2y + f2 * N2x*N2x);   // K[4][4] =f1 * (   N2y*N2y + f2*N2x*N2x)
/*____________________________________________________________________________

		 ELEMENT STIFFNESS MATRIX IN GLOBAL COORDINATES  (BATHE STYLE)
____________________________________________________________________________*/


		 sk [1][1] =  G11;
		 sk [1][2] =  G13;
		 sk [1][3] =  G12;
		 sk [1][4] =  G14;
		 sk [1][5] = -(G11+G12);
		 sk [1][6] = -(G13+G14);
		 sk [2][2] =  G33;
		 sk [2][3] =  G23;
		 sk [2][4] =  G34;
		 sk [2][5] = -(G13+G23);
		 sk [2][6] = -(G33+G34);
		 sk [3][3] =  G22;
		 sk [3][4] =  G24;
		 sk [3][5] = -(G12+G22);
		 sk [3][6] = -(G23+G24);
		 sk [4][4] =  G44;
		 sk [4][5] = -(G14+G24);
		 sk [4][6] = -(G34+G44);
		 sk [5][5] =  G11+2*G12+G22;
		 sk [5][6] =  G13+G14+G23+G24;
		 sk [6][6] =  G33+2*G34+G44;


/*___________________________________________________________________________

						 STORE STIFFNESS COEFFICIENTS Gij
____________________________________________________________________________*/

		i1 = 10*(ie-1);
		Gij[++i1] = G11;
		Gij[++i1] = G12;
		Gij[++i1] = G13;
		Gij[++i1] = G14;
		Gij[++i1] = G22;
		Gij[++i1] = G23;
		Gij[++i1] = G24;
		Gij[++i1] = G33;
		Gij[++i1] = G34;
		Gij[++i1] = G44;

/*___________________________________________________________________________

						 STORE STIFFNESS COEFFICIENTS Kij
____________________________________________________________________________*/

		i1 = 21*(ie-1);
		Kij[++i1] = sk[1][1];
      Kij[++i1] = sk[1][2];
      Kij[++i1] = sk[1][3];
      Kij[++i1] = sk[1][4];
      Kij[++i1] = sk[1][5];
      Kij[++i1] = sk[1][6];
      Kij[++i1] = sk[2][2];
      Kij[++i1] = sk[2][3];
      Kij[++i1] = sk[2][4];
      Kij[++i1] = sk[2][5];
      Kij[++i1] = sk[2][6];
      Kij[++i1] = sk[3][3];
      Kij[++i1] = sk[3][4];
      Kij[++i1] = sk[3][5];
      Kij[++i1] = sk[3][6];
      Kij[++i1] = sk[4][4];
      Kij[++i1] = sk[4][5];
      Kij[++i1] = sk[4][6];
      Kij[++i1] = sk[5][5];
      Kij[++i1] = sk[5][6];
      Kij[++i1] = sk[6][6];

/*___________________________________________________________________________

						 ACCUMULATES DIAGONAL STIFFNESS COEFFICIENTS
____________________________________________________________________________*/

		i1 = 2*(ni[1]-1);
		diag[++i1] += sk[1][1];
      diag[++i1] += sk[2][2];

		i1 = 2*(ni[2]-1);
		diag[++i1] += sk[3][3];
      diag[++i1] += sk[4][4];

		i1 = 2*(ni[3]-1);
		diag[++i1] += sk[5][5];
      diag[++i1] += sk[6][6];


/*___________________________________________________________________________

		EFFECTIVE STIFFNESS MATRIX GLOBAL COORDINATES
____________________________________________________________________________*/


		i1 = 0;
		for ( i=1 ; i<=nrw ; i++)
			for( j=i ; j<=nrw ; j++)  sv[++i1] = sk[i][j];

/*__________________________________________________________________________

		EQUATION NUMBERS OF ELEMENT DOF
___________________________________________________________________________*/

			i1 = ndof*(ni[1] - 1);
			lm[1] = id[++i1];
			lm[2] = id[++i1];

			i1 = ndof*(ni[2] - 1);
			lm[3] = id[++i1];
			lm[4] = id[++i1];

			i1 = ndof*(ni[3] - 1);
			lm[5] = id[++i1];
			lm[6] = id[++i1];


/*__________________________________________________________________________

		 ADDS TO GLOBAL STIFFNESS MATRIX
___________________________________________________________________________*/

		addban (sg, maxa, sv, lm, nrw);

/*________________________________________________________________________

		  END OF LOOP ON ELEMENTS
 ________________________________________________________________________*/

	 }
/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trimbc2dk() */

