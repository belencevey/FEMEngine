
/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: TRI3_K_SBS.C                                               #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #        Routines for assembling of global stiffness matrix for             #
 #        Linear Stress triangular elements.                                 #
 #        side by side assembling.                                           #
 #        Element code = 7.                                                  #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  TRI3_K_SBS()  : side by side computation of triangular element     #
 #                        stiffness matrix for linear stress triangle.       #
 #                                                                           #
 #****************************************************************************/


#include <math.h>

/*******************************************[ USER LIBRARY : TRIMEMPK.C ]*****
 #                                                                           #
 #  FUNCTION :  TRI3_K_SBS()                                                 #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of element stiffness matrices in global coordinates      #
 #      and assembling into the global stiffness matrix.                     #
 #      Side by side computations of adjacent elements stiffness matrices.   #
 #                                                                           #
 #      Linear Stress triangular elements.(side by side assembling)          #
 #      Element code = 7                                                     #
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

void TRI3_K_SBS (int sides    , int bsides     , int ndof,
               int nudim  , int *id      , int *side_nodes, int *side_bc,
               int xnsn,
			   double *xyz, int *side_mat, int *mat       , double *prop,
               int numpr , double *sg   , int *maxa      , int *lm,
               double *sv)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// int     sides;   /* number of sides.                                       */
// int    bsides;   /* number of boundary sides.                              */
// int     ndof ;   /* maximum number of nodal degrees of freedom.            */
// int     nudim;   /* number of spatial dimensions.                          */
// int    *id   ;   /* identification vector [numnp * ndof].                  */
// int *side_nodes; /* nodes of each side       [sides * xnsn]                */
// int  *side_bc;   /* side boundary conditions [bsides * xndof]             */
// int   xnsn   ;   /* maximum number of side nodes.                          */
// double *xyz  ;   /* nodal coordinates [numnp * nudim].                     */
// int *side_mat;   /* side material index of adjacent elements [2*sides[0]]  */
// int    *mat  ;   /* element material type vector [numel].                  */
// double *prop ;   /* material properties vector [numat * numpr].            */
// int     numpr;   /* maximum number of material properties.                 */
// double  *sg  ;   /* global stiffness matrix stored by skyline [nwk]        */
// int     *maxa;   /* pointer of diagonal elements for skyline [neq+1]       */
// int     *lm  ;   /* identification vector of element dof [nrw]             */
// double  *sv  ;   /* upper triangular element matrix [nrw *(nrw +1)/2]      */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

 void TRI3_KS_BS (int bsides     , int xndof,
                 int nudim  , int *id      , int *side_nodes, int *side_bc,
                 int xnsn,
				 double *xyz, int *side_mat, int *mat       , double *prop,
                 int xnumpr , double *sg   , int *maxa      , int *lm,
                 double *sv);

 void TRI3_KS_IS (int sides, int bsides     , int xndof,
                 int nudim  , int *id      , int *side_nodes, int *side_bc,
                 int xnsn,
			     double *xyz, int *side_mat, int *mat       , double *prop,
                 int xnumpr , double *sg   , int *maxa      , int *lm,
                 double *sv);
/*__________________________________________________________________________

				LOOP ON BOUNDARY SIDES
___________________________________________________________________________*/

  TRI3_KS_BS (bsides, ndof, nudim, id, side_nodes, side_bc,
             xnsn, xyz, side_mat, mat, prop, numpr, sg, maxa, lm, sv);

/*__________________________________________________________________________

				LOOP ON INTERIOR SIDES
___________________________________________________________________________*/

  TRI3_KS_IS (sides, bsides, ndof, nudim, id, side_nodes, side_bc,
             xnsn, xyz, side_mat, mat, prop, numpr, sg, maxa, lm, sv);


/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of TRI3_K_SBS() */



/*******************************************[ USER LIBRARY : TRIMEMPK.C ]*****
 #                                                                           #
 #  FUNCTION :  TRI3_KS_BS ()                                                 #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of element stiffness matrices in global coordinates      #
 #      and assembling into the global stiffness matrix.                     #
 #      Side by side computations of adjacent elements stiffness matrices.   #
 #                                                                           #
 #      ONLY BOUNDARY SIDES                                                  #
 #                                                                           #
 #      Linear Stress triangular elements.(side by side assembling)          #
 #      Element code = 7                                                     #
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

 void TRI3_KS_BS (int bsides     , int ndof,
                 int nudim  , int *id      , int *side_nodes, int *side_bc,
                 int xnsn,
				 double *xyz, int *side_mat, int *mat       , double *prop,
                 int numpr , double *sg   , int *maxa      , int *lm,
                 double *sv)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// int    bsides;   /* number of boundary sides.                              */
// int     ndof ;   /* maximum number of nodal degrees of freedom.            */
// int     nudim;   /* number of spatial dimensions.                          */
// int    *id   ;   /* identification vector [numnp * ndof].                  */
// int *side_nodes; /* nodes of each side       [sides * xnsn]                */
// int  *side_bc;   /* side boundary conditions [bsides * xndof]              */
// int   xnsn   ;   /* maximum number of side nodes.                          */
// double *xyz  ;   /* nodal coordinates [numnp * nudim].                     */
// int *side_mat;   /* side material index of adjacent elements [2*sides[0]]  */
// int    *mat  ;   /* element material type vector [numel].                  */
// double *prop ;   /* material properties vector [numat * numpr].            */
// int     numpr;   /* maximum number of material properties.                 */
// double  *sg  ;   /* global stiffness matrix stored by skyline [nwk]        */
// int     *maxa;   /* pointer of diagonal elements for skyline [neq+1]       */
// int     *lm  ;   /* identification vector of element dof [nrw]             */
// double  *sv  ;   /* upper triangular element matrix [nrw *(nrw +1)/2]      */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int     nrw  ;   /* number of matrix rows = 4 for boundary sides           */
                  /* number of matrix rows = 4 for interior sides           */
 int     ic   ;   /* side counter for loop on sides.                        */
 int     ni[5];   /* side node numbers.                                     */
 int     i, j ;   /* loop counters.                                         */
 int     i1   ;   /* auxiliar pointer for mesh arrays.                      */

/*____________________________________________________________________________

		  DECLARATION OF MATERIAL PROPERTIES
____________________________________________________________________________*/

 double  E    ;   /* elasticity modulus                                     */
 double  nu   ;   /* Poisson's ratio                                        */
 double  h    ;   /* thickness                                              */

/*____________________________________________________________________________

		        DECLARATION OF SIDE NODAL COORDINATES
____________________________________________________________________________*/

 double  xn1  ;   /* X coordinate of vertex node 1                          */
 double  xn2  ;   /* X coordinate of vertex node 2                          */
 double  xn3  ;   /* X coordinate of opposite node 3                        */
 double  yn1  ;   /* Y coordinate of vertex node 1                          */
 double  yn2  ;   /* Y coordinate of vertex node 2                          */
 double  yn3  ;   /* Y coordinate of opposite node 3                        */

/*____________________________________________________________________________

		           DECLARATION OF GEOMETRIC CONSTANTS
____________________________________________________________________________*/

 double  xn3n1, yn1n3, xn2n3, yn3n2;
 double  darea;   /* two times the element area                             */

/*____________________________________________________________________________

			DECLARATION OF STIFFNESS MATRIX CONSTANTS
____________________________________________________________________________*/

 double  f1;  /* 0.5*h/darea*E/(1-nu*nu)                                        */
 double  f2;  /* 0.5*(1-nu)                                                 */
 double  f3;  /* 0.5*(1+nu)                                                 */

/*____________________________________________________________________________

			        DECLARATION OF SIDE STIFFNESS MATRIX
____________________________________________________________________________*/

 double sk[5][5];


/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

 void addban (double *sg, int *maxa, double *sv, int *lm, int nrw);

/*__________________________________________________________________________

				LOOP ON BOUNDARY SIDES
___________________________________________________________________________*/

 for (ic=1; ic<=bsides ; ic++)
 {
/*________________________________________________________________________

			NODAL INCIDENCES OF BOUNDARY SIDE
 __________________________________________________________________________*/

		i1 = xnsn*(ic-1);
		ni[1] = side_nodes [++i1];
		ni[2] = side_nodes [++i1];
		ni[3] = side_nodes [++i1];

/*________________________________________________________________________

                        NODAL COORDINATES
 __________________________________________________________________________*/

        i1 = nudim*(ni[1] - 1);
	    xn1 = xyz[++i1];
	    yn1 = xyz[++i1];

	    i1 = nudim*(ni[2] - 1);
	    xn2 = xyz[++i1];
	    yn2 = xyz[++i1];

	    i1 = nudim*(ni[3] - 1);
	    xn3 = xyz[++i1];
	    yn3 = xyz[++i1];

/*__________________________________________________________________________

			 ELEMENT GEOMETRIC CONSTANTS
___________________________________________________________________________*/

		xn3n1 = xn3 - xn1;
		yn1n3 = yn1 - yn3;

		xn2n3 = xn2 - xn3;
		yn3n2 = yn3 - yn2;

		darea = xn2n3*yn1n3 - xn3n1*yn3n2;

/*__________________________________________________________________________

				MATERIAL PROPERTIES OF ELEMENT
___________________________________________________________________________*/

		i1 = numpr * (side_mat[2*ic-1] - 1);

		E  = prop[++i1];
		nu = prop[++i1];
		h  = prop[++i1];

/*__________________________________________________________________________

			 ELEMENT STIFFNESS CONSTANTS (PLANE STRESS)
___________________________________________________________________________*/

		f1 =  0.5*(h/darea)*E/(1-nu*nu);
		f2 =  0.5*(1-nu);
		f3 =  0.5*(1+nu);
/*__________________________________________________________________________

							1/2 MATRIX Kn1n1
___________________________________________________________________________*/

   sk [1][1] = 0.5 * f1 * (yn3n2*yn3n2 + f2*xn2n3*xn2n3);
   sk [1][2] = 0.5 * f1 * f3*xn2n3*yn3n2;
   sk [2][2] = 0.5 * f1 * (xn2n3*xn2n3 + f2*yn3n2*yn3n2);

/*__________________________________________________________________________

							 MATRIX Kn1n2
___________________________________________________________________________*/

   sk [1][3] = f1 * (yn3n2*yn1n3 + f2*xn2n3*xn3n1);
   sk [1][4] = f1 * (nu*yn3n2*xn3n1 + f2*xn2n3*yn1n3);
   sk [2][3] = f1 * (nu*xn2n3*yn1n3 + f2*xn3n1*yn3n2);
   sk [2][4] = f1 * (xn2n3*xn3n1 + f2*yn3n2*yn1n3);

/*__________________________________________________________________________

							1/2 MATRIX Kn2n2
___________________________________________________________________________*/

   sk [3][3] = 0.5 * f1 * (yn1n3*yn1n3 + f2*xn3n1*xn3n1);
   sk [3][4] = 0.5 * f1 * f3*xn3n1*yn1n3;
   sk [4][4] = 0.5 * f1 * (xn3n1*xn3n1 + f2*yn1n3*yn1n3);

/*___________________________________________________________________________

		EFFECTIVE STIFFNESS MATRIX GLOBAL COORDINATES
____________________________________________________________________________*/

   nrw  = 4;

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

/*__________________________________________________________________________

		 ADDS TO GLOBAL STIFFNESS MATRIX
___________________________________________________________________________*/

	 addban (sg, maxa, sv, lm, nrw);

/*________________________________________________________________________

		  END OF LOOP ON BOUNDARY SIDES
 ________________________________________________________________________*/

	 }
/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of TRI3_KS_BS() */






/*******************************************[ USER LIBRARY : TRIMEMPK.C ]*****
 #                                                                           #
 #  FUNCTION :  TRI3_KS_IS ()                                                 #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of element stiffness matrices in global coordinates      #
 #      and assembling into the global stiffness matrix.                     #
 #      Side by side computations of adjacent elements stiffness matrices.   #
 #                                                                           #
 #      ONLY INTERIOR SIDES                                                  #
 #                                                                           #
 #      Linear Stress triangular elements.(side by side assembling)          #
 #      Element code = 7                                                     #
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

 void TRI3_KS_IS (int sides, int bsides     , int ndof,
                 int nudim  , int *id      , int *side_nodes, int *side_bc,
                 int xnsn,
			     double *xyz, int *side_mat, int *mat       , double *prop,
                 int numpr , double *sg   , int *maxa      , int *lm,
                 double *sv)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// int     sides;   /* number of sides.                                       */
// int    bsides;   /* number of boundary sides.                              */
// int     ndof ;   /* maximum number of nodal degrees of freedom.            */
// int     nudim;   /* number of spatial dimensions.                          */
// int    *id   ;   /* identification vector [numnp * ndof].                  */
// int *side_nodes; /* nodes of each side       [sides * xnsn]                */
// int  *side_bc;   /* side boundary conditions [bsides * xndof]             */
// int   xnsn   ;   /* maximum number of side nodes.                          */
// double *xyz  ;   /* nodal coordinates [numnp * nudim].                     */
// int *side_mat;   /* side material index of adjacent elements [2*sides[0]]  */
// int    *mat  ;   /* element material type vector [numel].                  */
// double *prop ;   /* material properties vector [numat * numpr].            */
// int     numpr;   /* maximum number of material properties.                 */
// double  *sg  ;   /* global stiffness matrix stored by skyline [nwk]        */
// int     *maxa;   /* pointer of diagonal elements for skyline [neq+1]       */
// int     *lm  ;   /* identification vector of element dof [nrw]             */
// double  *sv  ;   /* upper triangular element matrix [nrw *(nrw +1)/2]      */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int     nrw  ;   /* number of matrix rows = 4 for boundary sides           */
                  /* number of matrix rows = 4 for interior sides           */
 int     ic   ;   /* side counter for loop on sides.                        */
 int     ni[5];   /* side node numbers.                                     */
 int     i, j ;   /* loop counters.                                         */
 int     i1   ;   /* auxiliar pointer for mesh arrays.                      */

/*____________________________________________________________________________

		  DECLARATION OF MATERIAL PROPERTIES
____________________________________________________________________________*/

 double  E1   ;   /* elasticity modulus element 1                           */
 double  nu1  ;   /* Poisson's ratio element 1                              */
 double  h1   ;   /* thickness element 1                                    */
 double  E2   ;   /* elasticity modulus element 2                           */
 double  nu2  ;   /* Poisson's ratio element 2                              */
 double  h2   ;   /* thickness element 2                                    */

/*____________________________________________________________________________

		        DECLARATION OF SIDE NODAL COORDINATES
____________________________________________________________________________*/

 double  xn1  ;   /* X coordinate of vertex node 1                          */
 double  xn2  ;   /* X coordinate of vertex node 2                          */
 double  xn3  ;   /* X coordinate of right opposite node 3                  */
 double  xn4  ;   /* X coordinate of left opposite node 4                   */
 double  yn1  ;   /* Y coordinate of vertex node 1                          */
 double  yn2  ;   /* Y coordinate of vertex node 2                          */
 double  yn3  ;   /* Y coordinate of right opposite node 3                  */
 double  yn4  ;   /* Y coordinate of left opposite node 4                   */

/*____________________________________________________________________________

		           DECLARATION OF GEOMETRIC CONSTANTS
____________________________________________________________________________*/

 double  xn3n1, yn1n3, xn2n3, yn3n2, xn1n4, yn4n1, xn4n2, yn2n4;
 double  darea1;  /* two times the element 1 area                           */
 double  darea2;  /* two times the element 2 area                           */

/*____________________________________________________________________________

			DECLARATION OF STIFFNESS MATRIX CONSTANTS
____________________________________________________________________________*/

 double  f1;  /* 0.5*h/darea*E/(1-nu*nu)                                    */
 double  f2;  /* 0.5*(1-nu)                                                 */
 double  f3;  /* 0.5*(1+nu)                                                 */

/*____________________________________________________________________________

			        DECLARATION OF SIDE STIFFNESS MATRIX
____________________________________________________________________________*/

 double sk[5][5];

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

 void addban (double *sg, int *maxa, double *sv, int *lm, int nrw);

/*__________________________________________________________________________

				LOOP ON INTERIOR SIDES
___________________________________________________________________________*/

 for (ic=bsides+1; ic<=sides ; ic++)
	 {
/*________________________________________________________________________

			NODAL INCIDENCES OF INTERIOR SIDE
 __________________________________________________________________________*/

		i1 = xnsn*(ic-1);
		ni[1] = side_nodes [++i1];
		ni[2] = side_nodes [++i1];
		ni[3] = side_nodes [++i1];
		ni[4] = side_nodes [++i1];

/*________________________________________________________________________

                        NODAL COORDINATES
 __________________________________________________________________________*/

        i1 = nudim*(ni[1] - 1);
	    xn1 = xyz[++i1];
	    yn1 = xyz[++i1];

	    i1 = nudim*(ni[2] - 1);
	    xn2 = xyz[++i1];
	    yn2 = xyz[++i1];

	    i1 = nudim*(ni[3] - 1);
	    xn3 = xyz[++i1];
	    yn3 = xyz[++i1];

	    i1 = nudim*(ni[4] - 1);
	    xn4 = xyz[++i1];
	    yn4 = xyz[++i1];

/*__________________________________________________________________________

			 ELEMENT GEOMETRIC CONSTANTS
___________________________________________________________________________*/

		xn3n1 = xn3 - xn1;
		yn1n3 = yn1 - yn3;

		xn2n3 = xn2 - xn3;
		yn3n2 = yn3 - yn2;

		darea1 = xn2n3*yn1n3 - xn3n1*yn3n2;

		xn4n2 = xn4 - xn2;
		yn2n4 = yn2 - yn4;

		xn1n4 = xn1 - xn4;
		yn4n1 = yn4 - yn1;

		darea2 = xn1n4*yn2n4 - xn4n2*yn4n1;

/*__________________________________________________________________________

				MATERIAL PROPERTIES OF ELEMENT 1
___________________________________________________________________________*/

		i1 = numpr * (side_mat[2*ic-1] - 1);

		E1  = prop[++i1];
		nu1 = prop[++i1];
		h1  = prop[++i1];

/*__________________________________________________________________________

			 ELEMENT 1 STIFFNESS CONSTANTS (PLANE STRESS)
___________________________________________________________________________*/

		f1 =  0.5*(h1/darea1)*E1/(1-nu1*nu1);
		f2 =  0.5*(1-nu1);
		f3 =  0.5*(1+nu1);

/*__________________________________________________________________________

							1/2 MATRIX Kn1n1
___________________________________________________________________________*/

   sk [1][1] = 0.5 * f1 * (yn3n2*yn3n2 + f2*xn2n3*xn2n3);
   sk [1][2] = 0.5 * f1 * f3*xn2n3*yn3n2;
   sk [2][2] = 0.5 * f1 * (xn2n3*xn2n3 + f2*yn3n2*yn3n2);

/*__________________________________________________________________________

							 MATRIX Kn1n2
___________________________________________________________________________*/

   sk [1][3] = f1 * (yn3n2*yn1n3 + f2*xn2n3*xn3n1);
   sk [1][4] = f1 * (nu1*yn3n2*xn3n1 + f2*xn2n3*yn1n3);
   sk [2][3] = f1 * (nu1*xn2n3*yn1n3 + f2*xn3n1*yn3n2);
   sk [2][4] = f1 * (xn2n3*xn3n1 + f2*yn3n2*yn1n3);

/*__________________________________________________________________________

							1/2 MATRIX Kn2n2
___________________________________________________________________________*/

   sk [3][3] = 0.5 * f1 * (yn1n3*yn1n3 + f2*xn3n1*xn3n1);
   sk [3][4] = 0.5 * f1 * f3*xn3n1*yn1n3;
   sk [4][4] = 0.5 * f1 * (xn3n1*xn3n1 + f2*yn1n3*yn1n3);

/*__________________________________________________________________________

				MATERIAL PROPERTIES OF ELEMENT 2
___________________________________________________________________________*/

		i1 = numpr * (side_mat[2*ic] - 1);

		E2  = prop[++i1];
		nu2 = prop[++i1];
		h2  = prop[++i1];

/*__________________________________________________________________________

			 ELEMENT 1 STIFFNESS CONSTANTS (PLANE STRESS)
___________________________________________________________________________*/

		f1 =  0.5*(h2/darea2)*E2/(1-nu2*nu2);
		f2 =  0.5*(1-nu2);
		f3 =  0.5*(1+nu2);

/*__________________________________________________________________________

							1/2 MATRIX Kn1n1
___________________________________________________________________________*/

   sk [1][1] += 0.5 * f1 * (yn2n4*yn2n4 + f2*xn4n2*xn4n2);
   sk [1][2] += 0.5 * f1 * f3*xn4n2*yn2n4;
   sk [2][2] += 0.5 * f1 * (xn4n2*xn4n2 + f2*yn2n4*yn2n4);

/*__________________________________________________________________________

							 MATRIX Kn1n2
___________________________________________________________________________*/

   sk [1][3] += f1 * (yn2n4*yn4n1 + f2*xn4n2*xn1n4);
   sk [1][4] += f1 * (nu2*yn2n4*xn1n4 + f2*xn4n2*yn4n1);
   sk [2][3] += f1 * (nu2*xn4n2*yn4n1 + f2*xn1n4*yn2n4);
   sk [2][4] += f1 * (xn1n4*xn4n2 + f2*yn2n4*yn4n1);

/*__________________________________________________________________________

							1/2 MATRIX Kn2n2
___________________________________________________________________________*/

   sk [3][3] += 0.5 * f1 * (yn4n1*yn4n1 + f2*xn1n4*xn1n4);
   sk [3][4] += 0.5 * f1 * f3*xn1n4*yn4n1;
   sk [4][4] += 0.5 * f1 * (xn1n4*xn1n4 + f2*yn4n1*yn4n1);

/*___________________________________________________________________________

		EFFECTIVE STIFFNESS MATRIX GLOBAL COORDINATES
____________________________________________________________________________*/

  nrw  = 4;

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

/*__________________________________________________________________________

		 ADDS TO GLOBAL STIFFNESS MATRIX
___________________________________________________________________________*/

		addban ( sg, maxa, sv, lm, nrw);

/*________________________________________________________________________

		  END OF LOOP ON INTERIOR SIDES
 ________________________________________________________________________*/

	 }
/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of TRI3_KS_IS() */

