/**********************[          Trimm2dk.cpp            ]********************/
#include"Trimm2dk.h"


/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: TRIMM2DK.C                                                 #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #        Routines for assembling of global stiffness matrix for             #
 #        Linear mixed triangular elements. Element code = 3.                #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  trimm2dk()    : computation of triangular element stiffness matrix.#
 #                                                                           #
 #****************************************************************************/



/*******************************************[ USER LIBRARY : TRIM2DK.C  ]*****
 #                                                                           #
 #  FUNCTION :  trimm2dk ()                                                  #
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
 #      Linear mixed triangular elements plane strain.                       #
 #      Element code = 3.                                                    #
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


 void trimm2dk (int sides, int bsides, int ndof,
                int nudim, int *id, int *side_nodes, int *side_bc,
                int xnsn,
			    double *xyz, int *side_mat, int *mat, double *prop,
                int numpr, double *sg, int *maxa, int *lm,
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

 void trimm2dbs (int bsides     , int xndof,
                 int nudim  , int *id      , int *side_nodes, int *side_bc,
                 int xnsn,
				 double *xyz, int *side_mat, int *mat       , double *prop,
                 int xnumpr , double *sg   , int *maxa      , int *lm,
                 double *sv);

 void trimm2dis (int sides, int bsides     , int xndof,
                 int nudim  , int *id      , int *side_nodes, int *side_bc,
                 int xnsn,
			     double *xyz, int *side_mat, int *mat       , double *prop,
                 int xnumpr , double *sg   , int *maxa      , int *lm,
                 double *sv);

/*__________________________________________________________________________

				LOOP ON BOUNDARY SIDES
___________________________________________________________________________*/

  trimm2dbs (bsides, ndof, nudim, id, side_nodes, side_bc,
             xnsn, xyz, side_mat, mat, prop, numpr, sg, maxa, lm, sv);

/*__________________________________________________________________________

				LOOP ON INTERIOR SIDES
___________________________________________________________________________*/

  trimm2dis (sides, bsides, ndof, nudim, id, side_nodes, side_bc,
             xnsn, xyz, side_mat, mat, prop, numpr, sg, maxa, lm, sv);

/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trimm2dk() */



/*******************************************[ USER LIBRARY : TRIM2DK.C  ]*****
 #                                                                           #
 #  FUNCTION :  trimm2dbs ()                                                 #
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
 #      Linear mixed triangular elements plane strain.                       #
 #      Element code = 3.                                                    #
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




 void trimm2dbs (int bsides, int ndof, int nudim, int *id, int *side_nodes,
                 int *side_bc, int xnsn, double *xyz, int *side_mat, int *mat,
				 double *prop, int numpr, double *sg, int *maxa, int *lm, double *sv)


/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

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

 int     nrw  ;   /* number of matrix rows = 6 for boundary sides           */
                  /* number of matrix rows = 8 for interior sides           */
 int     ic   ;   /* side counter for loop on sides.                        */
 int     ni[4];   /* side node numbers.                                     */
 int     i, j ;   /* loop counters.                                         */
 int     i1   ;   /* auxiliar pointer for mesh arrays.                      */

/*____________________________________________________________________________

		  DECLARATION OF MATERIAL PROPERTIES
____________________________________________________________________________*/

 double  E    ;   /* elasticity modulus                                     */
 double  nu   ;   /* Poisson's ratio                                        */
 double  h    ;   /* thickness                                              */

/*____________________________________________________________________________

			DECLARATION OF SIDE MATRICES
____________________________________________________________________________*/

 double sk[9][9];
 double Cl1e[4][3];
 double Cl2e[4][3];
 double Cl3e[4][3];
 double Dlle[4];
 double Cl1[4][3];
 double Cl2[4][3];
 double Cl3[4][3];
 double Dll[4];
 double Kij[3][3];

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

 void addban (double *sg, int *maxa, double *sv, int *lm, int nrw);
 void trimm2dCl (double Cl1[4][3], double Cl2[4][3], double Cl3[4][3], double Dll[4],
                 int n1, int n2, int n3, double E, double nu, double h, int nudim,
                 double *xyz);
 void trimmKij (double Kij[3][3], double Cli[4][3], double *Dll, double Clj[4][3]);

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

/*__________________________________________________________________________

				MATERIAL PROPERTIES OF ELEMENT
___________________________________________________________________________*/

		i1 = numpr * (side_mat[2*ic-1] - 1);

		E  = prop[++i1];
		nu = prop[++i1];
		h  = prop[++i1];

/*__________________________________________________________________________

									 MATRICES Cl1e, Cl2e, Cl3e, Dlle FOR ELEMENT 1
___________________________________________________________________________*/


  trimm2dCl (Cl1e, Cl2e, Cl3e, Dlle, ni[3], ni[2], ni[1], E, nu, h, nudim, xyz);


/*__________________________________________________________________________

									 ADD TO SIDE MATRICES Cl1, Cl2, Cl3, Dll FOR ELEMENT 1
___________________________________________________________________________*/

		for ( i=1 ; i<=3 ; i++)
			{
    Dll[i]=Dlle[i];
    for ( j=1 ; j<=2 ; j++)
       {
        Cl1[i][j]= Cl3e[i][j];
        Cl2[i][j]= Cl2e[i][j];
        Cl3[i][j]= Cl1e[i][j];
       }
    }
/*__________________________________________________________________________

									 INVERT MATRIX Dll (or zero if tractions are prescribed)
___________________________________________________________________________*/


		Dll [1] = 1/Dll[1];
		Dll [2] = (side_bc[2*ic  ]==1)? 1/Dll[2] : 0;
		Dll [3] = (side_bc[2*ic-1]==1)? 1/Dll[3] : 0;

/*__________________________________________________________________________

									 MATRIX K11
___________________________________________________________________________*/

  trimmKij (Kij, Cl1, Dll, Cl1);

		 sk [1][1] =  Kij[1][1];
		 sk [1][2] =  Kij[1][2];
		 sk [2][2] =  Kij[2][2];

/*__________________________________________________________________________

									 MATRIX K12
___________________________________________________________________________*/

  trimmKij (Kij, Cl1, Dll, Cl2);

		 sk [1][3] =  Kij[1][1];
		 sk [1][4] =  Kij[1][2];
   sk [2][3] =  Kij[2][1];
   sk [2][4] =  Kij[2][2];

/*__________________________________________________________________________

									 MATRIX K13
___________________________________________________________________________*/

  trimmKij (Kij, Cl1, Dll, Cl3);

		 sk [1][5] =  Kij[1][1];
		 sk [1][6] =  Kij[1][2];
   sk [2][5] =  Kij[2][1];
   sk [2][6] =  Kij[2][2];

/*__________________________________________________________________________

									 MATRIX K22
___________________________________________________________________________*/

  trimmKij (Kij, Cl2, Dll, Cl2);

		 sk [3][3] =  Kij[1][1];
		 sk [3][4] =  Kij[1][2];
		 sk [4][4] =  Kij[2][2];

/*__________________________________________________________________________

									 MATRIX K23
___________________________________________________________________________*/

  trimmKij (Kij, Cl2, Dll, Cl3);

		 sk [3][5] =  Kij[1][1];
		 sk [3][6] =  Kij[1][2];
   sk [4][5] =  Kij[2][1];
   sk [4][6] =  Kij[2][2];

/*__________________________________________________________________________

									 MATRIX K33
___________________________________________________________________________*/

  trimmKij (Kij, Cl3, Dll, Cl3);

		 sk [5][5] =  Kij[1][1];
		 sk [5][6] =  Kij[1][2];
		 sk [6][6] =  Kij[2][2];

/*___________________________________________________________________________

		EFFECTIVE STIFFNESS MATRIX GLOBAL COORDINATES
____________________________________________________________________________*/

  nrw  = 6;

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

		addban ( sg, maxa, sv, lm, nrw);

/*________________________________________________________________________

		  END OF LOOP ON BOUNDARY SIDES
 ________________________________________________________________________*/

	 }
/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trimm2dbs() */






/*******************************************[ USER LIBRARY : TRIM2DK.C  ]*****
 #                                                                           #
 #  FUNCTION :  trimm2dis ()                                                 #
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
 #      Linear mixed triangular elements plane strain.                       #
 #      Element code = 3.                                                    #
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


 void trimm2dis (int sides, int bsides, int ndof, int nudim, int *id, int *side_nodes,
                 int *side_bc, int xnsn, double *xyz, int *side_mat, int *mat,
                 double *prop, int numpr, double *sg, int *maxa,  int *lm, double *sv)


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

 int     nrw  ;   /* number of matrix rows = 6 for boundary sides           */
                  /* number of matrix rows = 8 for interior sides           */
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

			DECLARATION OF SIDE MATRICES
____________________________________________________________________________*/

 double sk[9][9];
 double Cl1e[4][3];
 double Cl2e[4][3];
 double Cl3e[4][3];
 double Dlle[4];
 double Cl1[4][3];
 double Cl2[4][3];
 double Cl3[4][3];
 double Cl4[4][3];
 double Dll[4];
 double Kij[3][3];

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

 void addban (double *sg, int *maxa, double *sv, int *lm, int nrw);
 void trimm2dCl (double Cl1[4][3], double Cl2[4][3], double Cl3[4][3], double Dll[4],
                 int n1, int n2, int n3, double E, double nu, double h, int nudim,
                 double *xyz);
 void trimmKij (double Kij[3][3], double Cli[4][3], double *Dll, double Clj[4][3]);

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

/*__________________________________________________________________________

				MATERIAL PROPERTIES OF ELEMENT 1
___________________________________________________________________________*/

		i1 = numpr * (side_mat[2*ic-1] - 1);

		E1  = prop[++i1];
		nu1 = prop[++i1];
		h1  = prop[++i1];

/*__________________________________________________________________________

									 MATRICES Cl1e, Cl2e, Cl3e FOR ELEMENT 1
___________________________________________________________________________*/

  trimm2dCl (Cl1e, Cl2e, Cl3e, Dlle, ni[3], ni[2], ni[1], E1, nu1, h1, nudim, xyz);

/*__________________________________________________________________________

									 ADD TO SIDE MATRICES Cl1, Cl2, Cl3, Dll FOR ELEMENT 1
___________________________________________________________________________*/

		for ( i=1 ; i<=3 ; i++)
			{
    Dll[i]=Dlle[i];
    for ( j=1 ; j<=2 ; j++)
       {
        Cl1[i][j]= Cl3e[i][j];
        Cl2[i][j]= Cl2e[i][j];
        Cl3[i][j]= Cl1e[i][j];
       }
    }

/*__________________________________________________________________________

				MATERIAL PROPERTIES OF ELEMENT 2
___________________________________________________________________________*/

		i1 = numpr * (side_mat[2*ic] - 1);

		E2  = prop[++i1];
		nu2 = prop[++i1];
		h2  = prop[++i1];

/*__________________________________________________________________________

									 MATRICES Cl1e, Cl2e, Cl3e FOR ELEMENT 2
___________________________________________________________________________*/

  trimm2dCl (Cl1e, Cl2e, Cl3e, Dlle, ni[4], ni[1], ni[2], E2, nu2, h2, nudim, xyz);

/*__________________________________________________________________________

									 ADD TO SIDE MATRICES Cl1, Cl2, Cl4 FOR ELEMENT 2
___________________________________________________________________________*/

		for ( i=1 ; i<=3 ; i++)
			{
     Dll[i]= Dll[i] + Dlle[i];
			  for ( j=1 ; j<=2 ; j++)
        {
         Cl1[i][j]= Cl1[i][j] + Cl2e[i][j];
         Cl2[i][j]= Cl2[i][j] + Cl3e[i][j];
         Cl4[i][j]= Cl1e[i][j];
        }
    }

/*__________________________________________________________________________

									 INVERT MATRIX Dll (or zero if restrained)
___________________________________________________________________________*/


		Dll [1] = 1/Dll[1];
		Dll [2] = 1/Dll[2];
		Dll [3] = 1/Dll[3];

/*__________________________________________________________________________

									 MATRIX K11
___________________________________________________________________________*/

  trimmKij (Kij, Cl1, Dll, Cl1);

		 sk [1][1] =  Kij[1][1];
		 sk [1][2] =  Kij[1][2];
		 sk [2][2] =  Kij[2][2];

/*__________________________________________________________________________

									 MATRIX K12
___________________________________________________________________________*/

  trimmKij (Kij, Cl1, Dll, Cl2);

		 sk [1][3] =  Kij[1][1];
		 sk [1][4] =  Kij[1][2];
   sk [2][3] =  Kij[2][1];
   sk [2][4] =  Kij[2][2];

/*__________________________________________________________________________

									 MATRIX K13
___________________________________________________________________________*/

  trimmKij (Kij, Cl1, Dll, Cl3);

		 sk [1][5] =  Kij[1][1];
		 sk [1][6] =  Kij[1][2];
   sk [2][5] =  Kij[2][1];
   sk [2][6] =  Kij[2][2];

/*__________________________________________________________________________

									 MATRIX K14
___________________________________________________________________________*/

  trimmKij (Kij, Cl1, Dll, Cl4);

		 sk [1][7] =  Kij[1][1];
		 sk [1][8] =  Kij[1][2];
   sk [2][7] =  Kij[2][1];
   sk [2][8] =  Kij[2][2];

/*__________________________________________________________________________

									 MATRIX K22
___________________________________________________________________________*/

  trimmKij (Kij, Cl2, Dll, Cl2);

		 sk [3][3] =  Kij[1][1];
		 sk [3][4] =  Kij[1][2];
		 sk [4][4] =  Kij[2][2];

/*__________________________________________________________________________

									 MATRIX K23
___________________________________________________________________________*/

  trimmKij (Kij, Cl2, Dll, Cl3);

		 sk [3][5] =  Kij[1][1];
		 sk [3][6] =  Kij[1][2];
   sk [4][5] =  Kij[2][1];
   sk [4][6] =  Kij[2][2];

/*__________________________________________________________________________

									 MATRIX K24
___________________________________________________________________________*/

  trimmKij (Kij, Cl2, Dll, Cl4);

		 sk [3][7] =  Kij[1][1];
		 sk [3][8] =  Kij[1][2];
   sk [4][7] =  Kij[2][1];
   sk [4][8] =  Kij[2][2];

/*__________________________________________________________________________

									 MATRIX K33
___________________________________________________________________________*/

  trimmKij (Kij, Cl3, Dll, Cl3);

		 sk [5][5] =  Kij[1][1];
		 sk [5][6] =  Kij[1][2];
		 sk [6][6] =  Kij[2][2];

/*__________________________________________________________________________

									 MATRIX K34
___________________________________________________________________________*/

  trimmKij (Kij, Cl3, Dll, Cl4);

		 sk [5][7] =  Kij[1][1];
		 sk [5][8] =  Kij[1][2];
   sk [6][7] =  Kij[2][1];
   sk [6][8] =  Kij[2][2];

/*__________________________________________________________________________

									 MATRIX K44
___________________________________________________________________________*/

  trimmKij (Kij, Cl4, Dll, Cl4);

		 sk [7][7] =  Kij[1][1];
		 sk [7][8] =  Kij[1][2];
		 sk [8][8] =  Kij[2][2];

/*___________________________________________________________________________

		EFFECTIVE STIFFNESS MATRIX GLOBAL COORDINATES
____________________________________________________________________________*/

  nrw  = 8;

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

			i1 = ndof*(ni[4] - 1);
			lm[7] = id[++i1];
			lm[8] = id[++i1];

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

}/* end of trimm2dis() */



/*******************************************[ USER LIBRARY : TRIM2DK.C  ]*****
 #                                                                           #
 #  FUNCTION :  trimm2dCl ()                                                 #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of side stiffness matrices C in global coordinates       #
 #      Only one element contribution.                                       #
 #                                                                           #
 #      Linear mixed triangular elements plane strain.                       #
 #      Element code = 3.                                                    #
 #                                                                           #
 #      Matrices are stored as submatrices Cli associted with node i         #
 #      Side direction clokwise from node 3 to 2. node 1 opposite node.      #
 #      Nodes are numbered anticlockwise.                                    #
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


 void trimm2dCl (double Cl1[4][3], double Cl2[4][3], double Cl3[4][3], double Dll[4],
                int n1, int n2, int n3, double E, double nu, double h, int nudim, double *xyz)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// double  Cl1[4][3];   /* submatrix Cl1 [3 x 2].                             */
// double  Cl2[4][3];   /* submatrix Cl2 [3 x 2].                             */
// double  Cl3[4][3];   /* submatrix Cl3 [3 x 2].                             */
// double  Dll[4];   /* submatrix diagonal [3 x 3].                        */
// int  n1, n2, n3;     /* node numbers of element                            */
// double   E   ;   /* elasticity modulus element                             */
// double   nu  ;   /* Poisson's ratio element                                */
// double   h   ;   /* thickness of element                                   */
// int     nudim;   /* number of spatial dimensions.                          */
// double *xyz  ;   /* nodal coordinates [numnp * nudim].                     */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

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

/*____________________________________________________________________________

		DECLARATION OF COMPUTED GEOMETRIC PROPERTIES
____________________________________________________________________________*/

 double  x31, y31, x23, y23, x12, y12;
 double  L1, L2, L3; /* length of sides                                     */
 double  darea   ;   /* two times the element area = x23*y31 - x31*y23 		   */
 double  c1, s1, c2, s2, c3, s3; /* director cosines of sides               */

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

 void trimmCij (double Cij[4][3], double Lj, double E, double nu, double h,
                double ci, double si, double cj, double sj);
 void trimmCii (double Cii[4][3], double Li, double ci, double si, double h);
 void trimmDll (double Dll[4], double darea, double E, double nu, double h);

/*____________________________________________________________________________

				 NODAL COORDINATES
____________________________________________________________________________*/

		i1 = nudim*(n1 - 1);
		x1 = xyz[++i1];
		y1 = xyz[++i1];

		i1 = nudim*(n2 - 1);
		x2 = xyz[++i1];
		y2 = xyz[++i1];

		i1 = nudim*(n3 - 1);
		x3 = xyz[++i1];
		y3 = xyz[++i1];

/*__________________________________________________________________________

			 ELEMENT CONSTANTS
___________________________________________________________________________*/

		x23 = x2 - x3;
		y23 = y2 - y3;

		x31 = x3 - x1;
		y31 = y3 - y1;

		x12 = x1 - x2;
		y12 = y1 - y2;

  L1 = sqrt(x23*x23 + y23*y23);
  L2 = sqrt(x31*x31 + y31*y31);
  L3 = sqrt(x12*x12 + y12*y12);

  c1 = x23/L1;
  s1 = y23/L1;
  c2 = x31/L2;
  s2 = y31/L2;
  c3 = x12/L3;
  s3 = y12/L3;

  darea = x23*y31 - x31*y23;

/*__________________________________________________________________________

									 MATRIX Dll
___________________________________________________________________________*/

  trimmDll (Dll, darea, E, nu, h);

/*__________________________________________________________________________

									 MATRIX Cl1
___________________________________________________________________________*/

  trimmCii (Cl1, L1, c1, s1, h);

/*__________________________________________________________________________

									 MATRIX Cl2
___________________________________________________________________________*/

  trimmCij (Cl2, L2, E, nu, h, c1, s1, c2, s2);

/*__________________________________________________________________________

									 MATRIX Cl3
___________________________________________________________________________*/

  trimmCij (Cl3, L3, E, nu, h, c1, s1, c3, s3);

/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trimm2dCl() */



/*******************************************[ USER LIBRARY : TRIM2DK.C  ]*****
 #                                                                           #
 #  FUNCTION :  trimmDll ()                                                  #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of diagonal submatrix Dll [3 x 3] in global coordinates. #
 #      Linear mixed triangular elements plane strain.                       #
 #      Element code = 3.                                                    #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED LIBRARY FUNCTIONS USED ]--------------------------------#
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

 void trimmDll (double Dll[4], double darea, double  E, double nu, double h)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Dll[4];  /* submatrix Cij [3 x 2].                                 */
// double   darea;  /* two times the element area 		                          */
// double       E;  /* elasticity modulus                                     */
// double      nu;  /* Poisson's ratio                                        */
// double   h   ;   /* thickness of element                                   */
/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*__________________________________________________________________________

			 DIAGONAL SUBMATRIX
___________________________________________________________________________*/

		Dll [1] = h*darea/6*E/(1-nu*nu);
		Dll [2] = h*darea/6*(1+nu)*(1-2*nu)/(1-nu)/E;
		Dll [3] = h*darea/6*2*(1+nu)/E;

/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trimmDll() */



/*******************************************[ USER LIBRARY : TRIM2DK.C  ]*****
 #                                                                           #
 #  FUNCTION :  trimmCij ()                                                  #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of stiffness submatrix Cij [3 x 2] in global coordinates.#
 #      ( for i<>j )                                                         #
 #      Linear mixed triangular elements plane strain.                       #
 #      Element code = 3.                                                    #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED LIBRARY FUNCTIONS USED ]--------------------------------#
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

 void trimmCij (double Cij[4][3], double Lj, double E, double nu, double h, double ci, double si, double cj, double sj)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Cij[4][3];   /* submatrix Cij [3 x 2].                                 */
// double     Lj;   /* length of side j                                       */
// double      E;   /* elasticity modulus                                     */
// double     nu;   /* Poisson's ratio                                        */
// double   h   ;   /* thickness of element                                   */
// double     ci;   /* cosine of angle alpha i of side i                      */
// double     si;   /* sine of angle alpha i of side i                        */
// double     cj;   /* cosine of angle alpha j of side j                      */
// double     sj;   /* sine of angle alpha j of side j                        */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 double  G1, G2, G3, G4; /* auxiliar coefficients                           */

/*____________________________________________________________________________

			DECLARATION OF MATRIX FACTORS AND COEFFICIENTS
____________________________________________________________________________*/

 double  f1;  /* factor Lj/6*E/(1-nu*nu)                                    */
 double  f2;  /* factor Lj/6/(1-nu)                                         */
 double  f3;  /* factor Lj/6                                                */

/*__________________________________________________________________________

			 ELEMENT CONSTANTS
___________________________________________________________________________*/

		f1 =  h*Lj/6*E/(1-nu*nu);
		f2 =  h*Lj/6/(1-nu);
		f3 =  h*Lj/6;

/*__________________________________________________________________________

									 AUXILIAR COEFFICIENTS
___________________________________________________________________________*/

		G1 = ci*sj-si*cj;
		G2 = si*sj+ci*cj;
		G3 = si*si-ci*ci;
		G4 = 2*ci*si;

/*____________________________________________________________________________

		 ELEMENT STIFFNESS SUBMATRIX MATRIX IN GLOBAL COORDINATES
____________________________________________________________________________*/


		 Cij[1][1] = f1*ci*G1;
		 Cij[1][2] = f1*si*G1;
		 Cij[2][1] = f2*( si*G2 - nu*( sj*G3 + cj*G4));
		 Cij[2][2] = f2*(-ci*G2 + nu*(-cj*G3 + sj*G4));
		 Cij[3][1] = f3*( cj*G3-sj*G4);
		 Cij[3][2] = f3*(-sj*G3-cj*G4);

/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trimmCij() */





/*******************************************[ USER LIBRARY : TRIM2DK.C  ]*****
 #                                                                           #
 #  FUNCTION :  trimmCii ()                                                  #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of stiffness submatrix Cii [3 x 2] in global coordinates.#
 #      Linear mixed triangular elements plane strain.                       #
 #      Element code = 3.                                                    #
 #                                                                           #
 #***************************************************************************/

void trimmCii (double Cii[4][3], double Li, double ci, double si, double h)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Cii[4][3];   /* submatrix Cii [3 x 2].                                 */
// double     Li;   /* length of side i                                       */
// double     ci;   /* cosine of angle alpha i of side i                      */
// double     si;   /* sine of angle alpha i of side i                        */
// double   h   ;   /* thickness of element                                   */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

/*____________________________________________________________________________

			DECLARATION OF MATRIX FACTORS AND COEFFICIENTS
____________________________________________________________________________*/

 double  f1;  /* factor Li/6                                                */

/*__________________________________________________________________________

			 ELEMENT CONSTANTS
___________________________________________________________________________*/

		f1 =  h*Li/6;

/*____________________________________________________________________________

		 ELEMENT STIFFNESS SUBMATRIX MATRIX IN GLOBAL COORDINATES
____________________________________________________________________________*/


		 Cii[1][1] = 0.0;
		 Cii[1][2] = 0.0;
		 Cii[2][1] = f1*si;
		 Cii[2][2] =-f1*ci;
		 Cii[3][1] =-f1*ci;
		 Cii[3][2] =-f1*si;

/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trimmCii() */


/*******************************************[ USER LIBRARY : TRIM2DK.C  ]*****
 #                                                                           #
 #  FUNCTION :  trimmKij ()                                                  #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of stiffness submatrix Kij [2 x 2] in global coordinates.#
 #                                                                           #
 #      Linear mixed triangular elements plane strain.                       #
 #      Element code = 3.                                                    #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED LIBRARY FUNCTIONS USED ]--------------------------------#
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

 void trimmKij (double Kij[3][3], double Cli[4][3], double *Dll, double Clj[4][3])

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Kij[3][3];   /* submatrix Kij [2 x 2].                                 */
// double  Cli[4][3];   /* submatrix Cli [3 x 2].                                 */
// double  Clj[4][3];   /* submatrix Clj [3 x 2].                                 */
// double   *Dll;   /* diadonal submatrix Dll [3 x 3].                        */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

		 ELEMENT STIFFNESS SUBMATRIX MATRIX IN GLOBAL COORDINATES
____________________________________________________________________________*/


		 Kij[1][1] = Cli[1][1]*Dll[1]*Clj[1][1] + Cli[2][1]*Dll[2]*Clj[2][1] + Cli[3][1]*Dll[3]*Clj[3][1];
		 Kij[1][2] = Cli[1][1]*Dll[1]*Clj[1][2] + Cli[2][1]*Dll[2]*Clj[2][2] + Cli[3][1]*Dll[3]*Clj[3][2];
		 Kij[2][1] = Cli[1][2]*Dll[1]*Clj[1][1] + Cli[2][2]*Dll[2]*Clj[2][1] + Cli[3][2]*Dll[3]*Clj[3][1];
		 Kij[2][2] = Cli[1][2]*Dll[1]*Clj[1][2] + Cli[2][2]*Dll[2]*Clj[2][2] + Cli[3][2]*Dll[3]*Clj[3][2];

/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trimmkKij() */




