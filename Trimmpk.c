
/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: TRIMMPK.C                                                  #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #        Routines for assembling of global stiffness matrix for             #
 #        Morley Plate triangular elements. Element code = 5.                #
 #        Side by side assembling                                            #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  trimmpk()    : computation of triangular element stiffness matrix. #
 #                                                                           #
 #****************************************************************************/


#include <math.h>

/*******************************************[ USER LIBRARY : TRIMMPK.C ]******
 #                                                                           #
 #  FUNCTION :  trimmpk ()                                                   #
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
 #      Morley Plate triangular elements.                                    #
 #      Element code = 5                                                     #
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

 void trimmpk (int sides, int bsides, int nudim, int *id, int *side_nodes,
               int *side_bc, int xnsn, double *xyz, int *side_mat, int *mat,
               double *prop, int numpr, double *sg, int *maxa, int *lm,
               double *sv)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// int     sides;   /* number of sides.                                       */
// int    bsides;   /* number of boundary sides.                              */
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

 void trimmpbs (int bsides , int nudim  , int *id      , int *side_nodes,
                 int *side_bc, int xnsn,
				 double *xyz, int *side_mat, int *mat       , double *prop,
                 int xnumpr , double *sg   , int *maxa      , int *lm,
                 double *sv);

 void trimmpis (int sides, int bsides     , int nudim  , int *id      ,
                 int *side_nodes, int *side_bc, int xnsn,
			     double *xyz, int *side_mat, int *mat       , double *prop,
                 int xnumpr , double *sg   , int *maxa      , int *lm,
                 double *sv);

/*__________________________________________________________________________

				LOOP ON BOUNDARY SIDES
___________________________________________________________________________*/

  trimmpbs (bsides, nudim, id, side_nodes, side_bc,
             xnsn, xyz, side_mat, mat, prop, numpr, sg, maxa, lm, sv);

/*__________________________________________________________________________

				LOOP ON INTERIOR SIDES
___________________________________________________________________________*/

  trimmpis (sides, bsides, nudim, id, side_nodes, side_bc,
             xnsn, xyz, side_mat, mat, prop, numpr, sg, maxa, lm, sv);

/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trimmpk() */



/*******************************************[ USER LIBRARY : TRIMMPK.C ]******
 #                                                                           #
 #  FUNCTION :  trimmpbs ()                                                  #
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
 #      Morley Plate triangular elements.                                    #
 #      Element code = 5                                                     #
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

 void trimmpbs (int bsides , int nudim  , int *id      , int *side_nodes,
                 int *side_bc, int xnsn,
				 double *xyz, int *side_mat, int *mat       , double *prop,
                 int numpr , double *sg   , int *maxa      , int *lm,
                 double *sv)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// int    bsides;   /* number of boundary sides.                              */
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
                  /* number of matrix rows = 9 for interior sides           */
 int     ic   ;   /* side counter for loop on sides.                        */
 int     ni[7];   /* side node numbers.                                     */
 int     isg[7];  /* side sign                                              */
 int     i, j ;   /* loop counters.                                         */
 int     i1   ;   /* auxiliar pointer for mesh arrays.                      */

/*____________________________________________________________________________

		  DECLARATION OF MATERIAL PROPERTIES
____________________________________________________________________________*/

 double  E    ;   /* elasticity modulus                                     */
 double  nu   ;   /* Poisson's ratio                                        */
 double  h    ;   /* thickness                                              */
 double  darea1  ;   /* two times the element area = x23*y31 - x31*y23      */

/*____________________________________________________________________________

			DECLARATION OF SIDE MATRICES
____________________________________________________________________________*/

 double sk[9][9];

 double Kw1e[4];
 double Kw2e[4];
 double Kw3e[4];
 double Krot4e[4];
 double Krot5e[4];
 double Krot6e[4];

 double Kw1[4];
 double Kw2[4];
 double Kw3[4];
 double Krot4[4];
 double Krot5[4];
 double Krot6[4];

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

 void addban (double *sg, int *maxa, double *sv, int *lm, int nrw);
 void trimmpKce (double Kw1e[4], double Kw2e[4], double Kw3e[4],
                 double Krot4e[4], double Krot5e[4], double Krot6e[4],
                 int n1, int n2, int n3, int nudim, double *xyz);
 double trimmpKij (double Kci[4], double E, double nu, double h, double Kcj[4]);
 double darea (int n1, int n2, int n3, int nudim, double* xyz);

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
		ni[4] = abs(side_nodes [++i1]);
		ni[5] = abs(side_nodes [++i1]);
		ni[6] = abs(side_nodes [++i1]);

        i1 = xnsn*(ic-1)+3;
        isg[4] = (side_nodes [++i1]>0)?1:-1;;
		isg[5] = (side_nodes [++i1]>0)?1:-1;;
		isg[6] = (side_nodes [++i1]>0)?1:-1;;

/*__________________________________________________________________________

				MATERIAL PROPERTIES OF ELEMENT
___________________________________________________________________________*/

		i1 = numpr * (side_mat[2*ic-1] - 1);

		E  = prop[++i1];
		nu = prop[++i1];
		h  = prop[++i1];
        darea1 = darea (ni[1], ni[2], ni[3],nudim, xyz);

/*__________________________________________________________________________

        MATRICES Kw1e, Kw2e, Kw3e, Krot4e, Krot5e, Krot6e FOR ELEMENT 1
___________________________________________________________________________*/


  trimmpKce (Kw1e, Kw2e, Kw3e, Krot4e, Krot5e, Krot6e, ni[3], ni[1], ni[2],
	         nudim, xyz);

/*__________________________________________________________________________

 GLOBAL CURVATURE VECTORS Kw1, Kw2, Kw3, Krot4, Krot5, Krot6 FOR ELEMENT 1
___________________________________________________________________________*/

	for ( i=1 ; i<=3 ; i++)
	{
      Kw1[i]= Kw2e[i];
      Kw2[i]= Kw3e[i];
      Kw3[i]= Kw1e[i];
      Krot4[i]= isg[4]*Krot4e[i];
      Krot5[i]= isg[5]*Krot6e[i];
      Krot6[i]= isg[6]*Krot5e[i];
    }
/*__________________________________________________________________________

			LOCAL STIFFNESS MATRIX SK (UPPER SYMMETRIC PART)
___________________________________________________________________________*/

   sk [1][1] =  trimmpKij (Kw1, E, nu, h, Kw1)*darea1/2;
   sk [1][2] =  trimmpKij (Kw1, E, nu, h, Kw2)*darea1/2;
   sk [1][3] =  trimmpKij (Kw1, E, nu, h, Kw3)*darea1/2;
   sk [1][4] =  trimmpKij (Kw1, E, nu, h, Krot4)*darea1/2;
   sk [1][5] =  trimmpKij (Kw1, E, nu, h, Krot5)*darea1/2;
   sk [1][6] =  trimmpKij (Kw1, E, nu, h, Krot6)*darea1/2;
   sk [2][2] =  trimmpKij (Kw2, E, nu, h, Kw2)*darea1/2;
   sk [2][3] =  trimmpKij (Kw2, E, nu, h, Kw3)*darea1/2;
   sk [2][4] =  trimmpKij (Kw2, E, nu, h, Krot4)*darea1/2;
   sk [2][5] =  trimmpKij (Kw2, E, nu, h, Krot5)*darea1/2;
   sk [2][6] =  trimmpKij (Kw2, E, nu, h, Krot6)*darea1/2;
   sk [3][3] =  trimmpKij (Kw3, E, nu, h, Kw3)*darea1/2;
   sk [3][4] =  trimmpKij (Kw3, E, nu, h, Krot4)*darea1/2;
   sk [3][5] =  trimmpKij (Kw3, E, nu, h, Krot5)*darea1/2;
   sk [3][6] =  trimmpKij (Kw3, E, nu, h, Krot6)*darea1/2;
   sk [4][4] =  trimmpKij (Krot4, E, nu, h, Krot4)*darea1/2;
   sk [4][5] =  trimmpKij (Krot4, E, nu, h, Krot5)*darea1/2;
   sk [4][6] =  trimmpKij (Krot4, E, nu, h, Krot6)*darea1/2;
   sk [5][5] =  trimmpKij (Krot5, E, nu, h, Krot5)*darea1/2;
   sk [5][6] =  trimmpKij (Krot5, E, nu, h, Krot6)*darea1/2;
   sk [6][6] =  trimmpKij (Krot6, E, nu, h, Krot6)*darea1/2;


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

	lm[1] = id[ni[1]];
    lm[2] = id[ni[2]];
    lm[3] = id[ni[3]];
    lm[4] = id[ni[4]];
    lm[5] = id[ni[5]];
    lm[6] = id[ni[6]];

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

}/* end of trimmpbs() */






/*******************************************[ USER LIBRARY : TRIMMPK.C ]******
 #                                                                           #
 #  FUNCTION :  trimmpis ()                                                  #
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
 #      Morley Plate triangular elements.                                    #
 #      Element code = 5                                                     #
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

 void trimmpis (int sides, int bsides     , int nudim  , int *id,
                 int *side_nodes, int *side_bc, int xnsn,
			     double *xyz, int *side_mat, int *mat       , double *prop,
                 int numpr , double *sg   , int *maxa      , int *lm,
                 double *sv)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// int     sides;   /* number of sides.                                       */
// int    bsides;   /* number of boundary sides.                              */
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
 int     ni[10];  /* side node numbers.                                     */
 int     isg[10]; /* side sign                                              */
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
 double  darea1  ;   /* two times the element1 area = x23*y31 - x31*y23     */
 double  darea2  ;   /* two times the element2 area = x23*y31 - x31*y23     */

/*____________________________________________________________________________

			DECLARATION OF SIDE MATRICES
____________________________________________________________________________*/

 double sk[10][10];

 double Kw1e[4];
 double Kw2e[4];
 double Kw3e[4];
 double Krot4e[4];
 double Krot5e[4];
 double Krot6e[4];

 double Kw1[4];
 double Kw2[4];
 double Kw3[4];
 double Kw7[4];
 double Krot4[4];
 double Krot5[4];
 double Krot6[4];
 double Krot8[4];
 double Krot9[4];

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

 void addban (double *sg, int *maxa, double *sv, int *lm, int nrw);
 void trimmpKce (double Kw1e[4], double Kw2e[4], double Kw3e[4],
                 double Krot4e[4], double Krot5e[4], double Krot6e[4],
                 int n1, int n2, int n3, int nudim, double *xyz);
 double trimmpKij (double Kci[4], double E, double nu, double h, double Kcj[4]);
 double darea (int n1, int n2, int n3, int nudim, double* xyz);
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
		ni[4] = abs(side_nodes [++i1]);
		ni[5] = abs(side_nodes [++i1]);
		ni[6] = abs(side_nodes [++i1]);
		ni[7] = side_nodes [++i1];
		ni[8] = abs(side_nodes [++i1]);
		ni[9] = abs(side_nodes [++i1]);

        i1 = xnsn*(ic-1)+3;
        isg[4] = (side_nodes [++i1]>0)?1:-1;
		isg[5] = (side_nodes [++i1]>0)?1:-1;
		isg[6] = (side_nodes [++i1]>0)?1:-1;
        isg[7] = (side_nodes [++i1]>0)?1:-1;
		isg[8] = (side_nodes [++i1]>0)?1:-1;
		isg[9] = (side_nodes [++i1]>0)?1:-1;

/*__________________________________________________________________________

				MATERIAL PROPERTIES OF ELEMENT 1
___________________________________________________________________________*/

		i1 = numpr * (side_mat[2*ic-1] - 1);

		E1  = prop[++i1];
		nu1 = prop[++i1];
		h1  = prop[++i1];

        darea1 = darea (ni[1], ni[2], ni[3], nudim, xyz);
        darea2 = darea (ni[1], ni[7], ni[2], nudim, xyz);

/*__________________________________________________________________________

        MATRICES Kw1e, Kw2e, Kw3e, Krot4e, Krot5e, Krot6e FOR ELEMENT 1
___________________________________________________________________________*/


  trimmpKce (Kw1e, Kw2e, Kw3e, Krot4e, Krot5e, Krot6e, ni[3], ni[1], ni[2],
	         nudim, xyz);

/*__________________________________________________________________________

 GLOBAL CURVATURE VECTORS Kw1, Kw2, Kw3, Krot4, Krot5, Krot6 FOR ELEMENT 1
___________________________________________________________________________*/

	for ( i=1 ; i<=3 ; i++)
	{
      Kw1[i]= Kw2e[i];
      Kw2[i]= Kw3e[i];
      Kw3[i]= Kw1e[i];
      Krot4[i]= isg[4]*Krot4e[i];
      Krot5[i]= isg[5]*Krot6e[i];
      Krot6[i]= isg[6]*Krot5e[i];
    }

/*__________________________________________________________________________

	 LOCAL STIFFNESS MATRIX SK (UPPER SYMMETRIC PART) FOR ELEMENT 1
___________________________________________________________________________*/

   sk [1][1] =  trimmpKij (Kw1, E1, nu1, h1, Kw1)*darea1/2;
   sk [1][2] =  trimmpKij (Kw1, E1, nu1, h1, Kw2)*darea1/2;
   sk [1][3] =  trimmpKij (Kw1, E1, nu1, h1, Kw3)*darea1/2;
   sk [1][4] =  trimmpKij (Kw1, E1, nu1, h1, Krot4)*darea1/2;
   sk [1][5] =  trimmpKij (Kw1, E1, nu1, h1, Krot5)*darea1/2;
   sk [1][6] =  trimmpKij (Kw1, E1, nu1, h1, Krot6)*darea1/2;
   sk [1][7] =  0;
   sk [1][8] =  0;
   sk [1][9] =  0;
   sk [2][2] =  trimmpKij (Kw2, E1, nu1, h1, Kw2)*darea1/2;
   sk [2][3] =  trimmpKij (Kw2, E1, nu1, h1, Kw3)*darea1/2;
   sk [2][4] =  trimmpKij (Kw2, E1, nu1, h1, Krot4)*darea1/2;
   sk [2][5] =  trimmpKij (Kw2, E1, nu1, h1, Krot5)*darea1/2;
   sk [2][6] =  trimmpKij (Kw2, E1, nu1, h1, Krot6)*darea1/2;
   sk [2][7] =  0;
   sk [2][8] =  0;
   sk [2][9] =  0;
   sk [3][3] =  trimmpKij (Kw3, E1, nu1, h1, Kw3)*darea1/2;
   sk [3][4] =  trimmpKij (Kw3, E1, nu1, h1, Krot4)*darea1/2;
   sk [3][5] =  trimmpKij (Kw3, E1, nu1, h1, Krot5)*darea1/2;
   sk [3][6] =  trimmpKij (Kw3, E1, nu1, h1, Krot6)*darea1/2;
   sk [3][7] =  0;
   sk [3][8] =  0;
   sk [3][9] =  0;
   sk [4][4] =  trimmpKij (Krot4, E1, nu1, h1, Krot4)*darea1/2;
   sk [4][5] =  trimmpKij (Krot4, E1, nu1, h1, Krot5)*darea1/2;
   sk [4][6] =  trimmpKij (Krot4, E1, nu1, h1, Krot6)*darea1/2;
   sk [4][7] =  0;
   sk [4][8] =  0;
   sk [4][9] =  0;
   sk [5][5] =  trimmpKij (Krot5, E1, nu1, h1, Krot5)*darea1/2;
   sk [5][6] =  trimmpKij (Krot5, E1, nu1, h1, Krot6)*darea1/2;
   sk [5][7] =  0;
   sk [5][8] =  0;
   sk [5][9] =  0;
   sk [6][6] =  trimmpKij (Krot6, E1, nu1, h1, Krot6)*darea1/2;
   sk [6][7] =  0;
   sk [6][8] =  0;
   sk [6][9] =  0;
   sk [7][7] =  0;
   sk [7][8] =  0;
   sk [7][9] =  0;
   sk [8][8] =  0;
   sk [8][9] =  0;
   sk [9][9] =  0;


/*__________________________________________________________________________

				MATERIAL PROPERTIES OF ELEMENT 2
___________________________________________________________________________*/

		i1 = numpr * (side_mat[2*ic] - 1);

		E2  = prop[++i1];
		nu2 = prop[++i1];
		h2  = prop[++i1];

/*__________________________________________________________________________

        MATRICES Kw1e, Kw2e, Kw3e, Krot4e, Krot5e, Krot6e FOR ELEMENT 2
___________________________________________________________________________*/


  trimmpKce (Kw1e, Kw2e, Kw3e, Krot4e, Krot5e, Krot6e, ni[7], ni[2], ni[1],
	         nudim, xyz);

/*__________________________________________________________________________

 GLOBAL CURVATURE VECTORS Kw1, Kw2, Kw7, Krot4, Krot8, Krot9 FOR ELEMENT 2
___________________________________________________________________________*/


  for ( i=1 ; i<=3 ; i++)
  {
      Kw1[i] = Kw3e[i];
      Kw2[i] = Kw2e[i];
      Kw7[i] = Kw1e[i];
      Krot4[i] = -isg[4]*Krot4e[i];
      Krot8[i] =  isg[8]*Krot6e[i];
      Krot9[i] =  isg[9]*Krot5e[i];
  }

/*__________________________________________________________________________

	   LOCAL STIFFNESS MATRIX SK (UPPER SYMMETRIC PART) for ELEMENT 2
___________________________________________________________________________*/

   sk [1][1] +=  trimmpKij (Kw1, E2, nu2, h2, Kw1)*darea2/2;
   sk [1][2] +=  trimmpKij (Kw1, E2, nu2, h2, Kw2)*darea2/2;
   sk [1][3] +=  0;
   sk [1][4] +=  trimmpKij (Kw1, E2, nu2, h2, Krot4)*darea2/2;
   sk [1][5] +=  0;
   sk [1][6] +=  0;
   sk [1][7] +=  trimmpKij (Kw1, E2, nu2, h2, Kw7)*darea2/2;
   sk [1][8] +=  trimmpKij (Kw1, E2, nu2, h2, Krot8)*darea2/2;
   sk [1][9] +=  trimmpKij (Kw1, E2, nu2, h2, Krot9)*darea2/2;
   sk [2][2] +=  trimmpKij (Kw2, E2, nu2, h2, Kw2)*darea2/2;
   sk [2][3] +=  0;
   sk [2][4] +=  trimmpKij (Kw2, E2, nu2, h2, Krot4)*darea2/2;
   sk [2][5] +=  0;
   sk [2][6] +=  0;
   sk [2][7] +=  trimmpKij (Kw2, E2, nu2, h2, Kw7)*darea2/2;
   sk [2][8] +=  trimmpKij (Kw2, E2, nu2, h2, Krot8)*darea2/2;
   sk [2][9] +=  trimmpKij (Kw2, E2, nu2, h2, Krot9)*darea2/2;
   sk [3][3] +=  0;
   sk [3][4] +=  0;
   sk [3][5] +=  0;
   sk [3][6] +=  0;
   sk [3][7] +=  0;
   sk [3][8] +=  0;
   sk [3][9] +=  0;
   sk [4][4] +=  trimmpKij (Krot4, E2, nu2, h2, Krot4)*darea2/2;
   sk [4][5] +=  0;
   sk [4][6] +=  0;
   sk [4][7] +=  trimmpKij (Krot4, E2, nu2, h2, Kw7)*darea2/2;
   sk [4][8] +=  trimmpKij (Krot4, E2, nu2, h2, Krot8)*darea2/2;
   sk [4][9] +=  trimmpKij (Krot4, E2, nu2, h2, Krot9)*darea2/2;
   sk [5][5] +=  0;
   sk [5][6] +=  0;
   sk [5][7] +=  0;
   sk [5][8] +=  0;
   sk [5][9] +=  0;
   sk [6][6] +=  0;
   sk [6][7] +=  0;
   sk [6][8] +=  0;
   sk [6][9] +=  0;
   sk [7][7] +=  trimmpKij (Kw7, E2, nu2, h2, Kw7)*darea2/2;
   sk [7][8] +=  trimmpKij (Kw7, E2, nu2, h2, Krot8)*darea2/2;
   sk [7][9] +=  trimmpKij (Kw7, E2, nu2, h2, Krot9)*darea2/2;
   sk [8][8] +=  trimmpKij (Krot8, E2, nu2, h2, Krot8)*darea2/2;
   sk [8][9] +=  trimmpKij (Krot8, E2, nu2, h2, Krot9)*darea2/2;
   sk [9][9] +=  trimmpKij (Krot9, E2, nu2, h2, Krot9)*darea2/2;

/*___________________________________________________________________________

		EFFECTIVE STIFFNESS MATRIX GLOBAL COORDINATES
____________________________________________________________________________*/

  nrw  = 9;

  i1 = 0;
  for ( i=1 ; i<=nrw ; i++)
   	 for( j=i ; j<=nrw ; j++)  sv[++i1] = sk[i][j];

/*__________________________________________________________________________

		EQUATION NUMBERS OF ELEMENT DOF
___________________________________________________________________________*/

  	lm[1] = id[ni[1]];
    lm[2] = id[ni[2]];
    lm[3] = id[ni[3]];
    lm[4] = id[ni[4]];
    lm[5] = id[ni[5]];
    lm[6] = id[ni[6]];
    lm[7] = id[ni[7]];
    lm[8] = id[ni[8]];
    lm[9] = id[ni[9]];

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

}/* end of trimmpis() */



/*******************************************[ USER LIBRARY : TRIMMPK.C ]******
 #                                                                           #
 #  FUNCTION :  trimmpKce()                                                  #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of curvature vectors Kw, Krot in global coordinates      #
 #      Only one element contribution.                                       #
 #                                                                           #
 #      Morley Plate triangular elements.                                    #
 #      Element code = 5                                                     #
 #                                                                           #
 #      Matrices are stored as submatrices Cli associted with node i         #
 #      Side direction clokwise from node v3 to v2. node v1 opposite node.   #
 #      Nodes are numbered anticlockwise.                                    #
 #      Interaction side v2 - v1 with triangle v2, v1, v3.
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

 void trimmpKce (double Kw1e[4], double Kw2e[4], double Kw3e[4],
                 double Krot4e[4], double Krot5e[4], double Krot6e[4],
                 int n1, int n2, int n3, int nudim, double *xyz)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// double Kw1e[4];   /* vector Kw1 [3 x 1].                                   */
// double Kw2e[4];   /* vector Kw2 [3 x 1].                                   */
// double Kw3e[4];   /* vector Kw3 [3 x 1].                                   */
// double Krot4e[4]; /* vector Krot4 [3 x 1].                                 */
// double Krot5e[4]; /* vector Krot5 [3 x 1].                                 */
// double Krot6e[4]; /* vector Krot6 [3 x 1].                                 */
// int  n1, n2, n3;     /* global node numbers of element vertexes            */
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

 double  L1, L2, L3; /* length of sides                                     */
 double  L12, L22, L32; /* square length of sides                           */
 double  darea   ;   /* two times the element area = x23*y31 - x31*y23      */
 double  b1, b2, b3, c1, c2, c3;           /* geometric constants           */
 double  b12, b22, b32, c12, c22, c32;     /* geometric constants           */
 double  bc12, bc13, bc23, bc1, bc2, bc3;  /* geometric constants           */

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

void trimmpKwi(double Kc[4], double b12, double b22, double b32, double c12,
                double c22, double c32, double bc12, double bc13, double bc1,
                double bc2, double bc3, double L22, double L32, double darea);

void trimmpKroti(double Kc[4], double b32, double c32, double bc3, double L3,
                  double darea);

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

	b1 = y2 - y3;
    b2 = y3 - y1;
    b3 = y1 - y2;

	c1 = x3 - x2;
	c2 = x1 - x3;
    c3 = x2 - x1;

    b12 = b1*b1;
    b22 = b2*b2;
    b32 = b3*b3;

    c12 = c1*c1;
    c22 = c2*c2;
    c32 = c3*c3;

    darea = b1*c2-b2*c1;

    L12 = c12 + b12;
    L22 = c22 + b22;
    L32 = c32 + b32;

    L1 = sqrt(L12);
    L2 = sqrt(L22);
    L3 = sqrt(L32);

    bc12 = b1*b2 + c1*c2;
    bc13 = b1*b3 + c1*c3;
    bc23 = b2*b3 + c2*c3;

    bc1 = b1*c1;
    bc2 = b2*c2;
    bc3 = b3*c3;

/*__________________________________________________________________________

									 MATRIX Kw1
___________________________________________________________________________*/

  trimmpKwi(Kw1e, b12, b22, b32, c12, c22, c32, bc12, bc13, bc1, bc2, bc3,
            L22, L32, darea);

/*__________________________________________________________________________

									 MATRIX Kw2
___________________________________________________________________________*/

  trimmpKwi(Kw2e, b22, b32, b12, c22, c32, c12, bc23, bc12, bc2, bc3, bc1,
             L32, L12, darea);

/*__________________________________________________________________________

									 MATRIX Kw3
___________________________________________________________________________*/

  trimmpKwi(Kw3e, b32, b12, b22, c32, c12, c22, bc13, bc23, bc3, bc1, bc2,
             L12, L22, darea);

/*__________________________________________________________________________

									 MATRIX Krot4
___________________________________________________________________________*/

  trimmpKroti(Krot4e, b12, c12, bc1, L1, darea);

/*__________________________________________________________________________

									 MATRIX Krot5
___________________________________________________________________________*/

  trimmpKroti(Krot5e, b22, c22, bc2, L2, darea);

/*__________________________________________________________________________

									 MATRIX Krot6
___________________________________________________________________________*/

  trimmpKroti(Krot6e, b32, c32, bc3, L3, darea);

/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trimmpKce() */



/*******************************************[ USER LIBRARY : TRIMEMPK.C  ]*****
 #                                                                           #
 #  FUNCTION :  trimmpKwi ()                                                 #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of global curvatures vector Kc [3] in global coordinates.#
 #      Only vertex displacements shape functions                            #
 #      Morley Plate triangular elements.                                    #
 #      Element code = 5                                                     #
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


void trimmpKwi(double  Kc[4], double  bi2, double bj2, double bk2,
                 double ci2, double cj2, double ck2, double  bcij,
                 double bcik, double bci, double bcj, double bck, double  Lj2,
                 double Lk2, double  darea)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Kc[4];        /* vector of global curvatures                       */
// double  bi2, bj2, bk2, ci2, cj2, ck2;     /* geometric constants           */
// double  bcij, bcik, bci, bcj, bck;        /* geometric constants           */
// double  Lj2, Lk2;   /* square length of sides                              */
// double  darea   ;   /* two times the element area = x23*y31 - x31*y23      */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

		        DISPLACEMENT SHAPE FUNCTION CURVATURES
____________________________________________________________________________*/


		 Kc[1] = 2*(bi2 + bj2*bcij/Lj2 + bk2*bcik/Lk2)/darea/darea;
         Kc[2] = 2*(ci2 + cj2*bcij/Lj2 + ck2*bcik/Lk2)/darea/darea;
         Kc[3] = 4*(bci + bcj*bcij/Lj2 + bck*bcik/Lk2)/darea/darea;

/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trimmpKwi() */




/*******************************************[ USER LIBRARY : TRIMEMPK.C  ]*****
 #                                                                           #
 #  FUNCTION :  trimempKroti ()                                              #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of global curvatures vector Kc [3] in global coordinates.#
 #      Only midside rotations shape functions                               #
 #      Enhanced Morley Plate triangular elements.                           #
 #      Element code = 4                                                     #
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

void trimmpKroti(double  Kc[4], double  bi2, double  ci2, double  bci, double  Li,
                 double  darea)


/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Kc[4];         /* vector of global curvatures                      */
// double  bi2, ci2, bci; /* geometric constants                              */
// double  Li;            /* length of side i                                 */
// double  darea   ;   /* two times the element area = x23*y31 - x31*y23      */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

		        DISPLACEMENT SHAPE FUNCTION CURVATURES
____________________________________________________________________________*/


		 Kc[1] = - 2*(bi2/Li)/darea;
         Kc[2] = - 2*(ci2/Li)/darea;
         Kc[3] = - 4*(bci/Li)/darea;


/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trimmpKroti() */


/******************************************[ USER LIBRARY : TRIMMPK.C  ]******
 #                                                                           #
 #  FUNCTION :  trimmpKij ()                                                 #
 #                                                                           #
 #  TYPE     :  double                                                       #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of stiffness coefficient Kij in global coordinates.      #
 #                                                                           #
 #      Morley Plate triangular elements.                                    #
 #      Element code = 5                                                     #
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
 double trimmpKij (double Kci[4], double E, double nu, double h, double Kcj[4])

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Kci[4];    /* submatrix Kci [3 x 1].                               */
// double  Kcj[4];    /* submatrix Kcj [3 x 1].                               */
// double  E    ;   /* elasticity modulus                                     */
// double  nu   ;   /* Poisson's ratio                                        */
// double  h    ;   /* thickness                                              */
/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 double  Kij ;   /* coefficient Kij                                        */
 double  D[4][4]; /* constitutive matrix                                   */
 double  Da;
/*____________________________________________________________________________

		             CONSTITUTIVE MATRIX
____________________________________________________________________________*/

 Da=E*h*h*h/(1-nu*nu)/12;
 D[1][1]=Da;
 D[1][2]=Da*nu;
 D[1][3]=0;
 D[2][1]=Da*nu;
 D[2][2]=Da;
 D[2][3]=0;
 D[3][1]=0;
 D[3][2]=0;
 D[3][3]=Da*(1-nu)/2;

/*____________________________________________________________________________

		 ELEMENT STIFFNESS SUBMATRIX MATRIX IN GLOBAL COORDINATES
____________________________________________________________________________*/


  Kij = Kci[1]*(D[1][1]*Kcj[1]+D[1][2]*Kcj[2]+D[1][3]*Kcj[3]) +
        Kci[2]*(D[2][1]*Kcj[1]+D[2][2]*Kcj[2]+D[2][3]*Kcj[3]) +
        Kci[3]*(D[3][1]*Kcj[1]+D[3][2]*Kcj[2]+D[3][3]*Kcj[3]);

  Kij = Kij/3;

/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

  return (Kij);

}/* end of trimmpKij() */

/*******************************************[ USER LIBRARY : TRIMMPK.C ]******
 #                                                                           #
 #  FUNCTION :  darea()                                                      #
 #                                                                           #
 #  TYPE     :  double                                                       #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of two times the area of triangle with nodes             #
 #      n1, n2, n3.                                                          #
 #                                                                           #
 #***************************************************************************/


 double darea (int n1, int n2, int n3,int nudim, double *xyz)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// int  n1, n2, n3;     /* global node numbers of element vertexes            */
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

 double  b1, b2, b3, c1, c2, c3;           /* geometric constants           */

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

	b1 = y2 - y3;
    b2 = y3 - y1;
    b3 = y1 - y2;

	c1 = x3 - x2;
	c2 = x1 - x3;
    c3 = x2 - x1;

    return (b1*c2-b2*c1);

/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of darea() */

