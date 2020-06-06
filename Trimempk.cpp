/**********************[          Trimempk.cpp            ]********************/
#include "Trimempk.h"


/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: TRIMEMPK.C                                                 #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #        Routines for assembling of global stiffness matrix for             #
 #        Enhanced Morley Plate triangular elements. Element code = 4.       #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  trimempk()    : computation of triangular element stiffness matrix.#
 #                                                                           #
 #****************************************************************************/


/*******************************************[ USER LIBRARY : TRIMEMPK.C ]*****
 #                                                                           #
 #  FUNCTION :  trimempk ()                                                  #
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
 #      Enhanced Morley Plate triangular elements.                           #
 #      Element code = 4                                                     #
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


  void trimempk (int sides, int bsides, int nudim, int *id, int *side_nodes,
                int *side_bc, int xnsn, double *xyz, int *side_mat, int *mat,
                double *prop, int numpr, double *sg, int *maxa, int *lm,
                double *sv, double *fp, double *fm)

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
// double *fp   ;   /* concentrated nodal loads coarse mesh only.             */
// double *fm   ;   /* concentrated mixed nodal loads coarse mesh only.       */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

 void trimempbs (int bsides , int nudim  , int *id      , int *side_nodes,
                 int *side_bc, int xnsn,
				 double *xyz, int *side_mat, int *mat       , double *prop,
                 int xnumpr , double *sg   , int *maxa      , int *lm,
                 double *sv, double *fp, double *fm);

 void trimempis (int sides, int bsides     , int nudim  , int *id      ,
                 int *side_nodes, int *side_bc, int xnsn,
			     double *xyz, int *side_mat, int *mat       , double *prop,
                 int xnumpr , double *sg   , int *maxa      , int *lm,
                 double *sv);
/*__________________________________________________________________________

				LOOP ON BOUNDARY SIDES
___________________________________________________________________________*/

  trimempbs (bsides, nudim, id, side_nodes, side_bc, xnsn, xyz, side_mat,
             mat, prop, numpr, sg, maxa, lm, sv, fp, fm);


/*__________________________________________________________________________

				LOOP ON INTERIOR SIDES
___________________________________________________________________________*/

  trimempis (sides, bsides, nudim, id, side_nodes, side_bc, xnsn,
             xyz, side_mat, mat, prop, numpr, sg, maxa, lm, sv);


/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trimempk() */



/*******************************************[ USER LIBRARY : TRIMEMPK.C ]*****
 #                                                                           #
 #  FUNCTION :  trimempbs ()                                                 #
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
 #      Enhanced Morley Plate triangular elements.                           #
 #      Element code = 4                                                     #
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


void trimempbs (int bsides, int nudim, int *id, int *side_nodes,
                 int *side_bc, int xnsn,
				 double *xyz, int *side_mat, int *mat, double *prop,
                 int numpr, double *sg, int *maxa, int *lm,
                 double *sv, double *fp, double *fm)

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
// double *fp   ;   /* concentrated nodal loads coarse mesh only.             */
// double *fm   ;   /* concentrated mixed nodal loads coarse mesh only.       */

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

/*____________________________________________________________________________

		  DECLARATION OF NODAL COORDINATES AND SIDE LENGHT
____________________________________________________________________________*/

 double  x1   ;   /* X coordinate of node 1.                                */
 double  x2   ;   /* X coordinate of node 2.                                */
 double  y1   ;   /* Y coordinate of node 1.                                */
 double  y2   ;   /* Y coordinate of node 2.                                */
 double  L    ;   /* lenght of boundary side between nodes 1 and 2.         */

/*____________________________________________________________________________

			DECLARATION OF SIDE MATRICES
____________________________________________________________________________*/

 double sk[9][9];
 double Cl1e[4];
 double Cl2e[4];
 double Cl3e[4];
 double Cl4e[4];
 double Cl5e[4];
 double Cl6e[4];
 double Ae[4];
 double Cl1[4];
 double Cl2[4];
 double Cl3[4];
 double Cl4[4];
 double Cl5[4];
 double Cl6[4];
 double All[4];

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

 void addban (double *sg, int *maxa, double *sv, int *lm, int nrw);
 void trimempCle (double Cl1e[4], double Cl2e[4], double Cl3e[4], double Cl4e[4],
	              double Cl5e[4], double Cl6e[4],double Ae[4],
                  int n1, int n2, int n3, double E, double nu, double h, int nudim,
                  double *xyz);
 double trimempKij (double Cli[4], double All[4], double Clj[4]);

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

/*________________________________________________________________________

                        LENGHT OF BOUNDARY SIDE
 __________________________________________________________________________*/

        i1 = nudim*(ni[1] - 1);
	    x1 = xyz[++i1];
	    y1 = xyz[++i1];

	    i1 = nudim*(ni[2] - 1);
	    x2 = xyz[++i1];
	    y2 = xyz[++i1];

         L = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
/*__________________________________________________________________________

				MATERIAL PROPERTIES OF ELEMENT
___________________________________________________________________________*/

		i1 = numpr * (side_mat[2*ic-1] - 1);

		E  = prop[++i1];
		nu = prop[++i1];
		h  = prop[++i1];

/*__________________________________________________________________________

        MATRICES Cl1e, Cl2e, Cl3e, Cl4e, Cl5e, Cl6e, Ae FOR ELEMENT 1
___________________________________________________________________________*/


  trimempCle (Cl1e, Cl2e, Cl3e, Cl4e, Cl5e, Cl6e, Ae, ni[3], ni[1], ni[2],
	          E, nu, h, nudim, xyz);

/*__________________________________________________________________________

     ADD TO SIDE MATRICES Cl1, Cl2, Cl3, Cl4, Cl5, Cl6, Dll FOR ELEMENT 1
___________________________________________________________________________*/

	for ( i=1 ; i<=3 ; i++)
	{
      All[i]= Ae[i];
      Cl1[i]= Cl2e[i];
      Cl2[i]= Cl3e[i];
      Cl3[i]= Cl1e[i];
      Cl4[i]= isg[4]*Cl4e[i];
      Cl5[i]= isg[5]*Cl6e[i];
      Cl6[i]= isg[6]*Cl5e[i];
    }

/*__________________________________________________________________________

     CORRECTIONS TO NODAL FORCES DUE TO PRESCRIBED MIXED VARIABLES
                (only NORMAL MOMENT considered)
___________________________________________________________________________*/

if (side_bc[2*ic]==0)
  {
    fm[ni[1]]-=Cl1[2]*isg[4]*fp[ni[4]]/L;
    fm[ni[2]]-=Cl2[2]*isg[4]*fp[ni[4]]/L;
    fm[ni[3]]-=Cl3[2]*isg[4]*fp[ni[4]]/L;
    fm[ni[4]]-=Cl4[2]*isg[4]*fp[ni[4]]/L;
    fm[ni[5]]-=Cl5[2]*isg[4]*fp[ni[4]]/L;
    fm[ni[6]]-=Cl6[2]*isg[4]*fp[ni[4]]/L;
  }

/*__________________________________________________________________________

	  INVERT MATRIX All (or zero if tractions are prescribed)

   w prescribed           ==> Ks  prescribed A[1]=0
   normal rot. prescribed ==> Mns prescribed A[3]=0
   w free                 ==> Mns prescribed A[3]=0
   normal rot. free       ==> Mn  prescribed A[2]=0

   apoyo simple (1,0)
     w  =0     ==> Ks =0  ==> A[1]=0
     Mn =0     ==> Mn =0  ==> A[2]=0
   apoyo empotrado (1,1)
     w  =0     ==> Ks =0  ==> A[1]=0
    rotn=0     ==> Mns=0  ==> A[3]=0
   apoyo simetrico (0,1)
    rotn=0     ==> Mns=0  ==> A[3]=0

    funcionando
    empotrada
              side_bc[2*ic-1]==1
              side_bc[2*ic]  ==1
              A11[1]=0, A11[2]=1/A11[2], A11[3]=1/A11[3]
               Ks=0   ,   Mn<>0,     Mns<>0

    funcionando
    simetria
              side_bc[2*ic-1]==0
              side_bc[2*ic]  ==1
              A11[1]=1/A11[1], A11[2]=1/A11[2], A11[3]=1/A11[3]
                  Ks<>0      ,   Mn<>0        ,     Mns<>0

    no funcionando
    SA:       A11[1]=0, A11[2]=0, A11[3]=1/A11[3]
               Ks=0   ,    Mn=0,     Mns<>0
___________________________________________________________________________*/

//    All [1] = 1/All[1];
//    All [2] = 1/All[2];
//    All [3] = 1/All[3];

    All [1] = (side_bc[2*ic-1]==1)? 0: 1/All[1];
    All [2] = (side_bc[2*ic  ]==0)? 0: 1/All[2];
    All [3] = (side_bc[2*ic  ]==1)? 0: 1/All[3];

/*__________________________________________________________________________

			LOCAL STIFFNESS MATRIX SK (UPPER SYMMETRIC PART)
___________________________________________________________________________*/

   sk [1][1] =  trimempKij (Cl1, All, Cl1);
   sk [1][2] =  trimempKij (Cl1, All, Cl2);
   sk [1][3] =  trimempKij (Cl1, All, Cl3);
   sk [1][4] =  trimempKij (Cl1, All, Cl4);
   sk [1][5] =  trimempKij (Cl1, All, Cl5);
   sk [1][6] =  trimempKij (Cl1, All, Cl6);
   sk [2][2] =  trimempKij (Cl2, All, Cl2);
   sk [2][3] =  trimempKij (Cl2, All, Cl3);
   sk [2][4] =  trimempKij (Cl2, All, Cl4);
   sk [2][5] =  trimempKij (Cl2, All, Cl5);
   sk [2][6] =  trimempKij (Cl2, All, Cl6);
   sk [3][3] =  trimempKij (Cl3, All, Cl3);
   sk [3][4] =  trimempKij (Cl3, All, Cl4);
   sk [3][5] =  trimempKij (Cl3, All, Cl5);
   sk [3][6] =  trimempKij (Cl3, All, Cl6);
   sk [4][4] =  trimempKij (Cl4, All, Cl4);
   sk [4][5] =  trimempKij (Cl4, All, Cl5);
   sk [4][6] =  trimempKij (Cl4, All, Cl6);
   sk [5][5] =  trimempKij (Cl5, All, Cl5);
   sk [5][6] =  trimempKij (Cl5, All, Cl6);
   sk [6][6] =  trimempKij (Cl6, All, Cl6);


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

}/* end of trimempbs() */






/*******************************************[ USER LIBRARY : TRIMEMPK.C ]*****
 #                                                                           #
 #  FUNCTION :  trimempis ()                                                 #
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
 #      Enhanced Morley Plate triangular elements.                           #
 #      Element code = 4                                                     #
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


 void trimempis (int sides, int bsides, int nudim, int *id,
                 int *side_nodes, int *side_bc, int xnsn,
			     double *xyz, int *side_mat, int *mat, double *prop,
                 int numpr, double *sg, int *maxa, int *lm,
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


/*____________________________________________________________________________

			DECLARATION OF SIDE MATRICES
____________________________________________________________________________*/

 double sk[10][10];
 double Cl1e[4];
 double Cl2e[4];
 double Cl3e[4];
 double Cl4e[4];
 double Cl5e[4];
 double Cl6e[4];
 double Ae[4];
 double Cl1[4];
 double Cl2[4];
 double Cl3[4];
 double Cl4[4];
 double Cl5[4];
 double Cl6[4];
 double Cl7[4];
 double Cl8[4];
 double Cl9[4];
 double All[4];

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

 void addban (double *sg, int *maxa, double *sv, int *lm, int nrw);
 void trimempCle (double Cl1e[4], double Cl2e[4], double Cl3e[4], double Cl4e[4],
	              double Cl5e[4], double Cl6e[4],double Ae[4],
                  int n1, int n2, int n3, double E, double nu, double h, int nudim,
                  double *xyz);
 double trimempKij (double Cli[4], double All[4], double Clj[4]);
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

/*__________________________________________________________________________

        MATRICES Cl1e, Cl2e, Cl3e, Cl4e, Cl5e, Cl6e, Ae FOR ELEMENT 1
___________________________________________________________________________*/


  trimempCle (Cl1e, Cl2e, Cl3e, Cl4e, Cl5e, Cl6e, Ae, ni[3], ni[1], ni[2],
	          E1, nu1, h1, nudim, xyz);

/*__________________________________________________________________________

     ADD TO SIDE MATRICES Cl1, Cl2, Cl3, Cl4, Cl5, Cl6, Dll FOR ELEMENT 1
___________________________________________________________________________*/

	for ( i=1 ; i<=3 ; i++)
	{
      All[i]=Ae[i];
      Cl1[i]= Cl2e[i];
      Cl2[i]= Cl3e[i];
      Cl3[i]= Cl1e[i];
      Cl4[i]= isg[4]*Cl4e[i];
      Cl5[i]= isg[5]*Cl6e[i];
      Cl6[i]= isg[6]*Cl5e[i];
    }
/*__________________________________________________________________________

				MATERIAL PROPERTIES OF ELEMENT 2
___________________________________________________________________________*/

		i1 = numpr * (side_mat[2*ic] - 1);

		E2  = prop[++i1];
		nu2 = prop[++i1];
		h2  = prop[++i1];

/*__________________________________________________________________________

        MATRICES Cl1e, Cl2e, Cl3e, Cl4e, Cl5e, Cl6e, Ae FOR ELEMENT 2
___________________________________________________________________________*/


  trimempCle (Cl1e, Cl2e, Cl3e, Cl4e, Cl5e, Cl6e, Ae, ni[7], ni[2], ni[1],
	          E2, nu2, h2, nudim, xyz);

/*__________________________________________________________________________

     ADD TO SIDE MATRICES Cl1, Cl2, Cl3, Cl4, Cl5, Cl6, Dll FOR ELEMENT 2
___________________________________________________________________________*/

  for ( i=1 ; i<=3 ; i++)
  {
      All[i]+= Ae[i];
      Cl1[i]+= Cl3e[i];
      Cl2[i]+= Cl2e[i];
      Cl4[i]-= isg[4]*Cl4e[i];
      Cl7[i] = Cl1e[i];
      Cl8[i] = isg[8]*Cl6e[i];
      Cl9[i] = isg[9]*Cl5e[i];
  }

/*__________________________________________________________________________

		 INVERT MATRIX Dll (assumed not restrained in interior sides)
___________________________________________________________________________*/


   All [1] = 1/All[1];
   All [2] = 1/All[2];
   All [3] = 1/All[3];

/*__________________________________________________________________________

			LOCAL STIFFNESS MATRIX SK (UPPER SYMMETRIC PART)
___________________________________________________________________________*/

   sk [1][1] =  trimempKij (Cl1, All, Cl1);
   sk [1][2] =  trimempKij (Cl1, All, Cl2);
   sk [1][3] =  trimempKij (Cl1, All, Cl3);
   sk [1][4] =  trimempKij (Cl1, All, Cl4);
   sk [1][5] =  trimempKij (Cl1, All, Cl5);
   sk [1][6] =  trimempKij (Cl1, All, Cl6);
   sk [1][7] =  trimempKij (Cl1, All, Cl7);
   sk [1][8] =  trimempKij (Cl1, All, Cl8);
   sk [1][9] =  trimempKij (Cl1, All, Cl9);
   sk [2][2] =  trimempKij (Cl2, All, Cl2);
   sk [2][3] =  trimempKij (Cl2, All, Cl3);
   sk [2][4] =  trimempKij (Cl2, All, Cl4);
   sk [2][5] =  trimempKij (Cl2, All, Cl5);
   sk [2][6] =  trimempKij (Cl2, All, Cl6);
   sk [2][7] =  trimempKij (Cl2, All, Cl7);
   sk [2][8] =  trimempKij (Cl2, All, Cl8);
   sk [2][9] =  trimempKij (Cl2, All, Cl9);
   sk [3][3] =  trimempKij (Cl3, All, Cl3);
   sk [3][4] =  trimempKij (Cl3, All, Cl4);
   sk [3][5] =  trimempKij (Cl3, All, Cl5);
   sk [3][6] =  trimempKij (Cl3, All, Cl6);
   sk [3][7] =  trimempKij (Cl3, All, Cl7);
   sk [3][8] =  trimempKij (Cl3, All, Cl8);
   sk [3][9] =  trimempKij (Cl3, All, Cl9);
   sk [4][4] =  trimempKij (Cl4, All, Cl4);
   sk [4][5] =  trimempKij (Cl4, All, Cl5);
   sk [4][6] =  trimempKij (Cl4, All, Cl6);
   sk [4][7] =  trimempKij (Cl4, All, Cl7);
   sk [4][8] =  trimempKij (Cl4, All, Cl8);
   sk [4][9] =  trimempKij (Cl4, All, Cl9);
   sk [5][5] =  trimempKij (Cl5, All, Cl5);
   sk [5][6] =  trimempKij (Cl5, All, Cl6);
   sk [5][7] =  trimempKij (Cl5, All, Cl7);
   sk [5][8] =  trimempKij (Cl5, All, Cl8);
   sk [5][9] =  trimempKij (Cl5, All, Cl9);
   sk [6][6] =  trimempKij (Cl6, All, Cl6);
   sk [6][7] =  trimempKij (Cl6, All, Cl7);
   sk [6][8] =  trimempKij (Cl6, All, Cl8);
   sk [6][9] =  trimempKij (Cl6, All, Cl9);
   sk [7][7] =  trimempKij (Cl7, All, Cl7);
   sk [7][8] =  trimempKij (Cl7, All, Cl8);
   sk [7][9] =  trimempKij (Cl7, All, Cl9);
   sk [8][8] =  trimempKij (Cl8, All, Cl8);
   sk [8][9] =  trimempKij (Cl8, All, Cl9);
   sk [9][9] =  trimempKij (Cl9, All, Cl9);

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

}/* end of trimempis() */



/*******************************************[ USER LIBRARY : TRIMEMPK.C ]*****
 #                                                                           #
 #  FUNCTION :  trimempCle()                                                 #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of side stiffness matrices C in global coordinates       #
 #      Only one element contribution.                                       #
 #                                                                           #
 #      Enhanced Morley Plate triangular elements.                           #
 #      Element code = 4                                                     #
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


 void trimempCle (double  Cl1[4], double Cl2[4], double Cl3[4], double Cl4[4],
                  double Cl5[4], double Cl6[4], double Ae[4], int n1, int n2,
				  int n3, double E, double nu, double h, int nudim, double *xyz)
/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// double  Cl1[4];  /* vector Cl1 [3 x 1].                                    */
// double  Cl2[4];  /* vector Cl2 [3 x 1].                                    */
// double  Cl3[4];  /* vector Cl3 [3 x 1].                                    */
// double  Cl4[4];  /* vector Cl4 [3 x 1].                                    */
// double  Cl5[4];  /* vector Cl5 [3 x 1].                                    */
// double  Cl6[4];  /* vector Cl6 [3 x 1].                                    */
// double  Ae[4];       /* submatrix diagonal [3 x 3].                        */
// int  n1, n2, n3;     /* global node numbers of element vertexes            */
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

 double  L1, L2, L3; /* length of sides                                     */
 double  L12, L22, L32; /* square length of sides                           */
 double  darea   ;   /* two times the element area = x23*y31 - x31*y23      */
 double  h3;         /* three times plate thickness                         */
 double  b1, b2, b3, c1, c2, c3;           /* geometric constants           */
 double  b12, b22, b32, c12, c22, c32;     /* geometric constants           */
 double  bc12, bc13, bc23, bc1, bc2, bc3;  /* geometric constants           */
 double  cos2, sin2, cos_sin_2;  /* director cosines of side 3              */
 double  Kc[4];        /* vector of global curvatures                       */

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

void trimempKwi(double Kc[4], double b12, double b22, double b32, double c12,
                double c22, double c32, double bc12, double bc13, double bc1,
                double bc2, double bc3, double L22, double L32, double darea);

void trimempKroti(double Kc[4], double b32, double c32, double bc3, double L3,
                  double darea);

void trimempCli (double Cli[4], double Kc[4], double cos2, double sin2,
                 double cos_sin_2, double E, double nu, double h3, double darea);

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

			 ANGLE OF SIDE Li (coincident with L1)
___________________________________________________________________________*/

    cos2 = b12/L12;
    sin2 = c12/L12;
    cos_sin_2 = 2*bc1/L12;

/*__________________________________________________________________________

			 DIAGONAL SUBMATRIX Ae
___________________________________________________________________________*/

        h3 = h*h*h;

        Ae [1] = darea/6*E*h3/12;
		Ae [2] = darea/6*12*(1-nu*nu)/E/h3;
		Ae [3] = darea/6*24*(1+nu)/E/h3;

/*__________________________________________________________________________

									 MATRIX Cl1 with w1
___________________________________________________________________________*/

  trimempKwi(Kc, b12, b22, b32, c12, c22, c32, bc12, bc13, bc1, bc2, bc3,
             L22, L32, darea);

  trimempCli (Cl1, Kc, cos2, sin2, cos_sin_2, E, nu, h3, darea);

/*__________________________________________________________________________

									 MATRIX Cl2 with w2
___________________________________________________________________________*/

  trimempKwi(Kc, b22, b32, b12, c22, c32, c12, bc23, bc12, bc2, bc3, bc1,
             L32, L12, darea);

  trimempCli (Cl2, Kc, cos2, sin2, cos_sin_2, E, nu, h3, darea);

/*__________________________________________________________________________

									 MATRIX Cl3 with w3
___________________________________________________________________________*/

  trimempKwi(Kc, b32, b12, b22, c32, c12, c22, bc13, bc23, bc3, bc1, bc2,
             L12, L22, darea);

  trimempCli (Cl3, Kc, cos2, sin2, cos_sin_2, E, nu, h3, darea);

/*__________________________________________________________________________

									 MATRIX Cl4 with rot4
___________________________________________________________________________*/

  trimempKroti(Kc, b12, c12, bc1, L1, darea);

  trimempCli (Cl4, Kc, cos2, sin2, cos_sin_2, E, nu, h3, darea);

/*__________________________________________________________________________

									 MATRIX Cl5 with rot5
___________________________________________________________________________*/

  trimempKroti(Kc, b22, c22, bc2, L2, darea);

  trimempCli (Cl5, Kc, cos2, sin2, cos_sin_2, E, nu, h3, darea);

/*__________________________________________________________________________

									 MATRIX Cl6 with rot6
___________________________________________________________________________*/

  trimempKroti(Kc, b32, c32, bc3, L3, darea);

  trimempCli (Cl6, Kc, cos2, sin2, cos_sin_2, E, nu, h3, darea);

/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trimempCle() */



/*******************************************[ USER LIBRARY : TRIMEMPK.C  ]*****
 #                                                                           #
 #  FUNCTION :  trimempKwi ()                                                #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of global curvatures vector Kc [3] in global coordinates.#
 #      Only vertex displacements shape functions                            #
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

void trimempKwi(double Kc[4], double bi2, double bj2, double bk2, double ci2,
                double cj2, double ck2, double bcij, double bcik, double bci,
				double bcj, double bck, double Lj2, double Lk2, double darea)

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

}/* end of trimempKwi() */




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

void trimempKroti(double Kc[4], double bi2, double ci2, double bci, double Li, double darea)

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

}/* end of trimempKroti() */




/*******************************************[ USER LIBRARY : TRIMEMPK.C  ]****
 #                                                                           #
 #  FUNCTION :  trimempCli()                                                 #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of stiffness vector Cij [3] in local coordinates.        #
 #      Enhanced Morley Plate triangular elements.                           #
 #      Element code = 4                                                     #
 #                                                                           #
 #***************************************************************************/

void trimempCli (double Cij[4], double Kc[4], double cos2, double sin2,
                 double cos_sin_2, double E, double nu, double h3, double darea)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Cij[4];  /* vector Cij [3]. (side i DOF j)                         */
// double   Kc[4];  /* vector of global curvatures of DOF j                   */
// double  cos2, sin2, cos_sin_2;  /* director cosines of side i              */
// double   E   ;   /* elasticity modulus element                             */
// double   nu  ;   /* Poisson's ratio element                                */
// double   h3  ;   /* cubic thickness of element                             */
// double  darea;   /* two times the element area = x23*y31 - x31*y23         */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

		 ELEMENT STIFFNESS SUBMATRIX MATRIX IN GLOBAL COORDINATES
____________________________________________________________________________*/


     Cij[1] =  darea*E*h3/72*(sin2*Kc[1] + cos2*Kc[2] - cos_sin_2*Kc[3]/2);
     Cij[2] = -darea/6*((cos2+nu*sin2)*Kc[1] + (sin2+nu*cos2)*Kc[2] +
                         cos_sin_2*(1-nu)*Kc[3]/2);
     Cij[3] = -darea/6*(cos_sin_2*(Kc[2]-Kc[1]) + 2*(cos2-sin2)*Kc[3]/2);


/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trimempCli() */


/******************************************[ USER LIBRARY : TRIMEMPK.C  ]*****
 #                                                                           #
 #  FUNCTION :  trimempKij ()                                                #
 #                                                                           #
 #  TYPE     :  double                                                       #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of stiffness coefficient Kij in global coordinates.      #
 #                                                                           #
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

 double trimempKij (double Cli[4], double All[4], double Clj[4])

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Cli[4];    /* submatrix Cli [3 x 1].                               */
// double  Clj[4];    /* submatrix Clj [3 x 1].                               */
// double  All[4];     /* diagonal submatrix All [3 x 3].                      */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 double  CKij ;   /* coefficient Kij                                        */

/*____________________________________________________________________________

		 ELEMENT STIFFNESS SUBMATRIX MATRIX IN GLOBAL COORDINATES
____________________________________________________________________________*/


  CKij = Cli[1]*All[1]*Clj[1] + Cli[2]*All[2]*Clj[2] + Cli[3]*All[3]*Clj[3];

/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

  return (CKij);

}/* end of trimempKij() */

