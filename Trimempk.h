/**********************[          Trimempk.h            ]********************/

#ifndef TRIMEMPK_H
#define TRIMEMPK_H

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


#include <cmath>
#include"SKYLINE.h"
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
                double *sv, double *fp, double *fm);

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

/*****************************************************************************/




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


void trimempbs (int bsides , int nudim  , int *id, int *side_nodes,
                 int *side_bc, int xnsn,
				 double *xyz, int *side_mat, int *mat, double *prop,
                 int numpr, double *sg, int *maxa, int *lm,
                 double *sv, double *fp, double *fm);

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

/******************************************************************************/





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
                 double *sv);

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

/******************************************************************************/




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
				  int n3, double E, double nu, double h, int nudim, double *xyz);

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

/*****************************************************************************/







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
				double bcj, double bck, double Lj2, double Lk2, double darea);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Kc[4];        /* vector of global curvatures                       */
// double  bi2, bj2, bk2, ci2, cj2, ck2;     /* geometric constants           */
// double  bcij, bcik, bci, bcj, bck;        /* geometric constants           */
// double  Lj2, Lk2;   /* square length of sides                              */
// double  darea   ;   /* two times the element area = x23*y31 - x31*y23      */

/****************************************************************************/



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

void trimempKroti(double Kc[4], double bi2, double ci2, double bci, double Li, double darea);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Kc[4];         /* vector of global curvatures                      */
// double  bi2, ci2, bci; /* geometric constants                              */
// double  Li;            /* length of side i                                 */
// double  darea   ;   /* two times the element area = x23*y31 - x31*y23      */

/**************************************************************************/




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
                 double cos_sin_2, double E, double nu, double h3, double darea);

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

/****************************************************************************/




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

 double trimempKij (double Cli[4], double All[4], double Clj[4]);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Cli[4];    /* submatrix Cli [3 x 1].                               */
// double  Clj[4];    /* submatrix Clj [3 x 1].                               */
// double  All[4];     /* diagonal submatrix All [3 x 3].                      */

/*******************************************************************************/


#endif
