/**********************[          Trimmpk.h            ]********************/
#ifndef TRIMMPK_H
#define TRIMMPK_H

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


#include <cmath>
#include"SKYLINE.h"

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


 void trimmpbs (int bsides, int nudim, int *id, int *side_nodes, int *side_bc,
                int xnsn, double *xyz, int *side_mat, int *mat, double *prop,
                int numpr, double *sg, int *maxa, int *lm, double *sv);

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

/*****************************************************************************/





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


 void trimmpis (int sides, int bsides, int nudim, int *id, int *side_nodes,
				 int *side_bc, int xnsn, double *xyz, int *side_mat, int *mat,
				 double *prop, int numpr, double *sg, int *maxa, int *lm, double *sv);

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

/***************************************************************************/



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


 void trimmpKce (double Kw1e[4], double Kw2e[4], double Kw3e[4], double Krot4e[4],
                 double Krot5e[4], double Krot6e[4], int n1, int n2, int n3,
                 int nudim, double *xyz);

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

/******************************************************************************/


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

void trimmpKwi(double Kc[4], double bi2, double bj2, double bk2, double ci2,
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

void trimmpKroti(double Kc[4], double bi2, double ci2, double bci, double Li, double darea);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Kc[4];         /* vector of global curvatures                      */
// double  bi2, ci2, bci; /* geometric constants                              */
// double  Li;            /* length of side i                                 */
// double  darea   ;   /* two times the element area = x23*y31 - x31*y23      */

/*****************************************************************************/





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
 double trimmpKij (double Kci[4], double E, double nu, double h, double Kcj[4]);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Kci[4];    /* submatrix Kci [3 x 1].                               */
// double  Kcj[4];    /* submatrix Kcj [3 x 1].                               */
// double  E    ;   /* elasticity modulus                                     */
// double  nu   ;   /* Poisson's ratio                                        */
// double  h    ;   /* thickness                                              */
/****************************************************************************/






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


 double darea (int n1, int n2, int n3, int nudim, double *xyz);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// int  n1, n2, n3;     /* global node numbers of element vertexes            */
// int     nudim;   /* number of spatial dimensions.                          */
// double *xyz  ;   /* nodal coordinates [numnp * nudim].                     */

/********************************************************************************/

#endif
