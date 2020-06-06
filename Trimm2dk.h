/**********************[          Trimm2dk.h            ]********************/

#ifndef TRIMM2DK_H
#define TRIMM2DK_H


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


#include <cmath>
#include"SKYLINE.h"
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
                double *sv);

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

/****************************************************************************/






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
				 double *prop, int numpr, double *sg, int *maxa, int *lm, double *sv);

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

/**************************************************************************/





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
                 double *prop, int numpr, double *sg, int *maxa,  int *lm, double *sv);



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

/****************************************************************************/





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
                int n1, int n2, int n3, double E, double nu, double h, int nudim, double *xyz);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// double  Cl1[4][3];   /* submatrix Cl1 [3 x 2].                             */
// double  Cl2[4][3];   /* submatrix Cl2 [3 x 2].                             */
// double  Cl3[4][3];   /* submatrix Cl3 [3 x 2].                             */
// double  Dll[4];   /* submatrix diagonal [3 x 3].                           */
// int  n1, int n2, int n3;     /* node numbers of element                    */
// double   E   ;   /* elasticity modulus element                             */
// double   nu  ;   /* Poisson's ratio element                                */
// double   h   ;   /* thickness of element                                   */
// int     nudim;   /* number of spatial dimensions.                          */
// double *xyz  ;   /* nodal coordinates [numnp * nudim].                     */

/*****************************************************************************/


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

 void trimmDll (double Dll[4], double darea, double  E, double nu, double h);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Dll[4];  /* submatrix Cij [3 x 2].                                 */
// double   darea;  /* two times the element area 		                    */
// double       E;  /* elasticity modulus                                     */
// double      nu;  /* Poisson's ratio                                        */
// double   h   ;   /* thickness of element                                   */
/*****************************************************************************/


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

 void trimmCij (double Cij[4][3], double Lj, double E, double nu, double h, double ci, double si, double cj, double sj);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Cij[4][3];   /* submatrix Cij [3 x 2].                             */
// double     Lj;   /* length of side j                                       */
// double      E;   /* elasticity modulus                                     */
// double     nu;   /* Poisson's ratio                                        */
// double   h   ;   /* thickness of element                                   */
// double     ci;   /* cosine of angle alpha i of side i                      */
// double     si;   /* sine of angle alpha i of side i                        */
// double     cj;   /* cosine of angle alpha j of side j                      */
// double     sj;   /* sine of angle alpha j of side j                        */

/*****************************************************************************/




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

void trimmCii (double Cii[4][3], double Li, double ci, double si, double h);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Cii[4][3];   /* submatrix Cii [3 x 2].                             */
// double     Li;   /* length of side i                                       */
// double     ci;   /* cosine of angle alpha i of side i                      */
// double     si;   /* sine of angle alpha i of side i                        */
// double   h   ;   /* thickness of element                                   */

/*******************************************************************************/




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

 void trimmKij (double Kij[3][3], double Cli[4][3], double Clj[4][3], double* Dll);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/


// double  Kij[3][3];   /* submatrix Kij [2 x 2].                               */
// double  Cli[4][3];   /* submatrix Cli [3 x 2].                               */
// double  Clj[4][3];   /* submatrix Clj [3 x 2].                               */
// double   *Dll;   /* diadonal submatrix Dll [3 x 3].                          */

/*********************************************************************************/


#endif
