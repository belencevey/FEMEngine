/**********************[          Trimm2dmf.h            ]********************/

#ifndef TRIMM2DMF_H
#define TRIMM2DMF_H


/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: TR2DMLMF.C                                                 #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #        Routines for computation of midside interpolated forces.           #
 #        All elements of current mesh for all levels.                       #
 #        Linear triangular elements. Element code = 1.                      #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  trimm2dmf() : Computation of midside forces linear mixed triangle. #
 #                                                                           #
 #****************************************************************************/


#include <cstdio>
#include <cmath>
#include <stdlib.h>
#include <memory.h>
//#include <mem.h>


/*******************************************[ USER LIBRARY : TRIMM2DMF.C ]****
 #                                                                           #
 #  FUNCTION :  t1m2dmf()                                                    #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of side stresses for linear mixed triangles.             #
 #      Linear variation for the stresses along the side is assumed.         #
 #      Average value is assigned for sides belonging to two elements.       #
 #      For boundary sides, external tractions are impossed.                 #
 #                                                                           #
 #      Linear triangular elements.                                          #
 #      Element code = 1.                                                    #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  MATH.H  :                                                                #
 #    double sqrt() : calculates square root.                                #
 #                                                                           #
 #  STDLIB.H:                                                                #
 #    int abs()     : macro that return the absolute value of integer.       #
 #                                                                           #
 #                                                                           #
 #***************************************************************************/


  void t1m2dmf (int sides, int bsides, int ndof,
                 int nudim, int *side_nodes, int *side_bc,
                 int xnsn, double *xyz, int *side_mat,
                 int *mat, double *prop, int numpr,
                 double *u, double *Fn, double *Fs,
                 double *bs_Fn, double *bs_Fs);

/*____________________________________________________________________________

							 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int     sides;   /* number of sides.                                       */
// int    bsides;   /* number of boundary sides.                              */
// int     ndof ;   /* maximum number of nodal degrees of freedom.            */
// int     nudim;   /* number of spatial dimensions.                          */
// int *side_nodes; /* nodes of each side       [sides * xnsn]                */
// int    *side_bc; /* side boundary conditions [Bsides[0] * xndof]           */
// int   xnsn   ;   /* maximum number of side nodes.                          */
// double *xyz  ;   /* nodal coordinates [numnp * nudim].                     */
// int *side_mat;   /* side material index of adjacent elements [2*sides[0]]  */
// int    *mat  ;   /* element material type vector [numel].                  */
// double *prop ;   /* material properties vector [numat * numpr].            */
// int     numpr;   /* maximum number of material properties.                 */
// double *u    ;   /* expanded vector of nodal displacements [xndof*numnp]   */
// double  *Fn   ;   /* side normal tractions [sides[curlev]]                  */
// double  *Fs   ;   /* side tangential tractons [sides[curlev]]               */
// double *bs_Fn;   /* midside normal boundary tractions current level        */
// double *bs_Fs;   /* midside tangential boundary tractions current level    */

/*****************************************************************************/



/*******************************************[ USER LIBRARY : TRIMM2DMF.C ]****
 #                                                                           #
 #  FUNCTION :  t1m2D_stress_bs()                                            #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of side stresses for linear mixed triangles.             #
 #      Linear variation for the stresses along the side is assumed.         #
 #      Average value is assigned for sides belonging to two elements.       #
 #      ONLY FOR BOUNDARY SIDES.                                             #
 #                                                                           #
 #      Linear mixed triangular elements.                                    #
 #      Element code = 3.                                                    #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  MATH.H  :                                                                #
 #    double sqrt() : calculates square root.                                #
 #                                                                           #
 #  STDLIB.H:                                                                #
 #    int abs()     : macro that return the absolute value of integer.       #
 #                                                                           #
 #                                                                           #
 #***************************************************************************/

   void t1m2D_stress_bs (int bsides, int ndof, int nudim, int *side_nodes,
						  int *side_bc, int xnsn, double *xyz, int *side_mat,
						  int *mat, double *prop, int numpr, double *u, double *Fn,
						  double *Fs, double *bs_Fn, double *bs_Fs);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// int    bsides;   /* number of boundary sides.                              */
// int     ndof ;   /* maximum number of nodal degrees of freedom.            */
// int     nudim;   /* number of spatial dimensions.                          */
// int *side_nodes; /* nodes of each side       [sides * xnsn]                */
// int    *side_bc; /* side boundary conditions [Bsides[0] * xndof]           */
// int   xnsn   ;   /* maximum number of side nodes.                          */
// double *xyz  ;   /* nodal coordinates [numnp * nudim].                     */
// int *side_mat;   /* side material index of adjacent elements [2*sides[0]]  */
// int    *mat  ;   /* element material type vector [numel].                  */
// double *prop ;   /* material properties vector [numat * numpr].            */
// int     numpr;   /* maximum number of material properties.                 */
// double *u    ;   /* expanded vector of nodal displacements [xndof*numnp]   */
// double  *Fn   ;   /* side normal tractions [sides[curlev]]                  */
// double  *Fs   ;   /* side tangential tractons [sides[curlev]]               */
// double *bs_Fn;   /* midside normal boundary tractions current level        */
// double *bs_Fs;   /* midside tangential boundary tractions current level    */

/*****************************************************************************/






/*******************************************[ USER LIBRARY : TRIMM2DMF.C ]****
 #                                                                           #
 #  FUNCTION :  t1m2D_stress_is()                                            #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of side stresses for linear mixed triangles.             #
 #      Linear variation for the stresses along the side is assumed.         #
 #      Average value is assigned for sides belonging to two elements.       #
 #      ONLY FOR INTERIROR SIDES.                                            #
 #                                                                           #
 #      Linear mixed triangular elements.                                    #
 #      Element code = 3.                                                    #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  MATH.H  :                                                                #
 #    double sqrt() : calculates square root.                                #
 #                                                                           #
 #  STDLIB.H:                                                                #
 #    int abs()     : macro that return the absolute value of integer.       #
 #                                                                           #
 #                                                                           #
 #***************************************************************************/

  void t1m2D_stress_is (int sides, int bsides, int ndof, int nudim, int *side_nodes,
                        int xnsn, double *xyz, int *side_mat, int *mat, double *prop,
						int numpr, double *u, double *Fn, double *Fs);


/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// int     sides;   /* number of sides.                                       */
// int    bsides;   /* number of boundary sides.                              */
// int     ndof ;   /* maximum number of nodal degrees of freedom.            */
// int     nudim;   /* number of spatial dimensions.                          */
// int *side_nodes; /* nodes of each side       [sides * xnsn]                */
// int   xnsn   ;   /* maximum number of side nodes.                          */
// double *xyz  ;   /* nodal coordinates [numnp * nudim].                     */
// int *side_mat;   /* side material index of adjacent elements [2*sides[0]]  */
// int    *mat  ;   /* element material type vector [numel].                  */
// double *prop ;   /* material properties vector [numat * numpr].            */
// int     numpr;   /* maximum number of material properties.                 */
// double *u    ;   /* expanded vector of nodal displacements [xndof*numnp]   */
// double  *Fn   ;   /* side normal tractions [sides[curlev]]                  */
// double  *Fs   ;   /* side tangential tractons [sides[curlev]]               */

/*****************************************************************************/




/*******************************************[ USER LIBRARY : TRIMM2DMF.C ]****
 #                                                                           #
 #  FUNCTION :  t1m2dSideStress()                                            #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of side stresses for linear mixed triangles.             #
 #      Planar Strain Stress.                                                #
 #      Only one element contribution.                                       #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  MATH.H  :                                                                #
 #    double sqrt() : calculates square root.                                #
 #                                                                           #
 #  STDLIB.H:                                                                #
 #    int abs()     : macro that return the absolute value of integer.       #
 #                                                                           #
 #                                                                           #
 #***************************************************************************/


  void t1m2dSideStress (double *Sn, double *Tsn, double *darea, int n1, int n2, int n3,
                        double E, double nu, int nudim, double *xyz, int xndof, double *u);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// double *Sn,*Tsn; /* midside tractions                                      */
// double  *darea;  /* two times the element area = x23*y31 - x31*y23 	    */
// int  n1, n2, n3; /* node numbers of element                                */
// double   E   ;   /* elasticity modulus element                             */
// double   nu  ;   /* Poisson's ratio element                                */
// int     nudim;   /* number of spatial dimensions.                          */
// double *xyz  ;   /* nodal coordinates [numnp * nudim].                     */
// int    xndof ;   /* maximum number of nodal degrees of freedom.            */
// double *u    ;   /* expanded vector of nodal displacements [xndof*numnp]   */

/******************************************************************************/



#endif
