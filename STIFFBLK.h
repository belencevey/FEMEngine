/**********************[          STIFFBLK.h            ]********************/

#ifndef STIFFBLK_H
#define STIFFBLK_H

/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: STIFFBLK.C                                                 #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #  void  stiffblk()   : computation of global stiffness matrix.             #
 #                                                                           #
 #****************************************************************************/

#include <cstdio>
#include <stdlib.h>
#include"TRIM2DK.h"

/*******************************************[ USER LIBRARY : STIFFBLK.C ]*****
 #                                                                           #
 #  FUNCTION :  stiffblk ()                                                  #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of global stiffness matrix for blocks of elements.       #
 #      Computation of coefficients Gij for right sided triangles.           #
 #                                                                           #
 #      Element codes for PLANE elements                                     #
 #                                                                           #
 #         1   : LINEAR TRIANGLE VERTEX NODES - PLANE STRESS                 #
 #         2   : QUADRATIC TRIANGLE RIGHT SIDES - PLANE STRESS               #
 #                                                                           #
 #                                                                           #
 #      Only the upper triangular part of the element matrix is computed     #
 #      and stored in vector form by rows. Then is assembled into the        #
 #      global stiffness matrix.                                             #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED LIBRARY FUNCTIONS USED ]--------------------------------#
 #                                                                           #
 #                                                                           #
 #  TRIM2DK.C :                                                              #
 #    void  trim2dk() : computation of stiffness matrix for linear triangular#
 #                      elements for plane states.                           #
 #                                                                           #
 #                                                                           #
 #***************************************************************************/


	void stiffblk (int ecode, int numel, int sides, int bsides, int xndof, int nudim, int *id, int *elm_sides,
                   int *side_nodes, int *side_bc, int xnsn, double *xyz, int *side_mat, int *mat,
                   double *prop, int xnumpr, double *sg, int *maxa, int *lm, double *sv, double *Gij,
                   double *Kij, double *diag, double *fp, double *fm);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// int     ecode;   /* element code number.                                   */
// int     numel;   /* number of elements in this block.                      */
// int     sides;   /* number of sides.                                       */
// int    bsides;   /* number of boundary sides.                              */
// int     xndof;   /* maximum number of nodal degrees of freedom.            */
// int     nudim;   /* number of spatial dimensions.                          */
// int    *id   ;   /* identification vector [numnp * ndof].                  */
// int *elm_sides ; /* sides of each element    [numel * xnes]                */
// int *side_nodes; /* nodes of each side       [sides * xnsn]                */
// int  *side_bc   ; /* side boundary conditions [bsides * xndof]             */
// int   xnsn   ;   /* maximum number of side nodes.                          */
// double *xyz  ;   /* nodal coordinates [numnp * nudim].                     */
// int    *side_mat; /* side material index of adjacent elements [2*sides[0]] */
// int    *mat  ;   /* element material type vector [numel].                  */
// double *prop ;   /* material properties vector [numat * numpr].            */
// int    xnumpr;   /* maximum number of material properties.                 */
// double  *sg  ;   /* global stiffness matrix stored by skyline [nwk]        */
// int     *maxa;   /* pointer of diagonal elements for skyline [neq+1]       */
// int     *lm  ;   /* identification vector of element dof [nrw]             */
// double  *sv  ;   /* upper triangular element matrix [nrw *(nrw +1)/2]      */
// double *Gij  ;   /* element stiffness coefficients [10*numel]              */
// double *Kij  ;   /* element stiffness coefficients [21*numel]              */
// double *diag ;   /* stiffness diagonal coefficients                        */
// double *fp   ;   /* concentrated nodal loads coarse mesh only.             */
// double *fm   ;   /* concentrated mixed nodal loads coarse mesh only.       */

/******************************************************************************/

#endif
