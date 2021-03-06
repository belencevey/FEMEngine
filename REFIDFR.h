/**********************[          REFINFR.h            ]********************/

#ifndef REFINFR_H
#define REFINFR_H


/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: REFIDFR.C                                                  #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #  void  refidfr()    : verification of element refinament                  #
 #                                                                           #
 #****************************************************************************/

#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include "Trimm2dmf.h"
#include "TR2DMLMF.h"
/*******************************************[ USER LIBRARY : REFID.C    ]*****
 #                                                                           #
 #  FUNCTION :  refidfr ()                                                   #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of local element refinament indicators.                  #
 #      Only for triangles.                                                  #
 #                                                                           #
 #                                                                           #
 #      Element codes for PLANE elements                                     #
 #                                                                           #
 #         1   : LINEAR TRIANGLE VERTEX NODES - PLANE STRESS                 #
 #         2   : QUADRATIC TRIANGLE RIGHT SIDES - PLANE STRESS               #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED LIBRARY FUNCTIONS USED ]--------------------------------#
 #                                                                           #
 #                                                                           #
 #  TRIM2DS.C :                                                              #
 #    void  trim2dms(): computation of interpolated midside stresses.        #
 #                      Linear triangular elements for plane states.         #
 #  TRIM2DS.C :                                                              #
 #    void  trim2dzz(): computation of refinament indicators for linear      #
 #                      triangular elements for plane states.                #
 #                      Zienckiewickz-Zhu procedure.                         #
 #                                                                           #
 #                                                                           #
 #***************************************************************************/

 void refidfr (int ecode      , int curlev,
                       int *acnumel, int *numel, int *nodes,
                       int *rnumel,
                       int *FRnumel, int **velm,  int xndof,
                       int nudim, int xnsn, double **xyz,
                       int **elm_sides , int **elm_nodes, int **elm_vlev,
                       int **parent_elm, int *sides    , int *Bsides,
                       int *BFRsides, int *FRsides, int *PRsides,
                       int **vside  ,  int **eside, int **bside,
                       int **tside, int **ssside,
                       int **side_nodes, int **side_vlev,
                       int **parent_side, int *side_bc, double **uold,
                       int *mat       , double *prop  , int xnumpr,
                       int *side_mat,
                       int **relm    , double **error,
                       double  **Fn   , double  **Fs,
                       double **Sx,  double **Sy,double **Sxy,
                       int *nse, int **idnode,
                       double **bs_Fn , double **bs_Fs,
                       double **a11, double **a12, double **a13,
                       double **a22, double **a23, double **a33,
                       double **b1x  , double **b2x  , double **b3x,
                       double **b1y  , double **b2y  , double **b3y,
                       double **b1xy , double **b2xy , double **b3xy,
                       int *refine,  int *levref, int ESTIMATOR,
                       long *TRnumel,
                       double *Uerr, double *Unew, double *Ufea,
                       double *errm , double **fi,
                       double **bn_nlength, double **bn_slength,
                       double **bn_qn, double **bn_qs,
                       double **bs_cs, double **bs_sn);

/*____________________________________________________________________________

				          PARAMETERS DECLARATION
____________________________________________________________________________*/

// int    ecode      ; /* element code number.                                */
// int   curlev      ; /* current level in analysis.                          */
// int   *acnumel    ; /* number of accumulated nodes in each level           */
// int   *numel      ; /* number of elements of each level.                   */
// int   *nodes      ; /* number of nodal points of each level.               */
// int   *Rnumel     ; /* number of refinable elements of each level.         */
// int   *FRnumel    ; /* number of refined elements (all in group rnumel).   */
// int  **velm       ; /* new element ordering after refinament [numel[level]]*/
// int    xndof      ; /* maximum number of nodal degrees of freedom.         */
// int   nudim       ; /* number of spatial dimensions.                       */
// int   xnsn        ; /* maximum number of side nodes.                       */
// double **xyz      ; /* nodal coordinates     [numnp[curlev] * nudim]       */
// int  **elm_sides  ; /* sides of each element of each level.                */
// int  **elm_nodes  ; /* nodes of each element of each level.                */
// int  **elm_vlev   ; /* vertex levels of each element of each level.        */
// int  **parent_elm ; /* parent element and location [3*numel[curlev]]       */
// int   *sides      ; /* number of sides of each level.                      */
// int   *Bsides     ; /* number of boundary sides of each level.             */
// int   *BFRsides   ; /* number of boundary sides not refined in each level. */
// int   *FRsides    ; /* number of refined sides in each level.              */
// int   *PRsides    ; /* number of transition sides in each level.           */
// int  **vside      ; /* new side ordering after refinament [sides[level]]   */
// int  **eside      ; /* inverse side ordering                               */
// int  **bside      ; /* list of boundary sides of each level.               */
// int  **tside      ; /* list of transition sides of each side.              */
// int  **ssside;
// int  **side_nodes ; /* nodes of each side       [sides[curlev] * xnsn]     */
// int  **side_vlev  ; /* levels of side vertex nodes  [xnev*numel[curlev]]   */
// int  **parent_side; /* parent side and location of each side of each level.*/
// int    *side_bc   ; /* side boundary conditions [bsides[0] * xndof]        */
// double **uold     ; /* vector of nodal displacements in previous levels    */
// int    *mat       ; /* element material type vector [numel].               */
// double *prop      ; /* material properties vector [numat * numpr].         */
// int     xnumpr    ; /* maximum number of material properties.              */
// int    *side_mat  ; /* side material index of adjacent elements [2*sides[0]] */
// int    **relm     ; /* element refinament indicator   [numel[level]]       */
// double  **error    ; /* element error energy norm      [rnumel[level]]      */
// double  **Fn       ; /* midside normal tractions of each level.             */
// double  **Fs       ; /* midside tangential tractions of each level.         */
// double  **Sx  ; /* nodal stress Sx                           [nodes[level]] */
// double  **Sy  ; /* nodal stress Sy                           [nodes[level]] */
// double  **Sxy ; /* nodal stress Sxy                          [nodes[level]] */
// int    *nse    ;  /* number of elements at each node of coarse mesh        */
// int    **idnode     ; /* identification of nodes on boundary sides.          */
// double **bs_Fn    ; /* midside normal boundary tractions of each level.    */
// double **bs_Fs    ; /* midside tangential boundary tractions of each level.*/
// double **a11      ; /* Z-Z coefficient matrix                              */
// double **a12      ; /* Z-Z coefficient matrix                              */
// double **a13      ; /* Z-Z coefficient matrix                              */
// double **a22      ; /* Z-Z coefficient matrix                              */
// double **a23      ; /* Z-Z coefficient matrix                              */
// double **a33      ; /* Z-Z coefficient matrix                              */
// double **b1x ; /* Z-Z right hand side component                            */
// double **b2x ; /* Z-Z right hand side component                            */
// double **b3x ; /* Z-Z right hand side component                            */
// double **b1y ; /* Z-Z right hand side component                            */
// double **b2y ; /* Z-Z right hand side component                            */
// double **b3y ; /* Z-Z right hand side component                            */
// double **b1xy; /* Z-Z right hand side component                            */
// double **b2xy; /* Z-Z right hand side component                            */
// double **b3xy; /* Z-Z right hand side component                            */
// int     *refine   ; /* refinement indicator.                               */
// int  *levref      ; /* level refinement indicator.                         */
// int ESTIMATOR     ; /* estimator type 1:nodal avg. 2:midside stress        */
// double *Uerr   ;    /* energy of error                                     */
// long   *TRnumel;    /* total number of element to be refined               */
// double *Unew   ;    /* energy of new displacements for all elements        */
// double *Ufea   ;    /* energy of finite element approximation              */
// double *errm   ;    /* mean error in energy norm                           */
// double **fi;
// double **bn_nlength; /* boundary node summ of adjacent lengths normal      */
// double **bn_slength; /* boundary node summ of adjacent lengths tangent     */
// double **bn_qn    ; /* boundary node normal reaction                       */
// double **bn_qs    ; /* boundary node tangent reaction                      */
// double **bs_cs    ; /* boundary side cosine of local system angle          */
// double **bs_sn    ; /* boundary side sine   of local system angle          */


/******************************************************************************/



#endif
