/**********************[          Trimm2dav.h            ]********************/

#ifndef TR2DMLAV_H
#define TR2DMLAV_H


/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: TR2DMLAV.C                                                 #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #        Routines for computation of local element refinament indicator.    #
 #        Multilevel computation for all elements.                           #
 #        Linear mixed triangular elements. Element code = 3.                #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  tr2dmlav()   : computation of refinament indicator.                #
 #                                                                           #
 #****************************************************************************/


#include <cstdio>
#include <cmath>
#include <stdlib.h>


/*******************************************[ USER LIBRARY : TR2DMLAV.C ]*****
 #                                                                           #
 #  FUNCTION :  tr1_msavg_err()                                              #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Multilevel Computation of element refinament indicator.              #
 #      Error estimated in energy norm between recovered stresses and        #
 #      computed stresses from finite element results.                       #
 #      Recovered stresses computed from averaged midside tractions.         #
 #      Equidistribution of error in energy norm.                            #
 #      Zienckiewicz's strategy mean error.                                  #
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

 void tr1_msavg_err (int curlev , int *numel     , int *FRnumel,
                      int **velm      , int xndof      , double **xyz,
                      int **elm_sides , int **elm_nodes, int **elm_vlev,
                      int **parent_elm, double **uold  , int *mat,
                      double *prop    , int xnumpr     ,
                      double  **Fn     , double  **Fs    ,
                      double **error  , double *Uerr   , double *Unew,
                      double *Ufea);

/*____________________________________________________________________________

							 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int   curlev      ; /* current level in analysis.                          */
// int  *numel       ; /* number of elements in each level.                   */
// int  *FRnumel     ; /* number of refined elements (all in group Rnumel).   */
// int  **velm       ; /* new element ordering after refinament [numel[level]]*/
// int    xndof      ; /* maximum number of nodal degrees of freedom.         */
// double **xyz      ; /* nodal coordinates     [nodes[curlev] * nudim]       */
// int  **elm_sides  ; /* sides incidences of each element [numel * xnes]     */
// int  **elm_nodes  ; /* nodes of each element    [numel[curlev]*xnen]       */
// int  **elm_vlev   ; /* level of element vertexs [numel * xnev]             */
// int  **parent_elm ; /* parent element and location [3*numel[curlev]]       */
// double **uold     ; /* vector of previous nodal displacements [xndof*nodes]*/
// int    *mat       ; /* element material type vector [numel].               */
// double *prop      ; /* material properties vector [numat * numpr].         */
// int     xnumpr    ; /* maximum number of material properties.              */
// double  **Fn       ; /* side normal tractions [sides[curlev]]               */
// double  **Fs       ; /* side tangential tractons [sides[curlev]]            */
// double **error    ; /* element error energy norm of each level             */
// double *Uerr   ;  /* energy of error                                       */
// double *Unew   ;  /* energy of new displacements for all elements          */
// double *Ufea   ;   /* energy of finite element approximation               */

/*****************************************************************************/


#endif
