/**********************[          Prtplt.h            ]********************/

#ifndef PRTPLT_H
#define PRTPLT_H

/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: PRTPLT.C                                                   #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #  void  prtplt()     : write geometry and stresses in output file.         #
 #  void  stress()     : write stresses in output file.                      #
 #                                                                           #
 #****************************************************************************/

#include <cstdio>
#include <stdlib.h>
#include <string>
#include<fstream>
#include<iostream>
#include<iomanip>

#include"TRIM2DS.h"
#include"Trimm2ds.h"


/*******************************************[ USER LIBRARY : PRTPLT.C   ]*****
 #                                                                           #
 #  FUNCTION :  prtplt()                                                     #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Write geometry data and stresses in output file for plot.            #
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
 #    void  trim2ds() : computation of stresses for linear triangular        #
 #                      elements for plane states.                           #
 #                                                                           #
 #                                                                           #
 #***************************************************************************/

   void prtplt (std::string *forname, int curlev, int curmesh, int ecode, int *numnp,
               int *acnodes, int *acnumel, int *numel, int *frnumel, int xndof,
			   int nudim, int **elm_sides, int **elm_nodes, int **elm_vlev,
			   int **parent_elm, int **velm, double **xyz, int *mat, double *prop,
			   int xnumpr, double **uold, double **Fn, double **Fs, double **Sx,
			   double **Sy, double **Sxy, double **error, double errm,
			   int ESTIMATOR, int SOLVER);

/*
  void prtplt (char *forname, int curlev, int curmesh, int ecode, int *numnp,
               int *acnodes, int *acnumel, int *numel, int *frnumel, int xndof,
			   int nudim, int **elm_sides, int **elm_nodes, int **elm_vlev,
			   int **parent_elm, int **velm, double **xyz, int *mat, double *prop,
			   int xnumpr, double **uold, double **Fn, double **Fs, double **Sx,
			   double **Sy, double **Sxy, double **error, double errm,
			   int ESTIMATOR, int SOLVER);*/


/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// char *forname   ; /* complete output file root name                        */
// int  curlev     ; /* fine level in analysis.                               */
// int  curmesh    ; /* current mesh in analysis.                             */
// int  ecode      ; /* element code number.                                  */
// int *numnp      ; /* number of nodal points                                */
// int *acnodes    ; /* number of accumulated nodes in previous levels        */
// int *acnumel    ; /* number of accumulated elements in previous levels     */
// int *numel      ; /* number of elements in each level.                     */
// int *frnumel    ; /* number of refined elements in each level.             */
// int  xndof      ; /* maximum number of nodal degrees of freedom.           */
// int  nudim      ; /* number of spatial dimensions.                         */
// int **elm_sides ; /* sides incidences of each element [numel * xnes]       */
// int **elm_nodes ; /* nodes of each element    [numel[curlev]*xnen]         */
// int **elm_vlev  ; /* levels of element vertex nodes  [xnev*numel[curlev]]  */
// int **parent_elm; /* parent element and location [3*numel[curlev]]         */
// int **velm      ; /* new element ordering after refinament[numel[curlev]]  */
// double **xyz    ; /* nodal coordinates [numnp * nudim].                    */
// int    *mat     ; /* element material type vector [numel].                 */
// double *prop    ; /* material properties vector [numat * numpr].           */
// int    xnumpr   ; /* maximum number of material properties.                */
// double **uold   ; /* vector of nodal displacements [xndof*numnp]           */
// double  **Fn     ; /* midside normal traction       [sides[level]]          */
// double  **Fs     ; /* midside tangential traction   [sides[level]]          */
// double  **Sx     ; /* nodal stress Sx                 [nodes[level]]        */
// double  **Sy     ; /* nodal stress Sy                 [nodes[level]]        */
// double  **Sxy    ; /* nodal stress Sxy                [nodes[level]]        */
// double **error  ; /* element error                 [numel[level]]          */
// double errm     ; /* mean error in energy norm of fine mesh                */
// int ESTIMATOR; /* estimator type 1: nodal avg.  2: midside stress  3: ZZ   */
// int  SOLVER  ;  /* type of iterative solver                                */
//                 /*    SOLVER = 1   hierarchical PCG                        */
//                 /*    SOLVER = 2   diagonal PCG                            */
//                 /*    SOLVER = 3   hierarchical MULTIGRID                  */
//                 /*    SOLVER = 4   diagonal MULTIGRID                      */

/****************************************************************************/







/*******************************************[ USER LIBRARY : PRTPLT.C   ]*****
 #                                                                           #
 #  FUNCTION :  stress()                                                     #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of local element stresses without interpolation.         #
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
 #    void  trim2ds() : computation of stresses for linear triangular        #
 #                      elements for plane states.                           #
 #                                                                           #
 #                                                                           #
 #***************************************************************************/


  void stress (FILE *fpo, FILE *fGiD, int ecode, int *acnodes, int elm1,
               int elm2, int xndof, int nudim, int *elm_sides, int *elm_nodes,
			   int *elm_vlev, int *parent_elm, int *velm, double **xyz,
			   int *mat,  double *prop, int xnumpr, double **uold,
			   double *Fn, double *Fs, double **Sx, double **Sy,
			   double **Sxy, double *error, double errm, int ESTIMATOR, int level);


/*____________________________________________________________________________

								PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE   *fpo    ; /* pointer to output file                                 */
// FILE   *fGiD   ; /* pointer to GiD output file                             */
// int     ecode  ; /* element code number.                                   */
// int *acnodes   ; /* number of accumulated nodes in previous levels        */
// int  elm1, elm2; /* first and last elements of this level.                 */
// int  xndof     ; /* maximum number of nodal degrees of freedom.            */
// int  nudim     ; /* number of spatial dimensions.                          */
// int *elm_sides ; /* sides incidences of each element [numel * xnes]        */
// int *elm_nodes ; /* nodes of each element    [numel[curlev]*2*xnen]        */
// int *elm_vlev  ; /* levels of element vertex nodes  [xnev*numel[curlev]]   */
// int *parent_elm; /* parent element and location [3*numel[curlev]]          */
// int *velm      ; /* new element ordering after refinament[numel[curlev]]   */
// double **xyz   ; /* nodal coordinates [numnp * nudim].                     */
// int    *mat    ; /* element material type vector [numel].                  */
// double *prop   ; /* material properties vector [numat * numpr].            */
// int    xnumpr  ; /* maximum number of material properties.                 */
// double **uold  ; /* vector of nodal displacements [xndof*numnp]            */
// double  *Fn     ; /* midside normal traction       [sides[level]]           */
// double  *Fs     ; /* midside tangential traction   [sides[level]]           */
// double  **Sx    ; /* nodal stress Sx                    [nodes[level]]      */
// double  **Sy    ; /* nodal stress Sy                    [nodes[level]]      */
// double  **Sxy   ; /* nodal stress Sxy                   [nodes[level]]      */
// double *error  ; /* element error of this level    [numel[level]]          */
// double errm    ; /* mean error in energy norm of fine mesh                 */
// int ESTIMATOR  ; /* estimator type 1:nodal avg. 2:midside stress           */
// int  level     ; /* current level in analysis.                             */

/*****************************************************************************/




/*******************************************[ USER LIBRARY : PRTPLT.C  ]******
 #                                                                           #
 #  FUNCTION :  elmdesc ()                                                   #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Description of element incidences for Gid output.                    #
 #                                                                           #
 #***************************************************************************/


  void elmdesc (FILE *fGiDmesh, int *acnodes, int elm1, int elm2, int *elm_nodes,
                int *elm_vlev, int *parent_elm, int *velm, int *mat);


/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE *fGiDmesh ; /* pointer to GiD mesh file                               */
// int *acnodes   ; /* number of accumulated nodes in previous levels         */
// int  elm1, elm2; /* first and last elements of this level.                 */
// int *elm_nodes ; /* nodes of each element    [numel[curlev]*xnen]          */
// int *elm_vlev  ; /* level of element vertexs [numel[curlev]*xnev]          */
// int *parent_elm; /* parent element and location [3*numel[curlev]]          */
// int *velm      ; /* new element ordering after refinament[numel[curlev]]   */
// int    *mat    ; /* element material type vector [numel].                  */

/****************************************************************************/



#endif
