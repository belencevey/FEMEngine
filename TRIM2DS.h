/**********************[          TRIM2DS.h            ]********************/

#ifndef TRIM2DS_H
#define TRIM2DS_H


/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: TRIM2DS.C                                                  #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #        Routines for computation of local element stresses for             #
 #        Linear triangular elements. Element code = 1.                      #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  trim2ds()    : computation of stresses in global directions.       #
 #                                                                           #
 #****************************************************************************/


#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include<fstream>
#include "Trimm2dmf.h"

/*******************************************[ USER LIBRARY : TRIM2DS.C  ]*****
 #                                                                           #
 #  FUNCTION :  trim2ds ()                                                   #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of element stresses in global coordinates.               #
 #                                                                           #
 #      Linear triangular elements.                                          #
 #      Element code = 1.                                                    #
 #                                                                           #
 #      Only the centroidal stresses are computed and stored in each node.   #
 #      Also the deformation energy is computed.                             #
 #                                                                           #
 #      Only the unrefined elements of each level are computed.              #
 #      (groups UNUMEL and URNUMEL).                                         #
 #                                                                           #
 #                                                                           #
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
void trim2ds (std::ofstream &fpo, std::ofstream &fGiD, int *acnodes, int  elm1, int elm2, int  xndof, int  nudim,
                 int *elm_sides, int *elm_nodes, int *elm_vlev, int *parent_elm, int *velm,
                 double **xyz, int    *mat, double *prop, int    xnumpr, double **uold,
                 double  *Fn, double  *Fs, double  **Sx, double  **Sy, double  **Sxy,
                 double *error, double errm, int ESTIMATOR, int  level);

 /* void trim2ds (FILE *fpo, FILE *fGiD, int *acnodes, int  elm1, int elm2, int  xndof, int  nudim,
                 int *elm_sides, int *elm_nodes, int *elm_vlev, int *parent_elm, int *velm,
                 double **xyz, int    *mat, double *prop, int    xnumpr, double **uold,
                 double  *Fn, double  *Fs, double  **Sx, double  **Sy, double  **Sxy,
                 double *error, double errm, int ESTIMATOR, int  level);*/


/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE   *fpo    ; /* pointer to output file                                 */
// FILE   *fGiD   ; /* pointer to GiD output file                             */
// int *acnodes   ; /* number of accumulated nodes in previous levels         */
// int  elm1, elm2; /* first and last elements of this level.                 */
// int  xndof     ; /* maximum number of nodal degrees of freedom.            */
// int  nudim     ; /* number of spatial dimensions.                          */
// int *elm_sides ; /* sides incidences of each element [numel * xnes]        */
// int *elm_nodes ; /* nodes of each element    [numel[curlev]*xnen]          */
// int *elm_vlev  ; /* level of element vertexs [numel[curlev]*xnev]          */
// int *parent_elm; /* parent element and location [3*numel[curlev]]          */
// int *velm      ; /* new element ordering after refinament[numel[curlev]]   */
// double **xyz   ; /* nodal coordinates [numnp * nudim].                     */
// int    *mat    ; /* element material type vector [numel].                  */
// double *prop   ; /* material properties vector [numat * numpr].            */
// int    xnumpr  ; /* maximum number of material properties.                 */
// double **uold  ; /* vector of nodal displacements [xndof*numnp]            */
// double  *Fn     ; /* midside normal traction       [sides[level]]           */
// double  *Fs     ; /* midside tangential traction   [sides[level]]           */
// double  **Sx  ; /* nodal stress Sx                           [nodes[level]] */
// double  **Sy  ; /* nodal stress Sy                           [nodes[level]] */
// double  **Sxy ; /* nodal stress Sxy                          [nodes[level]] */
// double *error  ; /* element error of this level   [numel[level]]           */
// double errm    ; /* mean error in energy norm of fine mesh                 */
// int ESTIMATOR; /* estimator type 1:nodal avg. 2:midside stress             */
// int  level     ; /* current level in analysis.                             */

/*********************************************************************************/


#endif
