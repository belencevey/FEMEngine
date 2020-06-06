// /**********************[          TR2DMLMF.h            ]********************/

#ifndef TR2DMLMF_H
#define TR2DMLMF_H

// /* *********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
// #                                                                           #
// #  LIBRARY NAME: TR2DMLMF.C                                                 #
// #                                                                           #
// *****************************************************************************
// #                                                                           #
// #                                                                           #
// #        Routines for computation of midside interpolated forces.           #
// #        All elements of current mesh for all levels.                       #
// #        Linear triangular elements. Element code = 1.                      #
// #                                                                           #
// #                                                                           #
// #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
// #                                                                           #
// #  void  tr2dmlmf()  : Multilevel computation of midside forces             #
// #                      linear triangle.                                     #
// #                                                                           #
// #*************************************************************************** */


#include <cstdio>
#include <cmath>
#include <stdlib.h>
#include <memory>
//#include <mem.h>

// /* ******************************************[ USER LIBRARY : TR2DMLMF.C ]*****
// #                                                                           #
// #  FUNCTION :  tr1_nodal_avg()                                              #
// #                                                                           #
// #  TYPE     :  void                                                         #
// #                                                                           #
// #                                                                           #
// #----[ DESCRIPTION ]--------------------------------------------------------#
// #                                                                           #
// #      Computation of side interpolated forces for linear triangles.        #
// #      Multilevel computations for ALL elements of the mesh.                #
// #      Linear variation for the stresses along the side is assumed.         #
// #      Average value is assigned for sides belonging to two elements.       #
// #      Special averaging for transtion sides.                               #
// #      For boundary sides, external tractions are impossed.                 #
// #                                                                           #
// #      Linear triangular elements.                                          #
// #      Element code = 1.                                                    #
// #                                                                           #
// #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
// #                                                                           #
// #                                                                           #
// #  MATH.H  :                                                                #
// #    double sqrt() : calculates square root.                                #
// #                                                                           #
// #  STDLIB.H:                                                                #
// #    int abs()     : macro that return the absolute value of integer.       #
// #                                                                           #
// #                                                                           #
// #************************************************************************** */


  void tr1_nodal_avg (int curlev, int *numel, int *nodes, int *FRnumel, int **velm,
                  int *BFRsides,int *PRsides, int **eside, int **bside, int **tside, int xndof,
                  double **xyz, int **elm_nodes, int **elm_vlev, int **parent_elm, int **side_nodes,
                  int **side_vlev, double **uold, int *mat, double *prop,int xnumpr,double **Sx,
                  double **Sy, double **Sxy, int *nse, int nudim);


/*____________________________________________________________________________

							 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int   curlev      ; /* current level in analysis.                          */
// int  *numel       ; /* number of elements in each level.                   */
// int   *nodes      ; /* number of nodal points of each level.               */
// int  *FRnumel     ; /* number of refined elements (all in group rnumel).   */
// int  **velm       ; /* new element ordering after refinament [numel[level]]*/
// int   *BFRsides   ; /* number of refined boundary sides in each level.     */
// int   *PRsides    ; /* number of transition sides in each level.           */
// int  **eside      ; /* inverse side ordering                               */
// int  **bside      ; /* list of boundary sides of each level.               */
// int  **tside      ; /* list of transition sides of each side.              */
// int    xndof      ; /* maximum number of nodal degrees of freedom.         */
// double **xyz      ; /* nodal coordinates     [numnp[curlev] * nudim]       */
// int  **elm_nodes  ; /* nodes of each element    [numel[curlev]*xnen]       */
// int  **elm_vlev   ; /* level of element vertexs [numel * xnev]             */
// int  **parent_elm ; /* parent element and location [3*numel[curlev]]       */
// int  **side_nodes ; /* nodes of each side       [sides[curlev] * xnsn]     */
// int  **side_vlev  ; /* levels of side vertex nodes  [xnev*numel[curlev]]   */
// double **uold     ; /* vector of previous nodal displacements [xndof*numnp]*/
// int    *mat       ; /* element material type vector [numel].               */
// double *prop      ; /* material properties vector [numat * numpr].         */
// int     xnumpr    ; /* maximum number of material properties.              */
// double  **Sx  ; /* nodal stress Sx                           [nodes[level]] */
// double  **Sy  ; /* nodal stress Sy                           [nodes[level]] */
// double  **Sxy ; /* nodal stress Sxy                          [nodes[level]] */
// int    *nse    ;  /* number of elements at each node of coarse mesh        */
// int   nudim;

/* *****************************************************************************/







/* ******************************************[ USER LIBRARY : TR2DMLMF.C ]*****
 #                                                                           #
 #  FUNCTION :  tr1_flux_avg()                                               #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of side interpolated forces for linear triangles.        #
 #      Multilevel computations for ALL elements of the mesh.                #
 #      Linear variation for the stresses along the side is assumed.         #
 #      Average value is assigned for sides belonging to two elements.       #
 #      Special averaging for transtion sides.                               #
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


  void tr1_flux_avg (int curlev, int *numel, int *FRnumel, int **velm, int *sides, int *Bsides,
                 int *BFRsides, int*FRsides, int *PRsides, int  **vside, int  **bside, int  **tside,
                 int  **ssside, int xndof, double **xyz, int nudim, int **elm_sides, int **elm_nodes,
                 int **elm_vlev, int **parent_elm, int **parent_side, int *side_bc, double **uold, int *mat,
                 double *prop, int xnumpr, double **Fn, double **Fs, double **bs_Fn, double **bs_Fs, double **fi,
                 double **bn_nlength, double **bn_slength, double **bn_qn, double **bn_qs, double **bs_cs,
                 double **bs_sn, int *nodes, int **side_nodes, int **side_vlev);


/*____________________________________________________________________________

							 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int   curlev      ; /* current level in analysis.                          */
// int  *numel       ; /* number of elements in each level.                   */
// int  *FRnumel     ; /* number of refined elements (all in group rnumel).   */
// int  **velm       ; /* new element ordering after refinament [numel[level]]*/
// int   *sides      ; /* number of sides in each level.                      */
// int   *Bsides     ; /* number of boundary sides in each level.             */
// int   *BFRsides   ; /* number of refined boundary sides in each level.     */
// int   *FRsides    ; /* number of refined sides in each level.              */
// int   *PRsides    ; /* number of transition sides in each level.           */
// int  **vside      ; /* new side ordering after refinament [sides[level]]   */
// int  **bside      ; /* list of boundary sides of each level.               */
// int  **tside      ; /* list of transition sides of each side.              */
// int **ssside      ;
// int    xndof      ; /* maximum number of nodal degrees of freedom.         */
// double **xyz      ; /* nodal coordinates     [numnp[curlev] * nudim]       */
// int    nudim      ; /* number of spatial dimensions (plus geometric weigths)   */
// int  **elm_sides  ; /* sides incidences of each element [numel * xnes]     */
// int  **elm_nodes  ; /* nodes of each element    [numel[curlev]*xnen]       */
// int  **elm_vlev   ; /* level of element vertexs [numel * xnev]             */
// int  **parent_elm ; /* parent element and location [3*numel[curlev]]       */
// int  **parent_side; /* parent side and location [2*sides[curlev]]          */
// int    *side_bc   ; /* side boundary conditions [Bsides[0] * xndof]        */
// double **uold     ; /* vector of previous nodal displacements [xndof*numnp]*/
// int    *mat       ; /* element material type vector [numel].               */
// double *prop      ; /* material properties vector [numat * numpr].         */
// int     xnumpr    ; /* maximum number of material properties.              */
// double  **Fn       ; /* side normal tractions [sides[curlev]]               */
// double  **Fs       ; /* side tangential tractons [sides[curlev]]            */
// double **bs_Fn    ; /* midside normal boundary tractions current level     */
// double **bs_Fs    ; /* midside tangential boundary tractions current level */
// double **fi;
// double **bn_nlength; /* boundary node summ of adjacent lengths normal      */
// double **bn_slength; /* boundary node summ of adjacent lengths tangent     */
// double **bn_qn    ; /* boundary node normal reaction                       */
// double **bn_qs    ; /* boundary node tangent reaction                      */
// double **bs_cs    ; /* boundary side cosine of local system angle          */
// double **bs_sn    ; /* boundary side sine   of local system angle          */
// int *nodes;
// int **side_nodes;
// int **side_vlev;

/* ******************************************************************************/





/* ******************************************[ USER LIBRARY : TR2DMLMF.C ]*****
 #                                                                           #
 #  FUNCTION :  tr1_zz_centsp()                                              #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of nodal estimated stresses for linear triangles.        #
 #      Linear variation for the stresses over elements is assumed.          #
 #      ZZ superconvegence VERTEX recovery CENTROIDAL sampling.              #
 #      Average value is assigned for boundary nodes.                        #
 #      Special averaging for transtion nodes.                               #
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


  void tr1_zz_centsp (int   curlev, int  *numel, int  *nodes, int  *FRnumel, int  **velm, int  *Bsides,
                 int *BFRsides, int *PRsides, int  **eside, int  **bside, int  **tside, int xndof,
                 double **xyz, int **elm_nodes, int **elm_vlev, int **parent_elm, int **side_nodes,
                 int **side_vlev, double **uold, int *mat, double *prop, int xnumpr, double **Sx,
                 double **Sy, double **Sxy, double **a11, double **a12, double **a13, double **a22,
                 double **a23, double **a33, double **b1x, double **b2x, double **b3x, double **b1y,
                 double **b2y, double **b3y, double **b1xy, double **b2xy, double **b3xy, int *nse,
                 int **idnode);


/*____________________________________________________________________________

							 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int   curlev      ; /* current level in analysis.                          */
// int  *numel       ; /* number of elements in each level.                   */
// int   *nodes      ; /* number of nodal points of each level.               */
// int  *FRnumel     ; /* number of refined elements (all in group rnumel).   */
// int  **velm       ; /* new element ordering after refinament [numel[level]]*/
// int   *Bsides     ; /* number of boundary sides of each level.             */
// int   *BFRsides   ; /* number of refined boundary sides in each level.     */
// int   *PRsides    ; /* number of transition sides in each level.           */
// int  **eside      ; /* inverse side ordering                               */
// int  **bside      ; /* list of boundary sides of each level.               */
// int  **tside      ; /* list of transition sides of each side.              */
// int    xndof      ; /* maximum number of nodal degrees of freedom.         */
// double **xyz      ; /* nodal coordinates     [numnp[curlev] * nudim]       */
// int  **elm_nodes  ; /* nodes of each element    [numel[curlev]*xnen]       */
// int  **elm_vlev   ; /* level of element vertexs [numel * xnev]             */
// int  **parent_elm ; /* parent element and location [3*numel[curlev]]       */
// int  **side_nodes ; /* nodes of each side       [sides[curlev] * xnsn]     */
// int  **side_vlev  ; /* levels of side vertex nodes  [xnev*numel[curlev]]   */
// double **uold     ; /* vector of previous nodal displacements [xndof*numnp]*/
// int    *mat       ; /* element material type vector [numel].               */
// double *prop      ; /* material properties vector [numat * numpr].         */
// int     xnumpr    ; /* maximum number of material properties.              */
// double  **Sx  ; /* nodal stress Sx                           [nodes[level]] */
// double  **Sy  ; /* nodal stress Sy                           [nodes[level]] */
// double  **Sxy ; /* nodal stress Sxy                          [nodes[level]] */
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
// int    *nse    ;  /* number of elements at each node of coarse mesh        */
// int **idnode  ; /* identification of boundary nodes.                       */

/* *******************************************************************************/




/* ******************************************[ USER LIBRARY : TR2DMLMF.C ]*****
 #                                                                           #
 #  FUNCTION :  tr1_zz_sidesp()                                              #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of nodal estimated stresses for linear triangles.        #
 #      Linear variation for the stresses over elements is assumed.          #
 #      ZZ superconvegence VERTEX recovery   SIDE   sampling.                #
 #      Average value is assigned for boundary nodes.                        #
 #      Special averaging for transtion nodes.                               #
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


  void tr1_zz_sidesp (int curlev, int *numel, int *nodes, int *FRnumel, int **velm, int *Bsides,
                  int *BFRsides, int *PRsides, int **eside, int **bside, int **tside, int xndof,
                  double **xyz, int **elm_nodes, int **elm_vlev, int **parent_elm, int **side_nodes,
                  int **side_vlev, double **uold, int *mat, double *prop, int xnumpr, double **Sx,
                  double **Sy, double **Sxy, double **a11, double **a12, double **a13, double **a22,
                  double **a23, double **a33, double **b1x, double **b2x, double **b3x, double **b1y,
                  double **b2y, double **b3y, double **b1xy, double **b2xy, double **b3xy, int *nse, int **idnode);


/*____________________________________________________________________________

							 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int   curlev      ; /* current level in analysis.                          */
// int  *numel       ; /* number of elements in each level.                   */
// int   *nodes      ; /* number of nodal points of each level.               */
// int  *FRnumel     ; /* number of refined elements (all in group rnumel).   */
// int  **velm       ; /* new element ordering after refinament [numel[level]]*/
// int   *Bsides     ; /* number of boundary sides of each level.             */
// int   *BFRsides   ; /* number of refined boundary sides in each level.     */
// int   *PRsides    ; /* number of transition sides in each level.           */
// int  **eside      ; /* inverse side ordering                               */
// int  **bside      ; /* list of boundary sides of each level.               */
// int  **tside      ; /* list of transition sides of each side.              */
// int    xndof      ; /* maximum number of nodal degrees of freedom.         */
// double **xyz      ; /* nodal coordinates     [numnp[curlev] * nudim]       */
// int  **elm_nodes  ; /* nodes of each element    [numel[curlev]*xnen]       */
// int  **elm_vlev   ; /* level of element vertexs [numel * xnev]             */
// int  **parent_elm ; /* parent element and location [3*numel[curlev]]       */
// int  **side_nodes ; /* nodes of each side       [sides[curlev] * xnsn]     */
// int  **side_vlev  ; /* levels of side vertex nodes  [xnev*numel[curlev]]   */
// double **uold     ; /* vector of previous nodal displacements [xndof*numnp]*/
// int    *mat       ; /* element material type vector [numel].               */
// double *prop      ; /* material properties vector [numat * numpr].         */
// int     xnumpr    ; /* maximum number of material properties.              */
// double  **Sx  ; /* nodal stress Sx                           [nodes[level]] */
// double  **Sy  ; /* nodal stress Sy                           [nodes[level]] */
// double  **Sxy ; /* nodal stress Sxy                          [nodes[level]] */
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
// int    *nse    ;  /* number of elements at each node of coarse mesh        */
// int **idnode  ; /* identification of boundary nodes.                       */

/* *****************************************************************************/







/* ******************************************[ USER LIBRARY : TR2DMLMF.C ]*****
 #                                                                           #
 #  FUNCTION :  tr1_side_avg()                                               #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of side interpolated forces for linear triangles.        #
 #      Multilevel computations for ALL elements of the mesh.                #
 #      Linear variation for the stresses along the side is assumed.         #
 #      Averaging of global side stresses.                                   #
 #      Special averaging for transtion nodes.                               #
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


  void tr1_side_avg (
                         int curlev, int  *numel, int  *FRnumel, int  **velm, int *sides,
                         int *Bsides, int *BFRsides, int *FRsides, int *PRsides, int  **vside,
                         int **bside, int **tside, int **ssside, int xndof, double **xyz, int **elm_sides,
                         int **elm_nodes, int **elm_vlev, int **parent_elm, double **uold, int *mat,
                         double *prop, int xnumpr, double **Sx, double  **Sy, double **Sxy);


/*____________________________________________________________________________

							 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int   curlev      ; /* current level in analysis.                          */
// int  *numel       ; /* number of elements in each level.                   */
// int  *FRnumel     ; /* number of refined elements (all in group rnumel).   */
// int  **velm       ; /* new element ordering after refinament [numel[level]]*/
// int   *sides      ; /* number of sides in each level.                      */
// int   *Bsides     ; /* number of boundary sides in each level.             */
// int   *BFRsides   ; /* number of refined boundary sides in each level.     */
// int   *FRsides    ; /* number of refined sides in each level.              */
// int   *PRsides    ; /* number of transition sides in each level.           */
// int  **vside      ; /* new side ordering after refinament [sides[level]]   */
// int  **bside      ; /* list of boundary sides of each level.               */
// int  **tside      ; /* list of transition sides of each side.              */
// int **ssside      ;
// int    xndof      ; /* maximum number of nodal degrees of freedom.         */
// double **xyz      ; /* nodal coordinates     [numnp[curlev] * nudim]       */
// int  **elm_sides  ; /* sides incidences of each element [numel * xnes]     */
// int  **elm_nodes  ; /* nodes of each element    [numel[curlev]*xnen]       */
// int  **elm_vlev   ; /* level of element vertexs [numel * xnev]             */
// int  **parent_elm ; /* parent element and location [3*numel[curlev]]       */
// double **uold     ; /* vector of previous nodal displacements [xndof*numnp]*/
// int    *mat       ; /* element material type vector [numel].               */
// double *prop      ; /* material properties vector [numat * numpr].         */
// int     xnumpr    ; /* maximum number of material properties.              */
// double  **Sx  ; /* nodal stress Sx                           [nodes[level]] */
// double  **Sy  ; /* nodal stress Sy                           [nodes[level]] */
// double  **Sxy ; /* nodal stress Sxy                          [nodes[level]] */

/* ******************************************************************************/

#endif
