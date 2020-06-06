/**********************[          TRIMBC2DK.h            ]********************/

#ifndef TRIM2DS_H
#define TRIM2DS_H


/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: TRIMBC2DK.C                                                #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #        Routines for assembling of global stiffness matrix for             #
 #        Linear triangular elements cubic Bezier shape.                     #
 #        Element code = 8.                                                  #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  trimbc2dk()  : computation of triangular element stiffness matrix. #
 #                                                                           #
 #****************************************************************************/


#include <cmath>
#include"SKYLINE.h"

/*******************************************[ USER LIBRARY : TRIM2DK.C  ]*****
 #                                                                           #
 #  FUNCTION :  trimbc2dk ()                                                 #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of element stiffness matrix in global coordinates        #
 #      and assembling into the global stiffness matrix.                     #
 #                                                                           #
 #      Linear triangular elements plane stress.                             #
 #      Element code = 8.                                                    #
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


  void trimbc2dk(int feblk  , int numel     , int ndof      , int nudim,
	             int *id    , int *elm_sides, int *side_nodes, int xnsn,
	             double *xyz, int *mat      , double *prop   , int numpr,
	             double *sg , int *maxa     , int *lm        , double *sv,
	             double *Gij, double *Kij   , double *diag);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

/* int     feblk;    first element number in this block.                    */
/* int     numel;    number of elements in this block.                      */
/* int     ndof ;    maximum number of nodal degrees of freedom.            */
/* int     nudim;    number of spatial dimensions.                          */
/* int    *id   ;    identification vector [numnp * ndof].                  */
/* int *side_nodes;  nodes of each side       [sides * xnsn]                */
/* int   xnsn   ;    maximum number of side nodes.                          */
/* int *elm_sides ;  sides of each element    [numel * xnes]                */
/* double *xyz  ;    nodal coordinates [numnp * nudim].                     */
/* int    *mat  ;    element material type vector [numel].                  */
/* double *prop ;    material properties vector [numat * numpr].            */
/* int     numpr;    maximum number of material properties.                 */
/* double  *sg  ;    global stiffness matrix stored by skyline [nwk]        */
/* int     *maxa;    pointer of diagonal elements for skyline [neq+1]       */
/* int     *lm  ;    identification vector of element dof [nrw]             */
/* double  *sv  ;    upper triangular element matrix [nrw *(nrw +1)/2]      */
/* double *Gij  ;    element stiffness coefficients [10*numel]              */
/* double *Kij  ;    element stiffness coefficients [21*numel]              */
/* double *diag ;    stiffness diagonal coefficients                        */

/********************************************************************************/

#endif
