/**********************[          SKYLINE.h            ]********************/

#ifndef SKYLINE_H
#define SKYLINE_H



/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME:  SKYLINE.C                                                 #
 #                                                                           #
 #****************************************************************************
 #                                                                           #
 #                                                                           #
 #       Routines for manipulation of skyline storage.                       #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  numeqn()   : assigns sequential equation numbers to each free dof. #
 #  void  expandv2D(): expand a compacted vector of vector into a global one.#
 #  void  expands()  : expand a compacted vector of scalar into a global one.#
 #  void  compressv2D(): compress global vector of vector into compacted one.#
 #  void  compresss(): compress global vector of scalars into compacted one. #
 #  void  profile()  : computes vector of diagonal elements of skyline.      #
 #  void  maxdiff()  : computes effective columns heights.                   #
 #  void  addban()   : assembling element matrix in compacted global matrix. #
 #                                                                           #
 #****************************************************************************/


#include <stdlib.h>


/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  numeqn ()                                                    #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Assigns sequential equation numbers to each free dof.                #
 #      The status of dof i is obtained from the identification              #
 #      vector id[i] following the next convention:                          #
 #                                                                           #
 #                  id[i]  = 0  free dof                                     #
 #                  id[i] != 0  restrained dof                               #
 #                                                                           #
 #      The number of equations associated with free dofs is                 #
 #      assigned to the parameter neq.                                       #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void numeqn (int numnp, int ndof, int *id, int *neq);

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int     numnp;   /* number of nodal points                                 */
// int     ndof ;   /* maximum number of nodal degrees of freedom             */
// int    *id   ;   /* identification vector [numnp * ndof]                   */
	              /* before processing dof status (0 = free, 1 = restrained)*/
		          /* after  processing numbering of free dof                */
// int    *neq  ;   /* number of equations                                    */

/*******************************************************************************/




/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  expandv2D()                                                  #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Expand a compacted vector dn of 2D vectors into a complete vector u. #
 #      The number of dof i is obtained from the identification              #
 #      vector id[i] following the next convention:                          #
 #                                                                           #
 #                  id[i]  = 0  u[i] = 0                                     #
 #                  id[i] != 0  u[i] = dn[id[i]]                             #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void expandv2D (double *dn, double *u, int *id, int numnp, int xndof, int *idb, double *bn_cs, double *bn_sn);

/*____________________________________________________________________________

								PARAMETERS DECLARATION
____________________________________________________________________________*/

// double *dn   ;  /* compacted vector to be expanded                         */
// double *u    ;  /* expanded vector                                         */
// int    *id   ;  /* identification vector of free dof [numnp * xndof]       */
// int     numnp;  /* number of nodal points                                  */
// int     xndof;  /* maximum number of nodal degrees of freedom              */
// int    *idb    ;  /* identification vector of boundary nodes [numnp[0]]    */
// double *bn_cs  ;  /* boundary node cosine of local axes rotation angle     */
// double *bn_sn  ;  /* boundary node sine   of local axes rotation angle     */

/*****************************************************************************/





/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  expands()                                                    #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Expand a compacted vector dn of scalars into a complete vector u.    #
 #      The number of dof i is obtained from the identification              #
 #      vector id[i] following the next convention:                          #
 #                                                                           #
 #                  id[i]  = 0  u[i] = 0                                     #
 #                  id[i] != 0  u[i] = dn[id[i]]                             #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void expands (double *dn, double *u, int *id, int numnp);

/*____________________________________________________________________________

								PARAMETERS DECLARATION
____________________________________________________________________________*/

// double *dn   ;  /* compacted vector to be expanded                         */
// double *u    ;  /* expanded vector                                         */
// int    *id   ;  /* identification vector of free dof [numnp * xndof]       */
// int     numnp;  /* number of nodal points                                  */

/******************************************************************************/






/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  compressv2D()                                                #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Compress a global vector u of 2D vectors into a compacted vector dn. #
 #      Rotates global components to local node system.                      #
 #      The number of dof i is obtained from the identification              #
 #      vector id[i] following the next convention:                          #
 #                                                                           #
 #                  id[i]  = 0  u[i] = 0                                     #
 #                  id[i] != 0  u[i] = dn[id[i]]                             #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void compressv2D (double *dn, double *u, int *id, int numnp, int xndof, int *idb, double *bn_cs, double *bn_sn);


/*____________________________________________________________________________

								PARAMETERS DECLARATION
____________________________________________________________________________*/

// double *dn   ;  /* compacted vector                                        */
// double *u    ;  /* global vector to be compressed                          */
// int    *id   ;  /* identification vector of free dof [numnp * xndof]       */
// int     numnp;  /* number of nodal points                                  */
// int     xndof;  /* maximum number of nodal degrees of freedom              */
// int    *idb    ;  /* identification vector of boundary nodes [numnp[0]]    */
// double *bn_cs  ;  /* boundary node cosine of local axes rotation angle     */
// double *bn_sn  ;  /* boundary node sine   of local axes rotation angle     */

/******************************************************************************/






/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  compresss()                                                  #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Compress a global vector u of scalars into a compacted vector dn.    #
 #      Rotates global components to local node system.                      #
 #      The number of dof i is obtained from the identification              #
 #      vector id[i] following the next convention:                          #
 #                                                                           #
 #                  id[i]  = 0  u[i] = 0                                     #
 #                  id[i] != 0  u[i] = dn[id[i]]                             #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void compresss (double *dn, double *u, int *id, int numnp);


/*____________________________________________________________________________

								PARAMETERS DECLARATION
____________________________________________________________________________*/

// double *dn   ;  /* compacted vector                                        */
// double *u    ;  /* global vector to be compressed                          */
// int    *id   ;  /* identification vector of free dof [numnp * xndof]       */
// int     numnp;  /* number of nodal points                                  */

/*****************************************************************************/






/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  profile ()                                                   #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computes effective column heights of skyline.                        #
 #      These heights are the maximum difference between the                 #
 #      equation number of each dof and the least element                    #
 #      equation number, for all the elements that contains                  #
 #      that dof. The height of the dof i is stored in the                   #
 #      location maxa[i+1].                                                  #
 #      Computes adresses of diagonal elements in skyline and                #
 #      stores them in the same array maxa[].                                #
 #      The elements of a column of the skyline are sequentially             #
 #      numbered from bottom to top, and the position number of              #
 #      of the diagonal elements of the matrix (bottom elements              #
 #      of columns) in the skyline is computed from the column               #
 #      heights as:                                                          #
 #                                                                           #
 #              maxa[1] = 1                                                  #
 #              maxa[2] = 2                                                  #
 #                                                                           #
 #              for i>2  maxa[i] = maxa[i] + maxa[i-1] + 1;                  #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


 void profile (int neq, int numel, int sides, int bsides, int ecode, int xndof, int *id, int *elm_sides,
               int *side_nodes, int *lm, int *nwk, int *maxa);

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int  neq       ;  /* number of equations                                   */
// int  numel     ;  /* number of elements                                    */
// int  sides     ;  /* number of sides.                                      */
// int  bsides    ;  /* number of boundary sides.                             */
// int  ecode     ;  /* element code number                                   */
// int  xndof     ;  /* maximum number of nodal degrees of freedom            */
// int *id        ;  /* equation number of free dof [numnp * ndof]            */
						             /* ( >0 eqn. number, =0 restrained)      */
// int *side_nodes;  /* nodes of each side       [sides * xnsn]               */
// int *elm_sides ;  /* sides of each element    [numel * xnes]               */
// int *lm        ;  /* equation number of element dof [nen*ndof]             */
// int *nwk       ;  /* number of matrix elements                             */
// int *maxa      ;  /* pointer of diagonal elements for skyline [neq+1]      */

/******************************************************************************/





/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  maxdiff ()                                                   #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computes effective column heights of skyline.                        #
 #      These heights are the maximum difference between the                 #
 #      equation number of each dof and the least element                    #
 #      equation number, for all the elements that contains                  #
 #      that dof. The height of the dof i is stored in the                   #
 #      location maxa[i+1].                                                  #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/

#include <cmath>

 void maxdiff (int neq, int numel, int sides, int bsides, int ecode, int xndof, int *id, int *elm_sides,
               int *side_nodes, int *lm, int *maxa);

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int  neq       ;  /* number of equations                                   */
// int  numel     ;  /* number of elements                                    */
// int  sides     ;  /* number of sides.                                      */
// int  bsides    ;  /* number of boundary sides.                             */
// int  ecode     ;  /* element code number                                   */
// int  xndof     ;  /* maximum number of nodal degrees of freedom            */
// int *id        ;  /* equation number of free dof [numnp * ndof]            */
						 /* ( >0 eqn. number, =0 restrained)                      */
// int *side_nodes;  /* nodes of each side       [sides * xnsn]               */
// int *elm_sides ;  /* sides of each element    [numel * xnes]               */
// int *lm        ;  /* equation number of element dof [nen*ndof] [nsn*ndof]  */
// int *maxa      ;  /* pointer of diagonal elements for skyline [neq+1]      */

/******************************************************************************/





/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  addban ()                                                    #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Assembling of an upper triangular element matrix into                #
 #      the compacted global matrix.                                         #
 #                                                                           #
 #       (translation to C language of subroutine ADDBAN from                #
 #        K.J.Bathe, "Finite Element Procedures in Engineering               #
 #        Analysis", Prentice Hall, 1982.)                                   #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void addban (double *A, int *maxa, double *sv, int *lm, int numrw);

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// double  *A   ;   /* global matrix stored by skyline [nwk]                  */
// int     *maxa;   /* pointer of diagonal elements for skyline [neq+1]       */
// double  *sv  ;   /* upper triagular element matrix stored by rows          */
// int     *lm  ;   /* equation numbers of element dof [nen*ndof]             */
// int     numrw;   /* number of element square matrix rows/columns           */

/******************************************************************************/



#endif
