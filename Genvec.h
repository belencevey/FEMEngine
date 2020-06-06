/**********************[             Genvec.h            ]********************/

#ifndef GENVEC_H
#define GENVEC_H

/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: GENVEC.C                                                   #
 #                                                                           #
 #****************************************************************************
 #                                                                           #
 #  void  nodvec()     : read/generates nodal vector data.                   #
 #  void  sidevec()    : read side vector surface load data.                 #
 #  void  elmvec()     : read element body load vector data.                 #
 #  void  sideload()   : computation of boundary nodal loads initial mesh.   #
 #  void  mlsideload() : multilevel computation of boundary nodal loads.     #
 #  void  elmld()      : computation of nodal loads over elements.           #
 #  void  genvec()     : generation of  nodal vector data.                   #
 #                                                                           #
 #****************************************************************************/


#include <cstdio>
#include <stdlib.h>
#include <cmath>
//#include <memory.h>
#include<string.h>
#include<fstream>
#include<iostream>

/*******************************************[ USER LIBRARY : GENVEC.C  ]******
 #                                                                           #
 #  FUNCTION :  nodvec()                                                     #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #    Input of nodal vector data.                                            #
 #    Optional data generation.                                              #
 #                                                                           #
 #    Each line of the input file is classified as:                          #
 #                                                                           #
 #     NODAL VECTOR DATA  LINE : contains vector data of one node.           #
 #     GENERATION COMMAND LINE : contains generation parameters.             #
 #                                                                           #
 #    The first number in a line of nodal vector data of the input file      #
 #    must be an integer value. If this value is NONZERO the line is         #
 #    assumed as a NODAL VECTOR DATA LINE  with the next syntax:             #
 #                                                                           #
 #          NNODE   F1  F2  F3  F4  ...  FXNDOF                              #
 #                                                                           #
 #    where                                                                  #
 #          NNODE : node number                                              #
 #             Fi : vector component in the direction of dof i applied       #
 #                  to node NNODE.                                           #
 #                                                                           #
 #                                                                           #
 #    the number of vector components must be equal to the maximum           #
 #    number of nodal degrees of freedom specified (parameter XNDOF).        #
 #    The program  accepts nodes with different numbers of dof, but          #
 #    for  nodes  with  less dof than  the maximum , the additional          #
 #    components MUST BE SPECIFIED with arbitrary values.                    #
 #                                                                           #
 #    If the first number in a line of nodal vector data of the input        #
 #    file is a ZERO integer value the line is assumed as a GENERATION       #
 #    COMMAND LINE with the next syntax:                                     #
 #                                                                           #
 #                                                                           #
 #          0  GCODE [generation parameters]                                 #
 #                                                                           #
 #    If GCODE is NONZERO a set of generation parameters must be given.      #
 #                                                                           #
 #                                                                           #
 #    TYPES OF NODAL VECTOR GENERATION                                       #
 #                                                                           #
 #    GCODE = 1 : assign a complete set of vectorial components.             #
 #    GCODE = 2 : increment a complete set of vectorial components.          #
 #                                                                           #
 #                                                                           #
 #    If GCODE is equal to ZERO, this line is interpreted as the last        #
 #    line of nodal coordinates data of the input file, and no generation    #
 #    parameters must be given.                                              #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED LIBRARY FUNCTIONS USED ]--------------------------------#
 #                                                                           #
 #                                                                           #
 #  GENVEC.C  :                                                              #
 #    void genvec() : generation of nodal vector data.                       #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H :                                                                #
 #    int  fscanf() : performs formatted input from a stream.                #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/

 void nodvec (std::ifstream &fpi, int xndof, int *id, double *f);
 //void nodvec (FILE *fpi, int xndof, int *id, double *f);

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE    *fpi ;   /* pointer input  file                                  */
// int     xndof;   /* maximum number of nodal degrees of freedom           */
// int     *id  ;   /* identification vector of numbered dof                */
// double  *f   ;   /* vector of nodal components [neq].                    */

/****************************************************************************/








/*******************************************[ USER LIBRARY : GENVEC.C  ]******
 #                                                                           #
 #  FUNCTION :  genvec ()                                                    #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #    Generation of nodal vector data.                                       #
 #                                                                           #
 #    The parameter GCODE defines the type of generation.                    #
 #                                                                           #
 #    TYPES OF NODAL VECTOR GENERATION                                       #
 #                                                                           #
 #    GCODE = 1 : ASSIGN A COMPLETE SET OF VECTORIAL COMPONENTS.             #
 #    GCODE = 2 : INCREMENT A COMPLETE SET OF VECTORIAL COMPONENTS.          #
 #                                                                           #
 #    PARAMETERS:  [ N1  N2  INC     F1  F2  F3  ...  FXNDOF ]               #
 #                                                                           #
 #             N1  : first node.                                             #
 #             N2  : last  node.                                             #
 #             INC : node number increment.                                  #
 #             Fi  : vector component in the direction of dof i applied      #
 #                   to this group of nodes.                                 #
 #                                                                           #
 #                                                                           #
 #    The number of  vector  components  must be equal to the maximum        #
 #    number of nodal degrees of freedom specified (parameter XNDOF).        #
 #                                                                           #
 #    Nodal vector components   F1, F2, F3, ..., FXNDOF, are assigned        #
 #    to a sequence of nodes beginning with node N1, and the subsequent      #
 #    nodes N1+INC, N1+2*INC, N1+3*INC,...., etc. Node N2 must be the        #
 #    last node of the sequence.                                             #
 #                                                                           #
 #    If N1>N2 then INC must be negative, but if it is given as a positive   #
 #    value, its sign is changed automatically.                              #
 #                                                                           #
 #                                                                           #
 #    If GCODE=1 then Fn  = Fi.                                              #
 #    If GCODE=2 then Fn += Fi.                                              #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H :                                                                #
 #    int  fscanf() : performs formatted input from a stream.                #
 #                                                                           #
 #  MATH.H  :                                                                #
 #    int    abs()  : gets the absolute value of an integer.                 #
 #    double fmod() : calculates the remainder of a division.                #
 #                                                                           #
 #****************************************************************************/

void genvec (std::ifstream &fpi, int gcode, int *id, int xndof, double *f);
 //void genvec (FILE *fpi, int gcode, int *id, int xndof, double *f);

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE    *fpi ;   /* pointer input  file                                  */
// int    gcode ;   /* generation code                                      */
// int     *id  ;   /* identification vector of dof status.                 */
//						/*      0 = free                                    */
//						/*      1 = restrained                              */
// int    xndof ;   /* maximum number of nodal degrees of freedom.          */
// double  *f   ;   /* vector of nodal components [neq].                    */

/****************************************************************************/




/*******************************************[ USER LIBRARY : GENVEC.C  ]******
 #                                                                           #
 #  FUNCTION :  sidevec()                                                    #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #    Input of side vector data.                                             #
 #    Optional data generation.                                              #
 #                                                                           #
 #    Each line of the input file is classified as:                          #
 #                                                                           #
 #     SIDE  VECTOR DATA  LINE : contains vector data of one node.           #
 #     GENERATION COMMAND LINE : contains generation parameters.             #
 #                                                                           #
 #    The first number in a line of nodal vector data of the input file      #
 #    must be an integer value. If this value is NONZERO the line is         #
 #    assumed as a SIDE VECTOR DATA LINE  with the next syntax:              #
 #                                                                           #
 #          NSIDE  DOF  F0  F1  F2  F3  ...  FXDF                            #
 #                                                                           #
 #    where                                                                  #
 #          NSIDE : node number                                              #
 #            DOF : direction of vector field                                #
 #                  DOF = 1 : global direction 1.                            #
 #                  DOF = 2 : global direction 2.                            #
 #             Fi : nodal values for a Lagrange's description of the         #
 #                  polynomial vector field over the side.                   #
 #                  Maximum degree of polynomial is XDF.                     #
 #                                                                           #
 #                                                                           #
 #    the number of nodal values must be equal to the maximum  degree        #
 #    of freedom specified plus one (parameter XDF + 1). These values        #
 #    correspond to nodes equally spaced along the side.                     #
 #    The program  accepts polynomial variations with less dof than the      #
 #    maximum, but the values on ALL the nodes loactions MUST BE SPECIFIED.  #
 #                                                                           #
 #    If the first number in a line of nodal vector data of the input        #
 #    file is a ZERO integer value the line is assumed as a GENERATION       #
 #    COMMAND LINE with the next syntax:                                     #
 #                                                                           #
 #                                                                           #
 #          0  GCODE [generation parameters]                                 #
 #                                                                           #
 #    If GCODE is NONZERO a set of generation parameters must be given.      #
 #                                                                           #
 #                                                                           #
 #    TYPES OF NODAL VECTOR GENERATION                                       #
 #                                                                           #
 #    GCODE = 1 : assign a complete set of vector components.                #
 #    GCODE = 2 : selective assignament of vector components.                #
 #                                                                           #
 #                                                                           #
 #    If GCODE is equal to ZERO, this line is interpreted as the last        #
 #    line of nodal coordinates data of the input file, and no generation    #
 #    parameters must be given.                                              #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED LIBRARY FUNCTIONS USED ]--------------------------------#
 #                                                                           #
 #                                                                           #
 #  GENVEC.C  :                                                              #
 #    void genvec() : generation of nodal vector data.                       #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H :                                                                #
 #    int  fscanf() : performs formatted input from a stream.                #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/

void sidevec (std::ifstream &fpi, int bsides, int xdf, int xnbdof, double **side_loads);
 //void sidevec (FILE *fpi, int bsides, int xdf, int xnbdof, double **side_loads);

/*____________________________________________________________________________

									  PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE    *fpi       ; /* pointer input  file                                */
// int      bsides    ; /* number of boundary sides.                          */
// int      xdf       ; /* maximum polynomial degree of external vector field.*/
// int      xnbdof    ; /* maximum number of boundary degrees of freedom      */
// double **side_loads; /* boundary side loads [xndof][bsides * (xdf+1)]      */

/***************************************************************************/






/*******************************************[ USER LIBRARY : GENVEC.C  ]******
 #                                                                           #
 #  FUNCTION :  sideload ()                                                  #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #    Computation of side nodal loads for initial level (coarse mesh).       #
 #                                                                           #
 #    Loads are assumed given in the tangential and normal directions of     #
 #    each boundary side.                                                    #
 #    Loads are transformed from the side local system to global x,y system  #
 #    and are transformed to the local system of each node.                  #
 #                                                                           #
 #    Polynomial description assumed for distributed boundary loads.         #
 #    The polynomial degree is XDF and the distribution is described         #
 #    by the nodal values at Lagrange points on the side.                    #
 #                                                                           #
 #    The maximum degree implemented is three.                               #
 #                                                                           #
 #    Nodal forces are computed on ALL the nodes of the side.                #
 #    These loads are stored in the load vector f[] and in compacted         #
 #    form (that is only on free dof) in the vector dn[].                    #
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


 void sideload (int ecode, int bsides, int xdf, int xndof, int xnsn,  int *side_nodes, double **side_loads,
				double *bs_Fs, double *bs_Fn, double *xyz, double *f, int *idb, double *bn_cs, double *bn_sn, double *bs_cs, double *bs_sn,
                int *side_mat, double *prop, int numpr, int nudim);

/*____________________________________________________________________________

									  PARAMETERS DECLARATION
____________________________________________________________________________*/

// int      ecode     ; /* element code number                                */
// int      bsides    ; /* number of boundary sides                           */
// int      xdf       ; /* maximum polynomial degree of external vector field */
// int      xndof     ; /* maximum number of nodal degrees of freedom         */
// int      xnsn      ; /* maximum number of side nodes                       */
// int     *side_nodes; /* nodes of each side initial level [sides * xnsn]    */
// double **side_loads; /* boundary side loads [xndof][bsides * (xdf+1)]      */
//                      /* 0:tangential 1:normal                              */
// double  *bs_Fs;      /* midside tangential boundary traction [Bsides]      */
// double  *bs_Fn;      /* midside normal boundary traction [Bsides]          */
// double  *xyz       ; /* nodal coordinates initial level [numnp * nudim]    */
// double  *f         ; /* vector of nodal forces [xndof * numnp]             */
// int    *idb    ;  /* identification vector of boundary nodes [numnp[0]]    */
// double *bn_cs  ;  /* boundary node cosine of local axes rotation angle     */
// double *bn_sn  ;  /* boundary node sine   of local axes rotation angle     */
// double *bs_cs  ;  /* boundary side cosine of local axes rotation angle     */
// double *bs_sn  ;  /* boundary side sine   of local axes rotation angle     */
// int    *side_mat; /* side material index of adjacent elements [2*sides[0]] */
// double *prop   ;  /* material properties   [numat * xnumpr]                */
// int     numpr;   /* maximum number of material properties.                 */
// int     nudim;  /* number of spatial dimensions (plus geometric weigths)   */

/***************************************************************************/







/*******************************************[ USER LIBRARY : GENVEC.C  ]******
 #                                                                           #
 #  FUNCTION :  mlsideload()                                                 #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #    Multilevel computation of boundary nodal loads.                        #
 #                                                                           #
 #    Polynomial description assumed for distributed boundary loads.         #
 #    The polynomial degree is XDF and the distribution is described         #
 #    by the nodal values at Lagrange points on the parent side.             #
 #                                                                           #
 #    The maximum degree implemented is three.                               #
 #                                                                           #
 #    Nodal forces are computed for all boundary nodes of all levels.        #
 #    These loads are stored in the load vector u1[].                        #
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


 void mlsideload (int ecode, int *nodes, int *Bsides, int *BFRsides, int xdf,
                  int xndof, int xnsn, int **bside, int **side_nodes, int **side_vlev, double **side_loads,
                  int **parent_side, int curlev, double **xyz, double **f, double **bs_cs, double **bs_sn);

/*____________________________________________________________________________

									  PARAMETERS DECLARATION
____________________________________________________________________________*/

// int       ecode    ; /* element code number                                */
// int      *nodes    ; /* number of nodes in each level                      */
// int      *Bsides   ; /* boundary sides of each level                       */
// int      *BFRsides ; /* number of refined boundary sides of group ursides. */
// int       xdf      ; /* maximum polynomial degree of external vector field */
// int       xndof    ; /* maximum number of nodal degrees of freedom         */
// int       xnsn     ; /* maximum number of side nodes                       */
// int    **bside     ; /* boundary side ordering bside[i]=old                */
// int    **side_nodes; /* nodes of each parent side [sides * xnsn]           */
// int    **side_vlev ; /* levels of side vertex nodes  [xnev*numel[curlev]]  */
// double **side_loads; /* boundary parent side loads [xndof][bsides*(xdf+1)] */
// int   **parent_side; /* parent side and location    [2*sides[curlev]]      */
// int      curlev    ; /* current level in analysis.                         */
// double **xyz       ; /* nodal coordinates of initial level.                */
// double **f         ; /* vector of nodal forces [xndof * numnp]             */
// double **bs_cs     ; /* boundary side cosine                               */
// double **bs_sn     ; /* boundary side sine                                 */

/****************************************************************************/


#endif
