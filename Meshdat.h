/**********************[          Meshdat.h            ]********************/

#ifndef MESHDAT_H
#define MESHDAT_H

/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: MESHDAT.C                                                  #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #                Routines for input of  mesh data.                          #
 #                Nodal coordinates. Side boundary conditions.               #
 #                Side and Element properties.                               #
 #                Only one type of element for the entire mesh.              #
 #                Linear static analysis.                                    #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  meshdat()    : read mesh data.                                     #
 #  void  coord()      : read nodal coordinates.                             #
 #  void  sideinc()    : read side incidences.                               #
 #  void  sidebc()     : read side boundary conditions.                      #
 #  void  matprop()    : read material properties.                           #
 #  void  elminc()     : read element incidences.                            #
 #                                                                           #
 #****************************************************************************/


#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include<fstream>
#include<iostream>

/*******************************************[ USER LIBRARY : MESHDAT.C ]******
 #                                                                           #
 #  FUNCTION :  meshdat ()                                                   #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Input of fine mesh data for two dimensional meshes.                  #
 #      Nodal coordinates and side boundary conditions.                      #
 #      Nodal incidences of sides.                                           #
 #      Element material type and incidences.                                #
 #      Material properties.                                                 #
 #                                                                           #
 #      There are not implemented generation commands.                       #
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
 #  MESHDAT.C :                                                              #
 #    void coord()       : read nodal coordinates.                           #
 #    void sideinc()     : read side incidences.                             #
 #    void sidebc()      : read side boundary conditions.                    #
 #    void matprop()     : read material properties.                         #
 #    void  elminc()     : read element incidences.                          #
 #    void  elmnod()     : computation of element nodal incidences.          #
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


  void  meshdat(std::ifstream &fpi      , int nodes     , int numel,
                       int sides      , int Bsides    , int nudim,
                       int ecode      , double *xyz   , double *xyzc,
                       int *id        , int xndof     , int xnbdof,
					   int xnsn       , int xnes      ,
                       int *side_nodes, int *side_vlev, int *elm_sides,
                       int *elm_nodes , int *elm_vlev , int *side_bc,
                       int *side_mat,
                       int numat      , int *mat      , double *prop,
                       int xnumpr     , int *idb      , double *bn_cs,
                       double *bn_sn  , int numcurv   , double *bs_cs,
                       double *bs_sn  ,
                       double *bn_nlength, double *bn_slength);


 /*
  void  meshdat(FILE *fpi      , int nodes     , int numel,
                       int sides      , int Bsides    , int nudim,
                       int ecode      , double *xyz   , double *xyzc,
                       int *id        , int xndof     , int xnbdof,
					   int xnsn       , int xnes      ,
                       int *side_nodes, int *side_vlev, int *elm_sides,
                       int *elm_nodes , int *elm_vlev , int *side_bc,
                       int *side_mat,
                       int numat      , int *mat      , double *prop,
                       int xnumpr     , int *idb      , double *bn_cs,
                       double *bn_sn  , int numcurv   , double *bs_cs,
                       double *bs_sn  ,
                       double *bn_nlength, double *bn_slength);
*/

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/


// FILE    *fpi   ; /* pointer input  file.                                   */
// int    numnp   ; /* number of nodal points.                                */
// int    numel   ; /* number of elements.                                    */
// int    sides   ; /* number of sides.                                       */
// int   bsides   ; /* number of boundary sides.                              */
// int    nudim   ; /* number of spatial dimensions.                          */
// int    ecode   ; /* element code number.                                   */
// double  *xyz   ; /* nodal coordinates [numnp * nudim].                     */
// double  *xyzc  ; /* coordinates of control points [numnp * nudim].         */
// int     *id    ; /* identification vector [numnp * xndof].                 */
// int    xndof   ; /* maximum number of nodal degrees of freedom.            */
// int    xnbdof  ; /* maximum number of boundary degrees of freedom.         */
// int     xnsn   ; /* maximum number of side nodes.                          */
// int     xnes   ; /* maximum number of element sides.                       */
// int *side_nodes; /* nodes of each side       [sides * xnsn]                */
// int *side_vlev ; /* level of side vertexs    [sides * xnsv]                */
// int *elm_sides ; /* sides of each element    [numel * xnes]                */
// int *elm_nodes ; /* nodes of each element    [numel * xnes]                */
// int *elm_vlev  ; /* level of element vertexs [numel * xnev]                */
// int *side_bc   ; /* side boundary conditions [bsides * xndof]              */
// int *side_mat  ; /* side material index of adjacent elements [2*sides[0]]  */
// int    numat   ; /* number of material sets.                               */
// int    *mat    ; /* element material type vector [numel].                  */
// double *prop   ; /* material properties vector [numat * numpr].            */
// int   xnumpr   ; /* maximum number of material properties.                 */
// int    *idb    ; /* identification vector of boundary nodes [numnp[0]]     */
// int numcurv    ;  /* number of curved boundary sides initial mesh          */
// double *bs_cs  ; /* boundary side cosine of local axes rotation angle      */
// double *bs_sn  ; /* boundary side sine   of local axes rotation angle      */
// double *bn_cs  ; /* boundary node cosine of local axes rotation angle      */
// double *bn_sn  ; /* boundary node sine   of local axes rotation angle      */
// double *bn_nlength; /* boundary node summ of adjacent lengths normal       */
// double *bn_slength; /* boundary node summ of adjacent lengths tangent      */

/****************************************************************************/





/*******************************************[ USER LIBRARY : MESHDAT.C ]******
 #                                                                           #
 #  FUNCTION :  coord ()                                                     #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #    Input of nodal coordinates.                                            #
 #                                                                           #
 #    A line of NODAL COORDINATES data has the next syntax:                  #
 #                                                                           #
 #          NNODE   Xi   Yi  [Zi]                                            #
 #                                                                           #
 #    where                                                                  #
 #          NNODE : node number                                              #
 #             Xi : coordinate x of node NNODE                               #
 #             Yi : coordinate y of node NNODE                               #
 #             Zi : coordinate z of node NNODE                               #
 #                                                                           #
 #    the brackets indicates that this data is condition dependent,          #
 #    that is, depending on the number of spatial dimensions specified       #
 #    (parameter NUDIM) the following data must be given:                    #
 #                                                                           #
 #           NUDIM=2 : only x, y coordinates must be given.                  #
 #           NUDIM=3 : the three x, y, z coordinates must be given.          #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H :                                                                #
 #    int  fscanf() : performs formatted input from a stream.                #
 #    int  printf() : formatted output to stdout.                            #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void coord (FILE *fpi, int  numnp, int nudim, double *xyz);

/*____________________________________________________________________________

									PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE    *fpi ;   /* pointer input  file                                    */
// int     numnp;   /* number of nodal points.                                */
// int     nudim;   /* number of spatial dimensions.                          */
// double  *xyz ;   /* nodal coordinates.                                     */

/****************************************************************************/




/*******************************************[ USER LIBRARY : MESHDAT.C ]******
 #                                                                           #
 #  FUNCTION :  nodebc ()                                                    #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #    Input of nodal boundary conditions for initial coarse mesh.            #
 #    Computation of cosine, sine of local angle of rotation anticlockwse.   #
 #                                                                           #
 #    A line of NODAL BOUNDARY CONDITIONS data has the next syntax:          #
 #                                                                           #
 #          BNODE   ANGLE   R1  R2  R3  R4  ...  RXNDOF                      #
 #                                                                           #
 #    where                                                                  #
 #          BNODE : boundary node number.                                    #
 #          ANGLE : rotation angle anticlockwise of node local system.       #
 #             Ri : restraint condition of degree of freedom (dof) i         #
 #                  of node BNODE in the local system.                       #
 #                                                                           #
 #    each restraint condition Ri must be an integer 0 or 1 which is         #
 #    interpreted as:                                                        #
 #                                                                           #
 #             Ri=1  : restrained.                                           #
 #             Ri=0  : free.                                                 #
 #                                                                           #
 #    The number of restraint conditions must be equal to the maximum        #
 #    number of nodal degrees of freedom specified (parameter XNDOF).        #
 #                                                                           #
 #    The program  accepts nodes with different number of dof,  but          #
 #    internally all the nodes are considered as having equal number         #
 #    of dofs. Then, for nodes with less dofs than the maximum , the         #
 #    additional dofs MUST GIVEN RESTRAINED.??                               #
 #                                                                           #
 #    Also, nodes with are only used for auxiliar purposes such as mesh      #
 #    generation, specification of axis orientation, etc, MUST BE            #
 #    COMPLETELY RESTRAINED.                                                 #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H :                                                                #
 #    int  fscanf() : performs formatted input from a stream.                #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void nodebc (FILE *fpi, int ecode, int numnp, int bnodes, int xndof, int *id, int *idb, double *bn_cs, double *bn_sn);

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE   *fpi    ; /* pointer input  file                                    */
// int    ecode   ; /* element code number.                                   */
// int     numnp  ; /* number of nodal points.                                */
// int     bnodes ; /* number of boundary nodes.                              */
// int     xndof  ; /* maximum number of nodal degrees of freedom.            */
// int    *id     ; /* vector of dof status (0=free, 1=restrained)            */
// int    *idb    ; /* identification vector of boundary nodes [bsides[0]]    */
// double *bn_cs  ; /* boundary node cosine of local axes rotation angle      */
// double *bn_sn  ; /* boundary node sine   of local axes rotation angle      */

/****************************************************************************/



/*******************************************[ USER LIBRARY : MESHDAT.C ]******
 #                                                                           #
 #  FUNCTION :  sideinc ()                                                   #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #    Input of side nodal incidences.                                        #
 #                                                                           #
 #    A line of SIDE INCIDENCES data has the next syntax:                    #
 #                                                                           #
 #          NSIDE   N1  N2  N3  N4  ...  NXNSN                               #
 #                                                                           #
 #    where                                                                  #
 #          NSIDE : side number.                                             #
 #             Ni : node i of side NSIDE.                                    #
 #                                                                           #
 #    The number of nodes must be equal to the maximum number of nodes       #
 #    per side specified (parameter XNSN). The program does not accepts      #
 #    sides with different number of nodes.                                  #
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


  void sideinc (FILE *fpi, int numnp, int sides, int xnsn, int *side_nodes, int *side_vlev, int ecode);

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE *fpi       ; /* pointer input  file                                   */
// int   numnp     ; /* number of nodal points.                               */
// int   sides     ; /* number of sides.                                      */
// int   xnsn      ; /* maximum number of side nodes.                         */
// int  *side_nodes; /* nodes of each side                                    */
// int  *side_vlev ; /* level of side vertexs                                 */
// int    ecode    ; /* element code number.                                  */

/****************************************************************************/





/*******************************************[ USER LIBRARY : MESHDAT.C ]******
 #                                                                           #
 #  FUNCTION :  sidebc ()                                                    #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #    Input of side boundary conditions, tangential and normal directions.   #
 #    Input of coordinates of control points of boundary sides.              #
 #    Computation of cosine, sine of local angle of rotation anticlockwse.   #
 #    Local system x-direction is tangential and y_direction is normal.      #
 #                                                                           #
 #    A line of SIDE BOUNDARY CONDITIONS data has the next syntax:           #
 #                                                                           #
 #          NSIDE   R1  R2  R3  R4  ...  RXNDOF  XC YC                       #
 #                                                                           #
 #    where                                                                  #
 #          NSIDE : side number                                              #
 #             Ri : restraint condition of degree of freedom (dof) i         #
 #                  of side NSIDE.                                           #
 #             XC : X coordinate of control point.                           #
 #             YC : Y coordinate of control point.                           #
 #                                                                           #
 #    each restraint condition Ri must be an integer 0 or 1 which is         #
 #    interpreted as:                                                        #
 #                                                                           #
 #             Ri=1  : restrained.                                           #
 #             Ri=0  : free.                                                 #
 #                                                                           #
 #    The number of restraint conditions must be equal to the maximum        #
 #    number of nodal degrees of freedom specified (parameter XNDOF).        #
 #                                                                           #
 #    The program  accepts sides with different number of dof,  but          #
 #    internally all the sides are considered as having equal number         #
 #    of dofs. Then, for sides with less dofs than the maximum , the         #
 #    additional dofs MUST BE RESTRAINED.                                    #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H :                                                                #
 #    int  fscanf() : performs formatted input from a stream.                #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void sidebc (FILE *fpi, int bsides, int xnsn, int xnbdof, int nudim, int *side_nodes, int *side_bc,
               double *xyz, double *xyzc, double *bs_cs, double *bs_sn, double *bn_nlength, double *bn_slength, int numcurv);

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE *fpi       ; /* pointer input  file                                   */
// int   bsides    ; /* number of boundary sides.                             */
// int   xnsn      ; /* maximum number of side nodes.                         */
// int   xnbdof    ; /* maximum number of boundary degrees of freedom.        */
// int   nudim     ; /* number of spatial dimensions.                         */
// int  *side_nodes; /* nodes of each side                                    */
// int  *side_bc   ; /* side boundary conditions [bsides * xndof]             */
// double  *xyz    ; /* nodal coordinates [numnp * nudim].                    */
// double  *xyzc   ; /* coordinates of control points [numnp * nudim].        */
// double *bs_cs   ; /* boundary side cosine of local axes rotation angle     */
// double *bs_sn   ; /* boundary side sine   of local axes rotation angle     */
// double *bn_nlength; /* boundary node summ of adjacent lengths normal       */
// double *bn_slength; /* boundary node summ of adjacent lengths tangent      */
// int numcurv     ; /* number of curved boundary sides initial mesh          */

/****************************************************************************/








/*******************************************[ USER LIBRARY : MESHDAT.C ]******
 #                                                                           #
 #  FUNCTION :  matprop ()                                                   #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #    Input of material properties for static analysis of solids.            #
 #                                                                           #
 #    First is read the control data of this block of elements               #
 #    (not implemented)                                                      #
 #                                                                           #
 #    Then the properties for each material set of this block are read       #
 #    from the input file.                                                   #
 #                                                                           #
 #    A line of MATERIAL PROPERTIES data has the next syntax:                #
 #                                                                           #
 #          MATNB   PROP[1]  PROP[2] .... PROP[XNUMPR]                       #
 #                                                                           #
 #    where                                                                  #
 #          MATNB  : material number                                         #
 #          PROP[i]: propertie i for material MATNB.                         #
 #                                                                           #
 #    The next list contains the codes for the elements implemented          #
 #    and its numbers of properties for static analysis:                     #
 #                                                                           #
 #             numpr   : number of material properties                       #
 #             prop[1] : material property  1                                #
 #             prop[2] : material property  2                                #
 #                :                                                          #
 #                :                                                          #
 #             prop[numpr] : material property  numpr                        #
 #                                                                           #
 #                                                                           #
 #    Element codes for PLANE elements                                       #
 #                                                                           #
 #    ECODE  1   : LINEAR TRIANGLE VERTEX NODES - PLANE STRESS               #
 #                                                                           #
 #                 numpr   : 3                                               #
 #                 prop[1] : Young's modulus.                                #
 #                 prop[2] : Poisson's ratio.                                #
 #                 prop[3] : thickness.                                      #
 #                                                                           #
 #    ECODE  2   : QUADRATIC TRIANGLE RIGHT SIDES - PLANE STRESS             #
 #                                                                           #
 #                 numpr   : 3                                               #
 #                 prop[1] : Young's modulus.                                #
 #                 prop[2] : Poisson's ratio.                                #
 #                 prop[3] : thickness.                                      #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H :                                                                #
 #    int  fscanf() : performs formatted input from a stream.                #
 #                                                                           #
 #                                                                           #
 #***************************************************************************/


  void matprop (FILE *fpi, int ecode, int numat, double *prop, int xnumpr);


/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE     *fpi;   /* pointer input  file.                                   */
// int     ecode;   /* element code number.                                   */
// int     numat;   /* number of material sets.                               */
// double  *prop;   /* material properties [numat * xnumpr].                  */
// int    xnumpr;   /* maximum number of material properties.                 */

/****************************************************************************/






/*******************************************[ USER LIBRARY : MESHDAT.C ]******
 #                                                                           #
 #  FUNCTION :  elminc ()                                                    #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #    Input of material type and side incidences for each element.           #
 #                                                                           #
 #    All elements of the block must be of the same type.                    #
 #                                                                           #
 #    A line of ELEMENT PROPERTIES data has the next syntax:                 #
 #                                                                           #
 #          NELMT   MAT  N1  N2  N3  .....  Nxnen                            #
 #                                                                           #
 #    where                                                                  #
 #          NELMT : element number                                           #
 #            MAT : material code for element NELMT                          #
 #             Ni : number of the i-th side of the element NELMT             #
 #                                                                           #
 #    The next list contains the codes for the elements implemented          #
 #    and its numbers of sides.                                              #
 #                                                                           #
 #    Element codes for PLANE TRIANGULAR elements - SIDES = 3.               #
 #                                                                           #
 #    ECODE  1   : LINEAR TRIANGLE VERTEX NODES - PLANE STRESS               #
 #    ECODE  2   : QUADRATIC TRIANGLE RIGHT SIDES - PLANE STRESS             #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H :                                                                #
 #    int  fscanf() : performs formatted input from a stream.                #
 #                                                                           #
 #                                                                           #
 #***************************************************************************/


  void elminc (FILE *fpi, int ecode, int numel, int sides, int bsides, int *side_mat,
               int numat, int *mat, int *elm_sides, int xnes);

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE       *fpi;   /* pointer input  file.                                 */
// int       ecode;   /* element code number.                                 */
// int       numel;   /* number of elements.                                  */
// int       sides;   /* number of sides.                                     */
// int      bsides;   /* number of bsides.                                    */
// int   *side_mat;   /* side material index of adjacent element              */
// int       numat;   /* number of material sets.                             */
// int       *mat ;   /* element material type [numel].                       */
// int  *elm_sides;   /* incidences vector  [numel * xnes].                   */
// int        xnes;   /* maximum number of element sides.                     */

/*****************************************************************************/





/*******************************************[ USER LIBRARY : MESHDAT.C ]******
 #                                                                           #
 #  FUNCTION :  elmnod()                                                     #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #    Computation of element nodal incidences.                               #
 #    Computation of side opposite nodal incidences.                         #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H :                                                                #
 #    int  fscanf() : performs formatted input from a stream.                #
 #                                                                           #
 #                                                                           #
 #***************************************************************************/


  void elmnod (int ecode, int numel, int xnsn, int *side_nodes, int *elm_sides, int *elm_nodes, int *elm_vlev);

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int       ecode;   /* element code number.                                 */
// int       numel;   /* number of elements.                                  */
// int        xnsn;   /* maximum number of side nodes.                        */
// int *side_nodes;   /* side nodes                                           */
// int  *elm_sides;   /* incidences vector  [numel * xnes].                   */
// int  *elm_nodes;   /* element nodes      [numel * xnen]                    */
// int  *elm_vlev ;   /* element vertex levels [numel * xnev]                 */

/****************************************************************************/


#endif
