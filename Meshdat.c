

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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>


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


  void  meshdat (FILE *fpi, int numnp, int numel, int sides, int bsides,
                 int nudim, int ecode, double *xyz, double *xyzc, int *id,
                 int xndof, int xnbdof, int xnsn, int xnes, int *side_nodes,
                 int *side_vlev, int *elm_sides, int *elm_nodes, int *elm_vlev,
                 int *side_bc, int *side_mat, int numat, int *mat, double *prop,
                 int xnumpr, int *idb, int numcurv, double *bs_cs, double *bs_sn,
                 double *bn_cs, double *bn_sn, double *bn_nlength, double *bn_slength)

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


/*____________________________________________________________________________

								  BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

							PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/


  void coord   (FILE*, int, int, double*);
  void nodebc  (FILE* fpi, int ecode, int numnp, int bsides, int xndof, int *id,
                int *idb, double *bn_cs, double *bn_sn);
  void sideinc (FILE*, int, int, int, int*, int*, int);
  void sidebc  (FILE* fpi, int bsides, int xnsn, int xnbdof, int nudim,
                int *side_nodes,
                int *side_bc, double *xyz, double *xyzc, double *bs_cs,
                double *bs_sn, double *bn_nlength, double *bn_slength,
                int numcurv);

  void matprop (FILE*, int, int, double*, int);
  void elminc  (FILE*, int, int, int, int, int*, int, int*, int*, int);
  void elmnod  (int ecode, int numel, int xnsn, int *side_nodes,
	            int *elm_sides, int *elm_nodes, int *elm_vlev);

/*____________________________________________________________________________

								  NODAL COORDINATES
____________________________________________________________________________*/


  coord (fpi, numnp, nudim, xyz);

/*____________________________________________________________________________

							 NODAL BOUNDARY CONDITIONS
___________________________________________________________________________ */


  nodebc (fpi, ecode, numnp, bsides, xndof, id, idb, bn_cs, bn_sn);

/*____________________________________________________________________________

								  SIDE INCIDENCES
___________________________________________________________________________ */


  sideinc (fpi, numnp, sides, xnsn, side_nodes, side_vlev, ecode);


/*____________________________________________________________________________

				 SIDE BOUNDARY CONDITIONS AND CONTROL POINTS
___________________________________________________________________________ */


  sidebc (fpi, bsides, xnsn, xnbdof, nudim, side_nodes, side_bc, xyz, xyzc,
          bs_cs, bs_sn, bn_nlength, bn_slength, numcurv);


/*____________________________________________________________________________

								MATERIAL PROPERTIES
___________________________________________________________________________ */


  matprop (fpi, ecode, numat, prop, xnumpr);


/*____________________________________________________________________________

								ELEMENT SIDE INCIDENCES
___________________________________________________________________________ */


  elminc (fpi, ecode, numel, sides, bsides, side_mat, numat, mat, elm_sides, xnes);


/*____________________________________________________________________________

								ELEMENT NODAL INCIDENCES
___________________________________________________________________________ */


  elmnod (ecode, numel, xnsn, side_nodes, elm_sides, elm_nodes, elm_vlev);


/*____________________________________________________________________________

				  END OF FUNCTION
____________________________________________________________________________*/

} /* end of meshdat() */





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


  void coord (FILE *fpi, int numnp, int nudim, double *xyz)

/*____________________________________________________________________________

									PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE    *fpi ;   /* pointer input  file                                    */
// int     numnp;   /* number of nodal points.                                */
// int     nudim;   /* number of spatial dimensions.                          */
// double  *xyz ;   /* nodal coordinates.                                     */

/*____________________________________________________________________________

									BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

							  LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int    node     ;   /* node number  (0 = generation command line)          */
 double coor     ;   /* nodal coordinate value                              */
 int    j1       ;   /* auxiliar pointer for vector xyz[]                   */
 int    i, j     ;   /* loop counter                                        */

/*____________________________________________________________________________

									  LOOP ON NODES
____________________________________________________________________________*/


 for (i=1 ; i<=numnp ; i++)
	 {

/*____________________________________________________________________________

									 READ NODE NUMBER
____________________________________________________________________________*/


	  fscanf(fpi,"%d",&node);
	  if (node>numnp || node<0)
		 {
			 printf ("\n");
			 printf ("FATAL ERROR: ");
			 printf ("ILLEGAL NODE NUMBER IN NODAL COORDINATES");
			 printf ("\n");
			 exit(1);
		 }

/*____________________________________________________________________________

								 READ NODAL COORDINATES
____________________________________________________________________________*/

	  j1 = nudim * (node - 1);
	  for (j=1 ; j<=nudim ; j++)
		  {
			 fscanf(fpi,"%lf",&coor);
			 xyz[j1 +  j] = coor;
		  }

/*____________________________________________________________________________

								 END OF LOOP ON NODES
____________________________________________________________________________*/

	 }
/*____________________________________________________________________________

									  END OF FUNCTION
____________________________________________________________________________*/

}/* end of coord() */




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


  void nodebc (FILE *fpi, int ecode, int numnp, int bnodes, int xndof, int *id,
               int  *idb, double *bn_cs, double *bn_sn)

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



/*____________________________________________________________________________

								  BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

							  LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int    node     ;   /* boundary node number                                */
 double angle    ;   /* local axes angle in degrees                         */
 double pi       ;   /* number pi                                           */
 int    rdof     ;   /* dof restraint condition (0:free , 1:restrained)     */
 int    j1       ;   /* auxiliar pointer for vector id[]                    */
 int    i, j, i1 ;   /* loop counters                                       */

/*____________________________________________________________________________

								      NUMBER PI
____________________________________________________________________________*/

//  pi = M_PI;
  pi=3.14159;

/*____________________________________________________________________________

								 LOOP ON BOUNDARY NODES
____________________________________________________________________________*/

 i1 = ((ecode==4) || (ecode==5)) ? 2*bnodes : bnodes;

 if (ecode == 8) i1 = numnp;

 for (i=1 ; i<=i1 ; i++)
	{

/*____________________________________________________________________________

									READ BOUNDARY NODE NUMBER
____________________________________________________________________________*/


		fscanf(fpi,"%d",&node);
		if (node>numnp || node<0)
		  {
			 printf ("\n");
			 printf ("FATAL ERROR: ");
			 printf ("ILLEGAL NODE NUMBER IN NODAL BOUNDARY CONDITIONS");
			 printf ("\n");
			 exit(1);
		  }


/*____________________________________________________________________________

							READ NODAL ANGLE IN DEGREES
____________________________________________________________________________*/


  if ((ecode!=4)&&(ecode!=5))
  {
     idb[node] = i;

     fscanf(fpi,"%lf",&angle);
		if (angle>360. || angle<-360.0)
        {
			 printf ("\n");
			 printf ("FATAL ERROR: ");
			 printf ("OUT OF RANGE NODAL ANGLE [-360, 360]");
			 printf ("\n");
			 exit(1);
        }

      angle *= pi/180.;
      bn_cs[i] = cos(angle);
      bn_sn[i] = sin(angle);
  }

/*____________________________________________________________________________

							READ NODAL BOUNDARY CONDITIONS
____________________________________________________________________________*/


		j1 = xndof * (node - 1);
		for (j=1 ; j<=xndof ; j++)
			{
			  fscanf (fpi,"%d", &rdof);
			  if (rdof>1 || rdof<0)
				 {
					printf ("\n");
					printf ("FATAL ERROR: ");
					printf ("ILLEGAL RESTRAINT IN NODAL BOUNDARY CONDITIONS");
					printf ("\n");
					exit(1);
				 }

           id[++j1] = rdof;

			}

/*____________________________________________________________________________

							END OF LOOP ON BOUNDARY NODES
____________________________________________________________________________*/

	 }
/*____________________________________________________________________________

									 END OF FUNCTION
____________________________________________________________________________*/

}/* end of nodebc() */




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


  void sideinc (FILE *fpi, int numnp, int sides, int xnsn, int *side_nodes,
                int *side_vlev, int ecode)

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
/*____________________________________________________________________________

								  BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

							  LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int   side     ;    /* side number  (0 = generation command line)          */
 int   node     ;    /* node number                                         */
 int   i1       ;    /* auxiliar pointer for vector side_nodes[]            */
 int   j1       ;    /* auxiliar                                            */
 int   i, j     ;    /* loop counters                                       */

/*____________________________________________________________________________

									  LOOP ON SIDES
____________________________________________________________________________*/


 for (i=1 ; i<=sides ; i++)
	 {

/*____________________________________________________________________________

									READ SIDE NUMBER
____________________________________________________________________________*/


		fscanf(fpi,"%d",&side);
		if (side>sides || side<0)
		  {
			 printf ("\n");
			 printf ("FATAL ERROR: ");
			 printf ("ILLEGAL SIDE NUMBER IN SIDE INCIDENCES");
			 printf ("\n");
			 exit(3);
		  }

/*____________________________________________________________________________

								READ SIDE INCIDENCES
____________________________________________________________________________*/


		i1 = xnsn * (side - 1);
		j1 = xnsn;
    for (j=1 ; j<=j1 ; j++)
			{
			  fscanf (fpi,"%d", &node);
			  if (node>numnp)
				 {
					printf ("\n");
					printf ("FATAL ERROR: ");
					printf ("ILLEGAL NODE NUMBER IN SIDE INCIDENCES");
					printf ("\n");
					exit(1);
				 }

			  side_nodes[i1 + j] = node;
			}


/*____________________________________________________________________________

					SIDE VERTEX LEVELS FOR INITIAL LEVEL (only two vertexs)
____________________________________________________________________________*/

		i1 = 2 * (side - 1);
		side_vlev[++i1] = 0;
		side_vlev[++i1] = 0;

/*____________________________________________________________________________

								 END OF LOOP ON SIDES
____________________________________________________________________________*/

	 }

/*____________________________________________________________________________

									 END OF FUNCTION
____________________________________________________________________________*/

}/* end of sideinc() */






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


  void sidebc (FILE *fpi, int bsides, int xnsn, int xnbdof, int nudim,
               int *side_nodes, int *side_bc, double *xyz, double *xyzc,
               double *bs_cs, double *bs_sn, double *bn_nlength, double *bn_slength,
               int numcurv)


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

/*____________________________________________________________________________

								  BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

							  LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int    side     ;  /* side number                                          */
 int    n1, n2   ;  /* node numbers of side nodes                           */
 int    rdof     ;  /* dof restraint condition (0:free , 1:restrained)      */
 int    i1       ;  /* auxiliar pointer for vectors side_bc[], side_nodes[] */
 int    i, j     ;  /* loop counters                                        */
 double x1, y1   ;  /* coordinates of first node                            */
 double x2, y2   ;  /* coordinates of second node                           */
 double x21, y21 ;  /* side projections on x, y axes.                       */
 double L21      ;  /* side lenght.                                         */

/*____________________________________________________________________________

								 LOOP ON BOUNDARY SIDES
____________________________________________________________________________*/


 for (i=1 ; i<=bsides ; i++)
	 {

/*____________________________________________________________________________

									READ SIDE NUMBER
____________________________________________________________________________*/


		fscanf(fpi,"%d",&side);
		if (side>bsides || side<0)
		  {
			 printf ("\n");
			 printf ("FATAL ERROR: ");
			 printf ("ILLEGAL SIDE NUMBER IN SIDE BOUNDARY CONDITIONS");
			 printf ("\n");
			 exit(1);
		  }


/*____________________________________________________________________________

							READ SIDE BOUNDARY CONDITIONS
____________________________________________________________________________*/


		i1 = xnbdof * (side - 1);
		for (j=1 ; j<=xnbdof ; j++)
			{
			  fscanf (fpi,"%d", &rdof);
			  if (rdof>1 || rdof<0)
				 {
					printf ("\n");
					printf ("FATAL ERROR: ");
					printf ("ILLEGAL RESTRAINT IN SIDE BOUNDARY CONDITIONS");
					printf ("\n");
					exit(1);
				 }

           side_bc[++i1] = rdof;
			}

/*____________________________________________________________________________

							 SIDE COSINE, SINE. (ONLY STRAIGHT SIDES)
____________________________________________________________________________*/


      i1 = xnsn * (side - 1);
      n1 = side_nodes[++i1];
      n2 = side_nodes[++i1];

      x1 = xyz[nudim*(n1-1)+1];
      y1 = xyz[nudim*(n1-1)+2];
      x2 = xyz[nudim*(n2-1)+1];
      y2 = xyz[nudim*(n2-1)+2];

      x21 = x2 - x1;
      y21 = y2 - y1;

      L21 = sqrt(x21*x21+y21*y21);

      bs_cs[side] = x21/L21;
      bs_sn[side] = y21/L21;

      i1 = xnbdof * (side - 1);
      if (side_bc[++i1])
      {
        bn_slength[n1] += L21;
        bn_slength[n2] += L21;
      }

     if (side_bc[++i1])
     {
       bn_nlength[n1] += L21;
       bn_nlength[n2] += L21;
     }

/*____________________________________________________________________________

					 CONTROL POINT COORDINATES (ASSUMED STRAIGHT SIDES)
____________________________________________________________________________*/


	  i1 = nudim * (side - 1);
      xyzc[++i1] = (x1+x2)/2;
      xyzc[++i1] = (y1+y2)/2;

/*____________________________________________________________________________

							END OF LOOP ON BOUNDARY SIDES
____________________________________________________________________________*/

	 }

/*____________________________________________________________________________

						   READ CONTROL POINT COORDINATES CURVED SIDES
____________________________________________________________________________*/

 if (numcurv)
 {
    for (i=1 ; i<=numcurv ; i++)
    {

/*____________________________________________________________________________

									READ BOUNDARY SIDE NUMBER
____________________________________________________________________________*/


		fscanf(fpi,"%d",&side);
		if (side>bsides || side<0)
        {
		  printf ("\n");
		  printf ("FATAL ERROR: ");
		  printf ("ILLEGAL BOUNDARY SIDE NUMBER IN CONTROL POINT COORDINATES");
		  printf ("\n");
		  exit(1);
        }

/*____________________________________________________________________________

							READ CONTROL POINT COORDINATES
____________________________________________________________________________*/


		i1 = nudim * (side - 1);
		for (j=1 ; j<=nudim ; j++)
		{
		  fscanf (fpi,"%lf", &x1);
          xyzc[++i1] = x1;
		}

/*____________________________________________________________________________

							END OF LOOP ON CURVED BOUNDARY SIDES
____________________________________________________________________________*/


    }

/*____________________________________________________________________________

					END OF READING CONTROL POINT COORDINATES
____________________________________________________________________________*/

 }

/*____________________________________________________________________________

									 END OF FUNCTION
____________________________________________________________________________*/

}/* end of sidebc() */








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


  void matprop (FILE *fpi, int ecode, int numat, double *prop, int xnumpr)


/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE     *fpi;   /* pointer input  file.                                   */
// int     ecode;   /* element code number.                                   */
// int     numat;   /* number of material sets.                               */
// double  *prop;   /* material properties [numat * xnumpr].                  */
// int    xnumpr;   /* maximum number of material properties.                 */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int    matnb;   /* material number                                         */
 double prp  ;   /* material property                                       */
 int    j1   ;   /* auxiliar pointer for vector prop[]                      */
 int    i, j ;   /* loop counters                                           */

/*____________________________________________________________________________

		SET ELEMENT NUMBER OF MATERIAL PROPERTIES FOR EACH CODE
			 STATIC ANALYSIS
____________________________________________________________________________*/

 int numpr[9] = {0, 3, 3, 3, 3, 3, 3, 3, 3};

/*____________________________________________________________________________

		 READ ELEMENT BLOCK DATA (NOT IMPLEMENTED)
			 bknumel, bkmumat, bkecode
____________________________________________________________________________*/
/*____________________________________________________________________________

			READ MATERIAL PROPERTIES OF EACH MATERIAL SET
____________________________________________________________________________*/

 for (i=1 ; i<=numat ; i++)
	 {
		fscanf(fpi,"%d",&matnb);

		j1 = xnumpr * (matnb-1);
		for (j=1 ; j<=numpr[ecode] ; j++)
			{
			  fscanf(fpi,"%lf",&prp);
			  prop[++j1] = prp;
			}
	 }

/*____________________________________________________________________________

				  END OF FUNCTION
____________________________________________________________________________*/

}/* end of matprop() */






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


  void elminc (FILE *fpi, int ecode, int numel, int sides, int bsides, int *side_mat, int numat,
               int *mat, int *elm_sides, int xnes)

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

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int    elm  ;   /* element number                                          */
 int    side ;   /* side number                                             */
 int    matnb;   /* material number                                         */
 int    j1   ;   /* auxiliar pointer for vector prop[]                      */
 int    i, j ;   /* loop counters                                           */

/*____________________________________________________________________________

					 SET ELEMENT NUMBER OF SIDES FOR EACH CODE
____________________________________________________________________________*/


 int nes[9] = {0, 3, 3, 3, 3, 3, 3, 3, 3};


/*____________________________________________________________________________

		 READ ELEMENT BLOCK DATA (NOT IMPLEMENTED)
			 bknumel, bkmumat, bkecode
____________________________________________________________________________*/
/*____________________________________________________________________________

									LOOP ON ELEMENTS
____________________________________________________________________________*/


 for (i=1 ; i<=numel ; i++)
	 {

/*____________________________________________________________________________

								  READ ELEMENT NUMBER
____________________________________________________________________________*/


		fscanf(fpi,"%d",&elm);

		if (elm>numel || elm<0)
		  {
			 printf ("\n");
			 printf ("FATAL ERROR: ");
			 printf ("ILLEGAL ELEMENT NUMBER IN ELEMENT INCIDENCES");
			 printf ("\n");
			 exit(2);
		  }


/*____________________________________________________________________________

					 READ ELEMENT MATERIAL TYPE AND INCIDENCES
____________________________________________________________________________*/

		fscanf(fpi,"%d",&matnb);

		if (matnb<1 || matnb>numat)
		  { printf ("\n");
			 printf ("FATAL ERROR: ");
			 printf ("ILLEGAL ELEMENT MATERIAL NUMBER");
			 printf ("\n");
			 exit(2);
		  }

		mat[elm] = matnb;

		j1 = xnes * (elm - 1);
		for (j=1 ; j<=nes[ecode] ; j++)
			{
			  fscanf (fpi,"%d", &side);
			  if (abs(side)>sides || side==0)
				 {
					printf ("\n");
					printf ("FATAL ERROR: ");
					printf ("ILLEGAL SIDE NUMBER IN ELEMENT INCIDENCES");
					printf ("\n");
					exit(1);
				 }

				elm_sides[j1 + j] = side;

              // propiedades de los lados (verifica primero si es de contorno)

                if (abs(side)>bsides)
                {
                  if (side>0)
                    {
                        side_mat[2*side]= matnb;
                    }
                  else
                    {
                        side_mat[-2*side-1]= matnb;
                    }
                }
                else
                {
                  side_mat[2*abs(side)-1]= matnb;
                }

	        }

/*____________________________________________________________________________

								END OF LOOP ON ELEMENTS
____________________________________________________________________________*/

	 }
/*____________________________________________________________________________

				  END OF FUNCTION
____________________________________________________________________________*/

}/* end of elminc() */






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


  void elmnod (int ecode, int numel, int xnsn, int *side_nodes, int *elm_sides,
               int *elm_nodes, int *elm_vlev)

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

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int L1, L2, L3; /* element sides                                           */
 int n1, n2, n3; /* element nodes                                           */
 int n4, n5, n6; /* element nodes                                           */
 int    i1   ;   /* auxiliar pointer for vector elm_nodes[]                 */
 int    i    ;   /* loop counter                                            */

/*____________________________________________________________________________

									LOOP ON ELEMENTS
____________________________________________________________________________*/

	switch (ecode)
	 {
	   case 1: /* linear triangle */
	   case 8: /* linear triangle cubic Bézier */
        for (i=1 ; i<=numel ; i++)
         	 {
            //----------  ELEMENT SIDES -------------------------------------
          		L1 = elm_sides[3*i-2];
          		L2 = elm_sides[3*i-1];
          		L3 = elm_sides[3*i];
            //----------  ELEMENT VERTEXS -----------------------------------
          		n1 = (L1<0) ? side_nodes[xnsn*(-L1-1)+2] : side_nodes[xnsn*(L1-1)+1];
          		n2 = (L2<0) ? side_nodes[xnsn*(-L2-1)+2] : side_nodes[xnsn*(L2-1)+1];
          		n3 = (L3<0) ? side_nodes[xnsn*(-L3-1)+2] : side_nodes[xnsn*(L3-1)+1];
          		i1 = 3*(i-1);
          		elm_nodes[++i1] = n1;
          		elm_nodes[++i1] = n2;
          		elm_nodes[++i1] = n3;
            //----------  ELEMENT VERTEX LEVELS -----------------------------
          		i1 = 3*(i-1);
          		elm_vlev[++i1] = 0;
          		elm_vlev[++i1] = 0;
          		elm_vlev[++i1] = 0;
	          }
       break;

	  case 2: /* quadratic triangle */
        for (i=1 ; i<=numel ; i++)
         	 {
            //----------  ELEMENT SIDES -------------------------------------
          		L1 = elm_sides[3*i-2];
          		L2 = elm_sides[3*i-1];
          		L3 = elm_sides[3*i];
            //----------  ELEMENT VERTEXS -----------------------------------
          		n1 = (L1<0) ? side_nodes[xnsn*(-L1-1)+2] : side_nodes[xnsn*(L1-1)+1];
          		n2 = (L2<0) ? side_nodes[xnsn*(-L2-1)+2] : side_nodes[xnsn*(L2-1)+1];
          		n3 = (L3<0) ? side_nodes[xnsn*(-L3-1)+2] : side_nodes[xnsn*(L3-1)+1];
          		n4 = (L1<0) ? side_nodes[xnsn*(-L1-1)+3] : side_nodes[xnsn*(L1-1)+3];
          		n5 = (L2<0) ? side_nodes[xnsn*(-L2-1)+3] : side_nodes[xnsn*(L2-1)+3];
          		n6 = (L3<0) ? side_nodes[xnsn*(-L3-1)+3] : side_nodes[xnsn*(L3-1)+3];
          		i1 = 6*(i-1);
          		elm_nodes[++i1] = n1;
          		elm_nodes[++i1] = n2;
          		elm_nodes[++i1] = n3;
          		elm_nodes[++i1] = n4;
          		elm_nodes[++i1] = n5;
          		elm_nodes[++i1] = n6;
            //----------  ELEMENT VERTEX LEVELS -----------------------------
          		i1 = 6*(i-1);
          		elm_vlev[++i1] = 0;
         		elm_vlev[++i1] = 0;
          		elm_vlev[++i1] = 0;
          		elm_vlev[++i1] = 0;
          		elm_vlev[++i1] = 0;
          		elm_vlev[++i1] = 0;
	          }
       break;

	   case 3: /* linear mixed triangle */
        for (i=1 ; i<=numel ; i++)
         	 {
            //----------  ELEMENT SIDES -------------------------------------
          		L1 = elm_sides[3*i-2];
          		L2 = elm_sides[3*i-1];
          		L3 = elm_sides[3*i];
            //----------  ELEMENT VERTEXS -----------------------------------
          		n1 = (L1<0) ? side_nodes[xnsn*(-L1-1)+2] : side_nodes[xnsn*(L1-1)+1];
          		n2 = (L2<0) ? side_nodes[xnsn*(-L2-1)+2] : side_nodes[xnsn*(L2-1)+1];
          		n3 = (L3<0) ? side_nodes[xnsn*(-L3-1)+2] : side_nodes[xnsn*(L3-1)+1];
          		i1 = 3*(i-1);
          		elm_nodes[++i1] = n1;
          		elm_nodes[++i1] = n2;
          		elm_nodes[++i1] = n3;
            //----------  ELEMENT VERTEX LEVELS -----------------------------
          		i1 = 3*(i-1);
          		elm_vlev[++i1] = 0;
          		elm_vlev[++i1] = 0;
          		elm_vlev[++i1] = 0;
	          }
       break;

	  case 4: /* mixed MORLEY triangle */
        for (i=1 ; i<=numel ; i++)
         	 {
            //----------  ELEMENT SIDES -------------------------------------
          		L1 = elm_sides[3*i-2];
          		L2 = elm_sides[3*i-1];
          		L3 = elm_sides[3*i];
            //----------  ELEMENT VERTEXS -----------------------------------
          		n1 = (L1<0) ? side_nodes[xnsn*(-L1-1)+2] : side_nodes[xnsn*(L1-1)+1];
          		n2 = (L2<0) ? side_nodes[xnsn*(-L2-1)+2] : side_nodes[xnsn*(L2-1)+1];
          		n3 = (L3<0) ? side_nodes[xnsn*(-L3-1)+2] : side_nodes[xnsn*(L3-1)+1];
          		n4 = (L1<0) ? side_nodes[xnsn*(-L1-1)+4] : side_nodes[xnsn*(L1-1)+4];
          		n5 = (L2<0) ? side_nodes[xnsn*(-L2-1)+4] : side_nodes[xnsn*(L2-1)+4];
          		n6 = (L3<0) ? side_nodes[xnsn*(-L3-1)+4] : side_nodes[xnsn*(L3-1)+4];
          		i1 = 6*(i-1);
          		elm_nodes[++i1] = n1;
          		elm_nodes[++i1] = n2;
          		elm_nodes[++i1] = n3;
          		elm_nodes[++i1] = n4;
          		elm_nodes[++i1] = n5;
          		elm_nodes[++i1] = n6;
            //----------  ELEMENT VERTEX LEVELS -----------------------------
          		i1 = 6*(i-1);
          		elm_vlev[++i1] = 0;
         		elm_vlev[++i1] = 0;
          		elm_vlev[++i1] = 0;
          		elm_vlev[++i1] = 0;
          		elm_vlev[++i1] = 0;
          		elm_vlev[++i1] = 0;
	          }
       break;

	   case 5: /* linear plate mixed triangle */
	   case 7: // linear triangle by side
        for (i=1 ; i<=numel ; i++)
         	 {
            //----------  ELEMENT SIDES -------------------------------------
          		L1 = elm_sides[3*i-2];
          		L2 = elm_sides[3*i-1];
          		L3 = elm_sides[3*i];
            //----------  ELEMENT VERTEXS -----------------------------------
          		n1 = (L1<0) ? side_nodes[xnsn*(-L1-1)+2] : side_nodes[xnsn*(L1-1)+1];
          		n2 = (L2<0) ? side_nodes[xnsn*(-L2-1)+2] : side_nodes[xnsn*(L2-1)+1];
          		n3 = (L3<0) ? side_nodes[xnsn*(-L3-1)+2] : side_nodes[xnsn*(L3-1)+1];
          		i1 = 3*(i-1);
          		elm_nodes[++i1] = n1;
          		elm_nodes[++i1] = n2;
          		elm_nodes[++i1] = n3;
            //----------  ELEMENT VERTEX LEVELS -----------------------------
          		i1 = 3*(i-1);
          		elm_vlev[++i1] = 0;
          		elm_vlev[++i1] = 0;
          		elm_vlev[++i1] = 0;
	          }
       break;


      }
/*____________________________________________________________________________

				  END OF FUNCTION
____________________________________________________________________________*/

}/* end of elmnod() */



