/**********************[          Trimm2dmf.cpp            ]********************/
#include"Trimm2dmf.h"


/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: TR2DMLMF.C                                                 #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #        Routines for computation of midside interpolated forces.           #
 #        All elements of current mesh for all levels.                       #
 #        Linear triangular elements. Element code = 1.                      #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  trimm2dmf() : Computation of midside forces linear mixed triangle. #
 #                                                                           #
 #****************************************************************************/





/*******************************************[ USER LIBRARY : TRIMM2DMF.C ]****
 #                                                                           #
 #  FUNCTION :  t1m2dmf()                                                    #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of side stresses for linear mixed triangles.             #
 #      Linear variation for the stresses along the side is assumed.         #
 #      Average value is assigned for sides belonging to two elements.       #
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


  void t1m2dmf (int sides, int bsides, int ndof,
                 int nudim, int *side_nodes, int *side_bc,
                 int xnsn, double *xyz, int *side_mat,
                 int *mat, double *prop, int numpr,
                 double *u , double *Fn, double *Fs,
                 double *bs_Fn, double *bs_Fs)

/*____________________________________________________________________________

							 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int     sides;   /* number of sides.                                       */
// int    bsides;   /* number of boundary sides.                              */
// int     ndof ;   /* maximum number of nodal degrees of freedom.            */
// int     nudim;   /* number of spatial dimensions.                          */
// int *side_nodes; /* nodes of each side       [sides * xnsn]                */
// int    *side_bc; /* side boundary conditions [Bsides[0] * xndof]           */
// int   xnsn   ;   /* maximum number of side nodes.                          */
// double *xyz  ;   /* nodal coordinates [numnp * nudim].                     */
// int *side_mat;   /* side material index of adjacent elements [2*sides[0]]  */
// int    *mat  ;   /* element material type vector [numel].                  */
// double *prop ;   /* material properties vector [numat * numpr].            */
// int     numpr;   /* maximum number of material properties.                 */
// double *u    ;   /* expanded vector of nodal displacements [xndof*numnp]   */
// double  *Fn   ;   /* side normal tractions [sides[curlev]]                  */
// double  *Fs   ;   /* side tangential tractons [sides[curlev]]               */
// double *bs_Fn;   /* midside normal boundary tractions current level        */
// double *bs_Fs;   /* midside tangential boundary tractions current level    */

/*____________________________________________________________________________

								 BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

 void t1m2D_stress_bs (int bsides   , int xndof      , int nudim,
                       int *side_nodes, int *side_bc,
                       int xnsn     , double *xyz    , int *side_mat,
                       int *mat     , double *prop   , int xnumpr,
                       double *u    , double  *Fn     , double  *Fs     ,
                       double *bs_Fn, double *bs_Fs);

 void t1m2D_stress_is (int sides  ,
                       int bsides , int xndof      , int nudim,
                       int *side_nodes,
                       int xnsn   , double *xyz    , int *side_mat,
                       int *mat   , double *prop   , int xnumpr,
                       double *u  , double  *Fn     , double  *Fs);

/*__________________________________________________________________________

		             SET TO ZERO ALL FORCE COMPONENTS
___________________________________________________________________________*/


   memset (Fn, 0, (sides+1)*sizeof(double ));
   memset (Fs, 0, (sides+1)*sizeof(double ));

/*__________________________________________________________________________

				LOOP ON BOUNDARY SIDES
___________________________________________________________________________*/

  t1m2D_stress_bs (bsides, ndof, nudim, side_nodes, side_bc, xnsn, xyz,
                   side_mat, mat, prop, numpr, u, Fn, Fs, bs_Fn, bs_Fs);

/*__________________________________________________________________________

				LOOP ON INTERIOR SIDES
___________________________________________________________________________*/

  t1m2D_stress_is (sides, bsides, ndof, nudim, side_nodes,
                   xnsn, xyz, side_mat, mat, prop, numpr, u, Fn, Fs);

/*________________________________________________________________________

				             END OF FUNCTION
 ________________________________________________________________________*/

}/* end of t1m2dmf() */




/*******************************************[ USER LIBRARY : TRIMM2DMF.C ]****
 #                                                                           #
 #  FUNCTION :  t1m2D_stress_bs()                                            #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of side stresses for linear mixed triangles.             #
 #      Linear variation for the stresses along the side is assumed.         #
 #      Average value is assigned for sides belonging to two elements.       #
 #      ONLY FOR BOUNDARY SIDES.                                             #
 #                                                                           #
 #      Linear mixed triangular elements.                                    #
 #      Element code = 3.                                                    #
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

   void t1m2D_stress_bs (int bsides, int ndof, int nudim, int *side_nodes,
						  int *side_bc, int xnsn, double *xyz, int *side_mat,
						  int *mat, double *prop, int numpr, double *u, double *Fn,
						  double *Fs, double *bs_Fn, double *bs_Fs)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// int    bsides;   /* number of boundary sides.                              */
// int     ndof ;   /* maximum number of nodal degrees of freedom.            */
// int     nudim;   /* number of spatial dimensions.                          */
// int *side_nodes; /* nodes of each side       [sides * xnsn]                */
// int    *side_bc; /* side boundary conditions [Bsides[0] * xndof]           */
// int   xnsn   ;   /* maximum number of side nodes.                          */
// double *xyz  ;   /* nodal coordinates [numnp * nudim].                     */
// int *side_mat;   /* side material index of adjacent elements [2*sides[0]]  */
// int    *mat  ;   /* element material type vector [numel].                  */
// double *prop ;   /* material properties vector [numat * numpr].            */
// int     numpr;   /* maximum number of material properties.                 */
// double *u    ;   /* expanded vector of nodal displacements [xndof*numnp]   */
// double  *Fn   ;   /* side normal tractions [sides[curlev]]                  */
// double  *Fs   ;   /* side tangential tractons [sides[curlev]]               */
// double *bs_Fn;   /* midside normal boundary tractions current level        */
// double *bs_Fs;   /* midside tangential boundary tractions current level    */

/*____________________________________________________________________________

								 BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

							LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int n1, n2, n3;  /* element node numbers.                                  */
 double  darea;   /* two times the element area                             */
 double Sn, Tsn;  /* midside tractions                                      */
 int i         ;  /* loop counter                                           */
 int i1        ;  /* auxiliar pointers for mesh arrays.                     */

/*____________________________________________________________________________

		               DECLARATION OF MATERIAL PROPERTIES
____________________________________________________________________________*/

 double  E    ;   /* elasticity modulus                                     */
 double  nu   ;   /* Poisson's ratio                                        */
 double  h    ;   /* thickness                                              */

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

  void t1m2dSideStress (double *Sn , double *Tsn, double  *darea,
                        int n1, int n2, int n3,
                        double E, double nu, int nudim,
                        double *xyz, int xndof, double *u);

/*__________________________________________________________________________

				LOOP ON BOUNDARY SIDES
___________________________________________________________________________*/

 for (i=1; i<=bsides ; i++)
	 {
/*________________________________________________________________________

			NODAL INCIDENCES OF BOUNDARY SIDE
 __________________________________________________________________________*/

		i1 = xnsn*(i-1);
		n1 = side_nodes [++i1];
		n2 = side_nodes [++i1];
		n3 = side_nodes [++i1];

/*__________________________________________________________________________

				MATERIAL PROPERTIES OF ELEMENT
___________________________________________________________________________*/

		i1 = numpr * (side_mat[2*i-1] - 1);

		E  = prop[++i1];
		nu = prop[++i1];
		h  = prop[++i1];

/*__________________________________________________________________________

									 STRESSES ON SIDE L1 (N3-N2)
___________________________________________________________________________*/


  t1m2dSideStress (&Sn, &Tsn, &darea, n3, n2, n1, E, nu, nudim, xyz, ndof, u);

  i1 = ndof*(i-1);

  Fs[i] = (side_bc[++i1]) ? Tsn : bs_Fs[i];
  Fn[i] = (side_bc[++i1]) ? Sn  : bs_Fn[i];

/*________________________________________________________________________

		  END OF LOOP ON BOUNDARY SIDES
 ________________________________________________________________________*/

	 }

/*________________________________________________________________________

				             END OF FUNCTION
 ________________________________________________________________________*/

}/* end of t1m2D_stress_bs() */




/*******************************************[ USER LIBRARY : TRIMM2DMF.C ]****
 #                                                                           #
 #  FUNCTION :  t1m2D_stress_is()                                            #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of side stresses for linear mixed triangles.             #
 #      Linear variation for the stresses along the side is assumed.         #
 #      Average value is assigned for sides belonging to two elements.       #
 #      ONLY FOR INTERIROR SIDES.                                            #
 #                                                                           #
 #      Linear mixed triangular elements.                                    #
 #      Element code = 3.                                                    #
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

  void t1m2D_stress_is (int sides, int bsides, int ndof, int nudim, int *side_nodes,
                        int xnsn, double *xyz, int *side_mat, int *mat, double *prop,
						int numpr, double *u, double *Fn, double *Fs)


/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// int     sides;   /* number of sides.                                       */
// int    bsides;   /* number of boundary sides.                              */
// int     ndof ;   /* maximum number of nodal degrees of freedom.            */
// int     nudim;   /* number of spatial dimensions.                          */
// int *side_nodes; /* nodes of each side       [sides * xnsn]                */
// int   xnsn   ;   /* maximum number of side nodes.                          */
// double *xyz  ;   /* nodal coordinates [numnp * nudim].                     */
// int *side_mat;   /* side material index of adjacent elements [2*sides[0]]  */
// int    *mat  ;   /* element material type vector [numel].                  */
// double *prop ;   /* material properties vector [numat * numpr].            */
// int     numpr;   /* maximum number of material properties.                 */
// double *u    ;   /* expanded vector of nodal displacements [xndof*numnp]   */
// double  *Fn   ;   /* side normal tractions [sides[curlev]]                  */
// double  *Fs   ;   /* side tangential tractons [sides[curlev]]               */

/*____________________________________________________________________________

								 BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

							LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int n1, n2, n3, n4;  /* side node numbers.                                 */
 double  darea1;  /* two times the element 1 area                           */
 double  darea2;  /* two times the element 2 area                           */
 double Sn1, Tsn1;  /* midside tractions element 1                          */
 double Sn2, Tsn2;  /* midside tractions element 2                          */
 int i         ;  /* loop counter                                           */
 int i1        ;  /* auxiliar pointers for mesh arrays.                     */

/*____________________________________________________________________________

		               DECLARATION OF MATERIAL PROPERTIES
____________________________________________________________________________*/

 double  E    ;   /* elasticity modulus                                     */
 double  nu   ;   /* Poisson's ratio                                        */
 double  h    ;   /* thickness                                              */

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

  void t1m2dSideStress (double *Sn , double *Tsn, double  *darea,
                        int n1     , int n2     , int n3,
                        double   E , double   nu, int nudim,
                        double *xyz, int    xndof, double *u);

/*__________________________________________________________________________

				LOOP ON INTERIOR SIDES
___________________________________________________________________________*/

 for (i=bsides+1; i<=sides ; i++)
	 {
/*________________________________________________________________________

			NODAL INCIDENCES OF BOUNDARY SIDE
 __________________________________________________________________________*/

		i1 = xnsn*(i-1);
		n1 = side_nodes [++i1];
		n2 = side_nodes [++i1];
 	n3 = side_nodes [++i1];
  n4 = side_nodes [++i1];

/*__________________________________________________________________________

				MATERIAL PROPERTIES OF ELEMENT 1
___________________________________________________________________________*/

		i1 = numpr * (side_mat[2*i-1] - 1);

		E  = prop[++i1];
		nu = prop[++i1];
		h  = prop[++i1];

/*__________________________________________________________________________

									 STRESSES ON SIDE L1 (N1-N2) DUE TO ELEMENT 1
___________________________________________________________________________*/


  t1m2dSideStress (&Sn1, &Tsn1, &darea1, n3, n2, n1, E, nu, nudim, xyz, ndof, u);

/*__________________________________________________________________________

				MATERIAL PROPERTIES OF ELEMENT 2
___________________________________________________________________________*/

		i1 = numpr * (side_mat[2*i] - 1);

		E  = prop[++i1];
		nu = prop[++i1];
		h  = prop[++i1];

/*__________________________________________________________________________

									 STRESSES ON SIDE L1 (N1-N2) DUE TO ELEMENT 2
___________________________________________________________________________*/


  t1m2dSideStress (&Sn2, &Tsn2, &darea2, n4, n1, n2, E, nu, nudim, xyz, ndof, u);

  Fs[i] = (darea1*Tsn1 + darea2*Tsn2)/(darea1 + darea2);
  Fn[i] = (darea1*Sn1 + darea2*Sn2)/(darea1 + darea2);

/*________________________________________________________________________

		  END OF LOOP ON INTERIOR SIDES
 ________________________________________________________________________*/

	 }

/*________________________________________________________________________

				             END OF FUNCTION
 ________________________________________________________________________*/

}/* end of t1m2D_stress_is() */





/*******************************************[ USER LIBRARY : TRIMM2DMF.C ]****
 #                                                                           #
 #  FUNCTION :  t1m2dSideStress()                                            #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of side stresses for linear mixed triangles.             #
 #      Planar Strain Stress.                                                #
 #      Only one element contribution.                                       #
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


  void t1m2dSideStress (double *Sn, double *Tsn, double *darea, int n1, int n2, int n3,
                        double E, double nu, int nudim, double *xyz, int xndof, double *u)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// double *Sn,*Tsn; /* midside tractions                                      */
// double  *darea;  /* two times the element area = x23*y31 - x31*y23 		      */
// int  n1, n2, n3; /* node numbers of element                                */
// double   E   ;   /* elasticity modulus element                             */
// double   nu  ;   /* Poisson's ratio element                                */
// int     nudim;   /* number of spatial dimensions.                          */
// double *xyz  ;   /* nodal coordinates [numnp * nudim].                     */
// int    xndof ;   /* maximum number of nodal degrees of freedom.            */
// double *u    ;   /* expanded vector of nodal displacements [xndof*numnp]   */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int     i1   ;   /* auxiliar pointer for mesh arrays.                      */

/*____________________________________________________________________________

			  DECLARATION OF NODAL COORDINATES
____________________________________________________________________________*/

 double  x1   ;   /* X coordinate of node 1.                                */
 double  x2   ;   /* X coordinate of node 2.                                */
 double  x3   ;   /* X coordinate of node 3.                                */
 double  y1   ;   /* Y coordinate of node 1.                                */
 double  y2   ;   /* Y coordinate of node 2.                                */
 double  y3   ;   /* Y coordinate of node 3.                                */

/*____________________________________________________________________________

		DECLARATION OF COMPUTED GEOMETRIC PROPERTIES
____________________________________________________________________________*/

 double  x13, y31, x32, y23;
 double  L12;    /* square length of side 1                                 */
 double  c1, s1; /* director cosines of sides                               */

/*____________________________________________________________________________

   DECLARATION OF MATERIAL COEFFICIENTS
____________________________________________________________________________*/

 double  f1   ;   /* E/(1+nu)                                               */
 double  f2   ;   /* nu/(1-2nu)                                             */
/*____________________________________________________________________________

				 DECLARATION OF ELEMENT DISPLACEMENTS COEFFICIENTS
____________________________________________________________________________*/

 double u1, u2, u3;
 double v1, v2, v3;
 double u13, u23, v13, v23;

/*____________________________________________________________________________

  				     DECLARATION OF CENTROIDAL STRAINS AND STRESSES
____________________________________________________________________________*/

 double Exc, Eyc, Exyc;
 double Sxc, Syc, Sxyc;

/*____________________________________________________________________________

						  VERTEX	NODAL DISPLACEMENTS
____________________________________________________________________________*/

		i1 = xndof*(n1 - 1);
		u1 = u[++i1];
		v1 = u[++i1];

		i1 = xndof*(n2 - 1);
		u2 = u[++i1];
		v2 = u[++i1];

		i1 = xndof*(n3 - 1);
		u3 = u[++i1];
		v3 = u[++i1];

/*____________________________________________________________________________

						  DISPLACEMENTS COEFFICIENTS
____________________________________________________________________________*/

	   u13 = u1-u3;
	   u23 = u2-u3;

	   v13 = v1-v3;
	   v23 = v2-v3;

/*____________________________________________________________________________

				 NODAL COORDINATES
____________________________________________________________________________*/

		i1 = nudim*(n1 - 1);
		x1 = xyz[++i1];
		y1 = xyz[++i1];

		i1 = nudim*(n2 - 1);
		x2 = xyz[++i1];
		y2 = xyz[++i1];

		i1 = nudim*(n3 - 1);
		x3 = xyz[++i1];
		y3 = xyz[++i1];

/*__________________________________________________________________________

			 ELEMENT CONSTANTS
___________________________________________________________________________*/

		x32 = x3 - x2;
		y23 = y2 - y3;

		x13 = x1 - x3;
		y31 = y3 - y1;

  *darea = x13*y23 - x32*y31;

/*__________________________________________________________________________

					MATERIAL PROPERTIES
___________________________________________________________________________*/

   f1 = E/(1+nu);
   f2 = nu/(1-2*nu);

/*____________________________________________________________________________

						       CENTROIDAL STRAINS
____________________________________________________________________________*/

		 Exc  = (y23*u13 + y31*u23)/(*darea);
		 Eyc  = (x32*v13 + x13*v23)/(*darea);
   Exyc = (x32*u13 + x13*u23 + y23*v13 + y31*v23)/(*darea);

/*____________________________________________________________________________

						       CENTROIDAL STRESSES (PLANE STRAIN)
____________________________________________________________________________*/

		 Sxc  = f1*(Exc + f2*(Exc+Eyc));
		 Syc  = f1*(Eyc + f2*(Exc+Eyc));
   Sxyc = 0.5*f1*Exyc;

/*____________________________________________________________________________

                   STRESSES ON SIDE L1
____________________________________________________________________________*/

  L12 = x32*x32 + y23*y23;

  c1 = -x32;
  s1 = y23;

  *Sn = (Sxc*s1*s1 + Syc*c1*c1 - 2*Sxyc*s1*c1)/L12;
  *Tsn = ((Syc-Sxc)*s1*c1 + Sxyc*(c1*c1-s1*s1))/L12;


/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of t1m2dSideStress() */





