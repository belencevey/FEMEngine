
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
 #  void  tr2dmlmf()  : Multilevel computation of midside forces             #
 #                      linear triangle.                                     #
 #                                                                           #
 #****************************************************************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
//#include <mem.h>

/*******************************************[ USER LIBRARY : TR2DMLMF.C ]*****
 #                                                                           #
 #  FUNCTION :  tr1_nodal_avg()                                              #
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

void tr1_nodal_avg (int curlev      , int *numel      , int *nodes,
                      int *FRnumel    , int **velm      , int *BFRsides,
                      int *PRsides    , int **eside     , int **bside,
                      int **tside     , int xndof       , double **xyz,
                      int **elm_nodes , int **elm_vlev  , int **parent_elm,
                      int **side_nodes, int ** side_vlev, double **uold,
                      int *mat        , double *prop    , int xnumpr,
                      double  **Sx     , double **Sy      , double **Sxy,
                      int *nse, int nudim)


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

/*____________________________________________________________________________

								 BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

							LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int ie        ;   /* element counter for loop on elements.                 */
 int ie1       ;   /* initial element for loop on elements.                 */
 int level     ;   /* level   counter for loop on levels.                   */
 int n1, n2, n3;   /* element node numbers.                                 */
 int lv1, lv2, lv3;/* element node levels.                                  */
 int pelm      ;   /* parent element                                        */
 int side      ;   /* side number                                           */
 int i1        ;   /* auxiliar pointers for mesh arrays.                    */
 int i         ;   /* loop counter                                          */

/*____________________________________________________________________________

				 DECLARATION OF ELEMENT COORDINATES COEFFICIENTS
____________________________________________________________________________*/

 double x1, x2, x3;
 double y1, y2, y3;
 double x13, x32, y31, y23;
 double  area;   /* two times the element area = x13*y23-x32*y31           */

/*____________________________________________________________________________

				 DECLARATION OF ELEMENT DISPLACEMENTS COEFFICIENTS
____________________________________________________________________________*/

 double u1, u2, u3;
 double v1, v2, v3;
 double u13, u23, v13, v23;

/*____________________________________________________________________________

				 DECLARATION OF ELEMENT DISPLACEMENTS DERIVATIVES
____________________________________________________________________________*/

 double ux, uy, vx, vy;

/*____________________________________________________________________________

		               DECLARATION OF MATERIAL PROPERTIES
____________________________________________________________________________*/

 double  E    ;   /* elasticity modulus                                     */
 double  nu   ;   /* Poisson's ratio                                        */

 double  f1   ;   /* E/(1-nu*nu)                                            */
 double  f2   ;   /* E/2/(1+nu)                                             */

/*____________________________________________________________________________

  				     DECLARATION OF CENTROIDAL STRAINS AND STRESSES
____________________________________________________________________________*/

 double Exc, Eyc, Exyc;
 double Sxc, Syc, Sxyc;

/*____________________________________________________________________________

  				     DECLARATION OF SIDE STRESSES
____________________________________________________________________________*/

 double Sx1, Sy1, Sxy1;
 double Sx2, Sy2, Sxy2;

/*____________________________________________________________________________

  				     DECLARATION OF SIDE TRACTIONS
____________________________________________________________________________*/

/*__________________________________________________________________________

		             SET TO ZERO ALL NODAL COMPONENTS
___________________________________________________________________________*/

    for (level=0; level<=curlev ; level++)
       {
         memset (Sx[level] , 0, (nodes[level]+1)*sizeof(double ));
         memset (Sy[level] , 0, (nodes[level]+1)*sizeof(double ));
         memset (Sxy[level], 0, (nodes[level]+1)*sizeof(double ));
       }


/*__________________________________________________________________________

			             LOOP ON ALL LEVELS
___________________________________________________________________________*/

 for (level=0; level<=curlev ; level++)
	 {

     ie1 = (level<curlev) ? FRnumel[level]+1 : 1;

/*__________________________________________________________________________

				  LOOP ON ALL UNREFINED ELEMENTS OF EACH PREVIOUS LEVEL
                    AND ON ALL ELEMENTS OF CURRENT LEVEL
___________________________________________________________________________*/

 for (ie=ie1; ie<=numel[level] ; ie++)
	 {

/*________________________________________________________________________

							  NODAL INCIDENCES
 __________________________________________________________________________*/

 	  i1 = 3*(velm[level][ie]-1);
 	  n1 = elm_nodes[level][++i1];
	   n2 = elm_nodes[level][++i1];
	   n3 = elm_nodes[level][++i1];

/*___________________________________________________________________________

								 VERTEX NODAL LEVELS
 __________________________________________________________________________*/


	   i1 = 3*(velm[level][ie]-1);
	   lv1= elm_vlev[level][++i1];
	   lv2= elm_vlev[level][++i1];
	   lv3= elm_vlev[level][++i1];

/*____________________________________________________________________________

              						  VERTEX	NODAL DISPLACEMENTS
____________________________________________________________________________*/


    i1 = xndof*(n1 - 1);
	   u1 = uold[lv1][++i1];
	   v1 = uold[lv1][++i1];

	   i1 = xndof*(n2 - 1);
	   u2 = uold[lv2][++i1];
    v2 = uold[lv2][++i1];

	   i1 = xndof*(n3 - 1);
    u3 = uold[lv3][++i1];
	   v3 = uold[lv3][++i1];

/*____________________________________________________________________________

						                 DISPLACEMENTS COEFFICIENTS
____________________________________________________________________________*/

	   u13 = u1-u3;
	   u23 = u2-u3;

	   v13 = v1-v3;
	   v23 = v2-v3;

/*____________________________________________________________________________

							               	VERTEX NODAL COORDINATES
____________________________________________________________________________*/


  	i1 = nudim*(n1 - 1);
	  x1 = xyz[lv1][++i1];
   y1 = xyz[lv1][++i1];

   i1 = nudim*(n2 - 1);
   x2 = xyz[lv2][++i1];
   y2 = xyz[lv2][++i1];

  	i1 = nudim*(n3 - 1);
  	x3 = xyz[lv3][++i1];
 	 y3 = xyz[lv3][++i1];

/*____________________________________________________________________________

						                COORDINATES COEFFICIENTS
____________________________________________________________________________*/


	  x13 = x1-x3;
	  x32 = x3-x2;

   y31 = y3-y1;
	  y23 = y2-y3;

   area = x13*y23 - x32*y31;

/*__________________________________________________________________________

         					MATERIAL PROPERTIES AND STIFFNESS COEFICIENTS
___________________________________________________________________________*/

		 pelm = parent_elm[level][3*velm[level][ie]-2];

		 i1 = xnumpr * (mat[pelm] - 1);
		 E  = prop[++i1];
		 nu = prop[++i1];

   f1 = E/(1-nu*nu);
   f2 = 0.5*E/(1+nu);

/*____________________________________________________________________________

               ELEMENT DISPLACEMENTS DERIVATIVES
____________________________________________________________________________*/


   ux = (y23*u13+y31*u23)/area;
   uy = (x32*u13+x13*u23)/area;
   vx = (y23*v13+y31*v23)/area;
   vy = (x32*v13+x13*v23)/area;

/*____________________________________________________________________________

						       CENTROIDAL STRAINS
____________________________________________________________________________*/

   Exc  = ux;
   Eyc  = vy;
   Exyc = uy+vx;

/*____________________________________________________________________________

						                     CENTROIDAL STRESSES
____________________________________________________________________________*/

   Sxc  = f1*(Exc+nu*Eyc);
   Syc  = f1*(nu*Exc+Eyc);
   Sxyc = f2*Exyc;

/*____________________________________________________________________________

						                    STRESSES ON NODE N1
____________________________________________________________________________*/


   Sx[lv1][n1]  += Sxc;
   Sy[lv1][n1]  += Syc;
   Sxy[lv1][n1] += Sxyc;

/*____________________________________________________________________________

						                   STRESSES ON NODE N2
____________________________________________________________________________*/

   Sx[lv2][n2]  += Sxc;
   Sy[lv2][n2]  += Syc;
   Sxy[lv2][n2] += Sxyc;

/*____________________________________________________________________________

						                   STRESSES ON NODE N3
____________________________________________________________________________*/

   Sx[lv3][n3]  += Sxc;
   Sy[lv3][n3]  += Syc;
   Sxy[lv3][n3] += Sxyc;

/*________________________________________________________________________

	        END OF LOOP ON UNREFINED ELEMENTS OF PREVIOUS LEVELS
                 AND ON ALL ELEMENTS OF CURRENT LEVEL
 ________________________________________________________________________*/

	 }

/*________________________________________________________________________

                     END OF LOOP ON LEVELS
 ________________________________________________________________________*/

 }

/*_________________________________________________________________________

          LOOP ON ALL NODES OF COARSE MESH
 __________________________________________________________________________*/

      for (i=1 ; i<=nodes[0] ; i++)
         {
           Sx[0][i] /= nse[i];
           Sy[0][i] /= nse[i];
           Sxy[0][i] /= nse[i];
         }


/*_________________________________________________________________________

              LOOP ON ALL INTERIOR NODES OF EACH LEVEL
            (Since it is not possible to identify the interior
             nodes, we loop first an all nodes, and then the
             boundary nodes are corrected)

               LOOP ON ALL BOUNDARY NODES OF EACH LEVEL
 __________________________________________________________________________*/

 for (level=1; level<=curlev ; level++)
	 {
      for (i=1 ; i<=nodes[level] ; i++)
         {
           Sx[level][i] /= 6.;
           Sy[level][i] /= 6.;
           Sxy[level][i] /= 6.;
         }

      for (i=1 ; i<=BFRsides[level-1] ; i++)
         {
           side = eside[level-1][bside[level-1][i]];
           Sx[level][side] *= 2.;
           Sy[level][side] *= 2.;
           Sxy[level][side] *= 2.;
         }
 }


/*__________________________________________________________________________

	            LOOP ON ALL TRANSITION SIDES OF PREVIOUS LEVELS

              (AVERAGE STRESS FROM NODES OF PARENT SIDE)
___________________________________________________________________________*/


for (level=0 ; level<curlev ; level++)
   {
     if (PRsides[level])
       {
        for (i=1; i<=PRsides[level]; i++)
           {

           	 i1 = 2*(tside[level][i]-1);
           	 n1 = side_nodes[level][++i1];
	            n2 = side_nodes[level][++i1];

           	 i1 = 2*(tside[level][i]-1);
		           lv1= side_vlev[level][++i1];
		           lv2= side_vlev[level][++i1];

             Sx1 = Sx[lv1][n1];
             Sy1 = Sy[lv1][n1];
             Sxy1 = Sxy[lv1][n1];

             Sx2 = Sx[lv2][n2];
             Sy2 = Sy[lv2][n2];
             Sxy2 = Sxy[lv2][n2];

             side = eside[level][tside[level][i]];
             Sx [level+1][side] = (Sx1+Sx2)/2.;
             Sy [level+1][side] = (Sy1+Sy2)/2.;
             Sxy[level+1][side] = (Sxy1+Sxy2)/2.;

	          }
       }
   }

/*________________________________________________________________________

        				             END OF FUNCTION
 ________________________________________________________________________*/

}/* end of tr1_nodal_avg() */








/*******************************************[ USER LIBRARY : TR2DMLMF.C ]*****
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
 void tr1_flux_avg (int curlev    , int *numel      , int *FRnumel,
                      int **velm      , int *sides       , int *Bsides,
                      int *BFRsides   , int *FRsides     , int *PRsides,
                      int **vside     , int **bside      , int **tside,
                      int **ssside    , int xndof        , double **xyz, int nudim,
                      int **elm_sides , int **elm_nodes  , int **elm_vlev,
                      int **parent_elm, int **parent_side, int *side_bc,
                      double **uold   , int *mat         , double *prop,
                      int xnumpr      , double  **Fn      , double  **Fs,
                      double **bs_Fn  , double **bs_Fs,
                      double **fi     ,
                      double **bn_nlength, double **bn_slength,
                      double **bn_qn, double **bn_qs,
                      double **bs_cs, double **bs_sn,
                      int *nodes, int **side_nodes, int **side_vlev)


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

/*____________________________________________________________________________

								 BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

							LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int ie        ;   /* element counter for loop on elements.                 */
 int ie1       ;   /* initial element for loop on elements.                 */
 int level     ;   /* level   counter for loop on levels.                   */
 int L1, L2, L3;   /* element side numbers.                                 */
 int n1, n2, n3;   /* element node numbers.                                 */
 int lv1, lv2, lv3;/* element node levels.                                  */
 int pelm      ;   /* parent element                                        */
 int pside     ;   /* parent side                                           */
 int side      ;   /* side number                                           */
 int side1     ;   /* side number                                           */
 int side2     ;   /* side number                                           */
 int i1,i2,i3  ;   /* auxiliar pointers for mesh arrays.                    */
 int i         ;   /* loop counter                                          */

/*____________________________________________________________________________

				 DECLARATION OF ELEMENT COORDINATES COEFFICIENTS
____________________________________________________________________________*/

 double x1, x2, x3;
 double y1, y2, y3;
 double x13, x32, y31, y23, x12, y12;
 double  area;   /* two times the element area = x13*y23-x32*y31           */
 double L12, L23, L31; /* side lengths                                     */
 double cos1, sin1; /* cosine and sine of side 1                           */
 double cos2, sin2; /* cosine and sine of side 2                           */
 double cos3, sin3; /* cosine and sine of side 3                           */

/*____________________________________________________________________________

				 DECLARATION OF ELEMENT DISPLACEMENTS COEFFICIENTS
____________________________________________________________________________*/

 double u1, u2, u3;
 double v1, v2, v3;
 double u13, u23, v13, v23;
/*____________________________________________________________________________

		               DECLARATION OF MATERIAL PROPERTIES
____________________________________________________________________________*/

 double  E    ;   /* elasticity modulus                                     */
 double  nu   ;   /* Poisson's ratio                                        */
 double  th   ;   /* thickness                                              */
 double  f1   ;   /* E/(1-nu*nu)                                            */
 double  f2   ;   /* E/2/(1+nu)                                             */

/*____________________________________________________________________________

  				     DECLARATION OF CENTROIDAL STRAINS AND STRESSES
____________________________________________________________________________*/

 double Exc, Eyc, Exyc;
 double Sxc, Syc, Sxyc;

/*____________________________________________________________________________

  				     DECLARATION OF SIDE TRACTIONS
____________________________________________________________________________*/

 double Fnj, Fsj, F1, F2, F3, F4;
/*____________________________________________________________________________

  				     DECLARATION OF SIDE AND NODE REACTIONS
____________________________________________________________________________*/

 double Rn1, Rs1, Rn2, Rs2, Fx1, Fx2, Fy1, Fy2;

 // eliminar nudim = 2;

/*__________________________________________________________________________

		             SET TO ZERO ALL FORCE COMPONENTS
___________________________________________________________________________*/


    for (level=0; level<=curlev ; level++)
       {
         memset (Fn[level], 0, (sides[level]+1)*sizeof(double ));
         memset (Fs[level], 0, (sides[level]+1)*sizeof(double ));
       }



/*__________________________________________________________________________

			             LOOP ON ALL LEVELS
___________________________________________________________________________*/

 for (level=0; level<=curlev ; level++)
	 {

     ie1 = (level<curlev) ? FRnumel[level]+1 : 1;

/*__________________________________________________________________________

				  LOOP ON ALL UNREFINED ELEMENTS OF EACH PREVIOUS LEVEL
                    AND ON ALL ELEMENTS OF CURRENT LEVEL
___________________________________________________________________________*/

 for (ie=ie1; ie<=numel[level] ; ie++)
	 {

/*________________________________________________________________________

							  NODAL INCIDENCES
 __________________________________________________________________________*/

		i1 = 3*(velm[level][ie]-1);
		L1 = elm_sides[level][++i1];
		L2 = elm_sides[level][++i1];
		L3 = elm_sides[level][++i1];

		i1 = 3*(velm[level][ie]-1);
		n1 = elm_nodes[level][++i1];
		n2 = elm_nodes[level][++i1];
		n3 = elm_nodes[level][++i1];

/*___________________________________________________________________________

								 VERTEX NODAL LEVELS
 __________________________________________________________________________*/


		i1 = 3*(velm[level][ie]-1);
		lv1= elm_vlev[level][++i1];
		lv2= elm_vlev[level][++i1];
		lv3= elm_vlev[level][++i1];

/*____________________________________________________________________________

						  VERTEX	NODAL DISPLACEMENTS
____________________________________________________________________________*/


		i1 = xndof*(n1 - 1);
		u1 = uold[lv1][++i1];
		v1 = uold[lv1][++i1];

		i1 = xndof*(n2 - 1);
		u2 = uold[lv2][++i1];
		v2 = uold[lv2][++i1];

		i1 = xndof*(n3 - 1);
		u3 = uold[lv3][++i1];
		v3 = uold[lv3][++i1];

/*____________________________________________________________________________

						  DISPLACEMENTS COEFFICIENTS
____________________________________________________________________________*/

	   u13 = u1-u3;
	   u23 = u2-u3;

	   v13 = v1-v3;
	   v23 = v2-v3;

/*____________________________________________________________________________

								VERTEX NODAL COORDINATES
____________________________________________________________________________*/


		 i1 = nudim*(n1 - 1);
		 x1 = xyz[lv1][++i1];
		 y1 = xyz[lv1][++i1];

		 i1 = nudim*(n2 - 1);
		 x2 = xyz[lv2][++i1];
		 y2 = xyz[lv2][++i1];

		 i1 = nudim*(n3 - 1);
		 x3 = xyz[lv3][++i1];
   y3 = xyz[lv3][++i1];

/*____________________________________________________________________________

						  COORDINATES COEFFICIENTS
____________________________________________________________________________*/


	  x13 = x1-x3;
	  x32 = x3-x2;

   y31 = y3-y1;
	  y23 = y2-y3;

	  x12 = x1-x2;
	  y12 = y1-y2;

   area = x13*y23 - x32*y31;

/*__________________________________________________________________________

					MATERIAL PROPERTIES AND STIFFNESS COEFICIENTS
___________________________________________________________________________*/

		 pelm = parent_elm[level][3*velm[level][ie]-2];

		 i1 = xnumpr * (mat[pelm] - 1);
		 E  = prop[++i1];
		 nu = prop[++i1];
   th = prop[++i1];

   f1 = E/(1-nu*nu);
   f2 = 0.5*E/(1+nu);

/*____________________________________________________________________________

						       CENTROIDAL STRAINS
____________________________________________________________________________*/

		 Exc  = (y23*u13+y31*u23)/area;
		 Eyc  = (x32*v13+x13*v23)/area;
   Exyc = (x32*u13+x13*u23+y23*v13+y31*v23)/area;

/*____________________________________________________________________________

						       CENTROIDAL STRESSES
____________________________________________________________________________*/

		 Sxc  = f1*(Exc+nu*Eyc);
		 Syc  = f1*(nu*Exc+Eyc);
   Sxyc = f2*Exyc;

/*____________________________________________________________________________

                   UNIT LENGTH FORCES ON SIDE L1
____________________________________________________________________________*/


    L12 = x12*x12 + y12*y12;

    cos1 = x12;
    sin1 = y12;

    Fnj = (Sxc*sin1*sin1+Syc*cos1*cos1-2*Sxyc*sin1*cos1)*th/L12;
    Fsj = ((Syc-Sxc)*sin1*cos1+Sxyc*(cos1*cos1-sin1*sin1))*th/L12;

		  i1 = abs(L1);
		  Fn[level][i1] += Fnj;
		  Fs[level][i1] += Fsj;

/*____________________________________________________________________________

                   UNIT LENGTH FORCES ON SIDE L2
____________________________________________________________________________*/

    L23 = x32*x32 + y23*y23;
    cos2 = -x32;
    sin2 =  y23;

    Fnj = (Sxc*sin2*sin2+Syc*cos2*cos2-2*Sxyc*sin2*cos2)*th/L23;
    Fsj = ((Syc-Sxc)*sin2*cos2+Sxyc*(cos2*cos2-sin2*sin2))*th/L23;

		  i1 = abs(L2);
		  Fn[level][i1] += Fnj;
		  Fs[level][i1] += Fsj;

/*____________________________________________________________________________

                   UNIT LENGTH FORCES ON SIDE L3
____________________________________________________________________________*/

    L31 = x13*x13 + y31*y31;
    cos3 = -x13;
    sin3 =  y31;

    Fnj = (Sxc*sin3*sin3+Syc*cos3*cos3-2*Sxyc*sin3*cos3)*th/L31;
    Fsj = ((Syc-Sxc)*sin3*cos3+Sxyc*(cos3*cos3-sin3*sin3))*th/L31;

	  	i1 = abs(L3);
		  Fn[level][i1] += Fnj;
		  Fs[level][i1] += Fsj;

/*________________________________________________________________________

	        END OF LOOP ON UNREFINED ELEMENTS OF PREVIOUS LEVELS
                 AND ON ALL ELEMENTS OF CURRENT LEVEL
 ________________________________________________________________________*/

	 }

/*________________________________________________________________________

                     END OF LOOP ON PREVIOUS LEVELS
 ________________________________________________________________________*/

 }


/*_________________________________________________________________________

          LOOP ON ALL UNREFINED INTERIOR SIDES OF PREVIOUS LEVELS
            (Since it is not possible to distinguish beetween
             interior and boundary sides, we loop first an all
             sides, and then the boundary sides are corrected)

          LOOP ON ALL UNREFINED BOUNDARY SIDES OF PREVIOUS LEVEL
                   BOUNDARY TRACTIONS ON FREE SIDES
 __________________________________________________________________________*/

 for (level=0; level<curlev ; level++)
	 {
    i1 = FRsides[level]+1;
    for (i=i1 ; i<=sides[level] ; i++)
       {
         side = vside[level][i];
         Fn[level][side] /= 2.;
         Fs[level][side] /= 2.;
       }

    i1 = BFRsides[level]+1;
    for (i=i1 ; i<=Bsides[level] ; i++)
       {
         side = bside[level][i];
         Fn[level][side] *= 2.;
         Fs[level][side] *= 2.;
       }


    i1 = BFRsides[level]+1;
    for (i=i1 ; i<=Bsides[level] ; i++)
       {
         pside = parent_side[level][2*i-1];
         i3 = xndof*(pside-1);
         side = bside[level][i];
         if (!side_bc[++i3]) Fs[level][side] = bs_Fs[level][i];
			if (!side_bc[++i3]) Fn[level][side] = bs_Fn[level][i];
       }

 }



/*__________________________________________________________________________

	            LOOP ON ALL INTERIOR SIDES OF CURRENT FINE LEVEL
                  (First loop on all sides and then
                    boundary sides are corrected)
___________________________________________________________________________*/


      for (i=1 ; i<=sides[curlev] ; i++)
         {
           Fn[curlev][i] /= 2.;
           Fs[curlev][i] /= 2.;
         }


      for (i=1 ; i<=Bsides[curlev] ; i++)
         {
           side = bside[curlev][i];
           Fn[curlev][side] *= 2.;
           Fs[curlev][side] *= 2.;
         }

/*__________________________________________________________________________

	            LOOP ON ALL BOUNDARY SIDES OF CURRENT LEVEL
                    BOUNDARY TRACTIONS ON FREE SIDES
___________________________________________________________________________*/


 for (i=1; i<=Bsides[curlev] ; i++)
    {
      pside = parent_side[curlev][2*i-1];
      i1 = xndof*(pside-1);
      side = bside[curlev][i];
    	if (!side_bc[++i1]) Fs[curlev][side] = bs_Fs[curlev][i];
   		if (!side_bc[++i1]) Fn[curlev][side] = bs_Fn[curlev][i];
    }


/*__________________________________________________________________________

       CORRECTION OF AVERAGED FORCES (ONLY ONE ELEMENT CONTRIBUTION)

                 ON SON SIDES OF TRANSITION SIDES

       (NOTE: TRANSITION SIDES ARE REFINED SIDES THEY AREN'T DIVIDED BY 2)
___________________________________________________________________________*/


 for (level=0; level<curlev ; level++)
	 {
       i2 = PRsides[level];
       for (i=1; i<=i2 ; i++)
         {
           side  = tside[level][i];
           side1 = ssside[level][side];
           side2 = side1+1;

           F1 = Fn[level+1][side1];
           F2 = Fn[level+1][side2];
           F3 = Fn[level][side];
           F4 = (F1+F2+F3)/2.;

           Fn[level][side] = F4;
           Fn[level+1][side1] = (5*F1+2*F3-F2)/4.;
           Fn[level+1][side2] = (5*F2+2*F3-F1)/4.;

           F1 = Fs[level+1][side1];
           F2 = Fs[level+1][side2];
           F3 = Fs[level][side];
           F4 = (F1+F2+F3)/2.;

           Fs[level][side] = F4;
           Fs[level+1][side1] = (5*F1+2*F3-F2)/4.;
           Fs[level+1][side2] = (5*F2+2*F3-F1)/4.;
         }
    }


/*
 for (level=0; level<curlev ; level++)
	 {
       i2 = PRsides[level];
       for (i=1; i<=i2 ; i++)
         {
           side  = tside[level][i];
           side1 = ssside[level][side];
           side2 = side1+1;

           Fn[level+1][side1] *= 2;
           Fn[level+1][side2] *= 2;

           Fs[level+1][side1] *= 2;
           Fs[level+1][side2] *= 2;
         }
    }
*/


/*_________________________________________________________________________

          LOOP ON ALL UNREFINED BOUNDARY SIDES OF PREVIOUS LEVEL
       COMPUTATION OF NODAL AVERAGED TANGENTIAL AND NORMAL REACTION

                 REACTION STORED IN VECTOR -fi
 __________________________________________________________________________*/

 for (level=0; level<=curlev ; level++)
   for (i=1; i<=nodes[level] ; i++)
      {
        bn_qs[level][i] = 0.0;
        bn_qn[level][i] = 0.0;
      }


 for (level=0; level<=curlev ; level++)
	 {

    i1 = (level==curlev)?1:BFRsides[level]+1;
    for (i=i1 ; i<=Bsides[level] ; i++)
       {
         side = bside[level][i];

         n1 = side_nodes[level][2*side-1];
         n2 = side_nodes[level][2*side];
         lv1 = side_vlev[level][2*side-1];
         lv2 = side_vlev[level][2*side];

         cos1 = bs_cs[level][i];
         sin1 = bs_sn[level][i];

         Fx1 = -fi[lv1][2*n1-1];
         Fy1 = -fi[lv1][2*n1];
         Fx2 = -fi[lv2][2*n2-1];
         Fy2 = -fi[lv2][2*n2];

         Rn1 = -Fx1*sin1+Fy1*cos1;
         Rs1 =  Fx1*cos1+Fy1*sin1;

         Rn2 = -Fx2*sin1+Fy2*cos1;
         Rs2 =  Fx2*cos1+Fy2*sin1;

         pside = parent_side[level][2*i-1];
         i3 = xndof*(pside-1);

         if (side_bc[++i3])
           {
             bn_qs[lv1][n1] += Rs1;
             bn_qs[lv2][n2] += Rs2;
           }

         if (side_bc[++i3])
           {
             bn_qn[lv1][n1] += Rn1;
             bn_qn[lv2][n2] += Rn2;
           }
       }
   }

/*__________________________________________________________________________

        	            LOOP ON ALL BOUNDARY SIDES
       COMPUTATION OF SIDE AVERAGED TANGENTIAL AND NORMAL FORCES
___________________________________________________________________________*/

/* for (level=0; level<=curlev ; level++)
	 {

    i1 = (level==curlev)?1:BFRsides[level]+1;
    for (i=i1 ; i<=Bsides[level] ; i++)
       {
         side = bside[level][i];

         n1 = side_nodes[level][2*side-1];
         n2 = side_nodes[level][2*side];
         lv1 = side_vlev[level][2*side-1];
         lv2 = side_vlev[level][2*side];

       	i2 = nudim*(n1 - 1);
		   x1 = xyz[lv1][++i2];
		   y1 = xyz[lv1][++i2];

      	i2 = nudim*(n2 - 1);
      	x2 = xyz[lv2][++i2];
      	y2 = xyz[lv2][++i2];

         L12 = sqrt(x12*x12 + y12*y12);

         Rn1 = bn_qn[lv1][n1];
         Rs1 = bn_qs[lv1][n1];

         Rn2 = bn_qn[lv2][n2];
         Rs2 = bn_qs[lv2][n2];

         pside = parent_side[level][2*i-1];
         i3 = xndof*(pside-1);

         if (side_bc[++i3])
           {
             Fs[level][side] =
               0.5*(Rs1/bn_slength[lv1][n1]+Rs2/bn_slength[lv2][n2]);
           }

         if (side_bc[++i3])
           {
             Fn[level][side] =
               0.5*(Rn1/bn_nlength[lv1][n1]+Rn2/bn_nlength[lv2][n2]);
           }
       }
   }
*/


/*________________________________________________________________________

				             END OF FUNCTION
 ________________________________________________________________________*/

}/* end of tr1_flux_avg() */






/*******************************************[ USER LIBRARY : TR2DMLMF.C ]*****
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

 void tr1_zz_centsp (int curlev      , int *numel     , int *nodes,
                        int *FRnumel    , int **velm     , int *Bsides,
                        int *BFRsides,
                        int *PRsides    , int **eside    , int **bside,
                        int **tside     , int xndof      , double **xyz,
                        int **elm_nodes , int **elm_vlev , int **parent_elm,
                        int **side_nodes, int **side_vlev, double **uold,
                        int *mat        , double *prop   , int xnumpr,
                        double  **Sx     , double  **Sy    , double **Sxy,
                        double **a11    , double **a12   , double **a13,
                        double **a22    , double **a23   , double **a33,
                        double **b1x    , double **b2x   , double **b3x,
                        double **b1y    , double **b2y   , double **b3y,
                        double **b1xy   , double **b2xy  , double **b3xy,
                        int *nse        , int **idnode)


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

/*____________________________________________________________________________

								 BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

							LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int ie        ;   /* element counter for loop on elements.                 */
 int ie1       ;   /* initial element for loop on elements.                 */
 int level     ;   /* level   counter for loop on levels.                   */
 int nudim     ;   /* number of spatial dimensions (two for plane problems) */
 int n1, n2, n3;   /* element node numbers.                                 */
 int lv1, lv2, lv3;/* element node levels.                                  */
 int pelm      ;   /* parent element                                        */
 int side      ;   /* side number                                           */
 int i1        ;   /* auxiliar pointers for mesh arrays.                    */
 int i         ;   /* loop counter                                          */

/*____________________________________________________________________________

				 DECLARATION OF ELEMENT COORDINATES COEFFICIENTS
____________________________________________________________________________*/

 double x1, x2, x3;
 double y1, y2, y3;
 double x13, x32, y31, y23;
 double area;    /* two times the element area = x13*y23-x32*y31           */
 double xc, yc;

/*____________________________________________________________________________

				 DECLARATION OF ELEMENT DISPLACEMENTS COEFFICIENTS
____________________________________________________________________________*/

 double u1, u2, u3;
 double v1, v2, v3;
 double u13, u23, v13, v23;

/*____________________________________________________________________________

				 DECLARATION OF ELEMENT DISPLACEMENTS DERIVATIVES
____________________________________________________________________________*/

 double ux, uy, vx, vy;

/*____________________________________________________________________________

		               DECLARATION OF MATERIAL PROPERTIES
____________________________________________________________________________*/

 double  E    ;   /* elasticity modulus                                     */
 double  nu   ;   /* Poisson's ratio                                        */

 double  f1   ;   /* E/(1-nu*nu)                                            */
 double  f2   ;   /* E/2/(1+nu)                                             */

/*____________________________________________________________________________

  				     DECLARATION OF CENTROIDAL STRAINS AND STRESSES
____________________________________________________________________________*/

 double Exc, Eyc, Exyc;
 double Sxc, Syc, Sxyc;

/*____________________________________________________________________________

  				     DECLARATION OF SIDE STRESSES
____________________________________________________________________________*/

 double Sx1, Sy1, Sxy1;
 double Sx2, Sy2, Sxy2;

/*____________________________________________________________________________

			     DECLARATION OF MATRIX COEFFICIENTS FOR ZZ SMOOTHER
____________________________________________________________________________*/

 double A11, A12, A13, A22, A23, A33, A01, A02, A03, det;

 nudim = 2;

/*____________________________________________________________________________

            	  IDENTIFICATION OF BOUNDARY NODES
____________________________________________________________________________*/


	for (level=0; level<=curlev ; level++)
		 memset (idnode[level], 0, (nodes[level]+1)*sizeof(int));

   for (i=1 ; i<=Bsides[0] ; i++)
      {
        idnode[0][side_nodes[0][2*i-1]] = 1;
        idnode[0][side_nodes[0][2*i]] = 1;
      }

   for (level=1; level<=curlev ; level++)
	   {
        for (i=1 ; i<=BFRsides[level-1] ; i++)
           {
             side = eside[level-1][bside[level-1][i]];
             idnode[level][side] = 1;
           }
      }

/*__________________________________________________________________________

		             SET TO ZERO ALL NODAL COMPONENTS
___________________________________________________________________________*/

    for (level=0; level<=curlev ; level++)
       {
         memset (Sx[level] , 0, (nodes[level]+1)*sizeof(double ));
         memset (Sy[level] , 0, (nodes[level]+1)*sizeof(double ));
         memset (Sxy[level], 0, (nodes[level]+1)*sizeof(double ));

         memset (a11[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (a12[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (a13[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (a22[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (a23[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (a33[level] , 0, (nodes[level]+1)*sizeof(double));

         memset (b1x[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b2x[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b3x[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b1y[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b2y[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b3y[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b1xy[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b2xy[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b3xy[level] , 0, (nodes[level]+1)*sizeof(double));
       }


/*__________________________________________________________________________

			             LOOP ON ALL LEVELS
___________________________________________________________________________*/

 for (level=0; level<=curlev ; level++)
	 {

     ie1 = (level<curlev) ? FRnumel[level]+1 : 1;

/*__________________________________________________________________________

				  LOOP ON ALL UNREFINED ELEMENTS OF EACH PREVIOUS LEVEL
                    AND ON ALL ELEMENTS OF CURRENT LEVEL
___________________________________________________________________________*/

 for (ie=ie1; ie<=numel[level] ; ie++)
	 {

/*________________________________________________________________________

							  NODAL INCIDENCES
 __________________________________________________________________________*/

	i1 = 3*(velm[level][ie]-1);
	n1 = elm_nodes[level][++i1];
	n2 = elm_nodes[level][++i1];
	n3 = elm_nodes[level][++i1];

/*___________________________________________________________________________

								 VERTEX NODAL LEVELS
 __________________________________________________________________________*/


	i1 = 3*(velm[level][ie]-1);
	lv1= elm_vlev[level][++i1];
	lv2= elm_vlev[level][++i1];
	lv3= elm_vlev[level][++i1];

/*____________________________________________________________________________

              						  VERTEX	NODAL DISPLACEMENTS
____________________________________________________________________________*/


	i1 = xndof*(n1 - 1);
	u1 = uold[lv1][++i1];
	v1 = uold[lv1][++i1];

	i1 = xndof*(n2 - 1);
	u2 = uold[lv2][++i1];
	v2 = uold[lv2][++i1];

	i1 = xndof*(n3 - 1);
	u3 = uold[lv3][++i1];
	v3 = uold[lv3][++i1];

/*____________________________________________________________________________

						   DISPLACEMENTS COEFFICIENTS
____________________________________________________________________________*/

	u13 = u1-u3;
	u23 = u2-u3;

	v13 = v1-v3;
	v23 = v2-v3;

/*____________________________________________________________________________

							VERTEX NODAL COORDINATES
____________________________________________________________________________*/


	  i1 = nudim*(n1 - 1);
	  x1 = xyz[lv1][++i1];
   y1 = xyz[lv1][++i1];

   i1 = nudim*(n2 - 1);
   x2 = xyz[lv2][++i1];
   y2 = xyz[lv2][++i1];

 	 i1 = nudim*(n3 - 1);
 	 x3 = xyz[lv3][++i1];
 	 y3 = xyz[lv3][++i1];

/*____________________________________________________________________________

						   COORDINATES COEFFICIENTS
____________________________________________________________________________*/


	  x13 = x1-x3;
	  x32 = x3-x2;

   y31 = y3-y1;
	  y23 = y2-y3;

   area = x13*y23 - x32*y31;

/*____________________________________________________________________________

						     CENTROID COORDINATES
____________________________________________________________________________*/

    xc = (x1 + x2 + x3)/3.;
    yc = (y1 + y2 + y3)/3.;

/*__________________________________________________________________________

         					MATERIAL PROPERTIES AND STIFFNESS COEFICIENTS
___________________________________________________________________________*/

  pelm = parent_elm[level][3*velm[level][ie]-2];

  i1 = xnumpr * (mat[pelm] - 1);
  E  = prop[++i1];
  nu = prop[++i1];

  f1 = E/(1-nu*nu);
  f2 = 0.5*E/(1+nu);

/*____________________________________________________________________________

               ELEMENT DISPLACEMENTS DERIVATIVES
____________________________________________________________________________*/


  ux = (y23*u13+y31*u23)/area;
  uy = (x32*u13+x13*u23)/area;
  vx = (y23*v13+y31*v23)/area;
  vy = (x32*v13+x13*v23)/area;

/*____________________________________________________________________________

						       CENTROIDAL STRAINS
____________________________________________________________________________*/

  Exc  = ux;
  Eyc  = vy;
  Exyc = uy+vx;

/*____________________________________________________________________________

						            CENTROIDAL STRESSES
____________________________________________________________________________*/

  Sxc  = f1*(Exc+nu*Eyc);
  Syc  = f1*(nu*Exc+Eyc);
  Sxyc = f2*Exyc;

/*____________________________________________________________________________

						        STRESSES ON NODE N1
____________________________________________________________________________*/


  Sx[lv1][n1]  += Sxc;
  Sy[lv1][n1]  += Syc;
  Sxy[lv1][n1] += Sxyc;

  a11[lv1][n1] += 1.;
  a12[lv1][n1] += xc-x1;
  a13[lv1][n1] += yc-y1;
  a22[lv1][n1] += (xc-x1)*(xc-x1);
  a23[lv1][n1] += (xc-x1)*(yc-y1);
  a33[lv1][n1] += (yc-y1)*(yc-y1);

  b1x[lv1][n1] += Sxc;
  b2x[lv1][n1] += Sxc*(xc-x1);
  b3x[lv1][n1] += Sxc*(yc-y1);

  b1y[lv1][n1] += Syc;
  b2y[lv1][n1] += Syc*(xc-x1);
  b3y[lv1][n1] += Syc*(yc-y1);

  b1xy[lv1][n1] += Sxyc;
  b2xy[lv1][n1] += Sxyc*(xc-x1);
  b3xy[lv1][n1] += Sxyc*(yc-y1);

/*____________________________________________________________________________

						        STRESSES ON NODE N2
____________________________________________________________________________*/

  Sx[lv2][n2]  += Sxc;
  Sy[lv2][n2]  += Syc;
  Sxy[lv2][n2] += Sxyc;

  a11[lv2][n2] += 1.;
  a12[lv2][n2] += xc-x2;
  a13[lv2][n2] += yc-y2;
  a22[lv2][n2] += (xc-x2)*(xc-x2);
  a23[lv2][n2] += (xc-x2)*(yc-y2);
  a33[lv2][n2] += (yc-y2)*(yc-y2);

  b1x[lv2][n2] += Sxc;
  b2x[lv2][n2] += Sxc*(xc-x2);
  b3x[lv2][n2] += Sxc*(yc-y2);

  b1y[lv2][n2] += Syc;
  b2y[lv2][n2] += Syc*(xc-x2);
  b3y[lv2][n2] += Syc*(yc-y2);

  b1xy[lv2][n2] += Sxyc;
  b2xy[lv2][n2] += Sxyc*(xc-x2);
  b3xy[lv2][n2] += Sxyc*(yc-y2);

/*____________________________________________________________________________

						       STRESSES ON NODE N3
____________________________________________________________________________*/

  Sx[lv3][n3]  += Sxc;
  Sy[lv3][n3]  += Syc;
  Sxy[lv3][n3] += Sxyc;

  a11[lv3][n3] += 1.;
  a12[lv3][n3] += xc-x3;
  a13[lv3][n3] += yc-y3;
  a22[lv3][n3] += (xc-x3)*(xc-x3);
  a23[lv3][n3] += (xc-x3)*(yc-y3);
  a33[lv3][n3] += (yc-y3)*(yc-y3);

  b1x[lv3][n3] += Sxc;
  b2x[lv3][n3] += Sxc*(xc-x3);
  b3x[lv3][n3] += Sxc*(yc-y3);

  b1y[lv3][n3] += Syc;
  b2y[lv3][n3] += Syc*(xc-x3);
  b3y[lv3][n3] += Syc*(yc-y3);

  b1xy[lv3][n3] += Sxyc;
  b2xy[lv3][n3] += Sxyc*(xc-x3);
  b3xy[lv3][n3] += Sxyc*(yc-y3);

/*________________________________________________________________________

	        END OF LOOP ON UNREFINED ELEMENTS OF PREVIOUS LEVELS
                 AND ON ALL ELEMENTS OF CURRENT LEVEL
 ________________________________________________________________________*/

	 }

/*________________________________________________________________________

                     END OF LOOP ON LEVELS
 ________________________________________________________________________*/

 }

/*_________________________________________________________________________

          LOOP ON ALL NODES OF COARSE MESH (only boundary nodes needed)
 __________________________________________________________________________*/

      for (i=1 ; i<=nodes[0] ; i++)
         {
           Sx[0][i] /= nse[i];
           Sy[0][i] /= nse[i];
           Sxy[0][i] /= nse[i];
         }


/*_________________________________________________________________________

               LOOP ON ALL BOUNDARY NODES OF EACH LEVEL
 __________________________________________________________________________*/

 for (level=1; level<=curlev ; level++)
	 {
      for (i=1 ; i<=BFRsides[level-1] ; i++)
         {
           side = eside[level-1][bside[level-1][i]];
           Sx[level][side] /= 3.;
           Sy[level][side] /= 3.;
           Sxy[level][side] /= 3.;
         }
    }

/*_________________________________________________________________________

               LOOP ON ALL INTERIOR NODES OF EACH LEVEL
 __________________________________________________________________________*/

 for (level=0; level<=curlev ; level++)
	 {
      for (i=1 ; i<=nodes[level] ; i++)
         {
           if (!idnode[level][i])
             {
               A11 = a11[level][i];
               A12 = a12[level][i];
               A13 = a13[level][i];
               A22 = a22[level][i];
               A23 = a23[level][i];
               A33 = a33[level][i];
               det = A11*A22*A33+2*A12*A13*A23-A22*A13*A13-A12*A12*A33-
                     A11*A23*A23;

               A01 = A22*A33-A23*A23;
               A02 = A13*A23-A12*A33;
               A03 = A12*A23-A13*A22;

               Sx[level][i] = (A01*b1x[level][i]+A02*b2x[level][i]+
                               A03*b3x[level][i])/det;
               Sy[level][i] = (A01*b1y[level][i]+A02*b2y[level][i]+
                               A03*b3y[level][i])/det;
               Sxy[level][i] = (A01*b1xy[level][i]+A02*b2xy[level][i]+
                                A03*b3xy[level][i])/det;
             }
         }
    }

/*__________________________________________________________________________

	            LOOP ON ALL TRANSITION SIDES OF ALL LEVELS

       CORRECTION OF AVERAGED STRESSES (ONLY ONE ELEMENT CONTRIBUTION)

                 ON SON SIDES OF TRANSITION SIDES
___________________________________________________________________________*/


for (level=0 ; level<curlev ; level++)
   {
     if (PRsides[level])
       {
        for (i=1; i<=PRsides[level]; i++)
           {

           	 i1 = 2*(tside[level][i]-1);
           	 n1 = side_nodes[level][++i1];
	            n2 = side_nodes[level][++i1];

           	 i1 = 2*(tside[level][i]-1);
		           lv1= side_vlev[level][++i1];
		           lv2= side_vlev[level][++i1];

             Sx1 = Sx[lv1][n1];
             Sy1 = Sy[lv1][n1];
             Sxy1 = Sxy[lv1][n1];

             Sx2 = Sx[lv2][n2];
             Sy2 = Sy[lv2][n2];
             Sxy2 = Sxy[lv2][n2];

             side = eside[level][tside[level][i]];
             Sx [level+1][side] = (Sx1+Sx2)/2.;
             Sy [level+1][side] = (Sy1+Sy2)/2.;
             Sxy[level+1][side] = (Sxy1+Sxy2)/2.;

	          }
       }
   }

/*________________________________________________________________________

				             END OF FUNCTION
 ________________________________________________________________________*/

}/* end of tr1_zz_centsp() */





/*******************************************[ USER LIBRARY : TR2DMLMF.C ]*****
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

 void tr1_zz_sidesp (int curlev      , int *numel     , int *nodes,
                        int *FRnumel    , int **velm     , int *Bsides,
                        int *BFRsides,
                        int *PRsides    , int **eside    , int **bside,
                        int **tside     , int xndof      , double **xyz,
                        int **elm_nodes , int **elm_vlev , int **parent_elm,
                        int **side_nodes, int **side_vlev, double **uold,
                        int *mat        , double *prop   , int xnumpr,
                        double  **Sx     , double  **Sy    , double **Sxy,
                        double **a11    , double **a12   , double **a13,
                        double **a22    , double **a23   , double **a33,
                        double **b1x    , double **b2x   , double **b3x,
                        double **b1y    , double **b2y   , double **b3y,
                        double **b1xy   , double **b2xy  , double **b3xy,
                        int *nse        , int **idnode)


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

/*____________________________________________________________________________

								 BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

							LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int ie        ;   /* element counter for loop on elements.                 */
 int ie1       ;   /* initial element for loop on elements.                 */
 int level     ;   /* level   counter for loop on levels.                   */
 int nudim     ;   /* number of spatial dimensions (two for plane problems) */
 int n1, n2, n3;   /* element node numbers.                                 */
 int lv1, lv2, lv3;/* element node levels.                                  */
 int pelm      ;   /* parent element                                        */
 int side      ;   /* side number                                           */
 int i1        ;   /* auxiliar pointers for mesh arrays.                    */
 int i         ;   /* loop counter                                          */

/*____________________________________________________________________________

				 DECLARATION OF ELEMENT COORDINATES COEFFICIENTS
____________________________________________________________________________*/

 double x1, x2, x3;
 double y1, y2, y3;
 double x13, x32, y31, y23;
 double area;    /* two times the element area = x13*y23-x32*y31           */
 double xL1, xL2, xL3;
 double yL1, yL2, yL3;

/*____________________________________________________________________________

				 DECLARATION OF ELEMENT DISPLACEMENTS COEFFICIENTS
____________________________________________________________________________*/

 double u1, u2, u3;
 double v1, v2, v3;
 double u13, u23, v13, v23;

/*____________________________________________________________________________

				 DECLARATION OF ELEMENT DISPLACEMENTS DERIVATIVES
____________________________________________________________________________*/

 double ux, uy, vx, vy;

/*____________________________________________________________________________

		               DECLARATION OF MATERIAL PROPERTIES
____________________________________________________________________________*/

 double  E    ;   /* elasticity modulus                                     */
 double  nu   ;   /* Poisson's ratio                                        */

 double  f1   ;   /* E/(1-nu*nu)                                            */
 double  f2   ;   /* E/2/(1+nu)                                             */

/*____________________________________________________________________________

  				     DECLARATION OF CENTROIDAL STRAINS AND STRESSES
____________________________________________________________________________*/

 double Exc, Eyc, Exyc;
 double Sxc, Syc, Sxyc;

/*____________________________________________________________________________

  				     DECLARATION OF SIDE STRESSES
____________________________________________________________________________*/

 double Sx1, Sy1, Sxy1;
 double Sx2, Sy2, Sxy2;

/*____________________________________________________________________________

			     DECLARATION OF MATRIX COEFFICIENTS FOR ZZ SMOOTHER
____________________________________________________________________________*/

 double A11, A12, A13, A22, A23, A33, A01, A02, A03, det;

 nudim = 2;

/*____________________________________________________________________________

            	  IDENTIFICATION OF BOUNDARY NODES
____________________________________________________________________________*/


	for (level=0; level<=curlev ; level++)
		 memset (idnode[level], 0, (nodes[level]+1)*sizeof(int));

   for (i=1 ; i<=Bsides[0] ; i++)
      {
        idnode[0][side_nodes[0][2*i-1]] = 1;
        idnode[0][side_nodes[0][2*i]] = 1;
      }

   for (level=1; level<=curlev ; level++)
	   {
        for (i=1 ; i<=BFRsides[level-1] ; i++)
           {
             side = eside[level-1][bside[level-1][i]];
             idnode[level][side] = 1;
           }
      }

/*__________________________________________________________________________

		             SET TO ZERO ALL NODAL COMPONENTS
___________________________________________________________________________*/

    for (level=0; level<=curlev ; level++)
       {
         memset (Sx[level] , 0, (nodes[level]+1)*sizeof(double ));
         memset (Sy[level] , 0, (nodes[level]+1)*sizeof(double ));
         memset (Sxy[level], 0, (nodes[level]+1)*sizeof(double ));

         memset (a11[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (a12[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (a13[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (a22[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (a23[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (a33[level] , 0, (nodes[level]+1)*sizeof(double));

         memset (b1x[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b2x[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b3x[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b1y[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b2y[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b3y[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b1xy[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b2xy[level] , 0, (nodes[level]+1)*sizeof(double));
         memset (b3xy[level] , 0, (nodes[level]+1)*sizeof(double));
       }


/*__________________________________________________________________________

			             LOOP ON ALL LEVELS
___________________________________________________________________________*/

 for (level=0; level<=curlev ; level++)
	 {

     ie1 = (level<curlev) ? FRnumel[level]+1 : 1;

/*__________________________________________________________________________

				  LOOP ON ALL UNREFINED ELEMENTS OF EACH PREVIOUS LEVEL
                    AND ON ALL ELEMENTS OF CURRENT LEVEL
___________________________________________________________________________*/

 for (ie=ie1; ie<=numel[level] ; ie++)
	 {

/*________________________________________________________________________

							  NODAL INCIDENCES
 __________________________________________________________________________*/

	i1 = 3*(velm[level][ie]-1);
	n1 = elm_nodes[level][++i1];
	n2 = elm_nodes[level][++i1];
	n3 = elm_nodes[level][++i1];

/*___________________________________________________________________________

								 VERTEX NODAL LEVELS
 __________________________________________________________________________*/


	i1 = 3*(velm[level][ie]-1);
	lv1= elm_vlev[level][++i1];
	lv2= elm_vlev[level][++i1];
	lv3= elm_vlev[level][++i1];

/*____________________________________________________________________________

              						  VERTEX	NODAL DISPLACEMENTS
____________________________________________________________________________*/


	i1 = xndof*(n1 - 1);
	u1 = uold[lv1][++i1];
	v1 = uold[lv1][++i1];

	i1 = xndof*(n2 - 1);
	u2 = uold[lv2][++i1];
	v2 = uold[lv2][++i1];

	i1 = xndof*(n3 - 1);
	u3 = uold[lv3][++i1];
	v3 = uold[lv3][++i1];

/*____________________________________________________________________________

						   DISPLACEMENTS COEFFICIENTS
____________________________________________________________________________*/

	u13 = u1-u3;
	u23 = u2-u3;

	v13 = v1-v3;
	v23 = v2-v3;

/*____________________________________________________________________________

							VERTEX NODAL COORDINATES
____________________________________________________________________________*/


	  i1 = nudim*(n1 - 1);
	  x1 = xyz[lv1][++i1];
   y1 = xyz[lv1][++i1];

   i1 = nudim*(n2 - 1);
   x2 = xyz[lv2][++i1];
   y2 = xyz[lv2][++i1];

 	 i1 = nudim*(n3 - 1);
 	 x3 = xyz[lv3][++i1];
 	 y3 = xyz[lv3][++i1];

/*____________________________________________________________________________

						   COORDINATES COEFFICIENTS
____________________________________________________________________________*/


	  x13 = x1-x3;
	  x32 = x3-x2;

   y31 = y3-y1;
	  y23 = y2-y3;

   area = x13*y23 - x32*y31;

/*____________________________________________________________________________

              					  MID-SIDE COORDINATES
____________________________________________________________________________*/

    xL1 = (x1+x2)/2.;
    xL2 = (x2+x3)/2.;
    xL3 = (x3+x1)/2.;

    yL1 = (y1+y2)/2.;
    yL2 = (y2+y3)/2.;
    yL3 = (y3+y1)/2.;

/*__________________________________________________________________________

         					MATERIAL PROPERTIES AND STIFFNESS COEFICIENTS
___________________________________________________________________________*/

  pelm = parent_elm[level][3*velm[level][ie]-2];

  i1 = xnumpr * (mat[pelm] - 1);
  E  = prop[++i1];
  nu = prop[++i1];

  f1 = E/(1-nu*nu);
  f2 = 0.5*E/(1+nu);

/*____________________________________________________________________________

               ELEMENT DISPLACEMENTS DERIVATIVES
____________________________________________________________________________*/


  ux = (y23*u13+y31*u23)/area;
  uy = (x32*u13+x13*u23)/area;
  vx = (y23*v13+y31*v23)/area;
  vy = (x32*v13+x13*v23)/area;

/*____________________________________________________________________________

						       CENTROIDAL STRAINS
____________________________________________________________________________*/

  Exc  = ux;
  Eyc  = vy;
  Exyc = uy+vx;

/*____________________________________________________________________________

						            CENTROIDAL STRESSES
____________________________________________________________________________*/

  Sxc  = f1*(Exc+nu*Eyc);
  Syc  = f1*(nu*Exc+Eyc);
  Sxyc = f2*Exyc;

/*____________________________________________________________________________

						        STRESSES ON NODE N1 (CONTRIBUTIONS OF SIDE 1 AND SIDE 3)
____________________________________________________________________________*/


  Sx[lv1][n1]  += Sxc;
  Sy[lv1][n1]  += Syc;
  Sxy[lv1][n1] += Sxyc;

  a11[lv1][n1] += 0.5;
  a12[lv1][n1] += 0.5*(xL1-x1);
  a13[lv1][n1] += 0.5*(yL1-y1);
  a22[lv1][n1] += 0.5*(xL1-x1)*(xL1-x1);
  a23[lv1][n1] += 0.5*(xL1-x1)*(yL1-y1);
  a33[lv1][n1] += 0.5*(yL1-y1)*(yL1-y1);

  b1x[lv1][n1] += 0.5*Sxc;
  b2x[lv1][n1] += 0.5*Sxc*(xL1-x1);
  b3x[lv1][n1] += 0.5*Sxc*(yL1-y1);

  b1y[lv1][n1] += 0.5*Syc;
  b2y[lv1][n1] += 0.5*Syc*(xL1-x1);
  b3y[lv1][n1] += 0.5*Syc*(yL1-y1);

  b1xy[lv1][n1] += 0.5*Sxyc;
  b2xy[lv1][n1] += 0.5*Sxyc*(xL1-x1);
  b3xy[lv1][n1] += 0.5*Sxyc*(yL1-y1);

  a11[lv1][n1] += 0.5;
  a12[lv1][n1] += 0.5*(xL3-x1);
  a13[lv1][n1] += 0.5*(yL3-y1);
  a22[lv1][n1] += 0.5*(xL3-x1)*(xL3-x1);
  a23[lv1][n1] += 0.5*(xL3-x1)*(yL3-y1);
  a33[lv1][n1] += 0.5*(yL3-y1)*(yL3-y1);

  b1x[lv1][n1] += 0.5*Sxc;
  b2x[lv1][n1] += 0.5*Sxc*(xL3-x1);
  b3x[lv1][n1] += 0.5*Sxc*(yL3-y1);

  b1y[lv1][n1] += 0.5*Syc;
  b2y[lv1][n1] += 0.5*Syc*(xL3-x1);
  b3y[lv1][n1] += 0.5*Syc*(yL3-y1);

  b1xy[lv1][n1] += 0.5*Sxyc;
  b2xy[lv1][n1] += 0.5*Sxyc*(xL3-x1);
  b3xy[lv1][n1] += 0.5*Sxyc*(yL3-y1);

/*____________________________________________________________________________

						        STRESSES ON NODE N2 (CONTRIBUTIONS OF SIDE 1 AND SIDE 2)
____________________________________________________________________________*/

  Sx[lv2][n2]  += Sxc;
  Sy[lv2][n2]  += Syc;
  Sxy[lv2][n2] += Sxyc;

  a11[lv2][n2] += 0.5;
  a12[lv2][n2] += 0.5*(xL1-x2);
  a13[lv2][n2] += 0.5*(yL1-y2);
  a22[lv2][n2] += 0.5*(xL1-x2)*(xL1-x2);
  a23[lv2][n2] += 0.5*(xL1-x2)*(yL1-y2);
  a33[lv2][n2] += 0.5*(yL1-y2)*(yL1-y2);

  b1x[lv2][n2] += 0.5*Sxc;
  b2x[lv2][n2] += 0.5*Sxc*(xL1-x2);
  b3x[lv2][n2] += 0.5*Sxc*(yL1-y2);

  b1y[lv2][n2] += 0.5*Syc;
  b2y[lv2][n2] += 0.5*Syc*(xL1-x2);
  b3y[lv2][n2] += 0.5*Syc*(yL1-y2);

  b1xy[lv2][n2] += 0.5*Sxyc;
  b2xy[lv2][n2] += 0.5*Sxyc*(xL1-x2);
  b3xy[lv2][n2] += 0.5*Sxyc*(yL1-y2);

  a11[lv2][n2] += 0.5;
  a12[lv2][n2] += 0.5*(xL2-x2);
  a13[lv2][n2] += 0.5*(yL2-y2);
  a22[lv2][n2] += 0.5*(xL2-x2)*(xL2-x2);
  a23[lv2][n2] += 0.5*(xL2-x2)*(yL2-y2);
  a33[lv2][n2] += 0.5*(yL2-y2)*(yL2-y2);

  b1x[lv2][n2] += 0.5*Sxc;
  b2x[lv2][n2] += 0.5*Sxc*(xL2-x2);
  b3x[lv2][n2] += 0.5*Sxc*(yL2-y2);

  b1y[lv2][n2] += 0.5*Syc;
  b2y[lv2][n2] += 0.5*Syc*(xL2-x2);
  b3y[lv2][n2] += 0.5*Syc*(yL2-y2);

  b1xy[lv2][n2] += 0.5*Sxyc;
  b2xy[lv2][n2] += 0.5*Sxyc*(xL2-x2);
  b3xy[lv2][n2] += 0.5*Sxyc*(yL2-y2);


/*____________________________________________________________________________

						       STRESSES ON NODE N3  (CONTRIBUTIONS OF SIDE 2 AND SIDE 3)
____________________________________________________________________________*/

  Sx[lv3][n3]  += Sxc;
  Sy[lv3][n3]  += Syc;
  Sxy[lv3][n3] += Sxyc;

  a11[lv3][n3] += 0.5;
  a12[lv3][n3] += 0.5*(xL2-x3);
  a13[lv3][n3] += 0.5*(yL2-y3);
  a22[lv3][n3] += 0.5*(xL2-x3)*(xL2-x3);
  a23[lv3][n3] += 0.5*(xL2-x3)*(yL2-y3);
  a33[lv3][n3] += 0.5*(yL2-y3)*(yL2-y3);

  b1x[lv3][n3] += 0.5*Sxc;
  b2x[lv3][n3] += 0.5*Sxc*(xL2-x3);
  b3x[lv3][n3] += 0.5*Sxc*(yL2-y3);

  b1y[lv3][n3] += 0.5*Syc;
  b2y[lv3][n3] += 0.5*Syc*(xL2-x3);
  b3y[lv3][n3] += 0.5*Syc*(yL2-y3);

  b1xy[lv3][n3] += 0.5*Sxyc;
  b2xy[lv3][n3] += 0.5*Sxyc*(xL2-x3);
  b3xy[lv3][n3] += 0.5*Sxyc*(yL2-y3);

  a11[lv3][n3] += 0.5;
  a12[lv3][n3] += 0.5*(xL1-x3);
  a13[lv3][n3] += 0.5*(yL1-y3);
  a22[lv3][n3] += 0.5*(xL1-x3)*(xL1-x3);
  a23[lv3][n3] += 0.5*(xL1-x3)*(yL1-y3);
  a33[lv3][n3] += 0.5*(yL1-y3)*(yL1-y3);

  b1x[lv3][n3] += 0.5*Sxc;
  b2x[lv3][n3] += 0.5*Sxc*(xL1-x3);
  b3x[lv3][n3] += 0.5*Sxc*(yL1-y3);

  b1y[lv3][n3] += 0.5*Syc;
  b2y[lv3][n3] += 0.5*Syc*(xL1-x3);
  b3y[lv3][n3] += 0.5*Syc*(yL1-y3);

  b1xy[lv3][n3] += 0.5*Sxyc;
  b2xy[lv3][n3] += 0.5*Sxyc*(xL1-x3);
  b3xy[lv3][n3] += 0.5*Sxyc*(yL1-y3);


/*________________________________________________________________________

	        END OF LOOP ON UNREFINED ELEMENTS OF PREVIOUS LEVELS
                 AND ON ALL ELEMENTS OF CURRENT LEVEL
 ________________________________________________________________________*/

	 }

/*________________________________________________________________________

                     END OF LOOP ON LEVELS
 ________________________________________________________________________*/

 }

/*_________________________________________________________________________

          LOOP ON ALL NODES OF COARSE MESH (only boundary nodes needed)
 __________________________________________________________________________*/

      for (i=1 ; i<=nodes[0] ; i++)
         {
           Sx[0][i] /= nse[i];
           Sy[0][i] /= nse[i];
           Sxy[0][i] /= nse[i];
         }


/*_________________________________________________________________________

               LOOP ON ALL BOUNDARY NODES OF EACH LEVEL
 __________________________________________________________________________*/

 for (level=1; level<=curlev ; level++)
	 {
      for (i=1 ; i<=BFRsides[level-1] ; i++)
         {
           side = eside[level-1][bside[level-1][i]];
           Sx[level][side] /= 3.;
           Sy[level][side] /= 3.;
           Sxy[level][side] /= 3.;
         }
    }

/*_________________________________________________________________________

               LOOP ON ALL INTERIOR NODES OF EACH LEVEL
 __________________________________________________________________________*/

 for (level=0; level<=curlev ; level++)
	 {
      for (i=1 ; i<=nodes[level] ; i++)
         {
           if (!idnode[level][i])
             {
               A11 = a11[level][i];
               A12 = a12[level][i];
               A13 = a13[level][i];
               A22 = a22[level][i];
               A23 = a23[level][i];
               A33 = a33[level][i];
               det = A11*A22*A33+2*A12*A13*A23-A22*A13*A13-A12*A12*A33-
                     A11*A23*A23;

               A01 = A22*A33-A23*A23;
               A02 = A13*A23-A12*A33;
               A03 = A12*A23-A13*A22;

               Sx[level][i] = (A01*b1x[level][i]+A02*b2x[level][i]+
                               A03*b3x[level][i])/det;
               Sy[level][i] = (A01*b1y[level][i]+A02*b2y[level][i]+
                               A03*b3y[level][i])/det;
               Sxy[level][i] = (A01*b1xy[level][i]+A02*b2xy[level][i]+
                                A03*b3xy[level][i])/det;
             }
         }
    }

/*__________________________________________________________________________

	            LOOP ON ALL TRANSITION SIDES OF ALL LEVELS

       CORRECTION OF AVERAGED STRESSES (ONLY ONE ELEMENT CONTRIBUTION)

                 ON SON SIDES OF TRANSITION SIDES
___________________________________________________________________________*/


for (level=0 ; level<curlev ; level++)
   {
     if (PRsides[level])
       {
        for (i=1; i<=PRsides[level]; i++)
           {

           	 i1 = 2*(tside[level][i]-1);
           	 n1 = side_nodes[level][++i1];
	            n2 = side_nodes[level][++i1];

           	 i1 = 2*(tside[level][i]-1);
		           lv1= side_vlev[level][++i1];
		           lv2= side_vlev[level][++i1];

             Sx1 = Sx[lv1][n1];
             Sy1 = Sy[lv1][n1];
             Sxy1 = Sxy[lv1][n1];

             Sx2 = Sx[lv2][n2];
             Sy2 = Sy[lv2][n2];
             Sxy2 = Sxy[lv2][n2];

             side = eside[level][tside[level][i]];
             Sx [level+1][side] = (Sx1+Sx2)/2.;
             Sy [level+1][side] = (Sy1+Sy2)/2.;
             Sxy[level+1][side] = (Sxy1+Sxy2)/2.;

	          }
       }
   }

/*________________________________________________________________________

				             END OF FUNCTION
 ________________________________________________________________________*/

}/* end of tr1_zz_sidesp() */




/*******************************************[ USER LIBRARY : TR2DMLMF.C ]*****
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

 void tr1_side_avg (int curlev    , int *numel      , int *FRnumel,
                      int **velm      , int *sides       , int *Bsides,
                      int *BFRsides   , int *FRsides     , int *PRsides,
                      int **vside     , int **bside      , int **tside,
                      int **ssside    , int xndof        , double **xyz,
                      int **elm_sides , int **elm_nodes  , int **elm_vlev,
                      int **parent_elm,
                      double **uold   , int *mat         , double *prop,
                      int xnumpr      ,
                      double  **Sx     , double **Sy      , double **Sxy)


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

/*____________________________________________________________________________

								 BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

							LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int ie        ;   /* element counter for loop on elements.                 */
 int ie1       ;   /* initial element for loop on elements.                 */
 int level     ;   /* level   counter for loop on levels.                   */
 int nudim     ;   /* number of spatial dimensions (two for plane problems) */
 int L1, L2, L3;   /* element side numbers.                                 */
 int n1, n2, n3;   /* element node numbers.                                 */
 int lv1, lv2, lv3;/* element node levels.                                  */
 int pelm      ;   /* parent element                                        */
 int side      ;   /* side number                                           */
 int side1     ;   /* side number                                           */
 int side2     ;   /* side number                                           */
 int i1,i2     ;   /* auxiliar pointers for mesh arrays.                    */
 int i         ;   /* loop counter                                          */

/*____________________________________________________________________________

				 DECLARATION OF ELEMENT COORDINATES COEFFICIENTS
____________________________________________________________________________*/

 double x1, x2, x3;
 double y1, y2, y3;
 double x13, x32, y31, y23;
 double  area;   /* two times the element area = x13*y23-x32*y31           */

/*____________________________________________________________________________

				 DECLARATION OF ELEMENT DISPLACEMENTS COEFFICIENTS
____________________________________________________________________________*/

 double u1, u2, u3;
 double v1, v2, v3;
 double u13, u23, v13, v23;
/*____________________________________________________________________________

		               DECLARATION OF MATERIAL PROPERTIES
____________________________________________________________________________*/

 double  E    ;   /* elasticity modulus                                     */
 double  nu   ;   /* Poisson's ratio                                        */
 double  f1   ;   /* E/(1-nu*nu)                                            */
 double  f2   ;   /* E/2/(1+nu)                                             */

/*____________________________________________________________________________

  				     DECLARATION OF CENTROIDAL STRAINS AND STRESSES
____________________________________________________________________________*/

 double Exc, Eyc, Exyc;
 double Sxc, Syc, Sxyc;

/*____________________________________________________________________________

  				     DECLARATION OF SIDE TRACTIONS
____________________________________________________________________________*/

 nudim = 2;

/*__________________________________________________________________________

		             SET TO ZERO ALL FORCE COMPONENTS
___________________________________________________________________________*/


    for (level=0; level<=curlev ; level++)
       {
         memset (Sx[level] , 0, (sides[level]+1)*sizeof(double ));
         memset (Sy[level] , 0, (sides[level]+1)*sizeof(double ));
         memset (Sxy[level], 0, (sides[level]+1)*sizeof(double ));
       }



/*__________________________________________________________________________

			             LOOP ON ALL LEVELS
___________________________________________________________________________*/

 for (level=0; level<=curlev ; level++)
	 {

     ie1 = (level<curlev) ? FRnumel[level]+1 : 1;

/*__________________________________________________________________________

				  LOOP ON ALL UNREFINED ELEMENTS OF EACH PREVIOUS LEVEL
                    AND ON ALL ELEMENTS OF CURRENT LEVEL
___________________________________________________________________________*/

 for (ie=ie1; ie<=numel[level] ; ie++)
	 {

/*________________________________________________________________________

							  NODAL INCIDENCES
 __________________________________________________________________________*/

		i1 = 3*(velm[level][ie]-1);
		L1 = elm_sides[level][++i1];
		L2 = elm_sides[level][++i1];
		L3 = elm_sides[level][++i1];

		i1 = 3*(velm[level][ie]-1);
		n1 = elm_nodes[level][++i1];
		n2 = elm_nodes[level][++i1];
		n3 = elm_nodes[level][++i1];

/*___________________________________________________________________________

								 VERTEX NODAL LEVELS
 __________________________________________________________________________*/


		i1 = 3*(velm[level][ie]-1);
		lv1= elm_vlev[level][++i1];
		lv2= elm_vlev[level][++i1];
		lv3= elm_vlev[level][++i1];

/*____________________________________________________________________________

						  VERTEX	NODAL DISPLACEMENTS
____________________________________________________________________________*/


		i1 = xndof*(n1 - 1);
		u1 = uold[lv1][++i1];
		v1 = uold[lv1][++i1];

		i1 = xndof*(n2 - 1);
		u2 = uold[lv2][++i1];
		v2 = uold[lv2][++i1];

		i1 = xndof*(n3 - 1);
		u3 = uold[lv3][++i1];
		v3 = uold[lv3][++i1];

/*____________________________________________________________________________

						  DISPLACEMENTS COEFFICIENTS
____________________________________________________________________________*/

	   u13 = u1-u3;
	   u23 = u2-u3;

	   v13 = v1-v3;
	   v23 = v2-v3;

/*____________________________________________________________________________

								VERTEX NODAL COORDINATES
____________________________________________________________________________*/


		 i1 = nudim*(n1 - 1);
		 x1 = xyz[lv1][++i1];
		 y1 = xyz[lv1][++i1];

		 i1 = nudim*(n2 - 1);
		 x2 = xyz[lv2][++i1];
		 y2 = xyz[lv2][++i1];

		 i1 = nudim*(n3 - 1);
		 x3 = xyz[lv3][++i1];
   y3 = xyz[lv3][++i1];

/*____________________________________________________________________________

						  COORDINATES COEFFICIENTS
____________________________________________________________________________*/


	  x13 = x1-x3;
	  x32 = x3-x2;

   y31 = y3-y1;
	  y23 = y2-y3;

   area = x13*y23 - x32*y31;

/*__________________________________________________________________________

					MATERIAL PROPERTIES AND STIFFNESS COEFICIENTS
___________________________________________________________________________*/

		 pelm = parent_elm[level][3*velm[level][ie]-2];

		 i1 = xnumpr * (mat[pelm] - 1);
		 E  = prop[++i1];
		 nu = prop[++i1];


   f1 = E/(1-nu*nu);
   f2 = 0.5*E/(1+nu);

/*____________________________________________________________________________

						       CENTROIDAL STRAINS
____________________________________________________________________________*/

		 Exc  = (y23*u13+y31*u23)/area;
		 Eyc  = (x32*v13+x13*v23)/area;
   Exyc = (x32*u13+x13*u23+y23*v13+y31*v23)/area;

/*____________________________________________________________________________

						       CENTROIDAL STRESSES
____________________________________________________________________________*/

		 Sxc  = f1*(Exc+nu*Eyc);
		 Syc  = f1*(nu*Exc+Eyc);
   Sxyc = f2*Exyc;

/*____________________________________________________________________________

						                ASSIGN STRESSES TO SIDES
____________________________________________________________________________*/


		  i1 = abs(L1);
		  Sx[level][i1] += Sxc;
		  Sy[level][i1] += Syc;
	  Sxy[level][i1] += Sxyc;

		  i1 = abs(L2);
		  Sx[level][i1] += Sxc;
		  Sy[level][i1] += Syc;
	  Sxy[level][i1] += Sxyc;

		  i1 = abs(L3);
		  Sx[level][i1] += Sxc;
		  Sy[level][i1] += Syc;
	  Sxy[level][i1] += Sxyc;


/*________________________________________________________________________

	        END OF LOOP ON UNREFINED ELEMENTS OF PREVIOUS LEVELS
                 AND ON ALL ELEMENTS OF CURRENT LEVEL
 ________________________________________________________________________*/

	 }

/*________________________________________________________________________

                     END OF LOOP ON PREVIOUS LEVELS
 ________________________________________________________________________*/

 }


/*_________________________________________________________________________

          LOOP ON ALL UNREFINED INTERIOR SIDES OF PREVIOUS LEVELS
            (Since it is not possible to distinguish beetween
             interior and boundary sides, we loop first an all
             sides, and then the boundary sides are corrected)

          LOOP ON ALL UNREFINED BOUNDARY SIDES OF PREVIOUS LEVEL
                   BOUNDARY TRACTIONS ON FREE SIDES
 __________________________________________________________________________*/

 for (level=0; level<curlev ; level++)
	 {
    i1 = FRsides[level]+1;
    for (i=i1 ; i<=sides[level] ; i++)
       {
         side = vside[level][i];
         Sx[level][side] /= 2.;
         Sy[level][side] /= 2.;
         Sxy[level][side] /= 2.;
       }

    i1 = BFRsides[level]+1;
    for (i=i1 ; i<=Bsides[level] ; i++)
       {
         side = bside[level][i];
         Sx[level][side] *= 2.;
         Sy[level][side] *= 2.;
         Sxy[level][side] *= 2.;
       }

/*________________________________________________________________________

                     END OF LOOP ON PREVIOUS LEVELS
 ________________________________________________________________________*/

 }



/*__________________________________________________________________________

	            LOOP ON ALL INTERIOR SIDES OF CURRENT FINE LEVEL
                  (First loop on all sides and then
                    boundary sides are corrected)
___________________________________________________________________________*/


      for (i=1 ; i<=sides[curlev] ; i++)
         {
           Sx[curlev][i] /= 2.;
           Sy[curlev][i] /= 2.;
          Sxy[curlev][i] /= 2.;
         }


      for (i=1 ; i<=Bsides[curlev] ; i++)
         {
           side = bside[curlev][i];
           Sx[curlev][side] *= 2.;
           Sy[curlev][side] *= 2.;
          Sxy[curlev][side] *= 2.;
         }


/*__________________________________________________________________________

	            LOOP ON ALL TRANSITION SIDES OF ALL LEVELS

       CORRECTION OF AVERAGED FORCES (ONLY ONE ELEMENT CONTRIBUTION)

                 ON SON SIDES OF TRANSITION SIDES
___________________________________________________________________________*/

/*
 for (level=0; level<curlev ; level++)
	 {
       i2 = PRsides[level];
       for (i=1; i<=i2 ; i++)
         {
           side  = tside[level][i];
           side1 = ssside[level][side];
           side2 = side1+1;

           F1 = Fn[level][side];

           Fn[level+1][side1] = F1;
           Fn[level+1][side2] = F1;

           F1 = Fs[level][side];

           Fs[level+1][side1] = F1;
           Fs[level+1][side2] = F1;
         }
    }
*/

 for (level=0; level<curlev ; level++)
	 {
       i2 = PRsides[level];
       for (i=1; i<=i2 ; i++)
         {
           side  = tside[level][i];
           side1 = ssside[level][side];
           side2 = side1+1;

           Sx[level+1][side1] *= 2;
           Sx[level+1][side2] *= 2;

           Sy[level+1][side1] *= 2;
           Sy[level+1][side2] *= 2;

          Sxy[level+1][side1] *= 2;
          Sxy[level+1][side2] *= 2;

         }
    }


/*________________________________________________________________________

				             END OF FUNCTION
 ________________________________________________________________________*/

}/* end of tr1_side_avg() */

