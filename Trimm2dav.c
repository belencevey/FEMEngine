
/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: TR2DMLAV.C                                                 #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #        Routines for computation of local element refinament indicator.    #
 #        Multilevel computation for all elements.                           #
 #        Linear mixed triangular elements. Element code = 3.                      #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  tr2dmlav()   : computation of refinament indicator.                #
 #                                                                           #
 #****************************************************************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>


/*******************************************[ USER LIBRARY : TR2DMLAV.C ]*****
 #                                                                           #
 #  FUNCTION :  tr1_msavg_err()                                              #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Multilevel Computation of element refinament indicator.              #
 #      Error estimated in energy norm between recovered stresses and        #
 #      computed stresses from finite element results.                       #
 #      Recovered stresses computed from averaged midside tractions.         #
 #      Equidistribution of error in energy norm.                            #
 #      Zienckiewicz's strategy mean error.                                  #
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

void tr1_msavg_err (int curlev , int *numel     , int *FRnumel,
                      int **velm      , int xndof      , double **xyz,
                      int **elm_sides , int **elm_nodes, int **elm_vlev,
                      int **parent_elm, double **uold  , int *mat,
                      double *prop    , int xnumpr     ,
                      double  **Fn     , double  **Fs    ,
                      double **error  , double *Uerr   , double *Unew,
                      double *Ufea)

/*____________________________________________________________________________

							 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int   curlev      ; /* current level in analysis.                          */
// int  *numel       ; /* number of elements in each level.                   */
// int  *FRnumel     ; /* number of refined elements (all in group Rnumel).   */
// int  **velm       ; /* new element ordering after refinament [numel[level]]*/
// int    xndof      ; /* maximum number of nodal degrees of freedom.         */
// double **xyz      ; /* nodal coordinates     [nodes[curlev] * nudim]       */
// int  **elm_sides  ; /* sides incidences of each element [numel * xnes]     */
// int  **elm_nodes  ; /* nodes of each element    [numel[curlev]*xnen]       */
// int  **elm_vlev   ; /* level of element vertexs [numel * xnev]             */
// int  **parent_elm ; /* parent element and location [3*numel[curlev]]       */
// double **uold     ; /* vector of previous nodal displacements [xndof*nodes]*/
// int    *mat       ; /* element material type vector [numel].               */
// double *prop      ; /* material properties vector [numat * numpr].         */
// int     xnumpr    ; /* maximum number of material properties.              */
// double  **Fn       ; /* side normal tractions [sides[curlev]]               */
// double  **Fs       ; /* side tangential tractons [sides[curlev]]            */
// double **error    ; /* element error energy norm of each level             */
// double *Uerr   ;  /* energy of error                                       */
// double *Unew   ;  /* energy of new displacements for all elements          */
// double *Ufea   ;   /* energy of finite element approximation               */

/*____________________________________________________________________________

								 BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

							LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int ie        ;   /* element counter for loop on elements.                 */
 int ie1       ;   /* initial element for loop on elements.                 */
 int ie2       ;   /* last    element for loop on elements.                 */
 int level     ;   /* level   counter for loop on levels.                   */
 int nudim     ;   /* number of spatial dimensions (two for plane problems) */
 int L1, L2, L3;   /* element side numbers.                                 */
 int n1, n2, n3;   /* element node numbers.                                 */
 int lv1, lv2, lv3;/* element node levels.                                  */
 int pelm      ;   /* parent element                                        */
 int i1        ;   /* auxiliar pointers for mesh arrays.                    */
 double Ueerr  ;   /* energy of error in each element                       */
 double Ufem   ;   /* energy of current approx. in each element             */

/*____________________________________________________________________________

				 DECLARATION OF ELEMENT COORDINATES COEFFICIENTS
____________________________________________________________________________*/

 double x1, x2, x3;
 double y1, y2, y3;
 double x12, x13, x32, y12, y31, y23;
 double  area;   /* two times the element area = x13*y23-x32*y31           */

 double L12, L23, L31; /* side lengths                                     */
 double cos1, sin1; /* cosine and sine of side 1                           */
 double cos2, sin2; /* cosine and sine of side 2                           */
 double cos3, sin3; /* cosine and sine of side 3                           */

 double dus, dvs;

/*____________________________________________________________________________

				 DECLARATION OF ELEMENT DISPLACEMENTS COEFFICIENTS
____________________________________________________________________________*/

 double u1, u2, u3;
 double v1, v2, v3;
 double u12, u13, u23, v12, v13, v23;

/*____________________________________________________________________________

		               DECLARATION OF MATERIAL PROPERTIES
____________________________________________________________________________*/

 double  E    ;   /* elasticity modulus                                     */
 double  nu   ;   /* Poisson's ratio                                        */
 double  th   ;   /* element thickness                                      */
 double  f1   ;   /* E/(1-nu*nu)                                            */
 double  f2   ;   /* E/2/(1+nu)                                             */

/*____________________________________________________________________________

  				     DECLARATION OF SIDE TRACTIONS
____________________________________________________________________________*/

 double FnL1, FsL1, FnL2, FsL2, FnL3, FsL3;

 /*____________________________________________________________________________

  				     DECLARATION OF MIDSIDE STRESSES
____________________________________________________________________________*/

 double SxL1, SyL1, SxyL1;
 double SxL2, SyL2, SxyL2;
 double SxL3, SyL3, SxyL3;

/*____________________________________________________________________________

  				     DECLARATION OF CENTROIDAL STRAINS AND STRESSES
____________________________________________________________________________*/

 double Exc, Eyc, Exyc;
 double Sxc, Syc, Sxyc;

/*____________________________________________________________________________

  				     DECLARATION OF MIDSIDE STRESSES AND SIDE STRAIN
____________________________________________________________________________*/

 double Ss, Sn, Sns, Es;


/*____________________________________________________________________________

  				     DECLARATION OF ERROR STRESSES AT INTEGRATION POINTS
____________________________________________________________________________*/

 double Isx1, Isy1, Isxy1;
 double Isx2, Isy2, Isxy2;
 double Isx3, Isy3, Isxy3;

/*____________________________________________________________________________

					NUMBER OF DIMENSION FOR NODAL COORDINATES
____________________________________________________________________________*/

 nudim = 2;

/*__________________________________________________________________________

         COMPUTATION OF MESH DEFORMATION ENERGY OF CURRENT LEVEL

              AND ESTIMATED ELEMENT ERROR IN ENERGY NORM
___________________________________________________________________________*/

 *Unew = 0;
 *Uerr = 0;
 *Ufea = 0;

/*__________________________________________________________________________

			                 LOOP ON ALL LEVELS
___________________________________________________________________________*/

 for (level=0; level<=curlev ; level++)
	 {

      ie1 = (level<curlev) ? FRnumel[level]+1 : 1;
      ie2 = numel[level];

/*__________________________________________________________________________

				  LOOP ON ALL UNREFINED ELEMENTS OF PREVIOUS LEVELS
                   AND ON ALL ELEMENTS OF CURRENT LEVEL
___________________________________________________________________________*/

 for (ie=ie1; ie<=ie2 ; ie++)
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

    u12 = u1-u2;
	   u13 = u1-u3;
	   u23 = u2-u3;

    v12 = v1-v2;
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

	   x12 = x1-x2;
	   y12 = y1-y2;

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

						       STRESSES ON SIDE L1
____________________________________________________________________________*/


    	i1 = abs(L1);
	   FnL1 = Fn[level][i1];
	   FsL1 = Fs[level][i1];

      L12 = x12*x12 + y12*y12;

      cos1 = x12;
      sin1 = y12;
      dus  = u12;
      dvs  = v12;

      Es = (dus*cos1 + dvs*sin1)/L12;

      Sns = FsL1/th;
      Sn  = FnL1/th;
      Ss  = nu*Sn + Es*E;

      SxL3  = (Ss*cos1*cos1 + Sn*sin1*sin1 - 2*Sns*sin1*cos1)/L12;
      SyL3  = (Ss*sin1*sin1 + Sn*cos1*cos1 + 2*Sns*sin1*cos1)/L12;
      SxyL3 = ((Ss - Sn)*sin1*cos1 + Sns*(cos1*cos1-sin1*sin1))/L12;


/*____________________________________________________________________________

						       STRESSES ON SIDE L2
____________________________________________________________________________*/


      i1 = abs(L2);
		FnL2 = Fn[level][i1];
		FsL2 = Fs[level][i1];

      L23 = x32*x32 + y23*y23;

      cos2 = -x32;
      sin2 =  y23;
      dus  =  u23;
      dvs  =  v23;

      Es = (dus*cos2 + dvs*sin2)/L23;

      Sns = FsL2/th;
      Sn  = FnL2/th;
      Ss  = nu*Sn + Es*E;

      SxL1  = (Ss*cos2*cos2 + Sn*sin2*sin2 - 2*Sns*sin2*cos2)/L23;
      SyL1  = (Ss*sin2*sin2 + Sn*cos2*cos2 + 2*Sns*sin2*cos2)/L23;
      SxyL1 = ((Ss - Sn)*sin2*cos2 + Sns*(cos2*cos2-sin2*sin2))/L23;

/*____________________________________________________________________________

						       STRESSES ON SIDE L3
____________________________________________________________________________*/


    	i1 = abs(L3);
		FnL3 = Fn[level][i1];
		FsL3 = Fs[level][i1];

      L31 = x13*x13 + y31*y31;

      cos3 = -x13;
      sin3 =  y31;
      dus  = -u13;
      dvs  = -v13;

      Es = (dus*cos3 + dvs*sin3)/L31;

      Sns = FsL3/th;
      Sn  = FnL3/th;
      Ss  = nu*Sn + Es*E;

      SxL2  = (Ss*cos3*cos3 + Sn*sin3*sin3 - 2*Sns*sin3*cos3)/L31;
      SyL2  = (Ss*sin3*sin3 + Sn*cos3*cos3 + 2*Sns*sin3*cos3)/L31;
      SxyL2 = ((Ss - Sn)*sin3*cos3 + Sns*(cos3*cos3-sin3*sin3))/L31;


/*__________________________________________________________________________

							          SQUARE ERROR ENERGY NORM
        (NUMERICAL INTEGRATED AT MIDSIDES, WEIGTHS = 1/3)
___________________________________________________________________________*/

      Isx1 = SxL1-Sxc;
      Isx2 = SxL2-Sxc;
      Isx3 = SxL3-Sxc;

      Isy1 = SyL1-Syc;
      Isy2 = SyL2-Syc;
      Isy3 = SyL3-Syc;

      Isxy1 = SxyL1-Sxyc;
      Isxy2 = SxyL2-Sxyc;
      Isxy3 = SxyL3-Sxyc;

      Ueerr = (th*area/6./E)*
              ((Isx1*Isx1-2*nu*Isx1*Isy1+Isy1*Isy1+2*(1+nu)*Isxy1*Isxy1)+
               (Isx2*Isx2-2*nu*Isx2*Isy2+Isy2*Isy2+2*(1+nu)*Isxy2*Isxy2)+
               (Isx3*Isx3-2*nu*Isx3*Isy3+Isy3*Isy3+2*(1+nu)*Isxy3*Isxy3));

/*__________________________________________________________________________

							 ENERGY OF FE APPROXIMATION
___________________________________________________________________________*/


      Ufem = th*(area/2)*(Sxc*Exc+Syc*Eyc+Sxyc*Exyc);

      *Ufea += Ufem;

/*__________________________________________________________________________

						ACCUMULATED	ENERGY OF NEW DISPLACEMENTS
___________________________________________________________________________*/

	   *Unew += Ufem + Ueerr;

/*__________________________________________________________________________

						   ENERGY ERROR OF CURRENT ELEMENT
___________________________________________________________________________*/

	   error[level][velm[level][ie]] = sqrt(Ueerr);

    *Uerr += Ueerr;

/*________________________________________________________________________

	        END OF LOOP ON UNREFINED ELEMENTS OF PREVIOUS LEVELS
                 AND ON ALL ELEMENTS OF CURRENT LEVEL
 ________________________________________________________________________*/

	 }

/*________________________________________________________________________

                     END OF LOOP ON ALL LEVELS
 ________________________________________________________________________*/

 }

/*________________________________________________________________________

				              END OF FUNCTION
 ________________________________________________________________________*/


}/* end of tr1_msavg_err() */


