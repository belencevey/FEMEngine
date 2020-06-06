
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

void trim2ds (FILE *fpo, FILE *fGiD    , int *acnodes,
                      int elm1, int elm2, int xndof, int nudim,
                      int *elm_sides , int *elm_nodes, int *elm_vlev,
                      int *parent_elm, int *velm     , double **xyz,
                      int *mat       , double *prop  , int xnumpr,
                      double **uold  , double  *Fn    , double  *Fs,
                      double **Sx     , double **Sy    , double **Sxy,
                      double *error  , double errm   ,int ESTIMATOR,
                      int level)


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


/*____________________________________________________________________________

			                BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

				          LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int ie        ;   /* element counter for loop on elements.                 */
 int L1, L2, L3;   /* element side numbers.                                 */
 int n1, n2, n3;   /* element node numbers.                                 */
 int lv1, lv2, lv3;/* element node levels.                                  */
 int pelm      ;   /* parent element                                        */
 int i1        ;   /* auxiliar pointers for mesh arrays.                    */

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

 double  x13, y31, x32, y23, x12, y12;

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

/* double  f1   ;    E/(1-nu*nu)                                            */
/* double  f2   ;    E/2/(1+nu)                                             */
/*double  area ;    two times the element area = x13*y23-x32*y31           */

/*____________________________________________________________________________

  				     DECLARATION OF SIDE TRACTIONS
____________________________________________________________________________*/

 double FnL1, FsL1, FnL2, FsL2, FnL3, FsL3;

 /*____________________________________________________________________________

  				     DECLARATION OF MIDSIDE STRESSES (GLOBAL COORDINATES)
____________________________________________________________________________*/

 double SxL1, SyL1, SxyL1;
 double SxL2, SyL2, SxyL2;
 double SxL3, SyL3, SxyL3;

/*____________________________________________________________________________

  				     DECLARATION OF CENTROIDAL STRAINS AND STRESSES
____________________________________________________________________________*/

/*
 double Exc, Eyc, Exyc;
 double Sxc, Syc, Sxyc;
*/

/*____________________________________________________________________________

  				     DECLARATION OF MIDSIDE STRESSES AND SIDE STRAIN
____________________________________________________________________________*/

 double Ss, Sn, Sns, Es;


/*____________________________________________________________________________

  				     DECLARATION OF RECOVERED VERTEX STRESSES
____________________________________________________________________________*/

 double Sx1, Sy1, Sxy1;
 double Sx2, Sy2, Sxy2;
 double Sx3, Sy3, Sxy3;

/*____________________________________________________________________________

  				     MIDSIDE ESTIMATOR
____________________________________________________________________________*/

// SIDE FLUX AVERAGING (DISCONTINUOUS) or SIDE AVERAGING (DISCONTINUOUS)
if (ESTIMATOR==2 || ESTIMATOR==5)
  {

/*__________________________________________________________________________

					LOOP ON UNREFINED ELEMENTS OF THIS LEVEL
___________________________________________________________________________*/

 for (ie=elm1; ie<=elm2 ; ie++)
	 {
/*________________________________________________________________________

							  NODAL INCIDENCES
 __________________________________________________________________________*/


/*________________________________________________________________________

							 SIDE AND NODAL INCIDENCES
 __________________________________________________________________________*/

		i1 = 3*(velm[ie]-1);
		L1 = elm_sides[++i1];
		L2 = elm_sides[++i1];
		L3 = elm_sides[++i1];

		i1 = 3*(velm[ie]-1);
		n1 = elm_nodes[++i1];
		n2 = elm_nodes[++i1];
		n3 = elm_nodes[++i1];

/*___________________________________________________________________________

								 VERTEX NODAL LEVELS
 __________________________________________________________________________*/


		i1 = 3*(velm[ie]-1);
		lv1= elm_vlev[++i1];
		lv2= elm_vlev[++i1];
		lv3= elm_vlev[++i1];


/*____________________________________________________________________________

								NODAL COORDINATES
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

							 VERTEX NODAL DISPLACEMENTS
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

						  COORDINATES COEFFICIENTS
____________________________________________________________________________*/

		x13 = x1 - x3;
		y31 = y3 - y1;

		x32 = x3 - x2;
		y23 = y2 - y3;

		x12 = x1 - x2;
		y12 = y1 - y2;

/*
  area = x13*y23 - x32*y31;
*/

/*__________________________________________________________________________

					MATERIAL PROPERTIES AND STIFFNESS COEFICIENTS
___________________________________________________________________________*/

		pelm = parent_elm[3*velm[ie]-2];

		i1 = xnumpr * (mat[pelm] - 1);
		E  = prop[++i1];
		nu = prop[++i1];
        th = prop[++i1];

/*
  f1 = E/(1-nu*nu);
  f2 = 0.5*E/(1+nu);
*/

/*____________________________________________________________________________

						       CENTROIDAL STRAINS
____________________________________________________________________________*/

/*
		Exc  = (y23*u13+y31*u23)/area;
		Eyc  = (x32*v13+x13*v23)/area;
        Exyc = (x32*u13+x13*u23+y23*v13+y31*v23)/area;
*/

/*____________________________________________________________________________

						       CENTROIDAL STRESSES
____________________________________________________________________________*/

/*
		Sxc  = f1*(Exc+nu*Eyc);
		Syc  = f1*(nu*Exc+Eyc);
        Sxyc = f2*Exyc;
*/

/*____________________________________________________________________________

						       STRESSES ON SIDE L1
____________________________________________________________________________*/

if (ESTIMATOR==2) // SIDE FLUX AVERAGING
  {
      i1 = abs(L1);
	  FnL1 = Fn[i1];
	  FsL1 = Fs[i1];

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
  }

if (ESTIMATOR==5) // GLOBAL SIDE FLUX AVERAGING
  {
    		i1 = abs(L1);
		    SxL3  =  Sx[level][i1];
		    SyL3  =  Sy[level][i1];
		    SxyL3 = Sxy[level][i1];
  }

/*____________________________________________________________________________

						       STRESSES ON SIDE L2
____________________________________________________________________________*/


if (ESTIMATOR==2) // SIDE FLUX AVERAGING
  {
    		i1 = abs(L2);
		    FnL2 = Fn[i1];
		    FsL2 = Fs[i1];

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
  }

if (ESTIMATOR==5) // GLOBAL SIDE FLUX AVERAGING
  {
    		i1 = abs(L2);
		    SxL1  =  Sx[level][i1];
		    SyL1  =  Sy[level][i1];
		    SxyL1 = Sxy[level][i1];
  }

/*____________________________________________________________________________

						       STRESSES ON SIDE L3
____________________________________________________________________________*/


if (ESTIMATOR==2) // SIDE FLUX AVERAGING
  {
    		i1 = abs(L3);
		    FnL3 = Fn[i1];
		    FsL3 = Fs[i1];

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
  }

if (ESTIMATOR==5)  // GLOBAL SIDE FLUX AVERAGING
  {
    		i1 = abs(L3);
		    SxL2  =  Sx[level][i1];
		    SyL2  =  Sy[level][i1];
		    SxyL2 = Sxy[level][i1];
  }


/*____________________________________________________________________________

						   RECOVERED NODAL STRESSES AT VERTEXS
____________________________________________________________________________*/


		Sx1  = -SxL1  + SxL2  + SxL3;
		Sy1  = -SyL1  + SyL2  + SyL3;
        Sxy1 = -SxyL1 + SxyL2 + SxyL3;

		Sx2  = SxL1  - SxL2  + SxL3;
		Sy2  = SyL1  - SyL2  + SyL3;
        Sxy2 = SxyL1 - SxyL2 + SxyL3;

		Sx3  = SxL1  + SxL2  - SxL3;
		Sy3  = SyL1  + SyL2  - SyL3;
        Sxy3 = SxyL1 + SxyL2 - SxyL3;


/*____________________________________________________________________________

          						 PRINT ELEMENT STRESSES AND ERROR
____________________________________________________________________________*/


	 fprintf (fpo,"%d  %d  %d  ", acnodes[lv1]+n1, acnodes[lv2]+n2,
             acnodes[lv3]+n3);

	 fprintf (fpo,"%+14.8e   %+14.8e   %+14.8e  ", Sx1, Sx2, Sx3);
	 fprintf (fpo,"%+14.8e   %+14.8e   %+14.8e  ", Sy1, Sy2, Sy3);
	 fprintf (fpo,"%+14.8e   %+14.8e   %+14.8e  ", Sxy1, Sxy2, Sxy3);
  fprintf (fpo,"%+14.8e ", error[velm[ie]]/errm);
	 fprintf (fpo,"\n");

/*____________________________________________________________________________

          						 PRINT ELEMENT STRESSES GiD
____________________________________________________________________________*/

  fprintf (fGiD,"\n %d %f %f %f",ie, SxL1, SyL1, SxyL1);
  fprintf (fGiD,"\n %f %f %f",SxL2, SyL2, SxyL2);
  fprintf (fGiD,"\n %f %f %f",SxL3, SyL3, SxyL3);

/*________________________________________________________________________

		                  END OF LOOP ON ELEMENTS
 ________________________________________________________________________*/

	 }

/*________________________________________________________________________

		                 END OF MIDSIDE ESTIMATOR
 ________________________________________________________________________*/

	 }



/*____________________________________________________________________________

  				     NODAL AVERAGE ESTIMATOR AND ZZ SUPERCONVERGENCE
____________________________________________________________________________*/


if (ESTIMATOR==1 || ESTIMATOR==3 || ESTIMATOR==4)
  {

/*__________________________________________________________________________

					LOOP ON UNREFINED ELEMENTS OF THIS LEVEL
___________________________________________________________________________*/

 for (ie=elm1; ie<=elm2 ; ie++)
	 {
/*________________________________________________________________________

							  NODAL INCIDENCES
 __________________________________________________________________________*/

		i1 = 3*(velm[ie]-1);
		n1 = elm_nodes[++i1];
		n2 = elm_nodes[++i1];
		n3 = elm_nodes[++i1];

/*___________________________________________________________________________

								 VERTEX NODAL LEVELS
 __________________________________________________________________________*/


		i1 = 3*(velm[ie]-1);
		lv1= elm_vlev[++i1];
		lv2= elm_vlev[++i1];
		lv3= elm_vlev[++i1];

/*____________________________________________________________________________

						       STRESSES ON NODE N1
____________________________________________________________________________*/

  Sx1 = Sx[lv1][n1];
  Sy1 = Sy[lv1][n1];
  Sxy1 = Sxy[lv1][n1];
/*____________________________________________________________________________

						       STRESSES ON NODE N2
____________________________________________________________________________*/

  Sx2 = Sx[lv2][n2];
  Sy2 = Sy[lv2][n2];
  Sxy2 = Sxy[lv2][n2];
/*____________________________________________________________________________

						       STRESSES ON NODE N3
____________________________________________________________________________*/

  Sx3 = Sx[lv3][n3];
  Sy3 = Sy[lv3][n3];
  Sxy3 = Sxy[lv3][n3];

/*____________________________________________________________________________

						 PRINT ELEMENT INCIDENCES AND STRESSES
____________________________________________________________________________*/


	 fprintf (fpo,"%d  %d  %d  ", acnodes[lv1]+n1, acnodes[lv2]+n2,
             acnodes[lv3]+n3);

	 fprintf (fpo,"%+14.8e   %+14.8e   %+14.8e  ", Sx1, Sx2, Sx3);
	 fprintf (fpo,"%+14.8e   %+14.8e   %+14.8e  ", Sy1, Sy2, Sy3);
	 fprintf (fpo,"%+14.8e   %+14.8e   %+14.8e  ", Sxy1, Sxy2, Sxy3);
     fprintf (fpo,"%+14.8e ", error[velm[ie]]/errm);
	 fprintf (fpo,"\n");

	 /*____________________________________________________________________________

	 PRINT ELEMENT STRESSES GiD
	 ____________________________________________________________________________*/

	 fprintf(fGiD, "\n %d %f %f %f", ie, Sx1, Sy1, Sxy1);
	 fprintf(fGiD, "\n    %f %f %f", Sx2, Sy2, Sxy2);
	 fprintf(fGiD, "\n    %f %f %f", Sx3, Sy3, Sxy3);


/*________________________________________________________________________

		  END OF LOOP ON ELEMENTS
 ________________________________________________________________________*/

	 }

/*________________________________________________________________________

		  END OF NODAL AVERAGE
 ________________________________________________________________________*/

	 }




/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of trim2ds() */
