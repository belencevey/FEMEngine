
/****************************[ MAIN FUNCTION ]********************************
 #                                                                           #
 #  PROGRAM NAME: SAMFE2D                                                    #
 #                                                                           #
 #  C LIBRARY: LSAFE.C                                                       #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #          Structural Analysis  by  Mixed Finite  Elements  in 2D           #
 #                                                                           #
 #                                                                           #
 #                S A M F E 2 D P (DOS version 1.0 - 30/07/02)               #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #                                                                           #
 #          Finite Element Structural Analysis for adaptive                  #
 #          solution of linear problems in two dimensions.                   #
 #          Stress and plate elements.                                       #
 #          Only one type de element for the entire mesh.                    #
 #          Multigrid adaptive solver.                                       #
 #          Direct Gauss solver implemented for skyline storage              #
 #          for initial coarse mesh.                                         #
 #          Energy norm used for error in adaptive refinament and            #
 #          for controlling the iterative solution on the fine mesh.         #
 #          Uniform refinament with irregular nodes, imcompatible            #
 #          meshes with constrained irregular nodes.                         #
 #          Curved boundary considered quadratic isoparametric description   #
 #          via control points.                                              #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED LIBRARY FUNCTIONS USED ]--------------------------------#
 #                                                                           #
 #                                                                           #
 #  INIT_LFE.C :                                                             #
 #    void   initp()     : initiates program. Display presentation screen.   #
 #    void   prtdis()    : print displacements in output file.               #
 #                                                                           #
 #  FILE_LFE.C :                                                             #
 #    void   iofile()    : open input/output files.                          #
 #                                                                           #
 #  MESHDAT.C  :                                                             #
 #    void   meshdat()   : read mesh data.                                   #
 #                                                                           #
 #  GENVEC.C   :                                                             #
 #    void  nodvec()     : read nodal vector data.                           #
 #    void  sidevec()    : read side vector surface load data.               #
 #    void  elmvec()     : read element body load vector data.               #
 #    void  sideload()   : computation of boundary nodal loads initial mesh. #
 #    void  mlsideload() : multilevel computation of boundary nodal loads.   #
 #    void  elmld()      : computation of nodal loads over elements.         #
 #    void  genvec()     : generation of  nodal vector data.                 #
 #                                                                           #
 #  SKYLINE.C  :                                                             #
 #    void   numeqn()    : sequential numeration of free dofs.               #
 #    void   profile()   : computation of pointers of skyline.               #
 #    void   expand()    : expand compacted vectors into global ones.        #
 #    void   compress()  : compress global vectors into compacted ones.      #
 #                                                                           #
 #  GAUSS  .C  :                                                             #
 #    void   decomp()    : Gauss factorization skyline storage.              #
 #    void   redbak()    : Gauss backsubstitution for skyline storage.       #
 #                                                                           #
 #  STIFFBLK.C :                                                             #
 #     void  stiffblk()  : assembling of global stiffness matrix.            #
 #                           for blocks of elements.                         #
 #                                                                           #
 #  AUX_LFE .C :                                                             #
 #     void  *mallocvf() : dynamic memory allocation.                        #
 #     void  *callocvf() : dynamic zeroing memory allocation.                #
 #                                                                           #
 #  PRTPLT  .C :                                                             #
 #     void  prtplt()    : computation of element stresses.                  #
 #                                                                           #
 #  DATALEV .C :                                                             #
 #     void  refelm()    : computation of number of elements to be refined.  #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO .H :                                                               #
 #    int  fscanf() : performs formatted input from a stream.                #
 #    int  fclose() : closes a stream.                                       #
 #    int  printf() : formatted output to standard output device.            #
 #                                                                           #
 #  MEM   .H :                                                               #
 #    void *memcpy(): copies a block of bytes.                               #
 #                                                                           #
 #  STDLIB.H :                                                               #
 #    void exit()   : terminates program.                                    #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/



#include <new>
//#include <cstdio>
//#include <stdlib.h>
#include <cstdlib>
//#include <memory.h>
//#include <mem.h>
#include <cmath>
#include <string>
///////////////////////////////

#include<fstream>
#include<iostream>


//Librerias PROPIAS//
	//#include "INIT_LFE.h"
	#include "FILE_LFE.h"
	#include "MESHDAT.h"
	#include "GENVEC.h"
	#include "SKYLINE.h"
	#include "GAUSS.h"
	#include "STIFFBLK.h"
	#include "AUX_LFE.h"
	#include "PRTPLT.h"


	#include "Inputest.h"
	#include "Inputpar.h"
	#include "REFIDFR.h"
	#include "TR2DMLMF.h"
	#include "TRI3_K_SBS.h"
	#include "TRIM2DK.h"
	#include "TRIM2DS.h"
	#include "TRIMBC2DK.h"
	#include "Trimempk.h"
	#include "Trimm2dav.h"
	#include "Trimm2dk.h"
	#include "Trimm2dmf.h"
	#include "Trimm2ds.h"
	#include "Trimmpk.h"

int main(int argc, char *argv[])
{
/*____________________________________________________________________________

            						 DECLARATION OF REFINAMENT PARAMETERS
____________________________________________________________________________*/

 int  maxlev ;  /* maximum number of levels.                                */
 int   curlev;  /* current level in analysis.                               */
 int  curmesh;  /* current mesh in analysis.                                */
 int  newlev ;  /* new level to be defined.                                 */
 double  rtol;  /* error tolerance for mesh refinament.                     */
 int   refine;  /* global refinement indicator.                             */
 int  *levref;  /* level refinement indicator.                              */
 int   numpl ;  /* number of protective layers (default is one).            */
 int ESTIMATOR; /* estimator type 1:nodal avg. 2:midside stress             */
// int    *iptr;  /* auxiliar integer pointer for reallocation                */
// double *dptr;  /* auxiliar double  pointer for reallocation                */
 double Uerr ;  /* energy of error                                          */
 double Unew ;  /* energy of new displacements for all elements             */
 double Uexc ;  /* energy of exact solution if non zero value given         */
 double Ufea ;  /* energy of finite element approximation                   */
 double errm ;  /* mean error in energy norm of fine mesh                   */
 double rtolene; /* tolerance for convergence check in energy.              */
 double rtolest; /* tolerance for convergence check in error estimator.     */
 long  TRnumel; /* total number of element to be refined                    */

/*____________________________________________________________________________

            						 DECLARATION OF RELAXATION PARAMETERS
____________________________________________________________________________*/


 int  *iter   ;  /* iterative loop counter for PCG iterations in each level */
 int  SOLVER  ;  /* type of iterative solver                                */
                 /*    SOLVER = 1   hierarchical PCG                        */
                 /*    SOLVER = 2   diagonal PCG                            */
                 /*    SOLVER = 3   hierarchical MULTIGRID                  */
                 /*    SOLVER = 4   diagonal MULTIGRID                      */
 int  ver_ener;  /* convergence check on energy increment    (1:yes, 0:no)  */
 int  ver_res ;  /* convergence check on residual norm       (1:yes, 0:no)  */
 int  ver_est ;  /* convergence check on error estimator     (1:yes, 0:no)  */
 int  ver_nst ;  /* check if maximum number of steps reached (1:yes, 0:no)  */
 int  v_cyc ;  /* counter of multigrid cycles                             */

/*____________________________________________________________________________

          				DECLARATION OF COMMON PARAMETERS FOR ALL MESHES
____________________________________________________________________________*/

 int   nudim ;  /* number of spatial dimensions.                            */
 int   ecode ;  /* element code number.                                     */
 int   xndof ;  /* maximum number of nodal degrees of freedom.              */
 int   xnbdof;  /* maximum number of boundary degrees of freedom.           */
 int   xnsn  ;  /* maximum number of incident side nodes.                   */
 int   xnsv  ;  /* maximum number of side vertex nodes.                     */
 int   xnes  ;  /* maximum number of element sides.                         */
 int   xnen  ;  /* maximum number of incident element nodes.                */
 int   xnev  ;  /* maximum number of element vertex nodes.                  */
 int   xnumpr;  /* maximum number of material properties.                   */
 int   xdf   ;  /* maximum polynomial degree of external vector fields.     */
 int   numat ;  /* number of material sets.                                 */

/*____________________________________________________________________________

			DECLARATION OF ARRAYS OF MESH PARAMETERS FOR EACH LEVEL
____________________________________________________________________________*/

 int *nodes   ; /* number of nodal points.                                  */
 int *acnodes ; /* number of accumulated nodes in each level                */
 int *acnumel ; /* number of accumulated nodes in each level                */
 int *numel   ; /* number of elements.                                      */
 int *Rnumel   ; /* number of refinable elements.                           */
 int *Unumel   ; /* number of unrefinable elements.                         */
 int *FRnumel  ; /* number of refined elements.                             */
 int *URnumel  ; /* number of refinable elements not refined.               */
 int *Lnumel  ; /* number of loaded elements.                               */
 int *plnumel ; /* number of elements in the protective layers.             */
 int *sides   ; /* number of sides.                                         */
 int *Rsides   ; /* number of refinable sides.                              */
 int *Bsides   ; /* number of boundary sides.                               */
 int *BFRsides ; /* number of refined boundary sides.                       */
 int *BURsides ; /* number of boundary sides not refined.                   */
 int *FRsides  ; /* number of refined sides.                                */
 int *PRsides  ; /* number of transition sides.                             */
 int *URsides  ; /* number of refinable sides not refined.                  */
 int *Usides   ; /* number of unrefinable sides.                            */
 int *CFRsides ; /* number of new refined active sides.                     */
 int *CPRsides ; /* number of new refined transition sides.                 */
 int *CURsides ; /* number of new unrefined sides of group rside.           */
 int *CBFRsides; /* number of new refined active boundary sides.            */
 int *RBURsides; /* number of new unrefined sides of group rside.           */
 int *RURsides ; /* new number of refinable sides not refined after refin.  */
 int *RUsides  ; /* new number of unrefinable sides after refinement.       */
 int *NPRsides ; /* new number of transition sides after refinement.        */
 int *CFRnumel ; /* number of new refined elements after refinement.        */
 int *RURnumel ; /* number of unrefined sides of group rside.               */
 int *CURnumel ; /* number of unrefinable elem. changed to refinable.       */
 int *RUnumel  ; /* number of unrefinable elem. not changed.                */
 int *NURnumel ; /* number of new generated unrefinable elements.           */
 int *NUnumel  ; /* number of new generated refinable elements.             */
 int *NURsides ; /* number of new generated unrefinable sides.              */
 int *NUsides  ; /* number of new generated refinable sides.                */
 int *NBsides  ; /* number of new generated boundary sides.                 */
 int *Nnodes   ; /* number of new nodes.                                    */
 int *Nnumel   ; /* number of new elements.                                 */
 int *Nsides   ; /* number of new sides.                                    */

/*____________________________________________________________________________

				DECLARATION OF MESH ARRAYS - ONLY COARSE MESH (LEVEL 0)
____________________________________________________________________________*/

   int    *id     ;  /* identification vector [nodes[0] * xndof]              */
   int    *idb    ;  /* identification vector of boundary nodes [nodes[0]]    */
   int    *side_bc;  /* side boundary conditions [Bsides[0] * xnbdof]         */
   int    *side_mat; /* side material index of adjacent elements [2*sides[0]] */
   double *bn_cs  ;  /* boundary node cosine of local axes rotation angle     */
   double *bn_sn  ;  /* boundary node sine   of local axes rotation angle     */
   int    *mat    ;  /* element material type [numel[0]]                      */
   double *prop   ;  /* material properties   [numat * xnumpr]                */
   double *Gij    ;  /* element stiffness coefficients [10*numel[0]]          */
   double *Kij    ;  /* element stiffness coefficients [21*numel[0]]          */
   double *fp     ;  /* concentrated nodal loads coarse mesh only.            */
   int numcurv    ;  /* number of curved boundary sides initial mesh          */
   int    *nse    ;  /* number of elements at each node of coarse mesh        */
   double *ddiag  ;  /* diagonal coefficients matrix A [3*nodes[0]]           */

   double *fm     ;  /* concentrated mixed nodal loads coarse mesh only.      */

/*____________________________________________________________________________

		  DECLARATION OF SKYLINE PARAMETERS AND ARRAYS - ONLY COARSE MESH
____________________________________________________________________________*/

 int      neq ;  /* number of equations                                     */
 int      nwk ;  /* number of matrix elements                               */
 double  *sg  ;  /*  global stiffness matrix stored by skyline [nwk]        */
 int     *maxa;  /*  pointer of diagonal elements  for skyline [neq+1]      */
 int      nrw ;  /*  maximum number of rows in element matrix (xnen*xndof)  */
 int     *lm  ;  /*  identification vector of element dof [nrw]             */
 double  *sv  ;  /*  upper triangular element matrix [nrw*(nrw+1)/2]        */

/*____________________________________________________________________________

			 DECLARATION OF POINTERS OF MESH ARRAYS FOR EACH LEVEL
____________________________________________________________________________*/

   double **xyz      ; /* nodal coordinates     [nodes[curlev] * nudim]       */
   double **xyzc     ; /* coordinates of control point [Bsides[curlev]*nudim] */
   double **bs_cs    ; /* boundary side cosine of local system angle          */
    double **bs_sn    ; /* boundary side sine   of local system angle          */
   double **bn_nlength; /* boundary node summ of adjacent lengths normal      */
   double **bn_slength; /* boundary node summ of adjacent lengths tangent     */
   double **bn_qn    ; /* boundary node normal reaction                       */
   double **bn_qs    ; /* boundary node tangent reaction                      */
   int  **side_nodes ; /* nodes of each side       [sides[curlev] * xnsn]     */
   int  **elm_sides  ; /* sides of each element    [numel[level] * xnes]      */
   int  **elm_nodes  ; /* nodes/level of each element [numel[level] * xnes]   */
   int  **parent_elm ; /* parent element and location [3*numel[level]         */
   int  **parent_side; /* parent side and location    [2*Bsides[level]]       */
   int  **elm_vlev   ; /* levels of element vertex nodes  [xnev*numel[curlev]]*/
   int  **side_vlev  ; /* levels of side vertex nodes  [xnev*numel[curlev]]   */
   int  **velm       ; /* new element ordering after refinament [numel[level]]*/
   int  **vside      ; /* new side ordering after refinament [sides[level]]   */
   int  **eside      ; /* side refinament indicator and old side ordering     */
   int  **bside      ; /* list of boundary sides of each level.               */
   int  **tside      ; /* list of transition sides of each side.              */
   int  **ssside     ; /* son sides of each side. [sides[level]]              */
   int  **seside     ; /* son sides of each element. [numel[level]]           */
   int  **idnode     ; /* identification of nodes in refined region.          */
   int  **relm       ; /* element refinament indicator  [numel[level]]        */
   double **error    ; /* element error                 [numel[level]]        */

/*____________________________________________________________________________

	DECLARATION OF ARRAYS FOR DISPLACEMENTS AND FORCES - ONLY COARSE MESH
____________________________________________________________________________*/

 double **side_loads;/* boundary side loads [xnbdof][Bsides * (xdf+1)]      */
 double **elm_loads ;/* element body loads [xndof][Lnumel*(xdf+1)*(xdf+2)/2]*/
 double  *dn      ; /* vector of compacted nodal forces/displacements [neq] */

/*____________________________________________________________________________

	DECLARATION OF ARRAYS FOR DISPLACEMENTS AND FORCES FOR EACH LEVEL
____________________________________________________________________________*/

   double **f   ; /* vector of nodal forces             [nodes[level]*xndof]  */
   double **u1  ; /* vector of previous nodal displacements  [nodes[level]]   */
   double **diag; /* vector of diagonal matrix coefficients  [nodes[level]]   */
 //double ***mdiag;/*vector of diagonal stiffness coeff ALL MESHES AND LEVELS */
 //double ***res; /* vector of multigrid residuals ALL MESHES AND LEVELS      */
// double ***du ; /* vector of multigrid corrections ALL MESHES AND LEVELS    */
   double **p   ; /* auxiliar vector for conjugate gradients [nodes[level]]   */
   double **r   ; /* auxiliar vector for conjugate gradients [nodes[level]]   */
   double **z   ; /* auxiliar vector for conjugate gradients [nodes[level]]   */
   double **fi  ; /* auxiliar vector for internal forces     [nodes[level]]   */
   double **vec ; /* auxiliar vector for conjugate gradients [nodes[level]]   */
   double  **Sx  ; /* nodal stress Sx                           [nodes[level]] */
   double  **Sy  ; /* nodal stress Sy                           [nodes[level]] */
   double  **Sxy ; /* nodal stress Sxy                          [nodes[level]] */
   double **Fn  ; /* midside normal traction                 [sides[level]]   */
   double **Fs  ; /* midside tangential traction             [sides[level]]   */
   double **bs_Fn; /* midside normal boundary traction       [Bsides[level]]  */
   double **bs_Fs;/* midside tangential boundary traction   [Bsides[level]]   */
   double **a11 ; /* Z-Z coefficient matrix                                   */
   double **a12 ; /* Z-Z coefficient matrix                                   */
   double **a13 ; /* Z-Z coefficient matrix                                   */
   double **a22 ; /* Z-Z coefficient matrix                                   */
   double **a23 ; /* Z-Z coefficient matrix                                   */
   double **a33 ; /* Z-Z coefficient matrix                                   */
   double **b1x ; /* Z-Z right hand side component                            */
   double **b2x ; /* Z-Z right hand side component                            */
   double **b3x ; /* Z-Z right hand side component                            */
   double **b1y ; /* Z-Z right hand side component                            */
   double **b2y ; /* Z-Z right hand side component                            */
   double **b3y ; /* Z-Z right hand side component                            */
   double **b1xy; /* Z-Z right hand side component                            */
   double **b2xy; /* Z-Z right hand side component                            */
   double **b3xy; /* Z-Z right hand side component                            */

/*____________________________________________________________________________

							  AUXILIAR VARIABLES
____________________________________________________________________________*/


 int  i      ;  /* auxiliar loop counter for allocation                     */
 int  level  ;  /* auxiliar loop counter for levels                         */

/*____________________________________________________________________________

							DECLARATION OF FILES POINTERS
____________________________________________________________________________*/

 //char    filename[160];  /* input filename                                  */
 std::string filename;
 //FILE    *fpi          ;  /* pointer input  file                             */
 std::ifstream fpi           ;  /* pointer input  file                             */
 //char      forname[170]   ;  /* complete output file root name                  */
 std::string forname;
 //FILE    *fpo          ;  /* pointer results  file                           */
 std::ofstream fpo           ;  /* pointer results  file                           */
 //FILE    *fpiter       ;  /* pointer iteration steps file                    */
 std::ofstream fpiter        ;  /* pointer iteration steps file                    */

/*____________________________________________________________________________

							DECLARATION OF ALLOCATION VARIABLES
____________________________________________________________________________*/

 unsigned size;
 unsigned size1;
 unsigned size2;
 unsigned size3;

 long unsigned msize;
 long unsigned msize1;
 long unsigned msize2;
 long unsigned msize3;


/*____________________________________________________________________________

						  PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/
/*
  extern void  initp (void);

  extern void  iofile (char *filename, FILE **fpi, char *forname, FILE **fpo,
                       FILE **fpiter);

  extern void  inputpar (int *SOLVER, int ecode);

  extern void  inputest (int *ESTIMATOR, double *Uexc, int ecode);

  extern void  meshdat(FILE *fpi      , int nodes     , int numel,
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

  extern void  numeqn(int nodes, int xndof, int *id, int *neq);

  extern void  profile(int neq, int numel, int sides, int Bsides,
                       int ecode, int xndof, int *id, int *elm_sides,
                       int *side_nodes, int *lm, int *nwk, int *maxa);

  extern void  expandv2D (double *dn, double *u, int *id, int nodes,
                          int xndof, int *idb, double *bn_cs, double *bn_sn);

  extern void  expands (double *dn, double *u, int *id, int nodes);

  extern void  compressv2D (double *dn, double *u, int *id, int nodes,
                            int xndof, int *idb, double *bn_cs, double *bn_sn);

  extern void  compresss (double *dn, double *u, int *id, int nodes);

  extern void nodvec (FILE *fpi, int xndof, int *id, double *f);

  extern void  sidevec(FILE *fpi, int Bsides, int xdf, int xnbdof,double **side_loads);

  extern void  sideload(int ecode, int Bsides, int xdf, int xndof, int xnsn,
				        int *side_nodes, double **side_loads, double *bs_Fs,
                        double *bs_Fn, double *xyz,
						double *f, int *idb, double *bn_cs, double *bn_sn,
                        double *bs_cs, double *bs_sn, int *side_mat,
                        double *prop, int xnumpr, int nudim);

  extern void  mlsideload(int ecode, int *nodes, int *Bsides, int *BFRsides,
                          int xdf, int xndof,
                          int xnsn, int **bside, int **side_nodes,
                          int **side_vlev, double **side_loads,
                          int **parent_side, int curlev, double **xyz,
                          double **f, double **bs_cs, double **bs_sn);

  extern void side2dbt (int ecode , int *Bsides ,      int **bside,
                        int xdf   , double **side_loads, int **parent_side,
                        int curlev, double **bs_Fn     , double **bs_Fs);

  extern void  stiffblk(int ecode, int numel, int sides, int bsides, int xndof,
                        int nudim, int *id, int *elm_sides, int *side_nodes,
                        int *side_bc,int xnsn, double *xyz, int *side_mat,
		                int *mat, double *prop, int xnumpr, double *sg,
		                int *maxa, int *lm, double *sv, double *Gij,
                        double *Kij, double *diag, double *fp, double *fm);

  extern void  decomp (double *sg, int *maxa, int neq);

  extern void  redbak (double *sg, double *dn, int *maxa, int neq);

  extern void  prtplt (char *forname   , int curlev     , int curmesh,
                       int ecode,
                       int *nodes      , int *acnodes   , int *acnumel,
                       int *numel      , int *FRnumel   , int xndof,
                       int nudim       , int **elm_sides, int **elm_nodes,
                       int **elm_vlev,
                       int **parent_elm, int **velm     , double **xyz,
                       int *mat        , double *prop   , int xnumpr,
                       double **u      , double  **Fn    , double  **Fs,
                       double **Sx, double **Sy, double **Sxy,
                       double **error, double errm,
                       int ESTIMATOR, int SOLVER);

  extern void *mallocvf (unsigned siz);

  extern void *callocvf (unsigned nel, unsigned elsiz);

  extern void tr1_k0v (double **vec    , int *numel     , int *FRnumel,
                       int xndof       , int **elm_nodes, int **elm_vlev,
                       int **parent_elm, int **velm     , double **p,
                       double **xyz    , int *mat       , double *prop,
                       int xnumpr      , int curlev);


  extern void refidfr (int ecode      , int curlev,
                       int *acnumel, int *numel, int *nodes,
                       int *rnumel,
                       int *FRnumel, int **velm,  int xndof,
                       int nudim, int xnsn, double **xyz,
                       int **elm_sides , int **elm_nodes, int **elm_vlev,
                       int **parent_elm, int *sides    , int *Bsides,
                       int *BFRsides, int *FRsides, int *PRsides,
                       int **vside  ,  int **eside, int **bside,
                       int **tside, int **ssside,
                       int **side_nodes, int **side_vlev,
                       int **parent_side, int *side_bc, double **uold,
                       int *mat       , double *prop  , int xnumpr,
                       int *side_mat,
                       int **relm    , double **error,
                       double  **Fn   , double  **Fs,
                       double **Sx,  double **Sy,double **Sxy,
                       int *nse, int **idnode,
                       double **bs_Fn , double **bs_Fs,
                       double **a11, double **a12, double **a13,
                       double **a22, double **a23, double **a33,
                       double **b1x  , double **b2x  , double **b3x,
                       double **b1y  , double **b2y  , double **b3y,
                       double **b1xy , double **b2xy , double **b3xy,
                       int *refine,  int *levref, int ESTIMATOR,
                       long *TRnumel,
                       double *Uerr, double *Unew, double *Ufea,
                       double *errm , double **fi,
                       double **bn_nlength, double **bn_slength,
                       double **bn_qn, double **bn_qs,
                       double **bs_cs, double **bs_sn);

  extern void datalev (int curlev   , int newlev    ,
                       int *nodes   , int *numel    , int *Rnumel,
                       int *Unumel  , int *FRnumel  , int *URnumel,
                       int *Rsides  , int *Bsides   , int *BFRsides,
                       int *BURsides, int *FRsides  , int *PRsides,
                       int *URsides , int *Usides   ,
                       int **relm   , int **velm    , int **vside,
                       int **bside  , int **eside   , int **tside,
                       int **ssside , int **seside,
                       int *sides   , int **idnode  , int **parent_elm,
                       int **parent_side, int **elm_sides, int **elm_nodes,
                       int **elm_vlev, int **side_nodes, int **side_vlev,
                       int numpl,
                       int *plnumel , int *CFRsides , int *CPRsides,
                       int *CURsides, int *CBFRsides, int *RBURsides,
                       int *RURsides, int *RUsides  , int *NPRsides,
                       int *CFRnumel, int *CURnumel , int *RURnumel,
                       int *RUnumel , int *NURnumel , int *NUnumel,
                       int *NURsides, int *NUsides  , int *NBsides,
                       int *Nnodes  , int *Nnumel   , int *Nsides,
                       double **xyzc, double **bs_cs, double **bs_sn);

 extern void newcoor (int curlev  , int *nodes      , int *Nnodes,
                      int *BFRsides, int *CBFRsides ,
                      int **vside , int **eside     , int **bside,
                      int **side_nodes, int **side_vlev,
                      double **xyz, double **xyzc,
                      double **bs_cs, double **bs_sn,
                      int **parent_side, int *side_bc,
                      double **bn_nlength, double **bn_slength, int xndof);


 extern void update_DG (int curlev    , int *numel      , int *FRnumel,
                        int xndof     , int **elm_sides , int **elm_nodes,
                        int **elm_vlev, int **parent_elm, int **velm,
                        int **eside   , int *nodes      , double **diag,
                        double **xyz  , int *mat        , double *prop,
                        int xnumpr);

 extern void update_HD (int curlev    , int *FRnumel    , int *CFRnumel,
                        int xndof     , int **elm_sides , int **elm_nodes,
                        int **elm_vlev, int **parent_elm, int **velm,
                        int **eside   , double **diag   , double **xyz,
                        int *mat      , double *prop    , int xnumpr);


  extern void extrapol(int curlev     , int xndof     , int **parent_side,
                       int *BFRsides  , int *CBFRsides, int *FRsides,
                       int *CFRsides  , int *CPRsides , int **side_nodes,
                       int **side_vlev, int *side_bc  , int **vside,
                       int **eside    , int **bside   , double **uold,
                       double **bs_cs   , double **bs_sn);

//  extern void solver (int SOLVER,
//                      int curlev     , int *nodes       , int *numel,
//                      int *FRnumel   , int xndof        ,
//                      int **elm_nodes, int **elm_vlev   , int **parent_elm,
//                      int **velm     , int **parent_side, int *BFRsides,
//                      int *FRsides   , int *PRsides     , int **side_nodes,
//                      int **side_vlev, int *side_bc     , int **vside,
//                      double **uold  , double **f       , double **diag,
//                      double ***mdiag,
//                      double *dn     , int *id,
//                      double *sg     , int *maxa        , int neq,
//                      double **p, double **r, double **z, double **fi,
//                      double **vec,
//                      double ***res, double ***du,
//                      int **eside, int **tside, int **bside, int *idb,
//                      double *bn_cs, double *bn_sn,
//                      double **bs_cs, double **bs_sn,
//                      int *iter, int *v_cyc, double **xyz,
//                      int *mat, double *prop, int xnumpr,
//                      int ver_ener, int ver_res, int ver_est,
//                      int ver_nst,
//                      FILE *fpiter, double rtolene, double rtolest,
//                      double Uexc,
//                      int ecode      , int *acnumel  , int *rnumel,
//                      int **elm_sides, int *sides    , int *Bsides,
//                      int **ssside   , double  rtol  , int **relm,
//                      double **error , double  **Fn   , double  **Fs,
//                      double **Sx     , double **Sy    , double **Sxy,
//                      int *nse       , int **idnode,
//                      double **bs_Fn , double **bs_Fs,
//                      double **a11   , double **a12  , double **a13,
//                      double **a22   , double **a23  , double **a33,
//                      double **b1x   , double **b2x  , double **b3x,
//                      double **b1y   , double **b2y  , double **b3y,
//                      double **b1xy  , double **b2xy , double **b3xy,
//                      int *levref   , int ESTIMATOR,
//                      double **bn_nlength, double **bn_slength,
//                      double **bn_qn, double **bn_qs);

*/

/*____________________________________________________________________________

							   VERIFY NUMBER OF COMMAND LINE ARGUMENTS
____________________________________________________________________________*/

/*
 if (argc==1)
   {
				printf ("\n");
				printf ("Usage: ");
				printf ("Samfe2d <Filename>");
				printf ("\n");
    exit(0);
   };
*/


 //strncpy(filename, argv[1], strlen(argv[1])+1);
 filename=std::string(argv[1]);

/*____________________________________________________________________________

							   DISPLAY PRESENTATION SCREEN
____________________________________________________________________________*/


  //initp ();

/*____________________________________________________________________________

						  	  OPEN INPUT/OUTPUT FILES
____________________________________________________________________________*/

    iofile (&filename, fpi, &forname, fpo, fpiter);
  //iofile (filename, &fpi, forname, &fpo, &fpiter);

/*____________________________________________________________________________

					 ALLOCATION FOR ARRAYS OF MESH PARAMETERS
____________________________________________________________________________*/

 ver_nst  = 0;
 maxlev = 1;


//  nodes     = (int *) callocvf ((maxlev+1) , sizeof (int));
    nodes=new int[maxlev+1]; for (int i = 0; i < maxlev + 1; i++)nodes[i] = 0;
//  acnodes   = (int *) callocvf ((maxlev+1) , sizeof (int));
    acnodes=new int[maxlev+1]; for (int i = 0; i < maxlev + 1; i++)acnodes[i] = 0;
//  acnumel   = (int *) callocvf ((maxlev+1) , sizeof (int));
    acnumel=new int[maxlev+1]; for (int i = 0; i < maxlev + 1; i++)acnumel[i] = 0;
//  numel     = (int *) callocvf ((maxlev+1) , sizeof (int));
    numel=new int[maxlev+1]; for (int i = 0; i < maxlev + 1; i++)numel[i] = 0;
//  Rnumel    = (int *) callocvf ((maxlev+1) , sizeof (int));
    Rnumel=new int[maxlev+1]; for (int i = 0; i < maxlev + 1; i++)Rnumel[i] = 0;
//  Unumel    = (int *) callocvf ((maxlev+1) , sizeof (int));
    Unumel=new int[maxlev+1]; for (int i = 0; i < maxlev + 1; i++)Unumel[i] = 0;
//  FRnumel   = (int *) callocvf ((maxlev+1) , sizeof (int));
    FRnumel=new int[maxlev+1]; for (int i = 0; i < maxlev + 1; i++)FRnumel[i] = 0;
//  URnumel   = (int *) callocvf ((maxlev+1) , sizeof (int));
    URnumel=new int[maxlev+1]; for (int i = 0; i < maxlev + 1; i++)URnumel[i] = 0;
//  Lnumel    = (int *) callocvf ((maxlev+1) , sizeof (int));
    Lnumel=new int[maxlev+1]; for (int i = 0; i < maxlev + 1; i++)Lnumel[i] = 0;
//  plnumel   = (int *) callocvf ((maxlev+1) , sizeof (int));
    plnumel=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)plnumel[i] = 0;
//  sides     = (int *) callocvf ((maxlev+1) , sizeof (int));
    sides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)sides[i] = 0;
//  Rsides    = (int *) callocvf ((maxlev+1) , sizeof (int));
    Rsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)Rsides[i] = 0;
//  Bsides    = (int *) callocvf ((maxlev+1) , sizeof (int));
    Bsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)Bsides[i] = 0;
//  BFRsides  = (int *) callocvf ((maxlev+1) , sizeof (int));
    BFRsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)BFRsides[i] = 0;
//  BURsides  = (int *) callocvf ((maxlev+1) , sizeof (int));
    BURsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)BURsides[i] = 0;
//  FRsides   = (int *) callocvf ((maxlev+1) , sizeof (int));
    FRsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)FRsides[i] = 0;
//  PRsides   = (int *) callocvf ((maxlev+1) , sizeof (int));
    PRsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)PRsides[i] = 0;
//  URsides   = (int *) callocvf ((maxlev+1) , sizeof (int));
    URsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)URsides[i] = 0;
//  Usides    = (int *) callocvf ((maxlev+1) , sizeof (int));
    Usides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)Usides[i] = 0;
//  CFRsides  = (int *) callocvf ((maxlev+1) , sizeof (int));
    CFRsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)CFRsides[i] = 0;
//  CPRsides  = (int *) callocvf ((maxlev+1) , sizeof (int));
    CPRsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)CPRsides[i] = 0;
//  CURsides  = (int *) callocvf ((maxlev+1) , sizeof (int));
    CURsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)CURsides[i] = 0;
//  CBFRsides = (int *) callocvf ((maxlev+1) , sizeof (int));
    CBFRsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)CBFRsides[i] = 0;
//  RBURsides = (int *) callocvf ((maxlev+1) , sizeof (int));
    RBURsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)RBURsides[i] = 0;
//  RURsides  = (int *) callocvf ((maxlev+1) , sizeof (int));
    RURsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)RURsides[i] = 0;
//  RUsides   = (int *) callocvf ((maxlev+1) , sizeof (int));
    RUsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)RUsides[i] = 0;
//  NPRsides  = (int *) callocvf ((maxlev+1) , sizeof (int));
    NPRsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)NPRsides[i] = 0;
//  CFRnumel  = (int *) callocvf ((maxlev+1) , sizeof (int));
    CFRnumel=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)CFRnumel[i] = 0;
//  RURnumel  = (int *) callocvf ((maxlev+1) , sizeof (int));
    RURnumel=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)RURnumel[i] = 0;
//  CURnumel  = (int *) callocvf ((maxlev+1) , sizeof (int));
    CURnumel=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)CURnumel[i] = 0;
//  RUnumel   = (int *) callocvf ((maxlev+1) , sizeof (int));
    RUnumel=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)RUnumel[i] = 0;
//  NURnumel  = (int *) callocvf ((maxlev+1) , sizeof (int));
    NURnumel=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)NURnumel[i] = 0;
//  NUnumel   = (int *) callocvf ((maxlev+1) , sizeof (int));
    NUnumel=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)NUnumel[i] = 0;
//  NURsides  = (int *) callocvf ((maxlev+1) , sizeof (int));
    NURsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)NURsides[i] = 0;
//  NUsides   = (int *) callocvf ((maxlev+1) , sizeof (int));
    NUsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)NUsides[i] = 0;
//  NBsides   = (int *) callocvf ((maxlev+1) , sizeof (int));
    NBsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)NBsides[i] = 0;
//  Nnodes    = (int *) callocvf ((maxlev+1) , sizeof (int));
    Nnodes=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)Nnodes[i] = 0;
//  Nnumel    = (int *) callocvf ((maxlev+1) , sizeof (int));
    Nnumel=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)Nnumel[i] = 0;
//  Nsides    = (int *) callocvf ((maxlev+1) , sizeof (int));
    Nsides=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)Nsides[i] = 0;
//  iter     = (int *) callocvf ((maxlev+1) , sizeof (int));
    iter=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)iter[i] = 0;
//  levref   = (int *) callocvf ((maxlev+1) , sizeof (int));
    levref=new int [maxlev+1]; for (int i = 0; i < maxlev + 1; i++)levref[i] = 0;

/*____________________________________________________________________________

					 ALLOCATION FOR POINTERS OF MESH ARRAYS
____________________________________________________________________________*/

//  side_nodes  = (int **) callocvf ((maxlev+1) , sizeof (int *));
    side_nodes  = new int* [maxlev+1];
//  side_vlev   = (int **) callocvf ((maxlev+1) , sizeof (int *));
    side_vlev   = new int* [maxlev+1];
//  elm_sides   = (int **) callocvf ((maxlev+1) , sizeof (int *));
    elm_sides   = new int* [maxlev+1];
//  elm_nodes   = (int **) callocvf ((maxlev+1) , sizeof (int *));
    elm_nodes   = new int* [maxlev+1];
//  elm_vlev    = (int **) callocvf ((maxlev+1) , sizeof (int *));
    elm_vlev    = new int* [maxlev+1];
//  parent_elm  = (int **) callocvf ((maxlev+1) , sizeof (int *));
    parent_elm  = new int* [maxlev+1];
//  parent_side = (int **) callocvf ((maxlev+1) , sizeof (int *));
    parent_side = new int* [maxlev+1];
//  velm        = (int **) callocvf ((maxlev+1) , sizeof (int *));
    velm        = new int* [maxlev+1];
//  vside       = (int **) callocvf ((maxlev+1) , sizeof (int *));
    vside       = new int* [maxlev+1];
//  eside       = (int **) callocvf ((maxlev+1) , sizeof (int *));
    eside       = new int* [maxlev+1];
//  bside       = (int **) callocvf ((maxlev+1) , sizeof (int *));
    bside       = new int* [maxlev+1];
//  tside       = (int **) callocvf ((maxlev+1) , sizeof (int *));
    tside       = new int* [maxlev+1];
//  ssside      = (int **) callocvf ((maxlev+1) , sizeof (int *));
    ssside      = new int* [maxlev+1];
//  seside      = (int **) callocvf ((maxlev+1) , sizeof (int *));
    seside      = new int* [maxlev+1];
//  idnode      = (int **) callocvf ((maxlev+1) , sizeof (int *));
    idnode      = new int* [maxlev+1];
//  relm        = (int **) callocvf ((maxlev+1) , sizeof (int *));
    relm        = new int* [maxlev+1];



  //xyz    = (double **) callocvf ((maxlev+1) , sizeof (double *));
    xyz    = new double* [maxlev+1];

  //xyzc   = (double **) callocvf ((maxlev+1) , sizeof (double *));
    xyzc   = new double* [maxlev+1];

  //bs_cs  = (double **) callocvf ((maxlev+1) , sizeof (double *));
	bs_cs = new double* [maxlev + 1];

  //bs_sn  = (double **) callocvf ((maxlev+1) , sizeof (double *));
    bs_sn  = new double* [maxlev+1];

  //bn_nlength = (double **) callocvf ((maxlev+1) , sizeof (double *));
    bn_nlength = new double* [maxlev+1];

  //bn_slength = (double **) callocvf ((maxlev+1) , sizeof (double *));
    bn_slength = new double* [maxlev+1];

  //bn_qn  = (double **) callocvf ((maxlev+1) , sizeof (double *));
    bn_qn = new double* [maxlev+1];

  //bn_qs  = (double **) callocvf ((maxlev+1) , sizeof (double *));
   bn_qs  =  new double* [maxlev+1];


  //f      = (double **) callocvf ((maxlev+1) , sizeof (double *));
    f      = new double* [maxlev+1];

  //u1     = (double **) callocvf ((maxlev+1) , sizeof (double *));
    u1     = new double* [maxlev+1];

  //diag   = (double **) callocvf ((maxlev+1) , sizeof (double *));
    diag   = new double* [maxlev+1];

  //p      = (double **) callocvf ((maxlev+1) , sizeof (double *));
    p      = new double* [maxlev+1];

  //r      = (double **) callocvf ((maxlev+1) , sizeof (double *));
    r      =  new double* [maxlev+1];

  //z      = (double **) callocvf ((maxlev+1) , sizeof (double *));
    z      = new double* [maxlev+1];

  //fi     = (double **) callocvf ((maxlev+1) , sizeof (double *));
    fi     = new double* [maxlev+1];

  //vec    = (double **) callocvf ((maxlev+1) , sizeof (double *));
    vec    = new double* [maxlev+1];

  //Sx     = (double  **) callocvf ((maxlev+1) , sizeof (double  *));
   Sx     =  new double* [maxlev+1];

  //Sy     = (double  **) callocvf ((maxlev+1) , sizeof (double  *));
    Sy     = new double* [maxlev+1];

  //Sxy    = (double  **) callocvf ((maxlev+1) , sizeof (double  *));
    Sxy    = new double* [maxlev+1];

  //Fn     = (double  **) callocvf ((maxlev+1) , sizeof (double  *));
    Fn     = new double* [maxlev+1];

  //Fs     = (double  **) callocvf ((maxlev+1) , sizeof (double  *));
    Fs     = new double* [maxlev+1];

  //bs_Fn  = (double **) callocvf ((maxlev+1) , sizeof (double *));
    bs_Fn  = new double* [maxlev+1];

  //bs_Fs  = (double **) callocvf ((maxlev+1) , sizeof (double *));
    bs_Fs  = new double* [maxlev+1];

  //a11    = (double **) callocvf ((maxlev+1) , sizeof (double *));
    a11    = new double* [maxlev+1];

  //a12    = (double **) callocvf ((maxlev+1) , sizeof (double *));
    a12    =  new double* [maxlev+1];

  //a13    = (double **) callocvf ((maxlev+1) , sizeof (double *));
    a13    = new double* [maxlev+1];

  //a22    = (double **) callocvf ((maxlev+1) , sizeof (double *));
    a22    = new double* [maxlev+1];

  //a23    = (double **) callocvf ((maxlev+1) , sizeof (double *));
    a23    = new double* [maxlev+1];

  //a33    = (double **) callocvf ((maxlev+1) , sizeof (double *));
    a33    = new double* [maxlev+1];

  //b1x    = (double **) callocvf ((maxlev+1) , sizeof (double *));
    b1x    = new double* [maxlev+1];

  //b2x    = (double **) callocvf ((maxlev+1) , sizeof (double *));
    b2x    = new double* [maxlev+1];

  //b3x    = (double **) callocvf ((maxlev+1) , sizeof (double *));
    b3x    = new double* [maxlev+1];

  //b1y    = (double **) callocvf ((maxlev+1) , sizeof (double *));
    b1y    = new double* [maxlev+1];

  //b2y    = (double **) callocvf ((maxlev+1) , sizeof (double *));
    b2y    = new double* [maxlev+1];

  //b3y    = (double **) callocvf ((maxlev+1) , sizeof (double *));
    b3y    = new double* [maxlev+1];

  //b1xy   = (double **) callocvf ((maxlev+1) , sizeof (double *));
    b1xy   = new double* [maxlev+1];

  //b2xy   = (double **) callocvf ((maxlev+1) , sizeof (double *));
    b2xy   = new double* [maxlev+1];

  //b3xy   = (double **) callocvf ((maxlev+1) , sizeof (double *));
    b3xy   = new double* [maxlev+1];

  //error  = (double **) callocvf ((maxlev+1) , sizeof (double *));
    error  = new double* [maxlev+1];


//  mdiag = (double ***) callocvf((maxlev+1), sizeof(double **));
//  double **mdiag [maxlev+1];
//  mdiag[0]= (double **) callocvf((maxlev+1)*(maxlev+2)/2, sizeof(double **));
//  for (i=1; i<=maxlev ; i++) mdiag[i] = mdiag[i-1] + i;

//  du = (double ***) callocvf((maxlev+1), sizeof(double **));
//  du[0]= (double **) callocvf((maxlev+1)*(maxlev+2)/2, sizeof(double **));
//  for (i=1; i<=maxlev ; i++) du[i] = du[i-1] + i;

//  res = (double ***) callocvf((maxlev+1), sizeof(double **));
//  res[0]= (double **) callocvf((maxlev+1)*(maxlev+2)/2, sizeof(double **));
//  for (i=1; i<=maxlev ; i++) res[i] = res[i-1] + i;


/*____________________________________________________________________________

									GLOBAL COARSE MESH DATA
____________________________________________________________________________*/


  //fscanf(fpi," %d %d %d %d %d", &nodes[0], &numel[0], &sides[0] , &Bsides[0],
  //                              &numcurv);
  fpi>>nodes[0]>>numel[0]>>sides[0]>>Bsides[0]>>numcurv;

  //fscanf(fpi," %d %d %d %d", &ecode   , &numat   , &Lnumel[0], &xdf      );
  fpi>>ecode>>numat>>Lnumel[0]>>xdf;

/*____________________________________________________________________________

								 ELEMENT PARAMETERS
____________________________________________________________________________*/

  switch (ecode)
		  {
			 case 1: /* TRI3 plane linear triangle */
				nudim = 2;
				xndof = 2;
				xnbdof= 2;
				xnsn  = 4;
				xnsv  = 4;
				xnes  = 3;
				xnen  = 3;
				xnev  = 3;
				xnumpr= 3;
			 break;

			 case 2: /* plane quadratic triangle - right sided */
				nudim = 2;
				xndof = 2;
				xnbdof= 2;
				xnsn  = 3;
				xnsv  = 2;
				xnes  = 3;
				xnen  = 6;
				xnev  = 6;
				xnumpr= 3;
			 break;

			 case 3: /* plane linear mixed triangle */
				nudim = 2;
				xndof = 2;
				xnbdof= 2;
				xnsn  = 4;
				xnsv  = 4;
				xnes  = 3;
				xnen  = 3;
				xnev  = 3;
				xnumpr= 3;
			 break;

			 case 4: /* Supeconvergent Morley triangle */
				nudim = 2;
				xndof = 1;
				xnbdof= 2;
				xnsn  = 9;
				xnsv  = 4;
				xnes  = 3;
				xnen  = 6;
				xnev  = 6;
				xnumpr= 3;
			 break;

             case 5: /* Morley triangle by side */
				nudim = 2;
				xndof = 1;
				xnbdof= 2;
				xnsn  = 9;
				xnsv  = 4;
				xnes  = 3;
				xnen  = 6;
				xnev  = 6;
				xnumpr= 3;
			 break;

			 case 6: /* Morley plate triangle */
				nudim = 2;
				xndof = 1;
				xnbdof= 2;
				xnsn  = 4;
				xnsv  = 4;
				xnes  = 3;
				xnen  = 3;
				xnev  = 3;
				xnumpr= 3;
			 break;


			 case 8: /* TRI3BC plane linear triangle - cubic Bézier shape */
				 nudim = 3;
				 xndof = 2;
				 xnbdof= 2;
				 xnsn  = 6;
				 xnsv  = 2;
				 xnes  = 3;
				 xnen  = 9;
				 xnev  = 3;
				 xnumpr= 3;
             break;

			 default:
                /*
				printf ("\n");
				printf ("FATAL ERROR: ");
				printf ("ELEMENT CODE NUMBER NOT RECOGNIZED");
				printf ("\n");
				*/
				std::cout<<std::endl;
				std::cout<<"FATAL ERROR: ";
				std::cout<<"ELEMENT CODE NUMBER NOT RECOGNIZED";
				std::cout<<std::endl;


				exit(1);

		  }

/*____________________________________________________________________________

			  ALLOCATION FOR GLOBAL VARIABLES OF COARSE MESH
____________________________________________________________________________*/


//  xyz[0]= (double *) callocvf ((nodes[0] * nudim  + 1) , sizeof (double));
		  xyz[0] = new double[nodes[0] * nudim + 1]; for (int i = 0; i < nodes[0] * nudim + 1; i++)xyz[0][i] = 0.0;


  //xyzc[0]=(double *) callocvf ((Bsides[0] * nudim  + 1) , sizeof (double));
   xyzc[0]= new double [Bsides[0] * nudim  + 1]; for (int i = 0; i < Bsides[0] * nudim + 1; i++)xyzc[0][i] = 0.0;

  //id    = (int    *) callocvf ((nodes[0] * xndof  + 1) , sizeof (int));
   id     = new int [nodes[0] * xndof  + 1]; for (int i = 0; i < nodes[0] * xndof + 1; i++)id[i] = 0;

  //mat   = (int    *) callocvf ((numel[0] + 1         ) , sizeof (int));
    mat   = new int [numel[0] + 1]; for (int i = 0; i < numel[0] + 1; i++)mat[i] = 0;

  //prop  = (double *) callocvf ((numat * xnumpr + 1) , sizeof (double));
    prop  = new double [numat * xnumpr + 1]; for (int i = 0; i < numat * xnumpr + 1; i++)prop[i] = 0.0;


  //side_bc = (int *) callocvf ((Bsides[0] * xnbdof + 1) , sizeof(int));
    side_bc = new int [Bsides[0] * xnbdof + 1];

  //side_mat      = (int *) callocvf ((sides[0] * 2 + 1), sizeof(int));
    side_mat      = new int [sides[0] * 2 + 1];

  //bs_cs[0]      = (double *) callocvf ((Bsides[0] + 1  ) , sizeof(double));
    bs_cs[0]      = new double [Bsides[0] + 1]; for (int i = 0; i < Bsides[0] + 1; i++) bs_cs[0][i] = 0.0;

  //bs_sn[0]      = (double *) callocvf ((Bsides[0] + 1  ) , sizeof(double));
    bs_sn[0]      = new double [Bsides[0] + 1]; for (int i = 0; i < Bsides[0] + 1; i++) bs_sn[0][i] = 0.0;

  //bn_nlength[0] = (double *) callocvf ((nodes[0] + 1  ) , sizeof(double));
    bn_nlength[0] = new double [nodes[0] + 1]; for (int i = 0; i < nodes[0] + 1; i++) bn_nlength[0][i] = 0.0;

  //bn_slength[0] = (double *) callocvf ((nodes[0] + 1  ) , sizeof(double));
    bn_slength[0] = new double [nodes[0] + 1]; for (int i = 0; i < nodes[0] + 1; i++) bn_slength[0][i] = 0.0;


  //idb           = (int *) callocvf ((nodes[0] + 1         ) , sizeof(int));
	idb = new int[nodes[0] + 1]; for (int i = 0; i < nodes[0] + 1; i++)idb[i]=0;


  //bn_cs = (double *)callocvf((nodes[0] + 1), sizeof(double));
   bn_cs = new double[nodes[0] + 1]; for (int i = 0; i < nodes[0] + 1; i++) bn_cs[i] = 0.0;

  //bn_sn = (double *)callocvf((nodes[0] + 1), sizeof(double));
    bn_sn = new double [nodes[0] + 1]; for (int i = 0; i < nodes[0] + 1; i++) bn_sn[i] = 0.0;


  //side_nodes[0] = (int *) callocvf ((sides[0]  * xnsn  + 1) , sizeof(int));
    side_nodes[0] =  new int [sides[0]  * xnsn  + 1]; for (int i = 0; i < sides[0] * xnsn + 1; i++) side_nodes[0][i] = 0;


  //side_vlev [0] = (int *) callocvf ((sides[0]  * xnsv  + 1) , sizeof(int));
    side_vlev [0] = new int [sides[0]  * xnsv  + 1]; for (int i = 0; i < sides[0] * xnsv + 1; i++) side_vlev[0][i] = 0;


  //elm_sides [0] = (int *) callocvf ((numel[0]  * xnes  + 1) , sizeof(int));
    elm_sides [0] =  new int [numel[0]  * xnes  + 1]; for (int i = 0; i < numel[0] * xnes + 1; i++) elm_sides[0][i] = 0;


  //elm_nodes [0] = (int *) callocvf ((numel[0]  * xnen  + 1) , sizeof(int));
    elm_nodes [0] =  new int [numel[0]  * xnen  + 1]; for (int i = 0; i < numel[0] * xnen + 1; i++) elm_nodes[0][i] = 0;


  //elm_vlev  [0] = (int *) callocvf ((numel[0]  * xnev  + 1) , sizeof(int));
    elm_vlev  [0] =  new int [numel[0]  * xnev  + 1]; for (int i = 0; i < numel[0] * xnev + 1; i++) elm_vlev[0][i] = 0;


  //idnode    [0] = (int *) callocvf ((nodes[0] + 1)          , sizeof(int));
    idnode    [0] =  new int [nodes[0] + 1]; for (int i = 0; i < nodes[0] + 1; i++) idnode[0][i] = 0;



  //Gij  = (double *) callocvf ((10*numel[0]+1), sizeof (double));
    Gij  = new double [10*numel[0]+1]; for (int i = 0; i < 10 * numel[0] + 1; i++) Gij[i] = 0.0;


  //Kij  = (double *) callocvf ((21*numel[0]+1), sizeof (double));
    Kij  = new double [21*numel[0]+1]; for (int i = 0; i < 21 * numel[0] + 1; i++) Kij[i] = 0.0;

  //fp   = (double *) callocvf ((xndof*nodes[0]+1), sizeof (double));
    fp   = new double [xndof*nodes[0]+1]; for (int i = 0; i < xndof * nodes[0] + 1; i++) fp[i] = 0.0;

  //fm   = (double *) callocvf ((xndof*nodes[0]+1), sizeof (double));
    fm   = new double [xndof*nodes[0]+1]; for (int i = 0; i < xndof * nodes[0] + 1; i++) fm[i] = 0.0;


/*  if (ecode==3)
    {
     ddiag  = (double *) callocvf ((3*nodes[0]+1), sizeof (double));
     C1ij  = (double *) callocvf ((6*sides[0]+1), sizeof (double));
     C2ij  = (double *) callocvf ((6*sides[0]+1), sizeof (double));
     C3ij  = (double *) callocvf ((6*sides[0]+1), sizeof (double));
     C4ij  = (double *) callocvf ((6*sides[0]+1), sizeof (double));
    }
*/

/*____________________________________________________________________________

                    INPUT RELAXATION PARAMETERS
____________________________________________________________________________*/

 //inputpar (&SOLVER, ecode);
 SOLVER=1;
 //inputest (&ESTIMATOR, &Uexc, ecode);
 ESTIMATOR=1;
 Uexc=1.0;




/*____________________________________________________________________________

                    OUTPUT RELAXATION PARAMETERS
____________________________________________________________________________*/

 //fprintf (fpo,"\n\n TYPE OF ANALYSIS: ADAPATIVE REFINEMENT");
 fpo<<"\n\n TYPE OF ANALYSIS: ADAPATIVE REFINEMENT"<<std::endl;
 //fprintf (fpo,"\n");

   switch (SOLVER)
     {
       case 1:
        //fprintf (fpo,"\n SOLVER TYPE : GAUSS ELIMINATION");
        fpo<<"\n SOLVER TYPE : GAUSS ELIMINATION";
       break;

       case 2:
        //fprintf (fpo,"\n SOLVER TYPE : ORTHOGONAL DECOMPOSITION ");
        fpo<<"\n SOLVER TYPE : ORTHOGONAL DECOMPOSITION ";
       break;
     }

    //fprintf (fpo,"\n");
    fpo<<std::endl;
     //fprintf (fpo,"\n ESTIMATOR TYPE                 :");
       fpo<<"\n ESTIMATOR TYPE                 :";
 switch (ESTIMATOR)
   {
     case 1:
      //fprintf (fpo," NODAL AVERAGING");
      fpo<<" NODAL AVERAGING";
     break;

     case 2:
      //fprintf (fpo," MIDSIDE AVERAGING");
      fpo<<" MIDSIDE AVERAGING";
     break;

     case 3:
      //fprintf (fpo," ZZ NODAL RECOVERY CENTROIDAL SAMPLING");
      fpo<<" ZZ NODAL RECOVERY CENTROIDAL SAMPLING";
     break;

     case 4:
      //fprintf (fpo," ZZ NODAL RECOVERY SIDE SAMPLING");
      fpo<<" ZZ NODAL RECOVERY SIDE SAMPLING";
     break;

     case 5:
      //fprintf (fpo," SIMPLE SIDE AVERAGING (WITHOUT EXTERNAL FLUXES)");
      fpo<<" SIMPLE SIDE AVERAGING (WITHOUT EXTERNAL FLUXES)";
     break;
   }

   if (Uexc)
      //fprintf (fpo,"\n\n EXACT SOLUTION ENERGY NORM: %lf ", Uexc);
      fpo<<"\n\n EXACT SOLUTION ENERGY NORM: %lf "<<Uexc;

   //fprintf(fpo,"\n\n MESH  NUMEL  DOF");
    fpo<<"\n\n MESH  NUMEL  DOF";

   if (Uexc)
   {
        //fprintf(fpo,"       Uh           Eh           Eex     THETA  Err* %%  ErrEX %%");
        fpo<<"       Uh           Eh           Eex     THETA  Err* %%  ErrEX %%";
   }



/*____________________________________________________________________________

									                  COARSE MESH DATA
____________________________________________________________________________*/

  meshdat (fpi, nodes[0], numel[0], sides[0], Bsides[0], nudim, ecode,
           xyz[0], xyzc[0], id, xndof, xnbdof, xnsn, xnes, side_nodes[0],
           side_vlev[0], elm_sides[0], elm_nodes[0], elm_vlev[0],
           side_bc, side_mat, numat,
           mat, prop, xnumpr, idb, bn_cs, bn_sn, numcurv,
           bs_cs[0], bs_sn[0], bn_nlength[0], bn_slength[0]);

/*____________________________________________________________________________

							  COMPUTES EQUATIONS NUMBERS
____________________________________________________________________________*/


  numeqn (nodes[0], xndof, id, &neq);


/*____________________________________________________________________________

	MEMORY ALLOCATION FOR SKYLINE POINTERS AND ELEMENT MATRIX ARRAYS
____________________________________________________________________________*/

  switch (ecode)
  {
  case 1: /* plane linear triangle */
    nrw = xnen * xndof;
    break;

  case 2: /* plane quadratic triangle - right sided */
    nrw = xnen * xndof;
    break;

  case 3: /* plane linear mixed triangle */
    nrw = xnsn * xndof;
    break;

  case 4: /* Supeconvergent Morley triangle */
    nrw = xnsn * xndof;
    break;

  case 5: /* Morley triangle by side */
    nrw = xnsn * xndof;
    break;

  case 7: /* TRI3 plane linear triangle by side */
    nrw = xnsn * xndof;
	break;

  case 8: /* TRI3BC plane linear triangle - cubic Bézier shape */
    nrw = xnen * xndof;
    break;

  default:
                /*
				printf ("\n");
				printf ("FATAL ERROR: ");
				printf ("ELEMENT CODE NUMBER NOT RECOGNIZED");
				printf ("\n");
				*/
				std::cout<<std::endl;
				std::cout<<"FATAL ERROR: ";
				std::cout<<"ELEMENT CODE NUMBER NOT RECOGNIZED";
				std::cout<<std::endl;
				exit(1);
    }

  //maxa = (int    *) callocvf ((neq + 2)            , sizeof (int));
  maxa = new int[neq + 2]; for (int i = 0; i < neq + 2; i++)  maxa[i] = 0;

  //lm   = (int    *) callocvf ((nrw + 1)            , sizeof (int));
  lm = new int[nrw + 1]; for (int i = 0; i < nrw + 1; i++)  lm[i] = 0;

  //sv   = (double *) callocvf ((nrw * (nrw+1)/2 + 1), sizeof (double));
  sv = new double[nrw * (nrw+1)/2 + 1]; for (int i = 0; i < nrw * (nrw + 1) / 2 + 1; i++)  sv[i] = 0.0;


/*____________________________________________________________________________

					  COMPUTES SKYLINE SIZE AND POINTERS
____________________________________________________________________________*/


  profile (neq, numel[0], sides[0], Bsides[0], ecode, xndof, id, elm_sides[0],
           side_nodes[0], lm, &nwk, maxa);

/*____________________________________________________________________________

				  MEMORY ALLOCATION FOR GLOBAL STIFFNESS MATRIX
____________________________________________________________________________*/


  //sg = (double *) callocvf ((nwk + 1), sizeof (double));
  sg = new double[nwk + 1]; for (int i = 0; i < nwk + 1; i++)  sv[i] = 0.0;

/*____________________________________________________________________________

				 MEMORY ALLOCATION FOR FORCES AND DISPLACEMENTS
____________________________________________________________________________*/


//  side_loads = (double **) callocvf (xndof, sizeof (double *));
side_loads = new double*[xndof];
  //side_loads[0] = (double *)  callocvf ((xndof * (Bsides[0] * (xdf+1)+1)),sizeof (double));
side_loads[0] = new double [xndof * (Bsides[0] * (xdf+1)+1)]; for (int i = 0; i <xndof * (Bsides[0] * (xdf + 1) + 1); i++)  side_loads[0][i] = 0.0;
//  for (i=1 ; i<=xndof-1 ; i++)
//				  side_loads[i] = side_loads[0] + i * (Bsides[0]*(xdf+1)+1);
for (i=1 ; i<=xndof-1 ; i++)  side_loads[i] = side_loads[0] + i * (Bsides[0]*(xdf+1)+1);




  //bs_Fn[0] = (double *) callocvf ((Bsides[0]+1), sizeof (double));
    bs_Fn[0] = new double [Bsides[0]+1];


  //bs_Fs[0] = (double *) callocvf ((Bsides[0]+1), sizeof (double));
    bs_Fs[0] = new double [Bsides[0]+1];


  if (Lnumel[0])
	 {
// 		elm_loads    = (double **) callocvf (xndof, sizeof (double *));
//	 	elm_loads[0] = (double *)callocvf (xndof * (Lnumel[0]*(xdf+1)*(xdf+2)/2+1), sizeof (double));
//		for (i=1 ; i<=xndof-1 ; i++)
//				 elm_loads[i] = elm_loads[0] + i*(Lnumel[0]*(xdf+1)*(xdf+2)/2+1);
	  elm_loads = new double* [xndof];
	  elm_loads[0]=new double [Lnumel[0]*(xdf+1)*(xdf+2)/2+1]; for (int i = 0; i < Lnumel[0] * (xdf + 1) * (xdf + 2) / 2 + 1; i++)elm_loads[0][i] = 0.0;
	  for (i=1 ; i<=xndof-1 ; i++) elm_loads[i] = elm_loads[0] + i*(Lnumel[0]*(xdf+1)*(xdf+2)/2+1);
	 }

  //dn = (double *) callocvf ((neq   + 1), sizeof (double));
    dn = new double [neq   + 1]; for (int i = 0; i < neq + 1; i++)  dn[i] = 0.0;

  //f[0]   = (double *) callocvf ((xndof*nodes[0] + 1), sizeof (double));
    f[0] = new double [xndof*nodes[0] + 1]; for (int i = 0; i < xndof * nodes[0] + 1; i++)  f[0][i] = 0.0;

  //diag[0]= (double *) callocvf ((xndof*nodes[0] + 1), sizeof (double));
    diag[0] = new double [xndof*nodes[0] + 1]; for (int i = 0; i < xndof * nodes[0] + 1; i++)  diag[0][i] = 0.0;

  //u1[0]  = (double *) callocvf ((xndof*nodes[0] + 1), sizeof (double));
    u1[0] = new double [xndof*nodes[0] + 1]; for (int i = 0; i < xndof * nodes[0] + 1; i++)  u1[0][i] = 0.0;


//mdiag = (double ***) callocvf((maxlev+1), sizeof(double **));
////double **mdiag [maxlev+1];
//mdiag[0]= (double **) callocvf((maxlev+1)*(maxlev+2)/2, sizeof(double **));
//for (i=1; i<=maxlev ; i++) mdiag[i] = mdiag[i-1] + i;
//mdiag[0][0]= (double *) callocvf ((xndof*nodes[0] + 1), sizeof (double));
//double mdiag[maxlev+1][(maxlev+2)/2][];

 //du[0][0]  = (double *) callocvf ((xndof*nodes[0] + 1), sizeof (double));
 //res[0][0]  = (double *) callocvf ((xndof*nodes[0] + 1), sizeof (double));

  //r[0]   = (double *) callocvf ((xndof*nodes[0] + 1), sizeof (double));
  //p[0]   = (double *) callocvf ((xndof*nodes[0] + 1), sizeof (double));
  //z[0]   = (double *) callocvf ((xndof*nodes[0] + 1), sizeof (double));
  //fi[0]  = (double *) callocvf ((xndof*nodes[0] + 1), sizeof (double));
    fi[0]  = new double[xndof*nodes[0] + 1]; for (int i = 0; i < xndof * nodes[0] + 1; i++)  fi[0][i] = 0.0;

//vec[0]   = (double *) callocvf ((xndof*nodes[0] + 1), sizeof (double));

/*____________________________________________________________________________

							  INPUT EXTERNAL NODAL FORCES
____________________________________________________________________________*/


  nodvec (fpi, xndof, id, fp);

/*____________________________________________________________________________

		 ASSEMBLING OF STIFFNESS MATRIX AND COMPUTATION OF Gij,S
____________________________________________________________________________*/


 stiffblk (ecode, numel[0], sides[0], Bsides[0], xndof, nudim, id,
           elm_sides[0], side_nodes[0], side_bc, xnsn, xyz[0], side_mat, mat,
           prop, xnumpr, sg, maxa, lm, sv, Gij, Kij, diag[0], fp, fm);

/*____________________________________________________________________________

					 ADDING MIXED NODAL FORCES
____________________________________________________________________________*/

  for (i=1; i<=xndof*nodes[0] ; i++) fp[i] += fm[i];

/*____________________________________________________________________________

							INPUT EXTERNAL SIDE LOADS DESCRIPTION
____________________________________________________________________________*/


  sidevec (fpi, Bsides[0], xdf, xnbdof, side_loads);

/*____________________________________________________________________________

			  COMPUTATION OF NODAL SIDE LOADS IN GLOBAL X,Y COORDINATES
____________________________________________________________________________*/


 sideload (ecode, Bsides[0], xdf, xndof, xnsn, side_nodes[0], side_loads,
           bs_Fs[0], bs_Fn[0], xyz[0], f[0], idb, bn_cs, bn_sn, bs_cs[0],
           bs_sn[0], side_mat, prop, xnumpr, nudim);

/*____________________________________________________________________________

      			  ADDITION OF CONCENTRATED LOADS IN COARSE MESH
____________________________________________________________________________*/

 for (i=1; i<=xndof*nodes[0] ; i++) f[0][i] += fp[i];

/*____________________________________________________________________________

						COMPRESS VECTOR FORCES f[0] IN VECTOR dn
____________________________________________________________________________*/

 if ((ecode==4) || (ecode==5))
   compresss (dn, f[0], id, nodes[0]);
 else
   compressv2D (dn, f[0], id, nodes[0], xndof, idb, bn_cs, bn_sn);

/*____________________________________________________________________________

								  GAUSS FACTORIZATION
____________________________________________________________________________*/


 decomp (sg, maxa, neq);

/*____________________________________________________________________________

							 DISPLACEMENTS COMPUTATION
____________________________________________________________________________*/


 redbak (sg, dn, maxa, neq);

/*____________________________________________________________________________

							EXPAND DISPLACEMENTS dn IN VECTOR u1[0]
____________________________________________________________________________*/

 if ((ecode==4) || (ecode==5))
   expands (dn, u1[0], id, nodes[0]);
 else
   expandv2D (dn, u1[0], id, nodes[0], xndof, idb, bn_cs, bn_sn);
/*____________________________________________________________________________

								 COARSE MESH PARAMETERS
____________________________________________________________________________*/


 acnodes[0] = 0;
 acnumel[0] = 0;

 FRnumel[0] = 0;
 URnumel[0] = numel[0];
 Rnumel[0]  = numel[0];
 Unumel[0]  = 0;

 plnumel[0] = 0;

 FRsides[0] = 0;
 URsides[0] = sides[0];
 Rsides[0]  = sides[0];
 Usides[0]  = 0;

 BFRsides[0] = 0;
 BURsides[0] = Bsides[0];


 //parent_side[0] = (int *) callocvf ((2*Bsides[0]+1), sizeof (int));
   parent_side[0]=new int [2*Bsides[0]+1]; for (int i = 0; i < 2 * Bsides[0] + 1; i++)  parent_side[0][i] = 0;


 //parent_elm[0]  = (int *) callocvf ((3*numel[0]+1), sizeof (int));
   parent_elm[0]=new int [2*Bsides[0]+1]; for (int i = 0; i < 2 * Bsides[0] + 1; i++)  parent_elm[0][i] = 0;

 //vside[0] = (int *) callocvf ((sides[0]+1), sizeof (int));
   vside[0]= new int [sides[0]+1]; for (int i = 0; i < sides[0] + 1; i++)  vside[0][i] = 0;

 //eside[0] = (int *) callocvf ((sides[0]+1), sizeof (int));
   eside[0] = new int [sides[0]+1]; for (int i = 0; i < sides[0] + 1; i++)  eside[0][i] = 0;

 //bside[0] = (int *) callocvf ((sides[0]+1), sizeof (int));
   bside[0] = new int [sides[0]+1]; for (int i = 0; i < sides[0] + 1; i++)  bside[0][i] = 0;

 //tside[0] = (int *) callocvf ((sides[0]+1), sizeof (int));
   tside[0] = new int [sides[0]+1]; for (int i = 0; i < sides[0] + 1; i++)  tside[0][i] = 0;

 //ssside[0] = (int *) callocvf ((sides[0]+1), sizeof (int));
   ssside[0] = new int [sides[0]+1]; for (int i = 0; i < sides[0] + 1; i++)  ssside[0][i] = 0;

 //seside[0] = (int *) callocvf ((numel[0]+1), sizeof (int));
   seside[0] = new int [numel[0]+1]; for (int i = 0; i < numel[0] + 1; i++)  seside[0][i] = 0;

 //velm [0] = (int *) callocvf ((numel[0]+1), sizeof (int));
   velm [0] = new int [numel[0]+1]; for (int i = 0; i < numel[0] + 1; i++)  seside[0][i] = 0;

 //relm [0] = (int *) callocvf ((numel[0]+1), sizeof (int));
   relm [0] = new int [numel[0]+1]; for (int i = 0; i < numel[0] + 1; i++)  relm[0][i] = 0;


 for (i=1; i<=sides[0]; i++)  vside[0][i]=i;


 for (i=1; i<=Bsides[0] ; i++)
	   {
      bside[0][i] = i;
		    parent_side[0][2*i-1] = i;
		    parent_side[0][2*i] = 1;
	   }

 for (i=1; i<=numel[0] ; i++)
	   {
      velm[0][i] = i;
		    parent_elm[0][3*i-2] = i;
		    parent_elm[0][3*i-1] = 1;
		    parent_elm[0][3*i]   = 1;
	   }

/*____________________________________________________________________________

		     ALLOCATION MEMORY FOR REFINAMENT INDICATORS OF INITIAL LEVEL
____________________________________________________________________________*/
      //double *Sx[maxlev+1];
        Sx[0]=new double [nodes[0]+1]; for (int i = 0; i < nodes[0] + 1; i++)  Sx[0][i] = 0.0;

      //Sy[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
        Sy[0]   = new double [nodes[0]+1]; for (int i = 0; i < nodes[0] + 1; i++)  Sy[0][i] = 0.0;

      //Sxy[0]  = (double  *) callocvf ((nodes[0]+1), sizeof (double));
        Sxy[0]  = new double [nodes[0]+1]; for (int i = 0; i < nodes[0] + 1; i++)  Sxy[0][i] = 0.0;

      //nse     = (int    *) callocvf ((nodes[0]+1), sizeof (int));
        nse     = new int [nodes[0]+1]; for (int i = 0; i < nodes[0] + 1; i++)  nse[i] = 0;


      for (i=1; i<=numel[0] ; i++)
         {
           ++nse[elm_nodes[0][3*i-2]];
           ++nse[elm_nodes[0][3*i-1]];
           ++nse[elm_nodes[0][3*i]];

         }

// switch (ESTIMATOR)
//   {
//    case 1: /* nodal averaging */
//      //Sx[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double Sx[maxlev+1][nodes[0]+1];
//
//      //Sx[0]=new double [nodes[0]+1];
//
//      //Sy[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double Sy[maxlev+1][nodes[0]+1];
//
//      //Sxy[0]  = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double Sxy[maxlev+1][nodes[0]+1];
//
//      //nse     = (int    *) callocvf ((nodes[0]+1), sizeof (int));
//      int nse[maxlev+1][nodes[0]+1];
//
//
//      for (i=1; i<=numel[0] ; i++)
//         {
//           ++nse[elm_nodes[0][3*i-2]];
//           ++nse[elm_nodes[0][3*i-1]];
//           ++nse[elm_nodes[0][3*i]];
//         }
//    break;
//
//    case 2: /* midside averaging */
//      //Fn[0]    = (double  *) callocvf (( sides[0]+1), sizeof (double ));
//      double Fn[maxlev+1][sides[0]+1];
//
//      //Fs[0]    = (double  *) callocvf (( sides[0]+1), sizeof (double ));
//      double Fs[maxlev+1][sides[0]+1];
//
//      //bn_qn[0] = (double *) callocvf (( nodes[0]+1) , sizeof(double));
//      double bn_qn[maxlev+1][sides[0]+1];
//
//      //bn_qs[0] = (double *) callocvf (( nodes[0]+1) , sizeof(double));
//      double bn_qs[maxlev+1][sides[0]+1];
//
//      //nse     = (int    *) callocvf ((nodes[0]+1), sizeof (int));
//      int nse[nodes[0]+1];
//
//    break;
//
//    case 3: /* ZZ smoother procedure */
//    case 4: /* ZZ smoother procedure */
//      //Sx[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double Sx[maxlev+1][nodes[0]+1];
//
//      //Sy[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double Sy[maxlev+1][nodes[0]+1];
//
//      //Sxy[0]  = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double Sxy[maxlev+1][nodes[0]+1];
//
//      nse     = (int    *) callocvf ((nodes[0]+1), sizeof (int));
//
//      for (i=1; i<=numel[0] ; i++)
//         {
//           ++nse[elm_nodes[0][3*i-2]];
//           ++nse[elm_nodes[0][3*i-1]];
//           ++nse[elm_nodes[0][3*i]];
//         }
//
//      //a11[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double a11[maxlev+1][nodes[0]+1];
//
//      //a12[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double a12[maxlev+1][nodes[0]+1];
//
//      //a13[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double a13[maxlev+1][nodes[0]+1];
//
//      //a22[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double a22[maxlev+1][nodes[0]+1];
//
//      //a23[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double a23[maxlev+1][nodes[0]+1];
//
//      //a33[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double a33[maxlev+1][nodes[0]+1];
//
//      //b1x[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double b1x[maxlev+1][nodes[0]+1];
//
//      //b2x[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double b2x[maxlev+1][nodes[0]+1];
//
//      //b3x[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double b3x[maxlev+1][nodes[0]+1];
//
//      //b1y[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double b1y[maxlev+1][nodes[0]+1];
//
//      //b2y[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double b2y[maxlev+1][nodes[0]+1];
//
//      //b3y[0]   = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double b3y[maxlev+1][nodes[0]+1];
//
//      //b1xy[0]  = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double b1xy[maxlev+1][nodes[0]+1];
//
//      //b2xy[0]  = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double b2xy[maxlev+1][nodes[0]+1];
//
//      //b3xy[0]  = (double  *) callocvf ((nodes[0]+1), sizeof (double));
//      double b3xy[maxlev+1][nodes[0]+1];
//
//    break;
//
//    case 5: /* ZZ side averaging */
//      //Sx[0]   = (double  *) callocvf ((sides[0]+1), sizeof (double));
//      double Sx[maxlev+1][sides[0]+1];
//
//      //Sy[0]   = (double  *) callocvf ((sides[0]+1), sizeof (double));
//      double Sy[maxlev+1][sides[0]+1];
//
//      //Sxy[0]  = (double  *) callocvf ((sides[0]+1), sizeof (double));
//      double Sxy[maxlev+1][sides[0]+1];
//
//    break;
//
//   }
//


  //idnode[0] = (int *) callocvf ((nodes[0]+1), sizeof (int));
//  int idnode[maxlev+1][nodes[0]+1];

  //error[0]  = (double *) callocvf ((numel[0]+1), sizeof (double));
error[0] = new double[numel[0] + 1]; for (int i = 0; i < numel[0] + 1; i++)error[0][i] = 0.0;



 curlev = 0;
 newlev = 1;
 curmesh = 0;
 iter[0] = 0;
 v_cyc = 0;

/*____________________________________________________________________________

			        COMPUTE STRESS FIELD
____________________________________________________________________________*/


   refidfr (ecode, curlev, acnumel, numel, nodes, Rnumel, FRnumel,
            velm, xndof, nudim, xnsn, xyz, elm_sides, elm_nodes,
            elm_vlev, parent_elm,
            sides, Bsides, BFRsides, FRsides, PRsides, vside,
            eside, bside, tside, ssside, side_nodes, side_vlev,
            parent_side, side_bc, u1, mat, prop, xnumpr, side_mat,
            relm, error, Fn, Fs, Sx, Sy, Sxy, nse, idnode,
            bs_Fn, bs_Fs,
            a11, a12, a13, a22, a23, a33, b1x, b2x, b3x, b1y, b2y, b3y,
            b1xy, b2xy, b3xy,
            &refine, levref, ESTIMATOR, &TRnumel, &Uerr, &Unew, &Ufea, &errm,
            fi, bn_nlength, bn_slength, bn_qn, bn_qs, bs_cs, bs_sn);

/*____________________________________________________________________________

    PRINT DISPLACEMENTS AND STRESSES RECENTLY COMPUTED IN OUTPUT FILE
                      OF LAST MESH ANALIZED
____________________________________________________________________________*/


  prtplt (&forname, curlev, curmesh, ecode, nodes, acnodes, acnumel, numel,
      	  FRnumel, xndof, nudim, elm_sides, elm_nodes, elm_vlev, parent_elm,
          velm, xyz, mat, prop, xnumpr, u1, Fn, Fs, Sx, Sy, Sxy, error,
          errm, ESTIMATOR, SOLVER);

/*____________________________________________________________________________

								  END OF MAIN PROGRAM
____________________________________________________________________________*/

  printf("\n\n  PROGRAM TERMINATED");
  std::cout<<"\n\n  PROGRAM TERMINATED";


  return 0;

} /* end of main */


