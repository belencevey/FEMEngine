/**********************[          Prtplt.cpp            ]********************/
#include"Prtplt.h"


/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: PRTPLT.C                                                   #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #  void  prtplt()     : write geometry and stresses in output file.         #
 #  void  stress()     : write stresses in output file.                      #
 #                                                                           #
 #****************************************************************************/


/*******************************************[ USER LIBRARY : PRTPLT.C   ]*****
 #                                                                           #
 #  FUNCTION :  prtplt()                                                     #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Write geometry data and stresses in output file for plot.            #
 #      Only for triangles.                                                  #
 #                                                                           #
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
 #  TRIM2DS.C :                                                              #
 #    void  trim2ds() : computation of stresses for linear triangular        #
 #                      elements for plane states.                           #
 #                                                                           #
 #                                                                           #
 #***************************************************************************/


  void prtplt (std::string *forname, int curlev, int curmesh, int ecode, int *numnp,
               int *acnodes, int *acnumel, int *numel, int *frnumel, int xndof,
			   int nudim, int **elm_sides, int **elm_nodes, int **elm_vlev,
			   int **parent_elm, int **velm, double **xyz, int *mat, double *prop,
			   int xnumpr, double **uold, double **Fn, double **Fs, double **Sx,
			   double **Sy, double **Sxy, double **error, double errm,
			   int ESTIMATOR, int SOLVER)
  /*void prtplt (char *forname, int curlev, int curmesh, int ecode, int *numnp,
               int *acnodes, int *acnumel, int *numel, int *frnumel, int xndof,
			   int nudim, int **elm_sides, int **elm_nodes, int **elm_vlev,
			   int **parent_elm, int **velm, double **xyz, int *mat, double *prop,
			   int xnumpr, double **uold, double **Fn, double **Fs, double **Sx,
			   double **Sy, double **Sxy, double **error, double errm,
			   int ESTIMATOR, int SOLVER)*/

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// char *forname   ; /* complete output file root name                        */
// int  curlev     ; /* fine level in analysis.                               */
// int  curmesh    ; /* current mesh in analysis.                             */
// int  ecode      ; /* element code number.                                  */
// int *numnp      ; /* number of nodal points                                */
// int *acnodes    ; /* number of accumulated nodes in previous levels        */
// int *acnumel    ; /* number of accumulated elements in previous levels     */
// int *numel      ; /* number of elements in each level.                     */
// int *frnumel    ; /* number of refined elements in each level.             */
// int  xndof      ; /* maximum number of nodal degrees of freedom.           */
// int  nudim      ; /* number of spatial dimensions.                         */
// int **elm_sides ; /* sides incidences of each element [numel * xnes]       */
// int **elm_nodes ; /* nodes of each element    [numel[curlev]*xnen]         */
// int **elm_vlev  ; /* levels of element vertex nodes  [xnev*numel[curlev]]  */
// int **parent_elm; /* parent element and location [3*numel[curlev]]         */
// int **velm      ; /* new element ordering after refinament[numel[curlev]]  */
// double **xyz    ; /* nodal coordinates [numnp * nudim].                    */
// int    *mat     ; /* element material type vector [numel].                 */
// double *prop    ; /* material properties vector [numat * numpr].           */
// int    xnumpr   ; /* maximum number of material properties.                */
// double **uold   ; /* vector of nodal displacements [xndof*numnp]           */
// double  **Fn     ; /* midside normal traction       [sides[level]]          */
// double  **Fs     ; /* midside tangential traction   [sides[level]]          */
// double  **Sx     ; /* nodal stress Sx                 [nodes[level]]        */
// double  **Sy     ; /* nodal stress Sy                 [nodes[level]]        */
// double  **Sxy    ; /* nodal stress Sxy                [nodes[level]]        */
// double **error  ; /* element error                 [numel[level]]          */
// double errm     ; /* mean error in energy norm of fine mesh                */
// int ESTIMATOR; /* estimator type 1: nodal avg.  2: midside stress  3: ZZ   */
// int  SOLVER  ;  /* type of iterative solver                                */
//                 /*    SOLVER = 1   hierarchical PCG                        */
//                 /*    SOLVER = 2   diagonal PCG                            */
//                 /*    SOLVER = 3   hierarchical MULTIGRID                  */
//                 /*    SOLVER = 4   diagonal MULTIGRID                      */

/*____________________________________________________________________________

                        BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

                     LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

//FILE   *fpmvw  ;  /* pointer to meshview file                               */
std::ofstream fpmvw;
//FILE   *fGiDmesh; /* pointer to GiD mesh file                               */
std::ofstream fGiDmesh;
//FILE   *fGiDres;  /* pointer to GiD results file                            */
std::ofstream fGiDres;
int  elm1, elm2;  /* first and last elements of each level.                 */
int  level   ;    /* loop counter on levels                                 */
int     i, j ;    /* loop counters                                          */
int   i1     ;    /* auxiliar pointer for mesh arrays                       */
int  idnode  ;    /* node number identificator                              */
//char foname[200];  /* complete output file name                              */
std::string foname;
//char meshname[200];
std::string meshname;

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

void stress (std::ofstream &fpmvw      , std::ofstream &fGiD, int ecode     , int *acnodes,
             int elm1       , int elm2      , int xndof,
             int nudim      , int *elm_sides, int *elm_nodes,
             int *elm_vlev  ,
             int *parent_elm, int *velm     , double **xyz,
             int *mat       , double *prop  , int xnumpr,
             double **uold  , double  *Fn, double  *Fs,
             double **Sx     , double **Sy, double **Sxy,
             double *error  , double errm, int ESTIMATOR,
             int level);
/*void stress (FILE *fpmvw      , FILE *fGiD, int ecode     , int *acnodes,
             int elm1       , int elm2      , int xndof,
             int nudim      , int *elm_sides, int *elm_nodes,
             int *elm_vlev  ,
             int *parent_elm, int *velm     , double **xyz,
             int *mat       , double *prop  , int xnumpr,
             double **uold  , double  *Fn, double  *Fs,
             double **Sx     , double **Sy, double **Sxy,
             double *error  , double errm, int ESTIMATOR,
             int level);*/

void elmdesc (std::ofstream &fGiDmesh, int *acnodes,
             int elm1       , int elm2      ,
             int *elm_nodes,
             int *elm_vlev  ,
             int *parent_elm, int *velm     , int *mat);
/*void elmdesc (FILE *fGiDmesh, int *acnodes,
             int elm1       , int elm2      ,
             int *elm_nodes,
             int *elm_vlev  ,
             int *parent_elm, int *velm     , int *mat);*/

//extern FILE *fopenvf (char *finame, char *mode); ////////////////////////////ESTO NO VA MAS
/*extern FILE *fopenvf (char *finame, char *mode);*/

/*____________________________________________________________________________

                    OPENING OF MESHVIEW OUTPUT FILES
____________________________________________________________________________*/

			//strcpy (foname, forname);
			//strcat (foname, ".mvw");
			foname=(*forname)+".mvw";

   //fpmvw = fopenvf (foname, "w");
   fpmvw.open(foname.c_str());
/*____________________________________________________________________________

                   OPENING OF GiD OUTPUT FILES
____________________________________________________________________________*/

			//strcpy (foname, forname);
			//strcat (foname, ".post.msh");
			foname=(*forname)+".post.msh";

   //fGiDmesh = fopenvf (foname, "w");
   fGiDmesh.open(foname.c_str());

			//strcpy (foname, forname);
			//strcat (foname, ".post.res");
			foname=(*forname)+".post.res";

   //fGiDres = fopenvf (foname, "w");
   fGiDres.open(foname.c_str());
   fGiDres << std::setiosflags(std::ios::fixed) << std::setiosflags(std::ios::showpoint) << std::setprecision(7); //Lo formateo 
/*____________________________________________________________________________

                  PRINT GLOBAL DATA AND NUMBER OF FUNCTIONS
              AND NUMBER OF FUNCTIONS VALUES PER ELEMENT FUNCTION
____________________________________________________________________________*/

  		 switch (ecode)
		 {
		   case 1:
		   case 3:
		   case 7:
		   case 8:
		     /*fprintf (fpmvw,"%d %d %d %d %d", acnodes[curlev]+numnp[curlev],
                   acnumel[curlev] + numel[curlev], 4, xndof, 4);*/
		     fpmvw<<acnodes[curlev]+numnp[curlev]<<" "<<acnumel[curlev]+numel[curlev] << " " << 4 << " " << xndof << " " << 4<<"\n";

			 switch (ecode)
			 {
				 case 1:
					 //strcpy(meshname, "Malla T1 ");
					 meshname="\"Malla T1 \"";
					 break;
				 case 8:
					 //strcpy(meshname, "Malla T8 ");
					 meshname="\"Malla T8 \"";
					 break;
				 default:
					 //strcpy(meshname, "Malla ");
					 meshname="\"Malla \"";
			 }

		     // fprintf (fGiDmesh,"MESH \"%s\" dimension %d ElemType Triangle Nnode %d", forname, nudim, 3);



			 //fprintf(fGiDmesh, "MESH \"%s\" dimension 2 ElemType Triangle Nnode %d", meshname, 3);
			 fGiDmesh<<"MESH \""<<meshname<<"\" dimension 2 ElemType Triangle Nnode "<<3;
			 //fprintf(fGiDmesh, "\ncoordinates");
			 fGiDmesh<<"\ncoordinates";


			 //fprintf(fGiDres, "GiD Post Results File 1.0\n\n");
			 fGiDres<<"GiD Post Results File 1.0\n\n";
			 //fprintf (fGiDres,"GaussPoints \"Board\" ElemType Triangle \"%s\"", meshname);
			 fGiDres<<"GaussPoints \"Board\" ElemType Triangle \""<<meshname<<"\"";
             //fprintf (fGiDres,"\n  Number of Gauss Points: 3");
             fGiDres<<"\n  Number of Gauss Points: 3" << "\n";
             if (ESTIMATOR==1)
			 {
				 /*
				 fprintf(fGiDres, "\n  Natural Coordinates: Given");
				 fprintf(fGiDres, "\n      0.0 0.0");
				 fprintf(fGiDres, "\n      1.0 0.0");
				 fprintf(fGiDres, "\n      0.0 1.0");
				 */
				 fGiDres<< "\n  Natural Coordinates: Given";
				 fGiDres<< "\n      0.0 0.0";
				 fGiDres<< "\n      1.0 0.0";
				 fGiDres<< "\n      0.0 1.0";

			 }
			 else
			 {
				 //fprintf(fGiDres, "\n  Natural Coordinates: Internal");
				 fGiDres<< "\n  Natural Coordinates: Internal";
			 }

             /*
             fprintf (fGiDres,"\nEnd gausspoints\n");
             fprintf (fGiDres,"\nResult \"Desplazamientos\" \"Benchmark\" 1 Vector OnNodes");
             fprintf (fGiDres,"\nComponentNames \"Desp - X\", \"Desp - Y\", \"Desp - Z\"");
             fprintf (fGiDres,"\nValues");
             */
             fGiDres<<"\nEnd gausspoints\n";
             fGiDres<<"\nResult \"Desplazamientos\" \"Benchmark\" 1 Vector OnNodes"<<"\n";
             fGiDres<<"\nComponentNames \"Desp - X\", \"Desp - Y\", \"Desp - Z\"";
             fGiDres<<"\nValues\n";

		   break;

		   case 4:
		     /*fprintf (fpmvw,"%d %d %d %d %d", acnodes[curlev]+numnp[curlev],
                   acnumel[curlev] + numel[curlev], 4, xndof, 4);*/
		     fpmvw<<acnodes[curlev]+numnp[curlev] << " " << acnumel[curlev]+numel[curlev] << " " << 4<<xndof << " " << 4 << "\n";

		     /*fprintf (fGiDmesh,"MESH \"%s\" dimension %d ElemType Triangle Nnode %d",
                forname, nudim, 3);*/
		     fGiDmesh<<"MESH \ "<< forname<<"\ dimension "<<nudim<< "ElemType Triangle Nnode "<<3 << "\n";

             /*
		     fprintf (fGiDres,"GaussPoints \"Triangles\" ElemType Triangle");
             fprintf (fGiDres,"\nNumber of Gauss Points: 3");
             fprintf (fGiDres,"\nNatural Coordinates: Internal");
             fprintf (fGiDres,"\nEnd gausspoints");
             fprintf (fGiDres,"\nDesplazamientos  4 1.0 Estatic 2 1 1");
             fprintf (fGiDres,"\nX");
             fprintf (fGiDres,"\nY");
             fprintf (fGiDres,"\nZ");
             */
             fGiDres<<"GaussPoints \"Triangles\" ElemType Triangle";
             fGiDres<<"\nNumber of Gauss Points: 3";
             fGiDres<<"\nNatural Coordinates: Internal";
             fGiDres<<"\nEnd gausspoints";
             fGiDres<<"\nDesplazamientos  4 1.0 Estatic 2 1 1";
             fGiDres<<"\nX";
             fGiDres<<"\nY";
             fGiDres<<"\nZ";

		   break;

		   case 5:
		     /*fprintf (fpmvw,"%d %d %d %d %d", acnodes[curlev]+numnp[curlev],
                   acnumel[curlev] + numel[curlev], 4, xndof, 4);*/
		     fpmvw<<acnodes[curlev]+numnp[curlev] << " " << acnumel[curlev] + numel[curlev] << " " << 4 << " " << xndof << " " << 4 << "\n";

		     /*fprintf (fGiDmesh,"MESH \"%s\" dimension %d ElemType Triangle Nnode %d",
                forname, nudim, 3);*/
		     fGiDmesh<<"MESH \ "<<forname<<" dimension "<<nudim<<" ElemType Triangle Nnode "<<3 << "\n";
             /*
		     fprintf (fGiDres,"GaussPoints \"Triangles\" ElemType Triangle");
             fprintf (fGiDres,"\nNumber of Gauss Points: 3");
             fprintf (fGiDres,"\nNatural Coordinates: Internal");
             fprintf (fGiDres,"\nEnd gausspoints");
             fprintf (fGiDres,"\nDesplazamientos  4 1.0 Estatic 2 1 1");
             fprintf (fGiDres,"\nX");
             fprintf (fGiDres,"\nY");
             fprintf (fGiDres,"\nZ");
             */
		     fGiDres<<"GaussPoints \"Triangles\" ElemType Triangle";
             fGiDres<<"\nNumber of Gauss Points: 3";
             fGiDres<<"\nNatural Coordinates: Internal";
             fGiDres<<"\nEnd gausspoints";
             fGiDres<<"\nDesplazamientos  4 1.0 Estatic 2 1 1";
             fGiDres<<"\nX";
             fGiDres<<"\nY";
             fGiDres<<"\nZ";

		   break;

           default :
               /*printf ("\n FATAL ERROR:  ELEMENT CODE NOT RECOGNIZED");
						exit(1);*/
               std::cout<<"\n FATAL ERROR:  ELEMENT CODE NOT RECOGNIZED";
						exit(1);
		 }

		  /*
		  fprintf (fpmvw,"\n");
		  fprintf (fpmvw,"%d %d %d %d", 3, 3, 3, 1);
		  fprintf (fpmvw,"\n");
		  */
		 fpmvw<<"\n";
		 fpmvw<<3 << " " << 3 << " " << 3 << " " << 1;
		 fpmvw<<"\n";

/*____________________________________________________________________________

                  PRINT GiD MESH
____________________________________________________________________________*/


   idnode=0;

/*____________________________________________________________________________

                     LOOP ON LEVELS FOR NODAL VALUES
____________________________________________________________________________*/


 for (level=0; level<=curlev ; level++)
	 {

/*____________________________________________________________________________

                    LOOP ON NODAL POINTS OF EACH LEVEL
____________________________________________________________________________*/

		for ( i=1 ; i<=numnp[level] ; i++)
			{
/*____________________________________________________________________________

		       PRINT MESHVIEW NODAL COORDINATES
____________________________________________________________________________*/

			  i1 = nudim * (i -1);

			  //fprintf (fpmvw,"%12.5e   ", xyz[level][i1 +1]);
			  fpmvw<<xyz[level][i1 +1] << " ";
			  //fprintf (fpmvw,"%12.5e   ", xyz[level][i1 +2]);
			  fpmvw<<xyz[level][i1 +2] << " ";
			  if (nudim==3)
				 //fprintf (fpmvw,"%12.5e   ", xyz[level][i1 +3]);
				 fpmvw<<xyz[level][i1 +3] << " " << "\n";
			  else
				 //fprintf (fpmvw,"%12.5e   ", 0.0);
				 fpmvw<< 0.0 << "\n";

/*____________________________________________________________________________

		  	 PRINT MESHVIEW NODAL DISPLACEMENTS
____________________________________________________________________________*/

			  i1 = xndof*(i - 1);
			  for (j=1 ; j<=xndof ; j++)
				  //fprintf (fpmvw,"%+14.8e   ", uold[level][i1+j]);
				  fpmvw<<uold[level][i1+j]<<" ";
				  //fprintf (fpmvw,"\n");
				  fpmvw<<"\n";

/*____________________________________________________________________________

		       PRINT GiD NODAL COORDINATES
____________________________________________________________________________*/

			  idnode++;

              i1 = nudim * (i -1);

              //fprintf (fGiDmesh,"\n %d ", idnode);
              fGiDmesh << "\n" <<idnode<<"  ";
			  //fprintf (fGiDmesh,"%f   ", xyz[level][i1 +1]);
			  fGiDmesh<<xyz[level][i1 +1] << " ";
			  //fprintf (fGiDmesh,"%f   ", xyz[level][i1 +2]);
			  fGiDmesh<<xyz[level][i1 +2] << " ";
			  if (nudim==3)
				 //fprintf (fGiDmesh,"%f   ", xyz[level][i1 +3]);
				 fGiDmesh<<xyz[level][i1 +3] << " " << "\n";
			  else
				 //fprintf (fGiDmesh,"%f   ", 0.0);
				 fGiDmesh<< 0.0 << " " << "\n";

/*____________________________________________________________________________

		       PRINT GiD NODAL DISPLACEMENTS
____________________________________________________________________________*/

          //fprintf (fGiDres,"\n %d ", idnode);
          fGiDres<<idnode<<" ";
          i1 = xndof*(i - 1);
		  for (j=1 ; j<=xndof ; j++)
             //fprintf (fGiDres,"%f   ", uold[level][i1+j]);
             fGiDres<< uold[level][i1+j] << " ";
	      //if (xndof<3) fprintf (fGiDres,"%f   ", 0.0);
	      if (xndof<3) fGiDres<<0.0 << " \n";

/*__________________________________________________________________________

		  	 END OF LOOP ON NODES
 ___________________________________________________________________________*/

			}

/*__________________________________________________________________________

		  	END OF LOOP ON  LEVELS
 ___________________________________________________________________________*/

	 }

   //fprintf (fGiDmesh,"\nend coordinates");
   fGiDmesh<<"\nend coordinates";
   //fprintf(fGiDres, "\nEnd Values\n");
   fGiDres<<"\nEnd Values\n";
/*__________________________________________________________________________

						 LOOP ON LEVELS FOR ELEMENT DESCRIPTION
___________________________________________________________________________*/

 //fprintf (fGiDmesh,"\n\nelements");
 fGiDmesh<<"\n\nelements";

 for (level=0 ; level<=curlev ; level++)
	 {
		if (level<curlev)
		  {
  			 if (frnumel[level]<numel[level])
	  	     {
                 elm1 = frnumel[level]+1;
		         elm2 = numel[level];
                 elmdesc (fGiDmesh, acnodes, elm1, elm2, elm_nodes[level],
                          elm_vlev[level], parent_elm[level], velm[level], mat);

             }
		  }
		else
		  {
			   elm1 = 1;
			   elm2 = numel[level];
               elmdesc (fGiDmesh, acnodes, elm1, elm2, elm_nodes[level],
                        elm_vlev[level], parent_elm[level], velm[level], mat);

		  }

	 }

//fprintf (fGiDmesh,"\nend elements\n");
fGiDmesh<<"\nend elements\n";

/*__________________________________________________________________________

						 LOOP ON LEVELS FOR ELEMENT VALUES
___________________________________________________________________________*/

       //fprintf (fGiDres,"\nResult \"Tensiones\" \"Benchmark\" 1 Vector OnGaussPoints \"Board\"");
       fGiDres<<"\nResult \"Tensiones\" \"Benchmark\" 1 Vector OnGaussPoints \"Board\"";
       //fprintf (fGiDres,"\nComponentNames \"Sx\" \"Sy\" \"Sxy\"");
       fGiDres<<"\nComponentNames \"Sx\" \"Sy\" \"Sxy\"";
	   //fprintf(fGiDres, "\nValues");
	   fGiDres<<"\nValues";


 for (level=0 ; level<=curlev ; level++)
	 {

		if (level<curlev)
		  {
  			 if (frnumel[level]<numel[level])
		    		{
				      elm1 = frnumel[level]+1;
				      elm2 = numel[level];
				      stress (fpmvw, fGiDres, ecode, acnodes, elm1, elm2, xndof, nudim,
                      elm_sides[level], elm_nodes[level], elm_vlev[level],
                      parent_elm[level], velm[level], xyz, mat, prop, xnumpr,
                      uold, Fn[level], Fs[level], Sx, Sy, Sxy, error[level],
                      errm, ESTIMATOR, level);
        }
		  }
		else
		  {
			   elm1 = 1;
			   elm2 = numel[level];
			   stress (fpmvw, fGiDres, ecode, acnodes, elm1, elm2, xndof, nudim,
              elm_sides[level], elm_nodes[level], elm_vlev[level],
              parent_elm[level], velm[level], xyz, mat, prop, xnumpr,
              uold, Fn[level], Fs[level], Sx, Sy, Sxy, error[level],
              errm, ESTIMATOR, level);
		  }
	 }

     //fprintf(fGiDres, "\nEnd Values");
     fGiDres<< "\nEnd Values";
/*____________________________________________________________________________

								  PRINT ANALYSIS PARAMETERS
____________________________________________________________________________*/

 //fprintf (fpmvw,"\n\n TYPE OF ANALYSIS: ADAPATIVE REFINEMENT");
 fpmvw<<"\n\n TYPE OF ANALYSIS: ADAPATIVE REFINEMENT";
 //fprintf (fpmvw,"\n CURRENT LEVEL OF REFINEMENT : %d", curlev);
 fpmvw<<"\n CURRENT LEVEL OF REFINEMENT :"<<curlev << std::endl;
 switch (ESTIMATOR)
   {
     case 1:
      //fprintf (fpmvw,"\n ESTIMATOR TYPE : NODAL AVERAGING");
      fpmvw<<"\n ESTIMATOR TYPE : NODAL AVERAGING";
     break;

     case 2:
      //fprintf (fpmvw,"\n ESTIMATOR TYPE : MIDSIDE AVERAGING");
      fpmvw<<"\n ESTIMATOR TYPE : MIDSIDE AVERAGING";
     break;

     case 3:
      //fprintf (fpmvw,"\n ESTIMATOR TYPE : ZZ NODAL RECOVERY CENTROIDAL SAMPLING");
      fpmvw<<"\n ESTIMATOR TYPE : ZZ NODAL RECOVERY CENTROIDAL SAMPLING";
     break;

     case 4:
      //fprintf (fpmvw,"\n ESTIMATOR TYPE : ZZ NODAL RECOVERY SIDE SAMPLING");
      fpmvw<<"\n ESTIMATOR TYPE : ZZ NODAL RECOVERY SIDE SAMPLING";
     break;

     case 5:
      /*fprintf (fpmvw,
        "\n ESTIMATOR TYPE : SIMPLE SIDE AVERAGING (WITHOUT EXTERNAL FLUXES)");*/
     fpmvw<<"\n ESTIMATOR TYPE : SIMPLE SIDE AVERAGING (WITHOUT EXTERNAL FLUXES)";
     break;
   }


  switch (SOLVER)
     {
       case 1:
        //fprintf (fpmvw,"\n SOLVER TYPE : HIERARCHICAL PCG");
        fpmvw<<"\n SOLVER TYPE : HIERARCHICAL PCG";
       break;

       case 2:
        //fprintf (fpmvw,"\n SOLVER TYPE : DIAGONAL PCG");
        fpmvw<<"\n SOLVER TYPE : DIAGONAL PCG";
       break;

       case 3:
        //fprintf (fpmvw,"\n SOLVER TYPE : HIERARCHICAL MULTIGRID");
        fpmvw<<"\n SOLVER TYPE : HIERARCHICAL MULTIGRID";
       break;

       case 4:
        //fprintf (fpmvw,"\n SOLVER TYPE : DIAGONAL MULTIGRID");
        fpmvw<<"\n SOLVER TYPE : DIAGONAL MULTIGRID";
       break;
     }

/*____________________________________________________________________________

								  CLOSE OUPUT FILE
____________________________________________________________________________*/

  /*
  fclose (fpmvw);
  fclose (fGiDmesh);
  fclose (fGiDres);
  */
  fpmvw.close();
  fGiDmesh.close();
  fGiDres.close();

/*____________________________________________________________________________

				END OF FUNCTION
____________________________________________________________________________*/

 } /* end of prtplt() */







/*******************************************[ USER LIBRARY : PRTPLT.C   ]*****
 #                                                                           #
 #  FUNCTION :  stress()                                                     #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computation of local element stresses without interpolation.         #
 #      Only for triangles.                                                  #
 #                                                                           #
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
 #  TRIM2DS.C :                                                              #
 #    void  trim2ds() : computation of stresses for linear triangular        #
 #                      elements for plane states.                           #
 #                                                                           #
 #                                                                           #
 #***************************************************************************/
  void stress (std::ofstream &fpo, std::ofstream &fGiD, int ecode, int *acnodes, int elm1,
               int elm2, int xndof, int nudim, int *elm_sides, int *elm_nodes,
			   int *elm_vlev, int *parent_elm, int *velm, double **xyz,
			   int *mat,  double *prop, int xnumpr, double **uold,
			   double *Fn, double *Fs, double **Sx, double **Sy,
			   double **Sxy, double *error, double errm, int ESTIMATOR, int level)
/*
  void stress (FILE *fpo, FILE *fGiD, int ecode, int *acnodes, int elm1,
               int elm2, int xndof, int nudim, int *elm_sides, int *elm_nodes,
			   int *elm_vlev, int *parent_elm, int *velm, double **xyz,
			   int *mat,  double *prop, int xnumpr, double **uold,
			   double *Fn, double *Fs, double **Sx, double **Sy,
			   double **Sxy, double *error, double errm, int ESTIMATOR, int level)

*/
/*____________________________________________________________________________

								PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE   *fpo    ; /* pointer to output file                                 */
// FILE   *fGiD   ; /* pointer to GiD output file                             */
// int     ecode  ; /* element code number.                                   */
// int *acnodes   ; /* number of accumulated nodes in previous levels        */
// int  elm1, elm2; /* first and last elements of this level.                 */
// int  xndof     ; /* maximum number of nodal degrees of freedom.            */
// int  nudim     ; /* number of spatial dimensions.                          */
// int *elm_sides ; /* sides incidences of each element [numel * xnes]        */
// int *elm_nodes ; /* nodes of each element    [numel[curlev]*2*xnen]        */
// int *elm_vlev  ; /* levels of element vertex nodes  [xnev*numel[curlev]]   */
// int *parent_elm; /* parent element and location [3*numel[curlev]]          */
// int *velm      ; /* new element ordering after refinament[numel[curlev]]   */
// double **xyz   ; /* nodal coordinates [numnp * nudim].                     */
// int    *mat    ; /* element material type vector [numel].                  */
// double *prop   ; /* material properties vector [numat * numpr].            */
// int    xnumpr  ; /* maximum number of material properties.                 */
// double **uold  ; /* vector of nodal displacements [xndof*numnp]            */
// double  *Fn     ; /* midside normal traction       [sides[level]]           */
// double  *Fs     ; /* midside tangential traction   [sides[level]]           */
// double  **Sx    ; /* nodal stress Sx                    [nodes[level]]      */
// double  **Sy    ; /* nodal stress Sy                    [nodes[level]]      */
// double  **Sxy   ; /* nodal stress Sxy                   [nodes[level]]      */
// double *error  ; /* element error of this level    [numel[level]]          */
// double errm    ; /* mean error in energy norm of fine mesh                 */
// int ESTIMATOR  ; /* estimator type 1:nodal avg. 2:midside stress           */
// int  level     ; /* current level in analysis.                             */

/*____________________________________________________________________________

								BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

/////////////////////////////////////////////////////////////////////////////////////////////////VER


/*
extern  void trim2ds (FILE *fpo      , FILE *fGiD    , int *acnodes,
                      int elm1, int elm2       , int xndof     , int nudim,
                      int *elm_sides , int *elm_nodes, int *elm_vlev,
                      int *parent_elm, int *velm     , double **xyz,
                      int *mat       , double *prop  , int xnumpr,
                      double **uold  , double  *Fn    , double  *Fs,
                      double **Sx     , double **Sy    , double **Sxy,
                      double *error  , double errm   ,int ESTIMATOR,
                      int level);

extern  void trimm2ds (FILE *fpo      , FILE *fGiD    ,  int *acnodes,
                      int elm1,        int elm2       , int xndof     , int nudim,
                      int *elm_sides , int *elm_nodes, int *elm_vlev,
                      int *parent_elm, int *velm     , double **xyz,
                      int *mat       , double *prop  , int xnumpr,
                      double **uold  , double  *Fn    , double  *Fs,
                      double **Sx     , double **Sy    , double **Sxy,
                      double *error  , double errm   ,
                      int level);
*/

/*__________________________________________________________________________

			  LOOP ON SEVERAL BLOCKS OF ELEMENTS (NOT IMPLEMENTED)
 ___________________________________________________________________________*/

/*__________________________________________________________________________

					 ELEMENT STRESSES IN GLOBAL COORDINATES
___________________________________________________________________________*/

	 switch (ecode)
		{
		  case  1: /* linear triangle */
		  case  7: /* linear triangle by side */
		  case  8: /* linear triangle cubic Bezier */
			 trim2ds (fpo,fGiD, acnodes, elm1, elm2, xndof, nudim, elm_sides,
                   elm_nodes, elm_vlev, parent_elm, velm, xyz, mat,
                   prop, xnumpr, uold, Fn, Fs, Sx, Sy, Sxy, error,
                   errm, ESTIMATOR, level);
		  break;

		  case  3: /* linear mixed triangle */
			 trimm2ds (fpo, fGiD, acnodes, elm1, elm2, xndof, nudim, elm_sides,
             elm_nodes, elm_vlev, parent_elm, velm, xyz, mat,
             prop, xnumpr, uold, Fn, Fs, Sx, Sy, Sxy, error,
             errm, level);
		  break;

          case 4:

          break;

          case 5:

          break;

		  default:
			 printf ("\n FATAL ERROR:  ELEMENT CODE NOT RECOGNIZED");
			 exit(1);
		}


/*____________________________________________________________________________

				END OF FUNCTION
____________________________________________________________________________*/

 } /* end of stress() */



/*******************************************[ USER LIBRARY : PRTPLT.C  ]******
 #                                                                           #
 #  FUNCTION :  elmdesc ()                                                   #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Description of element incidences for Gid output.                    #
 #                                                                           #
 #***************************************************************************/

  void elmdesc (std::ofstream &fGiDmesh, int *acnodes, int elm1, int elm2, int *elm_nodes,
                int *elm_vlev, int *parent_elm, int *velm, int *mat)
  /*void elmdesc (FILE *fGiDmesh, int *acnodes, int elm1, int elm2, int *elm_nodes,
                int *elm_vlev, int *parent_elm, int *velm, int *mat)*/


/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE *fGiDmesh ; /* pointer to GiD mesh file                               */
// int *acnodes   ; /* number of accumulated nodes in previous levels         */
// int  elm1, elm2; /* first and last elements of this level.                 */
// int *elm_nodes ; /* nodes of each element    [numel[curlev]*xnen]          */
// int *elm_vlev  ; /* level of element vertexs [numel[curlev]*xnev]          */
// int *parent_elm; /* parent element and location [3*numel[curlev]]          */
// int *velm      ; /* new element ordering after refinament[numel[curlev]]   */
// int    *mat    ; /* element material type vector [numel].                  */

/*____________________________________________________________________________

			                BEGINNING OF FUNCTION
____________________________________________________________________________*/

{
/*____________________________________________________________________________

				          LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int ie        ;   /* element counter for loop on elements.                 */
 int n1, n2, n3;   /* element node numbers.                                 */
 int lv1, lv2, lv3;/* element node levels.                                  */
 int pelm      ;   /* parent element                                        */
 int i1        ;   /* auxiliar pointers for mesh arrays.                    */

/*____________________________________________________________________________

		               DECLARATION OF MATERIAL PROPERTIES
____________________________________________________________________________*/

 int  idmat   ;   /* identificator number of material                       */

/*__________________________________________________________________________

					LOOP ON UNREFINED ELEMENTS OF THIS LEVEL
___________________________________________________________________________*/

 for (ie=elm1; ie<=elm2 ; ie++)
	 {

/*________________________________________________________________________

							 SIDE AND NODAL INCIDENCES
 __________________________________________________________________________*/

/*		i1 = 3*(velm[ie]-1);
		L1 = elm_sides[++i1];
		L2 = elm_sides[++i1];
		L3 = elm_sides[++i1];
*/

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

/*__________________________________________________________________________

					MATERIAL PROPERTIES AND STIFFNESS COEFICIENTS
___________________________________________________________________________*/

		pelm = parent_elm[3*velm[ie]-2];

		idmat = mat[pelm];

/*____________________________________________________________________________

          						 PRINT ELEMENT INCIDENCES AND MATERIAL
____________________________________________________________________________*/


	 /*fprintf (fGiDmesh,"\n %d  %d  %d  %d  ", ie, acnodes[lv1]+n1, acnodes[lv2]+n2,
             acnodes[lv3]+n3);*/
    fGiDmesh<<ie<<" "<<acnodes[lv1]+n1<<" "<<acnodes[lv2]+n2<<" "<<acnodes[lv3]+n3 << std::endl;

/*________________________________________________________________________

		                  END OF LOOP ON ELEMENTS
 ________________________________________________________________________*/

	 }



/*________________________________________________________________________

				END OF FUNCTION
 ________________________________________________________________________*/

}/* end of elmdesc() */
