/**********************[           Genvec.cpp            ]********************/
#include"Genvec.h"


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

 void nodvec (std::ifstream &fpi, int xndof, int *id, double *f)
 //void nodvec (FILE *fpi, int xndof, int *id, double *f)

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE    *fpi ;   /* pointer input  file                                    */
// int     xndof;   /* maximum number of nodal degrees of freedom             */
// int     *id  ;   /* identification vector of numbered dof                  */
// double  *f   ;   /* vector of nodal components [neq].                      */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int    iread    ;   /* reading parameter (1 = read , 0 = end reading)      */
 int    node     ;   /* node number  (0 = generation command line)          */
 int    gcode    ;   /* generation code                                     */
 double vec      ;   /* vectorial component of vector f[]                   */
 int    i1       ;   /* auxiliar pointer for vector id[]                    */
 int    j1       ;   /* auxiliar pointer for vector f[]                     */
 int    j        ;   /* loop counter                                        */

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

void genvec (std::ifstream &fpi, int gcode, int *id, int xndof, double *f);
//void genvec (FILE *fpi, int gcode, int *id, int xndof, double *f);
/*____________________________________________________________________________

			 LOOP WHILE FOR READING
____________________________________________________________________________*/

 iread = 1;

 while (iread)
	 {
/*____________________________________________________________________________

                           READ NODE NUMBER
____________________________________________________________________________*/


		//fscanf(fpi,"%d",&node);
		fpi>>node;

/*____________________________________________________________________________

				READ NODAL VECTOR COMPONENTS
____________________________________________________________________________*/

	  if (node>0)
		   {
		  	  i1 = xndof * (node - 1);
			    for ( j=1 ; j<=xndof ; j++)
			      	{
				        //fscanf(fpi,"%lf",&vec);
				        fpi>>vec;
				        j1 = id[i1 + j];
				        if (j1 > 0) f[i1+j] = vec;
				      }
		   }

/*____________________________________________________________________________

			  VERIFICATION END OF READING
____________________________________________________________________________*/

	  else
		  {
  	  //fscanf(fpi,"%d",&gcode);
  	  fpi>>gcode;
	    if (!gcode)
		     iread = 0;

/*____________________________________________________________________________

		  GENERATION OF NODAL VECTORS COMPONENTS
____________________________________________________________________________*/

    else
    {
        if (gcode<1 || gcode>2)
        {
            /*
            printf ("\n");
            printf ("FATAL ERROR: ");
            printf ("NODAL VECTOR GENERATION CODE NOT RECOGNIZED");
            printf ("\n");
            */
            std::cout<<std::endl;
            std::cout<<"FATAL ERROR: ";
            std::cout<<"NODAL VECTOR GENERATION CODE NOT RECOGNIZED";
            std::cout<<std::endl;
            exit(1);
        }
        genvec (fpi, gcode, id, xndof, f);
    }
/*____________________________________________________________________________

		  END OF VERIFICATION OF END OF READING
____________________________________________________________________________*/

		}
/*____________________________________________________________________________

			END OF WHILE FOR READING
____________________________________________________________________________*/

	 }
/*____________________________________________________________________________

				  END OF FUNCTION
____________________________________________________________________________*/

}/* end of nodvec() */









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

 void genvec (std::ifstream &fpi, int gcode, int *id, int xndof, double *f)
 //void genvec (FILE *fpi, int gcode, int *id, int xndof, double *f)

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE    *fpi ;   /* pointer input  file                                    */
// int    gcode ;   /* generation code                                        */
// int     *id  ;   /* identification vector of dof status.                   */
//						/*      0 = free                                      */
//						/*      1 = restrained                                */
// int    xndof ;   /* maximum number of nodal degrees of freedom.            */
// double  *f   ;   /* vector of nodal components [neq].                      */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int    node     ;   /* node number  (0 = generation command line)          */
 int    n1       ;   /* first node for generation                           */
 int    n2       ;   /* last  node for generation                           */
 int    inc      ;   /* node number increment for generation                */
 int    ns1      ;   /* number of intervals for generation                  */
 double vec      ;   /* vectorial component of vector f[]                   */
 int    i1       ;   /* auxiliar pointer for vector id[]                    */
 int    j1       ;   /* auxiliar pointer for vector f[]                     */
 int    i, j     ;   /* loop counters                                       */

/*____________________________________________________________________________

			  READING GENERATION PARAMETERS
____________________________________________________________________________*/

		  /*
		  fscanf(fpi,"%d",&n1);
		  fscanf(fpi,"%d",&n2);
		  fscanf(fpi,"%d",&inc);
		  */
		  fpi>>n1;
		  fpi>>n2;
		  fpi>>inc;

/*____________________________________________________________________________

			  VERIFICATION OF SUBDIVISIONS
____________________________________________________________________________*/


            if (!inc)
            {
                /*
                printf ("\n");
                printf ("FATAL ERROR: ");
                printf ("NULL INCREMENT IN NODAL VECTOR ASSIGNAMENT");
                printf ("\n");
				*/
                std::cout<<std::endl;
                std::cout<<"FATAL ERROR: ";
                std::cout<<"NULL INCREMENT IN NODAL VECTOR ASSIGNAMENT";
                std::cout<<std::endl;
				exit(1);

            }

		  ns1 = n2 - n1;
		  if (ns1<0) inc =-abs(inc);
		  if (fmod(abs(ns1),abs(inc)))
			 {
			     /*
			     printf ("\n");
				printf ("FATAL ERROR: ");
				printf ("ILLEGAL INCREMENT IN NODAL VECTOR ASSIGNAMENT");
				printf ("\n");
				*/
			    std::cout<<std::endl;
				std::cout<<"FATAL ERROR: ";
				std::cout<<"ILLEGAL INCREMENT IN NODAL VECTOR ASSIGNAMENT";
			    std::cout<<std::endl;

				exit(1);
			 }

		  ns1 /= inc;

/*____________________________________________________________________________

		  LOOP OVER THE NUMBER OF DEGREES OF FREEDOM
____________________________________________________________________________*/

		  for (i=1 ; i<=xndof ; i++)
		{
		  //fscanf (fpi,"%lf", &vec);
            fpi>>vec;
/*____________________________________________________________________________

		  SWITCH TO GENERATION ROUTINE
____________________________________________________________________________*/

		  switch (gcode)
			{
			  case 1 :
				 node = n1;
				 for (j=0 ; j<=ns1 ; j++)
					 {
						i1 = xndof * (node - 1) + i;
						j1 = id[i1];
						if (j1 > 0) f[j1] = vec;
						node += inc;
					 }
				 break;

			  case 2 :
				 if (vec)
					{
					  node = n1;
					  for (j=0 ; j<=ns1 ; j++)
						  {
							 i1 = xndof * (node - 1) + i;
							 j1 = id[i1];
							 if (j1 > 0) f[j1] += vec;
							 node += inc;
						  }
					}
				 break;
			  default :;
			}
/*____________________________________________________________________________

			 END OF LOOP OVER THE NUMBER OF DEGREES OF FREEDOM
____________________________________________________________________________*/

		}
/*____________________________________________________________________________

				  END OF FUNCTION
____________________________________________________________________________*/

}/* end of genvec() */





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

void sidevec (std::ifstream &fpi, int bsides, int xdf, int xnbdof, double **side_loads)
 //void sidevec (FILE *fpi, int bsides, int xdf, int xnbdof, double **side_loads)

/*____________________________________________________________________________

									  PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE    *fpi       ; /* pointer input  file                                */
// int      bsides    ; /* number of boundary sides.                          */
// int      xdf       ; /* maximum polynomial degree of external vector field.*/
// int      xnbdof    ; /* maximum number of boundary degrees of freedom      */
// double **side_loads; /* boundary side loads [xndof][bsides * (xdf+1)]      */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int    iread    ;   /* reading parameter (1 = read , 0 = end reading)      */
 int    side     ;   /* side number  (0 = generation command line)          */
 int    dof      ;   /* vector field component                              */
 int    gcode    ;   /* generation code                                     */
 double vec      ;   /* vectorial component of vector f[]                   */
 int    i1       ;   /* auxiliar pointer for vector side_loads[]            */
 int    i        ;   /* loop counter                                        */

/*____________________________________________________________________________

			PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/

/*void genvec (FILE *fpi, int gcode, int *id, int xndof, double *f);*/

/*____________________________________________________________________________

								LOOP WHILE FOR READING
____________________________________________________________________________*/

 iread = 1;

 while (iread)
	 {
/*____________________________________________________________________________

								  READ SIDE NUMBER
____________________________________________________________________________*/


	  //fscanf(fpi,"%d",&side);
	  fpi>>side;
	  if (side>bsides || side<0)
		 {
			/*
			printf ("\n");
			printf ("FATAL ERROR: ");
			printf ("ILLEGAL SIDE NUMBER IN SIDE LOADS");
			printf ("\n");
			*/
			std::cout<<std::endl;
			std::cout<<"FATAL ERROR: ";
			std::cout<<"ILLEGAL SIDE NUMBER IN SIDE LOADS";
			std::cout<<std::endl;
			exit(1);
		 }


/*____________________________________________________________________________

						  READ VECTOR COMPONENT DIRECTION
____________________________________________________________________________*/

	  if (side>0)
		 {
			//fscanf (fpi,"%d",&dof);
			fpi>>dof;
			if (dof>xnbdof || dof<0)
			  {
				 /*
				 printf ("\n");
				 printf ("FATAL ERROR: ");
				 printf ("ILLEGAL DOF COMPONENT IN SIDE LOADS");
				 printf ("\n");
				 */
			     std::cout<<std::endl;
				 std::cout<<"FATAL ERROR: ";
				 std::cout<<"ILLEGAL DOF COMPONENT IN SIDE LOADS";
		         std::cout<<std::endl;
				 exit(1);
			  }

/*____________________________________________________________________________

						  READ VECTOR FIELD COEFFICIENTS
____________________________________________________________________________*/


			i1 = (side-1) * (xdf+1);
			for ( i=1 ; i<=xdf+1 ; i++)
				{
				  //fscanf(fpi,"%lf", &vec);
				  fpi>>vec;
				  side_loads[dof-1][i1 + i] = vec;
				}

		 }
/*____________________________________________________________________________

									GENERATION CODE
____________________________________________________________________________*/

	  else
		 {
			//fscanf(fpi,"%d",&gcode);
			fpi>>gcode;
			if (!gcode)
				iread = 0;

/*____________________________________________________________________________

						GENERATION OF SIDE VECTOR COMPONENTS
____________________________________________________________________________*/

			else
			  {
				 /*if (gcode<1 || gcode>2)*/
					{
					  /*
					  printf ("\n");
					  printf ("FATAL ERROR: ");
					  printf ("SIDE LOADS GENERATION CODE NOT IMPLEMENTED");
					  printf ("\n");
					  */
					  std::cout<<std::endl;
					  std::cout<<"FATAL ERROR: ";
					  std::cout<<"SIDE LOADS GENERATION CODE NOT IMPLEMENTED";
					  std::cout<<std::endl;
					  exit(1);
					}

				 /*genvec (fpi, gcode, id, xndof, f);*/

			  }
		 }
/*____________________________________________________________________________

								END OF WHILE FOR READING
____________________________________________________________________________*/

	 }
/*____________________________________________________________________________

									END OF FUNCTION
____________________________________________________________________________*/

}/* end of sidevec() */







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
                int *side_mat, double *prop, int numpr, int nudim)

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
/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int    side     ;   /* loop counter on boundary sides.                     */
 int    n1, n2, n3;  /* node numbers of side nodes                          */
 double L        ;   /* length of side n1-n2                                */
 double LX, LY   ;   /* length projections on global coordinates            */
 int    i1       ;   /* auxiliar pointers for vectors id[], f[]             */
 int    k1       ;   /* auxiliar pointer for vector side_nodes[]            */
 double p0, p1   ,   /* nodal values for polinomial description of          */
	    	  p2, p3   ;   /* the load variation along the side                   */
 double f0s, f0n ;   /* nodal resultants on n1 in side or nodal directions  */
 double f1s, f1n ;   /* nodal resultants on n2 in side or nodal directions  */
 double f0x, f0y ;   /* nodal resultants on n1 in global x,y system.        */
 double f1x, f1y ;   /* nodal resultants on n2 in global x,y system.        */
 double f2x, f2y ;   /* nodal resultants on n3 in global x,y system.        */
 double cosi   ;   /* cosine of side angle of rotation.                     */
 double sini   ;   /* sine of side angle of rotation.                       */
 double cos1   ;   /* cosine of angle of rotation of first side node.       */
 double sin1   ;   /* sine of angle of rotation of first side node.         */
 double cos2   ;   /* cosine of angle of rotation of second side node.      */
 double sin2   ;   /* sine of angle of rotation of second side node.        */
// double cos3   ;   /* cosine of angle of rotation of opposite side node.    */
// double sin3   ;   /* sine of angle of rotation of opposite side node.      */
 double Sn, Tsn;   /* midside tractions on boundary sides                   */
 double  E    ;     /* elasticity modulus                                   */
 double  nu   ;     /* Poisson's ratio                                      */
 double xij,yij,xjk,yjk,xki,yki,Li,Lj,Lk,ci,si,cj,sj,ck,sk;
 double CijS, CikS, CiiS, CijT, CikT, CiiT;

/*____________________________________________________________________________

								LOOP ON BOUNDARY SIDES
____________________________________________________________________________*/


	 for (side=1 ; side<=bsides ; side++)
	    {

/*____________________________________________________________________________

					COMPUTATION OF NODAL VECTOR COMPONENTS
____________________________________________________________________________*/

			i1 = xnsn * (side - 1);
			switch (ecode)
            {
			  case 1: /* linear triangle */
			  case 7: /* linear triangle by side */
			  case 8: /* linear triangle cubic Bezier */
                cosi = bs_cs[side];
                sini = bs_sn[side];

                n1 = side_nodes[i1+1];
                cos1 = bn_cs[idb[n1]];
                sin1 = bn_sn[idb[n1]];

                n2 = side_nodes[i1+2];
                cos2 = bn_cs[idb[n2]];
                sin2 = bn_sn[idb[n2]];

                LX = xyz[nudim*(n1-1)+1] - xyz[nudim*(n2-1)+1];
				LY = xyz[nudim*(n1-1)+2] - xyz[nudim*(n2-1)+2];
				L = sqrt (LX*LX + LY*LY);

				switch (xdf)
				{
				  case 1: /* linear load variation */
                  k1  = 2 * (side - 1);
                  p0  = side_loads[0][k1+1];
                  p1  = side_loads[0][k1+2];
                  Tsn = (p0 + p1)/2;
                  f0s = L * (2*p0 + p1)/6;
                  f1s = L * (2*p1 + p0)/6;
                  p0  = side_loads[1][k1+1];
                  p1  = side_loads[1][k1+2];
                  Sn  = (p0 + p1)/2;
                  f0n = L * (2*p0 + p1)/6;
                  f1n = L * (2*p1 + p0)/6;
			      break;

			      case 2: /* quadratic load variation */
                  k1  = 3 * (side - 1);
                  p0  = side_loads[0][k1+1];
                  p1  = side_loads[0][k1+2];
                  p2  = side_loads[0][k1+3];
                  Tsn = (p0+4*p1+p2)/6;
                  f0s = L * (2*p1 + p0)/6;
                  f1s = L * (2*p1 + p2)/6;
                  p0  = side_loads[1][k1+1];
                  p1  = side_loads[1][k1+2];
                  p2  = side_loads[1][k1+3];
                  Sn  = (p0+4*p1+p2)/6;
                  f0n = L * (2*p1 + p0)/6;
                  f1n = L * (2*p1 + p2)/6;
			      break;

			      case 3: /* cubic load variation */
                  k1  = 4 * (side - 1);
                  p0  = side_loads[0][k1+1];
                  p1  = side_loads[0][k1+2];
                  p2  = side_loads[0][k1+3];
                  p3  = side_loads[0][k1+4];
                  Tsn = (p2+p3)/2; /* corregir */
                  f0s = L * (13*p0 + 36*p1 +  9*p2 +  2*p3)/120;
                  f1s = L * ( 2*p0 +  9*p1 + 36*p2 + 13*p3)/120;
                  p0  = side_loads[1][k1+1];
                  p1  = side_loads[1][k1+2];
                  p2  = side_loads[1][k1+3];
                  p3  = side_loads[1][k1+4];
                  Sn  = (p2+p3)/2; /* corregir; */
                  f0n = L * (13*p0 + 36*p1 +  9*p2 +  2*p3)/120;
                  f1n = L * ( 2*p0 +  9*p1 + 36*p2 + 13*p3)/120;
			      break;
                }

                f0x = f0s*cosi-f0n*sini;
                f0y = f0s*sini+f0n*cosi;
                f0s =  f0x*cos1+f0y*sin1;
                f0n = -f0x*sin1+f0y*cos1;
                i1 = xndof*(n1-1);
                f[++i1] += f0s;
                f[++i1] += f0n;

                f1x =  f1s*cosi-f1n*sini;
                f1y =  f1s*sini+f1n*cosi;
                f1s =  f1x*cos2+f1y*sin2;
                f1n = -f1x*sin2+f1y*cos2;
                i1 = xndof*(n2-1);
                f[++i1] += f1s;
                f[++i1] += f1n;

                bs_Fs[side]=Tsn;
                bs_Fn[side]=Sn;
				break;

		      case 2: /* quadratic triangle side - not implemented */
			  break;

			  case 3: /* linear mixed triangle side */
                cosi = bs_cs[side];
                sini = bs_sn[side];

                n1 = side_nodes[i1+1];
                cos1 = bn_cs[idb[n1]];
                sin1 = bn_sn[idb[n1]];

                n2 = side_nodes[i1+2];
                cos2 = bn_cs[idb[n2]];
                sin2 = bn_sn[idb[n2]];

				LX = xyz[nudim*(n1 - 1) + 1] - xyz[nudim*(n2 - 1) + 1];
				LY = xyz[nudim*(n1 - 1) + 2] - xyz[nudim*(n2 - 1) + 2];
				L = sqrt (LX*LX + LY*LY);

				switch (xdf)
                {
	              case 1: /* linear load variation */
                  k1  = 2 * (side - 1);
                  p0  = side_loads[0][k1+1];
                  p1  = side_loads[0][k1+2];
                  Tsn = (p0 + p1)/2;
                  f0s = L * (2*p0 + p1)/6;
                  f1s = L * (2*p1 + p0)/6;
                  p0  = side_loads[1][k1+1];
                  p1  = side_loads[1][k1+2];
                  Sn  = (p0 + p1)/2;
                  f0n = L * (2*p0 + p1)/6;
                  f1n = L * (2*p1 + p0)/6;
				  break;

				  case 2: /* quadratic load variation */
                  k1  = 3 * (side - 1);
                  p0  = side_loads[0][k1+1];
                  p1  = side_loads[0][k1+2];
                  p2  = side_loads[0][k1+3];
                  Tsn = (p0+4*p1+p2)/6;
                  f0s = L * (2*p1 + p0)/6;
                  f1s = L * (2*p1 + p2)/6;
                  p0  = side_loads[1][k1+1];
                  p1  = side_loads[1][k1+2];
                  p2  = side_loads[1][k1+3];
                  Sn  = (p0+4*p1+p2)/6;
                  f0n = L * (2*p1 + p0)/6;
                  f1n = L * (2*p1 + p2)/6;
				  break;

				  case 3: /* cubic load variation */
                  k1  = 4 * (side - 1);
                  p0  = side_loads[0][k1+1];
                  p1  = side_loads[0][k1+2];
                  p2  = side_loads[0][k1+3];
                  p3  = side_loads[0][k1+4];
                  f0s = L * (13*p0 + 36*p1 +  9*p2 +  2*p3)/120;
                  f1s = L * ( 2*p0 +  9*p1 + 36*p2 + 13*p3)/120;
                  Tsn = (p2+p3)/2; /*corregir*/
                  p0  = side_loads[1][k1+1];
                  p1  = side_loads[1][k1+2];
                  p2  = side_loads[1][k1+3];
                  p3  = side_loads[1][k1+4];
                  Sn  = (p2+p3)/2; /*corregir;*/
                  f0n = L * (13*p0 + 36*p1 +  9*p2 +  2*p3)/120;
                  f1n = L * ( 2*p0 +  9*p1 + 36*p2 + 13*p3)/120;
				  break;

				}

         f0x = f0s*cosi-f0n*sini;
         f0y = f0s*sini+f0n *cosi;
         f0s =  f0x*cos1+f0y*sin1;
         f0n = -f0x*sin1+f0y*cos1;
         i1 = xndof*(n1-1);
         f[++i1] += f0s;
         f[++i1] += f0n;

         f1x =  f1s*cosi-f1n*sini;
         f1y =  f1s*sini+f1n*cosi;
         f1s =  f1x*cos2+f1y*sin2;
         f1n = -f1x*sin2+f1y*cos2;
         i1 = xndof*(n2-1);
         f[++i1] += f1s;
         f[++i1] += f1n;

         /* efecto de tensiones impuestas en desplazamientos
         // si el lado esta descargado o esta con desplazamiento impuesto
         // no afecta al resultado pues las cargas son nulas.*/

	     i1 = numpr * (side_mat[2*side-1] - 1);

		 E  = prop[++i1];
		 nu = prop[++i1];

         i1 = xnsn * (side - 1);
         n3 = side_nodes[i1+3];

         xjk=xyz[nudim*(n2-1)+1] - xyz[nudim*(n1-1)+1]; /*/x21*/
         yjk=xyz[nudim*(n2-1)+2] - xyz[nudim*(n1-1)+2]; /*/y21*/
         xki=xyz[nudim*(n1-1)+1] - xyz[nudim*(n3-1)+1]; /*/x13*/
         yki=xyz[nudim*(n1-1)+2] - xyz[nudim*(n3-1)+2]; /*/y13*/
         xij=xyz[nudim*(n3-1)+1] - xyz[nudim*(n2-1)+1]; /*/x32*/
         yij=xyz[nudim*(n3-1)+2] - xyz[nudim*(n2-1)+2]; /*/y32*/

         Li = sqrt (xjk*xjk + yjk*yjk);
         Lj = sqrt (xki*xki + yki*yki);
         Lk = sqrt (xij*xij + yij*yij);

         si=yjk/Li;
         ci=xjk/Li;
         sj=yki/Lj;
         cj=xki/Lj;
         sk=yij/Lk;
         ck=xij/Lk;

         /*/  Fx debida a Sni, Tsni*/
         CijS=Lj/6/(1-nu)*(si*(si*sj+ci*cj)-
                          nu*(sj*(si*si-ci*ci)+2*cj*si*ci));
         CikS=Lk/6/(1-nu)*(si*(si*sk+ci*ck)-
                          nu*(sk*(si*si-ci*ci)+2*ck*si*ci));
         CiiS=Li/6*si;

         CijT=Lj/6*(cj*(si*si-ci*ci)-2*sj*si*ci);
         CikT=Lk/6*(ck*(si*si-ci*ci)-2*sk*si*ci);
         CiiT=Li/6*(-ci);

         /*/  Fx en nodos n1,n2,n3 (k,j,i)*/
         f0x=-CikS*Sn-CikT*Tsn;
         f1x=-CijS*Sn-CijT*Tsn;
         f2x=-CiiS*Sn-CiiT*Tsn;

         /*/  Fy debida a Sni, Tsni*/
         CijS=Lj/6/(1-nu)*(-ci*(si*sj+ci*cj)+
                          nu*(cj*(ci*ci-si*si)+2*sj*si*ci));
         CikS=Lk/6/(1-nu)*(-ci*(si*sk+ci*ck)+
                          nu*(ck*(ci*ci-si*si)+2*sk*si*ci));
         CiiS=Li/6*(-ci);

         CijT=Lj/6*(sj*(ci*ci-si*si)-2*cj*si*ci);
         CikT=Lk/6*(sk*(ci*ci-si*si)-2*ck*si*ci);
         CiiT=Li/6*(-si);

         /*/  Fy en nodos n1,n2,n3 (k,j,i)*/
         f0y=-CikS*Sn-CikT*Tsn;
         f1y=-CijS*Sn-CijT*Tsn;
         f2y=-CiiS*Sn-CiiT*Tsn;

         /*/ fuerzas en nodo n1*/
         f0s =  f0x*cos1+f0y*sin1;
         f0n = -f0x*sin1+f0y*cos1;
         i1 = xndof*(n1-1);
         f[++i1] += f0s;
         f[++i1] += f0n;

          /*/ fuerzas en nodo n2*/
         f1s =  f1x*cos2+f1y*sin2;
         f1n = -f1x*sin2+f1y*cos2;
         i1 = xndof*(n2-1);
         f[++i1] += f1s;
         f[++i1] += f1n;

          /*/ fuerzas en nodo n3	CORREGIR angulo en nodos se corrige despues*/
         i1 = xndof*(n3-1);
         f[++i1] += f2x;
         f[++i1] += f2y;

         bs_Fs[side]=Tsn;
         bs_Fn[side]=Sn;

		 break;

/*____________________________________________________________________________

							  END OF SWITCH ON ELEMENT CODES
____________________________________________________________________________*/

         }

/*____________________________________________________________________________

							  END OF LOOP ON BOUNDARY SIDES
____________________________________________________________________________*/

		 }


/*____________________________________________________________________________

									END OF FUNCTION
____________________________________________________________________________*/

}/* end of sideload() */








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
                  int **parent_side, int curlev, double **xyz, double **f, double **bs_cs, double **bs_sn)

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


/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 int    level    ;   /* loop counter on levels.                             */
 int    side     ;   /* loop counter on boundary sides.                     */
 int    side1    ;   /* intial side of loop.                                */
 int    side2    ;   /* final side of loop.                                 */
 int    pside    ;   /* parent side of current side.                        */
 int    n1, n2   ;   /* node numbers of side nodes                          */
 int    lv1, lv2 ;   /* level of vertex nodes of parent side                */
 double L        ;   /* length of side n1-n2                                */
 double x1, y1   ;   /* coordinates of first side node of parent side.      */
 double x2, y2   ;   /* coordinates of second side node of parent side.     */
 double LX, LY   ;   /* length projections on global coordinates            */
 int    i1       ;   /* auxiliar pointers for vectors u1[]                  */
 int    k1       ;   /* auxiliar pointer for vector side_nodes[]            */
 double p0, p1   ,   /* nodal values for polinomial description of          */
		  p2, p3   ;   /* the load variation along the side                   */
 double pb0, pb1 ,   /* nodal values for polinomial description of          */
		  pb2, pb3 ;   /* the load variation along the interior segment       */
 double f1s, f2s ;   /* nodal resultants of the side tangential loads       */
 double f1n, f2n ;   /* nodal resultants of the side normal loads           */
 double f1x, f2x ;   /* nodal resultants of loads in x direction.           */
 double f1y, f2y ;   /* nodal resultants of loads in y direction.           */
 double sl       ;   /* maximum number of side subdivions                   */
 int    m        ;   /* position in parent side                             */
 double e1, e2   ;   /* coordinates of segment extremes                     */
 double eb0, eb1 ,   /* coordinates of lagrangian points in segment         */
        eb2, eb3 ;
 double d1, d2, d3;  /* auxiliar coefficients                               */
 double cosi   ;   /* cosine of side angle of rotation.                     */
 double sini   ;   /* sine of side angle of rotation.                       */

/*____________________________________________________________________________

									CLEAR FORCE VECTOR
____________________________________________________________________________*/

  for (level=0; level<=curlev ; level++)
       memset (f[level], 0, (xndof*nodes[level]+1)*sizeof(double));

/*____________________________________________________________________________

									LOOP ON PREVIOUS LEVELS
____________________________________________________________________________*/

  for (level=0 ; level<=curlev ; level++)
   {

    sl = (double)(1<<level);

/*____________________________________________________________________________

					LOOP ON UNREFINED BOUNDARY SIDES FOR COARSE MESHES
               LOOP ON ALL BOUNDARY SIDES FOR FINE MESH
____________________________________________________________________________*/

    if (level<curlev)
      {
        side1 = BFRsides[level];
        side2 = Bsides[level];
      }
    else
      {
        side1 = 0;
        side2 = Bsides[level];
      }

	 for (side=side1+1 ; side<=side2 ; side++)
		 {

/*____________________________________________________________________________

					COMPUTATION OF NODAL VECTOR COMPONENTS
____________________________________________________________________________*/

 		pside = parent_side[level][2*side-1];
			m     = parent_side[level][2*side  ];

   cosi = bs_cs[level][side];
   sini = bs_sn[level][side];

			i1 = xnsn * (pside - 1);
			switch (ecode)
					{
					  case 1: /* linear triangle side */

						 n1 = side_nodes[0][i1+1];
						 n2 = side_nodes[0][i1+2];

       i1 = xndof*(n1-1);
       x1 = xyz[0][++i1];
       y1 = xyz[0][++i1];

       i1 = xndof*(n2-1);
       x2 = xyz[0][++i1];
       y2 = xyz[0][++i1];

						 LX = x2 - x1;
						 LY = y2 - y1;
						 L = sqrt (LX*LX + LY*LY);
       L /= sl;

       i1 = xnsn*(bside[level][side]-1);
						 n1 = side_nodes[level][i1+1];
						 n2 = side_nodes[level][i1+2];
						 lv1 = side_vlev[level][i1+1];
						 lv2 = side_vlev[level][i1+2];

       switch (xdf)
             {
               case 1: /* linear load variation */
            	  k1 = 2 * (pside - 1);
            	  e1 = (m-1)/sl;
               e2 = m/sl;

            	  p0 = side_loads[0][k1+1];
            	  p1 = side_loads[0][k1+2];
               pb0 = e1*(p1-p0)+p0;
               pb1 = e2*p1+(1-e2)*p0;
               f1s = L * (2*pb0 + pb1)/6;
               f2s = L * (2*pb1 + pb0)/6;

            	  p0 = side_loads[1][k1+1];
            	  p1 = side_loads[1][k1+2];
               pb0 = e1*(p1-p0)+p0;
               pb1 = e2*p1+(1-e2)*p0;
               f1n = L * (2*pb0 + pb1)/6;
               f2n = L * (2*pb1 + pb0)/6;

               break;

            	case 2: /* quadratic load variation */
               k1 = 3 * (pside - 1);
               e1 = (m-1)/sl;
               e2 = m/sl;
               eb0 = e1;
               eb1 = (e1+e2)/2.;
               eb2 = e2;

            	  p0 = side_loads[0][k1+1];
            	  p1 = side_loads[0][k1+2];
            	  p2 = side_loads[0][k1+3];
               d1 = p0-2*p1+p2;
               d2 = 3*p0-4*p1+p2;
               pb0 = 2*eb0*eb0*d1-eb0*d2+p0;
               pb1 = 2*eb1*eb1*d1-eb1*d2+p0;
               pb2 = 2*eb2*eb2*d1-eb2*d2+p0;
               f1s = L * (pb0 + 2*pb1)/6;
               f2s = L * (2*pb1 + pb2)/6;

            	  p0 = side_loads[1][k1+1];
            	  p1 = side_loads[1][k1+2];
            	  p2 = side_loads[1][k1+3];
               d1 = p0-2*p1+p2;
               d2 = 3*p0-4*p1+p2;
               pb0 = 2*eb0*eb0*d1-eb0*d2+p0;
               pb1 = 2*eb1*eb1*d1-eb1*d2+p0;
               pb2 = 2*eb2*eb2*d1-eb2*d2+p0;
               f1n = L * (pb0 + 2*pb1)/6;
               f2n = L * (2*pb1 + pb2)/6;

               break;

            	case 3: /* cubic load variation */
               k1 = 4 * (pside - 1);
               e1 = (m-1)/sl;
               e2 = m/sl;
               eb0 = e1;
               eb1 = e1+(e2-e1)/3.;
               eb2 = e1+2*(e2-e1)/3.;
               eb3 = e2;

            	  p0 = side_loads[0][k1+1];
            	  p1 = side_loads[0][k1+2];
            	  p2 = side_loads[0][k1+3];
               p3 = side_loads[0][k1+4];
               d1 = p3-3*p2+3*p1-p0;
               d2 = 2*p0-5*p1+4*p2-p3;
               d3 = 2*p3-9*p2+18*p1-11*p0;
               pb0 = 9/2*eb0*eb0*(eb0*d1+d2)+eb0/2*d3+p0;
               pb1 = 9/2*eb1*eb1*(eb1*d1+d2)+eb1/2*d3+p0;
               pb2 = 9/2*eb2*eb2*(eb2*d1+d2)+eb2/2*d3+p0;
               pb3 = 9/2*eb3*eb3*(eb3*d1+d2)+eb3/2*d3+p0;
               f1s = L * (13*pb0+36*pb1+9*pb2+2*pb3)/120;
               f2s = L * (2*pb0+9*pb1+36*pb2+13*pb3)/120;

            	  p0 = side_loads[1][k1+1];
            	  p1 = side_loads[1][k1+2];
            	  p2 = side_loads[1][k1+3];
               p3 = side_loads[1][k1+4];
               d1 = p3-3*p2+3*p1-p0;
               d2 = 2*p0-5*p1+4*p2-p3;
               d3 = 2*p3-9*p2+18*p1-11*p0;
               pb0 = 9/2*eb0*eb0*(eb0*d1+d2)+eb0/2*d3+p0;
               pb1 = 9/2*eb1*eb1*(eb1*d1+d2)+eb1/2*d3+p0;
               pb2 = 9/2*eb2*eb2*(eb2*d1+d2)+eb2/2*d3+p0;
               pb3 = 9/2*eb3*eb3*(eb3*d1+d2)+eb3/2*d3+p0;
               f1n = L * (13*pb0+36*pb1+9*pb2+2*pb3)/120;
               f2n = L * (2*pb0+9*pb1+36*pb2+13*pb3)/120;

               break;

               default:
                 printf ("\n");
				             printf ("FATAL ERROR: ");
				             printf ("TOO HIGHER DEGREE IN BOUNDARY LOADS");
				             printf ("\n");
				             exit(1);
             }

       f1x = f1s*cosi-f1n*sini;
       f1y = f1s*sini+f1n*cosi;
       i1 = xndof*(n1-1);
       f[lv1][++i1] += f1x;
       f[lv1][++i1] += f1y;

       f2x = f2s*cosi-f2n*sini;
       f2y = f2s*sini+f2n*cosi;
       i1 = xndof*(n2-1);
       f[lv2][++i1] += f2x;
       f[lv2][++i1] += f2y;
					  break;

					  case 2: /* quadratic triangle side - not implemented */
					  break;

					  case 3: /* linear mixed triangle side */

						 n1 = side_nodes[0][i1+1];
						 n2 = side_nodes[0][i1+2];

       i1 = xndof*(n1-1);
       x1 = xyz[0][++i1];
       y1 = xyz[0][++i1];

       i1 = xndof*(n2-1);
       x2 = xyz[0][++i1];
       y2 = xyz[0][++i1];

						 LX = x2 - x1;
						 LY = y2 - y1;
						 L = sqrt (LX*LX + LY*LY);
       L /= sl;

       i1 = xnsn*(bside[level][side]-1);
						 n1 = side_nodes[level][i1+1];
						 n2 = side_nodes[level][i1+2];
						 lv1 = side_vlev[level][i1+1];
						 lv2 = side_vlev[level][i1+2];

       switch (xdf)
             {
               case 1: /* linear load variation */
            	  k1 = 2 * (pside - 1);
            	  e1 = (m-1)/sl;
               e2 = m/sl;

            	  p0 = side_loads[0][k1+1];
            	  p1 = side_loads[0][k1+2];
               pb0 = e1*(p1-p0)+p0;
               pb1 = e2*p1+(1-e2)*p0;
               f1s = L * (2*pb0 + pb1)/6;
               f2s = L * (2*pb1 + pb0)/6;

            	  p0 = side_loads[1][k1+1];
            	  p1 = side_loads[1][k1+2];
               pb0 = e1*(p1-p0)+p0;
               pb1 = e2*p1+(1-e2)*p0;
               f1n = L * (2*pb0 + pb1)/6;
               f2n = L * (2*pb1 + pb0)/6;

               break;

            	case 2: /* quadratic load variation */
               k1 = 3 * (pside - 1);
               e1 = (m-1)/sl;
               e2 = m/sl;
               eb0 = e1;
               eb1 = (e1+e2)/2.;
               eb2 = e2;

            	  p0 = side_loads[0][k1+1];
            	  p1 = side_loads[0][k1+2];
            	  p2 = side_loads[0][k1+3];
               d1 = p0-2*p1+p2;
               d2 = 3*p0-4*p1+p2;
               pb0 = 2*eb0*eb0*d1-eb0*d2+p0;
               pb1 = 2*eb1*eb1*d1-eb1*d2+p0;
               pb2 = 2*eb2*eb2*d1-eb2*d2+p0;
               f1s = L * (pb0 + 2*pb1)/6;
               f2s = L * (2*pb1 + pb2)/6;

            	  p0 = side_loads[1][k1+1];
            	  p1 = side_loads[1][k1+2];
            	  p2 = side_loads[1][k1+3];
               d1 = p0-2*p1+p2;
               d2 = 3*p0-4*p1+p2;
               pb0 = 2*eb0*eb0*d1-eb0*d2+p0;
               pb1 = 2*eb1*eb1*d1-eb1*d2+p0;
               pb2 = 2*eb2*eb2*d1-eb2*d2+p0;
               f1n = L * (pb0 + 2*pb1)/6;
               f2n = L * (2*pb1 + pb2)/6;

               break;

            	case 3: /* cubic load variation */
               k1 = 4 * (pside - 1);
               e1 = (m-1)/sl;
               e2 = m/sl;
               eb0 = e1;
               eb1 = e1+(e2-e1)/3.;
               eb2 = e1+2*(e2-e1)/3.;
               eb3 = e2;

            	  p0 = side_loads[0][k1+1];
            	  p1 = side_loads[0][k1+2];
            	  p2 = side_loads[0][k1+3];
               p3 = side_loads[0][k1+4];
               d1 = p3-3*p2+3*p1-p0;
               d2 = 2*p0-5*p1+4*p2-p3;
               d3 = 2*p3-9*p2+18*p1-11*p0;
               pb0 = 9/2*eb0*eb0*(eb0*d1+d2)+eb0/2*d3+p0;
               pb1 = 9/2*eb1*eb1*(eb1*d1+d2)+eb1/2*d3+p0;
               pb2 = 9/2*eb2*eb2*(eb2*d1+d2)+eb2/2*d3+p0;
               pb3 = 9/2*eb3*eb3*(eb3*d1+d2)+eb3/2*d3+p0;
               f1s = L * (13*pb0+36*pb1+9*pb2+2*pb3)/120;
               f2s = L * (2*pb0+9*pb1+36*pb2+13*pb3)/120;

            	  p0 = side_loads[1][k1+1];
            	  p1 = side_loads[1][k1+2];
            	  p2 = side_loads[1][k1+3];
               p3 = side_loads[1][k1+4];
               d1 = p3-3*p2+3*p1-p0;
               d2 = 2*p0-5*p1+4*p2-p3;
               d3 = 2*p3-9*p2+18*p1-11*p0;
               pb0 = 9/2*eb0*eb0*(eb0*d1+d2)+eb0/2*d3+p0;
               pb1 = 9/2*eb1*eb1*(eb1*d1+d2)+eb1/2*d3+p0;
               pb2 = 9/2*eb2*eb2*(eb2*d1+d2)+eb2/2*d3+p0;
               pb3 = 9/2*eb3*eb3*(eb3*d1+d2)+eb3/2*d3+p0;
               f1n = L * (13*pb0+36*pb1+9*pb2+2*pb3)/120;
               f2n = L * (2*pb0+9*pb1+36*pb2+13*pb3)/120;

               break;

               default:
                 printf ("\n");
				             printf ("FATAL ERROR: ");
				             printf ("TOO HIGHER DEGREE IN BOUNDARY LOADS");
				             printf ("\n");
				             exit(1);
             }

       f1x = f1s*cosi-f1n*sini;
       f1y = f1s*sini+f1n*cosi;
       i1 = xndof*(n1-1);
       f[lv1][++i1] += f1x;
       f[lv1][++i1] += f1y;

       f2x = f2s*cosi-f2n*sini;
       f2y = f2s*sini+f2n*cosi;
       i1 = xndof*(n2-1);
       f[lv2][++i1] += f2x;
       f[lv2][++i1] += f2y;
					  break;

					}

/*____________________________________________________________________________

							  END OF LOOP ON BOUNDARY SIDES
____________________________________________________________________________*/

		 }

/*____________________________________________________________________________

							  END OF LOOP ON LEVELS
____________________________________________________________________________*/

   }


/*____________________________________________________________________________

									END OF FUNCTION
____________________________________________________________________________*/

}/* end of mlsideload() */
