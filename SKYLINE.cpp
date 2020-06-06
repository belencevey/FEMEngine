/**********************[          SKYLINE.h            ]********************/
#include"SKYLINE.h"


/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME:  SKYLINE.C                                                 #
 #                                                                           #
 #****************************************************************************
 #                                                                           #
 #                                                                           #
 #       Routines for manipulation of skyline storage.                       #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  numeqn()   : assigns sequential equation numbers to each free dof. #
 #  void  expandv2D(): expand a compacted vector of vector into a global one.#
 #  void  expands()  : expand a compacted vector of scalar into a global one.#
 #  void  compressv2D(): compress global vector of vector into compacted one.#
 #  void  compresss(): compress global vector of scalars into compacted one. #
 #  void  profile()  : computes vector of diagonal elements of skyline.      #
 #  void  maxdiff()  : computes effective columns heights.                   #
 #  void  addban()   : assembling element matrix in compacted global matrix. #
 #                                                                           #
 #****************************************************************************/






/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  numeqn ()                                                    #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Assigns sequential equation numbers to each free dof.                #
 #      The status of dof i is obtained from the identification              #
 #      vector id[i] following the next convention:                          #
 #                                                                           #
 #                  id[i]  = 0  free dof                                     #
 #                  id[i] != 0  restrained dof                               #
 #                                                                           #
 #      The number of equations associated with free dofs is                 #
 #      assigned to the parameter neq.                                       #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void numeqn (int numnp, int ndof, int *id, int *neq)

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int     numnp;   /* number of nodal points                                 */
// int     ndof ;   /* maximum number of nodal degrees of freedom             */
// int    *id   ;   /* identification vector [numnp * ndof]                   */
	              /* before processing dof status (0 = free, 1 = restrained)*/
		          /* after  processing numbering of free dof                */
// int    *neq  ;   /* number of equations                                    */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

  int    ne   ;   /* current equation number                                */
  int    numdf;   /* total number of degrees of freedom                     */
  int    i    ;   /* loop counter                                           */

/*____________________________________________________________________________

			  CONSTANTS DEFINITION
____________________________________________________________________________*/


			 numdf = numnp * ndof;

/*____________________________________________________________________________

		  ASSIGN EQUATION NUMBER TO EACH FREE DOF
____________________________________________________________________________*/


  ne = 0;

  for ( i=1 ; i<=numdf ; i++)
	  id[i] = (!id[i]) ? ++ne : 0;

/*____________________________________________________________________________

		  ASSIGN TOTAL NUMBER OF EQUATIONS TO PARAMETER NEQ
____________________________________________________________________________*/


				 *neq = ne;

/*____________________________________________________________________________

				  END OF FUNCTION
____________________________________________________________________________*/

} /* end of numeqn() */






/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  expandv2D()                                                  #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Expand a compacted vector dn of 2D vectors into a complete vector u. #
 #      The number of dof i is obtained from the identification              #
 #      vector id[i] following the next convention:                          #
 #                                                                           #
 #                  id[i]  = 0  u[i] = 0                                     #
 #                  id[i] != 0  u[i] = dn[id[i]]                             #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void expandv2D (double *dn, double *u, int *id, int numnp, int xndof, int *idb, double *bn_cs, double *bn_sn)

/*____________________________________________________________________________

								PARAMETERS DECLARATION
____________________________________________________________________________*/

// double *dn   ;  /* compacted vector to be expanded                         */
// double *u    ;  /* expanded vector                                         */
// int    *id   ;  /* identification vector of free dof [numnp * xndof]       */
// int     numnp;  /* number of nodal points                                  */
// int     xndof;  /* maximum number of nodal degrees of freedom              */
// int    *idb    ;  /* identification vector of boundary nodes [numnp[0]]    */
// double *bn_cs  ;  /* boundary node cosine of local axes rotation angle     */
// double *bn_sn  ;  /* boundary node sine   of local axes rotation angle     */

/*____________________________________________________________________________

								 BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

							  LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

  int    i     ;  /* loop counter                                           */
  int    n1    ;  /* boundary node number                                   */
  int    i1, i2;  /* auxiliar pointers for vector u[]                       */
  double cos1  ;  /* cosine of rotation of boundary node                    */
  double sin1  ;  /* sine of rotation of boundary node                      */
  double u1, u2;  /* vector components in local system of boundary node.    */


/*____________________________________________________________________________

						  EXPAND VECTOR dn INTO VECTOR u
____________________________________________________________________________*/


  for ( i=1 ; i<=numnp ; i++)
	  {
      n1 = idb[i];
      i2 = xndof*i;
      i1 = i2-1;

      if (n1)
        {
          cos1 = bn_cs[n1];
          sin1 = bn_sn[n1];

          u1 = (id[i1]) ? dn[id[i1]] : 0;
          u2 = (id[i2]) ? dn[id[i2]] : 0;

          u[i1] = u1*cos1 - u2*sin1;
          u[i2] = u1*sin1 + u2*cos1;
        }
      else
        {
          u[i1] = (id[i1]) ? dn[id[i1]] : 0;
          u[i2] = (id[i2]) ? dn[id[i2]] : 0;
        }

     }


/*____________________________________________________________________________

									END OF FUNCTION
____________________________________________________________________________*/

} /* end of expandv2D() */






/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  expands()                                                    #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Expand a compacted vector dn of scalars into a complete vector u.    #
 #      The number of dof i is obtained from the identification              #
 #      vector id[i] following the next convention:                          #
 #                                                                           #
 #                  id[i]  = 0  u[i] = 0                                     #
 #                  id[i] != 0  u[i] = dn[id[i]]                             #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void expands (double *dn, double *u, int *id, int numnp)

/*____________________________________________________________________________

								PARAMETERS DECLARATION
____________________________________________________________________________*/

// double *dn   ;  /* compacted vector to be expanded                         */
// double *u    ;  /* expanded vector                                         */
// int    *id   ;  /* identification vector of free dof [numnp * xndof]       */
// int     numnp;  /* number of nodal points                                  */

/*____________________________________________________________________________

								 BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

							  LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

  int    i     ;  /* loop counter                                           */

/*____________________________________________________________________________

						  EXPAND VECTOR dn INTO VECTOR u
____________________________________________________________________________*/


  for ( i=1 ; i<=numnp ; i++) u[i] = (id[i]) ? dn[id[i]] : 0;


/*____________________________________________________________________________

									END OF FUNCTION
____________________________________________________________________________*/

} /* end of expandvs() */







/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  compressv2D()                                                #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Compress a global vector u of 2D vectors into a compacted vector dn. #
 #      Rotates global components to local node system.                      #
 #      The number of dof i is obtained from the identification              #
 #      vector id[i] following the next convention:                          #
 #                                                                           #
 #                  id[i]  = 0  u[i] = 0                                     #
 #                  id[i] != 0  u[i] = dn[id[i]]                             #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void compressv2D (double *dn, double *u, int *id, int numnp, int xndof, int *idb, double *bn_cs, double *bn_sn)


/*____________________________________________________________________________

								PARAMETERS DECLARATION
____________________________________________________________________________*/

// double *dn   ;  /* compacted vector                                        */
// double *u    ;  /* global vector to be compressed                          */
// int    *id   ;  /* identification vector of free dof [numnp * xndof]       */
// int     numnp;  /* number of nodal points                                  */
// int     xndof;  /* maximum number of nodal degrees of freedom              */
// int    *idb    ;  /* identification vector of boundary nodes [numnp[0]]    */
// double *bn_cs  ;  /* boundary node cosine of local axes rotation angle     */
// double *bn_sn  ;  /* boundary node sine   of local axes rotation angle     */



/*____________________________________________________________________________

								 BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

							  LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

  int    i     ;  /* loop counter                                           */
  int    n1    ;  /* boundary node number                                   */
  int    i1, i2;  /* auxiliar pointers for vector u[]                       */
  double cos1  ;  /* cosine of rotation of boundary node                    */
  double sin1  ;  /* sine of rotation of boundary node                      */
  double ux, uy;  /* vector components in global coordinates                */
  double u1, u2;  /* vector components in local system of boundary node.    */

/*____________________________________________________________________________

						  COMPRESS VECTOR u INTO VECTOR dn
____________________________________________________________________________*/


  for (i = 1; i <= numnp; i++)
	  {
      n1 = idb[i];
      i2 = xndof*i;
      i1 = i2-1;

      if (n1)
      {
          cos1 = bn_cs[n1];
          sin1 = bn_sn[n1];

          ux = u[i1];
          uy = u[i2];

          u1 =  ux*cos1 + uy*sin1;
          u2 = -ux*sin1 + uy*cos1;
        }
      else
        {
          u1 = u[i1];
          u2 = u[i2];
        }

      if (id[i1]) dn[id[i1]] = u1;
      if (id[i2]) dn[id[i2]] = u2;

     }

/*____________________________________________________________________________

									END OF FUNCTION
____________________________________________________________________________*/

} /* end of compressv2D() */







/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  compresss()                                                  #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Compress a global vector u of scalars into a compacted vector dn.    #
 #      Rotates global components to local node system.                      #
 #      The number of dof i is obtained from the identification              #
 #      vector id[i] following the next convention:                          #
 #                                                                           #
 #                  id[i]  = 0  u[i] = 0                                     #
 #                  id[i] != 0  u[i] = dn[id[i]]                             #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void compresss (double *dn, double *u, int *id, int numnp)


/*____________________________________________________________________________

								PARAMETERS DECLARATION
____________________________________________________________________________*/

// double *dn   ;  /* compacted vector                                        */
// double *u    ;  /* global vector to be compressed                          */
// int    *id   ;  /* identification vector of free dof [numnp * xndof]       */
// int     numnp;  /* number of nodal points                                  */

/*____________________________________________________________________________

								 BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

							  LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

  int    i     ;  /* loop counter                                           */

/*____________________________________________________________________________

						  COMPRESS VECTOR u INTO VECTOR dn
____________________________________________________________________________*/


  for ( i=1 ; i<=numnp ; i++) if (id[i]) dn[id[i]] = u[i];

/*____________________________________________________________________________

									END OF FUNCTION
____________________________________________________________________________*/

} /* end of compresss() */







/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  profile ()                                                   #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computes effective column heights of skyline.                        #
 #      These heights are the maximum difference between the                 #
 #      equation number of each dof and the least element                    #
 #      equation number, for all the elements that contains                  #
 #      that dof. The height of the dof i is stored in the                   #
 #      location maxa[i+1].                                                  #
 #      Computes adresses of diagonal elements in skyline and                #
 #      stores them in the same array maxa[].                                #
 #      The elements of a column of the skyline are sequentially             #
 #      numbered from bottom to top, and the position number of              #
 #      of the diagonal elements of the matrix (bottom elements              #
 #      of columns) in the skyline is computed from the column               #
 #      heights as:                                                          #
 #                                                                           #
 #              maxa[1] = 1                                                  #
 #              maxa[2] = 2                                                  #
 #                                                                           #
 #              for i>2  maxa[i] = maxa[i] + maxa[i-1] + 1;                  #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


 void profile (int neq, int numel, int sides, int bsides, int ecode, int xndof, int *id, int *elm_sides,
               int *side_nodes, int *lm, int *nwk, int *maxa)

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int  neq       ;  /* number of equations                                   */
// int  numel     ;  /* number of elements                                    */
// int  sides     ;  /* number of sides.                                      */
// int  bsides    ;  /* number of boundary sides.                             */
// int  ecode     ;  /* element code number                                   */
// int  xndof     ;  /* maximum number of nodal degrees of freedom            */
// int *id        ;  /* equation number of free dof [numnp * ndof]            */
						             /* ( >0 eqn. number, =0 restrained)                      */
// int *side_nodes;  /* nodes of each side       [sides * xnsn]               */
// int *elm_sides ;  /* sides of each element    [numel * xnes]               */
// int *lm        ;  /* equation number of element dof [nen*ndof]             */
// int *nwk       ;  /* number of matrix elements                             */
// int *maxa      ;  /* pointer of diagonal elements for skyline [neq+1]      */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

  int    neq1 ;   /* number of equations plus one (neq + 1)                 */
  int    i    ;   /* loop counters                                          */

/*____________________________________________________________________________

	 FUNCTION DECLARATIONS
____________________________________________________________________________*/

  void maxdiff(int neq, int numel, int sides, int bsides, int ecode, int xndof,
               int* id,	int* elm_sides, int* side_nodes, int* lm, int* maxa);

/*____________________________________________________________________________

	 SEARCH FOR MAXIMUM DIFFERENCE IN DOF NUMBERING
____________________________________________________________________________*/

  maxdiff (neq, numel, sides, bsides, ecode, xndof, id, elm_sides,
           side_nodes, lm, maxa);

/*________________________________________  __________________________________

			  POINTER OF MAIN DIAGONAL
 ___________________________________________________________________________*/

  neq1 = neq + 1;

  maxa[1] = 1;
  maxa[2] = 2;

  if ( neq > 1 )
	 {
		for ( i=3 ; i<=neq1 ; i++)
			{
			  maxa[i] = maxa[i] + maxa[i-1] + 1;
			}
	 }

/*__________________________________________________________________________

	ASSIGN TOTAL NUMBER OF MATRIX ELEMENTS TO PARAMETER NWK
 ___________________________________________________________________________*/


  *nwk = maxa[neq1] - 1;

/*____________________________________________________________________________

				  END OF FUNCTION
____________________________________________________________________________*/

} /* end of profile() */






/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  maxdiff ()                                                   #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Computes effective column heights of skyline.                        #
 #      These heights are the maximum difference between the                 #
 #      equation number of each dof and the least element                    #
 #      equation number, for all the elements that contains                  #
 #      that dof. The height of the dof i is stored in the                   #
 #      location maxa[i+1].                                                  #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/

#include <math.h>

 void maxdiff (int neq, int numel, int sides, int bsides, int ecode, int xndof, int *id, int *elm_sides,
               int *side_nodes, int *lm, int *maxa)

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int  neq       ;  /* number of equations                                   */
// int  numel     ;  /* number of elements                                    */
// int  sides     ;  /* number of sides.                                      */
// int  bsides    ;  /* number of boundary sides.                             */
// int  ecode     ;  /* element code number                                   */
// int  xndof     ;  /* maximum number of nodal degrees of freedom            */
// int *id        ;  /* equation number of free dof [numnp * ndof]            */
						 /* ( >0 eqn. number, =0 restrained)                      */
// int *side_nodes;  /* nodes of each side       [sides * xnsn]               */
// int *elm_sides ;  /* sides of each element    [numel * xnes]               */
// int *lm        ;  /* equation number of element dof [nen*ndof] [nsn*ndof]  */
// int *maxa      ;  /* pointer of diagonal elements for skyline [neq+1]      */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

  int    neq1 ;   /* number of equations plus one (neq + 1)                 */
  int    nen  ;   /* number of element nodes                                */
  int    numed;   /* number of element degrees of freedom (nen * ndof)      */
  int    nsn  ;   /* number of side nodes                                   */
  int    numsd;   /* number of side degrees of freedom (nsn * ndof)         */
  int    ls   ;   /* least element dof number                               */
  int    dif  ;   /* difference between element dof number and least number */
  int    i,j,k;   /* loop counters                                          */
  int    j1   ;   /* auxiliar pointer for vector id[]                       */
  int    k1   ;   /* auxiliar pointer for vector lm[]                       */
  int    L[5] ;   /* element sides                                          */
  int    n[10] ;  /* element or side nodes                                  */


/*____________________________________________________________________________

									CONSTANTS DEFINITION
____________________________________________________________________________*/


 neq1 = neq + 1;

/*____________________________________________________________________________

	 SEARCH FOR MAXIMUM DIFFERENCE IN DOF NUMBERING
____________________________________________________________________________*/


 switch (ecode)
		 {
/*____________________________________________________________________________

	     ELEMENT BY ELEMENT SEARCH FOR LINEAR TRIANGLES
____________________________________________________________________________*/


   case 1: /* linear triangle */

			  nen = 3;
			  nsn = 4;
			  numed = nen * xndof;

			  /*----LOOP ON ELEMENTS----*/

			  for ( i=1 ; i<=numel ; i++)
				  {
					 ls = neq1;

					 L[1] = elm_sides [3*i-2];
					 L[2] = elm_sides [3*i-1];
					 L[3] = elm_sides [3*i];

					 n[1] = (L[1]>0) ? side_nodes[nsn*(L[1]-1)+1] : side_nodes[nsn*(-L[1]-1)+2];
					 n[2] = (L[2]>0) ? side_nodes[nsn*(L[2]-1)+1] : side_nodes[nsn*(-L[2]-1)+2];
					 n[3] = (L[3]>0) ? side_nodes[nsn*(L[3]-1)+1] : side_nodes[nsn*(-L[3]-1)+2];

					 k1 = 0;
					 for (j=1 ; j<=nen ; j++)
						 {
							j1 = (n[j]-1) * xndof;
							for (k=1 ; k<=xndof ; k++)
								{
								  lm[++k1] = id[j1 + k];
								  if (lm[k1] >0  && lm[k1] <ls) ls = lm[k1];
								}
						 }


					 /*----UPDATES MAXIMUM DIFFERENCE OF EACH DOF----*/

					 for (j=1 ; j<=numed ; j++)
						 {
							if (lm[j] > 0)
							  {
								 dif = lm[j] - ls;
								 if ( dif > maxa[lm[j]+1] ) maxa[lm[j]+1] = dif;
							  }
						 }

				  } /*----END OF LOOP ON LINEAR TRIANGLES---*/
			  break;

/*____________________________________________________________________________

		   ELEMENT BY ELEMENT SEARCH FOR LINEAR TRIANGLES CUBIC BEZIER
  ____________________________________________________________________________*/


   case 8: /* linear triangle cubic bezier */

	   nen = 3;
	   nsn = 6;
	   numed = nen * xndof;

	   /*----LOOP ON ELEMENTS----*/

	   for (i = 1; i <= numel; i++)
	   {
		   ls = neq1;

		   L[1] = elm_sides[3 * i - 2];
		   L[2] = elm_sides[3 * i - 1];
		   L[3] = elm_sides[3 * i];

		   n[1] = (L[1]>0) ? side_nodes[nsn*(L[1] - 1) + 1] : side_nodes[nsn*(-L[1] - 1) + 2];
		   n[2] = (L[2]>0) ? side_nodes[nsn*(L[2] - 1) + 1] : side_nodes[nsn*(-L[2] - 1) + 2];
		   n[3] = (L[3]>0) ? side_nodes[nsn*(L[3] - 1) + 1] : side_nodes[nsn*(-L[3] - 1) + 2];

		   k1 = 0;
		   for (j = 1; j <= nen; j++)
		   {
			   j1 = (n[j] - 1) * xndof;
			   for (k = 1; k <= xndof; k++)
			   {
				   lm[++k1] = id[j1 + k];
				   if (lm[k1] >0 && lm[k1] <ls) ls = lm[k1];
			   }
		   }


		   /*----UPDATES MAXIMUM DIFFERENCE OF EACH DOF----*/

		   for (j = 1; j <= numed; j++)
		   {
			   if (lm[j] > 0)
			   {
				   dif = lm[j] - ls;
				   if (dif > maxa[lm[j] + 1]) maxa[lm[j] + 1] = dif;
			   }
		   }

	   } /*----END OF LOOP ON LINEAR TRIANGLES CUBIC BEZIER---*/
	   break;

/*____________________________________________________________________________

       ELEMENT BY ELEMENT SEARCH FOR QUADRATIC TRIANGLES
____________________________________________________________________________*/



			case 2: /* quadratic triangle */

			  nen = 6;
			  numed = nen * xndof;

			  /*----LOOP ON ELEMENTS----*/

			  for ( i=1 ; i<=numel ; i++)
				  {
					 ls = neq1;

					 L[1] = elm_sides [3*i-2];
					 L[2] = elm_sides [3*i-1];
					 L[3] = elm_sides [3*i];

					 n[1] = (L[1]>0) ? side_nodes[3*L[1]-2] : side_nodes[-3*L[1]];
					 n[2] = (L[2]>0) ? side_nodes[3*L[2]-2] : side_nodes[-3*L[2]];
					 n[3] = (L[3]>0) ? side_nodes[3*L[3]-2] : side_nodes[-3*L[3]];
					 n[4] = (L[1]>0) ? side_nodes[3*L[1]-1] : side_nodes[-3*L[1]-1];
					 n[5] = (L[2]>0) ? side_nodes[3*L[2]-1] : side_nodes[-3*L[2]-1];
					 n[6] = (L[3]>0) ? side_nodes[3*L[3]-1] : side_nodes[-3*L[3]-1];


					 k1 = 0;
					 for (j=1 ; j<=nen ; j++)
						 {
							j1 = (n[j]-1) * xndof;
							for (k=1 ; k<=xndof ; k++)
								{
								  lm[++k1] = id[j1 + k];
								  if (lm[k1] >0  && lm[k1] <ls) ls = lm[k1];
								}
						 }


					 /*----UPDATES MAXIMUM DIFFERENCE OF EACH DOF----*/

					 for (j=1 ; j<=numed ; j++)
						 {
							if (lm[j] > 0)
							  {
								 dif = lm[j] - ls;
								 if ( dif > maxa[lm[j]+1] ) maxa[lm[j]+1] = dif;
							  }
						 }

				  } /*----END OF LOOP ON QUADRATIC TRIANGLES---*/
			  break;

/*____________________________________________________________________________

			  SIDE BY SIDE SEARCH FOR LINEAR MIXED TRIANGLES
____________________________________________________________________________*/

     case 3: /* linear mixed triangle */
	 case 7: /* TRI3 plane linear triangle by side */

			  /*----LOOP ON BOUNDARY SIDES----*/

			  nsn = 3;
			  numsd = nsn * xndof;

			  for ( i=1 ; i<=bsides ; i++)
			   {
			      ls = neq1;

                  n[1] = side_nodes[4*i-3];
		          n[2] = side_nodes[4*i-2];
		          n[3] = side_nodes[4*i-1];

                  k1 = 0;
			      for (j=1 ; j<=nsn ; j++)
			       {
		              j1 = (n[j]-1) * xndof;
				   	  for (k=1 ; k<=xndof ; k++)
			           {
				          lm[++k1] = id[j1 + k];
					      if (lm[k1] >0  && lm[k1] <ls) ls = lm[k1];
					   }
			  	   }

				 /*----UPDATES MAXIMUM DIFFERENCE OF EACH DOF----*/

			 	  for (j=1 ; j<=numsd ; j++)
				   {
					  if (lm[j] > 0)
					    {
					  	  dif = lm[j] - ls;
						  if ( dif > maxa[lm[j]+1] ) maxa[lm[j]+1] = dif;
					    }
				   }

				  } /*----END OF LOOP ON BOUNDARY SIDES---*/


			  /*----LOOP ON INTERIOR SIDES----*/

			  nsn = 4;
			  numsd = nsn * xndof;

			  for ( i=bsides+1 ; i<=sides ; i++)
			   {
	               ls = neq1;

				   n[1] = side_nodes[4*i-3];
				   n[2] = side_nodes[4*i-2];
				   n[3] = side_nodes[4*i-1];
                   n[4] = side_nodes[4*i];

				   k1 = 0;
				   for (j=1 ; j<=nsn ; j++)
				    {
				       j1 = (n[j]-1) * xndof;
					   for (k=1 ; k<=xndof ; k++)
						{
					       lm[++k1] = id[j1 + k];
						  if (lm[k1] >0  && lm[k1] <ls) ls = lm[k1];
						}
					 }

				   /*----UPDATES MAXIMUM DIFFERENCE OF EACH DOF----*/

				   for (j=1 ; j<=numsd ; j++)
				    {
				        if (lm[j] > 0)
						 {
						    dif = lm[j] - ls;
						    if ( dif > maxa[lm[j]+1] ) maxa[lm[j]+1] = dif;
					     }
			 	    }

				  } /*----END OF LOOP ON INTERIOR SIDES---*/

			  break;

/*____________________________________________________________________________

			  SIDE BY SIDE SEARCH FOR MIXED MORLEY TRIANGLES
____________________________________________________________________________*/

     case 4: /* mixed Morley triangle */
     case 5:

			  /*----LOOP ON BOUNDARY SIDES----*/

			  nsn = 6;
			  numsd = nsn * xndof;

			  for ( i=1 ; i<=bsides ; i++)
			   {
			      ls = neq1;

                  n[1] = side_nodes[9*i-8];
		          n[2] = side_nodes[9*i-7];
		          n[3] = side_nodes[9*i-6];
                  n[4] = abs(side_nodes[9*i-5]);
		          n[5] = abs(side_nodes[9*i-4]);
		          n[6] = abs(side_nodes[9*i-3]);

                  k1 = 0;
			      for (j=1 ; j<=nsn ; j++)
			       {
		              j1 = (n[j]-1) * xndof;
				   	  for (k=1 ; k<=xndof ; k++)
			           {
				          lm[++k1] = id[j1 + k];
					      if (lm[k1] >0  && lm[k1] <ls) ls = lm[k1];
					   }
			  	   }

				 /*----UPDATES MAXIMUM DIFFERENCE OF EACH DOF----*/

			 	  for (j=1 ; j<=numsd ; j++)
				   {
					  if (lm[j] > 0)
					    {
					  	  dif = lm[j] - ls;
						  if ( dif > maxa[lm[j]+1] ) maxa[lm[j]+1] = dif;
					    }
				   }

				  } /*----END OF LOOP ON BOUNDARY SIDES---*/


			  /*----LOOP ON INTERIOR SIDES----*/

			  nsn = 9;
			  numsd = nsn * xndof;

			  for ( i=bsides+1 ; i<=sides ; i++)
			   {
	               ls = neq1;

				   n[1] = side_nodes[9*i-8];
				   n[2] = side_nodes[9*i-7];
				   n[3] = side_nodes[9*i-6];
                   n[4] = abs(side_nodes[9*i-5]);
				   n[5] = abs(side_nodes[9*i-4]);
				   n[6] = abs(side_nodes[9*i-3]);
				   n[7] = side_nodes[9*i-2];
                   n[8] = abs(side_nodes[9*i-1]);
                   n[9] = abs(side_nodes[9*i]);

				   k1 = 0;
				   for (j=1 ; j<=nsn ; j++)
				    {
				       j1 = (n[j]-1) * xndof;
					   for (k=1 ; k<=xndof ; k++)
						{
					       lm[++k1] = id[j1 + k];
						  if (lm[k1] >0  && lm[k1] <ls) ls = lm[k1];
						}
					 }

				   /*----UPDATES MAXIMUM DIFFERENCE OF EACH DOF----*/

				   for (j=1 ; j<=numsd ; j++)
				    {
				        if (lm[j] > 0)
						 {
						    dif = lm[j] - ls;
						    if ( dif > maxa[lm[j]+1] ) maxa[lm[j]+1] = dif;
					     }
			 	    }

				  } /*----END OF LOOP ON INTERIOR SIDES---*/

			  break;

/*____________________________________________________________________________

			  SIDE BY SIDE SEARCH FOR LINEAR MIXED PLATE TRIANGLES
____________________________________________________________________________*/

     case 6: /* linear mixed plate triangle */

			  /*----LOOP ON BOUNDARY SIDES----*/

			  nsn = 3;
			  numsd = nsn * xndof;

			  for ( i=1 ; i<=bsides ; i++)
			   {
			      ls = neq1;

                  n[1] = side_nodes[4*i-3];
		          n[2] = side_nodes[4*i-2];
		          n[3] = side_nodes[4*i-1];

                  k1 = 0;
			      for (j=1 ; j<=nsn ; j++)
			       {
		              j1 = (n[j]-1) * xndof;
				   	  for (k=1 ; k<=xndof ; k++)
			           {
				          lm[++k1] = id[j1 + k];
					      if (lm[k1] >0  && lm[k1] <ls) ls = lm[k1];
					   }
			  	   }

				 /*----UPDATES MAXIMUM DIFFERENCE OF EACH DOF----*/

			 	  for (j=1 ; j<=numsd ; j++)
				   {
					  if (lm[j] > 0)
					    {
					  	  dif = lm[j] - ls;
						  if ( dif > maxa[lm[j]+1] ) maxa[lm[j]+1] = dif;
					    }
				   }

				  } /*----END OF LOOP ON BOUNDARY SIDES---*/


			  /*----LOOP ON INTERIOR SIDES----*/

			  nsn = 4;
			  numsd = nsn * xndof;

			  for ( i=bsides+1 ; i<=sides ; i++)
			   {
	               ls = neq1;

				   n[1] = side_nodes[4*i-3];
				   n[2] = side_nodes[4*i-2];
				   n[3] = side_nodes[4*i-1];
                   n[4] = side_nodes[4*i];

				   k1 = 0;
				   for (j=1 ; j<=nsn ; j++)
				    {
				       j1 = (n[j]-1) * xndof;
					   for (k=1 ; k<=xndof ; k++)
						{
					       lm[++k1] = id[j1 + k];
						  if (lm[k1] >0  && lm[k1] <ls) ls = lm[k1];
						}
					 }

				   /*----UPDATES MAXIMUM DIFFERENCE OF EACH DOF----*/

				   for (j=1 ; j<=numsd ; j++)
				    {
				        if (lm[j] > 0)
						 {
						    dif = lm[j] - ls;
						    if ( dif > maxa[lm[j]+1] ) maxa[lm[j]+1] = dif;
					     }
			 	    }

				  } /*----END OF LOOP ON INTERIOR SIDES---*/

			  break;

/*____________________________________________________________________________

				  END OF SELECT CASE ECODE
____________________________________________________________________________*/

         }

/*____________________________________________________________________________

				  END OF FUNCTION
____________________________________________________________________________*/

} /* end of maxdiff() */






/*******************************************[ USER LIBRARY : SKYLINE.C ]******
 #                                                                           #
 #  FUNCTION :  addban ()                                                    #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Assembling of an upper triangular element matrix into                #
 #      the compacted global matrix.                                         #
 #                                                                           #
 #       (translation to C language of subroutine ADDBAN from                #
 #        K.J.Bathe, "Finite Element Procedures in Engineering               #
 #        Analysis", Prentice Hall, 1982.)                                   #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void addban (double *A, int *maxa, double *sv, int *lm, int numrw)

/*____________________________________________________________________________

				 PARAMETERS DECLARATION
____________________________________________________________________________*/

// double  *A   ;   /* global matrix stored by skyline [nwk]                  */
// int     *maxa;   /* pointer of diagonal elements for skyline [neq+1]       */
// double  *sv  ;   /* upper triagular element matrix stored by rows          */
// int     *lm  ;   /* equation numbers of element dof [nen*ndof]             */
// int     numrw;   /* number of element square matrix rows/columns           */


/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

  int    ndi  ;   /* pointer of diagonal coefficient for element matrix     */
  int    i    ;   /* loop counter of element matrix rows                    */
  int    j    ;   /* loop counter of element matrix columns                 */
  int    ii   ;   /* global row ii associated to element matrix row i       */
  int    jj   ;   /* global column jj associated to element matrix column j */
  int    ij   ;   /* difference (ii - jj). ij>0  implies upper triangular   */
		  /* localization in global matrix                          */
  int    mi   ;   /* diagonal pointer for global column ii                  */
  int    kk   ;   /* pointer for global position (ii,jj) in skyline         */
  int    ks   ;   /* pointer for element matrix when i<j                    */
  int    kss  ;   /* pointer for local position (i,j) in element matrix     */


/*____________________________________________________________________________

			 LOOP ON ELEMENT MATRIX ROWS
____________________________________________________________________________*/

  ndi = 0;
  for ( i=1 ; i<=numrw ; i++)
	  {
/*____________________________________________________________________________

			VERIFICATION OF RESTRAINED ROW DOF
____________________________________________________________________________*/

	ii = lm[i];
	if (ii > 0)
	  {
/*____________________________________________________________________________

		SET GLOBAL DIAGONAL POINTER FOR COLUMN ii
____________________________________________________________________________*/

		 mi = maxa[ii];
		 ks = i;

/*____________________________________________________________________________

			 LOOP ON ELEMENT MATRIX COLUMNS
____________________________________________________________________________*/

		 for ( j=1 ; j<=numrw ; j++)
			 {
/*____________________________________________________________________________

		 VERIFICATION OF RESTRAINED COLUMN DOF
____________________________________________________________________________*/

		  jj = lm[j];
		  if (jj > 0)
			 {
/*____________________________________________________________________________

		VERIFICATION OF LOCALIZATION IN UPPER TRIANGULAR GLOBAL MATRIX
____________________________________________________________________________*/

				ij = ii - jj;
				if (ij >= 0)
			{
/*____________________________________________________________________________

	  ASSSEMBLES LOCAL COEFFICIENT (I,J) IN GLOBAL LOCALIZATION (II,JJ)
____________________________________________________________________________*/

			  kk = mi + ij;
			  kss = ks;
			  if (j >= i) kss = j + ndi;
			  A[kk] += sv[kss];
/*____________________________________________________________________________

	 END OF VERIFICATION OF LOCALIZATION IN UPPER TRIANGULAR GLOBAL MATRIX
____________________________________________________________________________*/

			}
/*____________________________________________________________________________

			 END OF VERIFICATION OF RESTRAINED COLUMN DOF
____________________________________________________________________________*/

			 }
		  ks += numrw - j;
/*____________________________________________________________________________

		END OF LOOP ON ELEMENT MATRIX COLUMNS
____________________________________________________________________________*/

			 }
/*____________________________________________________________________________

			 END OF VERIFICATION OF RESTRAINED ROW DOF
____________________________________________________________________________*/

	  }
/*____________________________________________________________________________

		SET POINTER OF DIAGONAL COEFFICIENTS TO A NEW ELEMENT MATRIX ROW
____________________________________________________________________________*/

	ndi += (numrw - i);

/*____________________________________________________________________________

		 END OF LOOP ON ELEMENT MATRIX ROWS
____________________________________________________________________________*/

	  }

} /* end of addban() */





