
/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: AUX_LFE.C                                                  #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #                Auxiliar routines for program LSAFE.                       #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  *callocvf()  : dynamic memory allocation storing zeros.            #
 #  void  *reallocvf() : dynamic memory reallocation.                        #
 #  void  *mallocvf()  : dynamic memory allocation.                          #
 #                                                                           #
 #****************************************************************************/


#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>


/*******************************************[ USER LIBRARY : AUX_LFE.C  ]*****
 #                                                                           #
 #  FUNCTION :  callocvf ()                                                  #
 #                                                                           #
 #  TYPE     :  void *                                                       #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Allocates space for a block of bytes and stores zero in              #
 #      the area. The block of bytes is defined by the number of             #
 #      elements of the block and the element size in bytes.                 #
 #                                                                           #
 #      Return a void * pointer to the newly allocated block or              #
 #      exits program if not enough space exists.                            #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  ALLOC.H  :                                                               #
 #    void *  calloc():   allocates main memory storing zero in the area.    #
 #                                                                           #
 #  STDIO.H  :                                                               #
 #    int     printf():   formatted output to stdout.                        #
 #                                                                           #
 #  STDLIB.H :                                                               #
 #    void    exit () :   terminates program.                                #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void *callocvf (unsigned nel, unsigned elsiz)

/*____________________________________________________________________________

			             PARAMETERS DECLARATION
____________________________________________________________________________*/

 //unsigned nel;
 //unsigned elsiz;

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 void   *p           ;  /* pointer to a newly block of memory               */

/*____________________________________________________________________________

			MEMORY ALLOCATION AND ZEROING WITH VERIFICATION
____________________________________________________________________________*/

/*	 if (((long)nel*(long)elsiz)>65535L)
		{
		  printf("\n MEMORY BLOCK TOO LARGE TO BE ALLOCATED....");
		  printf("\n\n PROGRAM ABORTED");
		  exit(1);
		}
*/
	 p = calloc (nel, elsiz);
	 if (!p)
		{
		  printf("\n INSUFFICIENT MEMORY.....");
		  exit(1);
		}

/*____________________________________________________________________________

		  RETURNS POINTER TO NEWLY MEMORY BLOCK
____________________________________________________________________________*/

		return p;

/*____________________________________________________________________________

			END OF FUNCTION
____________________________________________________________________________*/

}/* end of callocvf() */





/*******************************************[ USER LIBRARY : AUX_LFE.C  ]*****
 #                                                                           #
 #  FUNCTION :  reallocvf ()                                                 #
 #                                                                           #
 #  TYPE     :  void *                                                       #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Reallocates space for a block of bytes.                              #
 #      The block of bytes is defined by a pointer to first element          #
 #      of the block and the new block size in bytes.                        #
 #                                                                           #
 #      Return a void * pointer to the newly reallocated block or            #
 #      exits program if not enough space exists.                            #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  ALLOC.H  :                                                               #
 #    void * realloc():   reallocates a block of bytes in main memory.       #
 #                                                                           #
 #  STDIO.H  :                                                               #
 #    int     printf():   formatted output to stdout.                        #
 #                                                                           #
 #  STDLIB.H :                                                               #
 #    void    exit () :   terminates program.                                #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void *reallocvf (void *blk, unsigned nel, unsigned elsiz)

/*____________________________________________________________________________

			             PARAMETERS DECLARATION
____________________________________________________________________________*/

 //void *blk;
 //unsigned nel;     /* number of data units to be allocated                  */
 //unsigned elsiz;   /* size of each data unit                                */


/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 void   *p           ;  /* pointer to a newly block of memory               */

/*____________________________________________________________________________

			          MEMORY REALLOCATION WITH VERIFICATION
____________________________________________________________________________*/

	 if (((long)nel*(long)elsiz)>65535L)
		{
		  printf("\n MEMORY BLOCK TOO LARGE TO BE REALLOCATED....");
		  printf("\n\n PROGRAM ABORTED");
		  exit(1);
		}

	 if (nel==0)
		{
		  printf("\n NULL MEMORY BLOCK TO BE REALLOCATED....");
		  printf("\n\n PROGRAM ABORTED");
		  exit(1);
		}


    elsiz *= nel;
	 p = realloc (blk, elsiz);
	 if (!p)
		{
		  printf("\n MEMORY REALLOCATION FAILLED.....");
 		  printf("\n\n PROGRAM ABORTED");
		  exit(1);
		}

/*____________________________________________________________________________

		  RETURNS POINTER TO NEWLY MEMORY BLOCK
____________________________________________________________________________*/

		return p;

/*____________________________________________________________________________

			END OF FUNCTION
____________________________________________________________________________*/

}/* end of reallocvf() */






/*******************************************[ USER LIBRARY : AUX_LFE.C  ]*****
 #                                                                           #
 #  FUNCTION :  mallocvf ()                                                  #
 #                                                                           #
 #  TYPE     :  void *                                                       #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Allocates space for a block of bytes.                                #
 #                                                                           #
 #      Return a void * pointer to the newly allocated block or              #
 #      exits program if not enough space exists.                            #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  ALLOC.H  :                                                               #
 #    void *  malloc():   allocates main memory.                             #
 #                                                                           #
 #  STDIO.H  :                                                               #
 #    int     printf():   formatted output to stdout.                        #
 #                                                                           #
 #  STDLIB.H :                                                               #
 #    void    exit () :   terminates program.                                #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  void *mallocvf (unsigned siz)

/*____________________________________________________________________________

			             PARAMETERS DECLARATION
____________________________________________________________________________*/

 //unsigned siz;


/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

		      LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 void   *p           ;  /* pointer to a newly block of memory               */

/*____________________________________________________________________________

		  MEMORY ALLOCATION WITH VERIFICATION
____________________________________________________________________________*/


    p = malloc (siz);
    if (!p)
    {
	printf("\n INSUFFICIENT MEMORY.....");
	exit(1);
    }

/*____________________________________________________________________________

		  RETURNS POINTER TO NEWLY MEMORY BLOCK
____________________________________________________________________________*/

		return p;

/*____________________________________________________________________________

			END OF FUNCTION
____________________________________________________________________________*/

}/* end of mallocvf() */




/*******************************************[ USER LIBRARY : GENFL.C   ]******
 #                                                                           #
 #  FUNCTION :  dot ()                                                       #
 #                                                                           #
 #  TYPE     :  double                                                       #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #    Computation of dot product of two double vectors.                      #
 #    The result is returned as a double value.                              #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  double dot (double  *x, double  *sh, int n)


/*____________________________________________________________________________

			             PARAMETERS DECLARATION
____________________________________________________________________________*/


    //double  *x   ;   /*   first vector                                       */
	//double  *sh  ;   /*   second vector                                      */
	//int      n   ;   /*   number of elements of each vector                  */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 double d    ;   /* accumulator                                             */
 int    i    ;   /* loop counter                                            */

/*____________________________________________________________________________

				DOT PRODUCT
____________________________________________________________________________*/

		  d = 0.0;
		  for (i=1 ; i<=n ; i++) d += x[i] * sh[i];

/*____________________________________________________________________________

				RETURN DOT PRODUCT
____________________________________________________________________________*/

			return d;
/*____________________________________________________________________________

				  END OF FUNCTION
____________________________________________________________________________*/

}/* end of dot() */
