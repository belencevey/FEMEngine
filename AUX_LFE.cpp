/**********************[            AUX_LFE.cpp           ]********************/
#include "AUX_LFE.h"

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


void *callocvf (unsigned nel, unsigned elsiz){

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 void   *p;  /* pointer to a newly block of memory               */

/*____________________________________________________________________________

			MEMORY ALLOCATION AND ZEROING WITH VERIFICATION
____________________________________________________________________________*/

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

}







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


void *reallocvf (void *blk, unsigned nel, unsigned elsiz){
/* unsigned nel;      number of data units to be allocated                  */
/* unsigned elsiz;    size of each data unit                                */

    void   *p; // Pointer to a newly block of memory
    //MEMORY REALLOCATION WITH VERIFICATION
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

    return p;
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


void *mallocvf (unsigned siz){

    void   *p; //pointer to a newly block of memory

    //MEMORY ALLOCATION WITH VERIFICATION
    p = malloc (siz);
    if (!p)
    {
        printf("\n INSUFFICIENT MEMORY.....");
        exit(1);
    }

    return p;
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


double dot (double  *x, double *sh, int n){
//  double  *x   :first vector
//  double  *sh  :second vector
//  int      n   :number of elements of each vector

    double d; //Accumulator
    //DOT PRODUCT
    d = 0.0;
    for (int i=1 ; i<=n ; i++)
    {
        d = d+ x[i] * sh[i];
    }
    return d;
}/* end of dot() */


