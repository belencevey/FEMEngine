/**********************[            AUX_LFE.h           ]********************/

#ifndef AUX_LFE_H
#define AUX_LFE_H

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
#include <cstdio>
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


  void *callocvf (unsigned nel, unsigned elsiz);

/*____________________________________________________________________________

			                       PARAMETERS
____________________________________________________________________________*/

 /*unsigned nel;*/
 /*unsigned elsiz;*/

/****************************************************************************/






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


  void *reallocvf (void *blk, unsigned nel, unsigned elsiz);

/*____________________________________________________________________________

			                     PARAMETERS
____________________________________________________________________________*/

/* void *blk; */
/* unsigned nel;      number of data units to be allocated                  */
/* unsigned elsiz;    size of each data unit                                */

/****************************************************************************/






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


  void *mallocvf (unsigned siz);

/*____________________________________________________________________________

			                     PARAMETERS
____________________________________________________________________________*/

/* unsigned siz; */

/****************************************************************************/



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


  double dot (double  *x, double *sh, int n);

/*____________________________________________________________________________

			             PARAMETERS DECLARATION
____________________________________________________________________________*/


/*   double  *x   ;      first vector                                       */
/*   double  *sh  ;     second vector                                      */
/*   int      n   ;     number of elements of each vector                  */

/*****************************************************************************/




#endif
