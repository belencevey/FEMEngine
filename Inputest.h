/**********************[          Inputest.h            ]********************/

#ifndef INPUTEST_H
#define INPUTEST_H

/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: INPUTEST.C                                                 #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #                Routines for reading of estimator type                     #
 #                for stress smoother                                        #
 #                for program SAFE2DUR.                                      #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void   inputest()  : reading of input parameters.                        #
 #                                                                           #
 #****************************************************************************/



#include <cstdio>
#include <stdlib.h>
#include <string.h>



/*******************************************[ USER LIBRARY : INPUTPAR.C ]*****
 #                                                                           #
 #  FUNCTION :  inputest()                                                   #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Reading of estimator type for stress smoothin for program SAFE2DUR.  #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H  :                                                               #
 #    int     printf():   formatted output to stdout                         #
 #    char *  gets()  :   gets a string from stdin                           #
 #                                                                           #
 #  STRING.H :                                                               #
 #    size_t  strlen():   gets length of a string.                           #
 #    char *  strcpy():   copies  one string into another.                   #
 #    char *  strcat():   appends one string to another.                     #
 #    char *  strchr():   search the ocurrance of a character in a string.   #
 #    char *  strncpy():  copies n bytes of one string into another.         #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


	void inputest (int *ESTIMATOR, double *Uexc, int ecode);

/*____________________________________________________________________________

			                 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int *ESTIMATOR; /* estimator type for stress smoothing.                    */
//                 /* ESTIMATOR = 1   nodal averaging                         */
//                 /* ESTIMATOR = 2   side averaging                          */
//                 /* ESTIMATOR = 3   side plus nodal averaging               */
//                 /* ESTIMATOR = 4   ZZ stress smoother                      */
//
// double *Uexc ;  /* energy of exact solution if non zero value given        */
// int     ecode;  /* element code number.                                    */

/*****************************************************************************/


#endif
