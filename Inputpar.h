/**********************[          Inputpar.h            ]********************/

#ifndef INPUTPAR_H
#define INPUTPAR_H

/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: INPUTPAR.C                                                 #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #                Routines for reading of input parameters                   #
 #                for program SAFE2MIX.                                      #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void   inputpar()  : reading of input parameters.                        #
 #                                                                           #
 #****************************************************************************/



#include <cstdio>
#include <stdlib.h>
#include <string.h>



/*******************************************[ USER LIBRARY : INPUTPAR.C ]*****
 #                                                                           #
 #  FUNCTION :  inputpar()                                                   #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Reading of input parameters for program SAFE2DUR.                    #
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


	void inputpar (int *SOLVER, int ecode);

/*____________________________________________________________________________

			                 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int  *SOLVER ;  /* type of solver                                          */
//                 /*    SOLVER = 1   Gauss Eimination                        */
//                 /*    SOLVER = 2   Orthogonal Decomposition                */
//
// int    ecode ;  /* element code number.                                    */

/******************************************************************************/


#endif
