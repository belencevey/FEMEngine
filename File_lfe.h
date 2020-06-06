/**********************[           FILE_LFE.h           ]********************/

#ifndef FILE_LFE_H
#define FILE_LFE_H

/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: FILE_LFE.C                                                 #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #                Routines for manipulation of input/output                  #
 #                files of program SAFE2D.                                   #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void   iofile()    : opening of input/output files.                      #
 #  FILE  *fopenvf()   : opens a file with verification.                     #
 #  void   prtdis()    : print displacements in output file.                 #
 #                                                                           #
 #****************************************************************************/



#include <cstdio>
#include <stdlib.h>
#include <string>
#include<fstream>
#include<iostream>



/*******************************************[ USER LIBRARY : FILE_LFE.C ]*****
 #                                                                           #
 #  FUNCTION :  iofile ()                                                    #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Open input/output files for program SAFE2D.                          #
 #                                                                           #
 #      Gets input filename ######[.###] from command line arguments.        #
 #      The input file must be in the same directory of the program.         #
 #      The root  filename is  the name  of the input file  without any      #
 #      extension. A maximum size of six characters length is assumed.       #
 #      The extension is optional and is ignored in forming the root         #
 #      filename.                                                            #
 #                                                                           #
 #      The output filenames are formed by the root filename plus two        #
 #      digits associated with the corresponding level of discretization.    #
 #      The extension of the output files are related to the postprocessor   #
 #      program used.                                                        #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED LIBRARY FUNCTIONS USED ]--------------------------------#
 #                                                                           #
 #                                                                           #
 #  FILE_LFE.C :                                                             #
 #    FILE  *fopenvf()     : opens a file with verification.                 #
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
 #----[ STANDARD LIBRARY STRUCTURES USED ]-----------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H  :                                                               #
 #    struct  FILE {} :   File control structure for streams.                #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/
    void iofile (std::string *filename, std::ifstream &fpi, std::string *forname, std::ofstream &fpo, std::ofstream &fpiter);
	//void iofile (char *filename, FILE **fpi, char *forname, FILE **fpo, FILE **fpiter);

/*____________________________________________________________________________

			                 PARAMETERS DECLARATION
____________________________________________________________________________*/

// char   *filename   ;  /* complete input file name                        */
// FILE  **fpi        ;  /* pointer input  file                             */
// char   *forname    ;  /* complete output file root name                  */
// FILE  **fpo        ;  /* pointer results file                            */
// FILE  **fpiter     ;  /* pointer iteration steps file                    */


/*****************************************************************************/









/*******************************************[ USER LIBRARY : FILE_LFE.C ]*****
 #                                                                           #
 #  FUNCTION :  fopenvf ()                                                   #
 #                                                                           #
 #  TYPE     :  FILE *                                                       #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Opens a file with verification.                                      #
 #                                                                           #
 #      Return a FILE * pointer to the newly open file if successful,        #
 #      else exits program.                                                  #
 #                                                                           #
 #      The type of opening is controlled by the char parameter mode.        #
 #                                                                           #
 #      Codes for open types:                                                #
 #                                                                           #
 #                'r' : open for reading                                     #
 #                'w' : create for writing                                   #
 #                'a' : append, open for writing at end-of-file              #
 #                '+' : add symbol to allow read/write access                #
 #                'b' : open in binary mode                                  #
 #                't' : open in text mode                                    #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H  :                                                               #
 #    FILE *  fopen ():   opens a stream.                                    #
 #    int     printf():   formatted output to stdout                         #
 #                                                                           #
 #  STDLIB.H :                                                               #
 #    void    exit () :   terminates program.                                #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY STRUCTURES USED ]-----------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H  :                                                               #
 #    struct  FILE {} :   File control structure for streams.                #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  //FILE *fopenvf (char *finame, char *mode);

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// char   *finame      ;  /* name of the file to be open                    */
// char   *mode        ;  /* type of opening                                */


/****************************************************************************/








/*******************************************[ USER LIBRARY : FILE_LFE.C ]*****
 #                                                                           #
 #  FUNCTION :  prtdis ()                                                    #
 #                                                                           #
 #  TYPE     :  FILE *                                                       #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Print displacements in output file.                                  #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H  :                                                               #
 #    int    fprintf():   sends formatted output to a stream.                #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY STRUCTURES USED ]-----------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H  :                                                               #
 #    struct  FILE {} :   File control structure for streams.                #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/


  //void prtdis (FILE *fpo, int numnp, int ndof, int *id, double *dn);
    void prtdis (std::ifstream &fpo, int numnp, int ndof, int *id,double *dn);
/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE   *fpo  ;   /* pointer to output file                               */
// int     numnp;   /* number of nodal points                               */
// int     ndof ;   /* maximum number of nodal degrees of freedom           */
// int    *id   ;   /* equation number of free dof [numnp * ndof]           */
//		  /* ( >0 eqn. number, =0 restrained)                       */
// double *dn   ;   /*  vector of nodal displacements [neq]                 */

/*****************************************************************************/


#endif
