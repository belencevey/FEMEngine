

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



#include <stdio.h>
#include <stdlib.h>
#include <string.h>



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


	void iofile (char *filename, FILE **fpi, char *forname, FILE **fpo, FILE **fpiter)

/*____________________________________________________________________________

			                 PARAMETERS DECLARATION
____________________________________________________________________________*/

// char   *filename   ;  /* complete input file name                          */
// FILE  **fpi        ;  /* pointer input  file                               */
// char   *forname    ;  /* complete output file root name                    */
// FILE  **fpo        ;  /* pointer results file                              */
// FILE  **fpiter     ;  /* pointer iteration steps file                      */

/*____________________________________________________________________________

                           BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				             LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 char *ext_name     ;  /* extension of input file                           */
// char  filename[60] ;  /* input filename                                    */
 int   length       ;  /* length of root name string                        */

/*____________________________________________________________________________

			             PROTOTYPES FUNCTIONS DECLARATION
____________________________________________________________________________*/


 FILE *fopenvf (char *finame, char *mode);

/*____________________________________________________________________________

	                    VERIFY EXTENSION OCURRENCE
                        FILE ROOT NAME
____________________________________________________________________________*/

 ext_name = strrchr(filename, '.');

 if ((strlen(ext_name)>0) && (!strcmp(ext_name, ".2da") || !strcmp(ext_name, ".2DA")))
   {
     length = strlen(filename)-strlen(ext_name);
     strncpy (forname, filename, length);
     forname[length] = 0;
   }
 else
 {
		 printf("\n\n file extension must be .2da \n");
		 exit(2);
 }
   // strcpy (forname, filename);

/*____________________________________________________________________________

		              OPENING OF INPUT FILE
____________________________________________________________________________*/

			 *fpi = fopenvf (filename, "r");

/*____________________________________________________________________________

		              OPENING OF RESULTS FILE
____________________________________________________________________________*/

			strcpy (filename, forname);
			strcat (filename, ".out");

	 *fpo = fopenvf (filename, "w");

/*____________________________________________________________________________

		              OPENING OF ITERATION FILE
____________________________________________________________________________*/

			strcpy (filename, forname);
			strcat (filename, ".itr");

	 *fpiter = fopenvf (filename, "w");

/*____________________________________________________________________________

			END OF FUNCTION
____________________________________________________________________________*/

}/* end of iofile () */











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


  FILE *fopenvf (char *finame, char *mode)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// char   *finame      ;  /* name of the file to be open                      */
// char   *mode        ;  /* type of opening                                  */


/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 FILE   *fpi         ;  /* pointer to the newly open file                   */

/*____________________________________________________________________________

			  OPEN FILE WITH VERIFICATION
____________________________________________________________________________*/


		fpi = fopen (finame, mode);
		if (fpi == NULL)
	     {
		    printf("\n\n Can't open file %s \n",finame);
		    exit(2);
        }

/*____________________________________________________________________________

		  RETURNS POINTER TO NEWLY OPEN FILE
____________________________________________________________________________*/

		return fpi;

/*____________________________________________________________________________

			END OF FUNCTION
____________________________________________________________________________*/

}/* end of fopenvf() */









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


  void prtdis (FILE *fpo, int numnp, int ndof, int *id, double *dn)

/*____________________________________________________________________________

			PARAMETERS DECLARATION
____________________________________________________________________________*/

// FILE   *fpo  ;   /* pointer to output file                                 */
// int     numnp;   /* number of nodal points                                 */
// int     ndof ;   /* maximum number of nodal degrees of freedom             */
// int    *id   ;   /* equation number of free dof [numnp * ndof]             */
//		  /* ( >0 eqn. number, =0 restrained)                       */
// double *dn   ;   /*  vector of nodal displacements [neq]                   */

/*____________________________________________________________________________

			BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 double  disp ;   /* value of nodal displacement                            */
 int     j1   ;   /* auxiliar pointer for vector id[]                       */
 int     i    ;   /* loop counter on nodal points                           */
 int     j    ;   /* loop counter on nodal degrees of freedom               */

/*____________________________________________________________________________

			LOOP ON NODAL POINTS
____________________________________________________________________________*/

  for ( i=1 ; i<=numnp ; i++)
	  {
	 fprintf (fpo,"%5d  ", i);

/*____________________________________________________________________________

				PRINT NODAL DISPLACEMENTS
____________________________________________________________________________*/

	 j1 = ndof * (i -1);

	 for ( j=1 ; j<=ndof ; j++)
		 {
			 disp = (id[j1 + j] > 0) ? dn[id[j1 + j]] : 0.0;
			 fprintf (fpo,"%+12.5e   ", disp);
		 }

	 fprintf (fpo,"\n");

/*____________________________________________________________________________

			  END OF LOOP ON NODAL POINTS
____________________________________________________________________________*/

	  }
/*____________________________________________________________________________

			END OF FUNCTION
____________________________________________________________________________*/

}/* end of prtdis() */
