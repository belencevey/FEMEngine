/**********************[          Inputest.h            ]********************/
#include"Inputest.h"

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


	void inputest (int *ESTIMATOR, double *Uexc, int ecode)

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

/*____________________________________________________________________________

                           BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				             LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 char  name[60] ;  /* input filename                                        */

/*____________________________________________________________________________

         	          GETS ENERGY OF EXACT SOLUTION IF GIVEN
____________________________________________________________________________*/

        printf("\n\n  SQUARE ENERGY NORM OF EXACT SOLUTION: ");
        gets(name);
        *Uexc = atof(name);
        if (*Uexc<0)
          {
            printf("\n\n  MUST BE ZERO OR A POSITIVE NUMBER...");
            printf("\n\n  PROGRAM TERMINATED");
            exit(1);
          }

/*____________________________________________________________________________

         	          GETS ESTIMATOR TYPE
____________________________________________________________________________*/

   switch(ecode)
    {
      case 1:
      case 2:
	  case 7:
	  case 8:
        printf  ("\n\n  STRESS SMOOTHER TYPE 1 : NODAL AVERAGING ");
        printf  ("\n  STRESS SMOOTHER TYPE 2 : SIDE FLUX AVERAGING (DISCONTINUOUS)");
        printf  (
          "\n  STRESS SMOOTHER TYPE 3 : ZZ VERTEX RECOVERY CENTROIDAL SAMPLING");
        printf  (
          "\n  STRESS SMOOTHER TYPE 4 : ZZ VERTEX RECOVERY MID-SIDE SAMPLING");
        printf  (
          "\n  STRESS SMOOTHER TYPE 5 : SIDE AVERAGING (DISCONTINUOUS)");
        printf("\n\n  STRESS SMOOTHER TYPE : ");
        gets(name);
        *ESTIMATOR = atoi(name);
        if (*ESTIMATOR>5 || *ESTIMATOR<1)
          {
            printf("\n\n  ESTIMATOR OUT OF RANGE... ");
            printf("\n\n  PROGRAM TERMINATED");
            exit(1);
          }
       break;
       case 3:
         *ESTIMATOR=2;
       break;

       case 4:
         *ESTIMATOR=2;
       break;

       case 5:
         *ESTIMATOR=2;
       break;
     }

/*____________________________________________________________________________

		                       	END OF FUNCTION
____________________________________________________________________________*/

}/* end of inputest () */



