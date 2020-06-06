/**********************[         Inputpar.cpp           ]********************/
#include"Inputpar.h"


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


	void inputpar (int *SOLVER, int ecode)

/*____________________________________________________________________________

			                 PARAMETERS DECLARATION
____________________________________________________________________________*/

// int  *SOLVER ;  /* type of solver                                          */
//                 /*    SOLVER = 1   Gauss Eimination                        */
//                 /*    SOLVER = 2   Orthogonal Decomposition                */
//
// int    ecode ;  /* element code number.                                    */

/*____________________________________________________________________________

                           BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

/*____________________________________________________________________________

				             LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/

 char  SolverName[60] ;  /* solver type name                                */
 int  type;

/*____________________________________________________________________________

         	          GETS SOLVER TYPE
____________________________________________________________________________*/

   switch(ecode)
    {
      case 1:
      case 2:
	  case 7:
	  case 8:
        *SOLVER=1;
        break;

      case 3:
      case 4:
        printf("\n\n  SOLVER TYPE 1 : GAUSS ELIMINATION        ");
        printf  ("\n  SOLVER TYPE 2 : ORTHOGONAL DECOMPOSITION ");

        printf("\n\n  SOLVER TYPE : ");
        gets(SolverName);
        *SOLVER = atoi(SolverName);
        if (*SOLVER>2 || *SOLVER<1)
          {
            printf("\n\n  SOLVER OUT OF RANGE... ");
            printf("\n\n  PROGRAM TERMINATED");
            exit(1);
          }
        break;
      case 5:
        printf("\n\n  SOLVER TYPE 1 : GAUSS ELIMINATION        ");
        printf  ("\n  SOLVER TYPE 2 : ORTHOGONAL DECOMPOSITION ");

        printf("\n\n  SOLVER TYPE : ");
        gets(SolverName);
        *SOLVER = atoi(SolverName);
        if (*SOLVER>2 || *SOLVER<1)
          {
            printf("\n\n  SOLVER OUT OF RANGE... ");
            printf("\n\n  PROGRAM TERMINATED");
            exit(1);
          }
        break;

    }
/*____________________________________________________________________________

		                       	END OF FUNCTION
____________________________________________________________________________*/

}/* end of inputpar () */




