/**********************[           Init_lfe.cpp            ]********************/
#include"Init_lfe.h"


/**********************[ USER DEFINED LIBRARY FUNCTIONS ]*********************
 #                                                                           #
 #  LIBRARY NAME: INIT_LFE.C                                                 #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #                                                                           #
 #                Routines for display presentation                          #
 #                screen for program LSAFE.                                  #
 #                                                                           #
 #                                                                           #
 #----[ USER DEFINED FUNCTIONS IN THIS LIBRARY ]-----------------------------#
 #                                                                           #
 #  void  initp()      : displays presentation screen.                       #
 #                                                                           #
 #****************************************************************************/




/*******************************************[ USER LIBRARY : INIT_LFE.C ]*****
 #                                                                           #
 #  FUNCTION :  initp ()                                                     #
 #                                                                           #
 #  TYPE     :  void                                                         #
 #                                                                           #
 #                                                                           #
 #----[ DESCRIPTION ]--------------------------------------------------------#
 #                                                                           #
 #      Displays presentation screen for program LSAFE.                      #
 #      Uses extended  ASCII codes for some characters                       #
 #      displayed, so that special modifications could                       #
 #      be needed for operating systems other than DOS.                      #
 #                                                                           #
 #                                                                           #
 #----[ STANDARD LIBRARY FUNCTIONS USED ]------------------------------------#
 #                                                                           #
 #                                                                           #
 #  STDIO.H  :                                                               #
 #    int     printf():   formatted output to stdout                         #
 #                                                                           #
 #                                                                           #
 #****************************************************************************/

  void initp ()

/*____________________________________________________________________________

								 BEGINNING OF FUNCTION
____________________________________________________________________________*/

{

 printf("\n\n\n");
 printf("\n *============================================================");
 printf("=============* ");
 printf("\n *                                                           ");
 printf("              * ");
 printf("\n *                  MIXED BIDIMENSIONAL LINEAR STRUCTURAL    ");
 printf("              * ");
 printf("\n *                       ANALYSIS by FINITE ELEMENTS         ");
 printf("              * ");
 printf("\n *                                                           ");
 printf("              * ");
 printf("\n *                              S A F E 2 D - A              ");
 printf("              * ");
 printf("\n *                                                           ");
 printf("              * ");
 printf("\n *                                                           ");
 printf("              * ");
 printf("\n *           Dev. by  C.JOUGLARD - JUL/2002/FIUBA/ARGENTINA  ");
 printf("              * ");
 printf("\n *                                                           ");
 printf("              * ");
 printf("\n *===========================================================");
 printf("==============* ");

 printf("\n\n");


}/* end of initp() */


