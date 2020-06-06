/**********************[             GAUSS.h            ]********************/

#ifndef GAUSS_H
#define GAUSS_H

/******************************[ FUNCTION ]***********************************
 #                                                                           #
 #  NAME     : decomp ()                                                     #
 #                                                                           #
 #  TYPE     : void                                                          #
 #                                                                           #
 #  C LIBRARY: GAUSS.C                                                       #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #  DESCRIPTION : Gauss factorization for skyline storage saving diagonal    #
 #                pivots in array dg[].                                      #
 #                Execution breaking when a zero diagonal pivot appears.     #
 #                                                                           #
 #                (translation to C language of first part  of subroutine    #
 #                 COLSOL from  K.J.Bathe, "Finite Element Procedures in     #
 #                 Engineering Analysis", Prentice Hall, 1982.)              #
 #                                                                           #
 #****************************************************************************/


  void decomp (double *sg, int *maxa, int neq);

/*____________________________________________________________________________

                                  PARAMETERS
____________________________________________________________________________*/

//double *sg;
//int    *maxa;
//int    neq;

/****************************************************************************/





/******************************[ FUNCTION ]***********************************
 #                                                                           #
 #  NAME     : redbak ()                                                     #
 #                                                                           #
 #  TYPE     : void                                                          #
 #                                                                           #
 #  C LIBRARY: GAUSS.C                                                       #
 #                                                                           #
 *****************************************************************************
 #                                                                           #
 #  DESCRIPTION : Gauss backsubstitution for skyline storage.                #
 #                                                                           #
 #                (translation to C language of second part of subroutine    #
 #                 COLSOL from  K.J.Bathe, "Finite Element Procedures in     #
 #                 Engineering Analysis", Prentice Hall, 1982.)              #
 #                                                                           #
 #****************************************************************************/



 void redbak (double *sg, double *v, int *maxa, int neq);

// double *sg, *v;
// int *maxa, neq;

/*****************************************************************************/


#endif
