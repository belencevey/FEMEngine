
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


       void decomp (double *sg, int *maxa, int neq)


/*____________________________________________________________________________

                        PARAMETERS DECLARATION
____________________________________________________________________________*/


//double *sg;
//int    *maxa;
//int    neq;

{

/*____________________________________________________________________________

                       LOCAL VARIABLES DECLARATION
____________________________________________________________________________*/


 double  b, c;

 int   n, kn, kl, ku, kh, k, nnp,
       ic, klt, j, ki, nd, kk, l;

/*____________________________________________________________________________

             ZEROING COUNTER OF NEGATIVE PIVOTS
____________________________________________________________________________*/

                        nnp = 0;

/*____________________________________________________________________________

                 GAUSS FACTORIZATION
____________________________________________________________________________*/

   for (n = 1; n <= neq; n++)
      {
        kn = maxa[n];
        kl = kn + 1;
        ku = maxa [n + 1] - 1;
        kh = ku - kl;
        if (kh > 0)
          {
            k = n - kh;
            ic = 0;
            klt = ku;
            for ( j=1 ; j<=kh ; j++)
               {
                 ic++;
                 klt--;
                 ki = maxa[k];
                 nd = maxa[k + 1] - ki - 1;
                 if (nd > 0)
                   {
                     kk = (ic < nd) ? ic : nd;
                     c = 0.;
                     for ( l=1 ; l<=kk ; l++)  c += sg[ki + l] * sg[klt + l];
                     sg[klt] -= c;
                   }
                 k++;
               }
          }

        if (kh >= 0)
          {
            k = n;
            b = 0.0;
            for ( kk=kl ; kk<=ku ; kk++)
               {
                 k--;
                 ki = maxa[k];
                 c = sg[kk] / sg[ki];
                 b += c * sg[kk];
                 sg[kk] = c;
               }
            sg[kn] -= b;
          }

        if (sg[kn] < 0) nnp++;


      }


   }/* end of decomp() */






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



 void redbak (double *sg, double *v, int *maxa, int neq)

// double *sg, *v;
// int *maxa, neq;

 {
  double c;
  int    n, kl, k, ku, kk, l;

  for ( n=1 ; n<=neq ; n++)
     {
        kl = maxa[n] + 1;
        ku = maxa[n + 1] - 1;
        if (ku >= kl)
          {
            k = n;
            c = 0.;
            for ( kk=kl ; kk<=ku ; kk++)
               {
                  k--;
                  c += sg[kk] * v[k];
               }
            v[n] -= c;
          }
     }


    /* backsubstitution */

    for ( n=1 ; n<=neq ; n++)  v[n] /= sg[maxa[n]];
    if (neq > 1)
      {
        n = neq;
        for ( l=2 ; l<=neq ; l++)
           {
              kl = maxa[n] + 1;
              ku = maxa[n + 1] - 1;
              if ( ku >= kl)
                {
                   k = n;
                   for ( kk=kl ; kk<=ku ; kk++)
                      {
                         k--;
                         v[k] -= sg[kk] * v[n];
                      }
                 }
              n--;
           }
      }

  } /* end of redbak() */

