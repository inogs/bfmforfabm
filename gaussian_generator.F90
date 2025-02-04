        MODULE gaussian_generator
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     QUESTA FUNZIONE GENERA NUMERI GAUSSIANI
!     Definita a partire dall'articolo di Box-Muller 1958
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CONTAINS

      real function W(k,l)
          USE random_generator
          implicit none
          integer(4) :: k,l
          real       :: dummy
          real       :: r1,r2,pi,xx1,xx2
          r1 = ran2(k)
          r2 = ran2(l)
          pi = 3.141593
          dummy=1.0-r1
          xx1 = -2.0*log(dummy)
          xx2 = sqrt(xx1)*cos(2.0*pi*r2)
          W = med+var*xx2
          return
          end

         END MODULE
