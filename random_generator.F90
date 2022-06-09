       MODULE random_generator

       integer(4), parameter :: ntab=32
       real                  :: med=0.0, var=1.0 
       integer(4)            :: idum2=123456789,iv(ntab)=0,iy=0

!----------------------------------------------------------------------

      CONTAINS

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     QUESTA FUNZIONE GENERA NUMERI GAUSSIANI
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     QUESTA FUNZIONE GENERA NUMERI RANDOM (vedi NUMERICAL RECIPES 
!     di W.H. et al. a pag 272)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
          real function ran2(idum)
          implicit none
          integer(4)            :: idum,j,k
          integer(4)            :: idum2
          integer(4), parameter :: im1=2147483563
          integer(4), parameter :: im2=2147483399
          real, parameter       :: am=(1.0/im1)
          integer(4), parameter :: imm1=im1-1
          integer(4), parameter :: ia1=40014
          integer(4), parameter :: ia2=40692
          integer(4), parameter :: iq1=53668
          integer(4), parameter :: iq2=52774
          integer(4), parameter :: ir1=12211
          integer(4), parameter :: ir2=3791
          integer(4), parameter :: ndiv=1+(imm1/ntab)
          real, parameter       :: eps=1.2e-7
          real, parameter       :: rnmx=1.0-eps

          if (idum.le.0) then

             idum=max(-idum,1)
             idum2=idum

             do j=ntab+8,1,-1

                k=idum/iq1

                idum=ia1*(idum-k*iq1)-k*ir1

                if (idum.lt.0) idum=idum+im1

                if (j.le.ntab) iv(j)=idum

             enddo

             iy=iv(1)

          endif

          k=idum/iq1
          idum=ia1*(idum-k*iq1)-k*ir1

          if (idum.lt.0) idum=idum+im1

          k=idum2/iq2
          idum2=ia2*(idum2-k*iq2)-k*ir2

          if (idum2.lt.0) idum2=idum2+im2

          j=1+iy/ndiv
          iy=iv(j)-idum2
          iv(j)=idum

          if (iy.lt.1) iy=iy+imm1

          ran2=min((am*iy),rnmx)

          return

          end

      END MODULE
