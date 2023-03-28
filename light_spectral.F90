#include "fabm_driver.h"

module ogs_bfm_light_spectral

   use fabm_types
   use ogs_bfm_shared
   use adj_3stream, only: solve_direct

   implicit none

   private

   type,extends(type_base_model),public :: type_ogs_bfm_light_spectral
      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id)   :: id_par_P1,id_par_P2,id_par_P3,id_par_P4,id_par_P5,id_par_P6,id_par_P7,id_par_P8,id_par_P9
      type (type_diagnostic_variable_id)   :: id_PAR_tot
      type (type_diagnostic_variable_id)   :: id_P1chl,id_P2chl,id_P3chl,id_P4chl,id_P5chl,id_P6chl,id_P7chl,id_P8chl,id_P9chl
      type (type_diagnostic_variable_id)   :: id_P1c  ,id_P2c  ,id_P3c  ,id_P4c  ,id_P5c  ,id_P6c  ,id_P7c  ,id_P8c,  id_P9c      
!      type (type_diagnostic_variable_id)   :: id_anap450, id_aph450
!      type (type_diagnostic_variable_id)   :: id_acdom250, id_acdom325, id_acdom400, id_acdom425, id_acdom450
!      type (type_diagnostic_variable_id)   :: id_Scdom350_500, id_Scdom250_325       
!      type (type_diagnostic_variable_id)   :: id_bbp450, id_bbp550, id_bbp700

!      type (type_horizontal_diagnostic_variable_id) :: id_Rrs400, id_Rrs425, id_Rrs450, id_Rrs475
!      type (type_horizontal_diagnostic_variable_id) :: id_Rrs500, id_Rrs525, id_Rrs550, id_Rrs575, id_Rrs675
!      type (type_horizontal_diagnostic_variable_id) :: id_kd375, id_kd400, id_kd425, id_kd475, id_kd500
      
      type (type_dependency_id)            :: id_dz
      type (type_dependency_id)            :: id_aP1c, id_aP2c, id_aP3c, id_aP4c, id_aP5c, id_aP6c, id_aP7c, id_aP8c, id_aP9c
      type (type_dependency_id)            :: id_aP1chl,id_aP2chl,id_aP3chl,id_aP4chl,id_aP5chl,id_aP6chl,id_aP7chl,id_aP8chl,id_aP9chl
      type (type_state_variable_id)        :: id_R6c, id_X1c, id_X2c, id_X3c
      type (type_horizontal_dependency_id) :: id_zenithA

! BLOCK 1 python generated code see AUX_SCRIPTS/python_light_spectral.py
      type (type_horizontal_dependency_id) ::  id_Ed_0_0250, id_Ed_0_0325, id_Ed_0_0350, id_Ed_0_0375, id_Ed_0_0400
      type (type_horizontal_dependency_id) ::  id_Ed_0_0425, id_Ed_0_0450, id_Ed_0_0475, id_Ed_0_0500, id_Ed_0_0525
      type (type_horizontal_dependency_id) ::  id_Ed_0_0550, id_Ed_0_0575, id_Ed_0_0600, id_Ed_0_0625, id_Ed_0_0650
      type (type_horizontal_dependency_id) ::  id_Ed_0_0675, id_Ed_0_0700, id_Ed_0_0725, id_Ed_0_0775, id_Ed_0_0850
      type (type_horizontal_dependency_id) ::  id_Ed_0_0950, id_Ed_0_1050, id_Ed_0_1150, id_Ed_0_1250, id_Ed_0_1350
      type (type_horizontal_dependency_id) ::  id_Ed_0_1450, id_Ed_0_1550, id_Ed_0_1650, id_Ed_0_1750, id_Ed_0_1900
      type (type_horizontal_dependency_id) ::  id_Ed_0_2200, id_Ed_0_2900, id_Ed_0_3700
      type (type_horizontal_dependency_id) ::  id_Es_0_0250, id_Es_0_0325, id_Es_0_0350, id_Es_0_0375, id_Es_0_0400
      type (type_horizontal_dependency_id) ::  id_Es_0_0425, id_Es_0_0450, id_Es_0_0475, id_Es_0_0500, id_Es_0_0525
      type (type_horizontal_dependency_id) ::  id_Es_0_0550, id_Es_0_0575, id_Es_0_0600, id_Es_0_0625, id_Es_0_0650
      type (type_horizontal_dependency_id) ::  id_Es_0_0675, id_Es_0_0700, id_Es_0_0725, id_Es_0_0775, id_Es_0_0850
      type (type_horizontal_dependency_id) ::  id_Es_0_0950, id_Es_0_1050, id_Es_0_1150, id_Es_0_1250, id_Es_0_1350
      type (type_horizontal_dependency_id) ::  id_Es_0_1450, id_Es_0_1550, id_Es_0_1650, id_Es_0_1750, id_Es_0_1900
      type (type_horizontal_dependency_id) ::  id_Es_0_2200, id_Es_0_2900, id_Es_0_3700

! END BLOCK 1 python generated  code

      ! Parameters
      integer  :: nlt,npft
      real(rk) :: rd, rs, ru, vs, vu
      real(rk) :: SdomX1, X1coeff, SdomX2, X2coeff, SdomX3, X3coeff, lambda_aCDOM, Xmincoeff
      real(rk) :: Sapar, lambda_aPart, aparcoeff
      real(rk) :: Sbpar, lambda_bPart, bparcoeff, bb_to_b
      logical :: compute_acdom
      logical :: compute_anap
      real(rk) :: p_epsP1,p_epsP2,p_epsP3,p_epsP4,p_epsP5,p_epsP6,p_epsP7,p_epsP8,p_epsP9
      real(rk) :: p_bpsP1,p_bpsP2,p_bpsP3,p_bpsP4,p_bpsP5,p_bpsP6,p_bpsP7,p_bpsP8,p_bpsP9
      real(rk) :: p_bbrP1,p_bbrP2,p_bbrP3,p_bbrP4,p_bbrP5,p_bbrP6,p_bbrP7,p_bbrP8,p_bbrP9
      logical :: compute_aph, compute_bph, compute_bbc
      
   contains
!     Model procedures
      procedure :: initialize
      procedure :: do_column
   end type type_ogs_bfm_light_spectral

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ogs_bfm_light_spectral),intent(inout),target :: self
      integer,                     intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      real(rk) :: hc, hcoavo, rlamm, rlamm1, rlamm2, nl
      real(rk) :: n, dummy_p, cu_area, aph_mean, bph_mean

!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%nlt,    'nlt',  '-',   'number of wavelenghts', default=-1)
      call self%get_parameter(self%npft,   'npft', '-',   'number of PFT', default=-1)
      call self%get_parameter(self%rd,     'rd',   '-',   ' ', default=1.0_rk)
      call self%get_parameter(self%rs,     'rs',   '-',   ' ', default=1.5_rk)
      call self%get_parameter(self%ru,     'ru',   '-',   ' ', default=3.0_rk)
      call self%get_parameter(self%vs,     'vs',   '',    'avg cosine diffuse down', default=0.83_rk)
      call self%get_parameter(self%vu,     'vu',   '-',   'avg cosine diffuse up', default=0.4_rk)
      call self%get_parameter(self%SdomX1,        'SdomX1',        'nm-1',     'slope for aCDOM [X1c] wavelength dependence')
      call self%get_parameter(self%X1coeff,       'X1coeff',       'm2 mgC-1', 'specific absorption of X1c at lambda_aCDOM ')
      call self%get_parameter(self%SdomX2,        'SdomX2',        'nm-1',     'slope for aCDOM [X2c] wavelength dependence')
      call self%get_parameter(self%X2coeff,       'X2coeff',       'm2 mgC-1', 'specific absorption of X2c at lambda_aCDOM ')
      call self%get_parameter(self%SdomX3,        'SdomX3',        'nm-1',     'slope for aCDOM [X3c] wavelength dependence')
      call self%get_parameter(self%X3coeff,       'X3coeff',       'm2 mgC-1', 'specific absorption of X3c at lambda_aCDOM ')
      call self%get_parameter(self%lambda_aCDOM,  'lambda_aCDOM',  'nm',       'wavelength where reference aCDOM is given')
      call self%get_parameter(self%Xmincoeff,     'Xmincoeff',     'm-1',      'minimum aCDOM at 450nm')      
      call self%get_parameter(self%Sapar,         'Sapar',         'nm-1',     'slope parameter for aNAP wavelength dependence')
      call self%get_parameter(self%lambda_aPart,  'lambda_aPart',  'nm',       'wavelength where reference aNAP is given')
      call self%get_parameter(self%aparcoeff,     'aparcoeff',     'm2 mgC-1', 'specific absorption at lambda_aPart ')
      call self%get_parameter(self%Sbpar,         'Sbpar',         '-',        'exponent for bNAP wavelength dependence')
      call self%get_parameter(self%lambda_bPart,  'lambda_bPart',  'nm',       'wavelength where reference bNAP is given')
      call self%get_parameter(self%bparcoeff,     'bparcoeff',     'm2 mgC-1', 'specific scatter at lambda_bPart')
      call self%get_parameter(self%bb_to_b,       'bb_to_b',       '-',        'backscatter to total scatter ratio') 
      call self%get_parameter(self%compute_acdom, 'compute_acdom', '[T or F]', 'logical flag to compute acdom') 
      call self%get_parameter(self%compute_anap,  'compute_anap',  '[T or F]', 'logical flag to compute anap')
      call self%get_parameter(self%p_epsP1,       'p_epsP1',       'm2 mgChl-1',  'mean absorption coefficient from 400-700nm for P1', default=0.03_rk)
      call self%get_parameter(self%p_epsP2,       'p_epsP2',       'm2 mgChl-1',  'mean absorption coefficient from 400-700nm for P2', default=0.03_rk)
      call self%get_parameter(self%p_epsP3,       'p_epsP3',       'm2 mgChl-1',  'mean absorption coefficient from 400-700nm for P3', default=0.03_rk)
      call self%get_parameter(self%p_epsP4,       'p_epsP4',       'm2 mgChl-1',  'mean absorption coefficient from 400-700nm for P4', default=0.03_rk)
      call self%get_parameter(self%p_epsP5,       'p_epsP5',       'm2 mgChl-1',  'mean absorption coefficient from 400-700nm for P5', default=0.03_rk)
      call self%get_parameter(self%p_epsP6,       'p_epsP6',       'm2 mgChl-1',  'mean absorption coefficient from 400-700nm for P6', default=0.03_rk)
      call self%get_parameter(self%p_epsP7,       'p_epsP7',       'm2 mgChl-1',  'mean absorption coefficient from 400-700nm for P7', default=0.03_rk)
      call self%get_parameter(self%p_epsP8,       'p_epsP8',       'm2 mgChl-1',  'mean absorption coefficient from 400-700nm for P8', default=0.03_rk)
      call self%get_parameter(self%p_epsP9,       'p_epsP9',       'm2 mgChl-1',  'mean absorption coefficient from 400-700nm for P9', default=0.03_rk)
      call self%get_parameter(self%compute_aph,   'compute_aph',   '[T or F]',    'logical flag to scale aph')
      call self%get_parameter(self%p_bpsP1,       'p_bpsP1',       'm2 mgC-1',  'mean scattering coefficient from 400-700nm for P1')
      call self%get_parameter(self%p_bpsP2,       'p_bpsP2',       'm2 mgC-1',  'mean scattering coefficient from 400-700nm for P2')
      call self%get_parameter(self%p_bpsP3,       'p_bpsP3',       'm2 mgC-1',  'mean scattering coefficient from 400-700nm for P3')
      call self%get_parameter(self%p_bpsP4,       'p_bpsP4',       'm2 mgC-1',  'mean scattering coefficient from 400-700nm for P4')
      call self%get_parameter(self%p_bpsP5,       'p_bpsP5',       'm2 mgC-1',  'mean scattering coefficient from 400-700nm for P5')
      call self%get_parameter(self%p_bpsP6,       'p_bpsP6',       'm2 mgC-1',  'mean scattering coefficient from 400-700nm for P6')
      call self%get_parameter(self%p_bpsP7,       'p_bpsP7',       'm2 mgC-1',  'mean scattering coefficient from 400-700nm for P7')
      call self%get_parameter(self%p_bpsP8,       'p_bpsP8',       'm2 mgC-1',  'mean scattering coefficient from 400-700nm for P8')
      call self%get_parameter(self%p_bpsP9,       'p_bpsP9',       'm2 mgC-1',  'mean scattering coefficient from 400-700nm for P9')      
      call self%get_parameter(self%compute_bph,   'compute_bph',   '[T or F]',    'logical flag to scale bph')
      call self%get_parameter(self%p_bbrP1,       'p_bbrP1',       '-',           'backscattering to total scattering ratio for P1')
      call self%get_parameter(self%p_bbrP2,       'p_bbrP2',       '-',           'backscattering to total scattering ratio for P2')
      call self%get_parameter(self%p_bbrP3,       'p_bbrP3',       '-',           'backscattering to total scattering ratio for P3')
      call self%get_parameter(self%p_bbrP4,       'p_bbrP4',       '-',           'backscattering to total scattering ratio for P4')
      call self%get_parameter(self%p_bbrP5,       'p_bbrP5',       '-',           'backscattering to total scattering ratio for P5')
      call self%get_parameter(self%p_bbrP6,       'p_bbrP6',       '-',           'backscattering to total scattering ratio for P6')
      call self%get_parameter(self%p_bbrP7,       'p_bbrP7',       '-',           'backscattering to total scattering ratio for P7')
      call self%get_parameter(self%p_bbrP8,       'p_bbrP8',       '-',           'backscattering to total scattering ratio for P8')
      call self%get_parameter(self%p_bbrP9,       'p_bbrP9',       '-',           'backscattering to total scattering ratio for P9')
      call self%get_parameter(self%compute_bbc,   'compute_bbc',   '[T or F]',    'logical flag to compute bbc from bc*bbr')
      
      if (self%nlt>0) then
          allocate(lam(self%nlt));             lam(:)=huge(lam(1))
          allocate(lam1(self%nlt));            lam1(:)=huge(lam1(1))
          allocate(lam2(self%nlt));            lam2(:)=huge(lam2(1))
          allocate(aw(self%nlt));              aw(:)=huge(aw(1))
          allocate(bw(self%nlt));              bw(:)=huge(bw(1))
          allocate(bbw(self%nlt));             bbw(:)=huge(bbw(1))
          allocate(ac(self%npft,self%nlt));    ac(:,:)=huge(ac(1,1))
          allocate(ac_ps(self%npft,self%nlt)); ac_ps(:,:)=huge(ac_ps(1,1))
          allocate(bc(self%npft,self%nlt));    bc(:,:)=huge(bc(1,1))
          allocate(bbc(self%npft,self%nlt));   bbc(:,:)=huge(bbc(1,1))
          allocate(apoc(self%nlt));            apoc(:)=huge(apoc(1))
          allocate(bpoc(self%nlt));            bpoc(:)=huge(bpoc(1))
          allocate(bbpoc(self%nlt));           bbpoc(:)=huge(bbpoc(1))
          allocate(acdom_min(self%nlt));       acdom_min(:)=huge(acdom_min(1))
          allocate(acdom(3,self%nlt));         acdom(:,:)=huge(acdom(1,1))
          allocate(Ed_0(self%nlt));            Ed_0(:)=huge(Ed_0(1))
          allocate(Es_0(self%nlt));            Es_0(:)=huge(Es_0(1))
          allocate(WtoQ(self%nlt));            WtoQ(:)=huge(WtoQ(1))
!          allocate(equis(7));                  equis(:)=huge(equis(1))
!          allocate(ies(7));                    ies(:)=huge(ies(1))

          
!         load the IOP for the biogeochemical variables considered
          call lidata(self%nlt,self%npft)

      endif 

     hc = 1.0D0/(h_planck*c_light)
     hcoavo = hc*oavo


      do nl = 1,self%nlt
       rlamm = real(lam(nl),8)*1.0E-9      !lambda in m
       WtoQ(nl) = rlamm*hcoavo*1000000.0D0 !Watts to micro mol quanta conversion
       acdom_min(nl)= 0.0_rk
      enddo

!   with parameter in namelist     
      do nl = 1,self%nlt
       rlamm1 = real(lam1(nl),8)
       rlamm2 = real(lam2(nl),8)
       acdom_min(nl) = self%Xmincoeff*(exp(-self%SdomX2*(rlamm2-self%lambda_aCDOM))-exp(-self%SdomX2*(rlamm1-self%lambda_aCDOM)))/(-self%SdomX2*(rlamm2-rlamm1))
      enddo

!      do nl = 1,self%nlt
!        write(*,*) real(lam(nl),8), acdom_min(nl)
!      enddo
      
      
 !   CDOM absorption coefficients
      if (self%compute_acdom) then     
      do nl = 1,self%nlt
 !      rlamm = real(lam(nl),8)
       rlamm1 = real(lam1(nl),8)
       rlamm2 = real(lam2(nl),8)
 !      acdom(nl) = self%cdomcoeff * exp(-self%Sdom*(rlamm-self%lambda_aCDOM))
       acdom(1,nl) = self%X1coeff*(exp(-self%SdomX1*(rlamm2-self%lambda_aCDOM))-exp(-self%SdomX1*(rlamm1-self%lambda_aCDOM)))/(-self%SdomX1*(rlamm2-rlamm1))
       acdom(2,nl) = self%X2coeff*(exp(-self%SdomX2*(rlamm2-self%lambda_aCDOM))-exp(-self%SdomX2*(rlamm1-self%lambda_aCDOM)))/(-self%SdomX2*(rlamm2-rlamm1))
       acdom(3,nl) = self%X3coeff*(exp(-self%SdomX3*(rlamm2-self%lambda_aCDOM))-exp(-self%SdomX3*(rlamm1-self%lambda_aCDOM)))/(-self%SdomX3*(rlamm2-rlamm1))
      enddo
      endif

 !   POC absorption/scatter coefficients
      if (self%compute_anap) then
      do nl = 1,self%nlt
       rlamm = real(lam(nl),8)
       rlamm1 = real(lam1(nl),8)
       rlamm2 = real(lam2(nl),8)
 !      apoc(nl) = self%aparcoeff * exp(-self%Sapar*(rlamm-self%lambda_aPart))
       apoc(nl) = self%aparcoeff*(exp(-self%Sapar*(rlamm2-self%lambda_aPart))-exp(-self%Sapar*(rlamm1-self%lambda_aPart)))/(-self%Sapar*(rlamm2-rlamm1))     
       bpoc(nl) = self%bparcoeff * ((self%lambda_bPart/rlamm)**self%Sbpar)
       bbpoc(nl) = self%bparcoeff * ((self%lambda_bPart/rlamm)**self%Sbpar) * self%bb_to_b
      enddo
      endif

!       write(*,*) "acdom", acdom(1,7), acdom(2,7), acdom(3,7)
!       write(*,*) "apoc", apoc(7)
!       write(*,*) "bpoc", bpoc(7)
!       write(*,*) "bbpoc", bbpoc(7)

      do nl = 1,self%nlt
        write(*,*) real(lam(nl),8), acdom(1,nl), acdom(2,nl), acdom(3,nl)
      enddo

      do nl = 1,self%nlt
        write(*,*) real(lam(nl),8), apoc(nl), bpoc(nl), bbpoc(nl)
      enddo
       
 !   PHYTO absorption coefficients
      if (self%compute_aph) then
         do n = 1,self%npft
            if (n == 1) dummy_p = self%p_epsP1
            if (n == 2) dummy_p = self%p_epsP2
            if (n == 3) dummy_p = self%p_epsP3
            if (n == 4) dummy_p = self%p_epsP4
            if (n == 5) dummy_p = self%p_epsP5
            if (n == 6) dummy_p = self%p_epsP6
            if (n == 7) dummy_p = self%p_epsP7
            if (n == 8) dummy_p = self%p_epsP8
            if (n == 9) dummy_p = self%p_epsP9            
            ! find mean
            cu_area = 0.0_rk
            do nl = 5,17   ! indexes for 400-700nm
!                rlamm = real(lam(nl),8)
                rlamm1 = real(lam1(nl),8)
                rlamm2 = real(lam2(nl),8)
                cu_area = cu_area + ((rlamm2-rlamm1) * ac(n,nl))
            enddo
            aph_mean = cu_area / (real(lam2(17),8)-real(lam1(5),8)) ! total band width 400-700nm
            do nl = 1,self%nlt
                ac(n,nl) = ac(n,nl) * (dummy_p/aph_mean)
            enddo 
         enddo
      endif
       
 !   PHYTO scattering coefficients
      if (self%compute_bph) then
         do n = 1,self%npft
            if (n == 1) dummy_p = self%p_bpsP1
            if (n == 2) dummy_p = self%p_bpsP2
            if (n == 3) dummy_p = self%p_bpsP3
            if (n == 4) dummy_p = self%p_bpsP4
            if (n == 5) dummy_p = self%p_bpsP5
            if (n == 6) dummy_p = self%p_bpsP6
            if (n == 7) dummy_p = self%p_bpsP7
            if (n == 8) dummy_p = self%p_bpsP8
            if (n == 9) dummy_p = self%p_bpsP9
            ! find mean
            cu_area = 0.0_rk
            do nl = 5,17   ! indexes for 400-700nm
                rlamm1 = real(lam1(nl),8)
                rlamm2 = real(lam2(nl),8)
                cu_area = cu_area + ((rlamm2-rlamm1) * bc(n,nl))
            enddo
            bph_mean = cu_area / (real(lam2(17),8)-real(lam1(5),8)) ! total band width 400-700nm
            do nl = 1,self%nlt
                bc(n,nl) = bc(n,nl) * (dummy_p/bph_mean)
            enddo 
         enddo
      endif

 !   PHYTO backscattering coefficients       
      if (self%compute_bbc) then
       do nl = 1,19
          bbc(1,nl) = self%p_bbrP1
          bbc(2,nl) = self%p_bbrP2
          bbc(3,nl) = self%p_bbrP3
          bbc(4,nl) = self%p_bbrP4
          bbc(5,nl) = self%p_bbrP5
          bbc(6,nl) = self%p_bbrP6
          bbc(7,nl) = self%p_bbrP7
          bbc(8,nl) = self%p_bbrP8
          bbc(9,nl) = self%p_bbrP9
       enddo
      endif

!      do nl = 1,self%nlt
!        write(*,*) real(lam(nl),8), ac(1,nl), ac_ps(1,nl), bc(1,nl), bbc(1,nl)
!      enddo
!      do nl = 1,self%nlt
!        write(*,*) real(lam(nl),8), ac(2,nl), ac_ps(2,nl), bc(2,nl), bbc(2,nl)
!      enddo
!      do nl = 1,self%nlt
!        write(*,*) real(lam(nl),8), ac(3,nl), ac_ps(3,nl), bc(3,nl), bbc(3,nl)
!      enddo
!      do nl = 1,self%nlt
!        write(*,*) real(lam(nl),8), ac(4,nl), ac_ps(4,nl), bc(4,nl), bbc(4,nl)
!!      enddo
!      do nl = 1,self%nlt
!        write(*,*) real(lam(nl),8), ac(5,nl), ac_ps(5,nl), bc(5,nl), bbc(5,nl)
!      enddo
!      do nl = 1,self%nlt
!        write(*,*) real(lam(nl),8), ac(6,nl), ac_ps(6,nl), bc(6,nl), bbc(6,nl)
!      enddo
!      do nl = 1,self%nlt
!        write(*,*) real(lam(nl),8), ac(7,nl), ac_ps(7,nl), bc(7,nl), bbc(7,nl)
!      enddo
!      do nl = 1,self%nlt
!        write(*,*) real(lam(nl),8), ac(8,nl), ac_ps(8,nl), bc(8,nl), bbc(8,nl)
!      enddo
!      do nl = 1,self%nlt
!        write(*,*) real(lam(nl),8), ac(9,nl), ac_ps(9,nl), bc(9,nl), bbc(9,nl)
!      enddo
      
     ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_par_P1, 'PAR_P1', 'uE mgChl-1 d-1', 'PAR_diatoms', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_P2, 'PAR_P2', 'uE mgChl-1 d-1', 'PAR_flagellates', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_P3, 'PAR_P3', 'uE mgChl-1 d-1', 'PAR_picoeukaryotes', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_P4, 'PAR_P4', 'uE mgChl-1 d-1', 'PAR_dinoflagellates', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_P5, 'PAR_P5', 'uE mgChl-1 d-1', 'PAR_coccoloithophores', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_P6, 'PAR_P6', 'uE mgChl-1 d-1', 'PAR_prochlorococcus', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_P7, 'PAR_P7', 'uE mgChl-1 d-1', 'PAR_green1', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_P8, 'PAR_P8', 'uE mgChl-1 d-1', 'PAR_green2', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_P9, 'PAR_P9', 'uE mgChl-1 d-1', 'PAR_synechococcus', source=source_do_column)
      call self%register_diagnostic_variable(self%id_PAR_tot, 'PAR_tot',  'uE m-2 d-1 [400-700]','PAR_total', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_Scdom350_500, 'Scdom350_500', 'nm-1','visible spectral slope acdom', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_Scdom250_325, 'Scdom250_325', 'nm-1','UV spectral slope acdom', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_acdom250, 'acdom250', 'm-1', 'acdom in 250 nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_acdom325, 'acdom325', 'm-1', 'acdom in 325 nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_acdom400, 'acdom400', 'm-1', 'acdom in 400 nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_acdom425, 'acdom425', 'm-1', 'acdom in 425 nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_acdom450, 'acdom450', 'm-1', 'acdom in 450 nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_anap450, 'anap450',  'm-1', 'anap in 450 nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_aph450,  'aph450',   'm-1', 'aph in 450 nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_bbp450,  'bbp450',   'm-1', 'particle backscattering in 450 nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_bbp550,  'bbp550',   'm-1', 'particle backscattering in 550 nm band', source=source_do_column)      
!      call self%register_diagnostic_variable(self%id_bbp700,  'bbp700',   'm-1', 'particle backscattering in 700 nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_Rrs400,   'Rrs400',   '-',  'subsurface reflectance in 400nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_Rrs425,   'Rrs425',   '-',  'subsurface reflectance in 425nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_Rrs450,   'Rrs450',   '-',  'subsurface reflectance in 450nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_Rrs475,   'Rrs475',   '-',  'subsurface reflectance in 475nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_Rrs500,   'Rrs500',   '-',  'subsurface reflectance in 500nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_Rrs525,   'Rrs525',   '-',  'subsurface reflectance in 525nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_Rrs550,   'Rrs550',   '-',  'subsurface reflectance in 550nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_Rrs575,   'Rrs575',   '-',  'subsurface reflectance in 575nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_Rrs675,   'Rrs675',   '-',  'subsurface reflectance in 675nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_kd375,    'kd375',  'm-1',  'extinction coefficient in 375nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_kd400,    'kd400',  'm-1',  'extinction coefficient in 400nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_kd425,    'kd425',  'm-1',  'extinction coefficient in 425nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_kd475,    'kd475',  'm-1',  'extinction coefficient in 475nm band', source=source_do_column)
!      call self%register_diagnostic_variable(self%id_kd500,    'kd500',  'm-1',  'extinction coefficient in 500nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_P1chl, 'P1chl', 'mg m-3', 'chlorophyll_P1', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P2chl, 'P2chl', 'mg m-3', 'chlorophyll_P2', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P3chl, 'P3chl', 'mg m-3', 'chlorophyll_P3', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P4chl, 'P4chl', 'mg m-3', 'chlorophyll_P4', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P5chl, 'P5chl', 'mg m-3', 'chlorophyll_P5', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P6chl, 'P6chl', 'mg m-3', 'chlorophyll_P6', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P7chl, 'P7chl', 'mg m-3', 'chlorophyll_P7', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P8chl, 'P8chl', 'mg m-3', 'chlorophyll_P8', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P9chl, 'P9chl', 'mg m-3', 'chlorophyll_P9', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P1c, 'P1c', 'mg m-3', 'carbon_P1', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P2c, 'P2c', 'mg m-3', 'carbon_P2', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P3c, 'P3c', 'mg m-3', 'carbon_P3', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P4c, 'P4c', 'mg m-3', 'carbon_P4', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P5c, 'P5c', 'mg m-3', 'carbon_P5', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P6c, 'P6c', 'mg m-3', 'carbon_P6', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P7c, 'P7c', 'mg m-3', 'carbon_P7', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P8c, 'P8c', 'mg m-3', 'carbon_P8', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_P9c, 'P9c', 'mg m-3', 'carbon_P9', source=source_do_column)      
      
      ! Register dependencies on aggregated variables
      call self%register_dependency(self%id_aP1c, carbon_P1)
      call self%register_dependency(self%id_aP2c, carbon_P2)
      call self%register_dependency(self%id_aP3c, carbon_P3)
      call self%register_dependency(self%id_aP4c, carbon_P4)
      call self%register_dependency(self%id_aP5c, carbon_P5)
      call self%register_dependency(self%id_aP6c, carbon_P6)
      call self%register_dependency(self%id_aP7c, carbon_P7)
      call self%register_dependency(self%id_aP8c, carbon_P8)
      call self%register_dependency(self%id_aP9c, carbon_P9)

      call self%register_dependency(self%id_aP1chl,chlorophyll_P1)
      call self%register_dependency(self%id_aP2chl,chlorophyll_P2)
      call self%register_dependency(self%id_aP3chl,chlorophyll_P3)
      call self%register_dependency(self%id_aP4chl,chlorophyll_P4)
      call self%register_dependency(self%id_aP5chl,chlorophyll_P5)
      call self%register_dependency(self%id_aP6chl,chlorophyll_P6)
      call self%register_dependency(self%id_aP7chl,chlorophyll_P7)
      call self%register_dependency(self%id_aP8chl,chlorophyll_P8)
      call self%register_dependency(self%id_aP9chl,chlorophyll_P9)

      ! Register biogeochemical dependencies
      call self%register_state_dependency(self%id_R6c,'R6c','mg C/m^3', 'POC')
      call self%register_state_dependency(self%id_X1c,'X1c','mg C/m^3', 'labile CDOM')
      call self%register_state_dependency(self%id_X2c,'X2c','mg C/m^3', 'semi-labile CDOM')
      call self%register_state_dependency(self%id_X3c,'X3c','mg C/m^3', 'semi-refractory CDOM')

      ! Register environmental dependencies 
      call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
      call self%register_horizontal_dependency(self%id_zenithA, type_horizontal_standard_variable(name='zenith_angle'))

! BLOCK 2 python generate code see AUX_SCRIPTS/python_light_spectral.py
      call self%register_dependency(self%id_Ed_0_0250,type_surface_standard_variable(name='surf_direct_downward_irradiance_0250_nm'))
      call self%register_dependency(self%id_Ed_0_0325,type_surface_standard_variable(name='surf_direct_downward_irradiance_0325_nm'))
      call self%register_dependency(self%id_Ed_0_0350,type_surface_standard_variable(name='surf_direct_downward_irradiance_0350_nm'))
      call self%register_dependency(self%id_Ed_0_0375,type_surface_standard_variable(name='surf_direct_downward_irradiance_0375_nm'))
      call self%register_dependency(self%id_Ed_0_0400,type_surface_standard_variable(name='surf_direct_downward_irradiance_0400_nm'))
      call self%register_dependency(self%id_Ed_0_0425,type_surface_standard_variable(name='surf_direct_downward_irradiance_0425_nm'))
      call self%register_dependency(self%id_Ed_0_0450,type_surface_standard_variable(name='surf_direct_downward_irradiance_0450_nm'))
      call self%register_dependency(self%id_Ed_0_0475,type_surface_standard_variable(name='surf_direct_downward_irradiance_0475_nm'))
      call self%register_dependency(self%id_Ed_0_0500,type_surface_standard_variable(name='surf_direct_downward_irradiance_0500_nm'))
      call self%register_dependency(self%id_Ed_0_0525,type_surface_standard_variable(name='surf_direct_downward_irradiance_0525_nm'))
      call self%register_dependency(self%id_Ed_0_0550,type_surface_standard_variable(name='surf_direct_downward_irradiance_0550_nm'))
      call self%register_dependency(self%id_Ed_0_0575,type_surface_standard_variable(name='surf_direct_downward_irradiance_0575_nm'))
      call self%register_dependency(self%id_Ed_0_0600,type_surface_standard_variable(name='surf_direct_downward_irradiance_0600_nm'))
      call self%register_dependency(self%id_Ed_0_0625,type_surface_standard_variable(name='surf_direct_downward_irradiance_0625_nm'))
      call self%register_dependency(self%id_Ed_0_0650,type_surface_standard_variable(name='surf_direct_downward_irradiance_0650_nm'))
      call self%register_dependency(self%id_Ed_0_0675,type_surface_standard_variable(name='surf_direct_downward_irradiance_0675_nm'))
      call self%register_dependency(self%id_Ed_0_0700,type_surface_standard_variable(name='surf_direct_downward_irradiance_0700_nm'))
      call self%register_dependency(self%id_Ed_0_0725,type_surface_standard_variable(name='surf_direct_downward_irradiance_0725_nm'))
      call self%register_dependency(self%id_Ed_0_0775,type_surface_standard_variable(name='surf_direct_downward_irradiance_0775_nm'))
      call self%register_dependency(self%id_Ed_0_0850,type_surface_standard_variable(name='surf_direct_downward_irradiance_0850_nm'))
      call self%register_dependency(self%id_Ed_0_0950,type_surface_standard_variable(name='surf_direct_downward_irradiance_0950_nm'))
      call self%register_dependency(self%id_Ed_0_1050,type_surface_standard_variable(name='surf_direct_downward_irradiance_1050_nm'))
      call self%register_dependency(self%id_Ed_0_1150,type_surface_standard_variable(name='surf_direct_downward_irradiance_1150_nm'))
      call self%register_dependency(self%id_Ed_0_1250,type_surface_standard_variable(name='surf_direct_downward_irradiance_1250_nm'))
      call self%register_dependency(self%id_Ed_0_1350,type_surface_standard_variable(name='surf_direct_downward_irradiance_1350_nm'))
      call self%register_dependency(self%id_Ed_0_1450,type_surface_standard_variable(name='surf_direct_downward_irradiance_1450_nm'))
      call self%register_dependency(self%id_Ed_0_1550,type_surface_standard_variable(name='surf_direct_downward_irradiance_1550_nm'))
      call self%register_dependency(self%id_Ed_0_1650,type_surface_standard_variable(name='surf_direct_downward_irradiance_1650_nm'))
      call self%register_dependency(self%id_Ed_0_1750,type_surface_standard_variable(name='surf_direct_downward_irradiance_1750_nm'))
      call self%register_dependency(self%id_Ed_0_1900,type_surface_standard_variable(name='surf_direct_downward_irradiance_1900_nm'))
      call self%register_dependency(self%id_Ed_0_2200,type_surface_standard_variable(name='surf_direct_downward_irradiance_2200_nm'))
      call self%register_dependency(self%id_Ed_0_2900,type_surface_standard_variable(name='surf_direct_downward_irradiance_2900_nm'))
      call self%register_dependency(self%id_Ed_0_3700,type_surface_standard_variable(name='surf_direct_downward_irradiance_3700_nm'))
      call self%register_dependency(self%id_Es_0_0250,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0250_nm'))
      call self%register_dependency(self%id_Es_0_0325,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0325_nm'))
      call self%register_dependency(self%id_Es_0_0350,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0350_nm'))
      call self%register_dependency(self%id_Es_0_0375,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0375_nm'))
      call self%register_dependency(self%id_Es_0_0400,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0400_nm'))
      call self%register_dependency(self%id_Es_0_0425,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0425_nm'))
      call self%register_dependency(self%id_Es_0_0450,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0450_nm'))
      call self%register_dependency(self%id_Es_0_0475,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0475_nm'))
      call self%register_dependency(self%id_Es_0_0500,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0500_nm'))
      call self%register_dependency(self%id_Es_0_0525,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0525_nm'))
      call self%register_dependency(self%id_Es_0_0550,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0550_nm'))
      call self%register_dependency(self%id_Es_0_0575,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0575_nm'))
      call self%register_dependency(self%id_Es_0_0600,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0600_nm'))
      call self%register_dependency(self%id_Es_0_0625,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0625_nm'))
      call self%register_dependency(self%id_Es_0_0650,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0650_nm'))
      call self%register_dependency(self%id_Es_0_0675,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0675_nm'))
      call self%register_dependency(self%id_Es_0_0700,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0700_nm'))
      call self%register_dependency(self%id_Es_0_0725,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0725_nm'))
      call self%register_dependency(self%id_Es_0_0775,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0775_nm'))
      call self%register_dependency(self%id_Es_0_0850,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0850_nm'))
      call self%register_dependency(self%id_Es_0_0950,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0950_nm'))
      call self%register_dependency(self%id_Es_0_1050,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1050_nm'))
      call self%register_dependency(self%id_Es_0_1150,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1150_nm'))
      call self%register_dependency(self%id_Es_0_1250,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1250_nm'))
      call self%register_dependency(self%id_Es_0_1350,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1350_nm'))
      call self%register_dependency(self%id_Es_0_1450,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1450_nm'))
      call self%register_dependency(self%id_Es_0_1550,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1550_nm'))
      call self%register_dependency(self%id_Es_0_1650,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1650_nm'))
      call self%register_dependency(self%id_Es_0_1750,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1750_nm'))
      call self%register_dependency(self%id_Es_0_1900,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1900_nm'))
      call self%register_dependency(self%id_Es_0_2200,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_2200_nm'))
      call self%register_dependency(self%id_Es_0_2900,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_2900_nm'))
      call self%register_dependency(self%id_Es_0_3700,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_3700_nm'))
! END BLOCK2 python generated code
   end subroutine initialize

   subroutine do_column(self,_ARGUMENTS_VERTICAL_)
      class (type_ogs_bfm_light_spectral),intent(in) :: self
      _DECLARE_ARGUMENTS_VERTICAL_

      integer  :: kk,nlev,l
      real(rk) :: dz,zenithA,mud
      real(rk) :: phy_a,phy_b,phy_bb
      real(rk) :: cdom_a
      real(rk) :: tot_a,tot_b,tot_bb
      real(rk) :: R6c,X1c,X2c,X3c
      real(rk) :: P1c,P2c,P3c,P4c,P5c,P6c,P7c,P8c,P9c
      real(rk) :: P1chl,P2chl,P3chl,P4chl
      real(rk) :: P5chl,P6chl,P7chl,P8chl,P9chl
!      real(rk) :: aph450, anap450
!      real(rk) :: acdom250,acdom325,acdom400,acdom425,acdom450
!      real(rk) :: Scdom350_500, Scdom250_325
!      real(rk) :: bbp450, bbp550, bbp700
!      real(rk) :: equis(self%nlt),ies(self%nlt)
      real(rk) :: rlamm
      real(rk) :: zgrid(cache%n+1)
      real(rk) :: a_array(cache%n, self%nlt)
      real(rk) :: b_array(cache%n, self%nlt)
      real(rk) :: bb_array(cache%n, self%nlt)
      real(rk) :: vd(cache%n, self%nlt)
      real(rk) :: E(3, cache%n+1, self%nlt)
      real(rk) :: E_ave(3, cache%n, self%nlt)
      real(rk) :: rd, rs, ru, vs, vu
      real(rk) :: E_scalar(cache%n, self%nlt)
      real(rk) :: PAR_P1_array(cache%n)
      real(rk) :: PAR_P2_array(cache%n)
      real(rk) :: PAR_P3_array(cache%n)
      real(rk) :: PAR_P4_array(cache%n)
      real(rk) :: PAR_P5_array(cache%n)
      real(rk) :: PAR_P6_array(cache%n)
      real(rk) :: PAR_P7_array(cache%n)
      real(rk) :: PAR_P8_array(cache%n)
      real(rk) :: PAR_P9_array(cache%n)
      real(rk) :: PAR_scalar_array(cache%n)
      real(rk) :: PAR_P1,PAR_P4
      real(rk) :: PAR_P2,PAR_P5,PAR_P7,PAR_P8
      real(rk) :: PAR_P3,PAR_P6,PAR_P9
      real(rk) :: PAR_scalar

      _GET_HORIZONTAL_(self%id_zenithA,zenithA)   ! Zenith angle
      call getrmud(zenithA,mud) ! average cosine direct component in the water

!START BLOCK3
      _GET_SURFACE_(self%id_Ed_0_0250,Ed_0(1))
      _GET_SURFACE_(self%id_Ed_0_0325,Ed_0(2))
      _GET_SURFACE_(self%id_Ed_0_0350,Ed_0(3))
      _GET_SURFACE_(self%id_Ed_0_0375,Ed_0(4))
      _GET_SURFACE_(self%id_Ed_0_0400,Ed_0(5))
      _GET_SURFACE_(self%id_Ed_0_0425,Ed_0(6))
      _GET_SURFACE_(self%id_Ed_0_0450,Ed_0(7))
      _GET_SURFACE_(self%id_Ed_0_0475,Ed_0(8))
      _GET_SURFACE_(self%id_Ed_0_0500,Ed_0(9))
      _GET_SURFACE_(self%id_Ed_0_0525,Ed_0(10))
      _GET_SURFACE_(self%id_Ed_0_0550,Ed_0(11))
      _GET_SURFACE_(self%id_Ed_0_0575,Ed_0(12))
      _GET_SURFACE_(self%id_Ed_0_0600,Ed_0(13))
      _GET_SURFACE_(self%id_Ed_0_0625,Ed_0(14))
      _GET_SURFACE_(self%id_Ed_0_0650,Ed_0(15))
      _GET_SURFACE_(self%id_Ed_0_0675,Ed_0(16))
      _GET_SURFACE_(self%id_Ed_0_0700,Ed_0(17))
      _GET_SURFACE_(self%id_Ed_0_0725,Ed_0(18))
      _GET_SURFACE_(self%id_Ed_0_0775,Ed_0(19))
      _GET_SURFACE_(self%id_Ed_0_0850,Ed_0(20))
      _GET_SURFACE_(self%id_Ed_0_0950,Ed_0(21))
      _GET_SURFACE_(self%id_Ed_0_1050,Ed_0(22))
      _GET_SURFACE_(self%id_Ed_0_1150,Ed_0(23))
      _GET_SURFACE_(self%id_Ed_0_1250,Ed_0(24))
      _GET_SURFACE_(self%id_Ed_0_1350,Ed_0(25))
      _GET_SURFACE_(self%id_Ed_0_1450,Ed_0(26))
      _GET_SURFACE_(self%id_Ed_0_1550,Ed_0(27))
      _GET_SURFACE_(self%id_Ed_0_1650,Ed_0(28))
      _GET_SURFACE_(self%id_Ed_0_1750,Ed_0(29))
      _GET_SURFACE_(self%id_Ed_0_1900,Ed_0(30))
      _GET_SURFACE_(self%id_Ed_0_2200,Ed_0(31))
      _GET_SURFACE_(self%id_Ed_0_2900,Ed_0(32))
      _GET_SURFACE_(self%id_Ed_0_3700,Ed_0(33))
      _GET_SURFACE_(self%id_Es_0_0250,Es_0(1))
      _GET_SURFACE_(self%id_Es_0_0325,Es_0(2))
      _GET_SURFACE_(self%id_Es_0_0350,Es_0(3))
      _GET_SURFACE_(self%id_Es_0_0375,Es_0(4))
      _GET_SURFACE_(self%id_Es_0_0400,Es_0(5))
      _GET_SURFACE_(self%id_Es_0_0425,Es_0(6))
      _GET_SURFACE_(self%id_Es_0_0450,Es_0(7))
      _GET_SURFACE_(self%id_Es_0_0475,Es_0(8))
      _GET_SURFACE_(self%id_Es_0_0500,Es_0(9))
      _GET_SURFACE_(self%id_Es_0_0525,Es_0(10))
      _GET_SURFACE_(self%id_Es_0_0550,Es_0(11))
      _GET_SURFACE_(self%id_Es_0_0575,Es_0(12))
      _GET_SURFACE_(self%id_Es_0_0600,Es_0(13))
      _GET_SURFACE_(self%id_Es_0_0625,Es_0(14))
      _GET_SURFACE_(self%id_Es_0_0650,Es_0(15))
      _GET_SURFACE_(self%id_Es_0_0675,Es_0(16))
      _GET_SURFACE_(self%id_Es_0_0700,Es_0(17))
      _GET_SURFACE_(self%id_Es_0_0725,Es_0(18))
      _GET_SURFACE_(self%id_Es_0_0775,Es_0(19))
      _GET_SURFACE_(self%id_Es_0_0850,Es_0(20))
      _GET_SURFACE_(self%id_Es_0_0950,Es_0(21))
      _GET_SURFACE_(self%id_Es_0_1050,Es_0(22))
      _GET_SURFACE_(self%id_Es_0_1150,Es_0(23))
      _GET_SURFACE_(self%id_Es_0_1250,Es_0(24))
      _GET_SURFACE_(self%id_Es_0_1350,Es_0(25))
      _GET_SURFACE_(self%id_Es_0_1450,Es_0(26))
      _GET_SURFACE_(self%id_Es_0_1550,Es_0(27))
      _GET_SURFACE_(self%id_Es_0_1650,Es_0(28))
      _GET_SURFACE_(self%id_Es_0_1750,Es_0(29))
      _GET_SURFACE_(self%id_Es_0_1900,Es_0(30))
      _GET_SURFACE_(self%id_Es_0_2200,Es_0(31))
      _GET_SURFACE_(self%id_Es_0_2900,Es_0(32))
      _GET_SURFACE_(self%id_Es_0_3700,Es_0(33))

!END BLOCK3 python generated code
      
      kk=0
      zgrid(1)=0.0_rk

      _DOWNWARD_LOOP_BEGIN_
          kk = kk + 1
       _GET_(self%id_dz,dz)     ! Layer height (m)

         zgrid(kk+1)=zgrid(kk)+dz

       _GET_(self%id_aP1chl,P1chl)
       _GET_(self%id_aP2chl,P2chl)
       _GET_(self%id_aP3chl,P3chl)
       _GET_(self%id_aP4chl,P4chl)
       _GET_(self%id_aP5chl,P5chl)
       _GET_(self%id_aP6chl,P6chl)
       _GET_(self%id_aP7chl,P7chl)
       _GET_(self%id_aP8chl,P8chl)
       _GET_(self%id_aP9chl,P9chl)

       _GET_(self%id_aP1c,P1c)
       _GET_(self%id_aP2c,P2c)
       _GET_(self%id_aP3c,P3c)
       _GET_(self%id_aP4c,P4c)       
       _GET_(self%id_aP5c,P5c)
       _GET_(self%id_aP6c,P6c)
       _GET_(self%id_aP7c,P7c)
       _GET_(self%id_aP8c,P8c)       
       _GET_(self%id_aP9c,P9c)       
       
       _GET_(self%id_R6c,R6c)

       _GET_(self%id_X1c,X1c)
       _GET_(self%id_X2c,X2c)
       _GET_(self%id_X3c,X3c)

          _SET_DIAGNOSTIC_(self%id_P1chl, max(p_small,P1chl))
          _SET_DIAGNOSTIC_(self%id_P2chl, max(p_small,P2chl))       
          _SET_DIAGNOSTIC_(self%id_P3chl, max(p_small,P3chl))       
          _SET_DIAGNOSTIC_(self%id_P4chl, max(p_small,P4chl))       
          _SET_DIAGNOSTIC_(self%id_P5chl, max(p_small,P5chl))       
          _SET_DIAGNOSTIC_(self%id_P6chl, max(p_small,P6chl))       
          _SET_DIAGNOSTIC_(self%id_P7chl, max(p_small,P7chl))       
          _SET_DIAGNOSTIC_(self%id_P8chl, max(p_small,P8chl))       
          _SET_DIAGNOSTIC_(self%id_P9chl, max(p_small,P9chl))       
          _SET_DIAGNOSTIC_(self%id_P1c, max(p_small,P1c))
          _SET_DIAGNOSTIC_(self%id_P2c, max(p_small,P2c))       
          _SET_DIAGNOSTIC_(self%id_P3c, max(p_small,P3c))       
          _SET_DIAGNOSTIC_(self%id_P4c, max(p_small,P4c))       
          _SET_DIAGNOSTIC_(self%id_P5c, max(p_small,P5c))       
          _SET_DIAGNOSTIC_(self%id_P6c, max(p_small,P6c))       
          _SET_DIAGNOSTIC_(self%id_P7c, max(p_small,P7c))       
          _SET_DIAGNOSTIC_(self%id_P8c, max(p_small,P8c))       
          _SET_DIAGNOSTIC_(self%id_P9c, max(p_small,P9c))       
       
! Equations determining optical properties in relations to biogeochemical variables
          do l=1,self%nlt
             phy_a  = ac(1,l)*P1chl+ac(2,l)*P2chl+ac(3,l)*P3chl+ac(4,l)*P4chl +ac(2,l)*P5chl+ac(3,l)*P6chl+ac(2,l)*P7chl+ac(2,l)*P8chl+ac(3,l)*P9chl
             phy_b  = bc(1,l)*P1c + bc(2,l)*P2c + bc(3,l)*P3c + bc(4,l)*P4c   +bc(2,l)*P5c + bc(3,l)*P6c + bc(2,l)*P7c + bc(2,l)*P8c+ bc(3,l)*P9c
!             phy_b  = bc(1,l)*P1chl + bc(2,l)*P2chl + bc(3,l)*P3chl + bc(4,l)*P4chl
             phy_bb = bc(1,l)*bbc(1,l)*P1c + bc(2,l)*bbc(2,l)*P2c + bc(3,l)*bbc(3,l)*P3c + bc(4,l)*bbc(4,l)*P4c + bc(2,l)*bbc(2,l)*P5c + bc(3,l)*bbc(3,l)*P6c + bc(2,l)*bbc(2,l)*P7c + bc(2,l)*bbc(2,l)*P8c + bc(3,l)*bbc(3,l)*P9c
!             phy_bb = bc(1,l)*bbc(1,l)*P1chl + bc(2,l)*bbc(2,l)*P2chl + bc(3,l)*bbc(3,l)*P3chl + bc(4,l)*bbc(4,l)*P4chl
             cdom_a = acdom(1,l)*X1c  + acdom(2,l)*X2c  + acdom(3,l)*X3c 
             cdom_a = MAX(cdom_a, acdom_min(l))
             
! Need to add also cdom
             tot_a  =  aw(l) + phy_a  + apoc(l) * R6c + cdom_a
             tot_b  =  bw(l) + phy_b  + bpoc(l) * R6c 
             tot_bb = bbw(l) + phy_bb + bbpoc(l)* R6c 

             a_array(kk,l)  = tot_a
             b_array(kk,l)  = tot_b
             bb_array(kk,l) = tot_bb
          enddo

!! IOPs observed for diagnostics
!       acdom250 = MAX(acdom(1,1)*X1c + acdom(2,1)*X2c + acdom(3,1)*X3c, acdom_min(1))
!       acdom325 = MAX(acdom(1,2)*X1c + acdom(2,2)*X2c + acdom(3,2)*X3c, acdom_min(2))
!       acdom400 = MAX(acdom(1,5)*X1c + acdom(2,5)*X2c + acdom(3,5)*X3c, acdom_min(5))
!       acdom425 = MAX(acdom(1,6)*X1c + acdom(2,6)*X2c + acdom(3,6)*X3c, acdom_min(6))
!       acdom450 = MAX(acdom(1,7)*X1c + acdom(2,7)*X2c + acdom(3,7)*X3c, acdom_min(7))
!       aph450   = ac(1,7)*P1chl + ac(2,7)*P2chl + ac(3,7)*P3chl + ac(4,7)*P4chl + ac(2,7)*P5chl + ac(3,7)*P6chl + ac(2,7)*P7chl + ac(2,7)*P8chl + ac(3,7)*P9chl
!       anap450  = apoc(7) * R6c
!       bbp450 = (bc(1,7)*bbc(1,7)*P1c + bc(2,7)*bbc(2,7)*P2c + bc(3,7)*bbc(3,7)*P3c + bc(4,7)*bbc(4,7)*P4c + bc(2,7)*bbc(2,7)*P5c + bc(3,7)*bbc(3,7)*P6c + bc(2,7)*bbc(2,7)*P7c + bc(2,7)*bbc(2,7)*P8c + bc(3,7)*bbc(3,7)*P9c) + bbpoc(7)*R6c
!       bbp550 = (bc(1,11)*bbc(1,11)*P1c + bc(2,11)*bbc(2,11)*P2c + bc(3,11)*bbc(3,11)*P3c + bc(4,11)*bbc(4,11)*P4c + bc(2,11)*bbc(2,11)*P5c + bc(3,11)*bbc(3,11)*P6c + bc(2,11)*bbc(2,11)*P7c + bc(2,11)*bbc(2,11)*P8c + bc(3,11)*bbc(3,11)*P9c) + bbpoc(11)*R6c
!       bbp700 = (bc(1,17)*bbc(1,17)*P1c + bc(2,17)*bbc(2,17)*P2c + bc(3,17)*bbc(3,17)*P3c + bc(4,17)*bbc(4,17)*P4c + bc(2,17)*bbc(2,17)*P5c + bc(3,17)*bbc(3,17)*P6c + bc(2,17)*bbc(2,17)*P7c + bc(2,17)*bbc(2,17)*P8c + bc(3,17)*bbc(3,17)*P9c) + bbpoc(17)*R6c

!         _SET_DIAGNOSTIC_(self%id_acdom250, max(p_small,acdom250))            
!         _SET_DIAGNOSTIC_(self%id_acdom325, max(p_small,acdom325))            
!         _SET_DIAGNOSTIC_(self%id_acdom400, max(p_small,acdom400))            
!         _SET_DIAGNOSTIC_(self%id_acdom425, max(p_small,acdom425))            
!         _SET_DIAGNOSTIC_(self%id_acdom450, max(p_small,acdom450))
!         _SET_DIAGNOSTIC_(self%id_aph450, max(p_small,aph450))              
!         _SET_DIAGNOSTIC_(self%id_anap450, max(p_small,anap450))
!         _SET_DIAGNOSTIC_(self%id_bbp450, max(p_small,bbp450))
!         _SET_DIAGNOSTIC_(self%id_bbp550, max(p_small,bbp550))
!         _SET_DIAGNOSTIC_(self%id_bbp700, max(p_small,bbp700))
         
       ! linear fit of ln-transformed aCDOM(l) against wavelength:
       ! between 350 and 500 nm (Babin et al 2013, Organelli et al 2014) 
       ! Organelli 2014 uses non-linear least-squares fit!
       ! between 250 and 325 nm (Catalá et al 2018, Galletti et al 2019) 
!          do l=1,self%nlt
!             rlamm = real(lam(l),8)     
!             equis(l) = rlamm-self%lambda_aCDOM
!             ies(l) = LOG(acdom(1,l)*X1c+acdom(2,l)*X2c+acdom(3,l)*X3c)
!          enddo
!          call linear_regression(equis(3:9),ies(3:9),7,acdom450,Scdom350_500) ! indexes for 350-500nm
!          call linear_regression(equis(1:2),ies(1:2),7,acdom450,Scdom250_325) ! indexes for 250-325nm
!       ! acdom450=EXP(acdom450)
!         _SET_DIAGNOSTIC_(self%id_Scdom350_500, -Scdom350_500) 
!         _SET_DIAGNOSTIC_(self%id_Scdom250_325, -Scdom250_325) 
       
     _DOWNWARD_LOOP_END_
  
     rd      = self%rd
     rs      = self%rs
     ru      = self%ru
     vd(:,:) = mud
     vs      = self%vs
     vu      = self%vu

     call solve_direct(cache%n+1, zgrid, cache%n, zgrid, self%nlt, a_array, b_array, bb_array, rd, rs, ru, vd, vs, vu, Ed_0, Es_0, E, E_ave)


! Scalar irradiance
     E_scalar(:,:)=E_ave(1,:,:)/vd + E_ave(2,:,:)/vs + E_ave(3,:,:)/vu

     PAR_P1_array(:)       = 0.0_rk
     PAR_P2_array(:)       = 0.0_rk
     PAR_P3_array(:)       = 0.0_rk
     PAR_P4_array(:)       = 0.0_rk
     PAR_P5_array(:)       = 0.0_rk
     PAR_P6_array(:)       = 0.0_rk
     PAR_P7_array(:)       = 0.0_rk
     PAR_P8_array(:)       = 0.0_rk
     PAR_P9_array(:)       = 0.0_rk
     PAR_scalar_array(:)   = 0.0_rk

!     do l=1,self%nlt
     do l=5,17     
         PAR_P1_array(:) = PAR_P1_array(:) + (WtoQ(l) * ac_ps(1,l) * E_scalar(:,l)) * SEC_PER_DAY
         PAR_P2_array(:) = PAR_P2_array(:) + (WtoQ(l) * ac_ps(2,l) * E_scalar(:,l)) * SEC_PER_DAY
         PAR_P3_array(:) = PAR_P3_array(:) + (WtoQ(l) * ac_ps(3,l) * E_scalar(:,l)) * SEC_PER_DAY
         PAR_P4_array(:) = PAR_P4_array(:) + (WtoQ(l) * ac_ps(4,l) * E_scalar(:,l)) * SEC_PER_DAY
         PAR_P5_array(:) = PAR_P5_array(:) + (WtoQ(l) * ac_ps(5,l) * E_scalar(:,l)) * SEC_PER_DAY
         PAR_P6_array(:) = PAR_P6_array(:) + (WtoQ(l) * ac_ps(6,l) * E_scalar(:,l)) * SEC_PER_DAY
         PAR_P7_array(:) = PAR_P7_array(:) + (WtoQ(l) * ac_ps(7,l) * E_scalar(:,l)) * SEC_PER_DAY
         PAR_P8_array(:) = PAR_P8_array(:) + (WtoQ(l) * ac_ps(8,l) * E_scalar(:,l)) * SEC_PER_DAY
         PAR_P9_array(:) = PAR_P9_array(:) + (WtoQ(l) * ac_ps(9,l) * E_scalar(:,l)) * SEC_PER_DAY         
     enddo

     do l=5,17
         PAR_scalar_array(:)      = PAR_scalar_array(:) + (E_scalar(:,l) * WtoQ(l)) * SEC_PER_DAY
     enddo

! AOPs observed for diagnostics     
!    write(*,*) 'Rrs= ', E(3,1,7)/(Ed_0(7)+Es_0(7))

!     _HORIZONTAL_LOOP_BEGIN_
!     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs400,E(3,1,5)/max(p_small,(E(1,1,5)+E(2,1,5))))
!     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs425,E(3,1,6)/max(p_small,(E(1,1,6)+E(2,1,6))))
!     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs450,E(3,1,7)/max(p_small,(E(1,1,7)+E(2,1,7)))) 
!     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs475,E(3,1,8)/max(p_small,(E(1,1,8)+E(2,1,8))))
!     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs500,E(3,1,9)/max(p_small,(E(1,1,9)+E(2,1,9))))
!     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs525,E(3,1,10)/max(p_small,(E(1,1,10)+E(2,1,10)))) 
!     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs550,E(3,1,11)/max(p_small,(E(1,1,11)+E(2,1,11)))) 
!     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs575,E(3,1,12)/max(p_small,(E(1,1,12)+E(2,1,12)))) 
!     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs675,E(3,1,16)/max(p_small,(E(1,1,16)+E(2,1,16)))) 

!     write(*,*) 'Z9= ', zgrid(26)
     
!     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_kd375,-LOG(max(p_small,(E(1,26,4)+E(2,26,4))/(max(p_small,E(1,1,4)+E(2,1,4)))))/9.05_rk)
!     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_kd400,-LOG(max(p_small,(E(1,26,5)+E(2,26,5))/(max(p_small,E(1,1,5)+E(2,1,5)))))/9.05_rk)     
!     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_kd425,-LOG(max(p_small,(E(1,26,6)+E(2,26,6))/(max(p_small,E(1,1,6)+E(2,1,6)))))/9.05_rk)
!     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_kd475,-LOG(max(p_small,(E(1,26,8)+E(2,26,8))/(max(p_small,E(1,1,8)+E(2,1,8)))))/9.05_rk)
!     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_kd500,-LOG(max(p_small,(E(1,26,9)+E(2,26,9))/(max(p_small,E(1,1,9)+E(2,1,9)))))/9.05_rk)     
     
!      _HORIZONTAL_LOOP_END_

      kk=0

      _DOWNWARD_LOOP_BEGIN_

          kk = kk + 1

         _SET_DIAGNOSTIC_(self%id_par_P1, max(p_small,PAR_P1_array(kk)))                  
         _SET_DIAGNOSTIC_(self%id_par_P2, max(p_small,PAR_P2_array(kk)))            
         _SET_DIAGNOSTIC_(self%id_par_P3, max(p_small,PAR_P3_array(kk)))
         _SET_DIAGNOSTIC_(self%id_par_P4, max(p_small,PAR_P4_array(kk)))
         _SET_DIAGNOSTIC_(self%id_par_P5, max(p_small,PAR_P5_array(kk)))                  
         _SET_DIAGNOSTIC_(self%id_par_P6, max(p_small,PAR_P6_array(kk)))            
         _SET_DIAGNOSTIC_(self%id_par_P7, max(p_small,PAR_P7_array(kk)))
         _SET_DIAGNOSTIC_(self%id_par_P8, max(p_small,PAR_P8_array(kk)))          
         _SET_DIAGNOSTIC_(self%id_par_P9, max(p_small,PAR_P9_array(kk)))         
         _SET_DIAGNOSTIC_(self%id_PAR_tot, max(p_small,PAR_scalar_array(kk)))          

     _DOWNWARD_LOOP_END_

   end subroutine do_column

end module
