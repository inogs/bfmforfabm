#include "fabm_driver.h"

module ogs_bfm_light_spectral

   use fabm_types
   use ogs_bfm_shared
   use adj_3stream, only: solve_direct

   implicit none

   private

   type,extends(type_base_model),public :: type_ogs_bfm_light_spectral
      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id)   :: id_prova1, id_prova2
      type (type_diagnostic_variable_id)   :: id_tot_a, id_tot_b, id_tot_bb
      type (type_diagnostic_variable_id)   :: id_phy_a, id_phy_b, id_phy_bb
      type (type_dependency_id)            :: id_dz
      type (type_state_variable_id)        :: id_P1c, id_P2c, id_P3c, id_P4c
      type (type_state_variable_id)        :: id_P1chl, id_P2chl, id_P3chl, id_P4chl
      type (type_state_variable_id)        :: id_R6c, id_X1c, id_X2c, id_X3c
      type (type_diagnostic_variable_id)   :: id_par_dia, id_par_flag, id_par_pico, id_par_dino
      type (type_diagnostic_variable_id)   :: id_par_tot
! BLOCK 1 python generated  code
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
      integer  :: nlt
!     real(rk) :: EPSESSX,EPS0X,pEIR_eowX
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
      real(rk) :: hc, hcoavo, rlamm, nl
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%nlt,    'nlt',   '-',   'number of wavelenghts', default=-1)
       
      if (self%nlt>0) then
          allocate(lam(self%nlt)); lam(:)=huge(lam(1))
          allocate(aw(self%nlt));   aw(:)=huge(aw(1))
          allocate(bw(self%nlt));   bw(:)=huge(bw(1))
          allocate(bbw(self%nlt));  bbw(:)=huge(bbw(1))
          allocate(ac(4,self%nlt));   ac(4,:)=huge(ac(1,1))
          allocate(bc(4,self%nlt));   bc(4,:)=huge(bc(1,1))
          allocate(bbc(4,self%nlt));  bbc(4,:)=huge(bbc(1,1))
          allocate(apoc(self%nlt)); apoc(:)=huge(apoc(1))
          allocate(bpoc(self%nlt)); bpoc(:)=huge(bpoc(1))
          allocate(bbpoc(self%nlt)); bbpoc(:)=huge(bbpoc(1))
          allocate(Ed_0(self%nlt)); Ed_0(:)=huge(Ed_0(1))
          allocate(Es_0(self%nlt)); Es_0(:)=huge(Es_0(1))
          allocate(WtoQ(self%nlt)); WtoQ(:)=huge(WtoQ(1))

!         load the IOP for the biogeochemical variables considered
          call lidata(self%nlt,4)

      endif 

     hc = 1.0D0/(h_planck*c_light)
     hcoavo = hc*oavo


      do nl = 1,self%nlt
       rlamm = real(lam(nl),8)*1.0E-9      !lambda in m
       WtoQ(nl) = rlamm*hcoavo*1000000.0D0 !Watts to micro mol quanta conversion
      enddo

     ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_prova1,'prova1','W/m^2','prova1_paolo', source=source_do_column)
      call self%register_diagnostic_variable(self%id_tot_a,'tot_a','1/m','total absorption', source=source_do_column)
      call self%register_diagnostic_variable(self%id_tot_b,'tot_b','1/m','total scattering', source=source_do_column)
      call self%register_diagnostic_variable(self%id_tot_bb,'tot_bb','1/m','total back scattering', source=source_do_column)
      call self%register_diagnostic_variable(self%id_phy_a,'phy_a','1/m','phytoplankton absorption', source=source_do_column)
      call self%register_diagnostic_variable(self%id_phy_b,'phy_b','1/m','phytoplankton scattering', source=source_do_column)
      call self%register_diagnostic_variable(self%id_phy_bb,'phy_bb','1/m','phytoplankton back scattering', source=source_do_column)

      call self%register_diagnostic_variable(self%id_par_dia, 'PAR_dia','?????','PAR_diatoms', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_flag,'PAR_flag','?????','PAR_flagellates', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_pico,'PAR_pico','?????','PAR_picophytoplankton', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_dino,'PAR_dino','?????','PAR_dinoflagellates', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_tot,'PAR_tot','?????','PAR_total', source=source_do_column)

      ! Register biogeochemical dependencies 

      call self%register_state_dependency(self%id_P1c,'P1c','mg C/m^3', 'Diatoms carbon')
      call self%register_state_dependency(self%id_P2c,'P2c','mg C/m^3', 'Flagellates carbon')
      call self%register_state_dependency(self%id_P3c,'P3c','mg C/m^3', 'PicoPhytoplankton carbon')
      call self%register_state_dependency(self%id_P4c,'P4c','mg C/m^3', 'DinoFlagellates carbon')

      call self%register_state_dependency(self%id_P1chl,'P1chl','mg chl/m^3', 'Diatoms chlorophyll')
      call self%register_state_dependency(self%id_P2chl,'P2chl','mg chl/m^3', 'Flagellates chlorophyll')
      call self%register_state_dependency(self%id_P3chl,'P3chl','mg chl/m^3', 'PicoPhytoplankton chlorophyll')
      call self%register_state_dependency(self%id_P4chl,'P4chl','mg chl/m^3', 'DinoFlagellates chlorophyll')

      call self%register_state_dependency(self%id_R6c,'R6c','mg C/m^3', 'POC')
      call self%register_state_dependency(self%id_X1c,'X1c','mg C/m^3', 'labile CDOM')
      call self%register_state_dependency(self%id_X2c,'X2c','mg C/m^3', 'semi-labile CDOM')
      call self%register_state_dependency(self%id_X3c,'X3c','mg C/m^3', 'semi-refractory CDOM')

      ! Register environmental dependencies 
      call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
! BLOCK 1 python generate code see AUX_SCRIPTS/python_light_spectral.py
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
! END BLOCK1
   end subroutine initialize

   subroutine do_column(self,_ARGUMENTS_VERTICAL_)
      class (type_ogs_bfm_light_spectral),intent(in) :: self
      _DECLARE_ARGUMENTS_VERTICAL_

      integer  :: kk,nlev,l
      real(rk) :: buffer,dz,xEPS,xtnc,EIR,ESS
      real(rk) :: phy_a,phy_b,phy_bb
      real(rk) :: cdom_a
      real(rk) :: tot_a,tot_b,tot_bb
      real(rk) :: R6c,X1c,X2c,X3c
      real(rk) :: P1c,P2c,P3c,P4c
      real(rk) :: P1chl, P2chl, P3chl, P4chl 
      real(rk) :: zgrid(cache%n+1)
      real(rk) :: a_array(cache%n, self%nlt)
      real(rk) :: b_array(cache%n, self%nlt)
      real(rk) :: bb_array(cache%n, self%nlt)
      real(rk) :: vd(cache%n, self%nlt)
      real(rk) :: E(3, cache%n+1, self%nlt)
      real(rk) :: E_ave(3, cache%n, self%nlt)
      real(rk) :: rd, rs, ru, vs, vu 
      real(rk) :: E_scalar(cache%n, self%nlt)
      real(rk) :: PAR_diatoms_array(cache%n)
      real(rk) :: PAR_flagellates_array(cache%n)
      real(rk) :: PAR_picophytoplankton_array(cache%n)
      real(rk) :: PAR_dinoflagellates_array(cache%n)
      real(rk) :: PAR_scalar_array(cache%n)
      real(rk) :: PAR_diatoms
      real(rk) :: PAR_flagellates
      real(rk) :: PAR_picophytoplankton
      real(rk) :: PAR_dinoflagellates
      real(rk) :: PAR_scalar
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

!END BLOCK3
      

      write(*,*) "Ed_0_250", Ed_0(1)
      write(*,*) "Es_0_250", Es_0(1)
      write(*,*) "Ed_0_325", Ed_0(2)
      write(*,*) "Es_0_325", Es_0(2)
      write(*,*) " size(a_array) ", size(a_array)


!     if (buffer.lt.0._rk) buffer=0._rk
      kk=0
      zgrid(1)=0.0_rk
      _DOWNWARD_LOOP_BEGIN_
!     _VERTICAL_LOOP_BEGIN_
          kk = kk + 1
       _GET_(self%id_dz,dz)     ! Layer height (m)

         zgrid(kk+1)=zgrid(kk)+dz

       _GET_(self%id_P1chl,P1chl)
       _GET_(self%id_P2chl,P2chl)
       _GET_(self%id_P3chl,P3chl)
       _GET_(self%id_P4chl,P4chl)

       _GET_(self%id_R6c,R6c)

       _GET_(self%id_X1c,X1c)
       _GET_(self%id_X2c,X2c)
       _GET_(self%id_X3c,X3c)

! Equations determining optical properties in relations to biogeochemical variables
          do l=1,self%nlt
             phy_a  = ac(1,l)*P1chl + ac(2,l)*P2chl + ac(3,l)*P3chl + ac(4,l)*P4chl
             phy_b  = bc(1,l)*P1chl + bc(2,l)*P2chl + bc(3,l)*P3chl + bc(4,l)*P4chl
             phy_bb = bbc(1,l)*P1chl+ bbc(2,l)*P2chl+ bbc(3,l)*P3chl+ bbc(4,l)*P4chl
!            cdom_a = acdom(l)*X1c  + acdom(l)*X2c  + acdom(l)*X3c 

! Need to add also cdom
             tot_a  =  aw(l) + phy_a  + apoc(l) * R6c !+ cdom_a
             tot_b  =  bw(l) + phy_b  + bpoc(l) * R6c 
             tot_bb = bbw(l) + phy_bb + bbpoc(l)* R6c 

             a_array(kk,l)  = tot_a
             b_array(kk,l)  = tot_b
             bb_array(kk,l) = tot_bb
          enddo

     _DOWNWARD_LOOP_END_
  
     rd=1.0_rk
     rs=1.0_rk
     ru=1.0_rk
     vd(:,:)=1.0_rk
     vs=1.0_rk
     vu=1.0_rk

!integer, intent(in):: n,m                                             !n of layers and size of the vertical grid
!integer, intent(in):: nlt                                             !n of wavelenghts to be considered
!double precision, dimension(m), intent(in):: zz                       !must be some grid from z(1)=0 to z(n+1)=bottom
!double precision, dimension(n+1), intent(in):: z                      !layer boundaries (depth levels); z(1)=0 (must be), z(n+1) = bottom
!double precision, dimension(n,nlt), intent(in):: a, b, bb, vd         !input depth-dependent data
!double precision, intent(in):: rd, rs, ru, vs, vu                     !input parameters
!double precision, dimension(nlt),intent(in) :: EdOASIM, EsOASIM !boundary values
!double precision, dimension(3,m,nlt), intent(out):: E                     !the 3-stream solution on the zz grid

     bb_array(:,:)=0.001_rk
     write(*,*) 'a_array(:,1)',  a_array(:,1)
     write(*,*) 'b_array(:,1)',  b_array(:,1)
     write(*,*) 'bb_array(:,1)', bb_array(:,1)

     call solve_direct(cache%n+1, zgrid, cache%n, zgrid, self%nlt, a_array, b_array, bb_array, rd, rs, ru, vd, vs, vu, Ed_0, Es_0, E, E_ave)

     write(*,*) 'E(1,:,:)', E(1,:,1:2)
     write(*,*) 'E(2,:,:)', E(2,:,1:2)
     write(*,*) 'E(3,:,:)', E(3,:,1:2)
     write(*,*) 'E_ave(1,:,:)', E_ave(1,:,1:2)
     write(*,*) 'E_ave(2,:,:)', E_ave(2,:,1:2)
     write(*,*) 'E_ave(3,:,:)', E_ave(3,:,1:2)
! Scalar irradiance
     E_scalar(:,:)=E_ave(1,:,:)/vd + E_ave(2,:,:)/vs + E_ave(3,:,:)/vu

     PAR_diatoms_array(:)           = 0.0_rk
     PAR_flagellates_array(:)       = 0.0_rk
     PAR_picophytoplankton_array(:) = 0.0_rk
     PAR_dinoflagellates_array(:)   = 0.0_rk
     PAR_scalar_array(:)            = 0.0_rk

     do l=1,self%nlt
!    do l=5,17 or 19?  
         PAR_diatoms_array(:)           = PAR_diatoms_array(:)           + WtoQ(l) * ac(1,l) * E_scalar(:,l)
         PAR_flagellates_array(:)       = PAR_flagellates_array(:)       + WtoQ(l) * ac(2,l) * E_scalar(:,l)
         PAR_picophytoplankton_array(:) = PAR_picophytoplankton_array(:) + WtoQ(l) * ac(3,l) * E_scalar(:,l)
         PAR_dinoflagellates_array(:)   = PAR_dinoflagellates_array(:)   + WtoQ(l) * ac(4,l) * E_scalar(:,l)
     enddo

     do l=5,17
!    do l=5,17 or 19?  
         PAR_scalar_array(:)            = PAR_scalar_array(:) + E_scalar(:,l)
     enddo

      kk=0

      _DOWNWARD_LOOP_BEGIN_

          kk = kk + 1

!         PAR_diatoms           = PAR_diatoms_array(kk)
!         PAR_flagellates       = PAR_flagellates_array(kk)
!         PAR_picophytoplankton = PAR_picophytoplankton_array(kk)
!         PAR_dinoflagellates   = PAR_dinoflagellates_array(kk)
!         PAR_scalar            = PAR_scalar_array(kk)

         _SET_DIAGNOSTIC_(self%id_par_dia,PAR_diatoms_array(kk))                  
         _SET_DIAGNOSTIC_(self%id_par_flag,PAR_flagellates_array(kk))            
         _SET_DIAGNOSTIC_(self%id_par_pico,PAR_picophytoplankton_array(kk))
         _SET_DIAGNOSTIC_(self%id_par_dino,PAR_dinoflagellates_array(kk))           
         _SET_DIAGNOSTIC_(self%id_par_tot,PAR_scalar_array(kk))           

     _DOWNWARD_LOOP_END_

   STOP
   end subroutine do_column

end module
