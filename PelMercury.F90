
#include "fabm_driver.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelMercury
!
! DESCRIPTION
!       This process describes the dynamics of four Mercury species
!       (HgII, Hg0, MMHg, DMHg) in the water column.
!       Parameterized processes are:
!       - partitioning of HgII and MMHHg to POC and DOC
!       - methylation of Hg and MMHg
!       - demethylation of MMHg and DMHg
!       - redox transformations of HgII and Hg0
!       - photochemical redox transformations and demethylation
!
!       MMHg bioaccumulation in the food web is handled in
!      the subroutines Phyto, MicroZoo and MesoZoo
!
! !INTERFACE
 module bfm_PelMercury

   use fabm_types
   use ogs_bfm_shared
   use ogs_bfm_pelagic_base
!
!
! !AUTHORS
!   Original version by P. Ruardij and M. Vichi
!
!
!
! !REVISION_HISTORY
!
! COPYING
!
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij, M. Vichi
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!EOP
!-------------------------------------------------------------------------!
!BOC
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  private

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     type,extends(type_ogs_bfm_pelagic_base),public :: type_ogs_bfm_PelMercury
      ! NB: own state variables (c,n,p,s,f,chl) are added implicitly by deriving
      ! from type_ogs_bfm_pelagic_base!

      ! Identifiers for state variables of other models
    !  type (type_state_variable_id) :: id_O3c,id_O2o,id_O3h,id_O4n          !  dissolved inorganic carbon, oxygen, total alkalinity, N2
      type (type_state_variable_id) :: id_Hg2,id_Hg0,id_MMHg,id_DMHg   !  nutrients: phosphate, nitrate, ammonium, silicate, iron
      type (type_state_variable_id) :: id_R1c, id_R2c,id_R3c       !  dissolved organic carbon (R1: labile, R2: semi-labile, R3: semi-refractory)
      type (type_state_variable_id) :: id_R6c,id_R8c !,id_R6n,id_R6s          !  particulate organic carbon
    !  type (type_state_variable_id) :: id_O5c   !id_R1p,id_R1n,                   !  Free calcite (liths) - used by calcifiers only
    !  type (type_state_variable_id) :: id_X1c, id_X2c, id_X3c               !  CDOM
      type (type_dependency_id)     :: id_dz
      type(type_dependency_id) :: id_ruR6c,id_ruR8c
      ! Environmental dependencies
      type (type_dependency_id)    :: id_parEIR,id_ETW   ! PAR and temperature
!      type (type_dependency_id)    :: id_PAR_tot
      type (type_dependency_id)    :: id_isBen
      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id) :: id_skmet   ! Methylation of HgII
      type (type_diagnostic_variable_id) :: id_skmet2  ! Methylation of MMHg
      type (type_diagnostic_variable_id) :: id_skdem   ! Biotic Demethylation of MMHg
      type (type_diagnostic_variable_id) :: id_skdem2  ! Biotic Demethylation of MMHg
      type (type_diagnostic_variable_id) :: id_skphdem ! Photochemical Demethylation of MMHg
      type (type_diagnostic_variable_id) :: id_skpdem2 ! Photochemical Demethylation of DMHg
      type (type_diagnostic_variable_id) :: id_skbred  ! Biotic Reduction of HgII
      type (type_diagnostic_variable_id) :: id_skphr   ! Photochemical Reduction of HgII
      type (type_diagnostic_variable_id) :: id_skphox  ! Photochemical Oxidation of Hg0
      type (type_diagnostic_variable_id) :: id_skbox  ! Biotic Oxidation of Hg0
      type (type_diagnostic_variable_id) :: id_skox    ! Dark Oxidation of Hg0
      type (type_diagnostic_variable_id) :: id_HgCl2, id_MMHgCl
      type (type_diagnostic_variable_id) :: id_HgDOC, id_MMHgDOC
      type (type_diagnostic_variable_id) :: id_HgPOC6,id_MMHgPOC6,id_HgPOC8,id_MMHgPOC8

! variation of O3h for denitrification (-1 mole of NO3 (consumed) -> + 1 mole of alk)
      ! Parameters (described in subroutine initialize, below)
       real(rk) :: p_KdhDOC, p_KdmDOC
       real(rk) :: p_KdhPOC6,p_KdhPOC8, p_KdmPOC6,p_KdmPOC8
       real(rk) :: p_kmet,p_kmet2,p_kdem, p_kphdem, p_kdem2, p_kpdem2
       real(rk) :: p_kbred,p_kbox,p_kphr,p_kphox
       integer :: p_Esource,p_use_benthic
    contains

   ! Model procedures
      procedure :: initialize
      procedure :: do
    end type type_ogs_bfm_PelMercury

contains

  subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ogs_bfm_PelMercury),intent(inout),target :: self
      integer,                        intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      ! Set time unit to d-1
      ! This implies that all rates (sink/source terms, vertical velocities) are
      ! given in d-1.
!      self%dt = 86400._rk
!EOP
!-----------------------------------------------------------------------
!BOC
! Initialize pelagic base model (this also sets the time unit to per day,instead of the default per second)
     ! call self%initialize_bfm_PelMercury
      call self%initialize_bfm_base()

     ! Obtain the values of all model parameters from FABM.
      ! Specify the long name and units of the parameters, which could be used
      ! by FABM (or its host)
      ! to present parameters to the user for configuration (e.g., through a
      ! GUI)

      call self%get_parameter(self%p_KdhPOC6,    'p_KdhPOC6',  '[l/kg]', 'Partition coefficient of Hg to POC', default=1000000.0_rk)
      call self%get_parameter(self%p_KdhPOC8,    'p_KdhPOC8',  '[l/kg]', 'Partition coefficient of Hg to POC', default=1000000.0_rk)
      call self%get_parameter(self%p_KdhDOC,    'p_KdhDOC',  '[l/kg]', 'Partition coefficient of Hg to DOC', default=1000000.0_rk)
      call self%get_parameter(self%p_KdmPOC6,   'p_KdmPOC6', '[l/kg]', 'Partition coefficient of MMHg to POC', default=1000000.0_rk)
      call self%get_parameter(self%p_KdmPOC8,   'p_KdmPOC8', '[l/kg]', 'Partition coefficient of MMHg to POC', default=1000000.0_rk)
      call self%get_parameter(self%p_KdmDOC,    'p_KdmDOC',  '[l/kg]', 'Partition coefficient of MMHg to DOC', default=1000000.0_rk)
      call self%get_parameter(self%p_kmet,      'p_kmet',    '[1/d]',  'Rate constant for Hg methylation', default=0.114_rk)
      call self%get_parameter(self%p_kmet2,     'p_kmet2',   '[1/d]', 'Rate constant for MMHg methylation',default=0.0016_rk)
      call self%get_parameter(self%p_kdem,      'p_kdem',    '[1/d]', 'Rate constant for MMHg demethylation',default=0.0009504_rk)
      call self%get_parameter(self%p_kphdem,    'p_kphdem',  '[1/d]', 'Rate constant for MMHg photodemethylation',default=0.007_rk)
      call self%get_parameter(self%p_kdem2,     'p_kdem2',   '[1/d]', 'Rate constant for DMHg dark demethylation', default=0.001642_rk)    
      call self%get_parameter(self%p_kpdem2,    'p_kpdem2',  '[1/d]', 'Rate constant for DMHg photodemethylation',default=0.00032832_rk)
      call self%get_parameter(self%p_kbred,     'p_kbred',   '[1/d]', 'Rate constant for biotic HgII reduction',default=0.0860_rk)
      call self%get_parameter(self%p_kbox,      'p_kbox',    '[1/d]', 'Rate constant for biotic Hg0 oxidation',default=0.140_rk )
      call self%get_parameter(self%p_kphr,      'p_kphr',    '[1/d]', 'Rate constant for photochemical HgII reduction',default=0.138_rk)    
      call self%get_parameter(self%p_kphox,     'p_kphox',   '[1/d]', 'Rate constant for photochemical Hg0 oxidation',default=0.54_rk)                                      
        call self%get_parameter(self%p_Esource,   'p_Esource',  '1-6',   'source of light for Hg processes', default=6)
  !    call self%get_parameter(self%p_use_benthic, 'p_use_benthic','0-1',  'activate first-order benthic return', default=0)
  !    call self%get_parameter(self%p_TauC,        'p_TauC',       '[d]',  'remineralization of Carbon in detritus', default=365.0_rk)
  !    call self%get_parameter(self%p_TauN,        'p_TauN',       '[d]',  'remineralization of Nitrogen in detritus', default=365.0_rk)
  !    call self%get_parameter(self%p_TauP,        'p_TauP',       '[d]',  'remineralization of Phosphorus in detritus', default=365.0_rk)
!      call self%register_state_variable(self%id_MMHg, 'MMHg', 'nmol Hg/m^3', 'Methylmercury in water', missing_value=-1.0_rk)
      ! Register links to external variables
      call self%register_state_dependency(self%id_R3c,'R3c','mgC /m^3','semi-refractory DOC')
      call self%register_state_dependency(self%id_R1c,'R1c','mgC/m^3','labile DOC')
      call self%register_state_dependency(self%id_R2c,'R2c','mgC/m^3','semilabile DOC')
      call self%register_state_dependency(self%id_R6c,'R6c','mg C/m^3','detritus-C')
      call self%register_state_dependency(self%id_R8c,'R8c','mg C/m^3','small detritus-C')
      call self%register_dependency(self%id_ruR6c, 'ruR6c', 'mmol C m-3 s-1', 'small Particulate Organic Carbon remineralization rate')
      call self%register_dependency(self%id_ruR8c, 'ruR8c', 'mmol C m-3 s-1', 'large Particulate Organic Carbon remineralization rate')

      call self%register_state_dependency(self%id_MMHg,'MMHg','nmol Hg/m^3','monomethylmercury')      
      call self%register_state_dependency(self%id_DMHg,'DMHg','nmol Hg/m^3','dimethylmercury')      
      call self%register_state_dependency(self%id_Hg2,'Hg2','nmol Hg/m^3','oxidized mercury')      
      call self%register_state_dependency(self%id_Hg0,'Hg0','nmol Hg/m^3','elemental mercury')      

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_parEIR,standard_variables%downwelling_photosynthetic_radiative_flux)  !PAOLO EIR vs PAR?
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      ! Dependency from multispectral model
  !    call self%register_dependency(self%id_PAR_tot,type_bulk_standard_variable(name='PAR_tot'))
  !    call self%register_dependency(self%id_isBen,type_bulk_standard_variable(name='isBen'))

      call self%register_diagnostic_variable(self%id_skmet,   'skmet',  'nmolHg/m3/d',  'HgII methylation',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_skmet2,  'skmet2', 'nmolHg/m3/d',  'MMHg methylation',source=source_do_column, output=output_none)
      call self%register_diagnostic_variable(self%id_skdem,   'skdem' ,   'nmolHg/m3/d', 'MMHg demethylation',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_skdem2,  'skdem2', 'nmolHg/m3/d',  'DMHg dark demethylation',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_skphdem, 'skphdem', 'nmolHg/m3/d',  'MMHg photodemethylation',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_skpdem2, 'skpdem2', 'nmolHg/m3/d',  'DMHg photodemethylation',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_skbred,  'skbred',  'nmolHg/m3/d',  'HgII Biotic Reduction',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_skphr,   'skphr',   'nmolHg/m3/d',  'HgII Photoreduction',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_skphox,  'skphox',  'nmolHg/m3/d',  'Hg0 Photooxidation',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_skbox,   'skbox','nmolHg/m3/d','Hg0 Biotic Oxidation',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_HgCl2,   'HgCl2','nmolHg/m3','Dissolved Hg',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_MMHgCl,  'MMHgCl','nmolHg/m3','Dissolved MeHg',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_HgPOC6,   'HgPOC6','nmolHg/m3','Hg in small POC  ',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_MMHgPOC6, 'MMHgPOC6','nmolHg/m3','MMHg in small POC',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_HgPOC8,   'HgPOC8','nmolHg/m3','Hg in large POC  ',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_MMHgPOC8, 'MMHgPOC8','nmolHg/m3','MMHg in large POC',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_HgDOC,   'HgDOC','nmolHg/m3','Hg in DOC',source=source_do_column,output=output_none)
      call self%register_diagnostic_variable(self%id_MMHgDOC, 'MMHgDOC','nmolHg/m3','MMHg in DOC',source=source_do_column,output=output_none)
    !  call self%register_diagnostic_variable(self%id_remX3c,    'remX3c',     'mgC/m3/d',  'remineralization of cdom X3c')
    !  call self%register_diagnostic_variable(self%id_remR3c,    'remR3c',     'mgC/m3/d',  'remineralization of dom R3c')
!     call self%register_diagnostic_dependency(self%id_flPTN6r,'flPTN6r','mmolHS/m3/d','total rate of formation of reduction equivalent') ! from PelBac

  end subroutine

   subroutine do(self,_ARGUMENTS_DO_)

      class (type_ogs_bfm_PelMercury),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

   ! !LOCAL VARIABLES:
      real(rk) :: ETW, parEIR, parEIR_wm2, tempk
      real(rk) :: R6c,R8c,R3c,R2c,R1c
      real(rk) :: Hg2, MMHg, DMHg, Hg0
      real(rk) :: remPOCmmol, ruR6c, ruR8c
      real(rk) :: skmet, skmet2,skdem,skdem2, skphdem, skpdem2
      real(rk) :: skbred, skphr, skbox, skphox
!      real(rk) :: isBen
      real(rk) :: den1,den2,faqh,faqm, fdoch,fdocm,fpoc6h,fpoc6m,fpoc8h,fpoc8m         ! [-]
      real(rk) :: p_KdhPOC6,p_KdmPOC6,p_KdhDOC,p_KdmDOC           ! [kg l-1]
      real(rk) :: p_KdhPOC8,p_KdmPOC8                          ! [kg l-1]
      real(rk) :: DOC, POC6,POC8                                 ! [mg l-1]
      real(rk) :: HgCl2, HgPOC6,HgDOC,HgPOC8     !Hg species [nmol m-3]
      real(rk) :: MMHgCl, MMHgPOC6, MMHgDOC,MMHgPOC8 !Hg species [nmol m-3]
     ! Enter spatial loops (if any)
      _LOOP_BEGIN_

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocate local memory
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
         ! Retrieve ambient nutrient concentrations

         _GET_(self%id_R3c,R3c)
         _GET_(self%id_R2c,R2c)
         _GET_(self%id_R1c,R1c)
         _GET_(self%id_R6c,R6c)
         _GET_(self%id_R8c,R8c)
         _GET_(self%id_Hg2,Hg2)
         _GET_(self%id_MMHg,MMHg)
         _GET_(self%id_DMHg,DMHg)
         _GET_(self%id_Hg0,Hg0)
  !       _GET_(self%id_O2o,O2o)
  !       _GET_(self%id_isBen,    isBen)

         ! Retrieve environmental dependencies (water temperature,
         ! photosynthetically active radation)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_ruR6c,ruR6c)
         _GET_(self%id_ruR8c,ruR8c)

         ! From where to get the light
          ! Both parEIR and PAR_tot are in uE m-2 d-1, 6=parEIR from light, 5=PAR_tot from light_spectral
!         select case (self%p_Esource)
!         case (5)
!            _GET_(self%id_parEIR,   parEIR)   ! uE m-2 d-1
!!         case (6)
            _GET_(self%id_parEIR,    parEIR)   ! uE m-2 d-1
!         end select

  tempk=ETW+273.15d0
        parEIR_wm2 = parEIR/(WtoQuanta*SEC_PER_DAY)
   DOC = (R1c +R2c+R3c)/1000.0d0 !in mg/l
   POC6=R6c/1000.0d0
   POC8=R8c/1000.0d0
   remPOCmmol = (ruR6c+ruR8c)/12.0d0     !*1000.0d0)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Partition
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  den1=1.0d0+0.000001d0*(self%p_KdhPOC6*POC6+self%p_KdhPOC8*POC8+self%p_KdhDOC*DOC)
  den2=1.0d0+0.000001d0*(self%p_KdmPOC6*POC6+self%p_KdmPOC8*POC8+self%p_KdmDOC*DOC)
  faqh=1.0d0/den1
  faqm=1.0d0/den2
  fdoch=(0.000001d0*self%p_KdhDOC*DOC)/den1
  fdocm=(0.000001d0*self%p_KdmDOC*DOC)/den2 !(1.0d0+0.000001d0*(self%p_KmPOC*POC+self%p_KmDOC*DOC))
  fpoc6h=(0.000001d0*self%p_KdhPOC6*POC6)/den1 !(1.0d0+0.000001d0*(self%p_KhPOC*POC+self%p_KhDOC*DOC))
  fpoc6m=(0.000001d0*self%p_KdmPOC6*POC6)/den2 !(1.0d0+0.000001d0*(self%p_KmPOC*POC+self%p_KmDOC*DOC))
  fpoc8h=(0.000001d0*self%p_KdhPOC8*POC8)/den1 !(1.0d0+0.000001d0*(self%p_KhPOC*POC+self%p_KhDOC*DOC))
  fpoc8m=(0.000001d0*self%p_KdmPOC8*POC8)/den2 !(1.0d0+0.000001d0*(self%p_KmPOC*POC+self%p_KmDOC*DOC))

  HgCl2=Hg2*faqh
  HgPOC6=Hg2*fpoc6h
  HgPOC8=Hg2*fpoc8h
  HgDOC=Hg2*fdoch
  MMHgCl=MMHg*faqm
  MMHgPOC6=MMHg*fpoc6m
  MMHgPOC8=MMHg*fpoc8m
  MMHgDOC=MMHg*fdocm

    _SET_DIAGNOSTIC_(self%id_HgCl2,HgCl2)
    _SET_DIAGNOSTIC_(self%id_HgPOC6,HgPOC6)
    _SET_DIAGNOSTIC_(self%id_HgPOC8,HgPOC8)
    _SET_DIAGNOSTIC_(self%id_HgDOC,HgDOC)
    _SET_DIAGNOSTIC_(self%id_MMHgCl,MMHgCl)
    _SET_DIAGNOSTIC_(self%id_MMHgPOC6,MMHgPOC6)
    _SET_DIAGNOSTIC_(self%id_MMHgPOC8,MMHgPOC8)
    _SET_DIAGNOSTIC_(self%id_MMHgDOC,MMHgDOC)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Hg Methylation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   skmet  = self%p_kmet*remPOCmmol*(HgCl2+HgDOC) !+HgPOC)
   skmet2  =self%p_kmet2*(MMHgCl+MMHgDOC)
     _SET_DIAGNOSTIC_(self%id_skmet,skmet)
     _SET_DIAGNOSTIC_(self%id_skmet2,skmet2)
     ! flux_vector
     _SET_ODE_(self%id_MMHg,skmet)
     _SET_ODE_(self%id_Hg2, -skmet)
     _SET_ODE_(self%id_DMHg,skmet2)
     _SET_ODE_(self%id_MMHg, -skmet2)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! MMHg Demethylation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  skdem   = (self%p_kdem*exp(-5500.0d0*(1.0d0/tempk-1.0d0/293.150d0)))*(MMHgCl+MMHgDOC) !+MMHgPOC)
  skphdem = self%p_kphdem*parEIR_wm2*(MMHgCl+MMHgDOC)  ! convert PAR from uE m-2 d-1 to W/m2
  _SET_DIAGNOSTIC_(self%id_skdem,skdem)
  _SET_DIAGNOSTIC_(self%id_skphdem,skphdem)
! flux_vector
 _SET_ODE_(self%id_MMHg,-skdem)
 _SET_ODE_(self%id_Hg2, skdem)

 _SET_ODE_(self%id_MMHg,-skphdem)
 _SET_ODE_(self%id_Hg2,skphdem)
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 ! DMHg Demethylation
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  skpdem2  = (self%p_kpdem2*parEIR_wm2)*DMHg
  skdem2  = ((self%p_kdem2))*DMHg
  _SET_DIAGNOSTIC_(self%id_skpdem2,skpdem2)
  _SET_DIAGNOSTIC_(self%id_skdem2,skdem2)
  ! flux_vector
  _SET_ODE_(self%id_DMHg,-skdem2)
  _SET_ODE_(self%id_MMHg,skdem2)

  _SET_ODE_(self%id_DMHg,-skpdem2)
  _SET_ODE_(self%id_MMHg,skpdem2)
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 ! Hg2 Reduction
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   skbred= self%p_kbred*remPOCmmol*(HgCl2+HgDOC)
   _SET_DIAGNOSTIC_(self%id_skbred,skbred) ! impact of denitrification on reduction equivalent
   ! flux_vector
   _SET_ODE_(self%id_Hg2, -skbred)
   _SET_ODE_(self%id_Hg0,  skbred)

  skphr= self%p_kphr*parEIR_wm2*(HgCl2+HgDOC)
  _SET_DIAGNOSTIC_(self%id_skphr,skphr)
   ! flux_vector photo
   _SET_ODE_(self%id_Hg2, -skphr)
   _SET_ODE_(self%id_Hg0,  skphr)
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   ! Hg0 Oxidation
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    skbox=self%p_kbox*remPOCmmol*Hg0
    _SET_DIAGNOSTIC_(self%id_skbox,skbox) ! impact of denitrification on reduction equivalent
    skphox= self%p_kphox*parEIR_wm2*Hg0
    _SET_DIAGNOSTIC_(self%id_skphox,skphox) ! impact of denitrification on reduction equivalent

    ! flux_vector
    _SET_ODE_(self%id_Hg2, skbox)
    _SET_ODE_(self%id_Hg0, -skbox)
    _SET_ODE_(self%id_Hg2, skphox)
    _SET_ODE_(self%id_Hg0, -skphox)

! Check unit of measure of PAR here!   parEIR is in uE m-2 d-1, p_IXn in uE m-2 s-1
  !degX1c = X1c * ( self%p_bX1c * min(parEIR/(self%p_IX1*SEC_PER_DAY),1.0_rk) ) ! Eq 13
  !degX2c = X2c * ( self%p_bX2c * min(parEIR/(self%p_IX2*SEC_PER_DAY),1.0_rk) ) ! Eq 13
  !degX3c = X3c * ( self%p_bX3c * min(parEIR/(self%p_IX3*SEC_PER_DAY),1.0_rk) ) ! Eq 13

  !_SET_DIAGNOSTIC_(self%id_degX1c,degX1c) ! photodegradation of labile CDOM
  !_SET_DIAGNOSTIC_(self%id_degX2c,degX2c) ! photodegradation of semi-labile CDOM
  !_SET_DIAGNOSTIC_(self%id_degX3c,degX3c) ! photodegradation of semi-refractory CDOM

     ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do
end module

! GP !EOC

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
