#include "fabm_driver.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: MicroZoo
!
! DESCRIPTION
!  Microzooplankton dynamics
!
! !INTERFACE

  module bfm_MicroZoo

    use fabm_types
    use fabm_particle

    use ogs_bfm_shared
    use ogs_bfm_pelagic_base
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   use global_mem, ONLY:RLEN,ZERO,ONE,rnd_SEED
!   use constants,  ONLY:MW_C,C2ALK
! #ifdef NOPOINTERS
!   use mem
! #else
!   use mem, ONLY: D3STATE, O2o, R1c, R6c, R1n, R6n, &
!     R2c,R1p, R6p, N4n, N1p, PhytoPlankton, MicroZooPlankton, PelBacteria
!   use mem, ONLY: ppPelBacteria, ppO2o, ppR1c, ppR1l, ppR6c, ppR6s, Depth,&
!     ppR1n, ppR6n, ppR1p, ppR1l, ppR6p, ppN4n, ppN1p, ppPhytoPlankton, ppMicroZooPlankton, &
!     ETW, eO2mO2, qncPBA, qpcPBA, qncPPY, qpcPPY, qncMIZ, qpcMIZ, iiPelBacteria, &
!     qlcPPY, qscPPY, iiPhytoPlankton, iiMicroZooPlankton, iiC, iiN, iiP, iiL, iiS, &
!     NO_BOXES, iiBen, iiPel, flux_vector, quota_flux
! #ifdef INCLUDE_PELCO2
!   use mem, ONLY: ppO3c, ppO5c, ppO3h, qccPPY
! #endif
! #ifdef INCLUDE_PELFE
!   use mem, ONLY: iiF, qfcPPY, ppR6f
! #endif
! #endif
!   use mem_Param,     ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p, p_small
!   use bfm_error_msg, ONLY: bfm_error
!   use mem_MicroZoo
!   use standalone, only: delt,maxdelt
!   use random_generator
!   use gaussian_generator


!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   ! The following vector functions are used: eTq, MM, nutlim
!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   use mem_globalfun,   ONLY: eTq, MM, nutlim

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  implicit none

  private

!
! !AUTHORS
!   First ERSEM version by H. Baretta-Bekker and J.W. Baretta
!   Additional parametrizations by P. Ruardij and M. Vichi 
!   Dynamical allocation by G. Mattia 
!
! !REVISION_HISTORY
!   !
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
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  type,extends(type_ogs_bfm_pelagic_base),public :: type_ogs_bfm_microzoo
      ! NB: own state variables (c,n,p,s,f,chl) are added implicitly by deriving
      ! from type_ogs_bfm_pelagic_base!
      ! Identifiers for preys
  
      type (type_state_variable_id), allocatable,dimension(:) :: id_preyc
      type (type_state_variable_id), allocatable,dimension(:) :: id_preyn
      type (type_state_variable_id), allocatable,dimension(:) :: id_preyp
      ! type (type_state_variable_id),    allocatable,dimension(:) :: id_preys
      ! type (type_state_variable_id),    allocatable,dimension(:) :: id_preyf
      ! type (type_state_variable_id),    allocatable,dimension(:) :: id_preyl
      type (type_model_id),      allocatable,dimension(:) :: id_prey



      ! Identifiers for state variables of other models
      ! type (type_state_variable_id) :: id_O3c,id_O2o,id_O3h                  !  dissolved inorganic carbon, oxygen, total alkalinity
      type (type_state_variable_id)      :: id_O2o   !  oxygen

      ! Environmental dependencies
      type (type_dependency_id)          :: id_ETW   ! temperature

      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id) :: id_ETWd   ! temperature Celsius
      type (type_diagnostic_variable_id) :: id_et     ! temperature q10 factor
      type (type_diagnostic_variable_id) :: id_eO2    ! Oxygen limitation
      type (type_diagnostic_variable_id) :: id_rumc   ! total potential food
      type (type_diagnostic_variable_id) :: id_rugc   ! total food uptake rate (eq 38 Vichi et al. 2007)
      type (type_diagnostic_variable_id) :: id_sut    ! specific uptake rate considering potentially available food
      type (type_diagnostic_variable_id) :: id_rugn   ! tbd
      type (type_diagnostic_variable_id) :: id_rugp   ! tbd
      type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_preycd !prey c
      type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_preynd !prey n
      type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_preypd !prey p
      ! type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_preysd !prey s
      ! type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_preyfd !prey f
      ! type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_preydld !prey l

      ! Parameters (described in subroutine initialize, below)
      integer  :: nprey
      real(rk), allocatable :: p_pa(:)
      real(rk) :: p_q10, p_srs, p_sum, p_sdo, p_sd
      real(rk) :: p_pu, p_pu_ea, p_chro, p_chuc, p_minfood
      real(rk) :: p_pecaco3, p_qncMIZ, p_qpcMIZ


      ! Parameters (described in subroutine initialize, below)
  ! integer       :: i
  ! integer       :: ppzooc, ppzoon, ppzoop
  ! integer, save :: first =0
  ! integer(4)    :: local_seed
  ! integer       :: AllocStatus
  ! integer,dimension(NO_BOXES)  :: limit
  ! real(RLEN),allocatable,save,dimension(:) :: sut,et,eO2,rumc,  &
  ! rugc,rugn,rugp,runc,runn,runp, &
  ! rrsc,rrac,reac,rdc,rrtc,ruPBAc,ruPPYc,  &
  ! ruMIZc,rric,rr1c,rr6c,rr1p,rr1n, &
  ! rrip,rr6p,rep,rrin,zooc, tfluxC, tfluxN, tfluxP,MIZ_FLUCT
  ! real(RLEN),allocatable,save,dimension(:)    :: rr6n,ren,pu_ra,r
  ! real(RLEN),allocatable,save,dimension(:)    :: pe_N1p, pe_N4n, pe_R6c
  ! real(RLEN),allocatable,save,dimension(:,:)  :: PBAc,PPYc,MIZc
  ! #ifndef INCLUDE_PELCO2
  ! integer,parameter :: ppO3c = 0
  ! #endif
  contains
    ! Model procedures
    procedure :: initialize
    procedure :: do

  end type type_ogs_bfm_microzoo

contains

  subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ogs_bfm_microzoo),intent(inout),target :: self
      integer,                      intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      integer           :: iprey
      character(len=16) :: index
      real(rk) :: pippo1
!EOP
!-----------------------------------------------------------------------
!BOC
     ! Obtain the values of all model parameters from FABM.
      ! Specify the long name and units of the parameters, which could be used
      ! by FABM (or its host)
      ! to present parameters to the user for configuration (e.g., through a
      ! GUI)
      call self%get_parameter(self%p_q10,     'p_q10',     '-',         'Characteristic Q10 coefficient')
      call self%get_parameter(self%p_srs,     'p_srs',     '1/d',       'Respiration rate at 10 degrees Celsius')
      call self%get_parameter(self%p_sum,     'p_sum',     '1/d',       'Potential growth rate')
      call self%get_parameter(self%p_sdo,     'p_sdo',     '1/d',       'Mortality rate due to oxygen limitation')
      call self%get_parameter(self%p_sd,      'p_sd',      '1/d',       'Temperature independent mortality rate')
      call self%get_parameter(self%p_pu,      'p_pu',      '-',         'Assimilation efficiency')
      call self%get_parameter(self%p_pu_ea,   'p_pu_ea',   '-',         'Fraction of activity excretion')
      call self%get_parameter(self%p_chro,    'p_chro',    'mmol/m3',   'Half-saturation oxygen concentration')
      call self%get_parameter(self%p_chuc,    'p_chuc',    'mgC/m3',    'Half-saturation Food concentration for Type II')
      call self%get_parameter(self%p_minfood, 'p_minfood', 'mgC/m3',    'Half-saturation food concentration for preference factor')
      call self%get_parameter(self%p_pecaco3, 'p_pecaco3', '-',         'Portion of egested calcified shells during grazing')
      call self%get_parameter(self%p_qncMIZ,  'p_qncMIZ',  'mmolN/mgC', 'Maximum quotum P:C')
      call self%get_parameter(self%p_qpcMIZ,  'p_qpcMIZ',  'mmolN/mgC', 'Maximum quotum N:C')

      ! Register state variables (handled by type_bfm_pelagic_base)
      call self%add_constituent('c',1.e-4_rk)
      call self%add_constituent('n',1.26e-6_rk)
      call self%add_constituent('p',4.288e-8_rk)
      
      ! Determine number of prey types
      call self%get_parameter(self%nprey,'nprey','','number of prey types',default=0)
      ! Get prey-specific parameters.
      allocate(self%p_pa(self%nprey)) !Availability of nprey for microzoo group
      allocate(self%id_prey(self%nprey))
      allocate(self%id_preyc(self%nprey))
      allocate(self%id_preyn(self%nprey))
      allocate(self%id_preyp(self%nprey))
      
      do iprey=1,self%nprey
        write (index,'(i0)') iprey
        call self%get_parameter(self%p_pa(iprey),'suprey'//trim(index),'-','Availability for prey type '//trim(index))
        ! call self%register_state_dependency(self%id_preyc(iprey),'prey'//trim(index)//'','mmol C/m^3', 'prey '//trim(index)//' carbon')
        ! call self%register_state_dependency(self%id_preyc(iprey),'prey'//trim(index)//'c','mmol C/m^3', 'prey '//trim(index)//' carbon')
        ! call self%register_state_dependency(self%id_preyn(iprey),'prey'//trim(index)//'n','mmol N/m^3', 'prey '//trim(index)//' nitrogen')
        ! call self%register_state_dependency(self%id_preyp(iprey),'prey'//trim(index)//'p','mmol P/m^3', 'prey '//trim(index)//' phosphorous')
        ! call self%register_state_dependency(self%id_preys(iprey),'prey'//trim(index)//'s','mmol Si/m^3', 'prey '//trim(index)//' silicate')
        ! call self%register_state_dependency(self%id_preyl(iprey),'prey'//trim(index)//'l','mg C/m^3', 'prey '//trim(index)//' calcite')
        
        ! call self%register_diagnostic_variable(self%id_preycd(iprey),'prey'//trim(index)//'cd','mg C/m^3/d',   'ingestion of prey '//trim(index)//' carbon')
        ! call self%register_diagnostic_variable(self%id_preynd(iprey),'prey'//trim(index)//'nd','mmol N/m^3/d', 'ingestion of prey '//trim(index)//' nitrogen')
        ! call self%register_diagnostic_variable(self%id_preypd(iprey),'prey'//trim(index)//'pd','mmol P/m^3/d', 'ingestion of prey '//trim(index)//' phosphorus')
        ! call self%register_diagnostic_variable(self%id_preysd(iprey),'fprey'//trim(index)//'s','mmol Si/m^3/d','ingestion of prey '//trim(index)//' silicate')

        call self%register_state_dependency(self%id_preyc(iprey),'prey'//trim(index)//'c','mmol C/m^3', 'prey '//trim(index)//' carbon')
        call self%register_state_dependency(self%id_preyn(iprey),'prey'//trim(index)//'n','mmol n/m^3', 'prey '//trim(index)//' nitrogen')
        call self%register_state_dependency(self%id_preyp(iprey),'prey'//trim(index)//'p','mmol p/m^3', 'prey '//trim(index)//' phosphorous')

        call self%register_model_dependency(self%id_prey(iprey),'prey'//trim(index))
        call self%request_coupling_to_model(self%id_preyc(iprey),self%id_prey(iprey),'c')
        call self%request_coupling_to_model(self%id_preyn(iprey),self%id_prey(iprey),'n')
        call self%request_coupling_to_model(self%id_preyp(iprey),self%id_prey(iprey),'p')
      end do
      
      
      ! Register links to external nutrient pools.
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O2/m^3','dissolved oxygen')
      
      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      
      ! Register diagnostic variables (i.e., model outputs)
      call self%register_diagnostic_variable(self%id_ETWd, 'ETW',  'C',     'temperature Celsius')
      call self%register_diagnostic_variable(self%id_et,   'et',   '-',     'temperature factor')
      call self%register_diagnostic_variable(self%id_eO2,  'eO2',  '-',     'Oxygen limitation')
      call self%register_diagnostic_variable(self%id_rumc, 'rumc', 'mgC/m3',   'total potential food')
      call self%register_diagnostic_variable(self%id_rugc, 'rugc', 'mgC/m3/d', 'total food uptake rate')
      call self%register_diagnostic_variable(self%id_sut,  'sut',  '1/d',      'specific uptake rate')
      call self%register_diagnostic_variable(self%id_rugn, 'rugn', 'tbd',      'tbd')
      call self%register_diagnostic_variable(self%id_rugp, 'rugp', 'tbd',      'tbd')
      
      
    end subroutine

  subroutine do(self,_ARGUMENTS_DO_)
    
    class (type_ogs_bfm_microzoo),intent(in) :: self
    _DECLARE_ARGUMENTS_DO_
    
    !LOCAL VARIABLES:
    integer  :: iprey
    real(rk), dimension(self%nprey) :: preycP,preypP,preynP !,preysP, preylP
    real(rk), dimension(self%nprey) :: rupreyc
    real(rk) :: zooc, zoop, zoon
    real(rk) :: et,ETW,eO2
    real(rk) :: O2o
    real(rk) :: rumc,rugc,sut
    real(rk) :: rugn,rugp
    
    ! Enter spatial loops (if any)
    _LOOP_BEGIN_

    
    stop
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Allocate local memory
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! if (first==0) then
    !    ALLOCATE ( PBAc(NO_BOXES,iiPelBacteria),   PPYc(NO_BOXES,iiPhytoPlankton),  &
    !       &       MIZc(NO_BOXES,iiMicroZooPlankton),               &
    !       &       sut(NO_BOXES), et(NO_BOXES), eO2(NO_BOXES),      &
    !       &       rumc(NO_BOXES),  &
    !       &       rugc(NO_BOXES), rugn(NO_BOXES), rugp(NO_BOXES),  &
    !       &       runc(NO_BOXES), runn(NO_BOXES), runp(NO_BOXES),  &
    !       &       rrsc(NO_BOXES), rrac(NO_BOXES), reac(NO_BOXES),  &
    !       &       rdc(NO_BOXES) , rrtc(NO_BOXES), ruPBAc(NO_BOXES), ruPPYc(NO_BOXES), &
    !       &       ruMIZc(NO_BOXES), rric(NO_BOXES), rr1c(NO_BOXES), rr6c(NO_BOXES), &
    !       &       rr1p(NO_BOXES), rr1n(NO_BOXES), zooc(NO_BOXES), rrip(NO_BOXES), &
    !       &       rr6p(NO_BOXES), rep(NO_BOXES), rrin(NO_BOXES), rr6n(NO_BOXES), &
    !       &       ren(NO_BOXES), pu_ra(NO_BOXES), r(NO_BOXES), &
    !       &       tfluxC(NO_BOXES), tfluxN(NO_BOXES), tfluxP(NO_BOXES), &
    !       &       pe_N4n(NO_BOXES), pe_N1p(NO_BOXES), pe_R6c(NO_BOXES), MIZ_FLUCT(NO_BOXES), &
    !       &      STAT = AllocStatus )
    
    !    IF( AllocStatus /= 0 ) call bfm_error('MicroZooDynamics','Error allocating arrays')
    !    first=1
    ! end if
    
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    !  Copy  state var. object in local var
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Retrieve local biomass (carbon, phosphorus, nitrogen, chlorophyll
    ! concentrations).
    
    ! Concentrations excluding background (used in sink terms)
    ! _GET_(self%id_c,zooc)
    ! _GET_(self%id_n,zoon)
    ! _GET_(self%id_p,zoop)
    
    ! Retrieve ambient nutrient concentrations
    _GET_(self%id_O2o,O2o)
    
    ! Retrieve environmental dependencies (water temperature)
    _GET_(self%id_ETW,ETW)
    
    ! Get prey concentrations
    do iprey = 1, self%nprey
      _GET_(self%id_preyc(iprey), preycP(iprey))
      _GET_(self%id_preyn(iprey), preynP(iprey))
      _GET_(self%id_preyp(iprey), preynP(iprey))
      ! _GET_(self%idpreys(iprey), preysP(iprey))
      ! _GET_(self%idpreyl(iprey), preylP(iprey))
    enddo
    ! Prey carbon was returned in mmol (due to units of standard_variables%total_carbon); convert to mg
    ! MAYBE NO MORE NECESSARY preycP = preycP*CMass
    
    ! SKIP
    ! Quota collectors
    
    ! write(*,*) 'p_mez_sigma_rnd(zoo)', p_mez_sigma_rnd(zoo)
    ! write(*,*) 'delt ', delt
    ! write(*,*) 'dsqrt(p_mez_sigma_rnd(zoo) * delt)/delt ',
    ! dsqrt(p_mez_sigma_rnd(zoo) * delt)/delt
    ! local_seed = rnd_SEED
    ! MIZ_FLUCT(:)=ZERO
    ! do i=1,NO_BOXES
    ! SKIP      MIZ_FLUCT(i) = zooc(i) * dsqrt(p_miz_sigma_rnd(zoo)/86400.0D0 * delt)*86400.0D0/delt * W(local_seed,local_seed)
    ! enddo
    ! rnd_SEED = local_seed 
    ! TO HERE
    
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Temperature effect
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    et = eTq(ETW, self%p_q10)
    
    _SET_DIAGNOSTIC_(self%id_ETWd,ETW)
    _SET_DIAGNOSTIC_(self%id_et,et)
    
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Oxygen limitation
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    eO2 = min(ONE, MM(O2o, self%p_chro))
    
    _SET_DIAGNOSTIC_(self%id_eO2,eO2)
    
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Calculate total potential food given the non-dim prey availability
    ! and capture efficiency with loops over all LFGs.
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    rumc   = ZERO
    do iprey = 1, self%nprey
      rumc = rumc + self%p_pa(iprey)*preycP(iprey)* &
      MM(preycP(iprey), self%p_minfood)
    end do
    
    _SET_DIAGNOSTIC_(self%id_rumc,rumc)
    
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Calculate total food uptake rate (eq 38 Vichi et al. 2007) and 
    ! specific uptake rate considering potentially available food (sut)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! rugc = et*self%p_sum*MM(rumc, self%p_chuc)*zooc
    ! sut = rugc / (p_small + rumc)
    
    ! _SET_DIAGNOSTIC_(self%id_rugc,rugc)
    ! _SET_DIAGNOSTIC_(self%id_sut,sut)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Total Gross Uptakes from every LFG
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! ! Bacterioplankton
    ! rugn = ZERO
    ! rugp = ZERO
    
    _LOOP_END_
  end subroutine do
end module
    ! do iprey = 1, self%nprey
    !   rupreyc(iprey) = sut*preycP(iprey)
    !   rugn = rugn + rupreyc(iprey)*preynP(iprey)/preycP(iprey)
    !   !   do i = 1, iiPelBacteria
    !   !     ruPBAc = sut*PBAc(:,i)
    !   !     call quota_flux(iiPel, ppzooc, ppPelBacteria(i,iiC), ppzooc, ruPBAc            , tfluxC)
    !   !     call quota_flux(iiPel, ppzoon, ppPelBacteria(i,iiN), ppzoon, ruPBAc*qncPBA(i,:), tfluxN)
    !   !     call quota_flux(iiPel, ppzoop, ppPelBacteria(i,iiP), ppzoop, ruPBAc*qpcPBA(i,:), tfluxP)
    !   !     rugn = rugn + ruPBAc*qncPBA(i,:)
    !   !     rugp = rugp + ruPBAc*qpcPBA(i,:)
    !   ! SKIP MIZ_FLUCT!     call quota_flux(iiPel, ppzooc, ppPelBacteria(i,iiC), ppzooc, MIZ_FLUCT*ruPBAc/rugc, tfluxC)
    !   !     call quota_flux(iiPel, ppzoon, ppPelBacteria(i,iiN), ppzoon, MIZ_FLUCT*ruPBAc/rugc*qncPBA(i,:), tfluxN)
    !   !     call quota_flux(iiPel, ppzoop, ppPelBacteria(i,iiP), ppzoop, MIZ_FLUCT*ruPBAc/rugc*qpcPBA(i,:), tfluxP)
    !   ! to here
    ! end do
    !   ! Phytoplankton
    !   do i = 1, iiPhytoPlankton
    !     ruPPYc = sut*PPYc(:,i)
    !     call quota_flux(iiPel, ppzooc, ppPhytoPlankton(i,iiC), ppzooc, ruPPYc            , tfluxC)
    !     call quota_flux(iiPel, ppzoon, ppPhytoPlankton(i,iiN), ppzoon, ruPPYc*qncPPY(i,:), tfluxN)
    !     call quota_flux(iiPel, ppzoop, ppPhytoPlankton(i,iiP), ppzoop, ruPPYc*qpcPPY(i,:), tfluxP)
    !     rugn = rugn + ruPPYc*qncPPY(i,:)
    !     rugp = rugp + ruPPYc*qpcPPY(i,:)
    !     ! Chl is transferred to the infinite sink
    !     call flux_vector(iiPel, ppPhytoPlankton(i,iiL), &
    !                ppPhytoPlankton(i,iiL), -ruPPYc*qlcPPY(i,:))
    !     ! silicon constituent is transferred to biogenic silicate
    !     if ( ppPhytoPlankton(i,iiS) .gt. 0 ) & 
    !        call flux_vector(iiPel, ppPhytoPlankton(i,iiS), ppR6s, ruPPYc*qscPPY(i,:))
    
    ! #ifdef INCLUDE_PELFE
    !     ! Fe constituent is transferred to particulate iron
    !     if ( ppPhytoPlankton(i,iiF) .gt. 0 ) & 
!        call flux_vector(iiPel, ppPhytoPlankton(i,iiF), ppR6f, ruPPYc*qfcPPY(i,:))
! #endif

! SKIP #if defined INCLUDE_PELCO2
!     ! PIC (calcite/aragonite) production associated to the grazed biomass
!     ! The idea in PISCES is that the calcite flux exists only when associated
!     ! to a carbon release from phytoplankton (there is no calcite storage in phyto)
!     ! Use the realized rain ratio for each phytoplankton species and assume
!     ! that only a portion is egested
!     ! Calcite production is parameterized as a flux between DIC and PIC
!     ! that affects alkalinity
!     call flux_vector( iiPel, ppO3c,ppO5c, p_pecaco3(zoo)*ruPPYc*qccPPY(i,:))
!     call flux_vector( iiPel, ppO3h,ppO3h, -C2ALK*p_pecaco3(zoo)*ruPPYc*qccPPY(i,:))

!     call quota_flux(iiPel, ppzooc, ppPhytoPlankton(i,iiC), ppzooc, MIZ_FLUCT*ruPPYc/rugc, tfluxC)
!     call quota_flux(iiPel, ppzoon, ppPhytoPlankton(i,iiN), ppzoon, MIZ_FLUCT*ruPPYc/rugc*qncPPY(i,:), tfluxN)
!     call quota_flux(iiPel, ppzoop, ppPhytoPlankton(i,iiP), ppzoop, MIZ_FLUCT*ruPPYc/rugc*qpcPPY(i,:), tfluxP)
! TO HERE #endif
!   end do
!   ! Microzooplankton
!   do i = 1, iiMicroZooPlankton
!     ruMIZc = sut*MIZc(:,i)
!     ! Note that intra-group predation (cannibalism) is not added as a flux
!     if ( i/= zoo) then
!        call quota_flux(iiPel, ppzooc, ppMicroZooPlankton(i,iiC), ppzooc, ruMIZc            , tfluxC)
!        call quota_flux(iiPel, ppzoon, ppMicroZooPlankton(i,iiN), ppzoon, ruMIZc*qncMIZ(i,:), tfluxN)
!        call quota_flux(iiPel, ppzoop, ppMicroZooPlankton(i,iiP), ppzoop, ruMIZc*qpcMIZ(i,:), tfluxP)
! SKIP
!        call quota_flux(iiPel, ppzooc, ppMicroZooPlankton(i,iiC), ppzooc, MIZ_FLUCT*ruMIZc/rugc, tfluxC)
!        call quota_flux(iiPel, ppzoon, ppMicroZooPlankton(i,iiN), ppzoon, MIZ_FLUCT*ruMIZc/rugc*qncMIZ(i,:), tfluxN)
!        call quota_flux(iiPel, ppzoop, ppMicroZooPlankton(i,iiP), ppzoop, MIZ_FLUCT*ruMIZc/rugc*qpcMIZ(i,:), tfluxP)
! TO HERE
!     end if
!     rugn = rugn + ruMIZc*qncMIZ(i,:)
!     rugp = rugp + ruMIZc*qpcMIZ(i,:)
!   end do

!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!   ! Fluxes from microzooplankton
!   ! The metabolic balance is the following:
!   ! Ingestion = Growth + Excretion + Respiration
!   ! Assimilation efficiency p_pu = G/I
!   ! Excretion E = I*p_pu_ea
!   ! therefore R = (1-p_pu-p_pu_ea)*I
!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!   ! Rest, activity and total respiration fluxes
!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!   rrsc = p_srs(zoo)*et*zooc
!   ! the activity respiration is derived from the other constant parameters
!   rrac = rugc*(ONE - p_pu(zoo) - p_pu_ea(zoo))
!   rrtc = rrsc + rrac
!   call quota_flux(iiPel, ppzooc, ppzooc, ppO3c, rrtc, tfluxC)
!   call flux_vector(iiPel, ppO2o, ppO2o, -rrtc/MW_C)

!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!   ! Mortality (rdc) + Activity Excretion (reac)
!   ! and partitioning between particulate and dissolved
!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!   rdc  = ((ONE - eO2)*p_sdo(zoo) + p_sd(zoo))*zooc
!   reac = rugc*(ONE - p_pu(zoo))*p_pu_ea(zoo)
!   rric = reac + rdc
!   rr1c = rric*p_pe_R1c
!   rr6c = rric*(ONE - p_pe_R1c)
!   call quota_flux(iiPel, ppzooc, ppzooc, ppR1c, 0.98D0*rr1c, tfluxC)
!   call quota_flux(iiPel, ppzooc, ppzooc, ppR1l, 0.02D0*rr1c, tfluxC) ! To  CDOM
! ! call quota_flux(iiPel, ppzooc, ppzooc, ppR1c, rr1c, tfluxC)
!   call quota_flux(iiPel, ppzooc, ppzooc, ppR6c, rr6c, tfluxC)

!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!   !     Nutrient dynamics in microzooplankton
!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!   ! Organic Nitrogen dynamics
!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!   rrin = rugn*p_pu_ea(zoo) + rdc*qncMIZ(zoo,:)
!   rr1n = rrin*p_pe_R1n
!   rr6n = rrin - rr1n
!   call quota_flux(iiPel, ppzoon, ppzoon, ppR1n, rr1n, tfluxN)
!   call quota_flux(iiPel, ppzoon, ppzoon, ppR6n, rr6n, tfluxN)

!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!   ! Organic Phosphorus dynamics
!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!   rrip = rugp*p_pu_ea(zoo) + rdc*qpcMIZ(zoo,:)
!   rr1p = rrip*p_pe_R1p
!   rr6p = rrip - rr1p
!   call quota_flux(iiPel, ppzoop, ppzoop, ppR1p, rr1p, tfluxP)
!   call quota_flux(iiPel, ppzoop, ppzoop, ppR6p, rr6p, tfluxP)

!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!   ! Dissolved nutrient dynamics
!   ! Compare the quota of the net growth rates with the optimal quota
!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!   runc = max(ZERO, rugc*(ONE - p_pu_ea(zoo)) - rrac)
!   runn = max(ZERO, rugn*(ONE - p_pu_ea(zoo)) + rrsc*qncMIZ(zoo,:))
!   runp = max(ZERO, rugp*(ONE - p_pu_ea(zoo)) + rrsc*qpcMIZ(zoo,:))
!   ren  = max(ZERO, runn/(p_small + runc) - p_qncMIZ(zoo))* runc
!   rep  = max(ZERO, runp/(p_small + runc) - p_qpcMIZ(zoo))* runc
!   call quota_flux(iiPel, ppzoon, ppzoon, ppN4n, ren, tfluxN)
!   call quota_flux(iiPel, ppzoop, ppzoop, ppN1p, rep, tfluxP)

! SKIP FROM HERE (slow)   if ( ppzoon == 0 .or. ppzoop == 0 ) then
!      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!      ! Eliminate the excess of the non-limiting constituent under fixed quota
!      ! Determine whether C, P or N is limiting (Total Fluxes Formulation)
!      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!      limit = nutlim(tfluxc,tfluxn,tfluxp,qncMIZ(zoo,:),qpcMIZ(zoo,:),iiC,iiN,iiP)

!      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!      ! Compute the correction terms depending on the limiting constituent
!      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!      WHERE     ( limit == iiC )
!          pe_N1p = max(ZERO,tfluxp  - p_qpcMIZ(zoo)* tfluxc)
!          pe_N4n = max(ZERO,tfluxn  - p_qncMIZ(zoo)* tfluxc)
!          pe_R6c = ZERO
!      ELSEWHERE ( limit == iiP )
!          pe_N1p = ZERO
!          pe_N4n = max(ZERO, tfluxn  - tfluxp/p_qpcMIZ(zoo)*p_qncMIZ(zoo) )
!          pe_R6c = max(ZERO, tfluxc  - tfluxp/p_qpcMIZ(zoo))
!      ELSEWHERE ( limit == iiN )
!          pe_N1p = max(ZERO, tfluxp  - tfluxn/p_qncMIZ(zoo)*p_qpcMIZ(zoo))
!          pe_N4n = ZERO
!          pe_R6c = max(ZERO, tfluxc  - tfluxn/p_qncMIZ(zoo))
!      END WHERE

!      call flux_vector(iiPel, ppzooc, ppR6c, pe_R6c*(ONE-p_pe_R1c) )
!      call flux_vector(iiPel, ppzooc, ppR1c, 0.98D0*pe_R6c*(p_pe_R1c))
!      call flux_vector(iiPel, ppzooc, ppR1l, 0.02D0*pe_R6c*(p_pe_R1c)) ! To CDOM
! !    call flux_vector(iiPel, ppzooc, ppR1c, pe_R6c*(p_pe_R1c))
!      call flux_vector(iiPel, ppzoop, ppN1p, pe_N1p)
!      call flux_vector(iiPel, ppzoon, ppN4n, pe_N4n)

! #ifdef DEBUG
!      write(*,*) '+++++++++++++++'
!      if ( limit(1)==iiC ) then
!      write(*,*) 'tfluxp', tfluxp,'pe_N1p', tfluxp  - p_qpcMIZ(zoo)* tfluxc 
!      write(*,*) 'tfluxn', tfluxn,'pe_N4n', tfluxn  - p_qncMIZ(zoo)* tfluxc
!      write(*,*) 'tfluxc', tfluxc,'pe_R6c', ZERO 
!      write(*,*) 'ooooooooooooooo'
!      endif

!      if ( limit(1)==iiP ) then
!      write(*,*) 'tfluxp', tfluxp,'pe_N1p', ZERO 
!      write(*,*) 'tfluxn', tfluxn,'pe_N4n', tfluxn - tfluxp/p_qpcMIZ(zoo)*p_qncMIZ(zoo)
!      write(*,*) 'tfluxc', tfluxc,'pe_R6c', tfluxc - tfluxp/p_qpcMIZ(zoo)
!      write(*,*) 'ooooooooooooooo'
!      endif

!      if ( limit(1)==iiN ) then
!      write(*,*) 'tfluxp', tfluxp,'pe_N1p', tfluxp  - tfluxn/p_qncMIZ(zoo)*p_qpcMIZ(zoo)
!      write(*,*) 'tfluxn', tfluxn,'pe_N4n', ZERO
!      write(*,*) 'tfluxc', tfluxc,'pe_R6c', tfluxc  - tfluxn/p_qncMIZ(zoo) 
!      endif
!      write(*,*) '+++++++++++++++'
! #endif
! TO HERE
!   endif

!   end subroutine MicroZooDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
