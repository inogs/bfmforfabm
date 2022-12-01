#include "fabm_driver.h"

module ogs_bfm_ornstein_uhlenbeck

   use fabm_types
   use ogs_bfm_shared
   use gaussian_generator
   use ogs_bfm_pelagic_base

   implicit none

   private

   type,extends(type_ogs_bfm_pelagic_base),public :: type_ogs_bfm_ornstein_uhlenbeck
      ! Identifiers for diagnostic variables
!     type (type_state_variable_id)        :: id_FLUCf
!     type (type_state_variable_id)        :: id_FLUCs
      type (type_horizontal_diagnostic_variable_id)   :: id_diaFLUC 
      type (type_dependency_id)            :: id_dt

      ! Parameters
      real(rk) :: D_rnd,k_rnd,time_step,FLUCm
   contains
!     Model procedures
      procedure :: initialize
      procedure :: do_surface
   end type type_ogs_bfm_ornstein_uhlenbeck

contains

   subroutine initialize(self,configunit)


!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ogs_bfm_ornstein_uhlenbeck),intent(inout),target :: self
      integer,                     intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%D_rnd,  'D_rnd',   '[-]2*d',  'noise source intensity', default=0.0_rk)
      call self%get_parameter(self%k_rnd,  'k_rnd',  'd-1','inverse of correlation time', default=0.0_rk)
      call self%get_parameter(self%FLUCm,  'FLUCm',  '[-]','average fluctiation value', default=0.0_rk)
      call self%get_parameter(self%time_step,  'time_step', 's','fabm time step', default=3600.0_rk)

! Register state variables (handled by type_bfm_pelagic_base)
      call self%initialize_bfm_base() ! processes are considered in unit of days
      call self%add_constituent('fluc',0.0_rk)  ! NB this does nothing if iron support is disabled.
      call self%add_constituent('seed',0.0_rk)

      ! Register diagnostic variables
!     call self%register_horizontal_diagnostic_variable(self%id_diaFLUC,'diaFLUC','[-]','Random Fluctuation')
!     call self%register_dependency(self%id_diaFLUC,type_surface_standard_variable(name='Random_Fluctuation'))
      call self%register_diagnostic_variable(self%id_diaFLUC, 'diaFLUC',  '[-]',  'Random Fluctuation', source=source_do_surface)
              
   end subroutine

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_ogs_bfm_ornstein_uhlenbeck),intent(in) :: self
       _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: FLUC,D_rnd,k_rnd,FLUCm
      real(rk) :: dfluc,fSEED,dSEED
      real(rk) :: OUdt
      integer  :: local_seed,local_seed0

      _HORIZONTAL_LOOP_BEGIN_

        _GET_HORIZONTAL_(self%id_fluc,FLUC)
        _GET_HORIZONTAL_(self%id_seed,fSEED)

        OUdt=self%time_step

        local_seed=int(fSEED)
        local_seed0=local_seed

        D_rnd=self%D_rnd
        k_rnd=self%k_rnd
        FLUCm=self%FLUCm

!       write(*,*) 'local_seed1',local_seed

        dfluc=- k_rnd*(FLUC-FLUCm) &
                   + 86400.0_rk/OUdt*k_rnd*dsqrt(D_rnd*OUdt/86400.0_rk) *  REAL(W(local_seed,local_seed),8)

!       write(*,*) 'dfluc',dfluc
!       write(*,*) 'k_rnd',self%k_rnd
!       write(*,*) 'OUdt',OUdt
!       write(*,*) 'D_rnd',D_rnd
!       write(*,*) 'local_seed2',local_seed
!       write(*,*) '--------------------'

        dSEED=86400.0_rk/OUdt*REAL(local_seed-local_seed0,8)
!       dSEED=REAL(local_seed-local_seed0,8)

      _ADD_SURFACE_SOURCE_(self%id_fluc,dfluc)
      _ADD_SURFACE_SOURCE_(self%id_seed,dSEED)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_diaFLUC,dfluc)


      _HORIZONTAL_LOOP_END_

   end subroutine do_surface

end module
