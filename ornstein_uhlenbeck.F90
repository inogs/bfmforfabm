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
      type (type_diagnostic_variable_id)   :: id_diaFLUC 
      type (type_dependency_id)            :: id_dt

      ! Parameters
      real(rk) :: D_rnd,k_rnd,time_step
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
      call self%get_parameter(self%D_rnd,  'D_rnd',   '[-]2*s',  'noise source intensity', default=0.0_rk)
      call self%get_parameter(self%k_rnd,  'k_rnd',  's-1','inverse of correlation time', default=0.0_rk)
      call self%get_parameter(self%time_step,  'time_step', 's','fabm time step', default=3600.0_rk)

! Register state variables (handled by type_bfm_pelagic_base)
      call self%initialize_bfm_base()
      call self%add_constituent('fluc',0.0_rk)  ! NB this does nothing if iron support is disabled.
      call self%add_constituent('seed',0.0_rk)

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_diaFLUC,'diaFLUC','[-]','Random Fluctuation')
              
   end subroutine

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_ogs_bfm_ornstein_uhlenbeck),intent(in) :: self
       _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: FLUC,dt,D_rnd,k_rnd
      real(rk) :: dfluc,fSEED,dSEED
      real(rk) :: OUdt
      integer  :: local_seed,local_seed0

      _HORIZONTAL_LOOP_BEGIN_

        _GET_HORIZONTAL_(self%id_fluc,FLUC)
        _GET_HORIZONTAL_(self%id_seed,fSEED)

        OUdt=self%time_step

        local_seed=int(fSEED)
        local_seed0=local_seed

        dfluc=- self%k_rnd*FLUC &
                   + 1.0_rk/OUdt*k_rnd*dsqrt(D_rnd*OUdt) *  REAL(W(local_seed,local_seed),8)

        dSEED=1.0_rk/dt*REAL(local_seed-local_seed0)

!     _SET_ODE_(self%id_FLUC,dfluc)
!     _SET_ODE_(self%id_fSEED,dSEED)

      _SET_SURFACE_EXCHANGE_(self%id_f,dfluc)
      _SET_SURFACE_EXCHANGE_(self%id_s,dSEED)


      _HORIZONTAL_LOOP_END_

   end subroutine do_surface

end module
