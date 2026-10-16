#include "fabm_driver.h"
module bfm_PelagicGasHg

   use fabm_types
   use ogs_bfm_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ogs_bfm_PelagicGasHg
      type (type_state_variable_id) :: id_Hg0, id_Hg2
      type (type_dependency_id) :: id_ETW, id_ESW
      type (type_horizontal_dependency_id) :: id_EWIND
      type (type_horizontal_diagnostic_variable_id) :: id_fgHg, id_wnd_diag
      real(rk) :: p_Hgin

   contains
      procedure :: initialize
      procedure :: do_surface
   end type type_ogs_bfm_PelagicGasHg

contains

   subroutine initialize(self, configunit)
      class (type_ogs_bfm_PelagicGasHg), intent(inout), target :: self
      integer, intent(in) :: configunit

      call self%get_parameter(self%p_Hgin, 'p_Hgin', 'nmol Hg/m^2/d', &
         'constant Hg2 input to the water (opposite sign to evasion)', default=0.0_rk)
      call self%register_state_dependency(self%id_Hg0, 'Hg0', 'nmol Hg/m^3', 'elemental mercury')
      call self%register_state_dependency(self%id_Hg2, 'Hg2', 'nmol Hg/m^3', 'divalent mercury')
      call self%register_dependency(self%id_ETW, standard_variables%temperature)
      call self%register_dependency(self%id_ESW, standard_variables%practical_salinity)
      call self%register_dependency(self%id_EWIND, standard_variables%wind_speed)
      call self%register_diagnostic_variable(self%id_fgHg, 'fgHg0', 'nmol Hg/m^2/d', &
         'air-sea flux of Hg0', source=source_do_surface)
      call self%register_diagnostic_variable(self%id_wnd_diag, 'wind', 'm/s', &
         'surface windspeed', source=source_do_surface)

      self%dt = 86400._rk
   end subroutine initialize

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ogs_bfm_PelagicGasHg), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      logical :: bopen1
      real(rk) :: Hg0,Hg2, ETW, ESW, EWIND, wind
      real(rk) :: Hg0atm_nmol, v1, v2, v3, v4, v5, v6, Hlw
      real(rk) :: visc, kvis, rhow, diff, tk, kw, SchHg, ScCO2
      real(rk) :: temp, salt, fgHg, salt2, temp2, temp3, temp4, temp5
      real(rk), parameter :: Hg0atm = 1.6_rk       ! ng/m3 
      real(rk), parameter :: phiw = 2.26_rk
      real(rk), parameter :: mw = 18.0_rk
      real(rk), parameter :: molHg = 12.74_rk
      real(rk), parameter :: a = 0.0001529_rk, b = 0.000016826_rk
      real(rk), parameter :: c = 0.000000083885_rk, d = 1.026681642009_rk
      real(rk), parameter :: e = 0.0024727_rk, g = 0.000048429_rk
      real(rk), parameter :: h = 0.0000047172_rk, m = 0.000000075986_rk
      real(rk), parameter :: n = 0.0000060574_rk, o = 0.000000002676_rk
      real(rk), parameter :: p = 1.013253_rk

      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_Hg0, Hg0)
         _GET_(self%id_Hg2, Hg2)
         _GET_(self%id_ETW, ETW)
         _GET_(self%id_ESW, ESW)
         _GET_HORIZONTAL_(self%id_EWIND, EWIND)

         wind = max(EWIND, 0._rk)
         temp = ETW
         salt = ESW
         tk = temp - ZERO_KELVIN

         v1 = temp*(0.06144_rk-temp*(0.001451_rk-temp*b))
         v2 = a*p
         v3 = c*d
         v4 = e*salt
         v5 = (n*p-o*d)*temp
         v6 = ((temp*g)-temp*(h-temp*m))*salt
         visc = (1.791_rk-v1-v2+v3+v4+v5+v6)/1000._rk

         temp2 = temp*temp
         temp3 = temp2*temp
         temp4 = temp3*temp
         temp5 = temp4*temp
         salt2 = salt*salt
         rhow = 999.842594_rk + 6.793952e-2_rk*temp - 9.095290e-3_rk*temp2 &
              + 1.00168e-4_rk*temp3 - 1.120083e-6_rk*temp4 + 6.536332e-9_rk*temp5 &
              + (8.24493e-1_rk - 4.0899e-3_rk*temp + 7.6438e-5_rk*temp2 &
              - 8.2467e-7_rk*temp3 + 5.3875e-9_rk*temp4)*salt &
              + (-5.72466e-3_rk + 1.0227e-4_rk*temp - 1.6546e-6_rk*temp2)*salt**1.5_rk &
              + 4.8314e-4_rk*salt2

         kvis = (visc/rhow)*10000._rk
         diff = (7.4e-8_rk*((phiw*mw)**0.50_rk)*tk)/(visc*1000._rk*(molHg**0.6_rk))
         Hlw = exp((-2403.3_rk/tk)+6.92_rk)
         ScCO2 = 0.11_rk*temp2-6.16_rk*temp+644.7_rk
         SchHg = kvis/diff

         bopen1 = .true.
         Hg0atm_nmol = Hg0atm/200.59_rk
         if (bopen1) then
            kw = (0.222_rk*wind**2+0.333_rk*wind)*(SchHg/ScCO2)**(-0.5_rk)/100._rk*24._rk
         else
            kw = (0.1_rk+2.26_rk*wind*(SchHg/ScCO2)**(-0.5_rk))/100._rk*24._rk
         end if

         fgHg = kw*(Hg0-Hg0atm_nmol/Hlw)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_wnd_diag, wind)
         _SET_SURFACE_EXCHANGE_(self%id_Hg0, -fgHg)
         _SET_SURFACE_EXCHANGE_(self%id_Hg2,  + self%p_Hgin)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fgHg, fgHg)
      _HORIZONTAL_LOOP_END_
   end subroutine do_surface

end module bfm_PelagicGasHg
