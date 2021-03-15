module ogs_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use ogs_bfm_shared
   use ogs_bfm_pelagic_base
   use bfm_Phyto 
   use ogs_bfm_light




   ! Add use statements for new models here

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
!     procedure :: initialize
      procedure :: create
   end type

   type (type_factory), save, target, public :: ogs_model_factory

contains

!     subroutine initialize(self)
!     class (type_factory), intent(inout) :: self
!     call self%register_version('ERSEM',git_commit_id//' ('//git_branch_name//' branch)')
!      end subroutine initialize

   subroutine create(self, name, model)

      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         ! Add case statements for new models here
         case ('bfm_pelagic_base'); allocate(type_ogs_bfm_pelagic_base::model)
         case ('Phyto'); allocate(type_ogs_bfm_primary_producer::model)
         case ('light'); allocate(type_ogs_bfm_light::model)

      end select

   end subroutine create

end module
