# 1 "shr_kind_mod.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "shr_kind_mod.F90"
!===============================================================================
! SVN $Id: shr_kind_mod.F90 41285 2012-10-26 01:46:45Z sacks $
! SVN $URL: https:
!===============================================================================

MODULE shr_kind_mod

   !----------------------------------------------------------------------------
   ! precision/kind constants add data public
   !----------------------------------------------------------------------------
   public
   integer,parameter :: SHR_KIND_R8 = selected_real_kind(12) ! 8 byte real
   integer,parameter :: SHR_KIND_R4 = selected_real_kind( 6) ! 4 byte real
   integer,parameter :: SHR_KIND_RN = kind(1.0) ! native real
   integer,parameter :: SHR_KIND_I8 = selected_int_kind (13) ! 8 byte integer
   integer,parameter :: SHR_KIND_I4 = selected_int_kind ( 6) ! 4 byte integer
   integer,parameter :: SHR_KIND_IN = kind(1) ! native integer
   integer,parameter :: SHR_KIND_CS = 80 ! short char
   integer,parameter :: SHR_KIND_CL = 256 ! long char
   integer,parameter :: SHR_KIND_CX = 512 ! extra-long char
   integer,parameter :: SHR_KIND_CXX= 4096 ! extra-extra-long char

END MODULE shr_kind_mod