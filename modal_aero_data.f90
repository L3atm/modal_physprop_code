# 1 "modal_aero_data.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "modal_aero_data.F90"
      module modal_aero_data

!--------------------------------------------------------------
! ... Basic aerosol mode parameters and arrays
!--------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8
      !BSINGH - Actaul value of pcnst is not required her, so I am hardwiring it to 1 for now.
      !use constituents, only: pcnst
      use radconstants, only: nswbands, nlwbands

      implicit none
      save

     !BSINGH - Actaul value of pcnst is not required her, so I am hardwiring it to 1 for now.
     integer, parameter :: pcnst = 1
     integer, parameter :: maxd_aspectype = 14
    ! aerosol mode definitions
    !

    integer, parameter :: ntot_amode = 10






    !
    ! definitions for aerosol chemical components
    !

  integer, parameter :: ntot_aspectype = 7
  character(len=*),parameter :: specname_amode(ntot_aspectype) = (/ 'sulfate   ', 'nitrate   ', &
       'p-organic ', 's-organic ', 'black-c   ', &
       'seasalt   ', 'dust      ' /)
# 46 "modal_aero_data.F90"
    ! set specdens_amode from physprop files via rad_cnst_get_aer_props
    !specdens_amode(:ntot_aspectype) = (/1770.0,1770.0,1770.0, 1000.0, 1000.0, 1700.0,1900.0,2600.0 /)

    ! rce - 06-aug-2007 - changed specmw for almost everything to match mozart

    real(r8), parameter :: specmw_amode(ntot_aspectype) = (/ 96.0_r8, 62.0_r8, &
       12.0_r8, 12.0_r8, 12.0_r8, 58.5_r8, 135.0_r8 /)
# 62 "modal_aero_data.F90"
    ! input modename_amode, nspec_amode

    character(len=*), parameter :: modename_amode(ntot_amode) = (/ &
         'accum           ', &
         'aitken          ', &
         'primary_carbon  ', &
         'fine_seasalt    ', &
         'fine_dust       ', &
         'coarse_seasalt  ', &
         'coarse7_dust    ', &
         'coarse8_dust    ', &
         'coarse9_dust    ', &
         'coarse10_dust   '/)
# 92 "modal_aero_data.F90"
    !integer, parameter :: nspec_amode(ntot_amode) = (/ 5, 3, 2, 2, 2, 2, 2, 2, 2, 2 /) ! SS
    integer, parameter :: nspec_amode(ntot_amode) = (/ 5, 3, 2, 2, 9, 2, 9, 9, 9, 9 /) ! SS





    integer, parameter :: nspec_amode_max = 9!6
    ! input mprognum_amode, mdiagnum_amode, mprogsfc_amode, mcalcwater_amode

    integer, parameter :: mprognum_amode(ntot_amode) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/)
    integer, parameter :: mdiagnum_amode(ntot_amode) = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    integer, parameter :: mprogsfc_amode(ntot_amode) = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    integer, parameter :: mcalcwater_amode(ntot_amode) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/)
# 119 "modal_aero_data.F90"
    !BSINGH - Following are hardwired in this code as the user might need to change these properties


    ! size 600507
    real(r8), parameter :: dgnum_amode(ntot_amode) = (/ 0.1100e-6_r8, 0.0260e-6_r8, 0.050e-6_r8, 0.200e-6_r8, 0.560e-6_r8, 2.000e-6_r8, 2.200e-6_r8, 3.200e-6_r8, 8.600e-6_r8, 1.400e-5_r8 /)
    real(r8), parameter :: dgnumlo_amode(ntot_amode) = (/ 0.0535e-6_r8, 0.0087e-6_r8, 0.010e-6_r8, 0.050e-6_r8, 0.360e-6_r8, 1.000e-6_r8, 1.900e-6_r8, 3.000e-6_r8, 8.000e-6_r8, 1.000e-5_r8 /)
    real(r8), parameter :: dgnumhi_amode(ntot_amode) = (/ 0.4400e-6_r8, 0.0520e-6_r8, 0.100e-6_r8, 1.000e-6_r8, 0.760e-6_r8, 4.000e-6_r8, 2.600e-6_r8, 4.400e-6_r8, 9.500e-6_r8, 1.800e-5_r8 /)

    ! size 600505
    !real(r8), parameter :: dgnum_amode(ntot_amode) = (/ 0.1100e-6_r8, 0.0260e-6_r8, 0.050e-6_r8, 0.200e-6_r8, 0.560e-6_r8, 2.000e-6_r8, 2.000e-6_r8, 3.000e-6_r8, 8.600e-6_r8, 1.400e-5_r8 /)
    !real(r8), parameter :: dgnumlo_amode(ntot_amode) = (/ 0.0535e-6_r8, 0.0087e-6_r8, 0.010e-6_r8, 0.050e-6_r8, 0.360e-6_r8, 1.000e-6_r8, 1.700e-6_r8, 2.600e-6_r8, 8.000e-6_r8, 1.000e-5_r8 /)
    !real(r8), parameter :: dgnumhi_amode(ntot_amode) = (/ 0.4400e-6_r8, 0.0520e-6_r8, 0.100e-6_r8, 1.000e-6_r8, 0.760e-6_r8, 4.000e-6_r8, 2.400e-6_r8, 3.800e-6_r8, 9.500e-6_r8, 1.800e-5_r8 /)

    ! size 600504
    !real(r8), parameter :: dgnum_amode(ntot_amode) = (/ 0.1100e-6_r8, 0.0260e-6_r8, 0.050e-6_r8, 0.200e-6_r8, 0.560e-6_r8, 2.000e-6_r8, 2.000e-6_r8, 3.000e-6_r8, 8.600e-6_r8, 1.400e-5_r8 /)
    !real(r8), parameter :: dgnumlo_amode(ntot_amode) = (/ 0.0535e-6_r8, 0.0087e-6_r8, 0.010e-6_r8, 0.050e-6_r8, 0.360e-6_r8, 1.000e-6_r8, 1.700e-6_r8, 2.600e-6_r8, 7.600e-6_r8, 1.000e-5_r8 /)
    !real(r8), parameter :: dgnumhi_amode(ntot_amode) = (/ 0.4400e-6_r8, 0.0520e-6_r8, 0.100e-6_r8, 1.000e-6_r8, 0.760e-6_r8, 4.000e-6_r8, 2.400e-6_r8, 4.000e-6_r8, 9.000e-6_r8, 1.800e-5_r8 /)

    ! size 600503
    !real(r8), parameter :: dgnum_amode(ntot_amode) = (/ 0.1100e-6_r8, 0.0260e-6_r8, 0.050e-6_r8, 0.200e-6_r8, 0.560e-6_r8, 2.000e-6_r8, 2.000e-6_r8, 3.000e-6_r8, 8.600e-6_r8, 1.400e-5_r8 /)
    !real(r8), parameter :: dgnumlo_amode(ntot_amode) = (/ 0.0535e-6_r8, 0.0087e-6_r8, 0.010e-6_r8, 0.050e-6_r8, 0.360e-6_r8, 1.000e-6_r8, 1.700e-6_r8, 2.600e-6_r8, 8.300e-6_r8, 1.000e-5_r8 /)
    !real(r8), parameter :: dgnumhi_amode(ntot_amode) = (/ 0.4400e-6_r8, 0.0520e-6_r8, 0.100e-6_r8, 1.000e-6_r8, 0.760e-6_r8, 4.000e-6_r8, 2.400e-6_r8, 3.600e-6_r8, 1.000e-5_r8, 1.800e-5_r8 /)

    ! size 600502
    !real(r8), parameter :: dgnum_amode(ntot_amode) = (/ 0.1100e-6_r8, 0.0260e-6_r8, 0.050e-6_r8, 0.200e-6_r8, 0.560e-6_r8, 2.000e-6_r8, 2.000e-6_r8, 3.600e-6_r8, 8.600e-6_r8, 1.400e-5_r8 /)
    !real(r8), parameter :: dgnumlo_amode(ntot_amode) = (/ 0.0535e-6_r8, 0.0087e-6_r8, 0.010e-6_r8, 0.050e-6_r8, 0.360e-6_r8, 1.000e-6_r8, 1.700e-6_r8, 2.800e-6_r8, 8.300e-6_r8, 1.000e-5_r8 /)
    !real(r8), parameter :: dgnumhi_amode(ntot_amode) = (/ 0.4400e-6_r8, 0.0520e-6_r8, 0.100e-6_r8, 1.000e-6_r8, 0.760e-6_r8, 4.000e-6_r8, 2.400e-6_r8, 5.000e-6_r8, 1.000e-5_r8, 1.800e-5_r8 /)

    ! size 600501
    !real(r8), parameter :: dgnum_amode(ntot_amode) = (/ 0.1100e-6_r8, 0.0260e-6_r8, 0.050e-6_r8, 0.200e-6_r8, 0.560e-6_r8, 2.000e-6_r8, 2.000e-6_r8, 5.000e-6_r8, 8.600e-6_r8, 1.400e-5_r8 /)
    !real(r8), parameter :: dgnumlo_amode(ntot_amode) = (/ 0.0535e-6_r8, 0.0087e-6_r8, 0.010e-6_r8, 0.050e-6_r8, 0.360e-6_r8, 1.000e-6_r8, 1.700e-6_r8, 3.200e-6_r8, 8.300e-6_r8, 1.000e-5_r8 /)
    !real(r8), parameter :: dgnumhi_amode(ntot_amode) = (/ 0.4400e-6_r8, 0.0520e-6_r8, 0.100e-6_r8, 1.000e-6_r8, 0.760e-6_r8, 4.000e-6_r8, 2.400e-6_r8, 6.300e-6_r8, 1.000e-5_r8, 1.800e-5_r8 /)

    ! size 60040803, 60040804
    !real(r8), parameter :: dgnum_amode(ntot_amode) = (/ 0.1100e-6_r8, 0.0260e-6_r8, 0.050e-6_r8, 0.200e-6_r8, 0.560e-6_r8, 2.000e-6_r8, 2.000e-6_r8, 5.500e-6_r8, 8.600e-6_r8, 1.400e-5_r8 /)
    !real(r8), parameter :: dgnumlo_amode(ntot_amode) = (/ 0.0535e-6_r8, 0.0087e-6_r8, 0.010e-6_r8, 0.050e-6_r8, 0.360e-6_r8, 1.000e-6_r8, 1.700e-6_r8, 4.000e-6_r8, 8.300e-6_r8, 1.000e-5_r8 /)
    !real(r8), parameter :: dgnumhi_amode(ntot_amode) = (/ 0.4400e-6_r8, 0.0520e-6_r8, 0.100e-6_r8, 1.000e-6_r8, 0.760e-6_r8, 4.000e-6_r8, 2.400e-6_r8, 6.300e-6_r8, 1.000e-5_r8, 1.800e-5_r8 /)

    ! size 60040802
    !real(r8), parameter :: dgnum_amode(ntot_amode) = (/ 0.1100e-6_r8, 0.0260e-6_r8, 0.050e-6_r8, 0.200e-6_r8, 0.560e-6_r8, 2.000e-6_r8, 2.000e-6_r8, 5.500e-6_r8, 8.600e-6_r8, 1.800e-5_r8 /)
    !real(r8), parameter :: dgnumlo_amode(ntot_amode) = (/ 0.0535e-6_r8, 0.0087e-6_r8, 0.010e-6_r8, 0.050e-6_r8, 0.360e-6_r8, 1.000e-6_r8, 1.700e-6_r8, 4.000e-6_r8, 7.800e-6_r8, 1.200e-5_r8 /)
    !real(r8), parameter :: dgnumhi_amode(ntot_amode) = (/ 0.4400e-6_r8, 0.0520e-6_r8, 0.100e-6_r8, 1.000e-6_r8, 0.760e-6_r8, 4.000e-6_r8, 2.400e-6_r8, 6.300e-6_r8, 9.600e-6_r8, 2.500e-5_r8 /)
    ! size 600407
    !real(r8), parameter :: dgnum_amode(ntot_amode) = (/ 0.1100e-6_r8, 0.0260e-6_r8, 0.050e-6_r8, 0.200e-6_r8, 0.560e-6_r8, 2.000e-6_r8, 2.000e-6_r8, 5.800e-6_r8, 1.200e-5_r8, 2.500e-5_r8 /)
    !real(r8), parameter :: dgnumlo_amode(ntot_amode) = (/ 0.0535e-6_r8, 0.0087e-6_r8, 0.010e-6_r8, 0.050e-6_r8, 0.360e-6_r8, 1.000e-6_r8, 1.700e-6_r8, 5.000e-6_r8, 9.200e-6_r8, 2.000e-5_r8 /)
    !real(r8), parameter :: dgnumhi_amode(ntot_amode) = (/ 0.4400e-6_r8, 0.0520e-6_r8, 0.100e-6_r8, 1.000e-6_r8, 0.760e-6_r8, 4.000e-6_r8, 2.400e-6_r8, 6.600e-6_r8, 1.600e-5_r8, 3.500e-5_r8 /)
    ! size 6004, 600406
    !real(r8), parameter :: dgnum_amode(ntot_amode) = (/ 0.1100e-6_r8, 0.0260e-6_r8, 0.050e-6_r8, 0.200e-6_r8, 0.560e-6_r8, 2.000e-6_r8, 2.000e-6_r8, 5.800e-6_r8, 1.200e-6_r8, 3.000e-5_r8 /)
    !real(r8), parameter :: dgnumlo_amode(ntot_amode) = (/ 0.0535e-6_r8, 0.0087e-6_r8, 0.010e-6_r8, 0.050e-6_r8, 0.360e-6_r8, 1.000e-6_r8, 1.700e-6_r8, 5.000e-6_r8, 9.600e-6_r8, 2.500e-5_r8 /)
    !real(r8), parameter :: dgnumhi_amode(ntot_amode) = (/ 0.4400e-6_r8, 0.0520e-6_r8, 0.100e-6_r8, 1.000e-6_r8, 0.760e-6_r8, 4.000e-6_r8, 2.400e-6_r8, 6.600e-6_r8, 1.800e-5_r8, 4.000e-5_r8 /)

    ! size 6003
    !real(r8), parameter :: dgnum_amode(ntot_amode) = (/ 0.1100e-6_r8, 0.0260e-6_r8, 0.050e-6_r8, 0.200e-6_r8, 0.110e-6_r8, 2.000e-6_r8, 2.000e-6_r8, 6.000e-6_r8, 8.200e-6_r8, 1.200e-5_r8 /)
    !real(r8), parameter :: dgnumlo_amode(ntot_amode) = (/ 0.0535e-6_r8, 0.0087e-6_r8, 0.010e-6_r8, 0.050e-6_r8, 0.054e-6_r8, 1.000e-6_r8, 1.700e-6_r8, 5.400e-6_r8, 7.200e-6_r8, 1.000e-5_r8 /)
    !real(r8), parameter :: dgnumhi_amode(ntot_amode) = (/ 0.4400e-6_r8, 0.0520e-6_r8, 0.100e-6_r8, 1.000e-6_r8, 0.440e-6_r8, 4.000e-6_r8, 3.300e-6_r8, 7.000e-6_r8, 9.500e-6_r8, 2.000e-5_r8 /)
    ! size 6002
    !real(r8), parameter :: dgnum_amode(ntot_amode) = (/ 0.1100e-6_r8, 0.0260e-6_r8, 0.050e-6_r8, 0.200e-6_r8, 0.110e-6_r8, 2.000e-6_r8, 1.500e-6_r8, 4.800e-6_r8, 7.000e-6_r8, 1.000e-5_r8 /)
    !real(r8), parameter :: dgnumlo_amode(ntot_amode) = (/ 0.0535e-6_r8, 0.0087e-6_r8, 0.010e-6_r8, 0.050e-6_r8, 0.054e-6_r8, 1.000e-6_r8, 1.000e-6_r8, 3.000e-6_r8, 6.200e-6_r8, 9.000e-6_r8 /)
    !real(r8), parameter :: dgnumhi_amode(ntot_amode) = (/ 0.4400e-6_r8, 0.0520e-6_r8, 0.100e-6_r8, 1.000e-6_r8, 0.440e-6_r8, 4.000e-6_r8, 2.000e-6_r8, 5.800e-6_r8, 8.000e-6_r8, 1.600e-5_r8 /)

    ! size 4903,4904, 60, 6001
    !real(r8), parameter :: dgnum_amode(ntot_amode) = (/ 0.1100e-6_r8, 0.0260e-6_r8, 0.050e-6_r8, 0.200e-6_r8, 0.560e-6_r8, 2.000e-6_r8, 2.000e-6_r8, 5.800e-6_r8, 7.600e-6_r8, 1.200e-5_r8 /)
    !real(r8), parameter :: dgnumlo_amode(ntot_amode) = (/ 0.0535e-6_r8, 0.0087e-6_r8, 0.010e-6_r8, 0.050e-6_r8, 0.360e-6_r8, 1.000e-6_r8, 1.700e-6_r8, 4.800e-6_r8, 6.600e-6_r8, 9.600e-6_r8 /)
    !real(r8), parameter :: dgnumhi_amode(ntot_amode) = (/ 0.4400e-6_r8, 0.0520e-6_r8, 0.100e-6_r8, 1.000e-6_r8, 0.760e-6_r8, 4.000e-6_r8, 2.400e-6_r8, 6.200e-6_r8, 8.600e-6_r8, 1.800e-5_r8 /)

    ! size 4902
    !real(r8), parameter :: dgnum_amode(ntot_amode) = (/ 0.1100e-6_r8, 0.0260e-6_r8, 0.050e-6_r8, 0.200e-6_r8, 0.560e-6_r8, 2.000e-6_r8, 2.000e-6_r8, 6.000e-6_r8, 8.200e-6_r8, 1.500e-5_r8 /)
    !real(r8), parameter :: dgnumlo_amode(ntot_amode) = (/ 0.0535e-6_r8, 0.0087e-6_r8, 0.010e-6_r8, 0.050e-6_r8, 0.360e-6_r8, 1.000e-6_r8, 1.700e-6_r8, 5.000e-6_r8, 7.000e-6_r8, 1.000e-5_r8 /)
    !real(r8), parameter :: dgnumhi_amode(ntot_amode) = (/ 0.4400e-6_r8, 0.0520e-6_r8, 0.100e-6_r8, 1.000e-6_r8, 0.760e-6_r8, 4.000e-6_r8, 2.400e-6_r8, 6.500e-6_r8, 9.000e-6_r8, 2.000e-5_r8 /)

    ! size 49
    !real(r8), parameter :: dgnum_amode(ntot_amode) = (/ 0.1100e-6_r8, 0.0260e-6_r8, 0.050e-6_r8, 0.200e-6_r8, 0.560e-6_r8, 2.000e-6_r8, 2.000e-6_r8, 6.000e-6_r8, 1.200e-5_r8, 2.500e-5_r8 /)
    !real(r8), parameter :: dgnumlo_amode(ntot_amode) = (/ 0.0535e-6_r8, 0.0087e-6_r8, 0.010e-6_r8, 0.050e-6_r8, 0.360e-6_r8, 1.000e-6_r8, 1.700e-6_r8, 5.400e-6_r8, 1.000e-5_r8, 2.000e-5_r8 /)
    !real(r8), parameter :: dgnumhi_amode(ntot_amode) = (/ 0.4400e-6_r8, 0.0520e-6_r8, 0.100e-6_r8, 1.000e-6_r8, 0.760e-6_r8, 4.000e-6_r8, 2.400e-6_r8, 7.000e-6_r8, 1.500e-5_r8, 3.000e-5_r8 /)
# 203 "modal_aero_data.F90"
    ! input sigmag_amode, sigmaglo_amode, sigmaghi_amode


   ! size 60040803, 60040804
   real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 1.800_r8, 2.000_r8, 1.700_r8, 1.500_r8, 1.5000_r8, 1.500_r8 /)
   ! size 60040802
   !real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 1.550_r8, 2.000_r8, 1.550_r8, 1.500_r8, 1.400_r8, 1.500_r8 /)
   ! size 600408
   !real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 1.400_r8, 2.000_r8, 1.400_r8, 1.500_r8, 1.500_r8, 1.500_r8 /)
   ! size 600407
   !real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 1.400_r8, 2.000_r8, 1.400_r8, 1.400_r8, 1.400_r8, 1.400_r8 /)
   ! size 6004, 600406
   !real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 1.300_r8, 2.000_r8, 1.300_r8, 1.200_r8, 1.200_r8, 1.200_r8 /)
   ! size 6001
   !real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 1.300_r8, 2.000_r8, 1.300_r8, 1.300_r8, 1.300_r8, 1.300_r8 /)
   ! size 60
   !real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 1.500_r8, 2.000_r8, 1.500_r8, 1.500_r8, 1.500_r8, 1.500_r8 /)
   ! size 4905
   !real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 2.000_r8, 2.000_r8, 2.000_r8, 2.100_r8, 2.200_r8, 2.200_r8 /)
   ! size 4904
   !real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 2.300_r8, 2.000_r8, 2.000_r8, 2.100_r8, 1.800_r8, 1.600_r8 /)
   ! size 4903
   !real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 2.300_r8, 2.000_r8, 2.000_r8, 2.100_r8, 2.300_r8, 2.300_r8 /)
   !real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 2.300_r8, 2.000_r8, 2.000_r8, 2.100_r8, 2.600_r8, 2.600_r8 /)
   !real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 2.300_r8, 2.000_r8, 2.000_r8, 2.300_r8, 2.600_r8, 2.600_r8 /)
   !real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 2.600_r8, 2.000_r8, 2.300_r8, 2.300_r8, 2.600_r8, 2.600_r8 /)
   !real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 2.600_r8, 2.000_r8, 2.600_r8, 2.600_r8, 2.600_r8, 2.600_r8 /)
   !real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 2.000_r8, 2.000_r8, 2.100_r8, 1.800_r8, 2.100_r8, 2.200_r8 /)
   ! real(r8), parameter :: sigmag_amode(ntot_amode) = (/ 1.800_r8, 1.600_r8, 1.600_r8, 2.000_r8, 1.300_r8, 2.000_r8, 1.400_r8, 1.600_r8, 2.500_r8, 2.500_r8 /)






    ! input crystalization and deliquescence points

    real(r8), parameter :: rhcrystal_amode(ntot_amode) = (/ 0.350_r8, 0.350_r8, 0.350_r8, 0.350_r8, 0.350_r8, 0.350_r8, 0.350_r8, 0.350_r8, 0.350_r8, 0.350_r8 /)
    real(r8), parameter :: rhdeliques_amode(ntot_amode) = (/ 0.800_r8, 0.800_r8, 0.800_r8, 0.800_r8, 0.800_r8, 0.800_r8, 0.800_r8, 0.800_r8, 0.800_r8, 0.800_r8 /)







    !BSINGH -ENDS

    !BSINGH - Following is commented out as they are declared above as parameters
    ! input dgnum_amode, dgnumlo_amode, dgnumhi_amode (units = m)
    !real(r8) :: dgnum_amode(ntot_amode)
    !real(r8) :: dgnumlo_amode(ntot_amode)
    !real(r8) :: dgnumhi_amode(ntot_amode)

    ! input sigmag_amode
    !real(r8) :: sigmag_amode(ntot_amode)

    ! input crystalization and deliquescence points
    !real(r8) :: rhcrystal_amode(ntot_amode)
    !real(r8) :: rhdeliques_amode(ntot_amode)
    !BSINGH - ENDS


    integer :: msectional = -1


      integer & !
          lspectype_amode( maxd_aspectype, ntot_amode ), & !
          lmassptr_amode( maxd_aspectype, ntot_amode ), & !
          lmassptrcw_amode( maxd_aspectype, ntot_amode ), & !
          numptr_amode( ntot_amode ), & !
          numptrcw_amode( ntot_amode )


      real(r8) :: & !
          alnsg_amode( ntot_amode ), & !
          voltonumb_amode( ntot_amode ), & !
          voltonumblo_amode( ntot_amode ), & !
          voltonumbhi_amode( ntot_amode ), & !
          alnv2n_amode( ntot_amode ), & !
          alnv2nlo_amode( ntot_amode ), & !
          alnv2nhi_amode( ntot_amode ), & !
          specdens_amode( maxd_aspectype ), & !
          spechygro( maxd_aspectype )


      complex(r8) & !
          specrefndxsw( nswbands, maxd_aspectype ), & !
          specrefndxlw( nlwbands, maxd_aspectype )


      character(len=16) :: cnst_name_cw( pcnst )

      character(len=8) :: aodvisname(ntot_amode ), &
                          ssavisname(ntot_amode )
      character(len=48) :: aodvislongname(ntot_amode ), &
                           ssavislongname(ntot_amode )

      character(len=8) :: fnactname(ntot_amode ), &
                          fmactname(ntot_amode ), &
                          nactname(ntot_amode )
      character(len=48) :: fnactlongname(ntot_amode ), &
                           fmactlongname(ntot_amode ), &
                           nactlongname(ntot_amode )

      integer & !
          lptr_so4_a_amode(ntot_amode), lptr_so4_cw_amode(ntot_amode), & !
          lptr_msa_a_amode(ntot_amode), lptr_msa_cw_amode(ntot_amode), & !
          lptr_nh4_a_amode(ntot_amode), lptr_nh4_cw_amode(ntot_amode), & !
          lptr_no3_a_amode(ntot_amode), lptr_no3_cw_amode(ntot_amode), & !
          lptr_pom_a_amode(ntot_amode), lptr_pom_cw_amode(ntot_amode), & !
          lptr_soa_a_amode(ntot_amode), lptr_soa_cw_amode(ntot_amode), & !
          lptr_bc_a_amode(ntot_amode), lptr_bc_cw_amode(ntot_amode), & !
          lptr_nacl_a_amode(ntot_amode), lptr_nacl_cw_amode(ntot_amode),& !
          lptr_dust_a_amode(ntot_amode), lptr_dust_cw_amode(ntot_amode),& !
          modeptr_accum, modeptr_aitken, & !
          modeptr_ufine, modeptr_coarse, & !
          modeptr_pcarbon, & !
          modeptr_finedust, modeptr_fineseas, & !
          modeptr_coardust, modeptr_coarseas

      real(r8) :: &
          specmw_so4_amode, specdens_so4_amode, &
          specmw_nh4_amode, specdens_nh4_amode, &
          specmw_no3_amode, specdens_no3_amode, &
          specmw_pom_amode, specdens_pom_amode, &
          specmw_soa_amode, specdens_soa_amode, &
          specmw_bc_amode, specdens_bc_amode, &
          specmw_dust_amode, specdens_dust_amode, &
          specmw_seasalt_amode, specdens_seasalt_amode

 integer species_class(pcnst) ! indicates species class (
    ! cldphysics, aerosol, gas )

 integer spec_class_undefined
 parameter ( spec_class_undefined = 0 )
 integer spec_class_cldphysics
 parameter ( spec_class_cldphysics = 1 )
 integer spec_class_aerosol
 parameter ( spec_class_aerosol = 2 )
 integer spec_class_gas
 parameter ( spec_class_gas = 3 )
 integer spec_class_other
 parameter ( spec_class_other = 4 )


! threshold for reporting negatives from subr qneg3
      real(r8) :: qneg3_worst_thresh_amode(pcnst)

      integer, private :: qqcw(pcnst)=-1 ! Remaps modal_aero indices into pbuf
# 400 "modal_aero_data.F90"
      end module modal_aero_data

!----------------------------------------------------------------
!
! maxd_aspectype = maximum allowable number of chemical species
! in each aerosol mode
!
! ntot_amode = number of aerosol modes
! ( ntot_amode_gchm = number of aerosol modes in gchm
! ntot_amode_ccm2 = number of aerosol modes to be made known to ccm2
! These are temporary until multi-mode activation scavenging is going.
! Until then, ntot_amode is set to either ntot_amode_gchm or
! ntot_amode_ccm2 depending on which code is active )
!
! msectional - if positive, moving-center sectional code is utilized,
! and each mode is actually a section.
! msectional_concinit - if positive, special code is used to initialize
! the mixing ratios of all the sections.
!
! nspec_amode(m) = number of chemical species in aerosol mode m
! nspec_amode_ccm2(m) = . . . while in ccm2 code
! nspec_amode_gchm(m) = . . . while in gchm code
! nspec_amode_nontracer(m) = number of "non-tracer" chemical
! species while in gchm code
! lspectype_amode(l,m) = species type/i.d. for chemical species l
! in aerosol mode m. (1=sulfate, others to be defined)
! lmassptr_amode(l,m) = gchm r-array index for the mixing ratio
! (moles-x/mole-air) for chemical species l in aerosol mode m
! that is in clear air or interstitial air (but not in cloud water)
! lmassptrcw_amode(l,m) = gchm r-array index for the mixing ratio
! (moles-x/mole-air) for chemical species l in aerosol mode m
! that is currently bound/dissolved in cloud water
! lwaterptr_amode(m) = gchm r-array index for the mixing ratio
! (moles-water/mole-air) for water associated with aerosol mode m
! that is in clear air or interstitial air
! lkohlercptr_amode(m) = gchm r-array index for the kohler "c" parameter
! for aerosol mode m. This is defined on a per-dry-particle-mass basis:
! c = r(i,j,k,lkohlercptr_amode) * [rhodry * (4*pi/3) * rdry^3]
! numptr_amode(m) = gchm r-array index for the number mixing ratio
! (particles/mole-air) for aerosol mode m that is in clear air or
! interstitial are (but not in cloud water). If zero or negative,
! then number is not being simulated.
! ( numptr_amode_gchm(m) = same thing but for within gchm
! numptr_amode_ccm2(m) = same thing but for within ccm2
! These are temporary, to allow testing number in gchm before ccm2 )
! numptrcw_amode(m) = gchm r-array index for the number mixing ratio
! (particles/mole-air) for aerosol mode m
! that is currently bound/dissolved in cloud water
! lsfcptr_amode(m) = gchm r-array index for the surface area mixing ratio
! (cm^2/mole-air) for aerosol mode m that is in clear air or
! interstitial are (but not in cloud water). If zero or negative,
! then surface area is not being simulated.
! lsfcptrcw_amode(m) = gchm r-array index for the surface area mixing ratio
! (cm^2/mole-air) for aerosol mode m that is currently
! bound/dissolved in cloud water.
! lsigptr_amode(m) = gchm r-array index for sigmag for aerosol mode m
! that is in clear air or interstitial are (but not in cloud water).
! If zero or negative, then the constant sigmag_amode(m) is used.
! lsigptrcw_amode(m) = gchm r-array index for sigmag for aerosol mode m
! that is currently bound/dissolved in cloud water.
! If zero or negative, then the constant sigmag_amode(m) is used.
! lsigptrac_amode(m) = gchm r-array index for sigmag for aerosol mode m
! for combined clear-air/interstial plus bound/dissolved in cloud water.
! If zero or negative, then the constant sigmag_amode(m) is used.
!
! dgnum_amode(m) = geometric dry mean diameter (m) of the number
! distribution for aerosol mode m.
! (Only used when numptr_amode(m) is zero or negative.)
! dgnumlo_amode(m), dgnumhi_amode(m) = lower and upper limits on the
! geometric dry mean diameter (m) of the number distribution
! (Used when mprognum_amode>0, to limit dgnum to reasonable values)
! sigmag_amode(m) = geometric standard deviation for aerosol mode m
! sigmaglo_amode(m), sigmaghi_amode(m) = lower and upper limits on the
! geometric standard deviation of the number distribution
! (Used when mprogsfc_amode>0, to limit sigmag to reasonable values)
! alnsg_amode(m) = alog( sigmag_amode(m) )
! alnsglo_amode(m), alnsghi_amode(m) = alog( sigmaglo/hi_amode(m) )
! voltonumb_amode(m) = ratio of number to volume for mode m
! voltonumblo_amode(m), voltonumbhi_amode(m) = ratio of number to volume
! when dgnum = dgnumlo_amode or dgnumhi_amode, respectively
! voltosfc_amode(m), voltosfclo_amode(m), voltosfchi_amode(m) - ratio of
! surface to volume for mode m (like the voltonumb_amode s)
! alnv2n_amode(m), alnv2nlo_amode(m), alnv2nhi_amode(m) -
! alnv2n_amode(m) = alog( voltonumblo_amode(m) ), ...
! alnv2s_amode(m), alnv2slo_amode(m), alnv2shi_amode(m) -
! alnv2s_amode(m) = alog( voltosfclo_amode(m) ), ...
! rhcrystal_amode(m) = crystalization r.h. for mode m
! rhdeliques_amode(m) = deliquescence r.h. for mode m
! (*** these r.h. values are 0-1 fractions, not 0-100 percentages)
!
! mcalcwater_amode(m) - if positive, water content for mode m will be
! calculated and stored in rclm(k,lwaterptr_amode(m)). Otherwise, no.
! mprognum_amode(m) - if positive, number mixing-ratio for mode m will
! be prognosed. Otherwise, no.
! mdiagnum_amode(m) - if positive, number mixing-ratio for mode m will
! be diagnosed and put into rclm(k,numptr_amode(m)). Otherwise, no.
! mprogsfc_amode(m) - if positive, surface area mixing-ratio for mode m will
! be prognosed, and sigmag will vary temporally and spatially.
! Otherwise, sigmag is constant.
! *** currently surface area is not prognosed when msectional>0 ***
!
! ntot_aspectype = overall number of aerosol chemical species defined (over all modes)
! specdens_amode(l) = dry density (kg/m^3) of aerosol chemical species type l
! specmw_amode(l) = molecular weight (kg/kmol) of aerosol chemical species type l
! specname_amode(l) = name of aerosol chemical species type l
! specrefndxsw(l) = complex refractive index (visible wavelengths)
! of aerosol chemical species type l
! specrefndxlw(l) = complex refractive index (infrared wavelengths)
! of aerosol chemical species type l
! spechygro(l) = hygroscopicity of aerosol chemical species type l
!
! lptr_so4_a_amode(m), lptr_so4_cw_amode(m) = gchm r-array index for the
! mixing ratio for sulfate associated with aerosol mode m
! ("a" and "cw" phases)
! (similar for msa, oc, bc, nacl, dust)
!
! modename_amode(m) = character-variable name for mode m,
! read from mirage2.inp
! modeptr_accum - mode index for the main accumulation mode
! if modeptr_accum = 1, then mode 1 is the main accumulation mode,
! and modename_amode(1) = "accum"
! modeptr_aitken - mode index for the main aitken mode
! if modeptr_aitken = 2, then mode 2 is the main aitken mode,
! and modename_amode(2) = "aitken"
! modeptr_ufine - mode index for the ultrafine mode
! if modeptr_ufine = 3, then mode 3 is the ultrafine mode,
! and modename_amode(3) = "ufine"
! modeptr_coarseas - mode index for the coarse sea-salt mode
! if modeptr_coarseas = 4, then mode 4 is the coarse sea-salt mode,
! and modename_amode(4) = "coarse_seasalt"
! modeptr_coardust - mode index for the coarse dust mode
! if modeptr_coardust = 5, then mode 5 is the coarse dust mode,
! and modename_amode(5) = "coarse_dust"
!
! specdens_XX_amode = dry density (kg/m^3) of aerosol chemical species type XX
! where XX is so4, om, bc, dust, seasalt
! contains same values as the specdens_amode array
! allows values to be referenced differently
! specmw_XX_amode = molecular weight (kg/kmol) of aerosol chemical species type XX
! contains same values as the specmw_amode array
!
!-----------------------------------------------------------------------


!--------------------------------------------------------------
!
! ... aerosol size information for the current chunk
!
!--------------------------------------------------------------
!
! dgncur = current geometric mean diameters (cm) for number distributions
! dgncur_a - for unactivated particles, dry
! (in physics buffer as DGNUM)
! dgncur_awet - for unactivated particles, wet at grid-cell ambient RH
! (in physics buffer as DGNUMWET)
!
! the dgncur are computed from current mass and number
! mixing ratios in the grid cell, BUT are then adjusted to be within
! the bounds defined by dgnumlo/hi_amode
!
! v2ncur = current (number/volume) ratio based on dgncur and sgcur
! (volume in cm^3/whatever, number in particles/whatever)
! == 1.0 / ( pi/6 * dgncur**3 * exp(4.5*((log(sgcur))**2)) )
! v2ncur_a - for unactivated particles
! (currently just defined locally)
!
