# 1 "module_optics_support.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "module_optics_support.F90"
module module_optics_support

  use netcdf, only: nf90_open, NF90_NOWRITE, nf90_inq_varid, &
       nf90_get_var, nf90_noerr, nf90_strerror, nf90_inq_dimid, &
       nf90_inquire_dimension, nf90_close
  use shr_kind_mod, only: r8=>shr_kind_r8

  implicit none

  logical, parameter :: masterproc = .true.
  integer, parameter :: iulog = 6, cs1 = 256
  integer, parameter :: ncoef=5, prefr=7, prefi=10
  integer:: nswb, nlwb

  !local
  complex(r8), pointer :: refidx_aer_lw(:,:,:), refidx_aer_sw(:,:,:)

  contains

    subroutine endrun(msg)
      character(len=*), intent(in), optional :: msg

      write(iulog,*) msg
      !fortran stop
      stop
    end subroutine endrun

    subroutine read_modal_optics(nmodes, maxd_aspectype, nfiles, fnames, nswbands, nlwbands)

      implicit none

      integer, intent(in) :: nmodes, nfiles, nswbands, nlwbands, maxd_aspectype
      character(len=*), dimension(:), intent(in) :: fnames

      !local vars
      character(len=256) :: file, spec, tmpstr1, init_str, adv_str, rest_str
      character(len=256), dimension(nmodes,maxd_aspectype):: spec_in_mode
      integer :: ifp, nc_id, vid_im, vid_real, iw, pos, imode, iloop

      integer, dimension(nmodes) :: ispec !# of species in a particular mode

      real(r8), pointer :: ref_real(:), ref_im(:) ! tmp storage for components of complex index

      !Assign short and long wave band #s
      nswb = nswbands
      nlwb = nlwbands

      ispec(:) = 0
      allocate(refidx_aer_sw(nmodes,maxd_aspectype,nswbands))
      allocate(refidx_aer_lw(nmodes,maxd_aspectype,nlwbands))

      ifp = 0 !File counter
      iloop = 0 !Loop counter
      do while(.true.) !Loop should stop when ifp==nfiles
         iloop = iloop + 1
         !Following section of the code seprates out species name, mode number
         !and path to physprop files for that species. I am using 'index' function
         !to split the strings, which is rather cumbersome for this simple text
         !manupulation. This should be replaced with a better method in the future.

         init_str = fnames(iloop) !Store possible species+file name in a string
         adv_str = init_str(1:1) !adv_str should be 'A'[advected] for it to have a species and a physprop file path

         if(trim(adjustl(adv_str)) .eq. 'A') then !Extract species name and physprop file

            pos = index(init_str,':N:') !Find the position of first occurance of ":N:"
            rest_str = init_str(pos+1:) !Store rest of the string after ":N:"

            spec = init_str(3:pos-1) !This is the species name

            if(trim(adjustl(spec(1:3))) .ne. 'num') then !aerosol numbers (num_a1, num_a2 etc.) do not have physprop files
               !Find mode # of this species
               !write(tmpstr1,*)init_str(pos-1:pos-1)
               !read(tmpstr1,*)imode !This is the mode #

               !if(init_str(7:pos-1) .eq. 'a10') then
               if(init_str(8:pos-1) .eq. 'a10' .or. init_str(7:pos-1) .eq. 'a10') then
                  write(tmpstr1,*)init_str(pos-2:pos-1)
                  read(tmpstr1,*)imode !This is the mode #
               else
                  write(tmpstr1,*)init_str(pos-1:pos-1)
                  read(tmpstr1,*)imode !This is the mode #
               endif

               !Now find the file path
               pos = index(rest_str,':/') !Find the position of ":/"
               rest_str = rest_str(pos+1:)
               !Get rid of trailing ":+" if it exists
               pos = index(rest_str,':')
               if(pos == 0) then
                  file = rest_str
               else
                  file = rest_str(:pos-1)
               endif
               ifp = ifp + 1 !Increment file count
            else
               cycle ! For skipping aerosol numbers (num_a1, num_a2 etc.)
            endif
         else
            cycle !cycle as this string does not contain species name and physprop file
         endif

         ispec(imode) = ispec(imode) + 1 !Tracks how many species are there in each mode
         spec_in_mode(imode,ispec(imode)) = trim(spec) !Store species name as a function of mode and # of species in that mode

         !open file
         call check( nf90_open(trim(adjustl(file)), NF90_NOWRITE, nc_id) )

         !=== shortwave ===
         call check(nf90_inq_varid(nc_id, 'refindex_real_aer_sw', vid_real)) !real
         call check(nf90_inq_varid(nc_id, 'refindex_im_aer_sw', vid_im)) !imaginary

         !temporary allocation
         allocate(ref_real(nswbands), ref_im(nswbands))

         call check(nf90_get_var(nc_id, vid_real, ref_real),': ERROR reading refindex_real_aer_sw')
         call check(nf90_get_var(nc_id, vid_im, ref_im),': ERROR reading refindex_im_aer_sw')


         !successfully read refindex data -- set complex values in local vars
         do iw = 1, nswbands
            refidx_aer_sw(imode,ispec(imode),iw) = cmplx(ref_real(iw), abs(ref_im(iw)))
         end do

         !deallocate local vars
         deallocate(ref_real, ref_im)

         !=== longwave ===
         call check(nf90_inq_varid(nc_id, 'refindex_real_aer_lw', vid_real))
         call check(nf90_inq_varid(nc_id, 'refindex_im_aer_lw', vid_im))

         allocate(ref_real(nlwbands), ref_im(nlwbands))

         call check(nf90_get_var(nc_id, vid_real, ref_real),': ERROR reading refindex_real_aer_lw')
         call check(nf90_get_var(nc_id, vid_im, ref_im),': ERROR reading refindex_im_aer_lw')

         ! successfully read refindex data -- set complex value in physprop object
         do iw = 1, nlwbands
            refidx_aer_lw(imode,ispec(imode),iw) = cmplx(ref_real(iw), abs(ref_im(iw)))
         end do

         deallocate(ref_real, ref_im)

         if(ifp >= nfiles) exit
      enddo

    end subroutine read_modal_optics


    subroutine read_water_refindex(infilename, crefwsw, crefwlw)

      ! read water refractive index file and set module data

      character*(*), intent(in) :: infilename ! modal optics filename
      complex(r8), intent(out) :: crefwsw(nswb) ! complex refractive index for water visible
      complex(r8), intent(out) :: crefwlw(nlwb) ! complex refractive index for water infrared

      ! Local variables

      integer :: i, ierr
      integer :: ncid ! pio file handle
      integer :: did ! dimension ids
      integer :: dimlen ! dimension lengths
      integer :: vid ! variable ids
      real(r8) :: refrwsw(nswb), refiwsw(nswb) ! real, imaginary ref index for water visible
      real(r8) :: refrwlw(nlwb), refiwlw(nlwb) ! real, imaginary ref index for water infrared
      !----------------------------------------------------------------------------

      ! open file
      call check( nf90_open(trim(adjustl(infilename)), NF90_NOWRITE, ncid) )

      ! inquire dimensions. Check that file values match parameter values.

      call check(nf90_inq_dimid(ncid, 'lw_band', did))
      call check(nf90_inquire_dimension ( ncid, did, len = dimlen ))

      if (dimlen .ne. nlwb) then
         write(iulog,*) 'lw_band len=', dimlen, ' from ', infilename, ' ne nlwbands=', nlwb
         call endrun('read_modal_optics: bad lw_band value')
      endif

      call check(nf90_inq_dimid(ncid, 'sw_band', did))
      call check(nf90_inquire_dimension ( ncid, did, len = dimlen ))

      if (dimlen .ne. nswb) then
         write(iulog,*) 'sw_band len=', dimlen, ' from ', infilename, ' ne nswbands=', nswb
         call endrun('read_modal_optics: bad sw_band value')
      endif

      ! read variables
      call check(nf90_inq_varid(ncid, 'refindex_real_water_sw', vid))
      call check(nf90_get_var(ncid, vid, refrwsw),': ERROR reading refindex_real_water_sw')

      call check(nf90_inq_varid(ncid, 'refindex_im_water_sw', vid))
      call check(nf90_get_var(ncid, vid, refiwsw),': ERROR reading refindex_im_water_sw')

      call check(nf90_inq_varid(ncid, 'refindex_real_water_lw', vid))
      call check(nf90_get_var(ncid, vid, refrwlw),': ERROR reading refindex_real_water_lw')

      call check(nf90_inq_varid(ncid, 'refindex_im_water_lw', vid))
      call check(nf90_get_var(ncid, vid, refiwlw),': ERROR reading refindex_im_water_lw')

      ! set complex representation of refractive indices as module data
      do i = 1, nswb
         crefwsw(i) = cmplx(refrwsw(i), abs(refiwsw(i)))
      end do
      do i = 1, nlwb
         crefwlw(i) = cmplx(refrwlw(i), abs(refiwlw(i)))
      end do

      call check(nf90_close(ncid))

    end subroutine read_water_refindex



    subroutine rad_cnst_get_aer_props(list_idx, imode, ism, refindex_aer_sw, refindex_aer_lw)

      integer, intent(in) :: list_idx, imode, ism
      complex(r8), intent(out), optional :: refindex_aer_sw(:) ! species refractive index
      complex(r8), intent(out), optional :: refindex_aer_lw(:) ! species refractive index

      integer :: iw

      !Sanity chk for list_idx
      if(list_idx .ne. 0) then
         print*,'This code only support list_dix=0; list_idx is', list_idx
         stop "Stopped"
      endif

      if(present(refindex_aer_sw)) then
         do iw = 1, nswb
            refindex_aer_sw(iw) = refidx_aer_sw(imode,ism,iw)
         end do
      endif

      if(present(refindex_aer_lw)) then
         do iw = 1, nlwb
            refindex_aer_lw(iw) = refidx_aer_lw(imode,ism,iw)
         end do
      endif


    end subroutine rad_cnst_get_aer_props


    subroutine check(status, msg)
      character(len=*), intent(in), optional:: msg

      integer, intent (in) :: status

      if(status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         if(present(msg))print*,msg
         stop "Stopped"
      end if
    end subroutine check

end module module_optics_support
