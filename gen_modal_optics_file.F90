
program modal_optics

  !=======================================================================================
  !This code generates modal optics files for each mode which can then be used as an
  !input to the CESM (CAM) model. 
  !
  !Inputs: 
  !=======
  ! Files: Physprop files for each aerosol and water refindex file. User can copy-paste the
  !"mode_defs" section of atm_in (CESM1.2 or later) directly into the namelist 
  ! (namelist_optics.nml) for this code.
  !
  !Other inputs: 
  !
  !
  !Ouputs:Modal optics files with name: modal_optics_mode<mode #>.nc. For example, 
  !       modal_optics_mode1.nc, modal_optics_mode2.nc etc.
  !======
  !
  !
  !=======================================================================================

  !=== USE STATEMENTS === 
  use netcdf,                only: nf90_create, NF90_CLOBBER, nf90_def_dim, nf90_DOUBLE,          &
       nf90_def_var, nf90_char, nf90_put_att, nf90_enddef, nf90_put_var, nf90_GLOBAL, nf90_close

  use radconstants,          only: get_sw_spectral_boundaries, nswbands, nlwbands,                &
       wavenumber1_longwave
  use shr_kind_mod,          only: r8=>shr_kind_r8
  use module_optics_support, only: iulog, cs1, endrun, read_modal_optics, rad_cnst_get_aer_props, &
       ncoef, prefr, prefi, masterproc, read_water_refindex, check
  use modal_aero_data,       only: ntot_amode, nspec_amode, maxd_aspectype, alnsg_amode,          &
       sigmag_amode

  implicit none

  !=== LOCAL VARS ===
  integer, parameter  :: n_mode_str = 126!60       !Hardwired(It is hardwired in CESM aswell)
  character(len=cs1)  :: water_refindex_file   !File name for water ref indx
  character(len=cs1)  :: mode_defs(n_mode_str) !mode_defs must have a constant dimension as it is read from namelist
  real(r8), parameter :: rmmin = 0.01e-6_r8
  real(r8), parameter :: rmmax = 25.e-6_r8

  !Counters
  integer :: ns, nr, ni, ierr, nc, im
  integer :: m, l, ERR, nphysprop 

  real(r8) :: refrmin, refrmax, refimin, refimax, refr, refi
  
  real(r8) :: wavmaxsw(nswbands), wavminsw(nswbands), wavmidsw(nswbands)
  real(r8) :: wavmaxlw(nlwbands), wavminlw(nlwbands), wavmidlw(nlwbands)

  real(r8) :: refrtabsw_mo(prefr,nswbands) ! table of real refractive indices for aerosols
  real(r8) :: refitabsw_mo(prefi,nswbands) ! table of imag refractive indices for aerosols
  real(r8) :: refrtablw_mo(prefr,nlwbands) ! table of real refractive indices for aerosols
  real(r8) :: refitablw_mo(prefi,nlwbands) ! table of imag refractive indices for aerosols
  
  real(r8) :: extpsw_mo(ncoef,prefr,prefi,ntot_amode,nswbands) ! specific extinction
  real(r8) :: abspsw_mo(ncoef,prefr,prefi,ntot_amode,nswbands) ! specific absorption
  real(r8) :: asmpsw_mo(ncoef,prefr,prefi,ntot_amode,nswbands) ! asymmetry factor

  real(r8) :: extparam(ncoef,prefr,prefi,ntot_amode)
  real(r8) :: absparam(ncoef,prefr,prefi,ntot_amode)
  real(r8) :: asmparam(ncoef,prefr,prefi,ntot_amode)

  real(r8) :: extplw(ncoef,prefr,prefi,ntot_amode,nlwbands) ! specific extinction
  real(r8) :: absplw(ncoef,prefr,prefi,ntot_amode,nlwbands) ! specific absorption
  real(r8) :: asmplw(ncoef,prefr,prefi,ntot_amode,nlwbands) ! asymmetry factor
  
  real(r8) :: refrtab(prefr) ! table of real refractive indices for aerosols
  real(r8) :: refitab(prefi) ! table of imag refractive indices for aerosols

  complex(r8) :: crefwsw(nswbands) ! complex refractive index for water visible
  complex(r8) :: crefwlw(nlwbands) ! complex refractive index for water infrared

  complex(r8), allocatable :: specrefindex(:)     ! species refractive index

  !Read mode_defs and water ref indx from namelist
  namelist /optics_nl/ mode_defs, water_refindex_file
  open(101,file = 'namelist_optics.nml',status = 'old')
  read(101, nml=optics_nl, iostat=ierr)

  if (ierr /= 0) then
     call endrun('ERROR reading namelist')
  end if
  close(101)

  !Compute total # of physprop files(equals # of species) to read
  nphysprop = 0
  do im = 1, ntot_amode
     nphysprop = nphysprop + nspec_amode(im)
  end do
  
 print*, 'test01!!!' 
  !read refractive indicies from physprop files
  call read_modal_optics(ntot_amode,maxd_aspectype,nphysprop, mode_defs, nswbands, nlwbands)
 print*, 'test02!!!'

  !Read water refindex file to populate crefwsw and  crefwlw
  call read_water_refindex(water_refindex_file, crefwsw, crefwlw)

  !Compute alnsg_amode
  do im = 1, ntot_amode
     alnsg_amode(im) = log( sigmag_amode(im) )
  enddo
  

  write(iulog,*)'calculating aerosol optical coefficients'
  
  call get_sw_spectral_boundaries(wavminsw, wavmaxsw, 'm')
  
  allocate(specrefindex(nswbands))
  do  ns = 1, nswbands !short wave wavelength loop                    
     wavmidsw(ns) = 0.5*(wavminsw(ns) + wavmaxsw(ns))!wavelength (m)
     if (masterproc) write(iulog,*)'wavelength=',wavmidsw(ns)
     
     !first find min,max of real and imaginary parts of refractive index
     
     refrmin = real(crefwsw(ns))
     refrmax = real(crefwsw(ns)) 
     refimin = aimag(crefwsw(ns))
     refimax = aimag(crefwsw(ns))                  
     
     do m=1,ntot_amode    !aerosol mode loop                               
        do l=1,nspec_amode(m)! aerosol species loop 
           
           !real and imaginary parts of aerosol refractive index
           !BSINGH - Future work: Add 'ns' as an arg and avoid going through whole loop !BALLI
           call rad_cnst_get_aer_props(0, m, l, refindex_aer_sw=specrefindex)

           
           refr=real(specrefindex(ns))
           refi=aimag(specrefindex(ns))
           
           refrmin=min(refrmin,refr)
           refrmax=max(refrmax,refr)
           refimin=min(refimin,refi)
           refimax=max(refimax,refi)
           
        enddo
     enddo
     if (masterproc) then
        write(iulog,*)'sw band ',ns
        write(iulog,*)'refrmin=',refrmin
        write(iulog,*)'refrmax=',refrmax
        write(iulog,*)'refimin=',refimin
        write(iulog,*)'refimax=',refimax
     endif
     if(refimin<-1.e-3)then
        write(iulog,*)'negative refractive indices not allowed'
        call endrun('error in modal_aer_opt_init')
     endif
     
     call miefit(wavmidsw(ns),refrtab,refitab, extparam, absparam, asmparam )
     
     do nr=1,prefr
        refrtabsw_mo(nr,ns)=refrtab(nr)
     end do
     do ni=1,prefi
        refitabsw_mo(ni,ns)=refitab(ni)
     end do
     do  m=1,ntot_amode
        do  nr=1,prefr
           do  ni=1,prefi
              do  nc=1,ncoef
                 extpsw_mo(nc,nr,ni,m,ns) = extparam(nc,nr,ni,m)
                 abspsw_mo(nc,nr,ni,m,ns) = absparam(nc,nr,ni,m)
                 asmpsw_mo(nc,nr,ni,m,ns) = asmparam(nc,nr,ni,m)
              enddo
           enddo
        enddo
     enddo

  enddo !short wave wavelength loop ends
  
  !Deallocate specrefindex and again reallocate for the log wave
  deallocate(specrefindex)
  allocate(specrefindex(nlwbands))

  do ns = 1, nlwbands         
     wavmidlw(ns) = 0.5e-2*(1._r8/wavenumber1_longwave(ns) + 1._r8/wavenumber1_longwave(ns))
     
     !first find min,max of real and imaginary parts of infrared (10 micron) refractive index            
     !real and imaginary parts of water refractive index            
     refrmin = real(crefwlw(ns)) 
     refrmax = real(crefwlw(ns)) 
     refimin = aimag(crefwlw(ns))
     refimax = aimag(crefwlw(ns))

     do m = 1, ntot_amode  !aerosol mode loop                             
        do l = 1, nspec_amode(m)!aerosol species loop
           
           !real and imaginary parts of aerosol refractive index                  
           call rad_cnst_get_aer_props(0, m, l, refindex_aer_lw=specrefindex)
           
           refr    = real(specrefindex(ns))
           refi    = aimag(specrefindex(ns))
           refrmin = min(refrmin,refr)
           refrmax = max(refrmax,refr)
           refimin = min(refimin,refi)
           refimax = max(refimax,refi)
           
        enddo
     enddo
     
     if (masterproc) then
        write(iulog,*)'longwave band ',ns
        write(iulog,*)'refrmin=',refrmin
        write(iulog,*)'refrmax=',refrmax
        write(iulog,*)'refimin=',refimin
        write(iulog,*)'refimax=',refimax
     endif
     
     call miefit(wavmidlw(ns),refrtab,refitab, extparam, absparam, asmparam )
     
     do nr = 1, prefr
        refrtablw_mo(nr,ns) = refrtab(nr)
     end do
     
     do ni = 1, prefi
        refitablw_mo(ni,ns) = refitab(ni)
     end do

     do  m = 1, ntot_amode
        do  nr = 1, prefr
           do  ni = 1, prefi
              do  nc = 1, ncoef
                 extplw(nc,nr,ni,m,ns) = extparam(nc,nr,ni,m)
                 absplw(nc,nr,ni,m,ns) = absparam(nc,nr,ni,m)
                 asmplw(nc,nr,ni,m,ns) = asmparam(nc,nr,ni,m)
              enddo
           enddo
        enddo
     enddo

  enddo !long wave wavelength loop ends
  
  if (masterproc) then
     write(iulog,*) "modal_aer_opt_init: write modal optics files:"
     call write_modal_optics(refrtabsw_mo, refitabsw_mo, refrtablw_mo, refitablw_mo, sigmag_amode )
  endif

  deallocate(specrefindex)

contains

subroutine miefit(wavmid,refrtab,refitab, extparam,absparam,asmparam)
        
        use physconst,      only: rhoh2o
        
        implicit none
        integer nsiz ! number of wet particle sizes
        integer nlog ! number of log-normals modes to product mie fit
        parameter (nsiz=200,nlog=30)
        
        real(r8) rad(nsiz)
        real(r8) wavmid
        real(r8) extparam(ncoef,prefr,prefi,ntot_amode)
        real(r8) absparam(ncoef,prefr,prefi,ntot_amode)
        real(r8) asmparam(ncoef,prefr,prefi,ntot_amode)
        real(r8) refrtab(prefr),refitab(prefi)
        real(r8) size ! 2 pi radpart / waveleng = size parameter
        real(r8) qext(nsiz) ! array of extinction efficiencies
        real(r8) qsca(nsiz) ! array of scattering efficiencies
        real(r8) qabs(nsiz) ! array of absorption efficiencies
        real(r8) gqsc(nsiz) ! array of asymmetry factor * scattering efficiency
        real(r8) asymm(nsiz) ! array of asymmetry factor
        complex*16 sforw,sback,tforw(2),tback(2),refindx
        integer :: nmom,ipolzn,momdim,numang
        real(r8) pmom(0:1,1)
        logical :: perfct,prnt(2)
        complex*16 crefd
        real(r8) :: mimcut
        logical :: anyang
        
        real(r8) xmu(1)!BALLI
        complex*16 s1(1),s2(1)
        
        real(r8) dsdlogr(nsiz),dlogr
        real(r8) rmin,rmax ! min, max aerosol size bin
        real(r8) xr
        real(r8) drefr ! increment in real part of refractive index
        real(r8) drefi ! increment in imag part of refractive index
        integer m,n,nr,ni,nl
        
        real(r8) exparg
        real(r8) volwet ! sum of volume of wet aerosols
        real(r8) sumabs ! sum of specific absorption coefficients
        real(r8) sumsca ! sum of specific scattering coefficients
        real(r8) sumg   ! sum of asymmetry factors
        real(r8) rs(nlog) ! surface mode radius (cm)
        real(r8) specscat(nlog) ! specific scattering (m2/g)
        real(r8) specabs(nlog)  ! specific absorption (m2/g)
        real(r8) specext(nlog)  ! specfiic extinction (m2/g)
        real(r8) abs(nlog)      ! single scattering albedo
        real(r8) asym(nlog)     ! asymmetry factor
        real(r8) logr(nsiz)
        real(r8) bma,bpa
        real(r8) sq2pi
        real(r8) pie
        logical diag
        
        xmu(1) = 0.0_r8 !BALLI

        pie=4._r8*atan(1._r8)
        sq2pi=sqrt(2.*pie)
        
        nmom=0
        ipolzn=0
        momdim=1
        perfct=.false.
        prnt(1)=.false.
        prnt(2)=.false.
        mimcut=0.0_r8
        
        anyang=.false.
        numang=0
        if(masterproc.and.wavmid.gt.0.45e-6.and.wavmid.lt.0.6e-6)write(98,*)'wavmid=',wavmid
        
        drefr=(refrmax-refrmin)/(prefr-1)
        
        !size bins (m)
        
        rmin=0.001e-6
        rmax=100.e-6
        write(98,*)'wavmid=',wavmid
        dlogr=log(rmax/rmin)/(nsiz-1)
        logr(1)=log(rmin)
        do n=2,nsiz
           logr(n)=logr(n-1)+dlogr
        enddo
        
        do n=1,nsiz
           rad(n)=exp(logr(n))
        enddo
        
        !calibrate parameterization with range of refractive indices
        
        do  nr=1,prefr
           do  ni=1,prefi
              
              diag=masterproc.and.wavmid.gt.0.45e-6.and.wavmid.lt.0.6e-6.and.nr.eq.1.and.ni.eq.1
              
              refrtab(nr)=refrmin+(nr-1)*drefr
              refitab(ni)=refimax*0.3**(prefi-ni)
              if(ni.eq.1)refitab(ni)=0.
              crefd=dcmplx(refrtab(nr),refitab(ni))
              if(diag)write(98,*)'nr,ni,refr,refi=',nr,ni,refrtab(nr),refitab(ni)
              
              !mie calculations of optical efficiencies
              
              do n=1,nsiz
                 
                 !size parameter and weighted refractive index
                 
                 size=2.*pie*rad(n)/wavmid
                 size=min(size,400.d0)
                 
                 call miev0(size,crefd,perfct,mimcut,anyang,          &
                      numang,xmu,nmom,ipolzn,momdim,prnt,               &
                      qext(n),qsca(n),gqsc(n),pmom,sforw,sback,s1,      &
                      s2,tforw,tback )
                 
                 qsca(n)=min(qsca(n),qext(n))
                 qabs(n)=qext(n)-qsca(n)
                 asymm(n)=gqsc(n)/qsca(n)
                 if(diag)write(98,*)'rad=',rad(n),' qsca,qabs=',qsca(n),qabs(n)
              enddo
              
              !now consider a variety of lognormal size distributions
              !with surface mode radius rs (m) and log standard deviation alnsg_amode
              
              !aerosol mode loop
              
              do m=1,ntot_amode
                 
                 !size units are m
                 !bounds of surface mode radius
                 
                 bma=0.5*log(rmmax/rmmin)
                 bpa=0.5*log(rmmax*rmmin)
                 
                 do nl=1,nlog
                    
                    !pick rs to match zeroes of chebychev
                    
                    xr=cos(pie*(nl-0.5)/nlog)
                    rs(nl)=exp(xr*bma+bpa)
                    
                    !define wet size distribution
                    
                    !integrate over all size bins
                    
                    volwet=0
                    sumsca=0
                    sumabs=0
                    sumg=0
                    
                    do n=1,nsiz
                       exparg=log(rad(n)/rs(nl))/alnsg_amode(m)
                       dsdlogr(n)=exp(-0.5*exparg*exparg)
                       volwet=volwet+4./3.*rad(n)*dsdlogr(n)*dlogr
                       sumabs=sumabs+qabs(n)*dsdlogr(n)*dlogr
                       sumsca=sumsca+qsca(n)*dsdlogr(n)*dlogr
                       sumg=sumg+asymm(n)*qsca(n)*dsdlogr(n)*dlogr
                    enddo
                    
                    
                    !coefficients wrt wet mass, using water density
                    specscat(nl)=sumsca/(volwet*rhoh2o)
                    specabs(nl)=sumabs/(volwet*rhoh2o)
                    !coefficients wrt wet mass
                    specext(nl)=specscat(nl)+specabs(nl)
                    asym(nl)=sumg/sumsca
                    
                    if(masterproc)then
                       if(m.eq.2)then
                          write(98,'(a,e12.4,a,2f10.4,a,f10.4)')'rs=',rs(nl),' cref=',refrtab(nr),refitab(ni),' asym=',asym(nl)
                       endif
                    endif
                    
                 enddo
                 
                 call fitcurv(rs,specext,extparam(1,nr,ni,m),ncoef,nlog)
                 call fitcurvlin(rs,specabs,absparam(1,nr,ni,m),ncoef,nlog)
                 call fitcurvlin(rs,asym,asmparam(1,nr,ni,m),ncoef,nlog)
                 
              enddo !ntot_amode loop
           enddo !ni
        enddo !nr
        return
      end subroutine miefit
      !===============================================================================
      
      subroutine fitcurv(rs,yin,coef,ncoef,maxm)
        
        !fit y(x) using Chebyshev polynomials
        
        implicit none
        integer nmodes,maxm,m,ncoef
        parameter (nmodes=300)
        
        real(r8) rs(maxm) ! surface mode radius (cm)
        real(r8) yin(maxm),coef(ncoef)
        real(r8) x(nmodes),y(nmodes)
        real(r8) xmin,xmax
        
        if(maxm.gt.nmodes)then
           write(iulog,*)'nmodes too small in fitcurv',maxm
           call endrun()
        endif
        
        xmin=1.e20
        xmax=-1.e20
        do m=1,maxm
           x(m)=log(rs(m))
           xmin=min(xmin,x(m))
           xmax=max(xmax,x(m))
           y(m)=log(yin(m))
        enddo
        
        do  m=1,maxm
           x(m)=(2*x(m)-xmax-xmin)/(xmax-xmin)
        enddo
        
        call chebft(coef,ncoef,y,maxm)
        
        return
      end subroutine fitcurv

!===============================================================================

      subroutine fitcurvlin(rs,yin,coef,ncoef,maxm)
        
        !fit y(x) using Chebychev polynomials
        !calculates ncoef coefficients coef
        
        implicit none
        integer nmodes,m,ncoef
        parameter (nmodes=300)
        
        integer,intent(in) :: maxm
        real(r8),intent(in) :: rs(maxm),yin(maxm)
        real(r8) x(nmodes),y(nmodes)
        real(r8),intent(out) :: coef(ncoef)
        real(r8) xmin,xmax
        
        if(maxm.gt.nmodes)then
           write(iulog,*)'nmodes too small in fitcurv',maxm
           call endrun
        endif
        
        do  m=1,maxm
           x(m)=log(rs(m))
           y(m)=(yin(m))
        enddo
        
        call chebft(coef,ncoef,y,maxm)
        
        return
      end subroutine fitcurvlin

!===============================================================================

      subroutine chebft(c,m,f,n)
        !given a function f with values at zeroes x_k of Chebychef polynomial
        !T_n(x), calculate coefficients c_j such that
        !f(x)=sum(k=1,n) c_k t_(k-1)(y) - 0.5*c_1
        !where y=(x-0.5*(xmax+xmin))/(0.5*(xmax-xmin))
        !See Numerical Recipes, pp. 148-150.

        implicit none
        real(r8) pi
        parameter (pi=3.14159265)
        real(r8) c(m),f(n),fac,sum
        integer n,j,k,m
        
        fac=2./n
        do j=1,m
           sum=0
           do k=1,n
              sum=sum+f(k)*cos((pi*(j-1))*((k-0.5)/n))
           enddo
           c(j)=fac*sum
        enddo
        return
      end subroutine chebft

      !===============================================================================
      !BSINGH - Added following subroutines for writing out modal optics file
      subroutine write_modal_optics(refindex_real_sw,refindex_im_sw, refindex_real_lw, refindex_im_lw, sigma_logr_aer )
        
        use modal_aero_data, only: rhcrystal_amode, rhdeliques_amode, dgnum_amode, dgnumlo_amode, dgnumhi_amode
        implicit none
        
        character(len=256) outfilename
        ! error status return
        integer  ierr
        ! netCDF id
        integer ::  ncid
        ! dimension ids
        integer  lw_band_dim
        integer  sw_band_dim
        integer  refindex_real_dim
        integer  refindex_im_dim
        integer  mode_dim
        integer  spec_dim
        
        ! variable ids
        
        integer ::   lw_abs_id
        integer ::   sw_asm_id
        integer ::   sw_ext_id
        integer ::   sw_abs_id
        integer ::   sigma_logr_aer_id
        integer ::   dgnum_id, dgnumlo_id, dgnumhi_id
        integer ::   rhcrystal_id, rhdeliques_id
        integer ::   refindex_real_lw_id
        integer ::   refindex_im_lw_id
        integer ::   refindex_real_sw_id
        integer ::   refindex_im_sw_id
        integer ::   sigmag_id
        integer ::   opticsmethod_id
        ! rank (number of dimensions) for each variable
        integer  lw_abs_rank
        integer  sw_asm_rank
        integer  sw_ext_rank
        integer  sw_abs_rank
        integer refindex_rank
        parameter (lw_abs_rank = 5)
        parameter (sw_asm_rank = 5)
        parameter (sw_ext_rank = 5)
        parameter (sw_abs_rank = 5)
        parameter (refindex_rank=2)
        ! variable shapes
        integer  lw_abs_dims(lw_abs_rank)
        integer  sw_asm_dims(sw_asm_rank)
        integer  sw_ext_dims(sw_ext_rank)
        integer  sw_abs_dims(sw_abs_rank)
        integer refindex_dims(refindex_rank)
        integer coef_dim
        
        ! data variables
        
        real(r8), intent(in) ::  sigma_logr_aer(ntot_amode) ! geometric standard deviation of size distribution
        
        !     tables of real refractive indices for aerosols
        real(r8), intent(in) :: refindex_real_sw(prefr,nswbands) !
        real(r8), intent(in) :: refindex_im_sw(prefi,nswbands) !
        real(r8), intent(in) :: refindex_real_lw(prefr,nlwbands) !
        real(r8), intent(in) :: refindex_im_lw(prefi,nlwbands) !
        
        real(r8) :: extpswt(ncoef,prefr,prefi,1,nswbands) ! specific extinction
        real(r8) :: abspswt(ncoef,prefr,prefi,1,nswbands) ! specific absorption
        real(r8) :: asmpswt(ncoef,prefr,prefi,1,nswbands) ! asymmetry factor
        real(r8) :: absplwt(ncoef,prefr,prefi,1,nlwbands) ! specific absorption lw
        
        integer lenchr
        integer l,m
        integer nr,ni,nc,ns
        integer opticsmethod_len_dim
        integer, parameter :: opticsmethod_len=32
        character(len=opticsmethod_len)opticsmethod        
        
        integer len
        
        ! create netcdf file
        
        opticsmethod='modal'
        
        
        do m=1, ntot_amode

           if(m.lt.10)then
              write(outfilename,'(a17,i1,a3)')'modal_optics_mode',m,'.nc'
           else
              write(outfilename,'(a17,i2,a3)')'modal_optics_mode',m,'.nc'
           endif

           
           call check( nf90_create(outfilename, NF90_CLOBBER, ncid) )
           
           ! define dimensions
           
           call check( nf90_def_dim(ncid, 'lw_band', nlwbands, lw_band_dim) )
           call check( nf90_def_dim(ncid, 'sw_band', nswbands, sw_band_dim))
           call check( nf90_def_dim(ncid, 'refindex_real', prefr, refindex_real_dim))
           call check( nf90_def_dim(ncid, 'refindex_im', prefi, refindex_im_dim))
           call check( nf90_def_dim(ncid, 'mode', 1, mode_dim))
           call check( nf90_def_dim(ncid, 'coef_number', ncoef, coef_dim))
           call check( nf90_def_dim(ncid, 'opticsmethod_len', opticsmethod_len, opticsmethod_len_dim))
           
           ! define variables

           lw_abs_dims(5) = lw_band_dim
           lw_abs_dims(4) = mode_dim
           lw_abs_dims(3) = refindex_im_dim
           lw_abs_dims(2) = refindex_real_dim
           lw_abs_dims(1) = coef_dim
           call check( nf90_def_var(ncid, 'absplw', nf90_DOUBLE, lw_abs_dims, lw_abs_id ))
           
           sw_asm_dims(5) = sw_band_dim
           sw_asm_dims(4) = mode_dim
           sw_asm_dims(3) = refindex_im_dim
           sw_asm_dims(2) = refindex_real_dim
           sw_asm_dims(1) = coef_dim
           call check( nf90_def_var(ncid, 'asmpsw', nf90_DOUBLE, sw_asm_dims, sw_asm_id ))
           
           sw_ext_dims(5) = sw_band_dim
           sw_ext_dims(4) = mode_dim
           sw_ext_dims(3) = refindex_im_dim
           sw_ext_dims(2) = refindex_real_dim
           sw_ext_dims(1) = coef_dim
           call check( nf90_def_var(ncid, 'extpsw', nf90_DOUBLE, sw_ext_dims, sw_ext_id ))
           
           sw_abs_dims(5) = sw_band_dim
           sw_abs_dims(4) = mode_dim
           sw_abs_dims(3) = refindex_im_dim
           sw_abs_dims(2) = refindex_real_dim
           sw_abs_dims(1) = coef_dim
           call check( nf90_def_var(ncid, 'abspsw', nf90_DOUBLE,  sw_abs_dims, sw_abs_id ))
           
           refindex_dims(2) = sw_band_dim
           refindex_dims(1) = refindex_real_dim
           call check( nf90_def_var(ncid, 'refindex_real_sw', nf90_DOUBLE,  refindex_dims, refindex_real_sw_id))
           
           refindex_dims(2) = sw_band_dim
           refindex_dims(1) = refindex_im_dim
           call check( nf90_def_var(ncid, 'refindex_im_sw', nf90_DOUBLE,  refindex_dims, refindex_im_sw_id))
           
           refindex_dims(2) = lw_band_dim
           refindex_dims(1) = refindex_real_dim
           call check( nf90_def_var(ncid, 'refindex_real_lw', nf90_DOUBLE,  refindex_dims,  refindex_real_lw_id))
           
           refindex_dims(2) = lw_band_dim
           refindex_dims(1) = refindex_im_dim
           call check( nf90_def_var(ncid, 'refindex_im_lw', nf90_DOUBLE,  refindex_dims, refindex_im_lw_id))

           call check( nf90_def_var(ncid, 'sigmag', nf90_DOUBLE,  sigmag_id))
           call check( nf90_def_var(ncid, 'opticsmethod', nf90_char, opticsmethod_len_dim, opticsmethod_id))
           
           call check( nf90_def_var(ncid, 'dgnum', nf90_DOUBLE,  dgnum_id))
           call check( nf90_def_var(ncid, 'dgnumlo', nf90_DOUBLE,  dgnumlo_id))
           call check( nf90_def_var(ncid, 'dgnumhi', nf90_DOUBLE,  dgnumhi_id))
           
           call check( nf90_def_var(ncid, 'rhcrystal', nf90_DOUBLE,  rhcrystal_id))
           call check( nf90_def_var(ncid, 'rhdeliques', nf90_DOUBLE, rhdeliques_id))
                      
           ! assign attributes                      

           call check( nf90_put_att(ncid, lw_abs_id, 'long_name', 'coefficients of polynomial expression for longwave absorption'))

           call check( nf90_put_att(ncid, lw_abs_id, 'units', 'meter^2 kilogram^-1'))
           call check( nf90_put_att(ncid, sw_ext_id, 'long_name', 'coefficients of polynomial expression for shortwave extinction')) 
           call check( nf90_put_att(ncid, sw_ext_id, 'units',  'meter^2 kilogram^-1'))
           call check( nf90_put_att(ncid, sw_asm_id, 'long_name', 'coefficients of polynomial expression for shortwave asymmetry parameter'))
           call check( nf90_put_att(ncid, sw_asm_id, 'units',  'fraction'))
           call check( nf90_put_att(ncid, sw_abs_id, 'long_name', 'coefficients of polynomial expression for shortwave absorption'))
           call check( nf90_put_att(ncid, sw_abs_id, 'units',  'meter^2 kilogram^-1'))
           call check( nf90_put_att(ncid, refindex_real_lw_id, 'long_name', 'real refractive index of aerosol - longwave'))
           call check( nf90_put_att(ncid, refindex_im_lw_id,   'long_name',   'imaginary refractive index of aerosol - longwave'))
           call check( nf90_put_att(ncid, refindex_real_sw_id, 'long_name', 'real refractive index of aerosol - shortwave'))
           call check( nf90_put_att(ncid, refindex_im_sw_id,   'long_name', 'imaginary refractive index of aerosol - shortwave'))
           call check( nf90_put_att(ncid, sigmag_id, 'long_name', 'geometric standard deviation of aerosol'))
           call check( nf90_put_att(ncid, sigmag_id, 'units',  '-'))
           call check( nf90_put_att(ncid, dgnum_id, 'long_name', 'nominal geometric mean diameter for number distribution'))
           call check( nf90_put_att(ncid, dgnum_id, 'units',  'm'))
           call check( nf90_put_att(ncid, dgnumlo_id, 'long_name', 'lower bound on geometric mean diameter for number distribution'))
           call check( nf90_put_att(ncid, dgnumlo_id, 'units',  'm'))
           call check( nf90_put_att(ncid, dgnumhi_id, 'long_name', 'upper bound on geometric mean diameter for number distribution'))
           call check( nf90_put_att(ncid, dgnumhi_id, 'units',  'm'))
           call check( nf90_put_att(ncid, rhcrystal_id, 'long_name', 'crystalization relative humidity'))
           call check( nf90_put_att(ncid,rhcrystal_id, 'units',  '-'))
           call check( nf90_put_att(ncid, rhdeliques_id, 'long_name', 'deliquesence relative humidity'))
           call check( nf90_put_att(ncid,rhdeliques_id, 'units',  '-'))
           
           call check( nf90_put_att(ncid, nf90_GLOBAL, 'source',  'OPAC + miev0'))
           call check( nf90_put_att(ncid, nf90_GLOBAL, 'history',  'Ghan and Zaveri, JGR 2007; Longlei Li, 2023'))
           
           call check( nf90_enddef(ncid))
           
           call check( nf90_put_var(ncid, refindex_real_lw_id, refindex_real_lw))

           call check( nf90_put_var(ncid, refindex_im_lw_id, refindex_im_lw))
           call check( nf90_put_var(ncid, refindex_real_sw_id, refindex_real_sw))
           call check( nf90_put_var(ncid, refindex_im_sw_id, refindex_im_sw))
           
           do  nr=1,prefr
              do  ni=1,prefi
                 do  nc=1,ncoef
                    do ns=1, nswbands
                       extpswt(nc,nr,ni,1,ns)=extpsw_mo(nc,nr,ni,m,ns)
                       abspswt(nc,nr,ni,1,ns)=abspsw_mo(nc,nr,ni,m,ns)
                       asmpswt(nc,nr,ni,1,ns)=asmpsw_mo(nc,nr,ni,m,ns)
                    end do
                    do ns=1,nlwbands
                       absplwt(nc,nr,ni,1,ns)=absplw(nc,nr,ni,m,ns)
                    end do
                 enddo
              enddo
           enddo           
           
           call check( nf90_put_var(ncid, lw_abs_id, absplwt))
           call check( nf90_put_var(ncid, sw_ext_id, extpswt))
           call check( nf90_put_var(ncid, sw_abs_id, abspswt))
           call check( nf90_put_var(ncid, sw_asm_id, asmpswt))
           
           call check( nf90_put_var(ncid,  sigmag_id,  sigma_logr_aer(m)))
           call check( nf90_put_var(ncid,  opticsmethod_id, opticsmethod) )
           call check( nf90_put_var(ncid,  dgnum_id,  dgnum_amode(m)))
           call check( nf90_put_var(ncid,  dgnumlo_id,  dgnumlo_amode(m)))
           call check( nf90_put_var(ncid,  dgnumhi_id,  dgnumhi_amode(m)))
           call check( nf90_put_var(ncid,  rhcrystal_id,  rhcrystal_amode(m)))
           call check( nf90_put_var(ncid,  rhdeliques_id,  rhdeliques_amode(m)))
                                 
           ierr= nf90_close(ncid)
           print *,'modal optics file closed for mode ',m


        end do ! mode loop

      end subroutine write_modal_optics




end program modal_optics


