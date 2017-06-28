!> \file modcape.f90
!!   Dumps cross sections of dcape, CAPEmax, CINlower,CINupper, CINmax, ALE, W2max,
!!   QTcloudbase, THLcloudbase, Wcloudbase, THVcloudbase, QLcloudbase, LWP, RWP, cloud top
!
!>
!! Crosssections in the xy-plane
!! If netcdf is true, this module leads the dcape.myid.expnr.nc output

!!  \par Revision list
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!
module modcape

  use modglobal, only : longint,kmax

implicit none
private
PUBLIC :: initcape,docape,exitcape
save
!NetCDF variables
  integer,parameter :: nvar = 6
  integer :: ncid4 = 0
  integer :: nrec = 0
  character(80) :: fname = 'cape.xxxxyxxx.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname
  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  logical :: lcape = .false. !< switch for doing the crosssection (on/off)

contains

!> Initializing cape crossections. Read out the namelist, initializing the variables
  subroutine initcape
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical,cmyid
    use modglobal,only :imax,jmax,ifnamopt,fname_options,dtmax,dtav_glob,ladaptive,dt_lim,cexpnr,tres,btime
    use modstat_nc,only : lnetcdf,open_nc, define_nc, redefine_nc,ncinfo,writestat_dims_nc
   implicit none

    integer :: ierr

    namelist/NAMCAPE/ &
    lcape, dtav

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCAPE,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMCAPE'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMCAPE'
      endif
      write(6 ,NAMCAPE)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtav      ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lcape     ,1,MPI_LOGICAL,0,comm3d,mpierr)

    idtav = dtav/tres
    tnext   = idtav+btime
    if(.not.(lcape)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'cape: dtav should be a integer multiple of dtmax'
    end if
    if (lnetcdf) then
    fname(6:13) = cmyid
    fname(15:17) = cexpnr
    call ncinfo(tncname(1,:),'time','Time','s','time')
    call ncinfo(ncname( 1,:),'lwp','xy crosssections liquid water path','kg/m^2','tt0t')
    call ncinfo(ncname( 2,:),'rwp','xy crosssections rain water path','kg/m^2','tt0t')
    call ncinfo(ncname( 3,:),'twp','total water path','kg/m^2','tt0t')
    call ncinfo(ncname( 4,:),'surfprec','surface precipitation','-','tt0t')
    call ncinfo(ncname( 5,:),'mfcadv','vertical integral of horizontal qv advection','kg/m^2','tt0t')
    call ncinfo(ncname( 6,:),'mfccon','vertical integral of horizontal qv convergence','kg/m^2','tt0t')
    call open_nc(fname,  ncid4,nrec,n1=imax,n2=jmax)
    if (nrec==0) then
      call define_nc( ncid4, 1, tncname)
      call writestat_dims_nc(ncid4)
    end if
    call define_nc( ncid4, NVar, ncname)
    end if

  end subroutine initcape

!>Run crosssection.
  subroutine docape
    use modglobal, only : imax,jmax,i1,j1,k1,nsv,rk3step,timee,rtimee,dt_lim,&
    nsv,dzf,dxi5,dyi5
    use modfields, only : qt0,ql0,sv0,rhobf,u0,v0
    use modstat_nc, only : lnetcdf, writestat_nc
    use modmicrodata, only : iqr, precep, imicro
    use modmpi
    implicit none

    real, allocatable :: lwp(:,:),twp(:,:),rwp(:,:),sprec(:,:)
    ! for how to calc the following two terms, see Banacos, Peter: moisture flux convergence: its history...
    real, allocatable :: mfcadv(:,:) ! moisture flux convergence; advection term
    real, allocatable :: mfccon(:,:) ! moisture flux convergence; convergence term
    real, allocatable :: vars(:,:,:)

    ! LOCAL VARIABLES
    integer :: i,j,k
    real :: du0dx,dv0dy

    if (.not. lcape) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    allocate(lwp(2:i1,2:j1),rwp(2:i1,2:j1),twp(2:i1,2:j1))
    allocate(sprec(2:i1,2:j1))
    allocate(mfcadv(2:i1,2:j1),mfccon(2:i1,2:j1))

    ! loops over i,j,k
    lwp=0.
    twp=0.
    rwp=0.
    sprec=0.
    mfcadv=0.
    mfccon=0.
    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! water paths and surface precip
      lwp(i,j)=lwp(i,j)+rhobf(k)*ql0(i,j,k)*dzf(k)
      twp(i,j)=twp(i,j)+rhobf(k)*qt0(i,j,k)*dzf(k)
      if(nsv>1 .AND. imicro>0)then 
        rwp(i,j)=rwp(i,j)+rhobf(k)*sv0(i,j,k,iqr)*dzf(k)
        if(k==1)sprec(i,j)=precep(i,j,k)*rhobf(k)! correct for density to find total rain mass-flux
      end if
      ! vertical integral of moisture flux convergence
      ! finite central difference for "gradients"
      du0dx  = ( u0(i+1,j,k) - u0(i-1,j,k) ) * dxi5
      dv0dy  = ( v0(i,j+1,k) - v0(i,j-1,k) ) * dyi5
      ! calculate the terms
      mfcadv(i,j) = mfcadv(i,j) + rhobf(k)*( (-1)*u0(i,j,k)* &
      & (( (qt0(i+1,j,k) - ql0(i+1,j,k)) - (qt0(i-1,j,k) - ql0(i-1,j,k)) ) * dxi5) + &
      & (-1)*v0(i,j,k)*(( (qt0(i,j+1,k) - ql0(i,j+1,k)) - (qt0(i,j-1,k) - ql0(i,j-1,k)) ) * dyi5) ) * dzf(k)
      mfccon(i,j) = mfccon(i,j) + rhobf(k)*( (-1)*(qt0(i,j,k) - ql0(i,j,k))*(du0dx+dv0dy) ) * dzf(k)
    enddo
    enddo
    enddo

    if (lnetcdf) then
      allocate(vars(1:imax,1:jmax,nvar))
      vars(:,:,1) = lwp(2:i1,2:j1)
      vars(:,:,2) = rwp(2:i1,2:j1)
      vars(:,:,3) = twp(2:i1,2:j1)
      vars(:,:,4)= sprec(2:i1,2:j1)
      vars(:,:,5)= mfcadv(2:i1,2:j1)
      vars(:,:,6)= mfccon(2:i1,2:j1)
      call writestat_nc(ncid4,1,tncname,(/rtimee/),nrec,.true.)
      call writestat_nc(ncid4,nvar,ncname(1:nvar,:),vars,nrec,imax,jmax)
      deallocate(vars)
    end if

    deallocate(lwp,twp,rwp,sprec,mfcadv,mfccon)

  end subroutine docape

!> Clean up when leaving the run
  subroutine exitcape
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(lcape .and. lnetcdf) then
    call exitstat_nc(ncid4)
    end if

  end subroutine exitcape

end module modcape
