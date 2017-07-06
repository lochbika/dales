!> \file advec_52.f90
!!  Does advection with a 5th order upwind scheme and 2nd order in the vertical
!! \par Revision list
!! \par Authors
!! \see Wicker and Scamarock 2002
!!
!! By adding a small dissipative term to the sixth order flux, a fifth order
!! scheme is created that is nearly monotone:
!! \latexonly
!! \begin{eqnarray}
!!  F_{i-\frac{1}{2}}^{5th} &=& F_{i-\frac{1}{2}}^{6th} -
!! \left|\frac{\fav{u}_{i-\frac{1}{2}}}{60}\right|\left[10(\phi_i-\phi_{i-1})\right
!! . \nonumber\\\\
!! &&-\left.5(\phi_{i+1}-\phi_{i-2})+(\phi_{i+2}-\phi_{i-3})\right].
!! \end{eqnarray}
!! \endlatexonly
!!
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

!> Advection at cell center
subroutine advecc_52(putin, putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi,dyi,dzf
  use modfields, only : u0, v0, w0,rhobf
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the cell centered field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency
!  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rhoputin
!  real                                      :: rhobfinvk                      
  real                                       :: inv2dzfk, rhobf_p, rhobf_m

  integer :: i,j,k

  !if (leq) then

!  do k=1,k1
!    do j=2-jh,j1+jh
!      do i=2-ih,i1+ih
!      rhoputin(i,j,k)=rhobf(k)*putin(i,j,k)
!      end do
!    end do
!  end do

  do k=1,kmax
     !rhobfinvk = 1./rhobf(k)
     inv2dzfk = 1./(2. * dzf(k))
     rhobf_p = rhobf(k+1)/rhobf(k)
     if (k > 1) rhobf_m = rhobf(k-1)/rhobf(k)
    do j=2,j1
      do i=2,i1

        if(k==1) then

          putout(i,j,k)  = putout(i,j,k)- ( &
                ( &
                  u0(i+1,j,k)/60.&
                  *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                  -abs(u0(i+1,j,k))/60.&
                  *(10.*(putin(i+1,j,k)-putin(i,j,k))-5.*(putin(i+2,j,k)-putin(i-1,j,k))+(putin(i+3,j,k)-putin(i-2,j,k)))&
                  -u0(i,j,k)/60.&
                  *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                  +abs(u0(i,j,k))/60.&
                  *(10.*(putin(i,j,k)-putin(i-1,j,k))-5.*(putin(i+1,j,k)-putin(i-2,j,k))+(putin(i+2,j,k)-putin(i-3,j,k)))&
                  )*dxi&
                +(&
                  v0(i,j+1,k)/60.&
                  *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                  -abs(v0(i,j+1,k))/60.&
                  *(10.*(putin(i,j+1,k)-putin(i,j,k))-5.*(putin(i,j+2,k)-putin(i,j-1,k))+(putin(i,j+3,k)-putin(i,j-2,k)))&
                  -v0(i,j,k)/60.&
                  *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  +abs(v0(i,j,k))/60.&
                  *(10.*(putin(i,j,k)-putin(i,j-1,k))-5.*(putin(i,j+1,k)-putin(i,j-2,k))+(putin(i,j+2,k)-putin(i,j-3,k))) &
                  )* dyi &
                + ( &
                w0(i,j,k+1) * (rhobf_p * putin(i,j,k+1) + putin(i,j,k)) &
                ) * inv2dzfk  &
                )

        else
            putout(i,j,k)  = putout(i,j,k)- (  &
                  ( &
                      u0(i+1,j,k)/60.&
                      *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                      -abs(u0(i+1,j,k))/60.&
                      *(10.*(putin(i+1,j,k)-putin(i,j,k))-5.*(putin(i+2,j,k)-putin(i-1,j,k))+(putin(i+3,j,k)-putin(i-2,j,k)))&
                      -u0(i,j,k)/60.&
                      *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                      +abs(u0(i,j,k))/60.&
                      *(10.*(putin(i,j,k)-putin(i-1,j,k))-5.*(putin(i+1,j,k)-putin(i-2,j,k))+(putin(i+2,j,k)-putin(i-3,j,k)))&
                  )*dxi&
                +(&
                      v0(i,j+1,k)/60.&
                      *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                      -abs(v0(i,j+1,k))/60.&
                      *(10.*(putin(i,j+1,k)-putin(i,j,k))-5.*(putin(i,j+2,k)-putin(i,j-1,k))+(putin(i,j+3,k)-putin(i,j-2,k)))&
                      -v0(i,j,k)/60.&
                      *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                      +abs(v0(i,j,k))/60.&
                      *(10.*(putin(i,j,k)-putin(i,j-1,k))-5.*(putin(i,j+1,k)-putin(i,j-2,k))+(putin(i,j+2,k)-putin(i,j-3,k)))&
                  )* dyi &
                + ( &
                  w0(i,j,k+1) * (rhobf_p * putin(i,j,k+1) +  putin(i,j,k)) &
                  -w0(i,j,k)  * (rhobf_m * putin(i,j,k-1) +  putin(i,j,k)) &
                  ) * inv2dzfk &
                  )
        end if

      end do
    end do
  end do

end subroutine advecc_52


!> Advection at the u point.
subroutine advecu_52(putin,putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzf
  use modfields, only : u0, v0, w0,rhobf
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the u field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency
  !real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rhoputin

  integer :: i,j,k

  !if (leq) then

!  do k=1,k1
!    do j=2-jh,j1+jh
!      do i=2-ih,i1+ih
!      rhoputin(i,j,k)=rhobf(k)*putin(i,j,k)
!      end do
!    end do
!  end do

  do k=1,kmax
    do j=2,j1
      do i=2,i1

      if(k==1) then

        putout(i,j,k)  = putout(i,j,k)- ( &
              (&
                  (u0(i+1,j,k)+u0(i,j,k))/60.&
                  *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                  -sign(1.,(u0(i+1,j,k)+u0(i,j,k)))*(u0(i+1,j,k)+u0(i,j,k))/60.&
                  *(10.*(putin(i+1,j,k)-putin(i,j,k))-5.*(putin(i+2,j,k)-putin(i-1,j,k))+(putin(i+3,j,k)-putin(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i-1,j,k))/60.&
                  *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                  +sign(1.,(u0(i,j,k)+u0(i-1,j,k)))*(u0(i,j,k)+u0(i-1,j,k))/60.&
                  *(10.*(putin(i,j,k)-putin(i-1,j,k))-5.*(putin(i+1,j,k)-putin(i-2,j,k))+(putin(i+2,j,k)-putin(i-3,j,k)))&
              )*dxi5 &
            +(&
                  (v0(i,j+1,k)+v0(i-1,j+1,k))/60.&
                  *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                  -sign(1.,(v0(i,j+1,k)+v0(i-1,j+1,k)))*(v0(i,j+1,k)+v0(i-1,j+1,k))/60.&
                  *(10.*(putin(i,j+1,k)-putin(i,j,k))-5.*(putin(i,j+2,k)-putin(i,j-1,k))+(putin(i,j+3,k)-putin(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i-1,j,k))/60.&
                  *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  +sign(1.,(v0(i,j,k)+v0(i-1,j,k)))*(v0(i,j,k)+v0(i-1,j,k))/60.&
                  *(10.*(putin(i,j,k)-putin(i,j-1,k))-5.*(putin(i,j+1,k)-putin(i,j-2,k))+(putin(i,j+2,k)-putin(i,j-3,k)))&
              )* dyi5 &
            +(1./rhobf(k))*( &
              ( rhobf(k+1)*putin(i,j,k+1) + rhobf(k) * putin(i,j,k)) *(w0(i,j,k+1)+ w0(i-1,j,k+1)) &
              ) / (4.*dzf(k)) &
              )

      else
        putout(i,j,k)  = putout(i,j,k)- ( &
              ( &
                  (u0(i+1,j,k)+u0(i,j,k))/60.&
                  *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                  -sign(1.,(u0(i+1,j,k)+u0(i,j,k)))*(u0(i+1,j,k)+u0(i,j,k))/60.&
                  *(10.*(putin(i+1,j,k)-putin(i,j,k))-5.*(putin(i+2,j,k)-putin(i-1,j,k))+(putin(i+3,j,k)-putin(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i-1,j,k))/60.&
                  *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                  +sign(1.,(u0(i,j,k)+u0(i-1,j,k)))*(u0(i,j,k)+u0(i-1,j,k))/60.&
                  *(10.*(putin(i,j,k)-putin(i-1,j,k))-5.*(putin(i+1,j,k)-putin(i-2,j,k))+(putin(i+2,j,k)-putin(i-3,j,k)))&
              )*dxi5&
            +(&
                  (v0(i,j+1,k)+v0(i-1,j+1,k))/60.&
                  *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                  -sign(1.,(v0(i,j+1,k)+v0(i-1,j+1,k)))*(v0(i,j+1,k)+v0(i-1,j+1,k))/60.&
                  *(10.*(putin(i,j+1,k)-putin(i,j,k))-5.*(putin(i,j+2,k)-putin(i,j-1,k))+(putin(i,j+3,k)-putin(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i-1,j,k))/60.&
                  *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  +sign(1.,(v0(i,j,k)+v0(i-1,j,k)))*(v0(i,j,k)+v0(i-1,j,k))/60.&
                  *(10.*(putin(i,j,k)-putin(i,j-1,k))-5.*(putin(i,j+1,k)-putin(i,j-2,k))+(putin(i,j+2,k)-putin(i,j-3,k)))&
              )* dyi5 &
            +(1./rhobf(k))*( &
              (rhobf(k) * putin(i,j,k) + rhobf(k+1) * putin(i,j,k+1) )*(w0(i,j,k+1)+w0(i-1,j,k+1)) &
            -(rhobf(k) * putin(i,j,k) + rhobf(k-1) * putin(i,j,k-1) )*(w0(i,j,k  )+w0(i-1,j,k  )) &
              ) / (4. * dzf(k)) &
              )
      end if

      end do
    end do
  end do

end subroutine advecu_52


!> Advection at the v point.
subroutine advecv_52(putin, putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzf
  use modfields, only : u0, v0, w0,rhobf
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the v field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency
  !real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rhoputin

  integer :: i,j,k

  !if (leq) then

  !do k=1,k1
  !  do j=2-jh,j1+jh
  !    do i=2-ih,i1+ih
  !    rhoputin(i,j,k)=rhobf(k)*putin(i,j,k)
  !    end do
  !  end do
  !end do

  do k=1,kmax
    do j=2,j1
      do i=2,i1

      if(k==1) then

        putout(i,j,k)  = putout(i,j,k)- ( &
              ( &
                  (u0(i+1,j,k)+u0(i+1,j-1,k))/60.&
                  *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                  -sign(1.,(u0(i+1,j,k)+u0(i+1,j-1,k)))*(u0(i+1,j,k)+u0(i+1,j-1,k))/60.&
                  *(10.*(putin(i+1,j,k)-putin(i,j,k))-5.*(putin(i+2,j,k)-putin(i-1,j,k))+(putin(i+3,j,k)-putin(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i,j-1,k))/60.&
                  *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                  +sign(1.,(u0(i,j,k)+u0(i,j-1,k)))*(u0(i,j,k)+u0(i,j-1,k))/60.&
                  *(10.*(putin(i,j,k)-putin(i-1,j,k))-5.*(putin(i+1,j,k)-putin(i-2,j,k))+(putin(i+2,j,k)-putin(i-3,j,k)))&
                )*dxi5&
              +(&
                  (v0(i,j+1,k)+v0(i,j,k))/60.&
                  *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                  -sign(1.,(v0(i,j+1,k)+v0(i,j,k)))*(v0(i,j+1,k)+v0(i,j,k))/60.&
                  *(10.*(putin(i,j+1,k)-putin(i,j,k))-5.*(putin(i,j+2,k)-putin(i,j-1,k))+(putin(i,j+3,k)-putin(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i,j-1,k))/60.&
                  *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  +sign(1.,(v0(i,j,k)+v0(i,j-1,k)))*(v0(i,j,k)+v0(i,j-1,k))/60.&
                  *(10.*(putin(i,j,k)-putin(i,j-1,k))-5.*(putin(i,j+1,k)-putin(i,j-2,k))+(putin(i,j+2,k)-putin(i,j-3,k)))&
                )* dyi5 &
              +(1./rhobf(k))*( &
                (w0(i,j,k+1)+w0(i,j-1,k+1)) *(rhobf(k+1) * putin(i,j,k+1) + rhobf(k) * putin(i,j,k)) &
                ) / (4. * dzf(k)) &
                )

      else
        putout(i,j,k)  = putout(i,j,k)- ( &
              ( &
                  (u0(i+1,j,k)+u0(i+1,j-1,k))/60.&
                  *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                  -sign(1.,(u0(i+1,j,k)+u0(i+1,j-1,k)))*(u0(i+1,j,k)+u0(i+1,j-1,k))/60.&
                  *(10.*(putin(i+1,j,k)-putin(i,j,k))-5.*(putin(i+2,j,k)-putin(i-1,j,k))+(putin(i+3,j,k)-putin(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i,j-1,k))/60.&
                  *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                  +sign(1.,(u0(i,j,k)+u0(i,j-1,k)))*(u0(i,j,k)+u0(i,j-1,k))/60.&
                  *(10.*(putin(i,j,k)-putin(i-1,j,k))-5.*(putin(i+1,j,k)-putin(i-2,j,k))+(putin(i+2,j,k)-putin(i-3,j,k)))&
                )*dxi5&
              +(&
                  (v0(i,j+1,k)+v0(i,j,k))/60.&
                  *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                  -sign(1.,(v0(i,j+1,k)+v0(i,j,k)))*(v0(i,j+1,k)+v0(i,j+1,k))/60.&
                  *(10.*(putin(i,j+1,k)-putin(i,j,k))-5.*(putin(i,j+2,k)-putin(i,j-1,k))+(putin(i,j+3,k)-putin(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i,j-1,k))/60.&
                  *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  +sign(1.,(v0(i,j,k)+v0(i,j-1,k)))*(v0(i,j,k)+v0(i,j-1,k))/60.&
                  *(10.*(putin(i,j,k)-putin(i,j-1,k))-5.*(putin(i,j+1,k)-putin(i,j-2,k))+(putin(i,j+2,k)-putin(i,j-3,k)))&
                )* dyi5 &
              +(1./rhobf(k))*( &
                (w0(i,j,k+1)+w0(i,j-1,k+1))*(rhobf(k+1) * putin(i,j,k+1) + rhobf(k) * putin(i,j,k)) &
                -(w0(i,j,k) +w0(i,j-1,k))  *(rhobf(k-1) * putin(i,j,k-1) + rhobf(k) * putin(i,j,k)) &
                ) / (4. * dzf(k)) &
                )
      end if

      end do
    end do
  end do

end subroutine advecv_52


!> Advection at the w point.
subroutine advecw_52(putin, putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzh
  use modfields, only : u0, v0, w0,rhobh
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the w field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency
  !real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rhoputin

  integer :: i,j,k

  !if (leq) then

  !do k=1,k1
  !  do j=2-jh,j1+jh
  !    do i=2-ih,i1+ih
  !    rhoputin(i,j,k)=rhobh(k)*putin(i,j,k)
  !    end do
  !  end do
  !end do

  do k=2,kmax
    do j=2,j1
      do i=2,i1
          putout(i,j,k)  = putout(i,j,k)- ( &
                (&
                    (u0(i+1,j,k)+u0(i+1,j,k-1))/60.&
                    *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                    -sign(1.,(u0(i+1,j,k)+u0(i+1,j,k-1)))*(u0(i+1,j,k)+u0(i+1,j,k-1))/60.&
                    *(10.*(putin(i+1,j,k)-putin(i,j,k))-5.*(putin(i+2,j,k)-putin(i-1,j,k))+(putin(i+3,j,k)-putin(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i,j,k-1))/60.&
                    *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                    +sign(1.,(u0(i,j,k)+u0(i,j,k-1)))*(u0(i,j,k)+u0(i,j,k-1))/60.&
                    *(10.*(putin(i,j,k)-putin(i-1,j,k))-5.*(putin(i+1,j,k)-putin(i-2,j,k))+(putin(i+2,j,k)-putin(i-3,j,k)))&
                )*dxi5&
              + (&
                    (v0(i,j+1,k)+v0(i,j+1,k-1))/60.&
                    *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                    -sign(1.,(v0(i,j+1,k)+v0(i,j+1,k-1)))*(v0(i,j+1,k)+v0(i,j+1,k-1))/60.&
                    *(10.*(putin(i,j+1,k)-putin(i,j,k))-5.*(putin(i,j+2,k)-putin(i,j-1,k))+(putin(i,j+3,k)-putin(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i,j,k-1))/60.&
                    *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                    +sign(1.,(v0(i,j,k)+v0(i,j,k-1)))*(v0(i,j,k)+v0(i,j,k-1))/60.&
                    *(10.*(putin(i,j,k)-putin(i,j-1,k))-5.*(putin(i,j+1,k)-putin(i,j-2,k))+(putin(i,j+2,k)-putin(i,j-3,k)))&
                )* dyi5 &
              + (1./rhobh(k))*( &
              (rhobh(k) * putin(i,j,k) + rhobh(k+1) * putin(i,j,k+1) )*(w0(i,j,k) + w0(i,j,k+1)) &
                -(rhobh(k) * putin(i,j,k) + rhobh(k-1) * putin(i,j,k-1) )*(w0(i,j,k) + w0(i,j,k-1)) &
                )/ (4. * dzh(k)) &
                )
      end do
    end do
  end do

end subroutine advecw_52
