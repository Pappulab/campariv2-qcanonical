!--------------------------------------------------------------------------!
! LICENSE INFO:                                                            !
!--------------------------------------------------------------------------!
!    This file is part of CAMPARI.                                         !
!                                                                          !
!    Version 2.0                                                           !
!                                                                          !
!    Copyright (C) 2014, The CAMPARI development team (current and former  !
!                        contributors)                                     !
!                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang !
!                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     !
!                        Nicholas Lyle, Nicolas Bloechliger                !
!                                                                          !
!    Website: http://sourceforge.net/projects/campari/                     !
!                                                                          !
!    CAMPARI is free software: you can redistribute it and/or modify       !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!    CAMPARI is distributed in the hope that it will be useful,            !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
!    GNU General Public License for more details.                          !
!                                                                          !
!    You should have received a copy of the GNU General Public License     !
!    along with CAMPARI.  If not, see <http://www.gnu.org/licenses/>.      !
!--------------------------------------------------------------------------!
! AUTHORSHIP INFO:                                                         !
!--------------------------------------------------------------------------!
!                                                                          !
! MAIN AUTHOR:   Andreas Vitalis                                           !
! CONTRIBUTIONS: Xiaoling Wang, Rohit Pappu, Adam Steffen                  !
!                                                                          !
!--------------------------------------------------------------------------!
!
!
#include "macros.i"
!
!------------------------------------------------------------------------------
!
! ###########################################################
! #                                                         #
! #  this file contains auxiliary math routines, usually LA #
! #                                                         #
! ###########################################################
!
!--------------------------------------------------------------------------
!
! this function computes different types of means for two real numbers >= 0.0 (the L_x norms, where x is specified by i_genmu)
! for i_genmu <= 0, 0.0 is returned as result if either or both numbers are == 0.0
! supplying negative values for r1 and r2 runs counter to the intention of the fxn, and silently will return the arithmetic mean
! specially coded cases are 1, 0, -1, and 2; in addition, -999 and 999 serve as specific selectors for simple min/max
!
function genmu(r1,r2,i_genmu)
!
  implicit none
!
  RTYPE genmu,r1,r2
  integer i_genmu
!
  if ((r1.lt.0.0).OR.(r2.lt.0.0)) then
    genmu = 0.5*(r1 + r2)
  else if (i_genmu.gt.0) then
    if (i_genmu.eq.1) then
      genmu = 0.5*(r1 + r2)
    else if (i_genmu.eq.2) then
      genmu = sqrt(0.5*(r1*r1 + r2*r2))
    else if (i_genmu.eq.999) then
      genmu = max(r1,r2)
    else
      genmu = (0.5*(r1**(1.*i_genmu) + r2**(1.*i_genmu)))**(1./(1.*i_genmu))
    end if
  else
    if ((r1.le.0.0).OR.(r2.le.0.0)) then
      genmu = 0.0
    else if (i_genmu.eq.0) then
      genmu = sqrt(r1*r2)
    else if (i_genmu.eq.-1) then
      genmu = 1.0/(0.5*(1./r1 + 1./r2))
    else if (i_genmu.eq.-999) then
      genmu = min(r1,r2)
    else
      genmu = (0.5*(r1**(1.*i_genmu) + r2**(1.*i_genmu)))**(1./(1.*i_genmu))
    end if
  end if
!
end
!
!---------------------------------------------------------------
!
! makes cross product and its norm 
!
subroutine crossprod(v1,v2,vo,norm)
!
  implicit none
!
  RTYPE v1(3),v2(3),vo(3),norm
!
  vo(1) = v1(2)*v2(3) - v1(3)*v2(2)
  vo(2) = v1(3)*v2(1) - v1(1)*v2(3)
  vo(3) = v1(1)*v2(2) - v1(2)*v2(1)
  norm = sqrt(vo(1)*vo(1) + vo(2)*vo(2) + vo(3)*vo(3))
!
end
!
!---------------------------------------------------------------
!
! makes cross product and its squared norm
!
subroutine crossprod2(v1,v2,vo,norm2)
!
  implicit none
!
  RTYPE v1(3),v2(3),vo(3),norm2
!
  vo(1) = v1(2)*v2(3) - v1(3)*v2(2)
  vo(2) = v1(3)*v2(1) - v1(1)*v2(3)
  vo(3) = v1(1)*v2(2) - v1(2)*v2(1)
  norm2 = vo(1)*vo(1) + vo(2)*vo(2) + vo(3)*vo(3)
!
end
!
!---------------------------------------------------------------
!
! makes just cross product
!
subroutine crossprod3(v1,v2,vo)
!
  implicit none
!
  RTYPE v1(3),v2(3),vo(3)
!
  vo(1) = v1(2)*v2(3) - v1(3)*v2(2)
  vo(2) = v1(3)*v2(1) - v1(1)*v2(3)
  vo(3) = v1(1)*v2(2) - v1(2)*v2(1)
!
end
!
!-----------------------------------------------------------------------
!
! this routine inverts a matrix by means of Cholesky decomposition
! and simple linear solving
!
subroutine invmat(b,d,a)
!
  use iounit
!
  implicit none
!
  integer i,j,k,d
  integer icol,irow
  integer ipivot(d)
  integer indxc(d)
  integer indxr(d)
  RTYPE big,temp,pivot
  RTYPE a(d,d),b(d,d)
!
  do i=1,d
    do j=1,d
      a(i,j) = b(i,j)
    end do
  end do
!
  do i=1,d
    ipivot(i) = 0
  end do
  do i=1,d
    big = 0.0d0
    do j=1,d
      if (ipivot(j).ne.1) then
        do k=1,d
          if (ipivot(k).eq.0) then
            if (abs(a(j,k)).ge.big) then
              big = abs(a(j,k))
              irow = j
              icol = k
            end if
          else if (ipivot(k).gt.1) then
            write (ilog,*) 'Fatal. Encountered singular matrix in in&
 &vmat(...).' 
            call fexit()
          end if
        end do
      end if
    end do
    ipivot(icol) = ipivot(icol) + 1
    if (irow.ne.icol) then
      do j=1,d
        temp = a(irow,j)
        a(irow,j) = a(icol,j)
        a(icol,j) = temp
      end do
    end if
    indxr(i) = irow
    indxc(i) = icol
    if (a(icol,icol).eq.0.0) then
      write (ilog,*) 'Fatal. Encountered singular matrix in invmat(.&
 &..).' 
      call fexit()
    end if
    pivot = a(icol,icol)
    a(icol,icol) = 1.0d0
    do j=1,d
      a(icol,j) = a(icol,j) / pivot
    end do
    do j=1,d
      if (j.ne.icol) then
        temp = a(j,icol)
        a(j,icol) = 0.0d0
        do k=1,d
          a(j,k) = a(j,k) - a(icol,k)*temp
        end do
      end if
    end do
  end do
  do i=d,1,-1
    if (indxr(i).ne.indxc(i)) then
      do k=1,d
        temp = a(k,indxr(i))
        a(k,indxr(i)) = a(k,indxc(i))
        a(k,indxc(i)) = temp
      end do
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! a simple matrix multiplication for 2 3x3 matrices
!
subroutine max3max3(m1,m2,mp)
!
  implicit none
!
  RTYPE mp(3,3),m1(3,3),m2(3,3)
!
! mp = m1*m2
! note: m1*m2 ~= m2*m1
! 
  mp(1,1) = m1(1,1)*m2(1,1) + m1(1,2)*m2(2,1) + m1(1,3)*m2(3,1)
  mp(1,2) = m1(1,1)*m2(1,2) + m1(1,2)*m2(2,2) + m1(1,3)*m2(3,2)
  mp(1,3) = m1(1,1)*m2(1,3) + m1(1,2)*m2(2,3) + m1(1,3)*m2(3,3)
!
  mp(2,1) = m1(2,1)*m2(1,1) + m1(2,2)*m2(2,1) + m1(2,3)*m2(3,1)
  mp(2,2) = m1(2,1)*m2(1,2) + m1(2,2)*m2(2,2) + m1(2,3)*m2(3,2)
  mp(2,3) = m1(2,1)*m2(1,3) + m1(2,2)*m2(2,3) + m1(2,3)*m2(3,3)
!
  mp(3,1) = m1(3,1)*m2(1,1) + m1(3,2)*m2(2,1) + m1(3,3)*m2(3,1)
  mp(3,2) = m1(3,1)*m2(1,2) + m1(3,2)*m2(2,2) + m1(3,3)*m2(3,2)
  mp(3,3) = m1(3,1)*m2(1,3) + m1(3,2)*m2(2,3) + m1(3,3)*m2(3,3)
!
end
!  
!-----------------------------------------------------------------------
!
! an even simpler multiplication routine for 3x3 times 3x1
!
subroutine max3vet(m,v,vp)
!
  implicit none
!
  RTYPE vp(3),m(3,3),v(3)
!
! vp = m*v
! 
  vp(1) = m(1,1)*v(1) + m(1,2)*v(2) + m(1,3)*v(3)
  vp(2) = m(2,1)*v(1) + m(2,2)*v(2) + m(2,3)*v(3)
  vp(3) = m(3,1)*v(1) + m(3,2)*v(2) + m(3,3)*v(3)
!
end
!
!-------------------------------------------------------------
!
! a modified matrix multiplication (doubles) routine which expects AM flipped
! (such that dot product can be maximally efficient)
!
subroutine matmul2(sz1,sz2,sz3,AM,BM,CM)
!
  implicit none
!  
  integer sz1,sz2,sz3,i,j
  RTYPE AM(sz2,sz1)
  RTYPE BM(sz2,sz3)
  RTYPE CM(sz1,sz3)
!
  do i=1,sz1
    do j=1,sz3
      CM(i,j) = dot_product(AM(1:sz2,i),BM(1:sz2,j))
    end do
  end do
!
end
!  
!-----------------------------------------------------------------------
! 
! the core routine which does not perform sorting of eigenvalues/vectors
!
subroutine mat_diag2(mdim,inm,evals,evs)
!
  use iounit
!
  implicit none
!
  integer i,j,ii,kk,maxiters,mdim,niters
  RTYPE netsum,crit,f1,f2,f3,f4,f5,f6,f7
  RTYPE inm(mdim,mdim),evals(mdim),evs(mdim,mdim)
  RTYPE dum1(mdim),dum2(mdim)
!
! some checking and initialization
!
  maxiters = 1000
  niters = 0
!
  do ii=1,mdim
    do kk=1,mdim
      evs(ii,kk) = 0.0d0
    end do
    evs(ii,ii) = 1.0d0
  end do
  do ii=1,mdim
    dum1(ii)=inm(ii,ii)
    evals(ii)=dum1(ii)
    dum2(ii) = 0.0d0
  end do
!
  do i=1,maxiters
!   first check for all zeros in relevant half
    netsum = 0.0
    do ii=1,mdim-1
      netsum = netsum + sum(abs(inm(ii,(ii+1):mdim)))
    end do
!
    if (netsum.le.0.0) exit
!
    if (i.lt.4) then
      crit = 0.2*netsum/(mdim**2)
    else
      crit = 0.0d0
    end if

    do ii=1,mdim-1
      do kk=ii+1,mdim
        f1 = 100.0*abs(inm(ii,kk))
        if ((i.gt.4).AND.((abs(evals(ii))+f1).eq.abs(evals(ii)))&
 &          .AND.((abs(evals(kk))+f1).eq.abs(evals(kk)))) then
          inm(ii,kk) = 0.0d0
        else if (abs(inm(ii,kk)).gt.crit) then
          f2 = evals(kk) - evals(ii)
          if ((abs(f2)+f1).eq.abs(f2)) then
            f3 = inm(ii,kk)/f2
          else
            f4 = 0.5*f2/inm(ii,kk)
            f3 = 1.0/(abs(f4)+sqrt(1.0+f4**2))
            if (f4.lt.0.0)  f3 = -f3
          end if
          f5 = 1.0/sqrt(1.0+f3**2)
          f7 = f3*f5
          f6 = f7/(1.0+f5)
          f2 = f3 * inm(ii,kk)
          dum2(ii) = dum2(ii) - f2
          dum2(kk) = dum2(kk) + f2
          evals(ii) = evals(ii) - f2
          evals(kk) = evals(kk) + f2
          inm(ii,kk) = 0.0d0
          do j=1,ii-1
            f1 = inm(j,ii)
            f2 = inm(j,kk)
            inm(j,ii) = f1 - f7*(f2+f1*f6)
            inm(j,kk) = f2 + f7*(f1-f2*f6)
          end do
          do j=ii+1,kk-1
            f1 = inm(ii,j)
            f2 = inm(j,kk)
            inm(ii,j) = f1 - f7*(f2+f1*f6)
            inm(j,kk) = f2 + f7*(f1-f2*f6)
          end do
          do j=kk+1,mdim
            f1 = inm(ii,j)
            f2 = inm(kk,j)
            inm(ii,j) = f1 - f7*(f2+f1*f6)
            inm(kk,j) = f2 + f7*(f1-f2*f6)
          end do
          do j=1,mdim
            f1 = evs(j,ii)
            f2 = evs(j,kk)
            evs(j,ii) = f1 - f7*(f2+f1*f6)
            evs(j,kk) = f2 + f7*(f1-f2*f6)
          end do
          niters = niters + 1
        end if
      end do
    end do
    do ii=1,mdim
      dum1(ii) = dum1(ii) + dum2(ii)
      evals(ii) = dum1(ii)
      dum2(ii) = 0.0d0
    end do
  end do
!
  if (niters.eq.maxiters) then
    write(ilog,*) 'Warning. Matrix diagonalization not converged in &
 &mat_diag(...). Distrust any results.'
  end if
!
end
!
!-----------------------------------------------------------------------
!
! a wrapper which adds sort functionality
!
subroutine mat_diag(mdim,inm,evals,evs)
!
  use iounit
!
  implicit none
!
  integer i,k,j,mdim
  RTYPE crit
  RTYPE inm(mdim,mdim),evals(mdim),evs(mdim,mdim)
!
  call mat_diag2(mdim,inm,evals,evs)
!
  do i=1,mdim-1
    k = i
    crit = evals(i)
    do j=i+1,mdim
      if (evals(j).lt.crit) then
        k = j
        crit = evals(j)
      end if
    end do
    if (k.ne.i) then
      evals(k) = evals(i)
      evals(i) = crit
      do j=1,mdim
        crit = evs(j,i)
        evs(j,i) = evs(j,k)
        evs(j,k) = crit
      end do
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! matrix determinant (a of dimension n)
!
subroutine dtrm (a,n,d,indx)
! !
! ! Subroutine for evaluating the determinant of a matrix using 
! ! the partial-pivoting Gaussian elimination scheme.
! ! Copyright (c) Tao Pang 2001.
! !
   use iounit
! !
   implicit none
! !
   integer, INTENT (IN) :: n
   integer :: i,j,msgn
   integer, INTENT (OUT), DIMENSION (n) :: indx
   RTYPE, INTENT (OUT) :: d
   RTYPE, INTENT (INOUT), DIMENSION (n,n) :: A
! !
   call elgs(a,n,indx)
! !
   d = 1.0
   do i = 1, n
     d = d*a(indx(i),i)
   end do
! !
   msgn = 1
   do i = 1, n
     do while (indx(i).ne.i)
       msgn = -msgn
       j = indx(i)
       indx(i) = indx(j)
       indx(j) = j
       if(indx(i).eq.indx(j)) then
         write(ilog,*) 'Fatal. Determinant of bad matrix.'
         call fexit()
       end if
     end do
   end do
   d = msgn*d
end subroutine dtrm
! !
!-----------------------------------------------------------------------------
! !
subroutine elgs (a,n,indx)
! !
! ! Subroutine to perform the partial-pivoting Gaussian elimination.
! ! A(N,N) is the original matrix in the input and transformed matrix
! ! plus the pivoting element ratios below the diagonal in the output.
! ! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
! !
    implicit none
  !
    integer, INTENT (IN) :: n
    integer :: i,j,k,itmp
    integer, INTENT (OUT), DIMENSION (n) :: indx
    RTYPE c1,pi,pi1,pj
    RTYPE, INTENT (INOUT), DIMENSION (n,n) :: a
    RTYPE, DIMENSION (n) :: c
! !
! ! Initialize the index
! !
    
    do i = 1, n
      indx(i) = i
    end do
! !
! ! Find the rescaling factors, one from each row
! !
    do i = 1, n
      c1= 0.0
      do j = 1, n
        c1 = max(c1,abs(a(i,j)))
      end do
      c(i) = c1
    end do
! !
! ! Search the pivoting (largest) element from each column
! !
    do j = 1, n-1
      pi1 = 0.0
      do i = j, n
        pi = abs(a(indx(i),j))/c(indx(i))
        if (pi.gt.pi1) then
          pi1 = pi
          k   = i
        end if
      end do
! !
! ! Interchange the rows via INDX(N) to record pivoting order
! !
      itmp    = indx(j)
      indx(j) = indx(k) 
      indx(k) = itmp
      do i = j+1, n
        pj  = a(indx(i),j)/a(indx(j),j)
! !
! ! Record pivoting ratios below the diagonal
! !
        a(indx(i),j) = pj
! !
! ! Modify other elements accordingly
! !
        do k = j+1, n
          a(indx(i),k) = a(indx(i),k)-pj*a(indx(j),k)
        end do
      end do
    end do
! !
end subroutine elgs
!
!---------------------------------------------------------------------------------------------
!
! this subroutine finds the translation vector and quaternion for rotation of the set in olds to
! superpose optimally onto news
! all for the 3D case
!
subroutine align_3D(sz,olds,news,tvec,qrot,centroido,centroidn)
!
  use math
  use iounit
!
  implicit none
!
  integer sz,i
  RTYPE tvec(3),qrot(4),centroido(3),centroidn(3)
  RTYPE olds(3*sz),news(3*sz),mata(4,4),matb(4,4),matn(4,4),ntmp(4,4),evmat(4,4),ev(4)
!
  centroidn(:) = 0.0
  centroido(:) = 0.0
  do i=1,sz
    centroidn(1) = centroidn(1) + news(3*i-2)
    centroidn(2) = centroidn(2) + news(3*i-1)
    centroidn(3) = centroidn(3) + news(3*i)
    centroido(1) = centroido(1) + olds(3*i-2)
    centroido(2) = centroido(2) + olds(3*i-1)
    centroido(3) = centroido(3) + olds(3*i)
  end do
  centroidn(:) = centroidn(:)/(1.0*sz)
  centroido(:) = centroido(:)/(1.0*sz)
  tvec(:) = centroidn(:) - centroido(:)
  mata(:,:) = 0.0
  matb(:,:) = 0.0
  matn(:,:) = 0.0
  do i=1,sz
    mata(1,2) = -(olds(3*i-2)-centroido(1))
    mata(2,1) = olds(3*i-2)-centroido(1)
    mata(3,4) = olds(3*i-2)-centroido(1)
    mata(4,3) = -(olds(3*i-2)-centroido(1))
    mata(1,3) = -(olds(3*i-1)-centroido(2))
    mata(2,4) = -(olds(3*i-1)-centroido(2))
    mata(3,1) = olds(3*i-1)-centroido(2)
    mata(4,2) = olds(3*i-1)-centroido(2)
    mata(1,4) = -(olds(3*i)-centroido(3))
    mata(2,3) = olds(3*i)-centroido(3)
    mata(3,2) = -(olds(3*i)-centroido(3))
    mata(4,1) = olds(3*i)-centroido(3)
    matb(1,2) = -(news(3*i-2)-centroidn(1))
    matb(2,1) = news(3*i-2)-centroidn(1)
    matb(3,4) = -(news(3*i-2)-centroidn(1))
    matb(4,3) = news(3*i-2)-centroidn(1)
    matb(1,3) = -(news(3*i-1)-centroidn(2))
    matb(2,4) = news(3*i-1)-centroidn(2)
    matb(3,1) = news(3*i-1)-centroidn(2)
    matb(4,2) = -(news(3*i-1)-centroidn(2))
    matb(1,4) = -(news(3*i)-centroidn(3))
    matb(2,3) = -(news(3*i)-centroidn(3))
    matb(3,2) = news(3*i)-centroidn(3)
    matb(4,1) = news(3*i)-centroidn(3)
    ntmp(:,:) = matmul(transpose(mata),matb)
    matn(:,:) = matn(:,:) + ntmp(:,:)
  end do
  call mat_diag(4,matn,ev,evmat)
  qrot(:) = evmat(:,4)
!
end
!
!---------------------------------------------------------------------------------------------
!
! the same for the case where the centroid is meant to be a weighted average (center of mass for instance)
!
subroutine align_3D_wt(sz,olds,news,wts,tvec,qrot,centroido,centroidn)
!
  use math
  use iounit
!
  implicit none
!
  integer sz,i
  RTYPE tvec(3),qrot(4),centroido(3),centroidn(3)
  RTYPE olds(3*sz),news(3*sz),wts(sz),mata(4,4),matb(4,4),matn(4,4),ntmp(4,4),evmat(4,4),ev(4),netwt
!
  centroidn(:) = 0.0
  centroido(:) = 0.0
  do i=1,sz
    centroidn(1) = centroidn(1) + wts(i)*news(3*i-2)
    centroidn(2) = centroidn(2) + wts(i)*news(3*i-1)
    centroidn(3) = centroidn(3) + wts(i)*news(3*i)
    centroido(1) = centroido(1) + wts(i)*olds(3*i-2)
    centroido(2) = centroido(2) + wts(i)*olds(3*i-1)
    centroido(3) = centroido(3) + wts(i)*olds(3*i)
  end do
  netwt = sum(wts(1:sz))
  centroidn(:) = centroidn(:)/netwt
  centroido(:) = centroido(:)/netwt
  tvec(:) = centroidn(:) - centroido(:)
  mata(:,:) = 0.0
  matb(:,:) = 0.0
  matn(:,:) = 0.0
  do i=1,sz
    mata(1,2) = -(olds(3*i-2)-centroido(1))
    mata(2,1) = olds(3*i-2)-centroido(1)
    mata(3,4) = olds(3*i-2)-centroido(1)
    mata(4,3) = -(olds(3*i-2)-centroido(1))
    mata(1,3) = -(olds(3*i-1)-centroido(2))
    mata(2,4) = -(olds(3*i-1)-centroido(2))
    mata(3,1) = olds(3*i-1)-centroido(2)
    mata(4,2) = olds(3*i-1)-centroido(2)
    mata(1,4) = -(olds(3*i)-centroido(3))
    mata(2,3) = olds(3*i)-centroido(3)
    mata(3,2) = -(olds(3*i)-centroido(3))
    mata(4,1) = olds(3*i)-centroido(3)
    matb(1,2) = -(news(3*i-2)-centroidn(1))
    matb(2,1) = news(3*i-2)-centroidn(1)
    matb(3,4) = -(news(3*i-2)-centroidn(1))
    matb(4,3) = news(3*i-2)-centroidn(1)
    matb(1,3) = -(news(3*i-1)-centroidn(2))
    matb(2,4) = news(3*i-1)-centroidn(2)
    matb(3,1) = news(3*i-1)-centroidn(2)
    matb(4,2) = -(news(3*i-1)-centroidn(2))
    matb(1,4) = -(news(3*i)-centroidn(3))
    matb(2,3) = -(news(3*i)-centroidn(3))
    matb(3,2) = news(3*i)-centroidn(3)
    matb(4,1) = news(3*i)-centroidn(3)
    ntmp(:,:) = matmul(transpose(mata),matb)
    matn(:,:) = matn(:,:) + ntmp(:,:)
  end do
  call mat_diag(4,matn,ev,evmat)
  qrot(:) = evmat(:,4)
!
end
!
!------------------------------------------------------------------------------------------
!
subroutine centroidshf_3D(sz,olds,news)
!
  use math
  use iounit
!
  implicit none
!
  integer sz,i
  RTYPE centroido(3),centroidn(3),olds(3*sz),news(3*sz)
!
  centroidn(:) = 0.0
  centroido(:) = 0.0
  do i=1,sz
    centroidn(1) = centroidn(1) + news(3*i-2)
    centroidn(2) = centroidn(2) + news(3*i-1)
    centroidn(3) = centroidn(3) + news(3*i)
    centroido(1) = centroido(1) + olds(3*i-2)
    centroido(2) = centroido(2) + olds(3*i-1)
    centroido(3) = centroido(3) + olds(3*i)
  end do
  centroidn(:) = centroidn(:)/(1.0*sz)
  centroido(:) = centroido(:)/(1.0*sz)
  do i=1,sz
    news(3*i-2:3*i) = news(3*i-2:3*i) - centroidn(1:3)
    olds(3*i-2:3*i) = olds(3*i-2:3*i) - centroido(1:3)
  end do
!
end
!
!----------------------------------------------------------------------------------------
!
! provided single derivatives, cross-second derivative, and function values on a 2D grid of
! dimx,dimy, compute - for each grid point - the 16 parameters for the bicubic spline into
! pararray(1:dimx,1:dimy,1:16)
! periodicity can be chosen via wrapit
!
subroutine bicubic_spline_set_params(dimx,dimy,binx,biny,val,derivx,derivy,derivxy,pararray,wrapit)
!
  implicit none
!
  integer i,j,iip,jjp,dimx,dimy
  RTYPE binx,biny,binxy
  RTYPE val(dimx,dimy),derivx(dimx,dimy),derivy(dimx,dimy),derivxy(dimx,dimy)
  RTYPE pararray(dimx,dimy,16),grvals(4,4)
  logical wrapit
!
  binxy = binx*biny
  do i=1,dimx
    do j=1,dimy
      iip = i + 1
      if (wrapit.EQV..true.) then
        if (iip.gt.dimx) iip = iip - dimx
      else
        iip = i
      end if
      jjp = j + 1
      if (wrapit.EQV..true.) then
        if (jjp.gt.dimy) jjp = jjp - dimy
      else
        jjp = j
      end if
      grvals(1,2) = binx*derivx(i,j)
      grvals(2,2) = binx*derivx(iip,j)
      grvals(3,2) = binx*derivx(i,jjp)
      grvals(4,2) = binx*derivx(iip,jjp)
      grvals(1,3) = biny*derivy(i,j)
      grvals(2,3) = biny*derivy(iip,j)
      grvals(3,3) = biny*derivy(i,jjp)
      grvals(4,3) = biny*derivy(iip,jjp)
      grvals(1,4) = binxy*derivxy(i,j)
      grvals(2,4) = binxy*derivxy(iip,j)
      grvals(3,4) = binxy*derivxy(i,jjp)
      grvals(4,4) = binxy*derivxy(iip,jjp)
      grvals(1,1) = val(i,j)
      grvals(2,1) = val(iip,j)
      if (iip.eq.i) grvals(2,1) = val(i,j) + binx*derivx(i,j)
      grvals(3,1) = val(i,jjp)
      if (jjp.eq.j) grvals(3,1) = val(i,j) + biny*derivy(i,j)
      grvals(4,1) = val(iip,jjp)
      if ((iip.eq.i).AND.(jjp.eq.j)) then
        grvals(4,1) = val(i,j) + binx*derivx(i,j) + biny*derivy(i,j)
      else if (iip.eq.i) then
        grvals(4,1) = val(i,jjp) + binx*derivx(i,jjp)
      else if (jjp.eq.j) then
        grvals(4,1) = val(iip,j) + biny*derivy(iip,j)
      end if
      pararray(i,j,1) = grvals(1,1)
      pararray(i,j,2) = grvals(1,2)
      pararray(i,j,3) = -3.0*grvals(1,1) + 3.0*grvals(2,1) - 2.0*grvals(1,2) - grvals(2,2)
      pararray(i,j,4) = 2.0*grvals(1,1) - 2.0*grvals(2,1) + grvals(1,2) + grvals(2,2)
      pararray(i,j,5) = grvals(1,3)
      pararray(i,j,6) = grvals(1,4)
      pararray(i,j,7) = -3.0*grvals(1,3) + 3.0*grvals(2,3) - 2.0*grvals(1,4) - grvals(2,4)
      pararray(i,j,8) = 2.0*grvals(1,3) - 2.0*grvals(2,3) + grvals(1,4) + grvals(2,4)
      pararray(i,j,9) = -3.0*grvals(1,1) + 3.0*grvals(3,1) - 2.0*grvals(1,3) - grvals(3,3)
      pararray(i,j,10) = -3.0*grvals(1,2) + 3.0*grvals(3,2) - 2.0*grvals(1,4) - grvals(3,4)
      pararray(i,j,11)= 9.0*grvals(1,1) - 9.0*grvals(2,1) - 9.0*grvals(3,1) + 9.0*grvals(4,1) &
 &                        + 6.0*grvals(1,2) + 3.0*grvals(2,2) - 6.0*grvals(3,2) - 3.0*grvals(4,2) &
 &                        + 6.0*grvals(1,3) - 6.0*grvals(2,3) + 3.0*grvals(3,3) - 3.0*grvals(4,3) &
 &                        + 4.0*grvals(1,4) + 2.0*grvals(2,4) + 2.0*grvals(3,4) + 1.0*grvals(4,4)
      pararray(i,j,12)=-6.0*grvals(1,1) + 6.0*grvals(2,1) + 6.0*grvals(3,1) - 6.0*grvals(4,1) &
 &                        - 3.0*grvals(1,2) - 3.0*grvals(2,2) + 3.0*grvals(3,2) + 3.0*grvals(4,2) &
 &                        - 4.0*grvals(1,3) + 4.0*grvals(2,3) - 2.0*grvals(3,3) + 2.0*grvals(4,3) &
 &                        - 2.0*grvals(1,4) - 2.0*grvals(2,4) - 1.0*grvals(3,4) - 1.0*grvals(4,4)
      pararray(i,j,13) = 2.0*grvals(1,1) - 2.0*grvals(3,1) + grvals(1,3) + grvals(3,3)
      pararray(i,j,14) = 2.0*grvals(1,2) - 2.0*grvals(3,2) + grvals(1,4) + grvals(3,4)
      pararray(i,j,15)=-6.0*grvals(1,1) + 6.0*grvals(2,1) + 6.0*grvals(3,1) - 6.0*grvals(4,1) &
 &                        - 4.0*grvals(1,2) - 2.0*grvals(2,2) + 4.0*grvals(3,2) + 2.0*grvals(4,2) &
 &                        - 3.0*grvals(1,3) + 3.0*grvals(2,3) - 3.0*grvals(3,3) + 3.0*grvals(4,3) &
 &                        - 2.0*grvals(1,4) - 1.0*grvals(2,4) - 2.0*grvals(3,4) - 1.0*grvals(4,4)
      pararray(i,j,16)= 4.0*grvals(1,1) - 4.0*grvals(2,1) - 4.0*grvals(3,1) + 4.0*grvals(4,1) &
 &                        + 2.0*grvals(1,2) + 2.0*grvals(2,2) - 2.0*grvals(3,2) - 2.0*grvals(4,2) &
 &                        + 2.0*grvals(1,3) - 2.0*grvals(2,3) + 2.0*grvals(3,3) - 2.0*grvals(4,3) &
 &                        + 1.0*grvals(1,4) + 1.0*grvals(2,4) + 1.0*grvals(3,4) + 1.0*grvals(4,4)
    end do
  end do
!
end
!
!-------------------------------------------------------------------------------------------------
!
! smooth a single data vector using B-splines of supplied order and a dedicated output vector
! no sanity checks, this assumes nonperiodic variables (i.e., not circular variables) and a nonperiodic vector at boundaries
! choice of smoothing at boundaries is mean-preserving
!
subroutine vsmooth(dvec,vl,bso,dveco)
!
  implicit none
!
  integer jj,kk,vl,bsoby2,bso,kk2
  RTYPE cbspl(bso),cbspld(bso),tval
  RTYPE dvec(vl),dveco(vl)
!
  if (mod(bso,2).eq.0) then
    bsoby2 = bso/2
    tval = 0.0
    call cardBspline(tval,bso,cbspl,cbspld)
!
    dveco(:) = 0.0
    do jj=1,bsoby2
      kk = vl+jj-bsoby2
      kk2 = bsoby2-jj+1
      dveco(1:kk) = dveco(1:kk) + cbspl(jj)*dvec(kk2:vl)
      if (jj.lt.bsoby2) then
        dveco(kk2:vl) = dveco(kk2:vl) + cbspl(bso-jj)*dvec(1:kk)
!       corrections
        dveco((kk+1):vl) = dveco((kk+1):vl) + cbspl(jj)*dvec((kk+1):vl)
        dveco(1:(kk2-1)) = dveco(1:(kk2-1)) + cbspl(bso-jj)*dvec(1:(kk2-1))
      end if
    end do
!
  else
!
    bsoby2 = (bso+1)/2
    tval = 0.5
    call cardBspline(tval,bso,cbspl,cbspld)
!
    dveco(:) = 0.0
    do jj=1,bsoby2
      kk = vl+jj-bsoby2
      kk2 = bsoby2-jj+1
      dveco(1:kk) = dveco(1:kk) + cbspl(jj)*dvec(kk2:vl)
      if (jj.lt.bsoby2) then
        dveco(kk2:vl) = dveco(kk2:vl) + cbspl(bso-jj+1)*dvec(1:kk)
!       corrections
        dveco((kk+1):vl) = dveco((kk+1):vl) + cbspl(jj)*dvec((kk+1):vl)
        dveco(1:(kk2-1)) = dveco(1:(kk2-1)) + cbspl(bso-jj+1)*dvec(1:(kk2-1))
      end if
    end do
!
  end if
!
end
!
!----------------------------------------------------------------------------------------------
!
! evaluate autocorrelation fxn at fixed lag time 
!
subroutine vacf_fixtau(dvec,vl,tauv,acfv)
!
  implicit none
!
  integer vl,tauv
  RTYPE dvec(vl),acfv,cvar,cmean,icst
!
  acfv = 1.0
  icst = 1.0/(1.0*vl)
! mean
  cmean = icst*sum(dvec(1:vl))
! use the unbiased estimator (vl - 1) for normalization of variance
  cvar = sum((dvec(:) - cmean)**2)/(1.0*(vl-1))
  if (cvar.gt.0.0) then
!   for the ACF at fixed tau, we assume a single process with fixed mean and variance
    icst = 1.0/(1.0*(vl-tauv))
    acfv = icst*sum((dvec(1:(vl-tauv))-cmean)*(dvec((tauv+1):vl)-cmean))/cvar
  end if
!
end
!
!---------------------------------------------------------------------------------------------------
!
! simple 1D hist with automatic data range limits and specified no of bins
! note that hvec is not initialized, which must happen before calling this fxn
!
subroutine vautohist(dvec,vl,nbins,hvec,hlbnd,hspac)
!
  implicit none
!
  integer i,k,vl,nbins,hvec(nbins)
  RTYPE dvec(vl),dmin,dmax,hbnds(nbins+1),hspac,hlbnd
!
  dmin = minval(dvec(1:vl))
  dmax = maxval(dvec(1:vl))
! 
! set up h
  hspac = (dmax-dmin)*(nbins+1.0)/(1.0*nbins*nbins)
  do i=1,nbins+1
    hbnds(i) = dmin + (i-1.5)*hspac
  end do
  hlbnd = hbnds(1)
!
  do i=1,vl
    k = ceiling((dvec(i)-hbnds(1))/hspac)
    hvec(k) = hvec(k) + 1 
  end do
!
end
!
!---------------------------------------------------------------------------------------------------
!
! return a list of minima of a 1D data vector using strict (> rather than >=) criteria over
! a fixed window size (embargo)
!
subroutine vminima_int(dvec,vl,nminis,minilst,ml,embargo)
!
  implicit none
!
  integer vl,dvec(vl),i,nminis,ml,minilst(ml),embargo
!
  do i=embargo+2,vl-embargo-1
    if ((minval(dvec((i-embargo-1):(i-1))).gt.dvec(i)).AND.(minval(dvec((i+1):(i+embargo+1))).gt.dvec(i))) then
      nminis = nminis + 1
      minilst(nminis) = i
    end if
  end do
!
end
!
!---------------------------------------------------------------------------------------------------------
!
! compute a windowed MSF: strictly O(N), window size is made into an uneven number if necessary
! compared to an explicit two-pass computation, the accuracy is superior over an explicit computation
! via <X^2> - <X>^2 (tested for uncentered data up to a window size 2x10^5 and overall size 4x10^5)
!
subroutine vmsf_window(dvec,vl,wsize,msfv)
!
  implicit none
!
  integer vl,i,wsize,wszby2,clonk
  RTYPE msfv(vl),dvec(vl),iwsz,iwsz2,curm,curm2,curv,newm,newm2,newv,redsu
!
  msfv(:) = 0.0
  if (wsize.le.1) then
    return
  end if
!
  if (mod(wsize,2).eq.0) wsize = wsize + 1
  wszby2 = wsize/2
  iwsz = 1.0/(1.0*wsize)
  iwsz2 = 1.0/(1.0*wszby2)
!
  msfv(1:(wszby2+1)) = dvec((wszby2+1):wsize)
  if ((vl-wszby2).ge.(wszby2+2)) then
    msfv((wszby2+2):(vl-wszby2)) = dvec((wsize+1):vl) - dvec(1:(vl-wsize))
  end if
  msfv((vl-wszby2+1):vl) = -dvec((vl-wsize+1):(vl-wszby2-1))
  curm = sum(dvec(1:wszby2))
  curm2 = iwsz2*curm
  curv = sum((dvec(1:wszby2)-curm2)**2)
  do i=1,vl
    clonk = min(vl+wszby2+1-i,min(i+wszby2,wsize))
    iwsz2 = 1.0/(1.0*clonk)
    newm = curm + msfv(i) ! the sum across the window
    newm2 = iwsz2*newm
    newv = curv
    redsu = curm
    if (i.gt.(wszby2+1)) then
      newv = newv - (dvec(i-wszby2-1)-curm2)**2
      redsu = redsu - dvec(i-wszby2-1)
    end if
    if (i.gt.(vl-wszby2)) clonk = clonk + 1
    newv = newv + (clonk-1)*(newm2**2-curm2**2)-2.0*(newm2-curm2)*redsu ! the shift in mean
    if (i.le.(vl-wszby2)) newv = newv + (dvec(i+wszby2)-newm2)**2
    msfv(i) = iwsz2*newv
    curv = newv
    curm = newm
    curm2 = newm2
  end do
!
end
!
!----------------------------------------------------------------------------------------------------------
!
