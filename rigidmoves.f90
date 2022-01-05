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
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! the strategy is the following:
! for each molecule in the system maintain c.o.m as well as a the eigenvectors of the
! gyration tensor
! WARNING: use centroid / gyration tensor OR com / inertia tensor
!          currently mislabeled centroid as com!!!
!          mass-weighted quantities are used in the dynamics routines
!          (like update_rigidm) -> comm, rgpcsm, ...
! this allows construction of isotropic translation and rotation moves
! rejected moves are cancelled as usual (explicit xyz restore)
!
!
!----------------------------------------------------------------------
!
subroutine mcmove_trans(imol,mode)
!
  use iounit
  use molecule
  use cutoffs
!
  implicit none
!
  integer mode,j,imol
!
  if ((mode.eq.0).OR.(mode.eq.1)) then
    call makeref_formol(imol)
    call transxyz(imol,mode)
    if (use_mcgrid.EQV..true.) then
      do j=rsmol(imol,1),rsmol(imol,2)
        call updateresgp(j)
      end do
    end if
  else if (mode.eq.2) then
    call restore_com(imol)
    call getref_formol(imol)
    if (use_mcgrid.EQV..true.) then
      do j=rsmol(imol,1),rsmol(imol,2)
        call updateresgp(j)
      end do
    end if
  else
    write(ilog,*) 'Fatal. Called mcmove_trans(...) with unknown mode&
 &. Offending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!----------------------------------------------------------------------
!
subroutine mcmove_rot(imol,mode)
!
  use iounit
  use molecule
  use cutoffs
!
  implicit none
!
  integer mode,j,imol
!
  if ((mode.eq.0).OR.(mode.eq.1)) then
    call makeref_formol(imol)
    call rotxyz(imol,mode)
    if (use_mcgrid.EQV..true.) then
      do j=rsmol(imol,1),rsmol(imol,2)
        call updateresgp(j)
      end do
    end if
  else if (mode.eq.2) then
!    call restore_gyrten(imol)
    call getref_formol(imol)
    if (use_mcgrid.EQV..true.) then
      do j=rsmol(imol,1),rsmol(imol,2)
        call updateresgp(j)
      end do
    end if
  else
    write(ilog,*) 'Fatal. Called mcmove_rot(...) with unknown mode .&
 &Offending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!----------------------------------------------------------------------
!
subroutine mcmove_rigid(imol,mode)
!
  use iounit
  use molecule
  use cutoffs
! Martin : debug only  
!
  implicit none
!
  integer mode,j,imol
!
  if ((mode.eq.0).OR.(mode.eq.1)) then
    call makeref_formol(imol)
    call rotxyz(imol,mode)
    call transxyz(imol,mode)

    if (use_mcgrid.EQV..true.) then
      do j=rsmol(imol,1),rsmol(imol,2)
        call updateresgp(j)
      end do
    end if
  else if (mode.eq.2) then
    call restore_com(imol)
    call getref_formol(imol)
    if (use_mcgrid.EQV..true.) then
      do j=rsmol(imol,1),rsmol(imol,2)
        call updateresgp(j)
      end do
    end if
  else
    write(ilog,*) 'Fatal. Called mcmove_rigid(...) with unknown mode&
 &. Offending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!----------------------------------------------------------------------
!
subroutine mcmove_clurigid(imols,mmol,mode)
!
  use iounit
  use molecule
  use cutoffs
  use movesets
!
  implicit none
!
  integer mode,j,imols(nmol),mmol,imol,i
!
  if ((mode.eq.0).OR.(mode.eq.1)) then
    do i=1,mmol
      imol = imols(i)
      call makeref_formol(imol)
      call makeref_poly(imol)
    end do
    call rotxyz_some(imols,mmol,mode)
    call transxyz_some(imols,mmol,mode)
    if (use_mcgrid.EQV..true.) then
      do i=1,mmol
        imol = imols(i)
        do j=rsmol(imol,1),rsmol(imol,2)
          call updateresgp(j)
        end do
      end do
    end if
  else if (mode.eq.2) then
    do i=1,mmol
      imol = imols(i)
      call getref_formol(imol)
      call getref_poly(imol)
    end do
    if (use_mcgrid.EQV..true.) then
      do i=1,mmol
        imol = imols(i)
        do j=rsmol(imol,1),rsmol(imol,2)
          call updateresgp(j)
        end do
      end do
    end if
  else
    write(ilog,*) 'Fatal. Called mcmove_clurigid(...) with unknown m&
 &ode. Offending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------
!
! a subroutine which does rigid body translation of a chosen molecule
!
subroutine transxyz(imol,mode)
!
  use atoms
  use molecule
  use movesets
!
  implicit none
!
  integer i,mode,k,imol
  RTYPE tvec(3),random
!
! sample freshly
  if (mode.eq.0) then
!   get a scaled, uniform, BC-compliant displacement (includes com-translation)
    if (random().lt.rigid_randfreq) then
      call randomize_trans(imol,tvec)
    else
      call trans_bound(imol,tvec)
    end if
!   backup current translation
    do i=1,3
      cur_trans(i) = tvec(i)
    end do
!
! re-sample
  else if (mode.eq.1) then
    do i=1,3
      tvec(i) = cur_trans(i)
    end do
    do i=1,3
      com(imol,i) = com(imol,i) + tvec(i)
    end do
!
  end if
!
! now translate atoms belonging to molecule
  do k=atmol(imol,1),atmol(imol,2)
    x(k) = x(k) + tvec(1)
    y(k) = y(k) + tvec(2)
    z(k) = z(k) + tvec(3)
  end do
!
end
!
!---------------------------------------------------------------------------
!
! a subroutine which does rigid body translation of a chosen molecule
!
subroutine transxyz_some(imols,mmol,mode)
!
  use atoms
  use molecule
  use movesets
!
  implicit none
!
  integer j,mode,k,imols(nmol),mmol,imol
  RTYPE tvec(3),random,tlen,tvec2(3),comn(3)
!
! sample freshly
  if (mode.eq.0) then
!
!   first the rigid properties of the cluster need to be determined
!    cur_clcom(1) = 0.0d0
!    cur_clcom(2) = 0.0d0
!    cur_clcom(3) = 0.0d0
!    normz = 0.0
!    do j=1,mmol
!      imol = imols(j)
!      cur_clcom(:) = cur_clcom(:) + 
! &              dble(atmol(imol,2)-atmol(imol,1)+1)*com(imol,:)
!      normz = normz + dble(atmol(imol,2)-atmol(imol,1)+1)
!    end do
!    cur_clcom(:) = cur_clcom(:) / normz
!   get a scaled, uniform, BC-compliant displacement (includes com-translation)
    if (random().lt.rigid_randfreq) then
      comn(:) = cur_clcom(:)
      call randomize_trans2(comn,tvec)
    else
      call ranvec(tvec)
      tlen = random()*trans_stepsz
      tvec(:) = tlen*tvec(:)
    end if
    do j=1,mmol
      tvec2(:) = tvec(:)
      imol = imols(j)
      call trans_bound4(imol,tvec2)
      com(imol,:) = com(imol,:) + tvec2(:)
      cur_clurbt(j,:) = tvec2(:)
    end do
!   backup current translation
    cur_trans(:) = tvec(:)
!
! re-sample
  else if (mode.eq.1) then
    do j=1,mmol
      imol = imols(j)
      com(imol,:) = com(imol,:) + cur_clurbt(j,:)
    end do
!
  end if
!
! now translate atoms belonging to molecule
  do j=1,mmol
    imol = imols(j)
    do k=atmol(imol,1),atmol(imol,2)
      x(k) = x(k) + cur_clurbt(j,1)
      y(k) = y(k) + cur_clurbt(j,2)
      z(k) = z(k) + cur_clurbt(j,3)
    end do
  end do
!
end
!
!---------------------------------------------------------------------------
!
! a subroutine which does rigid body translation of a chosen molecule (c.o.m)
!
subroutine transxyzm(imol,mode,shiftv)
!
  use atoms
  use molecule
  use movesets
  use iounit
!
  implicit none
!
  integer i,mode,k,imol
  RTYPE tvec(3),shiftv(3)
!
! sample freshly
  if (mode.eq.0) then
!   get a scaled, uniform, BC-compliant displacement (includes com-translation)
    write(ilog,*) 'Fatal.'
!
! re-sample
  else if (mode.eq.1) then
    do i=1,3
      tvec(i) = cur_trans(i)
    end do
    shiftv(:) = 0.0
    call trans_bound2(imol,tvec,shiftv)
    do i=1,3
      comm(imol,i) = comm(imol,i) + tvec(i)
    end do
!
  end if
!
! now translate atoms belonging to molecule
  do k=atmol(imol,1),atmol(imol,2)
    x(k) = x(k) + tvec(1)
    y(k) = y(k) + tvec(2)
    z(k) = z(k) + tvec(3)
  end do
!
end
!
!---------------------------------------------------------------------------
!
! a subroutine to restore the com
! this seems a bit inconsistent (could be handled inside transxyz()), but
! is necessary since xyz-backup and restoration are handled outside of
! the actual sampler (since in general that is faster than the reverse transform)
!
subroutine restore_com(imol)
!
  use molecule
  use movesets
!
  implicit none
!
  integer i,imol
!
  do i=1,3
    com(imol,i) = com(imol,i) - cur_trans(i)
  end do
end
!
!---------------------------------------------------------------------------
!
! a subroutine to restore the gyration tensor principal axes
! this seems a bit inconsistent (could be handled inside rotxyz()), but
! is necessary since xyz-backup and restoration are handled outside of
! the actual sampler (since in general that is faster than the reverse transform)
!
subroutine restore_gyrten(imol)
!
  use molecule
  use movesets
!
  implicit none
!
  integer i,imol,j
  RTYPE qrot(3,4),q1(4),q2(4),q3(4),qp(4),qpp(4)
  RTYPE dum(3),dum2(3)
!
  if ((atmol(imol,2)-atmol(imol,1)).eq.0) return
!
  do i=1,3
!   inverse transform: principal gyration axis i
    qrot(i,1) = cos(-cur_rot(i)/2.0)
    do j=2,4
      qrot(i,j) = rgpcs(imol,i,j-1)*sin(-cur_rot(i)/2.0)
    end do
  end do
!
  do j=1,4
    q1(j) = qrot(1,j)
    q2(j) = qrot(2,j)
    q3(j) = qrot(3,j)
  end do
!
! now get the total transformation through quat. multiplication
! note that the quat-product is non-commutative
  call quat_product(q3,q2,qp)
  call quat_product(qp,q1,qpp)
!
! now apply the inv. transformation only to the principal axes (ref. frame)
  do i=1,3
    do j=1,3
      dum(j) = rgpcs(imol,i,j)
    end do
    call quat_conjugate(qpp,dum,dum2)
    do j=1,3
      rgpcs(imol,i,j) = dum2(j)
    end do
  end do
!
end
!
!---------------------------------------------------------------------------
!
! a subroutine which does rigid body rotation around gyration tensor PCs
! note that the quaternion product in this similar routines does not commute, but that
! this does not imply an orientional bias
! when using complete randomization, the mean position of any atom averaging a large number
! of rotations from the same(!) coordinates should indeed be the center of mass 
!
subroutine rotxyz(imol,mode)
!
  use molecule
  use math
  use movesets
  use atoms
!
  implicit none
!
  integer i,j,k,imol,mode
  RTYPE qrot(3,4),q1(4),q2(4),q3(4),qp(4),qpp(4)
  RTYPE random,rdr
!
  if ((atmol(imol,2)-atmol(imol,1)).eq.0) return
!
! setup the quaternions representing the transformation
!
  if (mode.eq.0) then
!   sample freshly
    rdr = random()
    do i=1,3
!     get a random rotational displacement angle
      if (rdr.lt.rigid_randfreq) then
        cur_rot(i) = random()*2.0*PI - PI
      else
        cur_rot(i) = (rot_stepsz/360.0)*(random()*2.0*PI - PI)
      end if
!     apply to principal gyration axis i
      qrot(i,1) = cos(cur_rot(i)/2.0)
      do j=2,4
        qrot(i,j) = rgpcs(imol,i,j-1)*sin(cur_rot(i)/2.0)
      end do
    end do
  else if (mode.eq.1) then
!   re-apply
    do i=1,3
      qrot(i,1) = cos(cur_rot(i)/2.0)
      do j=2,4
        qrot(i,j) = rgpcs(imol,i,j-1)*sin(cur_rot(i)/2.0)
      end do
    end do
  end if
!
  do j=1,4
    q1(j) = qrot(1,j)
    q2(j) = qrot(2,j)
    q3(j) = qrot(3,j)
  end do
!
! now get the total transformation through quat. multiplication
  call quat_product(q1,q2,qp)
  call quat_product(qp,q3,qpp)
!
! now apply the transformation
! note that this automatically takes care of the global vs. local frame issue
  do k=atmol(imol,1),atmol(imol,2)
    call quat_conjugate2(qpp,imol,k)
  end do
!
end
!
!---------------------------------------------------------------------------
!
! a subroutine which does rigid body rotation around molecule "cluster" gyration
! tensor axes
!
subroutine rotxyz_some(imols,mmol,mode)
!
  use molecule
  use iounit
  use math
  use movesets
  use atoms
  use mcsums
!
  implicit none
!
  integer i,j,k,imols(nmol),mmol,mode,imol
  RTYPE qrot(3,4),q1(4),q2(4),q3(4),qp(4),qpp(4),qinv(4)
  RTYPE random,rdr,normz,mmm,comn(3),shcom(3,mmol)
  RTYPE eval(3),evmat(3,3),rgten(3,3),sh(3)
!
! setup the quaternions representing the transformation
!
  if (mode.eq.0) then
!   sample freshly
!
!   first the rigid properties of the cluster need to be determined correctly observing BC!
!   initialize
    cur_clcom(1) = 0.0d0
    cur_clcom(2) = 0.0d0
    cur_clcom(3) = 0.0d0
    do i=1,3
      do j=1,3
        rgten(i,j) = 0.0
      end do
    end do
    normz = 0.0
!
!   determine cluster geometric center
    do j=1,mmol
      imol = imols(j)
!     first mol. in the cluster is our reference
      call shift_bound(imols(1),imol,sh)
      shcom(:,j) = sh(:)
!      write(*,*) j,imol,shcom(:,j)
      cur_clcom(:) = cur_clcom(:) + &
 &              dble(atmol(imol,2)-atmol(imol,1)+1)*shcom(:,j)
      normz = normz + dble(atmol(imol,2)-atmol(imol,1)+1)
    end do
    cur_clcom(:) = cur_clcom(:) / normz
!
!   determine cluster gyration tensor, note we're translating xyz into
!   the right frame here as well!
    do j=1,mmol
      imol = imols(j)
      sh(:) = shcom(:,j) - com(imol,:)
!      if (sum(sh).gt.0.1) write(*,*) 'hit'
      do i=atmol(imol,1),atmol(imol,2)
        x(i) = x(i) + sh(1)
        y(i) = y(i) + sh(2)
        z(i) = z(i) + sh(3)
        rgten(1,1) = rgten(1,1) + (x(i)-cur_clcom(1))**2
        rgten(1,2) = rgten(1,2) + (x(i)-cur_clcom(1))*&
 &                                (y(i)-cur_clcom(2))
        rgten(1,3) = rgten(1,3) + (x(i)-cur_clcom(1))*&
 &                                (z(i)-cur_clcom(3))
        rgten(2,2) = rgten(2,2) + (y(i)-cur_clcom(2))**2
        rgten(2,3) = rgten(2,3) + (y(i)-cur_clcom(2))*&
 &                                (z(i)-cur_clcom(3))
        rgten(3,3) = rgten(3,3) + (z(i)-cur_clcom(3))**2
      end do
    end do
    rgten(2,1) = rgten(1,2)
    rgten(3,1) = rgten(1,3)
    rgten(3,2) = rgten(2,3)
    do i=1,3
      do j=1,3
        rgten(i,j) = rgten(i,j)/normz
      end do
    end do
!
!   diagonalize cluster gyration tensor and store principal components
    call mat_diag(3,rgten,eval,evmat)
    do i=1,3
      do j=1,3
        cur_clrgpcs(i,j) = evmat(j,i)
      end do
    end do
!
    rdr = random()
    do i=1,3
!     get a random rotational displacement angle
      if (rdr.lt.rigid_randfreq) then
        cur_rot(i) = random()*2.0*PI - PI
      else
        cur_rot(i) = (rot_stepsz/360.0)*(random()*2.0*PI - PI)
      end if
!     apply to principal gyration axis i
      qrot(i,1) = cos(cur_rot(i)/2.0)
      do j=2,4
        qrot(i,j) = cur_clrgpcs(i,j-1)*sin(cur_rot(i)/2.0)
      end do
    end do
  else if (mode.eq.1) then
!   re-apply
    do j=1,mmol
      imol = imols(j)
!     first mol. in the cluster is our reference
      call shift_bound(imols(1),imol,sh)
      shcom(:,j) = sh(:)
      sh(:) = shcom(:,j) - com(imol,:)
      do i=atmol(imol,1),atmol(imol,2)
        x(i) = x(i) + sh(1)
        y(i) = y(i) + sh(2)
        z(i) = z(i) + sh(3)
      end do
    end do
    do i=1,3
      qrot(i,1) = cos(cur_rot(i)/2.0)
      do j=2,4
        qrot(i,j) = cur_clrgpcs(i,j-1)*sin(cur_rot(i)/2.0)
      end do
    end do
  end if
!
  do j=1,4
    q1(j) = qrot(1,j)
    q2(j) = qrot(2,j)
    q3(j) = qrot(3,j)
  end do
!
! now get the total transformation through quat. multiplication
  call quat_product(q1,q2,qp)
  call quat_product(qp,q3,qpp)
! construct the multiplicative inverse
  mmm = sum(qpp(:)*qpp(:))
  if (mmm.le.0.0) then
    write(ilog,*) 'Fatal error. Got a 0-quat. in rotxyz_some(...).'
    call fexit()
  end if
  qinv(1) = qpp(1)/mmm
  do i=2,4
    qinv(i) = -qpp(i)/mmm
  end do
!
! now apply the transformation
! note that this automatically takes care of the global vs. local frame issue
  do j=1,mmol
    imol = imols(j)
    call quat_conjugate1(qpp,qinv,shcom(:,j),comn,cur_clcom)
    com(imol,:) = comn(:)
    do k=atmol(imol,1),atmol(imol,2)
      call quat_conjugate3(qpp,qinv,k,cur_clcom)
    end do
  end do
  do j=1,mmol
    imol = imols(j)
    call update_image2(imol)
  end do
!
end
!
!---------------------------------------------------------------------------
!
! a similar routine for TMD:
! here, the elements in cur_rot are derived from angular velocities, hence the
! construction of the unit quaternion has to avoid all ambiguities 
! (noncommutative product approach would be inappropriate)
! note that this requires the angular velocities to be referring to the space-fixed
! cardinal axes (inertial frame)
!
subroutine rotxyzm(imol,mode)
!
  use molecule
  use math
  use movesets
  use iounit
  use atoms
  use forces
!
  implicit none
!
  integer k,imol,mode
  RTYPE qpp(4)
!
  if ((atmol(imol,2)-atmol(imol,1)).eq.0) return
!
! setup the quaternions representing the transformation
!
  if (mode.eq.0) then
!   sample freshly
    write(ilog,*) 'Fatal.'
    call fexit()
  else if (mode.eq.1) then
!   (re-)apply
    qpp(2:4) = sin(cur_rot(1:3)/2.0)
!   make a unit quaternion - note that asin is safe for small rotations
    qpp(1) = cos(asin(sqrt(dot_product(qpp(2:4),qpp(2:4)))))
  end if
!
! now apply the transformation
! note that this automatically takes care of the global vs. local frame issue
  do k=atmol(imol,1),atmol(imol,2)
    call quat_conjugate2m(qpp,imol,k)
  end do
!
end
!
!---------------------------------------------------------------------------
!
! a subroutine which updates the eigenvectors of the gyration tensor
! after a successful rigid body rotation (MC)
!
subroutine update_gyrten(imol)
!
  use iounit
  use movesets
  use molecule
!
  implicit none
!
  integer i,j,imol
  RTYPE qrot(3,4),q1(4),q2(4),q3(4),qp(4),qpp(4)
  RTYPE dum(3),dum2(3)
!
  if ((atmol(imol,2)-atmol(imol,1)).eq.0) return
!
! setup the quaternions representing the transformation
! with the current rotation move
  do i=1,3
    qrot(i,1) = cos(cur_rot(i)/2.0)
    do j=2,4
      qrot(i,j) = rgpcs(imol,i,j-1)*sin(cur_rot(i)/2.0)
    end do
  end do
!
  do j=1,4
    q1(j) = qrot(1,j)
    q2(j) = qrot(2,j)
    q3(j) = qrot(3,j)
  end do
!
! now get the total transformation through quat. multiplication
  call quat_product(q1,q2,qp)
  call quat_product(qp,q3,qpp)
!
! now apply the transformation to the gyration tensor eigenvectors
! (gyrational reference frame expressed in Cartesian coordinates)
  do i=1,3
    do j=1,3
      dum(j) = rgpcs(imol,i,j)
    end do
    call quat_conjugate(qpp,dum,dum2)
    do j=1,3
      rgpcs(imol,i,j) = dum2(j)
    end do 
  end do
!
! 26   format(i2,1x,3(g12.4,1x))
!
end
!
!---------------------------------------------------------------------------
!
! a subroutine which updates the rigid body descriptors (com, gyration tensor pcs)
!
subroutine update_rigid(imol)
!
  use iounit
  use molecule
  use atoms
!
  implicit none
!
  integer i,j,imol
  RTYPE rgten(3,3),eval(3),evmat(3,3)
!
! find the centroid of the atomic coordinates
!
! initialize
  com(imol,1) = 0.0d0
  com(imol,2) = 0.0d0
  com(imol,3) = 0.0d0
  do i=1,3
    do j=1,3
      rgten(i,j) = 0.0
    end do
  end do
!
! update com
  do i=atmol(imol,1),atmol(imol,2)
    com(imol,1) = com(imol,1) + x(i)
    com(imol,2) = com(imol,2) + y(i)
    com(imol,3) = com(imol,3) + z(i)
  end do
  com(imol,1) = com(imol,1) / dble(atmol(imol,2)-atmol(imol,1)+1)
  com(imol,2) = com(imol,2) / dble(atmol(imol,2)-atmol(imol,1)+1)
  com(imol,3) = com(imol,3) / dble(atmol(imol,2)-atmol(imol,1)+1)
!  write(*,*) 'com: ',imol,(com(imol,j),j=1,3)
!
  if ((atmol(imol,2)-atmol(imol,1)).eq.0) return
!
! update gyration tensor
  do i=atmol(imol,1),atmol(imol,2)
    rgten(1,1) = rgten(1,1) + (x(i)-com(imol,1))**2
    rgten(1,2) = rgten(1,2) + (x(i)-com(imol,1))*(y(i)-com(imol,2))
    rgten(1,3) = rgten(1,3) + (x(i)-com(imol,1))*(z(i)-com(imol,3))
    rgten(2,2) = rgten(2,2) + (y(i)-com(imol,2))**2
    rgten(2,3) = rgten(2,3) + (y(i)-com(imol,2))*(z(i)-com(imol,3))
    rgten(3,3) = rgten(3,3) + (z(i)-com(imol,3))**2
  end do
  rgten(2,1) = rgten(1,2)
  rgten(3,1) = rgten(1,3)
  rgten(3,2) = rgten(2,3)
  do i=1,3
    do j=1,3
      rgten(i,j) = rgten(i,j)/(dble(atmol(imol,2)-atmol(imol,1)+1))
    end do
!    write(*,*) i,(rgten(i,j),j=1,3)
  end do
  rgv(imol) = sqrt(rgten(1,1)+rgten(2,2)+rgten(3,3))
!
! diagonalize gyration tensor and store principal components
  call mat_diag(3,rgten,eval,evmat)
  do i=1,3
    do j=1,3
      rgpcs(imol,i,j) = evmat(j,i)
    end do
    rgevs(imol,i) = eval(i)
  end do
!
end
!
!---------------------------------------------------------------------------
!
! primitive wrapper
!
subroutine update_rigid_mc_all()
!
  use molecule
!
  implicit none
!
  integer imol
!
  do imol=1,nmol
    call update_rigid(imol)
    call update_image2(imol)
  end do
!
end
!
!------------------------------------------------------------------------------
!
! a subroutine which updates the rigid body descriptors (com, gyration tensor pcs)
! this is the mass-weighted adaptation for dynamics purposes.
! note, however, that the principal axes are currently not in use (functionality is there),
! so the "principal" axes are simply kept as the cardinal axes (tensor is not diagonal and computed here)
!
subroutine update_rigidm(imol)
!
  use iounit
  use molecule
  use atoms
  use forces
!
  implicit none
!
  integer i,j,imol
  RTYPE rgtenm(3,3),eval(3),evmat(3,3)
  RTYPE rgtenm2(3,3),dum1(3),hlp1,hlp2,hlp3,hlp4,hlp5,hlp6
  logical used(3,2),use_principal_axes
!
  use_principal_axes = .false.
!
  if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
    comm(imol,1) = x(atmol(imol,1))
    comm(imol,2) = y(atmol(imol,1))
    comm(imol,3) = z(atmol(imol,1))
    return
  end if
!
! find the centroid of the atomic coordinates
!
! initialize
  rgtenm(:,:) = 0.0
!
  comm(imol,1) = sum(mass(atmol(imol,1):atmol(imol,2))*x(atmol(imol,1):atmol(imol,2)))
  comm(imol,2) = sum(mass(atmol(imol,1):atmol(imol,2))*y(atmol(imol,1):atmol(imol,2)))
  comm(imol,3) = sum(mass(atmol(imol,1):atmol(imol,2))*z(atmol(imol,1):atmol(imol,2)))
  comm(imol,:) = comm(imol,:)/molmass(moltypid(imol))
!
! update inertia tensor
  do i=atmol(imol,1),atmol(imol,2)
    hlp1 = x(i)-comm(imol,1)
    hlp2 = y(i)-comm(imol,2)
    hlp3 = z(i)-comm(imol,3)
    hlp4 = hlp1*hlp1
    hlp5 = hlp2*hlp2
    hlp6 = hlp3*hlp3
    rgtenm(1,1) = rgtenm(1,1) + mass(i)*(hlp5+hlp6)
    rgtenm(1,2) = rgtenm(1,2) - mass(i)*hlp1*hlp2
    rgtenm(1,3) = rgtenm(1,3) - mass(i)*hlp1*hlp3
    rgtenm(2,2) = rgtenm(2,2) + mass(i)*(hlp4+hlp6)
    rgtenm(2,3) = rgtenm(2,3) - mass(i)*hlp2*hlp3
    rgtenm(3,3) = rgtenm(3,3) + mass(i)*(hlp4+hlp5)
  end do
  rgtenm(2,1) = rgtenm(1,2)
  rgtenm(3,1) = rgtenm(1,3)
  rgtenm(3,2) = rgtenm(2,3)
  rgvm(imol) = sqrt(rgtenm(1,1)+rgtenm(2,2)+rgtenm(3,3))
!
  if (use_principal_axes.EQV..true.) then
!   diagonalize inertia tensor and store principal components
    rgtenm2(:,:) = rgtenm(:,:)
    call mat_diag2(3,rgtenm2,eval,evmat)
    if (abs(dot_product(rgpcsm(imol,1,:),rgpcsm(imol,1,:))-1.0).gt.1.0e-3) then
!     not yet set, take default order
      do i=1,3
        rgpcsm(imol,i,:) = evmat(:,i)
        rgevsm(imol,i) = eval(i)
      end do
    else
!     preserve order (find former eigenvector < 45 deg away, should be unique)
      dum1(1) = 0.5*sqrt(2.) ! cos(pi/4)
      dum1(3) = sqrt(3.)/3.
      used(:,:) = .false.
      do i=1,3
        do j=1,3
          dum1(2) = dot_product(rgpcsm(imol,i,:),evmat(:,j))
          if ((abs(dum1(2)).ge.dum1(1)).AND.(used(j,2).EQV..false.)) then
            used(j,2) = .true.
            rgpcsm(imol,i,:) = nint(dum1(2))*evmat(:,j)
            rgevsm(imol,i) = eval(j)
            used(i,1) = .true.
          end if
        end do
      end do
!     in case of strong changes, find former eigenvector < ~54 deg away, may not be unique
      do i=1,3
        if (used(i,1).EQV..false.) then
          do j=1,3
            dum1(2) = dot_product(rgpcsm(imol,i,:),evmat(:,j))
            if ((abs(dum1(2)).ge.dum1(3)).AND.(used(j,2).EQV..false.)) then
              used(j,2) = .true.
              rgpcsm(imol,i,:) = nint(dum1(2))*evmat(:,j)
              rgevsm(imol,i) = eval(j)
              used(i,1) = .true.
            end if
          end do
        end if
      end do
      do i=1,3
        if (used(i,1).EQV..false.) then
          write(ilog,*) 'Fatal.'
          call fexit()
        end if
      end do
    end if
  else
!   note that rgevsm will not be set
    rgpcsm(imol,:,:) = 0.0
    rgpcsm(imol,1,1) = 1.0
    rgpcsm(imol,2,2) = 1.0
    rgpcsm(imol,3,3) = 1.0
  end if
!
end
!
!---------------------------------------------------------------------------
!
! a subroutine which computes and writes the current com to an argument vector
! without touching the global
!
subroutine get_com(imol,comn)
!
  use molecule
  use atoms
!
  implicit none
!
  integer i,imol
  RTYPE comn(3)
!
! find the centroid of the atomic coordinates
!
! initialize
  comn(1) = 0.0d0
  comn(2) = 0.0d0
  comn(3) = 0.0d0
!
! compute com
  do i=atmol(imol,1),atmol(imol,2)
    comn(1) = comn(1) + x(i)
    comn(2) = comn(2) + y(i)
    comn(3) = comn(3) + z(i)
  end do
  comn(1) = comn(1) / dble(atmol(imol,2)-atmol(imol,1)+1)
  comn(2) = comn(2) / dble(atmol(imol,2)-atmol(imol,1)+1)
  comn(3) = comn(3) / dble(atmol(imol,2)-atmol(imol,1)+1)
!
end
!
!----------------------------------------------------------------------------
!
! a subroutine which fully randomizes a molecule's orientation
!
subroutine randomize_rot(imol)
!
  use iounit
  use movesets
  use molecule
  use math
!
  implicit none
!
  integer i,j,k,imol
  RTYPE qrot(3,4),q1(4),q2(4),q3(4),qp(4),qpp(4),random
!
  if ((atmol(imol,2)-atmol(imol,1)).eq.0) return
!
! setup the quaternions representing the transformation
!
! randomize
  do i=1,3
    cur_rot(i) = random()*2*PI - PI
!   apply to principal gyration axis i
    qrot(i,1) = cos(cur_rot(i)/2.0)
    do j=2,4
      qrot(i,j) = rgpcs(imol,i,j-1)*sin(cur_rot(i)/2.0)
    end do
  end do
!
  do j=1,4
    q1(j) = qrot(1,j)
    q2(j) = qrot(2,j)
    q3(j) = qrot(3,j)
  end do
!
! now get the total transformation through quat. multiplication
  call quat_product(q1,q2,qp)
  call quat_product(qp,q3,qpp)
!
! now apply the transformation
! note that this automatically takes care of the global vs. local frame issue
  do k=atmol(imol,1),atmol(imol,2)
    call quat_conjugate2(qpp,imol,k)
  end do
!
end
!
!---------------------------------------------------------------------------
!
! a subroutine for quaternion multiplication 
!
subroutine quat_product(q1,q2,prodq)
!
  implicit none
!
  integer i
  RTYPE q1(4),q2(4),prodq(4),v1(3),v2(3),vp(3),nm
!
! first the scalar component 
  prodq(1) = q1(1)*q2(1) - (q1(2)*q2(2)+q1(3)*q2(3)+q1(4)*q2(4))
!
! now the vectorial component
  do i=2,4
    v1(i-1) = q1(i)
    v2(i-1) = q2(i)
  end do
  call crossprod(v1,v2,vp,nm)
  do i=2,4
    prodq(i) = q1(1)*q2(i) + q2(1)*q1(i) + vp(i-1)
  end do
!
end
!
!---------------------------------------------------------------------------
!
! a subroutine for the conjugation of a vector by a quaternion 
!
subroutine quat_conjugate(q,vin,vout)
!
  use iounit
!
  implicit none
!
  integer i
  RTYPE q(4),qinv(4),qf(4),qp(4),qpp(4),vin(3),vout(3),nm
!
! construct the multiplicative inverse
  nm = q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4)
  if (nm.le.0.0) then
    write(ilog,*) 'Fatal error. Called quat_conjugate with 0-quat.'
    call fexit()
  end if
  qinv(1) = q(1)/nm
  do i=2,4
    qinv(i) = -q(i)/nm
  end do
!
! make the input vector into a fake quaternion
  qf(1) = 0.0
  qf(2) = vin(1)
  qf(3) = vin(2)
  qf(4) = vin(3)
!
! now do the conjugation
  call quat_product(q,qf,qp)
  call quat_product(qp,qinv,qpp)
!
! extract the rotated vector
  vout(1) = qpp(2)
  vout(2) = qpp(3)
  vout(3) = qpp(4)
!
end
!
!---------------------------------------------------------------------------
!
! the same with a reference point 
!
subroutine quat_conjugatei(q,vin,vout,refv)
!
  use iounit
!
  implicit none
!
  integer i
  RTYPE q(4),qinv(4),qf(4),qp(4),qpp(4),vin(3),vout(3),nm,refv(3)
!
! construct the multiplicative inverse
  nm = q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4)
  if (nm.le.0.0) then
    write(ilog,*) 'Fatal error. Called quat_conjugate with 0-quat.'
    call fexit()
  end if
  qinv(1) = q(1)/nm
  do i=2,4
    qinv(i) = -q(i)/nm
  end do
!
! make the input vector into a fake quaternion
  qf(1) = 0.0
  qf(2) = vin(1) - refv(1)
  qf(3) = vin(2) - refv(2)
  qf(4) = vin(3) - refv(3)
!
! now do the conjugation
  call quat_product(q,qf,qp)
  call quat_product(qp,qinv,qpp)
!
! extract the rotated vector
  vout(1) = qpp(2) + refv(1)
  vout(2) = qpp(3) + refv(2)
  vout(3) = qpp(4) + refv(3)
!
end
!
!---------------------------------------------------------------------------
!
! same as quat_conjugate, only takes the inverse as input and uses ref. point
!
subroutine quat_conjugate1(q,qinv,vin,vout,refv)
!
  use iounit
!
  implicit none
!
  RTYPE q(4),qinv(4),qf(4),qp(4),qpp(4),vin(3),vout(3),refv(3)
!
! make the input vector into a fake quaternion
  qf(1) = 0.0
  qf(2) = vin(1) - refv(1)
  qf(3) = vin(2) - refv(2)
  qf(4) = vin(3) - refv(3)
!
! now do the conjugation
  call quat_product(q,qf,qp)
  call quat_product(qp,qinv,qpp)
!
! extract the rotated vector
  vout(1) = qpp(2) + refv(1)
  vout(2) = qpp(3) + refv(2)
  vout(3) = qpp(4) + refv(3)
!
end
!
!---------------------------------------------------------------------------
!
! a similar routine that works on an atom number instead
!
subroutine quat_conjugate2(q,imol,ati)
!
  use iounit
  use atoms
  use molecule
!
  implicit none
!
  integer i,ati,imol
  RTYPE q(4),qinv(4),qf(4),qp(4),qpp(4),nm
!
! construct the multiplicative inverse
  nm = q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4)
  if (nm.le.0.0) then
    write(ilog,*) 'Fatal error. Called quat_conjugate2 with 0-quat.'
    call fexit()
  end if
  qinv(1) = q(1)/nm
  do i=2,4
    qinv(i) = -q(i)/nm
  end do
!
! make the input vector into a fake quaternion
  if (ati.eq.1) then
!    write(*,*) 'com:',com(imol,1),com(imol,2),com(imol,3)
!    call update_rigid(imol)
  end if
  qf(1) = 0.0
  qf(2) = x(ati)-com(imol,1)
  qf(3) = y(ati)-com(imol,2)
  qf(4) = z(ati)-com(imol,3)
  if (ati.eq.1) then
!    write(*,*) qf(2),qf(3),qf(4)
  end if
!
! now do the conjugation
  call quat_product(q,qf,qp)
  call quat_product(qp,qinv,qpp)
!
! extract the rotated vector
  x(ati) = qpp(2)+com(imol,1)
  y(ati) = qpp(3)+com(imol,2)
  z(ati) = qpp(4)+com(imol,3)
  if (ati.eq.1) then
!    write(*,*) qpp(2),qpp(3),qpp(4)
  end if
!
end
!---------------------------------------------------------------------------
!
! identical to qc2, but uses inertia tensor instead
!
subroutine quat_conjugate2m(q,imol,ati)
!
  use iounit
  use atoms
  use molecule
!
  implicit none
!
  integer i,ati,imol
  RTYPE q(4),qinv(4),qf(4),qp(4),qpp(4),nm
!
! construct the multiplicative inverse
  nm = q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4)
  if (nm.le.0.0) then
    write(ilog,*) 'Fatal error. Called quat_conjugate2m() with 0-qua&
 &ternion.'
    call fexit()
  end if
  qinv(1) = q(1)/nm
  do i=2,4
    qinv(i) = -q(i)/nm
  end do
!
! make the input vector into a fake quaternion
  if (ati.eq.1) then
!    write(*,*) 'com:',com(imol,1),com(imol,2),com(imol,3)
!    call update_rigid(imol)
  end if
  qf(1) = 0.0
  qf(2) = x(ati)-comm(imol,1)
  qf(3) = y(ati)-comm(imol,2)
  qf(4) = z(ati)-comm(imol,3)
  if (ati.eq.1) then
!    write(*,*) qf(2),qf(3),qf(4)
  end if
!
! now do the conjugation
  call quat_product(q,qf,qp)
  call quat_product(qp,qinv,qpp)
!
! extract the rotated vector
  x(ati) = qpp(2)+comm(imol,1)
  y(ati) = qpp(3)+comm(imol,2)
  z(ati) = qpp(4)+comm(imol,3)
  if (ati.eq.1) then
!    write(*,*) qpp(2),qpp(3),qpp(4)
  end if
!
end
!
!---------------------------------------------------------------------------
!
! a similar routine that works on an atom number instead
! and uses a supplied reference vector
!
subroutine quat_conjugate3(q,qinv,ati,refv)
!
  use iounit
  use atoms
!
  implicit none
!
  integer ati
  RTYPE q(4),qinv(4),qf(4),qp(4),qpp(4),refv(3)
!
! make the input vector into a fake quaternion
  qf(1) = 0.0
  qf(2) = x(ati)-refv(1)
  qf(3) = y(ati)-refv(2)
  qf(4) = z(ati)-refv(3)
!
! now do the conjugation
  call quat_product(q,qf,qp)
  call quat_product(qp,qinv,qpp)
!
! extract the rotated vector
  x(ati) = qpp(2)+refv(1)
  y(ati) = qpp(3)+refv(2)
  z(ati) = qpp(4)+refv(3)
!
end
!
!-----------------------------------------------------------------------------
!
! a similar routine that works on an atom number instead and uses a reference vector but no inverse
!
subroutine quat_conjugate4(q,ati,refv)
!
  use iounit
  use atoms
  use molecule
!
  implicit none
!
  integer i,ati
  RTYPE q(4),qinv(4),qf(4),qp(4),qpp(4),nm,refv(3)
!
! construct the multiplicative inverse
  nm = q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4)
  if (nm.le.0.0) then
    write(ilog,*) 'Fatal error. Called quat_conjugate2 with 0-quat.'
    call fexit()
  end if
  qinv(1) = q(1)/nm
  do i=2,4
    qinv(i) = -q(i)/nm
  end do
!
! make the input vector into a fake quaternion
  if (ati.eq.1) then
!    write(*,*) 'com:',com(imol,1),com(imol,2),com(imol,3)
!    call update_rigid(imol)
  end if
  qf(1) = 0.0
  qf(2) = x(ati)-refv(1)
  qf(3) = y(ati)-refv(2)
  qf(4) = z(ati)-refv(3)
  if (ati.eq.1) then
!    write(*,*) qf(2),qf(3),qf(4)
  end if
!
! now do the conjugation
  call quat_product(q,qf,qp)
  call quat_product(qp,qinv,qpp)
!
! extract the rotated vector
  x(ati) = qpp(2)+refv(1)
  y(ati) = qpp(3)+refv(2)
  z(ati) = qpp(4)+refv(3)
  if (ati.eq.1) then
!    write(*,*) qpp(2),qpp(3),qpp(4)
  end if
!
end
!
!--------------------------------------------------------------------------
!
! this subroutine will re-calculate the positions of a subset of atoms based on
! a torsional change, for which axis, increment, and ref. point are all provided
! note that this routine is free of any opt. sanity checks: e.g., it does not check
! whether ati,atj are part of imol
!
subroutine alignC_onetor(ati,atj,refv,axl,dw)
!
  use iounit
  use atoms
  use molecule
  use fyoc
  use math
  use sequen
  use polypep
  use zmatrix
!
  implicit none
!
  integer k,ati,atj
  RTYPE q(4),axl(3),dw,refv(3),mmm,qinv(4)
!
  if (atj.lt.ati) then
    write(ilog,*) 'Fatal. Bad atom input in alignC_onetor(...).'
    call fexit()
  end if
!
! first normalize rotation axis
!
  mmm = sqrt(sum(axl(:)*axl(:)))
  axl(:) = axl(:)/mmm
!
! now construct the quaternion that corresponds to the inverse omega move
  q(1) = cos(-dw/2.0)
  q(2:4) = axl(1:3)*sin(-dw/2.0)
! construct the multiplicative inverse
  mmm = sum(q(:)*q(:))
  if (mmm.le.0.0) then
    write(ilog,*) 'Fatal error. Got a 0-quat. in alignC_omega(...).'
    call fexit()
  end if
  qinv(1) = q(1)/mmm
  qinv(2:4) = -q(2:4)/mmm
!
! now apply the transformation
! note that this automatically takes care of the global vs. local frame issue
  do k=ati,atj
    call quat_conjugate3(q,qinv,k,refv)
  end do
!
end
!
!---------------------------------------------------------------------------
!
subroutine alignC_omega(rs)
!
  use iounit
  use atoms
  use molecule
  use fyoc
  use math
  use sequen
  use polypep
  use zmatrix
!
  implicit none
!
  integer imol,i,k,rs,ll,p1,p2
  RTYPE q(4),axl(3),dw,refv(3),mmm,qinv(4)
!
  imol = molofrs(rs)
!
! first get the proper rotation axis, and an appropriate reference point
!
  if ((seqtyp(rs).eq.41).OR.(seqtyp(rs).eq.42)) then !NMF,NMA
    p1 = at(rs)%bb(3)
    p2 = at(rs)%bb(1)
  else
    p1 = ni(rs)
    p2 = ci(rs-1)
  end if
  axl(1) = x(p1) - x(p2)
  axl(2) = y(p1) - y(p2)
  axl(3) = z(p1) - z(p2)
  mmm = sqrt(axl(1)*axl(1) + axl(2)*axl(2) + axl(3)*axl(3))
  do i=1,3
    axl(i) = axl(i)/mmm
  end do
  refv(1) = x(p1)
  refv(2) = y(p1)
  refv(3) = z(p1)
!
! now construct the quaternion that corresponds to the inverse omega move
  ll = wline(rs)
  dw = (cur_omega - ztorpr(ll))/RADIAN
  q(1) = cos(-dw/2.0)
  do i=2,4
    q(i) = axl(i-1)*sin(-dw/2.0)
  end do
! construct the multiplicative inverse
  mmm = q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4)
  if (mmm.le.0.0) then
    write(ilog,*) 'Fatal error. Got a 0-quat. in alignC_omega(...).'
    call fexit()
  end if
  qinv(1) = q(1)/mmm
  do i=2,4
    qinv(i) = -q(i)/mmm
  end do
!
! now apply the transformation
! note that this automatically takes care of the global vs. local frame issue
  do k=atmol(imol,1),atmol(imol,2)
    call quat_conjugate3(q,qinv,k,refv)
  end do
!
end
!
!--------------------------------------------------------------------------------
!
subroutine alignC_pivot(rs)
!
  use iounit
  use atoms
  use molecule
  use fyoc
  use math
  use polypep
  use zmatrix
  use sequen
!
  implicit none
!
  integer imol,i,k,rs,ff,yy
  RTYPE qf(4),qy(4),q(4),axlf(3),axly(3),refv(3),mmm,df,dy,old_psi
  RTYPE qinv(4)
!
  imol = molofrs(rs)
!
! first get the proper rotation axes, and C-alpha as the appropriate reference point
!
  if (rs.ne.rsmol(imol,1)) then
    axlf(1) = x(cai(rs)) - x(ni(rs))
    axlf(2) = y(cai(rs)) - y(ni(rs))
    axlf(3) = z(cai(rs)) - z(ni(rs))
    mmm = sqrt(axlf(1)*axlf(1) + axlf(2)*axlf(2) + axlf(3)*axlf(3))
    do i=1,3
      axlf(i) = axlf(i)/mmm
    end do
    ff = fline(rs)
    df = (cur_phi - ztorpr(ff))/RADIAN
  end if
  axly(1) = x(ci(rs)) - x(cai(rs))
  axly(2) = y(ci(rs)) - y(cai(rs))
  axly(3) = z(ci(rs)) - z(cai(rs))
  mmm = sqrt(axly(1)*axly(1) + axly(2)*axly(2) + axly(3)*axly(3))
  do i=1,3
    axly(i) = axly(i)/mmm
  end do
  refv(1) = x(cai(rs))
  refv(2) = y(cai(rs))
  refv(3) = z(cai(rs))
!
! now construct the quaternion (prod.) that corresponds to the inverse pivot move
  yy = yline(rs)
  old_psi = ztorpr(yy)
  if (old_psi.lt.0.0d0) then
    old_psi = old_psi + 180.0d0
  else
    old_psi = old_psi - 180.0d0
  end if
  dy = (cur_psi - old_psi)/RADIAN
  if (rs.ne.rsmol(imol,1)) then
    qf(1) = cos(-df/2.0)
    do i=2,4
      qf(i) = axlf(i-1)*sin(-df/2.0)
    end do
    qy(1) = cos(-dy/2.0)
    do i=2,4
      qy(i) = axly(i-1)*sin(-dy/2.0)
    end do
    call quat_product(qf,qy,q)
  else
    q(1) = cos(-dy/2.0)
    do i=2,4
      q(i) = axly(i-1)*sin(-dy/2.0)
    end do
  end if
! construct the multiplicative inverse
  mmm = q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4)
  if (mmm.le.0.0) then
    write(ilog,*) 'Fatal error. Got a 0-quat. in alignC_pivot(...).'
    call fexit()
  end if
  qinv(1) = q(1)/mmm
  do i=2,4
    qinv(i) = -q(i)/mmm
  end do
!
! now apply the transformation
! note that this automatically takes care of the global vs. local frame issue
  do k=atmol(imol,1),atmol(imol,2)
    call quat_conjugate3(q,qinv,k,refv)
  end do
!
end
!
!-----------------------------------------------------------------------------
!
subroutine alignC_CR(rsi,rsf,dfy)
!
  use iounit
  use atoms
  use molecule
  use fyoc
  use math
  use polypep
  use sequen
  use movesets
!
  implicit none
!
  integer imol,i,k,rsi,rsf,rs
  RTYPE qf(4),qy(4),q(4),axlf(3),axly(3),refv(3),mmm
  RTYPE qinv(4),dfy(MAXCRDOF),dum
!
  imol = molofrs(rsi)
!
! a SJ-CR move is essentially a series of pivot moves, but we cannot construct the net quaternion
! by simply concatenating all the individual rotations, because the rotation axes share no common
! point
! this is tedious of course, as it means that we must rotate all coordinates several times
! the correct solution to this problem would of course be to make the CR algorithm exact so that
! no C-tail alignment is necessary - ever
!
  do rs=rsi,rsf-1
    refv(1) = x(cai(rs))
    refv(2) = y(cai(rs))
    refv(3) = z(cai(rs))
    axlf(1) = x(cai(rs)) - x(ni(rs))
    axlf(2) = y(cai(rs)) - y(ni(rs))
    axlf(3) = z(cai(rs)) - z(ni(rs))
    mmm =sqrt(axlf(1)*axlf(1) + axlf(2)*axlf(2) + axlf(3)*axlf(3))
    do i=1,3
      axlf(i) = axlf(i)/mmm
    end do
    axly(1) = x(ci(rs)) - x(cai(rs))
    axly(2) = y(ci(rs)) - y(cai(rs))
    axly(3) = z(ci(rs)) - z(cai(rs))
    mmm = sqrt(axly(1)*axly(1) + axly(2)*axly(2) + axly(3)*axly(3))
    do i=1,3
      axly(i) = axly(i)/mmm
    end do
    qf(1) = cos(-dfy(2*(rs-rsi)+1)/2.0)
    dum = sin(-dfy(2*(rs-rsi)+1)/2.0)
    do i=2,4
      qf(i) = axlf(i-1)*dum
    end do
    qy(1) = cos(-dfy(2*(rs-rsi)+2)/(2.0))
    dum = sin(-dfy(2*(rs-rsi)+2)/2.0)
    do i=2,4
      qy(i) = axly(i-1)*dum
    end do
    call quat_product(qf,qy,q)
    mmm = q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4) 
    if (mmm.le.0.0) then
      write(ilog,*) 'Fatal error. Got a 0-quat. in alignC_CR(...).'
      call fexit()
    end if
    qinv(1) = q(1)/mmm
    do i=2,4
      qinv(i) = -q(i)/mmm
    end do
    do k=atmol(imol,1),atmol(imol,2)
      call quat_conjugate3(q,qinv,k,refv)
    end do
  end do
!
end
!
!-----------------------------------------------------------------------------
!
subroutine alignC_nuc(rs)
!
  use iounit
  use atoms
  use molecule
  use fyoc
  use math
  use polypep
  use sequen
  use movesets
  use zmatrix
  use aminos
!
  implicit none
!
  integer imol,i,k,rs,otze,ff,yy
  RTYPE qf(4),qy(4),q(4),axl1(3),axl2(3),refv(3),mmm
  RTYPE qinv(4),dy,df
!
  imol = molofrs(rs)
!
! second two for 5'-terminal nucleosides
  if (seqflag(rs).eq.24) then
!   these are C5* and C4* here
    axl1(1) = x(nuci(rs,3)) - x(nuci(rs,2))
    axl1(2) = y(nuci(rs,3)) - y(nuci(rs,2))
    axl1(3) = z(nuci(rs,3)) - z(nuci(rs,2))
    mmm = sqrt(axl1(1)*axl1(1) + axl1(2)*axl1(2) + axl1(3)*axl1(3))
    do i=1,3
      axl1(i) = axl1(i)/mmm
    end do
    refv(1) = x(nuci(rs,2))
    refv(2) = y(nuci(rs,2))
    refv(3) = z(nuci(rs,2))
!   now construct the quaternion (prod.) that corresponds to the inverse pivot move
    ff = nucsline(2,rs)
    df = (cur_nucs(2) - ztorpr(ff))/RADIAN
!   now construct the quaternion (prod.) that corresponds to the inverse pivot move
    qf(1) = cos(-df/2.0)
    do i=2,4
      qf(i) = axl1(i-1)*sin(-df/2.0)
    end do
!   construct the multiplicative inverse
    mmm = qf(1)*qf(1) + qf(2)*qf(2) + qf(3)*qf(3) + qf(4)*qf(4)
    if (mmm.le.0.0) then
      write(ilog,*) 'Fatal error. Got a 0-quat. in alignC_nuc(...).'
      call fexit()
    end if
    qinv(1) = qf(1)/mmm
    do i=2,4
      qinv(i) = -qf(i)/mmm
    end do
!
!   now apply the transformation
!   note that this automatically takes care of the global vs. local frame issue
    do k=atmol(imol,1),atmol(imol,2)
      call quat_conjugate3(qf,qinv,k,refv)
    end do
!
    if (rs.eq.rsmol(imol,2)) then
      write(ilog,*) '3"-Nucleosides not yet supported by C-terminal &
 &alignment for nucleic acid moves.'
      call fexit()
    else
      otze = nuci(rs+1,1)
    end if
    axl2(1) = x(otze) - x(nuci(rs,4))
    axl2(2) = y(otze) - y(nuci(rs,4))
    axl2(3) = z(otze) - z(nuci(rs,4))
    mmm = sqrt(axl2(1)*axl2(1) + axl2(2)*axl2(2) + axl2(3)*axl2(3))
    do i=1,3
      axl2(i) = axl2(i)/mmm
    end do
    refv(1) = x(nuci(rs,4))
    refv(2) = y(nuci(rs,4))
    refv(3) = z(nuci(rs,4))
    yy = nucsline(3,rs)
    dy = (cur_nucs(3) - ztorpr(yy))/RADIAN
!   now construct the quaternion (prod.) that corresponds to the inverse pivot move
    qy(1) = cos(-dy/2.0)
    do i=2,4
      qy(i) = axl2(i-1)*sin(-dy/2.0)
    end do
!   construct the multiplicative inverse
    mmm = qy(1)*qy(1) + qy(2)*qy(2) + qy(3)*qy(3) + qy(4)*qy(4)
    if (mmm.le.0.0) then
      write(ilog,*) 'Fatal error. Got a 0-quat. in alignC_nuc(...).'
      call fexit()
    end if
    qinv(1) = qy(1)/mmm
    do i=2,4
      qinv(i) = -qy(i)/mmm
    end do
!
!   now apply the transformation
!   note that this automatically takes care of the global vs. local frame issue
    do k=atmol(imol,1),atmol(imol,2)
      call quat_conjugate3(qy,qinv,k,refv)
    end do
    return
  end if
!
! for other (full) nucleic acid residues this is even more tedious of course,
! as it means that we must rotate all coordinates several times, since there
! are 5+ backbone angles with no common point on their rotation axes
!
! first two
  if (rs.ne.rsmol(imol,1)) then
    axl1(1) = x(nuci(rs,2)) - x(nuci(rs,1))
    axl1(2) = y(nuci(rs,2)) - y(nuci(rs,1))
    axl1(3) = z(nuci(rs,2)) - z(nuci(rs,1))
    mmm = sqrt(axl1(1)*axl1(1) + axl1(2)*axl1(2) + axl1(3)*axl1(3))
    do i=1,3
      axl1(i) = axl1(i)/mmm
    end do
    ff = nucsline(1,rs)
    df = (cur_nucs(1) - ztorpr(ff))/RADIAN
  end if
  axl2(1) = x(nuci(rs,3)) - x(nuci(rs,2))
  axl2(2) = y(nuci(rs,3)) - y(nuci(rs,2))
  axl2(3) = z(nuci(rs,3)) - z(nuci(rs,2))
  mmm = sqrt(axl2(1)*axl2(1) + axl2(2)*axl2(2) + axl2(3)*axl2(3))
  do i=1,3
    axl2(i) = axl2(i)/mmm
  end do
  refv(1) = x(nuci(rs,2))
  refv(2) = y(nuci(rs,2))
  refv(3) = z(nuci(rs,2))
!
! now construct the quaternion (prod.) that corresponds to the inverse pivot move
  yy = nucsline(2,rs)
  dy = (cur_nucs(2) - ztorpr(yy))/RADIAN
  if (rs.ne.rsmol(imol,1)) then
    qf(1) = cos(-df/2.0)
    do i=2,4
      qf(i) = axl1(i-1)*sin(-df/2.0)
    end do
    qy(1) = cos(-dy/2.0)
    do i=2,4
      qy(i) = axl2(i-1)*sin(-dy/2.0)
    end do
    call quat_product(qf,qy,q)
  else
    q(1) = cos(-dy/2.0)
    do i=2,4
      q(i) = axl2(i-1)*sin(-dy/2.0)
    end do
  end if
! construct the multiplicative inverse
  mmm = q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4)
  if (mmm.le.0.0) then
    write(ilog,*) 'Fatal error. Got a 0-quat. in alignC_nuc(...).'
    call fexit()
  end if
  qinv(1) = q(1)/mmm
  do i=2,4
    qinv(i) = -q(i)/mmm
  end do
!
! now apply the transformation
! note that this automatically takes care of the global vs. local frame issue
  do k=atmol(imol,1),atmol(imol,2)
    call quat_conjugate3(q,qinv,k,refv)
  end do
!
! three and four
  axl1(1) = x(nuci(rs,4)) - x(nuci(rs,3))
  axl1(2) = y(nuci(rs,4)) - y(nuci(rs,3))
  axl1(3) = z(nuci(rs,4)) - z(nuci(rs,3))
  mmm = sqrt(axl1(1)*axl1(1) + axl1(2)*axl1(2) + axl1(3)*axl1(3))
  do i=1,3
    axl1(i) = axl1(i)/mmm
  end do
  ff = nucsline(3,rs)
  df = (cur_nucs(3) - ztorpr(ff))/RADIAN
  axl2(1) = x(nuci(rs,5)) - x(nuci(rs,4))
  axl2(2) = y(nuci(rs,5)) - y(nuci(rs,4))
  axl2(3) = z(nuci(rs,5)) - z(nuci(rs,4))
  mmm = sqrt(axl2(1)*axl2(1) + axl2(2)*axl2(2) + axl2(3)*axl2(3))
  do i=1,3
    axl2(i) = axl2(i)/mmm
  end do
  refv(1) = x(nuci(rs,4))
  refv(2) = y(nuci(rs,4))
  refv(3) = z(nuci(rs,4))
!
! now construct the quaternion (prod.) that corresponds to the inverse pivot move
  yy = nucsline(4,rs)
  dy = (cur_nucs(4) - ztorpr(yy))/RADIAN
  qf(1) = cos(-df/2.0)
  do i=2,4
    qf(i) = axl1(i-1)*sin(-df/2.0)
  end do
  qy(1) = cos(-dy/2.0)
  do i=2,4
    qy(i) = axl2(i-1)*sin(-dy/2.0)
  end do
  call quat_product(qf,qy,q)
! construct the multiplicative inverse
  mmm = q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4)
  if (mmm.le.0.0) then
    write(ilog,*) 'Fatal error. Got a 0-quat. in alignC_nuc(...).'
    call fexit()
  end if
  qinv(1) = q(1)/mmm
  do i=2,4
    qinv(i) = -q(i)/mmm
  end do
!
! now apply the transformation
! note that this automatically takes care of the global vs. local frame issue
  do k=atmol(imol,1),atmol(imol,2)
    call quat_conjugate3(q,qinv,k,refv)
  end do
!
! five (only if not C-terminal, otherwise let hydrogen spin regardless
  if (rs.eq.rsmol(imol,2)) then
    return
  else
    otze = nuci(rs+1,1)
  end if
  axl2(1) = x(otze) - x(nuci(rs,6))
  axl2(2) = y(otze) - y(nuci(rs,6))
  axl2(3) = z(otze) - z(nuci(rs,6))
  mmm = sqrt(axl2(1)*axl2(1) + axl2(2)*axl2(2) + axl2(3)*axl2(3))
  do i=1,3
    axl2(i) = axl2(i)/mmm
  end do
  refv(1) = x(nuci(rs,6))
  refv(2) = y(nuci(rs,6))
  refv(3) = z(nuci(rs,6))
!
  yy = nucsline(5,rs)
  dy = (cur_nucs(5) - ztorpr(yy))/RADIAN
  q(1) = cos(-dy/2.0)
  do i=2,4
    q(i) = axl2(i-1)*sin(-dy/2.0)
  end do
! construct the multiplicative inverse
  mmm = q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4)
  if (mmm.le.0.0) then
    write(ilog,*) 'Fatal error. Got a 0-quat. in alignC_nuc(...).'
    call fexit()
  end if
  qinv(1) = q(1)/mmm
  do i=2,4
    qinv(i) = -q(i)/mmm
  end do
!
! now apply the transformation
! note that this automatically takes care of the global vs. local frame issue
  do k=atmol(imol,1),atmol(imol,2)
    call quat_conjugate3(q,qinv,k,refv)
  end do
!
end
!
!-----------------------------------------------------------------------------
!
subroutine alignC_sugar(rs)
!
  use iounit
  use atoms
  use molecule
  use fyoc
  use math
  use sequen
  use polypep
  use zmatrix
!
  implicit none
!
  integer imol,i,k,rs,ll,p1,p2
  RTYPE q(4),axl1(3),d34,refv(3),mmm,qinv(4)
!
  imol = molofrs(rs)
!
! first get the proper rotation axes, and an appropriate reference point
!
  if (seqflag(rs).eq.24) then
    p1 = nuci(rs,4)
    p2 = nuci(rs,3)
  else
    p1 = nuci(rs,6)
    p2 = nuci(rs,5)
  end if
  axl1(1) = x(p1) - x(p2)
  axl1(2) = y(p1) - y(p2)
  axl1(3) = z(p1) - z(p2)
  mmm = sqrt(axl1(1)*axl1(1) + axl1(2)*axl1(2) + axl1(3)*axl1(3))
  do i=1,3
    axl1(i) = axl1(i)/mmm
  end do
  refv(1) = x(p1)
  refv(2) = y(p1)
  refv(3) = z(p1)
!
! now construct the quaternion that corresponds to the inverse sugar pucker move
  ll = nucsline(6,rs)
  d34 = (cur_pucks(7) - ztorpr(ll))/RADIAN
  q(1) = cos(-d34/2.0)
  do i=2,4
    q(i) = axl1(i-1)*sin(-d34/2.0)
  end do
!
! construct the multiplicative inverse
  mmm = q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4)
  if (mmm.le.0.0) then
    write(ilog,*) 'Fatal error. Got a 0-quat. in alignC_sugar(...).'
    call fexit()
  end if
  qinv(1) = q(1)/mmm
  do i=2,4
    qinv(i) = -q(i)/mmm
  end do
!
! now apply the transformation
! note that this automatically takes care of the global vs. local frame issue
  do k=atmol(imol,1),atmol(imol,2)
    call quat_conjugate3(q,qinv,k,refv)
  end do
!
  end
!
!----------------------------------------------------------------------------------------
!
subroutine alignC_pucker(rs)
!
  use iounit
  use atoms
  use molecule
  use fyoc
  use math
  use sequen
  use polypep
  use zmatrix
!
  implicit none
!
  integer imol,i,k,rs,ll,p1,p2
  RTYPE bufv(3),getztor,efftor,q(4),axl1(3),d34,refv(3),mmm,qinv(4)
!
  imol = molofrs(rs)
!
! first get the proper rotation axes, and an appropriate reference point
!
  p1 = cai(rs)
  p2 = ni(rs)
  axl1(1) = x(p1) - x(p2)
  axl1(2) = y(p1) - y(p2)
  axl1(3) = z(p1) - z(p2)
  mmm = sqrt(axl1(1)*axl1(1) + axl1(2)*axl1(2) + axl1(3)*axl1(3))
  do i=1,3
    axl1(i) = axl1(i)/mmm
  end do
  refv(1) = x(p1)
  refv(2) = y(p1)
  refv(3) = z(p1)
!
! now construct the quaternion that corresponds to the inverse pucker move
  if (rs.eq.rsmol(imol,1)) then
    ll = fline(rs)
    bufv(1) = x(at(rs)%bb(5))
    bufv(2) = y(at(rs)%bb(5))
    bufv(3) = z(at(rs)%bb(5))
    x(at(rs)%bb(5)) = xref(at(rs)%bb(5))
    y(at(rs)%bb(5)) = yref(at(rs)%bb(5))
    z(at(rs)%bb(5)) = zref(at(rs)%bb(5))
    efftor = getztor(ci(rs),cai(rs),ni(rs),at(rs)%bb(5))
    x(at(rs)%bb(5)) = bufv(1)
    y(at(rs)%bb(5)) = bufv(2)
    z(at(rs)%bb(5)) = bufv(3)
  else
    ll = ci(rs)
    efftor = cur_pucks(1)
  end if
  d34 = (efftor - ztorpr(ll))/RADIAN
  q(1) = cos(-d34/2.0)
  do i=2,4
    q(i) = axl1(i-1)*sin(-d34/2.0)
  end do
!
! construct the multiplicative inverse
  mmm = q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4)
  if (mmm.le.0.0) then
    write(ilog,*) 'Fatal error. Got a 0-quat. in alignC_pucker(...).'
    call fexit()
  end if
  qinv(1) = q(1)/mmm
  do i=2,4
    qinv(i) = -q(i)/mmm
  end do
!
! now apply the transformation
! note that this automatically takes care of the global vs. local frame issue
  do k=atmol(imol,1),atmol(imol,2)
    call quat_conjugate3(q,qinv,k,refv)
  end do
!
  end
!
!----------------------------------------------------------------------------------------
!
subroutine transfer_rbc(mode,rbcvec)
!
  use molecule
  use atoms
!
  implicit none
!
  integer imol,mode
  RTYPE rbcvec(9*nmol)
!
  do imol=1,nmol
    if (mode.eq.1) then
      rbcvec(9*imol-8) = x(atmol(imol,1))
      rbcvec(9*imol-7) = y(atmol(imol,1))
      rbcvec(9*imol-6) = z(atmol(imol,1))
    else
      x(atmol(imol,1)) = rbcvec(9*imol-8)
      y(atmol(imol,1)) = rbcvec(9*imol-7)
      z(atmol(imol,1)) = rbcvec(9*imol-6)
    end if
    if ((atmol(imol,2)-atmol(imol,1)).eq.0) cycle
    if (mode.eq.1) then
      rbcvec(9*imol-5) = x(atmol(imol,1)+1)
      rbcvec(9*imol-4) = y(atmol(imol,1)+1)
      rbcvec(9*imol-3) = z(atmol(imol,1)+1)
    else
      x(atmol(imol,1)+1) = rbcvec(9*imol-5)
      y(atmol(imol,1)+1) = rbcvec(9*imol-4)
      z(atmol(imol,1)+1) = rbcvec(9*imol-3)
    end if
    if ((atmol(imol,2)-atmol(imol,1)).eq.1) cycle
    if (mode.eq.1) then
      rbcvec(9*imol-2) = x(atmol(imol,1)+2)
      rbcvec(9*imol-1) = y(atmol(imol,1)+2)
      rbcvec(9*imol) =   z(atmol(imol,1)+2)
    else
      x(atmol(imol,1)+2) = rbcvec(9*imol-2)
      y(atmol(imol,1)+2) = rbcvec(9*imol-1)
      z(atmol(imol,1)+2) = rbcvec(9*imol)
    end if
  end do
!
end
!
!---------------------------------------------------------------------------------------
!
subroutine crosslink_follow(lk,imol,jmol,oldints)
!
  use sequen
  use molecule
  use iounit
  use atoms
  use polypep
  use system
!
  implicit none
!
  integer i,lk,imol,jmol,sz,azero,rs1,rs2,ptlst(5),shf
  RTYPE olds(9),oldints(9),news(9),qrot(4),tvec(3),roid1(3),roid2(3)
!
  azero = 0
  rs1 = crosslink(lk)%rsnrs(1)
  rs2 = crosslink(lk)%rsnrs(2)
!
  if (imol.eq.jmol) then
    write(ilog,*) 'Fatal. Encountered intramolecular crosslink in crosslink_follow(...).&
 & This is a bug.'
    call fexit()
  end if
!
  if (crosslink(lk)%itstype.eq.2) then
!   generate dummies for positions of SG, CB, and CA on rs2
    shf = 0
    if (ua_model.gt.0) shf = 1
    ptlst(1) = n+20-3
    call genxyz(ptlst(1),at(rs1)%sc(3-shf),oldints(3),at(rs1)%sc(2-shf),oldints(2),cai(rs1),oldints(1),azero)
    news(1) = x(ptlst(1))
    news(2) = y(ptlst(1))
    news(3) = z(ptlst(1))
    ptlst(2) = n+20-2
    call genxyz(ptlst(2),ptlst(1),oldints(6),at(rs1)%sc(3-shf),oldints(5),at(rs1)%sc(2-shf),oldints(4),azero)
    news(4) = x(ptlst(2))
    news(5) = y(ptlst(2))
    news(6) = z(ptlst(2))
    ptlst(3) = n+20-1
    call genxyz(ptlst(3),ptlst(2),oldints(9),ptlst(1),oldints(8),at(rs1)%sc(3-shf),oldints(7),azero)
    news(7) = x(ptlst(3))
    news(8) = y(ptlst(3))
    news(9) = z(ptlst(3))
    olds(1) = x(at(rs2)%sc(3-shf))
    olds(2) = y(at(rs2)%sc(3-shf))
    olds(3) = z(at(rs2)%sc(3-shf))
    olds(4) = x(at(rs2)%sc(2-shf))
    olds(5) = y(at(rs2)%sc(2-shf))
    olds(6) = z(at(rs2)%sc(2-shf))
    olds(7) = x(cai(rs2))
    olds(8) = y(cai(rs2))
    olds(9) = z(cai(rs2))
    sz = 3
!   obtain a translation vector and rotation quaternion to maintain correct linkage
    call align_3D(sz,olds,news,tvec,qrot,roid1,roid2)
!   rotate in place
    do i=atmol(molofrs(rs2),1),atmol(molofrs(rs2),2)
      call quat_conjugate4(qrot,i,roid1)
    end do
!   translate
    do i=atmol(molofrs(rs2),1),atmol(molofrs(rs2),2)
      x(i) = x(i) + tvec(1)
      y(i) = y(i) + tvec(2)
      z(i) = z(i) + tvec(3)
    end do
    olds(1) = x(at(rs2)%sc(3-shf))
    olds(2) = y(at(rs2)%sc(3-shf))
    olds(3) = z(at(rs2)%sc(3-shf))
    olds(4) = x(at(rs2)%sc(2-shf))
    olds(5) = y(at(rs2)%sc(2-shf))
    olds(6) = z(at(rs2)%sc(2-shf))
    olds(7) = x(cai(rs2))
    olds(8) = y(cai(rs2))
    olds(9) = z(cai(rs2))
  else
    write(ilog,*) 'Fatal. Encountered unsupported crosslink type in crosslink_follow(...).&
 & This is most likely an omission bug.'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------------------------
!

