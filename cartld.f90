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
! #######################################################
! ##                                                   ##
! ## subroutine ldmove -- Stochastic Dynamics moves    ##
! ##                                                   ##
! #######################################################
!     
! "ldmove" calls the LD integrator
!
!
subroutine cart_ldmove(istep,ndump,forcefirst)
!
  use iounit
  use math
  use atoms
  use molecule
  use forces
  use energies
  use movesets
  use system
  use units
  use mcsums
  use movesets
  use fyoc
  use torsn
  use cutoffs
  use shakeetal
!
  implicit none
!
  integer imol,j,istep,ndump,aone,ttc,rs,i
  RTYPE force3
  RTYPE ompl,ommi,intS,bold,rndnew,t1,t2,normal
  RTYPE exgdt2,exgdt
  RTYPE dsplxyz(3),boxovol
  logical afalse,forcefirst
!
  afalse = .false.
!
 57   format('Mol. ',i6,' Res. ',i8,' #',i2,' :',g14.7)
 58   format('Mol. ',i6,' Res. ',i8,' :',g14.7)
!
! set the constants for the Langevin
  exgdt = exp(-fric_ga*dyn_dt)
  exgdt2 = exp(-fric_ga*0.5*dyn_dt)
  ompl = (exgdt - 1.0 + fric_ga*dyn_dt)/&
 &    (fric_ga*dyn_dt*(1.0 - exgdt))
  ommi = 1.0 - ompl
  intS = 0.0
!
! missing initialization
!
  aone = 1
  nstep = nstep + 1
  mvcnt%nld = mvcnt%nld + 1
!
  if ((istep.eq.1).OR.(forcefirst.EQV..true.)) then
!    write(*,*) (1.0 - exgdt)/(fric_ga*exgdt2),exgdt2,ompl
    do imol=1,nmol
!     sanity check in first step against (partially) unsupported 2-atom molecules
      if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        write(ilog,*) 'Fatal. Two-atom molecules are not yet support&
 &ed in dynamics. Check back later.'
        call fexit()
      end if
!     make sure rigid body coordinates are up-to-date
      call update_rigidm(imol)
    end do
!   get the deterministic force
    call CPU_time(t1)
    esave = force3(esterms,esterms_tr,esterms_lr,forcefirst)
    if (no_shake.EQV..false.) call cart2cart_shake()
    ens%insU = esave
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
    ttc = 0
!   set the frictional parameters
    do i=1,n
      if (mass(i).le.0.0) cycle
      do j=1,3
        if (cart_frz(i,j).EQV..true.) cycle
        cart_ldp(i,j,1) = (u_dyn_kb*kelvin/mass(i))*(2.0*ompl*ompl*fric_ga*dyn_dt + ompl - ommi)
        cart_ldp(i,j,2) = (u_dyn_kb*kelvin/mass(i))*(2.0*ompl*ommi*fric_ga*dyn_dt - ompl + ommi)
        cart_ldp(i,j,3) = (u_dyn_kb*kelvin/mass(i))*(2.0*ommi*ommi*fric_ga*dyn_dt + ompl - ommi)
      end do
    end do
!
    call randomize_cart_velocities()
    cart_a(:,:) = cart_v(:,:) ! velocity back-up
!
!   now do the first step integration
!   now loop over all molecules
    do imol=1,nmol
!
      call makeref_formol(imol)
!
      do i=atmol(imol,1),atmol(imol,2)
        if (mass(i).le.0.0) cycle
!       gather the displacements in atomic coordinate space, update velocities    
        do j=1,3
          if (cart_frz(i,j).EQV..true.) then
            cart_v(i,j) = 0.0
            dsplxyz(j) = 0.0
            cycle
          end if
!         "increment" frictional/stochastic parameter to starting value
          cart_ldp(i,j,4) = sqrt(cart_ldp(i,j,1))
          cart_ldp(i,j,5) = normal()
          rndnew = cart_ldp(i,j,4)*cart_ldp(i,j,5)
!         increment velocity from time zero to 1/2dt
          cart_v(i,j) = exgdt2*(cart_v(i,j) + ompl*u_dyn_fconv*dyn_dt*cart_f(i,j)/mass(i) + rndnew)
!         increment position t to t + dt with staggered velocity at t + 1/2dt
          dsplxyz(j) = dyn_dt*cart_v(i,j)
        end do
        x(i) = x(i) + dsplxyz(1)
        y(i) = y(i) + dsplxyz(2)
        z(i) = z(i) + dsplxyz(3)
      end do
!
    end do
!
    if (no_shake.EQV..false.) call shake_wrap()
!
    do imol=1,nmol
!
!     from fresh xyz, we need to recompute internals (rigid-body and torsions)
      call update_rigidm(imol)
      call update_comv(imol)
      call genzmat(imol)
!     and shift molecule into central cell if necessary
      call update_image(imol)
!     update grid association if necessary
      if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
        do rs=rsmol(imol,1),rsmol(imol,2)
          call updateresgp(rs)
        end do
      end if
!
    end do
!
!   update pointer arrays such that analysis routines can work properly 
    call zmatfyc2()
!
  else
!
    cart_a(:,:) = cart_v(:,:) ! velocity back-up
!
!   get the deterministic force
    call CPU_time(t1)
    esave = force3(esterms,esterms_tr,esterms_lr,forcefirst)
    if (no_shake.EQV..false.) call cart2cart_shake()
    ens%insU = esave
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   set the frictional parameters
    do i=1,n
      if (mass(i).le.0.0) cycle
      do j=1,3
        if (cart_frz(i,j).EQV..true.) cycle
        cart_ldp(i,j,1) = (u_dyn_kb*kelvin/mass(i))*(2.0*ompl*ompl*fric_ga*dyn_dt + ompl - ommi)
        cart_ldp(i,j,2) = (u_dyn_kb*kelvin/mass(i))*(2.0*ompl*ommi*fric_ga*dyn_dt - ompl + ommi)
        cart_ldp(i,j,3) = (u_dyn_kb*kelvin/mass(i))*(2.0*ommi*ommi*fric_ga*dyn_dt + ompl - ommi)
      end do
    end do
!
!   now loop over all molecules
    do imol=1,nmol
      call makeref_formol(imol)
!
      do i=atmol(imol,1),atmol(imol,2)
        if (mass(i).le.0.0) cycle
!       gather the displacements in atomic coordinate space, update velocities    
        do j=1,3
          if (cart_frz(i,j).EQV..true.) then
            cart_v(i,j) = 0.0
            dsplxyz(j) = 0.0
            cycle
          end if
!         increment frictional/stochastic parameters to current values (trailing)
          bold = cart_ldp(i,j,2)/cart_ldp(i,j,4)
          cart_ldp(i,j,4) = sqrt(cart_ldp(i,j,1) + cart_ldp(i,j,3) - bold*bold)
          rndnew = bold*cart_ldp(i,j,5)
          cart_ldp(i,j,5) = normal()
          rndnew = rndnew + cart_ldp(i,j,4)*cart_ldp(i,j,5)
!         now increment velocity at time t - 1/2dt to t + 1/2dt
          cart_v(i,j) = exgdt2*(exgdt2*cart_v(i,j) + u_dyn_fconv*dyn_dt*cart_f(i,j)/mass(i) + rndnew)
!         and finally increment position to t + dt
          dsplxyz(j) = dyn_dt*cart_v(i,j)
        end do
        x(i) = x(i) + dsplxyz(1)
        y(i) = y(i) + dsplxyz(2)
        z(i) = z(i) + dsplxyz(3)
      end do
!
    end do
!
    if (no_shake.EQV..false.) call shake_wrap()
!
    do imol=1,nmol
!
!     from fresh xyz, we need to recompute internals (rigid-body and torsions)
      call update_rigidm(imol)
      call update_comv(imol)
      call genzmat(imol)
!     and shift molecule into central cell if necessary
      call update_image(imol)
!     update grid association if necessary
      if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
        do rs=rsmol(imol,1),rsmol(imol,2)
          call updateresgp(rs)
        end do
      end if
!
    end do
!
!   update pointer arrays such that analysis routines can work properly 
    call zmatfyc2()
!
  end if
!
  cart_a(1:n,:) = (cart_v(1:n,:)-cart_a(1:n,:))/dyn_dt ! finite diff. accelerations
  ens%insR(1) = 0.0
  do i=1,n
    if (mass(i).gt.0.0) then
      ens%insR(1) = ens%insR(1) + sum((mass(i)*cart_a(i,:)/u_dyn_fconv - cart_f(i,:))**2)/mass(i)
    end if
  end do
  call prt_curens(istep,boxovol,forcefirst)
!
end
!
!---------------------------------------------------------------------------
!
! this routine is only used to ensure that a run which was restarted from
! an MC restart file. this is necessary because the non-first step LD
! algorithm requires ldp(i,4) and ldp(i,5)
! of course, it required inertial masses to be set and up-to-date
!
subroutine init_cart_ldps()
!
  use iounit
  use forces
  use atoms
  use system
  use units
!
  implicit none
!
  integer i,j
  RTYPE ompl,ommi,normal
  RTYPE exgdt2,exgdt
!
! set the constants
  exgdt = exp(-fric_ga*dyn_dt)
  exgdt2 = exp(-fric_ga*0.5*dyn_dt)
  ompl = (exgdt - 1.0 + fric_ga*dyn_dt)/&
 &    (fric_ga*dyn_dt*(1.0 - exgdt))
  ommi = 1.0 - ompl
  exgdt = exp(-fric_ga*dyn_dt)
  do i=1,n
    if (mass(i).le.0.0) cycle
    do j=1,3
      cart_ldp(i,j,1) = (u_dyn_kb*kelvin/mass(i))*&
 &                  (2.0*ompl*ompl*fric_ga*dyn_dt + ompl - ommi)
      cart_ldp(i,j,2) = (u_dyn_kb*kelvin/mass(i))*&
 &                  (2.0*ompl*ommi*fric_ga*dyn_dt - ompl + ommi)
      cart_ldp(i,j,3) = (u_dyn_kb*kelvin/mass(i))*&
 &                  (2.0*ommi*ommi*fric_ga*dyn_dt + ompl - ommi)
      cart_ldp(i,j,4) = sqrt(cart_ldp(i,j,1))
      cart_ldp(i,j,5) = normal()
    end do
  end do
!
end
!
!-----------------------------------------------------------------------
!
! a routine that does nothing but assemble CLD parameters into provided arrays
!
subroutine manage_clds(vars,mode)
!
  use iounit
  use forces
  use atoms
!
  implicit none
!
  integer mode
  RTYPE vars(n,9)
!
  if (mode.eq.1) then
!   store in vars
    vars(1:n,1) = cart_v(1:n,1)
    vars(1:n,2) = cart_v(1:n,2)
    vars(1:n,3) = cart_v(1:n,3)
    vars(1:n,4) = cart_ldp(1:n,1,4)
    vars(1:n,5) = cart_ldp(1:n,1,5)
    vars(1:n,6) = cart_ldp(1:n,2,4)
    vars(1:n,7) = cart_ldp(1:n,2,5)
    vars(1:n,8) = cart_ldp(1:n,3,4)
    vars(1:n,9) = cart_ldp(1:n,3,5)
  else if (mode.eq.2) then
    cart_v(1:n,1) = vars(1:n,1)
    cart_v(1:n,2) = vars(1:n,2)
    cart_v(1:n,3) = vars(1:n,3)
    cart_ldp(1:n,1,4) = vars(1:n,4)
    cart_ldp(1:n,1,5) = vars(1:n,5)
    cart_ldp(1:n,2,4) = vars(1:n,6)
    cart_ldp(1:n,2,5) = vars(1:n,7)
    cart_ldp(1:n,3,4) = vars(1:n,8)
    cart_ldp(1:n,3,5) = vars(1:n,9)
  else
    write(ilog,*) 'Fatal. Called manage_clds(...) with unknown mode (identifier is ',&
  &mode,'). This is a bug.'
    call fexit()
  end if
!
  end
!
!-----------------------------------------------------------------------
!
