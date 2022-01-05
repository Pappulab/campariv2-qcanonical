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
! CONTRIBUTIONS: Jose Pulido                                               !
!                                                                          !
!--------------------------------------------------------------------------!
!     
!
#include "macros.i"
!
!-----------------------------------------------------------------------
!
! A collection of subroutines related to computing quantities
! relevant for characterizing the ensemble as a whole
!
!-----------------------------------------------------------------------
!
subroutine ensv_helper()
!
  use system
  use torsn
  use forces
  use atoms
  use units
  use molecule
  use math
!
  implicit none
!
  integer imol,boxdof,i,j,ttc
  RTYPE ttt,temp,boxovol,ttg(tstat%n_tgrps)
!
  ttt = 0.0
  ttg(:) = 0.0
!
  if (fycxyz.eq.1) then
    if (dyn_integrator_ops(1).gt.0) then ! we have estimated inertia at time t1.0
      do imol=1,nmol
        do j=1,3
          temp = dc_di(imol)%olddat(j,2)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        end do
        if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
          ttc = 3
        else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
          ttc = 5
!         WARNING: support missing
        else
          ttc = 6
          do j=4,6
            temp = (1.0/(RADIAN*RADIAN))*&
   &             dc_di(imol)%olddat(j,2)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
            ttt = ttt + temp
            ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
          end do    
        end if
        do j=ttc+1,ttc+ntormol(moltypid(imol))
          temp = (1.0/(RADIAN*RADIAN))*&
   &            dc_di(imol)%olddat(j,2)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        end do
      end do
    else
      do imol=1,nmol
        do j=1,3
          temp = dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        end do
        if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
          ttc = 3
        else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
          ttc = 5
!         WARNING: support missing
        else
          ttc = 6
          do j=4,6
            temp = (1.0/(RADIAN*RADIAN))*&
   &             dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
            ttt = ttt + temp
            ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
          end do    
        end if
        do j=ttc+1,ttc+ntormol(moltypid(imol))
          temp = (1.0/(RADIAN*RADIAN))*&
   &            dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        end do
      end do
    end if
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      if (pstat%flag.eq.1) then
        boxdof = 1
        temp = pstat%params(1)*bnd_v*bnd_v
        ttt = ttt + temp
        ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + temp
      else if (pstat%flag.eq.3) then
!       do nothing (extended ensemble particles not considered)
        boxdof = 0
      end if
    else
      boxdof = 0
    end if
    tstat%grpT(:) = ttg(:) / (tstat%grpdof(:)*u_dyn_kb)
    ens%insR(8) = ens%insK
    ens%insK = (1.0/u_dyn_fconv)*0.5*ttt
    ens%insT = ttt / ((totrbd+ndyntorsn+boxdof-n_constraints)*u_dyn_kb)
    ens%insR(9) = ens%insK2
    ttt = sum(mass(1:n)*(cart_v(1:n,1)*cart_v(1:n,1) + cart_v(1:n,2)*cart_v(1:n,2) + cart_v(1:n,3)*cart_v(1:n,3)))
    ens%insK2 = (1.0/u_dyn_fconv)*0.5*ttt
!
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then ! nonfunctional
      call virial(ens%insVirT)
      call ekinten(ens%insKinT)
      ens%insVirT(:,:) = 0.0
      ens%insP = 0.0
      ens%insP = (1./3.)*(u_dyn_virconv/boxovol)*((1.0/u_dyn_fconv)*&
   &            (ens%insKinT(1,1)+ens%insKinT(2,2)+ens%insKinT(3,3)) +&
   &            (ens%insVirT(1,1)+ens%insVirT(2,2)+ens%insVirT(3,3)) )
    end if
!
  else if (fycxyz.eq.2) then
!
    do imol=1,nmol
      do i=atmol(imol,1),atmol(imol,2)
        do j=1,3
          temp = mass(i)*cart_v(i,j)*cart_v(i,j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        end do
      end do
    end do
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      if (pstat%flag.eq.1) then
        boxdof = 1
        temp = pstat%params(1)*bnd_v*bnd_v
        ttt = ttt + temp
        ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + temp
      else if (pstat%flag.eq.2) then
        boxdof = 0
!       do nothing (Markov chain not considered)
      else if (pstat%flag.eq.3) then
        boxdof = 0
!       do nothing (extended ensemble particles not considered)
      end if
    else
      boxdof = 0
    end if
    tstat%grpT(:) = ttg(:) / (tstat%grpdof(:)*u_dyn_kb)
    ens%insR(8) = ens%insK
    ens%insK = (1.0/u_dyn_fconv)*0.5*ttt
    ens%insT = ttt / (dble(3*n+boxdof-n_constraints)*u_dyn_kb)
!
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then ! nonfunctional
      call virial(ens%insVirT)
      call ekinten(ens%insKinT)
      ens%insVirT(:,:) = 0.0
      ens%insP = 0.0
      ens%insP = (1./3.)*(u_dyn_virconv/boxovol)*((1.0/u_dyn_fconv)*&
   &        dble(3.0*n+boxdof)/dble(3*n+boxdof-n_constraints)*&
   &            (ens%insKinT(1,1)+ens%insKinT(2,2)+ens%insKinT(3,3)) +&
   &            (ens%insVirT(1,1)+ens%insVirT(2,2)+ens%insVirT(3,3)) )
    end if
  end if
end
!
!----------------------------------------------------------------------------------
!
subroutine get_ensv(boxovol)
!
  use iounit
  use system
  use molecule
  use forces
  use units
  use mcsums
  use math
  use torsn
  use energies
  use atoms
!
  implicit none
!
  integer imol,boxdof,i,j,ttc
  RTYPE ttt,temp,boxovol,ttg(tstat%n_tgrps)
!
  ttt = 0.0
  ttg(:) = 0.0
  if (dyn_integrator_ops(1).gt.0) then ! we have estimated inertia at time t1.0
    do imol=1,nmol
      do j=1,3
        temp = dc_di(imol)%olddat(j,2)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
        ttt = ttt + temp
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        if (nstep.gt.nequil) then
          dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp/u_dyn_kb
        end if
      end do
      if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
        ttc = 3
      else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        ttc = 5
!       WARNING: support missing
      else
        ttc = 6
        do j=4,6
          temp = (1.0/(RADIAN*RADIAN))*&
 &             dc_di(imol)%olddat(j,2)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
          if (nstep.gt.nequil) then
            dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp/u_dyn_kb
          end if
        end do    
      end if
      do j=ttc+1,ttc+ntormol(moltypid(imol))
        temp = (1.0/(RADIAN*RADIAN))*&
 &            dc_di(imol)%olddat(j,2)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
        ttt = ttt + temp
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        if (nstep.gt.nequil) then
          dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp/u_dyn_kb
        end if
      end do
    end do
  else
    do imol=1,nmol
      do j=1,3
        temp = dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
        ttt = ttt + temp
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        if (nstep.gt.nequil) then
          dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp/u_dyn_kb
        end if
      end do
      if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
        ttc = 3
      else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        ttc = 5
!       WARNING: support missing
      else
        ttc = 6
        do j=4,6
          temp = (1.0/(RADIAN*RADIAN))*&
 &             dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
          if (nstep.gt.nequil) then
            dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp/u_dyn_kb
          end if
        end do    
      end if
      do j=ttc+1,ttc+ntormol(moltypid(imol))
        temp = (1.0/(RADIAN*RADIAN))*&
 &            dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
        ttt = ttt + temp
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        if (nstep.gt.nequil) then
          dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp/u_dyn_kb
        end if
      end do
    end do
  end if
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      boxdof = 1
      temp = pstat%params(1)*bnd_v*bnd_v
      ttt = ttt + temp
      ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + temp
      if (nstep.gt.nequil) then
        bnd_avgT = bnd_avgT + temp/u_dyn_kb
      end if
    else if (pstat%flag.eq.3) then
!     do nothing (extended ensemble particles not considered)
      boxdof = 0
    else
      write(ilog,*) 'Fatal. Unsupported manostat in get_ensv(...). T&
 &his is most likely an omission bug.'
      call fexit()
    end if
  else
    boxdof = 0
  end if
  tstat%grpT(:) = ttg(:) / (tstat%grpdof(:)*u_dyn_kb)
  ens%insR(8) = ens%insK
  ens%insK = (1.0/u_dyn_fconv)*0.5*ttt
  ens%insT = ttt / ((totrbd+ndyntorsn+boxdof-n_constraints)*u_dyn_kb)
  ens%insR(9) = ens%insK2
  ttt = sum(mass(1:n)*(cart_v(1:n,1)*cart_v(1:n,1) + cart_v(1:n,2)*cart_v(1:n,2) + cart_v(1:n,3)*cart_v(1:n,3)))
  ens%insK2 = (1.0/u_dyn_fconv)*0.5*ttt
!
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then ! nonfunctional
    call virial(ens%insVirT)
    call ekinten(ens%insKinT)
    ens%insVirT(:,:) = 0.0
    ens%insP = 0.0
    ens%insP = (1./3.)*(u_dyn_virconv/boxovol)*((1.0/u_dyn_fconv)*&
 &            (ens%insKinT(1,1)+ens%insKinT(2,2)+ens%insKinT(3,3)) +&
 &            (ens%insVirT(1,1)+ens%insVirT(2,2)+ens%insVirT(3,3)) )
  end if
!  
  if (nstep.gt.nequil) then
    if (nstep.eq.(nequil+1)) then
      if (ens%flag.eq.2) then
        ens%insR(2) = ens%insK + ens%insU
        ens%insR(3) = ens%insK2 + ens%insU
      end if
    end if
    ens%avgcnt = ens%avgcnt + 1
    ens%avgK = ens%avgK + ens%insK
    ens%avgK2 = ens%avgK2 + ens%insK2
    ens%avgU = ens%avgU + ens%insU
    ens%avgT = ens%avgT + ens%insT
    ens%avgR(4) = ens%avgR(4) + abs(ens%insK2-ens%insK)
    ens%avgR(1) = ens%avgR(1) + ens%insR(1)
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
!      ens%avgV = ens%avgV + (nmol+1)*ens%insT*gasconst/
! &                              (boxovol*ens%insP/u_dyn_virconv)
      ens%avgV = ens%avgV + boxovol/1.0e6
      do i=1,3
        do j=1,3
          ens%avgVirT(i,j) = ens%avgVirT(i,j) +(u_dyn_virconv/boxovol)*ens%insVirT(i,j)
          ens%avgPT(i,j) = ens%avgPT(i,j) + (u_dyn_virconv/boxovol)*((1.0/u_dyn_fconv)*ens%insKinT(i,j) + ens%insVirT(i,j))
        end do
      end do
      ens%avgP = ens%avgP + ens%insP
    else
      ens%avgR(2) = ens%avgR(2) + (ens%insU + ens%insK)**2
      ens%avgR(3) = ens%avgR(3) + (ens%insU + ens%insK2)**2
      ens%avgR(5) = ens%avgR(5) + ens%insU**2
      ens%avgR(6) = ens%avgR(6) + ens%insK**2
      ens%avgR(7) = ens%avgR(7) + ens%insK2**2
    end if
  end if
!
end
!
!-----------------------------------------------------------------------
!
! same as get_ensv for straight Cartesian dynamics
!
subroutine get_cart_ensv(boxovol)
!
  use iounit
  use system
  use molecule
  use forces
  use units
  use mcsums
  use atoms
  use energies
!
  implicit none
!
  integer imol,boxdof,i,j
  RTYPE ttt,temp,boxovol,ttg(tstat%n_tgrps)
!
  ttg(:) = 0.0
  ttt = 0.0
  do imol=1,nmol
    do i=atmol(imol,1),atmol(imol,2)
      do j=1,3
        temp = mass(i)*cart_v(i,j)*cart_v(i,j)
        ttt = ttt + temp
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
!        if (nstep.gt.nequil) then
!          cart_avgT(i,j) = cart_avgT(i,j) + temp/u_dyn_kb
!        end if
      end do
    end do
  end do
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      boxdof = 1
      temp = pstat%params(1)*bnd_v*bnd_v
      ttt = ttt + temp
      ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + temp
      if (nstep.gt.nequil) then
        bnd_avgT = bnd_avgT + temp/u_dyn_kb
      end if
    else if (pstat%flag.eq.2) then
      boxdof = 0
!     do nothing (Markov chain not considered)
    else if (pstat%flag.eq.3) then
      boxdof = 0
!     do nothing (extended ensemble particles not considered)
    else
      write(ilog,*) 'Fatal. Unsupported manostat in get_cart_ensv(..&
 &.). This is most likely an omission bug.'
      call fexit()
    end if
  else
    boxdof = 0
  end if
  tstat%grpT(:) = ttg(:) / (tstat%grpdof(:)*u_dyn_kb)
  ens%insR(8) = ens%insK
  ens%insK = (1.0/u_dyn_fconv)*0.5*ttt
  ens%insT = ttt / (dble(3*n+boxdof-n_constraints)*u_dyn_kb)
!
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then ! nonfunctional
    call virial(ens%insVirT)
    call ekinten(ens%insKinT)
    ens%insVirT(:,:) = 0.0
    ens%insP = 0.0
    ens%insP = (1./3.)*(u_dyn_virconv/boxovol)*((1.0/u_dyn_fconv)*&
 &        dble(3.0*n+boxdof)/dble(3*n+boxdof-n_constraints)*&
 &            (ens%insKinT(1,1)+ens%insKinT(2,2)+ens%insKinT(3,3)) +&
 &            (ens%insVirT(1,1)+ens%insVirT(2,2)+ens%insVirT(3,3)) )
  end if
!  
  if (nstep.gt.nequil) then
    if (nstep.eq.(nequil+1)) then
      if (ens%flag.eq.2) then
        ens%insR(2) = ens%insK + ens%insU
      end if
    end if
    ens%avgcnt = ens%avgcnt + 1
    ens%avgK = ens%avgK + ens%insK
    ens%avgU = ens%avgU + ens%insU
    ens%avgT = ens%avgT + ens%insT
    ens%avgR(1) = ens%avgR(1) + ens%insR(1)
    ens%avgR(8) = ens%avgR(8) + ens%insR(6)
    ens%avgR(9) = ens%avgR(9) + ens%insR(7)
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
!      ens%avgV = ens%avgV + (nmol+1)*ens%insT*gasconst/
! &                              (boxovol*ens%insP/u_dyn_virconv)
      ens%avgV = ens%avgV + boxovol/1.0e6
      do i=1,3
        do j=1,3
       ens%avgVirT(i,j) = ens%avgVirT(i,j) +(u_dyn_virconv/boxovol)*ens%insVirT(i,j)
!        ens%avgPT(i,j) = ens%avgPT(i,j) + (u_dyn_virconv/boxovol)*(
! &                            (1.0/u_dyn_fconv)*ens%insKinT(i,j) +
! &                                              ens%insVirT(i,j) )
        end do
      end do
      ens%avgP = ens%avgP + ens%insP
    else
      ens%avgR(2) = ens%avgR(2) + (ens%insU + ens%insK)**2
      ens%avgR(5) = ens%avgR(5) + ens%insU**2
      ens%avgR(6) = ens%avgR(6) + ens%insK**2
    end if
  end if
!
end
!
!---------------------------------------------------------------------------
!
! the routine to initialize dof-velocities to a Boltzmann-distributed
! ensemble using target the target temperature and (current) masses
!
subroutine randomize_velocities()
!
  use system
  use forces
  use atoms
  use molecule
  use units
  use math
  use torsn
  use iounit
!
  implicit none
!
  integer imol,j,k,i,intdim,ttc,boxdof
  RTYPE vind,normal,ttt,fr1,fr2,temp
  RTYPE ttg(tstat%n_tgrps),corrfac(tstat%n_tgrps)
!
  do imol=1,nmol
!
!   first: center-of-mass translation
!
    do j=1,3
!     velocities in A/ps
      if (dc_di(imol)%frz(j).EQV..true.) then
        dc_di(imol)%v(j) = 0.0
        cycle
      end if
      vind = normal()
      vind = vind*sqrt(u_dyn_kb*kelvin/dc_di(imol)%im(j))
      dc_di(imol)%v(j) = vind
    end do
!
    if ((atmol(imol,2)-atmol(imol,1)).eq.0) cycle
!
!   second: rigid-body rotation

    if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
!
      ttc = 5
      write(ilog,*) 'Fatal. Diatomic molecules currently not support&
 &ed by dynamics implementation.'
      call fexit()
!     WARNING: support missing
!
    else
!
      ttc = 6
      do j=4,6
!       velocities in deg/ps
        if (dc_di(imol)%frz(j).EQV..true.) then
          dc_di(imol)%v(j) = 0.0
          cycle
        end if
        vind = normal()
        vind = vind*sqrt(u_dyn_kb*kelvin/dc_di(imol)%im(j))
        dc_di(imol)%v(j) = RADIAN*vind
      end do
!
    end if
!
!   third: torsional degrees of freedom
!
    do j=1,ntormol(moltypid(imol))
      if (dc_di(imol)%frz(ttc+j).EQV..true.) then
        dc_di(imol)%v(ttc+j) = 0.0
        cycle
      end if
      vind = normal()
      vind = vind*sqrt(u_dyn_kb*kelvin/dc_di(imol)%im(ttc+j))
      dc_di(imol)%v(ttc+j) = RADIAN*vind
    end do
!
  end do
!
! fourth: box degrees of freedom (if applicable)
!
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      vind = normal()
      vind = vind*sqrt(u_dyn_kb*kelvin/pstat%params(1))
      bnd_v = vind
    else if (pstat%flag.eq.3) then
      pstat%params(2) = 0.0
    else
      write(ilog,*) 'Fatal. Unsupported manostat in randomize_veloci&
 &ties(...). This is most likely an omission bug.'
      call fexit()
    end if
  end if
!
! fake-set ens%insK to avoid FPE
  ens%insK = 10.0
  call drift_removal(fr1,fr2)
!
! now re-scale to match exact temperature request
!
  ttt = 0.0
  ttg(:) = 0.0
  do imol=1,nmol
    do j=1,3
      temp = dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
      ttt = ttt + temp
      ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
    end do
    if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
      ttc = 3
    else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
      ttc = 5
!     WARNING: support missing
    else
      ttc = 6
      do j=4,6
        temp = (1.0/(RADIAN*RADIAN))*&
 &           dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
        ttt = ttt + temp
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
      end do    
    end if
    do j=ttc+1,ttc+ntormol(moltypid(imol))
      temp = (1.0/(RADIAN*RADIAN))*&
 &             dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
      ttt = ttt + temp
      ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
    end do
  end do
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      temp = pstat%params(1)*bnd_v*bnd_v
      ttt = ttt + temp
      ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + temp
      boxdof = 1
    else if (pstat%flag.eq.3) then
      boxdof = 0
    else
      write(ilog,*) 'Fatal. Unsupported manostat in randomize_veloci&
 &ties(...). This is most likely an omission bug.'
      call fexit()
    end if
  else
    boxdof = 0
  end if
  if ((totrbd+ndyntorsn+boxdof-n_constraints).le.0) then
    write(ilog,*) 'Attempting a calculation with no remaining degree&
 &s of freedom. This is fatal.'
    call fexit()
  end if
  if (ttt.le.0.0) then
    write(ilog,*) 'System is initially void of kinetic energy even t&
 &hough there are unconstrained degrees of freedom. This is either a&
 & bug or a numerical precision error due to an extremely small temp&
 &erature request. Fatal exit.'
    call fexit()
  end if
  do i=1,tstat%n_tgrps
    if (ttg(i).le.0.0) then
      write(ilog,*) 'At least on of the T-coupling groups is initial&
 &ly void of kinetic energy even though there are unconstrained degr&
 &ees of freedom. This is either a bug or a numerical precision erro&
 &r due to an extremely small temperature request. Fatal exit.'
      call fexit()
    end if
  end do
  ens%insK = (1.0/u_dyn_fconv)*0.5*ttt
  ttt = ttt / (dble(totrbd+ndyntorsn+boxdof-n_constraints)*u_dyn_kb)
  ttg(:) = ttg(:)/(u_dyn_kb*tstat%grpdof(:))
  corrfac(:) = sqrt(kelvin/ttg(:))
  do imol=1,nmol
    if (atmol(imol,1).eq.atmol(imol,2)) then
      intdim = 3
    else if (atmol(imol,1).eq.(atmol(imol,2)-1)) then
      intdim = 5
    else
      intdim = ntormol(moltypid(imol))+6
    end if
    do k=1,intdim
      dc_di(imol)%v(k) =corrfac(tstat%molgrp(imol))*dc_di(imol)%v(k)
    end do
  end do
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      bnd_v = corrfac(tstat%molgrp(nmol+1))*bnd_v
    end if
  end if
! 
end
!
!---------------------------------------------------------------------------
!
! the same for straight cartesian dynamics
!
subroutine randomize_cart_velocities()
!
  use system
  use forces
  use atoms
  use molecule
  use units
  use math
  use torsn
  use iounit
!
  implicit none
!
  integer imol,j,i,boxdof
  RTYPE vind,ttg(tstat%n_tgrps),normal,ttt,fr1,fr2
  RTYPE temp,corrfac(tstat%n_tgrps)
!
  do imol=1,nmol
!
!   straightforward
!
    do i=atmol(imol,1),atmol(imol,2)
      if (mass(i).le.0.0) cycle
      do j=1,3
        if (cart_frz(i,j).EQV..true.) then
          cart_v(i,j) = 0.0
          cycle
        end if
!       velocities in A/ps
        vind = normal()
        vind = vind*sqrt(u_dyn_kb*kelvin/mass(i))
        cart_v(i,j) = vind
      end do
    end do
!
    call update_comv(imol)
!
  end do
!
! box degrees of freedom (if applicable)
!
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      vind = normal()
      vind = vind*sqrt(u_dyn_kb*kelvin/pstat%params(1))
      bnd_v = vind
    else if (pstat%flag.eq.2) then
!     do nothing for Markov chain
    else if (pstat%flag.eq.3) then
      pstat%params(2) = 0.0
    else
      write(ilog,*) 'Fatal. Unsupported manostat in randomize_cart&
 &_velocities(...). This is most likely an omission bug.'
      call fexit()
    end if
  end if
!
! fake-set ens%insK to avoid FPE
  ens%insK = 10.0
  call drift_removal(fr1,fr2)
!
! now re-scale to match exact temperature request
!
  ttt = 0.0
  ttg(:) = 0.0
  do imol=1,nmol
    do i=atmol(imol,1),atmol(imol,2)
      do j=1,3
        temp = mass(i)*cart_v(i,j)*cart_v(i,j)
        ttt = ttt + temp
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
      end do
    end do
  end do
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      temp = pstat%params(1)*bnd_v*bnd_v
      ttt = ttt + temp
      ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + temp
      boxdof = 1
    else if (pstat%flag.eq.2) then
      boxdof = 0
    else if (pstat%flag.eq.3) then
      boxdof = 0
    else
      write(ilog,*) 'Fatal. Unsupported manostat in randomize_cart&
 &_velocities(...). This is most likely an omission bug.'
      call fexit()
    end if
  else
    boxdof = 0
  end if
  if ((3*n+boxdof-n_constraints).lt.0) then
    write(ilog,*) 3*n,boxdof,n_constraints
    write(ilog,*) 'Attempting a calculation with no remaining degree&
 &s of freedom. This is fatal.'
    call fexit()
  end if
  if (ttt.le.0.0) then
    write(ilog,*) 'System is initially void of kinetic energy even t&
 &hough there are unconstrained degrees of freedom. This is either a&
 & bug or a numerical precision error due to an extremely small temp&
 &erature request. Fatal exit.'
    call fexit()
  end if
  do i=1,tstat%n_tgrps
    if (ttg(i).le.0.0) then
      write(ilog,*) 'At least on of the T-coupling groups is initial&
 &ly void of kinetic energy even though there are unconstrained degr&
 &ees of freedom. This is either a bug or a numerical precision erro&
 &r due to an extremely small temperature request. Fatal exit.'
      call fexit()
    end if
  end do
  ens%insK = (1.0/u_dyn_fconv)*0.5*ttt
  ttt = ttt / (dble(3*n+boxdof-n_constraints)*u_dyn_kb)
  ttg(:) = ttg(:)/(u_dyn_kb*tstat%grpdof(:))
  corrfac(:) = sqrt(kelvin/ttg(:))
  do imol=1,nmol
    do i=atmol(imol,1),atmol(imol,2)
      do j=1,3
        if (cart_frz(i,j).EQV..true.) cycle
        cart_v(i,j) = corrfac(tstat%molgrp(imol))*cart_v(i,j)
      end do
    end do
  end do
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      bnd_v = corrfac(tstat%molgrp(nmol+1))*bnd_v
    end if
  end if
! 
end
!
!---------------------------------------------------------------------------
!
!  a subroutine to simply re-scale velocities by a T-increment
!
subroutine rescale_velocities(foreign_T)
!
  use system
  use forces
  use molecule
  use iounit
  use atoms
!
  implicit none
!
  integer imol,j,ttc
  RTYPE corrfac,foreign_T
!
  corrfac = sqrt(kelvin/foreign_T)
!
  do imol=1,nmol
!
!   first: center-of-mass translation
!
    do j=1,3
!     velocities in A/ps
      dc_di(imol)%v(j) = corrfac*dc_di(imol)%v(j)
    end do
!
    if ((atmol(imol,2)-atmol(imol,1)).eq.0) cycle
!
!   second: rigid-body rotation

    if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
!
      ttc = 5
      write(ilog,*) 'Fatal. Diatomic molecules currently not support&
 &ed by dynamics implementation.'
      call fexit()
!     WARNING: support missing
!
    else
!
      ttc = 6
      do j=4,6
!       velocities in deg/ps
        dc_di(imol)%v(j) = corrfac*dc_di(imol)%v(j)
      end do
!
    end if
!
!   third: torsional degrees of freedom
!
    do j=1,ntormol(moltypid(imol))
      dc_di(imol)%v(ttc+j) = corrfac*dc_di(imol)%v(ttc+j) 
    end do
!
  end do
!
! fourth: box degrees of freedom (if applicable)
!
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      bnd_v = corrfac*bnd_v
    else if (pstat%flag.eq.3) then
      pstat%params(2) = corrfac*pstat%params(2)
    else
      write(ilog,*) 'Fatal. Unsupported manostat in rescale_veloci&
 &ties(...). This is most likely an omission bug.'
      call fexit()
    end if
  end if
!
end
!
!-----------------------------------------------------------------------
!
! the same for straight cartesian dynamics
!
subroutine rescale_cart_velocities(foreign_T)
!
  use system
  use forces
  use molecule
  use iounit
  use atoms
!
  implicit none
!
  integer imol
  RTYPE foreign_T,corrfac
!
  corrfac = sqrt(kelvin/foreign_T)
!
  cart_v(1:n,1:3) = corrfac*cart_v(1:n,1:3)
  do imol=1,nmol
    call update_comv(imol)
  end do
!
! box degrees of freedom (if applicable)
!
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      bnd_v = corrfac*bnd_v
    else if (pstat%flag.eq.2) then
!     do nothing for Markov chain
    else if (pstat%flag.eq.3) then
      pstat%params(2) = corrfac*pstat%params(2)
    else
      write(ilog,*) 'Fatal. Unsupported manostat in rescale_cart&
 &_velocities(...). This is most likely an omission bug.'
      call fexit()
    end if
  end if
!
end
!
!-----------------------------------------------------------------------
!
! the subroutine to remove drift of the center of mass of the whole
! system (either as translation or rotation)
! not that drift is only computed over the intermolecular degrees of
! freedom
! in non-frictional dynamics, drift can easily absorb all of the kinetic
! energy of the system while effectively freezing the "real" degrees
! of freedom
!
subroutine drift_removal(fr1,fr2)
!
  use iounit
  use system
  use molecule
  use forces
  use mcsums
  use units
  use atoms
!
  implicit none
!
  integer j,imol,i
  RTYPE vavg(6),fr1,fr2,mte,ivmm,subm
!
 59   format("Fraction of drift translation: ",g11.5)
 60   format("Fraction of drift rotation   : ",g11.5)
!
  if (n.eq.1) then
    fr1 = 0.0
    fr2 = 0.0
    return
  end if
  if ((nmol.eq.1).AND.(fycxyz.eq.1)) then
    fr1 = 0.0
    fr2 = 0.0
    return
  end if
!
  do j=1,5
    vavg(j) = 0.0
  end do
  do imol=1,nmol
!
    subm = 0.0
    do j=atmol(imol,1),atmol(imol,2)
      if ((cart_frz(j,1).EQV..true.).AND.(cart_frz(j,2).EQV..true.).AND.(cart_frz(j,3).EQV..true.)) then
        subm = subm + mass(j)
      end if
    end do
    do j=1,3
      vavg(j) = vavg(j) + molmass(moltypid(imol))*dc_di(imol)%v(j)
    end do
    vavg(4) = vavg(4) + molmass(moltypid(imol)) - subm
!
  end do
  do j=1,3
    vavg(j) = vavg(j)/vavg(4)
  end do
  vavg(5) = vavg(1)*vavg(1) + vavg(2)*vavg(2) + vavg(3)*vavg(3)
  fr1 = 0.5*vavg(4)*vavg(5)/(u_dyn_fconv*ens%insK)
  if (ens%sysfrz.ge.2) then
    do imol=1,nmol
      if ((n.eq.1).OR.((nmol.eq.1).AND.(fycxyz.eq.1))) exit
      if (fycxyz.eq.1) then
        do j=1,3
          dc_di(imol)%v(j) = dc_di(imol)%v(j) - vavg(j)
        end do
      else if (fycxyz.eq.2) then
        ivmm = 1.0/molmass(moltypid(imol))
!       note that this is the same as writing it out over atoms in the first place instead
!       of taking the detour through c.o.m. velocities
        do i=atmol(imol,1),atmol(imol,2)
          mte = mass(i)*ivmm
          do j=1,3
            cart_v(i,j) = cart_v(i,j) - mte*vavg(j)
          end do
        end do
      end if
    end do
  end if
!
  if ((fycxyz.eq.1).AND.(nmol.gt.1)) then
    call rot_system_int(fr2)
  else if ((fycxyz.eq.2).AND.(n.gt.1)) then
    call rot_system_cart(fr2)
  else
    write(ilog,*) 'Fatal. Encountered unsupported choice for degrees&
 & of freedom in drift_removal(...). This is a bug.'
    call fexit()
  end if
  fr2 = fr2/ens%insK
!
end
!
!-----------------------------------------------------------------------
!
! note that this function cannot remove "internal" rotational drift
! potentially dangerous --> WARNING!!!!!!!!!!!
!
subroutine rot_system_int(erot)
!
  use iounit
  use forces
  use atoms
  use molecule
  use sequen
  use math
  use mcsums
  use system
  use units
!
  implicit none
!
  integer imol,i,j,athree
  RTYPE pos2(3),pos3(3),pos1(3),tst,xr,yr,zr,tem
  RTYPE or_pl(3),nopl,c2,totm,erot
  RTYPE jte(3,3),jti(3,3),vel(3),vel2(3)
  logical afalse
!
  afalse = .false.
  athree = 3
!
  do j=1,3
    pos3(j) = 0.0
    vel(j) = 0.0
    do i=1,3
      jte(j,i) = 0.0
    end do
  end do
  totm = 0.0
  do imol=1,nmol
    totm = totm + molmass(moltypid(imol))
    do j=1,3
      pos3(j) = pos3(j) + molmass(moltypid(imol))*comm(imol,j)
    end do
  end do
  do j=1,3
    pos3(j) = pos3(j)/totm
  end do
!
  if (nstep.eq.1) then
    tst = 0.5*dyn_dt
  else
    tst = dyn_dt
  end if
! get net angular momentum in A^2 * g / (mol * ps)
  do imol=1,nmol
    do j=1,3
      pos2(j) = molmass(moltypid(imol))*(comm(imol,j)-pos3(j))
    end do
    c2 = pos2(1)*pos2(1) + pos2(2)*pos2(2) + pos2(3)*pos2(3)
    do j=1,3
      pos1(j) = dc_di(imol)%v(j)
    end do
    call crossprod(pos2,pos1,or_pl,nopl)
    do j=1,3
      vel(j) = vel(j) + or_pl(j)
    end do
  end do
!  write(*,*) vel(1:3)
! now get inertial tensor in A^2 * g / mol
  if (nmol.gt.2) then
    do imol=1,nmol
      xr = comm(imol,1)-pos3(1)
      yr = comm(imol,2)-pos3(2)
      zr = comm(imol,3)-pos3(3)
      jte(1,1) = jte(1,1) + molmass(moltypid(imol))*(yr*yr + zr*zr)
      jte(1,2) = jte(1,2) - molmass(moltypid(imol))*xr*yr
      jte(1,3) = jte(1,3) - molmass(moltypid(imol))*xr*zr
      jte(2,2) = jte(2,2) + molmass(moltypid(imol))*(xr*xr + zr*zr)
      jte(2,3) = jte(2,3) - molmass(moltypid(imol))*yr*zr
      jte(3,3) = jte(3,3) + molmass(moltypid(imol))*(xr*xr + yr*yr)
    end do
    jte(3,2) = jte(2,3)
    jte(3,1) = jte(1,3)
    jte(2,1) = jte(1,2)
!   invert it
    call invmat(jte,athree,jti)  
!   derive angular velocity in 1 / ps
    do i=1,3
      vel2(i) = 0.0
      do j=1,3
        vel2(i) = vel2(i) + jti(i,j)*vel(j)
      end do
    end do
  else if (nmol.eq.2) then
    tem = (comm(1,1) - comm(2,1))**2 +&
 &             (comm(1,2) - comm(2,2))**2 +&
 &             (comm(1,3) - comm(2,3))**2 
    tem = tem*molmass(moltypid(1))*molmass(moltypid(2))/totm
    do j=1,3
      vel2(j) = vel(j)/tem
    end do
  end if
  erot = 0.0
  do i=1,3
     erot = erot + vel2(i)*vel(i)
  end do
  erot = 0.5*erot/u_dyn_fconv
! and finally remove by adjusting c.o.m.-velocties
  if (ens%sysfrz.eq.3) then
    do imol=1,nmol
      xr = comm(imol,1)-pos3(1)
      yr = comm(imol,2)-pos3(2)
      zr = comm(imol,3)-pos3(3)
      dc_di(imol)%v(1) = dc_di(imol)%v(1) - vel2(2)*zr + vel2(3)*yr
      dc_di(imol)%v(2) = dc_di(imol)%v(2) - vel2(3)*xr + vel2(1)*zr
      dc_di(imol)%v(3) = dc_di(imol)%v(3) - vel2(1)*yr + vel2(2)*xr
    end do
  end if
!
end
!
!-----------------------------------------------------------------------
!
! the same for full Cartesian description (simpler and slower)
!
subroutine rot_system_cart(erot)
!
  use iounit
  use forces
  use atoms
  use molecule
  use sequen
  use math
  use mcsums
  use system
  use units
!
  implicit none
!
  integer i,j,athree
  RTYPE pos2(3),pos3(3),pos1(3),tst,xr,yr,zr,tem
  RTYPE or_pl(3),nopl,c2,totm,erot
  RTYPE jte(3,3),jti(3,3),vel(3),vel2(3)
  logical afalse
!
  afalse = .false.
  athree = 3
!
  pos3(:) = 0.0
  vel(:) = 0.0
  jte(:,:) = 0.0
  totm = 0.0
!
  do i=1,n
    totm = totm + mass(i)
    pos3(1) = pos3(1) + mass(i)*x(i)
    pos3(2) = pos3(2) + mass(i)*y(i)
    pos3(3) = pos3(3) + mass(i)*z(i)
  end do
  pos3(:) = pos3(:)/totm
!
  if (nstep.eq.1) then
    tst = 0.5*dyn_dt
  else
    tst = dyn_dt
  end if
! get net angular momentum in A^2 * g / (mol * ps)
  do i=1,n
    pos2(1) = mass(i)*(x(i)-pos3(1))
    pos2(2) = mass(i)*(y(i)-pos3(2))
    pos2(3) = mass(i)*(z(i)-pos3(3))
    c2 = pos2(1)*pos2(1) + pos2(2)*pos2(2) + pos2(3)*pos2(3)
    do j=1,3
      pos1(j) = cart_v(i,j)
    end do
    call crossprod(pos2,pos1,or_pl,nopl)
    vel(:) = vel(:) + or_pl(:)
  end do
! now get inertial tensor in A^2 * g / mol
  if (n.gt.2) then
    do i=1,n
      xr = x(i)-pos3(1)
      yr = y(i)-pos3(2)
      zr = z(i)-pos3(3)
      jte(1,1) = jte(1,1) + mass(i)*(yr*yr + zr*zr)
      jte(1,2) = jte(1,2) - mass(i)*xr*yr
      jte(1,3) = jte(1,3) - mass(i)*xr*zr
      jte(2,2) = jte(2,2) + mass(i)*(xr*xr + zr*zr)
      jte(2,3) = jte(2,3) - mass(i)*yr*zr
      jte(3,3) = jte(3,3) + mass(i)*(xr*xr + yr*yr)
    end do
    jte(3,2) = jte(2,3)
    jte(3,1) = jte(1,3)
    jte(2,1) = jte(1,2)
!   invert it
    call invmat(jte,athree,jti)
!   derive angular velocity in 1 / ps
    do i=1,3
      vel2(i) = 0.0
      do j=1,3
        vel2(i) = vel2(i) + jti(i,j)*vel(j)
      end do
    end do
  else if (n.eq.2) then
    tem = (x(1) - x(2))**2 +&
 &             (y(1) - y(2))**2 +&
 &             (z(1) - z(2))**2 
    tem = tem*mass(1)*mass(2)/totm
    do j=1,3
      vel2(j) = vel(j)/tem
    end do
  end if
  erot = 0.0
  do i=1,3
     erot = erot + vel2(i)*vel(i)
  end do
  erot = 0.5*erot/u_dyn_fconv
! and finally remove by adjusting c.o.m.-velocties
  if (ens%sysfrz.eq.3) then
    do i=1,n
      xr = x(i)-pos3(1)
      yr = y(i)-pos3(2)
      zr = z(i)-pos3(3)
      cart_v(i,1) = cart_v(i,1) - vel2(2)*zr + vel2(3)*yr
      cart_v(i,2) = cart_v(i,2) - vel2(3)*xr + vel2(1)*zr
      cart_v(i,3) = cart_v(i,3) - vel2(1)*yr + vel2(2)*xr
    end do
  end if
!
end
!
!-----------------------------------------------------------------------
!
! here we compute the (internal, potential) pressure tensor from the virial     
!
subroutine virial(virpot)
!
  use iounit
  use system
  use forces
  use atoms
  use molecule
!
  implicit none
!
  integer i,j,k
  RTYPE virpot(3,3)
!
! basically, every force which acts along a distance or directly on a positional
! coordinate will affect the virial
!
  if ((bnd_type.eq.3).OR.(bnd_type.eq.4)) then
    if (bnd_shape.ne.2) then
      write(ilog,*) 'Fatal. Encountered unsupported box shape in virial calculation. This is an omission bug.'
      call fexit()
    end if
!   this is really simple assuming we have all the relevant forces collected in cart_f
!   note that this includes boundary forces
    do j=1,3
      do k=1,3
        virpot(k,j) = 0.0
      end do
    end do
    do i=1,nmol
      do j=1,3
!       note that we compute the virial over the previous configuration, since all forces
!       are computed for that
        virpot(j,1) = virpot(j,1) + dc_di(i)%f(j)*&
 & (commref(i,j)-bnd_params(1))
        virpot(j,2) = virpot(j,2) + dc_di(i)%f(j)*&
 & (commref(i,j)-bnd_params(2))
        virpot(j,3) = virpot(j,3) + dc_di(i)%f(j)*&
 & (commref(i,j)-bnd_params(3))
      end do
    end do
    do i=1,3
      do j=1,3
!       subtract out the boundary term (which is what maintains mechanical equilibrium and
!       averages the pressure to be zero)
        virpot(i,j) = virpot(i,j) + bnd_f(i,j)
      end do
    end do
  else if (bnd_type.eq.1) then
!   with periodic images it is much more complicated
!   if we can adopt the single-sum implementation as shown in GROMACS it might be ok;
!   otherwise collecting on-the-fly becomes inevitable (position vectors need to be
!   minimum image-corrected ...)
    write(ilog,*) 'Fatal. Pressure calculation currently does not support periodic boundary conditions. Check back later.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine gen_cartv(mshfs)
!
  use iounit
  use system
  use atoms
  use molecule
  use mcsums
  use forces
!
  implicit none
!
  integer imol,k
  RTYPE mshfs(3,nmol),tst
!
! PBC
  if (bnd_type.eq.1) then
    if (nstep.eq.1) then
!     reconstructed velocities are techniqually for the half-step t + 1/2dt and not for time t
!     (as would be desired by most ensemble routines)
      cart_v(1:n,1) = (x(1:n) - xref(1:n))/dyn_dt
      cart_v(1:n,2) = (y(1:n) - yref(1:n))/dyn_dt
      cart_v(1:n,3) = (z(1:n) - zref(1:n))/dyn_dt
      do imol=1,nmol
        do k=1,3
          tst = mshfs(k,imol)/dyn_dt
          cart_v(atmol(imol,1):atmol(imol,2),k) = cart_v(atmol(imol,1):atmol(imol,2),k) - tst
        end do
      end do
    else
!     if we have an old cart_v, we can estimate the velocity at timestep t by averaging
!     over the velocities at times t - 1/2dt and t + 1/2dt
      cart_v(1:n,1) = 0.5*(cart_v(1:n,1)+(x(1:n)-xref(1:n))/dyn_dt)
      cart_v(1:n,2) = 0.5*(cart_v(1:n,2)+(y(1:n)-yref(1:n))/dyn_dt)
      cart_v(1:n,3) = 0.5*(cart_v(1:n,3)+(z(1:n)-zref(1:n))/dyn_dt)
      do imol=1,nmol
        do k=1,3
          tst = 0.5*mshfs(k,imol)/dyn_dt
          cart_v(atmol(imol,1):atmol(imol,2),k) = cart_v(atmol(imol,1):atmol(imol,2),k) - tst
        end do
      end do
    end if
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
    if (nstep.eq.1) then
!     reconstructed velocities are techniqually for the half-step t + 1/2dt and not for time t
!     (as would be desired by most ensemble routines)
      cart_v(1:n,1) = (x(1:n) - xref(1:n))/dyn_dt
      cart_v(1:n,2) = (y(1:n) - yref(1:n))/dyn_dt
      cart_v(1:n,3) = (z(1:n) - zref(1:n))/dyn_dt
    else
!     if we have an old cart_v, we can estimate the velocity at timestep t by averaging
!     over the velocities at times t - 1/2dt and t + 1/2dt
      cart_v(1:n,1) = 0.5*(cart_v(1:n,1)+(x(1:n)-xref(1:n))/dyn_dt)
      cart_v(1:n,2) = 0.5*(cart_v(1:n,2)+(y(1:n)-yref(1:n))/dyn_dt)
      cart_v(1:n,3) = 0.5*(cart_v(1:n,3)+(z(1:n)-zref(1:n))/dyn_dt)
    end if
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition&
 & and/or box shape in gen_cartv(...). Please check back later.'
    call fexit()
  end if 
!
end
!
!-----------------------------------------------------------------------
!
! here we compute the (ideal, kinetic) pressure tensor from the (Cartesian) velocities
!     
subroutine ekinten(virkin)
!
  use iounit
  use system
  use forces
  use atoms
  use mcsums
  use molecule
!
  implicit none
!
  integer j,imol
  RTYPE virkin(3,3)
!
  if (fycxyz.eq.1) then
    virkin(:,:) = 0.0
    do j=1,3
      do imol=1,nmol
        virkin(j,1)=virkin(j,1) + dc_di(imol)%v(j)*dc_di(imol)%v(1)*&
 &                     molmass(moltypid(imol))
        virkin(j,2)=virkin(j,2) + dc_di(imol)%v(j)*dc_di(imol)%v(2)*&
 &                      molmass(moltypid(imol))
        virkin(j,3)=virkin(j,3) + dc_di(imol)%v(j)*dc_di(imol)%v(3)*&
 &                      molmass(moltypid(imol))
      end do
    end do
  else if (fycxyz.eq.2) then
!   these are half-step velocities and hence slightly incorrect
    do j=1,3
      virkin(j,1) = sum(cart_v(1:n,j)*cart_v(1:n,1)*mass(1:n))
      virkin(j,2) = sum(cart_v(1:n,j)*cart_v(1:n,2)*mass(1:n))
      virkin(j,3) = sum(cart_v(1:n,j)*cart_v(1:n,3)*mass(1:n))
    end do
  else
    write(ilog,*) 'Fatal. Encountered unsupported choice of degrees &
 &of freedom in ekinten(...). This is an omission bug.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine thermostat(tsc)
!
  use system
  use iounit
  use units
  use forces
  use system
  use atoms
  use molecule
  use math
!
  implicit none
!
  RTYPE tsc(tstat%n_tgrps),random,normal,vind,decayrate,expt,exp2t,tgtt,rn1,gn1,rndgamma
  integer j,imol,i,ttc
  logical didup
!
  decayrate = dyn_dt/tstat%params(1)
!
  if (tstat%flag.eq.1) then
!   Berendsen is really simple
    tsc(:) = sqrt( 1.0 + decayrate*&
 &                ((kelvin/tstat%grpT(:)) - 1.0) )
  else if (tstat%flag.eq.2) then
!   Andersen it also simple: note that we do not use tsc() here, but modify velocities
!   directly, and that we ignore T-groups (since they are irrelevant, but can still be used
!   for analysis!)
    tsc(:) = 1.0
    if (fycxyz.eq.1) then
!     in non-Cartesian dynamics we have the problem of vastly inhomogeneous degrees of freedom,
!     and we'll therefore couple each d.o.f. independently 
!     we also assume the inertial masses are xyz up-to-date (T-stat called soon after cart2int()!)
      do imol=1,nmol
!       first: center-of-mass translation
        do j=1,3
!         velocities in A/ps
          if (dc_di(imol)%frz(j).EQV..true.) then
            cycle
          end if
          if (random().lt.decayrate) then
            vind = normal()
            vind = vind*sqrt(u_dyn_kb*kelvin/dc_di(imol)%im(j))
            dc_di(imol)%v(j) = vind
          end if
        end do
        if ((atmol(imol,2)-atmol(imol,1)).eq.0) cycle
!       second: rigid-body rotation
        if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
          call fexit()
!         WARNING: support missing
        else
          ttc = 6
          do j=4,6
!           velocities in deg/ps
            if (dc_di(imol)%frz(j).EQV..true.) then
              cycle
            end if
            if (random().lt.decayrate) then
              vind = normal()
              vind = vind*sqrt(u_dyn_kb*kelvin/dc_di(imol)%im(j))
              dc_di(imol)%v(j) = RADIAN*vind
            end if
          end do
        end if
!       third: torsional degrees of freedom
        do j=1,ntormol(moltypid(imol))
          if (dc_di(imol)%frz(ttc+j).EQV..true.) then
            cycle
          end if
          if (random().lt.decayrate) then
            vind = normal()
            vind = vind*sqrt(u_dyn_kb*kelvin/dc_di(imol)%im(ttc+j))
            dc_di(imol)%v(ttc+j) = RADIAN*vind
          end if
        end do
!       fourth: box degrees of freedom (if applicable)
        if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
          if (random().lt.decayrate) then
            vind = normal()
            vind = vind*sqrt(u_dyn_kb*kelvin/pstat%params(1))
            bnd_v = vind
          end if
        end if
      end do
    else if (fycxyz.eq.2) then
!     in Cartesian dynamics we'll simply couple each atom individually
      do imol=1,nmol
        didup = .false.
        do i=atmol(imol,1),atmol(imol,2)
          if (mass(i).le.0.0) cycle
          if (random().lt.decayrate) then
            didup = .true.
            do j=1,3
!             velocities in A/ps
              vind = normal()
              vind = vind*sqrt(u_dyn_kb*kelvin/mass(i))
              cart_v(i,j) = vind
            end do
          end if
        end do
        if (didup.EQV..true.) then
          call update_comv(imol)
        end if
      end do
    end if
! for the Nose-Hoover analog method (by Stern), we'll simply set this to what is needed for the integrator
  else if (tstat%flag.eq.3) then
    tsc(:) = exp(-tstat%params(2)*0.5*dyn_dt)
! Bussi is really simple as well - note they strongly recommend using the exact (analytical) form derived in the appendix
! if all Gaussian random numbers are zero, and the gamma random # is ~ grpdof, analysis of leading terms recovers Berendsen
! conversely, if decayrate -> 0.0, NVE is recovered
  else if (tstat%flag.eq.4) then
    expt = exp(-decayrate)
    exp2t = 2.0*exp(-0.5*decayrate)
    do i=1,tstat%n_tgrps
      tgtt = kelvin/(tstat%grpT(i)*tstat%grpdof(i))
!     note that the factor of 2.0 in front of rndgamma is not immediately apparent from the reference
!     but is crucial for relating sums of squared Gaussian RNs to values pulled from the gamma-dist.
      if (tstat%grpdof(i).eq.1) then
        gn1 = 0.0
      else if (tstat%grpdof(i).eq.2) then
        rn1 = normal()
        gn1 = rn1*rn1
      else if (mod(floor(tstat%grpdof(i)+1.0e-7)-1,2).eq.0) then
        gn1 = 2.0*rndgamma((1.0*tstat%grpdof(i)-1.0)/2.0)
!     this split may not be necessary
      else
        rn1 = normal()
        gn1 = 2.0*rndgamma((1.0*tstat%grpdof(i)-2.0)/2.0) + rn1*rn1
      end if
      rn1 = normal()
      tsc(i) = sqrt(expt + tgtt*(1.0-expt)*(rn1*rn1 + gn1) + exp2t*rn1*sqrt(tgtt*(1.0-expt)))
    end do
  else
    write(ilog,*) 'Encountered unsupported thermostat flag. Offendin&
 &g code is ',tstat%flag,'. Please report this problem!'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine init_thermostat()
!
  use system
  use iounit
  use molecule
  use forces
  use torsn
  use atoms
  use shakeetal
  use sequen
!
  implicit none
!
  integer ttg(tstat%n_tgrps),i,ttf(tstat%n_tgrps),ttc,imol,j,boxdof
  integer drift_c,alldof
  logical samxyz(3)
!
  ttg(:) = 0
  ttf(:) = 0
  boxdof = 0
  drift_c = 0
! note that we have to spread out the drift-correction constraints onto
! the T-coupling groups
  if ((nmol.gt.1).AND.(fycxyz.eq.1)) then
    if (ens%sysfrz.eq.3) then
!     if there's exactly two molecules one rotational dof is missing
      if (nmol.eq.2) then
        drift_c = 5
      else
        drift_c = 6
      end if
    else if (ens%sysfrz.eq.2) then
      drift_c = 3
    else
      drift_c = 0
    end if
  else if (fycxyz.eq.1) then
    drift_c = 0 ! drift_c are lumped in dc_di%frz
  else
    samxyz(:) = .false.
    do i=1,n
      if (cart_frz(i,1).EQV..false.) samxyz(1) = .true.
      if (cart_frz(i,2).EQV..false.) samxyz(2) = .true.
      if (cart_frz(i,3).EQV..false.) samxyz(3) = .true.
    end do
    j = 0
    do i=1,3
      if (samxyz(i).EQV..false.) j = j + 1 ! j can never be 3
    end do
    if (ens%sysfrz.eq.3) then
!     if there's exactly two molecules one rotational dof is missing
      if (n.eq.2) then
        drift_c = 5 - j - min(j,1)*j
      else if (n.eq.1) then
        drift_c = 3 - j
      else
        drift_c = 6 - j - min(j,1)*(j+1)
      end if
    else if (ens%sysfrz.eq.2) then
      drift_c = 3 - j
    end if
  end if
  if (fycxyz.eq.1) then
    do imol=1,nmol
      do j=1,3
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + 1
        if (dc_di(imol)%frz(j).EQV..true.) then
          ttf(tstat%molgrp(imol)) = ttf(tstat%molgrp(imol)) + 1
        end if
      end do
      if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
        ttc = 3
      else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        ttc = 5
        write(ilog,*) 'Fatal. Diatomic molecules are not yet support&
 &ed in dynamics. Check back later.'
        call fexit()
!       WARNING: support missing
      else
        ttc = 6
        do j=4,6
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + 1
          if (dc_di(imol)%frz(j).EQV..true.) then
            ttf(tstat%molgrp(imol)) = ttf(tstat%molgrp(imol)) + 1
          end if
        end do    
      end if
      do j=ttc+1,ttc+ntormol(moltypid(imol))
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + 1
        if (dc_di(imol)%frz(j).EQV..true.) then
          ttf(tstat%molgrp(imol)) = ttf(tstat%molgrp(imol)) + 1
        end if
      end do
    end do
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      if (pstat%flag.eq.1) then
        ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + 1
        boxdof = 1
      end if
    end if
    alldof = totrbd+ndyntorsn-n_constraints+boxdof
    if ((sum(ttg)-sum(ttf)-drift_c).ne.alldof) then
      write(ilog,*) 'Fatal. Setup of degrees of freedom is inconsist&
 &ent. This is most certainly a bug.'
       call fexit()
    end if
  else if (fycxyz.eq.2) then
    do imol=1,nmol
      do i=atmol(imol,1),atmol(imol,2)
        do j=1,3
          if (cart_frz(i,j).EQV..false.) ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + 1
        end do
      end do
    end do
    do i=1,cart_cons_grps
      if (constraints(i)%nr.gt.0) then
        imol = molofrs(atmres(constraints(i)%idx(1,1)))
        ttf(tstat%molgrp(imol)) = ttf(tstat%molgrp(imol)) + constraints(i)%nr
      end if
      if (constraints(i)%nr3.gt.0) then
        imol = molofrs(atmres(constraints(i)%idx3(1,1))) ! the same
        ttf(tstat%molgrp(imol)) = ttf(tstat%molgrp(imol)) + constraints(i)%nr3
      end if
      if (constraints(i)%nr4.gt.0) then
        imol = molofrs(atmres(constraints(i)%idx4(1,1))) ! the same
        ttf(tstat%molgrp(imol)) = ttf(tstat%molgrp(imol)) + constraints(i)%nr4
      end if
    end do
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      if (pstat%flag.eq.1) then
        ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + 1
        boxdof = 1
      end if
    end if
    alldof = 3*n-n_constraints+boxdof
    if ((sum(ttg)-sum(ttf)-drift_c).ne.alldof) then
      write(ilog,*) 'Fatal. Setup of degrees of freedom is inconsist&
 &ent. This is most certainly a bug.'
      call fexit()
    end if
  else
    write(ilog,*) 'Fatal. Encountered unsupported choice of degrees &
 &of freedom in init_thermostat(). This is an omission bug.'
    call fexit()
  end if
!
  do i=1,tstat%n_tgrps
    tstat%grpT(i) = 0.0
    tstat%grpdof(i) = ttg(i) - ttf(i) - &
 &    dble(drift_c*(ttg(i) - ttf(i)))/dble(alldof+drift_c)
    if (tstat%grpdof(i).le.0.0) then
      write(ilog,*) 'Fatal. T-coupling group #',i,' has no active de&
 &grees of freedom. Check input files for unwanted constraints.'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine manostat(boxovol)
!
  use system
  use iounit
  use units
  use forces
  use system
  use energies
  use molecule
  use math
  use atoms
!
  implicit none
!
  RTYPE netf,deltaV,boxovol,oldbnd(3),evn(MAXENERGYTERMS)
  RTYPE random,scf,diff
  logical atrue,afalse
!
  atrue = .true.
  afalse = .false.
  bnd_params(8) = 1000.0
!
! in constant pressure simulations update box
  if (bnd_shape.eq.2) then
    if (pstat%flag.eq.1) then
      netf= bnd_fr - 4.0*PI*bnd_params(5)*extpress/u_dyn_virconv
      bnd_v= bnd_v + 0.5*u_dyn_fconv*dyn_dt*netf/pstat%params(1)
      bnd_params(4) = bnd_params(4) + dyn_dt*bnd_v
      call update_bound(afalse)
    else if (pstat%flag.eq.2) then
!     backup old radius
      oldbnd(1) = bnd_params(4)
      deltaV = (random()-0.5)*bnd_params(8)
      bnd_params(4) = ((3./(4.0*PI))*(boxovol+deltaV))**(1./3.)
      call update_bound(afalse)
      call e_boundary(evn)
      diff = 3.0*(evn(12)-esterms(12)) + &
 &                       3.0*deltaV*extpress/u_dyn_virconv
      invtemp = 1./(gasconst*ens%insT)
      if ((diff.le.0.0).OR.(random().lt.exp(-diff*invtemp))) then
        esterms(12) = evn(12)
      else
        bnd_params(4) = oldbnd(1)
        call update_bound(afalse)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported manostat for cho&
 &sen box and ensemble (manostat-flag: ',pstat%flag,').'
      call fexit()
    end if
  else if (bnd_shape.eq.1) then
    if (pstat%flag.eq.2) then
!     backup old radius
      oldbnd(1:3) = bnd_params(1:3)
      deltaV = (random()-0.5)*bnd_params(8)
      if ((boxovol+deltaV).le.0.0) then
        write(ilog,*) 'Fatal. Volume fluctuations yielded zero or ne&
 &gative volume in Markov chain. This indicates a poorly set up simu&
 &lation. Please check settings!'
        call fexit()
      end if
!     isotropic
      scf = ((boxovol+deltaV)**(1./3.))/bnd_params(1)
      bnd_params(1:3) = scf*bnd_params(1:3)
      call update_bound(afalse)
      diff = -3.0*deltaV*extpress/u_dyn_virconv
      invtemp = 1./(gasconst*ens%insT)
      if ((diff.le.0.0).OR.(random().lt.exp(-diff*invtemp))) then
        call rescale_xyz(scf,atrue)
      else
        bnd_params(1:3) = oldbnd(1:3)
        call update_bound(afalse)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported manostat for cho&
 &sen box and ensemble (manostat-flag: ',pstat%flag,').'
      call fexit()
    end if
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary shape for&
 & chosen ensemble (Flag: ',ens%flag,').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
