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
!
!------------------------------------------------------------------------
!
! routine for the smooth DSSP-based global potential energy term
! note that beta-bridges can very well stretch multiple molecules such
! that the DSSP-assignment is best done globally
! note that the DSSP-potential is written over positions, not torsions
! and therefore generates a Cartesian derivative (torsional derivative
! via cart2_int())
!
subroutine force_dssp_gl(evec,ca_f)
!
  use molecule
  use iounit
  use dssps
  use forces
  use energies
  use polypep
  use atoms
!
  implicit none
!
  integer i,ii
  RTYPE esc,wesc,hsc,whsc,evec(MAXENERGYTERMS),terme,termh
  RTYPE escm(npepmol),wescm(npepmol),hscm(npepmol),whscm(npepmol),ca_f(n,3)
  logical atrue
!
  atrue = .true.
  call update_dssp_der(esc,wesc,hsc,whsc,escm,wescm,hscm,whscm,atrue)
!
  evec(19) = evec(19) + scale_DSSP*&
 &            (par_DSSP(8)*(whsc - par_DSSP(7))**2 +&
 &             par_DSSP(10)*(wesc - par_DSSP(9))**2)
  terme = scale_DSSP*2.0*par_DSSP(10)*(wesc - par_DSSP(9))
  termh = scale_DSSP*2.0*par_DSSP(8)*(whsc - par_DSSP(7))
  do i=1,pep_sz
    ii = rspep(i)
    if (oi(ii).gt.0) ca_f(oi(ii),:) = ca_f(oi(ii),:) - termh*whsc_der(:,1,i) - terme*wesc_der(:,1,i)
    if (ci(ii).gt.0) ca_f(ci(ii),:) = ca_f(ci(ii),:) - termh*whsc_der(:,2,i) - terme*wesc_der(:,2,i)
    if (ni(ii).gt.0) ca_f(ni(ii),:) = ca_f(ni(ii),:) - termh*whsc_der(:,3,i) - terme*wesc_der(:,3,i)
    if (hni(ii).gt.0) ca_f(hni(ii),:) = ca_f(hni(ii),:) - termh*whsc_der(:,4,i) - terme*wesc_der(:,4,i)
  end do
!
end
!
!------------------------------------------------------------------------------
!
! forces for the smooth density restraint potential
!
subroutine force_emicro_gl(evec,ca_f)
!
  use ems
  use energies
  use units
  use atoms
  use mcsums
  use params
  use system
  use molecule
  use grandensembles
!
  implicit none
!
  RTYPE evec(MAXENERGYTERMS)
  RTYPE accum,normer,normer2
  RTYPE unitvol,soldens,temj,temk,teml,term0,ca_f(n,3)
  integer ints(3),intsm1(3),xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl,k,j,l,sh,imol,ii,cc,lidx
  logical afalse,doslice,doblocks,doline,ismember
!
  lidx = 1
  afalse = .false.
  doslice = .false.
  if (allocated(emcoarsegrid%lslc).EQV..true.) then
    doslice = .true.
  end if
  doline = .false.
  if (allocated(emcoarsegrid%llin).EQV..true.) then
    doline = .true.
  end if
  doblocks = .false.
  if (allocated(emcoarsegrid%lblk).EQV..true.) then
    doblocks = .true.
  end if
!
! independent of curmassm
  call em_pop_mv(emgrid,afalse)

  unitvol = emgrid%deltas(1)*emgrid%deltas(2)*emgrid%deltas(3)
  soldens = embgdensity/convdens ! in a.m.u per A^3
!
  if (empotmode.eq.1) then
    normer = convdens/unitvol
    emgrid%dens = normer*emgrid%mass
  else if (empotmode.eq.2) then
    normer = convdens/unitvol
    if (emgrid%cnt.le.0) then
      emgrid%dens = normer*emgrid%mass
    else
      normer = convdens/unitvol
      normer2 = (1.0-emiw)/(1.0*emgrid%cnt)
      emgrid%dens = normer*(normer2*emgrid%avgmass+emiw*emgrid%mass)
    end if
  end if
!
  ints(:) = nint(emcoarsegrid%deltas(:)/(emgrid%deltas(:)))
  intsm1(:) = ints(:) - 1
  normer = 1.0/(ints(1)*ints(2)*ints(3))
  accum = 0.0
  emcoarsegrid%deltadens(:,:,:) = 0. ! with heuristic, this will remain incomplete for parts with non-bg density but no atoms
!
  if (doblocks.EQV..true.) then
    do xk=1,emcoarsegrid%dimc(1)
      do yk=1,emcoarsegrid%dimc(2)
        do zk=1,emcoarsegrid%dimc(3)
          if (emcoarsegrid%lblk(xk,yk,zk,lidx).EQV..false.) cycle
          do yj=emcoarsegrid%blklms(yk,1,2),emcoarsegrid%blklms(yk,2,2)
            yi = yj*ints(2) - intsm1(2)
            do yl=0,ints(2)-1
              do zj=emcoarsegrid%blklms(zk,1,3),emcoarsegrid%blklms(zk,2,3)
                zi = zj*ints(3) - intsm1(3)
                do zl=0,ints(3)-1
                  do xj=emcoarsegrid%blklms(xk,1,1),emcoarsegrid%blklms(xk,2,1)
                    xi = xj*ints(1) - intsm1(1)
                    do xl=0,ints(1)-1
                      emcoarsegrid%deltadens(xj,yj,zj) = emcoarsegrid%deltadens(xj,yj,zj) + emgrid%dens(xi+xl,yi+yl,zi+zl)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    do xk=1,emcoarsegrid%dimc(1)
      do yk=1,emcoarsegrid%dimc(2)
        do zk=1,emcoarsegrid%dimc(3)
          if (emcoarsegrid%lblk(xk,yk,zk,lidx).EQV..false.) then
            accum = accum + emcoarsegrid%rblk(xk,yk,zk)
            cycle
          end if
          do yj=emcoarsegrid%blklms(yk,1,2),emcoarsegrid%blklms(yk,2,2)
            do zj=emcoarsegrid%blklms(zk,1,3),emcoarsegrid%blklms(zk,2,3)
              emcoarsegrid%deltadens(emcoarsegrid%blklms(xk,1,1):emcoarsegrid%blklms(xk,2,1),yj,zj) = &
 &                 normer*emcoarsegrid%deltadens(emcoarsegrid%blklms(xk,1,1):emcoarsegrid%blklms(xk,2,1),yj,zj) - &
 &                 emcoarsegrid%dens(emcoarsegrid%blklms(xk,1,1):emcoarsegrid%blklms(xk,2,1),yj,zj)
              accum = accum + sum(emcoarsegrid%deltadens(emcoarsegrid%blklms(xk,1,1):emcoarsegrid%blklms(xk,2,1),yj,zj)*&
 &                                emcoarsegrid%deltadens(emcoarsegrid%blklms(xk,1,1):emcoarsegrid%blklms(xk,2,1),yj,zj))
            end do
          end do
        end do
      end do
    end do
  else
    do yj=1,emcoarsegrid%dim(2)
      if (doslice.EQV..true.) then
        if (emcoarsegrid%lslc(yj,lidx).EQV..false.) cycle
      end if
      yi = yj*ints(2) - intsm1(2)
      do yk=0,ints(2)-1
        do zj=1,emcoarsegrid%dim(3)
          if (doline.EQV..true.) then
            if (emcoarsegrid%llin(yj,zj,lidx).EQV..false.) cycle
          end if
          zi = zj*ints(3) - intsm1(3)
          do zk=0,ints(3)-1
            do xj=1,emcoarsegrid%dim(1)
              xi = xj*ints(1) - intsm1(1)
              do xk=0,ints(1)-1
                emcoarsegrid%deltadens(xj,yj,zj) = emcoarsegrid%deltadens(xj,yj,zj) + emgrid%dens(xi+xk,yi+yk,zi+zk)
              end do
            end do
          end do
        end do
      end do
    end do
    do yj=1,emcoarsegrid%dim(2)
      if (doslice.EQV..true.) then
        if (emcoarsegrid%lslc(yj,lidx).EQV..false.) then
          accum = accum + emcoarsegrid%rslc(yj)
          cycle
        end if
      end if
      do zj=1,emcoarsegrid%dim(3)
        if (doline.EQV..true.) then
          if (emcoarsegrid%llin(yj,zj,lidx).EQV..false.) then
            accum = accum + emcoarsegrid%rlin(yj,zj)
            cycle
          end if
        end if
        emcoarsegrid%deltadens(:,yj,zj) = normer*emcoarsegrid%deltadens(:,yj,zj) - emcoarsegrid%dens(:,yj,zj)
        accum = accum + sum(emcoarsegrid%deltadens(:,yj,zj)*emcoarsegrid%deltadens(:,yj,zj))
      end do
    end do
  end if
!
  evec(21) = evec(21) + scale_EMICRO*accum
!
  unitvol = emgrid%deltas(1)*emgrid%deltas(2)*emgrid%deltas(3)
  if (empotmode.eq.1) then
    normer = scale_EMICRO*2.0*convdens/(ints(1)*ints(2)*ints(3)*unitvol)
  else if (empotmode.eq.2) then
    normer = scale_EMICRO*2.0*emiw*convdens/(ints(1)*ints(2)*ints(3)*unitvol)
  end if
  if (mod(emsplor,2).eq.0) then
    sh = emsplor/2 
  else
    sh = floor(emsplor/2.0)
  end if
  do imol=1,nmol
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
    do ii=atmol(imol,1),atmol(imol,2)
      if (empotprop.eq.1) then
        term0 = normer*(mass(ii) - soldens*atsavred(ii)*atvol(ii))
      else if (empotprop.eq.2) then
        term0 = normer*(1.0*lj_atnum(attyp(ii)) - soldens*atsavred(ii)*atvol(ii))
      end if
!     get the d/dxi (Mn(xi)): note that the grid axes are aligned with the coordinate axes,
!     so only diagonal terms contribute, i.e., d/dyi (Mn(xi)) is obviously zero
      do cc=1,3
        embsplbu(1:emsplor,1) = embspl(1:emsplor,cc,ii)
        embspl(1:emsplor,cc,ii) = embspld(1:emsplor,cc,ii)
        do j=1,emsplor
          temj = term0*embspl(j,1,ii)
          xi = emgfls(ii,1) + j - sh
          if (xi.gt.emgrid%dim(1)) xi = xi - emgrid%dim(1)
          if (xi.lt.1) xi = xi + emgrid%dim(1)
          xj = ceiling((xi*1.0-0.5)/(1.0*ints(1)))
          do k=1,emsplor
            yi = emgfls(ii,2) + k - sh
            if (yi.gt.emgrid%dim(2)) yi = yi - emgrid%dim(2)
            if (yi.lt.1) yi = yi + emgrid%dim(2)
            temk = temj*embspl(k,2,ii)
            yj = ceiling((yi*1.0-0.5)/(1.0*ints(2)))
            do l=1,emsplor
              zi = emgfls(ii,3) + l - sh
              if (zi.gt.emgrid%dim(3)) zi = zi - emgrid%dim(3)
              if (zi.lt.1) zi = zi + emgrid%dim(3)
              zj = ceiling((zi*1.0-0.5)/(1.0*ints(3)))
              teml = temk*embspl(l,3,ii)
              ca_f(ii,cc) = ca_f(ii,cc) - teml*emcoarsegrid%deltadens(xj,yj,zj)
            end do
          end do
        end do
        embspl(1:emsplor,cc,ii) = embsplbu(1:emsplor,1)
      end do
    end do
  end do
!
end
!
!-----------------------------------------------------------------------------
!
! routine for the polymer property-based global potential energy terms
! we are using an analytical eigenvalue formula derivable from transforming the eigenvalue problem
! into a root problem of a cubic polynomial and using Cardano's formula
! for a well-behaved matrix like the gyration tensor, this allows safe and accurate computation
! of derivatives 
!
subroutine force_poly_gl(imol,evec,ca_f)
!
  use atoms
  use molecule
  use energies
  use polyavg
!
  implicit none
!
  integer i,imol
  RTYPE t1,t2,t3,t4,t5,t6,t7,t8,h1,sh1,h2,st4,st6,pos(3),posr(3),nfac,sterm,cterm
  RTYPE drgtsdr(6,3),dt1dr(3),dt2dr(3),dt3dr(3),dt4dr(3),dt5dr(3)
  RTYPE dt6dr(3),dt7dr(3),dst4dr(3),drgevsdr(3,3),ddltsdr(3),dtpdr(3)
  RTYPE evec(MAXENERGYTERMS),dlts,tps,texp,tprefac,rgten(3,3),ca_f(n,3),dltsd,tpsd
!
  if (par_POLY2(imol).le.0) return
!
  if (atmol(imol,2).le.atmol(imol,1)) return
!
  nfac = 1.0/(dble(atmol(imol,2)-atmol(imol,1)+1))
  h1 = (1./3.)
  sh1 = sqrt(h1)
  com(imol,1) = nfac*sum(x(atmol(imol,1):atmol(imol,2)))
  com(imol,2) = nfac*sum(y(atmol(imol,1):atmol(imol,2)))
  com(imol,3) = nfac*sum(z(atmol(imol,1):atmol(imol,2)))
  rgten(:,:) = 0.0
  do i=atmol(imol,1),atmol(imol,2)
    rgten(1,1) = rgten(1,1) + (x(i)-com(imol,1))**2
    rgten(1,2) = rgten(1,2) + (x(i)-com(imol,1))*(y(i)-com(imol,2))
    rgten(1,3) = rgten(1,3) + (x(i)-com(imol,1))*(z(i)-com(imol,3))
    rgten(2,2) = rgten(2,2) + (y(i)-com(imol,2))**2
    rgten(2,3) = rgten(2,3) + (y(i)-com(imol,2))*(z(i)-com(imol,3))
    rgten(3,3) = rgten(3,3) + (z(i)-com(imol,3))**2
  end do
  rgten(:,:) = nfac*rgten(:,:)
  t1 = rgten(1,1)+rgten(2,2)+rgten(3,3)
  t2 = rgten(1,1)*rgten(2,2) + rgten(1,1)*rgten(3,3) + rgten(2,2)*rgten(3,3) - rgten(1,2)**2 - rgten(1,3)**2 - rgten(2,3)**2
  t3 = rgten(1,1)*rgten(2,3)**2 + rgten(2,2)*rgten(1,3)**2 + rgten(3,3)*rgten(1,2)**2 - rgten(1,1)*rgten(2,2)*rgten(3,3) - &
 &     2.0*rgten(1,2)*rgten(1,3)*rgten(2,3)
  t4 = t1**2 - 3.0*t2
  t5 = t1*(t4 - 1.5*t2) - 13.5*t3
  t6 = 27.0*( 0.25*(t2**2)*(t4 - t2) + t3*(t5 + 6.75*t3) )
  st6 = 0.0
  if (t6.ge.0.0) st6 = sqrt(t6)
  t7 = h1*atan2(st6,t5)
  sterm = sin(t7)
  cterm = cos(t7)
  st4 = 0.0
  if (t4.ge.0.0) st4 = sqrt(t4)
  h2 = sh1*st4*sin(t7)
  rgevs(imol,1) = h1*(t1-st4*cos(t7))
  rgevs(imol,2) = rgevs(imol,1) - h2
  rgevs(imol,3) = rgevs(imol,1) + h2
  rgevs(imol,1) = rgevs(imol,1) + st4*cos(t7)
  t8 = rgevs(imol,1)*rgevs(imol,2) + rgevs(imol,1)*rgevs(imol,3) + rgevs(imol,2)*rgevs(imol,3)
! order does not matter for delta* or t --> order later
  dlts = 1.0 - 3.0*t8/(t1**2)
  dltsd = scale_POLY*2.0*par_POLY(imol,4)*(dlts - par_POLY(imol,2))
  texp = (4.0/(dble(rsmol(imol,2)-rsmol(imol,1)+1))**tp_exp)
  tprefac = 2.5*((1.75/molcontlen(moltypid(imol)))**texp)
  tps = tprefac*rgv(imol)**texp
  tpsd = scale_POLY*2.0*par_POLY(imol,3)*(tps - par_POLY(imol,1))
  texp = 0.5*texp
  evec(14) = evec(14) + 0.5*( tpsd*(tps - par_POLY(imol,1)) + dltsd*(dlts - par_POLY(imol,2)) )
!
  nfac = 1.0/(dble(atmol(imol,2)-atmol(imol,1)+1))
  do i=atmol(imol,1),atmol(imol,2)
    pos(1) = x(i)
    pos(2) = y(i)
    pos(3) = z(i)
    posr(:) = pos(:) - com(imol,:)
    drgtsdr(:,:) = 0.0
!   diagonal elements
    drgtsdr(1,1) = 2.0*posr(1)*nfac
    drgtsdr(5,2) = 2.0*posr(2)*nfac
    drgtsdr(6,3) = 2.0*posr(3)*nfac
!   offdiagonal elements
    drgtsdr(2,1) = posr(2)*nfac ! xy
    drgtsdr(2,2) = posr(1)*nfac
    drgtsdr(3,1) = posr(3)*nfac ! xz
    drgtsdr(3,3) = posr(1)*nfac
    drgtsdr(4,2) = posr(3)*nfac ! yz
    drgtsdr(4,3) = posr(2)*nfac
!   t1 = Rg^2
    dt1dr(:) = drgtsdr(1,:)+drgtsdr(5,:)+drgtsdr(6,:)
!   t2 = Rxx*Ryy + Ryy*Rzz + Rxx*Rzz - Rxy^2 - Rxz^2 - Ryz^2
    dt2dr(:) =  drgtsdr(1,:)*rgten(2,2) + drgtsdr(5,:)*rgten(1,1) + &
 &              drgtsdr(5,:)*rgten(3,3) + drgtsdr(6,:)*rgten(2,2) + &
 &              drgtsdr(1,:)*rgten(3,3) + drgtsdr(6,:)*rgten(1,1) - &
 &              2.0*(drgtsdr(2,:)*rgten(1,2) + drgtsdr(3,:)*rgten(1,3) + drgtsdr(4,:)*rgten(2,3))
!   t3 = Rxx*Ryz^2 + Ryy*Rxz^2 + Rzz*Rxy^2 - Rxx*Ryy*Rzz  - 2*Rxy*Rxz*Ryz
    dt3dr(:) = drgtsdr(1,:)*rgten(2,3)**2 + 2.0*rgten(1,1)*rgten(2,3)*drgtsdr(4,:) + &
 &             drgtsdr(5,:)*rgten(1,3)**2 + 2.0*rgten(2,2)*rgten(1,3)*drgtsdr(3,:) + &
 &             drgtsdr(6,:)*rgten(1,2)**2 + 2.0*rgten(3,3)*rgten(1,2)*drgtsdr(2,:) - &
 &             rgten(2,2)*rgten(3,3)*drgtsdr(1,:) - rgten(1,1)*rgten(3,3)*drgtsdr(5,:) - rgten(1,1)*rgten(2,2)*drgtsdr(6,:) - &
 &        2.0*(rgten(1,3)*rgten(2,3)*drgtsdr(2,:) + rgten(1,2)*rgten(2,3)*drgtsdr(3,:) + rgten(1,2)*rgten(1,3)*drgtsdr(4,:) )
!   t4 = t1^2 - 3*t2
    dt4dr(:) = 2.0*t1*dt1dr(:) - 3.0*dt2dr(:)
!   t5 = t1*(t4 - 1.5*t2) - 13.5*t3
    dt5dr(:) = dt1dr(:)*(t4 - 1.5*t2) + t1*(dt4dr(:) - 1.5*dt2dr(:)) - 13.5*dt3dr(:)
!   t6 = 27*((1/4)*(t2**2)*(t4 - t2) + t3*(t5 + (27/4)*t3) )
    dt6dr(:) = 27.0*( dt3dr(:)*(t5 + 6.75*t3) + t3*(dt5dr(:) + 6.75*dt3dr(:)) + &
 &                    0.5*t2*dt2dr(:)*(t4-t2) + 0.25*(t2**2)*(dt4dr(:)-dt2dr(:)) )
!   t7 = (1/3)*atan2(sqrt(t6),t5)
    dt7dr(:) = 0.0
    if (t6.gt.0.0) dt7dr(:) = h1*(t5*dt6dr(:)/(2.0*st6) - st6*dt5dr(:))/(t6 + t5**2.0)
    dst4dr(:) = 0.0
    if (t4.gt.0.0) dst4dr(:) = dt4dr(:)/(2.0*st4)
    drgevsdr(1,:) = h1*(dt1dr(:)-dst4dr(:)*cterm + st4*sterm*dt7dr(:))
    drgevsdr(2,:) = drgevsdr(1,:) - sh1*st4*cterm*dt7dr(:) - sh1*dst4dr(:)*sterm
    drgevsdr(3,:) = drgevsdr(1,:) + sh1*st4*cterm*dt7dr(:) + sh1*dst4dr(:)*sterm
    drgevsdr(1,:) = drgevsdr(1,:) + dst4dr(:)*cterm - st4*sterm*dt7dr(:)
    ddltsdr(:) = (-3.0/t1**2)*( (rgevs(imol,1)*drgevsdr(2,:) + rgevs(imol,2)*drgevsdr(1,:) + rgevs(imol,1)*drgevsdr(3,:) + &
 &                               rgevs(imol,3)*drgevsdr(1,:) + rgevs(imol,2)*drgevsdr(3,:) + rgevs(imol,3)*drgevsdr(2,:)) - &
 &                               (2.0/t1)*t8*dt1dr(:) )
    dtpdr(:) = tprefac*texp*dt1dr(:)*(t1**(texp-1.0))
    ca_f(i,:) = ca_f(i,:) - ( tpsd*dtpdr(:) + dltsd*ddltsdr(:) )
  end do
!
! sort eigenvalues
  posr(:) = rgevs(imol,:)
  if (h2.gt.0.0) then
    rgevs(imol,1) = min(posr(1),posr(2))
    rgevs(imol,3) = max(posr(3),posr(1))
    rgevs(imol,2) = max(posr(2),min(posr(3),posr(1)))
  else
    rgevs(imol,1) = min(posr(1),posr(3))
    rgevs(imol,3) = max(posr(2),posr(1))
    rgevs(imol,2) = max(posr(3),min(posr(2),posr(1)))
  end if
  if ((rgevs(imol,1).gt.rgevs(imol,2)).OR.(rgevs(imol,2).gt.rgevs(imol,3)).OR.(rgevs(imol,1).gt.rgevs(imol,3))) call fexit()
!
end
!
!-----------------------------------------------------------------------------------------------
!
! a helper routine to get the outer derivative of the interpolation
! function
!
function fsdr_ipol(ati)
!
  use atoms
  use iounit
  use energies
!
  implicit none
!
  RTYPE fsdr_ipol,rat,expte
  integer ati
!
  rat = atsav(ati)/atbvol(ati)
  if (rat.gt.atsavmaxfr(ati)) then
    fsdr_ipol = 0.0
  else if (rat.gt.par_IMPSOLV(5)) then
    expte = exp(-(rat-atsavprm(ati,1))/par_IMPSOLV(3))
    fsdr_ipol = (atsavprm(ati,2)/par_IMPSOLV(3))*expte /&
 &        ((1.0 + expte)*(1.0 + expte))
  else
    fsdr_ipol = 0.0
  end if
!  if((ati.eq.3).OR.(ati.eq.11).OR.(ati.eq.10).OR.(ati.eq.17)) then
!    write(ilog,*) ati,' ',fs_ipol
!  end if
!
  return
!
end
!
!---------------------------------------------------------------------------
!
! this gets a bit messy (the function doesn't really show that, though),
! because of the coupling of so many degrees of freedom through atomic solvation states
!
subroutine force_freesolv(rs,evec,ca_f)
!
  use energies
  use polypep
  use iounit
  use sequen
  use atoms
  use forces
!
  implicit none
!
  integer i,rs,j,k,kk,ii,ssi
  RTYPE evec(MAXENERGYTERMS),ssav,fs_ipol,fsdr_ipol,fac,fos_Tdep
  RTYPE ca_f(n,3)
!
  do ssi=1,at(rs)%nfosgrps
!
    ssav = 0.0
!
    do i=1,at(rs)%fosgrp(ssi)%nats
!
      ii = at(rs)%fosgrp(ssi)%ats(i)
!
!     handle energy
      ssav = ssav + at(rs)%fosgrp(ssi)%wts(i)*fs_ipol(ii)
!
!     get the pre-factor (note that to normalize sav_dr/sisa_dr, atbvol
!     for atom ii is used)
      fac = fos_Tdep(at(rs)%fosgrp(ssi)%val(1:3))*at(rs)%fosgrp(ssi)%wts(i)*&
 &                      fsdr_ipol(ii)*scale_IMPSOLV/atbvol(ii)
!
!     handle dFOS_ssi/dxii
      do j=1,3
        ca_f(ii,j) = ca_f(ii,j) + fac*sav_dr(ii,j)
      end do
 
!     now loop over all the neighbors of atom ii (i.e., those terms, for which dSAVii/dxkk 
!     was not zero, and increment the respective cartesian derivatives
      do k=1,sisa_nr(ii)
        kk = sisa_i(k,ii)
        do j=1,3
          ca_f(kk,j) = ca_f(kk,j) + fac*sisa_dr(k,j,ii)
        end do
      end do
!
    end do
!    write(*,*) rs,ssav
!
!   finalize energy term
    evec(4) = evec(4) + scale_IMPSOLV*ssav*fos_Tdep(at(rs)%fosgrp(ssi)%val(1:3))
!
  end do
!
end
!
!----------------------------------------------------------------------------
!
!  THE ROUTINE TO SETUP THE DIPOLE-GROUP BASED SCREENED ATOMIC CHARGES 
!
!----------------------------------------------------------------------------
!
!
! a helper routine to get the outer derivative of the interpolation
! function
!
function sqdr_ipol(ati)
!
  use atoms
  use iounit
  use energies
!
  implicit none
!
  integer ati
  RTYPE sqdr_ipol,rat,expte
!
  rat = atsav(ati)/atbvol(ati)
!
  if (rat.ge.atsavmaxfr(ati)) then
    sqdr_ipol = 0.0
  else if (rat.gt.par_IMPSOLV(5)) then
    expte = exp(-(rat-atsavprm(ati,4))/par_IMPSOLV(4))

    sqdr_ipol = (atsavprm(ati,5)/par_IMPSOLV(4))*expte /&
 &        ((1.0 + expte)*(1.0 + expte))
  else
    sqdr_ipol = 0.0
  end if
!
  return
!
end
!
!-------------------------------------------------------------------------------
!
! now here it gets really nasty:
! we do in fact need to know the derivatives dscrqii/dxkk with arbitrary
! ii (any polar atom) and arbitrary kk (any other polar atom), since we'll
! have to take the derivative of coulomb's law later, which might have any
! type of usage of the scrqii's
! however, for every term in the Coulomb sum (ii,kk), we'll also need the non-zero
! derivatives of scrqii and scrqkk with respect to any arbitrary atom jj, since they'll
! clearly yield an (albeit small) term to -dU/dxjj
! the strategy is that we accumulate the dscrqii/dxkk similarly to the dsavii/dxkk
! in a 3D-array (sisq_dr), which is particularly complicated for screening models
! which are charge group-consistent 
!
subroutine force_setup_scrqs2()
!
  use iounit
  use aminos
  use sequen
  use polypep
  use energies
  use molecule
  use atoms
  use forces
!
  implicit none
!
  integer i,j,ii,rs,ll,k,kk,idx,m,mm,dpi,l,revidx(n)
  RTYPE sq_ipol,sqdr_ipol,fac,qss
  logical gotii
!
! in model #4 (pure distance-dependent dielectric), we don't need
! environmentally screened charges
  if (scrq_model.eq.4) then
!
!   do nothing
!
! in these models the screened charges are purely atom-based 
  else if ((scrq_model.eq.2).OR.(scrq_model.eq.6).OR.&
 &    (scrq_model.eq.8).OR.(scrq_model.eq.9)) then
!
    do rs=1,nseq
!
      do i=1,at(rs)%npol
!
        ii = at(rs)%pol(i)
!
!       handle actual screened charge
        scrq(ii) = (1.0 - coul_scr*sq_ipol(ii))
!
!       pre-factor
        fac = -coul_scr*sqdr_ipol(ii)/atbvol(ii)
!
!       handle dscrq_ii/dxii
!          write(*,*) 'now incrementing ',ii,fac*sav_dr(ii,3)
        do j=1,3
          scrq_dr(ii,j) = scrq_dr(ii,j) + fac*sav_dr(ii,j)
        end do
 
!       now loop over all the neighbors of atom ii (i.e., those terms, for which dSAVii/dxkk 
!       was not zero, and generate the corresponding derivatives dscrqii/dxkk
        do k=1,sisa_nr(ii)
          kk = sisa_i(k,ii)
!            write(*,*) 'now incr. ',kk,fac*sisa_dr(k,3,ii)
          do j=1,3
            sisq_dr(k,j,ii) = fac*sisa_dr(k,j,ii)
          end do
        end do
!
      end do
    end do
!
  else if ((scrq_model.eq.3).OR.(scrq_model.eq.5).OR.&
 &    (scrq_model.eq.1).OR.(scrq_model.eq.7)) then
!
    do rs=1,nseq
!
      do dpi=1,at(rs)%ndpgrps
!
        qss = 0.0
!
        do i=1,at(rs)%dpgrp(dpi)%nats
!
          ii = at(rs)%dpgrp(dpi)%ats(i)
!
!         handle actual screened charge
          qss = qss + abs(atq(ii))*sq_ipol(ii)
!
!         get the pre-factor (note that to normalize sav_dr/sisa_dr, atbvol
!         for atom ii is used)
          fac =-coul_scr*abs(atq(ii))*sqdr_ipol(ii)/&
 &                      (atbvol(ii)*at(rs)%dpgrp(dpi)%qnm)
!
!         handle dscrq_ii/dxii
!          write(*,*) 'now incrementing ',ii,fac*sav_dr(ii,3)
          do j=1,3
            scrq_dr(ii,j) = scrq_dr(ii,j) + fac*sav_dr(ii,j)
          end do
!
!         now handle dscrq_ll/dxii, for all ll which belong to the same charge group
          idx = which_dpg(ii)
          do l=1,at(rs)%dpgrp(idx)%nats
            ll = at(rs)%dpgrp(idx)%ats(l)
            if (ll.eq.ii) cycle
            gotii = .false.
            do m=1,sisa_nr(ll)
              mm = sisa_i(m,ll)
!             if ii is also already a neighbor of ll, we'll just use that position (m)
!             in the array for ll
              if (mm.eq.ii) then
                do j=1,3
                  sisq_dr(m,j,ll) = &
 &                      sisq_dr(m,j,ll)+fac*sav_dr(ii,j)
                end do
                gotii = .true.
                exit
              end if
            end do
!           if ii wasn't originally part of the neighbors of ll, we'll have to add it
!           to not screw up, we'll explicitly set sisa_dr for this to 0.0 (even though it
!           shouldn't be used hereafter)
            if (gotii.EQV..false.) then
              sisa_nr(ll) = sisa_nr(ll) + 1
              sisa_i(sisa_nr(ll),ll) = ii
              do j=1,3
                sisa_dr(sisa_nr(ll),j,ll) = 0.0
                sisq_dr(sisa_nr(ll),j,ll) = fac*sav_dr(ii,j)
              end do
            end if
          end do
!         now loop over all the neighbors of atom ii (i.e., those terms, for which dSAVii/dxkk 
!         was not zero), increment dscrqii/dxkk, and search for that same neighbor in the other
!         atoms of the charge group in order to increment the same term dscrqii/dxkk for those.
!         effectively, rather than xkk affecting a single screening factor strongly, in the group-
!         consistent modes, the derivative dscrqii/dxkk is splayed out over multiple charge
!         sites (which might be neither ii or kk!, but have to part of the same charge group
!         as atom ii)
          do j=1,3
            sisq_dr(1:sisa_nr(ii),j,ii) = sisq_dr(1:sisa_nr(ii),j,ii)+fac*sisa_dr(1:sisa_nr(ii),j,ii)
          end do
          revidx(:) = 0
          do l=1,at(rs)%dpgrp(idx)%nats
            ll = at(rs)%dpgrp(idx)%ats(l)
            if (ll.eq.ii) cycle
!           assemble the temporary reverse-indexed list
            do k=1,sisa_nr(ii)
              revidx(sisa_i(k,ii)) = k
            end do
            do m=1,sisa_nr(ll)
              mm = sisa_i(m,ll)
              if (revidx(mm).gt.0) then
!               if mm is also in the list of neighbors for ii (via revidx), we'll just use its position (m)
!               in the array for ll
                do j=1,3
                  sisq_dr(m,j,ll) = sisq_dr(m,j,ll)+fac*sisa_dr(revidx(mm),j,ii)
                end do
                revidx(mm) = 0
              end if
            end do
            do k=1,sisa_nr(ii)
              kk = sisa_i(k,ii)
              if (revidx(kk).gt.0) then
!               if kk wasn't originally part of the neighbors of ll, we'll have to add it
!               to not screw up, we'll explicitly set sisa_dr for this to 0.0 (even though it
!               shouldn't be used hereafter)
                sisa_nr(ll) = sisa_nr(ll) + 1
                sisa_i(sisa_nr(ll),ll) = kk
!                 write(*,*) 'now incr.2 ',ll,kk,fac*sisa_dr(k,3,ii)
                do j=1,3
                  sisa_dr(sisa_nr(ll),j,ll) = 0.0
                  sisq_dr(sisa_nr(ll),j,ll) = fac*sisa_dr(k,j,ii)
                end do
                revidx(kk) = 0
              end if
            end do
          end do
!
        end do
!
!       finalize actual screened charges
        qss = qss/at(rs)%dpgrp(dpi)%qnm
        do i=1,at(rs)%dpgrp(dpi)%nats
          scrq(at(rs)%dpgrp(dpi)%ats(i)) = (1.0-coul_scr*qss)
        end do
!
      end do
    end do
!
  end if
!
end
!
!---------------------------------------------------------------------------
!
! this helper routine's only job is to find the right derivative
! it is currently not in use
!
subroutine force_setup_scrqs3(ii,kk,which,sdr)
!
  use forces
!
  implicit none
!
  integer ii,kk,i,k,which,j
  RTYPE sdr(2,3)
!
! dscrqkk/dxii
  if (which.eq.1) then
    do k=1,sisa_nr(kk)
      if (ii.eq.sisa_i(k,kk)) then
        do j=1,3
          sdr(1,j) = sisq_dr(k,j,kk)
        end do
        return
      end if
    end do
    do j=1,3
      sdr(1,j) = 0.0
    end do
! dscrqii/dxkk
  else if (which.eq.2) then
    do i=1,sisa_nr(ii)
      if (kk.eq.sisa_i(i,ii)) then
        do j=1,3
          sdr(2,j) = sisq_dr(i,j,ii)
        end do
        return
      end if
    end do
    do j=1,3
      sdr(2,j) = 0.0
    end do
  end if
!
end
!
!---------------------------------------------------------------------------
!
! NOTE: this function is NOT obsolete unlike the other force_rsp_ ones 
!
! this is a specialized subroutine supposed to deal with the modified
! Coulomb term for Ewald summation / (Generalized) Reaction-Field
! currently this function is NOT supposed to work with the implicit solvent model
! support for some ghosting (FEG) calculations is only provided for (G)RF
!
subroutine force_rsp_long_mod(evec,rs1,rs2,cut,ca_f)
!
  use iounit
  use cutoffs
  use units
  use energies
  use atoms
  use polypep
  use inter
  use ewalds
  use forces
  use sequen
!
  implicit none
!
  RTYPE evec(MAXENERGYTERMS),term1,term0,terme,term2
  RTYPE dis,dis2,svec(3),dvec(3)
  RTYPE ca_f(n,3)
  integer rs,rs1,rs2,i,j,k,ii,kk
  logical cut
!
! first determine whether we need to jump out for FEG
  if (use_FEG.EQV..true.) then
    if ((fegcbmode.ne.1).OR.(fegmode.ne.2)) then
      write(ilog,*) 'Fatal. Wrong mode settings for Coulomb-scaling &
 &and/or FEG-mode itself for RF support.'
      call fexit()
    end if
    if (fegmode.eq.2) then
!     some sanity checks
!      if ((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..true.)&
! &.AND.(rs1.ne.rs2)) then
!        write(ilog,*) 'Fatal. Simple Cb-scaling in FEG is incompatib&
! &le with multiple ghosted residues. This is a bug.'
!        call fexit()
!      end if
      if((par_FEG(rs1).EQV..true.).OR.(par_FEG(rs2).EQV..true.))then
        call force_rsp_long_feg_mod(evec,rs1,rs2,cut,ca_f)
        return
      end if
    end if
  end if
!
! no intialization needed: note that on the calling side it requires a lot of
! care to handle fxns which "silently" increment arguments
!
! three different cases: intra-residue (rs1 == rs2), neighbors in seq., or others
!
  if (rs1.eq.rs2) then
!
    rs = rs1
!
!   loop over relevant intra-residue interactions, watch cutoff for sanity
    do i=1,nrpolintra(rs)
      ii = iaa(rs)%polin(i,1)
      kk = iaa(rs)%polin(i,2)
      dvec(1) = x(kk) - x(ii)
      dvec(2) = y(kk) - y(ii)
      dvec(3) = z(kk) - z(ii)
      dis2 = dvec(1)*dvec(1) + dvec(2)*dvec(2) + dvec(3)*dvec(3)
      dis = sqrt(dis2)
!
!     things are nice and simple in vacuo ...
      term0 = scale_POLAR*fudge(rs)%elin(i)*electric*&
 &                                 atq(ii)*atq(kk)
      term1 = term0/dis
      if (lrel_md.eq.2) then
        terme = term1*erfc(ewpm*dis)
        evec(6) = evec(6) + terme
        term2 = (1.0/dis2)*(terme + term1*ewpite*dis*exp(-ewpm2*dis2))
        do j=1,3
          ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
          ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
        end do
!     note that the RF energy correction is meant to be compatible with the
!     force correction (which is really the proper term)
      else if (lrel_md.eq.3) then
        terme = term1+term0*(par_POLAR(1)*dis2-par_POLAR(2))
        evec(6) = evec(6) + terme
        term2 = (1.0/dis2)*(term1 - 2.0*term0*dis2*par_POLAR(1))
        do j=1,3
          ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
          ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
        end do
      else
        call fexit()
      end if
!
!     extra correction for fudged terms (shows yet again how broken and undesirable fudges are)
      if (fudge(rs)%elin(i).lt.1.0) then
        term0 = scale_POLAR*electric*atq(ii)*atq(kk)*(1.0-fudge(rs)%elin(i))
        term1 = term0/dis
!
        if (lrel_md.eq.2) then
          terme = -term1*erf(ewpm*dis)
          evec(6) = evec(6) + terme
          term2 = (1.0/dis2)*(terme + term1*ewpite*dis*exp(-ewpm2*dis2))
          do j=1,3
            ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
            ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
          end do
        else if (lrel_md.eq.3) then
          terme = term0*(par_POLAR(1)*dis2-par_POLAR(2))
          evec(6) = evec(6) + terme
          term2 = -2.0*term0*par_POLAR(1)
          do j=1,3
            ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
            ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
          end do
        end if
      end if
    end do
!
!   for both Ewald and RF we have to go through the exclusion list and
!   explicitly add the corrections for those terms (some of them are constant
!   due to constraints, but it's easier, safer, and not much more expensive
!   to compute them instantaneously)
    do i=1,nrexpolin(rs)
      ii = iaa(rs)%expolin(i,1)
      kk = iaa(rs)%expolin(i,2)
      dvec(1) = x(kk) - x(ii)
      dvec(2) = y(kk) - y(ii)
      dvec(3) = z(kk) - z(ii)
      dis2 = dvec(1)*dvec(1) + dvec(2)*dvec(2) + dvec(3)*dvec(3)
      if ((cut.EQV..true.).AND.(dis2.gt.mcel_cutoff2)) cycle
      dis = sqrt(dis2)
!
      term0 = scale_POLAR*electric*atq(ii)*atq(kk)
      term1 = term0/dis
!
      if (lrel_md.eq.2) then
        terme = -term1*erf(ewpm*dis)
        evec(6) = evec(6) + terme
        term2 = (1.0/dis2)*(terme + term1*ewpite*dis*exp(-ewpm2*dis2))
        do j=1,3
          ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
          ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
        end do
      else if (lrel_md.eq.3) then
        terme = term0*(par_POLAR(1)*dis2-par_POLAR(2))
        evec(6) = evec(6) + terme
        term2 = -2.0*term0*par_POLAR(1)
        do j=1,3
          ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
          ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
        end do
      end if
    end do
!
! neighboring residues
!
  else if (abs(rs1-rs2).eq.1) then
!
!   there is no topological requirement for neighboring residues to be close (pure sequence)
!   so we have to check for BC 
!
    if (rs1.gt.rs2) then
      rs = rs2
      call dis_bound_rs(rs2,rs1,svec)
    else
      rs = rs1
      call dis_bound_rs(rs1,rs2,svec)
    end if
!
!   loop over relevant nearest-neighbor-residue interactions
    do i=1,nrpolnb(rs)
      ii = iaa(rs)%polnb(i,1)
      kk = iaa(rs)%polnb(i,2)
      call dis_bound5(ii,kk,svec,dis2,dvec)
      if ((cut.EQV..true.).AND.(dis2.gt.mcel_cutoff2)) cycle
      dis = sqrt(dis2)
!
!     things are nice and simple in vacuo ...
      term0 = scale_POLAR*fudge(rs)%elnb(i)*electric*&
 &                                 atq(ii)*atq(kk)
      term1 = term0/dis
      if (lrel_md.eq.2) then
        terme = term1*erfc(ewpm*dis)
        evec(6) = evec(6) + terme
        term2 = (1.0/dis2)*(terme + term1*ewpite*dis*exp(-ewpm2*dis2))
        do j=1,3
          ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
          ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
        end do
      else if (lrel_md.eq.3) then
        terme = term1+term0*(par_POLAR(1)*dis2-par_POLAR(2))
        evec(6) = evec(6) + terme
        term2 = (1.0/dis2)*(term1 - 2.0*term0*dis2*par_POLAR(1))
        do j=1,3
          ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
          ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
        end do
      else
        call fexit()
      end if
!
!     extra correction for fudged terms
      if (fudge(rs)%elnb(i).lt.1.0) then
        term0 = scale_POLAR*electric*atq(ii)*atq(kk)*(1.0-fudge(rs)%elnb(i))
        term1 = term0/dis
!
        if (lrel_md.eq.2) then
          terme = -term1*erf(ewpm*dis)
          evec(6) = evec(6) + terme
          term2 = (1.0/dis2)*(terme + term1*ewpite*dis*exp(-ewpm2*dis2))
          do j=1,3
            ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
            ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
          end do
        else if (lrel_md.eq.3) then
          terme = term0*(par_POLAR(1)*dis2-par_POLAR(2))
          evec(6) = evec(6) + terme
          term2 = -2.0*term0*par_POLAR(1)
          do j=1,3
            ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
            ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
          end do
        end if
      end if
    end do
!
!   for both Ewald and RF we have to go through the exclusion list and
!   explicitly add the corrections for those terms (some of them are constant
!   due to constraints, but it's easier, safer, and not much more expensive
!   to compute them instantaneously)
!   note that we do not use minimum image here, as excluded neighbor
!   interactions are ALWAYS on the same molecule and in that case dis_bound_rs
!   ALWAYS returns a zero shift vector (svec) 
    do i=1,nrexpolnb(rs)
      ii = iaa(rs)%expolnb(i,1)
      kk = iaa(rs)%expolnb(i,2)
      call dis_bound5(ii,kk,svec,dis2,dvec)
      dis = sqrt(dis2)
!
      term0 = scale_POLAR*electric*atq(ii)*atq(kk)
      term1 = term0/dis
      if (lrel_md.eq.2) then
        terme = -term1*erf(ewpm*dis)
        evec(6) = evec(6) + terme
        term2 = (1.0/dis2)*(terme + term1*ewpite*dis*exp(-ewpm2*dis2))
        do j=1,3
          ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
          ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
        end do
      else if (lrel_md.eq.3) then
        terme = term0*(par_POLAR(1)*dis2-par_POLAR(2))
        evec(6) = evec(6) + terme
        term2 = -2.0*term0*par_POLAR(1)
        do j=1,3
          ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
          ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
        end do
      end if
    end do
!
! all other residue pairs
!
  else
!
!   there is no topological relationship for remaining residues -> always check BC
!   since we have no exclusions, no special computation necessary
    call dis_bound_rs(rs1,rs2,svec)
!
    do i=1,at(rs1)%npol
!        write(*,*) sisa_nr(at(rs1)%pol(i))
      do k=1,at(rs2)%npol
        ii = at(rs1)%pol(i)
        kk = at(rs2)%pol(k)
        call dis_bound5(ii,kk,svec,dis2,dvec)
        if ((cut.EQV..true.).AND.(dis2.gt.mcel_cutoff2)) cycle
        dis = sqrt(dis2)
!
!       things are nice and simple in vacuo ...
        term0 = scale_POLAR*electric*atq(ii)*atq(kk)
        term1 = term0/dis
        if (lrel_md.eq.2) then
          terme = term1*erfc(ewpm*dis)
          evec(6) = evec(6) + terme
          term2 = (1.0/dis2)*(terme + term1*ewpite*dis*exp(-ewpm2*dis2))
          do j=1,3
            ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
            ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
          end do
        else if (lrel_md.eq.3) then
          terme = term1+term0*(par_POLAR(1)*dis2-par_POLAR(2))
          evec(6) = evec(6) + terme
          term2 = (1.0/dis2)*(term1 - 2.0*term0*dis2*par_POLAR(1))
          do j=1,3
            ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
            ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
          end do
        else
          call fexit()
        end if
      end do
    end do
!
  end if
!
end
!
!---------------------------------------------------------------------------
!
! NOTE: this function is NOT obsolete unlike the other force_rsp_ ones
!
! this variant exclusively supports RF electrostatics in conjunction with
! simple linear Cb-scaling for FEG
! the only way we can write a clean GRF treatment for scaled interactions, however,
! is to map the linear scaling factor to the charges which has consequences for
! ghost-ghost interactions which now need the scaling factor squared
! currently, this routine is mutually inconsistent with all other approaches
! toward ghost-ghost interactions(!!!)
!
subroutine force_rsp_long_feg_mod(evec,rs1,rs2,cut,ca_f)
!
  use iounit
  use energies
  use atoms
  use polypep
  use inter
  use units
  use tabpot
  use forces
  use cutoffs
!
  implicit none
!
  RTYPE evec(MAXENERGYTERMS),term1,dis2,terme
  RTYPE dis,svec(3),dvec(3),term0,term2,pfac2
  RTYPE ca_f(n,3)
  integer rs,rs1,rs2,i,j,k,ii,kk,aone,atwo
  logical cut
!
  aone = 1
  atwo = 2
!
! we have a potential override to cover in par_FEG3, which allows residues to be
! fully de-coupled at all times
  if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) return
!
!  pfac2 = par_FEG2(9)*par_FEG2(9)
  if ((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..true.)) then
    pfac2 = par_FEG2(9)*par_FEG2(9)
  else
    pfac2 = par_FEG2(9)
  end if
!
!
! three different cases: intra-residue (rs1 == rs2), neighbors in seq., or others
! note that we assume that this is the ghosted residue (!!!!)
!
  if (rs1.eq.rs2) then
!
    if (use_FEGS(6).EQV..true.) then
! 
      rs = rs1

!     loop over relevant intra-residue interactions, watch cutoff for sanity
      do i=1,nrpolintra(rs)
        ii = iaa(rs)%polin(i,1)
        kk = iaa(rs)%polin(i,2)
        dvec(1) = x(kk) - x(ii)
        dvec(2) = y(kk) - y(ii)
        dvec(3) = z(kk) - z(ii)
        dis2 = dvec(1)*dvec(1) + dvec(2)*dvec(2) + dvec(3)*dvec(3)
        if ((cut.EQV..true.).AND.(dis2.gt.mcel_cutoff2)) cycle
        dis = sqrt(dis2)
!
        term0 = pfac2*fudge(rs)%elin(i)*electric*atq(ii)*atq(kk)
        term1 = term0/dis
!       note that the RF energy correction is meant to be compatible with the
!       force correction (which is really the proper term)
        if (lrel_md.eq.3) then
          terme = term1+term0*(par_POLAR(1)*dis2-par_POLAR(2))
          evec(6) = evec(6) + terme
          term2 = (1.0/dis2)*(term1 - 2.0*term0*dis2*par_POLAR(1))
          do j=1,3
            ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
            ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
          end do
        else
          call fexit()
        end if
!
!       extra correction for fudged terms (shows yet again how broken and undesirable fudges are)
        if (fudge(rs)%elin(i).lt.1.0) then
          term0 = pfac2*electric*atq(ii)*atq(kk)*(1.0-fudge(rs)%elin(i))
          term1 = term0/dis
!
!          if (lrel_md.eq.2) then
!            terme = -term1*erf(ewpm*dis)
!            evec(6) = evec(6) + terme
!            term2 = (1.0/dis2)*(terme - &
! &                    term1*ewpite*dis*exp(-ewpm2*dis2))
!            do j=1,3
!              ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
!              ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
!            end do
          if (lrel_md.eq.3) then
            terme = term0*(par_POLAR(1)*dis2-par_POLAR(2))
            evec(6) = evec(6) + terme
            term2 = -2.0*term0*par_POLAR(1)
            do j=1,3
              ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
              ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
            end do
          end if
        end if
      end do
!
!     for RF we have to go through the exclusion list and
!     explicitly add the corrections for those terms (some of them are constant
!     due to constraints, but it's easier, safer, and not much more expensive
!     to compute them instantaneously)
      do i=1,nrexpolin(rs)
        ii = iaa(rs)%expolin(i,1)
        kk = iaa(rs)%expolin(i,2)
        dvec(1) = x(kk) - x(ii)
        dvec(2) = y(kk) - y(ii)
        dvec(3) = z(kk) - z(ii)
        dis2 = dvec(1)*dvec(1) + dvec(2)*dvec(2) + dvec(3)*dvec(3)
        dis = sqrt(dis2)
!
        term0 = pfac2*electric*atq(ii)*atq(kk)
        term1 = term0/dis
!
        if (lrel_md.eq.3) then
          terme = term0*(par_POLAR(1)*dis2-par_POLAR(2))
          evec(6) = evec(6) + terme
          term2 = -2.0*term0*par_POLAR(1)
          do j=1,3
            ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
            ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
          end do
        end if
      end do
!
    end if
!
!
! neighboring residues
! note that we assume implicitly here that only one of the two residues is ghosted(!!!!!)
!
  else if (abs(rs1-rs2).eq.1) then
!
    if (use_FEGS(6).EQV..true.) then
!
      if (rs1.gt.rs2) then
        rs = rs2
        call dis_bound_rs(rs2,rs1,svec)
      else
        rs = rs1
        call dis_bound_rs(rs1,rs2,svec)
      end if
!
!     loop over relevant nearest-neighbor-residue interactions
      do i=1,nrpolnb(rs)
        ii = iaa(rs)%polnb(i,1)
        kk = iaa(rs)%polnb(i,2)
        call dis_bound5(ii,kk,svec,dis2,dvec)
        if ((cut.EQV..true.).AND.(dis2.gt.mcel_cutoff2)) cycle
        dis = sqrt(dis2)
!
!       things are nice and simple in vacuo ...
        term0 = pfac2*fudge(rs)%elnb(i)*electric*atq(ii)*atq(kk)
        term1 = term0/dis
        if (lrel_md.eq.3) then
          terme = term1+term0*(par_POLAR(1)*dis2-par_POLAR(2))
          evec(6) = evec(6) + terme
          term2 = (1.0/dis2)*(term1 - 2.0*term0*dis2*par_POLAR(1))
          do j=1,3
            ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
            ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
          end do
        else
          call fexit()
        end if
!
!       extra correction for fudged terms
        if (fudge(rs)%elnb(i).lt.1.0) then
          term0 = pfac2*electric*atq(ii)*atq(kk)*(1.0-fudge(rs)%elnb(i))
          term1 = term0/dis
!
!          if (lrel_md.eq.2) then
!            terme = -term1*erf(ewpm*dis)
!            evec(6) = evec(6) + terme
!            term2 = (1.0/dis2)*(terme - &
! &                    term1*ewpite*dis*exp(-ewpm2*dis2))
!            do j=1,3
!              ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
!              ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
!            end do
          if (lrel_md.eq.3) then
            terme = term0*(par_POLAR(1)*dis2-par_POLAR(2))
            evec(6) = evec(6) + terme
            term2 = -2.0*term0*par_POLAR(1)
            do j=1,3
              ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
              ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
            end do
          end if
        end if
      end do
!
!     for RF we have to go through the exclusion list and
!     explicitly add the corrections for those terms (some of them are constant
!     due to constraints, but it's easier, safer, and not much more expensive
!     to compute them instantaneously)
!     note that we do not use minimum image here, as excluded neighbor
!     interactions are ALWAYS on the same molecule and in that case dis_bound_rs
!     ALWAYS returns a zero shift vector (svec) 
      do i=1,nrexpolnb(rs)
        ii = iaa(rs)%expolnb(i,1)
        kk = iaa(rs)%expolnb(i,2)
        call dis_bound5(ii,kk,svec,dis2,dvec)
        dis = sqrt(dis2)
!
        term0 = pfac2*electric*atq(ii)*atq(kk)
        term1 = term0/dis
        if (lrel_md.eq.3) then
          terme = term0*(par_POLAR(1)*dis2-par_POLAR(2))
          evec(6) = evec(6) + terme
          term2 = -2.0*term0*par_POLAR(1)
          do j=1,3
            ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
            ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
          end do
        end if
      end do
!
    end if
!
!
! all other residue pairs (again not we assume only one of the residues is ghosted(!!!!))
!
  else
!
    if (use_FEGS(6).EQV..true.) then
!
!     there is no topological relationship for remaining residues -> always check BC
      call dis_bound_rs(rs1,rs2,svec)
!
      do i=1,at(rs1)%npol
        do k=1,at(rs2)%npol
          ii = at(rs1)%pol(i)
          kk = at(rs2)%pol(k)
          call dis_bound5(ii,kk,svec,dis2,dvec)
          if ((cut.EQV..true.).AND.(dis2.gt.mcel_cutoff2)) cycle
          dis = sqrt(dis2)
!
!         things are nice and simple in vacuo ...
          term0 = pfac2*electric*atq(ii)*atq(kk)
          term1 = term0/dis
          if (lrel_md.eq.3) then
            terme = term1+term0*(par_POLAR(1)*dis2-par_POLAR(2))
            evec(6) = evec(6) + terme
            term2 = (1.0/dis2)*(term1 - 2.0*term0*dis2*par_POLAR(1))
            do j=1,3
              ca_f(ii,j) = ca_f(ii,j) - dvec(j)*term2
              ca_f(kk,j) = ca_f(kk,j) + dvec(j)*term2
            end do
          else
            call fexit()
          end if
        end do
      end do
    end if
!       
  end if
!
end
!
!---------------------------------------------------------------------------
!
subroutine genmu_dr(r1,r2,i_genmu,dr1,dr2)
!
  implicit none
!
  RTYPE dr1,dr2,r1,r2
  integer i_genmu
!
  if (i_genmu.eq.0) then
    dr1 = r2/(2.0*sqrt(r1*r2))
    dr2 = dr1*r1/r2
  else if (i_genmu.eq.1) then
    dr1 = 0.5
    dr2 = dr1
  else if (i_genmu.eq.-1) then
    dr1 = 1./(((r1*r1)*(1./r1 + 1./r2))*(0.5*(1./r1 + 1./r2)))
    dr2 = 1./(((r2*r2)*(1./r1 + 1./r2))*(0.5*(1./r1 + 1./r2)))
  else if (i_genmu.eq.2) then
    dr1 = r1/(2.0*sqrt(0.5*(r1*r1 + r2*r2)))
    dr2 = r2*dr1/r1
  else
    dr1 = ((0.5*(r1**(1.*i_genmu) + r2**(1.*i_genmu)))&
 &                                      **(1./(1.*i_genmu) - 1.0))*&
 &         (0.5*i_genmu*(r1**(1.*i_genmu - 1.0)))*(1./(1.*i_genmu))
    dr2 = ((0.5*(r1**(1.*i_genmu) + r2**(1.*i_genmu)))&
 &                                      **(1./(1.*i_genmu) - 1.0))*&
 &         (0.5*i_genmu*(r2**(1.*i_genmu - 1.0)))*(1./(1.*i_genmu))
  end if
!
end
!
!-----------------------------------------------------------------------
! 
! long-range Coulombic forces using neighbor lists: this uses two auxiliary functions:
! cgrp_ia(..) for monopoles reduced to atoms
! Vforce_rsp_long_lrel(..) for complete monopole-monopole and monopole-dipole interactions
!
subroutine force_P_LR_NBL(evec,ca_f,cycle_frz,sum_s)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use iounit
  use sequen
  use cutoffs
  use molecule
  use fyoc
!
  implicit none
!
  RTYPE svec(3),ca_f(n,3),evec(MAXENERGYTERMS)
  integer rs1,rs2,ii,i,j,dp2,imol,jmol
  logical docyc,cycle_frz
  RTYPE sum_s(n)
!
  if (lrel_md.eq.4) then
    do i=1,cglst%ncs
      ii = cglst%it(i)
      rs1 = atmres(cglst%it(i))
      imol = molofrs(rs1)
      if (rs1.lt.nseq) then
        docyc = .false.
        jmol = molofrs(rs1+1)
        if ((rsp_vec(rs1).eq.0).AND.(chgflag(rs1+1).EQV..true.)) then
          if (cycle_frz.EQV..true.) then
            if (imol.eq.jmol) then
              if (molfrzidx(imol).ge.3) docyc = .true.
              if ((molfrzidx(imol).eq.2).AND.(nchi(rs1).eq.0).AND.(nchi(rs1+1).eq.0)) docyc = .true.
            else
              if ((molfrzidx(imol).eq.4).AND.(molfrzidx(jmol).eq.4)) docyc = .true.
            end if
          end if
          if (docyc.EQV..false.) then
            call dis_bound_rs(rs1,rs1+1,svec)
            do dp2=1,at(rs1+1)%ndpgrps
              if (at(rs1+1)%dpgrp(dp2)%nc.eq.0) cycle
              call cgrp_ia(i,at(rs1+1)%dpgrp(dp2)%cgn,rs1,rs1+1,evec,ca_f,svec,sum_s)
            end do
          end if
        end if
!
        do j=1,rs_nbl(rs1)%nnblrs
          rs2 = rs_nbl(rs1)%nblr(j)
          jmol = molofrs(rs2)
          if (cycle_frz.EQV..true.) then
            if (imol.eq.jmol) then
              if (molfrzidx(imol).ge.3) cycle
              if ((molfrzidx(imol).eq.2).AND.(nchi(rs1).eq.0).AND.(nchi(rs1+1).eq.0)) cycle
            else
              if ((molfrzidx(imol).eq.4).AND.(molfrzidx(jmol).eq.4)) cycle
            end if
          end if
          do dp2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(dp2)%nc.eq.0) cycle
            svec(:) = rs_nbl(rs1)%lrsvec(j,:)
            call cgrp_ia(i,at(rs2)%dpgrp(dp2)%cgn,rs1,rs2,evec,ca_f,svec,sum_s)
          end do
        end do
!
        do j=1,rs_nbl(rs1)%ngnblrs
          rs2 = rs_nbl(rs1)%gnblr(j)
          jmol = molofrs(rs2)
          if (cycle_frz.EQV..true.) then
            if (imol.eq.jmol) then
              if (molfrzidx(imol).ge.3) cycle
              if ((molfrzidx(imol).eq.2).AND.(nchi(rs1).eq.0).AND.(nchi(rs1+1).eq.0)) cycle
            else
              if ((molfrzidx(imol).eq.4).AND.(molfrzidx(jmol).eq.4)) cycle
            end if
          end if
          do dp2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(dp2)%nc.eq.0) cycle
            svec(:) = rs_nbl(rs1)%glrsvec(j,:)
            call cgrp_ia(i,at(rs2)%dpgrp(dp2)%cgn,rs1,rs2,evec,ca_f,svec,sum_s)
          end do
        end do
      end if
    end do
!
  else if (lrel_md.eq.5) then
    do rs1=1,nseq
      imol = molofrs(rs1)
      if (rs1.lt.nseq) then
        docyc = .false.
        jmol = molofrs(rs1+1)
        if ((rsp_vec(rs1).eq.0).AND.((chgflag(rs1).EQV..true.).OR.(chgflag(rs1+1).EQV..true.))) then
          if (cycle_frz.EQV..true.) then
            if (imol.eq.jmol) then
              if (molfrzidx(imol).ge.3) docyc = .true.
              if ((molfrzidx(imol).eq.2).AND.(nchi(rs1).eq.0).AND.(nchi(rs1+1).eq.0)) docyc = .true.
            else
              if ((molfrzidx(imol).eq.4).AND.(molfrzidx(jmol).eq.4)) docyc = .true.
            end if
          end if
          if (docyc.EQV..false.) then
            call dis_bound_rs(rs1,rs1+1,svec)
            call Vforce_rsp_long_lrel(evec,rs1,rs1+1,svec,ca_f,sum_s)
          end if
        end if
!
        do j=1,rs_nbl(rs1)%nnblrs
          rs2 = rs_nbl(rs1)%nblr(j)
          jmol = molofrs(rs2)
          if (cycle_frz.EQV..true.) then
            if (imol.eq.jmol) then
              if (molfrzidx(imol).ge.3) cycle
              if ((molfrzidx(imol).eq.2).AND.(nchi(rs1).eq.0).AND.(nchi(rs1+1).eq.0)) cycle
            else
              if ((molfrzidx(imol).eq.4).AND.(molfrzidx(jmol).eq.4)) cycle
            end if
          end if
          call dis_bound_rs(rs1,rs2,svec)
          call Vforce_rsp_long_lrel(evec,rs1,rs2,svec,ca_f,sum_s)
        end do
!
        do j=1,rs_nbl(rs1)%ngnblrs
          rs2 = rs_nbl(rs1)%gnblr(j)
          jmol = molofrs(rs2)
          if (cycle_frz.EQV..true.) then
            if (imol.eq.jmol) then
              if (molfrzidx(imol).ge.3) cycle
              if ((molfrzidx(imol).eq.2).AND.(nchi(rs1).eq.0).AND.(nchi(rs1+1).eq.0)) cycle
            else
              if ((molfrzidx(imol).eq.4).AND.(molfrzidx(jmol).eq.4)) cycle
            end if
          end if
          call dis_bound_rs(rs1,rs2,svec)
          call Vforce_rsp_long_lrel(evec,rs1,rs2,svec,ca_f,sum_s)
        end do

      end if
    end do
  else
    write(ilog,*) 'Fatal. Called force_P_LR_NBL(...) with unsupported long-range treatment (got: ',lrel_md,'). This is a bug.'
    call fexit()
  end if
!
  end
!
!-----------------------------------------------------------------------
!
! this does not vectorize, so is listed here: interaction of two monopoles mapped to a single atom each
!
  subroutine cgrp_ia(i,j,rs1,rs2,evec,ca_f,svec,sum_s)
!
  use atoms
  use forces
  use energies
  use cutoffs
  use units
  use iounit
!
  implicit none
!
  integer ii,kk,l,i,j,rs1,rs2
  RTYPE dvec(3),dis2,svec(3),id2,d1,dd2,dd1,foin,id1,dum1
  RTYPE evec(MAXENERGYTERMS),ca_f(n,3)
  RTYPE term0,term1,term2,term3,term4,termq,terms,termr,termqi,termqk,pfac,pfac2,pfac3,pfac4
  RTYPE for_tmpi(3),for_tmpk(3),termf(3),tsdk(3),tsdi(3)
  RTYPE genmu
  logical dogh
  RTYPE sum_s(n)
!
  if (use_FEG.EQV..true.) then
    if (((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..false.)).OR.&
 &    ((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..true.).AND.(fegmode.eq.1))) then
      dogh = .false.
    else
      dogh = .true.
    end if
  else
    dogh = .false.
  end if
  ii = cglst%it(i)
  kk = cglst%it(j)
  dvec(1) = x(kk) - x(ii) + svec(1)
  dvec(2) = y(kk) - y(ii) + svec(2)
  dvec(3) = z(kk) - z(ii) + svec(3)
  dis2 = dvec(1)**2 + dvec(2)**2 + dvec(3)**2
  d1 = sqrt(dis2)
!  write(*,*) 'found ',rs1,rs2,sqrt(dis2)
  id1 = 1.0/d1
  id2 = id1*id1
  if (dogh.EQV..false.) then
    if (use_IMPSOLV.EQV..true.) then
!
!     a bit complicated
      pfac = electric*scale_POLAR
      termq = cglst%tc(i)*cglst%tc(j)
      terms = scrq(ii)*scrq(kk)
      termr = atr(ii)+atr(kk)
!
      if ((scrq_model.eq.1).OR.(scrq_model.eq.2)) then
        termqk = termq*scrq(ii)
        termqi = termq*scrq(kk)
        term0 = pfac*id1
        termqk = termqk*term0
        termqi = termqi*term0
        term0 = term0*termq
        for_tmpi(:) = term0*scrq(kk)*scrq_dr(ii,:)
        for_tmpk(:) = term0*scrq(ii)*scrq_dr(kk,:)
        term1 = term0*terms
        evec(6) = evec(6) + term1
        term2 = term1*id2
        termf(:) = term2*dvec(:)
        ca_f(ii,:) = ca_f(ii,:) - termf(:) + for_tmpi(:)
        ca_f(kk,:) = ca_f(kk,:) + termf(:) + for_tmpk(:)
        sum_s(ii) = sum_s(ii) + termqi
        sum_s(kk) = sum_s(kk) + termqk
!
      else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
        if ((scrq(ii).ne.0.0).AND.(scrq(kk).ne.0.0)) then
          terms = genmu(scrq(ii),scrq(kk),i_sqm)
          call genmu_dr(scrq(ii),scrq(kk),i_sqm,termqi,termqk)
          tsdi(:) = scrq_dr(ii,:)*termqi
          tsdk(:) = scrq_dr(kk,:)*termqk
          termqi = termqi*termq
          termqk = termqk*termq
          term0 = pfac*id1
          termqk = termqk*term0
          termqi = termqi*term0
          term0 = term0*termq
          for_tmpi(:) = term0*tsdi(:)
          for_tmpk(:) = term0*tsdk(:)
          term1 = term0*terms
          evec(6) = evec(6) + term1
          term2 = term1*id2
          termf(:) = term2*dvec(:)
          ca_f(ii,:) = ca_f(ii,:) - termf(:) + for_tmpi(:)
          ca_f(kk,:) = ca_f(kk,:) + termf(:) + for_tmpk(:)
          sum_s(ii) = sum_s(ii) + termqi
          sum_s(kk) = sum_s(kk) + termqk
        else 
!         do nothing (this should never happen of course)
        end if
!
      else if (scrq_model.eq.4) then
        pfac2 = par_IMPSOLV(8)*electric*scale_POLAR
!       contact regime
        if (d1.lt.termr) then
          term1 = 1.0/termr
          term0 = term1
!       distance-dependent regime
        else
          term1 = id1
          term0 = 2.0*id1
        end if
        term2 = pfac2*termq*id1
        termqi = term2*id2*term0
        term2 = term2*term1
        termf(:) = termqi*dvec(:)
        evec(6) = evec(6) + term2
        ca_f(ii,:) = ca_f(ii,:) - termf(:)
        ca_f(kk,:) = ca_f(kk,:) + termf(:) 
!
      else  !screening models #3,7,8,9
        pfac2 = par_IMPSOLV(8)*electric
        pfac3 = 1./par_IMPSOLV(1)
        pfac4 = par_IMPSOLV(9)*pfac3
        if ((scrq_model.eq.7).OR.(scrq_model.eq.8)) then
          if ((scrq(ii).ne.0.0).AND.(scrq(kk).ne.0.0)) then
            terms = genmu(scrq(ii),scrq(kk),i_sqm)
            call genmu_dr(scrq(ii),scrq(kk),i_sqm,termqi,termqk)
            tsdi(:) = scrq_dr(ii,:)*termqi
            tsdk(:) = scrq_dr(kk,:)*termqk
            termqi = termqi*termq
            termqk = termqk*termq
          else
!         
            tsdi(:) = 0.0
            tsdk(:) = 0.0
            termqi = 0.0
            termqk = 0.0
            termq = 0.0
            pfac = 0.0
          end if
        else
          termqk = termq*scrq(ii)
          termqi = termq*scrq(kk)
          tsdk(:) = scrq_dr(kk,:)*scrq(ii)
          tsdi(:) = scrq_dr(ii,:)*scrq(kk)
        end if
        term0 = pfac*id1
        term1 = electric*termq*terms
        term2 = pfac2*termq/termr
        for_tmpi(:) = 0.0
        if (abs(term1).gt.abs(term2)) then
          term4 = term1
          term3 = 1.0
        else
          dd1 = (d1 - termr)*pfac3
          dd2 = 1.0 - dd1
          if (dd1.lt.0.0) then
            term3 = (1.0 - par_IMPSOLV(9))
            term4 = term3*term1 + par_IMPSOLV(9)*term2
            term0 = term3*term0
          else if (dd2.gt.0.0) then
            dum1 = par_IMPSOLV(9)*dd2
            term3 = (1.0 - dum1)
            term4 = term3*term1 + dum1*term2
            term0 = term3*term0
            dum1 = -scale_POLAR*id2*(term2-term1)*pfac4
            for_tmpi(:) = dum1*dvec(:)
          else
            term4 = term1
            term3 = 1.0
          end if
        end if
        sum_s(ii) = sum_s(ii) + term0*termqi
        sum_s(kk) = sum_s(kk) + term0*termqk
        term0 = pfac*term3*termq
        for_tmpk(:) = -for_tmpi(:) + term0*tsdk(:)*id1
        for_tmpi(:) = for_tmpi(:) + term0*tsdi(:)*id1
        term1 = scale_POLAR*term4*id1
        evec(6) = evec(6) + term1
        term2 = term1*id2
        termf(:) = term2*dvec(:)
        ca_f(ii,:) = ca_f(ii,:) - termf(:) + for_tmpi(:)
        ca_f(kk,:) = ca_f(kk,:) + termf(:) + for_tmpk(:)
!
      end if ! choice of scrq_model
!
    else
!
!     much easier
      term1 = electric*cglst%tc(i)*cglst%tc(j)*scale_POLAR*id1
      evec(6) = evec(6) + term1
      do l=1,3
        foin = term1*dvec(l)*id2
        ca_f(ii,l) = ca_f(ii,l) - foin
        ca_f(kk,l) = ca_f(kk,l) + foin
      end do
    end if
!
  else
    if (use_IMPSOLV.EQV..true.) then
      write(ilog,*) 'Fatal. Missing support for combination of FEG and the ABSINTH&
 & implicit sovlation model in cgrp_ia. This is an omission bug.'
      call fexit()
    else
      term1 = electric*cglst%tc(i)*cglst%tc(j)*par_FEG2(9)/&
 & (par_FEG2(10) + d1)
      evec(6) = evec(6) + term1
      do l=1,3
        foin = term1*dvec(l)*id1/(par_FEG2(10) + d1)
        ca_f(ii,l) = ca_f(ii,l) - foin
        ca_f(kk,l) = ca_f(kk,l) + foin
      end do
    end if
  end if
!
  end
!
!-----------------------------------------------------------------------
!
! compute forces due to boundary
!
subroutine force_boundary_rs(rs,evec,mode,ca_f)
!
  use iounit
  use atoms
  use system
  use aminos
  use sequen
  use polypep
  use energies
  use forces
!
  implicit none
!
  integer rs,mode,i,ii,j
  RTYPE dis2,dis,evec(MAXENERGYTERMS),dvec(3),ddd,ca_f(n,3)
!
  if (bnd_type.eq.1) then
    if (bnd_shape.eq.3) then ! periodic cylinder w/ axis along z (continuous tube) using ASWBC
      do i=1,at(rs)%nbb+at(rs)%nsc
        if (i.le.at(rs)%nbb) then
          ii = at(rs)%bb(i)
        else
          ii = at(rs)%sc(i-at(rs)%nbb)
        end if
        dvec(1) = (x(ii)-bnd_params(1))
        dvec(2) = (y(ii)-bnd_params(2))
        dis2 = dvec(1)*dvec(1) + dvec(2)*dvec(2)
        dis = sqrt(dis2)
        ddd = dis - bnd_params(4)
        if (ddd.gt.0.0) then
          evec(12) = evec(12) + bnd_params(7)*ddd*ddd
          do j=1,2
            ca_f(ii,j) = ca_f(ii,j) - 2.0*bnd_params(7)*ddd*(dvec(j)/dis)
          end do
        end if
      end do
    else
!     do nothing
    end if
  else if (bnd_type.eq.3) then
!   since we use the reference atom for the particular residue to determine the boundary violation,
!   forces are trivial (note that that force will show up on the center of mass in internal coordinate
!   space of course)
!
    ii = refat(rs)
    if (bnd_shape.eq.1) then
      dvec(1) = (x(ii)-bnd_params(4))
      dvec(2) = (y(ii)-bnd_params(5))
      dvec(3) = (z(ii)-bnd_params(6))
      do j=1,3
        if (dvec(j).lt.0.0) then
          evec(12) = evec(12) + bnd_params(7)*dvec(j)*dvec(j)
          ca_f(ii,j) = ca_f(ii,j) - 2.0*bnd_params(7)*dvec(j)
        else if (dvec(j).gt.bnd_params(j)) then
          evec(12) = evec(12) + bnd_params(7)*(dvec(j)-bnd_params(j))**2
          ca_f(ii,j) = ca_f(ii,j) - 2.0*bnd_params(7)*(dvec(j)-bnd_params(j))
        end if
      end do
    else if (bnd_shape.eq.2) then
      dvec(1) = (x(ii)-bnd_params(1))
      dvec(2) = (y(ii)-bnd_params(2))
      dvec(3) = (z(ii)-bnd_params(3))
      dis2 = dvec(1)*dvec(1) + dvec(2)*dvec(2) + dvec(3)*dvec(3)
      dis = sqrt(dis2)
      ddd = dis - bnd_params(4)
      if (ddd.gt.0.0) then
        evec(12) = evec(12) + bnd_params(7)*ddd*ddd
        do j=1,3
          ca_f(ii,j) = ca_f(ii,j) - &
 &                    2.0*bnd_params(7)*ddd*(dvec(j)/dis)
          bnd_f(j,1) = bnd_f(j,1) + &
 &       2.0*bnd_params(7)*ddd*(dvec(j)/dis)*dvec(1)
          bnd_f(j,2) = bnd_f(j,2) + &
 &       2.0*bnd_params(7)*ddd*(dvec(j)/dis)*dvec(2)
          bnd_f(j,3) = bnd_f(j,3) + &
 &       2.0*bnd_params(7)*ddd*(dvec(j)/dis)*dvec(3)
        end do
        bnd_fr = bnd_fr + 2.0*bnd_params(7)*ddd!*bnd_params(4)
      end if
    else if (bnd_shape.eq.3) then ! cylinder w/ axis along z
      dvec(1) = (x(ii)-bnd_params(1))
      dvec(2) = (y(ii)-bnd_params(2))
      dvec(3) = (z(ii)-bnd_params(3)+0.5*bnd_params(6))
      dis2 = dvec(1)*dvec(1) + dvec(2)*dvec(2)
      dis = sqrt(dis2)
      ddd = dis - bnd_params(4)
      if (ddd.gt.0.0) then
        evec(12) = evec(12) + bnd_params(7)*ddd*ddd
        do j=1,2
          ca_f(ii,j) = ca_f(ii,j) - 2.0*bnd_params(7)*ddd*(dvec(j)/dis)
        end do
      end if
      if (dvec(3).lt.0.0) then
        evec(12) = evec(12) + bnd_params(7)*dvec(3)*dvec(3)
        ca_f(ii,3) = ca_f(ii,3) - 2.0*bnd_params(7)*dvec(3)
      else if (dvec(3).gt.bnd_params(6)) then
        evec(12) = evec(12) + bnd_params(7)*(dvec(3)-bnd_params(6))**2
        ca_f(ii,3) = ca_f(ii,3) - 2.0*bnd_params(7)*(dvec(3)-bnd_params(6))
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in force_b&
 &oundary_rs() (code # ',bnd_shape,').'
      call fexit()
    end if
  else if (bnd_type.eq.4) then
    if (bnd_shape.eq.1) then
      do i=1,at(rs)%nbb+at(rs)%nsc
        if (i.le.at(rs)%nbb) then
          ii = at(rs)%bb(i)
        else
          ii = at(rs)%sc(i-at(rs)%nbb)
        end if
        dvec(1) = (x(ii)-bnd_params(4))
        dvec(2) = (y(ii)-bnd_params(5))
        dvec(3) = (z(ii)-bnd_params(6))
        do j=1,3
          if (dvec(j).lt.0.0) then
            evec(12) = evec(12) + bnd_params(7)*dvec(j)*dvec(j)
            ca_f(ii,j) = ca_f(ii,j) - 2.0*bnd_params(7)*dvec(j)
          else if (dvec(j).gt.bnd_params(j)) then
            evec(12) = evec(12) + bnd_params(7)*(dvec(j)-bnd_params(j))**2
            ca_f(ii,j) = ca_f(ii,j) - 2.0*bnd_params(7)*(dvec(j)-bnd_params(j))
          end if
        end do
      end do
    else if (bnd_shape.eq.2) then
      do i=1,at(rs)%nbb+at(rs)%nsc
        if (i.le.at(rs)%nbb) then
          ii = at(rs)%bb(i)
        else
          ii = at(rs)%sc(i-at(rs)%nbb)
        end if
        dvec(1) = (x(ii)-bnd_params(1))
        dvec(2) = (y(ii)-bnd_params(2))
        dvec(3) = (z(ii)-bnd_params(3))
        dis2 = dvec(1)*dvec(1) + dvec(2)*dvec(2) + dvec(3)*dvec(3)
        dis = sqrt(dis2)
        ddd = dis - bnd_params(4)
        if (ddd.gt.0.0) then
          evec(12) = evec(12) + bnd_params(7)*ddd*ddd
          do j=1,3
            ca_f(ii,j) = ca_f(ii,j) - &
 &                    2.0*bnd_params(7)*ddd*(dvec(j)/dis)
            bnd_f(j,1) = bnd_f(j,1) + &
 &       2.0*bnd_params(7)*ddd*(dvec(j)/dis)*dvec(1)
            bnd_f(j,2) = bnd_f(j,2) + &
 &       2.0*bnd_params(7)*ddd*(dvec(j)/dis)*dvec(2)
            bnd_f(j,3) = bnd_f(j,3) + &
 &       2.0*bnd_params(7)*ddd*(dvec(j)/dis)*dvec(3)
          end do
          bnd_fr = bnd_fr + 2.0*bnd_params(7)*ddd
        end if
      end do
    else if (bnd_shape.eq.3) then ! cylinder w/ axis along z
      do i=1,at(rs)%nbb+at(rs)%nsc
        if (i.le.at(rs)%nbb) then
          ii = at(rs)%bb(i)
        else
          ii = at(rs)%sc(i-at(rs)%nbb)
        end if
        dvec(1) = (x(ii)-bnd_params(1))
        dvec(2) = (y(ii)-bnd_params(2))
        dvec(3) = (z(ii)-bnd_params(3)+0.5*bnd_params(6))
        dis2 = dvec(1)*dvec(1) + dvec(2)*dvec(2)
        dis = sqrt(dis2)
        ddd = dis - bnd_params(4)
        if (ddd.gt.0.0) then
          evec(12) = evec(12) + bnd_params(7)*ddd*ddd
          do j=1,2
            ca_f(ii,j) = ca_f(ii,j) - 2.0*bnd_params(7)*ddd*(dvec(j)/dis)
          end do
        end if
        if (dvec(3).lt.0.0) then
          evec(12) = evec(12) + bnd_params(7)*dvec(3)*dvec(3)
          ca_f(ii,3) = ca_f(ii,3) - 2.0*bnd_params(7)*dvec(3)
        else if (dvec(3).gt.bnd_params(6)) then
          evec(12) = evec(12) + bnd_params(7)*(dvec(3)-bnd_params(6))**2
          ca_f(ii,3) = ca_f(ii,3) - 2.0*bnd_params(7)*(dvec(3)-bnd_params(6))
        end if
      end do
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in force_b&
 &oundary_rs() (code # ',bnd_shape,').'
      call fexit()
    end if
  else if (bnd_type.eq.2) then
    write(ilog,*) 'Fatal. Force calculation does not support hard-wa&
 &ll boundary conditions at the moment.'
    call fexit()
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition&
 & in force_boundary_rs() (code # ',bnd_type,').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------------
!
! forces for atom-atom distance restraints
!
subroutine force_drest(evec,ca_f)
!
  use energies
  use distrest
  use atoms
  use forces
  use sequen
!
  implicit none
!
  RTYPE edr,dis,evec(MAXENERGYTERMS),svec(3),dvec(3),dis2,term0
  RTYPE ca_f(n,3)
  integer i,j,ii,kk
!
  edr = 0.0

  do i=1,ndrest
!
    ii = dresat(i,1)
    kk = dresat(i,2)
    if (molofrs(atmres(ii)).ne.molofrs(atmres(kk))) then
      call dis_bound_rs(atmres(ii),atmres(kk),svec)
      call dis_bound5(ii,kk,svec,dis2,dvec)
      dis = sqrt(dis2)
    else
      dvec(1) = x(kk)-x(ii)
      dvec(2) = y(kk)-y(ii)
      dvec(3) = z(kk)-z(ii)
      dis = sqrt(dvec(1)*dvec(1)+dvec(2)*dvec(2)+dvec(3)*dvec(3))
    end if
    if (drestyp(i).eq.1) then
      edr = edr + dresprm(i,2)*(dis-dresprm(i,1))*(dis-dresprm(i,1))
      term0 = scale_DREST*2.0*dresprm(i,2)*(dis-dresprm(i,1))/dis
      do j=1,3
        ca_f(ii,j) = ca_f(ii,j) + term0*dvec(j)
        ca_f(kk,j) = ca_f(kk,j) - term0*dvec(j)
      end do
    else if (drestyp(i).eq.2) then
      if (dis.gt.dresprm(i,1)) then
        edr = edr + dresprm(i,2)*(dis-dresprm(i,1))*(dis-dresprm(i,1))
        term0 = scale_DREST*2.0*dresprm(i,2)*(dis-dresprm(i,1))/dis
        do j=1,3
          ca_f(ii,j) = ca_f(ii,j) + term0*dvec(j)
          ca_f(kk,j) = ca_f(kk,j) - term0*dvec(j)
        end do
      end if
    else if (drestyp(i).eq.3) then
      if (dis.lt.dresprm(i,1)) then
        edr = edr + dresprm(i,2)*(dis-dresprm(i,1))*(dis-dresprm(i,1))
        term0 = scale_DREST*2.0*dresprm(i,2)*(dis-dresprm(i,1))/dis
        do j=1,3
          ca_f(ii,j) = ca_f(ii,j) + term0*dvec(j)
          ca_f(kk,j) = ca_f(kk,j) - term0*dvec(j)
        end do
      end if
    end if
!
  end do
!
  evec(10) = evec(10) + scale_DREST*edr
!
end
!
!-----------------------------------------------------------------------
!
subroutine force_corrector(rs,evec,ca_f)

  use energies
  use aminos
  use sequen
  use fyoc
  use movesets
  use iounit
  use forces
  use system
  use atoms
!
  implicit none
!
  RTYPE evec(MAXENERGYTERMS),fomega,fphenol,ca_f(n,3)
  integer rs
  character(3) resname
!
  resname = amino(seqtyp(rs))
!
! omega torsions (i.e., Z/E amide isomerization)
  if (wline(rs).gt.0) then
   call get_fomega(rs,resname,fomega,ca_f)
   evec(2) = evec(2) + scale_CORR*fomega
  end if
!
!! proline puckering
!  if (fycxyz.ne.2) then
!    if ((resname.eq.'PRO').OR.(resname.eq.'HYP')&
! &                   .OR.(resname.eq.'PCA')) then
!      if ((dc_di(molofrs(rs))%frz(chinr(1,rs)).EQV..false.).OR.&
! &        (dc_di(molofrs(rs))%frz(chinr(2,rs)).EQV..false.)) then
!      write(ilog,*) 'Fatal. Proline puckering not yet supported by f&
! &orce calculations. Use constraints to circumvent this issue.'
!        call fexit()
!       else
!!        do nothing
!       end if
!!      evec(2) = evec(2) + scale_CORR*epucker(rs)
!    end if
!  end if
!
! phenolic polar hydrogen Z/E
  if ((resname.eq.'PCR').OR.(resname.eq.'TYR')) then
    call get_fphenol(rs,resname,fphenol,ca_f)
    evec(2) = evec(2) + scale_CORR*fphenol
  end if
!  write(*,*) chi(3,rs),evec(2)
!
end
!
!-----------------------------------------------------------------------
!
! see eomega
!
subroutine get_fomega(rs,resname,fomega,ca_f)
!
  use fyoc
  use math
  use energies
  use sequen
  use polypep
  use forces
  use system
  use atoms
  use wl
!
  implicit none
!
  integer rs,imol,ito,j,shf
  character(3) resname
  RTYPE fomega,t1,t2,t3,t4,getztor,ca_f(n,3)
  RTYPE cdt,cd1(3),cd2(3),cd3(3),cd4(3),cost,si
  logical t2flag
!
! note that we'll do this only for secondary amides, primary/tertiary ones are a bit too redundant
!
  t2flag = .false.
  imol = molofrs(rs)
!  write(*,*) t1-t2,t3-t4,t2-t4,t1-t3
  if (resname.eq.'NMF') then
    t1=omega(rs)
    if (ua_model.gt.0) then
      t2=getztor(at(rs)%bb(2),at(rs)%bb(1),at(rs)%bb(3),at(rs)%bb(4))
      t2flag = .true.
    else
      t2=getztor(at(rs)%bb(2),at(rs)%bb(1),at(rs)%bb(3),at(rs)%bb(5))
      t3=getztor(at(rs)%bb(4),at(rs)%bb(1),at(rs)%bb(3),at(rs)%sc(1))
      t4=getztor(at(rs)%bb(4),at(rs)%bb(1),at(rs)%bb(3),at(rs)%bb(5))
     end if
!    write(*,*) t1-t2,t3-t4,t2-t4,t1-t3
  else if (resname.eq.'NMA') then
    t1=omega(rs)
    t2=getztor(at(rs)%bb(2),at(rs)%bb(1),at(rs)%bb(3),at(rs)%bb(4))
    t3=getztor(at(rs)%sc(2),at(rs)%bb(1),at(rs)%bb(3),at(rs)%sc(1))
    t4=getztor(at(rs)%sc(2),at(rs)%bb(1),at(rs)%bb(3),at(rs)%bb(4))
  else if (resname.eq.'PRO') then
    shf = 0
    if (ua_model.gt.0) shf = 1
    t1 = getztor(oi(rs-1),ci(rs-1),ni(rs),cai(rs))
    t2 = getztor(oi(rs-1),ci(rs-1),ni(rs),at(rs)%sc(4-shf))
    if ((ua_model.gt.0).AND.(seqtyp(rs-1).eq.28)) then
      fomega = (6.089/2.0)*(1.0 - cos(2.0*t1/RADIAN)) + &
 &             (2.300/2.0)*(1.0 + cos(1.0*t2/RADIAN)) + &
 &             (6.089/2.0)*(1.0 - cos(2.0*t2/RADIAN))
      t2flag = .true.
    else
      t3 = omega(rs)
      if (seqtyp(rs-1).eq.28) then
        t4 = getztor(at(rs-1)%bb(3),ci(rs-1),ni(rs),at(rs)%sc(4-shf))
      else
        t4 = getztor(cai(rs-1),ci(rs-1),ni(rs),at(rs)%sc(4-shf))
      end if
      fomega = (6.089/2.0)*(1.0 - cos(2.0*t1/RADIAN)) +&
 &         (6.089/2.0)*(1.0 - cos(2.0*t2/RADIAN)) +&
 &         (2.300/2.0)*(1.0 + cos(1.0*t3/RADIAN)) +&
 &         (6.089/2.0)*(1.0 - cos(2.0*t3/RADIAN)) +&
 &         (2.300/2.0)*(1.0 + cos(1.0*t4/RADIAN)) +&
 &         (6.089/2.0)*(1.0 - cos(2.0*t4/RADIAN))
    end if
!   note that the ti-dihedrals are all linearly related to omega, i.e.,
!   dti/dw = 1 throughout, since the bond is always defined the same direction
    if ((fycxyz.eq.1).AND.(hmjam%isin(2).EQV..false.).AND.(hmjam%isin(11).EQV..false.)) then
      if (t2flag.EQV..true.) then
        dc_di(imol)%f(wnr(rs)) = dc_di(imol)%f(wnr(rs)) + scale_CORR*(&
 & - (6.089/2.0)*(2.0*(sin(2.0*t1/RADIAN) + sin(2.0*t2/RADIAN))) &
 &                    + (2.300/2.0)*sin(t2/RADIAN))
      else
        dc_di(imol)%f(wnr(rs)) = dc_di(imol)%f(wnr(rs)) + scale_CORR*(&
 & - (6.089/2.0)*(2.0*(sin(2.0*t1/RADIAN) + sin(2.0*t2/RADIAN) +&
 &                     sin(2.0*t3/RADIAN) + sin(2.0*t4/RADIAN)))&
 & + (2.300/2.0)*(sin(t3/RADIAN) + sin(t4/RADIAN)))
      end if
    else
      call onetor_cosderiv(oi(rs-1),ci(rs-1),ni(rs),cai(rs),&
 &cost,cd1,cd2,cd3,cd4,si)
      cdt = -scale_CORR*6.089*2.0*cost
      do j=1,3
        ca_f(oi(rs-1),j)=ca_f(oi(rs-1),j) - cdt*cd1(j)
        ca_f(ci(rs-1),j)=ca_f(ci(rs-1),j) - cdt*cd2(j)
        ca_f(ni(rs),j)=ca_f(ni(rs),j) - cdt*cd3(j)
        ca_f(cai(rs),j)=ca_f(cai(rs),j) - cdt*cd4(j)
      end do
      call onetor_cosderiv(oi(rs-1),ci(rs-1),ni(rs),at(rs)%sc(4-shf),&
 &cost,cd1,cd2,cd3,cd4,si)
      if (t2flag.EQV..true.) then
        cdt = scale_CORR*((2.300/2.0) - 6.089*2.0*cost)
      else
        cdt = -scale_CORR*6.089*2.0*cost
      end if
      do j=1,3
        ca_f(oi(rs-1),j)=ca_f(oi(rs-1),j) - cdt*cd1(j)
        ca_f(ci(rs-1),j)=ca_f(ci(rs-1),j) - cdt*cd2(j)
        ca_f(ni(rs),j)=ca_f(ni(rs),j) - cdt*cd3(j)
        ca_f(at(rs)%sc(4-shf),j)=ca_f(at(rs)%sc(4-shf),j) - cdt*cd4(j)
      end do
      if (t2flag.EQV..true.) return
      if (seqtyp(rs-1).eq.28) then
        ito = at(rs-1)%bb(3)
      else 
        ito = cai(rs-1)
      end if
      call onetor_cosderiv(ito,ci(rs-1),ni(rs),cai(rs),&
 &cost,cd1,cd2,cd3,cd4,si)
      cdt = scale_CORR*((2.300/2.0) - 6.089*2.0*cost)
      do j=1,3
        ca_f(ito,j)=ca_f(ito,j) - cdt*cd1(j)
        ca_f(ci(rs-1),j)=ca_f(ci(rs-1),j) - cdt*cd2(j)
        ca_f(ni(rs),j)=ca_f(ni(rs),j) - cdt*cd3(j)
        ca_f(cai(rs),j)=ca_f(cai(rs),j) - cdt*cd4(j)
      end do
      call onetor_cosderiv(ito,ci(rs-1),ni(rs),at(rs)%sc(4-shf),&
 &cost,cd1,cd2,cd3,cd4,si)
      cdt = scale_CORR*((2.300/2.0) - 6.089*2.0*cost)
      do j=1,3
        ca_f(ito,j)=ca_f(ito,j) - cdt*cd1(j)
        ca_f(ci(rs-1),j)=ca_f(ci(rs-1),j) - cdt*cd2(j)
        ca_f(ni(rs),j)=ca_f(ni(rs),j) - cdt*cd3(j)
        ca_f(at(rs)%sc(4-shf),j)=ca_f(at(rs)%sc(4-shf),j) - cdt*cd4(j)
      end do
    end if
    return
  else
    t1 = getztor(oi(rs-1),ci(rs-1),ni(rs),cai(rs))
    t2 = getztor(oi(rs-1),ci(rs-1),ni(rs),hni(rs))
    if ((ua_model.gt.0).AND.(seqtyp(rs-1).eq.28)) then
      t2flag = .true.
    else
      t3 = omega(rs)
      if (seqtyp(rs-1).eq.28) then
        t4 = getztor(at(rs-1)%bb(3),ci(rs-1),ni(rs),hni(rs))
      else
        t4 = getztor(cai(rs-1),ci(rs-1),ni(rs),hni(rs))
      end if
    end if
  end if
!
  if (t2flag.EQV..true.) then
    fomega = (6.089/2.0)*(1.0 - cos(2.0*t1/RADIAN)) +&
 &         (2.300/2.0)*(1.0 + cos(1.0*t2/RADIAN)) +&
 &         (6.089/2.0)*(1.0 - cos(2.0*t2/RADIAN))
  else
    fomega = (6.089/2.0)*(1.0 - cos(2.0*t1/RADIAN)) +&
 &         (4.900/2.0)*(1.0 - cos(2.0*t2/RADIAN)) +&
 &         (2.300/2.0)*(1.0 + cos(1.0*t3/RADIAN)) +&
 &         (6.089/2.0)*(1.0 - cos(2.0*t3/RADIAN)) +&
 &         (4.900/2.0)*(1.0 - cos(2.0*t4/RADIAN))
  end if
! it's easy in torsional MD
  if ((fycxyz.eq.1).AND.(hmjam%isin(2).EQV..false.).AND.(hmjam%isin(11).EQV..false.)) then
    if (t2flag.EQV..true.) then
      dc_di(imol)%f(wnr(rs)) = dc_di(imol)%f(wnr(rs)) + scale_CORR*(&
 & - (6.089/2.0)*(2.0*(sin(2.0*t1/RADIAN) + sin(2.0*t2/RADIAN)))&
 & + (2.300/2.0)*(sin(t2/RADIAN)))
    else
      dc_di(imol)%f(wnr(rs)) = dc_di(imol)%f(wnr(rs)) + scale_CORR*(&
 & - (6.089/2.0)*(2.0*(sin(2.0*t1/RADIAN) + sin(2.0*t3/RADIAN)))&
 & + (2.300/2.0)*(sin(t3/RADIAN))&
 & - (4.900/2.0)*(2.0*(sin(2.0*t2/RADIAN) + sin(2.0*t4/RADIAN))))
    end if
! messy in Cartesian MD
  else
    if (resname.eq.'NMF') then
      call onetor_cosderiv(at(rs)%bb(2),at(rs)%bb(1),at(rs)%bb(3),&
 &at(rs)%sc(1),cost,cd1,cd2,cd3,cd4,si)
      cdt = -scale_CORR*6.089*2.0*cost
      do j=1,3
        ca_f(at(rs)%bb(2),j)=ca_f(at(rs)%bb(2),j) - cdt*cd1(j)
        ca_f(at(rs)%bb(1),j)=ca_f(at(rs)%bb(1),j) - cdt*cd2(j)
        ca_f(at(rs)%bb(3),j)=ca_f(at(rs)%bb(3),j) - cdt*cd3(j)
        ca_f(at(rs)%sc(1),j)=ca_f(at(rs)%sc(1),j) - cdt*cd4(j)
      end do
      if (t2flag.EQV..true.) then
        call onetor_cosderiv(at(rs)%bb(2),at(rs)%bb(1),at(rs)%bb(3),&
 &at(rs)%bb(4),cost,cd1,cd2,cd3,cd4,si)
        cdt = scale_CORR*((2.300/2.0) - 6.089*2.0*cost)
      else
        call onetor_cosderiv(at(rs)%bb(2),at(rs)%bb(1),at(rs)%bb(3),&
 &at(rs)%bb(5),cost,cd1,cd2,cd3,cd4,si)
        cdt = -scale_CORR*4.900*2.0*cost
      end if
      do j=1,3
        ca_f(at(rs)%bb(2),j)=ca_f(at(rs)%bb(2),j) - cdt*cd1(j)
        ca_f(at(rs)%bb(1),j)=ca_f(at(rs)%bb(1),j) - cdt*cd2(j)
        ca_f(at(rs)%bb(3),j)=ca_f(at(rs)%bb(3),j) - cdt*cd3(j)
        if (t2flag.EQV..true.) then
          ca_f(at(rs)%bb(4),j)=ca_f(at(rs)%bb(4),j) - cdt*cd4(j)
        else
          ca_f(at(rs)%bb(5),j)=ca_f(at(rs)%bb(5),j) - cdt*cd4(j)
        end if
      end do
      if (t2flag.EQV..true.) return
      call onetor_cosderiv(at(rs)%bb(4),at(rs)%bb(1),at(rs)%bb(3),&
 &at(rs)%sc(1),cost,cd1,cd2,cd3,cd4,si)
      cdt = scale_CORR*((2.300/2.0) - 6.089*2.0*cost)
      do j=1,3
        ca_f(at(rs)%bb(4),j)=ca_f(at(rs)%bb(4),j) - cdt*cd1(j)
        ca_f(at(rs)%bb(1),j)=ca_f(at(rs)%bb(1),j) - cdt*cd2(j)
        ca_f(at(rs)%bb(3),j)=ca_f(at(rs)%bb(3),j) - cdt*cd3(j)
        ca_f(at(rs)%sc(1),j)=ca_f(at(rs)%sc(1),j) - cdt*cd4(j)
      end do
      call onetor_cosderiv(at(rs)%bb(4),at(rs)%bb(1),at(rs)%bb(3),&
 &at(rs)%bb(5),cost,cd1,cd2,cd3,cd4,si)
      cdt = -scale_CORR*4.900*2.0*cost
      do j=1,3
        ca_f(at(rs)%bb(4),j)=ca_f(at(rs)%bb(4),j) - cdt*cd1(j)
        ca_f(at(rs)%bb(1),j)=ca_f(at(rs)%bb(1),j) - cdt*cd2(j)
        ca_f(at(rs)%bb(3),j)=ca_f(at(rs)%bb(3),j) - cdt*cd3(j)
        ca_f(at(rs)%bb(5),j)=ca_f(at(rs)%bb(5),j) - cdt*cd4(j)
      end do
    else if (resname.eq.'NMA') then
      call onetor_cosderiv(at(rs)%bb(2),at(rs)%bb(1),at(rs)%bb(3),&
 &at(rs)%sc(1),cost,cd1,cd2,cd3,cd4,si)
      cdt = -scale_CORR*6.089*2.0*cost
      do j=1,3
        ca_f(at(rs)%bb(2),j)=ca_f(at(rs)%bb(2),j) - cdt*cd1(j)
        ca_f(at(rs)%bb(1),j)=ca_f(at(rs)%bb(1),j) - cdt*cd2(j)
        ca_f(at(rs)%bb(3),j)=ca_f(at(rs)%bb(3),j) - cdt*cd3(j)
        ca_f(at(rs)%sc(1),j)=ca_f(at(rs)%sc(1),j) - cdt*cd4(j)
      end do
      call onetor_cosderiv(at(rs)%bb(2),at(rs)%bb(1),at(rs)%bb(3),&
 &at(rs)%bb(4),cost,cd1,cd2,cd3,cd4,si)
      cdt = -scale_CORR*4.900*2.0*cost
      do j=1,3
        ca_f(at(rs)%bb(2),j)=ca_f(at(rs)%bb(2),j) - cdt*cd1(j)
        ca_f(at(rs)%bb(1),j)=ca_f(at(rs)%bb(1),j) - cdt*cd2(j)
        ca_f(at(rs)%bb(3),j)=ca_f(at(rs)%bb(3),j) - cdt*cd3(j)
        ca_f(at(rs)%bb(4),j)=ca_f(at(rs)%bb(4),j) - cdt*cd4(j)
      end do
      call onetor_cosderiv(at(rs)%sc(2),at(rs)%bb(1),at(rs)%bb(3),&
 &at(rs)%sc(1),cost,cd1,cd2,cd3,cd4,si)
      cdt = scale_CORR*((2.300/2.0) - 6.089*2.0*cost)
      do j=1,3
        ca_f(at(rs)%sc(2),j)=ca_f(at(rs)%sc(2),j) - cdt*cd1(j)
        ca_f(at(rs)%bb(1),j)=ca_f(at(rs)%bb(1),j) - cdt*cd2(j)
        ca_f(at(rs)%bb(3),j)=ca_f(at(rs)%bb(3),j) - cdt*cd3(j)
        ca_f(at(rs)%sc(1),j)=ca_f(at(rs)%sc(1),j) - cdt*cd4(j)
      end do
      call onetor_cosderiv(at(rs)%sc(2),at(rs)%bb(1),at(rs)%bb(3),&
 &at(rs)%bb(4),cost,cd1,cd2,cd3,cd4,si)
      cdt = -scale_CORR*4.900*2.0*cost
      do j=1,3
        ca_f(at(rs)%sc(2),j)=ca_f(at(rs)%sc(2),j) - cdt*cd1(j)
        ca_f(at(rs)%bb(1),j)=ca_f(at(rs)%bb(1),j) - cdt*cd2(j)
        ca_f(at(rs)%bb(3),j)=ca_f(at(rs)%bb(3),j) - cdt*cd3(j)
        ca_f(at(rs)%bb(4),j)=ca_f(at(rs)%bb(4),j) - cdt*cd4(j)
      end do
    else 
      call onetor_cosderiv(oi(rs-1),ci(rs-1),ni(rs),cai(rs),&
 &cost,cd1,cd2,cd3,cd4,si)
      cdt = -scale_CORR*6.089*2.0*cost
      do j=1,3
        ca_f(oi(rs-1),j)=ca_f(oi(rs-1),j) - cdt*cd1(j)
        ca_f(ci(rs-1),j)=ca_f(ci(rs-1),j) - cdt*cd2(j)
        ca_f(ni(rs),j)=ca_f(ni(rs),j) - cdt*cd3(j)
        ca_f(cai(rs),j)=ca_f(cai(rs),j) - cdt*cd4(j)
      end do
      call onetor_cosderiv(oi(rs-1),ci(rs-1),ni(rs),hni(rs),&
 &cost,cd1,cd2,cd3,cd4,si)
      cdt = -scale_CORR*4.900*2.0*cost
      do j=1,3
        ca_f(oi(rs-1),j)=ca_f(oi(rs-1),j) - cdt*cd1(j)
        ca_f(ci(rs-1),j)=ca_f(ci(rs-1),j) - cdt*cd2(j)
        ca_f(ni(rs),j)=ca_f(ni(rs),j) - cdt*cd3(j)
        ca_f(hni(rs),j)=ca_f(hni(rs),j) - cdt*cd4(j)
      end do
      if (seqtyp(rs-1).eq.28) then
        ito = at(rs-1)%bb(3)
      else 
        ito = cai(rs-1)
      end if
      if (t2flag.EQV..true.) return
      call onetor_cosderiv(ito,ci(rs-1),ni(rs),cai(rs),&
 &cost,cd1,cd2,cd3,cd4,si)
      cdt = scale_CORR*((2.300/2.0) - 6.089*2.0*cost)
      do j=1,3
        ca_f(ito,j)=ca_f(ito,j) - cdt*cd1(j)
        ca_f(ci(rs-1),j)=ca_f(ci(rs-1),j) - cdt*cd2(j)
        ca_f(ni(rs),j)=ca_f(ni(rs),j) - cdt*cd3(j)
        ca_f(cai(rs),j)=ca_f(cai(rs),j) - cdt*cd4(j)
      end do
      call onetor_cosderiv(ito,ci(rs-1),ni(rs),hni(rs),&
 &cost,cd1,cd2,cd3,cd4,si)
      cdt = -scale_CORR*4.900*2.0*cost
      do j=1,3
        ca_f(ito,j)=ca_f(ito,j) - cdt*cd1(j)
        ca_f(ci(rs-1),j)=ca_f(ci(rs-1),j) - cdt*cd2(j)
        ca_f(ni(rs),j)=ca_f(ni(rs),j) - cdt*cd3(j)
        ca_f(hni(rs),j)=ca_f(hni(rs),j) - cdt*cd4(j)
      end do
    end if
  end if
!
end
!
!-----------------------------------------------------------------------
!
! see ephenol
! note that energy-scaling, but not force-scaling is done on the outside
! 
subroutine get_fphenol(rs,resname,fphenol,ca_f)
!
  use iounit
  use fyoc
  use math
  use forces
  use sequen
  use polypep
  use energies
  use system
  use atoms
  use wl
!
  implicit none
!
  integer rs,imol,j,shf,shf2,shf3
  RTYPE fphenol,t1,t2,getztor,si,ca_f(n,3)
  RTYPE cost,cd1(3),cd2(3),cd3(3),cd4(3),dt0
  character(3) resname
!
  fphenol = 0.0
  imol = molofrs(rs)
  shf = 0
  shf2 = 0
  shf3 = 0
  if (ua_model.gt.0) then
    shf = 1
    if (ua_model.eq.1) then
      shf2 = 3
      shf3 = 3
    else
      shf2 = 3
      shf3 = 7
    end if
  end if
!
  if (resname.eq.'TYR') then
    t1 = chi(3,rs)
    t2=getztor(at(rs)%sc(6-shf),at(rs)%sc(8-shf),at(rs)%sc(9-shf),at(rs)%sc(16-shf3))
    fphenol = 0.5*1.682*((1.0 - cos(2.0*t1/RADIAN)) + &
 &                   (1.0 - cos(2.0*t2/RADIAN)))
    if ((fycxyz.eq.1).AND.(hmjam%isin(2).EQV..false.).AND.(hmjam%isin(11).EQV..false.)) then
      dc_di(imol)%f(chinr(3,rs)) = dc_di(imol)%f(chinr(3,rs)) - &
 &                    scale_CORR*0.5*1.682*&
 &    (sin(2.0*t1/RADIAN)*2.0 + sin(2.0*t2/RADIAN)*2.0) !in kcal/(mol*rad)
!   in Cartesian dynamics we have to go through the tedious process of 
!   correctly mapping the torsional derivative to the xyz's of the four atoms
    else
      call onetor_cosderiv(at(rs)%sc(7-shf),at(rs)%sc(8-shf),at(rs)%sc(9-shf),&
 &at(rs)%sc(16-shf3),cost,cd1,cd2,cd3,cd4,si)
      dt0 = -scale_CORR*1.682*2.0*cost
      do j=1,3
        ca_f(at(rs)%sc(7-shf),j)=ca_f(at(rs)%sc(7-shf),j) - dt0*cd1(j)
        ca_f(at(rs)%sc(8-shf),j)=ca_f(at(rs)%sc(8-shf),j) - dt0*cd2(j)
        ca_f(at(rs)%sc(9-shf),j)=ca_f(at(rs)%sc(9-shf),j) - dt0*cd3(j)
        ca_f(at(rs)%sc(16-shf3),j)=ca_f(at(rs)%sc(16-shf3),j) - dt0*cd4(j)
      end do
      call onetor_cosderiv(at(rs)%sc(6-shf),at(rs)%sc(8-shf),at(rs)%sc(9-shf),&
 &at(rs)%sc(16-shf3),cost,cd1,cd2,cd3,cd4,si)
      dt0 = -scale_CORR*1.682*2.0*cost
      do j=1,3
        ca_f(at(rs)%sc(6-shf),j)=ca_f(at(rs)%sc(6-shf),j) - dt0*cd1(j)
        ca_f(at(rs)%sc(8-shf),j)=ca_f(at(rs)%sc(8-shf),j) - dt0*cd2(j)
        ca_f(at(rs)%sc(9-shf),j)=ca_f(at(rs)%sc(9-shf),j) - dt0*cd3(j)
        ca_f(at(rs)%sc(16-shf3),j)=ca_f(at(rs)%sc(16-shf3),j) - dt0*cd4(j)
      end do
    end if
  else if (resname.eq.'PCR') then
    t1 = chi(2-shf,rs)
    t2=getztor(at(rs)%bb(6),at(rs)%bb(7),at(rs)%bb(8),at(rs)%sc(4-shf2))
    fphenol = 0.5*1.682*((1.0 - cos(2.0*t1/RADIAN)) + &
 &                   (1.0 - cos(2.0*t2/RADIAN)))
    if ((fycxyz.eq.1).AND.(hmjam%isin(2).EQV..false.).AND.(hmjam%isin(11).EQV..false.)) then
      dc_di(imol)%f(chinr(2-shf,rs)) = dc_di(imol)%f(chinr(2-shf,rs)) - &
 &                    scale_CORR*0.5*1.682*&
 &    (sin(2.0*t1/RADIAN)*2.0 + sin(2.0*t2/RADIAN)*2.0) !in kcal/(mol*rad)
    else
      call onetor_cosderiv(at(rs)%bb(5),at(rs)%bb(7),at(rs)%bb(8),&
 &at(rs)%sc(4-shf2),cost,cd1,cd2,cd3,cd4,si)
      dt0 = -scale_CORR*1.682*2.0*cost
      do j=1,3
        ca_f(at(rs)%bb(5),j)=ca_f(at(rs)%bb(5),j) - dt0*cd1(j)
        ca_f(at(rs)%bb(7),j)=ca_f(at(rs)%bb(7),j) - dt0*cd2(j)
        ca_f(at(rs)%bb(8),j)=ca_f(at(rs)%bb(8),j) - dt0*cd3(j)
        ca_f(at(rs)%sc(4-shf2),j)=ca_f(at(rs)%sc(4-shf2),j) - dt0*cd4(j)
      end do
      call onetor_cosderiv(at(rs)%bb(6),at(rs)%bb(7),at(rs)%bb(8),&
 &at(rs)%sc(4-shf2),cost,cd1,cd2,cd3,cd4,si)
      dt0 = -scale_CORR*1.682*2.0*cost
      do j=1,3
        ca_f(at(rs)%bb(6),j)=ca_f(at(rs)%bb(6),j) - dt0*cd1(j)
        ca_f(at(rs)%bb(7),j)=ca_f(at(rs)%bb(7),j) - dt0*cd2(j)
        ca_f(at(rs)%bb(8),j)=ca_f(at(rs)%bb(8),j) - dt0*cd3(j)
        ca_f(at(rs)%sc(4-shf2),j)=ca_f(at(rs)%sc(4-shf2),j) - dt0*cd4(j)
      end do
    end if
  else
    write(ilog,*) 'Fatal. Called get_fphenol(...) with an inapplicable r&
 &esidue (',rs,': ',resname,').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
! these are routines meant to provide forces from bonded terms to be
! used in Cartesian dynamics (in this case bond length restraints)
! note that in all of these we're expanding the function computing
! the internal into this function to make it more efficient
!
subroutine force_bonds(rs,evec,ca_f)
!
  use iounit
  use params
  use energies
  use polypep
  use inter
  use atoms
  use forces
  use system
  use sequen
  use fyoc
!
  implicit none
!
  integer ibl1,ibl2,i,j,rs
  RTYPE dvec(nrsbleff(rs),3),d2(nrsbleff(rs)),d1(nrsbleff(rs))
  RTYPE globscal,svec(3),svecref(3)
  RTYPE id1(nrsbleff(rs)),incr,term0,term1,evec(MAXENERGYTERMS),ca_f(n,3)
!
  if (nrsbleff(rs).eq.0) return
!
  globscal = scale_BOND(1)
  if ((use_FEG.EQV..true.).AND.(fegmode.eq.2)) then
    if (par_FEG(rs).EQV..true.) then
      if (use_FEGS(15).EQV..false.) then
        return
      else
        globscal = scale_FEGS(15)
      end if
    end if
  end if
!
  if (disulf(rs).gt.0) then
    if (molofrs(disulf(rs)).ne.molofrs(rs)) then
      call dis_bound_rs(rs,disulf(rs),svecref)
    else
      svecref(:) = 0.0
    end if
  else
    svecref(:) = 0.0
  end if
  do i=1,nrsbleff(rs)
    ibl1 = iaa(rs)%bl(i,1)
    ibl2 = iaa(rs)%bl(i,2)
    svec(:) = 0.0
    if (atmres(ibl2).eq.disulf(rs)) svec(:) = svecref(:)
    dvec(i,1) = x(ibl2) - x(ibl1) + svec(1)
    dvec(i,2) = y(ibl2) - y(ibl1) + svec(2)
    dvec(i,3) = z(ibl2) - z(ibl1) + svec(3)
  end do
  d2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
  d1(:) = sqrt(d2(:))
  id1(:) = 1.0/d1(:)
!
  do i=1,nrsbleff(rs)
    ibl1 = iaa(rs)%bl(i,1)
    ibl2 = iaa(rs)%bl(i,2)
!   harmonic
    if (iaa(rs)%typ_bl(i).eq.1) then
      term0 = globscal*iaa(rs)%par_bl(i,1)
      evec(15) = evec(15)+ term0*((d1(i)-iaa(rs)%par_bl(i,2))**2)
      term0  = 2.0*term0*(d1(i)-iaa(rs)%par_bl(i,2))*id1(i)
      do j=1,3
        incr = term0*dvec(i,j)
        ca_f(ibl1,j) = ca_f(ibl1,j) + incr
        ca_f(ibl2,j) = ca_f(ibl2,j) - incr
!        if (pflag.EQV..true.) then
!          ens%insVirT(j,j) = ens%insVirT(j,j) - incr*dvec(i,j)
!        end if
      end do
!   Morse (anharmonic)
    else if (iaa(rs)%typ_bl(i).eq.2) then
      term0 = globscal*iaa(rs)%par_bl(i,3)
      term1 = exp(-iaa(rs)%par_bl(i,1)*(d1(i)-iaa(rs)%par_bl(i,2)))
      evec(15) = evec(15) + term0*((1.0 - term1)**2)
      term0 = term0*2.0*iaa(rs)%par_bl(i,1)*&
 &                               (1.0 - term1)*term1*id1(i)
      do j=1,3
        incr = term0*dvec(i,j)
        ca_f(ibl1,j) = ca_f(ibl1,j) + incr
        ca_f(ibl2,j) = ca_f(ibl2,j) - incr
!        if (pflag.EQV..true.) then
!          ens%insVirT(j,j) = ens%insVirT(j,j) - incr*dvec(i,j)
!        end if
      end do
!   GROMOS format (quartic)
    else if  (iaa(rs)%typ_bl(i).eq.3) then
      term0 = globscal*iaa(rs)%par_bl(i,1)
      evec(15) = evec(15)+ 0.25*term0*((d2(i)-iaa(rs)%par_bl(i,3))**2)
      term0  = term0*(d2(i)-iaa(rs)%par_bl(i,3))
      do j=1,3
        incr = term0*dvec(i,j)
        ca_f(ibl1,j) = ca_f(ibl1,j) + incr
        ca_f(ibl2,j) = ca_f(ibl2,j) - incr
      end do
    else
      write(ilog,*) 'Fatal. Encountered unsupported type of bond length&
 & potential in force_bonds(...) (offending type is ',&
 &iaa(rs)%typ_bl(i),').'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! computes dd/dx,y,z for a bond - does not support crosslinks
!
subroutine onebond_deriv(ibl1,ibl2,fvec1,fvec2)
!
  use atoms
!
  implicit none
!
  integer ibl1,ibl2
  RTYPE fvec1(3),fvec2(3)
  RTYPE dvec(3),id1,d2
!
  dvec(1) = x(ibl2) - x(ibl1)
  dvec(2) = y(ibl2) - y(ibl1)
  dvec(3) = z(ibl2) - z(ibl1)
  d2 = dvec(1)**2 + dvec(2)**2 + dvec(3)**2
  id1 = 1.0/sqrt(d2)
  fvec2(:) = id1*dvec(:)
  fvec1(:) = -fvec2(:)
!
end
!
!-------------------------------------------------------------------------------------------------
!
! note that for angles and higher terms the pairwise symmetry of 
! forces disappears of course (in favor of higher-order symmetry) ...
!
subroutine force_angles(rs,evec,ca_f)
!
  use iounit
  use params
  use energies
  use polypep
  use inter
  use atoms
  use forces
  use math
  use fyoc
  use sequen
!
  implicit none
!
  integer iba1,iba2,iba3,i,j,rs
  RTYPE dv12(nrsbaeff(rs),3),id12(nrsbaeff(rs))
  RTYPE dv32(nrsbaeff(rs),3),id32(nrsbaeff(rs)),dotp(nrsbaeff(rs))
  RTYPE rang(nrsbaeff(rs)),svec(3),svecref(3)
  RTYPE vecs(nrsbaeff(rs),3),nvecs(nrsbaeff(rs))
  RTYPE frang(nrsbaeff(rs)),dv31(nrsbaeff(rs),3)
  RTYPE fvec1(nrsbaeff(rs),3),fvec3(nrsbaeff(rs),3)
  RTYPE fvec2(nrsbaeff(rs),3),term2,si(nrsbaeff(rs))
  RTYPE incr,term0,evec(MAXENERGYTERMS),d2,d1,id1,globscal,ca_f(n,3),eps
!
  if (nrsbaeff(rs).eq.0) return
!
  eps = 1.0*(10.0**(-precision(eps)))
  globscal = scale_BOND(2)
  if ((use_FEG.EQV..true.).AND.(fegmode.eq.2)) then
    if (par_FEG(rs).EQV..true.) then
      if (use_FEGS(16).EQV..false.) then
        return
      else
        globscal = scale_FEGS(16)
      end if
    end if
  end if
!
  if (disulf(rs).gt.0) then
    if (molofrs(disulf(rs)).ne.molofrs(rs)) then
      call dis_bound_rs(rs,disulf(rs),svecref)
    else
      svecref(:) = 0.0
    end if
  else
    svecref(:) = 0.0
  end if
  do i=1,nrsbaeff(rs)
    iba1 = iaa(rs)%ba(i,1)
    iba2 = iaa(rs)%ba(i,2)
    iba3 = iaa(rs)%ba(i,3)
    svec(:) = 0.0
    if (atmres(iba2).eq.disulf(rs)) svec(:) = -svecref(:)
    dv12(i,1) = x(iba1) - x(iba2) + svec(1)
    dv12(i,2) = y(iba1) - y(iba2) + svec(2)
    dv12(i,3) = z(iba1) - z(iba2) + svec(3)
    svec(:) = 0.0
    if ((atmres(iba3).eq.disulf(rs)).AND.(atmres(iba2).ne.disulf(rs))) svec(:) = svecref(:)
    dv32(i,1) = x(iba3) - x(iba2) + svec(1)
    dv32(i,2) = y(iba3) - y(iba2) + svec(2)
    dv32(i,3) = z(iba3) - z(iba2) + svec(3)
!   add the 13-vector for Urey-Bradley
    if (iaa(rs)%typ_ba(i).eq.2) then
      svec(:) = 0.0
      if (atmres(iba3).eq.disulf(rs)) svec(:) = svecref(:)
      dv31(i,1) = x(iba3) - x(iba1) + svec(1)
      dv31(i,2) = y(iba3) - y(iba1) + svec(2)
      dv31(i,3) = z(iba3) - z(iba1) + svec(3)
    end if
  end do
  id12(:) = 1.0/sqrt(dv12(:,1)**2 + dv12(:,2)**2 + dv12(:,3)**2)
  id32(:) = 1.0/sqrt(dv32(:,1)**2 + dv32(:,2)**2 + dv32(:,3)**2)
  dotp(:) = dv12(:,1)*dv32(:,1) + dv12(:,2)*dv32(:,2) + &
 &          dv12(:,3)*dv32(:,3)
  do i=1,nrsbaeff(rs)
    if (dotp(i).lt.0.0) then
      si(i) = -1.0
    else
      si(i) = 1.0
    end if
  end do
! this half-angle formula has the advantage of being safe to use in all circumstances
! as 0.5*nvecs never gets near the -1/1 limits
! in addition there is considerable precision loss with the (simpler) acos-form
  vecs(:,1) = id32(:)*dv32(:,1) - si(:)*id12(:)*dv12(:,1)
  vecs(:,2) = id32(:)*dv32(:,2) - si(:)*id12(:)*dv12(:,2)
  vecs(:,3) = id32(:)*dv32(:,3) - si(:)*id12(:)*dv12(:,3)
  nvecs(:) = sqrt(vecs(:,1)**2 + vecs(:,2)**2 + vecs(:,3)**2)
  rang(:) = si(:)*2.0*RADIAN*asin(0.5*nvecs(:))
  frang(:) = si(:)*2.0*RADIAN/sqrt(1.0 - (0.5*nvecs(:))**2)
  do i=1,nrsbaeff(rs)
    if (dotp(i).lt.0.0) then
      rang(i) = rang(i) + RADIAN*PI
    end if
   nvecs(i) = max(nvecs(i),eps)
  end do
! there is a price, however, primarily in the derivatives
  fvec1(:,1) = 0.5*frang(:)*(1.0/nvecs(:))*(vecs(:,1)*&
 &  (-si(:)*id12(:) + si(:)*dv12(:,1)*(id12(:)**3)*dv12(:,1))&
 & + vecs(:,2)*(si(:)*dv12(:,2)*(id12(:)**3)*dv12(:,1))&
 & + vecs(:,3)*(si(:)*dv12(:,3)*(id12(:)**3)*dv12(:,1)) )
  fvec1(:,2) = 0.5*frang(:)*(1.0/nvecs(:))*(vecs(:,2)*&
 &  (-si(:)*id12(:) + si(:)*dv12(:,2)*(id12(:)**3)*dv12(:,2))&
 & + vecs(:,1)*(si(:)*dv12(:,1)*(id12(:)**3)*dv12(:,2))&
 & + vecs(:,3)*(si(:)*dv12(:,3)*(id12(:)**3)*dv12(:,2)) )
  fvec1(:,3) = 0.5*frang(:)*(1.0/nvecs(:))*(vecs(:,3)*&
 &  (-si(:)*id12(:) + si(:)*dv12(:,3)*(id12(:)**3)*dv12(:,3))&
 & + vecs(:,1)*(si(:)*dv12(:,1)*(id12(:)**3)*dv12(:,3))&
 & + vecs(:,2)*(si(:)*dv12(:,2)*(id12(:)**3)*dv12(:,3)) )
  fvec3(:,1) = 0.5*frang(:)*(1.0/nvecs(:))*(vecs(:,1)*&
 &  (id32(:) - dv32(:,1)*(id32(:)**3)*dv32(:,1))&
 & + vecs(:,2)*(-dv32(:,2)*(id32(:)**3)*dv32(:,1))&
 & + vecs(:,3)*(-dv32(:,3)*(id32(:)**3)*dv32(:,1)) )
  fvec3(:,2) = 0.5*frang(:)*(1.0/nvecs(:))*(vecs(:,2)*&
 &  (id32(:) - dv32(:,2)*(id32(:)**3)*dv32(:,2))&
 & + vecs(:,1)*(-dv32(:,1)*(id32(:)**3)*dv32(:,2))&
 & + vecs(:,3)*(-dv32(:,3)*(id32(:)**3)*dv32(:,2)) )
  fvec3(:,3) = 0.5*frang(:)*(1.0/nvecs(:))*(vecs(:,3)*&
 &  (id32(:) - dv32(:,3)*(id32(:)**3)*dv32(:,3))&
 & + vecs(:,1)*(-dv32(:,1)*(id32(:)**3)*dv32(:,3))&
 & + vecs(:,2)*(-dv32(:,2)*(id32(:)**3)*dv32(:,3)) )
  do j=1,3
    fvec2(:,j) = - (fvec1(:,j) + fvec3(:,j))
  end do
!
  do i=1,nrsbaeff(rs)
    iba1 = iaa(rs)%ba(i,1)
    iba2 = iaa(rs)%ba(i,2)
    iba3 = iaa(rs)%ba(i,3)
!   harmonic
    if (iaa(rs)%typ_ba(i).eq.1) then
      term0 = globscal*iaa(rs)%par_ba(i,1)
      evec(16) = evec(16) + term0*&
 &             ((rang(i)-iaa(rs)%par_ba(i,2))**2)
      term0  = 2.0*term0*(rang(i)-iaa(rs)%par_ba(i,2))
      do j=1,3
        ca_f(iba1,j) = ca_f(iba1,j) - term0*fvec1(i,j)
        ca_f(iba2,j) = ca_f(iba2,j) - term0*fvec2(i,j)
        ca_f(iba3,j) = ca_f(iba3,j) - term0*fvec3(i,j)
      end do
!   UB
    else if (iaa(rs)%typ_ba(i).eq.2) then
!     get the extra pseudo-bond term for Urey-Bradley
      d2 = dv31(i,1)**2 + dv31(i,2)**2 + dv31(i,3)**2
      d1 = sqrt(d2)
      id1 = 1.0/d1
      term0 = globscal*iaa(rs)%par_ba(i,1)
      term2 = globscal*iaa(rs)%par_ba(i,3)
      evec(16) = evec(16) + term0*&
 &             ((rang(i)-iaa(rs)%par_ba(i,2))**2) + &
 &         term2*((d1-iaa(rs)%par_ba(i,4))**2)
      term0 = 2.0*term0*(rang(i)-iaa(rs)%par_ba(i,2))
      term2 = 2.0*term2*(d1-iaa(rs)%par_ba(i,4))*id1
      do j=1,3
        incr = term2*dv31(i,j)
        ca_f(iba1,j) = ca_f(iba1,j) - term0*fvec1(i,j) + incr
        ca_f(iba2,j) = ca_f(iba2,j) - term0*fvec2(i,j)
        ca_f(iba3,j) = ca_f(iba3,j) - term0*fvec3(i,j) - incr
      end do
!   GROMOS harmonic cosine (note we're wasting the computational advantage of the cos-form here)
    else if (iaa(rs)%typ_ba(i).eq.3) then
      term0 = 0.5*globscal*iaa(rs)%par_ba(i,1)
      term2 = dotp(i)*(id12(i)*id32(i))
      evec(16) = evec(16) + term0*((term2 - iaa(rs)%par_ba(i,3))**2)
      term0 = -(2.0/RADIAN)*term0*(term2 - iaa(rs)%par_ba(i,3))*sin(rang(i)/RADIAN)
      do j=1,3
        ca_f(iba1,j) = ca_f(iba1,j) - term0*fvec1(i,j)
        ca_f(iba2,j) = ca_f(iba2,j) - term0*fvec2(i,j)
        ca_f(iba3,j) = ca_f(iba3,j) - term0*fvec3(i,j)
      end do
    else
      write(ilog,*) 'Fatal. Encountered unsupported type of bond angle &
 &potential in force_angles(...) (offending type is ',&
 &iaa(rs)%typ_ba(i),').'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! a routine to just get Cartesian derivatives dalpha/dx,y,z for a given bond angle
! does not work with crosslinks
!
subroutine oneang_deriv(iba1,iba2,iba3,rang,fvec1,fvec2,fvec3)
!
  use atoms
  use math
!
  implicit none
!
  integer iba1,iba2,iba3
  RTYPE dotp,si,dv12(3),dv32(3),id12,id32,vecs(3),fvec1(3),fvec2(3),fvec3(3),nvecs,rang,frang
!
  dv12(1) = x(iba1) - x(iba2)
  dv12(2) = y(iba1) - y(iba2)
  dv12(3) = z(iba1) - z(iba2)
  dv32(1) = x(iba3) - x(iba2)
  dv32(2) = y(iba3) - y(iba2)
  dv32(3) = z(iba3) - z(iba2)

  id12 = 1.0/sqrt(dv12(1)**2 + dv12(2)**2 + dv12(3)**2)
  id32 = 1.0/sqrt(dv32(1)**2 + dv32(2)**2 + dv32(3)**2)
  dotp = dv12(1)*dv32(1) + dv12(2)*dv32(2) + &
 &          dv12(3)*dv32(3)
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
! this half-angle formula has the advantage of being safe to use in all circumstances
! as 0.5*nvecs never gets near the -1/1 limits
! in addition there is considerable precision loss with the (simpler) acos-form
  vecs(1) = id32*dv32(1) - si*id12*dv12(1)
  vecs(2) = id32*dv32(2) - si*id12*dv12(2)
  vecs(3) = id32*dv32(3) - si*id12*dv12(3)
  nvecs = sqrt(vecs(1)**2 + vecs(2)**2 + vecs(3)**2)
  rang = si*2.0*RADIAN*asin(0.5*nvecs)
  frang = si*2.0*RADIAN/sqrt(1.0 - (0.5*nvecs)**2)
  if (dotp.lt.0.0) then
    rang = rang + RADIAN*PI
  end if
! there is a price, however, primarily in the derivatives
  fvec1(1) = 0.5*frang*(1.0/nvecs)*(vecs(1)*(-si*id12 + si*dv12(1)*(id12**3)*dv12(1))&
 & + vecs(2)*(si*dv12(2)*(id12**3)*dv12(1)) + vecs(3)*(si*dv12(3)*(id12**3)*dv12(1)) )
  fvec1(2) = 0.5*frang*(1.0/nvecs)*(vecs(2)*(-si*id12 + si*dv12(2)*(id12**3)*dv12(2))&
 & + vecs(1)*(si*dv12(1)*(id12**3)*dv12(2)) + vecs(3)*(si*dv12(3)*(id12**3)*dv12(2)) )
  fvec1(3) = 0.5*frang*(1.0/nvecs)*(vecs(3)*(-si*id12 + si*dv12(3)*(id12**3)*dv12(3))&
 & + vecs(1)*(si*dv12(1)*(id12**3)*dv12(3)) + vecs(2)*(si*dv12(2)*(id12**3)*dv12(3)) )
  fvec3(1) = 0.5*frang*(1.0/nvecs)*(vecs(1)*(id32 - dv32(1)*(id32**3)*dv32(1))&
 & + vecs(2)*(-dv32(2)*(id32**3)*dv32(1)) + vecs(3)*(-dv32(3)*(id32**3)*dv32(1)) )
  fvec3(2) = 0.5*frang*(1.0/nvecs)*(vecs(2)*(id32 - dv32(2)*(id32**3)*dv32(2))&
 & + vecs(1)*(-dv32(1)*(id32**3)*dv32(2)) + vecs(3)*(-dv32(3)*(id32**3)*dv32(2)) )
  fvec3(3) = 0.5*frang*(1.0/nvecs)*(vecs(3)*(id32 - dv32(3)*(id32**3)*dv32(3))&
 & + vecs(1)*(-dv32(1)*(id32**3)*dv32(3)) + vecs(2)*(-dv32(2)*(id32**3)*dv32(3)) )
  fvec2(:) = - (fvec1(:) + fvec3(:))
!
end
!
!-----------------------------------------------------------------------
!
subroutine numangtest(iba1,iba2,iba3,fv1,fv2,fv3)
!
  use atoms
!
  implicit none
!
  integer i,ttt(3),iba1,iba2,iba3
  RTYPE bla,getbang,blas(10),fv1(3),fv2(3),fv3(3)
  RTYPE tvv(3,3)
!
  bla = 0.0001
!
  ttt(1) = iba1
  ttt(2) = iba2
  ttt(3) = iba3
  tvv(1,:) = fv1(:)
  tvv(2,:) = fv2(:)
  tvv(3,:) = fv3(:)
  write(*,*) 'atoms: ',ttt(:)
  write(*,*) 'start: ',getbang(iba1,iba2,iba3)
  do i=1,3
    x(ttt(i)) = x(ttt(i)) + bla
    blas(2) = getbang(iba1,iba2,iba3)
    x(ttt(i)) = x(ttt(i)) - 2.0*bla
    blas(3) = getbang(iba1,iba2,iba3)
    x(ttt(i)) = x(ttt(i)) + bla
    write(*,*) tvv(i,1),(blas(3)-blas(2))/(2.0*bla)
    y(ttt(i)) = y(ttt(i)) + bla
    blas(2) = getbang(iba1,iba2,iba3)
    y(ttt(i)) = y(ttt(i)) - 2.0*bla
    blas(3) = getbang(iba1,iba2,iba3)
    y(ttt(i)) = y(ttt(i)) + bla
    write(*,*) tvv(i,2),(blas(3)-blas(2))/(2.0*bla)
    z(ttt(i)) = z(ttt(i)) + bla
    blas(2) = getbang(iba1,iba2,iba3)
    z(ttt(i)) = z(ttt(i)) - 2.0*bla
    blas(3) = getbang(iba1,iba2,iba3)
    z(ttt(i)) = z(ttt(i)) + bla
    write(*,*) tvv(i,3),(blas(3)-blas(2))/(2.0*bla)
  end do
  write(*,*) 'end: ',getbang(iba1,iba2,iba3)
!
end
!
!-----------------------------------------------------------------------
!
! note that for angles and higher terms the pairwise symmetry of 
! forces disappears of course ...
! the periodicity of dihedrals poses a particular difficulty, since
! the evaluation of the angle between the two planes becomes inaccurate
! around the break point which causes the gradient to diverge for
! non-periodic potentials (like harmonic)
! preferably use the Blondel version below
!
subroutine force_torsions2(rs,evec,ca_f)
!
  use iounit
  use params
  use energies
  use polypep
  use inter
  use atoms
  use forces
  use math
  use fyoc
  use sequen
!
  implicit none
!
  integer ito1,ito2,ito3,ito4,i,j,rs
  RTYPE dv21(nrsdieff(rs),3),dv32(nrsdieff(rs),3)
  RTYPE dv43(nrsdieff(rs),3),ca_f(n,3)
  RTYPE ncp123(nrsdieff(rs)),ncp234(nrsdieff(rs)),dotp(nrsdieff(rs))
  RTYPE termc(nrsdieff(rs)),rang,svec(3),svecref(3)
  RTYPE frang,inpm(nrsdieff(rs)),inpm2(nrsdieff(rs))
  RTYPE dotpd1(nrsdieff(rs),3)
  RTYPE dotpd2(nrsdieff(rs),3)
  RTYPE dotpd3(nrsdieff(rs),3),dotpd4(nrsdieff(rs),3)
  RTYPE cp123(nrsdieff(rs),3),cp234(nrsdieff(rs),3)
  RTYPE ncp123d1(nrsdieff(rs),3),ncp123d2(nrsdieff(rs),3)
  RTYPE ncp123d3(nrsdieff(rs),3),ncp123d4(nrsdieff(rs),3)
  RTYPE ncp234d1(nrsdieff(rs),3),ncp234d2(nrsdieff(rs),3)
  RTYPE ncp234d3(nrsdieff(rs),3),ncp234d4(nrsdieff(rs),3)
  RTYPE dtmd1(nrsdieff(rs),3),dtmd2(nrsdieff(rs),3)
  RTYPE dtmd3(nrsdieff(rs),3),globscal
  RTYPE dtmd4(nrsdieff(rs),3),si,vv,ctp2,ctp3,ctp4,ctp5,ctp6
  RTYPE term0,evec(MAXENERGYTERMS)
!
  if (nrsdieff(rs).eq.0) return
!
  globscal = scale_BOND(4)
  if ((use_FEG.EQV..true.).AND.(fegmode.eq.2)) then
    if (par_FEG(rs).EQV..true.) then
      if (use_FEGS(18).EQV..false.) then
        return
      else
        globscal = scale_FEGS(18)
      end if
    end if
  end if
!
  if (disulf(rs).gt.0) then
    if (molofrs(disulf(rs)).ne.molofrs(rs)) then
      call dis_bound_rs(rs,disulf(rs),svecref)
    else
      svecref(:) = 0.0
    end if
  else
    svecref(:) = 0.0
  end if
  do i=1,nrsdieff(rs)
    ito1 = iaa(rs)%di(i,1)
    ito2 = iaa(rs)%di(i,2)
    ito3 = iaa(rs)%di(i,3)
    ito4 = iaa(rs)%di(i,4)
    svec(:) = 0.0
    if (atmres(ito2).eq.disulf(rs)) svec(:) = svecref(:)
    dv21(i,1) = x(ito2) - x(ito1) + svec(1)
    dv21(i,2) = y(ito2) - y(ito1) + svec(2)
    dv21(i,3) = z(ito2) - z(ito1) + svec(3)
    svec(:) = 0.0
    if ((atmres(ito3).eq.disulf(rs)).AND.(atmres(ito2).ne.disulf(rs))) svec(:) = svecref(:)
    if ((atmres(ito3).ne.disulf(rs)).AND.(atmres(ito2).eq.disulf(rs))) svec(:) = -svecref(:)
    dv32(i,1) = x(ito3) - x(ito2) + svec(1)
    dv32(i,2) = y(ito3) - y(ito2) + svec(2)
    dv32(i,3) = z(ito3) - z(ito2) + svec(3)
    svec(:) = 0.0
    if ((atmres(ito4).eq.disulf(rs)).AND.(atmres(ito3).ne.disulf(rs))) svec(:) = svecref(:)
    if ((atmres(ito4).ne.disulf(rs)).AND.(atmres(ito3).eq.disulf(rs))) svec(:) = -svecref(:)
    dv43(i,1) = x(ito4) - x(ito3) + svec(1)
    dv43(i,2) = y(ito4) - y(ito3) + svec(2)
    dv43(i,3) = z(ito4) - z(ito3) + svec(3)
  end do
  cp123(:,1) = dv21(:,2)*dv32(:,3) - dv21(:,3)*dv32(:,2)
  cp123(:,2) = dv21(:,3)*dv32(:,1) - dv21(:,1)*dv32(:,3)
  cp123(:,3) = dv21(:,1)*dv32(:,2) - dv21(:,2)*dv32(:,1)
  cp234(:,1) = dv32(:,2)*dv43(:,3) - dv32(:,3)*dv43(:,2)
  cp234(:,2) = dv32(:,3)*dv43(:,1) - dv32(:,1)*dv43(:,3)
  cp234(:,3) = dv32(:,1)*dv43(:,2) - dv32(:,2)*dv43(:,1)
! norms and dot-product
  ncp123(:) = cp123(:,1)**2 + cp123(:,2)**2 + cp123(:,3)**2
  ncp234(:) = cp234(:,1)**2 + cp234(:,2)**2 + cp234(:,3)**2
  dotp(:) = cp123(:,1)*cp234(:,1) + cp123(:,2)*cp234(:,2) + &
 &          cp123(:,3)*cp234(:,3)
! dot product derivatives
  dotpd1(:,1) = dv32(:,3)*cp234(:,2) - dv32(:,2)*cp234(:,3)
  dotpd1(:,2) = -dv32(:,3)*cp234(:,1) + dv32(:,1)*cp234(:,3)
  dotpd1(:,3) = dv32(:,2)*cp234(:,1) - dv32(:,1)*cp234(:,2)
  dotpd2(:,1) =  (-dv21(:,3)-dv32(:,3))*cp234(:,2) +&
 &                cp123(:,2)*dv43(:,3) + &
 &               (dv32(:,2)+dv21(:,2))*cp234(:,3) -&
 &                cp123(:,3)*dv43(:,2)
  dotpd2(:,2) =  (dv32(:,3)+dv21(:,3))*cp234(:,1) -&
 &                cp123(:,1)*dv43(:,3) + &
 &               (-dv21(:,1)-dv32(:,1))*cp234(:,3) +&
 &                cp123(:,3)*dv43(:,1)
  dotpd2(:,3) =  (-dv21(:,2)-dv32(:,2))*cp234(:,1) +&
 &                cp123(:,1)*dv43(:,2) + &
 &               (dv32(:,1)+dv21(:,1))*cp234(:,2) -&
 &                cp123(:,2)*dv43(:,1)
  dotpd3(:,1) =  (-dv32(:,3)-dv43(:,3))*cp123(:,2) +&
 &                cp234(:,2)*dv21(:,3) + &
 &               (dv43(:,2)+dv32(:,2))*cp123(:,3) -&
 &                cp234(:,3)*dv21(:,2)
  dotpd3(:,2) =  (dv43(:,3)+dv32(:,3))*cp123(:,1) -&
 &                cp234(:,1)*dv21(:,3) + &
 &               (-dv32(:,1)-dv43(:,1))*cp123(:,3) +&
 &                cp234(:,3)*dv21(:,1)
  dotpd3(:,3) =  (-dv32(:,2)-dv43(:,2))*cp123(:,1) +&
 &                cp234(:,1)*dv21(:,2) + &
 &               (dv43(:,1)+dv32(:,1))*cp123(:,2) -&
 &                cp234(:,2)*dv21(:,1)
  dotpd4(:,1) = dv32(:,3)*cp123(:,2) - dv32(:,2)*cp123(:,3)
  dotpd4(:,2) = -dv32(:,3)*cp123(:,1) + dv32(:,1)*cp123(:,3)
  dotpd4(:,3) = dv32(:,2)*cp123(:,1) - dv32(:,1)*cp123(:,2) 
! quadratic norm derivatives
  ncp123d1(:,1) = 2.0*(cp123(:,2)*dv32(:,3) - cp123(:,3)*dv32(:,2))
  ncp123d1(:,2) = 2.0*(-cp123(:,1)*dv32(:,3) + cp123(:,3)*dv32(:,1))
  ncp123d1(:,3) = 2.0*(cp123(:,1)*dv32(:,2) - cp123(:,2)*dv32(:,1))
  ncp234d1(:,:) = 0.0
  ncp123d2(:,1) = 2.0*(cp123(:,2)*(-dv21(:,3)-dv32(:,3)) + &
 &                     cp123(:,3)*(dv32(:,2)+dv21(:,2)))
  ncp123d2(:,2) = 2.0*(cp123(:,1)*(dv32(:,3)+dv21(:,3)) + &
 &                     cp123(:,3)*(-dv21(:,1)-dv32(:,1)))
  ncp123d2(:,3) = 2.0*(cp123(:,1)*(-dv21(:,2)-dv32(:,2)) + &
 &                     cp123(:,2)*(dv32(:,1)+dv21(:,1)))
  ncp234d2(:,1) = 2.0*(cp234(:,2)*dv43(:,3) - cp234(:,3)*dv43(:,2))
  ncp234d2(:,2) = 2.0*(-cp234(:,1)*dv43(:,3) + cp234(:,3)*dv43(:,1))
  ncp234d2(:,3) = 2.0*(cp234(:,1)*dv43(:,2) - cp234(:,2)*dv43(:,1))
  ncp123d3(:,1) = 2.0*(cp123(:,2)*dv21(:,3) - cp123(:,3)*dv21(:,2))
  ncp123d3(:,2) = 2.0*(-cp123(:,1)*dv21(:,3) + cp123(:,3)*dv21(:,1))
  ncp123d3(:,3) = 2.0*(cp123(:,1)*dv21(:,2) - cp123(:,2)*dv21(:,1))
  ncp234d3(:,1) = 2.0*(cp234(:,2)*(-dv32(:,3)-dv43(:,3)) + &
 &                     cp234(:,3)*(dv43(:,2)+dv32(:,2)))
  ncp234d3(:,2) = 2.0*(cp234(:,1)*(dv43(:,3)+dv32(:,3)) + &
 &                     cp234(:,3)*(-dv32(:,1)-dv43(:,1)))
  ncp234d3(:,3) = 2.0*(cp234(:,1)*(-dv32(:,2)-dv43(:,2)) + &
 &                     cp234(:,2)*(dv43(:,1)+dv32(:,1)))
  ncp123d4(:,:) = 0.0
  ncp234d4(:,1) = 2.0*(cp234(:,2)*dv32(:,3) - cp234(:,3)*dv32(:,2))
  ncp234d4(:,2) = 2.0*(-cp234(:,1)*dv32(:,3) + cp234(:,3)*dv32(:,1))
  ncp234d4(:,3) = 2.0*(cp234(:,1)*dv32(:,2) - cp234(:,2)*dv32(:,1))
! filter for exceptions which would otherwise NaN
  do i=1,nrsdieff(rs)
!   colinear reference atoms
    if ((ncp123(i).le.0.0).OR.(ncp234(i).le.0.0)) then
!     this will set all derivatives to zero and the (now arbitrary) dihedral to zero deg.
      fo_wrncnt(4) = fo_wrncnt(4) + 1
      if (fo_wrncnt(4).eq.fo_wrnlmt(4)) then
        write(ilog,*) 'WARNING. Colinear reference atoms in computation o&
 &f torsional forces (acos). This indicates bad parameters, an unstable simulation, or a bug&
 &. Expect run to crash soon.'
        write(ilog,*) 'This was warning #',fo_wrncnt(4),' of this type not all of which may be displayed.'
        if (10.0*fo_wrnlmt(4).gt.0.5*HUGE(fo_wrnlmt(4))) then
          fo_wrncnt(4) = 0
        else
          fo_wrnlmt(4) = fo_wrnlmt(4)*10
        end if
      end if
      inpm2(i) = 0.0
    else
      inpm2(i) = 1.0/(ncp123(i)*ncp234(i))
    end if
  end do
! get cosine
  inpm(:) = sqrt(inpm2(:))
  termc(:) = dotp(:)*inpm(:)
! filter for divergence near periodicity point
  do i=1,nrsdieff(rs)
!   this isn't great, but keeps it simple 
    termc(i) = max(min(termc(i),1.0d0),-1.0d0)
  end do
! from here on it's relatively straightforward
  do j=1,3
    dtmd1(:,j) = inpm(:)*(dotpd1(:,j) - &
 &       0.5*dotp(:)*inpm2(:)*&
 &         (ncp123(:)*ncp234d1(:,j) + ncp234(:)*ncp123d1(:,j)) )
    dtmd2(:,j) = inpm(:)*(dotpd2(:,j) - &
 &       0.5*dotp(:)*inpm2(:)*&
 &         (ncp123(:)*ncp234d2(:,j) + ncp234(:)*ncp123d2(:,j)) )
    dtmd3(:,j) = inpm(:)*(dotpd3(:,j) - &
 &       0.5*dotp(:)*inpm2(:)*&
 &         (ncp123(:)*ncp234d3(:,j) + ncp234(:)*ncp123d3(:,j)) )
    dtmd4(:,j) = inpm(:)*(dotpd4(:,j) - &
 &       0.5*dotp(:)*inpm2(:)*&
 &         (ncp123(:)*ncp234d4(:,j) + ncp234(:)*ncp123d4(:,j)) )
  end do
!
  do i=1,nrsdieff(rs)
    if (iaa(rs)%typ_di(i).eq.0) cycle
    ito1 = iaa(rs)%di(i,1)
    ito2 = iaa(rs)%di(i,2)
    ito3 = iaa(rs)%di(i,3)
    ito4 = iaa(rs)%di(i,4)
!    call numdihtest(ito1,ito2,ito3,ito4,fvec1(i,:),fvec2(i,:),&
! &fvec3(i,:),fvec4(i,:))
!   Ryckaert-Bell
    if (iaa(rs)%typ_di(i).eq.1) then
      evec(18) = evec(18) + globscal*(iaa(rs)%par_di(i,1) + &
 &                   iaa(rs)%par_di(i,2)*termc(i) +&
 &                   iaa(rs)%par_di(i,3)*ctp2 +&
 &                   iaa(rs)%par_di(i,4)*ctp3 +&
 &                   iaa(rs)%par_di(i,5)*ctp4 +&
 &                   iaa(rs)%par_di(i,6)*ctp5 +&
 &                   iaa(rs)%par_di(i,7)*ctp6 )
      term0 = globscal*(iaa(rs)%par_di(i,2) + &
 &            iaa(rs)%par_di(i,3)*2.0*termc(i) +&
 &            iaa(rs)%par_di(i,4)*3.0*ctp2 +&
 &            iaa(rs)%par_di(i,5)*4.0*ctp3 +&
 &            iaa(rs)%par_di(i,6)*5.0*ctp4 +&
 &            iaa(rs)%par_di(i,7)*6.0*ctp5 )
      do j=1,3
        ca_f(ito1,j) = ca_f(ito1,j) - term0*dtmd1(i,j)
        ca_f(ito2,j) = ca_f(ito2,j) - term0*dtmd2(i,j)
        ca_f(ito3,j) = ca_f(ito3,j) - term0*dtmd3(i,j)
        ca_f(ito4,j) = ca_f(ito4,j) - term0*dtmd4(i,j)
      end do
!   harmonic
    else if (iaa(rs)%typ_di(i).eq.2) then

      si = dv21(i,1)*cp234(i,1) + dv21(i,2)*cp234(i,2) + &
 &        dv21(i,3)*cp234(i,3)
      rang = RADIAN*acos(termc(i))
      frang = -RADIAN/sqrt(1.0 - termc(i)**2)
      if (si.lt.0.0) then
        rang = -rang
        frang = -frang
      end if
      term0 = globscal*iaa(rs)%par_di(i,1)
!     make sure periodicity is accounted for (note that the shift is a linear
!     function with slope unity such that that there is no extra chain rule derivative)
      vv = rang - iaa(rs)%par_di(i,2)
      if (vv.lt.-180.0) vv = vv + 360.0
      if (vv.gt.180.0) vv = vv - 360.0
      evec(18) = evec(18) + 0.5*term0*(vv**2)
      term0  = term0*vv*frang
      do j=1,3
        ca_f(ito1,j) = ca_f(ito1,j) - term0*dtmd1(i,j)
        ca_f(ito2,j) = ca_f(ito2,j) - term0*dtmd2(i,j)
        ca_f(ito3,j) = ca_f(ito3,j) - term0*dtmd3(i,j)
        ca_f(ito4,j) = ca_f(ito4,j) - term0*dtmd4(i,j)
      end do
    else
      write(ilog,*) 'Fatal. Encountered unsupported type of torsional p&
 &otential in force_torsions2(...) (offending type is ',&
 &iaa(rs)%typ_di(i),').'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! the exact same routine only for impropers (do NOT use this routine -> see above)
!
subroutine force_impropers2(rs,evec,ca_f)
!
  use iounit
  use params
  use energies
  use polypep
  use inter
  use atoms
  use forces
  use math
  use sequen
  use fyoc
!
  implicit none
!
  integer ito1,ito2,ito3,ito4,i,j,rs
  RTYPE dv21(nrsimpteff(rs),3),dv32(nrsimpteff(rs),3)
  RTYPE dv43(nrsimpteff(rs),3),ca_f(n,3)
  RTYPE ncp123(nrsimpteff(rs)),ncp234(nrsimpteff(rs))
  RTYPE dotp(nrsimpteff(rs)),rang,svec(3),svecref(3)
  RTYPE termc(nrsimpteff(rs))
  RTYPE frang,inpm(nrsimpteff(rs))
  RTYPE inpm2(nrsimpteff(rs))
  RTYPE dotpd1(nrsimpteff(rs),3)
  RTYPE dotpd2(nrsimpteff(rs),3)
  RTYPE dotpd3(nrsimpteff(rs),3),dotpd4(nrsimpteff(rs),3)
  RTYPE cp123(nrsimpteff(rs),3),cp234(nrsimpteff(rs),3)
  RTYPE ncp123d1(nrsimpteff(rs),3),ncp123d2(nrsimpteff(rs),3)
  RTYPE ncp123d3(nrsimpteff(rs),3),ncp123d4(nrsimpteff(rs),3)
  RTYPE ncp234d1(nrsimpteff(rs),3),ncp234d2(nrsimpteff(rs),3)
  RTYPE ncp234d3(nrsimpteff(rs),3),ncp234d4(nrsimpteff(rs),3)
  RTYPE dtmd1(nrsimpteff(rs),3),dtmd2(nrsimpteff(rs),3)
  RTYPE dtmd3(nrsimpteff(rs),3),globscal
  RTYPE dtmd4(nrsimpteff(rs),3),si,vv,ctp2,ctp3,ctp4,ctp5,ctp6
  RTYPE term0,evec(MAXENERGYTERMS)
!
  if (nrsimpteff(rs).eq.0) return
!
  globscal = scale_BOND(3)
  if ((use_FEG.EQV..true.).AND.(fegmode.eq.2)) then
    if (par_FEG(rs).EQV..true.) then
      if (use_FEGS(17).EQV..false.) then
        return
      else
        globscal = scale_FEGS(17)
      end if
    end if
  end if
!
  if (disulf(rs).gt.0) then
    if (molofrs(disulf(rs)).ne.molofrs(rs)) then
      call dis_bound_rs(rs,disulf(rs),svecref)
    else
      svecref(:) = 0.0
    end if
  else
    svecref(:) = 0.0
  end if
  do i=1,nrsimpteff(rs)
    ito1 = iaa(rs)%impt(i,improper_conv(1))
    ito2 = iaa(rs)%impt(i,2)
    ito3 = iaa(rs)%impt(i,improper_conv(2))
    ito4 = iaa(rs)%impt(i,4)
    svec(:) = 0.0
    if (atmres(ito2).eq.disulf(rs)) svec(:) = svecref(:)
    dv21(i,1) = x(ito2) - x(ito1) + svec(1)
    dv21(i,2) = y(ito2) - y(ito1) + svec(2)
    dv21(i,3) = z(ito2) - z(ito1) + svec(3)
    svec(:) = 0.0
    if ((atmres(ito3).eq.disulf(rs)).AND.(atmres(ito2).ne.disulf(rs))) svec(:) = svecref(:)
    if ((atmres(ito3).ne.disulf(rs)).AND.(atmres(ito2).eq.disulf(rs))) svec(:) = -svecref(:)
    dv32(i,1) = x(ito3) - x(ito2) + svec(1)
    dv32(i,2) = y(ito3) - y(ito2) + svec(2)
    dv32(i,3) = z(ito3) - z(ito2) + svec(3)
    svec(:) = 0.0
    if ((atmres(ito4).eq.disulf(rs)).AND.(atmres(ito3).ne.disulf(rs))) svec(:) = svecref(:)
    if ((atmres(ito4).ne.disulf(rs)).AND.(atmres(ito3).eq.disulf(rs))) svec(:) = -svecref(:)
    dv43(i,1) = x(ito4) - x(ito3) + svec(1)
    dv43(i,2) = y(ito4) - y(ito3) + svec(2)
    dv43(i,3) = z(ito4) - z(ito3) + svec(3)
  end do
  cp123(:,1) = dv21(:,2)*dv32(:,3) - dv21(:,3)*dv32(:,2)
  cp123(:,2) = dv21(:,3)*dv32(:,1) - dv21(:,1)*dv32(:,3)
  cp123(:,3) = dv21(:,1)*dv32(:,2) - dv21(:,2)*dv32(:,1)
  cp234(:,1) = dv32(:,2)*dv43(:,3) - dv32(:,3)*dv43(:,2)
  cp234(:,2) = dv32(:,3)*dv43(:,1) - dv32(:,1)*dv43(:,3)
  cp234(:,3) = dv32(:,1)*dv43(:,2) - dv32(:,2)*dv43(:,1)
! norms and dot-product
  ncp123(:) = cp123(:,1)**2 + cp123(:,2)**2 + cp123(:,3)**2
  ncp234(:) = cp234(:,1)**2 + cp234(:,2)**2 + cp234(:,3)**2
  dotp(:) = cp123(:,1)*cp234(:,1) + cp123(:,2)*cp234(:,2) + &
 &          cp123(:,3)*cp234(:,3)
! dot product derivatives
  dotpd1(:,1) = dv32(:,3)*cp234(:,2) - dv32(:,2)*cp234(:,3)
  dotpd1(:,2) = -dv32(:,3)*cp234(:,1) + dv32(:,1)*cp234(:,3)
  dotpd1(:,3) = dv32(:,2)*cp234(:,1) - dv32(:,1)*cp234(:,2)
  dotpd2(:,1) =  (-dv21(:,3)-dv32(:,3))*cp234(:,2) +&
 &                cp123(:,2)*dv43(:,3) + &
 &               (dv32(:,2)+dv21(:,2))*cp234(:,3) -&
 &                cp123(:,3)*dv43(:,2)
  dotpd2(:,2) =  (dv32(:,3)+dv21(:,3))*cp234(:,1) -&
 &                cp123(:,1)*dv43(:,3) + &
 &               (-dv21(:,1)-dv32(:,1))*cp234(:,3) +&
 &                cp123(:,3)*dv43(:,1)
  dotpd2(:,3) =  (-dv21(:,2)-dv32(:,2))*cp234(:,1) +&
 &                cp123(:,1)*dv43(:,2) + &
 &               (dv32(:,1)+dv21(:,1))*cp234(:,2) -&
 &                cp123(:,2)*dv43(:,1)
  dotpd3(:,1) =  (-dv32(:,3)-dv43(:,3))*cp123(:,2) +&
 &                cp234(:,2)*dv21(:,3) + &
 &               (dv43(:,2)+dv32(:,2))*cp123(:,3) -&
 &                cp234(:,3)*dv21(:,2)
  dotpd3(:,2) =  (dv43(:,3)+dv32(:,3))*cp123(:,1) -&
 &                cp234(:,1)*dv21(:,3) + &
 &               (-dv32(:,1)-dv43(:,1))*cp123(:,3) +&
 &                cp234(:,3)*dv21(:,1)
  dotpd3(:,3) =  (-dv32(:,2)-dv43(:,2))*cp123(:,1) +&
 &                cp234(:,1)*dv21(:,2) + &
 &               (dv43(:,1)+dv32(:,1))*cp123(:,2) -&
 &                cp234(:,2)*dv21(:,1)
  dotpd4(:,1) = dv32(:,3)*cp123(:,2) - dv32(:,2)*cp123(:,3)
  dotpd4(:,2) = -dv32(:,3)*cp123(:,1) + dv32(:,1)*cp123(:,3)
  dotpd4(:,3) = dv32(:,2)*cp123(:,1) - dv32(:,1)*cp123(:,2) 
! quadratic norm derivatives
  ncp123d1(:,1) = 2.0*(cp123(:,2)*dv32(:,3) - cp123(:,3)*dv32(:,2))
  ncp123d1(:,2) = 2.0*(-cp123(:,1)*dv32(:,3) + cp123(:,3)*dv32(:,1))
  ncp123d1(:,3) = 2.0*(cp123(:,1)*dv32(:,2) - cp123(:,2)*dv32(:,1))
  ncp234d1(:,:) = 0.0
  ncp123d2(:,1) = 2.0*(cp123(:,2)*(-dv21(:,3)-dv32(:,3)) + &
 &                     cp123(:,3)*(dv32(:,2)+dv21(:,2)))
  ncp123d2(:,2) = 2.0*(cp123(:,1)*(dv32(:,3)+dv21(:,3)) + &
 &                     cp123(:,3)*(-dv21(:,1)-dv32(:,1)))
  ncp123d2(:,3) = 2.0*(cp123(:,1)*(-dv21(:,2)-dv32(:,2)) + &
 &                     cp123(:,2)*(dv32(:,1)+dv21(:,1)))
  ncp234d2(:,1) = 2.0*(cp234(:,2)*dv43(:,3) - cp234(:,3)*dv43(:,2))
  ncp234d2(:,2) = 2.0*(-cp234(:,1)*dv43(:,3) + cp234(:,3)*dv43(:,1))
  ncp234d2(:,3) = 2.0*(cp234(:,1)*dv43(:,2) - cp234(:,2)*dv43(:,1))
  ncp123d3(:,1) = 2.0*(cp123(:,2)*dv21(:,3) - cp123(:,3)*dv21(:,2))
  ncp123d3(:,2) = 2.0*(-cp123(:,1)*dv21(:,3) + cp123(:,3)*dv21(:,1))
  ncp123d3(:,3) = 2.0*(cp123(:,1)*dv21(:,2) - cp123(:,2)*dv21(:,1))
  ncp234d3(:,1) = 2.0*(cp234(:,2)*(-dv32(:,3)-dv43(:,3)) + &
 &                     cp234(:,3)*(dv43(:,2)+dv32(:,2)))
  ncp234d3(:,2) = 2.0*(cp234(:,1)*(dv43(:,3)+dv32(:,3)) + &
 &                     cp234(:,3)*(-dv32(:,1)-dv43(:,1)))
  ncp234d3(:,3) = 2.0*(cp234(:,1)*(-dv32(:,2)-dv43(:,2)) + &
 &                     cp234(:,2)*(dv43(:,1)+dv32(:,1)))
  ncp123d4(:,:) = 0.0
  ncp234d4(:,1) = 2.0*(cp234(:,2)*dv32(:,3) - cp234(:,3)*dv32(:,2))
  ncp234d4(:,2) = 2.0*(-cp234(:,1)*dv32(:,3) + cp234(:,3)*dv32(:,1))
  ncp234d4(:,3) = 2.0*(cp234(:,1)*dv32(:,2) - cp234(:,2)*dv32(:,1))
! filter for exceptions which would otherwise NaN
  do i=1,nrsimpteff(rs)
!   colinear reference atoms
    if ((ncp123(i).le.0.0).OR.(ncp234(i).le.0.0)) then
!     this will set all derivatives to zero and the (now arbitrary) dihedral to zero deg.
      inpm2(i) = 0.0
      fo_wrncnt(5) = fo_wrncnt(5) + 1
      if (fo_wrncnt(5).eq.fo_wrnlmt(5)) then
        write(ilog,*) 'WARNING. Colinear reference atoms in computation o&
 &f improper torsional forces (acos). This indicates bad parameters, an unstable simulation, or a bug&
 &. Expect run to crash soon.'
        write(ilog,*) 'This was warning #',fo_wrncnt(5),' of this type not all of which may be displayed.'
        if (10.0*fo_wrnlmt(5).gt.0.5*HUGE(fo_wrnlmt(5))) then
          fo_wrncnt(5) = 0
        else
          fo_wrnlmt(5) = fo_wrnlmt(5)*10
        end if
      end if
    else
      inpm2(i) = 1.0/(ncp123(i)*ncp234(i))
    end if
  end do
! get cosine
  inpm(:) = sqrt(inpm2(:))
  termc(:) = dotp(:)*inpm(:)
! filter for divergence near periodicity point
  do i=1,nrsimpteff(rs)
!   this isn't great, but keeps it simple 
    termc(i) = max(min(termc(i),1.0d0),-1.0d0)
  end do
! from here on it's relatively straightforward
  do j=1,3
    dtmd1(:,j) = inpm(:)*(dotpd1(:,j) - &
 &       0.5*dotp(:)*inpm2(:)*&
 &         (ncp123(:)*ncp234d1(:,j) + ncp234(:)*ncp123d1(:,j)) )
    dtmd2(:,j) = inpm(:)*(dotpd2(:,j) - &
 &       0.5*dotp(:)*inpm2(:)*&
 &         (ncp123(:)*ncp234d2(:,j) + ncp234(:)*ncp123d2(:,j)) )
    dtmd3(:,j) = inpm(:)*(dotpd3(:,j) - &
 &       0.5*dotp(:)*inpm2(:)*&
 &         (ncp123(:)*ncp234d3(:,j) + ncp234(:)*ncp123d3(:,j)) )
    dtmd4(:,j) = inpm(:)*(dotpd4(:,j) - &
 &       0.5*dotp(:)*inpm2(:)*&
 &         (ncp123(:)*ncp234d4(:,j) + ncp234(:)*ncp123d4(:,j)) )
  end do
!
  do i=1,nrsimpteff(rs)
    if (iaa(rs)%typ_impt(i).eq.0) cycle
    ito1 = iaa(rs)%impt(i,improper_conv(1))
    ito2 = iaa(rs)%impt(i,2)
    ito3 = iaa(rs)%impt(i,improper_conv(2))
    ito4 = iaa(rs)%impt(i,4)
!   Ryckaert-Bell (cosine powers)
    if (iaa(rs)%typ_impt(i).eq.1) then
      ctp2 = termc(i)*termc(i)
      ctp3 = ctp2*termc(i)
      ctp4 = ctp3*termc(i)
      ctp5 = ctp4*termc(i)
      ctp6 = ctp5*termc(i)
      evec(17) = evec(17) + globscal*(iaa(rs)%par_impt(i,1) + &
 &                iaa(rs)%par_impt(i,2)*termc(i) +&
 &                iaa(rs)%par_impt(i,3)*ctp2 +&
 &                iaa(rs)%par_impt(i,4)*ctp3 +&
 &                iaa(rs)%par_impt(i,5)*ctp4 +&
 &                iaa(rs)%par_impt(i,6)*ctp5 +&
 &                iaa(rs)%par_impt(i,7)*ctp6 )
      term0 = globscal*(iaa(rs)%par_impt(i,2) + &
 &            iaa(rs)%par_impt(i,3)*2.0*termc(i) +&
 &            iaa(rs)%par_impt(i,4)*3.0*ctp2 +&
 &            iaa(rs)%par_impt(i,5)*4.0*ctp3 +&
 &            iaa(rs)%par_impt(i,6)*5.0*ctp4 +&
 &            iaa(rs)%par_impt(i,7)*6.0*ctp5)
      do j=1,3
        ca_f(ito1,j) = ca_f(ito1,j) - term0*dtmd1(i,j)
        ca_f(ito2,j) = ca_f(ito2,j) - term0*dtmd2(i,j)
        ca_f(ito3,j) = ca_f(ito3,j) - term0*dtmd3(i,j)
        ca_f(ito4,j) = ca_f(ito4,j) - term0*dtmd4(i,j)
      end do
!   harmonic
    else if (iaa(rs)%typ_impt(i).eq.2) then
      si = dv21(i,1)*cp234(i,1) + dv21(i,2)*cp234(i,2) + &
 &        dv21(i,3)*cp234(i,3)
      rang = RADIAN*acos(termc(i))
      frang = -RADIAN/sqrt(1.0 - termc(i)**2)
      if (si.lt.0.0) then
        rang = -rang
        frang = -frang
      end if
      term0 = globscal*iaa(rs)%par_impt(i,1)
!     make sure periodicity is accounted for (note that the shift is a linear
!     function with slope unity such that there is no extra chain rule derivative)
      vv = rang - iaa(rs)%par_impt(i,2)
      if (vv.lt.-180.0) vv = vv + 360.0
      if (vv.gt.180.0) vv = vv - 360.0
      evec(17) = evec(17) + 0.5*term0*(vv**2)
      term0  = term0*vv*frang
      do j=1,3
        ca_f(ito1,j) = ca_f(ito1,j) - term0*dtmd1(i,j)
        ca_f(ito2,j) = ca_f(ito2,j) - term0*dtmd2(i,j)
        ca_f(ito3,j) = ca_f(ito3,j) - term0*dtmd3(i,j)
        ca_f(ito4,j) = ca_f(ito4,j) - term0*dtmd4(i,j)
      end do
    else
      write(ilog,*) 'Fatal. Encountered unsupported type of torsional p&
 &otential in force_impropers2(...) (offending type is ',&
 &iaa(rs)%typ_impt(i),').'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine numdihtest(iba1,iba2,iba3,iba4,fv1,fv2,fv3,fv4)
!
  use atoms
  use math
!
  implicit none
!
  integer i,ttt(4),iba1,iba2,iba3,iba4
  RTYPE bla,getztor,blas(10),fv1(3),fv2(3),fv3(3),fv4(3)
  RTYPE tvv(4,3)
!
  bla = 0.0001
!
  ttt(1) = iba1
  ttt(2) = iba2
  ttt(3) = iba3
  ttt(4) = iba4
  tvv(1,:) = RADIAN*fv1(:)
  tvv(2,:) = RADIAN*fv2(:)
  tvv(3,:) = RADIAN*fv3(:)
  tvv(4,:) = RADIAN*fv4(:)
  write(*,*) 'atoms: ',ttt(:)
  do i=1,4
    x(ttt(i)) = x(ttt(i)) + bla
    blas(2) = getztor(iba1,iba2,iba3,iba4)
    x(ttt(i)) = x(ttt(i)) - 2.0*bla
    blas(3) = getztor(iba1,iba2,iba3,iba4)
    x(ttt(i)) = x(ttt(i)) + bla
    if ((blas(3)-blas(2)).gt.180.0) blas(3) = blas(3) - 360.0
    if ((blas(3)-blas(2)).lt.-180.0) blas(3) = blas(3) + 360.0
    write(*,*) tvv(i,1),(blas(3)-blas(2))/(2.0*bla)
    y(ttt(i)) = y(ttt(i)) + bla
    blas(2) = getztor(iba1,iba2,iba3,iba4)
    y(ttt(i)) = y(ttt(i)) - 2.0*bla
    blas(3) = getztor(iba1,iba2,iba3,iba4)
    y(ttt(i)) = y(ttt(i)) + bla
    if ((blas(3)-blas(2)).gt.180.0) blas(3) = blas(3) - 360.0
    if ((blas(3)-blas(2)).lt.-180.0) blas(3) = blas(3) + 360.0
    write(*,*) tvv(i,2),(blas(3)-blas(2))/(2.0*bla)
    z(ttt(i)) = z(ttt(i)) + bla
    blas(2) = getztor(iba1,iba2,iba3,iba4)
    z(ttt(i)) = z(ttt(i)) - 2.0*bla
    blas(3) = getztor(iba1,iba2,iba3,iba4)
    z(ttt(i)) = z(ttt(i)) + bla
    if ((blas(3)-blas(2)).gt.180.0) blas(3) = blas(3) - 360.0
    if ((blas(3)-blas(2)).lt.-180.0) blas(3) = blas(3) + 360.0
    write(*,*) tvv(i,3),(blas(3)-blas(2))/(2.0*bla)
  end do
!
end
!
!-----------------------------------------------------------------------
!
! Blondel and Karplus derived the strategy below for getting around singularities
! it takes the (old) re-route through dphi/dr = dphi/dcosphi * dcosphi/dr = 
! -1/sinphi * dcosphi/dr
! however, instead of evaluating dcosphi/dr directly it circumvents the singularity
! in sinphi through merging the sin into the second derivative, in which in the 
! dangerous limiting cases continuity is preserved, and no trigonometric fxns are
! needed
! colinear reference atoms will still lead to the necessary singularity
! (as the dihedral angle does in fact become ill-defined in those cases)
!
subroutine force_torsions(rs,evec,ca_f)
!
  use iounit
  use params
  use energies
  use polypep
  use inter
  use atoms
  use forces
  use math
  use fyoc
  use sequen
!
  implicit none
!
  integer ito1,ito2,ito3,ito4,i,j,rs
  RTYPE dv12(nrsdieff(rs),3),dv23(nrsdieff(rs),3)
  RTYPE dv43(nrsdieff(rs),3),ndv23(nrsdieff(rs))
  RTYPE incp123(nrsdieff(rs)),incp234(nrsdieff(rs))
  RTYPE ncp123(nrsdieff(rs)),ncp234(nrsdieff(rs))
  RTYPE dotp(nrsdieff(rs)),dotpfg(nrsdieff(rs))
  RTYPE cp123(nrsdieff(rs),3),cp234(nrsdieff(rs),3)
  RTYPE dtmd1(nrsdieff(rs),3),dtmd2(nrsdieff(rs),3)
  RTYPE dtmd3(nrsdieff(rs),3),dotphg(nrsdieff(rs))
  RTYPE dtmd4(nrsdieff(rs),3),dotpah(nrsdieff(rs)),si,vv
  RTYPE termc(nrsdieff(rs)),terms(nrsdieff(rs))
  RTYPE inpm(nrsdieff(rs)),globscal,svec(3),svecref(3)
  RTYPE term0,evec(MAXENERGYTERMS),st2,ct2,cmin,smin
  RTYPE ctp2,ctp3,ctp4,ctp5,ctp6,ca_f(n,3)
!
  if (nrsdieff(rs).eq.0) return
!
  globscal = scale_BOND(4)
  if ((use_FEG.EQV..true.).AND.(fegmode.eq.2)) then
    if (par_FEG(rs).EQV..true.) then
      if (use_FEGS(18).EQV..false.) then
        return
      else
        globscal = scale_FEGS(18)
      end if
    end if
  end if
!
  if (disulf(rs).gt.0) then
    if (molofrs(disulf(rs)).ne.molofrs(rs)) then
      call dis_bound_rs(rs,disulf(rs),svecref)
    else
      svecref(:) = 0.0
    end if
  else
    svecref(:) = 0.0
  end if
  do i=1,nrsdieff(rs)
    ito1 = iaa(rs)%di(i,1)
    ito2 = iaa(rs)%di(i,2)
    ito3 = iaa(rs)%di(i,3)
    ito4 = iaa(rs)%di(i,4)
    svec(:) = 0.0
    if (atmres(ito2).eq.disulf(rs)) svec(:) = -svecref(:)
    dv12(i,1) = x(ito1) - x(ito2) + svec(1)
    dv12(i,2) = y(ito1) - y(ito2) + svec(2)
    dv12(i,3) = z(ito1) - z(ito2) + svec(3)
    svec(:) = 0.0
    if ((atmres(ito3).eq.disulf(rs)).AND.(atmres(ito2).ne.disulf(rs))) svec(:) = -svecref(:)
    if ((atmres(ito3).ne.disulf(rs)).AND.(atmres(ito2).eq.disulf(rs))) svec(:) = svecref(:)
    dv23(i,1) = x(ito2) - x(ito3) + svec(1)
    dv23(i,2) = y(ito2) - y(ito3) + svec(2)
    dv23(i,3) = z(ito2) - z(ito3) + svec(3)
    svec(:) = 0.0
    if ((atmres(ito4).eq.disulf(rs)).AND.(atmres(ito3).ne.disulf(rs))) svec(:) = svecref(:)
    if ((atmres(ito4).ne.disulf(rs)).AND.(atmres(ito3).eq.disulf(rs))) svec(:) = -svecref(:)
    dv43(i,1) = x(ito4) - x(ito3) + svec(1)
    dv43(i,2) = y(ito4) - y(ito3) + svec(2)
    dv43(i,3) = z(ito4) - z(ito3) + svec(3)
  end do
  cp123(:,1) = dv12(:,2)*dv23(:,3) - dv12(:,3)*dv23(:,2)
  cp123(:,2) = dv12(:,3)*dv23(:,1) - dv12(:,1)*dv23(:,3)
  cp123(:,3) = dv12(:,1)*dv23(:,2) - dv12(:,2)*dv23(:,1)
  cp234(:,1) = dv43(:,2)*dv23(:,3) - dv43(:,3)*dv23(:,2)
  cp234(:,2) = dv43(:,3)*dv23(:,1) - dv43(:,1)*dv23(:,3)
  cp234(:,3) = dv43(:,1)*dv23(:,2) - dv43(:,2)*dv23(:,1)
! norms and dot-product
  ncp123(:) = cp123(:,1)**2 + cp123(:,2)**2 + cp123(:,3)**2
  ncp234(:) = cp234(:,1)**2 + cp234(:,2)**2 + cp234(:,3)**2
  dotp(:) = cp123(:,1)*cp234(:,1) + cp123(:,2)*cp234(:,2) + &
 &          cp123(:,3)*cp234(:,3)
  dotpfg(:) = dv12(:,1)*dv23(:,1) + dv12(:,2)*dv23(:,2) +&
 &           dv12(:,3)*dv23(:,3)
  dotphg(:) = dv43(:,1)*dv23(:,1) + dv43(:,2)*dv23(:,2) +&
 &           dv43(:,3)*dv23(:,3)
  dotpah(:) = cp123(:,1)*dv43(:,1) + cp123(:,2)*dv43(:,2) +&
 &           cp123(:,3)*dv43(:,3)
  ndv23(:) = sqrt(dv23(:,1)**2 + dv23(:,2)**2 + dv23(:,3)**2)
!
! filter for exceptions which would otherwise NaN
  do i=1,nrsdieff(rs)
!   colinear reference atoms (note that this also detects the equally fatal |dv23| = 0)
    if ((ncp123(i).le.0.0).OR.(ncp234(i).le.0.0)) then
!     this will set all derivatives to zero and the (now arbitrary) dihedral to zero deg.
      fo_wrncnt(1) = fo_wrncnt(1) + 1
      if (fo_wrncnt(1).eq.fo_wrnlmt(1)) then
        write(ilog,*) 'WARNING. Colinear reference atoms in computation o&
 &f torsional forces (Blondel-Karplus). This indicates bad parameters, an unstable simulation, or a bug&
 &. Expect run to crash soon.'
        write(ilog,*) 'This was warning #',fo_wrncnt(1),' of this type not all of which may be displayed.'
        if (10.0*fo_wrnlmt(1).gt.0.5*HUGE(fo_wrnlmt(1))) then
          fo_wrncnt(1) = 0
        else
          fo_wrnlmt(1) = fo_wrnlmt(1)*10
        end if
      end if
      incp123(i) = 0.0
      incp234(i) = 0.0
    else
      incp123(i) = 1.0/ncp123(i)
      incp234(i) = 1.0/ncp234(i)
    end if
  end do
!
! equations 27 in Blondel and Karplus
  do j=1,3
    dtmd1(:,j) = -cp123(:,j)*incp123(:)*ndv23(:)

    dtmd2(:,j) = cp123(:,j)*incp123(:)*&
 &                  (ndv23(:) + dotpfg(:)/ndv23(:)) - &
 &               cp234(:,j)*dotphg(:)*incp234(:)/ndv23(:)

    dtmd3(:,j) = cp234(:,j)*incp234(:)*&
 &                  (dotphg(:)/ndv23(:) - ndv23(:)) - &
 &               cp123(:,j)*dotpfg(:)*incp123(:)/ndv23(:)
    dtmd4(:,j) = cp234(:,j)*incp234(:)*ndv23(:)
  end do
!
! now get the cosine and sine term (sine term from Blondel equ. 49)
  inpm(:) = sqrt(incp123(:)*incp234(:))
  termc(:) = dotp(:)*inpm(:)
  terms(:) = inpm(:)*ndv23(:)*dotpah(:)
!
  do i=1,nrsdieff(rs)
    if (iaa(rs)%typ_di(i).eq.0) cycle
    ito1 = iaa(rs)%di(i,1)
    ito2 = iaa(rs)%di(i,2)
    ito3 = iaa(rs)%di(i,3)
    ito4 = iaa(rs)%di(i,4)
!   Ryckaert-Bellemans
    if (iaa(rs)%typ_di(i).eq.1) then
      ctp2 = termc(i)*termc(i)
      ctp3 = ctp2*termc(i)
      ctp4 = ctp3*termc(i)
      ctp5 = ctp4*termc(i)
      ctp6 = ctp5*termc(i)
      evec(18) = evec(18) + globscal*(iaa(rs)%par_di(i,1) + &
 &                   iaa(rs)%par_di(i,2)*termc(i) +&
 &                   iaa(rs)%par_di(i,3)*ctp2 +&
 &                   iaa(rs)%par_di(i,4)*ctp3 +&
 &                   iaa(rs)%par_di(i,5)*ctp4 +&
 &                   iaa(rs)%par_di(i,6)*ctp5 +&
 &                   iaa(rs)%par_di(i,7)*ctp6 )
!     we generate the cosine-derivative (dtermc/dr) through -terms*(dphi/dr)
!     (standard chain rule)
      term0 = -terms(i)*globscal*(iaa(rs)%par_di(i,2) + &
 &            iaa(rs)%par_di(i,3)*2.0*termc(i) +&
 &            iaa(rs)%par_di(i,4)*3.0*ctp2 +&
 &            iaa(rs)%par_di(i,5)*4.0*ctp3 +&
 &            iaa(rs)%par_di(i,6)*5.0*ctp4 +&
 &            iaa(rs)%par_di(i,7)*6.0*ctp5 )
      do j=1,3
        ca_f(ito1,j) = ca_f(ito1,j) - term0*dtmd1(i,j)
        ca_f(ito2,j) = ca_f(ito2,j) - term0*dtmd2(i,j)
        ca_f(ito3,j) = ca_f(ito3,j) - term0*dtmd3(i,j)
        ca_f(ito4,j) = ca_f(ito4,j) - term0*dtmd4(i,j)
      end do
!   harmonic
    else if (iaa(rs)%typ_di(i).eq.2) then

      cmin = cos(iaa(rs)%par_di(i,2)/RADIAN) 
      smin = sin(iaa(rs)%par_di(i,2)/RADIAN)
!     get cos(phi-phi_0) and sin(phi-phi_0)
      ct2 = cmin*termc(i) + smin*terms(i)
      st2 = cmin*terms(i) - smin*termc(i)
      if (st2.lt.0.0) then
        si = -1.0
      else
        si = 1.0
      end if
      if (abs(ct2).gt.0.1) then ! safe regime for asin
        vv = RADIAN*asin(st2)
        if (ct2.lt.0.0) vv = 180.0 - vv ! safe for asin but angle very far away
      else ! switch to cosine (means angle is far away from target)
        vv = si*RADIAN*acos(ct2)
      end if
      term0 = globscal*iaa(rs)%par_di(i,1)
!     make sure periodicity is accounted for (note that the shift is a linear
!     function with slope unity such that that there is no extra chain rule derivative)
      if (vv.lt.-180.0) vv = vv + 360.0
      if (vv.gt.180.0) vv = vv - 360.0
      evec(18) = evec(18) + 0.5*term0*(vv**2)
      term0  = RADIAN*term0*vv
      do j=1,3
        ca_f(ito1,j) = ca_f(ito1,j) - term0*dtmd1(i,j)
        ca_f(ito2,j) = ca_f(ito2,j) - term0*dtmd2(i,j)
        ca_f(ito3,j) = ca_f(ito3,j) - term0*dtmd3(i,j)
        ca_f(ito4,j) = ca_f(ito4,j) - term0*dtmd4(i,j)
      end do
    else
      write(ilog,*) 'Fatal. Encountered unsupported type of torsional p&
 &otential in force_torsions(...) (offending type is ',&
 &iaa(rs)%typ_di(i),').'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! same for impropers
!
subroutine force_impropers(rs,evec,ca_f)
!
  use iounit
  use params
  use energies
  use polypep
  use inter
  use atoms
  use forces
  use math
  use sequen
  use fyoc
!
  implicit none
!
  integer ito1,ito2,ito3,ito4,i,j,rs
  RTYPE dv12(nrsimpteff(rs),3),dv23(nrsimpteff(rs),3)
  RTYPE dv43(nrsimpteff(rs),3),ndv23(nrsimpteff(rs))
  RTYPE incp123(nrsimpteff(rs)),incp234(nrsimpteff(rs))
  RTYPE ncp123(nrsimpteff(rs)),ncp234(nrsimpteff(rs))
  RTYPE dotp(nrsimpteff(rs)),dotpfg(nrsimpteff(rs))
  RTYPE cp123(nrsimpteff(rs),3),cp234(nrsimpteff(rs),3)
  RTYPE dtmd1(nrsimpteff(rs),3),dtmd2(nrsimpteff(rs),3)
  RTYPE dtmd3(nrsimpteff(rs),3),dotphg(nrsimpteff(rs))
  RTYPE dtmd4(nrsimpteff(rs),3),dotpah(nrsimpteff(rs)),si,vv
  RTYPE termc(nrsimpteff(rs)),terms(nrsimpteff(rs))
  RTYPE inpm(nrsimpteff(rs)),globscal,svec(3),svecref(3)
  RTYPE term0,evec(MAXENERGYTERMS),st2,ct2,cmin,smin
  RTYPE ctp2,ctp3,ctp4,ctp5,ctp6,ca_f(n,3)
!
  if (nrsimpteff(rs).eq.0) return
!
  globscal = scale_BOND(3)
  if ((use_FEG.EQV..true.).AND.(fegmode.eq.2)) then
    if (par_FEG(rs).EQV..true.) then
      if (use_FEGS(17).EQV..false.) then
        return
      else
        globscal = scale_FEGS(17)
      end if
    end if
  end if
!
  if (disulf(rs).gt.0) then
    if (molofrs(disulf(rs)).ne.molofrs(rs)) then
      call dis_bound_rs(rs,disulf(rs),svecref)
    else
      svecref(:) = 0.0
    end if
  else
    svecref(:) = 0.0
  end if
  do i=1,nrsimpteff(rs)
    ito1 = iaa(rs)%impt(i,improper_conv(1))
    ito2 = iaa(rs)%impt(i,2)
    ito3 = iaa(rs)%impt(i,improper_conv(2))
    ito4 = iaa(rs)%impt(i,4)
    svec(:) = 0.0
    if (atmres(ito2).eq.disulf(rs)) svec(:) = -svecref(:)
    dv12(i,1) = x(ito1) - x(ito2) + svec(1)
    dv12(i,2) = y(ito1) - y(ito2) + svec(2)
    dv12(i,3) = z(ito1) - z(ito2) + svec(3)
    svec(:) = 0.0
    if ((atmres(ito3).eq.disulf(rs)).AND.(atmres(ito2).ne.disulf(rs))) svec(:) = -svecref(:)
    if ((atmres(ito3).ne.disulf(rs)).AND.(atmres(ito2).eq.disulf(rs))) svec(:) = svecref(:)
    dv23(i,1) = x(ito2) - x(ito3) + svec(1)
    dv23(i,2) = y(ito2) - y(ito3) + svec(2)
    dv23(i,3) = z(ito2) - z(ito3) + svec(3)
    svec(:) = 0.0
    if ((atmres(ito4).eq.disulf(rs)).AND.(atmres(ito3).ne.disulf(rs))) svec(:) = svecref(:)
    if ((atmres(ito4).ne.disulf(rs)).AND.(atmres(ito3).eq.disulf(rs))) svec(:) = -svecref(:)
    dv43(i,1) = x(ito4) - x(ito3) + svec(1)
    dv43(i,2) = y(ito4) - y(ito3) + svec(2)
    dv43(i,3) = z(ito4) - z(ito3) + svec(3)
  end do
  cp123(:,1) = dv12(:,2)*dv23(:,3) - dv12(:,3)*dv23(:,2)
  cp123(:,2) = dv12(:,3)*dv23(:,1) - dv12(:,1)*dv23(:,3)
  cp123(:,3) = dv12(:,1)*dv23(:,2) - dv12(:,2)*dv23(:,1)
  cp234(:,1) = dv43(:,2)*dv23(:,3) - dv43(:,3)*dv23(:,2)
  cp234(:,2) = dv43(:,3)*dv23(:,1) - dv43(:,1)*dv23(:,3)
  cp234(:,3) = dv43(:,1)*dv23(:,2) - dv43(:,2)*dv23(:,1)
! norms and dot-product
  ncp123(:) = cp123(:,1)**2 + cp123(:,2)**2 + cp123(:,3)**2
  ncp234(:) = cp234(:,1)**2 + cp234(:,2)**2 + cp234(:,3)**2
  dotp(:) = cp123(:,1)*cp234(:,1) + cp123(:,2)*cp234(:,2) + &
 &          cp123(:,3)*cp234(:,3)
  dotpfg(:) = dv12(:,1)*dv23(:,1) + dv12(:,2)*dv23(:,2) +&
 &           dv12(:,3)*dv23(:,3)
  dotphg(:) = dv43(:,1)*dv23(:,1) + dv43(:,2)*dv23(:,2) +&
 &           dv43(:,3)*dv23(:,3)
  dotpah(:) = cp123(:,1)*dv43(:,1) + cp123(:,2)*dv43(:,2) +&
 &           cp123(:,3)*dv43(:,3)
  ndv23(:) = sqrt(dv23(:,1)**2 + dv23(:,2)**2 + dv23(:,3)**2)
!
! filter for exceptions which would otherwise NaN
  do i=1,nrsimpteff(rs)
!   colinear reference atoms (note that this also detects the equally fatal |dv23| = 0)
    if ((ncp123(i).le.0.0).OR.(ncp234(i).le.0.0)) then
!     this will set all derivatives to zero and the (now arbitrary) dihedral to zero deg.
      fo_wrncnt(2) = fo_wrncnt(2) + 1
      if (fo_wrncnt(2).eq.fo_wrnlmt(2)) then
        write(ilog,*) 'WARNING. Colinear reference atoms in computation o&
 &f improper torsional forces (Blondel-Karplus). This indicates bad parameters, an unstable simulation, or a bug&
 &. Expect run to crash soon.'
        write(ilog,*) 'This was warning #',fo_wrncnt(2),' of this type not all of which may be displayed.'
        if (10.0*fo_wrnlmt(2).gt.0.5*HUGE(fo_wrnlmt(2))) then
          fo_wrncnt(2) = 0
        else
          fo_wrnlmt(2) = fo_wrnlmt(2)*10
        end if
      end if
      incp123(i) = 0.0
      incp234(i) = 0.0
    else
      incp123(i) = 1.0/ncp123(i)
      incp234(i) = 1.0/ncp234(i)
    end if
  end do
!
! equations 27 in Blondel and Karplus
  do j=1,3
    dtmd1(:,j) = -cp123(:,j)*incp123(:)*ndv23(:)

    dtmd2(:,j) = cp123(:,j)*incp123(:)*&
 &                  (ndv23(:) + dotpfg(:)/ndv23(:)) - &
 &               cp234(:,j)*dotphg(:)*incp234(:)/ndv23(:)
    dtmd3(:,j) = cp234(:,j)*incp234(:)*&
 &                  (dotphg(:)/ndv23(:) - ndv23(:)) - &
 &               cp123(:,j)*dotpfg(:)*incp123(:)/ndv23(:)
    dtmd4(:,j) = cp234(:,j)*incp234(:)*ndv23(:)
  end do
!
! now get the cosine and sine term (sine term from Blondel equ. 49)
  inpm(:) = sqrt(incp123(:)*incp234(:))
  termc(:) = dotp(:)*inpm(:)
  terms(:) = inpm(:)*ndv23(:)*dotpah(:)
!
  do i=1,nrsimpteff(rs)
    if (iaa(rs)%typ_impt(i).eq.0) cycle
    ito1 = iaa(rs)%impt(i,improper_conv(1))
    ito2 = iaa(rs)%impt(i,2)
    ito3 = iaa(rs)%impt(i,improper_conv(2))
    ito4 = iaa(rs)%impt(i,4)
!   Ryckaert-Bellemans
    if (iaa(rs)%typ_impt(i).eq.1) then
      ctp2 = termc(i)*termc(i)
      ctp3 = ctp2*termc(i)
      ctp4 = ctp3*termc(i)
      ctp5 = ctp4*termc(i)
      ctp6 = ctp5*termc(i)
      evec(17) = evec(17) + globscal*(iaa(rs)%par_impt(i,1) + &
 &                iaa(rs)%par_impt(i,2)*termc(i) +&
 &                iaa(rs)%par_impt(i,3)*ctp2 +&
 &                iaa(rs)%par_impt(i,4)*ctp3 +&
 &                iaa(rs)%par_impt(i,5)*ctp4 +&
 &                iaa(rs)%par_impt(i,6)*ctp5 +&
 &                iaa(rs)%par_impt(i,7)*ctp6 )
!     we generate the cosine-derivative (dtermc/dr) through -terms*(dphi/dr)
!     (standard chain rule)
      term0 = -terms(i)*globscal*(iaa(rs)%par_impt(i,2) + &
 &            iaa(rs)%par_impt(i,3)*2.0*termc(i) +&
 &            iaa(rs)%par_impt(i,4)*3.0*ctp2 +&
 &            iaa(rs)%par_impt(i,5)*4.0*ctp3 +&
 &            iaa(rs)%par_impt(i,6)*5.0*ctp4 +&
 &            iaa(rs)%par_impt(i,7)*6.0*ctp5 )
      do j=1,3
        ca_f(ito1,j) = ca_f(ito1,j) - term0*dtmd1(i,j)
        ca_f(ito2,j) = ca_f(ito2,j) - term0*dtmd2(i,j)
        ca_f(ito3,j) = ca_f(ito3,j) - term0*dtmd3(i,j)
        ca_f(ito4,j) = ca_f(ito4,j) - term0*dtmd4(i,j)
      end do
!   harmonic
    else if (iaa(rs)%typ_impt(i).eq.2) then

      cmin = cos(iaa(rs)%par_impt(i,2)/RADIAN) 
      smin = sin(iaa(rs)%par_impt(i,2)/RADIAN)
!     get cos(phi-phi_0) and sin(phi-phi_0)
      ct2 = cmin*termc(i) + smin*terms(i)
      st2 = cmin*terms(i) - smin*termc(i)
      if (st2.lt.0.0) then
        si = -1.0
      else
        si = 1.0
      end if
      if (abs(ct2).gt.0.1) then ! safe regime for asin
        vv = RADIAN*asin(st2)
        if (ct2.lt.0.0) vv = 180.0 - vv ! safe for asin but angle very far away 
      else ! switch to cosine (means angle is far away from target)
        vv = si*RADIAN*acos(ct2)
      end if
      term0 = globscal*iaa(rs)%par_impt(i,1)
!     make sure periodicity is accounted for (note that the shift is a linear
!     function with slope unity such that that there is no extra chain rule derivative)
      if (vv.lt.-180.0) vv = vv + 360.0
      if (vv.gt.180.0) vv = vv - 360.0
      evec(17) = evec(17) + 0.5*term0*(vv**2)
      term0  = RADIAN*term0*vv
      do j=1,3
        ca_f(ito1,j) = ca_f(ito1,j) - term0*dtmd1(i,j)
        ca_f(ito2,j) = ca_f(ito2,j) - term0*dtmd2(i,j)
        ca_f(ito3,j) = ca_f(ito3,j) - term0*dtmd3(i,j)
        ca_f(ito4,j) = ca_f(ito4,j) - term0*dtmd4(i,j)
      end do
    else
      write(ilog,*) 'Fatal. Encountered unsupported type of torsional p&
 &otential in force_impropers(...) (offending type is ',&
 &iaa(rs)%typ_impt(i),').'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! WARNING: this function assumes that no crosslinked atoms contribute
!          to the CMAP dihedrals
!
subroutine force_CMAP(rs,evec,ca_f)
!
  use iounit
  use params
  use energies
  use inter
  use atoms
  use math
  use forces
!
  implicit none
!
  integer ito1,ito2,ito3,ito4,ito5,i,j,rs
  RTYPE dv12(nrscmeff(rs),3),dv23(nrscmeff(rs),3),dv54(nrscmeff(rs),3)
  RTYPE dv43(nrscmeff(rs),3),ndv23(nrscmeff(rs)),ndv43(nrscmeff(rs))
  RTYPE incp123(nrscmeff(rs)),incp234(nrscmeff(rs)),incp345(nrscmeff(rs))
  RTYPE ncp123(nrscmeff(rs)),ncp234(nrscmeff(rs)),ncp345(nrscmeff(rs))
  RTYPE dotp(nrscmeff(rs)),dotp2(nrscmeff(rs))
  RTYPE cp123(nrscmeff(rs),3),cp234(nrscmeff(rs),3),cp345(nrscmeff(rs),3)
  RTYPE dotpah(nrscmeff(rs)),dotpah2(nrscmeff(rs))
  RTYPE termc(nrscmeff(rs)),terms(nrscmeff(rs))
  RTYPE termc2(nrscmeff(rs)),terms2(nrscmeff(rs))
  RTYPE inpm(nrscmeff(rs)),inpm2(nrscmeff(rs))
  RTYPE dotpih(nrscmeff(rs)),dotpfg(nrscmeff(rs)),dotphg(nrscmeff(rs))
  RTYPE dtmd1(nrscmeff(rs),3),dtmd1b(nrscmeff(rs),3)
  RTYPE dtmd2(nrscmeff(rs),3),dtmd2b(nrscmeff(rs),3)
  RTYPE dtmd3(nrscmeff(rs),3),dtmd3b(nrscmeff(rs),3)
  RTYPE dtmd4(nrscmeff(rs),3),dtmd4b(nrscmeff(rs),3)
  RTYPE globscal,si,si2,vv,vv2,dvv,dvv2
  RTYPE term0,evec(MAXENERGYTERMS)
  RTYPE ca_f(n,3)
!
  if (nrscmeff(rs).eq.0) return
!
  globscal = scale_BOND(5)
  if ((use_FEG.EQV..true.).AND.(fegmode.eq.2)) then
    if (par_FEG(rs).EQV..true.) then
      if (use_FEGS(20).EQV..false.) then
        return
      else
        globscal = scale_FEGS(20)
      end if
    end if
  end if
!
  do i=1,nrscmeff(rs)
    ito1 = iaa(rs)%cm(i,1)
    ito2 = iaa(rs)%cm(i,2)
    ito3 = iaa(rs)%cm(i,3)
    ito4 = iaa(rs)%cm(i,4)
    ito5 = iaa(rs)%cm(i,5)
    dv12(i,1) = x(ito1) - x(ito2)
    dv12(i,2) = y(ito1) - y(ito2)
    dv12(i,3) = z(ito1) - z(ito2)
    dv23(i,1) = x(ito2) - x(ito3)
    dv23(i,2) = y(ito2) - y(ito3)
    dv23(i,3) = z(ito2) - z(ito3)
    dv43(i,1) = x(ito4) - x(ito3)
    dv43(i,2) = y(ito4) - y(ito3)
    dv43(i,3) = z(ito4) - z(ito3)
    dv54(i,1) = x(ito5) - x(ito4)
    dv54(i,2) = y(ito5) - y(ito4)
    dv54(i,3) = z(ito5) - z(ito4)
  end do
  cp123(:,1) = dv12(:,2)*dv23(:,3) - dv12(:,3)*dv23(:,2)
  cp123(:,2) = dv12(:,3)*dv23(:,1) - dv12(:,1)*dv23(:,3)
  cp123(:,3) = dv12(:,1)*dv23(:,2) - dv12(:,2)*dv23(:,1)
  cp234(:,1) = dv43(:,2)*dv23(:,3) - dv43(:,3)*dv23(:,2)
  cp234(:,2) = dv43(:,3)*dv23(:,1) - dv43(:,1)*dv23(:,3)
  cp234(:,3) = dv43(:,1)*dv23(:,2) - dv43(:,2)*dv23(:,1)
  cp345(:,1) = -dv54(:,2)*dv43(:,3) + dv54(:,3)*dv43(:,2)
  cp345(:,2) = -dv54(:,3)*dv43(:,1) + dv54(:,1)*dv43(:,3)
  cp345(:,3) = -dv54(:,1)*dv43(:,2) + dv54(:,2)*dv43(:,1)

! norms and dot-product
  ncp123(:) = cp123(:,1)**2 + cp123(:,2)**2 + cp123(:,3)**2
  ncp234(:) = cp234(:,1)**2 + cp234(:,2)**2 + cp234(:,3)**2
  ncp345(:) = cp345(:,1)**2 + cp345(:,2)**2 + cp345(:,3)**2

  dotp(:) = cp123(:,1)*cp234(:,1) + cp123(:,2)*cp234(:,2) + &
 &          cp123(:,3)*cp234(:,3)
  dotp2(:) = cp234(:,1)*cp345(:,1) + cp234(:,2)*cp345(:,2) + &
 &          cp234(:,3)*cp345(:,3)
  dotpah(:) = cp123(:,1)*dv43(:,1) + cp123(:,2)*dv43(:,2) +&
 &           cp123(:,3)*dv43(:,3)
  dotpah2(:) = cp234(:,1)*dv54(:,1) + cp234(:,2)*dv54(:,2) +&
 &           cp234(:,3)*dv54(:,3)

  dotpfg(:) = dv12(:,1)*dv23(:,1) + dv12(:,2)*dv23(:,2) +&
 &           dv12(:,3)*dv23(:,3)
  dotphg(:) = dv43(:,1)*dv23(:,1) + dv43(:,2)*dv23(:,2) +&
 &           dv43(:,3)*dv23(:,3)
  dotpih(:) = -dv54(:,1)*dv43(:,1) - dv54(:,2)*dv43(:,2) -&
 &           dv54(:,3)*dv43(:,3)
!
  ndv23(:) = sqrt(dv23(:,1)**2 + dv23(:,2)**2 + dv23(:,3)**2)
  ndv43(:) = sqrt(dv43(:,1)**2 + dv43(:,2)**2 + dv43(:,3)**2)
!
! filter for exceptions which would otherwise NaN
  do i=1,nrscmeff(rs)
!   colinear reference atoms (note that this also detects the equally fatal |dv23| = 0)
    if ((ncp123(i).le.0.0).OR.(ncp234(i).le.0.0)) then
!     this will set all derivatives to zero and the (now arbitrary) dihedral to zero deg.
      fo_wrncnt(3) = fo_wrncnt(3) + 1
      if (fo_wrncnt(3).eq.fo_wrnlmt(3)) then
        write(ilog,*) 'WARNING. Colinear reference atoms in computation o&
 &f CMAP forces (Blondel-Karplus). This indicates bad parameters, an unstable simulation, or a bug&
 &. Expect run to crash soon.'
        write(ilog,*) 'This was warning #',fo_wrncnt(3),' of this type not all of which may be displayed.'
        if (10.0*fo_wrnlmt(3).gt.0.5*HUGE(fo_wrnlmt(3))) then
          fo_wrncnt(3) = 0
        else
          fo_wrnlmt(3) = fo_wrnlmt(3)*10
        end if
      end if
      incp123(i) = 0.0
      incp234(i) = 0.0
    else
      incp123(i) = 1.0/ncp123(i)
      incp234(i) = 1.0/ncp234(i)
    end if
    if ((ncp234(i).le.0.0).OR.(ncp345(i).le.0.0)) then
!     this will set all derivatives to zero and the (now arbitrary) dihedral to zero deg.
      fo_wrncnt(3) = fo_wrncnt(3) + 1
      if (fo_wrncnt(3).eq.fo_wrnlmt(3)) then
        write(ilog,*) 'WARNING. Colinear reference atoms in computation o&
 &f CMAP forces (Blondel-Karplus). This indicates bad parameters, an unstable simulation, or a bug&
 &. Expect run to crash soon.'
        write(ilog,*) 'This was warning #',fo_wrncnt(3),' of this type not all of which may be displayed.'
        if (10.0*fo_wrnlmt(3).gt.0.5*HUGE(fo_wrnlmt(3))) then
          fo_wrncnt(3) = 0
        else
          fo_wrnlmt(3) = fo_wrnlmt(3)*10
        end if
      end if
      incp345(i) = 0.0
      incp234(i) = 0.0
    else
      incp345(i) = 1.0/ncp345(i)
      incp234(i) = 1.0/ncp234(i)
    end if
  end do
!
! equations 27 in Blondel and Karplus
  do j=1,3
    dtmd1(:,j) = -cp123(:,j)*incp123(:)*ndv23(:)
    dtmd2(:,j) = cp123(:,j)*incp123(:)*&
 &                  (ndv23(:) + dotpfg(:)/ndv23(:)) - &
 &               cp234(:,j)*dotphg(:)*incp234(:)/ndv23(:)
    dtmd3(:,j) = cp234(:,j)*incp234(:)*&
 &                  (dotphg(:)/ndv23(:) - ndv23(:)) - &
 &               cp123(:,j)*dotpfg(:)*incp123(:)/ndv23(:)
    dtmd4(:,j) = cp234(:,j)*incp234(:)*ndv23(:)
!
    dtmd1b(:,j) = -cp234(:,j)*incp234(:)*ndv43(:)
    dtmd2b(:,j) = cp234(:,j)*incp234(:)*&
 &                  (ndv43(:) - dotphg(:)/ndv43(:)) - &
 &               cp345(:,j)*dotpih(:)*incp345(:)/ndv43(:)
    dtmd3b(:,j) = cp345(:,j)*incp345(:)*&
 &                  (dotpih(:)/ndv43(:) - ndv43(:)) + &
 &               cp234(:,j)*dotphg(:)*incp234(:)/ndv43(:)
    dtmd4b(:,j) = cp345(:,j)*incp345(:)*ndv43(:)
  end do
!
! now get the cosine and sine terms (sine term from Blondel equ. 49)
  inpm(:) = sqrt(incp123(:)*incp234(:))
  inpm2(:) = sqrt(incp234(:)*incp345(:))
  termc(:) = dotp(:)*inpm(:)
  termc2(:) = dotp2(:)*inpm2(:)
  terms(:) = inpm(:)*ndv23(:)*dotpah(:)
  terms2(:) = inpm2(:)*ndv43(:)*dotpah2(:)
!
  do i=1,nrscmeff(rs)
    if (iaa(rs)%typ_cm(i).eq.0) cycle
!
    if (terms(i).lt.0.0) then
      si = -1.0
    else
      si = 1.0
    end if
    if (abs(termc(i)).gt.0.1) then ! safe regime for asin
      vv = RADIAN*asin(terms(i))
      if (termc(i).lt.0.0) vv = 180.0 - vv
    else 
      vv = si*RADIAN*acos(termc(i))
    end if
    if (terms2(i).lt.0.0) then
      si2 = -1.0
    else
      si2 = 1.0
    end if
    if (abs(termc2(i)).gt.0.1) then ! safe regime for asin
      vv2 = RADIAN*asin(terms2(i))
      if (termc2(i).lt.0.0) vv2 = 180.0 - vv2
    else ! switch to cosine (means angle is far away from zero)
      vv2 = si2*RADIAN*acos(termc2(i))
    end if
    if (vv.lt.-180.0) vv = vv + 360.0
    if (vv.gt.180.0) vv = vv - 360.0
    if (vv2.lt.-180.0) vv2 = vv2 + 360.0
    if (vv2.gt.180.0) vv2 = vv2 - 360.0

    if ((cm_typ(iaa(rs)%typ_cm(i),1).eq.1).OR.(cm_typ(iaa(rs)%typ_cm(i),1).eq.2)) then
      call cardinal_bspline2D(iaa(rs)%typ_cm(i),vv,vv2,term0,dvv,dvv2)
      evec(20) = evec(20) + globscal*term0
      dvv = globscal*RADIAN*dvv
      dvv2 = globscal*RADIAN*dvv2
      do j=1,3
        ca_f(ito1,j) = ca_f(ito1,j) - dvv*dtmd1(i,j)
        ca_f(ito2,j) = ca_f(ito2,j) - dvv*dtmd2(i,j) - dvv2*dtmd1b(i,j)
        ca_f(ito3,j) = ca_f(ito3,j) - dvv*dtmd3(i,j) - dvv2*dtmd2b(i,j)
        ca_f(ito4,j) = ca_f(ito4,j) - dvv*dtmd4(i,j) - dvv2*dtmd3b(i,j)
        ca_f(ito5,j) = ca_f(ito5,j) - dvv2*dtmd4b(i,j)
      end do
    else if ((cm_typ(iaa(rs)%typ_cm(i),1).eq.3).OR.(cm_typ(iaa(rs)%typ_cm(i),1).eq.4)) then
      call bicubic_spline(iaa(rs)%typ_cm(i),vv,vv2,term0,dvv,dvv2)
      evec(20) = evec(20) + globscal*term0
      dvv = globscal*RADIAN*dvv
      dvv2 = globscal*RADIAN*dvv2
      do j=1,3
        ca_f(ito1,j) = ca_f(ito1,j) - dvv*dtmd1(i,j)
        ca_f(ito2,j) = ca_f(ito2,j) - dvv*dtmd2(i,j) - dvv2*dtmd1b(i,j)
        ca_f(ito3,j) = ca_f(ito3,j) - dvv*dtmd3(i,j) - dvv2*dtmd2b(i,j)
        ca_f(ito4,j) = ca_f(ito4,j) - dvv*dtmd4(i,j) - dvv2*dtmd3b(i,j)
        ca_f(ito5,j) = ca_f(ito5,j) - dvv2*dtmd4b(i,j)
      end do
    else
      write(ilog,*) 'Fatal. Encountered unsupported type of CMAP p&
 &otential in force_CMAP(...) (offending type is ',&
 &cm_typ(iaa(rs)%typ_cm(i),1),').'
      call fexit()
    end if
  end do
!
end
!
!--------------------------------------------------------------------------------------------------
!
subroutine onetor_cosderiv(ito1,ito2,ito3,ito4,&
 &                         termc,dtmd1,dtmd2,dtmd3,dtmd4,si)
!
  use iounit
  use params
  use energies
  use polypep
  use inter
  use atoms
  use forces
  use math
!
  implicit none
!
  integer ito1,ito2,ito3,ito4,j
  RTYPE dv21(3),dv32(3)
  RTYPE dv43(3)
  RTYPE ncp123,ncp234,dotp,termc,inpm,inpm2
  RTYPE dotpd1(3),dotpd2(3),dotpd3(3),dotpd4(3)
  RTYPE cp123(3),cp234(3)
  RTYPE ncp123d1(3),ncp123d2(3)
  RTYPE ncp123d3(3),ncp123d4(3)
  RTYPE ncp234d1(3),ncp234d2(3)
  RTYPE ncp234d3(3),ncp234d4(3)
  RTYPE dtmd1(3),dtmd2(3),dtmd3(3),dtmd4(3)
  RTYPE si,sip
!
  dv21(1) = x(ito2) - x(ito1)
  dv21(2) = y(ito2) - y(ito1)
  dv21(3) = z(ito2) - z(ito1)
  dv32(1) = x(ito3) - x(ito2)
  dv32(2) = y(ito3) - y(ito2)
  dv32(3) = z(ito3) - z(ito2)
  dv43(1) = x(ito4) - x(ito3)
  dv43(2) = y(ito4) - y(ito3)
  dv43(3) = z(ito4) - z(ito3)
  cp123(1) = dv21(2)*dv32(3) - dv21(3)*dv32(2)
  cp123(2) = dv21(3)*dv32(1) - dv21(1)*dv32(3)
  cp123(3) = dv21(1)*dv32(2) - dv21(2)*dv32(1)
  cp234(1) = dv32(2)*dv43(3) - dv32(3)*dv43(2)
  cp234(2) = dv32(3)*dv43(1) - dv32(1)*dv43(3)
  cp234(3) = dv32(1)*dv43(2) - dv32(2)*dv43(1)
! norms and dot-product
  ncp123 = cp123(1)**2 + cp123(2)**2 + cp123(3)**2
  ncp234 = cp234(1)**2 + cp234(2)**2 + cp234(3)**2
  dotp = cp123(1)*cp234(1) + cp123(2)*cp234(2) + &
 &          cp123(3)*cp234(3)
! dot product derivatives
  dotpd1(1) = dv32(3)*cp234(2) - dv32(2)*cp234(3)
  dotpd1(2) = -dv32(3)*cp234(1) + dv32(1)*cp234(3)
  dotpd1(3) = dv32(2)*cp234(1) - dv32(1)*cp234(2)
  dotpd2(1) =  (-dv21(3)-dv32(3))*cp234(2) +&
 &                cp123(2)*dv43(3) + &
 &               (dv32(2)+dv21(2))*cp234(3) -&
 &                cp123(3)*dv43(2)
  dotpd2(2) =  (dv32(3)+dv21(3))*cp234(1) -&
 &                cp123(1)*dv43(3) + &
 &               (-dv21(1)-dv32(1))*cp234(3) +&
 &                cp123(3)*dv43(1)
  dotpd2(3) =  (-dv21(2)-dv32(2))*cp234(1) +&
 &                cp123(1)*dv43(2) + &
 &               (dv32(1)+dv21(1))*cp234(2) -&
 &                cp123(2)*dv43(1)
  dotpd3(1) =  (-dv32(3)-dv43(3))*cp123(2) +&
 &                cp234(2)*dv21(3) + &
 &               (dv43(2)+dv32(2))*cp123(3) -&
 &                cp234(3)*dv21(2)
  dotpd3(2) =  (dv43(3)+dv32(3))*cp123(1) -&
 &                cp234(1)*dv21(3) + &
 &               (-dv32(1)-dv43(1))*cp123(3) +&
 &                cp234(3)*dv21(1)
  dotpd3(3) =  (-dv32(2)-dv43(2))*cp123(1) +&
 &                cp234(1)*dv21(2) + &
 &               (dv43(1)+dv32(1))*cp123(2) -&
 &                cp234(2)*dv21(1)
  dotpd4(1) = dv32(3)*cp123(2) - dv32(2)*cp123(3)
  dotpd4(2) = -dv32(3)*cp123(1) + dv32(1)*cp123(3)
  dotpd4(3) = dv32(2)*cp123(1) - dv32(1)*cp123(2) 
! quadratic norm derivatives
  ncp123d1(1) = 2.0*(cp123(2)*dv32(3) - cp123(3)*dv32(2))
  ncp123d1(2) = 2.0*(-cp123(1)*dv32(3) + cp123(3)*dv32(1))
  ncp123d1(3) = 2.0*(cp123(1)*dv32(2) - cp123(2)*dv32(1))
  ncp234d1(:) = 0.0
  ncp123d2(1) = 2.0*(cp123(2)*(-dv21(3)-dv32(3)) + &
 &                     cp123(3)*(dv32(2)+dv21(2)))
  ncp123d2(2) = 2.0*(cp123(1)*(dv32(3)+dv21(3)) + &
 &                     cp123(3)*(-dv21(1)-dv32(1)))
  ncp123d2(3) = 2.0*(cp123(1)*(-dv21(2)-dv32(2)) + &
 &                     cp123(2)*(dv32(1)+dv21(1)))
  ncp234d2(1) = 2.0*(cp234(2)*dv43(3) - cp234(3)*dv43(2))
  ncp234d2(2) = 2.0*(-cp234(1)*dv43(3) + cp234(3)*dv43(1))
  ncp234d2(3) = 2.0*(cp234(1)*dv43(2) - cp234(2)*dv43(1))
  ncp123d3(1) = 2.0*(cp123(2)*dv21(3) - cp123(3)*dv21(2))
  ncp123d3(2) = 2.0*(-cp123(1)*dv21(3) + cp123(3)*dv21(1))
  ncp123d3(3) = 2.0*(cp123(1)*dv21(2) - cp123(2)*dv21(1))
  ncp234d3(1) = 2.0*(cp234(2)*(-dv32(3)-dv43(3)) + &
 &                     cp234(3)*(dv43(2)+dv32(2)))
  ncp234d3(2) = 2.0*(cp234(1)*(dv43(3)+dv32(3)) + &
 &                     cp234(3)*(-dv32(1)-dv43(1)))
  ncp234d3(3) = 2.0*(cp234(1)*(-dv32(2)-dv43(2)) + &
 &                     cp234(2)*(dv43(1)+dv32(1)))
  ncp123d4(:) = 0.0
  ncp234d4(1) = 2.0*(cp234(2)*dv32(3) - cp234(3)*dv32(2))
  ncp234d4(2) = 2.0*(-cp234(1)*dv32(3) + cp234(3)*dv32(1))
  ncp234d4(3) = 2.0*(cp234(1)*dv32(2) - cp234(2)*dv32(1))
! filter for exceptions which would otherwise NaN
! colinear reference atoms
  if ((ncp123.le.0.0).OR.(ncp234.le.0.0)) then
!   this will set all derivatives to zero and the (now arbitrary) dihedral to zero deg.
    inpm2 = 0.0
  else
    inpm2 = 1.0/(ncp123*ncp234)
  end if
! get cosine
  inpm = sqrt(inpm2)
  termc = dotp*inpm
! filter for divergence near periodicity point
  termc = max(min(termc,1.0d0),-1.0d0)
! from here on it's relatively straightforward
  do j=1,3
    dtmd1(j) = inpm*(dotpd1(j) - &
 &       0.5*dotp*inpm2*&
 &         (ncp123*ncp234d1(j) + ncp234*ncp123d1(j)) )
    dtmd2(j) = inpm*(dotpd2(j) - &
 &       0.5*dotp*inpm2*&
 &         (ncp123*ncp234d2(j) + ncp234*ncp123d2(j)) )
    dtmd3(j) = inpm*(dotpd3(j) - &
 &       0.5*dotp*inpm2*&
 &         (ncp123*ncp234d3(j) + ncp234*ncp123d3(j)) )
    dtmd4(j) = inpm*(dotpd4(j) - &
 &       0.5*dotp*inpm2*&
 &         (ncp123*ncp234d4(j) + ncp234*ncp123d4(j)) )
  end do
!
  sip = dv21(1)*cp234(1) + dv21(2)*cp234(2) + dv21(3)*cp234(3)
  if (sip.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
!
end
!
!-----------------------------------------------------------------------
!
! this routine uses the elegant Blondel formalism to get dphi/dr for
! a single torsion with dcosphi/dr and dsinphi/dr readily available through terms and termc
!
subroutine onetor_deriv(ito1,ito2,ito3,ito4,termc,terms,&
 &                         fvec1,fvec2,fvec3,fvec4)
!
  use math
  use atoms
  use iounit
  use forces
!
  implicit none
!
  integer j,ito1,ito2,ito3,ito4
  RTYPE dv12(3),dv23(3),dv43(3),cp123(3),cp234(3)
  RTYPE ncp123,ncp234,dotp,dotpfg,dotphg,dotpah,ndv23
  RTYPE incp123,incp234,inpm,termc,terms
  RTYPE fvec1(3),fvec2(3),fvec3(3),fvec4(3)
!
  dv12(1) = x(ito1) - x(ito2)
  dv12(2) = y(ito1) - y(ito2)
  dv12(3) = z(ito1) - z(ito2)
  dv23(1) = x(ito2) - x(ito3)
  dv23(2) = y(ito2) - y(ito3)
  dv23(3) = z(ito2) - z(ito3)
  dv43(1) = x(ito4) - x(ito3)
  dv43(2) = y(ito4) - y(ito3)
  dv43(3) = z(ito4) - z(ito3)
  cp123(1) = dv12(2)*dv23(3) - dv12(3)*dv23(2)
  cp123(2) = dv12(3)*dv23(1) - dv12(1)*dv23(3)
  cp123(3) = dv12(1)*dv23(2) - dv12(2)*dv23(1)
  cp234(1) = dv43(2)*dv23(3) - dv43(3)*dv23(2)
  cp234(2) = dv43(3)*dv23(1) - dv43(1)*dv23(3)
  cp234(3) = dv43(1)*dv23(2) - dv43(2)*dv23(1)
! norms and dot-product
  ncp123 = cp123(1)**2 + cp123(2)**2 + cp123(3)**2
  ncp234 = cp234(1)**2 + cp234(2)**2 + cp234(3)**2
  dotp = cp123(1)*cp234(1) + cp123(2)*cp234(2) + &
 &          cp123(3)*cp234(3)
  dotpfg = dv12(1)*dv23(1) + dv12(2)*dv23(2) +&
 &           dv12(3)*dv23(3)
  dotphg = dv43(1)*dv23(1) + dv43(2)*dv23(2) +&
 &           dv43(3)*dv23(3)
  dotpah = cp123(1)*dv43(1) + cp123(2)*dv43(2) +&
 &           cp123(3)*dv43(3)
  ndv23 = sqrt(dv23(1)**2 + dv23(2)**2 + dv23(3)**2)
!
! filter for exceptions which would otherwise NaN
! colinear reference atoms (note that this also detects the equally fatal |dv23| = 0)
  if ((ncp123.le.0.0).OR.(ncp234.le.0.0)) then
!   this will set all derivatives to zero and the (now arbitrary) dihedral to zero deg.
    fo_wrncnt(6) = fo_wrncnt(6) + 1
    if (fo_wrncnt(6).eq.fo_wrnlmt(6)) then
      write(ilog,*) 'WARNING. Colinear reference atoms in onetor_deriv(...)&
 & (Blondel-Karplus). This indicates bad parameters, an unstable simulation, or a bug&
 &. Expect run to crash soon.'
      write(ilog,*) 'This was warning #',fo_wrncnt(6),' of this type not all of which may be displayed.'
      if (10.0*fo_wrnlmt(6).gt.0.5*HUGE(fo_wrnlmt(6))) then
        fo_wrncnt(6) = 0
      else
        fo_wrnlmt(6) = fo_wrnlmt(6)*10
      end if
    end if
    incp123 = 0.0
    incp234 = 0.0
  else
    incp123 = 1.0/ncp123
    incp234 = 1.0/ncp234
  end if
!
! equations 27 in Blondel and Karplus in deg/A
  do j=1,3
    fvec1(j) = -RADIAN*cp123(j)*incp123*ndv23

    fvec2(j) = RADIAN*(cp123(j)*incp123*&
 &                  (ndv23 + dotpfg/ndv23) - &
 &               cp234(j)*dotphg*incp234/ndv23)
    fvec3(j) = RADIAN*(cp234(j)*incp234*&
 &                  (dotphg/ndv23 - ndv23) - &
 &               cp123(j)*dotpfg*incp123/ndv23)
    fvec4(j) = RADIAN*cp234(j)*incp234*ndv23
  end do
  inpm = sqrt(incp123*incp234)
  termc = dotp*inpm
  terms = inpm*ndv23*dotpah
!
end
!
!-----------------------------------------------------------------------
!
! the same operating on xref, yref, zref
!
subroutine onetor_deriv_ref(ito1,ito2,ito3,ito4,termc,terms,&
 &                         fvec1,fvec2,fvec3,fvec4)
!
  use math
  use atoms
  use iounit
  use forces
!
  implicit none
!
  integer j,ito1,ito2,ito3,ito4
  RTYPE dv12(3),dv23(3),dv43(3),cp123(3),cp234(3)
  RTYPE ncp123,ncp234,dotp,dotpfg,dotphg,dotpah,ndv23
  RTYPE incp123,incp234,inpm,termc,terms
  RTYPE fvec1(3),fvec2(3),fvec3(3),fvec4(3)
!
  dv12(1) = xref(ito1) - xref(ito2)
  dv12(2) = yref(ito1) - yref(ito2)
  dv12(3) = zref(ito1) - zref(ito2)
  dv23(1) = xref(ito2) - xref(ito3)
  dv23(2) = yref(ito2) - yref(ito3)
  dv23(3) = zref(ito2) - zref(ito3)
  dv43(1) = xref(ito4) - xref(ito3)
  dv43(2) = yref(ito4) - yref(ito3)
  dv43(3) = zref(ito4) - zref(ito3)
  cp123(1) = dv12(2)*dv23(3) - dv12(3)*dv23(2)
  cp123(2) = dv12(3)*dv23(1) - dv12(1)*dv23(3)
  cp123(3) = dv12(1)*dv23(2) - dv12(2)*dv23(1)
  cp234(1) = dv43(2)*dv23(3) - dv43(3)*dv23(2)
  cp234(2) = dv43(3)*dv23(1) - dv43(1)*dv23(3)
  cp234(3) = dv43(1)*dv23(2) - dv43(2)*dv23(1)
! norms and dot-product
  ncp123 = cp123(1)**2 + cp123(2)**2 + cp123(3)**2
  ncp234 = cp234(1)**2 + cp234(2)**2 + cp234(3)**2
  dotp = cp123(1)*cp234(1) + cp123(2)*cp234(2) + &
 &          cp123(3)*cp234(3)
  dotpfg = dv12(1)*dv23(1) + dv12(2)*dv23(2) +&
 &           dv12(3)*dv23(3)
  dotphg = dv43(1)*dv23(1) + dv43(2)*dv23(2) +&
 &           dv43(3)*dv23(3)
  dotpah = cp123(1)*dv43(1) + cp123(2)*dv43(2) +&
 &           cp123(3)*dv43(3)
  ndv23 = sqrt(dv23(1)**2 + dv23(2)**2 + dv23(3)**2)
!
! filter for exceptions which would otherwise NaN
! colinear reference atoms (note that this also detects the equally fatal |dv23| = 0)
  if ((ncp123.le.0.0).OR.(ncp234.le.0.0)) then
!   this will set all derivatives to zero and the (now arbitrary) dihedral to zero deg.
    fo_wrncnt(6) = fo_wrncnt(6) + 1
    if (fo_wrncnt(6).eq.fo_wrnlmt(6)) then
      write(ilog,*) 'WARNING. Colinear reference atoms in onetor_deriv(...)&
 & (Blondel-Karplus). This indicates bad parameters, an unstable simulation, or a bug&
 &. Expect run to crash soon.'
      write(ilog,*) 'This was warning #',fo_wrncnt(6),' of this type not all of which may be displayed.'
      if (10.0*fo_wrnlmt(6).gt.0.5*HUGE(fo_wrnlmt(6))) then
        fo_wrncnt(6) = 0
      else
        fo_wrnlmt(6) = fo_wrnlmt(6)*10
      end if
    end if
    incp123 = 0.0
    incp234 = 0.0
  else
    incp123 = 1.0/ncp123
    incp234 = 1.0/ncp234
  end if
!
! equations 27 in Blondel and Karplus in deg/A
  do j=1,3
    fvec1(j) = -RADIAN*cp123(j)*incp123*ndv23

    fvec2(j) = RADIAN*(cp123(j)*incp123*&
 &                  (ndv23 + dotpfg/ndv23) - &
 &               cp234(j)*dotphg*incp234/ndv23)
    fvec3(j) = RADIAN*(cp234(j)*incp234*&
 &                  (dotphg/ndv23 - ndv23) - &
 &               cp123(j)*dotpfg*incp123/ndv23)
    fvec4(j) = RADIAN*cp234(j)*incp234*ndv23
  end do
  inpm = sqrt(incp123*incp234)
  termc = dotp*inpm
  terms = inpm*ndv23*dotpah
!
end
!
!-----------------------------------------------------------------------
!
! usage discouraged: for dphi/dr use onetor_deriv instead
!
subroutine onetor_angderiv(ito1,ito2,ito3,ito4,rang,&
 &                         fvec1,fvec2,fvec3,fvec4)
!
  use iounit
  use math
!
  implicit none
!
  integer ito1,ito2,ito3,ito4,j
  RTYPE dtmd1(3),dtmd2(3),dtmd3(3),dtmd4(3)
  RTYPE fvec1(3),fvec2(3),fvec3(3),fvec4(3)
  RTYPE si,ct,rang,frang
!
! note that the cosine term is numerically adjusted in this call 
! to make the acos-call safe
  call onetor_cosderiv(ito1,ito2,ito3,ito4,&
 &                         ct,dtmd1,dtmd2,dtmd3,dtmd4,si)
  rang = RADIAN*acos(ct)
  frang = -RADIAN/sqrt(1.0 - ct**2)
! adjust handedness
  rang = si*rang
  frang = si*frang
!
! derivative drang/dri for all four atoms
  do j=1,3
    fvec1(j) = frang*dtmd1(j)
    fvec2(j) = frang*dtmd2(j)
    fvec3(j) = frang*dtmd3(j)
    fvec4(j) = frang*dtmd4(j)
  end do
!
end
!
!-----------------------------------------------------------------------
!
! note that this routine handles both torsional derivatives as well
! as Cartesian ones, whichever is needed(!)
!
subroutine force_torrs(evec,rs,ca_f)
!
  use iounit
  use energies
  use fyoc
  use forces
  use math
  use sequen
  use system
  use zmatrix
  use aminos
  use atoms
  use torsn
  use wl
!
  implicit none
!
  RTYPE evec(MAXENERGYTERMS),vv(MAXTORRES),expte,dt0,ca_f(n,3)
  RTYPE ad1(3),ad2(3),ad3(3),ad4(3),ctm,stm,getpuckertor,torincr
  integer k,rs,j,imol,vvz(MAXTORRES),zl,vvi(MAXTORRES)
!
  if (par_TOR2(rs).eq.0) then
    write(ilog,*) 'Fatal. Called force_torrs(rs) with residue, for w&
 &hich no torsional potential was set up.'
    call fexit()
  end if
!
!
  imol = molofrs(rs)
  vv(:) = 0.0
  vvz(:) = 0
  vvi(:) = 0
!
  if (seqpolty(rs).eq.'N') then
    do k=1,nnucs(rs)
      if (par_TOR(rs,MAXTORRES+k).gt.0) then
        vv(k) = nucs(k,rs) - par_TOR(rs,k)
        vvz(k) = nucsline(k,rs)
        vvi(k) = nucsnr(k,rs)
      end if
    end do
    vv(nnucs(rs)+1) = getpuckertor(rs,zl) - par_TOR(rs,nnucs(rs)+1)
    vvz(nnucs(rs)+1) = zl
  else if ((fline(rs).gt.0).OR.(yline(rs).gt.0).OR.(wline(rs).gt.0)) then
    if (par_TOR(rs,MAXTORRES+1).gt.0.0) then
      vv(1) = omega(rs) - par_TOR(rs,1)
      vvz(1) = wline(rs)
      vvi(1) = wnr(rs)
    end if
    if (par_TOR(rs,MAXTORRES+2).gt.0.0) then
      vv(2) = phi(rs) - par_TOR(rs,2)
      vvz(2) = fline(rs)
      vvi(2) = fnr(rs)
    end if
    if (par_TOR(rs,MAXTORRES+3).gt.0.0) then
      vv(3) = psi(rs) - par_TOR(rs,3)
      vvz(3) = yline(rs)
      vvi(3) = ynr(rs)
    end if
  end if
  do k=1,nchi(rs)
    if (par_TOR(rs,MAXTORRES+6+k).gt.0) then
      vv(6+k) = chi(k,rs) - par_TOR(rs,6+k)
      vvz(6+k) = chiline(k,rs)
      vvi(6+k) = chinr(k,rs)
    end if
  end do
! periodic boundary conditions (so to speak)
  do j=1,MAXTORRES
    if (vv(j).lt.-180.0) vv(j) = vv(j) + 360.0
    if (vv(j).gt.180.0) vv(j) = vv(j) - 360.0
  end do
  torincr = 0.0
  do j=1,MAXTORRES
    if ((par_TOR(rs,MAXTORRES+j).gt.0.0).AND.(vvz(j).gt.0)) then
      torincr = torincr + par_TOR(rs,MAXTORRES+j)*vv(j)*vv(j)
    end if
  end do
! the potential is harmonic  
  if ((par_TOR2(rs).eq.1).OR.(par_TOR2(rs).eq.3)) then
    evec(7) = evec(7) + scale_TOR*torincr
    if ((fycxyz.eq.1).AND.(hmjam%isin(7).EQV..false.).AND.(hmjam%isin(11).EQV..false.)) then
      do j=1,MAXTORRES
!       note that the vvi-check will cycle out the sugar in nucleotides
        if ((vvz(j).gt.0).AND.(vvi(j).gt.0)) then
          dc_di(imol)%f(vvi(j)) = dc_di(imol)%f(vvi(j)) - &
 &               scale_TOR*2.0*RADIAN*par_TOR(rs,MAXTORRES+j)*vv(j)
        end if
      end do
    else
      do j=1,MAXTORRES
        if (vvz(j).gt.0) then
          call onetor_deriv(iz(3,vvz(j)),iz(2,vvz(j)),iz(1,vvz(j)),vvz(j),&
 &                             ctm,stm,ad1,ad2,ad3,ad4)
          dt0 = scale_TOR*par_TOR(rs,MAXTORRES+j)*2.0*vv(j)
          do k=1,3
            ca_f(iz(3,vvz(j)),k) = ca_f(iz(3,vvz(j)),k) - dt0*ad1(k)
            ca_f(iz(2,vvz(j)),k) = ca_f(iz(2,vvz(j)),k) - dt0*ad2(k)
            ca_f(iz(1,vvz(j)),k) = ca_f(iz(1,vvz(j)),k) - dt0*ad3(k)
            ca_f(vvz(j),k)       = ca_f(vvz(j),k) - dt0*ad4(k)
          end do
        end if
      end do
    end if
! the potential is a Gaussian well
  else if ((par_TOR2(rs).eq.2).OR.(par_TOR2(rs).eq.4)) then
    expte = exp(-torincr)
    evec(7) = evec(7) - scale_TOR*expte
    if ((fycxyz.eq.1).AND.(hmjam%isin(7).EQV..false.).AND.(hmjam%isin(11).EQV..false.)) then
      do j=1,MAXTORRES
!       note that the vvi-check will cycle out the sugar in nucleotides
        if ((vvz(j).gt.0).AND.(vvi(j).gt.0)) then
          dc_di(imol)%f(vvi(j)) = dc_di(imol)%f(vvi(j)) - &
 &              expte*scale_TOR*2.0*RADIAN*par_TOR(rs,MAXTORRES+j)*vv(j)
        end if
      end do
    else
      do j=1,MAXTORRES
        if (vvz(j).gt.0) then
          call onetor_deriv(iz(3,vvz(j)),iz(2,vvz(j)),iz(1,vvz(j)),vvz(j),&
 &                             ctm,stm,ad1,ad2,ad3,ad4)
          dt0 = expte*scale_TOR*par_TOR(rs,MAXTORRES+j)*2.0*vv(j)
          do k=1,3
            ca_f(iz(3,vvz(j)),k) = ca_f(iz(3,vvz(j)),k) - dt0*ad1(k)
            ca_f(iz(2,vvz(j)),k) = ca_f(iz(2,vvz(j)),k) - dt0*ad2(k)
            ca_f(iz(1,vvz(j)),k) = ca_f(iz(1,vvz(j)),k) - dt0*ad3(k)
            ca_f(vvz(j),k)       = ca_f(vvz(j),k) - dt0*ad4(k)
          end do
        end if
      end do
    end if
  else
    write(ilog,*) 'Fatal. Called force_torrs(...) with unknown mode. Offen&
 &ding value is ',par_TOR2(rs),' (residue ',rs,').'
    call fexit()
  end if
!
!  imol = molofrs(rs)
!  ff = fline(rs)
!  yy = yline(rs)
!
!  if ((par_TOR2(rs).eq.1).OR.(par_TOR2(rs).eq.3)) then
!!   the potential is harmonic
!    vv(1) = phi(rs) - par_TOR(rs,1)
!    vv(2) = psi(rs) - par_TOR(rs,2)
!!   periodic boundary conditions (so to speak)
!    do j=1,2
!      if (vv(j).lt.-180.0) vv(j) = vv(j) + 360.0
!      if (vv(j).gt.180.0) vv(j) = vv(j) - 360.0
!    end do
!    evec(7) = evec(7) + scale_TOR*&
! &     (par_TOR(rs,3)*vv(1)**2 + par_TOR(rs,4)*vv(2)**2)
!    if (fycxyz.eq.1) then
!      dc_di(imol)%f(fnr(rs)) = dc_di(imol)%f(fnr(rs)) -&
! &    scale_TOR*2.0*RADIAN*par_TOR(rs,3)*vv(1)
!      dc_di(imol)%f(ynr(rs)) = dc_di(imol)%f(ynr(rs)) -&
! &    scale_TOR*2.0*RADIAN*par_TOR(rs,4)*vv(2)
!    else if (fycxyz.eq.2) then
!      call onetor_deriv(iz(3,ff),iz(2,ff),iz(1,ff),ff,&
! &                             ctm,stm,ad1,ad2,ad3,ad4)
!      dt0 = scale_TOR*par_TOR(rs,3)*2.0*vv(1)
!      do j=1,3
!        cart_f(iz(3,ff),j) = cart_f(iz(3,ff),j) - dt0*ad1(j)
!        cart_f(iz(2,ff),j) = cart_f(iz(2,ff),j) - dt0*ad2(j)
!        cart_f(iz(1,ff),j) = cart_f(iz(1,ff),j) - dt0*ad3(j)
!        cart_f(ff,j)       = cart_f(ff,j) - dt0*ad4(j)
!      end do
!      call onetor_deriv(iz(3,yy),iz(2,yy),iz(1,yy),yy,&
! &                             ctm,stm,ad1,ad2,ad3,ad4)
!      dt0 = scale_TOR*par_TOR(rs,4)*2.0*vv(2)
!      do j=1,3
!        cart_f(iz(3,yy),j) = cart_f(iz(3,yy),j) - dt0*ad1(j)
!        cart_f(iz(2,yy),j) = cart_f(iz(2,yy),j) - dt0*ad2(j)
!        cart_f(iz(1,yy),j) = cart_f(iz(1,yy),j) - dt0*ad3(j)
!        cart_f(yy,j)       = cart_f(yy,j) - dt0*ad4(j)
!      end do
!    end if
!  else if ((par_TOR2(rs).eq.2).OR.(par_TOR2(rs).eq.4)) then
!!   the potential is a Gaussian well
!    vv(1) = phi(rs) - par_TOR(rs,1)
!    vv(2) = psi(rs) - par_TOR(rs,2)
!!   periodic boundary conditions (so to speak)
!    do j=1,2
!      if (vv(j).lt.-180.0) vv(j) = vv(j) + 360.0
!      if (vv(j).gt.180.0) vv(j) = vv(j) - 360.0
!    end do
!    expte = exp(-par_TOR(rs,3)*vv(1)**2 - par_TOR(rs,4)*vv(2)**2)
!    evec(7) = evec(7) - scale_TOR*expte
!    if (fycxyz.eq.1) then
!      dc_di(imol)%f(fnr(rs)) = dc_di(imol)%f(fnr(rs)) -&
! &      scale_TOR*expte*2.0*RADIAN*par_TOR(rs,3)*vv(1)
!      dc_di(imol)%f(ynr(rs)) = dc_di(imol)%f(ynr(rs)) -&
! &      scale_TOR*expte*2.0*RADIAN*par_TOR(rs,4)*vv(2)
!    else if (fycxyz.eq.2) then
!      call onetor_deriv(iz(3,ff),iz(2,ff),iz(1,ff),ff,&
! &                             ctm,stm,ad1,ad2,ad3,ad4)
!      dt0 = expte*scale_TOR*par_TOR(rs,3)*2.0*vv(1)
!      do j=1,3
!        cart_f(iz(3,ff),j) = cart_f(iz(3,ff),j) - dt0*ad1(j)
!        cart_f(iz(2,ff),j) = cart_f(iz(2,ff),j) - dt0*ad2(j)
!        cart_f(iz(1,ff),j) = cart_f(iz(1,ff),j) - dt0*ad3(j)
!        cart_f(ff,j)       = cart_f(ff,j) - dt0*ad4(j)
!      end do
!      call onetor_deriv(iz(3,yy),iz(2,yy),iz(1,yy),yy,&
! &                             ctm,stm,ad1,ad2,ad3,ad4)
!      dt0 = expte*scale_TOR*par_TOR(rs,4)*2.0*vv(2)
!      do j=1,3
!        cart_f(iz(3,yy),j) = cart_f(iz(3,yy),j) - dt0*ad1(j)
!        cart_f(iz(2,yy),j) = cart_f(iz(2,yy),j) - dt0*ad2(j)
!        cart_f(iz(1,yy),j) = cart_f(iz(1,yy),j) - dt0*ad3(j)
!        cart_f(yy,j)       = cart_f(yy,j) - dt0*ad4(j)
!      end do
!    end if
!  else
!    write(ilog,*) 'Fatal. Called force_torrs(rs) with unknown mode. &
! &Offending value is ',par_TOR2(rs),' (residue ',rs,').'
!    call fexit()
!  end if
!
end
!
!---------------------------------------------------------------------------
!
! for ZSEC the forces are a little more tedious because of the collective
! nature of z_alpha/z_beta, but nonetheless straightforward for torsional
! MD
!
subroutine force_zsec_gl(imol,evec,ca_f)
!
  use energies
  use torsn
  use molecule
  use forces
  use fyoc
  use math
  use system
  use zmatrix
  use atoms
  use wl
!
  implicit none
!
  integer i,j,imol,nress,ff,yy
  RTYPE evec(MAXENERGYTERMS),dis,vv(2),ww(2),uu(2),ca_f(n,3)
  RTYPE dddf,dddy,dddw(2,2),dddu(2,2),expte
  RTYPE dzadf(rsmol(imol,2)-rsmol(imol,1)+1)
  RTYPE dzady(rsmol(imol,2)-rsmol(imol,1)+1)
  RTYPE dzbdf(rsmol(imol,2)-rsmol(imol,1)+1)
  RTYPE dzbdy(rsmol(imol,2)-rsmol(imol,1)+1)
  RTYPE dt0,ad1(3),ad2(3),ad3(3),ad4(3),ctm,stm
!
! remember has to take care of unsupported polymer types
  if (do_tors(moltypid(imol)).EQV..false.) return
!
  z_alpha(imol) = 0.0
  z_beta(imol) = 0.0
!
  nress = rsmol(imol,2)-rsmol(imol,1)-1
!
  do i=rsmol(imol,1)+1,rsmol(imol,2)-1
!
    if ((fline(i).le.0).OR.(yline(i).le.0)) cycle
!   alpha distance
    vv(1) = phi(i) - par_ZSEC2(1)
    vv(2) = psi(i) - par_ZSEC2(2)
!   periodic boundary conditions (so to speak)
    do j=1,2
      if (vv(j).lt.-180.0) vv(j) = vv(j) + 360.0
      if (vv(j).gt.180.0) vv(j) = vv(j) - 360.0
    end do
!   now the distance
    dis = sqrt(vv(1)**2 + vv(2)**2)
    dddf = (1.0/dis)*vv(1)
    dddy = (1.0/dis)*vv(2)
    if (dis.le.par_ZSEC2(3)) then
!     this residue is fully alpha (force is zero)
      z_alpha(imol) = z_alpha(imol) + 1.0
      dzadf(i-rsmol(imol,1)) = 0.0
      dzady(i-rsmol(imol,1)) = 0.0
    else
!     this residue is contributing less than a full count (most likely ~0.0)
      ww(1) = par_ZSEC2(1) + (par_ZSEC2(3)/dis)*vv(1)
!     both dis and vv(1) depend on phi
      dddw(1,1) = -dddf*vv(1)*par_ZSEC2(3)/(dis*dis) +&
 &                                      (par_ZSEC2(3)/dis)
      dddw(1,2) = -dddy*vv(1)*par_ZSEC2(3)/(dis*dis) 
      ww(2) = par_ZSEC2(2) + (par_ZSEC2(3)/dis)*vv(2)
      dddw(2,2) = -dddy*vv(2)*par_ZSEC2(3)/(dis*dis) +&
 &                                      (par_ZSEC2(3)/dis)
      dddw(2,1) = -dddf*vv(2)*par_ZSEC2(3)/(dis*dis)
!     PBC #2
      do j=1,2
        if (ww(j).lt.-180.0) ww(j) = ww(j) + 360.0
        if (ww(j).gt.180.0) ww(j) = ww(j) - 360.0
      end do
      uu(1) = phi(i) - ww(1)
      dddu(1,1) = 1.0 - dddw(1,1)
      dddu(1,2) = -dddw(1,2)
      uu(2) = psi(i) - ww(2)
      dddu(2,2) = 1.0 - dddw(2,2)
      dddu(2,1) = -dddw(2,1)
!     PBC #3
      do j=1,2
        if (uu(j).lt.-180.0) uu(j) = uu(j) + 360.0
        if (uu(j).gt.180.0) uu(j) = uu(j) - 360.0
      end do
      expte = exp(-par_ZSEC2(4)*(uu(1)**2+uu(2)**2))
      z_alpha(imol) = z_alpha(imol) + expte
      dzadf(i-rsmol(imol,1)) = -par_ZSEC2(4)*expte*&
 &         (2.0*uu(1)*dddu(1,1) + 2.0*uu(2)*dddu(2,1))/(1.0*nress)
      dzady(i-rsmol(imol,1)) = -par_ZSEC2(4)*expte*&
 &         (2.0*uu(1)*dddu(1,2) + 2.0*uu(2)*dddu(2,2))/(1.0*nress)
    end if
!   now beta
    vv(1) = phi(i) - par_ZSEC2(5)
    vv(2) = psi(i) - par_ZSEC2(6)
    do j=1,2
      if (vv(j).lt.-180.0) vv(j) = vv(j) + 360.0
      if (vv(j).gt.180.0) vv(j) = vv(j) - 360.0
    end do
    dis = sqrt(vv(1)**2 + vv(2)**2)
    dddf = (1.0/dis)*vv(1)
    dddy = (1.0/dis)*vv(2)
    if (dis.le.par_ZSEC2(7)) then
      z_beta(imol) = z_beta(imol) + 1.0
      dzbdf(i-rsmol(imol,1)) = 0.0
      dzbdy(i-rsmol(imol,1)) = 0.0
    else
      ww(1) = par_ZSEC2(5) + (par_ZSEC2(7)/dis)*vv(1)
!     both dis and vv(1) depend on phi
      dddw(1,1) = -dddf*vv(1)*par_ZSEC2(7)/(dis*dis) +&
 &                                      (par_ZSEC2(7)/dis)
      dddw(1,2) = -dddy*vv(1)*par_ZSEC2(7)/(dis*dis)
      ww(2) = par_ZSEC2(6) + (par_ZSEC2(7)/dis)*vv(2)
      dddw(2,2) = -dddy*vv(2)*par_ZSEC2(7)/(dis*dis) +&
 &                                      (par_ZSEC2(7)/dis)
      dddw(2,1) = -dddf*vv(2)*par_ZSEC2(7)/(dis*dis)
      do j=1,2
        if (ww(j).lt.-180.0) ww(j) = ww(j) + 360.0
        if (ww(j).gt.180.0) ww(j) = ww(j) - 360.0
      end do
      uu(1) = phi(i) - ww(1)
      dddu(1,1) = 1.0 - dddw(1,1)
      dddu(1,2) = -dddw(1,2)
      uu(2) = psi(i) - ww(2)
      dddu(2,2) = 1.0 - dddw(2,2)
      dddu(2,1) = -dddw(2,1)
      do j=1,2
        if (uu(j).lt.-180.0) uu(j) = uu(j) + 360.0
        if (uu(j).gt.180.0) uu(j) = uu(j) - 360.0
      end do
      expte = exp(-par_ZSEC2(8)*(uu(1)**2+uu(2)**2))
      z_beta(imol) = z_beta(imol) + expte
      dzbdf(i-rsmol(imol,1)) = -par_ZSEC2(8)*expte*&
 &         (2.0*uu(1)*dddu(1,1) + 2.0*uu(2)*dddu(2,1))/(1.0*nress)
      dzbdy(i-rsmol(imol,1)) = -par_ZSEC2(8)*expte*&
 &         (2.0*uu(1)*dddu(1,2) + 2.0*uu(2)*dddu(2,2))/(1.0*nress)
    end if
  end do
!
  z_alpha(imol) = z_alpha(imol)/(1.0*nress)
  z_beta(imol) = z_beta(imol)/(1.0*nress)
!
  evec(8) = evec(8) + scale_ZSEC*&
 &            (par_ZSEC(2)*(z_alpha(imol) - par_ZSEC(1))**2 +&
 &             par_ZSEC(4)*(z_beta(imol) - par_ZSEC(3))**2)
!
  do i=rsmol(imol,1)+1,rsmol(imol,2)-1
!
    if ((fline(i).le.0).OR.(yline(i).le.0)) cycle
    if ((fnr(i).le.0).OR.(ynr(i).le.0)) cycle
    if ((fycxyz.eq.1).AND.(hmjam%isin(8).EQV..false.).AND.(hmjam%isin(11).EQV..false.)) then
      dc_di(imol)%f(fnr(i)) = dc_di(imol)%f(fnr(i)) - scale_ZSEC*RADIAN*&
 &   (2.0*par_ZSEC(2)*(z_alpha(imol) - par_ZSEC(1))*dzadf(i-rsmol(imol,1)) +&
 &    2.0*par_ZSEC(4)*(z_beta(imol)  - par_ZSEC(3))*dzbdf(i-rsmol(imol,1)))

      dc_di(imol)%f(ynr(i)) = dc_di(imol)%f(ynr(i)) - scale_ZSEC*RADIAN*&
 &   (2.0*par_ZSEC(2)*(z_alpha(imol) - par_ZSEC(1))*dzady(i-rsmol(imol,1)) +&
 &    2.0*par_ZSEC(4)*(z_beta(imol)  - par_ZSEC(3))*dzbdy(i-rsmol(imol,1)))
    else 
!
      ff = fline(i)
      call onetor_deriv(iz(3,ff),iz(2,ff),iz(1,ff),ff,&
 &                             ctm,stm,ad1,ad2,ad3,ad4)
      dt0 = scale_ZSEC*&
 &   (2.0*par_ZSEC(2)*(z_alpha(imol) - par_ZSEC(1))*&
 &                dzadf(i-rsmol(imol,1)) +&
 &    2.0*par_ZSEC(4)*(z_beta(imol) - par_ZSEC(3))*&
 &                dzbdf(i-rsmol(imol,1)))
      do j=1,3
        ca_f(iz(3,ff),j) = ca_f(iz(3,ff),j) - dt0*ad1(j)
        ca_f(iz(2,ff),j) = ca_f(iz(2,ff),j) - dt0*ad2(j)
        ca_f(iz(1,ff),j) = ca_f(iz(1,ff),j) - dt0*ad3(j)
        ca_f(ff,j)       = ca_f(ff,j) - dt0*ad4(j)
      end do
      yy = yline(i)
      call onetor_deriv(iz(3,yy),iz(2,yy),iz(1,yy),yy,&
 &                             ctm,stm,ad1,ad2,ad3,ad4)
      dt0 = scale_ZSEC*&
 &   (2.0*par_ZSEC(2)*(z_alpha(imol) - par_ZSEC(1))*&
 &                dzady(i-rsmol(imol,1)) +&
 &    2.0*par_ZSEC(4)*(z_beta(imol) - par_ZSEC(3))*&
 &                dzbdy(i-rsmol(imol,1)))
      do j=1,3
        ca_f(iz(3,yy),j) = ca_f(iz(3,yy),j) - dt0*ad1(j)
        ca_f(iz(2,yy),j) = ca_f(iz(2,yy),j) - dt0*ad2(j)
        ca_f(iz(1,yy),j) = ca_f(iz(1,yy),j) - dt0*ad3(j)
        ca_f(yy,j)       = ca_f(yy,j) - dt0*ad4(j)
      end do
!
    end if
!
  end do
!
end
!
!-----------------------------------------------------------------------------
!
  subroutine gradient_test(epsil1,epsil2,epsil3)
!
  use forces
  use iounit
  use atoms
  use molecule
  use sequen
  use fyoc
  use movesets
  use math
  use energies
  use system
  use mpistuff
  use cutoffs
  use ems
!
  implicit none
!
  integer freeunit,i,j,k,ttc
  integer ingt1,ingt2,ii,jj
  RTYPE blas(1000),epsil1,epsil2,epsil3,force3,energy3
  RTYPE etd1(MAXENERGYTERMS),etd2(MAXENERGYTERMS),etd3(MAXENERGYTERMS)
  character(MAXSTRLEN) fn
  logical atrue,afalse,exists,rspmatthere
#ifdef ENABLE_MPI
  integer lext2
  character(3) xpont
#endif
!
  atrue = .true.
  afalse = .false.
  rspmatthere = allocated(rsp_mat)
!
 4467 format(180(g14.8,1x))
!
  if (use_dyn.EQV..false.) then
    write(ilog,*) 'Fatal. Called gradient_test(...) while forces are not in use. This is a bug.'
    call fexit()
  end if
!
#ifdef ENABLE_MPI
  lext2 = 3
  call int2str(myrank,xpont,lext2)
  fn =  'N_'//xpont(1:lext2)//'_NUM_GRAD_TEST_XYZ.dat'
#else
  fn = 'NUM_GRAD_TEST_XYZ.dat'
#endif
!
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    ingt1 = freeunit()
    open(unit=ingt1,file=fn(ii:jj),status='old')
    close(unit=ingt1,status='delete')
  end if
  ingt1=freeunit()
  open(unit=ingt1,file=fn(ii:jj),status='new')
!
  blas(1) = epsil1
!
! make sure everything is up-to-date RBC-wise
  do i=1,nmol
    call update_rigidm(i)
  end do
!
! make sure energy3 doesn't seg-fault
  if (rspmatthere.EQV..false.) then
    allocate(rsp_mat(nseq,nseq))
  end if
!
! first Cartesian derivatives
  ttc = 0
  do i=1,n
    if (use_POLY.EQV..true.) call update_rigid(molofrs(atmres(i)))
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(9) = force3(etd1,etd2,etd3,atrue)
    blas(12) = energy3(etd1,atrue)
    x(i) = x(i) + 0.5*blas(1)
    call genzmat(molofrs(atmres(i)))
    call zmatfyc()
    if (use_POLY.EQV..true.) call update_rigid(molofrs(atmres(i)))
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(10) = energy3(etd1,atrue)
    x(i) = x(i) - 1.0*blas(1)
    call genzmat(molofrs(atmres(i)))
    call zmatfyc()
    if (use_POLY.EQV..true.) call update_rigid(molofrs(atmres(i)))
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(11) = energy3(etd1,atrue)
    x(i) = x(i) + 0.5*blas(1)
    call genzmat(molofrs(atmres(i)))
    call zmatfyc()
    write(ingt1,4467) (blas(11)-blas(10))/blas(1),cart_f(i,1),blas(12),blas(9)
    ttc = ttc + 1
    y(i) = y(i) + 0.5*blas(1)
    call genzmat(molofrs(atmres(i)))
    call zmatfyc()
    if (use_POLY.EQV..true.) call update_rigid(molofrs(atmres(i)))
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(10) = energy3(etd1,atrue)
    y(i) = y(i) - 1.0*blas(1)
    call genzmat(molofrs(atmres(i)))
    call zmatfyc()
    if (use_POLY.EQV..true.) call update_rigid(molofrs(atmres(i)))
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(11) = energy3(etd1,atrue)
    y(i) = y(i) + 0.5*blas(1)
    call genzmat(molofrs(atmres(i)))
    call zmatfyc()
    write(ingt1,4467) (blas(11)-blas(10))/blas(1),cart_f(i,2),blas(12),blas(9)
    ttc = ttc + 1
    z(i) = z(i) + 0.5*blas(1)
    call genzmat(molofrs(atmres(i)))
    call zmatfyc()
    if (use_POLY.EQV..true.) call update_rigid(molofrs(atmres(i)))
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(10) = energy3(etd1,atrue)
    z(i) = z(i) - 1.0*blas(1)
    call genzmat(molofrs(atmres(i)))
    call zmatfyc()
    if (use_POLY.EQV..true.) call update_rigid(molofrs(atmres(i)))
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(11) = energy3(etd1,atrue)
    z(i) = z(i) + 0.5*blas(1)
    call genzmat(molofrs(atmres(i)))
    call zmatfyc()
    write(ingt1,4467) (blas(11)-blas(10))/blas(1),cart_f(i,3),blas(12),blas(9)
  end do
  close(unit=ingt1)
!
! exit if Cartesian
  if (fycxyz.ne.1) then
    if (rspmatthere.EQV..false.) deallocate(rsp_mat)
    return
  end if
!
! test internal forces
!
#ifdef ENABLE_MPI
  lext2 = 3
  call int2str(myrank,xpont,lext2)
  fn =  'N_'//xpont(1:lext2)//'_NUM_GRAD_TEST_INT.dat'
#else
  fn = 'NUM_GRAD_TEST_INT.dat'
#endif
!
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    ingt2 = freeunit()
    open(unit=ingt2,file=fn(ii:jj),status='old')
    close(unit=ingt2,status='delete')
  end if
  ingt2=freeunit()
  open(unit=ingt2,file=fn(ii:jj),status='new')
!
  do i=1,nmol
    call update_rigidm(i)
  end do
  do i=1,nmol
    ttc = 1
    blas(1) = epsil2
    if (use_POLY.EQV..true.) call update_rigid(i)
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(9) = force3(etd1,etd2,etd3,atrue)
    dc_di(i)%f(:) = 0.0
    call cart2int_f(skip_frz)
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(12) = energy3(etd1,atrue)
    do j=atmol(i,1),atmol(i,2)
      x(j) = x(j) + 0.5*blas(1)
    end do
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(10) = energy3(etd1,atrue)
    do j=atmol(i,1),atmol(i,2)
      x(j) = x(j) - 1.0*blas(1)
    end do
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(11) = energy3(etd1,atrue)
    do j=atmol(i,1),atmol(i,2)
      x(j) = x(j) + 0.5*blas(1)
    end do
    write(ingt1,4467) (blas(11)-blas(10))/blas(1),dc_di(i)%f(ttc),blas(12),blas(9),x(atmol(i,1))
    ttc = ttc + 1
    do j=atmol(i,1),atmol(i,2)
      y(j) = y(j) + 0.5*blas(1)
    end do
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(10) = energy3(etd1,atrue)
    do j=atmol(i,1),atmol(i,2)
      y(j) = y(j) - 1.0*blas(1)
    end do
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(11) = energy3(etd1,atrue)
    do j=atmol(i,1),atmol(i,2)
      y(j) = y(j) + 0.5*blas(1)
    end do
    write(ingt2,4467) (blas(11)-blas(10))/blas(1),dc_di(i)%f(ttc),blas(12),blas(9),y(atmol(i,1))
    ttc = ttc + 1
    do j=atmol(i,1),atmol(i,2)
      z(j) = z(j) + 0.5*blas(1)
    end do
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(10) = energy3(etd1,atrue)
    do j=atmol(i,1),atmol(i,2)
      z(j) = z(j) - 1.0*blas(1)
    end do
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(11) = energy3(etd1,atrue)
    do j=atmol(i,1),atmol(i,2)
      z(j) = z(j) + 0.5*blas(1)
    end do
    write(ingt2,4467) (blas(11)-blas(10))/blas(1),dc_di(i)%f(ttc),blas(12),blas(9),z(atmol(i,1))
    if (atmol(i,2).eq.atmol(i,1)) cycle
    ttc = ttc + 1
    blas(1) = epsil3
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(12) = energy3(etd1,atrue)
    cur_rot(1) = 0.5*blas(1)
    cur_rot(2) = 0.0
    cur_rot(3) = 0.0
    call rotxyzm(i,1)
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(10) = energy3(etd1,atrue)
    cur_rot(1) = -1.0*blas(1)
    call rotxyzm(i,1)
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(11) = energy3(etd1,atrue)
    cur_rot(1) = 0.5*blas(1)
    call rotxyzm(i,1)
    write(ingt2,4467) (blas(11)-blas(10))/blas(1),dc_di(i)%f(ttc),blas(12),blas(9),x(atmol(i,2))
    ttc = ttc + 1
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(12) = energy3(etd1,atrue)
    cur_rot(1) = 0.0
    cur_rot(2) = 0.5*blas(1)
    cur_rot(3) = 0.0
    call rotxyzm(i,1)
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(10) = energy3(etd1,atrue)
    cur_rot(2) = -1.0*blas(1)
    call rotxyzm(i,1)
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(11) = energy3(etd1,atrue)
    cur_rot(2) = 0.5*blas(1)
    call rotxyzm(i,1)
    write(ingt2,4467) (blas(11)-blas(10))/blas(1),dc_di(i)%f(ttc),blas(12),blas(9),x(atmol(i,2))
    ttc = ttc + 1
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(12) = energy3(etd1,atrue)
    cur_rot(1) = 0.0
    cur_rot(2) = 0.0
    cur_rot(3) = 0.5*blas(1)
    call rotxyzm(i,1)
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(10) = energy3(etd1,atrue)
    cur_rot(3) = -1.0*blas(1)
    call rotxyzm(i,1)
    if (use_EMICRO.EQV..true.) curmassm = .false.
    blas(11) = energy3(etd1,atrue)
    cur_rot(3) = 0.5*blas(1)
    call rotxyzm(i,1)
    write(ingt2,4467) (blas(11)-blas(10))/blas(1),dc_di(i)%f(ttc),blas(12),blas(9),x(atmol(i,2))
    ttc = 6
    do j=rsmol(i,1),rsmol(i,2)
      if (wline(j).gt.0) then
        ttc = ttc + 1
        dc_di(i)%f(:) = 0.0
        dc_di(i)%incr(:) = 0.0
        call fyczmat()
        call makexyz_formol(i)
        blas(8) = omega(j)
        if (use_POLY.EQV..true.) call update_rigid(i)
        if (use_EMICRO.EQV..true.) curmassm = .false.
        blas(9) = force3(etd1,etd2,etd3,atrue)
        blas(12) = energy3(etd1,atrue)
        call cart2int_f(skip_frz)
        omega(j) = omega(j) + 0.5*blas(1)
        call fyczmat()
        if (dc_di(i)%align(ttc).EQV..true.) then
          dc_di(i)%incr(ttc) = 0.5*blas(1)
          call IMD_prealign(i,afalse)
        end if
        call makexyz_formol(i)
        if (use_POLY.EQV..true.) call update_rigid(i)
        if (use_EMICRO.EQV..true.) curmassm = .false.
        blas(10) = energy3(etd1,atrue)
        omega(j) = omega(j) - 1.0*blas(1)
        call fyczmat()
        if (dc_di(i)%align(ttc).EQV..true.) then
          dc_di(i)%incr(ttc) = -1.0*blas(1)
          call IMD_prealign(i,afalse)
        end if
        call makexyz_formol(i)
        if (use_POLY.EQV..true.) call update_rigid(i)
        if (use_EMICRO.EQV..true.) curmassm = .false.
        blas(11) = energy3(etd1,atrue)
        omega(j) = blas(8)
        call fyczmat()
        if (dc_di(i)%align(ttc).EQV..true.) then
          dc_di(i)%incr(ttc) = 0.5*blas(1)
          call IMD_prealign(i,afalse)
        end if
        call makexyz_formol(i)
        write(ingt2,4467) RADIAN*(blas(11)-blas(10))/blas(1),dc_di(i)%f(ttc),blas(12),blas(9),blas(8)
      end if
      if ((fline(j).gt.0).AND.(seqflag(j).ne.5)) then
        ttc = ttc + 1
        dc_di(i)%f(:) = 0.0
        dc_di(i)%incr(:) = 0.0
        blas(8) = phi(j)
        if (use_POLY.EQV..true.) call update_rigid(i)
        if (use_EMICRO.EQV..true.) curmassm = .false.
        blas(9) = force3(etd1,etd2,etd3,atrue)
        blas(12) = energy3(etd1,atrue)
        call cart2int_f(skip_frz)
        phi(j) = phi(j) + 0.5*blas(1)
        call fyczmat()
        if (dc_di(i)%align(ttc).EQV..true.) then
          dc_di(i)%incr(ttc) = 0.5*blas(1)
          call IMD_prealign(i,afalse)
        end if
        call makexyz_formol(i)
        if (use_POLY.EQV..true.) call update_rigid(i)
        if (use_EMICRO.EQV..true.) curmassm = .false.
        blas(10) = energy3(etd1,atrue)
        phi(j) = phi(j) - 1.0*blas(1)  
        call fyczmat()
        if (dc_di(i)%align(ttc).EQV..true.) then
          dc_di(i)%incr(ttc) = -1.0*blas(1)
          call IMD_prealign(i,afalse)
        end if
        call makexyz_formol(i)
        if (use_POLY.EQV..true.) call update_rigid(i)
        if (use_EMICRO.EQV..true.) curmassm = .false.
        blas(11) = energy3(etd1,atrue)
        phi(j) = blas(8)
        call fyczmat()
        if (dc_di(i)%align(ttc).EQV..true.) then
          dc_di(i)%incr(ttc) = 0.5*blas(1)
          call IMD_prealign(i,afalse)
        end if
        call makexyz_formol(i)
        write(ingt2,4467) RADIAN*(blas(11)-blas(10))/blas(1),dc_di(i)%f(ttc),blas(12),blas(9),blas(8)
      else if (seqflag(j).eq.5) then
        ttc = ttc + 1
      end if
!
      if (yline(j).gt.0) then
        ttc = ttc + 1
        dc_di(i)%f(:) = 0.0
        dc_di(i)%incr(:) = 0.0
        blas(8) = psi(j)
        if (use_POLY.EQV..true.) call update_rigid(i)
        if (use_EMICRO.EQV..true.) curmassm = .false.
        blas(9) = force3(etd1,etd2,etd3,atrue)
        blas(12) = energy3(etd1,atrue)
        call cart2int_f(skip_frz)
        psi(j) = psi(j) + 0.5*blas(1)
        call fyczmat()
        if (dc_di(i)%align(ttc).EQV..true.) then
          dc_di(i)%incr(ttc) = 0.5*blas(1)
          call IMD_prealign(i,afalse)
        end if
        call makexyz_formol(i)
        if (use_POLY.EQV..true.) call update_rigid(i)
        if (use_EMICRO.EQV..true.) curmassm = .false.
        blas(10) = energy3(etd1,atrue)
        psi(j) = psi(j) - 1.0*blas(1)  
        call fyczmat()
        if (dc_di(i)%align(ttc).EQV..true.) then
          dc_di(i)%incr(ttc) = -1.0*blas(1)
          call IMD_prealign(i,afalse)
        end if
        call makexyz_formol(i)
        if (use_POLY.EQV..true.) call update_rigid(i)
        if (use_EMICRO.EQV..true.) curmassm = .false.
        blas(11) = energy3(etd1,atrue)
        psi(j) = blas(8)
        call fyczmat()
        if (dc_di(i)%align(ttc).EQV..true.) then
          dc_di(i)%incr(ttc) = 0.5*blas(1)
          call IMD_prealign(i,afalse)
        end if
        call makexyz_formol(i)
        write(ingt2,4467) RADIAN*(blas(11)-blas(10))/blas(1),dc_di(i)%f(ttc),blas(12),blas(9),blas(8)
      end if
      do k=1,nnucs(j)
        ttc = ttc + 1
        dc_di(i)%f(:) = 0.0
        dc_di(i)%incr(:) = 0.0
        blas(8) = nucs(k,j)
        if (use_POLY.EQV..true.) call update_rigid(i)
        if (use_EMICRO.EQV..true.) curmassm = .false.
        blas(9) = force3(etd1,etd2,etd3,atrue)
        blas(12) = energy3(etd1,atrue)
        call cart2int_f(skip_frz)
        nucs(k,j) = nucs(k,j) + 0.5*blas(1)
        call fyczmat()
        if (dc_di(i)%align(ttc).EQV..true.) then
          dc_di(i)%incr(ttc) = 0.5*blas(1)
          call IMD_prealign(i,afalse)
        end if
        call makexyz_formol(i)
        if (use_POLY.EQV..true.) call update_rigid(i)
        if (use_EMICRO.EQV..true.) curmassm = .false.
        blas(10) = energy3(etd1,atrue)
        nucs(k,j) = nucs(k,j) - 1.0*blas(1)  
        call fyczmat()
        if (dc_di(i)%align(ttc).EQV..true.) then
          dc_di(i)%incr(ttc) = -1.0*blas(1)
          call IMD_prealign(i,afalse)
        end if
        call makexyz_formol(i)
        if (use_POLY.EQV..true.) call update_rigid(i)
        if (use_EMICRO.EQV..true.) curmassm = .false.
        blas(11) = energy3(etd1,atrue)
        nucs(k,j) = blas(8)
        call fyczmat()
        if (dc_di(i)%align(ttc).EQV..true.) then
          dc_di(i)%incr(ttc) = 0.5*blas(1)
          call IMD_prealign(i,afalse)
        end if
        call makexyz_formol(i)
        write(ingt2,4467) RADIAN*(blas(11)-blas(10))/blas(1),dc_di(i)%f(ttc),blas(12),blas(9),blas(8)
      end do
      do k=1,nchi(j)
        call fyczmat()
        call makexyz_formol(i)
        ttc = ttc + 1
        dc_di(i)%f(:) = 0.0
        dc_di(i)%incr(:) = 0.0
        blas(8) = chi(k,j)
        if (use_POLY.EQV..true.) call update_rigid(i)
        if (use_EMICRO.EQV..true.) curmassm = .false.
        blas(9) = force3(etd1,etd2,etd3,atrue)
        blas(12) = energy3(etd1,atrue)
        call cart2int_f(skip_frz)
        chi(k,j) = chi(k,j) + 0.5*blas(1)
        call fyczmat()
        if (dc_di(i)%align(ttc).EQV..true.) then
          dc_di(i)%incr(ttc) = 0.5*blas(1)
          call IMD_prealign(i,afalse)
        end if
        call makexyz_formol(i)
        if (use_POLY.EQV..true.) call update_rigid(i)
        if (use_EMICRO.EQV..true.) curmassm = .false.
        blas(10) = energy3(etd1,atrue)
        chi(k,j) = chi(k,j) - 1.0*blas(1)  
        call fyczmat()
        if (dc_di(i)%align(ttc).EQV..true.) then
          dc_di(i)%incr(ttc) = -1.0*blas(1)
          call IMD_prealign(i,afalse)
        end if
        call makexyz_formol(i)
        if (use_POLY.EQV..true.) call update_rigid(i)
        if (use_EMICRO.EQV..true.) curmassm = .false.
        blas(11) = energy3(etd1,atrue)
        chi(k,j) = blas(8)
        call fyczmat()
        if (dc_di(i)%align(ttc).EQV..true.) then
          dc_di(i)%incr(ttc) = 0.5*blas(1)
          call IMD_prealign(i,afalse)
        end if
        call makexyz_formol(i)
        write(ingt2,4467) RADIAN*(blas(11)-blas(10))/blas(1),dc_di(i)%f(ttc),blas(12),blas(9),blas(8)
      end do
    end do
    dc_di(i)%f(:) = 0.0
    dc_di(i)%im(:) = 0.0
    dc_di(i)%incr(:) = 0.0
  end do
  close(unit=ingt2)
!
  if (rspmatthere.EQV..false.) deallocate(rsp_mat)
  return
!
end
!
!------------------------------------------------------------------------------------------
!
