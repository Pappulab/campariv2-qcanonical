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
!--------------------------------------------------------------------
!
subroutine setupmcgrid()
!
  use iounit
  use sequen
  use atoms
  use cutoffs
  use system
  use mcgrid
!
  implicit none
!
  integer i,j,which,gridnavi
  RTYPE pvd2,cutplusmargin2
  RTYPE refp(3)
!
! first use boundary condition to determine grid origin and deltas
! periodic BC (PBC)
  cutplusmargin2 = (mcel_cutoff + 2.0*mrrd)**2.0
  if (bnd_type.eq.1) then
!   rectangular box
    if (bnd_shape.eq.1) then
      do i=1,3
        grid%origin(i) = bnd_params(3+i)
        grid%deltas(i) = bnd_params(i)/(1.0*grid%dim(i))
      end do
!   cylinder
    else if (bnd_shape.eq.3) then
      grid%origin(1:2) = bnd_params(1:2) - 1.5*bnd_params(4)
      grid%deltas(1:2) = 3.0*bnd_params(4)/(1.0*grid%dim(1:2))
      grid%origin(3) = bnd_params(3) - 0.5*bnd_params(6)
      grid%deltas(3) = bnd_params(6)/(1.0*grid%dim(3))
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in set&
 &upmcgrid() (code # ',bnd_shape,').'
      call fexit()
    end if
! all-wall BCs
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   rectangular box
    if (bnd_shape.eq.1) then
      do i=1,3
        grid%origin(i) = bnd_params(3+i) - 0.25*bnd_params(i)
        grid%deltas(i) = (1.5*bnd_params(i))/(1.0*grid%dim(i))
      end do
!   spherical system
    else if (bnd_shape.eq.2) then
      do i=1,3
        grid%origin(i) = bnd_params(i)-1.5*bnd_params(4)
        grid%deltas(i) = 3.0*bnd_params(4)/(1.0*grid%dim(i))
      end do
!   cylinder
    else if (bnd_shape.eq.3) then
      grid%origin(1:2) = bnd_params(1:2) - 1.5*bnd_params(4)
      grid%deltas(1:2) = 3.0*bnd_params(4)/(1.0*grid%dim(1:2))
      grid%origin(3) = bnd_params(3) - 0.75*bnd_params(6)
      grid%deltas(3) = 1.5*bnd_params(6)/(1.0*grid%dim(3))
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in set&
 &upmcgrid() (code # ',bnd_shape,').'
      call fexit()
    end if
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition&
 & in setupmcgrid() (code # ',bnd_type,').'
    call fexit()
  end if
  grid%diags(1,1,1) = 0.0
  grid%diags(2,2,2) = sqrt(grid%deltas(1)*grid%deltas(1)&
 &                     + grid%deltas(2)*grid%deltas(2)&
 &                     + grid%deltas(3)*grid%deltas(3))
  grid%diags(2,2,1) = sqrt(grid%deltas(1)*grid%deltas(1)&
 &                     + grid%deltas(2)*grid%deltas(2))
  grid%diags(2,1,2) = sqrt(grid%deltas(1)*grid%deltas(1)&
 &                     + grid%deltas(3)*grid%deltas(3))
  grid%diags(1,2,2) = sqrt(grid%deltas(2)*grid%deltas(2)&
 &                     + grid%deltas(3)*grid%deltas(3))
  grid%diags(2,1,1) = grid%deltas(1)
  grid%diags(1,2,1) = grid%deltas(2)
  grid%diags(1,1,2) = grid%deltas(3)
  grid%pnts = grid%dim(1)*grid%dim(2)*grid%dim(3)
!
! now major allocation for NB lists (these arrays easily get very large)
  allocate(grid%nblist(grid%pnts+1,grid%maxnbs))
  allocate(grid%gpnblist(grid%pnts+1,grid%maxgpnbs))
  allocate(grid%gpnbmd(grid%pnts+1,grid%maxgpnbs))
  if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    allocate(grid%b_gpnblist(grid%pnts+1,grid%maxgpnbs))
    allocate(grid%b_gpnbmd(grid%pnts+1,grid%maxgpnbs))
  end if
!
  do i=1,grid%pnts
    do j=i+1,grid%pnts
      call gridmindiff(i,j,pvd2)
      if (pvd2.lt.cutplusmargin2) then
        grid%gpnbnr(i) = grid%gpnbnr(i) + 1
        if (grid%gpnbnr(i).gt.grid%maxgpnbs) then
          write(ilog,*) 'Fatal. Exceeded static neighbor number for grid points in setupmcgrid().&
 & Increase the value for FMCSC_GRIDMAXGPNB or choose a coarser grid.'
          call fexit()
        end if
        grid%gpnblist(i,grid%gpnbnr(i)) = j
        grid%gpnbmd(i,grid%gpnbnr(i)) = floor(sqrt(pvd2))
!       for MC, we need the complete double-mapped list
        if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
          grid%b_gpnbnr(j) = grid%b_gpnbnr(j) + 1
          if (grid%b_gpnbnr(j).gt.grid%maxgpnbs) then
            write(ilog,*) 'Fatal. Exceeded static neighbor number for grid points in setupmcgrid().&
 & Increase the value for FMCSC_GRIDMAXGPNB or choose a coarser grid.'
            call fexit()
          end if
          grid%b_gpnblist(j,grid%b_gpnbnr(j)) = i
          grid%b_gpnbmd(j,grid%b_gpnbnr(j)) = floor(sqrt(pvd2))
        end if    
      end if
    end do
  end do
!
! assign the various residues to 'their' grid points (via C-alpha or other reference atom)
!
  do i=1,nseq
    refp(1) = x(refat(i))
    refp(2) = y(refat(i))
    refp(3) = z(refat(i))
    which = gridnavi(refp)
    grid%nbnr(which) = grid%nbnr(which) + 1
    if (grid%nbnr(which).gt.grid%maxnbs) then
      write(ilog,*) 'Fatal. Exceeded static maximum for groups associated to a single grid point in setupmcgrid().&
 & Increase the value for FMCSC_GRIDMAXRSNB or choose a finer grid.'
      call fexit()
    end if
    grid%nblist(which,grid%nbnr(which)) = i
    grid%resgp(i) = which
  end do
!
end
!
!------------------------------------------------------------------------
!
function gridnavi(refp)
!
  use iounit
  use mcgrid
  use system
  use movesets
!
  implicit none
!
  integer xi,yi,zi,gridnavi,k(3),i
  RTYPE refp(3)
! see above
!
! get the grid indexing
!
  xi = floor(((refp(1)-grid%origin(1))/grid%deltas(1))) + 1
  yi = floor(((refp(2)-grid%origin(2))/grid%deltas(2))) + 1
  zi = floor(((refp(3)-grid%origin(3))/grid%deltas(3))) + 1
!
! periodic BC (PBC) (note that global shifts are molecule-wise and that refat(rs) may consequently
! hang "off" the grid at all times) -> shift back
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        k(i) = 0
        if (refp(i).lt.bnd_params(3+i)) then
          k(i) = floor((bnd_params(3+i)-refp(i))/bnd_params(i)) + 1
        else if(refp(i).gt.(bnd_params(3+i)+bnd_params(i))) then
          k(i) = floor((refp(i)-(bnd_params(3+i)+bnd_params(i)))/&
 &                                             bnd_params(i)) + 1
          k(i) = -k(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      k(:) = 0
      if (refp(3).lt.(bnd_params(3)-0.5*bnd_params(6))) then
        k(3) = floor((bnd_params(3)-0.5*bnd_params(6)-refp(3))/bnd_params(6)) + 1
      else if(refp(3).gt.(bnd_params(3)+0.5*bnd_params(6))) then
        k(3) = floor((refp(3)-(bnd_params(3)+0.5*bnd_params(6)))/bnd_params(6)) + 1
        k(3) = -k(3)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in gri&
 &dnavi() (code # ',bnd_shape,').'
      call fexit()
    end if
    xi = xi + k(1)*grid%dim(1)
    yi = yi + k(2)*grid%dim(2)
    zi = zi + k(3)*grid%dim(3)
  end if
!
 45 format(' x: ',g12.4,' y: ',g12.4,' z: ',g12.4)
  if ((xi.le.0).OR.(xi.gt.grid%dim(1))) then
    xi = max(xi,1)
    xi = min(xi,grid%dim(1))
!   in MC this happens -> lump into outermost grid-point
    if ((use_dyn.EQV..true.).OR.((bnd_type.eq.1).AND.(bnd_shape.eq.1))) then
      if ((bnd_type.eq.1).AND.(bnd_shape.eq.1)) then
        write(ilog,*) 'Warning. Position is off grid:'
        write(ilog,45) refp(1),refp(2),refp(3)
      else if (.NOT.(((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)).AND.(in_dyncyc.EQV..false.))) then
        write(ilog,*) 'Warning. Position is off grid:'
        write(ilog,45) refp(1),refp(2),refp(3)
        write(ilog,*) 'In non-periodic, soft boundaries, this will usually indicate that&
 & the system is poorly equilibrated and/or that the simulation has become unstable.'
      end if
    end if
  end if
  if ((yi.le.0).OR.(yi.gt.grid%dim(2))) then
    yi = max(yi,1)
    yi = min(yi,grid%dim(2))
!   in MC this happens -> lump into outermost grid-point
    if ((use_dyn.EQV..true.).OR.((bnd_type.eq.1).AND.(bnd_shape.eq.1))) then
      if ((bnd_type.eq.1).AND.(bnd_shape.eq.1)) then
        write(ilog,*) 'Warning. Position is off grid:'
        write(ilog,45) refp(1),refp(2),refp(3)
      else if (.NOT.(((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)).AND.(in_dyncyc.EQV..false.))) then
        write(ilog,*) 'Warning. Position is off grid:'
        write(ilog,45) refp(1),refp(2),refp(3)
        write(ilog,*) 'In non-periodic, soft boundaries, this will usually indicate that&
 & the system is poorly equilibrated and/or that the simulation has become unstable.'
      end if
    end if
  end if
  if ((zi.le.0).OR.(zi.gt.grid%dim(3))) then
    zi = max(zi,1)
    zi = min(zi,grid%dim(3))
!   in MC this happens -> lump into outermost grid-point
    if ((use_dyn.EQV..true.).OR.(bnd_type.eq.1)) then
      if (bnd_type.eq.1) then
        write(ilog,*) 'Warning. Position is off grid:'
        write(ilog,45) refp(1),refp(2),refp(3)
      else if (.NOT.(((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)).AND.(in_dyncyc.EQV..false.))) then
        write(ilog,*) 'Warning. Position is off grid:'
        write(ilog,45) refp(1),refp(2),refp(3)
        write(ilog,*) 'In non-periodic, soft boundaries, this will usually indicate that&
 & the system is poorly equilibrated and/or that the simulation has become unstable.'
      end if
    end if
  end if
!
  gridnavi = grid%dim(2)*grid%dim(3)*(xi-1) + grid%dim(3)*(yi-1)+ zi
  return
!
end
!
!-----------------------------------------------------------------------
!
subroutine gridpvec(w,pvec)
!
  use mcgrid
!
  implicit none
!
  integer w,xi,yi,zi
  RTYPE pvec(3)
!
  call getgridtriple(w,xi,yi,zi)
  pvec(1) = grid%origin(1) + (xi-0.5)*grid%deltas(1)
  pvec(2) = grid%origin(2) + (yi-0.5)*grid%deltas(2)
  pvec(3) = grid%origin(3) + (zi-0.5)*grid%deltas(3)
!
end
!
!-----------------------------------------------------------------------
!
subroutine gridmindiff(w1,w2,dis2)
!
  use iounit
  use mcgrid
  use system
!
  implicit none
!
  integer w1,w2,xi1,xi2,yi1,yi2,zi1,zi2,dxi,dyi,dzi
  RTYPE dvec(3),dis2
!
  call getgridtriple(w1,xi1,yi1,zi1)
  call getgridtriple(w2,xi2,yi2,zi2)
!
  dxi = abs(xi2 - xi1)
  dyi = abs(yi2 - yi1)
  dzi = abs(zi2 - zi1)
  if (bnd_type.eq.1) then
    if (bnd_shape.eq.1) then
      if (dxi.gt.(grid%dim(1)/2)) dxi = -dxi + grid%dim(1)
      if (dyi.gt.(grid%dim(2)/2)) dyi = -dyi + grid%dim(2)
      if (dzi.gt.(grid%dim(3)/2)) dzi = -dzi + grid%dim(3)
    else if (bnd_shape.eq.3) then
      if (dzi.gt.(grid%dim(3)/2)) dzi = -dzi + grid%dim(3)
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in grid&
 &diff(...) (code # ',bnd_shape,').'
      call fexit()
    end if
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition&
 & in griddiff(...) (code # ',bnd_type,').'
    call fexit()
  end if
  dvec(1) = max(0,dxi-1)*grid%deltas(1)
  dvec(2) = max(0,dyi-1)*grid%deltas(2)
  dvec(3) = max(0,dzi-1)*grid%deltas(3)
  dis2 = sum(dvec(:)*dvec(:))
!
end
!
!------------------------------------------------------------------------
!
subroutine getgridtriple(w,xi,yi,zi)
!
  use iounit
  use mcgrid
!
  implicit none
!
  integer buf,w,xi,yi,zi
!
  buf = w
  xi = int((w-1)/(grid%dim(2)*grid%dim(3))) + 1
  buf = buf - (xi-1)*grid%dim(2)*grid%dim(3)
!
  yi = int((buf-1)/grid%dim(3)) + 1
  zi = buf - (yi-1)*grid%dim(3)
!
  buf = grid%dim(2)*grid%dim(3)*(xi-1) + grid%dim(3)*(yi-1) + zi
  if (buf.ne.w) write(*,*) 'NOOOOOOOOOOOO'
!
end
!
!------------------------------------------------------------------------
!
subroutine updateresgp(rs)
!
  use iounit
  use sequen
  use atoms
  use mcgrid
!
  implicit none
!
  RTYPE refp(3)
  integer rs,gridnavi,wo,wn,i,j
!
  wo = grid%resgp(rs)
  !write(ilog,*) 'GRP# ',grpnr, 'OLD# ',wo, 'W/ ',grid%nbnr(wo)
  refp(1) = x(refat(rs))
  refp(2) = y(refat(rs))
  refp(3) = z(refat(rs))
  wn = gridnavi(refp)
  if ((wo.le.0).OR.(wo.gt.grid%pnts)) then
    write(ilog,*) 'Fatal. Grid-point was mis-assigned. This is either an exploded simulation or a bug.'
    call fexit()
  end if
!
  if (wn.eq.wo) then !nothing to do
    return 
  else !gotta switch group from old gp to new one
!    write(*,*) 'updating ',wo,' to ',wn,' for ',rs
    do i=1,grid%nbnr(wo)
      if (rs.eq.grid%nblist(wo,i)) then
        do j=i+1,grid%nbnr(wo)
          grid%nblist(wo,j-1) = grid%nblist(wo,j)
        end do
        exit
      end if
    end do
    grid%nbnr(wo) = grid%nbnr(wo) - 1
    grid%nbnr(wn) = grid%nbnr(wn) + 1
    if (grid%nbnr(wn).gt.grid%maxnbs) then
      write(ilog,*) 'Warning. Exceeded static maximum for groups associated to a single grid point in updateresgp().&
 & Increasing the value for FMCSC_GRIDMAXRSNB dynamically which may cause memory exceptions.'
      call gridnblrsz()
    end if
    grid%nblist(wn,grid%nbnr(wn)) = rs
    grid%resgp(rs) = wn
  end if
!
end
!
!-----------------------------------------------------------------
!
subroutine gridnblrsz()
!
  use mcgrid
  use iounit
!
  implicit none
!
  integer, ALLOCATABLE:: onelist(:,:)
  integer i
!
  allocate(onelist(grid%pnts,grid%maxnbs))
!
  do i=1,grid%pnts
    onelist(i,1:grid%maxnbs) = grid%nblist(i,1:grid%maxnbs)
  end do
  deallocate(grid%nblist)
  grid%maxnbs = ceiling(1.2*grid%maxnbs)
  allocate(grid%nblist(grid%pnts+1,grid%maxnbs))
  do i=1,grid%pnts
    grid%nblist(i,1:grid%maxnbs) = onelist(i,1:grid%maxnbs)
  end do
!
  deallocate(onelist)
!
  end
!
!-----------------------------------------------------------------
!
subroutine gridreport()
!
  use atoms
  use sequen
  use mcgrid
  use iounit
  use aminos
!
  implicit none
!
  integer i,j,k,kk,tg
  character(3) resname
!
 20   format('# ',i6,' (',a3,'): Ref. at ',f7.3,1x,f7.3,1x,f7.3)
 21   format('Gridpoint ',i6,' has ',i4,' residues.')
 22   format('There are ',i8,' residues on the grid.')
  write(ilog,*)
  write(ilog,*) '--- Summary of grid occupation ---'
  write(ilog,*)
!
  tg = 0
  do i=1,grid%dim(1)*grid%dim(2)*grid%dim(3)
    if (grid%nbnr(i).gt.0) then
      write(ilog,21) i,grid%nbnr(i)
    end if
    tg = tg + grid%nbnr(i)
    do j=1,grid%nbnr(i)

      k = grid%nblist(i,j)
      resname = amino(seqtyp(k))
      kk = refat(k)
      write(ilog,20) k,resname,x(kk),y(kk),z(kk)
    end do
  end do
  write(ilog,22) tg
  write(ilog,*)
!
end
!
!-----------------------------------------------------------------------------
!
! in both grd_respairs() and respairs() we're going to use one diagonal
! half-matrix (first index larger) for the standard IPP/LJ cutoff
! and the other one (first index smaller) for the EL cutoff
! additional terms are turned on on the EL side for charged species
!
subroutine grd_respairs(rs)
!
  use iounit
  use cutoffs
  use sequen
  use mcgrid
  use energies
  use system
  use atoms
!
  implicit none
!
  integer rs,rs2,nbi,ii,i,w
  RTYPE dis,dis2
!
  w = grid%resgp(rs)
!
! add groups associated with same grid-point (incl. self-term)
  do nbi=1,grid%nbnr(w)
    rs2 = grid%nblist(w,nbi)
!   with extra check
    call dis_bound(refat(rs),refat(rs2),dis2)
    dis = sqrt(dis2)
!   now check for the short-range cutoff
    if (dis.lt.(mcnb_cutoff+resrad(rs)+resrad(rs2))) then
      if (rs2.ge.rs) then
        rsp_mat(rs2,rs) = 1
      else
        rsp_mat(rs,rs2) = 1
      end if
    end if
    if ((use_POLAR.EQV..false.).AND.(use_TABUL.EQV..false.)) cycle
!   within standard "twin"-range second cutoff (TABUL, POLAR)
    if (dis.lt.(mcel_cutoff+resrad(rs)+resrad(rs2))) then
      if (rs2.ge.rs) then
        rsp_mat(rs,rs2) = 1
      else
        rsp_mat(rs2,rs) = 1
      end if 
    end if
!   this would be the version without extra check
!    rsp_mat(rs,rs2) = 1
!    rsp_mat(rs2,rs) = 1
  end do
! scan gp neighbors
  do i=1,grid%gpnbnr(w)
    ii = grid%gpnblist(w,i)
    do nbi=1,grid%nbnr(ii)
      rs2 = grid%nblist(ii,nbi)
!     with extra explicit check
      call dis_bound(refat(rs),refat(rs2),dis2)
      dis = sqrt(dis2)
!     now check for the short-range cutoff
      if (dis.lt.(mcnb_cutoff+resrad(rs)+resrad(rs2))) then
        if (rs2.ge.rs) then
          rsp_mat(rs2,rs) = 1
        else
          rsp_mat(rs,rs2) = 1
        end if
      end if
      if ((use_POLAR.EQV..false.).AND.(use_TABUL.EQV..false.)) cycle
!     within standard "twin"-range second cutoff (TABUL, POLAR)
      if (dis.lt.(mcel_cutoff+resrad(rs)+resrad(rs2))) then
        if (rs2.ge.rs) then
          rsp_mat(rs,rs2) = 1
        else
          rsp_mat(rs2,rs) = 1
        end if 
      end if
!     this would be the version without extra explicit check
!     this additional screen is expensive but crucial if the resrads in the simulation are very different
!     such as in the case of a macromolecule in a solvent box
!      chkdis = 1.0*grid%gpnbmd(w,i)-resrad(rs)-resrad(rs2)
!      if (chkdis.gt.mcel_cutoff) cycle
!      if (rs.ge.rs2) then
!        rsp_mat(rs2,rs) = 1
!      else
!        rsp_mat(rs,rs2) = 1
!      end if
!      if (chkdis.gt.mcnb_cutoff) cycle
!      if (rs.ge.rs2) then
!        rsp_mat(rs,rs2) = 1
!      else
!        rsp_mat(rs2,rs) = 1
!      end if
    end do
  end do
  do i=1,grid%b_gpnbnr(w)
    ii = grid%b_gpnblist(w,i)
    do nbi=1,grid%nbnr(ii)
      rs2 = grid%nblist(ii,nbi)
!     with extra explicit check
      call dis_bound(refat(rs),refat(rs2),dis2)
      dis = sqrt(dis2)
!     now check for the short-range cutoff
      if (dis.lt.(mcnb_cutoff+resrad(rs)+resrad(rs2))) then
        if (rs2.ge.rs) then
          rsp_mat(rs2,rs) = 1
        else
          rsp_mat(rs,rs2) = 1
        end if
      end if
      if ((use_POLAR.EQV..false.).AND.(use_TABUL.EQV..false.)) cycle
!     within standard "twin"-range second cutoff (TABUL, POLAR)
      if (dis.lt.(mcel_cutoff+resrad(rs)+resrad(rs2))) then
        if (rs2.ge.rs) then
          rsp_mat(rs,rs2) = 1
        else
          rsp_mat(rs2,rs) = 1
        end if 
      end if
!     this would be the version without extra explicit check
!     this additional screen is expensive but crucial if the resrads in the simulation are very different
!     such as in the case of a macromolecule in a solvent box
!      chkdis = 1.0*grid%b_gpnbmd(w,i)-resrad(rs)-resrad(rs2)
!      if (chkdis.gt.mcel_cutoff) cycle
!      if (rs.ge.rs2) then
!        rsp_mat(rs2,rs) = 1
!      else
!        rsp_mat(rs,rs2) = 1
!      end if
!      if (chkdis.gt.mcnb_cutoff) cycle
!      if (rs.ge.rs2) then
!        rsp_mat(rs,rs2) = 1
!      else
!        rsp_mat(rs2,rs) = 1
!      end if
    end do
  end do
!
  if (use_POLAR.EQV..false.) return
!
  if (lrel_mc.eq.1) then
    if (chgflag(rs).EQV..true.) then
      do rs2=1,rs
        if (rsp_mat(rs2,rs).eq.0) rsp_mat(rs2,rs) = 2
      end do
      do rs2=rs+1,nseq
        if (rsp_mat(rs,rs2).eq.0) rsp_mat(rs,rs2) = 2
      end do
    else
      do i=1,cglst%ncs
        ii = cglst%it(i)
        rs2 = atmres(cglst%it(i))
        if (chgflag(rs2).EQV..false.) cycle ! may be patched
        if (rs2.gt.rs) then
          if (rsp_mat(rs,rs2).eq.0) rsp_mat(rs,rs2) = 2
        else
          if (rsp_mat(rs2,rs).eq.0) rsp_mat(rs2,rs) = 2
        end if
      end do
    end if
  else if ((lrel_mc.eq.2).OR.(lrel_mc.eq.3)) then
    if (chgflag(rs).EQV..true.) then
      do i=1,cglst%ncs
        ii = cglst%it(i)
        rs2 = atmres(cglst%it(i))
        if (chgflag(rs2).EQV..false.) cycle ! may be patched
        if (rs2.gt.rs) then
          if (rsp_mat(rs,rs2).eq.0) rsp_mat(rs,rs2) = 2
        else
          if (rsp_mat(rs2,rs).eq.0) rsp_mat(rs2,rs) = 2
        end if
      end do
    end if
! do nothing for lrel_mc = 4
  end if
!
end
!
!-----------------------------------------------------------------------------
!
  function gridwrap(v1,d1)
!
  implicit none
!
  integer v1,d1,gridwrap
!
  if (v1.gt.d1) gridwrap = v1 - d1
  if (v1.lt.1) gridwrap = v1 + d1
!
  end 
!
!-----------------------------------------------------------------------------
!
! this is the dynamics-relevant neighbor list routine
! which creates (somewhat smart) neighbor list candidate index lists, neighbor lists themselves,
! and does all the memory allocation  
!
subroutine all_respairs_nbl(cycle_frz,cycle_tab,tpi)
!
  use iounit
  use cutoffs
  use sequen
  use mcgrid
  use energies
  use system
  use tabpot
  use mcsums
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer i,j,rs1,rs2,nbi,nbj,jj,rs,tpi,sta,sto,incr
  logical cycle_frz,cycle_tab
  RTYPE chkdis
#ifdef ENABLE_THREADS
  integer tpn,OMP_GET_NUM_THREADS,modtpn
#endif
!
#ifdef ENABLE_THREADS
  tpn = omp_get_num_threads()
  sta = tpi
  sto = nseq
  incr = tpn
  tpn = 1
#else
  tpi = 1
  sta = 1
  sto = nseq
  incr = 1
#endif
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) rs_nbl(:)%ntmpanb = 0
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP BARRIER
#endif
!
  if (use_cutoffs.EQV..true.) then
!
    if (use_mcgrid.EQV..true.) then
!  
      do i=1,grid%pnts
!       add all pairs within two pre-lists
        do nbi=1,grid%nbnr(i)
          do nbj=nbi+1,grid%nbnr(i)
            rs1 = grid%nblist(i,nbi)
            rs2 = grid%nblist(i,nbj)

            if (cycle_tab.EQV..true.) then
              if (tbp%rsmat(rs1,rs2).le.0) cycle
            end if
            if (rs1.gt.rs2) then
#ifdef ENABLE_THREADS
              modtpn = mod(rs2,tpn)
              if (modtpn.ne.(tpi-1)) cycle
#endif
              rs_nbl(rs2)%ntmpanb = rs_nbl(rs2)%ntmpanb + 1
              if (rs_nbl(rs2)%ntmpanb.gt.rs_nbl(rs2)%tmpalsz) then
                call nbl_resz(rs2,7,tpi)
              end if
              rs_nbl(rs2)%tmpanb(rs_nbl(rs2)%ntmpanb) = rs1
!              if (tpi.ne.1) write(*,*) rs2,rs_nbl(rs2)%ntmpanb,rs_nbl(rs2)%tmpanb(rs_nbl(rs2)%ntmpanb)
            else
#ifdef ENABLE_THREADS
              modtpn = mod(rs1,tpn)
              if (modtpn.ne.(tpi-1)) cycle
#endif
              rs_nbl(rs1)%ntmpanb = rs_nbl(rs1)%ntmpanb + 1
              if (rs_nbl(rs1)%ntmpanb.gt.rs_nbl(rs1)%tmpalsz) then
                call nbl_resz(rs1,7,tpi)
              end if
              rs_nbl(rs1)%tmpanb(rs_nbl(rs1)%ntmpanb) = rs2
!              if (tpi.ne.1) write(*,*) rs1,rs_nbl(rs1)%ntmpanb,rs_nbl(rs1)%tmpanb(rs_nbl(rs1)%ntmpanb)
            end if
          end do
        end do
        do j=1,grid%gpnbnr(i)
          jj = grid%gpnblist(i,j)
!         add all pairs between current point and proximal points
          do nbj=1,grid%nbnr(jj)
            do nbi=1,grid%nbnr(i)
              rs1 = grid%nblist(i,nbi)
              rs2 = grid%nblist(jj,nbj)
              if (cycle_tab.EQV..true.) then
                if (tbp%rsmat(rs1,rs2).le.0) cycle
              end if
!             this additional screen is expensive but crucial if the resrads in the simulation are very different
!             such as in the case of a macromolecule in a solvent box
              chkdis = 1.0*grid%gpnbmd(i,j)-resrad(rs1)-resrad(rs2)
              if (chkdis.gt.mcel_cutoff) cycle
              if (rs1.gt.rs2) then
#ifdef ENABLE_THREADS
                modtpn = mod(rs2,tpn)
                if (modtpn.ne.(tpi-1)) cycle
#endif
                rs_nbl(rs2)%ntmpanb = rs_nbl(rs2)%ntmpanb + 1
                if (rs_nbl(rs2)%ntmpanb.gt.rs_nbl(rs2)%tmpalsz) then
                  call nbl_resz(rs2,7,tpi)
                end if
                rs_nbl(rs2)%tmpanb(rs_nbl(rs2)%ntmpanb) = rs1
              else
#ifdef ENABLE_THREADS
                modtpn = mod(rs1,tpn)
                if (modtpn.ne.(tpi-1)) cycle
#endif
                rs_nbl(rs1)%ntmpanb = rs_nbl(rs1)%ntmpanb + 1
                if (rs_nbl(rs1)%ntmpanb.gt.rs_nbl(rs1)%tmpalsz) then
                  call nbl_resz(rs1,7,tpi)
                end if
                rs_nbl(rs1)%tmpanb(rs_nbl(rs1)%ntmpanb) = rs2
              end if
            end do
          end do
        end do
      end do
!
    else if (use_rescrit.EQV..true.) then
!
!     do nothing
!
    else
!
      write(ilog,*) 'Fatal. Encountered unsupported cutoff treatment in all_respairs_nbl(...).&
 & Please report this bug.'
      call fexit()
    end if
!
  else
!
    write(ilog,*) 'Fatal. Called all_respairs_nbl(...) with cutoffs turned off.'
    call fexit()
!
  end if
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
  do rs=sta,sto,incr ! 1,nseq
    call respairs_nbl(rs,cycle_frz,cycle_tab,tpi)
  end do
!
end
!
!-----------------------------------------------------------------------------
!
! this is the MC-relevant "neighbor-list" routine, which simply populates
! the residue-by-residue matrix rsp_mat(i,j) for i > j SR terms, and for j > i
! LR terms
!
subroutine respairs(rs)
!
  use iounit
  use sequen
  use energies
  use cutoffs
#if ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer i,ii,kk,rs,sta,sto
  RTYPE dis,dis2
#ifdef ENABLE_THREADS
  integer tpi,tpn,OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif
!
  kk = refat(rs)
!
#if ENABLE_THREADS
  call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpn,tpi,sta,sto,i,ii,dis,dis2) FIRSTPRIVATE(kk)
  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  sto = tpi*nseq/tpn
  sta = (tpi-1)*nseq/tpn + 1
#else
  sta=1
  sto=nseq
#endif
  
  do i=sta,sto
    ii = refat(i)
    call dis_bound(ii,kk,dis2)
    dis = sqrt(dis2)
!   now check for the short-range cutoff
    if (dis.lt.(mcnb_cutoff+resrad(rs)+resrad(i))) then

      if (i.ge.rs) then
        rsp_mat(i,rs) = 1
      else
        rsp_mat(rs,i) = 1
      end if
    end if

    if ((use_POLAR.EQV..false.).AND.(use_TABUL.EQV..false.)) cycle
!   within standard "twin"-range second cutoff (TABUL, POLAR)
    if (dis.lt.(mcel_cutoff+resrad(rs)+resrad(i))) then
      if (i.ge.rs) then
        rsp_mat(rs,i) = 1
      else
        rsp_mat(i,rs) = 1
      end if 
    else if (use_POLAR.EQV..true.) then
!     depending on LR electrostatics, turn on extra terms beyond that
      if (lrel_mc.eq.1) then
        if ((chgflag(rs).EQV..true.).OR.(chgflag(i).EQV..true.)) then
          if (i.ge.rs) then
            rsp_mat(rs,i) = 2
          else
            rsp_mat(i,rs) = 2
          end if
        end if
      else if ((lrel_mc.eq.2).OR.(lrel_mc.eq.3)) then
        if ((chgflag(rs).EQV..true.).AND.(chgflag(i).EQV..true.)) then
          if (i.ge.rs) then
            rsp_mat(rs,i) = 2
          else
            rsp_mat(i,rs) = 2
          end if
        end if
!     do nothing for lrel_mc = 4
      end if
    end if
  end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP END PARALLEL
#endif
!
end
!
!
!-----------------------------------------------------------------------------
!
! this extremely simple routine creates a later on static rs-rs neighbor list based
! on the presence of tabulated interactions (which may be sparse)
! note that this does not respect constraints (i.e, might be inefficient still)
!
  subroutine tab_respairs_nbl_pre(rs)
!
  use iounit
  use molecule
  use sequen
  use cutoffs
  use atoms
  use tabpot
!
  implicit none
!
  integer i,ii,rs
  integer, ALLOCATABLE:: inb(:)
!
  allocate(inb(nseq-rs))
!
  rs_nbl(rs)%ntabnbs = 0
  rs_nbl(rs)%ntabias = 0
!
  ii = refat(rs)
  do i=rs+1,nseq
    if (tbp%rsmat(rs,i).gt.0) then
      rs_nbl(rs)%ntabnbs = rs_nbl(rs)%ntabnbs + 1
      rs_nbl(rs)%ntabias = rs_nbl(rs)%ntabias + (tbp%rsmat(i,rs)-tbp%rsmat(rs,i)) + 1
      inb(rs_nbl(rs)%ntabnbs) = i
    end if
  end do
!  write(*,*) rs,rs_nbl(rs)%ntabnbs,rs_nbl(rs)%ntabias
!
  if (allocated(rs_nbl(rs)%tabnb).EQV..true.) then
    deallocate(rs_nbl(rs)%tabnb)
  end if
  if (rs_nbl(rs)%ntabnbs.gt.0) then
    allocate(rs_nbl(rs)%tabnb(rs_nbl(rs)%ntabnbs))
    rs_nbl(rs)%tabnb(1:rs_nbl(rs)%ntabnbs) = inb(1:rs_nbl(rs)%ntabnbs)
  end if
!
  deallocate(inb)
!
  end
!
!-----------------------------------------------------------------------------
!
! mode: 1 use all rsi > rs (assumes tmpanb not in use)
!       2 use those in rs_nbl(rs)%tmpanb
! cycle_frz: true exclude frozen ia.s from either
!            false do not exclude
! cycle_tab: true exclude non-tab-interacting ia.s from either
!            false do not exclude
!
subroutine setup_nblidx(rs,cycle_frz,cycle_tab,mode,tpi)
!
  use iounit
  use molecule
  use sequen
  use energies
  use cutoffs
  use polypep
  use system
  use atoms
  use tabpot
  use fyoc
!
  implicit none
!
  integer i,ii,rs,j,hi,mode,imol,jmol,tmpnb,tpi
  logical cycle_frz,cycle_tab
!
! initialize
  rs_nbl(rs)%nnbs = 0
  rs_nbl(rs)%nnbtrs = 0
  rs_nbl(rs)%nnblrs = 0
  rs_nbl(rs)%nnbats = 0
  rs_nbl(rs)%nnbtrats = 0
  rs_nbl(rs)%nnblrats = 0
  if (use_FEG.EQV..true.) then
    rs_nbl(rs)%ngnbs = 0
    rs_nbl(rs)%ngnbtrs = 0
    rs_nbl(rs)%ngnbats = 0
    rs_nbl(rs)%ngnbtrats = 0
    rs_nbl(rs)%ngnblrs = 0
    rs_nbl(rs)%ngnblrats = 0
  end if
!
  if (rs.lt.nseq) then
    rsp_vec(rs) = 0
  end if
!
  if (rs.eq.nseq) return
  if ((mode.eq.2).AND.(rs_nbl(rs)%ntmpanb.eq.0)) return
!
  imol = molofrs(rs)
  ii = refat(rs)
!
  if (mode.eq.1) then
!
    tmpnb = 0
!
    if (cycle_frz.EQV..true.) then
      if (cycle_tab.EQV..true.) then
        do j=1,rs_nbl(rs)%ntabnbs
          i = rs_nbl(rs)%tabnb(j)
          jmol = molofrs(i)
          if ((jmol.eq.imol).AND.(molfrzidx(imol).lt.3)) then
            if (molfrzidx(imol).eq.2) then
              if ((nchi(rs).eq.0).AND.(nchi(i).eq.0)) then
!               do nothing (backbone is frozen completely)
              else
                tmpnb = tmpnb + 1
                tmp_nbl(tpi)%dvec(tmpnb,1) = x(refat(i)) - x(ii)
                tmp_nbl(tpi)%dvec(tmpnb,2) = y(refat(i)) - y(ii)
                tmp_nbl(tpi)%dvec(tmpnb,3) = z(refat(i)) - z(ii)
                tmp_nbl(tpi)%idx(tmpnb) = i
              end if
            else
              tmpnb = tmpnb + 1
              tmp_nbl(tpi)%dvec(tmpnb,1) = x(refat(i)) - x(ii)
              tmp_nbl(tpi)%dvec(tmpnb,2) = y(refat(i)) - y(ii)
              tmp_nbl(tpi)%dvec(tmpnb,3) = z(refat(i)) - z(ii)
              tmp_nbl(tpi)%idx(tmpnb) = i
            end if
          else
            if ((molfrzidx(jmol).eq.4).AND.(molfrzidx(imol).eq.4)) then
!             do nothing (note there's no other way for the two molecules' interactions to never change)
            else
              tmpnb = tmpnb + 1
              tmp_nbl(tpi)%dvec(tmpnb,1) = x(refat(i)) - x(ii)
              tmp_nbl(tpi)%dvec(tmpnb,2) = y(refat(i)) - y(ii)
              tmp_nbl(tpi)%dvec(tmpnb,3) = z(refat(i)) - z(ii)
              tmp_nbl(tpi)%idx(tmpnb) = i
            end if
          end if
        end do
      else
        if (molfrzidx(imol).lt.3) then
          do i=rs+1,rsmol(imol,2)
            if (cycle_tab.EQV..true.) then
              if (tbp%rsmat(rs,i).le.0) cycle
            end if
            if (molfrzidx(imol).eq.2) then
              if ((nchi(rs).eq.0).AND.(nchi(i).eq.0)) then
!               do nothing (backbone is frozen completely)
              else
                tmpnb = tmpnb + 1
                tmp_nbl(tpi)%dvec(tmpnb,1) = x(refat(i)) - x(ii)
                tmp_nbl(tpi)%dvec(tmpnb,2) = y(refat(i)) - y(ii)
                tmp_nbl(tpi)%dvec(tmpnb,3) = z(refat(i)) - z(ii)
                tmp_nbl(tpi)%idx(tmpnb) = i
              end if
            else
              tmpnb = tmpnb + 1
              tmp_nbl(tpi)%dvec(tmpnb,1) = x(refat(i)) - x(ii)
              tmp_nbl(tpi)%dvec(tmpnb,2) = y(refat(i)) - y(ii)
              tmp_nbl(tpi)%dvec(tmpnb,3) = z(refat(i)) - z(ii)
              tmp_nbl(tpi)%idx(tmpnb) = i
            end if
          end do
        end if
        do jmol=imol+1,nmol
          if ((molfrzidx(jmol).eq.4).AND.(molfrzidx(imol).eq.4)) then
!           do nothing (note there's no other way for the two molecules' interactions to never change)
          else
            do i=rsmol(jmol,1),rsmol(jmol,2)
              if (cycle_tab.EQV..true.) then
                if (tbp%rsmat(rs,i).le.0) cycle
              end if
              tmpnb = tmpnb + 1
              tmp_nbl(tpi)%dvec(tmpnb,1) = x(refat(i)) - x(ii)
              tmp_nbl(tpi)%dvec(tmpnb,2) = y(refat(i)) - y(ii)
              tmp_nbl(tpi)%dvec(tmpnb,3) = z(refat(i)) - z(ii)
              tmp_nbl(tpi)%idx(tmpnb) = i
            end do
          end if
        end do
      end if
!
      hi = tmpnb
!
    else if (cycle_tab.EQV..true.) then ! but cycle_frz false
!
      tmpnb = 0
      do j=1,rs_nbl(rs)%ntabnbs
        i = rs_nbl(rs)%tabnb(j)
        if (tbp%rsmat(rs,i).le.0) cycle ! should never cycle
        tmpnb = tmpnb + 1
        tmp_nbl(tpi)%dvec(tmpnb,1) = x(refat(i)) - x(ii)
        tmp_nbl(tpi)%dvec(tmpnb,2) = y(refat(i)) - y(ii)
        tmp_nbl(tpi)%dvec(tmpnb,3) = z(refat(i)) - z(ii)
        tmp_nbl(tpi)%idx(tmpnb) = i
      end do
!
      hi = tmpnb
!
    else ! regular case
!
      hi = nseq-rs
      tmp_nbl(tpi)%dvec(1:hi,1) = x(refat(rs+1:nseq)) - x(ii)
      tmp_nbl(tpi)%dvec(1:hi,2) = y(refat(rs+1:nseq)) - y(ii)
      tmp_nbl(tpi)%dvec(1:hi,3) = z(refat(rs+1:nseq)) - z(ii)
      do i=rs+1,nseq
        tmp_nbl(tpi)%idx(i-rs) = i
      end do
!
    end if
!
  else if (mode.eq.2) then
!
    tmpnb = 0
!
    if (cycle_frz.EQV..true.) then
      do j=1,rs_nbl(rs)%ntmpanb
        i = rs_nbl(rs)%tmpanb(j)
        jmol = molofrs(i)
        if (cycle_tab.EQV..true.) then
          if (tbp%rsmat(rs,i).le.0) cycle
        end if
        if (imol.eq.jmol) then
          if (molfrzidx(imol).lt.3) then
            if (molfrzidx(imol).eq.2) then
!             the better approach here would be to scan only the segment between the two residues
!             for rigidity instead of relying in the unnecessarily stringent criterion in molfrzidx
              if ((nchi(rs).eq.0).AND.(nchi(i).eq.0)) then
!               do nothing (backbone is frozen completely)
              else
                tmpnb = tmpnb + 1
                tmp_nbl(tpi)%dvec(tmpnb,1) = x(refat(i)) - x(ii)
                tmp_nbl(tpi)%dvec(tmpnb,2) = y(refat(i)) - y(ii)
                tmp_nbl(tpi)%dvec(tmpnb,3) = z(refat(i)) - z(ii)
                tmp_nbl(tpi)%idx(tmpnb) = i
              end if
            else
              tmpnb = tmpnb + 1
              tmp_nbl(tpi)%dvec(tmpnb,1) = x(refat(i)) - x(ii)
              tmp_nbl(tpi)%dvec(tmpnb,2) = y(refat(i)) - y(ii)
              tmp_nbl(tpi)%dvec(tmpnb,3) = z(refat(i)) - z(ii)
              tmp_nbl(tpi)%idx(tmpnb) = i
            end if
          end if
        else      
          if ((molfrzidx(jmol).eq.4).AND.(molfrzidx(imol).eq.4)) then
!           do nothing (note there's no other way for the two molecules' interactions to never change)
          else
            tmpnb = tmpnb + 1
            tmp_nbl(tpi)%dvec(tmpnb,1) = x(refat(i)) - x(ii)
            tmp_nbl(tpi)%dvec(tmpnb,2) = y(refat(i)) - y(ii)
            tmp_nbl(tpi)%dvec(tmpnb,3) = z(refat(i)) - z(ii)
            tmp_nbl(tpi)%idx(tmpnb) = i
          end if
        end if
      end do
!
      hi = tmpnb
!
    else if (cycle_tab.EQV..true.) then ! but cycle_frz false (this may be inefficient)
!
      tmpnb = 0
!     note that the cutoff-aware pre-nb-list is also aware of tbp%rsmat
      do i=1,rs_nbl(rs)%ntmpanb
        if (tbp%rsmat(rs,rs_nbl(rs)%tmpanb(i)).le.0) cycle ! should never cycle
        tmpnb = tmpnb + 1
        tmp_nbl(tpi)%dvec(tmpnb,1) = x(refat(rs_nbl(rs)%tmpanb(i))) - x(ii)
        tmp_nbl(tpi)%dvec(tmpnb,2) = y(refat(rs_nbl(rs)%tmpanb(i))) - y(ii)
        tmp_nbl(tpi)%dvec(tmpnb,3) = z(refat(rs_nbl(rs)%tmpanb(i))) - z(ii)
        tmp_nbl(tpi)%idx(tmpnb) = rs_nbl(rs)%tmpanb(i)
      end do
      hi = tmpnb
!
    else ! regular case
!
      hi = rs_nbl(rs)%ntmpanb
      do i=1,hi
        tmp_nbl(tpi)%dvec(i,1) = x(refat(rs_nbl(rs)%tmpanb(i))) - x(ii)
        tmp_nbl(tpi)%dvec(i,2) = y(refat(rs_nbl(rs)%tmpanb(i))) - y(ii)
        tmp_nbl(tpi)%dvec(i,3) = z(refat(rs_nbl(rs)%tmpanb(i))) - z(ii)
        tmp_nbl(tpi)%idx(i) = rs_nbl(rs)%tmpanb(i)
      end do
!
    end if
!
  end if
!
! now take care of the shift vectors
!
! PBC
  if ((bnd_type.eq.1).AND.(hi.ge.1)) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        tmp_nbl(tpi)%k1(1:hi,j) = floor((-0.5*bnd_params(j)-tmp_nbl(tpi)%dvec(1:hi,j))/bnd_params(j))+1
        tmp_nbl(tpi)%k1(1:hi,j) = max(tmp_nbl(tpi)%k1(1:hi,j),0)
        tmp_nbl(tpi)%k2(1:hi,j) =  floor((tmp_nbl(tpi)%dvec(1:hi,j)-0.5*bnd_params(j))/bnd_params(j))+1
        tmp_nbl(tpi)%k2(1:hi,j) = max(tmp_nbl(tpi)%k2(1:hi,j),0)
        tmp_nbl(tpi)%svec(j,1:hi) = (tmp_nbl(tpi)%k1(1:hi,j) - tmp_nbl(tpi)%k2(1:hi,j))*bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      tmp_nbl(tpi)%svec(1,1:hi) = 0.0
      tmp_nbl(tpi)%svec(2,1:hi) = 0.0
      tmp_nbl(tpi)%k1(1:hi,3) = floor((-0.5*bnd_params(6)-tmp_nbl(tpi)%dvec(1:hi,3))/bnd_params(6))+1
      tmp_nbl(tpi)%k1(1:hi,3) = max(tmp_nbl(tpi)%k1(1:hi,3),0)
      tmp_nbl(tpi)%k2(1:hi,3) =  floor((tmp_nbl(tpi)%dvec(1:hi,3)-0.5*bnd_params(6))/bnd_params(6))+1
      tmp_nbl(tpi)%k2(1:hi,3) = max(tmp_nbl(tpi)%k2(1:hi,3),0)
      tmp_nbl(tpi)%svec(3,1:hi) = (tmp_nbl(tpi)%k1(1:hi,3) - tmp_nbl(tpi)%k2(1:hi,3))*bnd_params(6)
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in res&
 &pairs_nbl() (code # ',bnd_shape,').'
      call fexit()
    end if
    do i=1,hi
      if (molofrs(tmp_nbl(tpi)%idx(i)).eq.molofrs(rs)) then
        tmp_nbl(tpi)%svec(:,i) = 0.0
      else
        tmp_nbl(tpi)%dvec(i,1) = tmp_nbl(tpi)%dvec(i,1) + tmp_nbl(tpi)%svec(1,i)
        tmp_nbl(tpi)%dvec(i,2) = tmp_nbl(tpi)%dvec(i,2) + tmp_nbl(tpi)%svec(2,i)
        tmp_nbl(tpi)%dvec(i,3) = tmp_nbl(tpi)%dvec(i,3) + tmp_nbl(tpi)%svec(3,i)
      end if
    end do
  else
    tmp_nbl(tpi)%svec(:,1:hi) = 0.0
  end if
!
  tmp_nbl(tpi)%d1(1:hi) = sqrt(tmp_nbl(tpi)%dvec(1:hi,1)*tmp_nbl(tpi)%dvec(1:hi,1) + &
  &                            tmp_nbl(tpi)%dvec(1:hi,2)*tmp_nbl(tpi)%dvec(1:hi,2) + &
  &                            tmp_nbl(tpi)%dvec(1:hi,3)*tmp_nbl(tpi)%dvec(1:hi,3))
!
! now call the fxn which parses those arrays
  call populate_nbl(rs,hi,tmp_nbl(tpi)%svec(:,1:hi),tmp_nbl(tpi)%d1(1:hi),tmp_nbl(tpi)%idx(1:hi),tpi)
!
  end
!
!-----------------------------------------------------------------------------
!
  subroutine populate_nbl(rs,hi,svec,d1,idx,tpi)
!
  use iounit
  use molecule
  use sequen
  use energies
  use cutoffs
  use polypep
  use system
  use atoms
  use fyoc
!
  implicit none
!
  integer hi,rs,i,ii,tpi
  RTYPE svec(3,hi),d1(hi)
  integer idx(hi)
!
  if (use_FEG.EQV..false.) then
    do ii=1,hi
      i = idx(ii)
      if (i.eq.disulf(rs)) cycle ! exclude crosslinked pairs
      if (d1(ii).lt.(mcnb_cutoff+resrad(rs)+resrad(i))) then
        if (i.eq.(rs+1)) then
          rsp_vec(rs) = 1
        else
          rs_nbl(rs)%nnbs = rs_nbl(rs)%nnbs + 1
          if (rs_nbl(rs)%nnbs.gt.rs_nbl(rs)%nbalsz) call nbl_resz(rs,1,tpi)
          rs_nbl(rs)%nnbats = rs_nbl(rs)%nnbats + at(i)%nbb + at(i)%nsc
          rs_nbl(rs)%nb(rs_nbl(rs)%nnbs) = i
          if (molofrs(i).eq.molofrs(rs)) then
            rs_nbl(rs)%svec(rs_nbl(rs)%nnbs,:) = 0.0
          else
            rs_nbl(rs)%svec(rs_nbl(rs)%nnbs,:) = svec(:,ii)
          end if
        end if
      else if (d1(ii).lt.(mcel_cutoff+resrad(rs)+resrad(i))) then
        if (i.eq.(rs+1)) then
          rsp_vec(rs) = 2
        else
          rs_nbl(rs)%nnbtrs = rs_nbl(rs)%nnbtrs + 1
          if (rs_nbl(rs)%nnbtrs.gt.rs_nbl(rs)%tralsz) call nbl_resz(rs,2,tpi)
          rs_nbl(rs)%nnbtrats = rs_nbl(rs)%nnbtrats + at(i)%nbb + at(i)%nsc
          rs_nbl(rs)%nbtr(rs_nbl(rs)%nnbtrs) = i
          if (molofrs(i).eq.molofrs(rs)) then
            rs_nbl(rs)%trsvec(rs_nbl(rs)%nnbtrs,:) = 0.0
          else
            rs_nbl(rs)%trsvec(rs_nbl(rs)%nnbtrs,:) = svec(:,ii)
          end if
        end if
      else
        if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.5).OR.(lrel_md.eq.4))) then
          if ((chgflag(rs).EQV..true.).OR.(chgflag(i).EQV..true.)) then
            if (i.eq.(rs+1)) then
!             do nothing (detected as rsp_vec(i) being 0)
            else
              if ((lrel_md.eq.4).AND.((chgflag(rs).EQV..false.).OR.(chgflag(i).EQV..false.))) cycle
              rs_nbl(rs)%nnblrs = rs_nbl(rs)%nnblrs + 1
              if (rs_nbl(rs)%nnblrs.gt.rs_nbl(rs)%lralsz) call nbl_resz(rs,3,tpi)
              rs_nbl(rs)%nnblrats = rs_nbl(rs)%nnblrats + at(i)%nbb + at(i)%nsc
              rs_nbl(rs)%nblr(rs_nbl(rs)%nnblrs) = i
              if (molofrs(i).eq.molofrs(rs)) then
                rs_nbl(rs)%lrsvec(rs_nbl(rs)%nnblrs,:) = 0.0
              else
                rs_nbl(rs)%lrsvec(rs_nbl(rs)%nnblrs,:) = svec(:,ii)
              end if
            end if
          end if
        end if
      end if
    end do
  else
!   for calculations with ghosts we have (annoyingly) separate nb-lists
    do ii=1,hi
      i = idx(ii)
      if (i.eq.disulf(rs)) cycle ! exclude crosslinked pairs
      if (d1(ii).lt.(mcnb_cutoff+resrad(rs)+resrad(i))) then
        if (i.eq.(rs+1)) then
          rsp_vec(rs) = 1
        else
          if ((par_FEG3(rs).EQV..true.).OR.(par_FEG3(i).EQV..true.)) then
!           do nothing (at least one of the residues is fully de-coupled)
          else if (((par_FEG(rs).EQV..false.).AND.(par_FEG(i).EQV..false.)).OR.&
 &    ((par_FEG(rs).EQV..true.).AND.(par_FEG(i).EQV..true.).AND.(fegmode.eq.1))) then
            rs_nbl(rs)%nnbs = rs_nbl(rs)%nnbs + 1
            if (rs_nbl(rs)%nnbs.gt.rs_nbl(rs)%nbalsz) call nbl_resz(rs,1,tpi)
            rs_nbl(rs)%nnbats = rs_nbl(rs)%nnbats + at(i)%nbb + at(i)%nsc
            rs_nbl(rs)%nb(rs_nbl(rs)%nnbs) = i
            if (molofrs(i).eq.molofrs(rs)) then
              rs_nbl(rs)%svec(rs_nbl(rs)%nnbs,:) = 0.0
            else
              rs_nbl(rs)%svec(rs_nbl(rs)%nnbs,:) = svec(:,ii)
            end if
          else
            rs_nbl(rs)%ngnbs = rs_nbl(rs)%ngnbs + 1
            if (rs_nbl(rs)%ngnbs.gt.rs_nbl(rs)%gnbalsz) call nbl_resz(rs,4,tpi)
            rs_nbl(rs)%ngnbats = rs_nbl(rs)%ngnbats + at(i)%nbb + at(i)%nsc
            rs_nbl(rs)%gnb(rs_nbl(rs)%ngnbs) = i
            if (molofrs(i).eq.molofrs(rs)) then
              rs_nbl(rs)%gsvec(rs_nbl(rs)%ngnbs,:) = 0.0
            else
              rs_nbl(rs)%gsvec(rs_nbl(rs)%ngnbs,:) = svec(:,ii)
            end if
          end if
        end if
      else if (d1(ii).lt.(mcel_cutoff+resrad(rs)+resrad(i))) then
        if (i.eq.(rs+1)) then
          rsp_vec(rs) = 2
        else
          if ((par_FEG3(rs).EQV..true.).OR.(par_FEG3(i).EQV..true.)) then
!           do nothing (at least one of the residues is fully de-coupled)
          else if (((par_FEG(rs).EQV..false.).AND.(par_FEG(i).EQV..false.)).OR.&
 &    ((par_FEG(rs).EQV..true.).AND.(par_FEG(i).EQV..true.).AND.(fegmode.eq.1))) then
            rs_nbl(rs)%nnbtrs = rs_nbl(rs)%nnbtrs + 1
            if (rs_nbl(rs)%nnbtrs.gt.rs_nbl(rs)%tralsz) call nbl_resz(rs,2,tpi)
            rs_nbl(rs)%nnbtrats = rs_nbl(rs)%nnbtrats + at(i)%nbb + at(i)%nsc
            rs_nbl(rs)%nbtr(rs_nbl(rs)%nnbtrs) = i
            if (molofrs(i).eq.molofrs(rs)) then
              rs_nbl(rs)%trsvec(rs_nbl(rs)%nnbtrs,:) = 0.0
            else
              rs_nbl(rs)%trsvec(rs_nbl(rs)%nnbtrs,:) = svec(:,ii)
            end if
          else
            rs_nbl(rs)%ngnbtrs = rs_nbl(rs)%ngnbtrs + 1
            if (rs_nbl(rs)%ngnbtrs.gt.rs_nbl(rs)%gtralsz) call nbl_resz(rs,5,tpi)
            rs_nbl(rs)%ngnbtrats = rs_nbl(rs)%ngnbtrats + at(i)%nbb + at(i)%nsc
            rs_nbl(rs)%gnbtr(rs_nbl(rs)%ngnbtrs) = i
            if (molofrs(i).eq.molofrs(rs)) then
              rs_nbl(rs)%gtrsvec(rs_nbl(rs)%ngnbtrs,:) = 0.0
            else
              rs_nbl(rs)%gtrsvec(rs_nbl(rs)%ngnbtrs,:) = svec(:,ii)
            end if
          end if
        end if
      else
        if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.5).OR.(lrel_md.eq.4))) then
          if ((chgflag(rs).EQV..true.).OR.(chgflag(i).EQV..true.)) then
            if (i.eq.(rs+1)) then
!             do nothing (detected as rsp_vec(i) being 0)
            else
              if ((lrel_md.eq.4).AND.((chgflag(rs).EQV..false.).OR.(chgflag(i).EQV..false.))) cycle
              if ((par_FEG3(rs).EQV..true.).OR.(par_FEG3(i).EQV..true.)) then
!               do nothing (at least one of the residues is fully de-coupled)
              else if (((par_FEG(rs).EQV..false.).AND.(par_FEG(i).EQV..false.)).OR.&
 &    ((par_FEG(rs).EQV..true.).AND.(par_FEG(i).EQV..true.).AND.(fegmode.eq.1))) then
                rs_nbl(rs)%nnblrs = rs_nbl(rs)%nnblrs + 1
                if (rs_nbl(rs)%nnblrs.gt.rs_nbl(rs)%lralsz) call nbl_resz(rs,3,tpi)
                rs_nbl(rs)%nnblrats = rs_nbl(rs)%nnblrats + at(i)%nbb + at(i)%nsc
                rs_nbl(rs)%nblr(rs_nbl(rs)%nnblrs) = i
                if (molofrs(i).eq.molofrs(rs)) then
                  rs_nbl(rs)%lrsvec(rs_nbl(rs)%nnblrs,:) = 0.0
                else
                  rs_nbl(rs)%lrsvec(rs_nbl(rs)%nnblrs,:) = svec(:,ii)
                end if
              else
                rs_nbl(rs)%ngnblrs = rs_nbl(rs)%ngnblrs + 1
                if (rs_nbl(rs)%ngnblrs.gt.rs_nbl(rs)%glralsz) call nbl_resz(rs,6,tpi)
                rs_nbl(rs)%ngnblrats = rs_nbl(rs)%ngnblrats + at(i)%nbb + at(i)%nsc
                rs_nbl(rs)%gnblr(rs_nbl(rs)%ngnblrs) = i
                if (molofrs(i).eq.molofrs(rs)) then
                  rs_nbl(rs)%glrsvec(rs_nbl(rs)%ngnblrs,:) = 0.0
                else
                  rs_nbl(rs)%glrsvec(rs_nbl(rs)%ngnblrs,:) = svec(:,ii)
                end if
              end if
            end if
          end if
        end if
      end if
    end do
  end if
!
  end
!
!-----------------------------------------------------------------------------
!
! a subroutine to resize a neighbor list by a fixed increment
!
  subroutine nbl_resz(rs,which,tpi)
!
  use iounit
  use cutoffs
  use sequen
!
  implicit none
!
  integer rs,which,tpi
!
! normal close-range
  if (which.eq.1) then
!   backup, deallocate, increment, reallocate, copy back
    tmp_nbl(tpi)%nbl(1:rs_nbl(rs)%nbalsz) = rs_nbl(rs)%nb(1:rs_nbl(rs)%nbalsz)
    tmp_nbl(tpi)%svl(1:rs_nbl(rs)%nbalsz,:) = rs_nbl(rs)%svec(1:rs_nbl(rs)%nbalsz,:)
    deallocate(rs_nbl(rs)%nb)
    deallocate(rs_nbl(rs)%svec)
    rs_nbl(rs)%nbalsz = min(nseq,rs_nbl(rs)%nbalsz + 10)
    allocate(rs_nbl(rs)%nb(rs_nbl(rs)%nbalsz))
    allocate(rs_nbl(rs)%svec(rs_nbl(rs)%nbalsz,3))
    rs_nbl(rs)%nb(1:rs_nbl(rs)%nbalsz) = tmp_nbl(tpi)%nbl(1:rs_nbl(rs)%nbalsz)
    rs_nbl(rs)%svec(1:rs_nbl(rs)%nbalsz,:) = tmp_nbl(tpi)%svl(1:rs_nbl(rs)%nbalsz,:)
! normal twin-range
  else if (which.eq.2) then
    tmp_nbl(tpi)%nbl(1:rs_nbl(rs)%tralsz) = rs_nbl(rs)%nbtr(1:rs_nbl(rs)%tralsz)
    tmp_nbl(tpi)%svl(1:rs_nbl(rs)%tralsz,:) = rs_nbl(rs)%trsvec(1:rs_nbl(rs)%tralsz,:)
    deallocate(rs_nbl(rs)%nbtr)
    deallocate(rs_nbl(rs)%trsvec)
    rs_nbl(rs)%tralsz = min(nseq,rs_nbl(rs)%tralsz + 10)
    allocate(rs_nbl(rs)%nbtr(rs_nbl(rs)%tralsz))
    allocate(rs_nbl(rs)%trsvec(rs_nbl(rs)%tralsz,3))
    rs_nbl(rs)%nbtr(1:rs_nbl(rs)%tralsz) = tmp_nbl(tpi)%nbl(1:rs_nbl(rs)%tralsz)
    rs_nbl(rs)%trsvec(1:rs_nbl(rs)%tralsz,:) = tmp_nbl(tpi)%svl(1:rs_nbl(rs)%tralsz,:)
! normal long-range
  else if (which.eq.3) then
    tmp_nbl(tpi)%nbl(1:rs_nbl(rs)%lralsz) = rs_nbl(rs)%nblr(1:rs_nbl(rs)%lralsz)
    tmp_nbl(tpi)%svl(1:rs_nbl(rs)%lralsz,:) = rs_nbl(rs)%lrsvec(1:rs_nbl(rs)%lralsz,:)
    deallocate(rs_nbl(rs)%nblr)
    deallocate(rs_nbl(rs)%lrsvec)
    rs_nbl(rs)%lralsz = min(nseq,rs_nbl(rs)%lralsz + 10)
    allocate(rs_nbl(rs)%nblr(rs_nbl(rs)%lralsz))
    allocate(rs_nbl(rs)%lrsvec(rs_nbl(rs)%lralsz,3))
    rs_nbl(rs)%nblr(1:rs_nbl(rs)%lralsz) = tmp_nbl(tpi)%nbl(1:rs_nbl(rs)%lralsz)
    rs_nbl(rs)%lrsvec(1:rs_nbl(rs)%lralsz,:) = tmp_nbl(tpi)%svl(1:rs_nbl(rs)%lralsz,:)
! ghosted close-range
  else if (which.eq.4) then
    tmp_nbl(tpi)%nbl(1:rs_nbl(rs)%gnbalsz) = rs_nbl(rs)%gnb(1:rs_nbl(rs)%gnbalsz)
    tmp_nbl(tpi)%svl(1:rs_nbl(rs)%gnbalsz,:) = rs_nbl(rs)%gsvec(1:rs_nbl(rs)%gnbalsz,:)
    deallocate(rs_nbl(rs)%gnb)
    deallocate(rs_nbl(rs)%gsvec)
    rs_nbl(rs)%gnbalsz = min(nseq,rs_nbl(rs)%gnbalsz + 10)
    allocate(rs_nbl(rs)%gnb(rs_nbl(rs)%gnbalsz))
    allocate(rs_nbl(rs)%gsvec(rs_nbl(rs)%gnbalsz,3))
    rs_nbl(rs)%gnb(1:rs_nbl(rs)%gnbalsz) = tmp_nbl(tpi)%nbl(1:rs_nbl(rs)%gnbalsz)
    rs_nbl(rs)%gsvec(1:rs_nbl(rs)%gnbalsz,:) = tmp_nbl(tpi)%svl(1:rs_nbl(rs)%gnbalsz,:)
! ghosted twin-range
  else if (which.eq.5) then
    tmp_nbl(tpi)%nbl(1:rs_nbl(rs)%gtralsz) = rs_nbl(rs)%gnbtr(1:rs_nbl(rs)%gtralsz)
    tmp_nbl(tpi)%svl(1:rs_nbl(rs)%gtralsz,:) = rs_nbl(rs)%gtrsvec(1:rs_nbl(rs)%gtralsz,:)
    deallocate(rs_nbl(rs)%gnbtr)
    deallocate(rs_nbl(rs)%gtrsvec)
    rs_nbl(rs)%gtralsz = min(nseq,rs_nbl(rs)%gtralsz + 10)
    allocate(rs_nbl(rs)%gnbtr(rs_nbl(rs)%gtralsz))
    allocate(rs_nbl(rs)%gtrsvec(rs_nbl(rs)%gtralsz,3))
    rs_nbl(rs)%gnbtr(1:rs_nbl(rs)%gtralsz) = tmp_nbl(tpi)%nbl(1:rs_nbl(rs)%gtralsz)
    rs_nbl(rs)%gtrsvec(1:rs_nbl(rs)%gtralsz,:) = tmp_nbl(tpi)%svl(1:rs_nbl(rs)%gtralsz,:)
! ghosted long-range
  else if (which.eq.6) then
    tmp_nbl(tpi)%nbl(1:rs_nbl(rs)%glralsz) = rs_nbl(rs)%gnblr(1:rs_nbl(rs)%glralsz)
    tmp_nbl(tpi)%svl(1:rs_nbl(rs)%glralsz,:) = rs_nbl(rs)%glrsvec(1:rs_nbl(rs)%glralsz,:)
    deallocate(rs_nbl(rs)%gnblr)
    deallocate(rs_nbl(rs)%glrsvec)
    rs_nbl(rs)%glralsz = min(nseq,rs_nbl(rs)%glralsz + 10)
    allocate(rs_nbl(rs)%gnblr(rs_nbl(rs)%glralsz))
    allocate(rs_nbl(rs)%glrsvec(rs_nbl(rs)%glralsz,3))
    rs_nbl(rs)%gnblr(1:rs_nbl(rs)%glralsz) = tmp_nbl(tpi)%nbl(1:rs_nbl(rs)%glralsz)
    rs_nbl(rs)%glrsvec(1:rs_nbl(rs)%glralsz,:) = tmp_nbl(tpi)%svl(1:rs_nbl(rs)%glralsz,:)
! temporary index list
  else if (which.eq.7) then
    tmp_nbl(tpi)%nbl(1:rs_nbl(rs)%tmpalsz) = rs_nbl(rs)%tmpanb(1:rs_nbl(rs)%tmpalsz)
    deallocate(rs_nbl(rs)%tmpanb)
    rs_nbl(rs)%tmpalsz = min(nseq,rs_nbl(rs)%tmpalsz + 10)
    allocate(rs_nbl(rs)%tmpanb(rs_nbl(rs)%tmpalsz))
    rs_nbl(rs)%tmpanb(1:rs_nbl(rs)%tmpalsz) = tmp_nbl(tpi)%nbl(1:rs_nbl(rs)%tmpalsz)
  else
    write(ilog,*) 'Fatal. Called nbl_resz(...) with unkown list identifier (offending &
 &code is ',which,'). Please report this bug.'
    call fexit()
  end if
! 
  end
!
!-----------------------------------------------------------------------------
!
! with lrel_md set to 4 or 5, the entire local concept breaks down
! hence, we have to do a significant amount of clean-up work which may be prohibitively expensive
!
  subroutine complete_nbl(rs1,cycle_frz,tpi)
!
  use cutoffs
  use mcgrid
  use atoms
  use energies
  use sequen
  use polypep
  use fyoc
  use molecule
!
  implicit none
!
  integer imol,jmol,rs1,rs2,j,jj,tpi
  RTYPE svec(3)
  logical there,cycle_frz
!
  if (use_POLAR.EQV..false.) return
!
  imol = molofrs(rs1)
  if (lrel_md.eq.4) then
    if (chgflag(rs1).EQV..false.) return
    do j=1,cglst%ncs
      rs2 = atmres(cglst%it(j))
      if (rs2.eq.disulf(rs1)) cycle ! exclude crosslinked pairs
      jmol = molofrs(rs2)
      if (rs2.le.rs1+1) cycle
      if (cycle_frz.EQV..true.) then
        if (imol.eq.jmol) then
          if (molfrzidx(imol).ge.3) cycle
          if ((molfrzidx(imol).eq.2).AND.(nchi(rs1).eq.0).AND.(nchi(rs2).eq.0)) cycle
        else
          if ((molfrzidx(imol).eq.4).AND.(molfrzidx(jmol).eq.4)) cycle
        end if
      end if
      there = .false.
      do jj=1,rs_nbl(rs1)%nnblrs
        if (rs_nbl(rs1)%nblr(jj).eq.rs2) then
          there = .true.
          exit
        end if
      end do
      if (there.EQV..true.) cycle
      do jj=1,rs_nbl(rs1)%nnbtrs
        if (rs_nbl(rs1)%nbtr(jj).eq.rs2) then
          there = .true.
          exit
        end if
      end do
      if (there.EQV..true.) cycle
      do jj=1,rs_nbl(rs1)%nnbs
        if (rs_nbl(rs1)%nb(jj).eq.rs2) then
          there = .true.
          exit
        end if
      end do
      if (there.EQV..true.) cycle
      if (use_FEG.EQV..true.) then
        do jj=1,rs_nbl(rs1)%ngnblrs
          if (rs_nbl(rs1)%gnblr(jj).eq.rs2) then
            there = .true.
            exit
          end if
        end do
        if (there.EQV..true.) cycle
        do jj=1,rs_nbl(rs1)%ngnbtrs
          if (rs_nbl(rs1)%gnbtr(jj).eq.rs2) then
            there = .true.
            exit
          end if
        end do
        if (there.EQV..true.) cycle
        do jj=1,rs_nbl(rs1)%ngnbs
          if (rs_nbl(rs1)%gnb(jj).eq.rs2) then
            there = .true.
            exit
          end if
        end do
        if (there.EQV..true.) cycle
      end if
      if (use_FEG.EQV..false.) then
        rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
        if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3,tpi)
        rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + at(rs2)%nbb + at(rs2)%nsc
        rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
        if (molofrs(rs2).eq.molofrs(rs1)) then
          rs_nbl(rs1)%lrsvec(rs_nbl(rs1)%nnblrs,:) = 0.0
        else
          call dis_bound_rs(rs1,rs2,svec)
          rs_nbl(rs1)%lrsvec(rs_nbl(rs1)%nnblrs,:) = svec(:)
        end if
      else
        if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) then
!         do nothing (at least one of the residues is fully de-coupled)
        else if (((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..false.)).OR.&
 &    ((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..true.).AND.(fegmode.eq.1))) then
          rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
          if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3,tpi)
          rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + at(rs2)%nbb + at(rs2)%nsc
          rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
          if (molofrs(rs2).eq.molofrs(rs1)) then
            rs_nbl(rs1)%lrsvec(rs_nbl(rs1)%nnblrs,:) = 0.0
          else
            call dis_bound_rs(rs1,rs2,svec)
            rs_nbl(rs1)%lrsvec(rs_nbl(rs1)%nnblrs,:) = svec(:)
          end if
        else
          rs_nbl(rs1)%ngnblrs = rs_nbl(rs1)%ngnblrs + 1
          if (rs_nbl(rs1)%ngnblrs.gt.rs_nbl(rs1)%glralsz) call nbl_resz(rs1,6,tpi)
          rs_nbl(rs1)%ngnblrats = rs_nbl(rs1)%ngnblrats + at(rs2)%nbb + at(rs2)%nsc
          rs_nbl(rs1)%gnblr(rs_nbl(rs1)%ngnblrs) = rs2
          if (molofrs(rs2).eq.molofrs(rs1)) then
            rs_nbl(rs1)%glrsvec(rs_nbl(rs1)%ngnblrs,:) = 0.0
          else
            call dis_bound_rs(rs1,rs2,svec)
            rs_nbl(rs1)%glrsvec(rs_nbl(rs1)%ngnblrs,:) = svec(:)
          end if
        end if
      end if
    end do
!
  else if (lrel_md.eq.5) then
!   let's hope this is not a plasma
    if (chgflag(rs1).EQV..true.) then
      do rs2=rs1+2,nseq
        if (rs2.eq.disulf(rs1)) cycle ! exclude crosslinked pairs
        jmol = molofrs(rs2)
        if (cycle_frz.EQV..true.) then
          if (imol.eq.jmol) then
            if (molfrzidx(imol).ge.3) cycle
            if ((molfrzidx(imol).eq.2).AND.(nchi(rs1).eq.0).AND.(nchi(rs2).eq.0)) cycle
          else
            if ((molfrzidx(imol).eq.4).AND.(molfrzidx(jmol).eq.4)) cycle
          end if
        end if
        there = .false.
        do jj=1,rs_nbl(rs1)%nnblrs
          if (rs_nbl(rs1)%nblr(jj).eq.rs2) then
            there = .true.
            exit
          end if
        end do
        if (there.EQV..true.) cycle
        do jj=1,rs_nbl(rs1)%nnbtrs
          if (rs_nbl(rs1)%nbtr(jj).eq.rs2) then
            there = .true.
            exit
          end if
        end do
        if (there.EQV..true.) cycle
        do jj=1,rs_nbl(rs1)%nnbs
          if (rs_nbl(rs1)%nb(jj).eq.rs2) then
            there = .true.
            exit
          end if
        end do
        if (there.EQV..true.) cycle
        if (use_FEG.EQV..true.) then
          do jj=1,rs_nbl(rs1)%ngnblrs
            if (rs_nbl(rs1)%gnblr(jj).eq.rs2) then
              there = .true.
              exit
            end if
          end do
          if (there.EQV..true.) cycle
          do jj=1,rs_nbl(rs1)%ngnbtrs
            if (rs_nbl(rs1)%gnbtr(jj).eq.rs2) then
              there = .true.
              exit
            end if
          end do
          if (there.EQV..true.) cycle
          do jj=1,rs_nbl(rs1)%ngnbs
            if (rs_nbl(rs1)%gnb(jj).eq.rs2) then
              there = .true.
              exit
            end if
          end do
          if (there.EQV..true.) cycle
        end if
        if (use_FEG.EQV..false.) then
          rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
          if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3,tpi)
          rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + at(rs2)%nbb + at(rs2)%nsc
          rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
          if (molofrs(rs2).eq.molofrs(rs1)) then
            rs_nbl(rs1)%lrsvec(rs_nbl(rs1)%nnblrs,:) = 0.0
          else
            call dis_bound_rs(rs1,rs2,svec)
            rs_nbl(rs1)%lrsvec(rs_nbl(rs1)%nnblrs,:) = svec(:)
          end if
        else
          if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) then
!           do nothing (at least one of the residues is fully de-coupled)
          else if (((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..false.)).OR.&
 &      ((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..true.).AND.(fegmode.eq.1))) then
            rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
            if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3,tpi)
            rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + at(rs2)%nbb + at(rs2)%nsc
            rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
            if (molofrs(rs2).eq.molofrs(rs1)) then
              rs_nbl(rs1)%lrsvec(rs_nbl(rs1)%nnblrs,:) = 0.0
            else
              call dis_bound_rs(rs1,rs2,svec)
              rs_nbl(rs1)%lrsvec(rs_nbl(rs1)%nnblrs,:) = svec(:)
            end if
          else
            rs_nbl(rs1)%ngnblrs = rs_nbl(rs1)%ngnblrs + 1
            if (rs_nbl(rs1)%ngnblrs.gt.rs_nbl(rs1)%glralsz) call nbl_resz(rs1,6,tpi)
            rs_nbl(rs1)%ngnblrats = rs_nbl(rs1)%ngnblrats + at(rs2)%nbb + at(rs2)%nsc
            rs_nbl(rs1)%gnblr(rs_nbl(rs1)%ngnblrs) = rs2
            if (molofrs(rs2).eq.molofrs(rs1)) then
              rs_nbl(rs1)%glrsvec(rs_nbl(rs1)%ngnblrs,:) = 0.0
            else
              call dis_bound_rs(rs1,rs2,svec)
              rs_nbl(rs1)%glrsvec(rs_nbl(rs1)%ngnblrs,:) = svec(:)
            end if
          end if
        end if
      end do
    else ! residue rs1 is not charged
      do j=1,cglst%ncs
        rs2 = atmres(cglst%it(j))
        if (rs2.eq.disulf(rs1)) cycle ! exclude crosslinked pairs
        if (rs2.le.rs1+1) cycle
        jmol = molofrs(rs2)
        if (cycle_frz.EQV..true.) then
          if (imol.eq.jmol) then
            if (molfrzidx(imol).ge.3) cycle
            if ((molfrzidx(imol).eq.2).AND.(nchi(rs1).eq.0).AND.(nchi(rs2).eq.0)) cycle
          else
            if ((molfrzidx(imol).eq.4).AND.(molfrzidx(jmol).eq.4)) cycle
          end if
        end if
        there = .false.
        do jj=1,rs_nbl(rs1)%nnblrs
          if (rs_nbl(rs1)%nblr(jj).eq.rs2) then
            there = .true.
            exit
          end if
        end do
        if (there.EQV..true.) cycle
        do jj=1,rs_nbl(rs1)%nnbtrs
          if (rs_nbl(rs1)%nbtr(jj).eq.rs2) then
            there = .true.
            exit
          end if
        end do
        if (there.EQV..true.) cycle
        do jj=1,rs_nbl(rs1)%nnbs
          if (rs_nbl(rs1)%nb(jj).eq.rs2) then
            there = .true.
            exit
          end if
        end do
        if (there.EQV..true.) cycle
        if (use_FEG.EQV..true.) then
          do jj=1,rs_nbl(rs1)%ngnblrs
            if (rs_nbl(rs1)%gnblr(jj).eq.rs2) then
              there = .true.
              exit
            end if
          end do
          if (there.EQV..true.) cycle
          do jj=1,rs_nbl(rs1)%ngnbtrs
            if (rs_nbl(rs1)%gnbtr(jj).eq.rs2) then
              there = .true.
              exit
            end if
          end do
          if (there.EQV..true.) cycle
          do jj=1,rs_nbl(rs1)%ngnbs
            if (rs_nbl(rs1)%gnb(jj).eq.rs2) then
              there = .true.
              exit
            end if
          end do
          if (there.EQV..true.) cycle
        end if
        if (use_FEG.EQV..false.) then
          rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
          if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3,tpi)
          rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + at(rs2)%nbb + at(rs2)%nsc
          rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
          if (molofrs(rs2).eq.molofrs(rs1)) then
            rs_nbl(rs1)%lrsvec(rs_nbl(rs1)%nnblrs,:) = 0.0
          else
            call dis_bound_rs(rs1,rs2,svec)
            rs_nbl(rs1)%lrsvec(rs_nbl(rs1)%nnblrs,:) = svec(:)
          end if
        else
          if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) then
!           do nothing (at least one of the residues is fully de-coupled)
          else if (((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..false.)).OR.&
 &      ((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..true.).AND.(fegmode.eq.1))) then
            rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
            if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3,tpi)
            rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + at(rs2)%nbb + at(rs2)%nsc
            rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
            if (molofrs(rs2).eq.molofrs(rs1)) then
              rs_nbl(rs1)%lrsvec(rs_nbl(rs1)%nnblrs,:) = 0.0
            else
              call dis_bound_rs(rs1,rs2,svec)
              rs_nbl(rs1)%lrsvec(rs_nbl(rs1)%nnblrs,:) = svec(:)
            end if
          else
            rs_nbl(rs1)%ngnblrs = rs_nbl(rs1)%ngnblrs + 1
            if (rs_nbl(rs1)%ngnblrs.gt.rs_nbl(rs1)%glralsz) call nbl_resz(rs1,6,tpi)
            rs_nbl(rs1)%ngnblrats = rs_nbl(rs1)%ngnblrats + at(rs2)%nbb + at(rs2)%nsc
            rs_nbl(rs1)%gnblr(rs_nbl(rs1)%ngnblrs) = rs2
            if (molofrs(rs2).eq.molofrs(rs1)) then
              rs_nbl(rs1)%glrsvec(rs_nbl(rs1)%ngnblrs,:) = 0.0
            else
              call dis_bound_rs(rs1,rs2,svec)
              rs_nbl(rs1)%glrsvec(rs_nbl(rs1)%ngnblrs,:) = svec(:)
            end if
          end if
        end if
      end do
    end if
  end if
!
  end
!
!-----------------------------------------------------------------------------
!
  subroutine respairs_nbl(rs,cycle_frz,cycle_tab,tpi)
!
  use cutoffs
  use mcgrid
!
  implicit none
!
  integer mode,rs,tpi
  logical cycle_frz,cycle_tab
!
  mode = 1
  if (use_mcgrid.EQV..true.) mode = 2
  call setup_nblidx(rs,cycle_frz,cycle_tab,mode,tpi)
  if (use_mcgrid.EQV..true.) call complete_nbl(rs,cycle_frz,tpi)
!
  end
!
!-----------------------------------------------------------------------------
!
subroutine contactpairs(rs,cutdis)
!
  use iounit
  use sequen
  use cutoffs
!
  implicit none
!
  integer i,ii,kk,rs
  RTYPE dis,dis2,cutdis
!
  kk = refat(rs)
!
  do i=1,nseq
    ii = refat(i)
    call dis_bound(ii,kk,dis2)
    dis = sqrt(dis2)
    if (dis.lt.(cutdis+resrad(rs)+resrad(i))) then
      if (i.ge.rs) then
        rsp_mat(i,rs) = 1
      else
        rsp_mat(rs,i) = 1
      end if
    end if
  end do
!
end
!
!--------------------------------------------------------------------------------
!
subroutine clear_rsp()
!
  use cutoffs
  use sequen
!
  implicit none
!
  rsp_mat(:,:) = 0
!
end
!
!--------------------------------------------------------------------------------
!
subroutine clear_rsp2(rs)
!
  use cutoffs
  use sequen
!
  implicit none
!
  integer rs
!
  
  
  rsp_mat(rs,:) = 0
  rsp_mat(:,rs) = 0
!
end
!
!---------------------------------------------------------------------------------
!
subroutine cutoff_check()
!
  use iounit
  use atoms
  use mcgrid
  use polypep
  use sequen
  use cutoffs
!
  implicit none
!
  integer i,j,k,ji,ki,jj,kk
  RTYPE dis
!
  do i=1,nseq
    if (use_mcgrid.EQV..true.) then
      call grd_respairs(i)
    else if (use_rescrit.EQV..true.) then
      call respairs(i)
    end if
  end do
  do j=1,nseq
    do k=j,nseq
      if (rsp_mat(k,j).eq.0) then
        do ji=1,at(j)%nsc+at(j)%nbb
          if (ji.le.at(j)%nbb) then
            jj = at(j)%bb(ji)
          else
            jj = at(j)%sc(ji-at(j)%nbb)
          end if
          do ki=1,at(k)%nsc+at(k)%nbb
            if (ki.le.at(k)%nbb) then
              kk = at(k)%bb(ki)
            else
              kk = at(k)%sc(ki-at(k)%nbb)
            end if
            dis = (x(jj)-x(kk))**2 + (y(jj)-y(kk))**2 +&
 &                   (z(jj)-z(kk))**2
            if (dis.le.mcnb_cutoff2) then
              write(ilog,*) 'WARNING: Cutoff violation.'
              write(ilog,*) 'Residue: ',j,' Atom: ',ji,' /',jj
              write(ilog,*) 'Residue: ',k,' Atom: ',ki,' /',kk
              write(ilog,*) 'Distance: ',sqrt(dis),' (cutoff: ',&
 &                             mcnb_cutoff,')'
            end if
          end do
        end do
      end if
    end do
  end do
!
  
  
  call clear_rsp()
!
end
!
!---------------------------------------------------------------------------------
!
