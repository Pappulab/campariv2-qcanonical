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
! CONTRIBUTIONS: Rohit Pappu, Adam Steffen, Hoang Tran                     !
!                                                                          !
!--------------------------------------------------------------------------!
!
!-----------------------------------------------------------------------
!
! a set of trivial routines which copy subsets of atomic coordinates
! from or to backup arrays
!
!-----------------------------------------------------------------------
!
!
subroutine makeref_forrotlst(ati,alC,rsi,rsf)
!
  use atoms
  use zmatrix
  use sequen
  use molecule
  use atoms
!
  implicit none
!
  integer ati,i,ii,ik,k,imol,rsi,rsf
  logical alC
!
  if ((izrot(ati)%alsz.le.0).OR.(allocated(izrot(ati)%rotis).EQV..false.)) return
!
  rsf = 0
  rsi=HUGE(rsi)
  if (alC.EQV..true.) then
    imol = molofrs(atmres(ati))
    k = atmol(imol,1)
    i = 1
    do while (k.le.atmol(imol,2))
      if (i.le.izrot(ati)%alsz) then
        if ((k.ge.izrot(ati)%rotis(i,1)).AND.(k.le.izrot(ati)%rotis(i,2))) then
          k = izrot(ati)%rotis(i,2)+1
          i = i + 1
        else
          if ((k.ne.iz(1,ati)).AND.(k.ne.(iz(2,ati)))) then
            if (atmres(k).lt.rsi) rsi=atmres(k)
            if (atmres(k).gt.rsf) rsf=atmres(k)
          end if
          xref(k) = x(k)
          yref(k) = y(k)
          zref(k) = z(k)
          k = k + 1
        end if
      else
        if ((k.ne.iz(1,ati)).AND.(k.ne.(iz(2,ati)))) then
          if (atmres(k).lt.rsi) rsi=atmres(k)
          if (atmres(k).gt.rsf) rsf=atmres(k)
        end if
        xref(k) = x(k)
        yref(k) = y(k)
        zref(k) = z(k)
        k = k + 1
      end if
    end do
  else
    do k=1,izrot(ati)%alsz
      ii = izrot(ati)%rotis(k,1)
      ik = izrot(ati)%rotis(k,2)
      xref(ii:ik) = x(ii:ik)
      yref(ii:ik) = y(ii:ik)
      zref(ii:ik) = z(ii:ik)
      if (minval(atmres(ii:ik)).lt.rsi) rsi=minval(atmres(ii:ik))
      if (maxval(atmres(ii:ik)).gt.rsf) rsf=maxval(atmres(ii:ik))
    end do
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine makeref_forsc(rs)
!
  use atoms
  use polypep
!
  implicit none
!
  integer i,rs,ii
!
  do i=1,at(rs)%nbb+at(rs)%nsc
    if (i.le.at(rs)%nbb) then
      ii = at(rs)%bb(i)
    else
      ii = at(rs)%sc(i-at(rs)%nbb)
    end if
    xref(ii) = x(ii)
    yref(ii) = y(ii)
    zref(ii) = z(ii)
  end do
!
end
!
!---------------------------------------------------------------
!
!
subroutine makeref_forbb(imol,rs)
!
  use atoms
  use polypep
  use molecule
!
  implicit none
!
  integer i,rs,imol
!
! for a pivot move all atoms 'behind' active residue are moved as well
! note that this assumes that the first atom in the residue is a "backbone" atom 
  do i=at(rs)%bb(1),atmol(imol,2)
    xref(i) = x(i)
    yref(i) = y(i)
    zref(i) = z(i)
  end do
!
end
!
!---------------------------------------------------------------
!
subroutine makeref_formol(imol)
!
  use atoms
  use molecule
!
  implicit none
!
  integer i,imol
!
! all atoms in the molecule are backed up
  do i=atmol(imol,1),atmol(imol,2)
    xref(i) = x(i)
    yref(i) = y(i)
    zref(i) = z(i)
  end do
!
end
!
!---------------------------------------------------------------
!
subroutine makeref_poly(imol)
!
  use molecule
!
  implicit none
!
  integer i,imol,j
!
  rgvref(imol) = rgv(imol)
  do i=1,3
    rgevsref(imol,i) = rgevs(imol,i)
    comref(imol,i) = com(imol,i)
  end do
  do i=1,3
    do j=1,3
      rgpcsref(imol,i,j) = rgpcs(imol,i,j)
    end do
  end do
!
end
!
!---------------------------------------------------------------
!
subroutine makeref_polym(imol)
!
  use molecule
!
  implicit none
!
  integer imol
!
!  rgvmref(imol) = rgvm(imol)
!  rgevsmref(imol,1:3) = rgevsm(imol,1:3)
  commref(imol,1:3) = comm(imol,1:3)
!  do i=1,3
!    do j=1,3
!      rgpcsmref(imol,i,j) = rgpcsm(imol,i,j)
!    end do
!  end do
!
end
!
!---------------------------------------------------------------
!
subroutine makeref_zmat_gl()
!
  use atoms
  use zmatrix
!
  implicit none
!
  ztorpr(1:n) = ztor(1:n)
  bangpr(1:n) = bang(1:n)
!
end
!
!-----------------------------------------------------------------------
!
subroutine getref_forrotlst(ati,alC,rsi,rsf)
!
  use atoms
  use zmatrix
  use sequen
  use molecule
  use atoms
!
  implicit none
!
  integer ati,i,ii,ik,k,imol,rsi,rsf
  logical alC
!
  if ((izrot(ati)%alsz.le.0).OR.(allocated(izrot(ati)%rotis).EQV..false.)) return
!
  rsf = 0
  rsi=HUGE(rsi)
  if (alC.EQV..true.) then
    imol = molofrs(atmres(ati))
    k = atmol(imol,1)
    i = 1
    do while (k.le.atmol(imol,2))
      if (i.le.izrot(ati)%alsz) then
        if ((k.ge.izrot(ati)%rotis(i,1)).AND.(k.le.izrot(ati)%rotis(i,2))) then
          k = izrot(ati)%rotis(i,2)+1
          i = i + 1
        else
          if ((k.ne.iz(1,ati)).AND.(k.ne.(iz(2,ati)))) then
            if (atmres(k).lt.rsi) rsi=atmres(k)
            if (atmres(k).gt.rsf) rsf=atmres(k)
          end if
          x(k) = xref(k)
          y(k) = yref(k)
          z(k) = zref(k)
          k = k + 1
        end if
      else
        if ((k.ne.iz(1,ati)).AND.(k.ne.(iz(2,ati)))) then
          if (atmres(k).lt.rsi) rsi=atmres(k)
          if (atmres(k).gt.rsf) rsf=atmres(k)
        end if
        x(k) = xref(k)
        y(k) = yref(k)
        z(k) = zref(k)
        k = k + 1
      end if
    end do
  else
    do k=1,izrot(ati)%alsz
      ii = izrot(ati)%rotis(k,1)
      ik = izrot(ati)%rotis(k,2)
      x(ii:ik) = xref(ii:ik)
      y(ii:ik) = yref(ii:ik)
      z(ii:ik) = zref(ii:ik)
      if (minval(atmres(ii:ik)).lt.rsi) rsi=minval(atmres(ii:ik))
      if (maxval(atmres(ii:ik)).gt.rsf) rsf=maxval(atmres(ii:ik))
    end do
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine getref_forsc(rs)
!
  use atoms
  use polypep
!
  implicit none
!
  integer i,rs,ii
!
  do i=1,at(rs)%nbb+at(rs)%nsc
    if (i.le.at(rs)%nbb) then
      ii = at(rs)%bb(i)
    else
      ii = at(rs)%sc(i-at(rs)%nbb)
    end if
    x(ii) = xref(ii)
    y(ii) = yref(ii)
    z(ii) = zref(ii)
  end do
!
end
!
!-------------------------------------------------------------------
!
subroutine getref_forbb(imol,rs)
!
  use atoms
  use polypep
  use molecule
!
  implicit none
!
  integer i,rs,imol
!
  do i=at(rs)%bb(1),atmol(imol,2)
    x(i) = xref(i)
    y(i) = yref(i)
    z(i) = zref(i)
  end do
!
end
!
!-------------------------------------------------------------------
!
subroutine getref_forbb_nuccr(imol,rs)
!
  use atoms
  use polypep
  use molecule
!
  implicit none
!
  integer i,rs,imol
!
  do i=nuci(rs,1),atmol(imol,2)
    x(i) = xref(i)
    y(i) = yref(i)
    z(i) = zref(i)
  end do
!
end
!
!-------------------------------------------------------------------
!
subroutine getref_forbb_nuccr_pre(imol,rs)
!
  use atoms
  use polypep
  use molecule
  use system
!
  implicit none
!
  integer i,rs,imol
!
! rsf-1 should never be a 5'-cap residue
  do i=nuci(rs,6),atmol(imol,2)
    x(i) = xref(i)
    y(i) = yref(i)
    z(i) = zref(i)
  end do
!
end
!
!-------------------------------------------------------------------
!
subroutine getref_formol(imol)
!
  use atoms
  use molecule
!
  implicit none
!
  integer i,imol
!
! all atoms in the molecule are restored
  do i=atmol(imol,1),atmol(imol,2)
    x(i) = xref(i)
    y(i) = yref(i)
    z(i) = zref(i)
  end do
!
end
!
!---------------------------------------------------------------
!
subroutine getref_poly(imol)
!
  use molecule
!
  implicit none
!
  integer i,imol,j
!
  rgv(imol) = rgvref(imol)
  do i=1,3
    rgevs(imol,i) = rgevsref(imol,i)
    com(imol,i) = comref(imol,i)
  end do
  do i=1,3
    do j=1,3
      rgpcs(imol,i,j) = rgpcsref(imol,i,j)
    end do
  end do
!
end
!
!---------------------------------------------------------------
!
subroutine getref_zmat_gl()
!
  use atoms
  use zmatrix
!
  implicit none
!
  ztor(1:n) = ztorpr(1:n)
  bang(1:n) = bangpr(1:n)
!
end
!
!-------------------------------------------------------------------
!
subroutine getref_foruj(imol,rs)
!
  use atoms
  use iounit
  use polypep
  use molecule
  use sequen
!
  implicit none
!
  integer i,rs,imol
!
  if (seqpolty(rs).eq.'P') then
    do i=at(rs)%bb(3),atmol(imol,2)
      x(i) = xref(i)
      y(i) = yref(i)
      z(i) = zref(i)
    end do
  else
    write(ilog,*) 'Fatal. Called getref_foruj(...) with unsupported &
 &polymer or cap type. Please report this bug.'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------
!
subroutine getref_fordjo(imol,rs)
!
  use atoms
  use iounit
  use polypep
  use molecule
  use sequen
!
  implicit none
!
  integer i,rs,imol
!
  if (seqpolty(rs).eq.'P') then
    x(at(rs)%bb(4)) = xref(at(rs)%bb(4))
    y(at(rs)%bb(4)) = yref(at(rs)%bb(4))
    z(at(rs)%bb(4)) = zref(at(rs)%bb(4))
    do i=at(rs+1)%bb(1),atmol(imol,2)
      x(i) = xref(i)
      y(i) = yref(i)
      z(i) = zref(i)
    end do
  else
    write(ilog,*) 'Fatal. Called getref_fordjo(...) with unsupported &
 &polymer or cap type. Please report this bug.'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------
!
subroutine getref_fordo(imol,rs)
!
  use atoms
  use iounit
  use polypep
  use molecule
  use sequen
!
  implicit none
!
  integer i,rs,imol
!
  if (seqpolty(rs).eq.'P') then
    do i=at(rs)%bb(2),atmol(imol,2)
      x(i) = xref(i)
      y(i) = yref(i)
      z(i) = zref(i)
    end do
  else
    write(ilog,*) 'Fatal. Called getref_fordo(...) with unsupported &
 &polymer or cap type. Please report this bug.'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------
!

