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
!-----------------------------------------------------------------------
!
!
!             ####################################
!             #                                  #
!             #  A SET OF FORCE/ENERGY ROUTINES  #
!             #  WHICH WILL USUALLY TAKE UP THE  #
!             #  BULK OF CPU TIME                #
!             #                                  #
!             ####################################
!
!
!-----------------------------------------------------------------------
!
! this first 3 routines (force_PLJ_C, force_LJ_C, force_IPP_C) are basically
! to illustrate what the vectorizable routines further below are supposed to
! do
! currently not in use ...
!
subroutine force_PLJ_C(evec,rs1,rs2,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
!
  implicit none
!
  RTYPE dvec(3),dis2,svec(3),id2,eps12,eps6,term6,term12,foin,id1
  RTYPE evec(MAXENERGYTERMS),epsik,radik,term1,ca_f(n,3)
  integer rs1,rs2,ii,kk,j
!  
  call dis_bound_rs(rs1,rs2,svec)
! loop over complete set of atom-atom interactions in distant residues
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    do kk=at(rs2)%bb(1),at(rs2)%bb(1)+at(rs2)%nbb+at(rs2)%nsc-1
!
      dvec(1) = x(kk) - x(ii) + svec(1)
      dvec(2) = y(kk) - y(ii) + svec(2)
      dvec(3) = z(kk) - z(ii) + svec(3)
      dis2 = dvec(1)**2 + dvec(2)**2 + dvec(3)**2
      id2 = 1.0/dis2
!
      if (dis2.lt.mcnb_cutoff2) then
!
        epsik = 4.0*lj_eps(attyp(ii),attyp(kk))
        radik = lj_sig(attyp(ii),attyp(kk))
        if (epsik.gt.0.0) then
          eps6 = scale_attLJ*epsik
          eps12 = scale_IPP*epsik
          term6 = ((radik*id2)**3)
          term12 = term6*term6
          evec(1) = evec(1) + eps12*term12
          evec(3) = evec(3) - eps6*term6
          term12 = eps12*12.0*term12*id2
          term6 = eps6*6.0*term6*id2
          do j=1,3
            foin = dvec(j)*(term6 - term12)
            ca_f(ii,j) = ca_f(ii,j) + foin
            ca_f(kk,j) = ca_f(kk,j) - foin
          end do
        end if
!
      end if
!
      if ((atq(ii).eq.0.0).OR.(atq(kk).eq.0.0)) cycle
      id1 = sqrt(id2)
      term1 = electric*atq(ii)*atq(kk)*scale_POLAR*id1
      evec(6) = evec(6) + term1
      do j=1,3
        foin = term1*dvec(j)*id2
        ca_f(ii,j) = ca_f(ii,j) - foin
        ca_f(kk,j) = ca_f(kk,j) + foin
      end do
!
    end do
! 
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine force_LJ_C(evec,rs1,rs2,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
!
  implicit none
!
  RTYPE dvec(3),dis2,svec(3),id2,eps12,eps6,term6,term12,foin
  RTYPE evec(MAXENERGYTERMS),epsik,radik,ca_f(n,3)
  integer rs1,rs2,ii,kk,j
!  
  call dis_bound_rs(rs1,rs2,svec)
! loop over complete set of atom-atom interactions in distant residues
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    do kk=at(rs2)%bb(1),at(rs2)%bb(1)+at(rs2)%nbb+at(rs2)%nsc-1
!
      dvec(1) = x(kk) - x(ii) + svec(1)
      dvec(2) = y(kk) - y(ii) + svec(2)
      dvec(3) = z(kk) - z(ii) + svec(3)
      dis2 = dvec(1)**2 + dvec(2)**2 + dvec(3)**2
      
      if (dis2.lt.mcnb_cutoff2) then
!
        epsik = 4.0*lj_eps(attyp(ii),attyp(kk))
        if (epsik.gt.0.0) then
          id2 = 1.0/dis2
          radik = lj_sig(attyp(ii),attyp(kk))
          eps6 = scale_attLJ*epsik
          eps12 = scale_IPP*epsik
          term6 = ((radik*id2)**3)
          term12 = term6*term6
          evec(1) = evec(1) + eps12*term12
          evec(3) = evec(3) - eps6*term6
          term12 = eps12*12.0*term12*id2
          term6 = eps6*6.0*term6*id2
          do j=1,3
            foin = dvec(j)*(term6 - term12)
            ca_f(ii,j) = ca_f(ii,j) + foin
            ca_f(kk,j) = ca_f(kk,j) - foin
          end do
        end if
!
      end if
!
    end do
! 
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine force_IPP_C(evec,rs1,rs2,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
!
  implicit none
!
  RTYPE dvec(3),dis2,svec(3),id2,eps12,term6,term12,foin
  RTYPE evec(MAXENERGYTERMS),epsik,radik,ca_f(n,3)
  integer rs1,rs2,ii,kk,j
!  
  call dis_bound_rs(rs1,rs2,svec)
! loop over complete set of atom-atom interactions in distant residues
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    do kk=at(rs2)%bb(1),at(rs2)%bb(1)+at(rs2)%nbb+at(rs2)%nsc-1
!
      dvec(1) = x(kk) - x(ii) + svec(1)
      dvec(2) = y(kk) - y(ii) + svec(2)
      dvec(3) = z(kk) - z(ii) + svec(3)
      dis2 = dvec(1)**2 + dvec(2)**2 + dvec(3)**2
!
      if (dis2.lt.mcnb_cutoff2) then
!
        epsik = 4.0*lj_eps(attyp(ii),attyp(kk))
        if (epsik.gt.0.0) then
          id2 = 1.0/dis2
          radik = lj_sig(attyp(ii),attyp(kk))
          eps12 = scale_IPP*epsik
          term6 = ((radik*id2)**3)
          term12 = term6*term6
          evec(1) = evec(1) + eps12*term12
          term12 = eps12*12.0*term12*id2
          do j=1,3
            foin = -dvec(j)*term12
            ca_f(ii,j) = ca_f(ii,j) + foin
            ca_f(kk,j) = ca_f(kk,j) - foin
          end do
        end if
!
      end if
!
    end do
! 
  end do
!
end
!
!-----------------------------------------------------------------------
!
! the next four routines are vectorizable but don't use cutoffs, and
! consequently no neighbor lists
! they are easy to set up and read, but that useful
! (Vforce_PLJ, Vforce_LJ, Vforce_FEGPLJ, Vforce_IPP)
!
subroutine Vforce_PLJ(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use fyoc
!
  implicit none
!
  integer rs1,rs2,ii,hi,lo,lodi,hidi,lodi2,hidi2,k
  RTYPE sv(3)
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(n-at(rs1+2)%bb(1)+1)
  RTYPE term6(n-at(rs1+2)%bb(1)+1)
  RTYPE epsik(n-at(rs1+2)%bb(1)+1)
  RTYPE radik(n-at(rs1+2)%bb(1)+1)
  RTYPE term12(n-at(rs1+2)%bb(1)+1)
  RTYPE foin(n-at(rs1+2)%bb(1)+1,3)
  RTYPE dvec(n-at(rs1+2)%bb(1)+1,3)
  RTYPE svec(n-at(rs1+2)%bb(1)+1,3)
  RTYPE id2(n-at(rs1+2)%bb(1)+1)
  RTYPE dis2(n-at(rs1+2)%bb(1)+1)
  RTYPE id1(n-at(rs1+2)%bb(1)+1),ca_f(n,3)
!     
  lo = 1
  do rs2=rs1+2,nseq 
    hi = lo + at(rs2)%nsc + at(rs2)%nbb - 1
    call dis_bound_rs2(rs1,rs2,sv)
    svec(lo:hi,1) = sv(1)
    svec(lo:hi,2) = sv(2)
    svec(lo:hi,3) = sv(3)
    lo = hi + 1
  end do
  if (disulf(rs1).gt.rs1+1) then
    lodi = at(disulf(rs1))%bb(1) - at(rs1+2)%bb(1) + 1
    hidi = at(disulf(rs1))%bb(1) - at(rs1+2)%bb(1) + at(disulf(rs1))%nsc + at(disulf(rs1))%nbb
    lodi2 = at(disulf(rs1))%bb(1)
    hidi2 = at(disulf(rs1))%bb(1)+ at(disulf(rs1))%nsc + at(disulf(rs1))%nbb - 1
  end if
!
  hi = n!at(rs2)%bb(1)+n-at(rs1+1)%bb(1)+1-1
  lo = at(rs1+2)%bb(1)
!
! loop over complete set of atom-atom interactions in distant residues
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
!   use vector forms
    dvec(:,1) = x(lo:hi) - x(ii) + svec(:,1)
    dvec(:,2) = y(lo:hi) - y(ii) + svec(:,2)
    dvec(:,3) = z(lo:hi) - z(ii) + svec(:,3)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    epsik(:) = lj_eps(attyp(ii),attyp(lo:hi))
    radik(:) = lj_sig(attyp(ii),attyp(lo:hi))
    term6(:) = (radik(:)*id2(:))**3
    term12(:) = 4.0*scale_IPP*epsik(:)*term6(:)*term6(:)
    term6(:) = scale_attLJ*4.0*epsik(:)*term6(:)
    evec(1) = evec(1) + sum(term12(:))
    evec(3) = evec(3) - sum(term6(:))
    term12(:) = 12.0*term12(:)
    term6(:) = 6.0*term6(:)
    id1(:) = sqrt(id2(:))
    term1(:) = electric*atq(ii)*atq(lo:hi)*scale_POLAR*id1(:)
    evec(6) = evec(6) + sum(term1)
    foin(:,1) = dvec(:,1)*id2(:)*(term1 - term6(:) + term12(:))
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    ca_f(lo:hi,1) = ca_f(lo:hi,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*id2(:)*(term1 - term6(:) + term12(:))
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    ca_f(lo:hi,2) = ca_f(lo:hi,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*id2(:)*(term1 - term6(:) + term12(:))
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    ca_f(lo:hi,3) = ca_f(lo:hi,3) + foin(:,3)
    if (disulf(rs1).gt.rs1+1) then
      evec(1) = evec(1) - (1./12.)*sum(term12(lodi:hidi))
      evec(3) = evec(3) + (1./6.)*sum(term6(lodi:hidi))
      evec(6) = evec(6) - sum(term1(lodi:hidi))
      do k=1,3
        ca_f(ii,k) = ca_f(ii,k) + sum(foin(lodi:hidi,k))
        ca_f(lodi2:hidi2,k) = ca_f(lodi2:hidi2,k) - foin(lodi:hidi,k)
      end do
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! this routine deals exclusively with ghosted interactions
!
subroutine Vforce_FEGPLJ(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use fyoc
!
  implicit none
!
  integer rs1,rs2,ii,hi,lo,hidi,lodi,hidi2,lodi2,k
  RTYPE sv(3)
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(n-at(rs1+2)%bb(1)+1)
  RTYPE term6(n-at(rs1+2)%bb(1)+1)
  RTYPE term3(n-at(rs1+2)%bb(1)+1)
  RTYPE epsik(n-at(rs1+2)%bb(1)+1)
  RTYPE radik(n-at(rs1+2)%bb(1)+1)
  RTYPE term12(n-at(rs1+2)%bb(1)+1)
  RTYPE term12p(n-at(rs1+2)%bb(1)+1)
  RTYPE term6p(n-at(rs1+2)%bb(1)+1)
  RTYPE foin(n-at(rs1+2)%bb(1)+1,3)
  RTYPE dvec(n-at(rs1+2)%bb(1)+1,3)
  RTYPE svec(n-at(rs1+2)%bb(1)+1,3)
  RTYPE id2(n-at(rs1+2)%bb(1)+1)
  RTYPE dis2(n-at(rs1+2)%bb(1)+1)
  RTYPE id1(n-at(rs1+2)%bb(1)+1)
  RTYPE d1(n-at(rs1+2)%bb(1)+1),ca_f(n,3)
!     
  lo = 1
  do rs2=rs1+2,nseq 
    hi = lo + at(rs2)%nsc + at(rs2)%nbb - 1
    call dis_bound_rs2(rs1,rs2,sv)
    svec(lo:hi,1) = sv(1)
    svec(lo:hi,2) = sv(2)
    svec(lo:hi,3) = sv(3)
    lo = hi + 1
  end do
  if (disulf(rs1).gt.rs1+1) then
    lodi = at(disulf(rs1))%bb(1) - at(rs1+2)%bb(1) + 1
    hidi = at(disulf(rs1))%bb(1) - at(rs1+2)%bb(1) + at(disulf(rs1))%nsc + at(disulf(rs1))%nbb
    lodi2 = at(disulf(rs1))%bb(1)
    hidi2 = at(disulf(rs1))%bb(1)+ at(disulf(rs1))%nsc + at(disulf(rs1))%nbb - 1
  end if
!
  hi = n!at(rs2)%bb(1)+n-at(rs1+1)%bb(1)+1-1
  lo = at(rs1+2)%bb(1)
!
! loop over complete set of atom-atom interactions in distant residues
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
!   use vector forms: first the LJ term
    dvec(:,1) = x(lo:hi) - x(ii) + svec(:,1)
    dvec(:,2) = y(lo:hi) - y(ii) + svec(:,2)
    dvec(:,3) = z(lo:hi) - z(ii) + svec(:,3)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    epsik(:) = lj_eps(attyp(ii),attyp(lo:hi))
    radik(:) = lj_sig(attyp(ii),attyp(lo:hi))
    term3 = (dis2(:)/radik(:))**3
    term12p(:) = 1.0/(term3(:) + par_FEG2(2))
    term12(:) = 4.0*par_FEG2(1)*epsik(:)*term12p(:)*term12p(:)
    term6p(:) = 1.0/(term3(:) + par_FEG2(6))
    term6(:) = par_FEG2(5)*4.0*epsik(:)*term6p(:)
    evec(1) = evec(1) + sum(term12(:))
    evec(3) = evec(3) - sum(term6(:))
    term12p(:) = 12.0*term12(:)*term12p(:)*id2(:)*term3(:)
    term6p(:) = 6.0*term6(:)*term6p(:)*id2(:)*term3(:)
!   Coulomb term
    d1(:) = sqrt(dis2(:))
    id1(:) = 1.0/(d1(:) + par_FEG2(10))
    term1(:) = electric*atq(ii)*atq(lo:hi)*par_FEG2(9)*id1(:)
    evec(6) = evec(6) + sum(term1)
    term3(:) = term1(:)*id1(:)/d1(:)
!   Force increment
    foin(:,1) = dvec(:,1)*(term3 - term6p(:) + term12p(:))
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    ca_f(lo:hi,1) = ca_f(lo:hi,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*(term3 - term6p(:) + term12p(:))
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    ca_f(lo:hi,2) = ca_f(lo:hi,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*(term3 - term6p(:) + term12p(:))
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    ca_f(lo:hi,3) = ca_f(lo:hi,3) + foin(:,3)
    if (disulf(rs1).gt.rs1+1) then
      evec(1) = evec(1) - sum(term12(lodi:hidi))
      evec(3) = evec(3) + sum(term6(lodi:hidi))
      evec(6) = evec(6) - sum(term1(lodi:hidi))
      do k=1,3
        ca_f(ii,k) = ca_f(ii,k) + sum(foin(lodi:hidi,k))
        ca_f(lodi2:hidi2,k) = ca_f(lodi2:hidi2,k) - foin(lodi:hidi,k)
      end do
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine Vforce_LJ(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use fyoc
!
  implicit none
!
  integer rs1,rs2,ii,hi,lo,hidi,lodi,hidi2,lodi2,k
  RTYPE sv(3)
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(n-at(rs1+2)%bb(1)+1)
  RTYPE term6(n-at(rs1+2)%bb(1)+1)
  RTYPE epsik(n-at(rs1+2)%bb(1)+1)
  RTYPE radik(n-at(rs1+2)%bb(1)+1)
  RTYPE term12(n-at(rs1+2)%bb(1)+1)
  RTYPE foin(n-at(rs1+2)%bb(1)+1,3)
  RTYPE dvec(n-at(rs1+2)%bb(1)+1,3)
  RTYPE svec(n-at(rs1+2)%bb(1)+1,3)
  RTYPE id2(n-at(rs1+2)%bb(1)+1)
  RTYPE dis2(n-at(rs1+2)%bb(1)+1),ca_f(n,3)
!     
  lo = 1
  do rs2=rs1+2,nseq 
    hi = lo + at(rs2)%nsc + at(rs2)%nbb - 1
    call dis_bound_rs2(rs1,rs2,sv)
    svec(lo:hi,1) = sv(1)
    svec(lo:hi,2) = sv(2)
    svec(lo:hi,3) = sv(3)
    lo = hi + 1
  end do
  if (disulf(rs1).gt.rs1+1) then
    lodi = at(disulf(rs1))%bb(1) - at(rs1+2)%bb(1) + 1
    hidi = at(disulf(rs1))%bb(1) - at(rs1+2)%bb(1) + at(disulf(rs1))%nsc + at(disulf(rs1))%nbb
    lodi2 = at(disulf(rs1))%bb(1)
    hidi2 = at(disulf(rs1))%bb(1)+ at(disulf(rs1))%nsc + at(disulf(rs1))%nbb - 1
  end if
!
  hi = n!at(rs2)%bb(1)+n-at(rs1+1)%bb(1)+1-1
  lo = at(rs1+2)%bb(1)
!
! loop over complete set of atom-atom interactions in distant residues
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
!   use vector forms
    dvec(:,1) = x(lo:hi) - x(ii) + svec(:,1)
    dvec(:,2) = y(lo:hi) - y(ii) + svec(:,2)
    dvec(:,3) = z(lo:hi) - z(ii) + svec(:,3)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    epsik(:) = lj_eps(attyp(ii),attyp(lo:hi))
    radik(:) = lj_sig(attyp(ii),attyp(lo:hi))
    term6(:) = (radik(:)*id2(:))**3
    term12(:) = 4.0*scale_IPP*epsik(:)*term6(:)*term6(:)
    term6(:) = scale_attLJ*4.0*epsik(:)*term6(:)
    evec(1) = evec(1) + sum(term12(:))
    evec(3) = evec(3) - sum(term6(:))
    term1(:) = id2(:)*(12.0*term12(:) - 6.0*term6(:))
    foin(:,1) = dvec(:,1)*term1(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    ca_f(lo:hi,1) = ca_f(lo:hi,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*term1(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    ca_f(lo:hi,2) = ca_f(lo:hi,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*term1(:)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    ca_f(lo:hi,3) = ca_f(lo:hi,3) + foin(:,3)
    if (disulf(rs1).gt.rs1+1) then
      evec(1) = evec(1) - sum(term12(lodi:hidi))
      evec(3) = evec(3) + sum(term6(lodi:hidi))
      do k=1,3
        ca_f(ii,k) = ca_f(ii,k) + sum(foin(lodi:hidi,k))
        ca_f(lodi2:hidi2,k) = ca_f(lodi2:hidi2,k) - foin(lodi:hidi,k)
      end do
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine Vforce_IPP(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use fyoc
!
  implicit none
!
  integer rs1,rs2,ii,hi,lo,hidi,lodi,hidi2,lodi2,k
  RTYPE sv(3)
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term6(n-at(rs1+2)%bb(1)+1)
  RTYPE epsik(n-at(rs1+2)%bb(1)+1)
  RTYPE radik(n-at(rs1+2)%bb(1)+1)
  RTYPE term12(n-at(rs1+2)%bb(1)+1)
  RTYPE foin(n-at(rs1+2)%bb(1)+1,3)
  RTYPE dvec(n-at(rs1+2)%bb(1)+1,3)
  RTYPE svec(n-at(rs1+2)%bb(1)+1,3)
  RTYPE id2(n-at(rs1+2)%bb(1)+1)
  RTYPE dis2(n-at(rs1+2)%bb(1)+1),ca_f(n,3)

  lo = 1
  do rs2=rs1+2,nseq 
    hi = lo + at(rs2)%nsc + at(rs2)%nbb - 1
    call dis_bound_rs2(rs1,rs2,sv)
    svec(lo:hi,1) = sv(1)
    svec(lo:hi,2) = sv(2)
    svec(lo:hi,3) = sv(3)
    lo = hi + 1
  end do
  if (disulf(rs1).gt.rs1+1) then
    lodi = at(disulf(rs1))%bb(1) - at(rs1+2)%bb(1) + 1
    hidi = at(disulf(rs1))%bb(1) - at(rs1+2)%bb(1) + at(disulf(rs1))%nsc + at(disulf(rs1))%nbb
    lodi2 = at(disulf(rs1))%bb(1)
    hidi2 = at(disulf(rs1))%bb(1)+ at(disulf(rs1))%nsc + at(disulf(rs1))%nbb - 1
  end if
!
  hi = n!at(rs2)%bb(1)+n-at(rs1+1)%bb(1)+1-1
  lo = at(rs1+2)%bb(1)
!
! loop over complete set of atom-atom interactions in distant residues
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
!   use vector forms
    dvec(:,1) = x(lo:hi) - x(ii) + svec(:,1)
    dvec(:,2) = y(lo:hi) - y(ii) + svec(:,2)
    dvec(:,3) = z(lo:hi) - z(ii) + svec(:,3)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    epsik(:) = lj_eps(attyp(ii),attyp(lo:hi))
    radik(:) = lj_sig(attyp(ii),attyp(lo:hi))
    term6(:) = (radik(:)*id2(:))**3
    term12(:) = 4.0*scale_IPP*epsik(:)*term6(:)*term6(:)
    evec(1) = evec(1) + sum(term12(:))
    term6(:) = id2(:)*12.0*term12(:)
    foin(:,1) = dvec(:,1)*term6(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    ca_f(lo:hi,1) = ca_f(lo:hi,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*term6(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    ca_f(lo:hi,2) = ca_f(lo:hi,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*term6(:)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    ca_f(lo:hi,3) = ca_f(lo:hi,3) + foin(:,3)
    if (disulf(rs1).gt.rs1+1) then
      evec(1) = evec(1) - sum(term12(lodi:hidi))
      do k=1,3
        ca_f(ii,k) = ca_f(ii,k) + sum(foin(lodi:hidi,k))
        ca_f(lodi2:hidi2,k) = ca_f(lodi2:hidi2,k) - foin(lodi:hidi,k)
      end do
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! the following are the most important routines for standard MD/LD calculations
! these neighbor-list based fxns are not entirely consistent with the
! residue-based functions for the following reason:
! for short-range potentials the residue-based fxns can afford to 
! build in a second check for cutoffs, i.e., residue-consistency is
! sacrificed in favor of i) efficiency and ii) rigorous spherical cutoffs.
! here, we have no easy way of re-filtering the neighbor-list without building
! in conditionals which would partially break the vectorizability.
! hence, a larger subset of interactions is calculated here.
! note that for Coulomb interactions, this is currently not an issue, since
! the residue-based function do NOT build in secondary filtering for those.
! one solution would be purely atom-based neighbor-lists with the downside of
! having to maintain multiple and larger lists, and having to write a much
! more detailed neighbor-list routine (see fmcscgrid.f).
!
subroutine Vforce_PLJ_C(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use system
  use molecule
!
  implicit none
!
  integer rs1,rs2,ii,j,hi,lo,k,i,shii
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%nnbats)
  RTYPE term0(rs_nbl(rs1)%nnbats)
  RTYPE term6(rs_nbl(rs1)%nnbats)
  RTYPE term6p(rs_nbl(rs1)%nnbats)
  RTYPE term12(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE id1(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  RTYPE epsik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbats)
  RTYPE radik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbats)
  integer sh(rs_nbl(rs1)%nnbs),hira,lora
  integer k1(rs_nbl(rs1)%nnbs,3),k2(rs_nbl(rs1)%nnbs,3)
  RTYPE drav(rs_nbl(rs1)%nnbs,3),srav(rs_nbl(rs1)%nnbs,3)
  RTYPE ca_f(n,3)

! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
  ii = refat(rs1)
  drav(:,1) = x(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - x(ii)
  drav(:,2) = y(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - y(ii)
  drav(:,3) = z(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - z(ii)
  hira = rs_nbl(rs1)%nnbs
  lora = 1
!
! PBC
  if ((bnd_type.eq.1).AND.(lora.le.hira)) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(j)-drav(lora:hira,j))/bnd_params(j))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(j))/bnd_params(j))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      srav(lora:hira,1:2) = 0.0
      do j=3,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(6)-drav(lora:hira,j))/bnd_params(6))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(6))/bnd_params(6))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(6)
      end do
    end if
  else 
    srav(:,:) = 0.0
  end if
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))%nbb + at(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))%nsc - 1
  shii = at(rs1)%nbb + at(rs1)%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + sh(k)
    if ((rs2.gt.rsmol(molofrs(rs1),2)).OR.(rs2.lt.rsmol(molofrs(rs1),1))) then
      dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,1)
      dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,2)
      dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,3)
    else
      dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
      dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
      dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    end if
    term0(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    do ii=lo,hi
      epsik(:,ii) = lj_eps(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
      radik(:,ii) = lj_sig(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
    end do
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    id1(:) = sqrt(id2(:))
    term6p(:) = (radik(i,:)*id2(:))**3
    term12(:) = 4.0*scale_IPP*epsik(i,:)*term6p(:)*term6p(:)
    term6(:) = scale_attLJ*4.0*epsik(i,:)*term6p(:)
    evec(1) = evec(1) + sum(term12(:))
    evec(3) = evec(3) - sum(term6(:))
    term12(:) = 12.0*term12(:)*id2(:)
    term6(:) = 6.0*term6(:)*id2(:)
!   Coulomb potential
    term1(:) = electric*atq(ii)*term0(:)*scale_POLAR*id1(:)
    evec(6) = evec(6) + sum(term1)
    foin(:,1) = dvec(:,1)*(term1*id2(:) - term6(:) + term12(:))
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*(term1*id2(:) - term6(:) + term12(:))
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*(term1*id2(:) - term6(:) + term12(:))
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
!   virial
!    if (pflag.EQV..true.) then
!      ens%insVirT(1,1) = ens%insVirT(1,1) + sum(foin(:,1)*dvec(:,1))
!      ens%insVirT(1,2) = ens%insVirT(1,2) + sum(foin(:,1)*dvec(:,2))
!      ens%insVirT(1,3) = ens%insVirT(1,3) + sum(foin(:,1)*dvec(:,3))
!      ens%insVirT(2,1) = ens%insVirT(2,1) + sum(foin(:,2)*dvec(:,1))
!      ens%insVirT(2,2) = ens%insVirT(2,2) + sum(foin(:,2)*dvec(:,2))
!      ens%insVirT(2,3) = ens%insVirT(2,3) + sum(foin(:,2)*dvec(:,3))
!      ens%insVirT(3,1) = ens%insVirT(3,1) + sum(foin(:,3)*dvec(:,1))
!      ens%insVirT(3,2) = ens%insVirT(3,2) + sum(foin(:,3)*dvec(:,2))
!      ens%insVirT(3,3) = ens%insVirT(3,3) + sum(foin(:,3)*dvec(:,3))
!    end if
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + sh(k)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
!
! the same routine modified for Ewald electrostatics (the real-space part to be precise)
!
subroutine Vforce_PEWLJ_C(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use ewalds
  use sequen
  use system
  use molecule
!
  implicit none
!
  integer rs1,rs2,ii,j,hi,lo,k,i,shii
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%nnbats)
  RTYPE term2(rs_nbl(rs1)%nnbats)
  RTYPE term0(rs_nbl(rs1)%nnbats)
  RTYPE term6(rs_nbl(rs1)%nnbats)
  RTYPE term6p(rs_nbl(rs1)%nnbats)
  RTYPE term12(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE id1(rs_nbl(rs1)%nnbats)
  RTYPE d1(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  RTYPE epsik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbats)
  RTYPE radik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbats)
  integer sh(rs_nbl(rs1)%nnbs),hira,lora
  integer k1(rs_nbl(rs1)%nnbs,3),k2(rs_nbl(rs1)%nnbs,3)
  RTYPE drav(rs_nbl(rs1)%nnbs,3),srav(rs_nbl(rs1)%nnbs,3)
  RTYPE ca_f(n,3)

! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
  ii = refat(rs1)
  drav(:,1) = x(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - x(ii)
  drav(:,2) = y(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - y(ii)
  drav(:,3) = z(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - z(ii)
  hira = rs_nbl(rs1)%nnbs
  lora = 1
!
! PBC
  if ((bnd_type.eq.1).AND.(lora.le.hira)) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(j)-drav(lora:hira,j))/bnd_params(j))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(j))/bnd_params(j))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      srav(lora:hira,1:2) = 0.0
      do j=3,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(6)-drav(lora:hira,j))/bnd_params(6))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(6))/bnd_params(6))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(6)
      end do
    end if
  else 
    srav(:,:) = 0.0
  end if
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))%nbb + at(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))%nsc - 1
  shii = at(rs1)%nbb + at(rs1)%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + sh(k)
    if ((rs2.gt.rsmol(molofrs(rs1),2)).OR.(rs2.lt.rsmol(molofrs(rs1),1))) then
      dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,1)
      dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,2)
      dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,3)
    else
      dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
      dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
      dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    end if
    term0(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    do ii=lo,hi
      epsik(:,ii) = lj_eps(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
      radik(:,ii) = lj_sig(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
    end do
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    term6p(:) = (radik(i,:)*id2(:))**3
    term12(:) = 4.0*scale_IPP*epsik(i,:)*term6p(:)*term6p(:)
    term6(:) = scale_attLJ*4.0*epsik(i,:)*term6p(:)
    evec(1) = evec(1) + sum(term12(:))
    evec(3) = evec(3) - sum(term6(:))
    term12(:) = 12.0*term12(:)
    term6(:) = 6.0*term6(:)
!   Coulomb potential
    d1(:) = sqrt(dis2(:))
    id1(:) = 1.0/d1(:)
    term1(:) = electric*atq(ii)*term0(:)*scale_POLAR*id1(:)
    term2(:) = term1*(1.0 - erf(ewpm*d1(:)))
    evec(6) = evec(6) + sum(term2(:))
    term2(:) = id2(:)*(term2(:) - term6(:) + term12(:) + &
 &  term1(:)*ewpite*d1(:)*exp(-ewpm2*dis2(:)))
    foin(:,1) = dvec(:,1)*term2(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*term2(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*term2(:)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + sh(k)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
!
! the same routine for (G)RF electrostatics
!
subroutine Vforce_PRFLJ_C(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use math
  use system
  use sequen
  use molecule
!
  implicit none
!
  integer rs1,rs2,ii,j,hi,lo,k,i,shii
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%nnbats)
  RTYPE term2(rs_nbl(rs1)%nnbats)
  integer isout(rs_nbl(rs1)%nnbats)
  RTYPE term0(rs_nbl(rs1)%nnbats)
  RTYPE term6(rs_nbl(rs1)%nnbats)
  RTYPE term6p(rs_nbl(rs1)%nnbats)
  RTYPE term12(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE d2diff(rs_nbl(rs1)%nnbats)
  RTYPE id1(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  RTYPE epsik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbats)
  RTYPE radik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbats)
  integer sh(rs_nbl(rs1)%nnbs),hira,lora
  integer k1(rs_nbl(rs1)%nnbs,3),k2(rs_nbl(rs1)%nnbs,3)
  RTYPE drav(rs_nbl(rs1)%nnbs,3),srav(rs_nbl(rs1)%nnbs,3)
  RTYPE ca_f(n,3)

! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
  ii = refat(rs1)
  drav(:,1) = x(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - x(ii)
  drav(:,2) = y(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - y(ii)
  drav(:,3) = z(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - z(ii)
  hira = rs_nbl(rs1)%nnbs
  lora = 1
!
! PBC
  if ((bnd_type.eq.1).AND.(lora.le.hira)) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(j)-drav(lora:hira,j))/bnd_params(j))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(j))/bnd_params(j))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      srav(lora:hira,1:2) = 0.0
      do j=3,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(6)-drav(lora:hira,j))/bnd_params(6))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(6))/bnd_params(6))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(6)
      end do
    end if
  else 
    srav(:,:) = 0.0
  end if
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))%nbb + at(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))%nsc - 1
  shii = at(rs1)%nbb + at(rs1)%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + sh(k)
    if ((rs2.gt.rsmol(molofrs(rs1),2)).OR.(rs2.lt.rsmol(molofrs(rs1),1))) then
      dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,1)
      dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,2)
      dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,3)
    else
      dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
      dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
      dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    end if
    term0(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    do ii=lo,hi
      epsik(:,ii) = lj_eps(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
      radik(:,ii) = lj_sig(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
    end do
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    id1(:) = sqrt(id2(:))
    term6p(:) = (radik(i,:)*id2(:))**3
    term12(:) = 4.0*scale_IPP*epsik(i,:)*term6p(:)*term6p(:)
    term6(:) = scale_attLJ*4.0*epsik(i,:)*term6p(:)
    evec(1) = evec(1) + sum(term12(:))
    evec(3) = evec(3) - sum(term6(:))
    term12(:) = 12.0*term12(:)
    term6(:) = 6.0*term6(:)
!   Coulomb potential
    term1(:) = electric*atq(ii)*term0(:)*scale_POLAR
!   there is no general minimal cost trick, it entirely depends on how many 
!   atoms in the neighbor list are in fact outside of the cutoff (and hence
!   the cleanest fix would be to solve the issue in NBlist-generation, but
!   that's not recommended for many reasons).
!   here, we're going to zero out term1 for those beyond the cutoff
    d2diff(:) = max((mcel_cutoff2-dis2(:))*imcel2,0.0d0)
    isout(:) = ceiling(d2diff(:))
    term1(:) = term1(:)*dble(isout(:))
    term2(:) = term1(:)*(id1(:)+par_POLAR(1)*dis2(:)-par_POLAR(2))
    evec(6) = evec(6) + sum(term2)
    term2(:) = id2(:)*(term1(:)*(id1(:)-2.0*par_POLAR(1)*dis2(:))&
 &                                   - term6(:) + term12(:))
    foin(:,1) = dvec(:,1)*term2(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*term2(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*term2(:)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + sh(k)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
!
! the same routine for ghosted interactions
!
subroutine Vforce_FEGPLJ_C(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use system
  use sequen
  use molecule
!
  implicit none
!
  integer rs1,rs2,ii,j,hi,lo,k,i,shii
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%ngnbats)
  RTYPE term3(rs_nbl(rs1)%ngnbats)
  RTYPE term0(rs_nbl(rs1)%ngnbats)
  RTYPE term6(rs_nbl(rs1)%ngnbats)
  RTYPE term6p(rs_nbl(rs1)%ngnbats)
  RTYPE term12(rs_nbl(rs1)%ngnbats)
  RTYPE term12p(rs_nbl(rs1)%ngnbats)
  RTYPE d1(rs_nbl(rs1)%ngnbats)
  RTYPE foin(rs_nbl(rs1)%ngnbats,3)
  RTYPE dveci(rs_nbl(rs1)%ngnbats,3)
  RTYPE dvec(rs_nbl(rs1)%ngnbats,3)
  RTYPE id2(rs_nbl(rs1)%ngnbats)
  RTYPE dis2(rs_nbl(rs1)%ngnbats)
  RTYPE id1(rs_nbl(rs1)%ngnbats)
  RTYPE for_k(rs_nbl(rs1)%ngnbats,3)
  RTYPE epsik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%ngnbats)
  RTYPE radik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%ngnbats)
  integer sh(rs_nbl(rs1)%ngnbs),hira,lora
  integer k1(rs_nbl(rs1)%ngnbs,3),k2(rs_nbl(rs1)%ngnbs,3)
  RTYPE drav(rs_nbl(rs1)%ngnbs,3),srav(rs_nbl(rs1)%ngnbs,3)
  RTYPE ca_f(n,3)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
  ii = refat(rs1)
  drav(:,1) = x(refat(rs_nbl(rs1)%gnb(1:rs_nbl(rs1)%ngnbs))) - x(ii)
  drav(:,2) = y(refat(rs_nbl(rs1)%gnb(1:rs_nbl(rs1)%ngnbs))) - y(ii)
  drav(:,3) = z(refat(rs_nbl(rs1)%gnb(1:rs_nbl(rs1)%ngnbs))) - z(ii)
  hira = rs_nbl(rs1)%ngnbs
  lora = 1
!
! PBC
  if ((bnd_type.eq.1).AND.(lora.le.hira)) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(j)-drav(lora:hira,j))/bnd_params(j))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(j))/bnd_params(j))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      srav(lora:hira,1:2) = 0.0
      do j=3,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(6)-drav(lora:hira,j))/bnd_params(6))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(6))/bnd_params(6))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(6)
      end do
    end if
  else 
    srav(:,:) = 0.0
  end if
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%gnb(1:rs_nbl(rs1)%ngnbs))%nbb + at(rs_nbl(rs1)%gnb(1:rs_nbl(rs1)%ngnbs))%nsc - 1
  shii = at(rs1)%nbb + at(rs1)%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%ngnbs
    rs2 = rs_nbl(rs1)%gnb(k)
    hi = lo + sh(k)
    if ((rs2.gt.rsmol(molofrs(rs1),2)).OR.(rs2.lt.rsmol(molofrs(rs1),1))) then
      dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,1)
      dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,2)
      dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,3)
    else
      dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
      dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
      dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    end if
    term0(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    do ii=lo,hi
      epsik(:,ii) = lj_eps(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
      radik(:,ii) = lj_sig(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
    end do
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
!
!   ghosted LJ
    term3 = (dis2(:)/radik(i,:))**3
    term12p(:) = 1.0/(term3(:) + par_FEG2(2))
    term12(:) = 4.0*par_FEG2(1)*epsik(i,:)*term12p(:)*term12p(:)
    term6p(:) = 1.0/(term3(:) + par_FEG2(6))
    term6(:) = par_FEG2(5)*4.0*epsik(i,:)*term6p(:)
    evec(1) = evec(1) + sum(term12(:))
    evec(3) = evec(3) - sum(term6(:))
    term12(:) = 12.0*term12(:)*term12p(:)*id2(:)*term3(:)
    term6(:) = 6.0*term6(:)*term6p(:)*id2(:)*term3(:)
!
!   ghosted Coulomb potential
    d1(:) = sqrt(dis2(:))
    id1(:) = 1.0/(d1(:) + par_FEG2(10))
    term1(:) = electric*atq(ii)*term0(:)*par_FEG2(9)*id1(:)
    evec(6) = evec(6) + sum(term1)
    term1(:) = term1(:)*id1(:)/d1(:)
!
    foin(:,1) = dvec(:,1)*(term1 - term6(:) + term12(:))
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*(term1 - term6(:) + term12(:))
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*(term1 - term6(:) + term12(:))
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%ngnbs
    rs2 = rs_nbl(rs1)%gnb(k)
    hi = lo + sh(k)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
!
! the same routine for the combination of ghosted interactions and (G)RF electrostatics
! assuming simple linear scaling for the Cb-ghosting
!
subroutine Vforce_FEGPRFLJ_C(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use math
  use cutoffs
  use sequen
  use system
  use molecule
!
  implicit none
!
  integer rs1,rs2,ii,j,hi,lo,k,i,shii
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%ngnbats)
  RTYPE term2(rs_nbl(rs1)%ngnbats)
  integer isout(rs_nbl(rs1)%ngnbats)
  RTYPE term3(rs_nbl(rs1)%ngnbats)
  RTYPE term0(rs_nbl(rs1)%ngnbats)
  RTYPE term6(rs_nbl(rs1)%ngnbats)
  RTYPE term6p(rs_nbl(rs1)%ngnbats)
  RTYPE term12(rs_nbl(rs1)%ngnbats)
  RTYPE term12p(rs_nbl(rs1)%ngnbats)
  RTYPE d1(rs_nbl(rs1)%ngnbats)
  RTYPE foin(rs_nbl(rs1)%ngnbats,3)
  RTYPE dveci(rs_nbl(rs1)%ngnbats,3)
  RTYPE dvec(rs_nbl(rs1)%ngnbats,3)
  RTYPE id2(rs_nbl(rs1)%ngnbats)
  RTYPE dis2(rs_nbl(rs1)%ngnbats)
  RTYPE d2diff(rs_nbl(rs1)%ngnbats)
  RTYPE id1(rs_nbl(rs1)%ngnbats)
  RTYPE for_k(rs_nbl(rs1)%ngnbats,3)
  RTYPE epsik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%ngnbats)
  RTYPE radik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%ngnbats)
  integer sh(rs_nbl(rs1)%ngnbs),hira,lora
  integer k1(rs_nbl(rs1)%ngnbs,3),k2(rs_nbl(rs1)%ngnbs,3)
  RTYPE drav(rs_nbl(rs1)%ngnbs,3),srav(rs_nbl(rs1)%ngnbs,3)
  RTYPE ca_f(n,3)
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
  ii = refat(rs1)
  drav(:,1) = x(refat(rs_nbl(rs1)%gnb(1:rs_nbl(rs1)%ngnbs))) - x(ii)
  drav(:,2) = y(refat(rs_nbl(rs1)%gnb(1:rs_nbl(rs1)%ngnbs))) - y(ii)
  drav(:,3) = z(refat(rs_nbl(rs1)%gnb(1:rs_nbl(rs1)%ngnbs))) - z(ii)
  hira = rs_nbl(rs1)%ngnbs
  lora = 1
!
! PBC
  if ((bnd_type.eq.1).AND.(lora.le.hira)) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(j)-drav(lora:hira,j))/bnd_params(j))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(j))/bnd_params(j))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      srav(lora:hira,1:2) = 0.0
      do j=3,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(6)-drav(lora:hira,j))/bnd_params(6))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(6))/bnd_params(6))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(6)
      end do
    end if
  else 
    srav(:,:) = 0.0
  end if
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%gnb(1:rs_nbl(rs1)%ngnbs))%nbb + at(rs_nbl(rs1)%gnb(1:rs_nbl(rs1)%ngnbs))%nsc - 1
  shii = at(rs1)%nbb + at(rs1)%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%ngnbs
    rs2 = rs_nbl(rs1)%gnb(k)
    hi = lo + sh(k)
    if ((rs2.gt.rsmol(molofrs(rs1),2)).OR.(rs2.lt.rsmol(molofrs(rs1),1))) then
      dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,1)
      dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,2)
      dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,3)
    else
      dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
      dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
      dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    end if
    term0(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    do ii=lo,hi
      epsik(:,ii) = lj_eps(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
      radik(:,ii) = lj_sig(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
    end do
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
!
!   ghosted LJ
    term3 = (dis2(:)/radik(i,:))**3
    term12p(:) = 1.0/(term3(:) + par_FEG2(2))
    term12(:) = 4.0*par_FEG2(1)*epsik(i,:)*term12p(:)*term12p(:)
    term6p(:) = 1.0/(term3(:) + par_FEG2(6))
    term6(:) = par_FEG2(5)*4.0*epsik(i,:)*term6p(:)
    evec(1) = evec(1) + sum(term12(:))
    evec(3) = evec(3) - sum(term6(:))
    term12(:) = 12.0*term12(:)*term12p(:)*term3(:)
    term6(:) = 6.0*term6(:)*term6p(:)*term3(:)
!
!   ghosted Coulomb with RF potential: this only works correctly for fegcbmode = 1
    d1(:) = sqrt(dis2(:))
    id1(:) = 1.0/(d1(:))
    term1(:) = electric*atq(ii)*term0(:)*par_FEG2(9)
!   there is no general minimal cost trick, it entirely depends on how many 
!   atoms in the neighbor list are in fact outside of the cutoff (and hence
!   the cleanest fix would be to solve the issue in NBlist-generation, but
!   that's not recommended for many reasons).
!   here, we're going to zero out term1 for those beyond the cutoff
    d2diff(:) = max((mcel_cutoff2-dis2(:))*imcel2,0.0d0)
    isout(:) = ceiling(d2diff(:))
    term1(:) = term1(:)*dble(isout(:))
    term2(:) = term1(:)*(id1(:)+par_POLAR(1)*dis2(:)-par_POLAR(2))
    evec(6) = evec(6) + sum(term2)
    term2(:) = id2(:)*(term1(:)*(id1(:)-2.0*par_POLAR(1)*dis2(:))&
 &                                   - term6(:) + term12(:))
    foin(:,1) = dvec(:,1)*term2(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*term2(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*term2(:)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%ngnbs
    rs2 = rs_nbl(rs1)%gnb(k)
    hi = lo + sh(k)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
!
! just the LJ potential
!
subroutine Vforce_LJ_C(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use system
  use sequen
  use molecule
!
  implicit none
!
  integer rs1,rs2,ii,j,hi,lo,k,i,shii
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term6(rs_nbl(rs1)%nnbats)
  RTYPE term6p(rs_nbl(rs1)%nnbats)
  RTYPE term12(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  RTYPE epsik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbats)
  RTYPE radik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbats)
  integer sh(rs_nbl(rs1)%nnbs),hira,lora
  integer k1(rs_nbl(rs1)%nnbs,3),k2(rs_nbl(rs1)%nnbs,3)
  RTYPE drav(rs_nbl(rs1)%nnbs,3),srav(rs_nbl(rs1)%nnbs,3),ca_f(n,3)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
  ii = refat(rs1)
  drav(:,1) = x(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - x(ii)
  drav(:,2) = y(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - y(ii)
  drav(:,3) = z(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - z(ii)
  hira = rs_nbl(rs1)%nnbs
  lora = 1
!
! PBC
  if ((bnd_type.eq.1).AND.(lora.le.hira)) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(j)-drav(lora:hira,j))/bnd_params(j))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(j))/bnd_params(j))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      srav(lora:hira,1:2) = 0.0
      do j=3,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(6)-drav(lora:hira,j))/bnd_params(6))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(6))/bnd_params(6))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(6)
      end do
    end if
  else 
    srav(:,:) = 0.0
  end if
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh(:) = at(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))%nbb + at(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))%nsc - 1
  shii = at(rs1)%nbb + at(rs1)%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + sh(k)
    if ((rs2.gt.rsmol(molofrs(rs1),2)).OR.(rs2.lt.rsmol(molofrs(rs1),1))) then
      dveci(lo:hi,1) = x(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k))) + srav(k,1)
      dveci(lo:hi,2) = y(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k))) + srav(k,2)
      dveci(lo:hi,3) = z(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k))) + srav(k,3)
    else
      dveci(lo:hi,1) = x(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k)))
      dveci(lo:hi,2) = y(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k)))
      dveci(lo:hi,3) = z(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k)))
    end if
    do ii=lo,hi
      epsik(:,ii) = lj_eps(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
      radik(:,ii) = lj_sig(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
    end do
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    term6p(:) = (radik(i,:)*id2(:))**3
    term12(:) = 4.0*scale_IPP*epsik(i,:)*term6p(:)*term6p(:)
    term6(:) = scale_attLJ*4.0*epsik(i,:)*term6p(:)
    evec(1) = evec(1) + sum(term12(:))
    evec(3) = evec(3) - sum(term6(:))
    term6(:) = id2(:)*(12.0*term12(:) - 6.0*term6(:))
    foin(:,1) = dvec(:,1)*term6(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*term6(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*term6(:)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
!   virial
!    if (pflag.EQV..true.) then
!      ens%insVirT(1,1) = ens%insVirT(1,1) + sum(foin(:,1)*dvec(:,1))
!      ens%insVirT(2,2) = ens%insVirT(2,2) + sum(foin(:,2)*dvec(:,2))
!      ens%insVirT(3,3) = ens%insVirT(3,3) + sum(foin(:,3)*dvec(:,3))
!    end if
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + sh(k)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!
!-----------------------------------------------------------------------
!
! just excluded volume interactions
! 
subroutine Vforce_IPP_C(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use system
  use molecule
!
  implicit none
!
  integer rs1,rs2,ii,j,hi,lo,k,i,shii
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term6p(rs_nbl(rs1)%nnbats)
  RTYPE term12(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  RTYPE epsik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbats)
  RTYPE radik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbats)
  integer sh(rs_nbl(rs1)%nnbs),hira,lora
  integer k1(rs_nbl(rs1)%nnbs,3),k2(rs_nbl(rs1)%nnbs,3)
  RTYPE drav(rs_nbl(rs1)%nnbs,3),srav(rs_nbl(rs1)%nnbs,3),ca_f(n,3)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
  ii = refat(rs1)
  drav(:,1) = x(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - x(ii)
  drav(:,2) = y(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - y(ii)
  drav(:,3) = z(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - z(ii)
  hira = rs_nbl(rs1)%nnbs
  lora = 1
!
! PBC
  if ((bnd_type.eq.1).AND.(lora.le.hira)) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(j)-drav(lora:hira,j))/bnd_params(j))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(j))/bnd_params(j))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      srav(lora:hira,1:2) = 0.0
      do j=3,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(6)-drav(lora:hira,j))/bnd_params(6))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(6))/bnd_params(6))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(6)
      end do
    end if
  else 
    srav(:,:) = 0.0
  end if
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))%nbb + at(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))%nsc - 1
  shii = at(rs1)%nbb + at(rs1)%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + sh(k)
    if ((rs2.gt.rsmol(molofrs(rs1),2)).OR.(rs2.lt.rsmol(molofrs(rs1),1))) then
      dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,1)
      dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,2)
      dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + srav(k,3)
    else
      dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
      dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
      dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    end if
    do ii=lo,hi
      epsik(:,ii) = lj_eps(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
      radik(:,ii) = lj_sig(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
    end do
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!    write(*,*) dveci(1,1),x(ii),dveci(1,2)
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
!    write(*,*) 'res: ',dvec(1,1),dvec(1,2)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    term6p(:) = (radik(i,:)*id2(:))**3
    term12(:) = 4.0*scale_IPP*epsik(i,:)*term6p(:)*term6p(:)
!    do k=1,rs_nbl(rs1)%nnbats
!      evec(1) = evec(1) + term12(k)
!    end do
    evec(1) = evec(1) + sum(term12(:))
    term12(:) = id2(:)*12.0*term12(:)
    foin(:,1) = dvec(:,1)*term12(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*term12(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*term12(:)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + sh(k)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
! 
! specialized routine to compute standard interactions between TIP3P-style
! water molecules
!
subroutine Vforce_PLJ_C_T3P(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use system
!
  implicit none
!
  integer rs1,rs2,ii,j,hi,lo,k,nn,kkk
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%nnbats)
  RTYPE term0(rs_nbl(rs1)%nnbats)
  RTYPE term6(rs_nbl(rs1)%nnbats)
  RTYPE term12(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE id1(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  RTYPE epsik12,epsik6,radik
  integer hira,lora
  integer k1(rs_nbl(rs1)%nnbs,3),k2(rs_nbl(rs1)%nnbs,3)
  RTYPE drav(rs_nbl(rs1)%nnbs,3),srav(rs_nbl(rs1)%nnbs,3)
  RTYPE ca_f(n,3)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
  foin(:,:) = 0.0
!
  ii = refat(rs1)
  drav(:,1) = x(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - x(ii)
  drav(:,2) = y(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - y(ii)
  drav(:,3) = z(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - z(ii)
  hira = rs_nbl(rs1)%nnbs
  lora = 1
!
! PBC
  if ((bnd_type.eq.1).AND.(lora.le.hira)) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(j)-drav(lora:hira,j))/bnd_params(j))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(j))/bnd_params(j))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      srav(lora:hira,1:2) = 0.0
      do j=3,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(6)-drav(lora:hira,j))/bnd_params(6))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(6))/bnd_params(6))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(6)
      end do
    end if
  else 
    srav(:,:) = 0.0
  end if
!
  kkk = b_type(at(rs1)%bb(1))
  epsik12 = 4.0*scale_IPP*lj_eps(bio_ljtyp(kkk),bio_ljtyp(kkk))
  epsik6 = 4.0*scale_attLJ*lj_eps(bio_ljtyp(kkk),bio_ljtyp(kkk))
  radik = lj_sig(bio_ljtyp(kkk),bio_ljtyp(kkk)) 
  nn = rs_nbl(rs1)%nnbs
! we know that all neighbors are separate molecules
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + 2
    dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+2) + srav(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+2) + srav(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+2) + srav(k,3)
    term0(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+2)
    lo = hi + 1
  end do
!
  nn = rs_nbl(rs1)%nnbats
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+2
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    id1(:) = sqrt(id2(:))
!   LJ potential
!   note this is sub-optimal (a third of the terms is zero) but vectorizes cleanly
    if (ii.eq.at(rs1)%bb(1)) then
      term6(1:nn-2:3) = (radik*id2(1:nn-2:3))**3
      term12(1:nn-2:3) = epsik12*term6(1:nn-2:3)*term6(1:nn-2:3)
      term6(1:nn-2:3) = epsik6*term6(1:nn-2:3)
      evec(1) = evec(1) + sum(term12(1:nn-2:3))
      evec(3) = evec(3) - sum(term6(1:nn-2:3))
      term12(1:nn-2:3) = 12.0*term12(1:nn-2:3)
      term6(1:nn-2:3) = 6.0*term6(1:nn-2:3)
      foin(1:nn-2:3,1) = dvec(1:nn-2:3,1)*id2(1:nn-2:3)*&
 &(term12(1:nn-2:3) - term6(1:nn-2:3))
      foin(1:nn-2:3,2) = dvec(1:nn-2:3,2)*id2(1:nn-2:3)*&
 &(term12(1:nn-2:3) - term6(1:nn-2:3))
      foin(1:nn-2:3,3) = dvec(1:nn-2:3,3)*id2(1:nn-2:3)*&
 &(term12(1:nn-2:3) - term6(1:nn-2:3))
    else
      foin(:,:) = 0.0
    end if
!   Coulomb potential
    term1(:) = electric*atq(ii)*term0(:)*scale_POLAR*id1(:)
    evec(6) = evec(6) + sum(term1(:))
    foin(:,1) = foin(:,1) + dvec(:,1)*term1(:)*id2(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = foin(:,2) + dvec(:,2)*term1(:)*id2(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = foin(:,3) + dvec(:,3)*term1(:)*id2(:)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + 2
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
! 
! specialized routine to compute standard interactions between TIP4P-style
! water molecules
!
subroutine Vforce_PLJ_C_T4P(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use system
!
  implicit none
!
  integer rs1,rs2,ii,j,hi,lo,k,kkk,nn
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%nnbats)
  RTYPE term0(rs_nbl(rs1)%nnbats)
  RTYPE term6(rs_nbl(rs1)%nnbats)
  RTYPE term12(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE id1(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  RTYPE epsik12,epsik6,radik
  integer hira,lora
  integer k1(rs_nbl(rs1)%nnbs,3),k2(rs_nbl(rs1)%nnbs,3)
  RTYPE drav(rs_nbl(rs1)%nnbs,3),srav(rs_nbl(rs1)%nnbs,3)
  RTYPE ca_f(n,3)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
  ii = refat(rs1)
  drav(:,1) = x(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - x(ii)
  drav(:,2) = y(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - y(ii)
  drav(:,3) = z(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - z(ii)
  hira = rs_nbl(rs1)%nnbs
  lora = 1
!
! PBC
  if ((bnd_type.eq.1).AND.(lora.le.hira)) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(j)-drav(lora:hira,j))/bnd_params(j))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(j))/bnd_params(j))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      srav(lora:hira,1:2) = 0.0
      do j=3,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(6)-drav(lora:hira,j))/bnd_params(6))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(6))/bnd_params(6))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(6)
      end do
    end if
  else 
    srav(:,:) = 0.0
  end if
!
  kkk = b_type(at(rs1)%bb(1))
  epsik12 = 4.0*scale_IPP*lj_eps(bio_ljtyp(kkk),bio_ljtyp(kkk))
  epsik6 = 4.0*scale_attLJ*lj_eps(bio_ljtyp(kkk),bio_ljtyp(kkk))
  radik = lj_sig(bio_ljtyp(kkk),bio_ljtyp(kkk)) 
  nn = rs_nbl(rs1)%nnbs
!
! first the LJ term
  ii = at(rs1)%bb(1)
  dvec(1:nn,1) = x(refat(rs_nbl(rs1)%nb(1:nn))) - x(ii) + srav(:,1)
  dvec(1:nn,2) = y(refat(rs_nbl(rs1)%nb(1:nn))) - y(ii) + srav(:,2)
  dvec(1:nn,3) = z(refat(rs_nbl(rs1)%nb(1:nn))) - z(ii) + srav(:,3)
  dis2(1:nn) = dvec(1:nn,1)**2 + dvec(1:nn,2)**2 + dvec(1:nn,3)**2
  id2(1:nn) = 1.0/dis2(1:nn)
  term6(1:nn) = (radik*id2(1:nn))**3
  term12(1:nn) = epsik12*term6(1:nn)*term6(1:nn)
  term6(1:nn) = epsik6*term6(1:nn)
  evec(1) = evec(1) + sum(term12(1:nn))
  evec(3) = evec(3) - sum(term6(1:nn))
  term12(1:nn) = 12.0*term12(1:nn)
  term6(1:nn) = 6.0*term6(1:nn)
  foin(1:nn,1) = dvec(1:nn,1)*id2(1:nn)*(term6(1:nn) - term12(1:nn))
  foin(1:nn,2) = dvec(1:nn,2)*id2(1:nn)*(term6(1:nn) - term12(1:nn))
  foin(1:nn,3) = dvec(1:nn,3)*id2(1:nn)*(term6(1:nn) - term12(1:nn))
  ca_f(ii,1) = ca_f(ii,1) + sum(foin(1:nn,1))
  ca_f(ii,2) = ca_f(ii,2) + sum(foin(1:nn,2))
  ca_f(ii,3) = ca_f(ii,3) + sum(foin(1:nn,3))
  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),1) = &
 &  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),1) - foin(1:nn,1)
  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),2) = &
 &  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),2) - foin(1:nn,2)
  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),3) = &
 &  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),3) - foin(1:nn,3)
!
! now the Cb term
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + 2
    dveci(lo:hi,1) = x(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3) + srav(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3) + srav(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3) + srav(k,3)
    term0(lo:hi) = atq(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3)
    lo = hi + 1
  end do
!
  nn = 3*rs_nbl(rs1)%nnbs
  do ii=at(rs1)%bb(1)+1,at(rs1)%bb(1)+3
!   use vector forms as much as possible
!   first distance vectors
    dvec(1:nn,1) = dveci(1:nn,1) - x(ii)
    dvec(1:nn,2) = dveci(1:nn,2) - y(ii)
    dvec(1:nn,3) = dveci(1:nn,3) - z(ii)
    dis2(1:nn) = dvec(1:nn,1)**2 + dvec(1:nn,2)**2 + dvec(1:nn,3)**2
    id2(1:nn) = 1.0/dis2(1:nn)
    id1(1:nn) = sqrt(id2(1:nn))
!   Coulomb potential
    term1(1:nn) = electric*atq(ii)*term0(1:nn)*scale_POLAR*id1(1:nn)
    evec(6) = evec(6) + sum(term1(1:nn))
    foin(1:nn,1) = dvec(1:nn,1)*term1(1:nn)*id2(1:nn)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(1:nn,1))
    for_k(1:nn,1) = for_k(1:nn,1) + foin(1:nn,1)
    foin(1:nn,2) = dvec(1:nn,2)*term1(1:nn)*id2(1:nn)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(1:nn,2))
    for_k(1:nn,2) = for_k(1:nn,2) + foin(1:nn,2)
    foin(1:nn,3) = dvec(1:nn,3)*term1(1:nn)*id2(1:nn)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(1:nn,3))
    for_k(1:nn,3) = for_k(1:nn,3) + foin(1:nn,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + 2
    ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,1) = &
 &     ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,2) = &
 &     ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,3) = &
 &     ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
!
! modified for Ewald electrostatics
!
subroutine Vforce_PEWLJ_C_T3P(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use ewalds
  use sequen
  use system
!
  implicit none
!
  integer rs1,rs2,ii,j,hi,lo,k,nn,kkk
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%nnbats)
  RTYPE term2(rs_nbl(rs1)%nnbats)
  RTYPE term0(rs_nbl(rs1)%nnbats)
  RTYPE term6(rs_nbl(rs1)%nnbats)
  RTYPE term12(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE id1(rs_nbl(rs1)%nnbats)
  RTYPE d1(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  RTYPE epsik12,epsik6,radik
  integer hira,lora
  integer k1(rs_nbl(rs1)%nnbs,3),k2(rs_nbl(rs1)%nnbs,3)
  RTYPE drav(rs_nbl(rs1)%nnbs,3),srav(rs_nbl(rs1)%nnbs,3)
  RTYPE ca_f(n,3)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
  foin(:,:) = 0.0
!
  ii = refat(rs1)
  drav(:,1) = x(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - x(ii)
  drav(:,2) = y(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - y(ii)
  drav(:,3) = z(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - z(ii)
  hira = rs_nbl(rs1)%nnbs
  lora = 1
!
! PBC
  if ((bnd_type.eq.1).AND.(lora.le.hira)) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(j)-drav(lora:hira,j))/bnd_params(j))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(j))/bnd_params(j))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      srav(lora:hira,1:2) = 0.0
      do j=3,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(6)-drav(lora:hira,j))/bnd_params(6))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(6))/bnd_params(6))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(6)
      end do
    end if
  else 
    srav(:,:) = 0.0
  end if
!
  kkk = b_type(at(rs1)%bb(1))
  epsik12 = 4.0*scale_IPP*lj_eps(bio_ljtyp(kkk),bio_ljtyp(kkk))
  epsik6 = 4.0*scale_attLJ*lj_eps(bio_ljtyp(kkk),bio_ljtyp(kkk))
  radik = lj_sig(bio_ljtyp(kkk),bio_ljtyp(kkk)) 
  nn = rs_nbl(rs1)%nnbs
! we know that all neighbors are separate molecules
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + 2
    dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+2) + srav(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+2) + srav(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+2) + srav(k,3)
    term0(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+2)
    lo = hi + 1
  end do
!
  nn = rs_nbl(rs1)%nnbats
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+2
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
!   LJ potential
!   note this is sub-optimal (a third of the terms is zero) but vectorizes cleanly
    if (ii.eq.at(rs1)%bb(1)) then
      term6(1:nn-2:3) = (radik*id2(1:nn-2:3))**3
      term12(1:nn-2:3) = epsik12*term6(1:nn-2:3)*term6(1:nn-2:3)
      term6(1:nn-2:3) = epsik6*term6(1:nn-2:3)
      evec(1) = evec(1) + sum(term12(1:nn-2:3))
      evec(3) = evec(3) - sum(term6(1:nn-2:3))
      term12(1:nn-2:3) = 12.0*term12(1:nn-2:3)
      term6(1:nn-2:3) = 6.0*term6(1:nn-2:3)
      foin(1:nn-2:3,1) = dvec(1:nn-2:3,1)*id2(1:nn-2:3)*&
 &(term12(1:nn-2:3) - term6(1:nn-2:3))
      foin(1:nn-2:3,2) = dvec(1:nn-2:3,2)*id2(1:nn-2:3)*&
 &(term12(1:nn-2:3) - term6(1:nn-2:3))
      foin(1:nn-2:3,3) = dvec(1:nn-2:3,3)*id2(1:nn-2:3)*&
 &(term12(1:nn-2:3) - term6(1:nn-2:3))
    else
      foin(:,:) = 0.0
    end if
    d1(:) = sqrt(dis2(:))
    id1(:) = 1.0/d1(:)
    term1(:) = electric*atq(ii)*term0(:)*scale_POLAR*id1(:)
    term2(:) = term1*(1.0 - erf(ewpm*d1(:)))
    evec(6) = evec(6) + sum(term2(:))
    term2(:) = id2(:)*(term2(:) + &
 &          term1(:)*ewpite*d1(:)*exp(-ewpm2*dis2(:)))
    foin(:,1) = foin(:,1) + dvec(:,1)*term2(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = foin(:,2) + dvec(:,2)*term2(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = foin(:,3) + dvec(:,3)*term2(:)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + 2
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
!
! modified for Ewald electrostatics
! 
subroutine Vforce_PEWLJ_C_T4P(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use ewalds
  use sequen
  use system
!
  implicit none
!
  integer rs1,rs2,ii,j,hi,lo,k,kkk,nn
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%nnbats)
  RTYPE term2(rs_nbl(rs1)%nnbats)
  RTYPE term0(rs_nbl(rs1)%nnbats)
  RTYPE term6(rs_nbl(rs1)%nnbats)
  RTYPE term12(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE id1(rs_nbl(rs1)%nnbats)
  RTYPE d1(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  RTYPE epsik12,epsik6,radik
  integer hira,lora
  integer k1(rs_nbl(rs1)%nnbs,3),k2(rs_nbl(rs1)%nnbs,3)
  RTYPE drav(rs_nbl(rs1)%nnbs,3),srav(rs_nbl(rs1)%nnbs,3)
  RTYPE ca_f(n,3)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
  ii = refat(rs1)
  drav(:,1) = x(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - x(ii)
  drav(:,2) = y(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - y(ii)
  drav(:,3) = z(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - z(ii)
  hira = rs_nbl(rs1)%nnbs
  lora = 1
!
! PBC
  if ((bnd_type.eq.1).AND.(lora.le.hira)) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(j)-drav(lora:hira,j))/bnd_params(j))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(j))/bnd_params(j))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      srav(lora:hira,1:2) = 0.0
      do j=3,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(6)-drav(lora:hira,j))/bnd_params(6))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(6))/bnd_params(6))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(6)
      end do
    end if
  else 
    srav(:,:) = 0.0
  end if
!
  kkk = b_type(at(rs1)%bb(1))
  epsik12 = 4.0*scale_IPP*lj_eps(bio_ljtyp(kkk),bio_ljtyp(kkk))
  epsik6 = 4.0*scale_attLJ*lj_eps(bio_ljtyp(kkk),bio_ljtyp(kkk))
  radik = lj_sig(bio_ljtyp(kkk),bio_ljtyp(kkk)) 
  nn = rs_nbl(rs1)%nnbs
!
! first the LJ term
  ii = at(rs1)%bb(1)
  dvec(1:nn,1) = x(refat(rs_nbl(rs1)%nb(1:nn))) - x(ii) + srav(:,1)
  dvec(1:nn,2) = y(refat(rs_nbl(rs1)%nb(1:nn))) - y(ii) + srav(:,2)
  dvec(1:nn,3) = z(refat(rs_nbl(rs1)%nb(1:nn))) - z(ii) + srav(:,3)
  dis2(1:nn) = dvec(1:nn,1)**2 + dvec(1:nn,2)**2 + dvec(1:nn,3)**2
  id2(1:nn) = 1.0/dis2(1:nn)
  term6(1:nn) = (radik*id2(1:nn))**3
  term12(1:nn) = epsik12*term6(1:nn)*term6(1:nn)
  term6(1:nn) = epsik6*term6(1:nn)
  evec(1) = evec(1) + sum(term12(1:nn))
  evec(3) = evec(3) - sum(term6(1:nn))
  term12(1:nn) = 12.0*term12(1:nn)
  term6(1:nn) = 6.0*term6(1:nn)
  foin(1:nn,1) = dvec(1:nn,1)*id2(1:nn)*(term6(1:nn) - term12(1:nn))
  foin(1:nn,2) = dvec(1:nn,2)*id2(1:nn)*(term6(1:nn) - term12(1:nn))
  foin(1:nn,3) = dvec(1:nn,3)*id2(1:nn)*(term6(1:nn) - term12(1:nn))
  ca_f(ii,1) = ca_f(ii,1) + sum(foin(1:nn,1))
  ca_f(ii,2) = ca_f(ii,2) + sum(foin(1:nn,2))
  ca_f(ii,3) = ca_f(ii,3) + sum(foin(1:nn,3))
  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),1) = &
 &  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),1) - foin(1:nn,1)
  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),2) = &
 &  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),2) - foin(1:nn,2)
  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),3) = &
 &  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),3) - foin(1:nn,3)
!
! now the Cb term
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + 2
    dveci(lo:hi,1) = x(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3) + srav(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3) + srav(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3) + srav(k,3)
    term0(lo:hi) = atq(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3)
    lo = hi + 1
  end do
!
  nn = 3*rs_nbl(rs1)%nnbs
  do ii=at(rs1)%bb(1)+1,at(rs1)%bb(1)+3
!   use vector forms as much as possible
!   first distance vectors
    dvec(1:nn,1) = dveci(1:nn,1) - x(ii)
    dvec(1:nn,2) = dveci(1:nn,2) - y(ii)
    dvec(1:nn,3) = dveci(1:nn,3) - z(ii)
    dis2(1:nn) = dvec(1:nn,1)**2 + dvec(1:nn,2)**2 + dvec(1:nn,3)**2
    d1(1:nn) = sqrt(dis2(1:nn))
    id1(1:nn) = 1.0/d1(1:nn)
    id2(1:nn) = id1(1:nn)*id1(1:nn)
!   Coulomb potential
    term1(1:nn) = electric*atq(ii)*term0(1:nn)*scale_POLAR*id1(1:nn)
    term2(1:nn) = term1(1:nn)*(1.0 - erf(ewpm*d1(1:nn)))
    evec(6) = evec(6) + sum(term2(1:nn))
    term2(1:nn) = id2(1:nn)*(term2(1:nn) + &
 &          term1(1:nn)*ewpite*d1(1:nn)*exp(-ewpm2*dis2(1:nn)))
    foin(1:nn,1) = dvec(1:nn,1)*term2(1:nn)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(1:nn,1))
    for_k(1:nn,1) = for_k(1:nn,1) + foin(1:nn,1)
    foin(1:nn,2) = dvec(1:nn,2)*term2(1:nn)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(1:nn,2))
    for_k(1:nn,2) = for_k(1:nn,2) + foin(1:nn,2)
    foin(1:nn,3) = dvec(1:nn,3)*term2(1:nn)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(1:nn,3))
    for_k(1:nn,3) = for_k(1:nn,3) + foin(1:nn,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + 2
    ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,1) = &
 &     ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,2) = &
 &     ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,3) = &
 &     ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
! 
! the same routine modified for (G)RF electrostatics
!
subroutine Vforce_PRFLJ_C_T3P(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use math
  use system
!
  implicit none
!
  integer rs1,rs2,ii,j,hi,lo,k,nn,kkk
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%nnbats)
  RTYPE term2(rs_nbl(rs1)%nnbats)
  integer isout(rs_nbl(rs1)%nnbats)
  RTYPE term0(rs_nbl(rs1)%nnbats)
  RTYPE term6(rs_nbl(rs1)%nnbats)
  RTYPE term12(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE d2diff(rs_nbl(rs1)%nnbats)
  RTYPE id1(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  RTYPE epsik12,epsik6,radik,ca_f(n,3)
  integer hira,lora
  integer k1(rs_nbl(rs1)%nnbs,3),k2(rs_nbl(rs1)%nnbs,3)
  RTYPE drav(rs_nbl(rs1)%nnbs,3),srav(rs_nbl(rs1)%nnbs,3)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
  foin(:,:) = 0.0
!
  ii = refat(rs1)
  drav(:,1) = x(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - x(ii)
  drav(:,2) = y(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - y(ii)
  drav(:,3) = z(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - z(ii)
  hira = rs_nbl(rs1)%nnbs
  lora = 1
!
! PBC
  if ((bnd_type.eq.1).AND.(lora.le.hira)) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(j)-drav(lora:hira,j))/bnd_params(j))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(j))/bnd_params(j))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      srav(lora:hira,1:2) = 0.0
      do j=3,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(6)-drav(lora:hira,j))/bnd_params(6))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(6))/bnd_params(6))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(6)
      end do
    end if
  else 
    srav(:,:) = 0.0
  end if
!
  kkk = b_type(at(rs1)%bb(1))
  epsik12 = 4.0*scale_IPP*lj_eps(bio_ljtyp(kkk),bio_ljtyp(kkk))
  epsik6 = 4.0*scale_attLJ*lj_eps(bio_ljtyp(kkk),bio_ljtyp(kkk))
  radik = lj_sig(bio_ljtyp(kkk),bio_ljtyp(kkk)) 
  nn = rs_nbl(rs1)%nnbs
! we know that all neighbors are separate molecules
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + 2
    dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+2) + srav(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+2) + srav(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+2) + srav(k,3)
    term0(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+2)
    lo = hi + 1
  end do
!
  nn = rs_nbl(rs1)%nnbats
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+2
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    id1(:) = sqrt(id2(:))
!   LJ potential
!   note this is sub-optimal (a third of the terms is zero) but vectorizes cleanly
    if (ii.eq.at(rs1)%bb(1)) then
      term6(1:nn-2:3) = (radik*id2(1:nn-2:3))**3
      term12(1:nn-2:3) = epsik12*term6(1:nn-2:3)*term6(1:nn-2:3)
      term6(1:nn-2:3) = epsik6*term6(1:nn-2:3)
      evec(1) = evec(1) + sum(term12(1:nn-2:3))
      evec(3) = evec(3) - sum(term6(1:nn-2:3))
      term12(1:nn-2:3) = 12.0*term12(1:nn-2:3)
      term6(1:nn-2:3) = 6.0*term6(1:nn-2:3)
      foin(1:nn-2:3,1) = dvec(1:nn-2:3,1)*id2(1:nn-2:3)*&
 &(term12(1:nn-2:3) - term6(1:nn-2:3))
      foin(1:nn-2:3,2) = dvec(1:nn-2:3,2)*id2(1:nn-2:3)*&
 &(term12(1:nn-2:3) - term6(1:nn-2:3))
      foin(1:nn-2:3,3) = dvec(1:nn-2:3,3)*id2(1:nn-2:3)*&
 &(term12(1:nn-2:3) - term6(1:nn-2:3))
    else
      foin(:,:) = 0.0
    end if
!   Coulomb potential
    term1(:) = electric*atq(ii)*term0(:)*scale_POLAR
!   there is no general minimal cost trick, it entirely depends on how many 
!   atoms in the neighbor list are in fact outside of the cutoff (and hence
!   the cleanest fix would be to solve the issue in NBlist-generation, but
!   that's not recommended for many reasons).
!   here, we're going to zero out term1 for those beyond the cutoff
    d2diff(:) = max((mcel_cutoff2-dis2(:))*imcel2,0.0d0)
    isout(:) = ceiling(d2diff(:))
    term1(:) = term1(:)*dble(isout(:))
    term2(:) = term1(:)*(id1(:)+par_POLAR(1)*dis2(:)-par_POLAR(2))
    evec(6) = evec(6) + sum(term2(:))
    term2(:) = id2(:)*term1(:)*(id1(:)-2.0*par_POLAR(1)*dis2(:))
    foin(:,1) = foin(:,1) + dvec(:,1)*term2(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = foin(:,2) + dvec(:,2)*term2(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = foin(:,3) + dvec(:,3)*term2(:)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + 2
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+2,3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
! 
! the same modified for (G)RF electrostatics
!
subroutine Vforce_PRFLJ_C_T4P(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use math
  use system
!
  implicit none
!
  integer rs1,rs2,ii,j,hi,lo,k,kkk,nn
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%nnbats)
  RTYPE term2(rs_nbl(rs1)%nnbats)
  integer isout(rs_nbl(rs1)%nnbats)
  RTYPE term0(rs_nbl(rs1)%nnbats)
  RTYPE term6(rs_nbl(rs1)%nnbats)
  RTYPE term12(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE d2diff(rs_nbl(rs1)%nnbats)
  RTYPE id1(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  RTYPE epsik12,epsik6,radik,ca_f(n,3)
  integer hira,lora
  integer k1(rs_nbl(rs1)%nnbs,3),k2(rs_nbl(rs1)%nnbs,3)
  RTYPE drav(rs_nbl(rs1)%nnbs,3),srav(rs_nbl(rs1)%nnbs,3)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
  ii = refat(rs1)
  drav(:,1) = x(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - x(ii)
  drav(:,2) = y(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - y(ii)
  drav(:,3) = z(refat(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))) - z(ii)
  hira = rs_nbl(rs1)%nnbs
  lora = 1
!
! PBC
  if ((bnd_type.eq.1).AND.(lora.le.hira)) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(j)-drav(lora:hira,j))/bnd_params(j))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(j))/bnd_params(j))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      srav(lora:hira,1:2) = 0.0
      do j=3,3
        k1(lora:hira,j) = floor((-0.5*bnd_params(6)-drav(lora:hira,j))/bnd_params(6))+1
        k1(lora:hira,j) = max(k1(lora:hira,j),0)
        k2(lora:hira,j) = floor((drav(lora:hira,j)-0.5*bnd_params(6))/bnd_params(6))+1
        k2(lora:hira,j) = max(k2(lora:hira,j),0)
        srav(lora:hira,j) = (k1(lora:hira,j) - k2(lora:hira,j))*bnd_params(6)
      end do
    end if
  else 
    srav(:,:) = 0.0
  end if
!
  kkk = b_type(at(rs1)%bb(1))
  epsik12 = 4.0*scale_IPP*lj_eps(bio_ljtyp(kkk),bio_ljtyp(kkk))
  epsik6 = 4.0*scale_attLJ*lj_eps(bio_ljtyp(kkk),bio_ljtyp(kkk))
  radik = lj_sig(bio_ljtyp(kkk),bio_ljtyp(kkk)) 
  nn = rs_nbl(rs1)%nnbs
!
! first the LJ term
  ii = at(rs1)%bb(1)
  dvec(1:nn,1) = x(refat(rs_nbl(rs1)%nb(1:nn))) - x(ii) + srav(:,1)
  dvec(1:nn,2) = y(refat(rs_nbl(rs1)%nb(1:nn))) - y(ii) + srav(:,2)
  dvec(1:nn,3) = z(refat(rs_nbl(rs1)%nb(1:nn))) - z(ii) + srav(:,3)
  dis2(1:nn) = dvec(1:nn,1)**2 + dvec(1:nn,2)**2 + dvec(1:nn,3)**2
  id2(1:nn) = 1.0/dis2(1:nn)
  term6(1:nn) = (radik*id2(1:nn))**3
  term12(1:nn) = epsik12*term6(1:nn)*term6(1:nn)
  term6(1:nn) = epsik6*term6(1:nn)
  evec(1) = evec(1) + sum(term12(1:nn))
  evec(3) = evec(3) - sum(term6(1:nn))
  term12(1:nn) = 12.0*term12(1:nn)
  term6(1:nn) = 6.0*term6(1:nn)
  foin(1:nn,1) = dvec(1:nn,1)*id2(1:nn)*(term6(1:nn) - term12(1:nn))
  foin(1:nn,2) = dvec(1:nn,2)*id2(1:nn)*(term6(1:nn) - term12(1:nn))
  foin(1:nn,3) = dvec(1:nn,3)*id2(1:nn)*(term6(1:nn) - term12(1:nn))
  ca_f(ii,1) = ca_f(ii,1) + sum(foin(1:nn,1))
  ca_f(ii,2) = ca_f(ii,2) + sum(foin(1:nn,2))
  ca_f(ii,3) = ca_f(ii,3) + sum(foin(1:nn,3))
  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),1) = &
 &  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),1) - foin(1:nn,1)
  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),2) = &
 &  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),2) - foin(1:nn,2)
  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),3) = &
 &  ca_f(refat(rs_nbl(rs1)%nb(1:nn)),3) - foin(1:nn,3)
!
! now the Cb term
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + 2
    dveci(lo:hi,1) = x(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3) + srav(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3) + srav(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3) + srav(k,3)
    term0(lo:hi) = atq(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3)
    lo = hi + 1
  end do
!
  nn = 3*rs_nbl(rs1)%nnbs
  do ii=at(rs1)%bb(1)+1,at(rs1)%bb(1)+3
!   use vector forms as much as possible
!   first distance vectors
    dvec(1:nn,1) = dveci(1:nn,1) - x(ii)
    dvec(1:nn,2) = dveci(1:nn,2) - y(ii)
    dvec(1:nn,3) = dveci(1:nn,3) - z(ii)
    dis2(1:nn) = dvec(1:nn,1)**2 + dvec(1:nn,2)**2 + dvec(1:nn,3)**2
    id2(1:nn) = 1.0/dis2(1:nn)
    id1(1:nn) = sqrt(id2(1:nn))
!   Coulomb potential
    term1(1:nn) = electric*atq(ii)*term0(1:nn)*scale_POLAR
!   there is no general minimal cost trick, it entirely depends on how many 
!   atoms in the neighbor list are in fact outside of the cutoff (and hence
!   the cleanest fix would be to solve the issue in NBlist-generation, but
!   that's not recommended for many reasons).
!   here, we're going to zero out term1 for those beyond the cutoff
    d2diff(1:nn) = max((mcel_cutoff2-dis2(1:nn))*imcel2,0.0d0)
    isout(1:nn) = ceiling(d2diff(1:nn))
    term1(1:nn) = term1(1:nn)*dble(isout(1:nn))
    term2(1:nn) = term1(1:nn)*(id1(1:nn)+par_POLAR(1)*dis2(1:nn)-&
 &                                                   par_POLAR(2))
    evec(6) = evec(6) + sum(term2(1:nn))
    term2(1:nn) = id2(1:nn)*(term1(1:nn)*&
 &                (id1(1:nn)-2.0*par_POLAR(1)*dis2(1:nn)))
    foin(1:nn,1) = dvec(1:nn,1)*term2(1:nn)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(1:nn,1))
    for_k(1:nn,1) = for_k(1:nn,1) + foin(1:nn,1)
    foin(1:nn,2) = dvec(1:nn,2)*term2(1:nn)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(1:nn,2))
    for_k(1:nn,2) = for_k(1:nn,2) + foin(1:nn,2)
    foin(1:nn,3) = dvec(1:nn,3)*term2(1:nn)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(1:nn,3))
    for_k(1:nn,3) = for_k(1:nn,3) + foin(1:nn,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(k)
    hi = lo + 2
    ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,1) = &
 &     ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,2) = &
 &     ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,3) = &
 &     ca_f(at(rs2)%bb(1)+1:at(rs2)%bb(1)+3,3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
!
! this is exactly the same routine as Vforce_PLJ_C, only it uses the
! mid-range neighbor-list (between MCNBCUTOFF and MCELCUTOFF)
! note the use of rs_nbl(rs1)%trsvec for image shift vectors
!
subroutine Vforce_PLJ_TR(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
!
  implicit none
!
  integer rs1,rs2,ii,hi,lo,k,i,shii
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%nnbtrats)
  RTYPE term0(rs_nbl(rs1)%nnbtrats)
  RTYPE term6(rs_nbl(rs1)%nnbtrats)
  RTYPE term6p(rs_nbl(rs1)%nnbtrats)
  RTYPE term12(rs_nbl(rs1)%nnbtrats)
  RTYPE foin(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbtrats,3)
  RTYPE id2(rs_nbl(rs1)%nnbtrats)
  RTYPE dis2(rs_nbl(rs1)%nnbtrats)
  RTYPE id1(rs_nbl(rs1)%nnbtrats)
  RTYPE for_k(rs_nbl(rs1)%nnbtrats,3)
  RTYPE epsik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbtrats)
  RTYPE radik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbtrats)
  integer sh(rs_nbl(rs1)%nnbtrs)
  RTYPE ca_f(n,3)
!  
!  hi = n!at(rs2)%bb(1)+n-at(rs1+1)%bb(1)+1-1
!  lo = at(rs1+1)%bb(1)
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nbb + at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nsc - 1
  shii = at(rs1)%nbb + at(rs1)%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(k)
    hi = lo + sh(k)
    dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,3)
    term0(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    do ii=lo,hi
      epsik(:,ii) = lj_eps(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
      radik(:,ii) = lj_sig(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
    end do
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    id1(:) = sqrt(id2(:))
    term6p(:) = (radik(i,:)*id2(:))**3
    term12(:) = 4.0*scale_IPP*epsik(i,:)*term6p(:)*term6p(:)
    term6(:) = scale_attLJ*4.0*epsik(i,:)*term6p(:)
    evec(1) = evec(1) + sum(term12(:))
    evec(3) = evec(3) - sum(term6(:))
    term12(:) = 12.0*term12(:)*id2(:)
    term6(:) = 6.0*term6(:)*id2(:)
!   Coulomb potential
    term1(:) = electric*atq(ii)*term0(:)*scale_POLAR*id1(:)
    evec(6) = evec(6) + sum(term1(:))
    foin(:,1) = dvec(:,1)*(term1*id2(:) - term6(:) + term12(:))
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*(term1*id2(:) - term6(:) + term12(:))
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*(term1*id2(:) - term6(:) + term12(:))
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(k)
    hi = lo + sh(k)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1)+for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2)+for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3)+for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
!
! the same routine modified for Ewald electrostatics (never used)
! notice the use of rs_nbl(rs1)%trsvec for image shift vectors
!
subroutine Vforce_PEWLJ_TR(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use ewalds
  use units
  use cutoffs
!
  implicit none
!
  integer rs1,rs2,ii,hi,lo,k,i,shii
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%nnbtrats)
  RTYPE term2(rs_nbl(rs1)%nnbtrats)
  RTYPE term0(rs_nbl(rs1)%nnbtrats)
  RTYPE term6(rs_nbl(rs1)%nnbtrats)
  RTYPE term6p(rs_nbl(rs1)%nnbtrats)
  RTYPE term12(rs_nbl(rs1)%nnbtrats)
  RTYPE foin(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbtrats,3)
  RTYPE id2(rs_nbl(rs1)%nnbtrats)
  RTYPE dis2(rs_nbl(rs1)%nnbtrats)
  RTYPE id1(rs_nbl(rs1)%nnbtrats)
  RTYPE d1(rs_nbl(rs1)%nnbtrats)
  RTYPE for_k(rs_nbl(rs1)%nnbtrats,3)
  RTYPE epsik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbtrats)
  RTYPE radik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbtrats)
  integer sh(rs_nbl(rs1)%nnbtrs)
  RTYPE ca_f(n,3)
!  
!  hi = n!at(rs2)%bb(1)+n-at(rs1+1)%bb(1)+1-1
!  lo = at(rs1+1)%bb(1)
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nbb + at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nsc - 1
  shii = at(rs1)%nbb + at(rs1)%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(k)
    hi = lo + sh(k)
    dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,3)
    term0(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    do ii=lo,hi
      epsik(:,ii) = lj_eps(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
      radik(:,ii) = lj_sig(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
    end do
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    term6p(:) = (radik(i,:)*id2(:))**3
    term12(:) = 4.0*scale_IPP*epsik(i,:)*term6p(:)*term6p(:)
    term6(:) = scale_attLJ*4.0*epsik(i,:)*term6p(:)
    evec(1) = evec(1) + sum(term12(:))
    evec(3) = evec(3) - sum(term6(:))
    term12(:) = 12.0*term12(:)
    term6(:) = 6.0*term6(:)
!   Coulomb potential
    d1(:) = sqrt(dis2(:))
    id1(:) = 1.0/d1(:)
    term1(:) = electric*atq(ii)*term0(:)*scale_POLAR*id1(:)
    term2(:) = term1*(1.0 - erf(ewpm*d1(:)))
    evec(6) = evec(6) + sum(term2(:))
    term2(:) = id2(:)*(term2(:) - term6(:) + term12(:) + &
 &  term1(:)*ewpite*d1(:)*exp(-ewpm2*dis2(:)))
    foin(:,1) = dvec(:,1)*term2(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*term2(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*term2(:)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(k)
    hi = lo + sh(k)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1)+for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2)+for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3)+for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
!
! the same routine modified for RF electrostatics
! note the complication here due to particles outside(!) of the second cutoff
! also note the use of rs_nbl(rs1)%trsvec for image shift vectors
!
subroutine Vforce_PRFLJ_TR(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use math
!
  implicit none
!
  integer rs1,rs2,ii,hi,lo,k,i,shii
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%nnbtrats)
  integer isout(rs_nbl(rs1)%nnbtrats)
  RTYPE term2(rs_nbl(rs1)%nnbtrats)
  RTYPE term0(rs_nbl(rs1)%nnbtrats)
  RTYPE term6(rs_nbl(rs1)%nnbtrats)
  RTYPE term6p(rs_nbl(rs1)%nnbtrats)
  RTYPE term12(rs_nbl(rs1)%nnbtrats)
  RTYPE foin(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbtrats,3)
  RTYPE id2(rs_nbl(rs1)%nnbtrats)
  RTYPE dis2(rs_nbl(rs1)%nnbtrats)
  RTYPE d2diff(rs_nbl(rs1)%nnbtrats)
  RTYPE id1(rs_nbl(rs1)%nnbtrats)
  RTYPE for_k(rs_nbl(rs1)%nnbtrats,3)
  RTYPE epsik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbtrats)
  RTYPE radik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbtrats)
  RTYPE ca_f(n,3)
  integer sh(rs_nbl(rs1)%nnbtrs)
!  hi = n!at(rs2)%bb(1)+n-at(rs1+1)%bb(1)+1-1
!  lo = at(rs1+1)%bb(1)
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nbb + at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nsc - 1
  shii = at(rs1)%nbb + at(rs1)%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(k)
    hi = lo + sh(k)
    dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,3)
    term0(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    do ii=lo,hi
      epsik(:,ii) = lj_eps(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
      radik(:,ii) = lj_sig(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
    end do
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    id1(:) = sqrt(id2(:))
    term6p(:) = (radik(i,:)*id2(:))**3
    term12(:) = 4.0*scale_IPP*epsik(i,:)*term6p(:)*term6p(:)
    term6(:) = scale_attLJ*4.0*epsik(i,:)*term6p(:)
    evec(1) = evec(1) + sum(term12(:))
    evec(3) = evec(3) - sum(term6(:))
    term12(:) = 12.0*term12(:)
    term6(:) = 6.0*term6(:)
!   Coulomb potential
    term1(:) = electric*atq(ii)*term0(:)*scale_POLAR
!   there is no general minimal cost trick, it entirely depends on how many 
!   atoms in the neighbor list are in fact outside of the cutoff (and hence
!   the cleanest fix would be to solve the issue in NBlist-generation, but
!   that's not recommended for many reasons).
!   here, we're going to zero out term1 for those beyond the cutoff
    d2diff(:) = max((mcel_cutoff2-dis2(:))*imcel2,0.0d0)
    isout(:) = ceiling(d2diff(:))
    term1(:) = term1(:)*dble(isout(:))
    term2(:) = term1(:)*(id1(:)+par_POLAR(1)*dis2(:)-par_POLAR(2))
    evec(6) = evec(6) + sum(term2(:))
    term2(:) = id2(:)*(term1(:)*(id1(:)-2.0*par_POLAR(1)*dis2(:))&
 &                                   - term6(:) + term12(:))
    foin(:,1) = dvec(:,1)*term2(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*term2(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*term2(:)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(k)
    hi = lo + sh(k)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1)+for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2)+for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3)+for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
!
! the same routine for ghosted interactions
! notice the use of rs_nbl(rs1)%gtrsvec for image shift vectors
!
subroutine Vforce_FEGPLJ_TR(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
!
  implicit none
!
  integer rs1,rs2,ii,hi,lo,k,i,shii
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%ngnbtrats)
  RTYPE term3(rs_nbl(rs1)%ngnbtrats)
  RTYPE term0(rs_nbl(rs1)%ngnbtrats)
  RTYPE term6(rs_nbl(rs1)%ngnbtrats)
  RTYPE term6p(rs_nbl(rs1)%ngnbtrats)
  RTYPE term12(rs_nbl(rs1)%ngnbtrats)
  RTYPE term12p(rs_nbl(rs1)%ngnbtrats)
  RTYPE foin(rs_nbl(rs1)%ngnbtrats,3)
  RTYPE dveci(rs_nbl(rs1)%ngnbtrats,3)
  RTYPE dvec(rs_nbl(rs1)%ngnbtrats,3)
  RTYPE id2(rs_nbl(rs1)%ngnbtrats)
  RTYPE d1(rs_nbl(rs1)%ngnbtrats)
  RTYPE dis2(rs_nbl(rs1)%ngnbtrats)
  RTYPE id1(rs_nbl(rs1)%ngnbtrats)
  RTYPE for_k(rs_nbl(rs1)%ngnbtrats,3)
  RTYPE epsik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%ngnbtrats)
  RTYPE radik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%ngnbtrats)
  RTYPE ca_f(n,3)
  integer sh(rs_nbl(rs1)%ngnbtrs)
! 
!  hi = n!at(rs2)%bb(1)+n-at(rs1+1)%bb(1)+1-1
!  lo = at(rs1+1)%bb(1)
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%gnbtr(1:rs_nbl(rs1)%ngnbtrs))%nbb + at(rs_nbl(rs1)%gnbtr(1:rs_nbl(rs1)%ngnbtrs))%nsc-1
  shii = at(rs1)%nbb + at(rs1)%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%ngnbtrs
    rs2 = rs_nbl(rs1)%gnbtr(k)
    hi = lo + sh(k)
    dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%gtrsvec(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%gtrsvec(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%gtrsvec(k,3)
    term0(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    do ii=lo,hi
      epsik(:,ii) = lj_eps(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
      radik(:,ii) = lj_sig(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
    end do
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
!
!   ghosted LJ
    term3 = (dis2(:)/radik(i,:))**3
    term12p(:) = 1.0/(term3(:) + par_FEG2(2))
    term12(:) = 4.0*par_FEG2(1)*epsik(i,:)*term12p(:)*term12p(:)
    term6p(:) = 1.0/(term3(:) + par_FEG2(6))
    term6(:) = par_FEG2(5)*4.0*epsik(i,:)*term6p(:)
    evec(1) = evec(1) + sum(term12(:))
    evec(3) = evec(3) - sum(term6(:))
    term12(:) = 12.0*term12(:)*term12p(:)*id2(:)*term3(:)
    term6(:) = 6.0*term6(:)*term6p(:)*id2(:)*term3(:)
!
!   ghosted Coulomb potential
    d1(:) = sqrt(dis2(:))
    id1(:) = 1.0/(d1(:) + par_FEG2(10))
    term1(:) = electric*atq(ii)*term0(:)*par_FEG2(9)*id1(:)
    evec(6) = evec(6) + sum(term1)
    term1(:) = term1(:)*id1(:)/d1(:)
!
!   Force increment
    foin(:,1) = dvec(:,1)*(term1 - term6(:) + term12(:))
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*(term1 - term6(:) + term12(:))
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*(term1 - term6(:) + term12(:))
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%ngnbtrs
    rs2 = rs_nbl(rs1)%gnbtr(k)
    hi = lo + sh(k)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1)+for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2)+for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3)+for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
!
! the same routine for the combination of ghosted interactions and (G)RF electrostatics
! assuming simple linear scaling for the Cb-ghosting for the TR-regime
! notice the use of rs_nbl(rs1)%gtrsvec for image shift vectors
!
subroutine Vforce_FEGPRFLJ_TR(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use math
  use cutoffs
!
  implicit none
!
  integer rs1,rs2,ii,hi,lo,k,i,shii
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term1(rs_nbl(rs1)%ngnbtrats)
  RTYPE term2(rs_nbl(rs1)%ngnbtrats)
  integer isout(rs_nbl(rs1)%ngnbtrats)
  RTYPE term3(rs_nbl(rs1)%ngnbtrats)
  RTYPE term0(rs_nbl(rs1)%ngnbtrats)
  RTYPE term6(rs_nbl(rs1)%ngnbtrats)
  RTYPE term6p(rs_nbl(rs1)%ngnbtrats)
  RTYPE term12(rs_nbl(rs1)%ngnbtrats)
  RTYPE term12p(rs_nbl(rs1)%ngnbtrats)
  RTYPE d1(rs_nbl(rs1)%ngnbtrats)
  RTYPE foin(rs_nbl(rs1)%ngnbtrats,3)
  RTYPE dveci(rs_nbl(rs1)%ngnbtrats,3)
  RTYPE dvec(rs_nbl(rs1)%ngnbtrats,3)
  RTYPE id2(rs_nbl(rs1)%ngnbtrats)
  RTYPE dis2(rs_nbl(rs1)%ngnbtrats)
  RTYPE d2diff(rs_nbl(rs1)%ngnbtrats)
  RTYPE id1(rs_nbl(rs1)%ngnbtrats)
  RTYPE for_k(rs_nbl(rs1)%ngnbtrats,3)
  RTYPE epsik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%ngnbtrats)
  RTYPE radik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%ngnbtrats)
  integer sh(rs_nbl(rs1)%ngnbtrs)
  RTYPE ca_f(n,3)
!  
!  hi = n!at(rs2)%bb(1)+n-at(rs1+1)%bb(1)+1-1
!  lo = at(rs1+1)%bb(1)
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%gnbtr(1:rs_nbl(rs1)%ngnbtrs))%nbb+at(rs_nbl(rs1)%gnbtr(1:rs_nbl(rs1)%ngnbtrs))%nsc-1
  shii = at(rs1)%nbb + at(rs1)%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%ngnbtrs
    rs2 = rs_nbl(rs1)%gnbtr(k)
    hi = lo + sh(k)
    dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%gtrsvec(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%gtrsvec(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%gtrsvec(k,3)
    term0(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    do ii=lo,hi
      epsik(:,ii) = lj_eps(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
      radik(:,ii) = lj_sig(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
    end do
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
!
!   ghosted LJ
    term3 = (dis2(:)/radik(i,:))**3
    term12p(:) = 1.0/(term3(:) + par_FEG2(2))
    term12(:) = 4.0*par_FEG2(1)*epsik(i,:)*term12p(:)*term12p(:)
    term6p(:) = 1.0/(term3(:) + par_FEG2(6))
    term6(:) = par_FEG2(5)*4.0*epsik(i,:)*term6p(:)
    evec(1) = evec(1) + sum(term12(:))
    evec(3) = evec(3) - sum(term6(:))
    term12(:) = 12.0*term12(:)*term12p(:)*term3(:)
    term6(:) = 6.0*term6(:)*term6p(:)*term3(:)
!
!   ghosted Coulomb with RF potential: this only works correctly for fegcbmode = 1
    d1(:) = sqrt(dis2(:))
    id1(:) = 1.0/(d1(:))
    term1(:) = electric*atq(ii)*term0(:)*par_FEG2(9)
!   there is no general minimal cost trick, it entirely depends on how many 
!   atoms in the neighbor list are in fact outside of the cutoff (and hence
!   the cleanest fix would be to solve the issue in NBlist-generation, but
!   that's not recommended for many reasons).
!   here, we're going to zero out term1 for those beyond the cutoff
    d2diff(:) = max((mcel_cutoff2-dis2(:))*imcel2,0.0d0)
    isout(:) = ceiling(d2diff(:))
    term1(:) = term1(:)*dble(isout(:))
    term2(:) = term1(:)*(id1(:)+par_POLAR(1)*dis2(:)-par_POLAR(2))
    evec(6) = evec(6) + sum(term2)
    term2(:) = id2(:)*(term1(:)*(id1(:)-2.0*par_POLAR(1)*dis2(:))&
 &                                   - term6(:) + term12(:))
    foin(:,1) = dvec(:,1)*term2(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*term2(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*term2(:)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%ngnbtrs
    rs2 = rs_nbl(rs1)%gnbtr(k)
    hi = lo + sh(k)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end

!
!-----------------------------------------------------------------------
!
! the same routine for just the LJ potential
! notice the use of rs_nbl(rs1)%trsvec for image shift vectors
!
subroutine Vforce_LJ_TR(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use system
  use sequen
!
  implicit none
!
  integer rs1,rs2,ii,hi,lo,k,i,shii
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term6(rs_nbl(rs1)%nnbtrats)
  RTYPE term6p(rs_nbl(rs1)%nnbtrats)
  RTYPE term12(rs_nbl(rs1)%nnbtrats)
  RTYPE foin(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbtrats,3)
  RTYPE id2(rs_nbl(rs1)%nnbtrats)
  RTYPE dis2(rs_nbl(rs1)%nnbtrats)
  RTYPE for_k(rs_nbl(rs1)%nnbtrats,3)
  RTYPE epsik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbtrats)
  RTYPE radik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbtrats)
  integer sh(rs_nbl(rs1)%nnbtrs)
  RTYPE ca_f(n,3)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh(:) = at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nbb + at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nsc - 1
  shii = at(rs1)%nbb + at(rs1)%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(k)
    hi = lo + sh(k)
    dveci(lo:hi,1) = x(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k))) + rs_nbl(rs1)%trsvec(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k))) + rs_nbl(rs1)%trsvec(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k))) + rs_nbl(rs1)%trsvec(k,3)
    do ii=lo,hi
      epsik(:,ii) = lj_eps(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
      radik(:,ii) = lj_sig(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
    end do
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    term6p(:) = (radik(i,:)*id2(:))**3
    term12(:) = 4.0*scale_IPP*epsik(i,:)*term6p(:)*term6p(:)
    term6(:) = scale_attLJ*4.0*epsik(i,:)*term6p(:)
    evec(1) = evec(1) + sum(term12(:))
    evec(3) = evec(3) - sum(term6(:))
    term6(:) = id2(:)*(12.0*term12(:) - 6.0*term6(:))
    foin(:,1) = dvec(:,1)*term6(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*term6(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*term6(:)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
!   virial
!    if (pflag.EQV..true.) then
!      ens%insVirT(1,1) = ens%insVirT(1,1) + sum(foin(:,1)*dvec(:,1))
!      ens%insVirT(2,2) = ens%insVirT(2,2) + sum(foin(:,2)*dvec(:,2))
!      ens%insVirT(3,3) = ens%insVirT(3,3) + sum(foin(:,3)*dvec(:,3))
!    end if
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(k)
    hi = lo + sh(k)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!
!-----------------------------------------------------------------------
!
! the same routine for just excluded volume interactions
! notice the use of rs_nbl(rs1)%trsvec for image shift vectors
! 
subroutine Vforce_IPP_TR(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use system
  use molecule
!
  implicit none
!
  integer rs1,rs2,ii,hi,lo,k,i,shii
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term6p(rs_nbl(rs1)%nnbtrats)
  RTYPE term12(rs_nbl(rs1)%nnbtrats)
  RTYPE foin(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbtrats,3)
  RTYPE id2(rs_nbl(rs1)%nnbtrats)
  RTYPE dis2(rs_nbl(rs1)%nnbtrats)
  RTYPE for_k(rs_nbl(rs1)%nnbtrats,3)
  RTYPE epsik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbtrats)
  RTYPE radik(at(rs1)%nbb+at(rs1)%nsc,rs_nbl(rs1)%nnbtrats)
  integer sh(rs_nbl(rs1)%nnbtrs)
  RTYPE ca_f(n,3)
! 
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nbb + at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nsc - 1
  shii = at(rs1)%nbb + at(rs1)%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(k)
    hi = lo + sh(k)
    dveci(lo:hi,1) = x(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k))) + rs_nbl(rs1)%trsvec(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k))) + rs_nbl(rs1)%trsvec(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k))) + rs_nbl(rs1)%trsvec(k,3)
    do ii=lo,hi
      epsik(:,ii) = lj_eps(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
      radik(:,ii) = lj_sig(attyp(at(rs1)%bb(1):(at(rs1)%bb(1)+shii)),attyp(at(rs2)%bb(1)+ii-lo))
    end do
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    term6p(:) = (radik(i,:)*id2(:))**3
    term12(:) = 4.0*scale_IPP*epsik(i,:)*term6p(:)*term6p(:)
    evec(1) = evec(1) + sum(term12(:))
    term12(:) = id2(:)*12.0*term12(:)
    foin(:,1) = dvec(:,1)*term12(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*term12(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*term12(:)
    ca_f(ii,3) = ca_f(ii,3) - sum(foin(:,3))
    for_k(:,3) = for_k(:,3) + foin(:,3)
  end do
  lo = 1
  do k=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(k)
    hi = lo + sh(k)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1) + for_k(lo:hi,1)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2) + for_k(lo:hi,2)
    ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) = &
 &     ca_f(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3) + for_k(lo:hi,3)
    lo = hi + 1
  end do
!
end
!
!-------------------------------------------------------------------------------
!
! a specialized routine suitable for sparse tabulated potentials
! note that this is not crosslink-adjusted  and does not require 
! calls to the correction fxn in the wrapper construct (force_wrap)
!
subroutine Vforce_TAB(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use molecule
  use sequen
  use params
  use tabpot
  use units
  use math
  use cutoffs
!
  implicit none
!
  
  integer rs1,rs2,ii,kk,hi,lo,alcsz,k,i,jj,k2
  RTYPE sv(3),pfac,pfac2,hlp1,hlp2,hlp3,hlp4,dhlp1,dhlp2,dhlp3,evec(MAXENERGYTERMS)
  RTYPE dvec(rs_nbl(rs1)%ntabias,3)
  RTYPE d2(rs_nbl(rs1)%ntabias)
  RTYPE term0(rs_nbl(rs1)%ntabias),term1(rs_nbl(rs1)%ntabias),term2(rs_nbl(rs1)%ntabias),term3(rs_nbl(rs1)%ntabias)
  RTYPE id1(rs_nbl(rs1)%ntabias)
  RTYPE d1(rs_nbl(rs1)%ntabias),ca_f(n,3)
  integer tbin(rs_nbl(rs1)%ntabias)
!     
  lo = 1 
  do jj=1,rs_nbl(rs1)%ntabnbs
    rs2 = rs_nbl(rs1)%tabnb(jj)
    hi = lo + (tbp%rsmat(rs2,rs1)-tbp%rsmat(rs1,rs2))
    call dis_bound_rs2(rs1,rs2,sv)
    k = 0
    do i=tbp%rsmat(rs1,rs2),tbp%rsmat(rs2,rs1)
      ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
      kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
      dvec(lo+k,1) = x(kk) - x(ii)
      dvec(lo+k,2) = y(kk) - y(ii)
      dvec(lo+k,3) = z(kk) - z(ii)
      k = k + 1
    end do
    dvec(lo:hi,1) =  dvec(lo:hi,1) + sv(1)
    dvec(lo:hi,2) =  dvec(lo:hi,2) + sv(2)
    dvec(lo:hi,3) =  dvec(lo:hi,3) + sv(3)
    lo = hi + 1
  end do
!
  alcsz = hi
  d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
  d1(1:alcsz) = sqrt(d2(1:alcsz))
  pfac = 1.0/tbp%res
  pfac2 = pfac*scale_TABUL
  term0(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
  tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(term0(1:alcsz))))
  term1(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
  term2(1:alcsz) = term1(1:alcsz)*term1(1:alcsz) ! quadratic
  term3(1:alcsz) = term2(1:alcsz)*term1(1:alcsz) ! cubic
  id1(1:alcsz) = 1.0/d1(1:alcsz)
!  tbin(1:alcsz) = floor((d1(1:alcsz)-tbp%dis(1))/tbp%res) + 1
!  pfac = scale_TABUL*(1.0/tbp%res)
!
  lo = 1
  do jj=1,rs_nbl(rs1)%ntabnbs
    rs2 = rs_nbl(rs1)%tabnb(jj)
    if ((rs2-rs1).eq.1) then
!     essentially cycle since tabnb still contains sequence neighbors
      hi = lo + (tbp%rsmat(rs2,rs1)-tbp%rsmat(rs1,rs2))
      lo = hi + 1
      cycle
    end if
    hi = lo + (tbp%rsmat(rs2,rs1)-tbp%rsmat(rs1,rs2))
    k = 0
    do i=tbp%rsmat(rs1,rs2),tbp%rsmat(rs2,rs1)
      if (tbin(lo+k).ge.tbp%bins) then
        evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
      else if (tbin(lo+k).ge.1) then
        k2 = lo+k
        ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
        kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
        hlp1 = 2.0*term3(k2) - 3.0*term2(k2)
        hlp2 = term3(k2) - term2(k2)
        hlp3 = hlp2 - term2(k2) + term1(k2)
        dhlp1 = 6.0*term2(k2) - 6.0*term1(k2)
        dhlp2 = 3.0*term2(k2) - 2.0*term1(k2)
        dhlp3 = dhlp2 - 2.0*term1(k2) + 1.0
        hlp4 = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k2)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k2)+1) + &
 &               tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k2)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k2)+1))
!        term1 = (tbp%dis(tbin(k2)+1) - d1(k2))*tbp%pot(tbp%lst(i,3),tbin(k2))
!        term2 = (d1(k2) - tbp%dis(tbin(k2)))*tbp%pot(tbp%lst(i,3),tbin(k2)+1)
!        term3 = pfac*(term1+term2)
        evec(9) = evec(9) + scale_TABUL*hlp4
        hlp4 = pfac2*id1(k2)*(dhlp1*(tbp%pot(tbp%lst(i,3),tbin(k2)) - tbp%pot(tbp%lst(i,3),tbin(k2)+1)) + &
 &               tbp%res*(dhlp3*tbp%tang(tbp%lst(i,3),tbin(k2)) + dhlp2*tbp%tang(tbp%lst(i,3),tbin(k2)+1)))
!        term4 = id1(k2)*pfac*(tbp%pot(tbp%lst(i,3),tbin(k2)+1)-tbp%pot(tbp%lst(i,3),tbin(k2)))
        ca_f(ii,1:3) = ca_f(ii,1:3) + hlp4*dvec(k2,1:3)
        ca_f(kk,1:3) = ca_f(kk,1:3) - hlp4*dvec(k2,1:3)
      else
        evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),1)
      end if
      k = k + 1
    end do
    lo = hi + 1
  end do
!
  end
!
!-------------------------------------------------------------------------------
!
! the same things for standard NB-list
!
subroutine Vforce_TAB_C(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use molecule
  use sequen
  use params
  use tabpot
  use units
  use math
  use cutoffs
!
  implicit none
!
  
  integer rs1,rs2,ii,kk,hi,lo,alcsz,k,i,jj,k2
  RTYPE sv(3),pfac,pfac2,hlp1,hlp2,hlp3,hlp4,dhlp1,dhlp2,dhlp3,evec(MAXENERGYTERMS)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE d2(rs_nbl(rs1)%nnbats)
  RTYPE term0(rs_nbl(rs1)%ntabias),term1(rs_nbl(rs1)%ntabias),term2(rs_nbl(rs1)%ntabias),term3(rs_nbl(rs1)%ntabias)
  RTYPE id1(rs_nbl(rs1)%nnbats)
  RTYPE d1(rs_nbl(rs1)%nnbats),ca_f(n,3)
  integer tbin(rs_nbl(rs1)%nnbats)
!     
  lo = 1 
  do jj=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(jj)
    hi = lo + (tbp%rsmat(rs2,rs1)-tbp%rsmat(rs1,rs2))
    call dis_bound_rs2(rs1,rs2,sv)
    k = 0
    do i=tbp%rsmat(rs1,rs2),tbp%rsmat(rs2,rs1)
      ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
      kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
      dvec(lo+k,1) = x(kk) - x(ii)
      dvec(lo+k,2) = y(kk) - y(ii)
      dvec(lo+k,3) = z(kk) - z(ii)
      k = k + 1
    end do
    dvec(lo:hi,1) =  dvec(lo:hi,1) + sv(1)
    dvec(lo:hi,2) =  dvec(lo:hi,2) + sv(2)
    dvec(lo:hi,3) =  dvec(lo:hi,3) + sv(3)
    lo = hi + 1
  end do
!
  alcsz = hi
  d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
  d1(1:alcsz) = sqrt(d2(1:alcsz))
  pfac = 1.0/tbp%res
  pfac2 = pfac*scale_TABUL
  term0(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
  tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(term0(1:alcsz))))
  term1(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
  term2(1:alcsz) = term1(1:alcsz)*term1(1:alcsz) ! quadratic
  term3(1:alcsz) = term2(1:alcsz)*term1(1:alcsz) ! cubic
  id1(1:alcsz) = 1.0/d1(1:alcsz)
!
  lo = 1
  do jj=1,rs_nbl(rs1)%nnbs
    rs2 = rs_nbl(rs1)%nb(jj)
    hi = lo + (tbp%rsmat(rs2,rs1)-tbp%rsmat(rs1,rs2))
    k = 0
    do i=tbp%rsmat(rs1,rs2),tbp%rsmat(rs2,rs1)
      if (tbin(lo+k).ge.tbp%bins) then
        evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
      else if (tbin(lo+k).ge.1) then
        k2 = lo+k
        ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
        kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
        hlp1 = 2.0*term3(k2) - 3.0*term2(k2)
        hlp2 = term3(k2) - term2(k2)
        hlp3 = hlp2 - term2(k2) + term1(k2)
        dhlp1 = 6.0*term2(k2) - 6.0*term1(k2)
        dhlp2 = 3.0*term2(k2) - 2.0*term1(k2)
        dhlp3 = dhlp2 - 2.0*term1(k2) + 1.0
        hlp4 = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k2)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k2)+1) + &
 &               tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k2)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k2)+1))
        evec(9) = evec(9) + scale_TABUL*hlp4
        hlp4 = pfac2*id1(k2)*(dhlp1*(tbp%pot(tbp%lst(i,3),tbin(k2)) - tbp%pot(tbp%lst(i,3),tbin(k2)+1)) + &
 &               tbp%res*(dhlp3*tbp%tang(tbp%lst(i,3),tbin(k2)) + dhlp2*tbp%tang(tbp%lst(i,3),tbin(k2)+1)))
        ca_f(ii,1:3) = ca_f(ii,1:3) + hlp4*dvec(k2,1:3)
        ca_f(kk,1:3) = ca_f(kk,1:3) - hlp4*dvec(k2,1:3)
      else
        evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),1)
      end if
      k = k + 1
    end do
    lo = hi + 1
  end do
!
  end
!
!-------------------------------------------------------------------------------
!
! the same thing for the twin-range regime
!
subroutine Vforce_TAB_TR(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use molecule
  use sequen
  use params
  use tabpot
  use units
  use math
  use cutoffs
!
  implicit none
!
  integer rs1,rs2,ii,kk,hi,lo,alcsz,k,i,jj,k2
  RTYPE sv(3),pfac,pfac2,hlp1,hlp2,hlp3,hlp4,dhlp1,dhlp2,dhlp3,evec(MAXENERGYTERMS)
  RTYPE dvec(rs_nbl(rs1)%nnbtrats,3)
  RTYPE d2(rs_nbl(rs1)%nnbtrats)
  RTYPE term0(rs_nbl(rs1)%ntabias),term1(rs_nbl(rs1)%ntabias),term2(rs_nbl(rs1)%ntabias),term3(rs_nbl(rs1)%ntabias)
  RTYPE id1(rs_nbl(rs1)%nnbtrats)
  RTYPE d1(rs_nbl(rs1)%nnbtrats),ca_f(n,3)
  integer tbin(rs_nbl(rs1)%nnbtrats)
!   
  lo = 1 
  do jj=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(jj)
    hi = lo + (tbp%rsmat(rs2,rs1)-tbp%rsmat(rs1,rs2))
    call dis_bound_rs2(rs1,rs2,sv)
    k = 0
    do i=tbp%rsmat(rs1,rs2),tbp%rsmat(rs2,rs1)
      ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
      kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
      dvec(lo+k,1) = x(kk) - x(ii)
      dvec(lo+k,2) = y(kk) - y(ii)
      dvec(lo+k,3) = z(kk) - z(ii)
      k = k + 1
    end do
    dvec(lo:hi,1) =  dvec(lo:hi,1) + sv(1)
    dvec(lo:hi,2) =  dvec(lo:hi,2) + sv(2)
    dvec(lo:hi,3) =  dvec(lo:hi,3) + sv(3)
    lo = hi + 1
  end do
!
  alcsz = hi
  d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
  d1(1:alcsz) = sqrt(d2(1:alcsz))
  pfac = 1.0/tbp%res
  pfac2 = pfac*scale_TABUL
  term0(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
  tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(term0(1:alcsz))))
  term1(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
  term2(1:alcsz) = term1(1:alcsz)*term1(1:alcsz) ! quadratic
  term3(1:alcsz) = term2(1:alcsz)*term1(1:alcsz) ! cubic
  id1(1:alcsz) = 1.0/d1(1:alcsz)
!
  lo = 1
  do jj=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(jj)
    hi = lo + (tbp%rsmat(rs2,rs1)-tbp%rsmat(rs1,rs2))
    k = 0
    do i=tbp%rsmat(rs1,rs2),tbp%rsmat(rs2,rs1)
      if (tbin(lo+k).ge.tbp%bins) then
        evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
      else if (tbin(lo+k).ge.1) then
        k2 = lo+k
        ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
        kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
        hlp1 = 2.0*term3(k2) - 3.0*term2(k2)
        hlp2 = term3(k2) - term2(k2)
        hlp3 = hlp2 - term2(k2) + term1(k2)
        dhlp1 = 6.0*term2(k2) - 6.0*term1(k2)
        dhlp2 = 3.0*term2(k2) - 2.0*term1(k2)
        dhlp3 = dhlp2 - 2.0*term1(k2) + 1.0
        hlp4 = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k2)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k2)+1) + &
 &               tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k2)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k2)+1))
        evec(9) = evec(9) + scale_TABUL*hlp4
        hlp4 = pfac2*id1(k2)*(dhlp1*(tbp%pot(tbp%lst(i,3),tbin(k2)) - tbp%pot(tbp%lst(i,3),tbin(k2)+1)) + &
 &               tbp%res*(dhlp3*tbp%tang(tbp%lst(i,3),tbin(k2)) + dhlp2*tbp%tang(tbp%lst(i,3),tbin(k2)+1)))
        ca_f(ii,1:3) = ca_f(ii,1:3) + hlp4*dvec(k2,1:3)
        ca_f(kk,1:3) = ca_f(kk,1:3) - hlp4*dvec(k2,1:3)
      else
        evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),1)
      end if
      k = k + 1
    end do
    lo = hi + 1
  end do
!
  end
!
!-----------------------------------------------------------------------------------
!
