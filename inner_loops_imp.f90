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
!------------------------------------------------------------------------------
!
! the vectorized version of en_rsp
! note that the speed-up hinges crucially on the arrays being statically 
! allocated at call-time (for some reason dynamic allocation is terribly slow),
! and on defining the range -> 1:alcsz
!
subroutine Vforce_rsp(evec,rs1,rs2,cut,ca_f)
!
  use iounit
  use inter
  use atoms
  use params
  use energies
  use sequen
  use polypep
  use math
  use cutoffs
  use grandensembles
  use system
  use forces
!
  implicit none
!
  integer rs1,rs2,rs,i,j,k,ii,kk,imol1,imol2,alcsz
  RTYPE d2(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE d1(at(rs1)%na*at(rs2)%na),id2(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE termr(at(rs1)%na*at(rs2)%na)
  RTYPE terms(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  RTYPE term3(at(rs1)%na*at(rs2)%na)
  RTYPE idvn1(at(rs1)%na*at(rs2)%na,3)
  RTYPE tmaxd(at(rs1)%na*at(rs2)%na)
  integer idx2(at(rs1)%na*at(rs2)%na,2)
  RTYPE evec(MAXENERGYTERMS),incr
  RTYPE svec(3),efvoli,efvolk,datri,datrk,pfac,rt23,pwc2,pwc3,derinc
  logical cut,ismember,doid2
  RTYPE for_tmp(at(rs1)%na*at(rs2)%na,3)
  RTYPE ca_f(n,3)
!
! no intialization needed: note that on the calling side it requires a lot of
! care to handle fxns which "silently" increment arguments
!
  if (ideal_run.EQV..true.) then
    return
  end if
!  
! Prevent non-present molecules from interacting with each other
! or with present molecules.  Note that intra-nonpresent-molecule
! interactions are still allowed, in order to correctly simulate
! an infinite-dilution implicit solvent reference state particle
! bath -> partially redundant overlap with use_FEG
  if ((ens%flag.ge.5).AND.(ens%flag.le.6)) then
    imol1 = molofrs(rs1)
    imol2 = molofrs(rs2)
    if ((imol1.ne.imol2).AND.((.not.ismember(ispresent,imol1)).OR. &
 &   (.not.ismember(ispresent,imol2)))) return
  end if
!
! in FEG-1 we branch out if precisely one of the residues is ghosting
! always remember that the interaction of two ghost-residues (incl self!)
! is the full Hamiltonian, not the ghosted one
! in FEG-2 we branch out if either residue is ghosted
  if (use_FEG.EQV..true.) then
    if (fegmode.eq.1) then
   if(((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..false.)).OR.&
 & ((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..true.))) then
        call Vforce_rsp_feg(evec,rs1,rs2,cut,ca_f)
        return
      end if
    else if (fegmode.eq.2) then
      if((par_FEG(rs1).EQV..true.).OR.(par_FEG(rs2).EQV..true.))then
        call Vforce_rsp_feg(evec,rs1,rs2,cut,ca_f)
        return
      end if
    end if
  end if
!
  doid2 = .false.
  if  ((use_attLJ.EQV..true.).OR.&
 &     ((use_IPP.EQV..true.).AND.(use_hardsphere.EQV..false.)).OR.&
 &      (use_WCA.EQV..true.)) then
    doid2 = .true.
  end if
  for_tmp(:,:) = 0.0
!
! three different cases: intra-residue (rs1 == rs2), neighbors in seq., or others
!
  if (rs1.eq.rs2) then
!
!   intra-residue interactions are never subject to PBC (residue size should be within ~10A)
    rs = rs1
!   first set the necessary atom-specific parameters (mimimal)
    k = 1
    do i=1,nrsintra(rs)
      ii = iaa(rs)%atin(i,1)
      kk = iaa(rs)%atin(i,2)
      
      dvec(k,1) = x(kk) - x(ii)
      dvec(k,2) = y(kk) - y(ii)
      dvec(k,3) = z(kk) - z(ii)
      d2(k) = dvec(k,1)**2 + dvec(k,2)**2 + dvec(k,3)**2
      if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.&
 &      (cut.EQV..false.)) then
        term0(k) = fudge(rs)%rsin(i)*fudge(rs)%rsin_lje(i) !lj_eps(attyp(ii),attyp(kk))
        if ((term0(k).gt.0.0).OR.(use_IMPSOLV.EQV..true.)) then
          terms(k) = fudge(rs)%rsin_ljs(i) !lj_sig(attyp(ii),attyp(kk))
          if (use_IMPSOLV.EQV..true.) termr(k) = atr(ii)+atr(kk)
          idx2(k,1) = i
          k = k + 1
        end if
      end if
    end do
    k = k - 1
    alcsz = k ! MArtin : That is the length of the square distance vector d2(k)
    if (alcsz.gt.0) then
      if (doid2.EQV..true.) then
        id2(1:alcsz) = 1.0/d2(1:alcsz) ! Martin : gets the invert of d2 for some reason 
      end if
      if (use_IPP.EQV..true.) then
        term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**nhalf
        term2(1:alcsz) = term0(1:alcsz)*term1(1:alcsz)
        evec(1) = evec(1) + scale_IPP*4.0*&
 &                  sum(term2(1:alcsz))
        pfac = scale_IPP*4.0*nindx
        term3(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*pfac
        for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) - term3(1:alcsz)*dvec(1:alcsz,1)
        for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) - term3(1:alcsz)*dvec(1:alcsz,2)
        for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) - term3(1:alcsz)*dvec(1:alcsz,3)
      end if
      if (use_attLJ.EQV..true.) then
        term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**3
        term2(1:alcsz) = term0(1:alcsz)*term1(1:alcsz)
        evec(3) = evec(3) - scale_attLJ*4.0*&
 &                  sum(term2(1:alcsz))
        pfac = scale_attLJ*24.0
        term3(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*pfac
        for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) + term3(1:alcsz)*dvec(1:alcsz,1)
        for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) + term3(1:alcsz)*dvec(1:alcsz,2)
        for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) + term3(1:alcsz)*dvec(1:alcsz,3)
      end if
      if (use_IMPSOLV.EQV..true.) then
        d1(1:alcsz) = sqrt(d2(1:alcsz)) ! Martin : Actual distance
        idvn1(1:alcsz,1) = dvec(1:alcsz,1)/d1(1:alcsz)
        idvn1(1:alcsz,2) = dvec(1:alcsz,2)/d1(1:alcsz)
        idvn1(1:alcsz,3) = dvec(1:alcsz,3)/d1(1:alcsz)
        tmaxd(1:alcsz) = termr(1:alcsz) + par_IMPSOLV(1) ! Martin sum of atomic radii + solvation sphere thickness (Wait.... should be two par_IMSOLV)

        do k=1,alcsz
          ii = iaa(rs)%atin(idx2(k,1),1)
          kk = iaa(rs)%atin(idx2(k,1),2)
          if (d1(k).lt.tmaxd(k)) then
            efvoli = atsavred(ii)*atvol(ii)
            efvolk = atsavred(kk)*atvol(kk)
            datri = 2.0*atr(ii)
            datrk = 2.0*atr(kk)
            if (d1(k).gt.(tmaxd(k)-datrk)) then
              incr = -(tmaxd(k)-d1(k))/datrk
              derinc = -efvolk/datrk
            else if (d1(k).lt.termr(k)) then
              incr = -d1(k)/termr(k)
              derinc = efvolk/termr(k)
            else
              incr = -1.0
              derinc = 0.0
            end if
#ifdef ENABLE_THREADS
!$OMP CRITICAL(FORCE_SAV)
#endif
            svte(ii) = svte(ii) + incr*efvolk
            if (derinc.ne.0.0) then
              sisa_nr(ii) = sisa_nr(ii) + 1
              sisa_i(sisa_nr(ii),ii) = kk
              sisa_dr(sisa_nr(ii),1,ii) = idvn1(k,1)*derinc
              sisa_dr(sisa_nr(ii),2,ii) = idvn1(k,2)*derinc
              sisa_dr(sisa_nr(ii),3,ii) = idvn1(k,3)*derinc
              sav_dr(ii,1) = sav_dr(ii,1) - sisa_dr(sisa_nr(ii),1,ii)
              sav_dr(ii,2) = sav_dr(ii,2) - sisa_dr(sisa_nr(ii),2,ii)
              sav_dr(ii,3) = sav_dr(ii,3) - sisa_dr(sisa_nr(ii),3,ii)
            end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(FORCE_SAV)
#endif
            if (d1(k).gt.(tmaxd(k)-datri)) then
              incr = -(tmaxd(k)-d1(k))/datri
              derinc = efvoli/datri
            else if (d1(k).lt.termr(k)) then
              incr = -d1(k)/termr(k)
              derinc = -efvoli/termr(k)
            else
              incr = -1.0
              derinc = 0.0
            end if
#ifdef ENABLE_THREADS
!$OMP CRITICAL(FORCE_SAV)
#endif
            svte(kk) = svte(kk) + incr*efvoli
            if (derinc.ne.0.0) then
              sisa_nr(kk) = sisa_nr(kk) + 1
              sisa_i(sisa_nr(kk),kk) = ii
              sisa_dr(sisa_nr(kk),1,kk) = idvn1(k,1)*derinc
              sisa_dr(sisa_nr(kk),2,kk) = idvn1(k,2)*derinc
              sisa_dr(sisa_nr(kk),3,kk) = idvn1(k,3)*derinc
              sav_dr(kk,1) = sav_dr(kk,1) - sisa_dr(sisa_nr(kk),1,kk)
              sav_dr(kk,2) = sav_dr(kk,2) - sisa_dr(sisa_nr(kk),2,kk)
              sav_dr(kk,3) = sav_dr(kk,3) - sisa_dr(sisa_nr(kk),3,kk)
            end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(FORCE_SAV)
#endif
          end if
        end do
      end if
      if (use_WCA.EQV..true.) then
        pfac = scale_WCA*par_WCA(2)
        rt23 = ROOT26*ROOT26
        pwc2 = par_WCA(1)*par_WCA(1)
        pwc3 = 24.0*scale_WCA
        do k=1,alcsz
          if (term0(k).gt.0.0) then
            if (d2(k).lt.rt23*terms(k)) then
              term1(k) = (terms(k)*id2(k))**3
              term2(k) = term1(k)**2
              term3(k) = 4.0*(term2(k) - term1(k)) + 1.0
              evec(5) = evec(5) + scale_WCA*term0(k)*&
 &                               (term3(k) - par_WCA(2))
              term3(k) = -pwc3*(2.0*term2(k) - term1(k))*term0(k)*id2(k)
              for_tmp(k,1) = for_tmp(k,1) + term3(k)*dvec(k,1)
              for_tmp(k,2) = for_tmp(k,2) + term3(k)*dvec(k,2)
              for_tmp(k,3) = for_tmp(k,3) + term3(k)*dvec(k,3)
            else if (d2(k).lt.pwc2*terms(k)) then
              term1(k) = par_WCA(3)*d2(k)/terms(k) + par_WCA(4)
              term2(k) = 0.5*cos(term1(k)) - 0.5
              evec(5) = evec(5) + pfac*term0(k)*term2(k)
              term3(k) = -sin(term1(k))*pfac*term0(k)*par_WCA(3)/terms(k)
              for_tmp(k,1) = for_tmp(k,1) + term3(k)*dvec(k,1)
              for_tmp(k,2) = for_tmp(k,2) + term3(k)*dvec(k,2)
              for_tmp(k,3) = for_tmp(k,3) + term3(k)*dvec(k,3)
            end if
          end if
        end do
      end if
!     now increment net Cartesian force (note we take advantage of the fact that
!     all forces in here are exactly pairwise-decomposable
      do j=1,3
        do k=1,alcsz
          ii = iaa(rs)%atin(idx2(k,1),1)
          kk = iaa(rs)%atin(idx2(k,1),2)
          ca_f(ii,j) = ca_f(ii,j) + for_tmp(k,j)
          ca_f(kk,j) = ca_f(kk,j) - for_tmp(k,j)
        end do
      end do
    end if
!
!
!
  else if (abs(rs1-rs2).eq.1) then
!
!   there is no topological requirement for neighboring residues to be close (pure sequence)
!   so we have to check for BC 
    if (rs1.gt.rs2) then
      rs = rs2
      call dis_bound_rs(rs2,rs1,svec)
    else
      rs = rs1
      call dis_bound_rs(rs1,rs2,svec)
    end if
    k = 1
    do i=1,nrsnb(rs)
      ii = iaa(rs)%atnb(i,1)
      kk = iaa(rs)%atnb(i,2)
      dvec(k,1) = x(kk) - x(ii) + svec(1)
      dvec(k,2) = y(kk) - y(ii) + svec(2)
      dvec(k,3) = z(kk) - z(ii) + svec(3)
      d2(k) = dvec(k,1)**2 + dvec(k,2)**2 + dvec(k,3)**2
      if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.&
 &      (cut.EQV..false.)) then
        term0(k) = fudge(rs)%rsnb(i)*fudge(rs)%rsnb_lje(i) !lj_eps(attyp(ii),attyp(kk))
        if ((term0(k).gt.0.0).OR.(use_IMPSOLV.EQV..true.)) then
          terms(k) = fudge(rs)%rsnb_ljs(i) !lj_sig(attyp(ii),attyp(kk))
          if (use_IMPSOLV.EQV..true.) termr(k) = atr(ii)+atr(kk)
          idx2(k,1) = i
          k = k + 1
        end if
      end if
    end do
    k = k - 1
    alcsz = k
    if (alcsz.gt.0) then
      if (doid2.EQV..true.) then
        id2(1:alcsz) = 1.0/d2(1:alcsz)
      end if
      if (use_IPP.EQV..true.) then
        term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**nhalf
        term2(1:alcsz) = term0(1:alcsz)*term1(1:alcsz)
        evec(1) = evec(1) + scale_IPP*4.0*&
 &                  sum(term2(1:alcsz))
        pfac = scale_IPP*4.0*nindx
        term3(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*pfac
        for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) - term3(1:alcsz)*dvec(1:alcsz,1)
        for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) - term3(1:alcsz)*dvec(1:alcsz,2)
        for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) - term3(1:alcsz)*dvec(1:alcsz,3)
      end if
      if (use_attLJ.EQV..true.) then
        term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**3
        term2(1:alcsz) = term0(1:alcsz)*term1(1:alcsz)
        evec(3) = evec(3) - scale_attLJ*4.0*&
 &                  sum(term2(1:alcsz))
        pfac = scale_attLJ*24.0
        term3(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*pfac
        for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) + term3(1:alcsz)*dvec(1:alcsz,1)
        for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) + term3(1:alcsz)*dvec(1:alcsz,2)
        for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) + term3(1:alcsz)*dvec(1:alcsz,3)
      end if
      if (use_IMPSOLV.EQV..true.) then
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        idvn1(1:alcsz,1) = dvec(1:alcsz,1)/d1(1:alcsz)
        idvn1(1:alcsz,2) = dvec(1:alcsz,2)/d1(1:alcsz)
        idvn1(1:alcsz,3) = dvec(1:alcsz,3)/d1(1:alcsz)
        tmaxd(1:alcsz) = termr(1:alcsz) + par_IMPSOLV(1)
        do k=1,alcsz
          ii = iaa(rs)%atnb(idx2(k,1),1)
          kk = iaa(rs)%atnb(idx2(k,1),2)
          if (d1(k).lt.tmaxd(k)) then
            efvoli = atsavred(ii)*atvol(ii)
            efvolk = atsavred(kk)*atvol(kk)
            datri = 2.0*atr(ii)
            datrk = 2.0*atr(kk)
            if (d1(k).gt.(tmaxd(k)-datrk)) then
              incr = -(tmaxd(k)-d1(k))/datrk
              derinc = -efvolk/datrk
            else if (d1(k).lt.termr(k)) then !Martin : should this really ever happen ?
              incr = -d1(k)/termr(k)        ! Martin : Cause that's lie the definition of atoms overlapping 
              derinc = efvolk/termr(k)
            else
              incr = -1.0
              derinc = 0.0
            end if
#ifdef ENABLE_THREADS
!$OMP CRITICAL(FORCE_SAV)
#endif

            svte(ii) = svte(ii) + incr*efvolk
            if (derinc.ne.0.0) then
              sisa_nr(ii) = sisa_nr(ii) + 1
              sisa_i(sisa_nr(ii),ii) = kk
              sisa_dr(sisa_nr(ii),1,ii) = idvn1(k,1)*derinc
              sisa_dr(sisa_nr(ii),2,ii) = idvn1(k,2)*derinc
              sisa_dr(sisa_nr(ii),3,ii) = idvn1(k,3)*derinc
              sav_dr(ii,1) = sav_dr(ii,1) - sisa_dr(sisa_nr(ii),1,ii)
              sav_dr(ii,2) = sav_dr(ii,2) - sisa_dr(sisa_nr(ii),2,ii)
              sav_dr(ii,3) = sav_dr(ii,3) - sisa_dr(sisa_nr(ii),3,ii)
            end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(FORCE_SAV)
#endif
            if (d1(k).gt.(tmaxd(k)-datri)) then
              incr = -(tmaxd(k)-d1(k))/datri
              derinc = efvoli/datri
            else if (d1(k).lt.termr(k)) then
              incr = -d1(k)/termr(k)
              derinc = -efvoli/termr(k)
            else
              incr = -1.0
              derinc = 0.0
            end if
#ifdef ENABLE_THREADS
!$OMP CRITICAL(FORCE_SAV)
#endif
            svte(kk) = svte(kk) + incr*efvoli
            if (derinc.ne.0.0) then
              sisa_nr(kk) = sisa_nr(kk) + 1
              sisa_i(sisa_nr(kk),kk) = ii
              sisa_dr(sisa_nr(kk),1,kk) = idvn1(k,1)*derinc
              sisa_dr(sisa_nr(kk),2,kk) = idvn1(k,2)*derinc
              sisa_dr(sisa_nr(kk),3,kk) = idvn1(k,3)*derinc
              sav_dr(kk,1) = sav_dr(kk,1) - sisa_dr(sisa_nr(kk),1,kk)
              sav_dr(kk,2) = sav_dr(kk,2) - sisa_dr(sisa_nr(kk),2,kk)
              sav_dr(kk,3) = sav_dr(kk,3) - sisa_dr(sisa_nr(kk),3,kk)
            end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(FORCE_SAV)
#endif
          end if
        end do
      end if
      if (use_WCA.EQV..true.) then
        pfac = scale_WCA*par_WCA(2)
        rt23 = ROOT26*ROOT26
        pwc2 = par_WCA(1)*par_WCA(1)
        pwc3 = 24.0*scale_WCA
        do k=1,alcsz
          if (term0(k).gt.0.0) then
            if (d2(k).lt.rt23*terms(k)) then
              term1(k) = (terms(k)*id2(k))**3
              term2(k) = term1(k)**2
              term3(k) = 4.0*(term2(k) - term1(k)) + 1.0
              evec(5) = evec(5) + scale_WCA*term0(k)*&
 &                               (term3(k) - par_WCA(2))
              term3(k) = -pwc3*(2.0*term2(k) - term1(k))*term0(k)*id2(k)
              for_tmp(k,1) = for_tmp(k,1) + term3(k)*dvec(k,1)
              for_tmp(k,2) = for_tmp(k,2) + term3(k)*dvec(k,2)
              for_tmp(k,3) = for_tmp(k,3) + term3(k)*dvec(k,3)
            else if (d2(k).lt.pwc2*terms(k)) then
              term1(k) = par_WCA(3)*d2(k)/terms(k) + par_WCA(4)
              term2(k) = 0.5*cos(term1(k)) - 0.5
              evec(5) = evec(5) + pfac*term0(k)*term2(k)
              term3(k) = -sin(term1(k))*pfac*term0(k)*par_WCA(3)/terms(k)
              for_tmp(k,1) = for_tmp(k,1) + term3(k)*dvec(k,1)
              for_tmp(k,2) = for_tmp(k,2) + term3(k)*dvec(k,2)
              for_tmp(k,3) = for_tmp(k,3) + term3(k)*dvec(k,3)
            end if
          end if
        end do
      end if
!     now increment net Cartesian force (note we take advantage of the fact that
!     all forces in here are exactly pairwise-decomposable
      do j=1,3
        do k=1,alcsz
          ii = iaa(rs)%atnb(idx2(k,1),1)
          kk = iaa(rs)%atnb(idx2(k,1),2)
          ca_f(ii,j) = ca_f(ii,j) + for_tmp(k,j)
          ca_f(kk,j) = ca_f(kk,j) - for_tmp(k,j)
        end do
      end do
!
    end if
!
!
!
  else
!
!   there is no topological relationship for remaining residues -> always check BC
    call dis_bound_rs(rs1,rs2,svec)
!   loop over complete set of atom-atom interactions in distant residues
    k = 1
    do i=1,at(rs1)%na
      if (i.le.at(rs1)%nbb) then
        ii = at(rs1)%bb(i)
      else
        ii = at(rs1)%sc(i-at(rs1)%nbb)
      end if
      do j=1,at(rs2)%na
        if (j.le.at(rs2)%nbb) then
          kk = at(rs2)%bb(j)
        else
          kk = at(rs2)%sc(j-at(rs2)%nbb)
        end if
        dvec(k,1) = x(kk) - x(ii) + svec(1)
        dvec(k,2) = y(kk) - y(ii) + svec(2)
        dvec(k,3) = z(kk) - z(ii) + svec(3)
        d2(k) = dvec(k,1)**2 + dvec(k,2)**2 + dvec(k,3)**2
        if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.&
 &      (cut.EQV..false.)) then
          term0(k) = lj_eps(attyp(ii),attyp(kk))
          if ((term0(k).gt.0.0).OR.(use_IMPSOLV.EQV..true.)) then
            terms(k) = lj_sig(attyp(ii),attyp(kk))
            if (use_IMPSOLV.EQV..true.) termr(k) = atr(ii)+atr(kk)
            idx2(k,1) = ii
            idx2(k,2) = kk
            k = k + 1
          end if
        end if
      end do
    end do
    k = k - 1
    alcsz = k
    if (alcsz.gt.0) then
      if (doid2.EQV..true.) then
        id2(1:alcsz) = 1.0/d2(1:alcsz)
      end if
      if (use_IPP.EQV..true.) then
        term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**nhalf
        term2(1:alcsz) = term0(1:alcsz)*term1(1:alcsz)
        evec(1) = evec(1) + scale_IPP*4.0*&
 &                  sum(term2(1:alcsz))
        pfac = scale_IPP*4.0*nindx
        term3(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*pfac
        for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) - term3(1:alcsz)*dvec(1:alcsz,1)
        for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) - term3(1:alcsz)*dvec(1:alcsz,2)
        for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) - term3(1:alcsz)*dvec(1:alcsz,3)
      end if
      if (use_attLJ.EQV..true.) then
        term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**3
        term2(1:alcsz) = term0(1:alcsz)*term1(1:alcsz)
        evec(3) = evec(3) - scale_attLJ*4.0*&
 &                  sum(term2(1:alcsz))
        pfac = scale_attLJ*24.0
        term3(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*pfac
        for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) + term3(1:alcsz)*dvec(1:alcsz,1)
        for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) + term3(1:alcsz)*dvec(1:alcsz,2)
        for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) + term3(1:alcsz)*dvec(1:alcsz,3)
      end if
      if (use_IMPSOLV.EQV..true.) then
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        idvn1(1:alcsz,1) = dvec(1:alcsz,1)/d1(1:alcsz)
        idvn1(1:alcsz,2) = dvec(1:alcsz,2)/d1(1:alcsz)
        idvn1(1:alcsz,3) = dvec(1:alcsz,3)/d1(1:alcsz)
        tmaxd(1:alcsz) = termr(1:alcsz) + par_IMPSOLV(1)
        do k=1,alcsz
          ii = idx2(k,1)
          kk = idx2(k,2)
          if (d1(k).lt.tmaxd(k)) then
            efvoli = atsavred(ii)*atvol(ii)
            efvolk = atsavred(kk)*atvol(kk)
            datri = 2.0*atr(ii)
            datrk = 2.0*atr(kk)
            if (d1(k).gt.(tmaxd(k)-datrk)) then
              incr = -(tmaxd(k)-d1(k))/datrk
              derinc = -efvolk/datrk
            else if (d1(k).lt.termr(k)) then
              incr = -d1(k)/termr(k)
              derinc = efvolk/termr(k)
            else
              incr = -1.0
              derinc = 0.0
            end if
#ifdef ENABLE_THREADS
!$OMP CRITICAL(FORCE_SAV)
#endif
            svte(ii) = svte(ii) + incr*efvolk
            if (derinc.ne.0.0) then
              sisa_nr(ii) = sisa_nr(ii) + 1
              sisa_i(sisa_nr(ii),ii) = kk
              sisa_dr(sisa_nr(ii),1,ii) = idvn1(k,1)*derinc
              sisa_dr(sisa_nr(ii),2,ii) = idvn1(k,2)*derinc
              sisa_dr(sisa_nr(ii),3,ii) = idvn1(k,3)*derinc
              sav_dr(ii,1) = sav_dr(ii,1) - sisa_dr(sisa_nr(ii),1,ii)
              sav_dr(ii,2) = sav_dr(ii,2) - sisa_dr(sisa_nr(ii),2,ii)
              sav_dr(ii,3) = sav_dr(ii,3) - sisa_dr(sisa_nr(ii),3,ii)
            end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(FORCE_SAV)
#endif
            if (d1(k).gt.(tmaxd(k)-datri)) then
              incr = -(tmaxd(k)-d1(k))/datri
              derinc = efvoli/datri
            else if (d1(k).lt.termr(k)) then
              incr = -d1(k)/termr(k)
              derinc = -efvoli/termr(k)
            else
              incr = -1.0
              derinc = 0.0
            end if
#ifdef ENABLE_THREADS
!$OMP CRITICAL(FORCE_SAV)
#endif
            svte(kk) = svte(kk) + incr*efvoli
            if (derinc.ne.0.0) then
              sisa_nr(kk) = sisa_nr(kk) + 1
              sisa_i(sisa_nr(kk),kk) = ii
              sisa_dr(sisa_nr(kk),1,kk) = idvn1(k,1)*derinc
              sisa_dr(sisa_nr(kk),2,kk) = idvn1(k,2)*derinc
              sisa_dr(sisa_nr(kk),3,kk) = idvn1(k,3)*derinc
              sav_dr(kk,1) = sav_dr(kk,1) - sisa_dr(sisa_nr(kk),1,kk)
              sav_dr(kk,2) = sav_dr(kk,2) - sisa_dr(sisa_nr(kk),2,kk)
              sav_dr(kk,3) = sav_dr(kk,3) - sisa_dr(sisa_nr(kk),3,kk)
            end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(FORCE_SAV)
#endif
          end if
        end do
      end if
      if (use_WCA.EQV..true.) then
        pfac = scale_WCA*par_WCA(2)
        rt23 = ROOT26*ROOT26
        pwc2 = par_WCA(1)*par_WCA(1)
        pwc3 = 24.0*scale_WCA
        do k=1,alcsz
          ii = idx2(k,1)
          kk = idx2(k,2)
          if (term0(k).gt.0.0) then
            if (d2(k).lt.rt23*terms(k)) then
              term1(k) = (terms(k)*id2(k))**3
              term2(k) = term1(k)**2
              term3(k) = 4.0*(term2(k) - term1(k)) + 1.0
              evec(5) = evec(5) + scale_WCA*term0(k)*&
 &                               (term3(k) - par_WCA(2))
              term3(k) = -pwc3*(2.0*term2(k) - term1(k))*term0(k)*id2(k)
              for_tmp(k,1) = for_tmp(k,1) + term3(k)*dvec(k,1)
              for_tmp(k,2) = for_tmp(k,2) + term3(k)*dvec(k,2)
              for_tmp(k,3) = for_tmp(k,3) + term3(k)*dvec(k,3)
            else if (d2(k).lt.pwc2*terms(k)) then
              term1(k) = par_WCA(3)*d2(k)/terms(k) + par_WCA(4)
              term2(k) = 0.5*cos(term1(k)) - 0.5
              evec(5) = evec(5) + pfac*term0(k)*term2(k)
              term3(k) = -sin(term1(k))*pfac*term0(k)*par_WCA(3)/terms(k)
              for_tmp(k,1) = for_tmp(k,1) + term3(k)*dvec(k,1)
              for_tmp(k,2) = for_tmp(k,2) + term3(k)*dvec(k,2)
              for_tmp(k,3) = for_tmp(k,3) + term3(k)*dvec(k,3)
            end if
          end if
        end do
      end if
      do j=1,3
        do k=1,alcsz
          ii = idx2(k,1)
          kk = idx2(k,2)
          ca_f(ii,j) = ca_f(ii,j) + for_tmp(k,j)
          ca_f(kk,j) = ca_f(kk,j) - for_tmp(k,j)
        end do
      end do
!
    end if
!
!
  end if
!
end
!
!--------------------------------------------------------------------------
!
subroutine Vforce_rsp_feg(evec,rs1,rs2,cut,ca_f)
!
  use iounit
  use inter
  use atoms
  use params
  use energies
  use sequen
  use polypep
  use math
  use cutoffs
  use grandensembles
  use system
  use forces
!
  implicit none
!
  integer rs1,rs2,rs,i,j,k,ii,kk,alcsz
  RTYPE d2(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE id2(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE terms(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  RTYPE term3(at(rs1)%na*at(rs2)%na)
  RTYPE term4(at(rs1)%na*at(rs2)%na)
  RTYPE evec(MAXENERGYTERMS)
  RTYPE svec(3),pfac
  logical cut,doid2
  integer idx2(at(rs1)%na*at(rs2)%na,2)
  RTYPE for_tmp(at(rs1)%na*at(rs2)%na,3)
  RTYPE ca_f(n,3)
! we have a potential override to cover in par_FEG3, which allows residues to be
! fully de-coupled at all times
  if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) return
!
  doid2 = .false.
  if  ((use_FEGS(3).EQV..true.).OR.(use_FEGS(1).EQV..true.)) then
    doid2 = .true.
  else
!   currently (this might change), only these two terms supported
    return
  end if
  for_tmp(:,:) = 0.0
!
! three different cases: intra-residue (rs1 == rs2), neighbors in seq., or others
!
  if (rs1.eq.rs2) then
!
!   intra-residue interactions are never subject to PBC (residue size should be within ~10A)
    rs = rs1
!   first set the necessary atom-specific parameters (mimimal)
    k = 1
    do i=1,nrsintra(rs)
      ii = iaa(rs)%atin(i,1)
      kk = iaa(rs)%atin(i,2)
      dvec(k,1) = x(kk) - x(ii)
      dvec(k,2) = y(kk) - y(ii)
      dvec(k,3) = z(kk) - z(ii)
      d2(k) = dvec(k,1)**2 + dvec(k,2)**2 + dvec(k,3)**2
      if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.&
 &      (cut.EQV..false.)) then
        term0(k) = fudge(rs)%rsin(i)*fudge(rs)%rsin_lje(i) !lj_eps(attyp(ii),attyp(kk))
        if (term0(k).gt.0.0) then
          terms(k) = fudge(rs)%rsin_ljs(i) !lj_sig(attyp(ii),attyp(kk))
          idx2(k,1) = i
          k = k + 1
        end if
      end if
    end do
    k = k - 1
    alcsz = k
    if (alcsz.gt.0) then
      if (doid2.EQV..true.) then
        id2(1:alcsz) = 1.0/d2(1:alcsz)
        term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3
      end if
      if (use_FEGS(1).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(2))
        term3(1:alcsz) = term0(1:alcsz)*term2(1:alcsz)*term2(1:alcsz)
        evec(1) = evec(1) + par_FEG2(1)*4.0*sum(term3(1:alcsz))
        pfac = par_FEG2(1)*48.0
        term4(1:alcsz) = pfac*term3(1:alcsz)*term2(1:alcsz)*term1(1:alcsz)*id2(1:alcsz)
        for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) - term4(1:alcsz)*dvec(1:alcsz,1)
        for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) - term4(1:alcsz)*dvec(1:alcsz,2)
        for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) - term4(1:alcsz)*dvec(1:alcsz,3)
      end if
      if (use_FEGS(3).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(6))
        term3(1:alcsz) = term0(1:alcsz)*term2(1:alcsz)
        evec(3) = evec(3) - par_FEG2(5)*4.0*sum(term3(1:alcsz))
        pfac = par_FEG2(5)*24.0
        term4(1:alcsz) = pfac*term3(1:alcsz)*term2(1:alcsz)*term1(1:alcsz)*id2(1:alcsz)
        for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) + term4(1:alcsz)*dvec(1:alcsz,1)
        for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) + term4(1:alcsz)*dvec(1:alcsz,2)
        for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) + term4(1:alcsz)*dvec(1:alcsz,3)
      end if
!     now increment net Cartesian force (note we take advantage of the fact that
!     all forces in here are exactly pairwise-decomposable
      do j=1,3
        do k=1,alcsz
          ii = iaa(rs)%atin(idx2(k,1),1)
          kk = iaa(rs)%atin(idx2(k,1),2)
          ca_f(ii,j) = ca_f(ii,j) + for_tmp(k,j)
          ca_f(kk,j) = ca_f(kk,j) - for_tmp(k,j)
        end do
      end do
    end if
!
!   add correction terms (in torsional space)
    if (use_CORR.EQV..true.) then
      write(ilog,*) 'Fatal. This is a bug. In the FEG-mode 1, intram&
 &olecular contributions are to be omitted, and in FEG-mode 2 they o&
 &ught to be explicitly excluded. Still encountered computation of t&
 &orsional correction terms. Please report this problem.'
      call fexit()
    end if
!
!
!
  else if (abs(rs1-rs2).eq.1) then
!
!   there is no topological requirement for neighboring residues to be close (pure sequence)
!   so we have to check for BC 
    if (rs1.gt.rs2) then
      rs = rs2
      call dis_bound_rs(rs2,rs1,svec)
    else
      rs = rs1
      call dis_bound_rs(rs1,rs2,svec)
    end if
    k = 1
    do i=1,nrsnb(rs)
      ii = iaa(rs)%atnb(i,1)
      kk = iaa(rs)%atnb(i,2)
      dvec(k,1) = x(kk) - x(ii) + svec(1)
      dvec(k,2) = y(kk) - y(ii) + svec(2)
      dvec(k,3) = z(kk) - z(ii) + svec(3)
      d2(k) = dvec(k,1)**2 + dvec(k,2)**2 + dvec(k,3)**2
      if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.&
 &      (cut.EQV..false.)) then
        term0(k) = fudge(rs)%rsnb(i)*fudge(rs)%rsnb_lje(i) !lj_eps(attyp(ii),attyp(kk))
        if (term0(k).gt.0.0) then
          terms(k) = fudge(rs)%rsnb_ljs(i) !lj_sig(attyp(ii),attyp(kk))
          idx2(k,1) = i
          k = k + 1
        end if
      end if
    end do
    k = k - 1
    alcsz = k
    if (alcsz.gt.0) then
      if (doid2.EQV..true.) then
        id2(1:alcsz) = 1.0/d2(1:alcsz)
        term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3
      end if
      if (use_FEGS(1).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(2))
        term3(1:alcsz) = term0(1:alcsz)*term2(1:alcsz)*term2(1:alcsz)
        evec(1) = evec(1) + par_FEG2(1)*4.0*sum(term3(1:alcsz))
        pfac = par_FEG2(1)*48.0
        term4(1:alcsz) = pfac*term3(1:alcsz)*term2(1:alcsz)*term1(1:alcsz)*id2(1:alcsz)
        for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) - term4(1:alcsz)*dvec(1:alcsz,1)
        for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) - term4(1:alcsz)*dvec(1:alcsz,2)
        for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) - term4(1:alcsz)*dvec(1:alcsz,3)
      end if
      if (use_FEGS(3).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(6))
        term3(1:alcsz) = term0(1:alcsz)*term2(1:alcsz)
        evec(3) = evec(3) - par_FEG2(5)*4.0*sum(term3(1:alcsz))
        pfac = par_FEG2(5)*24.0
        term4(1:alcsz) = pfac*term3(1:alcsz)*term2(1:alcsz)*term1(1:alcsz)*id2(1:alcsz)
        for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) + term4(1:alcsz)*dvec(1:alcsz,1)
        for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) + term4(1:alcsz)*dvec(1:alcsz,2)
        for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) + term4(1:alcsz)*dvec(1:alcsz,3)
      end if
!     now increment net Cartesian force (note we take advantage of the fact that
!     all forces in here are exactly pairwise-decomposable
      do j=1,3
        do k=1,alcsz
          ii = iaa(rs)%atnb(idx2(k,1),1)
          kk = iaa(rs)%atnb(idx2(k,1),2)
          ca_f(ii,j) = ca_f(ii,j) + for_tmp(k,j)
          ca_f(kk,j) = ca_f(kk,j) - for_tmp(k,j)
        end do
      end do
!
    end if
!
!
!
  else
!
!   there is no topological relationship for remaining residues -> always check BC
    call dis_bound_rs(rs1,rs2,svec)
!   loop over complete set of atom-atom interactions in distant residues
    k = 1
    do i=1,at(rs1)%na
      if (i.le.at(rs1)%nbb) then
        ii = at(rs1)%bb(i)
      else
        ii = at(rs1)%sc(i-at(rs1)%nbb)
      end if
      do j=1,at(rs2)%na
        if (j.le.at(rs2)%nbb) then
          kk = at(rs2)%bb(j)
        else
          kk = at(rs2)%sc(j-at(rs2)%nbb)
        end if
        dvec(k,1) = x(kk) - x(ii) + svec(1)
        dvec(k,2) = y(kk) - y(ii) + svec(2)
        dvec(k,3) = z(kk) - z(ii) + svec(3)
        d2(k) = dvec(k,1)**2 + dvec(k,2)**2 + dvec(k,3)**2
        if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.&
 &      (cut.EQV..false.)) then
          term0(k) = lj_eps(attyp(ii),attyp(kk))
          if (term0(k).gt.0.0) then
            terms(k) = lj_sig(attyp(ii),attyp(kk))
            idx2(k,1) = ii
            idx2(k,2) = kk
            k = k + 1
          end if
        end if
      end do
    end do
    k = k - 1
    alcsz = k
    if (alcsz.gt.0) then
      if (doid2.EQV..true.) then
        id2(1:alcsz) = 1.0/d2(1:alcsz)
      end if
      if (doid2.EQV..true.) then
        id2(1:alcsz) = 1.0/d2(1:alcsz)
        term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3
      end if
      if (use_FEGS(1).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(2))
        term3(1:alcsz) = term0(1:alcsz)*term2(1:alcsz)*term2(1:alcsz)
        evec(1) = evec(1) + par_FEG2(1)*4.0*sum(term3(1:alcsz))
        pfac = par_FEG2(1)*48.0
        term4(1:alcsz) = pfac*term3(1:alcsz)*term2(1:alcsz)*term1(1:alcsz)*id2(1:alcsz)
        for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) - term4(1:alcsz)*dvec(1:alcsz,1)
        for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) - term4(1:alcsz)*dvec(1:alcsz,2)
        for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) - term4(1:alcsz)*dvec(1:alcsz,3)
      end if
      if (use_FEGS(3).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(6))
        term3(1:alcsz) = term0(1:alcsz)*term2(1:alcsz)
        evec(3) = evec(3) - par_FEG2(5)*4.0*sum(term3(1:alcsz))
        pfac = par_FEG2(5)*24.0
        term4(1:alcsz) = pfac*term3(1:alcsz)*term2(1:alcsz)*term1(1:alcsz)*id2(1:alcsz)
        for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) + term4(1:alcsz)*dvec(1:alcsz,1)
        for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) + term4(1:alcsz)*dvec(1:alcsz,2)
        for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) + term4(1:alcsz)*dvec(1:alcsz,3)
      end if
      do j=1,3
        do k=1,alcsz
          ii = idx2(k,1)
          kk = idx2(k,2)
          ca_f(ii,j) = ca_f(ii,j) + for_tmp(k,j)
          ca_f(kk,j) = ca_f(kk,j) - for_tmp(k,j)
        end do
      end do
!
    end if
!
!
  end if
!
end
!
!---------------------------------------------------------------------------------
!
! just the LJ potential with support for implicit solvent models
!
subroutine Vforce_LJIMP_C(evec,rs1,ca_f)
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
  use mcsums
  use molecule
!
  implicit none
!
  integer, VOLATILE:: oldsnrii,savnbs
  integer rs1,rs2,ii,kk,j,hi,lo,k,i,shii,kidx
  RTYPE efvoli,datri,idatri,irki,rki,idatrk,help1,incr,derinc
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term6(rs_nbl(rs1)%nnbats)
  RTYPE term6p(rs_nbl(rs1)%nnbats)
  RTYPE term12(rs_nbl(rs1)%nnbats)
  RTYPE tmaxd(rs_nbl(rs1)%nnbats)
  RTYPE tmaxd2(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE d1(rs_nbl(rs1)%nnbats)
  RTYPE atrk(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  RTYPE idvn1(rs_nbl(rs1)%nnbats,3)
  RTYPE epsik(rs_nbl(rs1)%nnbats,at(rs1)%nbb+at(rs1)%nsc)
  RTYPE radik(rs_nbl(rs1)%nnbats,at(rs1)%nbb+at(rs1)%nsc)
  integer savidx(rs_nbl(rs1)%nnbats)
  integer sh(rs_nbl(rs1)%nnbs),hira,lora
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
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh(:) = at(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))%nbb + at(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))%nsc - 1
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
    do ii=0,shii
      epsik(lo:hi,ii+1) = lj_eps(attyp(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k))),attyp(at(rs1)%bb(1)+ii))
      radik(lo:hi,ii+1) = lj_sig(attyp(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k))),attyp(at(rs1)%bb(1)+ii))
    end do
!    epsik(lo:hi,:) = &
! &  lj_eps(attyp(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)),attyp(at(rs1)%bb(1):at(rs1)%bb(1)+shii))
!    radik(lo:hi,:) = &
! &  lj_sig(attyp(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)),attyp(at(rs1)%bb(1):at(rs1)%bb(1)+shii))
    atrk(lo:hi) = atr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
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
    term6p(:) = (radik(:,i)*id2(:))**3
    term12(:) = 4.0*scale_IPP*epsik(:,i)*term6p(:)*term6p(:)
    term6(:) = scale_attLJ*4.0*epsik(:,i)*term6p(:)
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
!   now deal with solvent-accessible volume increments and forces
!   the problem is that here we cannot assume that the bulk of the atoms
!   in the neighbor list contributes -> re-parse to avoid unnecessary computation
    oldsnrii = sisa_nr(ii)
    help1 = atr(ii) + par_IMPSOLV(1)
    tmaxd(:) = atrk(:) + help1
    tmaxd2(:) = tmaxd(:)*tmaxd(:)
    savnbs = 0
    lo = 1
    do k=1,rs_nbl(rs1)%nnbs
      rs2 = rs_nbl(rs1)%nb(k)
      hi = lo + sh(k)
      do kk=lo,hi
        if (dis2(kk).lt.tmaxd2(kk)) then
          kidx = at(rs2)%bb(1)+kk-lo
          savnbs = savnbs + 1
          savidx(savnbs) = kidx
          term12(savnbs) = atvol(kidx)*atsavred(kidx) ! highjacked
          term6(savnbs) = atr(kidx) ! highjacked
          dis2(savnbs) = dis2(kk) ! shift down -> only safe if savnbs is rigorous subset and dis2 not used "normally" below
          dvec(savnbs,1) = dvec(kk,1)
          dvec(savnbs,2) = dvec(kk,2)
          dvec(savnbs,3) = dvec(kk,3) ! shift down -> only safe if savnbs is rigorous subset and dvec not used "normally" below
          tmaxd(savnbs) = tmaxd(kk) ! shift down -> only safe if savnbs is rigorous subset
        end if
      end do
      lo = hi + 1
    end do
!   to be safe (this should never happen of course -> wasteful otherwise)
    if (savnbs.eq.0) cycle
!
    d1(1:savnbs) = sqrt(dis2(1:savnbs))
    term6p(1:savnbs) = 1.0/d1(1:savnbs) ! highjacked
    idvn1(1:savnbs,1) = dvec(1:savnbs,1)*term6p(1:savnbs)
    idvn1(1:savnbs,2) = dvec(1:savnbs,2)*term6p(1:savnbs)
    idvn1(1:savnbs,3) = dvec(1:savnbs,3)*term6p(1:savnbs)
    efvoli = atsavred(ii)*atvol(ii)
    datri = 2.0*atr(ii)
    idatri = 1.0/datri
    do k=1,savnbs
      kk = savidx(k)
      idatrk = 0.5/term6(k)
      rki = (term6(k) + atr(ii))
      irki = 1.0/rki
      if (d1(k).gt.(tmaxd(k)-2.0*term6(k))) then
        incr = -(tmaxd(k)-d1(k))*idatrk
        derinc = -term12(k)*idatrk
      else if (d1(k).lt.rki) then
        incr = -d1(k)*irki
        derinc = term12(k)*irki
      else
        incr = -1.0
        derinc = 0.0
      end if
#ifdef ENABLE_THREADS
!$OMP CRITICAL(FORCE_SAV)
#endif
      svte(ii) = svte(ii) + incr*term12(k)
      if (derinc.ne.0.0) then
        sisa_nr(ii) = sisa_nr(ii) + 1
        sisa_i(sisa_nr(ii),ii) = kk
        sisa_dr(sisa_nr(ii),1,ii) = idvn1(k,1)*derinc
        sisa_dr(sisa_nr(ii),2,ii) = idvn1(k,2)*derinc
        sisa_dr(sisa_nr(ii),3,ii) = idvn1(k,3)*derinc
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(FORCE_SAV)
#endif
!     unfortunately, the terms incrementing the neighbor atoms are now
!     completely striped in memory -> no vectorization possible
      if (d1(k).gt.(tmaxd(k)-datri)) then
        incr = -(tmaxd(k)-d1(k))*idatri
        derinc = efvoli*idatri
      else if (d1(k).lt.rki) then
        incr = -d1(k)*irki
        derinc = -efvoli*irki
      else
        incr = -1.0
        derinc = 0.0
      end if
#ifdef ENABLE_THREADS
!$OMP CRITICAL(FORCE_SAV)
#endif
      svte(kk) = svte(kk) + incr*efvoli
      if (derinc.ne.0.0) then
        sisa_nr(kk) = sisa_nr(kk) + 1
        sisa_i(sisa_nr(kk),kk) = ii
        sisa_dr(sisa_nr(kk),1,kk) = idvn1(k,1)*derinc
        sisa_dr(sisa_nr(kk),2,kk) = idvn1(k,2)*derinc
        sisa_dr(sisa_nr(kk),3,kk) = idvn1(k,3)*derinc
        sav_dr(kk,1) = sav_dr(kk,1) - sisa_dr(sisa_nr(kk),1,kk)
        sav_dr(kk,2) = sav_dr(kk,2) - sisa_dr(sisa_nr(kk),2,kk)
        sav_dr(kk,3) = sav_dr(kk,3) - sisa_dr(sisa_nr(kk),3,kk)
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(FORCE_SAV)
#endif
    end do
!   for ii, at least sav_dr increment can be moved out
    if (sisa_nr(ii).gt.oldsnrii) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(FORCE_SAV)
#endif
      sav_dr(ii,1) = sav_dr(ii,1) - sum(sisa_dr(oldsnrii+1:sisa_nr(ii),1,ii))
      sav_dr(ii,2) = sav_dr(ii,2) - sum(sisa_dr(oldsnrii+1:sisa_nr(ii),2,ii))
      sav_dr(ii,3) = sav_dr(ii,3) - sum(sisa_dr(oldsnrii+1:sisa_nr(ii),3,ii))
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(FORCE_SAV)
#endif
    end if
  end do
!
! we'll always end up here to deal with the kk-forces
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
!---------------------------------------------------------------------------------
!
! just the LJ potential with support for implicit solvent models
!
subroutine Vforce_LJIMP_TR(evec,rs1,ca_f)
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
  integer rs1,rs2,ii,kk,hi,lo,k,i,shii,savnbs,kidx,oldsnrii
  RTYPE efvoli,datri,idatri,irki,rki,idatrk,help1,incr,derinc
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term6(rs_nbl(rs1)%nnbtrats)
  RTYPE term6p(rs_nbl(rs1)%nnbtrats)
  RTYPE term12(rs_nbl(rs1)%nnbtrats)
  RTYPE tmaxd(rs_nbl(rs1)%nnbtrats)
  RTYPE tmaxd2(rs_nbl(rs1)%nnbtrats)
  RTYPE foin(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbtrats,3)
  RTYPE id2(rs_nbl(rs1)%nnbtrats)
  RTYPE dis2(rs_nbl(rs1)%nnbtrats)
  RTYPE d1(rs_nbl(rs1)%nnbtrats)
  RTYPE atrk(rs_nbl(rs1)%nnbtrats)
  RTYPE for_k(rs_nbl(rs1)%nnbtrats,3)
  RTYPE idvn1(rs_nbl(rs1)%nnbtrats,3)
  RTYPE epsik(rs_nbl(rs1)%nnbtrats,at(rs1)%nbb+at(rs1)%nsc)
  RTYPE radik(rs_nbl(rs1)%nnbtrats,at(rs1)%nbb+at(rs1)%nsc)
  integer savidx(rs_nbl(rs1)%nnbtrats)
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
    dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,3)
    do ii=0,shii
      epsik(lo:hi,ii+1) = lj_eps(attyp(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k))),attyp(at(rs1)%bb(1)+ii))
      radik(lo:hi,ii+1) = lj_sig(attyp(at(rs2)%bb(1):(at(rs2)%bb(1)+sh(k))),attyp(at(rs1)%bb(1)+ii))
    end do
!    epsik(lo:hi,:) = &
! &  lj_eps(attyp(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)),attyp(at(rs1)%bb(1):at(rs1)%bb(1)+shii))
!    radik(lo:hi,:) = &
! &  lj_sig(attyp(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)),attyp(at(rs1)%bb(1):at(rs1)%bb(1)+shii))
    atrk(lo:hi) = atr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
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
    term6p(:) = (radik(:,i)*id2(:))**3
    term12(:) = 4.0*scale_IPP*epsik(:,i)*term6p(:)*term6p(:)
    term6(:) = scale_attLJ*4.0*epsik(:,i)*term6p(:)
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
!   now deal with solvent-accessible volume increments and forces
!   the problem is that here we cannot assume that the bulk of the atoms
!   in the neighbor list contributes -> re-parse to avoid unnecessary computation
    help1 = atr(ii) + par_IMPSOLV(1)
    tmaxd(:) = atrk(:) + help1
    tmaxd2(:) = tmaxd(:)*tmaxd(:)
    savnbs = 0
    lo = 1
    do k=1,rs_nbl(rs1)%nnbtrs
      rs2 = rs_nbl(rs1)%nbtr(k)
      hi = lo + sh(k)
      do kk=lo,hi
        if (dis2(kk).lt.tmaxd2(kk)) then
          kidx = at(rs2)%bb(1)+kk-lo
          savnbs = savnbs + 1
          savidx(savnbs) = kidx
          term12(savnbs) = atvol(kidx)*atsavred(kidx) ! highjacked
          term6(savnbs) = atr(kidx) ! highjacked
          dis2(savnbs) = dis2(kk) ! shift down -> only safe if savnbs is rigorous subset and dis2 not used "normally" below
          dvec(savnbs,1) = dvec(kk,1)
          dvec(savnbs,2) = dvec(kk,2)
          dvec(savnbs,3) = dvec(kk,3) ! shift down -> only safe if savnbs is rigorous subset and dvec not used "normally" below
          tmaxd(savnbs) = tmaxd(kk) ! shift down -> only safe if savnbs is rigorous subset
        end if
      end do
      lo = hi + 1
    end do
!   the cycle should always happen, since this is the mid-range function ->
!   support there to be safe (slight cost associated, though)
    if (savnbs.eq.0) cycle
!
    d1(1:savnbs) = sqrt(dis2(1:savnbs))
    term6p(1:savnbs) = 1.0/d1(1:savnbs) ! highjacked
    idvn1(1:savnbs,1) = dvec(1:savnbs,1)*term6p(1:savnbs)
    idvn1(1:savnbs,2) = dvec(1:savnbs,2)*term6p(1:savnbs)
    idvn1(1:savnbs,3) = dvec(1:savnbs,3)*term6p(1:savnbs)
    efvoli = atsavred(ii)*atvol(ii)
    datri = 2.0*atr(ii)
    idatri = 1.0/datri
    oldsnrii = sisa_nr(ii)
    do k=1,savnbs
      kk = savidx(k)
      idatrk = 0.5/term6(k)
      rki = (term6(k) + atr(ii))
      irki = 1.0/rki
      if (d1(k).gt.(tmaxd(k)-2.0*term6(k))) then
        incr = -(tmaxd(k)-d1(k))*idatrk
        derinc = -term12(k)*idatrk
      else if (d1(k).lt.rki) then
        incr = -d1(k)*irki
        derinc = term12(k)*irki
      else
        incr = -1.0
        derinc = 0.0
      end if
#ifdef ENABLE_THREADS
!$OMP CRITICAL(FORCE_SAV)
#endif
      svte(ii) = svte(ii) + incr*term12(k)
      if (derinc.ne.0.0) then
        sisa_nr(ii) = sisa_nr(ii) + 1
        sisa_i(sisa_nr(ii),ii) = kk
        sisa_dr(sisa_nr(ii),1,ii) = idvn1(k,1)*derinc
        sisa_dr(sisa_nr(ii),2,ii) = idvn1(k,2)*derinc
        sisa_dr(sisa_nr(ii),3,ii) = idvn1(k,3)*derinc
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(FORCE_SAV)
#endif
!     unfortunately, the terms incrementing the neighbor atoms are now
!     completely striped in memory -> no vectorization possible
      if (d1(k).gt.(tmaxd(k)-datri)) then
        incr = -(tmaxd(k)-d1(k))*idatri
        derinc = efvoli*idatri
      else if (d1(k).lt.rki) then
        incr = -d1(k)*irki
        derinc = -efvoli*irki
      else
        incr = -1.0
        derinc = 0.0
      end if
#ifdef ENABLE_THREADS
!$OMP CRITICAL(FORCE_SAV)
#endif
      svte(kk) = svte(kk) + incr*efvoli
      if (derinc.ne.0.0) then
        sisa_nr(kk) = sisa_nr(kk) + 1
        sisa_i(sisa_nr(kk),kk) = ii
        sisa_dr(sisa_nr(kk),1,kk) = idvn1(k,1)*derinc
        sisa_dr(sisa_nr(kk),2,kk) = idvn1(k,2)*derinc
        sisa_dr(sisa_nr(kk),3,kk) = idvn1(k,3)*derinc
        sav_dr(kk,1) = sav_dr(kk,1) - sisa_dr(sisa_nr(kk),1,kk)
        sav_dr(kk,2) = sav_dr(kk,2) - sisa_dr(sisa_nr(kk),2,kk)
        sav_dr(kk,3) = sav_dr(kk,3) - sisa_dr(sisa_nr(kk),3,kk)
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(FORCE_SAV)
#endif
    end do
!   for ii, at least sav_dr increment can be moved out
    if (sisa_nr(ii).gt.oldsnrii) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(FORCE_SAV)
#endif
      sav_dr(ii,1) = sav_dr(ii,1) - sum(sisa_dr(oldsnrii+1:sisa_nr(ii),1,ii))
      sav_dr(ii,2) = sav_dr(ii,2) - sum(sisa_dr(oldsnrii+1:sisa_nr(ii),2,ii))
      sav_dr(ii,3) = sav_dr(ii,3) - sum(sisa_dr(oldsnrii+1:sisa_nr(ii),3,ii))
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(FORCE_SAV)
#endif
    end if
  end do
!
! we'll always end up here to deal with the kk-forces
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
!-----------------------------------------------------------------------
!
! the vectorized version of force_rsp_long
! note that the speed-up hinges crucially on the arrays being statically 
! allocated at call-time (for some reason dynamic allocation is terribly slow),
! and on defining the range -> 1:alcsz
! note that unlike in en_rsp_long, long-range electrostatics is handled externally(!)
!
subroutine Vforce_rsp_long(evec,rs1,rs2,cut,ca_f,sum_s)
!
  use iounit
  use energies
  use atoms
  use polypep
  use inter
  use units
  use tabpot
  use molecule
  use sequen
  use cutoffs
  use system
  use grandensembles
  use forces
!
  implicit none
!
  integer rs,rs1,rs2,i,j,k,ii,kk,rhi,rlo,imol1,imol2,alcsz
  RTYPE sum_s(n)
  RTYPE evec(MAXENERGYTERMS)
  RTYPE svec(3),pfac,pfac2,pfac3,pfac4,dum1,dd1,dd2,hlp1,hlp2,hlp3,hlp4,dhlp1,dhlp2,dhlp3
  RTYPE d2(at(rs1)%na*at(rs2)%na),id2(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE d1(at(rs1)%na*at(rs2)%na),id1(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE termr(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  RTYPE term3(at(rs1)%na*at(rs2)%na)
  RTYPE term4(at(rs1)%na*at(rs2)%na)
  RTYPE termq(at(rs1)%na*at(rs2)%na),termqi(at(rs1)%na*at(rs2)%na),termqk(at(rs1)%na*at(rs2)%na)
  RTYPE terms(at(rs1)%na*at(rs2)%na),terms1(at(rs1)%na*at(rs2)%na),terms2(at(rs1)%na*at(rs2)%na)
  RTYPE for_tmpi(at(rs1)%na*at(rs2)%na,3)
  RTYPE for_tmpk(at(rs1)%na*at(rs2)%na,3)
  RTYPE tsdi(at(rs1)%na*at(rs2)%na,3)
  RTYPE tsdk(at(rs1)%na*at(rs2)%na,3)
  RTYPE termf(at(rs1)%na*at(rs2)%na,3)
  integer tbin(at(rs1)%na*at(rs2)%na)
  logical cut,ismember
  RTYPE ca_f(n,3)
!
  if ((ideal_run.EQV..true.).OR.&
 &    ((use_POLAR.EQV..false.).AND.(use_TABUL.EQV..false.))) then
    return
  end if
!  
! Prevent non-present molecules from interacting with each other
! or with present molecules.  Note that intra-nonpresent-molecule
! interactions are still allowed, in order to correctly simulate
! an infinite-dilution implicit solvent reference state particle
! bath -> partially redundant overlap with use_FEG
  if ((ens%flag.ge.5).AND.(ens%flag.le.6)) then
    imol1 = molofrs(rs1)
    imol2 = molofrs(rs2)
    if ((imol1.ne.imol2).AND.((.not.ismember(ispresent,imol1)).OR. &
 &   (.not.ismember(ispresent,imol2)))) return
  end if
!
! in FEG we branch out if precisely one of the residues is ghosting
! always remember that the interaction of two ghost-residues (incl self!)
! is the full Hamiltonian, not the ghosted one
  if (use_FEG.EQV..true.) then
    if (fegmode.eq.1) then
   if(((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..false.)).OR.&
 & ((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..true.))) then
        call Vforce_rsp_long_feg(evec,rs1,rs2,cut,ca_f,sum_s)
        return
      end if
    else if (fegmode.eq.2) then
      if((par_FEG(rs1).EQV..true.).OR.(par_FEG(rs2).EQV..true.))then
        call Vforce_rsp_long_feg(evec,rs1,rs2,cut,ca_f,sum_s)
        return
      end if
    end if
  end if
!
  svec(:) = 0.0
!
!
! three different cases: intra-residue (rs1 == rs2), neighbors in seq., or others
!
  if (rs1.eq.rs2) then
    rs = rs1
    if (use_POLAR.EQV..true.) then
!
!     first set the necessary atom-specific parameters (mimimal)
      k = 0
      pfac = electric*scale_POLAR
      alcsz = nrpolintra(rs)
      if (alcsz.gt.0) then
        if (use_IMPSOLV.EQV..true.) then
          if (scrq_model.le.2) then
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2) 
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
              terms(i) = scrq(ii)*scrq(kk)
              tsdk(i,1:3) = scrq_dr(kk,1:3)*scrq(ii)
              tsdi(i,1:3) = scrq_dr(ii,1:3)*scrq(kk)
              termqk(i) = termq(i)*scrq(ii)
              termqi(i) = termq(i)*scrq(kk)
            end do
          else if (scrq_model.eq.4) then
            pfac2 = par_IMPSOLV(8)*electric*scale_POLAR
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
              terms(i) = atr(ii)+atr(kk)
            end do
          else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
            pfac2 = par_IMPSOLV(8)*electric
            pfac3 = 1./par_IMPSOLV(1)
            pfac4 = par_IMPSOLV(9)*pfac3
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
              terms(i) = scrq(ii)*scrq(kk)
              tsdk(i,1:3) = scrq_dr(kk,1:3)*scrq(ii)
              tsdi(i,1:3) = scrq_dr(ii,1:3)*scrq(kk)
              termqk(i) = termq(i)*scrq(ii)
              termqi(i) = termq(i)*scrq(kk)
              termr(i) = atr(ii)+atr(kk)
            end do
          else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
              terms1(i) = scrq(ii)
              terms2(i) = scrq(kk)
              tsdk(i,1:3) = scrq_dr(kk,1:3)
              tsdi(i,1:3) = scrq_dr(ii,1:3)
            end do
          else ! 7,8
            pfac2 = par_IMPSOLV(8)*electric
            pfac3 = 1./par_IMPSOLV(1)
            pfac4 = par_IMPSOLV(9)*pfac3
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
              terms1(i) = scrq(ii)
              terms2(i) = scrq(kk)
              tsdk(i,1:3) = scrq_dr(kk,1:3)
              tsdi(i,1:3) = scrq_dr(ii,1:3)
              termr(i) = atr(ii)+atr(kk)
            end do
          end if
        else
          do i=1,nrpolintra(rs)
            ii = iaa(rs)%polin(i,1)
            kk = iaa(rs)%polin(i,2)
            dvec(i,1) = x(kk) - x(ii)
            dvec(i,2) = y(kk) - y(ii)
            dvec(i,3) = z(kk) - z(ii)
            termq(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
          end do
        end if
!       now perform the vectorizable bulk operations (maximal)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                    dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id2(1:alcsz) = id1(1:alcsz)*id1(1:alcsz)
!       finish up -> screening model specificity is high
        if (use_IMPSOLV.EQV..true.) then
          k = 0
!         with a minor if-statement, the normal models (1,2,5,6) can be handled together
          if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
              call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
              call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
              tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
              tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
              tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
              tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
              tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
              tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
              termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
              termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
            end if
            term0(1:alcsz) = pfac*id1(1:alcsz)
            termqk(1:alcsz) = termqk(1:alcsz)*term0(1:alcsz)
            termqi(1:alcsz) = termqi(1:alcsz)*term0(1:alcsz)
            term0(1:alcsz) = term0(1:alcsz)*termq(1:alcsz)
            for_tmpi(1:alcsz,1) = term0(1:alcsz)*tsdi(1:alcsz,1)
            for_tmpi(1:alcsz,2) = term0(1:alcsz)*tsdi(1:alcsz,2)
            for_tmpi(1:alcsz,3) = term0(1:alcsz)*tsdi(1:alcsz,3)
            for_tmpk(1:alcsz,1) = term0(1:alcsz)*tsdk(1:alcsz,1)
            for_tmpk(1:alcsz,2) = term0(1:alcsz)*tsdk(1:alcsz,2)
            for_tmpk(1:alcsz,3) = term0(1:alcsz)*tsdk(1:alcsz,3)
            term1(1:alcsz) = term0(1:alcsz)*terms(1:alcsz)
            evec(6) = evec(6) + sum(term1(1:alcsz))
            term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
            termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
              sum_s(ii) = sum_s(ii) + termqi(i)
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
              sum_s(kk) = sum_s(kk) + termqk(i)
              ca_f(ii,1:3) = ca_f(ii,1:3) - termf(i,1:3) + for_tmpi(i,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + termf(i,1:3) + for_tmpk(i,1:3)
            end do
!         the simplest one is always the straight distance-dependence (only complexity
!         comes from contact dielectric)
          else if (scrq_model.eq.4) then
            do i=1,nrpolintra(rs)
              if (d1(i).ge.terms(i)) then
                term1(i) = id1(i)
                term0(i) = 2.0*id1(i)
              else
                term1(i) = 1.0/terms(i)
                term0(i) = term1(i)
              end if
            end do
            term2(1:alcsz) = pfac2*termq(1:alcsz)*id1(1:alcsz)
            termqi(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*term0(1:alcsz)
            term2(1:alcsz) = term2(1:alcsz)*term1(1:alcsz)
            termf(1:alcsz,1) = termqi(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = termqi(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = termqi(1:alcsz)*dvec(1:alcsz,3)
            evec(6) = evec(6) + sum(term2(1:alcsz))
            k = 0
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              ca_f(ii,1:3) = ca_f(ii,1:3) - termf(i,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + termf(i,1:3)
            end do
!         the remaining models (3,7,8,9) have spliced-in distance-dependencies
!         this leads to bulkier-looking code but does introduce that much extra complexity compared
!         to the respective reference models (1,5,6,2 in that order)
          else
            if ((scrq_model.eq.7).OR.(scrq_model.eq.8)) then
              call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
              call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
              tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
              tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
              tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
              tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
              tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
              tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
              termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
              termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
            end if
            term0(1:alcsz) = pfac*id1(1:alcsz) ! cc4 without atq-terms
            term1(1:alcsz) = electric*termq(1:alcsz)*terms(1:alcsz) ! ff1
            term2(1:alcsz) = pfac2*termq(1:alcsz)/termr(1:alcsz) ! ff2
            for_tmpi(1:alcsz,1:3) = 0.0
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              if (abs(term1(i)).gt.abs(term2(i))) then
                term4(i) = term1(i)
                term3(i) = 1.0
              else
                dd1 = (d1(i) - termr(i))*pfac3
                dd2 = 1.0 - dd1
                if (dd1.lt.0.0) then
                  term3(i) = (1.0 - par_IMPSOLV(9))
                  term4(i) = term3(i)*term1(i) + par_IMPSOLV(9)*term2(i)
                  term0(i) = term3(i)*term0(i)
                else if (dd2.gt.0.0) then
                  dum1 = par_IMPSOLV(9)*dd2
                  term3(i) = (1.0 - dum1)
                  term4(i) = term3(i)*term1(i) + dum1*term2(i)
                  term0(i) = term3(i)*term0(i)
                  dum1 = -scale_POLAR*id2(i)*(term2(i)-term1(i))*pfac4
                  for_tmpi(i,1) = dum1*dvec(i,1)
                  for_tmpi(i,2) = dum1*dvec(i,2)
                  for_tmpi(i,3) = dum1*dvec(i,3)
                else
                  term4(i) = term1(i)
                  term3(i) = 1.0
                end if
              end if
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
              sum_s(ii) = sum_s(ii) + term0(i)*termqi(i)
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
              sum_s(kk) = sum_s(kk) + term0(i)*termqk(i)
            end do
            term0(1:alcsz) = pfac*term3(1:alcsz)*termq(1:alcsz) ! overwrite of cc4
            for_tmpk(1:alcsz,1) = -for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdk(1:alcsz,1)*id1(1:alcsz)
            for_tmpk(1:alcsz,2) = -for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdk(1:alcsz,2)*id1(1:alcsz)
            for_tmpk(1:alcsz,3) = -for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdk(1:alcsz,3)*id1(1:alcsz)
            for_tmpi(1:alcsz,1) = for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdi(1:alcsz,1)*id1(1:alcsz)
            for_tmpi(1:alcsz,2) = for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdi(1:alcsz,2)*id1(1:alcsz)
            for_tmpi(1:alcsz,3) = for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdi(1:alcsz,3)*id1(1:alcsz)
            term1(1:alcsz) = scale_POLAR*term4(1:alcsz)*id1(1:alcsz)
            evec(6) = evec(6) + sum(term1(1:alcsz))
            term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
            termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
            k = 0
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              ca_f(ii,1:3) = ca_f(ii,1:3) - termf(i,1:3) + for_tmpi(i,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + termf(i,1:3) + for_tmpk(i,1:3)
            end do
          end if
        else
          term1(1:alcsz) = pfac*termq(1:alcsz)*id1(1:alcsz)
          evec(6) = evec(6) + sum(term1(1:alcsz))
          term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
          do i=1,nrpolintra(rs)
            ii = iaa(rs)%polin(i,1)
            kk = iaa(rs)%polin(i,2)
            ca_f(ii,1:3) = ca_f(ii,1:3) - term2(i)*dvec(i,1:3)
            ca_f(kk,1:3) = ca_f(kk,1:3) + term2(i)*dvec(i,1:3)
          end do
        end if
      end if
!
    end if
!
!   note that this might be terribly inefficient if POLAR and TABUL are used concurrently
    if (use_TABUL.EQV..true.) then
!
      if (tbp%rsmat(rs1,rs2).gt.0) then
        rhi = max(rs1,rs2)
        rlo = min(rs1,rs2)
        k = 0
        alcsz = tbp%rsmat(rhi,rlo)-tbp%rsmat(rlo,rhi)+1
        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
          ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
          kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
          k = k + 1
          dvec(k,1) = x(kk) - x(ii)
          dvec(k,2) = y(kk) - y(ii)
          dvec(k,3) = z(kk) - z(ii)
        end do
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        pfac = 1.0/tbp%res
        term0(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
        tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(term0(1:alcsz))))
        term1(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
        term2(1:alcsz) = term1(1:alcsz)*term1(1:alcsz) ! quadratic
        term3(1:alcsz) = term2(1:alcsz)*term1(1:alcsz) ! cubic
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        pfac2 = pfac*scale_TABUL
        k = 0
        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
          k = k + 1
          if (tbin(k).ge.tbp%bins) then
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
          else if (tbin(k).ge.1) then
            ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
            kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
            hlp1 = 2.0*term3(k) - 3.0*term2(k)
            hlp2 = term3(k) - term2(k)
            hlp3 = hlp2 - term2(k) + term1(k)
            dhlp1 = 6.0*term2(k) - 6.0*term1(k)
            dhlp2 = 3.0*term2(k) - 2.0*term1(k)
            dhlp3 = dhlp2 - 2.0*term1(k) + 1.0
            hlp4 = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k)+1) + &
 &               tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1))
            evec(9) = evec(9) + scale_TABUL*hlp4
            hlp4 = pfac2*id1(k)*(dhlp1*(tbp%pot(tbp%lst(i,3),tbin(k)) - tbp%pot(tbp%lst(i,3),tbin(k)+1)) + &
 &               tbp%res*(dhlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + dhlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1)))
            ca_f(ii,1:3) = ca_f(ii,1:3) + hlp4*dvec(k,1:3)
            ca_f(kk,1:3) = ca_f(kk,1:3) - hlp4*dvec(k,1:3)
          else
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),1)
          end if
        end do
      end if
    end if
!
!
!
! neighboring residues
!
  else if (abs(rs1-rs2).eq.1) then
!
!   there is no topological requirement for neighboring residues to be close (pure sequence)
!   so we have to check for BC evtl.y
!
    if (rs1.gt.rs2) then
      rs = rs2
      call dis_bound_rs(rs2,rs1,svec)
    else
      rs = rs1
      call dis_bound_rs(rs1,rs2,svec)
    end if
!
    if (use_POLAR.EQV..true.) then
!
!     first set the necessary atom-specific parameters (mimimal)
      k = 0
      pfac = electric*scale_POLAR
      alcsz = nrpolnb(rs)
      if (alcsz.gt.0) then
        if (use_IMPSOLV.EQV..true.) then
          if (scrq_model.le.2) then
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2) 
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
              terms(i) = scrq(ii)*scrq(kk)
              tsdk(i,1:3) = scrq_dr(kk,1:3)*scrq(ii)
              tsdi(i,1:3) = scrq_dr(ii,1:3)*scrq(kk)
              termqk(i) = termq(i)*scrq(ii)
              termqi(i) = termq(i)*scrq(kk)
            end do
          else if (scrq_model.eq.4) then
            pfac2 = par_IMPSOLV(8)*electric*scale_POLAR
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
              terms(i) = atr(ii)+atr(kk)
            end do
          else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
            pfac2 = par_IMPSOLV(8)*electric
            pfac3 = 1./par_IMPSOLV(1)
            pfac4 = par_IMPSOLV(9)*pfac3
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
              terms(i) = scrq(ii)*scrq(kk)
              tsdk(i,1:3) = scrq_dr(kk,1:3)*scrq(ii)
              tsdi(i,1:3) = scrq_dr(ii,1:3)*scrq(kk)
              termqk(i) = termq(i)*scrq(ii)
              termqi(i) = termq(i)*scrq(kk)
              termr(i) = atr(ii)+atr(kk)
            end do
          else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
              terms1(i) = scrq(ii)
              terms2(i) = scrq(kk)
              tsdk(i,1:3) = scrq_dr(kk,1:3)
              tsdi(i,1:3) = scrq_dr(ii,1:3)
            end do
          else ! 7,8
            pfac2 = par_IMPSOLV(8)*electric
            pfac3 = 1./par_IMPSOLV(1)
            pfac4 = par_IMPSOLV(9)*pfac3
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
              terms1(i) = scrq(ii)
              terms2(i) = scrq(kk)
              tsdk(i,1:3) = scrq_dr(kk,1:3)
              tsdi(i,1:3) = scrq_dr(ii,1:3)
              termr(i) = atr(ii)+atr(kk)
            end do
          end if
        else
          do i=1,nrpolnb(rs)
            ii = iaa(rs)%polnb(i,1)
            kk = iaa(rs)%polnb(i,2)
            dvec(i,1) = x(kk) - x(ii)
            dvec(i,2) = y(kk) - y(ii)
            dvec(i,3) = z(kk) - z(ii)
            termq(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
          end do
        end if
!       now perform the vectorizable bulk operations (maximal)
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                    dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id2(1:alcsz) = id1(1:alcsz)*id1(1:alcsz)
!       finish up -> screening model specificity is high
        if (use_IMPSOLV.EQV..true.) then
          k = 0
          if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
              call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
              call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
              tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
              tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
              tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
              tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
              tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
              tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
              termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
              termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
            end if
            term0(1:alcsz) = pfac*id1(1:alcsz)
            termqk(1:alcsz) = termqk(1:alcsz)*term0(1:alcsz)
            termqi(1:alcsz) = termqi(1:alcsz)*term0(1:alcsz)
            term0(1:alcsz) = term0(1:alcsz)*termq(1:alcsz)
            for_tmpi(1:alcsz,1) = term0(1:alcsz)*tsdi(1:alcsz,1)
            for_tmpi(1:alcsz,2) = term0(1:alcsz)*tsdi(1:alcsz,2)
            for_tmpi(1:alcsz,3) = term0(1:alcsz)*tsdi(1:alcsz,3)
            for_tmpk(1:alcsz,1) = term0(1:alcsz)*tsdk(1:alcsz,1)
            for_tmpk(1:alcsz,2) = term0(1:alcsz)*tsdk(1:alcsz,2)
            for_tmpk(1:alcsz,3) = term0(1:alcsz)*tsdk(1:alcsz,3)
            term1(1:alcsz) = term0(1:alcsz)*terms(1:alcsz)
            evec(6) = evec(6) + sum(term1(1:alcsz))
            term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
            termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
              sum_s(ii) = sum_s(ii) + termqi(i)
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
              sum_s(kk) = sum_s(kk) + termqk(i)
              ca_f(ii,1:3) = ca_f(ii,1:3) + for_tmpi(i,1:3) - termf(i,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + for_tmpk(i,1:3) + termf(i,1:3)
            end do
          else if (scrq_model.eq.4) then
            do i=1,nrpolnb(rs)
              if (d1(i).ge.terms(i)) then
                term1(i) = id1(i)
                term0(i) = 2.0*id1(i)
              else
                term1(i) = 1.0/terms(i)
                term0(i) = term1(i)
              end if
            end do
            term2(1:alcsz) = pfac2*termq(1:alcsz)*id1(1:alcsz)
            termqi(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*term0(1:alcsz)
            term2(1:alcsz) = term2(1:alcsz)*term1(1:alcsz)
            termf(1:alcsz,1) = termqi(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = termqi(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = termqi(1:alcsz)*dvec(1:alcsz,3)
            evec(6) = evec(6) + sum(term2(1:alcsz))
            k = 0
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              ca_f(ii,1:3) = ca_f(ii,1:3) - termf(i,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + termf(i,1:3)
            end do
          else
            if ((scrq_model.eq.7).OR.(scrq_model.eq.8)) then
              call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
              call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
              tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
              tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
              tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
              tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
              tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
              tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
              termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
              termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
            end if
            term0(1:alcsz) = pfac*id1(1:alcsz) ! cc4 without atq-terms
            term1(1:alcsz) = electric*termq(1:alcsz)*terms(1:alcsz) ! ff1
            term2(1:alcsz) = pfac2*termq(1:alcsz)/termr(1:alcsz) ! ff2
            for_tmpi(1:alcsz,1:3) = 0.0
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              if (abs(term1(i)).gt.abs(term2(i))) then
                term4(i) = term1(i)
                term3(i) = 1.0
              else
                dd1 = (d1(i) - termr(i))*pfac3
                dd2 = 1.0 - dd1
                if (dd1.lt.0.0) then
                  term3(i) = (1.0 - par_IMPSOLV(9))
                  term4(i) = term3(i)*term1(i) + par_IMPSOLV(9)*term2(i)
                  term0(i) = term3(i)*term0(i)
                else if (dd2.gt.0.0) then
                  dum1 = par_IMPSOLV(9)*dd2
                  term3(i) = (1.0 - dum1)
                  term4(i) = term3(i)*term1(i) + dum1*term2(i)
                  term0(i) = term3(i)*term0(i)
                  dum1 = -scale_POLAR*id2(i)*(term2(i)-term1(i))*pfac4
                  for_tmpi(i,1) = dum1*dvec(i,1)
                  for_tmpi(i,2) = dum1*dvec(i,2)
                  for_tmpi(i,3) = dum1*dvec(i,3)
                else
                  term4(i) = term1(i)
                  term3(i) = 1.0
                end if
              end if
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
              sum_s(ii) = sum_s(ii) + term0(i)*termqi(i)
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
              sum_s(kk) = sum_s(kk) + term0(i)*termqk(i)
            end do
            term0(1:alcsz) = pfac*term3(1:alcsz)*termq(1:alcsz) ! overwrite of cc4
            for_tmpk(1:alcsz,1) = -for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdk(1:alcsz,1)*id1(1:alcsz)
            for_tmpk(1:alcsz,2) = -for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdk(1:alcsz,2)*id1(1:alcsz)
            for_tmpk(1:alcsz,3) = -for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdk(1:alcsz,3)*id1(1:alcsz)
            for_tmpi(1:alcsz,1) = for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdi(1:alcsz,1)*id1(1:alcsz)
            for_tmpi(1:alcsz,2) = for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdi(1:alcsz,2)*id1(1:alcsz)
            for_tmpi(1:alcsz,3) = for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdi(1:alcsz,3)*id1(1:alcsz)
            term1(1:alcsz) = scale_POLAR*term4(1:alcsz)*id1(1:alcsz)
            evec(6) = evec(6) + sum(term1(1:alcsz))
            term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
            termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
            k = 0
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              ca_f(ii,1:3) = ca_f(ii,1:3) - termf(i,1:3) + for_tmpi(i,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + termf(i,1:3) + for_tmpk(i,1:3)
            end do
          end if
        else
          term1(1:alcsz) = pfac*termq(1:alcsz)*id1(1:alcsz)
          evec(6) = evec(6) + sum(term1(1:alcsz))
          term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
          do i=1,nrpolnb(rs)
            ii = iaa(rs)%polnb(i,1)
            kk = iaa(rs)%polnb(i,2)
            ca_f(ii,1:3) = ca_f(ii,1:3) - term2(i)*dvec(i,1:3)
            ca_f(kk,1:3) = ca_f(kk,1:3) + term2(i)*dvec(i,1:3)
          end do
        end if
      end if
!
    end if
!
    if (use_TABUL.EQV..true.) then
      if (tbp%rsmat(rs1,rs2).gt.0) then
        rhi = max(rs1,rs2)
        rlo = min(rs1,rs2)
        k = 0
        alcsz = tbp%rsmat(rhi,rlo) - tbp%rsmat(rlo,rhi) + 1
        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
          ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
          kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
          k = k + 1
          dvec(k,1) = x(kk) - x(ii)
          dvec(k,2) = y(kk) - y(ii)
          dvec(k,3) = z(kk) - z(ii)
        end do
!       note we need to correct: the svec is adjusted such that the
!       sign is correct when using the polnb(rs)-array, but that TAB
!       assumes rs1,rs2 order
        if (rs.eq.rs2) then
          dvec(1:alcsz,1) = dvec(1:alcsz,1) - svec(1)
          dvec(1:alcsz,2) = dvec(1:alcsz,2) - svec(2)
          dvec(1:alcsz,3) = dvec(1:alcsz,3) - svec(3)
        else
          dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
          dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
          dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        end if
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        pfac = 1.0/tbp%res
        term0(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
        tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(term0(1:alcsz))))
        term1(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
        term2(1:alcsz) = term1(1:alcsz)*term1(1:alcsz) ! quadratic
        term3(1:alcsz) = term2(1:alcsz)*term1(1:alcsz) ! cubic
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        pfac2 = pfac*scale_TABUL
        k = 0
        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
          k = k + 1
          if (tbin(k).ge.tbp%bins) then
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
          else if (tbin(k).ge.1) then
            ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
            kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)

            hlp1 = 2.0*term3(k) - 3.0*term2(k)
            hlp2 = term3(k) - term2(k)
            hlp3 = hlp2 - term2(k) + term1(k)
            dhlp1 = 6.0*term2(k) - 6.0*term1(k)
            dhlp2 = 3.0*term2(k) - 2.0*term1(k)
            dhlp3 = dhlp2 - 2.0*term1(k) + 1.0
            hlp4 = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k)+1) + &
 &               tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1))
            evec(9) = evec(9) + scale_TABUL*hlp4
            hlp4 = pfac2*id1(k)*(dhlp1*(tbp%pot(tbp%lst(i,3),tbin(k)) - tbp%pot(tbp%lst(i,3),tbin(k)+1)) + &
 &               tbp%res*(dhlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + dhlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1)))
            ca_f(ii,1:3) = ca_f(ii,1:3) + hlp4*dvec(k,1:3)
            ca_f(kk,1:3) = ca_f(kk,1:3) - hlp4*dvec(k,1:3)
          else
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),1)
          end if
        end do
      end if
    end if
!
!
!
! all other residue pairs
!
 67   format(4g14.6)
  else
!
!   there is no topological relationship for remaining residues -> always check BC
    call dis_bound_rs(rs1,rs2,svec)
!
    if (use_POLAR.EQV..true.) then
!
!     first set the necessary atom-specific parameters (mimimal)
      k = 0
      pfac = electric*scale_POLAR
      alcsz = at(rs1)%npol*at(rs2)%npol
      if (alcsz.gt.0) then
        if (use_IMPSOLV.EQV..true.) then
          if (scrq_model.le.2) then
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)
                dvec(k,2) = y(kk) - y(ii)
                dvec(k,3) = z(kk) - z(ii)
                termq(k) = atq(ii)*atq(kk)
                terms(k) = scrq(ii)*scrq(kk)
                tsdk(k,1:3) = scrq_dr(kk,1:3)*scrq(ii)
                tsdi(k,1:3) = scrq_dr(ii,1:3)*scrq(kk)
                termqk(k) = termq(k)*scrq(ii)
                termqi(k) = termq(k)*scrq(kk)
              end do
            end do
          else if (scrq_model.eq.4) then
            pfac2 = par_IMPSOLV(8)*electric*scale_POLAR
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)
                dvec(k,2) = y(kk) - y(ii)
                dvec(k,3) = z(kk) - z(ii)
                termq(k) = atq(ii)*atq(kk)
                terms(k) = atr(ii)+atr(kk)
              end do
            end do
          else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
            pfac2 = par_IMPSOLV(8)*electric
            pfac3 = 1./par_IMPSOLV(1)
            pfac4 = par_IMPSOLV(9)*pfac3
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)
                dvec(k,2) = y(kk) - y(ii)
                dvec(k,3) = z(kk) - z(ii)
                termq(k) = atq(ii)*atq(kk)
                terms(k) = scrq(ii)*scrq(kk)
                tsdk(k,1:3) = scrq_dr(kk,1:3)*scrq(ii)
                tsdi(k,1:3) = scrq_dr(ii,1:3)*scrq(kk)
                termqk(k) = termq(k)*scrq(ii)
                termqi(k) = termq(k)*scrq(kk)
                termr(k) = atr(ii)+atr(kk)
              end do
            end do
          else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)
                dvec(k,2) = y(kk) - y(ii)
                dvec(k,3) = z(kk) - z(ii)
                termq(k) = atq(ii)*atq(kk)
                terms1(k) = scrq(ii)
                terms2(k) = scrq(kk)
                tsdk(k,1:3) = scrq_dr(kk,1:3)
                tsdi(k,1:3) = scrq_dr(ii,1:3)
              end do
            end do
          else ! 7,8
            pfac2 = par_IMPSOLV(8)*electric
            pfac3 = 1./par_IMPSOLV(1)
            pfac4 = par_IMPSOLV(9)*pfac3
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)
                dvec(k,2) = y(kk) - y(ii)
                dvec(k,3) = z(kk) - z(ii)
                termq(k) = atq(ii)*atq(kk)
                terms1(k) = scrq(ii)
                terms2(k) = scrq(kk)
                tsdk(k,1:3) = scrq_dr(kk,1:3)
                tsdi(k,1:3) = scrq_dr(ii,1:3)
                termr(k) = atr(ii)+atr(kk)
              end do
            end do
          end if
        else
          do i=1,at(rs1)%npol
            ii = at(rs1)%pol(i)
            do j=1,at(rs2)%npol
              kk = at(rs2)%pol(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              termq(k) = atq(ii)*atq(kk)
            end do
          end do
        end if
!       now perform the vectorizable bulk operations (maximal)
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                    dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id2(1:alcsz) = id1(1:alcsz)*id1(1:alcsz)
!       finish up -> screening model specificity is high
        if (use_IMPSOLV.EQV..true.) then
          k = 0
          if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
              call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
              call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
              tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
              tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
              tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
              tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
              tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
              tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
              termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
              termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
            end if
            term0(1:alcsz) = pfac*id1(1:alcsz)
            termqk(1:alcsz) = termqk(1:alcsz)*term0(1:alcsz)
            termqi(1:alcsz) = termqi(1:alcsz)*term0(1:alcsz)
            term0(1:alcsz) = term0(1:alcsz)*termq(1:alcsz)
            for_tmpi(1:alcsz,1) = term0(1:alcsz)*tsdi(1:alcsz,1)
            for_tmpi(1:alcsz,2) = term0(1:alcsz)*tsdi(1:alcsz,2)
            for_tmpi(1:alcsz,3) = term0(1:alcsz)*tsdi(1:alcsz,3)
            for_tmpk(1:alcsz,1) = term0(1:alcsz)*tsdk(1:alcsz,1)
            for_tmpk(1:alcsz,2) = term0(1:alcsz)*tsdk(1:alcsz,2)
            for_tmpk(1:alcsz,3) = term0(1:alcsz)*tsdk(1:alcsz,3)
            term1(1:alcsz) = term0(1:alcsz)*terms(1:alcsz)
            evec(6) = evec(6) + sum(term1(1:alcsz))
            term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
            termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
                sum_s(ii) = sum_s(ii) + termqi(k)
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
                sum_s(kk) = sum_s(kk) + termqk(k)
                ca_f(ii,1:3) = ca_f(ii,1:3) - termf(k,1:3) + for_tmpi(k,1:3)
                ca_f(kk,1:3) = ca_f(kk,1:3) + termf(k,1:3) + for_tmpk(k,1:3)
              end do
            end do
          else if (scrq_model.eq.4) then
            do i=1,at(rs1)%npol
              do j=1,at(rs2)%npol
                k = k + 1
                if (d1(k).ge.terms(k)) then
                  term1(k) = id1(k)
                  term0(k) = 2.0*id1(k)
                else
                  term1(k) = 1.0/terms(k)
                  term0(k) = term1(k)
                end if
              end do
            end do
            term2(1:alcsz) = pfac2*termq(1:alcsz)*id1(1:alcsz)
            termqi(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*term0(1:alcsz)
            term2(1:alcsz) = term2(1:alcsz)*term1(1:alcsz)
            termf(1:alcsz,1) = termqi(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = termqi(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = termqi(1:alcsz)*dvec(1:alcsz,3)
            evec(6) = evec(6) + sum(term2(1:alcsz))
            k = 0
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                ca_f(ii,1:3) = ca_f(ii,1:3) - termf(k,1:3)
                ca_f(kk,1:3) = ca_f(kk,1:3) + termf(k,1:3)
              end do
            end do
          else
            if ((scrq_model.eq.7).OR.(scrq_model.eq.8)) then
              call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
              call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
              tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
              tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
              tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
              tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
              tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
              tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
              termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
              termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
            end if
            term0(1:alcsz) = pfac*id1(1:alcsz) ! cc4 without atq-terms
            term1(1:alcsz) = electric*termq(1:alcsz)*terms(1:alcsz) ! ff1
            term2(1:alcsz) = pfac2*termq(1:alcsz)/termr(1:alcsz) ! ff2
            for_tmpi(1:alcsz,1:3) = 0.0
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                if (abs(term1(k)).gt.abs(term2(k))) then
                  term4(k) = term1(k)
                  term3(k) = 1.0
                else
                  dd1 = (d1(k) - termr(k))*pfac3
                  dd2 = 1.0 - dd1
                  if (dd1.lt.0.0) then
                    term3(k) = (1.0 - par_IMPSOLV(9))
                    term4(k) = term3(k)*term1(k) + par_IMPSOLV(9)*term2(k)
                    term0(k) = term3(k)*term0(k)
                  else if (dd2.gt.0.0) then
                    dum1 = par_IMPSOLV(9)*dd2
                    term3(k) = (1.0 - dum1)
                    term4(k) = term3(k)*term1(k) + dum1*term2(k)
                    term0(k) = term3(k)*term0(k)
                    dum1 = -scale_POLAR*id2(k)*(term2(k)-term1(k))*pfac4
                    for_tmpi(k,1) = dum1*dvec(k,1)
                    for_tmpi(k,2) = dum1*dvec(k,2)
                    for_tmpi(k,3) = dum1*dvec(k,3)
                  else
                    term4(k) = term1(k)
                    term3(k) = 1.0
                  end if
                end if
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
                sum_s(ii) = sum_s(ii) + term0(k)*termqi(k)
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
                sum_s(kk) = sum_s(kk) + term0(k)*termqk(k)
              end do
            end do
            term0(1:alcsz) = pfac*term3(1:alcsz)*termq(1:alcsz) ! overwrite of cc4
            for_tmpk(1:alcsz,1) = -for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdk(1:alcsz,1)*id1(1:alcsz)
            for_tmpk(1:alcsz,2) = -for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdk(1:alcsz,2)*id1(1:alcsz)
            for_tmpk(1:alcsz,3) = -for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdk(1:alcsz,3)*id1(1:alcsz)
            for_tmpi(1:alcsz,1) = for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdi(1:alcsz,1)*id1(1:alcsz)
            for_tmpi(1:alcsz,2) = for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdi(1:alcsz,2)*id1(1:alcsz)
            for_tmpi(1:alcsz,3) = for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdi(1:alcsz,3)*id1(1:alcsz)
            term1(1:alcsz) = scale_POLAR*term4(1:alcsz)*id1(1:alcsz)
            evec(6) = evec(6) + sum(term1(1:alcsz))
            term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
            termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
            k = 0
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                ca_f(ii,1:3) = ca_f(ii,1:3) - termf(k,1:3) + for_tmpi(k,1:3)
                ca_f(kk,1:3) = ca_f(kk,1:3) + termf(k,1:3) + for_tmpk(k,1:3)
              end do
            end do
          end if
        else
          term1(1:alcsz) = pfac*termq(1:alcsz)*id1(1:alcsz)
          evec(6) = evec(6) + sum(term1(1:alcsz))
          term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
          k = 0
          do i=1,at(rs1)%npol
            ii = at(rs1)%pol(i)
            do j=1,at(rs2)%npol
              kk = at(rs2)%pol(j)
              k = k + 1
              ca_f(ii,1:3) = ca_f(ii,1:3) - term2(k)*dvec(k,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + term2(k)*dvec(k,1:3)
            end do
          end do
        end if
      end if
!
    end if
!
    if (use_TABUL.EQV..true.) then
!
      if (tbp%rsmat(rs1,rs2).gt.0) then
        rhi = max(rs1,rs2)
        rlo = min(rs1,rs2)
        alcsz = tbp%rsmat(rhi,rlo) - tbp%rsmat(rlo,rhi) + 1
        k = 0
        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
          ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
          kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
          k = k + 1
          dvec(k,1) = x(kk) - x(ii)
          dvec(k,2) = y(kk) - y(ii)
          dvec(k,3) = z(kk) - z(ii)
        end do
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        pfac = 1.0/tbp%res
        term0(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
        tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(term0(1:alcsz))))
        term1(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
        term2(1:alcsz) = term1(1:alcsz)*term1(1:alcsz) ! quadratic
        term3(1:alcsz) = term2(1:alcsz)*term1(1:alcsz) ! cubic
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        pfac2 = pfac*scale_TABUL
!        id1(1:alcsz) = 1.0/d1(1:alcsz)
!        tbin(1:alcsz) = floor((d1(1:alcsz)-tbp%dis(1))/tbp%res) + 1
!        pfac = scale_TABUL*(1.0/tbp%res)
        k = 0
        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
          k = k + 1
          if (tbin(k).ge.tbp%bins) then
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
          else if (tbin(k).ge.1) then
            ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
            kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)

            hlp1 = 2.0*term3(k) - 3.0*term2(k)
            hlp2 = term3(k) - term2(k)
            hlp3 = hlp2 - term2(k) + term1(k)
            dhlp1 = 6.0*term2(k) - 6.0*term1(k)
            dhlp2 = 3.0*term2(k) - 2.0*term1(k)
            dhlp3 = dhlp2 - 2.0*term1(k) + 1.0

            hlp4 = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k)+1) + &
 &               tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1))

!            term1(k) = (tbp%dis(tbin(k)+1) - d1(k))*tbp%pot(tbp%lst(i,3),tbin(k))
!            term2(k) = (d1(k) - tbp%dis(tbin(k)))*tbp%pot(tbp%lst(i,3),tbin(k)+1)
!            term3(k) = pfac*(term1(k)+term2(k))
            evec(9) = evec(9) + scale_TABUL*hlp4
            hlp4 = pfac2*id1(k)*(dhlp1*(tbp%pot(tbp%lst(i,3),tbin(k)) - tbp%pot(tbp%lst(i,3),tbin(k)+1)) + &
 &               tbp%res*(dhlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + dhlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1)))
!            term4(k) = id1(k)*pfac*(tbp%pot(tbp%lst(i,3),tbin(k)+1)-tbp%pot(tbp%lst(i,3),tbin(k)))

            ca_f(ii,1:3) = ca_f(ii,1:3) + hlp4*dvec(k,1:3)
            ca_f(kk,1:3) = ca_f(kk,1:3) - hlp4*dvec(k,1:3)
         else
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),1)
          end if
        end do
      end if
    end if
!       
  end if
!
!
end
!
!--------------------------------------------------------------------------
!
! this subroutine computes all monopole-monopole and monopole-dipole, but not dipole-dipole terms
! between rs1 and rs2
!
subroutine Vforce_rsp_long_lrel(evec,rs1,rs2,svec,ca_f,sum_s)
!
  use atoms
  use polypep
  use energies
  use molecule
  use sequen
  use cutoffs
  use units
  use forces
  use iounit
!
  implicit none
!
  integer rs1,rs2,i,j,k,ii,kk,alcsz,g1,g2
  RTYPE sum_s(n)
  RTYPE evec(MAXENERGYTERMS)
  RTYPE svec(3),pfac,pfac2,pfac3,pfac4,dum1,dd1,dd2
  RTYPE d2(at(rs1)%na*at(rs2)%na),id2(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE d1(at(rs1)%na*at(rs2)%na),id1(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE termr(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  RTYPE term3(at(rs1)%na*at(rs2)%na)
  RTYPE term4(at(rs1)%na*at(rs2)%na)
  RTYPE termq(at(rs1)%na*at(rs2)%na),termqi(at(rs1)%na*at(rs2)%na),termqk(at(rs1)%na*at(rs2)%na)
  RTYPE terms(at(rs1)%na*at(rs2)%na),terms1(at(rs1)%na*at(rs2)%na),terms2(at(rs1)%na*at(rs2)%na)
  RTYPE for_tmpi(at(rs1)%na*at(rs2)%na,3)
  RTYPE for_tmpk(at(rs1)%na*at(rs2)%na,3)
  RTYPE tsdi(at(rs1)%na*at(rs2)%na,3)
  RTYPE tsdk(at(rs1)%na*at(rs2)%na,3)
  RTYPE termf(at(rs1)%na*at(rs2)%na,3)
  RTYPE ca_f(n,3)
  logical dogh
!
  if ((abs(rs2-rs1).eq.1).AND.(molofrs(rs1).eq.molofrs(rs2))) then
    write(ilog,*) 'Warning. Neighboring residues in the same molecule are beyond mid-range &
 &cutoff. This indicates an unstable simulation or a bug. Expect run to crash soon.'
  end if
!
  if (use_FEG.EQV..true.) then
    if (((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..false.)).OR.&
 &    ((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..true.).AND.(fegmode.eq.1))) then
      dogh = .false.
    else
      dogh = .true.
    end if
    if ((dogh.EQV..true.).AND.(use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
      write(ilog,*) 'Fatal. Missing support for combination of FEG and the ABSINTH&
 & implicit solvation model in Vforce_rsp_long_lrel(...). This is an omission bug.'
      call fexit()
    end if
  else
    dogh = .false.
  end if
!
  if ((use_POLAR.EQV..true.).AND.(lrel_md.eq.5)) then
!
!   first set the necessary atom-specific parameters (mimimal)
    k = 0
    pfac = electric*scale_POLAR
    alcsz = at(rs1)%npol*at(rs2)%npol
    if (alcsz.gt.0) then
      if (use_IMPSOLV.EQV..true.) then
        if (scrq_model.le.2) then
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle ! dipole-dipole
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              termq(k) = atq(ii)*atq(kk)
              terms(k) = scrq(ii)*scrq(kk)
              tsdk(k,1:3) = scrq_dr(kk,1:3)*scrq(ii)
              tsdi(k,1:3) = scrq_dr(ii,1:3)*scrq(kk)
              termqk(k) = termq(k)*scrq(ii)
              termqi(k) = termq(k)*scrq(kk)
            end do
            end do
          end do
          end do
        else if (scrq_model.eq.4) then
          pfac2 = par_IMPSOLV(8)*electric*scale_POLAR
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              termq(k) = atq(ii)*atq(kk)
              terms(k) = atr(ii)+atr(kk)
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
          pfac2 = par_IMPSOLV(8)*electric
          pfac3 = 1./par_IMPSOLV(1)
          pfac4 = par_IMPSOLV(9)*pfac3
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              termq(k) = atq(ii)*atq(kk)
              terms(k) = scrq(ii)*scrq(kk)
              tsdk(k,1:3) = scrq_dr(kk,1:3)*scrq(ii)
              tsdi(k,1:3) = scrq_dr(ii,1:3)*scrq(kk)
              termqk(k) = termq(k)*scrq(ii)
              termqi(k) = termq(k)*scrq(kk)
              termr(k) = atr(ii)+atr(kk)
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              termq(k) = atq(ii)*atq(kk)
              terms1(k) = scrq(ii)
              terms2(k) = scrq(kk)
              tsdk(k,1:3) = scrq_dr(kk,1:3)
              tsdi(k,1:3) = scrq_dr(ii,1:3)
            end do
            end do
          end do
          end do
        else ! 7,8
          pfac2 = par_IMPSOLV(8)*electric
          pfac3 = 1./par_IMPSOLV(1)
          pfac4 = par_IMPSOLV(9)*pfac3
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              termq(k) = atq(ii)*atq(kk)
              terms1(k) = scrq(ii)
              terms2(k) = scrq(kk)
              tsdk(k,1:3) = scrq_dr(kk,1:3)
              tsdi(k,1:3) = scrq_dr(ii,1:3)
              termr(k) = atr(ii)+atr(kk)
            end do
            end do
          end do
          end do
        end if
      else
        do g1=1,at(rs1)%ndpgrps
        do i=1,at(rs1)%dpgrp(g1)%nats
          ii = at(rs1)%dpgrp(g1)%ats(i)
          do g2=1,at(rs2)%ndpgrps
          if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
          do j=1,at(rs2)%dpgrp(g2)%nats
            kk = at(rs2)%dpgrp(g2)%ats(j)
            k = k + 1
            dvec(k,1) = x(kk) - x(ii)
            dvec(k,2) = y(kk) - y(ii)
            dvec(k,3) = z(kk) - z(ii)
            termq(k) = atq(ii)*atq(kk)
          end do
          end do
        end do
        end do
      end if
!     redefine alcsz
      alcsz = k
!     now perform the vectorizable bulk operations (maximal)
      dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
      dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
      dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
      d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                  dvec(1:alcsz,3)**2
      d1(1:alcsz) = sqrt(d2(1:alcsz))
      id1(1:alcsz) = 1.0/d1(1:alcsz)
      id2(1:alcsz) = id1(1:alcsz)*id1(1:alcsz)
!     finish up -> screening model specificity is high
      if (use_IMPSOLV.EQV..true.) then
        k = 0
        if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
            call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
            tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
            tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
            tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
            tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
            tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
            tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
            termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
            termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
          end if
          term0(1:alcsz) = pfac*id1(1:alcsz)
          termqk(1:alcsz) = termqk(1:alcsz)*term0(1:alcsz)
          termqi(1:alcsz) = termqi(1:alcsz)*term0(1:alcsz)
          term0(1:alcsz) = term0(1:alcsz)*termq(1:alcsz)
          for_tmpi(1:alcsz,1) = term0(1:alcsz)*tsdi(1:alcsz,1)
          for_tmpi(1:alcsz,2) = term0(1:alcsz)*tsdi(1:alcsz,2)
          for_tmpi(1:alcsz,3) = term0(1:alcsz)*tsdi(1:alcsz,3)
          for_tmpk(1:alcsz,1) = term0(1:alcsz)*tsdk(1:alcsz,1)
          for_tmpk(1:alcsz,2) = term0(1:alcsz)*tsdk(1:alcsz,2)
          for_tmpk(1:alcsz,3) = term0(1:alcsz)*tsdk(1:alcsz,3)
          term1(1:alcsz) = term0(1:alcsz)*terms(1:alcsz)
          evec(6) = evec(6) + sum(term1(1:alcsz))
          term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
          termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
          termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
          termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
              sum_s(ii) = sum_s(ii) + termqi(k)
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
              sum_s(kk) = sum_s(kk) + termqk(k)
              ca_f(ii,1:3) = ca_f(ii,1:3) - termf(k,1:3) + for_tmpi(k,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + termf(k,1:3) + for_tmpk(k,1:3)
            end do
            end do
          end do
          end do
        else if (scrq_model.eq.4) then
          do k=1,alcsz
            if (d1(k).ge.terms(k)) then
              term1(k) = id1(k)
              term0(k) = 2.0*id1(k)
            else
              term1(k) = 1.0/terms(k)
              term0(k) = term1(k)
            end if
          end do
          term2(1:alcsz) = pfac2*termq(1:alcsz)*id1(1:alcsz)
          termqi(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*term0(1:alcsz)
          term2(1:alcsz) = term2(1:alcsz)*term1(1:alcsz)
          termf(1:alcsz,1) = termqi(1:alcsz)*dvec(1:alcsz,1)
          termf(1:alcsz,2) = termqi(1:alcsz)*dvec(1:alcsz,2)
          termf(1:alcsz,3) = termqi(1:alcsz)*dvec(1:alcsz,3)
          evec(6) = evec(6) + sum(term2(1:alcsz))
          k = 0
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              ca_f(ii,1:3) = ca_f(ii,1:3) - termf(k,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + termf(k,1:3)
            end do
            end do
          end do
          end do
        else
          if ((scrq_model.eq.7).OR.(scrq_model.eq.8)) then
            call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
            call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
            tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
            tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
            tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
            tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
            tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
            tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
            termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
            termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
          end if
          term0(1:alcsz) = pfac*id1(1:alcsz) ! cc4 without atq-terms
          term1(1:alcsz) = electric*termq(1:alcsz)*terms(1:alcsz) ! ff1
          term2(1:alcsz) = pfac2*termq(1:alcsz)/termr(1:alcsz) ! ff2
          for_tmpi(1:alcsz,1:3) = 0.0
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              if (abs(term1(k)).gt.abs(term2(k))) then
                term4(k) = term1(k)
                term3(k) = 1.0
              else
                dd1 = (d1(k) - termr(k))*pfac3
                dd2 = 1.0 - dd1
                if (dd1.lt.0.0) then
                  term3(k) = (1.0 - par_IMPSOLV(9))
                  term4(k) = term3(k)*term1(k) + par_IMPSOLV(9)*term2(k)
                  term0(k) = term3(k)*term0(k)
                else if (dd2.gt.0.0) then
                  dum1 = par_IMPSOLV(9)*dd2
                  term3(k) = (1.0 - dum1)
                  term4(k) = term3(k)*term1(k) + dum1*term2(k)
                  term0(k) = term3(k)*term0(k)
                  dum1 = -scale_POLAR*id2(k)*(term2(k)-term1(k))*pfac4
                  for_tmpi(k,1) = dum1*dvec(k,1)
                  for_tmpi(k,2) = dum1*dvec(k,2)
                  for_tmpi(k,3) = dum1*dvec(k,3)
                else
                  term4(k) = term1(k)
                  term3(k) = 1.0
                end if
              end if
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
              sum_s(ii) = sum_s(ii) + term0(k)*termqi(k)
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
              sum_s(kk) = sum_s(kk) + term0(k)*termqk(k)
            end do
            end do
          end do
          end do
          term0(1:alcsz) = pfac*term3(1:alcsz)*termq(1:alcsz) ! overwrite of cc4
          for_tmpk(1:alcsz,1) = -for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdk(1:alcsz,1)*id1(1:alcsz)
          for_tmpk(1:alcsz,2) = -for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdk(1:alcsz,2)*id1(1:alcsz)
          for_tmpk(1:alcsz,3) = -for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdk(1:alcsz,3)*id1(1:alcsz)
          for_tmpi(1:alcsz,1) = for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdi(1:alcsz,1)*id1(1:alcsz)
          for_tmpi(1:alcsz,2) = for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdi(1:alcsz,2)*id1(1:alcsz)
          for_tmpi(1:alcsz,3) = for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdi(1:alcsz,3)*id1(1:alcsz)
          term1(1:alcsz) = scale_POLAR*term4(1:alcsz)*id1(1:alcsz)
          evec(6) = evec(6) + sum(term1(1:alcsz))
          term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
          termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
          termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
          termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
          k = 0
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              ca_f(ii,1:3) = ca_f(ii,1:3) - termf(k,1:3) + for_tmpi(k,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + termf(k,1:3) + for_tmpk(k,1:3)
            end do
            end do
          end do
          end do
        end if
      else
        if (dogh.EQV..true.) then
          term1(1:alcsz) = electric*par_FEG2(9)*termq(1:alcsz)/(par_FEG2(10) + d1(1:alcsz))
          evec(6) = evec(6) + sum(term1(1:alcsz))
          term2(1:alcsz) = term1(1:alcsz)*id1(1:alcsz)/(par_FEG2(10) + d1(1:alcsz))
        else
          term1(1:alcsz) = pfac*termq(1:alcsz)*id1(1:alcsz)
          evec(6) = evec(6) + sum(term1(1:alcsz))
          term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
        end if
        k = 0
        do g1=1,at(rs1)%ndpgrps
        do i=1,at(rs1)%dpgrp(g1)%nats
          ii = at(rs1)%dpgrp(g1)%ats(i)
          do g2=1,at(rs2)%ndpgrps
          if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
          do j=1,at(rs2)%dpgrp(g2)%nats
            kk = at(rs2)%dpgrp(g2)%ats(j)
            k = k + 1
            ca_f(ii,1:3) = ca_f(ii,1:3) - term2(k)*dvec(k,1:3)
            ca_f(kk,1:3) = ca_f(kk,1:3) + term2(k)*dvec(k,1:3)
          end do
          end do
        end do
        end do
      end if
    end if
!
  else if ((use_POLAR.EQV..true.).AND.(lrel_md.eq.4)) then
!
  end if
!
end
!
!---------------------------------------------------------------------------------------------
!
subroutine Vgenmu(r1,r2,i_genmu,ro,alcsz)
!
  implicit none
!
  integer i_genmu,alcsz
  RTYPE r1(alcsz),r2(alcsz),ro(alcsz),dum1
!
  if (i_genmu.eq.0) then
    ro(1:alcsz) = sqrt(r1(1:alcsz)*r2(1:alcsz))
  else if (i_genmu.eq.1) then
    ro(1:alcsz) = 0.5*(r1(1:alcsz) + r2(1:alcsz))
  else if (i_genmu.eq.-1) then
    ro(1:alcsz) = 1.0/(0.5*(1./r1(1:alcsz) + 1./r2(1:alcsz)))
  else if (i_genmu.eq.2) then
    ro(1:alcsz) = sqrt(0.5*(r1(1:alcsz)*r1(1:alcsz) + r2(1:alcsz)*r2(1:alcsz)))
  else
    dum1 = 1./(1.*i_genmu)
    ro(1:alcsz) = (0.5*(r1(1:alcsz)**(i_genmu) + r2(1:alcsz)**(i_genmu)))**(dum1)
  end if
!
end
!
!---------------------------------------------------------------------------
!
subroutine Vgenmu_dr(r1,r2,i_genmu,dr1,dr2,alcsz)
!
  implicit none
!
  integer i_genmu,alcsz
  RTYPE dr1(alcsz),dr2(alcsz),r1(alcsz),r2(alcsz),help(alcsz),dum1,dum3
!
  if (i_genmu.eq.0) then
    help(1:alcsz) = 1.0/(2.0*sqrt(r1(1:alcsz)*r2(1:alcsz)))
    dr1(1:alcsz) = help(1:alcsz)*r2(1:alcsz)
    dr2(1:alcsz) = dr1(1:alcsz)*r1(1:alcsz)/r2(1:alcsz)
  else if (i_genmu.eq.1) then
    dr1(1:alcsz) = 0.5
    dr2(1:alcsz) = dr1(1:alcsz)
  else if (i_genmu.eq.-1) then
    help(1:alcsz) = 1./r1(1:alcsz) + 1./r2(1:alcsz)
    dr1(1:alcsz) = 1./(((r1(1:alcsz)*r1(1:alcsz))*help(1:alcsz))*(0.5*help(1:alcsz)))
    dr2(1:alcsz) = 1./(((r2(1:alcsz)*r2(1:alcsz))*help(1:alcsz))*(0.5*help(1:alcsz)))
  else if (i_genmu.eq.2) then
    help(1:alcsz) = 1.0/(2.0*sqrt(0.5*(r1(1:alcsz)*r1(1:alcsz) + r2(1:alcsz)*r2(1:alcsz))))
    dr1(1:alcsz) = r1(1:alcsz)*help(1:alcsz)
    dr2(1:alcsz) = r2(1:alcsz)*help(1:alcsz) ! dr1(1:alcsz)/r1(1:alcsz)
  else
    help(1:alcsz) = 0.5*(r1(1:alcsz)**i_genmu + r2(1:alcsz)**i_genmu)
    dum1 = 1./(1.*i_genmu) - 1.0
    dum3 = 1.*i_genmu - 1.0
    dr1(1:alcsz) = 0.5*((help(1:alcsz))**dum1)*((r1(1:alcsz)**dum3))
    dr2(1:alcsz) = 0.5*((help(1:alcsz))**dum1)*((r2(1:alcsz)**dum3))
  end if
!
end
!
!-----------------------------------------------------------------------
!
! the vectorized version of force_rsp_long_feg
! this does not currently support IMPSOLV (note Ven_rsp_long does)
! support is very easy to add for screening models 1,2,5,6, but more annoying for the
! other ones
!
subroutine Vforce_rsp_long_feg(evec,rs1,rs2,cut,ca_f,sum_s)
!
  use iounit
  use energies
  use atoms
  use polypep
  use inter
  use units
  use tabpot
  use molecule
  use sequen
  use cutoffs
  use system
  use grandensembles
  use forces
!
  implicit none
!
  integer rs,rs1,rs2,i,j,k,ii,kk,alcsz
  RTYPE sum_s(n)
  RTYPE evec(MAXENERGYTERMS)
  RTYPE svec(3),pfac
  RTYPE d2(at(rs1)%na*at(rs2)%na),id1s(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE d1(at(rs1)%na*at(rs2)%na),id1(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  RTYPE termq(at(rs1)%na*at(rs2)%na)
  logical cut
  RTYPE ca_f(n,3)
!
! we have a potential override to cover in par_FEG3, which allows residues to be
! fully de-coupled at all times
  if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) return
!
  svec(:) = 0.0
!
!
! three different cases: intra-residue (rs1 == rs2), neighbors in seq., or others
!
  if (rs1.eq.rs2) then
    rs = rs1
    if (use_FEGS(6).EQV..true.) then
!
!     first set the necessary atom-specific parameters (mimimal)
      k = 0
      pfac = electric*par_FEG2(9)
      alcsz = nrpolintra(rs)
      if (alcsz.gt.0) then
        do i=1,nrpolintra(rs)
          ii = iaa(rs)%polin(i,1)
          kk = iaa(rs)%polin(i,2)
          dvec(i,1) = x(kk) - x(ii)
          dvec(i,2) = y(kk) - y(ii)
          dvec(i,3) = z(kk) - z(ii)
          termq(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
        end do
!       now perform the vectorizable bulk operations (maximal)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                    dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id1s(1:alcsz) = 1.0/(d1(1:alcsz) + par_FEG2(10))
        term1(1:alcsz) = pfac*termq(1:alcsz)*id1s(1:alcsz)
        evec(6) = evec(6) + sum(term1(1:alcsz))
        term2(1:alcsz) = term1(1:alcsz)*id1s(1:alcsz)*id1(1:alcsz)
        do i=1,nrpolintra(rs)
          ii = iaa(rs)%polin(i,1)
          kk = iaa(rs)%polin(i,2)
          ca_f(ii,1:3) = ca_f(ii,1:3) - term2(i)*dvec(i,1:3)
          ca_f(kk,1:3) = ca_f(kk,1:3) + term2(i)*dvec(i,1:3)
        end do
      end if
!
    end if
!
!
!
! neighboring residues
!
  else if (abs(rs1-rs2).eq.1) then
!
!   there is no topological requirement for neighboring residues to be close (pure sequence)
!   so we have to check for BC evtl.y
!
    if (rs1.gt.rs2) then
      rs = rs2
      call dis_bound_rs(rs2,rs1,svec)
    else
      rs = rs1
      call dis_bound_rs(rs1,rs2,svec)
    end if
!
    if (use_FEGS(6).EQV..true.) then
!
!     first set the necessary atom-specific parameters (mimimal)
      k = 0
      pfac = electric*par_FEG2(9)
      alcsz = nrpolnb(rs)
      if (alcsz.gt.0) then
        do i=1,nrpolnb(rs)
          ii = iaa(rs)%polnb(i,1)
          kk = iaa(rs)%polnb(i,2)
          dvec(i,1) = x(kk) - x(ii)
          dvec(i,2) = y(kk) - y(ii)
          dvec(i,3) = z(kk) - z(ii)
          termq(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
        end do
!       now perform the vectorizable bulk operations (maximal)
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                    dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id1s(1:alcsz) = 1.0/(d1(1:alcsz) + par_FEG2(10))
        term1(1:alcsz) = pfac*termq(1:alcsz)*id1s(1:alcsz)
        evec(6) = evec(6) + sum(term1(1:alcsz))
        term2(1:alcsz) = term1(1:alcsz)*id1s(1:alcsz)*id1(1:alcsz)
        do i=1,nrpolnb(rs)
          ii = iaa(rs)%polnb(i,1)
          kk = iaa(rs)%polnb(i,2)
          ca_f(ii,1:3) = ca_f(ii,1:3) - term2(i)*dvec(i,1:3)
          ca_f(kk,1:3) = ca_f(kk,1:3) + term2(i)*dvec(i,1:3)
        end do
      end if
!
    end if
!
!
!
! all other residue pairs
!
 67   format(4g14.6)
  else
!
!   there is no topological relationship for remaining residues -> always check BC
    call dis_bound_rs(rs1,rs2,svec)
!
    if (use_FEGS(6).EQV..true.) then
!
!     first set the necessary atom-specific parameters (mimimal)
      k = 0
      pfac = electric*par_FEG2(9)
      alcsz = at(rs1)%npol*at(rs2)%npol
      if (alcsz.gt.0) then
        do i=1,at(rs1)%npol
          ii = at(rs1)%pol(i)
          do j=1,at(rs2)%npol
            kk = at(rs2)%pol(j)
            k = k + 1
            dvec(k,1) = x(kk) - x(ii)
            dvec(k,2) = y(kk) - y(ii)
            dvec(k,3) = z(kk) - z(ii)
            termq(k) = atq(ii)*atq(kk)
          end do
        end do
!       now perform the vectorizable bulk operations (maximal)
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                    dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id1s(1:alcsz) = 1.0/(d1(1:alcsz) + par_FEG2(10))
        term1(1:alcsz) = pfac*termq(1:alcsz)*id1s(1:alcsz)
        evec(6) = evec(6) + sum(term1(1:alcsz))
        term2(1:alcsz) = term1(1:alcsz)*id1s(1:alcsz)*id1(1:alcsz)
        k = 0
        do i=1,at(rs1)%npol
          ii = at(rs1)%pol(i)
          do j=1,at(rs2)%npol
            kk = at(rs2)%pol(j)
            k = k + 1
            ca_f(ii,1:3) = ca_f(ii,1:3) - term2(k)*dvec(k,1:3)
            ca_f(kk,1:3) = ca_f(kk,1:3) + term2(k)*dvec(k,1:3)
          end do
        end do
      end if
!
    end if
!
!       
  end if
!
!
end
!
!-----------------------------------------------------------------------
!
subroutine Vforce_PSCRM12_C(evec,rs1,ca_f,sum_s)
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
  integer rs1,rs2,ii,j,hi,lo,k,i
  RTYPE en_incr
  RTYPE evec(MAXENERGYTERMS)
  RTYPE termc(rs_nbl(rs1)%nnbats)
  RTYPE term0(rs_nbl(rs1)%nnbats)
  RTYPE termq(rs_nbl(rs1)%nnbats)
  RTYPE termsq(rs_nbl(rs1)%nnbats)
  RTYPE termdii(rs_nbl(rs1)%nnbats,3)
  RTYPE termdkk(rs_nbl(rs1)%nnbats,3)
  RTYPE termsqdr(rs_nbl(rs1)%nnbats,3)
  RTYPE termcii(rs_nbl(rs1)%nnbats)
  RTYPE termckk(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE id1(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  integer sh(rs_nbl(rs1)%nnbs),hira,lora
  integer k1(rs_nbl(rs1)%nnbs,3),k2(rs_nbl(rs1)%nnbs,3)
  RTYPE drav(rs_nbl(rs1)%nnbs,3),srav(rs_nbl(rs1)%nnbs,3)
  RTYPE ca_f(n,3),sum_s(n)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
  termckk(:) = 0.0
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
    termq(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsq(lo:hi) = scrq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsqdr(lo:hi,1) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1)
    termsqdr(lo:hi,2) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2)
    termsqdr(lo:hi,3) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3)
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    if (atq(ii).eq.0.0) cycle
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    id1(:) = sqrt(id2(:))
!   screened Coulomb potential
    term0(:) = electric*termq(:)*atq(ii)
    termdkk(:,1) = term0(:)*scrq(ii)*termsqdr(:,1)
    termdkk(:,2) = term0(:)*scrq(ii)*termsqdr(:,2)
    termdkk(:,3) = term0(:)*scrq(ii)*termsqdr(:,3)
    termdii(:,1) = term0(:)*termsq(:)*scrq_dr(ii,1)
    termdii(:,2) = term0(:)*termsq(:)*scrq_dr(ii,2)
    termdii(:,3) = term0(:)*termsq(:)*scrq_dr(ii,3)
    termc(:) = term0(:)*scale_POLAR*id1(:)
    termcii(:) = termc(:)*termsq(:)
    termckk(:) = termckk(:) + termc(:)*scrq(ii)
    en_incr = sum(termcii(:))
    sum_s(ii) = sum_s(ii) + en_incr
    evec(6) = evec(6) + scrq(ii)*en_incr
!
    foin(:,1) = dvec(:,1)*scrq(ii)*termcii(:)*id2(:)
    ca_f(ii,1) = ca_f(ii,1) - (sum(foin(:,1)) - scale_POLAR*sum(termdii(:,1)*id1(:)))
    for_k(:,1) = for_k(:,1) + foin(:,1) + scale_POLAR*termdkk(:,1)*id1(:)
    foin(:,2) = dvec(:,2)*scrq(ii)*termcii(:)*id2(:)
    ca_f(ii,2) = ca_f(ii,2) - (sum(foin(:,2)) - scale_POLAR*sum(termdii(:,2)*id1(:)))
    for_k(:,2) = for_k(:,2) + foin(:,2) + scale_POLAR*termdkk(:,2)*id1(:)
    foin(:,3) = dvec(:,3)*scrq(ii)*termcii(:)*id2(:)
    ca_f(ii,3) = ca_f(ii,3) - (sum(foin(:,3)) - scale_POLAR*sum(termdii(:,3)*id1(:)))
    for_k(:,3) = for_k(:,3) + foin(:,3) + scale_POLAR*termdkk(:,3)*id1(:)
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
    sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) = &
 &     sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + termckk(lo:hi)
    lo = hi + 1
  end do
!
end
!
!------------------------------------------------------------------------------------
!
subroutine Vforce_PSCRM12_TR(evec,rs1,ca_f,sum_s)
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
  integer rs1,rs2,ii,hi,lo,k,i
  RTYPE en_incr
  RTYPE evec(MAXENERGYTERMS)
  RTYPE termc(rs_nbl(rs1)%nnbtrats)
  RTYPE term0(rs_nbl(rs1)%nnbtrats)
  RTYPE termq(rs_nbl(rs1)%nnbtrats)
  RTYPE termsq(rs_nbl(rs1)%nnbtrats)
  RTYPE termdii(rs_nbl(rs1)%nnbtrats,3)
  RTYPE termdkk(rs_nbl(rs1)%nnbtrats,3)
  RTYPE termsqdr(rs_nbl(rs1)%nnbtrats,3)
  RTYPE termcii(rs_nbl(rs1)%nnbtrats)
  RTYPE termckk(rs_nbl(rs1)%nnbtrats)
  RTYPE foin(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbtrats,3)
  RTYPE id2(rs_nbl(rs1)%nnbtrats)
  RTYPE dis2(rs_nbl(rs1)%nnbtrats)
  RTYPE id1(rs_nbl(rs1)%nnbtrats)
  RTYPE for_k(rs_nbl(rs1)%nnbtrats,3)
  integer sh(rs_nbl(rs1)%nnbtrs)
  RTYPE ca_f(n,3),sum_s(n)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
  termckk(:) = 0.0
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nbb + at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(k)
    hi = lo + sh(k)
    dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,3)
    termq(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsq(lo:hi) = scrq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsqdr(lo:hi,1) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1)
    termsqdr(lo:hi,2) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2)
    termsqdr(lo:hi,3) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3)
    lo = hi + 1
  end do
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    if (atq(ii).eq.0.0) cycle
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    id1(:) = sqrt(id2(:))
!   screened Coulomb potential
    term0(:) = electric*termq(:)*atq(ii)
    termdkk(:,1) = term0(:)*scrq(ii)*termsqdr(:,1)
    termdkk(:,2) = term0(:)*scrq(ii)*termsqdr(:,2)
    termdkk(:,3) = term0(:)*scrq(ii)*termsqdr(:,3)
    termdii(:,1) = term0(:)*termsq(:)*scrq_dr(ii,1)
    termdii(:,2) = term0(:)*termsq(:)*scrq_dr(ii,2)
    termdii(:,3) = term0(:)*termsq(:)*scrq_dr(ii,3)
    termc(:) = term0(:)*scale_POLAR*id1(:)
    termcii(:) = termc(:)*termsq(:)
    termckk(:) = termckk(:) + termc(:)*scrq(ii)
    en_incr = sum(termcii(:))
    sum_s(ii) = sum_s(ii) + en_incr
    evec(6) = evec(6) + scrq(ii)*en_incr
!
    foin(:,1) = dvec(:,1)*scrq(ii)*termcii(:)*id2(:)
    ca_f(ii,1) = ca_f(ii,1) - (sum(foin(:,1)) - scale_POLAR*sum(termdii(:,1)*id1(:)))
    for_k(:,1) = for_k(:,1) + foin(:,1) + scale_POLAR*termdkk(:,1)*id1(:)
    foin(:,2) = dvec(:,2)*scrq(ii)*termcii(:)*id2(:)
    ca_f(ii,2) = ca_f(ii,2) - (sum(foin(:,2)) - scale_POLAR*sum(termdii(:,2)*id1(:)))
    for_k(:,2) = for_k(:,2) + foin(:,2) + scale_POLAR*termdkk(:,2)*id1(:)
    foin(:,3) = dvec(:,3)*scrq(ii)*termcii(:)*id2(:)
    ca_f(ii,3) = ca_f(ii,3) - (sum(foin(:,3)) - scale_POLAR*sum(termdii(:,3)*id1(:)))
    for_k(:,3) = for_k(:,3) + foin(:,3) + scale_POLAR*termdkk(:,3)*id1(:)
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
    sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) = &
 &     sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + termckk(lo:hi)
    lo = hi + 1
  end do
!
end
!
!-------------------------------------------------------------------
!
subroutine Vforce_PSCRM4_C(evec,rs1,ca_f)
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
  integer rs1,rs2,ii,j,hi,lo,k,i
  RTYPE pfac2,dum1,dum3
  RTYPE evec(MAXENERGYTERMS)
  RTYPE termc(rs_nbl(rs1)%nnbats)
  RTYPE term0(rs_nbl(rs1)%nnbats)
  RTYPE term1(rs_nbl(rs1)%nnbats)
  RTYPE termq(rs_nbl(rs1)%nnbats)
  RTYPE termr(rs_nbl(rs1)%nnbats)
  RTYPE termh(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE id1(rs_nbl(rs1)%nnbats)
  RTYPE d1(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  integer sh(rs_nbl(rs1)%nnbs),hira,lora
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
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))%nbb + at(rs_nbl(rs1)%nb(1:rs_nbl(rs1)%nnbs))%nsc - 1
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
    termq(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termr(lo:hi) = atr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    lo = hi + 1
  end do
  pfac2 = par_IMPSOLV(8)*electric*scale_POLAR
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    if (atq(ii).eq.0.0) cycle
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    d1(:) = sqrt(dis2(:))
    id1(:) = 1.0/d1(:)
    id2(:) = id1(:)*id1(:)
    lo = 1
    do k=1,rs_nbl(rs1)%nnbs
      rs2 = rs_nbl(rs1)%nb(k)
      hi = lo + sh(k)
      do j=lo,hi
        dum1 = atr(ii) + termr(j)
        if (d1(j).ge.dum1) then
          term1(j) = id1(j)
          term0(j) = 2.0*id1(j)
        else
          term1(j) = 1.0/dum1
          term0(j) = term1(j)
        end if
      end do
      lo = hi + 1
    end do
!
    dum3 = pfac2*atq(ii)
    termc(:) = dum3*termq(:)*id1(:)
    termh(:) = termc(:)*id2(:)*term0(:)
    termc(:) = termc(:)*term1(:)
    evec(6) = evec(6) + sum(termc(:))
!
    foin(:,1) = dvec(:,1)*termh(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*termh(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*termh(:)
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
!-------------------------------------------------------------------
!
subroutine Vforce_PSCRM4_TR(evec,rs1,ca_f)
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
  integer rs1,rs2,ii,hi,lo,k,i,j
  RTYPE pfac2,dum1,dum3
  RTYPE evec(MAXENERGYTERMS)
  RTYPE termc(rs_nbl(rs1)%nnbtrats)
  RTYPE term0(rs_nbl(rs1)%nnbtrats)
  RTYPE term1(rs_nbl(rs1)%nnbtrats)
  RTYPE termq(rs_nbl(rs1)%nnbtrats)
  RTYPE termr(rs_nbl(rs1)%nnbtrats)
  RTYPE termh(rs_nbl(rs1)%nnbtrats)
  RTYPE foin(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbtrats,3)
  RTYPE id2(rs_nbl(rs1)%nnbtrats)
  RTYPE dis2(rs_nbl(rs1)%nnbtrats)
  RTYPE id1(rs_nbl(rs1)%nnbtrats)
  RTYPE d1(rs_nbl(rs1)%nnbtrats)
  RTYPE for_k(rs_nbl(rs1)%nnbtrats,3)
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
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(k)
    hi = lo + sh(k)
    dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,3)
    termq(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termr(lo:hi) = atr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    lo = hi + 1
  end do
  pfac2 = par_IMPSOLV(8)*electric*scale_POLAR
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    if (atq(ii).eq.0.0) cycle
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    d1(:) = sqrt(dis2(:))
    id1(:) = 1.0/d1(:)
    id2(:) = id1(:)*id1(:)
    lo = 1
    do k=1,rs_nbl(rs1)%nnbtrs
      rs2 = rs_nbl(rs1)%nbtr(k)
      hi = lo + sh(k)
      do j=lo,hi
        dum1 = atr(ii) + termr(j)
        if (d1(j).ge.dum1) then
          term1(j) = id1(j)
          term0(j) = 2.0*id1(j)
        else
          term1(j) = 1.0/dum1
          term0(j) = term1(j)
        end if
      end do
      lo = hi + 1
    end do
!
    dum3 = pfac2*atq(ii)
    termc(:) = dum3*termq(:)*id1(:)
    termh(:) = termc(:)*id2(:)*term0(:)
    termc(:) = termc(:)*term1(:)
    evec(6) = evec(6) + sum(termc(:))
!
    foin(:,1) = dvec(:,1)*termh(:)
    ca_f(ii,1) = ca_f(ii,1) - sum(foin(:,1))
    for_k(:,1) = for_k(:,1) + foin(:,1)
    foin(:,2) = dvec(:,2)*termh(:)
    ca_f(ii,2) = ca_f(ii,2) - sum(foin(:,2))
    for_k(:,2) = for_k(:,2) + foin(:,2)
    foin(:,3) = dvec(:,3)*termh(:)
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
!---------------------------------------------------------------------------
!
subroutine Vgenmu_all2(r1,r2,i_genmu,ro,dr1,dr2,alcsz)
!
  use iounit
!
  implicit none
!
  integer i_genmu,alcsz,i
  RTYPE dr1(alcsz),dr2(alcsz),r1,r2(alcsz),help(alcsz),dum1,dum2,dum3,ro(alcsz)
  RTYPE help2(alcsz)
!
  if ((r1.eq.0.0).AND.(i_genmu.ne.1)) then
    write(ilog,*) 'Fatal. Called Vgenmu_all2 with corrupted arguments.'
    call fexit()
  end if
!
  if (i_genmu.eq.0) then
    ro(1:alcsz) = sqrt(r1*r2(1:alcsz))
    do i=1,alcsz
      if (r2(i).eq.0.0) then
        dr1(i) = 0.0
        dr2(i) = 0.0
      else
        help(i) = 1.0/(2.0*ro(i))
        dr1(i) = help(i)*r2(i)
        dr2(i) = r1*dr1(i)/r2(i)
      end if
    end do
  else if (i_genmu.eq.1) then
    ro(1:alcsz) = 0.5*(r1 + r2(1:alcsz))
    dr1(1:alcsz) = 0.5
    dr2(1:alcsz) = dr1(1:alcsz)
  else if (i_genmu.eq.-1) then
    dum1 = 1./r1
    do i=1,alcsz
      if (r2(i).eq.0.0) then
        help(i) = 0.02
      else
        help(i) = dum1 + 1./r2(i)
      end if
    end do
    ro(1:alcsz) = 2.0/(help(1:alcsz))
    dum1 = r1*r1*0.5
    help2(1:alcsz) = help(1:alcsz)*help(1:alcsz)
    dr1(1:alcsz) = 1./(dum1*help2(1:alcsz))
    help2(1:alcsz) = (help2(1:alcsz)*r2(1:alcsz)*r2(1:alcsz))
    do i=1,alcsz
      if (r2(i).eq.0.0) then
        dr2(i) = 0.0
        dr1(i) = 0.0
        ro(i) = 0.0
      else
        dr2(i) = 2.0/help2(i)
      end if
    end do
  else if (i_genmu.eq.2) then
    dum1 = r1*r1
    ro(1:alcsz) = sqrt(0.5*(dum1 + r2(1:alcsz)*r2(1:alcsz)))
    dum2 = r1/2.0
    dr1(1:alcsz) = dum2/ro(1:alcsz)
    dum3 = 1./r1
    dr2(1:alcsz) = dum3*r2(1:alcsz)*dr1(1:alcsz)
  else
    dum3 = r1**i_genmu
    do i=1,alcsz
      if (r2(i).eq.0.0) then
        help(i) = 0.02
      else
        help(i) = 0.5*(dum3 + r2(i)**i_genmu)
      end if
    end do
    dum2 = 1./(1.*i_genmu)
    ro(1:alcsz) = help(1:alcsz)**(dum2)
    dum2 = dum2 - 1.0
    help2(1:alcsz) = (help(1:alcsz))**(dum2)
    dum1 = 1.*i_genmu - 1.0
    dr1(1:alcsz) = 0.5*help2(1:alcsz)*(r1**dum1)
    do i=1,alcsz
      if (r2(i).eq.0.0) then
        dr2(i) = 0.0
        dr1(i) = 0.0
        ro(i) = 0.0
      else
        dr2(i) = 0.5*help2(i)*(r2(i)**dum1)
      end if
    end do
  end if
!
end
!
!-------------------------------------------------------------------
!
subroutine Vforce_PSCRM56_C(evec,rs1,ca_f,sum_s)
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
  integer rs1,rs2,ii,j,hi,lo,k,i
  RTYPE en_incr,pfac
  RTYPE evec(MAXENERGYTERMS)
  RTYPE termc(rs_nbl(rs1)%nnbats)
  RTYPE term0(rs_nbl(rs1)%nnbats)
  RTYPE term1(rs_nbl(rs1)%nnbats)
  RTYPE term2(rs_nbl(rs1)%nnbats)
  RTYPE termh(rs_nbl(rs1)%nnbats)
  RTYPE termq(rs_nbl(rs1)%nnbats)
  RTYPE terms(rs_nbl(rs1)%nnbats)
  RTYPE termsq(rs_nbl(rs1)%nnbats)
  RTYPE termdii(rs_nbl(rs1)%nnbats,3)
  RTYPE termdkk(rs_nbl(rs1)%nnbats,3)
  RTYPE termsqdr(rs_nbl(rs1)%nnbats,3)
  RTYPE termcii(rs_nbl(rs1)%nnbats)
  RTYPE termckk(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE id1(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  integer sh(rs_nbl(rs1)%nnbs),hira,lora
  integer k1(rs_nbl(rs1)%nnbs,3),k2(rs_nbl(rs1)%nnbs,3)
  RTYPE drav(rs_nbl(rs1)%nnbs,3),srav(rs_nbl(rs1)%nnbs,3)
  RTYPE ca_f(n,3),sum_s(n)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
  termckk(:) = 0.0
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
    termq(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsq(lo:hi) = scrq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsqdr(lo:hi,1) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1)
    termsqdr(lo:hi,2) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2)
    termsqdr(lo:hi,3) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3)
    lo = hi + 1
  end do
  pfac = scale_POLAR*electric
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    if (atq(ii).eq.0.0) cycle
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    id1(:) = sqrt(id2(:))
!   screened Coulomb potential
    call Vgenmu_all2(scrq(ii),termsq,i_sqm,terms,term1,term2,rs_nbl(rs1)%nnbats)
    term0(:) = pfac*termq(:)*atq(ii)
    termh(:) = term0(:)*term2(:)
    termdkk(:,1) = termh(:)*termsqdr(:,1)
    termdkk(:,2) = termh(:)*termsqdr(:,2)
    termdkk(:,3) = termh(:)*termsqdr(:,3)
    termh(:) = term0(:)*term1(:)
    termdii(:,1) = termh(:)*scrq_dr(ii,1)
    termdii(:,2) = termh(:)*scrq_dr(ii,2)
    termdii(:,3) = termh(:)*scrq_dr(ii,3)
    termc(:) = term0(:)*id1(:)
    termcii(:) = termc(:)*term1(:)
    termckk(:) = termckk(:) + termc(:)*term2(:)
!   overwrite
    term0(:) = terms(:)*termc(:)
    en_incr = sum(term0(:))
    sum_s(ii) = sum_s(ii) + sum(termcii(:))
    evec(6) = evec(6) + en_incr
!
    foin(:,1) = dvec(:,1)*term0(:)*id2(:)
    ca_f(ii,1) = ca_f(ii,1) - (sum(foin(:,1)) - sum(termdii(:,1)*id1(:)))
    for_k(:,1) = for_k(:,1) + foin(:,1) + termdkk(:,1)*id1(:)
    foin(:,2) = dvec(:,2)*term0(:)*id2(:)
    ca_f(ii,2) = ca_f(ii,2) - (sum(foin(:,2)) - sum(termdii(:,2)*id1(:)))
    for_k(:,2) = for_k(:,2) + foin(:,2) + termdkk(:,2)*id1(:)
    foin(:,3) = dvec(:,3)*term0(:)*id2(:)
    ca_f(ii,3) = ca_f(ii,3) - (sum(foin(:,3)) - sum(termdii(:,3)*id1(:)))
    for_k(:,3) = for_k(:,3) + foin(:,3) + termdkk(:,3)*id1(:)
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
    sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) = &
 &     sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + termckk(lo:hi)
    lo = hi + 1
  end do
!
end
!
!-------------------------------------------------------------------
!
subroutine Vforce_PSCRM56_TR(evec,rs1,ca_f,sum_s)
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
  integer rs1,rs2,ii,hi,lo,k,i
  RTYPE en_incr,pfac
  RTYPE evec(MAXENERGYTERMS)
  RTYPE termc(rs_nbl(rs1)%nnbtrats)
  RTYPE term0(rs_nbl(rs1)%nnbtrats)
  RTYPE term1(rs_nbl(rs1)%nnbtrats)
  RTYPE term2(rs_nbl(rs1)%nnbtrats)
  RTYPE termh(rs_nbl(rs1)%nnbtrats)
  RTYPE termq(rs_nbl(rs1)%nnbtrats)
  RTYPE terms(rs_nbl(rs1)%nnbtrats)
  RTYPE termsq(rs_nbl(rs1)%nnbtrats)
  RTYPE termdii(rs_nbl(rs1)%nnbtrats,3)
  RTYPE termdkk(rs_nbl(rs1)%nnbtrats,3)
  RTYPE termsqdr(rs_nbl(rs1)%nnbtrats,3)
  RTYPE termcii(rs_nbl(rs1)%nnbtrats)
  RTYPE termckk(rs_nbl(rs1)%nnbtrats)
  RTYPE foin(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbtrats,3)
  RTYPE id2(rs_nbl(rs1)%nnbtrats)
  RTYPE dis2(rs_nbl(rs1)%nnbtrats)
  RTYPE id1(rs_nbl(rs1)%nnbtrats)
  RTYPE for_k(rs_nbl(rs1)%nnbtrats,3)
  integer sh(rs_nbl(rs1)%nnbtrs)
  RTYPE ca_f(n,3),sum_s(n)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
  termckk(:) = 0.0
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nbb + at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(k)
    hi = lo + sh(k)
    dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,3)
    termq(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsq(lo:hi) = scrq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsqdr(lo:hi,1) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1)
    termsqdr(lo:hi,2) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2)
    termsqdr(lo:hi,3) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3)
    lo = hi + 1
  end do
  pfac = scale_POLAR*electric
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    if (atq(ii).eq.0.0) cycle
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    id2(:) = 1.0/dis2(:)
    id1(:) = sqrt(id2(:))
!   screened Coulomb potential
    call Vgenmu_all2(scrq(ii),termsq,i_sqm,terms,term1,term2,rs_nbl(rs1)%nnbtrats)
    term0(:) = pfac*termq(:)*atq(ii)
    termh(:) = term0(:)*term2(:)
    termdkk(:,1) = termh(:)*termsqdr(:,1)
    termdkk(:,2) = termh(:)*termsqdr(:,2)
    termdkk(:,3) = termh(:)*termsqdr(:,3)
    termh(:) = term0(:)*term1(:)
    termdii(:,1) = termh(:)*scrq_dr(ii,1)
    termdii(:,2) = termh(:)*scrq_dr(ii,2)
    termdii(:,3) = termh(:)*scrq_dr(ii,3)
    termc(:) = term0(:)*id1(:)
    termcii(:) = termc(:)*term1(:)
    termckk(:) = termckk(:) + termc(:)*term2(:)
!   overwrite
    term0(:) = terms(:)*termc(:)
    en_incr = sum(term0(:))
    sum_s(ii) = sum_s(ii) + sum(termcii(:))
    evec(6) = evec(6) + en_incr
!
    foin(:,1) = dvec(:,1)*term0(:)*id2(:)
    ca_f(ii,1) = ca_f(ii,1) - (sum(foin(:,1)) - sum(termdii(:,1)*id1(:)))
    for_k(:,1) = for_k(:,1) + foin(:,1) + termdkk(:,1)*id1(:)
    foin(:,2) = dvec(:,2)*term0(:)*id2(:)
    ca_f(ii,2) = ca_f(ii,2) - (sum(foin(:,2)) - sum(termdii(:,2)*id1(:)))
    for_k(:,2) = for_k(:,2) + foin(:,2) + termdkk(:,2)*id1(:)
    foin(:,3) = dvec(:,3)*term0(:)*id2(:)
    ca_f(ii,3) = ca_f(ii,3) - (sum(foin(:,3)) - sum(termdii(:,3)*id1(:)))
    for_k(:,3) = for_k(:,3) + foin(:,3) + termdkk(:,3)*id1(:)
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
    sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) = &
 &     sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + termckk(lo:hi)
    lo = hi + 1
  end do
!
end
!
!------------------------------------------------------------------------------------
!
! for the mixed models, the fxn gets pretty messy, and hence provides a minor boost
! relative to Vforce_rsp_long(...)
!
subroutine Vforce_PSCRM78_C(evec,rs1,ca_f,sum_s)
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
  integer rs1,rs2,ii,j,hi,lo,k,i
  RTYPE en_incr,pfac,pfac2,pfac3,pfac4,dum1,dum2,dum3,dd1,dd2
  RTYPE evec(MAXENERGYTERMS)
  RTYPE termc(rs_nbl(rs1)%nnbats)
  RTYPE term0(rs_nbl(rs1)%nnbats)
  RTYPE term1(rs_nbl(rs1)%nnbats)
  RTYPE term2(rs_nbl(rs1)%nnbats)
  RTYPE term3(rs_nbl(rs1)%nnbats)
  RTYPE term4(rs_nbl(rs1)%nnbats)
  RTYPE termh1(rs_nbl(rs1)%nnbats)
  RTYPE termh2(rs_nbl(rs1)%nnbats)
  RTYPE termh(rs_nbl(rs1)%nnbats)
  RTYPE termhv(rs_nbl(rs1)%nnbats,3)
  RTYPE termq(rs_nbl(rs1)%nnbats)
  RTYPE terms(rs_nbl(rs1)%nnbats)
  RTYPE termr(rs_nbl(rs1)%nnbats)
  RTYPE termsq(rs_nbl(rs1)%nnbats)
  RTYPE termdii(rs_nbl(rs1)%nnbats,3)
  RTYPE termdkk(rs_nbl(rs1)%nnbats,3)
  RTYPE termsqdr(rs_nbl(rs1)%nnbats,3)
  RTYPE termckk(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE id1(rs_nbl(rs1)%nnbats)
  RTYPE d1(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  integer sh(rs_nbl(rs1)%nnbs),hira,lora
  integer k1(rs_nbl(rs1)%nnbs,3),k2(rs_nbl(rs1)%nnbs,3)
  RTYPE drav(rs_nbl(rs1)%nnbs,3),srav(rs_nbl(rs1)%nnbs,3)
  RTYPE ca_f(n,3),sum_s(n)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
  termckk(:) = 0.0
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
    termq(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termr(lo:hi) = atr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsq(lo:hi) = scrq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsqdr(lo:hi,1) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1)
    termsqdr(lo:hi,2) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2)
    termsqdr(lo:hi,3) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3)
    lo = hi + 1
  end do
  pfac = scale_POLAR*electric
  pfac2 = par_IMPSOLV(8)*electric
  pfac3 = 1./par_IMPSOLV(1)
  pfac4 = scale_POLAR*par_IMPSOLV(9)*pfac3
  dum2 = 1.0 - par_IMPSOLV(9)
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    if (atq(ii).eq.0.0) cycle
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    d1(:) = sqrt(dis2(:))
    id1(:) = 1.0/d1(:)
    id2(:) = id1(:)*id1(:)
!   screened Coulomb potential
    call Vgenmu_all2(scrq(ii),termsq,i_sqm,terms,term1,term2,rs_nbl(rs1)%nnbats)
    dum3 = electric*atq(ii)
    termh1(:) = dum3*termq(:)*terms(:)
    termh(:) = atr(ii)+termr(:)
    dum3 = atq(ii)*pfac2
    termh2(:) = dum3*termq(:)/termh(:)
    term0(:) = pfac*id1(:)
    lo = 1
    termhv(:,:) = 0.0
    do k=1,rs_nbl(rs1)%nnbs
      rs2 = rs_nbl(rs1)%nb(k)
      hi = lo + sh(k)
      do j=lo,hi
        if (abs(termh1(j)).gt.abs(termh2(j))) then
          term4(j) = termh1(j)
          term3(j) = 1.0
        else
          dd1 = (d1(j) - termh(j))*pfac3
          dd2 = 1.0 - dd1
          if (dd1.lt.0.0) then
            term3(j) = dum2
            term4(j) = term3(j)*termh1(j) + par_IMPSOLV(9)*termh2(j)
            term0(j) = term3(j)*term0(j)
          else if (dd2.gt.0.0) then
            dum1 = par_IMPSOLV(9)*dd2
            term3(j) = (1.0 - dum1)
            term4(j) = term3(j)*termh1(j) + dum1*termh2(j)
            term0(j) = term3(j)*term0(j)
            dum1 = -pfac4*id2(j)*(termh2(j)-termh1(j))
            termhv(j,1) = dum1*dvec(j,1)
            termhv(j,2) = dum1*dvec(j,2)
            termhv(j,3) = dum1*dvec(j,3)
          else
            term4(j) = termh1(j)
            term3(j) = 1.0
          end if
        end if
        dum3 = atq(ii)*term0(j)*termq(j)
        sum_s(ii) = sum_s(ii) + dum3*term1(j)
        termckk(j) = termckk(j) + dum3*term2(j)
      end do
      lo = hi + 1
    end do
!
    dum3 = pfac*atq(ii)
    term0(:) = dum3*termq(:)*id1(:)
    termh(:) = term0(:)*term3(:)*term2(:)
    termdkk(:,1) = -termhv(:,1) + termh(:)*termsqdr(:,1)
    termdkk(:,2) = -termhv(:,2) + termh(:)*termsqdr(:,2)
    termdkk(:,3) = -termhv(:,3) + termh(:)*termsqdr(:,3)
    termh(:) = term0(:)*term3(:)*term1(:)
    termdii(:,1) = termhv(:,1) + termh(:)*scrq_dr(ii,1)
    termdii(:,2) = termhv(:,2) + termh(:)*scrq_dr(ii,2)
    termdii(:,3) = termhv(:,3) + termh(:)*scrq_dr(ii,3)
    termc(:) = scale_POLAR*term4(:)*id1(:)
    en_incr = sum(termc(:))
    evec(6) = evec(6) + en_incr
!
    termh(:) = termc(:)*id2(:)
    foin(:,1) = dvec(:,1)*termh(:)
    ca_f(ii,1) = ca_f(ii,1) - (sum(foin(:,1)) - sum(termdii(:,1)))
    for_k(:,1) = for_k(:,1) + foin(:,1) + termdkk(:,1)
    foin(:,2) = dvec(:,2)*termh(:)
    ca_f(ii,2) = ca_f(ii,2) - (sum(foin(:,2)) - sum(termdii(:,2)))
    for_k(:,2) = for_k(:,2) + foin(:,2) + termdkk(:,2)
    foin(:,3) = dvec(:,3)*termh(:)
    ca_f(ii,3) = ca_f(ii,3) - (sum(foin(:,3)) - sum(termdii(:,3)))
    for_k(:,3) = for_k(:,3) + foin(:,3) + termdkk(:,3)
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
    sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) = &
 &     sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + termckk(lo:hi)
    lo = hi + 1
  end do
!
end
!
!-------------------------------------------------------------------
!
! for the mixed models, the fxn gets pretty messy, and hence provides a minor boost
! relative to Vforce_rsp_long(...)
!
subroutine Vforce_PSCRM39_C(evec,rs1,ca_f,sum_s)
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
  integer rs1,rs2,ii,j,hi,lo,k,i
  RTYPE en_incr,pfac,pfac2,pfac3,pfac4,dum1,dum2,dum3,dd1,dd2
  RTYPE evec(MAXENERGYTERMS)
  RTYPE termc(rs_nbl(rs1)%nnbats)
  RTYPE term0(rs_nbl(rs1)%nnbats)
  RTYPE term3(rs_nbl(rs1)%nnbats)
  RTYPE term4(rs_nbl(rs1)%nnbats)
  RTYPE termh1(rs_nbl(rs1)%nnbats)
  RTYPE termh2(rs_nbl(rs1)%nnbats)
  RTYPE termh(rs_nbl(rs1)%nnbats)
  RTYPE termhv(rs_nbl(rs1)%nnbats,3)
  RTYPE termq(rs_nbl(rs1)%nnbats)
  RTYPE termr(rs_nbl(rs1)%nnbats)
  RTYPE termsq(rs_nbl(rs1)%nnbats)
  RTYPE termdii(rs_nbl(rs1)%nnbats,3)
  RTYPE termdkk(rs_nbl(rs1)%nnbats,3)
  RTYPE termsqdr(rs_nbl(rs1)%nnbats,3)
  RTYPE termckk(rs_nbl(rs1)%nnbats)
  RTYPE foin(rs_nbl(rs1)%nnbats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbats,3)
  RTYPE id2(rs_nbl(rs1)%nnbats)
  RTYPE dis2(rs_nbl(rs1)%nnbats)
  RTYPE id1(rs_nbl(rs1)%nnbats)
  RTYPE d1(rs_nbl(rs1)%nnbats)
  RTYPE for_k(rs_nbl(rs1)%nnbats,3)
  integer sh(rs_nbl(rs1)%nnbs),hira,lora
  integer k1(rs_nbl(rs1)%nnbs,3),k2(rs_nbl(rs1)%nnbs,3)
  RTYPE drav(rs_nbl(rs1)%nnbs,3),srav(rs_nbl(rs1)%nnbs,3)
  RTYPE ca_f(n,3),sum_s(n)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
  termckk(:) = 0.0
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
    termq(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termr(lo:hi) = atr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsq(lo:hi) = scrq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsqdr(lo:hi,1) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1)
    termsqdr(lo:hi,2) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2)
    termsqdr(lo:hi,3) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3)
    lo = hi + 1
  end do
  pfac = scale_POLAR*electric
  pfac2 = par_IMPSOLV(8)*electric
  pfac3 = 1./par_IMPSOLV(1)
  pfac4 = scale_POLAR*par_IMPSOLV(9)*pfac3
  dum2 = 1.0 - par_IMPSOLV(9)
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    if (atq(ii).eq.0.0) cycle
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    d1(:) = sqrt(dis2(:))
    id1(:) = 1.0/d1(:)
    id2(:) = id1(:)*id1(:)
!   screened Coulomb potential
    dum3 = electric*atq(ii)*scrq(ii)
    termh1(:) = dum3*termq(:)*termsq(:)
    termh(:) = atr(ii)+termr(:)
    dum3 = atq(ii)*pfac2
    termh2(:) = dum3*termq(:)/termh(:)
    term0(:) = pfac*id1(:)
    lo = 1
    termhv(:,:) = 0.0
    do k=1,rs_nbl(rs1)%nnbs
      rs2 = rs_nbl(rs1)%nb(k)
      hi = lo + sh(k)
      do j=lo,hi
        if (abs(termh1(j)).gt.abs(termh2(j))) then
          term4(j) = termh1(j)
          term3(j) = 1.0
        else
          dd1 = (d1(j) - termh(j))*pfac3
          dd2 = 1.0 - dd1
          if (dd1.lt.0.0) then
            term3(j) = dum2
            term4(j) = term3(j)*termh1(j) + par_IMPSOLV(9)*termh2(j)
            term0(j) = term3(j)*term0(j)
          else if (dd2.gt.0.0) then
            dum1 = par_IMPSOLV(9)*dd2
            term3(j) = (1.0 - dum1)
            term4(j) = term3(j)*termh1(j) + dum1*termh2(j)
            term0(j) = term3(j)*term0(j)
            dum1 = -pfac4*id2(j)*(termh2(j)-termh1(j))
            termhv(j,1) = dum1*dvec(j,1)
            termhv(j,2) = dum1*dvec(j,2)
            termhv(j,3) = dum1*dvec(j,3)
          else
            term4(j) = termh1(j)
            term3(j) = 1.0
          end if
        end if
        dum3 = atq(ii)*term0(j)*termq(j)
        sum_s(ii) = sum_s(ii) + dum3*termsq(j)
        termckk(j) = termckk(j) + dum3*scrq(ii)
      end do
      lo = hi + 1
    end do
!
    dum3 = pfac*atq(ii)
    term0(:) = dum3*termq(:)*id1(:)
    termh(:) = term0(:)*term3(:)*scrq(ii)
    termdkk(:,1) = -termhv(:,1) + termh(:)*termsqdr(:,1)
    termdkk(:,2) = -termhv(:,2) + termh(:)*termsqdr(:,2)
    termdkk(:,3) = -termhv(:,3) + termh(:)*termsqdr(:,3)
    termh(:) = term0(:)*term3(:)*termsq(:)
    termdii(:,1) = termhv(:,1) + termh(:)*scrq_dr(ii,1)
    termdii(:,2) = termhv(:,2) + termh(:)*scrq_dr(ii,2)
    termdii(:,3) = termhv(:,3) + termh(:)*scrq_dr(ii,3)
    termc(:) = scale_POLAR*term4(:)*id1(:)
    en_incr = sum(termc(:))
    evec(6) = evec(6) + en_incr
!
    termh(:) = termc(:)*id2(:)
    foin(:,1) = dvec(:,1)*termh(:)
    ca_f(ii,1) = ca_f(ii,1) - (sum(foin(:,1)) - sum(termdii(:,1)))
    for_k(:,1) = for_k(:,1) + foin(:,1) + termdkk(:,1)
    foin(:,2) = dvec(:,2)*termh(:)
    ca_f(ii,2) = ca_f(ii,2) - (sum(foin(:,2)) - sum(termdii(:,2)))
    for_k(:,2) = for_k(:,2) + foin(:,2) + termdkk(:,2)
    foin(:,3) = dvec(:,3)*termh(:)
    ca_f(ii,3) = ca_f(ii,3) - (sum(foin(:,3)) - sum(termdii(:,3)))
    for_k(:,3) = for_k(:,3) + foin(:,3) + termdkk(:,3)
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
    sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) = &
 &     sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + termckk(lo:hi)
    lo = hi + 1
  end do
!
end
!
!-------------------------------------------------------------------
!
! for the mixed models, the fxn gets pretty messy, and hence provides a minor boost
! relative to Vforce_rsp_long(...)
!
subroutine Vforce_PSCRM78_TR(evec,rs1,ca_f,sum_s)
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
  integer rs1,rs2,ii,j,hi,lo,k,i
  RTYPE en_incr,pfac,pfac2,pfac3,pfac4,dum1,dum2,dum3,dd1,dd2
  RTYPE evec(MAXENERGYTERMS)
  RTYPE termc(rs_nbl(rs1)%nnbtrats)
  RTYPE term0(rs_nbl(rs1)%nnbtrats)
  RTYPE term1(rs_nbl(rs1)%nnbtrats)
  RTYPE term2(rs_nbl(rs1)%nnbtrats)
  RTYPE term3(rs_nbl(rs1)%nnbtrats)
  RTYPE term4(rs_nbl(rs1)%nnbtrats)
  RTYPE termh1(rs_nbl(rs1)%nnbtrats)
  RTYPE termh2(rs_nbl(rs1)%nnbtrats)
  RTYPE termh(rs_nbl(rs1)%nnbtrats)
  RTYPE termhv(rs_nbl(rs1)%nnbtrats,3)
  RTYPE termq(rs_nbl(rs1)%nnbtrats)
  RTYPE terms(rs_nbl(rs1)%nnbtrats)
  RTYPE termr(rs_nbl(rs1)%nnbtrats)
  RTYPE termsq(rs_nbl(rs1)%nnbtrats)
  RTYPE termdii(rs_nbl(rs1)%nnbtrats,3)
  RTYPE termdkk(rs_nbl(rs1)%nnbtrats,3)
  RTYPE termsqdr(rs_nbl(rs1)%nnbtrats,3)
  RTYPE termckk(rs_nbl(rs1)%nnbtrats)
  RTYPE foin(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbtrats,3)
  RTYPE id2(rs_nbl(rs1)%nnbtrats)
  RTYPE dis2(rs_nbl(rs1)%nnbtrats)
  RTYPE id1(rs_nbl(rs1)%nnbtrats)
  RTYPE d1(rs_nbl(rs1)%nnbtrats)
  RTYPE for_k(rs_nbl(rs1)%nnbtrats,3)
  integer sh(rs_nbl(rs1)%nnbtrs)
  RTYPE ca_f(n,3),sum_s(n)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
  termckk(:) = 0.0
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nbb + at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(k)
    hi = lo + sh(k)
    dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,3)
    termq(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termr(lo:hi) = atr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsq(lo:hi) = scrq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsqdr(lo:hi,1) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1)
    termsqdr(lo:hi,2) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2)
    termsqdr(lo:hi,3) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3)
    lo = hi + 1
  end do
  pfac = scale_POLAR*electric
  pfac2 = par_IMPSOLV(8)*electric
  pfac3 = 1./par_IMPSOLV(1)
  pfac4 = scale_POLAR*par_IMPSOLV(9)*pfac3
  dum2 = 1.0 - par_IMPSOLV(9)
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    if (atq(ii).eq.0.0) cycle
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    d1(:) = sqrt(dis2(:))
    id1(:) = 1.0/d1(:)
    id2(:) = id1(:)*id1(:)
!   screened Coulomb potential
    call Vgenmu_all2(scrq(ii),termsq,i_sqm,terms,term1,term2,rs_nbl(rs1)%nnbtrats)
    dum3 = electric*atq(ii)
    termh1(:) = dum3*termq(:)*terms(:)
    termh(:) = atr(ii)+termr(:)
    dum3 = atq(ii)*pfac2
    termh2(:) = dum3*termq(:)/termh(:)
    term0(:) = pfac*id1(:)
    lo = 1
    termhv(:,:) = 0.0
    do k=1,rs_nbl(rs1)%nnbtrs
      rs2 = rs_nbl(rs1)%nbtr(k)
      hi = lo + sh(k)
      do j=lo,hi
        if (abs(termh1(j)).gt.abs(termh2(j))) then
          term4(j) = termh1(j)
          term3(j) = 1.0
        else
          dd1 = (d1(j) - termh(j))*pfac3
          dd2 = 1.0 - dd1
          if (dd1.lt.0.0) then
            term3(j) = dum2
            term4(j) = term3(j)*termh1(j) + par_IMPSOLV(9)*termh2(j)
            term0(j) = term3(j)*term0(j)
          else if (dd2.gt.0.0) then
            dum1 = par_IMPSOLV(9)*dd2
            term3(j) = (1.0 - dum1)
            term4(j) = term3(j)*termh1(j) + dum1*termh2(j)
            term0(j) = term3(j)*term0(j)
            dum1 = -pfac4*id2(j)*(termh2(j)-termh1(j))
            termhv(j,1) = dum1*dvec(j,1)
            termhv(j,2) = dum1*dvec(j,2)
            termhv(j,3) = dum1*dvec(j,3)
          else
            term4(j) = termh1(j)
            term3(j) = 1.0
          end if
        end if
        dum3 = atq(ii)*term0(j)*termq(j)
        sum_s(ii) = sum_s(ii) + dum3*term1(j)
        termckk(j) = termckk(j) + dum3*term2(j)
      end do
      lo = hi + 1
    end do
!
    dum3 = pfac*atq(ii)
    term0(:) = dum3*termq(:)*id1(:)
    termh(:) = term0(:)*term3(:)*term2(:)
    termdkk(:,1) = -termhv(:,1) + termh(:)*termsqdr(:,1)
    termdkk(:,2) = -termhv(:,2) + termh(:)*termsqdr(:,2)
    termdkk(:,3) = -termhv(:,3) + termh(:)*termsqdr(:,3)
    termh(:) = term0(:)*term3(:)*term1(:)
    termdii(:,1) = termhv(:,1) + termh(:)*scrq_dr(ii,1)
    termdii(:,2) = termhv(:,2) + termh(:)*scrq_dr(ii,2)
    termdii(:,3) = termhv(:,3) + termh(:)*scrq_dr(ii,3)
    termc(:) = scale_POLAR*term4(:)*id1(:)
    en_incr = sum(termc(:))
    evec(6) = evec(6) + en_incr
!
    termh(:) = termc(:)*id2(:)
    foin(:,1) = dvec(:,1)*termh(:)
    ca_f(ii,1) = ca_f(ii,1) - (sum(foin(:,1)) - sum(termdii(:,1)))
    for_k(:,1) = for_k(:,1) + foin(:,1) + termdkk(:,1)
    foin(:,2) = dvec(:,2)*termh(:)
    ca_f(ii,2) = ca_f(ii,2) - (sum(foin(:,2)) - sum(termdii(:,2)))
    for_k(:,2) = for_k(:,2) + foin(:,2) + termdkk(:,2)
    foin(:,3) = dvec(:,3)*termh(:)
    ca_f(ii,3) = ca_f(ii,3) - (sum(foin(:,3)) - sum(termdii(:,3)))
    for_k(:,3) = for_k(:,3) + foin(:,3) + termdkk(:,3)
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
    sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) = &
 &     sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + termckk(lo:hi)
    lo = hi + 1
  end do
!
end
!
!-------------------------------------------------------------------
!
! for the mixed models, the fxn gets pretty messy, and hence provides a minor boost
! relative to Vforce_rsp_long(...)
!
subroutine Vforce_PSCRM39_TR(evec,rs1,ca_f,sum_s)
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
  integer rs1,rs2,ii,j,hi,lo,k,i
  RTYPE en_incr,pfac,pfac2,pfac3,pfac4,dum1,dum2,dum3,dd1,dd2
  RTYPE evec(MAXENERGYTERMS)
  RTYPE termc(rs_nbl(rs1)%nnbtrats)
  RTYPE term0(rs_nbl(rs1)%nnbtrats)
  RTYPE term3(rs_nbl(rs1)%nnbtrats)
  RTYPE term4(rs_nbl(rs1)%nnbtrats)
  RTYPE termh1(rs_nbl(rs1)%nnbtrats)
  RTYPE termh2(rs_nbl(rs1)%nnbtrats)
  RTYPE termh(rs_nbl(rs1)%nnbtrats)
  RTYPE termhv(rs_nbl(rs1)%nnbtrats,3)
  RTYPE termq(rs_nbl(rs1)%nnbtrats)
  RTYPE termr(rs_nbl(rs1)%nnbtrats)
  RTYPE termsq(rs_nbl(rs1)%nnbtrats)
  RTYPE termdii(rs_nbl(rs1)%nnbtrats,3)
  RTYPE termdkk(rs_nbl(rs1)%nnbtrats,3)
  RTYPE termsqdr(rs_nbl(rs1)%nnbtrats,3)
  RTYPE termckk(rs_nbl(rs1)%nnbtrats)
  RTYPE foin(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dveci(rs_nbl(rs1)%nnbtrats,3)
  RTYPE dvec(rs_nbl(rs1)%nnbtrats,3)
  RTYPE id2(rs_nbl(rs1)%nnbtrats)
  RTYPE dis2(rs_nbl(rs1)%nnbtrats)
  RTYPE id1(rs_nbl(rs1)%nnbtrats)
  RTYPE d1(rs_nbl(rs1)%nnbtrats)
  RTYPE for_k(rs_nbl(rs1)%nnbtrats,3)
  integer sh(rs_nbl(rs1)%nnbtrs)
  RTYPE ca_f(n,3),sum_s(n)
!  
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
  termckk(:) = 0.0
!
! loop over complete set of atom-atom interactions in distant residues
  lo = 1
  sh = at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nbb + at(rs_nbl(rs1)%nbtr(1:rs_nbl(rs1)%nnbtrs))%nsc - 1
! this loop does not vectorize cleanly -> as few operations as possible in here
  do k=1,rs_nbl(rs1)%nnbtrs
    rs2 = rs_nbl(rs1)%nbtr(k)
    hi = lo + sh(k)
    dveci(lo:hi,1) = x(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,1)
    dveci(lo:hi,2) = y(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,2)
    dveci(lo:hi,3) = z(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + rs_nbl(rs1)%trsvec(k,3)
    termq(lo:hi) = atq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termr(lo:hi) = atr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsq(lo:hi) = scrq(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k))
    termsqdr(lo:hi,1) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),1)
    termsqdr(lo:hi,2) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),2)
    termsqdr(lo:hi,3) = scrq_dr(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k),3)
    lo = hi + 1
  end do
  pfac = scale_POLAR*electric
  pfac2 = par_IMPSOLV(8)*electric
  pfac3 = 1./par_IMPSOLV(1)
  pfac4 = scale_POLAR*par_IMPSOLV(9)*pfac3
  dum2 = 1.0 - par_IMPSOLV(9)
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%nbb+at(rs1)%nsc-1
    if (atq(ii).eq.0.0) cycle
    i = ii - at(rs1)%bb(1) + 1
!   use vector forms as much as possible
!   first distance vectors
    dvec(:,1) = dveci(:,1) - x(ii)
    dvec(:,2) = dveci(:,2) - y(ii)
    dvec(:,3) = dveci(:,3) - z(ii)
    dis2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
    d1(:) = sqrt(dis2(:))
    id1(:) = 1.0/d1(:)
    id2(:) = id1(:)*id1(:)
!   screened Coulomb potential
    dum3 = electric*atq(ii)*scrq(ii)
    termh1(:) = dum3*termq(:)*termsq(:)
    termh(:) = atr(ii)+termr(:)
    dum3 = atq(ii)*pfac2
    termh2(:) = dum3*termq(:)/termh(:)
    term0(:) = pfac*id1(:)
    lo = 1
    termhv(:,:) = 0.0
    do k=1,rs_nbl(rs1)%nnbtrs
      rs2 = rs_nbl(rs1)%nbtr(k)
      hi = lo + sh(k)
      do j=lo,hi
        if (abs(termh1(j)).gt.abs(termh2(j))) then
          term4(j) = termh1(j)
          term3(j) = 1.0
        else
          dd1 = (d1(j) - termh(j))*pfac3
          dd2 = 1.0 - dd1
          if (dd1.lt.0.0) then
            term3(j) = dum2
            term4(j) = term3(j)*termh1(j) + par_IMPSOLV(9)*termh2(j)
            term0(j) = term3(j)*term0(j)
          else if (dd2.gt.0.0) then
            dum1 = par_IMPSOLV(9)*dd2
            term3(j) = (1.0 - dum1)
            term4(j) = term3(j)*termh1(j) + dum1*termh2(j)
            term0(j) = term3(j)*term0(j)
            dum1 = -pfac4*id2(j)*(termh2(j)-termh1(j))
            termhv(j,1) = dum1*dvec(j,1)
            termhv(j,2) = dum1*dvec(j,2)
            termhv(j,3) = dum1*dvec(j,3)
          else
            term4(j) = termh1(j)
            term3(j) = 1.0
          end if
        end if
        dum3 = atq(ii)*term0(j)*termq(j)
        sum_s(ii) = sum_s(ii) + dum3*termsq(j)
        termckk(j) = termckk(j) + dum3*scrq(ii)
      end do
      lo = hi + 1
    end do
!
    dum3 = pfac*atq(ii)
    term0(:) = dum3*termq(:)*id1(:)
    termh(:) = term0(:)*term3(:)*scrq(ii)
    termdkk(:,1) = -termhv(:,1) + termh(:)*termsqdr(:,1)
    termdkk(:,2) = -termhv(:,2) + termh(:)*termsqdr(:,2)
    termdkk(:,3) = -termhv(:,3) + termh(:)*termsqdr(:,3)
    termh(:) = term0(:)*term3(:)*termsq(:)
    termdii(:,1) = termhv(:,1) + termh(:)*scrq_dr(ii,1)
    termdii(:,2) = termhv(:,2) + termh(:)*scrq_dr(ii,2)
    termdii(:,3) = termhv(:,3) + termh(:)*scrq_dr(ii,3)
    termc(:) = scale_POLAR*term4(:)*id1(:)
    en_incr = sum(termc(:))
    evec(6) = evec(6) + en_incr
!
    termh(:) = termc(:)*id2(:)
    foin(:,1) = dvec(:,1)*termh(:)
    ca_f(ii,1) = ca_f(ii,1) - (sum(foin(:,1)) - sum(termdii(:,1)))
    for_k(:,1) = for_k(:,1) + foin(:,1) + termdkk(:,1)
    foin(:,2) = dvec(:,2)*termh(:)
    ca_f(ii,2) = ca_f(ii,2) - (sum(foin(:,2)) - sum(termdii(:,2)))
    for_k(:,2) = for_k(:,2) + foin(:,2) + termdkk(:,2)
    foin(:,3) = dvec(:,3)*termh(:)
    ca_f(ii,3) = ca_f(ii,3) - (sum(foin(:,3)) - sum(termdii(:,3)))
    for_k(:,3) = for_k(:,3) + foin(:,3) + termdkk(:,3)
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
    sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) = &
 &     sum_s(at(rs2)%bb(1):at(rs2)%bb(1)+sh(k)) + termckk(lo:hi)
    lo = hi + 1
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine Vforce_dcbdscrq(rs,ca_f)
!
  use polypep
  use energies
  use forces
  use atoms
!
  implicit none
!
  integer i,rs,ii,ll,l,j
  RTYPE ca_f(n,3)
!
  do i=1,at(rs)%npol
    ii = at(rs)%pol(i)
    do l=1,sisa_nr(ii)
      ll = sisa_i(l,ii)
      do j=1,3
!        write(*,*)  sum_scrcbs(ii)*sisq_dr(l,j,ii),cben_dr(l,j,ii)
        ca_f(ll,j) = ca_f(ll,j) + (sum_scrcbs(ii)+sum_scrcbs_tr(ii)+sum_scrcbs_lr(ii))*sisq_dr(l,j,ii)
      end do
    end do
  end do
!
end
!
!
!-------------------------------------------------------------------
!

