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
!------------------------------------------------------------------------------
!
! the vectorized version of en_rsp
! note that the speed-up hinges crucially on the arrays being statically 
! allocated at call-time (for some reason dynamic allocation is terribly slow),
! and on defining the range -> 1:alcsz
! note that for thread-safe execution we either have to specify multiple instances
! of svte(:) and pass them on as arguments or to organize the loop on the outside (difficult for
! pairwise sum). the latter would make this fxn very dangerous to use blindly.
! a third option is to ATOMICize the statements of concern (svte increments) which is
! slow but safe  ! Martin : Would be good to update this comment to reflect what has been chosen 
!
subroutine Ven_rsp(evec,rs1,rs2,cut)
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
  use fyoc
  ! Martin 
  !use mpistuff
              
  use mcsums
  use energies
 ! use martin_own ! This is for debuging 
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
  RTYPE tmaxd(at(rs1)%na*at(rs2)%na)
  integer idx2(at(rs1)%na*at(rs2)%na,2)
  RTYPE evec(MAXENERGYTERMS),incr
  RTYPE svec(3),efvoli,efvolk,datri,datrk,pfac,rt23,pwc2
  logical cut,ismember,doid2
  integer olap(at(rs1)%na*at(rs2)%na)
  RTYPE d2diff(at(rs1)%na*at(rs2)%na)
  integer term_softcore(at(rs1)%na*at(rs2)%na) ! martin : added
  RTYPE temp_pka_par(at(rs1)%na*at(rs2)%na)
  RTYPE temp_scale_IPP,temp_scale_attLJ
  RTYPE temp1,temp2
    
  character*10 rep
! no intialization needed: note that on the calling side it requires a lot of
! care to handle fxns which "silently" increment arguments
  term_softcore(:)=0
  temp_pka_par(:)=0.! This might need to be set to one

  
  
!  print *,"ABCDEFGHIJKL",rs1,rs2
!  print *,nrsnb(rs1),nrsnb(rs2)
!  temp1=evec(1)
!  temp2=evec(3)
  if (ideal_run.EQV..true.) then
    return
  end if
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
! for a crosslinked residue special exclusions apply -> we use a dirty
! trick here that's causing recursive execution of Ven_rsp
! we rely on: - crosslinked residues never being sequence neighbors
!             - all MC energy evaluations done via Ven_rsp and Ven_rsp_long
  if (disulf(rs1).eq.rs2) then
    call Ven_crosslink_corr(crlk_idx(rs1),evec,cut)
    return
  end if
! in FEG-1 we branch out if precisely one of the residues is ghosting
! always remember that the interaction of two ghost-residues (incl self!)
! is the full Hamiltonian, not the ghosted one
! in FEG-2 we branch out if either residue is ghosted
  if (use_FEG.EQV..true.) then
    if (fegmode.eq.1) then
        if(((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..false.)).OR.&
        & ((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..true.))) then
        call Ven_rsp_feg(evec,rs1,rs2,cut)
        return
        end if
    else if (fegmode.eq.2) then
      if((par_FEG(rs1).EQV..true.).OR.(par_FEG(rs2).EQV..true.))then
        call Ven_rsp_feg(evec,rs1,rs2,cut)
        return
      end if
    end if
  end if
  doid2 = .false.
  if  ((use_attLJ.EQV..true.).OR.&
 &     ((use_IPP.EQV..true.).AND.(use_hardsphere.EQV..false.))) then
    doid2 = .true.
  end if
  
  
  if (do_hs.eqv..true.) then 
      temp_scale_attLJ=scale_attLJ
      temp_scale_IPP=scale_IPP
      do i=1,size(hs_residues) ! This is specific to the LJ arrays, as the fudges cannot be set appart from one another.
          if ((hs_residues(i).eq.rs1).or.(hs_residues(i).eq.rs2)) then 
            scale_attLJ=hs_ATTLJ*par_HSQ(hs_residues(i),2)+(1-par_HSQ(hs_residues(i),2))*temp_scale_attLJ
            scale_IPP=hs_REPLJ*par_HSQ(hs_residues(i),2)+(1-par_HSQ(hs_residues(i),2))*temp_scale_IPP
            exit
          end if 
      end do 
  end if 
  
  
!  print *,"Here"
!  !print *,par_HSQ(:,2)  
!  print *,do_HSQ
!  print *,rs1,rs2
!  print *,cut
!  print *,use_softcore
!  
  
  
! three different cases: intra-residue (rs1 == rs2), neighbors in seq., or others
  if (rs1.eq.rs2) then
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
          if (use_IMPSOLV.EQV..true.) termr(k) = atr(ii)+atr(kk) ! Martin : atr is replica dependent, (as opposed to LJ_RAD))
          idx2(k,1) = i
          if (((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)).and.(use_softcore.eqv..true.)) then
            if ((abs(atr_limits(ii,1)).lt.1.0D-5).or.(abs(atr_limits(ii,2+spec_limit(atmres(ii)))).lt.1.0D-5)&
            &.or.(abs(atr_limits(kk,1)).lt.1.0D-5).or.&
            &(abs(atr_limits(kk,2+spec_limit(atmres(kk)))).lt.1.0D-5)) then ! Martin : if any of the end state is a singularity
              term_softcore(k)=1
              temp_pka_par(k)=1.
              !Correct for the directionality
              ! Note that for atoms that have no null limits, this factor is one in any case
              if (do_hsq.eqv..true.) then
                  if (abs(atr_limits(ii,1)).lt.1.0D-5) then 
                      temp_pka_par(k)=temp_pka_par(k)*par_HSQ(rs1,2)                  
                  end if 
                  if (abs(atr_limits(kk,1)).lt.1.0D-5) then 
                      temp_pka_par(k)=temp_pka_par(k)*par_HSQ(rs2,2)
                  end if 
                  if (abs(atr_limits(ii,2+spec_limit(atmres(ii)))).lt.1.0D-5) then 
                      temp_pka_par(k)=temp_pka_par(k)*(1-par_HSQ(rs1,2))
                  end if 
                  if (abs(atr_limits(kk,2+spec_limit(atmres(kk)))).lt.1.0D-5) then 
                      temp_pka_par(k)=temp_pka_par(k)*(1-par_HSQ(rs2,2))
                  end if 
              
              else 
                  if ((abs(atr_limits(ii,1)).lt.1.0D-5).or.(abs(atr_limits(kk,1)).lt.1.0D-5)) then 
                      temp_pka_par(k)=temp_pka_par(k)*par_pka(2)                  
                  end if 
                  if ((abs(atr_limits(ii,2+spec_limit(atmres(ii)))).lt.1.0D-5).or.&
                    &(abs(atr_limits(kk,2+spec_limit(atmres(kk)))).lt.1.0D-5)) then 
                      temp_pka_par(k)=temp_pka_par(k)*(1-par_pka(2))
                  end if 
              end if 
              
!                    print *,"ii",ii
!                    print *,atr_limits(ii,1),atr_limits(ii,2+spec_limit(atmres(ii))),par_HSQ(rs1,2)
!                    print *,"kk",kk
!                    print *,atr_limits(kk,1),atr_limits(kk,2+spec_limit(atmres(kk))),par_HSQ(rs2,2)
!                    print *,"temp",temp_pka_par(k)
                !@ MArtin : Don't see why I made this exclusive test, deleting
!                ! It actually has to be exclusive
!              if (do_hsq.eqv..true.) then 
!                  if ((seq_q_state(rs1).ne.0).and.(seq_q_state(rs2).ne.0)) then
!                      par_pka(2)=(par_hsq(rs1,2)+par_hsq(rs2,2))/2.
!                  else if (seq_q_state(rs1).ne.0) then 
!                      par_pka(2)=par_hsq(rs1,2)
!                  else if (seq_q_state(rs2).ne.0) then 
!                      par_pka(2)=par_hsq(rs2,2) 
!                  else
!                      par_pka(2)=0.
!                  end if 
!              end if 
!                if (((abs(atr_limits(ii,1)).lt.1.0D-5).and.(abs(atr_limits(kk,2+spec_limit(atmres(kk)))).lt.1.0D-5)).or.&
!                        &((abs(atr_limits(ii,2+spec_limit(atmres(ii)))).lt.1.0D-5).and.(abs(atr_limits(kk,1)).lt.1.0D-5))) then 
!                ! If they are both equal to zero (special case) at different ends
!                ! If one or more of the original state is equal to zero 
!                    if (do_hsq.eqv..true.) then 
!                        if (par_pka(2).le.0.5) then 
!                            temp_pka_par(k)=par_pka(2)
!                        else 
!                            temp_pka_par(k)=1-par_pka(2)
!                        end if 
!                    else 
!                        if (sqrt(par_hsq(rs1,2)*par_hsq(rs2,2)).le.0.5) then 
!                            temp_pka_par(k)=max(par_hsq(rs1,2),par_hsq(rs2,2))
!                        else 
!                            temp_pka_par(k)=1-max(par_hsq(rs1,2),par_hsq(rs2,2))
!                        end if                         
!                    end if 
!
!                ! If the first limit is the dummy atom
!                else if ((abs(atr_limits(ii,1)).lt.1.0D-5).or.(abs(atr_limits(kk,1)).lt.1.0D-5)) then  
!                    if (do_hsq.eqv..true.) then 
!                        temp_pka_par(k)=max(par_hsq(rs1,2),par_hsq(rs2,2))
!                    else
!                        temp_pka_par(k)=par_pka(2)
!                    end if 
!                ! If the second limit is the zero limit
!                else if ((abs(atr_limits(ii,2+spec_limit(atmres(ii)))).lt.1.0D-5).or.(abs(atr_limits(kk,2+spec_limit(atmres(kk))))&
!                    &.lt.1.0D-5)) then 
!                    if (do_hsq.eqv..true.) then 
!                        temp_pka_par(k)=1-max(par_hsq(rs1,2),par_hsq(rs2,2))
!                    else 
!                        temp_pka_par(k)=1-par_pka(2)
!                    end if 
!                else 
!                    print *,"Unaccounted for case"
!                    call fexit()
!                end if 
                
!                print *,"ii"
!                print *,atr_limits(ii,1),atr_limits(ii,2+spec_limit(atmres(ii))),par_HSQ(rs1,2)
!                print *,"kk"
!                print *,atr_limits(kk,1),atr_limits(kk,2+spec_limit(atmres(kk))),par_HSQ(rs2,2)
!                print *,"temp",temp_pka_par(k)
            else   
               term_softcore(k)=0
            end if
          else
              term_softcore(k)=0 
          end if 
          k = k + 1
        end if
      end if
    end do
    k = k - 1
    alcsz = k ! Martin : So that is the number of (potentially) interacting atoms for this residue

    if (alcsz.gt.0) then

      if (doid2.EQV..true.) then
        id2(1:alcsz) = 1.0/d2(1:alcsz)
      end if 
      if (use_IPP.EQV..true.) then
        if (use_hardsphere.EQV..true.) then
          d2diff(1:alcsz) = max((terms(1:alcsz)-d2(1:alcsz))/terms(1:alcsz),0.0d0)
          olap(1:alcsz) = ceiling(d2diff(1:alcsz))
          term1(1:alcsz) = dble(olap(1:alcsz))
          evec(1) = evec(1) + screenbarrier*sum(term1(1:alcsz)*(1-term_softcore(1:alcsz)))
          if (use_softcore.eqv..true.) then 
            print *,"!!! you should not be using hardspheres with the sofcore potential !!!"
            call fexit()
          end if 
        else 
          term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**nhalf! Martin nhalf is half of the default repulsve exponent (6))++
          
          evec(1) = evec(1) + scale_IPP*4.0*&
 &                  sum(term0(1:alcsz)*term1(1:alcsz)*(1-term_softcore(1:alcsz))) ! home : I changed that 
        end if
      end if
      if (use_attLJ.EQV..true.) then           
        term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**3! Martin : terms is the corresponding LJ sigma squared
        ! Martin : And id2 is the invert square distance 
        ! Martin : It appears the terms(1:alcsz) actually is the squared sigma
        evec(3) = evec(3) - scale_attLJ*4.0*&
 &                  sum(term0(1:alcsz)*term1(1:alcsz)*(1-term_softcore(1:alcsz))) ! Martin : So that is where the LJ potential is 
      end if
        ! Martin : Softcore potential : I will actually have to loop within the atoms to find those that need change
      if (use_IPP.EQV..true.) then ! -12 part 
        term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3! Martin nhalf is half of the default repulsve exponent (6))
        evec(1) = evec(1) + scale_IPP*4.0*sum(temp_pka_par(1:alcsz)*term0(1:alcsz)*(&
    &   (sc_alpha*(1-temp_pka_par(1:alcsz))+term1(1:alcsz))**(-2))*term_softcore(1:alcsz))
      end if 
      
      if (use_attLJ.EQV..true.) then ! -6 part 
        term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3 ! Martin : It appears the terms(1:alcsz) actually is the squared sigma
        evec(3) = evec(3) - scale_attLJ*4.0*&
    &   sum(temp_pka_par(1:alcsz)*term0(1:alcsz)*((sc_alpha*(1-temp_pka_par(1:alcsz))+term1(1:alcsz))**(-1))*term_softcore(1:alcsz))
      end if
      
      if (use_IMPSOLV.EQV..true.) then
        d1(1:alcsz) = sqrt(d2(1:alcsz))! MArtin : That is the actual distance 
        tmaxd(1:alcsz) = termr(1:alcsz) + par_IMPSOLV(1) !Martin : I still have a problem with that  
        do k=1,alcsz
          ii = iaa(rs)%atin(idx2(k,1),1)
          kk = iaa(rs)%atin(idx2(k,1),2)
          ! Martin : I need to put that back : This is just a test 
          if ((abs(atr(ii)).lt.1.0D-5).or.(abs(atr(kk)).lt.1.0D-5)) then 
                 cycle
          end if 
          
          if (d1(k).lt.tmaxd(k)) then
            efvoli = atsavred(ii)*atvol(ii) !
            efvolk = atsavred(kk)*atvol(kk) ! Martin : fraction times a volumes 
            datri = 2.0*atr(ii)
            datrk = 2.0*atr(kk)

            if (d1(k).gt.(tmaxd(k)-datrk)) then ! Martin : Shouldn't it be the radius here ? instead of the diametre
              incr = -(tmaxd(k)-d1(k))/datrk
            else if (d1(k).lt.termr(k)) then
              incr = -d1(k)/termr(k)
            else
              incr = -1.0
            end if
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
            svte(ii) = svte(ii) + incr*efvolk
            
            if (d1(k).gt.(tmaxd(k)-datri)) then
              incr = -(tmaxd(k)-d1(k))/datri
            else if (d1(k).lt.termr(k)) then
              incr = -d1(k)/termr(k)
            else
              incr = -1.0
            end if
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
            svte(kk) = svte(kk) + incr*efvoli
          end if
        end do
      end if
      
      if (use_WCA.EQV..true.) then
        pfac = scale_WCA*par_WCA(2)
        rt23 = ROOT26*ROOT26
        pwc2 = par_WCA(1)*par_WCA(1)
        do k=1,alcsz
          if (term0(k).gt.0.0) then
            if (d2(k).lt.rt23*terms(k)) then
              term1(k) = (terms(k)/d2(k))**3
              term1(k) = 4.0*(term1(k)**2 - term1(k)) + 1.0
              evec(5) = evec(5) + scale_WCA*term0(k)*&
 &                               (term1(k) - par_WCA(2))
            else if (d2(k).lt.pwc2*terms(k)) then
              term1(k) = par_WCA(3)*d2(k)/terms(k) + par_WCA(4)
              term1(k) = 0.5*cos(term1(k)) - 0.5
              evec(5) = evec(5) + pfac*term0(k)*term1(k)
            end if
          end if
        end do
      end if
    end if
!
    
!   add correction terms
    if (use_CORR.EQV..true.) then
      call e_corrector(rs,evec)
    end if

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
      dvec(k,2) = y(kk) - y(ii) + svec(2)! Martin : those vectors are zero in the case of our sphere
      dvec(k,3) = z(kk) - z(ii) + svec(3)
      
      d2(k) = dvec(k,1)**2 + dvec(k,2)**2 + dvec(k,3)**2
      
      if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.&
 &      (cut.EQV..false.)) then !Martin : as 
        term0(k) = fudge(rs)%rsnb(i)*fudge(rs)%rsnb_lje(i) !lj_eps(attyp(ii),attyp(kk))
        if ((term0(k).gt.0.0).OR.(use_IMPSOLV.EQV..true.)) then
          terms(k) = fudge(rs)%rsnb_ljs(i) !lj_sig(attyp(ii),attyp(kk))
          if (use_IMPSOLV.EQV..true.) termr(k) = atr(ii)+atr(kk)
          idx2(k,1) = i
          if (((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)).and.(use_softcore.eqv..true.)) then
            if ((abs(atr_limits(ii,1)).lt.1.0D-5).or.(abs(atr_limits(ii,2+spec_limit(atmres(ii)))).lt.1.0D-5)&
            &.or.(abs(atr_limits(kk,1)).lt.1.0D-5).or.&
            &(abs(atr_limits(kk,2+spec_limit(atmres(kk)))).lt.1.0D-5)) then 
                ! Martin : if any of the end state is a singularity
              term_softcore(k)=1
              temp_pka_par(k)=1.
              !Correct for the directionality
              ! Note that for atoms that hve no null limits, this factor is one in any case
              if (do_hsq.eqv..true.) then
                  if (abs(atr_limits(ii,1)).lt.1.0D-5) then 
                      temp_pka_par(k)=temp_pka_par(k)*par_HSQ(rs1,2)                  
                  end if 
                  if (abs(atr_limits(kk,1)).lt.1.0D-5) then 
                      temp_pka_par(k)=temp_pka_par(k)*par_HSQ(rs2,2)
                  end if 
                  if (abs(atr_limits(ii,2+spec_limit(atmres(ii)))).lt.1.0D-5) then 
                      temp_pka_par(k)=temp_pka_par(k)*(1-par_HSQ(rs1,2))
                  end if 
                  if (abs(atr_limits(kk,2+spec_limit(atmres(kk)))).lt.1.0D-5) then 
                      temp_pka_par(k)=temp_pka_par(k)*(1-par_HSQ(rs2,2))
                  end if 
              else 
                  if ((abs(atr_limits(ii,1)).lt.1.0D-5).or.(abs(atr_limits(kk,1)).lt.1.0D-5)) then 
                      temp_pka_par(k)=temp_pka_par(k)*par_pka(2)                  
                  end if 
                  if ((abs(atr_limits(ii,2+spec_limit(atmres(ii)))).lt.1.0D-5).or.&
                    &(abs(atr_limits(kk,2+spec_limit(atmres(kk)))).lt.1.0D-5)) then 
                      temp_pka_par(k)=temp_pka_par(k)*(1-par_pka(2))
                  end if 
              end if 
!                print *,"ii",ii
!                    print *,atr_limits(ii,1),atr_limits(ii,2+spec_limit(atmres(ii))),par_HSQ(rs1,2)
!                    print *,"kk",kk
!                    print *,atr_limits(kk,1),atr_limits(kk,2+spec_limit(atmres(kk))),par_HSQ(rs2,2)
!                    print *,"temp",temp_pka_par(k)
                
!                if (((abs(atr_limits(ii,1)).lt.1.0D-5).and.(abs(atr_limits(kk,2+spec_limit(atmres(kk)))).lt.1.0D-5)).or.&
!                    &((abs(atr_limits(ii,2+spec_limit(atmres(ii)))).lt.1.0D-5).and.(abs(atr_limits(kk,1)).lt.1.0D-5))) then 
!                        ! If they are both equal to zero (special case) at different ends
!                        ! 
!                    if (do_hsq.eqv..true.) then 
!                        if (sqrt(par_hsq(rs1,2)*par_hsq(rs2,2)).le.0.5) then 
!                            temp_pka_par(k)=max(par_hsq(rs1,2),par_hsq(rs2,2))
!                        else 
!                            temp_pka_par(k)=1-max(par_hsq(rs1,2),par_hsq(rs2,2))
!                        end if 
!                    else 
!                        if (par_pka(2).le.0.5) then 
!                            temp_pka_par(k)=par_pka(2)
!                        else 
!                            temp_pka_par(k)=1-par_pka(2)
!                        end if 
!                    end if 
!                ! If the first limit is the dummy atom
!                else if ((abs(atr_limits(ii,1)).lt.1.0D-5).or.(abs(atr_limits(kk,1)).lt.1.0D-5)) then  
!                    if (do_hsq.eqv..true.) then 
!!                        print *,"This can never happen."
!                        temp_pka_par(k)=max(par_hsq(rs1,2),par_hsq(rs2,2))
!                    else
!                        temp_pka_par(k)=par_pka(2)
!                    end if 
!                ! If the second limit is the zero limit
!                else if ((abs(atr_limits(ii,2+spec_limit(atmres(ii)))).lt.1.0D-5).or.(abs(atr_limits(kk,2+spec_limit(atmres(kk))))&
!                    &.lt.1.0D-5)) then 
!                    if (do_hsq.eqv..true.) then 
!                        temp_pka_par(k)=1-max(par_hsq(rs1,2),par_hsq(rs2,2))
!                    else 
!                        temp_pka_par(k)=1-par_pka(2)
!                    end if 
!                else 
!                    print *,"Unaccounted for case"
!                    call fexit()
!                end if 
!!                print *,"ii"
!!                print *,atr_limits(ii,1),atr_limits(ii,2+spec_limit(atmres(ii))),par_HSQ(rs1,2)
!!                print *,"kk"
!!                print *,atr_limits(kk,1),atr_limits(kk,2+spec_limit(atmres(kk))),par_HSQ(rs2,2)
!!                print *,"temp",temp_pka_par(k)
!!                
            else   
               term_softcore(k)=0
            end if
          else
              term_softcore(k)=0 
          end if 
          k = k + 1
        end if
      end if
    end do

    k = k - 1
    alcsz = k ! Martin : seriously that's a variable name ?
    ! Martin : k is the number of atoms in the interact array     
    if (alcsz.gt.0) then
      if (doid2.EQV..true.) then
        id2(1:alcsz) = 1.0/d2(1:alcsz)
      end if
      if (use_IPP.EQV..true.) then
        if (use_hardsphere.EQV..true.) then
          d2diff(1:alcsz) = max((terms(1:alcsz)-d2(1:alcsz))/terms(1:alcsz),0.0d0)
          olap(1:alcsz) = ceiling(d2diff(1:alcsz))
          term1(1:alcsz) = dble(olap(1:alcsz))
          evec(1) = evec(1) + screenbarrier*sum(term1(1:alcsz)*(1-term_softcore(1:alcsz)))
        else
          term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**nhalf! Martin nhalf is half of the default repulsive exponent (6))
          evec(1) = evec(1) + scale_IPP*4.0*&
 &                  sum(term0(1:alcsz)*term1(1:alcsz)*(1-term_softcore(1:alcsz)))
        end if
      end if
      if (use_attLJ.EQV..true.) then
        term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**3! Martin : terms is the corresponding LJ sigma
        !Martin : And id2 is the invert square distance 
        ! Martin : It appears the terms(1:alcsz) actually is the squared sigma
        evec(3) = evec(3) - scale_attLJ*4.0*&
 &                  sum(term0(1:alcsz)*term1(1:alcsz)*(1-term_softcore(1:alcsz))) ! Martin : So that is where the LJ potential is 
      end if
      
    ! Martin : Softcore potential : I will actually have to loop within the atoms to find those that need change
      if (use_IPP.EQV..true.) then ! -12 part
        term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3! Martin nhalf is half of the default repulsive exponent (6))
        evec(1) = evec(1) + scale_IPP*4.0*sum(temp_pka_par(1:alcsz)*term0(1:alcsz)*(&
    &   (sc_alpha*(1-temp_pka_par(1:alcsz))+term1(1:alcsz))**(-2))*term_softcore(1:alcsz))
      end if
      if (use_attLJ.EQV..true.) then ! -6 part 
        term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3 ! Martin :  terms(1:alcsz) actually is the squared sigma
        evec(3) = evec(3) - scale_attLJ*4.0*&
    &   sum(temp_pka_par(1:alcsz)*term0(1:alcsz)*((sc_alpha*(1-temp_pka_par(1:alcsz))+term1(1:alcsz))**(-1))*term_softcore(1:alcsz))
      end if
      
      if (use_IMPSOLV.EQV..true.) then
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        tmaxd(1:alcsz) = termr(1:alcsz) + par_IMPSOLV(1)
        do k=1,alcsz
          ii = iaa(rs)%atnb(idx2(k,1),1)
          kk = iaa(rs)%atnb(idx2(k,1),2)

              
          if ((abs(atr(ii)).lt.1.0D-5).or.(abs(atr(kk)).lt.1.0D-5)) then 
                cycle
          end if 
          if (d1(k).lt.tmaxd(k)) then
            efvoli = atsavred(ii)*atvol(ii)
            efvolk = atsavred(kk)*atvol(kk)
            datri = 2.0*atr(ii)
            datrk = 2.0*atr(kk)
            
            if (d1(k).gt.(tmaxd(k)-datrk)) then
              incr = -(tmaxd(k)-d1(k))/datrk
            else if (d1(k).lt.termr(k)) then
              incr = -d1(k)/termr(k)
            else
              incr = -1.0
            end if
            
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
            svte(ii) = svte(ii) + incr*efvolk
            if (d1(k).gt.(tmaxd(k)-datri)) then
              incr = -(tmaxd(k)-d1(k))/datri
            else if (d1(k).lt.termr(k)) then
              incr = -d1(k)/termr(k)
            else
              incr = -1.0
            end if
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
            
            svte(kk) = svte(kk) + incr*efvoli
          end if
        end do
      end if
      if (use_WCA.EQV..true.) then
        pfac = scale_WCA*par_WCA(2)
        rt23 = ROOT26*ROOT26
        pwc2 = par_WCA(1)*par_WCA(1)
        do k=1,alcsz
          if (term0(k).gt.0.0) then
            if (d2(k).lt.rt23*terms(k)) then
              term1(k) = (terms(k)/d2(k))**3
              term1(k) = 4.0*(term1(k)**2 - term1(k)) + 1.0
              evec(5) = evec(5) + scale_WCA*term0(k)*&
 &                               (term1(k) - par_WCA(2))
            else if (d2(k).lt.pwc2*terms(k)) then
              term1(k) = par_WCA(3)*d2(k)/terms(k) + par_WCA(4)
              term1(k) = 0.5*cos(term1(k)) - 0.5
              evec(5) = evec(5) + pfac*term0(k)*term1(k)
            end if
          end if
        end do
      end if
!
    end if
  else

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
        
        terms(k)=0.0
        term0(k)=0.0
        if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.&
 &      (cut.EQV..false.)) then
            if (do_hsq.eqv..true.) then 
!                if (epsrule.eq.2) then 
                    ! Adapted from the simga one 
                    term0(k) = lj_eps(attyp(ii),attyp(kk))*(1-par_hsq(rs1,2))*(1-par_hsq(rs2,2))+&! 1 1 
                    &lj_eps(attyp(ii),bio_ljtyp(transform_table(b_type(kk))))&
                    &*(1-par_hsq(rs1,2))*par_hsq(rs2,2)*(1-spec_limit(atmres(kk)))+&! 1 2
                    &lj_eps(bio_ljtyp(transform_table(b_type(ii))),attyp(kk))*&
                    &(par_hsq(rs1,2))*(1-par_hsq(rs2,2))*(1-spec_limit(atmres(ii)))+&! 2 1 
                    &lj_eps(bio_ljtyp(transform_table(b_type(ii))),bio_ljtyp(transform_table(b_type(kk))))*&
                    &(par_hsq(rs1,2))*(par_hsq(rs2,2))*(1-spec_limit(atmres(ii)))*(1-spec_limit(atmres(kk)))+&! 2 2
                    &lj_eps(attyp(ii),bio_ljtyp(his_eqv_table(b_type(kk))))*&
                    &(1-par_hsq(rs1,2))*par_hsq(rs2,2)*(spec_limit(atmres(kk)))+&!1 3 
                    &lj_eps(bio_ljtyp(his_eqv_table(b_type(ii))),attyp(kk))*&
                    &(1-par_hsq(rs2,2))*par_hsq(rs1,2)*(spec_limit(atmres(ii)))+&!3 1 
                    &lj_eps(bio_ljtyp(his_eqv_table(b_type(ii))),bio_ljtyp(transform_table(b_type(kk))))*&
                    &(par_hsq(rs2,2))*par_hsq(rs1,2)*(spec_limit(atmres(ii)))*(1-spec_limit(atmres(kk)))+&!3 2
                    &lj_eps(bio_ljtyp(transform_table(b_type(ii))),bio_ljtyp(his_eqv_table(b_type(kk))))*&
                    &(par_hsq(rs2,2))*par_hsq(rs1,2)*(1-spec_limit(atmres(ii)))*(spec_limit(atmres(kk)))+&!2 3
                    &lj_eps(bio_ljtyp(his_eqv_table(b_type(ii))),bio_ljtyp(his_eqv_table(b_type(kk))))*&
                    &(par_hsq(rs2,2))*par_hsq(rs1,2)*(spec_limit(atmres(ii)))*(spec_limit(atmres(kk)))!3 3
                    ! Old, working, but not for all 
!                    term0(k) = sqrt(lj_eps_limits(attyp(ii),attyp(ii),2+spec_limit(atmres(ii)))*&
!                    &lj_eps_limits(attyp(kk),attyp(kk),2+spec_limit(atmres(kk)))*par_hsq(rs1,2)*par_hsq(rs2,2))+&
!
!                    &sqrt(lj_eps_limits(attyp(ii),attyp(ii),1)*&                
!                    &lj_eps_limits(attyp(kk),attyp(kk),2+spec_limit(atmres(kk)))*(1-par_hsq(rs1,2))*par_hsq(rs2,2))+&
!
!                    &sqrt(lj_eps_limits(attyp(ii),attyp(ii),2+spec_limit(atmres(ii)))*&
!                    &lj_eps_limits(attyp(kk),attyp(kk),1)*par_hsq(rs1,2)*(1-par_hsq(rs2,2)))+&
!
!                    &sqrt(lj_eps_limits(attyp(ii),attyp(ii),1)*&
!                    &lj_eps_limits(attyp(kk),attyp(kk),1)*(1-par_hsq(rs1,2))*(1-par_hsq(rs2,2)))
                    ! Yes this is complicated, but it will have correct end states, which is all that matters
                    ! The reason for this is that in the HSQ the sig will actually vary, and combination between different sig 
                    ! will happen
!                else 
!                    print *,"EPSRULE that is not equal to 2 has not yet been implemented in HSQ"
!                end if 
            else     
              term0(k) = lj_eps(attyp(ii),attyp(kk))
            end if 
            
          if ((term0(k).gt.0.0).OR.(use_IMPSOLV.EQV..true.)) then
            if (do_hsq.eqv..true.) then 
!                if (sigrule.eq.1) then 
                    
                    
                    ! Need to acount for the third limit
                    !The following sum will be equal to only one of its 8 parts on the end states
                    ! This version take into account histidine tautomer.
                    !1 1 
                    !1 2
                    !2 1
                    !1 3
                    !3 1
                    !3 2
                    !2 3
                    !3 3
                
                  
                    ! Even simpler version would be with sigma itself 
                    ! ect...
                    ! WOuld actually stop the need for LJ limits to be setup.                       
                    terms(k) = lj_sig(attyp(ii),attyp(kk))*(1-par_hsq(rs1,2))*(1-par_hsq(rs2,2))+&! 1 1 
                    &lj_sig(attyp(ii),bio_ljtyp(transform_table(b_type(kk))))&
                    &*(1-par_hsq(rs1,2))*par_hsq(rs2,2)*(1-spec_limit(atmres(kk)))+&! 1 2
                    &lj_sig(bio_ljtyp(transform_table(b_type(ii))),attyp(kk))*&
                    &(par_hsq(rs1,2))*(1-par_hsq(rs2,2))*(1-spec_limit(atmres(ii)))+&! 2 1 
                    &lj_sig(bio_ljtyp(transform_table(b_type(ii))),bio_ljtyp(transform_table(b_type(kk))))*&
                    &(par_hsq(rs1,2))*(par_hsq(rs2,2))*(1-spec_limit(atmres(ii)))*(1-spec_limit(atmres(kk)))+&! 2 2
                    &lj_sig(attyp(ii),bio_ljtyp(his_eqv_table(b_type(kk))))*&
                    &(1-par_hsq(rs1,2))*par_hsq(rs2,2)*(spec_limit(atmres(kk)))+&!1 3 
                    &lj_sig(bio_ljtyp(his_eqv_table(b_type(ii))),attyp(kk))*&
                    &(1-par_hsq(rs2,2))*par_hsq(rs1,2)*(spec_limit(atmres(ii)))+&!3 1 
                    &lj_sig(bio_ljtyp(his_eqv_table(b_type(ii))),bio_ljtyp(transform_table(b_type(kk))))*&
                    &(par_hsq(rs2,2))*par_hsq(rs1,2)*(spec_limit(atmres(ii)))*(1-spec_limit(atmres(kk)))+&!3 2
                    &lj_sig(bio_ljtyp(transform_table(b_type(ii))),bio_ljtyp(his_eqv_table(b_type(kk))))*&
                    &(par_hsq(rs2,2))*par_hsq(rs1,2)*(1-spec_limit(atmres(ii)))*(spec_limit(atmres(kk)))+&!2 3
                    &lj_sig(bio_ljtyp(his_eqv_table(b_type(ii))),bio_ljtyp(his_eqv_table(b_type(kk))))*&
                    &(par_hsq(rs2,2))*par_hsq(rs1,2)*(spec_limit(atmres(ii)))*(spec_limit(atmres(kk)))!3 3
                    
!                    terms(k) = lj_sig_limits(attyp(ii),attyp(kk),1)*(1-par_hsq(rs1,2))*(1-par_hsq(rs2,2))+&! 1 1 
!                    &lj_sig_limits(attyp(ii),bio_ljtyp(transform_table(b_type(kk)),1)&
!                    &*(1-par_hsq(rs1,2))*par_hsq(rs2,2)*(1-spec_limit(atmres(kk)))&! 1 2
!                    &lj_sig_limits(attyp(ii),attyp(kk),2)*(par_hsq(rs1,2))*(1-par_hsq(rs2,2))*(1-spec_limit(atmres(ii)))&! 2 1 
!                    &lj_sig_limits(attyp(ii),bio_ljtyp(his_eqv_table(b_type(kk))),1)*&
!                    &(1-par_hsq(rs1,2))*par_hsq(rs2,2)*(spec_limit(atmres(kk)))&!1 3 
!                    &lj_sig_limits(attyp(ii),attyp(kk),3)*&
!                    &(1-par_hsq(rs2,2))*par_hsq(rs1,2)*(spec_limit(atmres(ii)))&!3 1 
!                    &lj_sig_limits(attyp(ii),bio_ljtyp(transform_table(b_type(kk)),3)*&
!                    &(par_hsq(rs2,2))*par_hsq(rs1,2)*(spec_limit(atmres(ii)))*(1-spec_limit(atmres(kk)))&!3 2
!                    &lj_sig_limits(attyp(ii),bio_ljtyp(his_eqv_table(b_type(kk))),2)*&
!                    &(par_hsq(rs2,2))*par_hsq(rs1,2)*(1-spec_limit(atmres(ii)))*(spec_limit(atmres(kk)))&!2 3
!                    &lj_sig_limits(attyp(ii),bio_ljtyp(his_eqv_table(b_type(kk))),3)*&
!                    &(par_hsq(rs2,2))*par_hsq(rs1,2)*(spec_limit(atmres(ii)))*(spec_limit(atmres(kk)))!3 3
                    
                    
!                    &0.5*lj_sig_limits(attyp(kk),attyp(ii),1)*(1-par_hsq(rs2,2))+&
!                    
!                    &0.5*lj_sig_limits(attyp(ii),attyp(kk),2+spec_limit(atmres(ii)))*(par_hsq(rs1,2))*(1-par_hsq(rs2,2))+&
!                    &0.5*lj_sig_limits(attyp(ii),bio_ljtyp(transform_table(b_type(kk))),2+spec_limit(atmres(ii)))&
!                    &*(par_hsq(rs1,2))*(par_hsq(rs2,2))+&
                    !&0.5*lj_sig_limits(attyp(kk),attyp(ii),2+spec_limit(atmres(kk)))*(par_hsq(rs2,2))
                    
                    ! Working, but maybe not for his
!                    terms(k) = 0.5*lj_sig_limits(attyp(ii),attyp(kk),1)*(1-par_hsq(rs1,2))+&
!                    &0.5*lj_sig_limits(attyp(kk),attyp(ii),1)*(1-par_hsq(rs2,2))+&
!                    
!                    &0.5*lj_sig_limits(attyp(ii),attyp(kk),2+spec_limit(atmres(ii)))*(par_hsq(rs1,2))+&
!                    &0.5*lj_sig_limits(attyp(kk),attyp(ii),2+spec_limit(atmres(kk)))*(par_hsq(rs2,2))
                    
                    ! Old not working 
!                    if ((spec_limit(atmres(ii)).eq.1).and.(spec_limit(atmres(kk)).eq.1))
!                    +&
!                    &lj_sig_limits(attyp(ii),attyp(kk),2)*0.5(par_hsq(rs1,2)+(par_hsq(rs2,2)))
!                    ! How do I ge the his limit ?
                    !1./2*(lj_sig_limits(attyp(ii),attyp(ii),2+spec_limit(atmres(ii)))+&
!                    &lj_sig_limits(attyp(kk),attyp(kk),2+spec_limit(atmres(kk))))*sqrt(par_hsq(rs1,2)*par_hsq(rs2,2))+&
!
!                    &1./2*(lj_sig_limits(attyp(ii),attyp(ii),1)+&
!                    &lj_sig_limits(attyp(kk),attyp(kk),2+spec_limit(atmres(kk))))*sqrt((1-par_hsq(rs1,2))*par_hsq(rs2,2))+&
!
!                    &1./2*(lj_sig_limits(attyp(ii),attyp(ii),2+spec_limit(atmres(ii)))+&
!                    &lj_sig_limits(attyp(kk),attyp(kk),1))*sqrt(par_hsq(rs1,2)*(1-par_hsq(rs2,2)))+&
!
!                    &1./2*(lj_sig_limits(attyp(ii),attyp(ii),1)+&
!                    &lj_sig_limits(attyp(kk),attyp(kk),1))*sqrt((1-par_hsq(rs1,2))*(1-par_hsq(rs2,2)))


!                else
!                    print *,"SIGRULE that is not equal to 1 has not yet been implemented in HSQ"
!                end if 
            else 
                terms(k) = lj_sig(attyp(ii),attyp(kk))
            end if 

            if (use_IMPSOLV.EQV..true.) termr(k) = atr(ii)+atr(kk)

            idx2(k,1) = ii
            idx2(k,2) = kk
            
            if (((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)).and.(use_softcore.eqv..true.)) then
                if ((abs(atr_limits(ii,1)).lt.1.0D-5).or.(abs(atr_limits(ii,2+spec_limit(atmres(ii)))).lt.1.0D-5)&
                &.or.(abs(atr_limits(kk,1)).lt.1.0D-5).or.&
                &(abs(atr_limits(kk,2+spec_limit(atmres(kk)))).lt.1.0D-5)) then ! Martin : if any of the end state is a singularity
                    term_softcore(k)=1
              temp_pka_par(k)=1.
              !Correct for the directionality
              ! Note that for atoms that have no null limits, this factor is one in any case
              if (do_hsq.eqv..true.) then
                  if (abs(atr_limits(ii,1)).lt.1.0D-5) then 
                      temp_pka_par(k)=temp_pka_par(k)*par_HSQ(rs1,2)                  
                  end if 
                  if (abs(atr_limits(kk,1)).lt.1.0D-5) then 
                      temp_pka_par(k)=temp_pka_par(k)*par_HSQ(rs2,2)
                  end if 
                  if (abs(atr_limits(ii,2+spec_limit(atmres(ii)))).lt.1.0D-5) then 
                      temp_pka_par(k)=temp_pka_par(k)*(1-par_HSQ(rs1,2))
                  end if 
                  if (abs(atr_limits(kk,2+spec_limit(atmres(kk)))).lt.1.0D-5) then 
                      temp_pka_par(k)=temp_pka_par(k)*(1-par_HSQ(rs2,2))
                  end if 
              
              else 
                  if ((abs(atr_limits(ii,1)).lt.1.0D-5).or.(abs(atr_limits(kk,1)).lt.1.0D-5)) then 
                      temp_pka_par(k)=temp_pka_par(k)*par_pka(2)                  
                  end if 
                  if ((abs(atr_limits(ii,2+spec_limit(atmres(ii)))).lt.1.0D-5).or.&
                    &(abs(atr_limits(kk,2+spec_limit(atmres(kk)))).lt.1.0D-5)) then 
                      temp_pka_par(k)=temp_pka_par(k)*(1-par_pka(2))
                  end if 
              end if 
              
              
              
!            print *,"ii",ii
!            print *,atr_limits(ii,1),atr_limits(ii,2+spec_limit(atmres(ii))),par_HSQ(rs1,2)
!            print *,"kk",kk
!            print *,atr_limits(kk,1),atr_limits(kk,2+spec_limit(atmres(kk))),par_HSQ(rs2,2)
!            print *,"temp",temp_pka_par(k)
!                    !@ MArtin : Don't see why I made this exclusive test, deleting 
!                    if (((abs(atr_limits(ii,1)).lt.1.0D-5).and.(abs(atr_limits(kk,2+spec_limit(atmres(kk)))).lt.1.0D-5)).or.&
!                        &((abs(atr_limits(ii,2+spec_limit(atmres(ii)))).lt.1.0D-5).and.(abs(atr_limits(kk,1)).lt.1.0D-5))) then 
!                        ! If they are both equal to zero (special case) at different ends
!                        ! If one or more of the original state is equal to zero 
!                         if (do_hsq.eqv..true.) then 
!                            if (sqrt(par_hsq(rs1,2)*par_hsq(rs2,2)).le.0.5) then 
!                                temp_pka_par(k)=max(par_hsq(rs1,2),par_hsq(rs2,2))
!                            else 
!                                temp_pka_par(k)=1-max(par_hsq(rs1,2),par_hsq(rs2,2))
!                            end if 
!                        else 
!                            if (par_pka(2).le.0.5) then 
!                                temp_pka_par(k)=par_pka(2)
!                            else 
!                                temp_pka_par(k)=1-par_pka(2)
!                            end if 
!                        end if 
!                    ! If the first limit is the dummy atom
!                    else if ((abs(atr_limits(ii,1)).lt.1.0D-5).or.(abs(atr_limits(kk,1)).lt.1.0D-5)) then  
!                        if (do_hsq.eqv..true.) then 
!                            temp_pka_par(k)=max(par_hsq(rs1,2),par_hsq(rs2,2))
!                        else
!                            temp_pka_par(k)=par_pka(2)
!                        end if 
!                    ! If the second limit is the zero limit
!                    else if ((abs(atr_limits(ii,2+spec_limit(atmres(ii)))).lt.1.0D-5).or.&
!                        &(abs(atr_limits(kk,2+spec_limit(atmres(kk)))).lt.1.0D-5)) then 
!                        if (do_hsq.eqv..true.) then 
!                            temp_pka_par(k)=1-max(par_hsq(rs1,2),par_hsq(rs2,2))
!                        else 
!                            temp_pka_par(k)=1-par_pka(2)
!                        end if 
!                    else 
!                        print *,"Unaccounted for case"
!                        call fexit()
!                    end if 
!                    print *,"ii"
!                    print *,atr_limits(ii,1),atr_limits(ii,2+spec_limit(atmres(ii))),par_HSQ(rs1,2)
!                    print *,"kk"
!                    print *,atr_limits(kk,1),atr_limits(kk,2+spec_limit(atmres(kk))),par_HSQ(rs2,2)
!                    print *,"temp",temp_pka_par(k)
                else   
                   term_softcore(k)=0
                end if
          else
              term_softcore(k)=0 
          end if  
!          if ((rs1.eq.3).and.(rs2.eq.5)) then
!            print *,temp_pka_par(k)
!            print *,ii,kk
!            print *,attyp(ii),attyp(kk)
!            print *,lj_sig(attyp(ii),attyp(kk))
!            print *,lj_sig(attyp(ii),attyp(ii))
!            print *,lj_sig(attyp(kk),attyp(kk))
!!            print *,lj_sig_limits(attyp(ii),attyp(kk),1)
!!            print *,lj_sig_limits(attyp(ii),attyp(ii),1)
!!            print *,lj_sig_limits(attyp(kk),attyp(kk),1)
!            print *,terms(k),term0(k)
!          end if 
          k = k + 1
          end if
        end if
      end do
    end do
    k = k - 1
    alcsz = k
    ! Martin L this is for debugging purposes only 
!    print *,"OK",use_attLJ,use_IPP
    if (alcsz.gt.0) then
      if (doid2.EQV..true.) then
        id2(1:alcsz) = 1.0/d2(1:alcsz)
      end if
      if (use_IPP.EQV..true.) then
        if (use_hardsphere.EQV..true.) then
          d2diff(1:alcsz) = max((terms(1:alcsz)-d2(1:alcsz))/terms(1:alcsz),0.0d0)
          olap(1:alcsz) = ceiling(d2diff(1:alcsz))
          term1(1:alcsz) = dble(olap(1:alcsz))
          evec(1) = evec(1) + screenbarrier*sum(term1(1:alcsz)*(1-term_softcore(1:alcsz)))
        else
          term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**nhalf! Martin nhalf is half of the default repulsive exponent (6))
          evec(1) = evec(1) + scale_IPP*4.0*&
 &                  sum(term0(1:alcsz)*term1(1:alcsz)*(1-term_softcore(1:alcsz)))
        end if
      end if
!
      if (use_attLJ.EQV..true.) then
        
        term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**3! Martin : terms is the corresponding LJ sigma
        ! Martin : And id2 is the invert square distance 
        ! Martin : It appears the terms(1:alcsz) actually is the squared sigma
        evec(3) = evec(3) - scale_attLJ*4.0*&
 &                  sum(term0(1:alcsz)*term1(1:alcsz)*(1-term_softcore(1:alcsz))) ! Martin : So that is where the LJ potential is 
      end if
    ! Martin : Softcore potential : I will actually have to loop within the atoms to find those that need change
      if (use_IPP.EQV..true.) then ! -12 part 
        term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3! Martin nhalf is half of the default repulsve exponent (6))
        evec(1) = evec(1) + scale_IPP*4.0*sum(temp_pka_par(1:alcsz)*term0(1:alcsz)*(& 
    &   (sc_alpha*(1-temp_pka_par(1:alcsz))+term1(1:alcsz))**(-2))*term_softcore(1:alcsz))
      end if

      if (use_attLJ.EQV..true.) then ! -6 part
        term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3 ! Martin : It appears the terms(1:alcsz) actually is the squared sigma
        evec(3) = evec(3) - scale_attLJ*4.0*&
    &   sum(temp_pka_par(1:alcsz)*term0(1:alcsz)*((sc_alpha*(1-temp_pka_par(1:alcsz))+term1(1:alcsz))**(-1))*term_softcore(1:alcsz))
      end if
      
      if (use_IMPSOLV.EQV..true.) then
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        tmaxd(1:alcsz) = termr(1:alcsz) + par_IMPSOLV(1)
        do k=1,alcsz
          ii = idx2(k,1)
          kk = idx2(k,2)
          if ((abs(atr(ii)).lt.1.0D-5).or.(abs(atr(kk)).lt.1.0D-5)) then 
                 cycle
          end if 

          if (d1(k).lt.tmaxd(k)) then
            
            efvoli = atsavred(ii)*atvol(ii) ! Martin : effective volume occupied by either atom 
            efvolk = atsavred(kk)*atvol(kk)
            datri = 2.0*atr(ii)! Martin : diameter 
            datrk = 2.0*atr(kk)! Martin : 
            if (d1(k).gt.(tmaxd(k)-datrk)) then ! if the distance between the atoms is more than the max distance minus 
                                                ! the diameter of the atom being considered 
              incr = -(tmaxd(k)-d1(k))/datrk ! then the increment is the distance of the overrlap dividied by the distance
              
                ! In esssence that makes incr the fraction of the current volume fraction that moves 
            else if (d1(k).lt.termr(k)) then
              incr = -d1(k)/termr(k)
            else
              incr = -1.0
            end if
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
            !Martin : So -1 is the totality of the atom is within the solvation shell : take this volume out
            ! The only problem is 
            svte(ii) = svte(ii) + incr*efvolk

            if (d1(k).gt.(tmaxd(k)-datri)) then
              incr = -(tmaxd(k)-d1(k))/datri
            else if (d1(k).lt.termr(k)) then
              incr = -d1(k)/termr(k)
            else
              incr = -1.0
            end if
#ifdef ENABLE_THREADS
!$OMP ATOMIC
#endif
            svte(kk) = svte(kk) + incr*efvoli
          end if
        end do
      end if

      if (use_WCA.EQV..true.) then
        pfac = scale_WCA*par_WCA(2)
        rt23 = ROOT26*ROOT26
        pwc2 = par_WCA(1)*par_WCA(1)
        do k=1,alcsz
          if (term0(k).gt.0.0) then
            if (d2(k).lt.rt23*terms(k)) then
              term1(k) = (terms(k)/d2(k))**3
              term1(k) = 4.0*(term1(k)**2 - term1(k)) + 1.0
              evec(5) = evec(5) + scale_WCA*term0(k)*&
 &                               (term1(k) - par_WCA(2))
            else if (d2(k).lt.pwc2*terms(k)) then
              term1(k) = par_WCA(3)*d2(k)/terms(k) + par_WCA(4)
              term1(k) = 0.5*cos(term1(k)) - 0.5
              evec(5) = evec(5) + pfac*term0(k)*term1(k)
            end if
          end if
        end do
      end if
!
    end if
    end if
    if (do_hs.eqv..true.) then 
        scale_attLJ=temp_scale_attLJ
        scale_IPP=temp_scale_IPP
  end if 
!  print *,evec(1)-temp1,evec(3)-temp2
  
end
!
!-----------------------------------------------------------------------
!
subroutine Ven_rsp_feg(evec,rs1,rs2,cut)
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
  RTYPE evec(MAXENERGYTERMS),svec(3)
  logical cut

 
  
! we have a potential override to cover in par_FEG3, which allows residues to be
! fully de-coupled at all times
! note that if fegmode is 1, only one of the two is ghosted, otherwise background Hamiltonian
! if fegmode is 2 and both are ghosted, we ensure the fully de-coupled ghost always remains
! fully ghosted so as to not pick up spurious interactions with increased coupling
  if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) &
 &  return
! three different cases: intra-residue (rs1 == rs2), neighbors in seq., or others

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
          k = k + 1
        end if
      end if
    end do
    k = k - 1
    alcsz = k
    if (alcsz.gt.0) then
      id2(1:alcsz) = 1.0/d2(1:alcsz)
      term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3 
      if (use_FEGS(1).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(2))
        evec(1) = evec(1) + par_FEG2(1)*4.0*&
 &      sum(term0(1:alcsz)*term2(1:alcsz)*term2(1:alcsz))
      end if
      if (use_FEGS(3).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(6))
        evec(3) = evec(3) - par_FEG2(5)*4.0*&
 &                  sum(term0(1:alcsz)*term2(1:alcsz))
      end if
!
    end if
!
!   add correction terms
    if (use_CORR.EQV..true.) then
      call e_corrector(rs,evec)
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
          k = k + 1
        end if
      end if
    end do
    k = k - 1
    alcsz = k
    if (alcsz.gt.0) then
      id2(1:alcsz) = 1.0/d2(1:alcsz)
      term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3 
      if (use_FEGS(1).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(2))
        evec(1) = evec(1) + par_FEG2(1)*4.0*&
 &      sum(term0(1:alcsz)*term2(1:alcsz)*term2(1:alcsz))
      end if
      if (use_FEGS(3).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(6))
        evec(3) = evec(3) - par_FEG2(5)*4.0*&
 &                  sum(term0(1:alcsz)*term2(1:alcsz))
      end if
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
            k = k + 1
          end if
        end if
      end do
    end do
    k = k - 1
    alcsz = k
    if (alcsz.gt.0) then
      id2(1:alcsz) = 1.0/d2(1:alcsz)
      term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3 
      if (use_FEGS(1).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(2))
        evec(1) = evec(1) + par_FEG2(1)*4.0*&
 &      sum(term0(1:alcsz)*term2(1:alcsz)*term2(1:alcsz))
      end if
      if (use_FEGS(3).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(6))
        evec(3) = evec(3) - par_FEG2(5)*4.0*&
 &                  sum(term0(1:alcsz)*term2(1:alcsz))
      end if
!
    end if
!
!
  end if
  
end
!
!-----------------------------------------------------------------------
!
! the vectorized version of en_rsp_long
! note that the speed-up hinges crucially on the arrays being statically 
! allocated at call-time (for some reason dynamic allocation is terribly slow),
! and on defining the range -> 1:alcsz
!
subroutine Ven_rsp_long(evec,rs1,rs2,cut)
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
  use fyoc
  use grandensembles
!
  implicit none
!
  integer rs,rs1,rs2,i,j,k,ii,kk,rhi,rlo,imol1,imol2,alcsz
  RTYPE evec(MAXENERGYTERMS)
  RTYPE genmu,svec(3),pfac,pfac2,hlp1,hlp2,hlp3
  RTYPE d2(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE d1(at(rs1)%na*at(rs2)%na),id1(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE termr(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  integer tbin(at(rs1)%na*at(rs2)%na)
  logical cut,ismember

  if ((ideal_run.EQV..true.).OR.&
 &    ((use_POLAR.EQV..false.).AND.(use_TABUL.EQV..false.))) then
    return
  end if

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
! for a crosslinked residue special exclusions apply -> we use a dirty
! trick here that's causing recursive execution of Ven_rsp
! we rely on: - crosslinked residues never being sequence neighbors
!             - all MC energy evaluations done via Ven_rsp and Ven_rsp_long
  if (disulf(rs1).eq.rs2) then
    call Ven_crosslink_corr_long(crlk_idx(rs1),evec,cut)
    return
  end if
!
!
! in FEG we branch out if precisely one of the residues is ghosting
! always remember that the interaction of two ghost-residues (incl self!)
! is the full Hamiltonian, not the ghosted one
  if (use_FEG.EQV..true.) then
    if (fegmode.eq.1) then
   if(((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..false.)).OR.&
 & ((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..true.))) then
        call Ven_rsp_long_feg(evec,rs1,rs2,cut)
        return
      end if
    else if (fegmode.eq.2) then
      if((par_FEG(rs1).EQV..true.).OR.(par_FEG(rs2).EQV..true.))then
        call Ven_rsp_long_feg(evec,rs1,rs2,cut)
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
              term0(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)*&
 &                            scrq(ii)*scrq(kk)

 
              
            end do
          else if (scrq_model.eq.4) then
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              term0(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
              term2(i) = 1.0/(atr(ii)+atr(kk))
            end do
          else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              term0(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)*&
 &                               scrq(ii)*scrq(kk)
              termr(i) = atr(ii)+atr(kk)
              term2(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)/termr(i)
            end do
          else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              term0(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)*&
 &                  genmu(scrq(ii),scrq(kk),i_sqm)
            end do
          else ! 7,8
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              term0(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)*&
 &                    genmu(scrq(ii),scrq(kk),i_sqm)
              termr(i) = atr(ii)+atr(kk)
              term2(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)/termr(i)
            end do
          end if
        else
          do i=1,nrpolintra(rs)
            ii = iaa(rs)%polin(i,1)
            kk = iaa(rs)%polin(i,2)
            dvec(i,1) = x(kk) - x(ii)
            dvec(i,2) = y(kk) - y(ii)
            dvec(i,3) = z(kk) - z(ii)
            term0(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
          end do
        end if
!       now perform the vectorizable bulk operations (maximal)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
!       finish up -> screening model specificity is high
        if (use_IMPSOLV.EQV..true.) then
          k = 0
          if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.&
 &            (scrq_model.eq.6)) then
 

            term1(1:alcsz) = term0(1:alcsz)*id1(1:alcsz)
            evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
            
          else if (scrq_model.eq.4) then
            do i=1,nrpolintra(rs)
              if (id1(i).le.term2(i)) then
                term1(i) = term0(i)*id1(i)
              else
                term1(i) = term0(i)*term2(i) 
              end if
            end do
            evec(6) = evec(6) + pfac*par_IMPSOLV(8)*&
 &                      sum(id1(1:alcsz)*term1(1:alcsz))
          else
            term2(1:alcsz) = term2(1:alcsz)*pfac2
            term0(1:alcsz) = term0(1:alcsz)*electric
            do i=1,nrpolintra(rs)
              if (abs(term0(i)).gt.abs(term2(i))) then
                term1(i) = term0(i)
              else
!               dissociated
                if ((termr(i)+par_IMPSOLV(1)).le.d1(i)) then
                  term1(i) = term0(i)
!               close contact
                else if (d1(i).lt.termr(i)) then
                  term1(i) = par_IMPSOLV(9)*(term2(i)-term0(i)) + &
 &                                      term0(i)
!               first-shell regime
                else
                  term1(i) = (termr(i)+par_IMPSOLV(1)-d1(i))/&
 &                           par_IMPSOLV(1)
                  term1(i)= (1.0-par_IMPSOLV(9)*term1(i))*term0(i)&
 &                 + par_IMPSOLV(9)*term1(i)*term2(i) 
                end if
              end if
            end do
            term1(1:alcsz) = term1(1:alcsz)*id1(1:alcsz)
            evec(6) = evec(6) + scale_POLAR*sum(term1(1:alcsz))
          end if
        else
          term1(1:alcsz) = term0(1:alcsz)*id1(1:alcsz)
          evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
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
          ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - &
 &        atmol(molofrs(atmres(tbp%lst(i,1))),1) 
          kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - &
 &        atmol(molofrs(atmres(tbp%lst(i,2))),1)
          k = k + 1
          dvec(k,1) = x(kk) - x(ii)
          dvec(k,2) = y(kk) - y(ii)
          dvec(k,3) = z(kk) - z(ii)
        end do
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        pfac = 1.0/tbp%res
        termr(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
        tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(termr(1:alcsz))))
        term0(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
        term1(1:alcsz) = term0(1:alcsz)*term0(1:alcsz) ! quadratic
        term2(1:alcsz) = term1(1:alcsz)*term0(1:alcsz) ! cubic
        k = 0
        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
          k = k + 1
          if (tbin(k).ge.tbp%bins) then
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
          else if (tbin(k).ge.1) then
            hlp1 = 2.0*term2(k) - 3.0*term1(k)
            hlp2 = term2(k) - term1(k)
            hlp3 = hlp2 - term1(k) + term0(k)
            termr(k) = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k)+1) + &
 &                   tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1))
            evec(9) = evec(9) + scale_TABUL*termr(k)
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
!   so we have to check for BC
!
!   code 2 indicates the residues are actually beyond the long cutoff-range, but are
!   included due to their charges, bail out into standard treatment if sequence neighbors
!   in the SAME molecule (since Ven_rsp_lrel does not handle fudges)
!   of course, the correctness of rsp_mat relies on cutoffs to be used -> check
    if (cut.EQV..true.) then
      if (molofrs(rs1).ne.molofrs(rs2)) then
        if (rsp_mat(min(rs1,rs2),max(rs1,rs2)).eq.2) then
          call dis_bound_rs(rs1,rs2,svec)
          call Ven_rsp_lrel(evec,rs1,rs2,svec)
          return
        end if
      end if
    end if
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
              dvec(i,1) = x(kk) - x(ii)! + svec(1)
              dvec(i,2) = y(kk) - y(ii)! + svec(2)
              dvec(i,3) = z(kk) - z(ii)! + svec(3)
              term0(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)*&
 &                            scrq(ii)*scrq(kk)
            end do
          else if (scrq_model.eq.4) then
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)! + svec(1)
              dvec(i,2) = y(kk) - y(ii)! + svec(2)
              dvec(i,3) = z(kk) - z(ii)! + svec(3)
              term0(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
              term2(i) = 1.0/(atr(ii)+atr(kk))
            end do
          else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)! + svec(1)
              dvec(i,2) = y(kk) - y(ii)! + svec(2)
              dvec(i,3) = z(kk) - z(ii)! + svec(3)
              term0(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)*&
 &                               scrq(ii)*scrq(kk)
              termr(i) = atr(ii)+atr(kk)
              term2(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)/termr(i)
            end do
          else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)! + svec(1)
              dvec(i,2) = y(kk) - y(ii)! + svec(2)
              dvec(i,3) = z(kk) - z(ii)! + svec(3)
              term0(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)*&
 &                  genmu(scrq(ii),scrq(kk),i_sqm)
            end do
          else ! 7,8
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)! + svec(1)
              dvec(i,2) = y(kk) - y(ii)! + svec(2)
              dvec(i,3) = z(kk) - z(ii)! + svec(3)
              term0(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)*&
 &                    genmu(scrq(ii),scrq(kk),i_sqm)
              termr(i) = atr(ii)+atr(kk)
              term2(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)/termr(i)
            end do
          end if
        else
          do i=1,nrpolnb(rs)
            ii = iaa(rs)%polnb(i,1)
            kk = iaa(rs)%polnb(i,2)
            dvec(i,1) = x(kk) - x(ii)! + svec(1)
            dvec(i,2) = y(kk) - y(ii)! + svec(2)
            dvec(i,3) = z(kk) - z(ii)! + svec(3)
            term0(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
          end do
        end if
!       now perform the vectorizable bulk operations (maximal)
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
!       finish up -> screening model specificity is high
        if (use_IMPSOLV.EQV..true.) then
          k = 0
          if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.&
 &            (scrq_model.eq.6)) then
            term1(1:alcsz) = term0(1:alcsz)*id1(1:alcsz)
            evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
          else if (scrq_model.eq.4) then
            do i=1,nrpolnb(rs)
              if (id1(i).le.term2(i)) then
                term1(i) = term0(i)*id1(i)
              else
                term1(i) = term0(i)*term2(i) 
              end if
            end do
            evec(6) = evec(6) + pfac*par_IMPSOLV(8)*&
 &                  sum(id1(1:alcsz)*term1(1:alcsz))

          else
            term2(1:alcsz) = term2(1:alcsz)*pfac2
            term0(1:alcsz) = term0(1:alcsz)*electric
            do i=1,nrpolnb(rs)
              if (abs(term0(i)).gt.abs(term2(i))) then
                term1(i) = term0(i)
              else
!               dissociated
                if ((termr(i)+par_IMPSOLV(1)).le.d1(i)) then
                  term1(i) = term0(i)
!               close contact
                else if (d1(i).lt.termr(i)) then
                  term1(i) = par_IMPSOLV(9)*(term2(i)-term0(i)) + &
 &                                      term0(i)
!               first-shell regime
                else
                  term1(i) = (termr(i)+par_IMPSOLV(1)-d1(i))/&
 &                           par_IMPSOLV(1)
                  term1(i)= (1.0-par_IMPSOLV(9)*term1(i))*term0(i)&
 &                 + par_IMPSOLV(9)*term1(i)*term2(i) 
                end if
              end if
            end do
            term1(1:alcsz) = term1(1:alcsz)*id1(1:alcsz)
            evec(6) = evec(6) + scale_POLAR*sum(term1(1:alcsz))

          end if
        else
          term1(1:alcsz) = term0(1:alcsz)*id1(1:alcsz)
          evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
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
        termr(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
        tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(termr(1:alcsz))))
        term0(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
        term1(1:alcsz) = term0(1:alcsz)*term0(1:alcsz) ! quadratic
        term2(1:alcsz) = term1(1:alcsz)*term0(1:alcsz) ! cubic
        k = 0
        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
          k = k + 1
          if (tbin(k).ge.tbp%bins) then
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
          else if (tbin(k).ge.1) then
            hlp1 = 2.0*term2(k) - 3.0*term1(k)
            hlp2 = term2(k) - term1(k)
            hlp3 = hlp2 - term1(k) + term0(k)
            termr(k) = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k)+1) + &
 &                   tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1))
            evec(9) = evec(9) + scale_TABUL*termr(k)
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
!   code 2 indicates the residues are actually beyond the long cutoff-range, but are
!   included due to their charges
!   of course, the correctness of rsp_mat relies on cutoffs to be used -> check
    if (cut.EQV..true.) then
      if (rsp_mat(min(rs1,rs2),max(rs1,rs2)).eq.2) then
        call Ven_rsp_lrel(evec,rs1,rs2,svec)
        return
      end if
    end if
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
                dvec(k,1) = x(kk) - x(ii)! + svec(1)
                dvec(k,2) = y(kk) - y(ii)! + svec(2)
                dvec(k,3) = z(kk) - z(ii)! + svec(3)
                term0(k) = atq(ii)*atq(kk)*scrq(ii)*scrq(kk)
              end do
            end do
          else if (scrq_model.eq.4) then
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)! + svec(1)
                dvec(k,2) = y(kk) - y(ii)! + svec(2)
                dvec(k,3) = z(kk) - z(ii)! + svec(3)
                term0(k) = atq(ii)*atq(kk)
                term2(k) = 1.0/(atr(ii)+atr(kk))
              end do
            end do
          else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)! + svec(1)
                dvec(k,2) = y(kk) - y(ii)! + svec(2)
                dvec(k,3) = z(kk) - z(ii)! + svec(3)
                term0(k) = atq(ii)*atq(kk)*scrq(ii)*scrq(kk)
                termr(k) = atr(ii)+atr(kk)
                term2(k) = atq(ii)*atq(kk)/termr(k)
              end do
            end do
          else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)! + svec(1)
                dvec(k,2) = y(kk) - y(ii)! + svec(2)
                dvec(k,3) = z(kk) - z(ii)! + svec(3)
                term0(k) = atq(ii)*atq(kk)*&
 &                  genmu(scrq(ii),scrq(kk),i_sqm)
              end do
            end do
          else ! 7,8
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)! + svec(1)
                dvec(k,2) = y(kk) - y(ii)! + svec(2)
                dvec(k,3) = z(kk) - z(ii)! + svec(3)
                term0(k) = atq(ii)*atq(kk)*&
 &                    genmu(scrq(ii),scrq(kk),i_sqm)
                termr(k) = atr(ii)+atr(kk)
                term2(k) = atq(ii)*atq(kk)/termr(k)
              end do
            end do
          end if
        else
          do i=1,at(rs1)%npol
            ii = at(rs1)%pol(i)
            do j=1,at(rs2)%npol
              kk = at(rs2)%pol(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)! + svec(1)
              dvec(k,2) = y(kk) - y(ii)! + svec(2)
              dvec(k,3) = z(kk) - z(ii)! + svec(3)
              term0(k) = atq(ii)*atq(kk)
            end do
          end do
        end if
!       now perform the vectorizable bulk operations (maximal)
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
!       finish up -> screening model specificity is high
        if (use_IMPSOLV.EQV..true.) then
          k = 0
          if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.&
 &            (scrq_model.eq.6)) then
            term1(1:alcsz) = term0(1:alcsz)*id1(1:alcsz)
            evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
          else if (scrq_model.eq.4) then
            do i=1,at(rs1)%npol
              do j=1,at(rs2)%npol
                k = k + 1
                if (id1(k).le.term2(k)) then
                  term1(k) = term0(k)*id1(k)
                else
                  term1(k) = term0(k)*term2(k) 
                end if
              end do
            end do
            evec(6) = evec(6) + par_IMPSOLV(8)*pfac*&
 &                       sum(id1(1:alcsz)*term1(1:alcsz))
          else
            term2(1:alcsz) = term2(1:alcsz)*pfac2
            term0(1:alcsz) = term0(1:alcsz)*electric
            do i=1,at(rs1)%npol
              do j=1,at(rs2)%npol
                k = k + 1
                if (abs(term0(k)).gt.abs(term2(k))) then
                  term1(k) = term0(k)
                else
!                 dissociated
                  if ((termr(k)+par_IMPSOLV(1)).le.d1(k)) then
                    term1(k) = term0(k)
!                 close contact
                  else if (d1(k).lt.termr(k)) then
                    term1(k) = par_IMPSOLV(9)*(term2(k)-term0(k)) + &
 &                                        term0(k)
!                 first-shell regime
                  else
                    term1(k) = (termr(k)+par_IMPSOLV(1)-d1(k))/&
 &                             par_IMPSOLV(1)
                    term1(k)= (1.0-par_IMPSOLV(9)*term1(k))*term0(k)&
 &                   + par_IMPSOLV(9)*term1(k)*term2(k) 
                  end if
                end if
              end do
            end do
            term1(1:alcsz) = term1(1:alcsz)*id1(1:alcsz)
            evec(6) = evec(6) + scale_POLAR*sum(term1(1:alcsz))
          end if
        else
          term1(1:alcsz) = term0(1:alcsz)*id1(1:alcsz)
          evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
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
          dvec(k,1) = x(kk) - x(ii)! + svec(1)
          dvec(k,2) = y(kk) - y(ii)! + svec(2)
          dvec(k,3) = z(kk) - z(ii)! + svec(3)
        end do
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        pfac = 1.0/tbp%res
        termr(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
        tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(termr(1:alcsz))))
        term0(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
        term1(1:alcsz) = term0(1:alcsz)*term0(1:alcsz) ! quadratic
        term2(1:alcsz) = term1(1:alcsz)*term0(1:alcsz) ! cubic
        k = 0
        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
          k = k + 1
          if (tbin(k).ge.tbp%bins) then
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
          else if (tbin(k).ge.1) then
            hlp1 = 2.0*term2(k) - 3.0*term1(k)
            hlp2 = term2(k) - term1(k)
            hlp3 = hlp2 - term1(k) + term0(k)
            termr(k) = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k)+1) + &
 &                   tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1))
            evec(9) = evec(9) + scale_TABUL*termr(k)
          else
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),1)
          end if
        end do
      end if
    end if
!       
    end if
    
  
end
!
!-----------------------------------------------------------------------
!
subroutine Ven_rsp_long_feg(evec,rs1,rs2,cut)
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
!
  implicit none
!
  integer rs,rs1,rs2,i,j,k,ii,kk,alcsz
  RTYPE evec(MAXENERGYTERMS)
  RTYPE genmu,svec(3),pfac,pfac2
  RTYPE d2(at(rs1)%na*at(rs2)%na)
  RTYPE id1s(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE d1(at(rs1)%na*at(rs2)%na),id1(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE termr(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  logical cut

! we have a potential override to cover in par_FEG3, which allows residues to be
! fully de-coupled at all times
! note that if fegmode is 1, only one of the two is ghosted, otherwise background Hamiltonian
! if fegmode is 2 and both are ghosted, we ensure the fully de-coupled ghost always remains
! fully ghosted so as to not pick up spurious interactions with increased coupling

  if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.))&
 &  return
!
  svec(:) = 0.0


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
        if (use_FEGS(4).EQV..true.) then
          if (scrq_model.le.2) then
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              term0(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)*&
 &                            scrq(ii)*scrq(kk)
            end do
          else if (scrq_model.eq.4) then
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              term0(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
              term2(i) = 1.0/(atr(ii)+atr(kk))
            end do
          else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              term0(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)*&
 &                               scrq(ii)*scrq(kk)
              termr(i) = atr(ii)+atr(kk)
              term2(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)/termr(i)
            end do
          else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              term0(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)*&
 &                  genmu(scrq(ii),scrq(kk),i_sqm)
            end do
          else ! 7,8
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              term0(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)*&
 &                    genmu(scrq(ii),scrq(kk),i_sqm)
              termr(i) = atr(ii)+atr(kk)
              term2(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)/termr(i)
            end do
          end if
        else
          do i=1,nrpolintra(rs)
            ii = iaa(rs)%polin(i,1)
            kk = iaa(rs)%polin(i,2)
            dvec(i,1) = x(kk) - x(ii)
            dvec(i,2) = y(kk) - y(ii)
            dvec(i,3) = z(kk) - z(ii)
            term0(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
          end do
        end if
!       now perform the vectorizable bulk operations (maximal)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id1s(1:alcsz) = 1.0/(d1(1:alcsz)+par_FEG2(10))
!       finish up -> screening model specificity is high
        if (use_FEGS(4).EQV..true.) then
          k = 0
          if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.&
 &            (scrq_model.eq.6)) then
            term1(1:alcsz) = term0(1:alcsz)*id1s(1:alcsz)
            evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
          else if (scrq_model.eq.4) then
            do i=1,nrpolintra(rs)
              if (id1(i).le.term2(i)) then
                term1(i) = term0(i)*id1(i)
              else
                term1(i) = term0(i)*term2(i) 
              end if
            end do
            evec(6) = evec(6) + pfac*par_IMPSOLV(8)*&
 &                      sum(id1s(1:alcsz)*term1(1:alcsz))
          else
            term2(1:alcsz) = term2(1:alcsz)*pfac2
            term0(1:alcsz) = term0(1:alcsz)*electric
            do i=1,nrpolintra(rs)
              if (abs(term0(i)).gt.abs(term2(i))) then
                term1(i) = term0(i)
              else
!               dissociated
                if ((termr(i)+par_IMPSOLV(1)).le.d1(i)) then
                  term1(i) = term0(i)
!               close contact
                else if (d1(i).lt.termr(i)) then
                  term1(i) = par_IMPSOLV(9)*(term2(i)-term0(i)) + &
 &                                      term0(i)
!               first-shell regime
                else
                  term1(i) = (termr(i)+par_IMPSOLV(1)-d1(i))/&
 &                           par_IMPSOLV(1)
                  term1(i)= (1.0-par_IMPSOLV(9)*term1(i))*term0(i)&
 &                 + par_IMPSOLV(9)*term1(i)*term2(i) 
                end if
              end if
            end do
            term1(1:alcsz) = term1(1:alcsz)*id1s(1:alcsz)
            evec(6) = evec(6) + par_FEG2(9)*sum(term1(1:alcsz))
          end if
        else
          term1(1:alcsz) = term0(1:alcsz)*id1s(1:alcsz)
          evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
        end if
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
!   code 2 indicates the residues are actually beyond the long cutoff-range, but are
!   included due to their charges, bail out into standard treatment if sequence neighbors
!   in the SAME molecule (since Ven_rsp_lrel does not handle fudges)
!   of course, the correctness of rsp_mat relies on cutoffs to be used -> check
!   (bail-out here is full minimum image convention)
    if (molofrs(rs1).ne.molofrs(rs2)) then
      if (cut.EQV..true.) then
        if (rsp_mat(min(rs1,rs2),max(rs1,rs2)).eq.2) then
          call dis_bound_rs(rs1,rs2,svec)
          call Ven_rsp_lrel_feg(evec,rs1,rs2,svec)
          return
        end if
      end if
    end if
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
        if (use_FEGS(4).EQV..true.) then
          if (scrq_model.le.2) then
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)! + svec(1)
              dvec(i,2) = y(kk) - y(ii)! + svec(2)
              dvec(i,3) = z(kk) - z(ii)! + svec(3)
              term0(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)*&
 &                            scrq(ii)*scrq(kk)
            end do
          else if (scrq_model.eq.4) then
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)! + svec(1)
              dvec(i,2) = y(kk) - y(ii)! + svec(2)
              dvec(i,3) = z(kk) - z(ii)! + svec(3)
              term0(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
              term2(i) = 1.0/(atr(ii)+atr(kk))
            end do
          else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)! + svec(1)
              dvec(i,2) = y(kk) - y(ii)! + svec(2)
              dvec(i,3) = z(kk) - z(ii)! + svec(3)
              term0(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)*&
 &                               scrq(ii)*scrq(kk)
              termr(i) = atr(ii)+atr(kk)
              term2(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)/termr(i)
            end do
          else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)! + svec(1)
              dvec(i,2) = y(kk) - y(ii)! + svec(2)
              dvec(i,3) = z(kk) - z(ii)! + svec(3)
              term0(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)*&
 &                  genmu(scrq(ii),scrq(kk),i_sqm)
            end do
          else ! 7,8
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)! + svec(1)
              dvec(i,2) = y(kk) - y(ii)! + svec(2)
              dvec(i,3) = z(kk) - z(ii)! + svec(3)
              term0(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)*&
 &                    genmu(scrq(ii),scrq(kk),i_sqm)
              termr(i) = atr(ii)+atr(kk)
              term2(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)/termr(i)
            end do
          end if
        else
          do i=1,nrpolnb(rs)
            ii = iaa(rs)%polnb(i,1)
            kk = iaa(rs)%polnb(i,2)
            dvec(i,1) = x(kk) - x(ii)! + svec(1)
            dvec(i,2) = y(kk) - y(ii)! + svec(2)
            dvec(i,3) = z(kk) - z(ii)! + svec(3)
            term0(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
          end do
        end if
!       now perform the vectorizable bulk operations (maximal)
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id1s(1:alcsz) = 1.0/(d1(1:alcsz)+par_FEG2(10))
!       finish up -> screening model specificity is high
        if (use_FEGS(4).EQV..true.) then
          k = 0
          if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.&
 &            (scrq_model.eq.6)) then
            term1(1:alcsz) = term0(1:alcsz)*id1s(1:alcsz)
            evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
          else if (scrq_model.eq.4) then
            do i=1,nrpolnb(rs)
              if (id1(i).le.term2(i)) then
                term1(i) = term0(i)*id1(i)
              else
                term1(i) = term0(i)*term2(i) 
              end if
            end do
            evec(6) = evec(6) + pfac*par_IMPSOLV(8)*&
 &                  sum(id1s(1:alcsz)*term1(1:alcsz))
          else
            term2(1:alcsz) = term2(1:alcsz)*pfac2
            term0(1:alcsz) = term0(1:alcsz)*electric
            do i=1,nrpolnb(rs)
              if (abs(term0(i)).gt.abs(term2(i))) then
                term1(i) = term0(i)
              else
!               dissociated
                if ((termr(i)+par_IMPSOLV(1)).le.d1(i)) then
                  term1(i) = term0(i)
!               close contact
                else if (d1(i).lt.termr(i)) then
                  term1(i) = par_IMPSOLV(9)*(term2(i)-term0(i)) + &
 &                                      term0(i)
!               first-shell regime
                else
                  term1(i) = (termr(i)+par_IMPSOLV(1)-d1(i))/&
 &                           par_IMPSOLV(1)
                  term1(i)= (1.0-par_IMPSOLV(9)*term1(i))*term0(i)&
 &                 + par_IMPSOLV(9)*term1(i)*term2(i) 
                end if
              end if
            end do
            term1(1:alcsz) = term1(1:alcsz)*id1s(1:alcsz)
            evec(6) = evec(6) + par_FEG2(9)*sum(term1(1:alcsz))
          end if
        else
          term1(1:alcsz) = term0(1:alcsz)*id1s(1:alcsz)
          evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
        end if
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
!   code 2 indicates the residues are actually beyond the long cutoff-range, but are
!   included due to their charges
!   of course, the correctness of rsp_mat relies on cutoffs to be used -> check
    if (cut.EQV..true.) then
      if (rsp_mat(min(rs1,rs2),max(rs1,rs2)).eq.2) then
        call Ven_rsp_lrel_feg(evec,rs1,rs2,svec)
        return
      end if
    end if
!
    if (use_FEGS(6).EQV..true.) then
!
!     first set the necessary atom-specific parameters (mimimal)
      k = 0
      pfac = electric*par_FEG2(9)
      alcsz = at(rs1)%npol*at(rs2)%npol
      if (alcsz.gt.0) then
        if (use_FEGS(4).EQV..true.) then
          if (scrq_model.le.2) then
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)! + svec(1)
                dvec(k,2) = y(kk) - y(ii)! + svec(2)
                dvec(k,3) = z(kk) - z(ii)! + svec(3)
                term0(k) = atq(ii)*atq(kk)*scrq(ii)*scrq(kk)
              end do
            end do
          else if (scrq_model.eq.4) then
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)! + svec(1)
                dvec(k,2) = y(kk) - y(ii)! + svec(2)
                dvec(k,3) = z(kk) - z(ii)! + svec(3)
                term0(k) = atq(ii)*atq(kk)
                term2(k) = 1.0/(atr(ii)+atr(kk))
              end do
            end do
          else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)! + svec(1)
                dvec(k,2) = y(kk) - y(ii)! + svec(2)
                dvec(k,3) = z(kk) - z(ii)! + svec(3)
                term0(k) = atq(ii)*atq(kk)*scrq(ii)*scrq(kk)
                termr(k) = atr(ii)+atr(kk)
                term2(k) = atq(ii)*atq(kk)/termr(k)
              end do
            end do
          else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)! + svec(1)
                dvec(k,2) = y(kk) - y(ii)! + svec(2)
                dvec(k,3) = z(kk) - z(ii)! + svec(3)
                term0(k) = atq(ii)*atq(kk)*&
 &                  genmu(scrq(ii),scrq(kk),i_sqm)
              end do
            end do
          else ! 7,8
            pfac2 = par_IMPSOLV(8)*electric
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)! + svec(1)
                dvec(k,2) = y(kk) - y(ii)! + svec(2)
                dvec(k,3) = z(kk) - z(ii)! + svec(3)
                term0(k) = atq(ii)*atq(kk)*&
 &                    genmu(scrq(ii),scrq(kk),i_sqm)
                termr(k) = atr(ii)+atr(kk)
                term2(k) = atq(ii)*atq(kk)/termr(k)
              end do
            end do
          end if
        else
          do i=1,at(rs1)%npol
            ii = at(rs1)%pol(i)
            do j=1,at(rs2)%npol
              kk = at(rs2)%pol(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)! + svec(1)
              dvec(k,2) = y(kk) - y(ii)! + svec(2)
              dvec(k,3) = z(kk) - z(ii)! + svec(3)
              term0(k) = atq(ii)*atq(kk)
            end do
          end do
        end if
!       now perform the vectorizable bulk operations (maximal)
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id1s(1:alcsz) = 1.0/(d1(1:alcsz)+par_FEG2(10))
!       finish up -> screening model specificity is high
        if (use_FEGS(4).EQV..true.) then
          k = 0
          if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.&
 &            (scrq_model.eq.6)) then
            term1(1:alcsz) = term0(1:alcsz)*id1s(1:alcsz)
            evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
          else if (scrq_model.eq.4) then
            do i=1,at(rs1)%npol
              do j=1,at(rs2)%npol
                k = k + 1
                if (id1(k).le.term2(k)) then
                  term1(k) = term0(k)*id1(k)
                else
                  term1(k) = term0(k)*term2(k) 
                end if
              end do
            end do
            evec(6) = evec(6) + par_IMPSOLV(8)*pfac*&
 &                       sum(id1s(1:alcsz)*term1(1:alcsz))
          else
            term2(1:alcsz) = term2(1:alcsz)*pfac2
            term0(1:alcsz) = term0(1:alcsz)*electric
            do i=1,at(rs1)%npol
              do j=1,at(rs2)%npol
                k = k + 1
                if (abs(term0(k)).gt.abs(term2(k))) then
                  term1(k) = term0(k)
                else
!                 dissociated
                  if ((termr(k)+par_IMPSOLV(1)).le.d1(k)) then
                    term1(k) = term0(k)
!                 close contact
                  else if (d1(k).lt.termr(k)) then
                    term1(k) = par_IMPSOLV(9)*(term2(k)-term0(k)) + &
 &                                        term0(k)
!                 first-shell regime
                  else
                    term1(k) = (termr(k)+par_IMPSOLV(1)-d1(k))/&
 &                             par_IMPSOLV(1)
                    term1(k)= (1.0-par_IMPSOLV(9)*term1(k))*term0(k)&
 &                   + par_IMPSOLV(9)*term1(k)*term2(k) 
                  end if
                end if
              end do
            end do
            term1(1:alcsz) = term1(1:alcsz)*id1s(1:alcsz)
            evec(6) = evec(6) + par_FEG2(9)*sum(term1(1:alcsz))
          end if
        else
          term1(1:alcsz) = term0(1:alcsz)*id1s(1:alcsz)
          evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
        end if
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
! note that this routine does NOT support crosslink corrections (we assume two crosslinked
! residues are never so far apart as to trigger the long-range corrections)
! it is ONLY called if the entry in rsp_mat is in fact 2
!
subroutine Ven_rsp_lrel(evec,rs1,rs2,svec)
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
!
  implicit none
!
  integer rs1,rs2,i,j,k,ii,kk,alcsz
  integer g1,g2
  RTYPE evec(MAXENERGYTERMS)
  RTYPE genmu,svec(3),pfac,pfac2
  RTYPE d2(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE d1(at(rs1)%na*at(rs2)%na),id1(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE termr(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)

  if (use_POLAR.EQV..true.) then
    pfac = electric*scale_POLAR
    if (lrel_mc.eq.1) then
      k = 0
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
              term0(k) = atq(ii)*atq(kk)*scrq(ii)*scrq(kk)
            end do
            end do
          end do
          end do
        else if (scrq_model.eq.4) then
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
              term0(k) = atq(ii)*atq(kk)
              term2(k) = 1.0/(atr(ii)+atr(kk))
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
          pfac2 = par_IMPSOLV(8)*electric
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
              term0(k) = atq(ii)*atq(kk)*scrq(ii)*scrq(kk)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          pfac2 = par_IMPSOLV(8)*electric
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
              term0(k) = atq(ii)*atq(kk)*genmu(scrq(ii),scrq(kk),i_sqm)
            end do
            end do
          end do
          end do
        else ! 7,8
          pfac2 = par_IMPSOLV(8)*electric
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
              term0(k) = atq(ii)*atq(kk)*genmu(scrq(ii),scrq(kk),i_sqm)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
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
            term0(k) = atq(ii)*atq(kk)
          end do
          end do
        end do
        end do
      end if
      alcsz = k
    else if (lrel_mc.eq.2) then
      k = 0
      if (use_IMPSOLV.EQV..true.) then
        if (scrq_model.le.2) then
          do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle ! dipole-dipole or dipole-monopole
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle ! dipole-dipole or dipole-monopole
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*atq(ii)*atq(kk)
            end do
            end do
          end do
          end do
        else if (scrq_model.eq.4) then
          do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = atq(ii)*atq(kk) 
              term2(k) = 1.0/(atr(ii)+atr(kk))
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*atq(ii)*atq(kk)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*atq(ii)*atq(kk)
            end do
            end do
          end do
          end do
        else ! 7,8
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*atq(ii)*atq(kk)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
            end do
            end do
          end do
          end do
        end if
      else
        do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          do i=1,at(rs1)%dpgrp(g1)%nats
          ii = at(rs1)%dpgrp(g1)%ats(i)
          do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
            kk = at(rs2)%dpgrp(g2)%ats(j)
            k = k + 1
            dvec(k,1) = x(kk) - x(ii)
            dvec(k,2) = y(kk) - y(ii)
            dvec(k,3) = z(kk) - z(ii)
            term0(k) = atq(ii)*atq(kk)
            end do
          end do
          end do
        end do
      end if
      alcsz = k
!   reduced monopole-monopole
    else if (lrel_mc.eq.3) then
      k = 0
      if (use_IMPSOLV.EQV..true.) then
        if (scrq_model.le.2) then
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle ! dipole-dipole or dipole-monopole
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle ! dipole-dipole or dipole-monopole
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        else if (scrq_model.eq.4) then
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
              term2(k) = 1.0/(atr(ii)+atr(kk))
            end do
          end do
        else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = (1.0/termr(k))*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        else ! 7,8
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = (1.0/termr(k))*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        end if
      else
        do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
          do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
            k = k + 1
            dvec(k,1) = x(kk) - x(ii)
            dvec(k,2) = y(kk) - y(ii)
            dvec(k,3) = z(kk) - z(ii)
            term0(k) = cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
          end do
        end do
      end if
      alcsz = k
    end if ! which lrel_mc
!
!   now perform the vectorizable bulk operations (maximal)
    dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
    dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
    dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
    d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
    d1(1:alcsz) = sqrt(d2(1:alcsz))
    id1(1:alcsz) = 1.0/d1(1:alcsz)
!
!   finish up -> screening model specificity is high
    if (use_IMPSOLV.EQV..true.) then
      k = 0
      if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.(scrq_model.eq.6)) then
        term1(1:alcsz) = term0(1:alcsz)*id1(1:alcsz)
        evec(6) = evec(6) + pfac*sum(term1(1:alcsz))

      else if (scrq_model.eq.4) then
        do k=1,alcsz
          if (id1(k).le.term2(k)) then
            term1(k) = term0(k)*id1(k)
          else
            term1(k) = term0(k)*term2(k) 
          end if
        end do           
        evec(6) = evec(6) + par_IMPSOLV(8)*pfac*sum(id1(1:alcsz)*term1(1:alcsz))
      else   ! 3,7,8,9
        term2(1:alcsz) = term2(1:alcsz)*pfac2
        term0(1:alcsz) = term0(1:alcsz)*electric
        do k=1,alcsz
          if (abs(term0(k)).gt.abs(term2(k))) then
            term1(k) = term0(k)
          else
!        dissociated
            if ((termr(k)+par_IMPSOLV(1)).le.d1(k)) then
              term1(k) = term0(k)
!           close contact
            else if (d1(k).lt.termr(k)) then
              term1(k) = par_IMPSOLV(9)*(term2(k)-term0(k)) + term0(k)
!           first-shell regime
            else
              term1(k) = (termr(k)+par_IMPSOLV(1)-d1(k))/par_IMPSOLV(1)
              term1(k)= (1.0-par_IMPSOLV(9)*term1(k))*term0(k) + par_IMPSOLV(9)*term1(k)*term2(k) 
            end if
          end if
        end do
        term1(1:alcsz) = term1(1:alcsz)*id1(1:alcsz)
        evec(6) = evec(6) + scale_POLAR*sum(term1(1:alcsz))
      end if
    else
      term1(1:alcsz) = term0(1:alcsz)*id1(1:alcsz)
      evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
    end if
  end if  
end
!
!-----------------------------------------------------------------------
!
! see above
!
subroutine Ven_rsp_lrel_feg(evec,rs1,rs2,svec)
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
!
  implicit none
!
  integer rs1,rs2,i,j,k,ii,kk,alcsz
  integer g1,g2
  RTYPE evec(MAXENERGYTERMS)
  RTYPE genmu,svec(3),pfac,pfac2
  RTYPE d2(at(rs1)%npol*at(rs2)%npol)
  RTYPE dvec(at(rs1)%npol*at(rs2)%npol,3)
  RTYPE d1(at(rs1)%npol*at(rs2)%npol),id1(at(rs1)%npol*at(rs2)%npol)
  RTYPE id1s(at(rs1)%npol*at(rs2)%npol)
  RTYPE term0(at(rs1)%npol*at(rs2)%npol)
  RTYPE termr(at(rs1)%npol*at(rs2)%npol)
  RTYPE term1(at(rs1)%npol*at(rs2)%npol)
  RTYPE term2(at(rs1)%npol*at(rs2)%npol)

  if (use_FEGS(6).EQV..true.) then
!
    pfac = electric*par_FEG2(9)
    if (lrel_mc.eq.1) then
      k = 0
      if (use_FEGS(4).EQV..true.) then
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
              term0(k) = atq(ii)*atq(kk)*scrq(ii)*scrq(kk)
            end do
            end do
          end do
          end do
        else if (scrq_model.eq.4) then
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
              term0(k) = atq(ii)*atq(kk)
              term2(k) = 1.0/(atr(ii)+atr(kk))
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
          pfac2 = par_IMPSOLV(8)*electric
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
              term0(k) = atq(ii)*atq(kk)*scrq(ii)*scrq(kk)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          pfac2 = par_IMPSOLV(8)*electric
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
              term0(k) = atq(ii)*atq(kk)*genmu(scrq(ii),scrq(kk),i_sqm)
            end do
            end do
          end do
          end do
        else ! 7,8
          pfac2 = par_IMPSOLV(8)*electric
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
              term0(k) = atq(ii)*atq(kk)*genmu(scrq(ii),scrq(kk),i_sqm)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
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
            term0(k) = atq(ii)*atq(kk)
          end do
          end do
        end do
        end do
      end if
      alcsz = k
    else if (lrel_mc.eq.2) then
      k = 0
      if (use_FEGS(4).EQV..true.) then
        if (scrq_model.le.2) then
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*atq(ii)*atq(kk)
              end do
            end do
            end do
          end do
        else if (scrq_model.eq.4) then
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = atq(ii)*atq(kk) 
              term2(k) = 1.0/(atr(ii)+atr(kk))
              end do
            end do
            end do
          end do
        else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*atq(ii)*atq(kk)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
              end do
            end do
            end do
          end do
        else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*atq(ii)*atq(kk)
              end do
            end do
            end do
          end do
        else ! 7,8
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*atq(ii)*atq(kk)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
              end do
            end do
            end do
          end do
        end if
      else
        do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          do i=1,at(rs1)%dpgrp(g1)%nats
          ii = at(rs1)%dpgrp(g1)%ats(i)
          do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
            kk = at(rs2)%dpgrp(g2)%ats(j)
            k = k + 1
            dvec(k,1) = x(kk) - x(ii)
            dvec(k,2) = y(kk) - y(ii)
            dvec(k,3) = z(kk) - z(ii)
            term0(k) = atq(ii)*atq(kk)
            end do
          end do
          end do
        end do
      end if
      alcsz = k
    else if (lrel_mc.eq.3) then
      k = 0
      if (use_FEGS(4).EQV..true.) then
        if (scrq_model.le.2) then
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        else if (scrq_model.eq.4) then
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
              term2(k) = 1.0/(atr(ii)+atr(kk))
            end do
          end do
        else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = (1.0/termr(k))*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        else ! 7,8
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = (1.0/termr(k))*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        end if
      else
        do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
          do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
            k = k + 1
            dvec(k,1) = x(kk) - x(ii)
            dvec(k,2) = y(kk) - y(ii)
            dvec(k,3) = z(kk) - z(ii)
            term0(k) = cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
          end do
        end do
      end if
      alcsz = k
    end if ! which lrel_mc
!
!   now perform the vectorizable bulk operations (maximal)
    dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
    dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
    dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
    d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
    d1(1:alcsz) = sqrt(d2(1:alcsz))
    id1(1:alcsz) = 1.0/d1(1:alcsz)
    id1s(1:alcsz) = 1.0/(d1(1:alcsz)+par_FEG2(10))
!   finish up -> screening model specificity is high
    if (use_FEGS(4).EQV..true.) then
      k = 0
      if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.(scrq_model.eq.6)) then
        term1(1:alcsz) = term0(1:alcsz)*id1s(1:alcsz)
        evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
      else if (scrq_model.eq.4) then
        do k=1,alcsz
          if (id1(k).le.term2(k)) then
            term1(k) = term0(k)*id1(k)
          else
           term1(k) = term0(k)*term2(k) 
          end if
        end do           
        evec(6) = evec(6) + par_IMPSOLV(8)*pfac*sum(id1s(1:alcsz)*term1(1:alcsz))
      else    !  #3,7,8,9
        term2(1:alcsz) = term2(1:alcsz)*pfac2
        term0(1:alcsz) = term0(1:alcsz)*electric
        do k=1,alcsz
          if (abs(term0(k)).gt.abs(term2(k))) then
            term1(k) = term0(k)
          else
!        dissociated
            if ((termr(k)+par_IMPSOLV(1)).le.d1(k)) then
              term1(k) = term0(k)
!           close contact
            else if (d1(k).lt.termr(k)) then
              term1(k) = par_IMPSOLV(9)*(term2(k)-term0(k)) + term0(k)
!           first-shell regime
            else
              term1(k) = (termr(k)+par_IMPSOLV(1)-d1(k))/par_IMPSOLV(1)
              term1(k)= (1.0-par_IMPSOLV(9)*term1(k))*term0(k) + par_IMPSOLV(9)*term1(k)*term2(k) 
            end if
          end if
        end do
        term1(1:alcsz) = term1(1:alcsz)*id1s(1:alcsz)
        evec(6) = evec(6) + par_FEG2(9)*sum(term1(1:alcsz))
      end if
    else
      term1(1:alcsz) = term0(1:alcsz)*id1s(1:alcsz)
      evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
    end if
!
  end if
!       
end
!
!--------------------------------------------------------------------------------------------
!
