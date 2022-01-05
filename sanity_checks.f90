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
! a routine which recomputes the global system energy
! and calls a sanity fxn for smart cutoffs (i.e., get list of atoms in range ->
! check all, are there any which were forgotten?)
!
subroutine energy_check()
!
  use iounit
  use energies
  use atoms
  use cutoffs
  use accept
!
  implicit none
!
  RTYPE energy,tsv,esav2
  RTYPE energy_debug
  integer i
  
!
  if (use_IMPSOLV.EQV..true.) then
    tsv = sum(atsav(1:n))
  end if
  write(ilog,*) 
  write(ilog,*) 'Current full energy   :  ',esave
  write(ilog,*) 'Current IPP/LJ/WCA    :  ',esterms(1) + esterms(3)+&
 &esterms(5)
  write(ilog,*) 'Current TAB/POLAR     :  ',esterms(6)+esterms(9)
  write(ilog,*) 'Current Bonded Terms  :  ',esterms(2) + esterms(15)&
 &                         + esterms(16) + esterms(17) + esterms(18) + esterms(20)
  write(ilog,*) 'Current Bias Terms    :  ',esterms(7) + esterms(8)&
 &          + esterms(13) + esterms(10) + esterms(11) + esterms(14)&
 &          + esterms(19) + esterms(21)
  esav2 = esave
  esave = energy(esterms)
  write(ilog,*) 'Reset of full energy  :  ',esave
  write(ilog,*) 'Full N^2 IPP/LJ/WCA   :  ',esterms(1) + esterms(3)+&
 &esterms(5)
  write(ilog,*) 'Full N^2 TAB/POLAR    :  ',esterms(6)+esterms(9)
  write(ilog,*) 'Accurate Bonded Terms :  ',esterms(2) + esterms(15)&
 &                         + esterms(16) + esterms(17) + esterms(18) + esterms(20)
  write(ilog,*) 'Accurate Bias Terms   :  ',esterms(7) + esterms(8)&
 &          + esterms(13) + esterms(10) + esterms(11) + esterms(14)&
 &          + esterms(19) + esterms(21)
!
  if ((use_cutoffs.EQV..true.).AND.((use_mcgrid.EQV..true.).OR.(use_rescrit.EQV..true.))) then
    if (do_n2loop.EQV..true.) then
      write(ilog,*) 'Performing sanity check for cutoffs:'
      call cutoff_check()
      write(ilog,*) 'End sanity check for cutoffs.'
    end if
  end if
  if (use_IMPSOLV.EQV..true.) then
    write(ilog,*) 'Current total SAV:    ',tsv
    tsv = sum(atsav(1:n))
    write(ilog,*) 'Accurate total SAV:   ',tsv
  end if
! uncomment to diagnose energy discrepancies
!  if (abs(esav2 -esave).gt.1.0e-9) then
!    write(*,*) esav2,esave,abs(esav2 -esave)
!    call fexit()
!  end if

!
end
!
!-----------------------------------------------------------------------------------------------
!
! a routine to correct sampling-relevant arrays for unsupported residues and set some
! parameters unsettable before
! requires rotation lists to be complete and Z-matrix to be final
!
subroutine correct_sampler()
!
  use movesets
  use sequen
  use fyoc
  use polypep
  use zmatrix
  use atoms
  use iounit
  use molecule
  use forces
!
  implicit none
!
  integer rs,shf2,k,imol,shf
  type(t_entlst):: bulst
  logical renter
!
  do rs=1,nseq
    imol = molofrs(rs)
    if (seqtyp(rs).ne.26) cycle
!
    if (seqpolty(rs).eq.'P') then
      seqflag(rs) = 2
      shf2 = 0
      do k=at(rs)%bb(1),at(rs)%bb(1)+at(rs)%nbb+at(rs)%nsc-1
        if ((iz(1,k).eq.cai(rs)).AND.(atnam(k)(1:1).eq.'H').AND.(n12(k).eq.1)) shf2 = shf2 + 1
      end do
      if ((yline(rs).gt.0).AND.(fline(rs).gt.0)) then
        if ((izrot(yline(rs))%alsz.gt.0).AND.(izrot(fline(rs))%alsz.gt.0).AND.(shf2.eq.0)) seqflag(rs) = 8
        if ((izrot(yline(rs))%alsz.gt.0).AND.(izrot(fline(rs))%alsz.le.0)) seqflag(rs) = 5
        if (((izrot(yline(rs))%alsz.gt.0).AND.(izrot(fline(rs))%alsz.gt.0)).OR.(seqflag(rs).eq.5)) then
          if (allocated(bulst%idx).EQV..true.) deallocate(bulst%idx)
          if (fylst%nr.gt.0) then
            allocate(bulst%idx(fylst%nr+1))
            bulst%idx(1:fylst%nr) = fylst%idx(:)
            deallocate(fylst%idx)
            deallocate(fylst%wt)
          end if
          allocate(fylst%idx(fylst%nr+1))
          allocate(fylst%wt(fylst%nr+1))
          if (fylst%nr.gt.0) then
            fylst%idx(1:fylst%nr) = bulst%idx(1:fylst%nr)
          end if
          fylst%nr = fylst%nr + 1
          fylst%idx(fylst%nr) = rs
          renter = .false.
          do k=0,fylst%nr-2
            if (k.eq.0) then
              if (fylst%idx(k+1).gt.rs) renter = .true.
            else
              if ((fylst%idx(k).lt.rs).AND.(fylst%idx(k+1).gt.rs)) renter = .true.
            end if
            if (renter.EQV..true.) then
              bulst%idx(1:(fylst%nr-k-1)) = fylst%idx((k+1):(fylst%nr-1))
              fylst%idx(k+1) = rs
              fylst%idx((k+2):fylst%nr) = bulst%idx(1:(fylst%nr-k-1))
              exit
            end if
          end do
        end if
      end if
      if (pucline(rs).gt.0) then
        if (allocated(bulst%idx).EQV..true.) deallocate(bulst%idx)
        if (puclst%nr.gt.0) then
          allocate(bulst%idx(puclst%nr+1))
          bulst%idx(1:puclst%nr) = puclst%idx(:)
          deallocate(puclst%idx)
          deallocate(puclst%wt)
        end if
        allocate(puclst%idx(puclst%nr+1))
        allocate(puclst%wt(puclst%nr+1))
        if (puclst%nr.gt.0) then
          puclst%idx(1:puclst%nr) = bulst%idx(1:puclst%nr)
        end if
        puclst%nr = puclst%nr + 1
        puclst%idx(puclst%nr) = rs
        renter = .false.
        do k=0,puclst%nr-2
          if (k.eq.0) then
            if (puclst%idx(k+1).gt.rs) renter = .true.
          else
            if ((puclst%idx(k).lt.rs).AND.(puclst%idx(k+1).gt.rs)) renter = .true.
          end if
          if (renter.EQV..true.) then
            bulst%idx(1:(puclst%nr-k-1)) = puclst%idx((k+1):(puclst%nr-1))
            puclst%idx(k+1) = rs
            puclst%idx((k+2):puclst%nr) = bulst%idx(1:(puclst%nr-k-1))
            exit
          end if
        end do
      end if
      if (wline(rs).gt.0) then
        if (izrot(wline(rs))%alsz.gt.0) then
          if (allocated(bulst%idx).EQV..true.) deallocate(bulst%idx)
          if (wlst%nr.gt.0) then
            allocate(bulst%idx(wlst%nr+1))
            bulst%idx(1:wlst%nr) = wlst%idx(:)
            deallocate(wlst%idx)
            deallocate(wlst%wt)
          end if
          allocate(wlst%idx(wlst%nr+1))
          allocate(wlst%wt(wlst%nr+1))
          if (wlst%nr.gt.0) then
            wlst%idx(1:wlst%nr) = bulst%idx(1:wlst%nr)
          end if
          wlst%nr = wlst%nr + 1
          wlst%idx(wlst%nr) = rs
          renter = .false.
          do k=0,wlst%nr-2
            if (k.eq.0) then
              if (wlst%idx(k+1).gt.rs) renter = .true.
            else
              if ((wlst%idx(k).lt.rs).AND.(wlst%idx(k+1).gt.rs)) renter = .true.
            end if
            if (renter.EQV..true.) then
              bulst%idx(1:(wlst%nr-k-1)) = wlst%idx((k+1):(wlst%nr-1))
              wlst%idx(k+1) = rs
              wlst%idx((k+2):wlst%nr) = bulst%idx(1:(wlst%nr-k-1))
              exit
            end if
          end do
        end if
      end if
    else if (seqpolty(rs).eq.'N') then
      seqflag(rs) = 22
      if ((nuci(rs,4).gt.0).AND.(nuci(rs,5).le.0)) seqflag(rs) = 24
      if (allocated(bulst%idx).EQV..true.) deallocate(bulst%idx)
      if (nuclst%nr.gt.0) then
        allocate(bulst%idx(nuclst%nr+1))
        bulst%idx(1:nuclst%nr) = nuclst%idx(:)
        deallocate(nuclst%idx)
        deallocate(nuclst%wt)
      end if
      allocate(nuclst%idx(nuclst%nr+1))
      allocate(nuclst%wt(nuclst%nr+1))
      if (nuclst%nr.gt.0) then
        nuclst%idx(1:nuclst%nr) = bulst%idx(1:nuclst%nr)
      end if
      nuclst%nr = nuclst%nr + 1
      nuclst%idx(nuclst%nr) = rs
      renter = .false.
      do k=0,nuclst%nr-2
        if (k.eq.0) then
          if (nuclst%idx(k+1).gt.rs) renter = .true.
        else
          if ((nuclst%idx(k).lt.rs).AND.(nuclst%idx(k+1).gt.rs)) renter = .true.
        end if
        if (renter.EQV..true.) then
          bulst%idx(1:(nuclst%nr-k-1)) = nuclst%idx((k+1):(nuclst%nr-1))
          nuclst%idx(k+1) = rs
          nuclst%idx((k+2):nuclst%nr) = bulst%idx(1:(nuclst%nr-k-1))
          exit
        end if
      end do
      if (nucsline(6,rs).gt.0) then
        if (allocated(bulst%idx).EQV..true.) deallocate(bulst%idx)
        if (nucpuclst%nr.gt.0) then
          allocate(bulst%idx(nucpuclst%nr+1))
          bulst%idx(1:nucpuclst%nr) = nucpuclst%idx(:)
          deallocate(nucpuclst%idx)
          deallocate(nucpuclst%wt)
        end if
        allocate(nucpuclst%idx(nucpuclst%nr+1))
        allocate(nucpuclst%wt(nucpuclst%nr+1))
        if (nucpuclst%nr.gt.0) then
          nucpuclst%idx(1:nucpuclst%nr) = bulst%idx(1:nucpuclst%nr)
        end if
        nucpuclst%nr = nucpuclst%nr + 1
        nucpuclst%idx(nucpuclst%nr) = rs
        renter = .false.
        do k=0,nucpuclst%nr-2
          if (k.eq.0) then
            if (nucpuclst%idx(k+1).gt.rs) renter = .true.
          else
            if ((nucpuclst%idx(k).lt.rs).AND.(nucpuclst%idx(k+1).gt.rs)) renter = .true.
          end if
          if (renter.EQV..true.) then
            bulst%idx(1:(nucpuclst%nr-k-1)) = nucpuclst%idx((k+1):(nucpuclst%nr-1))
            nucpuclst%idx(k+1) = rs
            nucpuclst%idx((k+2):nucpuclst%nr) = bulst%idx(1:(nucpuclst%nr-k-1))
            exit
          end if
        end do
      end if
    else if ((rs.ne.rsmol(imol,1)).AND.(rs.ne.rsmol(imol,2))) then
      seqflag(rs) = 50
    else if ((rs.eq.rsmol(imol,1)).AND.(rs.ne.rsmol(imol,2))) then
      if ((wline(rs+1).gt.0).AND.(seqpolty(rs+1).eq.'P')) then
        seqflag(rs) = 10
      else if ((wline(rs+1).le.0).AND.(seqpolty(rs+1).eq.'P')) then
        seqflag(rs) = 12
      else if (seqpolty(rs+1).eq.'N') then
        seqflag(rs) = 26
      else
        seqflag(rs) = 50
      end if
    else if ((rs.eq.rsmol(imol,2)).AND.(rs.ne.rsmol(imol,1))) then
      if ((wline(rs).gt.0).AND.(seqpolty(rs-1).eq.'P')) then
        seqflag(rs) = 11
      else if ((wline(rs).le.0).AND.(seqpolty(rs-1).eq.'P')) then
        seqflag(rs) = 13
      else if (seqpolty(rs-1).eq.'N') then
        seqflag(rs) = 28
      else
        seqflag(rs) = 50
      end if
    else if (rsmol(imol,1).eq.rsmol(imol,2)) then
      shf = 0
      shf2 = 0
      do k=atmol(imol,1)+1,atmol(imol,2)
        if (n14(k).gt.0) shf = shf + 1
        if (izrot(k)%alsz.gt.0) shf2 = shf2 + 1
      end do
      seqflag(rs) = 101
      if ((shf.gt.0).AND.(shf2.gt.0)) then
        seqflag(rs) = 103
      else if (shf.gt.0) then
        seqflag(rs) = 102
      end if
    end if
  end do
!
  do imol=1,nmol
    if ((ntormol(moltypid(imol)).gt.0).AND.(size(dc_di(imol)%frz(:)).ne.(ntormol(moltypid(imol))+6))) then
      write(ilog,*) 'Fatal. Setup of degrees of freedom indicates a bug in correct_sampler() for molecule ',imol,'. This &
 &is fatal.'
      call fexit()
    end if
  end do
!
!  write(*,*) 'ST',seqtyp(1:nseq),'SF',seqflag(1:nseq)
!  write(*,*) 'F',fline(1:nseq)
!  write(*,*) 'F2',fline2(1:nseq)
!  write(*,*) 'Y',yline(1:nseq)
!  write(*,*) 'Y2',yline2(1:nseq)
!  write(*,*) 'P',pucline(1:nseq)
!  write(*,*) 'W',wline(1:nseq)

!
end
!
!----------------------------------------------------------------------------
!
! fix potential problems with initial structure requests
!
subroutine strucinp_sanitychecks()
!
  use system
  use pdb
  use molecule
  use sequen
  use iounit
!
  implicit none
!
  if (pdb_analyze.EQV..true.) then
    globrandomize = 0 ! this would be utterly pointless
    fycinput = .false.
    if ((pdb_fileformat.ne.3).AND.(pdb_fileformat.ne.4).AND.(pdb_fileformat.ne.5).AND.(pdbinput.EQV..false.)) then
      write(ilog,*) 'Fatal. Input file for pdb-input in analysis mod&
 &e could not be opened. Check spec. for FMCSC_PDBFILE.'
      call fexit()
    else
      pdbinput = .false.
    end if
  else if (do_restart.EQV..true.) then
    pdbinput = .false.
    globrandomize = 0
    fycinput = .false.
  else
    if ((nmol.gt.1).AND.(globrandomize.eq.0).AND.(pdbinput.EQV..false.)) then
      write(ilog,*) 'Warning. Randomization of rigid-body coordinates is necessary for sy&
 &stems with multiple molecules and will be performed (alternatively, read in from PDB).'
      globrandomize = 2
    end if
    if ((nmol.gt.1).AND.(fycinput.EQV..true.)) then
      write(ilog,*) 'Warning. Requested fyc-read-in for a system wit&
 &h more than one molecule, which is not supported. Disabled.'
      fycinput = .false.
      if (pdbinput.EQV..false.) globrandomize = 1
    end if
    if ((fycinput.EQV..true.).AND.(pdbinput.EQV..true.)) then
      write(ilog,*) 'WARNING: Requested read-in of structure both from&
 & a torsional input file and a PDB-file. Defaulting to PDB.'
      fycinput = .false.
    end if
    if ((fycinput.EQV..true.).AND.((globrandomize.eq.1).OR.(globrandomize.eq.3))) then
      write(ilog,*) 'WARNING: Requested both read-in of a torsional input file and internal r&
 &andomization. Changing randomization settings to rigid-body only.'
      globrandomize = 2
    end if
    if ((globrandomize.eq.1).AND.(pdbinput.EQV..true.)) then
      write(ilog,*) 'Warning. Requested full randomization and read-in of PDB-file. Changing &
 &randomization settings to rigid-body and missing parts only.'
      globrandomize = 4
    end if
    if ((globrandomize.eq.3).AND.(pdbinput.EQV..false.)) then
      write(ilog,*) 'Warning. The requested randomization of just dihedral angles requires &
 &a PDB file to be read in for setting rigid-body coordinates. Defaulting to full ra&
 &ndomization.'
      globrandomize = 1
    end if
    if ((globrandomize.eq.0).AND.(nmol.gt.1).AND.(pdbinput.EQV..false.)) then
      write(ilog,*) 'Warning. Rigid-body randomization is mandatory for systems with multiple &
 &molecules and no external structure input. Enabled.'
      globrandomize = 2
    end if
  end if
!
end
!
!------------------------------------------------------------------
!
! note that this function must be called BEFORE user-level constraints are applied 
!
subroutine moveset_sanitychecks()
!
  use iounit
  use molecule
  use movesets
  use sequen
  use fyoc
  use torsn
  use grandensembles
  use system
  use ujglobals
  use mpistuff
  use aminos
  use wl
  use params
  use atoms
  use zmatrix
  use polypep
  use cutoffs
!
  implicit none
!
  integer dummy
  integer, ALLOCATABLE:: rslst(:),rslst2(:)
!
! first check for unsupported linkage types
!
  allocate(rslst(nseq))
  rslst(:) = 0
  do dummy=1,nseq
    if (dummy.eq.rsmol(molofrs(dummy),1)) cycle
    rslst(dummy) = atmres(iz(1,at(dummy)%bb(1)))
  end do
  do dummy=1,nseq
    if ((rslst(dummy).gt.0).AND.(rslst(dummy).ne.(dummy-1))) then
      if (pdb_analyze.EQV..false.) then
        write(ilog,*) 'Fatal. Branched polymers (i.e., a side chain (branch) using separate residues) are currently &
 &not supported in Monte Carlo runs. Lump all atoms in the shortest branch into the corresponding main chain residue &
 &to circumvent this error.'
        call fexit()
      end if
    end if
  end do
!
! first set up the residue lists for concerted rotation type moves
!
  allocate(rslst2(nseq))
  dummy = 0
  call setup_CR_eligible(dummy,rslst,rslst2)
  deallocate(rslst)
  deallocate(rslst2)
!
! a sanity check against redundant rigid-body moves
  if ((nmol.eq.1).AND.(ntorpuck.gt.0).AND.(rigidfreq.gt.0.0)) then
    write(ilog,*) 'Warning. Rigid-body moves are redundant with only one molecule in the system. Consider turning them off.'
  end if
!
  if ((particleflucfreq.gt.0).AND.(ens%flag.ne.5).AND.(ens%flag.ne.6)) then
    write(ilog,*) 'WARNING: Particle fluctuations moves are enabled,&
 & but neither the grand or semigrand ensembles are in use. Turning off.'
    particleflucfreq = 0.0
  end if
!
! a sanity check against internal moves
  if ((ntorpuck.eq.0).AND.(unklst%nr.eq.0).AND.(unslst%nr.eq.0).AND.(rigidfreq.lt.1.0).AND.(particleflucfreq.lt.1.0)) then
    write(ilog,*) 'WARNING: No internal degrees of freedom to sample. Turning off all corresponding moves.'
    rigidfreq = 1.0
    omegafreq = 0.0
    chifreq = 0.0
    crfreq = 0.0
    nucfreq = 0.0
    nuccrfreq = 0.0
    angcrfreq = 0.0
    torcrfreq = 0.0
    puckerfreq = 0.0
    nucpuckfreq = 0.0
    otherfreq = 0.0
  end if
!
! a sanity check against nonpucker internal moves
  if ((ntorsn.eq.0).AND.(unklst%nr.eq.0).AND.(unslst%nr.eq.0).AND.(rigidfreq.lt.1.0).AND.(particleflucfreq.lt.1.0)) then
    write(ilog,*) 'WARNING: No freely rotatable dihedral angles to sample. Turning off all corresponding moves.'
    omegafreq = 0.0
    chifreq = 0.0
    crfreq = 0.0
    nucfreq = 0.0
    nuccrfreq = 0.0
    angcrfreq = 0.0
    torcrfreq = 0.0
    otherfreq = 0.0
    puckerfreq = 1.0
    nucpuckfreq = 1.0
  end if
!
! a sanity check against backbone moves 
  dummy = sum(nchi(1:nseq))
  if ((rigidfreq.lt.1.0).AND.(particleflucfreq.lt.1.0).AND.(dummy.eq.ntorpuck).AND.(dummy.gt.0).AND.(otherfreq.le.0.0)) then
    chifreq = 1.0
    nucfreq = 0.0
    nuccrfreq = 0.0
    crfreq = 0.0
    angcrfreq = 0.0
    torcrfreq = 0.0
    omegafreq = 0.0
    puckerfreq = 0.0
  end if
!
! another one
  if ((rigidfreq.lt.1.0).AND.(particleflucfreq.lt.1.0).AND.(wlst%nr.gt.0).AND.(fylst%nr.eq.0).AND.&
 &    (chifreq.lt.1.0).AND.(nuclst%nr.eq.0).AND.(puclst%nr.eq.0).AND.(otherfreq.le.0.0)) then
    omegafreq = 1.0
    crfreq = 0.0
    angcrfreq = 0.0
    torcrfreq = 0.0
    puckerfreq = 0.0
    nuccrfreq = 0.0
    nucfreq = 0.0
    puckerfreq = 0.0
  end if
!
! another one
  if ((rigidfreq.lt.1.0).AND.(particleflucfreq.lt.1.0).AND.(wlst%nr.eq.0).AND.(fylst%nr.eq.0).AND.&
 &    (chifreq.lt.1.0).AND.(nuclst%nr.gt.0).AND.(otherfreq.le.0.0)) then
    omegafreq = 0.0
    crfreq = 0.0
    angcrfreq = 0.0
    torcrfreq = 0.0
    puckerfreq = 0.0
    nucfreq = 1.0
  end if
!
! another one - against alignment
  if ((align_NC.ne.1).AND.(nmol.eq.1).AND.(dyn_mode.eq.1)) then
    write(ilog,*) 'Warning. Within an MC simulation using a single molecular reference frame, it is probably inefficient &
 &to do (partial or rigorous) C-terminal alignment. Consider using N-terminal alignment instead.'
  end if
!
! a sanity check against RB and particule fluctuation rigid-body moves
  have_particlefluc = .true.
  if (particleflucfreq.le.0.0) have_particlefluc = .false.
  have_rigid = .true.
  if ((particleflucfreq.ge.1.0).AND.(rigidfreq.gt.0.0)) then
    write(ilog,*) 'WARNING: Particle insertion, transmutation, and deletion moves cannot alter the &
 &coordinates of molecules. If used exclusively, sampling is severely non-ergodic. Unless this is meant to &
 &emulate discretized space, results are physically meaningless.'
    write(ilog,*)
    have_rigid = .false.
  end if
  if (rigidfreq.le.0.0) have_rigid = .false.
!
! sanity checks against cluster RB moves
  have_clurb = .true.
  if ((have_rigid.EQV..true.).AND.(clurb_freq.gt.0.0).AND.(use_coupledrigid.EQV..true.)) then
    if (nmol.eq.1) then
      write(ilog,*) 'WARNING: Cluster rigid-body moves are impossibl&
 &e with only a single molecule in the system. Turning off.'
      clurb_freq = 0.0  
    else if (nmol.eq.2) then
      write(ilog,*) 'WARNING: Cluster rigid-body moves are possibly &
 & redundant with only two molecules in the system. Consider turning&
 & them off.'
    end if
    if (clurb_rdfreq.gt.0.0) then
      if (clurb_maxsz.gt.nmol) then
        clurb_maxsz = nmol
        write(ilog,*) 'WARNING: Max. size for random cluster RB move&
 &s exceeds number of molecules in the system. Lowering to ',&
 &clurb_maxsz,'.'
      end if
    end if
    if (clurb_rdfreq.lt.1.0) then
      write(ilog,*) 'Non-random cluster moves are currently not supp&
 &orte and will therefore be disabled. Please check back later.'
      clurb_rdfreq = 1.0
    end if
  else if ((have_rigid.EQV..true.).AND.(clurb_freq.gt.0.0).AND.(use_coupledrigid.EQV..false.)) then
    write(ilog,*) 'WARNING: Cluster rigid-body moves are only supported if FMCSC_COUPLERIGID is true. Turning off.'
    clurb_freq = 0.0  
  end if
  if ((have_rigid.EQV..false.).OR.(clurb_freq.le.0.0)) have_clurb = .false.
  if (have_clurb.EQV..true.) then
    allocate(clurblst%idx(rblst%nr))
    clurblst%idx(:) = rblst%idx
    clurblst%nr = rblst%nr
  end if
!
! a sanity check against sidechain moves
  have_chi = .true.
  if ((chifreq.gt.0.0).AND.(chilst%nr.eq.0)) then
    if ((rigidfreq.lt.1.0).AND.(particleflucfreq.lt.1.0)) then
      write(ilog,*) 'WARNING: No sidechains to sample! Adjusting frequency to 0.0.'
      write(ilog,*)
    end if
    chifreq = 0.0
  end if
  if ((chifreq.le.0.0).OR.(rigidfreq.ge.1.0).OR.(particleflucfreq.ge.1.0)) have_chi = .false.
!
! whether we have pH pseudo-moves
  have_ph = .false.
  if ((have_chi.EQV..true.).AND.(phfreq.gt.0.0)) have_ph = .true.
!
! sanity checks against non-NUC concerted rotation moves
  have_cr = .true.
  if ((crfreq.gt.0.0).AND.(ujcrlst%nr.eq.0).AND.(djcrlst%nr.eq.0).AND.(docrlst%nr.eq.0).AND.(sjcrlst%nr.eq.0)) then
    if ((rigidfreq.lt.1.0).AND.(particleflucfreq.lt.1.0).AND.(chifreq.lt.1.0)) then
      write(ilog,*) 'WARNING: No eligible polypeptide stretches for concerted rotation sampling &
 &of any kind! Adjusting frequency to 0.0.'
      write(ilog,*)
    end if
    crfreq = 0.0
  end if
  if ((crfreq.le.0.0).OR.(chifreq.ge.1.0).OR.(rigidfreq.ge.1.0).OR.(particleflucfreq.ge.1.0)) have_cr = .false.
  if (have_cr.EQV..true.) then
    have_torcr = .true.
!   first bond angle UJ method
    have_ujcr = .true.
    if ((have_cr.EQV..true.).AND.(ujcrlst%nr.eq.0)) then
      if (angcrfreq.gt.0.0) then
        write(ilog,*) 'WARNING: No suitable polypeptide chains for Ulmschneider-Jorgensen concerted rotation. &
   &Current settings require consecutive stretches of non-terminal residues that are not proline (or similar), which &
   &are at least ',ujminsz,' residues long.'
        write(ilog,*)
      end if
      angcrfreq = 0.0
    end if
    if ((have_cr.EQV..false.).OR.(angcrfreq.le.0.0)) have_ujcr = .false.
    if (have_ujcr.EQV..true.) then
      if (ujmaxsz.lt.ujminsz) then
        write(ilog,*) 'Warning. Maximum pre-rotation segment length &
   &in Ulmschneider-Jorgensen concerted rotation moves is shorter than minimum one. Adjusting to the same value.'
        ujmaxsz = ujminsz
      end if
      if (angcrfreq.ge.1.0) have_torcr = .false.
    end if
!   see if we have any torsional variants available, first DJ
    if ((djcrlst%nr.le.0).AND.(docrlst%nr.le.0).AND.(sjcrlst%nr.le.0)) have_torcr = .false.
    have_djcr = .true.
    if ((have_cr.EQV..true.).AND.(have_torcr.EQV..true.).AND.(djcrlst%nr.eq.0)) then
      if ((torcrfreq_omega.lt.1.0).AND.(torcrfreq.gt.0.0)) then
        write(ilog,*) 'WARNING: No suitable chains for Dinner-Ulmschneider polypeptide concerted rotation moves &
   &without omega sampling! Current settings require consecutive stretches of polypeptide residues that are at &
   &least 4 residues long (only N-terminal caps count).'
        write(ilog,*)
      end if
      if (torcrfreq.gt.0.0) then
        write(ilog,*) 'Tentatively turning requested Dinner-Ulmschneider polypeptide concerted rotation moves without &
   &omega sampling into those with omega sampling.'
        torcrfreq_omega = 1.0
      end if
    end if
    if ((have_torcr.EQV..false.).OR.(torcrfreq.le.0.0).OR.(torcrfreq_omega.ge.1.0)) have_djcr = .false.
    if (have_djcr.EQV..true.) then
      if (torcrmaxsz2.lt.torcrminsz2) then
        write(ilog,*) 'Warning. Maximum pre-rotation segment length &
   &in Dinner-Ulmschneider CR moves without omega sampling is shorter than preferred minimum one.'
        write(ilog,*) 'Adjusting to the same value.'
        torcrmaxsz2 = torcrminsz2
      end if
    end if
!   now DO
    have_docr = .true.
    if ((have_cr.EQV..true.).AND.(have_torcr.EQV..true.).AND.(docrlst%nr.eq.0)) then
      if ((torcrfreq_omega.gt.0.0).AND.(torcrfreq.gt.0.0)) then
        write(ilog,*) 'WARNING: No suitable chains for Dinner-Ulmschneider polypeptide concerted rotation moves &
   &with omega sampling! Current settings require consecutive stretches of polypeptide residues that are at &
   &least 3 residues long (caps are not eligible).'
        write(ilog,*)
      end if
      if (torcrfreq.gt.0.0) then
        write(ilog,*) 'Disabling all Dinner-Ulmschneider polypeptide concerted rotation moves.'
        torcrfreq = 0.0
      end if
    end if
    if ((have_torcr.EQV..false.).OR.(torcrfreq.le.0.0).OR.(torcrfreq_omega.le.0.0)) have_docr = .false.
    if (have_docr.EQV..true.) then
      if (torcrmaxsz.lt.torcrminsz) then
        write(ilog,*) 'Warning. Maximum pre-rotation segment length &
   &in Dinner-Ulmschneider CR moves with omega sampling is shorter than preferred minimum one.'
        write(ilog,*) 'Adjusting to the same value.'
        torcrmaxsz = torcrminsz
      end if
    end if
!   and finally SJ
    have_sjcr = .true.
    if ((have_torcr.EQV..true.).AND.(sjcrlst%nr.eq.0)) then
      if (torcrfreq.lt.1.0) then
        if (cr_mode.eq.1) then
          write(ilog,*) 'WARNING: No suitable polypeptide chains for Sjunnesson-Favrin concerted rot&
   &ation. Current settings require consecutive stretches of non-terminal residues that are not proline (or similar), which &
   &are at least ',nr_crres+2,' residues long.'
        else
          write(ilog,*) 'WARNING: No suitable polypeptide chains for Sjunnesson-Favrin concerted rot&
   &ation. Current settings require consecutive stretches of non-terminal residues that are not proline (or similar), which &
   &are at least 4 residues long.'
        end if
        write(ilog,*)
      end if
      torcrfreq = 1.0
    end if
    if ((have_cr.EQV..false.).OR.(have_torcr.EQV..false.).OR.(torcrfreq.ge.1.0)) have_sjcr = .false.
  else
    have_torcr = .false.
    have_ujcr = .false.
    have_djcr = .false.
    have_docr = .false.
    have_sjcr = .false.
  end if ! if (have_cr is true)
!
! a sanity check against omega moves
  have_omega = .true.
  if ((omegafreq.gt.0.0).AND.(wlst%nr.eq.0)) then
    if ((particleflucfreq.lt.1.0).AND.(rigidfreq.lt.1.0).AND.(chifreq.lt.1.0).AND.(crfreq.lt.1.0)) then
      write(ilog,*) 'WARNING: No omega angles to sample! Adjusting frequency to 0.0.'
      write(ilog,*)
    end if
    omegafreq = 0.0
  end if
  if ((omegafreq.le.0.0).OR.(rigidfreq.ge.1.0).OR.(particleflucfreq.ge.1.0).OR.(chifreq.ge.1.0).OR.&
 &    (crfreq.ge.1.0)) have_omega = .false.
!
! a sanity check against nuc.-moves
  have_nuc = .true.
  if ((nucfreq.gt.0.0).AND.(nuclst%nr.eq.0)) then
    if ((particleflucfreq.lt.1.0).AND.(rigidfreq.lt.1.0).AND.(chifreq.lt.1.0).AND.(omegafreq.lt.1.0).AND.(crfreq.lt.1.0)) then
      write(ilog,*) 'WARNING: No NUC-type backbone angles to sample! Adjusting frequency to 0.0.'
      write(ilog,*)
    end if
    nucfreq = 0.0
  end if
  if ((nucfreq.le.0.0).OR.(rigidfreq.ge.1.0).OR.(particleflucfreq.ge.1.0).OR.(chifreq.ge.1.0).OR.&
 &    (crfreq.ge.1.0).OR.(omegafreq.ge.1.0)) have_nuc = .false.
!
! a sanity check against NUCCR moves
  have_nuccr = .true.
  if ((have_nuc.EQV..true.).AND.(nuccrlst%nr.eq.0)) then
    if (nuccrfreq.gt.0.0) then
      write(ilog,*) 'WARNING: No suitable chains for nucleic acid concerted rotation moves to sample! This move type &
 &requires chains with 2 or 3 consecutive residues (depending on capping).'
      write(ilog,*)
    end if
    nuccrfreq = 0.0
  end if
  if ((have_nuc.EQV..false.).OR.(nuccrfreq.le.0.0)) have_nuccr = .false.
  if (have_nuccr.EQV..true.) then
    if (nuccrmaxsz.lt.nuccrminsz) then
      write(ilog,*) 'Warning. Maximum pre-rotation segment length &
 &in nucleic acid concerted rotation moves is shorter than minimum one. Adjusting to the same value.'
      nuccrmaxsz = nuccrminsz
    end if
  end if
!
! a sanity check against nuc. pucker moves
  have_nucpuck = .true.
  if ((have_nuc.EQV..true.).AND.(nuccrfreq.lt.1.0).AND.(nucpuclst%nr.eq.0)) then
    if (nucpuckfreq.gt.0.0) then
      write(ilog,*) 'WARNING: No nucleic acid pucker degrees of freedom to sample! Adjusting frequency to 0.0.'
      write(ilog,*)
    end if
    nucpuckfreq = 0.0
  end if
  if ((have_nuc.EQV..false.).OR.(nuccrfreq.ge.1.0).OR.(nucpuckfreq.le.0.0)) have_nucpuck = .false.
!
! a sanity check against pucker moves
  have_pucker = .true.
  if ((puckerfreq.gt.0.0).AND.(puclst%nr.eq.0)) then
    if ((rigidfreq.lt.1.0).AND.(chifreq.lt.1.0).AND.(nucfreq.lt.1.0).AND.&
 &      (crfreq.lt.1.0).AND.(omegafreq.lt.1.0)) then
      write(ilog,*) 'WARNING: No polypeptide pucker angles to sample! Adjusting f&
 &requency to 0.0.'
      write(ilog,*)
    end if
    puckerfreq = 0.0
  end if
  if ((puckerfreq.le.0.0).OR.(rigidfreq.ge.1.0).OR.(particleflucfreq.ge.1.0).OR.(chifreq.ge.1.0).OR.&
 &    (crfreq.ge.1.0).OR.(omegafreq.ge.1.0).OR.(nucfreq.ge.1.0)) have_pucker = .false.
!
! a sanity check against single dihedral moves
  have_other = .true.
  have_unkother = .true.
  have_natother = .true.
  have_unsother = .true.
  if ((otherfreq.gt.0.0).AND.(unklst%nr.eq.0).AND.(unslst%nr.eq.0).AND.(natlst%nr.eq.0)) then
    if ((rigidfreq.lt.1.0).AND.(chifreq.lt.1.0).AND.(nucfreq.lt.1.0).AND.&
 &      (crfreq.lt.1.0).AND.(omegafreq.lt.1.0).AND.(puckerfreq.lt.1.0)) then
      write(ilog,*) 'WARNING: No freely rotatable dihedral angles to sample! Adjusting f&
 &requency to 0.0.'
      write(ilog,*)
    end if
    otherfreq = 0.0
  else if (otherfreq.gt.0.0) then
    if ((other_unkfreq.gt.0.0).AND.(unklst%nr.eq.0)) then
      if ((rigidfreq.lt.1.0).AND.(chifreq.lt.1.0).AND.(nucfreq.lt.1.0).AND.&
 &        (crfreq.lt.1.0).AND.(omegafreq.lt.1.0).AND.(puckerfreq.lt.1.0)) then
        write(ilog,*) 'WARNING: No freely rotatable dihedral angles in unsupported residues! Adjusting f&
 &requency (FMCSC_OTHERUNKFREQ) to 0.0.'
        write(ilog,*)
      end if
      other_unkfreq = 0.0
    end if
    if ((other_natfreq.gt.0.0).AND.(natlst%nr.eq.0)) then
      if ((rigidfreq.lt.1.0).AND.(chifreq.lt.1.0).AND.(nucfreq.lt.1.0).AND.&
 &        (crfreq.lt.1.0).AND.(omegafreq.lt.1.0).AND.(puckerfreq.lt.1.0).AND.(other_unkfreq.lt.1.0)) then
        write(ilog,*) 'WARNING: No freely rotatable dihedral angles supported natively! Adjusting f&
 &requency (FMCSC_OTHERNATFREQ) to 0.0.'
        write(ilog,*)
      end if
      other_natfreq = 0.0
    end if
    if ((other_natfreq.lt.1.0).AND.(other_unkfreq.lt.1.0).AND.(unslst%nr.eq.0)) then
      if (other_natfreq.gt.0.0) then
        other_natfreq = 1.0
      else if (other_unkfreq.gt.0.0) then
        other_unkfreq = 1.0
      else
        write(ilog,*) 'WARNING: No freely rotatable dihedral angles not supported natively! Adjusting f&
 &requency for all corresponding moves (FMCSC_OTHERFREQ) to 0.0.'
        otherfreq = 0.0
      end if
    end if
  end if
!
  if ((otherfreq.le.0.0).OR.(rigidfreq.ge.1.0).OR.(particleflucfreq.ge.1.0).OR.(chifreq.ge.1.0).OR.&
 &    (crfreq.ge.1.0).OR.(omegafreq.ge.1.0).OR.(nucfreq.ge.1.0)) have_other = .false.
  if ((have_other.EQV..false.).OR.(other_unkfreq.le.0.0)) have_unkother = .false.
  if ((have_other.EQV..false.).OR.(other_unkfreq.ge.1.0).OR.(other_natfreq.le.0.0)) have_natother = .false.
  if ((have_other.EQV..false.).OR.(other_unkfreq.ge.1.0).OR.(other_natfreq.ge.1.0)) have_unsother = .false.
  if ((have_unsother.EQV..false.).AND.(have_natother.EQV..false.).AND.(have_unkother.EQV..false.)) have_other = .false.
!
! a sanity check against pivot moves
  have_pivot = .true.
  if ((otherfreq.lt.1.0).AND.(fylst%nr.eq.0)) then
    if ((rigidfreq.lt.1.0).AND.(chifreq.lt.1.0).AND.(nucfreq.lt.1.0).AND.&
 &      (crfreq.lt.1.0).AND.(omegafreq.lt.1.0).AND.(puckerfreq.lt.1.0)) then
      write(ilog,*) 'WARNING: No polypeptide backbone angles to sample! Adjusting f&
 &requency to 0.0.'
      write(ilog,*)
    end if
    if (otherfreq.gt.0.0) otherfreq = 1.0
    have_pivot = .false.
  end if
  if ((otherfreq.ge.1.0).OR.(puckerfreq.ge.1.0).OR.(rigidfreq.ge.1.0).OR.(particleflucfreq.ge.1.0).OR.(chifreq.ge.1.0).OR.&
 &    (crfreq.ge.1.0).OR.(omegafreq.ge.1.0).OR.(nucfreq.ge.1.0)) have_pivot = .false.
!
! a sanity check against an empty move set
  if ((have_particlefluc.EQV..false.).AND.(have_rigid.EQV..false.).AND.(have_chi.EQV..false.).AND.(have_other.EQV..false.).AND.&
 &(have_cr.EQV..false.).AND.(have_omega.EQV..false.).AND.(have_pucker.EQV..false.).AND.(have_pivot.EQV..false.).AND.&
 &(have_nuc.EQV..false.)) then
    if (pdb_analyze.EQV..false.) then
      write(ilog,*) 'Fatal. Chosen move set is empty (this is before considering possible extra constraints). &
 &Adjust selections and/or system accordingly.'
      call fexit()
    end if
  end if
!
! sanity checks against corrupted move sets
  if ((rigidfreq.lt.1.0).AND.(particleflucfreq.lt.1.0).AND.(crfreq.lt.1.0).AND.(chifreq.lt.1.0).AND.(omegafreq.lt.1.0).AND.&
 & (nucfreq.lt.1.0).AND.(puckerfreq.lt.1.0).AND.(otherfreq.lt.1.0).AND.(fylst%nr.eq.0)) then
    if (pdb_analyze.EQV..false.) then
      write(ilog,*) 'Fatal. Chosen move set is corrupt. Please disable moves that are unable to sample any degrees of &
 &freedom (in particular, FMCSC_OTHERFREQ and associated keywords).'
      call fexit()
    end if
  end if
!
#ifdef ENABLE_MPI
  if ((crfreq.gt.0.0).AND.(angcrfreq.gt.0.0)) then
    if (use_REMC.EQV..true.) then
      if (force_rexyz.EQV..false.) then
        write(ilog,*) 'Warning. When using Ulmschneider-Jorgensen concerted rotation moves &
 &with bond angle variations in conjunction with replica-exchange, it is necessary to force &
 &exchange of Cartesian coordinates (set FMCSC_REMC_DOXYZ to 1).'
        force_rexyz = .true.
      end if
    end if
  end if
#endif
!
! WL sanity checks
  if (mc_acc_crit.eq.3) then
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      write(ilog,*) 'Warning. The use of the Wang-Landau method in conjunction with&
 &fluctuating particle numbers will yield results that are not straightforward &
 &to analyze (if they are interpretable at all).'
    end if
    if ((wld%wl_mode.ne.2).AND.(wld%exrule.eq.3)) then
      write(ilog,*) 'Warning. Full dynamical extension of histograms in Wang-Landau &
 &calculations operating on energy histograms are likely to yield memory &
 &exceptions due to unbound upper range of energies. It is strongly recommended &
 &to use either 1 or 2 for FMCSC_WL_EXTEND.'
    end if
    if ((wld%wl_mode.eq.1).AND.(use_stericscreen.EQV..true.)) then
      if (pdb_analyze.EQV..false.) then
        write(ilog,*) 'Fatal. Using a steric screen for MC moves (FMCSC_USESCREEN) in conjunction &
 &with Wang-Landau sampling of the density of states is incorrect. This error can be overridden &
 &with FMCSC_UNSAFE.'
        if (be_unsafe.EQV..false.) call fexit()
      end if
    else if (use_stericscreen.EQV..true.) then
      write(ilog,*) 'Warning. Using a steric screen for MC moves (FMCSC_USESCREEN) in conjunction &
 &with Wang-Landau sampling is questionable due to the altered acceptance criterion.'
    end if
#ifdef ENABLE_MPI
    if ((use_MPIAVG.EQV..true.).AND.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
      if (pdb_analyze.EQV..false.) then
        write(ilog,*) 'Fatal. Parallel Wang-Landau runs (via FMCSC_MPIAVG) are currently not &
 &supported if the underlying sampler is a hybrid method.'
        call fexit()
      end if
    end if
#endif
  end if
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    if ((re_nbmode.eq.1).AND.(reol_all.EQV..false.).AND.(re_olcalc.le.nsim)&
  &.AND.(re_freq.le.nsim)) then
    write(ilog,*) 'Warning. Non-neighbor exchange is requested, but analysis only&
                  & collects/reports statistics for neighbor overlap.'
    end if
  end if
#endif
!
end
!
!-----------------------------------------------------------------------
!
subroutine hamiltonian_sanitychecks()
!
  use energies
  use iounit
  use sequen
  use cutoffs
!
  implicit none
!
  integer rs
!
! Sanity check for potential energy terms set up
!
  if ((use_WCA.EQV..true.).AND.(use_IPP.EQV..true.)) then
     write(ilog,*)'Warning: WCA and IPP CAN NOT be used together. WC&
 &A is turned off.'
     write(ilog,*)
     use_WCA = .false.
     scale_WCA = 0.0
#ifdef ENABLE_MPI
!    let's not take any risks of confusion here'
     write(ilog,*) 'This is FATAL in MPI-RE. Avoid confusion!'
     call fexit()
#endif
  end if
  if ((use_WCA.EQV..true.).AND.(use_attLJ.EQV..true.)) then
     write(ilog,*)'Warning: WCA and attractive LJ CANNOT be used tog&
 &ether. ATTLJ is turned off.'
     use_attLJ = .false.
     scale_attLJ = 0.0
     write(ilog,*)
#ifdef ENABLE_MPI
!    let's not take any risks of confusion here'
     write(ilog,*) 'This is FATAL in MPI-RE. Avoid confusion!'
     call fexit()
#endif
  end if
!
  if (((use_POLAR.EQV..true.).OR.(use_attLJ.EQV..true.))&
 &     .AND.((use_IPP.EQV..false.).AND.&
 &    (use_WCA.EQV..false.))) then
    write(ilog,*) 'Warning. Using net attractive potentials in the a&
 &bsence of excluded volume terms may lead to numerically unstable s&
 &imulations.'
  end if
!
  if ((use_FEG.EQV..true.).AND.(fegmode.eq.2)) then
    do rs=1,nseq
      if (par_FEG3(rs).EQV..true.) then
        if ((scale_FEGS(15).ne.1.0).OR.(scale_FEGS(16).ne.1.0).OR.&
 &          (scale_FEGS(17).ne.1.0).OR.(scale_FEGS(18).ne.1.0)) then
          write(ilog,*) 'Fatal. Scaling of ghost-ghost interactions &
 &while fully de-coupling parts of the system is not compatible with&
 & values other than unity for the FEG-bonded terms.'
          call fexit()
        end if
      end if
    end do
  end if
!
  if (use_TABUL.EQV..true.) then
    if (use_POLAR.EQV..true.) then
      write(ilog,*) 'Warning. Simultaneous usage of tabulated potentials and&
  & the polar potential is computationally inefficient.'
      if (lrel_md.ne.1) then
        write(ilog,*) 'Warning. All long-range electrostatics treatments outside of truncation&
  & render long-range treatment of tabulated potentials inconsistent.'
      end if
    else
      if (lrel_md.ne.1) then
        write(ilog,*) 'Warning. Long-range electrostatics treatment has no effect on tabulated&
  & potentials.'
      end if
    end if
  end if
!
  if ((par_IMPSOLV(11).ne.0.0).OR.(par_IMPSOLV(10).ne.0.0)) then
    write(ilog,*) 'Fatal. Gross solvent compressibility terms are &
 &currently not supported and are turned off.'
    par_IMPSOLV(10) = 0.0
    par_IMPSOLV(11) = 0.0
  end if

!
  if ((use_cutoffs.EQV..false.).AND.(do_n2loop.EQV..false.)) then
    write(ilog,*) 'Warning. Full pairwise sum must be calculated in runs without cutoffs (bad &
 &setting for FMCSC_N2LOOP).'
    do_n2loop = .true.
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine dynamics_sanitychecks()
!
  use iounit
  use system
  use energies
  use molecule
  use forces
  use atoms
  use polypep
  use zmatrix
  use movesets
  use sequen
  use aminos
  use fyoc
  use torsn
  use mini
  use grandensembles
  use cutoffs
  use mcsums
  use pdb
  use shakeetal
  use params
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer imol,i,k,rs,ttc
  integer, ALLOCATABLE:: rslst(:)
  RTYPE maxrrd,minsdl
  logical afalse,samxyz(3)
!
  afalse = .false.
!
! ensemble-related
  if (ens%flag.eq.1) then
!   currently there is no actual restrictions on the use of NVT
  else if (ens%flag.eq.2) then
    if (dyn_mode.ne.2) then
      write(ilog,*) 'Fatal. The NVE ensemble is only accessible in N&
 &ewtonian (full inertial) dynamics. All other approaches implemente&
 &d in CAMPARI automatically couple the system to a heat bath.'
     call fexit()
    end if
#ifdef ENABLE_MPI
    if ((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.)) then
      write(ilog,*) 'Warning. A PIGS calculation in the NVE ensemble is not generally meaningful.'
      if (re_velmode.ne.1) then
        write(ilog,*) 'Fatal. To enable this calculation anyway, velocity randomization upon reseeding must be &
 &turned on (FMCSC_RE_VELMODE).'
        call fexit()
      end if
    end if
#endif
    if (((add_settle.EQV..true.).OR.(cart_cons_mode.ne.1)).AND.(fycxyz.eq.2)) then
      write(ilog,*) 'Warning. First step of NVE simulation in the presence &
 &of holonomic constraints in Cartesian space will be erroneous and lead to temperature drop.'
    end if
  else if (ens%flag.eq.3) then
!   remove this top level sanity check once ready
    write(ilog,*) 'The NPT (isobaric-isothermal) ensemble is currently not supported. Please check back later.'
    call fexit()
    if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
      write(ilog,*) 'Manostats are not yet supported for MC calculat&
 &ions. Please check back later.'
      call fexit()
    end if
    if ((dyn_mode.eq.3).OR.(dyn_mode.eq.4).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8).OR.&
 &                             (use_IMPSOLV.EQV..true.)) then
      write(ilog,*) 'Manostats are not yet supported for calculation&
 &s involving an underlying continuum such as the ABSINTH implicit s&
 &olvation model or LD/BD approaches. Please check back later.'
      call fexit()
    else if ((dyn_mode.eq.2).OR.(dyn_mode.eq.5)) then
      if ((bnd_type.le.2).OR.(bnd_type.gt.4).OR.(bnd_shape.ne.1)) then
        write(ilog,*) 'Manostats are currently only supported for so&
 &ft-wall droplet boundary conditions and spherical simulation containers. Please check back later.'
        call fexit()
      end if
    end if
  else if (ens%flag.eq.4) then
!   remove this top level sanity check once ready
    write(ilog,*) 'The NPE (isobaric-isenthalpic) ensemble is currently not supported. Please check back later.'
    call fexit()
    if (dyn_mode.ne.2) then
      write(ilog,*) 'Fatal. The NPE (isobaric-isenthalpic) ensemble &
 &is not defined for approaches outside of Newtonian (full inertial)&
 & dynamics.'
      call fexit()
    else
      if ((bnd_type.le.2).OR.(bnd_type.gt.4).OR.(bnd_shape.ne.1)) then
        write(ilog,*) 'Manostats are currently only supported for so&
 &ft-wall droplet boundary conditions and spherical containers. Please check back later.'
        call fexit()
      end if
      write(ilog,*) 'Warning. The NPE (isobaric-isenthalpic) ensembl&
 &e currently outputs only internal energy, not enthalpy.'
    end if
  else if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
#ifdef ENABLE_MPI
    write(ilog,*) 'Fatal. Simulations in the (semi-)grand ensemble a&
 &re currently not supported for MPI runs. Check back later.'
    call fexit()
#endif
    if (dyn_mode.ne.1) then
      write(ilog,*) 'Fatal. The (semi-)grand ensemble is only access&
 &ible to standard (pure) Monte Carlo simulations at this point.'
      call fexit()
    end if
    if (ntorpuck.gt.0) then
      write(ilog,*) 'Warning. The existence of internal (torsional) &
 &degrees of freedom might lead to artifacts when using the (semi-)g&
 &rand ensembles.'
    end if
    if (fluclst%nr.le.0) then
      write(ilog,*) 'Fatal. The grand or semi-grand ensembles are us&
 &ed with no particle types allowed to fluctuate. Use NVT instead.'
      call fexit()
    end if
  end if
!
! thermostat/manostat related
  if ((dyn_mode.eq.2).OR.(dyn_mode.eq.5).OR.((mini_mode.eq.4).AND.(dyn_mode.eq.6))) then
    if (ens%flag.eq.1) then
      if (tstat%flag.eq.3) then
        write(ilog,*) 'Fatal. Extended ensemble thermostat (Nose-Hoover-like) not yet supported for NVT ensemble.'
        call fexit()
      end if
    else if (ens%flag.eq.3) then
      if (((tstat%flag.eq.3).AND.(pstat%flag.ne.3)).OR.((tstat%flag.ne.3).AND.(pstat%flag.eq.3))) then
        write(ilog,*) 'Fatal. Extended ensemble thermo- and manostat (Nose-Hoover-like) have to be used concurrently in the NPT&
 &ensemble.'
        call fexit()
      else if (pstat%flag.eq.2) then
        write(ilog,*) 'Fatal. Markov chain manostat is not yet supported. Please check back later.'
        call fexit()
      else if ((pstat%flag.eq.1).AND.((bnd_shape.ge.2).AND.(bnd_shape.le.4))) then
        write(ilog,*) 'Warning. Droplet manostat might give rise to artifacts in conjunction with standard thermostat (like &
 &Andersen or Berendsen) i.e., might sample an incorrect ensemble.'
      else if ((tstat%flag.eq.3).AND.(pstat%flag.eq.3)) then
        write(ilog,*) 'Fatal. Extended ensemble thermo- and manostat (Nose-Hoover-like) is not yet fully supported. Please &
 &check back later.'
        call fexit()
      end if
    else if (ens%flag.eq.4) then
      if (.NOT.((pstat%flag.eq.1).AND.((bnd_shape.ge.2).AND.(bnd_shape.le.4)))) then
        write(ilog,*) 'Fatal. Support for the NPE (isobaric-isenthalpic) ensemble is currently restricted to the droplet manostat.&
 & Please check back later.'
        call fexit()
      else if (pstat%flag.eq.2) then
        write(ilog,*) 'Fatal. Markov chain manostat is not yet supported. Please check back later.'
        call fexit()
      end if
    end if
    if (((tstat%flag.eq.1).OR.(tstat%flag.eq.4)).AND.((ens%flag.eq.1).OR.(ens%flag.eq.3))) then
      write(ilog,*) 'Warning. Irrespective of boundary conditions, global velocity rescaling thermostats can &
 &lead to severe equipartition / drift artifacts in dynamics. These effects are more common for heterogeneous &
 &systems with a mix of weak and strong energetic coupling between degrees of freedom, and are seen more easily in Cartesian &
 &rather than internal coordinate space dynamics.'
    end if
#ifdef ENABLE_MPI
    if ((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.)) then
      if (((tstat%flag.eq.1).OR.(tstat%flag.eq.3)).AND.((ens%flag.eq.1).OR.(ens%flag.eq.3)).AND.(re_velmode.ne.1)) then
        write(ilog,*) 'Fatal. For PIGS, the propagator should have a stochastic component. If not, velocity &
 &randomization after reseeding must be turned on (FMCSC_RE_VELMODE).'
        call fexit()
      end if
    end if
#endif
    if ((tstat%flag.eq.1).AND.((ens%flag.eq.1).OR.(ens%flag.eq.3))) then
      write(ilog,*) 'Warning. The Berendsen weak-coupling scheme does not provide sampling of the proper&
  & ensemble (NVT or NPT) by suppressing fluctuations (leading to incorrect heat capacities for instance).'
    end if
    if ((fycxyz.eq.2).AND.(ens%flag.ne.2).AND.(ens%flag.ne.4).AND.&
 &      (tstat%flag.eq.2).AND.((cart_cons_mode.gt.1).OR.(add_settle.EQV..true.))) then
      write(ilog,*) 'Fatal. The Andersen stochastic heat bath is incompatible with holonomic constraints &
  &in Cartesian dynamics since it decorrelates velocities within constraint groups. This leads to continuous &
  &energy dissipation upon velocity resets and subsequent constraint solving, and the resultant ensemble is ill-defined.'
      write(ilog,*) 'Use the Bussi-Parrinello thermostat instead, and/or switch to torsional dynamics.'
      if (be_unsafe.EQV..false.) call fexit()
    end if
  else if ((dyn_mode.eq.3).OR.(dyn_mode.eq.7)) then
    if ((fycxyz.eq.2).AND.((cart_cons_mode.gt.1).OR.(add_settle.EQV..true.))) then
      write(ilog,*) 'Warning. The Skeel-Izaguirre impulse integrator for Langevin dynamics has not been shown to &
  &sample the proper ensemble if holonomic constraints are present in Cartesian dynamics. Consult Chem. Phys. Lett. &
  &429 (2006), pp. 310-316 for details on this issue. Consider using Newtonian MD with the Bussi-Parrinello &
  &thermostat instead, and/or switching to torsional dynamics.'
    else if (fycxyz.eq.1) then
      write(ilog,*) 'Warning. The assumption of a uniform (flat) friction coefficient is highly questionable &
 &in Langevin (stochastic) dynamics when performed in complex, internal coordinates (flexible polymers). This implies that &
 &kinetics may be slowed critically for complex conformational transitions. Consider &
 &using torsional dynamics with the Bussi-Parrinello thermostat instead.'
    end if
  end if
#ifdef ENABLE_MPI
  if ((use_REMC.EQV..true.).AND.(re_velmode.eq.1).AND.(fycxyz.eq.2).AND.(re_freq.le.nsim).AND.&
 &    ((cart_cons_mode.gt.1).OR.(add_settle.EQV..true.))) then
    write(ilog,*) 'Warning. In REMC, randomization of velocities upon successful swaps is not recommended if holonomic &
 &constraints are in use in Cartesian space. Changing to re-scaling (option #2).'
    re_velmode = 2
  end if
#endif
!
! constraint-related
  if ((n_constraints.gt.0).AND.(do_frz.EQV..true.).AND.(use_IMPSOLV.EQV..true.).AND.&
 &(skip_frz.EQV..true.).AND.(fycxyz.eq.1)) then
    write(ilog,*) 'Warning. In torsional/rigid-body dynamics the skipping of constrained&
 & pairwise forces is illegal if forces are not pairwise decomposable. Disabling&
 & request for FMCSC_SKIPFRZ.'
    skip_frz = .false.
  end if
  if ((do_frz.EQV..false.).AND.(skip_frz.EQV..true.)) then
    skip_frz = .false.
  end if
  if (fycxyz.eq.2) then
    if (cart_cons_source.eq.3) then
      write(ilog,*) 'Warning. The use lof low precision input structures to define values for constraints &
 &is not recommended. Constraints will be slightly heterogeneous for identical chemical units, and differ &
 &if a different input structure of the same system is used.'
    else if ((use_dyn.EQV..true.).AND.(pdb_analyze.EQV..false.).AND.(do_restart.EQV..false.).AND.(cart_cons_grps.gt.0)) then
      write(ilog,*) 'Warning. Due to a structure-independent choice of constraint lengths (FMCSC_SHAKEFROM), the initial &
 &geometry is adjusted before the start of the simulation (see ',basename(1:bleng),'_START.*).'
    end if
    if ((cart_cons_mode.eq.4).AND.(cart_cons_method.ne.1)) then
      write(ilog,*) 'Warning. Explicit bond angle constraints are currently only supported with the simple SHAKE &
 &procedure. Adjusting algorithm choice accordingly.'
      cart_cons_method = 1
    else if (((cart_cons_mode.eq.4).OR.(cart_cons_mode.eq.5)).AND.(cart_cons_method.eq.4)) then
      write(ilog,*) 'Warning. Any type of highly coupled (bond angle) constraints will cause the LINCS procedure to crash.&
 & Adjusting algorithm choice accordingly (standard SHAKE).'
      cart_cons_method = 1
    end if
  end if
! c.o.m.-drift removal related
  if (ens%sysfrz.ge.2) then
    if ((.NOT.((bnd_type.eq.1).AND.(bnd_shape.eq.1))).AND.(dyn_mode.ne.1)) then
      write(ilog,*) 'Warning. Removing center-of-mass translation an&
 &d rotation can lead to misleading results for (partially) non-periodic boundar&
 &y conditions.'
    end if
    if ((fycxyz.eq.1).AND.(nmol.eq.1)) then
!     with just one molecule in the system we can just suppress rigid body
!     motion for that molecule
      do_frz = .true.
      if (((dyn_mode.ge.2).AND.(dyn_mode.le.5)).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
        do imol=1,nmol
          if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
            ttc = 3
          else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
            ttc = 5
          else
            ttc = 6
          end if
          do i=1,ttc
            if (dc_di(imol)%frz(i).EQV..false.) then
              dc_di(imol)%frz(i) = .true.
              n_constraints = n_constraints + 1
            end if
          end do
        end do
      end if
    else if (fycxyz.eq.1) then
!     with multiple molecules, we have to do system-wide corrections
      if (ens%sysfrz.eq.3) then
!       if there's exactly two molecules one rotational dof is missing
        if (nmol.eq.2) then
          n_constraints = n_constraints + 5
        else
          n_constraints = n_constraints + 6
        end if
      else
        n_constraints = n_constraints + 3
      end if
    else if (fycxyz.eq.2) then
!     in Cartesian MD all atom-based
      samxyz(:) = .false.
      do i=1,n
        if (cart_frz(i,1).EQV..false.) samxyz(1) = .true.
        if (cart_frz(i,2).EQV..false.) samxyz(2) = .true.
        if (cart_frz(i,3).EQV..false.) samxyz(3) = .true.
      end do
      k = 0
      do i=1,3
        if (samxyz(i).EQV..false.) k = k + 1 ! k can never be 3
      end do
      if (ens%sysfrz.eq.3) then
!       if there's exactly two molecules one rotational dof is missing
        if (n.eq.2) then
          n_constraints = n_constraints + 5 - k - min(k,1)*k
        else if (n.eq.1) then
          n_constraints = n_constraints + 3 - k
        else
          n_constraints = n_constraints + 6 - k - min(k,1)*(k+1)
        end if
      else
        n_constraints = n_constraints + 3 - k
      end if
    end if
  else if (ens%sysfrz.eq.1) then
    if ((bnd_type.eq.1).AND.(dyn_mode.ne.1).AND.((ens%flag.eq.2).OR.(ens%flag.eq.4).OR.&
 &   ((tstat%flag.ne.3).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.5))))) then
      write(ilog,*) 'Warning. The lack of removal of center-of-mass &
 &translation and rotation can sometimes lead to drift artifacts in (partially) periodic b&
 &oundary conditions.'
    end if
  end if
!
! unsupported residue-type related
  if (fycxyz.ne.2) then
    do rs=1,nseq
      if ((dyn_mode.ne.1).AND.((pucline(rs).gt.0).OR.(nucsline(6,rs).gt.0))) then
        write(ilog,*) 'Warning. Ring puckering is currently not supp&
 &orted for torsional space dynamics sampling (5-ring of residue ',rs,' is frozen).&
 & It may be possible to use hybrid MC/TMD sampling to address this issue.'
      end if
    end do
  end if
  if (n_pdbunk.gt.0) then
    allocate(rslst(nseq))
    rslst(:) = 0
    do i=1,nseq
      if (i.eq.rsmol(molofrs(i),1)) cycle
      rslst(i) = atmres(iz(1,at(i)%bb(1)))
    end do
    do i=1,nseq
      if ((rslst(i).gt.0).AND.(rslst(i).ne.(i-1))) then
        if (pdb_analyze.EQV..false.) then
          write(ilog,*) 'Fatal. Branched polymers (i.e., a side chain (branch) using separate residues) are only &
 &marginally supported in dynamics. The short-range interaction model will be erroneous (wrong exclusion rules, missing &
 &or incorrect bonded parameters). You can use FMCSC_UNSAFE to circumvent this error, but the resultant simulation &
 &(including analysis) is unlikely to be meaningful (especially when sampling in Cartesian space).'
        else
          write(ilog,*) 'Fatal. Branched polymers (i.e., a side chain (branch) using separate residues) are only &
 &marginally supported in trajectory analysis mode. Several analyses may fail and recomputed energies are likely &
 &to be meaningless. You can use FMCSC_UNSAFE to circumvent this error.'
        end if
        if (be_unsafe.EQV..false.) call fexit()
        exit
      end if
    end do
    deallocate(rslst)
  end if
!
! unsupported energy term
  if (dyn_mode.ne.1) then
    if ((par_IMPSOLV(11).ne.0.0).OR.(par_IMPSOLV(10).ne.0.0)) then
      write(ilog,*) 'Fatal. Gross solvent compressibility terms are &
 &currently not supported in dynamics runs. Check back later.'
      call fexit()
    end if
  end if
!
! Cartesian dynamics-related
  if (fycxyz.eq.2) then
    if (use_CORR.EQV..true.) then
      write(ilog,*) 'Warning. The use of structural correction terms&
 & in Cartesian dynamics is redundant with a standard parameter file&
 & accounting for bonded interactions of this type (if FMCSC_SC_BONDED_T is unity).'
    end if
    if (ens%flag.ge.3) then
      write(ilog,*) 'Fatal. Only the NVT and NVE ensembles are suppo&
 &rted for Cartesian dynamics at this time. Check back later.'
      call fexit()
    end if
    if ((dyn_mode.ne.6).OR.(mini_mode.eq.4)) then
      do i=1,n
        if (mass(i).le.0.0) then
          write(ilog,*) 'Warning. Massless (virtual) particles (sites) are fatal in unconstrained dynamics (first&
 & encountered for atom ',i,' of biotype ',b_type(i),'). Make sure that those sites &
 &are part of rigid constraint (SHAKE/SETTLE) groups.'
          exit
        end if
      end do
    end if
  end if
!
! IMD-related
  if ((fycxyz.eq.1).AND.(use_dyn.EQV..true.)) then
    do i=1,nmol
      call update_rigidm(i)
      dc_di(i)%f(:) = 0.0
      dc_di(i)%im(:) = 0.0
    end do
    call cart2int_I()
    do i=1,nmol
      dc_di(i)%im(:) = dc_di(i)%olddat(:,2)
      do k=1,size(dc_di(i)%im)
        if ((dc_di(i)%im(k).le.0.05).AND.(k.gt.6).AND.(dc_di(i)%frz(k).EQV..false.)) then
          write(ilog,*) 'Warning. Constraining an additional dihedral angle degree of freedom (#',k,') in molecule #',i,',&
 &because it has miniscule inertia (this can happen for colinear atoms, or if virtual (massless) &
 &sites are in use).'
          dc_di(i)%frz(k) = .true.
          n_constraints = n_constraints + 1
        end if
      end do
      dc_di(i)%f(:) = 0.0
    end do
  end if
!
! integrator-related
  if ((fycxyz.eq.1).AND.((dyn_mode.eq.3).OR.(dyn_mode.eq.7))) then
    write(ilog,*) "Warning. The TLD integrator based on Skeel's impulse integrator for solving uncoupled stochastic &
 &dynamics equations works only with the approximate, fully deterministic TMD scheme (in development). &
 &Better stability is probably achieved with variants of the TMD (Newtonian) integrator. Treat results with caution."
    dyn_integrator_ops(1) = 0
    dyn_integrator_ops(2) = 2
  end if
!
! minimizer related
  if (dyn_mode.eq.6) then
    if (ens%flag.ne.1) then
      write(ilog,*) 'Warning. Minimizers technically operate in cons&
 &tant volume conditions. Adjusting ensemble specification.'
      ens%flag = 1
    end if
#ifdef ENABLE_MPI
    write(ilog,*) 'Fatal. Minimizers are currently not supported in &
 &MPI runs. Check back later.'
    call fexit()
#endif
    if ((cart_cons_mode.ne.1).OR.(add_settle.EQV..true.)) then
      if (mini_mode.le.3) then
        write(ilog,*) 'Warning. Canonical minimizers are not compatible with holonomic constraints at the moment. &
 &It is recommended to use the stochastic option (FMCSC_MINI_MODE 4), switch to internal coordinates (FMCSC_CARTINT 1), &
 &or to simply perform dynamics at very low temperatures and small time steps.'
        add_settle = .false.
        cart_cons_mode = 1
      end if
    end if
  end if
!
! cutoff/PBC-related
  if (bnd_type.eq.1) then
    minsdl = huge(minsdl)
    if (bnd_shape.eq.1) then
      do i=1,3
        if (bnd_params(i).lt.minsdl) then
          minsdl = bnd_params(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      minsdl = bnd_params(6)
    end if
    maxrrd = tiny(maxrrd)
    do i=1,nseq
      if (resrad(i).gt.maxrrd) maxrrd = resrad(i)
    end do
    if ((use_IPP.EQV..true.).OR.(use_WCA.EQV..true.).OR.&
 &      (use_attLJ.EQV..true.).OR.(use_IMPSOLV.EQV..true.).OR.&
 &      (use_FEGS(1).EQV..true.).OR.(use_FEGS(3).EQV..true.)) then
      if ((mcnb_cutoff+2.0*maxrrd).gt.(0.5*minsdl)) then
        write(ilog,*) 'Warning. Short-range cutoff will lead to non-&
 &spherical symmetry due to system size. Consider increasing the sys&
 &tem size or decreasing cutoffs.'
      end if
    end if
    if ((use_TABUL.EQV..true.).OR.(use_POLAR.EQV..true.).OR.&
 &      (use_FEGS(6).EQV..true.)) then
      if ((mcel_cutoff+2.0*maxrrd).gt.(0.5*minsdl)) then
        write(ilog,*) 'Warning. Mid-(long)-range cutoff will lead to&
 & non-spherical symmetry due to system size. Consider increasing th&
 &e system size or decreasing cutoffs.'
      end if
    end if
    if ((dyn_mode.ne.1).AND.(lrel_md.ge.4).AND.(cglst%ncs.gt.0)) then
      write(ilog,*) 'Warning. In PBC, any non-spherical cutoff schemes can lead&
 & to significant artifacts. In the proposed calculation, the long-range&
 & electrostatics treatment (LREL_MD) creates such a non-spherical symmetry.'
    end if
    if (((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)).AND.&
 &(lrel_mc.ne.4).AND.(cglst%ncs.gt.0)) then
      write(ilog,*) 'Warning. In PBC, any non-spherical cutoff schemes can lead&
 & to significant artifacts. In the proposed calculation, the long-range&
 & electrostatics treatment (LREL_MC) creates such a non-spherical symmetry.'
    end if
    if (((lrel_mc.eq.4).OR.(lrel_md.eq.1)).AND.((cglst%ncs.gt.0).AND.((use_POLAR.EQV..true.).OR.&
 &      (use_FEGS(6).EQV..true.)))) then
      write(ilog,*) 'Warning. The simple truncation of long-range electrostatic interactions can &
 &lead to severe partitioning artifacts for almost all simulation setups.'
    end if
  else if ((bnd_type.gt.1).AND.(bnd_shape.ne.2)) then
    write(ilog,*) 'Warning. Should nonbonded interactions be in use, it must be kept in mind that &
 &the use of a nonspherical, aperiodic simulation container is likely to lead to major artifacts.'
  end if
!
! hybrid method-related
  if ((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    if (ens%flag.ne.1) then
      write(ilog,*) 'Fatal. Hybrid Monte Carlo / dynamics methods are only supported in&
 & the canonical (NVT) ensemble at the moment. Check back later for NPT support.'
      call fexit()
    end if
    if ((cart_cons_mode.gt.1).AND.(fycxyz.eq.2)) then
      write(ilog,*) 'Warning. In hybrid MC/MD sampling, velocities are normally randomly reassigned &
 &at the beginning of a dynamics segment. This can lead to temperature artifacts if the dynamics employ &
 &holonomic constraints (due to velocity decorrelation). Adjusting algorithm to remember velocities from &
 &previous segment (this will help with, but not eliminiate such artifacts). Consider using torsional &
 &dynamics instead.'
    end if    
    if ((dyn_mode.eq.5).AND.(tstat%flag.eq.1)) then
      write(ilog,*) 'Warning. It is incorrect to use hybrid MC/MD employing the Berendsen weak-coupling scheme&
 & to maintain constant temperature in the MD portion. This is because the latter does not rigorously&
 & sample the NVT ensemble, which the MC portion does.'
    end if
    if (use_POLAR.EQV..true.) then
      if ((lrel_md.eq.2).OR.(lrel_md.eq.3)) then
        write(ilog,*) 'Warning. Advanced long-range electrostatics treatments like (G)RF or Ewald&
 & are not supported in MC calculations, hence their use in hybrid methods will lead to&
 & inconsistencies and potential artifacts (see FMCSC_LREL_MD).'
      end if
      if (((lrel_md.eq.4).AND.(lrel_mc.ne.3)).OR.((lrel_md.eq.5).AND.(lrel_mc.ne.1))&
 &.OR.((lrel_md.eq.1).AND.(lrel_mc.ne.4))) then
        write(ilog,*) 'Warning. Simple approximations for long-range electrostatic&
 & interactions should be employed consistently for both parts of the hybrid MC/XD calculation.&
 & The current settings will give rise to potential artifacts (FMCSC_LREL_MD and FMCSC_LREL_MC).'
      end if
    end if
    if (mcel_cutoff.ne.mcnb_cutoff) then
      write(ilog,*) 'Warning. In hybrid runs, the most consistent choice of cutoffs is to&
 & eliminate the twin-range regime in the dynamics portion, i.e., to set both cutoffs to the&
 & same value. Otherwise, short-range interactions are truncated differently in the two&
 & sub-methods.'
    end if
    if (fycxyz.eq.2) then
      write(ilog,*) 'Warning. When using Cartesian dynamics within a hybrid approach, the MC&
 & moveset becomes severely non-ergodic. This implies that convergence to the correct ensemble&
 & averages will be significantly slowed due to the requirement of sampling a large number of&
 & MC segments. This is necessary to average out the temporary biases introduced by the&
 & non-ergodicity.'
    end if
    if (min_dyncyclen.gt.max_dyncyclen) then
      write(ilog,*) 'Warning. Bad specifications for FMCSC_CYCLE_DYN_MIN and FMCSC_CYCLE_DYN_MAX.&
 & Restoring to defaults.'
      min_dyncyclen = 500
      max_dyncyclen = 1000
    end if
    if (min_mccyclen.gt.max_mccyclen) then
      write(ilog,*) 'Warning. Bad specifications for FMCSC_CYCLE_MC_MIN and FMCSC_CYCLE_MC_MAX.&
 & Restoring to defaults.'
      min_mccyclen = 500
      max_mccyclen = 1000
    end if
!    if (first_mccyclen.
  end if
!
!  miscellanei
  if (use_dyn.EQV..true.) then
    if ((have_ph.EQV..true.).OR.(phout.lt.nsim)) then
      if ((have_ph.EQV..true.).AND.((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
        write(ilog,*) 'Warning. In hybrid MC/dynamics runs, pH-pseudo moves are not available. Turning off.'
      end if
      if (phout.lt.nsim) then
        write(ilog,*) 'Warning. Titration state output relies on pH-pseudo moves and is not available in&
 & any form of calculation employing dynamics. Turning off.'
      end if
      phout = nsim + 1
      phfreq = 0.0
    end if
  end if
!
end
!
#ifdef ENABLE_MPI
!
!------------------------------------------------------------------
!
! input file-dependent Hamiltonians pose a challenge, since bad input
! will yield a bad (or meaningless) RE-calculation
!
subroutine remc_sanitychecks()
!
  use iounit
  use mpistuff
  use energies
  use system
  use movesets
!
  implicit none
!
  integer k
!
  do k=1,re_conddim
    if ((re_types(k).eq.5).AND.(use_POLAR.EQV..false.)) then
      write(ilog,*) 'Fatal. Polar potential disabled after processin&
 &g. REMC calculation potentially meaningless.'
      call fexit()
    end if
    if ((re_types(k).eq.6).AND.(use_IMPSOLV.EQV..false.)) then
      write(ilog,*) 'Fatal. Implicit solvent potential disabled afte&
 &r processing. REMC calculation potentially meaningless.'
      call fexit()
    end if
    if ((re_types(k).eq.8).AND.(use_TOR.EQV..false.)) then
      write(ilog,*) 'Fatal. Torsional potential disabled after file &
 &processing. This might render this REMC calculation meaningless.'
      call fexit()
    end if
    if ((re_types(k).eq.20).AND.(use_FEGS(1).EQV..false.)) then
      write(ilog,*) 'Fatal. FEG-IPP potential disabled after file pr&
 &ocessing. This might render this REMC calculation meaningless.'
      call fexit()
    end if
    if ((re_types(k).eq.21).AND.(use_FEGS(3).EQV..false.)) then
      write(ilog,*) 'Fatal. FEG-attLJ potential disabled after file &
 &processing. This might render this REMC calculation meaningless.'
      call fexit()
    end if
    if ((re_types(k).eq.22).AND.(use_FEGS(6).EQV..false.)) then
      write(ilog,*) 'Fatal. FEG-POLAR potential disabled after file &
 &processing. This might render this REMC calculation meaningless.'
      call fexit()
    end if
    if ((re_types(k).eq.23).AND.(use_TABUL.EQV..false.)) then
      write(ilog,*) 'Fatal. Tabulated potential disabled after file &
 &processing. This might render this REMC calculation meaningless.'
      call fexit()
    end if
    if ((re_types(k).eq.24).AND.(use_POLY.EQV..false.)) then
      write(ilog,*) 'Fatal. Polymeric potential disabled after file &
 &processing. This might render this REMC calculation meaningless.'
      call fexit()
    end if
    if ((re_types(k).eq.25).AND.(use_DREST.EQV..false.)) then
      write(ilog,*) 'Fatal. Dis. restr. potential disabled after fil&
 &e processing. This might render REMC calculation meaningless.'
      call fexit()
    end if
!
  end do
!
  if ((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    write(ilog,*) 'Warning. Force routines are used exclusively for overlap&
 & calculations in this hybrid calculation whereas actual swaps use either&
 & force or energy routines depending on the sampling step they occur in.' 
  end if
!
  if (ens%flag.ne.1) then
    write(ilog,*) 'Fatal. Replica-exchange calculations are only supported in the&
 & canonical (NVT) ensemble at this point in time.'
    call fexit()
  end if
!
  if ((mc_acc_crit.eq.3).AND.(re_freq.le.nsim)) then
    write(ilog,*) 'Warning. Replica-exchange swap moves introduce canonical &
 &sampling moves to Wang-Landau simulations. This may hinder convergence and &
 &distort results.'
  end if
! 
end
!
!---------------------------------------------------------------------------
!
#endif
!
subroutine loop_checks()
!
  use energies
  use cutoffs
  use iounit
  use system
!
  implicit none
!
! if the Hamiltonian is "zero" we can save a lot of time ...
! note that ideal_run enforces only quitting the res.-res. or res.-nbl
! energy fxns prematurely (i.e., ZSEC, TOR, ... might still be used)
! other simplifications are similar (see energy and force routines and
! wrappers)
  if ((use_IPP.EQV..false.).AND.(use_attLJ.EQV..false.).AND.&
 &    (use_WCA.EQV..false.).AND.(use_CORR.EQV..false.).AND.&
 &    (use_IMPSOLV.EQV..false.).AND.(use_POLAR.EQV..false.).AND.&
 &    (use_TABUL.EQV..false.).AND.(use_FEG.EQV..false.)) then
    ideal_run = .true.
  else
    if ((nindx.ne.12).AND.(dyn_mode.ne.1)) then
      write(ilog,*) 'Warning. This calculation might run much slower than expected&
                   & due to non-standard exponent for LJ potential.'
!     do nothing
    else if ((use_IPP.EQV..true.).AND.(use_attLJ.EQV..true.).AND.&
 &    (use_WCA.EQV..false.).AND.(use_TABUL.EQV..false.).AND.&
 &    (use_IMPSOLV.EQV..false.).AND.(use_POLAR.EQV..false.).AND.&
 &    (use_FEG.EQV..false.)) then
      is_lj = .true.
    else if ((use_IPP.EQV..true.).AND.(use_attLJ.EQV..false.).AND.&
 &    (use_WCA.EQV..false.).AND.(use_TABUL.EQV..false.).AND.&
 &    (use_IMPSOLV.EQV..false.).AND.(use_POLAR.EQV..false.).AND.&
 &    (use_FEG.EQV..false.)) then
      is_ev = .true.
    else if ((use_IPP.EQV..true.).AND.(use_attLJ.EQV..true.).AND.&
 &    (use_WCA.EQV..false.).AND.(use_TABUL.EQV..false.).AND.&
 &    (use_IMPSOLV.EQV..false.).AND.(use_POLAR.EQV..true.).AND.&
 &    (use_FEG.EQV..false.)) then
      if ((lrel_md.eq.1).OR.(lrel_md.eq.4).OR.(lrel_md.eq.5)) then
        is_plj = .true.   ! TR cutoffs w/ wo/ extra terms
      else if (lrel_md.eq.2) then
        is_pewlj = .true. ! Ewald treatment 
      else if (lrel_md.eq.3) then
        is_prflj = .true. ! RF treatment
      end if
    else if ((use_IPP.EQV..true.).AND.(use_attLJ.EQV..true.).AND.&
 &    (use_WCA.EQV..false.).AND.(use_TABUL.EQV..true.).AND.&
 &    (use_IMPSOLV.EQV..false.).AND.(use_POLAR.EQV..false.).AND.&
 &    (use_FEG.EQV..false.)) then
      is_tablj = .true.
    else if ((use_IPP.EQV..false.).AND.(use_attLJ.EQV..false.).AND.&
 &    (use_WCA.EQV..false.).AND.(use_TABUL.EQV..true.).AND.&
 &    (use_IMPSOLV.EQV..false.).AND.(use_POLAR.EQV..false.).AND.&
 &    (use_FEG.EQV..false.)) then
      is_tab = .true.
    else if ((use_IPP.EQV..true.).AND.(use_attLJ.EQV..true.).AND.&
 &    (use_WCA.EQV..false.).AND.(use_TABUL.EQV..false.).AND.&
 &    (use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.).AND.&
 &    (use_FEG.EQV..false.)) then
      is_impljp = .true.
    else if ((use_IPP.EQV..true.).AND.(use_attLJ.EQV..true.).AND.&
 &    (use_WCA.EQV..false.).AND.(use_TABUL.EQV..false.).AND.&
 &    (use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..false.).AND.&
 &    (use_FEG.EQV..false.)) then
      is_implj = .true.
    else if ((use_IPP.EQV..true.).AND.(use_attLJ.EQV..true.).AND.&
 &    (use_WCA.EQV..false.).AND.(use_TABUL.EQV..false.).AND.&
 &    (use_IMPSOLV.EQV..false.).AND.(use_POLAR.EQV..true.).AND.&
 &    (use_FEG.EQV..true.)) then
      if ((lrel_md.eq.1).OR.(lrel_md.eq.4).OR.(lrel_md.eq.5)) then
        is_fegplj = .true.   ! TR cutoffs w/ wo/ extra terms
      else if ((lrel_md.eq.3).AND.(fegcbmode.eq.1)) then
        is_fegprflj = .true. ! RF treatment
      end if
    else if ((use_IPP.EQV..true.).AND.(use_attLJ.EQV..true.).AND.&
 &    (use_WCA.EQV..false.).AND.(use_TABUL.EQV..false.).AND.&
 &    (use_IMPSOLV.EQV..false.).AND.(use_POLAR.EQV..false.).AND.&
 &    (use_FEG.EQV..true.)) then
      is_feglj = .true.
    else if ((use_IPP.EQV..true.).AND.(use_attLJ.EQV..false.).AND.&
 &    (use_WCA.EQV..false.).AND.(use_TABUL.EQV..false.).AND.&
 &    (use_IMPSOLV.EQV..false.).AND.(use_POLAR.EQV..false.).AND.&
 &    (use_FEG.EQV..true.)) then
      is_fegev = .true.
    end if
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine analysis_checks()
!
  use torsn
  use polyavg
  use molecule
  use system
  use iounit
  use paircorr
  use dssps
  use contacts
  use energies
  use sequen
  use fyoc
  use dipolavg
  use diffrac
  use clusters
  use grandensembles
  use pdb
  use mcsums
  use movesets
  use mpistuff
  use fos
  use ems
!
  implicit none
!
  integer t1,t2,k,kk,i,j,modstep
  logical logdummy,turnoff(50)
  RTYPE ppp(4)
!
!! sanity check against FP (this is disabled due to speed and it being relevant for nonrepresentable
! values of nsim only)
!  if (nequil.lt.nsim) then
!    ppp(:) = 0.
!    do i=nequil+1,nsim
!      ppp(1) = ppp(1) + 1.0
!    end do
!    if (nint(ppp(1)).ne.(nsim-nequil)) then
!      write(ilog,*) 'Warning. If analyses are performed very frequently, the length of the &
! &simulation may lead to significant roundoff errors in accumulating averages and distribution functions.'
!    end if
!  end if
!
! sanity check against FP weights
  if (pdb_analyze.EQV..true.) then
    if (use_frame_weights.EQV..true.) then
      if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
        write(ilog,*) 'Fatal. The use of frame weights in trajectory analysis is incompatible with &
 &assumed ensembles with fluctuating particle numbers.'
        call fexit()
      end if
      ppp(4) = HUGE(ppp(4))
      do i=1,framecnt
        if ((framewts(i).gt.0.0).AND.(framewts(i).lt.ppp(4))) ppp(4) = framewts(i)
      end do
      ppp(1) = sum(framewts(:))+ppp(4)
      ppp(2) = sum(framewts(:))
      ppp(3) = ppp(1)-ppp(2)
!     WCS (single addition) with a 1% threshold
!     to be relevant, this assumes quantities multiplied with the weights have negligible variance
      if (abs(ppp(3)-ppp(4)).gt.1.0e-2*ppp(4)) then
        write(ilog,*) 'Warning. In using weights on snapshots in trajectory analysis mode, &
 &round-off errors on computed averages and distribution functions are likely to be significant &
 &due to spread of weight values.'
      end if
    end if
  end if
!
! some generic things
  if (pdb_analyze.EQV..true.) then
    ensout = nsim+1
    accout = nsim+1
    if (do_restart.EQV..true.) then
      write(ilog,*) 'Fatal. It is not possible to restart a trajectory analysis run.'
      call fexit()
    end if
    rstout = nsim+1
  end if
  turnoff(:) = .true.
  do kk=(polcalc*((nsim/polcalc)+2)),nequil,-polcalc
    if ((kk.le.nsim).AND.(kk.gt.nequil).AND.(mod(kk,polcalc).eq.0)) then
      turnoff(3) = .false.
      exit
    end if
  end do
  do kk=(holescalc*((nsim/holescalc)+2)),nequil,-holescalc
    if ((kk.le.nsim).AND.(kk.gt.nequil).AND.(mod(kk,holescalc).eq.0)) then
      turnoff(5) = .false.
      exit
    end if
  end do
  k = 0
  do kk=(savcalc*((nsim/savcalc)+2)),nequil,-savcalc
    if ((kk.le.nsim).AND.(kk.gt.nequil).AND.(mod(kk,savcalc).eq.0)) then
      turnoff(13) = .false.
      k = k + 1
      if (k.ge.savreq%instfreq) then
        turnoff(24) = .false.
        exit
      end if
    end if
  end do
  do kk=(covcalc*((nsim/covcalc)+2)),nequil,-covcalc
    if ((kk.le.nsim).AND.(kk.gt.nequil).AND.(mod(kk,covcalc).eq.0)) then
      turnoff(15) = .false.
      exit
    end if
  end do
  do kk=(phout*((nsim/phout)+2)),nequil,-phout
    if ((kk.le.nsim).AND.(kk.gt.nequil).AND.(mod(kk,phout).eq.0)) then
      turnoff(20) = .false.
      exit
    end if
  end do
  do kk=(xyzout*((nsim/xyzout)+2)),nequil,-xyzout
    if ((kk.le.nsim).AND.(kk.gt.nequil).AND.(mod(kk,xyzout).eq.0)) then
      turnoff(21) = .false.
      exit
    end if
  end do
  if ((particlenumcalc.le.nsim).AND.(ens%flag.ne.5).AND.(ens%flag.ne.6)) then
    write(ilog,*) 'Warning. Particle number analysis is only available &
 &if grand or semi-grand ensembles are used. Disabled.'
    particlenumcalc = nsim+1
  end if
  if ((turnoff(3).EQV..true.).AND.(polcalc.le.nsim)) then
    write(ilog,*) 'Warning. With chosen settings, no linear polymeric analysis &
 &is performed (',nsim,'/',nequil,'/',polcalc,'). Disabled.'
    polcalc = nsim+1
  end if
  if ((turnoff(5).EQV..true.).AND.(holescalc.le.nsim)) then
    write(ilog,*) 'Warning. With chosen settings, no void space analysis &
 &is performed (',nsim,'/',nequil,'/',holescalc,'). Disabled.'
    holescalc = nsim+1
  end if
  if ((turnoff(13).EQV..true.).AND.(savcalc.le.nsim)) then
    write(ilog,*) 'Warning. With chosen settings, no solvent accessibility analysis &
 &is performed (',nsim,'/',nequil,'/',savcalc,'). Disabled.'
    savcalc = nsim+1
  end if
  if ((turnoff(15).EQV..true.).AND.(covcalc.le.nsim)) then
    write(ilog,*) 'Warning. With chosen settings, no torsional signal trains &
 &are produced (',nsim,'/',nequil,'/',covcalc,'). Disabled.'
    covcalc = nsim+1
  end if
  if ((turnoff(20).EQV..true.).AND.(phout.le.nsim)) then
    write(ilog,*) 'Warning. With chosen settings, no titration state printout &
 &is produced (',nsim,'/',nequil,'/',phout,'). Disabled.'
    phout = nsim+1
  end if
  if ((turnoff(21).EQV..true.).AND.(xyzout.le.nsim)) then
    write(ilog,*) 'Warning. With chosen settings, no structural (trajectory) output &
 &is produced (',nsim,'/',nequil,'/',xyzout,'). Disabled.'
    xyzout = nsim+1
  end if
  if ((turnoff(24).EQV..true.).AND.(savcalc.le.nsim).AND.(savreq%instfreq*savcalc.le.nsim)) then
    write(ilog,*) 'Warning. With chosen settings, no instantaneous solvent-accessibility output &
 &is produced (',nsim,'/',nequil,'/',savcalc,'/',savreq%instfreq,'). Disabled.'
    savreq%instfreq = nsim+1
  end if
!
! sanity checks for available molecule types -> polymer stuff
  logdummy = .false.
  do i=1,nmoltyp
    if (do_pol(i).EQV..true.) then
      logdummy = .true.
      exit
    end if
  end do
  if (logdummy.EQV..false.) then
    write(ilog,*) 'WARNING: No molecule type with internal degrees o&
 &f freedom. Turning off all type-resolved analysis.'
    polcalc = nsim+1
    rhcalc = nsim+1
    holescalc = nsim+1
    segcalc = nsim+1
    covcalc = nsim+1
    angcalc = nsim+1
    torlccalc = nsim+1
    torout = nsim+1
  end if
! 
! sanity checks for available molecule types -> phi/psi stuff
  logdummy = .false.
  do i=1,nmoltyp
    if (do_tors(i).EQV..true.) then
      logdummy = .true.
      exit
    end if
  end do
  if (logdummy.EQV..false.) then
    write(ilog,*) 'WARNING: No molecule type with peptide-like torsi&
 &onal degrees of freedom. Turning off all dependent analysis.'
    segcalc = nsim+1
    angcalc = nsim+1
    dsspcalc = nsim+1
  end if
!
! sanity check for lack of solutes
  if (nsolutes.le.0) then
    write(ilog,*) 'Warning. No solute molecules present. Disabling all contact and cluster&
 & analyses.'
    contactcalc = nsim+1
  end if
!
! sanity check for structure alignment
  if ((pdb_analyze.EQV..false.).AND.(align%yes.EQV..true.).AND.(align%calc.le.nsim)) then
    write(ilog,*) 'Warning. Structure alignment (FMCSC_ALIGNCALC) is only supported in trajectory &
 &analysis runs. Disabled.'
    align%yes = .false.
    align%calc = nsim + 1
  end if
!
! sanity check against bad requested residues for torsional analysis
  if ((nrspecrama.gt.0).AND.(angcalc.le.nsim)) then
    k = 1
    do while (k.gt.0)
      k = 0
      do i=1,nrspecrama
        logdummy = .false.
        if ((specrama(i).le.nseq).AND.(specrama(i).gt.0)) then
          if((do_tors(moltypid(molofrs(specrama(i)))).EQV..false.)&
 & .OR.(fline(specrama(i)).le.0).OR.(yline(specrama(i)).le.0).OR.&
 & (specrama(i).eq.rsmol(molofrs(specrama(i)),1)).OR.&
 & (specrama(i).eq.rsmol(molofrs(specrama(i)),2))) then
            logdummy = .true.
          end if
        else
          logdummy = .true.
        end if
        if (logdummy.EQV..true.) then
          write(ilog,*)
          write(ilog,*) 'WARNING: Bad spec. for residue-specific Ram&
 &achandran analysis. Removing requested residue # ',specrama(i),'.'
          k = k + 1
          do j=i+1,nrspecrama
            specrama(j-1) = specrama(j)
          end do
          nrspecrama = nrspecrama - 1
          exit
        end if
      end do
    end do
  end if
!
! sanity check against bad requested molecule types
  if ((nrmolrama.gt.0).AND.(angcalc.le.nsim)) then
    k = 1
    do while (k.gt.0)
      k = 0
      do i=1,nrmolrama
        logdummy = .false.
        if ((molrama(i).gt.0).AND.(molrama(i).le.nangrps)) then
          if (do_tors(moltypid(molangr(molrama(i),1))).EQV..false.) then
            logdummy = .true.
          end if
        else
          logdummy = .true.
        end if
        if (logdummy.EQV..true.) then
          write(ilog,*)
          write(ilog,*) 'WARNING: Bad spec. for analysis group-specif&
 &ic Ramachandran analysis. Removing requ. group # ',molrama(i),'.'
          k = k + 1
          do j=i+1,nrmolrama
            molrama(j-1) = molrama(j)
          end do
          nrmolrama = nrmolrama - 1
          exit
        end if
      end do
    end do
  end if
!
! sanity check against bad requested contact analysis
  write(ilog,*)
  if ((contact_off.ge.nseq).AND.(contactcalc.le.nsim)) then
    write(ilog,*) 'WARNING: Bad spec. for offset for contact analys&
 &is (FMCSC_CONTACTOFF). In order to turn it off, it is &
 &recommended to adjust the frequency parameter.'
    write(ilog,*) 'Done that. (Ref.: Bad spec. is ',contact_off,')'
    contactcalc = nsim+1
  end if
!
! Sanity check against dipole request
!
  if ((dipcalc.le.nsim).AND.(use_POLAR.EQV..false.)) then
    write(ilog,*) 'WARNING: Dipole analysis requires polar potential&
 & to be turned on (assignment of charges). Turning off.'
    dipcalc = nsim+1
  end if
!
! Sanity check against bad diffraction calculation
!
  if (diffrcalc.le.nsim) then
    if (diffr_uax.EQV..true.) then
      if (sqrt(diffr_axis(1)**2 + &
 &             diffr_axis(2)**2 +&
 &             diffr_axis(3)**2).le.0.0) then
        write(ilog,*) 'Warning. Zero axis provided for diffraction ana&
 &lysis with constant axis. Turning off.'
        diffrcalc = nsim + 1
      end if
    end if
  end if
!
! Sanity check against pH write-out
!     
  if ((phfreq.le.0.0).AND.(phout.le.nsim)) then
    phout = nsim + 1
  end if
!
! Sanity check against holes calculation 
!
  if ((holescalc.le.nsim).AND.((bnd_type.eq.1).OR.(nmol.gt.1))) then
    holescalc = nsim+1
    write(ilog,*) 'WARNING: Void space analysis is not supported in periodic&
 & boundary conditions or for simulations with multiple molecules. Disabled.'
  end if
!
#ifdef LINK_LAPACK
! Sanity check against PCA
  if ((cstorecalc.le.nsim).AND.(pcamode.gt.1)) then
!
    if (cdis_crit.le.2) then
      write(ilog,*) 'Warning. Principal components are currently not supported&
 & for underlying data that are periodic (dihedral angles). Disabling PCA.'
      pcamode = 1
    end if
!
    if ((pcamode.eq.3).AND.(reduced_dim_clustering.gt.0)) then
      if ((reduced_dim_clustering.gt.clstsz).OR.&
 &        ((cdis_crit.eq.6).AND.(reduced_dim_clustering.ge.(cdofsbnds(2)-cdofsbnds(1)+1)))) then
        write(ilog,*) 'WARNING: Reduced dimensionality clustering after PCA implies&
 & that the chosen number of dimensions is smaller than the requested one (adjust&
 & clustering subset and/or setting for FMCSC_CREDUCEDIM). Turning off dimensionality&
 & reduction (PCA will still be performed).'
        reduced_dim_clustering = 0
      end if
    else
      reduced_dim_clustering = 0
    end if 
!
    if ((cdis_crit.eq.4).OR.(cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
      write(ilog,*) 'Warning. For PCA, locally adaptive weights are converted to static ones, which &
 &are then used to scale the data.'
      if (reduced_dim_clustering.gt.0) then
        write(ilog,*) 'This may lead to undesired effects for subsequent processing of data in reduced dimensionality.'
      end if
    end if
    if (((cdis_crit.eq.5).OR.(cdis_crit.eq.6)).AND.(align_for_clustering.EQV..true.)) then
      write(ilog,*) 'Warning. For PCA, the notion of pairwise alignment is lost for coordinate RMSD, i.e.,&
 & a trajectory prealigned to the last snapshot is used instead.'
      if (reduced_dim_clustering.gt.0) then
        write(ilog,*) 'This may lead to undesired effects for subsequent processing of data in reduced dimensionality.'
      end if
    end if
  end if
#endif
!
#ifdef ENABLE_MPI
!
! Sanity check against holes calculation 
!
  if ((holescalc.le.nsim).AND.(use_MPIAVG.EQV..true.)) then
    holescalc = nsim+1
    write(ilog,*) 'WARNING: Void space analysis is not supported in MP&
 &I averaging calculations. Disabled.'
  end if
!
! Sanity check against pH/MPI 
!
  if (phfreq.gt.0.0) then
    phfreq = 0.0
    phout = nsim + 1
    write(ilog,*) 'WARNING: pH moves are not supported in MPI calcul&
 &ations. Disabled.'
  end if
!
#endif
!
  if (select_frames.EQV..true.) then
    call strlims(frameidxfile,t1,t2)
    turnoff(:) = .true.
    do k=1,framecnt
      kk = framelst(k)
      modstep = mod(kk,enout)
      if (modstep.eq.0) turnoff(1) = .false.
      modstep = mod(kk,torout)
      if (modstep.eq.0) turnoff(2) = .false.
      modstep = mod(kk,polout)
      if (modstep.eq.0) turnoff(14) = .false.
      if (kk.le.nequil) cycle
      modstep = mod(kk,polcalc)
      if (modstep.eq.0) turnoff(3) = .false.
      modstep = mod(kk,particlenumcalc)
      if (modstep.eq.0) turnoff(4) = .false.
      modstep = mod(kk,holescalc)
      if (modstep.eq.0) turnoff(5) = .false.
      modstep = mod(kk,rhcalc)
      if (modstep.eq.0) turnoff(6) = .false.
      modstep = mod(kk,angcalc)
      if (modstep.eq.0) turnoff(7) = .false.
      modstep = mod(kk,intcalc)
      if (modstep.eq.0) turnoff(8) = .false.
      modstep = mod(kk,segcalc)
      if (modstep.eq.0) turnoff(9) = .false.
      modstep = mod(kk,dsspcalc)
      if (modstep.eq.0) turnoff(10) = .false.
      modstep = mod(kk,contactcalc)
      if (modstep.eq.0) turnoff(11) = .false.
      modstep = mod(kk,pccalc)
      if (modstep.eq.0) turnoff(12) = .false.
      modstep = mod(kk,savcalc)
      if (modstep.eq.0) turnoff(13) = .false.
      modstep = mod(kk,covcalc)
      if (modstep.eq.0) turnoff(15) = .false.
      modstep = mod(kk,dipcalc)
      if (modstep.eq.0) turnoff(16) = .false.
      modstep = mod(kk,diffrcalc)
      if (modstep.eq.0) turnoff(17) = .false.
      modstep = mod(kk,cstorecalc)
      if (modstep.eq.0) turnoff(18) = .false.
      modstep = mod(kk,torlccalc)
      if (modstep.eq.0) turnoff(19) = .false.
      modstep = mod(kk,phout)
      if (modstep.eq.0) turnoff(20) = .false.
      modstep = mod(kk,xyzout)
      if (modstep.eq.0) turnoff(21) = .false.
      modstep = mod(kk,rhcalc*sctcalc)
      if (modstep.eq.0) turnoff(22) = .false.
      modstep = mod(kk,contactcalc*clucalc)
      if (modstep.eq.0) turnoff(23) = .false.
      modstep = mod(kk,emcalc)
      if (modstep.eq.0) turnoff(24) = .false.
    end do
    if ((enout.le.nsim).AND.(turnoff(1).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling energy print-out due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_ENOUT.'
      enout = nsim+1
    end if
    if ((torout.le.nsim).AND.(turnoff(2).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling torsions print-out due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_TOROUT.'
      torout = nsim+1
    end if
    if ((polcalc.le.nsim).AND.(turnoff(3).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling simple polymeric analyses due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_POLCALC.'
      polcalc = nsim+1
    end if
    if ((particlenumcalc.le.nsim).AND.(turnoff(4).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling numbers histograms analysis due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_PARTICLENUMFREQ.'
      particlenumcalc = nsim+1
    end if
    if ((holescalc.le.nsim).AND.(turnoff(5).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling holes calculation due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_HOLESCALC.'
      holescalc = nsim+1
    end if
    if ((rhcalc.le.nsim).AND.(turnoff(6).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling quadratic polymeric analyses due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_RHCALC.'
      rhcalc = nsim+1
    end if
    if ((angcalc.le.nsim).AND.(turnoff(7).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling angular statistics collection due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_ANGCALC.'
      angcalc = nsim+1
    end if
    if ((intcalc.le.nsim).AND.(turnoff(8).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling internal coordinate analyses due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_INTCALC.'
      intcalc = nsim+1
    end if
    if ((segcalc.le.nsim).AND.(turnoff(9).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling torsional secondary structure analyses due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_SEGCALC.'
      segcalc = nsim+1
    end if
    if ((dsspcalc.le.nsim).AND.(turnoff(10).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling DSSP analyses due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_DSSPCALC.'
      dsspcalc = nsim+1
    end if
    if ((contactcalc.le.nsim).AND.(turnoff(11).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling contact analyses due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_CONTACTCALC.'
      contactcalc = nsim+1
    end if
    if ((pccalc.le.nsim).AND.(turnoff(12).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling all pair correlation analyses due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_PCCALC.'
      pccalc = nsim+1
    end if
    if ((savcalc.le.nsim).AND.(turnoff(13).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling solvent accessibility analyses due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_SAVCALC.'
      savcalc = nsim+1
    end if
    if ((polout.le.nsim).AND.(turnoff(14).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling inst. polymer print-out due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_POLOUT.'
      polout = nsim+1
    end if
    if ((covcalc.le.nsim).AND.(turnoff(15).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling torsional signal train print-out due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_COVCALC.'
      covcalc = nsim+1
    end if
    if ((dipcalc.le.nsim).AND.(turnoff(16).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling dipole analyses due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_DIPCALC.'
      dipcalc = nsim+1
    end if
    if ((diffrcalc.le.nsim).AND.(turnoff(17).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling fiber diffraction analysis due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_DIFFRCALC.'
      diffrcalc = nsim+1
    end if
    if ((cstorecalc.le.nsim).AND.(turnoff(18).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling structural clustering due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_CCOLLECT.'
      cstorecalc = nsim+1
    end if
    if ((torlccalc.le.nsim).AND.(turnoff(19).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling LCT analysis due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_TORLCCALC.'
      torlccalc = nsim+1
    end if
    if ((phout.le.nsim).AND.(turnoff(20).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling inst. protonation state print-out due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_PHOUT.'
      phout = nsim+1
    end if
    if ((xyzout.le.nsim).AND.(turnoff(21).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling trajectory output due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_XYZOUT.'
      xyzout = nsim+1
    end if
    if ((sctcalc*rhcalc.le.nsim).AND.(turnoff(22).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling Kratky analysis due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_SCATTERCALC.'
      sctcalc = nsim+1
    end if
    if ((clucalc*contactcalc.le.nsim).AND.(turnoff(23).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling molecular cluster analyses due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_CLUSTERCALC.'
      clucalc = nsim+1
    end if
    if ((emcalc.le.nsim).AND.(turnoff(24).EQV..true.)) then
      write(ilog,*) 'Warning. Disabling spatial density analysis due to lack of overlap &
 &between provided list of frames (',frameidxfile(t1:t2),') and setting for FMCSC_EMCALC.'
      clucalc = nsim+1
    end if
  end if
!
end
!
!--------------------------------------------------------------------------------------
!

