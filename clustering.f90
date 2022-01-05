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
! CONTRIBUTIONS: Nicolas Bloechliger                                       !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!-------------------------------------------------------------------------
!
subroutine store_for_clustering()
!
  use sequen
  use forces
  use clusters
  use fyoc
  use polypep
  use atoms
  use system
  use aminos
  use molecule
  use movesets
  use iounit
  use zmatrix
  use math
  use mcsums
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
#ifdef ENABLE_MPI
#include "mpif.h"
#endif
  integer kk,ttc,ttc2,i,rs,imol,kidx,cdofsz,cta
  logical afalse
  RTYPE d2v,mmtw,clutmp(calcsz)
#ifdef ENABLE_MPI
  integer clumpiavgtag,mstatus(MPI_STATUS_SIZE),masterrank,ierr
#endif
!
  afalse = .false.
!
  cstored = cstored + 1
  cta = cstored
#ifdef ENABLE_MPI
  masterrank = 0
  clumpiavgtag = 355
#endif
!
! transcribe breaks list (in-place)
  if ((ntbrks.ge.itbrklst).AND.(cstored.gt.1)) then
    if (nstep.gt.trbrkslst(itbrklst)) then
      ntbrks2 = ntbrks2 + 1
      trbrkslst(ntbrks2) = cstored - 1
      do while ((itbrklst.le.ntbrks).AND.(trbrkslst(itbrklst).lt.nstep))
        itbrklst = itbrklst + 1
        if (itbrklst.gt.ntbrks) exit
        if ((nstep.eq.nsim).AND.(trbrkslst(itbrklst).eq.nstep)) then
          ntbrks2 = ntbrks2 + 1
          trbrkslst(ntbrks2) = cstored
          itbrklst = itbrklst + 1
        end if
        if (itbrklst.gt.ntbrks) exit
      end do
    else if ((nstep.eq.nsim).AND.(trbrkslst(itbrklst).eq.nstep)) then
      ntbrks2 = ntbrks2 + 1
      trbrkslst(ntbrks2) = cstored
      itbrklst = itbrklst + 1
    end if
  end if
!
! torsional
  if ((cdis_crit.ge.1).AND.(cdis_crit.le.4)) then
    if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4)) then
!     may require backup 
      if ((dyn_integrator_ops(1).gt.0).AND.(use_dyn.EQV..true.).AND.(pdb_analyze.EQV..false.).AND.(fycxyz.eq.1)) then
        do imol=1,nmol
          dc_di(imol)%olddat(:,3) = dc_di(imol)%olddat(:,2)
        end do
      end if
      ! update com and populate olddat(:,2)
      do imol=1,nmol
        call update_rigidm(imol)
      end do
      call cart2int_I()
    end if
    cdofsz = calcsz/2
    if (cdis_crit.eq.1) cdofsz = calcsz
    if (cdis_crit.eq.4) cdofsz = calcsz/3
!
    clutmp(:) = 0.0
!   all torsions
    ttc = 0
    ttc2 = 1
    do rs=1,nseq
      imol = molofrs(rs)
!     omega
      if (wline(rs).gt.0) then
        ttc = ttc + 1
        mmtw = 1.0
        if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4)).AND.(wnr(rs).gt.0)) mmtw = dc_di(imol)%olddat(wnr(rs),2)
        if (cdofset(ttc2,1).eq.ttc) then
          if (cdis_crit.eq.1) then
            clutmp(ttc2) = omega(rs)
          else if (cdis_crit.eq.2) then
            clutmp(2*ttc2-1) = omega(rs)
            clutmp(2*ttc2) = mmtw
          else if (cdis_crit.eq.3) then
            clutmp(2*ttc2-1) = sin(omega(rs)/RADIAN)
            clutmp(2*ttc2) = cos(omega(rs)/RADIAN)
          else if (cdis_crit.eq.4) then
            clutmp(3*ttc2-2) = sin(omega(rs)/RADIAN)
            clutmp(3*ttc2-1) = cos(omega(rs)/RADIAN)
            clutmp(3*ttc2) = mmtw
          end if
          ttc2 = ttc2 + 1
          if (ttc2.gt.cdofsz) exit
        end if
      end if
!     phi
      if (fline(rs).gt.0) then
        ttc = ttc + 1
        mmtw = 1.0
        if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4)).AND.(fnr(rs).gt.0)) mmtw = dc_di(imol)%olddat(fnr(rs),2)
        if (cdofset(ttc2,1).eq.ttc) then
          if (cdis_crit.eq.1) then
            clutmp(ttc2) = phi(rs)
          else if (cdis_crit.eq.2) then
            clutmp(2*ttc2-1) = phi(rs)
            clutmp(2*ttc2) = mmtw
          else if (cdis_crit.eq.3) then
            clutmp(2*ttc2-1) = sin(phi(rs)/RADIAN)
            clutmp(2*ttc2) = cos(phi(rs)/RADIAN)
          else if (cdis_crit.eq.4) then
            clutmp(3*ttc2-2) = sin(phi(rs)/RADIAN)
            clutmp(3*ttc2-1) = cos(phi(rs)/RADIAN)
            clutmp(3*ttc2) = mmtw
          end if
          ttc2 = ttc2 + 1
          if (ttc2.gt.cdofsz) exit
        end if
      end if
!     psi
      if (yline(rs).gt.0) then
        ttc = ttc + 1
        mmtw = 1.0
        if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4)).AND.(ynr(rs).gt.0)) mmtw = dc_di(imol)%olddat(ynr(rs),2)
        if (cdofset(ttc2,1).eq.ttc) then
          if (cdis_crit.eq.1) then
            clutmp(ttc2) = psi(rs)
          else if (cdis_crit.eq.2) then
            clutmp(2*ttc2-1) = psi(rs)
            clutmp(2*ttc2) = mmtw
          else if (cdis_crit.eq.3) then
            clutmp(2*ttc2-1) = sin(psi(rs)/RADIAN)
            clutmp(2*ttc2) = cos(psi(rs)/RADIAN)
          else if (cdis_crit.eq.4) then
            clutmp(3*ttc2-2) = sin(psi(rs)/RADIAN)
            clutmp(3*ttc2-1) = cos(psi(rs)/RADIAN)
            clutmp(3*ttc2) = mmtw
          end if
          ttc2 = ttc2 + 1
          if (ttc2.gt.cdofsz) exit
        end if
      end if
!     nucleic acid bb angles
      do kk = 1,nnucs(rs)
        ttc = ttc + 1
        mmtw = 1.0
        if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4)).AND.(nucsnr(kk,rs).gt.0)) mmtw = dc_di(imol)%olddat(nucsnr(kk,rs),2)
        if (cdofset(ttc2,1).eq.ttc) then
          if (cdis_crit.eq.1) then
            clutmp(ttc2) = nucs(kk,rs)
          else if (cdis_crit.eq.2) then
            clutmp(2*ttc2-1) = nucs(kk,rs)
            clutmp(2*ttc2) = mmtw
          else if (cdis_crit.eq.3) then
            clutmp(2*ttc2-1) = sin(nucs(kk,rs)/RADIAN)
            clutmp(2*ttc2) = cos(nucs(kk,rs)/RADIAN)
          else if (cdis_crit.eq.4) then
            clutmp(3*ttc2-2) = sin(nucs(kk,rs)/RADIAN)
            clutmp(3*ttc2-1) = cos(nucs(kk,rs)/RADIAN)
            clutmp(3*ttc2) = mmtw
          end if
          ttc2 = ttc2 + 1
          if (ttc2.gt.cdofsz) exit
        end if
      end do
      if (ttc2.gt.cdofsz) exit
!     sugar bond in nucleotides
      if (seqpolty(rs).eq.'N') then
        kidx = nucsline(6,rs)
        ttc = ttc + 1
        mmtw = 1.0
        if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4)).AND.(nnucs(rs).gt.1)) then
          if (nucsnr(nnucs(rs)-1,rs).gt.0) mmtw = dc_di(imol)%olddat(nucsnr(nnucs(rs)-1,rs),2)
        end if
!       note the cheating on the inertial mass for the sugar bond
        if (cdofset(ttc2,1).eq.ttc) then
          if (cdis_crit.eq.1) then
            clutmp(ttc2) = ztor(kidx)
          else if (cdis_crit.eq.2) then
            clutmp(2*ttc2-1) = ztor(kidx)
            clutmp(2*ttc2) = mmtw
          else if (cdis_crit.eq.3) then
            clutmp(2*ttc2-1) = sin(ztor(kidx)/RADIAN)
            clutmp(2*ttc2) = cos(ztor(kidx)/RADIAN)
          else if (cdis_crit.eq.4) then
            clutmp(3*ttc2-2) = sin(ztor(kidx)/RADIAN)
            clutmp(3*ttc2-1) = cos(ztor(kidx)/RADIAN)
            clutmp(3*ttc2) = mmtw
          end if
          ttc2 = ttc2 + 1
          if (ttc2.gt.cdofsz) exit
        end if
      end if
!     chi angles
      do kk = 1,nchi(rs)
        ttc = ttc + 1
        mmtw = 1.0
        if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4)).AND.(chinr(kk,rs).gt.0)) mmtw = dc_di(imol)%olddat(chinr(kk,rs),2)
        if (cdofset(ttc2,1).eq.ttc) then
          if (cdis_crit.eq.1) then
            clutmp(ttc2) = chi(kk,rs)
          else if (cdis_crit.eq.2) then
            clutmp(2*ttc2-1) = chi(kk,rs)
            clutmp(2*ttc2) = mmtw
          else if (cdis_crit.eq.3) then
            clutmp(2*ttc2-1) = sin(chi(kk,rs)/RADIAN)
            clutmp(2*ttc2) = cos(chi(kk,rs)/RADIAN)
          else if (cdis_crit.eq.4) then
            clutmp(3*ttc2-2) = sin(chi(kk,rs)/RADIAN)
            clutmp(3*ttc2-1) = cos(chi(kk,rs)/RADIAN)
            clutmp(3*ttc2) = mmtw
          end if
          ttc2 = ttc2 + 1
          if (ttc2.gt.cdofsz) exit
        end if
      end do
      if (ttc2.gt.cdofsz) exit
    end do
    if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4)) then
!     may require restore
      if ((dyn_integrator_ops(1).gt.0).AND.(use_dyn.EQV..true.).AND.(pdb_analyze.EQV..false.).AND.(fycxyz.eq.1)) then
        do imol=1,nmol
          dc_di(imol)%olddat(:,2) = dc_di(imol)%olddat(:,3)
        end do
      end if
    end if
! xyz atom set
  else if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6).OR.(cdis_crit.eq.10)) then
    do i=1,clstsz/3
      clutmp(3*i-2) = x(cdofset(i,1))
      clutmp(3*i-1) = y(cdofset(i,1))
      clutmp(3*i) = z(cdofset(i,1))
    end do
    if (cdis_crit.eq.10) clutmp((clstsz+1):calcsz) = 1.0 ! currently no locally adaptive weight is collected live
! internal distance vector
  else if ((cdis_crit.ge.7).AND.(cdis_crit.le.9)) then
    do i=1,clstsz
      call dis_bound(cdofset(i,1),cdofset(i,2),d2v)
      clutmp(i) = sqrt(d2v)
    end do
    if (cdis_crit.eq.9) clutmp((clstsz+1):calcsz) = 1.0 ! currently no locally adaptive weight is collected live
  end if
!
#ifdef ENABLE_MPI
  if (use_MPIAVG.EQV..true.) then
    if (myrank.eq.masterrank) then
      cludata(:,cstored) = clutmp(:)
      do i=1,mpi_nodes-1
        call MPI_RECV(clutmp,calcsz,MPI_RTYPE,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.clumpiavgtag) then
          write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from master. Expected ',clumpiavgtag,'.'
          call fexit()
        end if
        cludata(:,i*cmaxsnaps+cstored) = clutmp(:)
      end do
    else
      call MPI_SEND(clutmp,calcsz,MPI_RTYPE,masterrank,clumpiavgtag,MPI_COMM_WORLD,ierr)
    end if
  else
    cludata(:,cstored) = clutmp(:)
  end if
#else
  cludata(:,cstored) = clutmp(:)
#endif
!
end
!
!--------------------------------------------------------------------------------------------------
!
! allow centering, variance normalization, and smoothing prior to clustering
!
subroutine preprocess_cludata()
!
  use clusters
  use iounit
!
  implicit none
!
  RTYPE cmean,icst,cvar
  integer effdd,dd
  RTYPE, ALLOCATABLE:: cvx(:),cvx2(:)
!
  if (csmoothie.ge.(cstored/2-1)) then
    write(ilog,*) 'Warning. Selected smoothing window size is too large (FMCSC_CSMOOTHORDER ',csmoothie,').'
    csmoothie = cstored/2-1
    write(ilog,*) 'Value is adjusted to (roughly) half the data set length (',csmoothie,').'
  end if
!
  if (cprepmode.le.0) return
!
  if ((cdis_crit.eq.1).OR.(cdis_crit.eq.2)) then
    write(ilog,*) 'Warning. Disabling data preprocessing for structural clustering due to use of &
 &explicit dihedral angles (periodic quantities with fixed range centered at 0.0).'
    return
  end if
!
  if ((cdis_crit.eq.3).OR.(cdis_crit.eq.4)) then
    write(ilog,*) 'Warning. Data preprocessing for structural clustering may not be particularly &
 &meaningful when performed on sine and cosine terms of dihedral angles.'
  end if
  if ((align_for_clustering.EQV..true.).AND.((cdis_crit.eq.5).OR.(cdis_crit.eq.6).OR.(cdis_crit.eq.10))) then
    if (cprepmode.gt.1) then
      write(ilog,*) 'Warning. Disabling data preprocessing for structural clustering except centering due to use of &
 &coordinate RMSD with alignment. Use a prealigned trajectory instead.'
      cprepmode = 1
    end if
  end if
!
  icst = 1.0/(1.0*cstored)
!
  effdd = clstsz
  if (cdis_crit.eq.4) effdd = 3*effdd/2
  if (cprepmode.ge.3) then
    allocate(cvx(cstored))
    allocate(cvx2(cstored))
  end if
  do dd=1,effdd
    if ((cdis_crit.eq.4).AND.(mod(dd,3).eq.0)) cycle ! skip weight dimension
    if (cprepmode.ne.3) then
      cmean = icst*sum(cludata(dd,1:cstored))
      cludata(dd,1:cstored) = cludata(dd,1:cstored) - cmean
    end if
    if ((cprepmode.eq.2).OR.(cprepmode.eq.5)) then
      cvar = 1.0/sqrt(sum(cludata(dd,1:cstored)**2)/(1.0*(cstored-1)))
      cludata(dd,1:cstored) = cvar*cludata(dd,1:cstored)
    end if
    if (cprepmode.ge.3) then
      cvx(:) = cludata(dd,1:cstored)
      call vsmooth(cvx,cstored,csmoothie,cvx2)
      cludata(dd,1:cstored) = cvx2(:)
    end if
  end do

!
  if (allocated(cvx).EQV..true.) deallocate(cvx)
  if (allocated(cvx2).EQV..true.) deallocate(cvx2)
!
end
!
!---------------------------------------------------------------------------------------------------
!
! allows weights to be repopulated with different types of information
! 0: leave as is
! 1: inverse standard deviation 
! 2: static weights using ACF at fixed lag
! 3: square root of inverse std dev and ACF at fixed lag
! 4: rate of transitions across global mean
! 5: same as 4 with smoothing for weights generation
! 6: sqrt(2*4)
! 7: sqrt(2*5)
! 8: rate of transitions across separators defined by histogram minima
! 9: same as 8 with smoothing for weights generation
!
subroutine repopulate_weights()
!
  use clusters
  use iounit
  use math
!
  implicit none
!
  integer dd,k,i,j,ii,ntrans,ncvs,nbins,nbinsm10,nbinsblk,effdd,wpos,lastk
  RTYPE cmean,cvar,icst,cacf,cm1,cm2,normer
  RTYPE, ALLOCATABLE:: cvx(:),cvx2(:),cvw(:),cvlst(:)
  integer, ALLOCATABLE:: cvi(:),ixlst(:)
!
 677 format(1000(g12.5,1x))
!
  if (cchangeweights.eq.0) return
!
  if (cstored.le.4) then
    write(ilog,*) 'Warning. Refusing to adjust weights for clustering due to insufficient data length.'
    return
  end if
!
  if (cchangeweights.eq.2) then
    if (cwacftau.ge.(cstored/2-1)) then
      write(ilog,*) 'Warning. Selected lag time for ACF weights is too large (FMCSC_CLAGTIME ',cwacftau,').'
      cwacftau = cstored/2-1
      write(ilog,*) 'Value is adjusted to roughly half the data set length (',cwacftau,').'
    end if
  end if
!
  nbins = 100
  nbinsm10 = nbins - 10
  nbinsblk = 3
  if (cchangeweights.ge.8) then
    allocate(cvi(nbins))
    allocate(cvlst(nbinsm10))
    allocate(ixlst(nbinsm10))
  end if
  allocate(cvx(cstored))
!
  if ((cchangeweights.eq.2).OR.(cchangeweights.eq.3).OR.(cchangeweights.eq.6).OR.(cchangeweights.eq.7)) then
    allocate(cvw(clstsz))
    cvw(:) = 0.0
  end if
  if ((cchangeweights.eq.1).OR.(cchangeweights.eq.3).OR.&
 &    (cchangeweights.eq.5).OR.(cchangeweights.eq.7).OR.(cchangeweights.eq.9)) allocate(cvx2(cstored))
!
  if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    if (cwwindowsz.ge.cstored) then
      write(ilog,*) 'Warning. Selected window size for locally adaptive weights is too large (FMCSC_CWINDOWSIZE ',cwwindowsz,').'
      cwwindowsz = cstored/2
      write(ilog,*) 'Value is adjusted to roughly half the data set length (',cwwindowsz,').'
    end if
  end if
  if (cchangeweights.eq.0) return
!
! locally adaptive weights of all types
  if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    lastk = 2
    if (cdis_crit.ge.9) lastk = 1
!   to avoid having to rewrite the same code several times, we pull the selector inside
    ii = cwwindowsz/2
    icst = 1.0/(1.0*cstored)
    effdd = clstsz
    if (cdis_crit.eq.4) effdd = effdd/2
    do dd=1,effdd
      do k=1,lastk
        if (cdis_crit.eq.4) then
          cvx(:) = cludata(3*dd-3+k,1:cstored)
          wpos = 3*dd
        else if (cdis_crit.eq.2) then
          if (k.eq.1) cvx(:) = sin(cludata(2*dd-1,1:cstored)/RADIAN)
          if (k.eq.2) cvx(:) = cos(cludata(2*dd-1,1:cstored)/RADIAN)
          wpos = 2*dd
        else
          cvx(:) = cludata(dd,1:cstored) 
          wpos = dd + clstsz
        end if
!       if so desired, populate temp vector of static weights with ACF values
        if ((cchangeweights.eq.2).OR.(cchangeweights.eq.3).OR.(cchangeweights.eq.6).OR.(cchangeweights.eq.7)) then
          call vacf_fixtau(cvx,cstored,cwacftau,cacf)
          cvw(dd) = max(cvw(dd),cacf)
        end if
!       do nothing more for static ACF weights
        if (cchangeweights.eq.2) then
!       locally adaptive measure based on variance possibly combined with ACF
        else if ((cchangeweights.eq.1).OR.(cchangeweights.eq.3)) then
          i = cwwindowsz
          call vmsf_window(cvx,cstored,i,cvx2)
          do i=1,cstored
            if (k.eq.1) cludata(wpos,i) = max(0.0,cvx2(i))
            if (k.eq.2) cludata(wpos,i) = cludata(wpos,i) + max(0.0,cvx2(i))
            if ((k.eq.lastk).AND.(cludata(wpos,i).gt.0.0)) cludata(wpos,i) = sqrt(1.0/(cludata(wpos,i)))
          end do
!       locally adaptive measure based on crossings of the global mean, possibly estimated from smoothed data,
!       possibly combined with ACF
        else if ((cchangeweights.ge.4).AND.(cchangeweights.le.7)) then
          if ((cchangeweights.eq.5).OR.(cchangeweights.eq.7)) then
            call vsmooth(cvx,cstored,csmoothie,cvx2)
            cvx(:) = cvx2(:)
          end if
          cmean = icst*sum(cvx)
          ntrans = 0
          do i=2,cwwindowsz
            if (((cvx(i-1).gt.cmean).AND.(cvx(i).lt.cmean)).OR.((cvx(i-1).lt.cmean).AND.(cvx(i).gt.cmean))) then
              ntrans = ntrans + 1
            end if
          end do
          if (k.eq.1) cludata(wpos,ii) = 1.0*ntrans+cdynwbuf
          if (k.eq.2) cludata(wpos,ii) = 1.0*min(cludata(wpos,ii),ntrans+cdynwbuf)
          do i=cwwindowsz+1,cstored
            if (((cvx(i-1).gt.cmean).AND.(cvx(i).lt.cmean)).OR.((cvx(i-1).lt.cmean).AND.(cvx(i).gt.cmean))) ntrans = ntrans + 1
            if (((cvx(i-cwwindowsz).gt.cmean).AND.(cvx(i-cwwindowsz+1).lt.cmean)).OR.&
 &              ((cvx(i-cwwindowsz).lt.cmean).AND.(cvx(i-cwwindowsz+1).gt.cmean))) ntrans = ntrans - 1
!
            if (k.eq.1) cludata(wpos,i-ii) = 1.0*ntrans+cdynwbuf
            if (k.eq.2) cludata(wpos,i-ii) = 1.0*min(cludata(wpos,i-ii),ntrans+cdynwbuf)
          end do
!       locally adaptive measure based on passes of separators, possibly of smoothed data
        else if ((cchangeweights.eq.8).OR.(cchangeweights.eq.9)) then
          cvi(:) = 0
          call vautohist(cvx,cstored,nbins,cvi,cm1,cm2)
          ncvs = 0
          ixlst(:) = 0
          call vminima_int(cvi,nbins,ncvs,ixlst,nbinsm10,nbinsblk)
          do i=1,ncvs
            cvlst(i) = cm1 + (ixlst(i)-0.5)*cm2
          end do
          if (ncvs.le.0) then
            ncvs = 1
            cvlst(ncvs) = icst*sum(cvx)
          end if
          if (cchangeweights.eq.9) then
            call vsmooth(cvx,cstored,csmoothie,cvx2)
            cvx(:) = cvx2(:)
          end if
!         we need to penalize low transition rates based on separators that are within a low likelihood tail
!         cm1 should correspond to the maximal value of the smaller fractional population for any considered separator
          cm1 = 0.0
          do j=1,ncvs
            if (ixlst(j).le.0) then
              cm1 = max(0.5,cm1)
            else
              cm1 = max(cm1,min(sum(cvi(1:ixlst(j)))/(1.0*cstored),sum(cvi((ixlst(j)+1):100))/(1.0*cstored)))
            end if
          end do
          ntrans = 0
          do i=2,cwwindowsz
            do j=1,ncvs
              if (((cvx(i-1).gt.cvlst(j)).AND.(cvx(i).lt.cvlst(j))).OR.((cvx(i-1).lt.cvlst(j)).AND.(cvx(i).gt.cvlst(j)))) then
                ntrans = ntrans + 1
                exit
              end if
            end do
          end do
          if (k.eq.1) cludata(wpos,ii) = (1.0*ntrans+cdynwbuf)/cm1
          if (k.eq.2) cludata(wpos,ii) = min(cludata(wpos,ii),(1.0*ntrans+cdynwbuf)/cm1)
          do i=cwwindowsz+1,cstored
            do j=1,ncvs
              if (((cvx(i-1).gt.cvlst(j)).AND.(cvx(i).lt.cvlst(j))).OR.((cvx(i-1).lt.cvlst(j)).AND.(cvx(i).gt.cvlst(j)))) then
                ntrans = ntrans + 1
                exit
              end if
            end do
            do j=1,ncvs
              if (((cvx(i-cwwindowsz).gt.cvlst(j)).AND.(cvx(i-cwwindowsz+1).lt.cvlst(j))).OR.&
 &                ((cvx(i-cwwindowsz).lt.cvlst(j)).AND.(cvx(i-cwwindowsz+1).gt.cvlst(j)))) then
                ntrans = ntrans - 1
                exit
              end if
            end do
            if (k.eq.1) cludata(wpos,i-ii) = (1.0*ntrans+cdynwbuf)/cm1
            if (k.eq.2) cludata(wpos,i-ii) = min(cludata(wpos,i-ii),(1.0*ntrans+cdynwbuf)/cm1)
          end do
        else
          write(ilog,*) 'Fatal. Weights of type ',cchangeweights,' are not supported at the moment. &
 &This is likely to be an omission bug.'
          call fexit()
        end if
      end do
    end do
    do dd=1,effdd
      if (cdis_crit.eq.4) then
        wpos = 3*dd
      else if (cdis_crit.eq.2) then
        wpos = 2*dd
      else
        wpos = dd+clstsz
      end if
!     complete those data with access to incomplete windows
      if (cchangeweights.gt.3) then  
        cludata(wpos,1:(ii-1)) = cludata(wpos,ii)
        cludata(wpos,(cstored-ii+1):cstored) = cludata(wpos,(cstored-ii))
      end if
      if ((cchangeweights.ge.4).AND.(cchangeweights.le.9)) then
        cludata(wpos,1:cstored) = 1.0/cludata(wpos,1:cstored)
      end if
      if ((cchangeweights.eq.3).OR.(cchangeweights.eq.6).OR.(cchangeweights.eq.7)) then
        cludata(wpos,1:cstored) = sqrt(cvw(dd)*cludata(wpos,1:cstored))
      else if (cchangeweights.eq.2) then
        cludata(wpos,1:cstored) = cvw(dd)
      end if
!      write(*,*) dd,icst*sum(cludata(wpos,1:cstored)),&
! &icst*sum((cludata(wpos,1:cstored)-icst*sum(cludata(wpos,1:cstored)))**2)
    end do
!    do i=1,cstored
!      write(0,677) cludata(clstsz,i),cludata(calcsz,i)
!    end do
! static weights are only for interatomic distances at the moment
  else if (cdis_crit.eq.8) then
    icst = 1.0/(1.0*cstored)
    if ((cchangeweights.eq.2).OR.(cchangeweights.eq.3).OR.(cchangeweights.eq.6).OR.(cchangeweights.eq.7)) then
      do dd=1,clstsz
        cvx(:) = cludata(dd,1:cstored)
        call vacf_fixtau(cvx,cstored,cwacftau,cacf)
        cvw(dd) = max(cvw(dd),cacf)
      end do
    end if
    if ((cchangeweights.eq.1).OR.(cchangeweights.eq.3)) then
      do dd=1,clstsz
        cvx(:) = cludata(dd,1:cstored)
        cmean = icst*sum(cvx)
        cvar = icst*sum((cvx(:)- cmean)**2)
        if (cvar.le.0.0) then
          cl_imvec(dd) = 0.0
          cycle
        end if
        cl_imvec(dd) = sqrt(1.0/cvar)
        if (cchangeweights.eq.3) cl_imvec(dd) = sqrt(cl_imvec(dd)*cvw(dd))
      end do
!   based on ACF
    else if (cchangeweights.eq.2) then
      cl_imvec(1:clstsz) = cvw(1:clstsz)
!   based on passes of global mean
    else if ((cchangeweights.ge.4).AND.(cchangeweights.le.7)) then
      do dd=1,clstsz
        cvx(:) = cludata(dd,1:cstored)
        if ((cchangeweights.eq.5).OR.(cchangeweights.eq.7)) then
          call vsmooth(cvx,cstored,csmoothie,cvx2)
          cvx(:) = cvx2(:)
        end if
        cmean = icst*sum(cvx)
        ntrans = 0
        do i=2,cstored
          if (((cvx(i-1).gt.cmean).AND.(cvx(i).lt.cmean)).OR.((cvx(i-1).lt.cmean).AND.(cvx(i).gt.cmean))) then
            ntrans = ntrans + 1
          end if
        end do
        cl_imvec(dd) = 1.0/(1.0*(ntrans+cdynwbuf))
        if ((cchangeweights.eq.6).OR.(cchangeweights.eq.7)) cl_imvec(dd) = sqrt(cl_imvec(dd)*cvw(dd))
      end do
!   based on passes of separators, possibly of smoothed data
    else if ((cchangeweights.eq.8).OR.(cchangeweights.eq.9)) then
      do dd=1,clstsz
        cvx(:) = cludata(dd,1:cstored)
        cvi(:) = 0
        call vautohist(cvx,cstored,nbins,cvi,cm1,cm2)
        ncvs = 0
        ixlst(:) = 0
        call vminima_int(cvi,nbins,ncvs,ixlst,nbinsm10,nbinsblk)
        do i=1,ncvs
          cvlst(i) = cm1 + (ixlst(i)-0.5)*cm2
        end do
        if (ncvs.le.0) then
          ncvs = 1
          cvlst(ncvs) = icst*sum(cvx)
        end if
        if (cchangeweights.eq.9) then
          call vsmooth(cvx,cstored,csmoothie,cvx2)
          cvx(:) = cvx2(:)
        end if
        ntrans = 0
        do i=2,cstored
          do k=1,ncvs
            if (((cvx(i-1).gt.cvlst(k)).AND.(cvx(i).lt.cvlst(k))).OR.((cvx(i-1).lt.cvlst(k)).AND.(cvx(i).gt.cvlst(k)))) then
              ntrans = ntrans + 1
              exit
            end if
          end do
          if (i.eq.cwwindowsz) ntrans = 0
        end do
!       we need to penalize low transition rates based on separators that are within a low likelihood tail
!       cm1 should correspond to the maximal value of the smaller fractional population for any considered separator
        cm1 = 0.0
        do k=1,ncvs
          if (ixlst(k).le.0) then
            cm1 = max(0.5,cm1)
          else
            cm1 = max(cm1,min(sum(cvi(1:ixlst(k)))/(1.0*cstored),sum(cvi((ixlst(k)+1):100))/(1.0*cstored)))
          end if
        end do
        cl_imvec(dd) = cm1/(1.0*(ntrans+cdynwbuf))
      end do
    else
      write(ilog,*) 'Fatal. Weights of type ',cchangeweights,' are not supported for interatomic distances. &
 &This is an omission bug.'
      call fexit()
    end if
    normer = sum(cl_imvec(:))
    if (normer.le.0.0) then ! no variance in data
      cl_imvec(:) = 1.0
    else
      cl_imvec(:) = sqrt(cl_imvec(:)/normer)
    end if
  end if
!
  if (allocated(cvi).EQV..true.) deallocate(cvi)
  if (allocated(cvx).EQV..true.) deallocate(cvx)
  if (allocated(cvx2).EQV..true.) deallocate(cvx2)
  if (allocated(ixlst).EQV..true.) deallocate(ixlst)
  if (allocated(cvlst).EQV..true.) deallocate(cvlst)
  if (allocated(cvw).EQV..true.) deallocate(cvw)
!
end
!
!------------------------------------------------------------------------------------
!
subroutine do_clustering()
!
  use clusters
  use iounit
!
  implicit none
!
  integer exitcode
  RTYPE t1,t2
!
  exitcode = 0
!
  if (cstored.le.0) return
!
! preprocess data if so desired
  if ((cchangeweights.gt.0).OR.(cprepmode.gt.0)) write(ilog,*)
  call CPU_time(t1)
  call preprocess_cludata()
  call repopulate_weights()
  call CPU_time(t2)
  if ((cchangeweights.gt.0).OR.(cprepmode.gt.0)) then
    write(ilog,*) 'Time elapsed for data preprocessing for structural clustering: ',t2-t1, ' [s]'
  end if
!
  if (pcamode.gt.1) then
#ifdef LINK_LAPACK
    if ((cdis_crit.ge.3)) then
      call get_principal_components(calcsz,cstored,clstsz)
    end if
#endif
  end if
!
  if ((cmode.eq.1).OR.(cmode.eq.2)) then
    call CPU_time(t1)
    call leader_clustering(cmode,nstruccls)
    call CPU_time(t2)
    write(ilog,*) 'Time elapsed for clustering: ',t2-t1, ' [s]'
  else if (cmode.eq.5) then
    call CPU_time(t1)
    call birch_clustering(cmode,nstruccls)
    call CPU_time(t2)
    write(ilog,*) 'Time elapsed for clustering: ',t2-t1, ' [s]'
  else
    if (cmode.eq.4) then
      call CPU_time(t1)
      if (cprogindex.eq.1) then
#ifdef LINK_NETCDF
        if (read_nbl_from_nc.EQV..true.) then
          call read_nbl_nc()
        else
          call gennbl_for_clustering()
        end if
#else
        call gennbl_for_clustering()
#endif
        call CPU_time(t2)
        write(ilog,*) 'Time elapsed for neighbor list creation/reading: ',t2-t1, ' [s]'
        call CPU_time(t1)
        call gen_MST_from_nbl()
        call do_prog_index()
        call CPU_time(t2)
        write(ilog,*) 'Time elapsed for generating progress index data: ',t2-t1, ' [s]'
      else
        call birch_clustering(cmode,nstruccls)
        call CPU_time(t2)
        write(ilog,*) 'Time elapsed for clustering: ',t2-t1, ' [s]'
        call CPU_time(t1)
        call gen_MST_from_treeclustering()
        call CPU_time(t2)
        write(ilog,*) 'Time elapsed for MST building: ',t2-t1, ' [s]'
        call CPU_time(t1)
        call do_prog_index()
        call CPU_time(t2)
        write(ilog,*) 'Time elapsed for generating progress index data: ',t2-t1, ' [s]'
      end if
    else if (cmode.eq.3) then
      call CPU_time(t1)
#ifdef LINK_NETCDF
      if (read_nbl_from_nc.EQV..true.) then
        call read_nbl_nc()
      else
        call gennbl_for_clustering()
      end if
#else
      call gennbl_for_clustering()
#endif
      call CPU_time(t2)
      write(ilog,*) 'Time elapsed for neighbor list creation/reading: ',t2-t1, ' [s]'
      call CPU_time(t1)
      call hierarchical_clustering(nstruccls)
      call CPU_time(t2)
      write(ilog,*) 'Time elapsed for clustering: ',t2-t1, ' [s]'
    else
      write(ilog,*) 'Fatal. Encountered unknown clustering mode in do_clustering(...). This&
 & is a bug.'
      call fexit()
    end if
    if (exitcode.ne.0) then
      write(ilog,*) 'Fatal. Clustering algorithm in do_clustering(...) exited with an error. This&
 & is either a bug or a corrupt setup (mode ',cmode,').'
      call fexit()
    end if
  end if
!
end
!
!--------------------------------------------------------------------------------------
!
subroutine gennbl_for_clustering()
!
  use clusters
  use iounit
  use interfaces
!
  implicit none
!
  integer i,j,ii,jj,kk,ll,k,l,mi,mj,maxalcsz,nsets,nzeros,nclalcsz
  integer(KIND=8) testcnt,testcnt2
  RTYPE rmsdtmp,rmsdctr,bucr,rmsdclcl
  RTYPE, ALLOCATABLE:: buf1(:),buf2(:),maxrads(:)
  integer, ALLOCATABLE:: bui1(:),bui2(:)
  logical atrue
!
  atrue = .true.
!
  allocate(cnblst(cstored))
  cnblst(:)%nbs = 0
  cnblst(:)%alsz = 4
  maxalcsz = 4
  do i=1,cstored
    allocate(cnblst(i)%idx(cnblst(i)%alsz))
    allocate(cnblst(i)%dis(cnblst(i)%alsz))
  end do
!
! local screening
  bucr = cradius 
  cradius = cmaxrad
  write(ilog,*)
  write(ilog,*) 'Now using truncated leader algorithm for pre-screening in neighbor list generation ...'
  nclalcsz = 10
  allocate(scluster(nclalcsz))
  scluster(:)%nmbrs = 0
  scluster(:)%alsz = 0
  scluster(:)%nb = 0
  scluster(:)%nchildren = 0
  scluster(:)%nbalsz = 0
  scluster(:)%chalsz = 0
  k = 0
  do i=1,cstored
    ii = -1
    do j=k,max(k-500,1),-1
      call snap_to_cluster_d(rmsdtmp,scluster(j),i)
      if (rmsdtmp.lt.cradius) then
        call cluster_addsnap(scluster(j),i,rmsdtmp)
        ii = j
        exit
      end if
    end do
    if (ii.eq.-1) then
      k = k + 1
      if (k.gt.nclalcsz) call scluster_resizelst(nclalcsz,scluster)
      call cluster_addsnap(scluster(k),i,rmsdtmp)
    end if
  end do
  nsets = k
!  call leader_clustering(1,nsets)
  write(ilog,*) '... done with initial cluster generation for neighbor list.'
  write(ilog,*)
!
  write(ilog,*) 'Now finding maximum radius of all ',nsets,' identified clusters ...'
  allocate(maxrads(nsets))
  do i=1,nsets
    maxrads(i) = 0.0
    do j=1,scluster(i)%nmbrs
      call snap_to_cluster_d(rmsdtmp,scluster(i),scluster(i)%snaps(j))
      if (rmsdtmp.gt.maxrads(i)) then
        maxrads(i) = rmsdtmp
      end if
    end do
  end do
  write(ilog,*) '... done.'
  write(ilog,*)
!
! now compare all blocks to each other (the slowest part) taking advantage of information
! generated previously (otherwise intractable) 
  ii = 0
  testcnt = 0
  testcnt2 = 0
  write(ilog,*) 'Now computing cutoff-assisted neighbor list ...'
  do i=1,nsets
    ii = ii + 1
    jj = ii 
    do kk=1,scluster(i)%nmbrs
      k = scluster(i)%snaps(kk)
      do ll=kk+1,scluster(i)%nmbrs
        l = scluster(i)%snaps(ll)
        call snap_to_snap_d(rmsdtmp,k,l)
        testcnt = testcnt + 1
        if (rmsdtmp.lt.chardcut) then
!         we'll store these redundantly
          testcnt2 = testcnt2 + 1
          cnblst(k)%nbs = cnblst(k)%nbs + 1
          if (cnblst(k)%nbs.gt.cnblst(k)%alsz) call cnbl_resz(k)
          if (cnblst(k)%alsz.gt.maxalcsz) maxalcsz = cnblst(k)%alsz
          cnblst(k)%idx(cnblst(k)%nbs) = l
          cnblst(k)%dis(cnblst(k)%nbs) = rmsdtmp
          cnblst(l)%nbs = cnblst(l)%nbs + 1
          if (cnblst(l)%nbs.gt.cnblst(l)%alsz) call cnbl_resz(l)
          if (cnblst(l)%alsz.gt.maxalcsz) maxalcsz = cnblst(l)%alsz
          cnblst(l)%idx(cnblst(l)%nbs) = k
          cnblst(l)%dis(cnblst(l)%nbs) = rmsdtmp
        end if
      end do
    end do
    do j=i+1,nsets
      jj = jj + 1
      call cluster_to_cluster_d(rmsdclcl,scluster(i),scluster(j))
      testcnt = testcnt + 1
      if ((rmsdclcl - maxrads(i) - maxrads(j)).gt.chardcut) cycle
      if ((scluster(i)%nmbrs.eq.1).AND.(scluster(j)%nmbrs.eq.1)) then
        if (rmsdclcl.lt.chardcut) then
          k = scluster(i)%snaps(1)
          l = scluster(j)%snaps(1)
!         we'll store these redundantly
          testcnt2 = testcnt2 + 1
          cnblst(k)%nbs = cnblst(k)%nbs + 1
          if (cnblst(k)%nbs.gt.cnblst(k)%alsz) call cnbl_resz(k)
          if (cnblst(k)%alsz.gt.maxalcsz) maxalcsz = cnblst(k)%alsz
          cnblst(k)%idx(cnblst(k)%nbs) = l
          cnblst(k)%dis(cnblst(k)%nbs) = rmsdclcl
          cnblst(l)%nbs = cnblst(l)%nbs + 1
          if (cnblst(l)%nbs.gt.cnblst(l)%alsz) call cnbl_resz(l)
          if (cnblst(l)%alsz.gt.maxalcsz) maxalcsz = cnblst(l)%alsz
          cnblst(l)%idx(cnblst(l)%nbs) = k
          cnblst(l)%dis(cnblst(l)%nbs) = rmsdclcl
        end if        
        cycle
      end if
!      if (maxrads(j).le.maxrads(i)) then
      if (scluster(i)%nmbrs.le.scluster(j)%nmbrs) then
        mj = j
        mi = i
      else
        mj = i
        mi = j
      end if
      do kk=1,scluster(mi)%nmbrs
        k = scluster(mi)%snaps(kk)
        call snap_to_cluster_d(rmsdctr,scluster(mj),k)
        testcnt = testcnt + 1
        if ((rmsdctr - maxrads(mj)).gt.chardcut) cycle
        do ll=1,scluster(mj)%nmbrs
          l = scluster(mj)%snaps(ll)
          if (scluster(mj)%nmbrs.eq.1) then
            rmsdtmp = rmsdctr
          else
            call snap_to_snap_d(rmsdtmp,k,l)
            testcnt = testcnt + 1
          end if
          if (rmsdtmp.lt.chardcut) then
!           we'll store these redundantly
            testcnt2 = testcnt2 + 1
            cnblst(k)%nbs = cnblst(k)%nbs + 1
            if (cnblst(k)%nbs.gt.cnblst(k)%alsz) call cnbl_resz(k)
            if (cnblst(k)%alsz.gt.maxalcsz) maxalcsz = cnblst(k)%alsz
            cnblst(k)%idx(cnblst(k)%nbs) = l
            cnblst(k)%dis(cnblst(k)%nbs) = rmsdtmp
            cnblst(l)%nbs = cnblst(l)%nbs + 1
            if (cnblst(l)%nbs.gt.cnblst(l)%alsz) call cnbl_resz(l)
            if (cnblst(l)%alsz.gt.maxalcsz) maxalcsz = cnblst(l)%alsz
            cnblst(l)%idx(cnblst(l)%nbs) = k
            cnblst(l)%dis(cnblst(l)%nbs) = rmsdtmp
          end if
        end do
      end do
    end do
  end do
  write(ilog,*) '... done after computing ',(100.0*testcnt)/(0.5*cstored*(cstored-1)),'% of &
 &possible terms with ',(100.0*testcnt2)/(1.0*testcnt),'% successful.'
  write(ilog,*)
!
  write(ilog,*) 'Now sorting neighbor lists ...'
 33 format(20000(i6,1x))
 34 format(20000(f10.3,1x))
  allocate(buf1(maxalcsz))
  allocate(buf2(maxalcsz))
  allocate(bui1(maxalcsz))
  allocate(bui2(maxalcsz))
  nzeros = 0
  do i=1,cstored
    if (cnblst(i)%nbs.le.0) then
      nzeros = nzeros + 1
!      write(ilog,*) 'Warning. Snapshot # ',i,' is without a neighbor (similar) structure. This &
! &may cause the clustering algorithm to crash or misbehave otherwise.'
    else
      if (cnblst(i)%nbs.gt.1) then
        buf1(1:cnblst(i)%nbs) = cnblst(i)%dis(1:cnblst(i)%nbs)
        do ii=1,cnblst(i)%nbs
          bui1(ii) = ii
        end do
        ii = 1
        jj = cnblst(i)%nbs
        call merge_sort(ldim=cnblst(i)%nbs,up=atrue,list=buf1(1:cnblst(i)%nbs),olist=buf2(1:cnblst(i)%nbs),ilo=ii,ihi=jj,&
 &                  idxmap=bui1(1:cnblst(i)%nbs),olist2=bui2(1:cnblst(i)%nbs))
        buf1(1:cnblst(i)%nbs) = cnblst(i)%dis(1:cnblst(i)%nbs)
        bui1(1:cnblst(i)%nbs) = cnblst(i)%idx(1:cnblst(i)%nbs)
        do j=1,cnblst(i)%nbs
          cnblst(i)%dis(j) = buf1(bui2(j))
          cnblst(i)%idx(j) = bui1(bui2(j))
        end do
      end if
      allocate(cnblst(i)%tagged(cnblst(i)%nbs))
      cnblst(i)%tagged(:) = .false.
    end if
  end do
  if (nzeros.gt.0) then
    write(ilog,*) 'Warning. ',nzeros,' snapshots are without a neighbor (similar) structure. This &
 &may in some cases cause the clustering algorithm to misbehave.'
  end if
  write(ilog,*) '... done.'
  write(ilog,*) 
!
#ifdef LINK_NETCDF
  write(ilog,*) 'Dumping to binary NetCDF-file ...'
    call dump_nbl_nc()
  write(ilog,*) '... done.'
  write(ilog,*)
#endif
!
  cradius = bucr
!
  deallocate(bui2)
  deallocate(bui1)
  deallocate(buf2)
  deallocate(buf1)
  deallocate(maxrads)
  do i=1,nsets
    if (allocated(scluster(i)%snaps).EQV..true.) deallocate(scluster(i)%snaps)
    if (allocated(scluster(i)%tmpsnaps).EQV..true.) deallocate(scluster(i)%tmpsnaps)
    if (allocated(scluster(i)%sums).EQV..true.) deallocate(scluster(i)%sums)
    if (allocated(scluster(i)%map).EQV..true.) deallocate(scluster(i)%map)
    if (allocated(scluster(i)%children).EQV..true.) deallocate(scluster(i)%children)
    if (allocated(scluster(i)%wghtsnb).EQV..true.) deallocate(scluster(i)%wghtsnb)
    if (allocated(scluster(i)%lstnb).EQV..true.) deallocate(scluster(i)%lstnb)
    if (allocated(scluster(i)%flwnb).EQV..true.) deallocate(scluster(i)%flwnb)
  end do
  deallocate(scluster)
!
end
!
!------------------------------------------------------------------------------
!
#ifdef LINK_NETCDF
!
subroutine dump_nbl_nc()
!
  use iounit
  use system
  use mcsums
  use clusters
  use netcdf
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer i,ncid,ii,jj,freeunit,xlen,istart
  integer, ALLOCATABLE:: helper(:)
!  integer nf90_redef,nf90_put_att,nf90_def_dim
  logical exists
  character(MAXSTRLEN) attstring,dumpfile
  real(KIND=4), ALLOCATABLE:: prthlp(:)
#ifdef ENABLE_MPI
  integer tl
  character(3) nod
#endif
!
#ifdef ENABLE_MPI
  tl = 3
  call int2str(myrank,nod,tl)
  if (use_REMC.EQV..true.) then
    dumpfile = 'N_'//nod(1:tl)//'_FRAMES_NBL.nc'
  else
    dumpfile = 'FRAMES_NBL.nc'
  end if
#else
  dumpfile = 'FRAMES_NBL.nc'
#endif
  call strlims(dumpfile,ii,jj)
  inquire(file=dumpfile(ii:jj),exist=exists)
  if (exists.EQV..true.) then
    ncid = freeunit()
    open(unit=ncid,file=dumpfile(ii:jj),status='old')
    close(unit=ncid,status='delete')
  end if
  ncid = freeunit()
  call check_fileionetcdf( nf90_create(path=dumpfile(ii:jj), cmode=IOR(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )
!
! enable definition
  do i=1,MAXSTRLEN
    attstring(i:i) = " "
  end do
  attstring(1:7) = "CAMPARI"
  call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "program", attstring(1:7)) )
  attstring(1:7) = "       "
  attstring(1:3) = "XXX"
  call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "programVersion", attstring(1:3)) )
  attstring(1:3) = "   "
  call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "title", basename(1:bleng)) )
! define dimensions
  attstring(1:10) = "framepairs"
  call check_fileionetcdf( nf90_def_dim(ncid, attstring(1:10), sum(cnblst(1:cstored)%nbs), cnc_ids(1)) )
  attstring(1:10) = "     " 
! define (not set) variables to hold distance and type of distance information
  attstring(1:9) = "snapshots"
  call check_fileionetcdf( nf90_def_var(ncid, attstring(1:9), NF90_INT, cnc_ids(1), cnc_ids(2)) )
  attstring(1:9) = "neighbors"
  call check_fileionetcdf( nf90_def_var(ncid, attstring(1:9), NF90_INT, cnc_ids(1), cnc_ids(3)) )
  attstring(1:9) = "distances"
  call check_fileionetcdf( nf90_def_var(ncid, attstring(1:9), NF90_FLOAT, cnc_ids(1), cnc_ids(4)) )
  attstring(1:9) = "         "
  if (cdis_crit.eq.1) then
    xlen = 14
    attstring(1:xlen) = "torsional RMSD"
  else if (cdis_crit.eq.2) then
    xlen = 37
    attstring(1:xlen) = "inertial-mass-weighted torsional RMSD"
  else if (cdis_crit.eq.3) then
    xlen = 31
    attstring(1:xlen) = "RMSD of torsional Fourier terms"
  else if (cdis_crit.eq.4) then
    xlen = 54
    attstring(1:xlen) = "inertial-mass-weighted RMSD of torsional Fourier terms"
  else if (cdis_crit.eq.5) then
    if (align_for_clustering.EQV..true.) then
      xlen = 19
      attstring(1:xlen) = "aligned atomic RMSD"
    else
      xlen = 21
      attstring(1:xlen) = "unaligned atomic RMSD"
    end if
  else if (cdis_crit.eq.6) then
    if (align_for_clustering.EQV..true.) then
      xlen = 42
      attstring(1:xlen) = "aligned atomic RMSD (may be separate sets)"
    else
      xlen = 21
      attstring(1:xlen) = "unaligned atomic RMSD"
    end if
  else if (cdis_crit.eq.7) then
    xlen = 29
    attstring(1:xlen) = "internal distance vector RMSD"
  else if (cdis_crit.eq.8) then
    xlen = 43
    attstring(1:xlen) = "mass-weighted internal distance vector RMSD"
  else if (cdis_crit.eq.9) then
    xlen = 46
    attstring(1:xlen) = "locally weighted internal distance vector RMSD" 
  else if (cdis_crit.eq.10) then
    xlen = 43
    attstring(1:xlen) = "locally weighted, unaligned coordinate RMSD"  
  else
    write(ilog,*) 'Fatal. Unsupported distance criterion in dump_nbl_nc(...). This is an omission bug.'
    call fexit()
  end if
  call check_fileionetcdf( nf90_put_att(ncid, cnc_ids(1), "type", attstring(1:xlen)) )
!  call check_fileionetcdf( nf90_def_var(ncid, attstring(1:xlen), NF90_FLOAT, cnc_ids(2), cnc_ids(5)) )
  attstring(1:5) = "units"
  if (cdis_crit.le.2) then
    attstring(6:12) = "degrees"
    call check_fileionetcdf( nf90_put_att(ncid, cnc_ids(1), attstring(1:5), attstring(6:12)) )
  else if ((cdis_crit.ge.3).AND.(cdis_crit.le.4)) then
    attstring(6:13) = "unitless"
    call check_fileionetcdf( nf90_put_att(ncid, cnc_ids(1), attstring(1:5), attstring(6:13)) )
  else
    attstring(6:13) = "angstrom"
    call check_fileionetcdf( nf90_put_att(ncid, cnc_ids(1), attstring(1:5), attstring(6:13)) )
  end if
  attstring(1:13) = "               "
! quit define mode
  call check_fileionetcdf( nf90_enddef(ncid) )
  call check_fileionetcdf( nf90_sync(ncid) )
!
! put the data
  allocate(helper(cstored))
  allocate(prthlp(cstored))
  istart = 1
  do i=1,cstored
    if (cnblst(i)%nbs.le.0) cycle
    helper(1:cnblst(i)%nbs) = i
    prthlp(1:cnblst(i)%nbs) = REAL(cnblst(i)%dis(1:cnblst(i)%nbs),KIND=4)
    call check_fileionetcdf( nf90_put_var(ncid, cnc_ids(2), helper(1:cnblst(i)%nbs), &
 &                                       start = (/ istart /), count = (/ cnblst(i)%nbs /)) )
    call check_fileionetcdf( nf90_put_var(ncid, cnc_ids(3), cnblst(i)%idx(1:cnblst(i)%nbs), &
 &                                       start = (/ istart /), count = (/ cnblst(i)%nbs /)) )
    call check_fileionetcdf( nf90_put_var(ncid, cnc_ids(4), REAL(prthlp(1:cnblst(i)%nbs),KIND=8), &
 &                                       start = (/ istart /), count = (/ cnblst(i)%nbs /)) )
    istart = istart + cnblst(i)%nbs
  end do
!
  deallocate(prthlp)
  deallocate(helper)
!
! close
  call check_fileionetcdf( nf90_close(ncid) )
! 
end
!
!---------------------------------------------
!
subroutine read_nbl_nc()
!
  use netcdf
  use clusters
  use iounit
!
  implicit none
!
  integer ncid,t1,t2,fndds,dimlen,nframes,i,ret,ilast,istart
  logical exists
  character(MAXSTRLEN) ucstr,trystr
  integer, ALLOCATABLE:: vnbs(:),vsnp(:)
  real(KIND=4), ALLOCATABLE:: vdis(:)
!
  call strlims(nblfilen,t1,t2)
  inquire (file=nblfilen(t1:t2),exist=exists)
! 
  if (exists.EQV..false.) then
    write(ilog,*) 'Fatal. Cannot open spec.d input file for reading &
 &from NetCDF (',nblfilen(t1:t2),') in setup_netcdftraj().'
    call fexit()
  end if
!
! open
 44 format('Warning. Ambiguous dimensions in NetCDF file. Encountered ',a,' twice but keeping&
 &only first.')
  call check_fileionetcdf( nf90_open(nblfilen(t1:t2), NF90_NOWRITE, ncid) )
!
! find the necessary dimensions: three are required
  fndds = 0
  nframes = 0
  do i=1,NF90_MAX_DIMS
    ret = nf90_inquire_dimension(ncid,i,trystr,dimlen)
    if (ret.eq.NF90_NOERR) then
      ucstr(1:15) = trystr(1:15)
      call toupper(ucstr(1:15))
      if (ucstr(1:10).eq.'FRAMEPAIRS') then
        if (fndds.eq.1) then 
          write(ilog,44) ucstr(1:10)
        else
          nframes = dimlen
          fndds = 1
          cnc_ids(1) = i
        end if
      end if
    else if (ret.eq.NF90_EBADDIM) then
!     do nothing
    else ! get us out of here
      call check_fileionetcdf( nf90_inquire_dimension(ncid,i,trystr,dimlen) )
    end if
  end do
!
  if (nframes.lt.1) then
    write(ilog,*) 'Fatal. NetCDF-file (',nblfilen(t1:t2),') has no neighbor &
 &data (empty containers).'
    call fexit()
  end if
!
! now find the necessary variables, only coordinates are required
  ret = nf90_inq_varid(ncid,"snapshots",cnc_ids(2))
  if (ret.eq.NF90_NOERR) then
!   do nothing
  else if (ret.eq.NF90_ENOTVAR) then
    write(ilog,*) 'Fatal. Variable "snapshots" not found in NetCDF-file (',&
 &nblfilen(t1:t2),'). Use NetCDFs ncdump utility to check file.'
    call fexit()
  else
    call check_fileionetcdf( nf90_inq_varid(ncid,"snapshots",cnc_ids(2)) )
  end if
  ret = nf90_inq_varid(ncid,"neighbors",cnc_ids(3))
  if (ret.eq.NF90_NOERR) then
!   do nothing
  else if (ret.eq.NF90_ENOTVAR) then
    write(ilog,*) 'Fatal. Variable "neighbors" not found in NetCDF-file (',&
 &nblfilen(t1:t2),'). Use NetCDFs ncdump utility to check file.'
    call fexit()
  else
    call check_fileionetcdf( nf90_inq_varid(ncid,"neighbors",cnc_ids(3)) )
  end if
  ret = nf90_inq_varid(ncid,"distances",cnc_ids(4))
  if (ret.eq.NF90_NOERR) then
!   do nothing
  else if (ret.eq.NF90_ENOTVAR) then
    write(ilog,*) 'Fatal. Variable "distances" not found in NetCDF-file (',&
 &nblfilen(t1:t2),'). Use NetCDFs ncdump utility to check file.'
    call fexit()
  else
    call check_fileionetcdf( nf90_inq_varid(ncid,"distances",cnc_ids(4)) )
  end if
!
  allocate(vnbs(nframes))
  allocate(vdis(nframes))
  allocate(vsnp(nframes))
!
! for some reason nf90_get_var chokes if the increment (count) is very large
  istart = 1
  ilast = min(nframes,10000)
  do while (istart.le.nframes)
    call check_fileionetcdf( nf90_get_var(ncid, cnc_ids(2), vsnp(istart:ilast), &
 &                                        start = (/ istart /) , count = (/ ilast-istart+1 /)) )
    call check_fileionetcdf( nf90_get_var(ncid, cnc_ids(3), vnbs(istart:ilast),&
 &                                        start = (/ istart /) , count = (/ ilast-istart+1 /)) )
    call check_fileionetcdf( nf90_get_var(ncid, cnc_ids(4), vdis(istart:ilast),&
 &                                        start = (/ istart /) , count = (/ ilast-istart+1 /)) )
    istart = ilast + 1
    ilast = min(nframes,istart+10000)
  end do
! 
  allocate(cnblst(cstored))
  cnblst(:)%nbs = 0
  i = 1
  ilast = 0
  do istart=vsnp(1),vsnp(nframes)
    if (istart.ne.vsnp(i)) cycle
    do while (vsnp(i).eq.istart)
      if (i.eq.nframes) exit
      if (vsnp(i+1).ne.vsnp(i)) exit
      i = i + 1
    end do
    cnblst(vsnp(i))%nbs = i - ilast
    cnblst(vsnp(i))%alsz = cnblst(vsnp(i))%nbs
    if (cnblst(vsnp(i))%nbs.gt.0) then
      allocate(cnblst(vsnp(i))%idx(cnblst(vsnp(i))%nbs))
      allocate(cnblst(vsnp(i))%dis(cnblst(vsnp(i))%nbs))
      allocate(cnblst(vsnp(i))%tagged(cnblst(vsnp(i))%nbs))
      cnblst(vsnp(i))%idx(:) = vnbs(ilast+1:i)
      cnblst(vsnp(i))%dis(:) = REAL(vdis(ilast+1:i),KIND=4)
      cnblst(vsnp(i))%tagged(:) = .false.
    else
      cnblst(vsnp(i))%nbs = 0
      cnblst(vsnp(i))%alsz = 2
      allocate(cnblst(vsnp(i))%idx(cnblst(vsnp(i))%alsz))
      allocate(cnblst(vsnp(i))%dis(cnblst(vsnp(i))%alsz))
      allocate(cnblst(vsnp(i))%tagged(cnblst(vsnp(i))%alsz))
    end if
    if (i.eq.nframes) exit
    ilast = i
    i = i  +1 
  end do
  deallocate(vsnp)
  deallocate(vnbs)
  deallocate(vdis)
!
end
!
#endif
!
!----------------------------------------------------------------------------
!
! set modei to 0 to obtain a silent version without post-processing
!
subroutine birch_clustering(modei,nnodes)
!
  use clusters
  use iounit
  use interfaces
!
  implicit none
!
  integer i,j,k,ii,jj,kk,ll,mm,atwo,fail,kkf,thekk,nlst1,nlst2,snapstart,snapend,snapinc,nnodes,modei
  integer, ALLOCATABLE:: kklst(:,:),errcnt(:),kkhistory(:)
  integer(KIND=8) cnt1,cnt2
  RTYPE, ALLOCATABLE:: scrcts(:)
  RTYPE rdv,mind,helper,maxd(2),normer(4),qualmet(4)
  logical atrue,afalse,notdone
!
 33 format('ERROR B: ',i6,' could be part of ',i5,' at ',g14.7,' (last d',g14.6,')')
 34 format('NEXT HIGHER: ',i6,' could have been part of ',i5,' at ',g14.7,'.')
 35 format('ERROR C: ',i6,' is meant to a child of ',i5,' at level ',i5,' but has distance ',g14.7,'.')
!
  atrue = .true.
  afalse = .false.
  atwo = 2
  if (cmaxrad.le.cradius) cmaxrad = 2.0*cradius
!
  if (cleadermode.le.2) then
    snapstart = 1
    snapend = cstored
    snapinc = 1
  else
    snapstart = cstored
    snapend = 1
    snapinc = -1
  end if
!
  allocate(kklst(cstored,2))
  allocate(kkhistory(c_nhier+1))
  allocate(errcnt(c_nhier+1))
  allocate(birchtree(c_nhier+1))
  do ii=1,c_nhier+1
    allocate(birchtree(ii)%cls(10))
    birchtree(ii)%ncls = 0
    birchtree(ii)%nclsalsz = 10
  end do
  allocate(scrcts(c_nhier+1))
  errcnt(:) = 0
  birchtree(1)%ncls = 1
  scrcts(c_nhier+1) = cradius
  if (c_nhier.gt.1) then
    scrcts(2) = cmaxrad
    do i=3,c_nhier
      scrcts(i) = cmaxrad - ((i-2.0)/(c_nhier-1.0))*(cmaxrad - cradius) ! linear so far
    end do
  end if
  do i=1,c_nhier+1
    birchtree(i)%cls(:)%nmbrs = 0
    birchtree(i)%cls(:)%nb = 0
    birchtree(i)%cls(:)%nchildren = 0
    birchtree(i)%cls(:)%alsz = 0
    birchtree(i)%cls(:)%nbalsz = 0
    birchtree(i)%cls(:)%chalsz = 0
    birchtree(i)%cls(:)%parent = 0
  end do

  allocate(birchtree(1)%cls(1)%snaps(2))
  birchtree(1)%cls(1)%alsz = 2
  allocate(birchtree(1)%cls(1)%children(2))
  birchtree(1)%cls(1)%chalsz = 2
!
  if (modei.gt.0) then
    write(ilog,*)
    write(ilog,*) 'Now performing tree-based clustering ...'
  end if
  cnt1 = 0
  cnt2 = 0
  do i=snapstart,snapend,snapinc
    kk = 1
    notdone = .true.
    fail = -1
    nlst1 = 1
    kklst(1,1) = kk
!    cnt1 = 0
    do ii=2,c_nhier
      jj = -1
      mind = HUGE(mind)
      nlst2 = 0
      do mm=1,nlst1
        kk = kklst(mm,1)
        do j=1,birchtree(ii-1)%cls(kk)%nchildren
          ll = birchtree(ii-1)%cls(kk)%children(j)
          cnt1 = cnt1 + 1
          if ((birchtree(ii)%cls(ll)%center.le.0).OR.(birchtree(ii)%cls(ll)%center.gt.cstored)) call fexit()
          call snap_to_cluster_d(rdv,birchtree(ii)%cls(ll),i)
          if (rdv.lt.mind) then
            mind = rdv
            jj = j
            thekk = kk
          end if
          if (rdv.lt.scrcts(ii)) then
            nlst2 = nlst2 + 1
            kklst(nlst2,2) = ll
          end if
        end do
      end do
!     store the path
      if (jj.eq.-1) then
        kkhistory(ii-1) = 1
      else
        kkhistory(ii-1) = thekk
      end if
      if ((ii.le.(c_nhier-1)).AND.(jj.gt.0)) then
        if (nlst2.le.0) then ! absolutely nothing nearby 
          nlst1 = 1
          kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj)
          if (fail.eq.-1) then
!            write(*,*) i,' failed at ',ii,' w/ ',mind
            fail = ii
            kkf = thekk
          end if
        else
          if (fail.gt.0) then
            do j=fail,ii-1
              birchtree(j)%ncls = birchtree(j)%ncls + 1
              if (birchtree(j)%ncls.gt.birchtree(j)%nclsalsz) call scluster_resizelst(birchtree(j)%nclsalsz,birchtree(j)%cls)
              kkhistory(j) = birchtree(j)%ncls
              if (j.gt.fail) then
                call cluster_addchild(birchtree(j-1)%cls(birchtree(j-1)%ncls),birchtree(j-1)%ncls,&
 &                                    birchtree(j)%cls(birchtree(j)%ncls),birchtree(j)%ncls)
                cnt2 = cnt2 + 1
              end if
            end do
            call cluster_addchild(birchtree(fail-1)%cls(kkf),kkf,birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)
            call cluster_addchild(birchtree(ii-1)%cls(birchtree(ii-1)%ncls),&
 &                                birchtree(ii)%cls(jj)%parent,birchtree(ii)%cls(jj),jj)
            fail = -1
          end if
          nlst1 = 1
          kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj)
        end if
        cycle
!     leaf
      else if ((ii.eq.c_nhier).AND.(mind.lt.scrcts(ii))) then
        kk = birchtree(ii-1)%cls(thekk)%children(jj)
        call cluster_addsnap(birchtree(ii)%cls(kk),i,rdv)
        do j=2,c_nhier
          call cluster_addsnap(birchtree(j-1)%cls(kkhistory(j-1)),i,rdv)
        end do
        notdone = .false.
      else
        if (fail.eq.-1) then
          fail = ii
          kkf = kk
        end if
        kk = kkf
        do j=fail,c_nhier
          birchtree(j)%ncls = birchtree(j)%ncls + 1
          if (birchtree(j)%ncls.gt.birchtree(j)%nclsalsz) call scluster_resizelst(birchtree(j)%nclsalsz,birchtree(j)%cls)
          call cluster_addsnap(birchtree(j)%cls(birchtree(j)%ncls),i,rdv)
          if (j.gt.fail) then
            call cluster_addchild(birchtree(j-1)%cls(birchtree(j-1)%ncls),birchtree(j-1)%ncls,&
 &                                birchtree(j)%cls(birchtree(j)%ncls),birchtree(j)%ncls)
            cnt2 = cnt2 + 1
          end if
        end do
        call cluster_addchild(birchtree(fail-1)%cls(kk),kk,birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)
        cnt2 = cnt2 + 1
        do j=2,fail
          call cluster_addsnap(birchtree(j-1)%cls(kkhistory(j-1)),i,rdv)
        end do
        call snap_to_cluster_d(maxd(1),birchtree(fail-1)%cls(kk),i)
        notdone = .false.
      end if
      if (notdone.EQV..false.) exit
    end do
  end do

  do i=snapstart,snapend,snapinc
    kk = 1
    notdone = .true.
    fail = -1
    nlst1 = 1
    kklst(1,1) = kk
    do ii=2,c_nhier+1
      jj = -1
      mind = HUGE(mind)
      nlst2 = 0
      do mm=1,nlst1
        kk = kklst(mm,1)
        do j=1,birchtree(ii-1)%cls(kk)%nchildren
          ll = birchtree(ii-1)%cls(kk)%children(j)
          cnt1 = cnt1 + 1
          if ((birchtree(ii)%cls(ll)%center.le.0).OR.(birchtree(ii)%cls(ll)%center.gt.cstored)) call fexit()
          call snap_to_cluster_d(rdv,birchtree(ii)%cls(ll),i)
          if (rdv.lt.mind) then
            mind = rdv
            jj = j
            thekk = kk
          end if
          if (rdv.lt.scrcts(ii)) then
            nlst2 = nlst2 + 1
            kklst(nlst2,2) = ll
          end if
        end do
      end do
!     store the path
      if (jj.eq.-1) then
        kkhistory(ii-1) = 1
      else
        kkhistory(ii-1) = thekk
      end if
      if ((ii.le.c_nhier).AND.(jj.gt.0)) then
        if (nlst2.le.0) then ! absolutely nothing nearby 
          nlst1 = 1
          kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj)
          if (fail.eq.-1) then
!            write(*,*) i,' failed at ',ii,' w/ ',mind
            fail = ii
            kkf = thekk
          end if
        else 
          nlst1 = 1
          kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj)
        end if
        cycle
!     leaf
      else if ((ii.eq.(c_nhier+1)).AND.(mind.lt.scrcts(ii))) then
        kk = birchtree(ii-1)%cls(thekk)%children(jj)
        call cluster_addsnap(birchtree(ii)%cls(kk),i,rdv)
        notdone = .false.
      else if (fail.gt.0) then
        kk = kkf
        do j=fail,c_nhier+1
          birchtree(j)%ncls = birchtree(j)%ncls + 1
          if (birchtree(j)%ncls.gt.birchtree(j)%nclsalsz) call scluster_resizelst(birchtree(j)%nclsalsz,birchtree(j)%cls)
          call cluster_addsnap(birchtree(j)%cls(birchtree(j)%ncls),i,rdv)
          if (j.gt.fail) then
            call cluster_addchild(birchtree(j-1)%cls(birchtree(j-1)%ncls),birchtree(j-1)%ncls,&
 &                                birchtree(j)%cls(birchtree(j)%ncls),birchtree(j)%ncls)
                     cnt2 = cnt2 + 1
          end if
        end do
        call cluster_addchild(birchtree(fail-1)%cls(kk),kk,birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)
        cnt2 = cnt2 + 1
        notdone = .false.
      else
        fail = ii
        j = c_nhier+1
        birchtree(j)%ncls = birchtree(j)%ncls + 1
        if (birchtree(j)%ncls.gt.birchtree(j)%nclsalsz) call scluster_resizelst(birchtree(j)%nclsalsz,birchtree(j)%cls)
        call cluster_addsnap(birchtree(j)%cls(birchtree(j)%ncls),i,rdv)
        call cluster_addchild(birchtree(fail-1)%cls(kk),kk,birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)
        notdone = .false.
      end if
      if (notdone.EQV..false.) exit
    end do
  end do

 66 format('Level    # Clusters     Threshold     Total Snaps    Total Children')
 67 format(i9,i10,1x,g14.4,4x,i12,4x,i12)
 68 format(i9,i10,5x,a7,7x,i12,4x,i12)
  if (modei.gt.0) then
    write(ilog,66)
    write(ilog,68) c_nhier+1,birchtree(1)%ncls,'MAXIMAL',sum(birchtree(1)%cls(1:birchtree(1)%ncls)%nmbrs),&
 &               sum(birchtree(1)%cls(1:birchtree(1)%ncls)%nchildren)
    do i=2,c_nhier+1
      write(ilog,67) c_nhier+2-i,birchtree(i)%ncls,scrcts(i),sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nmbrs),&
 &               sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nchildren)
    end do
    write(ilog,*) '---------------------------------------------------------------------'
    write(ilog,*)
    write(ilog,*) '... done after a total of ',cnt1,' distance evaluations.'
    write(ilog,*)
  end if
!
  do j=1,birchtree(c_nhier+1)%ncls
    call cluster_calc_params(birchtree(c_nhier+1)%cls(j),scrcts(c_nhier+1))
  end do
!
  if (refine_clustering.EQV..true.) then
    call quality_of_clustering(birchtree(c_nhier+1)%ncls,birchtree(c_nhier+1)%cls,scrcts(c_nhier+1),qualmet)
    if (modei.gt.0) then
      write(ilog,*) 'Now merging clusters that yield joint reduced average intracluster distance ...'
    end if
 77 format('Would join ',i5,' (',i6,') and ',i5,'(',i6,') from: ',/,'Diam: ',g10.4,' / ',g10.4,'; Rad.: ',&
 &g10.4,' / ',g10.4,' to ',g10.4,' / ',g10.4,'.')
    cnt1 = 0
    cnt2 = 0
    do i=1,birchtree(c_nhier)%ncls
      nlst1 = 0
      do kkf=1,birchtree(c_nhier)%ncls
        if (i.eq.kkf) cycle
        call cluster_to_cluster_d(rdv,birchtree(c_nhier)%cls(i),birchtree(c_nhier)%cls(kkf))
        cnt1 = cnt1 + 1
        if (rdv.lt.scrcts(c_nhier)) then
          nlst1 = nlst1 + 1
          kklst(nlst1,1) = kkf
        end if
      end do
      do j=1,birchtree(c_nhier)%cls(i)%nchildren
        jj = birchtree(c_nhier)%cls(i)%children(j)
        if (birchtree(c_nhier+1)%cls(jj)%nmbrs.le.0) cycle
        do kkf=1,nlst1
          do k=1,birchtree(c_nhier)%cls(kklst(kkf,1))%nchildren
            kk = birchtree(c_nhier)%cls(kklst(kkf,1))%children(k)
            if (birchtree(c_nhier+1)%cls(kk)%nmbrs.le.0) cycle
            if (jj.eq.kk) cycle
            call clusters_joint_diam(birchtree(c_nhier+1)%cls(jj),birchtree(c_nhier+1)%cls(kk),rdv,helper)
            cnt1 = cnt1 + 1
            normer(1) = 0.5*birchtree(c_nhier+1)%cls(jj)%nmbrs*(birchtree(c_nhier+1)%cls(jj)%nmbrs-1.0)
            normer(2) = 0.5*birchtree(c_nhier+1)%cls(kk)%nmbrs*(birchtree(c_nhier+1)%cls(kk)%nmbrs-1.0)
            normer(3) = 1.0/(1.0*(birchtree(c_nhier+1)%cls(jj)%nmbrs + birchtree(c_nhier+1)%cls(kk)%nmbrs))
            normer(4) = 0.0
            if (sum(normer(1:2)).gt.0) normer(4) = 1.0/sum(normer(1:2))
            maxd(1) =  normer(4)*(normer(1)*birchtree(c_nhier+1)%cls(jj)%diam + &
              &                   normer(2)*birchtree(c_nhier+1)%cls(kk)%diam)
            maxd(2) =  normer(3)*(birchtree(c_nhier+1)%cls(jj)%nmbrs*birchtree(c_nhier+1)%cls(jj)%radius + &
   &                              birchtree(c_nhier+1)%cls(kk)%nmbrs*birchtree(c_nhier+1)%cls(kk)%radius)
            if ((helper.le.maxd(2)).OR.(rdv.le.maxd(1))) then
              if (birchtree(c_nhier+1)%cls(jj)%nmbrs.gt.birchtree(c_nhier+1)%cls(kk)%nmbrs) then
                call join_clusters(birchtree(c_nhier+1)%cls(jj),birchtree(c_nhier+1)%cls(kk))
                birchtree(c_nhier+1)%cls(jj)%diam = rdv
                birchtree(c_nhier+1)%cls(jj)%radius = helper
                cnt2 = cnt2 + 1
              else
                call join_clusters(birchtree(c_nhier+1)%cls(kk),birchtree(c_nhier+1)%cls(jj))
                birchtree(c_nhier+1)%cls(kk)%diam = rdv
                birchtree(c_nhier+1)%cls(kk)%radius = helper
                cnt2 = cnt2 + 1
                exit
              end if
            end if
          end do
          if (birchtree(c_nhier+1)%cls(jj)%nmbrs.eq.0) exit
        end do
      end do
    end do
    if (modei.gt.0) then
      write(ilog,*) '... done after a total of ',cnt2,' merges requiring ',cnt1,' additional &
 &distance or joint size evaluations.'
      write(ilog,*)
    end if
  end if
!
! now shorten list and resort
  call clusters_shorten(birchtree(c_nhier+1)%cls,birchtree(c_nhier+1)%ncls)
  do j=1,birchtree(c_nhier+1)%ncls
    call cluster_calc_params(birchtree(c_nhier+1)%cls(j),scrcts(c_nhier+1))
  end do
  call clusters_sort(birchtree(c_nhier+1)%cls,birchtree(c_nhier+1)%ncls,afalse)
  do i=1,birchtree(c_nhier+1)%ncls
    call cluster_getcenter(birchtree(c_nhier+1)%cls(i))
  end do
! lastly, copy into global cluster array
  allocate(scluster(birchtree(c_nhier+1)%ncls))
  do i=1,birchtree(c_nhier+1)%ncls
    call copy_cluster(birchtree(c_nhier+1)%cls(i),scluster(i))
  end do
  nnodes = birchtree(c_nhier+1)%ncls
  do k=2,c_nhier+1
    do i=1,birchtree(k)%ncls
      call cluster_calc_params(birchtree(k)%cls(i),scrcts(k))
    end do
  end do
!
 63 format(i7,1x,i7,1x,i8,1000(1x,g12.5))
 64 format(1000(g12.5,1x))
!
  if (modei.gt.0) then
    write(ilog,*) '------------- CLUSTER SUMMARY ------------------'
    write(ilog,*) ' #       No.     "Center"  Diameter     Radius      '
    do i=1,birchtree(c_nhier+1)%ncls
      write(ilog,63) i,scluster(i)%nmbrs,scluster(i)%center,scluster(i)%diam,scluster(i)%radius
    end do
    write(ilog,*) '------------------------------------------------'
    write(ilog,*)
!
    atrue = .true.
    call quality_of_clustering(nnodes,scluster,scrcts(c_nhier+1),qualmet)
    call gen_graph_from_clusters(scluster,nnodes,atrue)
    call graphml_helper_for_clustering(scluster,nnodes)
    call vmd_helper_for_clustering(scluster,nnodes)
  end if
!
  deallocate(kklst)
  deallocate(kkhistory)
  deallocate(scrcts)
  deallocate(errcnt)  
!
end
!
!----------------------------------------------------------------------------
!
subroutine leader_clustering(mode,nnodes)
!
  use clusters
  use iounit
  use interfaces
!
  implicit none
!
  integer i,j,ii,jj,k,kk,l,nsets,mode,atwo,nsets_old,globi,nclalcsz
  integer snapstart,snapend,snapinc,nnodes
  integer(KIND=8) cnt1,cnt2
  integer, ALLOCATABLE:: iv1(:)
  RTYPE rdv,coreval,qualmet(4)
  logical atrue,afalse
!
  atrue = .true.
  afalse = .false.
  atwo = 2
  cnt1 = 0
  if (cleadermode.le.2) then
    snapstart = 1
    snapend = cstored
    snapinc = 1
  else
    snapstart = cstored
    snapend = 1
    snapinc = -1
  end if
  write(ilog,*)
!
  if (mode.eq.1) then
    write(ilog,*) 'Now performing LEADER-clustering ...'
    nclalcsz = 10
    allocate(scluster(nclalcsz))
    scluster(:)%nmbrs = 0
    scluster(:)%alsz = 0
    scluster(:)%nb = 0
    scluster(:)%nchildren = 0
    scluster(:)%nbalsz = 0
    scluster(:)%chalsz = 0
    scluster(:)%parent = 0
    k = 0 
    if ((cleadermode.eq.1).OR.(cleadermode.eq.3)) then
      do i=snapstart,snapend,snapinc
        ii = -1
        do j=k,1,-1
          cnt1 = cnt1 + 1
          call snap_to_snap_d(rdv,i,scluster(j)%center)
          if (rdv.lt.cradius) then
            call cluster_addsnap(scluster(j),i,rdv)
            ii = j
            exit
          end if
        end do
        if (ii.eq.-1) then
          k = k + 1
          if (k.gt.nclalcsz) call scluster_resizelst(nclalcsz,scluster)
          call cluster_addsnap(scluster(k),i,rdv)
        end if
      end do
   else
      do i=snapstart,snapend,snapinc
        ii = -1
        do j=1,k
          cnt1 = cnt1 + 1
          call snap_to_snap_d(rdv,i,scluster(j)%center)
          if (rdv.lt.cradius) then
            call cluster_addsnap(scluster(j),i,rdv)
            ii = j
            exit
          end if
        end do
        if (ii.eq.-1) then
          k = k + 1
          if (k.gt.nclalcsz) call scluster_resizelst(nclalcsz,scluster)
          call cluster_addsnap(scluster(k),i,rdv)
        end if
      end do
    end if
    write(ilog,*) 'Done after evaluating ',cnt1,' pairwise distances.'
    write(ilog,*)
!
    nsets = k
!
!   now sort
    call clusters_sort(scluster,nsets,afalse)
    do i=1,nsets
      call cluster_sortsnaps(scluster(i))
      call cluster_calc_params(scluster(i),cradius)
    end do
!
  else if (mode.eq.2) then
!
    write(ilog,*) 'Now performing first stage of modified LEADER-clustering ...'
    nclalcsz = 10
    allocate(scluster(nclalcsz))
    scluster(:)%nmbrs = 0
    scluster(:)%alsz = 0
    scluster(:)%nb = 0
    scluster(:)%nchildren = 0
    scluster(:)%nbalsz = 0
    scluster(:)%chalsz = 0
    scluster(:)%parent = 0
    k = 0
    if ((cleadermode.eq.1).OR.(cleadermode.eq.3)) then
      do i=snapstart,snapend,snapinc
        ii = -1
        do j=k,1,-1 !max(k-100,1),-1
          cnt1 = cnt1 + 1
          call snap_to_cluster_d(rdv,scluster(j),i)
          if (rdv.lt.cradius) then
            call cluster_addsnap(scluster(j),i,rdv)
            ii = j
            exit
          end if
        end do
        if (ii.eq.-1) then
          k = k + 1
          if (k.gt.nclalcsz) call scluster_resizelst(nclalcsz,scluster)
          call cluster_addsnap(scluster(k),i,rdv)
        end if
      end do
    else
      do i=snapstart,snapend,snapinc
        ii = -1
        do j=1,k
          cnt1 = cnt1 + 1
          call snap_to_cluster_d(rdv,scluster(j),i)
          if (rdv.lt.cradius) then
            call cluster_addsnap(scluster(j),i,rdv)
            ii = j
            exit
          end if
        end do
        if (ii.eq.-1) then
          k = k + 1
          if (k.gt.nclalcsz) call scluster_resizelst(nclalcsz,scluster)
          call cluster_addsnap(scluster(k),i,rdv)
        end if
      end do
    end if
    write(ilog,*) 'Done after evaluating ',cnt1,' pairwise distances.'
    write(ilog,*)
!
    nsets = k
    nsets_old = k
!
    call clusters_sort(scluster,nsets,afalse)
    do i=1,nsets
      call cluster_sortsnaps(scluster(i))
      call cluster_calc_params(scluster(i),cradius)
    end do
!
    cnt1 = 0
    cnt2 = 0
!
    if (refine_clustering.EQV..true.) then
      allocate(iv1(cstored))
!     the first block is to remove overlap
      write(ilog,*) 'Now removing cluster overlap by reassigning ambiguous frames to larger clusters ...'
!     we're doing this exactly twice - iteration would not be guaranteed to converge due to moving centers
      do globi=1,2
        do i=nsets,1,-1
          if (allocated(scluster(i)%snaps).EQV..false.) cycle
          do j=1,nsets
            if (i.eq.j) cycle
            if (allocated(scluster(j)%snaps).EQV..false.) cycle
            k = i
            l = j
            if (scluster(i)%nmbrs.gt.scluster(j)%nmbrs) then
              k = j
              l = i
            end if
            kk = 0
            call cluster_to_cluster_d(rdv,scluster(j),scluster(i))
            cnt1 = cnt1 + 1
            if (rdv.le.(scluster(i)%diam+scluster(j)%diam)) then
              do jj=1,scluster(k)%nmbrs
                if ((scluster(k)%snaps(jj).eq.scluster(k)%center).AND.(scluster(k)%nmbrs.gt.1)) cycle
                cnt2 = cnt2 + 1
                call snap_to_cluster_d(coreval,scluster(l),scluster(k)%snaps(jj))
                if (coreval.lt.cradius) then
                  kk = kk + 1
                  iv1(kk) = scluster(k)%snaps(jj)
                  if (coreval.gt.scluster(l)%diam) scluster(l)%diam = coreval
                end if
              end do
              if (kk.gt.0) then
                call cluster_transferframes(scluster(k),scluster(l),iv1(1:kk),kk)
              end if
            end if
            if (allocated(scluster(i)%snaps).EQV..false.) exit
          end do
        end do
      end do
!     now shorten list and resort
      call clusters_shorten(scluster(1:nsets),nsets)
      call clusters_sort(scluster(1:nsets),nsets,afalse)
      deallocate(iv1)
      write(ilog,*) 'Done with ',nsets,' remaining of ',nsets_old,' original clusters using ',&
   & cnt1+cnt2,' additional pairwise distance evaluations.'
      write(ilog,*)
    end if
  end if
!
!
 66 format(i7,1x,i7,1x,i8,1x,g12.5,1x,g12.5,1x,g12.5)
 64 format(1000(g12.5,1x))
!
  write(ilog,*) '------------- CLUSTER SUMMARY ------------------'
  write(ilog,*) ' #       No.     Origin    Diameter     Radius      '
  do i=1,nsets
    write(ilog,66) i,scluster(i)%nmbrs,scluster(i)%center,scluster(i)%diam,scluster(i)%radius
  end do
  write(ilog,*) '------------------------------------------------'
  write(ilog,*)
  nnodes = nsets
  call quality_of_clustering(nnodes,scluster(1:nnodes),cradius,qualmet)
  call gen_graph_from_clusters(scluster,nnodes,atrue)
  call graphml_helper_for_clustering(scluster,nnodes)
  call vmd_helper_for_clustering(scluster,nnodes)
!
  end
!
!---------------------------------------------------------------------------------------------
!
subroutine hierarchical_clustering(nnodes)
!
  use clusters
  use mpistuff
  use iounit
  use math
  use interfaces
!
  implicit none
!
  integer i,j,k,l,chosei,chosej,nclusters,nclalsz,aone,allnbs,globi,nnodes
  integer, ALLOCATABLE:: tryat(:),inset(:),nbsnr(:),alllnks(:,:),iv1(:),iv2(:,:),iv3(:)
  real(KIND=4), ALLOCATABLE:: alldiss(:),tmpv(:)
  logical notdone,candid,atrue,afalse
  RTYPE rdv,cdiameter,qualmet(4),rdvtmp
!
  atrue = .true.
  afalse = .false.
  aone = 1
!
  allocate(tryat(cstored))
  allocate(inset(cstored))
  allocate(nbsnr(cstored))
  nclalsz = 10
  allocate(scluster(nclalsz))
  scluster(:)%nmbrs = 0
  scluster(:)%alsz = 0
  scluster(:)%nb = 0
  scluster(:)%nchildren = 0
  scluster(:)%nbalsz = 0
  scluster(:)%chalsz = 0
  scluster(:)%parent = 0
!
  tryat(:) = 1
  inset(:) = 0
  notdone = .true.
  nclusters = 0
  cdiameter = 2.0*cradius
  candid = .false.
  nbsnr(:) = cnblst(1:cstored)%nbs
!
  write(ilog,*)
  write(ilog,*) 'Now creating global sorted list of neighbor pairs ...'
! this first step prunes the nb-list to our desired size
  do i=1,cstored
    do j=1,cnblst(i)%nbs
      if (cnblst(i)%dis(j).gt.cdiameter) exit
    end do
    cnblst(i)%nbs = j-1
  end do
!
  allnbs = sum(cnblst(1:cstored)%nbs)
  allocate(iv2(allnbs,2))
  allocate(alldiss(allnbs))
  allocate(tmpv(allnbs))
  allocate(iv1(allnbs))
  allocate(iv3(allnbs))
  j = 1
  do i=1,cstored
    if (cnblst(i)%nbs.gt.0) then
      iv2(j:j+cnblst(i)%nbs-1,1) = i
      iv2(j:j+cnblst(i)%nbs-1,2) = cnblst(i)%idx(1:cnblst(i)%nbs)
      tmpv(j:j+cnblst(i)%nbs-1) = cnblst(i)%dis(1:cnblst(i)%nbs)
      j = j + cnblst(i)%nbs
    end if
  end do
  do j=1,allnbs
    iv1(j) = j
  end do
  call merge_sort(ldim=allnbs,up=atrue,list=tmpv(1:allnbs),olist=alldiss(1:allnbs),&
 &                ilo=aone,ihi=allnbs,idxmap=iv1(1:allnbs),olist2=iv3(1:allnbs))
  deallocate(tmpv)
  deallocate(iv1)
  allocate(alllnks(allnbs,2))
  do i=1,allnbs
    alllnks(i,:) = iv2(iv3(i),:)
  end do
  deallocate(iv2)
  deallocate(iv3)
  write(ilog,*) '... done.'
  write(ilog,*)
!
  write(ilog,*) 'Now performing hierarchical clustering by considering shortest remaining link ...'
  globi = 1
  do while (notdone.EQV..true.)
    rdv = alldiss(globi)
    chosei = alllnks(globi,1)
    chosej = alllnks(globi,2)
    globi = globi + 1
    if (globi.gt.allnbs) exit
    if (rdv.gt.cdiameter) exit
    if ((inset(chosei).gt.0).AND.(inset(chosej).gt.0)) then
!     merge inset(chosei) and inset(chosej) if permissible
      if (inset(chosei).ne.inset(chosej)) then
        if (clinkage.eq.1) then ! maximum linkage
          candid = .true.
          do i=1,scluster(inset(chosei))%nmbrs
            do j=1,scluster(inset(chosej))%nmbrs
              k = 0
              do l=1,cnblst(scluster(inset(chosej))%snaps(j))%nbs
                if (inset(cnblst(scluster(inset(chosej))%snaps(j))%idx(l)).eq.inset(chosei)) k = k + 1
                if (k.eq.scluster(inset(chosei))%nmbrs) exit
              end do
              if (k.lt.scluster(inset(chosei))%nmbrs) then
                candid = .false.
                exit
              end if
            end do
            if (candid.EQV..false.) exit
          end do
        else if (clinkage.eq.2) then ! minimum linkage
          candid = .false.
          do i=1,scluster(inset(chosei))%nmbrs
            do j=1,scluster(inset(chosej))%nmbrs
              do l=1,cnblst(scluster(inset(chosej))%snaps(j))%nbs
                if (inset(cnblst(scluster(inset(chosej))%snaps(j))%idx(l)).eq.inset(chosei)) then
                  candid = .true.
                  exit
                end if
              end do
              if (candid.EQV..true.) exit
            end do
            if (candid.EQV..true.) exit
          end do
        else if (clinkage.eq.3) then ! mean linkage
          candid = .false.
          call cluster_to_cluster_d(rdvtmp,scluster(inset(chosei)),scluster(inset(chosej)))
          if (rdvtmp.le.cradius) candid = .true.
        end if
        if (candid.EQV..true.) then
          if (scluster(inset(chosei))%nmbrs.gt.scluster(inset(chosej))%nmbrs) then
            k = inset(chosei)
            l = inset(chosej)
          else
            k = inset(chosej)
            l = inset(chosei)
          end if
          inset(scluster(l)%snaps(1:scluster(l)%nmbrs)) = k
          call join_clusters(scluster(k),scluster(l))
!         transfer last to l
          if (l.lt.nclusters) then
            call copy_cluster(scluster(nclusters),scluster(l))
            deallocate(scluster(nclusters)%snaps)
            scluster(nclusters)%alsz = 0
            scluster(nclusters)%nmbrs = 0
            nclusters = nclusters - 1
            do j=1,scluster(l)%nmbrs
              if (inset(scluster(l)%snaps(j)).ne.nclusters+1) call fexit()
              inset(scluster(l)%snaps(j)) = l
            end do
          else if (nclusters.eq.l) then
            nclusters = nclusters - 1
          end if
        end if
      end if
!   append inset(chosei) if permissible
    else if (inset(chosei).gt.0) then
      if (clinkage.eq.1) then ! maximum linkage
        k = 0
        candid = .false.
        do i=1,cnblst(chosej)%nbs
          if (inset(cnblst(chosej)%idx(i)).eq.inset(chosei)) k = k + 1
          if (k.eq.scluster(inset(chosei))%nmbrs) exit
        end do
        if (k.eq.scluster(inset(chosei))%nmbrs) candid = .true.
      else if (clinkage.eq.2) then ! minimum
        candid = .false.
!       note that tryat(chosej) has to be 1, otherwise it inset(chosej) cannot be zero
        do i=1,cnblst(chosej)%nbs
          if (inset(cnblst(chosej)%idx(i)).eq.inset(chosei)) then
            candid = .true.
            exit
          end if
        end do
      else if (clinkage.eq.3) then ! mean linkage
        candid = .false.
        call snap_to_cluster_d(rdvtmp,scluster(inset(chosei)),chosej)
        if (rdvtmp.le.cradius) candid = .true.
      end if
      if (candid.EQV..true.) then
        inset(chosej) = inset(chosei)
        call cluster_addsnap(scluster(inset(chosei)),chosej,rdv)
      end if
!   append inset(chosej) if permissible
    else if (inset(chosej).gt.0) then
      if (clinkage.eq.1) then ! maximum linkage
        k = 0
        candid = .false.
        do i=1,cnblst(chosei)%nbs
          if (inset(cnblst(chosei)%idx(i)).eq.inset(chosej)) k = k + 1
          if (k.eq.scluster(inset(chosej))%nmbrs) exit
        end do
        if (k.eq.scluster(inset(chosej))%nmbrs) candid = .true.
      else if (clinkage.eq.2) then ! minimum
        candid = .false.
!       note that tryat(chosei) has to be 1, otherwise it inset(chosei) cannot be zero
        do i=1,cnblst(chosei)%nbs
          if (inset(cnblst(chosei)%idx(i)).eq.inset(chosej)) then
            candid = .true.
            exit
          end if
        end do
      else if (clinkage.eq.3) then ! mean linkage
        candid = .false.
        call snap_to_cluster_d(rdvtmp,scluster(inset(chosej)),chosei)
        if (rdvtmp.le.cradius) candid = .true.
      end if
      if (candid.EQV..true.) then
        inset(chosei) = inset(chosej)
        call cluster_addsnap(scluster(inset(chosej)),chosei,rdv)
      end if 
!   create new cluster of size 2
    else
      nclusters = nclusters + 1
      if (nclusters.gt.nclalsz) call scluster_resizelst(nclalsz,scluster)
      call cluster_addsnap(scluster(nclusters),chosei,rdv)
      call cluster_addsnap(scluster(nclusters),chosej,rdv)
      inset(chosei) = nclusters
      inset(chosej) = nclusters
    end if
!
    tryat(chosei) = tryat(chosei) + 1
!   update nb-list counter
    if (tryat(chosei).le.cnblst(chosei)%nbs) then
      do while ((inset(chosei).eq.inset(cnblst(chosei)%idx(tryat(chosei)))).AND.(inset(chosei).gt.0))
        tryat(chosei) = tryat(chosei) + 1
        if (tryat(chosei).gt.cnblst(chosei)%nbs) exit
      end do
    end if
    if (tryat(chosej).le.cnblst(chosej)%nbs) then
      if (cnblst(chosej)%idx(tryat(chosej)).eq.chosei) tryat(chosej) = tryat(chosej) + 1
      if (tryat(chosej).le.cnblst(chosej)%nbs) then
        do while ((inset(chosej).eq.inset(cnblst(chosej)%idx(tryat(chosej))).AND.inset(chosej).gt.0))
          tryat(chosej) = tryat(chosej) + 1
          if (tryat(chosej).gt.cnblst(chosej)%nbs) exit
        end do
      end if
    end if
  end do
  write(ilog,*) '... done.'
  deallocate(alldiss)
  deallocate(alllnks)
  do i=1,cstored
    if (inset(i).le.0) then
      nclusters = nclusters + 1
      if (nclusters.gt.nclalsz) call scluster_resizelst(nclalsz,scluster)
      call cluster_addsnap(scluster(nclusters),i,rdv)
    end if
  end do
!
! shorten list
  call clusters_shorten(scluster(1:nclusters),nclusters)
  call clusters_sort(scluster(1:nclusters),nclusters,afalse)
!
 66 format(i7,1x,i7,1x,i8,1x,g12.5,1x,g12.5)
 65 format(1000(g12.5,1x))
!
  write(ilog,*) '------------- CLUSTER SUMMARY ------------------'
  write(ilog,*) ' #       No.     "Center"  Diameter     Radius      '
  do i=1,nclusters
!   get radius etc.
    call cluster_calc_params(scluster(i),cradius)
!   re-determine central structures
    call cluster_getcenter(scluster(i))
    write(ilog,66) i,scluster(i)%nmbrs,scluster(i)%center,scluster(i)%diam,scluster(i)%radius
  end do
  write(ilog,*) '------------------------------------------------'
  write(ilog,*)
!
  call quality_of_clustering(nclusters,scluster(1:nclusters),cradius,qualmet)
!
  nnodes = nclusters
  if (nnodes.gt.0) call gen_graph_from_clusters(scluster,nclusters,atrue)
  if (nnodes.gt.0) call vmd_helper_for_clustering(scluster,nclusters)
  if (nnodes.gt.0) call graphml_helper_for_clustering(scluster,nclusters)
!
  cnblst(1:cstored)%nbs = nbsnr(:)
!
  deallocate(nbsnr)
  deallocate(tryat)
  deallocate(inset)
!
end
!
!---------------------------------------------------------------------------------
!
! this subroutine generates an exact MST assuming it is provided with a nb-list
! object that holds all the necessary edges
! this routine is very memory-intensive due to the duplication of the already large nb-list object
! it is highly related to hierarchical clustering with minimum linkage and max threshold
!
subroutine gen_MST_from_nbl()
!
  use clusters
  use iounit
  use interfaces
!
  implicit none
!
  integer allnbs,i,j,k,aone,globi,chosei,chosej,ntrees,nlnks,forsz
  real(KIND=4), ALLOCATABLE:: alldiss(:),tmpv(:)
  RTYPE rdv
  integer, ALLOCATABLE:: iv2(:,:),iv1(:),iv3(:),alllnks(:,:)
  type(t_scluster), ALLOCATABLE:: it(:)
  logical atrue,notdone
!
  write(ilog,*)
  write(ilog,*) 'Now creating global sorted list of neighbor pairs ...'
!
  aone = 1
  atrue = .true.
  notdone = .true.
  allnbs = sum(cnblst(1:cstored)%nbs)/2
  allocate(iv2(allnbs,2))
  allocate(alldiss(allnbs))
  allocate(tmpv(allnbs))
  allocate(iv1(allnbs))
  allocate(iv3(allnbs))
  j = 0
  do i=1,cstored
    do k=1,cnblst(i)%nbs
      if (cnblst(i)%idx(k).gt.i) then
        j = j + 1
        iv2(j,1) = i
        iv2(j,2) = cnblst(i)%idx(k)
        tmpv(j) = cnblst(i)%dis(k)
      end if
    end do
  end do
  do j=1,allnbs
    iv1(j) = j
  end do
  call merge_sort(ldim=allnbs,up=atrue,list=tmpv(1:allnbs),olist=alldiss(1:allnbs),&
 &                ilo=aone,ihi=allnbs,idxmap=iv1(1:allnbs),olist2=iv3(1:allnbs))
  deallocate(tmpv)
  deallocate(iv1)
  allocate(alllnks(2,allnbs))
  do i=1,allnbs
    alllnks(:,i) = iv2(iv3(i),:)
!    write(*,*) alllnks(1:2,i),alldiss(i)
  end do
  deallocate(iv2)
  deallocate(iv3)
  write(ilog,*) '... done.'
  write(ilog,*)
!
  allocate(iv1(cstored))
  iv1(:) = 0
  write(ilog,*) 'Now generating MST by considering shortest remaining link and merging ...'
  globi = 1
  nlnks = 0
  ntrees = 0
  forsz = 10
  allocate(it(forsz))
  do i=1,forsz
    it(i)%alsz = 0
    it(i)%nmbrs = 0
  end do
  do while (notdone.EQV..true.)
    rdv = alldiss(globi)
    chosei = alllnks(1,globi)
    chosej = alllnks(2,globi)
    globi = globi + 1
    if (globi.gt.allnbs) exit
    if ((iv1(chosei).le.0).AND.(iv1(chosej).le.0)) then
      nlnks = nlnks + 1
      ntrees = ntrees + 1
      iv1(chosei) = ntrees
      iv1(chosej) = ntrees
      if (ntrees.gt.forsz) call scluster_resizelst(forsz,it)
      call cluster_addsnap(it(iv1(chosei)),chosei,0.0)
      call cluster_addsnap(it(iv1(chosei)),chosej,0.0)
    else if ((iv1(chosei).gt.0).AND.(iv1(chosej).gt.0)) then
      if (iv1(chosei).eq.iv1(chosej)) cycle
      nlnks = nlnks + 1
      if (it(iv1(chosei))%nmbrs.gt.it(iv1(chosej))%nmbrs) then
        j = iv1(chosej)
        do i=1,it(j)%nmbrs
          iv1(it(j)%snaps(i)) = iv1(chosei)
        end do
        call join_clusters(it(iv1(chosei)),it(j))
      else
        j = iv1(chosei)
        do i=1,it(j)%nmbrs
          iv1(it(j)%snaps(i)) = iv1(chosej)
        end do
        call join_clusters(it(iv1(chosej)),it(j))
      end if
    else if (iv1(chosei).gt.0) then
      nlnks = nlnks + 1
      iv1(chosej) = iv1(chosei)
      call cluster_addsnap(it(iv1(chosei)),chosej,0.0)
    else
      nlnks = nlnks + 1
      iv1(chosei) = iv1(chosej)
      call cluster_addsnap(it(iv1(chosej)),chosei,0.0)
    end if
    alllnks(1,nlnks) = chosei
    alllnks(2,nlnks) = chosej
    alldiss(nlnks) = rdv
    if (nlnks.eq.(cstored-1)) exit
  end do
!
  if (nlnks.ne.(cstored-1)) then
    write(ilog,*) 'Fatal. Neighbor list is insufficient to create minimum spanning tree. &
 &Increase relevant thresholds.'
    call fexit()
  end if
!
  deallocate(iv1)
  do i=1,forsz
    if (allocated(it(i)%snaps).EQV..true.) deallocate(it(i)%snaps)
    if (allocated(it(i)%tmpsnaps).EQV..true.) deallocate(it(i)%tmpsnaps)
    if (allocated(it(i)%sums).EQV..true.) deallocate(it(i)%sums)
    if (allocated(it(i)%map).EQV..true.) deallocate(it(i)%map)
    if (allocated(it(i)%children).EQV..true.) deallocate(it(i)%children)
    if (allocated(it(i)%wghtsnb).EQV..true.) deallocate(it(i)%wghtsnb)
    if (allocated(it(i)%lstnb).EQV..true.) deallocate(it(i)%lstnb)
    if (allocated(it(i)%flwnb).EQV..true.) deallocate(it(i)%flwnb)
  end do
  deallocate(it)
!
  allocate(approxmst(cstored))
  approxmst(1:cstored)%deg = 0
 587 format(' Weight of Minimum Spanning Tree: ',1x,g12.5,a)
  write(ilog,*)
  if (cdis_crit.le.2) then
    write(ilog,587) sum(alldiss(1:cstored-1)),' degrees'
  else if (cdis_crit.le.4) then
    write(ilog,587) sum(alldiss(1:cstored-1)),' '
  else if (cdis_crit.le.10) then
    write(ilog,587) sum(alldiss(1:cstored-1)),' Angstrom'
  end if
  write(ilog,*)
!
  call gen_MST(alllnks(:,1:(cstored-1)),alldiss(1:(cstored-1)),approxmst)
!
  deallocate(alldiss)
  deallocate(alllnks)
!
  write(ilog,*) '... done.'
  write(ilog,*)
!
end
!
!-----------------------------------------------------------------------------------------
!
! this subroutine generates a set of links and their lengths that constitute an approximate
! MST based on results from tree-based clustering in the birchtree object
! it then transcribes this into an adjacency list object
! the accuracy of the approximate MST depends on the number of guesses (cprogindrmax) and the properties
! of the clustering
!
subroutine gen_MST_from_treeclustering()
!
  use clusters
  use iounit
!
  implicit none
!
  integer N_NEARS
  parameter(N_NEARS=5)
!
  integer i,j,k,r,e,i1,i2,l,m,mm,ixx,ixx2,cursnp
  integer cntloop,cntloop2,cntguess,tlstsz,lastcl,startcl
  logical skippy
  integer ntrees,oldntrees ! number of active trees
  type(t_progindextree), ALLOCATABLE:: trees(:)
  integer, ALLOCATABLE:: tmpix(:,:)
  RTYPE, ALLOCATABLE:: tmpdis(:,:)
  integer, ALLOCATABLE:: snap2tree(:) ! index of its tree for each snapshot
  RTYPE random ! random number
  RTYPE jkdist ! distance between snapshots i and k
  RTYPE t1,t2
  integer forestsize ! number of active trees
  integer, ALLOCATABLE:: treeptr(:) ! used to merge trees
  integer boruvkasteps ! counts how many boruvka steps are needed
  integer, ALLOCATABLE::  edgelst(:,:) ! array holding the selected edges of the current Boruvka step.
                                       ! Indices: edgenumber (1-nedges), endpoint (1-2)
  real(KIND=4), ALLOCATABLE::  ledgelst(:) ! array holding the lengths of the selected edges of the current Boruvka step
  integer nedges ! number of edges in edgelst
  integer nmstedges ! number of edges in the approximate MST
  integer, ALLOCATABLE:: snap2clus(:,:) ! map (snapshot, tree level) -> cluster based on birchtree
  integer, ALLOCATABLE:: snap2lst(:)  ! map (snapshot) -> list entry in list of unique clusters "testlst"
  integer, ALLOCATABLE:: testlst(:,:) ! temporary list of eligible candidates for random search with annotations
  
  type(t_adjlist), ALLOCATABLE:: snplst(:) ! for each tree, the hierarchy level to start from
  logical, ALLOCATABLE:: entered(:) ! membership array
  integer a,b ! numbers of snapshots, used for sophisticated nearest neighbor guessing
  logical emptylevel ! used for sophisticated nearest neighbor guessing
  integer kk ! local snapshot index within cluster, used for sophisticated nearest neighbor guessing
  integer, ALLOCATABLE:: mstedges(:,:) ! contains all the edges in the approximate MST
                                       ! Indices: edgenumber (1-nedges), endpoint (1-2)
  real(KIND=4), ALLOCATABLE:: lmstedges(:) ! contains the lengths of all the edges in the approximate MST
  integer(KIND=8) testcnt
!
  write(ilog,*)
  write(ilog,*) 'Now generating approximate MST based on tree-based clustering ...'
!
!
  ntrees = cstored
  allocate(trees(ntrees))
  allocate(treeptr(ntrees))
  allocate(snap2tree(cstored))
  do i=1,ntrees
    allocate(trees(i)%snaps(1))
  end do
  allocate(snap2clus(cstored,c_nhier+1))
  allocate(mstedges(2,cstored-1))
  allocate(lmstedges(cstored-1))
  allocate(tmpix(cstored,N_NEARS))
  allocate(tmpdis(cstored,N_NEARS)) ! N_NEARS+1))
!
  trees(:)%nsnaps = 1
  trees(:)%mine(1) = -1
  trees(:)%mine(2) = -1
  trees(:)%nsibalsz = 0
  forestsize = ntrees
  boruvkasteps = 0
  nmstedges = 0
  cntloop = 0
  cntguess = 0
  testcnt = 0
  a = maxval(birchtree(1:(c_nhier+1))%ncls)
  allocate(snap2lst(max(cstored,a)))
  allocate(entered(a))
  allocate(snplst(a))
  allocate(testlst(a,3))
  snap2lst(:) = 0
  snplst(:)%deg = 0
  snplst(:)%alsz = 0
  entered(:) = .false.
!
! Generate map (snapshot, tree level) -> cluster based on birchtree; in the process remove double entries for non-terminal levels
! arrange trees such that they initially follow clustering structure
  mm = 0
  do l=c_nhier+1,1,-1
    snap2clus(:,l) = 0
    do k=1,birchtree(l)%ncls
      do j=1,birchtree(l)%cls(k)%nmbrs
        if (snap2clus(birchtree(l)%cls(k)%snaps(j),l).gt.0) then
          if (birchtree(l)%cls(snap2clus(birchtree(l)%cls(k)%snaps(j),l))%nmbrs.gt.birchtree(l)%cls(k)%nmbrs) then
            call cluster_removesnap(birchtree(l)%cls(k),birchtree(l)%cls(k)%snaps(j))
          else
            call cluster_removesnap(birchtree(l)%cls(snap2clus(birchtree(l)%cls(k)%snaps(j),l)),birchtree(l)%cls(k)%snaps(j))
            snap2clus(birchtree(l)%cls(k)%snaps(j),l) = k
          end if
        else
          snap2clus(birchtree(l)%cls(k)%snaps(j),l) = k
        end if
        if (l.eq.(c_nhier+1)) then
          mm = mm + 1
          trees(mm)%snaps(1) = birchtree(l)%cls(k)%snaps(j)
          snap2tree(birchtree(l)%cls(k)%snaps(j)) = mm
        end if
      end do
    end do
  end do
  tmpdis(:,:) = HUGE(tmpdis(1,1)) ! initialize
  tmpix(:,:) = 0                  ! initialize
!
  do while (forestsize.ge.2)
    call CPU_time(t1)
    nedges = 0
    allocate(edgelst(forestsize,2))
    allocate(ledgelst(forestsize))
    do i=1,ntrees
      lastcl = 0
      if (trees(i)%nsnaps.gt.0) then ! tree i is still active
! determine shortest edge from tree i to another tree:
        l = c_nhier + 1 ! start on the leaf level
        if ((boruvkasteps.eq.0).AND.(i.gt.1)) then
          lastcl = snap2clus(trees(i-1)%snaps(1),l)
        end if
        emptylevel = .true.
        do while (emptylevel.EQV..true.)
! first assemble a list of unique clusters that the current tree spans into
          tlstsz = 0
          do j=1,trees(i)%nsnaps
            mm = abs(snap2clus(trees(i)%snaps(j),l))
            if (entered(mm).EQV..false.) then
              entered(mm) = .true.
              tlstsz = tlstsz + 1
              testlst(tlstsz,1) = mm
              testlst(tlstsz,2) = j
              snap2lst(mm) = tlstsz
              testlst(tlstsz,3) = 0
            end if
            testlst(snap2lst(mm),3) = testlst(snap2lst(mm),3) + 1
          end do
          emptylevel = .false.
! scan the list of unique clusters to see which ones are usable
          do j=1,tlstsz
            a = birchtree(l)%cls(testlst(j,1))%nmbrs
            if (testlst(j,3).ge.a) then
              emptylevel = .true.
              entered(testlst(j,1)) = .false.
              snplst(j)%deg = 0
              cycle
            end if
            if (snplst(j)%alsz.lt.a) then ! increase size if needed
              if (allocated(snplst(j)%adj).EQV..true.) deallocate(snplst(j)%adj)
              allocate(snplst(j)%adj(a))
              snplst(j)%alsz = a
            end if
            if ((boruvkasteps.gt.0).OR.(a.eq.1).OR.(i.eq.1).OR.(lastcl.ne.snap2clus(trees(i)%snaps(1),l)).OR.&
 &              (l.lt.(c_nhier+1))) then
              snplst(j)%deg = 0
              startcl = i
              if ((2*testlst(j,3).lt.a).AND.((a-testlst(j,3)).gt.cprogindrmax).AND.(l.le.c_nhier)) then
!               do nothing if we know that the cluster is not leaf, has capacity for random drawing, and is less than
!               half-occupied by the current tree
              else
                do m=1,a
                  if (snap2tree(birchtree(l)%cls(testlst(j,1))%snaps(m)).ne.i) then
                    snplst(j)%deg = snplst(j)%deg + 1
                    snplst(j)%adj(snplst(j)%deg) = m
                  end if
                end do
              end if
            else
              snplst(j)%adj(i-startcl) = i-startcl
            end if
          end do
!         descend toward root if all clusters at this level are unusable
          skippy = .true.
          do j=1,tlstsz
            if (entered(testlst(j,1)).EQV..true.) then
              skippy = .false.
              exit
            end if
          end do
          if (skippy.EQV..true.) then
            l = l - 1
            cycle
          end if
!         for those snaps that do have a finite search space, perform search
          do j=1,trees(i)%nsnaps
            cursnp = trees(i)%snaps(j)
            if (snap2clus(cursnp,l).lt.0) cycle
            if (entered(snap2clus(cursnp,l)).EQV..false.) cycle
            a = snap2lst(snap2clus(cursnp,l))
            if ((snplst(a)%deg.gt.0).AND.(snplst(a)%deg.lt.cprogindrmax)) then
!   Loop over all eligible snapshots in cluster of snapshot j of tree i:
              cntloop = cntloop + 1
              cntloop2 = 0
              do mm=1,snplst(a)%deg
                m = snplst(a)%adj(mm) 
                cntloop2 = cntloop2 + 1
                if (cntloop2.gt.cprogindrmax) then
                  write(ilog,*) 'NICO: Bug! Loop too big.'
                  write(ilog,*) 'Tree=',i
                  write(ilog,*) 'j=',j
                  write(ilog,*) 'Level=',l
                  write(ilog,*) 'Cluster=',snap2clus(cursnp,l)
                  write(ilog,*) 'cntloop2=',cntloop2
                  write(ilog,*) 'a=',a
                  write(ilog,*) 'trees(i)%nsnaps=',trees(i)%nsnaps
                  write(ilog,*) 'b=',b
                  write(ilog,*) 'In cluster: (snapshot, tree)'
                  do kk=1,snplst(a)%deg
                    k = snplst(a)%adj(kk)
                    write(ilog,*) birchtree(l)%cls(snap2clus(cursnp,l))%snaps(k),&
&                                 snap2tree(birchtree(l)%cls(snap2clus(cursnp,l))%snaps(k))
                  end do
                  write(ilog,*) 'In tree: (snapshot, cluster)'
                  do k=1,trees(i)%nsnaps
                    write(ilog,*) trees(i)%snaps(k),snap2clus(trees(i)%snaps(k),l)
                  end do
                  call fexit()
                end if
                k = birchtree(l)%cls(snap2clus(cursnp,l))%snaps(m) ! global snapshot
!   compute length of edge (cursnp,k):
                call snap_to_snap_d(jkdist,cursnp,k)
                testcnt = testcnt + 1
!   update shortest eligible edges per snap
                do ixx=1,N_NEARS
                  if (tmpix(cursnp,ixx).eq.k) exit
                  if (jkdist.lt.tmpdis(cursnp,ixx)) then
                    do ixx2=N_NEARS,ixx+1,-1 ! shift
                      tmpix(cursnp,ixx2) = tmpix(cursnp,ixx2-1)
                      tmpdis(cursnp,ixx2) = tmpdis(cursnp,ixx2-1)
                    end do
                    tmpix(cursnp,ixx) = k
                    tmpdis(cursnp,ixx) = jkdist 
                    exit
                  end if
                end do
                do ixx=1,N_NEARS
                  if (tmpix(k,ixx).eq.cursnp) exit
                  if (jkdist.lt.tmpdis(k,ixx)) then
                    do ixx2=N_NEARS,ixx+1,-1 ! shift
                      tmpix(k,ixx2) = tmpix(k,ixx2-1)
                      tmpdis(k,ixx2) = tmpdis(k,ixx2-1)
                    end do
                    tmpix(k,ixx) = cursnp
                    tmpdis(k,ixx) = jkdist 
                    exit
                  end if
                end do
              end do
            else
!   Guess nearest neighbor of snapshot j of tree i in its cluster cprogindrmax times:
              cntguess = cntguess + 1
!   this is the variant where the list is not assembled (not cost efficient)
              if (snplst(a)%deg.eq.0) then
                r = 0
                do while (r.le.cprogindrmax)
                  kk = ceiling(random()*birchtree(l)%cls(snap2clus(trees(i)%snaps(j),l))%nmbrs)
                  k = birchtree(l)%cls(snap2clus(trees(i)%snaps(j),l))%snaps(kk) ! nearest neighbor candidate
                  if (snap2tree(k).eq.i) cycle
!   compute length of edge (trees(i)%snaps(j),k):
                  call snap_to_snap_d(jkdist,trees(i)%snaps(j),k)
                  testcnt = testcnt + 1
                  r = r + 1
!   update shortest eligible edges per snap
                  do ixx=1,N_NEARS
                    if (tmpix(cursnp,ixx).eq.k) exit
                    if (jkdist.lt.tmpdis(cursnp,ixx)) then
                      do ixx2=N_NEARS,ixx+1,-1 ! shift
                        tmpix(cursnp,ixx2) = tmpix(cursnp,ixx2-1)
                        tmpdis(cursnp,ixx2) = tmpdis(cursnp,ixx2-1)
                      end do
                      tmpix(cursnp,ixx) = k
                      tmpdis(cursnp,ixx) = jkdist 
                      exit
                    end if
                  end do
                  do ixx=1,N_NEARS
                    if (tmpix(k,ixx).eq.cursnp) exit
                    if (jkdist.lt.tmpdis(k,ixx)) then
                      do ixx2=N_NEARS,ixx+1,-1 ! shift
                        tmpix(k,ixx2) = tmpix(k,ixx2-1)
                        tmpdis(k,ixx2) = tmpdis(k,ixx2-1)
                      end do
                      tmpix(k,ixx) = cursnp
                      tmpdis(k,ixx) = jkdist 
                      exit
                    end if
                  end do
                end do
              else
!   this is the variant where the list is assembled and hence no checks are required and no delays are possible
                do r=1,cprogindrmax
                  kk = ceiling(random()*snplst(a)%deg)
                  k = birchtree(l)%cls(snap2clus(trees(i)%snaps(j),l))%snaps(snplst(a)%adj(kk)) ! nearest neighbor candidate
!   compute length of edge (trees(i)%snaps(j),k):
                  call snap_to_snap_d(jkdist,trees(i)%snaps(j),k)
                  testcnt = testcnt + 1
!   update shortest eligible edges per snap:
                  do ixx=1,N_NEARS
                    if (tmpix(cursnp,ixx).eq.k) exit
                    if (jkdist.lt.tmpdis(cursnp,ixx)) then
                      do ixx2=N_NEARS,ixx+1,-1 ! shift
                        tmpix(cursnp,ixx2) = tmpix(cursnp,ixx2-1)
                        tmpdis(cursnp,ixx2) = tmpdis(cursnp,ixx2-1)
                      end do
                      tmpix(cursnp,ixx) = k
                      tmpdis(cursnp,ixx) = jkdist 
                      exit
                    end if
                  end do
                  do ixx=1,N_NEARS
                    if (tmpix(k,ixx).eq.cursnp) exit
                    if (jkdist.lt.tmpdis(k,ixx)) then
                      do ixx2=N_NEARS,ixx+1,-1 ! shift
                        tmpix(k,ixx2) = tmpix(k,ixx2-1)
                        tmpdis(k,ixx2) = tmpdis(k,ixx2-1)
                      end do
                      tmpix(k,ixx) = cursnp
                      tmpdis(k,ixx) = jkdist 
                      exit
                    end if
                  end do
                end do
              end if
            end if
            snap2clus(trees(i)%snaps(j),1:l) = -snap2clus(trees(i)%snaps(j),1:l)
          end do
          do j=1,tlstsz
            entered(testlst(j,1)) = .false.
          end do
          l = max(1,l - 1)
        end do
        do j=1,trees(i)%nsnaps
          if (minval(snap2clus(trees(i)%snaps(j),:)).gt.0) then
            write(ilog,*) 'Fatal. This means that at least one snapshot during progress index-based clustering was not &
 &processed. Please report this bug.'
            call fexit()
          end if
          snap2clus(trees(i)%snaps(j),:) = abs(snap2clus(trees(i)%snaps(j),:))
        end do
      end if
    end do
!    do j=1,cstored
!      tmpdis(j,N_NEARS+1) = 0.
!      ixx2 = 0
!      do ixx=1,N_NEARS
!        if (tmpix(j,ixx).le.0) exit
!        ixx2 = ixx2 + 1
!        tmpdis(j,N_NEARS+1) = tmpdis(j,N_NEARS+1) + tmpdis(j,ixx)
!      end do
!      tmpdis(j,N_NEARS+1) = tmpdis(j,N_NEARS+1)/(1.0*ixx2)
!    end do
    do i=1,ntrees
!     per tree, identify the edge we want to use for the spanning tree
      nedges = nedges + 1
      ledgelst(nedges) = HUGE(ledgelst(nedges))
      do j=1,trees(i)%nsnaps
        cursnp = trees(i)%snaps(j)
        do ixx=1,N_NEARS
          jkdist = tmpdis(cursnp,1)  !  + tmpdis(cursnp,N_NEARS+1) + tmpdis(tmpix(cursnp,ixx),N_NEARS+1)
          if (jkdist.lt.ledgelst(nedges)) then
            edgelst(nedges,1) = cursnp
            edgelst(nedges,2) = tmpix(cursnp,1)
            ledgelst(nedges) = jkdist ! sum(tmpdis(cursnp,1:5))
          end if
        end do
      end do
    end do
!
! merge trees (& update approximate minimum spanning tree):
    do j=1,ntrees
      treeptr(j) = j
    end do
    trees(1:ntrees)%nsiblings = 0
    do e=1,nedges
      i1 = treeptr(snap2tree(edgelst(e,1))) ! tree index of first edge endpoint
      i2 = treeptr(snap2tree(edgelst(e,2))) ! tree index of second edge endpoint
      if (i1.ne.i2) then ! edge does not introduce cycle
! add corresponding edge to approximate minimum spanning tree:
        nmstedges = nmstedges+1
        if (nmstedges.ge.cstored) then
          write(ilog,*) 'Fatal. The number of edges in the approximate minimum spanning tree is not correct. Please report &
 &this bug.'
          call fexit()
        end if
        mstedges(1,nmstedges) = edgelst(e,1)
        mstedges(2,nmstedges) = edgelst(e,2)
        lmstedges(nmstedges) = ledgelst(e)
! update pointer for tree and all its siblings obtained through prior merge operations
        treeptr(i1) = i2
        trees(i2)%nsiblings = trees(i2)%nsiblings + 1
        if (trees(i2)%nsiblings.gt.trees(i2)%nsibalsz) call pidxtree_growsiblings(trees(i2))
        trees(i2)%siblings(trees(i2)%nsiblings) = i1
        do j=1,trees(i1)%nsiblings
          mm = trees(i1)%siblings(j)
          trees(i2)%nsiblings = trees(i2)%nsiblings + 1
          if (trees(i2)%nsiblings.gt.trees(i2)%nsibalsz) call pidxtree_growsiblings(trees(i2))
          trees(i2)%siblings(trees(i2)%nsiblings) = mm
          treeptr(mm) = i2
        end do
! decrease forest size
        forestsize = forestsize-1
      end if
    end do
! reassign snap2tree based on pointer
    do j=1,cstored
      snap2tree(j) = treeptr(snap2tree(j))
    end do
! now construct a new set of trees from the up-to-date snap2tree
    mm = 0
    treeptr(1:cstored) = 0
    oldntrees = ntrees
    ntrees = 0
    do j=1,cstored
      if (treeptr(snap2tree(j)).eq.0) then
        ntrees = ntrees + 1
        trees(ntrees)%mine(1) = snap2tree(j) ! hijack
        trees(snap2tree(j))%mine(2) = ntrees ! hijack
      end if
      treeptr(snap2tree(j)) = treeptr(snap2tree(j)) + 1
    end do
    do mm=1,ntrees
      j = trees(mm)%mine(1)
      if (allocated(trees(mm)%snaps).EQV..true.) deallocate(trees(mm)%snaps)
      allocate(trees(mm)%snaps(treeptr(j)))
      trees(mm)%nsnaps = 0
    end do
    do j=1,cstored
      mm = trees(snap2tree(j))%mine(2)
      trees(mm)%nsnaps = trees(mm)%nsnaps + 1
      trees(mm)%snaps(trees(mm)%nsnaps) = j
      snap2tree(j) = trees(snap2tree(j))%mine(2)
    end do
    do j=ntrees+1,oldntrees
      if (allocated(trees(j)%snaps).EQV..true.) deallocate(trees(j)%snaps)
      if (allocated(trees(j)%siblings).EQV..true.) deallocate(trees(j)%siblings)
    end do
!
! manage stored guesses
    do j=1,cstored
      ixx2 = 0
      do ixx=1,N_NEARS
        if (tmpix(j,ixx).gt.0) then
          if (snap2tree(j).ne.snap2tree(tmpix(j,ixx))) then
            ixx2 = ixx2 + 1
            tmpix(j,ixx2) = tmpix(j,ixx)
            tmpdis(j,ixx2) = tmpdis(j,ixx)
          end if
        end if
      end do
      do ixx=max(1,ixx2),N_NEARS
        tmpix(j,ixx) = 0
        tmpdis(j,ixx) = HUGE(tmpdis(j,ixx))
      end do
    end do
!
! reset edgelist
    deallocate(edgelst)
    deallocate(ledgelst)
    boruvkasteps = boruvkasteps+1
    call CPU_time(t2)
 567 format('... time for Boruvka stage ',i4,': ',g11.3,'s ...')
    write(ilog,567) boruvkasteps,t2-t1 
  end do
!
!  write(ilog,*) 'NICO:   boruvkasteps = ',boruvkasteps
  if (cstored-1.ne.nmstedges) then
    write(ilog,*) 'Fatal. The number of edges in the approximate minimum spanning tree is not correct. Please report &
 &this bug.'
    call fexit()
  end if
!
  allocate(approxmst(cstored))
  approxmst(1:cstored)%deg = 0
 587 format(' Weight of Short Spanning Tree: ',1x,g12.5,a)
  write(ilog,*)
  if (cdis_crit.le.2) then
    write(ilog,587) sum(lmstedges(1:cstored-1)),' degrees'
  else if (cdis_crit.le.4) then
    write(ilog,587) sum(lmstedges(1:cstored-1)),' '
  else if (cdis_crit.le.10) then
    write(ilog,587) sum(lmstedges(1:cstored-1)),' Angstrom'
  end if
  write(ilog,*)
  call gen_MST(mstedges,lmstedges,approxmst)
!
  do i=1,ntrees
    if (allocated(trees(i)%snaps).EQV..true.) deallocate(trees(i)%snaps)
    if (allocated(trees(i)%siblings).EQV..true.) deallocate(trees(i)%siblings)
  end do
  deallocate(entered)
  deallocate(testlst)
  do i=1,size(snplst)
    if (allocated(snplst(i)%adj).EQV..true.) deallocate(snplst(i)%adj)
  end do
  deallocate(tmpdis)
  deallocate(tmpix)
  deallocate(snplst)
  deallocate(treeptr)
  deallocate(trees)
  deallocate(lmstedges)
  deallocate(mstedges)
  deallocate(snap2tree)
  deallocate(snap2lst)
  deallocate(snap2clus)
!
 77 format(a,20(i18,1x))
  write(ilog,*) '... done after ',testcnt,' additional distance evaluations.'
  write(ilog,*)
!
end
!
!--------------------------------------------------------------------------------
!
! this routines transcribes a list of edges with their distances, that effectively 
! describe the MST, into an array of adjacency list objects
!
subroutine gen_MST(mstedges,lmstedges,mst)
!
  use clusters
!
  implicit none
!
  integer e,v
  type(t_adjlist) mst(cstored)
  integer mstedges(2,cstored-1)
  REAL(kind=4) lmstedges(cstored-1)
!
! transform edgelist of MST (mstedges) to adjacencylist:
  do e=1,cstored
    allocate(mst(e)%adj(1))
    allocate(mst(e)%dist(1))
  end do
  do e=1,cstored-1
    do v=1,2
      if (mst(mstedges(v,e))%deg.gt.0) call extend_adjlst_byone(mst(mstedges(v,e)))
      mst(mstedges(v,e))%deg = mst(mstedges(v,e))%deg + 1
      mst(mstedges(v,e))%adj(mst(mstedges(v,e))%deg) = mstedges(3-v,e)
      mst(mstedges(v,e))%dist(mst(mstedges(v,e))%deg) = lmstedges(e)
    end do
  end do
!
end
!
!------------------------------------------------------------------------------------------------------
!
! a simple helper to grow as conservatively (and slowly) as possible an adjacency list object
!
subroutine extend_adjlst_byone(mstnode)
!
  use clusters
!
  implicit none
!
  integer, ALLOCATABLE:: itmp1(:)
  real(KIND=4), ALLOCATABLE:: rtmp1(:)
  type(t_adjlist) mstnode
!
  allocate(itmp1(mstnode%deg))
  allocate(rtmp1(mstnode%deg))
!
  itmp1(1:mstnode%deg) = mstnode%adj(1:mstnode%deg)
  rtmp1(1:mstnode%deg) = mstnode%dist(1:mstnode%deg)
  deallocate(mstnode%adj)
  deallocate(mstnode%dist)
!
  allocate(mstnode%adj(mstnode%deg+1))
  allocate(mstnode%dist(mstnode%deg+1))
  mstnode%adj(1:mstnode%deg) = itmp1(1:mstnode%deg)
  mstnode%dist(1:mstnode%deg) = rtmp1(1:mstnode%deg)
!
  deallocate(itmp1)
  deallocate(rtmp1)
!
end
!
!-----------------------------------------------------------------------------------------------------
!
! a subroutine to contract terminal vertices into their parent by artificially setting the distance 
! to a negative value
!
subroutine contract_mst(alst,nrnds)
!
  use clusters
  use iounit, ONLY: ilog
!
  implicit none
!
  integer, INTENT(IN):: nrnds
!
  integer i,j,kk,thej,ll
  type(t_adjlist) alst(cstored)
!
  logical, ALLOCATABLE:: terminal(:)
!
  allocate(terminal(cstored))
  terminal(:) = .false.
!
  do i=1,cstored
    if (alst(i)%deg.eq.1) terminal(i) = .true.
  end do
!
  do kk=1,nrnds
    do i=1,cstored
      if (terminal(i).EQV..true.) then
        do j=1,alst(i)%deg
          if (alst(i)%dist(j).ge.0.0) then
            thej = alst(i)%adj(j)
            if (alst(i)%dist(j).eq.0.0) then
              alst(i)%dist(j) = -10.0*TINY(alst(i)%dist(j))
            else
              alst(i)%dist(j) = -1.0/alst(i)%dist(j) 
            end if
            exit
          end if
        end do
        do ll=1,alst(thej)%deg
          if (alst(thej)%adj(ll).eq.i) then
            alst(thej)%dist(ll) = alst(i)%dist(j)
            exit
          end if
        end do
      end if
    end do
!   reset terminal status
    terminal(:) = .false.
    ll = 0
    do i=1,cstored
      thej = 0
      do j=1,alst(i)%deg
        if (alst(i)%dist(j).ge.0.0) thej = thej + 1
      end do
      if (thej.eq.1) then
        terminal(i) = .true.
        ll = ll + 1
      end if
    end do
    if (ll.eq.0) exit
  end do
!
  write(ilog,*) 
  write(ilog,104) cstored-ll,100.0*(cstored-ll)/(1.0*cstored)
 104 format('The number of promoted (folded) edges in the MST is ',i10,' (',f8.3,' %).')
  deallocate(terminal)
!
end
!
!
!------------------------------------------------------------------------------------------------------
!
! using PRIM's algorithm, this routine will trace the MST adjacency list to derive the progress index per se
! (along with the minimal distance to the current set (distv), the inverse map (invvec) and the parent/source vector (iv2))
!
subroutine gen_progind_from_adjlst(alst,starter,progind,distv,invvec,iv2)
!
  use clusters
  use iounit
!
  implicit none
!
  type(t_adjlist) alst(cstored)
  logical, ALLOCATABLE:: added(:), inprogind(:)
  integer, ALLOCATABLE:: heap(:),hsource(:)
  integer heapsize,lprogind,invvec(cstored+2),progind(cstored),iv2(cstored),j,starter
  real(KIND=4), ALLOCATABLE:: key(:)
  real(KIND=4) distv(cstored)
!
  allocate(inprogind(cstored))
  allocate(added(cstored))
  allocate(key(cstored))
  allocate(heap(cstored))
  allocate(hsource(cstored))
  added(:) = .false.
  inprogind(:) = .false.
  distv(:) = 0.0
  key(:) = 0.0
  hsource(:) = 0
!
! add first snapshot to progress index:
  if (starter.le.cstored) then
    progind(1) = starter
  else
    progind(1) = 1
    write(ilog,*) 'Warning. The snapshot index requested is not available (there are ',cstored,' snapshots in memory). &
 &Using first one instead.'
  end if
  lprogind = 1
  distv(lprogind) = 0.0
  added(progind(lprogind)) = .true.
  inprogind(progind(lprogind)) = .true.
  invvec(progind(lprogind)+1) = 1
  iv2(lprogind) = progind(lprogind)
! build heap:
  heapsize = alst(progind(lprogind))%deg
  heap(1:heapsize) = alst(progind(lprogind))%adj(:)
  key(1:heapsize) = alst(progind(lprogind))%dist(:)
  hsource(1:heapsize) = progind(lprogind)
  call hbuild(heap(1:heapsize),heapsize,key(1:heapsize),hsource(1:heapsize))
  do j=1,alst(progind(lprogind))%deg
    added(alst(progind(lprogind))%adj(j)) = .true.
  end do
  do while (lprogind.lt.cstored)
    do j=1,alst(progind(lprogind))%deg
! add neighbors of last snapshot in progress index to heap
      if (added(alst(progind(lprogind))%adj(j)).EQV..false.) then
        call hinsert(heap(1:(heapsize+1)),heapsize,key(1:(heapsize+1)),hsource(1:(heapsize+1)),&
 &                   alst(progind(lprogind))%adj(j),alst(progind(lprogind))%dist(j),progind(lprogind))
        added(alst(progind(lprogind))%adj(j)) = .true.
      end if
    end do
! append next snapshot to progind():
    lprogind = lprogind+1
    progind(lprogind) = heap(1)
    iv2(lprogind) = hsource(1)
    distv(lprogind) = key(1)
    invvec(heap(1)+1) = lprogind
    inprogind(progind(lprogind)) = .true.
! remove added snapshot from heap:
    call hremovemin(heap(1:heapsize),heapsize,key(1:heapsize),hsource(1:heapsize))
  end do
!
  do j=1,cstored
    if (distv(j).lt.0.0) then
      if (distv(j).ge.-10.0*TINY(distv(j))) then
        distv(j) = 0.0
      else
        distv(j) = -1.0/distv(j)
      end if
    end if
  end do
!
  deallocate(key)
  deallocate(added)
  deallocate(inprogind)
  deallocate(heap)
  deallocate(hsource)
!
end
!
!------------------------------------------------------------------------------------------------------
!
subroutine do_prog_index()
!
  use clusters
  use iounit
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer refi,i,kk,aone,idxc,ii,jj
  RTYPE random
  integer, ALLOCATABLE:: progind(:),invvec(:),cutv(:),iv2(:),proflst(:)
  real(KIND=4), ALLOCATABLE:: distv(:)
  logical candid,profile_min,exists
  integer tl2,iu,freeunit
  character(12) nod2
  character(MAXSTRLEN) fn
#ifdef ENABLE_MPI
  integer tl
  character(3) nod
#endif
!
  write(ilog,*)
  write(ilog,*) 'Now deriving progress index and computing MFPT-related quantities ...'
!
  allocate(distv(cstored))
  allocate(progind(cstored))
  allocate(invvec(cstored+2))
  allocate(iv2(cstored))
  allocate(proflst(max(1,cstored/csivmin)))
  if (csivmin.ge.cstored) then
    write(ilog,*) 'Warning. Setting for FMCSC_CBASINMAX is inappropriate (too large). Adjusting to '&
 &,max(cstored/10,1),'.'
    csivmin = max(1,cstored/10)
  end if
  if (((cprogindstart.eq.-1).OR.(cprogindstart.eq.-2)).AND.(cprogindex.eq.2)) then
    write(ilog,*) 'Using representative snapshot of largest cluster (#',birchtree(c_nhier+1)%cls(1)%center,&
 &') as reference for progress index (setting -1/-2).'
    cprogindstart = birchtree(c_nhier+1)%cls(1)%center ! center of largest cluster in hierarchical tree
  end if
  if (cprogindstart.eq.-3) then
    write(ilog,*) 'Using first snapshot as reference for progress index (setting -3).'
    cprogindstart = 1 ! first snapshot
  end if
  if (cprogindstart.gt.cstored) then
    write(ilog,*) 'Warning. Setting for FMCSC_CPROGINDSTART requests a snapshot that does not exist. Remember &
 &that numbering is with respect to the snapshots currently in memory. Using first snapshot instead.'
    cprogindstart = 1
  end if
  if (cprogpwidth.ge.cstored/2) then
    write(ilog,*) 'Warning. Setting for FMCSC_CPROGINDWIDTH is inappropriate (too large). Adjusting to '&
 &,max(cstored/10,1),'.'
    cprogpwidth = max(1,cstored/10)
  end if
!
  aone = 1
!
  if (cprogfold.gt.0) then
    call contract_mst(approxmst,cprogfold)
  end if
!
  if (cprogindstart.eq.0) then
    allocate(cutv(cstored))
    invvec(1) = 3*cstored
    invvec(cstored+2) = 6*cstored
    refi = ceiling(random()*cstored)
!
    call gen_progind_from_adjlst(approxmst,refi,progind,distv,invvec,iv2)
!
!   generate N(A->B) + N(B->A) with boundary condition allB
    cutv(1) = 2 
    do i=2,cstored
      kk = progind(i) + 1
      if ((invvec(kk+1).lt.i).AND.(invvec(kk-1).lt.i)) then
        cutv(i) = cutv(i-1) - 2
      else if ((invvec(kk+1).gt.i).AND.(invvec(kk-1).gt.i)) then
        cutv(i) = cutv(i-1) + 2
      else
        cutv(i) = cutv(i-1)
      end if
    end do
    idxc = 0
    do i=csivmin+1,cstored-csivmin
      candid = profile_min(cutv(1:cstored),cstored,i,csivmin,csivmax,aone,kk)
      if (candid.EQV..true.) then
        idxc = idxc + 1
        proflst(idxc) = i
      end if
    end do
    if (idxc.le.0) then
      write(ilog,*) 'Warning. Automatic identification of starting snapshots for profile &
 &generation failed. Creating a single profile from snapshot #1.'
      idxc = 1
      proflst(idxc) = 1
    end if
  else if (cprogindstart.gt.0) then
    idxc = 1
    proflst(idxc) = cprogindstart
  end if
!
  do i=1,idxc
    invvec(1) = 3*cstored
    invvec(cstored+2) = 6*cstored
!
    call gen_progind_from_adjlst(approxmst,proflst(i),progind,distv,invvec,iv2)
    tl2 = 12
    call int2str(proflst(i),nod2,tl2)
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      tl = 3
      call int2str(myrank,nod,tl)
      fn = 'N_'//nod(1:tl)//'_PROGIDX_'//nod2(1:tl2)//'.dat'
    else if (use_MPIAVG.EQV..true.) then
      fn = 'PROGIDX_'//nod2(1:tl2)//'.dat'
    end if
#else
    fn = 'PROGIDX_'//nod2(1:tl2)//'.dat'
#endif
    call strlims(fn,ii,jj)
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      iu = freeunit()
      open(unit=iu,file=fn(ii:jj),status='old',position='append')
      close(unit=iu,status='delete')
    end if
    iu=freeunit()
    open(unit=iu,file=fn(ii:jj),status='new')
    call gen_manycuts(progind,distv,invvec,iv2,iu,cprogpwidth)
    close(unit=iu)
  end do
!
  deallocate(proflst)
  deallocate(iv2)
  deallocate(invvec)
  deallocate(distv)
  deallocate(progind)
!
  write(ilog,*) '... done.'
  write(ilog,*)
!
end
!
!----------------------------------------------------------------------------------------------------
!
subroutine gen_manycuts(setis,distv,invvec,ivec2,iu,pwidth)
!
  use clusters
  use iounit
!
  implicit none
!
  integer tmat(3,3),vv1(3),vv2(3),vv3(3)
  integer pvo(5),pvn(5),iu,i,j,k,pwidth,ii,jj,kk,ll
  integer invvec(cstored+2),ivec2(cstored),setis(cstored)
  real(KIND=4) distv(cstored)
  integer, ALLOCATABLE:: cutv(:)
  logical, ALLOCATABLE:: ibrkx(:)
!
  allocate(cutv(cstored))
  allocate(ibrkx(cstored+1))
  ibrkx(:) = .false.
  do i=1,ntbrks2
    ibrkx(trbrkslst(i)+1) = .true.
  end do
!
  cutv(1) = 2
  if (ibrkx(cstored+1).EQV..true.) cutv(1) = 1 
  do i=2,cstored
    kk = setis(i) + 1
    if ((invvec(kk+1).lt.i).AND.(invvec(kk-1).lt.i)) then
      cutv(i) = cutv(i-1) - 2
      if (ibrkx(kk-1).EQV..true.) cutv(i) = cutv(i) + 1
      if (ibrkx(kk).EQV..true.) cutv(i) = cutv(i) + 1
    else if ((invvec(kk+1).gt.i).AND.(invvec(kk-1).gt.i)) then
      cutv(i) = cutv(i-1) + 2
      if (ibrkx(kk-1).EQV..true.) cutv(i) = cutv(i) - 1
      if (ibrkx(kk).EQV..true.) cutv(i) = cutv(i) - 1
    else
      cutv(i) = cutv(i-1)
      if (ibrkx(kk-1).EQV..true.) then
        if (invvec(setis(i)).gt.i) then
          cutv(i) = cutv(i) - 1
        else
          cutv(i) = cutv(i) + 1
        end if
      end if
      if (ibrkx(kk).EQV..true.) then
        if (invvec(setis(i)+2).gt.i) then
          cutv(i) = cutv(i) - 1
        else
          cutv(i) = cutv(i) + 1
        end if
      end if
    end if
  end do
!
  tmat(:,:) = 0
  tmat(3,3) = cstored+1-ntbrks2
  do i=1,pwidth
    kk = setis(i) + 1
    pvo(2) = 3
    if (invvec(kk+1).gt.i) then
      pvo(3) = 3
    else
      pvo(3) = 2
    end if
    if (invvec(kk-1).gt.i) then
      pvo(1) = 3
    else
      pvo(1) = 2
    end if
    pvn(1:3) = pvo(1:3)
    pvn(2) = 2
    do k=1,2
      if (ibrkx(kk+k-2).EQV..true.) cycle
      tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
      tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
    end do
  end do
  do i=1,cstored
    if ((i.gt.pwidth).AND.((cstored-i).ge.pwidth)) then
      if (((abs(setis(i+pwidth)-setis(i-pwidth)).le.1).AND.(abs(setis(i)-setis(i-pwidth)).le.1)).OR.&
 &        ((abs(setis(i+pwidth)-setis(i-pwidth)).le.1).AND.(abs(setis(i)-setis(i+pwidth)).le.1)).OR.&
 &        ((abs(setis(i)-setis(i+pwidth)).le.1).AND.(abs(setis(i)-setis(i-pwidth)).le.1))) then
!       this should be excessively rare
        j = min(setis(i+pwidth),setis(i-pwidth))
        j = min(j,setis(i))
        pvo(1:5) = 3
        if (invvec(j).ge.i) then
          if (abs(invvec(j)-i).lt.pwidth) pvo(1) = 2
        else
          if (abs(invvec(j)-i).le.pwidth) pvo(1) = 1
        end if
        if (invvec(j+4).ge.i) then
          if (abs(invvec(j+4)-i).lt.pwidth) pvo(5) = 2
        else
          if (abs(invvec(j+4)-i).le.pwidth) pvo(5) = 1
        end if
        pvn(1:5) = pvo(1:5)
        do k=j,j+2
          if (k.eq.setis(i)) then
            pvo(k-j+2) = 2
            pvn(k-j+2) = 1
          else if (k.eq.setis(i-pwidth)) then
            pvo(k-j+2) = 1
            pvn(k-j+2) = 3
          else if (k.eq.setis(i+pwidth)) then
            pvo(k-j+2) = 3
            pvn(k-j+2) = 2
          end if
        end do
        do k=1,4
          if (ibrkx(j+k-1).EQV..true.) cycle
          tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
          tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
        end do
      else if ((abs(setis(i)-setis(i+pwidth)).le.1).OR.(abs(setis(i)-setis(i-pwidth)).le.1).OR.&
 &(abs(setis(i+pwidth)-setis(i-pwidth)).le.1)) then
        if (abs(setis(i+pwidth)-setis(i-pwidth)).le.1) then
          ll = 0
          ii = min(setis(i+pwidth),setis(i-pwidth)) + 1
          kk = max(setis(i+pwidth),setis(i-pwidth)) + 1
          if (setis(i+pwidth).gt.setis(i-pwidth)) then
            pvo(2) = 1
            pvo(3) = 3
          else
            pvo(2) = 3
            pvo(3) = 1
          end if
          pvo(4) = 3
          if (invvec(kk+1).ge.i) then
            if (abs(invvec(kk+1)-i).lt.pwidth) pvo(4) = 2
          else
            if (abs(invvec(kk+1)-i).le.pwidth) pvo(4) = 1
          end if
          pvo(1) = 3
          if (invvec(ii-1).ge.i) then
            if (abs(invvec(ii-1)-i).lt.pwidth) pvo(1) = 2
          else
            if (abs(invvec(ii-1)-i).le.pwidth) pvo(1) = 1
          end if
          pvn(1:4) = pvo(1:4)
          if (setis(i+pwidth).gt.setis(i-pwidth)) then
            pvn(2) = 3
            pvn(3) = 2
          else
            pvn(2) = 2
            pvn(3) = 3
          end if
          do k=1,3
            if (ibrkx(ii+k-2).EQV..true.) cycle
            tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
            tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
          end do
        else if (abs(setis(i)-setis(i-pwidth)).le.1) then
          ll = 1
          ii = min(setis(i),setis(i-pwidth)) + 1
          kk = max(setis(i),setis(i-pwidth)) + 1
          if (setis(i).gt.setis(i-pwidth)) then
            pvo(2) = 1
            pvo(3) = 2
          else
            pvo(2) = 2
            pvo(3) = 1
          end if
          pvo(4) = 3
          if (invvec(kk+1).gt.i) then
            if (abs(invvec(kk+1)-i).lt.pwidth) pvo(4) = 2
          else
            if (abs(invvec(kk+1)-i).le.pwidth) pvo(4) = 1
          end if
          pvo(1) = 3
          if (invvec(ii-1).gt.i) then
            if (abs(invvec(ii-1)-i).lt.pwidth) pvo(1) = 2
          else
            if (abs(invvec(ii-1)-i).le.pwidth) pvo(1) = 1
          end if
          pvn(1:4) = pvo(1:4)
          if (setis(i).gt.setis(i-pwidth)) then
            pvn(2) = 3
            pvn(3) = 1
          else
            pvn(2) = 1
            pvn(3) = 3
          end if
          do k=1,3
            if (ibrkx(ii+k-2).EQV..true.) cycle
            tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
            tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
          end do
        else if (abs(setis(i)-setis(i+pwidth)).le.1) then
          ll = -1
          ii = min(setis(i),setis(i+pwidth)) + 1
          kk = max(setis(i),setis(i+pwidth)) + 1
          if (setis(i).gt.setis(i+pwidth)) then
            pvo(2) = 3
            pvo(3) = 2
          else
            pvo(2) = 2
            pvo(3) = 3
          end if
          pvo(4) = 3
          if (invvec(kk+1).gt.i) then
            if (abs(invvec(kk+1)-i).lt.pwidth) pvo(4) = 2
          else
            if (abs(invvec(kk+1)-i).le.pwidth) pvo(4) = 1
          end if
          pvo(1) = 3
          if (invvec(ii-1).gt.i) then
            if (abs(invvec(ii-1)-i).lt.pwidth) pvo(1) = 2
          else
            if (abs(invvec(ii-1)-i).le.pwidth) pvo(1) = 1
          end if
          pvn(1:4) = pvo(1:4)
          if (setis(i).gt.setis(i+pwidth)) then
            pvn(2) = 2
            pvn(3) = 1
          else
            pvn(2) = 1
            pvn(3) = 2
          end if
          do k=1,3
            if (ibrkx(ii+k-2).EQV..true.) cycle
            tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
            tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
          end do
        end if
        do j=ll,ll
          kk = setis(i+j*pwidth) + 1
          pvo(2) = j+2
          pvo(3) = 3
          if (invvec(kk+1).gt.i) then
            if (abs(invvec(kk+1)-i).lt.pwidth) pvo(3) = 2
          else
            if (abs(invvec(kk+1)-i).le.pwidth) pvo(3) = 1
          end if
          pvo(1) = 3
          if (invvec(kk-1).gt.i) then
            if (abs(invvec(kk-1)-i).lt.pwidth) pvo(1) = 2
          else
            if (abs(invvec(kk-1)-i).le.pwidth) pvo(1) = 1
          end if
          pvn(1:3) = pvo(1:3)
          pvn(2) = j+1
          if (pvn(2).eq.0) pvn(2) = 3
          do k=1,2
            if (ibrkx(kk+k-2).EQV..true.) cycle
            tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
            tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
          end do
        end do
      else
        do j=-1,1
          kk = setis(i+j*pwidth) + 1
          pvo(2) = j+2
          pvo(3) = 3
          if (invvec(kk+1).gt.i) then
            if (abs(invvec(kk+1)-i).lt.pwidth) pvo(3) = 2
          else
            if (abs(invvec(kk+1)-i).le.pwidth) pvo(3) = 1
          end if
          pvo(1) = 3
          if (invvec(kk-1).gt.i) then
            if (abs(invvec(kk-1)-i).lt.pwidth) pvo(1) = 2
          else
            if (abs(invvec(kk-1)-i).le.pwidth) pvo(1) = 1
          end if
          pvn(1:3) = pvo(1:3)
          pvn(2) = j+1
          if (pvn(2).eq.0) pvn(2) = 3
          do k=1,2
            if (ibrkx(kk+k-2).EQV..true.) cycle
            tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
            tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
          end do
        end do
      end if
    else if (i.le.pwidth) then
      if (abs(setis(i)-setis(i+pwidth)).le.1) then
        ii = min(setis(i),setis(i+pwidth)) + 1
        kk = max(setis(i),setis(i+pwidth)) + 1
        if (setis(i).gt.setis(i+pwidth)) then
          pvo(2) = 3
          pvo(3) = 2
        else
          pvo(2) = 2
          pvo(3) = 3
        end if
        pvo(4) = 3
        if (invvec(kk+1).gt.i) then
          if (abs(invvec(kk+1)-i).lt.pwidth) pvo(4) = 2
        else
          if (abs(invvec(kk+1)-i).le.pwidth) pvo(4) = 1
        end if
        pvo(1) = 3
        if (invvec(ii-1).gt.i) then
          if (abs(invvec(ii-1)-i).lt.pwidth) pvo(1) = 2
        else
          if (abs(invvec(ii-1)-i).le.pwidth) pvo(1) = 1
        end if
        pvn(1:4) = pvo(1:4)
        if (setis(i).gt.setis(i+pwidth)) then
          pvn(2) = 2
          pvn(3) = 1
        else
          pvn(2) = 1
          pvn(3) = 2
        end if
        do k=1,3
          if (ibrkx(ii+k-2).EQV..true.) cycle
          tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
          tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
        end do
      else
        do j=0,1
          kk = setis(i+j*pwidth) + 1
          pvo(2) = j+2
          pvo(3) = 3
          if (invvec(kk+1).gt.i) then
            if (abs(invvec(kk+1)-i).lt.pwidth) pvo(3) = 2
          else
            if (abs(invvec(kk+1)-i).le.pwidth) pvo(3) = 1
          end if
          pvo(1) = 3
          if (invvec(kk-1).gt.i) then
            if (abs(invvec(kk-1)-i).lt.pwidth) pvo(1) = 2
          else
            if (abs(invvec(kk-1)-i).le.pwidth) pvo(1) = 1
          end if
          pvn(1:3) = pvo(1:3)
          pvn(2) = j+1
          if (pvn(2).eq.0) pvn(2) = 3
          do k=1,2
            if (ibrkx(kk+k-2).EQV..true.) cycle
            tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
            tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
          end do
        end do
      end if
    else
      if (abs(setis(i)-setis(i-pwidth)).le.1) then
        ii = min(setis(i),setis(i-pwidth)) + 1
        kk = max(setis(i),setis(i-pwidth)) + 1
        if (setis(i).gt.setis(i-pwidth)) then
          pvo(2) = 1
          pvo(3) = 2
        else
          pvo(2) = 2
          pvo(3) = 1
        end if
        pvo(4) = 3
        if (invvec(kk+1).gt.i) then
          if (abs(invvec(kk+1)-i).lt.pwidth) pvo(4) = 2
        else
          if (abs(invvec(kk+1)-i).le.pwidth) pvo(4) = 1
        end if
        pvo(1) = 3
        if (invvec(ii-1).gt.i) then
          if (abs(invvec(ii-1)-i).lt.pwidth) pvo(1) = 2
        else
          if (abs(invvec(ii-1)-i).le.pwidth) pvo(1) = 1
        end if
        pvn(1:4) = pvo(1:4)
        if (setis(i).gt.setis(i-pwidth)) then
          pvn(2) = 3
          pvn(3) = 1
        else
          pvn(2) = 1
          pvn(3) = 3
        end if
        do k=1,3
          if (ibrkx(ii+k-2).EQV..true.) cycle
          tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
          tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
        end do
      else
        do j=-1,0
          kk = setis(i+j*pwidth) + 1
          pvo(2) = j+2
          pvo(3) = 3
          if (invvec(kk+1).gt.i) then
            if (abs(invvec(kk+1)-i).lt.pwidth) pvo(3) = 2
          else
            if (abs(invvec(kk+1)-i).le.pwidth) pvo(3) = 1
          end if
          pvo(1) = 3
          if (invvec(kk-1).gt.i) then
            if (abs(invvec(kk-1)-i).lt.pwidth) pvo(1) = 2
          else
            if (abs(invvec(kk-1)-i).le.pwidth) pvo(1) = 1
          end if
          pvn(1:3) = pvo(1:3)
          pvn(2) = j+1
          if (pvn(2).eq.0) pvn(2) = 3
          do k=1,2
            if (ibrkx(kk+k-2).EQV..true.) cycle
            tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
            tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
          end do
        end do
      end if
    end if
    ii = min(i,pwidth)
    jj = min(cstored-i,pwidth)
    kk = max(0,cstored - ii - jj)
    vv1(:) = tmat(1,:)
    vv2(:) = tmat(2,:)
    vv3(:) = tmat(3,:)
    write(iu,666) i,cstored-i,setis(i),cutv(i),distv(i),ivec2(i),invvec(ivec2(i)+1),cutv(invvec(ivec2(i)+1)),vv1,vv2,vv3,ii,jj,kk
  end do
 666 format(4(i10,1x),g12.5,1x,1000(i10,1x))
!
  deallocate(cutv)
!
end
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
