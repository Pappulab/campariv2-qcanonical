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
!------------------------------------------------------------------------
!
!
!             #################################
!             #                               #
!             # FIRST: WRAPPER ROUTINES       #
!             #        FOR GLOBAL ENERGIES    #
!             #                               #
!             #################################
!
!
!------------------------------------------------------------------------
!
! "energy" a simple service routine to calculate the total
! potential energy
! uses the residue-based function directly
! the commented block allow testing of the vectorized fxns against
! non-vectorized ones
!
function energy(evec)
!
  use sequen
  use energies
  use molecule
  use atoms
  ! Martin : debug only
  use iounit
  
#ifdef ENABLE_THREADS
  use threads
#endif
  use wl
  
  use movesets ! for hsq
!
  implicit none
!
  integer i,j,imol,rs,sta(2),sto(2),incr,tpi,azero
  RTYPE evec(MAXENERGYTERMS),energy
  logical sayno
#ifdef ENABLE_THREADS
  integer tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
  maxtpn = 0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
  evec(:) = 0.0
  sayno = .false.
  azero = 0
!
  if (ideal_run.EQV..false.) then
    call clear_rsp()
    if (use_IMPSOLV.EQV..true.) call init_svte(3)
  end if
  
!    print *,"here1"
!  print *,use_BOND(:)
!  flush(6)
  
!  call CPU_time(t1)
!  do i=1,nseq
!    do j=i,nseq
!      call Ven_rsp(evec,i,j,sayno)
!    end do
!  end do
!  call CPU_time(t2)
!  write(*,*) 'SR-novec: ',t2-t1,evec(1)+evec(3),sum(svte(:))
!  evec(1) = 0.0
!  evec(2) = 0.0
!  evec(3) = 0.0
!  if (use_IMPSOLV.EQV..true.) call init_svte(3)
!  call CPU_time(t1)
#ifdef ENABLE_THREADS
  call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,incr,tpn,tpi,i,j,rs,imol)
  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  sto(1) = nseq
  sto(2) = nmol
  sta(1:2) = tpi
  incr = tpn
#else
  sta(1:2) = 1
  sto(1) = nseq
  sto(2) = nmol
  tpi = 1
  incr = 1
#endif
  do i=sta(1),sto(1),incr !1,nseq
    do j=i,nseq
      call Ven_rsp(evec_thr(:,tpi),i,j,sayno)
    end do
  end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif 
!  call CPU_time(t2)
!  write(*,*) 'SR-Vect.: ',t2-t1,evec(1)+evec(3),sum(svte(:))
!
!  call CPU_time(t1)
  if (use_IMPSOLV.EQV..true.) then
    call init_svte(6)

    if (use_POLAR.EQV..true.) then
      call setup_scrqs2()
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP BARRIER
#endif
!  call CPU_time(t2)
!  write(*,*) 'setup 1: ',t2-t1
! 
!  call CPU_time(t1)
!  do i=1,nseq
!    do j=i,nseq
!      call en_rsp_long(evec,i,j,sayno)
!    end do
!  end do
!  call CPU_time(t2)
!  write(*,*) 'LR-novec: ',t2-t1,evec(6)+evec(9)
!  evec(6) = 0.0
!  evec(9) = 0.0
!  call CPU_time(t1)
  do i=sta(1),sto(1),incr !1,nseq
      
    do j=i,nseq
      call Ven_rsp_long(evec_thr(:,tpi),i,j,sayno)
    end do
  end do
!  call CPU_time(t2)
!  write(*,*) 'LR-Vect.: ',t2-t1,evec(6)+evec(9)
!  call CPU_time(t1)
  if (use_IMPSOLV.EQV..true.) then
    do i=sta(1),sto(1),incr !1,nseq
      call en_freesolv(evec_thr(:,tpi),i)
    end do
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    call en_totsav(evec_thr(:,tpi))
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
  end if
!  call CPU_time(t2)
!  write(*,*) 'FOS: ',t2-t1
!
!  call CPU_time(t1)
  do rs=sta(1),sto(1),incr !1,nseq
    if (use_BOND(1).EQV..true.) then
      call en_bonds(rs,evec_thr(:,tpi))
    end if
    if (use_BOND(2).EQV..true.) then
      call en_angles(rs,evec_thr(:,tpi))
    end if
    if (use_BOND(3).EQV..true.) then
      call en_impropers(rs,evec_thr(:,tpi))
    end if
    if (use_BOND(4).EQV..true.) then
      call en_torsions(rs,evec_thr(:,tpi))
    end if
    if (use_BOND(5).EQV..true.) then
      call en_cmap(rs,evec_thr(:,tpi))
    end if
  end do

!  call CPU_time(t2)
!  write(*,*) 'bonded: ',t2-t1
!
!  call CPU_time(t1)
  if (use_TOR.EQV..true.) then
    do i=sta(1),sto(1),incr !1,nseq
      if (par_TOR2(i).gt.0) then
        call en_torrs(evec_thr(:,tpi),i)
      end if
    end do
  end if
!
  do rs=sta(1),sto(1),incr !1,nseq
    call e_boundary_rs(rs,evec_thr(:,tpi),1)
  end do
!
  if (use_ZSEC.EQV..true.) then
    do imol=sta(2),sto(2),incr !1,nmol
      call en_zsec_gl(imol,evec_thr(:,tpi))
    end do
  end if
!
  if (use_POLY.EQV..true.) then
    do imol=sta(2),sto(2),incr !1,nmol
      call en_poly_gl(imol,evec_thr(:,tpi),0)
    end do
  end if
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  if (use_DSSP.EQV..true.) then
    call en_dssp_gl(evec_thr(:,tpi))
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  if (use_LCTOR.EQV..true.) then
    call en_lc_tor(evec_thr(:,tpi))
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  if (use_EMICRO.EQV..true.) then
    call en_emicro_gl(evec_thr(:,tpi),azero)
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  if (use_DREST.EQV..true.) then
    call edrest(evec_thr(:,tpi))
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP BARRIER
!$OMP SINGLE
  if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
!
  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif
!
  if (do_accelsim.EQV..true.) then
    call els_manage_justE(evec,1)
    evec(:) = evec(:) + hmjam%boosts(:,1)
  end if
  
  if ((do_pka.eqv..true.).or.(do_hsq.eqv..true.)) then 
    evec(13)=BATH_ENER
  else if (FEG_FOS_OFFSET.eqv..true.) then 
    evec(13)=(scale_FEGS(4)**(FEG_FOS_OFFSET_FACT))*FEG_FOS_OFFSET_VAL 
  end if 
  energy = sum(evec(:))
!
  return
!
end
!
!-----------------------------------------------------------------------
!
! same as energy, but uses all cutoffs (simple) turned on by the user
! we have to be careful, however, about clearing out the current rsp_mat(,)
! it is not advisable to leave 
!
function energy3(evec,nbup)
!
  use sequen
  use energies
  use molecule
  use cutoffs
  use system
  use atoms
  use energies
  use martin_own
#ifdef ENABLE_THREADS
  use threads
#endif
  use wl
  
  use movesets ! for hsq
!
  implicit none

  integer i,j,imol,rs,sta(2),sto(2),incr,tpi,azero
  RTYPE energy3,evec(MAXENERGYTERMS)
  logical nbup
#ifdef ENABLE_THREADS
  integer tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
  maxtpn = 0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
  energy3 = 0.0
  evec(:) = 0.0
  azero = 0

  if (ideal_run.EQV..false.) then
!
    if (use_IMPSOLV.EQV..true.) call init_svte(3)
!
    if (use_cutoffs.EQV..true.) then
      if (nbup.EQV..true.) then
        call clear_rsp()
        if (use_mcgrid.EQV..true.) then
          do rs=1,nseq
            call grd_respairs(rs)
          end do
        else if (use_rescrit.EQV..true.) then
          do rs=1,nseq
            call respairs(rs)
          end do
        end if
      end if
    end if
  end if
  
#ifdef ENABLE_THREADS
  call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,incr,tpn,tpi,i,j,rs,imol)
  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  sto(1) = nseq
  sto(2) = nmol
  sta(1:2) = tpi
  incr = tpn
#else
  sta(1:2) = 1
  sto(1) = nseq
  sto(2) = nmol
  tpi = 1
  incr = 1
#endif
  if (ideal_run.EQV..false.) then
    do i=sta(1),sto(1),incr !1,nseq
      do j=i,nseq
        if ((rsp_mat(j,i).eq.0).AND.(use_cutoffs)) cycle
        call Ven_rsp(evec_thr(:,tpi),i,j,use_cutoffs)
      end do
    end do
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
    if (use_IMPSOLV.EQV..true.) then
      call init_svte(6)
      if (use_POLAR.EQV..true.) then
        call setup_scrqs2()
      end if
    end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP BARRIER
#endif


    do i=sta(1),sto(1),incr !1,nseq
      do j=i,nseq
        if ((rsp_mat(i,j).eq.0).AND.(use_cutoffs)) cycle
        call Ven_rsp_long(evec_thr(:,tpi),i,j,use_cutoffs)
      end do
    end do
!
    if (use_IMPSOLV.EQV..true.) then
      !en_fos_memory(:,:)=0.0!martin ; need to reset at every frame
       
      !sav_memory(:,:,:)=0.0
      do i=sta(1),sto(1),incr !1,nseq ! martin goes through all the residues 
        call en_freesolv(evec_thr(:,tpi),i)
      end do
      
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
      call en_totsav(evec_thr(:,tpi))
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
    end if
  end if ! if ideal_run
!
  do rs=sta(1),sto(1),incr !1,nseq
      
    if (use_BOND(1).EQV..true.) then
      call en_bonds(rs,evec_thr(:,tpi))
    end if
    if (use_BOND(2).EQV..true.) then
      call en_angles(rs,evec_thr(:,tpi))
    end if
    if (use_BOND(3).EQV..true.) then
      call en_impropers(rs,evec_thr(:,tpi))
    end if
    if (use_BOND(4).EQV..true.) then
      call en_torsions(rs,evec_thr(:,tpi))
    end if
    if (use_BOND(5).EQV..true.) then
      call en_cmap(rs,evec_thr(:,tpi))
    end if
  end do
!
  if (use_TOR.EQV..true.) then
    do i=sta(1),sto(1),incr !1,nseq
      if (par_TOR2(i).gt.0) then
        call en_torrs(evec_thr(:,tpi),i)
      end if
    end do
  end if
!
  if ((bnd_type.ne.1).OR.(bnd_shape.ne.1)) then
    do rs=sta(1),sto(1),incr !1,nseq
      call e_boundary_rs(rs,evec_thr(:,tpi),1)
    end do
  end if
!
  if (use_ZSEC.EQV..true.) then
    do imol=sta(2),sto(2),incr !1,nmol
      call en_zsec_gl(imol,evec_thr(:,tpi))
    end do
  end if
!
  if (use_POLY.EQV..true.) then
    do imol=sta(2),sto(2),incr !1,nmol
      call en_poly_gl(imol,evec_thr(:,tpi),0)
    end do
  end if
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  if (use_DSSP.EQV..true.) then
    call en_dssp_gl(evec_thr(:,tpi))
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  if (use_LCTOR.EQV..true.) then
    call en_lc_tor(evec_thr(:,tpi))
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  if (use_EMICRO.EQV..true.) then
    call en_emicro_gl(evec_thr(:,tpi),azero)
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  if (use_DREST.EQV..true.) then
    call edrest(evec_thr(:,tpi))
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP BARRIER
!$OMP SINGLE
  if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
!
  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif
!
  if (do_accelsim.EQV..true.) then
    call els_manage_justE(evec,1)
    evec(:) = evec(:) + hmjam%boosts(:,1)
  end if

  if ((do_pka.eqv..true.).or.(do_hsq.eqv..true.)) then 
    evec(13)=BATH_ENER
  else if (FEG_FOS_OFFSET.eqv..true.) then 
    evec(13)=(scale_FEGS(4)**(10/12.))*FEG_FOS_OFFSET_VAL

  end if 
  
  energy3 = sum(evec(:))
!
  return
!
end
!
!------------------------------------------------------------------------
!
!
!             #################################
!             #                               #
!             # SECOND: WRAPPER ROUTINES      #
!             #   FOR INCREMENTAL ENERGIES    #
!             #                               #
!             #################################
!
!
! a simple wrapper routine that handles short-range energy calculations for sidechain
! moves (only relevant terms are computed ...)
!
subroutine chi_energy_short(rs,evec,cut,mode)
!
  use sequen
  use iounit
  use energies
  use cutoffs
  use molecule
  use fyoc
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rs,j,mode,imol,sto,sta,tpi,incr
  RTYPE evec(MAXENERGYTERMS)
  logical cut
#ifdef ENABLE_THREADS
  integer cvl,tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
! the boundary terms depend only on coordinates -> can be taken care of immediately
  call e_boundary_rs(rs,evec,mode)
!
! bonded terms may change but are restricted to a truly local effect (note -1/+1-padding, though)
  imol = molofrs(rs)
#ifdef ENABLE_THREADS
  cvl = min(rsmol(imol,2),rs+1) - max(rsmol(imol,1),rs-1) + 1
  call omp_set_num_threads(min(cvl,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,incr)
  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  sto = min(rsmol(imol,2),rs+1)
  sta = tpi - 1 + max(rsmol(imol,1),rs-1)
  incr = tpn
#else
  sta = max(rsmol(imol,1),rs-1)
  sto = min(rsmol(imol,2),rs+1)
  tpi = 1
  incr = 1
#endif
  do j=sta,sto,incr
    if (use_BOND(1).EQV..true.) call en_bonds(j,evec_thr(:,tpi))
    if (use_BOND(2).EQV..true.) call en_angles(j,evec_thr(:,tpi))
    if (use_BOND(3).EQV..true.) call en_impropers(j,evec_thr(:,tpi))
    if (use_BOND(4).EQV..true.) call en_torsions(j,evec_thr(:,tpi))
    if (use_BOND(5).EQV..true.) call en_cmap(j,evec_thr(:,tpi))
    if (disulf(j).gt.0) then
      if ((disulf(j).lt.sta).OR.(disulf(j).gt.sto)) then
        if (use_BOND(1).EQV..true.) call en_bonds(disulf(j),evec_thr(:,tpi))
        if (use_BOND(2).EQV..true.) call en_angles(disulf(j),evec_thr(:,tpi))
        if (use_BOND(3).EQV..true.) call en_impropers(disulf(j),evec_thr(:,tpi))
        if (use_BOND(4).EQV..true.) call en_torsions(disulf(j),evec_thr(:,tpi))
        if (use_BOND(5).EQV..true.) call en_cmap(disulf(j),evec_thr(:,tpi))
      end if
    end if
  end do
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
  if (tpn.gt.maxtpn) maxtpn = tpn
#endif
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     zero out svte and svbu (temporary arrays for changes in SAV)
      call init_svte(0)
    else if (mode.eq.1) then
!     store svte (the changes induced by the sampling) in svbu, zero out svte again
      call init_svte(1)
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP BARRIER
!$OMP END PARALLEL
#endif
! let's only do all this setup work if necessary
  if ((use_IMPSOLV.EQV..true.).OR.(use_WCA.EQV..true.).OR.(use_IPP.EQV..true.).OR.&
 &    (use_attLJ.EQV..true.).OR.(use_CORR.EQV..true.).OR.(use_FEGS(1).EQV..true.).OR.&
 &    (use_FEGS(3).EQV..true.)) then
!
    if (cut.EQV..true.) then
      if (use_mcgrid.EQV..true.) then
        call grd_respairs(rs)
      else if (use_rescrit.EQV..true.) then
        call respairs(rs)
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in chi_energy_shor&
 &t(...).'
        call fexit()
      end if
    end if
#ifdef ENABLE_THREADS
    call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,incr)
    tpn = omp_get_num_threads()
    incr = tpn
    tpi = omp_get_thread_num() + 1
    sto = nseq
    sta = tpi
#else
    sta = 1
    sto = nseq
    incr = 1
#endif

    if (cut.EQV..true.) then
      do j=sta,sto,incr ! 1,nseq
        if (rsp_mat(max(rs,j),min(rs,j)).ne.0) then
          call Ven_rsp(evec_thr(:,tpi),min(rs,j),max(rs,j),use_cutoffs)
        end if
      end do
    else
      do j=sta,sto,incr ! 1,nseq
        call Ven_rsp(evec_thr(:,tpi),rs,j,cut)
      end do
    end if
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
    if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    if (cut.EQV..false.) then
      call clear_rsp2(rs)
    end if
!
  else
!   do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif 
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.1) then
!     get the difference in SAV
      call init_svte(2)
!     now use that difference to determine which dependent interactions we have to re-compute
    end if
  end if
!
end
!c
!--------------------------------------------------------------------------
!
! a simple wrapper routine that handles long-range energy calculations for sidechain
! moves (only relevant terms are computed ...)
!
subroutine chi_energy_long(rs,evec,cut,mode)
!
  use iounit
  use sequen
  use energies
  use cutoffs
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rs,j,mode,i,k,l,nrsl,rsl(nseq),tpi,sta(2),sto(2),incr
  RTYPE evec(MAXENERGYTERMS)
  logical cut,inrsl(nseq)
#ifdef ENABLE_THREADS
  integer tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
  maxtpn = 0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    call setup_scrqs2()
  end if
!
  nrsl = 0
!     
  if (use_TOR.EQV..true.) then
    if (par_TOR2(rs).gt.0) then
      call en_torrs(evec,rs)
    end if
  end if
! 
  if ((use_IMPSOLV.EQV..true.).AND.((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.))) then
!
!   the problematic part is that the interaction of two residues can
!   indeed be changed by a third residue, which is changing conformation.
!   it is henceforth not safe to compute just the terms that change
!   by virtue of the conformational change.
!   to detect the implicitly affected terms we're using the temporary SAV
!   array (i.e., any atom whose SAV was affected by the conformational
!   change potentially interacts differently with all other atoms now)
!   through the rs_vec-array which is created during the last init_svte-call
!   in the short-range energy computation
!   note, however, that just like for every other interaction, the residue
!   based cutoff criterion still has to be fulfilled

    if (cut.EQV..true.) then
      inrsl(:) = .false.
      do j=1,rs-1
       if (rs_vec(j).eq.1) then
          nrsl = nrsl + 1
          rsl(nrsl) = j
          inrsl(j) = .true.
        end if
      end do
      do j=rs,rs
        nrsl = nrsl + 1
        rsl(nrsl) = j
        inrsl(j) = .true.
      end do
      do j=rs+1,nseq
        if (rs_vec(j).eq.1) then
          nrsl = nrsl + 1
          rsl(nrsl) = j
          inrsl(j) = .true.
        end if
      end do
      if (use_mcgrid.EQV..true.) then
        do i=1,nrsl
          call grd_respairs(rsl(i))
        end do
      else if (use_rescrit.EQV..true.) then
        do i=1,nrsl
          call respairs(rsl(i))
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in chi_energy_lo&
 &ng(...).'
        call fexit()
      end if
!
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,i,j,k,l,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
      sto(2) = nrsl
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      sto(2) = nrsl
      incr = 1
#endif
      do i=1,nrsl
        j = rsl(i)
        do k=sta(1),sto(1),incr ! 1,nseq
!         we avoid double counting through the use of the inverse map inrsl
          if (inrsl(k).EQV..true.) cycle
!
          if (((j.ge.k).AND.(rsp_mat(k,j).ne.0)).OR.&
 &            ((k.gt.j).AND.(rsp_mat(j,k).ne.0))) then
            call Ven_rsp_long(evec_thr(:,tpi),j,k,use_cutoffs)
          end if
        end do
!       and take care of the intra-rsl terms here
        sta(2) = i + tpi - 1
        do l=sta(2),sto(2),incr!i,nrsl
          k = rsl(l)
          if (((j.ge.k).AND.(rsp_mat(k,j).ne.0)).OR.&
 &            ((k.gt.j).AND.(rsp_mat(j,k).ne.0))) then
            call Ven_rsp_long(evec_thr(:,tpi),j,k,use_cutoffs)
          end if
        end do
      end do
!
      do i=sta(1),sto(1),incr ! 1,nseq
        if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
          call en_freesolv(evec_thr(:,tpi),i)
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
#endif
      do i=1,nrsl
        call clear_rsp2(rsl(i))
      end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,i,j,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
      do i=sta(1),sto(1),incr !1,nseq
        do j=i,nseq
          call Ven_rsp_long(evec_thr(:,tpi),i,j,cut)
        end do
        call en_freesolv(evec_thr(:,tpi),i)
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
!
  else if (use_IMPSOLV.EQV..true.) then
#ifdef ENABLE_THREADS
    call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,i,j,incr)
    tpn = omp_get_num_threads()
    incr = tpn
    tpi = omp_get_thread_num() + 1
    sto(1) = nseq
    sta(1) = tpi
#else
    tpi = 1
    sta(1) = 1
    sto(1) = nseq
    incr = 1
#endif
    do i=sta(1),sto(1),incr !1,nseq
      call en_freesolv(evec_thr(:,tpi),i)
    end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
!
  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.)) then
    if (cut.EQV..true.) then
      if (use_mcgrid.EQV..true.) then
        call grd_respairs(rs)
      else if (use_rescrit.EQV..true.) then
        call respairs(rs)
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in chi_energy_lo&
 &ng(...).'
        call fexit()
      end if
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
      do j=sta(1),sto(1),incr
        if (rsp_mat(min(j,rs),max(rs,j)).ne.0) then
          call Ven_rsp_long(evec_thr(:,tpi),min(j,rs),max(rs,j),cut)
         end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
      call clear_rsp2(rs)
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
      do j=sta(1),sto(1),incr !1,nseq
        call Ven_rsp_long(evec_thr(:,tpi),rs,j,cut)
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     now update the atsav-array using svte (to be ready for posterior long-range calculation)
      call init_svte(5)
    else if (mode.eq.1) then
!     clear out the residue-svte-monitoring-vector
      do i=1,nseq
        rs_vec(i) = 0
        rs_vec2(i) = 0
      end do
    end if
  end if
!
end
!
!------------------------------------------------------------------------
!
! a simple wrapper routine that handles short-range energy calculations for pivot
! moves (only relevant terms are computed ...)

subroutine pivot_energy_short(rs,evec,cut,mode,ct)
!
  use iounit
  use sequen
  use cutoffs
  use energies
  use molecule
  use fyoc
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rs,j,mode,i,k,irs,frs,sta(2),sto(2),incr,tpi,imol
  logical ct
  RTYPE evec(MAXENERGYTERMS)
  logical cut
#ifdef ENABLE_THREADS
  integer tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
  maxtpn = 0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
! set residue range (lever arm) based on whether N- or C-aligned (the moving parts change!)
  if (ct.EQV..true.) then
    irs = rsmol(molofrs(rs),1)
    frs = rs
  else
    irs = rs
    frs = rsmol(molofrs(rs),2)
  end if
!
! the boundary terms depend only on coordinates -> can be taken care of immediately
#ifdef ENABLE_THREADS
  call omp_set_num_threads(min(frs-irs+1,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,incr)
  
  

  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  sto(1) = frs
  sta(1) = irs + tpi - 1
  sta(2) = max(rsmol(molofrs(rs),1),irs-1) + tpi - 1
  sto(2) = min(rsmol(molofrs(rs),2),frs+1)
  incr = tpn
#else
  sta(1) = irs
  sto(1) = frs
  sta(2) = max(rsmol(molofrs(rs),1),irs-1)
  sto(2) = min(rsmol(molofrs(rs),2),frs+1)
  incr = 1
  tpi = 1
#endif
  do j=sta(1),sto(1),incr
    call e_boundary_rs(j,evec_thr(:,tpi),mode)
  end do

  if (n_crosslinks.gt.0) then
    do j=sta(2),sto(2),incr
      if (disulf(j).gt.0) then
        if (abs(rs-j).le.1) then
          if (use_BOND(1).EQV..true.) call en_bonds(disulf(j),evec_thr(:,tpi))
          if (use_BOND(2).EQV..true.) call en_angles(disulf(j),evec_thr(:,tpi))
          if (use_BOND(3).EQV..true.) call en_impropers(disulf(j),evec_thr(:,tpi))
          if (use_BOND(4).EQV..true.) call en_torsions(disulf(j),evec_thr(:,tpi))
          if (use_BOND(5).EQV..true.) call en_cmap(disulf(j),evec_thr(:,tpi))
        else if ((disulf(j).gt.sto(2)).OR.(disulf(j).lt.sta(2))) then
          if (use_BOND(1).EQV..true.) call en_bonds(j,evec_thr(:,tpi))
          if (use_BOND(2).EQV..true.) call en_angles(j,evec_thr(:,tpi))
          if (use_BOND(3).EQV..true.) call en_impropers(j,evec_thr(:,tpi))
          if (use_BOND(4).EQV..true.) call en_torsions(j,evec_thr(:,tpi))
          if (use_BOND(5).EQV..true.) call en_cmap(j,evec_thr(:,tpi))
          if (use_BOND(1).EQV..true.) call en_bonds(disulf(j),evec_thr(:,tpi))
          if (use_BOND(2).EQV..true.) call en_angles(disulf(j),evec_thr(:,tpi))
          if (use_BOND(3).EQV..true.) call en_impropers(disulf(j),evec_thr(:,tpi))
          if (use_BOND(4).EQV..true.) call en_torsions(disulf(j),evec_thr(:,tpi))
          if (use_BOND(5).EQV..true.) call en_cmap(disulf(j),evec_thr(:,tpi))
        end if
      end if
    end do
  end if
#ifdef ENABLE_THREADS
!$OMP SINGLE
  if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP BARRIER
!$OMP END PARALLEL
#endif
!
! bonded terms may change but are restricted to a truly local effect (note -1/+1-padding, though)
  imol = molofrs(rs)
#ifdef ENABLE_THREADS
  call omp_set_num_threads(min(min(rsmol(imol,2),rs+1)-max(rsmol(imol,1),rs-1)+1,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,incr) FIRSTPRIVATE(imol)
  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  sto(1) = min(rsmol(imol,2),rs+1)
  sta(1) = tpi - 1 + max(rsmol(imol,1),rs-1)
  incr = tpn
#else
  sta(1) = max(rsmol(imol,1),rs-1)
  sto(1) = min(rsmol(imol,2),rs+1)
  tpi = 1
  incr = 1
#endif

  do j=sta(1),sto(1),incr
    if (use_BOND(1).EQV..true.) call en_bonds(j,evec_thr(:,tpi))
    if (use_BOND(2).EQV..true.) call en_angles(j,evec_thr(:,tpi))
    if (use_BOND(3).EQV..true.) call en_impropers(j,evec_thr(:,tpi))
    if (use_BOND(4).EQV..true.) call en_torsions(j,evec_thr(:,tpi))
    if (use_BOND(5).EQV..true.) call en_cmap(j,evec_thr(:,tpi))
  end do

#ifdef ENABLE_THREADS
!$OMP SINGLE
  if (tpn.gt.maxtpn) maxtpn = tpn
#endif
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     zero out svte and svbu (temporary arrays for changes in SAV)
      call init_svte(0)
    else if (mode.eq.1) then
!     store svte (the changes induced by the sampling) in svbu, zero out svte again
      call init_svte(1)
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP BARRIER
!$OMP END PARALLEL
#endif

! let's only do all this setup work if necessary
  if ((use_IMPSOLV.EQV..true.).OR.(use_WCA.EQV..true.).OR.(use_IPP.EQV..true.).OR.&
 &    (use_attLJ.EQV..true.).OR.(use_CORR.EQV..true.).OR.(use_FEGS(1).EQV..true.).OR.&
 &    (use_FEGS(3).EQV..true.)) then
!
    if (cut.EQV..true.) then
      if (use_mcgrid.EQV..true.) then
        do i=irs,frs
          call grd_respairs(i)
        end do
      else if (use_rescrit.EQV..true.) then
        do i=irs,frs
          call respairs(i)
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in pivot_energy_&
 &long(...).'
        call fexit()
      end if
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
      do j=irs,frs
        do k=sta(1),sto(1),incr !1,irs-1 and frs+1,nseq
          if ((k.ge.irs).AND.(k.le.frs)) cycle 
          if (rsp_mat(max(j,k),min(j,k)).ne.0) then
            call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
      end do
      sto(2) = frs
      sta(2) = irs + tpi - 1
      do j=sta(2),sto(2),incr
        if (rsp_mat(max(j,rs),min(j,rs)).ne.0) then
          call Ven_rsp(evec_thr(:,tpi),min(rs,j),max(rs,j),use_cutoffs)
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
      do i=irs,frs
        call clear_rsp2(i)
      end do
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
!      do j=sta(1),sto(1),incr !1,nseq
      do j=irs,frs
        do k=sta(1),sto(1),incr !1,irs-1 and frs+1,nseq
          if ((k.ge.irs).AND.(k.le.frs)) cycle
            call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),cut)
        end do
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
        call Ven_rsp(evec_thr(:,tpi),rs,j,cut)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
!
  else
!   do nothing
  end if


  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif

  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.1) then
      call init_svte(2)
    end if
  end if

  end
!
!--------------------------------------------------------------------------
!
! a simple wrapper routine that handles long-range energy calculations for pivot
! moves (only relevant terms are computed ...)
!
subroutine pivot_energy_long(rs,evec,cut,mode,ct)
!
  use iounit
  use sequen
  use cutoffs
  use energies
  use molecule
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rs,j,mode,i,k,l,rsl(nseq),nrsl,irs,frs,sta(2),sto(2),incr,tpi
  logical ct
  RTYPE evec(MAXENERGYTERMS)
  logical cut,inrsl(nseq)
#ifdef ENABLE_THREADS
  integer tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
  maxtpn = 0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    call setup_scrqs2()
  end if
!
  nrsl = 0
!     
  if (use_TOR.EQV..true.) then
    if (par_TOR2(rs).gt.0) then
      call en_torrs(evec,rs)
    end if
  end if
!
! set residue range (lever arm) based on whether N- or C-aligned (the moving parts change!)
  if (ct.EQV..true.) then
    irs = rsmol(molofrs(rs),1)
    frs = rs
  else
    irs = rs
    frs = rsmol(molofrs(rs),2)
  end if
!
  if ((use_IMPSOLV.EQV..true.).AND.((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.))) then
!
!   the problematic part is that the interaction of two residues can
!   indeed be changed by a third residue, which is changing conformation.
!   it is henceforth not safe to compute just the terms that change
!   by virtue of the conformational change.
!   to detect the implicitly affected terms we're using the temporary SAV
!   array (i.e., any atom whose SAV was affected by the conformational
!   change potentially interacts differently with all other atoms now)
!   through the rs_vec-array which is created during the last init_svte-call
!   in the short-range energy computation
!   note, however, that just like for every other interaction, the residue
!   based cutoff criterion still has to be fulfilled
!
    if (cut.EQV..true.) then
      inrsl(:) = .false.
      do j=1,irs-1
       if (rs_vec(j).eq.1) then
          nrsl = nrsl + 1
          rsl(nrsl) = j
          inrsl(j) = .true.
        end if
      end do
      do j=irs,frs
        nrsl = nrsl + 1
        rsl(nrsl) = j
        inrsl(j) = .true.
      end do
      do j=frs+1,nseq
        if (rs_vec(j).eq.1) then
          nrsl = nrsl + 1
          rsl(nrsl) = j
          inrsl(j) = .true.
        end if
      end do
      if (use_mcgrid.EQV..true.) then
        do i=1,nrsl
          call grd_respairs(rsl(i))
        end do
      else if (use_rescrit.EQV..true.) then
        do i=1,nrsl
          call respairs(rsl(i))
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in pivot_energy_&
 &long(...).'
        call fexit()
      end if
!
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,i,j,k,l,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
      sto(2) = nrsl
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      sto(2) = nrsl
      incr = 1
#endif
      do i=1,nrsl
        j = rsl(i)
        do k=sta(1),sto(1),incr !1,nseq
!         to avoid double-counting, we'll utilize the inverse map
          if (inrsl(k).EQV..true.) cycle
!
          if (((j.ge.k).AND.(rsp_mat(k,j).ne.0)).OR.&
 &            ((k.gt.j).AND.(rsp_mat(j,k).ne.0))) then
            call Ven_rsp_long(evec_thr(:,tpi),j,k,use_cutoffs)
          end if
        end do
        sta(2) = i + tpi - 1
!       and take care of the rest here
        do l=sta(2),sto(2),incr !do l=i,nrsl
          k = rsl(l)
          if (((j.ge.k).AND.(rsp_mat(k,j).ne.0)).OR.&
 &            ((k.gt.j).AND.(rsp_mat(j,k).ne.0))) then
            call Ven_rsp_long(evec_thr(:,tpi),j,k,use_cutoffs)
          end if
        end do
      end do
!
      do i=sta(1),sto(1),incr ! 1,nseq
        if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
          call en_freesolv(evec_thr(:,tpi),i)
        end if
      end do
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
#endif
      do i=1,nrsl
        call clear_rsp2(rsl(i))
      end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
!
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
      do j=sta(1),sto(1),incr !1,nseq
        do k=j,nseq
          call Ven_rsp_long(evec_thr(:,tpi),j,k,cut)
        end do
        call en_freesolv(evec_thr(:,tpi),j)
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
!
  else if (use_IMPSOLV.EQV..true.) then
#ifdef ENABLE_THREADS
    call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,incr)
    tpn = omp_get_num_threads()
    incr = tpn
    tpi = omp_get_thread_num() + 1
    sto(1) = nseq
    sta(1) = tpi
#else
    tpi = 1
    sta(1) = 1
    sto(1) = nseq
    incr = 1
#endif
    do j=sta(1),sto(1),incr !1,nseq
      call en_freesolv(evec_thr(:,tpi),j)
    end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
    if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
!
  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.)) then
    if (cut.EQV..true.) then
      if (use_mcgrid.EQV..true.) then
        do i=irs,frs
          call grd_respairs(i)
        end do
      else if (use_rescrit.EQV..true.) then
        do i=irs,frs
          call respairs(i)
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in pivot_energy_&
 &long(...).'
        call fexit()
      end if
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
      do j=irs,frs
        do k=sta(1),sto(1),incr !1,irs-1 and frs+1,nseq
          if ((k.ge.irs).AND.(k.le.frs)) cycle 
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
      end do
      sto(2) = frs
      sta(2) = irs + tpi - 1
      do j=sta(2),sto(2),incr
        if (rsp_mat(min(j,rs),max(j,rs)).ne.0) then
          call Ven_rsp_long(evec_thr(:,tpi),min(rs,j),max(rs,j),use_cutoffs)
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
      do i=irs,frs
        call clear_rsp2(i)
      end do
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
!      do j=sta(1),sto(1),incr !1,nseq
      do j=irs,frs
        do k=sta(1),sto(1),incr !1,irs-1 and frs+1,nseq
          if ((k.ge.irs).AND.(k.le.frs)) cycle
          call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
        end do
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
        call Ven_rsp_long(evec_thr(:,tpi),rs,j,cut)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
!
  else
!  do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     now update the atsav-array using svte (to be ready for posterior long-range calculation)
      call init_svte(5)
    else if (mode.eq.1) then
!     clear out the residue-svte-monitoring-vector
      do i=1,nseq
        rs_vec(i) = 0
        rs_vec2(i) = 0
      end do
    end if
  end if
!
end
!
!------------------------------------------------------------------------
!
! this function is almost exactly identical to pivot_energy_short except that
! it takes residue boudns as arguments (would be much less code to have a wrapper)

subroutine other_energy_short(rs,irs,frs,evec,cut,mode)
!
  use iounit
  use sequen
  use cutoffs
  use energies
  use molecule
  use fyoc
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer j,mode,i,k,irs,frs,rs,sta(2),sto(2),incr,tpi,imol
  RTYPE evec(MAXENERGYTERMS)
  logical cut
#ifdef ENABLE_THREADS
  integer tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
  maxtpn = 0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
  imol = molofrs(rs)
!
! the boundary terms depend only on coordinates -> can be taken care of immediately
#ifdef ENABLE_THREADS
  call omp_set_num_threads(min(frs-irs+1,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,incr)
  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  sto(1) = frs
  sta(1) = irs + tpi - 1
  sta(2) = max(rsmol(imol,1),irs-1) + tpi - 1
  sto(2) = min(rsmol(imol,2),frs+1)
  incr = tpn
#else
  sta(1) = irs
  sto(1) = frs
  sta(2) = max(rsmol(imol,1),irs-1)
  sto(2) = min(rsmol(imol,2),frs+1)
  incr = 1
  tpi = 1
#endif
  do j=sta(1),sto(1),incr
    call e_boundary_rs(j,evec_thr(:,tpi),mode)
  end do
  if (n_crosslinks.gt.0) then
    do j=sta(2),sto(2),incr
      if (disulf(j).gt.0) then
        if (abs(rs-j).le.1) then
          if (use_BOND(1).EQV..true.) call en_bonds(disulf(j),evec_thr(:,tpi))
          if (use_BOND(2).EQV..true.) call en_angles(disulf(j),evec_thr(:,tpi))
          if (use_BOND(3).EQV..true.) call en_impropers(disulf(j),evec_thr(:,tpi))
          if (use_BOND(4).EQV..true.) call en_torsions(disulf(j),evec_thr(:,tpi))
          if (use_BOND(5).EQV..true.) call en_cmap(disulf(j),evec_thr(:,tpi))
        else if ((disulf(j).gt.sto(2)).OR.(disulf(j).lt.sta(2))) then
          if (use_BOND(1).EQV..true.) call en_bonds(j,evec_thr(:,tpi))
          if (use_BOND(2).EQV..true.) call en_angles(j,evec_thr(:,tpi))
          if (use_BOND(3).EQV..true.) call en_impropers(j,evec_thr(:,tpi))
          if (use_BOND(4).EQV..true.) call en_torsions(j,evec_thr(:,tpi))
          if (use_BOND(5).EQV..true.) call en_cmap(j,evec_thr(:,tpi))
          if (use_BOND(1).EQV..true.) call en_bonds(disulf(j),evec_thr(:,tpi))
          if (use_BOND(2).EQV..true.) call en_angles(disulf(j),evec_thr(:,tpi))
          if (use_BOND(3).EQV..true.) call en_impropers(disulf(j),evec_thr(:,tpi))
          if (use_BOND(4).EQV..true.) call en_torsions(disulf(j),evec_thr(:,tpi))
          if (use_BOND(5).EQV..true.) call en_cmap(disulf(j),evec_thr(:,tpi))
        end if
      end if
    end do
  end if

#ifdef ENABLE_THREADS
!$OMP SINGLE
  if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP BARRIER
!$OMP END PARALLEL
#endif
!
! bonded terms may change but are restricted to a truly local effect (note -1/+1-padding, though)
#ifdef ENABLE_THREADS
  call omp_set_num_threads(min(min(rsmol(imol,2),rs+1)-max(rsmol(imol,1),rs-1)+1,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,incr) FIRSTPRIVATE(imol)
  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  sto(1) = min(rsmol(imol,2),rs+1)
  sta(1) = tpi - 1 + max(rsmol(imol,1),rs-1)
  incr = tpn
#else
  sta(1) = max(rsmol(imol,1),rs-1)
  sto(1) = min(rsmol(imol,2),rs+1)
  tpi = 1
  incr = 1
#endif
  do j=sta(1),sto(1),incr
    if (use_BOND(1).EQV..true.) call en_bonds(j,evec_thr(:,tpi))
    if (use_BOND(2).EQV..true.) call en_angles(j,evec_thr(:,tpi))
    if (use_BOND(3).EQV..true.) call en_impropers(j,evec_thr(:,tpi))
    if (use_BOND(4).EQV..true.) call en_torsions(j,evec_thr(:,tpi))
    if (use_BOND(5).EQV..true.) call en_cmap(j,evec_thr(:,tpi))
  end do
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
  if (tpn.gt.maxtpn) maxtpn = tpn
#endif
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     zero out svte and svbu (temporary arrays for changes in SAV)
      call init_svte(0)
    else if (mode.eq.1) then
!     store svte (the changes induced by the sampling) in svbu, zero out svte again
      call init_svte(1)
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP BARRIER
!$OMP END PARALLEL
#endif
!
! let's only do all this setup work if necessary
  if ((use_IMPSOLV.EQV..true.).OR.(use_WCA.EQV..true.).OR.(use_IPP.EQV..true.).OR.&
 &    (use_attLJ.EQV..true.).OR.(use_CORR.EQV..true.).OR.(use_FEGS(1).EQV..true.).OR.&
 &    (use_FEGS(3).EQV..true.)) then
!
    if (cut.EQV..true.) then
      if (use_mcgrid.EQV..true.) then
        do i=irs,frs
          call grd_respairs(i)
        end do
      else if (use_rescrit.EQV..true.) then
        do i=irs,frs
          call respairs(i)
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in pivot_energy_&
 &long(...).'
        call fexit()
      end if
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
      do j=irs,frs
        do k=sta(1),sto(1),incr !1,irs-1 and frs+1,nseq
          if ((k.ge.irs).AND.(k.le.frs)) cycle 
          if (rsp_mat(max(j,k),min(j,k)).ne.0) then
            call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
      end do
      sto(2) = frs
      sta(2) = irs + tpi - 1
      do j=sta(2),sto(2),incr
        if (rsp_mat(max(j,rs),min(j,rs)).ne.0) then
          call Ven_rsp(evec_thr(:,tpi),min(rs,j),max(rs,j),use_cutoffs)
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
      do i=irs,frs
        call clear_rsp2(i)
      end do
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
!      do j=sta(1),sto(1),incr !1,nseq
      do j=irs,frs
        do k=sta(1),sto(1),incr !1,irs-1 and frs+1,nseq
          if ((k.ge.irs).AND.(k.le.frs)) cycle
          call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),cut)
        end do
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
        call Ven_rsp(evec_thr(:,tpi),rs,j,cut)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
!
  else
!   do nothing
  end if
!
!
  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.1) then
      call init_svte(2)
    end if
  end if
!
  end
!
!--------------------------------------------------------------------------
!
! exactly identical to pivot_energy_long except for arguments
!
subroutine other_energy_long(rs,irs,frs,evec,cut,mode)
!
  use iounit
  use sequen
  use cutoffs
  use energies
  use molecule
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rs,j,mode,i,k,l,rsl(nseq),nrsl,irs,frs,sta(2),sto(2),incr,tpi
  RTYPE evec(MAXENERGYTERMS)
  logical cut,inrsl(nseq)
#ifdef ENABLE_THREADS
  integer tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
  maxtpn = 0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    call setup_scrqs2()
  end if
!
  nrsl = 0
!     
  if (use_TOR.EQV..true.) then
    if (par_TOR2(rs).gt.0) then
      call en_torrs(evec,rs)
    end if
  end if
!
  if ((use_IMPSOLV.EQV..true.).AND.((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.))) then
!
!   the problematic part is that the interaction of two residues can
!   indeed be changed by a third residue, which is changing conformation.
!   it is henceforth not safe to compute just the terms that change
!   by virtue of the conformational change.
!   to detect the implicitly affected terms we're using the temporary SAV
!   array (i.e., any atom whose SAV was affected by the conformational
!   change potentially interacts differently with all other atoms now)
!   through the rs_vec-array which is created during the last init_svte-call
!   in the short-range energy computation
!   note, however, that just like for every other interaction, the residue
!   based cutoff criterion still has to be fulfilled
!
    if (cut.EQV..true.) then
      inrsl(:) = .false.
      do j=1,irs-1
       if (rs_vec(j).eq.1) then
          nrsl = nrsl + 1
          rsl(nrsl) = j
          inrsl(j) = .true.
        end if
      end do
      do j=irs,frs
        nrsl = nrsl + 1
        rsl(nrsl) = j
        inrsl(j) = .true.
      end do
      do j=frs+1,nseq
        if (rs_vec(j).eq.1) then
          nrsl = nrsl + 1
          rsl(nrsl) = j
          inrsl(j) = .true.
        end if
      end do
      if (use_mcgrid.EQV..true.) then
        do i=1,nrsl
          call grd_respairs(rsl(i))
        end do
      else if (use_rescrit.EQV..true.) then
        do i=1,nrsl
          call respairs(rsl(i))
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in pivot_energy_&
 &long(...).'
        call fexit()
      end if
!
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,i,j,k,l,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
      sto(2) = nrsl
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      sto(2) = nrsl
      incr = 1
#endif
      do i=1,nrsl
        j = rsl(i)
        do k=sta(1),sto(1),incr !1,nseq
!         to avoid double-counting, we'll utilize the inverse map
          if (inrsl(k).EQV..true.) cycle
!
          if (((j.ge.k).AND.(rsp_mat(k,j).ne.0)).OR.&
 &            ((k.gt.j).AND.(rsp_mat(j,k).ne.0))) then
            call Ven_rsp_long(evec_thr(:,tpi),j,k,use_cutoffs)
          end if
        end do
        sta(2) = i + tpi - 1
!       and take care of the rest here
        do l=sta(2),sto(2),incr !do l=i,nrsl
          k = rsl(l)
          if (((j.ge.k).AND.(rsp_mat(k,j).ne.0)).OR.&
 &            ((k.gt.j).AND.(rsp_mat(j,k).ne.0))) then
            call Ven_rsp_long(evec_thr(:,tpi),j,k,use_cutoffs)
          end if
        end do
      end do
!
      do i=sta(1),sto(1),incr ! 1,nseq
        if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
          call en_freesolv(evec_thr(:,tpi),i)
        end if
      end do
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
#endif
      do i=1,nrsl
        call clear_rsp2(rsl(i))
      end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
!
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
      do j=sta(1),sto(1),incr !1,nseq
        do k=j,nseq
          call Ven_rsp_long(evec_thr(:,tpi),j,k,cut)
        end do
        call en_freesolv(evec_thr(:,tpi),j)
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
!
  else if (use_IMPSOLV.EQV..true.) then
#ifdef ENABLE_THREADS
    call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,incr)
    tpn = omp_get_num_threads()
    incr = tpn
    tpi = omp_get_thread_num() + 1
    sto(1) = nseq
    sta(1) = tpi
#else
    tpi = 1
    sta(1) = 1
    sto(1) = nseq
    incr = 1
#endif
    do j=sta(1),sto(1),incr !1,nseq
      call en_freesolv(evec_thr(:,tpi),j)
    end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
    if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
!
  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.)) then
    if (cut.EQV..true.) then
      if (use_mcgrid.EQV..true.) then
        do i=irs,frs
          call grd_respairs(i)
        end do
      else if (use_rescrit.EQV..true.) then
        do i=irs,frs
          call respairs(i)
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in pivot_energy_&
 &long(...).'
        call fexit()
      end if
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
      do j=irs,frs
        do k=sta(1),sto(1),incr !1,irs-1 and frs+1,nseq
          if ((k.ge.irs).AND.(k.le.frs)) cycle 
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
      end do
      sto(2) = frs
      sta(2) = irs + tpi - 1
      do j=sta(2),sto(2),incr
        if (rsp_mat(min(j,rs),max(j,rs)).ne.0) then
          call Ven_rsp_long(evec_thr(:,tpi),min(rs,j),max(rs,j),use_cutoffs)
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
      do i=irs,frs
        call clear_rsp2(i)
      end do
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
!      do j=sta(1),sto(1),incr !1,nseq
      do j=irs,frs
        do k=sta(1),sto(1),incr !1,irs-1 and frs+1,nseq
          if ((k.ge.irs).AND.(k.le.frs)) cycle
          call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
        end do
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
        call Ven_rsp_long(evec_thr(:,tpi),rs,j,cut)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
!
  else
!  do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     now update the atsav-array using svte (to be ready for posterior long-range calculation)
      call init_svte(5)
    else if (mode.eq.1) then
!     clear out the residue-svte-monitoring-vector
      do i=1,nseq
        rs_vec(i) = 0
        rs_vec2(i) = 0
      end do
    end if
  end if
!
end
!
!------------------------------------------------------------------------
!
! a simple wrapper routine that handles short-range energy calculations for pivot
! moves (only relevant terms are computed ...)
!
subroutine cr_energy_short(rsi,rsf,evec,cut,mode,ct)
!
  use sequen
  use iounit
  use energies
  use cutoffs
  use molecule
  use fyoc
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rsi,rsf,j,k,mode,irs,frs,imol,sta,sto,sta2,sto2,incr,tpi
  RTYPE evec(MAXENERGYTERMS)
  logical cut,ct
#ifdef ENABLE_THREADS
  integer tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,pvl
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
  maxtpn = 0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
!
! set residue range (lever arm) based on whether N- or C-aligned (the moving parts change!)
  if (ct.EQV..true.) then
    irs = rsmol(molofrs(rsi),1)
    frs = rsf
  else
    irs = rsi
    frs = rsmol(molofrs(rsi),2)
  end if
!
! the boundary terms depend only on coordinates -> can be taken care of immediately
#ifdef ENABLE_THREADS
  pvl = frs - irs + 1
  call omp_set_num_threads(min(pvl,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,incr)
  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  sto = frs
  sta = tpi - 1 + irs
  sta2 = max(rsmol(molofrs(rsi),1),irs-1) + tpi - 1
  sto2 = min(rsmol(molofrs(rsi),2),frs+1)
  incr = tpn
#else
  sta = irs
  sto = frs
  sta2 = max(rsmol(molofrs(rsi),1),irs-1)
  sto2 = min(rsmol(molofrs(rsi),2),frs+1)
  incr = 1
  tpi = 1
#endif
  do j=sta,sto,incr
    call e_boundary_rs(j,evec_thr(:,tpi),mode)
  end do
  if (n_crosslinks.gt.0) then
    do j=sta2,sto2,incr
      if (disulf(j).gt.0) then
        if ((j.le.(rsf+1)).AND.(j.ge.(rsi-1)).AND.&
 &          ((disulf(j).gt.(rsf+1)).OR.(disulf(j).lt.(rsi-1)))) then
          if (use_BOND(1).EQV..true.) call en_bonds(disulf(j),evec_thr(:,tpi))
          if (use_BOND(2).EQV..true.) call en_angles(disulf(j),evec_thr(:,tpi))
          if (use_BOND(3).EQV..true.) call en_impropers(disulf(j),evec_thr(:,tpi))
          if (use_BOND(4).EQV..true.) call en_torsions(disulf(j),evec_thr(:,tpi))
          if (use_BOND(5).EQV..true.) call en_cmap(disulf(j),evec_thr(:,tpi))
        else if ((disulf(j).gt.sto2).OR.(disulf(j).lt.sta2)) then
          if (use_BOND(1).EQV..true.) call en_bonds(j,evec_thr(:,tpi))
          if (use_BOND(2).EQV..true.) call en_angles(j,evec_thr(:,tpi))
          if (use_BOND(3).EQV..true.) call en_impropers(j,evec_thr(:,tpi))
          if (use_BOND(4).EQV..true.) call en_torsions(j,evec_thr(:,tpi))
          if (use_BOND(5).EQV..true.) call en_cmap(j,evec_thr(:,tpi))
          if (use_BOND(1).EQV..true.) call en_bonds(disulf(j),evec_thr(:,tpi))
          if (use_BOND(2).EQV..true.) call en_angles(disulf(j),evec_thr(:,tpi))
          if (use_BOND(3).EQV..true.) call en_impropers(disulf(j),evec_thr(:,tpi))
          if (use_BOND(4).EQV..true.) call en_torsions(disulf(j),evec_thr(:,tpi))
          if (use_BOND(5).EQV..true.) call en_cmap(disulf(j),evec_thr(:,tpi))
        end if
      end if
    end do
  end if
#ifdef ENABLE_THREADS
!$OMP SINGLE
  if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP BARRIER
!$OMP END PARALLEL
#endif
!
! bonded terms may change but are restricted to a truly local effect (note -1/+1-padding, though)
  imol = molofrs(rsi)
#ifdef ENABLE_THREADS
  pvl = min(rsmol(imol,2),rsf+1) - max(rsmol(imol,1),rsi-1) + 1
  call omp_set_num_threads(min(pvl,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,incr,tpn,tpi,j)
  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  sto = min(rsmol(imol,2),rsf+1)
  sta = max(rsmol(imol,1),rsi-1) + tpi - 1
  incr = tpn
#else
  sta = max(rsmol(imol,1),rsi-1)
  sto = min(rsmol(imol,2),rsf+1)
  tpi = 1
  incr = 1
#endif
  do j=sta,sto,incr
    if (use_BOND(1).EQV..true.) call en_bonds(j,evec_thr(:,tpi))
    if (use_BOND(2).EQV..true.) call en_angles(j,evec_thr(:,tpi))
    if (use_BOND(3).EQV..true.) call en_impropers(j,evec_thr(:,tpi))
    if (use_BOND(4).EQV..true.) call en_torsions(j,evec_thr(:,tpi))
    if (use_BOND(5).EQV..true.) call en_cmap(j,evec_thr(:,tpi))
  end do
#ifdef ENABLE_THREADS
!$OMP SINGLE
  if (tpn.gt.maxtpn) maxtpn = tpn
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     zero out svte and svbu (temporary arrays for changes in SAV)
      call init_svte(0)
    else if (mode.eq.1) then
!     store svte (the changes induced by the sampling) in svbu, zero out svte again
      call init_svte(1)
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP BARRIER
!$OMP END PARALLEL
#endif
!
! let's only do all this setup work if necessary
  if ((use_IMPSOLV.EQV..true.).OR.(use_WCA.EQV..true.).OR.(use_IPP.EQV..true.).OR.&
 &    (use_attLJ.EQV..true.).OR.(use_CORR.EQV..true.).OR.(use_FEGS(1).EQV..true.).OR.&
 &    (use_FEGS(3).EQV..true.)) then
!
    if (cut.EQV..true.) then
      if (use_mcgrid.EQV..true.) then
        do j=irs,frs
          call grd_respairs(j)
        end do
      else if (use_rescrit.EQV..true.) then
        do j=irs,frs
          call respairs(j)
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in cr_energy_short&
 &(...).'
        call fexit()
      end if
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto = nseq
      sta = tpi
#else
      tpi = 1
      sta = 1
      sto = nseq
      incr = 1
#endif
      do k=sta,sto,incr
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            if (rsp_mat(max(j,k),min(j,k)).ne.0) then
              call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),cut)
            end if
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            if (rsp_mat(max(j,k),min(j,k)).ne.0) then
              call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),cut)
            end if
          end do
        else 
!         CR-stretch to rest of lever arm
          do j=rsi,rsf
            if (rsp_mat(max(j,k),min(j,k)).ne.0) then
              call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),cut)
            end if
          end do
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
      do j=irs,frs
        call clear_rsp2(j)
      end do
!
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto = nseq
      sta = tpi
#else
      tpi = 1
      sta = 1
      sto = nseq
      incr = 1
#endif
      do k=sta,sto,incr
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),cut)
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),cut)
          end do
        else 
!         CR-stretch to rest of lever arm
          do j=rsi,rsf
            call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),cut)
          end do
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
!
  else
!   do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.1) then
!     get the difference in SAV
      call init_svte(2)
    end if
  end if
!
end
!
!--------------------------------------------------------------------------
!
! a simple wrapper routine that handles long-range energy calculations for pivot
! moves (only relevant terms are computed ...)
!
subroutine cr_energy_long(rsi,rsf,evec,cut,mode,ct)
!
  use iounit
  use sequen
  use cutoffs
  use energies
  use molecule
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rsi,rsf,j,mode,i,k,l,rs,nrsl,rsl(nseq),frs,irs,sta(2),sto(2),incr,tpi
  RTYPE evec(MAXENERGYTERMS)
  logical cut,ct,inrsl(nseq)
#ifdef ENABLE_THREADS
  integer tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,pvl
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
  maxtpn = 0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
  if (ct.EQV..true.) then
    irs = rsmol(molofrs(rsi),1)
    frs = rsf
  else
    irs = rsi
    frs = rsmol(molofrs(rsi),2)
  end if
!
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    call setup_scrqs2()
  end if
!
  nrsl = 0
!
  if (use_TOR.EQV..true.) then
#ifdef ENABLE_THREADS
    pvl = rsf - rsi + 1
    call omp_set_num_threads(min(pvl,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,i,incr)
    tpn = omp_get_num_threads()
    tpi = omp_get_thread_num() + 1
    sto(1) = tpi*pvl/tpn + rsi - 1
    sta(1) = (tpi-1)*pvl/tpn + rsi
    incr = tpn
#else
    sta(1) = rsi
    sto(1) = rsf
    incr = 1
    tpi = 1
#endif
    do i=sta(1),sto(1),incr !rsi,rsf
      if (par_TOR2(i).gt.0) then
        call en_torrs(evec_thr(:,tpi),i)
      end if
    end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
    if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
  end if
!
  if ((use_IMPSOLV.EQV..true.).AND.((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.))) then
!
!   the problematic part is that the interaction of two residues can
!   indeed be changed by a third residue, which is changing conformation.
!   it is henceforth not safe to compute just the terms that change
!   by virtue of the conformational change.
!   to detect the implicitly affected terms we're using the temporary SAV
!   array (i.e., any atom whose SAV was affected by the conformational
!   change potentially interacts differently with all other atoms now)
!   through the rs_vec-array which is created during the last init_svte-call
!   in the short-range energy computation
!   note, however, that just like for every other interaction, the residue
!   based cutoff criterion still has to be fulfilled
!
    if (cut.EQV..true.) then
      inrsl(:) = .false.
      do rs=1,irs-1
       if (rs_vec(rs).eq.1) then
          nrsl = nrsl + 1
          rsl(nrsl) = rs
          inrsl(rs) = .true.
        end if
      end do
      do rs=irs,frs
        nrsl = nrsl + 1
        rsl(nrsl) = rs
        inrsl(rs) = .true.
      end do
      do rs=frs+1,nseq
        if (rs_vec(rs).eq.1) then
          nrsl = nrsl + 1
          rsl(nrsl) = rs
          inrsl(rs) = .true.
        end if
      end do
      if (use_mcgrid.EQV..true.) then
        do i=1,nrsl
          call grd_respairs(rsl(i))
        end do
      else if (use_rescrit.EQV..true.) then
        do i=1,nrsl
          call respairs(rsl(i))
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in cr_energy_lon&
 &g(...).'
        call fexit()
      end if
!
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,i,j,k,l,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
      sto(2) = nrsl
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      sto(2) = nrsl
      incr = 1
#endif
      do i=1,nrsl
        j = rsl(i)
        do k=sta(1),sto(1),incr ! 1,nseq
!         we use the inverse map to circumvent problems with double counting
          if (inrsl(k).EQV..true.) cycle
!
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
!       and do the rest here
        sta(2) = i + tpi - 1
        do l=sta(2),sto(2),incr !i,nrsl
          k = rsl(l)
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
      end do
!
      do i=sta(1),sto(1),incr !1,nseq
        if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
          call en_freesolv(evec_thr(:,tpi),i)
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
#endif
      do i=1,nrsl
        call clear_rsp2(rsl(i))
      end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
!
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
      do j=sta(1),sto(1),incr !1,nseq
        do k=j,nseq
          call Ven_rsp_long(evec_thr(:,tpi),j,k,cut)
        end do
        call en_freesolv(evec_thr(:,tpi),j)
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
!
  else if (use_IMPSOLV.EQV..true.) then
#ifdef ENABLE_THREADS
    call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,incr)
    tpn = omp_get_num_threads()
    incr = tpn
    tpi = omp_get_thread_num() + 1
    sto(1) = nseq
    sta(1) = tpi
#else
    tpi = 1
    sta(1) = 1
    sto(1) = nseq
    incr = 1
#endif
    do j=sta(1),sto(1),incr !1,nseq
      call en_freesolv(evec_thr(:,tpi),j)
    end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
    if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
!
  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.)) then
    if (cut.EQV..true.) then
      if (use_mcgrid.EQV..true.) then
        do i=irs,frs
          call grd_respairs(i)
        end do
      else if (use_rescrit.EQV..true.) then
        do i=irs,frs
          call respairs(i)
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in pivot_energy_&
 &long(...).'
        call fexit()
      end if
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
      do k=sta(1),sto(1),incr ! 1nseq
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            if (rsp_mat(min(j,k),max(j,k)).ne.0) then
              call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
            end if
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            if (rsp_mat(min(j,k),max(j,k)).ne.0) then
              call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
            end if
          end do
        else 
!         CR-stretch to rest of lever arm
          do j=rsi,rsf
            if (rsp_mat(min(j,k),max(j,k)).ne.0) then
              call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
            end if
          end do
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
      do j=irs,frs
        call clear_rsp2(j)
      end do
!
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
      do k=sta(1),sto(1),incr ! 1nseq
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
          end do
        else 
!         CR-stretch to rest of lever arm
          do j=rsi,rsf
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
          end do
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
!
    end if
!
  else
!   do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     now update the atsav-array using svte (to be ready for posterior long-range calculation)
      call init_svte(5)
    else if (mode.eq.1) then
!     clear out the residue-svte-monitoring-vector
      do i=1,nseq
        rs_vec(i) = 0
        rs_vec2(i) = 0
      end do
    end if
  end if
!
end
!
!
!------------------------------------------------------------------------
!
! a simple wrapper routine that handles short-range energy calculations for truly
! local CR moves (only relevant terms are computed ...)
! for long peptides, the complexity is significantly reduced when compared to pivot
! moves
!
subroutine ujcr_energy_short(rsi,rsf,evec,cut,mode)
!
  use sequen
  use iounit
  use energies
  use cutoffs
  use molecule
  use fyoc
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rsi,rsf,j,k,mode,irs,frs,imol,sta,sto,incr,tpi
  RTYPE evec(MAXENERGYTERMS)
  logical cut
#ifdef ENABLE_THREADS
  integer tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,pvl
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
  maxtpn = 0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
  irs = rsi
  frs = rsf
!
! the boundary terms depend only on coordinates -> can be taken care of immediately
#ifdef ENABLE_THREADS
  pvl = frs - irs + 1
  call omp_set_num_threads(min(pvl,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,incr)
  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  sto = rsf
  sta = tpi + irs - 1
  incr = tpn
#else
  sta = irs
  sto = frs
  tpi = 1
  incr = 1
#endif
  do j=sta,sto,incr
    call e_boundary_rs(j,evec_thr(:,tpi),mode)
  end do
#ifdef ENABLE_THREADS
!$OMP SINGLE
  if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP BARRIER
!$OMP END PARALLEL
#endif
!
! bonded terms may change but are restricted to a truly local effect (note -1/+1-padding, though)
  imol = molofrs(rsi)
#ifdef ENABLE_THREADS
  pvl = min(rsmol(imol,2),rsf+1) - max(rsmol(imol,1),rsi-1) + 1
  call omp_set_num_threads(min(pvl,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,incr)
  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  sto = min(rsmol(imol,2),rsf+1)
  sta = tpi - 1 + max(rsmol(imol,1),rsi-1)
  incr = 1
#else
  sta = max(rsmol(imol,1),rsi-1)
  sto = min(rsmol(imol,2),rsf+1)
  tpi = 1
  incr = 1
#endif
  do j=sta,sto,incr
    if (use_BOND(1).EQV..true.) call en_bonds(j,evec_thr(:,tpi))
    if (use_BOND(2).EQV..true.) call en_angles(j,evec_thr(:,tpi))
    if (use_BOND(3).EQV..true.) call en_impropers(j,evec_thr(:,tpi))
    if (use_BOND(4).EQV..true.) call en_torsions(j,evec_thr(:,tpi))
    if (use_BOND(5).EQV..true.) call en_cmap(j,evec_thr(:,tpi))
    if (disulf(j).gt.0) then
      if ((disulf(j).lt.sta).OR.(disulf(j).gt.sto)) then
        if (use_BOND(1).EQV..true.) call en_bonds(disulf(j),evec_thr(:,tpi))
        if (use_BOND(2).EQV..true.) call en_angles(disulf(j),evec_thr(:,tpi))
        if (use_BOND(3).EQV..true.) call en_impropers(disulf(j),evec_thr(:,tpi))
        if (use_BOND(4).EQV..true.) call en_torsions(disulf(j),evec_thr(:,tpi))
        if (use_BOND(5).EQV..true.) call en_cmap(disulf(j),evec_thr(:,tpi))
      end if
    end if
  end do
#ifdef ENABLE_THREADS
!$OMP SINGLE
  if (tpn.gt.maxtpn) maxtpn = tpn
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     zero out svte and svbu (temporary arrays for changes in SAV)
      call init_svte(0)
    else if (mode.eq.1) then
!     store svte (the changes induced by the sampling) in svbu, zero out svte again
      call init_svte(1)
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP BARRIER
!$OMP END PARALLEL
#endif
!
! let's only do all this setup work if necessary
  if ((use_IMPSOLV.EQV..true.).OR.(use_WCA.EQV..true.).OR.(use_IPP.EQV..true.).OR.&
 &    (use_attLJ.EQV..true.).OR.(use_CORR.EQV..true.).OR.(use_FEGS(1).EQV..true.).OR.&
 &    (use_FEGS(3).EQV..true.)) then
!
    if (cut.EQV..true.) then
      if (use_mcgrid.EQV..true.) then
        do j=irs,frs
          call grd_respairs(j)
        end do
      else if (use_rescrit.EQV..true.) then
        do j=irs,frs
          call respairs(j)
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in ujcr_energy_sho&
 &rt(...).'
        call fexit()
      end if
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto = nseq
      sta = tpi
#else
      tpi = 1
      sta = 1
      sto = nseq
      incr = 1
#endif
      do k=sta,sto,incr ! 1nseq
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            if (rsp_mat(max(j,k),min(j,k)).ne.0) then
              call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),cut)
            end if
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            if (rsp_mat(max(j,k),min(j,k)).ne.0) then
              call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),cut)
            end if
          end do
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
      do j=irs,frs
        call clear_rsp2(j)
      end do
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto = nseq
      sta = tpi
#else
      tpi = 1
      sta = 1
      sto = nseq
      incr = 1
#endif
      do k=sta,sto,incr ! 1nseq
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),cut)
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),cut)
          end do
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
!
  else
!   do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.1) then
!     get the difference in SAV
      call init_svte(2)
    end if
  end if
!
end
!
!--------------------------------------------------------------------------
!
! a simple wrapper routine that handles long-range energy calculations for pivot
! moves (only relevant terms are computed ...)
!
subroutine ujcr_energy_long(rsi,rsf,evec,cut,mode)
!
  use iounit
  use sequen
  use cutoffs
  use energies
  use molecule
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rsi,rsf,j,mode,i,k,l,rs,nrsl,rsl(nseq),frs,irs,sta(2),sto(2),incr,tpi
  RTYPE evec(MAXENERGYTERMS)
  logical cut,inrsl(nseq)
#ifdef ENABLE_THREADS
  integer tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,pvl
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
  maxtpn = 0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
  frs = rsf
  irs = rsi
!
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    call setup_scrqs2()
  end if
!
  nrsl = 0
!
  if (use_TOR.EQV..true.) then
#ifdef ENABLE_THREADS
    pvl = rsf - rsi + 1
    call omp_set_num_threads(min(pvl,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,i,incr)
    tpn = omp_get_num_threads()
    tpi = omp_get_thread_num() + 1
    sto(1) = tpi*pvl/tpn + rsi - 1
    sta(1) = (tpi-1)*pvl/tpn + rsi
    incr = tpn
#else
    sta(1) = rsi
    sto(1) = rsf
    incr = 1
    tpi = 1
#endif
    do i=sta(1),sto(1),incr !rsi,rsf
      if (par_TOR2(i).gt.0) then
        call en_torrs(evec_thr(:,tpi),i)
      end if
    end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
    if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
  end if
!
  if (use_IMPSOLV.EQV..true.) then
!
!   the problematic part is that the interaction of two residues can
!   indeed be changed by a third residue, which is changing conformation.
!   it is henceforth not safe to compute just the terms that change
!   by virtue of the conformational change.
!   to detect the implicitly affected terms we're using the temporary SAV
!   array (i.e., any atom whose SAV was affected by the conformational
!   change potentially interacts differently with all other atoms now)
!   through the rs_vec-array which is created during the last init_svte-call
!   in the short-range energy computation
!   note, however, that just like for every other interaction, the residue
!   based cutoff criterion still has to be fulfilled
!
    if (cut.EQV..true.) then
      inrsl(:) = .false.
      do rs=1,irs-1
       if (rs_vec(rs).eq.1) then
          nrsl = nrsl + 1
          rsl(nrsl) = rs
          inrsl(rs) = .true.
        end if
      end do
      do rs=irs,frs
        nrsl = nrsl + 1
        rsl(nrsl) = rs
        inrsl(rs) = .true.
      end do
      do rs=frs+1,nseq
        if (rs_vec(rs).eq.1) then
          nrsl = nrsl + 1
          rsl(nrsl) = rs
          inrsl(rs) = .true.
        end if
      end do
      if (use_mcgrid.EQV..true.) then
        do i=1,nrsl
          call grd_respairs(rsl(i))
        end do
      else if (use_rescrit.EQV..true.) then
        do i=1,nrsl
          call respairs(rsl(i))
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in ujcr_energy_l&
 &ong(...).'
        call fexit()
      end if
!
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,i,j,k,l,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
      sto(2) = nrsl
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      sto(2) = nrsl
      incr = 1
#endif
      do i=1,nrsl
        j = rsl(i)
        do k=sta(1),sto(1),incr ! 1,nseq
!         we solve the issue of double counting by using the inverse map
          if (inrsl(k).EQV..true.) cycle
!
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
!       and deal with the rest here
        sta(2) = i + tpi - 1
        do l=sta(2),sto(2),incr ! i,nrsl
          k = rsl(l)
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
      end do
      do i=sta(1),sto(1),incr ! 1,nseq
        if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
          call en_freesolv(evec_thr(:,tpi),i)
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
#endif
      do i=1,nrsl
        call clear_rsp2(rsl(i))
      end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
!
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
      do j=sta(1),sto(1),incr ! 1,nseq
        do k=j,nseq
          call Ven_rsp_long(evec_thr(:,tpi),j,k,cut)
        end do
        call en_freesolv(evec_thr(:,tpi),j)
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
!
  else if (use_IMPSOLV.EQV..true.) then

  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.)) then
    if (cut.EQV..true.) then
      if (use_mcgrid.EQV..true.) then
        do i=irs,frs
          call grd_respairs(i)
        end do
      else if (use_rescrit.EQV..true.) then
        do i=irs,frs
          call respairs(i)
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in ujcr_energy_l&
 &ong(...).'
        call fexit()
      end if
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
      do k=sta(1),sto(1),incr
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            if (rsp_mat(min(j,k),max(j,k)).ne.0) then
              call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
            end if
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            if (rsp_mat(min(j,k),max(j,k)).ne.0) then
              call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
            end if
          end do
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
      do j=irs,frs
        call clear_rsp2(j)
      end do
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sta,sto,tpn,tpi,j,k,incr)
      tpn = omp_get_num_threads()
      incr = tpn
      tpi = omp_get_thread_num() + 1
      sto(1) = nseq
      sta(1) = tpi
#else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
#endif
      do k=sta(1),sto(1),incr
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
          end do
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
!
  else
!   do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     now update the atsav-array using svte (to be ready for posterior long-range calculation)
      call init_svte(5)
    else if (mode.eq.1) then
!     clear out the residue-svte-monitoring-vector
      do i=1,nseq
        rs_vec(i) = 0
        rs_vec2(i) = 0
      end do
    end if
  end if
!
end
!
!--------------------------------------------------------------------------
!
! a simple wrapper routine that handles short-range energy calculations for rigid body
! moves (only relevant terms are computed ...)
!
subroutine rigid_energy_short(imol,evec,cut,mode)
!
  use sequen
  use iounit
  use energies
  use molecule
  use cutoffs
  use fyoc
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rs,imol,j,mode,tpi,sta,sto,incr
  RTYPE evec(MAXENERGYTERMS)
  logical cut
#ifdef ENABLE_THREADS
  integer tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,rvl
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
  maxtpn = 0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
! the boundary terms depend only on coordinates -> can be taken care of immediately
#ifdef ENABLE_THREADS
  rvl = rsmol(imol,2) - rsmol(imol,1) + 1
  call omp_set_num_threads(min(rvl,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rs,j,tpi,tpn,sto,sta,incr)
  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  sta = tpi + rsmol(imol,1) - 1
  sto = rsmol(imol,2)
  incr = tpn
#else
  sta = rsmol(imol,1)
  sto = rsmol(imol,2)
  incr = 1
  tpi = 1
#endif

  do j=sta,sto,incr ! rsmol(imol,1),rsmol(imol,2)
    call e_boundary_rs(j,evec_thr(:,tpi),mode)
    if (disulf(j).gt.0) then
      if ((disulf(j).lt.sta).OR.(disulf(j).gt.sto)) then
        if (use_BOND(1).EQV..true.) call en_bonds(j,evec_thr(:,tpi))
        if (use_BOND(2).EQV..true.) call en_angles(j,evec_thr(:,tpi))
        if (use_BOND(3).EQV..true.) call en_impropers(j,evec_thr(:,tpi))
        if (use_BOND(4).EQV..true.) call en_torsions(j,evec_thr(:,tpi))
        if (use_BOND(5).EQV..true.) call en_cmap(j,evec_thr(:,tpi))
        !
        if (use_BOND(1).EQV..true.) call en_bonds(disulf(j),evec_thr(:,tpi))
        if (use_BOND(2).EQV..true.) call en_angles(disulf(j),evec_thr(:,tpi))
        if (use_BOND(3).EQV..true.) call en_impropers(disulf(j),evec_thr(:,tpi))
        if (use_BOND(4).EQV..true.) call en_torsions(disulf(j),evec_thr(:,tpi))
        if (use_BOND(5).EQV..true.) call en_cmap(disulf(j),evec_thr(:,tpi))
      end if
    end if
  end do
#ifdef ENABLE_THREADS
!$OMP SINGLE
  if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP BARRIER
!$OMP END PARALLEL
#endif
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     zero out svte and svbu (temporary arrays for changes in SAV)
      call init_svte(0)
    else if (mode.eq.1) then
!     store svte (the changes induced by the sampling) in svbu, zero out svte again
      call init_svte(1)
    end if
  end if
! let's only do all this setup work if necessary
  if ((use_IMPSOLV.EQV..true.).OR.(use_WCA.EQV..true.).OR.(use_IPP.EQV..true.).OR.&
 &    (use_attLJ.EQV..true.).OR.(use_CORR.EQV..true.).OR.(use_FEGS(1).EQV..true.).OR.&
 &    (use_FEGS(3).EQV..true.)) then
!
    if ((cut.EQV..true.).AND.((use_rescrit.EQV..true.).OR.&
 &         (use_mcgrid.EQV..true.))) then
      if (use_mcgrid.EQV..true.) then
        do rs=rsmol(imol,1),rsmol(imol,2)
          call grd_respairs(rs)
        end do
      else if (use_rescrit.EQV..true.) then
        do rs=rsmol(imol,1),rsmol(imol,2)
          call respairs(rs)
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in rigid_energy_sh&
 &ort(...).'
        call fexit()
      end if
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rs,j,tpi,tpn,sto,sta,incr)
      tpn = omp_get_num_threads()
      tpi = omp_get_thread_num() + 1
      sta = tpi
      sto = nseq
      incr = tpn
#else
      sta = 1
      incr = 1
      sto = nseq
      tpi = 1
#endif

      do rs=rsmol(imol,1),rsmol(imol,2)
        do j=sta,sto,incr ! 1,rsmol(imol,1)-1,rsmol(imol,2)+1,nseq
          if ((j.ge.rsmol(imol,1)).AND.(j.le.rsmol(imol,2))) cycle
          if (rsp_mat(max(j,rs),min(j,rs)).ne.0) then
            call Ven_rsp(evec_thr(:,tpi),min(j,rs),max(j,rs),use_cutoffs)
          end if
        end do
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
      do rs=rsmol(imol,1),rsmol(imol,2)
        call clear_rsp2(rs)
      end do
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rs,j,tpi,tpn,sto,sta,incr)
      tpn = omp_get_num_threads()
      tpi = omp_get_thread_num() + 1
      sta = tpi
      sto = nseq
      incr = tpn
#else
      sta = 1
      incr = 1
      sto = nseq
      tpi = 1
#endif
      do rs=rsmol(imol,1),rsmol(imol,2)
        do j=sta,sto,incr ! 1,nseq
          if ((j.gt.rsmol(imol,2)).OR.(j.lt.rsmol(imol,1))) then 
            call Ven_rsp(evec_thr(:,tpi),rs,j,cut)
          end if
        end do
      end do

#ifdef ENABLE_THREADS
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP BARRIER
!$OMP END PARALLEL
#endif
    end if
!
  else
!   do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.1) then
!     get the difference in SAV
      call init_svte(2)
!     now use that difference to determine which dependent interactions we have to re-compute
    end if
  end if
!
end
!
!--------------------------------------------------------------------------
!
! the same for clusters of molecules
!
subroutine clurb_energy_short(imols,mmol,fmols,evec,cut,mode)
!
  use sequen
  use iounit
  use energies
  use molecule
  use cutoffs
  use system
  use fyoc
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rs,imol,j,mode,i,mmol,imols(nmol),fmols(nmol),k,sta,sto,incr,tpi,rsl(nseq),nrsl
  RTYPE evec(MAXENERGYTERMS)
  logical cut
#ifdef ENABLE_THREADS
  integer tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,sta2,sto2
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
  maxtpn = 0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
! the boundary terms depend only on coordinates -> can be taken care of immediately
#ifdef ENABLE_THREADS
  nrsl = 0
  do i=1,mmol 
    imol = imols(i)
    do j=rsmol(imol,1),rsmol(imol,2)
      nrsl = nrsl + 1
      rsl(nrsl) = j
    end do
  end do
  call omp_set_num_threads(min(nrsl,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,rs,tpi,tpn,sto,sta,incr)
  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  sta = tpi
  sto = nrsl
  sta2 = tpi
  sto2 = n_crosslinks
  incr = tpn
  do j=sta,sto,incr
    rs = rsl(j)
    call e_boundary_rs(rs,evec_thr(:,tpi),mode)
  end do
  do j=sta2,sto2,incr
    if (use_BOND(1).EQV..true.) call en_bonds(crosslink(j)%rsnrs(1),evec_thr(:,tpi))
    if (use_BOND(2).EQV..true.) call en_angles(crosslink(j)%rsnrs(1),evec_thr(:,tpi))
    if (use_BOND(3).EQV..true.) call en_impropers(crosslink(j)%rsnrs(1),evec_thr(:,tpi))
    if (use_BOND(4).EQV..true.) call en_torsions(crosslink(j)%rsnrs(1),evec_thr(:,tpi))
    if (use_BOND(5).EQV..true.) call en_cmap(crosslink(j)%rsnrs(1),evec_thr(:,tpi))
    if (use_BOND(1).EQV..true.) call en_bonds(crosslink(j)%rsnrs(2),evec_thr(:,tpi))
    if (use_BOND(2).EQV..true.) call en_angles(crosslink(j)%rsnrs(2),evec_thr(:,tpi))
    if (use_BOND(3).EQV..true.) call en_impropers(crosslink(j)%rsnrs(2),evec_thr(:,tpi))
    if (use_BOND(4).EQV..true.) call en_torsions(crosslink(j)%rsnrs(2),evec_thr(:,tpi))
    if (use_BOND(5).EQV..true.) call en_cmap(crosslink(j)%rsnrs(2),evec_thr(:,tpi))
  end do
!$OMP BARRIER
!$OMP SINGLE
  if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#else
  do i=1,mmol 
    imol = imols(i)
    do j=rsmol(imol,1),rsmol(imol,2)
      call e_boundary_rs(j,evec,mode)
    end do
  end do
! be safe
  tpi = 1
  do j=1,n_crosslinks
    if (use_BOND(1).EQV..true.) call en_bonds(crosslink(j)%rsnrs(1),evec_thr(:,tpi))
    if (use_BOND(2).EQV..true.) call en_angles(crosslink(j)%rsnrs(1),evec_thr(:,tpi))
    if (use_BOND(3).EQV..true.) call en_impropers(crosslink(j)%rsnrs(1),evec_thr(:,tpi))
    if (use_BOND(4).EQV..true.) call en_torsions(crosslink(j)%rsnrs(1),evec_thr(:,tpi))
    if (use_BOND(5).EQV..true.) call en_cmap(crosslink(j)%rsnrs(1),evec_thr(:,tpi))
    if (use_BOND(1).EQV..true.) call en_bonds(crosslink(j)%rsnrs(2),evec_thr(:,tpi))
    if (use_BOND(2).EQV..true.) call en_angles(crosslink(j)%rsnrs(2),evec_thr(:,tpi))
    if (use_BOND(3).EQV..true.) call en_impropers(crosslink(j)%rsnrs(2),evec_thr(:,tpi))
    if (use_BOND(4).EQV..true.) call en_torsions(crosslink(j)%rsnrs(2),evec_thr(:,tpi))
    if (use_BOND(5).EQV..true.) call en_cmap(crosslink(j)%rsnrs(2),evec_thr(:,tpi))
  end do
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     zero out svte and svbu (temporary arrays for changes in SAV)
      call init_svte(0)
    else if (mode.eq.1) then
!     store svte (the changes induced by the sampling) in svbu, zero out svte again
      call init_svte(1)
    end if
  end if
!
! let's only do all this setup work if necessary
  if ((use_IMPSOLV.EQV..true.).OR.(use_WCA.EQV..true.).OR.(use_IPP.EQV..true.).OR.&
 &    (use_attLJ.EQV..true.).OR.(use_CORR.EQV..true.).OR.(use_FEGS(1).EQV..true.).OR.&
 &    (use_FEGS(3).EQV..true.)) then
!
!   assemble list (makes code below easier)
    nrsl = 0
    do imol=1,nmol
      if (fmols(imol).eq.1) then
        if (bnd_type.eq.1) then
          do rs=rsmol(imol,1),rsmol(imol,2)
            nrsl = nrsl + 1
            rsl(nrsl) = rs
          end do
        end if
      else
        do rs=rsmol(imol,1),rsmol(imol,2)
          nrsl = nrsl + 1
          rsl(nrsl) = rs
        end do
      end if
    end do
!
    if ((cut.EQV..true.).AND.((use_rescrit.EQV..true.).OR.&
 &         (use_mcgrid.EQV..true.))) then
      if (use_mcgrid.EQV..true.) then
        do i=1,mmol 
          imol = imols(i)
          do rs=rsmol(imol,1),rsmol(imol,2)
            call grd_respairs(rs)
          end do
        end do
      else if (use_rescrit.EQV..true.) then
        do i=1,mmol 
          imol = imols(i)
          do rs=rsmol(imol,1),rsmol(imol,2)
            call respairs(rs)
          end do
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in clurb_energy_sh&
 &ort(...).'
        call fexit()
      end if
!
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nrsl,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,rs,sta,sto,incr,tpi,tpn,imol)
      tpn = omp_get_num_threads()
      tpi = omp_get_thread_num() + 1
      sta = tpi
      sto = nrsl
      incr = tpn
#else
      sta = 1
      sto = nrsl
      incr = 1
      tpi = 1
#endif
      do i=1,mmol 
        imol = imols(i)
        do rs=rsmol(imol,1),rsmol(imol,2)
          do j=sta,sto,incr
            k = rsl(j)
            if (fmols(molofrs(k)).eq.1) then
              if (imol.ge.molofrs(k)) cycle
            end if
            if (rsp_mat(max(rs,k),min(rs,k)).ne.0) then
              call Ven_rsp(evec_thr(:,tpi),min(rs,k),max(rs,k),use_cutoffs)
            end if
          end do
        end do
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
      do i=1,mmol 
        imol = imols(i)
        do rs=rsmol(imol,1),rsmol(imol,2)
          call clear_rsp2(rs)
        end do
      end do
    else
!
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nrsl,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,rs,sta,sto,incr,tpi,tpn,imol)
      tpn = omp_get_num_threads()
      tpi = omp_get_thread_num() + 1
      sta = tpi
      sto = nrsl
      incr = tpn
#else
      sta = 1
      sto = nrsl
      incr = 1
      tpi = 1
#endif
      do i=1,mmol 
        imol = imols(i)
        do rs=rsmol(imol,1),rsmol(imol,2)
          do j=sta,sto,incr
            k = rsl(j)
            if (fmols(molofrs(k)).eq.1) then
              if (imol.ge.molofrs(k)) cycle
            end if
            call Ven_rsp(evec_thr(:,tpi),min(rs,k),max(rs,k),use_cutoffs)
          end do
        end do
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
!
  else
!   do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.1) then
!     get the difference in SAV
      call init_svte(2)
!     now use that difference to determine which dependent interactions we have to re-compute
    end if
  end if
!
end
!
!------------------------------------------------------------------------
!
! a simple wrapper routine that handles long-range energy calculations for rigid body
! moves (only relevant terms are computed ...)
!
subroutine rigid_energy_long(imol,evec,cut,mode)
!
  use iounit
  use sequen
  use energies
  use molecule
  use cutoffs
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer imol,rs,j,mode,i,k,l,rsl(nseq),nrsl,sta(2),sto(2),incr,tpi
  RTYPE evec(MAXENERGYTERMS)
  logical cut,inrsl(nseq)
#ifdef ENABLE_THREADS
  integer tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
  maxtpn = 0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    call setup_scrqs2()
  end if
!
  nrsl = 0
!
  if ((use_IMPSOLV.EQV..true.).AND.((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.))) then
!
!   the problematic part is that the interaction of two residues can
!   indeed be changed by a third residue, which is changing conformation.
!   it is henceforth not safe to compute just the terms that change
!   by virtue of the conformational change.
!   to detect the implicitly affected terms we're using the temporary SAV
!   array (i.e., any atom whose SAV was affected by the conformational
!   change potentially interacts differently with all other atoms now)
!   through the rs_vec-array which is created during the last init_svte-call
!   in the short-range energy computation
!   note, however, that just like for every other interaction, the residue
!   based cutoff criterion still has to be fulfilled
!
    if ((cut.EQV..true.).AND.((use_rescrit.EQV..true.).OR.&
 &       (use_mcgrid.EQV..true.))) then
      inrsl(:) = .false.
      do rs=1,rsmol(imol,1)-1
        if (rs_vec(rs).eq.1) then ! Martin : that is if the sav has changed
          nrsl = nrsl + 1
          rsl(nrsl) = rs
          inrsl(rs) = .true.
        end if
      end do
      do rs=rsmol(imol,1),rsmol(imol,2)
        nrsl = nrsl + 1
        rsl(nrsl) = rs
        inrsl(rs) = .true.
       end do
      do rs=rsmol(imol,2)+1,nseq
        if (rs_vec(rs).eq.1) then
          nrsl = nrsl + 1
          rsl(nrsl) = rs
          inrsl(rs) = .true.
        end if
      end do
      if (use_mcgrid.EQV..true.) then
        do i=1,nrsl
          call grd_respairs(rsl(i))
        end do
      else if (use_rescrit.EQV..true.) then
        do i=1,nrsl
          call respairs(rsl(i))
        end do
      end if
!
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l,tpi,tpn,sto,sta,incr)
      tpn = omp_get_num_threads()
      tpi = omp_get_thread_num() + 1
      sta(1) = tpi
      sto(1) = nseq
      sto(2) = nrsl
      incr = tpn
#else
      sta(1) = 1
      incr = 1
      sto(1) = nseq
      sto(2) = nrsl
      tpi = 1
#endif
      do i=1,nrsl
        j = rsl(i)
        do k=sta(1),sto(1),incr !1,nseq
!         we take care of double-counting by using the inverse map inrsl
          if (inrsl(k).EQV..true.) cycle
!
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
        sta(2) = i + tpi - 1
!       and do the rsl-internal terms separately
        do l=sta(2),sto(2),incr !i,nrsl
          k = rsl(l)
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
      end do
!
      do i=sta(1),sto(1),incr ! 1,nseq
        if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
          call en_freesolv(evec_thr(:,tpi),i)
        end if
      end do
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
#endif
      do i=1,nrsl
        call clear_rsp2(rsl(i))
      end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
!
    else
!     this a disaster of course, use of cutoffs essential when doing IMPSOLV
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,k,tpi,tpn,sto,sta,incr)
      tpn = omp_get_num_threads()
      tpi = omp_get_thread_num() + 1
      sta(1) = tpi
      sto(1) = nseq
      incr = tpn
#else
      sta(1) = 1
      incr = 1
      sto(1) = nseq
      tpi = 1
#endif
      do j=sta(1),sto(1),incr !1,nseq
        do k=j,nseq
          call Ven_rsp_long(evec_thr(:,tpi),j,k,cut)
        end do
        call en_freesolv(evec_thr(:,tpi),j)
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
!
  else if (use_IMPSOLV.EQV..true.) then
#ifdef ENABLE_THREADS
    call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,tpi,tpn,sto,sta,incr)
    tpn = omp_get_num_threads()
    tpi = omp_get_thread_num() + 1
    sta(1) = tpi
    sto(1) = nseq
    incr = tpn
#else
    sta(1) = 1
    incr = 1
    sto(1) = nseq
    tpi = 1
#endif
    do j=sta(1),sto(1),incr
      call en_freesolv(evec_thr(:,tpi),j)
    end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
    if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
!
  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.)) then
    if ((cut.EQV..true.).AND.(use_rescrit.EQV..true.).OR.&
 &                      (use_mcgrid.EQV..true.)) then
      if (use_mcgrid.EQV..true.) then
        do rs=rsmol(imol,1),rsmol(imol,2)
          call grd_respairs(rs)
        end do
      else if (use_rescrit.EQV..true.) then
        do rs=rsmol(imol,1),rsmol(imol,2)
          call respairs(rs)
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in rigid_energy_&
 &long(...).'
        call fexit()
      end if
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rs,j,tpi,tpn,sto,sta,incr)
      tpn = omp_get_num_threads()
      tpi = omp_get_thread_num() + 1
      sta(1) = tpi
      sto(1) = nseq
      incr = tpn
#else
      sta(1) = 1
      incr = 1
      sto(1) = nseq
      tpi = 1
#endif
      do rs=rsmol(imol,1),rsmol(imol,2)
        do j=sta(1),sto(1),incr !1,rsmol(imol,1)-1 and rsmol(imol,2)+1,nseq
          if ((j.ge.rsmol(imol,1)).AND.(j.le.rsmol(imol,2))) cycle
          if (rsp_mat(min(j,rs),max(j,rs)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,rs),max(j,rs),use_cutoffs)
          end if
        end do
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
      do rs=rsmol(imol,1),rsmol(imol,2)
        call clear_rsp2(rs)
      end do
    else
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rs,j,tpi,tpn,sto,sta,incr)
      tpn = omp_get_num_threads()
      tpi = omp_get_thread_num() + 1
      sta(1) = tpi
      sto(1) = nseq
      incr = tpn
#else
      sta(1) = 1
      incr = 1
      sto(1) = nseq
      tpi = 1
#endif
      do rs=rsmol(imol,1),rsmol(imol,2)
        do j=sta(1),sto(1),incr !1,nseq
          if ((j.gt.rsmol(imol,2)).OR.(j.lt.rsmol(imol,1))) then
            call Ven_rsp_long(evec_thr(:,tpi),rs,j,cut)
          end if
        end do
      end do
#ifdef ENABLE_THREADS
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP BARRIER
!$OMP END PARALLEL
#endif
    end if
!
  else
!   do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     now update the atsav-array using svte (to be ready for posterior long-range calculation)
      call init_svte(5)
    else if (mode.eq.1) then
!     clear out the residue-svte-monitoring-vector
      do i=1,nseq
        rs_vec(i) = 0
        rs_vec2(i) = 0
      end do
    end if
  end if
!
end
!
!------------------------------------------------------------------------
!
! the same for clusters of molecules
!
subroutine clurb_energy_long(imols,mmol,fmols,evec,cut,mode)
!
  use iounit
  use sequen
  use energies
  use molecule
  use cutoffs
  use system
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer imol,rs,j,mode,i,k,l,rsl(nseq),nrsl,sta(2),sto(2),incr,tpi
  integer mmol,imols(nmol),fmols(nmol)
  RTYPE evec(MAXENERGYTERMS)
  logical cut,inrsl(nseq)
#ifdef ENABLE_THREADS
  integer tpn,maxtpn,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
  RTYPE evec_thr(MAXENERGYTERMS,thrdat%maxn)
!
  evec_thr(:,:) = 0.0
  maxtpn = 0
#else
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
#endif
!
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    call setup_scrqs2()
  end if
!
  nrsl = 0
!
  if ((use_IMPSOLV.EQV..true.).AND.((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.))) then
!
!   the problematic part is that the interaction of two residues can
!   indeed be changed by a third residue, which is changing conformation.
!   it is henceforth not safe to compute just the terms that change
!   by virtue of the conformational change.
!   to detect the implicitly affected terms we're using the temporary SAV
!   array (i.e., any atom whose SAV was affected by the conformational
!   change potentially interacts differently with all other atoms now)
!   through the rs_vec-array which is created during the last init_svte-call
!   in the short-range energy computation
!   note, however, that just like for every other interaction, the residue
!   based cutoff criterion still has to be fulfilled
!
    if ((cut.EQV..true.).AND.((use_rescrit.EQV..true.).OR.&
 &       (use_mcgrid.EQV..true.))) then
      inrsl(:) = .false.
      do imol=1,nmol
        if (fmols(imol).eq.1) then
          do rs=rsmol(imol,1),rsmol(imol,2)
            nrsl = nrsl + 1
            rsl(nrsl) = rs
            inrsl(rs) = .true.
          end do
        else
          do rs=rsmol(imol,1),rsmol(imol,2)
            if (rs_vec(rs).eq.1) then
              nrsl = nrsl + 1
              rsl(nrsl) = rs
              inrsl(rs) = .true.
            end if
          end do
        end if
      end do
      if (use_mcgrid.EQV..true.) then
        do i=1,nrsl
          call grd_respairs(rsl(i))
        end do
      else if (use_rescrit.EQV..true.) then
        do i=1,nrsl
          call respairs(rsl(i))
        end do
      end if
!
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l,sta,sto,incr,tpi,tpn)
      tpn = omp_get_num_threads()
      tpi = omp_get_thread_num() + 1
      sta(1) = tpi
      sto(1) = nseq
      sto(2) = nrsl
      incr = tpn
#else
      sta(1) = 1
      sto(1) = nseq
      sto(2) = nrsl
      incr = 1
      tpi = 1
#endif
      do i=1,nrsl
        j = rsl(i)
        do k=sta(1),sto(1),incr !1,nseq
!         in order to avoid double counting, we use the inverse map to skip out if both residues are in rsl
          if (inrsl(k).EQV..true.) cycle
!
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
!       and handle the interactions within rsl here
        sta(2) = i + tpi - 1
        do l=sta(2),sto(2),incr !i,nrsl
          k = rsl(l)
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
      end do
      do i=sta(1),sto(1),incr !1,nseq
        if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
          call en_freesolv(evec_thr(:,tpi),i)
        end if
      end do
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
      do i=1,nrsl
        call clear_rsp2(rsl(i))
      end do
!
    else ! meaning no cutoffs
!     this a disaster of course, use of cutoffs essential when doing IMPSOLV
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,k,sta,sto,incr,tpi,tpn)
      tpn = omp_get_num_threads()
      tpi = omp_get_thread_num() + 1
      sta(1) = tpi
      sto(1) = nseq
      incr = tpn
#else
      sta(1) = 1
      sto(1) = nseq
      incr = 1
      tpi = 1
#endif
      do j=sta(1),sto(1),incr !1,nseq
        do k=j,nseq
          call Ven_rsp_long(evec_thr(:,tpi),j,k,cut)
        end do
        call en_freesolv(evec_thr(:,tpi),j)
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL 
#endif
    end if
!
  else if (use_IMPSOLV.EQV..true.) then
#ifdef ENABLE_THREADS
    call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,sta,sto,incr,tpi,tpn)
    tpn = omp_get_num_threads()
    tpi = omp_get_thread_num() + 1
    sta(1) = tpi
    sto(1) = nseq
    incr = tpn
#else
    sta(1) = 1
    sto(1) = nseq
    incr = 1
    tpi = 1
#endif
    do i=sta(1),sto(1),incr !1,nseq
      if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
        call en_freesolv(evec_thr(:,tpi),i)
      end if
    end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL 
#endif
!
  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.)) then
    if ((cut.EQV..true.).AND.(use_rescrit.EQV..true.).OR.&
 &                      (use_mcgrid.EQV..true.)) then
      if (use_mcgrid.EQV..true.) then
        do i=1,mmol 
          imol = imols(i)
          do rs=rsmol(imol,1),rsmol(imol,2)
            call grd_respairs(rs)
          end do
        end do
      else if (use_rescrit.EQV..true.) then
        do i=1,mmol 
          imol = imols(i)
          do rs=rsmol(imol,1),rsmol(imol,2)
            call respairs(rs)
          end do
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in clurb_energy_&
 &long(...).'
        call fexit()
      end if
!     we'll re-interpret and re-utilize rsl here
!     this includes the potential pitfall of having to recompute internal cluster energies
!     if PBC are used
      nrsl = 0
      do imol=1,nmol
        if (fmols(imol).eq.1) then
          if (bnd_type.eq.1) then
            do rs=rsmol(imol,1),rsmol(imol,2)
              nrsl = nrsl + 1
              rsl(nrsl) = rs
            end do
          end if
        else
          do rs=rsmol(imol,1),rsmol(imol,2)
            nrsl = nrsl + 1
            rsl(nrsl) = rs
          end do
        end if
      end do
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nrsl,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,rs,sta,sto,incr,tpi,tpn,imol)
      tpn = omp_get_num_threads()
      tpi = omp_get_thread_num() + 1
      sta(1) = tpi
      sto(1) = nrsl
      incr = tpn
#else
      sta(1) = 1
      sto(1) = nrsl
      incr = 1
      tpi = 1
#endif
      do i=1,mmol 
        imol = imols(i)
        do rs=rsmol(imol,1),rsmol(imol,2)
          do j=sta(1),sto(1),incr
            k = rsl(j)
            if (fmols(molofrs(k)).eq.1) then
              if (imol.ge.molofrs(k)) cycle
            end if
            if (rsp_mat(min(rs,k),max(rs,k)).ne.0) then
              call Ven_rsp_long(evec_thr(:,tpi),min(rs,k),max(rs,k),use_cutoffs)
            end if
          end do
        end do
      end do
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
      do i=1,mmol 
        imol = imols(i)
        do rs=rsmol(imol,1),rsmol(imol,2)
          call clear_rsp2(rs)
        end do
      end do
    else
!     ditto
      nrsl = 0
      do imol=1,nmol
        if (fmols(imol).eq.1) then
          if (bnd_type.eq.1) then
            do rs=rsmol(imol,1),rsmol(imol,2)
              nrsl = nrsl + 1
              rsl(nrsl) = rs
            end do
          end if
        else
          do rs=rsmol(imol,1),rsmol(imol,2)
            nrsl = nrsl + 1
            rsl(nrsl) = rs
          end do
        end if
      end do
#ifdef ENABLE_THREADS
      call omp_set_num_threads(min(nrsl,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,rs,sta,sto,incr,tpi,tpn,imol)
      tpn = omp_get_num_threads()
      tpi = omp_get_thread_num() + 1
      sta(1) = tpi
      sto(1) = nrsl
      incr = tpn
#else
      sta(1) = 1
      sto(1) = nrsl
      incr = 1
      tpi = 1
#endif
      do i=1,mmol 
        imol = imols(i)
        do rs=rsmol(imol,1),rsmol(imol,2)
          do j=sta(1),sto(1),incr
            k = rsl(j)
            if (fmols(molofrs(k)).eq.1) then
              if (imol.ge.molofrs(k)) cycle
            end if
            call Ven_rsp_long(evec_thr(:,tpi),min(rs,k),max(rs,k),use_cutoffs)
          end do
        end do
      end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      if (tpn.gt.maxtpn) maxtpn = tpn
!$OMP END SINGLE
!$OMP END PARALLEL
#endif
    end if
!
  else
!   do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
#ifdef ENABLE_THREADS
  do j=2,maxtpn
    evec(:) = evec(:) + evec_thr(:,j)
  end do
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     now update the atsav-array using svte (to be ready for posterior long-range calculation)
      call init_svte(5)
    else if (mode.eq.1) then
!     clear out the residue-svte-monitoring-vector
      do i=1,nseq
        rs_vec(i) = 0
        rs_vec2(i) = 0
      end do
    end if
  end if
!
end
!
!------------------------------------------------------------------------
!
!             #################################
!             #                               #
!             # THIRD: OTHER ROUTINES         #
!             #                               #
!             #################################
!
!------------------------------------------------------------------------
!
! this subroutine is a wrapper mimicking steric/SAV interactions between crosslinked
! residues as a neighbor relation -> currently this relies on all energy calculations
! using Vforce_rsp and Vforce_rsp_long but no neighbor list-based routines (see force_wrap.f90)
!
subroutine Ven_crosslink_corr(lk,evec,cut)
!
  use sequen
  use inter
  use energies
  use atoms
  use params
  use polypep
!
  implicit none
!
  integer i,rs1,rs2,bun(4),lk,ii,jj,kk,mm
  RTYPE evec(MAXENERGYTERMS)
  integer bulst1(nrsnb(crosslink(lk)%rsnrs(1)),2)
  RTYPE bulst3(nrsnb(crosslink(lk)%rsnrs(1)),3)
  logical cut,foundit,bul
!
  if (ideal_run.EQV..true.) then
    return
  end if
!
  rs1 = crosslink(lk)%rsnrs(1)
  rs2 = crosslink(lk)%rsnrs(2)
  if (rs2.lt.rs1) then
    rs2 = crosslink(lk)%rsnrs(1)
    rs1 = crosslink(lk)%rsnrs(2)
  end if
  bun(1) = nrsnb(rs1)
  bun(2) = at(rs1+1)%na
  bun(3) = molofrs(rs1+1)
  bun(4) = refat(rs1+1)
  if (use_FEG.EQV..true.) then
    bul = par_FEG(rs1+1)
  end if
  bulst1(1:nrsnb(rs1),1) = iaa(rs1)%atnb(1:nrsnb(rs1),1)
  bulst1(1:nrsnb(rs1),2) = iaa(rs1)%atnb(1:nrsnb(rs1),2)
  bulst3(1:nrsnb(rs1),1) = fudge(rs1)%rsnb_ljs(1:nrsnb(rs1))
  bulst3(1:nrsnb(rs1),2) = fudge(rs1)%rsnb_lje(1:nrsnb(rs1))
  bulst3(1:nrsnb(rs1),3) = fudge(rs1)%rsnb(1:nrsnb(rs1))
! overwrite
  nrsnb(rs1) = 0
  do ii=at(rs1)%bb(1),at(rs1)%bb(1)+at(rs1)%na-1
    do kk=at(rs2)%bb(1),at(rs2)%bb(1)+at(rs2)%na-1
      foundit = .false.
      do i=1,crosslink(lk)%nrsin
        jj = crosslink(lk)%exclin(i,1)
        mm = crosslink(lk)%exclin(i,2)
        if (((ii.eq.jj).AND.(kk.eq.mm)).OR.&
 &          ((ii.eq.mm).AND.(kk.eq.jj))) then
          if (crosslink(lk)%is14in(i).EQV..true.) then
            nrsnb(rs1) = nrsnb(rs1) + 1
            iaa(rs1)%atnb(nrsnb(rs1),1) = ii
            iaa(rs1)%atnb(nrsnb(rs1),2) = kk
            fudge(rs1)%rsnb_ljs(nrsnb(rs1)) = lj_sig_14(attyp(ii),attyp(kk))
            fudge(rs1)%rsnb_lje(nrsnb(rs1)) = lj_eps_14(attyp(ii),attyp(kk))
            fudge(rs1)%rsnb(nrsnb(rs1)) = fudge_st_14
            foundit = .true.
          else
            foundit = .true.
          end if
       end if
        if (foundit.EQV..true.) exit
      end do
      if (foundit.EQV..false.) then
        nrsnb(rs1) = nrsnb(rs1) + 1
        iaa(rs1)%atnb(nrsnb(rs1),1) = ii
        iaa(rs1)%atnb(nrsnb(rs1),2) = kk
        fudge(rs1)%rsnb_ljs(nrsnb(rs1)) = lj_sig(attyp(ii),attyp(kk))
        fudge(rs1)%rsnb_lje(nrsnb(rs1)) = lj_eps(attyp(ii),attyp(kk))
        fudge(rs1)%rsnb(nrsnb(rs1)) = 1.0
      end if
    end do
  end do
  at(rs1+1)%na = at(rs2)%na
  molofrs(rs1+1) = molofrs(rs2)
  refat(rs1+1) = refat(rs2)
  if (use_FEG.EQV..true.) then
    par_FEG(rs1+1) = par_FEG(rs2)
  end if
! now call the appropriate routine
  call Ven_rsp(evec,rs1,rs1+1,cut)
! and restore what we destroyed
  nrsnb(rs1) = bun(1)
  at(rs1+1)%na = bun(2)
  molofrs(rs1+1) = bun(3)
  refat(rs1+1) = bun(4)
  if (use_FEG.EQV..true.) then
    par_FEG(rs1+1) = bul
  end if
  fudge(rs1)%rsnb_ljs(1:nrsnb(rs1)) = bulst3(1:nrsnb(rs1),1)
  fudge(rs1)%rsnb_lje(1:nrsnb(rs1)) = bulst3(1:nrsnb(rs1),2)
  fudge(rs1)%rsnb(1:nrsnb(rs1)) = bulst3(1:nrsnb(rs1),3)
  iaa(rs1)%atnb(1:nrsnb(rs1),1) = bulst1(1:nrsnb(rs1),1) 
  iaa(rs1)%atnb(1:nrsnb(rs1),2) = bulst1(1:nrsnb(rs1),2)
!
end
!
!-----------------------------------------------------------------------
!
! this subroutine is a wrapper mimicking polar/tab. interactions between crosslinked
! residues as a neighbor relation; this also supports Ewald/GRF
!
subroutine Ven_crosslink_corr_long(lk,evec,cut)
!
  use sequen
  use inter
  use energies
  use atoms
  use params
  use polypep
  use cutoffs
!
  implicit none
!
  integer i,rs1,rs2,bun(4),lk,ii,jj,kk,iii,kkk,mm
  RTYPE evec(MAXENERGYTERMS)
  integer bulst1(nrsnb(crosslink(lk)%rsnrs(1)),2)
  RTYPE bulst3(nrsnb(crosslink(lk)%rsnrs(1)))
  logical cut,foundit,bul,bul2
!
  if ((ideal_run.EQV..true.).OR.&
 &    ((use_POLAR.EQV..false.).AND.(use_TABUL.EQV..false.))) then
    return
  end if
!
  rs1 = crosslink(lk)%rsnrs(1)
  rs2 = crosslink(lk)%rsnrs(2)
  if (rs2.lt.rs1) then
    rs2 = crosslink(lk)%rsnrs(1)
    rs1 = crosslink(lk)%rsnrs(2)
  end if
  if (use_POLAR.EQV..false.) then
!   for just the tabulated potential does not require any trickery (all interactions computed
!   irrespective of topology)
    call Ven_rsp_long(evec,rs1,rs2,cut)
    return
  end if
  bun(1) = nrpolnb(rs1)
  bun(2) = at(rs1+1)%na
  bun(3) = molofrs(rs1+1)
  bun(4) = refat(rs1+1)
  if (use_FEG.EQV..true.) then
    bul = par_FEG(rs1+1)
  end if
  bul2 = use_TABUL
  bulst1(1:nrpolnb(rs1),1) = iaa(rs1)%polnb(1:nrpolnb(rs1),1)
  bulst1(1:nrpolnb(rs1),2) = iaa(rs1)%polnb(1:nrpolnb(rs1),2)
  bulst3(1:nrpolnb(rs1)) = fudge(rs1)%elnb(1:nrpolnb(rs1))
  nrpolnb(rs1) = 0
  do iii=1,at(rs1)%npol
    ii = at(rs1)%pol(iii)
    do kkk=1,at(rs2)%npol
      kk = at(rs2)%pol(kkk)
      foundit = .false.
      do i=1,crosslink(lk)%nrspol
        jj = crosslink(lk)%exclpol(i,1)
        mm = crosslink(lk)%exclpol(i,2)
        if (((ii.eq.jj).AND.(kk.eq.mm)).OR.&
 &          ((ii.eq.mm).AND.(kk.eq.jj))) then
          if (crosslink(lk)%is14pol(i).EQV..true.) then
            nrpolnb(rs1) = nrpolnb(rs1) + 1
            iaa(rs1)%polnb(nrpolnb(rs1),1) = ii
            iaa(rs1)%polnb(nrpolnb(rs1),2) = kk
            fudge(rs1)%elnb(nrpolnb(rs1)) = fudge_el_14
            foundit = .true.
          else
            foundit = .true.
          end if
        end if
        if (foundit.EQV..true.) exit
      end do
      if (foundit.EQV..false.) then
        nrpolnb(rs1) = nrpolnb(rs1) + 1
        iaa(rs1)%polnb(nrpolnb(rs1),1) = ii
        iaa(rs1)%polnb(nrpolnb(rs1),2) = kk
        fudge(rs1)%elnb(nrpolnb(rs1)) = 1.0
      end if
    end do
  end do
  at(rs1+1)%na = at(rs2)%na
  molofrs(rs1+1) = molofrs(rs2)
  refat(rs1+1) = refat(rs2)
  if (use_FEG.EQV..true.) then
    par_FEG(rs1+1) = par_FEG(rs2)
  end if
  use_TABUL = .false. ! see below
! now call the appropriate routine
  call Ven_rsp_long(evec,rs1,rs1+1,cut)
! and restore what we destroyed
  nrpolnb(rs1) = bun(1)
  at(rs1+1)%na = bun(2)
  molofrs(rs1+1) = bun(3)
  refat(rs1+1) = bun(4)
  if (use_FEG.EQV..true.) then
    par_FEG(rs1+1) = bul
  end if
  use_TABUL = bul2
  fudge(rs1)%elnb(1:nrpolnb(rs1)) = bulst3(1:nrpolnb(rs1))
  iaa(rs1)%polnb(1:nrpolnb(rs1),1) = bulst1(1:nrpolnb(rs1),1) 
  iaa(rs1)%polnb(1:nrpolnb(rs1),2) = bulst1(1:nrpolnb(rs1),2)
  if (use_TABUL.EQV..true.) then
    bul2 = use_POLAR
    use_POLAR = .false.
    call Ven_rsp_long(evec,rs1,rs2,cut)
    use_POLAR = bul2
  end if
!
end
!
!-------------------------------------------------------------------------------
!
! this simple wrapper routine sets the system's Hamiltonian to different
! conditions allow swap moves in REMC ...
! note that lamenergy is implictily threaded (if so compiled) through the call to the
! main time-consuming function energy3
!
subroutine lamenergy(rve,fve)
!
  use iounit
  use system
  use energies
  use units
  use mpistuff
  use molecule
  use dssps
  use ems
!
!  martin debug 
!  use atoms 
!  use sequen
  
  
  implicit none
!
  RTYPE rve(mpi_nodes),fve(mpi_nodes)
#ifdef ENABLE_MPI
  integer i,j,k,imol,which,i_start,i_end
  RTYPE energy3,evec(MAXENERGYTERMS)
  RTYPE vbu(MAXREDIMS),eva,dum,dum2
  logical needsavup,needemup,badflg(MAXREDIMS),atrue,afalse
  integer vbui(MAXREDIMS)
  logical needpka
  
  
  
  !Martin added
  logical ee
  integer rs1,rs2
  
  
  
  
#endif
!
#ifdef ENABLE_MPI
!
  needsavup = .false.
  needemup = .false.
  needpka = .false.
  atrue = .true.
  afalse = .false.
  dum2 = 0.0
!  
  
  do i=1,re_conddim
    if (re_types(i).eq.1) vbu(i)=kelvin
    if (re_types(i).eq.2) vbu(i)=scale_IPP
    if (re_types(i).eq.3) vbu(i)=scale_attLJ
    if (re_types(i).eq.4) vbu(i)=scale_WCA
    if (re_types(i).eq.5) vbu(i)=scale_POLAR
    if (re_types(i).eq.6) vbu(i)=scale_IMPSOLV
    if (re_types(i).eq.7) vbu(i)=par_IMPSOLV(2)
    if (re_types(i).eq.8) vbu(i)=scale_TOR
    if (re_types(i).eq.9) vbu(i)=scale_ZSEC
    if (re_types(i).eq.10) vbu(i)=par_ZSEC(1)
    if (re_types(i).eq.11) vbu(i)=par_ZSEC(3)
    if (re_types(i).eq.12) vbui(i)=scrq_model
    if (re_types(i).eq.13) vbu(i)=par_IMPSOLV(3)
    if (re_types(i).eq.14) vbu(i)=par_IMPSOLV(4)
    if (re_types(i).eq.15) vbu(i)=par_IMPSOLV(6)
    if (re_types(i).eq.16) vbu(i)=par_IMPSOLV(7)
    if (re_types(i).eq.17) vbu(i)=par_IMPSOLV(8)
    if (re_types(i).eq.18) vbui(i)=i_sqm
    if (re_types(i).eq.19) vbu(i)=par_IMPSOLV(9)
    if (re_types(i).eq.20) vbu(i)=scale_FEGS(1)
    if (re_types(i).eq.21) vbu(i)=scale_FEGS(3)
    if (re_types(i).eq.22) vbu(i)=scale_FEGS(6)
    if (re_types(i).eq.23) vbu(i)=scale_TABUL
    if (re_types(i).eq.24) vbu(i)=scale_POLY
    if (re_types(i).eq.25) vbu(i)=scale_DREST
    if (re_types(i).eq.26) vbu(i)=scale_FEGS(15)
    if (re_types(i).eq.27) vbu(i)=scale_FEGS(16)
    if (re_types(i).eq.28) vbu(i)=scale_FEGS(17)
    if (re_types(i).eq.29) vbu(i)=scale_FEGS(18)
    if (re_types(i).eq.30) vbu(i)=par_DSSP(9)
    if (re_types(i).eq.31) vbu(i)=par_DSSP(7)
!   do nothing for 32
    if (re_types(i).eq.33) vbu(i)=scale_EMICRO
    if (re_types(i).eq.34) vbu(i)=emthreshdensity
    if (re_types(i).eq.35) vbu(i)=par_pka(1)
    if (re_types(i).eq.36) vbu(i)=par_pka(2)
    if (re_types(i).eq.37) vbu(i)=par_pka(3)! Martin : added
    if (re_types(i).eq.38) vbu(i)=par_pka(4)! Martin : added again
    if (re_types(i).eq.39) vbu(i)=scale_FEGS(4)! Martin : added
  end do
!
! note that MPI_REMaster(...) temporarily changes reol_all for the
! the actual exchange overlap calculations if need be (i.e., if re_nbmode is 1)
  if (reol_all.EQV..true.) then
    i_start = 1
    i_end = re_conditions
  else
    i_start = max((myrank+1)-1,1)
    i_end = min((myrank+1)+1,re_conditions)
  end if
!
  do i=i_start,i_end
!   set the different conditions
!   remember that we don't change the use_XX-flags
!   this implies, however, that for scale_XX = 0.0, we compute
!   a bunch of terms all multiplied by 0.0. in order to preserve
!   this information (to compute derivatives with respect to scale_XX),
!   we therefore use a little detour (set to 1.0, subtract out)
    do j=1,re_conddim
      badflg(j) = .false.
      if (re_types(j).eq.1) then
        kelvin = re_mat(i,j)
        invtemp = 1.0/(gasconst*kelvin)
      else if (re_types(j).eq.2) then
        scale_IPP   = re_mat(i,j)
        if (scale_IPP.le.0.0) then
          badflg(j) = .true.
          scale_IPP = 1.0
        end if
      else if (re_types(j).eq.3) then
        scale_attLJ = re_mat(i,j)
        if (scale_attLJ.le.0.0) then
          badflg(j) = .true.
          scale_attLJ = 1.0
        end if
      else if (re_types(j).eq.4) then
        scale_WCA   = re_mat(i,j)
        if (scale_WCA.le.0.0) then
          badflg(j) = .true.
          scale_WCA = 1.0
        end if
      else if (re_types(j).eq.5) then
        scale_POLAR = re_mat(i,j)
        if (scale_POLAR.le.0.0) then
          badflg(j) = .true.
          scale_POLAR = 1.0
        end if
      else if (re_types(j).eq.6) then
        scale_IMPSOLV = re_mat(i,j)
        if (scale_IMPSOLV.le.0.0) then
          badflg(j) = .true.
          scale_IMPSOLV = 1.0
        end if
      else if (re_types(j).eq.7) then
!       note that the relevance of this relies entirely on whether use_IMPSOLV is true
        par_IMPSOLV(2) = re_mat(i,j)
        if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
          coul_scr = 1.0 - 1.0/par_IMPSOLV(2)
        else
          coul_scr = 1.0 - 1.0/sqrt(par_IMPSOLV(2))
        end if
      else if (re_types(j).eq.8) then
        scale_TOR = re_mat(i,j)
        if (scale_TOR.le.0.0) then
          badflg(j) = .true.
          scale_TOR = 1.0
        end if
      else if (re_types(j).eq.9) then
        scale_ZSEC = re_mat(i,j)
        if (scale_ZSEC.le.0.0) then
          badflg(j) = .true.
          scale_ZSEC = 1.0
        end if
      else if (re_types(j).eq.10) then
!       note that the relevance of this hinges on scale_ZSEC
        par_ZSEC(1) = re_mat(i,j)
      else if (re_types(j).eq.11) then
!       note that the relevance of this hinges on scale_ZSEC
        par_ZSEC(3) = re_mat(i,j)
      else if (re_types(j).eq.12) then
!       note that the relevance of this hinges on scale_IMPSOLV
        scrq_model = nint(re_mat(i,j))
        if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
          coul_scr = 1.0 - 1.0/par_IMPSOLV(2)
        else
          coul_scr = 1.0 - 1.0/sqrt(par_IMPSOLV(2))
        end if
      else if (re_types(j).eq.13) then
!       note that the relevance of this hinges on scale_IMPSOLV
        needsavup = .true.
        par_IMPSOLV(3) = re_mat(i,j)
      else if (re_types(j).eq.14) then
!       note that the relevance of this hinges on scale_IMPSOLV
        needsavup = .true.
        par_IMPSOLV(4) = re_mat(i,j)
      else if (re_types(j).eq.15) then
!       note that the relevance of this hinges on scale_IMPSOLV
        needsavup = .true.
        par_IMPSOLV(6) = re_mat(i,j)
      else if (re_types(j).eq.16) then
!       note that the relevance of this hinges on scale_IMPSOLV
        needsavup = .true.
        par_IMPSOLV(7) = re_mat(i,j)
      else if (re_types(j).eq.17) then
!       note that the relevance of this hinges on scale_IMPSOLV and scrq_model
        par_IMPSOLV(8) = 1./re_mat(i,j)
      else if (re_types(j).eq.18) then
!       note that the relevance of this hinges on scale_IMPSOLV and scrq_model
        i_sqm = nint(re_mat(i,j))
      else if (re_types(j).eq.19) then
!       note that the relevance of this hinges on scale_IMPSOLV and scrq_model
        par_IMPSOLV(9) = re_mat(i,j)
      else if (re_types(j).eq.20) then
        scale_FEGS(1) = re_mat(i,j)
        call setup_parFEG(1)
!        if (scale_FEGS(1).le.0.0) then
!          badflg(j) = .true.
!          scale_FEGS(1) = 1.0
!        end if
      else if (re_types(j).eq.21) then
        scale_FEGS(3) = re_mat(i,j)
        call setup_parFEG(3)
!        if (scale_FEGS(3).le.0.0) then
!          badflg(j) = .true.
!          scale_FEGS(3) = 1.0
!        end if
      else if (re_types(j).eq.22) then
        scale_FEGS(6) = re_mat(i,j)
        call setup_parFEG(6)
!        if (scale_FEGS(6).le.0.0) then
!          badflg(j) = .true.
!          scale_FEGS(6) = 1.0
!        end if
      else if (re_types(j).eq.23) then
        scale_TABUL = re_mat(i,j)
        if (scale_TABUL.le.0.0) then
          badflg(j) = .true.
          scale_TABUL = 1.0
        end if
      else if (re_types(j).eq.24) then
        scale_POLY = re_mat(i,j)
        if (scale_POLY.le.0.0) then
          badflg(j) = .true.
          scale_POLY = 1.0
        end if
      else if (re_types(j).eq.25) then
        scale_DREST = re_mat(i,j)
        if (scale_DREST.le.0.0) then
          badflg(j) = .true.
          scale_DREST = 1.0
        end if
      else if (re_types(j).eq.26) then
        scale_FEGS(15) = re_mat(i,j)
      else if (re_types(j).eq.27) then
        scale_FEGS(16) = re_mat(i,j)
      else if (re_types(j).eq.28) then
        scale_FEGS(17) = re_mat(i,j)
      else if (re_types(j).eq.29) then
        scale_FEGS(18) = re_mat(i,j)
      else if (re_types(j).eq.30) then
!       note that the relevance of this hinges on scale_DSSP
        par_DSSP(9) = re_mat(i,j)
      else if (re_types(j).eq.31) then
!       note that the relevance of this hinges on scale_DSSP
        par_DSSP(7) = re_mat(i,j)
!     do nothing for 32
      else if (re_types(j).eq.33) then
        scale_EMICRO = re_mat(i,j)
        if (scale_EMICRO.le.0.0) then
          badflg(j) = .true.
          scale_EMICRO = 1.0
        end if
      else if (re_types(j).eq.34) then
        emthreshdensity = re_mat(i,j)
        needemup = .true.
      else if (re_types(j).eq.35) then
        par_pka(1) = re_mat(i,j)
        needpka = .true.
      else if (re_types(j).eq.36) then
        par_pka(2) = re_mat(i,j)
        !needsavup=.true. ! Martin added
        needpka = .true.
      else if (re_types(j).eq.37) then
        par_pka(3) = re_mat(i,j)
        needpka = .true.
      else if (re_types(j).eq.38) then
        par_pka(4) = re_mat(i,j)
        needpka = .true.
      else if (re_types(j).eq.39) then!Martin added for Moses
        scale_FEGS(4) = re_mat(i,j)
      end if
!
    end do
    
    if (needsavup.EQV..true.) call absinth_savprm()
    
    
    if (needpka.eqv..true.) then !Martin added
      call scale_stuff_pka()
    end if 

    if (needemup.EQV..true.) then
      call scale_emmap(emthreshdensity,dum,dum2)
      call precompute_diff_emmap()
    end if
!   now compute the current energy at this condition with
!   our structure
    if (i.eq.i_start) then
      eva = invtemp*energy3(evec,atrue)
    else
      eva = invtemp*energy3(evec,afalse)
    end if
!   the computation of derivatives is really only meaningful when
!   noTI is false
    fve(i) = 0.0
    do j=1,re_conddim ! goes through the number of replica exchange conditions, in this part at least no reference to the pk titration is being made
      if (badflg(j).EQV..true.) then
        if (re_types(j).eq.2) then
          fve(i) = fve(i) + evec(1)
          evec(1) = 0.0
        else if (re_types(j).eq.3) then
          fve(i) = fve(i) + evec(3)
          evec(3) = 0.0
        else if (re_types(j).eq.4) then
          fve(i) = fve(i) + evec(5)
          evec(5) = 0.0
        else if (re_types(j).eq.5) then
          fve(i) = fve(i) + evec(6)
          evec(6) = 0.0
        else if (re_types(j).eq.6) then
          fve(i) = fve(i) + evec(4)
          evec(4) = 0.0
        else if (re_types(j).eq.8) then
          fve(i) = fve(i) + evec(7)
          evec(7) = 0.0
        else if (re_types(j).eq.9) then
          fve(i) = fve(i) + evec(8)
          evec(8) = 0.0
        else if (re_types(j).eq.23) then
          fve(i) = fve(i) + evec(9)
          evec(9) = 0.0
        else if (re_types(j).eq.24) then
          fve(i) = fve(i) + evec(14)
          evec(14) = 0.0
        else if (re_types(j).eq.25) then
          fve(i) = fve(i) + evec(10)
          evec(10) = 0.0
        else if (re_types(j).eq.33) then
          fve(i) = fve(i) + evec(21)
          evec(21) = 0.0
        end if
      else
        if (re_types(j).eq.2) then
          fve(i) = fve(i) + evec(1)/scale_IPP
        else if (re_types(j).eq.3) then
          fve(i) = fve(i) + evec(3)/scale_attLJ
        else if (re_types(j).eq.4) then
          fve(i) = fve(i) + evec(5)/scale_WCA
        else if (re_types(j).eq.5) then
          fve(i) = fve(i) + evec(6)/scale_POLAR
        else if (re_types(j).eq.6) then
          fve(i) = fve(i) + evec(4)/scale_IMPSOLV
        else if (re_types(j).eq.8) then
          fve(i) = fve(i) + evec(7)/scale_TOR
        else if (re_types(j).eq.9) then
          fve(i) = fve(i) + evec(8)/scale_ZSEC
        else if (re_types(j).eq.10) then
          which = 1
          do imol=1,nmol
            call der_zsec_gl(imol,fve(i),which)
          end do
        else if (re_types(j).eq.11) then
          which = 2
          do imol=1,nmol
            call der_zsec_gl(imol,fve(i),which)
          end do
        else if (re_types(j).eq.23) then
          fve(i) = fve(i) + evec(9)/scale_TABUL
        else if (re_types(j).eq.24) then
          fve(i) = fve(i) + evec(14)/scale_POLY
        else if (re_types(j).eq.25) then
          fve(i) = fve(i) + evec(10)/scale_DREST
!        else if (re_types(j).eq.30) then
!          which = 2
!          do imol=1,nmol
!            call der_dssp_gl(imol,fve(i),which)
!          end do
!        else if (re_types(j).eq.31) then
!          which = 1
!          do imol=1,nmol
!            call der_dssp_gl(imol,fve(i),which)
!          end do
        else if (re_types(j).eq.33) then
          fve(i) = fve(i) + evec(21)/scale_EMICRO
        end if
      end if
    end do
    
!   now recover the actual energy at this condition (see above)
    rve(i) = 0.0
    do k=1,MAXENERGYTERMS
      rve(i) = rve(i) + evec(k)
    end do
    rve(i) = rve(i)*invtemp
    needsavup = .false.
    needemup = .false.
    needpka = .false.
  end do
  

!
  do i=1,re_conddim
    if (re_types(i).eq.1) then 
      kelvin = vbu(i)
      invtemp = 1.0/(gasconst*kelvin)
    else if (re_types(i).eq.2) then
      scale_IPP   = vbu(i)
    else if (re_types(i).eq.3) then
      scale_attLJ = vbu(i)
    else if (re_types(i).eq.4) then
      scale_WCA   = vbu(i)
    else if (re_types(i).eq.5) then
      scale_POLAR = vbu(i)
    else if (re_types(i).eq.6) then
      scale_IMPSOLV = vbu(i)
    else if (re_types(i).eq.7) then
      par_IMPSOLV(2) = vbu(i)
      if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
        coul_scr = 1.0 - 1.0/par_IMPSOLV(2)
      else
        coul_scr = 1.0 - 1.0/sqrt(par_IMPSOLV(2))
      end if
    else if (re_types(i).eq.8) then
      scale_TOR = vbu(i)
    else if (re_types(i).eq.9) then
      scale_ZSEC = vbu(i)
    else if (re_types(i).eq.10) then
      par_ZSEC(1) = vbu(i)
    else if (re_types(i).eq.11) then
      par_ZSEC(3) = vbu(i)
    else if (re_types(i).eq.12) then
      scrq_model = vbui(i)
      if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
        coul_scr = 1.0 - 1.0/par_IMPSOLV(2)
      else
        coul_scr = 1.0 - 1.0/sqrt(par_IMPSOLV(2))
      end if
    else if (re_types(i).eq.13) then
      needsavup = .true.
      par_IMPSOLV(3) = vbu(i)
    else if (re_types(i).eq.14) then
      needsavup = .true.
      par_IMPSOLV(4) = vbu(i)
    else if (re_types(i).eq.15) then
      needsavup = .true.
      par_IMPSOLV(6) = vbu(i)
    else if (re_types(i).eq.16) then
      needsavup = .true.
      par_IMPSOLV(7) = vbu(i)
    else if (re_types(i).eq.17) then
      par_IMPSOLV(8) = vbu(i)
    else if (re_types(i).eq.18) then
      i_sqm = vbui(i)
    else if (re_types(i).eq.19) then
      par_IMPSOLV(9) = vbu(i)
    else if (re_types(i).eq.20) then
      scale_FEGS(1) = vbu(i)
      call setup_parFEG(1) 
    else if (re_types(i).eq.21) then
      scale_FEGS(3) = vbu(i)
      call setup_parFEG(3)
    else if (re_types(i).eq.22) then
      scale_FEGS(6) = vbu(i)
      call setup_parFEG(6)
    else if (re_types(i).eq.23) then
      scale_TABUL = vbu(i)
    else if (re_types(i).eq.24) then
      scale_POLY  = vbu(i)
    else if (re_types(i).eq.25) then
      scale_DREST = vbu(i)
    else if (re_types(i).eq.26) then
      scale_FEGS(15) = vbu(i)
    else if (re_types(i).eq.27) then
      scale_FEGS(16) = vbu(i)
    else if (re_types(i).eq.28) then
      scale_FEGS(17) = vbu(i)
    else if (re_types(i).eq.29) then
      scale_FEGS(18) = vbu(i)
    else if (re_types(i).eq.30) then
      par_DSSP(9) = vbu(i)
    else if (re_types(i).eq.31) then
      par_DSSP(7) = vbu(i)
!   do nothing for 32
    else if (re_types(i).eq.33) then
      scale_EMICRO = vbu(i)
    else if (re_types(i).eq.34) then
      emthreshdensity = vbu(i)
      needemup = .true.
    else if (re_types(i).eq.35) then
      par_pka(1) = vbu(i)
      needpka = .true.
    else if (re_types(i).eq.36) then
      par_pka(2) = vbu(i)
      !needsavup= .true. ! Martin : added becasue the parameters do need to be recomputed
      needpka = .true.
    else if (re_types(i).eq.37) then
      par_pka(3) = vbu(i)
      needpka = .true.
    else if (re_types(i).eq.38) then
      par_pka(4) = vbu(i)
      needpka = .true.
    else if (re_types(i).eq.39) then ! Martin added for moses
      scale_FEGS(4) = vbu(i)
    end if
  end do
  
    if (needpka.eqv..true.) then !Martin added
      call scale_stuff_pka()
      eva = invtemp*energy3(evec,afalse)
    end if 
  
  if (needsavup.EQV..true.) call absinth_savprm()
  if (needemup.EQV..true.) then
    call scale_emmap(emthreshdensity,dum,dum2)
    call precompute_diff_emmap()
  end if
  
#else
!
  write(ilog,*) 'Fatal. Called lamenergy(...) in non MPI-calculation&
 &. This is most definitely a bug.'
  call fexit()
#endif
end
!
!-----------------------------------------------------------------------------
!
subroutine setup_parFEG(which)
!
  use energies
  use iounit 
!
  implicit none
!
  integer which

  if (which.eq.1) then
    if (scale_FEGS(1).ge.1.0) then
      par_FEG2(1) = scale_FEGS(1)
      par_FEG2(2) = 0.0
    else
      if (fegljmode.eq.1) then
        par_FEG2(1) = scale_FEGS(1)
        par_FEG2(2) = 0.0
      else if (fegljmode.eq.2) then
        par_FEG2(1) = scale_FEGS(1)**par_FEG2(8)
        par_FEG2(2) = par_FEG2(3)*(1.0-scale_FEGS(1)**par_FEG2(7))
      else if (fegljmode.eq.3) then
        par_FEG2(1) = (1.0-exp(-par_FEG2(8)*scale_FEGS(1)))/&
 &                                 (1.0-exp(-par_FEG2(8)))
        par_FEG2(2) = par_FEG2(3)*(1.0-scale_FEGS(1))**par_FEG2(7)
      else
        write(ilog,*) 'Fatal. Got bad scaling mode for FEG-LJ (',&
 &fegljmode,').'
        call fexit()
      end if
    end if
  else if (which.eq.3) then
    if (scale_FEGS(3).ge.1.0) then
      par_FEG2(5) = scale_FEGS(3)
      par_FEG2(6) = 0.0
    else
      if (fegljmode.eq.1) then
        par_FEG2(5) = scale_FEGS(3)
        par_FEG2(6) = 0.0
      else if (fegljmode.eq.2) then
        par_FEG2(5) = scale_FEGS(3)**par_FEG2(8)
        par_FEG2(6) = par_FEG2(3)*(1.0-scale_FEGS(3)**par_FEG2(7))
      else if (fegljmode.eq.3) then
        par_FEG2(5) = (1.0-exp(-par_FEG2(8)*scale_FEGS(3)))/&
 &                                 (1.0-exp(-par_FEG2(8)))
        par_FEG2(6) = par_FEG2(3)*(1.0-scale_FEGS(3))**par_FEG2(7)
      else
        write(ilog,*) 'Fatal. Got bad scaling mode for FEG-LJ (',&
 &fegljmode,').'
        call fexit()
      end if
    end if
    ! Don't actually need, going t o use scale directly 
!  else if (which.eq.4) then ! Martin added for moses
!    if (scale_FEGS(4).ge.1.0) then
!      par_FEG2(4) = scale_FEGS(4)
!    else
!      if (fegljmode.eq.1) then
!        par_FEG2(4) = scale_FEGS(4)
!      else
!        write(ilog,*) 'Fatal. Got bad scaling mode for FEG-FOS (',&
! &fegljmode,').'
!        call fexit()
!      end if
!    end if
    
  else if (which.eq.6) then
    if (scale_FEGS(6).ge.1.0) then
      par_FEG2(9) = scale_FEGS(6)
      par_FEG2(10) = 0.0
    else
      if (fegcbmode.eq.1) then
        par_FEG2(9) = scale_FEGS(6)
        par_FEG2(10) = 0.0
      else if (fegcbmode.eq.2) then
        par_FEG2(9) = scale_FEGS(6)**par_FEG2(12)
        par_FEG2(10) = par_FEG2(4)*(1.0-scale_FEGS(3)**par_FEG2(11))
      else
        write(ilog,*) 'Fatal. Got bad scaling mode for FEG-Cb (',&
 &fegcbmode,').'
        call fexit()
      end if
    end if
  else
    write(ilog,*) 'Fatal. Called setup_parFEG(...) for unsupported p&
 &tential term (got ',which,').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
! a simple wrapper to summarize torsional potentials applied to individual
! backbone torsions (handled through function en_torrs)
! currently NOT IN USE
!
function etorgl()
!
  use sequen
  use energies
!
  implicit none
!
  integer i
  RTYPE etorgl,ev(MAXENERGYTERMS)
!
  etorgl = 0.0d0
!
  do i = 1,nseq
    if (par_TOR2(i).gt.0) then
      call en_torrs(ev,i)
      etorgl = etorgl + ev(7)
    end if
  end do
!
end
!
!------------------------------------------------------------------------------------
!
#ifdef ENABLE_THREADS
!
! with threaded loop execution it is of course difficult to set perfect bounds if the "dimensions" of the
! loop continuously change at runtime
! this routine attempts to set per-thread bounds smartly using simple heuristics
!
  subroutine get_thread_loop_bounds(mode,rsi,rsf,cut,tpn,stas,stos)
!
  use threads
  use sequen
  use cutoffs
  use polypep
!
  implicit none
!
  integer i,j,rsi,rsf,mode,totia,stas(thrdat%maxn),stos(thrdat%maxn),tpn
  logical cut
!
! short-range chi
  if (mode.eq.1) then
    if (nseq.le.thrdat%maxn) then
      tpn = nseq
!       stas(1:
      return
    end if
    if (cut.EQV..true.) then
      do j=1,rsi
        if (rsp_mat(rsi,j).ne.0) totia = totia + at(rsi)%na*at(j)%na
      end do
      do j=rsi+1,nseq
        if (rsp_mat(j,rsi).ne.0) totia = totia + at(rsi)%na*at(j)%na
      end do
    end if
  else
  end if
!
  end
!
#endif

