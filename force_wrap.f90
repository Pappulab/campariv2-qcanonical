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
!             #        FOR CARTESIAN FORCES   #
!             #        AND GLOBAL ENERGIES    #
!             #                               #
!             #################################
!
!
!------------------------------------------------------------------------
!
function force1(evec)
!
  use sequen
  use energies
  use forces
  use molecule
  use atoms
  use fyoc
  use system
  use iounit
  use polypep
  use wl
!
  implicit none
!
  integer i,j,imol,rs,aone,atwo
  RTYPE evec(MAXENERGYTERMS),force1,edum
  logical sayno,sayyes,pflag
!
  aone = 1
  atwo = 2
  evec(:) = 0.0
  bnd_f(:,:) = 0.0
  bnd_fr = 0.0
  sayno = .false.
  sayyes = .true.
! for later use: pass on to inner_loops fxns (virial is non-trivial for PBC!)
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    pflag = .true.
  else
    pflag = .false.
  end if
  ens%insVirT(:,:) = 0.0
  cart_f(:,:) = 0.0
  if (do_accelsim.EQV..true.) hmjam%ca_f(:,:) = 0.0
  if (use_IMPSOLV.EQV..true.) then
    sisa_nr(:) = 0
    sum_scrcbs(:) = 0.0
    sum_scrcbs_tr(:) = 0.0
    sum_scrcbs_lr(:) = 0.0
    sav_dr(:,:) = 0.0
    scrq_dr(:,:) = 0.0
    sisa_dr(:,:,:) = 0.0
    sisq_dr(:,:,:) = 0.0
  end if
!
  if (use_IMPSOLV.EQV..true.) call init_svte(3)
!
  do i=1,nseq
    do j=i,nseq
      if (disulf(i).eq.j) cycle
      call Vforce_rsp(evec,i,j,sayno,cart_f)
    end do
  end do
! crosslinks will always be treated as short-range 
  do i=1,n_crosslinks
    call Vforce_crosslink_corr(i,evec,sayno,cart_f)
  end do
  if ((do_accelsim.EQV..true.).AND.(hmjam%isin(1).EQV..true.)) then ! only IPP/ATTLJ/WCA in cart_f* so far
    hmjam%ca_f(:,:) = cart_f(:,:)
    i = 1
    edum = evec(1)+evec(3)+evec(5)
    call els_manage_one(i,edum,cart_f,atwo)
    if (use_WCA.EQV..true.) then
      evec(5) = evec(5) + hmjam%boosts(1,1)
    else if (use_IPP.EQV..true.) then
      evec(1) = evec(1) + hmjam%boosts(1,1)
    else 
      evec(3) = evec(3) + hmjam%boosts(1,1)
    end if
  end if
!
  if (use_IMPSOLV.EQV..true.) then
    call init_svte(6)
    if (use_POLAR.EQV..true.) then
      call force_setup_scrqs2()
    end if
  end if
!
  if ((do_accelsim.EQV..true.).AND.(hmjam%isin(6).EQV..true.)) then ! joint TABUL/POLAR
    do i=1,nseq
      do j=i,nseq
        if (disulf(i).eq.j) cycle
        call Vforce_rsp_long(evec,i,j,sayno,hmjam%ca_f,sum_scrcbs)
      end do
    end do
!   crosslinks will always be treated as short-range 
    do i=1,n_crosslinks
      call Vforce_crosslink_corr_long(i,evec,sayno,hmjam%ca_f,sum_scrcbs)
    end do
!
    if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
      do i=1,nseq
        call Vforce_dcbdscrq(i,hmjam%ca_f)
      end do
    end if
    edum = evec(6) + evec(9)
    i = 6
    call els_manage_one(i,edum,cart_f,aone)
    if (use_POLAR.EQV..true.) then
      evec(6) = evec(6) + hmjam%boosts(6,1)
    else
      evec(9) = evec(9) + hmjam%boosts(6,1)
    end if
  else
    do i=1,nseq
      do j=i,nseq
        if (disulf(i).eq.j) cycle
        call Vforce_rsp_long(evec,i,j,sayno,cart_f,sum_scrcbs)
      end do
    end do
!   crosslinks will always be treated as short-range 
    do i=1,n_crosslinks
      call Vforce_crosslink_corr_long(i,evec,sayno,cart_f,sum_scrcbs)
    end do
!
    if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
      do i=1,nseq
        call Vforce_dcbdscrq(i,cart_f)
      end do
    end if
  end if
!
  if (use_IMPSOLV.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(4).EQV..true.)) then
      do i=1,nseq
        call force_freesolv(i,evec,hmjam%ca_f)
      end do
      i = 4
      call els_manage_one(i,evec(4),cart_f,aone)
    else
      do i=1,nseq
        call force_freesolv(i,evec,cart_f)
      end do
    end if
  end if
!
  if (use_BOND(1).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(15).EQV..true.)) then
      do rs=1,nseq
        call force_bonds(rs,evec,hmjam%ca_f)
      end do
      i = 15
      call els_manage_one(i,evec(15),cart_f,aone)
    else
      do rs=1,nseq
        call force_bonds(rs,evec,cart_f)
      end do
    end if
  end if
  if (use_BOND(2).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(16).EQV..true.)) then
      do rs=1,nseq
        call force_angles(rs,evec,hmjam%ca_f)
      end do
      i = 16
      call els_manage_one(i,evec(16),cart_f,aone)
    else
      do rs=1,nseq
        call force_angles(rs,evec,cart_f)
      end do
    end if
  end if
  if (use_BOND(3).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(17).EQV..true.)) then
      do rs=1,nseq
        call force_impropers(rs,evec,hmjam%ca_f)
      end do
      i = 17
      call els_manage_one(i,evec(17),cart_f,aone)
    else
      do rs=1,nseq
        call force_impropers(rs,evec,cart_f)
      end do
    end if
  end if
  if (use_BOND(4).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(18).EQV..true.)) then
      do rs=1,nseq
        call force_torsions(rs,evec,hmjam%ca_f)
      end do
      i = 18
      call els_manage_one(i,evec(18),cart_f,aone)
    else
      do rs=1,nseq
        call force_torsions(rs,evec,cart_f)
      end do
    end if
  end if
  if (use_BOND(5).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(20).EQV..true.)) then
      do rs=1,nseq
        call force_cmap(rs,evec,hmjam%ca_f)
      end do
      i = 20
      call els_manage_one(i,evec(20),cart_f,aone)
    else
      do rs=1,nseq
        call force_cmap(rs,evec,cart_f)
      end do
    end if
  end if
  if (use_CORR.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(2).EQV..true.)) then
      do rs=1,nseq
        call force_corrector(rs,evec,hmjam%ca_f)
      end do
      i = 2
      call els_manage_one(i,evec(2),cart_f,aone)
    else
      do rs=1,nseq
        call force_corrector(rs,evec,cart_f)
      end do
    end if
  end if
!
  if (use_TOR.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(7).EQV..true.)) then
      do i=1,nseq
        if (par_TOR2(i).gt.0) then
          call force_torrs(evec,i,hmjam%ca_f)
        end if
      end do
      i = 7
      call els_manage_one(i,evec(7),cart_f,aone)
    else
      do i=1,nseq
        if (par_TOR2(i).gt.0) then
          call force_torrs(evec,i,cart_f)
        end if
      end do
    end if
  end if
!
  if ((do_accelsim.EQV..true.).AND.(hmjam%isin(12).EQV..true.)) then
    do rs=1,nseq
      call force_boundary_rs(rs,evec,aone,hmjam%ca_f)
    end do
    i = 12
    call els_manage_one(i,evec(12),cart_f,aone)
  else
    do rs=1,nseq
      call force_boundary_rs(rs,evec,aone,cart_f)
    end do
  end if
!
  if (use_ZSEC.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(8).EQV..true.)) then
      do imol=1,nmol
        call force_zsec_gl(imol,evec,hmjam%ca_f)
      end do
      i = 8
      call els_manage_one(i,evec(8),cart_f,aone)
    else
      do imol=1,nmol
        call force_zsec_gl(imol,evec,cart_f)
      end do
    end if
  end if
!
  if (use_DSSP.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(19).EQV..true.)) then
      call force_dssp_gl(evec,hmjam%ca_f)
      i = 19
      call els_manage_one(i,evec(19),cart_f,aone)
    else
      call force_dssp_gl(evec,cart_f)
    end if
  end if
!
  if (use_EMICRO.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(21).EQV..true.)) then
      call force_emicro_gl(evec,hmjam%ca_f)
      i = 21
      call els_manage_one(i,evec(21),cart_f,aone)
    else
      call force_emicro_gl(evec,cart_f)
    end if
  end if
!
  if (use_POLY.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(14).EQV..true.)) then
      do imol=1,nmol
        call force_poly_gl(imol,evec,hmjam%ca_f)
      end do
      i = 14
      call els_manage_one(i,evec(14),cart_f,aone)
    else
      do imol=1,nmol
        call force_poly_gl(imol,evec,cart_f)
      end do
    end if
  end if
!
  if (use_DREST.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(10).EQV..true.)) then
      call force_drest(evec,hmjam%ca_f)
      i = 10
      call els_manage_one(i,evec(10),cart_f,aone)
    else
      call force_drest(evec,cart_f)
    end if
  end if
!
  if (use_IMPSOLV.EQV..true.) then
    do i=1,n
      if (sisa_nr(i).gt.MAXNB) then
        write(ilog,*) 'WARNING. Maximum number of neighbor atoms exc&
 &eeded for atom ',i,'. It is very likely this run will crash soon. &
 &Ignore any of the results obtained.'
         call fexit()
      end if
    end do
  end if
!
  force1 = sum(evec(:))
!
  if ((do_accelsim.EQV..true.).AND.(hmjam%isin(11).EQV..true.)) then
    hmjam%ca_f(:,:) = cart_f(:,:)
    edum = force1
    i = 11
    call els_manage_one(i,edum,cart_f,atwo)
    evec(11) = hmjam%boosts(11,1)
    force1 = force1 + evec(11)
  end if
!
end
!
!-----------------------------------------------------------------------
!
function force3(evec,evec_tr,evec_lr,forceflag)
!
  use sequen
  use energies
  use forces
  use molecule
  use atoms
  use iounit
  use cutoffs
  use mcsums
  use inter
  use polypep
  use units
  use fyoc
  use system
  use movesets
#ifdef ENABLE_MPI
  use mpistuff
#endif
  use wl
!
  implicit none
!
#ifdef ENABLE_MPIBLOED
#include "mpif.h"
#endif
  integer i,j,imol,rs,tpi,aone,atwo
  RTYPE evec(MAXENERGYTERMS),force3
  RTYPE evec_tr(MAXENERGYTERMS),evec_lr(MAXENERGYTERMS)
  logical sayno,sayyes,forceflag,freshnbl,pflag,elsflag
  RTYPE evec_thr(MAXENERGYTERMS,3,2),edum
#ifdef ENABLE_MPIBLOED
  integer mstatus(MPI_STATUS_SIZE),ierr,Reqlist(mpi_nodes),mstlst(MPI_STATUS_SIZE,mpi_nodes),req
  RTYPE somex(n)
#endif
!
#ifdef ENABLE_THREADS
  RTYPE force3_threads
  force3 = force3_threads(evec,evec_tr,evec_lr,forceflag)
  return
#endif
!
#ifdef ENABLE_MPIBLOED
  if (use_MPIAVG.EQV..true.) then
!    call CPU_time(t3)
    k = 44
    if (myrank.eq.0) then
      do i=1,mpi_nodes-1
        call MPI_ISEND(x,n,MPI_RTYPE,i,k,MPI_COMM_WORLD,Reqlist(i+1),ierr)
      end do
!      call MPI_WAITALL(mpi_nodes-1,Reqlist(2:mpi_nodes),mstlst(:,2:mpi_nodes),ierr)
    else
      call MPI_IRECV(somex,n,MPI_RTYPE,0,MPI_ANY_TAG,MPI_COMM_WORLD,req,ierr)
      call MPI_WAIT(req,mstatus,ierr)
    end if
!    call CPU_time(t4)
!    write(*,*) myrank,': ',t4-t3
  end if
#endif
!
! first determine whether we need to update the NB-list
  freshnbl = .false.
  if ((mod(nstep,nbl_up).eq.0).OR.(nstep.le.1).OR.&
 &    (forceflag.EQV..true.)) then
    freshnbl = .true.
  end if
!
! init
  aone = 1
  atwo = 2
  elsflag = .false.
  evec_thr(:,:,1) = 0.0
  tpi = 1
  evec(:) = 0.0
  bnd_fr = 0.0
  sayno = .false.
  sayyes = .true.
! for later use: pass on to inner_loops fxns (virial is non-trivial for PBC!)
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    pflag = .true.
  else
    pflag = .false.
  end if
  ens%insVirT(:,:) = 0.0
  bnd_f(:,:) = 0.0
  cart_f(:,:) = 0.0
! the mid- and long-range interaction arrays (energies and cartesian forces) are only updated
! every so many steps
  if ((freshnbl.EQV..true.)) then
    evec_tr(:) = 0.0
    cart_f_tr(:,:) = 0.0
    if (use_IMPSOLV.EQV..true.) sum_scrcbs_tr(:) = 0.0
  else
    if (lrel_md.eq.2) then
      cart_f_lr(:,:) = 0.0
      evec_lr(:) = 0.0
    end if
  end if
  if ((freshnbl.EQV..true.)) then
    evec_lr(:) = 0.0
    cart_f_lr(:,:) = 0.0
    if (use_IMPSOLV.EQV..true.) sum_scrcbs_lr(:) = 0.0
  end if
  if (do_accelsim.EQV..true.) hmjam%ca_f(:,:) = 0.0
  if (use_IMPSOLV.EQV..true.) then
    sisa_nr(:) = 0
    sum_scrcbs(:) = 0.0
    sav_dr(:,:) = 0.0
    scrq_dr(:,:) = 0.0
    sisa_dr(:,:,:) = 0.0
    sisq_dr(:,:,:) = 0.0
  end if
!
!
  if (is_plj.EQV..true.) then
!
!   first calculations which use the standard NB gas phase potential (Coulomb+LJ)
!   note that many additional potentials can still be handled including EXTRA (CORR)
!   these routines are strictly incompatible with IMPSOLV, TABUL, and WCA, but will
!   also seg-fault when called WITHOUT the Coulomb term
!
    if (use_cutoffs.EQV..true.) then
!
      if (freshnbl.EQV..true.) then
        call all_respairs_nbl(skip_frz,sayno,tpi)
      end if
!
      if (use_waterloops.EQV..true.) then
        do i=1,rsw1-1
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
                call Vforce_rsp_long(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          call Vforce_rsp_long(evec_thr(:,1,tpi),i,i,sayno,cart_f,sum_scrcbs)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
              call Vforce_rsp_long(evec_thr(:,1,tpi),i,i+1,sayno,cart_f,sum_scrcbs)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
        end do
        do i=rsw1,nseq
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
                call Vforce_rsp_long(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          call Vforce_rsp_long(evec_thr(:,1,tpi),i,i,sayno,cart_f,sum_scrcbs)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
              call Vforce_rsp_long(evec_thr(:,1,tpi),i,i+1,sayno,cart_f,sum_scrcbs)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            if ((seqtyp(rsw1).eq.45).OR.(seqtyp(rsw1).eq.103)) then
              call Vforce_PLJ_C_T4P(evec_thr(:,1,tpi),i,cart_f)
            else
              call Vforce_PLJ_C_T3P(evec_thr(:,1,tpi),i,cart_f)
            end if
          end if
        end do
!       add long-range electrostatics corrections
        call force_P_LR(evec_thr(:,3,tpi),cart_f_lr,sum_scrcbs_lr,freshnbl)
      else
        do i=1,nseq
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
                call Vforce_rsp_long(evec_thr(:,2,tpi),i,i+1,use_cutoffs,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          call Vforce_rsp_long(evec_thr(:,1,tpi),i,i,sayno,cart_f,sum_scrcbs)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
              call Vforce_rsp_long(evec_thr(:,1,tpi),i,i+1,use_cutoffs,cart_f,sum_scrcbs)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
        end do
!       add long-range electrostatics corrections
        call force_P_LR(evec_thr(:,3,tpi),cart_f_lr,sum_scrcbs_lr,freshnbl)
      end if
!
    else
!     no cutoffs and no simplifying assumptions
      do i=1,nseq
        call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
        call Vforce_rsp_long(evec_thr(:,1,tpi),i,i,sayno,cart_f,sum_scrcbs)
        if (i.lt.nseq) then
          call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
          call Vforce_rsp_long(evec_thr(:,1,tpi),i,i+1,sayno,cart_f,sum_scrcbs)
        end if
        if (i.lt.(nseq-1)) then
          call Vforce_PLJ(evec_thr(:,1,tpi),i,cart_f)
        end if
      end do
    end if
!   crosslinks will always be treated as short-range 
    if (n_crosslinks.gt.0) call all_crosslink_corrs(evec_thr(:,1,tpi),sayno,cart_f,sum_scrcbs)
!
  else if (is_pewlj.EQV..true.) then
!
!   this is the modified tree with all routines being able to handle Ewald sums
!
    if (use_cutoffs.EQV..true.) then
!     note that we enforce during parsing that there is no twin-regime, so we
!     don't have to check it here.
!     also note that this is potentially a weak assumption for LJ interactions,
!     but that Ewald is generally safer if more density is in the real-space 
!     component using a slightly larger SR cutoff anyway
      if (freshnbl.EQV..true.) then
        call all_respairs_nbl(skip_frz,sayno,tpi)
      end if
!
      if (use_waterloops.EQV..true.) then
        do i=1,rsw1-1
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          call force_rsp_long_mod(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
              call force_rsp_long_mod(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PEWLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
        end do
        do i=rsw1,nseq
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          call force_rsp_long_mod(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
              call force_rsp_long_mod(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            if ((seqtyp(rsw1).eq.45).OR.(seqtyp(rsw1).eq.103)) then
              call Vforce_PEWLJ_C_T4P(evec_thr(:,1,tpi),i,cart_f)
            else
              call Vforce_PEWLJ_C_T3P(evec_thr(:,1,tpi),i,cart_f)
            end if
          end if
        end do
!       add reciprocal sum/forces and constant term to energy
        call force_P_LR(evec_thr(:,3,tpi),cart_f_lr,sum_scrcbs_lr,freshnbl)
      else
        do i=1,nseq
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          call force_rsp_long_mod(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
              call force_rsp_long_mod(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PEWLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
        end do
!       add reciprocal sum/forces and constant term to energy
        call force_P_LR(evec_thr(:,3,tpi),cart_f_lr,sum_scrcbs_lr,freshnbl)
      end if
!     crosslinks will always be treated as short-range 
      if (n_crosslinks.gt.0) call all_crosslink_corrs(evec_thr(:,1,tpi),sayno,cart_f,sum_scrcbs)
!
    else
!     Ewald must have a real-space cutoff
      write(ilog,*) 'Fatal. Ewald summation requires the use of a cutoff.'
      call fexit()
    end if
!
  else if (is_prflj.EQV..true.) then
!
!   this is the modified tree with all routines being able to handle (G)RF electrostatics
!
    if (use_cutoffs.EQV..true.) then
!
      if (freshnbl.EQV..true.) then
        call all_respairs_nbl(skip_frz,sayno,tpi)
      end if
!
      if (use_waterloops.EQV..true.) then
        do i=1,rsw1-1
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
                call force_rsp_long_mod(evec_thr(:,2,tpi),i,i+1,use_cutoffs,cart_f_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PRFLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          call force_rsp_long_mod(evec_thr(:,1,tpi),i,i,use_cutoffs,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
              call force_rsp_long_mod(evec_thr(:,1,tpi),i,i+1,use_cutoffs,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PRFLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
        end do
        do i=rsw1,nseq
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
                call force_rsp_long_mod(evec_thr(:,2,tpi),i,i+1,use_cutoffs,cart_f_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PRFLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          call force_rsp_long_mod(evec_thr(:,1,tpi),i,i,use_cutoffs,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
              call force_rsp_long_mod(evec_thr(:,1,tpi),i,i+1,use_cutoffs,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            if ((seqtyp(rsw1).eq.45).OR.(seqtyp(rsw1).eq.103)) then
              call Vforce_PRFLJ_C_T4P(evec_thr(:,1,tpi),i,cart_f)
            else
              call Vforce_PRFLJ_C_T3P(evec_thr(:,1,tpi),i,cart_f)
            end if
          end if
        end do
!       add long-range electrostatics corrections
        call force_P_LR(evec_thr(:,3,tpi),cart_f_lr,sum_scrcbs_lr,freshnbl)
      else
        do i=1,nseq
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
                call force_rsp_long_mod(evec_thr(:,2,tpi),i,i+1,use_cutoffs,cart_f_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PRFLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          call force_rsp_long_mod(evec_thr(:,1,tpi),i,i,use_cutoffs,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
              call force_rsp_long_mod(evec_thr(:,1,tpi),i,i+1,use_cutoffs,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PRFLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
        end do
!       add long-range electrostatics corrections
        call force_P_LR(evec_thr(:,3,tpi),cart_f_lr,sum_scrcbs_lr,freshnbl)
      end if
!     crosslinks will always be treated as short-range 
      if (n_crosslinks.gt.0) call all_crosslink_corrs(evec_thr(:,1,tpi),sayno,cart_f,sum_scrcbs)
!
    else
!     RF must have a real-space cutoff
      write(ilog,*) 'Fatal. (G)RF corrections require the use of a cutoff.'
      call fexit()
    end if
!
  else if (is_fegplj.EQV..true.) then
!
!   this is the modified tree with all routines being able to handle ghosted Hamiltonians
!   essentially, separate neighbor lists are set up for ghosted interactions
!
    if (use_cutoffs.EQV..true.) then
!
      if (freshnbl.EQV..true.) then
        call all_respairs_nbl(skip_frz,sayno,tpi)
      end if
!
      if (use_waterloops.EQV..true.) then
        do i=1,rsw1-1
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
                call Vforce_rsp_long(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
            if (rs_nbl(i)%ngnbtrs.gt.0) then
              call Vforce_FEGPLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          call Vforce_rsp_long(evec_thr(:,1,tpi),i,i,sayno,cart_f,sum_scrcbs)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
              call Vforce_rsp_long(evec_thr(:,1,tpi),i,i+1,sayno,cart_f,sum_scrcbs)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
          if (rs_nbl(i)%ngnbs.gt.0) then
            call Vforce_FEGPLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
        end do
        do i=rsw1,nseq
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
                call Vforce_rsp_long(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
            if (rs_nbl(i)%ngnbtrs.gt.0) then
              call Vforce_FEGPLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          call Vforce_rsp_long(evec_thr(:,1,tpi),i,i,sayno,cart_f,sum_scrcbs)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
              call Vforce_rsp_long(evec_thr(:,1,tpi),i,i+1,sayno,cart_f,sum_scrcbs)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            if ((seqtyp(rsw1).eq.45).OR.(seqtyp(rsw1).eq.103)) then
              call Vforce_PLJ_C_T4P(evec_thr(:,1,tpi),i,cart_f)
            else
              call Vforce_PLJ_C_T3P(evec_thr(:,1,tpi),i,cart_f)
            end if
          end if
          if (rs_nbl(i)%ngnbs.gt.0) then
            call Vforce_FEGPLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
        end do
!       add long-range electrostatics corrections
        call force_P_LR(evec_thr(:,3,tpi),cart_f_lr,sum_scrcbs_lr,freshnbl)
      else
        do i=1,nseq
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
                call Vforce_rsp_long(evec_thr(:,2,tpi),i,i+1,use_cutoffs,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
            if (rs_nbl(i)%ngnbtrs.gt.0) then
              call Vforce_FEGPLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          call Vforce_rsp_long(evec_thr(:,1,tpi),i,i,sayno,cart_f,sum_scrcbs)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
              call Vforce_rsp_long(evec_thr(:,1,tpi),i,i+1,sayno,cart_f,sum_scrcbs)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
          if (rs_nbl(i)%ngnbs.gt.0) then
            call Vforce_FEGPLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
        end do
!       add long-range electrostatics corrections
        call force_P_LR(evec_thr(:,3,tpi),cart_f_lr,sum_scrcbs_lr,freshnbl)
      end if
!
    else
!     no cutoffs, no simplifying assumptions, and no speed-up
      do i=1,nseq
        do j=i,nseq
          if (disulf(i).eq.j) cycle
          call Vforce_rsp(evec_thr(:,1,tpi),i,j,sayno,cart_f)
          call Vforce_rsp_long(evec_thr(:,1,tpi),i,j,sayno,cart_f,sum_scrcbs)
        end do
      end do
    end if
!   crosslinks will always be treated as short-range 
    if (n_crosslinks.gt.0) call all_crosslink_corrs(evec_thr(:,1,tpi),sayno,cart_f,sum_scrcbs)
!
  else if (is_fegprflj.EQV..true.) then
!
!   this is the modified tree with all routines being able to handle ghosted Hamiltonians + RF
!   for simple Cb-scaling (linear) only
!
    if (use_cutoffs.EQV..true.) then
!
      if (freshnbl.EQV..true.) then
        call all_respairs_nbl(skip_frz,sayno,tpi)
      end if
!
      if (use_waterloops.EQV..true.) then
        do i=1,rsw1-1
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
                call force_rsp_long_mod(evec_thr(:,2,tpi),i,i+1,use_cutoffs,cart_f_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PRFLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
            if (rs_nbl(i)%ngnbtrs.gt.0) then
              call Vforce_FEGPRFLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          call force_rsp_long_mod(evec_thr(:,1,tpi),i,i,use_cutoffs,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
              call force_rsp_long_mod(evec_thr(:,1,tpi),i,i+1,use_cutoffs,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PRFLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
          if (rs_nbl(i)%ngnbs.gt.0) then
            call Vforce_FEGPRFLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
        end do
        do i=rsw1,nseq
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
                call force_rsp_long_mod(evec_thr(:,2,tpi),i,i+1,use_cutoffs,cart_f_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PRFLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
            if (rs_nbl(i)%ngnbtrs.gt.0) then
              call Vforce_FEGPRFLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          call force_rsp_long_mod(evec_thr(:,1,tpi),i,i,use_cutoffs,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
              call force_rsp_long_mod(evec_thr(:,1,tpi),i,i+1,use_cutoffs,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            if ((seqtyp(rsw1).eq.45).OR.(seqtyp(rsw1).eq.103)) then
              call Vforce_PRFLJ_C_T4P(evec_thr(:,1,tpi),i,cart_f)
            else
              call Vforce_PRFLJ_C_T3P(evec_thr(:,1,tpi),i,cart_f)
            end if
          end if
          if (rs_nbl(i)%ngnbs.gt.0) then
            call Vforce_FEGPRFLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
        end do
!       add long-range electrostatics corrections
        call force_P_LR(evec_thr(:,3,tpi),cart_f_lr,sum_scrcbs_lr,freshnbl)
      else
        do i=1,nseq
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
                call force_rsp_long_mod(evec_thr(:,2,tpi),i,i+1,use_cutoffs,cart_f_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PRFLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
            if (rs_nbl(i)%ngnbtrs.gt.0) then
              call Vforce_FEGPRFLJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          call force_rsp_long_mod(evec_thr(:,1,tpi),i,i,use_cutoffs,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
              call force_rsp_long_mod(evec_thr(:,1,tpi),i,i+1,use_cutoffs,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PRFLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
          if (rs_nbl(i)%ngnbs.gt.0) then
            call Vforce_FEGPRFLJ_C(evec_thr(:,1,tpi),i,cart_f)
          end if
        end do
!       add long-range electrostatics corrections
        call force_P_LR(evec_thr(:,3,tpi),cart_f_lr,sum_scrcbs_lr,freshnbl)
      end if
!     crosslinks will always be treated as short-range 
      if (n_crosslinks.gt.0) call all_crosslink_corrs(evec_thr(:,1,tpi),sayno,cart_f,sum_scrcbs)
    else
!     RF must have a real-space cutoff
      write(ilog,*) 'Fatal. (G)RF corrections require the use of a cutoff.'
      call fexit()
    end if
!
  else if ((is_lj.EQV..true.).OR.(is_ev.EQV..true.)) then
!
    if (use_cutoffs.EQV..true.) then
!
      if (freshnbl.EQV..true.) then
        call all_respairs_nbl(skip_frz,sayno,tpi)
      end if
!
      do i=1,nseq
        if ((freshnbl.EQV..true.)) then
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.2) then
              call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
            end if
          end if
          if (rs_nbl(i)%nnbtrs.gt.0) then
            if (is_lj.EQV..true.) call Vforce_LJ_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            if (is_ev.EQV..true.) call Vforce_IPP_TR(evec_thr(:,2,tpi),i,cart_f_tr)
          end if
        end if
        call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
        if (i.lt.nseq) then
          if (rsp_vec(i).eq.1) then
            call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
          end if
        end if
        if (rs_nbl(i)%nnbs.gt.0) then
          if (is_lj.EQV..true.) call Vforce_LJ_C(evec_thr(:,1,tpi),i,cart_f)
          if (is_ev.EQV..true.) call Vforce_IPP_C(evec_thr(:,1,tpi),i,cart_f)
        end if
      end do
!
    else
!     no cutoffs and no simplifying assumptions
      do i=1,nseq
        call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f) !cart_f)
        if (i.lt.nseq) then
          call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f) !cart_f)
        end if
        if (i.lt.(nseq-1)) then
          if (is_lj.EQV..true.) call Vforce_LJ(evec_thr(:,1,tpi),i,cart_f)
          if (is_ev.EQV..true.) call Vforce_IPP(evec_thr(:,1,tpi),i,cart_f)
        end if
      end do
    end if
!   crosslinks will always be treated as short-range 
    do i=1,n_crosslinks
      call Vforce_crosslink_corr(i,evec_thr(:,1,tpi),sayno,cart_f)
    end do
!
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(1).EQV..true.)) then ! only IPP/ATTLJ in cart_f* so far
      hmjam%ca_f(:,:) = cart_f(:,:) + cart_f_tr(:,:)
      i = 1
      edum = evec_thr(1,1,tpi) + evec_thr(1,2,tpi) + evec_thr(3,1,tpi) + evec_thr(3,2,tpi)
      call els_manage_one(i,edum,cart_f,atwo)
      evec_thr(1,1,tpi) = evec_thr(1,1,tpi) + hmjam%boosts(1,1)
    end if
!
  else if (is_tab.EQV..true.) then
! 
!   a modified tree designed to handle calculations with just tabulated potentials
!   including very sparse assignments
!
    if (use_cutoffs.EQV..true.) then
!
      if (freshnbl.EQV..true.) then
        call all_respairs_nbl(skip_frz,sayyes,tpi)
      end if
!
      do i=1,nseq
        if ((freshnbl.EQV..true.)) then
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.2) then
!              call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
              call Vforce_rsp_long(evec_thr(:,2,tpi),i,i+1,use_cutoffs,cart_f_tr,sum_scrcbs_tr)
            end if
          end if
          if (rs_nbl(i)%nnbtrs.gt.0) then
            call Vforce_TAB_TR(evec_thr(:,2,tpi),i,cart_f_tr)
          end if
        end if
!        call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
        call Vforce_rsp_long(evec_thr(:,1,tpi),i,i,sayno,cart_f,sum_scrcbs)
        if (i.lt.nseq) then
          if (rsp_vec(i).eq.1) then
!            call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
            call Vforce_rsp_long(evec_thr(:,1,tpi),i,i+1,use_cutoffs,cart_f,sum_scrcbs)
          end if
        end if
        if (rs_nbl(i)%nnbs.gt.0) then
          call Vforce_TAB_C(evec_thr(:,1,tpi),i,cart_f)
        end if
      end do
!     crosslinks will always be treated as short-range 
      do i=1,n_crosslinks
        call Vforce_crosslink_corr_long(i,evec_thr(:,1,tpi),sayno,cart_f,sum_scrcbs)
      end do
!
    else ! no cutoffs
!
      do i=1,nseq
!        call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
        call Vforce_rsp_long(evec_thr(:,1,tpi),i,i,sayno,cart_f,sum_scrcbs)
        if (i.lt.nseq) then
!          call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
          call Vforce_rsp_long(evec_thr(:,1,tpi),i,i+1,sayno,cart_f,sum_scrcbs)
        end if
        if (rs_nbl(i)%ntabnbs.gt.0) then
          call Vforce_TAB(evec_thr(:,1,tpi),i,cart_f)
        end if
      end do
!     note that Vforce_TAB does not skip crrosslinked interactions and that no
!     corrections are needed
    end if
!
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(6).EQV..true.)) then ! only TABUL in cart_f* so far
      hmjam%ca_f(:,:) = cart_f(:,:) + cart_f_tr(:,:)
      i = 6
      edum = evec_thr(9,1,tpi) + evec_thr(9,2,tpi)
      call els_manage_one(i,edum,cart_f,atwo)
      evec_thr(9,1,tpi) = evec_thr(9,1,tpi) + hmjam%boosts(6,1)
    end if
!
!   missing
!
  else if (ideal_run.EQV..true.) then
!
!   do nothing
!
  else
!
!   the standard calculation with support for arbitrary Hamiltonians
!   we'll do some subdistinctions for specialized routines, but keep the framework general 
!   any one-shot Hamiltonians should probably be separated out into a construct as above,
!   since distances will be calculated twice ...
!
    if (use_IMPSOLV.EQV..true.) call init_svte(3)
!
    if (use_cutoffs.EQV..true.) then
!
      if (freshnbl.EQV..true.) then
        call all_respairs_nbl(skip_frz,sayno,tpi)
      end if
!
!     first the short-range terms (ATTLJ, IPP, WCA, FOS)
      if ((is_impljp.EQV..true.).OR.(is_implj.EQV..true.)) then
!       this is the specialized set of fxns for LJ 6/12 with ABSINTH implicit solvent model
        do i=1,nseq
!         if a nb-list update is due we recompute all short-range terms within second cutoff
!         note that this should never yield any extra terms in sisa_dr
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_LJIMP_TR(evec_thr(:,2,tpi),i,cart_f_tr)
            end if
          end if
!         now compute true short-range terms
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_LJIMP_C(evec_thr(:,1,tpi),i,cart_f)
          end if
        end do
      else
!       the general case
        do i=1,nseq
!         if a nb-list update is due we recompute all short-range terms within second cutoff
!         note that this should never yield any extra terms in sisa_dr
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr)
              end if
            end if
            do j=1,rs_nbl(i)%nnbtrs
              call Vforce_rsp(evec_thr(:,2,tpi),i,rs_nbl(i)%nbtr(j),sayno,cart_f_tr)
            end do
            if (use_FEG.EQV..true.) then
              do j=1,rs_nbl(i)%ngnbtrs
                call Vforce_rsp(evec_thr(:,2,tpi),i,rs_nbl(i)%gnbtr(j),sayno,cart_f_tr)
              end do
            end if
          end if
!         now compute true short-range terms
          call Vforce_rsp(evec_thr(:,1,tpi),i,i,sayno,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1,tpi),i,i+1,sayno,cart_f)
            end if
          end if
          do j=1,rs_nbl(i)%nnbs
            call Vforce_rsp(evec_thr(:,1,tpi),i,rs_nbl(i)%nb(j),sayno,cart_f)
          end do
          if (use_FEG.EQV..true.) then
            do j=1,rs_nbl(i)%ngnbs
              call Vforce_rsp(evec_thr(:,1,tpi),i,rs_nbl(i)%gnb(j),sayno,cart_f)
            end do
          end if
        end do
      end if
!     crosslinks will always be treated as short-range 
      do i=1,n_crosslinks
        call Vforce_crosslink_corr(i,evec_thr(:,1,tpi),sayno,cart_f)
      end do
!
      if ((do_accelsim.EQV..true.).AND.(hmjam%isin(1).EQV..true.)) then ! only ATTLJ/IPP or WCA in cart_f* so far
        hmjam%ca_f(:,:) = cart_f(:,:) + cart_f_tr(:,:)
        i = 1
        edum = evec_thr(1,1,tpi) + evec_thr(1,2,tpi) + evec_thr(3,1,tpi) + evec_thr(3,2,tpi) + evec_thr(5,1,tpi) + evec_thr(5,2,tpi)
        call els_manage_one(i,edum,cart_f,atwo)
        if (use_WCA.EQV..true.) then
          evec_thr(5,1,tpi) = evec_thr(5,1,tpi) + hmjam%boosts(1,1)
        else if (use_IPP.EQV..true.) then
          evec_thr(1,1,tpi) = evec_thr(1,1,tpi) + hmjam%boosts(1,1)
        else 
          evec_thr(3,1,tpi) = evec_thr(3,1,tpi) + hmjam%boosts(1,1)
        end if
      end if
!
!     setup the screened charge terms (which have the same short-range
!     as the sisa_dr terms unless large charge groups are used with coupled models:
!     in that latter case larger range is always accounted for fully(!))
      if (use_IMPSOLV.EQV..true.) then
        call init_svte(6)
        if (use_POLAR.EQV..true.) then
          call force_setup_scrqs2()
        end if
      end if
!
      if ((do_accelsim.EQV..true.).AND.(hmjam%isin(6).EQV..true.)) then
        elsflag = .true.
        hmjam%ca_f(:,:) = cart_f(:,:)
        cart_f(:,:) = 0.0
        hmjam%ca_f_tr(:,:) = cart_f_tr(:,:)
        cart_f_tr(:,:) = 0.0
      end if
!
      if (is_implj.EQV..true.) then
!       do nothing here
      else if (is_impljp.EQV..true.) then
!       this is the specialized set of fxns for POLAR with ABSINTH implicit solvent model
        if (scrq_model.le.2) then
          do i=1,nseq
!           compute the long-range terms within short cutoff
            call Vforce_rsp_long(evec_thr(:,1,tpi),i,i,sayno,cart_f,sum_scrcbs)
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.1) then
                call Vforce_rsp_long(evec_thr(:,1,tpi),i,i+1,sayno,cart_f,sum_scrcbs)
              end if
            end if
!
            if (rs_nbl(i)%nnbs.gt.0) then
              call Vforce_PSCRM12_C(evec_thr(:,1,tpi),i,cart_f,sum_scrcbs)
            end if
!           and now add the long-range terms within twin-range
!           note that this has the funny effect that the terms in sum_scrcbs are
!           or are not incremented further on account of charge-charge interactions
!           between particles relatively far apart, i.e., it's a pseudo-cutoff
!           since the atom will still be close to one of the two charges 
            if ((freshnbl.EQV..true.)) then
              if (i.lt.nseq) then
                if (rsp_vec(i).eq.2) then
                  call Vforce_rsp_long(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr,sum_scrcbs_tr)
                end if
              end if
!
              if (rs_nbl(i)%nnbtrs.gt.0) then
                call Vforce_PSCRM12_TR(evec_thr(:,2,tpi),i,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
          end do
        else if (scrq_model.eq.4) then
          do i=1,nseq
            call Vforce_rsp_long(evec_thr(:,1,tpi),i,i,sayno,cart_f,sum_scrcbs)
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.1) then
                call Vforce_rsp_long(evec_thr(:,1,tpi),i,i+1,sayno,cart_f,sum_scrcbs)
              end if
            end if
            if (rs_nbl(i)%nnbs.gt.0) then
              call Vforce_PSCRM4_C(evec_thr(:,1,tpi),i,cart_f)
            end if
            if ((freshnbl.EQV..true.)) then
              if (i.lt.nseq) then
                if (rsp_vec(i).eq.2) then
                  call Vforce_rsp_long(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr,sum_scrcbs_tr)
                end if
              end if
              if (rs_nbl(i)%nnbtrs.gt.0) then
                call Vforce_PSCRM4_TR(evec_thr(:,2,tpi),i,cart_f_tr)
              end if
            end if
          end do
        else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          do i=1,nseq
            call Vforce_rsp_long(evec_thr(:,1,tpi),i,i,sayno,cart_f,sum_scrcbs)
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.1) then
                call Vforce_rsp_long(evec_thr(:,1,tpi),i,i+1,sayno,cart_f,sum_scrcbs)
              end if
            end if
            if (rs_nbl(i)%nnbs.gt.0) then
              call Vforce_PSCRM56_C(evec_thr(:,1,tpi),i,cart_f,sum_scrcbs)
            end if
            if ((freshnbl.EQV..true.)) then
              if (i.lt.nseq) then
                if (rsp_vec(i).eq.2) then
                  call Vforce_rsp_long(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr,sum_scrcbs_tr)
                end if
              end if
              if (rs_nbl(i)%nnbtrs.gt.0) then
                call Vforce_PSCRM56_TR(evec_thr(:,2,tpi),i,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
          end do
        else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
          do i=1,nseq
            call Vforce_rsp_long(evec_thr(:,1,tpi),i,i,sayno,cart_f,sum_scrcbs)
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.1) then
                call Vforce_rsp_long(evec_thr(:,1,tpi),i,i+1,sayno,cart_f,sum_scrcbs)
              end if
            end if
            if (rs_nbl(i)%nnbs.gt.0) then
              call Vforce_PSCRM39_C(evec_thr(:,1,tpi),i,cart_f,sum_scrcbs)
            end if
            if ((freshnbl.EQV..true.)) then
              if (i.lt.nseq) then
                if (rsp_vec(i).eq.2) then
                  call Vforce_rsp_long(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr,sum_scrcbs_tr)
                end if
              end if
              if (rs_nbl(i)%nnbtrs.gt.0) then
                call Vforce_PSCRM39_TR(evec_thr(:,2,tpi),i,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
          end do
        else if ((scrq_model.eq.7).OR.(scrq_model.eq.8)) then
          do i=1,nseq
            call Vforce_rsp_long(evec_thr(:,1,tpi),i,i,sayno,cart_f,sum_scrcbs)
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.1) then
                call Vforce_rsp_long(evec_thr(:,1,tpi),i,i+1,sayno,cart_f,sum_scrcbs)
              end if
            end if
            if (rs_nbl(i)%nnbs.gt.0) then
              call Vforce_PSCRM78_C(evec_thr(:,1,tpi),i,cart_f,sum_scrcbs)
            end if
            if ((freshnbl.EQV..true.)) then
              if (i.lt.nseq) then
                if (rsp_vec(i).eq.2) then
                  call Vforce_rsp_long(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr,sum_scrcbs_tr)
                end if
              end if
              if (rs_nbl(i)%nnbtrs.gt.0) then
                call Vforce_PSCRM78_TR(evec_thr(:,2,tpi),i,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
          end do
        else
          write(ilog,*) 'Fatal. Encountered unsupported screening model in force3(...).&
                       & This is most certainly an omission bug.'
          call fexit()
        end if
      else
!       the general case
        do i=1,nseq
!         now compute the long-range terms within short cutoff
          call Vforce_rsp_long(evec_thr(:,1,tpi),i,i,sayno,cart_f,sum_scrcbs)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp_long(evec_thr(:,1,tpi),i,i+1,sayno,cart_f,sum_scrcbs)
            end if
          end if
!
          do j=1,rs_nbl(i)%nnbs
            call Vforce_rsp_long(evec_thr(:,1,tpi),i,rs_nbl(i)%nb(j),sayno,cart_f,sum_scrcbs)
          end do
          if (use_FEG.EQV..true.) then
            do j=1,rs_nbl(i)%ngnbs
              call Vforce_rsp_long(evec_thr(:,1,tpi),i,rs_nbl(i)%gnb(j),sayno,cart_f,sum_scrcbs)
            end do
          end if
!         and now add the long-range terms within twin-range
!         note that this has the funny effect that the terms in sum_scrcbs are
!         or are not incremented further on account of charge-charge interactions
!         between particles relatively far apart, i.e., it's a pseudo-cutoff
!         since the atom will still be close to one of the two charges 
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp_long(evec_thr(:,2,tpi),i,i+1,sayno,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
!
            do j=1,rs_nbl(i)%nnbtrs
              call Vforce_rsp_long(evec_thr(:,2,tpi),i,rs_nbl(i)%nbtr(j),sayno,cart_f_tr,sum_scrcbs_tr)
            end do
            if (use_FEG.EQV..true.) then
              do j=1,rs_nbl(i)%ngnbtrs
                call Vforce_rsp_long(evec_thr(:,2,tpi),i,rs_nbl(i)%gnbtr(j),sayno,cart_f_tr,sum_scrcbs_tr)
              end do
            end if
          end if
        end do
      end if
!     finally, add long-range electrostatics corrections (note that these include
!     the next neighbor problem explicitly)
!     remember that IMPSOLV and FEG are mutually exclusive, if that ever changes
!     the construct inside the fxn below has to be changed as well!
      if (use_POLAR.EQV..true.) then
        call force_P_LR(evec_thr(:,3,tpi),cart_f_lr,sum_scrcbs_lr,freshnbl)
      end if
!     crosslinks will always be treated as short-range 
      do i=1,n_crosslinks
        call Vforce_crosslink_corr_long(i,evec_thr(:,1,tpi),sayno,cart_f,sum_scrcbs)
      end do
!
    else
!
      do i=1,nseq
        do j=i,nseq
           if (j.eq.disulf(i)) cycle
           call Vforce_rsp(evec_thr(:,1,tpi),i,j,sayno,cart_f)
        end do
      end do
!     crosslinks will always be treated as short-range 
      do i=1,n_crosslinks
        call Vforce_crosslink_corr(i,evec_thr(:,1,tpi),sayno,cart_f)
      end do
!
      if ((do_accelsim.EQV..true.).AND.(hmjam%isin(1).EQV..true.)) then ! only ATTLJ/IPP or WCA in cart_f* so far
        hmjam%ca_f(:,:) = cart_f(:,:)
        i = 1
        edum = evec_thr(1,1,tpi) + evec_thr(3,1,tpi) + evec_thr(5,1,tpi)
        call els_manage_one(i,edum,cart_f,atwo)
        if (use_WCA.EQV..true.) then
          evec_thr(5,1,tpi) = evec_thr(5,1,tpi) + hmjam%boosts(1,1)
        else if (use_IPP.EQV..true.) then
          evec_thr(1,1,tpi) = evec_thr(1,1,tpi) + hmjam%boosts(1,1)
        else 
          evec_thr(3,1,tpi) = evec_thr(3,1,tpi) + hmjam%boosts(1,1)
        end if
      end if
!
      if (use_IMPSOLV.EQV..true.) then
        call init_svte(6)
        if (use_POLAR.EQV..true.) then
          call force_setup_scrqs2()
        end if
      end if
!
      if ((do_accelsim.EQV..true.).AND.(hmjam%isin(6).EQV..true.)) then
        elsflag = .true.
        hmjam%ca_f(:,:) = cart_f(:,:)
        cart_f(:,:) = 0.0
      end if
!
      do i=1,nseq
        do j=i,nseq
          if (j.eq.disulf(i)) cycle
          call Vforce_rsp_long(evec_thr(:,1,tpi),i,j,sayno,cart_f,sum_scrcbs)
        end do
      end do
!     crosslinks will always be treated as short-range 
      do i=1,n_crosslinks
        call Vforce_crosslink_corr_long(i,evec_thr(:,1,tpi),sayno,cart_f,sum_scrcbs)
      end do
!
    end if
!
  end if
!
! the rest are cutoff-independent terms
!
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    do i=1,nseq
      call Vforce_dcbdscrq(i,cart_f)
    end do
  end if
!
  if ((do_accelsim.EQV..true.).AND.(hmjam%isin(6).EQV..true.).AND.(elsflag.EQV..true.)) then
!   at this point, we have only contributions from POLAR and TABUL in cart_f, etc.
!   this is due to the prior backup and zeroing
    i = 6
    rs = 3
    edum = evec_thr(6,1,tpi) + evec_thr(9,1,tpi) + evec_thr(6,2,tpi) + evec_thr(9,2,tpi) + evec_thr(6,3,tpi) + evec_thr(9,3,tpi)
    if (use_cutoffs.EQV..true.) then
      hmjam%ca_f_bu(:,:) = cart_f(:,:) + cart_f_tr(:,:) + cart_f_lr(:,:) ! total force for terms 6/9
    else
      hmjam%ca_f_bu(:,:) = cart_f(:,:)
    end if
    call els_manage_one(i,edum,hmjam%ca_f_bu,rs) ! this transforms hmjam%ca_f_bu to hold the entire contribution (incl. bias)
    if (use_cutoffs.EQV..true.) then
!     note that we subtract out only the unbiased _lr and _tr contributions, i.e., the bias is mapped entirely onto cart_f
      cart_f(:,:) = hmjam%ca_f(:,:) + hmjam%ca_f_bu(:,:) - cart_f_tr(:,:) - cart_f_lr(:,:)
      cart_f_tr(:,:) = hmjam%ca_f_tr(:,:) + cart_f_tr(:,:)
    else
      cart_f(:,:) = hmjam%ca_f(:,:) + hmjam%ca_f_bu(:,:)
    end if
    if (use_POLAR.EQV..true.) then
      evec_thr(6,1,tpi) = evec_thr(6,1,tpi) + hmjam%boosts(6,1)
    else 
      evec_thr(9,1,tpi) = evec_thr(9,1,tpi) + hmjam%boosts(6,1)
    end if
  end if
!
  if (use_IMPSOLV.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(4).EQV..true.)) then
      do i=1,nseq
        call force_freesolv(i,evec_thr(:,1,tpi),hmjam%ca_f)
      end do
      i = 4
      call els_manage_one(i,evec_thr(4,1,tpi),cart_f,aone)
    else
      do i=1,nseq
        call force_freesolv(i,evec_thr(:,1,tpi),cart_f)
      end do
    end if
  end if
!
  if (use_BOND(1).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(15).EQV..true.)) then
      do rs=1,nseq
        call force_bonds(rs,evec_thr(:,1,tpi),hmjam%ca_f)
      end do
      i = 15
      call els_manage_one(i,evec_thr(15,1,tpi),cart_f,aone)
    else
      do rs=1,nseq
        call force_bonds(rs,evec_thr(:,1,tpi),cart_f)
      end do
    end if
  end if
  if (use_BOND(2).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(16).EQV..true.)) then
      do rs=1,nseq
        call force_angles(rs,evec_thr(:,1,tpi),hmjam%ca_f)
      end do
      i = 16
      call els_manage_one(i,evec_thr(16,1,tpi),cart_f,aone)
    else
      do rs=1,nseq
        call force_angles(rs,evec_thr(:,1,tpi),cart_f)
      end do
    end if
  end if
  if (use_BOND(3).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(17).EQV..true.)) then
      do rs=1,nseq
        call force_impropers(rs,evec_thr(:,1,tpi),hmjam%ca_f)
      end do
      i = 17
      call els_manage_one(i,evec_thr(17,1,tpi),cart_f,aone)
    else
      do rs=1,nseq
        call force_impropers(rs,evec_thr(:,1,tpi),cart_f)
      end do
    end if
  end if
  if (use_BOND(4).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(18).EQV..true.)) then
      do rs=1,nseq
        call force_torsions(rs,evec_thr(:,1,tpi),hmjam%ca_f)
      end do
      i = 18
      call els_manage_one(i,evec_thr(18,1,tpi),cart_f,aone)
    else
      do rs=1,nseq
        call force_torsions(rs,evec_thr(:,1,tpi),cart_f)
      end do
    end if
  end if
  if (use_BOND(5).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(20).EQV..true.)) then
      do rs=1,nseq
        call force_cmap(rs,evec_thr(:,1,tpi),hmjam%ca_f)
      end do
      i = 20
      call els_manage_one(i,evec_thr(20,1,tpi),cart_f,aone)
    else
      do rs=1,nseq
        call force_cmap(rs,evec_thr(:,1,tpi),cart_f)
      end do
    end if
  end if
  if (use_CORR.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(2).EQV..true.)) then
      do rs=1,nseq
        call force_corrector(rs,evec_thr(:,1,tpi),hmjam%ca_f)
      end do
      i = 2
      call els_manage_one(i,evec_thr(2,1,tpi),cart_f,aone)
    else
      do rs=1,nseq
        call force_corrector(rs,evec_thr(:,1,tpi),cart_f)
      end do
    end if
  end if
!
  if (use_TOR.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(7).EQV..true.)) then
      do i=1,nseq
        if (par_TOR2(i).gt.0) then
          call force_torrs(evec_thr(:,1,tpi),i,hmjam%ca_f)
        end if
      end do
      i = 7
      call els_manage_one(i,evec_thr(7,1,tpi),cart_f,aone)
    else
      do i=1,nseq
        if (par_TOR2(i).gt.0) then
          call force_torrs(evec_thr(:,1,tpi),i,cart_f)
        end if
      end do
    end if
  end if
!
  if ((do_accelsim.EQV..true.).AND.(hmjam%isin(12).EQV..true.)) then
    do rs=1,nseq
      call force_boundary_rs(rs,evec_thr(:,1,tpi),1,hmjam%ca_f)
    end do
    i = 12
    call els_manage_one(i,evec_thr(12,1,tpi),cart_f,aone)
  else
    do rs=1,nseq
      call force_boundary_rs(rs,evec_thr(:,1,tpi),1,cart_f)
    end do
  end if
!
  if (use_ZSEC.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(8).EQV..true.)) then
      do imol=1,nmol
        call force_zsec_gl(imol,evec_thr(:,1,tpi),hmjam%ca_f)
      end do
      i = 8
      call els_manage_one(i,evec_thr(8,1,tpi),cart_f,aone)
    else
      do imol=1,nmol
        call force_zsec_gl(imol,evec_thr(:,1,tpi),cart_f)
      end do
    end if
  end if
!       call CPU_time(t2)
!        write(*,*) 'boundary, zsec, tor',t2-t1
!
!       call CPU_time(t1)
!
!
  if (use_DSSP.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(19).EQV..true.)) then
      call force_dssp_gl(evec_thr(:,1,tpi),hmjam%ca_f)
      i = 19
      call els_manage_one(i,evec_thr(19,1,tpi),cart_f,aone)
    else
      call force_dssp_gl(evec_thr(:,1,tpi),cart_f)
    end if
  end if
!
  if (use_EMICRO.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(21).EQV..true.)) then
      call force_emicro_gl(evec_thr(:,1,tpi),hmjam%ca_f)
      i = 21
      call els_manage_one(i,evec_thr(21,1,tpi),cart_f,aone)
    else
      call force_emicro_gl(evec_thr(:,1,tpi),cart_f)
    end if
  end if
!
  if (use_POLY.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(14).EQV..true.)) then
      do imol=1,nmol
        call force_poly_gl(imol,evec_thr(:,1,tpi),hmjam%ca_f)
      end do
      i = 14
      call els_manage_one(i,evec_thr(14,1,tpi),cart_f,aone)
    else
      do imol=1,nmol
        call force_poly_gl(imol,evec_thr(:,1,tpi),cart_f)
      end do
    end if
  end if
!
  if (use_DREST.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(10).EQV..true.)) then
      call force_drest(evec_thr(:,1,tpi),hmjam%ca_f)
      i = 10
      call els_manage_one(i,evec_thr(10,1,tpi),cart_f,aone)
    else
      call force_drest(evec_thr(:,1,tpi),cart_f)
    end if
  end if
!
  if (use_IMPSOLV.EQV..true.) then
    do i=1,n
      if (sisa_nr(i).gt.MAXNB) then
      write(ilog,*) 'WARNING. Maximum number of neighbor atoms excee&
 &ded for atom ',i,'. It is very likely this run will crash soon. Ig&
 &nore any of the results obtained.'
        call fexit()
      end if
    end do
  end if
!
!  write(*,*) invtemp*sum(evec_thr(:,1,tpi))
  if ((freshnbl.EQV..true.)) evec_tr(:) = evec_thr(:,2,1)
  if ((freshnbl.EQV..true.).OR.(lrel_md.eq.2)) evec_lr(:) = evec_thr(:,3,1)
  evec(:) = evec_thr(:,1,1) + evec_tr(:) + evec_lr(:)
  cart_f(:,:) = cart_f(:,:) + cart_f_tr(:,:) + cart_f_lr(:,:)
!
  force3 = sum(evec(:))
!
  if ((do_accelsim.EQV..true.).AND.(hmjam%isin(11).EQV..true.)) then
    hmjam%ca_f(:,:) = cart_f(:,:)
    edum = force3
    i = 11
    call els_manage_one(i,edum,cart_f,atwo)
    evec(11) = hmjam%boosts(11,1)
    force3 = force3 + evec(11)
  end if
!
  return
!
end
!
!------------------------------------------------------------------------
!
! a wrapper which calls the "correct" long-range electrostatics subroutine
! note that "long-range" really stands for multiple, fundamentally different
! things
!
subroutine force_P_LR(evec,ca_f,sum_s,freshnbl)
!
  use energies
  use atoms
  use cutoffs
  use ewalds
  use iounit
  use mcsums
  use movesets
  use sequen
!
  implicit none
!
  logical freshnbl
  RTYPE evec(MAXENERGYTERMS),ca_f(n,3)
  RTYPE sum_s(n)
!
  if (use_POLAR.EQV..false.) return
!
  if ((lrel_md.eq.4).OR.(lrel_md.eq.5)) then
    if (freshnbl.EQV..true.) then
      ca_f(:,:) = 0.0
      evec(:) = 0.0
      call force_P_LR_NBL(evec,ca_f,skip_frz,sum_s)
    end if
  else if (lrel_md.eq.2) then
!   the reciprocal space sum must be calculated at every step
    ca_f(:,:) = 0.0
    evec(:) = 0.0
    if (ewald_mode.eq.1) then
      call force_pme(evec,ca_f)
    else if (ewald_mode.eq.2) then
      call force_ewald(evec,ca_f)
    end if
  else if (lrel_md.eq.3) then
!   do nothing except increment energy by self-term (RF has no LR or quasi-LR terms)
    if (freshnbl.EQV..true.) then
      evec(:) = 0.0
      evec(6) = evec(6) + rfcnst
    end if
  else if (lrel_md.eq.1) then
!   really do nothing
  else
    write(ilog,*) 'Fatal. Called force_P_LR(...) with unsupported &
 &LR electrostatics model (offending mode is ',lrel_md,'). Please re&
 &port this bug.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
! this subroutine is a wrapper mimicking steric/SAV interactions between crosslinked
! residues as a neighbor relation -> dirty trick that works well as long
! as crosslinked residues are reliably excluded in neighbor list setup
!
subroutine Vforce_crosslink_corr(lk,evec,cut,ca_f)
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
  RTYPE evec(MAXENERGYTERMS),ca_f(n,3)
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
  call Vforce_rsp(evec,rs1,rs1+1,cut,ca_f)
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
subroutine Vforce_crosslink_corr_long(lk,evec,cut,ca_f,sum_s)
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
  integer i,rs1,rs2,bun(4),bun2,lk,ii,jj,kk,iii,kkk,mm
  RTYPE evec(MAXENERGYTERMS),ca_f(n,3),sum_s(n)
  integer bulst1(nrsnb(crosslink(lk)%rsnrs(1)),2)
  integer bulst2(nrsnb(crosslink(lk)%rsnrs(1)),2)
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
!   case for just the tabulated potential does not require any trickery (all interactions computed
!   irrespective of topology)
    call Vforce_rsp_long(evec,rs1,rs2,cut,ca_f,sum_s)
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
! for Ewald and RF also back up exclusion lists
  if ((lrel_md.eq.2).OR.(lrel_md.eq.3)) then
    bun2 = nrexpolnb(rs1)
    bulst2(1:nrexpolnb(rs1),1) = iaa(rs1)%expolnb(1:nrexpolnb(rs1),1)
    bulst2(1:nrexpolnb(rs1),2) = iaa(rs1)%expolnb(1:nrexpolnb(rs1),2)
    nrexpolnb(rs1) = 0
  end if
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
            if ((lrel_md.eq.2).OR.(lrel_md.eq.3)) then
              nrexpolnb(rs1) = nrexpolnb(rs1) + 1
              iaa(rs1)%expolnb(nrexpolnb(rs1),1) = ii
              iaa(rs1)%expolnb(nrexpolnb(rs1),2) = kk
            end if
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
  if ((lrel_md.eq.2).OR.(lrel_md.eq.3)) then
    call force_rsp_long_mod(evec,rs1,rs1+1,cut,ca_f)
  else
    call Vforce_rsp_long(evec,rs1,rs1+1,cut,ca_f,sum_s)
  end if
! and restore what we destroyed
  if ((lrel_md.eq.2).OR.(lrel_md.eq.3)) then
    nrexpolnb(rs1) = bun2
    iaa(rs1)%expolnb(1:nrexpolnb(rs1),1) = bulst2(1:nrexpolnb(rs1),1)
    iaa(rs1)%expolnb(1:nrexpolnb(rs1),2) = bulst2(1:nrexpolnb(rs1),2)
  end if
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
    call Vforce_rsp_long(evec,rs1,rs2,cut,ca_f,sum_s)
    use_POLAR = bul2
  end if
!
end
!
!-------------------------------------------------------------------------------
!
! wrapper wrapper
!
subroutine all_crosslink_corrs(evec,cut,ca_f,sum_s)
!
  use sequen
  use energies
  use atoms
!
  implicit none
!
  logical cut
  RTYPE evec(MAXENERGYTERMS),ca_f(n,3),sum_s(n)
  integer lk
!
  do lk=1,n_crosslinks
    call Vforce_crosslink_corr(lk,evec,cut,ca_f)
    call Vforce_crosslink_corr_long(lk,evec,cut,ca_f,sum_s)
  end do
!
end
!
!------------------------------------------------------------------------------------
!
! this is a subroutine providing support for foreign energy calculations
! in REMD runs
! care has to be taken such that forces from the original Hamiltonian
! are retained in principle (this routine does NOT deal with swaps per se)
!
subroutine lamforce(rve,fve)
!
  use energies
  use atoms
  use iounit
  use molecule
  use mpistuff
  use system
  use forces
  use units
  use cutoffs
  use dssps
  use ems
!
  implicit none
!
  RTYPE rve(mpi_nodes),fve(mpi_nodes)
#ifdef ENABLE_MPI
  integer i,j,jj,k,imol,which,i_start,i_end
  RTYPE force3,evec(MAXENERGYTERMS)
  RTYPE vbu(MAXREDIMS),eva,dum,dum2
  RTYPE evec_tr(MAXENERGYTERMS),evec_lr(MAXENERGYTERMS)
  logical needforcebu,needsavup,needemup,badflg(MAXREDIMS),atrue
  integer vbui(MAXREDIMS)
  RTYPE, ALLOCATABLE:: cfbu(:,:),ifbu(:)
#endif
!
#ifdef ENABLE_MPI
!
  needforcebu = .false.
  needsavup = .false.
  needemup = .false.
  atrue = .true.
  rve(:) = 0.0
  fve(:) = 0.0
  dum2 = 0.0
!
  do i=1,re_conddim
!    if ((re_types(i).eq.7)).OR.&
! &      ((re_types(i).ge.12).AND.(re_types(i).le.19))) then
!      write(ilog,*) 'Fatal. REMD with the chosen condition (',&
! &re_types(i),') is not yet supported. Check back later.'
!      call fexit() 
!    end if
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
!   if 32 do nothing
    if (re_types(i).eq.33) vbu(i)=scale_EMICRO
    if (re_types(i).eq.34) vbu(i)=emthreshdensity
  end do
!
  if (needforcebu.EQV..true.) then
    allocate(cfbu(n,3))
    if (fycxyz.eq.2) then
      cfbu(1:n,:) = cart_f(1:n,:)
    else if (fycxyz.eq.1) then
      cfbu(1:n,:) = cart_f(1:n,:)
      k = 0
      do imol=1,nmol
        k = k + size(dc_di(imol)%f)
      end do
      allocate(ifbu(k))
      jj = 0
      do imol=1,nmol
        do j=1,size(dc_di(imol)%f)
          jj = jj + 1
          ifbu(jj) = dc_di(imol)%f(j)
        end do
      end do
    end if
  end if
!
  if (reol_all.EQV..true.) then
    i_start = 1
    i_end = re_conditions
  else
    i_start = max((myrank+1)-1,1)
    i_end = min((myrank+1)+1,re_conditions)
  end if
  do i=i_start,i_end
!  do i=1,re_conditions
!   set the other conditions (as requested: either neighbors only or all)
!   remember that we don't change the use_XX-flags
!   this implies, however, that for scale_XX = 0.0, we compute
!   a bunch of terms all multiplied by 0.0. in order to preserve
!   this information (to compute derivatives with respect to scale_XX),
!   we therefore use a little detour (set to 1.0, subtract out)
    eva = 0.0
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
        if (lrel_md.eq.3) then
          call setup_rfcnst2()
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
      else if (re_types(j).eq.21) then
        scale_FEGS(3) = re_mat(i,j)
        call setup_parFEG(3)
      else if (re_types(j).eq.22) then
        scale_FEGS(6) = re_mat(i,j)
        call setup_parFEG(6)
        if (lrel_md.eq.3) then
          call setup_rfcnst()
        end if
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
        par_DSSP(9) = re_mat(i,j)
      else if (re_types(j).eq.31) then
        par_DSSP(7) = re_mat(i,j)
!     if 32 do nothing
      else if (re_types(j).eq.33) then
        scale_EMICRO = re_mat(i,j)
        if (scale_EMICRO.le.0.0) then
          badflg(j) = .true.
          scale_EMICRO = 1.0
        end if
      else if (re_types(j).eq.34) then
        emthreshdensity = re_mat(i,j)
        needemup = .true.
      end if
    end do
    if (needsavup.EQV..true.) call absinth_savprm()
    if (needemup.EQV..true.) then
      call scale_emmap(emthreshdensity,dum,dum2)
      call precompute_diff_emmap()
    end if
!   now compute the current energy at this condition with
!   our structure
!   note that we correct for energy/structure mismatch for our home node
!   (artifact of the leapfrog algorithm) by re-computing ourselves as well(!)
    eva = invtemp*force3(evec,evec_tr,evec_lr,atrue)
!   the computation of derivatives is really only meaningful when
!   noTI is false
    fve(i) = 0.0
    do j=1,re_conddim
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
    rve(i) = sum(evec(:))
    rve(i) = rve(i)*invtemp
    needsavup = .false.
    needemup = .false.
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
      if (lrel_md.eq.3) then
        call setup_rfcnst2()
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
      if (lrel_md.eq.3) then
        call setup_rfcnst()
      end if
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
!   if 32 do nothing
    else if (re_types(i).eq.33) then
      scale_EMICRO = vbu(i)
    else if (re_types(i).eq.34) then
      emthreshdensity = vbu(i)
      needemup = .true.
    end if
  end do
  if (needsavup.EQV..true.) call absinth_savprm()
  if (needemup.EQV..true.) then
    call scale_emmap(emthreshdensity,dum,dum2)
    call precompute_diff_emmap()
  end if
!
  if (needforcebu.EQV..true.) then
    if (fycxyz.eq.2) then
      cart_f(1:n,:) = cfbu(1:n,:)
    else if (fycxyz.eq.1) then
      cart_f(1:n,:) = cfbu(1:n,:)
      jj = 0
      do imol=1,nmol
        do j=1,size(dc_di(imol)%f)
          jj = jj + 1
          dc_di(imol)%f(j) = ifbu(jj)
        end do
      end do
      deallocate(ifbu)
    end if
    deallocate(cfbu)
  end if
#else
!
  write(ilog,*) 'Fatal. Called lamforce(...) in non MPI-calculation.&
 & This is most definitely a bug.'
  call fexit()
#endif
!
end
!
!---------------------------------------------------------------------cc
!
! this routine assumes that there is only a single ghosted
! residue with effective scaled charges of par_FEG2(9)
! we also assume this residue is a dipole residue, not a group
! with a net charge
!
subroutine setup_rfcnst()
!
  use energies
  use iounit
  use cutoffs
  use sequen
  use polypep
  use atoms
  use system
  use math
  use units
!
  implicit none
!
  integer rs,i,k
  RTYPE resQ,netQ,istr,t1,t2,kap,epsr
!
  if (use_FEG.EQV..false.) return
!
  istr= 0.0
!
! we can only assume these are homogeneously distributed (see comments in
! grf_setup() in polar.f)
  do i=1,cglst%ncs
    if (par_FEG(atmres(cglst%it(i))).EQV..true.) then
      write(ilog,*) 'Fatal. FEG in conjunction with (G)RF treatment &
 & is only allowed for net-neutral residues.'
      call fexit()
    end if
    istr = istr + cglst%tc(i)*cglst%tc(i)
  end do
  epsr = par_IMPSOLV(2)
  netQ = 0.0
  k = 0
  do rs=1,nseq
    resQ = 0.0
    if (par_FEG3(rs).EQV..true.) then
!     remember that fully de-coupled residues are treated as if they had zero charge
      cycle
    else if (par_FEG(rs).EQV..true.) then
      k = k + 1
      do i=1,at(rs)%npol
        resQ = resQ + atq(at(rs)%pol(i))*atq(at(rs)%pol(i))
      end do
      netQ = netQ + par_FEG2(9)*par_FEG2(9)*resQ
    else
      do i=1,at(rs)%npol
        resQ = resQ + atq(at(rs)%pol(i))*atq(at(rs)%pol(i))
      end do
      netQ = netQ + scale_POLAR*resQ
    end if
  end do
!
! standard RF (GRF limiting case with zero ionic strength)
  if (rf_mode.eq.2) then
    par_POLAR(1) = (epsr-1.0)/((2.0*epsr+1.0)*(mcel_cutoff**3.0))
    par_POLAR(2) = 1.0/mcel_cutoff + mcel_cutoff2*par_POLAR(1)
! GRF
  else
    kap = sqrt(4.0*PI*electric*invtemp*istr/ens%insV)
    t1 = 1.0 + kap*mcel_cutoff
    t2 = kap*kap*mcel_cutoff2
    par_POLAR(1) = (t1*(epsr-1.0) + 0.5*epsr*t2)/&
 &        (((2.0*epsr+1.0)*t1 + epsr*t2)*(mcel_cutoff**3.0))
    par_POLAR(2) = 1.0/mcel_cutoff + mcel_cutoff2*par_POLAR(1)
  end if
  rfcnst = -0.5*electric*netQ*par_POLAR(2)
!
end
!
!-----------------------------------------------------------------------
!
subroutine setup_rfcnst2()
!
  use energies
  use iounit
  use cutoffs
  use sequen
  use polypep
  use atoms
  use system
  use math
  use units
!
  implicit none
!
  integer rs,i
  RTYPE netQ,istr,t1,t2,kap,epsr
!
  if (use_FEG.EQV..true.) then
    call setup_rfcnst()
    return
  end if
! 
  istr= 0.0
!
! we can only assume these are homogeneously distributed (see comments in
! grf_setup() in polar.f)
  do i=1,cglst%ncs
    istr = istr + cglst%tc(i)*cglst%tc(i)
  end do
  epsr = par_IMPSOLV(2)
  netQ = 0.0
  do rs=1,nseq
    do i=1,at(rs)%npol
      netQ = netQ + atq(at(rs)%pol(i))*atq(at(rs)%pol(i))
    end do
  end do
!
! standard RF (GRF limiting case with zero ionic strength)
  if (rf_mode.eq.2) then
    par_POLAR(1) = (epsr-1.0)/((2.0*epsr+1.0)*(mcel_cutoff**3.0))
    par_POLAR(2) = 1.0/mcel_cutoff + mcel_cutoff2*par_POLAR(1)
! GRF
  else
    kap = sqrt(4.0*PI*electric*invtemp*istr/ens%insV)
    t1 = 1.0 + kap*mcel_cutoff
    t2 = kap*kap*mcel_cutoff2
    par_POLAR(1) = (t1*(epsr-1.0) + 0.5*epsr*t2)/&
 &        (((2.0*epsr+1.0)*t1 + epsr*t2)*(mcel_cutoff**3.0))
    par_POLAR(2) = 1.0/mcel_cutoff + mcel_cutoff2*par_POLAR(1)
  end if
  rfcnst = -0.5*electric*scale_POLAR*netQ*par_POLAR(2)
!
end
!
!-----------------------------------------------------------------------
