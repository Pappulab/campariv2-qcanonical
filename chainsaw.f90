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
!
#include "macros.i"
!
!
! ############################################
! ##                                        ##
! ## THIS IS THE CAMPARI MASTER SOURCE FILE ##
! ##                                        ##
! ############################################
!
!
program chainsaw
!
  use iounit
  use torsn
  use mcsums
  use energies
  use atoms
  use sequen
  use mcgrid
  use pdb
  use molecule
  use polyavg
  use movesets
  use paircorr
  use mpistuff
  use cutoffs
  use system
  use contacts
  use dipolavg
  use forces
  use diffrac
  use mini
  use dssps
  use clusters
  use shakeetal
  use grandensembles
  use ems
  use martin_own
  use polypep
  use inter
  use fyoc !martin added
  use params
  use aminos
  use accept
  use zmatrix
  
  implicit none
!
#ifdef ENABLE_MPI
#include "mpif.h"
#endif
  character*1000 WW ! martin add
  character*1000 WW2 ! martin add
  character*1000 WW3 ! martin add
  integer urr
  
  integer dt(8)
  character (len=12) rc(3)
  integer rs,aone,imol
  integer istep,ndump,m1,ifos
  integer firststep,i,j,k
  RTYPE ee,ee1,ee2,random,t1,t2,eit
  RTYPE energy,energy3,edum(MAXENERGYTERMS),force1,force3
  RTYPE, ALLOCATABLE:: blas(:)
  RTYPE ed1(MAXENERGYTERMS),ed2(MAXENERGYTERMS)
  RTYPE tt
  integer alcsz,temp_var_int1(3),temp_int! Martin added
  integer, allocatable:: memory(:) ! Martin added 
  integer, allocatable:: temp_var3(:)
  integer, allocatable:: temp_var5(:,:),temp_var6(:,:),temp_var7(:),dpgrp_check(:,:,:)
  integer, allocatable:: take(:)
  integer, allocatable:: polin_added1(:,:),polin_added2(:,:),polnb_added1(:,:),polnb_added2(:,:)!Martin added
  integer limit,nsalt
  RTYPE temp_sum,temp_1,temp_2
  
  integer sub,at1,at2,dp1,dp2,ii,kk,rs1,rs2,maxlim
  character WWW*200000,WW_tmp*1000,my_format*10
  character*1 tab
  
!
#ifdef ENABLE_THREADS
  integer omp_get_num_threads
#endif
  character(MAXSTRLEN) intfile
  logical logdummy,atrue,afalse,logd2,logbuff,good
#ifdef ENABLE_MPI
  integer masterrank,modstep,ierr
#endif
!


DEBUG_ENER_DRIFT=0.0

#ifdef ENABLE_MPI
  masterrank = 0
#endif
!
!  allocate(lj_att_memory(200,200))
!  allocate(lj_rep_memory(200,200))

  allocate(blas(3))
  atrue = .true.
  afalse = .false.
  aone = 1
  m1 = -1
#ifdef ENABLE_MPI
!
! MPI: obviously some things need to be different ...
! this routine will have the master distribute the simulation information
! and run the setup procedures for all CPUs
  time_comm = 0.0
  call CPU_time(ee1)
  call MPI_StartMC()
  call CPU_time(t2)
  time_comm = time_comm + t2 - ee1
#else
!
! when we start
  call Date_and_Time(rc(1),rc(2),rc(3),dt)
  write(*,*)
  call write_license()
  write(*,*) 'Execution started at ',rc(2)(1:2),':',rc(2)(3:4),&
 &' on ',rc(1)(5:6),'/',rc(1)(7:8),'/',rc(1)(1:4),'.'
  call CPU_time(ee1)
  write(*,*)
!     
! all initializations for start up
  call initial()
  logdummy = .false.
!
! we need to read (not yet parse) the key (input) file
  call getkey()
!
! level 1 input reader
  call parsekey(1)
! level 2 input reader
  call parsekey(2)
! level 3 input reader
  call parsekey(3)
#endif
!
  if ((pdb_analyze.eqv..false.).and.(print_det.eqv..true.)) then 
      print *,"Solvation energy details analysis can only be done in analysis runs. Turning off."
      print_det=.false.
  end if 
  
  
  logbuff = use_trajidx
  use_trajidx = .false.
!
! set the derived box variables
  call update_bound(atrue)
! build the system 
  call makepept()

  ! Even if this new residue is not in the sequence
! resolve problems with requests for structural input
  call strucinp_sanitychecks()

  if (use_cutoffs.EQV..false.) then ! MARTIN : Overall empty statement, unecessary
!   overwrite checkfreq
    if (dyn_mode.eq.1) then
!      nsancheck = nsim + 1
    end if
  end if 
  
  if (do_hsq.eqv..false.) then 
      hsq_freq=0.0
      enout=min(enout,xyzout)
      xyzout=min(enout,xyzout)
  else if (hsq_freq.eq.0.0) then 
      do_hsq=.false.
  end if 
  
  if ((do_pka.eqv..true.).or.(do_hsq.eqv..true.)) then 
     good=.true.
     do i=1,n_biotyp
      if ((bio_ljtyp(i).lt.31).and.(bio_ljtyp(i)-bio_ljtyp(transform_table(i)).ne.0)) then 
          print *,"Wrong reference to LJ type for biotype",i,"(transforms to ",transform_table(i),", LJ from ",bio_ljtyp(i)," to "&
          &,bio_ljtyp(transform_table(i)),")"
          good=.false.
      end if
    end do 
    if (good.eqv..false.) then 
      print *,"One or more titrable biotype that has changing LJ type does not refer to a titrable LJ type. This will break ALL&
      & residues using that biotype. Fix parameters file before continuing."
      call fexit()
    end if 
    if (par_top.eq.0) then 
        print *,"You need to be using the USE_PARAM_TOP options with the DO_PKA or DO_HSQ options. Turning on."
        par_top=1
    end if 
    if (use_BOND(1).EQV..true.) then 
        print *,"You should not be using the bond length potential with the DO_PKA or DO_HSQ options. Turning off."
        use_BOND(1)=.false.
    end if      
  end if 
  
  if (((do_pka.eqv..true.).or.(do_hsq.eqv..true.)).and.(use_POLAR.EQV..true.)) then
    !First, checking a few thing relative to the pKa calculations
    
      ! Martin : this is a failsafe in case the newly added biotypes do not work, which would results in all residue potentially 
    !breaking
 
    ! Martin : still have to align the arrays : 
    !   -polin
    !   -polnb
    !   -expolin
    !   -expolnb
    ! I think that's done now
    !!!! Test for the new subroutine to handle all of this
    

    allocate(temp_which_dpg(n,3))
    allocate(temp_nrpolnb(nseq,3))
    allocate(temp_nrpolintra(nseq,3))
    allocate(temp_nrexpolin(nseq,3))
    allocate(temp_nrexpolnb(nseq,3))
    allocate(ats_temp(nseq,3))
    
    temp_which_dpg(:,:)=0
    temp_nrpolnb(:,:)=0
    temp_nrpolintra(:,:)=0
    temp_nrexpolin(:,:)=0
    temp_nrexpolnb(:,:)=0

    if ((nhis.ne.0).and.(do_hsq.eqv..true.)) then
        call qcann_polar_group(3)
    end if
    

    
    call qcann_polar_group(2)

    call qcann_polar_group(1)
    
    ! Martin : keep the biggest size as the effective size, to make sure we  iterate through all elements

    if (nhis.ne.0) then
        do i=1,n
            if (atq_limits(i,1).eq.atq_limits(i,4)) then 
                atq_limits(i,2)=atq_limits(i,4)
            else 
                atq_limits(i,2)=0.0
            end if 
        end do 
    end if 

    do i=1,n
        if (atq_limits(i,1).eq.atq_limits(i,3)) then 
            atq_limits(i,2)=atq_limits(i,3)
        else 
            atq_limits(i,2)=0.0
        end if 
    end do 

    do limit=1,3
        do rs=1,nseq 
            
            if (temp_nrpolnb(rs,limit).gt.nrpolnb(rs)) then 
                nrpolnb(rs)=temp_nrpolnb(rs,limit)
            end if 

            if (temp_nrpolintra(rs,limit).gt.nrpolintra(rs)) then 
                nrpolintra(rs)=temp_nrpolintra(rs,limit)
            end if 

            if (temp_nrexpolin(rs,limit).gt.nrexpolin(rs)) then 
                nrexpolin(rs)=temp_nrexpolin(rs,limit)
            end if 

            if (temp_nrexpolnb(rs,limit).gt.nrexpolnb(rs)) then 
                nrexpolnb(rs)=temp_nrexpolnb(rs,limit)
            end if 
        end do 
    end do 

    do rs=1,nseq
        do i=1,at(rs)%ndpgrps
            allocate(at(rs)%dpgrp(i)%qnm_limits(4))
            allocate(at(rs)%dpgrp(i)%tc_limits(4))         
!            ! Since the net charge is not interpolable (integer)), I will give the maximum value
!            ! This is fine since it seems it is only used for cycling through loops
        end do
    end do
      allocate(cglst%tc_limits(cglst%ncs,3))
    if (nhis.ne.0) then 

        call setup_charge_last_step(1,3)
    end if 
    call setup_charge_last_step(1,2)
        
    deallocate(temp_which_dpg)
    deallocate(temp_nrpolnb)
    deallocate(temp_nrpolintra)
    deallocate(temp_nrexpolin)
    deallocate(temp_nrexpolnb)
    deallocate(ats_temp)
  else if (use_POLAR.EQV..true.) then 
    call polar_groups()
  end if
  
  if (use_IMPSOLV.EQV..true.) then
    call setup_freesolv()
    call solvation_groups()
      if ((do_pka.eqv..true.).or.(do_hsq.eqv..true.)) then
        do imol=1,nmol ! Martin : this bit is taken from unbound's firtst subroutine
            do rs=rsmol(imol,1),rsmol(imol,2)
                do ifos=1,at(rs)%nfosgrps ! Goes though all the FOS of the groups composing the residue
                    if (allocated(at(rs)%fosgrp(ifos)%wts_limits).eqv..false.) then 
                        allocate(at(rs)%fosgrp(ifos)%wts_limits(at(rs)%fosgrp(ifos)%nats,2))
                    end if 
                    if (at(rs)%fosgrp(ifos)%val_limits(1,1).eq.at(rs)%fosgrp(ifos)%val_limits(1,2)) then 
                        at(rs)%fosgrp(ifos)%wts_limits(:,1)=at(rs)%fosgrp(ifos)%wts(:)
                        at(rs)%fosgrp(ifos)%wts_limits(:,2)=at(rs)%fosgrp(ifos)%wts(:)
                        at(rs)%fosgrp(ifos)%val_limits(:,1)=at(rs)%fosgrp(ifos)%val(:)
                        at(rs)%fosgrp(ifos)%val_limits(:,2)=at(rs)%fosgrp(ifos)%val(:)
                    end if 
                    if (nhis.ne.0) then 
                        if (ANY(his_state(:,1)==rs)) then
                            allocate(temp_var5(at(rs)%fosgrp(ifos)%nats,2))
                            temp_var5(:,:)=at(rs)%fosgrp(ifos)%wts_limits(:,:)
                            deallocate(at(rs)%fosgrp(ifos)%wts_limits)
                            allocate(at(rs)%fosgrp(ifos)%wts_limits(at(rs)%fosgrp(ifos)%nats,3))
                            at(rs)%fosgrp(ifos)%wts_limits(:,1:2)=temp_var5(:,:)
                            at(rs)%fosgrp(ifos)%wts_limits(:,1:2)=temp_var5(:,:)
                            at(rs)%fosgrp(ifos)%val_limits(:,3)=at(rs)%fosgrp(ifos)%val_limits(:,2)
                            deallocate(temp_var5)
                        end if
                    end if
                end do
            end do
        end do
      end if
  end if

  if (use_TABUL.EQV..true.) then
    call read_tabfiles()
  end if
!
  if (use_POLY.EQV..true.) then
    call read_polfile()
  end if
!
  if (use_DREST.EQV..true.) then
    call read_drestfile()
  end if

  if ((do_pka.eqv..true.).or.(do_hsq.eqv..true.)) then
    ! Martin : these loops are to correct the array so that the order of the bonds are the same in both limits
    allocate(temp_nrsdi(nseq,3))
    allocate(temp_nrsba(nseq,3))
    allocate(temp_nrsbl(nseq,3))
    allocate(temp_nrsimpt(nseq,3))
    if ((nhis.ne.0).and.(do_hsq.eqv..true.)) then 
        call qcann_assign_bndtprms(3)
    end if 
    call qcann_assign_bndtprms(2)
    call qcann_assign_bndtprms(1)

    iaa(:)=iaa_limits(:,1)
 
    if (nhis.ne.0) then 
        maxlim=3
    else 
        maxlim=2
    end if 
    do limit=1,maxlim
        do rs=1,nseq
            if (nrsdieff(rs).lt.temp_nrsdi(rs,limit)) then 
                nrsdieff(rs)=temp_nrsdi(rs,limit)
            end if 
            if (nrsbleff(rs).lt.temp_nrsbl(rs,limit)) then 
                nrsbleff(rs)=temp_nrsbl(rs,limit)
            end if 
            if (nrsbaeff(rs).lt.temp_nrsba(rs,limit)) then 
                nrsbaeff(rs)=temp_nrsba(rs,limit)
            end if 
            if (nrsimpteff(rs).lt.temp_nrsimpt(rs,limit)) then 
                nrsimpteff(rs)=temp_nrsimpt(rs,limit)
            end if 
        end do    
    end do 
    
    ! Put that back to have the corrected atoms
    iaa(:)=iaa_limits(:,2)
    
    call qcann_iaa_index_corrections(1,2)    

    deallocate(temp_nrsdi)
    deallocate(temp_nrsba)
    deallocate(temp_nrsbl)
    deallocate(temp_nrsimpt)

    iaa(:)=iaa_limits(:,2)
    fudge(:)=fudge_limits(:,2)
  else 
    call assign_bndtprms()
  end if 

  if (do_pka.eqv..true.) then 
    call scale_stuff_pka()
  end if 

  if (do_hsq.eqv..true.) then 

      if (hsq_mode.eq.2) then 
        call initialize_charge_state_cste_FOS()
      else if (hsq_mode.eq.1) then 
        call initialize_charge_state_cste_net_Q()
      end if 
    
    call scale_stuff_hsq()

    call switch_bonds()

  end if

  if (phfreq.gt.0.0) then
     call ionize_setup()
  end if
!
  if (do_restart.EQV..true.) then
!   do nothing (see below)
  else if (pdbinput.EQV..true.) then
    if (pdb_readmode.eq.1) then
      call FMSMC_readpdb()
    else
      call FMSMC_readpdb3()
    end if
  else if (fycinput.EQV..true.) then
    call readfyc()
  end if

 call randomize_bb() 
! torsional setup requires all structure manipulation to be complete
  if (use_TOR.EQV..true.) then
    call read_torfile()
  end if
!
! FEG requires all other energy setup to be complete
  if (use_FEG.EQV..true.) then
    call read_fegfile()
  end if
!
! Particle fluctuation setup also requires setup to be complete
  if ((ens%flag.ge.5).AND.(ens%flag.le.6)) then
    call read_particleflucfile()
  end if
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call remc_sanitychecks()
  end if
#endif
!  debug_pka=.true.
  if (do_hsq.eqv..true.) then
      if (debug_pka.eqv..true.) then
        call print_limits_summary()
      end if 
  end if 

  if (do_hs.eqv..true.) then
      write(ilog,*) 'Using the HS using mode ',HS_mode,' and freq ',hs_freq
      
      
      if (HS_mode.eq.1) then 
            call read_HS_residues()
            
      else if (HS_mode.eq.2) then
          allocate(hs_residues(nseq))
          do i=1,nseq
            hs_residues(i)=i
          end do 
      else if (HS_mode.eq.3) then
          nsalt=0
          do i=1,nseq
              if ((amino(seqtyp(i)).eq.'NA+').or.(amino(seqtyp(i)).eq.'CL-')) then
                  nsalt=nsalt+1
              end if 
          end do 
          
          allocate(hs_residues(nsalt))
          nsalt=0
          do i=1,nseq
              if ((amino(seqtyp(i)).eq.'NA+').or.(amino(seqtyp(i)).eq.'CL-')) then 
                  nsalt=nsalt+1
                  hs_residues(nsalt)=i
              end if 
          end do 
        end if 
    ! This is to simplify and use the same convention as HSQ
    ! Note that while all do have to be initialized, most will contain redundant information
    ! The only things that are to be modified in the second limit are 
    !       -the fudge factors (not the charge or the LJ param)
    !       -the FOS weights
    allocate(fudge_limits(nseq,2))
    allocate(iaa_limits(nseq,2))
    allocate(par_hsq(nseq,4))
    
    do rs=1,nseq
        do i=1,at(rs)%ndpgrps
            allocate(at(rs)%dpgrp(i)%qnm_limits(3))
            allocate(at(rs)%dpgrp(i)%tc_limits(3))
            
            at(rs)%dpgrp(i)%qnm_limits(:)=at(rs)%dpgrp(i)%qnm
            at(rs)%dpgrp(i)%tc_limits(:)=at(rs)%dpgrp(i)%tc
        end do 
        
        allocate(fudge_limits(rs,1)%elin(size(fudge(rs)%elin)))
        allocate(fudge_limits(rs,2)%elin(size(fudge(rs)%elin)))
        
        allocate(fudge_limits(rs,1)%elnb(size(fudge(rs)%elnb)))
        allocate(fudge_limits(rs,2)%elnb(size(fudge(rs)%elnb)))

        if (nrsintra(rs).ne.0) then  
            allocate(fudge_limits(rs,1)%rsin(size(fudge(rs)%rsin(:))))
            allocate(fudge_limits(rs,2)%rsin(size(fudge(rs)%rsin(:))))
            
            allocate(fudge_limits(rs,1)%rsin_ljs(size(fudge(rs)%rsin_ljs)))
            allocate(fudge_limits(rs,2)%rsin_ljs(size(fudge(rs)%rsin_ljs)))

            allocate(fudge_limits(rs,1)%rsin_lje(size(fudge(rs)%rsin_lje)))
            allocate(fudge_limits(rs,2)%rsin_lje(size(fudge(rs)%rsin_lje)))

            fudge_limits(rs,1)%rsin_ljs(:)=fudge(rs)%rsin_ljs(:)
            fudge_limits(rs,1)%rsin_lje(:)=fudge(rs)%rsin_lje(:)
            
            fudge_limits(rs,2)%rsin_ljs(:)=fudge(rs)%rsin_ljs(:)
            fudge_limits(rs,2)%rsin_lje(:)=fudge(rs)%rsin_lje(:)
            
            fudge_limits(rs,1)%rsin(:)=fudge(rs)%rsin(:)
            fudge_limits(rs,2)%rsin(:)=1.0 ! Has to be 1. The scaling is going to be done by changing the scale, in order to control 
                                            ! LJ_rep and LJ_att independently 
        end if 

        if (nrsnb(rs).ne.0) then  
            allocate(fudge_limits(rs,1)%rsnb(size(fudge(rs)%rsnb(:))))
            allocate(fudge_limits(rs,2)%rsnb(size(fudge(rs)%rsnb(:))))
            
            allocate(fudge_limits(rs,1)%rsnb_ljs(size(fudge(rs)%rsnb_ljs)))
            allocate(fudge_limits(rs,2)%rsnb_ljs(size(fudge(rs)%rsnb_ljs)))
        
            allocate(fudge_limits(rs,1)%rsnb_lje(size(fudge(rs)%rsnb_lje)))
            allocate(fudge_limits(rs,2)%rsnb_lje(size(fudge(rs)%rsnb_lje)))
            
            fudge_limits(rs,1)%rsnb_ljs(:)=fudge(rs)%rsnb_ljs(:)
            fudge_limits(rs,1)%rsnb_lje(:)=fudge(rs)%rsnb_lje(:)

            fudge_limits(rs,2)%rsnb_ljs(:)=fudge(rs)%rsnb_ljs(:)
            fudge_limits(rs,2)%rsnb_lje(:)=fudge(rs)%rsnb_lje(:)
            
            fudge_limits(rs,1)%rsnb(:)=1.0
            fudge_limits(rs,2)%rsnb(:)=1.0! Has to be 1. The scaling is going to be done by changing the scale, in order to control 
                                            ! LJ_rep and LJ_att independently 
        end if 
        
        if (nrpolintra(rs).gt.0) then 
            fudge_limits(rs,1)%elin(:)=fudge(rs)%elin(:)
            fudge_limits(rs,2)%elin(:)=0.
        end if
        if (nrpolnb(rs).gt.0) then 
            fudge_limits(rs,1)%elnb(:)=fudge(rs)%elnb(:)
            fudge_limits(rs,2)%elnb(:)=0.
        end if 

        allocate(iaa_limits(rs,1)%par_bl(nrsbl(rs),MAXBOPAR+1))
        allocate(iaa_limits(rs,2)%par_bl(nrsbl(rs),MAXBOPAR+1))
        
        allocate(iaa_limits(rs,1)%par_ba(nrsba(rs),MAXBAPAR+1))
        allocate(iaa_limits(rs,2)%par_ba(nrsba(rs),MAXBAPAR+1))
        
        allocate(iaa_limits(rs,1)%par_di(nrsdi(rs),MAXDIPAR))
        allocate(iaa_limits(rs,2)%par_di(nrsdi(rs),MAXDIPAR))
        
        allocate(iaa_limits(rs,1)%par_impt(size(iaa(rs)%par_impt(:,1)),MAXDIPAR))
        allocate(iaa_limits(rs,2)%par_impt(size(iaa(rs)%par_impt(:,1)),MAXDIPAR))
        
        iaa_limits(rs,1)%par_bl(:,:)=iaa(rs)%par_bl(:,:)
        iaa_limits(rs,2)%par_bl(:,:)=iaa(rs)%par_bl(:,:)
        
        iaa_limits(rs,1)%par_ba(:,:)=iaa(rs)%par_ba(:,:)
        iaa_limits(rs,2)%par_ba(:,:)=iaa(rs)%par_ba(:,:)
        
        iaa_limits(rs,1)%par_di(:,:)=iaa(rs)%par_di(:,:)
        iaa_limits(rs,2)%par_di(:,:)=iaa(rs)%par_di(:,:)
        
        iaa_limits(rs,1)%par_impt(:,:)=iaa(rs)%par_impt(:,:)
        iaa_limits(rs,2)%par_impt(:,:)=iaa(rs)%par_impt(:,:)
        
        do ifos=1,at(rs)%nfosgrps
            allocate(at(rs)%fosgrp(ifos)%wts_limits(at(rs)%fosgrp(ifos)%nats,2))
            
            at(rs)%fosgrp(ifos)%wts_limits(:,1)=at(rs)%fosgrp(ifos)%wts(:)
            at(rs)%fosgrp(ifos)%wts_limits(:,2)=0.
            
            at(rs)%fosgrp(ifos)%val_limits(:,1)=at(rs)%fosgrp(ifos)%val(:)
            at(rs)%fosgrp(ifos)%val_limits(:,2)=at(rs)%fosgrp(ifos)%val(:)
        end do 
    end do 
    allocate(atq_limits(size(atq),3))
    atq_limits(:,1)=atq(:)
    atq_limits(:,2)=atq(:)
    atq_limits(:,3)=atq(:)
    
    allocate(cglst%tc_limits(size(cglst%tc),3))
    cglst%tc_limits(:,1)=cglst%tc(:)
    cglst%tc_limits(:,2)=0.
    cglst%tc_limits(:,3)=0.
    
    allocate(atr_limits(size(atr),2))
    atr_limits(:,1)=atr(:)
    atr_limits(:,2)=atr(:)
    
    allocate(atbvol_limits(size(atbvol),2))
    atbvol_limits(:,1)=atbvol(:)
    atbvol_limits(:,2)=atbvol(:)
    
    allocate(atvol_limits(size(atvol),2))
    atvol_limits(:,1)=atvol(:)
    atvol_limits(:,2)=atvol(:)
    
    allocate(atsavred_limits(size(atsavred),2))
    atsavred_limits(:,1)=atsavred(:)
    atsavred_limits(:,2)=atsavred(:)
    
    allocate(atsavmaxfr_limits(size(atsavmaxfr),2))
    atsavmaxfr_limits(:,1)=atsavmaxfr(:)
    atsavmaxfr_limits(:,2)=atsavmaxfr(:)
    
    allocate(atstv_limits(size(atstv),2))
    atstv_limits(:,1)=atstv(:)
    atstv_limits(:,2)=atstv(:)
    
    allocate(atsavprm_limits(n,size(atsavprm(1,:)),2))
    atsavprm_limits(:,:,1)=atsavprm(:,:)
    atsavprm_limits(:,:,2)=atsavprm(:,:)
    
    
    allocate(blen_limits(n,2))
    blen_limits(:,1)=blen(:)
    blen_limits(:,2)=blen(:)
    
    allocate(bang_limits(n,2))
    bang_limits(:,1)=bang(:)
    bang_limits(:,2)=bang(:)
  end if 


  if (do_pka.eqv..true.) then 
      deallocate(x_limits)
      deallocate(y_limits)
      deallocate(z_limits)
      
  end if 

! if equilibration run manually turn off analysis, so we don't do
! unnecessary and potentially ill-defined setup work
  if (nequil.ge.nsim) then! Martin shouldn't it be turn off analysis if nsim==0 and nequil !=0 ?cause that is a real 
                            !manual equilibration
    pccalc = nsim+1
    covcalc = nsim+1
    segcalc = nsim+1
    torlccalc = nsim+1
    polcalc = nsim+1
    rhcalc = nsim+1
    sctcalc = nsim+1
    holescalc = nsim+1
    savcalc = nsim+1
    contactcalc = nsim+1
    particlenumcalc = nsim+1
    clucalc = nsim+1
    angcalc = nsim+1
    dipcalc = nsim+1
    intcalc = nsim+1
    dsspcalc = nsim+1
    cstorecalc = nsim+1
    xyzout = nsim+1
    phout = nsim+1
    diffrcalc = nsim+1
    emcalc = nsim+1
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then!MARTIN : Not sur this is necessary since this is for manual equilibration runs, it probably shouldnt
!need any RE stuff, but I need to look into what re_olcalc is
      re_olcalc = nsim+1
    end if
#endif
    if (pdb_analyze.EQV..true.) then
      write(ilog,*) 'Warning. No analysis performed in PDB analysis &
 &mode due to equilibration period exceeding snapshot number.'
      write(ilog,*)
    end if
  end if! MARTIN end of the manual equilibration part
!
  if (pdb_analyze.EQV..true.) then
    if (use_dyn.EQV..true.) then
      write(ilog,*) 'Warning. Reconstruction of dynamical variables in PDB an&
 &alysis mode not yet supported. Nonetheless using force routines to compute energies (if &
 &applicable).'
      write(ilog,*)
    end if
    if (ens%flag.ne.1) then
      write(ilog,*) 'Fatal. Ensembles other than NVT are currently n&
 &ot supported by PDB analysis mode.'
      write(ilog,*)
      call fexit()
    end if
  end if
!
  if (just_solutes.EQV..true.) then
    if (nsolutes.le.0) then
      write(ilog,*) 'Warning. Suppressing structural output for a system &
 &composed exclusively of solvent is not supported. Turning off suppression&
 & flag.'
      just_solutes = .false.
      pdbeffn = n
    end if
  end if
!
! temporary shutdown of LCT stuff
  if (use_LCTOR.EQV..true.) then
    write(ilog,*) 'Warning. LCT potential term is currently not supp&
 &orted. Turning off.'
    use_LCTOR = .false.
    scale_LCTOR = 0.0
  end if
  if (use_lctmoves.EQV..true.) then
    write(ilog,*) 'Warning. LCT moves are currently not supported. T&
 &urning off.'
    use_lctmoves = .false.
  end if
  if (torlccalc.le.nsim) then
    write(ilog,*) 'Warning. LCT analysis is currently not supported.&
 & Turning off.'
    torlccalc = nsim + 1
  end if
!
  
  
  
! pdb analysis mode requires more setup work
  if (pdb_analyze.EQV..true.) then
!    call randomize_bb()
    if (align%yes.EQV..true.) then
      call read_alignfile(align)
      if (pdb_fileformat.le.2) then
        if (use_pdb_template.EQV..true.) then
          call read_pdb_template()
        end if
      end if
    end if
    if (pdb_fileformat.eq.1) call setup_pdbtraj()
    if (pdb_fileformat.eq.2) call setup_pdbsnaps()
#ifdef LINK_XDR
    if (pdb_fileformat.eq.3) call setup_xtctraj()
#endif
    if (pdb_fileformat.eq.4) call setup_dcdtraj()
#ifdef LINK_NETCDF
    if (pdb_fileformat.eq.5) call setup_netcdftraj()
#endif
#ifdef ENABLE_MPI
    if (select_frames.EQV..true.) then
      write(ilog,*) 'Warning. Trajectory subset analysis of user-selected framesfile is currently not &
 &supported in MPI trajectory analysis mode. Disabled.'
      select_frames = .false.
    end if
#else
    if (select_frames.EQV..true.) call read_framesfile()
#endif
  end if
!
  if (do_restart.EQV..false.) then
!   dump out a pdb- and Z matrix-file for sanity check
    if ((pdb_analyze.EQV..false.).OR.(n_pdbunk.gt.0)) then
      intfile = 'START '
      call FMCSC_dump(intfile)
    end if
  end if
!
! read in file containing information about regions in
! phi/psi space to sample (center of boxes, boxsizes)
  if (use_stericgrids.EQV..true.) then
    call allocate_stericgrid(aone)
    call readgrid()
  end if
!
! allocate torsional analysis arrays (up here because of ZSEC-
! coupling)
  call allocate_torsn(aone)
!
! setup for LC torsional analysis
  if ((covcalc.le.nsim).AND.(covmode.eq.3)) then
    call setup_lc_tor()
  end if
!
! setup work for DSSP analysis (here because of DSSP-restraint coupling)
  if ((dsspcalc.le.nsim).OR.(use_DSSP.EQV..true.)) then
    call setup_dssp()
  end if
!
#ifdef LINK_NETCDF
!
! read in file-based EM restraint map
  if (use_EMICRO.EQV..true.) then
    call read_netcdf3Dmap()
  end if
!
#endif
!
!
! now we proceed through some sanity checks for the completed settings
!
#ifdef ENABLE_MPI
!
! Sanity check against pdb-analyzer (not supported in MPI)
!
  if ((pdb_analyze.EQV..true.).AND.(use_REMC.EQV..false.)) then
    write(ilog,*) 'Fatal. Requested trajectory analysis mode in an MPI setting without using &
 &replica exchange setup. This is not supported (please enable FMCSC_REMC).'
    call fexit()
  end if
#endif
!
! setup for domain grids
  if ((do_restart.EQV..false.).AND.(use_mcgrid.EQV..true.)) then
    call allocate_mcgrid(aone)
    call setupmcgrid()
    if (grid%report.EQV..true.) then
      call gridreport()
    end if
  end if
!
! setup work for rigid-body moves
  rblst%nr = nmol
  allocate(rblst%idx(rblst%nr))
  do imol=1,nmol
    rblst%idx(imol) = imol
  end do
!
! Hamiltonian sanity checks
  call hamiltonian_sanitychecks()
!
! a sanity check against a meaningless calculation
  if ((fycxyz.eq.1).AND.(ntorpuck.eq.0).AND.(nmol.eq.1)) then
    write(ilog,*) 'Warning. Attempting a calculation with no relevan&
 &t degrees of freedom (1 rigid molecule in a box). Input error?'
  end if
!
! big sanity check routine for move set
  if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    call moveset_sanitychecks()
  end if
!
! setup work for torsional dynamics (must precede read_frzfile())
  if (((dyn_mode.ne.1).AND.(fycxyz.eq.1)).OR.((cstorecalc.le.nsim).AND.&
 &((cdis_crit.eq.2).OR.(cdis_crit.eq.4)))) then
    call set_IMDalign()
  end if
!
! handling of constraints for internal coordinate spaces and adjustment of MC sampling weights
  if (do_frz.EQV..true.) then
    call read_frzfile()
  end if
  call preferential_sampling_setup()
!
  call dynamics_sanitychecks()
  call loop_checks()
!
! setup work for constraints in Cartesian dynamics
  if ((dyn_mode.ne.1).AND.(fycxyz.eq.2)) then
    call shake_setup()
  end if
!
! T-coupling groups for thermostatting (if nothing else, at least for analysis)
  if (dyn_mode.ne.1) then
    call read_tgrpfile()
    call init_thermostat()
  end if
!
! allocation for movesets
  if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    call allocate_moveset(aone)
  end if
!
! setup work for Ewald sums (might need loop_checks to check on Hamiltonian)
! and (G)RF electrostatics
  if (dyn_mode.ne.1) then
    if (lrel_md.eq.2) then
      call setup_ewald()
    else if (lrel_md.eq.3) then
      call grf_setup()
      if (use_FEG.EQV..true.) then
        call setup_rfcnst()
      end if
    end if
  end if
!
  if (dyn_mode.ne.1) then 
    if (use_TABUL.EQV..true.) then
      do rs=1,nseq
        call tab_respairs_nbl_pre(rs)
      end do
    end if
  end if
!
! setup work for WL sampling
  if (mc_acc_crit.eq.3) call wl_init()
!
! setup work for energy landscape sculpting
  call els_init()
! 
! monitor initial setup time
  call CPU_time(t2)
  eit = t2 - ee1
  
  
  
  !Martin : allocates the enfresolv memory array befroe the first energy computation
  !allocate (en_fos_memory(nseq,4))
  !allocate (sav_memory(nseq,4,20))
  

  !Martin : creates a pointer to access nseq from a subroutine
  !stock=nseq
  !nseq_pointer => stock
  
  
  ! monitor time for initial energy/force calculation
  if ((do_restart.EQV..false.).AND.(pdb_analyze.EQV..false.)) then
!
    call CPU_time(t1)
!
!    evtl.y perform a gradient check
    if ((use_dyn.EQV..true.).AND.(grad_check.EQV..true.)) then
!     for now, just use these hard-coded parameters
      blas(1) = 0.0001
      blas(2) = 0.0001
      blas(3) = 0.0001
      call gradient_test(blas(1),blas(2),blas(3))
    end if
!
!   use force routines unless in straight Monte Carlo !!Martin : This is the first energy assesement
    if (dyn_mode.eq.1) then
      if (do_n2loop.EQV..true.) then
          
          esave=energy(esterms)
          
      end if
      if (use_cutoffs.EQV..true.) then
        esavec = energy3(edum,atrue)! martin : label 1
        if (do_n2loop.EQV..false.) then
          write(ilog,*) 'Warning. Using cutoff-based initial energy also as "reference" energy (due to FMCSC_N2LOOP).'
          esave = esavec
          
          esterms(:) = edum(:)
        end if
      end if
    else 
      if (do_n2loop.EQV..true.) then
        esave = force1(esterms)
      end if
      if (use_cutoffs.EQV..true.) then
        esavec = force3(edum,ed1,ed2,atrue)
        if (do_n2loop.EQV..false.) then
          write(ilog,*) 'Warning. Using cutoff-based initial energy also as "reference" energy (due to FMCSC_N2LOOP).'
          esave = esavec
          esterms(:) = edum(:)
        end if
      end if
    end if
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
    call CPU_time(t1) 
!
!   setup work for rigid body moves (not sure whether or why this would be needed)
    if (rigidfreq.gt.0.0) then
      do imol=1,nmol
        call update_rigid(imol)
      end do
    end if
!
!   cycle initialization for hybrid methods (this is overwritten in restart cases)
    if ((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
      in_dyncyc = .false.
      curcyc_end = first_mccyclen
      curcyc_start = 1
      curcyc_end = min(curcyc_end,nsim)
#ifdef ENABLE_MPI
      if (use_REMC.EQV..true.) then
        if (mod(curcyc_end,re_freq).eq.0) then
          curcyc_end = curcyc_end + 1
        end if
      end if
#endif
    end if
  else
    call CPU_time(t1)
  end if


! setup for pc analysis
  if (pccalc.lt.nsim) then
    call read_gpcfile()
  end if

! sanity checks and setup for structural clustering
  if (cstorecalc.le.nsim) then
    call read_clusteringfile()
  end if

! if still relevant, read breaks file
  if (cstorecalc.le.nsim) then
    call read_trajbrksfile()
  end if

! setup work for internal coordinate histograms
  if (intcalc.le.nsim) then
    call setup_inthists()
  end if

! additional SAV requests
  if (savcalc.le.nsim) then
    call read_savreqfile()
  end if

! various sanity checks for requested analyses
  call analysis_checks()

! set up file I/O (get handles/open files)
  call makeio(aone)

! setup for indexed trajectory output
  use_trajidx = logbuff
  if (use_trajidx.EQV..true.) call read_trajidxfile()

! setup for covariance analysis (needs I/O to be done)
  if (covcalc.le.nsim) then
    call setup_covar_int()
  end if

! setup for pH calculations (mostly informational output) (needs I/O to be done)
  if (phfreq.gt.0.0) then
     call ph_setup()
  end if

! setup for diffraction calculations
  if (diffrcalc.le.nsim) then
    call bessel_read()
  end if

! setup for torsional output (needs I/O to be done)
  if (torout.le.nsim) then
    call torsion_header()
  end if

! before we proceed to the actual calculation we take care of
! all memory allocation which hasn't happened yet
  call allocate_rest()
! with memory allocated, get volume element for PC
  if (pccalc.le.nsim)  call get_voli_pc()
!
! print summary before proceeding
  call summary()
!     
! set local counter to zero
  ndump = 0
!
  tt = 0.0

! switch back for restart purposes
  if (do_restart.EQV..true.) then
    call CPU_time(t2)
    eit = eit + t2 - t1
    call CPU_time(t1)
!   grid can only be initialized properly after reading restart
!   so use residue-based cut in first force/energy calculation
    if (use_cutoffs.EQV..true.) then
      logdummy = use_mcgrid
      logd2 = use_rescrit
      use_mcgrid = .false.
      use_rescrit = .true.   
    end if
    call read_restart()
    if (use_cutoffs.EQV..true.) then
      use_mcgrid = logdummy
      use_rescrit = logd2
    end if
    if (use_mcgrid.EQV..true.) then
      call allocate_mcgrid(aone)
      call setupmcgrid()
      if (grid%report.EQV..true.) then
        call gridreport()
      end if
    end if
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
    call CPU_time(t1)
  end if
  
  if (pccalc.le.nsim) call get_voli_pc()
  firststep = nstep + 1
!
! deal with hybrid method cycle counters
  if ((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    if (nstep.lt.nsim) then
      if (in_dyncyc.EQV..false.) then
        mvcnt%nmcseg = mvcnt%nmcseg + 1
        mvcnt%avgmcseglen = mvcnt%avgmcseglen + 1.0*(curcyc_end-max(nstep+1,curcyc_start)+1)
      else
        mvcnt%ndynseg = mvcnt%ndynseg + 1
        mvcnt%avgdynseglen = mvcnt%avgdynseglen + 1.0*(curcyc_end-max(nstep+1,curcyc_start)+1)
      end if
    end if
  end if
!
! finalize time counter for initial setup
!
  call CPU_time(t2)
  eit = eit + t2 - t1

  curframe = 0


  if ((do_pka.eqv..true.).or.(do_hsq.eqv..true.)) then 
      if (debug_pka.eqv..true.) then 

        call write_pka_details(1)
        !call print_limits_summary()
        
      end if 
    call print_values_summary()
  end if 
!  call print_values_summary()
!  call print_limits_summary()
!  flush(6)
!  flush(0)
  done_dets=.false.
    
  do istep = firststep,nsim ! MARTIN : here is the main loop      
    if (print_det.eqv..true.) then 
        if (mod(istep,sav_det_freq).eq.0) then 
            if (sav_det_mode.eq.1) then 
                call write_detailed_energy_1(en_fos_memory,istep)
            else 
                call write_detailed_energy_2(sav_memory,istep)
            end if 
        end if 
    end if 
!   in the pdb-analysis mode a "move" is the next structure
    if (pdb_analyze.EQV..true.) then ! strictly for analysis using a pre determined trajectory
      nstep = nstep + 1
      if (pdb_fileformat.eq.1) call FMSMC_readpdb2(istep)
      if (pdb_fileformat.eq.2) call FMSMC_readpdb2(istep)
#ifdef LINK_XDR
      if (pdb_fileformat.eq.3) call FMSMC_readxtc(istep)! MARTIN : this is the common pdb+xtc mode
#endif
      if (pdb_fileformat.eq.4) call FMSMC_readdcd(istep)
#ifdef LINK_NETCDF
      if (pdb_fileformat.eq.5) call FMSMC_readnetcdf(istep)
#endif
      if (select_frames.EQV..true.) then
        if (istep.eq.framelst(curframe+1)) then
          curframe = curframe + 1
        else
          cycle
        end if
      end if
#ifdef ENABLE_MPI
      if (re_aux(3).eq.1) call MPI_RETrajMode(istep) ! un-REX trajectories
#endif
      call CPU_time(t1)
      if (use_dyn.EQV..true.) then
        esave = force3(esterms,esterms_tr,esterms_lr,atrue)
      else
        esave = energy3(esterms,atrue)
      end if
      call CPU_time(t2)
      time_energy = time_energy + t2 - t1
      if (align%yes.EQV..true.) then
        if (mod(istep,align%calc).eq.0) call struct_align()
      end if
      call CPU_time(t1)
      call mcstat(istep,ndump)
      call CPU_time(t2)
      time_analysis = time_analysis + t2 - t1
      if ((mod(istep,nsancheck).eq.0).AND.((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
        call CPU_time(t1)
        call energy_check()
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
      end if

      !call write_detailed_energy(en_fos_memory,istep)
      cycle
    end if!MARTIN end of the loop for the analysis run
!
    if (use_dyn.EQV..true.) then ! martin ; this is only relevent for md
!
      if (dyn_mode.eq.2) then
#ifdef ENABLE_MPI
        if (use_REMC.EQV..true.) then
          modstep = mod(istep,re_freq)
          if ((modstep.eq.0).and.(HS_in_progress.EQV..false.)) then
            call MPI_REMaster(istep,ndump)
            cycle
          end if
        end if
#endif
        if (fycxyz.eq.1) then
          call int_mdmove(istep,ndump,afalse)
        else
          call cart_mdmove(istep,ndump,afalse)
        end if
!
      else if (dyn_mode.eq.3) then
#ifdef ENABLE_MPI
        if (use_REMC.EQV..true.) then
          modstep = mod(istep,re_freq)
          if ((modstep.eq.0).and.(HS_in_progress.EQV..false.)) then
            call MPI_REMaster(istep,ndump)
            cycle
          end if
        end if
#endif
        if (fycxyz.eq.1) then
          call int_ldmove(istep,ndump,afalse)
        else
          call cart_ldmove(istep,ndump,afalse)
        end if
!
      else if (dyn_mode.eq.6) then
!       this is a bit ugly, as the minimizers internalize the loop
!       hence, exit immediately when finished 
        if ((mini_mode.gt.0).AND.(mini_mode.lt.4)) then
          call minimize(nsim,mini_econv,mini_stepsize,mini_mode)
          exit
        else if (mini_mode.eq.4) then
          call stochastic_min(nsim)
          exit
        end if 
!
!     hybrid method requires special adjustment
      else if ((dyn_mode.eq.5).OR.(dyn_mode.eq.7)) then
!       if new cycle, switch
!       note that RE-moves should never be placed on curcyc_start or curcyc_end
        if ((istep.eq.curcyc_start).AND.(istep.ne.firststep)) then
          if (in_dyncyc.EQV..true.) then
            in_dyncyc = .false.
          else
            in_dyncyc = .true.
          end if
        end if
#ifdef ENABLE_MPI
        if (use_REMC.EQV..true.) then
          modstep = mod(istep,re_freq)
          if ((modstep.eq.0).and.(HS_in_progress.EQV..false.)) then
            call MPI_REMaster(istep,ndump)
            if (in_dyncyc.EQV..false.) then
              if (mod(istep,nsancheck).eq.0) then
                call CPU_time(t1)
                call energy_check()
                call CPU_time(t2)
                time_energy = time_energy + t2 - t1
              end if
            end if
            cycle
          end if
        end if
#endif
        if (in_dyncyc.EQV..true.) then
          logdummy = .false.
          if (istep.eq.curcyc_start) then
            logdummy = .true.
            if ((fycxyz.eq.2).AND.(cart_cons_mode.gt.1).AND.(istep.gt.(first_mccyclen+1))) logdummy = .false.
          end if
          if (dyn_mode.eq.5) then
            if (fycxyz.eq.1) then
              call int_mdmove(istep,ndump,logdummy)
            else
              call cart_mdmove(istep,ndump,logdummy)
            end if
          else if (dyn_mode.eq.7) then
            if (fycxyz.eq.1) then
              call int_ldmove(istep,ndump,logdummy)
            else
              call cart_ldmove(istep,ndump,logdummy)
            end if
          end if
          if (istep.eq.curcyc_end) then
            curcyc_start = istep + 1
            curcyc_end = istep + min_mccyclen + floor(random()*(max_mccyclen-min_mccyclen) + 0.5)
            curcyc_end = min(curcyc_end,nsim)
#ifdef ENABLE_MPI
            if (use_REMC.EQV..true.) then
!             important: re_freq is never less than 2
              if (mod(curcyc_start,re_freq).eq.0) then
                curcyc_start = curcyc_start + 1
              end if
              if (mod(curcyc_end,re_freq).eq.0) then
                curcyc_end = curcyc_end + 1
              end if
!             publish the new cycle to all nodes and receive on slaves (note that this will override
!             the cycle settings just obtained on the slave node)
              call MPI_SyncHybridCycle(istep)
            else if (use_MPIAVG.EQV..true.) then
              call MPI_SyncHybridCycle(istep)
            end if

#endif
            if (istep.lt.nsim) then
              mvcnt%nmcseg = mvcnt%nmcseg + 1
              mvcnt%avgmcseglen = mvcnt%avgmcseglen + 1.0*(curcyc_end-curcyc_start+1)
            end if
          end if
        else
!         re-set energy to MC
          if (istep.eq.curcyc_start) then
            logdummy = .true.
!           this corrects for image mismatches due to comm / com - differences in PBC
            call update_rigid_mc_all()
            esave = energy3(esterms,logdummy)
          end if
!         pick an elementary, non-RE move type we want
          call select_mcmove_tree()
!         displace, evaluate and deal with it
          call mcmove(istep,ndump)
!
          if (mod(istep,nsancheck).eq.0) then
            call CPU_time(t1)
            call energy_check()
            call CPU_time(t2)
            time_energy = time_energy + t2 - t1
          end if
          if (istep.eq.curcyc_end) then
            curcyc_start = istep + 1
            curcyc_end = istep + min_dyncyclen + floor(random()*(max_dyncyclen-min_dyncyclen) + 0.5)
            curcyc_end = min(curcyc_end,nsim)
#ifdef ENABLE_MPI
            if (use_REMC.EQV..true.) then
!             important: re_freq is never less than 2
              if (mod(curcyc_start,re_freq).eq.0) then
                curcyc_start = curcyc_start + 1
              end if
              if (mod(curcyc_end,re_freq).eq.0) then
                curcyc_end = curcyc_end + 1
              end if
!             publish the new cycle to all nodes and receive on slaves (note that this will override
!             the cycle settings just obtained on the slave node)
              call MPI_SyncHybridCycle(istep)
            else if (use_MPIAVG.EQV..true.) then
              call MPI_SyncHybridCycle(istep)
            end if
#endif
            if (istep.lt.nsim) then
              mvcnt%ndynseg = mvcnt%ndynseg + 1
              mvcnt%avgdynseglen = mvcnt%avgdynseglen + 1.0*(curcyc_end-curcyc_start+1)
            end if
          end if
        end if
!
      end if
!
      call CPU_time(t1)

      call mcstat(istep,ndump)
      call CPU_time(t2)
      time_analysis = time_analysis + t2 - t1
!
#ifdef ENABLE_MPI
      if ((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.).AND.(istep.gt.1)) then
        modstep = mod(istep,re_freq)
        if ((modstep.eq.0).AND.(istep.lt.nsim)) then
          call MPI_ASMaster(istep) ! does not increment step -> do not cycle
        end if
      end if
#endif
!
    else
!
!   determine what movetype we want
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then ! martin : in case of the use replica exchange, and actual run 
      modstep = mod(istep,re_freq)
      if ((modstep.eq.0).and.(HS_in_progress.EQV..false.)) then ! Martin : check if it is time to try and exahcnge replicas
        call MPI_REMaster(istep,ndump)
        if (mod(istep,nsancheck).eq.0) then ! CHekc if it is time to do a sanity check 
          call CPU_time(t1)
          call energy_check()
          call CPU_time(t2)
          time_energy = time_energy + t2 - t1
        end if
        cycle ! Martin : Cycles through the main loop, because if a RE move is accepted, no other type of move should be tested this step
      end if ! Martin :: although I could argue that since the move is not necessarly acccepted, we should probably still the iteration if it is not 
    end if
#endif
!   pick an elementary, non-RE move type we want
    call select_mcmove_tree()
    do while (hsq_move.eqv..true.) ! if the move selected was an HSQ move, execute the HSQ bit, and try
        call do_HSQ_walk(istep)!another move until you have some other type of move
        call select_mcmove_tree()
    end do 
    do while (hs_move.eqv..true.) ! if the move selected was an HSQ move, execute the HSQ bit, and try
        call do_HS_walk(istep)!another move until you have some other type of move
        call select_mcmove_tree()
    end do 
!  write(*,*) pivot_move,fyc_move,omega_move,chi_move,ph_move,pucker_move
!  write(*,*) particleinsertion_move,particledeletion_move,particleidentity_move
!  write(*,*) UJconrot_move,SJconrot_move,torcr_move,torcr_move_omega
!  write(*,*) nuc_move,nucpuck_move,nuccr_move
!  write(*,*) rigid_move,clurb_move,trans_move,rot_move,hsq_move
!
    call mcmove(istep,ndump)
!
!   write out current statistics
    
    call CPU_time(t1)
    call mcstat(istep,ndump)
    call CPU_time(t2)
    time_analysis = time_analysis + t2 - t1

#ifdef ENABLE_MPI
    if ((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.).AND.(istep.gt.1)) then
      modstep = mod(istep,re_freq)
      if ((modstep.eq.0).AND.(istep.lt.nsim)) then
        call MPI_ASMaster(istep) ! does not increment step -> do not cycle
      end if
    end if
#endif
!
    if ((mod(istep,nsancheck).eq.0)) then!@.and.(do_hsq.eqv..false.).and.(do_pka.eqv..false.).and.(do_hs.eqv..false.)) then ! I am getting rid of this because of 
    ! lack of support for the cutoff version 
    ! That should no longer be the case, correcting
    
      call CPU_time(t1)

      call energy_check()
      call CPU_time(t2)
      time_energy = time_energy + t2 - t1
    end if
!
    end if
  end do ! martin : end of the main loop : the simulation/ analysis is over when this point is reached
!
  intfile = 'END '
  call FMCSC_dump(intfile)
!
! if in pdb-analysis mode, we might have to close the trajectory file
  if (pdb_analyze.EQV..true.) then
    if (pdb_fileformat.eq.1) call close_pdbtraj()
#ifdef LINK_XDR
    if (pdb_fileformat.eq.3) call close_xtctraj()
#endif
    if (pdb_fileformat.eq.4) call close_dcdtraj()
#ifdef LINK_NETCDF
    if (pdb_fileformat.eq.5) call close_netcdftraj()
#endif
  end if
!
! finalize I/O (close files)
  call makeio(2)
  
!
! summarize performance information
  call CPU_time(ee2)
  ee = ee2 - ee1
  write(ilog,*)
  write(ilog,*) 'Total CPU time elapsed [s]     : ',ee
  write(ilog,*) 'Fraction for energy functions  : ',100.0*time_energy/ee,'%'
  write(ilog,*) 'Fraction for analysis routines : ',100.0*time_analysis/ee,'%'
  if (time_struc.gt.0.0) write(ilog,*) 'Fraction for chain closure     : ',100.0*time_struc/ee,'%'
  if (time_holo.gt.0.0) write(ilog,*) 'Fraction for constraint solvers: ',100.0*time_holo/ee,'%'
  if (time_ph.gt.0.0) write(ilog,*) 'Fraction for pH routines       : ',100.0*time_ph/ee,'%'
  if (time_pme.gt.0.0) write(ilog,*) 'Fraction for PME reciprocal sum: ',100.0*time_pme/ee,'%'
  if (time_nbl.gt.0.0) write(ilog,*) 'Fraction for neighbor lists    : ',100.0*time_nbl/ee,'%'
  write(ilog,*) 'Fraction for initial setup     : ',100.0*eit/ee,'%'
#ifdef ENABLE_MPI
  write(ilog,*) 'Fraction for communication     : ',100.0*time_comm/ee,'%'
  write(ilog,*) 'Fraction of remainder          : ',100.0*(1.0 - time_energy/ee - time_analysis/ee - time_struc/ee -&
 &time_comm/ee - time_ph/ee - time_holo/ee - time_nbl/ee - time_pme/ee - eit/ee),'%'
#else
  write(ilog,*) 'Fraction of remainder          : ',100.0*(1.0 - time_energy/ee - time_analysis/ee - time_struc/ee -&
 & time_ph/ee - time_holo/ee - time_nbl/ee - time_pme/ee - eit/ee),'%'
#endif
 if (do_hs.eqv..true.) then 
     write(ilog,*) 'Number of tries for     HS : ',    mvcnt%nhs
     write(ilog,*) 'Acceptance  for     HS : ',    acc%nhs
     write(ilog,*) 'Acceptance ratio for     HS : ',    acc%nhs/mvcnt%nhs

     
     
 end if 
 
 
 if (do_hsq.eqv..true.) then
    write(ilog,*) 'Number of HSQ moves tried',mvcnt%nhsq
    write(ilog,*) 'Number of HSQ moves accepted',acc%nhsq
    if (mvcnt%nhsq.ne.0) then 
        temp_1=acc%nhsq
        temp_2=mvcnt%nhsq
        temp_sum=temp_1/temp_2
        write(ilog,*) 'Fraction of HSQ moves accepted',temp_sum
    end if 
    write(ilog,*) 'The folowing is the charge states transition matrix number of tries, with the line and column being the &
    &charge state index as indicated in the input charge states file. The value is the number of tries from line index to column&
    & index. The sum over the comlumn for each states should be significant, meaning they have been proposed a high enough number&
    &of times'
    
    do i=1,n_charge_states
        WWW=''
        do j=1,n_charge_states
            WW_tmp=''
            write(WW_tmp,'(I0,A1)') acc%HSQ_state_tries(i,j),char(9)
            WWW=trim(WWW)//trim(WW_tmp)
        end do 
        write(my_format,'(A2,I0,A1)')'(A',len(trim(WWW)),')'
        write(HSQ_tries,my_format) WWW
    end do
    
    write(ilog,*) 'The following is the charge states transition matrix number of acceptances.'
    do i=1,n_charge_states
        WWW=''
        
        do j=1,n_charge_states
            WW_tmp=''
            write(WW_tmp,'(I0,A1)') acc%HSQ_state_acc(i,j),char(9)
            WWW=trim(WWW)//trim(WW_tmp)
        end do 
        write(my_format,'(A2,I0,A1)')'(A',len(trim(WWW)),')'
        write(HSQ_acc,my_format) WWW
    end do
    

    write(ilog,*) 'The following is the charge states transition matrix acceptance ratio.'
    do i=1,n_charge_states
        WWW=''
        do j=1,n_charge_states
            WW_tmp=''
            if (acc%HSQ_state_tries(i, j).ne.0) then 
                temp_1=acc%HSQ_state_acc(i,j)
                temp_2=acc%HSQ_state_tries(i,j)
                temp_sum=temp_1/temp_2
                write(WW_tmp,'(F6.5,A1)') temp_sum,char(9)
            else 
                write(WW_tmp,'(F6.5,A1)') 0.0,char(9)
            end if 
            WWW=trim(WWW)//trim(WW_tmp)
        end do 
        write(my_format,'(A2,I0,A1)')'(A',len(trim(WWW)),')'
        write(HSQ_ratio,my_format) WWW
    end do
    
    WWW=''
    write(ilog,*) 'For a simpler visualization : this is over the column sum of tries and acceptance )'
        do j=1,n_charge_states
            WW_tmp=''
            if (acc%HSQ_state_tries(i, j).ne.0) then 
                temp_1=SUM(acc%HSQ_state_acc(:,j))
                temp_2=SUM(acc%HSQ_state_tries(:,j))
                temp_sum=temp_1/temp_2
                write(WW_tmp,'(I0,I0,F6.5,A1)') SUM(acc%HSQ_state_acc(:,j)),SUM(acc%HSQ_state_tries(:,j)),temp_sum,char(9)
            else 
                write(WW_tmp,'(I0,I0,F6.5,A1)') SUM(acc%HSQ_state_acc(:,j)),SUM(acc%HSQ_state_tries(:,j)),temp_sum,char(9)
            end if 
            WWW=trim(WWW)//trim(WW_tmp)
        end do 
        write(my_format,'(A2,I0,A1)')'(A',len(trim(WWW)),')'
        write(HSQ_summary,my_format) WWW
 end if
 
  write(ilog,*)
! and clean-up of course ....
  call deallocate_all()
  
  deallocate(blas)

  call Date_and_Time(rc(1),rc(2),rc(3),dt)
  write(ilog,*) 'Execution ended at ',rc(2)(1:2),':',rc(2)(3:4),&
 &' on ',rc(1)(5:6),'/',rc(1)(7:8),'/',rc(1)(1:4),'.'
!
#ifdef ENABLE_MPI
  call makelogio(2)
#endif
!
! in parallel mode, need to shutdown MPI universe
#ifdef ENABLE_MPI
  call MPI_StopMC()
  call MPI_FINALIZE(ierr)
#endif
!   
end
!
!------------------------------------------------------------------------------
!
subroutine write_license()
!
  write(*,*) '--------------------------------------------------------------------------'
  write(*,*)
  write(*,*) '                  ---        Welcome to CAMPARI       ---'
  write(*,*)
  write(*,*) '                  ---       Released Version 2.0      ---'
  write(*,*)
  write(*,*) '    Copyright (C) 2014, CAMPARI Development Team'
  write(*,*)
  write(*,*) '    Websites: http://campari.sourceforge.net'
  write(*,*) '              http://sourceforge.net/projects/campari/'
  write(*,*) '              http://pappulab.wustl.edu/campari'
  write(*,*)
  write(*,*) '    CAMPARI is free software: you can redistribute it and/or modify'
  write(*,*) '    it under the terms of the GNU General Public License as published by'
  write(*,*) '    the Free Software Foundation, either version 3 of the License, or'
  write(*,*) '    (at your option) any later version.'
  write(*,*)
  write(*,*) '    CAMPARI is distributed in the hope that it will be useful,'
  write(*,*) '    but WITHOUT ANY WARRANTY; without even the implied warranty of'
  write(*,*) '    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
  write(*,*) '    GNU General Public License for more details. '
  write(*,*)
  write(*,*) '    You should have received a copy of the GNU General Public License' 
  write(*,*) '    along with CAMPARI. If not, see <http://www.gnu.org/licenses/>.'
  write(*,*)
  write(*,*)
  write(*,*) '--------------------------------------------------------------------------'
  write(*,*)
  write(*,*)
!
end
!
!------------------------------------------------------------------------------
!

