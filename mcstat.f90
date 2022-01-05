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
! CONTRIBUTIONS: Rohit Pappu, Hoang Tran, Albert Mao                       !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! ##################################################
! ##                                              ##
! ## analysis master routine (also prints final   ##
! ## simulation statistics to log-output)         ## 
! ##                                              ##
! ##################################################
!
!
subroutine mcstat(istep,ndump)
!
  use system
  use iounit
  use mcsums
  use energies
  use pdb
  use polyavg
  use torsn
  use movesets
  use paircorr
  use mpistuff
  use accept
  use contacts
  use sequen
  use molecule
  use aminos
  use dipolavg
  use atoms
  use diffrac
  use forces
  use fyoc
  use grandensembles
  use dssps
  use ujglobals
  use shakeetal
  use clusters
  use ems
  use wl
  use units
  
  ! Martin added 
  use martin_own
#ifdef LINK_NETCDF
  use netcdf
#endif
!
  implicit none
!
  integer freeunit,idump,ndump,lext,i,j,ttc,ihelp,ii,jj
  integer istep,modstep,btc,acnts(10),ret,rs,imol,aone
  RTYPE prtout(MAXCHI+6)
  character(1) charman
  character(7) string
  character(3) resname
  character(MAXSTRLEN) dumpfile,helpfile
  character(10) heady(10)
  logical exists,sayyes,sayno,dodcdhd
  ! Martin added
  character WWW*200000,WW_tmp*1000,my_format*10
  character*1 tab
#ifdef ENABLE_MPI
  character(3) nod
  integer tl
#endif
  tab = char(9)
  sayyes = .true.
  sayno = .false.
  aone = 1
! don't get confused with istep vs. nstep, they should agree at all times,
! the only difference being that nstep is global, while istep is passed on
!
! use output/analysis frequencies to determine what needs to be done
!
! energy output: equilibration-independent
  
  modstep = mod(nstep,xyzout)  
  if ((nstep.gt.nequil).AND.(modstep.eq.0).AND.((do_hsq.eqv..true.))) then
  !Martin : added : every time you print the structure, print the charge state
  
        WWW=''
        WW_tmp=''    
        write(WW_tmp,*) nstep
        WWW=trim(WWW)//trim(WW_tmp)


        WW_tmp=''
        write(WW_tmp,*) tab,seq_q_state(1:n_res) ! Changed to limit the size of the output filewill bot ouput the salt anymore
        WWW=trim(WWW)//trim(WW_tmp)


        write(my_format,'(A2,I0,A1)')'(A',len(trim(WWW)),')'
        if (len(trim(WWW)).ge.200000) then 
            print *,'The size of the output is higher than the container'
            stop
        end if
        write(q_state_f,my_format) trim(WWW)
        flush(q_state_f)
   end if 
  
  modstep = mod(nstep,enout)
  if (((modstep.eq.0).and.(do_hsq.eqv..false.)).or.((do_hsq.eqv..true.).and.(modstep.eq.0).and.(nstep.gt.nequil))) then
    if (do_accelsim.EQV..true.) then
      esterms(:) = esterms(:) - hmjam%boosts(:,1)
      esave = sum(esterms)
      if ((use_dyn.EQV..true.).AND.(hmjam%isin(6).EQV..true.).AND.(use_POLAR.EQV..false.)) then
        esterms(6) = esterms(6) + hmjam%boosts(6,1)
        esterms(9) = esterms(9) - hmjam%boosts(6,1)
      end if
      if ((use_dyn.EQV..true.).AND.(hmjam%isin(6).EQV..true.).AND.(use_IPP.EQV..false.)) then
        if (use_WCA.EQV..true.) then
          esterms(1) = esterms(1) + hmjam%boosts(1,1)
          esterms(5) = esterms(5) - hmjam%boosts(1,1)
        else if (use_attLJ.EQV..true.) then
          esterms(1) = esterms(1) + hmjam%boosts(1,1)
          esterms(3) = esterms(3) - hmjam%boosts(1,1)
        end if
      end if
    end if
 14     format(i12,1x,22(g14.7,1x))
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      write(iene,14) istep,esave,esterms(1),esterms(2),esterms(3),&
 &                 esterms(4),esterms(5),esterms(6),esterms(7),&
 &                 esterms(8),esterms(9),esterms(10),esterms(11),&
 &                 esterms(12),esterms(13),esterms(14),esterms(15),&
 &                 esterms(16),esterms(17),esterms(18),esterms(19),&
 &                 esterms(20),esterms(21)
    else if (use_MPIAVG.EQV..true.) then
      if (myrank.eq.0) then
        call MPI_AVGWriteEn()
      else
        call MPI_AVGSendEn()
      end if
      if (mpi_cnt_en.eq.(mpi_nodes-1)) then
        mpi_cnt_en = 0
      else
        mpi_cnt_en = mpi_cnt_en + 1
      end if
    end if
#else

    write(iene,14) istep,esave,esterms(1),esterms(2),esterms(3),&
 &                 esterms(4),esterms(5),esterms(6),esterms(7),&
 &                 esterms(8),esterms(9),esterms(10),esterms(11),&
 &                 esterms(12),esterms(13),esterms(14),esterms(15),&
 &                 esterms(16),esterms(17),esterms(18),esterms(19),&
 &                 esterms(20),esterms(21)
#endif
    if (do_accelsim.EQV..true.) then
      esterms(:) = esterms(:) + hmjam%boosts(:,1)
      esave = sum(esterms)
      if ((use_dyn.EQV..true.).AND.(hmjam%isin(6).EQV..true.).AND.(use_POLAR.EQV..false.)) then
        esterms(6) = esterms(6) - hmjam%boosts(6,1)
        esterms(9) = esterms(9) + hmjam%boosts(6,1)
      end if
      if ((use_dyn.EQV..true.).AND.(hmjam%isin(6).EQV..true.).AND.(use_IPP.EQV..false.)) then
        if (use_WCA.EQV..true.) then
          esterms(1) = esterms(1) - hmjam%boosts(1,1)
          esterms(5) = esterms(5) + hmjam%boosts(1,1)
        else if (use_attLJ.EQV..true.) then
          esterms(1) = esterms(1) - hmjam%boosts(1,1)
          esterms(3) = esterms(3) + hmjam%boosts(1,1)
        end if
      end if
    end if
  end if
!
! acceptance output: equilibration-independent
  if (((dyn_mode.eq.1).OR.(((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))&
 &.AND.(in_dyncyc.EQV..false.))).AND.(pdb_analyze.EQV..false.)) then
    modstep = mod(istep,accout)
    if (modstep.eq.0) then
 24     format(i12,1x,19i11)
#ifdef ENABLE_MPI
      if (use_REMC.EQV..true.) then
        write(iacc,24) istep,&
 &acc%naccept,acc%nfy,acc%nomega,acc%ndjcr,acc%ndocr,acc%nsjcr,acc%nujcr,acc%nchi,&
 &acc%nnuc,acc%npucker,acc%nnuccr,&
 &acc%nrb,acc%nrot,acc%ntrans,acc%nclurb,&
 &acc%ninsert,acc%ndelete,acc%nidentitychange,&
 &acc%nre
      else if (use_MPIAVG.EQV..true.) then
        if (myrank.eq.0) then
          call MPI_AVGWriteAcc()
        else
          call MPI_AVGSendAcc()
        end if
      end if
#else
        write(iacc,24) istep,&
 &acc%naccept,acc%nfy,acc%nomega,acc%ndjcr,acc%ndocr,acc%nsjcr,acc%nujcr,acc%nchi,&
 &acc%nnuc,acc%npucker,acc%nnuccr,&
 &acc%nrb,acc%nrot,acc%ntrans,acc%nclurb,&
 &acc%ninsert,acc%ndelete,acc%nidentitychange,&
 &acc%nre
#endif
    end if
  end if
!
! ensemble variables output: equilibration-independent
  if ((pdb_analyze.EQV..false.).AND.(((dyn_mode.ge.2).AND.(dyn_mode.le.4)).OR.&
 &(((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)).AND.(in_dyncyc.EQV..true.)))) then
    modstep = mod(istep,ensout)
    if (modstep.eq.0) then
#ifdef ENABLE_MPI
      if (use_REMC.EQV..true.) then
        if ((ens%flag.eq.1).OR.(ens%flag.eq.2)) then
          if (fycxyz.eq.1) then
            write(iens,14) istep,ens%insU,ens%insK,ens%insK+ens%insU,ens%insT,ens%insR(1),ens%insK2,ens%insK2+ens%insU
          else
            write(iens,14) istep,ens%insU,ens%insK,ens%insK+ens%insU,ens%insT,ens%insR(1)
          end if
        else if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
          if (fycxyz.eq.1) then
            write(iens,14) istep,ens%insU,ens%insK,ens%insK+ens%insU,ens%insT,ens%insR(1),ens%insK2,ens%insK2+ens%insU,&
 &ens%insK+ens%insU+bnd_pV,ens%insP,ens%insV,bnd_pV
          else
            write(iens,14) istep,ens%insU,ens%insK,ens%insK+ens%insU,ens%insT,ens%insR(1),&
 &ens%insK+ens%insU+bnd_pV,ens%insP,ens%insV,bnd_pV
          end if
        end if
      else if (use_MPIAVG.EQV..true.) then
        if (myrank.eq.0) then
          call MPI_AVGWriteEns()
        else
          call MPI_AVGSendEns()
        end if
      end if
#else
      if ((ens%flag.eq.1).OR.(ens%flag.eq.2)) then
        if (fycxyz.eq.1) then
          write(iens,14) istep,ens%insU,ens%insK,ens%insK+ens%insU,ens%insT,ens%insR(1),ens%insK2,ens%insK2+ens%insU
        else
          write(iens,14) istep,ens%insU,ens%insK,ens%insK+ens%insU,ens%insT,ens%insR(1)
        end if
      else if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
        if (fycxyz.eq.1) then
          write(iens,14) istep,ens%insU,ens%insK,ens%insK+ens%insU,ens%insT,ens%insR(1),ens%insK2,ens%insK2+ens%insU,&
 &ens%insK+ens%insU+bnd_pV,ens%insP,ens%insV,bnd_pV
        else
          write(iens,14) istep,ens%insU,ens%insK,ens%insK+ens%insU,ens%insT,ens%insR(1),&
 &ens%insK+ens%insU+bnd_pV,ens%insP,ens%insV,bnd_pV
        end if
      end if
#endif
    end if
  end if
!
! bla
!  call TRPZ_analysis()!DPDP_analysis()
!  call bla(ndump)
!
! if doing Wang-Landau, keep relevant RC values up-to-date
  if (do_wanglandau.EQV..true.) then
    call wl_keep_rc_updated()
  end if
!
! instantaneous torsions: equilibration-independent
  modstep = mod(nstep,torout)
  if (modstep .eq. 0) then
    call torsion()
  end if
! polymer statistics
  modstep = mod(nstep,polcalc)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
    call do_poly()
  end if
! particle fluctuation statistics
  if ((ens%flag.ge.5).AND.(ens%flag.le.6)) then
    modstep = mod(nstep,particlenumcalc)
    if ((istep.gt.nequil).AND.(modstep.eq.0)) then
      call particleflucstat()
    end if
  end if
! perform holes calculation
  modstep = mod(nstep,holescalc)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
    call holes()
  end if
! hydrodynamic radius 
  modstep = mod(nstep,rhcalc)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
    call hydroradavg()
  end if
! Ramachandran analysis
  modstep = mod(nstep,angcalc)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
    call do_rama()
  end if
! internal coordinate analysis
  modstep = mod(nstep,intcalc)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
    call do_inthists()
  end if
! backbone segment analysis
  modstep = mod(nstep,segcalc)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
     call bb_segments()
  end if
! DSSP analysis
  modstep = mod(nstep,dsspcalc)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
     call do_dssp()
  end if
! contact analysis
  modstep = mod(nstep,contactcalc)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
     call rescontacts()
  end if
! pair correlation analysis
  modstep = mod(nstep,pccalc)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
     call do_rbc_pc()
     if (do_amidepc.EQV..true.) call do_amid_pc()
     if (gpc%nos.gt.0) call do_general_pc()
  end if
! SAV analysis (note that instantaneous output is coupled to averaging
! so we don't allow it pre-equil)
  modstep = mod(nstep,savcalc)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
    call do_sav()
  end if
! instantaneous system polymer stuff: equilibration-independent
  modstep = mod(nstep,polout)
  if (modstep.eq.0) then
    call polymer()
  end if
! covariance analysis
  modstep = mod(nstep,covcalc)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
    call do_covar_int()
  end if
! dipole analysis
  modstep = mod(nstep,dipcalc)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
    call do_dipoles()
  end if
! diffraction analysis
  modstep = mod(nstep,diffrcalc)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
    call do_diffraction()
  end if
! EM analysis
  modstep = mod(nstep,emcalc)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
    call do_em_density_map()
  end if
! data collection for structural clustering
  modstep = mod(nstep,cstorecalc)
#ifdef ENABLE_MPI
  if ((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.)) then
    if (modstep.eq.0) then
      call store_for_clustering()
    end if
  else
    if ((istep.gt.nequil).AND.(modstep.eq.0)) then
      call store_for_clustering()
    end if
  end if
#else
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
    call store_for_clustering()
  end if
#endif
! LC of torsions analysis
  modstep = mod(nstep,torlccalc)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
    call collect_lc_tor()
  end if
! pH (titration state) output
  modstep = mod(nstep,phout)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
    call prt_phstate()
  end if
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
!   RE overlap analysis
    modstep = mod(nstep,re_olcalc)
    if ((istep.gt.nequil).AND.(modstep.eq.0)) then
      call MPI_REOverlap()
    end if
  end if
#endif
! restart files are equilibration-independent (of course), but do not work for trajectory analysis runs
  if (pdb_analyze.EQV..false.) then
    modstep = mod(nstep,rstout)
#ifdef ENABLE_MPI
    if ((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.)) then
      if (modstep.eq.0) then
        if ((mod(nstep,re_freq).eq.0).AND.(nstep.lt.nsim)) modstep = 1
      end if
    end if
#endif
    if (modstep.eq.0) then
      call prt_restart()
    end if
  end if
! frame weights for ELS
  if (do_accelsim.EQV..true.) then
    modstep = mod(nstep,hmjam%prtfrmwts)
    if (modstep.eq.0) then
      call els_prt_wts()
    end if
  end if
!
! structure output (note that this instantaneous output is still
! equilibration-dependent (since it's most likely for post-processing)
! might change in the future, though
  modstep = mod(nstep,xyzout)
  if ((istep.gt.nequil).AND.(modstep.eq.0)) then
!
#ifdef ENABLE_MPI
     if ((use_REMC.EQV..true.).OR.((use_MPIAVG.EQV..true.).AND.(force_singlexyz.EQV..true.))) then
       ndump = ndump + 1
       lext = 5
       tl = 3
       call int2str(ndump,string,lext)
       call int2str(myrank,nod,tl)
!
       if (xyzmode.eq.1) then
         dumpfile = 'N_'//nod(1:tl)//'_'//basename(1:bleng)//'_'//string(1:lext)//'.arc'
         call strlims(dumpfile,ii,jj)
         inquire(file=dumpfile(ii:jj),exist=exists)
         if (exists.EQV..true.) then
           idump = freeunit()
           open(unit=idump,file=dumpfile(ii:jj),status='old')
           close(unit=idump,status='delete')
         end if
         idump=freeunit()
         open(unit=idump,file=dumpfile(ii:jj),status='new')
         if (ndump.eq.1) then
           helpfile = 'N_'//nod(1:tl)//'_'//basename(1:bleng)//'_VIS.vmd'
           call strlims(helpfile,ii,jj)
           inquire(file=helpfile(ii:jj),exist=exists)
           if (exists.EQV..true.) then
             ihelp = freeunit()
             open(unit=ihelp,file=helpfile(ii:jj),status='old')
             close(unit=ihelp,status='delete')
             open(unit=ihelp,file=helpfile(ii:jj),status='new')
             write(ilog,*) 'WARNING: VMD-visualization file already &
 &exists (',helpfile(ii:jj),'). Overwriting (might not work properly when s&
 &imulation appended).'
           else
             ihelp=freeunit()
             open(unit=ihelp,file=helpfile(ii:jj),status='new')
           end if
           call FMSMC_pdb_helper(ihelp,ndump)
           close(unit=ihelp)
         end if
         call tinkerxyz(idump)
         close(unit=idump)
       end if
!
       if (xyzmode.eq.2) then
         if (pdb_writemode.eq.1) then
           dumpfile = 'N_'//nod(1:tl)//'_'//basename(1:bleng)//'_'//string(1:lext)//'.pdb'
           call strlims(dumpfile,ii,jj)
           inquire(file=dumpfile(ii:jj),exist=exists)
           if (exists.EQV..true.) then
             idump = freeunit()
             open(unit=idump,file=dumpfile(ii:jj),status='old')
             close(unit=idump,status='delete')
           end if
           idump=freeunit()
           open(unit=idump,file=dumpfile(ii:jj),status='new')
         else
           dumpfile = 'N_'//nod(1:tl)//'_'//basename(1:bleng)//'_traj.pdb'
           call strlims(dumpfile,ii,jj)
           inquire(file=dumpfile(ii:jj),exist=exists)
           if (exists.EQV..true.) then
             idump = freeunit()
             open(unit=idump,file=dumpfile(ii:jj),status='old',position='append')
             if (ndump.eq.1) then
               write(ilog,*) 'WARNING: PDB-trajectory file already e&
 &xists (',dumpfile(ii:jj),'). Appending to previous simulation?!'
             end if
           else
             if (ndump.ne.1) then
               write(ilog,*) 'WARNING: PDB-trajectory file has been &
 &moved or removed (step ',nstep,'). Creating new file.'
             end if
             idump=freeunit()
             open(unit=idump,file=dumpfile(ii:jj),status='new')
           end if
         end if
         if (ndump.eq.1) then
           helpfile = 'N_'//nod(1:tl)//'_'//basename(1:bleng)//'_VIS.vmd'
           call strlims(helpfile,ii,jj)
           inquire(file=helpfile(ii:jj),exist=exists)
           if (exists.EQV..true.) then
             ihelp = freeunit()
             open(unit=ihelp,file=helpfile(ii:jj),status='old')
             close(unit=ihelp,status='delete')
             open(unit=ihelp,file=helpfile(ii:jj),status='new')
             write(ilog,*) 'WARNING: VMD-visualization file already &
 &exists (',helpfile(ii:jj),'). Overwriting (might not work properly when s&
 &imulation appended).'
           else
             ihelp=freeunit()
             open(unit=ihelp,file=helpfile(ii:jj),status='new')
           end if
           call FMSMC_pdb_helper(ihelp,ndump)
           close(unit=ihelp)
         end if
         call FMSMC_pdb(idump,ndump)
         close(unit=idump)
       end if
!
       if (xyzmode.eq.3) then
         dumpfile = 'N_'//nod(1:tl)//'_'//basename(1:bleng)//'_traj.dcd'
         call strlims(dumpfile,ii,jj)
         inquire(file=dumpfile(ii:jj),exist=exists)
         if (exists.EQV..true.) then
           dodcdhd = .false.
           idump = freeunit()
           open(unit=idump,file=dumpfile(ii:jj),status='old',form='unformatted',position='append')
           if (ndump.eq.1) then
             write(ilog,*) 'WARNING: DCD-trajectory file already exi&
 &sts (',dumpfile(ii:jj),'). Appending to previous simulation?!'
           end if
         else
           dodcdhd = .true.
           if (ndump.ne.1) then
             write(ilog,*) 'WARNING: DCD-trajectory file has been mo&
 &ved or removed (step ',nstep,'). Creating new file.'
           end if
           idump = freeunit()
           open(unit=idump,file=dumpfile(ii:jj),status='new',form='unformatted')
         end if
         if (ndump.eq.1) then
           helpfile = 'N_'//nod(1:tl)//'_'//basename(1:bleng)//'_VIS.vmd'
           call strlims(helpfile,ii,jj)
           inquire(file=helpfile(ii:jj),exist=exists)
           if (exists.EQV..true.) then
             ihelp = freeunit()
             open(unit=ihelp,file=helpfile(ii:jj),status='old')
             close(unit=ihelp,status='delete')
             open(unit=ihelp,file=helpfile(ii:jj),status='new')
             write(ilog,*) 'WARNING: VMD-visualization file already &
 &exists (',helpfile(ii:jj),'). Overwriting (might not work properly when s&
 &imulation appended).'
           else
             ihelp=freeunit()
             open(unit=ihelp,file=helpfile(ii:jj),status='new')
           end if
           call FMSMC_pdb_helper(ihelp,ndump)
           close(unit=ihelp)
         end if
         call FMSMC_prtdcd(idump,ndump,dodcdhd)
         close(unit=idump)
       end if
!
#ifdef LINK_XDR
       if (xyzmode.eq.4) then
         dumpfile = 'N_'//nod(1:tl)//'_'//basename(1:bleng)//'_traj.xtc'
         call strlims(dumpfile,ii,jj)
         inquire(file=dumpfile(ii:jj),exist=exists)
         if (exists.EQV..true.) then
           charman = "a"
           call xdrfopen(idump,dumpfile(ii:jj),charman,ret)
           if (ndump.eq.1) then
             write(ilog,*) 'WARNING: XTC-trajectory file already exi&
 &sts (',dumpfile(ii:jj),'). Appending to previous simulation?!'
           end if
         else
           if (ndump.ne.1) then
             write(ilog,*) 'WARNING: XTC-trajectory file has been mo&
 &ved or removed (step ',nstep,'). Creating new file.'
           end if
           charman = "w"
           call xdrfopen(idump,dumpfile(ii:jj),charman,ret)
         end if
         if (ndump.eq.1) then
           helpfile = 'N_'//nod(1:tl)//'_'//basename(1:bleng)//'_VIS.vmd'
           call strlims(helpfile,ii,jj)
           inquire(file=helpfile(ii:jj),exist=exists)
           if (exists.EQV..true.) then
             ihelp = freeunit()
             open(unit=ihelp,file=helpfile(ii:jj),status='old')
             close(unit=ihelp,status='delete')
             open(unit=ihelp,file=helpfile(ii:jj),status='new')
             write(ilog,*) 'WARNING: VMD-visualization file already &
 &exists (',helpfile(ii:jj),'). Overwriting (might not work properly when s&
 &imulation appended).'
           else
             ihelp=freeunit()
             open(unit=ihelp,file=helpfile(ii:jj),status='new')
           end if
           call FMSMC_pdb_helper(ihelp,ndump)
           close(unit=ihelp)
         end if
         call FMSMC_prtxtc(idump,ndump)
         call xdrfclose(idump,ret)
       end if
#endif
!
#ifdef LINK_NETCDF
       if (xyzmode.eq.5) then
         dumpfile = 'N_'//nod(1:tl)//'_'//basename(1:bleng)//'_traj.nc'
         call strlims(dumpfile,ii,jj)
         inquire(file=dumpfile(ii:jj),exist=exists)
         if (exists.EQV..true.) then
           dodcdhd = .false.
           idump = freeunit()
           call check_fileionetcdf( nf90_open(dumpfile(ii:jj), NF90_WRITE, idump) )
           if (ndump.eq.1) then
             write(ilog,*) 'WARNING: NetCDF-trajectory file already exi&
 &sts (',dumpfile(ii:jj),'). Appending to previous simulation?!'
             call check_netcdfappend(idump,dodcdhd,dumpfile(ii:jj))
           end if
         else
           dodcdhd = .true.
           if (ndump.ne.1) then
             write(ilog,*) 'WARNING: NetCDF-trajectory file has been mo&
 &ved or removed (step ',nstep,'). Creating new file.'
           end if
           idump = freeunit()
           call check_fileionetcdf( nf90_create(dumpfile(ii:jj), NF90_CLOBBER, idump) )
         end if
         if (ndump.eq.1) then
           helpfile = 'N_'//nod(1:tl)//'_'//basename(1:bleng)//'_VIS.vmd'
           call strlims(helpfile,ii,jj)
           inquire(file=helpfile(ii:jj),exist=exists)
           if (exists.EQV..true.) then
             ihelp = freeunit()
             open(unit=ihelp,file=helpfile(ii:jj),status='old')
             close(unit=ihelp,status='delete')
             open(unit=ihelp,file=helpfile(ii:jj),status='new')
             write(ilog,*) 'WARNING: VMD-visualization file already &
 &exists (',helpfile(ii:jj),'). Overwriting (might not work properly when s&
 &imulation appended).'
           else
             ihelp=freeunit()
             open(unit=ihelp,file=helpfile(ii:jj),status='new')
           end if
           call FMSMC_pdb_helper(ihelp,ndump)
           close(unit=ihelp)
         end if
         call FMSMC_prtnetcdf(idump,ndump,dodcdhd)
         call check_fileionetcdf( nf90_close(idump) )
       end if
#endif
!
     else if ((use_MPIAVG.EQV..true.).AND.(force_singlexyz.EQV..false.)) then
       ndump = ndump + 1
       lext = 5
       call int2str(ndump,string,lext)
       if (myrank.eq.0) then
         call MPI_AVGWaitStr()
       else
         call MPI_AVGPingStr()
       end if
!
       if (mpi_cnt_xyz.eq.myrank) then
!
         if (xyzmode.eq.1) then
           dumpfile = basename(1:bleng)//'_'//string(1:lext)//'.arc'
           call strlims(dumpfile,ii,jj)
           inquire(file=dumpfile(ii:jj),exist=exists)
           if (exists.EQV..true.) then
             idump = freeunit()
             open(unit=idump,file=dumpfile(ii:jj),status='old')
             close(unit=idump,status='delete')
           end if
           idump=freeunit()
           open(unit=idump,file=dumpfile(ii:jj),status='new')
           if (ndump.eq.1) then
             helpfile = basename(1:bleng)//'_VIS.vmd'
             call strlims(helpfile,ii,jj)
             inquire(file=helpfile(ii:jj),exist=exists)
             if (exists.EQV..true.) then
               ihelp = freeunit()
               open(unit=ihelp,file=helpfile(ii:jj),status='old')
               close(unit=ihelp,status='delete')
               open(unit=ihelp,file=helpfile(ii:jj),status='new')
               write(ilog,*) 'WARNING: VMD-visualization file alread&
 &y exists (',helpfile(ii:jj),'). Overwriting (might not work properly when&
 & simulation appended).'
             else
               ihelp=freeunit()
               open(unit=ihelp,file=helpfile(ii:jj),status='new')
             end if
             call FMSMC_pdb_helper(ihelp,ndump)
             close(unit=ihelp)
           end if
           call tinkerxyz(idump)
           close(unit=idump)
         end if
!
         if (xyzmode.eq.2) then
           if (pdb_writemode.eq.1) then
             dumpfile = basename(1:bleng)//'_'//string(1:lext)//'.pdb'
             call strlims(dumpfile,ii,jj)
             inquire(file=dumpfile(ii:jj),exist=exists)
             if (exists.EQV..true.) then
               idump = freeunit()
               open(unit=idump,file=dumpfile(ii:jj),status='old')
               close(unit=idump,status='delete')
             end if
             idump=freeunit()
             open(unit=idump,file=dumpfile(ii:jj),status='new')
           else
             dumpfile = basename(1:bleng)//'_traj.pdb'
             call strlims(dumpfile,ii,jj)
             inquire(file=dumpfile(ii:jj),exist=exists)
             if (exists.EQV..true.) then
               idump = freeunit()
!               do while (1.eq.1)
!                 inquire(file=dumpfile(ii:jj),access=testy)
!                 if (testy(1:9).eq.'UNDEFINED') exit
!               end do
               open(unit=idump,file=dumpfile(ii:jj),status='old',position='append')
               if (ndump.eq.1) then
                 write(ilog,*) 'WARNING: PDB-trajectory file already&
 & exists (',dumpfile(ii:jj),'). Appending to previous simulation?!'
               end if
             else
               if (ndump.ne.1) then
                 write(ilog,*) 'WARNING: PDB-trajectory file has bee&
 &n moved or removed (step ',nstep,'). Creating new file.'
               end if
               idump=freeunit()
               open(unit=idump,file=dumpfile(ii:jj),status='new')
             end if
           end if
           if (ndump.eq.1) then
             helpfile = basename(1:bleng)//'_VIS.vmd'
             call strlims(helpfile,ii,jj)
             inquire(file=helpfile(ii:jj),exist=exists)
             if (exists.EQV..true.) then
               ihelp = freeunit()
               open(unit=ihelp,file=helpfile(ii:jj),status='old')
               close(unit=ihelp,status='delete')
               open(unit=ihelp,file=helpfile(ii:jj),status='new')
               write(ilog,*) 'WARNING: VMD-visualization file alread&
 &y exists (',helpfile(ii:jj),'). Overwriting (might not work properly when&
 & simulation appended).'
             else
               ihelp=freeunit()
               open(unit=ihelp,file=helpfile(ii:jj),status='new')
             end if
             call FMSMC_pdb_helper(ihelp,ndump)
             close(unit=ihelp)
           end if
           call FMSMC_pdb(idump,ndump)
           close(unit=idump)
         end if
!
         if (xyzmode.eq.3) then
           dumpfile = basename(1:bleng)//'_traj.dcd'
           call strlims(dumpfile,ii,jj)
           inquire(file=dumpfile(ii:jj),exist=exists)
           if (exists.EQV..true.) then
             dodcdhd = .false.
             idump = freeunit()
             open(unit=idump,file=dumpfile(ii:jj),status='old',form='unformatted',position='append')
             if (ndump.eq.1) then
               write(ilog,*) 'WARNING: DCD-trajectory file already e&
 &xists (',dumpfile(ii:jj),'). Appending to previous simulation?!'
             end if
           else
             dodcdhd = .true.
             if (ndump.ne.1) then
               write(ilog,*) 'WARNING: DCD-trajectory file has been &
 &moved or removed (step ',nstep,'). Creating new file.'
             end if
             idump = freeunit()
             open(unit=idump,file=dumpfile(ii:jj),status='new',form='unformatted')
           end if
           if (ndump.eq.1) then
             helpfile = basename(1:bleng)//'_VIS.vmd'
             call strlims(helpfile,ii,jj)
             inquire(file=helpfile(ii:jj),exist=exists)
             if (exists.EQV..true.) then
               ihelp = freeunit()
               open(unit=ihelp,file=helpfile(ii:jj),status='old')
               close(unit=ihelp,status='delete')
               open(unit=ihelp,file=helpfile(ii:jj),status='new')
               write(ilog,*) 'WARNING: VMD-visualization file alread&
 &y exists (',helpfile(ii:jj),'). Overwriting (might not work properly when&
 & simulation appended).'
             else
               ihelp=freeunit()
               open(unit=ihelp,file=helpfile(ii:jj),status='new')
             end if
             call FMSMC_pdb_helper(ihelp,ndump)
             close(unit=ihelp)
           end if
           call FMSMC_prtdcd(idump,ndump,dodcdhd)
           close(unit=idump)
         end if
!
#ifdef LINK_XDR
         if (xyzmode.eq.4) then
           dumpfile = basename(1:bleng)//'_traj.xtc'
           call strlims(dumpfile,ii,jj)
           inquire(file=dumpfile(ii:jj),exist=exists)
           if (exists.EQV..true.) then
             charman = "a"
             call xdrfopen(idump,dumpfile(ii:jj),charman,ret)
             if (ndump.eq.1) then
               write(ilog,*) 'WARNING: XTC-trajectory file already e&
 &xists (',dumpfile(ii:jj),'). Appending to previous simulation?!'
             end if
           else
             if (ndump.ne.1) then
               write(ilog,*) 'WARNING: XTC-trajectory file has been &
 &moved or removed (step ',nstep,'). Creating new file.'
             end if
             charman = "w"
             call xdrfopen(idump,dumpfile(ii:jj),charman,ret)
           end if
           if (ndump.eq.1) then
             helpfile = basename(1:bleng)//'_VIS.vmd'
             call strlims(helpfile,ii,jj)
             inquire(file=helpfile(ii:jj),exist=exists)
             if (exists.EQV..true.) then
               ihelp = freeunit()
               open(unit=ihelp,file=helpfile(ii:jj),status='old')
               close(unit=ihelp,status='delete')
               open(unit=ihelp,file=helpfile(ii:jj),status='new')
               write(ilog,*) 'WARNING: VMD-visualization file alread&
 &y exists (',helpfile(ii:jj),'). Overwriting (might not work properly when&
 & simulation appended).'
             else
               ihelp=freeunit()
               open(unit=ihelp,file=helpfile(ii:jj),status='new')
             end if
             call FMSMC_pdb_helper(ihelp,ndump)
             close(unit=ihelp)
           end if
           call FMSMC_prtxtc(idump,ndump)
           call xdrfclose(idump,ret)
         end if
#endif
!
#ifdef LINK_NETCDF
         if (xyzmode.eq.5) then
           dumpfile = basename(1:bleng)//'_traj.nc'
           call strlims(dumpfile,ii,jj)
           inquire(file=dumpfile(ii:jj),exist=exists)
           if (exists.EQV..true.) then
             dodcdhd = .false.
             idump = freeunit()
             call check_fileionetcdf( nf90_open(dumpfile(ii:jj), NF90_WRITE, idump) )
             if (ndump.eq.1) then
               write(ilog,*) 'WARNING: NetCDF-trajectory file already e&
 &xists (',dumpfile(ii:jj),'). Appending to previous simulation?!'
               call check_netcdfappend(idump,dodcdhd,dumpfile(ii:jj))
             end if
           else
             dodcdhd = .true.
             if (ndump.ne.1) then
               write(ilog,*) 'WARNING: NetCDF-trajectory file has been &
 &moved or removed (step ',nstep,'). Creating new file.'
             end if
             idump = freeunit()
             call check_fileionetcdf( nf90_create(dumpfile(ii:jj), NF90_CLOBBER, idump) )
           end if
           if (ndump.eq.1) then
             helpfile = basename(1:bleng)//'_VIS.vmd'
             call strlims(helpfile,ii,jj)
             inquire(file=helpfile(ii:jj),exist=exists)
             if (exists.EQV..true.) then
               ihelp = freeunit()
               open(unit=ihelp,file=helpfile(ii:jj),status='old')
               close(unit=ihelp,status='delete')
               open(unit=ihelp,file=helpfile(ii:jj),status='new')
               write(ilog,*) 'WARNING: VMD-visualization file alread&
 &y exists (',helpfile(ii:jj),'). Overwriting (might not work properly when&
 & simulation appended).'
             else
               ihelp=freeunit()
               open(unit=ihelp,file=helpfile(ii:jj),status='new')
             end if
             call FMSMC_pdb_helper(ihelp,ndump)
             close(unit=ihelp)
           end if
           call FMSMC_prtnetcdf(idump,ndump,dodcdhd)
           call check_fileionetcdf( nf90_close(idump) )
         end if
#endif
!
       end if
!
#ifdef LINK_NETCDF
       if ((xyzmode.eq.5).AND.(netcdf_ids(40).ne.-1)) then
         call MPI_AVGSyncNetCDF()
       end if
#endif
!
       if (myrank.eq.0) then
         call MPI_AVGWaitStr()
       else
         call MPI_AVGPingStr()
       end if
!
       if (mpi_cnt_xyz.eq.(mpi_nodes-1)) then
         mpi_cnt_xyz = 0
       else
         mpi_cnt_xyz = mpi_cnt_xyz + 1
       end if
!
     end if
#else
     ndump = ndump + 1
     lext = 5
     call int2str(ndump,string,lext)
!
     if (xyzmode.eq.1) then
       dumpfile = basename(1:bleng)//'_'//string(1:lext)//'.arc'
       call strlims(dumpfile,ii,jj)
       inquire(file=dumpfile(ii:jj),exist=exists)
       if (exists.EQV..true.) then
         idump = freeunit()
         open(unit=idump,file=dumpfile(ii:jj),status='old')
         close(unit=idump,status='delete')
       end if
       idump=freeunit()
       open(unit=idump,file=dumpfile(ii:jj),status='new')
       if (ndump.eq.1) then
         helpfile = basename(1:bleng)//'_VIS.vmd'
         call strlims(helpfile,ii,jj)
         inquire(file=helpfile(ii:jj),exist=exists)
         if (exists.EQV..true.) then
           ihelp = freeunit()
           open(unit=ihelp,file=helpfile(ii:jj),status='old')
           close(unit=ihelp,status='delete')
           open(unit=ihelp,file=helpfile(ii:jj),status='new')
           write(ilog,*) 'WARNING: VMD-visualization file already ex&
 &ists (',helpfile(ii:jj),'). Overwriting (might not work properly when sim&
 &ulation appended).'
         else
           ihelp=freeunit()
           open(unit=ihelp,file=helpfile(ii:jj),status='new')
         end if
         call FMSMC_pdb_helper(ihelp,ndump)
         close(unit=ihelp)
       end if
       call tinkerxyz(idump)
       close(unit=idump)
      end if
!
     if (xyzmode.eq.2) then
       if (pdb_writemode.eq.1) then
         dumpfile = basename(1:bleng)//'_'//string(1:lext)//'.pdb'
         call strlims(dumpfile,ii,jj)
         inquire(file=dumpfile(ii:jj),exist=exists)
         if (exists.EQV..true.) then
           idump = freeunit()
           open(unit=idump,file=dumpfile(ii:jj),status='old')
           close(unit=idump,status='delete')
         end if
         idump=freeunit()
         open(unit=idump,file=dumpfile(ii:jj),status='new')
       else
         dumpfile = basename(1:bleng)//'_traj.pdb'
         call strlims(dumpfile,ii,jj)
         inquire(file=dumpfile(ii:jj),exist=exists)
         if (exists.EQV..true.) then
           idump = freeunit()
           open(unit=idump,file=dumpfile(ii:jj),status='old',position='append')
           if (ndump.eq.1) then
             write(ilog,*) 'WARNING: PDB-trajectory file already exi&
 &sts (',dumpfile(ii:jj),'). Appending to previous simulation?!'
           end if
         else
           if (ndump.ne.1) then
             write(ilog,*) 'WARNING: PDB-trajectory file has been mo&
 &ved or removed (step ',nstep,'). Creating new file.'
           end if
           idump=freeunit()
           open(unit=idump,file=dumpfile(ii:jj),status='new')
         end if
       end if
       if (ndump.eq.1) then
         helpfile = basename(1:bleng)//'_VIS.vmd'
         call strlims(helpfile,ii,jj)
         inquire(file=helpfile(ii:jj),exist=exists)
         if (exists.EQV..true.) then
           ihelp = freeunit()
           open(unit=ihelp,file=helpfile(ii:jj),status='old')
           close(unit=ihelp,status='delete')
           open(unit=ihelp,file=helpfile(ii:jj),status='new')
           write(ilog,*) 'WARNING: VMD-visualization file already ex&
 &ists (',helpfile(ii:jj),'). Overwriting (might not work properly when sim&
 &ulation appended).'
         else
           ihelp=freeunit()
           open(unit=ihelp,file=helpfile(ii:jj),status='new')
         end if
         call FMSMC_pdb_helper(ihelp,ndump)
         close(unit=ihelp)
       end if
       call FMSMC_pdb(idump,ndump)
       close(unit=idump)
     end if
!
     if (xyzmode.eq.3) then
       dumpfile = basename(1:bleng)//'_traj.dcd'
       call strlims(dumpfile,ii,jj)
       inquire(file=dumpfile(ii:jj),exist=exists)
       if (exists.EQV..true.) then
         dodcdhd = .false.
         idump = freeunit()
         open(unit=idump,file=dumpfile(ii:jj),status='old',form='unformatted',position='append')
         if (ndump.eq.1) then
           write(ilog,*) 'WARNING: DCD-trajectory file already exist&
 &s (',dumpfile(ii:jj),'). Appending to previous simulation?!'
         end if
       else
         dodcdhd = .true.
         if (ndump.ne.1) then
           write(ilog,*) 'WARNING: DCD-trajectory file has been move&
 &d or removed (step ',nstep,'). Creating new file.'
         end if
         idump = freeunit()
         open(unit=idump,file=dumpfile(ii:jj),status='new',form='unformatted')
       end if
       if (ndump.eq.1) then
         helpfile = basename(1:bleng)//'_VIS.vmd'
         call strlims(helpfile,ii,jj)
         inquire(file=helpfile(ii:jj),exist=exists)
         if (exists.EQV..true.) then
           ihelp = freeunit()
           open(unit=ihelp,file=helpfile(ii:jj),status='old')
           close(unit=ihelp,status='delete')
           open(unit=ihelp,file=helpfile(ii:jj),status='new')
           write(ilog,*) 'WARNING: VMD-visualization file already ex&
 &ists (',helpfile(ii:jj),'). Overwriting (might not work properly when sim&
 &ulation appended).'
         else
           ihelp=freeunit()
           open(unit=ihelp,file=helpfile(ii:jj),status='new')
         end if
         call FMSMC_pdb_helper(ihelp,ndump)
         close(unit=ihelp)
       end if
       call FMSMC_prtdcd(idump,ndump,dodcdhd)
       close(unit=idump)
     end if
!
#ifdef LINK_XDR
     if (xyzmode.eq.4) then
       dumpfile = basename(1:bleng)//'_traj.xtc'
       call strlims(dumpfile,ii,jj)
       inquire(file=dumpfile(ii:jj),exist=exists)
       if (exists.EQV..true.) then
         charman = "a"
         call xdrfopen(idump,dumpfile(ii:jj),charman,ret)
         if (ndump.eq.1) then
           write(ilog,*) 'WARNING: XTC-trajectory file already exist&
 &s (',dumpfile(ii:jj),'). Appending to previous simulation?!'
         end if
       else
         if (ndump.ne.1) then
           write(ilog,*) 'WARNING: XTC-trajectory file has been move&
 &d or removed (step ',nstep,'). Creating new file.'
         end if
         charman = "w"
         call xdrfopen(idump,dumpfile(ii:jj),charman,ret)
       end if
       if (ndump.eq.1) then
         helpfile = basename(1:bleng)//'_VIS.vmd'
         call strlims(helpfile,ii,jj)
         inquire(file=helpfile(ii:jj),exist=exists)
         if (exists.EQV..true.) then
           ihelp = freeunit()
           open(unit=ihelp,file=helpfile(ii:jj),status='old')
           close(unit=ihelp,status='delete')
           open(unit=ihelp,file=helpfile(ii:jj),status='new')
           write(ilog,*) 'WARNING: VMD-visualization file already ex&
 &ists (',helpfile(ii:jj),'). Overwriting (might not work properly when sim&
 &ulation appended).'
         else
           ihelp=freeunit()
           open(unit=ihelp,file=helpfile(ii:jj),status='new')
         end if
         call FMSMC_pdb_helper(ihelp,ndump)
         close(unit=ihelp)
       end if
       call FMSMC_prtxtc(idump,ndump)
       call xdrfclose(idump,ret)
     end if
#endif
!
#ifdef LINK_NETCDF
     if (xyzmode.eq.5) then
       dumpfile = basename(1:bleng)//'_traj.nc'
       call strlims(dumpfile,ii,jj)
       inquire(file=dumpfile(ii:jj),exist=exists)
       if (exists.EQV..true.) then
         dodcdhd = .false.
         idump = freeunit()
         call check_fileionetcdf( nf90_open(dumpfile(ii:jj), NF90_WRITE, idump) )
         if (ndump.eq.1) then
           write(ilog,*) 'WARNING: NetCDF-trajectory file already exist&
 &s (',dumpfile(ii:jj),'). Appending to previous simulation?!'
           call check_netcdfappend(idump,dodcdhd,dumpfile(ii:jj))
         end if
       else
         dodcdhd = .true.
         if (ndump.ne.1) then
           write(ilog,*) 'WARNING: NetCDF-trajectory file has been move&
 &d or removed (step ',nstep,'). Creating new file.'
         end if
         idump = freeunit()
         call check_fileionetcdf( nf90_create(dumpfile(ii:jj), NF90_CLOBBER, idump) )
       end if
       if (ndump.eq.1) then
         helpfile = basename(1:bleng)//'_VIS.vmd'
         call strlims(helpfile,ii,jj)
         inquire(file=helpfile(ii:jj),exist=exists)
         if (exists.EQV..true.) then
           ihelp = freeunit()
           open(unit=ihelp,file=helpfile(ii:jj),status='old')
           close(unit=ihelp,status='delete')
           open(unit=ihelp,file=helpfile(ii:jj),status='new')
           write(ilog,*) 'WARNING: VMD-visualization file already ex&
 &ists (',helpfile(ii:jj),'). Overwriting (might not work properly when sim&
 &ulation appended).'
         else
           ihelp=freeunit()
           open(unit=ihelp,file=helpfile(ii:jj),status='new')
         end if
         call FMSMC_pdb_helper(ihelp,ndump)
         close(unit=ihelp)
       end if
       call FMSMC_prtnetcdf(idump,ndump,dodcdhd)
       call check_fileionetcdf( nf90_close(idump) )
     end if
#endif
#endif
!
  end if
!
! analysis in final step: first the parts that are
! only meaningful when nequil < nsim
  if ((istep.eq.nsim).AND.(nsim.gt.nequil)) then
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      if(angcalc .le. nsim) then
        call printrama()
      end if
      if (intcalc.le.nsim) call prt_inthists()
      if (polcalc.le.nsim) call prt_rdhist()
      if (segcalc.le.nsim) call prt_bb_segments()
      if (dsspcalc.le.nsim) call prt_dssp()
      if (polcalc.le.nsim) call prt_pers()
      if (nsim.ge.rhcalc) then
        call prt_hydrorad()
      end if
      if (nsim.ge.contactcalc) then
        call prt_rescontacts()
      end if
      if (nsim.ge.savcalc) then
        call prt_sav()
      end if
      if (nsim.ge.pccalc) then
        if (do_amidepc.EQV..true.) call prt_amid_pc()
        call prt_rbc_pc()
        if (gpc%nos.gt.0) call prt_general_pc()
      end if
      if (nsim.ge.dipcalc) then
        call prt_dipoles()
      end if
      if (nsim.ge.torlccalc) then
        call prt_lc_tor()
      end if
      if (nsim.ge.diffrcalc) then
        call prt_diffraction()
      end if
      if (nsim.ge.emcalc) then
        call prt_emmap(emgrid,aone)
      end if
      if (nsim.ge.cstorecalc) then
        call do_clustering()
      end if
      if (nsim.ge.re_olcalc) then
        call MPI_WriteREOlap()
      end if
    else if (use_MPIAVG.EQV..true.) then
      if (myrank.eq.0) then
        call MPI_AVGCollect()
      else
        call MPI_AVGSend()
      end if
    end if
#else
    if(angcalc .le. nsim) then
       call printrama()
    end if
    if (intcalc.le.nsim) call prt_inthists()
    if (segcalc.le.nsim) call prt_bb_segments()
    if (dsspcalc.le.nsim) call prt_dssp()
    if (polcalc.le.nsim) call prt_pers()
    if (polcalc.le.nsim) call prt_rdhist()
    if (((ens%flag.eq.5).OR.(ens%flag.eq.6)).AND.&
 &      (nsim.ge.particlenumcalc)) then
      call prt_particlenumhistogram()
    end if
    if (nsim.ge.rhcalc) then
      call prt_hydrorad()
    end if
    if (nsim.ge.contactcalc) then
      call prt_rescontacts()
    end if
    if (nsim.ge.pccalc) then
      call prt_rbc_pc()
      if (do_amidepc.EQV..true.) call prt_amid_pc()
      if (gpc%nos.gt.0) call prt_general_pc()
    end if
    if (nsim.ge.dipcalc) then
      call prt_dipoles()
    end if
    if (nsim.ge.savcalc) then
      call prt_sav()
    end if
    if (nsim.ge.torlccalc) then
      call prt_lc_tor()
    end if
    if (nsim.ge.diffrcalc) then
      call prt_diffraction()
    end if
    if (nsim.ge.emcalc) then
      call prt_emmap(emgrid,aone)
    end if
    if (nsim.ge.cstorecalc) then
      call do_clustering()
    end if
#endif
  end if
! and now the parts that are always worth printing out
  if (istep.eq.nsim) then
    write(ilog,*)
    write(ilog,*)
    write(ilog,*) ' ------------------ Final Processing ------------------'
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      if ((myrank.eq.0).AND.(nsim.ge.re_freq)) then
        call MPI_REWritePB()
      end if
    end if
#endif
 28     format(a3,':',i4,10(1x,i10))
 281    format(8x,10(1x,a10))
 282    format(22x,10(1x,a10))
 29     format(i4,': Res. ',i5,'-',i5,1x,i10,1x,i10)
 283    format(a20,': ',i10,' / ',i10,' = ',g12.5,'%')
 284    format(a35,': ',g14.7,' Steps (from ',i10,' Segments)')
 35     format('Frac. of Ulmschneider-Jorgensen moves w/o closure  : ',g10.4,'%')
 36     format('Frac. of non-w Dinner-Ulmschneider moves w/o clos. : ',g10.4,'%')
 37     format('Frac. of w-Dinner-Ulmschneider moves w/o closure   : ',g10.4,'%')
 38     format('Frac. of nucl. acid CR moves without closure       : ',g10.4,'%')
 39     format('Avg. # of closure attmpts. per U.-J. CR-moves      : ',g10.4)
 40     format('Avg. # of closure attmpts. per non-w-D.-U. CR-moves: ',g10.4)
 41     format('Avg. # of closure attmpts. per w-D.-U. CR-moves    : ',g10.4)
 42     format('Avg. # of closure attmpts. per nucl. acid CR moves : ',g10.4)

!
    if (((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)).AND.(pdb_analyze.EQV..false.)) then
      write(ilog,*)
      write(ilog,*) 'Bulk acceptance rates:'
      if (mvcnt%nrb.gt.0) write(ilog,283) 'RIGID-BODY',acc%nrb,mvcnt%nrb,100.0*acc%nrb/(1.0*mvcnt%nrb)
      if (mvcnt%nclurb.gt.0) write(ilog,283) 'CLUSTER',acc%nclurb,mvcnt%nclurb,100.0*acc%nclurb/(1.0*mvcnt%nclurb)
      if (mvcnt%ntrans.gt.0) write(ilog,283) 'RIGID TRANSLATION',acc%ntrans,mvcnt%ntrans,100.0*acc%ntrans/(1.0*mvcnt%ntrans)
      if (mvcnt%nrot.gt.0) write(ilog,283) 'RIGID ROTATION',acc%nrot,mvcnt%nrot,100.0*acc%nrot/(1.0*mvcnt%nrot)
      if (mvcnt%nchi.gt.0) write(ilog,283) 'SIDECHAIN',acc%nchi,mvcnt%nchi,100.0*acc%nchi/(1.0*mvcnt%nchi)
      if (mvcnt%nph.gt.0) write(ilog,283) 'TITRATION (PH)',acc%nph,mvcnt%nph,100.0*acc%nph/(1.0*mvcnt%nph)
      if (mvcnt%nomega.gt.0) write(ilog,283) 'OMEGA',acc%nomega,mvcnt%nomega,100.0*acc%nomega/(1.0*mvcnt%nomega)
      if (mvcnt%nfy.gt.0) write(ilog,283) 'PHI/PSI (PIVOT)',acc%nfy,mvcnt%nfy,100.0*acc%nfy/(1.0*mvcnt%nfy)
      if (mvcnt%nfyc.gt.0) write(ilog,283) 'PHI/PSI/CHI',acc%nfyc,mvcnt%nfyc,100.0*acc%nfyc/(1.0*mvcnt%nfyc)
      if (mvcnt%npucker.gt.0) write(ilog,283) '5-RING PUCKERING',acc%npucker,mvcnt%npucker,100.0*acc%npucker/(1.0*mvcnt%npucker)
      if (mvcnt%nnuc.gt.0) write(ilog,283) 'NUCLEIC ACID',acc%nnuc,mvcnt%nnuc,100.0*acc%nnuc/(1.0*mvcnt%nnuc)
      if (mvcnt%nother.gt.0) write(ilog,283) 'SINGLE PIVOT (OTHER)',acc%nother,mvcnt%nother,100.0*acc%nother/(1.0*mvcnt%nother)
      if (mvcnt%nsjcr.gt.0) write(ilog,283) 'CONCERTED ROT. (SJ)',acc%nsjcr,mvcnt%nsjcr,100.0*acc%nsjcr/(1.0*mvcnt%nsjcr)
      if (mvcnt%nujcr.gt.0) write(ilog,283) 'CONCERTED ROT. (UJ)',acc%nujcr,mvcnt%nujcr,100.0*acc%nujcr/(1.0*mvcnt%nujcr)
      if (mvcnt%ndjcr.gt.0) write(ilog,283) 'CONCERTED ROT. (DJ)',acc%ndjcr,mvcnt%ndjcr,100.0*acc%ndjcr/(1.0*mvcnt%ndjcr)
      if (mvcnt%ndocr.gt.0) write(ilog,283) 'CONCERTED ROT. (DO)',acc%ndocr,mvcnt%ndocr,100.0*acc%ndocr/(1.0*mvcnt%ndocr)
      if (mvcnt%nnuccr.gt.0) write(ilog,283) 'CONCERTED ROT. (NC)',acc%nnuccr,mvcnt%nnuccr,100.0*acc%nnuccr/(1.0*mvcnt%nnuccr)
      if (mvcnt%ninsert.gt.0) write(ilog,283) 'INSERTION',acc%ninsert,mvcnt%ninsert,100.0*acc%ninsert/(1.0*mvcnt%ninsert)
      if (mvcnt%ndelete.gt.0) write(ilog,283) 'DELETION',acc%ndelete,mvcnt%ndelete,100.0*acc%ndelete/(1.0*mvcnt%ndelete)
      if (mvcnt%nidentitychange.gt.0) write(ilog,283) 'ID-CHANGE',acc%nidentitychange,&
 & mvcnt%nidentitychange,100.0*acc%nidentitychange/(1.0*mvcnt%nidentitychange)
      if (mvcnt%nre.gt.0) write(ilog,283) 'REPLICA EXCHANGE',acc%nre,mvcnt%nre,100.0*acc%nre/(1.0*mvcnt%nre)
      btc = mvcnt%nrb + mvcnt%nrot + mvcnt%nclurb + mvcnt%ntrans + mvcnt%nchi&
 & + mvcnt%nre + mvcnt%nidentitychange + mvcnt%ndelete + mvcnt%ninsert + mvcnt%nsjcr&
 & + mvcnt%nujcr + mvcnt%ndjcr + mvcnt%ndocr  + mvcnt%nfyc + mvcnt%nfy + mvcnt%nomega&
 & + mvcnt%nnuc + mvcnt%nnuccr + mvcnt%nph + mvcnt%npucker + mvcnt%nother
      if (dyn_mode.eq.5) then
        btc = btc + mvcnt%nmd
        if (btc.gt.0) then
          write(ilog,*)
          write(ilog,*) 'Hybrid MC/MD statistics:'
          write(ilog,283) 'DYNAMICS MOVES',mvcnt%nmd,btc,100.0*mvcnt%nmd/(1.0*btc)
          write(ilog,283) 'MONTE CARLO MOVES',btc-mvcnt%nmd,btc,100.0*(btc-mvcnt%nmd)/(1.0*btc)
          write(ilog,*)
          if (mvcnt%ndynseg.gt.0) write(ilog,284) 'AVERAGE DYNAMICS SEGMENT LENGTH',&
 & mvcnt%avgdynseglen/(1.0*mvcnt%ndynseg),mvcnt%ndynseg
          if (mvcnt%nmcseg.gt.0) write(ilog,284) 'AVERAGE MONTE CARLO SEGMENT LENGTH',&
 & mvcnt%avgmcseglen/(1.0*mvcnt%nmcseg),mvcnt%nmcseg
          write(ilog,*)
        end if
      else if (dyn_mode.eq.7) then
        btc = btc + mvcnt%nld
        if (btc.gt.0) then
          write(ilog,*)
          write(ilog,*) 'Hybrid MC/LD statistics:'
          write(ilog,283) 'LANGEVIN MOVES',mvcnt%nld,btc,100.0*mvcnt%nld/(1.0*btc)
          write(ilog,283) 'MONTE CARLO MOVES',btc-mvcnt%nld,btc,100.0*(btc-mvcnt%nld)/(1.0*btc)
          write(ilog,*)
         if (mvcnt%ndynseg.gt.0) write(ilog,284) 'AVERAGE LANGEVIN SEGMENT LENGTH',&
 & mvcnt%avgdynseglen/(1.0*mvcnt%ndynseg),mvcnt%ndynseg
          if (mvcnt%nmcseg.gt.0) write(ilog,284) 'AVERAGE MONTE CARLO SEGMENT LENGTH',&
 & mvcnt%avgmcseglen/(1.0*mvcnt%nmcseg),mvcnt%nmcseg
          write(ilog,*)
        end if
      else if (dyn_mode.eq.8) then
        btc = btc + mvcnt%nbd
        if (btc.gt.0) then
          write(ilog,*)
          write(ilog,*) 'Hybrid MC/BD statistics:'
          write(ilog,283) 'BROWNIAN MOVES',mvcnt%nbd,btc,100.0*mvcnt%nbd/(1.0*btc)
          write(ilog,283) 'MONTE CARLO MOVES',btc-mvcnt%nbd,btc,100.0*(btc-mvcnt%nbd)/(1.0*btc)
          write(ilog,*)
         if (mvcnt%ndynseg.gt.0) write(ilog,284) 'AVERAGE BROWNIAN SEGMENT LENGTH',&
 & mvcnt%avgdynseglen/(1.0*mvcnt%ndynseg),mvcnt%ndynseg
          if (mvcnt%nmcseg.gt.0) write(ilog,284) 'AVERAGE MONTE CARLO SEGMENT LENGTH',&
 & mvcnt%avgmcseglen/(1.0*mvcnt%nmcseg),mvcnt%nmcseg
          write(ilog,*)
        end if
      end if
      if ((wld%stepnum.gt.0).AND.(mc_acc_crit.eq.3)) then
        call wl_sim_finished()
      end if
      if ((rigidfreq.lt.1.0).AND.(particleflucfreq.lt.1.0)) then
        write(ilog,*)
        write(ilog,*) 'Residue-resolved acceptance counts (int. moves):'
        btc = 0
        if ((use_globmoves.EQV..true.).OR.((have_chi.EQV..true.).AND.(phfreq.lt.1.0))) then
          btc = btc + 1
          heady(btc) = ' SIDECHAIN'
        end if
        if ((have_cr.EQV..true.).OR.(have_nuccr.EQV..true.)) then
          btc = btc + 1
          heady(btc) = ' CON. ROT.'
        end if
        if (have_omega.EQV..true.) then
          btc = btc + 1
          heady(btc) = '     OMEGA'
        end if
        if ((have_nuc.EQV..true.).AND.(nuccrfreq.lt.1.0).AND.(nucpuckfreq.lt.1.0)) then
          btc = btc + 1
          heady(btc) = ' NUC. ACID'
        end if
        if ((have_pucker.EQV..true.).OR.(have_nucpuck.EQV..true.)) then
          btc = btc + 1
          heady(btc) = ' PUCKERING'
        end if
        if (have_pivot.EQV..true.) then
          btc = btc + 1
          heady(btc) = '     PIVOT'
        end if
        if (have_other.EQV..true.) then
          btc = btc + 1
          heady(btc) = '     OTHER'
        end if
        write(ilog,281) (heady(j),j=1,btc)
        do i=1,nseq
          btc = 0
          if ((use_globmoves.EQV..true.).OR.((have_chi.EQV..true.).AND.(phfreq.lt.1.0))) then
            btc = btc + 1
            acnts(btc) = acc%chi(i)
          end if
          if ((have_cr.EQV..true.).OR.(have_nuccr.EQV..true.)) then
            btc = btc + 1
            acnts(btc) = acc%cr(i)
          end if
          if (have_omega.EQV..true.) then
            btc = btc + 1
            acnts(btc) = acc%omega(i)
          end if
          if ((have_nuc.EQV..true.).AND.(nuccrfreq.lt.1.0).AND.(nucpuckfreq.lt.1.0)) then
            btc = btc + 1
            acnts(btc) = acc%nuc(i)
          end if
          if ((have_pucker.EQV..true.).OR.(have_nucpuck.EQV..true.)) then
            btc = btc + 1
            acnts(btc) = acc%pucker(i)
          end if
          if (have_pivot.EQV..true.) then
            btc = btc + 1
            acnts(btc) = acc%fy(i)
          end if
          if (have_other.EQV..true.) then
            btc = btc + 1
            acnts(btc) = acc%other(i)
          end if
          resname = amino(seqtyp(i))
!        if (use_globmoves.EQV..true.) then
!          write(ilog,28) resname,i,fyresacc(i)
!        else
          write(ilog,28) resname,i,(acnts(j),j=1,btc)
!        end if
        end do
        write(ilog,*)
        if (mvcnt%nujcr.gt.0) then
          write(ilog,35) (100.0*uj_deadcnt)/(1.0*mvcnt%nujcr)
!          write(ilog,39) (uj_totcnt)/(1.0*mvcnt%nujcr)
          write(ilog,*)
        end if
        if (mvcnt%ndjcr.gt.0) then
          write(ilog,36) (100.0*dj_deadcnt)/(1.0*mvcnt%ndjcr)
          if (torcrmode.eq.2) write(ilog,40) (dj_totcnt)/(1.0*mvcnt%ndjcr)
          write(ilog,*)
        end if
        if (mvcnt%ndocr.gt.0) then
          write(ilog,37) (100.0*do_deadcnt)/(1.0*mvcnt%ndocr)
          if (torcrmode.eq.2) write(ilog,41) (do_totcnt)/(1.0*mvcnt%ndocr)
          write(ilog,*)
        end if
        if (mvcnt%nnuccr.gt.0) then
          write(ilog,38) (100.0*nc_deadcnt)/(1.0*mvcnt%nnuccr)
          if (torcrmode.eq.2) write(ilog,42) (nc_totcnt)/(1.0*mvcnt%nnuccr)
          write(ilog,*)
        end if
      end if
      if ((rigidfreq.gt.0.0).OR.(particleflucfreq.gt.0.0)) then
        write(ilog,*)
        write(ilog,*) 'Molecule-resolved acceptance counts:'
        btc = 0
        if (have_particlefluc.EQV..true.) then
          btc = btc + 1
          heady(btc) = ' (S)GCMC  '
        end if
        if (have_rigid.EQV..true.) then
          btc = btc + 1
          heady(btc) = 'RIGID-BODY'
        end if
        write(ilog,282) (heady(j),j=1,btc)
        do i=1,nmol
          btc = 0
          if (have_particlefluc.EQV..true.) then
            btc = btc + 1
            acnts(btc) = acc%insert(i)+acc%delete(i)+acc%permute(i)
          end if
          if (have_rigid.EQV..true.) then
            btc = btc + 1
            acnts(btc) = acc%rigid(i)
          end if
          write(ilog,29) i,rsmol(i,1),rsmol(i,2),(acnts(j),j=1,btc)
        end do
      end if
    end if
!
    if ((nequil.ge.nsim).OR.(dyn_mode.eq.6).OR.(dyn_mode.eq.1).OR.(ens%avgcnt.eq.0)) then
      write(ilog,*)
!     nothing else to do
      return
    else
 57   format(2x,a,g14.7)
 61   format(2x,a,i10)
 58   format(a3,a,10(g12.5,1x))
 59   format(i5,1x,i5,a,3(g12.5,1x))
 62   format(i10,1x,a,100000(g12.5,1x))
 60   format(a,3(g12.5,1x))
!
      write(ilog,*)
      write(ilog,*) 'Average thermal variables:'
      write(ilog,57) 'Temperature in K            :  ',ens%avgT/(1.0*ens%avgcnt)
      write(ilog,57) 'Kinetic Energy in kcal/mol  :  ',ens%avgK/(1.0*ens%avgcnt)
      write(ilog,57) 'K.E. Stand. Dev. in kcal/mol:  ',sqrt(ens%avgR(6)/(1.0*ens%avgcnt) - (ens%avgK/(1.0*ens%avgcnt))**2)
      if (fycxyz.eq.1) then
        write(ilog,57) 'Cart. kin. en. in kcal/mol  :  ',ens%avgK2/(1.0*ens%avgcnt)
        write(ilog,57) 'CKE Stand. Dev. in kcal/mol :  ',sqrt(ens%avgR(7)/(1.0*ens%avgcnt) - (ens%avgK2/(1.0*ens%avgcnt))**2)
        write(ilog,57) 'MUD of K.E.s in kcal/mol    :  ',ens%avgR(4)/(1.0*ens%avgcnt)
      else if (ens%avgR(8).gt.0.0) then
        write(ilog,57) 'Int. trans. temp. in K      :  ',ens%avgR(8)/(1.0*ens%avgcnt)
        write(ilog,57) 'Int. rot. temp. in K        :  ',ens%avgR(9)/(1.0*ens%avgcnt)
      end if
      write(ilog,57) 'Potential energy in kcal/mol:  ',ens%avgU/(1.0*ens%avgcnt)
      write(ilog,57) 'P.E. Stand. Dev. in kcal/mol:  ',sqrt(ens%avgR(5)/(1.0*ens%avgcnt) - (ens%avgU/(1.0*ens%avgcnt))**2)
      write(ilog,57) 'T.E. Stand. Dev. in kcal/mol:  ',sqrt(ens%avgR(2)/(1.0*ens%avgcnt) - &
 &                                                          ((ens%avgU+ens%avgK)/(1.0*ens%avgcnt))**2)
      if (fycxyz.eq.1) then
        write(ilog,57) 'CTE Stand. Dev. in kcal/mol :  ',sqrt(ens%avgR(3)/(1.0*ens%avgcnt) - &
 &                                                            ((ens%avgU+ens%avgK2)/(1.0*ens%avgcnt))**2)
      end if
      if (ens%flag.eq.2) then
        write(ilog,57) 'T.E.-Drift in kcal/mol/ps   :  ',(ens%insU+ens%insK-ens%insR(2))/(dyn_dt*ens%avgcnt)
        if (fycxyz.eq.1) then
          write(ilog,57) 'CTE-Drift in kcal/mol/ps    :  ',(ens%insU+ens%insK2-ens%insR(3))/(dyn_dt*ens%avgcnt)
        end if
      end if
      write(ilog,57) 'GLCF in gA^2/mol/ps^4       :  ',ens%avgR(1)/(1.0*ens%avgcnt)
      if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
        write(ilog,57) 'Pressure in bar             :  ',ens%avgP/(1.0*ens%avgcnt)
        write(ilog,57) 'Volume in nm^3              :  ',ens%avgV/(1.0*ens%avgcnt)
      end if
      if ((cart_cons_mode.gt.1).AND.(no_shake.EQV..false.).AND.&
 &(cart_cons_grps.gt.(settle_tip3ps+settle_tip4ps+settle_tip5ps+settle_tip4pes+settle_spcs+settle_rest))) then
        write(ilog,*)
        write(ilog,*) 'Statistics of constraint solver (averaged over all appl. constraints groups):'
        if (cart_cons_method.eq.4) then
          write(ilog,57) 'Average order of LINCS matrix inverse expansion :   ',&
 &sum((1.0*constraints(nonsettlelst(1:shake_cnt))%iters)/(1.0*constraints(nonsettlelst(1:shake_cnt))%itercnt))/&
 &(1.0*shake_cnt)
        else
          write(ilog,57) 'Average number of steps for iterative solver    :  ',&
 &sum((1.0*constraints(nonsettlelst(1:shake_cnt))%iters)/(1.0*constraints(nonsettlelst(1:shake_cnt))%itercnt))/&
 &(1.0*shake_cnt)
        end if
      end if
#ifdef ENABLE_MPI
      if (use_REMC.EQV..true.) then
        write(ilog,61) 'Net number of acc. RE-Swaps: ',acc%nre
      end if
#endif
      if (fycxyz.eq.1) then
       write(ilog,*)
       write(ilog,*) 'Molecule-resolved average temperatures for RB d.o.f.:'
       write(ilog,*) 'Mol. |Typ.|       | Compon. 1  | Compon. 2  | Compon. 3  |'
       do imol=1,nmol
        write(ilog,59) imol,moltypid(imol),' TRANS.: ',(dc_di(imol)%avgT(i)/(1.0*ens%avgcnt),i=1,3)
        if (atmol(imol,2).gt.(atmol(imol,1)+1)) then
          write(ilog,59) imol,moltypid(imol),' ROT.:   ',(dc_di(imol)%avgT(i)/(1.0*ens%avgcnt),i=4,6)
        end if
       end do
       ttc = 0 
       do imol=1,nmol
        if (othidxmol(moltypid(imol)).gt.0) then
         if (ttc.eq.0) then
           write(ilog,*)
           write(ilog,*) 'Average temperatures for OTHER internal d.o.f. for entire molecule:'
           write(ilog,*) 'Molecule # | Various angles ...'
           ttc = 1
         end if
         write(ilog,62) imol,' :',(dc_di(imol)%avgT(i)/(1.0*ens%avgcnt),&
 &        i=othidxmol(moltypid(imol)),ntormol(moltypid(imol))+6)
        end if
       end do
       if (((ens%flag.eq.3).OR.(ens%flag.eq.4)).AND.&
 &        ((bnd_type.ge.3).AND.(bnd_type.le.4))) then
        write(ilog,60) 'BOUNDARY PARTICLE :',bnd_avgT/(1.0*ens%avgcnt)
       end if
       if (ntorsn.gt.0) then
        write(ilog,*)
        write(ilog,*) 'Residue-resolved average temperatures for internal d.o.f.:'
       end if
       do imol=1,nmol
        do rs=rsmol(imol,1),rsmol(imol,2)
          resname = amino(seqtyp(rs))
          if (wnr(rs).gt.0) then
            write(ilog,58) resname,'->OMEGA:',dc_di(imol)%avgT(wnr(rs))/(1.0*ens%avgcnt)
          end if
          if (fnr(rs).gt.0) then
            write(ilog,58) resname,'->PHI  :',dc_di(imol)%avgT(fnr(rs))/(1.0*ens%avgcnt)
          end if
          if (ynr(rs).gt.0) then
            write(ilog,58) resname,'->PSI  :',dc_di(imol)%avgT(ynr(rs))/(1.0*ens%avgcnt)
          end if
          ttc = 0
          do i=1,nnucs(rs)
            if (nucsnr(i,rs).gt.0) then
              ttc = ttc + 1
              prtout(ttc) = dc_di(imol)%avgT(nucsnr(i,rs))/(1.0*ens%avgcnt)
            end if
          end do
          if (ttc.gt.0) write(ilog,58) resname,'->NUC.S:',prtout(1:ttc)
          ttc = 0
          do i=1,nchi(rs)
            if (chinr(i,rs).gt.0) then
              ttc = ttc + 1
              prtout(ttc) = dc_di(imol)%avgT(chinr(i,rs))/(1.0*ens%avgcnt)
            end if
          end do
          if (ttc.gt.0) write(ilog,58) resname,'->CHIS :',prtout(1:ttc)
        end do
       end do
      end if
    end if
  end if
!
  return
!
  end 
!
!--------------------------------------------------------------------------
!

