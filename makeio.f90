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
! MAIN AUTHOR:   Rohit Pappu                                               !
! CONTRIBUTIONS: Andreas Vitalis, Hoang Tran                               !
!                                                                          !
!--------------------------------------------------------------------------!
!
!
! #############################################################
! ##                                                         ##
! ## subroutine makeio -- open or close desired output files ##
! ##                                                         ##
! #############################################################
!
! "makeio" opens or closes output files for Monte Carlo simulations
! mode = 1: open
! mode = 2: or anything else close
!
! for parallel computing filehandles are a mess, here we open all files in the
! same (running) directory (REMC) or we suppress output (averaging)
! WARNING:
! this means that ALL GENERAL CHANGES HAVE TO BE MADE AT LEAST TWICE 
!
subroutine makeio(mode)
!
  use iounit
  use mcsums
  use polyavg
  use system
  use molecule
  use torsn
  use movesets
  use mpistuff
  use energies
  use dssps
  use fos
  use wl
  use paircorr
  use clusters

!
  implicit none
!
  integer freeunit,lext,mode,mt,imol
  logical exists,doperser
  character(3) mtstr
  character(60) dihedfile,polmrfile
  character(60) holesfile,accfile,trcvfile
  character(60) persfile,sasafile,enfile,dsspfile,sav_det_file,val_det_file,sav_det_file_map
  character(60) pka_test_file_fos,pka_test_file_fos_weights,pka_test_file_bond_length ! MArtin , I should move that for a declaration only in MPI version
  
  character(60) HSQ_tries_file,HSQ_acc_file,HSQ_acc_ratio,HSQ_acc_summary
  
  character(60) pka_test_file_atm_vol,pka_test_file_atm_id,pka_test_file_sav,pka_test_file_charge
  character(60) pka_test_file_LJ_rad,pka_test_file_LJ_eps,pka_test_file_LJ_sig
  character(60) pka_test_file_atr,pka_test_file_biotype
  character(60) pka_charge2_file, pka_bond_length, pka_bond_angle
  character(60) pka_test_file_fos_id
  character(60) pka_atbvol,keep_track_q ! Martin hsq
  character(60) phtfile,ensfile,elsfile,pcfile,rmsdfile
  logical fycopen,polymopen,ensopen,dsspopen,&
 &accopen,enopen,persopen,holesopen,sasaopen,phtopen,elsopen,pcopen,rmsdopen,&
 &sav_detailopen,val_det_file_open,state_q_open
 logical test_file1open,test_file2open,test_file3open,test_file4open,test_file5open
 logical test_file6open,test_file7open,test_file8open,test_file9open,test_file10open
 logical test_file11open, test_file12open,test_file13open,test_file14open,test_file15open,test_file16open,test_file17open
 logical test_file18open, test_file19open,test_file20open,test_file21open
 logical inst_der_ham_out_open,sav_detail_map_open

 character(60) pka_test_file_limits


 
 
#ifdef ENABLE_MPI
  integer masterrank,t1,t2
  character(3) xpont
  character(60) remcfile,retrfile
  logical remcopen,retropen
#endif
!
! the save/data construct ensures that the logicals are treated like subroutine-globals
! including initialization (once(!)) 
! Martin : the save statement seems actually quite unuseful, since the variable are never used in other part of the program
  data fycopen/.false./,&
 &polymopen/.false./,&
 &accopen/.false./,dsspopen/.false./,&
 &enopen/.false./,persopen/.false./,ensopen/.false./,&
 &holesopen/.false./,sasaopen/.false./,phtopen/.false./,elsopen/.false./,pcopen/.false./,rmsdopen/.false./,&
 &sav_detailopen/.false./test_file1open/.false./,test_file2open/.false./,test_file3open/.false./,test_file4open/.false./&
 &,test_file5open/.false./,test_file6open/.false./,test_file7open/.false./,test_file8open/.false./,test_file9open/.false./,&
 &test_file10open/.false./,test_file11open/.false./, test_file12open/.false./,test_file13open/.false./,test_file14open/.false./,&
 &test_file15open/.false./,test_file16open/.false./,val_det_file_open/.false./,state_q_open/.false./, test_file17open/.false./,&
 &sav_detail_map_open/.false./,test_file18open/.false./,test_file19open/.false./,test_file20open/.false./,test_file21open/.false./

#ifdef ENABLE_MPI
  data remcopen/.false./,retropen/.false./
#endif
!
  save fycopen,polymopen,ensopen,&
 &accopen,enopen,persopen,holesopen,sasaopen,&
 &phtopen,dsspopen,elsopen,pcopen,rmsdopen,&
 &sav_detailopen
#ifdef ENABLE_MPI
  save remcopen,retropen,test_file3open
#endif
!
  doperser = .false.
  do imol=1,nmol
    if (do_pers(moltypid(imol)).EQV..true.) then
      doperser = .true.
      exit
    end if
  end do
 
#ifdef ENABLE_MPI
!
  masterrank = 0
  if (use_REMC.EQV..true.) then
!   in this case we have to open the output files for each node individually
    if (mode .eq. 1) then
      lext = 3
      call int2str(myrank,xpont,lext)
      
      dihedfile  = 'N_'//xpont(1:lext)//'_FYC.dat'
      polmrfile  = 'N_'//xpont(1:lext)//'_POLYMER.dat'
      holesfile  = 'N_'//xpont(1:lext)//'_HOLES.dat'
      enfile  = 'N_'//xpont(1:lext)//'_ENERGY.dat'
      if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
        accfile = 'N_'//xpont(1:lext)//'_ACCEPTANCE.dat'
      end if
      if (((dyn_mode.ge.2).AND.(dyn_mode.le.5)).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
        ensfile = 'N_'//xpont(1:lext)//'_ENSEMBLE.dat'
      end if
      persfile = 'N_'//xpont(1:lext)//'_PERSISTENCE.dat'
      sasafile = 'N_'//xpont(1:lext)//'_SAV.dat'
      dsspfile = 'N_'//xpont(1:lext)//'_DSSP_RUNNING.dat'
      phtfile = 'N_'//xpont(1:lext)//'_PHTIT.dat'
      remcfile = 'N_'//xpont(1:lext)//'_EVEC.dat'
      retrfile = 'N_'//xpont(1:lext)//'_REXTRACE.dat'
      elsfile = 'N_'//xpont(1:lext)//'_ELS_WFRAMES.dat'
      pcfile = 'N_'//xpont(1:lext)//'_GENERAL_DIS.dat'
      rmsdfile = 'N_'//xpont(1:lext)//'_RMSD.dat'
      if (print_det.eqv..true.) then 
        sav_det_file='T_'//xpont(1:lext)//'_SAV_detail.dat'!martin
        sav_det_file_map='T_'//xpont(1:lext)//'_SAV_detail_mapping.dat'
      end if 
      if (do_pka.eqv..true.) then 
        if (debug_pka.eqv..true.) then 

            pka_test_file_charge='SN_'//xpont(1:lext)//'_charge.dat'!martin 
            pka_test_file_fos='SN_'//xpont(1:lext)//'_FOS.dat'!martin 
            
            pka_test_file_fos_id='SN_'//xpont(1:lext)//'_FOS_atoms_id.dat'!martin 
            
            pka_test_file_fos_weights='SN_'//xpont(1:lext)//'_FOS_weights.dat'!martin 
            pka_test_file_bond_length='SN_'//xpont(1:lext)//'_bond_length.dat'!martin 
            pka_test_file_atm_vol='SN_'//xpont(1:lext)//'_atom_volume.dat'!martin 
            pka_test_file_atm_id='SN_'//xpont(1:lext)//'_atom_id.dat'!martin 
            pka_test_file_sav='SN_'//xpont(1:lext)//'_SAV.dat'!martin 
            pka_test_file_LJ_rad='SN_'//xpont(1:lext)//'_LJ_radius.dat'!martin
            pka_test_file_LJ_eps='SN_'//xpont(1:lext)//'_LJ_eps.dat'!martin
            pka_test_file_LJ_sig='SN_'//xpont(1:lext)//'_LJ_sig.dat'!martin !
            pka_test_file_atr='SN_'//xpont(1:lext)//'_atr.dat'!martin !
            pka_test_file_biotype='SN_'//xpont(1:lext)//'_biotype.dat'!martin !
            pka_test_file_limits='X_N_'//xpont(1:lext)//'_limits.dat' ! MArtin 
            pka_bond_angle='SN_'//xpont(1:lext)//'_bond_angle.dat' ! MArtin 
            pka_atbvol='SN_'//xpont(1:lext)//'_atbvol'!martin 
            val_det_file='SN_'//xpont(1:lext)//'_values'!martin 
        else
            val_det_file='SN_'//xpont(1:lext)//'_values'!martin 
        end if 
      end if 

      if (print_det.eqv..true.) then 

          inquire(file=sav_det_file,exist=exists)
          if (exists.EQV..true.) then
            isav_det = freeunit()
            open (unit=isav_det,file=sav_det_file,status='old')
            close(unit=isav_det,status='delete')
          end if
          if ((sav_det_freq.gt.0.0).AND.(sav_det_freq.le.nsim)) then
            isav_det = freeunit()
            open (unit=isav_det,file=sav_det_file,status='new') 
            sav_detailopen = .true.
          end if 
          
          inquire(file=sav_det_file_map,exist=exists)
          if (exists.EQV..true.) then
            isav_det_map = freeunit()
            open (unit=isav_det_map,file=sav_det_file_map,status='old')
            close(unit=isav_det_map,status='delete')
          end if
          if ((sav_det_freq.gt.0.0).AND.(sav_det_freq.le.nsim)) then
            isav_det_map = freeunit()
            open (unit=isav_det_map,file=sav_det_file_map,status='new') 
            sav_detail_map_open = .true.
          end if 
          
      end if 
      if (do_pka.eqv..true.) then 
          if (debug_pka.eqv..true.) then 
            inquire(file=pka_test_file_limits,exist=exists)
            if (exists.EQV..true.) then
                file_limits  = freeunit()
                open (unit=file_limits ,file=pka_test_file_limits,status='old')
                close(unit=file_limits ,status='delete')
            end if
            if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                file_limits = freeunit()
                open (unit=file_limits,file=pka_test_file_limits,status='new') 
                test_file13open = .true.
            end if 
            
              inquire(file=pka_test_file_charge,exist=exists)
              if (exists.EQV..true.) then
                test_file1 = freeunit()
                open (unit=test_file1,file=pka_test_file_charge,status='old')
                close(unit=test_file1,status='delete')
              end if
              if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                test_file1 = freeunit()
                open (unit=test_file1,file=pka_test_file_charge,status='new') 
                test_file1open = .true.
              end if 
              
              inquire(file=pka_test_file_fos,exist=exists)
              if (exists.EQV..true.) then
                test_file2 = freeunit()
                open (unit=test_file2,file=pka_test_file_fos,status='old')
                close(unit=test_file2,status='delete')
              end if
              if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                test_file2 = freeunit()
                open (unit=test_file2,file=pka_test_file_fos,status='new') 
                test_file2open = .true.
              end if 

              inquire(file=pka_test_file_fos_weights,exist=exists)
              if (exists.EQV..true.) then
                test_file3 = freeunit()
                open (unit=test_file3,file=pka_test_file_fos_weights,status='old')
                close(unit=test_file3,status='delete')
              end if
              if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                test_file3 = freeunit()
                open (unit=test_file3,file=pka_test_file_fos_weights,status='new') 
                test_file3open = .true.
              end if  

              inquire(file=pka_test_file_bond_length,exist=exists)
              if (exists.EQV..true.) then
                test_file4 = freeunit()
                open (unit=test_file4,file=pka_test_file_bond_length,status='old')
                close(unit=test_file4,status='delete')
              end if
              if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                test_file4 = freeunit()
                open (unit=test_file4,file=pka_test_file_bond_length,status='new') 
                test_file4open = .true.
              end if  

              inquire(file=pka_test_file_atm_vol,exist=exists)
              if (exists.EQV..true.) then
                test_file5 = freeunit()
                open (unit=test_file5,file=pka_test_file_atm_vol,status='old')
                close(unit=test_file5,status='delete')
              end if
              if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                test_file5 = freeunit()
                open (unit=test_file5,file=pka_test_file_atm_vol,status='new') 
                test_file5open = .true.
              end if  

              inquire(file=pka_test_file_atm_id,exist=exists)
              if (exists.EQV..true.) then
                test_file6 = freeunit()
                open (unit=test_file6,file=pka_test_file_atm_id,status='old')
                close(unit=test_file6,status='delete')
              end if
              if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                test_file6 = freeunit()
                open (unit=test_file6,file=pka_test_file_atm_id,status='new') 
                test_file6open = .true.
              end if  

              inquire(file=pka_test_file_sav,exist=exists)
              if (exists.EQV..true.) then
                test_file7 = freeunit()
                open (unit=test_file7,file=pka_test_file_sav,status='old')
                close(unit=test_file7,status='delete')
              end if
              if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                test_file7 = freeunit()
                open (unit=test_file7,file=pka_test_file_sav,status='new') 
                test_file7open = .true.
              end if  

              inquire(file=pka_test_file_LJ_rad,exist=exists)
              if (exists.EQV..true.) then
                test_file8 = freeunit()
                open (unit=test_file8,file=pka_test_file_LJ_rad,status='old')
                close(unit=test_file8,status='delete')
              end if
              if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                test_file8 = freeunit()
                open (unit=test_file8,file=pka_test_file_LJ_rad,status='new') 
                test_file8open = .true.
              end if

              inquire(file=pka_test_file_LJ_eps,exist=exists)
              if (exists.EQV..true.) then
                test_file9 = freeunit()
                open (unit=test_file9,file=pka_test_file_LJ_eps,status='old')
                close(unit=test_file9,status='delete')
              end if
              if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                test_file9 = freeunit()
                open (unit=test_file9,file=pka_test_file_LJ_eps,status='new') 
                test_file9open = .true.
              end if  

              inquire(file=pka_test_file_LJ_sig,exist=exists)
              if (exists.EQV..true.) then
                test_file10 = freeunit()
                open (unit=test_file10,file=pka_test_file_LJ_sig,status='old')
                close(unit=test_file10,status='delete')
              end if
              if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                test_file10 = freeunit()
                open (unit=test_file10,file=pka_test_file_LJ_sig,status='new') 
                test_file11open = .true.
              end if  

              inquire(file=pka_test_file_atr,exist=exists)
              if (exists.EQV..true.) then
                test_file11 = freeunit()
                open (unit=test_file11,file=pka_test_file_atr,status='old')
                close(unit=test_file11,status='delete')
              end if
              if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                test_file11 = freeunit()
                open (unit=test_file11,file=pka_test_file_atr,status='new') 
                test_file11open = .true.
              end if  

              inquire(file=pka_test_file_biotype,exist=exists)
              if (exists.EQV..true.) then
                test_file12 = freeunit()
                open (unit=test_file12,file=pka_test_file_biotype,status='old')
                close(unit=test_file12,status='delete')
              end if
              if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                test_file12 = freeunit()
                open (unit=test_file12,file=pka_test_file_biotype,status='new') 
                test_file12open = .true.
              end if

            inquire(file=pka_bond_angle,exist=exists)
            if (exists.EQV..true.) then
                test_file15  = freeunit()
                open (unit=test_file15 ,file=pka_bond_angle,status='old')
                close(unit=test_file15 ,status='delete')
            end if
            if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                test_file15 = freeunit()
                open (unit=test_file15,file=pka_bond_angle,status='new') 
                test_file15open = .true.
            end if  
                            
              

            inquire(file=pka_atbvol,exist=exists)
            if (exists.EQV..true.) then
                test_file14  = freeunit()
                open (unit=test_file14 ,file=pka_atbvol,status='old')
                close(unit=test_file14 ,status='delete')
            end if
            if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                test_file14 = freeunit()
                open (unit=test_file14,file=pka_atbvol,status='new') 
                test_file14open = .true.                    
            end if      


              
            inquire(file=pka_test_file_fos_id,exist=exists)
            if (exists.EQV..true.) then
                test_file17  = freeunit()
                open (unit=test_file17 ,file=pka_test_file_fos_id,status='old')
                close(unit=test_file17 ,status='delete')
            end if
            if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                test_file17 = freeunit()
                open (unit=test_file17,file=pka_test_file_fos_id,status='new') 
                test_file17open = .true.                    
            end if   
            
            inquire(file=val_det_file,exist=exists)
              if (exists.EQV..true.) then
                val_det = freeunit()
                open (unit=val_det,file=val_det_file,status='old')
                close(unit=val_det,status='delete')
              end if

              val_det = freeunit()
              open (unit=val_det,file=val_det_file,status='new') 
              val_det_file_open = .true.

          else 
              inquire(file=val_det_file,exist=exists)
              if (exists.EQV..true.) then
                val_det = freeunit()
                open (unit=val_det,file=val_det_file,status='old')
                close(unit=val_det,status='delete')
              end if

              val_det = freeunit()
              open (unit=val_det,file=val_det_file,status='new') 
              val_det_file_open = .true.

          end if  
      end if 
      inquire(file=dihedfile,exist=exists)
      if (exists.EQV..true.) then
        idihed = freeunit()
        open (unit=idihed,file=dihedfile,status='old')
        close(unit=idihed,status='delete')
      end if
      if (torout.le.nsim) then
        idihed = freeunit()
        open (unit=idihed,file=dihedfile,status='new')
        fycopen = .true.
      end if
! setup filehandle for polymeric outputs
      inquire(file=polmrfile,exist=exists)
      if (exists.EQV..true.) then
        ipolmr = freeunit()
        open (unit=ipolmr,file=polmrfile,status='old')
        close(unit=ipolmr,status='delete')
      end if
      if (polout .le. nsim) then
        ipolmr = freeunit()
        open (unit=ipolmr,file=polmrfile,status='new') 
        polymopen = .true.
      end if 
!   
! setup filehandle for holes outputs
      inquire(file=holesfile,exist=exists)
      if (exists.EQV..true.) then
        iholes = freeunit()
        open (unit=iholes,file=holesfile,status='old')
        close(unit=iholes,status='delete')
      end if
      if (holescalc.le.nsim) then
        iholes = freeunit()
        open (unit=iholes,file=holesfile,status='new') 
        holesopen = .true.
      end if
! 
! setup filehandle for pH titration outputs
      inquire(file=phtfile,exist=exists)
      if (exists.EQV..true.) then
        ipht = freeunit()
        open (unit=ipht,file=phtfile,status='old')
        close(unit=ipht,status='delete')
      end if
      if ((phfreq.gt.0.0).AND.(phout.le.nsim)) then
        ipht = freeunit()
        open (unit=ipht,file=phtfile,status='new') 
        phtopen = .true.
      end if 
      

!   
! setup filehandle for energy outputs
      inquire(file=enfile,exist=exists)
      if (exists.EQV..true.) then
        iene = freeunit()
        open (unit=iene,file=enfile,status='old')
        close(unit=iene,status='delete')
      end if
      if (enout.le.nsim) then
        iene = freeunit()
        open (unit=iene,file=enfile,status='new')     
        enopen = .true.
      end if
!
! setup filehandle for acceptance outputs
      if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
        inquire(file=accfile,exist=exists)
        if (exists.EQV..true.) then
          iacc = freeunit()
          open (unit=iacc,file=accfile,status='old')
          close(unit=iacc,status='delete')
        end if
        if (accout.le.nsim) then
          iacc = freeunit()
          open (unit=iacc,file=accfile,status='new')
          accopen = .true.
          end if
      end if
!
! setup filehandle for ensemble outputs
      if ((pdb_analyze.EQV..false.).AND.(((dyn_mode.ge.2).AND.(dyn_mode.le.5)).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
        inquire(file=ensfile,exist=exists)
        if (exists.EQV..true.) then
          iens = freeunit()
          open (unit=iens,file=ensfile,status='old')
          close(unit=iens,status='delete')
        end if
        if (ensout.le.nsim) then
          iens = freeunit()
          open (unit=iens,file=ensfile,status='new')
          ensopen = .true.
        end if
      end if
!
! setup filehandle for persistence length calculation output file
      inquire(file=persfile,exist=exists)
      if (exists.EQV..true.) then
        ipers = freeunit()
        open (unit=ipers,file=persfile,status='old')
        close(unit=ipers,status='delete')
      end if
      if ((polcalc.le.nsim).AND.(doperser.EQV..true.)) then
        ipers = freeunit()
        open (unit=ipers,file=persfile,status='new')
        persopen = .true.
      end if
!
! setup filehandle for on-the-fly SAV output
      inquire(file=sasafile,exist=exists)
      if (exists.EQV..true.) then
        isasa = freeunit()
        open (unit=isasa,file=sasafile,status='old')
        close(unit=isasa,status='delete')
      end if
      if (savreq%instfreq.le.nsim) then
        isasa = freeunit()
        open (unit=isasa,file=sasafile,status='new')
        sasaopen = .true.
      end if
!
! setup filehandle for on-the-fly DSSP output
      inquire(file=dsspfile,exist=exists)
      if (exists.EQV..true.) then
        idssp = freeunit()
        open (unit=idssp,file=dsspfile,status='old')
        close(unit=idssp,status='delete')
      end if
      if ((dsspcalc.le.nsim).AND.(inst_dssp.EQV..true.)) then
        idssp = freeunit()
        open (unit=idssp,file=dsspfile,status='new')
        dsspopen = .true.
      end if
!
! setup filehandle for on-the-fly distance output
      inquire(file=pcfile,exist=exists)
      if (exists.EQV..true.) then
        igpc = freeunit()
        open (unit=igpc,file=pcfile,status='old')
        close(unit=igpc,status='delete')
      end if
      if ((pccalc.le.nsim).AND.(inst_gpc.gt.0).AND.(inst_gpc.le.(nsim/pccalc)).AND.(gpc%nos.gt.0)) then
        igpc = freeunit()
        open (unit=igpc,file=pcfile,status='new')
        pcopen = .true.
      end if
!
! setup filehandle for on-the-fly RMSD output
      inquire(file=rmsdfile,exist=exists)
      if (exists.EQV..true.) then
        irmsd = freeunit()
        open (unit=irmsd,file=rmsdfile,status='old')
        close(unit=irmsd,status='delete')
      end if
      if ((align%calc.le.nsim).AND.(align%yes.EQV..true.).AND.(align%instrmsd.EQV..true.)) then
        irmsd = freeunit()
        open (unit=irmsd,file=rmsdfile,status='new')
        rmsdopen = .true.
      end if
!
! setup filehandle for on-the-fly cross-node energy output
      inquire(file=remcfile,exist=exists)
      if (exists.EQV..true.) then
        iremc = freeunit()
        open (unit=iremc,file=remcfile,status='old')
        close(unit=iremc,status='delete')
      end if
      if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
        iremc = freeunit()
        open (unit=iremc,file=remcfile,status='new')
        remcopen = .true.
      end if
!
! setup filehandle for on-the-fly RE trace (original structure vs. condition)
      if ((myrank.eq.masterrank).AND.(pdb_analyze.EQV..false.)) then
        inquire(file=retrfile,exist=exists)
        if (exists.EQV..true.) then
          iretr = freeunit()
          open(unit=iretr,file=retrfile,status='old',position='append')
          if (do_restart.EQV..false.) close(unit=iretr,status='delete')
        end if
        inquire(file=retrfile,exist=exists)
        if ((re_freq.le.nsim).AND.(inst_retr.EQV..true.).AND.((do_restart.EQV..false.).OR.(exists.EQV..false.))) then
          iretr = freeunit()
          open (unit=iretr,file=retrfile,status='new')
          retropen = .true.
        else if ((re_freq.le.nsim).AND.(inst_retr.EQV..true.).AND.(do_restart.EQV..true.).AND.(exists.EQV..true.)) then
          retropen = .true.
          write(ilog,*) 'Warning. Trace file for replica exchange run exists and will be appended (this may be inappropriate).'
          close(unit=iretr)
          open(unit=iretr,file=retrfile,status='old',action='read')
          do while (.true.)
            read(iretr,*,iostat=imol) mt,mpi_lmap(1:mpi_nodes,1)
            if (imol.ne.0) exit
            mpi_lmap(1:mpi_nodes,2) = mpi_lmap(1:mpi_nodes,1)
          end do
          close(unit=iretr)
          open(unit=iretr,file=retrfile,status='old',position='append')
        else if ((do_restart.EQV..true.).AND.(exists.EQV..true.)) then
          close(unit=iretr)
        end if
      else if ((myrank.eq.masterrank).AND.(pdb_analyze.EQV..true.).AND.(re_aux(3).eq.1)) then
        call strlims(re_traceinfile,t1,t2)
        inquire(file=re_traceinfile(t1:t2),exist=exists)
        if (exists.EQV..true.) then
          iretr = freeunit()
          open(unit=iretr,file=re_traceinfile(t1:t2),status='old',action='read')
        else
          write(ilog,*) 'Warning. In MPI trajectory analysis mode, requested un"scrambling" of replica exchange trajectories &
 &cannot be performed because trace file cannot be opened (got: ',re_traceinfile(t1:t2),').'
          re_aux(3) = 0
        end if
      else if ((pdb_analyze.EQV..true.).AND.(re_aux(3).eq.1)) then
        call strlims(re_traceinfile,t1,t2)
        inquire(file=re_traceinfile(t1:t2),exist=exists)
        if (exists.EQV..false.) then
          write(ilog,*) 'Warning. In MPI trajectory analysis mode, requested un"scrambling" of replica exchange trajectories &
 &cannot be performed because trace file cannot be opened (got: ',re_traceinfile(t1:t2),').'
          re_aux(3) = 0
        end if
      end if
!
! setup filehandle for temporary TRCV (covariance of torsions) file
      if (covcalc.le.nsim) then
        do mt=1,nangrps
          call int2str(mt,mtstr,lext)
          trcvfile = 'N_'//xpont(1:lext)//'_TRCV_'//mtstr(1:lext)&
 &                                                //'.tmp'
          inquire(file=trcvfile,exist=exists)
          if (exists.EQV..true.) then
            itrcv(mt) = freeunit()
            open (unit=itrcv(mt),file=trcvfile,status='old')
            close(unit=itrcv(mt),status='delete')
          end if
!
          if (do_pol(moltypid(molangr(mt,1))).EQV..false.) cycle
          itrcv(mt) = freeunit()
          open (unit=itrcv(mt),file=trcvfile,status='new')
          trcvopen(mt) = .true.
        end do
      end if
!
! setup filehandle for frame weights from ELS (appended for restarts)
      inquire(file=elsfile,exist=exists)
      if (exists.EQV..true.) then
        hmjam%iwtsfile = freeunit()
        open(unit=hmjam%iwtsfile,file=elsfile,status='old',position='append')
        if ((do_restart.EQV..false.).AND.(hmjam%prtfrmwts.le.nsim)) then
          close(unit=hmjam%iwtsfile,status='delete')
        else
          close(unit=hmjam%iwtsfile)
        end if
      end if
      inquire(file=elsfile,exist=exists)
      if ((hmjam%prtfrmwts.le.nsim).AND.(exists.EQV..false.)) then
        hmjam%iwtsfile = freeunit()
        open (unit=hmjam%iwtsfile,file=elsfile,status='new')
        elsopen = .true.
      else if (hmjam%prtfrmwts.le.nsim) then
        hmjam%iwtsfile = freeunit()
        open(unit=hmjam%iwtsfile,file=elsfile,status='old',position='append')
        elsopen = .true.
      end if
!
!   close files and free up filehandles
    else
      if (fycopen) close(unit=idihed)
      if (polymopen) close(unit=ipolmr)
      if (enopen) close(unit=iene)
      if (accopen) close(unit=iacc)
      if (ensopen) close(unit=iens)       
      if (persopen) close(unit=ipers)
      if (sasaopen) close(unit=isasa)
      if (dsspopen) close(unit=idssp)
      if (phtopen) close(unit=ipht)
      if (sav_detailopen) close(unit=isav_det)
      if (test_file1open) close(unit=test_file1) ! Martin : this part needs to be checked, but isn't crucial
      if (test_file2open) close(unit=test_file2)
      if (test_file3open) close(unit=test_file3)
      if (test_file4open) close(unit=test_file4)
      if (test_file5open) close(unit=test_file5)
      if (test_file6open) close(unit=test_file6)
      if (test_file7open) close(unit=test_file7)
      if (test_file8open) close(unit=test_file8)
      if (test_file9open) close(unit=test_file9)
      if (test_file10open) close(unit=test_file10)
      if (test_file11open) close(unit=test_file11)
      if (test_file12open) close(unit=test_file12)
      if (test_file13open) close(unit=file_limits)
      
      

      if (remcopen) close(unit=iremc)
      if (retropen) close(unit=iretr)
      if (covcalc.le.nsim) then
        do mt=1,nangrps
          if (trcvopen(mt)) close(unit=itrcv(mt))
        end do
      end if
      if (elsopen) close(unit=hmjam%iwtsfile)
      if (pcopen) close(unit=igpc)
      if (rmsdopen) close(unit=irmsd)
    end if
!
  else if (use_MPIAVG.EQV..true.) then
!
    if (mode .eq. 1) then
      lext = 3
      call int2str(myrank,xpont,lext)
      dihedfile  = 'FYC.dat'
      polmrfile  = 'POLYMER.dat'
      enfile  = 'ENERGY.dat'
      if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
        accfile  = 'ACCEPTANCE.dat'
      end if
      if (((dyn_mode.ge.2).AND.(dyn_mode.le.5)).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
        ensfile = 'ENSEMBLE.dat'
      end if
      persfile = 'PERSISTENCE.dat'
      sasafile = 'SAV.dat'
      phtfile = 'PHTIT.dat'
      retrfile = 'N_'//xpont(1:lext)//'_PIGSTRACE.dat'
      elsfile = 'N_'//xpont(1:lext)//'_ELS_WFRAMES.dat'
      pcfile = 'N_'//xpont(1:lext)//'_GENERAL_DIS.dat'
      
      
      if (print_det.eqv..true.) then 
        sav_det_file = 'SAV_detail.dat' !martin
        sav_det_file_map='SAV_detail_mapping.dat'
      end if 
        
        if (myrank.eq.masterrank) then
!   setup filehandle for dihedral angle output
        inquire(file=dihedfile,exist=exists)
        if (exists.EQV..true.) then
          idihed = freeunit()
          open (unit=idihed,file=dihedfile,status='old')
          close(unit=idihed,status='delete')
        end if
        if (torout .le. nsim) then
          idihed = freeunit()
          open (unit=idihed,file=dihedfile,status='new')
          fycopen = .true.
        end if
!
!   setup filehandle for polymeric outputs
        inquire(file=polmrfile,exist=exists)
        if (exists.EQV..true.) then
          ipolmr = freeunit()
          open (unit=ipolmr,file=polmrfile,status='old')
          close(unit=ipolmr,status='delete')
        end if
        if (polout .le. nsim) then
          ipolmr = freeunit()
          open (unit=ipolmr,file=polmrfile,status='new')
          polymopen = .true.
        end if
!
!   setup filehandle for energy outputs
        inquire(file=enfile,exist=exists)
        if (exists.EQV..true.) then
          iene = freeunit()
          open (unit=iene,file=enfile,status='old')
          close(unit=iene,status='delete')
        end if
        if (enout.le.nsim) then
          iene = freeunit()
          open (unit=iene,file=enfile,status='new')
          enopen = .true.
        end if
!
!   setup filehandle for acceptance outputs
        if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
          inquire(file=accfile,exist=exists)
          if (exists.EQV..true.) then
            iacc = freeunit()
            open (unit=iacc,file=accfile,status='old')
            close(unit=iacc,status='delete')
          end if
          if (accout.le.nsim) then
            iacc = freeunit()
            open (unit=iacc,file=accfile,status='new')
            accopen = .true.
          end if
        end if
!
!   setup filehandle for ensemble outputs
        if ((pdb_analyze.EQV..false.).AND.(((dyn_mode.ge.2).AND.(dyn_mode.le.5)).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
          inquire(file=ensfile,exist=exists)
          if (exists.EQV..true.) then
            iens = freeunit()
            open (unit=iens,file=ensfile,status='old')
            close(unit=iens,status='delete')
          end if
          if (ensout.le.nsim) then
            iens = freeunit()
            open (unit=iens,file=ensfile,status='new')
            ensopen = .true.
          end if
        end if
!
!   setup filehandle for persistence length calculation output file
        inquire(file=persfile,exist=exists)
        if (exists.EQV..true.) then
          ipers = freeunit()
          open (unit=ipers,file=persfile,status='old')
          close(unit=ipers,status='delete')
        end if
        if ((polcalc.le.nsim).AND.(doperser.EQV..true.)) then
         ipers = freeunit()
         open (unit=ipers,file=persfile,status='new')
         persopen = .true.
        end if
!
!   setup filehandle for pH titration outputs
        inquire(file=phtfile,exist=exists)
        if (exists.EQV..true.) then
           ipht = freeunit()
           open (unit=ipht,file=phtfile,status='old')
           close(unit=ipht,status='delete')
        end if
        if ((phfreq.gt.0.0d0).AND.(phout.le.nsim)) then
           ipht = freeunit()
           open (unit=ipht,file=phtfile,status='new') 
           phtopen = .true.
        end if 
!
!   setup filehandle for on-the-fly SAV output
        inquire(file=sasafile,exist=exists)
        if (exists.EQV..true.) then
          isasa = freeunit()
          open (unit=isasa,file=sasafile,status='old')
          close(unit=isasa,status='delete')
        end if
        if (savreq%instfreq.le.nsim) then
          isasa = freeunit()
          open (unit=isasa,file=sasafile,status='new')
          sasaopen = .true.
        end if
!
      end if ! if master
!
! setup filehandle for on-the-fly distance output
      inquire(file=pcfile,exist=exists)
      if (exists.EQV..true.) then
        igpc = freeunit()
        open (unit=igpc,file=pcfile,status='old')
        close(unit=igpc,status='delete')
      end if
      if ((pccalc.le.nsim).AND.(inst_gpc.gt.0).AND.(inst_gpc.le.(nsim/pccalc)).AND.(gpc%nos.gt.0)) then
        igpc = freeunit()
        open (unit=igpc,file=pcfile,status='new')
        pcopen = .true.
      end if
!
! setup filehandle for temporary TRCV (covariance of torsions) file
      lext = 3
      if ((covcalc.le.nsim).AND.(myrank.eq.masterrank)) then
        do mt=1,nangrps
          call int2str(mt,mtstr,lext)
          trcvfile = 'TRCV_'//mtstr(1:lext)//'.tmp'
          inquire(file=trcvfile,exist=exists)
          if (exists.EQV..true.) then
            itrcv(mt) = freeunit()
            open (unit=itrcv(mt),file=trcvfile,status='old')
            close(unit=itrcv(mt),status='delete')
          end if
!
          if (do_pol(moltypid(molangr(mt,1))).EQV..false.) cycle
          itrcv(mt) = freeunit()
          open (unit=itrcv(mt),file=trcvfile,status='new')
          trcvopen(mt) = .true.
        end do
      end if
!
! setup filehandle for on-the-fly RE trace (original structure vs. condition)
      if ((myrank.eq.masterrank).AND.(pdb_analyze.EQV..false.).AND.(use_MPIMultiSeed.EQV..true.)) then
        inquire(file=retrfile,exist=exists)
        if (exists.EQV..true.) then
          iretr = freeunit()
          open(unit=iretr,file=retrfile,status='old',position='append')
          if (do_restart.EQV..false.) close(unit=iretr,status='delete')
        end if
        inquire(file=retrfile,exist=exists)
        if ((re_freq.le.nsim).AND.(inst_retr.EQV..true.).AND.((do_restart.EQV..false.).OR.(exists.EQV..false.))) then
          iretr = freeunit()
          open (unit=iretr,file=retrfile,status='new')
          retropen = .true.
        else if ((re_freq.le.nsim).AND.(inst_retr.EQV..true.).AND.(do_restart.EQV..true.).AND.(exists.EQV..true.)) then
          retropen = .true.
          write(ilog,*) 'Warning. Trace file for PIGS run exists and will be appended (this may be inappropriate).'
          close(unit=iretr)
          open(unit=iretr,file=retrfile,status='old',position='append')
        else if ((do_restart.EQV..true.).AND.(exists.EQV..true.)) then
          close(unit=iretr)
        end if
      end if
!
! setup filehandle for frame weights from ELS (appended)
      inquire(file=elsfile,exist=exists)
      if (exists.EQV..true.) then
        hmjam%iwtsfile = freeunit()
        open(unit=hmjam%iwtsfile,file=elsfile,status='old',position='append')
        if ((do_restart.EQV..false.).AND.(hmjam%prtfrmwts.le.nsim)) then
          close(unit=hmjam%iwtsfile,status='delete')
        else
          close(unit=hmjam%iwtsfile)
        end if
      end if
      inquire(file=elsfile,exist=exists)
      if ((hmjam%prtfrmwts.le.nsim).AND.(exists.EQV..false.)) then
        hmjam%iwtsfile = freeunit()
        open (unit=hmjam%iwtsfile,file=elsfile,status='new')
        elsopen = .true.
      else if (hmjam%prtfrmwts.le.nsim) then
        hmjam%iwtsfile = freeunit()
        open(unit=hmjam%iwtsfile,file=elsfile,status='old',position='append')
        elsopen = .true.
      end if
!
    else
      if (fycopen) close(unit=idihed)
      if (polymopen) close(unit=ipolmr)
      if (enopen) close(unit=iene)
      if (accopen) close(unit=iacc)
      if (ensopen) close(unit=iens)
      if (persopen) close(unit=ipers)
      if (sasaopen) close(unit=isasa)
      if (phtopen) close(unit=ipht)
      if (sav_detailopen) close(unit=isav_det)
      if (state_q_open) close(unit=q_state_f)
      
      if (covcalc.le.nsim) then
        do mt=1,nangrps
          if (trcvopen(mt)) close(unit=itrcv(mt))
        end do
      end if
      if (elsopen) close(unit=hmjam%iwtsfile)
      if (pcopen) close(unit=igpc)
    end if
!
  end if
#else
  if(mode .eq. 1) then
    dihedfile  = 'FYC.dat'
    polmrfile  = 'POLYMER.dat'
    holesfile = 'HOLES.dat'
    enfile  = 'ENERGY.dat'
    if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
      accfile  = 'ACCEPTANCE.dat'
    end if
    if (((dyn_mode.ge.2).AND.(dyn_mode.le.5)).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
      ensfile = 'ENSEMBLE.dat'
    end if
    persfile = 'PERSISTENCE.dat'
    sasafile = 'SAV.dat'
    dsspfile = 'DSSP_RUNNING.dat'
    pcfile = 'GENERAL_DIS.dat'
    trcvfile = 'TRCV.tmp'
    phtfile = 'PHTIT.dat'
    elsfile = 'ELS_WFRAMES.dat'
    rmsdfile = 'RMSD.dat'
    
    
    if (print_det.eqv..true.) then 
        sav_det_file='SAV_detail.dat' !martin 
        sav_det_file_map='SAV_detail_mapping.dat'
    end if 

    if (do_hsq.eqv..true.) then 
        val_det_file='values_detail.dat' !martin 
        keep_track_q='charge_state.dat' 
        pka_test_file_limits='values_limits.dat'
        HSQ_tries_file="HSQ_tries_summary.dat"
        HSQ_acc_file="HSQ_moves_acc_summary.dat"
        HSQ_acc_ratio="HSQ_moves_acc_ratio.dat"
        HSQ_acc_summary="HSQ_moves_summary.dat"
                
        
        inquire(file=HSQ_tries_file,exist=exists)
        if (exists.EQV..true.) then
            HSQ_tries = freeunit()
            open (unit=HSQ_tries,file=HSQ_tries_file,status='old')
            close(unit=HSQ_tries,status='delete')
        end if
        HSQ_tries  = freeunit()
        open (unit=HSQ_tries,file=HSQ_tries_file,status='new') 
        test_file18open = .true.

        
        inquire(file=HSQ_acc_file,exist=exists)
        if (exists.EQV..true.) then
            HSQ_acc = freeunit()
            open (unit=HSQ_acc,file=HSQ_acc_file,status='old')
            close(unit=HSQ_acc,status='delete')
        end if
        HSQ_acc = freeunit()
        open (unit=HSQ_acc,file=HSQ_acc_file,status='new') 
        test_file19open = .true.

        
        inquire(file=HSQ_acc_ratio,exist=exists)
        if (exists.EQV..true.) then
            HSQ_ratio = freeunit()
            open (unit=HSQ_ratio,file=HSQ_acc_ratio,status='old')
            close(unit=HSQ_ratio,status='delete')
        end if
        HSQ_ratio = freeunit()
        open (unit=HSQ_ratio,file=HSQ_acc_ratio,status='new') 
        test_file20open = .true.

        
        inquire(file=HSQ_acc_summary,exist=exists)
        if (exists.EQV..true.) then
            HSQ_summary = freeunit()
            open (unit=HSQ_summary,file=HSQ_acc_summary,status='old')
            close(unit=HSQ_summary,status='delete')
        end if
        HSQ_summary = freeunit()
        open (unit=HSQ_summary,file=HSQ_acc_summary,status='new') 
        test_file21open = .true.

        
        
        inquire(file=pka_test_file_limits,exist=exists)
            if (exists.EQV..true.) then
                file_limits  = freeunit()
                open (unit=file_limits ,file=pka_test_file_limits,status='old')
                close(unit=file_limits ,status='delete')
            end if
            if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
                file_limits = freeunit()
                open (unit=file_limits,file=pka_test_file_limits,status='new') 
                test_file13open = .true.
            end if 
    
        inquire(file=keep_track_q,exist=exists)
        if (exists.EQV..true.) then
           q_state_f = freeunit()
           open (unit=q_state_f,file=keep_track_q,status='old',position='append')
           close(unit=q_state_f,status='delete')
        end if
        
        if (qs_freq .le. nsim) then
           q_state_f = freeunit()
           open (unit=q_state_f,file=keep_track_q,status='new')
           state_q_open = .true.
        end if
        
        inquire(file=val_det_file,exist=exists)
        if (exists.EQV..true.) then
        val_det = freeunit()
        open (unit=val_det,file=val_det_file,status='old')
        close(unit=val_det,status='delete')
        end if
        val_det = freeunit()
        open (unit=val_det,file=val_det_file,status='new') 
        val_det_file_open = .true.
    end if 
! setup filehandle for dihedral angle output
    inquire(file=dihedfile,exist=exists)
    if (exists.EQV..true.) then
       idihed = freeunit()
       open (unit=idihed,file=dihedfile,status='old')
       close(unit=idihed,status='delete')
    end if
    if (torout .le. nsim) then
       idihed = freeunit()
       open (unit=idihed,file=dihedfile,status='new')
       fycopen = .true.
    end if
!
! setup filehandle for polymeric outputs
    inquire(file=polmrfile,exist=exists)
    if (exists.EQV..true.) then
       ipolmr = freeunit()
       open (unit=ipolmr,file=polmrfile,status='old')
       close(unit=ipolmr,status='delete')
    end if
    if (polout .le. nsim) then
       ipolmr = freeunit()
       open (unit=ipolmr,file=polmrfile,status='new') 
       polymopen = .true.
    end if  
!   
! setup filehandle for holes outputs
    inquire(file=holesfile,exist=exists)
    if (exists.EQV..true.) then
      iholes = freeunit()
      open (unit=iholes,file=holesfile,status='old')
      close(unit=iholes,status='delete')
    end if
    if (holescalc.le.nsim) then
      iholes = freeunit()
      open (unit=iholes,file=holesfile,status='new') 
      holesopen = .true.
    end if 
!   
! setup filehandle for energy outputs
    inquire(file=enfile,exist=exists)
    if (exists.EQV..true.) then
       iene = freeunit()
       open (unit=iene,file=enfile,status='old')
       close(unit=iene,status='delete')
    end if
    if (enout.le.nsim) then
      iene = freeunit()
      open (unit=iene,file=enfile,status='new')     
      enopen = .true.
    end if
!
! setup filehandle for acceptance outputs
    if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
      inquire(file=accfile,exist=exists)
      if (exists.EQV..true.) then
        iacc = freeunit()
        open (unit=iacc,file=accfile,status='old')
        close(unit=iacc,status='delete')
      end if
      if (accout.le.nsim) then
        iacc = freeunit()
        open (unit=iacc,file=accfile,status='new')
        accopen = .true.
      end if
    end if
!
! setup filehandle for ensemble outputs
    if ((pdb_analyze.EQV..false.).AND.(((dyn_mode.ge.2).AND.(dyn_mode.le.5)).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
      inquire(file=ensfile,exist=exists)
      if (exists.EQV..true.) then
        iens = freeunit()
        open (unit=iens,file=ensfile,status='old')
        close(unit=iens,status='delete')
      end if
      if (ensout.le.nsim) then
        iens = freeunit()
        open (unit=iens,file=ensfile,status='new')
        ensopen = .true.
      end if
    end if
!
! setup filehandle for persistence length calculation output file
    inquire(file=persfile,exist=exists)
    if (exists.EQV..true.) then
       ipers = freeunit()
       open (unit=ipers,file=persfile,status='old')
       close(unit=ipers,status='delete')
    end if
    if ((polcalc.le.nsim).AND.(doperser.EQV..true.)) then
      ipers = freeunit()
      open (unit=ipers,file=persfile,status='new')
      persopen = .true.
    end if
!
! setup filehandle for pH titration outputs
    
    if ((phfreq.gt.0.0).AND.(phout.le.nsim)) then
       inquire(file=phtfile,exist=exists)
       if (exists.EQV..true.) then
          ipht = freeunit()
          open (unit=ipht,file=phtfile,status='old')
          close(unit=ipht,status='delete')
       end if
       ipht = freeunit()
       open (unit=ipht,file=phtfile,status='new') 
       phtopen = .true.
    end if
    
      ! martin : try to create my own file
      
      if (print_det.eqv..true.) then 

        inquire(file=sav_det_file,exist=exists)
        if (exists.EQV..true.) then
            isav_det = freeunit()
            open (unit=isav_det,file=sav_det_file,status='old')
            close(unit=isav_det,status='delete')
        end if
      
        isav_det = freeunit()
        open (unit=isav_det,file=sav_det_file,status='new') 
        sav_detailopen = .true.
        
        inquire(file=sav_det_file_map,exist=exists)
        if (exists.EQV..true.) then
            isav_det_map = freeunit()
            open (unit=isav_det_map,file=sav_det_file_map,status='old')
            close(unit=isav_det_map,status='delete')
        end if
      
        isav_det_map = freeunit()
        open (unit=isav_det_map,file=sav_det_file_map,status='new') 
        sav_detail_map_open = .true.
        
      end if 
      
      

    
! setup filehandle for on-the-fly SAV output
    inquire(file=sasafile,exist=exists)
    if (exists.EQV..true.) then
       isasa = freeunit()
       open (unit=isasa,file=sasafile,status='old')
       close(unit=isasa,status='delete')
    end if
    if (savreq%instfreq.le.nsim) then
      isasa = freeunit()
      open (unit=isasa,file=sasafile,status='new')
      sasaopen = .true.
    end if
!
! setup filehandle for on-the-fly DSSP output
    inquire(file=dsspfile,exist=exists)
    if (exists.EQV..true.) then
      idssp = freeunit()
      open (unit=idssp,file=dsspfile,status='old')
      close(unit=idssp,status='delete')
    end if
    if ((dsspcalc.le.nsim).AND.(inst_dssp.EQV..true.)) then
      idssp = freeunit()
      open (unit=idssp,file=dsspfile,status='new')
      dsspopen = .true.
    end if
!
! setup filehandle for on-the-fly distance output
    inquire(file=pcfile,exist=exists)
    if (exists.EQV..true.) then
      igpc = freeunit()
      open (unit=igpc,file=pcfile,status='old')
      close(unit=igpc,status='delete')
    end if
    if ((pccalc.le.nsim).AND.(inst_gpc.gt.0).AND.(inst_gpc.le.(nsim/pccalc)).AND.(gpc%nos.gt.0)) then
      igpc = freeunit()
      open (unit=igpc,file=pcfile,status='new')
      pcopen = .true.
    end if
!
! setup filehandle for on-the-fly RMSD output
    inquire(file=rmsdfile,exist=exists)
    if (exists.EQV..true.) then
      irmsd = freeunit()
      open (unit=irmsd,file=rmsdfile,status='old')
      close(unit=irmsd,status='delete')
    end if
    if ((align%calc.le.nsim).AND.(align%yes.EQV..true.).AND.(align%instrmsd.EQV..true.)) then
      irmsd = freeunit()
      open (unit=irmsd,file=rmsdfile,status='new')
      rmsdopen = .true.
    end if
!
! setup filehandle for temporary TRCV (covariance of torsions) file
    lext = 3
    if (covcalc.le.nsim) then 
      do mt=1,nangrps
        call int2str(mt,mtstr,lext)
        trcvfile = 'TRCV_'//mtstr(1:lext)//'.tmp'
        inquire(file=trcvfile,exist=exists)
        if (exists.EQV..true.) then
          itrcv(mt) = freeunit()
          open (unit=itrcv(mt),file=trcvfile,status='old')
          close(unit=itrcv(mt),status='delete')
        end if
!
        if (do_pol(moltypid(molangr(mt,1))).EQV..false.) cycle
        itrcv(mt) = freeunit()
        open (unit=itrcv(mt),file=trcvfile,status='new')
        trcvopen(mt) = .true.
      end do
    end if
!
! setup filehandle for frame weights from ELS (appended)
    inquire(file=elsfile,exist=exists)
    if (exists.EQV..true.) then
      hmjam%iwtsfile = freeunit()
      open(unit=hmjam%iwtsfile,file=elsfile,status='old',position='append')
      if ((do_restart.EQV..false.).AND.(hmjam%prtfrmwts.le.nsim)) then
        close(unit=hmjam%iwtsfile,status='delete')
      else
        close(unit=hmjam%iwtsfile)
      end if
    end if
    inquire(file=elsfile,exist=exists)
    if ((hmjam%prtfrmwts.le.nsim).AND.(exists.EQV..false.)) then
      hmjam%iwtsfile = freeunit()
      open (unit=hmjam%iwtsfile,file=elsfile,status='new')
      elsopen = .true.
    else if (hmjam%prtfrmwts.le.nsim) then 
      hmjam%iwtsfile = freeunit()
      open(unit=hmjam%iwtsfile,file=elsfile,status='old',position='append')
      elsopen = .true.
    end if
!
! close files and free up filehandles
  else
    if (fycopen) close(unit=idihed)
    if (polymopen) close(unit=ipolmr)
    if (enopen) close(unit=iene)
    if (accopen) close(unit=iacc)
    if (persopen) close(unit=ipers)
    if (ensopen) close(unit=iens)
    if (sasaopen) close(unit=isasa)
    if (dsspopen) close(unit=idssp)
    if (phtopen) close (unit=ipht)
    if (sav_detailopen) close(unit=isav_det)
    if (sav_detail_map_open) close(unit=isav_det)
    if (val_det_file_open) close(unit=val_det)
    if (state_q_open) close(unit=q_state_f)

    if (covcalc.le.nsim) then
      do mt=1,nangrps
        if (trcvopen(mt).EQV..true.) close(unit=itrcv(mt))
      end do
    end if
    if (elsopen) close(unit=hmjam%iwtsfile)
    if (pcopen) close(unit=igpc)
    if (rmsdopen) close(unit=irmsd)
  end if
!
#endif
!
end
!
!
#ifdef ENABLE_MPI
!
subroutine makelogio(mode)
!
  use iounit
  use mcsums
  use mpistuff
!
  implicit none
!
  logical exists
  integer freeunit,lext,mode
  character(60) mclogfile
  character(3) xpont
!
  if (mode.eq.1) then
!
    lext = 3
    call int2str(myrank,xpont,lext) 
    mclogfile ='N_'//xpont(1:lext)//'.log' 
!
!   setup filehandle for log-file
    inquire(file=mclogfile,exist=exists)
    if (exists.EQV..true.) then
      ilog = freeunit()
      open (unit=ilog,file=mclogfile,status='old')
      close(unit=ilog,status='delete')
    end if
    ilog = freeunit()
    open (unit=ilog,file=mclogfile,status='new')
!
  else
!
    inquire(unit=ilog,opened=exists)
    if (exists.EQV..true.) close(unit=ilog)
!
  end if
end
!
#endif
!--------------------------------------------------------------------
!
!
! ####################################################
! ##                                                ##
! ## subroutine contours -- file I/O for dipeptides ##
! ##                                                ##
! ####################################################
!
! "contours" opens or closes output files generating 
! dipeptide maps
! mode = 1: open
! mode = 2: or anything else close
!
!  subroutine contours(mode)
!  implicit none
!  include 'files.i"
!  integer mode,freeunit
!  integer ien,ignorm,is,iphi,ipsi
!  character*60 enfile,gnormfile
!  character*60 phifile,psifile
!  character*60 basinfile
!  logical exists
!c     
!c     open output files for explicit or averaged energies
!  if(mode .eq. 1) then
!     enfile = filename(1:leng)//'En.dat'
!     inquire(file=enfile,exist=exists)
!     if (exists.EQV..true.) then
!        ien = freeunit()
!        open (unit=ien,file=enfile,status='old')
!        close(unit=ien,status='delete')
!     end if
!     ien = freeunit()
!     open (unit=ien,file=enfile,status='new')
!c     
!c     open output files to store gradient norms
!     gnormfile = filename(1:leng)//'Gnorm.dat'
!     inquire(file=gnormfile,exist=exists)
!     if (exists.EQV..true.) then
!        ignorm = freeunit()
!        open (unit=ignorm,file=gnormfile,status='old')
!        close(unit=ignorm,status='delete')
!     end if
!     ignorm = freeunit()
!     open (unit=ignorm,file=gnormfile,status='new')
!c     
!c     open file to store phi angles
!     phifile = filename(1:leng)//'Phi.dat'
!     inquire(file=phifile,exist=exists)
!     if (exists.EQV..true.) then
!        iphi = freeunit()
!        open (unit=iphi,file=phifile,status='old')
!        close(unit=iphi,status='delete')
!     end if
!     iphi = freeunit()
!     open (unit=iphi,file=phifile,status='new')
!c         
!c     open file to store psi angles
!     psifile = filename(1:leng)//'Psi.dat'
!     inquire(file=psifile,exist=exists)
!     if (exists.EQV..true.) then
!c           ipsi = freeunit()
!        open (unit=ipsi,file=psifile,status='old')
!        close(unit=ipsi,status='delete')
!     end if
!     ipsi = freeunit()
!     open (unit=ipsi,file=psifile,status='new')
!c         
!c     open file to store inherent structures
!     basinfile = filename(1:leng)//'IS.dat'
!     inquire(file=basinfile,exist=exists)
!     if (exists.EQV..true.) then
!        is = freeunit()
!        open (unit=is,file=basinfile,status='old')
!        close(unit=is,status='delete')
!     end if
!     is = freeunit()
!     open (unit=is,file=basinfile,status='new')
!  else
!     close(unit=iphi)
!     close(unit=ipsi)
!     close(unit=ien)
!     close(unit=ignorm)
!     close(unit=is)
!  end if
!  return
!  end
!
