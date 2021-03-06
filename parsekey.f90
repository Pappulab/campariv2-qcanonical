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
! CONTRIBUTIONS: Adam Steffen, Albert Mao, Rohit Pappu, Nicolas Bloechliger!
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! the master routine cycles through the key-file extracts keywords and arguments
! and passes them to the parsing-routine (organized loosely by module)
!
!-----------------------------------------------------------------------
!
subroutine parsekey(which)
! 
  use iounit
  use keys
  
!
  implicit none
!
  integer j,next,which,jjj
  character(MAXKWLEN) keyword
  character(MAXKEYLEN) nexten,kline,wo1,wo2,wo3,rest
  integer wt11,wt12,wt21,wt22,wt31,wt32,next2,kt1,kt2
!
! generic format: make sure this matches MAXKEYLEN
   79 format(a200)
!
! parse keywords according to selection     

  do j = 1,nkey
    wo1(:) = ' '
    wo2(:) = ' '
    wo3(:) = ' '
    rest(:) = ' '
    next = 1
    do jjj=1,size(key(j)%line)
      kline(jjj:jjj) = key(j)%line(jjj)
    end do
    do jjj=size(key(j)%line)+1,MAXKEYLEN
      kline(jjj:jjj) = ' '
    end do
    call extract_str(kline,keyword,next)
    call strlims_quiet(keyword,kt1,kt2)
    call toupper(keyword)
    nexten(1:MAXKEYLEN-next+1) = kline(next:MAXKEYLEN)
    next2 = 1
    call extract_str(nexten,wo1,next2)
    call strlims(wo1,wt11,wt12)
    if (((wt11.eq.wt12).AND.(wo1(wt11:wt11).eq.' ')).OR.(wo1(wt11:wt11).eq.'#').OR.(next2.le.1)) then
      if (keyword(1:6).eq.'FMCSC_') then
        if (which.eq.1) write(ilog,*) 'Warning. Keyword ',keyword(kt1:kt2),' has no specification (empty or commented). Ignored.'
      end if
      cycle
    end if
    rest(1:MAXKEYLEN-next2+1) = nexten(next2:MAXKEYLEN)
    call extract_str(nexten,wo2,next2)
    call strlims_quiet(wo2,wt21,wt22)
    call extract_str(nexten,wo3,next2)
    call strlims_quiet(wo3,wt31,wt32)
!    write(*,*) keyword,wo1(wt11:wt12),wo2(wt21:wt22),wo3(wt31:wt32)
!
    call key_params(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 &rest,which,keyword)
    call key_system(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 &rest,which,keyword)
    call key_sequen(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 &rest,which,keyword)
    call key_energies(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 &rest,which,keyword)
    call key_cutoffs(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 &rest,which,keyword)
    call key_movesets(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 &rest,which,keyword)
    call key_mpi(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 &rest,which,keyword)
    call key_pdb(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 &rest,which,keyword)
    call key_torsn(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 &rest,which,keyword)
    call key_polyavg(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 &rest,which,keyword)
    call key_analysis(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 &rest,which,keyword)
!
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine dum2int(dum,choice)
!
  implicit none
!
  RTYPE dum
  integer choice
!
  if (dum.gt.(1.0*HUGE(choice))) then
    choice = HUGE(choice) - 10
  else if (dum.lt.(-1.0*HUGE(choice))) then
    choice = 10 - HUGE(choice)
  else
    choice = INT(dum)
  end if
!
end
!
!-----------------------------------------------------------------------
!
! next are sub-subroutines which read in groups of keywords -> compiles
! more easily and quickly
!
!-----------------------------------------------------------------------
!
! covers: energies,distrest,tabpot,ewalds,dssps(some)
!
subroutine key_energies(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 & rest,which,keyword)
!
  use iounit
  use energies
  use keys
  use params
  use dssps
  use wl  
  use sequen
  use mcsums! Martin added 
  use movesets ! Martin added
!
  implicit none
!
  character(MAXKEYLEN) wo1,wo2,wo3,rest
  character(MAXKWLEN) keyword
  character(MAXKEYLEN) wo1rest
  integer which,choice,wt11,wt12,wt21,wt22,wt31,wt32,hlp1
  logical exists
  
  
  RTYPE dum
!

  exists=.false.
  
  if (which.eq.1) then
!
!
!   scaling factor for repulsive inverse power potential
!
    if (keyword(1:13) .eq. 'FMCSC_SC_IPP ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_IPP = .false.
        scale_IPP = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_IPP. Defaulting &
 &to ',scale_IPP,'.'
        end if
      else
        scale_IPP = dum
        use_IPP = .true.
      end if
!
!   scaling factor for attractive LJ interactions
!
    else if (keyword(1:15) .eq. 'FMCSC_SC_ATTLJ ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_attLJ = .false.
        scale_attLJ = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_ATTLJ. Defaultin&
 &g to ',scale_attLJ,'.'
        end if
      else
        scale_attLJ = dum
        use_attLJ = .true.
      end if
!
!   scaling factor for polar interactions
!
    else if (keyword(1:15) .eq. 'FMCSC_SC_POLAR ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_POLAR = .false.
        scale_POLAR = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_POLAR. Defaultin&
 &g to ',scale_POLAR,'.'
        end if
      else
        scale_POLAR = dum
        use_POLAR = .true.
      end if
!
!   scaling factor for DMFI of ABSINTH model (also turns on charge screening)
!
    else if (keyword(1:16) .eq. 'FMCSC_SC_IMPSOLV ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_IMPSOLV = .false.
        scale_IMPSOLV = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_IMPSOLV. Default&
 &ing to ',scale_IMPSOLV,'.'
        end if
      else
        scale_IMPSOLV = dum
        use_IMPSOLV = .true.
      end if
!
!   scaling factor for WCA potential
!
    else if (keyword(1:12) .eq. 'FMCSC_SC_WCA ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_WCA = .false.
        scale_WCA = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_WCA. Defaulting &
 &to ',scale_WCA,'.'
        end if
      else
        scale_WCA = dum
        use_WCA = .true.
      end if
!
!   scaling factor for hard-coded structural correction terms (discouraged)
!
    else if (keyword(1:15) .eq. 'FMCSC_SC_EXTRA ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_CORR = .false.
        scale_CORR = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_EXTRA. Defaultin&
 &g to ',scale_CORR,'.'
        end if
      else
        scale_CORR = dum
        use_CORR = .true.
      end if
!
!   whether to guess bonded terms from initial setup
!
    else if (keyword(1:19) .eq. 'FMCSC_GUESS_BONDED ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        guess_bonded = .true.
      else
        guess_bonded = .false.
      end if
!
!   scaling factors for bonded terms (Cartesian sampling)
!
    else if (keyword(1:18) .eq. 'FMCSC_SC_BONDED_B ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_BOND(1) = .false.
        scale_BOND(1) = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_BONDED_B. Defaul&
 &ting to ',scale_BOND(1),'.'
        end if
      else
        scale_BOND(1) = dum
        use_BOND(1) = .true.
      end if
    else if (keyword(1:18) .eq. 'FMCSC_SC_BONDED_A ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_BOND(2) = .false.
        scale_BOND(2) = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_BONDED_A. Defaul&
 &ting to ',scale_BOND(2),'.'
        end if
      else
        scale_BOND(2) = dum
        use_BOND(2) = .true.
      end if
    else if (keyword(1:18) .eq. 'FMCSC_SC_BONDED_I ') then
      read (wo1(wt11:wt12),*) dum

      if (dum .le. 0.0d0)  then
        use_BOND(3) = .false.
        scale_BOND(3) = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_BONDED_I. Defaul&
 &ting to ',scale_BOND(3),'.'
        end if
      else
        scale_BOND(3) = dum
        use_BOND(3) = .true.
      end if

    else if (keyword(1:18) .eq. 'FMCSC_SC_BONDED_T ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_BOND(4) = .false.
        scale_BOND(4) = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_BONDED_T. Defaul&
 &ting to ',scale_BOND(4),'.'
        end if
      else
        scale_BOND(4) = dum
        use_BOND(4) = .true.
      end if 
    else if (keyword(1:18) .eq. 'FMCSC_SC_BONDED_M ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_BOND(5) = .false.
        scale_BOND(5) = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_BONDED_M. Defaul&
 &ting to ',scale_BOND(5),'.'
        end if
      else
        scale_BOND(5) = dum
        use_BOND(5) = .true.
      end if
!
!   Directory with CMAP files
!
    else if (keyword(1:14) .eq. 'FMCSC_CMAPDIR ') then
      if ((wt12.lt.MAXKEYLEN).AND.(wo1(wt12:wt12).ne.SLASHCHAR)) then
        wo1(wt12+1:wt12+1)=SLASHCHAR
        wt12 = wt12 + 1
      end if
      cm_dir = wo1(wt11:wt12)
!
!   if using B-splines to "inter"polate CMAPS, this sets order
!
    else if (keyword(1:16) .eq. 'FMCSC_CMAPORDER ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.lt.2).OR.(choice.gt.10)) then
        cm_splor = 4
      else
        cm_splor = dum
      end if
!
!   scaling factor for torsional potentials
!
    else if (keyword(1:13) .eq. 'FMCSC_SC_TOR ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_TOR = .false.
        scale_TOR = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_TOR. Defaulting &
 &to ',scale_TOR,'.'
        end if
      else
        scale_TOR = dum
        use_TOR = .true.
      end if
!
!   scaling factor for global secondary structure biasing potential
!
    else if (keyword(1:14) .eq. 'FMCSC_SC_ZSEC ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_ZSEC = .false.
        scale_ZSEC = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_ZSEC. Defaulting&
 & to ',scale_ZSEC,'.'
        end if
      else
        scale_ZSEC = dum
        use_ZSEC = .true.
      end if
!
!   scaling factor for global DSSP biasing potential
!
    else if (keyword(1:14) .eq. 'FMCSC_SC_DSSP ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_DSSP = .false.
        scale_DSSP = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_DSSP. Defaulting&
 & to ',scale_DSSP,'.'
        end if
      else
        scale_DSSP = dum
        use_DSSP = .true.
      end if
!
!   scaling factor for polymeric biasing potentials
!
    else if (keyword(1:14) .eq. 'FMCSC_SC_POLY ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_POLY = .false.
        scale_POLY = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_POLY. Defaulting&
 & to ',scale_POLY,'.'
        end if
      else
        scale_POLY = dum
        use_POLY = .true.
      end if
!
!   scaling factor for tabulated potentials
!
    else if (keyword(1:15) .eq. 'FMCSC_SC_TABUL ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_TABUL = .false.
        scale_TABUL = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_TABUL. Defaultin&
 &g to ',scale_TABUL,'.'
        end if
      else
        scale_TABUL = dum
        use_TABUL = .true.
      end if
!
!   scaling factor for distance restraint potentials
!
    else if (keyword(1:15) .eq. 'FMCSC_SC_DREST ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_DREST = .false.
        scale_DREST = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_DREST. Defaultin&
 &g to ',scale_DREST,'.'
        end if
      else
        scale_DREST = dum
        use_DREST = .true.
      end if
!
#ifdef LINK_NETCDF
!
!   scaling factor for EM restraint potentials
!
    else if (keyword(1:16) .eq. 'FMCSC_SC_EMICRO ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_EMICRO = .false.
        scale_EMICRO = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_EMICRO. Defaultin&
 &g to ',scale_EMICRO,'.'
        end if
      else
        scale_EMICRO = dum
        use_EMICRO = .true.
      end if
!
#endif
!
!   scaling factor for potentials based on linear combinations of sin/cos terms of system torsions, DISABLED FUNCTIONALITY
!
    else if (keyword(1:15) .eq. 'FMCSC_SC_LCTOR ') then
      read (wo1(wt11:wt12),*) dum
      if (dum .le. 0.0d0)  then
        use_LCTOR = .false.
        scale_LCTOR = 0.0
        if (dum.lt.0.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_SC_LCTOR. Defaultin&
 &g to ',scale_LCTOR,'.'
        end if
      else
        scale_LCTOR = dum
        use_LCTOR = .true.
      end if
!
!    whether or not to (potentially) do growth calculation with ghost particles
!
     else if (keyword(1:12) .eq. 'FMCSC_GHOST ') then
       read (wo1(wt11:wt12),*) dum
       call dum2int(dum,choice)
       if (choice.eq.1) then
         use_FEG = .true.
       else
         use_FEG = .false.
       end if
!
!    what dimensions to resculpt energy landscape for (a la accelerated molecular dynamics)
!
     else if (keyword(1:13) .eq. 'FMCSC_SCULPT ') then
       if (allocated(hmjam%lst).EQV..false.) allocate(hmjam%lst(MAXENERGYTERMS))
       if (allocated(hmjam%thresh).EQV..false.) allocate(hmjam%thresh(MAXENERGYTERMS,2))
       hmjam%thresh(:,1) = -1.0*(HUGE(dum)/1000.0)
       hmjam%thresh(:,2) = HUGE(dum)
       if (allocated(hmjam%alpha).EQV..false.) allocate(hmjam%alpha(MAXENERGYTERMS,2))
       hmjam%alpha(:,:) = -1.0
       wo1rest = wo1(wt11:wt12)//' '//rest(1:(MAXKEYLEN-wt12+wt11-2))
       hlp1 = MAXENERGYTERMS
       choice = 1
       call get_ints_from_str(wo1rest,hlp1,hmjam%lst,choice,hmjam%nlst)
       if (hmjam%nlst.le.0) then
         write(ilog,*) 'Invalid or missing spec. for FMCSC_SCULPT. Turning off energy landscape sculpting.'
         do_accelsim = .false.
       else
         do_accelsim = .true.
       end if
!
!   location of beta-circle for ZSEC analysis
!
    else if (keyword(1:15) .eq. 'FMCSC_ZS_POS_B ') then
      read (wo1(wt11:wt12),*) par_ZSEC2(5)
      read (wo2(wt21:wt22),*) par_ZSEC2(6)
      if ((par_ZSEC2(5).lt.-180.0d0).OR.&
 &                   (par_ZSEC2(5).ge.180.0)) then
        par_ZSEC2(5) = -155.0
        write(ilog,*) 'Invalid spec. for FMCSC_ZS_POS_B (phi). Defau&
 &lting to ',par_ZSEC2(5),'.'
      end if
      if ((par_ZSEC2(6).lt.-180.0d0).OR.&
 &                   (par_ZSEC2(6).ge.180.0)) then
        par_ZSEC2(6) = 160.0
        write(ilog,*) 'Invalid spec. for FMCSC_ZS_POS_B (psi). Defau&
 &lting to ',par_ZSEC2(6),'.'
      end if
!
!   radius of beta-circle
!
    else if (keyword(1:15) .eq. 'FMCSC_ZS_RAD_B ') then
      read (wo1(wt11:wt12),*) par_ZSEC2(7)
      if ((par_ZSEC2(7).lt.0.0d0).OR.&
 &                     (par_ZSEC2(7).gt.180.0)) then
        par_ZSEC2(7) = 35.0
        write(ilog,*) 'Invalid spec. for FMCSC_ZS_RAD_B. Defaulting &
 &to ',par_ZSEC2(7),'.'
      end if
!
!   decay parameter around beta-circle
!
    else if (keyword(1:15) .eq. 'FMCSC_ZS_STP_B ') then
      read (wo1(wt11:wt12),*) par_ZSEC2(8)
      if ((par_ZSEC2(8).le.0.0d0).OR.&
 &                     (par_ZSEC2(8).gt.10.0)) then
        par_ZSEC2(8) = 0.002
        write(ilog,*) 'Invalid spec. for FMCSC_ZS_STP_B. Defaulting &
 &to ',par_ZSEC2(8),'.'
      end if
!
!   location of alpha-circle
!
    else if (keyword(1:15) .eq. 'FMCSC_ZS_POS_A ') then
      read (wo1(wt11:wt12),*) par_ZSEC2(1)
      read (wo2(wt21:wt22),*) par_ZSEC2(2)
      if ((par_ZSEC2(1).lt.-180.0d0).OR.&
 &                   (par_ZSEC2(1).ge.180.0)) then
        par_ZSEC2(1) = -60.0
        write(ilog,*) 'Invalid spec. for FMCSC_ZS_POS_A (phi). Defau&
 &lting to ',par_ZSEC2(1),'.'
      end if
      if ((par_ZSEC2(2).lt.-180.0d0).OR.&
 &                   (par_ZSEC2(2).ge.180.0)) then
        par_ZSEC2(2) = -50.0
        write(ilog,*) 'Invalid spec. for FMCSC_ZS_POS_A (psi). Defau&
 &lting to ',par_ZSEC2(2),'.'
      end if
!
!   radius of beta-circle
!
    else if (keyword(1:15) .eq. 'FMCSC_ZS_RAD_A ') then
      read (wo1(wt11:wt12),*) par_ZSEC2(3)
      if ((par_ZSEC2(3).lt.0.0d0).OR.&
 &                   (par_ZSEC2(3).gt.180.0)) then
        par_ZSEC2(3) = 35.0
        write(ilog,*) 'Invalid spec. for FMCSC_ZS_RAD_A. Defaulting &
 &to ',par_ZSEC2(3),'.'
      end if
!
!   decay parameter around alpha-circle
!
    else if (keyword(1:15) .eq. 'FMCSC_ZS_STP_A ') then
      read (wo1(wt11:wt12),*) par_ZSEC2(4)
      if ((par_ZSEC2(4).le.0.0d0).OR.&
 &                     (par_ZSEC2(4).gt.10.0)) then
        par_ZSEC2(4) = 0.002
        write(ilog,*) 'Invalid spec. for FMCSC_ZS_STP_A. Defaulting &
 &to ',par_ZSEC2(4),'.'
      end if
!
!   DSSP minimum hydrogen bond energy (cutoff to discourage fusion)
!
    else if (keyword(1:17) .eq. 'FMCSC_DSSP_MINHB ') then
      read (wo1(wt11:wt12),*) par_DSSP(3)
      if ((par_DSSP(3).lt.-10.0).OR.(par_DSSP(3).gt.-4.0)) then
        par_DSSP(3) = -4.0
        write(ilog,*) 'Invalid spec. for FMCSC_DSSP_MINHB. Defaultin&
 &g to ',par_DSSP(3),'.'
      end if
!
!   perfect (good) hydrogen bond energy
!
    else if (keyword(1:18) .eq. 'FMCSC_DSSP_GOODHB ') then
      read (wo1(wt11:wt12),*) par_DSSP(1)
      if ((par_DSSP(1).gt.-1.0).OR. (par_DSSP(1).lt.-4.0)) then
        par_DSSP(1) = -2.5
        write(ilog,*) 'Invalid spec. for FMCSC_DSSP_GOODHB. Defaulti&
 &g to ',par_DSSP(1),'.'
      end if
!
!   DSSP maximum hydrogen bond energy (cutoff to ensure finite range)
!
    else if (keyword(1:17) .eq. 'FMCSC_DSSP_MAXHB ') then
      read (wo1(wt11:wt12),*) par_DSSP(2)
      if ((par_DSSP(2).lt.-1.0).OR.(par_DSSP(2).ge.0.0)) then
        par_DSSP(2) = -0.5
        write(ilog,*) 'Invalid spec. for FMCSC_DSSP_MAXHB. Defaultin&
 &g to ',par_DSSP(2),'.'
      end if
!
!   DSSP distance cutoff (CA-CA) for HB search
!
    else if (keyword(1:15) .eq. 'FMCSC_DSSP_CUT ') then
      read (wo1(wt11:wt12),*) par_DSSP(4)
      if (par_DSSP(4).lt.5.0) then
        par_DSSP(4) = 10.0
        write(ilog,*) 'Invalid spec. for FMCSC_DSSP_CUT. Defaulting &
 &to ',par_DSSP(4),'.'
      end if
!
!   exponent for HB-energy weight term
!
    else if (keyword(1:15) .eq. 'FMCSC_DSSP_EXP ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.lt.1).OR.(choice.gt.20)) then
        par_DSSP2(1) = 3
        write(ilog,*) 'Invalid spec. for FMCSC_DSSP_EXP. Defaulting &
 &to ',par_DSSP2(1),'.'
     else
        par_DSSP2(1) = choice
      end if
!
!   which mode to use for weighing E- and H-fractions with HB-energies:
!   1) adjust every HB individually to DSSP_GOODHB, 2) adjust net HB-score to min(it,1), 3) adjust 
!   net E/H-score to min(it,1)
!
    else if (keyword(1:16) .eq. 'FMCSC_DSSP_MODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.lt.1).OR.(choice.gt.3)) then
        par_DSSP2(2) = 3
        write(ilog,*) 'Invalid spec. for FMCSC_DSSP_MODE. Defaulting&
 &to ',par_DSSP2(2),'.'
      else
        par_DSSP2(2) = choice
      end if
!
!   which mode of torsional LCTs: 1) = Fourier  2) = direct, DISABLED FUNCTIONALITY 
!
    else if (keyword(1:16) .eq. 'FMCSC_TORLCMODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.2) then
        torlcmode = 2
      else
        torlcmode = 1
        write(ilog,*) 'Invalid spec. for FMCSC_TORLCMODE. Defaulting&
 & to ',torlcmode,'.'
      end if
!
!   dielectric for implicit solvent model and for RF electrostatics
!
    else if (keyword(1:14) .eq. 'FMCSC_IMPDIEL ') then
      read (wo1(wt11:wt12),*) par_IMPSOLV(2)
      if (par_IMPSOLV(2).lt.1.0) then
        par_IMPSOLV(2) = 78.0
        write(ilog,*) 'Invalid spec. for FMCSC_IMPDIEL. Defaulting t&
 &o ',par_IMPSOLV(2),'.'
      end if
!
!   calculation frequency for solvent-accessible volume (anaylsis purposes)
!
    else if (keyword(1:14) .eq. 'FMCSC_SAVCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        savcalc = 500
        write(ilog,*) 'Invalid spec. for FMCSC_SAVCALC. Defaulting t&
 &o ',savcalc,'.'
      else
        savcalc = choice
      end if
!
!   probe radius for SAV (and IMPSOLV)
!
    else if (keyword(1:15) .eq. 'FMCSC_SAVPROBE ') then
      read (wo1(wt11:wt12),*) savprobe
      if (savprobe.lt.0.0) then
        savprobe = 2.5
        par_IMPSOLV(1) = 5.0
        write(ilog,*) 'Invalid spec. for FMCSC_SAVPROBE. Defaulting &
 &to ',savprobe,'.'
      else
        par_IMPSOLV(1) = 2.0*savprobe
      end if
!
!   whether to calculate N**2 loop initially (time saver for big systems)
!
    else if (keyword(1:13) .eq. 'FMCSC_N2LOOP ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        do_n2loop = .true.
      else
        do_n2loop = .false.
      end if
    ! Martin : keyword to force the program to use the parameters for the param file for topology 
    else if (keyword(1:20) .eq. 'FMCSC_USE_PARAM_TOP ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((dum.ge.0).and.(dum.le.2)) then 
          par_top = dum
        else 
          par_top=0
        end if

    else if (keyword(1:15) .eq. 'FMCSC_FOS_CORR ') then !Martin : should be 
        read (wo1(wt11:wt12),*) dum
        if (dum.ne.0.0) then
          fos_corr = .true.
        else
          fos_corr=.false.
        end if
        
    else if (keyword(1:17) .eq. 'FMCSC_DPGRP_CORR ') then !Martin : should be 
        read (wo1(wt11:wt12),*) dum
        if (dum.ne.0.0) then
          dpgrp_corr = .true.
        else
          dpgrp_corr=.false.
        end if
        
    else if (keyword(1:21) .eq. 'FMCSC_FEG_FOS_OFFSET ') then
        read (wo1(wt11:wt12),*) dum
        if (dum.ne.0.0) then
         FEG_FOS_OFFSET = .true.
         FEG_FOS_OFFSET_VAL = dum
         
        else
            FEG_FOS_OFFSET=.false.
        end if
        
    else if (keyword(1:24) .eq. 'FMCSC_FOS_OFFSET_FACTOR ') then
        read (wo1(wt11:wt12),*) dum
        if (dum.ge.0.0) then
            FEG_FOS_OFFSET_FACT = dum
        else
            FEG_FOS_OFFSET_FACT=10./12
        end if
      
        
    else if (keyword(1:14) .eq. 'FMCSC_HS_FREQ ') then !Frequency at which to apply the hamiltonian switch part
        read (wo1(wt11:wt12),*) dum
        if (dum.gt.0.0) then
            hs_freq = dum
            do_hs=.true.
            hs_mode=3
        else
            hs_freq=0.0
            do_hs=.false.
        end if
    else if (keyword(1:21) .eq. 'FMCSC_Q_SWITCH_STEPS ') then !      
          read (wo1(wt11:wt12),*) dum
          call dum2int(dum,choice)

          if (choice.lt.0) then
            hsq_steps = 50
          else
            hsq_steps = choice
          end if
 
        
!Martin 
    else if (keyword(1:19) .eq. 'FMCSC_SAV_DET_FREQ ') then ! 
     
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        sav_det_freq = 0
        print_det=.false.
     else
        sav_det_freq = choice
        print_det=.true.
      end if
  
      else if (keyword(1:19) .eq. 'FMCSC_SAV_DET_MODE ') then ! 
     
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        sav_det_mode = 1
     else
        sav_det_mode = choice
      end if
      
  else if (keyword(1:13) .eq. 'FMCSC_DO_HSQ ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      
      if (choice.eq.1) then
        do_hsq=.true.
      else
        do_hsq=.false.
      end if
      
! #tyler pka: flag if the simulation is doing a pka free energy simulation
    else if (keyword(1:13) .eq. 'FMCSC_DO_PKA ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice .eq. 1) then
        do_pka = .true.
        do_pka_2 = .true.
        do_pka_3 = .true.
      else
        do_pka = .false.
        do_pka_2 = .false.
        do_pka_3 = .false.
      end if

!    else if (keyword(1:17) .eq. 'FMCSC_PKA_SIGMAS ') then
!        read (wo1(wt11:wt12),*) par_PKA(1)
!        if ((par_PKA(1).lt.0.0d0).OR.(par_PKA(1).gt.1.0)) then
!          par_PKA(1) = 0.0
!          write(ilog,*) 'Invalid spec. for FMCSC_PKA_SIGMAS. Defaulting&
! & to ',par_PKA(1),'.'
!        end if
!
!    else if (keyword(1:17) .eq. 'FMCSC_PKA_POLARS ') then
!        read (wo1(wt11:wt12),*) par_PKA(2)
!        if ((par_PKA(2).lt.0.0d0).OR.(par_PKA(2).gt.1.0)) then
!          par_PKA(2) = 0.0
!          write(ilog,*) 'Invalid spec. for FMCSC_PKA_POLARS. Defaulting&
! & to ',par_PKA(2),'.'
!        end if
!    else if (keyword(1:17) .eq. 'FMCSC_PKA_FOS ') then
!        read (wo1(wt11:wt12),*) par_PKA(3)
!        if ((par_PKA(2).lt.0.0d0).OR.(par_PKA(3).gt.1.0)) then
!          par_PKA(3) = 0.0
!          write(ilog,*) 'Invalid spec. for FMCSC_PKA_POLARS. Defaulting&
! & to ',par_PKA(2),'.'
!        end if
    end if
!
  else if (which.eq.2) then
!
     call key_energies_l2(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 & rest,which,keyword)
!
!
!
  else if (which.eq.3) then
!
     call key_energies_l3(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 & rest,which,keyword)
!
!
!
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! covers: level 2-parts of energies (incl. distrest,tabpot)
!
subroutine key_energies_l2(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,&
 &wt32,rest,which,keyword)
!
  use iounit
  use energies
  use keys
  use mcgrid
  use distrest
  use tabpot
  use mpistuff
  use movesets
  use system
  use torsn
  use math
  use dssps
  use params
  use ems
  use fos
  use wl
  use mcsums! Martin added 
  use sequen!Martin added
  
  implicit none
!
  character(MAXKEYLEN) wo1,wo2,wo3,rest,wo1rest
  character(MAXKWLEN) keyword
  integer which,choice,wt11,wt12,wt21,wt22,wt31,wt32,hlp1,hlp2
  RTYPE dum
  logical exists
!
  hlp2 = 1
!
  if (which.eq.2) then
      
      if ((do_hs.eqv..true.).or.(do_hsq.eqv..true.)) then 
        if (keyword(1:21) .eq. 'FMCSC_Q_SWITCH_STEPS ') then !      
          read (wo1(wt11:wt12),*) dum
          call dum2int(dum,choice)

          if (choice.lt.0) then
            hsq_steps = 0
          else
            hsq_steps = choice
          end if
        end if 
      end if 
      
      
      if (do_hs.eqv..true.) then 
       if (keyword(1:14) .eq. 'FMCSC_HS_MODE ') then !HS mode : 1 is from file, 2 is everything, 3 is NA+ CL- only
            !hamiltonian applies on
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.gt.0).and.(choice.le.3)) then 
            hs_mode =  dum 
        else 
            hs_mode =  3
        end if 
        
    else if (keyword(1:14) .eq. 'FMCSC_HS_FILE ') then !pointer to the file listing the residue that the alternate 
            !hamiltonian applies on
      inquire(file=wo1(wt11:wt12),exist=exists)
      if (exists.EQV..true.) then
            hs_seq_file =  wo1(wt11:wt12)
        end if
        
    else if (keyword(1:14) .eq. 'FMCSC_HS_FOS ') then !Potential in the alternate state (only applies to specified residues)
      read (wo1(wt11:wt12),*) dum
        if ((dum.ge.0).and.(dum.le.1)) then
            hs_FOS = dum
        else 
            hs_FOS=1.
        end if
        
        
        
    else if (keyword(1:14) .eq. 'FMCSC_HS_ATTLJ ') then !Potential in the alternate state (only applies to specified residues)
    
      read (wo1(wt11:wt12),*) dum
        if ((dum.ge.0).and.(dum.le.1)) then
            hs_ATTLJ = dum
        else 
            hs_ATTLJ=1.
        end if
        
        
    else if (keyword(1:14) .eq. 'FMCSC_HS_REPLJ ') then !Potential in the alternate state (only applies to specified residues)
      read (wo1(wt11:wt12),*) dum
        if ((dum.ge.0).and.(dum.le.1)) then
            hs_REPLJ = dum
        else 
            hs_REPLJ=1.
        end if
        
        
    else if (keyword(1:14) .eq. 'FMCSC_HS_POLAR ') then !Potential in the alternate state (only applies to specified residues)
      read (wo1(wt11:wt12),*) dum
        if ((dum.ge.0).and.(dum.le.1)) then
            hs_POLAR = dum
        else 
            hs_POLAR=1.
        end if
    end if 
  end if
      
!   some spec.s related to energy terms to be used 
    if ((do_pka.eqv..true.).or.(do_hsq.eqv..true.)) then


        if (keyword(1:19) .eq. 'FMCSC_PKA_SOFTCORE ') then
          read (wo1(wt11:wt12),*) dum
          call dum2int(dum,choice)
          
          if (choice .eq. 1) then
            use_softcore = .true.
          else
            use_softcore = .false.
          end if
          
        end if 

        if (keyword(1:16) .eq. 'FMCSC_PKA_DEBUG ') then
          read (wo1(wt11:wt12),*) dum
          call dum2int(dum,choice)
          if (choice .eq. 1) then
            debug_pka = .true.
          else
            debug_pka = .false.
          end if
        end if 
        
        if (do_hsq.eqv..true.) then 
           if (keyword(1:15) .eq. 'FMCSC_HSQ_MODE ') then
              read (wo1(wt11:wt12),*) dum
              call dum2int(dum,choice)

              if (choice.eq.1) then
                hsq_mode=choice
              else if (choice.eq.2) then
                hsq_mode=choice
              else 
                hsq_mode=1
              end if

            else if (keyword(1:22) .eq. 'FMCSC_HSQ_READ_CHARGE ') then !      
              read (wo1(wt11:wt12),*) dum
              call dum2int(dum,choice)

              if (choice.ne.1) then
                read_charge_state = .false.
              else 
                  read_charge_state = .true.
              end if    

            else if (keyword(1:20) .eq. 'FMCSC_Q_SWITCH_FREQ ') then ! 
              read (wo1(wt11:wt12),*) dum
              if (dum.lt.0) then
                qs_freq = 0
                do_hsq=.false.
              else
                qs_freq = dum
                do_hsq=.true.
                do_pka=.false.
              end if
            else if (keyword(1:13) .eq. 'FMCSC_BATH_E ') then ! 
              read (wo1(wt11:wt12),*) dum
              if (dum.lt.0) then
                ener_bath = 0
              else
                ener_bath= dum
              end if

            else if (keyword(1:16) .eq. 'FMCSC_HSQ_NET_Q ') then 

              read (wo1(wt11:wt12),*) dum
              call dum2int(dum,choice)
              tot_QX = choice

            else if (keyword(1:15) .eq. 'FMCSC_Q_TOT_UP ') then !      Basically, this is like FMCSC_HSQ_NET_Q , but for constant FOS
              read (wo1(wt11:wt12),*) dum
              call dum2int(dum,choice)

                if (choice.lt.0) then
                    TOT_X = 0.
                else
                    TOT_X = choice
                end if    

            else if (keyword(1:19) .eq. 'FMCSC_HSQ_N_CHANGE ') then !      
              read (wo1(wt11:wt12),*) dum
              call dum2int(dum,choice)

                if (choice.lt.0) then
                    MAX_SWITCH = 0
                else
                    MAX_SWITCH = choice
                end if    

            else if (keyword(1:15) .eq. 'FMCSC_HSQ_FREQ ') then !      
              read (wo1(wt11:wt12),*) hsq_freq

              if (hsq_freq.lt.0) then
                hsq_freq = 0
              end if 
           else if (keyword(1:19) .eq. 'FMCSC_HSQ_ALT_MODE ') then !      
              read (wo1(wt11:wt12),*) dum
              call dum2int(dum,choice)
              
              if ((choice.lt.0).or.(choice.gt.4)) then
                  alt_mode = 1
              else 
                alt_mode = choice
              end if 
            end if 
        end if 

        if (do_pka.eqv..true.) then 
            if (keyword(1:18) .eq. 'FMCSC_PKA_REP_POT ') then
                read (wo1(wt11:wt12),*) dum
                call dum2int(dum,choice)
                if (dum.ne.0.0) then
                  E_rep_offset_val= choice
                  E_rep_offset= .true.
                else
                  E_rep_offset=.false.
                end if
            end if 
        end if
     end if 
     
    if ((use_IPP.EQV..true.).OR.(use_REMC.EQV..true.)) then  !if REMC we don't fully know what energy terms we get
!
!     exponent for IPP
!
      if (keyword(1:13) .eq. 'FMCSC_IPPEXP ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.le.0)  then
          nindx = 12
          write(ilog,*) 'Invalid spec. for FMCSC_IPPEXP. Defaulting to ',nindx,'.'
          nhalf = nindx/2
        else if (choice.gt.100)  then
          nindx = 2
          nhalf = 1
          use_hardsphere = .true.
        else
          nindx = choice
          if (mod(nindx,2).ne.0) then
            nindx = nindx + 1
            write(ilog,*) 'Warning. FMCSC_IPPEXP needs to be even. Increasing by 1 to give ',nindx,'.'
          end if
          nhalf = nindx/2
        end if
!
      end if
!
    end if !(end if for use_IPP)
!   some spec.s related to energy terms to be used 
!
    if ((use_BOND(1).EQV..true.).OR.(use_BOND(2).EQV..true.).OR.(use_BOND(3).EQV..true.).OR.&
 &      (use_BOND(4).EQV..true.).OR.(use_BOND(5).EQV..true.).OR.(use_REMC.EQV..true.)) then 
!
!     patch file for bonded interactions
!
      if (keyword(1:17) .eq. 'FMCSC_BPATCHFILE ') then
        bpatchfile = wo1(wt11:wt12)
!
      end if
!
    end if !(end if for any bonded term)
!
    if ((use_BOND(3).EQV..true.).OR.(use_REMC.EQV..true.)) then  !if REMC we don't fully know what energy terms we get
!
!     force field convention-based overwrite option to more easily conform with impropers
!
      if (keyword(1:20) .eq. 'FMCSC_IMPROPER_CONV ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (dum.eq.2)  then
          improper_conv(1) = 3
          improper_conv(2) = 1
        end if
!
      end if
!
    end if !(end if for use_BOND(3))
!
    if ((use_POLAR.EQV..true.).OR.(use_REMC.EQV..true.)) then  !if REMC we don't fully know what energy terms we get
!
!     a way to override chemical equivalence assumed in biotype spec. for primary amides (universal)
!
      if (keyword(1:15) .eq. 'FMCSC_AMIDEPOL ') then
        read (wo1(wt11:wt12),*) dum
        primamide_cis_chgshft = dum
!
!     tolerance setting for neutral groups
!
      else if (keyword(1:13) .eq. 'FMCSC_POLTOL ') then
        read (wo1(wt11:wt12),*) dum
        if (dum.lt.0.0) then
          dpgrp_neut_tol = 0.0
          write(ilog,*) 'Invalid spec. for FMCSC_POLTOL. Defaulting to precision limit.'
        else
          dpgrp_neut_tol = dum
        end if
!
!     a way to override residue-level charge flags (used in certain cutoff treatments)
!
      else if (keyword(1:18) .eq. 'FMCSC_NCPATCHFILE ') then
        ncpatchfile = wo1(wt11:wt12)
!
!     a way to override charges at atomic resolution
!
      else if (keyword(1:17) .eq. 'FMCSC_CPATCHFILE ') then
        cpatchfile = wo1(wt11:wt12)
!
      end if
!
    end if !(end if for use_POLAR)
!
    if ((use_WCA.EQV..true.).OR.(use_REMC.EQV..true.)) then
!
!     scaling-lambda for attraction spliced in within WCA
!
      if (keyword(1:14) .eq. 'FMCSC_ATT_WCA ') then
        read (wo1(wt11:wt12),*) dum
        if ((dum .lt. 0.0d0)) then! .or.(dum .gt. 1.0d0))  then
           par_WCA(2) = 0.0d0
           write(ilog,*) 'Invalid spec. for FMCSC_ATT_WCA. Defa&
 &ulting to ',par_WCA(2),'.'
        else
           par_WCA(2) = dum
        end if
!
!     cutoff for WCA
!
      else if (keyword(1:14) .eq. 'FMCSC_CUT_WCA ') then
        read (wo1(wt11:wt12),*) dum
        if (dum .lt. 1.5d0)  then
           par_WCA(1) = 1.5d0
           write(ilog,*) 'Invalid spec. for FMCSC_CUT_WCA. Defa&
 &ulting to ',par_WCA(1),'.'
        else
           par_WCA(1) = dum
        end if
!       update the two dependent parameters 
        par_WCA(3) = pi/(par_WCA(1)**2 - ROOT26*ROOT26) !attr. param. 1
        par_WCA(4) = pi - par_WCA(3)*ROOT26*ROOT26      !attr. param. 2
      end if
!
    end if !(end if for use_WCA)
!
    if ((use_IMPSOLV.EQV..true.).OR.(use_REMC.EQV..true.)) then
!
!     contact dielectric for (partially) distance-dependent screening models
!
      if (keyword(1:18) .eq. 'FMCSC_CONTACTDIEL ') then
        read (wo1(wt11:wt12),*) par_IMPSOLV(8)
        if (par_IMPSOLV(8).lt.1.0) then
          par_IMPSOLV(8) = 5.0
          write(ilog,*) 'Invalid spec. for FMCSC_CONTACTDIEL. Defaul&
 &ting to ',par_IMPSOLV(8),'.'
          par_IMPSOLV(8) = 1./par_IMPSOLV(8)
        else
          par_IMPSOLV(8) = 1./par_IMPSOLV(8)
        end if
!
!     steepness for sigmoidal interpolation for FOS-term (the smaller the steeper)
!
      else if (keyword(1:13) .eq. 'FMCSC_FOSTAU ') then
        read (wo1(wt11:wt12),*) par_IMPSOLV(3)
        if ((par_IMPSOLV(3).gt.100.0).OR.(par_IMPSOLV(3).le.0.0))&
 &then
          par_IMPSOLV(3) = 0.25
          write(ilog,*) 'Invalid spec. for FMCSC_FOSTAU. Defaulting &
 &to ',par_IMPSOLV(3),'.'
        end if
!
!     center for sigmoidal interpolation for FOS term 
!
      else if (keyword(1:13) .eq. 'FMCSC_FOSMID ') then
        read (wo1(wt11:wt12),*) par_IMPSOLV(6)
        if ((par_IMPSOLV(6).gt.1.0).OR.(par_IMPSOLV(6).lt.0.0))&
 &then
          par_IMPSOLV(6) = 0.1
          write(ilog,*) 'Invalid spec. for FMCSC_FOSMID. Defaulting &
 &to ',par_IMPSOLV(6),'.'
        end if
!
!     a toggle to enable T-dependent rFOS values
!
      else if (keyword(1:14) .eq. 'FMCSC_FOSMODE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.gt.1).OR.(choice.lt.0)) then
           use_Tdepfosvals = .false.
          write(ilog,*) 'Invalid spec. for FMCSC_FOSMODE. Defaulting to fixed reference free energy of solvation values.'
        else
          if (choice.eq.1) use_Tdepfosvals = .true.
        end if
!
!     reference temperature to assume for the rFOS values 
!
      else if (keyword(1:14) .eq. 'FMCSC_FOSREFT ') then
        read (wo1(wt11:wt12),*) fos_Tdepref
        if ((fos_Tdepref.gt.1000.0).OR.(fos_Tdepref.lt.10.0)) then
          fos_Tdepref = 298.
          write(ilog,*) 'Invalid spec. for FMCSC_FOSREFT. Defaulting to ',fos_Tdepref,'.'
        end if
!
!     which model to use for charge-screening:
!     pure environmental: 1) group-consistent 2) purely atom-based, 3) #1 + distance-dependence, 4) pure distance-dependent
!     5) generalizd group-consistent 6) generalized atom-based, 7) #4 + d.d., 8) #6 + d.d., 9) #2 + d.d.
!
      else if (keyword(1:15) .eq. 'FMCSC_SCRMODEL ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.gt.9).OR.(choice.lt.1)) then
          scrq_model = 2
          write(ilog,*) 'Invalid spec. for FMCSC_SCRMODEL. Defaultin&
 &g to ',scrq_model,'.'
        else
          scrq_model = choice
        end if
!
!     steepness parameter for sigmoidal interpolation (the smaller the steeper) for q-SCR term
!
      else if (keyword(1:13) .eq. 'FMCSC_SCRTAU ') then
        read (wo1(wt11:wt12),*) par_IMPSOLV(4)
        if ((par_IMPSOLV(4).gt.100.0).OR.(par_IMPSOLV(4).le.0.0))&
 &then
          par_IMPSOLV(4) = 0.5
          write(ilog,*) 'Invalid spec. for FMCSC_SCRTAU. Defaulting &
 &to ',par_IMPSOLV(4),'.'
        end if
!
!     center for sigmoidal interpolation for q-SCR term
!
      else if (keyword(1:13) .eq. 'FMCSC_SCRMID ') then
        read (wo1(wt11:wt12),*) par_IMPSOLV(7)
        if ((par_IMPSOLV(7).gt.1.0).OR.(par_IMPSOLV(7).lt.0.0))&
 &then
          par_IMPSOLV(7) = 0.9
          write(ilog,*) 'Invalid spec. for FMCSC_SCRMID. Defaulting &
 &to ',par_IMPSOLV(7),'.'
        end if
!
!     for the mixed model, the impact the distance-dependence can have
!     ranges from 0.0 (identical to model 1) to 1.0 (maximum impact, i.e.,
!     well-depth at contact minimum is independent of environment in
!     most cases
!
      else if (keyword(1:13) .eq. 'FMCSC_SCRMIX ') then
        read (wo1(wt11:wt12),*) par_IMPSOLV(9)
        if ((par_IMPSOLV(9).gt.1.0).OR.(par_IMPSOLV(9).le.0.0))&
 &then
          par_IMPSOLV(9) = 0.5
          write(ilog,*) 'Invalid spec. for FMCSC_SCRMIX. Defaulting &
 &to ',par_IMPSOLV(9),'.'
        end if
!
!     for the generalized screening models, a true mean is taken for the two
!     atomic, inverse pseudo-dielectrics to compute every interaction
!     this mean can have any kind of order m, including 0 for the geometric mean
!     mean(r1,r2) = (0.5*(r1^m + r2^m))^(1./m)
!
      else if (keyword(1:11) .eq. 'FMCSC_ISQM ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.gt.10).OR.(choice.lt.-10)) then
          i_sqm = 0
          write(ilog,*) 'Specified value for FMCSC_ISQM is out of range. Defau&
 &lting to ',i_sqm,' (restrictions apply, because this is used in the inner loops).'
        else
          i_sqm = choice
        end if
!
!     a way to override solvation group assignments, weights, and ref. values
!
      else if (keyword(1:19) .eq. 'FMCSC_FOSPATCHFILE ') then
        fospatchfile = wo1(wt11:wt12)
!
!     a way to override volume reduction factors
!
      else if (keyword(1:19) .eq. 'FMCSC_ASRPATCHFILE ') then
        ardpatchfile = wo1(wt11:wt12)
!
!     a way to override maximum solvent-accessible volume fractions
!
      else if (keyword(1:19) .eq. 'FMCSC_SAVPATCHFILE ') then
        asmpatchfile = wo1(wt11:wt12)
!
!     DISABLED FUNCTIONALITY
!
      else if (keyword(1:16) .eq. 'FMCSC_COMPRESS1 ') then
!        read (wo1(wt11:wt12),*) par_IMPSOLV(10)
!
!     DISABLED FUNCTIONALITY
!
      else if (keyword(1:16) .eq. 'FMCSC_COMPRESS2 ') then
!        read (wo1(wt11:wt12),*) par_IMPSOLV(11)
      end if
!
    end if !(end if for use_IMPSOLV)
!
    if ((use_TOR.EQV..true.).OR.(use_REMC.EQV..true.)) then
!
!     input file for torsional bias potentials
!
      if (keyword(1:14) .eq. 'FMCSC_TORFILE ') then
        inquire(file=wo1(wt11:wt12),exist=exists)
        if (exists.EQV..true.) then
          torfile = wo1(wt11:wt12)
        else
          write(ilog,*) 'Warning. Cannot open spec.d input file (',&
 &wo1(wt11:wt12),') for FMCSC_TORFILE. Turning off file-based torsional bias terms (FMCSC_SC_TOR).'
          use_TOR = .false.
          scale_TOR = 0.0
        end if
!
!     whether to write out an overview of the torsional bias potentials in the system
!
      else if (keyword(1:16) .eq. 'FMCSC_TORREPORT ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.1) then
          tor_report = .true.
        else
          tor_report = .false.
        end if
!
      end if
!
    end if !(end if for use_TOR)
!
    if ((use_ZSEC.EQV..true.).OR.(use_REMC.EQV..true.)) then
!
!     target content for alpha
!
      if (keyword(1:14) .eq. 'FMCSC_ZS_FR_A ') then
        read (wo1(wt11:wt12),*) par_ZSEC(1)
        if ((par_ZSEC(1).lt.0.0d0).OR.(par_ZSEC(1).gt.1.0)) then
          par_ZSEC(1) = 0.45
          write(ilog,*) 'Invalid spec. for FMCSC_ZS_FR_A. Defaulting&
 & to ',par_ZSEC(1),'.'
        end if
!
!     harmonic force constant for potential on z_alpha (kcal/mol)
!
      else if (keyword(1:15) .eq. 'FMCSC_ZS_FR_KA ') then
        read (wo1(wt11:wt12),*) par_ZSEC(2)
        if (par_ZSEC(2).lt.0.0d0) then
          par_ZSEC(2) = 0.0
          write(ilog,*) 'Invalid spec. for FMCSC_ZS_FR_KA. Defaultin&
 &g to ',par_ZSEC(2),'.'
        end if
!
!     Minimum length for z_alpha (number of residues) #tyler
!
      else if (keyword(1:15) .eq. 'FMCSC_ZS_FR_LA ') then
        read (wo1(wt11:wt12),*) par_ZSEC_L(1)
        use_ZSEC_MinL(1)=.true.
        if (par_ZSEC_L(1).lt.2) then
          use_ZSEC_MinL(1)=.false.
          write(ilog,*) 'Invalid spec. for FMCSC_ZS_FR_LA. Defaultin&
 &g to turning Minimum Length Off'
        end if
!
!     target content for beta
!
      else if (keyword(1:14) .eq. 'FMCSC_ZS_FR_B ') then
        read (wo1(wt11:wt12),*) par_ZSEC(3)
        if ((par_ZSEC(3).lt.0.0d0).OR.(par_ZSEC(3).gt.1.0)) then
          par_ZSEC(3) = 0.45
          write(ilog,*) 'Invalid spec. for FMCSC_ZS_FR_B. Defaulting&
 & to ',par_ZSEC(3),'.'
        end if
!
!     harmonic force constant for potential on z_beta (kcal/mol)
!
      else if (keyword(1:15) .eq. 'FMCSC_ZS_FR_KB ') then
        read (wo1(wt11:wt12),*) par_ZSEC(4)
        if (par_ZSEC(4).lt.0.0d0) then
          par_ZSEC(4) = 0.0
          write(ilog,*) 'Invalid spec. for FMCSC_ZS_FR_KB. Defaultin&
 &g to ',par_ZSEC(4),'.'
        end if
!
!     Minimum length for z_alpha (number of residues) #tyler
!
      else if (keyword(1:15) .eq. 'FMCSC_ZS_FR_LB ') then
        read (wo1(wt11:wt12),*) par_ZSEC_L(2)
        use_ZSEC_MinL(2)=.true.
        if (par_ZSEC_L(2).lt.2) then
          use_ZSEC_MinL(2)=.false.
          write(ilog,*) 'Invalid spec. for FMCSC_ZS_FR_LB. Defaultin&
 &g to turning Minimum Length Off'
        end if
!
      end if
!
    end if !(end if for use_ZSEC)
!
    if ((use_DSSP.EQV..true.).OR.(use_REMC.EQV..true.)) then
!
!     target content for H (alpha)
!
      if (keyword(1:15) .eq. 'FMCSC_DSSP_HSC ') then
        read (wo1(wt11:wt12),*) par_DSSP(7)
        if ((par_DSSP(7).lt.0.0d0).OR.(par_DSSP(7).gt.1.0)) then
          par_DSSP(7) = 0.45
          write(ilog,*) 'Invalid spec. for FMCSC_DSSP_HSC. Defaultin&
 &g to ',par_DSSP(7),'.'
        end if
!
!     harmonic force constant for potential on H-score (kcal/mol)
!
      else if (keyword(1:17) .eq. 'FMCSC_DSSP_HSC_K ') then
        read (wo1(wt11:wt12),*) par_DSSP(8)
        if (par_DSSP(8).lt.0.0d0) then
          par_DSSP(8) = 0.0
          write(ilog,*) 'Invalid spec. for FMCSC_DSSP_HSC_K. Default&
 &ing to ',par_DSSP(8),'.'
        end if
!
!     target content for E (beta)
!
      else if (keyword(1:15) .eq. 'FMCSC_DSSP_ESC ') then
        read (wo1(wt11:wt12),*) par_DSSP(9)
        if ((par_DSSP(9).lt.0.0d0).OR.(par_DSSP(9).gt.1.0)) then
          par_DSSP(9) = 0.45
          write(ilog,*) 'Invalid spec. for FMCSC_DSSP_ESC. Defaultin&
 &g to ',par_DSSP(9),'.'
        end if
!
!     harmonic force constant for potential on E-score (kcal/mol)
!
      else if (keyword(1:17) .eq. 'FMCSC_DSSP_ESC_K ') then
        read (wo1(wt11:wt12),*) par_DSSP(10)
        if (par_DSSP(10).lt.0.0d0) then
          par_DSSP(10) = 0.0
          write(ilog,*) 'Invalid spec. for FMCSC_DSSP_ESC_K. Default&
 &ing to ',par_DSSP(10),'.'
        end if
!
      end if
!
    end if !(end if for use_DSSP)
!
    if ((use_POLY.EQV..true.).OR.(use_REMC.EQV..true.)) then
!
!     input file for polymeric bias potentials
!
      if (keyword(1:15) .eq. 'FMCSC_POLYFILE ') then
        inquire(file=wo1(wt11:wt12),exist=exists)
        if (exists.EQV..true.) then
          polfile = wo1(wt11:wt12)
        else
          write(ilog,*) 'Warning. Cannot open spec.d input file (',&
 &wo1(wt11:wt12),') for FMCSC_POLYFILE. Turning off polymeric biasing terms (FMCSC_SC_POLY).'
          use_POLY = .false.
          scale_POLY = 0.0
        end if
!
!     whether to write out an overview of the polymeric bias potentials in the system
!
      else if (keyword(1:17) .eq. 'FMCSC_POLYREPORT ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.1) then
          poly_report = .true.
        else
          poly_report = .false.
        end if
!
      end if
    end if !(end if for use_POLY)
!
    if ((use_TABUL.EQV..true.).OR.(use_REMC.EQV..true.)) then
!
!     path to file indicating which tab. potential to use for which atom pairs
!
      if (keyword(1:18) .eq. 'FMCSC_TABCODEFILE ') then
!       will check for existence in the reader
        tabcodefile = wo1(wt11:wt12)
!
!     file with actual potential input, includes distance column
!
      else if (keyword(1:17) .eq. 'FMCSC_TABPOTFILE ') then
!       will check for existence in the reader
        tabpotfile = wo1(wt11:wt12)
!
!     file with derivatives (tangents), but no distance column
!
      else if (keyword(1:18) .eq. 'FMCSC_TABTANGFILE ') then
!       will check for existence in the reader
        tabtangfile = wo1(wt11:wt12)
!
!     whether to provide a summary of the tabulated interactions
!
      else if (keyword(1:16) .eq. 'FMCSC_TABREPORT ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.1) then
          tabul_report = .true.
        else
          tabul_report = .false.
        end if
!
!     tightness parameter for simplified Kochanek-Bartels spline (when not reading tangents from file)
!
      else if (keyword(1:16) .eq. 'FMCSC_TABITIGHT ') then
        read (wo1(wt11:wt12),*) tbp%ipprm(1)
        if ((tbp%ipprm(1).lt.-1.0).OR.(tbp%ipprm(1).gt.1.0)) then
          tbp%ipprm(1) = 0.0
          write(ilog,*) 'Invalid spec. for FMCSC_TABITIGHT. Defaulti&
 &ng to ',tbp%ipprm(1),'.'
        end if
!
!     left/right bias parameter for simplified Kochanek-Bartels spline (when not reading tangents from file)
!
      else if (keyword(1:15) .eq. 'FMCSC_TABIBIAS ') then
        read (wo1(wt11:wt12),*) tbp%ipprm(2)
        if ((tbp%ipprm(2).lt.-1.0).OR.(tbp%ipprm(2).gt.1.0)) then
          tbp%ipprm(2) = 0.0
          write(ilog,*) 'Invalid spec. for FMCSC_TABIBIAS. Defaulti&
 &ng to ',tbp%ipprm(2),'.'
        end if

!
      end if
!
    end if !(end if for use_TAB)
!
    if ((use_DREST.EQV..true.).OR.(use_REMC.EQV..true.)) then
!
!     input file for distance restraint terms
!
      if (keyword(1:16) .eq. 'FMCSC_DRESTFILE ') then
        inquire(file=wo1(wt11:wt12),exist=exists)
        if (exists.EQV..true.) then
          drestfile = wo1(wt11:wt12)
        else
          write(ilog,*) 'Warning. Cannot open spec.d input file (',&
 &wo1(wt11:wt12),') for FMCSC_DRESTFILE. Turning off file-based distance restraint terms (FMCSC_SC_DREST).'
          use_DREST = .false.
          scale_DREST = 0.0
        end if
!
!     whether to write out an overview of the distance restraint terms in the system
!
      else if (keyword(1:18) .eq. 'FMCSC_DRESTREPORT ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.1) then
          drest_report = .true.
        else
          drest_report = .false.
        end if
!
      end if
    end if !(end if for use_DREST)
!
    if ((use_FEG.EQV..true.).OR.(use_REMC.EQV..true.)) then
!
!     input file for particles to be ghosted
!
      if (keyword(1:14) .eq. 'FMCSC_FEGFILE ') then
        inquire(file=wo1(wt11:wt12),exist=exists)
        if (exists.EQV..true.) then
          fegfile = wo1(wt11:wt12)
        else
          write(ilog,*) 'Warning. Cannot open spec.d input file (',&
 &wo1(wt11:wt12),') for FMCSC_FEGFILE. Turning off ghosting (FMCSC_GHOST).'
          use_FEG = .false.
        end if
!
!     whether to write out an overview of the ghosted particles in the system
!
      else if (keyword(1:16) .eq. 'FMCSC_FEGREPORT ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.1) then
          feg_report = .true.
        else
          feg_report = .false.
        end if
!
!     whether or not (default) to include ghost-ghost interactions in the scaling
!
      else if (keyword(1:15) .eq. 'FMCSC_FEG_MODE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,fegmode)
        if ((fegmode.lt.1).OR.(fegmode.gt.2)) then
          fegmode = 1
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_MODE. Defaultin&
 &g to ',fegmode,'.'
        end if
!
!     mode for LJ scaling
!
      else if (keyword(1:17) .eq. 'FMCSC_FEG_LJMODE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,fegljmode)
        if ((fegljmode.gt.3).OR.(fegljmode.le.0)) then
          fegljmode = 2
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_LJMODE. Default&
 &ing to ',fegljmode,'.'
        end if
!
!     soft-core radius for LJ terms
!
      else if (keyword(1:16) .eq. 'FMCSC_FEG_LJRAD ') then
        read (wo1(wt11:wt12),*) par_FEG2(3)
        if (par_FEG2(3).lt.0.0) then
          par_FEG2(3) = 0.5
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_LJRAD. Defaulti&
 &ng to ',par_FEG2(3),'.'
        end if
!
!     inner exponent for soft-core LJ
!
      else if (keyword(1:18) .eq. 'FMCSC_FEG_LJSCEXP ') then
        read (wo1(wt11:wt12),*) par_FEG2(7)
        if (par_FEG2(7).le.0.0) then
          par_FEG2(7) = 2.0
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_LJSCEXP. Defaul&
 &ting to ',par_FEG2(7),'.'
        end if
!
!     pre-(polynomial)-exponent for FEG-LJ scaling
!
      else if (keyword(1:16) .eq. 'FMCSC_FEG_LJEXP ') then
        read (wo1(wt11:wt12),*) par_FEG2(8)
        if (par_FEG2(8).le.0.0) then
          par_FEG2(8) = 2.0
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_LJEXP. Defaulti&
 &ng to ',par_FEG2(8),'.'
        end if
!
!     mode for Cb scaling
!
      else if (keyword(1:17) .eq. 'FMCSC_FEG_CBMODE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,fegcbmode)
        if ((fegcbmode.gt.2).OR.(fegcbmode.le.0)) then
          fegcbmode = 1
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_CBMODE. Default&
 &ing to ',fegcbmode,'.'
        end if
!
!     soft-core radius for Cb term
!
      else if (keyword(1:16) .eq. 'FMCSC_FEG_CBRAD ') then
        read (wo1(wt11:wt12),*) par_FEG2(4)
        if (par_FEG2(4).lt.0.0) then
          par_FEG2(4) = 0.5
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_CBRAD. Defaulti&
 &ng to ',par_FEG2(4),'.'
        end if
!
!     inner exponent for soft-core Cb
!
      else if (keyword(1:18) .eq. 'FMCSC_FEG_CBSCEXP ') then
        read (wo1(wt11:wt12),*) par_FEG2(11)
        if (par_FEG2(11).le.0.0) then
          par_FEG2(11) = 1.0
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_CBSCEXP. Defaul&
 &ting to ',par_FEG2(11),'.'
        end if
!
!     pre-(polynomial)-exponent for FEG-Cb scaling
!
      else if (keyword(1:16) .eq. 'FMCSC_FEG_CBEXP ') then
        read (wo1(wt11:wt12),*) par_FEG2(12)
        if (par_FEG2(12).le.0.0) then
          par_FEG2(12) = 2.0
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_CBEXP. Defaulti&
 &ng to ',par_FEG2(12),'.'
        end if
!
      end if
!
    end if ! if use_FEG
!
!
! Settings for EM restraints

    if ((use_EMICRO.EQV..true.).OR.(use_REMC.EQV..true.)) then
!
!     threshold signal for input map to consider solute
!
      if (keyword(1:18) .eq. 'FMCSC_EMTHRESHOLD ') then
        read (wo1(wt11:wt12),*) emthreshdensity
!
!     at what level to truncate (homogenize) values below threshold (optional)
!
      else if (keyword(1:17) .eq. 'FMCSC_EMTRUNCATE ') then
        read (wo1(wt11:wt12),*) emtruncate
!
!     what to assume for the optical background (optional, otherwise histogram peak used)
!
      else if (keyword(1:19) .eq. 'FMCSC_EMBACKGROUND ') then
        read (wo1(wt11:wt12),*) emoptbg
!
#ifdef LINK_NETCDF
!
!     input file for restraint EM map
!
      else if (keyword(1:17) .eq. 'FMCSC_EMMAPFILE ') then
        inquire(file=wo1(wt11:wt12),exist=exists)
        if (exists.EQV..false.) then
          write(ilog,*) 'Fatal. Cannot open spec.d file (',wo1(wt11:wt12),') for FMCSC_EMMAPFILE (FMCSC_SC_EMICRO).'
          call fexit()
        else
          emmapfile = wo1(wt11:wt12)
        end if
!
#endif
!
!     whether to flatten the map at the top (the sanity has to be checked elsewhere)
!
      else if (keyword(1:16) .eq. 'FMCSC_EMFLATTEN ') then
        read (wo1(wt11:wt12),*) emflatval
!
!     assumed total mass
!
      else if (keyword(1:16) .eq. 'FMCSC_EMTOTMASS ') then
        read (wo1(wt11:wt12),*) emtotalmass
        if (emtotalmass.le.0.0) then
          emtotalmass = -1.0
          write(ilog,*) 'Invalid spec. for FMCSC_EMTOTMASS. Defaulting&
 & to computed total mass.'
        end if
!
!     mode selector for type of restraint potential (1: instant., 2: ensemble)
!
      else if (keyword(1:13) .eq. 'FMCSC_EMMODE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.2)) then
          empotmode = 1
          write(ilog,*) 'Invalid spec. for FMCSC_EMMODE. Defaulting&
 & to ',empotmode,'.'
        else
          empotmode = choice
        end if
!
!     weight for average split in ensemble-averaged density restraint
!
      else if (keyword(1:16) .eq. 'FMCSC_EMIWEIGHT ') then
        read (wo1(wt11:wt12),*) emiw
        if ((emiw.le.0.0).OR.(emiw.ge.1.0)) then
          emiw = 0.1
          write(ilog,*) 'Invalid spec. for FMCSC_EMIWEIGHT. Defaulting&
 & to ',emiw,'.'
        end if
!
!     mode selector for efficiency heuristics for energy/force calculations (0: off, 1: slice, 2: line, 3: block)
!
      else if (keyword(1:18) .eq. 'FMCSC_EMHEURISTIC ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.0).OR.(choice.gt.3)) then
          emheuristic = 3
          write(ilog,*) 'Invalid spec. for FMCSC_EMHEURISTIC. Defaulting&
 & to ',emheuristic,'.'
        else
          emheuristic = choice
        end if
!
!     increased deltas for input map to lower resolution
!
      else if (keyword(1:15) .eq. 'FMCSC_EMREDUCE ') then
        read (wo1(wt11:wt12),*) fdlts(1)
        read (wo2(wt21:wt22),*) fdlts(2)
        read (wo3(wt31:wt32),*) fdlts(3)
!
      end if
!
    end if ! if use_EMICRO
!
!   Settings for analysis/use torsional linear combinations, currently DISABLED FUNCTIONALITY
!
    if ((torlccalc.le.nsim).OR.(use_REMC.EQV..true.).OR.&
 &                                 (use_LCTOR.EQV..true.)) then
!
!     read-in file with the coefficients for the torsional LCs
! 
      if (keyword(1:16) .eq. 'FMCSC_TORLCFILE ') then
        torlcfile = wo1(wt11:wt12)
        inquire(file=wo1(wt11:wt12),exist=exists)
        if (exists.EQV..false.) then
          write(ilog,*) 'Coefficent file for torsional LCs does not &
 &seem to exist. Turning off analysis/pot. (',wo1(wt11:wt12),').'
          ntorlcs = 0
          torlccalc = nsim + 1
          use_LCTOR = .false.
          scale_LCTOR = 0.0
        end if
!
!     read-in file with the coefficients for the inverse transform
!
      else if (keyword(1:17) .eq. 'FMCSC_TORLCFILE2 ') then
        torlcfile2 = wo1(wt11:wt12)
        use_lctmoves = .true.
        inquire(file=wo1(wt11:wt12),exist=exists)
        if (exists.EQV..false.) then
          write(ilog,*) 'Coefficent file for torsional LCs (inverse)&
 & does not seem to exist. (',wo1(wt11:wt12),').'
           use_lctmoves = .false.
!          ntorlcs = 0
!          torlccalc = nsim + 1
!          use_LCTOR = .false.
!          scale_LCTOR = 0.0
        end if
!
!     resolution for torsional LCs
!
      else if (keyword(1:15) .eq. 'FMCSC_TORLCRES ') then
        read (wo1(wt11:wt12),*) torlc_params(2)
        if (torlc_params(2) .le. 0.0d0) then
          torlc_params(2) = 0.1
          write(ilog,*) 'Invalid spec. for FMCSC_TORLCRES. Defaultin&
 &g to ',torlc_params(2),'.'
        end if
!
!     minimum expected value for torsional LCs
!
      else if (keyword(1:15) .eq. 'FMCSC_TORLCMIN ') then
        read (wo1(wt11:wt12),*) torlc_params(1)
        if (torlc_params(1) .ge. 0.0d0) then
          torlc_params(1) = -5.0
          write(ilog,*) 'Invalid spec. for FMCSC_TORLCMIN. Defaultin&
 &g to ',torlc_params(1),'.'
        end if
!
!     number of LCTs of interest
!
      else if (keyword(1:14) .eq. 'FMCSC_NRTORLC ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.lt.0) then
          ntorlcs = 0
          write(ilog,*) 'Invalid spec. for FMCSC_NRTORLC. Defaulting&
 & to ',ntorlcs,'.'
          torlccalc = nsim + 1
          use_LCTOR = .false.
          scale_LCTOR = 0.0
        else
          ntorlcs = choice
        end if
!
!     minimum wall height (so that LCT moves can identify correct range for LCTs)
!
      else if (keyword(1:19) .eq. 'FMCSC_LCTTHRESHOLD ') then
        read (wo1(wt11:wt12),*) par_LCTOR(3)
!    
      end if
!
    end if !(end if for use of torsional LCs (analysis and potential))
!
    if ((use_REMC.EQV..true.).OR.(use_LCTOR.EQV..true.)) then
!
!    read-in file with the coefficients for the torsional LCs, currently DISABLED FUNCTIONALITY
! 
      if (keyword(1:14) .eq. 'FMCSC_LCTBINS ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.le.0) then
          par_LCTOR2(1) = 100
          write(ilog,*) 'Invalid spec. for FMCSC_NRTORLC. Defaulting&
 & to ',par_LCTOR2(1),'.'
        else
          par_LCTOR2(1) = choice   
        end if
!    
      end if
!
    end if !(end if for torsional LCs (LCT) potential (exclusively))
!
    if (do_accelsim.EQV..true.) then
!
!     lower threshold energies for landscape sculpting (basin filling)
! 
      if (keyword(1:16) .eq. 'FMCSC_ELS_FILLS ') then
        wo1rest = wo1(wt11:wt12)//' '//rest(1:(MAXKEYLEN-wt12+wt11-2))
        hlp1 = MAXENERGYTERMS
!       sanity checks are applied later
        call get_reals_from_str(wo1rest,hlp1,hmjam%thresh(:,1),hlp2,choice)
!
!     upper threshold energies for landscape sculpting (barrier shaving)
!
      else if (keyword(1:17) .eq. 'FMCSC_ELS_SHAVES ') then
        wo1rest = wo1(wt11:wt12)//' '//rest(1:(MAXKEYLEN-wt12+wt11-2))
        hlp1 = MAXENERGYTERMS
!       sanity checks are applied later
        call get_reals_from_str(wo1rest,hlp1,hmjam%thresh(:,2),hlp2,choice)
!
!     alpha parameters for landscape sculpting (basin filling) in kcal/mol
! 
      else if (keyword(1:18) .eq. 'FMCSC_ELS_ALPHA_F ') then
        wo1rest = wo1(wt11:wt12)//' '//rest(1:(MAXKEYLEN-wt12+wt11-2))
        hlp1 = MAXENERGYTERMS
!       sanity checks are applied later
        call get_reals_from_str(wo1rest,hlp1,hmjam%alpha(:,1),hlp2,choice)
!
!     alpha parameters for landscape sculpting (barrier shaving) in kcal/mol
!
      else if (keyword(1:18) .eq. 'FMCSC_ELS_ALPHA_S ') then
        wo1rest = wo1(wt11:wt12)//' '//rest(1:(MAXKEYLEN-wt12+wt11-2))
        hlp1 = MAXENERGYTERMS
!       sanity checks are applied later
        call get_reals_from_str(wo1rest,hlp1,hmjam%alpha(:,2),hlp2,choice)
!    
!     whether to automatically write a file with frames with nonunity weights
!
      else if (keyword(1:24) .eq. 'FMCSC_ELS_PRINT_WEIGHTS ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.gt.0) then
          hmjam%prtfrmwts = choice
        else
          hmjam%prtfrmwts = HUGE(hmjam%prtfrmwts)
          write(ilog,*) 'Invalid spec. for FMCSC_ELS_PRINT_WEIGHTS. Output file ELS_WFRAMES.dat is not produced.'
        end if
!
!     a threshold to discard frames with extremely low weight
!
      else if (keyword(1:20) .eq. 'FMCSC_ELS_THRESHOLD ') then
        read (wo1(wt11:wt12),*) dum
        if (dum.ge.0.0) then
          hmjam%threshwt = dum
        else
          hmjam%threshwt = 0.0
          write(ilog,*) 'Invalid spec. for FMCSC_ELS_THRESHOLD. Defaulting to ',hmjam%threshwt,'.'
        end if
!
      end if
!
    end if !(end if for energy landscape sculpting is true)
!
!
!
  else
!
    write(ilog,*) 'Fatal. Called key_energies_l2(...) with level oth&
 &er than 2. This is most definitely a bug.'
    call fexit()
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! covers: level 3-parts of energies (incl. ewalds)
!
subroutine key_energies_l3(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,&
 &wt32,rest,which,keyword)
!
  use iounit
  use energies
  use keys
  use mpistuff
  use ewalds
  use torsn
  use cutoffs
  use params
!
  implicit none
!
  character(MAXKEYLEN) wo1,wo2,wo3,rest
  character(MAXKWLEN) keyword
  integer which,choice,wt11,wt12,wt21,wt22,wt31,wt32,i,k,iunit
  integer freeunit
  RTYPE dum
  logical exists
!
  if (which.eq.3) then
!
! FEG scaling factors
!
      
    if ((do_pka.eqv..true.).or.(do_hsq.eqv..true.)) then
        if (use_softcore.eqv..true.) then 
            if (keyword(1:25) .eq. 'FMCSC_PKA_SOFTCORE_ALPHA ') then
                read (wo1(wt11:wt12),*) sc_alpha
                if (sc_alpha.lt.0.0) then
                  sc_alpha = 1.0
                  write(ilog,*) 'Invalid spec. for FMCSC_PKA_SOFTCORE_ALPHA. Defaulting t&
            &o ',sc_alpha,'.'
                end if
      
        !Martin : Lambda exponent for the softcore potential. Default 1
            else if (keyword(1:21) .eq. 'FMCSC_PKA_SOFTCORE_N ') then
                read (wo1(wt11:wt12),*) sc_n
                if (sc_n.le.0.0) then 
                  sc_n = 1.0
                  write(ilog,*) 'Invalid spec. for FMCSC_PKA_SOFTCORE_N. Defaulting t&
            &o ',sc_n,'.'
                end if
            end if
            
        end if 
        if (read_charge_state.eqv..true.) then
            if (keyword(1:22) .eq. 'FMCSC_HSQ_CHARGE_FILE ') then
              exists=.false.
              inquire(file=wo1(wt11:wt12),exist=exists)
              if (exists.EQV..false.) then
                write(ilog,*) 'Fatal. Cannot open spec.d HSQ_CHARGE file (',&
         &wo1(wt11:wt12),') for FMCSC_SEQFILE.'
                 read_charge_state = .false.
                call fexit()
              else
                charge_state_file = wo1(wt11:wt12)
              end if
            end if 

        end if 
    end if 
    
    
    if (use_FEG.EQV..true.) then
!
!     opacity for IPP
!
      if (keyword(1:14) .eq. 'FMCSC_FEG_IPP ') then
        read (wo1(wt11:wt12),*) scale_FEGS(1)
        if (scale_FEGS(1).lt.0.0) then
          scale_FEGS(1) = 0.0
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_IPP. Defaulting&
 & to ',scale_FEGS(1),'.'
        else if (scale_FEGS(1).gt.0.0) then
          use_FEGS(1) = .true.
          if (use_IPP.EQV..false.) then
            write(ilog,*) 'Warning. Standard IPP-term has to be enab&
 &led if ghosting with IPP is requested. Turning term on with zero s&
 &caling factor (FMCSC_FEG_IPP).'
            use_IPP = .true.
            scale_IPP = 0.0
          end if
        end if
!
!     opacity for attractive LJ
!
      else if (keyword(1:16) .eq. 'FMCSC_FEG_ATTLJ ') then
        read (wo1(wt11:wt12),*) scale_FEGS(3)
        if (scale_FEGS(3).lt.0.0) then
          scale_FEGS(3) = 0.0
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_ATTLJ. Defaulti&
 &ng to ',scale_FEGS(3),'.'
        else if (scale_FEGS(3).gt.0.0) then
          use_FEGS(3) = .true.
          if (use_attLJ.EQV..false.) then
            write(ilog,*) 'Warning. Standard ATTLJ-term has to be en&
 &abled if ghosting with ATTLJ is requested. Turning term on with ze&
 &ro scaling factor (FMCSC_FEG_ATTLJ).'
            use_attLJ = .true.
            scale_attLJ = 0.0
          end if
        end if
!
!     opacity for polar term
!
      else if (keyword(1:16) .eq. 'FMCSC_FEG_POLAR ') then
        read (wo1(wt11:wt12),*) scale_FEGS(6)
        if (scale_FEGS(6).lt.0.0) then
          scale_FEGS(6) = 0.0
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_POLAR. Defaulti&
 &ng to ',scale_FEGS(6),'.'
        else if (scale_FEGS(6).gt.0.0) then
          use_FEGS(6) = .true.
          if (use_POLAR.EQV..false.) then
            write(ilog,*) 'Warning. Standard POLAR-term has to be en&
 &abled if ghosting with POLAR is requested. Turning term on with ze&
 &ro scaling factor (FMCSC_FEG_POLAR).'
            use_POLAR = .true.
            scale_POLAR = 0.0
          end if
        end if
!
!     opacity for torsional term(s)
!
      else if (keyword(1:19) .eq. 'FMCSC_FEG_BONDED_T ') then
        read (wo1(wt11:wt12),*) scale_FEGS(18)
        if (scale_FEGS(18).lt.0.0) then
          scale_FEGS(18) = 0.0
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_BONDED_T. Defau&
 &lting to ',scale_FEGS(18),'.'
        else if (scale_FEGS(18).gt.0.0) then
          use_FEGS(18) = .true.
          if (use_BOND(4).EQV..false.) then
            write(ilog,*) 'Warning. Standard BONDED_T-term has to be&
 & enabled if ghosting with BONDED_T is requested. Turning term on w&
 &ith zero scaling factor (FMCSC_FEG_BONDED_T).'
            use_BOND(4) = .true.
            scale_BOND(4) = 0.0
          end if
        end if
!
!     opacity for improper dihedral term(s)
!
      else if (keyword(1:19) .eq. 'FMCSC_FEG_BONDED_I ') then
        read (wo1(wt11:wt12),*) scale_FEGS(17)
        if (scale_FEGS(17).lt.0.0) then
          scale_FEGS(17) = 0.0
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_BONDED_I. Defau&
 &lting to ',scale_FEGS(17),'.'
        else if (scale_FEGS(17).gt.0.0) then
          use_FEGS(17) = .true.
          if (use_BOND(3).EQV..false.) then
            write(ilog,*) 'Warning. Standard BONDED_I-term has to be&
 & enabled if ghosting with BONDED_I is requested. Turning term on w&
 &ith zero scaling factor (FMCSC_FEG_BONDED_I).'
            use_BOND(3) = .true.
            scale_BOND(3) = 0.0
          end if
        end if
!
!     opacity for bond angle term(s) -> iffy 
!
      else if (keyword(1:19) .eq. 'FMCSC_FEG_BONDED_A ') then
        read (wo1(wt11:wt12),*) scale_FEGS(16)
        if (scale_FEGS(16).lt.0.0) then
          scale_FEGS(16) = 0.0
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_BONDED_A. Defau&
 &lting to ',scale_FEGS(16),'.'
        else if (scale_FEGS(16).gt.0.0) then
          use_FEGS(16) = .true.
          if (use_BOND(2).EQV..false.) then
            write(ilog,*) 'Warning. Standard BONDED_A-term has to be&
 & enabled if ghosting with BONDED_A is requested. Turning term on w&
 &ith zero scaling factor (FMCSC_FEG_BONDED_A).'
            use_BOND(2) = .true.
            scale_BOND(2) = 0.0
          end if
        end if
!
!     opacity for bond length term(s) -> Ig nobel prize in the works?
!
      else if (keyword(1:19) .eq. 'FMCSC_FEG_BONDED_B ') then
        read (wo1(wt11:wt12),*) scale_FEGS(15)
        if (scale_FEGS(15).lt.0.0) then
          scale_FEGS(15) = 0.0
          write(ilog,*) 'Invalid spec. for FMCSC_FEG_BONDED_B. Defau&
 &lting to ',scale_FEGS(15),'.'
        else if (scale_FEGS(15).gt.0.0) then
          use_FEGS(15) = .true.
          if (use_BOND(1).EQV..false.) then
            write(ilog,*) 'Warning. Standard BONDED_B-term has to be&
 & enabled if ghosting with BONDED_B is requested. Turning term on w&
 &ith zero scaling factor (FMCSC_FEG_BONDED_B).'
            use_BOND(1) = .true.
            scale_BOND(1) = 0.0
          end if
        end if
!
      end if
!
    end if

!
! Ewald spec.s
!
    if (lrel_md.eq.2) then
!
!     order of B-spline interpolation for PME
!
      if (keyword(1:14) .eq. 'FMCSC_BSPLINE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.le.1).OR.(mod(choice,2).ne.0).OR.&
 &                  (choice.gt.20)) then
          splor = 6
          write(ilog,*) 'Invalid spec. for FMCSC_BSPLINE. Defaulting&
 & to ',splor,'.'
        else
          splor = choice
        end if
!
!     1) PME, 2) standard
!
      else if (keyword(1:12) .eq. 'FMCSC_EWALD ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.le.0).OR.(choice.gt.2)) then
          ewald_mode = 1
          write(ilog,*) 'Invalid spec. for FMCSC_EWALD. Defaulting t&
 &o ',ewald_mode,'.'
        else
          ewald_mode = choice
        end if
        if (ewald_mode.eq.1) then
          ewald_mode = 2
#ifdef LINK_FFTW
          ewald_mode = 1
#endif
          if (ewald_mode.eq.2) then
            write(ilog,*) 'PME is only supported while linking to ex&
 &ternal FFTW libraries. Switching to standard Ewald.'
          end if
        end if
!
!     the Fourier spacing for Ewald sums (interpreted differently depending on mode)
!
      else if (keyword(1:14) .eq. 'FMCSC_EWFSPAC ') then
        read (wo1(wt11:wt12),*) dum
        if (dum .le. 0.0d0) then
          ewfspac = 1.0
          write(ilog,*) 'Invalid spec. for FMCSC_EWFSPAC. Defaulting to ',ewfspac,'.'
        else
          ewfspac = dum
        end if
!
!
!     a user override for the Ewald parameter
!
      else if (keyword(1:12) .eq. 'FMCSC_EWPRM ') then
        read (wo1(wt11:wt12),*) dum
        if (dum .le. 0.0d0) then
          ewpm_pre = -0.1 ! negative value indicates that value needs to be determined by CAMPARI
          write(ilog,*) 'Invalid spec. for FMCSC_EWPRM. Defaulting to automatic determination.'
        else
          ewpm_pre = dum
        end if
!
      end if
!
    end if
!
!   RF spec.s
!     
    if (lrel_md.eq.3) then
!
!     mode for reaction-field treatment (1 is standard, 2 is GRF)
!
      if (keyword(1:13) .eq. 'FMCSC_RFMODE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.le.0).OR.(choice.gt.2)) then
          rf_mode = 1
          write(ilog,*) 'Invalid spec. for FMCSC_RFMODE. Defaulting &
 &to ',rf_mode,'.'
        else
          rf_mode = choice
        end if
!
      end if
!
    end if
!
    if ((use_LCTOR.EQV..true.).OR.(use_REMC.EQV..true.)) then
!
!     DISABLED FUNCTIONALITY 
!
      if (keyword(1:17) .eq. 'FMCSC_LCTPOTFILE ') then
        inquire(file=wo1(wt11:wt12),exist=exists)
        if (exists.EQV..true.) then
          lctpotfile = wo1(wt11:wt12)
          iunit = freeunit()
          open (unit=iunit,file=wo1(wt11:wt12),status='old')
          do i=1,par_LCTOR2(1)
            read(iunit,*) lct_val(i),&
 &                   (lct_pot(k,i),k=1,ntorlcs)
          end do
          par_LCTOR(2) = lct_val(1)
          par_LCTOR(1) = lct_val(2)-lct_val(1)
        else
          write(ilog,*) 'Warning. Cannot open spec.d input file (',&
 &wo1(wt11:wt12),') for FMCSC_LCTPOTFILE. Turning off LCT potential.'
          use_LCTOR = .false.
          scale_LCTOR = 0.0
        end if
!
      end if
!
    end if
!
!
!
  else
!
    write(ilog,*) 'Fatal. Called key_energies_l3(...) with level oth&
 &er than 3. This is most definitely a bug.'
    call fexit()
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! covers: torsn,dssps(some)
!
subroutine key_torsn(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 & rest,which,keyword)
!
  use iounit
  use torsn
  use keys
  use dssps
  use system
!
  implicit none
!
  character(MAXKEYLEN) wo1,wo2,wo3,rest
  character(MAXKWLEN) keyword
  integer which,choice,wt11,wt12,wt21,wt22,wt31,wt32,i
  RTYPE dum,dum4(4)
  logical exists
!
  if (which.eq.1) then
!
!
!   output frequency for system dihedrals (large file!), negative value switches mode
!
    if (keyword(1:13) .eq. 'FMCSC_TOROUT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.0) then
        write(ilog,*) 'Invalid spec. for FMCSC_TOROUT. Default.'
        torout = 10000
        toroutmode = 1
      else
        if (choice.lt.0) toroutmode = 2
        torout = abs(choice)
      end if
!
!   generic resolution for things controlled by FMCSC_ANGCALC 
!
    else if (keyword(1:14) .eq. 'FMCSC_ANGRES ') then
      read (wo1(wt11:wt12),*) fyres
      if ((fyres.lt.1.0).OR.(fyres.gt.180.0)) then
        fyres = 10.0
        write(ilog,*) 'Invalid spec. for FMCSC_ANGRES. Defaulting&
 & to ',fyres,'.'
      end if
!
!   calculation frequency for peptide-specific angular (dihedral) distribution functions
!
    else if (keyword(1:14) .eq. 'FMCSC_ANGCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        write(ilog,*) 'Invalid spec. for FMCSC_ANGCALC. Default.'
        angcalc = 50
      else
        angcalc = choice
      end if
!
!   calculation frequency for internal coordinate space variable distribution functions
!
    else if (keyword(1:14) .eq. 'FMCSC_INTCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        write(ilog,*) 'Invalid spec. for FMCSC_INTCALC. Default.'
        intcalc = 250
      else
        intcalc = choice
      end if
!
!   which classes of internal coordinate space variables to analyze
!
    else if (keyword(1:15) .eq. 'FMCSC_WHICHINT ') then
      read (wo1(wt11:wt12),*) dum4(1)
      call dum2int(dum4(1),choice)
      if (choice.eq.1) then
        do_ints(1) = .true.
      else
        do_ints(1) = .false.
      end if
      read (rest,*) (dum4(i),i=2,4)
      do i=2,4
        call dum2int(dum4(i),choice)
        if (choice.eq.1) then 
          do_ints(i) = .true.
        else
          do_ints(i) = .false.
        end if
      end do
!
!   residues for specific backbone distributions (ANGCALC)
!
    else if (keyword(1:14) .eq. 'FMCSC_RAMARES ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,nrspecrama)
      if ((nrspecrama.lt.0).OR.(nrspecrama.gt.MAXFYMAPS)) then
        nrspecrama = 0
        write(ilog,*) 'Invalid spec. for FMCSC_RAMARES. Defaulting t&
 &o ',nrspecrama,'.'
      else
        read(rest,*) (specrama(i),i=1,nrspecrama)
!        do i=1,nrspecrama
!          specrama(i) = specrama(i)
!        end do
      end if
!
!   molecules for specific backbone distributions (ANGCALC)
!
    else if (keyword(1:14) .eq. 'FMCSC_RAMAMOL ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,nrmolrama)
      if ((nrmolrama.lt.0).OR.(nrmolrama.gt.MAXFYMAPS)) then
        nrmolrama = 0
        write(ilog,*) 'Invalid spec. for FMCSC_RAMAMOL. Defaulting t&
 &o ',nrmolrama,'.'
      else
        read(rest,*) (molrama(i),i=1,nrmolrama)
!        do i=1,nrmolrama
!          molrama(i) = molrama(i)
!        end do
      end if
!
!   calculation frequency for segment distributions
!
    else if (keyword(1:14) .eq. 'FMCSC_SEGCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        segcalc = 50
        write(ilog,*) 'Invalid spec. for FMCSC_SEGCALC. Defaulting t&
 &o ',segcalc,'.'
      else
        segcalc = choice
      end if
!
!   calculation frequency for DSSP segment distributions
!
    else if (keyword(1:15) .eq. 'FMCSC_DSSPCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        dsspcalc = 1000
        write(ilog,*) 'Invalid spec. for FMCSC_DSSPCALC. Defaulting &
 &to ',dsspcalc,'.'
      else
        dsspcalc = choice
      end if
!
!   whether to provide classical (character-based) instantaneous DSSP-writeouts
!
    else if (keyword(1:15) .eq. 'FMCSC_INSTDSSP ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        inst_dssp = .true.
      end if
!
!   calculation frequency for covariance analysis
!
    else if (keyword(1:14) .eq. 'FMCSC_COVCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        covcalc = HUGE(covcalc)
        write(ilog,*) 'Invalid spec. for FMCSC_COVCALC. Defaulting t&
 &o ',covcalc,'.'
      else
        covcalc = choice
      end if
!
!   what type of covariance analysis to perform
!
    else if (keyword(1:14) .eq. 'FMCSC_COVMODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,covmode)
      if ((covmode.lt.1).OR.(covmode.gt.3)) then
        covmode = 1
        write(ilog,*) 'Invalid spec. for FMCSC_COVMODE. Defaulting t&
 &o ',covmode,'.'
      end if
!
!   calculation frequency for torsional linear combinations analysis, currently DISABLED FUNCTIONALITY
!
    else if (keyword(1:16) .eq. 'FMCSC_TORLCCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        torlccalc = 50
        write(ilog,*) 'Invalid spec. for FMCSC_TORLCCALC. Defaulting&
 & to ',torlccalc,'.'
      else
        torlccalc = choice
      end if
!
    end if
!
!
!
  else if (which.eq.2) then
!
    if (segcalc.le.nsim) then
!
!     backbone segments data input file
!
      if (keyword(1:16) .eq. 'FMCSC_BBSEGFILE ') then
        inquire(file=wo1(wt11:wt12),exist=exists)
        if (exists.EQV..false.) then
          write(ilog,*) 'Warning. Cannot open spec.d file with polypeptide secondary structure annotation (',&
 &wo1(wt11:wt12),') for FMCSC_BBSEGFILE (relates to FMCSC_SEGCALC). Disabling analysis.'
          segcalc = nsim + 1
        else
          bbsegfile = wo1(wt11:wt12)
          call bb_segment_reader()
        end if
      end if
!
    end if ! if segcalc.le.nsim

!
  else if (which.eq.3) then
!
  else
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! covers: system, grandensembles, shakeetal
!
subroutine key_system(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 & rest,which,keyword)
!
  use iounit
  use system
  use units
  use keys
  use grandensembles
  use forces
  use movesets
  use shakeetal
!
  implicit none
!
  character(MAXKEYLEN) wo1,wo2,wo3,rest
  character(MAXKWLEN) keyword
  integer which,choice,wt11,wt12,wt21,wt22,wt31,wt32,t1,t2
  RTYPE dum
!
  if (which.eq.1) then
!
!
!   whether run is simply analysis of pre-existing trajectory
!
    if (keyword(1:17) .eq. 'FMCSC_PDBANALYZE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        pdb_analyze = .true.
      else
        pdb_analyze = .false.
      end if
!
!   whether run uses united atoms
!
    else if (keyword(1:14) .eq. 'FMCSC_UAMODEL ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.le.2).AND.(choice.ge.0)) then
        ua_model = choice
      else
        ua_model = 0
        write(ilog,*) 'Invalid spec. for FMCSC_UAMODEL. De&
 &faulting to all-atom model (mode ',ua_model,').'
      end if
!
!   integer as to what kind of sampling methodology to use
!   1) MC (default), 2) MD, 3) LD, 4) BD (not yet), 5) mix MC/MD, 6) minimizer, 7) mix MC/LD
!
    else if (keyword(1:15) .eq. 'FMCSC_DYNAMICS ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (((choice.gt.1).AND.(choice.le.3)).OR.((choice.ge.5).AND.(choice.le.7))) then
        use_dyn = .true.
        dyn_mode = choice
      else if (choice.eq.1) then
        use_dyn = .false.
        dyn_mode = choice
      else
        dyn_mode = 1
        write(ilog,*) 'Invalid spec. for FMCSC_DYNAMICS. De&
 &faulting to MC simulation (mode ',dyn_mode,').'
      end if
!
!   choice of boundary condition
!
    else if (keyword(1:15) .eq. 'FMCSC_BOUNDARY ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.gt.4).OR.(choice.lt.1).OR.(choice.eq.2)) then
        bnd_type = 1
        write(ilog,*) 'Invalid spec. for FMCSC_BOUNDARY. De&
 &faulting to periodic boundaries (mode ',bnd_type,').'
      else
        bnd_type = choice
      end if
!
!   equilibration steps
!
    else if (keyword(1:12) .eq. 'FMCSC_EQUIL ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,nequil)
      if (nequil .lt. 0) then
        nequil = 10000
        write(ilog,*) 'Invalid spec. for FMCSC_EQUIL. Defaulting to &
 &',nequil,'.'
      end if
!
!   overall number of steps
!
    else if (keyword(1:14) .eq. 'FMCSC_NRSTEPS ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,nsim)
      if (nsim .le. 0) then
        nsim = 100000
        write(ilog,*) 'Invalid spec. for FMCSC_NRSTEPS. Defaulting t&
 &o ',nsim,'.'
      end if
!
!   dynamics timestep in ps
!
    else if (keyword(1:15) .eq. 'FMCSC_TIMESTEP ') then
      read (wo1(wt11:wt12),*) dyn_dt
      if (dyn_dt .le. 0.0d0)  then
        dyn_dt = 0.002 ! 2fs
        write(ilog,*) 'Invalid spec. for FMCSC_TIMESTEP. Defaulting &
 &to ',dyn_dt,'.'
      end if
!
!   working ensemble (1: NVT, 2: NVE, 3: NPT(not yet), 4: NPE (not yet), 5:uVT, 6:(du/N)VT
!
    else if (keyword(1:15) .eq. 'FMCSC_ENSEMBLE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,ens%flag)
      if ((ens%flag.lt.1).OR.(ens%flag.gt.6)) then
        ens%flag = 1
        write(ilog,*) 'Invalid spec. for FMCSC_ENSEMBLE. Defaulting &
 &to NVT (',ens%flag,').'
      end if
!
!   working temperature (thermostat target in NVT-MD, initial T in NVE-MD)
!
    else if (keyword(1:11) .eq. 'FMCSC_TEMP ') then
      read (wo1(wt11:wt12),*) kelvin
      if (kelvin .le. 0.0d0)  then
        kelvin = 298.0d0
        write(ilog,*) 'Invalid spec. for FMCSC_TEMP. Defaulting to '&
 &,kelvin,'.'
      end if
      invtemp = 1.0d0 /(kelvin*gasconst)
!
!   working external pressure for XPX-ensembles, currently DISABLED FUNCTIONALITY
!
    else if (keyword(1:15) .eq. 'FMCSC_PRESSURE ') then
      read (wo1(wt11:wt12),*) extpress
      if (extpress .le. 0.0d0)  then
        extpress = 1.0d0
        write(ilog,*) 'Invalid spec. for FMCSC_PRESSURE. Defaulting &
 &to ',extpress,'.'
      end if
!
!   calculation frequency for molecule type number histograms
!
    else if (keyword(1:14) .eq. 'FMCSC_NUMCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        write(ilog,*) 'Invalid spec. for FMCSC_NUMCALC. Default.'
        particlenumcalc = 100
      else
        particlenumcalc = choice
      end if
!
!   chemical potential input for GCMC
!        
    else if (keyword(1:23) .eq. 'FMCSC_PARTICLEFLUCFILE ') then
      particleflucfile = wo1(wt11:wt12)
!
!   whether to write a summary of the fluctuating molecule types in the system
!
    else if (keyword(1:18) .eq. 'FMCSC_GRANDREPORT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        sgc_report = .true.
      else
        sgc_report = .false.
      end if
!
!   whether to use absolute chemical potentials or excess mu + expected number
!
    else if (keyword(1:16) .eq. 'FMCSC_GRANDMODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.lt.1).OR.(choice.gt.2)) then
        write(ilog,*) 'Invalid spec. for FMCSC_GRANDMODE. Defaulting to the use of &
 &excess chemical potentials and expected particle numbers.'
        gc_mode = 2
      else
        gc_mode = choice
      end if
!
!   Base filename used for trajectories and structural files ...
!
    else if (keyword(1:15) .eq. 'FMCSC_BASENAME ') then
      basename = wo1(wt11:wt12)
      call strlims(basename,t1,t2)
      bleng = t2-t1+1
!
!   whether to restart run from the latest restart point
!
    else if (keyword(1:14) .eq. 'FMCSC_RESTART ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        do_restart = .true.
      else
        do_restart = .false.
      end if
!
!   whether to allow skipping of reading dynamical variables from restart
!   file, and to re-assign them instead (this only works for restart files
!   generated by MC runs!)
!
    else if (keyword(1:16) .eq. 'FMCSC_RST_MC2MD ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        mc_compat_flag = 1
      end if
!
!   whether to fudge masses in inertial dynamics 
!
    else if (keyword(1:16) .eq. 'FMCSC_FUDGE_DYN ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        fudge_mass = .true.
      end if
!
!   choice of specific TMD integrator variant
!
    else if (keyword(1:21) .eq. 'FMCSC_TMD_INTEGRATOR ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.lt.1).OR.(choice.gt.2)) then
        dyn_integrator = 2
        write(ilog,*) 'Warning. Invalid spec. for FMCSC_TMD_INTEGRATOR. Def&
 &aulting to ',dyn_integrator,'.'
      else
        dyn_integrator = choice
      end if
!
!   incremental velocity update option for TMD integrator variants
!
    else if (keyword(1:17) .eq. 'FMCSC_TMD_INT2UP ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.lt.0) then
        dyn_integrator_ops(1) = 0
        write(ilog,*) 'Warning. Invalid spec. for FMCSC_TMD_INT2UP. Def&
 &aulting to ',dyn_integrator_ops(1),'.'
      else
        dyn_integrator_ops(1) = choice
      end if
!
!    flag to decide what to do with natively unsupported d.o.f.s in dynamics
!    0) sample none, 1) sample UNK, 2) sample unsupported in native, 3) sample all
!
    else if (keyword(1:18) .eq. 'FMCSC_TMD_UNKMODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.lt.0).OR.(choice.gt.3)) then
        dyn_integrator_ops(10) = 1
        write(ilog,*) 'Warning. Invalid spec. for FMCSC_TMD_UNKMODE. Def&
 &aulting to ',dyn_integrator_ops(10),'.'
      else
        dyn_integrator_ops(10) = choice
      end if

!
    end if
!
!
!
  else if (which.eq.2) then
!
!     
!   choice of sampling methodology
!
    if (use_dyn.EQV..true.) then
!
!     what kind of system motion removal to do
!
      if (keyword(1:13) .eq. 'FMCSC_SYSFRZ ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.3)) then
          ens%sysfrz = 1
          write(ilog,*) 'Invalid spec. for FMCSC_SYSFRZ. De&
 &faulting to ',ens%sysfrz,'.'
        else
          ens%sysfrz = choice
        end if
!
!     thermostat type
!
      else if (keyword(1:12) .eq. 'FMCSC_TSTAT ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.4)) then
          tstat%flag = 4
          write(ilog,*) 'Invalid spec. for FMCSC_TSTAT. Def&
 &aulting to ',tstat%flag,' (Bussi-Parrinello velocity rescaling).'
        else
          tstat%flag = choice
        end if
!
!     thermostat decay constant
!
      else if (keyword(1:16) .eq. 'FMCSC_TSTAT_TAU ') then
        read (wo1(wt11:wt12),*) tstat%params(1)
        if (tstat%params(1).lt.dyn_dt)  then
          tstat%params(1) = dyn_dt
          write(ilog,*) 'Invalid spec. for FMCSC_TSTAT_TAU. Defaulti&
 &ng to ',tstat%params(1),'.'
        end if
!
!     input file for T-coupling group requests
!
      else if (keyword(1:17) .eq. 'FMCSC_TSTAT_FILE ') then
        tstat%fnam = wo1(wt11:wt12)
!
!     manostat type, currently DISABLED FUNCTIONALITY
!
      else if (keyword(1:12) .eq. 'FMCSC_PSTAT ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.3)) then
          pstat%flag = 1
          write(ilog,*) 'Invalid spec. for FMCSC_PSTAT. Def&
 &aulting to ',pstat%flag,'.'
        else
          pstat%flag = choice
        end if
!
!     manostat mass, currently DISABLED FUNCTIONALITY
!
      else if (keyword(1:17) .eq. 'FMCSC_PSTAT_MASS ') then
        read (wo1(wt11:wt12),*) pstat%params(1)
        if (pstat%params(1).le.0.0)  then
          pstat%params(1) = 10.0
          write(ilog,*) 'Invalid spec. for FMCSC_PSTAT_MASS. Default&
 &ing to ',pstat%params(1),'.'
        end if
!
!     generalized friction parameter (LD)
!
      else if (keyword(1:15) .eq. 'FMCSC_FRICTION ') then
        read (wo1(wt11:wt12),*) fric_ga
        if (fric_ga.le.0.0)  then
          fric_ga = 1.0
          write(ilog,*) 'Invalid spec. for FMCSC_FRICTION. Defaultin&
 &g to ',fric_ga,'.'
        end if
!
!     general report flag for dynamics-related things
!
      else if (keyword(1:16) .eq. 'FMCSC_DYNREPORT ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.1) then
          dyn_report = .true.
        else
          dyn_report = .false.
        end if
!
!     whether to perform a numerical check of the gradients
!
      else if (keyword(1:16) .eq. 'FMCSC_CHECKGRAD ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.1) then
          grad_check = .true.
        else
          grad_check = .false.
        end if
!
!     cycle lengths for hybrid methods
!
      else if (keyword(1:19) .eq. 'FMCSC_CYCLE_MC_MIN ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.ge.nsim)) then
          min_mccyclen = 200
          write(ilog,*) 'Invalid spec. for FMCSC_CYCLE_MC_MIN. Default&
 &ing to ',min_mccyclen,'.'
        else
          min_mccyclen = choice
        end if
      else if (keyword(1:19) .eq. 'FMCSC_CYCLE_MC_MAX ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.ge.nsim)) then
          max_mccyclen = 500
          write(ilog,*) 'Invalid spec. for FMCSC_CYCLE_MC_MAX. Default&
 &ing to ',max_mccyclen,'.'
        else
          max_mccyclen = choice
        end if
      else if (keyword(1:20) .eq. 'FMCSC_CYCLE_DYN_MIN ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.ge.nsim)) then
          min_dyncyclen = 500
          write(ilog,*) 'Invalid spec. for FMCSC_CYCLE_DYN_MIN. Default&
 &ing to ',min_dyncyclen,'.'
        else
          min_dyncyclen = choice
        end if
      else if (keyword(1:20) .eq. 'FMCSC_CYCLE_DYN_MAX ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.ge.nsim)) then
          max_dyncyclen = 1000
          write(ilog,*) 'Invalid spec. for FMCSC_CYCLE_DYN_MAX. Default&
 &ing to ',max_dyncyclen,'.'
        else
          max_dyncyclen = choice
        end if
      else if (keyword(1:21) .eq. 'FMCSC_CYCLE_MC_FIRST ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.ge.nsim)) then
          first_mccyclen = 50000
          write(ilog,*) 'Invalid spec. for FMCSC_CYCLE_MC_FIRST. Default&
 &ing to ',first_mccyclen,'.'
        else
          first_mccyclen = choice
        end if
!
      end if
!
    end if ! if use_dyn.EQV..true.
!
!   what to consider as the degrees of freedom
!
    if (keyword(1:14) .eq. 'FMCSC_CARTINT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.lt.1).OR.(choice.gt.2)) then
        fycxyz = 1
        write(ilog,*) 'Invalid spec. for FMCSC_CARTINT. Def&
 &aulting to mixed rigid-body / torsional space simulation (mode ',fycxyz,').'
      else
        if (choice.eq.2) then
          if ((dyn_mode.ne.2).AND.(dyn_mode.ne.3).AND.(dyn_mode.ne.5).AND.(dyn_mode.ne.7).AND.(dyn_mode.ne.6)) then
            write(ilog,*) 'Warning. Cartesian space sampling is not supported in pure Monte Carlo &
 &runs. Switching to mixed rigid-body / torsional space sampling.'
            choice = 1
          end if
        end if
        fycxyz = choice
      end if
    end if
! 
!   choice of simulation container shape (as a fxn of boundary condition)
!
    if (keyword(1:12) .eq. 'FMCSC_SHAPE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
!     in case of periodic boundary conditions
      if (bnd_type.eq.1) then
        if (choice.ne.1) then ! .AND.(choice.ne.3)) then
          bnd_shape = 1
          write(ilog,*) 'Invalid spec. for FMCSC_SHAPE with&
 & boundary condition of type ',bnd_type,'. Defaulting to a cubic bo&
 &x (mode ',bnd_shape,').'
        else
          bnd_shape = choice
        end if
      else if (bnd_type.eq.2) then
        if (choice.ne.2) then
          bnd_shape = 2
          write(ilog,*) 'Invalid spec. for FMCSC_SHAPE with&
 & boundary condition of type ',bnd_type,'. Defaulting to a spherica&
 &l system (mode ',bnd_shape,').'
        else
          bnd_shape = choice
        end if
      else if (bnd_type.eq.3) then
        if (choice.gt.2) then
          bnd_shape = 2
          write(ilog,*) 'Invalid spec. for FMCSC_SHAPE with&
 & boundary condition of type ',bnd_type,'. Defaulting to a spherica&
 &l system (mode ',bnd_shape,').'
        else
          bnd_shape = choice
        end if
      else if (bnd_type.eq.4) then
        if (choice.gt.2) then
          bnd_shape = 2
          write(ilog,*) 'Invalid spec. for FMCSC_SHAPE with&
 & boundary condition of type ',bnd_type,'. Defaulting to a spherica&
 &l system (mode ',bnd_shape,').'
        else
          bnd_shape = choice
        end if
      end if
!
    end if
!
! spec's for specific boundary condition
!
    if ((bnd_type.eq.3).OR.(bnd_type.eq.4)) then
!
!   force constant for soft wall BC
!
      if (keyword(1:15) .eq. 'FMCSC_SOFTWALL ') then
        read (wo1(wt11:wt12),*) bnd_params(7)
        if (bnd_params(7).le.0.0d0)  then
          bnd_params(7) = 1.0
            write(ilog,*) 'Invalid spec. for FMCSC_SOFTWALL. Default&
 &ing to ',bnd_params(7),' kcal/(mol*A*A).'
        end if
!
      end if
!
    end if
!
!
!
  else if (which.eq.3) then
!
!
!  third order spec.s for simulation system
!
!   choice of origin for simulation system
!
    if (keyword(1:13) .eq. 'FMCSC_ORIGIN ') then
!     in case of a cubic system
      if (bnd_shape.eq.1) then
        read (wo1(wt11:wt12),*) bnd_params(4)
        read (wo2(wt21:wt22),*) bnd_params(5)
        read (wo3(wt31:wt32),*) bnd_params(6)
      else if ((bnd_shape.eq.2).OR.(bnd_shape.eq.3)) then
        read (wo1(wt11:wt12),*) bnd_params(1)
        read (wo2(wt21:wt22),*) bnd_params(2)
        read (wo3(wt31:wt32),*) bnd_params(3)
      end if
!
!    as well as dimensions 
!
    else if (keyword(1:11) .eq. 'FMCSC_SIZE ') then
!     in case of a cubic system
      if (bnd_shape.eq.1) then
        read (wo1(wt11:wt12),*) bnd_params(1)
        read (wo2(wt21:wt22),*) bnd_params(2)
        read (wo3(wt31:wt32),*) bnd_params(3)
      else if (bnd_shape.eq.2) then
        read (wo1(wt11:wt12),*) bnd_params(4)
        bnd_params(5) = bnd_params(4)**2
      else if (bnd_shape.eq.3) then
        read (wo1(wt11:wt12),*) bnd_params(4)
        read (wo2(wt21:wt22),*) bnd_params(6)
        bnd_params(5) = bnd_params(4)**2
      end if
!
    end if
!
!   SHAKE parameters for Cartesian integrator
!
    if ((fycxyz.eq.2).AND.(use_dyn.EQV..true.)) then 
!
!     relative tolerance for bond lengths (otherwise unit-dependent)
!
      if (keyword(1:15) .eq. 'FMCSC_SHAKETOL ') then
        read (wo1(wt11:wt12),*) shake_tol
        if ((shake_tol.le.0.0d0).OR.(shake_tol.ge.1.0e-2))  then
          shake_tol = 1.0e-4
            write(ilog,*) 'Invalid spec. for FMCSC_SHAKETOL. Default&
 &ing to ',shake_tol,'.'
        end if

!     absolute tolerance for cosine of bond angles (normalized)
!
      else if (keyword(1:16) .eq. 'FMCSC_SHAKEATOL ') then
        read (wo1(wt11:wt12),*) shake_atol
        if ((shake_atol.le.0.0d0).OR.(shake_atol.ge.1.0e-2))  then
          shake_atol = 1.0e-4
            write(ilog,*) 'Invalid spec. for FMCSC_SHAKEATOL. Default&
 &ing to ',shake_atol,'.'
        end if
!
!     which bonds/angles to constrain
!
      else if (keyword(1:15) .eq. 'FMCSC_SHAKESET ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.6)) then
          cart_cons_mode = 1
          write(ilog,*) 'Invalid spec. for FMCSC_SHAKESET. Default&
 &ing to the use of no constraints (option ',cart_cons_mode,').'
        else
          cart_cons_mode = choice
        end if
!
!     how to obtain the required input lengths/angles
!
      else if (keyword(1:16) .eq. 'FMCSC_SHAKEFROM ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.3)) then
          cart_cons_source = 1
          write(ilog,*) 'Invalid spec. for FMCSC_SHAKEFROM. Default&
 &ing to the use of database information (CAMPARI default geometry, i.e., option ',cart_cons_source,').'
        else
          cart_cons_source = choice
        end if
!
!     add settle specifically for water molecules?
!
      else if (keyword(1:16) .eq. 'FMCSC_SETTLEH2O ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.1) then
          add_settle = .true.
        else
          add_settle = .false.
        end if
!
!     effective tolerance (accuracy) setting for LINCS
!
      else if (keyword(1:17) .eq. 'FMCSC_LINCSORDER ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.2).OR.(choice.gt.16)) then
          lincs_order = 8
          write(ilog,*) 'Invalid spec. for FMCSC_LINCSORDER. Default&
 &ing to ',lincs_order,'.'
        else
          lincs_order = choice
        end if
!
!     which method to use for holonomic constraints
!
      else if (keyword(1:18) .eq. 'FMCSC_SHAKEMETHOD ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
#ifdef LINK_LAPACK
        if ((choice.lt.1).OR.(choice.gt.4)) then
#else
        if ((choice.ne.1).AND.(choice.ne.4)) then
#endif
          cart_cons_method = 1
          write(ilog,*) 'Invalid spec. for FMCSC_SHAKEMETHOD. Default&
 &ing to standard SHAKE (method ',cart_cons_method,').'
        else
          cart_cons_method = choice
        end if
!
!     how many iterations to accept before termination in SHAKE
!
      else if (keyword(1:19) .eq. 'FMCSC_SHAKEMAXITER ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.100).OR.(choice.gt.10000000)) then
          shake_maxiter = 1000
          write(ilog,*) 'Invalid spec. for FMCSC_SHAKEMAXITER. Default&
 &ing to a limit of ',shake_maxiter,' iterations.'
        else
          shake_maxiter = choice
        end if
!
!     optional file input as to which atom-atom distances to constrain
!
      else if (keyword(1:16) .eq. 'FMCSC_SHAKEFILE ') then
        cart_cons_file = wo1(wt11:wt12)
!
      end if

!
    end if ! if Cartesian dynamics
!
!
  else
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! covers: movesets, grids, ionize, minimize, ujglobals, wl
!
subroutine key_movesets(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 & rest,which,keyword)
!
  use iounit
  use movesets
  use keys
  use math
  use ionize
  use grids
  use mini
  use ujglobals
  use wl
  use system
!
  implicit none
!
  character(MAXKEYLEN) wo1,wo2,wo3,rest
  character(MAXKWLEN) keyword
  character(MAXKEYLEN+20) str2
  integer which,choice,wt11,wt12,wt21,wt22,wt31,wt32
  RTYPE dum
  logical exists
!
  if (which.eq.1) then
!
!
!   pick the acceptance criterion to use
!
    if (keyword(1:17) .eq. 'FMCSC_MC_ACCEPT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.le.0).OR.(choice.gt.3)) then
        mc_acc_crit = 1
        write(ilog,*) 'Invalid spec. for FMCSC_MC_ACCEPT. Defaulting to ',mc_acc_crit,' (standard Metropolis).'
        do_wanglandau = .false.
      else
        mc_acc_crit = choice
        if (mc_acc_crit.eq.3) do_wanglandau = .true.
      end if
!
!   which WL mode to use (1 = energy, always 1D, 2 = reaction coordinate (Rg), incl 2D, 3 = energy PMF, incl 2D)
!
    else if (keyword(1:14) .eq. 'FMCSC_WL_MODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.le.0).OR.(choice.gt.3)) then
        wld%wl_mode = 1
        write(ilog,*) 'Invalid spec. for FMCSC_WL_MODE. Defaulting to ',wld%wl_mode,'.'
      else
        wld%wl_mode = choice
      end if
!
!   the very general input file to tune the picking probabilities for move types
!
    else if (keyword(1:14) .eq. 'FMCSC_PSWFILE ') then
      prefsamplingfile = wo1(wt11:wt12)
!
!   whether to write a summary of the adjusted sampling weights in the system
!
    else if (keyword(1:16) .eq. 'FMCSC_PSWREPORT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        psw_report = .true.
      else
        psw_report = .false.
      end if
!
!   choice of pivot biased technique: none, grids or HCs
!
    else if (keyword(1:16) .eq. 'FMCSC_PIVOTMODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        use_stericgrids = .false.
      else if (choice.eq.2) then
        use_stericgrids = .true.
      else
        write(ilog,*) 'Invalid spec. for FMCSC_PIVOTMODE. Defaulting&
 & to simple (blind) mode.'
        use_stericgrids = .false.
      end if
!
!   choice of frequency of pivot moves which totally randomize a pair of bb-angles
!   (vs. those that slightly perturb a pair of bb-angles)
!
    else if (keyword(1:18) .eq. 'FMCSC_PIVOTRDFREQ ') then
      read (wo1(wt11:wt12),*) pivot_randfreq
      if ((pivot_randfreq.lt.0.0d0).OR.(pivot_randfreq.gt.1.0))&
 & then
        pivot_randfreq = 0.1
        write(ilog,*) 'Invalid spec. for FMCSC_PIVOTRDFREQ. Defaulti&
 &ng to ',pivot_randfreq,'.'
      end if
!
!   step size (in degrees) for the local brand of pivot moves
!
    else if (keyword(1:18) .eq. 'FMCSC_PIVOTSTEPSZ ') then
      read (wo1(wt11:wt12),*) pivot_stepsz
      if (pivot_stepsz.le.0.0d0) then
        pivot_stepsz = 2.0
        write(ilog,*) 'Invalid spec. for FMCSC_PIVOTSTEPSZ. Defaulti&
 &ng to ',pivot_stepsz,'.'
      end if
!
!   choice of frequency of omega moves amongst total single residue (pivot) moves
!
    else if (keyword(1:16) .eq. 'FMCSC_OMEGAFREQ ') then
      read (wo1(wt11:wt12),*) omegafreq
      if ((omegafreq.lt.0.0d0).OR.(omegafreq.gt.1.0))&
 & then
        omegafreq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_OMEGAFREQ. Defaulting&
 & to ',omegafreq,'.'
      end if
!
!   choice of frequency of omega moves which totally randomize an individual omega angle
!
    else if (keyword(1:18) .eq. 'FMCSC_OMEGARDFREQ ') then
      read (wo1(wt11:wt12),*) omega_randfreq
      if ((omega_randfreq.lt.0.0d0).OR.(omega_randfreq.gt.1.0))&
 & then
        omega_randfreq = 0.4
        write(ilog,*) 'Invalid spec. for FMCSC_OMEGARDFREQ. Defaulti&
 &ng to ',omega_randfreq,'.'
      end if
!
!   step size (in degrees) for the local brand of omega moves
!
    else if (keyword(1:18) .eq. 'FMCSC_OMEGASTEPSZ ') then
      read (wo1(wt11:wt12),*) omega_stepsz
      if (omega_stepsz.le.0.0d0) then
        omega_stepsz = 2.0
        write(ilog,*) 'Invalid spec. for FMCSC_OMEGASTEPSZ. Defaulti&
 &ng to ',omega_stepsz,'.'
      end if
!
!   how to deal with chain building direction: normally the N-terminus stays put while the
!   C-terminus flaps around. this keyword allow changing that default behavior
!   1) always N-aligned (def.) 2) always C-aligned 3) shorter tail always aligned
!   4) shorter tail is preferably aligned
!
    else if (keyword(1:12) .eq. 'FMCSC_ALIGN ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.gt.4).OR.(choice.lt.1)) then
        align_NC = 3
        write(ilog,*) 'Invalid spec. for FMCSC_ALIGN. Defaulting to &
 & mandatory N-terminal alignemnt (mode ',align_NC,').'
      else
        align_NC = choice
      end if
!
!   coupling of backbone and sidechain torsional degrees of freedom
!
    else if (keyword(1:13) .eq. 'FMCSC_COUPLE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        use_globmoves = .true.
      end if
!
!   how many chi angles to sample in an elementary chi-move (note
!   that large values are meaningless since the max. is restricted by
!   the number of chis in an individual sidechain)
!
    else if (keyword(1:12) .eq. 'FMCSC_NRCHI ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.lt.1) then
        nrchis = 2
        write(ilog,*) 'Invalid spec. for FMCSC_NRCHI. Defaulting to &
 &',nrchis,'.'
      else
        nrchis = choice
      end if
!
!   fraction of total torsional moves to be sidechain moves
!
    else if (keyword(1:14) .eq. 'FMCSC_CHIFREQ ') then
      read (wo1(wt11:wt12),*) chifreq
      if ((chifreq.lt.0.0).OR.(chifreq.gt.1.0)) then
        chifreq = 0.5
        write(ilog,*) 'Invalid spec. for FMCSC_CHIFREQ. Defaulting t&
 &o ',chifreq,'.'
      end if
!
!   choice of frequency of chi-moves which totally randomize a number of chi-angles
!   (vs. those that slightly perturb a number of chi-angles)
!
    else if (keyword(1:16) .eq. 'FMCSC_CHIRDFREQ ') then
      read (wo1(wt11:wt12),*) chi_randfreq
      if ((chi_randfreq.lt.0.0d0).OR.(chi_randfreq.gt.1.0))&
 & then
        chi_randfreq = 0.1
        write(ilog,*) 'Invalid spec. for FMCSC_CHIRDFREQ. Defaulting&
 & to ',chi_randfreq,'.'
      end if
!
!   step size (in degrees) for the local brand of sidechain moves
!
    else if (keyword(1:16) .eq. 'FMCSC_CHISTEPSZ ') then
      read (wo1(wt11:wt12),*) chi_stepsz
      if (chi_stepsz.le.0.0d0) then
        chi_stepsz = 10.0
        write(ilog,*) 'Invalid spec. for FMCSC_CHISTEPSZ. Defaulting&
 & to ',chi_stepsz,'.'
      end if
!
!   choice of frequency for polypeptide concerted rotation moves
!
    else if (keyword(1:13) .eq. 'FMCSC_CRFREQ ') then
      read (wo1(wt11:wt12),*) crfreq
      if ((crfreq.lt.0.0).OR.(crfreq.gt.1.0)) then
        crfreq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_CRFREQ. Defaulting to ',crfreq,'.' 
      end if
!
!   amongst fraction of CR moves, fraction to include bond-angle degrees of freedom
!
    else if (keyword(1:16) .eq. 'FMCSC_ANGCRFREQ ') then
      read (wo1(wt11:wt12),*) angcrfreq
      if ((angcrfreq.lt.0.0).OR.(angcrfreq.gt.1.0)) then
        angcrfreq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_ANGCRFREQ. Defaulting to ',angcrfreq,'.'
      end if
!
!   amongst fraction of non-angular CR backbone moves, fraction to be UJ-tor-CR
!
    else if (keyword(1:16) .eq. 'FMCSC_TORCRFREQ ') then
      read (wo1(wt11:wt12),*) torcrfreq
      if ((torcrfreq.lt.0.0).OR.(torcrfreq.gt.1.0)) then
        torcrfreq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_TORCRFREQ. Defaulting to ',torcrfreq,'.'
      end if
!
!   whether to use an algorithm rigorously satisfying detailed balance in theory (1) or
!   one that is computationally much more efficient (2)
!
    else if (keyword(1:16) .eq. 'FMCSC_TORCRMODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.lt.1).OR.(choice.gt.2)) then
        torcrmode = 1
        write(ilog,*) 'Invalid spec. for FMCSC_TORCRMODE. Defaulting to ',torcrmode,'.'
      else
        torcrmode = choice
      end if
!
!   amongst fraction of Dinner-Ulmschneider CR backbone moves, fraction to include omega 
!
    else if (keyword(1:17) .eq. 'FMCSC_TORCROFREQ ') then
      read (wo1(wt11:wt12),*) torcrfreq_omega
      if ((torcrfreq_omega.lt.0.0).OR.(torcrfreq_omega.gt.1.0)) then
        torcrfreq_omega = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_TORCROFREQ. Defaulting to ',torcrfreq_omega,'.'
      end if
!
!   choice of frequency for polypeptide ring pucker moves
!
    else if (keyword(1:14) .eq. 'FMCSC_PKRFREQ ') then
      read (wo1(wt11:wt12),*) puckerfreq
      if ((puckerfreq.lt.0.0).OR.(puckerfreq.gt.1.0)) then
        puckerfreq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_PKRFREQ.  Defaulting to ',puckerfreq,'.'
      end if
!
!   keyword for ratio of pucker moves that are inversion moves
!  
    else if (keyword(1:16) .eq. 'FMCSC_PKRRDFREQ ') then
      read (wo1(wt11:wt12),*) puckerrdfreq
      if ((puckerrdfreq.lt.0.0).OR.(puckerrdfreq.gt.1.0)) then
        puckerrdfreq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_PKRRDFREQ.  Defaulting to ',puckerrdfreq,'.'
      end if
!
!   choice of frequency for nucleic acid backbone moves
!
    else if (keyword(1:14) .eq. 'FMCSC_NUCFREQ ') then
      read (wo1(wt11:wt12),*) nucfreq
      if ((nucfreq.lt.0.0).OR.(nucfreq.gt.1.0)) then
        nucfreq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_NUCFREQ. Defaulting to ',nucfreq,'.'
      end if
!
!   how many generic backbone angles to sample in a generic backbone-move (note
!   that large values are meaningless since the max. is restricted by
!   the number of backbone angles in an individual nucleotide residue)
!
    else if (keyword(1:12) .eq. 'FMCSC_NRNUC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.lt.1) then
        nrnucim = 2
        write(ilog,*) 'Invalid spec. for FMCSC_NRNUC. Defaulting to ',nrchis,'.'
      else
        nrnucim = choice
      end if
!
!   choice of frequency of generic backbone-moves which totally randomize a number of backbone angles
!
    else if (keyword(1:16) .eq. 'FMCSC_NUCRDFREQ ') then
      read (wo1(wt11:wt12),*) nuc_randfreq
      if ((nuc_randfreq.lt.0.0d0).OR.(nuc_randfreq.gt.1.0)) then
        nuc_randfreq = 0.1
        write(ilog,*) 'Invalid spec. for FMCSC_NUCRDFREQ. Defaulting to ',nuc_randfreq,'.'
      end if
!
!   step size (in degrees) for the local brand of generic backbone moves for nucleotides
!
    else if (keyword(1:16) .eq. 'FMCSC_NUCSTEPSZ ') then
      read (wo1(wt11:wt12),*) nuc_stepsz
      if (nuc_stepsz.le.0.0d0) then
        nuc_stepsz = 10.0
        write(ilog,*) 'Invalid spec. for FMCSC_NUCSTEPSZ. Defaulting to ',nuc_stepsz,'.'
      end if
!
!   amongst fraction of nucleic acid moves, what is fraction of concerted rotation (n.a.) moves
!
    else if (keyword(1:16) .eq. 'FMCSC_NUCCRFREQ ') then
      read (wo1(wt11:wt12),*) nuccrfreq
      if ((nuccrfreq.lt.0.0).OR.(nuccrfreq.gt.1.0)) then
        nuccrfreq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_NUCCRFREQ. Defaulting to ',nuccrfreq,'.'
      end if
!
!   amongst fraction of non-CR nucleic acid moves, what is fraction of sugar pucker moves
!
    else if (keyword(1:16) .eq. 'FMCSC_SUGARFREQ ') then
      read (wo1(wt11:wt12),*) nucpuckfreq
      if ((nucpuckfreq.lt.0.0).OR.(nucpuckfreq.gt.1.0)) then
        nucpuckfreq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_SUGARFREQ. Defaulting to ',nucpuckfreq,'.'
      end if
!
!   amongst fraction fraction of sugar pucker moves, what fraction of moves is of ring inversion type
!
    else if (keyword(1:18) .eq. 'FMCSC_SUGARRDFREQ ') then
      read (wo1(wt11:wt12),*) nucpuckrdfreq
      if ((nucpuckrdfreq.lt.0.0).OR.(nucpuckrdfreq.gt.1.0)) then
        nucpuckrdfreq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_SUGARRDFREQ. Defaulting to ',nucpuckrdfreq,'.'
      end if
!
!   step size parameter for torsional changes in random perturbations of pucker state (both peptide and n.a.)
!
    else if (keyword(1:20) .eq. 'FMCSC_PUCKERSTEP_DI ') then
      read (wo1(wt11:wt12),*) pucker_distp
      if ((pucker_distp.lt.0.0).OR.(pucker_distp.gt.180.0)) then
        pucker_distp = 5.0
        write(ilog,*) 'Invalid spec. for FMCSC_PUCKERSTEP_DI. Defaulting to ',pucker_distp,'.'
      end if
!
!   step size parameter for angle changes in random perturbations of pucker state (both peptide and n.a.)
!
    else if (keyword(1:20) .eq. 'FMCSC_PUCKERSTEP_AN ') then
      read (wo1(wt11:wt12),*) pucker_anstp
      if ((pucker_anstp.lt.0.0).OR.(pucker_anstp.gt.20.0)) then
        pucker_anstp = 2.0
        write(ilog,*) 'Invalid spec. for FMCSC_PUCKERSTEP_AN. Defaulting to ',pucker_anstp,'.'
      end if
!
!   choice of frequency for single dihedral moves on all d.o.f.s
!
    else if (keyword(1:16) .eq. 'FMCSC_OTHERFREQ ') then
      read (wo1(wt11:wt12),*) otherfreq
      if ((otherfreq.lt.0.0).OR.(otherfreq.gt.1.0)) then
        otherfreq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_OTHERFREQ. Defaulting to ',otherfreq,'.'
      end if
!
!   the corresponding full randomization freq.
!
    else if (keyword(1:18) .eq. 'FMCSC_OTHERRDFREQ ') then
      read (wo1(wt11:wt12),*) other_randfreq
      if ((other_randfreq.lt.0.0).OR.(other_randfreq.gt.1.0)) then
        other_randfreq = 0.2
        write(ilog,*) 'Invalid spec. for FMCSC_OTHERRDFREQ. Defaulting to ',other_randfreq,'.'
      end if
!
!   the corresponding step size parameter
!
    else if (keyword(1:18) .eq. 'FMCSC_OTHERSTEPSZ ') then
      read (wo1(wt11:wt12),*) other_stepsz
      if ((other_stepsz.lt.0.0).OR.(other_stepsz.gt.180.0)) then
        other_stepsz = 20.0
        write(ilog,*) 'Invalid spec. for FMCSC_OTHERSTEPSZ. Defaulting to ',other_stepsz,'.'
      end if
!
!   amongst single dihedral moves on all d.o.f.s, what is fraction operating on d.o.f.s in unsupported residues
!
    else if (keyword(1:19) .eq. 'FMCSC_OTHERUNKFREQ ') then
      read (wo1(wt11:wt12),*) other_unkfreq
      if ((other_unkfreq.lt.0.0).OR.(other_unkfreq.gt.1.0)) then
        other_unkfreq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_OTHERUNKFREQ. Defaulting to ',other_unkfreq,'.'
      end if
!
!   amongst single dihedral moves on all d.o.f.s not operating on d.o.f.s in unsupported residues, what is fraction
!   operating on native and supported d.o.f.s
!
    else if (keyword(1:19) .eq. 'FMCSC_OTHERNATFREQ ') then
      read (wo1(wt11:wt12),*) other_natfreq
      if ((other_natfreq.lt.0.0).OR.(other_natfreq.gt.1.0)) then
        other_natfreq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_OTHERNATFREQ. Defaulting to ',other_natfreq,'.'
      end if
!
!   choice of frequency for rigid-body moves
!
    else if (keyword(1:16) .eq. 'FMCSC_RIGIDFREQ ') then
      read (wo1(wt11:wt12),*) rigidfreq
      if ((rigidfreq.lt.0.0).OR.(rigidfreq.gt.1.0)) then
        rigidfreq = 0.5
        write(ilog,*) 'Invalid spec. for FMCSC_RIGIDFREQ. Defaulting&
 & to ',rigidfreq,'.' 
      end if
!
!   choice of frequency of rigid-body moves which fully randomize the orientation
!   and position of a molecule in the box 
!
    else if (keyword(1:18) .eq. 'FMCSC_RIGIDRDFREQ ') then
      read (wo1(wt11:wt12),*) rigid_randfreq
      if ((rigid_randfreq.lt.0.0d0).OR.(rigid_randfreq.gt.1.0))&
 & then
        rigid_randfreq = 0.1
        write(ilog,*) 'Invalid spec. for FMCSC_RIGIDRDFREQ. Defaulti&
 &ng to ',rigid_randfreq,'.'
      end if
!
!   whether or not to translate and rotate simultaneously in single-molecule RB moves
!
    else if (keyword(1:18) .eq. 'FMCSC_COUPLERIGID ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        use_coupledrigid = .true.
      else
        use_coupledrigid = .false.
      end if
!
!   amongst single-mol. rigid-body moves, fraction to be rotation moves (only uncoupled)
!
    else if (keyword(1:14) .eq. 'FMCSC_ROTFREQ ') then
      read (wo1(wt11:wt12),*) rotfreq
      if ((rotfreq.lt.0.0).OR.(rotfreq.gt.1.0)) then
        rotfreq = 0.5
        write(ilog,*) 'Invalid spec. for FMCSC_ROTFREQ. Defaulting t&
 &o ',rotfreq,'.' 
      end if
!
!   step size (in A) for all rigid-body translations
!   (interpreted as drawing random vector with length uniform from [0;stepsz])
!
    else if (keyword(1:18) .eq. 'FMCSC_TRANSSTEPSZ ') then
      read (wo1(wt11:wt12),*) trans_stepsz
      if (trans_stepsz.le.0.0d0) then
        trans_stepsz = 5.0
        write(ilog,*) 'Invalid spec. for FMCSC_TRANSSTEPSZ. Defaulti&
 &ng to ',trans_stepsz,'.'
      end if
!
!   step size (in degrees) for all rigid-body rotations
!   (interpreted as drawing rotation angle for each axis from [0;stepsz])
!
    else if (keyword(1:16) .eq. 'FMCSC_ROTSTEPSZ ') then
      read (wo1(wt11:wt12),*) rot_stepsz
      if (rot_stepsz.le.0.0d0) then
        rot_stepsz = 60.0
        write(ilog,*) 'Invalid spec. for FMCSC_ROTSTEPSZ. Defaulting&
 & to ',rot_stepsz,'.'
      end if
!
!   frequency of cluster rigid-body moves amongst total RB moves (-> rigidfreq) 
!
    else if (keyword(1:16) .eq. 'FMCSC_CLURBFREQ ') then
      read (wo1(wt11:wt12),*) clurb_freq
      if ((clurb_freq.lt.0.0).OR.(clurb_freq.gt.1.0)) then
        clurb_freq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_CLURBFREQ. Defaulting&
 & to ',clurb_freq,'.'
      end if
!
!   frequency of random cluster moves amongst total cluster RB moves (-> clurb_freq) 
!
    else if (keyword(1:18) .eq. 'FMCSC_CLURBRDFREQ ') then
      read (wo1(wt11:wt12),*) clurb_rdfreq
      if ((clurb_rdfreq.lt.0.0).OR.(clurb_rdfreq.gt.1.0)) then
        clurb_rdfreq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_CLURBRDFREQ. Defaulti&
 &ng to ',clurb_rdfreq,'.'
      end if
!
!   maximum cluster size for random cluster rigid-body moves  
!
    else if (keyword(1:15) .eq. 'FMCSC_CLURBMAX ') then
     read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.lt.2) then
        clurb_maxsz = 2
        write(ilog,*) 'Invalid spec. for FMCSC_CLURBMAX. Defaulting &
 &to ',clurb_maxsz,'.'
      else
        clurb_maxsz = choice
      end if
!
!   reset frequency for cluster distr. for structural cluster moves, currently DISABLED FUNCTIONALITY
!
    else if (keyword(1:17) .eq. 'FMCSC_CLURBRESET ') then
     read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.lt.2) then
        clurb_reset = 500
        write(ilog,*) 'Invalid spec. for FMCSC_CLURBRESET. Defaultin&
 &g to ',clurb_reset,'.'
      else
        clurb_reset = choice
      end if
!
!   bail-out frequency into random cluster move for structural cluster
!   moves if molecules form single cluster, currently DISABLED FUNCTIONALITY
!
    else if (keyword(1:16) .eq. 'FMCSC_CLURBBAIL ') then
      read (wo1(wt11:wt12),*) clurb_strbail
      if ((clurb_strbail.lt.0.0).OR.(clurb_strbail.gt.1.0)) then
        clurb_strbail = 0.99
        write(ilog,*) 'Invalid spec. for FMCSC_CLURBBAIL. Defaulting&
 & to ',clurb_strbail,'.'
      end if
!
!   cluster definition criterion in A as minimum atom-atom distance, currently DISABLED FUNCTIONALITY
!
    else if (keyword(1:16) .eq. 'FMCSC_CLURBCRIT ') then
      read (wo1(wt11:wt12),*) clurb_discrit
      if (clurb_discrit.lt.0.0) then
        clurb_discrit = 4.0
        write(ilog,*) 'Invalid spec. for FMCSC_CLURBCRIT. Defaulting&
 & to ',clurb_discrit,'.'
      end if
!
!   frequency of particle fluctuation moves 
!      
    else if (keyword(1:23) .eq. 'FMCSC_PARTICLEFLUCFREQ ') then
      read (wo1(wt11:wt12),*) particleflucfreq
      if ((particleflucfreq.lt.0.0).OR.(particleflucfreq.gt.1.0))&
 & then
        particleflucfreq = 0.0
        write(ilog,*) 'Invalid specification for FMCSC_PARTICLEFLUCF&
 &REQ. Defaulting to ', particleflucfreq, '.'
      end if
!
!   input file for constraints requests
!
    else if (keyword(1:14) .eq. 'FMCSC_FRZFILE ') then
      frzfile = wo1(wt11:wt12)
      do_frz = .true.
!
!   whether to write a summary of the constraints in the system
!
    else if (keyword(1:16) .eq. 'FMCSC_FRZREPORT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        frz_report = .true.
      else
        frz_report = .false.
      end if
!
!   whether to exclude fully constrained interactions from force evaluations in dynamics
!
    else if (keyword(1:14) .eq. 'FMCSC_SKIPFRZ ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        skip_frz = .true.
      else
        skip_frz = .false.
      end if
!
!   pH pseudo-move calculations
!
    else if (keyword(1:9) .eq. 'FMCSC_PH ') then
      read (wo1(wt11:wt12),*) phs
      if ((phs.lt.0.0).OR.(phs.gt.14.0))  then
        phs = 7.0
        write(ilog,*) 'Invalid spec. for FMCSC_PH. pH must be between 0 and 14. Defaul&
 &ting to ',phs,'.'
      end if
!
!   pH pseudo-move choice of frequency
!
    else if (keyword(1:13) .eq. 'FMCSC_PHFREQ ') then
      read (wo1(wt11:wt12),*) phfreq
      if ((phfreq.lt.0.0).OR.(phfreq.gt.1.0)) then
        phfreq = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_PHFREQ. Defa&
 &ulting to ',phfreq,'.'
      end if
!
!   ionic strength for pH analysis
!
    else if (keyword(1:15) .eq. 'FMCSC_IONICSTR ') then
      read (wo1(wt11:wt12),*) ionicstr
      if (ionicstr.lt.0.0)  then
        ionicstr = 0.0
        write(ilog,*) 'Invalid spec. for FMCSC_IONICSTR (must be greater or equal to ze&
 &o). Defaulting to ',ionicstr,'.'
      end if
!
!   gradient-based convergence criterion for minimizers
!
    else if (keyword(1:16) .eq. 'FMCSC_MINI_GRMS ') then
      read (wo1(wt11:wt12),*) mini_econv
      if (mini_econv .le. 0.0) then
        mini_econv = 0.1
        write(ilog,*) 'Invalid spec. for FMCSC_MINI_GRMS. Defaulting&
 & to ',mini_econv,'.'
      end if
!
!   uphill tolerance for BFGS minimizer
!
    else if (keyword(1:17) .eq. 'FMCSC_MINI_UPTOL ') then
      read (wo1(wt11:wt12),*) mini_uphill
      if (mini_uphill.lt.0.0) then
        mini_uphill = 0.1
        write(ilog,*) 'Invalid spec. for FMCSC_MINI_UPTOL. Defaulting&
 & to ',mini_uphill,'.'
      end if
! 
!   memory length of BFGS minimizer
!
    else if (keyword(1:18) .eq. 'FMCSC_MINI_MEMORY ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.lt.5).OR.(choice.gt.2000)) then
        mini_mem = 10
        write(ilog,*) 'Invalid spec. for FMCSC_MINI_MEMORY. Defaulting&
 & to ',mini_mem,'.'
      else
        mini_mem = choice
      end if
!          
!   step-sizes for minimizers
!
    else if (keyword(1:20) .eq. 'FMCSC_MINI_STEPSIZE ') then
      read (wo1(wt11:wt12),*) mini_stepsize
      if (mini_stepsize .le. 0.0) then
        mini_stepsize = 0.5
        write(ilog,*) 'Invalid spec. for FMCSC_MINI_STEPSIZE. Defaul&
 &ting to ',mini_econv,'.'
      end if
!         
    else if (keyword(1:24) .eq. 'FMCSC_MINI_XYZ_STEPSIZE ') then
      read (wo1(wt11:wt12),*) mini_xyzstep
      if (mini_xyzstep .le. 0.0) then
        mini_xyzstep = 0.1
        write(ilog,*) 'Invalid spec. for FMCSC_MINI_XYZ_STEPSIZE. Defaul&
 &ting to ',mini_xyzstep,'.'
      end if
!
    else if (keyword(1:24) .eq. 'FMCSC_MINI_ROT_STEPSIZE ') then
      read (wo1(wt11:wt12),*) mini_rotstep
      if (mini_rotstep .le. 0.0) then
        mini_rotstep = 0.5
        write(ilog,*) 'Invalid spec. for FMCSC_MINI_ROT_STEPSIZE. Defaul&
 &ting to ',mini_rotstep,'.'
      end if
!
    else if (keyword(1:24) .eq. 'FMCSC_MINI_INT_STEPSIZE ') then
      read (wo1(wt11:wt12),*) mini_intstep
      if (mini_intstep .le. 0.0) then
        mini_intstep = 0.5
        write(ilog,*) 'Invalid spec. for FMCSC_MINI_INT_STEPSIZE. Defaul&
 &ting to ',mini_intstep,'.'
      end if
! 
!   type of minimizer
!
    else if (keyword(1:16) .eq. 'FMCSC_MINI_MODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.le.0).OR.(choice.gt.4)) then
        mini_mode = 1
        write(ilog,*) 'Invalid spec. for FMCSC_MINI_MODE. Defaulting&
 & to ',mini_mode,'.'
      else
        mini_mode = choice
      end if
!
!   number of SD steps prior to stochastic minimization / annealing
!
    else if (keyword(1:22) .eq. 'FMCSC_MINI_SC_SDSTEPS ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        mini_sc_sdsteps = 0
        write(ilog,*) 'Invalid spec. for FMCSC_MINI_SC_SDSTEPS. Defaulting&
 & to ',mini_sc_sdsteps,'.'
      else
        mini_sc_sdsteps = choice
      end if
!
!   fraction of stochastic minimizer / annealing protocol spent in hot phase
!
    else if (keyword(1:19) .eq. 'FMCSC_MINI_SC_HEAT ') then
      read (wo1(wt11:wt12),*) mini_sc_heat
      if (mini_sc_heat .le. 0.001) then
        mini_sc_heat = 0.001
        write(ilog,*) 'Invalid spec. for FMCSC_MINI_SC_HEAT. Defaul&
 &ting to ',mini_sc_heat,'.'
      end if
      if (mini_sc_heat .ge. 0.999) then
        mini_sc_heat = 0.999
        write(ilog,*) 'Invalid spec. for FMCSC_MINI_SC_HEAT. Defaul&
 &ting to ',mini_sc_heat,'.'
      end if
! 
!   target temperature for quenching of stochastic minimizer / annealing protocol
!
    else if (keyword(1:20) .eq. 'FMCSC_MINI_SC_TBATH ') then
      read (wo1(wt11:wt12),*) mini_sc_tbath
      if (mini_sc_tbath .le. 0.0) then
        mini_sc_tbath = 0.1
        write(ilog,*) 'Invalid spec. for FMCSC_MINI_SC_HEAT. Defaul&
 &ting to ',mini_sc_tbath,'.'
      end if
!
    end if
!
!
  else if (which.eq.2) then
!
!
!   Sjunesson, Ulmschneider, etc. concerted rotation spec.s 
!
    if ((nucfreq.gt.0.0).OR.(crfreq.gt.0.0)) then
!
!     the variance or width of the distribution we're pulling from in SJCR
!
      if (keyword(1:14) .eq. 'FMCSC_CRWIDTH ') then
        read (wo1(wt11:wt12),*) dum
        if ((dum .le. 0.0d0).OR.(dum .gt. (PI/2))) then
          cr_a = 1.0/(PI/80.0)
          write(ilog,*) 'Invalid spec. for FMCSC_CRWIDTH. Defaulting to ',cr_a,'.'
        else
          cr_a = 1.0/dum
        end if
!
!     the scaling factor for the biasing term to be put in (co-regulates width!!!) in SJCR
!
      else if (keyword(1:13) .eq. 'FMCSC_CRBIAS ') then
        read (wo1(wt11:wt12),*) dum
        if (dum .le. 0.0d0) then
          cr_b = 1.0
          write(ilog,*) 'Invalid spec. for FMCSC_CRBIAS. Defaulting to ',cr_b,'.'
        else
          cr_b = dum
        end if
!
!     SJCR mode: 1) CNC of next (first unsampled) residue as reference, derivatives via eff. lever arms
!                2) CCO of final (last sampled) residue as reference, derivatives via nested rotation matrices
!
      else if (keyword(1:13) .eq. 'FMCSC_CRMODE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.le.0).OR.(choice.gt.2)) then
          cr_mode = 1
          write(ilog,*) 'Invalid spec. for FMCSC_CRMODE. Deefaulting to ',cr_mode,'.'
        else
          cr_mode = choice
        end if
!
!     the number of backbone degrees of freedom used to approximate the solution of the loop closure problem for SJCR
!
      else if (keyword(1:12) .eq. 'FMCSC_CRDOF ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice .lt. 7) then
          nr_crdof = 8
          write(ilog,*) 'Invalid spec. for FMCSC_CRDOF. Defaulting to ',nr_crdof,'.'
        else
          nr_crdof = choice
        end if
        if (mod(nr_crdof,2).eq.0) then
          nr_crres = nr_crdof/2
        else
          nr_crres = nr_crdof/2 + 1
        end if
!
!     the variance or width of the distribution we're pulling from for exact CR
!
      else if (keyword(1:15) .eq. 'FMCSC_UJCRBIAS ') then
        read (wo1(wt11:wt12),*) dum
        if (dum.lt.0.0) then
          UJ_params(2) = 8.0
          write(ilog,*) 'Invalid spec. for FMCSC_UJCRBIAS. Defaulting to ',UJ_params(2),'.'
        else
          UJ_params(2) = dum
        end if
!
!     the net adjusting factor for the width of the resultant distribution in the pre-rotation segment for exact CR
!
      else if (keyword(1:16) .eq. 'FMCSC_UJCRWIDTH ') then
        read (wo1(wt11:wt12),*) dum
        if ((dum.le.0.0).OR.(dum.gt.90.0)) then
          UJ_params(1) = 2.0
          write(ilog,*) 'Invalid spec. for FMCSC_UJCRWIDTH. Defaulting to ',UJ_params(1),'.'
        else
          UJ_params(1) = dum
        end if
        UJ_params(1) = UJ_params(1)/RADIAN
        UJ_params(1) = 1.0/(UJ_params(1)**2.0)
!
!     the relative adjustment factor for scaling down the magnitude of the
!     (biased) fluctuations on the angles vs. the dihedrals in the pre-rotation segment of UJ-CR moves
!
      else if (keyword(1:16) .eq. 'FMCSC_UJCRSCANG ') then
        read (wo1(wt11:wt12),*) dum
        if (dum.le.0.0) then
          UJ_params(3) = 20.0
          write(ilog,*) 'Invalid spec. for FMCSC_UJCRSCANG. Defaulting to ',1.0/UJ_params(3),'.'
        else
          UJ_params(3) = 1.0/dum
        end if
!
!     the same for pre-rotation omega angles in exact torsional CR methods
!
      else if (keyword(1:17) .eq. 'FMCSC_TORCRSCOME ') then
        read (wo1(wt11:wt12),*) dum
        if (dum.le.0.0) then
          UJ_params(6) = 1.0
          write(ilog,*) 'Invalid spec. for FMCSC_TORCRSCOME. Defaulting to ',1.0/UJ_params(6),'.'
        else
          UJ_params(6) = 1.0/dum
        end if
!
!     the maximum length of the pre-rotation segment for bond angle UJCR
!
      else if (keyword(1:14) .eq. 'FMCSC_UJCRMAX ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.18)) then
          ujmaxsz = 8
          write(ilog,*) 'Invalid spec. for FMCSC_UJCRMAX. Defaulting to ',ujmaxsz-2,'.'
        else
          ujmaxsz = choice + 2
        end if
!
!     the minimum length of the pre-rotation segment for bond angle UJCR
!
      else if (keyword(1:14) .eq. 'FMCSC_UJCRMIN ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.18)) then
          ujminsz = 3
          write(ilog,*) 'Invalid spec. for FMCSC_UJCRMIN. Defaulting to ',ujminsz-2,'.'
        else
          ujminsz = choice + 2
        end if
!
!     the search stepsize for the 1-D root search for exact CR
!
      else if (keyword(1:17) .eq. 'FMCSC_UJCRSTEPSZ ') then
        read (wo1(wt11:wt12),*) dum
        if ((dum.lt.0.001).OR.(dum.gt.50.0)) then
          UJ_params(4) = 1.0
          write(ilog,*) 'Invalid spec. for FMCSC_UJCRSTEPSZ. Defaulting to ',UJ_params(4),'.'
        else
          UJ_params(4) = dum
        end if
!
!     interval size restriction for root search for bond angle UJCR
!
      else if (keyword(1:19) .eq. 'FMCSC_UJCRINTERVAL ') then
        read (wo1(wt11:wt12),*) dum
        if ((dum.lt.3.0).OR.(dum.gt.90.0)) then
          UJ_params(5) = 20.0
          write(ilog,*) 'Invalid spec. for FMCSC_UJCRINTERVAL. Defaulting to ',UJ_params(5),'.'
        else
          UJ_params(5) = dum
        end if
!
!     the number of iterations before bailing out for exact torsional CR
!
      else if (keyword(1:17) .eq. 'FMCSC_UJMAXTRIES ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.lt.1) then
          maxcrtries = 1
          write(ilog,*) 'Invalid spec. for FMCSC_UJMAXTRIES. Defaulting to ',maxcrtries,'.'
        else
          maxcrtries = choice
        end if
!
!     pre-rotation segment length settings for exact torsional CR methods
! 
      else if (keyword(1:18) .eq. 'FMCSC_TORCRMAX_DO ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.MAXUJDOF)) then
          torcrmaxsz = 12
          write(ilog,*) 'Invalid spec. for FMCSC_TORCRMAX. Defaulting to ',torcrmaxsz,'.'
        else
          torcrmaxsz = choice
        end if
      else if (keyword(1:18) .eq. 'FMCSC_TORCRMAX_DJ ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.MAXUJDOF)) then
          torcrmaxsz2 = 12
          write(ilog,*) 'Invalid spec. for FMCSC_TORCRMAX_DJ. Defaulting to ',torcrmaxsz2,'.'
        else
          torcrmaxsz2 = choice
        end if
      else if (keyword(1:15) .eq. 'FMCSC_NUCCRMAX ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.18)) then
          nuccrmaxsz = 10
          write(ilog,*) 'Invalid spec. for FMCSC_NUCCRMAX. Defaulting to ',nuccrmaxsz,'.'
        else
          nuccrmaxsz = choice
        end if
      else if (keyword(1:18) .eq. 'FMCSC_TORCRMIN_DO ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.MAXUJDOF)) then
          torcrminsz = 1
          write(ilog,*) 'Invalid spec. for FMCSC_TORCRMIN. Defaulting to ',torcrminsz,'.'
        else
          torcrminsz = choice
        end if
      else if (keyword(1:18) .eq. 'FMCSC_TORCRMIN_DJ ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.MAXUJDOF)) then
          torcrminsz2 = 1
          write(ilog,*) 'Invalid spec. for FMCSC_TORCRMIN_DJ. Defaulting to ',torcrminsz2,'.'
        else
          torcrminsz2 = choice
        end if
      else if (keyword(1:15) .eq. 'FMCSC_NUCCRMIN ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.MAXUJDOF)) then
          nuccrminsz = 1
          write(ilog,*) 'Invalid spec. for FMCSC_NUCCRMIN. Defaulting to ',nuccrminsz,'.'
        else
          nuccrminsz = choice
        end if
!
      end if
!
    end if !(end if for (crfreq.gt.0))
!
!   Use of steric exclusion grids (dipeptide derived)
!
    if (use_stericgrids.EQV..true.) then
!
!     grid spacing (half window size)
! 
      if (keyword(1:16) .eq. 'FMCSC_GRDWINDOW ') then
        read (wo1(wt11:wt12),*) stgr%halfwindow
        if ((stgr%halfwindow.le.0.0).OR.(stgr%halfwindow.gt.180.0)) then
          stgr%halfwindow = 5.0
          write(ilog,*) 'Invalid spec. for FMCSC_GRDWINDOW. Defaulting to ',stgr%halfwindow,'.'
        end if
!
!     Directory with steric grids
!
      else if (keyword(1:14) .eq. 'FMCSC_GRIDDIR ') then
        if ((wt12.lt.MAXKEYLEN).AND.(wo1(wt12:wt12).ne.SLASHCHAR)) then
          wo1(wt12+1:wt12+1)=SLASHCHAR
          wt12 = wt12 + 1
        end if
        str2 = wo1(wt11:wt12)//'ala_grid.dat'
        inquire(file=str2,exist=exists)
        if (exists.EQV..false.) then
          write(ilog,*) 'Cannot open file ala_grid.dat in directory &
 &',wo1(wt11:wt12),'. This is fatal. Check FMCSC_GRIDDIR!'
          call fexit()
        else
          griddir = wo1(wt11:wt12)
        end if
!
      end if
!
    end if !(end if for use_stericgrids)
!
!   Use of Wang-Landau sampling
!
    if (do_wanglandau.EQV..true.) then
!
!     whether to dynamically grow histograms (1 = no, 2 = only low end, 3 = both ends)
!
      if (keyword(1:16) .eq. 'FMCSC_WL_EXTEND ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.le.0).OR.(choice.gt.3)) then
          wld%exrule = 1
          write(ilog,*) 'Invalid spec. for FMCSC_WL_EXTEND. Defaulting to ',wld%exrule,'.'
        else
          wld%exrule = choice
        end if
!
!     how many buffer steps for f-updates and whether to freeze the histogram
!
      else if (keyword(1:16) .eq. 'FMCSC_WL_FREEZE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.0) then
          wld%buffer = max(1,nsim/10)
          write(ilog,*) 'Invalid spec. for FMCSC_WL_FREEZE. Defaulting to ',wld%buffer,'.'
        else
          if (choice.le.0) then
            wld%buffer = -choice
            wld%freeze = .false.
          else
            wld%buffer = choice
          end if
        end if
!            
!     whether to print out intermediate quantities
!
      else if (keyword(1:15) .eq. 'FMCSC_WL_DEBUG ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.1) debug_wanglandau = .true.
!            
!     Number of steps to check if reference histogram has enough visits
!
      else if (keyword(1:19) .eq. 'FMCSC_WL_FLATCHECK ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.le.0).OR.(choice.gt.(nsim/10))) then
          wld%gh_flatcheck_freq = 10000
          write(ilog,*) 'Invalid spec. for FMCSC_WL_FLATCHECK. Defaulting to ',wld%gh_flatcheck_freq,'.'
        else
          wld%gh_flatcheck_freq = choice
        end if
!            
!     Update histograms only after this many steps
!
      else if (keyword(1:16) .eq. 'FMCSC_WL_HUFREQ ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.le.0).OR.(choice.gt.(nsim/10))) then
          wld%hufreq = 10
          write(ilog,*) 'Invalid spec. for FMCSC_WL_HUFREQ. Defaulting to ',wld%hufreq,'.'
        else
          wld%hufreq = choice
        end if
!            
!     Initial log(f) value
!
      else if (keyword(1:12) .eq. 'FMCSC_WL_F0 ') then
        read (wo1(wt11:wt12),*) wld%fval
        if (wld%fval.le.0.0) then
          wld%fval = 1.0
          write(ilog,*) 'Invalid spec. for FMCSC_WL_F0. Defaulting to ',wld%fval,'.'
        end if
!
!     How often do we need to visit each energy bin before adjusting the F value?
!
      else if (keyword(1:16) .eq. 'FMCSC_WL_HVMODE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.le.0).OR.(choice.gt.2)) then
          wld%hvmode = 1
          write(ilog,*) 'Invalid spec. for FMCSC_WL_HVMODE. Defaulting to ',wld%hvmode,'.'
        else
          wld%hvmode = choice
        end if
!
!     input file to intitialize g(E) 
!
      else if(keyword(1:19) .eq. 'FMCSC_WL_GINITFILE ') then
        wld%ginitfile = wo1(wt11:wt12)
        wld%use_ginitfile = .true.
!            
!     which molecule to use for reaction coordinate-based WL (modes 2, 3 only)
!
      else if (keyword(1:13) .eq. 'FMCSC_WL_MOL ') then
        if (wld%wl_mode.ge.2) then
          read (wo1(wt11:wt12),*) dum
          call dum2int(dum,choice)
          if (choice.le.0) then
            if (wld%wl_mode.eq.3) wld%dimensionality = 1
            write(ilog,*) 'Invalid spec. for FMCSC_WL_MOL (for first molecule). This means that either the &
 &default is used (mol. ',wld%wl_mol(1),' for FMCSC_WL_MODE 2) or that the 2D request is ignored (for FMCSC_WL_MODE 3).'
          else
            wld%wl_mol(1) = choice
            if (wld%wl_mode.eq.3) wld%dimensionality = 2
            if (wld%wl_mode.eq.2) wld%dimensionality = 1
          end if
        end if
        if ((wld%wl_mode.eq.2).AND.(wt22.ge.wt21)) then
          if ((wo2(wt21:wt21).ne.' ').AND.(wo2(wt21:wt21).ne.'#')) then
            read (wo2(wt21:wt22),*) dum
            call dum2int(dum,choice)
            if (choice.le.0) then
              wld%dimensionality = 1
              write(ilog,*) 'Invalid spec. for FMCSC_WL_MOL for second molecule. This means that the request for 2D is &
 &ignored (for FMCSC_WL_MODE 2).'
            else
              wld%dimensionality = 2
              wld%wl_mol(2) = choice
            end if
          end if
        else if (wld%wl_mode.eq.2) then
          wld%dimensionality = 1 
        end if
!            
!     what reaction coordinate(s) to use for reaction coordinate-based or mixed WL
!
      else if (keyword(1:12) .eq. 'FMCSC_WL_RC ') then
        if (wld%wl_mode.ge.2) then
          read (wo1(wt11:wt12),*) dum
          call dum2int(dum,choice)
          if ((choice.le.0).OR.(choice.gt.3)) then
            write(ilog,*) 'Invalid spec. for FMCSC_WL_RC (for first molecule). Using default (radius of gyration). This may &
 &produce an error later.'
          else
            wld%wl_rc(1) = choice
          end if
        end if
        if ((wld%wl_mode.eq.2).AND.(wt22.ge.wt21)) then
          if ((wo2(wt21:wt21).ne.' ').AND.(wo2(wt21:wt21).ne.'#')) then
            read (wo2(wt21:wt22),*) dum
            call dum2int(dum,choice)
            if ((choice.le.0).OR.(choice.gt.3)) then
              wld%wl_rc(2) = 1
              write(ilog,*) 'Invalid spec. for FMCSC_WL_RC for second dimension. Using default (radius of gyration). This may &
 &produce an error later.'
            else
              wld%wl_rc(2) = choice
            end if
          end if
        end if
!
      end if
!
    end if ! if do_wanglandau
!
!
  else if (which.eq.3) then
!
    if (do_wanglandau.EQV..true.) then
!
!     WL maximum allowable value (energy or order parameter)
!
      if (keyword(1:13) .eq. 'FMCSC_WL_MAX ') then
        if (wld%dimensionality.eq.1) then
          read (wo1(wt11:wt12),*) wld%g_max
        else if (wld%dimensionality.eq.2) then
          read (wo1(wt11:wt12),*) wld%g_max2d(1)
          wld%g_max = wld%g_max2d(1)
          if (wt22.ge.wt21) then
            if ((wo2(wt21:wt21).ne.' ').AND.(wo2(wt21:wt21).ne.'#')) then
              read (wo2(wt21:wt22),*) wld%g_max2d(2)
            end if
          end if
        end if
!
!     WL bin size (energy or order parameter)
!   
      else if (keyword(1:15).eq.'FMCSC_WL_BINSZ ') then
        if (wld%dimensionality.eq.1) then
          read (wo1(wt11:wt12),*) wld%g_binsz
          if (wld%g_binsz.le.0.0) then
            wld%g_binsz = 1.0
            write(ilog,*) 'Invalid spec. for FMCSC_WL_BINSZ. Defaulting to ',wld%g_binsz,'.'
          end if
        else if (wld%dimensionality.eq.2) then
          read (wo1(wt11:wt12),*) wld%g_binsz2d(1)
          if (wt22.ge.wt21) then
            if ((wo2(wt21:wt21).ne.' ').AND.(wo2(wt21:wt21).ne.'#')) then
              read (wo2(wt21:wt22),*) wld%g_binsz2d(2)
            end if
          end if
          if (wld%g_binsz2d(1).le.0.0) then
            wld%g_binsz2d(1) = 0.1
            write(ilog,*) 'Invalid spec. for FMCSC_WL_BINSZ (first dimension). Defaulting to ',wld%g_binsz2d(1),'.'
          end if
          if (wld%g_binsz2d(2).le.0.0) then
            wld%g_binsz2d(2) = 0.1
            write(ilog,*) 'Invalid spec. for FMCSC_WL_BINSZ (second dimension). Defaulting to ',wld%g_binsz2d(2),'.'
          end if
        end if
!
      end if
!
    end if ! if do_wanglandau
!
  else
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! covers: params, inter, fos
!
subroutine key_params(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 & rest,which,keyword)
!
  use iounit
  use params
  use fos
  use keys
  use inter
!
  implicit none
!
  character(MAXKEYLEN) wo1,wo2,wo3,rest
  character(MAXKWLEN) keyword
  integer which,choice,wt11,wt12,wt21,wt22,wt31,wt32
  RTYPE dum
  logical exists
!
  if (which.eq.1) then
!
!
!   parameter file
!
    if (keyword(1:11) .eq. 'PARAMETERS ') then
      inquire(file=wo1(wt11:wt12),exist=exists)
      if (exists.EQV..false.) then
        write(ilog,*) 'Fatal. Cannot open spec.d parameter file (',wo1(wt11:wt12),') for PARAMETERS.'
        call fexit()
      else
        paramfile = wo1(wt11:wt12)
      end if
!
!   patch file for biotype (re)assignments
!
    else if (keyword(1:23) .eq. 'FMCSC_BIOTYPEPATCHFILE ') then
      biotpatchfile = wo1(wt11:wt12)
!
!   patch file for atomic masses
!
    else if (keyword(1:17) .eq. 'FMCSC_MPATCHFILE ') then
      masspatchfile = wo1(wt11:wt12)
!
!   patch file for atomic radii
!
    else if (keyword(1:17) .eq. 'FMCSC_RPATCHFILE ') then
      radpatchfile = wo1(wt11:wt12)
!
!   patch file for LJ types ("atom" entries in parameter file)
!
    else if (keyword(1:18) .eq. 'FMCSC_LJPATCHFILE ') then
      ljpatchfile = wo1(wt11:wt12)
!
!   LJ Sigma combination rule
!
    else if (keyword(1:14) .eq. 'FMCSC_SIGRULE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.gt.3).OR.(choice.lt.1)) then
        sigrule = 1
        write(ilog,*) 'Invalid spec. for FMCSC_SIGRULE. Defaulting t&
 &o arithmetic rule (',sigrule,').'
      else
        sigrule = choice
      end if
!
!   LJ Epsilon combination rule
!
    else if (keyword(1:14) .eq. 'FMCSC_EPSRULE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.gt.3).OR.(choice.lt.1)) then
        epsrule = 2
        write(ilog,*) 'Invalid spec. for FMCSC_EPSRULE. Defaulting t&
 &o geometric rule (',epsrule,').'
      else
        epsrule = choice
      end if
!
!   HS scaling factor
!
    else if (keyword(1:14) .eq. 'FMCSC_HSSCALE ') then
      read (wo1(wt11:wt12),*) reduce
      if (reduce .le. 0.0d0)  then
        reduce = 1.0d0
      end if
!
!   whether to write out an overview of the vdW parameters
!
    else if (keyword(1:16) .eq. 'FMCSC_VDWREPORT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        vdW_report = .true.
      else
        vdW_report = .false.
      end if
!
!   Mode for NB-SR interactions: 1 = all atoms separated by at least 3 bonds, and at least
!                                    one relevant rotatable bond interact. the omega bond in amides
!                                    is not considered a relevant, rotatable bond!
!                                2 = all atoms separated by at least 3 bonds interact
!                                3 = special GROMOS exclusion rules
!
    else if (keyword(1:17) .eq. 'FMCSC_INTERMODEL ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.gt.3).OR.(choice.lt.1)) then
        nbsr_model = 1
        write(ilog,*) 'Invalid spec. for FMCSC_INTERMODEL. Defaultin&
 &g to ',nbsr_model,'.'
      else
        nbsr_model = choice
      end if
!
!   Mode for 1-4 interactions: 1 = only true 14 are treated 14)
!                              2 = all pairs connected through a single rotatable
!                                  bond are treated 14. besides aromatic rings, gdn, 
!                                  this affects amides (omega always assumed pseudo-rigid)
!
    else if (keyword(1:14) .eq. 'FMCSC_MODE_14 ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,mode_14)
      if ((mode_14.lt.1).OR.(mode_14.gt.2)) then
        mode_14 = 1
        write(ilog,*) 'Invalid spec. for FMCSC_MODE_14. Defaulting t&
 &o ',mode_14,'.'
      end if
!
!   Scaling factor for 1-4 specifically for attLJ, WCA, and IPP terms
!
    else if (keyword(1:18) .eq. 'FMCSC_FUDGE_ST_14 ') then
      read (wo1(wt11:wt12),*) fudge_st_14
      if (fudge_st_14.lt.0.0d0)  then
        fudge_st_14 = 1.0
        write(ilog,*) 'Invalid spec. for FMCSC_FUDGE_ST_14. Defaulti&
 &ng to ',fudge_st_14,'.'
      end if
!
!   Scaling factor for 1-4 specifically for POLAR (but not for TABUL!) terms
!
    else if (keyword(1:18) .eq. 'FMCSC_FUDGE_EL_14 ') then
      read (wo1(wt11:wt12),*) fudge_el_14
      if (fudge_el_14.lt.0.0d0)  then
        fudge_el_14 = 1.0
        write(ilog,*) 'Invalid spec. for FMCSC_FUDGE_EL_14. Defaulti&
 &ng to ',fudge_el_14,'.'
      end if
!
!   whether to write out an overview of the non-bonded interactions in the system
!
    else if (keyword(1:18) .eq. 'FMCSC_INTERREPORT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        ia_report = .true.
      else
        ia_report = .false.
      end if
!
!   whether to write out an overview of the electrostatic interactions in the system
!
    else if (keyword(1:17) .eq. 'FMCSC_ELECREPORT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        elec_report = .true.
      else
        elec_report = .false.
      end if
!
!   a report flag whether to write out a summary of found and missing parameters for bonded terms
!
    else if (keyword(1:17) .eq. 'FMCSC_BONDREPORT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        bonded_report = .true.
      else
        bonded_report = .false.
      end if 
!
!    whether to write out an overview of the FOS reference values
!
    else if (keyword(1:16) .eq. 'FMCSC_FOSREPORT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        fos_report = .true.
      else
        fos_report = .false.
      end if
!
!   whether to write out an overview of the automatically determined dipole groups (VMD-file)
!
    else if (keyword(1:16) .eq. 'FMCSC_DIPREPORT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        dip_report = .true.
      else
        dip_report = .false.
      end if
!
!   which short-range electrostatic base model to use: 1) MMFF 2) sane charge group-based approach
!
    else if (keyword(1:16) .eq. 'FMCSC_ELECMODEL ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.gt.2).OR.(choice.lt.1)) then
        elec_model = 1
        write(ilog,*) 'Invalid spec. for FMCSC_ELECMODEL. Defaulting&
 & to ',elec_model,'.'
      else
        elec_model = choice
      end if
!
!   whether to turn a few fatal sanity checks into warnings instead
!
    else if (keyword(1:13) .eq. 'FMCSC_UNSAFE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        be_unsafe = .true.
      else
        be_unsafe = .false.
      end if
!
    end if
!
!
  else if (which.eq.2) then
!
  else if (which.eq.3) then
!
  else
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! covers: polyavg, diffrac, ems
!
subroutine key_polyavg(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 & rest,which,keyword)
!
  use iounit
  use polyavg
  use diffrac
  use ems
  use keys
!
  implicit none
!
  character(MAXKEYLEN) wo1,wo2,wo3,rest
  character(MAXKWLEN) keyword
  integer which,choice,wt11,wt12,wt21,wt22,wt31,wt32
  RTYPE dum
! we'll use an old "size"-variable here to prevent users from
! asking for arbitrary (and ridicously unnecessary) memory
  integer MAXQVEC
  parameter (MAXQVEC=1000)
!
  if (which.eq.1) then
!
!
!   output frequency for system-wide polymer measures (rarely useful)
!
    if (keyword(1:13) .eq. 'FMCSC_POLOUT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        write(ilog,*) 'Invalid spec. for FMCSC_POLOUT. Default.'
        polout = 500
      else
        polout = choice
      end if
!
!   calculation frequency for averaged polymer prop.s
!
    else if (keyword(1:14) .eq. 'FMCSC_POLCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        polcalc = 50
        write(ilog,*) 'Invalid spec. for FMCSC_POLCALC. Defaulting t&
 &o ',polcalc,'.'
      else
        polcalc = choice
      end if
!
!   maximum number of Rg-bins: this affects RGHIST, RETEHIST, and DENSPROF
!
    else if (keyword(1:16) .eq. 'FMCSC_POLRGBINS ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.1) then
        maxrgbins = 1000
        write(ilog,*) 'Invalid spec. for FMCSC_POLRGBINS. Defaulting&
 & to ',maxrgbins,'.'
      else
        maxrgbins = choice
      end if
!
!   subfrequency of calculation of scattering intensities (of FMCSC_RHCALC)
!
    else if (keyword(1:18) .eq. 'FMCSC_SCATTERCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        sctcalc = 10
        write(ilog,*) 'Invalid spec. for FMCSC_SCATTERCALC. Defaulti&
 &ng to ',sctcalc,'.'
      else
        sctcalc = choice
      end if
!
!   number of wave vectors to use
!
    else if (keyword(1:18) .eq. 'FMCSC_SCATTERVECS ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.le.0).OR.(choice.gt.MAXQVEC)) then
        nqv = 20
        write(ilog,*) 'Invalid spec. for FMCSC_SCATTERVECS. Defaulti&
 &ng to ',sctcalc,'.'
      else
        nqv = choice
      end if
!
!  maximum wave vectors length to cover
!
    else if (keyword(1:17) .eq. 'FMCSC_SCATTERRES ') then
      read (wo1(wt11:wt12),*) qv_res
      if (qv_res.le.0) then
        qv_res = 0.025
        write(ilog,*) 'Invalid spec. for FMCSC_SCATTERRES. Defaultin&
 &g to ',qv_res,'.'
      end if
!
!   frequency of calculation of internal void space (holes)
!
    else if (keyword(1:16) .eq. 'FMCSC_HOLESCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        holescalc = HUGE(holescalc)
        write(ilog,*) 'Invalid spec. for FMCSC_HOLESCALC. Defaulting&
 & to ',holescalc,'.'
      else
        holescalc = choice
      end if
!
!   whether to output instantaneous values of cosines (for calculation of persistence length) instead of just a final average
!
    else if (keyword(1:16) .eq. 'FMCSC_INST_PERS ') then
      inst_pers = .true.
!     
!   calculation frequency for averaged inverse of hydrodynamic radius 
!
    else if (keyword(1:13) .eq. 'FMCSC_RHCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        write(ilog,*) 'Invalid spec. for FMCSC_RHCALC. Default.'
        rhcalc = 500
      else
        rhcalc = choice
      end if
!
!   bin-size for Rg axis of Rg vs. delta* histograma and for Ree-histograms
!
    else if (keyword(1:16) .eq. 'FMCSC_RGBINSIZE ') then
      read (wo1(wt11:wt12),*) rg_binsz
      if (rg_binsz.le.0.0) then
        rg_binsz = 0.2
        write(ilog,*) 'Invalid spec. for FMCSC_RGBINSIZE. Defaulting&
 & to ',rg_binsz,'.'
      end if
!
!   frequency of calculation of spatial density maps
!
    else if (keyword(1:13) .eq. 'FMCSC_EMCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        emcalc = huge(emcalc)
        write(ilog,*) 'Invalid spec. for FMCSC_EMCALC. Defaulting&
 & to ',emcalc,'.'
      else
        emcalc = choice
      end if
!
!   grid resolution vector for EM analysis
!
    else if (keyword(1:15) .eq. 'FMCSC_EMDELTAS ') then
      read (wo1(wt11:wt12),*) emgrid%deltas(1)
      read (wo2(wt21:wt22),*) emgrid%deltas(2)
      read (wo3(wt31:wt32),*) emgrid%deltas(3)
!
!   assumed density of background
!
    else if (keyword(1:18) .eq. 'FMCSC_EMBGDENSITY ') then
      read (wo1(wt11:wt12),*) embgdensity
      if (embgdensity.lt.0.0) then
        embgdensity = 1.0
        write(ilog,*) 'Invalid spec. for FMCSC_EMBGDENSITY. Defaulting&
 & to ',embgdensity,'.'
      end if
!
!   order for cardinal B-splines
!
    else if (keyword(1:16) .eq. 'FMCSC_EMBSPLINE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.lt.1) then
        emsplor = 3
        write(ilog,*) 'Invalid spec. for FMCSC_EMBSPLINE. Defaulting&
 & to ',emsplor,'.'
      else
        emsplor = choice
      end if
!
!   property selector for density restraint potential
!
    else if (keyword(1:17) .eq. 'FMCSC_EMPROPERTY ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.lt.1).OR.(choice.gt.3)) then
        empotprop = 1
        write(ilog,*) 'Invalid spec. for FMCSC_EMPROPERTY. Defaulting&
 & to ',empotprop,'.'
      else
        empotprop = choice
      end if
!
!   input file for user-defined property
!
    else if (keyword(1:21) .eq. 'FMCSC_EMPROPERTYFILE ') then
      read (wo1(wt11:wt12),*) empropfile
!
!   frequency of calculation of diffraction patterns
!
    else if (keyword(1:16) .eq. 'FMCSC_DIFFRCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        diffrcalc = HUGE(diffrcalc)
        write(ilog,*) 'Invalid spec. for FMCSC_DIFFRCALC. Defaulting&
 & to ',diffrcalc,'.'
      else
        diffrcalc = choice
      end if
!
!   resolution of diffraction map in axial coordinate
!
    else if (keyword(1:16) .eq. 'FMCSC_DIFFRZRES ') then
      read (wo1(wt11:wt12),*) zz_res
      if (zz_res.le.0) then
        zz_res = 0.1
        write(ilog,*) 'Invalid spec. for FMCSC_DIFFRZRES. Defaulting&
 & to ',zz_res,'.'
      end if
!
!   resolution of diffraction map in radial coordinate
!
    else if (keyword(1:16) .eq. 'FMCSC_DIFFRRRES ') then
      read (wo1(wt11:wt12),*) rr_res
      if (rr_res.le.0) then
        rr_res = 0.1
        write(ilog,*) 'Invalid spec. for FMCSC_DIFFRRRES. Defaulting&
 & to ',rr_res,'.'
      end if
!
!   number of diffraction bins in axial coordinate
!
    else if (keyword(1:16) .eq. 'FMCSC_DIFFRZMAX ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice .le. 1) then
        zz_max = 100
        write(ilog,*) 'Invalid spec. for FMCSC_DIFFRZMAX. Defaulting&
 & to ',zz_max,'.'
      else
        zz_max = choice
      end if
!
!   number of diffraction bins in radial coordinate
!
    else if (keyword(1:16) .eq. 'FMCSC_DIFFRRMAX ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice .le. 1) then
        rr_max = 100
        write(ilog,*) 'Invalid spec. for FMCSC_DIFFRRMAX. Defaulting&
 & to ',rr_max,'.'
      else
        rr_max = choice
      end if
!
!   highest order of bessel fxn to use
!
    else if (keyword(1:16) .eq. 'FMCSC_DIFFRJMAX ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice .lt. 0) then
        bes_max = 10
        write(ilog,*) 'Invalid spec. for FMCSC_DIFFRJMAX. Defaulting&
 & to ',bes_max,'.'
      else
        bes_max = choice
      end if
!
!   if this keyword is present, the user requests a custom, constant axis
!   for diffaction calculations
!
    else if (keyword(1:16) .eq. 'FMCSC_DIFFRAXIS ') then
      diffr_uax = .true.
      read (wo1(wt11:wt12),*) diffr_axis(1)
      read (wo2(wt21:wt22),*) diffr_axis(2)
      read (wo3(wt31:wt32),*) diffr_axis(3)
!
!   this keyword is redundant if DIFFRAXIS is not given
!
    else if (keyword(1:16) .eq. 'FMCSC_DIFFRAXON ') then
      read (wo1(wt11:wt12),*) diffr_axon(1)
      read (wo2(wt21:wt22),*) diffr_axon(2)
      read (wo3(wt31:wt32),*) diffr_axon(3)
!
!   input file for tabulated Bessel fxns
!
    else if (keyword(1:17) .eq. 'FMCSC_BESSELFILE ') then
      besselfile = wo1(wt11:wt12)
!
    end if
!
!
!
  else if (which.eq.2) then
!
  else if (which.eq.3) then
!
  else
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! covers: cutoffs, mcgrid
!
subroutine key_cutoffs(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 & rest,which,keyword)
!
  use iounit
  use cutoffs
  use keys
  use mcgrid
!
  implicit none
!
  character(MAXKEYLEN) wo1,wo2,wo3,rest
  character(MAXKWLEN) keyword
  integer which,choice,wt11,wt12,wt21,wt22,wt31,wt32,k
  RTYPE dum
!
  if (which.eq.1) then
!
!
!   choice of cutoff implementation: none (1), topology-assisted (4), or grid-based (3)
!
    if (keyword(1:16) .eq. 'FMCSC_CUTOFFMODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        use_cutoffs = .false.
        use_mcgrid = .false.
        use_rescrit = .false.
!      else if (choice.eq.2) then
!        use_cutoffs = .true.
!        use_mcgrid = .false.
!        use_rescrit = .false.
      else if (choice.eq.3) then
        use_cutoffs = .true.
        use_mcgrid = .true.
        use_rescrit = .false.
      else if (choice.eq.4) then
        use_cutoffs = .true.
        use_mcgrid = .false.
        use_rescrit = .true.
      else
        write(ilog,*) 'Invalid spec. for FMCSC_CUTOFFMODE. Defaulting to no cutoffs (probably undesired).'
        use_cutoffs = .false.
        use_mcgrid = .false.
        use_rescrit = .false.
      end if
!
!   printout for dynamics, energy check for MC
!
    else if (keyword(1:16) .eq. 'FMCSC_CHECKFREQ ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,nsancheck)
      if (nsancheck.lt.1) then
        nsancheck = 10000
        write(ilog,*) 'Invalid spec. for FMCSC_CHECKFREQ. Defaulting to ',nsancheck,'.'
      end if
!
!   whether to use a steric pre-screen for faster energy evaluation in MC moves 
!
    else if (keyword(1:16) .eq. 'FMCSC_USESCREEN ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice .eq. 1)  then
        use_stericscreen = .true. 
      end if
!
!   the height of the screen (also used as HS penalty) 
!
    else if (keyword(1:14) .eq. 'FMCSC_BARRIER ') then
      read (wo1(wt11:wt12),*) screenbarrier
      if (screenbarrier.lt.1000.0) then
        screenbarrier = 10000000.0
        write(ilog,*) 'Invalid spec. for FMCSC_BARRIER. Defaulting to ',screenbarrier,'.'
      end if
!
    end if
!
!
!
  else if (which.eq.2) then
!
!
!   cutoff settings (including long-range electrostatics base setting)
!
    if (use_cutoffs.EQV..true.) then
!
!     this is the general non-bonded cutoff for IPP and LJ
!
      if (keyword(1:15) .eq. 'FMCSC_NBCUTOFF ') then
        read (wo1(wt11:wt12),*) mcnb_cutoff
        if (mcnb_cutoff .le. 0.0d0)  then
          mcnb_cutoff = 10.0
          write(ilog,*) 'Invalid spec. for FMCSC_NBCUTOFF. Defaulting to ',mcnb_cutoff,'.'
        end if
        mcnb_cutoff2 = mcnb_cutoff**2
!
!     a cutoff for long-range attractions (1/r^3)
!
      else if (keyword(1:15) .eq. 'FMCSC_ELCUTOFF ') then
        read (wo1(wt11:wt12),*) mcel_cutoff
        if (mcel_cutoff .le. 0.0d0)  then
          mcel_cutoff = 20.0
          write(ilog,*) 'Invalid spec. for FMCSC_ELCUTOFF. Defaulting to ',mcel_cutoff,'.'
        end if
        mcel_cutoff2 = mcel_cutoff**2
        imcel2 = 1./mcel_cutoff2
!
!     how to deal with long-range electrostatics in MC
!     1) C-C and C-Dip fully, 2) just C-C fully
!     3) just C-C in single-site monopole resolution, 4) not at all
!
      else if (keyword(1:14) .eq. 'FMCSC_LREL_MC ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.le.0).OR.(choice.gt.4)) then
          lrel_mc = 3
          write(ilog,*) 'Invalid spec. for FMCSC_LREL_MC. Defaulting to ',lrel_mc,'.'
        else
          lrel_mc = choice
        end if
!
!     how to deal with long-range electrostatics in MD/BD/LD
!     1) truncation, 2) Ewald, 3) RF, 4) as third option in LREL_MC, 5) as first option in LREL_MC
!
      else if (keyword(1:14) .eq. 'FMCSC_LREL_MD ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.le.0).OR.(choice.gt.5)) then
          lrel_md = 1
          write(ilog,*) 'Invalid spec. for FMCSC_LREL_MD. Defaulting to ',lrel_md,'.'
        else
          lrel_md = choice
        end if
!
      else if (keyword(1:13) .eq. 'FMCSC_NBL_UP ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.10000)) then
          nbl_up = 10
          write(ilog,*) 'Invalid spec. for FMCSC_NBL_UP. Defaulting to ',nbl_up,'.'
        else
          nbl_up = choice
        end if
!
      end if
!
    end if !(end if for use_cutoffs)
!   
!   MC-grid spec.s
!
    if (use_mcgrid.EQV..true.) then
!
      if (keyword(1:14) .eq. 'FMCSC_GRIDDIM ') then
        read (wo1(wt11:wt12),*) grid%dim(1)
        read (wo2(wt21:wt22),*) grid%dim(2)
        read (wo3(wt31:wt32),*) grid%dim(3)
        do k=1,3
          if (grid%dim(k).le.1) then
            grid%dim(k) = 10
          end if
          grid%bu(k) = grid%dim(k)
        end do
!
      else if (keyword(1:18) .eq. 'FMCSC_GRIDMAXGPNB ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.lt.1) then
          grid%maxgpnbs = 125
          write(ilog,*) 'Invalid spec. for FMCSC_GRIDMAXGPNB. Defaulting to ',grid%maxgpnbs,'.'
        else
          grid%maxgpnbs = choice
        end if
!
      else if (keyword(1:18) .eq. 'FMCSC_GRIDMAXRSNB ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.lt.1) then
          grid%maxnbs = 200
          write(ilog,*) 'Invalid spec. for FMCSC_GRIDMAXRSNB. Defaulting to ',grid%maxnbs,'.'
        else
          grid%maxnbs = choice
        end if
!
!       whether to provide a summary of initial grid occupation statistics 
!
      else if (keyword(1:17) .eq. 'FMCSC_GRIDREPORT ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.1)  then
          grid%report = .true. 
        end if
!
      end if
!
    end if !(end if for use_mcgrid)

!
  else if (which.eq.3) then
!
  else
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! covers: paircorr, mcsums, contacts, dipolavg, clusters
!
subroutine key_analysis(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 & rest,which,keyword)
!
  use iounit
  use paircorr
  use keys
  use mcsums
  use contacts
  use dipolavg
  use system
  use molecule
  use clusters
  use pdb
  use fos
  use energies
!
  implicit none
!
  character(MAXKEYLEN) wo1,wo2,wo3,rest
  character(MAXKWLEN) keyword
  integer which,choice,wt11,wt12,wt21,wt22,wt31,wt32
  RTYPE dum
!
  if (which.eq.1) then
!
!   analysis group file
!
    if (keyword(1:16) .eq. 'FMCSC_ANGRPFILE ') then
!     we'll check for existence in the reader
      angrpfile = wo1(wt11:wt12)
!
!   output frequency for energy values
!
    else if (keyword(1:12) .eq. 'FMCSC_ENOUT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        enout = 500
        write(ilog,*) 'Invalid spec. for FMCSC_ENOUT. Defaulting to ',enout,'.'
      else
        enout = choice
      end if

!   output frequency for ensemble values (only written if using dynamics to sample)
!
    else if (keyword(1:13) .eq. 'FMCSC_ENSOUT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        ensout = 500
        write(ilog,*) 'Invalid spec. for FMCSC_ENSOUT. Defaulting to ',ensout,'.'
      else
        ensout = choice
      end if
!
!   output frequency for acceptance rates
!
    else if (keyword(1:13) .eq. 'FMCSC_ACCOUT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        accout = 2000
        write(ilog,*) 'Invalid spec. for FMCSC_ACCOUT. Defaulting to ',accout,'.'
      else
        accout = choice
      end if
!
! output frequency for xyz(pdb)-structure files
!
    else if (keyword(1:13) .eq. 'FMCSC_XYZOUT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        xyzout = 10000
        write(ilog,*) 'Invalid spec. for FMCSC_XYZOUT. Defaulting to ',xyzout,'.'
      else
        xyzout = choice
      end if


!
!   whether to write arc (1), pdb (2=default), dcd (3), xtc (4) (if LINK_XDR), or NetCDF (nc, if LINK_NETCDF)
!
    else if (keyword(1:13) .eq. 'FMCSC_XYZPDB ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
#ifdef LINK_XDR
#ifdef LINK_NETCDF
      if ((choice.le.0).OR.(choice.gt.5)) then
        write(ilog,*) 'Invalid spec. for FMCSC_XYZPDB. Defaulting to pdb-trajectory.'
        xyzmode = 2
      else
        xyzmode = choice
      end if
#else
      if ((choice.le.0).OR.(choice.gt.4)) then
        write(ilog,*) 'Invalid spec. for FMCSC_XYZPDB. Defaulting to pdb-trajectory.'
        xyzmode = 2
      else
        xyzmode = choice
      end if
#endif
#else
#ifdef LINK_NETCDF
      if ((choice.le.0).OR.((choice.gt.3).AND.(choice.ne.5))) then
        write(ilog,*) 'Invalid spec. for FMCSC_XYZPDB. Defaulting to pdb-trajectory.'
        xyzmode = 2
      else
        xyzmode = choice
      end if
#else
      if ((choice.le.0).OR.(choice.gt.3)) then
        write(ilog,*) 'Invalid spec. for FMCSC_XYZPDB. Defaulting to pdb-trajectory.'
        xyzmode = 2
      else
        xyzmode = choice
      end if
#endif
#endif
!
!   frequency of output of titration state (pH calculations)
!
    else if (keyword(1:12) .eq. 'FMCSC_PHOUT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice .le. 0) then
        phout = HUGE(phout)
        write(ilog,*) 'Invalid spec. for FMCSC_PHOUT. Analysis turned off (by default).'
      else
        phout = choice
      end if
!
!   calculation frequency for the contact analysis
!
    else if (keyword(1:18) .eq. 'FMCSC_CONTACTCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        contactcalc = 500
        write(ilog,*) 'Invalid spec. for FMCSC_CONTACTCALC. Defaulting to ',contactcalc,'.'
      else
        contactcalc = choice
      end if
!
!   calculation frequency for cluster analysis (from within contactcalc)
!
    else if (keyword(1:18) .eq. 'FMCSC_CLUSTERCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        clucalc = 100
        write(ilog,*) 'Invalid spec. for FMCSC_CLUSTERCALC. Defaulting to ',clucalc,'.'
      else
        clucalc = choice
      end if
!
!   offset for contact analysis, i.e., do not count sequence neighbors ...
!
    else if (keyword(1:17) .eq. 'FMCSC_CONTACTOFF ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.lt.0) then
        contact_off = 2
        write(ilog,*) 'Invalid spec. for FMCSC_CONTACTOFF. Defaulting to ',contact_off,'.'
      else
        contact_off = choice
      end if
!
!   cutoff (i.e., contact def.) for residue-c.o.m-based contact
!
    else if (keyword(1:17) .eq. 'FMCSC_CONTACTMIN ') then
      read (wo1(wt11:wt12),*) contact_cuts(1)
      if (contact_cuts(1).le.0.0) then
        contact_cuts(1) = 5.0
        write(ilog,*) 'Invalid spec. for FMCSC_CONTACTMIN. Defaulting to ',contact_cuts(1),'.'
      end if
      contact_cuts(1) = contact_cuts(1)**2
!
!   cutoff (i.e., contact def.) for residue-c.o.m-based contact
!
    else if (keyword(1:17) .eq. 'FMCSC_CONTACTCOM ') then
      read (wo1(wt11:wt12),*) contact_cuts(2)
      if (contact_cuts(2).le.0.0) then
        contact_cuts(2) = 5.0
        write(ilog,*) 'Invalid spec. for FMCSC_CONTACTCOM. Defaulting to ',contact_cuts(2),'.'
      end if
      contact_cuts(2) = contact_cuts(2)**2
!
!   calculation frequency for pair correlation analysis
!
    else if (keyword(1:13) .eq. 'FMCSC_PCCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        pccalc = 500
        write(ilog,*) 'Invalid spec. for FMCSC_PCCALC. Defaulting to ',pccalc,'.'
      else
        pccalc = choice
      end if
!
!   whether to obtain non-PBC distance distributions for polyamides
!
    else if (keyword(1:17) .eq. 'FMCSC_DO_AMIDEPC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        do_amidepc = .true.
      end if
!
!   bin-size for PC analysis
!
    else if (keyword(1:16) .eq. 'FMCSC_PCBINSIZE ') then
      read (wo1(wt11:wt12),*) pc_binsz
      if (pc_binsz.le.0.0) then
        pc_binsz = 0.2
        write(ilog,*) 'Invalid spec. for FMCSC_PCBINSIZE. Defaulting to ',pc_binsz,'.'
      end if
!
!   calculation frequency for instantaneous distance output
!
    else if (keyword(1:14) .eq. 'FMCSC_INSTGPC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.lt.0) then
        inst_gpc = 0
        write(ilog,*) 'Invalid spec. for FMCSC_INSTGPC. Instantaneous distance output turned off by default.'
      else
        inst_gpc = choice
      end if
!
!   calculation frequency for averaged dipole moments
!
    else if (keyword(1:14) .eq. 'FMCSC_DIPCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        write(ilog,*) 'Invalid spec. for FMCSC_DIPCALC. Analysis turned off (by default).'
        dipcalc = HUGE(dipcalc) 
      else
        dipcalc = choice
      end if
!
!   "calculation" frequency for aligning a read-in trajectory
!
    else if (keyword(1:16) .eq. 'FMCSC_ALIGNCALC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      align%yes = .true.
      if (choice.eq.0) then
        align%calc = 1
        write(ilog,*) 'Invalid spec. for FMCSC_ALIGNCALC. Defaulting to ',align%calc,'.'
      else if (choice.lt.0) then
        align%calc = -choice
        align%instrmsd = .true.
      else
        align%calc = choice
      end if
!
!   "collection" frequency for structural clustering 
!
    else if (keyword(1:15) .eq. 'FMCSC_CCOLLECT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        cstorecalc = HUGE(cstorecalc)
        write(ilog,*) 'Invalid spec. for FMCSC_CCOLLECT. Analysis turned off (by default)'
      else
        cstorecalc = choice
      end if
!
!   whether to align structures before applying RMSD distance metrics in clustering
!
    else if (keyword(1:13) .eq. 'FMCSC_CALIGN ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        align_for_clustering = .true.
      else
        align_for_clustering = .false.
      end if
!
!   what type of clustering algorithm to use
!
    else if (keyword(1:12) .eq. 'FMCSC_CMODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.le.0).OR.(choice.gt.5)) then
        cmode = 1
        write(ilog,*) 'Invalid spec. for FMCSC_CMODE. Defaulting to leader-based clustering (mode ',cmode,').'
      else
        cmode = choice
      end if
!
!   type of distance metric to use in structural clustering
!
    else if (keyword(1:16) .eq. 'FMCSC_CDISTANCE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.le.0).OR.(choice.gt.10)) then
        cdis_crit = 5
        write(ilog,*) 'Invalid spec. for FMCSC_CDISTANCE. Defaulting to ',cdis_crit,'.'
      else
        cdis_crit = choice
      end if
!
!   direction flag
!
    else if (keyword(1:14) .eq. 'FMCSC_CLEADER ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.le.0).OR.(choice.gt.4)) then
        cleadermode = 1
        write(ilog,*) 'Invalid spec. for FMCSC_CLEADER. Defaulting to forward/backward (mode ',cleadermode,').'
      else
        cleadermode = choice
      end if
!
!   whether to refine results
!
    else if (keyword(1:14) .eq. 'FMCSC_CREFINE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.0) then
        refine_clustering = .false.
      else
        refine_clustering = .true.
      end if
!
!   cutoff to use in neighbor list generation for advanced clustering/cFEP methods
!
    else if (keyword(1:14) .eq. 'FMCSC_CCUTOFF ') then
      read (wo1(wt11:wt12),*) chardcut
      if (chardcut.le.0.0) then
        chardcut = 3.0
        write(ilog,*) 'Invalid spec. for FMCSC_CCUTOFF. Defaulting to ',chardcut,'.'
      end if
!
!   linkage criterion for hierarchical clustering
!
    else if (keyword(1:14) .eq. 'FMCSC_CLINKAGE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.le.0).OR.(choice.ge.4)) then
        clinkage = 1 ! maximum
        write(ilog,*) 'Invalid spec. for FMCSC_CLINKAGE. Defaulting to maximum linkage.'
      else
        clinkage = choice
      end if
!
!   how many snapshots (roughly) are minimally required to call something a basin
!   used currently only for automatic starter generation for PROGIND method
!
    else if (keyword(1:16) .eq. 'FMCSC_CBASINMIN ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.3) then
        csivmin = 19
        write(ilog,*) 'Invalid spec. for FMCSC_CBASINMIN. Defaulting to ',csivmin,'.'
      else
        csivmin = choice
        if (mod(csivmin,2).eq.0) csivmin = csivmin + 1
      end if
!
!   how to compute sequence for one-shot processing (progress index)
!
    else if (keyword(1:19) .eq. 'FMCSC_CPROGINDMODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.ge.1).AND.(choice.le.2)) then
        cprogindex = choice
      else
        cprogindex = 2 
        write(ilog,*) 'Invalid spec. for FMCSC_CPROGINDMODE. Defaulting to approximate scheme.'
      end if
!
!   number of nearest neighbor guessings for progress index (CPROGINDMODE)
!
    else if (keyword(1:19) .eq. 'FMCSC_CPROGINDRMAX ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.gt.1) then
        cprogindrmax = choice
      else
        cprogindrmax = 50
        write(ilog,*) 'Invalid spec. for FMCSC_CPROGINDRMAX. Defaulting to ',cprogindrmax,'.'
      end if
!
!   A/B width (in numbers of snapshots) for ABC partitioning 
!
    else if (keyword(1:20) .eq. 'FMCSC_CPROGINDWIDTH ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.ge.2) then
        cprogpwidth = choice
      else
        cprogpwidth = 1000
        write(ilog,*) 'Invalid spec. for FMCSC_CPROGINDWIDTH. Defaulting to ',cprogpwidth,'.'
      end if
!
!   fold oder (outlier impact reducer)
!
    else if (keyword(1:19) .eq. 'FMCSC_CPROGMSTFOLD ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.ge.0) then
        cprogfold = choice
      else
        cprogfold = 0
        write(ilog,*) 'Invalid spec. for FMCSC_CPROGMSTFOLD. Defaulting to ',cprogfold,'.'
      end if
!
!   set cluster radius for applicable clustering methods
!
    else if (keyword(1:14) .eq. 'FMCSC_CRADIUS ') then
      read (wo1(wt11:wt12),*) cradius
      if (cradius.le.0.0) then
        cradius = 2.0 ! Angstrom for default RMSD distance
        write(ilog,*) 'Invalid spec. for FMCSC_CRADIUS. Defaulting to ',cradius,'.'
      end if
!
!   how many hierarchies to use in BIRCH-like clustering
!
    else if (keyword(1:18) .eq. 'FMCSC_BIRCHHEIGHT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.le.0) then
        c_nhier = 8
        write(ilog,*) 'Invalid spec. for FMCSC_BIRCHHEIGHT. Defaulting to ',c_nhier,'.'
      else
        c_nhier = choice
      end if
!
!   filename for custom clustering set requests
!
    else if (keyword(1:12) .eq. 'FMCSC_CFILE ') then
!     we'll check for existence in the reader
      cfilen = wo1(wt11:wt12)
!
!   filename for pre-generated NBL
!
    else if (keyword(1:14) .eq. 'FMCSC_NBLFILE ') then
!     we'll check for existence in the reader
      read_nbl_from_nc = .true.
      nblfilen = wo1(wt11:wt12)
!
!   filename for trajectory breaks
!
    else if (keyword(1:21) .eq. 'FMCSC_TRAJBREAKSFILE ') then
!     we'll check for existence in the reader
      tbrkfilen = wo1(wt11:wt12)
!
    end if
!
!
!
  else if (which.eq.2) then
!
!   output frequency for restart-files
!
    if (keyword(1:13) .eq. 'FMCSC_RSTOUT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (pdb_analyze.EQV..true.) then
        if ((choice.le.nsim).AND.(choice.gt.0)) then
          write(ilog,*) 'Warning. Disabling writing of restart files (requested by FMCSC_RSTOUT) for trajectory analysis runs.'
        end if
        rstout = nsim + 1
      else
        if (choice.le.0) then
          rstout = 10000
          write(ilog,*) 'Invalid spec. for FMCSC_RSTOUT. Defaulting to ',rstout,'.'
        else
          rstout = choice
        end if
      end if
    end if
!
    if (pccalc.le.nsim) then
!
!     input file for GPC analysis
!
      if (keyword(1:17) .eq. 'FMCSC_PCCODEFILE ') then
!       we'll check for existence in the reader
        pccodefile = wo1(wt11:wt12)
!
!     whether to write out a summary of the requested GPC analysis to a file
!
      else if (keyword(1:16) .eq. 'FMCSC_GPCREPORT ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.1) then
          gpc_report = .true.
        else
          gpc_report = .false.
        end if
!
      end if
!
    end if !(end if for pccalc<=nsim)
!
!
    if (savcalc.le.nsim) then
!
!     input file for atom-based SAV distributions
!
      if (keyword(1:18) .eq. 'FMCSC_SAVATOMFILE ') then
!       we'll check for existence in the reader
        savreq%filen = wo1(wt11:wt12)
!
!     whether to also print instantaneous SAV values for selected atoms
!
      else if (keyword(1:14) .eq. 'FMCSC_INSTSAV ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice*savcalc.gt.nsim)) then
          savreq%instfreq = nsim + 1
        else 
          savreq%instfreq = choice
        end if
!
      end if
!
    end if !(end if for savcalc<=nsim)
!
!
    if (cstorecalc.le.nsim) then
!
!     set top-level hierarchy radius for applicable clustering methods
!
      if (keyword(1:14) .eq. 'FMCSC_CMAXRAD ') then
        read (wo1(wt11:wt12),*) cmaxrad
        if (cmaxrad.le.0.0) then
          cmaxrad = 2.0*cradius
          write(ilog,*) 'Invalid spec. for FMCSC_CMAXRAD. Defaulting to ',cmaxrad,'.'
        else if ((cmaxrad.le.cradius).AND.(cmode.ne.3).AND.(cmode.ne.4)) then
          cmaxrad = 2.0*cradius ! Angstrom for default RMSD distance
          write(ilog,*) 'Invalid spec. for FMCSC_CMAXRAD. Defaulting to ',cmaxrad,'.'
        end if
!
!     for certain distance functions (CDISTANCE), this gives options for data preprocessing
!     centering, scaling, smoothing 
!
      else if (keyword(1:16) .eq. 'FMCSC_CPREPMODE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.ge.0).AND.(choice.le.5)) then
          cprepmode = choice
        else
          cprepmode = 0
          write(ilog,*) 'Invalid spec. for FMCSC_CPREPMODE. Defaulting to ',cprepmode,'.'
        end if
!
!     for certain distance functions (CDISTANCE), this gives option to alter the weights used by default
!
      else if (keyword(1:18) .eq. 'FMCSC_CMODWEIGHTS ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.ge.0).AND.(choice.le.9)) then
          cchangeweights = choice
        else
          cchangeweights = 0
          write(ilog,*) 'Invalid spec. for FMCSC_CMODWEIGHTS. Defaulting to ',cchangeweights,'.'
        end if
!
!     how to combine locally adaptive weights
!
      else if (keyword(1:20) .eq. 'FMCSC_CWCOMBINATION ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        cwcombination = choice
!
!     lag time for ACF weights at fixed lag time (in # of steps)
!
      else if (keyword(1:15) .eq. 'FMCSC_CLAGTIME ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.gt.0) then
          cwacftau = choice
        else
          cwacftau = 100
          write(ilog,*) 'Invalid spec. for FMCSC_CLAGTIME. Defaulting to ',cwacftau,'.'
        end if
!
!     window size for locally adaptive weights (in # of steps)
!
      else if (keyword(1:18) .eq. 'FMCSC_CWINDOWSIZE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.gt.0) then
          cwwindowsz = choice
        else
          cwwindowsz = 1000
          write(ilog,*) 'Invalid spec. for FMCSC_CWINDOWSIZE. Defaulting to ',cwwindowsz,'.'
        end if

!     when taking the inverse for some types of locally adaptive weights, buffer parameter to avoid infinity
!
      else if (keyword(1:16) .eq. 'FMCSC_CTRANSBUF ') then
        read (wo1(wt11:wt12),*) dum
        if (dum.gt.0.0) then
          cdynwbuf = dum
        else
          cdynwbuf = 1.0
          write(ilog,*) 'Invalid spec. for FMCSC_CTRANSBUF. Defaulting to ',cdynwbuf,'.'
        end if
!
!     order for B-splines for smoothing of clustering data
!
      else if (keyword(1:19) .eq. 'FMCSC_CSMOOTHORDER ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.gt.0) then
          csmoothie = choice
        else
          csmoothie = 10
          write(ilog,*) 'Invalid spec. for FMCSC_CSMOOTHORDER. Defaulting to ',csmoothie,'.'
        end if
!
!     in progress index method: which snapshot to start from (for 0, set automatically to multiple points,
!     for -1 with CPROGINDMODE = 2 set to center of largest cluster)
!
      else if (keyword(1:20) .eq. 'FMCSC_CPROGINDSTART ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.ge.0).OR.((choice.ge.-3).AND.(cmode.ne.4).OR.(cprogindex.eq.2))) then
          cprogindstart = choice
        else
          write(ilog,*) 'Invalid spec. for FMCSC_CPROGINDSTART. Defaulting to ',cprogindstart,'.'
          cprogindstart = 1
        end if
!
!     compute some sort of cut-based pseudo free energy profile 
!
      else if (keyword(1:15) .eq. 'FMCSC_CMSMCFEP ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.ge.1).AND.(choice.le.7)) then
          ccfepmode = choice
        else
          write(ilog,*) 'Invalid spec. for FMCSC_CMSMCFEP. Disabling (default).'
          ccfepmode = 0
        end if
!
!     treat SCCs independently in some analyses and reequilibrate separately 
!
      else if (keyword(1:19) .eq. 'FMCSC_CEQUILIBRATE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.1) then
          cequil = choice
        else if (cequil.ne.0) then
          write(ilog,*) 'Invalid spec. for FMCSC_CEQUILIBRATE. Disabling (default).'
          cequil = 0
        end if
!
!     parameter related to CBASINMIN in min-search
!
      else if (keyword(1:16) .eq. 'FMCSC_CBASINMAX ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.lt.csivmin) then
          csivmax = max(csivmin,49)
          write(ilog,*) 'Invalid spec. for FMCSC_CBASINMAX. Defaulting to ',csivmax,'.'
        else
          csivmax = choice
          if (mod(csivmax,2).eq.0) csivmax = csivmax + 1
        end if
!
#ifdef LINK_LAPACK
!
!     whether to do PCA and whether to print
!
      else if (keyword(1:14) .eq. 'FMCSC_PCAMODE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.le.0).OR.(choice.gt.3)) then
          pcamode = 1
          write(ilog,*) 'Invalid spec. for FMCSC_PCAMODE. Defaulting to ',pcamode,' (no PCA performed).'
        else
          pcamode = choice
        end if
!
!     PCA-based dimensionality reduction
!
      else if (keyword(1:17) .eq. 'FMCSC_CREDUCEDIM ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.le.0)) then
          reduced_dim_clustering = 0
          write(ilog,*) 'Invalid spec. for FMCSC_CREDUCEDIM. Performing original clustering.'
        else
          reduced_dim_clustering = choice
        end if
!
#endif
      end if
!
    end if ! (end if for cstorecalc<=nsim)
!
!
    if (xyzout.le.nsim) then
!
!     input file for atom-wise control over trajectory output
!
      if (keyword(1:18) .eq. 'FMCSC_TRAJIDXFILE ') then
!       we'll check for existence in the reader
        trajidxfile = wo1(wt11:wt12)
        use_trajidx = .true.
!
      end if
!
    end if !(end if for xyzout<=nsim)
!
!
    if (pdb_analyze.EQV..true.) then
!
!     input file for frame selection and weighting for trajectory analysis
!
      if (keyword(1:17) .eq. 'FMCSC_FRAMESFILE ') then
!       we'll check for existence in the reader
        frameidxfile = wo1(wt11:wt12)
        select_frames = .true.
!
      end if
!
    end if !(end if for pdb-analysis mode)
!
!
!
  else if (which.eq.3) then
!
    if (((align%calc.le.nsim).AND.(align%yes.EQV..true.)).OR.((cdis_crit.eq.6).AND.(align_for_clustering.EQV..true.))) then
!
!     input file for alignment procedure (defines alignment set, not structure)
!
      if (keyword(1:17) .eq. 'FMCSC_ALIGNFILE ') then
!       we'll check for existence in the reader
        align%filen = wo1(wt11:wt12)
!
      end if
!
    end if !(end if for aligncalc<=nsim)
!
  else
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! covers: sequen
!
subroutine key_sequen(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 & rest,which,keyword)
!
  use iounit
  use sequen
  use keys
  use pdb
  use system
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  character(MAXKEYLEN) wo1,wo2,wo3,rest
  character(MAXKWLEN) keyword
  integer which,choice,wt11,wt12,wt21,wt22,wt31,wt32
  RTYPE dum
  logical exists
!
  if (which.eq.1) then
!
!
!   initial structure randomization: 1) full, 2) RB only, 3) int only (requires pdb), 0) none (requires pdb or fyc)
!
    if (keyword(1:16) .eq. 'FMCSC_RANDOMIZE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.ge.0).OR.(choice.le.3))  then
        globrandomize = choice
      else
        globrandomize = 0
        write(ilog,*) 'Invalid spec. for FMCSC_RANDOMIZE. Defaulting to ',globrandomize,' &
 &(no structural randomization performed).'
      end if
!
!   a global maximum for attempts during initial structure randomization
!
    else if (keyword(1:17) .eq. 'FMCSC_RANDOMATTS ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.ge.3)  then
        globrdmatts = choice
      else
        globrdmatts = 300
        write(ilog,*) 'Invalid spec. for FMCSC_RANDOM_ATTS. Defaulting to ',globrdmatts,'.'
      end if
!
!   a universal IPP (and boundary) energy threshold for initial structure randomization
!
    else if (keyword(1:19) .eq. 'FMCSC_RANDOMTHRESH ') then
      read (wo1(wt11:wt12),*) globrdmthresh
      if ((globrdmthresh.le.1.0).OR.(globrdmthresh.gt.0.25*HUGE(dum))) then
        globrdmthresh = 10000.0
        write(ilog,*) 'Invalid spec. for FMCSC_RANDOMTHRESH. Defaulting to ',globrdmthresh,'.'
      end if
!
!   make the peptide chain cyclic (red herring)
!
    else if (keyword(1:13) .eq. 'FMCSC_CYCLIC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice .eq. 1)  then
        globcyclic = .true. 
      end if
!
!   sequence input file
!
    else if (keyword(1:14) .eq. 'FMCSC_SEQFILE ') then
      inquire(file=wo1(wt11:wt12),exist=exists)
      if (exists.EQV..false.) then
        write(ilog,*) 'Fatal. Cannot open spec.d sequence file (',&
 &wo1(wt11:wt12),') for FMCSC_SEQFILE.'
        call fexit()
      else
        seqfile = wo1(wt11:wt12)
      end if
!
!   whether to write a summary of the features of the input sequence
!
    else if (keyword(1:16) .eq. 'FMCSC_SEQREPORT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        seq_report = .true.
      else
        seq_report = .false.
      end if
!
!   FYC input file (deprecated)
!
    else if (keyword(1:14) .eq. 'FMCSC_FYCFILE ') then
      inquire(file=wo1(wt11:wt12),exist=exists)
      if (exists.EQV..true.) then
        fycinput = .true.
        fycfile = wo1(wt11:wt12)
      else
        write(ilog,*) 'Warning. Cannot open spec.d input file (',&
 &wo1(wt11:wt12),') for FMCSC_FYCFILE. Omitting request.'
        fycinput = .false.
      end if
!
    end if
!
!
!
  else if (which.eq.2) then
!
#ifdef ENABLE_MPI
!
    if (pdb_analyze.EQV..false.) then
!
!     whether to read several starting structure pdbs for certain parallel runs
!
      if (keyword(1:18) .eq. 'FMCSC_PDB_MPIMANY ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.1) then
          pdb_mpimany = .true.
        else
          pdb_mpimany = .false.
        end if
      end if
!
    end if
!
#endif
!
!
!
  else if (which.eq.3) then
!
!
!   spec's for pdb read-in-and-analyze or pdb starting structure mode
!
!
!   PDB input file
!
    if (keyword(1:14) .eq. 'FMCSC_PDBFILE ') then
      inquire(file=wo1(wt11:wt12),exist=exists)
      if (exists.EQV..true.) then
        pdbinput = .true.
        pdbinfile = wo1(wt11:wt12)
      else
#ifdef ENABLE_MPI
        if ((pdb_analyze.EQV..true.).OR.(pdb_mpimany.EQV..true.)) then
          pdbinfile = wo1(wt11:wt12)
          pdbinput = .true. ! check elsewhere
        else
          write(ilog,*) 'Warning. Cannot open spec.d input file (', &
 &wo1(wt11:wt12),') for FMCSC_PDBFILE. Omitting request.'
          pdbinput = .false.
        end if
#else
        write(ilog,*) 'Warning. Cannot open spec.d input file (', &
 &wo1(wt11:wt12),') for FMCSC_PDBFILE. Omitting request.'
        pdbinput = .false.
#endif
      end if
    end if
!
!
!
    if (pdb_analyze.EQV..true.) then
!
#ifdef LINK_XDR
!
!   XTC input file
!
      if (pdb_fileformat.eq.3) then
        if (keyword(1:14) .eq. 'FMCSC_XTCFILE ') then
          xtcinfile = wo1(wt11:wt12)
        end if
      end if
!
#endif
!
#ifdef LINK_NETCDF
!
!   NetCDF input file
!
      if (pdb_fileformat.eq.5) then
        if (keyword(1:17) .eq. 'FMCSC_NETCDFFILE ') then
          netcdfinfile = wo1(wt11:wt12)
        end if
      end if
!
#endif
!
!   DCD input file
!
      if (pdb_fileformat.eq.4) then
        if (keyword(1:14) .eq. 'FMCSC_DCDFILE ') then
          dcdinfile = wo1(wt11:wt12)
        end if
      end if
!
    end if !(end if for pdb_analyze)
!
!
!
  else
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! covers: pdb
!
subroutine key_pdb(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 & rest,which,keyword)
!
  use iounit
  use pdb
  use keys
  use system
  use sequen
  use clusters
!
  implicit none
!
  character(MAXKEYLEN) wo1,wo2,wo3,rest
  character(MAXKWLEN) keyword
  integer which,choice,wt11,wt12,wt21,wt22,wt31,wt32
  RTYPE dum
!
  if (which.eq.1) then
!
!
!   precision multiplier for XTC-files (the larger the higher precision)
!
#ifdef LINK_XDR
    if (keyword(1:14) .eq. 'FMCSC_XTCPREC ') then
      read (wo1(wt11:wt12),*) xtc_prec
      if (xtc_prec.le.100.0) then
        xtc_prec = 1000.0d0
        write(ilog,*) 'Invalid spec. for FMCSC_XTCPREC. Defaulting t&
 &o ',xtc_prec,'.'
      end if
    end if
#endif
!
!   mode: 1) use CAMPARI-conv. for O3* 2) use PDB-conv. for O3*
!
    if (keyword(1:18) .eq. 'FMCSC_PDB_NUCMODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,pdb_nucmode)
      if ((pdb_nucmode.gt.2).OR.(pdb_nucmode.lt.1)) then
        pdb_nucmode = 1
        write(ilog,*) 'Invalid spec. for FMCSC_PDB_NUCMODE. Defaulti&
 &ng to ',pdb_nucmode,'.'
      end if
!
!   mode: 1) use CAMPARI-conv. 2) use GROMACS-conv., 3) use CHARMM/NAMD-conv., 4) use AMBER conv.
!
    else if (keyword(1:17) .eq. 'FMCSC_PDB_W_CONV ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,pdb_convention(1))
      if ((pdb_convention(1).gt.4).OR.(pdb_convention(1).lt.1)) then
        pdb_convention(1) = 1
        write(ilog,*) 'Invalid spec. for FMCSC_PDB_W_CONV. Defaulti&
 &ng to ',pdb_convention(1),'.'
      end if
    else if (keyword(1:17) .eq. 'FMCSC_PDB_R_CONV ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,pdb_convention(2))
      if ((pdb_convention(2).gt.4).OR.(pdb_convention(2).lt.1)) then
        pdb_convention(2) = 1
        write(ilog,*) 'Invalid spec. for FMCSC_PDB_R_CONV. Defaulti&
 &ng to ',pdb_convention(2),'.'
      end if
!
!   whether to write to a single file or to a series of files
!
    else if (keyword(1:14) .eq. 'FMCSC_XYZMODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.le.0).OR.(choice.gt.2)) then
        pdb_writemode = 2
        write(ilog,*) 'Invalid spec. for FMCSC_XYZMODE. Defaulting to ',pdb_writemode,'.'
      else
        pdb_writemode = choice
      end if
!
!   tolerance settings for read-in structures
!
    else if (keyword(1:22) .eq. 'FMCSC_PDB_TOLERANCE_B ') then
      read (wo1(wt11:wt12),*) pdb_tolerance(1)
      read (wo2(wt21:wt22),*) pdb_tolerance(2)
      if ((pdb_tolerance(2).le.pdb_tolerance(1)).OR.(pdb_tolerance(1).lt.0.0).OR.(pdb_tolerance(2).le.0.0)) then
        pdb_tolerance(1) = 0.8
        pdb_tolerance(2) = 1.25
        write(ilog,*) 'Invalid spec.s for FMCSC_PDB_TOLERANCE_B. Defaulting t&
 &o ',100.0*pdb_tolerance(1),'% and ',100.0*pdb_tolerance(2),'%.'
      end if
!
    else if (keyword(1:22) .eq. 'FMCSC_PDB_TOLERANCE_A ') then
      read (wo1(wt11:wt12),*) dum
      if (dum.lt.0.0) then
        pdb_tolerance(3) = 20.0
        write(ilog,*) 'Invalid spec. for FMCSC_PDB_TOLERANCE_A. Defaulting t&
 &o ',pdb_tolerance(3),'.'
      else
        pdb_tolerance(3) = dum
      end if
!
!   whether to omit records for molecules tagged as solvent
!
    else if (keyword(1:18) .eq. 'FMCSC_XYZ_SOLVENT ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (dum.eq.1) then
        just_solutes = .false.
      else
        just_solutes = .true.
      end if
!
!   whether to break up molecules into portions of different images for PBC simulation trajectory output
!
    else if (keyword(1:19) .eq. 'FMCSC_XYZ_FORCEBOX ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (dum.eq.1) then
        pdb_force_box = .true.
      else
        pdb_force_box = .false.
      end if
!
!   pdb template (used for various things)
!
    else if (keyword(1:19) .eq. 'FMCSC_PDB_TEMPLATE ') then
      use_pdb_template = .true.
      pdbtmplfile = wo1(wt11:wt12)
!
!   mode: 1) build CAMPARI covalent geometry, extract only torsions; 2) extract all xyz, adopt PDB geometry
!
    else if (keyword(1:19) .eq. 'FMCSC_PDB_READMODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,pdb_readmode)
      if ((pdb_readmode.gt.2).OR.(pdb_readmode.lt.1)) then
        pdb_readmode = 1
        write(ilog,*) 'Invalid spec. for FMCSC_PDB_READMODE. Defau&
 &lting to ',pdb_readmode,'.'
      end if
!
    end if
!
!
!
  else if (which.eq.2) then
!
!
!   spec's for pdb read-in-and-analyze mode
!
    if (pdb_analyze.EQV..true.) then
!
!     format identifier (whether pdb input is a single pdb-trajectory (model) file or a series
!     of numbered pdb-files or an xtc-trajectory file or a dcd-trajectory file or a netcdf-trajectory file) 
!
      if (keyword(1:17) .eq. 'FMCSC_PDB_FORMAT ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,pdb_fileformat)
#ifdef LINK_XDR
#ifdef LINK_NETCDF
        if ((pdb_fileformat.gt.5).OR.(pdb_fileformat.lt.1)) then
          pdb_fileformat = 1
          write(ilog,*) 'Invalid spec. for FMCSC_PDB_FORMAT. Default&
 &ing to ',pdb_fileformat,'.'
        end if
#else
        if ((pdb_fileformat.gt.4).OR.(pdb_fileformat.lt.1)) then
          pdb_fileformat = 1
          write(ilog,*) 'Invalid spec. for FMCSC_PDB_FORMAT. Default&
 &ing to ',pdb_fileformat,'.'
        end if
#endif
#else
#ifdef LINK_NETCDF
        if ((pdb_fileformat.gt.5).OR.(pdb_fileformat.lt.1).OR.(pdb_fileformat.eq.3)) then
          pdb_fileformat = 1
          write(ilog,*) 'Invalid spec. for FMCSC_PDB_FORMAT. Default&
 &ing to ',pdb_fileformat,'.'
        end if
#else
        if ((pdb_fileformat.gt.4).OR.(pdb_fileformat.lt.1).OR.(pdb_fileformat.eq.3)) then
          pdb_fileformat = 1
          write(ilog,*) 'Invalid spec. for FMCSC_PDB_FORMAT. Default&
 &ing to ',pdb_fileformat,'.'
        end if
#endif
#endif
!
      end if
!
    end if !(end if for pdb_analyze)
!
!
!
  else if (which.eq.3) then
!
!   spec.s for starting a simulation from a pdb-structure
!
    if ((pdb_readmode.eq.2).OR.(pdb_analyze.EQV..true.))  then
!
!     mode: 1) never and 2) always re-build hydrogens
!       (since hydrogens are often ill-defined in structural input)
!
      if (keyword(1:16) .eq. 'FMCSC_PDB_HMODE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,pdb_hmode)
        if ((pdb_hmode.gt.2).OR.(pdb_hmode.lt.1)) then
          pdb_hmode = 1
          write(ilog,*) 'Invalid spec. for FMCSC_PDB_HMODE. Defaulti&
 &ng to ',pdb_hmode,'.'
        end if
!
      end if
!
    end if !(end if for pdb-stuff)
!
!
!
  else
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! covers: mpistuff and threads
!
subroutine key_mpi(wo1,wt11,wt12,wo2,wt21,wt22,wo3,wt31,wt32,&
 & rest,which,keyword)
!
  use iounit
  use keys
  use mpistuff
  use threads
!
  implicit none
!
  character(MAXKEYLEN) wo1,wo2,wo3,rest
  character(MAXKWLEN) keyword
  integer which,choice,wt11,wt12,wt21,wt22,wt31,wt32
  RTYPE dum
#ifdef ENABLE_MPI
  logical exists
#endif
!
  if (which.eq.1) then
!
!    flag, whether to use RE
!
    if (keyword(1:11) .eq. 'FMCSC_REMC ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice .eq. 1) then
        use_REMC = .true.
      else
        use_REMC = .false.
      end if
!
#ifdef ENABLE_MPI
!   flag, whether to use MPI averaging
!
    else if (keyword(1:13) .eq. 'FMCSC_MPIAVG ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice .eq. 1) then
        use_MPIAVG = .true.
      else
        use_MPIAVG = .false.
      end if
!
!   frequency for attempting RE (swap) moves
!
    else if (keyword(1:13) .eq. 'FMCSC_REFREQ ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,re_freq)
      if (re_freq .le. 1) then
        re_freq = 1000
        write(ilog,*) 'Invalid spec. for FMCSC_REFREQ. Defaulting &
 &to ',re_freq,'.'
      end if
!
!   whether or not to force exchange of xyz-coordinates in REMC (does NOT work in REMD)
    else if (keyword(1:17) .eq. 'FMCSC_REMC_DOXYZ ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        force_rexyz = .true.
      else
        force_rexyz = .false.
      end if
!
!   how to deal with velocities in REMD after successful exchange 
!   1) re-randomize to current T, 2) keep, 3) re-scale by T-offset (identical to 2) if no T-exchange)
!
    else if (keyword(1:17) .eq. 'FMCSC_RE_VELMODE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.le.0).OR.(choice.gt.3)) then
        re_velmode = 1
        write(ilog,*) 'Invalid spec. for FMCSC_RE_VELMODE. Defaulting to ',re_velmode,'.'
      else
        re_velmode = choice
      end if
!
!   whether or not to write out instantaneous trace of structure vs. condition
!
    else if (keyword(1:14) .eq. 'FMCSC_RETRACE ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.eq.1) then
        inst_retr = .true.
      else 
        inst_retr = .false.
      end if
!
!   choice for granularity in designing communication flow
!
    else if (keyword(1:19) .eq. 'FMCSC_MPIGRANULESZ ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice.lt.1) then
        mpi_granularity = 1
        write(ilog,*) 'Invalid spec. for FMCSC_MPIGRANULESZ. Defaulting to ',mpi_granularity,'.'
      else
        mpi_granularity = choice
      end if
!
!   flag, whether to use MPI standard collective communication routines in place of custom
!   ones
!
    else if (keyword(1:15) .eq. 'FMCSC_MPICOLLS ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice .eq. 1) then
        use_MPIcolls = .true.
      else
        use_MPIcolls = .false.
      end if
!
!   flag, whether to use reseeding MPIAVG feature
!
    else if (keyword(1:15) .eq. 'FMCSC_MPI_PIGS ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if (choice .eq. 1) then
        use_MPIMultiSeed = .true.
      else
        use_MPIMultiSeed = .false.
      end if
!
!   how many replicas to protect from reseeding for PIGS
!
    else if (keyword(1:19) .eq. 'FMCSC_MPI_GOODPIGS ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.lt.1).OR.(choice.gt.mpi_nodes)) then
        re_aux(7) = ceiling(0.5*mpi_nodes)
      else
        re_aux(7) = choice
      end if
!
#endif
!
#ifdef ENABLE_THREADS
!   max number of threads to allow 
!
    else if (keyword(1:16) .eq. 'FMCSC_NRTHREADS ') then
      read (wo1(wt11:wt12),*) dum
      call dum2int(dum,choice)
      if ((choice.lt.1).OR.(choice.gt.MAXTHREADS)) then
        thrdat%maxn = 2
        write(ilog,*) 'Invalid spec. for FMCSC_NRTHREADS. Defaulting &
 &to ',thrdat%maxn,'.'
      else
        thrdat%maxn = dum
      end if
#endif
!
    end if
!
!
!
  else if (which.eq.2) then
!
!   REMC settings
!
    if (use_REMC.EQV..true.) then
!
#ifdef ENABLE_MPI
!
!    the options for the parallel executable are conditionally compiled
!    (i.e., there's no need for a "i am parallel"-flag)
!
!     the RE input file
!
      if (keyword(1:13) .eq. 'FMCSC_REFILE ') then
        inquire(file=wo1(wt11:wt12),exist=exists)
        if (exists.EQV..true.) then
          re_infile = wo1(wt11:wt12)
        else
          write(ilog,*) 'Fatal. Cannot open spec.d input file (',wo1(wt11:wt12),') for FMCSC_REFILE.'
          call fexit()
        end if
!
!     in trajectory analysis, offer read-in of trace file for unREXing trajectories
!
      else if (keyword(1:16) .eq. 'FMCSC_TRACEFILE ') then
        re_traceinfile = wo1(wt11:wt12)
        re_aux(3) = 1
!
!     number of replicas to use (equals the number of nodes in the mpirun-command)
!     there has to be a corresponding input file with the appropriate number of conditions
!
      else if (keyword(1:15) .eq. 'FMCSC_REPLICAS ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,re_conditions)
        if (re_conditions .lt. 2) then
          re_conditions = 2
          write(ilog,*) 'Invalid spec. for FMCSC_REPLICAS. Defaulting to ',re_conditions,'.'
        end if
!
!     dimensionality of conditions (there has to be an appropriate number of columns in the
!     input file)
!
      else if (keyword(1:12) .eq. 'FMCSC_REDIM ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,re_conddim)
        if (re_conddim .le. 0) then
          re_conddim = 1
          write(ilog,*) 'Invalid spec. for FMCSC_REDIM. Defaulting to ',re_conddim,'.'
        end if
!
!     number of swap attempts during each exchange block
!
      else if (keyword(1:14) .eq. 'FMCSC_RESWAPS ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,re_tryswap)
        if (re_tryswap .le. 0) then
          re_tryswap = 10
          write(ilog,*) 'Invalid spec. for FMCSC_RESWAPS. Defaulting&
 & to ',re_tryswap,'.'
        end if
!
!     neighbor mode (whether to exchange only with (linear) neighbors
!     or with the whole set of conditions)
!     2) swap only neighbors       1) swap all (biased)
!
      else if (keyword(1:15) .eq. 'FMCSC_RENBMODE ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if ((choice.lt.1).OR.(choice.gt.2)) then
          re_nbmode = 2
          write(ilog,*) 'Invalid spec. for FMCSC_RENBMODE. Defaultin&
 &g to ',re_nbmode,'.'
        else
          re_nbmode = choice
        end if
!
!     how often to calculate overlap with other Hamiltonians
!
      else if (keyword(1:15) .eq. 'FMCSC_REOLCALC ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.le.0) then
          re_olcalc = 1000
          write(ilog,*) 'Invalid spec. for FMCSC_REOLCALC. Defaultin&
 &g to ',re_olcalc,'.'
        else
          re_olcalc = choice
        end if
!
!     whether or not to write out inst. values for overlap comput.
!
      else if (keyword(1:15) .eq. 'FMCSC_REOLINST ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.gt.0) then
          inst_reol = choice
        else if (choice.lt.0) then
          inst_reol = 0
          write(ilog,*) 'Invalid spec. for FMCSC_REOLINST. Defaultin&
 &g to ',inst_reol,'.'
        end if
!
!     whether or not to compute a full overlap matrix
!
      else if (keyword(1:14) .eq. 'FMCSC_REOLALL ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.1) then
          reol_all = .true.
        else
          reol_all = .false.
        end if
!
!     in MPI trajectory analysis, what trajectory output frequency to assume for unREXing
!      
      else if (keyword(1:17) .eq. 'FMCSC_RE_TRAJOUT ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.le.0) then
          re_aux(4) = 10000
          write(ilog,*) 'Invalid spec. for FMCSC_RE_TRAJOUT. Defaulting to ',re_aux(4),'.'
        else
          re_aux(4) = choice
        end if
!
!     in MPI trajectory analysis, what equilibration phase in trace file to assume for unREXing
!      
      else if (keyword(1:18) .eq. 'FMCSC_RE_TRAJSKIP ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.lt.0) then
          re_aux(5) = 10000
          write(ilog,*) 'Invalid spec. for FMCSC_RE_TRAJSKIP. Defaulting to ',re_aux(5),'.'
        else
          re_aux(5) = choice
        end if
!
      end if
!
#endif
!
    end if !(end if for use_REMC)
!
!   MPIAVG settings
!
    if (use_MPIAVG.EQV..true.) then
!
#ifdef ENABLE_MPI
!
!     whether to write individual or collective trajectory data
!
      if (keyword(1:17) .eq. 'FMCSC_MPIAVG_XYZ ') then
        read (wo1(wt11:wt12),*) dum
        call dum2int(dum,choice)
        if (choice.eq.1) then
          force_singlexyz = .true.
        else
          force_singlexyz = .false.
        end if
      end if
!
#endif
!
    end if !(end if for use_MPIAVG)
!
!
  else if (which.eq.3) then
!
  else
!
  end if
!
end
!
!-----------------------------------------------------------------------
