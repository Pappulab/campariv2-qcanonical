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
! CONTRIBUTIONS: Rohit Pappu, Hoang Tran                                   !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!-----------------------------------------------------------------------------
!
! this routine parses the sequence and then sets up the simulated molecules
!
!------------------------------------------------------------------------------
!
subroutine makepept()
!
  use fyoc
  use iounit
  use sequen
  use molecule
  use aminos
  use zmatrix
  use atoms
  use polypep
  use system
  use pdb
  use params ! Martin added 
  use inter ! Martin  : added 
  use martin_own ! Martin added 
  use energies ! Martin added
  use movesets
  use torsn
  
  implicit none
!
  integer rs,imol,jj,iomessage
  integer freeunit,iseqfile,aone,atwo,sickle,athree
  integer j,k,next,length,i,t1,t2,p1,p2
  RTYPE rr,random!,getbang,getztor,getblen
  character(60) resname
  character(MAXSTRLEN) seqlin
  logical done,startmol,foundN,foundC,foundM,foundit
  !integer, ALLOCATABLE :: temp_nrsin(:),temp_nrsnb(:)
  !integer, ALLOCATABLE :: temp_atin(:,:)
  integer temp_var1(2),temp
  RTYPE temp_var2
  integer mem(2)
  logical test
  integer limit
!  RTYPE  HS_bond_E(:) ! This is for special bond energy for HSQ
  
  
  
  integer alcsz,ntorsn_save ! Martin 
!
! a little bit of init
  aone = 1
  atwo = 2
  athree = 3
!
! read the parameter file
  call readprm()
  
!
! open sequence input file
  iseqfile = freeunit()
  call strlims(seqfile,t1,t2)
  open(unit=iseqfile,file=seqfile(t1:t2),status='old')
!
! first do a dry run to get the limits for the calculation
! also check for errors
!
! get the peptide sequence
  i = 0
  n_pdbunk = 0
  done = .false.
  startmol = .false.
  sickle = -1
  n_crosslinks = 0
  
  

  do while (done.EQV..false.)
!
     foundN = .false.
     foundC = .false.
     foundM = .false.
     if (sickle.ne.-1) then ! last res. was CYS
       jj = 0
       read(seqlin(next:MAXSTRLEN),49,iostat=iomessage) jj
       if ((iomessage.eq.0).AND.(jj.gt.0).AND.(jj.ne.i)) then
         if (i.lt.jj) n_crosslinks = n_crosslinks + 1
       end if
       do k=1,MAXSTRLEN
         seqlin(k:k) = ' '
       end do
     end if

     sickle = -1
     foundit = .false.
!
 49      format(i100)
 50      format(A)
!
     read(iseqfile,50,iostat=iomessage) seqlin
     if (iomessage.eq.-1) then
       close(unit=iseqfile)
       exit
     else if (iomessage.eq.2) then
       write(ilog,*) 'Fatal. File I/O error while processing sequenc&
 &e input file (got: ',seqfile(t1:t2),').'
       call fexit()
     end if
     call toupper(seqlin)
     next = 1
     call extract_str(seqlin,resname,next)
     call strlims(resname,p1,p2)
     length = p2-p1+1
!
!    first let's check for 3-letter-or-less codes
     if ((length.eq.3).OR.(length.eq.2)) then
        if (resname.eq.'END') then
          if (startmol.EQV..true.) then
            write(ilog,*) 'Fatal. Encountered end of sequence input in s&
 &equence file without end of previous molecule.'
            call fexit()
          end if
          done = .true.
          exit
        end if
        if (length.eq.2) then
          resname(3:3) = ' '
        end if
!here
        do j = 1,maxamino
          if(resname(1:3) .eq. amino(j)) then
            i = i + 1
            sickle = j
            foundit = .true.
          end if
        end do
        
!       N-terminal caps
        if ((foundit.EQV..true.).AND.((sickle.eq.27).OR.(sickle.eq.28).OR.&
 &          ((sickle.ge.76).AND.(sickle.le.87)))) then
          foundN = .true.
!       C-terminal caps
        else if ((foundit.EQV..true.).AND.((sickle.eq.29).OR.(sickle.eq.30))) then
          foundC = .true.
!       standard in-chain residues
        else if ((foundit.EQV..true.).AND.(((sickle.ge.1).AND.(sickle.le.25)).OR.&
 &               ((sickle.ge.31).AND.(sickle.le.35)).OR.&
 &               (sickle.eq.51).OR.&
 &               ((sickle.ge.64).AND.(sickle.le.75)).OR.&
 &               ((sickle.ge.108).AND.(sickle.le.115)).or.(sickle.ge.119).or.(sickle.eq.116))) then
          foundM = .true.
!       small molecules and atoms
        else if ((foundit.EQV..true.).AND.(((sickle.ge.36).AND.(sickle.le.50)).OR.&
 &               ((sickle.ge.52).AND.(sickle.le.63)).OR.&
 &               ((sickle.ge.88).AND.(sickle.le.107)).OR.&
 &               ((sickle.ge.117).AND.(sickle.le.118)))) then
          foundN = .true.  
          foundC = .true.
!       unknown residues (have to always be explicitly N/C labeled)
        else ! if (pdb_analyze.EQV..true.) then
          foundM = .true.
          i = i + 1
        end if
!
!       now see what we've got
        if ((foundN.EQV..true.).AND.(foundC.EQV..true.)) then
          if (startmol.EQV..true.) then
            write(ilog,*) 'Fatal. Found start of new molecule in seq&
 &uence file without end of previous molecule (',resname(1:3),').'
            call fexit()
          end if
          nmol = nmol + 1
        else if (foundN.EQV..true.) then
          if (startmol.EQV..true.) then
            write(ilog,*) 'Fatal. Found start of new molecule in seq&
 &uence file without end of previous molecule (',resname(1:3),').'
            call fexit()
          end if
          nmol = nmol + 1
          startmol = .true.
        else if (foundC.EQV..true.) then
          if (startmol.EQV..false.) then
            write(ilog,*) 'Fatal. Found end of current molecule in s&
 &equence file without ever starting it (',resname(1:3),').'
            call fexit()
          end if
          startmol = .false.
        else if (foundM.EQV..true.) then
          if (startmol.EQV..false.) then
            write(ilog,*) 'Fatal. Found intra-chain residue without &
 &ever starting a chain in seq. file (',resname(1:3),').'
            call fexit()
          end if
        else
          write(ilog,*) 'Fatal. Found an unknown residue in seq. fil&
 &e, which is currently not supported (',resname(1:3),').'
          call fexit()
        end if
!
!
!    the 5-letter codes are really 3-letter codes plus an appendix (_D, _C, _N)
     else if (length.eq.5) then
        do j = 1,maxamino
          if(resname(1:3) .eq. amino(j)) then
            i = i + 1
            sickle = j
            foundit = .true.
          end if
        end do
!
!       parse for appendix
        if ((resname(4:4).ne.'_').OR.((resname(5:5).ne.'D').AND.&
 &           (resname(5:5).ne.'N').AND.(resname(5:5).ne.'C'))) then
          write(ilog,*) 'Fatal. Found an unknown residue in seq. fil&
 &e, which is currently not supported (',resname(1:5),').'
          call fexit()
        end if
!
!       chirality request
        if (resname(5:5).eq.'D') then
!         D-aa in-chain residues
          if ((foundit.EQV..true.).AND.(((sickle.ge.2).AND.(sickle.le.23)).OR.&
 &            ((sickle.ge.31).AND.(sickle.le.35)).OR.((sickle.ge.108).AND.(sickle.le.115)).OR.&
 &            (sickle.eq.25).OR.(sickle.eq.51).or.(sickle.ge.119).or.(sickle.eq.116))) then
            foundM = .true.              
          end if
!
          if (foundM.EQV..true.) then
            if (startmol.EQV..false.) then
              write(ilog,*) 'Fatal. Found intra-chain residue withou&
 &t ever starting a chain in seq. file (',resname(1:5),').'
              call fexit()
            end if
          else
            write(ilog,*) 'Fatal. Found an unknown res. or incompreh&
 &ensible chirality request in seq. file (',resname(1:5),').'
            call fexit()
          end if
!
!       N-terminal residue request
        else if (resname(5:5).eq.'N') then
!         N-terminal standard residues
          if ((foundit.EQV..true.).AND.(((sickle.ge.1).AND.(sickle.le.25)).OR.&
 &            ((sickle.ge.31).AND.(sickle.le.35)).OR.&
 &            (sickle.eq.51).OR.&
 &            ((sickle.ge.64).AND.(sickle.le.75)).OR.((sickle.ge.108).AND.(sickle.le.115))&
 &            .or.(sickle.ge.119).or.(sickle.eq.116))) then
            foundN = .true.
!         unknown residues (have to always be explicitly N/C labeled)
          else if (foundit.EQV..false.) then
            foundN = .true.
            i = i + 1
!         allow redundant N-terminus spec.
          else if ((sickle.eq.27).OR.(sickle.eq.28).OR.((sickle.ge.76).AND.(sickle.le.87))) then
            foundN = .true.          
          end if
!
          if (foundN.EQV..true.) then
            if (startmol.EQV..true.) then
              write(ilog,*) 'Fatal. Found start of new molecule in s&
 &equence file without end of previous molecule (',resname(1:5),').'
              call fexit()
            end if
            nmol = nmol + 1
            startmol = .true.
          else
            write(ilog,*) 'Fatal. Found an unknown res. or incompreh&
 &ensible request for N-term. in seq. file (',resname(1:5),').'
            call fexit()
          end if
!
!       C-terminal residue request
        else if (resname(5:5).eq.'C') then
!         C-terminal standard residues
          if ((foundit.EQV..true.).AND.(((sickle.ge.1).AND.(sickle.le.25)).OR.&
 &            ((sickle.ge.31).AND.(sickle.le.35)).OR.&
 &            (sickle.eq.51).OR.(sickle.ge.119).or.(sickle.eq.116).OR.&
 &            ((sickle.ge.64).AND.(sickle.le.75)).OR.((sickle.ge.108).AND.(sickle.le.115)))) then
            foundC = .true.
!         unknown residues (have to always be explicitly N/C labeled)
          else if (foundit.EQV..false.) then
            foundC = .true.
            i = i + 1
          else if ((sickle.eq.29).OR.(sickle.eq.30)) then ! allow redundant C-terminus spec.
            foundC = .true.
          end if
!
          if (foundC.EQV..true.) then
            if (startmol.EQV..false.) then
              write(ilog,*) 'Fatal. Found end of current molecule in&
 & sequence file without ever starting it (',resname(1:5),').'
              call fexit()
            end if
            startmol = .false.
          else
            write(ilog,*) 'Fatal. Found an unknown res. or incompreh&
 &ensible request for C-term. in seq. file (',resname(1:5),').'
            call fexit()
          end if
!
        end if
!
!    the 7-letter codes are really 3-letter codes plus 2 appendices (_D, _C, _N)
!
     else if (length.eq.7) then
        do j = 1,maxamino
          if(resname(1:3) .eq. amino(j)) then
            i = i + 1
            sickle = j
            foundit = .true.
          end if
        end do
!
!       parse for appendix
        if ((resname(4:4).ne.'_').OR.(resname(6:6).ne.'_').OR.&
 &          ((resname(5:5).ne.'D').AND.&
 &           (resname(5:5).ne.'N').AND.(resname(5:5).ne.'C')).OR.&
 &          ((resname(7:7).ne.'D').AND.&
 &           (resname(7:7).ne.'N').AND.(resname(7:7).ne.'C'))) then
          write(ilog,*) 'Fatal. Found an unknown residue in seq. fil&
 &e, which is currently not supported (',resname(1:7),').'
          call fexit()
        end if
!
!       chirality request and N-term. request
        if (((resname(5:5).eq.'D').AND.(resname(7:7).eq.'N')).OR.&
 &          ((resname(5:5).eq.'N').AND.(resname(7:7).eq.'D'))) then
!         D-aa in-chain residues
          if ((foundit.EQV..true.).AND.(((sickle.ge.2).AND.(sickle.le.23)).OR.&
 &            ((sickle.ge.31).AND.(sickle.le.35)).OR.((sickle.ge.108).AND.(sickle.le.115)).OR.&
 &            (sickle.eq.25).OR.(sickle.eq.51).or.(sickle.ge.119).or.(sickle.eq.116))) then
            foundN = .true.              
          end if
!
          if (foundN.EQV..true.) then
            if (startmol.EQV..true.) then
              write(ilog,*) 'Fatal. Found start of new molecule in s&
 &equence file without end of previous molecule (',resname(1:7),').'
              call fexit()
            end if
            nmol = nmol + 1
            startmol = .true.
          else
            write(ilog,*) 'Fatal. Found an unknown res. or incompreh&
 &ensible N/D-request in seq. file (',resname(1:7),').'
            call fexit()
          end if
!
!       chirality request and C-term. request
        else if (((resname(5:5).eq.'D').AND.(resname(7:7).eq.'C'))&
 &       .OR.((resname(5:5).eq.'C').AND.(resname(7:7).eq.'D'))) then
!         C-terminal standard residues
          if ((foundit.EQV..true.).AND.(((sickle.ge.1).AND.(sickle.le.23)).OR.&
 &            ((sickle.ge.31).AND.(sickle.le.35)).OR.((sickle.ge.108).AND.(sickle.le.115)).OR.&
 &            (sickle.eq.51).OR.(sickle.eq.25).or.(sickle.ge.119).or.(sickle.eq.116))) then
            foundC = .true.              
          end if
!
          if (foundC.EQV..true.) then
            if (startmol.EQV..false.) then
              write(ilog,*) 'Fatal. Found end of current molecule in&
 & sequence file without ever starting it (',resname(1:7),').'
              call fexit()
            end if
            startmol = .false.
          else
            write(ilog,*) 'Fatal. Found an unknown res. or incompreh&
 &ensible C/D-request in seq. file (',resname(1:7),').'
            call fexit()
          end if
!
!       N/C-terminal (free amino acid) request
        else if (((resname(5:5).eq.'N').AND.(resname(7:7).eq.'C'))&
 &       .OR.((resname(5:5).eq.'C').AND.(resname(7:7).eq.'N'))) then
!         N/C-terminal standard residues
          if ((foundit.EQV..true.).AND.(((sickle.ge.1).AND.(sickle.le.25)).OR.&
 &            ((sickle.ge.31).AND.(sickle.le.35)).OR.((sickle.ge.108).AND.(sickle.le.115)).OR.&
 &            (sickle.eq.51).or.(sickle.ge.119).or.(sickle.eq.116).OR.&
 &            ((sickle.ge.64).AND.(sickle.le.75)))) then
            foundC = .true.
            foundN = .true.
!         unknown residues (have to always be explicitly N/C labeled)
          else if (foundit.EQV..false.) then
            foundN = .true.
            foundC = .true.
            i = i + 1
          end if
!
          if ((foundC.EQV..true.).AND.(foundN.EQV..true.)) then
            if (startmol.EQV..true.) then
              write(ilog,*) 'Fatal. Found start of new molecule in s&
 &equence file without end of previous molecule (',resname(1:7),').'
              call fexit()
            end if
            nmol = nmol + 1
          else
            write(ilog,*) 'Fatal. Found an unknown res. or incompreh&
 &ensible request for N/C-term. in seq. file (',resname(1:7),').'
            call fexit()
          end if
!
        end if
!
!    the 9-letter codes are really 3-letter codes plus all 3 appendices (_D, _C, _N)
!
!     else if (length.eq.9) then
!        do j = 1,maxamino
!          if(resname(1:3) .eq. amino(j)) then
!            seqtyp(i) = j
!            chiral(i) = 1
!          end if
!        end do
!
     end if
!
     if (sickle.eq.-1) then
       n_pdbunk = n_pdbunk + 1
     end if
!
  end do
!
  if (nmol.eq.0) then
    write(ilog,*) 'Fatal error. Requested calculation without any mo&
 &lecules. Sane, but empty sequence file?'
    call fexit()
  end if
!
  if (done.EQV..false.) then
    write(ilog,*) 'Fatal error. Please terminate sequence file with &
 &END line.'
    call fexit()
  end if
!
! we now have the full residue and molecule information
  nseq = i
!
!  write(*,*) 'allocating pdbunk'
  if (n_pdbunk.gt.0) allocate(pdb_unknowns(n_pdbunk))
!  write(*,*) 'allocating seq'
  call allocate_sequen(aone)
!  write(*,*) 'allocating mol.'
  call allocate_molecule(aone)
!  write(*,*) 'allocating fyoc'
  call allocate_fyoc(aone)
!  write(*,*) 'allocating done'
!
  rewind(unit=iseqfile)
!  
! get the peptide sequence
  i = 0
  n_pdbunk = 0
  nmol = 0
  done = .false.
  startmol = .false.
  n_crosslinks = 0
  seqtyp(:) = -1
  chiral(:) = 1
  
  do while (done.EQV..false.)
!
     if (i.gt.0) then
       jj = 0
       read(seqlin(next:MAXSTRLEN),49,iostat=iomessage) jj
       if ((iomessage.eq.0).AND.(jj.gt.0).AND.(jj.ne.i).AND.(jj.le.nseq)) then
         if (i.lt.jj) then
           n_crosslinks = n_crosslinks + 1
           crosslink(n_crosslinks)%rsnrs(1) = i
           crosslink(n_crosslinks)%rsnrs(2) = jj
         end if
       end if
       do k=1,MAXSTRLEN
         seqlin(k:k) = ' '
       end do
     end if

     foundN = .false.
     foundC = .false.
     foundM = .false.
     foundit = .false.
!
     read(iseqfile,50,iostat=iomessage) seqlin
     
     if (iomessage.eq.-1) then
       close(unit=iseqfile)
       exit
     else if (iomessage.eq.2) then
       write(ilog,*) 'Fatal. File I/O error while processing sequenc&
 &e input file (got: ',seqfile(t1:t2),').'
       call fexit()
     end if
     call toupper(seqlin)
     next = 1
     call extract_str(seqlin,resname,next)
     call strlims(resname,p1,p2)
     length = p2-p1+1
!
!    first let's check for N-terminal residues
     if ((length.eq.3).OR.(length.eq.2)) then
        if (resname.eq.'END') then
          done = .true.
          exit
        end if
        if (length.eq.2) then
          resname(3:3) = ' '
        end if
        do j = 1,maxamino
          if(resname(1:3) .eq. amino(j)) then
            i = i + 1
            
            seqtyp(i) = j
            chiral(i) = 1
            foundit = .true.
          end if
        end do
        if (foundit.EQV..true.) then
!         N-terminal caps
          if ((seqtyp(i).eq.27).OR.(seqtyp(i).eq.28).OR.&
 &                 ((seqtyp(i).ge.76).AND.(seqtyp(i).le.87))) then
            foundN = .true.
!         C-terminal caps
          else if ((seqtyp(i).eq.29).OR.(seqtyp(i).eq.30)) then
            foundC = .true.
!         standard in-chain residues
          else if (((seqtyp(i).ge.1).AND.(seqtyp(i).le.25)).OR.&
 &                 ((seqtyp(i).ge.31).AND.(seqtyp(i).le.35)).OR.&
 &                 (seqtyp(i).eq.51).OR.(seqtyp(i).eq.116).or.(seqtyp(i).ge.119).or.&
 &                 ((seqtyp(i).ge.64).AND.(seqtyp(i).le.75)).OR.&
 &                 ((seqtyp(i).ge.108).AND.(seqtyp(i).le.115))) then
            foundM = .true.
!         small molecules and atoms
          else if (((seqtyp(i).ge.36).AND.(seqtyp(i).le.50)).OR.&
 &                 ((seqtyp(i).ge.52).AND.(seqtyp(i).le.63)).OR.&
 &                 ((seqtyp(i).ge.88).AND.(seqtyp(i).le.107)).OR.&
 &                 ((seqtyp(i).ge.117).AND.(seqtyp(i).le.118))) then
            foundN = .true.  
            foundC = .true.
!          else
!            write(ilog,*) 'Fatal. Unable to parse residue class in makepept(). This is a bug.'
!            call fexit()
          end if
!       unknown residues (have to always be explicitly N/C labeled)
        else ! if (pdb_analyze.EQV..true.) then
          foundM = .true.
          i = i + 1
        end if
!
!       now see what we've got
        if ((foundN.EQV..true.).AND.(foundC.EQV..true.)) then
          nmol = nmol + 1
          rsmol(nmol,1) = i
          rsmol(nmol,2) = i
          phi(i) = 0.0d0
          psi(i) = 0.0d0
          omega(i) = 0.0d0
          do j = 1, 4
            chi(j,i) = 0.0d0
          end do
        else if (foundN.EQV..true.) then
          nmol = nmol + 1
          rsmol(nmol,1) = i
          phi(i) = 0.0d0
          psi(i) = 0.0d0
          omega(i) = 0.0d0
          do j = 1, 4
            chi(j,i) = 0.0d0
          end do
          startmol = .true.
        else if (foundC.EQV..true.) then
          rsmol(nmol,2) = i
          phi(i) = 0.0d0
          psi(i) = 0.0d0
          omega(i) = 0.0d0
          do j = 1, 4
            chi(j,i) = 0.0d0
          end do
          startmol = .false.
        else if (foundM.EQV..true.) then
          phi(i) = 0.0d0
          psi(i) = 0.0d0
          omega(i) = 0.0d0
          do j = 1, 4
            chi(j,i) = 0.0d0
          end do
        end if
!
!    the 5-letter codes are really 3-letter codes plus an appendix (_D, _C, _N)
     else if (length.eq.5) then
        do j = 1,maxamino
          if(resname(1:3) .eq. amino(j)) then
            i = i + 1
            seqtyp(i) = j
            chiral(i) = 1
            foundit = .true.
          end if
        end do
!
!       chirality request
        if (resname(5:5).eq.'D') then
!         D-aa in-chain residues
          if (foundit.EQV..true.) then
            if (((seqtyp(i).ge.2).AND.(seqtyp(i).le.23)).OR.&
 &              ((seqtyp(i).ge.31).AND.(seqtyp(i).le.35)).OR.&
 &              (seqtyp(i).eq.25).OR.(seqtyp(i).eq.51).OR.(seqtyp(i).eq.116).or.(seqtyp(i).ge.119).or.&
 &              ((seqtyp(i).ge.108).AND.(seqtyp(i).le.115))) then
              foundM = .true.              
            end if
          end if
!
          if (foundM.EQV..true.) then
            chiral(i) = -1
            phi(i) = 0.0d0
            psi(i) = 0.0d0
            omega(i) = 0.0d0
            do j = 1, 4
              chi(j,i) = 0.0d0
            end do
          end if
!
!       N-terminal residue request
        else if (resname(5:5).eq.'N') then
!         unknown residues (have to always be explicitly N/C labeled)
          if (foundit.EQV..false.) then
            foundN = .true.
            i = i + 1
!         N-terminal standard residues
          else if (((seqtyp(i).ge.1).AND.(seqtyp(i).le.25)).OR.&
 &              ((seqtyp(i).ge.31).AND.(seqtyp(i).le.35)).OR.&
 &              (seqtyp(i).eq.51).OR.(seqtyp(i).eq.116).or.(seqtyp(i).ge.119).or.&
 &              ((seqtyp(i).ge.64).AND.(seqtyp(i).le.75)).OR.&
 &              ((seqtyp(i).ge.108).AND.(seqtyp(i).le.115))) then
            foundN = .true.
!         allow redundant N-terminus spec.
          else if ((seqtyp(i).eq.27).OR.(seqtyp(i).eq.28).OR.((seqtyp(i).ge.76).AND.(seqtyp(i).le.87))) then
            foundN = .true.
          end if
!
          if (foundN.EQV..true.) then
            nmol = nmol + 1
            rsmol(nmol,1) = i
            phi(i) = 0.0d0
            psi(i) = 0.0d0
            omega(i) = 0.0d0
            do j = 1, 4
              chi(j,i) = 0.0d0
            end do
            moltermid(nmol,1) = 1
            startmol = .true.
          end if
!
!       C-terminal residue request
        else if (resname(5:5).eq.'C') then
!         unknown residues (have to always be explicitly N/C labeled)
          if (foundit.EQV..false.) then
            foundC = .true.
            i = i + 1
!         C-terminal standard residues
          else if (((seqtyp(i).ge.1).AND.(seqtyp(i).le.25)).OR.&
 &             ((seqtyp(i).ge.31).AND.(seqtyp(i).le.35)).OR.&
 &             (seqtyp(i).eq.51).OR.(seqtyp(i).eq.116).or.(seqtyp(i).ge.119).OR.&
 &             ((seqtyp(i).ge.64).AND.(seqtyp(i).le.75)).OR.&
 &             ((seqtyp(i).ge.108).AND.(seqtyp(i).le.115))) then
            foundC = .true.
          else if ((seqtyp(i).eq.29).OR.(seqtyp(i).eq.30)) then ! allow redundant C-terminus spec.
            foundC = .true.
          end if
!
          if (foundC.EQV..true.) then
            rsmol(nmol,2) = i
            phi(i) = 0.0d0
            psi(i) = 0.0d0
            omega(i) = 0.0d0
            do j = 1, 4
              chi(j,i) = 0.0d0
            end do
            moltermid(nmol,2) = 1
            startmol = .false.
          end if
!
        end if
!
!    the 7-letter codes are really 3-letter codes plus 2 appendices (_D, _C, _N)
!
     else if (length.eq.7) then
        do j = 1,maxamino
          if(resname(1:3) .eq. amino(j)) then
            i = i + 1
            seqtyp(i) = j
            chiral(i) = 1
            foundit = .true.
          end if
        end do
!
!       chirality request and N-term. request
        if (((resname(5:5).eq.'D').AND.(resname(7:7).eq.'N')).OR.&
 &          ((resname(5:5).eq.'N').AND.(resname(7:7).eq.'D'))) then
!         D-aa in-chain residues
          if (foundit.EQV..true.) then
            if (((seqtyp(i).ge.2).AND.(seqtyp(i).le.23)).OR.&
 &              ((seqtyp(i).ge.31).AND.(seqtyp(i).le.35)).OR.&
 &              (seqtyp(i).eq.25).OR.(seqtyp(i).eq.51).OR.(seqtyp(i).eq.116).or.(seqtyp(i).ge.119).OR.&
 &              ((seqtyp(i).ge.108).AND.(seqtyp(i).le.115))) then
              foundN = .true.              
            end if
          end if
!
          if (foundN.EQV..true.) then
            nmol = nmol + 1
            rsmol(nmol,1) = i
            chiral(i) = -1
            phi(i) = 0.0d0
            psi(i) = 0.0d0
            omega(i) = 0.0d0
            do j = 1, 4
              chi(j,i) = 0.0d0
            end do
            moltermid(nmol,1) = 1
            startmol = .true.
          end if
!
!       chirality request and C-term. request
        else if (((resname(5:5).eq.'D').AND.(resname(7:7).eq.'C'))&
 &       .OR.((resname(5:5).eq.'C').AND.(resname(7:7).eq.'D'))) then
!         C-terminal standard residues
          if (foundit.EQV..true.) then
            if (((seqtyp(i).ge.1).AND.(seqtyp(i).le.23)).OR.&
 &              ((seqtyp(i).ge.31).AND.(seqtyp(i).le.35)).OR.&
 &              (seqtyp(i).eq.51).OR.(seqtyp(i).eq.25).OR.(seqtyp(i).eq.116).or.(seqtyp(i).ge.119).OR.&
 &              ((seqtyp(i).ge.108).AND.(seqtyp(i).le.115))) then
              foundC = .true.              
            end if
          end if
!
          if (foundC.EQV..true.) then
            rsmol(nmol,2) = i
            chiral(i) = -1
            phi(i) = 0.0d0
            psi(i) = 0.0d0
            omega(i) = 0.0d0
            do j = 1, 4
              chi(j,i) = 0.0d0
            end do
            moltermid(nmol,2) = 1
            startmol = .false.
          end if
!
!       N/C-terminal (free amino acid) request
        else if (((resname(5:5).eq.'N').AND.(resname(7:7).eq.'C'))&
 &       .OR.((resname(5:5).eq.'C').AND.(resname(7:7).eq.'N'))) then
!         unknown residues (have to always be explicitly N/C labeled)
          if (foundit.EQV..false.) then
            foundN = .true.
            foundC = .true.
            i = i + 1 
!         N/C-terminal standard residues
          else if (((seqtyp(i).ge.1).AND.(seqtyp(i).le.25)).OR.&
 &              ((seqtyp(i).ge.31).AND.(seqtyp(i).le.35)).OR.&
 &              (seqtyp(i).eq.51).OR.(seqtyp(i).eq.116).or.(seqtyp(i).ge.119).OR.&
 &              ((seqtyp(i).ge.64).AND.(seqtyp(i).le.75)).OR.&
 &              ((seqtyp(i).ge.108).AND.(seqtyp(i).le.115))) then
            foundC = .true.
            foundN = .true.
          end if
!
          if ((foundC.EQV..true.).AND.(foundN.EQV..true.)) then
            nmol = nmol + 1
            rsmol(nmol,1) = i
            rsmol(nmol,2) = i
            phi(i) = 0.0d0
            psi(i) = 0.0d0
            omega(i) = 0.0d0
            do j = 1, 4
              chi(j,i) = 0.0d0
            end do
            moltermid(nmol,1) = 1
            moltermid(nmol,2) = 1
          end if
!
        end if
!
!    the 9-letter codes are really 3-letter codes plus all 3 appendices (_D, _C, _N)
!
!     else if (length.eq.9) then
!        do j = 1,maxamino
!          if(resname(1:3) .eq. amino(j)) then
!            seqtyp(i) = j
!            chiral(i) = 1
!          end if
!        end do
!
     end if
!
     if (seqtyp(i).eq.-1) then
       n_pdbunk = n_pdbunk + 1
       pdb_unknowns(n_pdbunk)(1:3) = resname(1:3)
     end if
!
  end do
!
  if (done.EQV..false.) then
    write(ilog,*) 'Fatal error. Please terminate sequence file with &
 &END line.'
    call fexit()
  end if
!
  nseq = i
! for any unsupported residues, grab as much info as possible from PDB
  if (n_pdbunk.gt.0) then
    call infer_from_pdb(1)
  else if ((use_pdb_template.EQV..true.).AND.(pdb_analyze.EQV..false.)) then
    write(ilog,*) 'Warning. Ignoring superfluous request for use of pdb-template (no unsupported residues &
 &and not in trajectory analysis mode).'
    use_pdb_template = .false.
  end if
!
! now let's initialize some more arrays
  call setup_resrad()
!  write(*,*) 'allocating polypep.'
  call allocate_polypep(aone)
!  write(*,*) 'allocating atomestls.'
  call allocate_atomestls(aone)
!  write(*,*) 'allocating done'
  if (n_crosslinks.gt.0) then
    call setup_crosslinks()
  end if
!
  do imol=1,nmol
    do rs = rsmol(imol,1),rsmol(imol,2)
      if (((seqtyp(rs).ge.1).AND.(seqtyp(rs).le.25)).OR.&
 &        ((seqtyp(rs).ge.31).AND.(seqtyp(rs).le.35)).OR.(seqtyp(rs).eq.116).or.(seqtyp(rs).ge.119).OR.&
 &        (seqtyp(rs).eq.51).OR.((seqtyp(rs).ge.108).AND.(seqtyp(rs).le.115))) then
        phi(rs) = 179.5d0!-57.0
        psi(rs) = 179.5d0!-47.0
        omega(rs) = 180.0d0
        if (ua_model.gt.0) then
          if (rs.eq.rsmol(imol,1)+1) then
            if (seqtyp(rs-1).eq.28) then
              omega(rs) = 0.0
            end if
          end if
        end if
        do jj = 1,MAXCHI
          chi(jj,rs) = 180.0d0
!         override for P
          if (seqtyp(rs).eq.9) chi(jj,rs) = 0.0d0
        end do
!       chi3-override for Q/E
        if (seqtyp(rs).eq.19) chi(3,rs) = 90.0d0
        if (seqtyp(rs).eq.18) chi(3,rs) = 90.0d0
        if (seqtyp(rs).eq.108) chi(3,rs) = 90.0d0
        if (seqtyp(rs).eq.113) chi(3,rs) = 90.0d0
        if (seqtyp(rs).eq.116) chi(3,rs) = 120.0d0
!       the hack of initializing pucker states randomly is made defunct (rr is always <= 1.0)
!       the chi-hijack still works to communicate the choice of pucker state to
!       sidechain(...)
!       do not change unless you understand the consequences
        if(seqtyp(rs).eq.9) then
          if (rs.gt.rsmol(imol,1)) then
            rr = 1.0 ! random()
            if (rr.le.1.5) then
! exo (C-gamma far from Ci; phi +/-60.8)
              if (chiral(rs).eq.1) then
                chi(1,rs) = -28.0d0
                chi(2,rs) =  4.3d0  ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = -60.8d0
              else
                chi(1,rs) = 30.5d0  ! 28.0d0
                chi(2,rs) = -10.0d0 ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = 60.8d0
              end if
            else
! endo (C-gamma closer to Ci; phi -/+71.0)
              if (chiral(rs).eq.1) then
                chi(1,rs) = 31.0d0   ! 32.0d0
                chi(2,rs) =  -11.9d0 ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = -71.0d0
              else
                chi(1,rs) = -29.0d0  ! -31.0d0
                chi(2,rs) = 6.4d0 ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = 71.0d0
              end if
            end if
          else
            rr = 1.0 ! random()
            if (rr.le.1.5) then
              if (chiral(rs).eq.1) then
                chi(1,rs) = -28.0d0
                chi(2,rs) =  4.3d0  ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = -0.8d0
              else
                chi(1,rs) = 30.5d0  ! 28.0d0
                chi(2,rs) = -10.0d0 ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = 120.8d0
              end if
            else
              if (chiral(rs).eq.1) then
                chi(1,rs) = 31.0d0   ! 32.0d0
                chi(2,rs) =  -11.9d0 ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = -11.0d0
              else
                chi(1,rs) = -29.0d0  ! -31.0d0
                chi(2,rs) = 6.4d0 ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = 131.0d0
              end if
            end if
          end if
!       dealing with hydroxyproline
        else if (seqtyp(rs).eq.32) then
          phi(rs) = -60.0d0
          psi(rs) = 110.0d0
          chi(1,rs) = -29.4d0
          chi(2,rs) =  39.6d0 
        end if
!       some small molecules have residual degrees of freedom
      else if (((seqtyp(rs).ge.41).AND.(seqtyp(rs).le.44)).OR.&
 &             (seqtyp(rs).eq.48).OR.(seqtyp(rs).eq.50).OR.&
 &             ((seqtyp(rs).ge.88).AND.(seqtyp(rs).le.99))) then
        do jj = 1,MAXCHI
          chi(jj,rs) = 180.0d0
        end do
!       override for PPA
        if (seqtyp(rs).eq.44) then
          if (random().gt.0.5) then
            chi(1,rs) = 120.0d0
          else
            chi(1,rs) = -120.0d0
          end if
        end if
      else if ((seqtyp(rs).ge.29).AND.(seqtyp(rs).le.30)) then
        omega(rs) = 180.0d0
      else if ((seqtyp(rs).ge.41).AND.(seqtyp(rs).le.42)) then
        omega(rs) = 0.0d0
!      nucleic acids have all kinds of degrees of freedom
      else if ((seqtyp(rs).ge.64).AND.(seqtyp(rs).le.87)) then
        do jj = 1,6
          nucs(jj,rs) = 180.0d0
        end do
        do jj = 1,MAXCHI
          chi(jj,rs) = 180.0d0
        end do
!       sugar pucker missing
      end if
    end do
  end do
!     
! set up polypeptide chain geometry
  call proteus_init(aone)
  do imol=1,nmol
   call proteus(imol)
  end do
  call proteus_init(atwo)
!
  if (n_pdbunk.gt.0) call infer_from_pdb(2)
!
! top level patches first -> biotypes (these potentially reset everything beyond topology)
  call read_biotpatchfile()
!
! potential mass patches a.s.a.p.
  call read_atpatchfile(athree)
!


   allocate(spec_limit(nseq))
   spec_limit(:)=0
  if ((do_hsq.eqv..true.).or.(do_hs.eqv..true.)) then 
    nhis=0
    do rs=1,nseq
        if ((amino(seqtyp(rs)).eq."HIX").and.(do_hsq.eqv..true.)) then
            nhis=nhis+1
        else if ((do_hsq.eqv..true.).and.((amino(seqtyp(rs)).eq."HDX").or.(amino(seqtyp(rs)).eq."HEX"))) then 
            print *,"FATAL : You need HIX in the sequence file if you are going to use Histidine in the protonation walk (HSQ))"
            call fexit()  

!        else if ((do_pka.eqv..true.).and.(amino(seqtyp(rs)).eq."HDX")) then 
!            print *,"WARNING : It is recommended to use the HEX intead of HDX, but this is not fatal"
!        else if ((amino(seqtyp(rs)).eq."HIX").and.(do_pka.eqv..true.)) then
!            print *,"FATAL : The use of HIX is reserved for the protonation walk. Please change your sequence to HEX"
!            call fexit()  
        end if 
    end do
    if (nhis.ne.0) then 
        allocate(his_state(nhis,2))
        allocate(his_state_save(nhis,2))
        nhis=0
        do rs=1,nseq
            if ((amino(seqtyp(rs)).eq."HIX")) then 
                nhis=nhis+1
                his_state(nhis,1)=rs
                his_state(nhis,2)=0

    !            if (do_hsq.eqv..true.) then ! if it is HSQ, we are setting the tautomer to the one we want
    !                his_state(nhis,2)=0 !zero or one, can be directly addedd to the limit number
    !            else if (amino(seqtyp(rs)).eq."HDX") then 
    !                his_state(nhis,2)=1
    !            else if (amino(seqtyp(rs)).eq."HEX") then 
    !                his_state(nhis,2)=0
    !            end if 
            end if 
        end do 
    end if 
  end if 
  

  if (par_top.eq.1) then
    call replace_bond_values()
    if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then 
        call get_alternative_coord()
    end if 
    
  end if 
  
! from the fully assembled Z-matrix, build backbone coordinates
! and populate 1-2-topology arrays
  
  call makexyz2()

  call setup_connectivity_2()

! now regenerate internals from the newly found Cartesian coordinates
! this serves to remove intrinsic incompatibilities in the Z-matrix
! (like ill-defined bond angles) which were adjusted by genxyz(...) called
! from within makexyz2()
! also set up 1-3 and 1-4 topology arrays

  do imol=1,nmol
    call genzmat(imol)
  end do

  call makexyz2()

  call setup_connectivity_34()

! potential changes to LJ parameters must occur early due to reliance of other settings on them
  call read_ljpatchfile()

! now allocate memory for atomic parameters and assign their values
!  write(*,*) 'allocating atomprms.'
  call allocate_atomprms(aone)

!  write(*,*) 'done'
!
! check atomic valencies
  call valence_check()

! some high-level residue-wise parameters
  call absinth_residue()
  
! populate short-range interaction  and more atomic arrays: note that
! setup_srinter has to be called before setup_savol (within absinth_atom)
! can be called!
!  write(*,*) 'allocating inter.'
  call allocate_inter(aone)
!  write(*,*) 'done'

  if ((do_pka_3.eqv..true.).or.(do_hsq.eqv..true.)) then 
      
 
    call setup_pka_prm() ! Martin : the goal of this is that the correct epsilon values for the limits are set BEFORE setupt_inter is 
    
    ! Martin : next bit is to implement histidine three limits case

    allocate(fudge_limits(nseq,3))
    allocate(iaa_limits(nseq,3))
    allocate(temp_nrsin(nseq,3))
    allocate(temp_nrsnb(nseq,3))
    allocate(temp_atin(nseq,2,3))
    
    ! So Now I want all the histidine to all be in the HE state for the neutral (that will be the second limit)
    if ((nhis.ne.0).and.(do_hsq.eqv..true.)) then 
        call qcann_setup_srinter(3)
    end if 
    call qcann_setup_srinter(2)
    call qcann_setup_srinter(1)

    if ((nhis.ne.0).and.(do_hsq.eqv..true.)) then 
        temp=3 
    else 
        temp=2  
    end if  

    do limit=1,temp
        do i=1,nseq             
        !    This part is to make sure that future iteration goes through the entire array in all cases
            if (nrsintra(i).lt.temp_nrsin(i,limit)) then 
                nrsintra(i)=temp_nrsin(i,limit)
            end if 
            if (nrsnb(i).lt.temp_nrsnb(i,limit)) then 
                nrsnb(i)=temp_nrsnb(i,limit)
            end if 
        end do
    end do 
    
    if ((nhis.ne.0).and.(do_hsq.eqv..true.)) then 
        call align_and_extend_arrays1(1,3)
    end if 

    ! it may be necessary to 
    call align_and_extend_arrays1(1,2)

    deallocate(temp_nrsin)
    deallocate(temp_nrsnb)
    
    !!!!!! Next !!!!!!
    
    do rs=1,nseq
        do i=1,nrsnb(rs)
            iaa(rs)%atnb(i,:)=iaa_limits(rs,1)%atnb(i,:)
        end do 
    end do 
    
    do rs=1,nseq
        do i=1,nrsintra(rs)
            iaa(rs)%atin(i,:)=iaa_limits(rs,1)%atin(i,:)
        end do 
    end do 
    do rs=1,nseq
        do i=1,nrsnb(rs)
            if ((iaa_limits(rs,1)%atnb(i,1).ne.iaa_limits(rs,2)%atnb(i,1)).or.&
            &(iaa_limits(rs,1)%atnb(i,2).ne.iaa_limits(rs,2)%atnb(i,2)).or.&
            &(iaa_limits(rs,1)%atnb(i,1).ne.iaa(rs)%atnb(i,1)).or.&
            &(iaa_limits(rs,1)%atnb(i,2).ne.iaa(rs)%atnb(i,2))) then
                print *,"nb",iaa_limits(rs,1)%atnb(i,:),iaa_limits(rs,2)%atnb(i,:),iaa(rs)%atnb(i,:)
            end if 
        end do 
    end do 
    do rs=1,nseq
        do i=1,nrsintra(rs)
            if ((iaa_limits(rs,1)%atin(i,1).ne.iaa_limits(rs,2)%atin(i,1)).or.&
            &(iaa_limits(rs,1)%atin(i,2).ne.iaa_limits(rs,2)%atin(i,2)).or.&
            &(iaa_limits(rs,1)%atin(i,1).ne.iaa(rs)%atin(i,1)).or.&
            &(iaa_limits(rs,1)%atin(i,2).ne.iaa(rs)%atin(i,2))) then
                print *,"intra",iaa_limits(rs,1)%atin(i,:),iaa_limits(rs,2)%atin(i,:),iaa(rs)%atin(i,:)
            end if 
        end do 
    end do 
!
    if (nhis.ne.0) then
        call qcann_absinth_atom(3) 
    end if 
    
    call qcann_absinth_atom(2) 
    call qcann_absinth_atom(1) 

    ! And now retrieve all the recently formed components

    do i=1,n
        do j=1,3
            if (atr_limits(i,j).eq.0) then 
                atstv_limits(i,j)=0.0
            end if 
        end do 
    end do 

    !end new part ())
  else   
    call setup_srinter()! Martin : sets up the fudge factors 
    call absinth_atom() 
  end if 
    
  
  
  if (par_top.eq.2) then 
    ! Note that this is not needed for the histidines since the bond type do not change in between the tautomers
    call rewrite_params()
  end if 
  
  
  

  
! Martin : made major change, but that should be the easiest : getting the same workflow as for the dermination of the other limits
  !call absinth_atom() ! Martin : Seting up of the pka params here 

  if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then 
      !That part is to get rid of the E-310 zeros
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            do i=1,size(iaa(rs)%par_impt(:,1))
                do j=1,size(iaa(rs)%par_impt(i,:))
                        iaa_limits(rs,1)%par_impt(i,j)=0.
                        iaa_limits(rs,2)%par_impt(i,j)=0.
                end do 
            end do 
            do i=1,size(iaa(rs)%par_di(:,1))
                do j=1,size(iaa(rs)%par_di(i,:))
                    iaa_limits(rs,1)%par_di(i,j)=0.
                    iaa_limits(rs,2)%par_di(i,j)=0.
                end do 
            end do 
        end do 
    end do  
    
    if (do_pka_2.eqv..true.) then ! Why is that only for pKa ? 
        call scale_stuff_pka()
    end if 

    
!    if (do_hsq.eqv..true.) then
!        call scale_stuff_hsq()
!    end if 
  end if 

  call makexyz2()
! assign residue- and molecule-wise parameters
  
  call parse_sequence()
  call absinth_molecule()


  

  if (n_crosslinks.gt.0) then
    call setup_crosslinks2()
  end if

  do imol=1,nmol
    call find_rotlsts(imol)
    call parse_rotlsts(imol)
  end do
  call trans_rotlstdof()

  call find_ntorsn()
!  
!  call find_ntorsn()
! forces requires find_ntorsn() to be done
!  write(*,*) 'allocating forces.'
  call allocate_forces(aone)
!  write(*,*) 'done'
  call correct_srinter()
  call correct_sampler()
!
! populate molecular RBC-arrays
  do imol=1,nmol
    call update_rigid(imol)
  end do
!
! set some analysis flags and check water loops (note that BT/LJ patches have already been applied)
  call parse_moltyp()
!
! make a back up of Z matrix
  bangpr(1:n) = bang(1:n)
  blenpr(1:n) = blen(1:n)

!  do imol=1,nmol! Goes through all molecules
!        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
!            print *,"Residue",rs
!            do i=1,size(iaa(rs)%par_di(:,1))
!                print *,iaa(rs)%par_di(i,:)
!
!            end do 
!        end do
!  end do 

 
end
!
!----------------------------------------------------------------------
!
subroutine setup_crosslinks()
!
  use sequen
  use aminos
  use iounit
  use fyoc
  use molecule
!
  implicit none
!
  integer i,j,rs1,rs2,imol
  logical samemol
!
  do i=1,n_crosslinks
    do j=i+1,n_crosslinks
      if ((crosslink(j)%rsnrs(1).eq.crosslink(i)%rsnrs(1)).OR.&
 &        (crosslink(j)%rsnrs(1).eq.crosslink(i)%rsnrs(2)).OR.&
 &        (crosslink(j)%rsnrs(2).eq.crosslink(i)%rsnrs(1)).OR.&
 &        (crosslink(j)%rsnrs(2).eq.crosslink(i)%rsnrs(2))) then
        write(ilog,*) 'Fatal. Only one type of chemical cross-link&
 & is possible for a given residue. Check sequence input.'
        call fexit()
      end if
    end do
    rs1 = crosslink(i)%rsnrs(1)
    rs2 = crosslink(i)%rsnrs(2)
    if ((seqtyp(rs1).eq.8).AND.(seqtyp(rs2).eq.8)) then
      samemol = .true.
      do imol=1,nmol
        if ((rsmol(imol,1).le.rs1).AND.(rsmol(imol,2).ge.rs1)) then
          if ((rsmol(imol,1).gt.rs2).OR.(rsmol(imol,2).lt.rs2)) then
            samemol = .false.
            exit
          end if
        end if
      end do
      if (samemol.EQV..true.) then
        crosslink(i)%itstype = 1 ! intramol. disulfide
      else
        crosslink(i)%itstype = 2 ! intermol. disulfide
      end if
      disulf(rs1) = rs2
      disulf(rs2) = rs1
    else
      write(ilog,*) 'Fatal. Chemical cross-links are currently not supported&
 & between residues ',amino(seqtyp(rs1)),' and ',amino(seqtyp(rs2)),'. Check back later.'
      call fexit()
    end if
  end do
!
end
!
!------------------------------------------------------------------------------------------
!
subroutine setup_crosslinks2()
!
  use sequen
  use iounit
  use molecule
!
  implicit none
!
  integer i,j,k,l,kk,ends,cnt,nlks,imol,jmol,last,first
  integer, ALLOCATABLE:: lklst(:,:), counts(:), idx(:), counts2(:)
  logical, ALLOCATABLE:: donewith(:)
  logical notdone
!
  allocate(lklst(nseq,2))
  allocate(counts(nmol))
  allocate(counts2(nmol))
  allocate(donewith(nseq))
  allocate(idx(n_crosslinks+3))
!
! try to arrange intermolecular crosslinks such that they can be observed sequentially (inside-out)
  donewith(:) = .false.
  counts(:) = 0
!
  do imol=1,nmol
    do j=1,n_crosslinks
      if (donewith(j).EQV..true.) cycle
      if (molofrs(crosslink(j)%rsnrs(1)).ne.molofrs(crosslink(j)%rsnrs(2))) then
        if ((molofrs(crosslink(j)%rsnrs(1)).eq.imol).OR.(molofrs(crosslink(j)%rsnrs(2)).eq.imol)) then
          counts(imol) = counts(imol) + 1
        end if
      end if
    end do
  end do
!
  idx(:) = 0
  first = 1
  last = n_crosslinks
  do while (first.le.last)
!
    counts2(:) = counts(:)
    do j=1,n_crosslinks
      imol = molofrs(crosslink(j)%rsnrs(1))
      jmol = molofrs(crosslink(j)%rsnrs(2))
      if (donewith(j).EQV..true.) cycle
      if (imol.ne.jmol) then
        if ((counts(imol).eq.1).OR.(counts(jmol).eq.1)) then
          idx(j) = last
          last = last - 1
          counts2(imol) = counts(imol) - 1 
          counts2(jmol) = counts(jmol) - 1
          if (counts(jmol).gt.1) then
            kk = crosslink(j)%rsnrs(1)
            crosslink(j)%rsnrs(1) = crosslink(j)%rsnrs(2)
            crosslink(j)%rsnrs(2) = kk
          end if
          donewith(j) = .true.
        end if
      else
        idx(j) = first
        first = first + 1
        donewith(j) = .true.
      end if
    end do
 !
    counts(:) = counts2(:)
  end do
!
  do j=1,n_crosslinks
    if (j.eq.idx(j)) cycle
    do i=1,n_crosslinks
      if (idx(i).eq.j) exit
    end do
    do k=1,2
      kk = crosslink(j)%rsnrs(k)
      crosslink(j)%rsnrs(k) = crosslink(i)%rsnrs(k)
      crosslink(i)%rsnrs(k) = kk
    end do
    kk = crosslink(j)%itstype
    crosslink(j)%itstype = crosslink(i)%itstype
    crosslink(i)%itstype = kk
    kk = idx(j)
    idx(j) = idx(i)
    idx(i) = kk
  end do
!
! set coupling sets for crosslinks
  nlks = 0
  donewith(:) = .false.
  do i=1,n_crosslinks
    if (donewith(i).EQV..true.) cycle
    if (molofrs(crosslink(i)%rsnrs(2)).ne.molofrs(crosslink(i)%rsnrs(1))) then
      nlks = 1
      lklst(1,1) = i
      notdone = .true.
      imol = molofrs(crosslink(i)%rsnrs(1))
      jmol = molofrs(crosslink(i)%rsnrs(2))
      counts(1) = imol
      counts(2) = jmol
      donewith(i) = .true.
      cnt = 2
      ends = 2
      do while (notdone.EQV..true.)
        notdone = .false.
        do j=1,n_crosslinks
          if (donewith(j).EQV..true.) cycle
          if (molofrs(crosslink(j)%rsnrs(2)).ne.molofrs(crosslink(j)%rsnrs(1))) then
            do k=1,cnt
              if (molofrs(crosslink(j)%rsnrs(2)).eq.counts(k)) then
                donewith(j) = .true.
                notdone = .true.
                kk = 0
                do l=1,cnt
                  if (l.eq.k) cycle
                  if (counts(l).eq.molofrs(crosslink(j)%rsnrs(1))) then
                    ends = ends - 2
                    kk = 1
                    write(ilog,*) 'Fatal. Intermolecular crosslinks which produce&
 & ring topologies are currently not supported. Please check back later.'
                    call fexit()
                    exit
                  end if
                end do
                if (kk.eq.0) then
                  cnt = cnt + 1
                  counts(cnt) = molofrs(crosslink(j)%rsnrs(1)) 
                  nlks = nlks + 1
                  lklst(nlks,1) = j
                  exit
                end if
              else if (molofrs(crosslink(j)%rsnrs(1)).eq.counts(k)) then
                donewith(j) = .true.
                notdone = .true.
                kk = 0
                do l=1,cnt
                  if (l.eq.k) cycle
                  if (counts(l).eq.molofrs(crosslink(j)%rsnrs(2))) then
                    ends = ends - 2
                    kk = 1
                    write(ilog,*) 'Fatal. Intermolecular crosslinks which produce&
 & ring topologies are currently not supported. Please check back later.'
                    call fexit()
                    exit
                  end if
                end do
                if (kk.eq.0) then
                  cnt = cnt + 1
                  counts(cnt) = molofrs(crosslink(j)%rsnrs(2))
                  nlks = nlks + 1
                  lklst(nlks,1) = j
                  exit
                end if
              end if
            end do
          end if
        end do
        if (ends.eq.0) notdone = .false.
      end do
      do k=1,nlks
        crosslink(lklst(k,1))%nolks = nlks-1
        kk = 0
        do l=1,nlks
          if (l.eq.k) cycle
          kk = kk + 1
          crosslink(lklst(k,1))%olks(kk) = lklst(l,1)
        end do
!        write(ilog,*) 'Found ',crosslink(k)%olks(1:crosslink(k)%nolks),' for ',lklst(k,1)
      end do
    end if
  end do
!
  deallocate(lklst)
  deallocate(donewith)
  deallocate(counts)
  deallocate(counts2)
  deallocate(idx)
!
  do i=1,n_crosslinks
    crlk_idx(crosslink(i)%rsnrs(1)) = i
    crlk_idx(crosslink(i)%rsnrs(2)) = i
  end do
  call crosslink_excludes()
!
end
!
!--------------------------------------------------------------------------
!
! a subroutine to parse the sequence (molecule-wise) to set some analysis
! flags and to allow pooling of the data for identical molecules in some
! cases
!
subroutine parse_sequence()
!
  use iounit
  use sequen
  use molecule
  use polypep
  use pdb
  use atoms
  use fyoc
!
  implicit none
!
  integer imol,rs,jmol,k,kkk,aone,cl
  integer, ALLOCATABLE:: dummy(:)
  logical newmt
!
  allocate(dummy(nmol))
!
! first a dry run to establish nmoltyp
  nmoltyp = 1
  dummy(1) = 1
!
  do imol=2,nmol
    newmt = .true.
    do jmol=1,nmoltyp
      if ((rsmol(dummy(jmol),2)-rsmol(dummy(jmol),1)).eq.&
 &        (rsmol(imol,2)-rsmol(imol,1))) then
        newmt = .false.
        do rs=1,rsmol(imol,2)-rsmol(imol,1)+1
          if (seqtyp(rs+rsmol(imol,1)-1).ne.&
 &            seqtyp(rs+rsmol(dummy(jmol),1)-1)) then
            newmt = .true.
            exit
          end if
          if (chiral(rs+rsmol(imol,1)-1).ne.&
 &            chiral(rs+rsmol(dummy(jmol),1)-1)) then
            newmt = .true.
            exit
          end if
          if (seqtyp(rs+rsmol(imol,1)-1).eq.26) then ! always declare new for unknowns unless biotype seq is identical
            cl = at(rs+rsmol(dummy(jmol),1)-1)%bb(1)
            if ((at(rs+rsmol(dummy(jmol),1)-1)%nbb.ne.at(rs+rsmol(imol,1)-1)%nbb).OR.&
 &              (at(rs+rsmol(dummy(jmol),1)-1)%nsc.ne.at(rs+rsmol(imol,1)-1)%nsc)) then
              newmt = .true.
              exit
            end if
            do k=at(rs+rsmol(imol,1)-1)%bb(1),at(rs+rsmol(imol,1)-1)%bb(1)+at(rs+rsmol(imol,1)-1)%nbb+at(rs+rsmol(imol,1)-1)%nsc-1
              if (b_type(k).ne.b_type(cl)) then
                newmt = .true.
                exit
              else
                cl = cl + 1
              end if
            end do
            if (newmt.EQV..true.) exit
          end if
        end do
        if (newmt.EQV..false.) then
!         always declare a new molecule type for crosslinked molecules (this may be overkill but so what)
          do cl=1,n_crosslinks
            if (((crosslink(cl)%rsnrs(1).ge.rsmol(imol,1)).AND.(crosslink(cl)%rsnrs(1).le.rsmol(imol,2))).OR.&
 &((crosslink(cl)%rsnrs(1).ge.rsmol(dummy(jmol),1)).AND.(crosslink(cl)%rsnrs(1).le.rsmol(dummy(jmol),2)))) then
              newmt = .true.
            end if
            if (((crosslink(cl)%rsnrs(2).ge.rsmol(imol,1)).AND.(crosslink(cl)%rsnrs(2).le.rsmol(imol,2))).OR.&
 &((crosslink(cl)%rsnrs(2).ge.rsmol(dummy(jmol),1)).AND.(crosslink(cl)%rsnrs(2).le.rsmol(dummy(jmol),2)))) then
              newmt = .true.
            end if
          end do
        end if
        if (newmt.EQV..false.) exit
      end if
    end do
    if (newmt.EQV..true.) then
      nmoltyp = nmoltyp + 1
      dummy(nmoltyp) = imol
    end if
  end do
!
! now allocate memory
  aone = 1
  call allocate_molecule_type(aone)
!
! the first molecule defines the first type and adds one count to that
  kkk = 1
  moltyp(1,1) = 1
  moltyp(1,2) = 1
  moltypid(1) = 1
!
  do imol=2,nmol
!
    newmt = .true.
!
    do jmol=1,kkk
      if ((rsmol(moltyp(jmol,1),2)-rsmol(moltyp(jmol,1),1)).eq.&
 &        (rsmol(imol,2)-rsmol(imol,1))) then
        newmt = .false.
        do rs=1,rsmol(imol,2)-rsmol(imol,1)+1
          if (seqtyp(rs+rsmol(imol,1)-1).ne.&
 &            seqtyp(rs+rsmol(moltyp(jmol,1),1)-1)) then
            newmt = .true.
            exit
          end if
          if (chiral(rs+rsmol(imol,1)-1).ne.&
 &            chiral(rs+rsmol(moltyp(jmol,1),1)-1)) then
            newmt = .true.
            exit
          end if
          if (seqtyp(rs+rsmol(imol,1)-1).eq.26) then ! always declare new for unknowns
            cl = at(rs+rsmol(moltyp(jmol,1),1)-1)%bb(1)
            if ((at(rs+rsmol(moltyp(jmol,1),1)-1)%nbb.ne.at(rs+rsmol(imol,1)-1)%nbb).OR.&
 &              (at(rs+rsmol(moltyp(jmol,1),1)-1)%nsc.ne.at(rs+rsmol(imol,1)-1)%nsc)) then
              newmt = .true.
              exit
            end if
            do k=at(rs+rsmol(imol,1)-1)%bb(1),at(rs+rsmol(imol,1)-1)%bb(1)+at(rs+rsmol(imol,1)-1)%nbb+at(rs+rsmol(imol,1)-1)%nsc-1
              if (b_type(k).ne.b_type(cl)) then
                newmt = .true.
                exit
              else
                cl = cl + 1
              end if
            end do
            if (newmt.EQV..true.) exit
          end if
        end do
      end if
!
      if (newmt.EQV..false.) then
!       always declare a new molecule type for crosslinked molecules (this may be overkill but so what)
        do cl=1,n_crosslinks
          if (((crosslink(cl)%rsnrs(1).ge.rsmol(imol,1)).AND.(crosslink(cl)%rsnrs(1).le.rsmol(imol,2))).OR.&
 &((crosslink(cl)%rsnrs(1).ge.rsmol(dummy(jmol),1)).AND.(crosslink(cl)%rsnrs(1).le.rsmol(dummy(jmol),2)))) then
            newmt = .true.
          end if
          if (((crosslink(cl)%rsnrs(2).ge.rsmol(imol,1)).AND.(crosslink(cl)%rsnrs(2).le.rsmol(imol,2))).OR.&
 &((crosslink(cl)%rsnrs(2).ge.rsmol(dummy(jmol),1)).AND.(crosslink(cl)%rsnrs(2).le.rsmol(dummy(jmol),2)))) then
            newmt = .true.
          end if
        end do
      end if
!
      if (newmt.EQV..false.) then
        moltyp(jmol,2) = moltyp(jmol,2) + 1
        moltypid(imol) = jmol
        exit
      end if
    end do
! 
    if (newmt.EQV..true.) then
      kkk = kkk + 1
      moltyp(kkk,1) = imol
      moltyp(kkk,2) = 1
      moltypid(imol) = kkk
    end if
!
  end do
!
! set total number of rigid-body degrees of freedom, assign analysis groups (default), and set solvent/solute
  totrbd = 0
  nressolute = 0
  nsolutes = 0
  nangrps = nmoltyp
  do imol=1,nmol
    an_grp_mol(imol) = moltypid(imol)
    if ((rsmol(imol,2)-rsmol(imol,1)).eq.0) is_solvent(imol) = .true.
    if (atmol(imol,2)-atmol(imol,1).eq.0) then 
      totrbd = totrbd + 3
    else if (atmol(imol,2)-atmol(imol,1).eq.1) then
      totrbd = totrbd + 5
    else
      totrbd = totrbd + 6
    end if
  end do
!
! read (and eventually overwrite analysis group spec.s)
  call read_analysisgrpfile()
  call alloc_angrps(an_grp_mol,is_solvent)
  if (just_solutes.EQV..true.) then
    pdbeffn = 0
    do imol=1,nmol
      if (is_solvent(imol).EQV..false.) pdbeffn = pdbeffn + (atmol(imol,2)-atmol(imol,1)+1)
    end do
  else
    pdbeffn = n
  end if
!
  deallocate(dummy)
!
end
!
!---------------------------------------------------------------------------
!
! this routine sets analysis flags and checks for water loop support
! also provides sequence report (formerly in parse_sequence)
!
subroutine parse_moltyp()
!
  use iounit
  use sequen
  use molecule
  use atoms
  use polypep
  use zmatrix
  use cutoffs
  use params
  use aminos
!
  implicit none
!
  integer jmol,ati,rs,k,j,rs2,ii
  RTYPE, ALLOCATABLE:: epsik(:)
  logical dropout,foundit
!
! fix resrad for polynucleotides with 5' phosphate
  do rs=1,nseq
    if ((seqflag(rs).eq.22).AND.(rs.eq.rsmol(molofrs(rs),1)).AND.(amino(seqtyp(rs)).ne.'R5P').AND.&
 &      (amino(seqtyp(rs)).ne.'D5P')) then
      resrad(rs) = resrad(rs) + 0.8
    end if
  end do
!
! go through the molecule types and analyze what can be done with them
  do jmol=1,nmoltyp
!
    do ati=atmol(moltyp(jmol,1),1)+1,atmol(moltyp(jmol,1),2)
      if (izrot(ati)%alsz.gt.0) then
        do_pol(jmol) = .true.
        exit
      end if
    end do
!
    if (rsmol(moltyp(jmol,1),2).gt.rsmol(moltyp(jmol,1),1)) then
!     do_pers requires two specific reference atoms per residues (caps are allowed to misbehave) and
!     at least two eligible residues
      do_pers(jmol) = .true.
      k = 0
      do rs=rsmol(moltyp(jmol,1),1),rsmol(moltyp(jmol,1),2)
        if ((ci(rs).gt.0).AND.(ni(rs).gt.0)) then
          k = k + 1
        else if ((nuci(rs,2).gt.0).AND.(nuci(rs,6).gt.0)) then
          k = k + 1
        else 
          if ((rs.gt.rsmol(moltyp(jmol,1),1)).AND.(rs.lt.rsmol(moltyp(jmol,1),2))) then
            do_pers(jmol) = .false.
            exit
          end if
        end if
      end do
      if (k.lt.2) do_pers(jmol) = .false.
!     terminal residues are never suited for traditional polypeptide torsional
!     analysis as they have at most one informative backbone angle
!     note that this has to change should one do generic torsional analysis
      do rs=rsmol(moltyp(jmol,1),1)+1,rsmol(moltyp(jmol,1),2)-1
        if (seqpolty(rs).eq.'P') then
          do_tors(jmol) = .true.
          exit
        end if
      end do
    end if
  end do
!
! finally, find if there is a way to speed up a prototypical calculation
! through dedicated loops
! the assumed architecture of the sequence file is:
! arbitrary residues / molecules first, terminated by water solvent of
! single type
  foundit = .false.
  use_waterloops = .false.
  rsw1 = nseq + 1
  k = 0
  do rs=1,nseq
    if (((seqtyp(rs).eq.39).OR.(seqtyp(rs).eq.40).OR.&
 &       (seqtyp(rs).eq.45).OR.(seqtyp(rs).eq.103)).AND.(foundit.EQV..false.)) then
      foundit = .true.
      k = seqtyp(rs)
      rsw1 = rs
      do rs2=rs+1,nseq
        if (seqtyp(rs2).ne.k) then
          rsw1 = nseq + 1
          foundit = .false.
          if ((seqtyp(rs2).ne.39).AND.(seqtyp(rs2).ne.40).AND.(seqtyp(rs2).ne.45).AND.&
 &            (seqtyp(rs2).ne.103)) then
            write(ilog,*) 'Warning. This calculation might be able to run &
 &faster if water molecules are placed at the end of sequence input.'
          else
            write(ilog,*) 'Warning. Concurrent use of different water models &
 &means that specialized water loops are no longer usable.'
          end if
          write(ilog,*)
          exit
        end if
      end do
      if (rsw1.eq.nseq) then
        foundit = .false.
      end if
      if (foundit.EQV..false.) then
        use_waterloops = .false.
        exit
      end if
    end if
    if (foundit.EQV..true.) then
      use_waterloops = .true.
      exit
    end if
  end do
! check for parameter conformity and homogeneity
  if (use_waterloops.EQV..true.) then
    dropout = .false.
    allocate(epsik(n-at(rsw1)%bb(1)+1))
    ii = b_type(at(rsw1)%bb(1))
    do rs=rsw1,nseq
      if ((attyp(at(rs)%bb(1)).ne.bio_ljtyp(b_type(at(rs)%bb(1)))).OR.(b_type(at(rs)%bb(1)).ne.ii)) then
        write(ilog,*) 'Warning. Applied biotype or LJ patches imply that for chosen rigid water model optimized &
 &loops are not available.'
        use_waterloops = .false.
        exit
      end if
      j = n-at(rs)%bb(1)+1
      if ((seqtyp(rs).eq.39).OR.(seqtyp(rs).eq.40)) then  
        do k=at(rs)%bb(2),at(rs)%bb(3)
          epsik(:) = lj_eps(attyp(k),attyp(at(rsw1)%bb(1):n))
          if (maxval(epsik).gt.0.0) then
            write(ilog,*) 'Warning. LJ parameters imply that for chosen 3-site water model optimized &
 &loops are not available.'
            dropout = .true.
            use_waterloops = .false.
            exit
          end if
        end do
        if (dropout.EQV..true.) exit
      else if ((seqtyp(rs).eq.45).OR.(seqtyp(rs).eq.103)) then
        do k=at(rs)%bb(2),at(rs)%bb(4)
          epsik(:) = lj_eps(attyp(k),attyp(at(rsw1)%bb(1):n))
          if (maxval(epsik).gt.0.0) then
            write(ilog,*) 'Warning. LJ parameters imply that for chosen 4-site water model optimized &
 &loops are not available.'
            dropout = .true.
            use_waterloops = .false.
            exit
          end if
        end do
        if (dropout.EQV..true.) exit
      else
        write(ilog,*) 'Fatal. Error in setting up specialized water loops. This &
 &is most likely an omission bug.'
        call fexit()
      end if
    end do
    deallocate(epsik)
  end if
!
! provide summary
!
 23   format('Molecule type ',i3,' has ',i6,' residues and a mass of ',g16.8,' D. ')
 24   format('   The first example molecule is # ',i6,' (residues ',i6,'-',&
 &i6,').')
 25   format('   There are ',i6,' molecules of this type.')
 26   format('Analysis group # ',i6,' is equivalent to molecule type # ',i6,'.')
 27   format('   These molecule(s) are (is) considered solvent.')
 28   format('   These molecule(s) are (is) considered solute(s).')
 29   format('Analysis group # ',i6,' is a subset of molecule type # ',i6,' and constitutes:')
 30   format('   Molecule # ',i6)
!
  if (seq_report.EQV..true.) then
    write(ilog,*)
    write(ilog,*) '--- Overview of sequence ---'
    do k=1,nmoltyp
      write(ilog,*)
      write(ilog,23) k,rsmol(moltyp(k,1),2)-rsmol(moltyp(k,1),1)+1,sum(mass(atmol(moltyp(k,1),1):atmol(moltyp(k,1),2)))
      write(ilog,24) moltyp(k,1),rsmol(moltyp(k,1),1),&
 &                                 rsmol(moltyp(k,1),2)
      write(ilog,25) moltyp(k,2)
      if (do_pol(k).EQV..false.) then
        write(ilog,*) '  Polymer analysis is turned off for this molec&
 &ule type.'
      else
        write(ilog,*) '  Polymer analysis is turned on for this molecu&
 &le type.'
      end if
      if (do_pers(k).EQV..false.) then
        write(ilog,*) '  Turns/angular statistics are turned off for t&
 &his molecule type.'
      else
        write(ilog,*) '  Turns/angular statistics are turned on for th&
 &is molecule type.'
      end if
      if (do_tors(k).EQV..false.) then
        write(ilog,*) '  Torsional analysis is turned off for this mol&
 &ecule type.'
      else
        write(ilog,*) '  Torsional analysis is turned on for this mole&
 &cule type.'
      end if
    end do
    write(ilog,*)
    write(ilog,*)
    write(ilog,*) '--- Overview of analysis groups ---'
    write(ilog,*)
    do k=1,nangrps
      if (moltyp(moltypid(molangr(k,1)),2).eq.molangr(k,2)) then
        write(ilog,26) k,moltypid(molangr(k,1))
        if (is_solvent(molangr(k,1)).EQV..true.) then
          write(ilog,27)
        else
          write(ilog,28)
        end if
      else
        write(ilog,29) k,moltypid(molangr(k,1))
        do j=1,nmol
          if (an_grp_mol(j).eq.k) then
            write(ilog,30) j
          end if
        end do
        if (is_solvent(molangr(k,1)).EQV..true.) then
          write(ilog,27)
        else
          write(ilog,28)
        end if
      end if
    end do
    write(ilog,*)
  end if
!
end
!
!-------------------------------------------------------------------------------
!
! initial structure randomization (may be silent --> globrandomize)
!
subroutine randomize_bb()
!
  use iounit
  use sequen
  use energies
  use atoms
  use molecule
  use cutoffs
  use aminos
  use ujglobals
  use polypep
  use pdb
!
  implicit none

  RTYPE et(MAXENERGYTERMS),scale_1,cut,dis,dis2,random
  RTYPE dum3(3),ccmat(6),oldvals(12),solmat(MAXSOLU*4,6),intsvec(9)
  integer i,j,cnt,azero,imol,jmol,aone,nlks,solucount,lastrs,firstrs,rs,rsi,k,kk,rs1,rs2,frs
  integer i1,i2,i3,i4,i5,i6,i7,i8,i9
  integer, ALLOCATABLE:: lklst(:,:)
  logical log_1,log_2,log_3,log_4,atrue,log_5,log_6,log_7,log_hs,log_8,notdone,dorand
  logical, ALLOCATABLE:: donewith(:)
!
  if (n_pdbunk.gt.0) call zmatfyc2() ! pointer arrays may not be up2date
!
  azero = 0
  aone = 1
  atrue = .true.
  scale_1 = scale_IPP
  log_hs = use_hardsphere
  scale_IPP = 1.0 ! note that this is only relevant if use_IPP is true, i.e., it's completely random if FMCSC_SC_IPP = 0.0
  use_hardsphere = .false.
  cut = mcnb_cutoff 
  mcnb_cutoff = 5.0
  log_1 = use_attLJ
  log_2 = use_CORR
  log_3 = use_IMPSOLV
  log_4 = use_WCA
  log_5 = use_POLAR
  log_6 = use_ZSEC
  log_7 = use_TOR
  log_8 = use_FEG
  use_attLJ = .false.
  use_CORR = .false.  
  use_IMPSOLV = .false.
  use_WCA = .false.
  use_POLAR = .false.
  use_ZSEC = .false.
  use_TOR = .false.
  use_FEG = .false.
  et(:) = 0.0
  if (pdb_ihlp.le.0) then
    frs = nseq
  else
    frs = pdb_ihlp
  end if
!
  do imol=1,nmol
!   for no internal randomization exit right away
    if ((globrandomize.eq.0).OR.(globrandomize.eq.2)) exit
    if ((globrandomize.eq.4).AND.(rsmol(imol,2).le.frs)) cycle
!
    allocate(lklst(rsmol(imol,2)-rsmol(imol,1)+1,2))
    nlks = 0
    firstrs = rsmol(imol,2)
    lastrs = rsmol(imol,1)
    do i=1,n_crosslinks
      if ((molofrs(crosslink(i)%rsnrs(1)).eq.imol).AND.&
 &        (molofrs(crosslink(i)%rsnrs(2)).eq.imol)) then
        nlks = nlks + 1
        lklst(nlks,1) = crosslink(i)%rsnrs(1)
        lklst(nlks,2) = crosslink(i)%rsnrs(2)
        if ((minval(lklst(nlks,1:2)).le.frs).AND.(maxval(lklst(nlks,1:2)).gt.frs)) then
          write(ilog,*) 'Warning. Intramolecular crosslink is ignored during partial &
 &structure randomization (residues ',lklst(nlks,1),' and ',lklst(nlks,2),').'
          nlks = nlks - 1
          cycle
        end if
        if (lklst(nlks,1).lt.firstrs) firstrs = lklst(nlks,1)
        if (lklst(nlks,1).gt.lastrs) lastrs = lklst(nlks,1)
        if (lklst(nlks,2).lt.firstrs) firstrs = lklst(nlks,2)
        if (lklst(nlks,2).gt.lastrs) lastrs = lklst(nlks,2)
        do j=1,nlks-1
          if (lklst(nlks,1).lt.lklst(j,1)) then
            write(ilog,*) 'Fatal. List of crosslinks is not in expected order. This&
 & is a bug.'
            call fexit()
          end if
        end do
      end if
    end do
    if ((rsmol(imol,1).le.frs).AND.(rsmol(imol,2).gt.frs)) then
      rs1 = frs + 1
    else
      rs1 = rsmol(imol,1)
    end if
!   building randomized conformations gets extremely messy for crosslinked peptides
!   due to the implied loop closures; without them, it's trivial
    if (nlks.eq.0) then
      call stretch_randomize(imol,rs1,rsmol(imol,2))
    else
      call stretch_randomize(imol,rs1,firstrs-1)
      k = 0
      allocate(donewith(nlks))
      donewith(:) = .false.
      do while (k.lt.nlks)
!       pick a crosslink whose 2nd residue is within another crosslink stretch
        notdone = .true.
        do i=1,nlks
          if (donewith(i).EQV..true.) cycle
          do j=1,nlks
            if ((lklst(i,2).lt.lklst(j,2)).AND.(lklst(i,2).gt.lklst(j,1))) then
              notdone = .false.
              if (donewith(j).EQV..true.) then
                write(ilog,*) 'Fatal. Crosslink topology is too complex for CAMPARI&
 & to close all loops consistently. Please check back later.'
                call fexit()
              end if
              exit
            end if
          end do
          if (notdone.EQV..false.) exit
        end do
        if (notdone.EQV..true.) then
          do i=1,nlks
            if (donewith(i).EQV..false.) exit
          end do
        end if
!
!       now sample whatever we can randomly, then close crosslink via concerted rotation
        notdone = .true.
        kk = 0
        do while (notdone.EQV..true.)
          kk = kk + 1
          do rs=lklst(i,1),lklst(i,2)
            dorand = .true.
            do j=1,nlks
              if (donewith(j).EQV..false.) cycle
              if ((rs.ge.lklst(j,1)).AND.(rs.le.lklst(j,2))) then
                dorand = .false.
                exit
              end if
            end do
            if (dorand.EQV..false.) cycle
            call stretch_randomize(imol,rs,rs)
          end do
          call crosslink_indices(i,i1,i2,i3,i4,i5,i6,i7,i8,i9)
          call crosslink_refvals(i,i1,i2,i3,i4,i5,i6,i7,i8,i9,ccmat,oldvals)
          call torchainclose(i1,i2,i3,i7,i8,i9,ccmat,solmat,oldvals,solucount)
          if (solucount.ge.1) then
            notdone = .false.
            call crosslink_picksolu(i,i1,i2,i3,i4,i5,i6,i7,i8,i9,solucount,solmat,ccmat,oldvals)
          end if
          if (kk.gt.max(globrdmatts,10000)) then
            write(ilog,*) 'Fatal. Randomly determined crosslinks led to&
 & non-closable topology within ',max(globrdmatts,10000),' attempts. Simplify crosslinks or try again.'
            write(ilog,*) 'This was the crosslink between residues ',lklst(i,1),' and ',&
 &lklst(i,2),'.'
            call fexit()
          end if
        end do
        donewith(i) = .true.
        k = k + 1
      end do
      deallocate(donewith)
      rs = firstrs+1
      do while (rs.lt.lastrs-1)
        dorand = .true.
        do i=1,nlks
          if ((rs.ge.lklst(i,1)).AND.(rs.le.lklst(i,2))) then
            dorand = .false.
            exit
          end if
        end do
        if (dorand.EQV..true.) then
          rsi = rs
          notdone = .true.
          do while (notdone.EQV..true.)
            rs = rs + 1
            if (rs.eq.lastrs-1) notdone = .false.
            do i=1,nlks
              if ((rs.ge.lklst(i,1)).AND.(rs.le.lklst(i,2))) then
                notdone = .false.
                rs = rs - 1
              end if
            end do
          end do
          dorand = .false.
          call stretch_randomize(imol,rsi,rs)
          rs = rs + 1
        else
          rs = rs + 1
        end if
      end do
      call stretch_randomize(imol,lastrs+1,rsmol(imol,2))
    end if
    deallocate(lklst)
  end do
!
  do imol=1,nmol
    call update_rigid(imol)
  end do
!
! now change the rigid-body degrees of freedom to eliminate intermolecular clashes
  do imol=1,nmol
!
    if (((globrandomize.eq.0).OR.(globrandomize.eq.3)).AND.(imol.le.molofrs(frs))) cycle ! in general, frs < nmol
    cnt = 0
    et(1) = 0.25*HUGE(et(1))
    et(12) = 0.25*HUGE(et(12))
!
    do while ((et(1)+et(12)).gt.globrdmthresh)
      cnt = cnt + 1
!     obviously we got into a mess: clear-up
      if (cnt.gt.globrdmatts) then
        write(ilog,*) 'WARNING: Could not arrange molecules &
 &with non-dramatic overlap (for molecule ',imol,').'
        exit
      end if
!     randomize pos for molecule imol
      call randomize_trans(imol,dum3)
      do i=atmol(imol,1),atmol(imol,2)
        x(i) = x(i) + dum3(1)
        y(i) = y(i) + dum3(2)
        z(i) = z(i) + dum3(3)
      end do
      call randomize_rot(imol)
!     make sure intermolecular energies (only backward) are sane
      et(1) = 0.0
      et(12) = 0.0
      do i=rsmol(imol,1),rsmol(imol,2)
        call e_boundary_rs(i,et,azero)
        do j=1,rsmol(imol,1)-1
          if (j.lt.rsmol(imol,1)) then
            call dis_bound(refat(i),refat(j),dis2) 
            dis = sqrt(dis2)
            if (dis.lt.(mcnb_cutoff+resrad(i)+resrad(j))) then
              call Ven_rsp(et,j,i,atrue)
            end if
          end if
        end do
      end do
    end do
    call update_rigid(imol)
  end do
!
  do i=1,n_crosslinks
    if ((globrandomize.eq.0).OR.(globrandomize.eq.3)) exit
    rs1 = crosslink(i)%rsnrs(1)
    imol = molofrs(rs1)
    rs2 = crosslink(i)%rsnrs(2)
    jmol = molofrs(rs2)
    if (imol.eq.jmol) cycle
    if (crosslink(i)%itstype.eq.2) then
      intsvec(3) = 2.03 ! h.c. -S-S
      intsvec(6) = 1.822 ! h.c. - match with sidechain
      intsvec(9) = 1.53  ! h.c. - match with sidechain
      intsvec(2) = 103.0 ! h.c. -CB-S-S
      intsvec(5) = 103.0 ! h.c. -S-S-CB
      intsvec(8) = 114.4 ! h.c. - match with sidechain
    else
      write(ilog,*) 'Fatal. Encountered unsupported crosslink type in randomize_bb(...).&
 & This is most likely an omission bug.'
      call fexit()
    end if
    et(1) = 0.25*HUGE(et(1))
    et(12) = 0.25*HUGE(et(12))
    cnt = 0
!
    do while ((et(1)+et(12)).gt.globrdmthresh)
      cnt = cnt + 1
      if (cnt.gt.globrdmatts) then
        write(ilog,*) 'WARNING: Could not arrange crosslinked chains &
 &with non-dramatic overlap (molecule ',imol,').'
        exit
      end if
      if (crosslink(i)%itstype.eq.2) then
        intsvec(1) = random()*360.0 - 180.0
        intsvec(4) = random()*360.0 - 180.0
        intsvec(7) = random()*360.0 - 180.0
      else
        write(ilog,*) 'Fatal. Encountered unsupported crosslink type in randomize_bb(...).&
 & This is most likely an omission bug.'
        call fexit()
      end if
      call crosslink_follow(i,imol,jmol,intsvec)
!     make sure intermolecular energies (only backward) are sane
      et(1) = 0.0
      et(12) = 0.0
      do k=rsmol(imol,1),rsmol(imol,2)
        call e_boundary_rs(k,et,azero)
        do j=1,rsmol(imol,1)-1
          if (j.lt.rsmol(imol,1)) then
            call dis_bound(refat(k),refat(j),dis2) 
            dis = sqrt(dis2)
            if (dis.lt.(mcnb_cutoff+resrad(k)+resrad(j))) then
              call Ven_rsp(et,j,k,atrue)
            end if
          end if
        end do
      end do
    end do
  end do
!
! finalize a bit
  use_attLJ = log_1
  use_CORR = log_2
  use_IMPSOLV = log_3
  use_WCA = log_4
  use_POLAR = log_5
  use_ZSEC = log_6
  use_TOR = log_7
  use_FEG = log_8
  scale_IPP = scale_1
  use_hardsphere = log_hs
  mcnb_cutoff = cut
!
end
!
!-----------------------------------------------------------------
!
subroutine stretch_randomize(imol,rsi,rsf)
!
  use aminos
  use iounit
  use energies
  use sequen
  use molecule
  use polypep
  use fyoc
  use cutoffs
  use zmatrix
!
  implicit none
!
  RTYPE et(MAXENERGYTERMS),dis,dis2,etbu
  RTYPE, ALLOCATABLE:: curvals(:,:),potstr(:,:)
  integer imol,rsi,rsf,i,j,npots,cnt,azero,aone,atwo,athree,dofcnt,dofcntbu,dofcntbu2,alcsz,phs1,phs2,k
  logical afalse,atrue
!
  if (rsi.gt.rsf) return
!
  alcsz = at(rsf)%bb(1)-at(rsi)%bb(1)+at(rsf)%nbb+at(rsf)%nsc
  allocate(curvals(alcsz,2))
  afalse = .false.
  atrue = .true.
  azero = 0
  aone = 1
  atwo = 2
  athree = 3
  dofcnt = 0
  dofcntbu2 = 0
  phs1 = 1*globrdmatts/3 ! at least 1
  phs2 = 2*globrdmatts/3 ! at least 2
!
! setup up an extra potential if this is a randomization within an intramolecular crosslink loop
  if (n_crosslinks.gt.0) allocate(potstr(n_crosslinks,4))
  npots = 0
  do i=1,n_crosslinks
    if ((molofrs(crosslink(i)%rsnrs(1)).eq.imol).AND.&
 &      (molofrs(crosslink(i)%rsnrs(2)).eq.imol)) then
      if ((rsi.ge.minval(crosslink(i)%rsnrs(1:2))).AND.(rsf.le.maxval(crosslink(i)%rsnrs(1:2))).AND.&
 &        (crosslink(i)%itstype.eq.1)) then
        npots = npots + 1
        potstr(npots,1) = 1.0*cai(crosslink(i)%rsnrs(1))
        potstr(npots,2) = 1.0*cai(crosslink(i)%rsnrs(2))
        potstr(npots,3) = 3.0 + 2.0*(maxval(crosslink(i)%rsnrs(1:2))-rsf) ! max dis before penalty sets in
        potstr(npots,4) = 0.1*globrdmthresh ! kcal/molA^2
      end if
    end if
  end do
! 
  do i=rsi,rsf
    et(1) = HUGE(etbu)
    etbu = et(1)
    cnt = 0
    if (i.gt.rsi) dofcntbu2 = dofcntbu
    dofcntbu = dofcnt
    do while (etbu.gt.globrdmthresh)
      cnt = cnt + 1
!     obviously we got into a mess: clear-up
      if (cnt.gt.globrdmatts) exit
!
      dofcnt = dofcntbu
      if ((cnt.gt.phs2).AND.(i.gt.rsi)) then
!       randomize bb-angles for residue i-1 (this becomes expensive quickly)
        dofcnt = dofcntbu2
        call setrandom_forres(i-1,alcsz,curvals,dofcnt,athree,atrue)
        if (dofcnt.ne.dofcntbu) then
          write(ilog,*) 'Fatal. D.o.f counter is corrupt in stretch_randomize(...). This is a bug.'
          call fexit()
        end if
      else if (cnt.gt.phs1) then
!       add side chains if still failing
        call setrandom_forres(i,alcsz,curvals,dofcnt,atwo,atrue)
        dofcnt = dofcntbu
      end if
!     randomize bb-angles for residue i
      call setrandom_forres(i,alcsz,curvals,dofcnt,aone,atrue)
!     make sure backward-in-chain-energies are sane
      et(1) = 0.0
      do j=rsmol(imol,1),min(i+1,rsf)
        call dis_bound(refat(i),refat(j),dis2)
        dis = sqrt(dis2)
        if (dis.lt.(mcnb_cutoff+resrad(i)+resrad(j))) then
          call Ven_rsp(et,min(j,i),max(j,i),atrue)
        end if
        if ((cnt.gt.phs2).AND.(j.ne.i).AND.(i.gt.rsi)) then
          call dis_bound(refat(i-1),refat(j),dis2) 
          dis = sqrt(dis2)
          if (dis.lt.(mcnb_cutoff+resrad(i-1)+resrad(j))) then
            call Ven_rsp(et,min(j,i-1),max(j,i-1),atrue)
          end if
        end if
      end do
      
      
      do j=1,npots
        call dis_bound(nint(potstr(j,1)),nint(potstr(j,2)),dis2)
        dis = sqrt(dis2)
        if (dis.gt.potstr(j,3)) then
          et(1) = et(1) + potstr(j,4)*(-potstr(j,3))*(dis-potstr(j,3))
        end if
      end do

      if ((et(1).lt.etbu).AND.((dofcnt.gt.dofcntbu).OR.((cnt.gt.phs2).AND.(i.gt.rsi).AND.(dofcnt.gt.dofcntbu2)))) then
        if ((cnt.gt.phs2).AND.(i.gt.rsi)) then
          curvals((dofcntbu2+1):dofcnt,2) = curvals((dofcntbu2+1):dofcnt,1)
        else
          curvals((dofcntbu+1):dofcnt,2) = curvals((dofcntbu+1):dofcnt,1)
        end if
        etbu = et(1)
      end if
      
      if ((dofcnt.eq.dofcntbu).AND.(cnt.lt.phs1)) then
!       turn on side chain addition right away
        cnt = phs1
      else if ((dofcnt.eq.dofcntbu).AND.(cnt.lt.phs2)) then
!       turn on panic mode right away
        cnt = phs2
      end if
      if (cnt.eq.globrdmatts) then ! still here means trouble - restore best guess
        if (i.gt.rsi) then
          curvals((dofcntbu2+1):dofcnt,1) = curvals((dofcntbu2+1):dofcnt,2)
          dofcnt = dofcntbu2
          call setrandom_forres(i-1,alcsz,curvals,dofcnt,athree,afalse)
          call setrandom_forres(i,alcsz,curvals,dofcnt,athree,afalse)
        else
          curvals((dofcntbu+1):dofcnt,1) = curvals((dofcntbu+1):dofcnt,2)
          dofcnt = dofcntbu
          call setrandom_forres(i,alcsz,curvals,dofcnt,athree,afalse)
        end if
      end if
    end do
  end do
  deallocate(curvals)
  if (n_crosslinks.gt.0) deallocate(potstr)
!
end
!
!----------------------------------------------------------------------------------
!
! a small helper for stretch_randomize
!
subroutine setrandom_forres(i,alcsz,curvals,curcur,mode,dosam)
!
  use fyoc
  use polypep
  use molecule
  use atoms
  use zmatrix
  use sequen
!
  implicit none
!
  integer alcsz,curcur,mode,j,ati,i,azero
  RTYPE curvals(alcsz,2),random
  logical isnat,dosam
!
  azero = 0
!
  if ((mode.eq.1).OR.(mode.eq.3)) then
    if (fline(i).gt.0) then
      if (izrot(fline(i))%alsz.gt.0) then
        if (seqflag(i).ne.5) then
          if (dosam.EQV..true.) curvals(curcur+1,1) = random()*360.0 - 180.0
        else
          if (dosam.EQV..true.) curvals(curcur+1,1) = phi(i)
        end if
      else
        if (dosam.EQV..true.) curvals(curcur+1,1) = phi(i)
      end if
    end if
    if (yline(i).gt.0) then
      if (izrot(yline(i))%alsz.gt.0) then
        if (dosam.EQV..true.) curvals(curcur+2,1) = random()*360.0 - 180.0
      else
        if (dosam.EQV..true.) curvals(curcur+2,1) = psi(i)
      end if
    end if
  end if
  if ((fline(i).gt.0).OR.(yline(i).gt.0)) then
    if ((mode.eq.1).OR.(mode.eq.3)) call setfy(i,curvals(curcur+1,1),curvals(curcur+2,1),azero)
    curcur = curcur + 2
  end if
  do ati=at(i)%bb(1),at(i)%bb(1)+at(i)%nbb+at(i)%nsc-1
    if (izrot(ati)%alsz.gt.0) then
      if ((ati.eq.fline(i)).OR.(ati.eq.yline(i)).OR.(ati.eq.wline(i))) cycle
      isnat = .false.
      if ((ati.eq.nucsline(1,i)).OR.(ati.eq.nucsline(2,i)).OR.(ati.eq.nucsline(3,i)).OR.&
 &        (ati.eq.nucsline(4,i))) isnat = .true.
      if (i.gt.1) then
        if ((ati.eq.nucsline(5,i-1)).OR.(ati.eq.nucsline(3,i-1))) isnat = .true.
      end if
      if ((izrot(ati)%rotis(izrot(ati)%alsz,2).gt.(at(i)%bb(1)+at(i)%nbb+at(i)%nsc-1)).OR.(isnat.EQV..true.)) then
        curcur = curcur + 1
        if ((mode.eq.1).OR.(mode.eq.3)) then
          if (dosam.EQV..true.) curvals(curcur,1) = random()*360.0 - 180.0
          call setother(ati,curvals(curcur,1),azero,isnat)
        end if
      end if
    end if
  end do
  do ati=at(i)%bb(1),at(i)%bb(1)+at(i)%nbb+at(i)%nsc-1
    if (izrot(ati)%alsz.gt.0) then
      if ((ati.eq.fline(i)).OR.(ati.eq.yline(i)).OR.(ati.eq.wline(i))) cycle
      if ((ati.eq.nucsline(1,i)).OR.(ati.eq.nucsline(2,i)).OR.(ati.eq.nucsline(3,i)).OR.&
 &        (ati.eq.nucsline(4,i))) cycle
      if (i.gt.1) then
        if ((ati.eq.nucsline(5,i-1)).OR.(ati.eq.nucsline(3,i-1))) cycle
      end if
      isnat = .false.
      do j=1,nchi(i)
        if (ati.eq.chiline(j,i)) isnat = .true.
      end do
      if ((seqtyp(i).ne.26).AND.(isnat.EQV..false.)) cycle  ! don't perturb those in unslst
      if (atmres(izrot(ati)%rotis(izrot(ati)%alsz,2)).eq.atmres(ati)) then
        curcur = curcur + 1
        if ((mode.eq.2).OR.(mode.eq.3)) then
          if (dosam.EQV..true.) curvals(curcur,1) = random()*360.0 - 180.0
          call setother(ati,curvals(curcur,1),azero,isnat)
        end if
      end if
    end if
  end do
  call makexyz_forbb(i)
!
end
!
!--------------------------------------------------------------------------------------------------------
!
