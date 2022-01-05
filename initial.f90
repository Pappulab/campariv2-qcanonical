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
! CONTRIBUTIONS: Rohit Pappu, Adam Steffen, Albert Mao                     !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!-----------------------------------------------------------------------
!
subroutine initialize_charge_state_cste_net_Q() ! The advantage of the constant net charge is that the number of free proton is kept 
    !constant
    ! Of course, the FOS will vary a lot.
    use atoms
    use energies 
    use movesets
    use sequen
    use aminos ! Martin : debug only 
    use zmatrix
    use accept 
    
    
    implicit none
    integer pos_ion(nseq),pi_c
    integer neg_ion(nseq),ni_c
    integer pos_res(nseq),pr_c
    integer neg_res(nseq),nr_c
    integer i,rs,net_q,tot_up,j
    integer salt_n,salt_p,res_n,res_p

    RTYPE temp,random,temp_proba,temp_proba_pos,temp_proba_neg
    
    allocate(seq_q_state(nseq))
    allocate(seq_q_state_save(nseq))

    allocate(seq_q_H_dir(nseq))
    allocate(q_switch(nseq))
    allocate(conjugate(nseq))
    allocate(par_hsq(nseq,4))
    allocate(par_hsq_save(nseq,4))
    allocate(hsq_sign(nseq))
    allocate(blen_save(size(blen)))
    allocate(bang_save(size(bang)))
    
    allocate(ztor_save(size(ztor)))

    ! Ok so I am not going to worry about this now, but I will not test for overall neutral system.
    ! This should be done before anyway
    ! Note that in htis routine, I start from the charged states, unlike in the constnat FOS one
    
    print *,"Executing HSQ in the constant charge mode"
    
    pos_ion(:)=0
    neg_ion(:)=0
    pos_res(:)=0
    neg_res(:)=0
    
    salt_n=0 
    salt_p=0 
    res_n=0 
    res_p=0 
    
    pi_c=0 
    ni_c=0 
    pr_c=0 
    nr_c=0 
    
    tot_up=0
    
    selec_proba=0
    
    if (MAX_SWITCH.eq.0) then 
        MAX_SWITCH=1
    end if 

    ! the problem is that I need a good way of finding atoms and conjugate 
    seq_q_state(:)=999
    seq_q_H_dir(:)=0
    pos_res(:)=0 ! Note that the neg ion work differently from the neg res array
    
    salt_n=0
    
    do rs=1,nseq
        if ((amino(seqtyp(rs)).eq.'NAX').or.(amino(seqtyp(rs)).eq.'CLX')) then 
            print *,"You should not have titrable ions in a constant charge HSQ, aborting"
            call fexit()

        else if ((amino(seqtyp(rs)).eq.'GLX').or.(amino(seqtyp(rs)).eq.'ASX').or.(amino(seqtyp(rs)).eq.'TYX')&
            .or.(amino(seqtyp(rs)).eq.'TXP').or.(amino(seqtyp(rs)).eq.'SXP').or.(amino(seqtyp(rs)).eq.'YXP')) then 
            nr_c=nr_c+1
            neg_res(rs)=1
            seq_q_state(rs)=1
            tot_up=tot_up+1
            seq_q_H_dir(rs)=2 ! That means the atom appears 
            hsq_sign(rs)=-1
            
        else if ((amino(seqtyp(rs)).eq.'LYX').or.(amino(seqtyp(rs)).eq.'HDX').or.(amino(seqtyp(rs)).eq.'HEX')&
            &.or.(amino(seqtyp(rs)).eq.'HIX')) then
            pr_c=pr_c+1
            pos_res(rs)=1
            seq_q_state(rs)=1
            hsq_sign(rs)=+1
            tot_up=tot_up+1
            seq_q_H_dir(rs)=1 ! That means the atom disapears 
            
        else if ((amino(seqtyp(rs)).eq.'NA+')) then 
            salt_p=salt_p +1
            seq_q_state(rs)=0
            hsq_sign(rs)=0 ! That is as far as the titration is concerned
            
        else if ((amino(seqtyp(rs)).eq.'CL-'))  then 
            salt_n=salt_n +1
            seq_q_state(rs)=0
            hsq_sign(rs)=0
            
        else if ((amino(seqtyp(rs)).eq.'LYS').or.(amino(seqtyp(rs)).eq.'HIP').or.(amino(seqtyp(rs)).eq.'ARG')) then 
            res_p=res_p+1
            seq_q_state(rs)=0
            hsq_sign(rs)=0 
            
        else if ((amino(seqtyp(rs)).eq.'ASP').or.(amino(seqtyp(rs)).eq.'GLU')) then 
            res_n=res_n+1
            seq_q_state(rs)=0
            hsq_sign(rs)=0
            
        else if ((amino(seqtyp(rs)).eq.'T1P').or.(amino(seqtyp(rs)).eq.'S1P').or.(amino(seqtyp(rs)).eq.'Y1P')) then 
            res_n=res_n+1
            seq_q_state(rs)=0
            hsq_sign(rs)=0
            
        else if ((amino(seqtyp(rs)).eq.'T2P').or.(amino(seqtyp(rs)).eq.'S2P').or.(amino(seqtyp(rs)).eq.'Y2P')) then 
            res_n=res_n+1
            seq_q_state(rs)=0
            hsq_sign(rs)=0

        end if 
    end do

    ! The total charge of the protein is supposed to be tot_QX (I know the name is misleading)
    ! 
    
    if ((salt_p-salt_n).ne.(-tot_QX)) then 
        print *,salt_p,salt_n,tot_QX
        print *,"The system is not net neutral, check your sequence. You have an ion net charge of ",salt_p-salt_n,&
        &" a protein net charge target of ",(tot_QX)," with ",res_p,"non ionazible positive residues, and "&
        &,res_n,"non ionazible negative residues "
        call fexit()
    end if 
    
    if (read_charge_state.eqv..true.) then 
        n_res=nseq-(salt_n+salt_p)
        call read_charge_states()
        if (n_charge_states.eq.1) then 
            hsq_freq=0.0
        end if 
    end if
    
    if (n_charge_states.eq.1) then !This is to make sure no infinite loop will happen
        hsq_freq=0.0
    end if 
    
    n_titrable=pr_c+nr_c
    
    allocate(titrable_res_index(n_titrable))

    ! Martin : Truned that off to avoid preferential protonation that is due to the more stable ion.
    ! Does not matter given that the walk is in a net charge dimension.
    if ((nr_c.eq.0).and.(pr_c.eq.0)) then 
        print *,"No point in using HSQ with no titrable residues, please turn off"
        call fexit()
    end if 
    if (nr_c.ne.0.0) then 
        selec_proba_neg=(1.0)/nr_c
    else 
        selec_proba_neg=0.0
    end if 
    if (pr_c.ne.0.0) then 
        selec_proba_pos=(1.0)/pr_c
    else 
        selec_proba_pos=0.0
    end if 
    
    selec_proba=1./(nr_c+pr_c)!selec_proba_neg+selec_proba_pos
    ! I think I need to add the 'inert' residues here. We will see, but I think it currently will generate ifinity loops
    ! Martin : Ok I am adding it, but not going to test it, so if not working, may be this
    net_q=-nr_c+pr_c+res_p-res_n
    if (nr_c.ne.0) then 
        temp_proba_neg=1./nr_c !tot_X*selec_proba_neg
    else 
        temp_proba_neg=0
    end if 
    if (pr_c.ne.0) then 
        temp_proba_pos=1./pr_c !tot_X*selec_proba_pos
    else 
        temp_proba_pos=0
    end if 

    pr_c=1
    nr_c=1 
    ! okay why is that one again ? : because it is an index of the actual array
    ! First : how do I decide which charge state to start with 
!    do while (tot_up.lt.tot_X) 
!    Having a constant net charge means that the state of the groups will have to be changed simultenaously

!    do while (net_q.ne.tot_QX) ! All this isn't very efficient, but it is only intitilization anyway
!        if (net_q.lt.tot_QX) then ! If the charge is less than the prescribed charge : uncharge a negative residue
!            do i=1,nseq
!                temp=random()
!                if ((temp.lt.temp_proba_neg).and.(net_q.lt.tot_QX).and.(seq_q_state(i).eq.1)) then 
!                    if (1.eq.neg_res(i)) then
!                        seq_q_state(i)=2
!                        neg_res(i)=0
!                        tot_up=tot_up-1
!                        nr_c=nr_c+1 ! That just means don't try it again
!                        net_q=net_q+1
!                    end if 
!                else if (net_q.ge.tot_QX) then   
!                    exit
!                end if 
!            end do 
!        else if (net_q.gt.tot_QX) then 
!            do i=1,nseq
!                temp=random()
!
!                
!!                print *,(temp.lt.temp_proba_pos),(net_q.gt.tot_QX),seq_q_state(i),i,nseq
!                if ((temp.lt.temp_proba_pos).and.(net_q.gt.tot_QX).and.(seq_q_state(i).eq.1)) then        
!                    if (1.eq.pos_res(i)) then
!                        seq_q_state(i)=2
!                        pos_res(i)=0
!                        tot_up=tot_up-1
!                        pr_c=pr_c+1
!                        net_q=net_q-1
!                    end if 
!                else if (net_q.le.tot_QX) then   
!                    exit
!                end if 
!            end do 
!        else if (net_q.eq.tot_QX) then   
!                exit
!        end if 
!    end do 

    j=1
    do i=1,nseq
        if ((seq_q_state(i).eq.1).or.(seq_q_state(i).eq.2).or.(seq_q_state(i).eq.3)) then
            titrable_res_index(j)=i
            j=j+1
        end if 
    end do 
    
    if (read_charge_state.eqv..true.) then 
        i=CEILING(random()*(n_charge_states))
        acc%curr_tries(:)=i
        seq_q_state(1:n_res)=q_states(i,:)
    else      ! dont think anything of this is useful 
        do while (net_q.ne.tot_QX) ! All this isn't very efficient, but it is only initialization anyway
            i=CEILING(random()*(n_res))
            if ((net_q.lt.tot_QX).and.(1.eq.neg_res(i)).and.(seq_q_state(i).eq.1)) then 
                ! If the charge is less than the prescribed charge : uncharge a negative residue
                seq_q_state(i)=2
                neg_res(i)=0
                tot_up=tot_up-1
                nr_c=nr_c+1 ! That just means don't try it again
                net_q=net_q+1
                
            else if ((net_q.gt.tot_QX).and.(1.eq.pos_res(i)).and.(seq_q_state(i).eq.1)) then 
                seq_q_state(i)=2
                pos_res(i)=0
                tot_up=tot_up-1
                pr_c=pr_c+1
                net_q=net_q-1
                
            else if (net_q.eq.tot_QX) then   
                exit
            end if 
        end do 
    end if 
!        print *,"bababab"
    do i=1,nseq
!        print *,rs
        if (seq_q_state(i).eq.1) then
            par_hsq(i,:)=0
        else if ((seq_q_state(i).eq.2).or.(seq_q_state(i).eq.3)) then 
            par_hsq(i,:)=1
        else if ((seq_q_state(i).eq.0).or.(seq_q_state(i).eq.-1)) then 
            par_hsq(i,:)=0
        end if 
!        print *,par_hsq(i,:)
!        print *,seq_q_state(i)
    end do 
    
    do i=1,nhis 
        if (seq_q_state(his_state(i,1)).eq.3) then 
            his_state(i,2)=1
        else 
            his_state(i,2)=0
        end if 
    end do 
    
    print *,"Starting with charge state :",seq_q_state(:)

end subroutine

subroutine initialize_charge_state_cste_FOS()
    use atoms
    use energies 
    use movesets
    use sequen
    use zmatrix
    use aminos ! Martin : debug only 
    
    
    implicit none 
    integer pos_ion(nseq),pi_c
    integer neg_ion(nseq),ni_c
    integer pos_res(nseq),pr_c
    integer neg_res(nseq),nr_c

    
    integer i,rs, tot_up
    
    RTYPE temp,random,temp_proba
    
    allocate(seq_q_state(nseq))
    allocate(seq_q_state_save(nseq))
    allocate(q_switch(nseq))
    allocate(conjugate(nseq))
    allocate(par_hsq(nseq,4))
    allocate(par_hsq_save(nseq,4))
    allocate(ztor_save(n+20))
    
    allocate(blen_save(size(blen)))
    allocate(bang_save(size(bang)))
    
    print *,"HSQ modes other than constant Q are not yet fully implemented, aborting"
    call fexit()
    
    print *,"Executing HSQ in the constant FOS mode"
    
     

    
    pos_ion(:)=0
    neg_ion(:)=0
    pos_res(:)=0
    neg_res(:)=0
        
    pi_c=0 
    ni_c=0 
    pr_c=0 
    nr_c=0 
    
    selec_proba=0

    if (MAX_SWITCH.eq.0) then 
        MAX_SWITCH=TOT_X
    end if 
    
!   the problem is that I need a good way of finding atoms and conjugate 
   
    seq_q_state(:)=999
    
   ! I also need a conjugate ion array 
    do rs=1,nseq
        if ((seqtyp(rs).eq.117).or.(seqtyp(rs).eq.118)) then
            seq_q_state(rs)=-1
            if (seqtyp(rs).eq.117) then 
                ni_c=ni_c+1
                neg_ion(ni_c)=rs
            else
                pi_c=pi_c+1
                pos_ion(pi_c)=rs
            end if 
            
        else if (seqtyp(rs).eq.116) then 
            nr_c=nr_c+1
            neg_res(nr_c)=rs
            seq_q_state(rs)=2
            !seq_q_H_dir(rs)=2 ! That means the atom appears 
        else if (seqtyp(rs).eq.119) then 
            pr_c=pr_c+1
            pos_res(pr_c)=rs
            seq_q_state(rs)=2
            !seq_q_H_dir(rs)=1 ! That means the atom disapears 
        else
            seq_q_state(rs)=0
            !seq_q_H_dir(rs)=0
        end if 
    end do

    if ((nr_c.ne.pi_c).or.(pr_c.ne.ni_c)) then 
        print *,"Imbalance in the titrable charge groups for charge change hamiltonian switch. Aborting."
        call fexit()
    end if 
    
    ! That is when I iterate through the molcecule, what is the probability of me getting sellecting the current residue

    selec_proba=(1.0)/(nr_c+pr_c)

    temp_proba=tot_X*selec_proba

    pr_c=1
    nr_c=1

    do i=1,nseq
        if (i.eq.pos_res(pr_c)) then 
            conjugate(i)=neg_ion(pr_c)
            pr_c=pr_c+1
        else if (i.eq.neg_res(nr_c)) then
            conjugate(i)=pos_ion(nr_c)
            nr_c=nr_c+1
        else 
            conjugate(i)=0
        end if 
            
    end do
    
    tot_up=0
    
    pr_c=1
    nr_c=1
    
    do while (tot_up.lt.tot_X) 
        do i=1,nseq
            temp=random()
            if ((temp.lt.temp_proba).and.(tot_up.lt.tot_X).and.(seq_q_state(i).eq.2)) then 
                if (i.eq.pos_res(pr_c)) then
                    
                    seq_q_state(i)=1
                    tot_up=tot_up+1
                    pr_c=pr_c+1
                    
                else if (i.eq.neg_res(nr_c)) then
                    
                    seq_q_state(i)=1
                    tot_up=tot_up+1
                    nr_c=nr_c+1

                end if 
                
            else if (tot_up.eq.tot_X) then   
                exit
            end if 
        end do 
    end do 
    
    do i=1,nseq
        if (seq_q_state(i).eq.1) then
            par_hsq(i,:)=0
            par_hsq(conjugate(i),:)=0
        else if (seq_q_state(i).eq.2) then 
            par_hsq(i,:)=1
            par_hsq(conjugate(i),:)=1
        else if (seq_q_state(i).eq.0) then 
            par_hsq(i,:)=0
        end if 
    end do 
end subroutine

! large-scale initialization to be able to support "lazy"
! key-files
subroutine initial()
!
  use commline
  use atoms
  use distrest
  use inter
  use iounit
  use keys
  use mcsums
  use math
  use pdb
  use polyavg
  use sequen
  use grids
  use mcgrid
  use units
  use polypep
  use fyoc
  use molecule
  use params
  use torsn
  use ionize
  use zmatrix
  use movesets
  use paircorr
  use mpistuff
  use cutoffs
  use accept
  use contacts
  use energies
  use tabpot
  use system
  use fos
  use dipolavg
  use diffrac
  use forces
  use ewalds
  use mini
  use grandensembles
  use ujglobals
  use dssps
  use threads
  use shakeetal
  use clusters
  use ems
  use wl
!
  implicit none
!
  integer i,j,t1,t2
!
! THIS CONSTRUCT IS OBSOLETE: NON-STANDARD FORTRAN DUE TO CONVERSION LOGICAL->INT
! IT WILL WORK WITH MOST MODERN COMPILERS, HOWEVER
!!     see how the compiler deals with logical-to-integer conversions and set
!!     appropriate parameters such that LOGMULTI*(int(logical) + LOGSHIFT)
!!     always gives unity for a .true. and zero for a .false.
!  istrue = .true.
!  isfalse = .false.
!  TRUE2INT = 0
!  FALSE2INT = 0
!  TRUE2INT = istrue
!  FALSE2INT = isfalse
!  if (TRUE2INT.eq.FALSE2INT) then
!    write(*,*) 'Fatal. Compiler does not support intrinsic conversio&
! &n of logical variables to integers.'
!    call fexit()
!  else
!    LOGSHIFT = -FALSE2INT
!    LOGMULTI = 1.0/dble(TRUE2INT + LOGSHIFT)
!  end if
!
! I/O, command line interface, etc. ...
!
  ilog = std_f_out
  nargs = 0
  ndepravg = 0.0
  nkey = 0
!
#ifdef ENABLE_THREADS
! thread (OpenMP) settings
  thrdat%maxn = 2
#endif

  

  
!
! number of atoms in the system
!
  n = 0
!
! number of molecules in the system
!
  nmol = 0
!
  ! Martin Added 
  
  HSQ_in_progress=.false.
  
! force field parameters 
!
  n_biotyp = 0
  n_ctyp =0 
  n_ljtyp = 0
  epsrule = 2
  sigrule = 1
  reduce = 1.0
  nindx = 12
  nhalf = nindx/2
  primamide_cis_chgshft = 0.0
  dpgrp_neut_tol = 0.0
  be_unsafe = .false.
  ua_model = 0
  bt_patched = .false.
  lj_patched = .false.
  charge_patched = .false.
  fos_patched = .false.
  bonded_patched = .false.
  asm_patched = .false.
  ard_patched = .false.
  nc_patched = .false.
  mass_patched = .false.
  rad_patched = .false.
  improper_conv(1) = 1
  improper_conv(2) = 3
  FOS_PEP_BB(:) = 0.0
  FOS_PEP_PROBB(:) = 0.0
  FOS_PEP_BB_FOR(:) = 0.0
  FOS_PEP_BB_NH2(:) = 0.0
  FOS_PEP_PROBB_FOR(:) = 0.0
  FOS_PEP_CCT(:) = 0.0
  FOS_PEP_UCT(:) = 0.0
  FOS_PEP_CNT(:) = 0.0
  FOS_PEP_PROCNT(:) = 0.0
  FOS_PEP_UNT(:) = 0.0
  FOS_PEP_PROUNT(:) = 0.0
  FOS_ALA(:) = 0.0
  FOS_VAL(:) = 0.0
  FOS_LEU(:) = 0.0
  FOS_ILE(:) = 0.0
  FOS_HIE(:) = 0.0
  FOS_HID(:) = 0.0
  FOS_MET(:) = 0.0
  FOS_SER(:) = 0.0
  FOS_THR(:) = 0.0
  FOS_CYS(:) = 0.0
  FOS_NVA(:) = 0.0
  FOS_NLE(:) = 0.0
  FOS_PHE(:) = 0.0
  FOS_TYR(:) = 0.0
  FOS_TRP(:) = 0.0
  FOS_ASN(:) = 0.0
  FOS_PRO(:) = 0.0
  FOS_GAM(:) = 0.0
  FOS_GLN(:) = 0.0
  FOS_ASP(:) = 0.0
  FOS_GLU(:) = 0.0
  FOS_ARG(:) = 0.0
  FOS_ORN(:) = 0.0
  FOS_DAB(:) = 0.0
  FOS_LYS(:) = 0.0
  FOS_HYP(:) = 0.0
  FOS_ABA(:) = 0.0
  FOS_PCA(:) = 0.0
  FOS_AIB(:) = 0.0
  FOS_HIP(:) = 0.0
  FOS_LINK_CC(:) = 0.0
  FOS_ACA(:) = 0.0
  FOS_NMF(:) = 0.0
  FOS_NMA(:) = 0.0
  FOS_URE(:) = 0.0
  FOS_DMA(:) = 0.0
  FOS_FOA(:) = 0.0
  FOS_PPA(:) = 0.0
  FOS_SPC(:) = 0.0
  FOS_T3P(:) = 0.0
  FOS_T4P(:) = 0.0
  FOS_T5P(:) = 0.0
  FOS_T4E(:) = 0.0
  FOS_CH4(:) = 0.0
  FOS_MOH(:) = 0.0
  FOS_PCR(:) = 0.0
  FOS_NA(:) = 0.0
  FOS_CL(:) = 0.0
  FOS_K(:) = 0.0
  FOS_BR(:) = 0.0
  FOS_I(:) = 0.0
  FOS_CS(:) = 0.0
  FOS_O2(:) = 0.0
  FOS_NH4(:) = 0.0
  FOS_AC(:) = 0.0
  FOS_GDN(:) = 0.0
  FOS_CYT(:) = 0.0
  FOS_URA(:) = 0.0
  FOS_THY(:) = 0.0
  FOS_GUA(:) = 0.0
  FOS_ADE(:) = 0.0
  FOS_PUR(:) = 0.0
  FOS_PRP(:) = 0.0
  FOS_IBU(:) = 0.0
  FOS_NBU(:) = 0.0
  FOS_MSH(:) = 0.0
  FOS_EOH(:) = 0.0
  FOS_EMT(:) = 0.0
  FOS_BEN(:) = 0.0
  FOS_NAP(:) = 0.0
  FOS_TOL(:) = 0.0
  FOS_IMD(:) = 0.0
  FOS_IME(:) = 0.0
  FOS_MIN(:) = 0.0
  FOS_1MN(:) = 0.0
  FOS_2MN(:) = 0.0
  FOS_LCP(:) = 0.0
  FOS_NO3(:) = 0.0
  FOS_NUC_PO4M(:) = 0.0
  FOS_NUC_PO4D(:) = 0.0
  FOS_NUC_RIBO3(:) = 0.0
  FOS_NUC_RIBO2(:) = 0.0
  FOS_NUC_RIBO1(:) = 0.0
  FOS_NUC_MTHF(:) = 0.0
  FOS_XPU(:) = 0.0
  FOS_XPT(:) = 0.0
  FOS_XPC(:) = 0.0
  FOS_XPG(:) = 0.0
  FOS_XPA(:) = 0.0
!
! number of residues, residues with chi-sidechains, and chains in biopolymer sequence
!
  nseq = 0
  chilst%nr = 0 ! see proteus_init
  fylst%nr = 0
  wlst%nr = 0
  nuclst%nr = 0
  nucpuclst%nr = 0
  puclst%nr = 0
  n_constraints = 0
  crlk_mode = 1
!
! # of ring closures for Z-matrix
!
  nadd = 0
!
! control flags for MC engine
!
  do_restart = .false.
  mc_compat_flag = 0
  pdb_analyze = .false.
  pdb_ihlp = 0
  gc_mode = 2
  framecnt = 1
  select_frames = .false.
  use_frame_weights = .false.
  align%yes = .false.
  align%calc = 1
  align%refset = .false.
  align%haveset = .false.
  align%instrmsd = .false.
  align%nr = 0
  align%mmol = 0
  pdb_fileformat = 1
  pdb_mpimany = .false.
  use_pdb_template = .false.
  just_solutes = .false.
  use_trajidx = .false.
  n_pdbunk = 0
  pdb_format(:) = 0
  pdb_writemode = 2
  pdb_readmode = 1
  pdb_hmode = 1
  pdb_nucmode = 1 !CAMPARI-compatible
  pdb_convention(1) = 1 !CAMPARI-style writing
  pdb_convention(2) = 1 !CAMPARI-style reading
#ifdef LINK_XDR
  xtc_prec = 1000.0d0
#endif
#ifdef LINK_NETCDF
  netcdf_ids(:) = 0
  netcdf_fr1 = 0
#endif
  dcd_withbox = .true.
  pdb_force_box = .false.
  pdb_tolerance(1) = 0.8
  pdb_tolerance(2) = 1.25
  pdb_tolerance(3) = 20.0
  nbsr_model = 1
  use_14 = .true.
  mode_14 = 1  ! only real 14 interactions are treated 14
  fudge_el_14 = 1.0
  fudge_st_14 = 1.0
  elec_model = 1
  elec_report = .false.
  tor_report = .false.
  tabul_report = .false.
  ia_report = .false.
  seq_report = .false.
  vdW_report = .false.
  dip_report = .false.
  poly_report = .false.
  drest_report = .false.
  gpc_report = .false.
  fos_report = .false.
  feg_report = .false.
  frz_report = .false.
  psw_report = .false.
  bonded_report = .false.
  sgc_report = .false.
  dyn_report = .false.
  ideal_run = .false.
  is_ev = .false.
  is_lj = .false.
  is_tab = .false.
  is_tablj = .false.
  is_plj = .false.
  is_pewlj = .false.
  is_prflj = .false.
  is_impljp = .false.
  is_implj = .false.
  is_fegplj = .false.
  is_fegprflj = .false.
  is_feglj = .false.
  is_fegev = .false.
  do_n2loop = .true.
  use_IPP = .true.
  use_hardsphere = .false.
  use_IMPSOLV = .false.
  use_attLJ = .false.
  use_CORR = .false.
  use_WCA = .false.
  use_POLAR = .false.
  use_TOR = .false.
  use_POLY = .false.
  use_ZSEC = .false.
! set default to false
  use_ZSEC_MinL = .false.
  use_TABUL = .false.
  use_DREST = .false.
  use_LCTOR = .false.
  use_FEG = .false.
  use_DSSP = .false.

  
  do i=1,4
    use_BOND(i) = .false.
  end do
  guess_bonded = .false.
  use_FEGS_FOS=.false.
  do i=1,MAXENERGYTERMS
    use_FEGS(i) = .false.
  end do
  use_EMICRO = .false.
  use_cutoffs = .false.
  use_mcgrid = .false.
  use_REMC = .false.
#ifdef ENABLE_MPI
  use_MPIAVG = .false.
  use_MPIMultiSeed = .false.
  re_aux(7) = ceiling(0.5*mpi_nodes)
#endif
  fycinput = .false.
!
! boundary conditions
!
  bnd_wrnlmt(:) = 1
  bnd_wrncnt(:) = 0
  bnd_type = 1 !PBC
  bnd_shape = 1 !cube
  bnd_params(1) = 40.0
  bnd_params(2) = 40.0
  bnd_params(3) = 40.0 !side lengths
  bnd_params(4) = -20.0
  bnd_params(5) = -20.0
  bnd_params(6) = -20.0 !lower bnds
  bnd_params(7) = 1.0   !force constant for soft walls
  ens%insV = bnd_params(1)*bnd_params(2)*bnd_params(3)
  bnd_pV = 0.0
!
! computation frequencies/cycles and moveset parameters
!
  mc_acc_crit = 1 ! Metropolis
  nrchis = 2
  particleflucfreq = 0.0
  chifreq = 0.5
  crfreq = 0.0
  torcrfreq_omega = 0.0
  angcrfreq = 0.0
  torcrfreq = 0.0
  puckerfreq = 0.0
  puckerrdfreq = 0.0
  nuccrfreq = 0.0
  maxcrtries = 1
  uj_deadcnt = 0
  dj_deadcnt = 0
  do_deadcnt = 0
  nc_deadcnt = 0
  dj_wrncnt(:) = 0
  do_wrncnt(:) = 0
  nc_wrncnt(:) = 0
  pc_wrncnt(:) = 0
  gc_wrncnt(:) = 0
  fo_wrncnt(:) = 0
  sj_wrncnt(:) = 0
  dj_wrnlmt(:) = 1
  do_wrnlmt(:) = 1
  nc_wrnlmt(:) = 1
  pc_wrnlmt(:) = 1
  gc_wrnlmt(:) = 1
  fo_wrnlmt(:) = 1
  sj_wrnlmt(:) = 1
  uj_totcnt = 0
  dj_totcnt = 0
  do_totcnt = 0
  nc_totcnt = 0
  ujminsz = 3
  ujmaxsz = 8
  UJ_params(1) = 2.0
  UJ_params(2) = 8.0
  UJ_params(3) = 20.0
  UJ_params(4) = 1.0
  UJ_params(5) = 20.0
  UJ_params(6) = 1.0
  torcrmode = 1
  torcrminsz = 1
  torcrminsz2 = 1
  torcrmaxsz = 12
  torcrmaxsz2 = 12
  nuccrminsz = 1
  nuccrmaxsz = 12
  omegafreq = 0.0
  pivot_randfreq = 0.1
  pivot_stepsz = 2.0  ! in deg
  chi_randfreq = 0.1
  chi_stepsz = 10.0   ! in deg
  rigidfreq = 0.0
  rigid_randfreq = 0.1
  clurb_freq = 0.0
  clurb_rdfreq = 1.0   ! NOT the fraction of full randomization amongst clurb-moves!!! (that's covered by rigid_rdfreq)
  clurb_maxsz = 2
  clurb_reset = 500
  clurb_strbail = 0.99
  clurb_discrit = 4.0
  omega_randfreq = 0.4
  omega_stepsz = 2.0  ! in deg
  nucfreq = 0.0
  nrnucim = 2
  nuc_randfreq = 0.4
  nuc_stepsz = 2.0
  nucpuckfreq = 0.0
  nucpuckrdfreq = 0.0
  pucker_distp = 5.0
  pucker_anstp = 2.0
  align_NC = 3   ! shorter tail alignment
  rotfreq = 0.5
  phfreq = 0.0
  trans_stepsz = 5.0  ! in A
  rot_stepsz = 60.0   ! in deg
  otherfreq = 0.0
  other_stepsz = 20.0
  other_randfreq = 0.2
  other_unkfreq = 0.0
  other_natfreq = 0.0
  use_stericgrids = .false.
  use_globmoves = .false.
  use_lctmoves = .false.
  use_coupledrigid = .true.
  do i=1,MAXCHI
    cur_chiflag(i) = .false.
  end do
  do_frz = .false.
  skip_frz = .false.
  in_dyncyc = .false.
  curcyc_start = 1
  curcyc_end = 10000
  max_dyncyclen = 1000
  min_dyncyclen = 500
  max_mccyclen = 500
  min_mccyclen = 200
  first_mccyclen = 50000
  mvcnt%ninsert = 0
  mvcnt%ndelete = 0
  mvcnt%nidentitychange = 0
  mvcnt%nrot = 0
  mvcnt%ntrans = 0
  mvcnt%nrb = 0
  mvcnt%nclurb = 0
  mvcnt%nfyc = 0
  mvcnt%nfy = 0
  mvcnt%nomega = 0
  mvcnt%nre = 0
  mvcnt%nnuc = 0
  mvcnt%nchi = 0
  mvcnt%nnuccr = 0
  mvcnt%npucker = 0
  mvcnt%nujcr = 0
  mvcnt%nsjcr = 0
  mvcnt%ndjcr = 0
  mvcnt%ndocr = 0
  mvcnt%nph = 0
  mvcnt%nother = 0
  mvcnt%nld = 0
  mvcnt%nbd= 0
  mvcnt%nmd = 0
  mvcnt%nlct = 0
  mvcnt%ndynseg = 0
  mvcnt%nmcseg = 0
  mvcnt%nhsq = 0
  mvcnt%avgdynseglen = 0.0
  mvcnt%avgmcseglen = 0.0

  acc%nhsq = 0
  
  HSQ_in_progress=.false.
!
! output and analysis frequencies
!
  enout = 500
  accout = 2000
  ensout = 500
  polcalc = 50
  rhcalc = 500
  angcalc = 50
  intcalc = 250
  dsspcalc = 1000
  segcalc = 50
  sctcalc = 10000
  holescalc = huge(holescalc)
  pccalc = 500
  savcalc = 500
  savreq%nats = 0 
  savreq%instfreq = huge(savcalc)
  particlenumcalc = 100
  polout = 500
  torout = 10000
  toroutmode = 1
  xyzout = 10000
  rstout = 10000
  xyzmode = 2
  covcalc = huge(covcalc)
  dipcalc = huge(dipcalc)
  torlccalc = huge(torlccalc)
  torlcmode = 1
  phout = huge(phout)
!
! molecular volume variables
  sysvol(1) = 0.0
!
! character input
  paramfile = '/usr/local/campari/params/abs3.2_opls.prm'
  seqfile = '/usr/local/campari/examples/test.in'
  prefsamplingfile = '/usr/local/campari/examples/yourcalc.pfs'
  angrpfile = '/usr/local/campari/examples/test_angrps.in'
  bbsegfile = '/usr/local/campari/data/bbseg.dat'
  besselfile = '/usr/local/campari/data/BesselFunctions.dat'
  griddir = '/usr/local/campari/data/grids/'
  cm_dir = '/usr/local/campari/data/'
  pdb_suff(1:4) = '.pdb'
  pdb_pref(1:9) = 'yourcalc_'
  pdbinfile = '/usr/local/campari/examples/yourcalc_test.pdb'
  xtcinfile = '/usr/local/campari/examples/yourcalc_traj.xtc'
  dcdinfile = '/usr/local/campari/examples/yourcalc_traj.dcd'
  netcdfinfile = '/usr/local/campari/examples/yourcalc_traj.nc'
  pdbtmplfile = '/usr/local/campari/examples/yourcalc_tmpl.pdb'
  torfile = '/usr/local/campari/examples/test_tor.in'
  pccodefile = '/usr/local/campari/examples/test.idx'
  tabcodefile = '/usr/local/campari/examples/test.idx'
  tabpotfile = '/usr/local/campari/examples/test_pot.dat'
  tabtangfile = '/usr/local/campari/examples/test_tang.dat'
  torlcfile = '/usr/local/campari/examples/test.cof'
  torlcfile2 = '/usr/local/campari/examples/invtest.cof'
  lctpotfile = '/usr/local/campari/examples/test_lct.dat'
  polfile = '/usr/local/campari/examples/test_poly.in'
  drestfile = '/usr/local/campari/examples/test_drest.in'
  frzfile = '/usr/local/campari/examples/test_frz.in'
  align%filen = '/usr/local/campari/examples/test_align.in'
  fegfile = '/usr/local/campari/examples/test_feg.in'
  particleflucfile = '/usr/local/campari/examples/test_pfluc.in'
  tstat%fnam = '/usr/local/campari/examples/test_tgrps.in'
  cart_cons_file = '/usr/local/campari/examples/test_shake.in'
  cpatchfile = '/usr/local/campari/examples/yourcalc.cpi'
  nblfilen = '/usr/local/campari/examples/yourcalc_nbl.nc'
  emmapfile = '/usr/local/campari/examples/test_input_density_map.nc'
  cfilen = '/usr/local/campari/examples/clustering.idx'
  tbrkfilen = '/usr/local/campari/examples/yourcalc.tbr'
  fospatchfile = '/usr/local/campari/examples/yourcalc.fpi'
  asmpatchfile = '/usr/local/campari/examples/yourcalc.spi'
  ardpatchfile = '/usr/local/campari/examples/yourcalc.xpi'
  bpatchfile = '/usr/local/campari/examples/yourcalc.bpi'
  ljpatchfile = '/usr/local/campari/examples/yourcalc.lpi'
  masspatchfile = '/usr/local/campari/examples/yourcalc.mpi'
  radpatchfile = '/usr/local/campari/examples/yourcalc.rpi'
  biotpatchfile = '/usr/local/campari/examples/yourcalc.tpi'
  savreq%filen = '/usr/local/campari/examples/yourcalc.sri'
  basename = 'yourcalc'
  call strlims(basename,t1,t2)
  bleng = t2-t1+1
!     
! step+term counters
!
  nstep = 0
  dihed_wrncnt = 0
  dihed_wrnlmt = 1
  dblba_wrncnt = 0
  dblba_wrnlmt = 1
  acc%nfyc = 0
  acc%nsjcr = 0
  acc%nujcr = 0
  acc%nchi = 0
  acc%nfy = 0
  acc%ndjcr = 0
  acc%ndocr = 0
  acc%nre = 0
  acc%nlct = 0
  acc%nrb = 0
  acc%nclurb = 0
  acc%nrot = 0
  acc%ntrans = 0
  acc%nomega = 0
  acc%nnuc = 0
  acc%nnuccr = 0
  acc%ninsert = 0
  acc%ndelete = 0
  acc%nidentitychange = 0
  acc%npucker = 0
  acc%naccept = 0
  acc%nph = 0
  acc%nother = 0
  npucmoves = 0
!
  nsavavg = 0
  time_analysis = 0.0
  time_energy = 0.0
  time_struc = 0.0
  time_comm = 0.0
  time_ph = 0.0
  time_holo = 0.0
  time_pme = 0.0
  time_nbl = 0.0
#ifdef ENABLE_MPI
  mpi_granularity = 1
  use_MPIcolls = .false.
  mpi_cnt_en = 0
  mpi_cnt_tor = 0
  mpi_cnt_pol = 0
  mpi_cnt_xyz = 0
  mpi_cnt_sav = 0
  mpi_cnt_trcv = 0
  mpi_cnt_ens = 0
#endif
!
! energy variables
!
  esave = 0.0
  esaver = 0.0
  esavec = 0.0
  do i=1,MAXENERGYTERMS
    esterms(i) = 0.0
  end do
!
! general settings
!
  nsim = 100000
  nequil = 10000
  fycxyz = 1
  dyn_dt = 2.0e-3
  cart_cons_grps = 0
  cart_cons_method = 1
  cart_cons_mode = 1
  cart_cons_source = 1 ! from Engh-Huber et al.
  no_shake = .false.
  add_settle = .true.
  shake_tol = 1.0e-4
  shake_atol = 1.0e-4
  lincs_order = 8
  shake_maxiter = 1000
  cs_wrncnt(:) = 0
  cs_wrnlmt(:) = 1
  use_dyn = .false.
  grad_check = .false.
  fudge_mass = .false.
  dyn_mode = 1
  dyn_integrator = 2
  dyn_integrator_ops(:) = 0
  dyn_integrator_ops(10) = 1 ! what to do with unsupported d.o.f.s
  kelvin = 298.0
  fos_Tdepref = kelvin
  use_Tdepfosvals = .false.
  invtemp = 1.0d0 /(kelvin*gasconst)
  fric_ga = 1.0  ! in 1/ps
  ens%avgP = 0.0
  ens%avgT = 0.0
  ens%avgK = 0.0
  ens%avgK2 = 0.0
  ens%avgU = 0.0
  ens%avgV = 0.0
  ens%avgPT(:,:) = 0.0
  ens%avgVirT(:,:) = 0.0
  ens%avgP = 0.0
  ens%avgR(:) = 0.0
  extpress = 1.0 ! in bar
  invcompress = 23.0 ! vacuum
  ens%avgcnt = 0
  ens%flag = 1 ! NVT
  tstat%flag = 4 ! Bussi-Parrinello
  tstat%params(1) = 1.0 ! in ps
  tstat%n_tgrps = 1
  pstat%flag = 1 ! droplet boundary
  pstat%params(1) = 10.0 ! in g/mol
  ens%sysfrz = 1 ! no com translation removal or rotation removal
  ionicstr = 0.0d0
  invdhlen = 0.0d0
  nsancheck = 10000
  mcnb_cutoff = 10.0
  mcel_cutoff = 12.0
  mcnb_cutoff2 = mcnb_cutoff*mcnb_cutoff
  mcel_cutoff2 = mcel_cutoff*mcel_cutoff
  imcel2 = 1.0/mcel_cutoff2
  nbl_up = 5
  nbl_maxnbs = 1000
  lr_up = 10
  lrel_md = 4
  rf_mode = 1
  ewpm_pre = -0.1
  ewpm = 0.1
  ewpm2 = ewpm*ewpm
  ewpite = 2.0*ewpm/sqrt(PI)
  ewfspac = 1.0
  ewald_mode = 1 
  splor = 8
  ewcnst = 0.0
  lrel_mc = 3
  use_stericscreen = .false.
  screenbarrier = 10000000
  scale_attLJ = 0.0
  scale_IMPSOLV = 0.0
  scale_CORR = 0.0
  scale_IPP = 1.0
  scale_WCA = 0.0
  scale_POLAR = 0.0
  scale_TOR = 0.0
  scale_POLY = 0.0
  scale_ZSEC = 0.0
  scale_TABUL = 0.0
  scale_DREST = 0.0
  scale_LCTOR = 0.0
  scale_DSSP = 0.0
  do i=1,MAXENERGYTERMS
    scale_FEGS(i) = 0.0
  end do
  do i=1,4
    scale_BOND(i) = 0.0
  end do
  scale_EMICRO = 0.0
!
  par_IMPSOLV(1) = 5.0 !diameter of first solvation shell
  par_IMPSOLV(2) = 78.0 !dielectric of water
!  par_IMPSOLV(3) = 1.0  !fraction of water molecule used to get atsavminfr
!  par_IMPSOLV(4) = 1.0  !fraction of water molecule used to get atsavminfr2
  par_IMPSOLV(3) = 0.25  !interpolation constant for sigmoidal function (FOS)
  par_IMPSOLV(4) = 0.5  !interpolation constant for sigmoidal function (q-SCR)
  scrq_model = 2        !choice of model for charge-screening
  i_sqm = 0             !choice of mean for generalized charge-screening (modes 5/6)
  par_IMPSOLV(6) = 0.1  !center for sigmoidal function (relative to allowed interval) (FOS)
  par_IMPSOLV(7) = 0.9  !center for sigmoidal function (relative to allowed interval) (q-SCR) 
  par_IMPSOLV(5) = 0.2595
  coul_scr = 1.0 - 1.0/sqrt(par_IMPSOLV(2)) !the resulting screening factor (see energy.f)
  par_IMPSOLV(8) = 1./5. !contact dielectric (only relevant for screening model 3/4)
  par_IMPSOLV(9) = 0.5  !how much distance-dependence to mix in (from 0.0 to 1.0) -> only relevant for model 3
  par_IMPSOLV(10) = 0.0 ! SAV-response compressibility 
  par_IMPSOLV(11) = 0.0 ! SEXV-response compressibility
!
  par_WCA(1) = 1.5 !cutoff in sigma
  par_WCA(2) = 0.0 !attractive lambda
  par_WCA(3) = pi/(par_WCA(1)**2 - ROOT26*ROOT26) !attr. param. 1
  par_WCA(4) = pi - par_WCA(3)*ROOT26*ROOT26      !attr. param. 2
  tbp%ipprm(:) = 0.0 ! Catmull Rom
!
  polbiasmode = 1
!
! #tyler pka
  ! MArtin : 
  do_pka = .false.
  do_pka_2 = .false.
  do_pka_3 = .false.
  do_hs=.false.
  use_softcore= .false.
  sc_alpha = 0.5
  sc_n = 1
  par_top=.false.
  print_det=.false.
  sav_det_freq=0
  sav_det_mode=1
    !Martin added
  FEG_FOS_OFFSET=.false.
  FEG_FOS_OFFSET_FACT=10./12
  debug_pka=.false.
  HS_in_progress=.false.
  HSQ_steps =50
  hs_mode=3
 ! pka_limits_scaled = .false.
!
  do i=1,MAXENERGYPARAMS
    par_ZSEC(i) = 0.0
  end do
! #tyler set the default required length for zsec
  par_ZSEC_L = 1

  par_ZSEC2(1) = -60.0     !center of alpha-circle
  par_ZSEC2(2) = -50.0     !
  par_ZSEC2(3) = 35.0      !radius of alpha-circle
  par_ZSEC2(4) = 0.002     !stiffness around alpha-plateau
  par_ZSEC2(5) = -155.0    !center of beta-circle
  par_ZSEC2(6) = 160.0     !
  par_ZSEC2(7) = 35.0      !radius of beta-circle
  par_ZSEC2(8) = 0.002     !stiffness around beta-plateau
  par_DSSP(7) = 0.45       !target H-score
  par_DSSP(9) = 0.45       !target E-score
  par_DSSP(8) = 0.0        !harmonic H-restraint
  par_DSSP(10) = 0.0       !harmonic E-restraint
!
  fegcbmode = 1            !linear (2 = vanilla soft-core)
  fegljmode = 2            !vanilla soft-core (1 = linear, 3 = exp. soft-core)
  fegmode = 1
  par_FEG2(3) = 0.5        !the soft-core limit for IPP/LJ
  par_FEG2(4) = 0.5        !the soft-core limit for Cb
  par_FEG2(7) = 2.0        !exponent in the SC-term (LJ)
  par_FEG2(8) = 1.0        !exponent for global polynomial/exp. scaling (LJ)
  par_FEG2(11) = 1.0       !exponent in the SC-term (Cb)
  par_FEG2(12) = 1.0       !exponent for global polynomial/exp. scaling (Cb)
  par_LCTOR(1) = 0.1       !default resolution
  par_LCTOR(2) = -4.95     !default value for first bin
  par_LCTOR(3) = 245.0     !default value for artificial wall height
  par_LCTOR2(1) = 100      !default number of bins
  do i=1,MAXTORLC
    do j=1,MAXTORLCBINS
      lct_pot(i,j) = 0.0
    end do
  end do
  harappa = .false.
!
  globcyclic = .false.
  globrandomize = 0
  globrdmatts = 300
  globrdmthresh = 1.0e4
  mini_mem = 10
  mini_mode = 1
  mini_econv = 0.5
  mini_stepsize = 0.1
  mini_uphill = 0.1
  mini_xyzstep = 0.1 ! in A
  mini_rotstep = 0.5 ! in deg
  mini_intstep = 0.5 ! in deg
  mini_sc_sdsteps = 0
  mini_sc_tbath = 0.1 ! Kelvin
  mini_sc_heat = 0.1
  cm_splor = 4
!
! steric grids and HC (biased MC) settings
!
  stgr%halfwindow = 5.0
!
! grid settings
!
  grid%dim(1) = 10
  grid%dim(2) = 10
  grid%dim(3) = 10
  grid%origin(1) = -20.0
  grid%origin(2) = -20.0
  grid%origin(3) = -20.0
  grid%deltas(1) = 4.0
  grid%deltas(2) = 4.0
  grid%deltas(3) = 4.0
  grid%maxgpnbs = 125
  grid%maxnbs = 200
  grid%report = .false.
!
! concerted rotation settings
!
  cr_mode = 1
  cr_a = 1.0/(PI/80.0)
  nr_crdof = 8
  nr_crres = 4
  nr_crfit = 9
  cr_b = 1.0
  call consmax()
!
! WL settings
!
  do_wanglandau = .false.
  debug_wanglandau = .false.
  finit_wanglandau = .false.
  wld%wl_mode = 1 ! 1 = generate g(E), 2 = generate g(RG), 3 = generate g(E) but with temperature bias
  wld%dimensionality = 1
  wld%exrule = 1 ! do not grow histograms dynamically
  wld%buffer = max(1,nsim/10)
  wld%freeze = .true.
  wld%wl_mol(:) = 1 ! only relevant if wld%wl_mode == 2
  wld%wl_rc(:) = 1 
  wld%fval = 1.0  !WL simulation convergence value initial value (LOG(f))
  wld%gh_flatcheck_freq = 10000 !how often to check if reference energy histogram is flat
  wld%hvmode = 1 !how many times to visit each energy bin in order to update the F value 
  wld%t1on = .false. !activate 1/t stage of the algorithm
  wld%stepnum = 0    !total actual WL moves (not = istep)
  wld%stepnumbuf = 0
  wld%flevel = 0 ! number of "stages" (gh resets) before using 1/t algorithm
  wld%hufreq = 10  !Decorrelation: update histograms with this frequency
! Various WL statistics
  wld%overflow = 0  !counts that went over max bin window
  wld%underflow = 0  !counts that went under bin window
  wld%accepted = 0 ! Total WL moves accepted
  wld%maxb = 1
  wld%maxb2d(:) = 1
  wld%nbins = 0
  wld%nbins2d(:) = 0
  wld%g_binsz = 0.1
  wld%g_binsz2d(:) = 0.1
  wld%g_max = 0.0
  wld%g_max2d(:) = 0.0
! fileIO
  wld%use_ginitfile = .false.
!
! ELS settings
  do_accelsim = .false.
  hmjam%prtfrmwts = HUGE(hmjam%prtfrmwts)
  hmjam%threshwt = 0.0
  hmjam%nlst = 0
!
! RE settings
!
#ifdef ENABLE_MPI
  re_nbmode = 2
  re_tryswap = 1
  re_freq = 1000
  inst_retr = .true.
  re_conditions = 2
  re_conddim = 1
  do j=1,MAXREDIMS
    re_types(j) = 1
  end do
  nreolap = 0
  re_olcalc = 1000
  inst_reol = 0
  reol_all = .false.
  force_rexyz = .false.
  force_singlexyz = .true.
  re_velmode = 1
  re_aux(:) = 0
  re_aux(4:5) = 10000
  re_infile = '/usr/local/campari/examples/yourcalc.rex'
  re_traceinfile = '/usr/local/campari/examples/yourcalc.rxt'
#endif
!
! DSSP analysis
!
  bdgtb%nrs = 0
!
! SASA analysis
!
  savprobe = 1.5 ! = par_IMPSOLV(1)/2.0 (see above)
  savavg = 0.0
!
! contact analysis
!
  contactcalc = 500
  contact_off = 0
  clucalc = 100
  contact_cuts(1) = 5.0*5.0
  contact_cuts(2) = 5.0*5.0
!
! angular maps, segment distributions, general torsion stuff
!
  nrspecrama = 0
  fyres = 10.0
  ntorlcs = 0
  torlc_params(1) = -5.0
  torlc_params(2) = 0.1
  do_ints(:) = .false.
  do_ints(4) = .true.
  intres(1) = 0.01
  intres(2) = 1.0
  intres(3) = 1.0
  intres(4) = 3.6
!
! Rg-delta* and PC analysis
!
  rg_binsz = 0.2
  tp_binsz = 0.01
  dlst_binsz = 0.01
  maxrgbins = 1000
  pc_binsz = 0.2
  tp_exp = 1./3.    ! universal scaling exponent for the globule in 3D
  do_amidepc = .false.
!
! Diffraction analysis
!
  diffrcalc = huge(diffrcalc)
  bes_max = 10
  rr_max = 100
  zz_max = 100
  rr_res = 0.1
  zz_res = 0.1
  bes_res = 0.25
  bes_bins = 2001
  bes_num = 501
  diffr_cnt = 0
  diffr_uax = .false.
  do i=1,3
    diffr_axon(i) = 0.0
  end do
!
! GPC settings
!
  inst_gpc = 0
  n_pccalls = 0
  gpc%mode = -1
  gpc%atms = 0
  gpc%nos = 0
  gpc%nlst = 0
!
! some more analysis settings
!
! scatter analysis
  qv_res = 0.025
  nqv = 20
!
! initialize whether holes calculation has been done
  firsthole = .true.
!
! initialize whether to output instantaneous values of cosines (for 
! calculation of persistence length) instead of just a final
! average
  inst_pers = .false.
!
! DSSP-settings
  inst_dssp = .false.
  par_DSSP(1) = -2.5
  par_DSSP(2) = -0.5
  par_DSSP(3) = -4.0
  par_DSSP(4) = 10.0
  par_DSSP2(1) = 3
  par_DSSP2(2) = 3
!
! initialize distance restraint parameters
  ndrest = 0
  allndr = 0
!
! covariance stuff
  covmode = 1
!
! structural clustering
  cmode = 1 ! LEADER clustering
  cdis_crit = 5 ! standard all-atom RMSD with alignment 
  cstored = 0
  cstorecalc = HUGE(cstorecalc)
  cchangeweights = 0
  cprepmode = 0
  cwcombination = 1 ! arithmetic mean
  cwacftau = 100
  csmoothie = 10
  cwwindowsz = 1000
  cdynwbuf = 1.0
  cradius = 2.0
  cmaxrad = 2.0*cradius
  chardcut = 3.0
  csivmin = 19
  csivmax = 49
  read_nbl_from_nc = .false.
  align_for_clustering = .true.
  c_nhier = 8
  clinkage = 1
  refine_clustering = .false.
  cleadermode = 1
  pcamode = 1
  reduced_dim_clustering = 0
  cprogindex = 2 ! use approximate scheme
  cprogindrmax = 50 ! max guesses for approximate scheme
  cprogindstart = 1 ! starter snap for prog index
  cprogpwidth = 1000 ! width for localized cut
  ccfepmode = 0 ! disabled
  cprogfold = 0 ! disabled
  cequil = 0 ! do not split and reequilibrate MSM
  itbrklst = 1
  ntbrks = 0
  ntbrks2 = 0
!
! EM analysis
  emcalc = HUGE(emcalc)
  emgrid%dim(:) = 40
  emgrid%origin(:) = -20.0
  emgrid%deltas(:) = 0.0
  emgrid%cnt = 0
  embgdensity = 1.0
  emheuristic = 3
  emsplor = 3
  emthreshdensity = 0.0
  emtotalmass = -1.0
  emflatval = HUGE(emflatval)
  empotmode = 1 ! instant., 2=time averaged
  emiw = 0.1 ! weight of instantaneous term in time-averaged method
  empotprop = 1 ! mass, 2 = atomic number
  emtruncate =-1.0*(HUGE(emtruncate))/1000.0
  emoptbg = HUGE(emoptbg)
  emmcacc = 0
  curmassm = .false.
  emgrid%idx(:,1) = 1
  emgrid%idx(:,2) = emgrid%dim(:)
  fdlts(:) = 0.0
!
! grab command line arguments
  call grab_args()
!
end
!
!---------------------------------------------------------------------------------
