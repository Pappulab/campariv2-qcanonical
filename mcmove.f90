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
! CONTRIBUTIONS: Hoang Tran, Xiaoling Wang, Rohit Pappu, Adam Steffen      !
!                Albert Mao                                                !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!-----------------------------------------------------------------------------
! 
!subroutine keep_track_charge()
!
!end subroutine 
!  
subroutine do_HS_walk(istep) 
    ! This is a reimplementation of the HS method, since I haven't been able to find the keywords associated with it 
    ! A few notes : I think I should be making it so the moves can only be applied the the residue that have a different 
    ! Hamiltoninan. This should be as an option, since it only really makes sense for the salt 
    ! I can still use the HS_walk method
    ! The other thing is whether I should add an option to force the acceptance by trying to get out pf the auxilary hamiltonian 
    ! IN that case I could get away with it as long as I set a meximum length for the auxilary markov chain 
    use energies
    use movesets
    use atoms
    use polypep
    use sequen
    use system
    use mcsums
    use accept
    use zmatrix 
    use units
    use iounit
    use mpistuff
    
    
    
    
    implicit none
    
    character WWW*200000,WW_tmp*1000,my_format*10
    character*1 tab
    integer rs,i
    RTYPE energy3,energy,edum(MAXENERGYTERMS),temp,expterm,random,save_es(MAXENERGYTERMS),kelvin_save
    RTYPE sqrt
    RTYPE evec_11(4),evec_12(4),evec_21(4),evec_22(4)
    RTYPE, ALLOCATABLE::  x_save(:),y_save(:), z_save(:)
    
    integer, intent(in) :: istep
    logical moveok
    
    mvcnt%nhs=mvcnt%nhs+1
    evec_11(1)=esterms(1)
    evec_11(2)=esterms(3)
    evec_11(3)=esterms(4)
    evec_11(4)=esterms(6)
    
    allocate(x_save(n))
    allocate(y_save(n))
    allocate(z_save(n))
    
    allocate(ztor_save(size(ztor)))
    allocate(blen_save(size(blen)))
    allocate(bang_save(size(bang)))
    
    x_save(:)=x(:)
    y_save(:)=y(:)
    z_save(:)=z(:)
    
    ztor_save(:)=ztor(:)
    blen_save(:)=blen(:)
    bang_save(:)=bang(:)
    
    save_es(:)=esterms(:) 
    save_E_HS(1,1)=esave 
    par_hsq(:,:)=0. ! Just to be safe
    do i=1,size(hs_residues)
        par_hsq(hs_residues(i),1)=(1.-hs_POLAR)
        par_hsq(hs_residues(i),2)=1.! Because this is taken care of in the inner loop routine
        par_hsq(hs_residues(i),3)=(1.-hs_FOS)
        par_hsq(hs_residues(i),4)=0. ! We don't need to use the third limit
    end do 
    call scale_stuff_hsq()
    save_E_HS(2,1)=energy(esterms)

    HS_in_progress=.true.
    
    call HS_walk() ! Now if this move is true, then I should be able
    
    HS_in_progress=.false.
    save_E_HS(2,2)=energy(esterms) !
    par_hsq(:,:)=0.
    call scale_stuff_hsq()

    save_E_HS(1,2)=energy(esterms) ! Get and save the new energy
    evec_12(1)=esterms(1)
    evec_12(2)=esterms(3)
    evec_12(3)=esterms(4)
    evec_12(4)=esterms(6)
    expterm = exp(-(save_E_HS(1,2)-save_E_HS(1,1))*invtemp)*exp(-(save_E_HS(2,1)-save_E_HS(2,2))*invtemp)! Evalutate
    
    temp=random()
    
      if (expterm.gt.temp) then
        moveok = .true.
        acc%nhs=acc%nhs+1
        BATH_ENER=0
        esave=energy(esterms)

      else
        moveok = .false.
        x(:)=x_save
        y(:)=y_save
        z(:)=z_save
        ztor(:)=ztor_save(:)
        blen(:)=blen_save(:)
        bang(:)=bang_save(:)
        call scale_stuff_hsq()
        call switch_bonds()
        esave=energy(esterms)  ! Resets the energy to a no cut-off value
        esterms(:)=save_es(:) ! restore energy vectors
!        print *,"Rejected HS"
        
      end if
    par_hsq(:,:)=0.0
    deallocate(ztor_save)
    deallocate(blen_save)
    deallocate(bang_save)
      
    deallocate(x_save)
    deallocate(y_save)
    deallocate(z_save)
      
      
end subroutine
!
!
!
subroutine do_HSQ_walk(istep) !This subroutine will contain all the steps of the HSQ 
!    First, find which groups to switch.
    use energies
    use movesets
    use atoms
    use polypep
    use sequen
    use system
    use mcsums
    use accept
    use zmatrix 
    use units
    use iounit
    use inter
    
    use mpistuff
    
    implicit none
    
    character WWW*200000,WW_tmp*1000,my_format*10
    character*1 tab
    integer rs,i
    RTYPE energy3,energy,edum(MAXENERGYTERMS),temp,expterm,random,save_es(MAXENERGYTERMS),kelvin_save
    RTYPE sqrt
    
    RTYPE evec_11(4),evec_12(4),evec_21(4),evec_22(4)
    
    integer, intent(in) :: istep
    logical moveok
    
    tab = char(9)
    edum(:)=0.0
    
    save_es(:)=esterms(:) 

!    If I just set esterms to 0, the the energy is recomputed from scratch
    !hsq_move=.false. ! Put that here for noral execution of mc move in 
    ! From here, decided which residue to change the state of 
    ! Martin : here I reuse thelimits variables to save the state before the HS walk
    ! Do here I need to save :
!       -the par_pka
!       -the energy (to save time) 
!       -the coordinates
    save_E_HS(1,1)=esave !energy3(edum,.true.) ! gets the original energy
    
    evec_11(1)=esterms(1)
    evec_11(2)=esterms(3)
    evec_11(3)=esterms(4)
    evec_11(4)=esterms(6)
    
    x_limits(:,1)=x(:)! This is just reusing the xlimits array
    y_limits(:,1)=y(:)
    z_limits(:,1)=z(:)
    
    ztor_save(:)=ztor(:)
    blen_save(:)=blen(:)
    bang_save(:)=bang(:)
    
    par_hsq_save(:,:)=par_hsq(:,:)
    seq_q_state_save(:)=seq_q_state(:)
    
    if (nhis.ne.0) then 
        his_state_save(:,:)=his_state(:,:)
    end if 
    if (hsq_mode.eq.1) then 
        call HS_q_move_tree_cste_Q()
    else if (hsq_mode.eq.2) then 
        call HS_q_move_tree_cste_FOS() !This mode is not currently working 
    else 
        print *,"Fatal. You have to specify a mode if you are going to use HSQ."
        call fexit()
    end if 
    
    acc%HSQ_state_tries(acc%curr_tries(1),acc%curr_tries(2))=acc%HSQ_state_tries(acc%curr_tries(1),acc%curr_tries(2))+1
    
    
    
    seq_q_state(:)=seq_q_state_save(:)+q_switch(:)

    !Selection of the alternative Markov Chain Hamiltonian parameters
    do rs=1,nseq ! Here my goal is to enhance the acceptance ratio V3
        if ((q_switch(rs).ne.0).and.((seq_q_state(rs).ne.3).or.((seq_q_state_save(rs).ne.3)))) then 
            ! What I need to do is to find a state in which the residue in which the proton is grown is solvent exposed. 
            ! But also get the charges 
            !q_switch(i)*hsq_sign(i)*-1 ! This is equal to one if proton is grown, -1 if it is lost, 0 if nothing.
            call select_par_hsq(alt_mode,rs)
        end if
   !     par_hsq(rs,3)=par_hsq_save(rs,3)+q_switch(rs) ! Sets the FOS to the final value,(try to get charges less close)
   !     par_hsq(rs,4)=par_hsq_save(rs,4)+q_switch(rs)
    end do 
    
    ! Histidine correction 
    do i=1,nhis
            ! if the residue is in fact changing 
        if (seq_q_state(his_state(i,1)).eq.3) then 
            his_state(i,2)=1  
        else if (seq_q_state(his_state(i,1)).eq.2) then 
            his_state(i,2)=0                
        else if (seq_q_state(his_state(i,1)).eq.1) then             
            his_state(i,2)=0                
        end if 
        call select_par_hsq(alt_mode,his_state(i,1))
    end do 
    
    BATH_ENER=ener_bath 
    
    call scale_stuff_hsq()
    call switch_bonds()

    save_E_HS(2,1)=energy(esterms) ! Get and save the new energy
    
    esave=save_E_HS(2,1)
    
    evec_21(1)=esterms(1)
    evec_21(2)=esterms(3)
    evec_21(3)=esterms(4)
    evec_21(4)=esterms(6)
    
    HSQ_in_progress=.true.
    
    call HS_walk() ! Now if this move is true, then I should be able
    
    HSQ_in_progress=.false.
    
    save_E_HS(2,2)=energy(esterms) ! Get the final alternative hamiltonian energy
    evec_22(1)=esterms(1)
    evec_22(2)=esterms(3)
    evec_22(3)=esterms(4)
    evec_22(4)=esterms(6)
    
    do rs=1,nseq
        if (seq_q_state(rs).eq.1) then 
            par_hsq(rs,1)=0.
            par_hsq(rs,2)=0.
            par_hsq(rs,3)=0.
            par_hsq(rs,4)=0.
        else if ((seq_q_state(rs).eq.2).or.(seq_q_state(rs).eq.3)) then 
            par_hsq(rs,1)=1.
            par_hsq(rs,2)=1.
            par_hsq(rs,3)=1.
            par_hsq(rs,4)=1.
        end if 
    end do

    do i=1,nhis
        if (seq_q_state(his_state(i,1)).eq.1) then 
            his_state(i,2)=0
        else if (seq_q_state(his_state(i,1)).eq.2) then 
            his_state(i,2)=0
        else if (seq_q_state(his_state(i,1)).eq.3) then  
            his_state(i,2)=1
        end if 
    end do 

    BATH_ENER=0
    
    call scale_stuff_hsq()
    call switch_bonds()

    save_E_HS(1,2)=energy(esterms) ! Get and save the new energy
    evec_12(1)=esterms(1)
    evec_12(2)=esterms(3)
    evec_12(3)=esterms(4)
    evec_12(4)=esterms(6)
    expterm = exp(-(save_E_HS(1,2)-save_E_HS(1,1))*invtemp)*exp(-(save_E_HS(2,1)-save_E_HS(2,2))*invtemp)! Evalutate
    
    temp=random()
    
      if (expterm.gt.temp) then
        moveok = .true.
        acc%nhsq=acc%nhsq+1

        BATH_ENER=0
        esave=energy(esterms)
        
        acc%HSQ_state_acc(acc%curr_tries(1),acc%curr_tries(2))=acc%HSQ_state_acc(acc%curr_tries(1),acc%curr_tries(2))+1
         
        
        if (debug_pka.eqv..true.) then 
            print *,"HSQ Accepted at ",istep
        end if 
      else
        moveok = .false.
        
        par_hsq(:,:)=par_hsq_save(:,:)
        seq_q_state(:)=seq_q_state_save(:)
        if (nhis.ne.0) then 
            his_state(:,:)=his_state_save(:,:)
        end if 
        acc%curr_tries(2)=acc%curr_tries(1)
        
        x(:)=x_limits(:,1)
        y(:)=y_limits(:,1)
        z(:)=z_limits(:,1)
        
        ztor(:)=ztor_save(:)
        blen(:)=blen_save(:)
        bang(:)=bang_save(:)
        
        call scale_stuff_hsq()
        call switch_bonds()
        
        esave=energy(esterms)  ! Resets the energy to a no cut-off value
        
        esterms(:)=save_es(:) ! restore energy vectors
        if (debug_pka.eqv..true.) then 
            print *,"HSQ rejected at ",istep
        end if 
      end if
      
    
      

    if ((mvcnt%nhsq.ne.0))then
!        write(ilog,*) "so far the acceptance ratio",FLOAT(acc%nhsq)/mvcnt%nhsq
!        write(ilog,*) "ratio_alt ",exp(-(save_E_HS(2,1)-save_E_HS(2,2))*invtemp)
!        write(ilog,*) "detail_alt",save_E_HS(2,1),save_E_HS(2,2)
!        write(ilog,*) "det 21",evec_21(:)
!        write(ilog,*) "det 22",evec_22(:)
!        write(ilog,*) "ratio_rel",exp(-(save_E_HS(1,2)-save_E_HS(1,1))*invtemp)
!        write(ilog,*) "detail_rel",save_E_HS(1,1),save_E_HS(1,2)
!        write(ilog,*) "det 11",evec_11(:)
!        write(ilog,*) "det 12",evec_12(:)
!        flush(ilog)
    end if
    flush(6)
    
end subroutine

subroutine select_par_hsq(alt_mode,rs) 
    use energies
    use atoms
    integer alt_mode,rs
    ! Mode 1 : Half charge, atom fully grown, full charged FOS
    ! Mode 2 : No charge, atom fully grown, full charged FOS
    
    if (alt_mode.eq.1) then 
        if ((q_switch(rs)*hsq_sign(rs)*-1).eq.1) then
            par_hsq(rs,1)=1.
            par_hsq(rs,2)=0.5
            par_hsq(rs,3)=0.
            par_hsq(rs,4)=1.
        else if ((q_switch(rs)*hsq_sign(rs)*-1).eq.-1) then 
            par_hsq(rs,1)=0.
            par_hsq(rs,2)=0.5
            par_hsq(rs,3)=1.
            par_hsq(rs,4)=0.
        end if

    else if (alt_mode.eq.2) then 
        if ((q_switch(rs)*hsq_sign(rs)*-1).eq.1) then
            par_hsq(rs,1)=1.
            par_hsq(rs,2)=0.
            par_hsq(rs,3)=0.
            par_hsq(rs,4)=1.
        else if ((q_switch(rs)*hsq_sign(rs)*-1).eq.-1) then 
            par_hsq(rs,1)=0.
            par_hsq(rs,2)=0.
            par_hsq(rs,3)=1.
            par_hsq(rs,4)=0.
        end if 

    else if (alt_mode.eq.3) then 
        if ((q_switch(rs)*hsq_sign(rs)*-1).eq.-1) then 
            par_hsq(rs,1)=1.
            par_hsq(rs,2)=1.
            par_hsq(rs,3)=0.
            par_hsq(rs,4)=1.
        else if ((q_switch(rs)*hsq_sign(rs)*-1).eq.1) then 
            par_hsq(rs,1)=0.
            par_hsq(rs,2)=0.
            par_hsq(rs,3)=1.
            par_hsq(rs,4)=0.
        end if 

    else if (alt_mode.eq.4) then 
        if ((q_switch(rs)*hsq_sign(rs)*-1).eq.-1) then
            par_hsq(rs,1)=1.
            par_hsq(rs,2)=0.
            par_hsq(rs,3)=0.
            par_hsq(rs,4)=1.
        else if ((q_switch(rs)*hsq_sign(rs)*-1).eq.1) then 
            par_hsq(rs,1)=0.
            par_hsq(rs,2)=1.
            par_hsq(rs,3)=1.
            par_hsq(rs,4)=0.
        end if 
        
        
        
        
    else 
        print *,"FATAL : You need to select a mode for the HSQ (protonation walk) parameters"
        call fexit()
    end if

end subroutine 


subroutine HS_walk()
    use energies
    use movesets
    use atoms
    use system
    use mcsums
    use iounit
    
    implicit none

    logical done

    integer i 
    RTYPE energy3,edum(MAXENERGYTERMS)
    !print *,"HSQ_step",hsq_steps
    do i=1,hsq_steps
!        print *,"in loop",i
        !flush(ilog)
        ! Here i will try at each step to move the par_pka one step further. If by the end of the loop, I have not 
        ! been able to reach the end state, then I will go back to the original conformation
        ! and try a different change in protonation state
        ! For a first version, I am going to try and simply change everything at once.
        ! In later version, i could try high dimensional HS, but I do not have the energy now  
        ! From there if the HS move is accepted, we go into the HS walk    
        call select_mcmove_tree()    
        
        !Proceeds to the move and evaluation of the energy and acceptance ratio
        call mcmove(nstep,0)
    end do 
    nstep=nstep-hsq_steps ! we should not be counting those steps
    
end subroutine 
! I think the optimization of the HS procedure I should put the 

subroutine switch_bonds() ! Martin : this routine was adapted from makexyz_formol
!
  use iounit
  use sequen
  use polypep
  use zmatrix
  use molecule
  
  
  
  use atoms
!
  implicit none
!
  integer i,i2,i3,i4,chiral,imol
  RTYPE bl,ba,baodi
 ! Martin : so te goal here is to get the alorithm tpo change the bond length dynamically when the hamiltonian is changed
! 
    do imol=1,nmol
      do i=atmol(imol,1)+3,atmol(imol,2)
!        if ((blen_save(i).ne.blen(i)).or.(bang_save(i).ne.bang(i))) then 
            i2 = iz(1,i)
            i3 = iz(2,i)
            i4 = iz(3,i)
            
            chiral = iz(4,i)

            bl = blen(i)
            ba = bang(i)
            baodi = ztor(i)

            call regenxyz(i,i2,bl,i3,ba,i4,baodi,chiral)
!        end if 
      end do
    end do 
end

subroutine regenxyz(at,i2,bl,i3,ba,i4,baodi,chirality) ! Martin : this routine was adapted from  genxyz
!
  use atoms
  use iounit
  use math
  use zmatrix
!
  implicit none
!
  integer at,i2,i3,i4,chirality,i
  RTYPE bl,ba,bl2,baodi,dv3(3),dv2(3),stm1,ctm1,dv4(3)
  RTYPE trx,trz,stm2,ctm2,s1,c1,s2,c2,bl3,cp1(3),np1,cp2(3)
! Martin : so in this routine instead of the usual, we are only going to check if the bond parameters have changed, and only change
!  the coordinates if they have
!
 67   format(a,i6,1x,i6,1x,i6,a,i7,a)
 68   format('Fatal. Ill-defined dihedral angle for atom ',i7,' (1-2:',&
 &i7,': ',f9.4,' A, 1-3:',i7,': ',f9.4,' deg., 1-4:',i7,': ',f9.4,&
 &' deg.)')
 69   format('Warning. Ill-defined bond angles for atom ',i7,' (1-2:',&
 &i7,': ',f9.4,' A, 1-3:',i7,': ',f9.4,' deg., 1-3:',i7,': ',f9.4,&
 &' deg.)')
 70   format('Fatal. Colinear reference atoms for atom ',i7,' (1-2:',&
 &i7,': ',f9.4,' A, 1-3:',i7,': ',f9.4,' deg., 1-3:',i7,': ',f9.4,&
 &' deg.)')
! sanity check
  if ((i2.lt.0).OR.(i3.lt.0).OR.(i4.lt.0)) then
    write(ilog,67) 'Fatal. Bad atom index in xyzgen (',&
 & i2,i3,i4,') for atom ',at,'. This is a bug.'
    call fexit()
  end if
!
! first atom in a molecule always at the origin (note that
! for re-building later these are skipped!), second on +z-axis,
! third in (+)xz-plane



    if (chirality.eq.0) then

!
!     with a torsional angle
!
      stm1 = ba/RADIAN
      stm2 = baodi/RADIAN
      ctm1 = cos(stm1)
      ctm2 = cos(stm2)
      stm2 = sin(stm2)
      stm1 = sin(stm1)
      dv2(1) = x(i2) - x(i3)
      dv2(2) = y(i2) - y(i3)
      dv2(3) = z(i2) - z(i3)
      bl2 = sqrt(sum(dv2(:)**2))
      dv2(:) = dv2(:)/bl2

      if (stm1.le.0.0) then  ! at,i2,i3 colinear, torsion is irrelevant
        dv4(:) = -dv2(:)*ctm1
        x(at) = x(i2) + bl*dv4(1)
        y(at) = y(i2) + bl*dv4(2)
        z(at) = z(i2) + bl*dv4(3)
        return
      end if

      dv3(1) = x(i3) - x(i4)
      dv3(2) = y(i3) - y(i4)
      dv3(3) = z(i3) - z(i4)
      bl3 = sqrt(sum(dv3(:)**2))
      dv3(:) = dv3(:)/bl3
      
!     get the cross product of the normed bond vectors (cp1)

      call crossprod3(dv3,dv2,cp1)

!     get the dot product (acos) of the normed bond vectors
      np1 = sum(dv2(:)*dv3(:))

      if (abs(np1).ge.1.0) then

        write(ilog,*) 'Fatal. This indicates an instable simulation,&
 & a bug, or bad input in genxyz(...). Please report if necessary.'
        write(ilog,*) at,i2,bl,i3,ba,i4,baodi

        call fexit()
      end if


      s1 = sqrt(1.0 - np1**2)
      cp1(:) = cp1(:)/s1
!     get the cross product of the first bond vector with cp1 
      call crossprod3(cp1,dv2,cp2)

!     build the displacement vector and generate coordinates
      dv4(:) = cp2(:)*stm1*ctm2 + cp1(:)*stm1*stm2 - dv2(:)*ctm1
      x(at) = x(i2) + bl*dv4(1)
      y(at) = y(i2) + bl*dv4(2)
      z(at) = z(i2) + bl*dv4(3)

! 
    else if ((chirality.eq.1).OR.(chirality.eq.-1)) then

!     with two bond angles
!
      stm1 = ba/RADIAN
      stm2 = baodi/RADIAN
      ctm1 = cos(stm1)
      ctm2 = cos(stm2)
      dv2(1) = x(i3) - x(i2)
      dv2(2) = y(i3) - y(i2)
      dv2(3) = z(i3) - z(i2)
      bl2 = sqrt(sum(dv2(:)**2))
      dv2(:) = dv2(:)/bl2
      dv3(1) = x(i2) - x(i4)
      dv3(2) = y(i2) - y(i4)
      dv3(3) = z(i2) - z(i4)
      bl3 = sqrt(sum(dv3(:)**2))
      dv3(:) = dv3(:)/bl3
      call crossprod3(dv3,dv2,cp1)
      np1 = sum(dv2(:)*dv3(:))

      if (np1.ge.1.0) then
!       i2, i3, and i4 are colinear, which is fatal as orientation around
!       axis becomes ill-defined (in other words, this case HAS to be expressed
!       as a torsional problem)
        write(ilog,70) at,i2,bl,i3,ba,i4,baodi
        call fexit()
      end if
      s2 = 1.0/(1.0 - np1**2)
      
      ! MArtin : this can be infinity, which means that np1 can be 1 (actually it is -1))
      
      c1 = (-ctm2 - np1*ctm1)*s2
      c2 = (ctm1 + np1*ctm2)*s2
      s1 = (c1*ctm2 - c2*ctm1 + 1.0)*s2
      if (s1.gt.0.0) then
        s1 = chirality*sqrt(s1)
      else if (s1.lt.0.0) then
!       the two bond angle setup can potentially be incompatible at planar centers:
!       warn eventually ...
        if (s1.lt.(-sqrt(1.0*(10.0**(-precision(s1)))))) then
          dblba_wrncnt = dblba_wrncnt + 1
          if (dblba_wrncnt.eq.dblba_wrnlmt) then
            write(ilog,69) at,i2,bl,i3,ba,i4,baodi
            write(ilog,*) 'This is warning number #',dblba_wrncnt,' of this type not all of&
 & which may be displayed.'
            if (10.0*dblba_wrnlmt.gt.0.5*HUGE(dblba_wrnlmt)) then
              dblba_wrncnt = 0 ! reset
            else
              dblba_wrnlmt = dblba_wrnlmt*10
            end if
          end if
        end if
        s1 = sqrt(sum((c1*dv3(:)+c2*dv2(:))**2))
        c1 = c1/s1
        c2 = c2/s1
        s1 = 0.0
      end if
      
      dv4(:) = dv3(:)*c1 + dv2(:)*c2 + cp1(:)*s1
      x(at) = x(i2) + bl*dv4(1)
      y(at) = y(i2) + bl*dv4(2)
      z(at) = z(i2) + bl*dv4(3)
      
    end if
end

subroutine HS_q_move_tree_cste_Q()
    
    use atoms
    use energies 
    use movesets
    use sequen
    use accept
    
    implicit none    
!    This routine is the to get residue to change while conserving the total number of ups 
    integer i,tot_change        

    integer tot_sign,index
    RTYPE random
    logical found_one
    
    q_switch(:)=0
    tot_sign=0 ! Here tot sign only represent the residues that are changing in the moves selected
    tot_change=0    
    
    ! Now I have to go through the residues, and for each pos res I invert, 
    ! I need to either invert another same sign one in the opposite direction, 
    ! or a different sign in the same direction
    
    if (read_charge_state.eqv..false.) then 
        print *,"WARNING : the use of q-cannonical without predetermined charge state is not recommended,&
        & and will fail with any sequence containing Histidines. "
        do while (((tot_sign.eq.0).and.(tot_change.ne.0)).eqv..false.)
            i=CEILING(random()*(n_res))
            if ((seq_q_state(i).ne.-1).and.(seq_q_state(i).ne.0).and.(tot_change.lt.MAX_SWITCH*2).and.(q_switch(i).eq.0)) then 
                if (((((hsq_sign(i).eq.-1).and.((seq_q_state(i).eq.2).or.(seq_q_state(i).eq.3)))&
                &.or.((hsq_sign(i).eq.1).and.(seq_q_state(i).eq.1)))).and.((tot_sign.eq.0).or.(tot_sign.eq.1))) then!That is the &
                !&case where the charge is made more negative
                    tot_sign=tot_sign-1
                    tot_change=tot_change+1
                    if (seq_q_state(i).eq.1) then 
                        q_switch(i)=1
                        seq_q_state(i)=2
                    else
                        q_switch(i)=-1
                        seq_q_state(i)=1
                    end if
                else if (((hsq_sign(i).eq.-1).and.(seq_q_state(i).eq.1))& ! That is the case where the charge is made more positive
                &.or.((hsq_sign(i).eq.1).and.((seq_q_state(i).eq.2).or.(seq_q_state(i).eq.3))).and.((tot_sign.eq.0).or.&
                &(tot_sign.eq.-1))) then 
                    tot_sign=tot_sign+1
                    tot_change=tot_change+1
                    if (seq_q_state(i).eq.1) then 
                        q_switch(i)=1
                        seq_q_state(i)=2
                    else
                        q_switch(i)=-1
                        seq_q_state(i)=1
                    end if
                end if 
            end if 
        end do 
    else 
        found_one=.false.
!    New version relying on predetermined charged states
        do while(found_one.eqv..false.)
            index=CEILING(random()*(n_charge_states))
            if (ALL(q_states(index,:).eq.seq_q_state(1:n_res),n_res).eqv..false.) then 
                q_switch(1:n_res)=q_states(index,:)-seq_q_state(1:n_res)
                seq_q_state(1:n_res)=q_states(index,:)
                found_one=.true.
                acc%curr_tries(1)=acc%curr_tries(2)
                acc%curr_tries(2)=index
                exit 
            end if 
        end do 
    end if 
    
end subroutine

subroutine HS_q_move_tree_cste_FOS()
    use atoms
    use energies 
    use movesets
    use sequen
    implicit none    
!    This routine is the to get residue to change while conserving the total number of ups 
    integer i,tot_change        
    integer tot_sign
    RTYPE random
    
    q_switch(:)=0
    tot_sign=0
    tot_change=0
    
    do while (((tot_sign.eq.0).and.(tot_change.ne.0)).eqv..false.)
        do i=1,nseq
            if ((random().lt.selec_proba).and.(seq_q_state(i).ne.-1).and.(seq_q_state(i).ne.0).and.(tot_change.lt.MAX_SWITCH*2)&
            &.and.(q_switch(i).eq.0)) then 
                if ((seq_q_state(i).eq.1).and.((tot_sign.eq.0).or.(tot_sign.eq.-1))) then
                    tot_change=tot_change+1
                    q_switch(i)=1
                    q_switch(conjugate(i))=1
                    seq_q_state(i)=2
                    tot_sign=tot_sign+1
                else if ((seq_q_state(i).eq.2).and.((tot_sign.eq.0).or.(tot_sign.eq.1))) then
                    tot_change=tot_change+1
                    q_switch(i)=-1
                    q_switch(conjugate(i))=-1
                    seq_q_state(i)=1
                    tot_sign=tot_sign-1
                end if 
             else if (seq_q_state(i).eq.0) then
                    q_switch(i)=0
            end if 
        end do 
    end do 

end subroutine

! this routine picks an elementary Metropolis MC move based on the tree structure
! and increments the appropriate counter

subroutine select_mcmove_tree()
  use movesets
  use system
  use mcsums
  ! Martin added
  use energies
!
  implicit none
!
  RTYPE random
!
  particleinsertion_move = .false.
  particledeletion_move = .false.
  particleidentity_move = .false. 
  pivot_move = .false.
  chi_move = .false.
  fyc_move = .false.
  pucker_move = .false.
  SJconrot_move = .false.
  UJconrot_move = .false.
  lct_move = .false.
  rigid_move = .false.
  rot_move = .false.
  trans_move = .false.
  ph_move = .false.
  omega_move = .false.
  nuc_move = .false.
  nucpuck_move = .false.
  clurb_move = .false.
  torcr_move = .false.
  torcr_move_omega = .false.
  nuccr_move = .false.
  other_move = .false.
  hsq_move= .false.
  hs_move= .false.
! Martin in all these sections, checks the frequency of the different types of move as describe in the parameter file   
if ((do_hsq.eqv..true.).and.(HSQ_in_progress.eqv..false.)) then !If w are in an hsq chain, do not try this 
    if (random().lt.hsq_freq) then ! Do an hsq move (note that hsq will be zero if the hsq mode is off), although making a boolean 
    !test may be more efficient
        hsq_move=.true.
        mvcnt%nhsq=mvcnt%nhsq+1
    end if 
end if 

if ((do_hs.eqv..true.).and.(HS_in_progress.eqv..false.)) then !If w are in an hsq chain, do not try this 
    if (random().lt.hs_freq) then ! Do an hs move (note that hsq will be zero if the hsq mode is off), although making a boolean 
        !test may be more efficient
        hs_move=.true.
        mvcnt%nhs=mvcnt%nhs+1
    end if 
end if 
if ((hsq_move.eqv..false.).and.(hs_move.eqv..false.)) then
  if (random().lt.particleflucfreq) then !
    if (ens%flag.eq.5) then
      if (random().lt.0.5) then
        particleinsertion_move = .true.
        mvcnt%ninsert = mvcnt%ninsert + 1
      else
        particledeletion_move = .true.
        mvcnt%ndelete = mvcnt%ndelete + 1
      end if
    else if (ens%flag.eq.6) then
      particleidentity_move = .true.
      mvcnt%nidentitychange = mvcnt%nidentitychange + 1
    end if
  else
    if (random().lt.(rigidfreq)) then
      if (use_coupledrigid.EQV..true.) then
        if (random().lt.clurb_freq) then
          clurb_move = .true.
          mvcnt%nclurb = mvcnt%nclurb + 1
        else
          rigid_move = .true.
          mvcnt%nrb = mvcnt%nrb + 1
        end if
      else if (random().le.rotfreq) then
        rot_move = .true.
        mvcnt%nrot = mvcnt%nrot + 1
      else
        trans_move = .true.
        mvcnt%ntrans = mvcnt%ntrans + 1
      end if
    else
      if (use_lctmoves.EQV..true.) then
        lct_move = .true.
        mvcnt%nlct = mvcnt%nlct + 1
      else if (use_globmoves.EQV..true.) then
        if (random().le.chifreq) then 
          if (random() .lt. phfreq) then
            ph_move = .true.
            mvcnt%nph = mvcnt%nph + 1
          else 
            chi_move = .true.
            mvcnt%nchi = mvcnt%nchi + 1
          end if
        else
          if (random().le.crfreq) then
            if (random().le.angcrfreq) then
              UJconrot_move = .true.
              mvcnt%nujcr = mvcnt%nujcr + 1
            else
              if (random().le.torcrfreq) then
                if (random().le.torcrfreq_omega) then
                  torcr_move_omega = .true.
                  mvcnt%ndocr = mvcnt%ndocr + 1
                else
                  torcr_move = .true.
                  mvcnt%ndjcr = mvcnt%ndjcr + 1
                end if
              else
                SJconrot_move = .true.
                mvcnt%nsjcr = mvcnt%nsjcr + 1
              end if
            end if
          else
            if (random().le.omegafreq) then
              omega_move = .true.
              mvcnt%nomega = mvcnt%nomega + 1
            else
              if (random().le.nucfreq) then
                if(random().le.nuccrfreq) then
                  nuccr_move = .true.
                  mvcnt%nnuccr = mvcnt%nnuccr + 1
                else
                  if (random().le.nucpuckfreq) then
                    nucpuck_move = .true.
                    mvcnt%npucker = mvcnt%npucker + 1
                  else
                    nuc_move = .true.
                    mvcnt%nnuc = mvcnt%nnuc + 1
                  end if
                end if
              else
                if (random().le.puckerfreq) then
                  pucker_move = .true.
                  mvcnt%npucker = mvcnt%npucker + 1
                else
                  if (random().le.otherfreq) then
                    other_move = .true.
                    mvcnt%nother = mvcnt%nother + 1
                  else
                    fyc_move = .true.
                    mvcnt%nfyc = mvcnt%nfyc + 1
                  end if
                end if
              end if
            end if
          end if
        end if
      else
        if (random().le.chifreq) then 
          if (random() .lt. phfreq) then
            ph_move = .true.
            mvcnt%nph = mvcnt%nph + 1
          else 
            chi_move = .true.
            mvcnt%nchi = mvcnt%nchi + 1
          end if
        else
          if (random().le.crfreq) then
            if (random().le.angcrfreq) then
              UJconrot_move = .true.
              mvcnt%nujcr = mvcnt%nujcr + 1
            else
              if (random().le.torcrfreq) then
                if (random().le.torcrfreq_omega) then
                  torcr_move_omega = .true.
                  mvcnt%ndocr = mvcnt%ndocr + 1
                else
                  torcr_move = .true.
                  mvcnt%ndjcr = mvcnt%ndjcr + 1
                end if
              else
                SJconrot_move = .true.
                mvcnt%nsjcr = mvcnt%nsjcr + 1
              end if
            end if
          else
            if (random().le.omegafreq) then
              omega_move = .true.
              mvcnt%nomega = mvcnt%nomega + 1
            else
              if (random().le.nucfreq) then
                if (random().le.nuccrfreq) then
                  nuccr_move = .true.
                  mvcnt%nnuccr = mvcnt%nnuccr + 1
                else
                 if (random().le.nucpuckfreq) then
                    nucpuck_move = .true.
                    mvcnt%npucker = mvcnt%npucker + 1
                  else
                    nuc_move = .true.
                    mvcnt%nnuc = mvcnt%nnuc + 1
                  end if
                end if
              else
                if (random().le.puckerfreq) then
                  pucker_move = .true.
                  mvcnt%npucker = mvcnt%npucker + 1
                else
                  if (random().le.otherfreq) then
                    other_move = .true.
                    mvcnt%nother = mvcnt%nother + 1
                  else
                    pivot_move = .true.
                    mvcnt%nfy = mvcnt%nfy + 1
                  end if
                end if
              end if
            end if
          end if
        end if
      end if 
    end if
  end if
end if 
!
!  write(*,*) pivot_move,fyc_move,omega_move,chi_move,ph_move,pucker_move
!  write(*,*) particleinsertion_move,particledeletion_move,particleidentity_move
!  write(*,*) UJconrot_move,SJconrot_move,torcr_move,torcr_move_omega
!  write(*,*) nuc_move,nucpuck_move,nuccr_move
!  write(*,*) rigid_move,clurb_move,trans_move,rot_move,hsq_move,hs_move
!
end
!
!-----------------------------------------------------------------------------------------
!
! this is the core routine for individual (elementary) MC moves, but on its own
! it mostly just provides the calling construct for the actual routines
! (for the bulk of those, see setconf.f90, energy_wrap.f90, inner_loops_en.f90, energy.f90)
! all MC moves here use the Metropolis algorithm
!
subroutine mcmove(istep,ndump)
!
  use iounit
  use sequen
  use molecule
  use torsn
  use energies
  use system
  use movesets
  use accept
  use mcsums
  use cutoffs
  use fyoc
  use grandensembles
  use math
  use atoms
  use ujglobals
  use polypep
  use zmatrix
  use dssps
  use aminos
  use mpistuff
  use interfaces
  
  ! MArtin : debug only 
  use martin_own
  
  
!
  implicit none
!
  integer i,j,k,kk,istep,rs,rsi,rsf,imol,jmol,which,ndump
  integer aone,azero,atwo,afour,mmol,athree
  logical moveok,sayyes,sayno,yep,chidone,ct,flipset,foundok,dsspyes
  integer cardinality,currentnum,currentnum2,exactcrtries
  integer randomfluctype,random_wt_fluctype,randommoloftype,typetoinsert,typetodelete
  RTYPE diff,evecp(MAXENERGYTERMS),eveca(MAXENERGYTERMS),difft
  RTYPE evecd(MAXENERGYTERMS),energy3
  RTYPE expterm,random,rcheck,t1,t2,fugacity,normal,blas(40)
! CR and pucker variables
  logical findsol,prolineres
  integer nresidue,nDOF,dihed_displace,dof,solucount,solcnt2,nome,omepos(MAXUJDOF)
  RTYPE oldvals(4),oldvals2(12),oldvals3(18),oldvals_nuc(14),puckdofs(7)
  RTYPE dsquare,jac_a,jac_b,dsquare_star,det_L,det_L_star,prob_move,jacv(MAXSOLU*4),jacvp(MAXSOLU*4)
  RTYPE ccmat(6),ccmat_nuc(7),jacmat(MAXCRDOF,MAXCRDOF),dfy(MAXCRDOF),dgv(MAXCRDOF)
  RTYPE lmat(MAXCRDOF,MAXCRDOF),gmat(MAXCRDOF,MAXCRDOF),bt
  RTYPE bias_fwd,bias_bwd,dum,dumfwd,dumbwd,amat(MAXCRDOF,MAXCRDOF)
  RTYPE ujmat(MAXUJDOF,MAXUJDOF),nuccrmat(MAXUJDOF,MAXUJDOF)
  RTYPE cdecomp_upper_star(MAXUJDOF,MAXUJDOF)
  RTYPE pos_def(MAXUJDOF,MAXUJDOF),pos_def_star(MAXUJDOF,MAXUJDOF)
  RTYPE cdecomp_lower(MAXUJDOF,MAXUJDOF),gvec(MAXUJDOF),curpksp(MAXSOLU*4,7,3)
  RTYPE cdecomp_upper(MAXUJDOF,MAXUJDOF),dphi(MAXUJDOF+6)
  RTYPE cdecomp_lower_star(MAXUJDOF,MAXUJDOF),curpks(MAXSOLU*4,7,3),torcrmat(MAXUJDOF,MAXUJDOF),dccmat(6)
  RTYPE gvec_star(MAXUJDOF),dpfmat(MAXUJDOF+6),dpfmat2(MAXUJDOF+6),solmat(MAXSOLU*4,6),solmatp(MAXSOLU*4,6)
!
! some initializations
  sayyes = .true.
  sayno = .false.
  if (align_NC.eq.2) then
    ct = .true.
  else ! see below in individual move-sets for further modification
    ct = .false.
  end if
  aone = 1
  azero = 0
  atwo = 2
  athree = 3
  afour = 4
  evecp(:) = 0.0
  eveca(:) = 0.0
  prolineres = .false.
  prob_move = 1.0
!
  nstep = nstep + 1
  if (particleinsertion_move.EQV..true.) then
!
!   here we will sample the number of particles of a certain allowed type
!   as part of the grand-canonical ensemble
!   insertion first
!
    typetoinsert = random_wt_fluctype()
    imol = randommoloftype(typetoinsert, 2)
    if (imol.le.0) then
      gc_wrncnt(1) = gc_wrncnt(1) + 1
      if (gc_wrncnt(1).eq.gc_wrnlmt(1)) then
        write(ilog,*) 'WARNING: An insertion move was attempted &
 &but there were no remaining reserve non-present molec&
 &ules of type ', typetoinsert, '. The Markov chain is broken at thi&
 &s point, and it is recommended that all results from this step onw&
 &ards are excluded from any analysis. To prevent this problem, incr&
 &ease the number of reserve particles specified in the sequence fil&
 &e.'
        write(ilog,*) 'This was warning #',gc_wrncnt(1),' of this type not all of which may be displayed.'
        if (10.0*gc_wrnlmt(1).gt.0.5*HUGE(gc_wrnlmt(1))) then
          gc_wrncnt(1) = 0
        else
          gc_wrnlmt(1) = gc_wrnlmt(1)*10
        end if
      end if
      call mc_acceptance_bailout(mc_acc_crit,istep)
      return
    end if
    dsspyes = .false.
    if (use_DSSP.EQV..true.) then
      if (pepmol(imol).gt.0) dsspyes = .true.
    end if
    moveok = .false.
!   Acceptance criteria must be computed using the number of particles
!   BEFORE the insertion
    currentnum = cardinality(typesets(typetoinsert,1))
!   Tentatively inserts the particle and calculates the
!   short-range energy for steric screening purposes
!   We still need to add a rigid-body translation to somewhere 
!   compatible with the box boundaries
    call insertmol(imol)
    call CPU_time(t1)
    call rigid_energy_short(imol, eveca, use_cutoffs, azero)
    if ((use_DSSP.EQV..true.).AND.(dsspyes.EQV..true.)) then
      call en_dssp_gl(eveca)
    end if
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!   Removes the particle and screens for steric clashes
    call deletemol(imol) 
    if (use_stericscreen.EQV..true.) then
      difft = sum(eveca)
      if ((difft.gt.2.0*esave).AND.(difft.gt.screenbarrier)) then
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!   Calculates the original short and long-range energies
    call CPU_time(t1)
    call rigid_energy_short(imol, evecp, use_cutoffs, aone)
    if ((use_DSSP.EQV..true.).AND.(dsspyes.EQV..true.)) then
      call en_dssp_gl(evecp)
    end if
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
    call rigid_energy_long(imol, evecp, use_cutoffs, azero)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!   Re-inserts the particle and calculates the long-range
!   energy only; the short-range energy has been done above
    call insertmol(imol)
    call CPU_time(t1)
    call rigid_energy_long(imol, eveca, use_cutoffs, aone)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   get energy difference and determine acceptance
    fugacity = exp(invtemp*chempot(typetoinsert))
    if (gc_mode.eq.1) then
      prob_move = fugacity*ens%insV/(dble(currentnum+1)*thermalvolume(typetoinsert))
    else
      prob_move = fugacity*eqnum(typetoinsert)/dble(currentnum+1)
    end if

    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
    if (moveok.EQV..true.) then
      acc%ninsert = acc%ninsert + 1
      acc%insert(imol) = acc%insert(imol) + 1
      esave = esave + diff
      esterms(:) = esterms(:) + evecd(:)
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
    else
      call deletemol(imol)
      if (use_IMPSOLV) call init_svte(afour)
    end if
!    write (*,*) 'Current number: ',cardinality(typesets(typetoinsert
! &,1))
!
  else if (particledeletion_move.EQV..true.) then
!
!   now deletion
! 
    typetodelete = random_wt_fluctype()
    imol = randommoloftype(typetodelete, 1)
!   Aborts the move if there are no molecules of the chosen type
!   to delete.  Note: One might also consider only choosing from
!   types which have at least one molecule present, but this would
!   make preserving detailed balance more complicated.
    if (imol.le.0) then
      call mc_acceptance_bailout(mc_acc_crit,istep)
      return
    end if
    dsspyes = .false.
    if (use_DSSP.EQV..true.) then
      if (pepmol(imol).gt.0) dsspyes = .true.
    end if
    moveok = .false.
!   Acceptance criteria must be computed using the number of particles
!   BEFORE the deletion
    currentnum = cardinality(typesets(typetodelete,1))
!   Tentatively deletes the particle and calculates the short-
!   range energy.  There is no need to do steric screening.
    call deletemol(imol)
    call CPU_time(t1)
    call rigid_energy_short(imol, eveca, use_cutoffs, azero)
    if ((use_DSSP.EQV..true.).AND.(dsspyes.EQV..true.)) then
      call en_dssp_gl(eveca)
    end if
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!   Re-inserts the particle and calculates the short and long-
!   range energies.
    call insertmol(imol)
    call CPU_time(t1)
    call rigid_energy_short(imol, evecp, use_cutoffs, aone)
    if ((use_DSSP.EQV..true.).AND.(dsspyes.EQV..true.)) then
      call en_dssp_gl(evecp)
    end if
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
    call rigid_energy_long(imol, evecp, use_cutoffs, azero)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!   Re-deletes the particle and calculates the long-range energy
!   only; the short-range energy has been done above
    call deletemol(imol)
    call CPU_time(t1)
    call rigid_energy_long(imol, eveca, use_cutoffs, aone)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
    fugacity = exp(-invtemp*chempot(typetodelete))
    
    if (gc_mode.eq.1) then
      prob_move = dble(currentnum)*fugacity*thermalvolume(typetodelete)/ens%insV
    else
      prob_move = dble(currentnum)*fugacity/eqnum(typetodelete)
    end if
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
    if (moveok.EQV..true.) then
      acc%ndelete = acc%ndelete + 1
      acc%delete(imol) = acc%delete(imol) + 1
      esave = esave + diff
      esterms(:) = esterms(:) + evecd(:)
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
    else
      call insertmol(imol)
      if (use_IMPSOLV) call init_svte(afour)
    end if
!    write (*,*) 'Current number: ',cardinality(typesets(typetodelete
! &,1))
!
  else if (particleidentity_move.EQV..true.) then
!
!   here we will sample the identity of two different molecule types (identity swap)
!   as part of the semi-grand canonical ensemble
!   to avoid redundant moves, types are forced to be different
!
    typetodelete = randomfluctype()
    imol = randommoloftype(typetodelete, 1)
    if ((imol.le.0).OR.(fluclst%nr.le.1)) then
      call mc_acceptance_bailout(mc_acc_crit,istep)
      return
    end if
    if (fluclst%nr.eq.2) then
      if (typetodelete.eq.1) then
        typetoinsert = 2
      else
        typetoinsert = 1
      end if
    else
      typetoinsert = typetodelete
      do while (typetoinsert.eq.typetodelete)
        typetoinsert = randomfluctype()
      end do
    end if
    jmol = randommoloftype(typetoinsert, 2)
    if (jmol.le.0) then
      write(ilog,*) 'FATAL: A move was attempted to change a particl&
 &e of type ', typetodelete, ' to a particle of type ', typetoinsert&
 &, ', but there were no remaining reserve non-present particles of &
 &type ', typetoinsert, '.  This does not ever need to happen; pleas&
 &e rerun the simulation with sufficient reserve particles specified&
 & in the sequence file.'
      call fexit()
    end if
    dsspyes = .false.
    if (use_DSSP.EQV..true.) then
      if (pepmol(imol).gt.0) dsspyes = .true.
      if (pepmol(jmol).gt.0) dsspyes = .true.
    end if
    moveok = .false.
    currentnum = cardinality(typesets(typetodelete,1))
    currentnum2 = cardinality(typesets(typetoinsert,1))
!   Tentatively changes the particle and calculates the
!   short-range energy for steric screening purposes
    call transmute(imol, jmol)
    call CPU_time(t1)
!   Double counting does not happen here because imol always has a
!   zero interaction with jmol
    call rigid_energy_short(imol, eveca, use_cutoffs, azero)
    call rigid_energy_short(jmol, eveca, use_cutoffs, azero)
    if ((use_DSSP.EQV..true.).AND.(dsspyes.EQV..true.)) then
      call en_dssp_gl(eveca)
    end if
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!   Restores the original particle and screens for steric clashes
    call transmute(jmol, imol)
    if (use_stericscreen.EQV..true.) then
      difft = sum(eveca)
      if ((difft.gt.2.0*esave).AND.(difft.gt.screenbarrier)) then
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!   Calculates the original short and long-range energies
    call CPU_time(t1)
    call rigid_energy_short(imol, evecp, use_cutoffs, aone)
    call rigid_energy_short(jmol, evecp, use_cutoffs, aone)
    if ((use_DSSP.EQV..true.).AND.(dsspyes.EQV..true.)) then
      call en_dssp_gl(evecp)
    end if
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
    call rigid_energy_long(imol, evecp, use_cutoffs, azero)
    call rigid_energy_long(jmol, evecp, use_cutoffs, azero)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!   Re-inserts the particle and calculates the long-range
!   energy only; the short-range energy has been done above
    call transmute(imol, jmol)
    call CPU_time(t1)
    call rigid_energy_long(imol, eveca, use_cutoffs, aone)
    call rigid_energy_long(jmol, eveca, use_cutoffs, aone)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
    fugacity = exp(invtemp*(chempot(typetoinsert)-chempot(typetodelete)))
    if (gc_mode.eq.1) then
      prob_move = fugacity*dble(currentnum)*thermalvolume(typetodelete)/(dble(currentnum2+1)*thermalvolume(typetoinsert))
    else
      prob_move = fugacity*dble(currentnum)*eqnum(typetoinsert)/(eqnum(typetodelete)*dble(currentnum2+1))
    end if
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
    if (moveok.EQV..true.) then
      acc%nidentitychange = acc%nidentitychange + 1
      acc%permute(imol) = acc%permute(imol) + 1
      acc%permute(jmol) = acc%permute(jmol) + 1
      esave = esave + diff
      esterms(:) = esterms(:) + evecd(:)
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
    else
      call transmute(jmol, imol)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
    end if
!
!
!
!
  else if ((rigid_move.EQV..true.).OR.(rot_move.EQV..true.).OR.&
 &                             (trans_move.EQV..true.)) then
!   here we will sample only rigid-body degrees of freedom
    moveok = .false.
!   randomly pick a molecule
    rcheck = random()
    call binary_search(rblst%nr,rblst%wt(1:rblst%nr),rcheck,imol)
    imol = rblst%idx(min(rblst%nr,imol+1))
    dsspyes = .false.
    if (use_DSSP.EQV..true.) then
      if (pepmol(imol).gt.0) dsspyes = .true.
    end if
!   ensure that for mono-atomic molecules no rotations are attempted
    if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
      if (rigid_move.EQV..true.) then
        rigid_move = .false.
        trans_move = .true.
!       correction needed
        mvcnt%nrb = mvcnt%nrb - 1
        mvcnt%ntrans = mvcnt%ntrans + 1
      else if (rot_move.EQV..true.) then
        rot_move = .false.
        trans_move = .true.
!       correction needed
        mvcnt%nrot = mvcnt%nrot - 1
        mvcnt%ntrans = mvcnt%ntrans + 1
      end if
    end if
!   we will sample first to be able to pre-screen the new conformation for
!   steric overlap without even bothering to compute the prior energy
!   backup reference coordinates, perturb coord.s and re-assign grid-stuff
!    write(*,*) '0',x(15),x(30)
    if (rigid_move.EQV..true.) then
      call mcmove_rigid(imol,azero)
    else if (rot_move.EQV..true.) then
      call mcmove_rot(imol,azero)
    else if (trans_move.EQV..true.) then
      call mcmove_trans(imol,azero)
    end if ! Martin : Since this is to get rid of the steric overlap
!   compute posterior short-range (steric) energy
    call CPU_time(t1)
!   first let's compute the posterior short-range contributions
    call rigid_energy_short(imol,eveca,use_cutoffs,azero)
!
!   strictly global energy terms are handled outside of the residue-based functions
    if (use_DREST.EQV..true.) call edrest(eveca)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
    if ((use_DSSP.EQV..true.).AND.(dsspyes.EQV..true.)) then
      call en_dssp_gl(eveca)
    end if
!   now make use of pre-screening (if so desired)
    if (use_stericscreen.EQV..true.) then
      difft = 0.0
      do i=1,MAXENERGYTERMS
        difft = difft + eveca(i)
      end do
!
      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then ! If that happens, conformation is rejected, bailout is called
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
        if (rigid_move.EQV..true.) then
          call mcmove_rigid(imol,atwo)
        else if (rot_move.EQV..true.) then
          call mcmove_rot(imol,atwo)
        else if (trans_move.EQV..true.) then
          call mcmove_trans(imol,atwo)
        end if
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
    

    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   restore the coordinates to the prior state 
    if (rigid_move.EQV..true.) then
      call mcmove_rigid(imol,atwo)
    else if (rot_move.EQV..true.) then
      call mcmove_rot(imol,atwo)
    else if (trans_move.EQV..true.) then
      call mcmove_trans(imol,atwo)
    end if
!    write(*,*) '2',x(15),x(30)
!   now compute the prior energies, the short-range part will already
!   yield the full SAV-information (if relevant)
    call CPU_time(t1)
    call rigid_energy_short(imol,evecp,use_cutoffs,aone)
!   strictly global energy terms are handled outside of the residue-based functions
    if (use_DREST.EQV..true.) call edrest(evecp)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
    if ((use_DSSP.EQV..true.).AND.(dsspyes.EQV..true.)) then
      call en_dssp_gl(evecp)
    end if
    
!   the long-range part has a different cutoff and is computed separately
!   also, it will use the old atsav's, but compute only those pair terms that changed
!   at the end the fxn updates the atsav's to the sampled conformation (generated in the next step)
    call rigid_energy_long(imol,evecp,use_cutoffs,azero)

    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   now re-apply the sampling move
    if (rigid_move.EQV..true.) then
      call mcmove_rigid(imol,aone)
    else if (rot_move.EQV..true.) then
      call mcmove_rot(imol,aone)
    else if (trans_move.EQV..true.) then
      call mcmove_trans(imol,aone)
    end if
!    write(*,*) '3',x(15),x(30)
!   time to get the long-range posterior energies (with updated atsav's)
    call CPU_time(t1)
    call rigid_energy_long(imol,eveca,use_cutoffs,aone)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!   get energy difference and determine acceptance
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!   update or restore
    if (moveok.EQV..true.) then
      if (rigid_move.EQV..true.) then
        acc%nrb = acc%nrb + 1
      else if (rot_move.EQV..true.) then
        acc%nrot = acc%nrot + 1
      else if (trans_move.EQV..true.) then
        acc%ntrans = acc%ntrans + 1
      end if
      acc%rigid(imol) = acc%rigid(imol) + 1
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = esterms(i) + evecd(i)
      end do
      if ((rigid_move.EQV..true.).OR.(rot_move.EQV..true.)) then
        call update_gyrten(imol)
      end if
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
    else
!     restore the coordinates to the prior state
      if (rigid_move.EQV..true.) then
        call mcmove_rigid(imol,atwo)
      else if (rot_move.EQV..true.) then
        call mcmove_rot(imol,atwo)
      else if (trans_move.EQV..true.) then
        call mcmove_trans(imol,atwo)
      end if
      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
    end if
!    write(*,*) '4',diff,x(15),x(30),moveok
!    write(*,*) 'F:',rgpcs(1,1,1),rgpcs(1,1,2),rgpcs(1,1,3)
!    write(*,*) 'F:',com(imol,1),com(imol,2),com(imol,3)
  else if ((clurb_move.EQV..true.)) then
! 
!   here we will sample rigid-body degrees of freedom for several molecules simultaneously
!
    moveok = .false.
    clurbf(:) = 0
    mmol = 0
!
!   initially we pick whether to pick a cluster structurally (proximity)
!   or randomly
    if (random().ge.clurb_rdfreq) then  ! structural
!     determine whether to re-determine cluster distribution
      if (clurb_incr.eq.0) then
!        call molclusters2(clurb_clun,clurb_cluszs,clurb_clus)
        clurb_incr = clurb_incr + 1
      else
        clurb_incr = clurb_incr + 1
        if (clurb_incr.eq.clurb_reset) clurb_incr = 0
      end if
!     we need to accumulate the forward and backward biasing probability!!1 WARNING
      if (clurb_clun.eq.0) then
!       do nothing: we bail out into the random assignment
      else if ((clurb_clun.eq.1).AND.(clurb_cluszs(1).eq.nmol).AND.&
 &             (random().lt.clurb_strbail)) then
!       do nothing with proper bail-out frequency (otherwise, whole system is moved)
      else
        j = floor(random()*clurb_clun)+1
        mmol = clurb_cluszs(j)
        do i=1,mmol
          clurbi(i) = clurb_clus(j,i)
          clurbf(clurbi(i)) = 1
        end do
      end if
    end if
!
    if (mmol.eq.0) then ! random
!     randomly pick a cluster-size (note this uniform prob. corresponds neither to real prob. nor
!     to combinatorial prob.)
      mmol = floor(random()*(clurb_maxsz-1)) + 2
!     randomly pick first molecule
      rcheck = random()
      call binary_search(clurblst%nr,clurblst%wt(1:clurblst%nr),rcheck,kk)
      clurbi(1) = clurblst%idx(min(clurblst%nr,kk+1))
      clurbf(clurbi(1)) = 1
      jac_a = 1.0 - clurblst%wt(kk+1)
!     find further molecules -> list is truncated for unpleasant cases
      do i=2,mmol
        foundok = .false.
        k = 0
        do while ((foundok.EQV..false.).AND.((k.le.100).OR.(jac_a.gt.0.1).OR.(i.le.2)))
          k = k + 1
          rcheck = random()
          call binary_search(clurblst%nr,clurblst%wt(1:clurblst%nr),rcheck,kk)
          clurbi(i) = clurblst%idx(min(clurblst%nr,kk+1))
          foundok = .true.
          do j=1,i-1
            if (clurbi(i).eq.clurbi(j)) foundok = .false.
          end do
        end do
        if (k.gt.100) then
          mmol = i -1
          exit
        end if
        clurbf(clurbi(i)) = 1
        jac_a = 1.0 - clurblst%wt(kk+1)
      end do
    end if
    dsspyes = .false.
    if (use_DSSP.EQV..true.) then
      do i=1,mmol
        if (pepmol(clurbi(i)).gt.0) then
          dsspyes = .true.
          exit
        end if
      end do
    end if
!
!   we will sample first to be able to pre-screen the new conformation for
!   steric overlap without even bothering to compute the prior energy
!
!   backup reference coordinates, perturb coord.s and re-assign grid-stuff
!    write(*,*) '0',x(1),x(5)
    call mcmove_clurigid(clurbi,mmol,azero)
!    write(*,*) '1',x(1),x(5)
!
!   compute posterior short-range (steric) energy
    call CPU_time(t1)
!
!   first let's compute the posterior short-range contributions
    call clurb_energy_short(clurbi,mmol,clurbf,eveca,use_cutoffs,&
 &azero)
!   strictly global energy terms are handled outside of the residue-based functions
    if (use_DREST.EQV..true.) call edrest(eveca)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
    if ((use_DSSP.EQV..true.).AND.(dsspyes.EQV..true.)) then
      call en_dssp_gl(eveca)
    end if
!
!   now make use of pre-screening (if so desired)
    if (use_stericscreen.EQV..true.) then
      difft = 0.0
      do i=1,MAXENERGYTERMS
        difft = difft + eveca(i)
      end do
      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
        call mcmove_clurigid(clurbi,mmol,atwo)
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   restore the coordinates to the prior state
    call mcmove_clurigid(clurbi,mmol,atwo)
!    write(*,*) '2',x(1),x(5)
!
!   now compute the prior energies, the short-range part will already
!   yield the full SAV-information (if relevant)
    call CPU_time(t1)
    call clurb_energy_short(clurbi,mmol,clurbf,evecp,use_cutoffs,&
 &aone)
!   strictly global energy terms are handled outside of the residue-based functions
    if (use_DREST.EQV..true.) call edrest(evecp)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
    if ((use_DSSP.EQV..true.).AND.(dsspyes.EQV..true.)) then
      call en_dssp_gl(evecp)
    end if
!    write(*,*) evecp(1)
!
!   the long-range part has a different cutoff and is computed separately
!   also, it will use the old atsav's, but compute only those pair terms that changed
!   at the end the fxn updates the atsav's to the sampled conformation (generated in the next step)
    call clurb_energy_long(clurbi,mmol,clurbf,evecp,use_cutoffs,&
 &azero)

    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   now re-apply the sampling move
    call mcmove_clurigid(clurbi,mmol,aone)
!    write(*,*) '3',x(1),x(5)
!
!   time to get the long-range posterior energies (with updated atsav's)
    call CPU_time(t1)
    call clurb_energy_long(clurbi,mmol,clurbf,eveca,use_cutoffs,&
 &aone)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   get energy difference and determine acceptance
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
!   update or restore
    if (moveok.EQV..true.) then
      acc%nclurb = acc%nclurb + 1
      do i=1,mmol
        acc%rigid(clurbi(i)) = acc%rigid(clurbi(i)) + 1
      end do
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = esterms(i) + evecd(i)
      end do
      if ((rigid_move.EQV..true.).OR.(rot_move.EQV..true.)) then
        do i=1,mmol
          call update_gyrten(clurbi(i))
        end do
      end if
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
    else
!     restore the coordinates to the prior state
      call mcmove_clurigid(clurbi,mmol,atwo)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour) ! Martin : That is an important step 
    end if
!    write(*,*) '4',x(1),x(5),diff,moveok
!    write(*,*) 'F:',rgpcs(1,1,1),rgpcs(1,1,2),rgpcs(1,1,3)
!    write(*,*) 'F:',com(imol,1),com(imol,2),com(imol,3)
!
!
!
!
  else if (chi_move.EQV..true.) then
! 
!   here we will sample only internal sidechain degrees of freedom
!
    moveok = .false.
!
!   randomly pick a sidechain
    rcheck = random()
    call binary_search(chilst%nr,chilst%wt(1:chilst%nr),rcheck,rs)
    rs = chilst%idx(min(chilst%nr,rs+1))
!   now randomly pick the appropriate number of chi-angles to-be-sampled
    chidone = .false.
    
    
    do while (chidone.EQV..false.)
      if (nchi(rs).le.nrchis) then
        do i=1,nchi(rs)
          cur_chiflag(i) = .true.
        end do
        chidone = .true.
      else
        do i=1,nchi(rs)
          if (random().lt.((1.0*nrchis)/nchi(rs))) then
            cur_chiflag(i) = .true.
            chidone = .true.
          end if
        end do
      end if
    end do
!
!   we will sample first to be able to pre-screen the new conformation for
!   steric overlap without even bothering to compute the prior energy
!
!   backup reference coordinates, perturb coord.s and re-assign grid-stuff
!    write(*,*) '0',chi(1,2)
    call mcmove_chi(rs,azero)
!    write(*,*) '1',chi(1,2)
!
!   compute posterior short-range (steric) energy
    call CPU_time(t1)
!
!   first let's compute the posterior short-range contributions
    call chi_energy_short(rs,eveca,use_cutoffs,azero)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),eveca,azero)
    if (use_DREST.EQV..true.) call edrest(eveca)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
!   now make use of pre-screening (if so desired)
    if (use_stericscreen.EQV..true.) then
      difft = 0.0
      do i=1,MAXENERGYTERMS
        difft = difft + eveca(i)
      end do
      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
        call mcmove_chi(rs,atwo)
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   restore the coordinates to the prior state
    call mcmove_chi(rs,atwo)
!    write(*,*) '2',chi(1,2)
!
!   now compute the prior energies, the short-range part will already
!   yield the full SAV-information (if relevant)
    call CPU_time(t1)
    call chi_energy_short(rs,evecp,use_cutoffs,aone)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),evecp,aone)
    if (use_DREST.EQV..true.) call edrest(evecp)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
!   the long-range part has a different cutoff and is computed separately
!   also, it will use the old atsav's, but compute only those pair terms that changed
!   at the end the fxn updates the atsav's to the sampled conformation (generated in the next step)
    call chi_energy_long(rs,evecp,use_cutoffs,azero)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   now re-apply the sampling move
    call mcmove_chi(rs,aone)
!    write(*,*) '3',chi(1,2)
!
!   time to get the long-range posterior energies (with updated atsav's)
    call CPU_time(t1)
    call chi_energy_long(rs,eveca,use_cutoffs,aone)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1            
!
!   get energy difference and determine acceptance
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
!   update or restore
    if (moveok.EQV..true.) then
      acc%nchi = acc%nchi + 1
      acc%chi(rs) = acc%chi(rs) + 1
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = esterms(i) + evecd(i)
      end do
      call internal_bound(molofrs(rs))
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
    else
!     restore the coordinates to the prior state
      call mcmove_chi(rs,atwo)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
    end if
!
!   finally, restore the active-chi array 
    do i=1,nchi(rs)
      cur_chiflag(i) = .false.
    end do
!    write(*,*) '4',chi(1,1),chi(2,1),moveok
!
!
!
!
!
  else if (pivot_move.EQV..true.) then
!

!   here we will sample only a single set of backbone angles (pivot move)
!
    moveok = .false.
!
!   randomly pick a residue
    rcheck = random()
    call binary_search(fylst%nr,fylst%wt(1:fylst%nr),rcheck,rs)
    rs = fylst%idx(min(fylst%nr,rs+1))
    if ((rsmol(molofrs(rs),1).eq.rsmol(molofrs(rs),2))) then
      ct = .false.
    else if (align_NC.eq.3) then
      if ((rsmol(molofrs(rs),2)-rs).gt.(rs-rsmol(molofrs(rs),1)))&
 & then
        ct = .true.
      else if (((rsmol(molofrs(rs),2)-rs)&
 &.eq.(rs-rsmol(molofrs(rs),1))).AND.(random().gt.0.5)) then
        ct = .true.
      end if
    else if (align_NC.eq.4) then
      bt = 1.0*(rsmol(molofrs(rs),2)-rs+1)/&
 &         (rsmol(molofrs(rs),2)-rsmol(molofrs(rs),1)+2)
      if (random().le.bt) ct = .true.
    end if    
!   backup reference coordinates, perturb coord.s and re-assign grid-stuff
!    write(*,*) '0',ztor(yline(53)),psi(53),z(822)
    call mcmove_pivot(rs,azero,ct)
!    write(*,*) '1',ztor(yline(4)),psi(4)
!   compute posterior short-range (steric) energy
    call CPU_time(t1)
!
!   first let's compute the posterior short-range contributions

    call pivot_energy_short(rs,eveca,use_cutoffs,azero,ct)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rs),eveca)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),eveca,azero)
    if (use_DSSP.EQV..true.) call en_dssp_gl(eveca)
    if (use_DREST.EQV..true.) call edrest(eveca)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
!   now make use of pre-screening (if so desired)
    if (use_stericscreen.EQV..true.) then
      difft = 0.0
      do i=1,MAXENERGYTERMS
        difft = difft + eveca(i)
      end do
      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
        call mcmove_pivot(rs,atwo,ct)
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
    
!   restore the coordinates to the prior state
    call mcmove_pivot(rs,atwo,ct)
!    write(*,*) '2',ztor(yline(4)),psi(4)

!   now compute the prior energies, the short-range part will already
!   yield the full SAV-information (if relevant)
    call CPU_time(t1)

    call pivot_energy_short(rs,evecp,use_cutoffs,aone,ct)
   
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rs),evecp)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),evecp,aone)
    if (use_DSSP.EQV..true.) call en_dssp_gl(evecp)
    if (use_DREST.EQV..true.) call edrest(evecp)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
!   the long-range part has a different cutoff and is computed separately
!   also, it will use the old atsav's, but compute only those pair terms that changed
!   at the end the fxn updates the atsav's to the sampled conformation (generated in the next step)
    call pivot_energy_long(rs,evecp,use_cutoffs,azero,ct)
    
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   now re-apply the sampling move
    call mcmove_pivot(rs,aone,ct)
!    write(*,*) '3',ztor(yline(4)),psi(4)
!
!   time to get the long-range posterior energies (with updated atsav's)
    call CPU_time(t1)
    call pivot_energy_long(rs,eveca,use_cutoffs,aone,ct)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1            
!   get energy difference and determine acceptance
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
!   update or restore
    if (moveok.EQV..true.) then
      acc%nfy = acc%nfy + 1
      acc%fy(rs) = acc%fy(rs) + 1
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = esterms(i) + evecd(i)
      end do
      if (use_POLY.EQV..false.) call internal_bound(molofrs(rs))
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
    else
      call mcmove_pivot(rs,atwo,ct)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
    end if
!
!
!
!
!
  else if (nuc_move.EQV..true.) then
!
!   here we will sample a (sub)set of nucleic acid (or generic) backbone angles on a single
!   residue
!
    moveok = .false.
!
!   randomly pick a residue
    rcheck = random()
    call binary_search(nuclst%nr,nuclst%wt(1:nuclst%nr),rcheck,rs)
    rs = nuclst%idx(min(nuclst%nr,rs+1))
!
!   randomly pick the individual degrees of freedom
    chidone = .false.
    do while (chidone.EQV..false.)
      if (nnucs(rs).le.nrnucim) then
        do i=1,nnucs(rs)
          cur_nucflag(i) = .true.
        end do
        chidone = .true.
      else
        do i=1,nnucs(rs)
          if (random().lt.((1.0*nrnucim)/nnucs(rs))) then
            cur_nucflag(i) = .true.
            chidone = .true.
          end if
        end do
      end if
    end do
!
!   handle alignment
    if (align_NC.eq.3) then
      if ((rsmol(molofrs(rs),2)-rs).gt.(rs-rsmol(molofrs(rs),1)))&
 & then
        ct = .true.
      else if (((rsmol(molofrs(rs),2)-rs)&
 &.eq.(rs-rsmol(molofrs(rs),1))).AND.(random().gt.0.5)) then
        ct = .true.
      end if
    else if (align_NC.eq.4) then
      bt = 1.0*(rsmol(molofrs(rs),2)-rs+1)/&
 &         (rsmol(molofrs(rs),2)-rsmol(molofrs(rs),1)+2)
      if (random().le.bt) ct = .true.
    end if
!
!   backup reference coordinates, perturb coord.s and re-assign grid-stuff
!    write(*,*) '0',ztor(nucsline(4,2)),nucs(4,2)
    call mcmove_nuc(rs,azero,ct)
!    write(*,*) '1',ztor(nucsline(4,2)),nucs(4,2)
!
!   compute posterior short-range (steric) energy
    call CPU_time(t1)
!
!   first let's compute the posterior short-range contributions
    call pivot_energy_short(rs,eveca,use_cutoffs,azero,ct)
!   strictly global energy fxns are handled outside of the residue-based fxns
!    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rs),eveca)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),eveca,azero)
    if (use_DREST.EQV..true.) call edrest(eveca)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
!   now make use of pre-screening (if so desired)
    if (use_stericscreen.EQV..true.) then
      difft = 0.0
      do i=1,MAXENERGYTERMS
        difft = difft + eveca(i)
      end do
      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
        call mcmove_nuc(rs,atwo,ct)
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   restore the coordinates to the prior state
    call mcmove_nuc(rs,atwo,ct)
!    write(*,*) '2',ztor(nucsline(4,2)),nucs(4,2)
!
!   now compute the prior energies, the short-range part will already
!   yield the full SAV-information (if relevant)
    call CPU_time(t1)
    call pivot_energy_short(rs,evecp,use_cutoffs,aone,ct)
!   strictly global energy fxns are handled outside of the residue-based fxns
!    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rs),evecp)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),evecp,aone)
    if (use_DREST.EQV..true.) call edrest(evecp)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
!   the long-range part has a different cutoff and is computed separately
!   also, it will use the old atsav's, but compute only those pair terms that changed
!   at the end the fxn updates the atsav's to the sampled conformation (generated in the next step)
    call pivot_energy_long(rs,evecp,use_cutoffs,azero,ct)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   now re-apply the sampling move
    call mcmove_nuc(rs,aone,ct)
!    write(*,*) '3',ztor(nucsline(4,2)),nucs(4,2)
!
!   time to get the long-range posterior energies (with updated atsav's)
    call CPU_time(t1)
    call pivot_energy_long(rs,eveca,use_cutoffs,aone,ct)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   get energy difference and determine acceptance
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
!   update or restore
    if (moveok.EQV..true.) then
      acc%nnuc = acc%nnuc + 1
      acc%nuc(rs) = acc%nuc(rs) + 1
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = esterms(i) + evecd(i)
      end do
      if (use_POLY.EQV..false.) call internal_bound(molofrs(rs))
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
!      
    else
      call mcmove_nuc(rs,atwo,ct)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
    end if
!    write(*,*) 'F',ztor(nucsline(4,2)),nucs(4,2),moveok
!
!
!
!
!
  else if (nucpuck_move.EQV..true.) then
!
!   here we will sample a single sugar pucker "angle" via two methods
!   ensuring proper closure, maintenance, and sampling of reasonable 5-ring geometries
!
    moveok = .false.
!
!   randomly pick a residue
    rcheck = random()
    call binary_search(nucpuclst%nr,nucpuclst%wt(1:nucpuclst%nr),rcheck,rs)
    rs = nucpuclst%idx(min(nucpuclst%nr,rs+1))
!
!   handle alignment
    if (align_NC.eq.3) then
      if ((rsmol(molofrs(rs),2)-rs).gt.(rs-rsmol(molofrs(rs),1)))&
 & then
        ct = .true.
      else if (((rsmol(molofrs(rs),2)-rs)&
 &.eq.(rs-rsmol(molofrs(rs),1))).AND.(random().gt.0.5)) then
        ct = .true.
      end if
    else if (align_NC.eq.4) then
      bt = 1.0*(rsmol(molofrs(rs),2)-rs+1)/&
 &         (rsmol(molofrs(rs),2)-rsmol(molofrs(rs),1)+2)
      if (random().le.bt) ct = .true.
    end if
!
!   get reference jacobian
    call jacobian_pucker(rs,jac_a)
!
!   backup reference coordinates, perturb driver coord.s, close ring and re-assign grid-stuff
!    write(*,*) '0',ztor(nucsline(4,2)),nucs(4,2)
    call mcmove_sugar(rs,puckdofs,azero,ct)
!    write(*,*) '1',ztor(nucsline(4,2)),nucs(4,2)
!
!   compute posterior short-range (steric) energy
    call CPU_time(t1)
!
!   first let's compute the posterior short-range contributions
    call pivot_energy_short(rs,eveca,use_cutoffs,azero,ct)
!   strictly global energy fxns are handled outside of the residue-based fxns
!    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rs),eveca)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),eveca,azero)
    if (use_DREST.EQV..true.) call edrest(eveca)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
!   now make use of pre-screening (if so desired)
    if (use_stericscreen.EQV..true.) then
      difft = 0.0
      do i=1,MAXENERGYTERMS
        difft = difft + eveca(i)
      end do
      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
        call mcmove_sugar(rs,puckdofs,atwo,ct)
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   get posterior jacobian
    call jacobian_pucker(rs,jac_b)
!
!   restore the coordinates to the prior state
    call mcmove_sugar(rs,puckdofs,atwo,ct)
!    write(*,*) '2',ztor(nucsline(4,2)),nucs(4,2)
!
!   now compute the prior energies, the short-range part will already
!   yield the full SAV-information (if relevant)
    call CPU_time(t1)
    call pivot_energy_short(rs,evecp,use_cutoffs,aone,ct)
!   strictly global energy fxns are handled outside of the residue-based fxns
!    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rs),evecp)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),evecp,aone)
    if (use_DREST.EQV..true.) call edrest(evecp)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
!   the long-range part has a different cutoff and is computed separately
!   also, it will use the old atsav's, but compute only those pair terms that changed
!   at the end the fxn updates the atsav's to the sampled conformation (generated in the next step)
    call pivot_energy_long(rs,evecp,use_cutoffs,azero,ct)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   now re-apply the sampling move
    call mcmove_sugar(rs,puckdofs,aone,ct)
!    write(*,*) '3',ztor(nucsline(4,2)),nucs(4,2)
!
!   time to get the long-range posterior energies (with updated atsav's)
    call CPU_time(t1)
    call pivot_energy_long(rs,eveca,use_cutoffs,aone,ct)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   get energy difference and determine acceptance
    prob_move = jac_b/jac_a
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
!   update or restore
    if (moveok.EQV..true.) then
      acc%npucker = acc%npucker + 1
      acc%pucker(rs) = acc%pucker(rs) + 1
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = esterms(i) + evecd(i)
      end do
      if (use_POLY.EQV..false.) call internal_bound(molofrs(rs))
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
!      
    else
      call mcmove_sugar(rs,puckdofs,atwo,ct)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
    end if
!    if (moveok) write(*,*) '+',getztor(nuci(2,4),nuci(2,5),nuci(2,6),nuci(3,1))-ztorpr(nuci(3,1))
!    if (moveok) write(*,*) '++',getztor(nuci(2,3),nuci(2,4),nuci(2,5),nuci(2,6))-ztorpr(nuci(2,6))
!    write(*,*) 'F',ztor(nucsline(4,2)),nucs(4,2),moveok
!
!
!
!
!
  else if (pucker_move.EQV..true.) then
!
!   here we will sample a single sugar pucker "angle" via two algorithms
!   ensuring proper closure, maintenance, and sampling of reasonable 5-ring geometries
!
    moveok = .false.
!
!   randomly pick a proline (or evtl. HYP, PCA) residue
    rcheck = random()
    call binary_search(puclst%nr,puclst%wt(1:puclst%nr),rcheck,rs)
    rs = puclst%idx(min(puclst%nr,rs+1))
!
!   handle alignment
    if ((rsmol(molofrs(rs),1).eq.rsmol(molofrs(rs),2))) then
      ct = .false.
    else if (align_NC.eq.3) then
      if ((rsmol(molofrs(rs),2)-rs).gt.(rs-rsmol(molofrs(rs),1)))&
 & then
        ct = .true.
      else if (((rsmol(molofrs(rs),2)-rs)&
 &.eq.(rs-rsmol(molofrs(rs),1))).AND.(random().gt.0.5)) then
        ct = .true.
      end if
    else if (align_NC.eq.4) then
      bt = 1.0*(rsmol(molofrs(rs),2)-rs+1)/&
 &         (rsmol(molofrs(rs),2)-rsmol(molofrs(rs),1)+2)
      if (random().le.bt) ct = .true.
    end if
!
!   get reference jacobian
    call jacobian_pucker(rs,jac_a)
!
!   backup reference coordinates, perturb driver coord.s, close ring and re-assign grid-stuff
    call mcmove_pucker(rs,puckdofs,azero,ct)
!
!   compute posterior short-range (steric) energy
    call CPU_time(t1)
!
!   first let's compute the posterior short-range contributions
    blas(17) = phi(rs)
    phi(rs) = ztor(fline(rs))
    call pivot_energy_short(rs,eveca,use_cutoffs,azero,ct)
    phi(rs) = blas(17)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),eveca,azero)
!   make sure pointer array is right for posterior energy
    if (use_ZSEC.EQV..true.) then
      blas(17) = phi(rs)
      phi(rs) = ztor(fline(rs))
      call en_zsec_gl(molofrs(rs),eveca)
      phi(rs) = blas(17)
    end if
    if (use_DSSP.EQV..true.) call en_dssp_gl(eveca)
    if (use_DREST.EQV..true.) call edrest(eveca)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
!   now make use of pre-screening (if so desired)
    if (use_stericscreen.EQV..true.) then
      difft = 0.0
      do i=1,MAXENERGYTERMS
        difft = difft + eveca(i)
      end do
      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
        call mcmove_pucker(rs,puckdofs,atwo,ct)
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   get posterior jacobian
    call jacobian_pucker(rs,jac_b)
!
!   restore the coordinates to the prior state
    call mcmove_pucker(rs,puckdofs,atwo,ct)
!
!   now compute the prior energies, the short-range part will already
!   yield the full SAV-information (if relevant)
    call CPU_time(t1)
    call pivot_energy_short(rs,evecp,use_cutoffs,aone,ct)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),evecp,aone)
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rs),evecp)
    if (use_DSSP.EQV..true.) call en_dssp_gl(evecp)
    if (use_DREST.EQV..true.) call edrest(evecp)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
!   the long-range part has a different cutoff and is computed separately
!   also, it will use the old atsav's, but compute only those pair terms that changed
!   at the end the fxn updates the atsav's to the sampled conformation (generated in the next step)
    call pivot_energy_long(rs,evecp,use_cutoffs,azero,ct)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   now re-apply the sampling move
    call mcmove_pucker(rs,puckdofs,aone,ct)
!
!   time to get the long-range posterior energies (with updated atsav's)
    call CPU_time(t1)
    blas(17) = phi(rs)
    phi(rs) = ztor(fline(rs))
    call pivot_energy_long(rs,eveca,use_cutoffs,aone,ct)
    phi(rs) = blas(17)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   get energy difference and determine acceptance
    prob_move = jac_b/jac_a
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
!   update or restore
    if (moveok.EQV..true.) then
      acc%npucker = acc%npucker + 1
      acc%pucker(rs) = acc%pucker(rs) + 1
!     update pointer array (important!)
      phi(rs) = ztor(fline(rs))
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = esterms(i) + evecd(i)
      end do
      if (use_POLY.EQV..false.) call internal_bound(molofrs(rs))
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
!      
    else
      call mcmove_pucker(rs,puckdofs,atwo,ct)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
    end if
!
!
!
!
!
  else if (omega_move.EQV..true.) then
!
!   here we will sample only a single omega angle
!
    moveok = .false.
!
!   randomly pick a residue
    rcheck = random()
    call binary_search(wlst%nr,wlst%wt(1:wlst%nr),rcheck,rs)
    rs = wlst%idx(min(wlst%nr,rs+1))
    if (align_NC.eq.3) then
      if ((rsmol(molofrs(rs),2)-rs).gt.(rs-rsmol(molofrs(rs),1)))&
 & then
        ct = .true.
      else if (((rsmol(molofrs(rs),2)-rs)&
 &.eq.(rs-rsmol(molofrs(rs),1))).AND.(random().gt.0.5)) then
        ct = .true.
      end if
    else if (align_NC.eq.4) then
      bt = 1.0*(rsmol(molofrs(rs),2)-rs+1)/&
 &         (rsmol(molofrs(rs),2)-rsmol(molofrs(rs),1)+2)
      if (random().le.bt) ct = .true.
    end if
!
!   backup reference coordinates, perturb coord.s and re-assign grid-stuff
!    write(*,*) '0',ztor(yline(4)),psi(4)
    call mcmove_omega(rs,azero,ct)
!    write(*,*) '1',ztor(yline(4)),psi(4)
!
!   compute posterior short-range (steric) energy
    call CPU_time(t1)
!
!   first let's compute the posterior short-range contributions
    call pivot_energy_short(rs,eveca,use_cutoffs,azero,ct)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rs),eveca)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),eveca,azero)
    if (use_DSSP.EQV..true.) call en_dssp_gl(eveca)
    if (use_DREST.EQV..true.) call edrest(eveca)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
!   now make use of pre-screening (if so desired)
    if (use_stericscreen.EQV..true.) then
      difft = 0.0
      do i=1,MAXENERGYTERMS
        difft = difft + eveca(i)
      end do
      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
        call mcmove_omega(rs,atwo,ct)
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   restore the coordinates to the prior state
    call mcmove_omega(rs,atwo,ct)
!    write(*,*) '2',ztor(wline(3)),omega(3)
!
!   now compute the prior energies, the short-range part will already
!   yield the full SAV-information (if relevant)
    call CPU_time(t1)
    call pivot_energy_short(rs,evecp,use_cutoffs,aone,ct)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rs),evecp)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),evecp,aone)
    if (use_DSSP.EQV..true.) call en_dssp_gl(evecp)
    if (use_DREST.EQV..true.) call edrest(evecp)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
!   the long-range part has a different cutoff and is computed separately
!   also, it will use the old atsav's, but compute only those pair terms that changed
!   at the end the fxn updates the atsav's to the sampled conformation (generated in the next step)
    call pivot_energy_long(rs,evecp,use_cutoffs,azero,ct)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   now re-apply the sampling move
    call mcmove_omega(rs,aone,ct)
!    write(*,*) '3',ztor(wline(3)),omega(3)
!
!   time to get the long-range posterior energies (with updated atsav's)
    call CPU_time(t1)
    call pivot_energy_long(rs,eveca,use_cutoffs,aone,ct)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1            
!
!   get energy difference and determine acceptance
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
!   update or restore
    if (moveok.EQV..true.) then
      acc%nomega = acc%nomega + 1
      acc%omega(rs) = acc%omega(rs) + 1
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = esterms(i) + evecd(i)
      end do
      call internal_bound(molofrs(rs))
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
    else
      call mcmove_omega(rs,atwo,ct)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour) ! Martin : So theree is the restore 
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
    end if
!    write(*,*) moveok,'4',ztor(wline(3)),omega(3)
!
!
!
!
!
  else if (other_move.EQV..true.) then
!
!   here we will sample a single dihedral angle 
!
    moveok = .false.
!
!   randomly pick a d.o.f
    yep = .false.
    if (random().le.other_unkfreq) then
      rcheck = random()
      call binary_search(unklst%nr,unklst%wt(1:unklst%nr),rcheck,which)
      which = unklst%idx(min(unklst%nr,which+1))
    else
      if (random().le.other_natfreq) then
        rcheck = random()
        call binary_search(natlst%nr,natlst%wt(1:natlst%nr),rcheck,which)
        which = natlst%idx(min(natlst%nr,which+1))
        yep = .true.
      else
        rcheck = random()
        call binary_search(unslst%nr,unslst%wt(1:unslst%nr),rcheck,which)
        which = unslst%idx(min(unslst%nr,which+1))
      end if
    end if
    rs = atmres(which)
    imol = molofrs(rs) ! robust
!   the rotation set may already be flipped, detect this with the help of iz(3,which)
    flipset = .false.
    rsi = HUGE(rsi)
    do k=1,izrot(which)%alsz
      j = izrot(which)%rotis(k,1)
      kk = izrot(which)%rotis(k,2)
      if ((iz(3,which).ge.j).AND.(iz(3,which).le.kk)) then
        flipset = .true.
        exit
      end if
    end do
    if (align_NC.ne.1) then
      nome = 0
      dum = 0.0
      if (izrot(which)%alsz.gt.0) currentnum = minval(atmres(izrot(which)%rotis(1,1:2)))
      if (rsmol(imol,2).ne.rsmol(imol,1)) ct = .false.
      do i=1,izrot(which)%alsz
        nome = nome + izrot(which)%rotis(i,2) - izrot(which)%rotis(i,1) + 1
        dum = dum + sum(mass(izrot(which)%rotis(i,1):izrot(which)%rotis(i,2)))
      end do
      if (align_NC.eq.2) then
        if (flipset.EQV..false.) then ! we have to check (default rotation set)
          do i=1,izrot(which)%alsz
            if ((atmres(izrot(which)%rotis(i,2)).ne.atmres(izrot(which)%rotis(i,1))).OR.&
 &            (maxval(atmres(izrot(which)%rotis(i,1:2))).ne.currentnum)) then
              ct = .true. ! all unflipped rotation lists spanning multiple residues are eligible
            end if
          end do
          if (izrot(which)%alsz.gt.0) then
            if ((izrot(which)%treevs(2).eq.1).AND.(izrot(which)%treevs(5).gt.0)) ct = .true.
          end if
        else ! the set has already been flipped which indicates C-alignment in hybrid run
          ct = .true.
        end if
      end if
      if (align_NC.eq.2) then
        if (dum.le.11.5) ct = .false.
      else if (align_NC.eq.3) then
        if (flipset.EQV..true.) then ! C-terminal rotation set
          if (nome.le.(atmol(imol,2)-atmol(imol,1)-1-nome)) ct = .true.
        else
          if (nome.gt.(atmol(imol,2)-atmol(imol,1)-1-nome)) ct = .true.
        end if
      else if (align_NC.eq.4) then
        if (flipset.EQV..true.) then
          bt = (1.0*(atmol(imol,2)-atmol(imol,1)-1-nome))/(1.0*(atmol(imol,2)-atmol(imol,1)-1))
        else
          bt = (1.0*nome)/(1.0*(atmol(imol,2)-atmol(imol,1)-1))
        end if
        if (random().le.bt) ct = .true.
      end if
    end if
    if (ct.EQV..true.) then
      if (flipset.EQV..true.) then ! flipset true means a C-terminal rot list
        flipset = .false.
      else
        flipset = .true.
      end if
    end if
!
!   backup reference coordinates, perturb coord.s and re-assign grid-stuff
    call mcmove_other(which,azero,ct,flipset,rsi,rsf,yep)
!   double check and adjustment needed to make sure that effective rotation set matches rs
    if (rs.gt.rsf) rs = rsf
    if (rs.lt.rsi) rs = rsi
!   this condition most likely means that a sidechain angle is C-terminally aligned, which is
!   common in mode 4 -> ignore since it works with other_energy_*
!    if ((rs.ne.rsf).AND.(rs.ne.rsi)) write(*,*) 'WEIRD',rs,rsi,rsf
!
!   compute posterior short-range (steric) energy
    call CPU_time(t1)
!
!   first let's compute the posterior short-range contributions
    call other_energy_short(rs,rsi,rsf,eveca,use_cutoffs,azero)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if ((use_ZSEC.EQV..true.).AND.(seqpolty(rs).eq.'P')) call en_zsec_gl(imol,eveca)
    if (use_POLY.EQV..true.) call en_poly_gl(imol,eveca,azero)
    if ((use_DSSP.EQV..true.).AND.(seqpolty(rs).eq.'P')) call en_dssp_gl(eveca)
    if (use_DREST.EQV..true.) call edrest(eveca)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
!   now make use of pre-screening (if so desired)
    if (use_stericscreen.EQV..true.) then
      difft = 0.0
      do i=1,MAXENERGYTERMS
        difft = difft + eveca(i)
      end do
      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
        call mcmove_other(which,atwo,ct,flipset,rsi,rsf,yep)
        if (use_POLY.EQV..true.) call restore_poly(imol)
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   restore the coordinates to the prior state
    call mcmove_other(which,atwo,ct,flipset,rsi,rsf,yep)
!
!   now compute the prior energies, the short-range part will already
!   yield the full SAV-information (if relevant)
    call CPU_time(t1)
    call other_energy_short(rs,rsi,rsf,evecp,use_cutoffs,aone)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if ((use_ZSEC.EQV..true.).AND.(seqpolty(rs).eq.'P')) call en_zsec_gl(imol,evecp)
    if (use_POLY.EQV..true.) call en_poly_gl(imol,evecp,aone)
    if ((use_DSSP.EQV..true.).AND.(seqpolty(rs).eq.'P')) call en_dssp_gl(evecp)
    if (use_DREST.EQV..true.) call edrest(evecp)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
!   the long-range part has a different cutoff and is computed separately
!   also, it will use the old atsav's, but compute only those pair terms that changed
!   at the end the fxn updates the atsav's to the sampled conformation (generated in the next step)
    call other_energy_long(rs,rsi,rsf,evecp,use_cutoffs,azero)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   now re-apply the sampling move
    call mcmove_other(which,aone,ct,flipset,rsi,rsf,yep)
!
!   time to get the long-range posterior energies (with updated atsav's)
    call CPU_time(t1)
    call other_energy_long(rs,rsi,rsf,eveca,use_cutoffs,aone)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1            
!
!   get energy difference and determine acceptance
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
!   update or restore
    if (moveok.EQV..true.) then
      acc%nother = acc%nother + 1
      acc%other(rs) = acc%other(rs) + 1
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = esterms(i) + evecd(i)
      end do
      call internal_bound(imol)
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
    else
      call mcmove_other(which,atwo,ct,flipset,rsi,rsf,yep)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
      if (use_POLY.EQV..true.) call restore_poly(imol)
    end if
!
!
!
!
!
  else if (UJconrot_move.EQV..true.) then
!
!   this is the true local CR implementation by Ulmschneider / Jorgensen employing
!   flexible bond lengths 
!
    moveok = .false.
!
!   get a suitable residue for closure (last one)
    rcheck = random()
    call binary_search(ujcrlst%nr,ujcrlst%wt(1:ujcrlst%nr),rcheck,rs)
    rsf = ujcrlst%idx(min(ujcrlst%nr,rs+1))
    imol = molofrs(rsf)
    nresidue = min(floor(random()*(ujmaxsz-ujminsz+1)) + ujminsz,rsf-ujcrlst%idx2(min(ujcrlst%nr,rs+1))+1)
    cur_ujsz = nresidue - 2
    rsi = rsf - cur_ujsz - 1
!
!   number of DOF for residue chain = 5(i-2)
    nDOF = 5*cur_ujsz + 6
    dof = 5*cur_ujsz
!
!   testing
    blas(1) = x(ci(rsf))
    blas(2) = y(ci(rsf))
    blas(3) = z(ci(rsf))
!
!   ----------- Pre-Rotation -----------
!
!   build the DOF(i) matrix for the value of midpoint a(phi,alph...)
!   returned value is square matrix ujmat = I = da/dphi * da/phi
    call buildUJdof(rsi,rsf,dof,ujmat)
!
!   calculate J = c1*(1+c2*I)
    do j=1,dof
      do k=1,dof
         pos_def(j,k) = ujmat(j,k)*UJ_params(1)*UJ_params(2)
      end do
      pos_def(j,j) = pos_def(j,j)+UJ_params(1)
    end do
!
!   find Lt via Cholesky decomp, J=LLt
    call cdecomp(pos_def,MAXUJDOF,dof,cdecomp_upper,.true.)
    call cdecomp(pos_def,MAXUJDOF,dof,cdecomp_lower,.false.)
!
!   weigh Lt by angles vs. dihedral(columns dihe_displace by c3)
    dihed_displace = cur_ujsz*2 + 1
    do j=1,dof
      do k=dihed_displace,dof
        cdecomp_upper(k,j) = cdecomp_upper(k,j)*UJ_params(3)
        cdecomp_lower(j,k) = cdecomp_lower(j,k)*UJ_params(3)
      end do
    end do
!
!   generate random vector dchi, gaussian dist
    do j=1,dof
      gvec(j) = normal()
    end do
!
!   calculate d**2 = gvec_trans*gvec
    dsquare = dot_product(gvec(1:dof),gvec(1:dof))
!
!   calculate delta phi from equation (L^t * delta_phi = gve)
    dphi(:) = 0
    call trisolv(cdecomp_upper,dphi,gvec,MAXUJDOF,dof)
!
!   get reference jacobian
    call jacobian(rsf,jac_a)
!
!   ! store the old values needed for chain closure prior to moving chain
    call uj_refvals(rsf,ccmat,oldvals)
!
!   move prerotation segment
    call mcmove_uj(rsi,rsf,dphi,azero)
!
!   repeat steps ABOVE for new atomic coordinates
    call buildUJdof(rsi,rsf,dof,ujmat)
!
    do j=1,dof
      do k=1,dof
        pos_def_star(j,k)=ujmat(j,k)*UJ_params(1)*UJ_params(2)
      end do
      pos_def_star(j,j) = pos_def_star(j,j)+UJ_params(1)
    end do
    call cdecomp(pos_def_star,MAXUJDOF,dof,cdecomp_upper_star,.true.)
    call cdecomp(pos_def_star,MAXUJDOF,dof,cdecomp_lower_star,.false.)
    do j=1,dof
      do k=dihed_displace,dof
        cdecomp_upper_star(k,j) = cdecomp_upper_star(k,j)*UJ_params(3)
        cdecomp_lower_star(j,k) = cdecomp_lower_star(j,k)*UJ_params(3)
      end do
    end do
!   calculate set of gvec phis based on (gvec_star = L_star^t * delta_phi)
    do j=1,dof
      gvec_star(j) = 0.0
      do k=1,dof
        gvec_star(j) = gvec_star(j) + &
 &          cdecomp_upper_star(j,k)*dphi(k)
      end do
    end do
!   cal d_star**2
    dsquare_star = dot_product(gvec_star(1:dof),gvec_star(1:dof))
!
!   now calculate the biasing probability for the move and reverse move
    det_L = 1.0d0
    do j=1,dof
      det_L = det_L*cdecomp_lower(j,j)
    end do
    det_L_star = 1.0d0
    do j=1,dof
      det_L_star = det_L_star*cdecomp_lower_star(j,j)
    end do
!! ! WARNING : WARNING:: WARNING:: WARNING:: 
!    prob_move = (det_L_star/det_L)*exp(-(dsquare_star-dsquare))
!    prob_move = (det_L_star/det_L)*exp(-(dsquare_star-dsquare))*(jac_b/jac_a)
!
!    prob_move = det_L*exp(-dsquare)
!    prob_back = det_L_star*exp(-dsquare_star)
!
!   restore xyz of segment post prerotation
    call mcmove_uj(rsi,rsf,dphi,afour)
!
!   ----------- Chain Close -----------
!
!   find the displacement angles (dccmat) needed to close chain
    call CPU_time(blas(4))
    call UJchainclose(rsf,ccmat,dccmat,oldvals,findsol)
    call CPU_time(blas(5))
    time_struc = time_struc + blas(5) - blas(4)
!
    do j=1,dof
      dpfmat(j) = dphi(j)
    end do
    do j=dof+1,nDOF
      dpfmat(j) = dccmat(j-dof)
    end do
!
!   restore pre-rotation segment
    call mcmove_uj(rsi,rsf,dphi,athree)
!
    if(findsol.EQV..false.) then
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
      uj_deadcnt = uj_deadcnt + 1
      call mc_acceptance_bailout(mc_acc_crit,istep)
      return
    end if
!
!   move entire internal chain
    call mcmove_uj(rsi,rsf,dpfmat,aone)
!
!   and test for proper closure (testing)
    if ((abs(x(ci(rsf))-blas(1))+abs(y(ci(rsf))-blas(2))+abs(z(ci(rsf))-blas(3)))&
 &        .gt.1.0e-4) then
      blas(4) = 1000.0
!      write(*,*) abs(x(ci(rsf))-blas(1))+abs(y(ci(rsf))-blas(2))+abs(z(ci(rsf))-blas(3))
      write(ilog,*) 'Inexact closure during chain closure algorithm.&
 & This is potentially a bug.'
!      write(ilog,*) findsol,istep
!      write(ilog,*) 'ccmat  : ',ccmat
!      write(ilog,*) 'dccmat : ',dccmat
!      write(ilog,*) 'ovals  : ',oldvals
      call fexit()
    end if
!
! ----------- MC Accept -----------
!
    call jacobian(rsf,jac_b)
!
!   compute posterior short-range (steric) energy
    call CPU_time(t1)
!
!   first lets compute the posterior short-range contributions
    call ujcr_energy_short(rsi,rsf,eveca,use_cutoffs,azero)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rsi),eveca)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rsi),eveca,azero)
    if (use_DSSP.EQV..true.) call en_dssp_gl(eveca)
    if (use_DREST.EQV..true.) call edrest(eveca)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
!   now make use of pre-screening (if so desired)
    if (use_stericscreen.EQV..true.) then
      difft = 0.0
      do i=1,MAXENERGYTERMS
        difft = difft + eveca(i)
      end do
      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
        call mcmove_uj(rsi,rsf,dpfmat,atwo)
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   restore the coordinates to the prior state
    call mcmove_uj(rsi,rsf,dpfmat,atwo)
!
!   now compute the prior energies, the short-range part will already
!   yield the full SAV-information (if relevant)
    call CPU_time(t1)
    call ujcr_energy_short(rsi,rsf,evecp,use_cutoffs,aone)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rsi),evecp)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rsi),evecp,aone)
    if (use_DSSP.EQV..true.) call en_dssp_gl(evecp)
    if (use_DREST.EQV..true.) call edrest(evecp)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
!   the long-range part has a different cutoff and is computed separately
!   also, it will use the old atsavs, but compute only those pair terms that changed
!   at the end the fxn updates the atsavs to the sampled conformation (generated in the next step)
    call ujcr_energy_long(rsi,rsf,evecp,use_cutoffs,azero)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   now re-apply the sampling move
    call mcmove_uj(rsi,rsf,dpfmat,aone)
!
!   time to get the long-range posterior energies (with updated atsav's)
    call CPU_time(t1)
    call ujcr_energy_long(rsi,rsf,eveca,use_cutoffs,aone)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   get energy difference and determine acceptance
    prob_move = (det_L_star/det_L)*exp(-(dsquare_star-dsquare))*(jac_b/jac_a)
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
!   update or restore
    if (moveok.EQV..true.) then
      acc%nujcr = acc%nujcr + 1
      do i=rsi,rsf ! note that in rsf only(!) the phi-angle is changed
        acc%cr(i) = acc%cr(i) + 1
      end do
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = esterms(i) + evecd(i)
      end do
      call internal_bound(molofrs(rsi))
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
    else
      call mcmove_uj(rsi,rsf,dpfmat,atwo)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
    end if
!
!
!
!
  else if (SJconrot_move.EQV..true.) then
!
!   here we will sample a few residues in a biased way to minimize lever-arm effects
!   effectively, this is a more local move, which is somewhat useful, but suffers from
!   non-local energy evaluations
!
    moveok = .false.
!
!   get a set of consecutive residues
    rcheck = random()
    call binary_search(sjcrlst%nr,sjcrlst%wt(1:sjcrlst%nr),rcheck,rs)
    rsf = sjcrlst%idx(min(sjcrlst%nr,rs+1)) + 1
    imol = molofrs(rsf)
    rsi = rsf - nr_crres
!
    if (align_NC.eq.3) then
      if ((rsmol(molofrs(rsi),2)-rsf+1)&
 &.gt.(rsi-rsmol(molofrs(rsi),1))) then
        ct = .true.
      else if (((rsmol(molofrs(rsi),2)-rsf+1)&
 &.eq.(rsi-rsmol(molofrs(rsi),1))).AND.(random().gt.0.5)) then
        ct = .true.
      end if
    else if (align_NC.eq.4) then
      bt = 1.0*(rsmol(molofrs(rsi),2)-0.5*(rsi+rsf-1)+1)/&
 &         (rsmol(molofrs(rsi),2)-rsmol(molofrs(rsi),1)+2)
      if (random().le.bt) ct = .true.
    end if
!
    do i=1,MAXCRDOF
      do j=1,MAXCRDOF
        jacmat(i,j) = 0.0
        gmat(i,j) = 0.0
      end do
    end do
!
!   build the (product) derivative-matrix (target positions vs. DOF) G
    if (cr_mode.eq.2) then
      call buildGMat2(rsi,rsf,gmat)
    else
      call buildGMat(rsi,rsf,gmat)
    end if
!
!   transform to matrix A (A = a/2*(I+bG)
    do i=1,nr_crdof
      do j=1,nr_crdof
        amat(i,j) = 0.5*cr_a*cr_b*gmat(i,j)
        if (i.eq.j) then
          amat(i,j) = amat(i,j) + 0.5*cr_a
        end if
      end do
    end do
!
!   Cholesky-decompose matrix A 
    call cdecomp(amat,MAXCRDOF,nr_crdof,lmat,sayyes)
!
!   now solve a triangular set of equations to get a displacement vector
    call trisolv_sj(lmat,dfy,dgv,MAXCRDOF,nr_crdof)
!
!   get the foward biasing probability
    bias_fwd = lmat(1,1)
!   determinant of A via L
    do i=2,nr_crdof
      bias_fwd = bias_fwd*lmat(i,i)
    end do
!   phiT A phi via psi or explicitly
    dumfwd = 0.0
    do i=1,nr_crdof
      dumfwd = dumfwd + dgv(i)*dgv(i)
    end do
!
!   backup reference coordinates and apply the perturbation given by dfy
!    write(*,*) '0',x(33),x(46)
    call mcmove_cr(rsi,rsf,dfy,azero,ct)
!    write(*,*) '1',x(33),x(46)!phi(rsi+1),psi(rsi+1)
!
!   compute posterior short-range (steric) energy
    call CPU_time(t1)
!
!   first let's compute the posterior short-range contributions
    call cr_energy_short(rsi,rsf,eveca,use_cutoffs,azero,ct)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rsi),eveca)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rsi),eveca,azero)
    if (use_DSSP.EQV..true.) call en_dssp_gl(eveca)
    if (use_DREST.EQV..true.) call edrest(eveca)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
!   now make use of pre-screening (if so desired)
    if (use_stericscreen.EQV..true.) then
      difft = 0.0
      do i=1,MAXENERGYTERMS
        difft = difft + eveca(i)
      end do
      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
        call mcmove_cr(rsi,rsf,dfy,atwo,ct)
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   restore the coordinates to the prior state
    call mcmove_cr(rsi,rsf,dfy,atwo,ct)
!    write(*,*) '2',x(33),x(46)!phi(rsi+1),psi(rsi+1)
!
!   now compute the prior energies, the short-range part will already
!   yield the full SAV-information (if relevant)
    call CPU_time(t1)
    call cr_energy_short(rsi,rsf,evecp,use_cutoffs,aone,ct)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rsi),evecp)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rsi),evecp,aone)
    if (use_DSSP.EQV..true.) call en_dssp_gl(evecp)
    if (use_DREST.EQV..true.) call edrest(evecp)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
!   the long-range part has a different cutoff and is computed separately
!   also, it will use the old atsav's, but compute only those pair terms that changed
!   at the end the fxn updates the atsav's to the sampled conformation (generated in the next step)
    call cr_energy_long(rsi,rsf,evecp,use_cutoffs,azero,ct)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   now re-apply the sampling move
    call mcmove_cr(rsi,rsf,dfy,aone,ct)
!    write(*,*) '3',x(33),x(46)!phi(rsi+1),psi(rsi+1)
!
!   time to get the long-range posterior energies (with updated atsav's)
    call CPU_time(t1)
    call cr_energy_long(rsi,rsf,eveca,use_cutoffs,aone,ct)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   in order to be able to evaluate the acceptance prob. we need to to know the reverse bias as well
!   unfortunately there is no indirect way to get this, so again get G-matrix (of the perturbed state)
    if (cr_mode.eq.2) then
      call buildGMat2(rsi,rsf,gmat)
    else
      call buildGMat(rsi,rsf,gmat)
    end if
!
!   convert to A 
    do i=1,nr_crdof
      do j=1,nr_crdof
       amat(i,j) = 0.5*cr_a*cr_b*gmat(i,j)
        if (i.eq.j) then
          amat(i,j) = amat(i,j) + 0.5*cr_a
        end if
      end do
    end do
!
!   Cholesky-decompose matrix A (for calculation of determinant)
    call cdecomp(amat,MAXCRDOF,nr_crdof,lmat,sayyes)
!
!   we know the displacement vector already (-dphi), so need for LE-solving 
    do i=1,nr_crdof
      dfy(i) = -dfy(i)
    end do
!
!   get the backward biasing probability, now the term in exponential is more complicated
    bias_bwd = lmat(1,1)
!   determinant of A via L
    do i=2,nr_crdof
      bias_bwd = bias_bwd*lmat(i,i)
    end do
!   phiT A phi explicitly (matrix products)
    dumbwd = 0.0
    do i=1,nr_crdof
      dum = 0.0
      do j=1,nr_crdof
        dum = dum + dfy(j)*amat(i,j)
      end do
      dumbwd = dumbwd + dfy(i)*dum
    end do
!
!   get energy difference and determine acceptance
    prob_move = (bias_bwd/bias_fwd)*exp(-(dumbwd-dumfwd))
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
!   update or restore
    if (moveok.EQV..true.) then
      acc%nsjcr = acc%nsjcr + 1
      do i=rsi,rsf-1
        acc%cr(i) = acc%cr(i) + 1
      end do
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = esterms(i) + evecd(i)
      end do
      call internal_bound(molofrs(rsi))
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
    else
      call mcmove_cr(rsi,rsf,dfy,atwo,ct)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
    end if
!
!
!
!
  else if (torcr_move.EQV..true.) then
!
!   concerted rotation with exact chain closure, constant bond angles (skip omega torsions)
!
! DOES THIS USE OMEGA IN PREROTATION OR NOT? CHECK PARTICLEFLUC
    moveok = .false.
!
!   get a set of consecutive residues
    rcheck = random()
    call binary_search(djcrlst%nr,djcrlst%wt(1:djcrlst%nr),rcheck,rs)
    rsf = djcrlst%idx(min(djcrlst%nr,rs+1))
    imol = molofrs(rsf)
    dof = min(floor(random()*(torcrmaxsz2-torcrminsz2+1)) + torcrminsz2,3*(rsf-djcrlst%idx2(min(djcrlst%nr,rs+1))+1)-9)
    rsi = rsf - 2 - floor((dof+1.5)/3.0)
    if (rsi.eq.rsmol(imol,1)) then
      if (mod(dof,3).eq.0) dof = dof - 1 ! for any terminal polypeptide residue phi is missing
      if (seqpolty(rsi).ne.'P') then
        if (mod(dof,3).eq.2) dof = dof - 1 ! for cap, psi is missing, too
        rsi = rsi + 1
      end if
    end if
!   we kick phi-angles out for cyclic residues (note that the total always remains >= 1)
    if (rsi.lt.(rsf-2)) then
      if (seqflag(rsi).eq.5) then
        if (mod(dof,3).le.1) dof = dof - 1
      end if
    end if
    do rs=rsi+1,rsf-3
      if (seqflag(rs).eq.5) dof = dof - 1
    end do
!
!    write(*,*) 'in',z(at(rsf)%sc(4)),z(at(rsf-1)%sc(4)),z(at(rsf-2)%sc(4))
!    write(*,*) 'jn',z(at(rsf)%sc(3)),z(at(rsf-1)%sc(3)),z(at(rsf-2)%sc(3))
!    write(*,*) imol,rsi,rsf,dof,torcrmaxsz2,torcrminsz2
!
!   number of DOF for residue chain = 3(i-2)
    nDOF = dof + 6
!
!   test for closure and prerotation
    blas(1) = x(ci(rsf))
    blas(2) = y(ci(rsf))
    blas(3) = z(ci(rsf))
!
    if (use_POLY.EQV..true.) then
      call makeref_poly(molofrs(rsi))
    end if
!
!   ----------- Pre-Rotation -----------
!
!   build the DOF(i) matrix for the value of midpoint a(phi,alph...)
!   returned value is square matrix torcrmat = I = da/dphi * da/phi
    call buildtorcrdof_dj(rsi,rsf,dof,torcrmat,omepos,nome)
!
!   calculate J = c1*(1+c2*I)
    do j=1,dof
      do k=1,dof
        pos_def(j,k) = torcrmat(j,k)*UJ_params(1)*UJ_params(2)
      end do
      pos_def(j,j) = pos_def(j,j)+UJ_params(1)
    end do

!
!   find Lt via Cholesky decomp, J=LLt
    call cdecomp(pos_def,MAXUJDOF,dof,cdecomp_upper,.true.)
    call cdecomp(pos_def,MAXUJDOF,dof,cdecomp_lower,.false.)
!
!   re-weight omegas to reduce impact
    dihed_displace = 3
    do j=1,dof
!      do k=omepos,dof,dihed_displace
      do kk=1,nome
        k = omepos(kk)
        cdecomp_upper(k,j) = cdecomp_upper(k,j)*UJ_params(6)
        cdecomp_lower(j,k) = cdecomp_lower(j,k)*UJ_params(6)
      end do
    end do
!
!   get reference jacobian
    call jacobian_torcr_dj(rsf,jac_a)
!
!   store the old values needed for chain closure prior to moving chain
    call torcr_refvals_dj(rsf,ccmat,oldvals3)
!
!   in the "safe" mode, assemble half-set of solutions with no pre-rotation (includes original
!   conformation for detailed balance)
    if (torcrmode.eq.1) then
!
      dphi(1:dof) = 0.0
!     empty move to get backup arrays populated correctly
      call mcmove_torcr_dj(rsi,rsf,dof,dphi,azero)
      call mcmove_torcr_dj(rsi,rsf,dof,dphi,afour)
      solmatp(:,:) = 0.0
      call CPU_time(blas(10))
      call torcrchainclose_dj(rsi,rsf,dof,dphi,ccmat,solmatp,oldvals3,solcnt2,jacvp)
      call proline_cr(rsf,puckdofs,solmatp,jacvp,solcnt2,curpksp,aone)
      if (solcnt2.eq.0) then
        dj_wrncnt(4) = dj_wrncnt(4) + 1
        if (dj_wrncnt(4).eq.dj_wrnlmt(4)) then
          write(ilog,*) 'Warning. Did not find original solution for Dinner-Ulmschneider concerted rotation&
 & algorithm without omega bond sampling. This indicates that the root search is not failsafe. &
 &If this problem occurs frequently, the simulation is most likely producing&
 & biased or otherwise corrupted results.'
          write(ilog,*) 'This was warning #',dj_wrncnt(4),' of this type not all of which may be displayed.'
          if (10.0*dj_wrnlmt(4).gt.0.5*HUGE(dj_wrnlmt(4))) then
            dj_wrncnt(4) = 0
          else
            dj_wrnlmt(4) = dj_wrnlmt(4)*10
          end if
        end if
      end if
      call CPU_time(blas(11))
      time_struc = time_struc + blas(11) - blas(10)
    end if
!
    solucount = 0
    exactcrtries = 0

!
    do while ((solucount.eq.0).AND.(exactcrtries.lt.MAXCRTRIES))
!
      exactcrtries = exactcrtries + 1
!
!     generate random vector dchi, gaussian dist
      do j=1,dof
        gvec(j) = normal()
      end do
!
!     calculate d**2 = gvec_trans*gvec
      dsquare = 0.0
      do j=1,dof
        dsquare = dsquare + gvec(j)*gvec(j)
      end do
!
!     calculate delta phi from equation (L^t * delta_phi = gve)
      dphi(:) = 0
      call trisolv(cdecomp_upper,dphi,gvec,MAXUJDOF,dof)
!
!     move prerotation segment
      call mcmove_torcr_dj(rsi,rsf,dof,dphi,azero)
!      if (dof.eq.10) write(*,*) abs(x(cai(rsf-2))-blas(21))+abs(y(cai(rsf-2))-blas(22))+abs(z(cai(rsf-2))-blas(23))
!
!     repeat steps ABOVE for new atomic coordinates
      call buildtorcrdof_dj(rsi,rsf,dof,torcrmat,omepos,nome)
!
      do j=1,dof
        do k=1,dof
          pos_def_star(j,k)=torcrmat(j,k)*UJ_params(1)*UJ_params(2)
        end do
        pos_def_star(j,j) = pos_def_star(j,j)+UJ_params(1)
      end do
      call cdecomp(pos_def_star,MAXUJDOF,dof,cdecomp_upper_star,.true.)
      call cdecomp(pos_def_star,MAXUJDOF,dof,cdecomp_lower_star,.false.)
!
!     re-weight omegas to reduce impact
      dihed_displace = 3
      do j=1,dof
!      do k=omepos,dof,dihed_displace
        do kk=1,nome
          k = omepos(kk)
          cdecomp_upper_star(k,j) = cdecomp_upper_star(k,j)*UJ_params(6)
          cdecomp_lower_star(j,k) = cdecomp_lower_star(j,k)*UJ_params(6)
        end do
      end do
!
!     calculate set of gvec phis based on (gvec_star = L_star^t * delta_phi)
      do j=1,dof
        gvec_star(j) = 0.0
        do k=1,dof
          gvec_star(j) = gvec_star(j) + cdecomp_upper_star(j,k)*dphi(k)
        end do
      end do
!     cal d_star**2
      dsquare_star = 0.0
      do j=1,dof
        dsquare_star = dsquare_star + gvec_star(j)*gvec_star(j)
      end do
!
!     now calculate the biasing probability for the move and reverse move
      det_L = 1.0
      do j=1,dof
        det_L = det_L*cdecomp_lower(j,j) ! remember given pos.def. A=L(trans)*L, then det(A)=(Prod.(diagonal(L)))^2
      end do
      det_L_star = 1.0
      do j=1,dof
        det_L_star = det_L_star*cdecomp_lower_star(j,j)
      end do
      prob_move = (det_L_star/det_L)*exp(-(dsquare_star-dsquare))
!      call trisolv(cdecomp_upper_star,dum1,gvec_star,MAXUJDOF,dof)
!      if (sum(abs(dphi(1:dof)-dum1(1:dof))).ge.1.0e-9) write(*,*) dof,prob_move,prob_back
!
!     restore xyz of segment post prerotation
      call mcmove_torcr_dj(rsi,rsf,dof,dphi,afour)
!
!     ----------- Chain Close -----------
!
!     find the displacement angles (dpfmat) needed to close chain
      call CPU_time(blas(10))
      solmat(:,:) = 0.0
      call torcrchainclose_dj(rsi,rsf,dof,dphi,ccmat,solmat,oldvals3,solucount,jacv)
      call proline_cr(rsf,puckdofs,solmat,jacv,solucount,curpks,aone)
      if (torcrmode.eq.1) then
        if (solcnt2.gt.0) then
          jacv(solucount+1:solucount+solcnt2) = jacvp(1:solcnt2)
          solmat(solucount+1:solucount+solcnt2,:) = solmatp(1:solcnt2,:)
          curpks(solucount+1:solucount+solcnt2,:,:) = curpksp(1:solcnt2,:,:)
        end if
        if (solucount+solcnt2.gt.0) then
          call picksolu_dj(rsi,rsf,dof,dphi,solmat,curpks,solucount,solcnt2,jacv,dpfmat2,which,1.0,1.0)!det_L,det_L_star)
        end if
      else
        if (solucount.gt.0) then
          call picksolu_dj(rsi,rsf,dof,dphi,solmat,curpks,solucount,solcnt2,jacv,dpfmat2,which,det_L,det_L_star)
          jac_b = jacv(which)
        end if
      end if
      call CPU_time(blas(11))
      time_struc = time_struc + blas(11) - blas(10)
!
      if (solucount.eq.0) then
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
      end if
!     do not recycle for "safe" mode
      if (torcrmode.eq.1) exit

    end do ! while solucount is 0 and exactcrtries less than MAXCRTRIES

!
    dj_totcnt = dj_totcnt + exactcrtries
!
    if (torcrmode.eq.1) then
      if ((solucount.eq.0).AND.(solcnt2.eq.0)) then
        dj_wrncnt(1) = dj_wrncnt(1) + 1
        if (dj_wrncnt(1).eq.dj_wrnlmt(1)) then
          write(ilog,*) 'Warning. Did not find any valid solution for Dinner-Ulmschneider concerted rotation&
 & algorithm without omega bond sampling. This is indicative of a bug or of bad search settings for the&
 & closure algorithm (decrease stepsize). If this problem occurs frequently, the simulation is probably&
 & producing biased or otherwise corrupted results.'
          write(ilog,*) 'This was warning #',dj_wrncnt(1),' of this type not all of which may be displayed.'
          if (10.0*dj_wrnlmt(1).gt.0.5*HUGE(dj_wrnlmt(1))) then
            dj_wrncnt(1) = 0
          else
            dj_wrnlmt(1) = dj_wrnlmt(1)*10
          end if
        end if
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
        dj_deadcnt = dj_deadcnt + 1
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    else
      if ((exactcrtries.eq.MAXCRTRIES).AND.(solucount.eq.0)) then
        dj_wrncnt(1) = dj_wrncnt(1) + 1
        if (dj_wrncnt(1).eq.dj_wrnlmt(1)) then
          write(ilog,*) 'Warning. Did not find a valid solution for Dinner-Ulmschneider concerted rotation&
 & algorithm without omega bond sampling after ',MAXCRTRIES,' pre-rotation attempts.'
          write(ilog,*) 'This is indicative of a bug, of too large pre-rotation stepsizes or of bad search settings &
 & for the closure algorithm (decrease stepsize). If this problem occurs frequently, the simulation is probably&
 & producing biased or otherwise corrupted results.'
          write(ilog,*) 'This was warning #',dj_wrncnt(1),' of this type not all of which may be displayed.'
          if (10.0*dj_wrnlmt(1).gt.0.5*HUGE(dj_wrnlmt(1))) then
            dj_wrncnt(1) = 0
          else
            dj_wrnlmt(1) = dj_wrnlmt(1)*10
          end if
        end if
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
        dj_deadcnt = dj_deadcnt + 1
!        mvcnt%nfy = mvcnt%nfy + 1
!        j = floor(random()*3.0)
!        call pivot_bailout(rsf-j,istep)
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
!   and test for proper closure (testing)
    if ((abs(x(ci(rsf))-blas(1))+abs(y(ci(rsf))-blas(2))+abs(z(ci(rsf))-blas(3)))&
 &          .gt.1.0e-4) then
      write(ilog,*) 'Inexact closure during chain closure algorithm.&
 & This is potentially a bug. Please record following information.'
      write(ilog,*) 'Step #',istep,' and final res. ',rsf
      call fexit()
    end if
!
 444 format(i2,1x,3(g12.6,1x))
!   evaluate pre-rotation biases
    if (torcrmode.eq.1) then
      if (which.le.solucount) then
        prob_move = (det_L_star/det_L)*exp(-(dsquare_star-dsquare))
      else
        prob_move = 1.0
      end if
!      prob_back = exp(-dsquare_star)
    else
      prob_move = (det_L_star/det_L)*exp(-(dsquare_star-dsquare))*(jac_b/jac_a)
!      prob_back = det_L_star*exp(-dsquare_star)
    end if
!
!   compute posterior short-range (steric) energy
    call CPU_time(t1)
!
!   first lets compute the posterior short-range contributions
    if (use_POLY.EQV..true.) then
      call update_rigid(molofrs(rsi))
    end if
    call ujcr_energy_short(rsi,rsf,eveca,use_cutoffs,azero)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rsi),eveca)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rsi),eveca,azero)
    if (use_DSSP.EQV..true.) call en_dssp_gl(eveca)
    if (use_DREST.EQV..true.) call edrest(eveca)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
!   now make use of pre-screening (if so desired)
    if (use_stericscreen.EQV..true.) then
      difft = 0.0
      do i=1,MAXENERGYTERMS
        difft = difft + eveca(i)
      end do
      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
        call mcmove_torcr_dj(rsi,rsf,dof,dpfmat2,atwo)
        call mcmove_torcr_aux(rsi,rsf,dof,curpks,atwo,aone)
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   restore the coordinates to the prior state
    call mcmove_torcr_dj(rsi,rsf,dof,dpfmat2,atwo)
    call mcmove_torcr_aux(rsi,rsf,which,curpks,atwo,aone)
!
!   now compute the prior energies, the short-range part will already
!   yield the full SAV-information (if relevant)
    call CPU_time(t1)
    call ujcr_energy_short(rsi,rsf,evecp,use_cutoffs,aone)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rsi),evecp)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rsi),evecp,aone)
    if (use_DSSP.EQV..true.) call en_dssp_gl(evecp)
    if (use_DREST.EQV..true.) call edrest(evecp)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
!   the long-range part has a different cutoff and is computed separately
!   also, it will use the old atsavs, but compute only those pair terms that changed
!   at the end the fxn updates the atsavs to the sampled conformation (generated in the next step)
    call ujcr_energy_long(rsi,rsf,evecp,use_cutoffs,azero)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   now re-apply the sampling move
    call mcmove_torcr_dj(rsi,rsf,dof,dpfmat2,aone)
    call mcmove_torcr_aux(rsi,rsf,which,curpks,aone,aone)
!
!   time to get the long-range posterior energies (with updated atsav's)
    call CPU_time(t1)
    call ujcr_energy_long(rsi,rsf,eveca,use_cutoffs,aone)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   get energy difference and determine acceptance
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
!   update or restore
    if (moveok.EQV..true.) then
      acc%ndjcr = acc%ndjcr + 1
      do i=rsi,rsf
        acc%cr(i) = acc%cr(i) + 1
      end do
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = esterms(i) + evecd(i)
      end do
      call internal_bound(molofrs(rsi))
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
    else
      call mcmove_torcr_dj(rsi,rsf,dof,dpfmat2,atwo)
      call mcmove_torcr_aux(rsi,rsf,which,curpks,atwo,aone)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
    end if
!    rcheck = getbang(at(rsf)%sc(1),at(rsf)%bb(2),at(rsf)%sc(2))
!    if ((rcheck.gt.150.0).OR.(rcheck.lt.80.0)) write(*,*) rcheck
!
!
!
!
  else if (torcr_move_omega.EQV..true.) then
!
!   torsion conrot, include omega angles
    moveok = .false.
!
!   get a set of consecutive residues
    rcheck = random()
    call binary_search(docrlst%nr,docrlst%wt(1:docrlst%nr),rcheck,rs)
    rsf = docrlst%idx(min(docrlst%nr,rs+1))
    imol = molofrs(rsf)
    dof = min(floor(random()*(torcrmaxsz-torcrminsz+1)) + torcrminsz,3*(rsf-docrlst%idx2(min(docrlst%nr,rs+1))+1)-6)
    rsi = rsf - 1 - floor((dof+2.5)/3.0)
    if ((rsi.eq.rsmol(imol,1)).AND.(mod(dof,3).eq.0)) then ! omega and phi are missing
      dof = dof - 2
    else if ((rsi.eq.rsmol(imol,1)).AND.(mod(dof,3).eq.2)) then ! phi is missing
      dof = dof - 1
    end if
!   we kick phi-angles out for cyclic residues (note that the total always remains >= 1)
    do rs=rsi+1,rsf-2
      if (seqflag(rs).eq.5) dof = dof - 1
    end do
!
!    write(*,*) imol,rsi,rsf,dof,istep
!   number of DOF for residue chain = 3(i-2)
    nDOF = dof + 6
!
!   test for closure and prerotation
    blas(1) = x(ci(rsf))
    blas(2) = y(ci(rsf))
    blas(3) = z(ci(rsf))
!
    if (use_POLY.EQV..true.) then
      call makeref_poly(molofrs(rsi))
    end if
!
!   build the DOF(i) matrix for the value of midpoint a(phi,alph...)
!   returned value is square matrix torcrmat = I = da/dphi * da/phi
    call buildtorcrdof_do(rsi,rsf,dof,torcrmat,omepos,nome)
!
!   calculate J = c1*(1+c2*I)
    do j=1,dof
      do k=1,dof
         pos_def(j,k) = torcrmat(j,k)*UJ_params(1)*UJ_params(2)
      end do
      pos_def(j,j) = pos_def(j,j)+UJ_params(1)
    end do
!
!   find Lt via Cholesky decomp, J=LLt
    call cdecomp(pos_def,MAXUJDOF,dof,cdecomp_upper,.true.)
    call cdecomp(pos_def,MAXUJDOF,dof,cdecomp_lower,.false.)
!
!   weigh omega angles by c3
    dihed_displace = 3
    do j=1,dof
!      do k=omepos,dof,dihed_displace
      do kk=1,nome
        k = omepos(kk)
        cdecomp_upper(k,j) = cdecomp_upper(k,j)*UJ_params(6)
        cdecomp_lower(j,k) = cdecomp_lower(j,k)*UJ_params(6)
      end do
    end do
!
!   get reference jacobian
    call jacobian_torcr_do(rsf,jac_a)
!
!   store the old values needed for chain closure prior to moving chain
    call torcr_refvals_do(rsf,ccmat,oldvals2)
!
    if (torcrmode.eq.1) then
      dphi(1:dof) = 0.0
!     empty move to get backup arrays populated correctly
      call mcmove_torcr_do(rsi,rsf,dof,dphi,azero)
      call mcmove_torcr_do(rsi,rsf,dof,dphi,afour)
      solmatp(:,:) = 0.0
      call CPU_time(blas(10))
      call torcrchainclose_do(rsi,rsf,dof,dphi,ccmat,solmatp,oldvals2,solcnt2,jacvp)
      call proline_cr(rsf,puckdofs,solmatp,jacvp,solcnt2,curpksp,atwo)
!
      if (solcnt2.eq.0) then
        do_wrncnt(4) = do_wrncnt(4) + 1
        if (do_wrncnt(4).eq.do_wrnlmt(4)) then
          write(ilog,*) 'Warning. Did not find original solution for Dinner-Ulmschneider concerted rotation&
 & algorithm with omega bond sampling. This indicates that the root search is not failsafe. &
 &If this problem occurs frequently, the simulation is most likely producing&
 & biased or otherwise corrupted results.'
          write(ilog,*) 'This was warning #',do_wrncnt(4),' of this type not all of which may be displayed.'
          if (10.0*do_wrnlmt(4).gt.0.5*HUGE(do_wrnlmt(4))) then
            do_wrncnt(4) = 0
          else
            do_wrnlmt(4) = do_wrnlmt(4)*10
          end if
        end if
      end if
      call CPU_time(blas(11))
      time_struc = time_struc + blas(11) - blas(10)
    end if
!
    solucount = 0
    exactcrtries = 0
!
    do while ((solucount.eq.0).AND.(exactcrtries.lt.MAXCRTRIES))
!
      exactcrtries = exactcrtries + 1

!     generate random vector dchi, gaussian dist
      do j=1,dof
        gvec(j) = normal()
      end do
!
!     calculate d**2 = gvec_trans*gvec
      dsquare = 0.0d0
      do j=1,dof
        dsquare = dsquare + gvec(j)*gvec(j)
      end do
!
!     calculate delta phi from equation (L^t * delta_phi = gve)
      dphi(:) = 0.0
      call trisolv(cdecomp_upper,dphi,gvec,MAXUJDOF,dof)
!      if (dof.ge.8) write(*,*) dof,dphi(1:dof)/gvec(1:dof),(RADIAN/(1.0*dof))*sum(abs(dphi(1:dof)))
!
!     move prerotation segment
      call mcmove_torcr_do(rsi,rsf,dof,dphi,azero)
!
!     repeat steps ABOVE for new atomic coordinates
      call buildtorcrdof_do(rsi,rsf,dof,torcrmat,omepos,nome)
!
      do j=1,dof
        do k=1,dof
          pos_def_star(j,k)=torcrmat(j,k)*UJ_params(1)*UJ_params(2)
        end do
        pos_def_star(j,j) = pos_def_star(j,j)+UJ_params(1)
      end do
      call cdecomp(pos_def_star,MAXUJDOF,dof,cdecomp_upper_star,.true.)
      call cdecomp(pos_def_star,MAXUJDOF,dof,cdecomp_lower_star,.false.)
      do j=1,dof
!        do k=omepos,dof,dihed_displace
        do kk=1,nome
          k = omepos(kk)
          cdecomp_upper_star(k,j) = cdecomp_upper_star(k,j)*UJ_params(6)
          cdecomp_lower_star(j,k) = cdecomp_lower_star(j,k)*UJ_params(6)
        end do
      end do
!
!     calculate set of gvec phis based on (gvec_star = L_star^t * delta_phi)
      do j=1,dof
        gvec_star(j) = 0.0
        do k=1,dof
          gvec_star(j) = gvec_star(j) + cdecomp_upper_star(j,k)*dphi(k)
        end do
      end do
!     cal d_star**2
      dsquare_star = 0.0d0
      do j=1,dof
        dsquare_star = dsquare_star + gvec_star(j)*gvec_star(j)
      end do
!
!     now calculate the biasing probability for the move and reverse move
      det_L = 1.0d0
      do j=1,dof
        det_L = det_L*cdecomp_lower(j,j)
      end do
      det_L_star = 1.0d0
      do j=1,dof
        det_L_star = det_L_star*cdecomp_lower_star(j,j)
      end do
!
!     restore xyz of segment post prerotation
      call mcmove_torcr_do(rsi,rsf,dof,dphi,afour)
!
!     ----------- Chain Close -----------
!
!     find the displacement angles (dpfmat) needed to sample and close the chain
      call CPU_time(blas(10))
      solmat(:,:) = 0.0
      call torcrchainclose_do(rsi,rsf,dof,dphi,ccmat,solmat,oldvals2,solucount,jacv)
      call proline_cr(rsf,puckdofs,solmat,jacv,solucount,curpks,atwo)
      if (torcrmode.eq.1) then
        if (solcnt2.gt.0) then
          jacv(solucount+1:solucount+solcnt2) = jacvp(1:solcnt2)
          solmat(solucount+1:solucount+solcnt2,:) = solmatp(1:solcnt2,:)
          curpks(solucount+1:solucount+solcnt2,:,:) = curpksp(1:solcnt2,:,:)
        end if
        if (solucount+solcnt2.gt.0) then
          call picksolu_do(rsi,rsf,dof,dphi,solmat,curpks,solucount,solcnt2,jacv,dpfmat2,which,1.0,1.0)!,det_L,det_L_star)
        end if
      else
        if (solucount.gt.0) then
          call picksolu_do(rsi,rsf,dof,dphi,solmat,curpks,solucount,solcnt2,jacv,dpfmat2,which,1.0,1.0)!,det_L,det_L_star)
          jac_b = jacv(which)
        end if
      end if
      call CPU_time(blas(11))
      time_struc = time_struc + blas(11) - blas(10)
!
      if (solucount.eq.0) then
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
      end if
!
!     do not recycle for "safe" mode
      if (torcrmode.eq.1) exit

    end do ! while solucount is 0 and exactcrtries less than MAXCRTRIES
!
    do_totcnt = do_totcnt + exactcrtries
!
    if (torcrmode.eq.1) then
      if ((solucount+solcnt2).eq.0) then
        do_wrncnt(1) = do_wrncnt(1) + 1
        if (do_wrncnt(1).eq.do_wrnlmt(1)) then
          write(ilog,*) 'Warning. Did not find any valid solution for Dinner-Ulmschneider concerted rotation&
 & algorithm with omega bond sampling. This is indicative of a bug or of bad search settings for the closure&
 & algorithm (decrease stepsize). If this problem occurs frequently, the simulation is probably producing&
 & biased or otherwise corrupted results.'
          write(ilog,*) 'This was warning #',do_wrncnt(1),' of this type not all of which may be displayed.'
          if (10.0*do_wrnlmt(1).gt.0.5*HUGE(do_wrnlmt(1))) then
            do_wrncnt(1) = 0
          else
            do_wrnlmt(1) = do_wrnlmt(1)*10
          end if
        end if
        do_deadcnt = do_deadcnt + 1
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    else
     if ((solucount.eq.0).AND.(exactcrtries.eq.MAXCRTRIES)) then
        do_wrncnt(1) = do_wrncnt(1) + 1
        if (do_wrncnt(1).eq.do_wrnlmt(1)) then
          write(ilog,*) 'Warning. Did not find a valid solution for Dinner-Ulmschneider concerted rotation&
 & algorithm with omega bond sampling after ',MAXCRTRIES,' pre-rotation attempts.'
          write(ilog,*) 'This is indicative of a bug, of too large pre-rotation stepsizes or of bad search &
 &settings for the closure algorithm (decrease stepsize). If this problem occurs frequently, the simulation is&
 & probably producing biased or otherwise corrupted results.'
          write(ilog,*) 'This was warning #',do_wrncnt(1),' of this type not all of which may be displayed.'
          if (10.0*do_wrnlmt(1).gt.0.5*HUGE(do_wrnlmt(1))) then
            do_wrncnt(1) = 0
          else
            do_wrnlmt(1) = do_wrnlmt(1)*10
          end if
        end if
        do_deadcnt = do_deadcnt + 1
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
!        call mcmove_torcr_do(rsi,rsf,dof,dphi,atwo)
!        k = floor(random()*3.0)
!        if (k.eq.2) then
!          mvcnt%nomega = mvcnt%nomega + 1
!          j = floor(random()*2.0)
!          call omega_bailout(rsf-j,istep)
!        else
!          mvcnt%nfy = mvcnt%nfy + 1
!          j = floor(random()*2.0)
!          call pivot_bailout(rsf-j,istep)
!        end if
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
!   and test for proper closure (testing)
    if ((abs(x(ci(rsf))-blas(1))+abs(y(ci(rsf))-blas(2))+abs(z(ci(rsf))-blas(3)))&
 &          .gt.1.0e-4) then
      write(ilog,*) 'Inexact closure during chain closure algorithm.&
 & This is potentially a bug. Please record following information.'
      write(ilog,*) 'Step #',istep,' and final res. ',rsf
      call fexit()
    end if
!
!   evaluate pre-rotation biases
    if (torcrmode.eq.1) then
      if (which.le.solucount) then
        prob_move = (det_L_star/det_L)*exp(-(dsquare_star-dsquare))
      else
        prob_move = 1.0
      end if
!      prob_back = exp(-dsquare_star)
    else
      prob_move = (det_L_star/det_L)*exp(-(dsquare_star-dsquare))*(jac_b/jac_a)
!      prob_back = det_L_star*exp(-dsquare_star)
    end if
!
!   compute posterior short-range (steric) energy
    call CPU_time(t1)
!
!   first lets compute the posterior short-range contributions
    if (use_POLY.EQV..true.) then
      call update_rigid(molofrs(rsi))
    end if
    call ujcr_energy_short(rsi,rsf,eveca,use_cutoffs,azero)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rsi),eveca)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rsi),eveca,azero)
    if (use_DSSP.EQV..true.) call en_dssp_gl(eveca)
    if (use_DREST.EQV..true.) call edrest(eveca)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
!   now make use of pre-screening (if so desired)
    if (use_stericscreen.EQV..true.) then
      difft = 0.0
      do i=1,MAXENERGYTERMS
        difft = difft + eveca(i)
      end do
      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
        call mcmove_torcr_do(rsi,rsf,dof,dpfmat2,atwo)
        call mcmove_torcr_aux(rsi,rsf,which,curpks,atwo,atwo)
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   restore the coordinates to the prior state
    call mcmove_torcr_do(rsi,rsf,dof,dpfmat2,atwo)
    call mcmove_torcr_aux(rsi,rsf,which,curpks,atwo,atwo)
!
!   now compute the prior energies, the short-range part will already
!   yield the full SAV-information (if relevant)
    call CPU_time(t1)
    call ujcr_energy_short(rsi,rsf,evecp,use_cutoffs,aone)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rsi),evecp)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rsi),evecp,aone)
    if (use_DSSP.EQV..true.) call en_dssp_gl(evecp)
    if (use_DREST.EQV..true.) call edrest(evecp)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
!   the long-range part has a different cutoff and is computed separately
!   also, it will use the old atsavs, but compute only those pair terms that changed
!   at the end the fxn updates the atsavs to the sampled conformation (generated in the next step)
    call ujcr_energy_long(rsi,rsf,evecp,use_cutoffs,azero)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   now re-apply the sampling move
    call mcmove_torcr_do(rsi,rsf,dof,dpfmat2,aone)
    call mcmove_torcr_aux(rsi,rsf,which,curpks,aone,atwo)
!
!   time to get the long-range posterior energies (with updated atsav's)
    call CPU_time(t1)
    call ujcr_energy_long(rsi,rsf,eveca,use_cutoffs,aone)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   get energy difference and determine acceptance
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
!   update or restore
    if (moveok.EQV..true.) then
      acc%ndocr = acc%ndocr + 1
      do i=rsi,rsf
        acc%cr(i) = acc%cr(i) + 1
      end do
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = esterms(i) + evecd(i)
      end do
      call internal_bound(molofrs(rsi))
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
    else
      call mcmove_torcr_do(rsi,rsf,dof,dpfmat2,atwo)
      call mcmove_torcr_aux(rsi,rsf,which,curpks,atwo,atwo)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
    end if
!
!
!
!
  else if (nuccr_move.EQV..true.) then
!
    moveok = .false.
!
!   get a set of consecutive residues
    rcheck = random()
    call binary_search(nuccrlst%nr,nuccrlst%wt(1:nuccrlst%nr),rcheck,rs)
    rsf = nuccrlst%idx(min(nuccrlst%nr,rs+1))
    imol = molofrs(rsf)
    dof = min(floor(random()*(nuccrmaxsz-nuccrminsz+1)) + nuccrminsz,5*(rsf-nuccrlst%idx2(min(nuccrlst%nr,rs+1))+1)-7)
    rsi = rsf - 1 - floor((dof+1.5)/5.0)
    if ((seqflag(rsi).eq.24).AND.(mod(dof,5).lt.4)) then
      dof = 5*(dof/5) ! integer division
    else if (rsi.eq.rsmol(imol,1)) then
      if (mod(dof,5).eq.3) dof = dof - 1
    end if
    nDOF = dof + 6
!
    blas(1) = x(nuci(rsf,6))
    blas(2) = y(nuci(rsf,6))
    blas(3) = z(nuci(rsf,6))
!
    if (use_POLY.EQV..true.) then
      call makeref_poly(molofrs(rsi))
    end if
!
!   --------- Pre-Rotation ----------
!
    call buildnuccrdof(rsi,rsf,dof,nuccrmat)
!
!   calculate J = c1*(1+c2*I)
    do j=1,dof
      do k=1,dof
        pos_def(j,k) = nuccrmat(j,k)*UJ_params(1)*UJ_params(2)
      end do
      pos_def(j,j) = pos_def(j,j)+UJ_params(1)
    end do
!
!   find Lt via Cholesky decomp, J=LLt
    call cdecomp(pos_def,MAXUJDOF,dof,cdecomp_upper,.true.)
    call cdecomp(pos_def,MAXUJDOF,dof,cdecomp_lower,.false.)
!
!   get reference jacobian
    call jacobian_nuccr(rsf,jac_a)
!
!   store the old values needed for chain closure prior to moving chain
    call nuccr_refvals(rsf,ccmat_nuc,oldvals_nuc)
!
    if (torcrmode.eq.1) then
      dphi(1:dof) = 0.0
!     empty move to get backup arrays populated correctly
      call mcmove_nuccr(rsi,rsf,dof,dphi,azero)
      call mcmove_nuccr(rsi,rsf,dof,dphi,afour)
      solmatp(:,:) = 0.0
      solmat(:,:) = 0.0
      call CPU_time(blas(10))
      call nuccrchainclose(rsi,rsf,dof,dphi,ccmat_nuc,solmatp,oldvals_nuc,solcnt2,jacvp)
!
      if (solcnt2.eq.0) then
        nc_wrncnt(4) = nc_wrncnt(4) + 1
        if (nc_wrncnt(4).eq.nc_wrnlmt(4)) then
          write(ilog,*) 'Warning. Did not find original solution for Dinner-Ulmschneider concerted rotation&
 & algorithm for polynucleotides. This indicates that the root search is not failsafe. &
 &If this problem occurs frequently, the simulation is most likely producing&
 & biased or otherwise corrupted results.'
          write(ilog,*) 'This was warning #',nc_wrncnt(4),' of this type not all of which may be displayed.'
          if (10.0*nc_wrnlmt(4).gt.0.5*HUGE(nc_wrnlmt(4))) then
            nc_wrncnt(4) = 0
          else
            nc_wrnlmt(4) = nc_wrnlmt(4)*10
          end if
        end if
      end if
      call CPU_time(blas(11))
      time_struc = time_struc + blas(11) - blas(10)
    end if
!
    solucount = 0
    exactcrtries = 0
!
    do while ((solucount.eq.0).AND.(exactcrtries.lt.MAXCRTRIES))
!
      exactcrtries = exactcrtries + 1
!
!     generate random vector dchi, gaussian dist
      do j=1,dof
        gvec(j) = normal()
      end do
!
!     calculate d**2 = gvec_trans*gvec
      dsquare = 0.0d0
      do j=1,dof
        dsquare = dsquare + gvec(j)*gvec(j)
      end do
!
!     calculate delta phi from equation (L^t * delta_phi = gvec)
      dphi(:) = 0.0
      call trisolv(cdecomp_upper,dphi,gvec,MAXUJDOF,dof)
!
!     move prerotation segment
      call mcmove_nuccr(rsi,rsf,dof,dphi,azero)
!
!     repeat steps ABOVE for new atomic coordinates
      call buildnuccrdof(rsi,rsf,dof,nuccrmat)
!
      do j=1,dof
        do k=1,dof
          pos_def_star(j,k)=nuccrmat(j,k)*UJ_params(1)*UJ_params(2)
        end do
        pos_def_star(j,j) = pos_def_star(j,j)+UJ_params(1)
      end do
      call cdecomp(pos_def_star,MAXUJDOF,dof,cdecomp_upper_star,.true.)
      call cdecomp(pos_def_star,MAXUJDOF,dof,cdecomp_lower_star,.false.)
! 
!     calculate set of gvec phis based on (gvec_star = L_star^t * delta_phi)
      do j=1,dof
        gvec_star(j) = 0.0
        do k=1,dof
          gvec_star(j) = gvec_star(j) + cdecomp_upper_star(j,k)*dphi(k)
        end do
      end do
!     cal d_star**2
      dsquare_star = 0.0d0
      do j=1,dof
        dsquare_star = dsquare_star + gvec_star(j)*gvec_star(j)
      end do
!
!     now calculate the biasing probability for the move and reverse move
      det_L = 1.0d0
      do j=1,dof
        det_L = det_L*cdecomp_lower(j,j)
      end do
      det_L_star = 1.0d0
      do j=1,dof
        det_L_star = det_L_star*cdecomp_lower_star(j,j)
      end do
!
!     restore xyz of segment post prerotation
      call mcmove_nuccr(rsi,rsf,dof,dphi,afour)
!
!     ------------- Chain Close -------------
!
!     find the displacement angles (dpfmat) needed to close chain
      call CPU_time(blas(4))
      call nuccrchainclose(rsi,rsf,dof,dphi,ccmat_nuc,solmat,oldvals_nuc,solucount,jacv)
      if (torcrmode.eq.1) then
        if (solcnt2.gt.0) then
          jacv(solucount+1:solucount+solcnt2) = jacvp(1:solcnt2)
          solmat(solucount+1:solucount+solcnt2,:) = solmatp(1:solcnt2,:)
        end if
        if (solucount+solcnt2.gt.0) then
          call picksolu_nc(rsi,rsf,dof,dphi,solmat,solucount,solcnt2,jacv,dpfmat,which,1.0,1.0)!det_L,det_L_star)
        end if
      else
        if (solucount.gt.0) then
          call picksolu_nc(rsi,rsf,dof,dphi,solmat,solucount,solcnt2,jacv,dpfmat,which,1.0,1.0)!det_L,det_L_star)
          jac_b = jacv(which)
        end if
      end if
      call CPU_time(blas(5))
      time_struc = time_struc + blas(5) - blas(4)
!
      if (solucount.eq.0) then
!       restore pre-rotation segment (initial state fully restored)
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
      end if
!
!     do not recycle for "safe" mode
      if (torcrmode.eq.1) exit

    end do ! while solucount is 0 and exactcrtries less than MAXCRTRIES
!
    nc_totcnt = nc_totcnt + exactcrtries
!
    if (torcrmode.eq.1) then
      if ((solucount+solcnt2).eq.0) then
        nc_wrncnt(1) = nc_wrncnt(1) + 1
        if (nc_wrncnt(1).eq.nc_wrnlmt(1)) then
          write(ilog,*) 'Warning. Did not find any valid solution for Dinner-Ulmschneider concerted rotation&
 & algorithm for polynucleotides. This is indicative of a bug or of bad search settings for the closure&
 & algorithm (decrease stepsize). If this problem occurs frequently, the simulation is probably producing&
 & biased or otherwise corrupted results.'
          write(ilog,*) 'This was warning #',nc_wrncnt(1),' of this type not all of which may be displayed.'
          if (10.0*nc_wrnlmt(1).gt.0.5*HUGE(nc_wrnlmt(1))) then
            nc_wrncnt(1) = 0
          else
            nc_wrnlmt(1) = nc_wrnlmt(1)*10
          end if
        end if
!       restore pre-rotation segment (initial state fully restored)
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
        nc_deadcnt = nc_deadcnt + 1
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    else
      if ((solucount.eq.0).AND.(exactcrtries.eq.MAXCRTRIES)) then
        nc_wrncnt(1) = nc_wrncnt(1) + 1
        if (nc_wrncnt(1).eq.nc_wrnlmt(1)) then
          write(ilog,*) 'Warning. Did not find a valid solution for Dinner-Ulmschneider concerted rotation&
 & algorithm for polynucleotides after ',MAXCRTRIES,' pre-rotation attempts.'
          write(ilog,*) 'This is indicative of a bug, of too large pre-rotation stepsizes or of bad search &
 &settings for the closure algorithm (decrease stepsize). If this problem occurs frequently, the simulation is&
 & probably producing biased or otherwise corrupted results.'
          write(ilog,*) 'This was warning #',nc_wrncnt(1),' of this type not all of which may be displayed.'
          if (10.0*nc_wrnlmt(1).gt.0.5*HUGE(nc_wrnlmt(1))) then
            nc_wrncnt(1) = 0
          else
            nc_wrnlmt(1) = nc_wrnlmt(1)*10
          end if
        end if
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
        nc_deadcnt = nc_deadcnt + 1
!        mvcnt%nnuc = mvcnt%nnuc + 1
!        j = 0
!        if (random().gt.(5./6.)) j = 1
!        call nuc_bailout(rsf-j,istep)
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
!   and test for proper closure (testing)
    if ((abs(x(nuci(rsf,6))-blas(1))+abs(y(nuci(rsf,6))-blas(2))+abs(z(nuci(rsf,6))-blas(3)))&
 &        .gt.1.0e-4) then
!      write(*,*) (abs(x(nuci(rsf,6))-blas(1))+abs(y(nuci(rsf,6))-blas(2))+abs(z(nuci(rsf,6))-blas(3)))
      write(ilog,*) 'Inexact closure during chain closure algorithm.&
 & This is potentially a bug. Please record following information.'
      write(ilog,*) 'Step #',istep,' and final res. ',rsf
      call fexit()
    end if
!
!   evaluate pre-rotation biases
    if (torcrmode.eq.1) then
      if (which.le.solucount) then
        prob_move = (det_L_star/det_L)*exp(-(dsquare_star-dsquare))
      else
        prob_move = 1.0
      end if
!      prob_back = exp(-dsquare_star)
    else
      prob_move = (det_L_star/det_L)*exp(-(dsquare_star-dsquare))*(jac_b/jac_a)
!      prob_back = det_L_star*exp(-dsquare_star)
    end if
!
!   compute posterior short-range (steric) energy
    call CPU_time(t1)
!
!   first lets compute the posterior short-range contributions
    if (use_POLY.EQV..true.) then
      call update_rigid(molofrs(rsi))
    end if
    call ujcr_energy_short(rsi,rsf,eveca,use_cutoffs,azero)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rsi),eveca)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rsi),eveca,azero)
    if (use_DSSP.EQV..true.) call en_dssp_gl(eveca)
    if (use_DREST.EQV..true.) call edrest(eveca)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
!   now make use of pre-screening (if so desired)
    if (use_stericscreen.EQV..true.) then
      difft = 0.0
      do i=1,MAXENERGYTERMS
        difft = difft + eveca(i)
      end do
      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
        call mcmove_nuccr(rsi,rsf,dof,dpfmat,atwo)
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   restore the coordinates to the prior state
    call mcmove_nuccr(rsi,rsf,dof,dpfmat,atwo)
!
!   now compute the prior energies, the short-range part will already
!   yield the full SAV-information (if relevant)
    call CPU_time(t1)
    call ujcr_energy_short(rsi,rsf,evecp,use_cutoffs,aone)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rsi),evecp)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rsi),evecp,aone)
    if (use_DSSP.EQV..true.) call en_dssp_gl(evecp)
    if (use_DREST.EQV..true.) call edrest(evecp)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
!   the long-range part has a different cutoff and is computed separately
!   also, it will use the old atsavs, but compute only those pair terms that changed
!   at the end the fxn updates the atsavs to the sampled conformation (generated in the next step)
    call ujcr_energy_long(rsi,rsf,evecp,use_cutoffs,azero)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   now re-apply the sampling move
    call mcmove_nuccr(rsi,rsf,dof,dpfmat,aone)
!
!   time to get the long-range posterior energies (with updated atsav's)
    call CPU_time(t1)
    call ujcr_energy_long(rsi,rsf,eveca,use_cutoffs,aone)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   get energy difference and determine acceptance
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
!   update or restore
    if (moveok.EQV..true.) then
      acc%nnuccr = acc%nnuccr + 1
      do i=rsi,rsf
        acc%cr(i) = acc%cr(i) + 1
      end do
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = esterms(i) + evecd(i)
      end do
      if (use_POLY.EQV..false.) call internal_bound(molofrs(rsi))
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
    else
      call mcmove_nuccr(rsi,rsf,dof,dpfmat,atwo)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rsi))
    end if
!  else if (alpha_move.EQV..true.) then
!c
!c       here we will sample a few residues in a biased way to minimize lever-arm effects
!c       effectively, this is a more local move, which is very useful for collapsed chains
!c
!    moveok = .false.
!c
!c       get a suitable molecule
!    k = 1!random()*nmolcr + 1
!    imol = 1!molcr(k)
!c       get a set of consecutive residues
!    if ((rsmol(imol,2)-rsmol(imol,1)-2).eq.1) then
!      rsi = rsmol(imol,1) + 1
!      rsf = rsmol(imol,2) - 1
!    else
!      rsi = floor(random()*(rsmol(imol,2)-rsmol(imol,1)-2))
! &                            + rsmol(imol,1) + 1
!      rsf = rsi + 1
!    end if
!    do i=1,MAXCRDOF
!      dfy(i) = 0.0d0
!    end do
!c
!c       get the forward move and the forward biasing probability
!    do i=rsi,rsf
!      dfy(2*(i-rsi)+1) = (random()*90.0 - 90.0 - phi(i))/RADIAN
!      dfy(2*(i-rsi)+2) = (random()*90.0 - 90.0 - psi(i))/RADIAN
!    end do
!    bias_fwd = 16**2
!c
!c       backup reference coordinates and apply the perturbation given by dfy
!c        write(*,*) '0',x(33),x(46)
!    call mcmove_cr(rsi,rsf+1,dfy,0)
!c        write(*,*) '1',x(33),x(46)!phi(rsi+1),psi(rsi+1)
!c
!c       compute posterior short-range (steric) energy
!    call CPU_time(t1)
!c
!c       first let's compute the posterior short-range contributions
!   call cr_energy_short(rsi,rsf,eveca,use_cutoffs,0)
!c       strictly global energy fxns are handled outside of the residue-based fxns
!    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rsi),eveca)
!c
!c       now make use of pre-screening (if so desired)
!    if (use_stericscreen.EQV..true.) then
!      difft = 0.0
!      do i=1,MAXENERGYTERMS
!        difft = difft + eveca(i)
!      end do
!      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
!        call CPU_time(t2)
!        time_energy = time_energy + t2 - t1
!        call mcmove_cr(rsi,rsf+1,dfy,2)
!        return
!      end if
!    end if
!c
!    call CPU_time(t2)
!    time_energy = time_energy + t2 - t1
!c
!c       restore the coordinates to the prior state
!    call mcmove_cr(rsi,rsf+1,dfy,2)
!c        write(*,*) '2',x(33),x(46)!phi(rsi+1),psi(rsi+1)
!c
!c       now compute the prior energies, the short-range part will already
!c       yield the full SAV-information (if relevant)
!    call CPU_time(t1)
!    call cr_energy_short(rsi,rsf,evecp,use_cutoffs,1)
!c       strictly global energy fxns are handled outside of the residue-based fxns
!    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rsi),evecp)
!c
!c       the long-range part has a different cutoff and is computed separately
!c       also, it will use the old atsav's, but compute only those pair terms that changed
!c       at the end the fxn updates the atsav's to the sampled conformation (generated in the next step)
!    call cr_energy_long(rsi,rsf,evecp,use_cutoffs,0)
!    call CPU_time(t2)
!    time_energy = time_energy + t2 - t1
!c
!c       now re-apply the sampling move
!    call mcmove_cr(rsi,rsf+1,dfy,1)
!c        write(*,*) '3',x(33),x(46)!phi(rsi+1),psi(rsi+1)
!c
!c       time to get the long-range posterior energies (with updated atsav's)
!    call CPU_time(t1)
!    call cr_energy_long(rsi,rsf,eveca,use_cutoffs,1)
!    call CPU_time(t2)
!    time_energy = time_energy + t2 - t1
!c
!c       get the total energy difference
!    diff = 0.0
!    do i=1,MAXENERGYTERMS
!      evecd(i) = eveca(i) - evecp(i)
!      diff = diff + evecd(i)
!    end do
!c
!c       get the backward biasing probability
!    bias_bwd = 1.0
!c
!c       finally do the metropolis
!    bt = (bias_bwd/bias_fwd)
!    expterm = bt*exp(-diff*invtemp)
!    rcheck = random()
!    if(expterm .gt. rcheck) then
!      moveok = .true.
!    end if
!c
!c       update or restore
!    if (moveok.EQV..true.) then
!      ncraccept = ncraccept + 1
!      do i=rsi,rsf
!        fyresacc(i) = fyresacc(i) + 1
!      end do
!      esave = esave + diff
!      do i=1,MAXENERGYTERMS
!        esterms(i) = esterms(i) + evecd(i)
!      end do
!      call internal_bound(molofrs(rsi))
!    else
!      call mcmove_cr(rsi,rsf+1,dfy,2)
!      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
!    end if
!c        write(*,*) '4',x(33),x(46),moveok!phi(rsi+1),psi(rsi+1),moveok
!c
!c
!c
!c
!c
  else if (fyc_move.EQV..true.) then
!
!   here we will sample a single set of backbone angles (pivot part)
!   and the set of chi-angles for that particular residue simultaneously
!
    moveok = .false.
!
!   randomly pick a residue
    rcheck = random()
    call binary_search(fylst%nr,fylst%wt(1:fylst%nr),rcheck,rs)
    rs = fylst%idx(min(fylst%nr,rs+1))
    do i=1,nchi(rs)
      cur_chiflag(i) = .true.
    end do
!
    if (align_NC.eq.3) then
      if ((rsmol(molofrs(rs),2)-rs).gt.(rs-rsmol(molofrs(rs),1)))&
 & then
        ct = .true.
      else if (((rsmol(molofrs(rs),2)-rs)&
 &.eq.(rs-rsmol(molofrs(rs),1))).AND.(random().gt.0.5)) then
        ct = .true.
      end if
    else if (align_NC.eq.4) then
      bt = 1.0*(rsmol(molofrs(rs),2)-rs+1)/&
 &         (rsmol(molofrs(rs),2)-rsmol(molofrs(rs),1)+2)
      if (random().le.bt) ct = .true.
    end if
!
!   backup reference coordinates, perturb coord.s and re-assign grid-stuff
    call mcmove_fyc(rs,azero,ct)
!
!   compute posterior short-range (steric) energy
    call CPU_time(t1)
!
!   first let's compute the posterior short-range contributions
    call pivot_energy_short(rs,eveca,use_cutoffs,azero,ct)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rs),eveca)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),eveca,azero)
    if (use_DSSP.EQV..true.) call en_dssp_gl(eveca)
    if (use_DREST.EQV..true.) call edrest(eveca)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
!   now make use of pre-screening (if so desired)
    if (use_stericscreen.EQV..true.) then
      difft = 0.0
      do i=1,MAXENERGYTERMS
        difft = difft + eveca(i)
      end do
      if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
        call CPU_time(t2)
        time_energy = time_energy + t2 - t1
        call mcmove_fyc(rs,atwo,ct)
        if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
        call mc_acceptance_bailout(mc_acc_crit,istep)
        return
      end if
    end if
!
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   restore the coordinates to the prior state
    call mcmove_fyc(rs,atwo,ct)
!
!   now compute the prior energies, the short-range part will already
!   yield the full SAV-information (if relevant)
    call CPU_time(t1)
    call pivot_energy_short(rs,evecp,use_cutoffs,aone,ct)
!   strictly global energy fxns are handled outside of the residue-based fxns
    if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rs),evecp)
    if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),evecp,aone)
    if (use_DSSP.EQV..true.) call en_dssp_gl(evecp)
    if (use_DREST.EQV..true.) call edrest(evecp)
    if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
!   the long-range part has a different cutoff and is computed separately
!   also, it will use the old atsav's, but compute only those pair terms that changed
!   at the end the fxn updates the atsav's to the sampled conformation (generated in the next step)
    call pivot_energy_long(rs,evecp,use_cutoffs,azero,ct)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   now re-apply the sampling move
    call mcmove_fyc(rs,aone,ct)
!
!   time to get the long-range posterior energies (with updated atsav's)
    call CPU_time(t1)
    call pivot_energy_long(rs,eveca,use_cutoffs,aone,ct)
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   get energy difference and determine acceptance
    call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
!   update or restore
    if (moveok.EQV..true.) then
      acc%naccept = acc%naccept + 1
      acc%nfyc = acc%nfyc + 1
      acc%fy(rs) = acc%fy(rs) + 1
      if (nchi(rs).gt.0) acc%chi(rs) = acc%chi(rs) + 1
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = esterms(i) + evecd(i)
      end do
      if (use_POLY.EQV..false.) call internal_bound(molofrs(rs))
      if (use_EMICRO.EQV..true.) then
        call em_transfer_m()
        esave = esave + evecp(21) - (esterms(21) - evecd(21))
        esterms(21) = eveca(21)
      end if
    else
      call mcmove_fyc(rs,atwo,ct)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
    end if
!
!
!
!
  else if (lct_move.EQV..true.) then
!
!   here we will sample only a single LCT (uniformly in LCT space)
!
    moveok = .false.
!
!   randomly pick an LCT
    rs = floor(random()*ntorlcs)+1
!
!   backup reference coordinates, perturb coord.s and re-assign grid-stuff
    call mcmove_lct(rs,azero)
!
!   compute full posterior energy
    call CPU_time(t1)
!
!   everything is global in LCT moves (which is a waste if only LCT-potential is
!   used OR if a significant amount of entries in separatrix are zero)
    diff = energy3(eveca,sayyes) - esave
!
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
!
!   evaluate metropolis
    if(diff .le. 0.0d0) then
      moveok = .true.
    else
      expterm = exp(-diff*invtemp)
      rcheck = random()
      if(expterm .gt. rcheck) then
        moveok = .true.
      end if
    end if
    if (harappa.EQV..true.) then
      moveok = .false.
      harappa = .false.
    end if
!    write(*,*) lct(1),lct(2),moveok
!
!   update or restore
    if (moveok.EQV..true.) then
      acc%nlct = acc%nlct + 1
      esave = esave + diff
      do i=1,MAXENERGYTERMS
        esterms(i) = eveca(i)
      end do
    else
      call mcmove_lct(rs,atwo)
      if (use_IMPSOLV.EQV..true.) call init_svte(afour)
    end if
!    write(*,*) 'A:',lct(1),lct(2)
!    write(*,*) 'B:',phi(2),psi(2)
!
  else if (ph_move .EQV. .true.) then
     call CPU_time(t1)
     call protmov()
     call CPU_time(t2)
     time_ph = time_ph + t2 - t1
!
  else
     write(ilog,*) 'Fatal. Entered unsupported move type. Check with&
 & the developer.'
     nstep = nstep - 1
     call fexit()

!   WARNING: old version - currently not supported
!
!    moveok = .false.
!
!    eprior = esave
!    call makeref
!    call pickres(ires)
!
!  if (.not. use_elscp) then
!     eafter = energy()
!     diff = eafter - eprior
!  else
!     call elscp(ires,earmafter)
!     diff = earmafter - earmsave
!  end if
!
!  if(diff .le. 0.0d0) then
!     moveok = .true.
!  else
!     expterm = exp(-diff*invtemp)
!     rcheck = random()
!     if(expterm .gt. rcheck) then
!        moveok = .true.
!     end if
!  end if
!
! process results of proposed move based on accept / reject 
! criterion
!
!  if(moveok.EQV..true.) then
!     esave = eafter
!     esatt = eattract
!     naccept = naccept + 1
!  else
!     if(use_hc .and. (hccount .eq. 0)) then
!        reshc = reshc_old
!        hcnum = hcnum_old
!     end if
!     call restore(ires)
!     call getref
!  end if
!c
  end if
  
!  write(*,*) pivot_move,fyc_move,omega_move,chi_move,ph_move,pucker_move
!  write(*,*) particleinsertion_move,particledeletion_move,particleidentity_move
!  write(*,*) UJconrot_move,SJconrot_move,torcr_move,torcr_move_omega
!  write(*,*) nuc_move,nucpuck_move,nuccr_move
!  write(*,*) rigid_move,clurb_move,trans_move,rot_move,hsq_move
end
!
!
!---------------------------------------------------------------------------
!
! emat: different conditions in first and second index, however structures
!       only preserved along second index
!       contains beta*E 
! lmap: the mapping vector for possible swaps between conditions
!
subroutine RE_swap(emat,lmap)
!
  use iounit
  use system
  use mpistuff
!
  implicit none
!
  integer lmap(mpi_nodes),i,ii,jj,ik,jk,k,dummy
  RTYPE emat(mpi_nodes,mpi_nodes),diff,expterm,random
!
  do k=1,re_conditions
    lmap(k) = k
  end do
!
  do i=1,re_tryswap
    ii = floor(random()*re_conditions)+1
    if (re_nbmode.eq.1) then
!     we try to exchange with arbitrary conditions
      jj = ii
      do while (ii.eq.jj)
        jj = floor(random()*re_conditions)+1
      end do
      re_trans(ii,jj) = re_trans(ii,jj) + 1
      re_trans(jj,ii) = re_trans(jj,ii) + 1
      ik = lmap(ii)
      jk = lmap(jj)
      diff = emat(ii,jk)+emat(jj,ik)-emat(ii,ik)-emat(jj,jk)
      expterm = exp(-diff)
      if (expterm.gt.random()) then
        dummy = lmap(ii)
        lmap(ii) = lmap(jj)
        lmap(jj) = dummy
      end if
      if (diff.le.0.0) then
        re_probs(ii,jj) = re_probs(ii,jj) + 1.0
        re_probs(jj,ii) = re_probs(jj,ii) + 1.0
      else
        re_probs(ii,jj) = re_probs(ii,jj) + expterm
        re_probs(jj,ii) = re_probs(jj,ii) + expterm
      end if
    else if (re_nbmode.eq.2) then
!     we only exchange with neighboring conditions
!     note that users always supply a vector regardless of dimensionality
!     hence, the neighboring idea isn't all that well-defined for dimensions > 1
      if (ii.eq.1) then
        jj = ii + 1
      else if (ii.eq.re_conditions) then
        jj = ii - 1
      else
        if (random().gt.0.5) then
          jj = ii + 1
        else
          jj = ii - 1
        end if
      end if
      re_trans(ii,jj) = re_trans(ii,jj) + 1
      re_trans(jj,ii) = re_trans(jj,ii) + 1
      ik = lmap(ii)
      jk = lmap(jj)
      if (abs(ik-jk).gt.1) cycle
      diff = emat(ii,jk)+emat(jj,ik)-emat(ii,ik)-emat(jj,jk)
      expterm = exp(-diff*invtemp)
      if (expterm.gt.random()) then
        dummy = lmap(ii)
        lmap(ii) = lmap(jj)
        lmap(jj) = dummy
      end if
      if (diff.le.0.0) then
        re_probs(ii,jj) = re_probs(ii,jj) + 1.0
        re_probs(jj,ii) = re_probs(jj,ii) + 1.0
      else
        re_probs(ii,jj) = re_probs(ii,jj) + expterm
        re_probs(jj,ii) = re_probs(jj,ii) + expterm
      end if
    end if
  end do
!
end
!
!---------------------------------------------------------------------------------------
!
subroutine pivot_bailout(rs,istep)
!
  use iounit
  use movesets
  use energies
  use sequen
  use molecule
  use accept
  use system
  use mcsums
  use cutoffs
  use atoms
  use polypep
!
  implicit none
!
  integer rs,i,istep
  integer aone,azero,atwo,afour
  logical moveok,sayyes,sayno,ct
  RTYPE diff,evecp(MAXENERGYTERMS),eveca(MAXENERGYTERMS),difft
  RTYPE evecd(MAXENERGYTERMS)
  RTYPE random,t1,t2,bt,prob_move
!
  moveok = .false.
  azero = 0
  aone = 1
  atwo = 2
  afour = 4
  sayyes = .true.
  sayno = .false.
  evecp(:) = 0.0
  eveca(:) = 0.0
  prob_move = 1.0
!
  if (align_NC.eq.2) then
    ct = .true.
  else ! see below in individual move-sets for further modification
    ct = .false.
  end if
! 
  if (align_NC.eq.3) then
    if ((rsmol(molofrs(rs),2)-rs).gt.(rs-rsmol(molofrs(rs),1)))&
 & then
      ct = .true.
    else if (((rsmol(molofrs(rs),2)-rs)&
 &.eq.(rs-rsmol(molofrs(rs),1))).AND.(random().gt.0.5)) then
      ct = .true.
    end if
  else if (align_NC.eq.4) then
    bt = 1.0*(rsmol(molofrs(rs),2)-rs+1)/&
 &         (rsmol(molofrs(rs),2)-rsmol(molofrs(rs),1)+2)
    if (random().le.bt) ct = .true.
  end if
!
  call mcmove_pivot(rs,azero,ct)
!
!   compute posterior short-range (steric) energy
  call CPU_time(t1)
!
! first let's compute the posterior short-range contributions
  call pivot_energy_short(rs,eveca,use_cutoffs,azero,ct)
! strictly global energy fxns are handled outside of the residue-based fxns
  if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rs),eveca)
  if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),eveca,azero)
  if (use_DSSP.EQV..true.) call en_dssp_gl(eveca)
  if (use_DREST.EQV..true.) call edrest(eveca)
  if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
! now make use of pre-screening (if so desired)
  if (use_stericscreen.EQV..true.) then
    difft = 0.0
    do i=1,MAXENERGYTERMS
      difft = difft + eveca(i)
    end do
    if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
      call CPU_time(t2)
      time_energy = time_energy + t2 - t1
  call mcmove_pivot(rs,atwo,ct)
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
      call mc_acceptance_bailout(mc_acc_crit,istep)
      return
    end if
  end if
!
  call CPU_time(t2)
  time_energy = time_energy + t2 - t1
!
! restore the coordinates to the prior state
  call mcmove_pivot(rs,atwo,ct)
!    write(*,*) '2',ztor(yline(4)),psi(4)
!
! now compute the prior energies, the short-range part will already
! yield the full SAV-information (if relevant)
  call CPU_time(t1)
  call pivot_energy_short(rs,evecp,use_cutoffs,aone,ct)
! strictly global energy fxns are handled outside of the residue-based fxns
  if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rs),evecp)
  if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),evecp,aone)
  if (use_DSSP.EQV..true.) call en_dssp_gl(evecp)
  if (use_DREST.EQV..true.) call edrest(evecp)
  if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
! the long-range part has a different cutoff and is computed separately
! also, it will use the old atsav's, but compute only those pair terms that changed
! at the end the fxn updates the atsav's to the sampled conformation (generated in the next step)
  call pivot_energy_long(rs,evecp,use_cutoffs,azero,ct)
  call CPU_time(t2)
  time_energy = time_energy + t2 - t1
! now re-apply the sampling move
  call mcmove_pivot(rs,aone,ct)
!
! time to get the long-range posterior energies (with updated atsav's)
  call CPU_time(t1)
  call pivot_energy_long(rs,eveca,use_cutoffs,aone,ct)
  call CPU_time(t2)
  time_energy = time_energy + t2 - t1
!
! get energy difference and determine acceptance
  call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
! update or restore
  if (moveok.EQV..true.) then
    acc%nfy = acc%nfy + 1
    acc%fy(rs) = acc%fy(rs) + 1
    esave = esave + diff
    do i=1,MAXENERGYTERMS
      esterms(i) = esterms(i) + evecd(i)
    end do
    if (use_POLY.EQV..false.) call internal_bound(molofrs(rs))
    if (use_EMICRO.EQV..true.) then
      call em_transfer_m()
      esave = esave + evecp(21) - (esterms(21) - evecd(21))
      esterms(21) = eveca(21)
    end if
!      
  else
    call mcmove_pivot(rs,atwo,ct)
    if (use_IMPSOLV.EQV..true.) call init_svte(afour)
    if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
  end if
!
end
!
!------------------------------------------------------------------------------------
!
subroutine nuc_bailout(rs,istep)
!
  use iounit
  use movesets
  use energies
  use sequen
  use molecule
  use accept
  use system
  use mcsums
  use cutoffs
  use fyoc
  use polypep
!
  implicit none
!
  integer rs,i,istep
  integer aone,azero,atwo,afour
  logical moveok,sayyes,sayno,ct,chidone
  RTYPE diff,evecp(MAXENERGYTERMS),eveca(MAXENERGYTERMS),difft
  RTYPE evecd(MAXENERGYTERMS)
  RTYPE random,t1,t2,bt,prob_move
!
  moveok = .false.
  azero = 0
  aone = 1
  atwo = 2
  afour = 4
  sayyes = .true.
  sayno = .false.
  evecp(:) = 0.0
  eveca(:) = 0.0
  prob_move = 1.0
!
! randomly pick the individual degrees of freedom
  chidone = .false.
  do while (chidone.EQV..false.)
    if (nnucs(rs).le.nrnucim) then
      do i=1,nnucs(rs)
        cur_nucflag(i) = .true.
      end do
      chidone = .true.
    else
      do i=1,nnucs(rs)
        if (random().lt.((1.0*nrnucim)/nnucs(rs))) then
          cur_nucflag(i) = .true.
          chidone = .true.
        end if
      end do
    end if
  end do
!
! handle alignment
  if (align_NC.eq.3) then
    if ((rsmol(molofrs(rs),2)-rs).gt.(rs-rsmol(molofrs(rs),1)))&
 & then
      ct = .true.
    else if (((rsmol(molofrs(rs),2)-rs)&
 &.eq.(rs-rsmol(molofrs(rs),1))).AND.(random().gt.0.5)) then
      ct = .true.
    end if
  else if (align_NC.eq.4) then
    bt = 1.0*(rsmol(molofrs(rs),2)-rs+1)/&
 &         (rsmol(molofrs(rs),2)-rsmol(molofrs(rs),1)+2)
    if (random().le.bt) ct = .true.
  end if
!
! backup reference coordinates, perturb coord.s and re-assign grid-stuff
  call mcmove_nuc(rs,azero,ct)
!
! compute posterior short-range (steric) energy
  call CPU_time(t1)
!
! first let's compute the posterior short-range contributions
  call pivot_energy_short(rs,eveca,use_cutoffs,azero,ct)
! strictly global energy fxns are handled outside of the residue-based fxns
  if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),eveca,azero)
  if (use_DREST.EQV..true.) call edrest(eveca)
  if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
! now make use of pre-screening (if so desired)
  if (use_stericscreen.EQV..true.) then
    difft = 0.0
    do i=1,MAXENERGYTERMS
      difft = difft + eveca(i)
    end do
    if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
      call CPU_time(t2)
      time_energy = time_energy + t2 - t1
      call mcmove_nuc(rs,atwo,ct)
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
      call mc_acceptance_bailout(mc_acc_crit,istep)
      return
    end if
  end if
!
  call CPU_time(t2)
  time_energy = time_energy + t2 - t1
!
! restore the coordinates to the prior state
  call mcmove_nuc(rs,atwo,ct)
!
! now compute the prior energies, the short-range part will already
! yield the full SAV-information (if relevant)
  call CPU_time(t1)
  call pivot_energy_short(rs,evecp,use_cutoffs,aone,ct)
! strictly global energy fxns are handled outside of the residue-based fxns
  if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),evecp,aone)
  if (use_DREST.EQV..true.) call edrest(evecp)
  if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
! the long-range part has a different cutoff and is computed separately
! also, it will use the old atsav's, but compute only those pair terms that changed
! at the end the fxn updates the atsav's to the sampled conformation (generated in the next step)
  call pivot_energy_long(rs,evecp,use_cutoffs,azero,ct)
  call CPU_time(t2)
  time_energy = time_energy + t2 - t1
!
! now re-apply the sampling move
  call mcmove_nuc(rs,aone,ct)
!
! time to get the long-range posterior energies (with updated atsav's)
  call CPU_time(t1)
  call pivot_energy_long(rs,eveca,use_cutoffs,aone,ct)
  call CPU_time(t2)
  time_energy = time_energy + t2 - t1
!
! get energy difference and determine acceptance
  call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
! update or restore
  if (moveok.EQV..true.) then
    acc%nnuc = acc%nnuc + 1
    acc%nuc(rs) = acc%nuc(rs) + 1
    esave = esave + diff
    do i=1,MAXENERGYTERMS
      esterms(i) = esterms(i) + evecd(i)
    end do
    if (use_POLY.EQV..false.) call internal_bound(molofrs(rs))
    if (use_EMICRO.EQV..true.) then
      call em_transfer_m()
      esave = esave + evecp(21) - (esterms(21) - evecd(21))
      esterms(21) = eveca(21)
    end if
!      
  else
    call mcmove_nuc(rs,atwo,ct)
    if (use_IMPSOLV.EQV..true.) call init_svte(afour)
    if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
  end if
!
end
!
!------------------------------------------------------------------------------
!
subroutine omega_bailout(rs,istep)
!
  use iounit
  use movesets
  use energies
  use sequen
  use molecule
  use accept
  use system
  use mcsums
  use cutoffs
  use polypep
!
  implicit none
!
  integer rs,i,istep
  integer aone,azero,atwo,afour
  logical moveok,sayyes,sayno,ct
  RTYPE diff,evecp(MAXENERGYTERMS),eveca(MAXENERGYTERMS),difft
  RTYPE evecd(MAXENERGYTERMS)
  RTYPE random,t1,t2,bt,prob_move
!
  moveok = .false.
  azero = 0
  aone = 1
  atwo = 2
  afour = 4
  sayyes = .true.
  sayno = .false.
  evecp(:) = 0.0
  eveca(:) = 0.0
  prob_move = 1.0
!
!
  if (align_NC.eq.3) then
    if ((rsmol(molofrs(rs),2)-rs).gt.(rs-rsmol(molofrs(rs),1)))&
 & then
      ct = .true.
    else if (((rsmol(molofrs(rs),2)-rs)&
 &.eq.(rs-rsmol(molofrs(rs),1))).AND.(random().gt.0.5)) then
      ct = .true.
    end if
  else if (align_NC.eq.4) then
    bt = 1.0*(rsmol(molofrs(rs),2)-rs+1)/&
 &       (rsmol(molofrs(rs),2)-rsmol(molofrs(rs),1)+2)
    if (random().le.bt) ct = .true.
  end if
!
! backup reference coordinates, perturb coord.s and re-assign grid-stuff
  call mcmove_omega(rs,azero,ct)
!
! compute posterior short-range (steric) energy
  call CPU_time(t1)
!
! first let's compute the posterior short-range contributions
  call pivot_energy_short(rs,eveca,use_cutoffs,azero,ct)
! strictly global energy fxns are handled outside of the residue-based fxns
  if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rs),eveca)
  if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),eveca,azero)
  if (use_DSSP.EQV..true.) call en_dssp_gl(eveca)
  if (use_DREST.EQV..true.) call edrest(eveca)
  if (use_EMICRO.EQV..true.) call en_emicro_gl(eveca,atwo)
!
! now make use of pre-screening (if so desired)
  if (use_stericscreen.EQV..true.) then
    difft = 0.0
    do i=1,MAXENERGYTERMS
      difft = difft + eveca(i)
    end do
    if ((difft.gt.(2*esave)).AND.(difft.gt.screenbarrier)) then
      call CPU_time(t2)
      time_energy = time_energy + t2 - t1
      call mcmove_omega(rs,atwo,ct)
      if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
      call mc_acceptance_bailout(mc_acc_crit,istep)
      return
    end if
  end if
!
  call CPU_time(t2)
  time_energy = time_energy + t2 - t1
!
! restore the coordinates to the prior state
  call mcmove_omega(rs,atwo,ct)
!
! now compute the prior energies, the short-range part will already
! yield the full SAV-information (if relevant)
  call CPU_time(t1)
  call pivot_energy_short(rs,evecp,use_cutoffs,aone,ct)
! strictly global energy fxns are handled outside of the residue-based fxns
  if (use_ZSEC.EQV..true.) call en_zsec_gl(molofrs(rs),evecp)
  if (use_POLY.EQV..true.) call en_poly_gl(molofrs(rs),evecp,aone)
  if (use_DSSP.EQV..true.) call en_dssp_gl(evecp)
  if (use_DREST.EQV..true.) call edrest(evecp)
  if (use_EMICRO.EQV..true.) call en_emicro_gl(evecp,aone)
!
! the long-range part has a different cutoff and is computed separately
! also, it will use the old atsav's, but compute only those pair terms that changed
! at the end the fxn updates the atsav's to the sampled conformation (generated in the next step)
  call pivot_energy_long(rs,evecp,use_cutoffs,azero,ct)
  call CPU_time(t2)
  time_energy = time_energy + t2 - t1
!
! now re-apply the sampling move
  call mcmove_omega(rs,aone,ct)
!
! time to get the long-range posterior energies (with updated atsav's)
  call CPU_time(t1)
  call pivot_energy_long(rs,eveca,use_cutoffs,aone,ct)
  call CPU_time(t2)
  time_energy = time_energy + t2 - t1
!
! get energy difference and determine acceptance
  call mc_acceptance(diff,esave,evecp,eveca,evecd,mc_acc_crit,moveok,istep,prob_move)
!
! update or restore
  if (moveok.EQV..true.) then
    acc%nomega = acc%nomega + 1
    acc%omega(rs) = acc%omega(rs) + 1
    esave = esave + diff
    do i=1,MAXENERGYTERMS
      esterms(i) = esterms(i) + evecd(i)
    end do
    call internal_bound(molofrs(rs))
    if (use_EMICRO.EQV..true.) then
      call em_transfer_m()
      esave = esave + evecp(21) - (esterms(21) - evecd(21))
      esterms(21) = eveca(21)
    end if
  else
    call mcmove_omega(rs,atwo,ct)
    if (use_IMPSOLV.EQV..true.) call init_svte(afour)
    if (use_POLY.EQV..true.) call restore_poly(molofrs(rs))
  end if
end
!
!------------------------------------------------------------------------------------
!
! a generic MC acceptance routine
! acctyp = 1: Metropolis
! acctyp = 2: Fermi
! acctyp = 3: WL (requires auxiliary constructs in module wl)
!
subroutine mc_acceptance(diff,ecurr,evecp,eveca,evecd,acctyp,moveok,istep,biasterm)!Martin : Here is the actual acceptance 
    !subroutine
!
  use iounit
  use energies
  use wl
  use molecule
  use system
  use torsn
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer istep,idx_new,idx_old,idx_new2,idx_old2,acctyp,kk,atwo
  RTYPE biasterm,ecurr,evecp(MAXENERGYTERMS),eveca(MAXENERGYTERMS),evecd(MAXENERGYTERMS),diff,expterm,random
  RTYPE els(MAXENERGYTERMS)
  logical moveok,inrange
!
  atwo = 2
  evecd(:) = eveca(:)-evecp(:)

  if (do_accelsim.EQV..true.) then
    els(:) = esterms(:) + evecd(:) - hmjam%boosts(:,1)
    call els_manage_justE(els,atwo)
    eveca(:) = eveca(:) + hmjam%boosts(:,2)
    evecp(:) = evecp(:) + hmjam%boosts(:,1)
    evecd(:) = eveca(:)-evecp(:)
  end if
!
  diff = sum(evecd)

  
  
  
! Metropolis-Rosenbluth-Hastings
  if (acctyp.eq.1) then
!   depending on move, bias term may be 1.0 and diff may be 0.0 even for nonzero Ham
    expterm = biasterm*exp(-diff*invtemp)
    if (expterm.gt.1.000001) then ! saves a random while avoiding issues with exact reproducibility (see above)
      moveok = .true.
    else
      if (expterm.gt.random()) then
        moveok = .true.
      else
        moveok = .false.
      end if
    end if
! Fermi
  else if (acctyp.eq.2) then
    expterm = biasterm/(exp(diff*invtemp) + biasterm)
    if (expterm.gt.random()) then
      moveok = .true.
    else
      moveok = .false.
    end if
! Metropolis bail-out for Wang-Landau
  else if ((acctyp.eq.3).AND.(istep.le.nequil)) then
    expterm = biasterm*exp(-diff*invtemp)
    if (expterm.ge.1.0) then
      moveok = .true.
    else
      if (expterm.gt.random()) then
        moveok = .true.
      else
        moveok = .false.
      end if
    end if
! Wang-Landau
  else if (acctyp.eq.3) then
!   If you are here, we are now only doing WL -> increment counter
    wld%stepnum = wld%stepnum + 1
    call wl_get_binpos(inrange,ecurr,diff,idx_old,idx_new,idx_old2,idx_new2)
!
    if (inrange.EQV..false.) then
!     this includes cases where allocation size is violated and where fixed interval boundaries are violated
      moveok = .false.
      idx_new = idx_old
      if (wld%dimensionality.eq.2) idx_new2 = idx_old2
    else
!     use WL criterion if in-range, otherwise reject (technically incorrect!)
      if (wld%dimensionality.eq.1) then
        if (wld%wl_mode.eq.1) then
          expterm = biasterm*exp(wld%g(idx_old) - wld%g(idx_new))
        else if ((wld%wl_mode.eq.2).OR.(wld%wl_mode.eq.3)) then
          expterm = biasterm*exp(-diff*invtemp + (wld%g(idx_old) - wld%g(idx_new)))
        end if
      else if (wld%dimensionality.eq.2) then
        if ((wld%wl_mode.eq.2).OR.(wld%wl_mode.eq.3)) then
          expterm = biasterm*exp(-diff*invtemp + (wld%g2d(idx_old,idx_old2) - wld%g2d(idx_new,idx_new2)))
        end if
      end if
      if (expterm.ge.1.0) then 
        moveok = .true.
        wld%accepted = wld%accepted + 1
      else
        if (expterm.gt.random()) then
          moveok = .true.
          wld%accepted = wld%accepted + 1
        else
          moveok = .false.
          idx_new = idx_old  !stay where you are
          if (wld%dimensionality.eq.2) idx_new2 = idx_old2
        end if
      end if
!     extend current range counters
      if (wld%dimensionality.eq.1) then
        kk = max(1,min(wld%nbins,idx_new))
        wld%maxb = max(wld%maxb,kk)
        wld%minb = min(wld%minb,kk)
      else if (wld%dimensionality.eq.2) then
        kk = max(1,min(wld%nbins2d(1),idx_new))
        wld%maxb2d(1) = max(wld%maxb2d(1),kk)
        wld%minb2d(1) = min(wld%minb2d(1),kk)
        kk = max(1,min(wld%nbins2d(2),idx_new2))
        wld%maxb2d(2) = max(wld%maxb2d(2),kk)
        wld%minb2d(2) = min(wld%minb2d(2),kk)
      end if
    end if
!   udpate histograms
    call wl_update_histos(istep,idx_new,idx_new2)
  else
    write(ilog,*) 'Fatal. Encountered unsupported acceptance criterion flag in &
 &mc_acceptance(...). This is a bug.'
    call fexit()
  end if
!
! may integrate more default operations into this block
  if ((moveok.EQV..true.).AND.(do_accelsim.EQV..true.)) then
    hmjam%boosts(:,1) = hmjam%boosts(:,2)
  end if
!
end
!
!------------------------------------------------------------------------------------
!
! for WL runs, premature exits in the various branches of mcmove need histograms
! to be updated regardless
!
subroutine mc_acceptance_bailout(acctyp,istep)
!
  use wl
  use energies
  use system
!
  implicit none
!
  RTYPE diff
  integer istep,idx_old,idx_old2,idx_new,idx_new2,acctyp
  logical inrange
!
  if (acctyp.eq.3) then
    if (istep.le.nequil) return ! irrelevant
!
!   for actual WL increment counter and update bins
    wld%stepnum = wld%stepnum + 1
    diff = 0.0
    call wl_get_binpos(inrange,esave,diff,idx_old,idx_new,idx_old2,idx_new2)
!
    idx_new = idx_old
    if (wld%dimensionality.eq.2) idx_new2 = idx_old2
!   udpate histograms
    call wl_update_histos(istep,idx_new,idx_new2)
  end if
!
end
!
!------------------------------------------------------------------------------------
!
