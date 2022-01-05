!--------------------------------------------------------------------------!
!    This file is part of martin's  personal hacked version of CAMPARI.    !
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
! AUTHORSHIP INFO: Added standard CAMPARI by Martin                        !
!--------------------------------------------------------------------------!
!                                                                          !
! MAIN AUTHOR:   Martin Fossat                                             !
!                                                                          !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"

module martin_own
!
  RTYPE, ALLOCATABLE:: en_fos_memory(:,:)! This array has the dimension : (nbfr,res,nsolgrp)
  RTYPE, ALLOCATABLE:: sav_memory(:,:,:)
  RTYPE, ALLOCATABLE :: lj_rep_memory(:,:),lj_att_memory(:,:)
  integer, ALLOCATABLE :: memory_1_4(:,:)
  logical done_dets
  
  public test_combi ,fos_Tdep_duplicate   
  
  contains
  
  function test_combi(array1,array2,size) result(res)
    ! Martin : this is a simple function that I may need in several parts of the program
    ! It checks that there is s single equivalent for each element of a list in another list, and returns the mapping
    implicit none
    integer, intent(in) :: size
    integer, dimension(size), intent(in) :: array1,array2
    integer i,j
    integer, dimension(size):: res,local_mem

    do i=1,size
        local_mem(i)=0
        res(i)=0
    end do 
    
    do i=1,size
        do j=1,size
            if (((array1(i).eq.array2(j)).and.(res(i).eq.0)).and.(local_mem(j).ne.1)) then
                res(i)=j
                local_mem(j)=1
            end if
        end do
    end do
    
    return 
    
end function test_combi

function fos_Tdep_duplicate(tdvs) result(fos_Tdep)
!
  use system
  use fos
!
  implicit none
!


  RTYPE tdvs(3),fos_Tdep,Trat
!
  if (use_Tdepfosvals.EQV..true.) then
    Trat = kelvin/fos_Tdepref
    fos_Tdep = Trat*(tdvs(1)-tdvs(2)) + tdvs(2) + 0.01*(kelvin*(1.0-log(Trat))-fos_Tdepref)
  else
    fos_Tdep = tdvs(1)
  end if

  return
  
end function fos_Tdep_duplicate

subroutine  align_and_extend_arrays1(limit1,limit2)
    use atoms 
    use inter
    use sequen
    use params

    logical test,reverse
    integer limit1,limit2,temp_int

    integer rs,i,j,k,imol
    integer mem(2)
    integer temp_var1(2)
    RTYPE temp_var2
    integer, allocatable:: temp_var3(:)
    
    
! The next part is the alignement of the iaa and fudges arrays for both limits. It also adds non present contact to the shorter 
! list if the number of entries is not the same
   ! For the HIS case, I can do a three way check (1vs2,1vs3,2vs1 )
    do rs=1,nseq
        do i=1,nrsintra(rs)
            test=.true.
            ! Test if the combination is the same (i.e. 14 15 is the same as 15 14)  
            mem=test_combi(iaa_limits(rs,limit1)%atin(i,:),iaa_limits(rs,limit2)%atin(i,:),2)
            do k=1,2
                if (mem(k).eq.0) then 
                    test=.false.
                end if 
            end do 
            !If not the same, check the rest of the array to see if it is somewhere else 
            if (test.eqv..false.) then 
                do j=1,nrsintra(rs)-i
                    test=.true.
                    mem=test_combi(iaa_limits(rs,limit1)%atin(i+j,:),iaa_limits(rs,limit2)%atin(i,:),2)
                    do k=1,2
                        if (mem(k).eq.0) then 
                            test=.false.
                        end if 
                    end do 
                    ! If it is somewhere else, switch values to end up with the same order
                    if (test.eqv..true.) then 
                        temp_var1(:)=iaa_limits(rs,limit1)%atin(i,:)
                        iaa_limits(rs,limit1)%atin(i,:)=iaa_limits(rs,limit1)%atin(i+j,:)
                        iaa_limits(rs,limit1)%atin(i+j,:)=temp_var1(:)

                        temp_var2=fudge_limits(rs,limit1)%rsin_ljs(i)
                        fudge_limits(rs,limit1)%rsin_ljs(i)=fudge_limits(rs,limit1)%rsin_ljs(i+j)
                        fudge_limits(rs,limit1)%rsin_ljs(i+j)=temp_var2

                        temp_var2=fudge_limits(rs,limit1)%rsin_lje(i)
                        fudge_limits(rs,limit1)%rsin_lje(i)=fudge_limits(rs,limit1)%rsin_lje(i+j)
                        fudge_limits(rs,limit1)%rsin_lje(i+j)=temp_var2
                        exit
                    end if 
                end do  
                ! If it is not somewhere else (that is you went through the entire array without finding it)
                if ((test.eqv..false.)) then ! If you didn't find a match
                    if (temp_nrsin(rs,limit1).le.temp_nrsin(rs,limit2)) then 
                    ! Not a problem since the fudge is way over allocated
                        do j=1,nrsintra(rs)-i
                            iaa_limits(rs,limit1)%atin(nrsintra(rs)+1-j,:)=iaa_limits(rs,limit1)%atin(nrsintra(rs)-j,:)! Shifts the 
                            !array by one
                            fudge_limits(rs,limit1)%rsin_ljs(nrsintra(rs)+1-j)=fudge_limits(rs,limit1)%rsin_ljs(nrsintra(rs)-j)
                            fudge_limits(rs,limit1)%rsin_lje(nrsintra(rs)+1-j)=fudge_limits(rs,limit1)%rsin_lje(nrsintra(rs)-j)
                        end do 

                        iaa_limits(rs,limit1)%atin(i,:)=iaa_limits(rs,limit2)%atin(i,:)
                        fudge_limits(rs,limit1)%rsin_ljs(i)=fudge_limits(rs,limit2)%rsin_ljs(i)
                        fudge_limits(rs,limit1)%rsin_lje(i)=fudge_limits(rs,limit2)%rsin_lje(i)

                    else 
                        do j=1,nrsintra(rs)-i
                            iaa_limits(rs,limit2)%atin(nrsintra(rs)+1-j,:)=iaa_limits(rs,limit2)%atin(nrsintra(rs)-j,:)
                            fudge_limits(rs,limit2)%rsin_ljs(nrsintra(rs)+1-j)=fudge_limits(rs,limit2)%rsin_ljs(nrsintra(rs)-j)
                            fudge_limits(rs,limit2)%rsin_lje(nrsintra(rs)+1-j)=fudge_limits(rs,limit2)%rsin_lje(nrsintra(rs)-j)
                        end do 
                        iaa_limits(rs,limit2)%atin(i,:)=iaa_limits(rs,limit1)%atin(i,:)
                        fudge_limits(rs,limit2)%rsin_ljs(i)=fudge_limits(rs,limit1)%rsin_ljs(i)
                        fudge_limits(rs,limit2)%rsin_lje(i)=fudge_limits(rs,limit1)%rsin_lje(i)
                        
                    end if ! Martin : I should not need to make any limits zero, as the softcore should take care of all of it 
                end if 
            end if 
        end do 
    end do

    do rs=1,nseq
        do i=1,nrsnb(rs)
            test=.true.
            mem=test_combi(iaa_limits(rs,limit1)%atnb(i,:),iaa_limits(rs,limit2)%atnb(i,:),2)
            do k=1,2
                if (mem(k).eq.0) then 
                    test=.false.
                end if 
            end do 

            if (test.eqv..false.) then ! If the array line is not equivalent
                do j=1,nrsnb(rs)-i
                    test=.true.
                    mem=test_combi(iaa_limits(rs,limit1)%atnb(i+j,:),iaa_limits(rs,limit2)%atnb(i,:),2)
                    do k=1,2
                        if (mem(k).eq.0) then 
                            test=.false.
                        end if 
                    end do 
                    if (test.eqv..true.) then 
                        temp_var1(:)=iaa_limits(rs,limit1)%atnb(i,:)
                        iaa_limits(rs,limit1)%atnb(i,:)=iaa_limits(rs,limit1)%atnb(i+j,:)
                        iaa_limits(rs,limit1)%atnb(i+j,:)=temp_var1(:)

                        temp_var2=fudge_limits(rs,limit1)%rsnb_ljs(i)
                        fudge_limits(rs,limit1)%rsnb_ljs(i)=fudge_limits(rs,limit1)%rsnb_ljs(i+j)
                        fudge_limits(rs,limit1)%rsnb_ljs(i+j)=temp_var2

                        temp_var2=fudge_limits(rs,limit1)%rsnb_lje(i)
                        fudge_limits(rs,limit1)%rsnb_lje(i)=fudge_limits(rs,limit1)%rsnb_lje(i+j)
                        fudge_limits(rs,limit1)%rsnb_lje(i+j)=temp_var2
                        exit
                    end if 
                end do  
                
                if ((test.eqv..false.)) then ! If you didn't find a match
                    if (temp_nrsnb(rs,limit1).le.temp_nrsnb(rs,limit2)) then ! if limits two is the bigger array
                        do j=1,nrsnb(rs)-i
                            iaa_limits(rs,limit1)%atnb(nrsnb(rs)+1-j,:)=iaa_limits(rs,limit1)%atnb(nrsnb(rs)-j,:)
                            fudge_limits(rs,limit1)%rsnb_ljs(nrsnb(rs)+1-j)=fudge_limits(rs,limit1)%rsnb_ljs(nrsnb(rs)-j)
                            fudge_limits(rs,limit1)%rsnb_lje(nrsnb(rs)+1-j)=fudge_limits(rs,limit1)%rsnb_lje(nrsnb(rs)-j)
                        end do 
                        iaa_limits(rs,limit1)%atnb(i,:)=iaa_limits(rs,limit2)%atnb(i,:)
                        fudge_limits(rs,limit1)%rsnb_ljs(i)=fudge_limits(rs,limit2)%rsnb_ljs(i)
                        fudge_limits(rs,limit1)%rsnb_lje(i)=fudge_limits(rs,limit2)%rsnb_lje(i)
                    else 
                        do j=1,nrsnb(rs)-i
                            iaa_limits(rs,limit2)%atnb(nrsnb(rs)+1-j,:)=iaa_limits(rs,limit2)%atnb(nrsnb(rs)-j,:)
                            fudge_limits(rs,limit2)%rsnb_ljs(nrsnb(rs)+1-j)=fudge_limits(rs,limit2)%rsnb_ljs(nrsnb(rs)-j)
                            fudge_limits(rs,limit2)%rsnb_lje(nrsnb(rs)+1-j)=fudge_limits(rs,limit2)%rsnb_lje(nrsnb(rs)-j)
                        end do 
                        
                        iaa_limits(rs,limit2)%atnb(i,:)=iaa_limits(rs,limit1)%atnb(i,:)
                        fudge_limits(rs,limit2)%rsnb_ljs(i)=fudge_limits(rs,limit1)%rsnb_ljs(i)
                        fudge_limits(rs,limit2)%rsnb_lje(i)=fudge_limits(rs,limit1)%rsnb_lje(i)
                    end if ! Martin : I should not need to make any limits zero, as the softcore should take care of all of it 
                end if 
            end if 
            
            test=.true.
            mem=test_combi(iaa_limits(rs,limit1)%atnb(i,:),iaa_limits(rs,limit2)%atnb(i,:),2)
            
            do k=1,2
                if (mem(k).eq.0) then 
                    test=.false.
                end if 
            end do
        end do 
    end do
end subroutine 

subroutine qcann_setup_srinter(limit)
    use params
    use inter
    use sequen
    use atoms
    use fyoc
    use molecule
    use polypep
    
    integer rs,alcsz,i,imol,k,j
    integer limit,temp_int
    
    ! Note that the third limit for the HIS is not going to be necessary, since changing the biotypes alone will be sufficient. 
    lj_eps(:,:)=lj_eps_limits(:,:,limit)
    lj_eps_14(:,:)=lj_eps_14_limits(:,:,limit)
    lj_sig(:,:)=lj_sig_limits(:,:,limit)
    lj_sig_14(:,:)=lj_sig_14_limits(:,:,limit)
    
    if (limit.eq.2) then 
        temp_save(:)=attyp(:)
        temp_bio_save(:)=b_type(:)
        do i=1,n 
            b_type(i)=transform_table(b_type(i))
            attyp(i)=bio_ljtyp(b_type(i))
        end do 
    else if (limit.eq.3) then 
        temp_save(:)=attyp(:)
        temp_bio_save(:)=b_type(:)
        
        ! Place holder for the biotype transform switch
        do i=1,nhis
            do j=1,at(his_state(i,1))%nsc
                
                
                b_type(at(his_state(i,1))%sc(j))=transform_table(his_eqv_table(b_type(at(his_state(i,1))%sc(j))))
                ! Martin : added after, though it was missing 
                attyp(at(his_state(i,1))%sc(j))=bio_ljtyp(b_type(at(his_state(i,1))%sc(j)))
            end do 
        end do 
    end if 

    call setup_srinter()
    
    do imol=1,nmol
        do rs=rsmol(imol,1),rsmol(imol,2)
            if (natres(rs).gt.1) then
                allocate(fudge_limits(rs,limit)%rsin((natres(rs)*(natres(rs)-1))/2))
                allocate(fudge_limits(rs,limit)%rsin_ljs((natres(rs)*(natres(rs)-1))/2))
                allocate(fudge_limits(rs,limit)%rsin_lje((natres(rs)*(natres(rs)-1))/2))
                allocate(iaa_limits(rs,limit)%atin((natres(rs)*(natres(rs)-1))/2,2))
            end if
            if (rs.lt.nseq) then
                alcsz = natres(rs)*natres(rs+1)
                if (disulf(rs).gt.0) then
                    alcsz = max(alcsz,natres(rs)*natres(disulf(rs)))
                end if

                allocate(fudge_limits(rs,limit)%rsnb(alcsz))
                allocate(fudge_limits(rs,limit)%rsnb_ljs(alcsz))
                allocate(fudge_limits(rs,limit)%rsnb_lje(alcsz))

            end if
        ! Martin : I should to need any of the things that are not ACTUALLY PARAMS, the rest will cause issues anyway

            allocate(iaa_limits(rs,limit)%bl(nrsbl(rs),2))
            allocate(iaa_limits(rs,limit)%par_bl(nrsbl(rs),MAXBOPAR+1))
            allocate(iaa_limits(rs,limit)%typ_bl(nrsbl(rs)))
            
            allocate(iaa_limits(rs,limit)%ba(nrsba(rs),3))
            allocate(iaa_limits(rs,limit)%par_ba(nrsba(rs),MAXBAPAR+1))
            allocate(iaa_limits(rs,limit)%typ_ba(nrsba(rs)))

            allocate(iaa_limits(rs,limit)%di(nrsdi(rs),4))
            allocate(iaa_limits(rs,limit)%par_di(nrsdi(rs),MAXDIPAR))
            allocate(iaa_limits(rs,limit)%typ_di(nrsdi(rs)))
            
            allocate(iaa_limits(rs,limit)%impt(3*nrsimpt(rs),4))
            allocate(iaa_limits(rs,limit)%par_impt(3*nrsimpt(rs),MAXDIPAR))
            allocate(iaa_limits(rs,limit)%typ_impt(3*nrsimpt(rs)))
            ! Martin ": those should be allocated or just above where we actually get the correct alcsz
            ! Martin : I will try to allocate everything, even the things that are not useuful, because
            ! it may resolves my seg fault issue
        end do 
    end do

    fudge_limits(:,limit)=fudge(:)
    iaa_limits(:,limit)=iaa(:)
    
    temp_nrsin(:,limit)=nrsintra(:)
    temp_nrsnb(:,limit)=nrsnb(:)
    
    if (limit.ne.1) then 
        attyp(:)=temp_save(:)
        b_type(:)=temp_bio_save(:)
        call deallocate_iaa()
    end if 
    
end subroutine 

subroutine qcann_iaa_index_corrections(limit1,limit2)
    use inter
    use sequen
    use energies
    use aminos
    use atoms
    use params
    
    logical good,reverse
    RTYPE, allocatable:: temp_var1(:),  temp_var2(:)
    
    integer, allocatable:: temp_var5(:,:),temp_var6(:,:)
    integer, allocatable:: memory(:) 
    integer limit1,limit2
    integer i,j,rs,k,temp_int
    
    allocate(temp_var1(3))
    allocate(temp_var2(5))
    
    do rs=1,nseq
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Bond angles
        do i=1,nrsba(rs) ! This part is to get the right order
            allocate(memory(size(iaa(rs)%ba(i,:))))
            good=.true.
            memory(:)=test_combi(iaa_limits(rs,limit1)%ba(i,:),iaa_limits(rs,limit2)%ba(i,:),size(iaa(rs)%ba(i,:)))            

            do k=1,size(memory)
                if (memory(k).eq.0) then 
                    good=.false.
                end if 
            end do 
            if (good.eqv..false.) then 
             ! So if we find a residue array that is not strictly indentical, we go in the reste of the array to make a permutation 
                if (temp_nrsba(rs,limit1).lt.temp_nrsba(rs,limit2)) then ! Thta is if the first array is shorter
                    do j=1,nrsba(rs)-i
                        memory(:)=test_combi(iaa_limits(rs,limit1)%ba(i+j,:),iaa_limits(rs,limit2)%ba(i,:),3)

                        if ((memory(2).eq.2).and.&
                        &(((memory(3).eq.3).and.(memory(1).eq.1)).or.((memory(3).eq.1).and.(memory(1).eq.3)))) then                       

                            temp_var1(:)=iaa_limits(rs,limit1)%ba(i,:)
                            iaa_limits(rs,limit1)%ba(i,:)=iaa_limits(rs,limit1)%ba(i+j,:)
                            iaa_limits(rs,limit1)%ba(i+j,:)=temp_var1(:)

                            temp_var2(:)=iaa_limits(rs,limit1)%par_ba(i,:)
                            iaa_limits(rs,limit1)%par_ba(i,:)=iaa_limits(rs,limit1)%par_ba(i+j,:)
                            iaa_limits(rs,limit1)%par_ba(i+j,:)=temp_var2(:)

                            temp_int =iaa_limits(rs,limit1)%typ_ba(i)
                            iaa_limits(rs,limit1)%typ_ba(i)=iaa_limits(rs,limit1)%typ_ba(i+j)
                            iaa_limits(rs,limit1)%typ_ba(i+j)=temp_int 
                            exit 
                        end if   
                    end do 
                else    
                    do j=1,nrsba(rs)-i
                        memory(:)=test_combi(iaa_limits(rs,limit2)%ba(i+j,:),iaa_limits(rs,limit1)%ba(i,:),3)
                        if ((memory(2).eq.2).and.&
                        &(((memory(3).eq.3).and.(memory(1).eq.1)).or.((memory(3).eq.1).and.(memory(1).eq.3)))) then                       

                            temp_var1(:)=iaa_limits(rs,limit2)%ba(i,:)
                            iaa_limits(rs,limit2)%ba(i,:)=iaa_limits(rs,limit2)%ba(i+j,:)
                            iaa_limits(rs,limit2)%ba(i+j,:)=temp_var1(:)

                            temp_var2(:)=iaa_limits(rs,limit2)%par_ba(i,:)
                            iaa_limits(rs,limit2)%par_ba(i,:)=iaa_limits(rs,limit2)%par_ba(i+j,:)
                            iaa_limits(rs,limit2)%par_ba(i+j,:)=temp_var2(:)

                            temp_int =iaa_limits(rs,limit2)%typ_ba(i)
                            iaa_limits(rs,limit2)%typ_ba(i)=iaa_limits(rs,limit2)%typ_ba(i+j)
                            iaa_limits(rs,limit2)%typ_ba(i+j)=temp_int
                            exit 
                        end if   
                    end do 
                end if 
            end if  
            deallocate(memory)
        end do 
        ! This part is to correct for changes of types
        ! I must start with 2 cause that is where all the eff angles are
        do i=1,nrsba(rs)
            
            if ((iaa_limits(rs,limit1)%typ_ba(i).eq.2).and.(iaa_limits(rs,limit2)%typ_ba(i).eq.1)) then 
                iaa_limits(rs,limit2)%typ_ba(i)=2
                iaa_limits(rs,limit2)%par_ba(i,3)=0.
                iaa_limits(rs,limit2)%par_ba(i,4)=iaa_limits(rs,limit1)%par_ba(i,4)
                
            else if ((iaa_limits(rs,limit1)%typ_ba(i).eq.1).and.(iaa_limits(rs,limit2)%typ_ba(i).eq.2)) then 
                iaa_limits(rs,limit1)%typ_ba(i)=2
                iaa_limits(rs,limit1)%par_ba(i,3)=0.
                iaa_limits(rs,limit1)%par_ba(i,4)=iaa_limits(rs,limit2)%par_ba(i,4)
                
            else if ((iaa_limits(rs,limit1)%typ_ba(i).eq.3).or.(iaa_limits(rs,limit2)%typ_ba(i).eq.3)) then 
                print *,"Wrong angle type for atom ",iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),iaa(rs)%ba(i,3)," : RETI does&
                & not support cosine based harmonic potential"
                call fexit()
                
            else if ((iaa_limits(rs,limit1)%typ_ba(i).eq.1).and.(iaa_limits(rs,limit2)%typ_ba(i).eq.0)) then 

                iaa_limits(rs,limit2)%typ_ba(i)=1
                iaa_limits(rs,limit2)%par_ba(i,:)=iaa_limits(rs,limit1)%par_ba(i,:)

            else if ((iaa_limits(rs,limit1)%typ_ba(i).eq.2).and.(iaa_limits(rs,limit2)%typ_ba(i).eq.0)) then 
                iaa_limits(rs,limit2)%typ_ba(i)=2
                iaa_limits(rs,limit2)%par_ba(i,:)=iaa_limits(rs,limit1)%par_ba(i,:)
            end if
            
        end do
        
        do i=1,nrsba(rs) ! This par is to zero out whatever should not be here because of dummy
            if ((abs(atr_limits(iaa_limits(rs,limit1)%ba(i,1),1)).lt.1.0D-5).or.&
            &(abs(atr_limits(iaa_limits(rs,limit1)%ba(i,2),1)).lt.1.0D-5)&
            &.or.(abs(atr_limits(iaa_limits(rs,limit1)%ba(i,3),1)).lt.1.0D-5)) then 
            
                iaa_limits(rs,limit1)%par_ba(i,1)=0.
                iaa_limits(rs,limit1)%par_ba(i,3)=0.
                
            else if ((abs(atr_limits(iaa_limits(rs,limit1)%ba(i,1),2)).lt.1.0D-5).or.&
                &(abs(atr_limits(iaa_limits(rs,limit1)%ba(i,2),2)).lt.1.0D-5).or.&
                &(abs(atr_limits(iaa_limits(rs,limit1)%ba(i,3),2)).lt.1.0D-5)) then
                
                iaa_limits(rs,limit2)%par_ba(i,1)=0.
                iaa_limits(rs,limit2)%par_ba(i,3)=0.
                
            end if
        end do
    end do 
    
    deallocate(temp_var1)
    deallocate(temp_var2)

    allocate(temp_var1(2))
    allocate(temp_var2(4))
    
    ! MARTIN :: super important, I need to align the smallest array to the biggest one ! Because 
    do rs=1,nseq
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Bond length
        do i=1,nrsbl(rs) ! This part is to get the right order
            allocate(memory(size(iaa(rs)%bl(i,:))))
            good=.true.
            memory(:)=test_combi(iaa_limits(rs,limit1)%bl(i,:),iaa_limits(rs,limit2)%bl(i,:),size(iaa(rs)%bl(i,:)))            
            
            do k=1,size(memory)
                if (memory(k).eq.0) then 
                    good=.false.
                end if 
            end do 
            if (good.eqv..false.) then 
                if (temp_nrsbl(rs,limit1).lt.temp_nrsbl(rs,limit2)) then ! Thta is if the first array is shorter
                 ! So if we find a residue array that is not strictly indentical,
                ! we go in the reste of the array to make a permutation 
                    do j=1,nrsbl(rs)-i
                        memory(:)=test_combi(iaa_limits(rs,limit1)%bl(i+j,:),iaa_limits(rs,limit2)%bl(i,:),size(iaa(rs)%bl(i,:)))

                        if (((memory(2).eq.2).and.(memory(1).eq.1)).or.((memory(1).eq.2).and.(memory(2).eq.1))) then                       

                            temp_var1(:)=iaa_limits(rs,limit1)%bl(i,:)
                            iaa_limits(rs,limit1)%bl(i,:)=iaa_limits(rs,limit1)%bl(i+j,:)
                            iaa_limits(rs,limit1)%bl(i+j,:)=temp_var1(:)

                            temp_var2(:)=iaa_limits(rs,limit1)%par_bl(i,:)
                            iaa_limits(rs,limit1)%par_bl(i,:)=iaa_limits(rs,limit1)%par_bl(i+j,:)
                            iaa_limits(rs,limit1)%par_bl(i+j,:)=temp_var2(:)

                            temp_int =iaa_limits(rs,limit1)%typ_bl(i)
                            iaa_limits(rs,limit1)%typ_bl(i)=iaa_limits(rs,limit1)%typ_bl(i+j)
                            iaa_limits(rs,limit1)%typ_bl(i+j)=temp_int 
                            exit 

                        end if   
                    end do 
                else 
                    do j=1,nrsbl(rs)-i
                        memory(:)=test_combi(iaa_limits(rs,limit2)%bl(i+j,:),iaa_limits(rs,limit1)%bl(i,:),size(iaa(rs)%bl(i,:)))

                        if (((memory(2).eq.2).and.(memory(1).eq.1)).or.((memory(1).eq.2).and.(memory(2).eq.1))) then                       

                            temp_var1(:)=iaa_limits(rs,limit2)%bl(i,:)
                            iaa_limits(rs,limit2)%bl(i,:)=iaa_limits(rs,2)%bl(i+j,:)
                            iaa_limits(rs,limit2)%bl(i+j,:)=temp_var1(:)

                            temp_var2(:)=iaa_limits(rs,limit2)%par_bl(i,:)
                            iaa_limits(rs,limit2)%par_bl(i,:)=iaa_limits(rs,limit2)%par_bl(i+j,:)
                            iaa_limits(rs,limit2)%par_bl(i+j,:)=temp_var2(:)

                            temp_int =iaa_limits(rs,limit2)%typ_bl(i)
                            iaa_limits(rs,limit2)%typ_bl(i)=iaa_limits(rs,limit2)%typ_bl(i+j)
                            iaa_limits(rs,limit2)%typ_bl(i+j)=temp_int 
                            exit 

                        end if   
                    end do 
                end if 
            end if  
            deallocate(memory)
        end do 
        ! This part is to correct for changes of types
        ! I must start with 2 cause that is where all the eff angles are
        do i=1,nrsbl(rs)
            if (iaa_limits(rs,limit1)%typ_bl(i).ne.iaa_limits(rs,limit2)%typ_bl(i)) then 
                ! The following part is because the interaction is not defined for the dyummy atom. 
                ! Thus I set it to that of the other side 
                if (iaa_limits(rs,limit1)%typ_bl(i).eq.0) then 
                    iaa_limits(rs,limit1)%typ_bl(i)=iaa_limits(rs,limit2)%typ_bl(i)
                else if (iaa_limits(rs,limit2)%typ_bl(i).eq.0) then 
                    iaa_limits(rs,limit2)%typ_bl(i)=iaa_limits(rs,limit1)%typ_bl(i)
                else 
                    
                    print *,"Unfortunately, RETI does not support varying bond length types(",iaa_limits(rs,limit1)%typ_bl(i),&
                    &iaa_limits(rs,limit2)%typ_bl(i),"). Exiting."
                    print *,"The previous message corresponds to atoms :",iaa(rs)%bl(i,:),&
                    &iaa_limits(rs,limit1)%bl(i,:),iaa_limits(rs,limit2)%bl(i,:)
                    flush(6)
                    call fexit()
                end if 
            end if
        end do
        
        do i=1,nrsbl(rs) ! This par is to zero out whatever should not be here because of dummy
            if ((abs(atr_limits(iaa(rs)%bl(i,1),limit1)).lt.1.0D-5).or.(abs(atr_limits(iaa(rs)%bl(i,2),limit1)).lt.1.0D-5)) then 
                iaa_limits(rs,1)%par_bl(i,1)=0.
            else if ((abs(atr_limits(iaa(rs)%bl(i,1),limit2)).lt.1.0D-5).or.(abs(atr_limits(iaa(rs)%bl(i,2),limit2)).lt.1.0D-5)) then 
                iaa_limits(rs,2)%par_bl(i,1)=0.
            end if 
        end do
    end do 
    deallocate(temp_var1)
    deallocate(temp_var2)

    allocate(temp_var1(4))! save atoms
    allocate(temp_var2(9))! save params
    do rs=1,nseq
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Bond dihedral
        do i=1,nrsdi(rs) ! This part is to get the right order
            allocate(memory(size(iaa(rs)%di(i,:))))
            good=.true.
            memory(:)=test_combi(iaa_limits(rs,limit1)%di(i,:),iaa_limits(rs,limit2)%di(i,:),size(iaa(rs)%di(i,:)))            
            do k=1,size(memory)
                if (memory(k).eq.0) then 
                    good=.false.
                end if 
            end do 
            if (good.eqv..false.) then 
             ! So if we find a residue array that is not strictly indentical, we go in the reste of the array to make a permutation 
                if (temp_nrsdi(rs,limit1).lt.temp_nrsdi(rs,limit2)) then ! Thta is if the first array is shorter
                    do j=1,nrsdi(rs)-i
                        memory(:)=test_combi(iaa_limits(rs,limit1)%di(i+j,:),iaa_limits(rs,limit2)%di(i,:),size(iaa(rs)%di(i,:)))
                        good=.true. ! Reusing variables
                        reverse=.true.
                        do k=1,size(iaa(rs)%di(i,:))
                            if ((memory(k).ne.k).and.(memory(k).ne.size(iaa(rs)%di(i,:))-k+1)) then 
                                good=.false.
                            end if 
                            if (memory(k).ne.size(iaa(rs)%di(i,:))-k+1) then 
                                reverse=.false.
                            end if 
                        end do 

                        if (good.eqv..true.) then                       
                            temp_var1(:)=iaa_limits(rs,limit1)%di(i,:)
                            iaa_limits(rs,limit1)%di(i,:)=iaa_limits(rs,limit1)%di(i+j,:)
                            iaa_limits(rs,limit1)%di(i+j,:)=temp_var1(:)

                            temp_var2(:)=iaa_limits(rs,limit1)%par_di(i,:)
                            iaa_limits(rs,limit1)%par_di(i,:)=iaa_limits(rs,limit1)%par_di(i+j,:)
                            iaa_limits(rs,limit1)%par_di(i+j,:)=temp_var2(:)

                            temp_int =iaa_limits(rs,limit1)%typ_di(i)
                            iaa_limits(rs,limit1)%typ_di(i)=iaa_limits(rs,limit1)%typ_di(i+j)
                            iaa_limits(rs,limit1)%typ_di(i+j)=temp_int 
                            exit 
                        end if   
                    end do 
                else    
                    
                    
                    do j=1,nrsdi(rs)-i
                        memory(:)=test_combi(iaa_limits(rs,limit2)%di(i+j,:),iaa_limits(rs,limit1)%di(i,:),size(iaa(rs)%di(i,:)))
                        good=.true. ! Reusing variables
                        do k=1,size(iaa(rs)%di(i,:))
                            if ((memory(k).ne.k).and.(memory(k).ne.size(iaa(rs)%di(i,:))-k+1)) then 
                                good=.false.
                            end if 
                        end do 
                        if (good.eqv..true.) then                        

                            temp_var1(:)=iaa_limits(rs,limit2)%di(i,:)
                            iaa_limits(rs,limit2)%di(i,:)=iaa_limits(rs,limit2)%di(i+j,:)
                            iaa_limits(rs,limit2)%di(i+j,:)=temp_var1(:)

                            temp_var2(:)=iaa_limits(rs,limit2)%par_di(i,:)
                            iaa_limits(rs,limit2)%par_di(i,:)=iaa_limits(rs,limit2)%par_di(i+j,:)
                            iaa_limits(rs,limit2)%par_di(i+j,:)=temp_var2(:)

                            temp_int =iaa_limits(rs,limit2)%typ_di(i)
                            iaa_limits(rs,limit2)%typ_di(i)=iaa_limits(rs,limit2)%typ_di(i+j)
                            iaa_limits(rs,limit2)%typ_di(i+j)=temp_int 
                            exit 
                        end if   
                    end do 
                end if 
            end if  
            deallocate(memory)
        end do 
        allocate(memory(4))
        
        do i=1,nrsdi(rs)
            ! Here, we have already determined that the same index is the same dihedral, 
            ! but it may still have a difference in the order of the atom, so reindexing just in case
            ! Thus we need to realign one or the other 
            ! Since at this point the atoms in iaa are the same as in iaalim2, we will switch lim1
            memory(:)=test_combi(iaa_limits(rs,limit1)%di(i,:),iaa_limits(rs,limit2)%di(i,:),size(iaa(rs)%di(i,:)))
            
            if ((memory(1).eq.0).or.(memory(2).eq.0).or.(memory(3).eq.0).or.(memory(4).eq.0)) then 
                ! If this prints, then you have exclusive contacts in both lists. Should not happen.
                print *,"Problem with the dihedral interpolation arrays."
                print *,rs,i
                print *,memory(:)
                print *,iaa_limits(rs,limit1)%di(i,:)
                print *,iaa_limits(rs,limit2)%di(i,:)
                flush(6)
                call fexit()
            end if 
            
            if (iaa_limits(rs,limit1)%par_di(i,1).lt.1.0D-5) then 
                iaa_limits(rs,limit1)%di(i,:)=iaa_limits(rs,limit2)%di(i,:)
                iaa_limits(rs,limit1)%par_di(i,:)=iaa_limits(rs,limit2)%par_di(i,:)

            else if (iaa_limits(rs,limit2)%par_di(i,1).lt.1.0D-5) then
                iaa_limits(rs,limit2)%di(i,:)=iaa_limits(rs,limit1)%di(i,memory(:))
                iaa_limits(rs,limit2)%par_di(i,:)=iaa_limits(rs,limit1)%par_di(i,:)
                
            end if 
        end do
        deallocate(memory)
        ! This part is to correct for changes of types
        ! I must start with 2 cause that is where all the eff angles are
        do i=1,nrsdi(rs)
            if (iaa_limits(rs,limit1)%typ_di(i).ne.iaa_limits(rs,limit2)%typ_di(i)) then
                if (iaa_limits(rs,limit1)%typ_di(i).eq.0) then
                    iaa_limits(rs,limit1)%typ_di(i)=iaa_limits(rs,limit2)%typ_di(i)
                else if (iaa_limits(rs,limit2)%typ_di(i).eq.0) then
                    iaa_limits(rs,limit2)%typ_di(i)=iaa_limits(rs,limit1)%typ_di(i)
                else
                    print *,iaa_limits(rs,limit1)%di(i,:)
                    print *,iaa_limits(rs,limit2)%di(i,:)
                    print *,"Unfortunately, RETI does not support varying dihedral types. Exiting."
                    call fexit()
                end if 
            end if
        end do
        
        do i=1,nrsdi(rs) ! This par is to zero out whatever should not be here because of dummy
            if ((abs(atr_limits(iaa(rs)%di(i,1),limit1)).lt.1.0D-5).or.(abs(atr_limits(iaa(rs)%di(i,limit2),1)).lt.1.0D-5).or.&
            &(abs(atr_limits(iaa(rs)%di(i,3),limit1)).lt.1.0D-5).or.(abs(atr_limits(iaa(rs)%di(i,4),limit1)).lt.1.0D-5)) then 
                iaa_limits(rs,limit1)%par_di(i,1)=0.
                iaa_limits(rs,limit1)%par_di(i,2)=0.
                iaa_limits(rs,limit1)%par_di(i,3)=0.
                iaa_limits(rs,limit1)%par_di(i,4)=0.
            else if ((abs(atr_limits(iaa(rs)%di(i,1),limit2)).lt.1.0D-5).or.(abs(atr_limits(iaa(rs)%di(i,2),limit2)).lt.1.0D-5).or.&
            &(abs(atr_limits(iaa(rs)%di(i,3),limit2)).lt.1.0D-5).or.(abs(atr_limits(iaa(rs)%di(i,4),limit1)).lt.1.0D-5)) then 
                iaa_limits(rs,limit2)%par_di(i,1)=0.
                iaa_limits(rs,limit2)%par_di(i,2)=0.
                iaa_limits(rs,limit2)%par_di(i,3)=0.
                iaa_limits(rs,limit2)%par_di(i,4)=0.
            end if 
        end do
    end do 
    deallocate(temp_var1)
    deallocate(temp_var2)
    
    allocate(temp_var1(4))
    allocate(temp_var2(MAXDIPAR))

    
    
    do rs=1,nseq
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Bond dihedral
        do i=1,nrsimpt(rs) ! This part is to get the right order

            allocate(memory(size(iaa(rs)%impt(i,:))))
            good=.true.
            memory(:)=test_combi(iaa_limits(rs,limit1)%impt(i,:),iaa_limits(rs,limit2)%impt(i,:),size(iaa(rs)%impt(i,:)))
            do k=1,size(memory)
                if (memory(k).eq.0) then 
                    good=.false.
                end if 
            end do 
            
            
            if (good.eqv..false.) then 
             ! So if we find a residue array that is not strictly indentical, we go in the rest of the array to make a permutation 
                if (temp_nrsimpt(rs,limit1).lt.temp_nrsimpt(rs,limit2)) then ! That is if the first array is shorter
                    temp_int = 999 ! Just for additional debugging purposes
                    do j=1,nrsimpt(rs) -i
                        memory(:)=test_combi(iaa_limits(rs,limit1)%impt(i+j,:),iaa_limits(rs,limit2)%impt(i,:),&
                        &size(iaa(rs)%impt(i,:)))
                        good=.true. ! Reusing variables 
                        do k=1,size(iaa(rs)%impt(i,:))
                            if ((memory(k).ne.k).and.(memory(k).ne.size(iaa(rs)%impt(i,:))-k+1)) then 
                                good=.false.
                            end if 
                        end do 
                        if (good.eqv..true.) then 
                            temp_var1(:)=iaa_limits(rs,limit1)%impt(i,:)
                            iaa_limits(rs,limit1)%impt(i,:)=iaa_limits(rs,limit1)%impt(i+j,:)
                            iaa_limits(rs,limit1)%impt(i+j,:)=temp_var1(:)

                            temp_var2(:)=iaa_limits(rs,limit1)%par_impt(i,:)
                            iaa_limits(rs,limit1)%par_impt(i,:)=iaa_limits(rs,limit1)%par_impt(i+j,:)
                            iaa_limits(rs,limit1)%par_impt(i+j,:)=temp_var2(:)

                            temp_int =iaa_limits(rs,limit1)%typ_impt(i)
                            iaa_limits(rs,limit1)%typ_impt(i)=iaa_limits(rs,limit1)%typ_impt(i+j)
                            iaa_limits(rs,limit1)%typ_impt(i+j)=temp_int
                            exit 
                        end if   
                    end do
                    ! IN the case where the "contact" was NOT found, just add it 
                    if (good.eqv..false.) then
                        allocate(temp_var5(size(iaa_limits(rs,limit1)%impt(:,1)),4))
                        allocate(temp_var6(size(iaa_limits(rs,limit1)%par_impt(:,1)),MAXDIPAR))
                        

                        temp_var5(:,:)=iaa_limits(rs,limit1)%impt(:,:)
                        temp_var6(:,:)=iaa_limits(rs,limit1)%par_impt(:,:)
                        
                        iaa_limits(rs,limit1)%impt(i,:)=iaa_limits(rs,limit2)%impt(i,:)
                        iaa_limits(rs,limit1)%par_impt(i,:)=iaa_limits(rs,limit2)%par_impt(i,:)
                        
                        ! Shift all elements by one in order to introduce new "contact"
                        do j=i+1,nrsimpt(rs)
                            iaa_limits(rs,limit1)%impt(j,:)=temp_var5(j-1,:)
                            iaa_limits(rs,limit1)%par_impt(j,:)=temp_var6(j-1,:)
                        end do 

                        deallocate(temp_var5)
                        deallocate(temp_var6)
                    end if 
                else    
                    do j=1,nrsimpteff(rs) -i
                        memory(:)=test_combi(iaa_limits(rs,limit2)%impt(i+j,:),iaa_limits(rs,limit1)%impt(i,:),&
                        &size(iaa(rs)%impt(i,:)))
                        good=.true. ! Reusing variables
                        do k=1,size(iaa(rs)%impt(i,:))
                            if ((memory(k).ne.k).and.(memory(k).ne.size(iaa(rs)%impt(i,:))-k+1)) then 
                                good=.false.
                            end if 
                        end do 
                        if (good.eqv..true.) then
                            temp_var1(:)=iaa_limits(rs,limit2)%impt(i,:)
                            iaa_limits(rs,limit2)%impt(i,:)=iaa_limits(rs,limit2)%impt(i+j,:)
                            iaa_limits(rs,limit2)%impt(i+j,:)=temp_var1(:)

                            temp_var2(:)=iaa_limits(rs,limit2)%par_impt(i,:)
                            iaa_limits(rs,limit2)%par_impt(i,:)=iaa_limits(rs,limit2)%par_impt(i+j,:)
                            iaa_limits(rs,limit2)%par_impt(i+j,:)=temp_var2(:)

                            temp_int =iaa_limits(rs,limit2)%typ_impt(i)
                            iaa_limits(rs,limit2)%typ_impt(i)=iaa_limits(rs,limit2)%typ_impt(i+j)
                            iaa_limits(rs,limit2)%typ_impt(i+j)=temp_int 
                            exit 
                        end if   
                    end do
                    ! This means that the contact was not found in the "other list"
                    if (good.eqv..false.) then 
                        allocate(temp_var5(size(iaa_limits(rs,limit2)%impt(:,1)),4))
                        allocate(temp_var6(size(iaa_limits(rs,limit2)%par_impt(:,1)),MAXDIPAR))
                        
                        temp_var5(:,:)=iaa_limits(rs,limit2)%impt(:,:)
                        temp_var6(:,:)=iaa_limits(rs,limit2)%par_impt(:,:)
                        
                        iaa_limits(rs,limit2)%impt(i,:)=iaa_limits(rs,limit1)%impt(i,:)
                        iaa_limits(rs,limit2)%par_impt(i,:)=iaa_limits(rs,limit1)%par_impt(i,:)
                        
                        ! Shift all elements by one in order to introduce new "contact"
                        do j=i+1,nrsimpt(rs)
                            iaa_limits(rs,limit2)%impt(j,:)=temp_var5(j-1,:)
                            iaa_limits(rs,limit2)%par_impt(j,:)=temp_var6(j-1,:)
                        end do 
                        
                        deallocate(temp_var5)
                        deallocate(temp_var6)
                    end if 
                end if 
            end if  
            deallocate(memory)
        end do 
        ! This part is to correct for changes of types
        ! I must start with 2 cause that is where all the eff angles are
        do i=1,nrsimpt(rs)
            if (iaa_limits(rs,limit1)%typ_impt(i).eq.0) then 
                iaa_limits(rs,limit1)%typ_impt(i)=iaa_limits(rs,limit2)%typ_impt(i)
            else if (iaa_limits(rs,limit2)%typ_impt(i).eq.0) then 
                iaa_limits(rs,limit2)%typ_impt(i)=iaa_limits(rs,limit1)%typ_impt(i)
            else if (iaa_limits(rs,limit1)%typ_impt(i).ne.iaa_limits(rs,limit2)%typ_impt(i)) then 
                print *,"Unfortunately, RETI does not support varying improper dihedral types. Exiting."
                call fexit() ! Martin : just seeing further to catch the problem 
                
            end if
        end do
        
        do i=1,nrsimpt(rs) ! This par is to zero out whatever should not be here because of dummy
            if ((abs(atr_limits(iaa(rs)%impt(i,1),limit1)).lt.1.0D-5).or.(abs(atr_limits(iaa(rs)%impt(i,limit2),1)).lt.1.0D-5).or.&
            &(abs(atr_limits(iaa(rs)%impt(i,3),limit1)).lt.1.0D-5).or.(abs(atr_limits(iaa(rs)%impt(i,4),limit1)).lt.1.0D-5)) then 
                iaa_limits(rs,limit1)%par_impt(i,1)=0.
                iaa_limits(rs,limit1)%par_impt(i,2)=0.
                iaa_limits(rs,limit1)%par_impt(i,3)=0.
                iaa_limits(rs,limit1)%par_impt(i,4)=0.
                
            else if ((abs(atr_limits(iaa(rs)%impt(i,1),limit2)).lt.1.0D-5).or.(abs(atr_limits(iaa(rs)%impt(i,2),limit2)).lt.1.0D-5)&
                &.or.(abs(atr_limits(iaa(rs)%impt(i,3),limit2)).lt.1.0D-5).or.(abs(atr_limits(iaa(rs)%impt(i,4),limit1))&
                &.lt.1.0D-5)) then 
                iaa_limits(rs,limit2)%par_impt(i,1)=0.
                iaa_limits(rs,limit2)%par_impt(i,2)=0.
                iaa_limits(rs,limit2)%par_impt(i,3)=0.
                iaa_limits(rs,limit2)%par_impt(i,4)=0.
            end if 
        end do
        
        allocate(memory(4))
        do i=1,nrsimpt(rs)
            ! Here, we have already determined that the same index is the same dihedral, 
            ! but it may still have a difference in the order of the atom, so reindexing just in case
            ! Thus we need to realign one or the other 
            ! Since at this point the atoms in iaa are the same as in iaalim2, we will switch lim1
            memory(:)=test_combi(iaa_limits(rs,limit1)%impt(i,:),iaa_limits(rs,limit2)%impt(i,:),size(iaa(rs)%impt(i,:)))
            
            if ((memory(1).eq.0).or.(memory(2).eq.0).or.(memory(3).eq.0).or.(memory(4).eq.0)) then 
                ! If this prints, then you have exclusive contacts in both lists. Should not happen
                print *,"impt","If this prints, then you have exclusive contacts in both lists. Should not happen"
                print *,rs,i
                print *,"Problem is with limits ",limit1,limit2
                print *,memory(:)
                print *,iaa_limits(rs,limit1)%impt(i,:)
                print *,iaa_limits(rs,limit2)%impt(i,:)
                flush(6)
                call fexit()
            end if 
            if (iaa_limits(rs,limit1)%par_impt(i,1).lt.1.0D-5) then 
                iaa_limits(rs,limit1)%impt(i,:)=iaa_limits(rs,limit2)%impt(i,:)
                iaa_limits(rs,limit1)%par_impt(i,:)=iaa_limits(rs,limit2)%par_impt(i,:)
            else if (iaa_limits(rs,limit2)%par_impt(i,1).lt.1.0D-5) then

                iaa_limits(rs,limit2)%impt(i,:)=iaa_limits(rs,limit1)%impt(i,memory(:))
                iaa_limits(rs,limit2)%par_impt(i,:)=iaa_limits(rs,limit1)%par_impt(i,:)
            end if 
        end do
        deallocate(memory)

    end do 
    deallocate(temp_var1)
    deallocate(temp_var2)

    allocate(memory(4))
    
     ! Here, we have already determined that the same index is the same dihedral, 
     ! but it may still have a difference in the order of the atom, which can influence the energy. 
     ! Thus we need to realign one or the other 
     ! Since at this point the atoms in iaa are the same as in iaalim2, we will switch lim1
     ! allocate(temp_var1(4))
    do rs=1,nseq
        do i=1,nrsdi(rs)
            memory(:)=test_combi(iaa_limits(rs,limit1)%di(i,:),iaa_limits(rs,limit2)%di(i,:),size(iaa(rs)%di(i,:)))     
            do j=1,size(memory(:))
                if (memory(j).ne.j) then  
                    print *,"Got a problem with DI post",memory(:),rs,i
                    print *,"Problem is with limits ",limit1,limit2
                    print *,"lim1",iaa_limits(rs,limit1)%di(i,:)
                    print *,"lim2",iaa_limits(rs,limit2)%di(i,:)
                    exit
                end if 
            end do
        end do 
    end do 
    
    deallocate(memory)
    allocate(memory(4))
    do rs=1,nseq
        do i=1,nrsimpt(rs)
            memory(:)=test_combi(iaa_limits(rs,limit1)%impt(i,:),iaa_limits(rs,limit2)%impt(i,:),size(iaa(rs)%impt(i,:)))            
            do j=1,size(memory(:))
                if (memory(j).ne.j) then 
                    print *,"Got a problem with IMPT",memory(:),rs,i
                    
                    print *,"Problem is with limits ",limit1,limit2
                    print *,"lim",limit1,iaa_limits(rs,limit1)%impt(i,1),iaa_limits(rs,limit1)%impt(i,2),&
                    iaa_limits(rs,limit1)%impt(i,3),&
                    &iaa_limits(rs,1)%impt(i,4)
                    print *,"lim",limit2,iaa_limits(rs,limit2)%impt(i,1),iaa_limits(rs,limit2)%impt(i,2),&
                    &iaa_limits(rs,limit2)%impt(i,3),&
                    &iaa_limits(rs,2)%impt(i,4)
                    print *,"memory",memory
                    print *,"the index is ",i
                    print *,"the residue is ",rs
                    exit
                end if 
            end do
        end do 
    end do 
    deallocate(memory)

    iaa(:)=iaa_limits(:,limit2)
    fudge(:)=fudge_limits(:,limit2)
    allocate(memory(4))
    
    do rs=1,nseq
        do i=1,nrsdi(rs)
            if ((iaa_limits(rs,limit1)%di(i,1).ne.iaa(rs)%di(i,1)).or.&
            &(iaa_limits(rs,limit1)%di(i,2).ne.iaa(rs)%di(i,2)).or.&
            &(iaa_limits(rs,limit1)%di(i,3).ne.iaa(rs)%di(i,3)).or.&
            &(iaa_limits(rs,limit1)%di(i,4).ne.iaa(rs)%di(i,4)).or.&
            &(iaa_limits(rs,limit2)%di(i,1).ne.iaa(rs)%di(i,1)).or.&
            &(iaa_limits(rs,limit2)%di(i,2).ne.iaa(rs)%di(i,2)).or.&
            &(iaa_limits(rs,limit2)%di(i,3).ne.iaa(rs)%di(i,3)).or.&
            &(iaa_limits(rs,limit2)%di(i,4).ne.iaa(rs)%di(i,4))) then
                iaa(rs)%di(i,:)=iaa_limits(rs,limit2)%di(i,:)
            end if 
        end do 
    end do 

    
    do rs=1,nseq
        do i=1,nrsimpt(rs)
            if ((iaa_limits(rs,limit1)%impt(i,1).ne.iaa(rs)%impt(i,1)).or.&
            &(iaa_limits(rs,limit1)%impt(i,2).ne.iaa(rs)%impt(i,2)).or.&
            &(iaa_limits(rs,limit1)%impt(i,3).ne.iaa(rs)%impt(i,3)).or.&
            &(iaa_limits(rs,limit1)%impt(i,4).ne.iaa(rs)%impt(i,4)).or.&
            
            &(iaa_limits(rs,limit2)%impt(i,1).ne.iaa(rs)%impt(i,1)).or.&
            &(iaa_limits(rs,limit2)%impt(i,2).ne.iaa(rs)%impt(i,2)).or.&
            &(iaa_limits(rs,limit2)%impt(i,3).ne.iaa(rs)%impt(i,3)).or.&
            &(iaa_limits(rs,limit2)%impt(i,4).ne.iaa(rs)%impt(i,4))) then
                iaa(rs)%impt(i,:)=iaa_limits(rs,limit2)%impt(i,:)
            end if 
        end do 
    end do 
    deallocate(memory)
    
    do rs=1,nseq
        do i=1,nrsba(rs)
            if ((iaa_limits(rs,limit1)%ba(i,1).ne.iaa_limits(rs,limit2)%ba(i,1)).or.&
            &(iaa_limits(rs,limit1)%ba(i,2).ne.iaa_limits(rs,limit2)%ba(i,2)).or.&
            &(iaa_limits(rs,limit1)%ba(i,3).ne.iaa_limits(rs,limit2)%ba(i,3)).or.&
            
            &(iaa_limits(rs,limit1)%ba(i,1).ne.iaa(rs)%ba(i,1)).or.&
            &(iaa_limits(rs,limit1)%ba(i,2).ne.iaa(rs)%ba(i,2)).or.&
            &(iaa_limits(rs,limit1)%ba(i,3).ne.iaa(rs)%ba(i,3))) then
                print *,"ba",iaa(rs)%ba(i,:)
                print *,"ba_lim1",iaa_limits(rs,limit1)%ba(i,:)
                print *,"ba_lim2",iaa_limits(rs,limit2)%ba(i,:)

            end if 
        end do
    end do 
    
    do rs=1,nseq
        do i=1,nrsbl(rs)
            if ((iaa_limits(rs,limit1)%bl(i,1).ne.iaa_limits(rs,limit2)%bl(i,1)).or.&
            &(iaa_limits(rs,limit1)%bl(i,2).ne.iaa_limits(rs,limit2)%bl(i,2)).or.&
            &(iaa_limits(rs,limit1)%bl(i,1).ne.iaa(rs)%bl(i,1)).or.&
            &(iaa_limits(rs,limit1)%bl(i,2).ne.iaa(rs)%bl(i,2))) then
                print *,"bl",iaa(rs)%bl(i,:)
                print *,"bl_lim1",iaa_limits(rs,1)%bl(i,:)
                print *,"bl_lim2",iaa_limits(rs,2)%bl(i,:)
            end if 
        end do 
    end do 
    
    if (use_BOND(2).EQV..true.) then 
        ! Martin : this is essentially to have the option use_bond(2) (angle potential), be on only in the case of prolines.
        do rs=1,nseq 
            if (amino(seqtyp(rs)).eq.'PRO') cycle
                nrsbaeff(rs)=0
        end do 
    else if (use_BOND(2).EQV..false.) then 
        do rs=1,nseq 
            if (amino(seqtyp(rs)).eq.'PRO') then ! Martin I think this test exist somewhere else, but better safe than sorry 
                print *,"You must turn on the SC_BONDED_A option if you are going to have prolines in the sequence ! &
                &Stopping execution."
                call fexit()
            end if 
        end do 
    end if 
    
end subroutine


subroutine qcann_assign_bndtprms(limit)
    use params
    use sequen
    use atoms
    use inter
    use polypep
    
    
    integer limit
    integer i,j,rs
    
    iaa(:)=iaa_limits(:,limit) 
    
    ! Martin : these loops are to correct the array so that the order of the bonds are the same in both limits
    if (limit.eq.2) then 
        temp_save(:)=attyp(:)
        temp_bio_save(:)=b_type(:)
        do i=1,n 
            b_type(i)=transform_table(b_type(i))
            attyp(i)=bio_ljtyp(b_type(i))
        end do 
    else if (limit.eq.3) then 
        temp_save(:)=attyp(:)
        temp_bio_save(:)=b_type(:)
        
        ! Place holder for the biotype transform switch
        do i=1,nhis
            do j=1,at(his_state(i,1))%nsc
                b_type(at(his_state(i,1))%sc(j))=transform_table(his_eqv_table(b_type(at(his_state(i,1))%sc(j))))
            end do 
        end do 
    end if 

    call assign_bndtprms()
    
    iaa_limits(:,limit)=iaa(:)

    temp_nrsdi(:,limit)=nrsdieff(:)
    temp_nrsba(:,limit)=nrsbaeff(:)
    temp_nrsbl(:,limit)=nrsbleff(:)
    temp_nrsimpt(:,limit)=nrsimpteff(:)

    if (limit.ne.1) then 
        b_type(:)=temp_bio_save(:)
        attyp(:)=temp_save(:)
    end if 
    
end subroutine

subroutine qcann_absinth_atom(limit) 
    use params
    use sequen
    use atoms
    use energies
    use polypep
    
    integer limit
    integer i,kk,imol,j
    
    lj_eps(:,:)=lj_eps_limits(:,:,limit)
    lj_eps_14(:,:)=lj_eps_14_limits(:,:,limit)
    lj_sig(:,:)=lj_sig_limits(:,:,limit)
    lj_sig_14(:,:)=lj_sig_14_limits(:,:,limit)

    ! Martin : tranform all the needed quantities into the end state values
    
    if (limit.eq.2) then 
        temp_save(:)=attyp(:)
        temp_bio_save(:)=b_type(:)
        do i=1,n 
            b_type(i)=transform_table(b_type(i))
            attyp(i)=bio_ljtyp(b_type(i))
        end do 
    else if (limit.eq.3) then 
        temp_save(:)=attyp(:)
        temp_bio_save(:)=b_type(:)
        
        ! Place holder for the biotype transform switch
        do i=1,nhis
            do j=1,at(his_state(i,1))%nsc
                b_type(at(his_state(i,1))%sc(j))=transform_table(his_eqv_table(b_type(at(his_state(i,1))%sc(j))))
            end do 
        end do 
    end if 

    x(:)=x_limits(:,limit)
    y(:)=y_limits(:,limit)
    z(:)=z_limits(:,limit)
    
    atr(:)=atr_limits(:,limit)
    
    ! I do not need to do that with atvol and atbvol, because they are already in the function absinth atom
    call absinth_atom() 
    
    atsavred_limits(:,limit)=atsavred(:)
    atsavmaxfr_limits(:,limit)=atsavmaxfr(:)
    if (use_IMPSOLV.EQV..true.) then 
        atsavprm_limits(:,:,limit)=atsavprm(:,:)
    end if 
    atstv_limits(:,limit)=atstv(:)
    
    if (limit.ne.1) then 
        b_type(:)=temp_bio_save(:)
        attyp(:)=temp_save(:)
    end if 

end subroutine 

subroutine qcann_polar_group(limit)
  use atoms
  use sequen
  use molecule
  use cutoffs
  use polypep
  use inter
  use fyoc !martin added
  use params
  

  integer limit,rs,dp,alcsz,i,kk,imol,j
  
    iaa(:)=iaa_limits(:,limit)
    fudge(:)=fudge_limits(:,limit)

    
    if (limit.eq.2) then 
        temp_save(:)=attyp(:)
        temp_bio_save(:)=b_type(:)
        do i=1,n 
            b_type(i)=transform_table(b_type(i))
            attyp(i)=bio_ljtyp(b_type(i))
        end do 
    else if (limit.eq.3) then 
        temp_save(:)=attyp(:)
        temp_bio_save(:)=b_type(:)
        

        do i=1,nhis
            do j=1,at(his_state(i,1))%nsc
                b_type(at(his_state(i,1))%sc(j))=transform_table(his_eqv_table(b_type(at(his_state(i,1))%sc(j))))
                attyp(at(his_state(i,1))%sc(j))=bio_ljtyp(transform_table(his_eqv_table(b_type(at(his_state(i,1))%sc(j)))))
            end do 
        end do 
    end if 
    
    call polar_groups()
    
    kk=0
    do rs=1,nseq
        do i=1,at(rs)%ndpgrps
            if (at(rs)%dpgrp(i)%nc.ne.0) kk = kk + 1
        end do
    end do

    do rs=1,nseq
        allocate(ats_temp(rs,limit)%dpgrp(at(rs)%ndpgrps))
        ats_temp(rs,limit)%ndpgrps=at(rs)%ndpgrps
        do dp=1,at(rs)%ndpgrps
            allocate(ats_temp(rs,limit)%dpgrp(dp)%ats(at(rs)%dpgrp(dp)%nats))
            ats_temp(rs,limit)%dpgrp(dp)%nats=at(rs)%dpgrp(dp)%nats
            ats_temp(rs,limit)%dpgrp(dp)%cgn=at(rs)%dpgrp(dp)%cgn
            ats_temp(rs,limit)%dpgrp(dp)%ats(:)=at(rs)%dpgrp(dp)%ats(:)
        end do 
    end do 

    do imol=1,nmol ! Martin : this bit is taken from unbound's first subroutine
        do rs=rsmol(imol,1),rsmol(imol,2)
            if (natres(rs).gt.1) then
                allocate(iaa_limits(rs,limit)%expolin((natres(rs)*(natres(rs)-1))/2,2))
            end if 
            
            allocate(fudge_limits(rs,limit)%elin((at(rs)%npol*(at(rs)%npol-1))/2))
            allocate(iaa_limits(rs,limit)%polin((at(rs)%npol*(at(rs)%npol-1))/2,2))

            if (rs.lt.nseq) then
                alcsz = natres(rs)*natres(rs+1)
                if (disulf(rs).gt.0) alcsz = max(alcsz,natres(rs)*natres(disulf(rs)))

                allocate(iaa_limits(rs,limit)%expolnb(alcsz,2))
                
                alcsz = at(rs)%npol*at(rs+1)%npol
                if (disulf(rs).gt.0) alcsz = max(alcsz,at(rs)%npol*at(disulf(rs))%npol)

                allocate(fudge_limits(rs,limit)%elnb(alcsz))
                allocate(iaa_limits(rs,limit)%polnb(alcsz,2))
                
            end if
        end do
    end do
    
    temp_which_dpg(:,limit)=which_dpg(:)
    temp_nrpolnb(:,limit)=nrpolnb(:)
    temp_nrpolintra(:,limit)=nrpolintra(:)
    temp_nrexpolin(:,limit)=nrexpolin(:)
    temp_nrexpolnb(:,limit)=nrexpolnb(:)
        
    do rs=1,nseq
        if (rs.lt.nseq) then
            fudge_limits(rs,limit)%elnb(:)=fudge(rs)%elnb(:)
            iaa_limits(rs,limit)%polnb(:,:)=iaa(rs)%polnb(:,:)
            iaa_limits(rs,limit)%expolnb(:,:)=iaa(rs)%expolnb(:,:)
        end if 
        fudge_limits(rs,limit)%elin(:)=fudge(rs)%elin(:)
        iaa_limits(rs,limit)%polin(:,:)=iaa(rs)%polin(:,:)
        if (natres(rs).gt.1) then
            iaa_limits(rs,limit)%expolin(:,:)=iaa(rs)%expolin(:,:)
        end if 
    end do 
    if (limit.eq.2) then 
        atq_limits(:,3)=atq(:)
    else if (limit.eq.1) then 
        atq_limits(:,1)=atq(:)
    else if (limit.eq.3) then 
        atq_limits(:,4)=atq(:)
    end if 
    ! Return to the actual biotypes
    if (limit.ne.1) then 
        b_type(:)=temp_bio_save(:)
        do rs=1,nseq
            at(rs)%npol=0
            deallocate(fudge(rs)%elin)
            deallocate(iaa(rs)%polin)
            if (allocated(fudge(rs)%elnb)) then 
                deallocate(iaa(rs)%polnb)
                deallocate(fudge(rs)%elnb)
            end if 
            if (allocated(iaa(rs)%expolnb)) then 
                deallocate(iaa(rs)%expolnb)
            end if 
            if (natres(rs).gt.1) then
                deallocate(iaa(rs)%expolin)
            end if
            if (at(rs)%ndpgrps.eq.0) cycle
            deallocate(at(rs)%dpgrp)
        end do 
        deallocate(cglst%it)
        deallocate(cglst%nc)
        deallocate(cglst%tc)
    end if 
    if (limit.ne.1) then 
        b_type(:)=temp_bio_save(:)
        attyp(:)=temp_save(:)
    end if 
end subroutine

subroutine setup_charge_last_step(limit1,limit2)
    use inter
    use sequen
    use atoms
    use polypep
    use cutoffs
    use iounit
    
    integer, allocatable:: memory(:) ! Martin added 
    integer, allocatable:: polin_added1(:,:),polin_added2(:,:),polnb_added1(:,:),polnb_added2(:,:)!Martin added
    RTYPE, allocatable:: temp_var1(:),  temp_var2(:)
    integer, allocatable:: temp_var3(:)
    integer, allocatable:: temp_var5(:,:),temp_var6(:,:),temp_var7(:)
    integer limit1,limit2,rs,i,dp,ii,k,j,limit
    RTYPE temp_sum
    logical good
    !type(t_at), ALLOCATABLE :: ats_temp(:) ! Martin : added this 
    type(t_cglst):: cglst_temp
    integer lim2_tmp
    
!    Ok once here what we need to do is : 
    !   -find out which of the limits has the dummy atom.
    !   -lower the number of dipolar groups the limits that contains it, if necessary which rearanging the whichdp array
    !   -fix the intras list to contain none of the intra in dipolar group interaction
    !!!!!!WARNING : the algorithm here makes the assumption that the dummy atoms will ALWAYS
    !!!!!! be assigned to the last dipolar group. Should always be true, but may not be
    !!!!!! In essence it will take the atoms that is a dummy one and put it in the second to last dipolar group, then remove
    !!!!!! The last dipolar group from the list
    ! I should fix the dipolar groups first, then get the limits for the different qnm, tc,nc, 
    !Therefor I should not worry about deallocating aall the dpgrp
    ! I will however need to get a lmiit for this array, so I can replace with th appropriate one
    !will need temperary arrays for the folowing 
!        -Fixed arrays (won't be interpolated)
!        -at(rs)%dpgrp(dp1)%ats(i)
!        -at(rs)%dpgrp(dp1)%cgn ! number of atom per dipolar group
!        -at(rs)%dpgrp(dp1)%nats 
!        -at(rs)%dpgrp(dp1)%nats
!        
!        -Arrays to interpolate
!            -at(rs)%dpgrp(dp1)%nc
!             -at(rs)%dpgrp(dp1)%tc
!             -at(rs)%dpgrp(dp1)%qnm
! Additonal problem comes from the fact that the interpolation is in three and not two steps
!    or just for :
!    The logic behind this is going to be to check if a given AA has a dummy atom.
!    If it does, then take the value of the limit in which the atom is present for %ats,at%dpgrp%nats,at%ndpgrps
!    I think at that point, what I should do is only focus on %ats, %nats, maybe cgn, and from there on only, 
!    deduce the %nc%tc%qnm, and correct the which_dpgrp array 

    do rs=1,nseq
        good=.true.
        ! If there is no dpgrp in one but there is in the other, we have a disapearing ion. Take the highest number of dipolar group
        if (((ats_temp(rs,limit1)%ndpgrps.ne.ats_temp(rs,limit2)%ndpgrps)).and.&
        &((ats_temp(rs,limit1)%ndpgrps.eq.0).or.(ats_temp(rs,limit2)%ndpgrps.eq.0))) then
            if (ats_temp(rs,limit1)%ndpgrps.eq.0) then
                
                ats_temp(rs,limit1)%ndpgrps=ats_temp(rs,limit2)%ndpgrps
                ! Martin : I may need to come back to this 
!                deallocate(at(rs)%dpgrp)
!                allocate(at(rs)%dpgrp(at(rs)%ndpgrps))
                deallocate(ats_temp(rs,limit1)%dpgrp)
                allocate(ats_temp(rs,limit1)%dpgrp(ats_temp(rs,limit1)%ndpgrps))
                do dp=1,ats_temp(rs,limit1)%ndpgrps
                    ats_temp(rs,limit1)%dpgrp(dp)%nats=ats_temp(rs,limit2)%dpgrp(dp)%nats
                    deallocate(ats_temp(rs,limit1)%dpgrp(dp)%ats)
                    allocate(ats_temp(rs,limit1)%dpgrp(dp)%ats(ats_temp(rs,limit1)%dpgrp(dp)%nats))
                    ats_temp(rs,limit1)%dpgrp(dp)%ats(:)=ats_temp(rs,limit2)%dpgrp(dp)%ats(:)
                    ats_temp(rs,limit1)%dpgrp(dp)%cgn=ats_temp(rs,limit2)%dpgrp(dp)%cgn
                    !at(rs)%dpgrp(dp)%nc=ats_temp(rs)%dpgrp(dp)%nc
                end do
            else
                ats_temp(rs,limit2)%ndpgrps=ats_temp(rs,limit1)%ndpgrps
                ! Martin : I may need to come back to this 
!                deallocate(at(rs)%dpgrp)
!                allocate(at(rs)%dpgrp(at(rs)%ndpgrps))
                
                deallocate(ats_temp(rs,limit2)%dpgrp)
                allocate(ats_temp(rs,limit2)%dpgrp(ats_temp(rs,limit2)%ndpgrps))
                do dp=1,ats_temp(rs,limit2)%ndpgrps
                    ats_temp(rs,limit2)%dpgrp(dp)%nats=ats_temp(rs,limit1)%dpgrp(dp)%nats
!                    deallocate(ats_temp(rs,limit2)%dpgrp(dp)%ats)
                    allocate(ats_temp(rs,limit2)%dpgrp(dp)%ats(ats_temp(rs,limit2)%dpgrp(dp)%nats))
                    ats_temp(rs,limit2)%dpgrp(dp)%ats(:)=ats_temp(rs,limit1)%dpgrp(dp)%ats(:)
                    ats_temp(rs,limit2)%dpgrp(dp)%cgn=ats_temp(rs,limit1)%dpgrp(dp)%cgn
                    !at(rs)%dpgrp(dp)%nc=ats_temp(rs)%dpgrp(dp)%nc
                end do
            end if 
        end if 
        ! Ok that takes care of the ions
        ! Now for the dummy atoms
        ! If any of the atoms has a null radius for one limits, then takes the other FOR THE ENTIRE RESIDUE
        do dp=1,ats_temp(rs,limit1)%ndpgrps
            do i=1,ats_temp(rs,limit1)%dpgrp(dp)%nats
                ii=ats_temp(rs,limit1)%dpgrp(dp)%ats(i)
                if (abs(atr_limits(ii,limit1)).lt.1.0D-5) then
                    good=.false.
                end if 
            end do 
        end do 
!       if the second limits has the zero, then do nothing, cause its already good, because the actual array is already based on the 
!       first limit
        if (good.eqv..false.) then 
            do dp=1,ats_temp(rs,limit1)%ndpgrps
                deallocate(ats_temp(rs,limit1)%dpgrp(dp)%ats) 
            end do 
            ats_temp(rs,limit1)%ndpgrps=ats_temp(rs,limit2)%ndpgrps
            deallocate(ats_temp(rs,limit1)%dpgrp)
            allocate(ats_temp(rs,limit1)%dpgrp(ats_temp(rs,limit1)%ndpgrps))
            do dp=1,ats_temp(rs,limit1)%ndpgrps
                ats_temp(rs,limit1)%dpgrp(dp)%nats=ats_temp(rs,limit2)%dpgrp(dp)%nats
                allocate(ats_temp(rs,limit1)%dpgrp(dp)%ats(ats_temp(rs,limit1)%dpgrp(dp)%nats))
                ats_temp(rs,limit1)%dpgrp(dp)%ats(:)=ats_temp(rs,limit2)%dpgrp(dp)%ats(:)
                ats_temp(rs,limit1)%dpgrp(dp)%cgn=ats_temp(rs,limit2)%dpgrp(dp)%cgn
            end do 
        end if 
    !
    end do     
    !
    !deallocate(ats_temp)
    !
    ! Now to rebuild which_dpgrp
    do rs=1,nseq
        do dp=1,at(rs)%ndpgrps
            do i=1,at(rs)%dpgrp(dp)%nats
                which_dpg(at(rs)%dpgrp(dp)%ats(i))=dp
            end do 
        end do 
    end do     
    !
!        do i=1,n
!        if (which_dpg(i).ne.temp_which_dpg(i)) then 
!            print *,amino(seqtyp(atmres(i)))," is not correctly implemented. Exiting."
!            ! If you receive this message, it is most likely due to the dipolar group being different in each end states.
!            ! I saw no other option than to hardcode the dipolar group for those residues, which is done in get_dipgrps()
!            call fexit()
!        end if 
!    end do     
!             -Arrays to interpolate
!             -at(rs)%dpgrp(dp1)%nc
!             -at(rs)%dpgrp(dp1)%tc
!             -at(rs)%dpgrp(dp1)%qnm
    
    
    
!    !Martin : commmentned 15/6/20
!    do rs=1,nseq
!        do dp=1,at(rs)%ndpgrps
!            allocate(at(rs)%dpgrp(dp)%qnm_limits(4))
!            allocate(at(rs)%dpgrp(dp)%tc_limits(4))         
!!            ! Since the net charge is not interpolable (integer)), I will give the maximum value
!!            ! This is fine since it seems it is only used for cycling through loops
!        end do
!    end do
!    Now I know this is a bit confusiong, but we are going to change atq limits again. This time, all atoms that are part of a 
!    dipolar group that has at least one residue that changes start will be interpolated to zero in the middle state
! Needs to be done in two part because of the way I constructued stuff. Rewritting the all thing may ead to simpler expressions, but 
    ! This is only executed once at the beginning i.e. I don't really care, nit worth the time
  do rs=1,nseq
    do i=1,at(rs)%ndpgrps
        good=.true.
        ! This is for the his limits, has the same net charge among both tautomers
        at(rs)%dpgrp(i)%qnm_limits(4)=at(rs)%dpgrp(i)%qnm_limits(3)
        at(rs)%dpgrp(i)%tc_limits(4)=at(rs)%dpgrp(i)%tc_limits(3)
        
        do k=1,at(rs)%dpgrp(i)%nats
            if (atq_limits(at(rs)%dpgrp(i)%ats(k),1).ne.atq_limits(at(rs)%dpgrp(i)%ats(k),3)) then 
                good=.false.
            end if 
        end do 
        if (good.eqv..false.) then 
                atq_limits(at(rs)%dpgrp(i)%ats(:),2)=0
        end if 
    end do 
  end do 
    
  do rs=1,nseq
    do i=1,at(rs)%ndpgrps
      at(rs)%dpgrp(i)%qnm_limits(:) = 0.0
      at(rs)%dpgrp(i)%tc_limits(:)=0.0

      do k=1,at(rs)%dpgrp(i)%nats
          do j=1,3
            at(rs)%dpgrp(i)%qnm_limits(j) = at(rs)%dpgrp(i)%qnm_limits(j) + abs(atq_limits(at(rs)%dpgrp(i)%ats(k),j))
            at(rs)%dpgrp(i)%tc_limits(j) = at(rs)%dpgrp(i)%tc_limits(j) + atq_limits(at(rs)%dpgrp(i)%ats(k),j)
          end do 
          at(rs)%dpgrp(i)%nc=at(rs)%dpgrp(i)%nc+atq(at(rs)%dpgrp(i)%ats(k))
      end do
      ! In this next part, we assign the nc charge as being the the one of the end state with the highest absolute charge
      if (abs(at(rs)%dpgrp(i)%tc_limits(1)).le.abs(at(rs)%dpgrp(i)%tc_limits(2))) then
          at(rs)%dpgrp(i)%nc=nint(at(rs)%dpgrp(i)%tc_limits(2))
      else if  (abs(at(rs)%dpgrp(i)%tc_limits(2)).lt.abs(at(rs)%dpgrp(i)%tc_limits(1))) then
          at(rs)%dpgrp(i)%nc=nint(at(rs)%dpgrp(i)%tc_limits(1))
      end if 
    end do
  end do
! The cglst needs to be done after interpolation 
! So I need to put the monopole array into the 
! To make life simpler, and because the interpolation is ALWAYS going to be from the charged to the uncharged state
! We then know that the starting state array is going to have all the monopoles from the end state one, plus the titratable groups

! Martin : moved 15/6/20
! allocate(cglst%tc_limits(cglst%ncs,3))

  ! Way overallocated, but will on;ly be used for a short time
  allocate(polin_added1(nseq,1000))
  allocate(polin_added2(nseq,1000))
  allocate(polnb_added1(nseq,1000))
  allocate(polnb_added2(nseq,1000))

  polin_added1(:,:)=0
  polin_added2(:,:)=0
  polnb_added1(:,:)=0
  polnb_added2(:,:)=0
  
  do i=1,cglst%ncs
    temp_sum=0.
    do k=1,at(atmres(cglst%it(i)))%dpgrp(which_dpg(cglst%it(i)))%nats
        temp_sum=temp_sum+atq_limits(at(atmres(cglst%it(i)))%dpgrp(which_dpg(cglst%it(i)))%ats(k),1)
    end do 
    cglst%tc_limits(i,1)=temp_sum
    temp_sum=0.
    do k=1,at(atmres(cglst%it(i)))%dpgrp(which_dpg(cglst%it(i)))%nats
        temp_sum=temp_sum+atq_limits(at(atmres(cglst%it(i)))%dpgrp(which_dpg(cglst%it(i)))%ats(k),2)
    end do 
    cglst%tc_limits(i,2)=temp_sum
    temp_sum=0.
    do k=1,at(atmres(cglst%it(i)))%dpgrp(which_dpg(cglst%it(i)))%nats
        temp_sum=temp_sum+atq_limits(at(atmres(cglst%it(i)))%dpgrp(which_dpg(cglst%it(i)))%ats(k),3)
    end do 
    cglst%tc_limits(i,3)=temp_sum
  end do 

    ! Now the next part should unecessary 
    ! Martin : first get the size of the arrays to match
    do rs=1,nseq
        if (rs.lt.nseq) then
            if (size(iaa_limits(rs,1)%polnb(:,1)).ne.size(iaa_limits(rs,2)%polnb(:,1))) then 
                if (size(iaa_limits(rs,1)%polnb(:,1)).gt.size(iaa_limits(rs,2)%polnb(:,1))) then 
                    !!!!polnb
                    allocate(temp_var5(size(iaa_limits(rs,2)%polnb(:,1)),2))
                    temp_var5(:,:)=iaa_limits(rs,2)%polnb(:,:)
                    deallocate(iaa_limits(rs,2)%polnb)
                    allocate(iaa_limits(rs,2)%polnb(size(iaa_limits(rs,1)%polnb(:,1)),2))
                    iaa_limits(rs,2)%polnb(:,:)=0
                    do i=1,size(temp_var5(:,1))
                        iaa_limits(rs,2)%polnb(i,:)=temp_var5(i,:)
                    end do 
                    
                    !!!!!!!%expolnb
                    allocate(temp_var6(size(temp_var5(:,1)),2))
                    temp_var6(:,:)=iaa_limits(rs,2)%expolnb(:,:)
                    deallocate(iaa_limits(rs,2)%expolnb)
                    allocate(iaa_limits(rs,2)%expolnb(size(iaa_limits(rs,1)%expolnb(:,1)),2))
                    iaa_limits(rs,2)%expolnb(:,:)=0.
                    do i=1,size(temp_var5(:,1))
                        iaa_limits(rs,2)%expolnb(i,:)=temp_var5(i,:)
                    end do 
                    

                    !!!!!)%elnb
                    allocate(temp_var7(size(temp_var5(:,1))))
                    temp_var7(:)=fudge_limits(rs,2)%elnb(:)
                    
                    deallocate(fudge_limits(rs,2)%elnb)
                    allocate(fudge_limits(rs,2)%elnb(size(fudge_limits(rs,1)%elnb(:))))
                    
                    fudge_limits(rs,2)%elnb(:)=0.
                    do i=1,size(temp_var5(:,1))
                        fudge_limits(rs,2)%elnb(i)=temp_var7(i)
                    end do 
                    
                    deallocate(temp_var5)
                    deallocate(temp_var6)
                    deallocate(temp_var7)
                else 
                    !!!!polnb
                    allocate(temp_var5(size(iaa_limits(rs,1)%polnb(:,1)),2))
                    temp_var5(:,:)=iaa_limits(rs,1)%polnb(:,:)
                    deallocate(iaa_limits(rs,1)%polnb)
                    allocate(iaa_limits(rs,1)%polnb(size(iaa_limits(rs,2)%polnb(:,1)),2))
                    iaa_limits(rs,1)%polnb(:,:)=0
                    do i=1,size(temp_var5(:,1))
                        iaa_limits(rs,1)%polnb(i,:)=temp_var5(i,:)
                    end do 
                    !!!!!!!%expolnb
                    allocate(temp_var6(size(temp_var5(:,1)),2))
                    temp_var6(:,:)=iaa_limits(rs,1)%expolnb(:,:)
                    deallocate(iaa_limits(rs,1)%expolnb)
                    
                    allocate(iaa_limits(rs,1)%expolnb(size(iaa_limits(rs,2)%expolnb(:,1)),2))
                    iaa_limits(rs,1)%expolnb(:,:)=0.
                    do i=1,size(temp_var5(:,1))
                        iaa_limits(rs,1)%expolnb(i,:)=temp_var5(i,:)
                    end do 
                    !!!!!)%elnb
                    allocate(temp_var7(size(temp_var5(:,1))))
                    temp_var7(:)=fudge_limits(rs,1)%elnb(:)
                    deallocate(fudge_limits(rs,1)%elnb)
                    
                    allocate(fudge_limits(rs,1)%elnb(size(fudge_limits(rs,2)%elnb(:))))
                    fudge_limits(rs,1)%elnb(:)=0.
                    do i=1,size(temp_var5(:,1))
                        fudge_limits(rs,1)%elnb(i)=temp_var7(i)
                    end do 
                    
                    deallocate(temp_var5)
                    deallocate(temp_var6)
                    deallocate(temp_var7)
                end if 
            end if 
        end if 
        
    ! and now atin 
        
        if (size(iaa_limits(rs,1)%polin(:,1)).ne.size(iaa_limits(rs,2)%polin(:,1))) then 
                if (size(iaa_limits(rs,1)%polin(:,1)).gt.size(iaa_limits(rs,2)%polin(:,1))) then 
                    !!!!polin
                    allocate(temp_var5(size(iaa_limits(rs,2)%polin(:,1)),2))
                    temp_var5(:,:)=iaa_limits(rs,2)%polin(:,:)
                    deallocate(iaa_limits(rs,2)%polin)
                    allocate(iaa_limits(rs,2)%polin(size(iaa_limits(rs,1)%polin(:,1)),2))
                    iaa_limits(rs,2)%polin(:,:)=0
                    do i=1,size(temp_var5(:,1))
                        iaa_limits(rs,2)%polin(i,:)=temp_var5(i,:)
                    end do 
                    
                    !!!!!!!%expolin
                    if (natres(rs).gt.1) then
                        allocate(temp_var6(size(temp_var5(:,1)),2))
                        temp_var6(:,:)=iaa_limits(rs,2)%expolin(:,:)
                        deallocate(iaa_limits(rs,2)%expolin)
                        allocate(iaa_limits(rs,2)%expolin(size(iaa_limits(rs,1)%expolin(:,1)),2))
                        iaa_limits(rs,2)%expolin(:,:)=0.
                        do i=1,size(temp_var5(:,1))
                            iaa_limits(rs,2)%expolin(i,:)=temp_var5(i,:)
                        end do 
                    end if 

                    !!!!!)%elin
                    allocate(temp_var7(size(temp_var5(:,1))))
                    temp_var7(:)=fudge_limits(rs,2)%elin(:)
                    
                    deallocate(fudge_limits(rs,2)%elin)
                    allocate(fudge_limits(rs,2)%elin(size(fudge_limits(rs,1)%elin(:))))
                    
                    fudge_limits(rs,2)%elin(:)=0.
                    do i=1,size(temp_var5(:,1))
                        fudge_limits(rs,2)%elin(i)=temp_var7(i)
                    end do 

                    deallocate(temp_var5)
                    deallocate(temp_var6)
                    deallocate(temp_var7)
                else 
                    !!!!polin
                    allocate(temp_var5(size(iaa_limits(rs,1)%polin(:,1)),2))
                    temp_var5(:,:)=iaa_limits(rs,1)%polin(:,:)
                    deallocate(iaa_limits(rs,1)%polin)
                    allocate(iaa_limits(rs,1)%polin(size(iaa_limits(rs,2)%polin(:,1)),2))
                    iaa_limits(rs,1)%polin(:,:)=0
                    do i=1,size(temp_var5(:,1))
                        iaa_limits(rs,1)%polin(i,:)=temp_var5(i,:)
                    end do 
                    !!!!!!!%expolin
                    if (natres(rs).gt.1) then
                        allocate(temp_var6(size(temp_var5(:,1)),2))
                        temp_var6(:,:)=iaa_limits(rs,1)%expolin(:,:)
                        deallocate(iaa_limits(rs,1)%expolin)

                        allocate(iaa_limits(rs,1)%expolin(size(iaa_limits(rs,2)%expolin(:,1)),2))
                        iaa_limits(rs,1)%expolin(:,:)=0.
                        do i=1,size(temp_var5(:,1))
                            iaa_limits(rs,1)%expolin(i,:)=temp_var5(i,:)
                        end do 
                    end if 
                    !!!!!)%elin
                    allocate(temp_var7(size(temp_var5(:,1))))
                    temp_var7(:)=fudge_limits(rs,1)%elin(:)
                    deallocate(fudge_limits(rs,1)%elin)
                    
                    allocate(fudge_limits(rs,1)%elin(size(fudge_limits(rs,2)%elin(:))))
                    fudge_limits(rs,1)%elin(:)=0.
                    do i=1,size(temp_var5(:,1))
                        fudge_limits(rs,1)%elin(i)=temp_var7(i)
                    end do 
                    
                    deallocate(temp_var5)
                    deallocate(temp_var6)
                    deallocate(temp_var7)
                end if 
            end if 
    end do 
    ! Martin : ok so after that I still need to align the values
    allocate(temp_var3(2))
    allocate(temp_var2(1))
    allocate(temp_var7(2))
    allocate(memory(2))

!    write(ilog,*)"iaa_limits(3,1)%polnb(:,:)",nrpolnb(3)
!    write(ilog,*) iaa_limits(3,1)%polnb(:,:)
!    
!    write(ilog,*) "iaa_limits(3,2)%polnb(:,:)",nrpolnb(3)
!    write(ilog,*) iaa_limits(3,2)%polnb(:,:)
!    
!    write(ilog,*)"iaa_limits(3,1)%polin(:,:)",nrpolintra(3)
!    write(ilog,*) iaa_limits(3,1)%polin(:,:)
!    
!    write(ilog,*) "iaa_limits(3,2)%polin(:,:)",nrpolintra(3)
!    write(ilog,*) iaa_limits(3,2)%polin(:,:)
    
    do rs=1,nseq
        if (rs.ge.nseq) cycle 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! nrpolnb
        do i=1,nrpolnb(rs) ! This part is to get the right order
            
            good=.true.
            memory(:)=test_combi(iaa_limits(rs,1)%polnb(i,:),iaa_limits(rs,2)%polnb(i,:),size(iaa(rs)%polnb(i,:)))            

            do k=1,size(memory)
                if (memory(k).eq.0) then 
                    good=.false.
                end if 
            end do 
            if (good.eqv..false.) then 
             ! So if we find a residue array that is not strictly indentical, we go in the reste of the array to make a permutation 
                if (temp_nrpolnb(rs,limit1).le.temp_nrpolnb(rs,limit2)) then !That is if the first array is shorter
                    do j=1,nrpolnb(rs)-i
                        good=.true.
                        memory(:)=test_combi(iaa_limits(rs,1)%polnb(i+j,:),iaa_limits(rs,2)%polnb(i,:),size(iaa(rs)%polnb(i,:)))
                        do k=1,size(memory)
                            if (memory(k).eq.0) then 
                                good=.false.
                            end if 
                        end do 

                        if (good.eqv..true.) then                       
                            temp_var3(:)=iaa_limits(rs,1)%polnb(i,:)
                            iaa_limits(rs,1)%polnb(i,:)=iaa_limits(rs,1)%polnb(i+j,:)
                            iaa_limits(rs,1)%polnb(i+j,:)=temp_var3(:)

                            temp_var2(1)=fudge_limits(rs,1)%elnb(i)
                            fudge_limits(rs,1)%elnb(i)=fudge_limits(rs,1)%elnb(i+j)
                            fudge_limits(rs,1)%elnb(i+j)=temp_var2(1)
                            
                            ! Even though the expoling seems to only be used in the force using routines,
                            ! I will do it for consistency's sake
                            temp_var7(:)=iaa_limits(rs,1)%expolnb(i,:)
                            iaa_limits(rs,1)%expolnb(i,:)=iaa_limits(rs,1)%expolnb(i+j,:)
                            iaa_limits(rs,1)%expolnb(i+j,:)=temp_var7(:)
                            
                            exit 
                        end if   
                    end do 
                    ! If you did not find it, then assign the other limit values, and remember to turn
                    ! Fudge to 0 in a latter step 
                    if (good.eqv..false.) then !that is supposed to be if you didn't find it 
                        polnb_added1(rs,i)=1
                        iaa_limits(rs,2)%polnb(i,:)=iaa_limits(rs,1)%polnb(i,:)
                    else 
                        polnb_added1(rs,i)=0
                    end if 
                else    
                    do j=1,nrpolnb(rs)-i
                        good=.true.
                        memory(:)=test_combi(iaa_limits(rs,2)%polnb(i+j,:),iaa_limits(rs,1)%polnb(i,:),2)
                        do k=1,size(memory)
                            if (memory(k).eq.0) then 
                                good=.false.
                            end if 
                        end do 

                        if (good.eqv..true.) then                       
                            temp_var3(:)=iaa_limits(rs,2)%polnb(i,:)
                            iaa_limits(rs,2)%polnb(i,:)=iaa_limits(rs,2)%polnb(i+j,:)
                            iaa_limits(rs,2)%polnb(i+j,:)=temp_var3(:)

                            temp_var2(1)=fudge_limits(rs,2)%elnb(i)
                            fudge_limits(rs,2)%elnb(i)=fudge_limits(rs,2)%elnb(i+j)
                            fudge_limits(rs,2)%elnb(i+j)=temp_var2(1)
                            
                            temp_var7(:)=iaa_limits(rs,2)%expolnb(i,:)
                            iaa_limits(rs,2)%expolnb(i,:)=iaa_limits(rs,2)%expolnb(i+j,:)
                            iaa_limits(rs,2)%expolnb(i+j,:)=temp_var7(:)
                            exit 
                        end if   
                    end do 
                    ! If you did not find it, then assign the other limit values, and remember to turn
                    ! Fudge to 0 in a latter step 
                    if (good.eqv..false.) then !that is supposed to be if you didn't find it 
                        polnb_added2(rs,i)=1
                        iaa_limits(rs,2)%polnb(i,:)=iaa_limits(rs,1)%polnb(i,:)
                    else 
                        polnb_added2(rs,i)=0
                    end if 
                end if 
            end if  
            
        end do 

    end do 
    deallocate(memory)
    
    do rs=1,nseq
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! nrpolin
        do i=1,nrpolintra(rs) ! This part is to get the right order
            allocate(memory(size(iaa(rs)%polin(i,:))))
            good=.true.
            memory(:)=test_combi(iaa_limits(rs,1)%polin(i,:),iaa_limits(rs,2)%polin(i,:),size(iaa(rs)%polin(i,:)))
            do k=1,size(memory)
                if (memory(k).eq.0) then
                    good=.false.
                end if 
            end do 
            if (good.eqv..false.) then 
            ! So if we find a residue array that is not strictly identical, we go in the rest of the array to make a permutation 
                if (temp_nrpolintra(rs,limit1).le.temp_nrpolintra(rs,limit2)) then ! That is if the first array is shorter
                    do j=1,nrpolintra(rs)-i
                        good=.true.
                        memory(:)=test_combi(iaa_limits(rs,1)%polin(i+j,:),iaa_limits(rs,2)%polin(i,:),size(iaa(rs)%polin(i,:)))
                        do k=1,size(memory)
                            if (memory(k).eq.0) then 
                                good=.false.
                            end if 
                        end do 
                        if (good.eqv..true.) then                       
                            temp_var3(:)=iaa_limits(rs,1)%polin(i,:)
                            iaa_limits(rs,1)%polin(i,:)=iaa_limits(rs,1)%polin(i+j,:)
                            iaa_limits(rs,1)%polin(i+j,:)=temp_var3(:)

                            temp_var2(1)=fudge_limits(rs,1)%elin(i)
                            fudge_limits(rs,1)%elin(i)=fudge_limits(rs,1)%elin(i+j)
                            fudge_limits(rs,1)%elin(i+j)=temp_var2(1)
                            
                            temp_var7(:)=iaa_limits(rs,1)%expolin(i,:)
                            iaa_limits(rs,1)%expolin(i,:)=iaa_limits(rs,1)%expolin(i+j,:)
                            iaa_limits(rs,1)%expolin(i+j,:)=temp_var7(:)
                            exit 
                        end if   
                    end do 
                    ! If you did not find it, then assign the other limit values, and remember to turn
                    ! Fudge to 0 in a latter step 
                    if (good.eqv..false.) then 
                        polin_added1(rs,i)=1  
                        iaa_limits(rs,1)%polin(i,:)=iaa_limits(rs,2)%polin(i,:)
                    else 
                        polin_added1(rs,i)=0
                    end if 
                else
                    do j=1,nrpolintra(rs)-i
                        good=.true.
                        memory(:)=test_combi(iaa_limits(rs,2)%polin(i+j,:),iaa_limits(rs,1)%polin(i,:),size(iaa(rs)%polin(i,:)))
                        do k=1,size(memory)
                            if (memory(k).eq.0) then 
                                good=.false.
                            end if 
                        end do 

                        if (good.eqv..true.) then                       
                            temp_var3(:)=iaa_limits(rs,2)%polin(i,:)
                            iaa_limits(rs,2)%polin(i,:)=iaa_limits(rs,2)%polin(i+j,:)
                            iaa_limits(rs,2)%polin(i+j,:)=temp_var3(:)

                            temp_var2(1)=fudge_limits(rs,2)%elin(i)
                            fudge_limits(rs,2)%elin(i)=fudge_limits(rs,2)%elin(i+j)
                            fudge_limits(rs,2)%elin(i+j)=temp_var2(1)

                            temp_var7(:)=iaa_limits(rs,2)%expolin(i,:)
                            iaa_limits(rs,2)%expolin(i,:)=iaa_limits(rs,2)%expolin(i+j,:)
                            iaa_limits(rs,2)%expolin(i+j,:)=temp_var7(:)
                            exit 
                        end if   
                    end do 
                    ! If you did not find it, then assign the other limit values, and remember to turn
                    ! Fudge to 0 in a latter step 
                    if (good.eqv..false.) then !that is supposed to be if you didn't find it 
                        polin_added2(rs,i)=1
                        iaa_limits(rs,2)%polin(i,:)=iaa_limits(rs,1)%polin(i,:)
                    else 
                        polin_added2(rs,i)=0
                    end if 
                end if 
            end if  
            deallocate(memory)
        end do 
    end do 
    
    deallocate(temp_var3)
    deallocate(temp_var2)
    deallocate(temp_var7)
    lim2_tmp=limit2+1
    do i=1,n
        rs=atmres(i)! Actually the (atq_limits(i,1).eq.0) is unuseful, since charges always disapear
        if ((temp_which_dpg(i,limit2).ne.temp_which_dpg(i,limit1)).and.&
        &((atq_limits(i,1).eq.0).or.(atq_limits(i,lim2_tmp).eq.0))) then 
            if (temp_which_dpg(i,limit1).eq.0) then! That only happens for disapearing ions
                
                temp_which_dpg(i,limit1)=1
                at(rs)%ndpgrps=at(rs)%ndpgrps+1
                
                allocate(at(rs)%dpgrp(which_dpg(i)))
                at(rs)%dpgrp(temp_which_dpg(i,limit1))%nats=1
                at(rs)%dpgrp(temp_which_dpg(i,limit1))%ats(at(rs)%dpgrp(temp_which_dpg(i,limit1))%nats)=i
                
            else if (temp_which_dpg(i,limit2).eq.0) then! That only happens for disapearing ions
                cycle
                
            else if (temp_which_dpg(i,limit1).gt.temp_which_dpg(i,limit2)) then ! Ok so some assumption are wrong. But, let's automaticaly assign the 
                ! dummy atoms to whatever groups it is directly connected to
                ! Assumptions here are : the dummy atom is part of the last dipolar group of the residue
                ! That the dummy atoms is connected to only one atom
                
                do j=1,nrsbl(rs)! So here : try to find the atom that it is directly connected to 
                    if ((i.eq.iaa(rs)%bl(j,1))) then 
                        k=iaa(rs)%bl(j,2)
                        exit
                    else if ((i.eq.iaa(rs)%bl(j,2))) then 
                        k=iaa(rs)%bl(j,1)
                        exit
                    end if 
                end do     
                
                !which_dpg(i)=which_dpg(k)
                temp_which_dpg(i,limit1)=temp_which_dpg(k,limit1)
                
                at(rs)%ndpgrps=at(rs)%ndpgrps-1
                
                allocate(temp_var7(at(rs)%dpgrp(temp_which_dpg(i,limit1))%nats))
                temp_var7(:)=at(rs)%dpgrp(temp_which_dpg(i,limit1))%ats(:) ! Save all the atoms from the group
   
                
                at(rs)%dpgrp(temp_which_dpg(i,limit1))%nats=at(rs)%dpgrp(temp_which_dpg(i,limit1))%nats+1
                
                deallocate(at(rs)%dpgrp(temp_which_dpg(i,limit1))%ats)
                allocate(at(rs)%dpgrp(temp_which_dpg(i,limit1))%ats(at(rs)%dpgrp(temp_which_dpg(i,limit1))%nats))
                
                
                do j=1,at(rs)%dpgrp(temp_which_dpg(i,limit1))%nats-1
                    at(rs)%dpgrp(temp_which_dpg(i,limit1))%ats(j)=temp_var7(j)
                end do
                at(rs)%dpgrp(temp_which_dpg(i,limit1))%ats(at(rs)%dpgrp(temp_which_dpg(i,limit1))%nats)=i
                deallocate(temp_var7)
            end if 
        end if 
    end do
    
!!Martin : commented out to see if it solves problems
!    allocate(memory_1_4(n,n))
!    memory_1_4(:,:)=0
!
!    call get_all_1_4(memory_1_4)
!    ! Now I have to go through all the dipolar groups and identify pairs those that have one or more 
!    ! more atoms that are connected with less than 3 bonds
!    allocate(dpgrp_check(nseq,10,10)) ! Martin : overallocated
!    dpgrp_check(:,:,:)=1
!    do rs=1,nseq
!        do dp1=1,at(rs)%ndpgrps
!          do dp2=dp1+1,at(rs)%ndpgrps
!            do i=1,at(rs)%dpgrp(dp1)%nats
!              do k=1,at(rs)%dpgrp(dp2)%nats
!                  at1 = at(rs)%dpgrp(dp1)%ats(i)
!                  at2 = at(rs)%dpgrp(dp2)%ats(k)
!                  if ((memory_1_4(at1,at2).lt.3)) then !.and.(dpgrp_check(rs,dp1,dp2).eq.1)) then 
!                    dpgrp_check(rs,dp1,dp2)=0 ! is zero if any of the atoms in the dpgrp is less than 4 bonds away
!                    dpgrp_check(rs,dp2,dp1)=0
!                    ! Which means 0 if the dipolar group should not interact
!                  end if 
!              end do 
!            end do 
!          end do 
!        end do    
!    end do 
!    
!!    !Martin : Not sure what that what suppose to do, but it definitly causes problems
!    do rs=1,nseq
!        sub=0
!        do i=1,nrpolintra(rs)
!            at1=iaa_limits(rs,1)%polin(i,1)
!            at2=iaa_limits(rs,1)%polin(i,2)
!            if ((atmres(at1).eq.atmres(at2)).and.(dpgrp_check(rs,which_dpg(at1),which_dpg(at2)).eq.0)) then 
!                sub=sub+1
!                do j=i,nrpolintra(rs)-1
!                    iaa_limits(rs,1)%polin(j,:)=iaa_limits(rs,1)%polin(j+1,:)
!                    fudge_limits(rs,1)%elin(j)=fudge_limits(rs,1)%elin(j+1)
!                    iaa_limits(rs,1)%expolin(j,:)=iaa_limits(rs,1)%expolin(j+1,:)
!
!                    iaa_limits(rs,2)%polin(j,:)=iaa_limits(rs,2)%polin(j+1,:)
!                    fudge_limits(rs,2)%elin(j)=fudge_limits(rs,2)%elin(j+1)
!                    iaa_limits(rs,2)%expolin(j,:)=iaa_limits(rs,2)%expolin(j+1,:)
!
!                end do
!            end if
!        end do
!        nrpolintra(rs)=nrpolintra(rs)-sub
!    end do
! Ok so now, that the array have been aligned, and rectified in size and content, I still need to : 
!-make sure that the arrays are consistent (since the none limit array is not allign, I can just reasssign)
!-make sure none of the fudges have been set to zero, unless it is because the interaction only existed in one of the end states


    do rs=1,nseq
        do i=1,nrpolnb(rs)
            if (polnb_added1(rs,i).eq.1) then 
                fudge_limits(rs,1)%elnb(i)=0.0
            else if (polnb_added1(rs,i).eq.0) then 
                fudge_limits(rs,1)%elnb(i)=1.0
            else 
                print *,rs,i,polnb_added1(rs,i)
                print *,"This should not be happening under any circumstance (nb1)"
                call fexit()
            end if 
        end do 
        
        do i=1,nrpolnb(rs)
            if (polnb_added2(rs,i).eq.1) then 
                fudge_limits(rs,2)%elnb(i)=0.0
            else if (polnb_added2(rs,i).eq.0) then 
                fudge_limits(rs,2)%elnb(i)=1.0
            else 
                print *,rs,i,polnb_added2(rs,i)
                print *,"This should not be happening under any circumstance (nb2)"
                call fexit()
            end if 
        end do 
        
        do i=1,nrpolintra(rs)
            if (polin_added1(rs,i).eq.1) then 
                fudge_limits(rs,1)%elin(i)=0.0
            else if (polin_added1(rs,i).eq.0) then 
                fudge_limits(rs,1)%elin(i)=1.0
            else 
                print *,rs,i,polin_added1(rs,i)
                print *,"This should not be happening under any circumstance (in1)"
                call fexit()
            end if 
        end do 
        
        do i=1,nrpolintra(rs)
            if (polin_added2(rs,i).eq.1) then 
                fudge_limits(rs,2)%elin(i)=0.0
            else if (polin_added2(rs,i).eq.0) then  
                fudge_limits(rs,2)%elin(i)=1.0
            else 
                print *,rs,i,polin_added2(rs,i)
                print *,"This should not be happening under any circumstance (in2)"
                call fexit()
            end if 
        end do 
        
!        write(ilog,*) "nrpolintra_lim1",rs,nrpolintra(rs)
!        do i=1,nrpolintra(rs)
!            write(ilog,*) iaa_limits(rs,1)%polin(i,:),fudge_limits(rs,1)%elin(i)
!        end do 
!        
!        write(ilog,*) "nrpolintra_lim2",rs,nrpolintra(rs)
!        do i=1,nrpolintra(rs)
!            write(ilog,*) iaa_limits(rs,2)%polin(i,:),fudge_limits(rs,2)%elin(i)
!        end do 
!        
!        write(ilog,*) "nrpolnb_lim1",rs,nrpolnb(rs)
!        do i=1,nrpolnb(rs)
!            write(ilog,*) iaa_limits(rs,1)%polnb(i,:),fudge_limits(rs,1)%elnb(i)
!        end do 
!        
!        write(ilog,*) "nrpolnb_lim2",rs,nrpolnb(rs)
!        do i=1,nrpolnb(rs)
!            write(ilog,*) iaa_limits(rs,2)%polnb(i,:),fudge_limits(rs,2)%elnb(i)
!        end do 
!        
!        write(ilog,*) "polnb_added",rs
!        do i=1,nrpolnb(rs)
!            write(ilog,*) polnb_added1(rs,i),polnb_added2(rs,i)
!        end do 
!        
!        write(ilog,*) "polin_added",rs
!        do i=1,nrpolintra(rs)
!            write(ilog,*) polin_added1(rs,i),polin_added2(rs,i)
!        end do 

    end do 
    
    deallocate(polin_added1)
    deallocate(polin_added2)
    deallocate(polnb_added1)
    deallocate(polnb_added2)
    
    iaa(:)=iaa_limits(:,1) 
    fudge(:)=fudge_limits(:,1) 
    
end subroutine 
subroutine get_all_1_4(memory_1_4) 
    ! The goal here is to create a array of n*n that contains binary information on whether two atoms are
    ! separated by four bonds or less
    ! the goal is to have the number of bonds separating each atom from each other atom as a value in the array
    
    use inter
    use atoms
    use sequen
    
    implicit none
    integer at1,at2,rs,i,j,k,l
    integer, intent(inout) :: memory_1_4(:,:)

    ! So I can try to just put the distance in the index.
    
    do rs=1,nseq! First look at all the first neighbors
        do i=1,nrsbl(rs)
            at1=iaa(rs)%bl(i,1)
            at2=iaa(rs)%bl(i,2)
            memory_1_4(at1,at2)=1
            memory_1_4(at2,at1)=1
        end do  
    end do 

    
    do l=1,n ! Okay, that may be an overkill, but it must work in that case
        do i=1,n
            do j=1,n
                if ((memory_1_4(i,j).ne.1)) then ! If you haven't found a first neighbor contact
                    do k=1,n
                        if (((memory_1_4(i,k).ne.0).and.(memory_1_4(j,k).ne.0)).and.&
                        &(((memory_1_4(i,k)+memory_1_4(j,k)).lt.memory_1_4(i,j)).or.(memory_1_4(i,j).eq.0))) then
                            memory_1_4(i,j)=memory_1_4(i,k)+memory_1_4(j,k)
                            memory_1_4(j,i)=memory_1_4(i,k)+memory_1_4(j,k)
                        end if 
                    end do 
                end if 
            end do 
        end do 
    end do 

end subroutine 




subroutine sort_2_at_lists(input_list) 
    ! The goal here is to create a array of n*n that contains binary information on whether two atoms are
    ! separated by four bonds or less
    ! the goal is to have the number of bonds separating each atom from each other atom as a value in the array
    
    use inter
    use atoms
    use sequen
    
    implicit none
    integer at1,at2,rs,i,j,k,l,at
    integer temp,temp_arr(2)
    integer, intent(inout) :: input_list(:,:)
    
    !First, put larger number on the left column
    do at=1,size(input_list(:,1))
        if (input_list(at,1).gt.input_list(at,2)) then 
            temp =input_list(at,2)
            
            input_list(at,2)=input_list(at,1)
            input_list(at,1)=temp 
            
        end if 
    end do 
    
    !Then sort in order
    do at1=1,size(input_list(:,1))
        do at2=at1,size(input_list(:,1))
            
            if((input_list(at1,1).eq.0).or.(input_list(at1,2).eq.0).or.(input_list(at2,1).eq.0).or.(input_list(at2,2).eq.0)) cycle
            
            
            if (input_list(at2,1).lt.input_list(at1,1)) then 
                    temp_arr=input_list(at1,:)
                    input_list(at1,:)=input_list(at2,:)
                    input_list(at2,:)=temp_arr
                    
            else if (input_list(at2,1).eq.input_list(at1,1)) then 
                if (input_list(at2,2).lt.input_list(at1,2)) then 
                    temp=input_list(at1,2)
                    input_list(at1,2)=input_list(at2,2)
                    input_list(at2,2)=temp
                end if 
            end if
        end do 
    end do 

end subroutine 




subroutine print_iz()
    use zmatrix
    use mcsums
    
    
    implicit none
    integer i,j,k
    character WWW*200000,WW_tmp*1000,my_format*10
    
    
    do i=1,size(iz(1,:))
        WW_tmp=''
        write(WW_tmp,*) iz(1,i),char(9),iz(2,i),char(9),iz(3,i),char(9),iz(4,i),char(10)
        WWW=trim(WWW)//trim(WW_tmp)
!    
    end do 
!    
    write(my_format,'(A2,I0,A1)')'(A',len(trim(WWW)),')'
    write(val_det,my_format) trim(WWW)
    
    
end subroutine 

  subroutine write_detailed_energy_1(en_fos_memory,istep)! I may have to pas the istep too 
        
        use mcsums
        
        implicit none
        integer i,j,k
        integer, intent(in) :: istep
        character WWW*200000,WW_tmp*1000,my_format*10
        character*1 tab
        RTYPE, intent(in) :: en_fos_memory(:,:)
        WWW=''
        tab = char(9)


        do i=1 , size(en_fos_memory,dim=1)
            WW_tmp=''
            write(WW_tmp,'(I0,A1)') istep,tab
            WWW=trim(WWW)//trim(WW_tmp)

            do j=1 , size(en_fos_memory,dim=2)
                WW_tmp=''
                write(WW_tmp,'(A4)') 'res_'
                WWW=trim(WWW)//trim(WW_tmp)
                WW_tmp=''
                write(WW_tmp,'(I0,A1)') i,tab
                WWW=trim(WWW)//trim(WW_tmp)
            
                WW_tmp=''
                write(WW_tmp,'(A4,I0,A1)') "grp_",j,tab
                WWW=trim(WWW)//trim(WW_tmp)
                
                WW_tmp=''
                write(WW_tmp,*) en_fos_memory(i,j)
                WWW=trim(WWW)//trim(WW_tmp)
                    
                WWW=trim(WWW)//char(10)
            end do 

        end do

        write(my_format,'(A2,I0,A1)')'(A',len(trim(WWW)),')'
        if (len(trim(WWW)).ge.200000) then 
            print *,'The length of the output is higher than the container'
            stop
        end if 
        write(isav_det,my_format) trim(WWW)

  end subroutine write_detailed_energy_1
  
  subroutine write_detailed_energy_2(en_fos_memory,istep)! I may have to pas the istep too 
        
        use mcsums
        use sequen
        use polypep
        use atoms
        implicit none
        integer i,j,k,rs
        integer, intent(in) :: istep
        character WWW*200000,WW_tmp*1000,my_format*10
        character*1 tab
        RTYPE, intent(in) :: en_fos_memory(:,:,:)
        WWW=''
        tab = char(9)
        if (done_dets.eqv..false.) then
            done_dets=.true.
            WW_tmp='This is mapping to the atoms reflecting the &
            &structure of the file SAV details file. Below is the atomic number for&
            & each entry as a function of their resiue and solcation group. Second &
            &is the weights for each atoms used in the simulaton,&
            & third is the total volume of the hydration layer &
            &(ignoring overlaps), and fourth is the maximum fraction of the layer &
            & sovlation volume that does not overlap with neighboring atom, &
            &and finally the effective temperature depedent (if applies)&
            & free energy of sovlation for the residue (for each atom, redundant).&
            & Thus multiplying each of these values with the &
            & corresponding values in the SAV details file will give you the&
            & effective free energy of solvation contribution for the &
            & correspoonding atom, in kcal/mol.'

            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

            WW_tmp='Atomic numbers'
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)


            !first, print the corresponding atoms 
            do rs=1,nseq
                do i=1,at(rs)%nfosgrps
                    WW_tmp=''
                    write(WW_tmp,'(A4)') 'res_'
                    WWW=trim(WWW)//trim(WW_tmp)
                    WW_tmp=''
                    write(WW_tmp,'(I0,A1)') rs,tab
                    WWW=trim(WWW)//trim(WW_tmp)

                    WW_tmp=''
                    write(WW_tmp,'(A4,I0,A1)') "grp_",i,tab
                    WWW=trim(WWW)//trim(WW_tmp)

                    do j=1,at(rs)%fosgrp(i)%nats
                        WW_tmp=''
                        write(WW_tmp,'(I0,A1)') ,at(rs)%fosgrp(i)%ats(j),tab
                        WWW=trim(WWW)//trim(WW_tmp)
                    end do 
                    WWW=trim(WWW)//char(10)
                end do 

            end do 


            WW_tmp='Atomic FOS weights'
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
            do rs=1,nseq
                do i=1,at(rs)%nfosgrps
                    WW_tmp=''
                    write(WW_tmp,'(A4)') 'res_'
                    WWW=trim(WWW)//trim(WW_tmp)
                    WW_tmp=''
                    write(WW_tmp,'(I0,A1)') rs,tab
                    WWW=trim(WWW)//trim(WW_tmp)

                    WW_tmp=''
                    write(WW_tmp,'(A4,I0,A1)') "grp_",i,tab
                    WWW=trim(WWW)//trim(WW_tmp)

                    do j=1,at(rs)%fosgrp(i)%nats
                        WW_tmp=''
                        write(WW_tmp,*) at(rs)%fosgrp(i)%wts(j),tab
                        WWW=trim(WWW)//trim(WW_tmp)
                    end do 
                    WWW=trim(WWW)//char(10)
                end do 
            end do 

            WW_tmp='Atomic solvation layers volume'
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

            do rs=1,nseq
                do i=1,at(rs)%nfosgrps
                    WW_tmp=''
                    write(WW_tmp,'(A4)') 'res_'
                    WWW=trim(WWW)//trim(WW_tmp)
                    WW_tmp=''
                    write(WW_tmp,'(I0,A1)') rs,tab
                    WWW=trim(WWW)//trim(WW_tmp)

                    WW_tmp=''
                    write(WW_tmp,'(A4,I0,A1)') "grp_",i,tab
                    WWW=trim(WWW)//trim(WW_tmp)

                    do j=1,at(rs)%fosgrp(i)%nats
                        WW_tmp=''
                        write(WW_tmp,*) atbvol(at(rs)%fosgrp(i)%ats(j)),tab
                        WWW=trim(WWW)//trim(WW_tmp)
                    end do 
                    WWW=trim(WWW)//char(10)
                end do 
            end do 

            WW_tmp= 'Atomic maximum fraction of accessible solvation volume'
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

            do rs=1,nseq
                do i=1,at(rs)%nfosgrps
                    WW_tmp=''
                    write(WW_tmp,'(A4)') 'res_'
                    WWW=trim(WWW)//trim(WW_tmp)
                    WW_tmp=''
                    write(WW_tmp,'(I0,A1)') rs,tab
                    WWW=trim(WWW)//trim(WW_tmp)

                    WW_tmp=''
                    write(WW_tmp,'(A4,I0,A1)') "grp_",i,tab
                    WWW=trim(WWW)//trim(WW_tmp)

                    do j=1,at(rs)%fosgrp(i)%nats
                        WW_tmp=''
                        write(WW_tmp,*) atsavmaxfr(at(rs)%fosgrp(i)%ats(j)),tab
                        WWW=trim(WWW)//trim(WW_tmp)
                    end do 
                    WWW=trim(WWW)//char(10)
                end do 
            end do 

            WW_tmp= 'Residue temperature dependent reference FOS'
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

            do rs=1,nseq
                do i=1,at(rs)%nfosgrps
                    WW_tmp=''
                    write(WW_tmp,'(A4)') 'res_'
                    WWW=trim(WWW)//trim(WW_tmp)
                    WW_tmp=''
                    write(WW_tmp,'(I0,A1)') rs,tab
                    WWW=trim(WWW)//trim(WW_tmp)

                    WW_tmp=''
                    write(WW_tmp,'(A4,I0,A1)') "grp_",i,tab
                    WWW=trim(WWW)//trim(WW_tmp)

                    do j=1,at(rs)%fosgrp(i)%nats
                        WW_tmp=''
                        write(WW_tmp,*) fos_Tdep_duplicate(at(rs)%fosgrp(i)%val(1:3)),tab
                        WWW=trim(WWW)//trim(WW_tmp)
                    end do 
                    WWW=trim(WWW)//char(10)
                end do 
            end do 

            write(my_format,'(A2,I0,A1)')'(A',len(trim(WWW)),')'

            if (len(trim(WWW)).ge.200000) then 
                print *,'The length of the output is higher than the container'
                stop
            end if 
            write(isav_det_map,my_format) trim(WWW)
        end if 
        
        
        
        WWW=''
        WW_tmp=''
        write(WW_tmp,'(I0,A1)') istep,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
        do i=1 ,nseq !  size(en_fos_memory,dim=1)
            do j=1 , at(i)%nfosgrps !size(en_fos_memory,dim=2)
                WW_tmp=''
                write(WW_tmp,'(A4)') 'res_'
                WWW=trim(WWW)//trim(WW_tmp)
                WW_tmp=''
                write(WW_tmp,'(I0,A1)') i,tab
                WWW=trim(WWW)//trim(WW_tmp)
            
                WW_tmp=''
                write(WW_tmp,'(A4,I0,A1)') "grp_",j,tab
                WWW=trim(WWW)//trim(WW_tmp)
                
                do k=1,at(i)%fosgrp(j)%nats!size(en_fos_memory,dim=3)
                    WW_tmp=''
                    write(WW_tmp,*) en_fos_memory(i,j,k)
                    WWW=trim(WWW)//trim(WW_tmp)
                    
                end do 
                WWW=trim(WWW)//char(10)
            end do 

        end do

        write(my_format,'(A2,I0,A1)')'(A',len(trim(WWW)),')'
        if (len(trim(WWW)).ge.200000) then 
            print *,'The length of the output is higher than the container'
            stop
        end if 
        write(isav_det,my_format) trim(WWW)
            

  end subroutine write_detailed_energy_2
  
  
  
  subroutine write_pka_details(istep)

    use molecule
    use polypep
    use atoms
    use zmatrix
    use params
    use mcsums
    use atoms !martin 
    use inter
    
    implicit none 
    integer, intent(in) :: istep
    integer i,j, iqg, icat,k
    
    character WWW*200000, WW_tmp*1000 ,my_format*10
    character WWW1*200000,WW_tmp1*1000,my_format1*10
    character WWW2*200000,WW_tmp2*1000,my_format2*10
    character WWW3*200000,WW_tmp3*1000,my_format3*10
    character WWW4*200000,WW_tmp4*1000,my_format4*10
    character WWW5*200000,WW_tmp5*1000,my_format5*10
    character WWW6*200000,WW_tmp6*1000,my_format6*10
    character WWW7*200000,WW_tmp7*1000,my_format7*10
    character WWW8*200000,WW_tmp8*1000,my_format8*10
    character WWW9*200000,WW_tmp9*1000,my_format9*10
    character WWW10*200000,WW_tmp10*1000,my_format10*10
    character WWW11*200000,WW_tmp11*1000,my_format11*10
    character WWW12*200000,WW_tmp12*1000,my_format12*10
    character WWW14*200000,WW_tmp14*1000,my_format14*10
    character WWW15*200000,WW_tmp15*1000,my_format15*10
    character WWW16*200000,WW_tmp16*1000,my_format16*10
    character WWW17*200000,WW_tmp17*1000,my_format17*10
    character WWW18*200000,WW_tmp18*1000,my_format18*10
    character WWW19*200000,WW_tmp19*1000,my_format19*10
    character WWW20*200000,WW_tmp20*1000,my_format20*10
                        
    character*1 tab
    integer rs, imol, ifos, iat
    
    !character(60) out_limits
    

    tab = char(9)
    
    WWW=''
    WWW1=''
    WWW2=''
    WWW3=''
    WWW4=''
    WWW5=''
    WWW6=''
    WWW7=''
    WWW8=''
    WWW9=''
    WWW10=''
    WWW11=''
    WWW12=''
    WWW14=''
    WWW15=''
    WWW16=''
    WWW17=''
    WWW18=''
    WWW19=''
    WWW20=''
    
    WW_tmp=''
    write(WW_tmp,'(I0,A1)') istep,char(10)
    WWW1=trim(WWW1)//trim(WW_tmp)
    WWW2=trim(WWW2)//trim(WW_tmp)
    WWW3=trim(WWW3)//trim(WW_tmp)
    WWW4=trim(WWW4)//trim(WW_tmp)
    WWW5=trim(WWW5)//trim(WW_tmp)
    WWW6=trim(WWW6)//trim(WW_tmp)
    WWW7=trim(WWW7)//trim(WW_tmp)
    WWW8=trim(WWW8)//trim(WW_tmp)
    WWW9=trim(WWW9)//trim(WW_tmp)
    WWW10=trim(WWW10)//trim(WW_tmp)
    WWW11=trim(WWW11)//trim(WW_tmp)
    WWW12=trim(WWW12)//trim(WW_tmp)
    WWW14=trim(WWW14)//trim(WW_tmp)
    WWW15=trim(WWW15)//trim(WW_tmp)
    WWW16=trim(WWW16)//trim(WW_tmp)
    WWW17=trim(WWW17)//trim(WW_tmp)
    WWW18=trim(WWW18)//trim(WW_tmp)
    WWW19=trim(WWW19)//trim(WW_tmp)
    WWW20=trim(WWW20)//trim(WW_tmp)
                    
    WW_tmp=''
    WW_tmp1=''
    WW_tmp2=''
    WW_tmp3=''
    WW_tmp4=''
    WW_tmp5=''
    WW_tmp6=''
    WW_tmp7=''
    WW_tmp8=''
    WW_tmp9=''
    WW_tmp10=''
    WW_tmp11=''
    WW_tmp12=''
    WW_tmp14=''
    WW_tmp15=''
    WW_tmp16=''
    WW_tmp17=''
    WW_tmp18=''
    WW_tmp19=''
    WW_tmp20=''
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
            WWW1=trim(WWW1)//trim(WW_tmp)
            WWW2=trim(WWW2)//trim(WW_tmp)
            WWW3=trim(WWW3)//trim(WW_tmp)
            WWW4=trim(WWW4)//trim(WW_tmp)
            WWW5=trim(WWW5)//trim(WW_tmp)
            WWW6=trim(WWW6)//trim(WW_tmp)
            WWW7=trim(WWW7)//trim(WW_tmp)
            WWW8=trim(WWW8)//trim(WW_tmp)
            WWW9=trim(WWW9)//trim(WW_tmp)
            WWW10=trim(WWW10)//trim(WW_tmp)
            WWW11=trim(WWW11)//trim(WW_tmp)
            WWW12=trim(WWW12)//trim(WW_tmp)
            WWW14=trim(WWW14)//trim(WW_tmp)
            WWW15=trim(WWW15)//trim(WW_tmp)
            WWW16=trim(WWW16)//trim(WW_tmp)
            WWW17=trim(WWW17)//trim(WW_tmp)
            WWW18=trim(WWW18)//trim(WW_tmp)
            WWW19=trim(WWW19)//trim(WW_tmp)
            WWW20=trim(WWW20)//trim(WW_tmp)

            WW_tmp3=''
            if (allocated(at(rs)%fosgrp).eqv..true.) then! If this residue has an allocated free energy of solvation 
                do ifos=1,at(rs)%nfosgrps
                    WW_tmp3=''

                    write(WW_tmp3,*) at(rs)%fosgrp(ifos)%wts
                    WWW3=trim(WWW3)//trim(WW_tmp3)
                    WW_tmp2=''
                    write(WW_tmp2,*) at(rs)%fosgrp(ifos)%val(1) ! Martin : The 1 is to get only the actual free energy
                    WWW2=trim(WWW2)//trim(WW_tmp2)
                    
                    WW_tmp20=''
                    write(WW_tmp20,*) at(rs)%fosgrp(ifos)%ats
                    WWW20=trim(WWW20)//trim(WW_tmp20)
                    WWW20=trim(WWW20)//tab
                    
                end do
            end if
            WWW3=trim(WWW3)//char(10)
            WWW2=trim(WWW2)//char(10)
            WWW20=trim(WWW20)//char(10)
            if (allocated(at(rs)%dpgrp).eqv..true.) then
                do iqg=1,at(rs)%ndpgrps
                    do icat=1, at(rs)%dpgrp(iqg)%nats
                        
                        WW_tmp1=''
                        write(WW_tmp1,*) atq(at(rs)%dpgrp(iqg)%ats(icat))
                        WWW1=trim(WWW1)//trim(WW_tmp1)
                        
!                        WW_tmp4=''
!                        write(WW_tmp4,*) blen(at(rs)%dpgrp(iqg)%ats(icat))
!                        WWW4=trim(WWW4)//trim(WW_tmp4)
                        
                        WW_tmp5=''
                        write(WW_tmp5,*) atvol(at(rs)%dpgrp(iqg)%ats(icat))
                        WWW5=trim(WWW5)//trim(WW_tmp5)
                        
                        WW_tmp6=''
                        write(WW_tmp6,*) at(rs)%dpgrp(iqg)%ats(icat)
                        WWW6=trim(WWW6)//trim(WW_tmp6)
                        
                        WW_tmp7=''
                        write(WW_tmp7,*) atsav(at(rs)%dpgrp(iqg)%ats(icat))
                        WWW7=trim(WWW7)//trim(WW_tmp7)
                        
                        WW_tmp8=''
                        write(WW_tmp8,*) lj_rad(attyp(at(rs)%dpgrp(iqg)%ats(icat)))
                        WWW8=trim(WWW8)//trim(WW_tmp8)
                        
                        WW_tmp9=''
                        write(WW_tmp9,*) lj_eps(attyp(at(rs)%dpgrp(iqg)%ats(icat)),attyp(at(rs)%dpgrp(iqg)%ats(icat)))
                        WWW9=trim(WWW9)//trim(WW_tmp9)
                        
!                        WW_tmp10=''
!                        write(WW_tmp10,*) at,attyp(at(rs)%dpgrp(iqg)%ats(icat)  
!                        WWW10=trim(WWW10)//trim(WW_tmp10)
                        

                        WW_tmp10=''
                        write(WW_tmp10,*) lj_sig(attyp(at(rs)%dpgrp(iqg)%ats(icat)),attyp(at(rs)%dpgrp(iqg)%ats(icat)))           
                        WWW10=trim(WWW10)//trim(WW_tmp10)
                        
                        WW_tmp11=''
                        write(WW_tmp11,*) atr(at(rs)%dpgrp(iqg)%ats(icat))            
                        WWW11=trim(WWW11)//trim(WW_tmp11)
                        
                        WW_tmp12=''
                        write(WW_tmp12,*) b_type(at(rs)%dpgrp(iqg)%ats(icat))           
                        WWW12=trim(WWW12)//trim(WW_tmp12)
                        
!                        WW_tmp14=''
!                        write(WW_tmp14,*) bang(at(rs)%dpgrp(iqg)%ats(icat))
!                        WWW14=trim(WWW14)//trim(WW_tmp14)

                        WW_tmp15=''
                        write(WW_tmp15,*) atbvol(at(rs)%dpgrp(iqg)%ats(icat))
                        WWW15=trim(WWW15)//trim(WW_tmp15)

!                        WW_tmp16=''
!                        write(WW_tmp16,*) iaa(rs)%par_bl(:,:)
!                        WWW16=trim(WWW16)//trim(WW_tmp16)
!                        
!                        WW_tmp17=''
!                        write(WW_tmp17,*) atbvol(at(rs)%dpgrp(iqg)%ats(icat))
!                        WWW15=trim(WWW17)//trim(WW_tmp17)
!                        
!                        WW_tmp18=''
!                        write(WW_tmp18,*) atbvol(at(rs)%dpgrp(iqg)%ats(icat))
!                        WWW18=trim(WWW18)//trim(WW_tmp18)
!                        
!                        WW_tmp19=''
!                        write(WW_tmp19,*) atbvol(at(rs)%dpgrp(iqg)%ats(icat))
!                        WWW19=trim(WWW19)//trim(WW_tmp19)
                    end do 
                    WWW1=trim(WWW1)//tab
                    WWW4=trim(WWW4)//tab
                    WWW5=trim(WWW5)//tab
                    WWW6=trim(WWW6)//tab
                    WWW7=trim(WWW7)//tab
                    WWW8=trim(WWW8)//tab
                    WWW9=trim(WWW9)//tab
                    WWW10=trim(WWW10)//tab
                    WWW11=trim(WWW11)//tab
                    WWW12=trim(WWW12)//tab
                    WWW14=trim(WWW14)//tab
                    WWW15=trim(WWW15)//tab
                    WWW16=trim(WWW16)//tab
                    WWW17=trim(WWW17)//tab
                    WWW18=trim(WWW18)//tab
                    WWW19=trim(WWW19)//tab
                end do
            end if 
            do i=1,nrsbl(rs)
                WW_tmp4=''
                write(WW_tmp4,*) i,tab
                WWW4=trim(WWW4)//trim(WW_tmp4)
                WWW4=trim(WWW4)//char(10)
                WW_tmp4=''
                write(WW_tmp4,*) iaa(rs)%par_bl(i,1),iaa(rs)%par_bl(i,2),iaa(rs)%par_bl(i,3),iaa(rs)%par_bl(i,4),char(10)
                WWW4=trim(WWW4)//trim(WW_tmp4)
                WWW4=trim(WWW4)//char(10)
            end do
            do i=1,nrsba(rs)
                WW_tmp14=''
                write(WW_tmp14,*) i,tab
                WWW14=trim(WWW14)//trim(WW_tmp14)
                WW_tmp14=''
                write(WW_tmp14,*) iaa(rs)%par_ba(i,1),iaa(rs)%par_ba(i,2),iaa(rs)%par_ba(i,3),iaa(rs)%par_ba(i,4),char(10)
                WWW14=trim(WWW14)//trim(WW_tmp14)
            end do 

            WWW1=trim(WWW1)//char(10)
            WWW4=trim(WWW4)//char(10)
            WWW5=trim(WWW5)//char(10)
            WWW6=trim(WWW6)//char(10)
            WWW7=trim(WWW7)//char(10)
            WWW8=trim(WWW8)//char(10)
            WWW9=trim(WWW9)//char(10)
            WWW10=trim(WWW10)//char(10)
            WWW11=trim(WWW11)//char(10)
            WWW12=trim(WWW12)//char(10)
            WWW14=trim(WWW14)//char(10)
            WWW15=trim(WWW15)//char(10)
            WWW15=trim(WWW15)//char(10)
            WWW16=trim(WWW16)//char(10)
            WWW17=trim(WWW17)//char(10)
            WWW18=trim(WWW18)//char(10)
            WWW19=trim(WWW19)//char(10)
        end do
    end do

    write(my_format1,'(A2,I0,A1)')'(A',len(trim(WWW1)),')'
    write(my_format2,'(A2,I0,A1)')'(A',len(trim(WWW2)),')'
    write(my_format3,'(A2,I0,A1)')'(A',len(trim(WWW3)),')'
    write(my_format4,'(A2,I0,A1)')'(A',len(trim(WWW4)),')'
    write(my_format5,'(A2,I0,A1)')'(A',len(trim(WWW5)),')'
    write(my_format6,'(A2,I0,A1)')'(A',len(trim(WWW6)),')'
    write(my_format7,'(A2,I0,A1)')'(A',len(trim(WWW7)),')'
    write(my_format8,'(A2,I0,A1)')'(A',len(trim(WWW8)),')'
    write(my_format9,'(A2,I0,A1)')'(A',len(trim(WWW9)),')'
    write(my_format10,'(A2,I0,A1)')'(A',len(trim(WWW10)),')'
    write(my_format11,'(A2,I0,A1)')'(A',len(trim(WWW11)),')'
    write(my_format12,'(A2,I0,A1)')'(A',len(trim(WWW12)),')'

    write(my_format14,'(A2,I0,A1)')'(A',len(trim(WWW14)),')'
    write(my_format15,'(A2,I0,A1)')'(A',len(trim(WWW15)),')'
    write(my_format16,'(A2,I0,A1)')'(A',len(trim(WWW16)),')'
    write(my_format17,'(A2,I0,A1)')'(A',len(trim(WWW17)),')'
    write(my_format18,'(A2,I0,A1)')'(A',len(trim(WWW18)),')'
    write(my_format19,'(A2,I0,A1)')'(A',len(trim(WWW19)),')'
    write(my_format20,'(A2,I0,A1)')'(A',len(trim(WWW20)),')'
    
    write(test_file1,my_format1) trim(WWW1)
    write(test_file2,my_format2) trim(WWW2)
    write(test_file3,my_format3) trim(WWW3)
    write(test_file4,my_format4) trim(WWW4)
    write(test_file5,my_format5) trim(WWW5)
    write(test_file6,my_format6) trim(WWW6)
    write(test_file7,my_format7) trim(WWW7)
    write(test_file8,my_format8) trim(WWW8)
    write(test_file9,my_format9) trim(WWW9)
    write(test_file10,my_format10) trim(WWW10)
    write(test_file11,my_format11) trim(WWW11)
    write(test_file12,my_format12) trim(WWW12)
    write(test_file15,my_format14) trim(WWW14)
    write(test_file14,my_format15) trim(WWW15)
!    write(test_file17,my_format16) trim(WWW16)
!    write(test_file18,my_format17) trim(WWW17)
!    write(test_file19,my_format18) trim(WWW18)
    write(test_file17,my_format20) trim(WWW20)
    
    call print_limits_summary

  end subroutine write_pka_details
  
  subroutine print_limits_summary()
     
    use mcsums
    use molecule
    use polypep
    use atoms
    use zmatrix
    use params
    use inter
    use sequen
    use energies
    use cutoffs

    implicit none 
    integer i,j, iqg, icat, k
    character WWW*2000000, WW_tmp*1000 ,my_format*10
    !character WWW13*200000,WW_tmp13*1000,my_format13*10
    character*1 tab
    integer rs, imol, ifos, iat 
!    character(60) out_limits

    
    WWW=''
    tab = char(9)
    WW_tmp=''
    write(WW_tmp,*) "Please note that for the LJ limits, some may be 0 in both. This is normal, and simply means no biotype in the &
    &parameter file use this LJ type",char(10)
    WWW=trim(WWW)//trim(WW_tmp)

    
    
    !WWW=''
    !tab = char(9)
    WW_tmp=''
    write(WW_tmp,*) "eps limits",char(10)
    WWW=trim(WWW)//trim(WW_tmp)

    WW_tmp=''
    do i=1, n_ljtyp
        write(WW_tmp,*) i,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WW_tmp=''
        write(WW_tmp,*) lj_eps_limits(i,i,1),tab,lj_eps_limits(i,i,2),char(10)
        WWW=trim(WWW)//trim(WW_tmp)
    end do 
    
    if (use_POLAR.eqv..true.) then 
        tab = char(9)
        WW_tmp=''
        write(WW_tmp,*) "cglst%tc_limits",char(10)
        WWW=trim(WWW)//trim(WW_tmp)
        do i=1,cglst%ncs
            write(WW_tmp,*) cglst%it(i),cglst%nc(i),tab
            WWW=trim(WWW)//trim(WW_tmp)
            WW_tmp=''
            write(WW_tmp,*) cglst%tc_limits(i,1),tab,cglst%tc_limits(i,2),tab,cglst%tc_limits(i,3),char(10)
            WWW=trim(WWW)//trim(WW_tmp)
        end do

        tab = char(9)    
        WW_tmp=''
        write(WW_tmp,*) "cglst%nc_limits",char(10)
        WWW=trim(WWW)//trim(WW_tmp)
        do i=1,cglst%ncs
            write(WW_tmp,*) cglst%it(i),tab
            WWW=trim(WWW)//trim(WW_tmp)
            WW_tmp=''
            write(WW_tmp,*) cglst%tc_limits(i,1),tab,cglst%tc_limits(i,2),tab,cglst%tc_limits(i,3),char(10)
            WWW=trim(WWW)//trim(WW_tmp)
        end do   
    end if 
!    tab = char(9)
!    WW_tmp=''
!    write(WW_tmp,*) "atom_groups",char(10)
!    WWW=trim(WWW)//trim(WW_tmp)
!      do rs=1,nseq
!        do i=1,at(rs)%ndpgrps
!          do k=1,at(rs)%dpgrp(i)%nats
!            
!            write(WW_tmp,*) i,tab
!            WWW=trim(WWW)//trim(WW_tmp)
!            WW_tmp=''
!            write(WW_tmp,*) "res",rs,"dp",i,"at",at(rs)%dpgrp(i)%ats(k),char(10)
!            WWW=trim(WWW)//trim(WW_tmp)
!        end do    
!    end do 
!      end do 
!    
    
      WW_tmp=''  
      write(WW_tmp,*)"LJ limits per biotype",char(10)
      WWW=trim(WWW)//trim(WW_tmp)
      do i=1,n_biotyp
          WW_tmp='' 
          write(WW_tmp,*) i,bio_ljtyp(i),bio_ljtyp(transform_table(i)),bio_ljtyp(transform_table(his_eqv_table(i))),char(10)
          WWW=trim(WWW)//trim(WW_tmp)
      end do 
    
    tab = char(9)
    WW_tmp=''
    write(WW_tmp,*) "at(rs)%dpgrp(i)%qnm_limits",char(10)
    WWW=trim(WWW)//trim(WW_tmp)

    do rs=1,nseq
        do i=1,at(rs)%ndpgrps
            write(WW_tmp,*) rs,i,tab
            WWW=trim(WWW)//trim(WW_tmp)
            WW_tmp=''
            write(WW_tmp,*) at(rs)%dpgrp(i)%qnm_limits(1),tab,at(rs)%dpgrp(i)%qnm_limits(2)&
            &,tab,at(rs)%dpgrp(i)%qnm_limits(3),char(10)
            WWW=trim(WWW)//trim(WW_tmp)
        end do 
    end do
    

!    tab = char(9)
!    WW_tmp=''
!    write(WW_tmp,*) "at(rs)%dpgrp(i)%tc_limits",char(10)
!    WWW=trim(WWW)//trim(WW_tmp)
!
!    do rs=1,nseq
!        do i=1,at(rs)%ndpgrps
!            write(WW_tmp,*) rs,i,tab
!            WWW=trim(WWW)//trim(WW_tmp)
!            WW_tmp=''
!            write(WW_tmp,*) at(rs)%dpgrp(i)%tc_limits(1),tab,at(rs)%dpgrp(i)%tc_limits(2)&
!            &,tab,at(rs)%dpgrp(i)%tc_limits(3),char(10)
!            WWW=trim(WWW)//trim(WW_tmp)
!            
!
!        end do 
!    end do
    
    
    

    WW_tmp=''
    write(WW_tmp,*) , char(10)
    WWW=trim(WWW)//trim(WW_tmp)
    
    WW_tmp=''
    write(WW_tmp,*) "sig limits",char(10)
    WWW=trim(WWW)//trim(WW_tmp)
    do i=1, n_ljtyp
        write(WW_tmp,*) i,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WW_tmp=''
        write(WW_tmp,*) lj_sig_limits(i,i,1),tab,lj_sig_limits(i,i,2),char(10)
        WWW=trim(WWW)//trim(WW_tmp)
    end do 
    WW_tmp=''
    write(WW_tmp,*), char(10)
    WWW=trim(WWW)//trim(WW_tmp)

    WW_tmp=''
    write(WW_tmp,*) "atr limits",char(10)
    WWW=trim(WWW)//trim(WW_tmp)
    do i=1, n
        WW_tmp=''
        write(WW_tmp,*) i,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WW_tmp=''
        write(WW_tmp,*) atr_limits(i,1),tab,atr_limits(i,2),tab,atr_limits(i,3),char(10)
        WWW=trim(WWW)//trim(WW_tmp)
    end do 
    
    WW_tmp=''
    write(WW_tmp,*) "atvol limits",char(10)
    WWW=trim(WWW)//trim(WW_tmp)
    do i=1, n
        WW_tmp=''
        write(WW_tmp,*) i,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WW_tmp=''
        write(WW_tmp,*) atvol_limits(i,1),tab,atvol_limits(i,2),char(10)
        WWW=trim(WWW)//trim(WW_tmp)
    end do  
    
    
    WW_tmp=''
    write(WW_tmp,*) "atbvol limits",char(10)
    WWW=trim(WWW)//trim(WW_tmp)
    do i=1, n
        WW_tmp=''
        write(WW_tmp,*) i,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WW_tmp=''
        write(WW_tmp,*) atbvol_limits(i,1),tab,atbvol_limits(i,2),char(10)
        WWW=trim(WWW)//trim(WW_tmp)
    end do     
    
        WW_tmp=''
    write(WW_tmp,*) "atstv limits",char(10)
    WWW=trim(WWW)//trim(WW_tmp)
    do i=1, n
        WW_tmp=''
        write(WW_tmp,*) i,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WW_tmp=''
        write(WW_tmp,*) atstv_limits(i,1),tab,atstv_limits(i,2),char(10)
        WWW=trim(WWW)//trim(WW_tmp)
    end do     
    
    
    if (use_POLAR.EQV..true.) then
    
    WW_tmp=''
    write(WW_tmp,*) , char(10)
    WWW=trim(WWW)//trim(WW_tmp)
    WW_tmp=''
    write(WW_tmp,*) "atq limits",char(10)
    WWW=trim(WWW)//trim(WW_tmp)
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            WW_tmp=''
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,char(10)
            WWW=trim(WWW)//trim(WW_tmp)
            
            if (allocated(at(rs)%dpgrp).eqv..true.) then
                do iqg=1,at(rs)%ndpgrps
                    do icat=1, at(rs)%dpgrp(iqg)%nats
                        WW_tmp=''
                        write(WW_tmp,*) at(rs)%dpgrp(iqg)%ats(icat),tab
                        WWW=trim(WWW)//trim(WW_tmp)
                        
                        WW_tmp=''
                        write(WW_tmp,*) atq_limits(at(rs)%dpgrp(iqg)%ats(icat),1),tab,atq_limits(at(rs)%dpgrp(iqg)%ats(icat),3)&
                        &,tab,atq_limits(at(rs)%dpgrp(iqg)%ats(icat),4),char(10)
                        WWW=trim(WWW)//trim(WW_tmp)
                    end do
                end do 
            end if 
        WWW=trim(WWW)//char(10)    
        end do 
    end do 
    

    
    WW_tmp=''
    write(WW_tmp,'(A4,A1)') 'elin',char(10)
    WWW=trim(WWW)//trim(WW_tmp)

    do imol=1,nmol ! Martin : this bit is taken from unbound's first subroutine
        do rs=rsmol(imol,1),rsmol(imol,2)
            WW_tmp=''
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,char(10)
            WWW=trim(WWW)//trim(WW_tmp)
            
            do i=1,nrpolintra(rs)       

                WW_tmp=''
                write(WW_tmp,*) 'atoms1',iaa_limits(rs,1)%polin(i,:),char(10)
                WWW=trim(WWW)//trim(WW_tmp)

                WW_tmp=''
                write(WW_tmp,*) 'atoms2',iaa_limits(rs,2)%polin(i,:),char(10)
                WWW=trim(WWW)//trim(WW_tmp)
                
                WW_tmp=''
                write(WW_tmp,*) 'atoms3',iaa(rs)%polin(i,:),char(10)
                WWW=trim(WWW)//trim(WW_tmp)


                WW_tmp=''
                write(WW_tmp,*) ,fudge_limits(rs,1)%elin(i),tab,fudge_limits(rs,2)%elin(i),fudge(rs)%elin(i),char(10)
                WWW=trim(WWW)//trim(WW_tmp)
            end do 
        end do
    end do
    
    WW_tmp=''
    write(WW_tmp,'(A7,A1)') 'expolin',char(10)
    WWW=trim(WWW)//trim(WW_tmp)

    do imol=1,nmol ! Martin : this bit is taken from unbound's first subroutine
        do rs=rsmol(imol,1),rsmol(imol,2)
            WW_tmp=''
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,char(10)
            WWW=trim(WWW)//trim(WW_tmp)
            
            do i=1,nrpolintra(rs)   
                WW_tmp=''
                write(WW_tmp,*) 'atoms1',iaa_limits(rs,1)%polin(i,:),char(10)
                WWW=trim(WWW)//trim(WW_tmp)

                WW_tmp=''
                write(WW_tmp,*) 'atoms2',iaa_limits(rs,2)%polin(i,:),char(10)
                WWW=trim(WWW)//trim(WW_tmp)

                WW_tmp=''
                write(WW_tmp,*) ,iaa_limits(rs,1)%expolin(i,:),tab,iaa_limits(rs,2)%expolin(i,:),char(10)
                WWW=trim(WWW)//trim(WW_tmp)
            end do 
        end do
    end do
    
    WW_tmp=''
    write(WW_tmp,'(A4,A1)') 'elnb',char(10)
    WWW=trim(WWW)//trim(WW_tmp)

    do imol=1,nmol ! Martin : this bit is taken from unbound's first subroutine
        do rs=rsmol(imol,1),rsmol(imol,2)
            WW_tmp=''
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,char(10)
            WWW=trim(WWW)//trim(WW_tmp)
            
            do i=1,nrpolnb(rs)       

                WW_tmp=''
                write(WW_tmp,*) 'atoms1',iaa_limits(rs,1)%polnb(i,:),char(10)
                WWW=trim(WWW)//trim(WW_tmp)

                WW_tmp=''
                write(WW_tmp,*) 'atoms2',iaa_limits(rs,2)%polnb(i,:),char(10)
                WWW=trim(WWW)//trim(WW_tmp)
                
                WW_tmp=''
                write(WW_tmp,*) 'atoms3',iaa(rs)%polnb(i,:),char(10)
                WWW=trim(WWW)//trim(WW_tmp)


                WW_tmp=''
                write(WW_tmp,*) ,fudge_limits(rs,1)%elnb(i),tab,fudge_limits(rs,2)%elnb(i),fudge(rs)%elnb(i),char(10)
                WWW=trim(WWW)//trim(WW_tmp)
            end do 
        end do
    end do
    
    WW_tmp=''
    write(WW_tmp,'(A7,A1)') 'expolnb',char(10)
    WWW=trim(WWW)//trim(WW_tmp)

    do imol=1,nmol ! Martin : this bit is taken from unbound's first subroutine
        do rs=rsmol(imol,1),rsmol(imol,2)
            WW_tmp=''
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,char(10)
            WWW=trim(WWW)//trim(WW_tmp)
            
            do i=1,nrpolnb(rs)   
                WW_tmp=''
                write(WW_tmp,*) 'atoms1',iaa_limits(rs,1)%polnb(i,:),char(10)
                WWW=trim(WWW)//trim(WW_tmp)

                WW_tmp=''
                write(WW_tmp,*) 'atoms2',iaa_limits(rs,2)%polnb(i,:),char(10)
                WWW=trim(WWW)//trim(WW_tmp)

                WW_tmp=''
                write(WW_tmp,*) ,iaa_limits(rs,1)%expolnb(i,:),tab,iaa_limits(rs,2)%expolnb(i,:),iaa(rs)%expolnb(i,:),char(10)
                WWW=trim(WWW)//trim(WW_tmp)
                
            end do 
        end do
    end do
    end if 
    
    WW_tmp=''
    write(WW_tmp,*) "eps limits",char(10)
    WWW=trim(WWW)//trim(WW_tmp)
    do i=1, n_ljtyp
        do j=1,n_ljtyp
            WW_tmp=''
            write(WW_tmp,*) "type",i,tab,j
            WWW=trim(WWW)//trim(WW_tmp)

            WW_tmp=''
            write(WW_tmp,*) lj_eps_limits(i,j,1),lj_eps_limits(i,j,2),char(10)
            WWW=trim(WWW)//trim(WW_tmp)
        end do 
    end do 
    
        WW_tmp=''
    write(WW_tmp,*) "sig limits",char(10)
    WWW=trim(WWW)//trim(WW_tmp)
    do i=1, n_ljtyp
        do j=1,n_ljtyp
            WW_tmp=''
            write(WW_tmp,*) "type",i,tab,j
            WWW=trim(WWW)//trim(WW_tmp)

            WW_tmp=''
            write(WW_tmp,*) lj_sig_limits(i,j,1),lj_sig_limits(i,j,2),char(10)
            WWW=trim(WWW)//trim(WW_tmp)
        end do 
    end do 
    
    
    
    
    WW_tmp=''
    write(WW_tmp,*) "Fudges Limits intra res rsin"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)

    do rs=1,nseq
        WW_tmp=''
        write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
        
        do i=1,nrsintra(rs)
            WW_tmp=''
            write(WW_tmp,*) 'atoms1',iaa(rs)%atin(i,1),iaa(rs)%atin(i,2)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

            write(WW_tmp,*) 'atoms2',iaa_limits(rs,1)%atin(i,1),iaa_limits(rs,1)%atin(i,2)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

            write(WW_tmp,*) 'atoms3',iaa_limits(rs,2)%atin(i,1),iaa_limits(rs,2)%atin(i,2)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

            WW_tmp=''
            write(WW_tmp,*) fudge_limits(rs,1)%rsin(i),tab,fudge_limits(rs,2)%rsin(i)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

        end do     
    end do 
   WW_tmp=''
    write(WW_tmp,*) "Fudges Limits intra res rsnb"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)

    do rs=1,nseq
        WW_tmp=''
        write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
        
        do i=1,nrsnb(rs)
            WW_tmp=''
            write(WW_tmp,*) 'atoms',iaa(rs)%atnb(i,1),iaa(rs)%atnb(i,2)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
            
            write(WW_tmp,*) 'atoms1',iaa_limits(rs,1)%atnb(i,1),iaa_limits(rs,1)%atnb(i,2)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

            write(WW_tmp,*) 'atoms2',iaa_limits(rs,2)%atnb(i,1),iaa_limits(rs,2)%atnb(i,2)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)


            WW_tmp=''
            write(WW_tmp,*) fudge_limits(rs,1)%rsnb(i),tab,fudge_limits(rs,2)%rsnb(i)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

        end do     
    end do 

   WW_tmp=''
    write(WW_tmp,*) "Fudges Limits intra res rsin_ljs"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)

    do rs=1,nseq
        WW_tmp=''
        write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
        
        do i=1,nrsintra(rs)
            WW_tmp=''
            write(WW_tmp,*) 'atoms',iaa(rs)%atin(i,1),iaa(rs)%atin(i,2)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
            
!            WW_tmp=''
!            write(WW_tmp,*) 'LJ_type1',bio_ljtyp(b_type(iaa(rs)%atin(i,1))),bio_ljtyp(b_type(iaa(rs)%atin(i,2)))
!            WWW=trim(WWW)//trim(WW_tmp)
!            WWW=trim(WWW)//char(10)
!                        
!            WW_tmp=''
!            write(WW_tmp,*) 'LJ_type2',bio_ljtyp(transform_table(b_type(iaa(rs)%atin(i,1)))),&
!            &transform_table(bio_ljtyp(b_type(iaa(rs)%atin(i,2))))
!            WWW=trim(WWW)//trim(WW_tmp)
!            WWW=trim(WWW)//char(10)

            WW_tmp=''
            write(WW_tmp,*) fudge_limits(rs,1)%rsin_ljs(i),tab,fudge_limits(rs,2)%rsin_ljs(i)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

        end do     
    end do 
      WW_tmp=''
    write(WW_tmp,*) "Fudges Limits intra res rsnb_ljs"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)

    do rs=1,nseq
        WW_tmp=''
        write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
        
        do i=1,nrsnb(rs)
            WW_tmp=''
            write(WW_tmp,*) 'atoms',iaa(rs)%atnb(i,1),iaa(rs)%atnb(i,2),char(10)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
            
!            WW_tmp=''
!            write(WW_tmp,*) 'LJ_type1',bio_ljtyp(b_type(iaa(rs)%atnb(i,1))),bio_ljtyp(b_type(iaa(rs)%atnb(i,2)))
!            WWW=trim(WWW)//trim(WW_tmp)
!            WWW=trim(WWW)//char(10)
!                        
!            WW_tmp=''
!            write(WW_tmp,*) 'LJ_type2',bio_ljtyp(transform_table(b_type(iaa(rs)%atnb(i,1))))&
!            &,transform_table(bio_ljtyp(b_type(iaa(rs)%atnb(i,2))))
!            WWW=trim(WWW)//trim(WW_tmp)
!            WWW=trim(WWW)//char(10)

            WW_tmp=''
            write(WW_tmp,*) fudge_limits(rs,1)%rsnb_ljs(i),tab,fudge_limits(rs,2)%rsnb_ljs(i)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

        end do     
    end do 
    
    write(WW_tmp,*) "Fudges Limits intra res rsin_lje"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)

    do rs=1,nseq
        WW_tmp=''
        write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
        
        do i=1,nrsintra(rs)
            WW_tmp=''
            write(WW_tmp,*) 'atoms',iaa(rs)%atin(i,1),iaa(rs)%atin(i,2)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
                        
!            WW_tmp=''
!            write(WW_tmp,*) 'LJ_type1',bio_ljtyp(b_type(iaa(rs)%atin(i,1))),bio_ljtyp(b_type(iaa(rs)%atin(i,2)))
!            WWW=trim(WWW)//trim(WW_tmp)
!            WWW=trim(WWW)//char(10)
!                        
!            WW_tmp=''
!            write(WW_tmp,*) 'LJ_type2',bio_ljtyp(transform_table(b_type(iaa(rs)%atin(i,1)))),&
!            &transform_table(bio_ljtyp(b_type(iaa(rs)%atin(i,2))))
!            WWW=trim(WWW)//trim(WW_tmp)
!            WWW=trim(WWW)//char(10)


            WW_tmp=''
            write(WW_tmp,*) fudge_limits(rs,1)%rsin_lje(i),tab,fudge_limits(rs,2)%rsin_lje(i)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

        end do     
    end do 
    write(WW_tmp,*) "Fudges Limits intra res rsnb_lje"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)

    do rs=1,nseq
        WW_tmp=''
        write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
        
        do i=1,nrsnb(rs)
            WW_tmp=''
            write(WW_tmp,*) 'atoms',iaa(rs)%atnb(i,1),iaa(rs)%atnb(i,2),char(10)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
                        
!            WW_tmp=''
!            write(WW_tmp,*) 'LJ_type1',bio_ljtyp(b_type(iaa(rs)%atnb(i,1))),bio_ljtyp(b_type(iaa(rs)%atnb(i,2)))
!            WWW=trim(WWW)//trim(WW_tmp)
!            WWW=trim(WWW)//char(10)
!                        
!            WW_tmp=''
!            write(WW_tmp,*) 'LJ_type2',bio_ljtyp(transform_table(b_type(iaa(rs)%atnb(i,1)))),&
!            &transform_table(bio_ljtyp(b_type(iaa(rs)%atnb(i,2))))
!            WWW=trim(WWW)//trim(WW_tmp)
!            WWW=trim(WWW)//char(10)


            WW_tmp=''
            write(WW_tmp,*) fudge_limits(rs,1)%rsnb_lje(i),tab,fudge_limits(rs,2)%rsnb_lje(i)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

        end do     
    end do 
    
    write(WW_tmp,*) "Bond length limits"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)

    do i=1,n
        WW_tmp=''
        write(WW_tmp,'(A5,I0,A1)') 'atom_',i,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WW_tmp=''
        write(WW_tmp,*) blen_limits(i,1),tab,blen_limits(i,2)
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)

   
    end do 
    
    write(WW_tmp,*) "Bond angle limits"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)

    do i=1,n
        WW_tmp=''
        write(WW_tmp,'(A5,I0,A1)') 'atom_',i,tab
        WWW=trim(WWW)//trim(WW_tmp)

        WW_tmp=''
        write(WW_tmp,*) bang_limits(i,1),tab,bang_limits(i,2)
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
    end do 
    write(WW_tmp,*) "atsavred limits ",char(10)
    WWW=trim(WWW)//trim(WW_tmp)
!    WWW=trim(WWW)//char(10)
    do i=1,n
        WW_tmp=''
        write(WW_tmp,'(A5,I0,A1)') 'atom_',i,tab
        WWW=trim(WWW)//trim(WW_tmp)

        WW_tmp=''
        write(WW_tmp,*) atsavred_limits(i,1),tab,atsavred_limits(i,2)
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
    end do 

    
    write(WW_tmp,*) "atsavmaxfr"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do i=1,n
        WW_tmp=''
        write(WW_tmp,'(A5,I0,A1)') 'atom_',i,tab
        WWW=trim(WWW)//trim(WW_tmp)

        WW_tmp=''
        write(WW_tmp,*) atsavmaxfr_limits(i,1),tab,atsavmaxfr_limits(i,2)
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
    end do 
    
  if (use_IMPSOLV.EQV..true.) then
    write(WW_tmp,*) "fos_value"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
            WWW=trim(WWW)//trim(WW_tmp)
                do ifos=1,at(rs)%nfosgrps
                     WW_tmp=''

                    write(WW_tmp,*) at(rs)%fosgrp(ifos)%val_limits(1,1),tab,at(rs)%fosgrp(ifos)%val_limits(1,2) 
                    ! Martin : The 1 is to get only the acutual free energy
                    WWW=trim(WWW)//trim(WW_tmp)
                end do
            WWW=trim(WWW)//char(10)
        end do 
    end do 
    write(WW_tmp,*) "fos_weights"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
                do ifos=1,at(rs)%nfosgrps
                     WW_tmp=''
                    write(WW_tmp,*) at(rs)%fosgrp(ifos)%wts_limits(1,1),tab,at(rs)%fosgrp(ifos)%wts_limits(1,2) 
                    ! Martin : The 1 is to get only the acutual free energy
                    WWW=trim(WWW)//trim(WW_tmp)
                end do
            WWW=trim(WWW)//char(10)
        end do 
    end do 
  end if 
    write(WW_tmp,*) "par_bl"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,char(10)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
                do i=1,size(iaa(rs)%typ_bl(:))
                    WW_tmp=''
                    write(WW_tmp,*) "atoms",iaa(rs)%bl(i,1),iaa(rs)%bl(i,2)
                    WWW=trim(WWW)//trim(WW_tmp)
                    WWW=trim(WWW)//char(10)
                    do j=1,size(iaa(rs)%par_bl(i,:))
   
                        
                        WW_tmp=''
                        write(WW_tmp,*) iaa_limits(rs,1)%par_bl(i,j),tab,iaa_limits(rs,2)%par_bl(i,j) 
                        WWW=trim(WWW)//trim(WW_tmp)
                        WWW=trim(WWW)//char(10)
                    end do 
                end do 
            WWW=trim(WWW)//char(10)
        end do 
    end do 
    
    write(WW_tmp,*) "par_ba"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,char(10)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
                do i=1,nrsba(rs)
                    WW_tmp=''
                    write(WW_tmp,*) "atoms",iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),iaa(rs)%ba(i,3) 
                    WWW=trim(WWW)//trim(WW_tmp)
                    WWW=trim(WWW)//char(10)
                    do j=1,size(iaa(rs)%par_ba(i,:))
                        WW_tmp=''
                        write(WW_tmp,*) iaa_limits(rs,1)%par_ba(i,j),tab,iaa_limits(rs,2)%par_ba(i,j) 
                        WWW=trim(WWW)//trim(WW_tmp)
                        WWW=trim(WWW)//char(10)
                    end do 
                end do 
            WWW=trim(WWW)//char(10)
        end do 
    end do 
    
    write(WW_tmp,*) "par_di"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            
            if (natres(rs).gt.1) then 
                write(WW_tmp,'(A4,I0,A1)') 'res_',rs,char(10)
                WWW=trim(WWW)//trim(WW_tmp)
                WWW=trim(WWW)//char(10)
                    do i=1,nrsdieff(rs)
                        WW_tmp=''
                        write(WW_tmp,*) "atoms",iaa(rs)%di(i,1),iaa(rs)%di(i,2),iaa(rs)%di(i,3),iaa(rs)%di(i,4)
                        WWW=trim(WWW)//trim(WW_tmp)
                        WWW=trim(WWW)//char(10)
                        do j=1,size(iaa(rs)%par_di(i,:))

                            WW_tmp=''
                            write(WW_tmp,*) iaa_limits(rs,1)%par_di(i,j),tab,iaa_limits(rs,2)%par_di(i,j) 
                            WWW=trim(WWW)//trim(WW_tmp)
                            WWW=trim(WWW)//char(10)
                        end do 
                    end do 
                WWW=trim(WWW)//char(10)
            end if 
        end do 
    end do 
    
    write(WW_tmp,*) "par_impt"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,char(10)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
            do i=1,nrsimpteff(rs)
                write(WW_tmp,*) "atoms",iaa(rs)%impt(i,1),iaa(rs)%impt(i,2),iaa(rs)%impt(i,3),iaa(rs)%impt(i,4)
                WWW=trim(WWW)//trim(WW_tmp)
                WWW=trim(WWW)//char(10)
                do j=1,size(iaa(rs)%par_impt(i,:))
                    WW_tmp=''
                    write(WW_tmp,*) iaa_limits(rs,1)%par_impt(i,j),tab,iaa_limits(rs,2)%par_impt(i,j) 
                    WWW=trim(WWW)//trim(WW_tmp)
                    WWW=trim(WWW)//char(10)
                end do 
            end do 
            WWW=trim(WWW)//char(10)
        end do 
    end do 
    
    
    if (use_POLAR.EQV..true.) then
        do rs=1,nseq
            WW_tmp=''
            write(WW_tmp,*) "nrpolnb",nrpolnb(rs)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
            do i=1,nrpolnb(rs)
                WW_tmp=''
                write(WW_tmp,*) "fudgenb",rs,iaa(rs)%polnb(i,1),iaa(rs)%polnb(i,2),&
                &fudge_limits(rs,1)%elnb(i),fudge_limits(rs,2)%elnb(i)
                WWW=trim(WWW)//trim(WW_tmp)
                WWW=trim(WWW)//char(10)

            end do 
        end do 

        do rs=1,nseq
            WW_tmp=''
            write(WW_tmp,*) "nrsintra",nrpolintra(rs)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

            do i=1,nrpolintra(rs)
                WW_tmp=''
                write(WW_tmp,*) "fudgein",rs,iaa(rs)%polin(i,1),iaa(rs)%polin(i,2),&
                &fudge_limits(rs,1)%elin(i),fudge_limits(rs,2)%elin(i)
                WWW=trim(WWW)//trim(WW_tmp)
                WWW=trim(WWW)//char(10)

            end do 
        end do 
    end if 
!    do i=1,n
!        WW_tmp=''
!        write(WW_tmp,*) "scrq",i,scrq(i)
!        WWW=trim(WWW)//trim(WW_tmp)
!        WWW=trim(WWW)//char(10)
!    end do
    do j=1,6
    do i=1,n
        WW_tmp=''
        write(WW_tmp,*) "atsavprm",i,atsavprm_limits(i,j,1),atsavprm_limits(i,j,2)
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
    end do 
    end do 

    write(my_format,'(A2,I0,A1)')'(A',len(trim(WWW)),')'
    write(file_limits,my_format) trim(WWW)
    
  end subroutine print_limits_summary
  
  
  
  subroutine print_values_summary()
     
    use mcsums
    use molecule
    use polypep
    use atoms
    use zmatrix
    use params
    use inter
    use sequen
    use energies 

    implicit none 
    integer i,j, iqg, icat, k
    character WWW*200000, WW_tmp*1000 ,my_format*10
    !character WWW13*200000,WW_tmp13*1000,my_format13*10
    character*1 tab
    integer rs, imol, ifos, iat 
   ! character(60) out_limits

    WWW=''
    tab = char(9)
    WW_tmp=''
    write(WW_tmp,*) "eps",char(10)
    WWW=trim(WWW)//trim(WW_tmp)

    WW_tmp=''
    do i=1, n_ljtyp
        write(WW_tmp,*) i,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WW_tmp=''
        write(WW_tmp,*) lj_eps(i,i),char(10)
        WWW=trim(WWW)//trim(WW_tmp)
    end do 

    WW_tmp=''
    write(WW_tmp,*) , char(10)
    WWW=trim(WWW)//trim(WW_tmp)
    
    WW_tmp=''
    write(WW_tmp,*) "sig",char(10)
    WWW=trim(WWW)//trim(WW_tmp)
    do i=1, n_ljtyp
        write(WW_tmp,*) i,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WW_tmp=''
        write(WW_tmp,*) lj_sig(i,i),char(10)
        WWW=trim(WWW)//trim(WW_tmp)
    end do 
    WW_tmp=''
    write(WW_tmp,*), char(10)
    WWW=trim(WWW)//trim(WW_tmp)

    WW_tmp=''
    write(WW_tmp,*) "atr",char(10)
    WWW=trim(WWW)//trim(WW_tmp)
    do i=1, n
        WW_tmp=''
        write(WW_tmp,*) i,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WW_tmp=''
        write(WW_tmp,*) atr(i),char(10)
        WWW=trim(WWW)//trim(WW_tmp)
    end do 
    

    
    
    WW_tmp=''
    write(WW_tmp,*) , char(10)
    WWW=trim(WWW)//trim(WW_tmp)
    WW_tmp=''
    write(WW_tmp,*) "atq",char(10)
    WWW=trim(WWW)//trim(WW_tmp)
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            WW_tmp=''
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
            WWW=trim(WWW)//trim(WW_tmp)
            if (allocated(at(rs)%dpgrp).eqv..true.) then
                do iqg=1,at(rs)%ndpgrps
                    do icat=1,at(rs)%dpgrp(iqg)%nats
                        WW_tmp=''
                        write(WW_tmp,*) atq(at(rs)%dpgrp(iqg)%ats(icat)),char(10)
                        WWW=trim(WWW)//trim(WW_tmp)
                    end do
                end do 
            end if 
        WWW=trim(WWW)//char(10)    
        end do 
    end do 
!    

    WW_tmp=''
    write(WW_tmp,*) "Fudges Limits intra res rsin"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)

    do rs=1,nseq
        WW_tmp=''
        write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
        
        do i=1,nrsintra(rs)
            WW_tmp=''
            write(WW_tmp,*) 'atoms',iaa(rs)%atin(i,1),iaa(rs)%atin(i,2)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

            WW_tmp=''
            write(WW_tmp,*) fudge(rs)%rsin(i)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

        end do     
    end do 
 WW_tmp=''
    write(WW_tmp,*) "Fudges intra res rsnb"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)

    do rs=1,nseq
        WW_tmp=''
        write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
        
        do i=1,nrsnb(rs)
            WW_tmp=''
            write(WW_tmp,*) 'atoms',iaa(rs)%atnb(i,1),iaa(rs)%atnb(i,2)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

            WW_tmp=''
            write(WW_tmp,*) fudge(rs)%rsnb(i)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

        end do     
    end do 

    WW_tmp=''
    write(WW_tmp,*) "Fudges intra res rsin_ljs"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)

    do rs=1,nseq
        WW_tmp=''
        write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
        
        do i=1,nrsintra(rs)
            WW_tmp=''
            write(WW_tmp,*) 'atoms',iaa(rs)%atin(i,1),iaa(rs)%atin(i,2)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
            
            WW_tmp=''
            write(WW_tmp,*) fudge(rs)%rsin_ljs(i)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

        end do     
    end do 
    
    WW_tmp=''
    write(WW_tmp,*) "Fudges intra res rsnb_ljs"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)

    do rs=1,nseq
        WW_tmp=''
        write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
        
        do i=1,nrsnb(rs)
            WW_tmp=''
            write(WW_tmp,*) 'atoms',iaa(rs)%atnb(i,1),iaa(rs)%atnb(i,2)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
            
            WW_tmp=''
            write(WW_tmp,*) fudge(rs)%rsnb_ljs(i)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

        end do     
    end do 
    
write(WW_tmp,*) "Fudges intra res rsin_lje"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)

    do rs=1,nseq
        
        WW_tmp=''
        write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
        
        do i=1,nrsintra(rs)
            WW_tmp=''
            write(WW_tmp,*) 'atoms',iaa(rs)%atin(i,1),iaa(rs)%atin(i,2)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

            WW_tmp=''
            write(WW_tmp,*) fudge(rs)%rsin_lje(i)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

        end do     
    end do 
    
    
    write(WW_tmp,*) "Fudges intra res rsnb_lje"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)

    do rs=1,nseq
        WW_tmp=''
        write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
        
        do i=1,nrsnb(rs)
            WW_tmp=''
            write(WW_tmp,*) 'atoms',iaa(rs)%atnb(i,1),iaa(rs)%atnb(i,2)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
            
            WW_tmp=''
            write(WW_tmp,*) fudge(rs)%rsnb_lje(i)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)

        end do     
    end do 
    
    write(WW_tmp,*) "Bond length"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)

    do i=1,n
        WW_tmp=''
        write(WW_tmp,'(A5,I0,A1)') 'atom_',i,tab
        WWW=trim(WWW)//trim(WW_tmp)
        WW_tmp=''
        write(WW_tmp,*) blen(i)
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)

   
    end do 
    
    write(WW_tmp,*) "Bond angle"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do i=1,n
        WW_tmp=''
        write(WW_tmp,'(A5,I0,A1)') 'atom_',i,tab
        WWW=trim(WWW)//trim(WW_tmp)

        WW_tmp=''
        write(WW_tmp,*) bang(i)
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
    end do 

    write(WW_tmp,*) "Bond torsion"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do i=1,n
        WW_tmp=''
        write(WW_tmp,'(A5,I0,A1)') 'atom_',i,tab
        WWW=trim(WWW)//trim(WW_tmp)

        WW_tmp=''
        write(WW_tmp,*) ztor(i)
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
    end do 

    write(WW_tmp,*) "atsavred "
    WWW=trim(WWW)//trim(WW_tmp)
!    WWW=trim(WWW)//char(10)
    do i=1,n
        WW_tmp=''
        write(WW_tmp,'(A5,I0,A1)') 'atom_',i,tab
        WWW=trim(WWW)//trim(WW_tmp)

        WW_tmp=''
        write(WW_tmp,*) atsavred(i)
        WWW=trim(WWW)//WW_tmp
        WWW=trim(WWW)//char(10)
    end do 

    
    write(WW_tmp,*) "atsavmaxfr"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do i=1,n
        WW_tmp=''
        write(WW_tmp,'(A5,I0,A1)') 'atom_',i,tab
        WWW=trim(WWW)//trim(WW_tmp)

        WW_tmp=''
        write(WW_tmp,*) atsavmaxfr(i)
        WWW=trim(WWW)//trim(WW_tmp)
        WWW=trim(WWW)//char(10)
    end do 
    
    write(WW_tmp,*) "Fos_value"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
            WWW=trim(WWW)//trim(WW_tmp)
                do ifos=1,at(rs)%nfosgrps
                     WW_tmp=''
                    write(WW_tmp,*) at(rs)%fosgrp(ifos)%val(1) ! Martin : The 1 is to get only the actual free energy
                    WWW=trim(WWW)//trim(WW_tmp)
                end do
            WWW=trim(WWW)//char(10)
        end do 
    end do 
    
    write(WW_tmp,*) "Fos_weights"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,tab
            WWW=trim(WWW)//trim(WW_tmp)
                do ifos=1,at(rs)%nfosgrps
                    write(WW_tmp,'(A4,I0,A1)') 'fosgrp',ifos,tab
                    WWW=trim(WWW)//trim(WW_tmp)
                    do i=1,at(rs)%fosgrp(ifos)%nats
                        WW_tmp=''
                        write(WW_tmp,*) at(rs)%fosgrp(ifos)%ats(i),at(rs)%fosgrp(ifos)%wts(i) 
                        WWW=trim(WWW)//trim(WW_tmp)
                        WWW=trim(WWW)//char(10)
                    end do
                end do
            WWW=trim(WWW)//char(10)
        end do
    end do 
    
    write(WW_tmp,*) "par_bl"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,char(10)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
                do i=1,size(iaa(rs)%typ_bl(:))
                    WW_tmp=''
                    write(WW_tmp,*) "atoms",iaa(rs)%bl(i,1),iaa(rs)%bl(i,2)
                    WWW=trim(WWW)//trim(WW_tmp)
                    WWW=trim(WWW)//char(10)
                    do j=1,size(iaa(rs)%par_bl(i,:))
                        WW_tmp=''
                        write(WW_tmp,*) iaa(rs)%par_bl(i,j)
                        WWW=trim(WWW)//trim(WW_tmp)
                        WWW=trim(WWW)//char(10)
                    end do 
                end do 
            WWW=trim(WWW)//char(10)
        end do 
    end do 
    
    write(WW_tmp,*) "par_ba"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,char(10)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
                do i=1,nrsba(rs)
                    WW_tmp=''
                    write(WW_tmp,*) "atoms",iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),iaa(rs)%ba(i,3) 
                    WWW=trim(WWW)//trim(WW_tmp)
                    WWW=trim(WWW)//char(10)
                    do j=1,size(iaa(rs)%par_ba(i,:))
                        WW_tmp=''
                        write(WW_tmp,*) iaa(rs)%par_ba(i,j)
                        WWW=trim(WWW)//trim(WW_tmp)
                        WWW=trim(WWW)//char(10)
                    end do 
                end do 
            WWW=trim(WWW)//char(10)
        end do 
    end do 
    
    write(WW_tmp,*) "par_di"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            if (natres(rs).gt.1) then 
                write(WW_tmp,'(A4,I0,A1)') 'res_',rs,char(10)
                WWW=trim(WWW)//trim(WW_tmp)
                WWW=trim(WWW)//char(10)
                    do i=1,nrsdieff(rs)
                        WW_tmp=''
                        write(WW_tmp,*) "atoms",iaa(rs)%di(i,1),iaa(rs)%di(i,2),iaa(rs)%di(i,3),iaa(rs)%di(i,4)
                        WWW=trim(WWW)//trim(WW_tmp)
                        WWW=trim(WWW)//char(10)
                        
                        do j=1,size(iaa(rs)%par_di(i,:))
                            WW_tmp=''
                            write(WW_tmp,*) iaa(rs)%par_di(i,j)
                            WWW=trim(WWW)//trim(WW_tmp)
                            WWW=trim(WWW)//char(10)
                        end do 
                        
                    end do 
                WWW=trim(WWW)//char(10)
            end if 
        end do 
    end do 
    
    
    write(WW_tmp,*) "par_impt"
    WWW=trim(WWW)//trim(WW_tmp)
    WWW=trim(WWW)//char(10)
    do imol=1,nmol! Goes through all molecules
        do rs=rsmol(imol,1),rsmol(imol,2)! Goes through all residues
            write(WW_tmp,'(A4,I0,A1)') 'res_',rs,char(10)
            WWW=trim(WWW)//trim(WW_tmp)
            WWW=trim(WWW)//char(10)
            do i=1,nrsimpteff(rs)
                write(WW_tmp,*) "atoms",iaa(rs)%impt(i,1),iaa(rs)%impt(i,2),iaa(rs)%impt(i,3),iaa(rs)%impt(i,4)
                WWW=trim(WWW)//trim(WW_tmp)
                WWW=trim(WWW)//char(10)
                do j=1,size(iaa(rs)%par_impt(i,:))
                    WW_tmp=''
                    write(WW_tmp,*) iaa(rs)%par_impt(i,j)
                    WWW=trim(WWW)//trim(WW_tmp)
                    WWW=trim(WWW)//char(10)
                end do 
            end do 
            WWW=trim(WWW)//char(10)
        end do 
    end do 
    
    if (use_IMPSOLV.EQV..true.) then
        do j=1,6
            do i=1,n
                WW_tmp=''
                write(WW_tmp,*) "atsavprm",i,atsavprm(i,j)
                WWW=trim(WWW)//trim(WW_tmp)
                WWW=trim(WWW)//char(10)
            end do 
        end do 
    end if
    
    write(my_format,'(A2,I0,A1)')'(A',len(trim(WWW)),')'
    write(val_det,my_format) trim(WWW)
    
    end subroutine print_values_summary
end module martin_own
