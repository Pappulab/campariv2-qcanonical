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
! CONTRIBUTIONS: Rohit Pappu                                               !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!-----------------------------------------------------------------------
!
!
!             ##################################
!             #                                #
!             #  A SET OF SUBROUTINES WHICH    #
!             #  DEAL WITH INTERNAL COORDINATE #
!             #  SPACE (Z-MATRIX) AND MOLEC.   #
!             #  TOPOLOGY IN GENERAL           #
!             #                                #
!             ##################################
!
!
!-----------------------------------------------------------------------------
!
! this subroutine appends the Z-matrix by adding the atom specified
! according to the input coordinate values as well as the input ref.
! atoms
subroutine get_alternative_coord()
    ! Martin : this subroutine goal is to get alternative distances for the Hamiltonian  exchange 
    ! Algorithm, such that the distance used for the limit cases of solvation volume accessible fraction 
    ! , ammong others, can be correctly evaluated regardless of the replica. It thus aims at producing for each replicas 
    ! two set of coordinates corresponding to the initial and final state of the transformation.
    use atoms
    use zmatrix
    use molecule
    use params !(denbug)
    
    implicit none
    
    integer i2,i3,i4,chiral ,i,imol
    RTYPE bl,ba,baodi
    
    allocate(x_limits(size(x(:)),3))
    allocate(y_limits(size(y(:)),3))
    allocate(z_limits(size(z(:)),3))
   

    do imol=1,nmol
      do i=atmol(imol,1),atmol(imol,2) !Goes trouh all atoms
    
        i2 = iz(1,i)
        i3 = iz(2,i)
        i4 = iz(3,i)
        chiral = iz(4,i)
        bl = blen_limits(i,1)
        ba = bang_limits(i,1)
        
!        print *,bio_botyp(b_type(i)),bio_botyp(b_type(iz(1,i))),bio_botyp(b_type(iz(2,i)))
!        print *,i,bl,ba
        
        baodi = ztor(i)
        
        call genxyz(i,i2,bl,i3,ba,i4,baodi,chiral)
        
        x_limits(i,1)=x(i)
        y_limits(i,1)=y(i)
        z_limits(i,1)=z(i)
        
        x(:)=x_limits(:,2)
        y(:)=y_limits(:,2)
        z(:)=z_limits(:,2)
        
        bl = blen_limits(i,2)
        ba = bang_limits(i,2)
        baodi = ztor(i)
        
        call genxyz(i,i2,bl,i3,ba,i4,baodi,chiral)
        
        x_limits(i,2)=x(i)
        y_limits(i,2)=y(i)
        z_limits(i,2)=z(i)
        
        ! This is valid because the bond length and angle type do not change in HIS
        x_limits(i,3)=x(i)
        y_limits(i,3)=y(i)
        z_limits(i,3)=z(i)
        
      end do
    end do
end  
!
subroutine z_add(bti,z1,z2,z3,iz1,iz2,iz3,iz4)
!
  use atoms
  use iounit
  use params
  use zmatrix
!
  implicit none
!
  integer bti,iz1,iz2,iz3,iz4
  RTYPE z1,z2,z3
  
  if (bti.gt.n_biotyp) then
!
    write(ilog,*) 'Fatal. Requested unknown biotype in z_add(...)&
 &or sidechain(...). Offending type is ',bti,' (atom #',n,').'
    call fexit()
!
  else if (bti.gt.0) then
!
    attyp(n) = bio_ljtyp(bti)
    b_type(n) = bti
    if ((attyp(n).lt.0).OR.(attyp(n).gt.n_ljtyp)) then
      write(ilog,*) 'Fatal. Encountered unknown LJ-type in z_add(...&
 &). Offending type is ',attyp(n),' for atom #',n,'.'
      call fexit()
    else if (attyp(n).eq.0) then
      write(ilog,*) 'Fatal. Parameter file lacks support for biotype&
 & ',b_type(n),' (atom #',n,').'
      call fexit()
    else
      atnam(n) = lj_symbol(attyp(n))
      mass(n) = lj_weight(attyp(n))
    end if
!
    blen(n) = z1
    bang(n) = z2
    ztor(n) = z3
    ! So what I need to do is to check 
    
    if (ztor(n).lt.-180.0) then
      ztor(n) = ztor(n) + 360.0
    else if (ztor(n).gt.180.0) then
      ztor(n) = ztor(n) - 360.0
    end if
    
    iz(1,n) = iz1 ! This actually refers to the atoms number that are connected to n
    iz(2,n) = iz2 ! Normaly the first one should be self, but not the case here 
    iz(3,n) = iz3 ! 
    iz(4,n) = iz4 ! 
!   we have to increase the number of atoms by 1
    n = n + 1
!
! we use the -1-signal to indicate a ring closure
!
  else if (bti.eq.-1) then
     nadd = nadd + 1
     iadd(1,nadd) = iz1
     iadd(2,nadd) = iz2
!
  else

    write(ilog,*) 'Fatal. Received unsupported (less or equal to zero excepting -1) biotype in&
 & z_add(...). This is a bug (z_add).'
    call fexit()
!
  end if
!
end

subroutine  rewrite_params() ! The goal of this routine is to replace the value of the bonds angle, dihedrals and length
    ! from the ones hard wired into campari to the ones present in the parameter files
    ! I could also try to make that for the case where the hardwired parameters are used,.
    !but their is the issue  of hardcoding each values, which is bad practice as well as confusing and lots of work to implement 
    ! I do not think it is necessary at all, so I will just keep it as one of the things to do later. 
  use atoms
  use iounit
  use params
  ! Martin ; needed for sure :
  use zmatrix
  use inter
  use sequen
  use martin_own, only: test_combi
  use energies
  use movesets ! For hsq
  ! Martin : check there is probably some unuseful extras
  integer i,j,k,l,rs,h,m,r,last,skipped,reused,bo_curr,old
  integer added_ba,added_bo,mem_ba(1000,3),mem_bo(1000,3),changed_ba,changed_bo,added_temp,temp_new,temp_new2
  integer i1,i2,i3,i4,k1,k2,k3,k4,a
  integer, allocatable  :: mem_all(:,:)
  
  integer, allocatable :: array1(:),array2(:),bio_replace_pka(:,:)
  integer, allocatable :: used_ba(:),used_bo(:),indexes(:)
  integer, allocatable :: memory(:),mem_b_type_used(:,:)
  RTYPE, allocatable :: angle_target(:),bond_target(:),angle_current(:),bond_current(:)
  RTYPE, allocatable :: mem_added_angle(:,:),mem_added_bond(:,:)
  integer, allocatable :: mem_added_bo_type(:,:),mem_changed_b_type(:,:),angle_index(:),bond_index(:)
  integer, allocatable :: mem_btype_switch(:,:),used_b_type(:),bo_type_used(:),indexes_reused(:,:)
  integer, allocatable :: temp(:),divergence_bo(:,:)
  logical found
  integer ADDED_ANGLES,ADDED_BONDS,skip,added_di,added_tor,b,n_replace
  character my_format*100
  integer, allocatable:: b_types_used(:),bo_type_switch(:,:)
  
  integer, allocatable:: first_bond(:),first_angle(:),add_di(:,:),add_tor(:,:)

 
  
  RTYPE, allocatable :: test_temp_array1(:),test_temp_array2(:)
  RTYPE mem_new_ba(1000),mem_new_bo(1000),temp_angle,temp_bond
  logical test,need_fix
  
! The goal of all of this is to rewrite the parameter file bond definitions such that the bond angles and bond length potential are 
! Identical in the param file to the hard coded values.
  
! Currently, CAMPARI uses hardcoded values to construct the molecule, but parameter files values to compute the bond and angle energies
! This is not recomended for q-cannonical measurments, because the contribution to the free energy can result in non zero contributions
! for the bond energy. Thus to make sure the bond energy was not included in the calculation we set bond energies to zero by matching
! the length with the equilibrium position of the harmonic potential. 
! Note that for the moment, we assume that the bond energy between likely conformations in a meostates are similar. 
! This is likely to be a good assumption since ASP and GLU have similar bond, so do the phosphogroups, and the tautoers of HIS.
! This assumption also implies that the difference in pKa observed is solely due to other factors (mainly FOS, but also electroctatics and LJ)
! There may come a time when this assumption is proven wrong for some residues, but I doubt it.
! First, we identify all of the incorrect bonded type definitions (by comparing parameter file bond length and angles to topoloogy values),
! creating a new bo type for every new definitions, and create new anlge and bond values when needed while trying to reuse as much as possible
! Since the topology is unaffected by the dihedrals, we just copy exsiting definitions with new BOtype (only thos that apply to the system, as to not break them)
! There may be ways to improve this, but it does the job.
  
  
  
! This array is to replace the botype of interpolable residue based on the value of their "normal" counterpart
  n_replace=64
  allocate(bio_replace_pka(n_replace,2))
  
  bio_replace_pka(1,:)=(/1208,282 /)
  bio_replace_pka(2,:)=(/1209,1087 /)
  bio_replace_pka(3,:)=(/1210 ,239 /)
  bio_replace_pka(4,:)=(/ 1211,239/)
  bio_replace_pka(5,:)=(/1212 ,236 /)
  bio_replace_pka(6,:)=(/1213 ,238 /)
  bio_replace_pka(7,:)=(/ 1217,279/)
  bio_replace_pka(8,:)=(/1218, 280 /)
  bio_replace_pka(9,:)=(/1219,281/)
  bio_replace_pka(10,:)=(/ 1220,282/)
  bio_replace_pka(11,:)=(/ 1221 ,282/)
  bio_replace_pka(12,:)=(/1222,276 /)
  bio_replace_pka(13,:)=(/ 1224  ,  212/)
  bio_replace_pka(14,:)=(/1225  ,  214/)
  bio_replace_pka(15,:)=(/1226  ,  215/)
  bio_replace_pka(16,:)=(/1227  ,  215 /)
  bio_replace_pka(17,:)=(/ 1228  ,  1099/)
  bio_replace_pka(18,:)=(/1229  , 163 /)
  bio_replace_pka(19,:)=(/ 1230 ,   164/)
  bio_replace_pka(20,:)=(/1231  ,  165 /)
  bio_replace_pka(21,:)=(/1232  ,  166 /)
  bio_replace_pka(22,:)=(/ 1233 ,   167/)
  bio_replace_pka(23,:)=(/1234  ,  168 /)
  bio_replace_pka(24,:)=(/ 1235  ,  169/)
  bio_replace_pka(25,:)=(/1236  ,  170 /)
  bio_replace_pka(26,:)=(/1237  ,  171/)
  bio_replace_pka(27,:)=(/1238  ,  172/)
  bio_replace_pka(28,:)=(/ 1239   , 173/)
  bio_replace_pka(29,:)=(/1240  ,  163 /)
  bio_replace_pka(30,:)=(/ 1241  ,  164/)
  bio_replace_pka(31,:)=(/1242  ,  165  /)
  bio_replace_pka(32,:)=(/ 1243 ,   166/)
  bio_replace_pka(33,:)=(/ 1244   , 167 /)
  bio_replace_pka(34,:)=(/ 1245  ,  168/)
  bio_replace_pka(35,:)=(/ 1246  ,  169 /)
  bio_replace_pka(36,:)=(/ 1247 ,   170/)
  bio_replace_pka(37,:)=(/1248  ,  171/)
  bio_replace_pka(38,:)=(/1249  ,  172/)
  bio_replace_pka(39,:)=(/1250   , 173/)
  bio_replace_pka(40,:)=(/1292  ,  1274 /)
  bio_replace_pka(41,:)=(/1293  ,  1275/)
  bio_replace_pka(42,:)=(/1294 ,   1276 /)
  bio_replace_pka(43,:)=(/1295  ,  1277 /)
  bio_replace_pka(44,:)=(/1296  ,  1278 /)
  bio_replace_pka(45,:)=(/ 1297 ,  1278/)
  bio_replace_pka(46,:)=(/1298   , 1171 /)
  bio_replace_pka(47,:)=(/1299   , 1264 /)
  bio_replace_pka(48,:)=(/1300   , 1265 /)
  bio_replace_pka(49,:)=(/1301  ,  1266 /)
  bio_replace_pka(50,:)=(/1302  ,  1267/)
  bio_replace_pka(51,:)=(/1303   , 1267 /)
  bio_replace_pka(52,:)=(/ 1304  ,  1158/)
  bio_replace_pka(53,:)=(/1305  ,  1285/)
  bio_replace_pka(54,:)=(/1306  ,  1286/)
  bio_replace_pka(55,:)=(/1307  ,  1287 /)
  bio_replace_pka(56,:)=(/1308  ,  1290 /)
  bio_replace_pka(57,:)=(/ 1309  ,  1291/)
  bio_replace_pka(58,:)=(/1310 ,   1186 /)
  bio_replace_pka(59,:)=(/1311   , 1186 /)
  bio_replace_pka(60,:)=(/1312  ,  1136 /)
  bio_replace_pka(61,:)=(/ 1313  ,  1137/)
  bio_replace_pka(62,:)=(/ 1314  ,  1138/)
  bio_replace_pka(63,:)=(/1315 , 1139 /)
  bio_replace_pka(64,:)=(/1316 , 133 /)

  test=.true.
  
  need_fix=.false.
  found=.false.
  allocate(memory(3))
  allocate(array1(3))
  allocate(array2(3))
  allocate(test_temp_array1(size(blen(:))))
  allocate(test_temp_array2(size(bang(:))))
  
  allocate(mem_added_bo_type(100000,5))
  allocate(mem_added_angle(100000,3))
  allocate(mem_added_bond(100000,3))
  allocate(mem_changed_b_type(100000,3))

  allocate(mem_btype_switch(n_biotyp,3))

  mem_btype_switch(:,:)=0
  
  test_temp_array1(:)=blen(:)
  test_temp_array2(:)=bang(:)
  
  added_ba=0
  added_bo=0
  
  changed_ba=0
  changed_bo=0
  
  mem_ba(:,:)=0
  mem_bo(:,:)=0
  
  mem_new_ba(:)=0
  mem_new_bo(:)=0
  
  
  allocate(mem_all(n,3))
  allocate(angle_target(n))
  allocate(angle_index(n))
  allocate(bond_index(n))
  allocate(mem_b_type_used(n,3))

  mem_all(:,:)=0
  angle_target(:)=0.
  angle_index(:)=0
  bond_index(:)=0
  
 mem_b_type_used(:,:)=0
  
  skipped=0
  ! So the revised goal is : 
  ! I identify biotype that need a change in 

  do i=1,n
    if ((iz(1,i).gt.0).and.(iz(2,i).gt.0)) then 
        array1(1)=bio_botyp(b_type(i))
        array1(2)=bio_botyp(b_type(iz(1,i)))
        array1(3)=bio_botyp(b_type(iz(2,i)))
    else
        cycle 
    end if 

    do k=1,ba_lstsz 
        array2(1)=ba_lst(k,1)
        array2(2)=ba_lst(k,2)
        array2(3)=ba_lst(k,3)

        test=.true.
        memory(:)=test_combi(array1,array2,size(memory(:)))

        do j=1,size(memory(:))
            if (memory(j).eq.0) then
                test=.false.
            end if 
        end do 
        if (test.eqv..true.) then 

            if (array2(2).ne.array1(2)) then 
                test=.false.
            end if 
        end if

        
        if (test.eqv..true.) then

            angle_index(i)=ba_lst(k,4)
            exit
        end if 
    end do 
  end do 
  deallocate(memory)
  deallocate(array1)
  deallocate(array2)

  allocate(memory(2))
  allocate(array1(2))
  allocate(array2(2))

! Second part, bond length
  do i=1,n
    if (iz(1,i).gt.0) then 
        array1(1)=bio_botyp(b_type(i))
        array1(2)=bio_botyp(b_type(iz(1,i)))
    else
        cycle 
    end if 

    do k=1,bo_lstsz      
        array2(:)=bo_lst(k,1:2)
        test=.true.
        memory(:)=test_combi(array1,array2,size(memory(:)))
        do j=1,size(memory(:))
            if (memory(j).eq.0) then
                test=.false.
            end if 
        end do 
        if (test.eqv..true.) then
            bond_index(i)=bo_lst(k,3)          
           exit 
        end if 
    end do    
  end do 
! Bond parameter check
allocate(used_b_type(n_biotyp))

used_b_type(:)=0
! This is a 
old=0
added_bo=0
added_angles=0
added_bonds=0
mem_added_bo_type(:,:)=0 !@ CHanged : ger the two connected atoms, then the bond type, then angle type

!mem_added_ba_type=(:,:)=0! Added bond_type_angle definition
mem_added_angle(:,:)=0. ! Added corresponding angle
mem_added_bond(:,:)=0. !added corresponding bond
mem_changed_b_type(:,:)=0 ! Memory of the modification to be made to biotype

deallocate(memory)
deallocate(array1)
deallocate(array2)

allocate(memory(3))
allocate(array1(3))
allocate(array2(3))

allocate(first_bond(n))
allocate(first_angle(n))
first_bond(:)=0
first_angle(:)=0


!Note that I am ALWAYS going to add a bo_type if the angle is incorrect
! Last change : I cannot use a similar bo_type for things that were precendently different 
do i=1,n
    if ((bond_index(i).eq.0)) then ! .or.(angle_index(i).eq.0)) then 
        temp_bond =0
        
    else 
        temp_bond =bo_par(bond_index(i),2)
    end if 
!    if (mem_btype_switch(b_type(i),1).ne.0) then 
!        cycle  
!    end if    
    bo_curr=0
    found=.False.
    do k=1,added_bo
        if (mem_changed_b_type(k,1).eq.b_type(i)) then 
            found=.true.
            bo_curr=mem_changed_b_type(k,2)
            exit
        end if 
    end do 

!    if (found) then
!        cycle
!    end if 
    ! This is to avoid unallocated references
    if (angle_index(i).eq.0) then 
        temp_angle=0.
    else 
        temp_angle=ba_par(angle_index(i),2)
    end if 
    
        ! That means something needs to be changed
    if (((abs(blen(i)-temp_bond).gt.0.00001).or.((abs(bang(i)-temp_angle).gt.0.00001))).and.(iz(1,i).ne.0)) then !.and.((angle_index(i).ne.0))))then 
!        print *,"For atom",i 
!        print *,blen(i),temp_bond
!        print *,bang(i),temp_angle
        array1(1)=bio_botyp(b_type(i))
        
        array1(2)=bio_botyp(b_type(iz(1,i)))
!        if (iz(2,i).ne.0) then 
!            array1(3)=bio_botyp(b_type(iz(2,i)))
!            print *,i,iz(1,i),iz(2,i)
!            print *,b_type(i),b_type(iz(1,i)),b_type(iz(2,i))
!        else 
!            array1(3)=0
!            print *,i,iz(1,i),iz(2,i)
!            print *,b_type(i),b_type(iz(1,i)),0
!        end if 
!        
!        print *,array1(:)
!        print *,blen(i),temp_bond,bang(i),temp_angle
        

        ! Now check if any of the newly created types has the right values.
        if (added_bo.ne.0) then 
            found=.false.
            do j=1,added_bo
                    ! If angles and bonds are the same, reuse
                    if (((abs(bang(i)-mem_added_angle(j,3)).lt.0.00001).and.(abs(blen(i)-mem_added_bond(j,3)).lt.0.00001))) then                         
                        ! If not already changed
!                        print *,bang(i),mem_added_angle(j,3)
!                        print *,blen(i),mem_added_bond(j,3)
!                        print *,"GGGG",mem_changed_b_type(j,:)
!                        print *,mem_btype_switch(mem_changed_b_type(j,1),1),bio_botyp(b_type(i))
                        added_bo=added_bo+1
                        ! This is to make sure previously different botypes get different new BOtype, even with same angles 
                        !(important for duplication of dihedrals))
                        ! Should make an exception for  if no definitions overlap
                        
                        if (mem_btype_switch(mem_changed_b_type(j,1),1).eq.bio_botyp(b_type(i))) then 
                        
                            mem_added_bo_type(added_bo,:)=mem_added_bo_type(j,:)
                            
                            mem_added_angle(added_bo,:)=mem_added_angle(j,:)
                            mem_added_bond(added_bo,:)=mem_added_bond(j,:)
                            
                            
                            mem_btype_switch(b_type(i),1)=bio_botyp(b_type(i))
                            mem_btype_switch(b_type(i),2)=mem_changed_b_type(j,2)                   
                            
                            mem_changed_b_type(added_bo,1)=b_type(i)
                            mem_changed_b_type(added_bo,2)=mem_changed_b_type(j,2)
                            skipped=skipped+1
                            
                            mem_b_type_used(added_bo,1)=b_type(i)                     
                            mem_b_type_used(added_bo,2)=b_type(iz(1,i))
                            
                            if (iz(2,i).ne.0) then 
                                mem_b_type_used(added_bo,3)=b_type(iz(2,i))
                            end if 
                            
                            found=.true.
                            exit
                            
                        else 
                           
                            !mem_added_bo_type(added_bo,1)=n_botyp+added_bo-skipped
                            mem_added_bo_type(added_bo,2)=bio_botyp(b_type(iz(1,i)))
                            mem_added_bo_type(added_bo,4)=mem_added_bo_type(j,4)
                            if (iz(2,i).ne.0) then 
                                mem_added_bo_type(added_bo,3)=bio_botyp(b_type(iz(2,i)))
                                mem_added_bo_type(added_bo,5)=mem_added_bo_type(j,5)
                            end if 
                            
                            
                            
                            mem_added_angle(added_bo,:)=mem_added_angle(j,:)
                            mem_added_bond(added_bo,:)=mem_added_bond(j,:)

                            if (bo_curr.eq.0) then 
                                mem_added_bo_type(added_bo,1)=n_botyp+added_bo-skipped
                                mem_btype_switch(b_type(i),1)=bio_botyp(b_type(i))
                                mem_btype_switch(b_type(i),2)=n_botyp+added_bo-skipped

                                mem_changed_b_type(added_bo,1)=b_type(i)
                                mem_changed_b_type(added_bo,2)=n_botyp+added_bo-skipped
                            else 
                                mem_added_bo_type(added_bo,1)=bo_curr
                                mem_btype_switch(b_type(i),1)=bio_botyp(b_type(i))
                                mem_btype_switch(b_type(i),2)=bo_curr

                                mem_changed_b_type(added_bo,1)=b_type(i)
                                mem_changed_b_type(added_bo,2)=bo_curr
                            end if 
                            
                            mem_b_type_used(added_bo,1)=b_type(i)
                            mem_b_type_used(added_bo,2)=b_type(iz(1,i))
                            if (iz(2,i).ne.0) then 
                                mem_b_type_used(added_bo,3)=b_type(iz(2,i))
                            end if 
                            found=.true.
                            exit
                            
                            
                            
                        end if 
              
                    end if 
!                end if 
            end do
            !if you haven't found a match, create one 
            if (found.eqv..false.) then 
                if (bo_curr.eq.0) then 
                    added_bo=added_bo+1

                    mem_added_bo_type(added_bo,1)=n_botyp+added_bo-skipped
                else 
                    added_bo=added_bo+1

                    mem_added_bo_type(added_bo,1)=bo_curr
                    skipped=skipped+1
                end if 

                mem_added_bo_type(added_bo,2)=bio_botyp(b_type(iz(1,i)))
                if (iz(2,i).ne.0) then 
                mem_added_bo_type(added_bo,3)=bio_botyp(b_type(iz(2,i)))
                end if 
                
                mem_b_type_used(added_bo,1)=b_type(i)
                mem_b_type_used(added_bo,2)=b_type(iz(1,i))
                if (iz(2,i).ne.0) then 
                    mem_b_type_used(added_bo,3)=b_type(iz(2,i))
                end if 
                found=.false.
                ! Now check if we can refer to an already existing angle
                do k=1,n_angltyp
                    if (abs(ba_par(k,2)-bang(i)).lt.0.0001) then 
                                                    
                        mem_added_angle(added_bo,1)=k
                        mem_added_angle(added_bo,2)=100.0
                        mem_added_angle(added_bo,3)=ba_par(k,2)
                        !mem_added_angle(added_angles,1)=k
                        found=.true.
                        exit
                    end if 
                end do 
                
                ! That is one of the bond values corresponds to the desired value
                if (found) then 

                    mem_added_bo_type(added_bo,5)=k
                else ! if none
                !    print *,"DD"                    
                    do k=1,added_bo
                        if ((abs(mem_added_angle(k,3)-bang(i)).lt.0.0001).and.(bang(i).ne.0.)) then 

                            mem_added_bo_type(added_bo,5)=mem_added_bo_type(k,5)
                            
                            mem_added_angle(added_bo,1)=mem_added_angle(k,1)
                            mem_added_angle(added_bo,2)=mem_added_angle(k,2)
                            mem_added_angle(added_bo,3)=mem_added_angle(k,3)
                            exit 
                        end if 
                    end do 
                    
                    if ((k.eq.added_bo+1).and.(bang(i).ne.0.)) then 

                        added_angles=added_angles+1
                        mem_added_bo_type(added_bo,5)=n_angltyp+added_angles
                        first_angle(added_bo)=1
                        mem_added_angle(added_bo,1)=n_angltyp+added_angles
                        mem_added_angle(added_bo,2)=100.00 ! The spring constant, pretty irelevant for our purpose
                        mem_added_angle(added_bo,3)=bang(i)
                    end if 
                end if 
                
                ! Now for the bonds
                ! Now check if we can refer to an already existing angle
                do k=1,n_bondtyp 
                    if (abs(bo_par(k,2)-blen(i)).lt.0.0001) then 
                        mem_added_bond(added_bo,1)=k
                        mem_added_bond(added_bo,2)=100.0
                        mem_added_bond(added_bo,3)=bo_par(k,2)

                        exit
                    end if 
                end do 
                
                ! That is one of the bond values corresponds to the desired value
                if (k.ne.n_bondtyp+1) then 

                    mem_added_bo_type(added_bo,4)=k
                else ! if none
                    
                    do k=1,added_bo
                        if (abs(mem_added_bond(k,3)-blen(i)).lt.0.0001) then 
                            mem_added_bond(added_bo,1)=mem_added_bond(k,1)
                            mem_added_bond(added_bo,2)=mem_added_bond(k,2)
                            mem_added_bond(added_bo,3)=mem_added_bond(k,3)
                            mem_added_bo_type(added_bo,4)=mem_added_bo_type(k,4)
                            exit 
                        end if 
                    end do 
                    
                    if (k.eq.added_bo+1) then 

                        added_bonds=added_bonds+1
                        first_bond(added_bo)=1
                        mem_added_bo_type(added_bo,4)=n_bondtyp+added_bonds
                        
                        mem_added_bond(added_bo,1)=n_bondtyp+added_bonds
                        mem_added_bond(added_bo,2)=50 ! The spring constant, pretty irelevant for our purpose
                        mem_added_bond(added_bo,3)=blen(i)
                    end if 

                end if 
                if (bo_curr.eq.0) then 
                    mem_btype_switch(b_type(i),1)=bio_botyp(b_type(i))
                    mem_btype_switch(b_type(i),2)=n_botyp+added_bo-skipped
                
                    mem_changed_b_type(added_bo,1)=b_type(i)
                    mem_changed_b_type(added_bo,2)=n_botyp+added_bo-skipped
                else 
                    mem_btype_switch(b_type(i),1)=bio_botyp(b_type(i))
                    mem_btype_switch(b_type(i),2)=bo_curr
                
                    mem_changed_b_type(added_bo,1)=b_type(i)
                    mem_changed_b_type(added_bo,2)=bo_curr
                end if 
            end if 
        else ! This is to make sure we do not add anything that unecessary
            
            ! That is only for the first element added
            added_bo=added_bo+1

            mem_added_bo_type(added_bo,1)=n_botyp+added_bo-skipped
            mem_added_bo_type(added_bo,2)=bio_botyp(b_type(iz(1,i)))
            
            mem_b_type_used(added_bo,1)=b_type(i)
            mem_b_type_used(added_bo,2)=b_type(iz(1,i))
            if (iz(2,i).ne.0) then 
            mem_b_type_used(added_bo,3)=b_type(iz(2,i))
            end if 
            if (temp_angle.ne.0.) then 
                mem_added_bo_type(added_bo,3)=bio_botyp(b_type(iz(2,i)))
            else 
                mem_added_bo_type(added_bo,3)=0
            end if 
            found=.false.
            ! Now check if we can refer to an already existing angle
            do k=1,n_angltyp
                if ((abs(ba_par(k,2)-bang(i)).lt.0.0001).and.(bang(i).ne.0.)) then 
                        mem_added_angle(added_bo,1)=k
                        mem_added_angle(added_bo,2)=100.0
                        mem_added_angle(added_bo,3)=ba_par(k,2)
                    found=.true.
                    exit
                end if 
            end do 
            ! That is one of the bond values corresponds to the desired value
            if (found) then 


                mem_added_bo_type(added_bo,5)=k
            else if (bang(i).ne.0.) then ! if none
                added_angles=added_angles+1
                mem_added_bo_type(added_bo,5)=n_angltyp+added_angles
                
                mem_added_angle(added_bo,1)=n_angltyp+added_angles
                first_angle(added_bo)=1
                if (temp_angle.ne.0.) then 
                    mem_added_angle(added_bo,2)=100.
                else 
                    mem_added_angle(added_bo,2)=0.
                end if 
                mem_added_angle(added_bo,3)=bang(i)
            end if 

            ! Now for the bonds
            ! Now check if we can refer to an already existing angle
            do k=1,n_bondtyp
                if (abs(bo_par(k,2)-blen(i)).lt.0.0001) then 
                    mem_added_bond(added_bo,1)=k
                    mem_added_bond(added_bo,2)=100.0
                    mem_added_bond(added_bo,3)=bo_par(k,2)

                    exit
                end if 
            end do 
            ! That is one of the bond values corresponds to the desired value
            if (k.ne.n_bondtyp+1) then 


                mem_added_bo_type(added_bo,4)=k
            else ! if none


                added_bonds=added_bonds+1
                first_bond(added_bo)=1

                mem_added_bo_type(added_bo,4)=n_bondtyp+added_bonds
                mem_added_bond(added_bo,1)=n_bondtyp+added_bonds
                mem_added_bond(added_bo,2)=50.
                mem_added_bond(added_bo,3)=blen(i)
            end if 
                if (bo_curr.eq.0) then 
                    mem_btype_switch(b_type(i),1)=bio_botyp(b_type(i))
                    mem_btype_switch(b_type(i),2)=n_botyp+added_bo-skipped
                
                    mem_changed_b_type(added_bo,1)=b_type(i)
                    mem_changed_b_type(added_bo,2)=n_botyp+added_bo-skipped
                else 
                    mem_btype_switch(b_type(i),1)=bio_botyp(b_type(i))
                    mem_btype_switch(b_type(i),2)=bo_curr
                
                    mem_changed_b_type(added_bo,1)=b_type(i)
                    mem_changed_b_type(added_bo,2)=bo_curr
                end if 
            
        end if   

    else if (iz(2,i).ne.0) then ! That is second atom of zmat for molecule
        old=old+1
        added_bo=added_bo+1
        mem_added_bo_type(added_bo,1)=bio_botyp(b_type(i))
        mem_added_bo_type(added_bo,2)=bio_botyp(b_type(iz(1,i)))

        
        mem_b_type_used(added_bo,1)=b_type(i)
        mem_b_type_used(added_bo,2)=b_type(iz(1,i))
        
        
        if (iz(2,i).ne.0) then 
            mem_b_type_used(added_bo,3)=b_type(iz(2,i))
            mem_added_bo_type(added_bo,3)=bio_botyp(b_type(iz(2,i)))
        else 
            mem_b_type_used(added_bo,3)=0
            mem_added_bo_type(added_bo,3)=0
        end if 
        mem_added_bo_type(added_bo,4)=bond_index(i)
        mem_added_bo_type(added_bo,5)=angle_index(i)
        
        
        array1(1)=bio_botyp(b_type(i))
        array1(2)=bio_botyp(b_type(iz(1,i)))

    else  if (bio_botyp(b_type(i)).ne.0) then ! As to not apply to salt
        ! This only applies to and iz1=0 and iz2 =0
        ! For multiple chains, we check this btype has not been used before
        if (mem_btype_switch(b_type(i),2).eq.0) then         
            added_bo=added_bo+1
            old=old+1
            mem_btype_switch(b_type(i),1)=bio_botyp(b_type(i))
            mem_btype_switch(b_type(i),2)=n_botyp+added_bo-skipped
        end if 

    end if  
    ! At that point we already have changed what needs change
    !print *,mem_added_bo_type(added_bo,:)
    used_b_type(b_type(i))=1
end do

!do i=1,added_bo
!    do j=1,added_bo
!        if ((mem_changed_b_type(i,1).eq.mem_changed_b_type(j,1)).and.(mem_changed_b_type(i,2).ne.mem_changed_b_type(j,2))) then 
!            print *,"Problem",i
!            print *,mem_changed_b_type(i,1),mem_changed_b_type(j,1),mem_changed_b_type(i,2),mem_changed_b_type(j,2)
!            
!        end if 
!    end do 
!end do 





allocate(bo_type_used(n_botyp+added_bo))
bo_type_used(:)=0
last=0
k1=0
allocate(temp(n_biotyp))
temp(:)=bio_botyp(:)

  ! Check each biotypes
do k=1,n_biotyp
    do i=1,added_bo
        ! If right biotype, switch bo_type, mark as used
        if (mem_changed_b_type(i,1).eq.k) then
            bio_botyp(k)=mem_changed_b_type(i,2)
            if (bio_botyp(k).gt.last) then 
                
                last=bio_botyp(k)
                
            end if 
!            if (mem_changed_b_type(i,2).ne.0) then 
!                bo_type_used(mem_changed_b_type(i,2))=1
!            end if
            exit
        end if 
    end do 
    
    ! Now go throught old bo_type, and mark if used
!    do i=1,n_botyp
!        if (i.eq.bio_botyp(k)) then
!            bo_type_used(i)=1
!            exit
!        end if 
!    end do 
      
!    if  (bio_botyp(k).ne.0) then 
!        bo_type_used(bio_botyp(k))=1
!    end if 
end do 

! THis is just to automatically adapt the pKa interpolable versio opf the residues
do i=1,n_replace    
    bio_botyp(bio_replace_pka(i,1))=bio_botyp(bio_replace_pka(i,2))
end do 



!do k=1,n_biotyp
!    do i=1,added_bo
!        !This is just to be in the right order
!        if (mem_changed_b_type(i,1).eq.k) then
!            last=last+1
!            bio_botyp(k)=mem_changed_b_type(i,2)
!            
!            if (mem_changed_b_type(i,2).ne.0) then 
!                bo_type_used(mem_changed_b_type(i,2)-n_botyp)=1
!            end if
!        end if 
!    end do 
!    
!    if (last.ne.k) then 
!        last=last+1
!        
!
!        if  (bio_botyp(k).ne.0) then 
!            bo_type_used(bio_botyp(k))=1
!        end if 
!    end if
!end do 
444   format(A8,I4,A,A,A,A,A,I3,A,I3,A,I3)
  
!  allocate(b_types_used(n_botyp+added_bo+100))
!  b_types_used(:)=0

  do i=1,n_biotyp
      if (mem_btype_switch(i,1).ne.0) then
          bo_type_used(mem_btype_switch(i,2))=1
      else 
          if (bio_botyp(i).ne.0) then 
                bo_type_used(bio_botyp(i))=1
          end if 
      end if 
  end do 
  
  

 ! I should try reorganizing the bonded types, to avoid making dummies

  reused=0
  skip=0
  k1=0 
   445   format(A8,I4,A,I3,A,I3,A,I3)
   ! Check how many types are used

   do i=1,last!n_botyp+added_bo-skipped
       if (bo_type_used(i).eq.0) then 
           skip=skip+1
       end if 
   end do 
   ! Now we know we need to change the last-k1 highest used botype to lower ones
   !skip=last-k1
!   allocate(indexes_reused(skip))
 
   allocate(indexes_reused(skip,2))

   indexes_reused(:,:)=0

   ! Get which used biotypes are the skip highest
   reused=0
   do i=last,1,-1
        if ((bo_type_used(i).eq.1).and.(reused.lt.skip)) then 
            reused=reused+1
            indexes_reused(reused,1)=i
        end if 
   end do 
     reused=0
     k1=0
      ! Get which unused biotypes are the skip lowest 
   do i=1,last
       if ((bo_type_used(i).eq.0).and.(reused.lt.skip)) then 
            if (indexes_reused(reused+1,1).ge.i) then ! (i.le.n_botyp)) then 
                reused=reused+1
                indexes_reused(reused,2)=i
            end if 
       end if 
   end do 
   ! Needed 
   !reused=reused+1
   !indexes_reused(reused+1,2)=i
   do i=1,last
       print *,i,bo_type_used(i)
   end do 
!   
!   print *,n_botyp,added_bo,skipped
!   print *,"skip",skip
!   print *,last,k1
   
   ! Now indexes should have the original and final botypes. It should also delete any mention of the biotypes in previous definitions

!       !print *,i,b_types_used(i)
!!       if (b_types_used(i).eq.0) then 
!           ! This
!           if (i.gt.n_botyp+added_bo-reused-skip+1) then 
!               cycle 
!           end if 
!           reused=reused+1
!
!           do while (b_types_used(n_botyp+added_bo-reused-skip+1).eq.0)
!
!               skip=skip+1
!           end do
!           if (n_botyp+added_bo-reused-skip+1.gt.i) then 
!                indexes_reused(reused,1)=n_botyp+added_bo-reused-skip+1
!                indexes_reused(reused,2)=i
!           end if 
!!       end if 
!   end do 
!  
!  allocate(indexes_reused(added_bo,2))
!  indexes_reused(:,:)=0
!  reused=0
!  skip=0
!  last =n_biotyp
!   445   format(A8,I4,A,I3,A,I3,A,I3)
!   do i=1,n_botyp+added_bo
!       !print *,i,b_types_used(i)
!       if (b_types_used(i).eq.0) then 
!           ! This
!           if (i.gt.n_botyp+added_bo-reused-skip+1) then 
!               cycle 
!           end if 
!           reused=reused+1
!
!           do while (b_types_used(n_botyp+added_bo-reused-skip+1).eq.0)
!
!               skip=skip+1
!           end do
!           if (n_botyp+added_bo-reused-skip+1.gt.i) then 
!                indexes_reused(reused,1)=n_botyp+added_bo-reused-skip+1
!                indexes_reused(reused,2)=i
!           end if 
!       end if 
!   end do 

 do i=1,n_replace    
     do j=1,reused
        if (bio_botyp(bio_replace_pka(i,2)).eq.indexes_reused(j,1)) then 
            bio_botyp(bio_replace_pka(i,1))=  indexes_reused(j,2)   ! =bio_botyp(bio_replace_pka(i,2))
        end if 
     end do 
 end do  
   
   
   do i=1,added_bo
       do j=1,reused
           if (indexes_reused(j,1).eq.mem_added_bo_type(i,1)) then 
               mem_added_bo_type(i,1)=indexes_reused(j,2)
           end if
           if (indexes_reused(j,1).eq.mem_added_bo_type(i,2)) then 
               mem_added_bo_type(i,2)=indexes_reused(j,2)
           end if
           if (indexes_reused(j,1).eq.mem_added_bo_type(i,3)) then 
               mem_added_bo_type(i,3)=indexes_reused(j,2)
           end if
       end do 
   end do 

   do i=1,added_bo
       do j=1,reused
           if (indexes_reused(j,1).eq.mem_changed_b_type(i,2)) then 
               mem_changed_b_type(i,2)=indexes_reused(j,2)
           end if
       end do 
   end do 
   
   do i=1,n_biotyp
       do j=1,reused
           if (indexes_reused(j,1).eq.mem_btype_switch(i,2)) then 
               mem_btype_switch(i,2)=indexes_reused(j,2)
           end if
       end do 
   end do 
   
   print *,"Start of reused"
   do j=1,skip
    print *,indexes_reused(j,:)
   end do 

   
   
!bo_type_used(:)=0
last=0

do k=1,n_biotyp
    if (mem_btype_switch(k,2).ne.0) then 
        write(ilog, 444 ) "biotype ",k,"  ",bio_code(k),'  "',bio_res(k),'"  ',&
        &bio_ljtyp(k),"  ",bio_ctyp(k),"  ",mem_btype_switch(k,2)
    else 
        write(ilog, 444 ) "biotype ",k,"  ",bio_code(k),'  "',bio_res(k),'"  ',&
        &bio_ljtyp(k),"  ",bio_ctyp(k),"  ",bio_botyp(k)
    end if
end do 

! Now I need to change the angle definition for each types that were charged
! For angles, we need a double loop though biotypes 
added_temp=added_bo

deallocate(memory)
deallocate(array1)
deallocate(array2)

allocate(memory(3))
allocate(array1(3))
allocate(array2(3))

!n_biotyp= last 

! Here just changing old definitions to new ones

do i=1,added_bo
    do j=1,3
        if (mem_b_type_used(i,j).eq.0) cycle
        if (mem_btype_switch(mem_b_type_used(i,j),1).ne.0) then 
            mem_added_bo_type(i,j)=mem_btype_switch(mem_b_type_used(i,j),2)
        end if
    end do 
end do 

!Here I need to duplicate 
allocate(divergence_bo(n_botyp,50)) ! First dimension is the number of divergent botype
divergence_bo(:,:)=0


found=.false.
do i=1,n_botyp
    do j=1,n_biotyp
        if (mem_btype_switch(j,1).eq.i) then 
            do k=2,100
                if (divergence_bo(i,k).eq.mem_btype_switch(j,2)) then ! If already found
                    exit
                else if (divergence_bo(i,k).eq.0) then
                    divergence_bo(i,1)=divergence_bo(i,1)+1
                    
                    divergence_bo(i,k)=mem_btype_switch(j,2)
                    exit
                end if 
            end do 
        end if 
    end do 
end do 

deallocate(array1)
deallocate(array2)
deallocate(memory)

allocate(array1(2))
allocate(array2(2))
allocate(memory(2))


do i=1,bo_lstsz ! Also, should also be 
    do j=2,divergence_bo(bo_lst(i,1),1)+1
        do k=2,divergence_bo(bo_lst(i,2),1)+1
            
            if ((divergence_bo(bo_lst(i,1),1).eq.1).and.(divergence_bo(bo_lst(i,2),1).eq.1)) then 
                !  Not a real change 
                 old =old+1
            end if 
            added_bo=added_bo+1
            mem_added_bo_type(added_bo,1)=divergence_bo(bo_lst(i,1),j)
            mem_added_bo_type(added_bo,2)=divergence_bo(bo_lst(i,2),k)
            mem_added_bo_type(added_bo,4)=bo_lst(i,3)
        end do 
    end do 
end do 

!do i=1,bo_lstsz
!    do j=1,added_bo
!        found=.false.
!        array1(:)=bo_lst(i,1:2)
!        array2(:)=mem_added_bo_type(j,1:2)
!        memory(:)=test_combi(array1,array2,size(memory(:)))
!        if ((memory(1).ne.0).and.(memory(2).ne.0)) then
!          
!            found=.true.
!            exit
!        end if 
!    end do 
!    if (found.eqv..false.) then 
!        print *,"bonded_type_bond",bo_lst(i,:)
!    end if 
!end do 
   
do i=1,bo_lstsz
    found=.false.
    do j=1,reused
        if ((bo_lst(i,1).eq.indexes_reused(j,2)).or.(bo_lst(i,2).eq.indexes_reused(j,2))) then 
            found=.true.
            exit
        end if  
    end do 
    if (found.eqv..false.) then 
        print *," bonded_type_bond",bo_lst(i,:)
    else ! This is just to get a memory of the printed
        bo_lst(i,:)=0
    end if 
end do 

 do i=1,added_bo
     !
     i1=mem_added_bo_type(i,1)
     i2=mem_added_bo_type(i,2)
     found=.false.
     do k=1,bo_lstsz
     

        if ((((i1.eq.bo_lst(k,1)).and.(i2.eq.bo_lst(k,2))).OR.&
          &((i2.eq.bo_lst(k,1)).and.(i1.eq.bo_lst(k,2))))) then 
          found=.true.
          exit
          end if 
     end do
     if (found.eqv..true.) cycle
     if (mem_added_bo_type(i,4).eq.0) cycle 
     found=.false.
     do j=1,i-1
         k1=mem_added_bo_type(j,1)
         k2=mem_added_bo_type(j,2)
        if (mem_added_bo_type(j,4).eq.0) cycle 
          if ((((i1.eq.k1).and.(i2.eq.k2)).OR.&
          &((i1.eq.k2).and.(i2.eq.k1)))) then 
            found=.true.
            exit
          end if 
     end do
     if (found.eqv..false.) then 
        print *,"bonded_type_bond",mem_added_bo_type(i,1:2),mem_added_bo_type(i,4)
     end if 
 end do 


deallocate(array1)
deallocate(array2)
deallocate(memory)

allocate(array1(3))
allocate(array2(3))
allocate(memory(3))


!added_di=0
!Finally, I need to report the changes in biotypes for all other pair not used in topology, but only useed for potentials

do i=1,ba_lstsz
    do j=2,divergence_bo(ba_lst(i,1),1)+1
        do k=2,divergence_bo(ba_lst(i,2),1)+1
            do l=2,divergence_bo(ba_lst(i,3),1)+1
                if ((divergence_bo(ba_lst(i,1),1).eq.1).and.(divergence_bo(ba_lst(i,2),1).eq.1)&
                &.and.(divergence_bo(ba_lst(i,3),1).eq.1))  then 
                    !  Not a real change 
                    old =old+1
                end if 
                added_bo=added_bo+1
                mem_added_bo_type(added_bo,1)=divergence_bo(ba_lst(i,1),j)
                mem_added_bo_type(added_bo,2)=divergence_bo(ba_lst(i,2),k)
                mem_added_bo_type(added_bo,3)=divergence_bo(ba_lst(i,3),l)
                mem_added_bo_type(added_bo,5)=ba_lst(i,4)
            end do      
        end do 
    end do 
end do 

!do i=1,ba_lstsz
!    found=.false.
!    do j=1,added_bo
!        array1(:)=ba_lst(i,1:3)
!        array2(:)=mem_added_bo_type(j,1:3)
!        memory(:)=test_combi(array1,array2,size(memory(:)))
!        if ((memory(1).ne.0).and.(array2(2).eq.array1(2)).and.(memory(3).ne.0)) then
!            found=.true.
!            exit
!        end if 
!    end do 
!    if (found.eqv..false.) then 
!        print *,"bonded_type_angle",ba_lst(i,:)
!    end if  
!end do 

do i=1,ba_lstsz
    found=.false.
    do j=1,reused
        if ((ba_lst(i,1).eq.indexes_reused(j,2)).or.(ba_lst(i,2).eq.indexes_reused(j,2)).or.&
        &(ba_lst(i,3).eq.indexes_reused(j,2))) then 
            found=.true.
            exit
        end if  
    end do 
    if (found.eqv..false.) then 
        print *,"bonded_type_angle",ba_lst(i,:)
    else  
        ba_lst(i,:)=0
    end if 
end do 

 do i=1,added_bo
     i1=mem_added_bo_type(i,1)
     i2=mem_added_bo_type(i,2)
     i3=mem_added_bo_type(i,3)
     found=.false.
     do k=1,ba_lstsz
          if ((((i1.eq.ba_lst(k,1)).and.(i2.eq.ba_lst(k,2)).and.(i3.eq.ba_lst(k,3))).OR.&
          &((i1.eq.ba_lst(k,3)).and.(i2.eq.ba_lst(k,2)).and.(i3.eq.ba_lst(k,1))))) then 
            found=.true.
            exit
          end if 
     end do 
     if (found.eqv..true.) cycle
    if (mem_added_bo_type(i,5).eq.0) cycle 
     found=.false.
     do j=1,i-1
         k1=mem_added_bo_type(j,1)
         k2=mem_added_bo_type(j,2)
         k3=mem_added_bo_type(j,3)
         if (mem_added_bo_type(j,5).eq.0) cycle 
          if ((((i1.eq.k1).and.(i2.eq.k2).and.(i3.eq.k3)).OR.&
          &((i1.eq.k3).and.(i2.eq.k2).and.(i3.eq.k1)))) then 
            found=.true.
            exit
          end if 
     end do
     if (found.eqv..false.) then 
        print *,"bonded_type_angle",mem_added_bo_type(i,1:3),mem_added_bo_type(i,5)
     end if 
 end do 
 
added_di=0

!Finally, I need to report the changes in biotypes for all other pair not used in topology, but only useed for potentials

!!!!!!!!!!!!!!
! So that problem as it manifest itself curently is not with the definition per sey it seems. (except for the CAP-glycine)



allocate(add_di(10000000,5))

add_di(:,:)=0
do rs=1,nseq
    do a=1,nrsdi(rs)
        i1 = temp(b_type(iaa(rs)%di(a,1)))
        i2 = temp(b_type(iaa(rs)%di(a,2)))
        i3 = temp(b_type(iaa(rs)%di(a,3)))
        i4 = temp(b_type(iaa(rs)%di(a,4)))
        do i=1,di_lstsz
            k1 = di_lst(i,1)
            k2 = di_lst(i,2)
            k3 = di_lst(i,3)
            k4 = di_lst(i,4)

           ! This is to not add a definition if there is not in the normal verison
            if (.NOT.((((i2.eq.k2).AND.(i3.eq.k3)).AND.&
 &                 ((i1.eq.k1).AND.(i4.eq.k4))).OR.&
 &                (((i2.eq.k3).AND.(i3.eq.k2)).AND.&
 &                 ((i1.eq.k4).AND.(i4.eq.k1))))) cycle
            do j=2,divergence_bo(di_lst(i,1),1)+1
                do k=2,divergence_bo(di_lst(i,2),1)+1
                    do l=2,divergence_bo(di_lst(i,3),1)+1
                        do m=2,divergence_bo(di_lst(i,4),1)+1
                            added_di=added_di+1
                            add_di(added_di,1)=divergence_bo(di_lst(i,1),j)
                            add_di(added_di,2)=divergence_bo(di_lst(i,2),k)
                            add_di(added_di,3)=divergence_bo(di_lst(i,3),l)
                            add_di(added_di,4)=divergence_bo(di_lst(i,4),m)
                            add_di(added_di,5)=di_lst(i,5)
                        end do 
                    end do    
                end do 
            end do 
        end do 
    end do 
end do 


! do k=1,di_lstsz
!    print *,"bonded_type_torsion",di_lst(k,:)
! end do 
do i=1,di_lstsz
    found=.false.
    do j=1,reused
        if ((di_lst(i,1).eq.indexes_reused(j,2)).or.(di_lst(i,2).eq.indexes_reused(j,2)).or.&
        &(di_lst(i,3).eq.indexes_reused(j,2)).or.(di_lst(i,4).eq.indexes_reused(j,2))) then 
            found=.true.
            exit
        end if  
    end do 
    if (found.eqv..false.) then 
        print *,"bonded_type_torsion",di_lst(i,:)
     else 
         di_lst(i,:)=0

    end if 
end do 

 
 
 
 do i=1,added_di
     !
     i1=add_di(i,1)
     i2=add_di(i,2)
     i3=add_di(i,3)
     i4=add_di(i,4)
     found=.false.
     do k=1,di_lstsz
         if ((((i1.eq.di_lst(k,1)).and.(i2.eq.di_lst(k,2)).and.(i3.eq.di_lst(k,3)).and.(i4.eq.di_lst(k,4))).OR.&
          &((i1.eq.di_lst(k,4)).and.(i2.eq.di_lst(k,3)).and.(i3.eq.di_lst(k,2)).and.(i4.eq.di_lst(k,1))))) then 
            found=.true.
            exit
          end if 
     end do 
     
     if (found.eqv..true.) cycle 
     
     
     found=.false.
     do j=1,i-1
         k1=add_di(j,1)
         k2=add_di(j,2)
         k3=add_di(j,3)
         k4=add_di(j,4)
         
          if ((((i1.eq.k1).and.(i2.eq.k2).and.(i3.eq.k3).and.(i4.eq.k4)).OR.&
          &((i1.eq.k4).and.(i2.eq.k3).and.(i3.eq.k2).and.(i4.eq.k1)))) then 
            found=.true.
            exit
          end if 
     end do
     if (found.eqv..false.) then 
        print *,"bonded_type_torsion",add_di(i,:)
     end if 
 end do 
 
 

deallocate(add_di)

added_di=0
!Finally, I need to report the changes in biotypes for all other pairs not used in topology, but only useed for potentials
allocate(add_di(10000000,5))
add_di(:,:)=0
do rs=1,nseq
    do a=1,nrsimpt(rs)
        i1 = temp(b_type(iaa(rs)%impt(a,1)))
        i2 = temp(b_type(iaa(rs)%impt(a,2)))
        i3 = temp(b_type(iaa(rs)%impt(a,3)))
        i4 = temp(b_type(iaa(rs)%impt(a,4)))
        do i=1,impt_lstsz
            k1 = impt_lst(i,1)
            k2 = impt_lst(i,2)
            k3 = impt_lst(i,3)
            k4 = impt_lst(i,4)
            if (.NOT.((((i1.eq.k1).AND.(i2.eq.k2)).AND.&
 &                 ((i3.eq.k3).AND.(i4.eq.k4))).OR.&
 &                (((i1.eq.k1).AND.(i2.eq.k3)).AND.&
 &                 ((i3.eq.k2).AND.(i4.eq.k4))))) cycle

            do j=2,divergence_bo(impt_lst(i,1),1)+1
                do k=2,divergence_bo(impt_lst(i,2),1)+1
                    do l=2,divergence_bo(impt_lst(i,3),1)+1
                        do m=2,divergence_bo(impt_lst(i,4),1)+1
                            added_di=added_di+1
                            add_di(added_di,1)=divergence_bo(impt_lst(i,1),j)
                            add_di(added_di,2)=divergence_bo(impt_lst(i,2),k)
                            add_di(added_di,3)=divergence_bo(impt_lst(i,3),l)
                            add_di(added_di,4)=divergence_bo(impt_lst(i,4),m)
                            add_di(added_di,5)=impt_lst(i,5)
                        end do 
                    end do      
                end do 
            end do 
        end do 
    end do 
end do 

! do k=1,impt_lstsz
!    print *,"bonded_type_imptors",impt_lst(k,:)
! end do 
! 



do i=1,impt_lstsz
    found=.false.
    do j=1,reused
        if ((impt_lst(i,1).eq.indexes_reused(j,2)).or.(impt_lst(i,2).eq.indexes_reused(j,2)).or.&
        &(impt_lst(i,3).eq.indexes_reused(j,2)).or.(impt_lst(i,4).eq.indexes_reused(j,2))) then 
            found=.true.
            exit
        end if  
    end do 
    if (found.eqv..false.) then 
        print *,"bonded_type_imptors",impt_lst(i,:)
    else 
        impt_lst(i,:)=0
    end if 
end do 



 do i=1,added_di
     i1=add_di(i,1)
     i2=add_di(i,2)
     i3=add_di(i,3)
     i4=add_di(i,4)
     
     found=.false.
     do k=1,impt_lstsz
         if ((((i1.eq. impt_lst(k,1)).AND.(i2.eq. impt_lst(k,2))).AND.&
 &                 ((i3.eq. impt_lst(k,3)).AND.(i4.eq. impt_lst(k,4)))).OR.&
 &                (((i1.eq. impt_lst(k,1)).AND.(i2.eq.impt_lst(k,3))).AND.&
 &                 ((i3.eq.impt_lst(k,2)).AND.(i4.eq. impt_lst(k,4))))) then 
             found=.true.
          end if 
     end do
     if (found.eqv..true.) cycle 
     do j=1,i-1
         k1=add_di(j,1)
         k2=add_di(j,2)
         k3=add_di(j,3)
         k4=add_di(j,4)
         
            if ((((i1.eq.k1).AND.(i2.eq.k2)).AND.&
 &                 ((i3.eq.k3).AND.(i4.eq.k4))).OR.&
 &                (((i1.eq.k1).AND.(i2.eq.k3)).AND.&
 &                 ((i3.eq.k2).AND.(i4.eq.k4)))) then 

            found=.true.
            exit
          end if 
     end do
     if (found.eqv..false.) then 
        print *,"bonded_type_imptors",add_di(i,:)
     end if 
 end do 
 

do i=1,added_temp
    
    if (first_bond(i).ne.0) then 
        print *,"bond",int(mem_added_bond(i,1)),1,mem_added_bond(i,2:3)
    end if 
end do 

do i=1,added_temp
    if (first_angle(i).ne.0) then 
        print *,"angle",int(mem_added_angle(i,1)),1,mem_added_angle(i,2:3)
    end if 
end do 

! Martin : This is the special dummy atom case
    if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then 
      do i=1,n
          if ((blen_limits(i,2).eq.0.).and.(blen_limits(i,1).ne.0.)) then 
            blen_limits(i,2)=blen_limits(i,1)
          else if ((blen_limits(i,1).eq.0.).and.(blen_limits(i,2).ne.0.)) then 
            blen_limits(i,1)=blen_limits(i,2)
          end if 
          if ((blen_limits(i,3).eq.0.).and.(blen_limits(i,1).ne.0.)) then 
            blen_limits(i,3)=blen_limits(i,1)
          end if 
      end do 
    end if 
    ! Martin : scale the values now 
  if ((do_pka.eqv..true.).or.(do_hsq.eqv..true.)) then 
      blen(:)= blen_limits(:,2)*par_pka(2) + blen_limits(:,1)*(1.0-par_pka(2))
      bang(:) = bang_limits(:,2)*par_pka(2) + bang_limits(:,1)*(1.0-par_pka(2))
  end if 

  if ((added_bo.ne.old).or.(need_fix)) then 
    print *,old,added_bo,need_fix
    print *,"Please fix the parameter file before continuing"
    call fexit()
  end if 
  
end  
subroutine  replace_bond_values2_save() ! The goal of this routine is to replace the value of the bonds angle, dihedrals and length
    ! from the ones hard wired into campari to the ones present in the parameter files
    ! I could also try to make that for the case where the hardwired parameters are used,.
    !but their is the issue  of hardcoding each values, which is bad practice as well as confusing and lots of work to implement 
    ! I do not think it is necessary at all, so I will just keep it as one of the things to do later. 
  use atoms
  use iounit
  use params
  ! Martin ; needed for sure :
  use zmatrix
  use inter
  use sequen
  use martin_own, only: test_combi
  use energies
  use movesets ! For hsq
  ! Martin : check there is probably some unuseful extras
  integer i,j,k,l,rs,h,m,r,last,skipped,reused,bo_curr,old
  integer added_ba,added_bo,mem_ba(1000,3),mem_bo(1000,3),changed_ba,changed_bo,added_temp,temp_new,temp_new2
  
  integer, allocatable  :: mem_all(:,:)
  
  integer, allocatable :: array1(:),array2(:)
  integer, allocatable :: used_ba(:),used_bo(:)
  integer, allocatable :: memory(:),mem_b_type_used(:,:)
  RTYPE, allocatable :: angle_target(:),bond_target(:),angle_current(:),bond_current(:)
  RTYPE, allocatable :: mem_added_angle(:,:),mem_added_bond(:,:)
  integer, allocatable :: mem_added_bo_type(:,:),mem_changed_b_type(:,:),angle_index(:),bond_index(:)
  integer, allocatable :: mem_btype_switch(:,:),used_b_type(:),bo_type_used(:),indexes_reused(:,:)
  logical found
  integer ADDED_ANGLES,ADDED_BONDS,skip
  character my_format*100
  integer, allocatable:: b_types_used(:),bo_type_switch(:,:)
  
  integer, allocatable:: first_bond(:),first_angle(:)

  RTYPE, allocatable :: test_temp_array1(:),test_temp_array2(:)
  RTYPE mem_new_ba(1000),mem_new_bo(1000),temp_angle,temp_bond
  logical test,need_fix
  
  test=.true.
  need_fix=.false.
  found=.false.
  allocate(memory(3))
  allocate(array1(3))
  allocate(array2(3))
  allocate(test_temp_array1(size(blen(:))))
  allocate(test_temp_array2(size(bang(:))))
  
  allocate(mem_added_bo_type(100000,5))
  allocate(mem_added_angle(100000,3))
  allocate(mem_added_bond(100000,3))
  allocate(mem_changed_b_type(100000,3))

  allocate(mem_btype_switch(n_biotyp,3))

  mem_btype_switch(:,:)=0
  
  test_temp_array1(:)=blen(:)
  test_temp_array2(:)=bang(:)
  
  added_ba=0
  added_bo=0
  
  changed_ba=0
  changed_bo=0
  
  mem_ba(:,:)=0
  mem_bo(:,:)=0
  
  mem_new_ba(:)=0
  mem_new_bo(:)=0
  
  
  allocate(mem_all(n,3))
  allocate(angle_target(n))
  allocate(angle_index(n))
  allocate(bond_index(n))
  allocate(mem_b_type_used(n,3))

  mem_all(:,:)=0
  angle_target(:)=0.
  angle_index(:)=0
  bond_index(:)=0
  

  
  skipped=0
  ! So the revised goal is : 
  ! I identify biotype that need a change in 

  do i=1,n
    if ((iz(1,i).gt.0).and.(iz(2,i).gt.0)) then 
        array1(1)=bio_botyp(b_type(i))
        array1(2)=bio_botyp(b_type(iz(1,i)))
        array1(3)=bio_botyp(b_type(iz(2,i)))
    else
        cycle 
    end if 

    do k=1,ba_lstsz 
        array2(1)=ba_lst(k,1)
        array2(2)=ba_lst(k,2)
        array2(3)=ba_lst(k,3)

        test=.true.
        memory(:)=test_combi(array1,array2,size(memory(:)))

        do j=1,size(memory(:))
            if (memory(j).eq.0) then
                test=.false.
            end if 
        end do 
        if (test.eqv..true.) then 

            if (array2(2).ne.array1(2)) then 
                test=.false.
            end if 
        end if

        
        if (test.eqv..true.) then

            angle_index(i)=ba_lst(k,4)
            exit
        end if 
    end do 
  end do 
  deallocate(memory)
  deallocate(array1)
  deallocate(array2)

  allocate(memory(2))
  allocate(array1(2))
  allocate(array2(2))

! Second part, bond length
  do i=1,n
    if (iz(1,i).gt.0) then 
        array1(1)=bio_botyp(b_type(i))
        array1(2)=bio_botyp(b_type(iz(1,i)))
    else
        cycle 
    end if 

    do k=1,bo_lstsz      
        array2(:)=bo_lst(k,1:2)
        test=.true.
        memory(:)=test_combi(array1,array2,size(memory(:)))
        do j=1,size(memory(:))
            if (memory(j).eq.0) then
                test=.false.
            end if 
        end do 
        if (test.eqv..true.) then
            bond_index(i)=bo_lst(k,3)          
           exit 
        end if 
    end do    
  end do 
! Bond parameter check
allocate(used_b_type(n_biotyp))

used_b_type(:)=0
! This is a 
old=0
added_bo=0
added_angles=0
added_bonds=0
mem_added_bo_type(:,:)=0 !@ CHanged : ger the two connected atoms, then the bond type, then angle type

!mem_added_ba_type=(:,:)=0! Added bond_type_angle definition
mem_added_angle(:,:)=0. ! Added corresponding angle
mem_added_bond(:,:)=0. !added corresponding bond
mem_changed_b_type(:,:)=0 ! Memory of the modification to be made to biotype

deallocate(memory)
deallocate(array1)
deallocate(array2)

allocate(memory(3))
allocate(array1(3))
allocate(array2(3))

allocate(first_bond(n))
allocate(first_angle(n))
first_bond(:)=0
first_angle(:)=0


!Note that I am ALWAYS going to add a bo_type if the angle is incorrect
do i=1,n
    if ((bond_index(i).eq.0)) then ! .or.(angle_index(i).eq.0)) then 
        temp_bond =0
        
    else 
        temp_bond =bo_par(bond_index(i),2)
    end if 
!    if (mem_btype_switch(b_type(i),1).ne.0) then 
!        cycle  
!    end if    
    bo_curr=0
    found=.False.
    do k=1,added_bo
        if (mem_changed_b_type(k,1).eq.b_type(i)) then 
            found=.true.
            bo_curr=mem_changed_b_type(k,2)
            exit
        end if 
    end do 

!    if (found) then
!        cycle
!    end if 
    ! This is to avoid unallocated references
    if (angle_index(i).eq.0) then 
        temp_angle=0.
    else 
        temp_angle=ba_par(angle_index(i),2)
    end if 
    
        ! That means something needs to be changed
    if ((blen(i).ne.temp_bond).or.((bang(i).ne.temp_angle))) then !.and.((angle_index(i).ne.0))))then 
    
        array1(1)=bio_botyp(b_type(i))
        array1(2)=bio_botyp(b_type(iz(1,i)))
!        if (iz(2,i).ne.0) then 
!            array1(3)=bio_botyp(b_type(iz(2,i)))
!            print *,i,iz(1,i),iz(2,i)
!            print *,b_type(i),b_type(iz(1,i)),b_type(iz(2,i))
!        else 
!            array1(3)=0
!            print *,i,iz(1,i),iz(2,i)
!            print *,b_type(i),b_type(iz(1,i)),0
!        end if 
!        
!        print *,array1(:)
!        print *,blen(i),temp_bond,bang(i),temp_angle
        

        ! Now check if any of the newly created types has the right values.
        if (added_bo.ne.0) then 
            found=.false.
            do j=1,added_bo 

                    
                ! If angles and bonds are the same, reuse
                    if (((bang(i).eq.mem_added_angle(j,3)).and.(blen(i).eq.mem_added_bond(j,3)))) then                         
                        ! If not already changed

                            
                            added_bo=added_bo+1
                            mem_added_bo_type(added_bo,:)=mem_added_bo_type(j,:)
                            
                            mem_added_angle(added_bo,:)=mem_added_angle(j,:)
                            mem_added_bond(added_bo,:)=mem_added_bond(j,:)
                            
                            
                            mem_btype_switch(b_type(i),1)=bio_botyp(b_type(i))
                            mem_btype_switch(b_type(i),2)=mem_changed_b_type(j,2)                   
                            
                            mem_changed_b_type(added_bo,1)=b_type(i)
                            mem_changed_b_type(added_bo,2)=mem_changed_b_type(j,2)
                            skipped=skipped+1
                            
                            mem_b_type_used(added_bo,1)=b_type(i)
                            
                            mem_b_type_used(added_bo,2)=b_type(iz(1,i))
                            if (iz(2,i).ne.0) then 
                                mem_b_type_used(added_bo,3)=b_type(iz(2,i))
                            end if 
                            found=.true.
                            exit

                        
                    end if 
!                end if 
            end do
            !if you haven't found a match, create one 
            if (found.eqv..false.) then 
                if (bo_curr.eq.0) then 
                    added_bo=added_bo+1

                    mem_added_bo_type(added_bo,1)=n_botyp+added_bo-skipped
                else 
                    added_bo=added_bo+1

                    mem_added_bo_type(added_bo,1)=bo_curr
                    skipped=skipped+1
                end if 

                mem_added_bo_type(added_bo,2)=bio_botyp(b_type(iz(1,i)))
                mem_added_bo_type(added_bo,3)=bio_botyp(b_type(iz(2,i)))
                
                
                mem_b_type_used(added_bo,1)=b_type(i)
                mem_b_type_used(added_bo,2)=b_type(iz(1,i))
                if (iz(2,i).ne.0) then 
                    mem_b_type_used(added_bo,3)=b_type(iz(2,i))
                end if 
                found=.false.
                ! Now check if we can refer to an already existing angle
                do k=1,n_angltyp
                    if (ba_par(k,2).eq.bang(i)) then 
                                                    
                        mem_added_angle(added_bo,1)=k
                        mem_added_angle(added_bo,2)=100.0
                        mem_added_angle(added_bo,3)=ba_par(k,2)
                        !mem_added_angle(added_angles,1)=k
                        found=.true.
                        exit
                    end if 
                end do 
                
                ! That is one of the bond values corresponds to the desired value
                if (found) then 

                    mem_added_bo_type(added_bo,5)=k
                else ! if none
                !    print *,"DD"                    
                    do k=1,added_bo
                        if ((mem_added_angle(k,3).eq.bang(i)).and.(bang(i).ne.0.)) then 

                            mem_added_bo_type(added_bo,5)=mem_added_bo_type(k,5)
                            
                            mem_added_angle(added_bo,1)=mem_added_angle(k,1)
                            mem_added_angle(added_bo,2)=mem_added_angle(k,2)
                            mem_added_angle(added_bo,3)=mem_added_angle(k,3)
                            exit 
                        end if 
                    end do 
                    
                    if ((k.eq.added_bo+1).and.(bang(i).ne.0.)) then 

                        added_angles=added_angles+1
                        mem_added_bo_type(added_bo,5)=n_angltyp+added_angles
                        first_angle(added_bo)=1
                        mem_added_angle(added_bo,1)=n_angltyp+added_angles
                        mem_added_angle(added_bo,2)=100.00 ! The spring constant, pretty irelevant for our purpose
                        mem_added_angle(added_bo,3)=bang(i)
                    end if 
                end if 
                
                ! Now for the bonds
                ! Now check if we can refer to an already existing angle
                do k=1,n_bondtyp 
                    if (bo_par(k,2).eq.blen(i)) then 
                        mem_added_bond(added_bo,1)=k
                        mem_added_bond(added_bo,2)=100.0
                        mem_added_bond(added_bo,3)=bo_par(k,2)

                        exit
                    end if 
                end do 
                
                ! That is one of the bond values corresponds to the desired value
                if (k.ne.n_bondtyp+1) then 

                    mem_added_bo_type(added_bo,4)=k
                else ! if none
                    
                    do k=1,added_bo
                        if (mem_added_bond(k,3).eq.blen(i)) then 
                            mem_added_bond(added_bo,1)=mem_added_bond(k,1)
                            mem_added_bond(added_bo,2)=mem_added_bond(k,2)
                            mem_added_bond(added_bo,3)=mem_added_bond(k,3)
                            mem_added_bo_type(added_bo,4)=mem_added_bo_type(k,4)
                            exit 
                        end if 
                    end do 
                    
                    if (k.eq.added_bo+1) then 

                        added_bonds=added_bonds+1
                        first_bond(added_bo)=1
                        mem_added_bo_type(added_bo,4)=n_bondtyp+added_bonds
                        
                        mem_added_bond(added_bo,1)=n_bondtyp+added_bonds
                        mem_added_bond(added_bo,2)=50 ! The spring constant, pretty irelevant for our purpose
                        mem_added_bond(added_bo,3)=blen(i)
                    end if 

                end if 
                mem_btype_switch(b_type(i),1)=bio_botyp(b_type(i))
                mem_btype_switch(b_type(i),2)=n_botyp+added_bo-skipped
                
                mem_changed_b_type(added_bo,1)=b_type(i)
                mem_changed_b_type(added_bo,2)=n_botyp+added_bo-skipped

            end if 
        else ! This is to make sure we do not add anything that unecessary

            ! That is only for the first element added
            added_bo=added_bo+1

            mem_added_bo_type(added_bo,1)=n_botyp+added_bo-skipped
            mem_added_bo_type(added_bo,2)=bio_botyp(b_type(iz(1,i)))
            
            mem_b_type_used(added_bo,1)=b_type(i)
            mem_b_type_used(added_bo,2)=b_type(iz(1,i))
            if (iz(2,i).ne.0) then 
            mem_b_type_used(added_bo,3)=b_type(iz(2,i))
            end if 
            if (temp_angle.ne.0.) then 
                mem_added_bo_type(added_bo,3)=bio_botyp(b_type(iz(2,i)))
            else 
                mem_added_bo_type(added_bo,3)=0
            end if 
            found=.false.
            ! Now check if we can refer to an already existing angle
            do k=1,n_angltyp
                if ((ba_par(k,2).eq.bang(i)).and.(bang(i).ne.0.)) then 

                        mem_added_angle(added_bo,1)=k
                        mem_added_angle(added_bo,2)=100.0
                        mem_added_angle(added_bo,3)=ba_par(k,2)
                    found=.true.
                    exit
                end if 
            end do 
            ! That is one of the bond values corresponds to the desired value
            if (found) then 


                mem_added_bo_type(added_bo,5)=k
            else if (bang(i).ne.0.) then ! if none
                added_angles=added_angles+1
                mem_added_bo_type(added_bo,5)=n_angltyp+added_angles
                
                mem_added_angle(added_bo,1)=n_angltyp+added_angles
                first_angle(added_bo)=1
                if (temp_angle.ne.0.) then 
                    mem_added_angle(added_bo,2)=100.
                else 
                    mem_added_angle(added_bo,2)=0.
                end if 
                mem_added_angle(added_bo,3)=bang(i)
            end if 

            ! Now for the bonds
            ! Now check if we can refer to an already existing angle
            do k=1,n_bondtyp
                if (bo_par(k,2).eq.blen(i)) then 
                    mem_added_bond(added_bo,1)=k
                    mem_added_bond(added_bo,2)=100.0
                    mem_added_bond(added_bo,3)=bo_par(k,2)

                    exit
                end if 
            end do 
            ! That is one of the bond values corresponds to the desired value
            if (k.ne.n_bondtyp+1) then 


                mem_added_bo_type(added_bo,4)=k
            else ! if none


                added_bonds=added_bonds+1
                first_bond(added_bo)=1

                mem_added_bo_type(added_bo,4)=n_bondtyp+added_bonds
                mem_added_bond(added_bo,1)=n_bondtyp+added_bonds
                mem_added_bond(added_bo,2)=50.
                mem_added_bond(added_bo,3)=blen(i)
            end if 
            mem_btype_switch(b_type(i),1)=bio_botyp(b_type(i))
            mem_btype_switch(b_type(i),2)=n_botyp+added_bo-skipped
            
            mem_changed_b_type(added_bo,1)=b_type(i)
            mem_changed_b_type(added_bo,2)=n_botyp+added_bo-skipped
            
        end if   

    else if (iz(2,i).ne.0) then  
        old=old+1
        added_bo=added_bo+1
        mem_added_bo_type(added_bo,1)=bio_botyp(b_type(i))
        mem_added_bo_type(added_bo,2)=bio_botyp(b_type(iz(1,i)))

        mem_b_type_used(added_bo,1)=b_type(i)
        mem_b_type_used(added_bo,2)=b_type(iz(1,i))
        if (iz(2,i).ne.0) then 
            mem_b_type_used(added_bo,3)=b_type(iz(2,i))
            mem_added_bo_type(added_bo,3)=bio_botyp(b_type(iz(2,i)))
        else 
            mem_b_type_used(added_bo,3)=0
            mem_added_bo_type(added_bo,3)=0
        end if 
        mem_added_bo_type(added_bo,4)=bond_index(i)
        mem_added_bo_type(added_bo,5)=angle_index(i)
        
        
        array1(1)=bio_botyp(b_type(i))
        array1(2)=bio_botyp(b_type(iz(1,i)))

    else 
        
        added_bo=added_bo+1
        mem_btype_switch(b_type(i),1)=bio_botyp(b_type(i))
        mem_btype_switch(b_type(i),2)=n_botyp+added_bo-skipped
        

    end if  
    ! At that point we already have changed what needs change
    used_b_type(b_type(i))=1
end do


allocate(bo_type_used(n_botyp+added_bo))
bo_type_used(:)=0
last=0
do k=1,n_biotyp
    do i=1,added_bo
        !This is just to be in the right order
        if (mem_changed_b_type(i,1).eq.k) then
            last=last+1

            bio_botyp(k)=mem_changed_b_type(i,2)

            if (mem_changed_b_type(i,2).ne.0) then 
                bo_type_used(mem_changed_b_type(i,2)-n_botyp)=1
            end if
        end if 
    end do 
    
    if (last.ne.k) then 
        last=last+1
        

        if  (bio_botyp(k).ne.0) then 
            bo_type_used(bio_botyp(k))=1
        end if 
    end if
end do 
444   format(A8,I4,A,A,A,A,A,I3,A,I3,A,I3)
    
  allocate(b_types_used(n_botyp+added_bo+100))
  b_types_used(:)=0
  do i=1,n_biotyp
      if (bio_botyp(i).ne.0) then 
            b_types_used(bio_botyp(i))=1
      end if 
  end do 
  
  

 ! I should try reorganizing the bonded types, to avoid making dummies

  allocate(indexes_reused(added_bo,2))
  indexes_reused(:,:)=0
  reused=0
  skip=0
  last =n_biotyp
   445   format(A8,I4,A,I3,A,I3,A,I3)
   do i=1,n_botyp+added_bo
       !print *,i,b_types_used(i)
       if (b_types_used(i).eq.0) then
           ! This
           if (i.gt.n_botyp+added_bo-reused-skip+1) then 
               cycle 
           end if 
           reused=reused+1

           do while (b_types_used(n_botyp+added_bo-reused-skip+1).eq.0)

               skip=skip+1
           end do
           if (n_botyp+added_bo-reused-skip+1.gt.i) then 
                indexes_reused(reused,1)=n_botyp+added_bo-reused-skip+1
                
                indexes_reused(reused,2)=i
                
           end if 

       end if 
   end do 

   do i=1,added_bo
       do j=1,reused 
           if (indexes_reused(j,1).eq.mem_added_bo_type(i,1)) then 
               !print *,"unused1",indexes_reused(j,:)
               mem_added_bo_type(i,1)=indexes_reused(j,2)
           end if
           if (indexes_reused(j,1).eq.mem_added_bo_type(i,2)) then 
               !print *,"unused2",indexes_reused(j,:)
               mem_added_bo_type(i,2)=indexes_reused(j,2)
           end if
           if (indexes_reused(j,1).eq.mem_added_bo_type(i,3)) then 
               !print *,"unused3",indexes_reused(j,:)
               mem_added_bo_type(i,3)=indexes_reused(j,2)
           end if
       end do 
   end do 

   do i=1,added_bo

       do j=1,reused 

           if (indexes_reused(j,1).eq.mem_changed_b_type(i,2)) then 
               mem_changed_b_type(i,2)=indexes_reused(j,2)
           end if
       end do 
   end do 
   
   do i=1,n_biotyp
       do j=1,reused 
           if (indexes_reused(j,1).eq.mem_btype_switch(i,2)) then 
               print *,i
               print *,indexes_reused(j,:)
               
               mem_btype_switch(i,2)=indexes_reused(j,2)
           end if
       end do 
   end do 
   print *,"Start of reused"
   do j=1,reused 
    print *,indexes_reused(j,:)
   end do 

   
   
bo_type_used(:)=0
last=0

do k=1,n_biotyp
    if (mem_btype_switch(k,2).ne.0) then 
            write(ilog, 444 ) "biotype ",k,"  ",bio_code(k),'  "',bio_res(k),'"  ',&
            &bio_ljtyp(k),"  ",bio_ctyp(k),"  ",mem_btype_switch(k,2)

    else 
        write(ilog, 444 ) "biotype ",k,"  ",bio_code(k),'  "',bio_res(k),'"  ',&
            &bio_ljtyp(k),"  ",bio_ctyp(k),"  ",bio_botyp(k)
            
            
    end if
end do 

! Now I need to change the angle definition for each types that were charged
! For angles, we need a double loop though biotypes 
added_temp=added_bo

deallocate(memory)
deallocate(array1)
deallocate(array2)

allocate(memory(3))
allocate(array1(3))
allocate(array2(3))

n_biotyp= last 

! Here just addding 

do i=1,added_bo
    do j=1,3
        if (mem_b_type_used(i,j).eq.0) cycle
        if (mem_btype_switch(mem_b_type_used(i,j),1).ne.0) then 
            mem_added_bo_type(i,j)=mem_btype_switch(mem_b_type_used(i,j),2)

        end if
    end do 
end do 

!Here I need to duplicate 


deallocate(array1)
deallocate(array2)
deallocate(memory)

allocate(array1(2))
allocate(array2(2))
allocate(memory(2))


do i=1,bo_lstsz
    do j=1,added_bo
        found=.false.
        array1(:)=bo_lst(i,1:2)
        array2(:)=mem_added_bo_type(j,1:2)
        memory(:)=test_combi(array1,array2,size(memory(:)))
        if ((memory(1).ne.0).and.(memory(2).ne.0)) then
          
            found=.true.
            exit
        end if 
    end do 
   
    if (found.eqv..false.) then 
        print *,"bonded_type_bond",bo_lst(i,:)
    end if 
end do 
deallocate(array1)
deallocate(array2)
deallocate(memory)

allocate(array1(2))
allocate(array2(2))
allocate(memory(2))

do i=1,added_bo
    found=.false.
    ! This is to not double print a similar bond definition
    if (mem_added_bo_type(i,4).eq.0) then
        cycle 
    end if
    
    do j=1,i-1
        array1=mem_added_bo_type(i,1:2)
        array2=mem_added_bo_type(j,1:2)
        
        memory(:)=test_combi(array1,array2,size(memory(:)))
        if ((memory(1).ne.0).and.(memory(2).ne.0)) then
            found=.true.
            exit
        end if 
    end do 
!    print *,i
!    print *,array1
!    print *,array2
    if (found.eqv..false.) then 
        
        print *,"bonded_type_bond",mem_added_bo_type(i,1:2),mem_added_bo_type(i,4)
    end if 
end do 

deallocate(array1)
deallocate(array2)
deallocate(memory)

allocate(array1(3))
allocate(array2(3))
allocate(memory(3))

do i=1,ba_lstsz
    found=.false.
    do j=1,added_bo
        array1(:)=ba_lst(i,1:3)
        array2(:)=mem_added_bo_type(j,1:3)
        memory(:)=test_combi(array1,array2,size(memory(:)))
        if ((memory(1).ne.0).and.(array2(2).eq.array1(2)).and.(memory(3).ne.0)) then
            
            found=.true.
            exit
        end if 
    end do 
    if (found.eqv..false.) then 
        print *,"bonded_type_angle",ba_lst(i,:)
    end if  
end do 

deallocate(array1)
deallocate(array2)
deallocate(memory)

allocate(array1(3))
allocate(array2(3))
allocate(memory(3))
do i=1,added_bo
    if (mem_added_bo_type(i,5).eq.0) then
        cycle 
    end if
    found=.false.
    do j=1,i-1
        array1=mem_added_bo_type(i,1:3) ! Changed to 1:3
        array2=mem_added_bo_type(j,1:3)
        
        memory(:)=test_combi(array1,array2,size(memory(:)))
!        if ((((array1(1).eq.295).and.(array1(2).eq.295).and.(array1(3).eq.309)).or.&
!        &((array1(1).eq.309).and.(array1(2).eq.295).and.(array1(3).eq.295))).and.((array2(1).eq.295).or.(array2(1).eq.309))&
!        &) then 
!            
!            print *,"mem"
!            PRINT *,array1
!            PRINT *,array2
!            PRINT *,memory
!            
!        end if 
        if (memory(2).ne.0) then 
            if ((memory(1).ne.0).and.(array2(2).eq.array1(2)).and.(memory(3).ne.0)) then
                found=.true.
                exit
            end if 
        end if 
    end do 
    
    if (found.eqv..false.) then 
        print *,"bonded_type_angle",mem_added_bo_type(i,1:3),mem_added_bo_type(i,5)
    end if 
end do


do i=1,added_bo
    if (first_bond(i).ne.0) then 
        print *,"bond",int(mem_added_bond(i,1)),1,mem_added_bond(i,2:3)
    end if 
end do 

do i=1,added_bo
    if (first_angle(i).ne.0) then 
        print *,"angle",int(mem_added_angle(i,1)),1,mem_added_angle(i,2:3)
    end if 
end do 




! Martin : This is the special dummy atom case
    if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then 
      do i=1,n
          if ((blen_limits(i,2).eq.0.).and.(blen_limits(i,1).ne.0.)) then 
            blen_limits(i,2)=blen_limits(i,1)
          else if ((blen_limits(i,1).eq.0.).and.(blen_limits(i,2).ne.0.)) then 
            blen_limits(i,1)=blen_limits(i,2)
          end if 
          if ((blen_limits(i,3).eq.0.).and.(blen_limits(i,1).ne.0.)) then 
            blen_limits(i,3)=blen_limits(i,1)
          end if 
      end do 
    end if 
    ! Martin : scale the values now 
  if ((do_pka.eqv..true.).or.(do_hsq.eqv..true.)) then 
      blen(:)= blen_limits(:,2)*par_pka(2) + blen_limits(:,1)*(1.0-par_pka(2))
      bang(:) = bang_limits(:,2)*par_pka(2) + bang_limits(:,1)*(1.0-par_pka(2))
  end if 

  if ((added_bo.ne.old).or.(need_fix)) then 
    print *,old,added_bo,need_fix
    print *,"Please fix the parameter file before continuing"
    call fexit()
  end if 
  
end  
! This is now unused, 
subroutine  replace_bond_values() ! The goal of this routine is to replace the value of the bonds angle, dihedrals and length
    ! from the ones hard wired into campari to the ones present in the parameter files
    ! I could also try to make that for the case where the hardwired parameters are used,.
    !but their is the issue  of hardcoding each values, which is bad practice as well as confusing and lots of work to implement 
    ! I do not think it is necessary at all, so I will just keep it as one of the things to do later. 
  use atoms
  use iounit
  use params
  ! Martin ; needed for sure :
  use zmatrix
  use inter
  use sequen
  use martin_own, only: test_combi
  use energies
  use movesets ! For hsq
  ! Martin : check there is probably some unuseful extras
  integer i,j,k,l,rs,h
  integer added_ba,added_bo,mem_ba(1000,3),mem_bo(10000,3),changed_ba,changed_bo
  
  integer, allocatable  :: mem_all(:,:)
  
  integer, allocatable :: array1(:),array2(:)
  integer, allocatable :: used_ba(:),used_bo(:)
  integer, allocatable :: memory(:)
  RTYPE, allocatable:: angle_target(:),bond_target(:),angle_current(:),bond_current(:),angle_index(:),bond_index(:)
  logical found
  integer i1,i2,i3,k1,k2,k3
  
  
  
  RTYPE, allocatable :: test_temp_array1(:),test_temp_array2(:)
  RTYPE mem_new_ba(10000),mem_new_bo(10000)
  logical test,need_fix
  test=.true.
  need_fix=.false.
  found=.false.
  allocate(memory(3))
  allocate(array1(3))
  allocate(array2(3))
 
   do i=1,n
    if ((iz(1,i).gt.0).and.(iz(2,i).gt.0)) then 
        i1=bio_botyp(b_type(i))
        i2=bio_botyp(b_type(iz(1,i)))
        i3=bio_botyp(b_type(iz(2,i)))
    else
        cycle 
    end if 
    found=.false.
    do k=1,ba_lstsz 
        
        k1=ba_lst(k,1)
        k2=ba_lst(k,2)
        k3=ba_lst(k,3)

        test=.true.
!        memory(:)=test_combi(array1,array2,size(memory(:)))
!
!        do j=1,size(memory(:))
!            if (memory(j).eq.0) then
!                test=.false.
!            end if 
!        end do 
!        
!        if (test.eqv..true.) then 
!            if (array2(memory(2)).ne.array1(2)) then 
!                test=.false.
!            end if 
!        end if
        
        if (.not.(((i1.eq.k1).and.(i2.eq.k2).and.(i3.eq.k3)).or.((i1.eq.k3).and.(i2.eq.k2).and.(i3.eq.k1)))) cycle
           
        

        if (bang(i).ne.ba_par(ba_lst(k,4),2)) then 

            if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then
                bang_limits(i,1)=ba_par(ba_lst(k,4),2)
            else 
                bang(i)=ba_par(ba_lst(k,4),2)
            end if 
        
        end if 
        exit 
        ! This is to avoid a lack of definition

    end do    
!    if (found.eqv..false.) then 
!        print *,"Could not find angle corresponding to "
!        print *,array1(:)
!    end if 

    if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then 
        if ((iz(1,i).gt.0).and.(iz(2,i).gt.0)) then 
            i1=bio_botyp(transform_table(b_type(i)))
            i2=bio_botyp(transform_table(b_type(iz(1,i))))
            i3=bio_botyp(transform_table(b_type(iz(2,i))))
        else 
            cycle
        end if 
        
        do k=1,ba_lstsz      
            k1=ba_lst(k,1)
            k2=ba_lst(k,2)
            k3=ba_lst(k,3)
            if (.not.(((i1.eq.k1).and.(i2.eq.k2).and.(i3.eq.k3)).or.((i1.eq.k3).and.(i2.eq.k2).and.(i3.eq.k1)))) cycle
            bang_limits(i,2)=ba_par(ba_lst(k,4),2)
            exit

        end do
        
        if (nhis.ne.0) then ! if there is a histidine (HSQ only)

            if ((iz(1,i).gt.0).and.(iz(2,i).gt.0)) then 
                i1=bio_botyp(transform_table(his_eqv_table(b_type(i))))
                i2=bio_botyp(transform_table(his_eqv_table(b_type(iz(1,i)))))
                i3=bio_botyp(transform_table(his_eqv_table(b_type(iz(2,i)))))
            else 
                cycle
            end if 

            do k=1,ba_lstsz      
                k1=ba_lst(k,1)
                k2=ba_lst(k,2)
                k3=ba_lst(k,3)
                if (.not.(((i1.eq.k1).and.(i2.eq.k2).and.(i3.eq.k3)).or.((i1.eq.k3).and.(i2.eq.k2).and.(i3.eq.k1)))) cycle

                bang_limits(i,3)=ba_par(ba_lst(k,4),2)
                exit
            end do    
        end if 
    end if 
  end do 
! Martin : This is the special dummy atom case
  if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then 
      do i=1,n
          if ((bang_limits(i,2).eq.0.).and.(bang_limits(i,1).ne.0.)) then 
              bang_limits(i,2)=bang_limits(i,1)
          else if ((bang_limits(i,1).eq.0.).and.(bang_limits(i,2).ne.0.)) then 
              bang_limits(i,1)=bang_limits(i,2)
          end if 
          if ((bang_limits(i,3).eq.0.).and.(bang_limits(i,1).ne.0.)) then 
              bang_limits(i,3)=bang_limits(i,1)
          end if 
      end do 
    end if 
   
 

  deallocate(memory)
  deallocate(array1)
  deallocate(array2)

  allocate(memory(2))
  allocate(array1(2))
  allocate(array2(2))

! Second part, bond length
  do i=1,n
    if (iz(1,i).gt.0) then 
        i1=bio_botyp(b_type(i))
        i2=bio_botyp(b_type(iz(1,i)))
    else
        cycle 
    end if 

    do k=1,bo_lstsz      
        k1=bo_lst(k,1)
        k2=bo_lst(k,2)
        if (.not.(((i1.eq.k1).and.(i2.eq.k2)).or.((i1.eq.k2).and.(i2.eq.k1)))) cycle

        if (blen(i).ne.bo_par(bo_lst(k,3),2)) then 
            if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then 
                blen_limits(i,1)=bo_par(bo_lst(k,3),2)
            else 
                blen(i)=bo_par(bo_lst(k,3),2)
            end if
        end if 
        exit 
            

    end do 

    if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then 
        if (iz(1,i).gt.0) then 
            i1=bio_botyp(transform_table(b_type(i)))
            i2=bio_botyp(transform_table(b_type(iz(1,i))))
        else
            cycle ! Martin : this is not necessary because it already does it in the previous loop
        end if 

        do k=1,bo_lstsz      
            k1=bo_lst(k,1)
            k2=bo_lst(k,2)
            
            
            if (.not.(((i1.eq.k1).and.(i2.eq.k2)).or.((i1.eq.k2).and.(i2.eq.k1)))) cycle
            

                
            blen_limits(i,2)=bo_par(bo_lst(k,3),2)

            exit
        end do  
        if (nhis.ne.0) then 
            if (iz(1,i).gt.0) then 
                i1=bio_botyp(transform_table(his_eqv_table(b_type(i))))
                i2=bio_botyp(transform_table(his_eqv_table(b_type(iz(1,i)))))
            else
                cycle ! Martin : this is not necessary because it already does it in the previous loop
            end if 

            do k=1,bo_lstsz      
                k1=bo_lst(k,1)
                k2=bo_lst(k,2)

                if (.not.(((i1.eq.k1).and.(i2.eq.k2)).or.((i1.eq.k2).and.(i2.eq.k1)))) cycle
                
                blen_limits(i,3)=bo_par(bo_lst(k,3),2)
                exit

            end do  
        end if 
    end if 
end do 
! Bond Length parameter check

! Martin : This is the special dummy atom case
    if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then 
      do i=1,n
          if ((blen_limits(i,2).eq.0.).and.(blen_limits(i,1).ne.0.)) then 
            blen_limits(i,2)=blen_limits(i,1)
          else if ((blen_limits(i,1).eq.0.).and.(blen_limits(i,2).ne.0.)) then 
            blen_limits(i,1)=blen_limits(i,2)
          end if 
          if ((blen_limits(i,3).eq.0.).and.(blen_limits(i,1).ne.0.)) then 
            blen_limits(i,3)=blen_limits(i,1)
          end if 
      end do 
    end if 
    ! Martin : scale the values now 
  if ((do_pka.eqv..true.).or.(do_hsq.eqv..true.)) then 
      blen(:)= blen_limits(:,2)*par_pka(2) + blen_limits(:,1)*(1.0-par_pka(2))
      bang(:) = bang_limits(:,2)*par_pka(2) + bang_limits(:,1)*(1.0-par_pka(2))
  end if 


end  
!subroutine  replace_bond_values3() ! The goal of this routine is to replace the value of the bonds angle, dihedrals and length
!    ! from the ones hard wired into campari to the ones present in the parameter files
!    ! I could also try to make that for the case where the hardwired parameters are used,.
!    !but their is the issue  of hardcoding each values, which is bad practice as well as confusing and lots of work to implement 
!    ! I do not think it is necessary at all, so I will just keep it as one of the things to do later. 
!  use atoms
!  use iounit
!  use params
!  ! Martin ; needed for sure :
!  use zmatrix
!  use inter
!  use sequen
!  use martin_own, only: test_combi
!  use energies
!  use movesets ! For hsq
!  ! Martin : check there is probably some unuseful extras
!  integer i,j,k,l,rs,h
!  integer added_ba,added_bo,mem_ba(1000,3),mem_bo(1000,3),changed_ba,changed_bo
!  
!  integer, allocatable  :: mem_all(:,:)
!  
!  integer, allocatable :: array1(:),array2(:)
!  integer, allocatable :: used_ba(:),used_bo(:)
!  integer, allocatable :: memory(:)
!  RTYPE, allocatable:: angle_target(:),bond_target(:),angle_current(:),bond_current(:),angle_index(:),bond_index(:)
!  logical found
!  
!  
!  RTYPE, allocatable :: test_temp_array1(:),test_temp_array2(:)
!  RTYPE mem_new_ba(1000),mem_new_bo(1000)
!  logical test,need_fix
!  test=.true.
!  need_fix=.false.
!  found=.false.
!  allocate(memory(3))
!  allocate(array1(3))
!  allocate(array2(3))
!  allocate(test_temp_array1(size(blen(:))))
!  allocate(test_temp_array2(size(bang(:))))
!  
!
!  allocate(used_ba(ba_lstsz))
!  allocate(used_bo(bo_lstsz))
!  
!  used_ba(:)=0
!  used_bo(:)=0
!  
!  test_temp_array1(:)=blen(:)
!  test_temp_array2(:)=bang(:)
!  
!  added_ba=0
!  added_bo=0
!  
!  changed_ba=0
!  changed_bo=0
!  
!  mem_ba(:,:)=0
!  mem_bo(:,:)=0
!  
!  mem_new_ba(:)=0
!  mem_new_bo(:)=0
!  
!  
!  allocate(mem_all(n,3))
!  allocate(angle_target(n))
!  allocate(angle_index(n))
!  allocate(bond_index(n))
!  ! So the revised goal is : 
!  ! I dentify biotype that need a change in 
!
!  
!  
!  
!  
!  
!  
!  do i=1,n
!    if ((iz(1,i).gt.0).and.(iz(2,i).gt.0)) then 
!        array1(1)=bio_botyp(b_type(i))
!        array1(2)=bio_botyp(b_type(iz(1,i)))
!        array1(3)=bio_botyp(b_type(iz(2,i)))
!    else
!        cycle 
!    end if 
!
!    do k=1,ba_lstsz 
!        array2(1)=ba_lst(k,1)
!        array2(2)=ba_lst(k,2)
!        array2(3)=ba_lst(k,3)
!
!        test=.true.
!        memory(:)=test_combi(array1,array2,size(memory(:)))
!
!        do j=1,size(memory(:))
!            if (memory(j).eq.0) then
!                test=.false.
!            end if 
!        end do 
!        if (test.eqv..true.) then 
!            if (array2(memory(2)).ne.array1(2)) then 
!                test=.false.
!            end if 
!        end if
!        
!        if (test.eqv..true.) then
!
!            
!            if (bang(i).ne.ba_par(ba_lst(k,4),2)) then 
!
!                changed_ba=changed_ba+1
!                
!                mem_ba(changed_ba,1)=i
!                mem_ba(changed_ba,2)=k
!                
!                ! There I am going to see if that required value already exist. If it does not, readily print the lines that need modifying
!                do h=1,n_angltyp
!                    if (bang(i).eq.ba_par(h,2)) then
!                        ! The first dimension is going to be the number of the 
!                        !write(*,*) ba_lst(k,1:3),h
!                        mem_ba(changed_ba,3)=h
!                        !need_fix=.true.
!                        exit
!                    end if 
!                end do 
!                
!                if (h-1.eq.n_angltyp) then 
!                    print *,"Could not find an already existing angle value."
!                    print *,"Please add the following line to the angle definition"
!                    write(*,*) "angle",n_angltyp+added_ba+1,ba_typ(ba_lst(k,4)),ba_par(ba_lst(k,4),1),bang(i)
!                    write(*,*) "bonded_type_angle",array2(:),n_angltyp+added_ba+1
!                    print *,"And replace bond angle ",ba_lst(k,4),"by",n_angltyp+added_ba+1
!                    added_ba=added_ba+1
!                end if 
!            else 
!                used_ba(k)=1
!            end if 
!            
!            if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then
!                bang_limits(i,1)=ba_par(ba_lst(k,4),2)
!            else 
!                bang(i)=ba_par(ba_lst(k,4),2)
!            end if 
!            exit 
!            
!        else if (k.eq.ba_lstsz) then 
!            ! If there is not a bond definition, create one.
!            !print *,"k",k,ba_lst(k,:)
!            do h=1,n_angltyp
!                if (bang(i).eq.ba_par(h,2)) then
!                    
!                    !print *,"missing bond definition."
!                    !print *,"please add the following line"
!                    print *,"bonded_type_angle",array1,h
!                    test=.true.
!                    exit
!                end if 
!            end do 
!            if (test.eqv..false.) then 
!                print *,"no angles found for this value",bang(i)
!            end if 
!        end if 
!    end do    
!
!
!    if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then 
!        if ((iz(1,i).gt.0).and.(iz(2,i).gt.0)) then 
!            array1(1)=bio_botyp(transform_table(b_type(i)))
!            array1(2)=bio_botyp(transform_table(b_type(iz(1,i))))
!            array1(3)=bio_botyp(transform_table(b_type(iz(2,i))))
!        else 
!            cycle
!        end if 
!        
!        do k=1,ba_lstsz      
!            array2(1)=ba_lst(k,1)
!            array2(2)=ba_lst(k,2)
!            array2(3)=ba_lst(k,3)
!
!            test=.true.
!            memory(:)=test_combi(array1,array2,size(memory(:)))
!
!            do j=1,size(memory(:))
!                if (memory(j).eq.0) then
!                    test=.false.
!                end if 
!            end do 
!            if (test.eqv..true.) then 
!                if (array2(memory(2)).ne.array1(2)) then 
!                    test=.false.
!                end if 
!            end if 
!            if (test.eqv..true.) then 
!                bang_limits(i,2)=ba_par(ba_lst(k,4),2)
!                exit
!            end if 
!        end do
!        
!        if (nhis.ne.0) then ! if there is a histidine (HSQ only)
!
!            if ((iz(1,i).gt.0).and.(iz(2,i).gt.0)) then 
!                array1(1)=bio_botyp(transform_table(his_eqv_table(b_type(i))))
!                array1(2)=bio_botyp(transform_table(his_eqv_table(b_type(iz(1,i)))))
!                array1(3)=bio_botyp(transform_table(his_eqv_table(b_type(iz(2,i)))))
!            else 
!                cycle
!            end if 
!
!            do k=1,ba_lstsz      
!                array2(1)=ba_lst(k,1)
!                array2(2)=ba_lst(k,2)
!                array2(3)=ba_lst(k,3)
!
!                test=.true.
!                memory(:)=test_combi(array1,array2,size(memory(:)))
!
!                do j=1,size(memory(:))
!                    if (memory(j).eq.0) then
!                        test=.false.
!                    end if 
!                end do 
!                if (test.eqv..true.) then 
!                    if (array2(memory(2)).ne.array1(2)) then 
!                        test=.false.
!                    end if 
!                end if
!                
!                if (test.eqv..true.) then 
!                    bang_limits(i,3)=ba_par(ba_lst(k,4),2)
!                    exit
!                end if 
!            end do    
!        end if 
!    end if 
!  end do 
!! Martin : This is the special dummy atom case
!  if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then 
!      do i=1,n
!          if ((bang_limits(i,2).eq.0.).and.(bang_limits(i,1).ne.0.)) then 
!              bang_limits(i,2)=bang_limits(i,1)
!          else if ((bang_limits(i,1).eq.0.).and.(bang_limits(i,2).ne.0.)) then 
!              bang_limits(i,1)=bang_limits(i,2)
!          end if 
!          if ((bang_limits(i,3).eq.0.).and.(bang_limits(i,1).ne.0.)) then 
!              bang_limits(i,3)=bang_limits(i,1)
!          end if 
!      end do 
!    end if 
!   
! do j=1,changed_ba 
!    do k=1, changed_ba
!        ! Here i need to identify how many values share this line
!        if ((mem_ba(j,3).eq.mem_ba(k,3)).and.(bang(mem_ba(j,1)).ne.bang(mem_ba(k,1)))) then 
!            found=.false.
!            do l=1,added_ba
!                if ((mem_new_ba(l).eq.bang(mem_ba(j,1))).or.(mem_new_ba(l).ne.bang(mem_ba(k,1)))) then 
!                    found=.true.
!                    exit
!                end if 
!            end do 
!            if (found.eqv..true.) then 
!                exit
!            end if 
!            ! If this one is used already, it needs two more (one for each new target values)
!            if (used_ba(mem_ba(j,2)).eq.1) then 
!                ! Change each of the bond types, and add one for each
!                write(*,*) "bonded_type_angle",array2(1:2),n_botyp+added_ba+1,n_angltyp+added_ba+1
!                write(*,*) "angle",n_angltyp+added_ba+1,ba_typ(ba_lst(mem_ba(j,2),4)),ba_par(ba_lst(mem_ba(j,2),4),1),&
!                &bang(mem_ba(j,1))
!
!                mem_new_ba(added_ba+1)=bang(mem_ba(j,1))
!                print *,"replace bond type to",n_botyp+added_ba+1,"in biotype",b_type(mem_ba(j,1))
!                
!                write(*,*) "bonded_type_angle",array2(1:2),n_botyp+added_ba+2,n_angltyp+added_ba+2
!                write(*,*) "angle",n_angltyp+added_ba+2,ba_typ(ba_lst(mem_ba(k,2),4)),ba_par(ba_lst(mem_ba(k,2),4),1),&
!                &bang(mem_ba(k,1))
!                mem_new_ba(added_ba+2)=bang(mem_ba(k,1))
!                print *,"replace bond type to",n_botyp+added_ba+2,"in biotype",b_type(mem_ba(k,1))
!                added_ba=added_ba+2
!                
!            else
!                print *,"replace bond type to",n_botyp+added_ba+1,"in biotype",b_type(mem_ba(j,1))
!                write(*,*) "bonded_type_angle",array2(1:2),n_botyp+added_ba+1,n_angltyp+added_ba+1
!                write(*,*) "angle",n_angltyp+added_ba+1,ba_typ(ba_lst(k,4)),ba_par(ba_lst(k,4),1),bang(mem_ba(k,1))
!                mem_new_ba(added_ba+1)=bang(mem_ba(k,1))
!                added_ba=added_ba+1
!                
!            end if 
!        end if 
!    end do 
!    
!    
!end do     
!  if ((added_ba.ne.0).or.(added_bo.ne.0).or.(need_fix)) then 
!    print *,"Please fix the parameter file bond angles before continuing"
!    call fexit()
!  end if 
!  deallocate(memory)
!  deallocate(array1)
!  deallocate(array2)
!
!  allocate(memory(2))
!  allocate(array1(2))
!  allocate(array2(2))
!
!! Second part, bond length
!  do i=1,n
!    if (iz(1,i).gt.0) then 
!        array1(1)=bio_botyp(b_type(i))
!        array1(2)=bio_botyp(b_type(iz(1,i)))
!    else
!        cycle 
!    end if 
!
!    do k=1,bo_lstsz      
!        array2(1)=bo_lst(k,1)
!        array2(2)=bo_lst(k,2)
!        test=.true.
!        memory(:)=test_combi(array1,array2,size(memory(:)))
!        do j=1,size(memory(:))
!            if (memory(j).eq.0) then
!                test=.false.
!            end if 
!        end do 
!        if (test.eqv..true.) then 
!
!            bond_index(i)=k
!            
!            if (blen(i).ne.bo_par(bo_lst(k,3),2)) then 
!                changed_bo=changed_bo+1
!                
!                mem_bo(changed_bo,1)=i
!                mem_bo(changed_bo,2)=k
!                
!                ! There I am going to see if that required value already exist. If it does not, readily print the lines that need modifying
!                do h=1,n_bondtyp
!                    if (blen(i).eq.bo_par(h,2)) then
!                        mem_bo(changed_bo,3)=h
!                        exit
!                    end if 
!                end do 
!                
!                if (h-1.eq.n_bondtyp) then 
!                    print *,"Could not find an already existing bond value."
!                    print *,"Please add the following line to the bond definition"
!                    
!                    write(*,*) "bond",n_bondtyp+added_bo+1,bo_typ(bo_lst(k,3)),bo_par(bo_lst(k,3),1),blen(i)
!                    write(*,*) "bonded_type_bond",array2(:),n_bondtyp+added_bo+1
!                    print *,"And replace bond  ",bo_lst(k,3),"by",n_bondtyp+added_bo+1
!                    added_bo=added_bo+1
!                end if 
!            else 
!                used_bo(k)=1
!            end if 
!            if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then 
!                blen_limits(i,1)=bo_par(bo_lst(k,3),2)
!            else 
!                blen(i)=bo_par(bo_lst(k,3),2)
!            end if 
!            exit 
!        end if 
!    end do    
!
!    if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then 
!        if (iz(1,i).gt.0) then 
!            array1(1)=bio_botyp(transform_table(b_type(i)))
!            array1(2)=bio_botyp(transform_table(b_type(iz(1,i))))
!        else
!            cycle ! Martin : this is not necessary because it already does it in the previous loop
!        end if 
!
!        do k=1,bo_lstsz      
!            array2(1)=bo_lst(k,1)
!            array2(2)=bo_lst(k,2)
!
!            test=.true.
!            memory(:)=test_combi(array1,array2,size(memory(:)))
!
!            do j=1,size(memory(:))
!                if (memory(j).eq.0) then
!                    test=.false.
!                end if 
!            end do 
!            if (test.eqv..true.) then 
!                blen_limits(i,2)=bo_par(bo_lst(k,3),2)
!                exit
!            end if 
!        end do  
!        if (nhis.ne.0) then 
!            if (iz(1,i).gt.0) then 
!                array1(1)=bio_botyp(transform_table(his_eqv_table(b_type(i))))
!                array1(2)=bio_botyp(transform_table(his_eqv_table(b_type(iz(1,i)))))
!            else
!                cycle ! Martin : this is not necessary because it already does it in the previous loop
!            end if 
!
!            do k=1,bo_lstsz      
!                array2(1)=bo_lst(k,1)
!                array2(2)=bo_lst(k,2)
!
!                test=.true.
!                memory(:)=test_combi(array1,array2,size(memory(:)))
!
!                do j=1,size(memory(:))
!                    if (memory(j).eq.0) then
!                        test=.false.
!                    end if 
!                end do 
!                if (test.eqv..true.) then 
!                    blen_limits(i,3)=bo_par(bo_lst(k,3),2)
!                    exit
!                end if 
!            end do  
!        end if 
!    end if 
!end do 
!do j=1,changed_bo 
!    do k=1, changed_bo
!        ! Here i need to identify how many values share this line
!        if ((mem_bo(j,3).eq.mem_bo(k,3)).and.(blen(mem_bo(j,1)).ne.blen(mem_bo(k,1)))) then 
!            found=.false.
!            do l=1,added_bo
!                if ((mem_new_bo(l).eq.blen(mem_bo(j,1))).or.(mem_new_bo(l).ne.blen(mem_bo(k,1)))) then 
!                    found=.true.
!                    exit
!                end if 
!            end do 
!            if (found.eqv..true.) then 
!                exit
!            end if 
!            ! If this one is used already, it needs two more (one for each new target values)
!            if (used_bo(mem_bo(j,2)).eq.1) then 
!                ! Change each of the bond types, and add one for each
!                write(*,*) "bonded_type_bond",array2(1),n_botyp+added_bo+1,n_bondtyp+added_bo+1
!                write(*,*) "bond",n_bondtyp+added_bo+1,bo_typ(bo_lst(mem_bo(j,2),3)),bo_par(bo_lst(mem_bo(j,2),3),1),&
!                &blen(mem_bo(j,1))
!
!                mem_new_bo(added_bo+1)=blen(mem_bo(j,1))
!                print *,"replace bond type to",n_botyp+added_bo+1,"in biotype",b_type(mem_bo(j,1))
!                
!                write(*,*) "bonded_type_bond",array2(1),n_botyp+added_bo+2,n_angltyp+added_bo+2
!                write(*,*) "bond",n_bondtyp+added_bo+2,bo_typ(bo_lst(mem_bo(k,2),3)),bo_par(bo_lst(mem_bo(k,2),3),1),&
!                &blen(mem_bo(k,1))
!                mem_new_bo(added_bo+2)=blen(mem_bo(k,1))
!                print *,"replace bond type to",n_botyp+added_bo+2,"in biotype",b_type(mem_bo(k,1))
!                added_bo=added_bo+2
!                
!            else
!                print *,"replace bond type to",n_botyp+added_bo+1,"in biotype",b_type(mem_bo(j,1))
!                write(*,*) "bonded_type_bond",array2(1),n_botyp+added_bo+1,n_bondtyp+added_bo+1
!                write(*,*) "bond",n_bondtyp+added_bo+1,bo_typ(bo_lst(k,3)),bo_par(bo_lst(k,3),1),blen(mem_bo(k,1))
!                mem_new_bo(added_bo+1)=blen(mem_bo(k,1))
!                added_bo=added_bo+1
!                
!            end if 
!        end if 
!    end do 
!end do     
!! Martin : This is the special dummy atom case
!    if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then 
!      do i=1,n
!          if ((blen_limits(i,2).eq.0.).and.(blen_limits(i,1).ne.0.)) then 
!            blen_limits(i,2)=blen_limits(i,1)
!          else if ((blen_limits(i,1).eq.0.).and.(blen_limits(i,2).ne.0.)) then 
!            blen_limits(i,1)=blen_limits(i,2)
!          end if 
!          if ((blen_limits(i,3).eq.0.).and.(blen_limits(i,1).ne.0.)) then 
!            blen_limits(i,3)=blen_limits(i,1)
!          end if 
!      end do 
!    end if 
!    ! Martin : scale the values now 
!  if ((do_pka.eqv..true.).or.(do_hsq.eqv..true.)) then 
!      blen(:)= blen_limits(:,2)*par_pka(2) + blen_limits(:,1)*(1.0-par_pka(2))
!      bang(:) = bang_limits(:,2)*par_pka(2) + bang_limits(:,1)*(1.0-par_pka(2))
!  end if 
!
!end  
!!
!!------------------------------------------------------------------------------------
!
! the following two subroutines populate the Z-matrix differently: by extracting
! the internals from the Cartesian coordinates
!
subroutine genzmat(imol)
!
  use atoms
  use iounit
  use zmatrix
  use molecule
!
  implicit none
!
  integer i,i1,i2,i3,i4,i5,imol
  integer iz0(0:n),iz1(n)
  RTYPE getbang,getztor,dasign
!
!
! zero out the values of the local defining atoms
!
  i1 = 0
  i2 = 0
  i3 = 0
  i4 = 0
  i5 = 0
  iz0(0) = 0
  do i = atmol(imol,1),atmol(imol,2)
    iz0(i) = 0
    iz1(i) = 0
  end do
!
! zero out the final internal coordinate values

  do i = atmol(imol,1),atmol(imol,2)
    blen(i) = 0.0d0
    bang(i) = 0.0d0
    ztor(i) = 0.0d0
  end do
!
! first, decide which of the atoms to define next
!
  do i = atmol(imol,1),atmol(imol,2)
    i1 = i
!
! define the bond length for the current atom
!
    if (i.ge.(atmol(imol,1)+1)) then
      i2 = iz(1,i1)
      blen(i1) = sqrt((x(i1)-x(i2))**2 + &
 &                    (y(i1)-y(i2))**2 +&
 &                    (z(i1)-z(i2))**2)
    end if

! define the bond angle for the current atom
!
    if (i.ge.(atmol(imol,1)+2)) then
      i3 = iz(2,i1)
      bang(i1) = getbang(i1,i2,i3)
    end if

!
! decide whether to use a dihedral or second bond angle;
! then find the value of the angle
!
    if (i.ge.(atmol(imol,1)+3)) then
      if (iz(4,i1) .eq. 0) then
        i4 = iz(3,i1)
        i5 = 0
        ztor(i1) = getztor(i1,i2,i3,i4)
      else
        i4 = iz(3,i1)
        ztor(i1) = getbang(i1,i2,i4)
        i5 = 1
        dasign = getztor(i1,i2,i3,i4)
        if (dasign .gt. 0.0d0)  i5 = -1
      end if
    end if

! transfer defining atoms to the permanent array iz
!
    iz(1,i1) = iz0(i2)
    iz(2,i1) = iz0(i3)
    iz(3,i1) = iz0(i4)
    iz(4,i1) = i5
    iz0(i1) = i
    iz1(i1) = i2
  end do


end
!
!
!------------------------------------------------------------------------------------
!
! this modified routine is only used in readpdb.f, for the pdb-analysis mode, which
! has to be able to correctly reflect varying geometries ...
! it also employs a bunch of sanity checks which allow it to probe a system for
! which bond lengths and angles are ill-defined initially, and to replace the latter
! with CAMPARI default values (useful for partial re-building)
!
subroutine genzmat_f(imol,steppi)
!
  use atoms
  use iounit
  use zmatrix
  use molecule
  use pdb
  use system
  use fyoc
  use params
!
  implicit none
!
  integer i,i1,i2,i3,i4,i5,imol,steppi,rs,i5bu
  integer iz0(0:n),iz1(n)
  RTYPE getbang_safe,getztor,dasign,blbu,babu,dhbu,anv,bll,blh
  logical flagi1,flagi1f
!
! some initialization
  if (steppi.eq.1) then
    anv = pdb_tolerance(3)
    bll = pdb_tolerance(1)
    blh = pdb_tolerance(2)
! if we're referencing back to the last structures, the backup-
! values won't necessarily hold values close to the minimum -> bump up tolerance
  else if (pdb_analyze.EQV..true.) then
    anv = 2.0*pdb_tolerance(3)
    bll = max(1.0 - 2.0*(1.0-pdb_tolerance(1)),0.0d0)
    blh = max(1.0 - 2.0*(1.0-pdb_tolerance(2)),0.0d0)
  end if
!
! zero out the values of the local defining atoms
!
  i1 = 0
  i2 = 0
  i3 = 0
  i4 = 0
  i5 = 0
  iz0(0) = 0
  do i = atmol(imol,1),atmol(imol,2)
    iz0(i) = 0
    iz1(i) = 0
  end do
!
! first, decide which of the atoms to define next
!
  do i = atmol(imol,1),atmol(imol,2)
    blbu = blen(i)
    blen(i) = 0.0
    babu = bang(i)
    bang(i) = 0.0
    dhbu = ztor(i)
    ztor(i) = 0.0
    i1 = i
    flagi1 = .false.
    flagi1f = .false.
    i5bu = iz(4,i1)
    i5 = 0
    rs = atmres(i1)
!
! define the bond length for the current atom
!
    if (i.ge.(atmol(imol,1)+1)) then
      i2 = iz(1,i1)
      blen(i1) = sqrt((x(i1)-x(i2))**2 + (y(i1)-y(i2))**2 + (z(i1)-z(i2))**2)
      if ((steppi.eq.1).OR.(pdb_analyze.EQV..true.)) then
        if(((blen(i1)/blbu).gt.blh).OR.((blen(i1)/blbu).lt.bll))then
          write(ilog,*) 'Bond length exception at atoms ',i1,'and ',i2,'.'
          if (i.ge.(atmol(imol,1)+3)) then
            flagi1 = .true.
          else
!            write(ilog,*) 'WARNING: Atom ',i1,' will NOT BE ADJUSTED!'
          end if
        end if
        if ((pdb_hmode.eq.2).AND.(mass(i).lt.3.5)) then
          flagi1 = .true.
        end if
      end if
    end if
!
! define the bond angle for the current atom
!
    if (i.ge.(atmol(imol,1)+2)) then
      i3 = iz(2,i1)
      bang(i1) = getbang_safe(i1,i2,i3)
      if ((steppi.eq.1).OR.(pdb_analyze.EQV..true.)) then
        if ((abs(bang(i1)-babu).gt.anv).OR.(bang(i1).lt.0.0)) then
          bang(i1) = max(bang(i1),0.0)
          write(ilog,*) 'Bond angle exception at atoms ',i1,'and ',i2,'and ',i3,'.'
          if (i.ge.(atmol(imol,1)+3)) then
            flagi1 = .true.
          else
!            write(ilog,*) 'WARNING: Atom ',i1,' will NOT BE ADJUSTED!'
          end if
        end if
        if ((pdb_hmode.eq.2).AND.(mass(i).lt.3.5)) then
          flagi1 = .true.
        end if
      else if (bang(i1).lt.0.0) then
        bang(i1) = max(bang(i1),0.0)
      end if
    end if
!
! decide whether to use a dihedral or second bond angle;
! then find the value of the angle
!
    if (i.ge.(atmol(imol,1)+3)) then
      if (iz(4,i1).eq.0) then
        i4 = iz(3,i1)
        i5 = 0
        ztor(i1) = getztor(i1,i2,i3,i4)
        if ((steppi.eq.1).OR.(pdb_analyze.EQV..true.)) then
          if (((pdb_hmode.eq.2).OR.(flagi1.EQV..true.)).AND.(mass(i).lt.3.5)) then
            if (fline2(rs).eq.i1) then
              flagi1f = .true.
              if (ztor(fline(rs)).gt.0.0) then
                ztor(i1) = ztor(fline(rs)) - 180.0
              else
                ztor(i1) = ztor(fline(rs)) + 180.0
              end if
            else
              flagi1 = .true.
            end if
          end if
!
!         a problem with the y(f)line/y(f)line2 setup is that an actual dihedral
!         which is NOT a degree of freedom will be used
!         if it is ever requested to correct oxygen positions, use construct below (untested)
!
!          if (whatever) then
!            else if ((rs.gt.rsmol(imol,1)).AND.
! &                          (yline2(rs-1).eq.i1)) then
!              write(*,*) 'in with ',i1
!              if (ztor(yline(rs-1)).gt.0.0) then
!                ztor(i1) = ztor(yline(rs-1)) - 180.0
!              else
!                ztor(i1) = ztor(yline(rs-1)) + 180.0
!              end if
!            end if
!          end if
        end if
      else
        i4 = iz(3,i1)
        ztor(i1) = getbang_safe(i1,i2,i4)
        i5 = 1
        dasign = getztor(i1,i2,i3,i4)
        if (dasign .gt. 0.0d0)  i5 = -1
        if ((steppi.eq.1).OR.(pdb_analyze.EQV..true.)) then
          if ((abs(ztor(i1)-dhbu).gt.anv).OR.(ztor(i1).lt.0.0)) then
            ztor(i1) = max(ztor(i1),0.0)
            write(ilog,*) 'Bond angle exception at atoms ',i1,'and ',i2,'and ',i4,'.'
            if (i.ge.(atmol(imol,1)+3)) then
              flagi1 = .true.
            else
!              write(ilog,*)'WARNING: Atom ',i1,' will NOT BE ADJUSTED!'
            end if
          end if
          if ((pdb_hmode.eq.2).AND.(mass(i).lt.3.5)) then
            flagi1 = .true.
          end if
        else if (ztor(i1).lt.0.0) then
          ztor(i1) = max(ztor(i1),0.0)
        end if
      end if
    end if
!
!   for a re-builder, restore every(!)thing except fline2 torsion (which is already set)
    if (flagi1.EQV..true.) then
      blen(i1) = blbu
      bang(i1) = babu
      if (flagi1f.EQV..false.) ztor(i1) = dhbu
      i5 = i5bu
    end if
!
!   force re-building message for phi-case
    if (flagi1f.EQV..true.) flagi1 = .true.
!
! transfer defining atoms to the permanent array iz
!
    iz(1,i1) = iz0(i2)
    iz(2,i1) = iz0(i3)
    iz(3,i1) = iz0(i4)
    iz(4,i1) = i5
    iz0(i1) = i
    iz1(i1) = i2
    if ((flagi1.EQV..true.).AND.(steppi.eq.1)) then
      write(ilog,*) 'Re-building atom: ',i1,' (ty: ',b_type(i),' = '&
 &,bio_code(b_type(i)),').'
      
!     note we'll do that below (call to makexyz_formol)
    end if
  end do
!
  if ((steppi.eq.1).OR.(pdb_analyze.EQV..true.)) then
    do rs=rsmol(imol,1),rsmol(imol,2)
!     because of the way, phi- and psi-angles are coded up, there are two associated
!     ztor()-entries for each
!     in case the geometry comes from pdb we cannot assume what their (static) relationship
!     is and need to determine it explicitly
      if ((fline(rs).gt.0).AND.(fline2(rs).gt.0)) then
        phish(rs) = ztor(fline2(rs)) - ztor(fline(rs))
        if (phish(rs).gt.180.0) phish(rs) = phish(rs) - 360.0
        if (phish(rs).le.-180.0) phish(rs) = phish(rs) + 360.0
      end if
      if ((yline(rs).gt.0).AND.(yline2(rs).gt.0)) then 
        if (ztor(yline(rs)).lt.0.0) then
          psish(rs) = ztor(yline2(rs)) - (ztor(yline(rs)) + 180.0)
        else
          psish(rs) = ztor(yline2(rs)) - (ztor(yline(rs)) - 180.0)
        end if
        if (psish(rs).gt.180.0) psish(rs) = psish(rs) - 360.0
        if (psish(rs).le.-180.0) psish(rs) = psish(rs) + 360.0
      end if
    end do
  end if
!
  call zmatfyc()
  call makexyz_formol(imol)
!
end
!
!----------------------------------------------------------------------------------
!
! a finally a routine which prints the Z-matrix 
!
subroutine prtzmat(iu)
!
  use atoms
  use molecule
  use zmatrix
  use iounit
!
  implicit none
!
  integer i,k,iu,nnn,kk,imol
  logical isopen
!
  inquire (unit=iu,opened=isopen)
  if (isopen.EQV..false.) then
    write(ilog,*) 'Fatal. Got a stale filehandle in prtzmat(...).'
    call fexit()
  end if
!
! these formattings are chosen to be compatible with TINKER (as much as possible
! since molecule handling is different) 
!
 51   format (i6)
 52   format (i6,2x,a3,i5)
 53   format (i6,2x,a3,i5,i6,f10.5)
 54   format (i6,2x,a3,i5,i6,f10.5,i6,f10.4)
 55   format (i6,2x,a3,i5,i6,f10.5,i6,f10.4,i6,f10.4,i6)
 56   format (2i6)
!
  write(iu,51) n
  do imol=1,nmol
!
!   first three are special cases
!
    nnn = atmol(imol,2)-atmol(imol,1)+1
    if (nnn.ge.1) then
      kk = atmol(imol,1)
      write(iu,52)  kk,atnam(kk),attyp(kk)
    end if
    if (nnn.ge.2) then
      kk = atmol(imol,1) + 1
      write(iu,53)  kk,atnam(kk),attyp(kk),iz(1,kk),blen(kk)
    end if
    if (nnn.ge.3) then
      kk = atmol(imol,1) + 2
      write(iu,54)  kk,atnam(kk),attyp(kk),iz(1,kk),blen(kk),&
 &                    iz(2,kk),bang(kk)
    end if
!
!   now standard rest
!
    do i=atmol(imol,1)+3,atmol(imol,2)
      write(iu,55)  i,atnam(i),attyp(i),iz(1,i),blen(i),&
 &                    iz(2,i),bang(i),iz(3,i),ztor(i),iz(4,i)
    end do
!
  end do
!
! ring closures
!
  if (nadd.ne.0)  write(iu,*)
  do i=1,nadd
    write(iu,56) (iadd(k,i),k=1,2)
  end do
!
end
!
!------------------------------------------------------------------------------
!
! a subroutine to build a working (and - if possible - chemically reasonable) Z-matrix from PDB for an
! unknown entity; also sets valency, atom types (through guess_types) and detects ring bonds
! works for intramolecular linkages as well
!
subroutine infer_topology(i,nl,unkbnd,xyzdum,d2max,pdbnm,rsnam)
!
  use pdb
  use zmatrix
  use sequen
  use molecule
  use math
  use params
  use atoms
  use polypep
  use iounit
  use aminos
!
  implicit none
!
  integer, INTENT(IN) :: i,nl
  integer, INTENT(INOUT) :: unkbnd(n_pdbunk,4)
  RTYPE xyzdum(3,nl),d2max(nl)
  character(4), INTENT(IN) :: pdbnm(nl)
  character(3), INTENT(IN) :: rsnam(nl)
  integer j,k,ii,kk,rs,dihk,knty,rcnt,startj
  integer, ALLOCATABLE:: invpdbmap(:),btyh(:)
  RTYPE d2,d2t
  logical dihthere,superwarn
!
 77 format('Map extracted from PDB is incomplete.'/'Coordinates for atom ',i10,' of biotype ',a3,' in residue ',i8,' (',a3,') &
 &could not be read.'/'This is fatal for simulations with unsupported residues.')
  allocate(invpdbmap(n))
  invpdbmap(:) = 0
  do j=1,n
    if (pdbmap(j).gt.0) then
      invpdbmap(pdbmap(j)) = j
    else
!     there are certain cases where we can tolerate pdbmap(j) being 0
      if (atmres(j).le.(maxval(unkbnd(:,3)))) then
        write(ilog,77) j,bio_code(b_type(j)),atmres(j),amino(seqtyp(atmres(j)))
        call fexit()
      else if ((atmres(j).eq.(maxval(unkbnd(:,3))+1)).AND.(molofrs(atmres(j)).eq.molofrs(min(nseq,maxval(unkbnd(:,3))+1)))) then
        write(ilog,77) j,bio_code(b_type(j)),atmres(j),amino(seqtyp(atmres(j)))
        call fexit()
      end if
    end if
  end do
  superwarn = .false.
  do j=1,n
    if ((invpdbmap(j).le.0).AND.(superwarn.EQV..false.)) then
      write(ilog,*) 'Warning. Superfluous atoms in pdb input for atoms in supported residues are ignored.'
      superwarn = .true.
!      invpdbmap(j) = j
    end if
  end do
  rs = unkbnd(i,3)
 !
  knty = 0
  do j=1,i-1
    if (rsnam(unkbnd(i,1))(1:3).eq.rsnam(unkbnd(j,1))(1:3)) then
      knty = j
      unkbnd(i,4) = unkbnd(knty,3)
      exit
    end if
  end do
!
! sanity
  if (knty.gt.0) then
    k = 0
    if ((unkbnd(i,2)-unkbnd(i,1)).ne.(unkbnd(knty,2)-unkbnd(knty,1))) then
      write(ilog,*) 'Fatal. When using multiple unsupported residues with the same residue name, &
 &it is mandatory that all those residues possess identical numbers and types of atoms in identical &
 &order.'
      call fexit()
    end if
    if (((unkbnd(knty,3).eq.rsmol(molofrs(unkbnd(knty,3)),1)).AND.(unkbnd(i,3).ne.rsmol(molofrs(unkbnd(i,3)),1))).OR.&
 &      ((unkbnd(knty,3).eq.rsmol(molofrs(unkbnd(knty,3)),2)).AND.(unkbnd(i,3).ne.rsmol(molofrs(unkbnd(i,3)),2))).OR.&
 &      ((unkbnd(knty,3).ne.rsmol(molofrs(unkbnd(knty,3)),2)).AND.(unkbnd(i,3).eq.rsmol(molofrs(unkbnd(i,3)),2))).OR.&
 &      ((unkbnd(knty,3).ne.rsmol(molofrs(unkbnd(knty,3)),2)).AND.(unkbnd(i,3).eq.rsmol(molofrs(unkbnd(i,3)),2)))) then
      write(ilog,*) 'Fatal. When using multiple unsupported residues with the same residue name, &
 &it is mandatory that these residues do not occur in both terminal and non-terminal positions.'
      call fexit()
    end if
    do j=unkbnd(i,1),unkbnd(i,2)
      if (pdbnm(j).ne.pdbnm(unkbnd(knty,1)+k)) then
        write(ilog,*) 'Fatal. When using multiple unsupported residues with the same residue name, &
 &it is mandatory that all those residues possess identical numbers and types of atoms in identical &
 &order.'
        call fexit()
      end if
      k = k + 1
    end do
  end if
!
  allocate(btyh(n_biotyp+unkbnd(i,2)-unkbnd(i,1)+1))
  btyh(:) = 0
  rcnt = 0
  do j=unkbnd(i,1),unkbnd(i,2)
!   note that necessarily j == invpdbmap(j) for this atom range
    if ((j.ne.invpdbmap(j)).OR.(j.ne.pdbmap(j))) then
      write(ilog,*) 'Fatal. This is a bug in infer_topology(...).'
      call fexit()
    end if
    d2 = HUGE(d2)
!   find the nearest atom listed prior that is part of the same molecule (partner for Z-matrix bond)
    startj = 0
    do k=j-1,1,-1
      if (invpdbmap(k).le.0) cycle
      if (atmres(invpdbmap(k)).eq.rsmol(molofrs(unkbnd(i,3)),1)) then
        startj = k
      else if (startj.gt.0) then
        exit
      end if
    end do
    if (startj.eq.0) startj = 1
    ii = 0
    do k=startj,j-1
      if (invpdbmap(k).le.0) cycle
      call dis_bound2(xyzdum(:,j),xyzdum(:,k),d2t)
      if (d2t.lt.d2) then 
        d2 = d2t
        ii = k
      end if
    end do
!   get an estimate of the net valence (this is highly approximate of course and makes very many assumptions
!   about atomistic systems)
    kk = 0
    do k=1,nl
      if (k.eq.j) cycle
      if (invpdbmap(k).le.0) cycle
      call dis_bound2(xyzdum(:,j),xyzdum(:,k),d2t)
      if ((k.ge.unkbnd(i,1)).AND.(k.le.unkbnd(i,2))) then
        if (d2t.lt.d2max(j)*d2max(k)) then
          kk = kk + 1
        end if
      else if (attyp(invpdbmap(k)).le.0) then
        if (d2t.lt.d2max(j)*d2max(k)) then
          kk = kk + 1
        end if
      else
        if (d2t.lt.d2max(j)*0.5*sqrt(lj_sig(attyp(invpdbmap(k)),attyp(invpdbmap(k))))) then
          kk = kk + 1
        end if
      end if
    end do
!
    if (knty.gt.0) then
      b_type(j) = b_type(unkbnd(knty,1)+rcnt)
      attyp(j) = attyp(unkbnd(knty,1)+rcnt)
      mass(j) = mass(unkbnd(knty,1)+rcnt)
      rcnt = rcnt + 1
    else
      dihthere = .false.
      do k=unkbnd(i,1),j-1
!       try to infer (stringently) equivalent atom types
        if ((pdbnm(j)(2:4).eq.pdbnm(k)(2:4)).AND.(iz(1,k).eq.invpdbmap(ii)).AND.&
 &((pdbnm(j)(1:1).eq.'1').OR.(pdbnm(j)(1:1).eq.'2').OR.(pdbnm(j)(1:1).eq.'3').OR.(pdbnm(j)(1:1).eq.'4')).AND.&
 &((pdbnm(k)(1:1).eq.'1').OR.(pdbnm(k)(1:1).eq.'2').OR.(pdbnm(k)(1:1).eq.'3').OR.(pdbnm(k)(1:1).eq.'4')).AND.&
 &(pdbnm(j)(1:1).ne.pdbnm(k)(1:1))) then
          b_type(j) = b_type(k)
          btyh(b_type(j)) = btyh(b_type(j)) + 1
          attyp(j) = attyp(k)
          mass(j) = mass(k)
          dihthere = .true.
          exit
        end if
      end do
      if (dihthere.EQV..false.) then
        call guess_types(kk,mass(j),pdbnm(j),attyp(j),b_type(j))
        btyh(b_type(j)) = btyh(b_type(j)) + 1
!       mass should be matched to atom type (info about guessed mass would be lost -> undesirable)
        mass(j) = lj_weight(attyp(j))
      end if
    end if
!
!    write(*,*) attyp(j),': ',lj_symbol(attyp(j)),' and ',b_type(j),': ',bio_code(b_type(j)),' | ',kk
    iz(1:4,j) = 0
    if ((j.eq.unkbnd(i,1)).AND.(rs.eq.rsmol(molofrs(rs),1))) then
      iz(1:4,j) = 0
    else if ((j.eq.unkbnd(i,1)+1).AND.(rs.eq.rsmol(molofrs(rs),1))) then
      iz(1,j) = invpdbmap(ii) ! i = i
      blen(j) = sqrt(d2)
      iz(2:4,j) = 0
    else
!     the angle reference is set up without any particular considerations (input atom order critical)
      iz(1,j) = invpdbmap(ii)
      blen(j) = sqrt(d2)
      if (iz(1,iz(1,j)).gt.0) then
        iz(2,j) = iz(1,iz(1,j))
        call bondang2(xyzdum(:,pdbmap(iz(2,j))),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,j),bang(j))
        bang(j) = bang(j)*RADIAN
      else
        do k=1,j-1
          if (invpdbmap(k).le.0) cycle
          if (iz(1,k).eq.iz(1,j)) exit
        end do
        iz(2,j) = invpdbmap(k)
        call bondang2(xyzdum(:,pdbmap(iz(2,j))),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,j),bang(j))
        bang(j) = bang(j)*RADIAN
      end if
      if ((j.eq.unkbnd(i,1)+2).AND.(rs.eq.rsmol(molofrs(rs),1))) then
        iz(3:4,j) = 0
      else
!       we should scan whether we already have an atom utilizing the same bond iz(1,j),iz(2,j)
!       if not, this will become a putative d.o.f-type dihedral, otherwise create indirect reference to fix
!       covalent geometry
        dihthere = .false.
        do k=1,j-1
          if (invpdbmap(k).le.0) cycle
          if (((iz(1,k).eq.iz(2,j)).AND.(iz(2,k).eq.iz(1,j))).OR.&
 &            ((iz(1,k).eq.iz(1,j)).AND.(iz(2,k).eq.iz(2,j)))) then
            dihk = k
            dihthere = .true.
            exit
          end if
        end do
        if (dihthere.EQV..false.) then
          if ((iz(1,iz(2,j)).gt.0).AND.(iz(1,iz(2,j)).ne.iz(1,j))) then
            iz(3,j) = iz(1,iz(2,j))
            call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,j),ztor(j))
            ztor(j) = ztor(j)*RADIAN
            iz(4,j) = 0
          else
            do k=1,j-1
              if (invpdbmap(k).le.0) cycle
              if ((iz(1,k).eq.iz(2,j)).AND.(k.ne.iz(1,j))) exit
            end do
            iz(3,j) = invpdbmap(k)
            call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,j),ztor(j))
            ztor(j) = ztor(j)*RADIAN
            iz(4,j) = 0
          end if
          if (iz(3,j).eq.0) then
            do k=1,j-1
              if (invpdbmap(k).le.0) cycle
              if (((iz(1,k).eq.iz(1,j)).OR.(iz(1,iz(1,j)).eq.k)).AND.(k.ne.iz(2,j))) then
                iz(3,j) = invpdbmap(k)
                call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),&
 &                         xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                iz(4,j) = -1
                if (ztor(j).lt.0.0) iz(4,j) = 1
                call bondang2(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                ztor(j) = ztor(j)*RADIAN
                exit
              end if
            end do
          end if
        else
          if ((iz(1,j).eq.iz(1,dihk)).AND.(iz(2,j).eq.(iz(2,dihk))).AND.(iz(3,dihk).gt.0)) then
            iz(3,j) = dihk
            iz(4,j) = 0
!           basically an improper (more stable)
            call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,j),ztor(j))
            ztor(j) = ztor(j)*RADIAN
          end if
          if (iz(3,j).eq.0) then
            do k=1,j-1
              if (invpdbmap(k).le.0) cycle
              if (((iz(1,k).eq.iz(1,j)).OR.(iz(1,iz(1,j)).eq.k)).AND.(k.ne.iz(2,j))) then
                iz(3,j) = k
                call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,j),ztor(j))
                if ((iz(1,k).eq.iz(1,j))) then
                  iz(4,j) = 0
                else            
                  iz(4,j) = -1
                  if (ztor(j).lt.0.0) iz(4,j) = 1
                  call bondang2(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                end if
                ztor(j) = ztor(j)*RADIAN
                exit
              end if
            end do
          end if
          if (iz(3,j).eq.0) then
            if ((iz(1,iz(2,j)).gt.0).AND.(iz(1,iz(2,j)).ne.iz(1,j))) then
              iz(3,j) = iz(1,iz(2,j))
              call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,j),ztor(j))
              ztor(j) = ztor(j)*RADIAN
              iz(4,j) = 0
            else
              do k=1,j-1
                if (invpdbmap(k).le.0) cycle
                if ((iz(1,k).eq.iz(2,j)).AND.(k.ne.iz(1,j))) exit
              end do
              iz(3,j) = k
              call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,j),ztor(j))
              ztor(j) = ztor(j)*RADIAN
              iz(4,j) = 0
            end if
          end if
        end if
      end if
    end if
!    write(*,*) iz(1,j),blen(j),iz(2,j),bang(j),iz(3,j),ztor(j)    
  end do
!
  if (knty.gt.0) then
    do j=unkbnd(i,1),unkbnd(i,2)
      atnam(j)(3:3) = atnam(unkbnd(knty,1)+j-unkbnd(i,1))(3:3)
    end do
  else
    do j=1,n_biotyp
      if (btyh(j).le.1) then
        btyh(j) = -100
      else
        btyh(j) = 0
      end if
    end do
    k = 1
    do j=unkbnd(i,1),unkbnd(i,2)
      if (btyh(b_type(j)).ge.0) then
        btyh(b_type(j)) = btyh(b_type(j)) + 1
        call int2str(btyh(b_type(j)),atnam(j)(3:3),k)
      end if
    end do
  end if
  deallocate(btyh)
!
! the linkage to the next residue also needs to be fixed;
! only do this if next residue is natively supported and not in a new molecule
  if (rs.lt.rsmol(molofrs(rs),2)) then
    if (seqtyp(rs+1).gt.0) then
!     these must be processed in their native order
      do j=unkbnd(i,2)+1,unkbnd(i,2)+at(rs+1)%nbb+at(rs+1)%nsc
!       bond lengths and angles across the linkage are handled pretty rigidly again
        if (iz(1,j).eq.0) then
          d2 = HUGE(d2)
          do k=1,j-1
            if (invpdbmap(k).le.0) cycle
            call dis_bound2(xyzdum(:,pdbmap(j)),xyzdum(:,pdbmap(k)),d2t)
            if (d2t.lt.d2) then 
              d2 = d2t
              ii = k
            end if
          end do
          iz(1,j) = ii
          blen(j) = sqrt(d2)
        end if
        if (iz(2,j).eq.0) then
          if (iz(1,iz(1,j)).gt.0) then
            iz(2,j) = iz(1,iz(1,j))
            call bondang2(xyzdum(:,pdbmap(iz(2,j))),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),bang(j))
            bang(j) = bang(j)*RADIAN
          else
            do k=1,j-1
             if (invpdbmap(k).le.0) cycle
              if (iz(1,k).eq.iz(1,j)) then
                iz(2,j) = invpdbmap(k)
                call bondang2(xyzdum(:,pdbmap(iz(2,j))),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),bang(j))
                bang(j) = bang(j)*RADIAN
                exit
              end if
            end do
          end if
        end if
!        
        if (iz(3,j).eq.0) then
          dihthere = .false.
          do k=1,j-1
            if (invpdbmap(k).le.0) cycle
            if (((iz(1,k).eq.iz(2,j)).AND.(iz(2,k).eq.iz(1,j))).OR.&
 &              ((iz(1,k).eq.iz(1,j)).AND.(iz(2,k).eq.iz(2,j)))) then
              dihthere = .true.
              dihk = k
              exit
            end if
          end do
          if (dihthere.EQV..false.) then
            if ((iz(1,iz(2,j)).gt.0).AND.(iz(1,iz(2,j)).ne.iz(1,j))) then
              iz(3,j) = iz(1,iz(2,j))
              call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),&
 &                       xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
              ztor(j) = ztor(j)*RADIAN
              iz(4,j) = 0
            else
              do k=1,j-1
                if (invpdbmap(k).le.0) cycle
                if ((iz(1,k).eq.iz(2,j)).AND.(k.ne.iz(1,j))) then
                  iz(3,j) = invpdbmap(k)
                  call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),&
 &                           xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                  ztor(j) = ztor(j)*RADIAN
                  iz(4,j) = 0
                  exit
                end if
              end do
            end if
            if (iz(3,j).eq.0) then
              do k=1,j-1
                if (invpdbmap(k).le.0) cycle
                if (((iz(1,k).eq.iz(1,j)).OR.(iz(1,iz(1,j)).eq.k)).AND.(k.ne.iz(2,j))) then
                  iz(3,j) = invpdbmap(k)
                  call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),&
 &                           xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                  iz(4,j) = -1
                  if (ztor(j).lt.0.0) iz(4,j) = 1
                  call bondang2(xyzdum(:,pdbmap(k)),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                  ztor(j) = ztor(j)*RADIAN
                  exit
                end if
              end do
            end if
          else
            if ((iz(1,j).eq.iz(1,dihk)).AND.(iz(2,j).eq.(iz(2,dihk))).AND.(iz(3,dihk).gt.0)) then
              iz(3,j) = dihk
              iz(4,j) = 0
!             basically an improper (more stable)
              call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),xyzdum(:,pdbmap(iz(1,j))),&
 &                      xyzdum(:,pdbmap(j)),ztor(j))
              ztor(j) = ztor(j)*RADIAN
            end if
            do k=1,j-1
              if (invpdbmap(k).le.0) cycle
              if (((iz(1,k).eq.iz(1,j)).OR.(iz(1,iz(1,j)).eq.k)).AND.(k.ne.iz(2,j))) then
                iz(3,j) = invpdbmap(k)
                call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),&
 &                         xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                if ((iz(1,k).eq.iz(1,j))) then
                  iz(4,j) = 0
                else
                  ztor(j) = -1
                  if (ztor(j).lt.0.0) iz(4,j) = 1
                  call bondang2(xyzdum(:,pdbmap(k)),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                end if
                ztor(j) = ztor(j)*RADIAN
                exit
              end if
            end do
            if (iz(3,j).eq.0) then
              if ((iz(1,iz(2,j)).gt.0).AND.(iz(1,iz(2,j)).ne.iz(1,j))) then
                iz(3,j) = iz(1,iz(2,j))
                call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),&
 &                         xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                ztor(j) = ztor(j)*RADIAN
                iz(4,j) = 0
              else
                do k=1,j-1
                  if (invpdbmap(k).le.0) cycle
                  if ((iz(1,k).eq.iz(2,j)).AND.(k.ne.iz(1,j))) then
                    iz(3,j) = invpdbmap(k)
                    call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),&
 &                             xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                    ztor(j) = ztor(j)*RADIAN
                    iz(4,j) = 0
                    exit
                  end if
                end do
              end if
            end if
          end if
        end if
!    write(*,*) iz(1,j),blen(j),iz(2,j),bang(j),iz(3,j),ztor(j)
      end do
!     correct for 
      do j=unkbnd(i,2)+1,unkbnd(i,2)+at(rs+1)%nbb+at(rs+1)%nsc
        if (iz(4,j).eq.0) then
          if ((iz(1,j).gt.unkbnd(i,2)).AND.(iz(2,j).gt.unkbnd(i,2)).AND.(iz(3,j).le.unkbnd(i,2))) then
            do k=j+1,unkbnd(i,2)+at(rs+1)%nbb+at(rs+1)%nsc
              if (invpdbmap(k).le.0) cycle
              if (((iz(1,k).eq.iz(2,j)).AND.(iz(2,k).eq.iz(1,j))).OR.&
 &                ((iz(1,k).eq.iz(1,j)).AND.(iz(2,k).eq.iz(2,j)).AND.(iz(4,k).eq.0).AND.(iz(3,k).ne.j))) then
                if (iz(1,k).eq.iz(2,j)) then
                  iz(3,k) = iz(3,j)
                else if (iz(1,k).eq.iz(1,j)) then
                  iz(3,k) = j                    
                end if
                call dihed(xyzdum(:,pdbmap(iz(3,k))),xyzdum(:,pdbmap(iz(2,k))),&
 &                         xyzdum(:,pdbmap(iz(1,k))),xyzdum(:,pdbmap(k)),ztor(k))
                iz(4,k) = 0
                ztor(k) = ztor(k)*RADIAN         
                exit
              end if
            end do
            exit
          end if
        end if
      end do
    end if
  end if
!
  do j=unkbnd(i,1),unkbnd(i,2)
    do k=j+1,unkbnd(i,2)
      call dis_bound2(xyzdum(:,j),xyzdum(:,k),d2t)
      if ((iz(1,j).ne.k).AND.(iz(1,k).ne.j).AND.(d2t.lt.d2max(j)*d2max(k))) then
        nadd = nadd + 1
        iadd(1,nadd) = j
        iadd(2,nadd) = k
      end if
    end do
  end do
!
  deallocate(invpdbmap)
!
end
!
!-----------------------------------------------------------------------------
!
! a subroutine to correct linkages if an unknown residue was identified as a supported polymer type
! also appends eligible residue lists for move types
!
subroutine correct_topology(i,nl,unkbnd,xyzdum)
!
  use pdb
  use zmatrix
  use sequen
  use molecule
  use math
  use polypep
  use fyoc
  use aminos
  use params
  use movesets
  use atoms
  use system
  use iounit
!
  implicit none
!
  integer rs,imol,shf,shf2,k,cbi,cgi,cdi,c2i,c1i,o4i,c3i,c4i
  integer, INTENT(IN) :: i,nl,unkbnd(n_pdbunk,4)
  RTYPE xyzdum(3,nl),blenp
  logical pckyes,renter
!
  pckyes = .false.
  rs = unkbnd(i,3)
  imol = molofrs(rs)
!
  if (seqpolty(rs).eq.'P') then
    if (rs.gt.rsmol(imol,1)) then
      if ((ci(rs-1).gt.0).AND.(cai(rs-1).gt.0)) then
        iz(2,cai(rs)) = ci(rs-1)
        iz(3,cai(rs)) = cai(rs-1)
        wline(rs) = cai(rs)
        call bondang2(xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(ci(rs-1))),bang(cai(rs)))
        bang(cai(rs)) = RADIAN*bang(cai(rs))
        call dihed(xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(ci(rs-1))),&
 &                                                                        xyzdum(:,pdbmap(cai(rs-1))),ztor(cai(rs)))
        ztor(cai(rs)) = RADIAN*ztor(cai(rs))
      end if
      if (((seqtyp(rs-1).eq.26).AND.(seqpolty(rs-1).eq.'P')).OR.(aminopolty(seqtyp(rs-1)).eq.'P')) then
        if (yline(rs-1).gt.0) then
          if (iz(1,yline(rs-1)).eq.iz(1,ni(rs))) then
            yline2(rs-1) = ni(rs)
          end if
        end if
      end if
      if (ci(rs-1).gt.0) then ! should be redundant check
        iz(2,ci(rs)) = ni(rs)
        iz(3,ci(rs)) = ci(rs-1)
        fline(rs) = ci(rs)
        call bondang2(xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(ni(rs))),bang(ci(rs)))
        bang(ci(rs)) = RADIAN*bang(ci(rs))
        call dihed(xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(ni(rs))),&
 &                                                                        xyzdum(:,pdbmap(ci(rs-1))),ztor(ci(rs)))
        ztor(ci(rs)) = RADIAN*ztor(ci(rs))
      end if
      if (oi(rs).gt.0) then ! should also be redundant
        iz(2,oi(rs)) = cai(rs)
        iz(3,oi(rs)) = ni(rs)
        yline(rs) = oi(rs)
        call bondang2(xyzdum(:,pdbmap(oi(rs))),xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),bang(oi(rs)))
        bang(oi(rs)) = RADIAN*bang(oi(rs))
        call dihed(xyzdum(:,pdbmap(oi(rs))),xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),&
 &                                                                        xyzdum(:,pdbmap(ni(rs))),ztor(oi(rs)))
        ztor(oi(rs)) = RADIAN*ztor(oi(rs))
      end if
      if (hni(rs).gt.0) then
        iz(1,hni(rs)) = ni(rs)
        iz(2,hni(rs)) = cai(rs)
        iz(3,hni(rs)) = ci(rs)
        if ((fline(rs).ne.hni(rs)).AND.(fline(rs).gt.0)) fline2(rs) = hni(rs)
        call dis_bound2(xyzdum(:,pdbmap(hni(rs))),xyzdum(:,pdbmap(ni(rs))),blen(hni(rs)))
        blen(hni(rs)) = sqrt(blen(hni(rs)))
        call bondang2(xyzdum(:,pdbmap(hni(rs))),xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(cai(rs))),bang(hni(rs)))
        bang(hni(rs)) = RADIAN*bang(hni(rs))
        call dihed(xyzdum(:,pdbmap(hni(rs))),xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(cai(rs))),&
 &                                                                        xyzdum(:,pdbmap(ci(rs))),ztor(hni(rs)))
        ztor(hni(rs)) = RADIAN*ztor(hni(rs))
      end if
    end if
    shf = 0
    if (ua_model.ne.0) shf = 1
!   detect if there is a pucker to consider for sampling (exclude N-terminal residues for this)
    if ((at(rs)%nsc.ge.(4-shf)).AND.(rs.gt.rsmol(imol,1))) then
      if ((iz(1,at(rs)%sc(3-shf)).eq.at(rs)%sc(2-shf)).AND.(iz(1,at(rs)%sc(2-shf)).eq.cai(rs)).AND.&
 &        ((iz(1,at(rs)%sc(4-shf)).eq.at(rs)%sc(3-shf)).OR.(iz(1,at(rs)%sc(4-shf)).eq.ni(rs))).AND.(ci(rs-1).gt.0)) then
        call dis_bound2(xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(at(rs)%sc(4-shf))),blenp)
        if ((bio_code(b_type(at(rs)%sc(2-shf)))(1:2).eq.'CB').AND.(bio_code(b_type(at(rs)%sc(3-shf)))(1:2).eq.'CG').AND.&
 &          (bio_code(b_type(at(rs)%sc(4-shf)))(1:2).eq.'CD').AND.(blenp.le.1.7)) then
!         set Z-matrix to necessary format
          cbi = at(rs)%sc(2-shf)
          cgi = at(rs)%sc(3-shf)
          cdi = at(rs)%sc(4-shf)
          iz(1,cbi) = cai(rs)
          iz(2,cbi) = ni(rs)
          iz(3,cbi) = ci(rs-1)
          iz(4,cbi) = 0
          iz(1,cgi) = cbi
          iz(2,cgi) = cai(rs)
          iz(3,cgi) = ni(rs)
          iz(4,cgi) = 0
          iz(1,cdi) = ni(rs)
          iz(2,cdi) = cai(rs)
          iz(3,cdi) = cbi
          iz(4,cdi) = 0
          call dis_bound2(xyzdum(:,pdbmap(cbi)),xyzdum(:,pdbmap(cai(rs))),blen(cbi))
          blen(cbi) = sqrt(blen(cbi))
          call dis_bound2(xyzdum(:,pdbmap(cgi)),xyzdum(:,pdbmap(cbi)),blen(cgi))
          blen(cgi) = sqrt(blen(cgi))
          call dis_bound2(xyzdum(:,pdbmap(cdi)),xyzdum(:,pdbmap(ni(rs))),blen(cdi))
          blen(cdi) = sqrt(blen(cdi))
          call bondang2(xyzdum(:,pdbmap(cbi)),xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(ni(rs))),bang(cbi))
          bang(cbi) = RADIAN*bang(cbi)
          call bondang2(xyzdum(:,pdbmap(cgi)),xyzdum(:,pdbmap(cbi)),xyzdum(:,pdbmap(cai(rs))),bang(cgi))
          bang(cgi) = RADIAN*bang(cgi)
          call bondang2(xyzdum(:,pdbmap(cdi)),xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(cai(rs))),bang(cdi))
          bang(cdi) = RADIAN*bang(cdi)
          call dihed(xyzdum(:,pdbmap(cbi)),xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(ci(rs-1))),ztor(cbi))
          ztor(cbi) = RADIAN*ztor(cbi)
          call dihed(xyzdum(:,pdbmap(cgi)),xyzdum(:,pdbmap(cbi)),xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(ni(rs))),ztor(cgi))
          ztor(cgi) = RADIAN*ztor(cgi)
          call dihed(xyzdum(:,pdbmap(cdi)),xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(cbi)),ztor(cdi))
          ztor(cdi) = RADIAN*ztor(cdi)
          pucline(rs) = ci(rs) 
          do k=1,nadd
            if ((iadd(1,k).eq.cbi).OR.(iadd(1,k).eq.cdi).OR.(iadd(2,k).eq.cbi).OR.(iadd(2,k).eq.cdi)) then
              iadd(1,k) = cgi
              iadd(2,k) = cdi
              exit
            end if
          end do
!         improve setup for atoms directly attached (primarily to make auxiliary bond angle variations consistent between systems)
          do k=at(rs)%bb(1),at(rs)%bb(1)+at(rs)%nsc+at(rs)%nbb
!           find atoms directly attached to the ring with a dihedral Z-matrix setup
            if (((iz(1,k).eq.cbi).OR.(iz(1,k).eq.cgi).OR.(iz(1,k).eq.cdi).OR.&
 &                (iz(1,k).eq.cai(rs)).OR.(iz(1,k).eq.ni(rs))).AND.(iz(4,k).eq.0).AND.(k.ne.ni(rs)).AND.(k.ne.cai(rs)).AND.&
 &               (k.ne.cbi).AND.(k.ne.cgi).AND.(k.ne.cdi).AND.(k.ne.ci(rs))) then
              if ((iz(1,k).eq.ni(rs)).AND.(iz(2,k).eq.cai(rs)).AND.(iz(3,k).ne.cdi).AND.(k.gt.cdi)) iz(3,k) = cdi
              if ((iz(1,k).eq.ni(rs)).AND.(iz(3,k).ne.cai(rs)).AND.(iz(2,k).eq.cdi).AND.(k.gt.cai(rs))) iz(3,k) = cai(rs)
              if ((iz(1,k).eq.cai(rs)).AND.(iz(2,k).eq.ni(rs)).AND.(iz(3,k).ne.cbi).AND.(k.gt.cbi)) iz(3,k) = cbi
              if ((iz(1,k).eq.cai(rs)).AND.(iz(3,k).ne.ni(rs)).AND.(iz(2,k).eq.cbi).AND.(k.gt.ni(rs))) iz(3,k) = ni(rs)
              if ((iz(1,k).eq.cbi).AND.(iz(2,k).eq.cai(rs)).AND.(iz(3,k).ne.cgi).AND.(k.gt.cgi)) iz(3,k) = cgi
              if ((iz(1,k).eq.cbi).AND.(iz(3,k).ne.cai(rs)).AND.(iz(2,k).eq.cgi).AND.(k.gt.cai(rs))) iz(3,k) = cai(rs)
              if ((iz(1,k).eq.cgi).AND.(iz(2,k).eq.cbi).AND.(iz(3,k).ne.cdi).AND.(k.gt.cdi)) iz(3,k) = cdi
              if ((iz(1,k).eq.cgi).AND.(iz(3,k).ne.cbi).AND.(iz(2,k).eq.cdi).AND.(k.gt.cbi)) iz(3,k) = cbi
              if ((iz(1,k).eq.cdi).AND.(iz(2,k).eq.cgi).AND.(iz(3,k).ne.ni(rs)).AND.(k.gt.ni(rs))) iz(3,k) = ni(rs)
              if ((iz(1,k).eq.cdi).AND.(iz(3,k).ne.cgi).AND.(iz(2,k).eq.ni(rs)).AND.(k.gt.cgi)) iz(3,k) = cgi
              call dihed(xyzdum(:,pdbmap(k)),xyzdum(:,pdbmap(iz(1,k))),xyzdum(:,pdbmap(iz(2,k))),xyzdum(:,pdbmap(iz(3,k))),&
 &                       ztor(k))
              ztor(k) = RADIAN*ztor(k)
            end if
          end do
        end if
      end if
    end if
    if (rs.lt.rsmol(imol,2)) then
      if (((seqtyp(rs+1).eq.26).AND.(seqpolty(rs+1).eq.'P')).OR.(aminopolty(seqtyp(rs+1)).eq.'P')) then
        iz(1,ni(rs+1)) = ci(rs)
        iz(2,ni(rs+1)) = cai(rs)
        iz(3,ni(rs+1)) = ni(rs)
        iz(4,ni(rs+1)) = 0
        iz(2,cai(rs+1)) = ci(rs)
        iz(3,cai(rs+1)) = cai(rs)
        iz(4,cai(rs+1)) = 0
        iz(3,ci(rs+1)) = ci(rs)
        iz(4,ci(rs+1)) = 0
        call dis_bound2(xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),blen(ni(rs+1)))
        blen(ni(rs+1)) = sqrt(blen(ni(rs+1)))
        call bondang2(xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),bang(ni(rs+1)))
        bang(ni(rs+1)) = RADIAN*bang(ni(rs+1))
        call bondang2(xyzdum(:,pdbmap(cai(rs+1))),xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),bang(cai(rs+1)))
        bang(cai(rs+1)) = RADIAN*bang(cai(rs+1))
        call dihed(xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),&
 &                                                                     xyzdum(:,pdbmap(ni(rs))),ztor(ni(rs+1)))
        ztor(ni(rs+1)) = RADIAN*ztor(ni(rs+1))
        if (yline(rs).gt.0) yline2(rs) = ni(rs+1)
        call dihed(xyzdum(:,pdbmap(cai(rs+1))),xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),&
 &                                                                        xyzdum(:,pdbmap(cai(rs))),ztor(cai(rs+1)))
        ztor(cai(rs+1)) = RADIAN*ztor(cai(rs+1))
        wline(rs+1) = cai(rs+1)
        call dihed(xyzdum(:,pdbmap(ci(rs+1))),xyzdum(:,pdbmap(cai(rs+1))),xyzdum(:,pdbmap(ni(rs+1))),&
 &                                                                        xyzdum(:,pdbmap(ci(rs))),ztor(ci(rs+1)))
        ztor(ci(rs+1)) = RADIAN*ztor(ci(rs+1))
        fline(rs+1) = ci(rs+1)
        if (hni(rs+1).gt.0) then
          iz(1,hni(rs+1)) = ni(rs+1)
          iz(2,hni(rs+1)) = cai(rs+1)
          iz(3,hni(rs+1)) = ci(rs+1)
          iz(4,hni(rs+1)) = 0
          fline2(rs+1) = hni(rs+1)
          call dihed(xyzdum(:,pdbmap(hni(rs+1))),xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(cai(rs+1))),&
 &                                                                          xyzdum(:,pdbmap(ci(rs+1))),ztor(hni(rs+1)))
          ztor(hni(rs+1)) = RADIAN*ztor(hni(rs+1))
        end if
      else if ((cai(rs+1).le.0).OR.(ni(rs+1).le.0)) then
!       do nothing
      else if (iz(1,cai(rs+1)).eq.ni(rs+1)) then
        iz(1,ni(rs+1)) = ci(rs)
        iz(2,ni(rs+1)) = cai(rs)
        iz(3,ni(rs+1)) = ni(rs)
        iz(4,ni(rs+1)) = 0
        iz(2,cai(rs+1)) = ci(rs)
        iz(3,cai(rs+1)) = cai(rs)
        iz(4,cai(rs+1)) = 0
        call dis_bound2(xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),blen(ni(rs+1)))
        blen(ni(rs+1)) = sqrt(blen(ni(rs+1)))
        call bondang2(xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),bang(ni(rs+1)))
        bang(ni(rs+1)) = RADIAN*bang(ni(rs+1))
        call bondang2(xyzdum(:,pdbmap(cai(rs+1))),xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),bang(cai(rs+1)))
        bang(cai(rs+1)) = RADIAN*bang(cai(rs+1))
        call dihed(xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),&
 &                                                                     xyzdum(:,pdbmap(ni(rs))),ztor(ni(rs+1)))
        ztor(ni(rs+1)) = RADIAN*ztor(ni(rs+1))
        if (yline(rs).gt.0) yline2(rs) = ni(rs+1)
        call dihed(xyzdum(:,pdbmap(cai(rs+1))),xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),&
 &                                                                        xyzdum(:,pdbmap(cai(rs))),ztor(cai(rs+1)))
        ztor(cai(rs+1)) = RADIAN*ztor(cai(rs+1))
        wline(rs+1) = cai(rs+1)
      end if
    end if
  else if (seqpolty(rs).eq.'N') then
    shf = 0
    if ((nuci(rs,5).le.0).AND.(nuci(rs,6).le.0).AND.(nuci(rs,4).gt.0)) shf = 2
!   make sure own linkage to supported residue prior in sequence is ok (note that d.o.f setup is not handled)
    if (rs.gt.rsmol(imol,1)) then
      if ((seqpolty(rs-1).eq.'N').OR.(aminopolty(seqtyp(rs-1)).eq.'N')) nucsline(6,rs-1) = nuci(rs,1)
      if ((seqtyp(rs-1).ne.26).AND.(aminopolty(seqtyp(rs-1)).eq.'N')) then
        shf2 = 0
        if ((seqtyp(rs-1).ge.76).AND.(seqtyp(rs-1).le.87)) shf2 = 2
        iz(2,nuci(rs,1)) = nuci(rs-1,5-shf2)
        iz(3,nuci(rs,1)) = nuci(rs-1,4-shf2)
        iz(4,nuci(rs,1)) = 0
        iz(2,nuci(rs,2)) = nuci(rs-1,6-shf2)
        iz(3,nuci(rs,2)) = nuci(rs-1,5-shf2)
        iz(4,nuci(rs,2)) = 0
        call bondang2(xyzdum(:,pdbmap(nuci(rs,1))),xyzdum(:,pdbmap(nuci(rs-1,6-shf2))),xyzdum(:,pdbmap(nuci(rs-1,5-shf2))),&
 &                    bang(nuci(rs,1)))
        bang(nuci(rs,1)) = RADIAN*bang(nuci(rs,1))
        call dihed(xyzdum(:,pdbmap(nuci(rs,1))),xyzdum(:,pdbmap(nuci(rs-1,6-shf2))),xyzdum(:,pdbmap(nuci(rs-1,5-shf2))),&
 &                                                                     xyzdum(:,pdbmap(nuci(rs-1,4-shf2))),ztor(nuci(rs,1)))
        ztor(nuci(rs,1)) = RADIAN*ztor(nuci(rs,1))
        call bondang2(xyzdum(:,pdbmap(nuci(rs,2))),xyzdum(:,pdbmap(nuci(rs,1))),xyzdum(:,pdbmap(nuci(rs-1,6-shf2))),&
 &                    bang(nuci(rs,2)))
        bang(nuci(rs,2)) = RADIAN*bang(nuci(rs,2))
        call dihed(xyzdum(:,pdbmap(nuci(rs,2))),xyzdum(:,pdbmap(nuci(rs,1))),xyzdum(:,pdbmap(nuci(rs-1,6-shf2))),&
 &                                                                     xyzdum(:,pdbmap(nuci(rs-1,5-shf2))),ztor(nuci(rs,2)))
        ztor(nuci(rs,2)) = RADIAN*ztor(nuci(rs,2))
        nucsline(5-shf2,rs-1) = nuci(rs,2)
      end if
    end if
!   detect if the sugar is sampleable
    if (at(rs)%nsc.ge.3) then
      if ((iz(1,at(rs)%sc(2)).eq.at(rs)%sc(1)).AND.(iz(1,at(rs)%sc(1)).eq.nuci(rs,6-shf)).AND.&
 &        ((iz(1,at(rs)%sc(3)).eq.at(rs)%sc(2)).OR.(iz(1,at(rs)%sc(3)).eq.nuci(rs,5-shf)))) then
        call dis_bound2(xyzdum(:,pdbmap(nuci(rs,5-shf))),xyzdum(:,pdbmap(at(rs)%sc(3))),blenp)
        if ((bio_code(b_type(at(rs)%sc(1)))(1:2).eq.'C2').AND.(bio_code(b_type(at(rs)%sc(2)))(1:2).eq.'C1').AND.&
 &          (bio_code(b_type(at(rs)%sc(3)))(1:2).eq.'O4').AND.(blenp.le.1.7)) then
!         set Z-matrix to necessary format
          c2i = at(rs)%sc(1)
          c1i = at(rs)%sc(2)
          o4i = at(rs)%sc(3)
          c4i = nuci(rs,5-shf)
          c3i = nuci(rs,6-shf)
          iz(1,c2i) = c3i
          iz(2,c2i) = c4i
          iz(3,c2i) = nuci(rs,4-shf) ! C5*
          iz(4,c2i) = 0
          iz(1,c1i) = c2i
          iz(2,c1i) = c3i
          iz(3,c1i) = c4i
          iz(4,c1i) = 0
          iz(1,o4i) = c4i
          iz(2,o4i) = c3i
          iz(3,o4i) = c2i
          iz(4,o4i) = 0
          call dis_bound2(xyzdum(:,pdbmap(c2i)),xyzdum(:,pdbmap(c3i)),blen(c2i))
          blen(c2i) = sqrt(blen(c2i))
          call dis_bound2(xyzdum(:,pdbmap(c1i)),xyzdum(:,pdbmap(c2i)),blen(c1i))
          blen(c1i) = sqrt(blen(c1i))
          call dis_bound2(xyzdum(:,pdbmap(o4i)),xyzdum(:,pdbmap(c4i)),blen(o4i))
          blen(o4i) = sqrt(blen(o4i))
          call bondang2(xyzdum(:,pdbmap(c2i)),xyzdum(:,pdbmap(c3i)),xyzdum(:,pdbmap(c4i)),bang(c2i))
          bang(c2i) = RADIAN*bang(c2i)
          call bondang2(xyzdum(:,pdbmap(c1i)),xyzdum(:,pdbmap(c2i)),xyzdum(:,pdbmap(c3i)),bang(c1i))
          bang(c1i) = RADIAN*bang(c1i)
          call bondang2(xyzdum(:,pdbmap(o4i)),xyzdum(:,pdbmap(c4i)),xyzdum(:,pdbmap(c3i)),bang(o4i))
          bang(o4i) = RADIAN*bang(o4i)
          call dihed(xyzdum(:,pdbmap(c2i)),xyzdum(:,pdbmap(c3i)),xyzdum(:,pdbmap(c4i)),xyzdum(:,pdbmap(nuci(rs,4-shf))),ztor(c2i))
          ztor(c2i) = RADIAN*ztor(c2i)
          call dihed(xyzdum(:,pdbmap(c1i)),xyzdum(:,pdbmap(c2i)),xyzdum(:,pdbmap(c3i)),xyzdum(:,pdbmap(c4i)),ztor(c1i))
          ztor(c1i) = RADIAN*ztor(c1i)
          call dihed(xyzdum(:,pdbmap(o4i)),xyzdum(:,pdbmap(c4i)),xyzdum(:,pdbmap(c3i)),xyzdum(:,pdbmap(c2i)),ztor(o4i))
          ztor(o4i) = RADIAN*ztor(o4i)
          pckyes = .true.
          do k=1,nadd
            if ((iadd(1,k).eq.c2i).OR.(iadd(1,k).eq.o4i).OR.(iadd(2,k).eq.c2i).OR.(iadd(2,k).eq.o4i)) then
              iadd(1,k) = c1i
              iadd(2,k) = o4i
              exit
            end if
          end do
          if (rs.eq.rsmol(imol,2)) then
            renter = .false.
            do k=at(rs)%bb(1),at(rs)%bb(1)+at(rs)%nbb+at(rs)%nsc-1
              if ((iz(1,k).eq.c3i.AND.(iz(2,k).eq.c4i)).AND.(iz(4,k).eq.0).AND.(iz(1,iz(3,k)).ne.c3i)) then
                nucsline(6,rs) = k
                renter = .true.
                exit
              end if
            end do
            if (renter.EQV..false.) pckyes = .false.
          end if
          if (pckyes.EQV..true.) then
!           improve setup for atoms directly attached (primarily to make auxiliary bond angle variations consistent between systems)
            do k=at(rs)%bb(1),at(rs)%bb(1)+at(rs)%nsc+at(rs)%nbb
!             find atoms directly attached to the ring with a dihedral Z-matrix setup
              if (((iz(1,k).eq.c2i).OR.(iz(1,k).eq.c1i).OR.(iz(1,k).eq.o4i).OR.&
 &                  (iz(1,k).eq.c3i).OR.(iz(1,k).eq.c4i)).AND.(iz(4,k).eq.0).AND.(k.ne.c4i).AND.(k.ne.c3i).AND.&
 &                 (k.ne.c2i).AND.(k.ne.c1i).AND.(k.ne.o4i).AND.(k.ne.nuci(rs,4-shf))) then
                if ((iz(1,k).eq.c4i).AND.(iz(2,k).eq.c3i).AND.(iz(3,k).ne.o4i).AND.(k.gt.o4i)) iz(3,k) = o4i
                if ((iz(1,k).eq.c4i).AND.(iz(3,k).ne.c3i).AND.(iz(2,k).eq.o4i).AND.(k.gt.c3i)) iz(3,k) = c3i
                if ((iz(1,k).eq.c3i).AND.(iz(2,k).eq.c4i).AND.(iz(3,k).ne.c2i).AND.(k.gt.c2i)) iz(3,k) = c2i
                if ((iz(1,k).eq.c3i).AND.(iz(3,k).ne.c4i).AND.(iz(2,k).eq.c2i).AND.(k.gt.c4i)) iz(3,k) = c4i
                if ((iz(1,k).eq.c2i).AND.(iz(2,k).eq.c3i).AND.(iz(3,k).ne.c1i).AND.(k.gt.c1i)) iz(3,k) = c1i
                if ((iz(1,k).eq.c2i).AND.(iz(3,k).ne.c3i).AND.(iz(2,k).eq.c1i).AND.(k.gt.c3i)) iz(3,k) = c3i
                if ((iz(1,k).eq.c1i).AND.(iz(2,k).eq.c2i).AND.(iz(3,k).ne.o4i).AND.(k.gt.o4i)) iz(3,k) = o4i
                if ((iz(1,k).eq.c1i).AND.(iz(3,k).ne.c2i).AND.(iz(2,k).eq.o4i).AND.(k.gt.c2i)) iz(3,k) = c2i
                if ((iz(1,k).eq.o4i).AND.(iz(2,k).eq.c1i).AND.(iz(3,k).ne.c4i).AND.(k.gt.c4i)) iz(3,k) = c4i
                if ((iz(1,k).eq.o4i).AND.(iz(3,k).ne.c1i).AND.(iz(2,k).eq.c4i).AND.(k.gt.c1i)) iz(3,k) = c1i
                call dihed(xyzdum(:,pdbmap(k)),xyzdum(:,pdbmap(iz(1,k))),xyzdum(:,pdbmap(iz(2,k))),xyzdum(:,pdbmap(iz(3,k))),&
 &                         ztor(k))
                ztor(k) = RADIAN*ztor(k)
              end if
            end do
          end if
        end if
      end if
    end if
    if (rs.lt.rsmol(imol,2)) then
      if ((nuci(rs+1,1).gt.0).AND.(nuci(rs+1,2).gt.0)) then
        iz(1,nuci(rs+1,1)) = nuci(rs,6-shf)
        iz(2,nuci(rs+1,1)) = nuci(rs,5-shf)
        iz(3,nuci(rs+1,1)) = nuci(rs,4-shf)
        iz(4,nuci(rs+1,1)) = 0
        iz(2,nuci(rs+1,2)) = nuci(rs,6-shf)
        iz(3,nuci(rs+1,2)) = nuci(rs,5-shf)
        iz(4,nuci(rs+1,2)) = 0
        if (nuci(rs+1,3).gt.0) then
          iz(3,nuci(rs+1,3)) = nuci(rs,6-shf)
          iz(4,nuci(rs+1,3)) = 0
        end if
        call dis_bound2(xyzdum(:,pdbmap(nuci(rs+1,1))),xyzdum(:,pdbmap(nuci(rs,6-shf))),blen(nuci(rs+1,1)))
        blen(nuci(rs+1,1)) = sqrt(blen(nuci(rs+1,1)))
        call bondang2(xyzdum(:,pdbmap(nuci(rs+1,1))),xyzdum(:,pdbmap(nuci(rs,6-shf))),xyzdum(:,pdbmap(nuci(rs,5-shf))),&
 &                    bang(nuci(rs+1,1)))
        bang(nuci(rs+1,1)) = RADIAN*bang(nuci(rs+1,1))
        call bondang2(xyzdum(:,pdbmap(nuci(rs+1,2))),xyzdum(:,pdbmap(nuci(rs+1,1))),xyzdum(:,pdbmap(nuci(rs,6-shf))),&
 &                    bang(nuci(rs+1,2)))
        bang(nuci(rs+1,2)) = RADIAN*bang(nuci(rs+1,2))
        call dihed(xyzdum(:,pdbmap(nuci(rs+1,1))),xyzdum(:,pdbmap(nuci(rs,6-shf))),xyzdum(:,pdbmap(nuci(rs,5-shf))),&
 &                                                                     xyzdum(:,pdbmap(nuci(rs,4-shf))),ztor(nuci(rs+1,1)))
        ztor(nuci(rs+1,1)) = RADIAN*ztor(nuci(rs+1,1))
        call dihed(xyzdum(:,pdbmap(nuci(rs+1,2))),xyzdum(:,pdbmap(nuci(rs+1,1))),xyzdum(:,pdbmap(nuci(rs,6-shf))),&
 &                                                                        xyzdum(:,pdbmap(nuci(rs,5-shf))),ztor(nuci(rs+1,2)))
        ztor(nuci(rs+1,2)) = RADIAN*ztor(nuci(rs+1,2))
        nnucs(rs) = nnucs(rs) + 1
        nucsline(nnucs(rs),rs) = nuci(rs+1,2)
        if (nuci(rs+1,3).gt.0) then
          call dihed(xyzdum(:,pdbmap(nuci(rs+1,3))),xyzdum(:,pdbmap(nuci(rs+1,2))),xyzdum(:,pdbmap(nuci(rs+1,1))),&
 &                                                                        xyzdum(:,pdbmap(nuci(rs,6-shf))),ztor(nuci(rs+1,3)))
          ztor(nuci(rs+1,3)) = RADIAN*ztor(nuci(rs+1,3))
          if (nnucs(rs+1).gt.0) nucsline(1,rs+1) = nuci(rs+1,3)
        end if
      end if
    end if
  else if (rs.lt.rsmol(imol,2)) then
    if ((seqtyp(rs+1).ne.26).AND.(aminopolty(seqtyp(rs+1)).eq.'P')) then
      if ((ci(rs).le.0).OR.(iz(1,at(rs+1)%bb(1)).ne.ci(rs))) then
        write(ilog,*) 'Warning. Due to unsupported residue of unsupported polymer type being &
 &immediately N-terminal, the meaning of the phi angle in polypeptide residue ',rs+1,' may be &
 &altered.'
      end if
      if ((ci(rs).le.0).OR.(iz(1,at(rs+1)%bb(1)).ne.ci(rs)).OR.&
 &        (oi(rs).le.0)) then
        write(ilog,*) 'Warning. Due to unsupported residue of unsupported polymer type being &
 &immediately N-terminal, the meaning of the omega angle in polypeptide residue ',rs+1,' may be &
 &altered.'
      else if (iz(1,oi(rs)).ne.ci(rs)) then
        write(ilog,*) 'Warning. Due to unsupported residue of unsupported polymer type being &
 &immediately N-terminal, the meaning of the omega angle in polypeptide residue ',rs+1,' may be &
 &altered.'
      end if
    else if ((seqtyp(rs+1).ne.26).AND.(aminopolty(seqtyp(rs+1)).eq.'N')) then
      if (bio_code(b_type(iz(1,at(rs+1)%bb(1))))(1:2).ne.'C3') then
        write(ilog,*) 'Warning. Due to unsupported residue of unsupported polymer type being &
 &immediately N-terminal, the meaning of the first two nucleic acid backbone angles in polynucleotide &
 &residue ',rs+1,' may be altered.'
      end if
    end if
  end if
!
end
!
!-----------------------------------------------------------------------------
!
! this routine generates the position of atom at based on the Cartesian
! coordinates of the reference atoms (lower index number) and the passed
! internals
! obviously has to be called in order (atom-wise)
! chirality indicates the chirality with respect to the ordering of ref.
! atoms (1 or -1), if set to 0 two bond angles are used instead (generally
! discouraged due to numerical in"stability")
! arguments are passed in distance units and degrees, respectively
! for the first few atoms, not all reference atoms might be provided
! -> build in a standard frame
!
subroutine genxyz(at,i2,bl,i3,ba,i4,baodi,chirality)
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
!
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

  if (i2.eq.0) then
      ! third in (+)xz-plane

    x(at) = 0.0d0
    y(at) = 0.0d0
    z(at) = 0.0d0
  else if (i3.eq.0) then

    x(at) = x(i2)
    y(at) = y(i2)
    z(at) = z(i2) + bl
  else if (i4.eq.0) then

    stm1 = sin(ba/RADIAN)
    ctm1 = cos(ba/RADIAN)
    dv2(1) = x(i2) - x(i3)
    dv2(2) = y(i2) - y(i3)
    dv2(3) = z(i2) - z(i3)
    bl2 = sqrt(sum(dv2(:)**2))
    dv2(:) = dv2(:)/bl2
!   rel. xz in the assumed frame
    trx = bl*stm1
    trz = bl2 - bl*ctm1
!   shift/rotate into real frame 
    s1 = sqrt(sum(dv2(1:2)**2))
    if (s1.gt.0.0) then
      c1 = dv2(2)/s1
      s2 = dv2(1)/s1
    else
      c1 = 1.0
      s2 = 0.0
    end if
    x(at) = x(i3) + trx*c1 + trz*s2*s1
    y(at) = y(i3) - trx*s2 + trz*c1*s1
    z(at) = z(i3) + trz*dv2(3)
!
! now the general (and most often used) cases
!
  else

    if (chirality.eq.0) then
!     with a torsional angle
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

!
    end if
!
  end if
!
end
!
!----------------------------------------------------------------
!
! next are a few wrappers for genxyz(...)
!
!----------------------------------------------------------------
!
subroutine makexyz_forsc(rs)
!
  use iounit
  use sequen
  use polypep
  use zmatrix
  use molecule
!
  implicit none
!
  integer i,firstat,rs,i2,i3,i4,chiral,imol
  RTYPE bl,ba,baodi
!
! even though this is slightly redundant, we'll rebuild the whole residue
! regardless
! note that the definition of chi-angles and "sidechain" implies that a part
! of the molecule is branched off and that this is what is rebuilt here
! (i.e., this would NOT work for generic branched polymers)
! it is very important NOT to have parts in backbone arrays that move as a
! result of chi-perturbation
!
  if (at(rs)%nsc.eq.0) then
    write(ilog,*) 'Fatal. Called makexyz_forsc() with ineligible residue (zero sidechain&
 & atoms).'
    call fexit()
  end if
  imol = molofrs(rs)
  firstat = at(rs)%bb(1)
  if (firstat.lt.atmol(imol,1)+3) firstat = atmol(imol,1)+3
  do i=firstat,at(rs)%sc(at(rs)%nsc)
    i2 = iz(1,i)
    i3 = iz(2,i)
    i4 = iz(3,i)
    chiral = iz(4,i)
    bl = blen(i)
    ba = bang(i)
    baodi = ztor(i)
    call genxyz(i,i2,bl,i3,ba,i4,baodi,chiral)
  end do
!
end
!
!---------------------------------------------------------------------
!
subroutine makexyz_forbb(rs)
!
  use iounit
  use sequen
  use polypep
  use zmatrix
  use molecule
!
  implicit none
!
  integer i,rs,i2,i3,i4,chirali,firstat,imol
  RTYPE bl,ba,baodi
!
! rebuild all atoms from first one of pivot residue
!
  imol = molofrs(rs)
  firstat = at(rs)%bb(1)
  if (firstat.lt.atmol(imol,1)+3) firstat = atmol(imol,1)+3
  do i=firstat,atmol(imol,2)
    i2 = iz(1,i)
    i3 = iz(2,i)
    i4 = iz(3,i)
    chirali = iz(4,i)
    bl = blen(i)
    ba = bang(i)
    baodi = ztor(i)
    call genxyz(i,i2,bl,i3,ba,i4,baodi,chirali)
  end do
!
end
!
!---------------------------------------------------------------------
!  
subroutine makexyz_formol(imol)
!
  use iounit
  use sequen
  use polypep
  use zmatrix
  use molecule
!
  implicit none
!
  integer i,i2,i3,i4,chiral,imol
  RTYPE bl,ba,baodi
!
! never build for first three atoms of a molecule (since those are only affected
! by rigid body moves ("feature" of the N->C building paradigm)) 
! alternative solution if this not feasible: maintaining a ref. array/frame
! for each molecule
! 
  do i=atmol(imol,1)+3,atmol(imol,2)
    i2 = iz(1,i)
    i3 = iz(2,i)
    i4 = iz(3,i)
    chiral = iz(4,i)
    bl = blen(i)
    ba = bang(i)
    baodi = ztor(i)
    call genxyz(i,i2,bl,i3,ba,i4,baodi,chiral)
  end do
!
end
!
!---------------------------------------------------------------------
!  
subroutine makexyz_formol2(imol)
!
  use iounit
  use sequen
  use polypep
  use zmatrix
  use molecule
!
  implicit none
!
  integer i,i2,i3,i4,chiral,imol
  RTYPE bl,ba,baodi
!
! initially, the first three atoms have to be built of course, which is done by this fxn
! 
  do i=atmol(imol,1),atmol(imol,2)
    i2 = iz(1,i)
    i3 = iz(2,i)
    i4 = iz(3,i)
    chiral = iz(4,i)
    bl = blen(i)
    ba = bang(i)
    baodi = ztor(i)
    call genxyz(i,i2,bl,i3,ba,i4,baodi,chiral)

  end do
!
end

!
!---------------------------------------------------------------------
!
! simplest of all wrappers
!   
subroutine makexyz2()
!
  use iounit
  use molecule
!
  implicit none
!
  integer imol
!
  do imol=1,nmol
    call makexyz_formol2(imol)
  end do
!
end
!
!-----------------------------------------------------------------------
!
! this concludes the wrappers for genxyz(...)
!
!-----------------------------------------------------------------------
!
! the next routine is related in that it also builds coordinates
! rather than following the Z-matrix hierarchy, however, it assumes only a single
! rotation takes place and applies this by means of a quaternion
!
subroutine quatxyz_forrotlst(ati,alC,dw)
!
  use iounit
  use sequen
  use zmatrix
  use molecule
  use atoms
  use system
  use forces
!
  implicit none
!
  integer i,k,ati,imol,ii,ik
  RTYPE axl(3),mmm,q(4),dw,qinv(4),refv(3)
  logical alC,setalC
!
  if ((izrot(ati)%alsz.le.0).OR.(allocated(izrot(ati)%rotis).EQV..false.)) return
!
  refv(1) = x(iz(1,ati))
  refv(2) = y(iz(1,ati))
  refv(3) = z(iz(1,ati))
  axl(1) = x(iz(1,ati)) - x(iz(2,ati))
  axl(2) = y(iz(1,ati)) - y(iz(2,ati))
  axl(3) = z(iz(1,ati)) - z(iz(2,ati))
  if (alC.EQV..true.) axl(:) = -1.0*axl(:)
!
  mmm = sqrt(sum(axl(:)*axl(:)))
  axl(:) = axl(:)/mmm
!
! construct the quaternion
  q(1) = cos(dw/2.0)
  q(2:4) = axl(1:3)*sin(dw/2.0)
!
! the rotation set may already be flipped, detect this with the help of iz(3,ati)
  setalC = .false.
  do k=izrot(ati)%alsz,1,-1
    ii = izrot(ati)%rotis(k,1)
    ik = izrot(ati)%rotis(k,2)
    if ((iz(3,ati).ge.ii).AND.(iz(3,ati).le.ik)) then
      setalC = .true.
      exit
    end if
  end do
  if (alC.EQV..true.) then
    if (setalC.EQV..true.) then ! setalC true means a C-terminal rot list
      setalC = .false.
    else
      setalC = .true.
    end if
  end if
!  
! construct the multiplicative inverse
  mmm = sum(q(:)*q(:))
  if (mmm.le.0.0) then
    write(ilog,*) 'Fatal error. Got a 0-quat. in quatxyz_forrotlst(...).'
    call fexit()
  end if
  qinv(1) = q(1)/mmm
  qinv(2:4) = -q(2:4)/mmm

  if (setalC.EQV..true.) then
    imol = molofrs(atmres(ati))
    k = atmol(imol,1)
    i = 1
    do while (k.le.atmol(imol,2))
      if (i.le.izrot(ati)%alsz) then
        if ((k.ge.izrot(ati)%rotis(i,1)).AND.(k.le.izrot(ati)%rotis(i,2))) then
          k = izrot(ati)%rotis(i,2)+1
          i = i + 1
        else
          if ((k.ne.iz(1,ati)).AND.(k.ne.iz(2,ati))) then
            call quat_conjugate3(q,qinv,k,refv)
          end if
          k = k + 1
        end if
      else
        if ((k.ne.iz(1,ati)).AND.(k.ne.iz(2,ati))) then
          call quat_conjugate3(q,qinv,k,refv)
        end if
        k = k + 1
      end if
    end do
  else
    do k=1,izrot(ati)%alsz
      ii = izrot(ati)%rotis(k,1)
      ik = izrot(ati)%rotis(k,2)
      do i=ii,ik
        call quat_conjugate3(q,qinv,i,refv)
      end do
    end do
  end if
!
end
!
!---------------------------------------------------------------------
!
! for this next routine, see description below
! note this routine relies on (rather unfortunate) explicit listing
! of residue types via name -> needs to be appended if new residues
! are created
!
subroutine find_ntorsn()
!
  use fyoc
  use iounit
  use sequen
  use torsn
  use zmatrix
  use aminos
  use molecule
  use movesets
  use atoms
!
  implicit none
!
  integer rs
  integer kk,mty,nmt,curmol,buk
  integer, ALLOCATABLE:: molnrs(:)
!
  ntorsn = 0
  ntorpuck = 0
  ndyntorsn = 0
  allocate(molnrs(nmol))
!
! go through the backbone and identify rotatable bonds
! now this was originally meant for a more general torsional implementation,
! i.e., to have arbitrary defined torsions
! since, however, the MC sampler needs to distinguish (efficiencies, movesets,
! etc., all depend on type of torsion) them anyway, we use more descript pointer
! arrays (see fyoc.i) directly into the relevant z-matrix array
! -> this function only serves to establish "ntorsn" and "ntorpuck" and nothing else
!
  curmol = -1
!
  do rs=1,nseq
!
    if (molofrs(rs).ne.curmol) then
      curmol = molofrs(rs)
      if (atmol(molofrs(rs),1).eq.atmol(molofrs(rs),2)) then
        nmt = 3 ! will cycle o'TOR',ut anyway
      else if ((atmol(molofrs(rs),2)-atmol(molofrs(rs),2)).eq.1)&
 &                                                             then
        nmt = 5 ! will cycle out anyway (just for clarity)
      else
        nmt = 6
      end if
    end if
    mty = moltypid(molofrs(rs))
    buk = ntorpuck
!
    if (wline(rs).gt.0) then
      if (izrot(wline(rs))%alsz.gt.0) then
        ntorsn = ntorsn + 1
        ntorpuck = ntorpuck + 1
        nmt = nmt + 1
        wnr(rs) = nmt
        if (moltyp(mty,1).eq.molofrs(rs)) then
          ntormol(mty) = ntormol(mty) + 1
        end if
      end if
    end if
!
    if (fline(rs).gt.0) then
      if (izrot(fline(rs))%alsz.gt.0) then
        ntorsn = ntorsn + 1
        ntorpuck = ntorpuck + 1
        nmt = nmt + 1
        fnr(rs) = nmt
        if (moltyp(mty,1).eq.molofrs(rs)) then
          ntormol(mty) = ntormol(mty) + 1
        end if
!      else
!        ndyntorsn = ndyntorsn + 1
      end if
    end if
!
    if ((pucline(rs).gt.0).OR.(nucsline(6,rs).gt.0)) then
      ntorpuck = ntorpuck + 7
    end if
!
    if (yline(rs).gt.0) then
      if (izrot(yline(rs))%alsz.gt.0) then
        ntorsn = ntorsn + 1
        ntorpuck = ntorpuck + 1
        nmt = nmt + 1
        ynr(rs) = nmt
        if (moltyp(mty,1).eq.molofrs(rs)) then
          ntormol(mty) = ntormol(mty) + 1
        end if
      end if
!
    end if
!
    do kk=1,nnucs(rs)
      if (izrot(nucsline(kk,rs))%alsz.le.0) cycle
      ntorsn = ntorsn + 1
      ntorpuck = ntorpuck + 1
      nmt = nmt + 1
      nucsnr(kk,rs) = nmt
      if (moltyp(mty,1).eq.molofrs(rs)) then
        ntormol(mty) = ntormol(mty) + 1
      end if
    end do
!
    do kk=1,nchi(rs)
      if (izrot(chiline(kk,rs))%alsz.le.0) cycle
      ntorsn = ntorsn + 1
      ntorpuck = ntorpuck + 1
      nmt = nmt + 1
      chinr(kk,rs) = nmt
      if (moltyp(mty,1).eq.molofrs(rs)) then
        ntormol(mty) = ntormol(mty) + 1
      end if
    end do
!
    if (buk.eq.ntorpuck) notors(rs) = .true.
    molnrs(molofrs(rs)) = nmt
!
  end do
!
  do kk=1,unslst%nr
    rs = atmres(iz(1,unslst%idx(kk)))
    notors(rs) = .false.
    ntorsn = ntorsn + 1
    ntorpuck = ntorpuck + 1
    molnrs(molofrs(rs)) = molnrs(molofrs(rs)) + 1
    if (moltyp(moltypid(molofrs(rs)),1).eq.molofrs(rs)) then
      if (othidxmol(moltypid(molofrs(rs))).eq.0) othidxmol(moltypid(molofrs(rs))) = molnrs(molofrs(rs))
      ntormol(moltypid(molofrs(rs))) = ntormol(moltypid(molofrs(rs))) + 1
    end if
  end do
!
  do kk=1,unklst%nr
    rs = atmres(iz(1,unklst%idx(kk)))
    notors(rs) = .false.
    ntorsn = ntorsn + 1
    ntorpuck = ntorpuck + 1
    molnrs(molofrs(rs)) = molnrs(molofrs(rs)) + 1
    if (moltyp(moltypid(molofrs(rs)),1).eq.molofrs(rs)) then
      if (othidxmol(moltypid(molofrs(rs))).eq.0) othidxmol(moltypid(molofrs(rs))) = molnrs(molofrs(rs))
      ntormol(moltypid(molofrs(rs))) = ntormol(moltypid(molofrs(rs))) + 1
    end if
  end do
!
  ndyntorsn = ndyntorsn + ntorsn ! includes proline phi
!
  deallocate(molnrs)
!
end
!
!-----------------------------------------------------------------------
!
! transfer 
!
subroutine transfer_torsn(mode,tvec)
!
  use fyoc
  use iounit
  use sequen
  use zmatrix
  use system
  use polypep
  use torsn
  use movesets
!
  implicit none
!
  integer rs,ttc,mode,shf,shfbu
  integer kk
  RTYPE tvec(ntorpuck)
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
  ttc = 0
!
  do rs=1,nseq
!
    if (wline(rs).gt.0) then
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = omega(rs)
      if (mode.eq.2) omega(rs) = tvec(ttc)
    end if
!
    if ((fline(rs).gt.0).AND.(seqflag(rs).ne.5)) then
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = phi(rs)
      if (mode.eq.2) phi(rs) = tvec(ttc)
    end if
!
    if (yline(rs).gt.0) then
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = psi(rs)
      if (mode.eq.2) psi(rs) = tvec(ttc)
    end if
!
    do kk=1,nnucs(rs)
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = nucs(kk,rs)
      if (mode.eq.2) nucs(kk,rs) = tvec(ttc)
    end do
!
    do kk=1,nchi(rs)
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = chi(kk,rs)
      if (mode.eq.2) chi(kk,rs) = tvec(ttc)
    end do
!
    if (((seqflag(rs).eq.5).AND.(fline(rs).gt.0)).OR.((seqpolty(rs).eq.'N').AND.(nucsline(6,rs).gt.0))) then
      if (seqpolty(rs).eq.'N') then
        shfbu = shf
        shf = 1
      end if
!
      if (seqpolty(rs).eq.'N') then
        ttc = ttc + 1
        if (mode.eq.1) tvec(ttc) = ztor(nucsline(6,rs))
        if (mode.eq.2) ztor(nucsline(6,rs)) = tvec(ttc)
      else
        ttc = ttc + 1
        if (mode.eq.1) tvec(ttc) = ztor(fline(rs))
        if (mode.eq.2) then
          ztor(fline(rs)) = tvec(ttc)
          phi(rs) = tvec(ttc)
        end if
      end if
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = ztor(at(rs)%sc(2-shf))
      if (mode.eq.2) ztor(at(rs)%sc(2-shf)) = tvec(ttc)
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = ztor(at(rs)%sc(3-shf))
      if (mode.eq.2) ztor(at(rs)%sc(3-shf)) = tvec(ttc)
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = ztor(at(rs)%sc(4-shf))
      if (mode.eq.2) ztor(at(rs)%sc(4-shf)) = tvec(ttc)
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = bang(at(rs)%sc(2-shf))
      if (mode.eq.2) bang(at(rs)%sc(2-shf)) = tvec(ttc)
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = bang(at(rs)%sc(3-shf))
      if (mode.eq.2) bang(at(rs)%sc(3-shf)) = tvec(ttc)
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = bang(at(rs)%sc(4-shf))
      if (mode.eq.2) bang(at(rs)%sc(4-shf)) = tvec(ttc)
!
      if (seqpolty(rs).eq.'N') shf = shfbu
    end if
!
  end do
!
  do rs=1,unslst%nr
    ttc = ttc + 1
    if (mode.eq.1) tvec(ttc) = ztor(unslst%idx(rs))
    if (mode.eq.2) ztor(unslst%idx(rs)) = tvec(ttc)
  end do
!
  do rs=1,unklst%nr
    ttc = ttc + 1
    if (mode.eq.1) tvec(ttc) = ztor(unklst%idx(rs))
    if (mode.eq.2) ztor(unklst%idx(rs)) = tvec(ttc)
  end do
!
end
!
!-----------------------------------------------------------------------------
!
! this subroutine populates the 1-2 bonded arrays which form the basis
! for the 1-3 and 1-4 ones (see below)
!
! strategy outline: every molecule effectively maintains an independent Z-matrix,
!                   i.e., would one write out the Z-matrix, the first three atoms
!                   of each molecule would indeed not have a well-defined reference
!                   this is more meaningful than the alternative solution to simply
!                   chain the Z-matrix through molecular boundaries leading to
!                   quite meaningless internals which are difficult to maintain
!
!
subroutine setup_connectivity_2()
!
  use zmatrix
  use iounit
  use molecule
  use atoms
!
  implicit none
!
  integer imol,at,at2,i,j,k
!
  i12(:,:) = 0
  n12(:) = 0
!
  do imol=1,nmol
    do at=atmol(imol,1)+1,atmol(imol,2)
      at2 = iz(1,at)
!     sanity, should be handled by the +1-offset
      if (at2.le.0) then
        write(ilog,*) 'Warning! Missing reference atom for second atom&
 & in molecule ',imol,' (atom #',at,' returned ',at2,'). This is mos&
 &t likely a bug. Please report.'
        cycle
      end if
      n12(at) = n12(at) + 1
      n12(at2) = n12(at2) + 1
      i12(n12(at),at) = at2
      i12(n12(at2),at2) = at
    end do
  end do
!
! cyclic systems require the inclusion of extra 1-2 pairs not apparent
! from the simple Z-matrix term
! these were set up previously in the z-matrix generating routines
!
  do i=1,nadd
    do j=1,2
      k = iadd(j,i)
      n12(k) = n12(k) + 1
      i12(n12(k),k) = iadd(3-j,i)
    end do
  end do
!
! finally sort the list
!
  do at=1,n
    call sort(n12(at),i12(1,at))
  end do
!
end
!
!-----------------------------------------------------------------------
!
! this subroutine populates the 1-3 and 1-4 bonded arrays
!
subroutine setup_connectivity_34()
!
  use atoms
!
  implicit none
!
  integer at,i,ii,j,jj,k,kk,l
  logical docyc
!
! first 1-3
!
  i13(:,:) = 0
  do at=1,n
    n13(at) = 0
    do i=1,n12(at)
      ii = i12(i,at)
      do j=1,n12(ii)
        jj = i12(j,ii)
        docyc = .false.
        if (jj.eq.at) cycle
        do k=1,n12(at)
          if (jj.eq.i12(k,at)) then
            docyc = .true.
            exit
          end if
        end do
        if (docyc.EQV..true.) cycle
        n13(at) = n13(at) + 1
        i13(n13(at),at) = jj
      end do
    end do
    call sort(n13(at),i13(1,at))
  end do
!
! now 1-4
!
  i14(:,:) = 0
  do at=1,n
    n14(at) = 0
    do i=1,n12(at)
      ii = i12(i,at)
      do j=1,n12(ii)
        jj = i12(j,ii)
        do k=1,n12(jj)
          kk = i12(k,jj)
          docyc = .false.
          if (kk.eq.at) cycle
          do l=1,n12(at)
            if (kk.eq.i12(l,at)) then
              docyc = .true.
              exit
            end if
          end do
          do l=1,n13(at)
            if (kk.eq.i13(l,at)) then
              docyc = .true.
              exit
            end if
          end do
          if (docyc.EQV..true.) cycle
          n14(at) = n14(at) + 1
          i14(n14(at),at) = kk
        end do
      end do
    end do
    call sort(n14(at),i14(1,at))
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine valence_check()
!
  use atoms
  use iounit
  use params
!
  implicit none
!
  integer at,k,kk,t1,t2
  logical, ALLOCATABLE:: chkd(:)
!
! since the parameter file has valence information, check whether
! any exceptions occur: good bug-check against builder+parameters
!    
 75   format('Fatal. Atom ',i7,' (biotype ',i5,'), which was assigned LJ&
 &-type ',i4,', has ',i2,' bound partner atoms instead of the ',&
 &i2,' defined for it in the parameter file (',a,'). This is the first error of this type&
 & for LJ-type ',i4,'. Further warnings omitted.')
!
  allocate(chkd(n_ljtyp))
  chkd(:) = .false.
  call strlims(paramfile,t1,t2)
!
  do at=1,n
    k = attyp(at)
    kk = bio_ljtyp(b_type(at))
    if (n12(at).ne.lj_val(k)) then
      if ((kk.ne.k).AND.(lj_patched.EQV..true.).AND.(lj_val(kk).eq.n12(at))) cycle
      if (chkd(k).EQV..false.) then
        write(ilog,75) at,b_type(at),k,n12(at),lj_val(k),&
 &paramfile(t1:t2),k
        chkd(k) = .true.
      end if
      if (be_unsafe.EQV..false.) call fexit()
    end if
  end do
!
  deallocate(chkd)
!
end
!
!-------------------------------------------------------------------------
!
! subroutine to determine whether an atom is part of a ring backbone
! it expects an integer array of size 10 to hold up to 10 rings
!
  subroutine determine_upto6ring(ati,nrgs,rnglst)
!
  use atoms
!
  implicit none
!
  integer ati,nrgs,rnglst(10),i,j,k,l,m
  logical isok
!
  nrgs = 0
!
! 3 
  do i=1,n12(ati)
    do j=i+1,n12(ati)
      do k=1,n12(i12(i,ati))
        if (i12(k,i12(i,ati)).eq.i12(j,ati)) then
          nrgs = nrgs + 1
          rnglst(nrgs) = 3
        end if
      end do
    end do
  end do
! 4 - fails for bicyclo [1.1.0]butane
  do i=1,n13(ati)
    do j=i+1,n13(ati)
      if (i13(i,ati).eq.i13(j,ati)) then
        nrgs = nrgs + 1
        rnglst(nrgs) = 4
      end if
    end do
  end do
! 5 - filters out atoms connected to 3-rings
  do i=1,n13(ati)
    do j=i+1,n13(ati)
      isok = .false.
      do k=1,n12(i13(i,ati))
        if (i12(k,i13(i,ati)).eq.i13(j,ati)) then
          isok = .true.
          exit
        end if
      end do
      if (isok.EQV..false.) cycle
      do k=1,n12(i13(i,ati))
        do l=1,n12(i13(j,ati))
          if (i12(k,i13(i,ati)).eq.i12(l,i13(j,ati))) isok = .false.
        end do
      end do
      if (isok.EQV..false.) cycle
      nrgs = nrgs + 1
      rnglst(nrgs) = 5
    end do
  end do
! 6 - filters out atoms connected to 4-rings
  do i=1,n13(ati)
    do j=i+1,n13(ati)
      isok = .false.
      do k=1,n13(i13(i,ati))
        if (i13(k,i13(i,ati)).eq.i13(j,ati)) then
          isok = .true.
          exit
        end if
      end do
      if (isok.EQV..false.) cycle
      do k=1,n12(ati)
        do l=1,n12(i13(i,ati))
          do m=1,n12(i13(j,ati))
            if ((i12(k,ati).eq.i12(l,i13(i,ati))).AND.(i12(k,ati).eq.i12(m,i13(j,ati)))) isok = .false.
          end do
        end do
      end do
      if (isok.EQV..false.) cycle
      nrgs = nrgs + 1
      rnglst(nrgs) = 6
    end do
  end do
!
  end 
!
!---------------------------------------------------------------------------
!
! this routine scans all bonds in a molecule, identifies those that are freely rotatable, picks a single
! atom to describe rotation around this bond (f,yline/f,yline2 exceptions apply), and allocates memory and
! stores lists of atoms rotating (default building direction) upon changes of the corresponding dihedral Z-matrix entry
! lists are stored as stretches of consecutive atoms
! this information can also be useful in identifying angles that should be eligible, but
! fail due to Z-matrix restrictions (usually caused by setting up multiple rotating atoms with iz(4,i) being 0)
! for this to work with hard crosslinks (crlk_mode is 2), crosslinks would have to first be established, intermolecular
! ones coperturbed (crosslink_follow), and checks as per genints(...) may have to be altered
!
subroutine find_rotlsts(imol)
!
  use zmatrix
  use sequen
  use molecule
  use atoms
  use polypep 
  use fyoc
  use inter
  use aminos
  use params
  use iounit
!
  implicit none
!
  integer imol,rs,tnat,ati,k,i,j,l,zid,zid2
  RTYPE incr,epsa,epsl,zbu,zbu2
  integer, ALLOCATABLE:: blst(:,:)
  RTYPE, ALLOCATABLE:: coords(:,:),lvec(:,:),ivec(:,:),tvec(:,:),avec(:,:)
  logical elig,isnew,proflag
!
  izrot(atmol(imol,1):atmol(imol,2))%alsz = 0
  izrot(atmol(imol,1):atmol(imol,2))%alsz2 = 0
!
 95 format('Bond in question: ',a4,'(',i10,')--',a4,'(',i10,')',/,'Terminal atoms in question: ',a4,'(',i10,') / ',a4,'(',i10,')')
! first make a list of all non-terminal bonds
  incr = 27.
  epsl = 1.0e-6
  epsa = 1.0e-4
  tnat = atmol(imol,2) - atmol(imol,1) + 1
  allocate(blst(tnat,MAXVLNC+1))
  blst(:,1) = 0
!
  do ati=atmol(imol,1),atmol(imol,2)
    if (n12(ati).le.1) cycle
    do k=1,n12(ati)
      if ((ati.lt.i12(k,ati)).AND.(n12(i12(k,ati)).gt.1)) then
        blst(ati-atmol(imol,1)+1,1) = blst(ati-atmol(imol,1)+1,1) + 1
        blst(ati-atmol(imol,1)+1,1+blst(ati-atmol(imol,1)+1,1)) = i12(k,ati)
      end if
    end do
  end do
!
! for each bond, check whether it is rotatable and whether it is already in a pointer array
  allocate(coords(tnat,3))
 
  !print *,"GH"
  !print *,atmol(imol,1),atmol(imol,2)
  coords(1:tnat,1) = x(atmol(imol,1):atmol(imol,2))
  coords(1:tnat,2) = y(atmol(imol,1):atmol(imol,2))
  coords(1:tnat,3) = z(atmol(imol,1):atmol(imol,2))
  allocate(lvec(max(1,sum(nrsbl(rsmol(imol,1):rsmol(imol,2)))),2))
  allocate(avec(max(1,sum(nrsba(rsmol(imol,1):rsmol(imol,2)))),2))
  allocate(ivec(max(1,sum(nrsimpt(rsmol(imol,1):rsmol(imol,2)))),2))
  allocate(tvec(max(1,sum(nrsdi(rsmol(imol,1):rsmol(imol,2)))),2))
  lvec(:,:) = 0.0
  avec(:,:) = 0.0
  ivec(:,:) = 0.0
  tvec(:,:) = 0.0
  call genints(imol,lvec(:,2),avec(:,2),ivec(:,2),tvec(:,2))
!
  do i=1,tnat
    do k=1,blst(i,1)
      zid = 0
      zid2 = 0
      isnew = .false.
      proflag = .false.
      do rs=max(rsmol(imol,1),atmres(blst(i,1+k))-2),min(rsmol(imol,2),atmres(blst(i,1+k))+2)
        if (wline(rs).gt.0) then
          if (((iz(1,wline(rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(2,wline(rs)).eq.blst(i,1+k))).OR.&
 &            ((iz(2,wline(rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(1,wline(rs)).eq.blst(i,1+k)))) then
            zid = wline(rs)
          end if
        end if
        if (fline(rs).gt.0) then
          if (((iz(1,fline(rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(2,fline(rs)).eq.blst(i,1+k))).OR.&
 &            ((iz(2,fline(rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(1,fline(rs)).eq.blst(i,1+k)))) then
            zid = fline(rs)
            if ((seqtyp(rs).eq.9).OR.(seqtyp(rs).eq.25).OR.(seqtyp(rs).eq.32)) then
              proflag = .true.
            end if
            zid2 = fline2(rs)
          end if
        end if
        if (yline(rs).gt.0) then
          if (((iz(1,yline(rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(2,yline(rs)).eq.blst(i,1+k))).OR.&
 &            ((iz(2,yline(rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(1,yline(rs)).eq.blst(i,1+k)))) then
            zid = yline(rs)
            zid2 = yline2(rs)
          end if
        end if
        do j=1,nchi(rs)
          if (((iz(1,chiline(j,rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(2,chiline(j,rs)).eq.blst(i,1+k))).OR.&
 &            ((iz(2,chiline(j,rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(1,chiline(j,rs)).eq.blst(i,1+k)))) then
            zid = chiline(j,rs)
          end if
        end do
        do j=1,nnucs(rs)
          if (((iz(1,nucsline(j,rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(2,nucsline(j,rs)).eq.blst(i,1+k))).OR.&
 &            ((iz(2,nucsline(j,rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(1,nucsline(j,rs)).eq.blst(i,1+k)))) then
            zid = nucsline(j,rs)
          end if
        end do
      end do
      if (zid.eq.0) then
        do j=1,n12(blst(i,1+k))
          if ((iz(1,i12(j,blst(i,1+k))).eq.blst(i,1+k)).AND.(iz(2,i12(j,blst(i,1+k))).eq.(i+atmol(imol,1)-1))) then
            elig = .true.
            do l=1,n12(blst(i,1+k))
              if (i12(l,blst(i,1+k)).eq.iz(3,i12(j,blst(i,1+k)))) then
                elig = .false.
                exit ! improper setup
              end if
            end do
!           further conditions are dihedral setup and having all reference atoms defined
            if ((elig.EQV..true.).AND.(iz(4,i12(j,blst(i,1+k))).eq.0).AND.(minval(iz(1:3,i12(j,blst(i,1+k)))).gt.0)) then
              zid = i12(j,blst(i,1+k))
              isnew = .true.
              exit
            end if
          end if
        end do
        do j=1,n12(i+atmol(imol,1)-1)
          if ((iz(1,i12(j,(i+atmol(imol,1)-1))).eq.(i+atmol(imol,1)-1)).AND.(iz(2,i12(j,(i+atmol(imol,1)-1))).eq.blst(i,1+k))) then
            elig = .true.
            do l=1,n12((i+atmol(imol,1)-1))
              if (i12(l,(i+atmol(imol,1)-1)).eq.iz(3,i12(j,(i+atmol(imol,1)-1)))) then
                elig = .false.
                exit ! improper setup
              end if
            end do
!           further conditions are dihedral setup and having all reference atoms defined
            if ((elig.EQV..true.).AND.(iz(4,i12(j,(i+atmol(imol,1)-1))).eq.0).AND.&
 &              (minval(iz(1:3,i12(j,(i+atmol(imol,1)-1)))).gt.0)) then
              zid2 = i12(j,(i+atmol(imol,1)-1))
              isnew = .true.
              exit
            end if
          end if
        end do
      end if
      if ((zid.le.0).AND.(zid2.le.0)) cycle
!      write(*,*) i,k,i+atmol(imol,1)-1,blst(i,1+k),isnew
!      if (zid.gt.0) write(*,*) 'trying ',zid,bio_code(b_type(zid)),atmres(zid),amino(seqtyp(atmres(zid))),' for rot ',isnew
!      if (zid2.gt.0) write(*,*) 'also trying ',zid2,bio_code(b_type(zid2)),atmres(zid2),amino(seqtyp(atmres(zid2))),' for rot '
      if (zid.gt.0) then
        zbu = ztor(zid)
        ztor(zid) = ztor(zid) + incr
        if (ztor(zid).gt.180.) ztor(zid) = ztor(zid) - 360.
        if (ztor(zid).le.-180.) ztor(zid) = ztor(zid) + 360.
        call makexyz_formol(imol)
        call genints(imol,lvec(:,1),avec(:,1),ivec(:,1),tvec(:,1))
        elig = .true.
        if ((maxval(abs(lvec(:,1)-lvec(:,2))).gt.epsl).OR.(maxval(abs(avec(:,1)-avec(:,2))).gt.epsa).OR.&
 &          (maxval(abs(ivec(:,1)-ivec(:,2))).gt.epsa)) then
          elig = .false.
        end if
        if ((isnew.EQV..false.).AND.(zid2.gt.0).AND.(elig.EQV..true.)) then
          do rs=max(rsmol(imol,1),atmres(blst(i,1+k))-2),min(rsmol(imol,2),atmres(blst(i,1+k))+2)
            if (zid2.eq.fline2(rs)) then
              fline2(rs) = 0
              exit
            else if (zid2.eq.yline2(rs)) then
              yline2(rs) = 0
              exit
            end if
          end do
          if (seqtyp(rs).ne.26) then
            if ((rs.lt.rsmol(molofrs(rs),2)).AND.(seqtyp(min(rs+1,nseq)).eq.26)) then
              write(ilog,*) 'Warning. Existing and rotatable dihedral angle with double setup does not &
 &require second Z-matrix entry for proper rotation. This is likely to be caused by linkage to an unsupported &
 &residue and may change the meaning of phi/psi angles in the affected polypeptide residue.'
              write(ilog,95) bio_code(b_type(i+atmol(imol,1)-1)),i+atmol(imol,1)-1,bio_code(b_type(blst(i,1+k))),blst(i,1+k),&
 &                         bio_code(b_type(zid)),zid,bio_code(b_type(zid2)),zid2
            else
              write(ilog,*) 'Warning. Existing and rotatable dihedral angle with double setup does not &
 &require second Z-matrix entry for proper rotation. This is likely an inconsistency in setting up the Z-matrix or &
 &the result of linkage to a residue of unsupported polymer type.'
              write(ilog,95) bio_code(b_type(i+atmol(imol,1)-1)),i+atmol(imol,1)-1,bio_code(b_type(blst(i,1+k))),blst(i,1+k),&
 &                       bio_code(b_type(zid)),zid,bio_code(b_type(zid2)),zid2
            end if
          end if
          zid2 = 0
        else if ((isnew.EQV..true.).AND.(elig.EQV..true.)) then
          zid2 = 0
        end if

        if (elig.EQV..true.) call gen_rotlst(zid,zid2,imol,tnat,coords(:,:),isnew)
        
!       restore
        x(atmol(imol,1):atmol(imol,2)) = coords(1:tnat,1)
        y(atmol(imol,1):atmol(imol,2)) = coords(1:tnat,2)
        z(atmol(imol,1):atmol(imol,2)) = coords(1:tnat,3)
        ztor(zid) = zbu
        if ((elig.EQV..false.).AND.((zid2.le.0).OR.(isnew.EQV..true.))) zid = 0
      end if
!
      if ((zid2.gt.0).AND.(isnew.EQV..true.)) then
        zbu = ztor(zid2)
        ztor(zid2) = ztor(zid2) + incr
        if (ztor(zid2).gt.180.) ztor(zid2) = ztor(zid2) - 360.
        if (ztor(zid2).le.-180.) ztor(zid2) = ztor(zid2) + 360.
        call makexyz_formol(imol)
        call genints(imol,lvec(:,1),avec(:,1),ivec(:,1),tvec(:,1))
        elig = .true.
        if ((maxval(abs(lvec(:,1)-lvec(:,2))).gt.epsl).OR.(maxval(abs(avec(:,1)-avec(:,2))).gt.epsa).OR.&
 &          (maxval(abs(ivec(:,1)-ivec(:,2))).gt.epsa)) then
          elig = .false.
        end if
        if (elig.EQV..true.) then
          zid = zid2
        end if

        if (elig.EQV..true.) call gen_rotlst(zid,zid2,imol,tnat,coords(:,:),isnew)
        x(atmol(imol,1):atmol(imol,2)) = coords(1:tnat,1)
        y(atmol(imol,1):atmol(imol,2)) = coords(1:tnat,2)
        z(atmol(imol,1):atmol(imol,2)) = coords(1:tnat,3)
        ztor(zid2) = zbu
        if ((elig.EQV..false.).AND.((zid.le.0).OR.(isnew.EQV..true.))) zid2 = 0 
      end if
      if ((zid2.gt.0).AND.(zid.gt.0).AND.(isnew.EQV..false.)) then
        zbu = ztor(zid)
        zbu2 = ztor(zid2)
        ztor(zid) = ztor(zid) + incr
        ztor(zid2) = ztor(zid2) + incr
        if (ztor(zid).gt.180.) ztor(zid) = ztor(zid) - 360.
        if (ztor(zid2).gt.180.) ztor(zid2) = ztor(zid2) - 360.
        if (ztor(zid).le.-180.) ztor(zid) = ztor(zid) + 360.
        if (ztor(zid2).le.-180.) ztor(zid2) = ztor(zid2) + 360.
        call makexyz_formol(imol)
        call genints(imol,lvec(:,1),avec(:,1),ivec(:,1),tvec(:,1))
        elig = .true.
        if ((maxval(abs(lvec(:,1)-lvec(:,2))).gt.epsl).OR.(maxval(abs(avec(:,1)-avec(:,2))).gt.epsa).OR.&
 &          (maxval(abs(ivec(:,1)-ivec(:,2))).gt.epsa)) then
          elig = .false.
        end if
        if (elig.EQV..true.) call gen_rotlst(zid,zid2,imol,tnat,coords(:,:),isnew)
        
        x(atmol(imol,1):atmol(imol,2)) = coords(1:tnat,1)
        y(atmol(imol,1):atmol(imol,2)) = coords(1:tnat,2)
        z(atmol(imol,1):atmol(imol,2)) = coords(1:tnat,3)
        ztor(zid) = zbu
        ztor(zid2) = zbu2
        if (elig.EQV..false.) then
          zid = 0
          zid2 = 0
        end if
      end if
      if ((isnew.EQV..false.).AND.(zid.le.0).AND.(proflag.EQV..false.)) then
        write(ilog,*) 'Fatal. A native dihedral angle appears to not be rotatable freely. This is probably &
 &caused by an unusual linkage between an unsupported and either a supported or an unsupported residue. &
 &The problem may be solved by masking the supported polymer type of unsupported residues, or by renaming &
 &additional residues to be considered as unsupported. It could also indicate a bug in setting up &
 &the Z-matrix.'
        call fexit()
      end if
    end do
  end do
!
  deallocate(tvec)
  deallocate(ivec)
  deallocate(avec)
  deallocate(lvec)
  deallocate(coords)
  deallocate(blst)
!
end 
!
!------------------------------------------------------------------------
!
! compute all internal coordinates in the complete lists (not Z-matrix!)
! this is used to test preservation of covalent geometry 
!
subroutine genints(imol,vl,vb,vi,vt)
!
  use inter
  use molecule
  use sequen
  use fyoc
  use atoms
  use iounit
!
  implicit none
!
  integer imol
  RTYPE vl(max(1,sum(nrsbl(rsmol(imol,1):rsmol(imol,2)))))
  RTYPE vb(max(1,sum(nrsba(rsmol(imol,1):rsmol(imol,2)))))
  RTYPE vi(max(1,sum(nrsimpt(rsmol(imol,1):rsmol(imol,2)))))
  RTYPE vt(max(1,sum(nrsdi(rsmol(imol,1):rsmol(imol,2)))))
  integer rs,cn(4),i,rs1,rs2
  RTYPE getblen,getbang,getztor
!
  cn(:) = 0
  do rs=rsmol(imol,1),rsmol(imol,2)
    if ((disulf(rs).gt.0).AND.(crlk_mode.eq.1)) then ! crosslink treated as restraints has to be skipped
      rs1 = min(rs,disulf(rs))
      rs2 = max(rs,disulf(rs))
      do i=1,nrsbl(rs)
        cn(1) = cn(1) + 1
        if ((maxval(atmres(iaa(rs)%bl(i,1:2))).eq.rs2).AND.(minval(atmres(iaa(rs)%bl(i,1:2))).eq.rs1)) then
          vl(cn(1)) = 0.0
        else
          vl(cn(1)) = getblen(iaa(rs)%bl(i,1),iaa(rs)%bl(i,2))
        end if
      end do
      do i=1,nrsba(rs)
        cn(2) = cn(2) + 1
        if ((maxval(atmres(iaa(rs)%ba(i,1:3))).eq.rs2).AND.(minval(atmres(iaa(rs)%ba(i,1:3))).eq.rs1)) then
          vb(cn(2)) = 0.0
        else
          vb(cn(2)) = getbang(iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),iaa(rs)%ba(i,3))
        end if
      end do
      do i=1,nrsimpt(rs)
        cn(3) = cn(3) + 1
        if ((maxval(atmres(iaa(rs)%impt(i,1:4))).eq.rs2).AND.(minval(atmres(iaa(rs)%impt(i,1:4))).eq.rs1)) then
          vi(cn(3)) = 0.0
        else
          vi(cn(3)) = getztor(iaa(rs)%impt(i,1),iaa(rs)%impt(i,2),iaa(rs)%impt(i,3),iaa(rs)%impt(i,4))
        end if
      end do
      do i=1,nrsdi(rs)
        cn(4) = cn(4) + 1
        if ((maxval(atmres(iaa(rs)%di(i,1:4))).eq.rs2).AND.(minval(atmres(iaa(rs)%di(i,1:4))).eq.rs1)) then
          vt(cn(4)) = 0.0
        else
          vt(cn(4)) = getztor(iaa(rs)%di(i,1),iaa(rs)%di(i,2),iaa(rs)%di(i,3),iaa(rs)%di(i,4))
        end if
      end do
    else if ((crlk_mode.eq.2).AND.(disulf(rs).gt.0)) then
      write(ilog,*) 'Fatal. Hard crosslink constraints are not yet supported in genints(...). This is an &
 &omission bug.'
      call fexit()
    else
      do i=1,nrsbl(rs)
        cn(1) = cn(1) + 1
        vl(cn(1)) = getblen(iaa(rs)%bl(i,1),iaa(rs)%bl(i,2))
      end do
      do i=1,nrsba(rs)
        cn(2) = cn(2) + 1
        vb(cn(2)) = getbang(iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),iaa(rs)%ba(i,3))
      end do
      do i=1,nrsimpt(rs)
        cn(3) = cn(3) + 1
        vi(cn(3)) = getztor(iaa(rs)%impt(i,1),iaa(rs)%impt(i,2),iaa(rs)%impt(i,3),iaa(rs)%impt(i,4))
      end do
      do i=1,nrsdi(rs)
        cn(4) = cn(4) + 1
        vt(cn(4)) = getztor(iaa(rs)%di(i,1),iaa(rs)%di(i,2),iaa(rs)%di(i,3),iaa(rs)%di(i,4))
      end do
    end if
  end do
!
end
!
!----------------------------------------------------------------------
!
! populate rotation list from coordinate set differences
! this assumes the natural building direction (N->C, branches outward)
! the complement set (sans axis atoms) is obtained in straightforward manner
!
subroutine gen_rotlst(z1,z2,imol,tnat,cds,isnew)
!
  use molecule
  use zmatrix
  use atoms
  use iounit
!
  implicit none
!
  integer i,z1,z2,zf,imol,nsgs,tnat
  RTYPE cds(tnat,3),epsx
  integer, ALLOCATABLE:: tmplst(:,:)
  logical inseg,isnew
!
  allocate(tmplst(tnat,2))
  epsx = 1.0e-7
  nsgs = 0
  inseg = .false.
!

  do i=1,tnat

    if ((abs(cds(i,1)-x(i+atmol(imol,1)-1)).gt.epsx).OR.(abs(cds(i,2)-y(i+atmol(imol,1)-1)).gt.epsx).OR.&
        (abs(cds(i,3)-z(i+atmol(imol,1)-1)).gt.epsx)) then

      if (inseg.EQV..false.) then
        nsgs = nsgs + 1
        tmplst(nsgs,1) = i+atmol(imol,1)-1
        inseg = .true.
      end if
    else
      if (inseg.EQV..true.) then
        tmplst(nsgs,2) = i+atmol(imol,1)-2
        inseg = .false.
      end if
    end if
  end do
  if (inseg.EQV..true.) then ! terminate last stretch if needed
    tmplst(nsgs,2) = atmol(imol,2)
  end if
  if ((z1.gt.0).AND.(z2.gt.0)) then
    zf = z1
  else
    zf = z1
    if (z2.gt.0) zf = z2
  end if
!
  if (nsgs.gt.0) then
    allocate(izrot(zf)%rotis(nsgs,2))
    izrot(zf)%alsz = nsgs
    izrot(zf)%rotis(:,:) = tmplst(1:nsgs,:)
  else if (isnew.EQV..false.) then
    print *,z1,z2,zf
    write(ilog,*) 'Fatal. A dihedral angle explicitly defined as a degree of freedom does not rotate any atoms &
 &(Z-matrix line ',zf,'). This is usually the result of introducing colinear geometries, or could be indicative of a bug.'
    call fexit()
  end if
!
  deallocate(tmplst)
!  
end
!
!----------------------------------------------------------------------------
!
subroutine trans_rotlstdof()
!
  use atoms
  use zmatrix
  use movesets
  use fyoc
  use molecule
  use sequen
!
  implicit none
!
  integer i,rs,rsi,j,imol
  logical skip
!
  unklst%nr = 0
  unslst%nr = 0
  natlst%nr = 0
  do i=1,n
    if (izrot(i)%alsz.gt.0) then
      rs = atmres(i)
      imol = molofrs(rs)
      skip = .false.
      do rsi=max(rs-1,rsmol(imol,1)),min(rs+1,rsmol(imol,2))
        if ((i.eq.wline(rsi)).OR.((i.eq.fline(rsi)).AND.(fline2(rsi).le.0)).OR.&
 &          ((i.eq.yline(rsi)).AND.(yline2(rsi).le.0))) then
          natlst%nr = natlst%nr + 1
          skip = .true.
          exit
        else if ((i.eq.fline(rsi)).OR.(i.eq.yline(rsi))) then
          skip=.true.
          exit
        end if
        do j=1,nchi(rsi)
          if (i.eq.chiline(j,rsi)) then
            natlst%nr = natlst%nr + 1
            skip = .true.
            exit
          end if
        end do
        if (skip.EQV..true.) exit
        do j=1,nnucs(rsi)
          if (i.eq.nucsline(j,rsi)) then
            natlst%nr = natlst%nr + 1
            skip = .true.
            exit
          end if
        end do
        if (skip.EQV..true.) exit 
      end do
      if (skip.EQV..false.) then
        if (seqtyp(atmres(i)).eq.26) then
          unklst%nr = unklst%nr + 1
        else
          unslst%nr = unslst%nr + 1
        end if
      end if
    end if
  end do
  if (natlst%nr.gt.0) then
    allocate(natlst%idx(natlst%nr))
    allocate(natlst%wt(natlst%nr))
  end if
  if (unklst%nr.gt.0) then
    allocate(unklst%idx(unklst%nr))
    allocate(unklst%wt(unklst%nr))
  end if
  if (unslst%nr.gt.0) then
    allocate(unslst%idx(unslst%nr))
    allocate(unslst%wt(unslst%nr))
  end if
!
  unklst%nr = 0
  unslst%nr = 0
  natlst%nr = 0
  do i=1,n
    if (izrot(i)%alsz.gt.0) then
      rs = atmres(i)
      imol = molofrs(rs)
      if (i.eq.atmol(imol,1)) cycle
      skip = .false.
      do rsi=max(rs-1,rsmol(imol,1)),min(rs+1,rsmol(imol,2))
        if ((i.eq.wline(rsi)).OR.((i.eq.fline(rsi)).AND.(fline2(rsi).le.0)).OR.&
 &          ((i.eq.yline(rsi)).AND.(yline2(rsi).le.0))) then
          natlst%nr = natlst%nr + 1
          natlst%idx(natlst%nr) = i
          skip = .true.
          exit
        else if ((i.eq.fline(rsi)).OR.(i.eq.yline(rsi))) then
          skip=.true.
          exit
        end if
        do j=1,nchi(rsi)
          if (i.eq.chiline(j,rsi)) then
            natlst%nr = natlst%nr + 1
            natlst%idx(natlst%nr) = i
            skip = .true.
            exit
          end if
        end do
        if (skip.EQV..true.) exit
        do j=1,nnucs(rsi)
          if (i.eq.nucsline(j,rsi)) then
            natlst%nr = natlst%nr + 1
            natlst%idx(natlst%nr) = i
            skip = .true.
            exit
          end if
        end do
        if (skip.EQV..true.) exit 
      end do
      if (skip.EQV..false.) then
        if (seqtyp(atmres(i)).eq.26) then
          unklst%nr = unklst%nr + 1
          unklst%idx(unklst%nr) = i
        else
          unslst%nr = unslst%nr + 1
          unslst%idx(unslst%nr) = i
        end if
      end if
    end if
  end do
!  write(*,*) 'UNS'
!  write(*,*)  unslst%idx(1:unslst%nr)
!  write(*,*) 'NAT'
!  write(*,*)  natlst%idx(1:natlst%nr)
!  write(*,*) 'UNK'
!  write(*,*)  unklst%idx(1:unklst%nr)
!
end
!
!------------------------------------------------------------------------------
!
! this routine parses the rotation set lists for a molecule to identify putative
! overlap relationships between them that are crucial for recursive computation
! of properties
!
subroutine parse_rotlsts(imol)
!
  use molecule
  use zmatrix
  use sequen
  use fyoc
  use atoms
  use interfaces
  use forces
  use iounit
!
  implicit none
!
  integer imol,i,k,l,nelg,kk,ll,rs,tmp(8),kid
  integer, ALLOCATABLE:: elgis(:,:),hlp(:,:),tmpl(:,:)
  logical fndit
!
  k = 0
  do i=atmol(imol,1)+3,atmol(imol,2)
    if (allocated(izrot(i)%rotis).EQV..true.) k = k + 1
  end do
!
  if (k.le.0) then
    izrot(atmol(imol,1))%alsz = 1
    if (allocated(izrot(atmol(imol,1))%rotis).EQV..true.) deallocate(izrot(atmol(imol,1))%rotis)
    allocate(izrot(atmol(imol,1))%rotis(1,2))
    izrot(atmol(imol,1))%rotis(1,1) = atmol(imol,1)
    izrot(atmol(imol,1))%rotis(1,2) = atmol(imol,2)
    if (allocated(izrot(atmol(imol,1))%diffis).EQV..true.) deallocate(izrot(atmol(imol,1))%diffis)
    izrot(atmol(imol,1))%alsz2 = 1
    allocate(izrot(atmol(imol,1))%diffis(izrot(atmol(imol,1))%alsz2,2))
    izrot(atmol(imol,1))%diffis(:,:) = izrot(atmol(imol,1))%rotis(:,:)
    if (allocated(izrot(atmol(imol,1))%treevs).EQV..true.) deallocate(izrot(atmol(imol,1))%treevs)
    allocate(izrot(atmol(imol,1))%treevs(5))
    izrot(atmol(imol,1))%treevs(:) = 0
    izrot(atmol(imol,1))%treevs(2) = 1
    izrot(atmol(imol,1))%treevs(4) = 1
    return
  end if
!
  allocate(elgis(k,9))
  elgis(:,:) = 0
  k = 0
  do i=atmol(imol,2),atmol(imol,1)+3,-1
    if (allocated(izrot(i)%rotis).EQV..true.) then
      k = k + 1
      elgis(k,1) = i
!     set total number of atoms for each group
      do kk=1,izrot(elgis(k,1))%alsz
        elgis(k,2) = elgis(k,2) + izrot(elgis(k,1))%rotis(kk,2) - izrot(elgis(k,1))%rotis(kk,1) + 1
      end do
    end if
  end do
!
! for each group, find the smallest larger group that completely contains it
  nelg = k
  do k=1,nelg
    do l=1,nelg
      if (k.eq.l) cycle
      if (elgis(k,2).gt.elgis(l,2)) cycle ! discard all smaller groups
      if (elgis(k,4).gt.0) then
        if (elgis(l,2).gt.elgis(elgis(k,4),2)) cycle ! discard all groups larger than the one currently found
      end if
      do kk=1,izrot(elgis(k,1))%alsz
        fndit = .false.
        do ll=1,izrot(elgis(l,1))%alsz
          if ((izrot(elgis(l,1))%rotis(ll,2).ge.izrot(elgis(k,1))%rotis(kk,2)).AND.&
 &            (izrot(elgis(l,1))%rotis(ll,1).le.izrot(elgis(k,1))%rotis(kk,1))) then
            fndit = .true.
            exit
          end if
        end do
        if (fndit.EQV..false.) exit
      end do
      if (fndit.EQV..false.) cycle
      elgis(k,4) = l
    end do
!    write(*,*) elgis(k,1),elgis(elgis(k,4),1)
    if (elgis(k,4).gt.0) elgis(elgis(k,4),3) = 1 ! flag the minimum container as a collection point
  end do
!
  l = 0
  do k=1,nelg
    if (elgis(k,3).gt.0) elgis(k,6) = -1
  end do
  do while (1.eq.1)
    fndit = .false.
    do k=1,nelg
      if (elgis(k,6).eq.l) then
        if (elgis(k,4).gt.0) then
          elgis(elgis(k,4),6) = l + 1 ! depth level for access
          fndit = .true.
        end if
      end if
    end do
    l = l + 1
    if (fndit.EQV..false.) exit
  end do
!
! sort according to depth
  do i=1,nelg
    elgis(i,5) = i
  end do
  kk = 1
  ll = nelg
  k = nelg
  fndit = .true.
  call merge_sort(k,fndit,elgis(:,6),elgis(:,7),kk,ll,elgis(:,5),elgis(:,8))
  do i=1,nelg
    elgis(elgis(i,8),5) = i
  end do
  elgis(:,6) = elgis(:,7)
  do i=1,nelg
    if (elgis(i,4).gt.0) elgis(i,4) = elgis(elgis(i,4),5)
  end do
  do i=1,nelg
    if (i.eq.elgis(i,8)) cycle
    tmp(1:4) = elgis(i,1:4)
    elgis(i,1:4) = elgis(elgis(i,8),1:4)
    elgis(elgis(i,8),1:4) = tmp(1:4)
    kk = elgis(i,5)
    ll = elgis(i,8)
    elgis(i,5) = elgis(ll,5)
    elgis(ll,5) = kk
    elgis(kk,8) = elgis(i,8)
  end do
!
! now index branches
  elgis(:,7:9) = 0
  kk = 1
  ll = 1
  do while (elgis(kk,6).eq.0)
    elgis(kk,7) = ll
    ll = ll + 1
    kk = kk + 1
    if (kk.gt.nelg) exit
  end do
  do k=1,nelg
    if (elgis(k,4).eq.0) cycle
    if (elgis(elgis(k,4),7).eq.0) then
      elgis(elgis(k,4),7) = elgis(k,7)
    else
!     this indicates branch merging
      if (elgis(elgis(k,4),7).lt.elgis(k,7)) then
        elgis(k,8) = elgis(elgis(k,4),7) ! the target container already has lower ID, end current branch
      else
        elgis(elgis(k,4),7) = elgis(k,7) ! the target container has higher ID, override and end relevant branches
        do i=1,k-1
          if (elgis(i,4).eq.elgis(k,4)) then
           elgis(i,8) = elgis(k,7)
          end if
        end do
      end if
    end if
  end do
!
! for all mergers, we need to choose the right place in hierarchy to perform merging
! first, populate a counter for how many branches are merging in at a given dihedral
  do k=1,nelg
    if (elgis(k,8).gt.0) then
      elgis(elgis(k,4),9) = elgis(elgis(k,4),9) + 1 ! elgis(k,7)
    end if
  end do  
!
! determine difference sets for rotation lists
  allocate(hlp(nelg,2))
  hlp(:,1) = 0
  do k=1,nelg
    if (elgis(k,4).gt.0) hlp(elgis(k,4),1) = hlp(elgis(k,4),1) + 1
  end do
  ll = maxval(hlp(:,1))
  deallocate(hlp)
  allocate(hlp(nelg,ll+1))
  hlp(:,1) = 0
  do k=1,nelg
    if (elgis(k,4).gt.0) then
      hlp(elgis(k,4),1) = hlp(elgis(k,4),1) + 1
      hlp(elgis(k,4),1+hlp(elgis(k,4),1)) = k
    end if
  end do
  allocate(tmpl(sum(izrot(atmol(imol,1):atmol(imol,2))%alsz),2))
  do k=1,nelg
    if (allocated(izrot(elgis(k,1))%diffis).EQV..true.) deallocate(izrot(elgis(k,1))%diffis)
    if (hlp(k,1).gt.0) then
      tmpl(1:izrot(elgis(k,1))%alsz,:) = izrot(elgis(k,1))%rotis(1:izrot(elgis(k,1))%alsz,:)
      izrot(elgis(k,1))%alsz2 = izrot(elgis(k,1))%alsz
      do rs=2,hlp(k,1)+1
        l = hlp(k,rs)
        do ll=1,izrot(elgis(l,1))%alsz
          kk = 1
          do while (kk.le.izrot(elgis(k,1))%alsz2)
!           note it can never happen that sets need to be joined because this would imply adding atoms
            if ((izrot(elgis(l,1))%rotis(ll,2).lt.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).gt.tmpl(kk,1))) then
              izrot(elgis(k,1))%alsz2 = izrot(elgis(k,1))%alsz2 + 1
              tmpl(izrot(elgis(k,1))%alsz2,1) = izrot(elgis(l,1))%rotis(ll,2) + 1
              tmpl(izrot(elgis(k,1))%alsz2,2) = tmpl(kk,2)
              tmpl(kk,2) = izrot(elgis(l,1))%rotis(ll,1) - 1 ! existing segment shortened
              exit
            else if ((izrot(elgis(l,1))%rotis(ll,2).eq.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).gt.tmpl(kk,1))) then
              tmpl(kk,2) = izrot(elgis(l,1))%rotis(ll,1) - 1 ! existing segment shortened
              exit
            else if ((izrot(elgis(l,1))%rotis(ll,2).lt.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).eq.tmpl(kk,1))) then
              tmpl(kk,1) = izrot(elgis(l,1))%rotis(ll,2) + 1 ! existing segment shortened
              exit
            else if ((izrot(elgis(l,1))%rotis(ll,2).ge.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).le.tmpl(kk,1))) then
              ! delete segment
              tmpl(kk:(izrot(elgis(k,1))%alsz2-1),:) = tmpl((kk+1):izrot(elgis(k,1))%alsz2,:)
              izrot(elgis(k,1))%alsz2 = izrot(elgis(k,1))%alsz2 - 1
            else
              kk = kk + 1
            end if
          end do
          do kk=1,izrot(elgis(k,1))%alsz2
            if ((izrot(elgis(l,1))%rotis(ll,2).gt.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).le.tmpl(kk,2))) then
              tmpl(kk,2) = izrot(elgis(l,1))%rotis(ll,1) - 1 ! existing segment shortened
            else if ((izrot(elgis(l,1))%rotis(ll,2).ge.tmpl(kk,1)).AND.(izrot(elgis(l,1))%rotis(ll,1).lt.tmpl(kk,1))) then
              tmpl(kk,1) = izrot(elgis(l,1))%rotis(ll,2) + 1 ! existing segment shortened
            end if
          end do
        end do
      end do
      if (izrot(elgis(k,1))%alsz2.le.0) then
        write(ilog,*) 'Fatal. In computing difference lists of rotating atoms, list associated with atom ',elgis(k,1),' turns &
 &out empty. This is a bug.'
        call fexit()
      else
        allocate(izrot(elgis(k,1))%diffis(izrot(elgis(k,1))%alsz2,2))
        izrot(elgis(k,1))%diffis(:,:) = tmpl(1:izrot(elgis(k,1))%alsz2,:)
      end if
    else
      izrot(elgis(k,1))%alsz2 = 0 ! this is needed because a torsion is now terminal, but may not have been previously 
    end if
  end do
! sanity check
  deallocate(tmpl)
  deallocate(hlp)
  allocate(tmpl(atmol(imol,2)-atmol(imol,1)+1,2))
  tmpl(:,:) = 0
  do l=1,nelg
    if (elgis(l,6).eq.0) then
      do ll=1,izrot(elgis(l,1))%alsz
        do k=izrot(elgis(l,1))%rotis(ll,1),izrot(elgis(l,1))%rotis(ll,2)
          kk = k - atmol(imol,1) + 1
          if (tmpl(kk,1).eq.0) then
            tmpl(kk,1) = elgis(l,1)
          else
            write(ilog,*) 'Fatal. In setting up difference lists of rotating atoms, atom ',k,' is part of multiple lists.&
 & This is a bug.'
            call fexit()
          end if
        end do
      end do
    else
      do ll=1,izrot(elgis(l,1))%alsz2
        do k=izrot(elgis(l,1))%diffis(ll,1),izrot(elgis(l,1))%diffis(ll,2)
          kk = k - atmol(imol,1) + 1
          if (tmpl(kk,1).eq.0) then
            tmpl(kk,1) = elgis(l,1)
          else
            write(ilog,*) 'Fatal. In setting up difference lists of rotating atoms, atom ',k,' is part of multiple lists.&
 & This is a bug.'
            call fexit()
          end if
        end do
      end do
    end if
  end do
  deallocate(tmpl)
!
! allocate and transfer
  do k=1,nelg
     if (allocated(izrot(elgis(k,1))%treevs).EQV..true.) deallocate(izrot(elgis(k,1))%treevs)
     allocate(izrot(elgis(k,1))%treevs(5+elgis(k,9)))
  end do
  elgis(:,9) = 0
  do k=1,nelg
    if (elgis(k,4).gt.0) then
      izrot(elgis(k,1))%treevs(1) = elgis(elgis(k,4),1) ! minimal container
    else
      izrot(elgis(k,1))%treevs(1) = 0
    end if
    izrot(elgis(k,1))%treevs(2:3) = elgis(k,7:8) ! branch and merger indicator
    izrot(elgis(k,1))%treevs(3) = elgis(k,9)
    izrot(elgis(k,1))%treevs(4) = k ! rank
    izrot(elgis(k,1))%treevs(5) = elgis(k,6) ! depth
    if (elgis(k,8).gt.0) then
      elgis(elgis(k,4),9) = elgis(elgis(k,4),9) + 1
      izrot(elgis(elgis(k,4),1))%treevs(5+elgis(elgis(k,4),9)) = elgis(k,7)
    end if
  end do
!
! set up a rigid-body motion list pinned to the first atom
  izrot(atmol(imol,1))%alsz = 1
  if (allocated(izrot(atmol(imol,1))%rotis).EQV..true.) deallocate(izrot(atmol(imol,1))%rotis)
  allocate(izrot(atmol(imol,1))%rotis(1,2))
  izrot(atmol(imol,1))%rotis(1,1) = atmol(imol,1)
  izrot(atmol(imol,1))%rotis(1,2) = atmol(imol,2)
! difference list also
  allocate(tmpl(atmol(imol,2)-atmol(imol,1)+1,2))
  kid = 1
  tmpl(1,:) =  izrot(atmol(imol,1))%rotis(1,:)
  do l=1,nelg
    do ll=1,izrot(elgis(l,1))%alsz
!     subsets, supersets
      kk = 1
      do while (kk.le.kid)
        if ((izrot(elgis(l,1))%rotis(ll,2).lt.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).gt.tmpl(kk,1))) then
          kid = kid + 1
          tmpl(kid,1) = izrot(elgis(l,1))%rotis(ll,2) + 1
          tmpl(kid,2) = tmpl(kk,2)
          tmpl(kk,2) = izrot(elgis(l,1))%rotis(ll,1) - 1 ! existing segment shortened
          exit
        else if ((izrot(elgis(l,1))%rotis(ll,2).eq.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).gt.tmpl(kk,1))) then
          tmpl(kk,2) = izrot(elgis(l,1))%rotis(ll,1) - 1 ! existing segment shortened
          exit
        else if ((izrot(elgis(l,1))%rotis(ll,2).lt.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).eq.tmpl(kk,1))) then
          tmpl(kk,1) = izrot(elgis(l,1))%rotis(ll,2) + 1 ! existing segment shortened
          exit
        else if ((izrot(elgis(l,1))%rotis(ll,2).ge.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).le.tmpl(kk,1))) then
          ! delete segment
          tmpl(kk:(kid-1),:) = tmpl((kk+1):kid,:)
          kid = kid - 1
        else
          kk = kk + 1
        end if
      end do
!     overlaps
      do kk=1,kid
        if ((izrot(elgis(l,1))%rotis(ll,2).gt.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).le.tmpl(kk,2))) then
          tmpl(kk,2) = izrot(elgis(l,1))%rotis(ll,1) - 1 ! existing segment shortened
        else if ((izrot(elgis(l,1))%rotis(ll,2).ge.tmpl(kk,1)).AND.(izrot(elgis(l,1))%rotis(ll,1).lt.tmpl(kk,1))) then
          tmpl(kk,1) = izrot(elgis(l,1))%rotis(ll,2) + 1 ! existing segment shortened
        end if
      end do
    end do
  end do
  if (allocated(izrot(atmol(imol,1))%diffis).EQV..true.) deallocate(izrot(atmol(imol,1))%diffis)
  izrot(atmol(imol,1))%alsz2 = kid
  allocate(izrot(atmol(imol,1))%diffis(kid,2))
  izrot(atmol(imol,1))%diffis(:,:) = tmpl(1:kid,:)
! finally, set treevs for base atom
  kk = 0
  do k=1,nelg
    if (elgis(k,4).le.0) kk = kk + 1
  end do
  if (allocated(izrot(atmol(imol,1))%treevs).EQV..true.) deallocate(izrot(atmol(imol,1))%treevs)
  allocate(izrot(atmol(imol,1))%treevs(5+kk))
  izrot(atmol(imol,1))%treevs(1) = 0
  izrot(atmol(imol,1))%treevs(2) = 1
  izrot(atmol(imol,1))%treevs(3) = kk
  izrot(atmol(imol,1))%treevs(4) = nelg + 1
  izrot(atmol(imol,1))%treevs(5) = maxval(elgis(1:nelg,6)) + 1
  kk = 0
  do k=1,nelg
    if (elgis(k,4).le.0) then
      kk = kk + 1
      izrot(atmol(imol,1))%treevs(5+kk) = elgis(k,7)
    end if
  end do
!
  deallocate(tmpl)
  deallocate(elgis)
!
end
!
!--------------------------------------------------------------------------------------------------------------
!
! a subroutine to determine the interaction status of two atoms via standard topology and rotation lists (including effective 14)
! this function has to deal with all specifics for both intra-residue and next-residue terms
!
subroutine ia_rotlsts(ii,kk,doia,are14)
!
  use atoms
  use sequen
  use molecule
  use zmatrix
  use inter
  use polypep
  use fyoc
  use system
!
  implicit none
!
  logical doia,are14
  integer j,jj,rs,rs1,rs2,ii,kk,shf,iij,kkj,rb1,rb2
  logical are14bu
!
  doia = .true.
  are14 = .false.
!
  if (molofrs(atmres(ii)).ne.(molofrs(atmres(kk)))) return
  if (abs(atmres(ii)-atmres(kk)).gt.1) return
!
  do j=1,n12(ii)
    if (i12(j,ii).eq.kk) doia = .false.
  end do
  do j=1,n13(ii)
    if (i13(j,ii).eq.kk) doia = .false.
  end do
  do j=1,n14(ii)
    if (use_14.EQV..false.) then
      if (i14(j,ii).eq.kk) doia = .false.
    else
      if (i14(j,ii).eq.kk) are14 = .true.
    end if
  end do
  are14bu = are14
!
! rotation list-based checks
  if (doia.EQV..true.) then
    iij = 0
    kkj = 0
    rs1 = min(atmres(ii),atmres(kk))
    rs2 = max(atmres(ii),atmres(kk))
    rs = rs1
    if (rs.gt.rsmol(molofrs(atmres(ii)),1)) rs = rs-1
    do j=at(rs)%bb(1),at(rs2)%bb(1)+at(rs2)%nbb+at(rs2)%nsc-1
      if (izrot(j)%alsz.le.0) cycle
      if (izrot(j)%treevs(5).eq.0) then
        do jj=1,izrot(j)%alsz
          if ((izrot(j)%rotis(jj,1).le.ii).AND.(izrot(j)%rotis(jj,2).ge.ii)) iij = j
          if ((izrot(j)%rotis(jj,1).le.kk).AND.(izrot(j)%rotis(jj,2).ge.kk)) kkj = j
        end do
      else
        do jj=1,izrot(j)%alsz2
          if ((izrot(j)%diffis(jj,1).le.ii).AND.(izrot(j)%diffis(jj,2).ge.ii)) iij = j
          if ((izrot(j)%diffis(jj,1).le.kk).AND.(izrot(j)%diffis(jj,2).ge.kk)) kkj = j
        end do
      end if
    end do
!
    if (iij.eq.kkj) then
      doia = .false.
    end if
!   need to check for corrections caused by axis atoms and set extended 14
    if (iij.gt.0) then
      if ((izrot(iij)%treevs(1).eq.kkj).AND.(mode_14.eq.2)) are14 = .true.
      if ((iz(1,iij).eq.kk).OR.(iz(2,iij).eq.kk)) doia = .false.
      if (izrot(iij)%treevs(1).gt.0) then
        if (((iz(1,izrot(iij)%treevs(1)).eq.kk).OR.(iz(2,izrot(iij)%treevs(1)).eq.kk)).AND.(mode_14.eq.2)) are14 = .true.
      end if
    end if
    if (kkj.gt.0) then
      if ((izrot(kkj)%treevs(1).eq.iij).AND.(mode_14.eq.2)) are14 = .true.
      if ((iz(1,kkj).eq.ii).OR.(iz(2,kkj).eq.ii)) doia = .false.
      if (izrot(kkj)%treevs(1).gt.0) then
        if (((iz(1,izrot(kkj)%treevs(1)).eq.ii).OR.(iz(2,izrot(kkj)%treevs(1)).eq.ii)).AND.(mode_14.eq.2)) are14 = .true.
      end if
    end if
!
!   we need an override if rotation lists are not appropriate to determine rigidity: this is the case
!   for flexible rings
    if ((seqpolty(rs1).eq.'N').AND.(nucsline(6,rs1).gt.0).AND.(at(rs1)%nsc.ge.3)) then
      rb1 = 0
      rb2 = 0 
      shf = 0
      if (seqflag(rs1).eq.24) shf = 2
      do j=1,n12(ii)
        if (((i12(j,ii).ge.at(rs1)%sc(1)).AND.(i12(j,ii).le.at(rs1)%sc(3))).OR.(i12(j,ii).eq.nuci(rs1,5-shf)).OR.&
 &          (i12(j,ii).eq.nuci(rs1,6-shf))) then
          rb1 = i12(j,ii)
          exit
        end if
      end do
      do j=1,n12(kk)
        if (((i12(j,kk).ge.at(rs1)%sc(1)).AND.(i12(j,kk).le.at(rs1)%sc(3))).OR.(i12(j,kk).eq.nuci(rs1,5-shf)).OR.&
 &          (i12(j,kk).eq.nuci(rs1,6-shf))) then
          rb2 = i12(j,kk)
          exit
        end if
      end do
      if ((rb1.gt.0).AND.(rb2.gt.0)) then
        doia = .true.
        are14 = are14bu
      end if
    else if ((seqpolty(rs1).eq.'P').OR.(seqpolty(rs2).eq.'P')) then
      shf = 0
      if (ua_model.ne.0) shf = 1
      do jj=rs1,rs2
        rb1 = 0
        rb2 = 0
        if ((at(jj)%nsc.lt.(4-shf)).OR.(seqpolty(jj).ne.'P')) cycle
        do j=1,n12(ii)
          if (((i12(j,ii).ge.at(jj)%sc(2-shf)).AND.(i12(j,ii).le.at(jj)%sc(4-shf))).OR.(i12(j,ii).eq.ni(jj)).OR.&
   &          (i12(j,ii).eq.cai(jj))) then
            rb1 = i12(j,ii)
            exit
          end if
        end do
        do j=1,n12(kk)
          if (((i12(j,kk).ge.at(jj)%sc(2-shf)).AND.(i12(j,kk).le.at(jj)%sc(4-shf))).OR.(i12(j,kk).eq.ni(jj)).OR.&
   &          (i12(j,kk).eq.cai(jj))) then
            rb2 = i12(j,kk)
            exit
          end if
        end do
        if ((rb1.gt.0).AND.(rb2.gt.0)) then
          doia = .true.
          are14 = are14bu
          exit
        end if
      end do
    end if      
  end if
!
end
!
!
!-----------------------------------------------------------------------
!
! the following are a few subroutines/fxns to compute individual
! z-matrix terms (lengths, angles, dihedrals) from atomic xyzs
!
!-----------------------------------------------------------------------
!
!
! a function which gives the bond length
! returned as output for two atomic indices (i1,i2)
!
function getblen(i1,i2)
!
  use atoms
!
  implicit none
!
  integer i1,i2
  RTYPE c1(3),getblen
!
  c1(1) = x(i1) - x(i2)
  c1(2) = y(i1) - y(i2)
  c1(3) = z(i1) - z(i2)
  getblen = sqrt(sum(c1(:)**2.0))
!
end
!
!------------------------------------------------------------------
!
! the same operating on alternative coordinates
!
function getblen_ref(i1,i2)
!
  use atoms
!
  implicit none
!
  integer i1,i2
  RTYPE c1(3),getblen_ref
!
  c1(1) = xref(i1) - xref(i2)
  c1(2) = yref(i1) - yref(i2)
  c1(3) = zref(i1) - zref(i2)
  getblen_ref = sqrt(sum(c1(:)**2.0))
!
end
!
!-----------------------------------------------------------------------
!
! a function which gives bond angle mapped onto interval [0,180.0]
! returned as output for three atom indices (i1,i2,i3 in that order)
!
function getbang(i1,i2,i3)
!
  use atoms
  use math
!
  implicit none
!
  integer i1,i2,i3
  RTYPE c1(3),c2(3),c3(3),cc1,cc2,cc3,getbang,dotp,si
!
  c1(1) = x(i1) - x(i2)
  c1(2) = y(i1) - y(i2)
  c1(3) = z(i1) - z(i2)
  c2(1) = x(i3) - x(i2)
  c2(2) = y(i3) - y(i2)
  c2(3) = z(i3) - z(i2)
  cc1 = sqrt(sum(c1(:)**2.0))
  cc2 = sqrt(sum(c2(:)**2.0))
  dotp = sum(c1(:)*c2(:))
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
  c3(:) = c2(:)/cc2 - si*c1(:)/cc1
  cc3 = sqrt(sum(c3(:)**2.0))
!
  getbang = si*2.0*RADIAN*asin(0.5*cc3)
  if (dotp.lt.0.0) then
    getbang = getbang + RADIAN*PI
  end if
!
end
!
!-----------------------------------------------------------------------
!
! a function which gives bond angle mapped onto interval [0,180.0]
! returned as output for three atom indices (i1,i2,i3 in that order)
! differently from getbang, it explicitly checks for superimposed atoms and returns -1.0
! to indicate this condition
!
function getbang_safe(i1,i2,i3)
!
  use atoms
  use math
!
  implicit none
!
  integer i1,i2,i3
  RTYPE c1(3),c2(3),c3(3),cc1,cc2,cc3,getbang_safe,dotp,si
!
  c1(1) = x(i1) - x(i2)
  c1(2) = y(i1) - y(i2)
  c1(3) = z(i1) - z(i2)
  c2(1) = x(i3) - x(i2)
  c2(2) = y(i3) - y(i2)
  c2(3) = z(i3) - z(i2)
  cc1 = sqrt(sum(c1(:)**2.0))
  cc2 = sqrt(sum(c2(:)**2.0))
  if ((cc1.le.0.0).OR.(cc2.le.0.0)) then
    getbang_safe = -1.0
    return
  end if
  dotp = sum(c1(:)*c2(:))
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
  c3(:) = c2(:)/cc2 - si*c1(:)/cc1
  cc3 = sqrt(sum(c3(:)**2.0))
!
  getbang_safe = si*2.0*RADIAN*asin(0.5*cc3)
  if (dotp.lt.0.0) then
    getbang_safe = getbang_safe + RADIAN*PI
  end if
!
end
!
!-----------------------------------------------------------------------
!
! the same operating on alternative coordinates
!
function getbang_ref(i1,i2,i3)
!
  use atoms
  use math
!
  implicit none
!
  integer i1,i2,i3
  RTYPE c1(3),c2(3),c3(3),cc1,cc2,cc3,getbang_ref,dotp,si
!
  c1(1) = xref(i1) - xref(i2)
  c1(2) = yref(i1) - yref(i2)
  c1(3) = zref(i1) - zref(i2)
  c2(1) = xref(i3) - xref(i2)
  c2(2) = yref(i3) - yref(i2)
  c2(3) = zref(i3) - zref(i2)
  cc1 = sqrt(sum(c1(:)**2.0))
  cc2 = sqrt(sum(c2(:)**2.0))
  dotp = sum(c1(:)*c2(:))
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
  c3(:) = c2(:)/cc2 - si*c1(:)/cc1
  cc3 = sqrt(sum(c3(:)**2.0))
!
  getbang_ref = si*2.0*RADIAN*asin(0.5*cc3)
  if (dotp.lt.0.0) then
    getbang_ref = getbang_ref + RADIAN*PI
  end if
!
end
!
!-----------------------------------------------------------------------
!
! a function which gives torsional angle in degrees [-180:180] returned
! as output for four atom indices (i1,i2,i3,i4 in that order)
!
function getztor(i1,i2,i3,i4)
!
  use atoms
  use math
!
  implicit none
!
  integer i1,i2,i3,i4
  RTYPE getztor,c1(3),c2(3),c3(3),c4(3)
  RTYPE cp1(3),np1,cp2(3),np2,np4,si,dotp
!
  c1(1) = x(i2) - x(i1)
  c1(2) = y(i2) - y(i1)
  c1(3) = z(i2) - z(i1)
  c2(1) = x(i3) - x(i2)
  c2(2) = y(i3) - y(i2)
  c2(3) = z(i3) - z(i2)
  c3(1) = x(i4) - x(i3)
  c3(2) = y(i4) - y(i3)
  c3(3) = z(i4) - z(i3)
  call crossprod(c1,c2,cp1,np1)
  call crossprod(c2,c3,cp2,np2)
!
  if (np1*np2 .gt. 0.0d0) then
    dotp = sum(cp1(:)*cp2(:))
    if (dotp.lt.0.0) then
      si = -1.0
    else
      si = 1.0
    end if
    c4(:) = cp2(:)/np2 - si*cp1(:)/np1
    np4 = sqrt(sum(c4(:)**2.0))
    getztor = si*2.0*RADIAN*asin(0.5*np4)
    if (dotp.lt.0.0) then
      getztor = getztor + RADIAN*PI
    end if
! 
    si = sum(c1(:)*cp2(:))
    if (si.lt.0.0) getztor = -getztor
  else
!   torsion is ill-defined, left at 0.0 (not necessarily fatal, though)
    getztor = 0.0
  end if
!
end
!
!-------------------------------------------------------------------------
!
! a simple wrapper to substitute for a (missing) pointer array
!
function getpuckertor(rs,zline)
!
  use polypep
  use sequen
  use molecule
  use iounit
  use aminos
  use fyoc
!
  implicit none
!
  integer rs,zline,shf
  RTYPE getztor,getpuckertor
!
  if (seqpolty(rs).eq.'N') then
    shf = 0
    if (((seqtyp(rs).ge.76).AND.(seqtyp(rs).le.87)).OR.((seqtyp(rs).eq.26).AND.(nuci(rs,5).le.0))) shf = 2
    getpuckertor = getztor(nuci(rs,4-shf),nuci(rs,5-shf),nuci(rs,6-shf),nucsline(6,rs))
    zline = nucsline(6,rs)
  else if ((seqtyp(rs).eq.9).OR.(seqtyp(rs).eq.25).OR.(seqtyp(rs).eq.32)) then
    if (rs.eq.rsmol(molofrs(rs),1)) then
      getpuckertor = getztor(ci(rs),cai(rs),ni(rs),fline(rs))
      zline = fline(rs)
    else
      getpuckertor = getztor(ci(rs),cai(rs),ni(rs),ci(rs-1))
      zline = ci(rs)
    end if
  else
     write(ilog,*) 'Fatal. Called getpuckertor(...) with unsupported residue.'
     call fexit()
  end if
!
  return
!
end
!
!-----------------------------------------------------------------------
!
! gives bond angle in radian
! angle is always mapped onto interval [0,PI/2] (smaller inside angle only)
!
subroutine bondang(p1,p2,p3,ang)
!
  implicit none
!
  integer k
  RTYPE p1(3),p2(3),p3(3),c1(3),c2(3),c3(3),cc1,cc2,cc3,dotp,ang,si
!
  cc1 = 0.0
  cc2 = 0.0
  do k=1,3
    c1(k) = p1(k) - p2(k)
    cc1 = cc1 + c1(k)**2
    c2(k) = p2(k) - p3(k)
    cc2 = cc2 + c2(k)**2
  end do
  cc1 = sqrt(cc1)
  cc2 = sqrt(cc2)
  dotp = sum(c1(:)*c2(:))
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
  c3(:) = c2(:)/cc2 - si*c1(:)/cc1
  cc3 = sqrt(sum(c3(:)**2.0))
!
  ang = si*2.0*asin(0.5*cc3)
  if (dotp.lt.0.0) then
    ang = -ang
  end if
!
end
!
!------------------------------------------------------------
!
! gives bond angle in radian, assumes atoms p1,p2,p3 connected
! in this order
! angle is always mapped onto interval [0,PI] (inside angles only)
!
subroutine bondang2(p1,p2,p3,ang)
!
  use math
!
  implicit none
!
  integer k
  RTYPE p1(3),p2(3),p3(3),c1(3),c2(3),c3(3),cc1,cc2,cc3,si,dotp,ang
!
  cc1 = 0.0
  cc2 = 0.0
  do k=1,3
    c1(k) = p1(k) - p2(k)
    cc1 = cc1 + c1(k)**2
    c2(k) = p3(k) - p2(k)
    cc2 = cc2 + c2(k)**2
  end do
  cc1 = sqrt(cc1)
  cc2 = sqrt(cc2)
  dotp = sum(c1(:)*c2(:))
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
  c3(:) = c2(:)/cc2 - si*c1(:)/cc1
  cc3 = sqrt(sum(c3(:)**2.0))
!
  ang = si*2.0*asin(0.5*cc3)
  if (dotp.lt.0.0) then
    ang = ang + PI
  end if
!
end
!
!-----------------------------------------------------------------------
!
! a similar subroutine which also gives the derivatives with respect to p1,p2,p3
!
subroutine bondang2_dcart(p1,p2,p3,ang,fvec1,fvec2,fvec3)
!
  use iounit
  use polypep
  use atoms
  use math
!
  implicit none
!
  RTYPE dv12(3),id12,p1(3),p2(3),p3(3),ang
  RTYPE dv32(3),id32,dotp
  RTYPE vecs(3),nvecs,frang
  RTYPE fvec1(3),fvec3(3),fvec2(3),si
!
  dv12(:) = p1(:) - p2(:)
  dv32(:) = p3(:) - p2(:)
  id12 = 1.0/sqrt(dv12(1)**2 + dv12(2)**2 + dv12(3)**2)
  id32 = 1.0/sqrt(dv32(1)**2 + dv32(2)**2 + dv32(3)**2)
  dotp = dv12(1)*dv32(1) + dv12(2)*dv32(2) + &
 &          dv12(3)*dv32(3)
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
! this half-angle formula has the advantage of being safe to use in all circumstances
! as 0.5*nvecs never gets near the -1/1 limits
! in addition there is considerable precision loss with the (simpler) acos-form
  vecs(1) = id32*dv32(1) - si*id12*dv12(1)
  vecs(2) = id32*dv32(2) - si*id12*dv12(2)
  vecs(3) = id32*dv32(3) - si*id12*dv12(3)
  nvecs = sqrt(vecs(1)**2 + vecs(2)**2 + vecs(3)**2)
  ang = si*2.0*asin(0.5*nvecs)
  frang = si*2.0/sqrt(1.0 - (0.5*nvecs)**2)
  if (dotp.lt.0.0) then
    ang = ang + PI
  end if
! there is a price, however, primarily in the derivatives
  fvec1(1) = 0.5*frang*(1.0/nvecs)*(vecs(1)*&
 &  (-si*id12 + si*dv12(1)*(id12**3)*dv12(1))&
 & + vecs(2)*(si*dv12(2)*(id12**3)*dv12(1))&
 & + vecs(3)*(si*dv12(3)*(id12**3)*dv12(1)) )
  fvec1(2) = 0.5*frang*(1.0/nvecs)*(vecs(2)*&
 &  (-si*id12 + si*dv12(2)*(id12**3)*dv12(2))&
 & + vecs(1)*(si*dv12(1)*(id12**3)*dv12(2))&
 & + vecs(3)*(si*dv12(3)*(id12**3)*dv12(2)) )
  fvec1(3) = 0.5*frang*(1.0/nvecs)*(vecs(3)*&
 &  (-si*id12 + si*dv12(3)*(id12**3)*dv12(3))&
 & + vecs(1)*(si*dv12(1)*(id12**3)*dv12(3))&
 & + vecs(2)*(si*dv12(2)*(id12**3)*dv12(3)) )
  fvec3(1) = 0.5*frang*(1.0/nvecs)*(vecs(1)*&
 &  (id32 - dv32(1)*(id32**3)*dv32(1))&
 & + vecs(2)*(-dv32(2)*(id32**3)*dv32(1))&
 & + vecs(3)*(-dv32(3)*(id32**3)*dv32(1)) )
  fvec3(2) = 0.5*frang*(1.0/nvecs)*(vecs(2)*&
 &  (id32 - dv32(2)*(id32**3)*dv32(2))&
 & + vecs(1)*(-dv32(1)*(id32**3)*dv32(2))&
 & + vecs(3)*(-dv32(3)*(id32**3)*dv32(2)) )
  fvec3(3) = 0.5*frang*(1.0/nvecs)*(vecs(3)*&
 &  (id32 - dv32(3)*(id32**3)*dv32(3))&
 & + vecs(1)*(-dv32(1)*(id32**3)*dv32(3))&
 & + vecs(2)*(-dv32(2)*(id32**3)*dv32(3)) )
  fvec2(:) = - (fvec1(:) + fvec3(:))
!
end
!
!------------------------------------------------------------
!
! gives angle between bond two vectors
! angle is always mapped onto interval [0,PI] (inside angle only)
!
subroutine bondang3(c1,c2,ang)
!
  use math
!
  implicit none
!
  RTYPE c1(3),c2(3),c3(3),cc1,cc2,cc3,dotp,ang,si
!
  cc1 = sqrt(sum(c1(:)**2.0))
  cc2 = sqrt(sum(c2(:)**2.0))
  dotp = sum(c1(:)*c2(:))
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
  c3(:) = c2(:)/cc2 - si*c1(:)/cc1
  cc3 = sqrt(sum(c3(:)**2.0))
!
  ang = si*2.0*asin(0.5*cc3)
  if (dotp.lt.0.0) then
    ang = ang + PI
  end if
!
end
!
!----------------------------------------------------------------------------
!
! gives torsion in radian, note sign convention
!
subroutine dihed(p1,p2,p3,p4,tor)
!
  use iounit
  use math
  use zmatrix
!
  implicit none
!
  integer k
  RTYPE p1(3),p2(3),p3(3),p4(3),c1(3),c2(3),c3(3),c4(3),np4
  RTYPE pl_1(3),pl_2(3),np1,np2,si,sign_conv,tor,dotp
!
  tor = 0.0
!
  do k=1,3
    c1(k) = p2(k) - p1(k)
    c2(k) = p3(k) - p2(k)
    c3(k) = p4(k) - p3(k)
  end do
  call crossprod(c1,c2,pl_1,np1)
  call crossprod(c2,c3,pl_2,np2)
!
  if (np1*np2 .gt. 0.0d0) then
    dotp = sum(pl_1(:)*pl_2(:))
    if (dotp.lt.0.0) then
      si = -1.0
    else
      si = 1.0
    end if
    c4(:) = pl_2(:)/np2 - si*pl_1(:)/np1
    np4 = sqrt(sum(c4(:)**2.0))
    tor = si*2.0*asin(0.5*np4)
    if (dotp.lt.0.0) then
      tor = tor + PI
    end if
! 
    sign_conv = c1(1)*pl_2(1) + c1(2)*pl_2(2) + c1(3)*pl_2(3)
    if (sign_conv .lt. 0.0d0)  tor = -tor
  else
!   torsion is ill-defined, left at 0.0 (not necessarily fatal, though)
    dihed_wrncnt = dihed_wrncnt + 1
    if (dihed_wrncnt.eq.dihed_wrnlmt) then
      write(ilog,*) 'Warning: Call to dihed(...) with bad vectors.'
      write(ilog,*) 'This is warning number #',dihed_wrncnt,' of this type not all of&
 & which may be displayed.'
      if (10.0*dihed_wrnlmt.gt.0.5*HUGE(dihed_wrnlmt)) then
        dihed_wrncnt = 0 ! reset
      else
        dihed_wrnlmt = dihed_wrnlmt*10
      end if
    end if
  end if
!
end
!
!-----------------------------------------------------------------------
!
! a similar subroutine which lso gives derivatives with respect to
! p1-p4
!
subroutine dihed_dcart(p1,p2,p3,p4,tor,dtmd1,dtmd2,dtmd3,dtmd4)
!
  use iounit
  use polypep
  use atoms
  use math
!
  implicit none
!
  integer j
  RTYPE dv12(3),dv23(3),p1(3),p2(3),p3(3),p4(3)
  RTYPE dv43(3),ndv23,incp123,incp234,ncp123,ncp234
  RTYPE dotp,dotpfg,cp123(3),cp234(3),dtmd1(3),dtmd2(3)
  RTYPE dtmd3(3),dotphg,dtmd4(3),dotpah,si
  RTYPE tor,c4(3),np4
!
  dv12(:) = p1(:) - p2(:) 
  dv23(:) = p2(:) - p3(:)
  dv43(:) = p4(:) - p3(:)

  cp123(1) = dv12(2)*dv23(3) - dv12(3)*dv23(2)
  cp123(2) = dv12(3)*dv23(1) - dv12(1)*dv23(3)
  cp123(3) = dv12(1)*dv23(2) - dv12(2)*dv23(1)
  cp234(1) = dv43(2)*dv23(3) - dv43(3)*dv23(2)
  cp234(2) = dv43(3)*dv23(1) - dv43(1)*dv23(3)
  cp234(3) = dv43(1)*dv23(2) - dv43(2)*dv23(1)
! norms and dot-product
  ncp123 = cp123(1)**2 + cp123(2)**2 + cp123(3)**2
  ncp234 = cp234(1)**2 + cp234(2)**2 + cp234(3)**2
  dotp = cp123(1)*cp234(1) + cp123(2)*cp234(2) + &
 &          cp123(3)*cp234(3)
  dotpfg = dv12(1)*dv23(1) + dv12(2)*dv23(2) +&
 &           dv12(3)*dv23(3)
  dotphg = dv43(1)*dv23(1) + dv43(2)*dv23(2) +&
 &           dv43(3)*dv23(3)
  dotpah = cp123(1)*dv43(1) + cp123(2)*dv43(2) +&
 &           cp123(3)*dv43(3)
  ndv23 = sqrt(dv23(1)**2 + dv23(2)**2 + dv23(3)**2)
!
! filter for exceptions which would otherwise NaN
! colinear reference atoms (note that this also detects the equally fatal |dv23| = 0)
  if ((ncp123.le.0.0).OR.(ncp234.le.0.0)) then
!   this will set all derivatives to zero and the (now arbitrary) dihedral to zero deg.
    write(ilog,*) 'WARNING. Colinear reference atoms in computation of &
 &torsional derivative (dihed_dcart(...)). This indicates an unstabl&
 &e simulation or a bug. Expect run to crash soon.'
    incp123 = 0.0
    incp234 = 0.0
  else
    incp123 = 1.0/ncp123
    incp234 = 1.0/ncp234
  end if
!
! equations 27 in Blondel and Karplus
  do j=1,3
    dtmd1(j) = -cp123(j)*incp123*ndv23

    dtmd2(j) = cp123(j)*incp123*&
 &                  (ndv23 + dotpfg/ndv23) - &
 &               cp234(j)*dotphg*incp234/ndv23

    dtmd3(j) = cp234(j)*incp234*&
 &                  (dotphg/ndv23 - ndv23) - &
 &               cp123(j)*dotpfg*incp123/ndv23
    dtmd4(j) = cp234(j)*incp234*ndv23
  end do
!
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
  c4(:) = cp234(:)*sqrt(incp234) - si*cp123(:)*sqrt(incp123)
  np4 = sqrt(sum(c4(:)**2.0))
  tor = si*2.0*asin(0.5*np4)
  if (dotp.lt.0.0) then
    tor = tor + PI
  end if
! 
  si = sum(-dv12(:)*cp234(:))
  if (si.lt.0.0) tor = -tor
!
end
!
!-----------------------------------------------------------------------
!
! a routine to check for colinear reference atoms for torsional tuples
!
function check_colinear(i1,i2,i3,i4)
!
  implicit none
!
  integer i1,i2,i3,i4
  logical check_colinear
  RTYPE ang,getbang
!
  check_colinear = .false.
!
  ang = getbang(i1,i2,i3)
  if ((abs(ang).le.0.01).OR.(abs(180.0-ang).le.0.01)) then
    check_colinear = .true.
    return
  end if
!
  ang = getbang(i2,i3,i4)
  if ((abs(ang).le.0.01).OR.(abs(180.0-ang).le.0.01)) check_colinear = .true.
!
  return
!
end
!
!--------------------------------------------------------------------------
!
