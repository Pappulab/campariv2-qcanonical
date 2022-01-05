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
!-----------------------------------------------------------------
!
! fill some fmsmc-specific atom arrays
!
subroutine setup_pka_prm() ! Martin Added
  use atoms
  use params
  use sequen
  use polypep
  use energies
  use aminos
  use math
  use system
  use fyoc 
  
  use movesets ! For hsq
    
  integer i,j

  atr(1:n) = lj_rad(attyp(1:n)) ! Martin : changed that, it should work 

  if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.))then
    
    do i=1,n
        atr_limits(i,1)=lj_rad(bio_ljtyp(b_type(i)))! Martin : I could go back to the more compact form here 
        atr_limits(i,2)=lj_rad(bio_ljtyp(transform_table(b_type(i))))
        atr_limits(i,3)=lj_rad(bio_ljtyp(transform_table(his_eqv_table(b_type(i)))))
    end do 
    atvol_limits(:,1)=(atr_limits(:,1)**3)*(4./3.)*PI
    atvol_limits(:,2)=(atr_limits(:,2)**3)*(4./3.)*PI
    atvol_limits(:,3)=(atr_limits(:,3)**3)*(4./3.)*PI
    atbvol_limits(:,1) = (4./3.)*PI*(atr_limits(:,1)+par_IMPSOLV(1))**3 - atvol_limits(:,1)
    atbvol_limits(:,2) = (4./3.)*PI*(atr_limits(:,2)+par_IMPSOLV(1))**3 - atvol_limits(:,2)
    atbvol_limits(:,3) = (4./3.)*PI*(atr_limits(:,3)+par_IMPSOLV(1))**3 - atvol_limits(:,3)
! Update :
! So I need to keep the zeros for the atr, but replace them later so they don't change value.
! Martin : Here I need a smater check, because the overide of radius is still happening 
    
    do i=1,n_biotyp ! home : why is that again 
        ! Martin : Goes through the atoms i, if find one that is equal to zero, replaces it by the new value
        ! The reason is that we needed to keep the zero for accurate setting of the limits of atr
        ! But the sigma can not be equal to 0
        ! MArtin : added test to avoid zero index error 
        if (bio_ljtyp(i).ne.0) then 
            if (((abs(lj_sig(bio_ljtyp(transform_table(i)),bio_ljtyp(transform_table(i)))).lt.1.0D-5).or.&
            &(abs(lj_sig(bio_ljtyp(i),bio_ljtyp(i))).lt.1.0D-5)).and.&
            &((abs(lj_sig(bio_ljtyp(transform_table(i)),bio_ljtyp(transform_table(i)))).lt.1.0D-5).neqv.&
            &(abs(lj_sig(bio_ljtyp(i),bio_ljtyp(i))).lt.1.0D-5))) then 

                lj_sig(bio_ljtyp(i),bio_ljtyp(i))=DMAX1(lj_sig(bio_ljtyp(transform_table(i)),bio_ljtyp(transform_table(i))),&  
                &lj_sig(bio_ljtyp(i),bio_ljtyp(i)))
                lj_sig(bio_ljtyp(transform_table(i)),bio_ljtyp(transform_table(i)))=lj_sig(bio_ljtyp(i),bio_ljtyp(i))
                ! Martin if you find one, replace by the other value, then execute the clean readik command again to clean the
!                ! combination rules
            else if (((i.ne.1216).and.(i.ne.1223).and.(i.ne.1208)).and.&
            &((abs(lj_sig(bio_ljtyp(transform_table(i)),bio_ljtyp(transform_table(i)))).lt.1.0D-5).and.&
            &(abs(lj_sig(bio_ljtyp(i),bio_ljtyp(i))).lt.1.0D-5))) then
                print *,"!!! ERROR : cannot not have an atom going from sigma=0 to sigma=0"
                print *,lj_sig(bio_ljtyp(transform_table(i)),bio_ljtyp(transform_table(i))),lj_sig(bio_ljtyp(i),bio_ljtyp(i))
                print *,i,transform_table(i)
                call fexit()
            end if
        end if 
    end do 

    do i=1, n_ljtyp
        do j=1, n_ljtyp
            if (lj_sig(i,j).ne.lj_sig(lj_transform_table(i),lj_transform_table(j))) then ! This indicates that this contact should
            ! have a special rule
                lj_sig(i,j)=lj_sig(lj_transform_table(i),lj_transform_table(j))
                lj_sig_14(i,j)=lj_sig_14(lj_transform_table(i),lj_transform_table(j))
                
            end if 
            if (lj_eps(i,j).ne.lj_eps(lj_transform_table(i),lj_transform_table(j))) then ! This indicates that this contact should
            ! have a special rule
                lj_eps(i,j)=lj_eps(lj_transform_table(i),lj_transform_table(j))
                lj_eps_14(i,j)=lj_eps_14(lj_transform_table(i),lj_transform_table(j))
                
            end if 
        end do 
    end do
    
    call clean_radik()
    call clean_epsik()
    
    do i=1,n_biotyp
        do j=1,n_biotyp
            if ((bio_ljtyp(i).ne.0).and.(bio_ljtyp(j).ne.0)) then ! Martin : I don't test for transformed,
            !becauseI know that it can't happen
                lj_eps_limits(bio_ljtyp(i),bio_ljtyp(j),3)=lj_eps(bio_ljtyp(transform_table(his_eqv_table(i))),&
                &bio_ljtyp(transform_table(his_eqv_table(j))))
                lj_eps_limits(bio_ljtyp(i),bio_ljtyp(j),2)=lj_eps(bio_ljtyp(transform_table(i)),bio_ljtyp(transform_table(j)))
                lj_eps_limits(bio_ljtyp(i),bio_ljtyp(j),1)=lj_eps(bio_ljtyp(i),bio_ljtyp(j))
                !transform_table(his_eqv_table(i)) could be replace by jsut his_eqv_table(i)
                lj_eps_14_limits(bio_ljtyp(i),bio_ljtyp(j),3)=lj_eps_14(bio_ljtyp(transform_table(his_eqv_table(i))),&
                &bio_ljtyp(transform_table(his_eqv_table(j))))
                lj_eps_14_limits(bio_ljtyp(i),bio_ljtyp(j),2)=lj_eps_14(bio_ljtyp(transform_table(i)),bio_ljtyp(transform_table(j)))
                lj_eps_14_limits(bio_ljtyp(i),bio_ljtyp(j),1)=lj_eps_14(bio_ljtyp(i),bio_ljtyp(j))
                
                lj_sig_limits(bio_ljtyp(i),bio_ljtyp(j),3)=lj_sig(bio_ljtyp(transform_table(his_eqv_table(i))),&
                &bio_ljtyp(transform_table(his_eqv_table(j))))
                lj_sig_limits(bio_ljtyp(i),bio_ljtyp(j),2)=lj_sig(bio_ljtyp(transform_table(i)),bio_ljtyp(transform_table(j)))
                lj_sig_limits(bio_ljtyp(i),bio_ljtyp(j),1)=lj_sig(bio_ljtyp(i),bio_ljtyp(j))
                
                lj_sig_14_limits(bio_ljtyp(i),bio_ljtyp(j),3)=lj_sig_14(bio_ljtyp(transform_table(his_eqv_table(i))),&
                &bio_ljtyp(transform_table(his_eqv_table(j))))
                lj_sig_14_limits(bio_ljtyp(i),bio_ljtyp(j),2)=lj_sig_14(bio_ljtyp(transform_table(i)),bio_ljtyp(transform_table(j)))
                lj_sig_14_limits(bio_ljtyp(i),bio_ljtyp(j),1)=lj_sig_14(bio_ljtyp(i),bio_ljtyp(j))
            end if 
        end do 
    end do 
end if
end subroutine setup_pka_prm


subroutine absinth_atom()
!
  use atoms
  use params
  use sequen
  use polypep
  use energies
  use aminos
  use zmatrix
  use math
  use system
  use iounit
  use fyoc
  
  
  use movesets ! for hsq
!
  implicit none
!
  integer i,j,k,rs,ati,shf,shf2,shf3,atwo,afour, lim
  logical sayyes
  RTYPE dum,dis,d1,d2
  logical saidit,saidit2

  atwo = 2
  afour = 4
  sayyes = .true.
  saidit = .false.
 55 format('Here, atom ',i7,' of type ',i4,' has radius ',g10.3,'.')
 56 format('Warning. Atom ',i7,' of biotype ',i5,' and LJ-type ',i5,' has ill-defined mass. Further warnings omitted.')

 if ((do_pka.eqv..false.).and.(do_hsq.eqv..false.)) then 
     atr(1:n) = lj_rad(attyp(1:n))
 end if 
 
 
! potential patch overrides for radii
  call read_atpatchfile(afour)
!
! now checks
  do i=1,n
    if ((atr(i).le.0.0).AND.((use_IMPSOLV.EQV..true.).OR.(savcalc.le.nsim)).AND.(do_pka_2.eqv..false.).and.&
    &(do_hsq.eqv..false.)) then
        write(ilog,*) 'Fatal. Atomic radii have to be finite for ABSINTH implicit solvent model&
 & or solvent-accessible volume analysis to work.'
      write(ilog,55) i,attyp(i),atr(i)
      write(ilog,*) 'Use the RADIUS keyword in the parameter file or a patch (FMCSC_RPATCHFILE) to solve this problem.'
      call fexit()
    end if
    atvol(i) = (4./3.)*PI*atr(i)**3
    atbvol(i) = (4./3.)*PI*(atr(i)+par_IMPSOLV(1))**3 - atvol(i)
    
    if ((mass(i).le.0.0).AND.(saidit.EQV..false.)) then
      saidit = .true.
      write(ilog,56) i,b_type(i),bio_ljtyp(b_type(i))
    end if
  end do
  
  
  do i=1,n
    atsavred(i) = atvol(i)
  end do
  
  do i=1,n
    k = iz(1,i)
!   disconnected atoms have to be excluded
!   note that even though the first atom in a multi-atom molecule
!   will always have k=0, the connectivity will still be handled
!   through the atom(s) actually connected to the first one
    if (k.le.0) cycle ! Martin :
    dis = sqrt((x(k) - x(i))**2&
 &           + (y(k) - y(i))**2&
 &           + (z(k) - z(i))**2)

    if (dis.lt.(atr(i)+atr(k))) then ! MArtin : If the spheres touch
      if (atr(i).ge.atr(k)) then ! Martin : if i is the largest atom
        if (dis.lt.(atr(i)-atr(k))) then !Martin : If the canter of the smaller sphere is inside the bigger one 
          atsavred(i) = atsavred(i) - 0.5*atvol(k)! Martin : Only problem is that this is not really correct....
          atsavred(k) = atsavred(k) - 0.5*atvol(k)! Rohit says does not matter, so it does not matter
        else
          atsavred(i) = atsavred(i) - 0.5*atvol(k)*&
 &                       ((atr(i)+atr(k)-dis)/(2.0*atr(k)))
          atsavred(k) = atsavred(k) - 0.5*atvol(k)*&
 &                       ((atr(i)+atr(k)-dis)/(2.0*atr(k)))
        end if
      else
        if (dis.lt.(atr(k)-atr(i))) then
          atsavred(i) = atsavred(i) - 0.5*atvol(i)
          atsavred(k) = atsavred(k) - 0.5*atvol(i)

        else
          atsavred(i) = atsavred(i) - 0.5*atvol(i)*&
 &                       ((atr(k)+atr(i)-dis)/(2.0*atr(i)))
          atsavred(k) = atsavred(k) - 0.5*atvol(i)*&
 &                       ((atr(k)+atr(i)-dis)/(2.0*atr(i)))
        end if
      end if
    end if        
  end do

! rings have additional bonds not showing up through iz - recover through nadd/iadd
  shf = 0
  shf2 = 0
  shf3 = 0

  if (ua_model.gt.0) then
    shf = 1
    shf2 = 2
    shf3 = 3
  end if
  
  do rs=1,nadd
    i = iadd(1,rs)
    k = iadd(2,rs)
    dis = sqrt((x(k) - x(i))**2&
 &           + (y(k) - y(i))**2&
 &           + (z(k) - z(i))**2)
!   distance override for crosslinks may be necessary
    do j=1,n_crosslinks
      if ((crosslink(j)%itstype.eq.1).OR.(crosslink(j)%itstype.eq.2)) then
        if (((i.eq.at(crosslink(j)%rsnrs(1))%sc(3-shf)).AND.(k.eq.at(crosslink(j)%rsnrs(2))%sc(3-shf))).OR.&
 &          ((k.eq.at(crosslink(j)%rsnrs(1))%sc(3-shf)).AND.(i.eq.at(crosslink(j)%rsnrs(2))%sc(3-shf)))) then
          dis = 2.03 ! make consistent with elsewhere
!              -> note we cannot infer the distance at this point since structure does not satisfy crosslink yet
          exit
        end if
      else
        write(ilog,*) 'Fatal. Encountered unsupported crosslink type in absinth_atom().&
 & This is most likely an omission bug.'
        call fexit()
      end if
    end do
    
    if (dis.lt.(atr(i)+atr(k))) then ! Martin : if atoms touch 
      if (atr(i).ge.atr(k)) then
        if (dis.lt.(atr(i)-atr(k))) then
          atsavred(i) = atsavred(i) - 0.5*atvol(k)
          atsavred(k) = atsavred(k) - 0.5*atvol(k)
        else
          atsavred(i) = atsavred(i) - 0.5*atvol(k)*&
 &                       ((atr(i)+atr(k)-dis)/(2.0*atr(k)))
          atsavred(k) = atsavred(k) - 0.5*atvol(k)*&
 &                       ((atr(i)+atr(k)-dis)/(2.0*atr(k)))
        end if
      else
        if (dis.lt.(atr(k)-atr(i))) then
          atsavred(i) = atsavred(i) - 0.5*atvol(i)
          atsavred(k) = atsavred(k) - 0.5*atvol(i)
        else
          atsavred(i) = atsavred(i) - 0.5*atvol(i)*&
 &                       ((atr(k)+atr(i)-dis)/(2.0*atr(i)))
          atsavred(k) = atsavred(k) - 0.5*atvol(i)*&
 &                       ((atr(k)+atr(i)-dis)/(2.0*atr(i)))
        end if
      end if
   end if
  end do
!
  do i=1,n
    if (atvol(i).gt.0.0) then
      atsavred(i) = atsavred(i)/atvol(i)
    else
      atsavred(i) = 0.0 
    end if
  end do

  call read_atpatchfile(atwo)

  call setup_savol(sayyes,dum)
  !
  if ((use_IMPSOLV.EQV..true.).OR.(savcalc.le.nsim)) then
    saidit = .false.
    saidit2 = .false.
    do ati=1,n
      if (atsavmaxfr(ati).gt.1.0) then
         write(ilog,*) 'Fatal. Setup of solvation parameters is incorrect for&
 & atom ',ati,' (type: ',b_type(ati),'). This is a bug.'
         call fexit()
      else if (atsavmaxfr(ati).le.0.0) then
         if (saidit.EQV..false.) then
           saidit = .true.
           write(ilog,*) 'Fatal. Maximum solvent-accessible volume fraction for&
 & atom ',ati,' (type: ',b_type(ati),') is zero. This could be a bug or indicates &
 &radii and bond lengths that are incompatible with underlying assumptions. Check &
 &FMCSC_SAVPATCHFILE, FMCSC_ASRPATCHFILE, and FMCSC_RPATCHFILE for workarounds or use FMCSC_UNSAFE to skip this error.&
 & Further warnings omitted.'
         end if
         
         if (be_unsafe.EQV..false.) call fexit()
      end if
      if ((atsavred(ati).le.0.0).AND.(do_pka_2.eqv..false.).and.(do_hsq.eqv..false.)) then
         if (saidit2.EQV..false.) then
           saidit2 = .true.
           write(ilog,*) 'Fatal. Atomic volume reduction factor for atom ',ati,' (type: ',b_type(ati),') &
 &is zero or negative. This indicates assumed radii that are incompatible with underlying assumptions. Check &
 &FMCSC_ASRPATCHFILE and FMCSC_RPATCHFILE for workarounds or use FMCSC_UNSAFE to skip this error. Further warnings omitted.'
         end if
         if (be_unsafe.EQV..false.) call fexit()
      end if
      atsavprm(ati,1) = par_IMPSOLV(6)*atsavmaxfr(ati) + &
 &                    (1.0-par_IMPSOLV(6))*par_IMPSOLV(5)
      d1 = 1.0/(1.0 + exp(-(par_IMPSOLV(5)-atsavprm(ati,1))/&
 &par_IMPSOLV(3)))
      d2 = 1.0/(1.0 + exp(-(atsavmaxfr(ati)-atsavprm(ati,1))/&
 &par_IMPSOLV(3)))
      atsavprm(ati,2) = 1.0/(d2-d1)
      atsavprm(ati,3) = 1.0 - atsavprm(ati,2)*(d2-0.5)
      atsavprm(ati,4) = par_IMPSOLV(7)*atsavmaxfr(ati) + &
 &                    (1.0-par_IMPSOLV(7))*par_IMPSOLV(5)
      d1 = 1.0/(1.0 + exp(-(par_IMPSOLV(5)-atsavprm(ati,4))/&
 &par_IMPSOLV(4)))
      d2 = 1.0/(1.0 + exp(-(atsavmaxfr(ati)-atsavprm(ati,4))/&
 &par_IMPSOLV(4)))
      atsavprm(ati,5) = 1.0/(d2-d1)
      atsavprm(ati,6) = 1.0 - atsavprm(ati,5)*(d2-0.5)
    end do
  end if
end
!
!-------------------------------------------------------------------
!
subroutine absinth_savprm()
!
  use energies
  use atoms
!
  implicit none
!
  integer ati
  RTYPE d1,d2
!
  do ati=1,n
    atsavprm(ati,1) = par_IMPSOLV(6)*atsavmaxfr(ati) + &
 &                    (1.0-par_IMPSOLV(6))*par_IMPSOLV(5)
    d1 = 1.0/(1.0 + exp(-(par_IMPSOLV(5)-atsavprm(ati,1))/&
 &par_IMPSOLV(3)))
    d2 = 1.0/(1.0 + exp(-(atsavmaxfr(ati)-atsavprm(ati,1))/&
 &par_IMPSOLV(3)))
    atsavprm(ati,2) = 1.0/(d2-d1)
    atsavprm(ati,3) = 1.0 - atsavprm(ati,2)*(d2-0.5)
    atsavprm(ati,4) = par_IMPSOLV(7)*atsavmaxfr(ati) + &
 &                    (1.0-par_IMPSOLV(7))*par_IMPSOLV(5)
    d1 = 1.0/(1.0 + exp(-(par_IMPSOLV(5)-atsavprm(ati,4))/&
 &par_IMPSOLV(4)))
    d2 = 1.0/(1.0 + exp(-(atsavmaxfr(ati)-atsavprm(ati,4))/&
 &par_IMPSOLV(4)))
    atsavprm(ati,5) = 1.0/(d2-d1)
    atsavprm(ati,6) = 1.0 - atsavprm(ati,5)*(d2-0.5)
  end do
!
end
!
!-----------------------------------------------------------------------------------
!
! sorry, this function probably is not expected to be here,
! but i couldn't find a better place
! this one sets up molecular volume parameters from atomic ones
! has to be called much later than absinth_atom, though
!
! quantities: molvol(mt,1): net vdW-volume occupied by molecules of type mt
!             molvol(mt,2): net SAV (not overlap-, only topology-corrected)
!             resvol(rs,1): net vdW-volume occupied by residue rs
!             resvol(rs,2): net heavy atom SAV (not overlap- only topology-corrected) for res. rs
!             resvol(rs,3): number of heavy atoms in residue rs
!             sysvol(3)   : net vdW-volume
! 
subroutine absinth_molecule()
!
  use iounit
  use energies
  use molecule
  use atoms
  use polypep
  use sequen
  use math
  use aminos
  use torsn
!
  implicit none
!
  integer mt,i,rs,ii
  RTYPE d1,getblen
!
  longestmol = 1
  do mt=1,nmoltyp
    do i=atmol(moltyp(mt,1),1),atmol(moltyp(mt,1),2)
      molvol(mt,1) = molvol(mt,1) + atsavred(i)*atvol(i)
      molvol(mt,2) = molvol(mt,2) + atsavmaxfr(i)*atbvol(i)
      molmass(mt) = molmass(mt) + mass(i)
    end do
    do rs=rsmol(moltyp(mt,1),1),rsmol(moltyp(mt,1),2)
      if (seqtyp(rs).ne.26) then
        molcontlen(mt) = molcontlen(mt) + rescontlen(seqtyp(rs))
      else if (seqpolty(rs).eq.'P') then
        molcontlen(mt) = molcontlen(mt) + rescontlen(1)
      else if ((seqpolty(rs).eq.'N').AND.((nuci(rs,6).gt.0).OR.(rs.gt.rsmol(moltyp(mt,1),1)))) then
        molcontlen(mt) = molcontlen(mt) + rescontlen(64)
      else if (seqpolty(rs).eq.'N') then
        molcontlen(mt) = molcontlen(mt) + rescontlen(76)
      else if ((rs.gt.rsmol(moltyp(mt,1),1)).AND.(rs.lt.rsmol(moltyp(mt,1),2))) then
        molcontlen(mt) = molcontlen(mt) + getblen(at(rs)%bb(1),at(rs+1)%bb(1))
      else
        molcontlen(mt) = molcontlen(mt) + resrad(rs)
      end if
    end do
    if ((rsmol(moltyp(mt,1),2)-rsmol(moltyp(mt,1),1)+1).gt.longestmol) then
      longestmol = rsmol(moltyp(mt,1),2)-rsmol(moltyp(mt,1),1)+1
    end if
  end do
  maxseglen = longestmol
  sysvol(3) = 0.0
!
  do rs=1,nseq
    resvol(rs,:) = 0.0
    do i=1,at(rs)%nbb+at(rs)%nsc
      if (i.le.at(rs)%nbb) then
        ii = at(rs)%bb(i)
      else
        ii = at(rs)%sc(i-at(rs)%nbb)
      end if
      resvol(rs,1) = resvol(rs,1) + atsavred(ii)*atvol(ii)
      sysvol(1) = sysvol(1) + atsavred(ii)*atvol(ii)
      if (mass(ii).gt.5.0) then
        resvol(rs,2) = resvol(rs,2) + atbvol(ii)
        resvol(rs,3) = resvol(rs,3) + 1.0
      end if
    end do
    sysvol(3) = sysvol(3) + resvol(rs,2)
  end do
!
  d1 = ((3./(4.*PI))*(sysvol(1)))**(1./3.)
  sysvol(2) = 1.3504*(4./3.)*PI*((d1 + par_IMPSOLV(1))**3 - d1**3)
!  d1 =  sysvol(1)/(1.5*1.5*PI)
!  sysvol(3) = 1.5*PI*(d1+par_IMPSOLV(1))*((1.5+par_IMPSOLV(1))**2.0)
!
end
!
!----------------------------------------------------------------------------
!
! set a couple of reference arrays for residues
! ditto for expected location of fxn
!
subroutine absinth_residue()
!
  use iounit
  use polypep
  use sequen
  use atoms
  use molecule
  use aminos
!
  implicit none
!
  integer rs
  character(3) resname
!
  do rs=1,nseq
!
    at(rs)%na = at(rs)%nbb + at(rs)%nsc
    resname = amino(seqtyp(rs))
    if (resname.ne.'UNK') seqpolty(rs) = aminopolty(seqtyp(rs))
!
!   first set the reference atom for each residue
    if ((resname.eq.'NA+').OR.(resname.eq.'CL-').OR.&
 &      (resname.eq.'K+ ').OR.(resname.eq.'BR-').OR.&
 &      (resname.eq.'CS+').OR.(resname.eq.'I- ').OR.&
 &      (resname.eq.'NH4').OR.(resname.eq.'AC-').OR.&
 &      (resname.eq.'T3P').OR.(resname.eq.'SPC').OR.&
 &      (resname.eq.'URE').OR.(resname.eq.'FOR').OR.&
 &      (resname.eq.'NH2').OR.(resname.eq.'ACA').OR.&
 &      (resname.eq.'PPA').OR.(resname.eq.'NMA').OR.&
 &      (resname.eq.'CH4').OR.(resname.eq.'FOA').OR.&
 &      (resname.eq.'T4P').OR.(resname.eq.'DMA').OR.&
 &      (resname.eq.'MOH').OR.(resname.eq.'GDN').OR.&
 &      (resname.eq.'EOH').OR.(resname.eq.'MSH').OR.&
 &      (resname.eq.'1MN').OR.(resname.eq.'LCP').OR.&
 &      (resname.eq.'NO3').OR.(resname.eq.'T5P').OR.&
 &      (resname.eq.'T4E').OR.(resname.eq.'BEN').OR.&
 &      (resname.eq.'NAP').OR.(resname.eq.'O2 ').OR.&
 &      (resname.eq.'NAX').OR.(resname.eq.'CLX')) then
       refat(rs) = at(rs)%bb(1)
    else if (resname.eq.'NMF') then
       refat(rs) = at(rs)%bb(3)
    else if ((resname.eq.'PCR').OR.(resname.eq.'TOL').OR.&
 &           (resname.eq.'IMD').OR.(resname.eq.'IME').OR.&
 &           (resname.eq.'PRP').OR.(resname.eq.'IBU').OR.&
 &           (resname.eq.'EMT').OR.(resname.eq.'NBU').OR.&
 &           (resname.eq.'2MN')) then
       refat(rs) = at(rs)%bb(2)
    else if (resname.eq.'THY') then
       refat(rs) = at(rs)%bb(8)
    else if ((resname.eq.'URA').OR.(resname.eq.'MIN')) then
       refat(rs) = at(rs)%bb(4)
    else if (resname.eq.'CYT') then
       refat(rs) = at(rs)%bb(7)
    else if (resname.eq.'ADE') then
       refat(rs) = at(rs)%bb(4)
    else if ((resname.eq.'GUA').OR.(resname.eq.'PUR')) then
       refat(rs) = at(rs)%bb(2)
    else if ((resname.eq.'R5P').OR.(resname.eq.'D5P')) then
!      the C5*
       refat(rs) = nuci(rs,4)
    else if ((resname.eq.'RIB').OR.(resname.eq.'DIB')) then
!      the C4*
       refat(rs) = nuci(rs,3)
    else if ((resname.eq.'RPU').OR.(resname.eq.'RPC').OR.&
 &           (resname.eq.'RPT').OR.(resname.eq.'RPA').OR.&
 &           (resname.eq.'RPG').OR.(resname.eq.'DPU').OR.&
 &           (resname.eq.'DPC').OR.(resname.eq.'DPT').OR.&
 &           (resname.eq.'DPA').OR.(resname.eq.'DPG').OR.&
 &           (resname.eq.'RIU').OR.(resname.eq.'RIC').OR.&
 &           (resname.eq.'RIT').OR.(resname.eq.'RIA').OR.&
 &           (resname.eq.'RIG').OR.(resname.eq.'DIU').OR.&
 &           (resname.eq.'DIC').OR.(resname.eq.'DIT').OR.&
 &           (resname.eq.'DIA').OR.(resname.eq.'DIG')) then
       refat(rs) = at(rs)%sc(2)
    else if ((resname.eq.'GLY').OR.(resname.eq.'PRO').OR.&
 &      (resname.eq.'ALA').OR.(resname.eq.'ABA').OR.&
 &      (resname.eq.'NVA').OR.(resname.eq.'VAL').OR.&
 &      (resname.eq.'LEU').OR.(resname.eq.'ILE').OR.&
 &      (resname.eq.'MET').OR.(resname.eq.'PHE').OR.&
 &      (resname.eq.'NLE').OR.(resname.eq.'SER').OR.&
 &      (resname.eq.'THR').OR.(resname.eq.'CYS').OR.&
 &      (resname.eq.'HIE').OR.(resname.eq.'HID').OR.&
 &      (resname.eq.'TYR').OR.(resname.eq.'TRP').OR.&
 &      (resname.eq.'ASP').OR.(resname.eq.'GLU').OR.&
 &      (resname.eq.'LYS').OR.(resname.eq.'ARG').OR.&
 &      (resname.eq.'GLN').OR.(resname.eq.'ASN').OR.&
 &      (resname.eq.'HYP').OR.(resname.eq.'ACE').OR.&
 &      (resname.eq.'NME').OR.(resname.eq.'AIB').OR.&
 &      (resname.eq.'DAB').OR.(resname.eq.'ORN').OR.&
 &      (resname.eq.'GAM').OR.(resname.eq.'PCA').OR.&
 &      (resname.eq.'GLH').OR.(resname.eq.'ASH').OR.&
 &      (resname.eq.'TYO').OR.(resname.eq.'CYX').OR.&
 &      (resname.eq.'LYD').OR.(resname.eq.'Y1P').OR.&
 &      (resname.eq.'S1P').OR.(resname.eq.'T1P').OR.&
 &      (resname.eq.'HIP').OR.(resname.eq.'GLX').OR.&
 &      (resname.eq.'LYX').OR.(resname.eq.'ASX').OR.&
 &      (resname.eq.'HEX').OR.(resname.eq.'HDX').OR.&
 &      (resname.eq.'Y2P').OR.(resname.eq.'S2P').OR.&
 &      (resname.eq.'T2P').OR.(resname.eq.'SXP').OR.&
 &      (resname.eq.'TXP').OR.(resname.eq.'YXP').OR.&
 &      (resname.eq.'TYX').OR.(resname.eq.'HIX')) then
      refat(rs) = cai(rs)
    else if (resname.eq.'UNK') then
      cycle
    else
      write(ilog,*) 'Fatal. Reference atom undefined for residue ',&
 &resname,'. Check back later.'
      call fexit()
    end if
!
    if (refat(rs).le.0) then
      write(ilog,*) 'Fatal. Reference atom is not set for residue ',&
 &rs,' (',resname,'). This is an omission bug.'
      call fexit()
    end if
!
!   seqflag is a way to replace checks with simpler ones
!   other types are: 12 other N-cap for peptides, 13 other C-cap for peptides;
!                    26 other 5'-cap for nucleotides; 28 3'-cap for nucleotides; 50 - unknown polymer; 
    if (((seqtyp(rs).ge.1).AND.(seqtyp(rs).le.23).AND.(seqtyp(rs).ne.9)).OR.((seqtyp(rs).ge.33).AND.(seqtyp(rs).le.35)).OR.&
 &      (seqtyp(rs).eq.51).OR.(seqtyp(rs).eq.31).OR.((seqtyp(rs).ge.108).AND.(seqtyp(rs).le.116)).OR.&
 &      (seqtyp(rs).ge.119)) then
      seqflag(rs) = 2
    else if ((seqtyp(rs).eq.9).OR.(seqtyp(rs).eq.25).OR.(seqtyp(rs).eq.32)) then
      seqflag(rs) = 5
    else if (seqtyp(rs).eq.24) then
      seqflag(rs) = 8
    else if ((seqtyp(rs).eq.27).OR.(seqtyp(rs).eq.28)) then
      seqflag(rs) = 10 ! polypeptide N-cap supporting omega in rs+1 only
    else if (seqtyp(rs).eq.29) then
      seqflag(rs) = 11 ! polypeptide C-cap supporting omega on itself only
    else if (seqtyp(rs).eq.30) then
      seqflag(rs) = 13 ! polypeptide C-cap not supporting omega
    else if ((seqtyp(rs).ge.76).AND.(seqtyp(rs).le.87)) then 
      seqflag(rs) = 24 ! nucleoside residue with or without sugar sampling (5'-terminal, the latter via nucsline(6,rs))
    else if ((seqtyp(rs).ge.64).AND.(seqtyp(rs).le.75)) then 
      seqflag(rs) = 22 ! nucleotide residue with or without sugar sampling
    else if ((seqtyp(rs).eq.36).OR.(seqtyp(rs).eq.37).OR.(seqtyp(rs).eq.39).OR.(seqtyp(rs).eq.40).OR.&
 &           (seqtyp(rs).eq.45).OR.(seqtyp(rs).eq.46).OR.((seqtyp(rs).ge.52).AND.(seqtyp(rs).le.56)).OR.&
 &           ((seqtyp(rs).ge.100).AND.(seqtyp(rs).le.103)).OR.(seqtyp(rs).eq.117).OR.(seqtyp(rs).eq.118)) then
      seqflag(rs) = 101 ! mandatory small molecules with no possibility or intent for any torsional degrees of freedom
    else if ((seqtyp(rs).eq.60).OR.((seqtyp(rs).ge.104).AND.(seqtyp(rs).le.107))) then
      seqflag(rs) = 102 ! mandatory small molecules that are topologically completely rigidified in a torsional sense
    else if ((seqtyp(rs).eq.38).OR.(seqtyp(rs).eq.57).OR.(seqtyp(rs).eq.58).OR.(seqtyp(rs).eq.59).OR.&
 &           ((seqtyp(rs).ge.61).AND.(seqtyp(rs).le.63)).OR.((seqtyp(rs).ge.41).AND.(seqtyp(rs).le.44)).OR.&
 &           ((seqtyp(rs).ge.47).AND.(seqtyp(rs).le.50)).OR.((seqtyp(rs).ge.88).AND.(seqtyp(rs).le.99))) then
      seqflag(rs) = 103 ! mandatory small molecules with torsional degrees of freedom (whether active by default or not)
    end if
  end do
!
end
!
!--------------------------------------------------------------------
!
! a subroutine to assign parameters for bonded interactions from
! the parsed parameter file
!

subroutine assign_bndtprms()
!
  use params
  use iounit
  use polypep
  use sequen
  use aminos
  use inter
  use atoms
  use math
  use system
  use energies
  ! MArtin : temporary 
  use zmatrix
  ! end martin
!
  implicit none
!
  integer rs,i,j,jj,k,kk,i1,i2,i3,i4,k1,k2,k3,k4,i5,k5,t1,t2
  integer, ALLOCATABLE:: dum(:,:),dum2(:,:)
  RTYPE, ALLOCATABLE:: pams(:,:)
  RTYPE testimp,getztor
  logical missing,reverse,allmissing(5),check_colinear
!
 777 format('BLXX ',a4,a4,a4,a4,' ',6(f10.5,1x))
 778  format('BLBA ',a4,a4,a4,' ',2(f10.5,1x))
 76   format('Matching bond length potential (#',i4,') for atoms ',&
 &i7,' and ',i7,' in residue # ',i6)
 86   format('Guessed bond length potential for atoms ',&
 &i7,' and ',i7,' in residue # ',i6)
 88   format('Cannot guess bond length potential for atoms ',&
 &i7,' and ',i7,' in residue # ',i6)
 87   format('Type: ',i3,' | Length: ',f7.2,' | Parameter 1: ',f7.2)
 77   format('Missing bond length potential specification for atoms ',&
 &i7,' and ',i7,' in residue # ',i6)
 78   format('--> Biotypes ',i4,' (',a4,') and ',i4,' (',a4,') in residu&
 &e type ',a3,'.')
 67   format('Missing bond angle potential specification for atoms ',&
 &i7,',',i7,', and ',i7,' in residue # ',i6)
 66   format('Matching bond angle potential (#',i4,') for atoms ',&
 &i7,',',i7,', and ',i7,' in residue # ',i6)
 96   format('Guessed bond angle potential for atoms ',&
 &i7,',',i7,', and ',i7,' in residue # ',i6)
 98   format('Cannot guess bond angle potential for atoms ',&
 &i7,',',i7,', and ',i7,' in residue # ',i6)
 97   format('Type: ',i3,' | Angle : ',f7.2,' | Parameter 1: ',f7.2)
 68   format('--> Biotypes ',i4,' (',a4,'), ',i4,' (',a4,'), and ',i4,' &
 &(',a4,') in residue type ',a3,'.')
 56   format('Matching torsional potential (#',i4,') for atoms ',&
 &i7,',',i7,',',i7,', and ',i7,' in residue # ',i6)
 57   format('No matching torsional potential specification for atoms ',&
 &i7,',',i7,',',i7,', and ',i7,' in residue # ',i6)
 58   format('--> Biotypes ',i4,' (',a4,'), ',i4,' (',a4,'), ',i4,&
 &' (',a4,'), and ',i4,' (',a4,') in residue type ',a3,'.')
 36   format('Matching CMAP potential (#',i4,') for atoms ',&
 &i7,',',i7,',',i7,',',i7,', and ',i7,' in residue # ',i6)
 37   format('No matching CMAP potential specification for atoms ',&
 &i7,',',i7,',',i7,',',i7,', and ',i7,' with start residue # ',i6)
 38   format('--> Biotypes ',i4,' (',a4,'), ',i4,' (',a4,'), ',i4,&
 &' (',a4,'), ',i4,' (',a4,'), and ',i4,' (',a4,') -> residue type of middle atom: ',a3,'.')
 46   format('Matching improper dihedral potential (#',i4,') for &
 &atoms ',i7,',',i7,',',i7,', and ',i7,' in residue # ',i6)
 47   format('No matching improper dihedral potential specification for &
 &atoms ',i7,',',i7,',',i7,', and ',i7,' in residue # ',i6)
 48   format('--> Biotypes ',i4,' (',a4,'), ',i4,' (',a4,'), ',i4,&
 &' (',a4,'), and ',i4,' (',a4,') in residue type ',a3,'.')
 27   format('Proper or improper dihedral potentials disabled for &
 &atoms ',i7,',',i7,',',i7,', and ',i7,' in residue # ',i6,' due to colinear reference atoms.')
!    
  call strlims(paramfile,t1,t2)
  allmissing(:) = .true.
  if ((bonded_report.EQV..true.).AND.(use_BOND(1).EQV..true.)) then
    write(ilog,*) 
    write(ilog,*) '--------    Bond Length Terms    --------'
  end if
  do rs=1,nseq
    if (use_BOND(1).EQV..false.) exit
!   bond length terms (mandatory)
    do i=1,nrsbl(rs)
      missing = .true.
      i1 = bio_botyp(b_type(iaa(rs)%bl(i,1)))
      i2 = bio_botyp(b_type(iaa(rs)%bl(i,2)))
      do k=1,bo_lstsz

        k1 = bo_lst(k,1)
        k2 = bo_lst(k,2)
        if (((i1.eq.k1).AND.(i2.eq.k2)).OR.& ! Martin : goes through all the bond type to find the two right ones 
 &          ((i2.eq.k1).AND.(i1.eq.k2))) then 
          allmissing(1) = .false.
          missing = .false.
          kk = bo_lst(k,3)
          
          iaa(rs)%typ_bl(i) = bo_typ(kk) ! Martin : here is the assignment 
          if (bonded_report.EQV..true.) then
            write(ilog,76) kk,iaa(rs)%bl(i,1),iaa(rs)%bl(i,2),rs
            write(ilog,78) b_type(iaa(rs)%bl(i,1)),&
 &bio_code(b_type(iaa(rs)%bl(i,1))),b_type(iaa(rs)%bl(i,2)),&
 &bio_code(b_type(iaa(rs)%bl(i,2))),amino(seqtyp(rs))
          end if
          do j=1,MAXBOPAR
            iaa(rs)%par_bl(i,j) = bo_par(kk,j)
          end do
          if (iaa(rs)%typ_bl(i).eq.3) then ! GROMOS quartic: store squared equ length
            iaa(rs)%par_bl(i,3) = iaa(rs)%par_bl(i,2)*iaa(rs)%par_bl(i,2)
          end if
        end if
      end do
      if (missing.EQV..true.) then
        
        print *,bio_botyp(b_type(iaa(rs)%bl(i,1))),bio_botyp(b_type(iaa(rs)%bl(i,2)))
        if (iaa(rs)%typ_bl(i).le.0) then
          write(ilog,77) iaa(rs)%bl(i,1),iaa(rs)%bl(i,2),rs
          write(ilog,78) b_type(iaa(rs)%bl(i,1)),&
 &bio_code(b_type(iaa(rs)%bl(i,1))),b_type(iaa(rs)%bl(i,2)),&
 &bio_code(b_type(iaa(rs)%bl(i,2))),amino(seqtyp(rs))
        else
          allmissing(1) = .false.
          write(ilog,86) iaa(rs)%bl(i,1),iaa(rs)%bl(i,2),rs
          write(ilog,78) b_type(iaa(rs)%bl(i,1)),&
 &bio_code(b_type(iaa(rs)%bl(i,1))),b_type(iaa(rs)%bl(i,2)),&
 &bio_code(b_type(iaa(rs)%bl(i,2))),amino(seqtyp(rs))
          write(ilog,87) iaa(rs)%typ_bl(i),iaa(rs)%par_bl(i,2),iaa(rs)%par_bl(i,1)
        end if
      end if
    end do
  end do
  
  
  
  
 
  if ((allmissing(1).EQV..true.).AND.(use_BOND(1).EQV..true.).AND.&
 &    (fycxyz.ne.2)) then
    write(ilog,*) 'No applicable assignments for bond length potentials found in parameter file ',paramfile(t1:t2),'.'
  end if
!
! bond angle terms (mandatory)
  if ((bonded_report.EQV..true.).AND.(use_BOND(2).EQV..true.)) then
    write(ilog,*)
    write(ilog,*) '--------    Bond Angle Terms     --------'
  end if

  do rs=1,nseq
    if (use_BOND(2).EQV..false.) exit
    do i=1,nrsba(rs)
      missing = .true.
      i1 = bio_botyp(b_type(iaa(rs)%ba(i,1)))
      i2 = bio_botyp(b_type(iaa(rs)%ba(i,2)))
      i3 = bio_botyp(b_type(iaa(rs)%ba(i,3)))
      
      do k=1,ba_lstsz
        k1 = ba_lst(k,1)
        k2 = ba_lst(k,2)
        k3 = ba_lst(k,3)
        if (i2.ne.k2) cycle
        if (((i1.eq.k1).AND.(i3.eq.k3)).OR.&
 &          ((i3.eq.k1).AND.(i1.eq.k3))) then
          missing = .false.
          allmissing(2) = .false.
          kk = ba_lst(k,4)
          iaa(rs)%typ_ba(i) = ba_typ(kk)
          if (bonded_report.EQV..true.) then
               if (((do_pka.eqv..true.).or.(do_hsq.eqv..true.)).and.((atr(iaa(rs)%ba(i,1)).ne.0).and.(atr(iaa(rs)%ba(i,2)).ne.0)&
               &.and.(atr(iaa(rs)%ba(i,3)).ne.0))) then 
                  write(ilog,66) kk,iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),&
         &iaa(rs)%ba(i,3),rs
                  write(ilog,68) b_type(iaa(rs)%ba(i,1)),&
         &bio_code(b_type(iaa(rs)%ba(i,1))),b_type(iaa(rs)%ba(i,2)),&
         &bio_code(b_type(iaa(rs)%ba(i,2))),b_type(iaa(rs)%ba(i,3)),&
         &bio_code(b_type(iaa(rs)%ba(i,3))),amino(seqtyp(rs))
               end if 
          end if
          do j=1,MAXBAPAR
            iaa(rs)%par_ba(i,j) = ba_par(kk,j)
          end do
          if ((ba_typ(kk).ge.1).AND.(ba_typ(kk).le.2)) then ! note that unit conversion is not necessary for GROMOS cos-harmonic
            iaa(rs)%par_ba(i,1)=iaa(rs)%par_ba(i,1)/(RADIAN*RADIAN)
          end if
          if (iaa(rs)%typ_ba(i).eq.3) then ! GROMOS cos-harmonic: store cosine of equ angle
            iaa(rs)%par_ba(i,3) = cos(iaa(rs)%par_ba(i,2)/RADIAN)
          end if
        end if
      end do
      if (missing.EQV..true.) then
        if (iaa(rs)%typ_ba(i).le.0) then
          write(ilog,67) iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),&
 &iaa(rs)%ba(i,3),rs
          write(ilog,68) b_type(iaa(rs)%ba(i,1)),&
 &bio_code(b_type(iaa(rs)%ba(i,1))),b_type(iaa(rs)%ba(i,2)),&
 &bio_code(b_type(iaa(rs)%ba(i,2))),b_type(iaa(rs)%ba(i,3)),&
 &bio_code(b_type(iaa(rs)%ba(i,3))),amino(seqtyp(rs))
        else
          allmissing(2) = .false.
          write(ilog,96) iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),iaa(rs)%ba(i,3),rs
          write(ilog,68) b_type(iaa(rs)%ba(i,1)),&
 &bio_code(b_type(iaa(rs)%ba(i,1))),b_type(iaa(rs)%ba(i,2)),&
 &bio_code(b_type(iaa(rs)%ba(i,2))),b_type(iaa(rs)%ba(i,3)),&
 &bio_code(b_type(iaa(rs)%ba(i,3))),amino(seqtyp(rs))
          write(ilog,97) iaa(rs)%typ_ba(i),iaa(rs)%par_ba(i,2),iaa(rs)%par_ba(i,1)
        end if
      end if
    end do
  end do
  if ((allmissing(2).EQV..true.).AND.(use_BOND(2).EQV..true.).AND.&
 &    (fycxyz.ne.2)) then
    write(ilog,*) 'No applicable assignments for bond angle potentials found in parameter file ',paramfile(t1:t2),'.'
  end if
!
! improper dihedral terms (optional)
  if ((bonded_report.EQV..true.).AND.(use_BOND(3).EQV..true.)) then
    write(ilog,*)
    write(ilog,*) '-------- Improper Dihedral Terms --------'
  end if
  do rs=1,nseq
    if (use_BOND(3).EQV..false.) exit
    do i=1,nrsimpt(rs),3
      missing = .true.
      if (check_colinear(iaa(rs)%impt(i,improper_conv(1)),iaa(rs)%impt(i,2),iaa(rs)%impt(i,improper_conv(2)),&
 &         iaa(rs)%impt(i,4)).EQV..true.) then
        write(ilog,27) iaa(rs)%impt(i,1:4),rs
        cycle
      end if
      i1 = bio_botyp(b_type(iaa(rs)%impt(i,1)))
      i2 = bio_botyp(b_type(iaa(rs)%impt(i,2)))
      i3 = bio_botyp(b_type(iaa(rs)%impt(i,3)))
      i4 = bio_botyp(b_type(iaa(rs)%impt(i,4)))
!     all three same (triple degeneracy)
      if ((i2.eq.i3).AND.(i2.eq.i4)) then
        jj = -1
        do k=1,impt_lstsz
          k1 = impt_lst(k,1)
          k2 = impt_lst(k,2)
          k3 = impt_lst(k,3)
          k4 = impt_lst(k,4)
          if (.NOT.((k4.eq.k3).AND.(k4.eq.k2).AND.(k4.eq.i2))) cycle
          jj = jj + 1
          missing = .false.
          allmissing(3) = .false.
          kk = impt_lst(k,5)
          iaa(rs)%typ_impt(i+jj) = di_typ(kk)
          if (bonded_report.EQV..true.) then
            write(ilog,46) kk,iaa(rs)%impt(i+jj,1),iaa(rs)%impt(i+jj,2),&
 &iaa(rs)%impt(i+jj,3),iaa(rs)%impt(i+jj,4),rs
            write(ilog,48) b_type(iaa(rs)%impt(i+jj,1)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,1))),b_type(iaa(rs)%impt(i+jj,2)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,2))),b_type(iaa(rs)%impt(i+jj,3)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,3))),b_type(iaa(rs)%impt(i+jj,4)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,4))),amino(seqtyp(rs))
          end if
          do j=1,MAXDIPAR
            iaa(rs)%par_impt(i+jj,j) = di_par(kk,j)
          end do
!         convert assumed peptide convention RB to internal polymer convention RB
          if (di_typ(kk).eq.3) then
            iaa(rs)%par_impt(i+jj,2) = -iaa(rs)%par_impt(i+jj,2)
            iaa(rs)%par_impt(i+jj,4) = -iaa(rs)%par_impt(i+jj,4)
            iaa(rs)%par_impt(i+jj,6) = -iaa(rs)%par_impt(i+jj,6)
            iaa(rs)%par_impt(i+jj,8) = -iaa(rs)%par_impt(i+jj,8)
!         convert units from kcal/(mol*deg^2) to kcal/(mol*rad^2)
          else if (di_typ(kk).eq.2) then
            iaa(rs)%par_impt(i+jj,1)=iaa(rs)%par_impt(i+jj,1)/(RADIAN*RADIAN)
          end if
!         re-flag polymer convention Ryckaert-Bellemans
          if (di_typ(kk).eq.3) then
            iaa(rs)%typ_impt(i+jj) = 1
          end if
        end do
!     single degeneracy
      else if (((i2.ne.i3).OR.(i2.ne.i4).OR.(i3.ne.i4)).AND.&
 &        ((i2.eq.i3).OR.(i2.eq.i4).OR.(i3.eq.i4))) then
!       in this case, the first entry will be the unique one with no relevant degeneracy (see unbond.f90)
!       note that here the biotype permutation XYZY vs. XZYY is in fact relevant since the Y's are degenerate
!       hence the redudancy check is simply omitted for a quick+dirty solution
        do jj=0,2
          i2 = bio_botyp(b_type(iaa(rs)%impt(i+jj,2)))
          i3 = bio_botyp(b_type(iaa(rs)%impt(i+jj,3)))
          i4 = bio_botyp(b_type(iaa(rs)%impt(i+jj,4)))
          do k=1,impt_lstsz
            k1 = impt_lst(k,1)
            k2 = impt_lst(k,2)
            k3 = impt_lst(k,3)
            k4 = impt_lst(k,4)
            if (.NOT.(((i1.eq.k1).AND.(i2.eq.k2)).AND.&
 &                 ((i3.eq.k3).AND.(i4.eq.k4)))) cycle
            missing = .false.
            allmissing(3) = .false.
            kk = impt_lst(k,5)
            iaa(rs)%typ_impt(i+jj) = di_typ(kk)
            if (bonded_report.EQV..true.) then
              write(ilog,46) kk,iaa(rs)%impt(i+jj,1),iaa(rs)%impt(i+jj,2),&
 &iaa(rs)%impt(i+jj,3),iaa(rs)%impt(i+jj,4),rs
              write(ilog,48) b_type(iaa(rs)%impt(i+jj,1)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,1))),b_type(iaa(rs)%impt(i+jj,2)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,2))),b_type(iaa(rs)%impt(i+jj,3)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,3))),b_type(iaa(rs)%impt(i+jj,4)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,4))),amino(seqtyp(rs))
            end if
            do j=1,MAXDIPAR
              iaa(rs)%par_impt(i+jj,j) = di_par(kk,j)
            end do
!           convert assumed peptide convention RB to internal polymer convention RB
            if (di_typ(kk).eq.3) then
              iaa(rs)%par_impt(i+jj,2) = -iaa(rs)%par_impt(i+jj,2)
              iaa(rs)%par_impt(i+jj,4) = -iaa(rs)%par_impt(i+jj,4)
              iaa(rs)%par_impt(i+jj,6) = -iaa(rs)%par_impt(i+jj,6)
              iaa(rs)%par_impt(i+jj,8) = -iaa(rs)%par_impt(i+jj,8)
!           convert units from kcal/(mol*deg^2) to kcal/(mol*rad^2)
            else if (di_typ(kk).eq.2) then
              iaa(rs)%par_impt(i+jj,1)=iaa(rs)%par_impt(i+jj,1)/(RADIAN*RADIAN)
            end if
!           re-flag polymer convention Ryckaert-Bellemans
            if (di_typ(kk).eq.3) then
              iaa(rs)%typ_impt(i+jj) = 1
            end if
          end do
        end do
!     no (relevant) degeneracy
      else
        do jj=0,2
          i2 = bio_botyp(b_type(iaa(rs)%impt(i+jj,2)))
          i3 = bio_botyp(b_type(iaa(rs)%impt(i+jj,3)))
          i4 = bio_botyp(b_type(iaa(rs)%impt(i+jj,4)))
          do k=1,impt_lstsz
            k1 = impt_lst(k,1)
            k2 = impt_lst(k,2)
            k3 = impt_lst(k,3)
            k4 = impt_lst(k,4)
            if (.NOT.((((i1.eq.k1).AND.(i2.eq.k2)).AND.&
 &                 ((i3.eq.k3).AND.(i4.eq.k4))).OR.&
 &                (((i1.eq.k1).AND.(i2.eq.k3)).AND.&
 &                 ((i3.eq.k2).AND.(i4.eq.k4)))) ) cycle
            reverse = .true.
            if ((i2.eq.k2).AND.(i3.eq.k3)) reverse = .false.
            if ((improper_conv(1).ne.1).AND.(reverse.EQV..true.)) then
              j = iaa(rs)%impt(i+jj,2)
              iaa(rs)%impt(i+jj,2) = iaa(rs)%impt(i+jj,3)
              iaa(rs)%impt(i+jj,3) = j
            end if
            missing = .false.
            allmissing(3) = .false.
            kk = impt_lst(k,5)
            iaa(rs)%typ_impt(i+jj) = di_typ(kk)
            if (bonded_report.EQV..true.) then
              write(ilog,46) kk,iaa(rs)%impt(i+jj,1),iaa(rs)%impt(i+jj,2),&
 &iaa(rs)%impt(i+jj,3),iaa(rs)%impt(i+jj,4),rs
              write(ilog,48) b_type(iaa(rs)%impt(i+jj,1)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,1))),b_type(iaa(rs)%impt(i+jj,2)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,2))),b_type(iaa(rs)%impt(i+jj,3)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,3))),b_type(iaa(rs)%impt(i+jj,4)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,4))),amino(seqtyp(rs))
            end if
            do j=1,MAXDIPAR
              iaa(rs)%par_impt(i+jj,j) = di_par(kk,j)
            end do
!           convert assumed peptide convention RB to internal polymer convention RB
            if (di_typ(kk).eq.3) then
              iaa(rs)%par_impt(i+jj,2) = -iaa(rs)%par_impt(i+jj,2)
              iaa(rs)%par_impt(i+jj,4) = -iaa(rs)%par_impt(i+jj,4)
              iaa(rs)%par_impt(i+jj,6) = -iaa(rs)%par_impt(i+jj,6)
              iaa(rs)%par_impt(i+jj,8) = -iaa(rs)%par_impt(i+jj,8)
!           convert units from kcal/(mol*rad^2) to kcal/(mol*deg^2) and flip sign of 
!           equilbirium position based on crude heuristic (otherwise, we would have to
!           maintain six independent impropers for each site)
            else if (di_typ(kk).eq.2) then
              iaa(rs)%par_impt(i+jj,1)=iaa(rs)%par_impt(i+jj,1)/(RADIAN*RADIAN)
              if (iaa(rs)%par_impt(i+jj,2).ne.0.0) then
                testimp = getztor(iaa(rs)%impt(i+jj,1),iaa(rs)%impt(i+jj,2),&
 &iaa(rs)%impt(i+jj,3),iaa(rs)%impt(i+jj,4))
                if (testimp*iaa(rs)%par_impt(i+jj,2).lt.0.0) then
                  iaa(rs)%par_impt(i+jj,2) = -iaa(rs)%par_impt(i+jj,2)
                end if
              end if
            end if
!           re-flag polymer convention Ryckaert-Bellemans
            if (di_typ(kk).eq.3) then
              iaa(rs)%typ_impt(i+jj) = 1
            end if
          end do
        end do
      end if
      if (missing.EQV..true.) then
        if (bonded_report.EQV..true.) then
          write(ilog,47) iaa(rs)%impt(i,1),iaa(rs)%impt(i,2),&
 &iaa(rs)%impt(i,3),iaa(rs)%impt(i,4),rs
          write(ilog,48) b_type(iaa(rs)%impt(i,1)),&
 &bio_code(b_type(iaa(rs)%impt(i,1))),b_type(iaa(rs)%impt(i,2)),&
 &bio_code(b_type(iaa(rs)%impt(i,2))),b_type(iaa(rs)%impt(i,3)),&
 &bio_code(b_type(iaa(rs)%impt(i,3))),b_type(iaa(rs)%impt(i,4)),&
 &bio_code(b_type(iaa(rs)%impt(i,4))),amino(seqtyp(rs))
        end if
      end if
    end do
  end do
  if ((allmissing(3).EQV..true.).AND.(use_BOND(3).EQV..true.)) then
    write(ilog,*) 'No applicable assignments for improper dihedral potentials found in parameter file ',paramfile(t1:t2),'.'
  end if
!
  
! torsional terms (optional)
  if ((bonded_report.EQV..true.).AND.(use_BOND(4).EQV..true.)) then
    write(ilog,*)
    write(ilog,*) '--------     Torsional Terms     --------'
  end if

  do rs=1,nseq
    if (use_BOND(4).EQV..false.) exit
    do i=1,nrsdi(rs)
      if (check_colinear(iaa(rs)%di(i,1),iaa(rs)%di(i,2),iaa(rs)%di(i,3),iaa(rs)%di(i,4)).EQV..true.) then
        write(ilog,27) iaa(rs)%di(i,1:4),rs
        cycle
      end if
      missing = .true.
      
      i1 = bio_botyp(b_type(iaa(rs)%di(i,1)))
      i2 = bio_botyp(b_type(iaa(rs)%di(i,2)))
      i3 = bio_botyp(b_type(iaa(rs)%di(i,3)))
      i4 = bio_botyp(b_type(iaa(rs)%di(i,4)))
      
      do k=1,di_lstsz
        k1 = di_lst(k,1)
        k2 = di_lst(k,2)
        k3 = di_lst(k,3)
        k4 = di_lst(k,4)
        if (.NOT.((((i2.eq.k2).AND.(i3.eq.k3)).AND.&
 &                 ((i1.eq.k1).AND.(i4.eq.k4))).OR.&
 &                (((i2.eq.k3).AND.(i3.eq.k2)).AND.&
 &                 ((i1.eq.k4).AND.(i4.eq.k1))))) cycle
        reverse = .true.
        if ((i2.eq.k2).AND.(i3.eq.k3)) reverse = .false.
        missing = .false.
        allmissing(4) = .false.
        kk = di_lst(k,5)
        iaa(rs)%typ_di(i) = di_typ(kk)
        if (bonded_report.EQV..true.) then
          write(ilog,56) kk,iaa(rs)%di(i,1),iaa(rs)%di(i,2),&
 &iaa(rs)%di(i,3),iaa(rs)%di(i,4),rs
          write(ilog,58) b_type(iaa(rs)%di(i,1)),&
 &bio_code(b_type(iaa(rs)%di(i,1))),b_type(iaa(rs)%di(i,2)),&
 &bio_code(b_type(iaa(rs)%di(i,2))),b_type(iaa(rs)%di(i,3)),&
 &bio_code(b_type(iaa(rs)%di(i,3))),b_type(iaa(rs)%di(i,4)),&
 &bio_code(b_type(iaa(rs)%di(i,4))),amino(seqtyp(rs))
             write(ilog,*) bio_botyp(b_type(iaa(rs)%di(i,1))),bio_botyp(b_type(iaa(rs)%di(i,2))),bio_botyp(b_type(iaa(rs)%di(i,3))),&
            &bio_botyp(b_type(iaa(rs)%di(i,4)))
!              write(ilog,777) bio_code(b_type(iaa(rs)%di(i,1))),&
! &bio_code(b_type(iaa(rs)%di(i,2))),bio_code(b_type(iaa(rs)%di(i,3))),&
! &bio_code(b_type(iaa(rs)%di(i,4))),di_par(kk,1:4)

        end if
!        write(*,*) di_par(kk,1:6)
        do j=1,MAXDIPAR
          iaa(rs)%par_di(i,j) = di_par(kk,j)
        end do
!       convert Ryckaert-Bellemans from assumed polymer convention to internal peptide convention
        if (di_typ(kk).eq.3) then
          iaa(rs)%par_di(i,2) = -iaa(rs)%par_di(i,2)
          iaa(rs)%par_di(i,4) = -iaa(rs)%par_di(i,4)
          iaa(rs)%par_di(i,6) = -iaa(rs)%par_di(i,6)
          iaa(rs)%par_di(i,8) = -iaa(rs)%par_di(i,8)
!       convert units from kcal/(mol*deg^2) to kcal/(mol*rad^2)
        else if (di_typ(kk).eq.2) then
          iaa(rs)%par_di(i,1)=iaa(rs)%par_di(i,1)/(RADIAN*RADIAN)
        end if
!       re-flag polymer convention Ryckaert-Bellemans
        if (di_typ(kk).eq.3) then
          iaa(rs)%typ_di(i) = 1
        end if
      end do
      if (missing.EQV..true.) then
        if (bonded_report.EQV..true.) then
          write(ilog,57) iaa(rs)%di(i,1),iaa(rs)%di(i,2),&
 &iaa(rs)%di(i,3),iaa(rs)%di(i,4),rs
          write(ilog,58) b_type(iaa(rs)%di(i,1)),&
 &bio_code(b_type(iaa(rs)%di(i,1))),b_type(iaa(rs)%di(i,2)),&
 &bio_code(b_type(iaa(rs)%di(i,2))),b_type(iaa(rs)%di(i,3)),&
 &bio_code(b_type(iaa(rs)%di(i,3))),b_type(iaa(rs)%di(i,4)),&
 &bio_code(b_type(iaa(rs)%di(i,4))),amino(seqtyp(rs))
            write(ilog,*) bio_botyp(b_type(iaa(rs)%di(i,1))),bio_botyp(b_type(iaa(rs)%di(i,2))),bio_botyp(b_type(iaa(rs)%di(i,3))),&
            &bio_botyp(b_type(iaa(rs)%di(i,4)))
        end if
      end if
    end do
  end do
  if ((allmissing(4).EQV..true.).AND.(use_BOND(4).EQV..true.)) then
    write(ilog,*) 'No applicable assignments for torsional angle potentials found in parameter file ',paramfile(t1:t2),'.'
  end if
!
! CMAP terms (utterly optional)
  if ((bonded_report.EQV..true.).AND.(use_BOND(5).EQV..true.)) then
    write(ilog,*)
    write(ilog,*) '--------       CMAP terms         --------'
  end if
  do rs=1,nseq
    if (use_BOND(5).EQV..false.) exit
    do i=1,nrscm(rs)
      missing = .true.
      i1 = bio_botyp(b_type(iaa(rs)%cm(i,1)))
      i2 = bio_botyp(b_type(iaa(rs)%cm(i,2)))
      i3 = bio_botyp(b_type(iaa(rs)%cm(i,3)))
      i4 = bio_botyp(b_type(iaa(rs)%cm(i,4)))
      i5 = bio_botyp(b_type(iaa(rs)%cm(i,5)))
      do k=1,cm_lstsz
        k1 = cm_lst(k,1)
        k2 = cm_lst(k,2)
        k3 = cm_lst(k,3)
        k4 = cm_lst(k,4)
        k5 = cm_lst(k,5)
!       match needs to be exact
        if (.NOT.(((i2.eq.k2).AND.(i3.eq.k3)).AND.&
 &                ((i1.eq.k1).AND.(i4.eq.k4)).AND.(i5.eq.k5)) ) cycle
        missing = .false.
        allmissing(5) = .false.
        kk = cm_lst(k,6)
        iaa(rs)%typ_cm(i) = kk
        if (bonded_report.EQV..true.) then
          write(ilog,36) iaa(rs)%typ_cm(i),iaa(rs)%cm(i,1),iaa(rs)%cm(i,2),&
 &iaa(rs)%cm(i,3),iaa(rs)%cm(i,4),iaa(rs)%cm(i,5),atmres(iaa(rs)%cm(i,3))
          write(ilog,38) b_type(iaa(rs)%cm(i,1)),&
 &bio_code(b_type(iaa(rs)%cm(i,1))),b_type(iaa(rs)%cm(i,2)),&
 &bio_code(b_type(iaa(rs)%cm(i,2))),b_type(iaa(rs)%cm(i,3)),&
 &bio_code(b_type(iaa(rs)%cm(i,3))),b_type(iaa(rs)%cm(i,4)),&
 &bio_code(b_type(iaa(rs)%cm(i,4))),b_type(iaa(rs)%cm(i,5)),&
 &bio_code(b_type(iaa(rs)%cm(i,5))),amino(seqtyp(atmres(iaa(rs)%cm(i,3))))
        end if
      end do
!     never report missing CMAPs 
    end do
  end do
  if ((allmissing(5).EQV..true.).AND.(use_BOND(5).EQV..true.)) then
    write(ilog,*) 'No applicable assignments for CMAP correction potentials found in parameter file ',paramfile(t1:t2),'.'
  end if
!
  

  
  
  
! eventually patch the default assignments from parameter file
  call read_bondedpatchfile(allmissing)
!
  do rs=1,nseq
!   sort terms such that ones with energy terms assigned are first
!   use second counter variable to mimick smaller array size
    if (use_BOND(1).EQV..true.) then
      k1 = 0
      k2 = 0
      allocate(pams(nrsbl(rs),MAXBOPAR+1))
      allocate(dum(nrsbl(rs),3))
      allocate(dum2(nrsbl(rs),2))
      do i=1,nrsbl(rs)
        if (iaa(rs)%typ_bl(i).gt.0) then ! MArtin  there might be something with this 
        ! because if the bonds are reshuffled(nad they are) the [presence of the dummy early on could explain yhe trouble 
          k1 = k1 + 1
          do j=1,2
            dum(k1,j) = iaa(rs)%bl(i,j)
          end do
          dum(k1,3) = iaa(rs)%typ_bl(i)
          do j=1,MAXBOPAR+1
              
            pams(k1,j) = iaa(rs)%par_bl(i,j)
          end do
        else
          k2 = k2 + 1
          do j=1,2
            dum2(k2,j) = iaa(rs)%bl(i,j)
          end do
        end if
      end do
      do i=1,k1
        do j=1,2
          iaa(rs)%bl(i,j) = dum(i,j)
        end do
        iaa(rs)%typ_bl(i) = dum(i,3)
        do j=1,MAXBOPAR+1
          iaa(rs)%par_bl(i,j) = pams(i,j)
        end do
      end do
      do i=k1+1,k1+k2
        do j=1,2
          iaa(rs)%bl(i,j) = dum2(i-k1,j)
        end do
        iaa(rs)%typ_bl(i) = 0
      end do
      nrsbleff(rs) = k1
      deallocate(pams)
      deallocate(dum)
      deallocate(dum2)
    end if
    if (use_BOND(2).EQV..true.) then
      k1 = 0
      k2 = 0
      allocate(pams(nrsba(rs),MAXBAPAR+1))
      allocate(dum(nrsba(rs),4))
      allocate(dum2(nrsba(rs),3))
      do i=1,nrsba(rs)
        if (iaa(rs)%typ_ba(i).gt.0) then
          k1 = k1 + 1
          do j=1,3
            dum(k1,j) = iaa(rs)%ba(i,j)
          end do
          dum(k1,4) = iaa(rs)%typ_ba(i)
          do j=1,MAXBAPAR+1
            pams(k1,j) = iaa(rs)%par_ba(i,j)
          end do
        else
          k2 = k2 + 1
          do j=1,3
            dum2(k2,j) = iaa(rs)%ba(i,j)
          end do
        end if
      end do
      do i=1,k1
        do j=1,3
          iaa(rs)%ba(i,j) = dum(i,j)
        end do
        iaa(rs)%typ_ba(i) = dum(i,4)
        do j=1,MAXBAPAR+1
          iaa(rs)%par_ba(i,j) = pams(i,j)
        end do
      end do
      do i=k1+1,k1+k2
        do j=1,3
          iaa(rs)%ba(i,j) = dum2(i-k1,j)
        end do
        iaa(rs)%typ_ba(i) = 0
      end do
      nrsbaeff(rs) = k1
      deallocate(pams)
      deallocate(dum)
      deallocate(dum2)
    end if
    if (use_BOND(3).EQV..true.) then
      k1 = 0
      k2 = 0
      allocate(pams(nrsimpt(rs),MAXDIPAR))
      allocate(dum(nrsimpt(rs),5))
      allocate(dum2(nrsimpt(rs),4))
      do i=1,nrsimpt(rs)
        if (iaa(rs)%typ_impt(i).gt.0) then
          k1 = k1 + 1
          do j=1,4
            dum(k1,j) = iaa(rs)%impt(i,j)
          end do
          dum(k1,5) = iaa(rs)%typ_impt(i)
          do j=1,MAXDIPAR
            pams(k1,j) = iaa(rs)%par_impt(i,j)
          end do
        else
          k2 = k2 + 1
          do j=1,4
            dum2(k2,j) = iaa(rs)%impt(i,j)
          end do
        end if
      end do
      do i=1,k1
        do j=1,4
          iaa(rs)%impt(i,j) = dum(i,j)
        end do
        iaa(rs)%typ_impt(i) = dum(i,5)
        do j=1,MAXDIPAR
          iaa(rs)%par_impt(i,j) = pams(i,j)
        end do
      end do
      do i=k1+1,k1+k2
        do j=1,4
          iaa(rs)%impt(i,j) = dum2(i-k1,j)
        end do
        iaa(rs)%typ_impt(i) = 0
      end do
      nrsimpteff(rs) = k1
      deallocate(pams)
      deallocate(dum)
      deallocate(dum2)
    end if
    if (use_BOND(4).EQV..true.) then
      k1 = 0
      k2 = 0
      allocate(pams(nrsdi(rs),MAXDIPAR))
      allocate(dum(nrsdi(rs),5))
      allocate(dum2(nrsdi(rs),4))
      do i=1,nrsdi(rs)
        if (iaa(rs)%typ_di(i).gt.0) then
          k1 = k1 + 1
          do j=1,4
            dum(k1,j) = iaa(rs)%di(i,j)
          end do
          dum(k1,5) = iaa(rs)%typ_di(i)
          do j=1,MAXDIPAR
            pams(k1,j) = iaa(rs)%par_di(i,j)
          end do
        else
          k2 = k2 + 1
          do j=1,4
            dum2(k2,j) = iaa(rs)%di(i,j)
          end do
        end if
      end do
      do i=1,k1
        do j=1,4
          iaa(rs)%di(i,j) = dum(i,j)
        end do
        iaa(rs)%typ_di(i) = dum(i,5)
        do j=1,MAXDIPAR
          iaa(rs)%par_di(i,j) = pams(i,j)
        end do
      end do
      do i=k1+1,k1+k2
        do j=1,4
          iaa(rs)%di(i,j) = dum2(i-k1,j)
        end do
        iaa(rs)%typ_di(i) = 0
      end do
      nrsdieff(rs) = k1
      deallocate(pams)
      deallocate(dum)
      deallocate(dum2)
    end if
    if (use_BOND(5).EQV..true.) then
      k1 = 0
      k2 = 0
      allocate(dum(nrscm(rs),6))
      allocate(dum2(nrscm(rs),5))
      do i=1,nrscm(rs)
        if (iaa(rs)%typ_cm(i).gt.0) then
          k1 = k1 + 1
          do j=1,5
            dum(k1,j) = iaa(rs)%cm(i,j)
          end do
          dum(k1,6) = iaa(rs)%typ_cm(i)
        else
          k2 = k2 + 1
          do j=1,5
            dum2(k2,j) = iaa(rs)%cm(i,j)
          end do
        end if
      end do
      do i=1,k1
        do j=1,5
          iaa(rs)%cm(i,j) = dum(i,j)
        end do
        iaa(rs)%typ_cm(i) = dum(i,6)
      end do
      do i=k1+1,k1+k2
        do j=1,5
          iaa(rs)%cm(i,j) = dum2(i-k1,j)
        end do
        iaa(rs)%typ_cm(i) = 0
      end do
      nrscmeff(rs) = k1
      deallocate(dum)
      deallocate(dum2)
    end if
  end do
!
  if ((allmissing(1).EQV..true.).AND.(use_BOND(1).EQV..true.).AND.&
 &    (fycxyz.ne.2)) then
    use_BOND(1) = .false.
    scale_BOND(1) = 0.0
  end if
  if ((allmissing(2).EQV..true.).AND.(use_BOND(2).EQV..true.).AND.&
 &    (fycxyz.ne.2)) then
    use_BOND(2) = .false.
    scale_BOND(2) = 0.0
  end if
  if ((allmissing(3).EQV..true.).AND.(use_BOND(3).EQV..true.)) then
    use_BOND(3) = .false.
    scale_BOND(3) = 0.0
  end if
  if ((allmissing(4).EQV..true.).AND.(use_BOND(4).EQV..true.)) then
    use_BOND(4) = .false.
    scale_BOND(4) = 0.0
  end if
  if ((allmissing(5).EQV..true.).AND.(use_BOND(5).EQV..true.)) then
    use_BOND(5) = .false.
    scale_BOND(5) = 0.0
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine guess_types(inval,inmass,innam,ljty,boty)
!
  use params
!
  implicit none
!
  character(4) innam
  integer inval,ljty,boty,k
  RTYPE inmass,thresh
  
  
  
  ljty = -1
  thresh = 0.01
  do while (thresh.lt.HUGE(thresh)/100.0) 
    do k=1,n_ljtyp
      if ((abs(inmass-lj_weight(k)).le.thresh).AND.(lj_val(k).eq.inval)) then
        ljty = k
        exit
      end if
    end do
    if (ljty.eq.-1) then ! has to trigger eventually
      do k=1,n_ljtyp
        if (abs(inmass-lj_weight(k)).le.thresh) then
          ljty = k
          exit
        end if
      end do
    end if
    if (ljty.eq.-1) then
      do k=1,n_ljtyp
        if (inval.eq.lj_val(k)) then
          ljty = k
          exit
        end if
      end do
    end if
    if (ljty.gt.0) exit
    thresh = thresh*10.0
  end do
!
  n_biotyp = n_biotyp + 1
  boty = n_biotyp
  bio_code(boty) = innam(2:4)
  bio_ljtyp(boty) = ljty
  bio_ctyp(boty) = 0
  bio_botyp(boty) = 0
!
end
!
!-----------------------------------------------------------------------------------------------------
