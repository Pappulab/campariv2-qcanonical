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
!-----------------------------------------------------------------------
!
! THE MAIN SETUP ROUTINE FOR ASSIGNMENT OF PARTIAL CHARGES
!
!-----------------------------------------------------------------------
!
subroutine polar_groups()
!
  use iounit
  use polypep
  use sequen
  use inter
  use aminos
  use atoms
  use energies
  use params
  use molecule
  use system
  use fyoc
  use zmatrix
  
  
  use energies ! Martin added just for TOTQX
!
  
  use mpistuff ! Martin : added for debug only   
  
  implicit none
!
  integer i,j,k,l,ii,kk,rs,imol,cty,shf,shf3,shf5
  integer, ALLOCATABLE:: skip(:),tmpvec(:)
  type(t_at), ALLOCATABLE :: at_temp(:)
  integer nccat,ncan,npcat,npan,dp2,dp1,alcsz,offset
  logical proceed,foundit,its14,docycle
  character(3) resname,resname1
  RTYPE tc,eps
  
  
  ! Martin : added 
  integer number
  RTYPE chg
  logical test1,test2

  shf = 0
  shf3 = 0
  shf5 = 0
  if (ua_model.gt.0) then
    shf = 1
    shf3 = 3
    shf5 = 5
  end if
  eps = 1.0*(10.0**(-precision(eps)))
!
  foundit = .false.
!
  do i=1,n
    atq(i) = 0.0
    cty = bio_ctyp(b_type(i))
    if (cty.eq.0) then
      if (seqtyp(atmres(i)).eq.26) then
        if (foundit.EQV..false.) write(ilog,*) 'Warning. Parameter file cannot be used to set partial charge for atoms&
 & in unsupported residues (first for atom #',i,'). This run will terminate unless a suitable charge patch is found.&
 & Omitting further warnings of this type.'
        foundit = .true.
      else

        write(ilog,*) 'Warning. Parameter file has no partial charge for biotype&
 & ',b_type(i),' (',bio_code(b_type(i)),') (encountered for atom&
 & #',i,'). This run will terminate unless a suitable charge patch is found.'
      end if
    end if
  end do
!
! since these parameters are assigned automatically, it is crucial that
! we leave ourselves an option to check. after all, while it's easy to come
! up with arbitrary charges, it is certainly not easy to make sure that the
! interaction arrays and charge sets are all sane and fully supported at all
! levels

! loop over all molecules in the system
  do imol=1,nmol
!
!   and over all residues within each molecule
    do rs=rsmol(imol,1),rsmol(imol,2)
!
      resname = amino(seqtyp(rs))

!     first deal with unsupported residues (charges may be available through biotype patch)
      if (resname.eq.'UNK') then
        at(rs)%npol = 0
        do i=1,at(rs)%nbb+at(rs)%nsc
          if (i.le.at(rs)%nbb) then
            ii = at(rs)%bb(i)
          else
            ii = at(rs)%sc(i-at(rs)%nbb)
          end if
          
          cty = bio_ctyp(b_type(ii))
          !! Ok so I need to reformulate that so the charges are taken into account in the case of a partial charge, but also 
          !! in the case of a monopole. One way could be to assign fakes charges to the atoms in the 
          !! changing partial charges, to "simulate" a charge, that keeps the overall charge of the charged state
          !! I also want to get rid of the special amino acid name thingy
          !! better solution : at the end of the charge assignement, go through residues, find those that are 
          !! 0 in the second limits, find how many pol in the dipolar group, and assign same charge, but equivalently to all at
          if (cty.eq.0) cycle
          if ((abs(c_charge(cty)).gt.eps).or.(cty.eq.275))  then
            at(rs)%npol = at(rs)%npol + 1
            at(rs)%pol(at(rs)%npol) = ii
            atq(ii) = c_charge(cty)
          end if
        end do
        cycle
      end if
!
      if (rsmol(imol,1).eq.rsmol(imol,2)) then
!
!       deal with small molecules first
        if ((seqflag(rs).ge.101).AND.(seqflag(rs).le.103)) then
          at(rs)%npol = 0
          do i=1,at(rs)%nbb+at(rs)%nsc
            if (i.le.at(rs)%nbb) then
              ii = at(rs)%bb(i)
            else
              ii = at(rs)%sc(i-at(rs)%nbb)
            end if
            cty = bio_ctyp(b_type(ii))
            if (cty.eq.0) cycle
            if ((abs(c_charge(cty)).gt.eps).or.(cty.eq.275))  then
                at(rs)%npol = at(rs)%npol + 1
                at(rs)%pol(at(rs)%npol) = ii
                atq(ii) = c_charge(cty)
            end if
          end do
!
!       now onto free amino acids
        else if ((seqflag(rs).eq.2).OR.(seqflag(rs).eq.8)) then
          at(rs)%npol = 0
          do i=1,at(rs)%nbb+at(rs)%nsc
            if (i.le.at(rs)%nbb) then
              ii = at(rs)%bb(i)
            else
              ii = at(rs)%sc(i-at(rs)%nbb)
            end if
            cty = bio_ctyp(b_type(ii))
            if (cty.eq.0) cycle
            if ((abs(c_charge(cty)).gt.eps).or.(cty.eq.275))  then
                at(rs)%npol = at(rs)%npol + 1
                at(rs)%pol(at(rs)%npol) = ii
                atq(ii) = c_charge(cty)
            end if
          end do
          if ((moltermid(imol,1).ne.1).AND.(moltermid(imol,1).ne.2)) then
            write(ilog,*) 'Fatal. Encountered unsupported terminus t&
 &ype for molecule ',imol,' in polar_groups() (N:',moltermid(imol,1),').'
            call fexit()
          end if
          if (moltermid(imol,2).ne.1) then
            write(ilog,*) 'Fatal. Encountered unsupported terminus t&
 &ype for molecule ',imol,' in polar_groups() (C:',moltermid(imol,2),').'
            call fexit()
          end if
!
        else if (seqflag(rs).eq.5) then
          at(rs)%npol = 0
          do i=1,at(rs)%nbb+at(rs)%nsc
            if (i.le.at(rs)%nbb) then
              ii = at(rs)%bb(i)
            else
              ii = at(rs)%sc(i-at(rs)%nbb)
            end if
            cty = bio_ctyp(b_type(ii))
            if (cty.eq.0) cycle
            if ((abs(c_charge(cty)).gt.eps).or.(cty.eq.275))  then
                at(rs)%npol = at(rs)%npol + 1
                at(rs)%pol(at(rs)%npol) = ii
                atq(ii) = c_charge(cty)
            end if
          end do
          if (moltermid(imol,1).eq.2) then
            write(ilog,*) 'Fatal. Polar potential does not support a&
 &n uncharged terminus with residue ',resname,' at the moment.'
            call fexit()
          else if (moltermid(imol,1).ne.1) then
            write(ilog,*) 'Fatal. Encountered unsupported terminus t&
 &ype for molecule ',imol,' in polar_groups() (N:',moltermid(imol,1),').'
            call fexit()
          end if
          if (moltermid(imol,2).ne.1) then
            write(ilog,*) 'Fatal. Encountered unsupported terminus t&
 &ype for molecule ',imol,' in polar_groups() (C:',moltermid(imol,2),').'
            call fexit()
          end if
!
        else if (seqflag(rs).eq.22) then
          at(rs)%npol = 0
          do i=1,at(rs)%nbb+at(rs)%nsc
            if (i.le.at(rs)%nbb) then
              ii = at(rs)%bb(i)
            else
              ii = at(rs)%sc(i-at(rs)%nbb)
            end if
            cty = bio_ctyp(b_type(ii))
            if (cty.eq.0) cycle
            if ((abs(c_charge(cty)).gt.eps).or.(cty.eq.275))  then
                at(rs)%npol = at(rs)%npol + 1
                at(rs)%pol(at(rs)%npol) = ii
                atq(ii) = c_charge(cty)
            end if
          end do
          if (moltermid(imol,1).ne.1) then
            write(ilog,*) 'Fatal. Encountered unsupported terminus t&
 &ype for molecule ',imol,' in polar_groups() (5":',moltermid(imol,1),').'
            call fexit()
          end if
          if (moltermid(imol,2).ne.1) then
            write(ilog,*) 'Fatal. Encountered unsupported terminus t&
 &ype for molecule ',imol,' in polar_groups() (3":',moltermid(imol,2),').'
            call fexit()
          end if
        else
          write(ilog,*) 'Fatal. Polar potential does not support res&
 &idue ',resname,' (free) at the moment. Check back later.'
          call fexit()
        end if
!
!
!     N-terminal residues
!
      else if (rs.eq.rsmol(imol,1)) then
!
!       these controls are for terminus conditions and to warn developers when adding unusual stuff
!
        if ((seqflag(rs).eq.10).OR.(seqflag(rs).eq.24)) then
          at(rs)%npol = 0
          do i=1,at(rs)%nbb+at(rs)%nsc
            if (i.le.at(rs)%nbb) then
              ii = at(rs)%bb(i)
            else
              ii = at(rs)%sc(i-at(rs)%nbb)
            end if
            cty = bio_ctyp(b_type(ii))
            if (cty.eq.0) cycle
            if ((abs(c_charge(cty)).gt.eps).or.(cty.eq.275))  then
                at(rs)%npol = at(rs)%npol + 1
                at(rs)%pol(at(rs)%npol) = ii
                atq(ii) = c_charge(cty)
            end if
          end do
!
        else if ((seqflag(rs).eq.2).OR.(seqflag(rs).eq.8)) then
          at(rs)%npol = 0
          do i=1,at(rs)%nbb+at(rs)%nsc
            if (i.le.at(rs)%nbb) then
              ii = at(rs)%bb(i)
            else
              ii = at(rs)%sc(i-at(rs)%nbb)
            end if
            cty = bio_ctyp(b_type(ii))
            if (cty.eq.0) cycle
            if ((abs(c_charge(cty)).gt.eps).or.(cty.eq.275))  then
                at(rs)%npol = at(rs)%npol + 1
                at(rs)%pol(at(rs)%npol) = ii
                atq(ii) = c_charge(cty)
            end if
          end do
          if ((moltermid(imol,1).ne.1).AND.(moltermid(imol,1).ne.2)) then
            write(ilog,*) 'Fatal. Encountered unsupported terminus t&
 &ype for molecule ',imol,' in polar_groups() (N:',moltermid(imol,1),').'
            call fexit()
          end if
!
        else if (seqflag(rs).eq.5) then
          at(rs)%npol = 0
          do i=1,at(rs)%nbb+at(rs)%nsc
            if (i.le.at(rs)%nbb) then
              ii = at(rs)%bb(i)
            else
              ii = at(rs)%sc(i-at(rs)%nbb)
            end if
            cty = bio_ctyp(b_type(ii))
            if (cty.eq.0) cycle
            if ((abs(c_charge(cty)).gt.eps).or.(cty.eq.275))  then
                at(rs)%npol = at(rs)%npol + 1
                at(rs)%pol(at(rs)%npol) = ii
                atq(ii) = c_charge(cty)
            end if
          end do
          if (moltermid(imol,1).eq.2) then
            write(ilog,*) 'Fatal. Polar potential does not support a&
 &n uncharged N-terminus with residue ',resname,' at the moment.'
            call fexit()
          else if (moltermid(imol,1).ne.1) then
            write(ilog,*) 'Fatal. Encountered unsupported terminus t&
 &ype for molecule ',imol,' in polar_groups() (N:',moltermid(imol,1),').'
            call fexit()
          end if
!
        else if (seqflag(rs).eq.22) then
          at(rs)%npol = 0
          do i=1,at(rs)%nbb+at(rs)%nsc
            if (i.le.at(rs)%nbb) then
              ii = at(rs)%bb(i)
            else
              ii = at(rs)%sc(i-at(rs)%nbb)
            end if
            cty = bio_ctyp(b_type(ii))
            if (cty.eq.0) cycle
            if ((abs(c_charge(cty)).gt.eps).or.(cty.eq.275))  then
                at(rs)%npol = at(rs)%npol + 1
                at(rs)%pol(at(rs)%npol) = ii
                atq(ii) = c_charge(cty)
            end if
          end do
          if (moltermid(imol,1).ne.1) then
            write(ilog,*) 'Fatal. Encountered unsupported terminus t&
 &ype for molecule ',imol,' in polar_groups() (5":',moltermid(imol,1),').'
            call fexit()
          end if
!
        else
          write(ilog,*) 'Fatal. Polar potential does not support N-t&
 &erminal residue ',resname,' (?) at the moment. Check back later.'
          call fexit()
        end if
!
!     C-terminal residues
      else if (rs.eq.rsmol(imol,2)) then
!
!       these controls are for terminus conditions and to warn developers when adding unusual stuff
!
        if ((seqflag(rs).eq.11).OR.(seqflag(rs).eq.13)) then
          at(rs)%npol = 0
          do i=1,at(rs)%nbb+at(rs)%nsc
            if (i.le.at(rs)%nbb) then
              ii = at(rs)%bb(i)
            else
              ii = at(rs)%sc(i-at(rs)%nbb)
            end if
            cty = bio_ctyp(b_type(ii))
            if (cty.eq.0) cycle
            if ((abs(c_charge(cty)).gt.eps).or.(cty.eq.275))  then
                at(rs)%npol = at(rs)%npol + 1
                at(rs)%pol(at(rs)%npol) = ii
                atq(ii) = c_charge(cty)
            end if
          end do
!
        else if ((seqflag(rs).eq.2).OR.(seqflag(rs).eq.5).OR.(seqflag(rs).eq.8)) then
          at(rs)%npol = 0
          do i=1,at(rs)%nbb+at(rs)%nsc
            if (i.le.at(rs)%nbb) then
              ii = at(rs)%bb(i)
            else
              ii = at(rs)%sc(i-at(rs)%nbb)
            end if
            cty = bio_ctyp(b_type(ii))
            if (cty.eq.0) cycle
            if ((abs(c_charge(cty)).gt.eps).or.(cty.eq.275))  then
                at(rs)%npol = at(rs)%npol + 1
                at(rs)%pol(at(rs)%npol) = ii
                atq(ii) = c_charge(cty)
            end if
          end do
          if (moltermid(imol,2).ne.1) then
            write(ilog,*) 'Fatal. Encountered unsupported terminus t&
 &ype for molecule ',imol,' in polar_groups() (C:',moltermid(imol,2),').'
            call fexit()
          end if
!
        else if (seqflag(rs).eq.22) then
          at(rs)%npol = 0
          do i=1,at(rs)%nbb+at(rs)%nsc
            if (i.le.at(rs)%nbb) then
              ii = at(rs)%bb(i)
            else
              ii = at(rs)%sc(i-at(rs)%nbb)
            end if
            cty = bio_ctyp(b_type(ii))
            if (cty.eq.0) cycle
            if ((abs(c_charge(cty)).gt.eps).or.(cty.eq.275))  then
                at(rs)%npol = at(rs)%npol + 1
                at(rs)%pol(at(rs)%npol) = ii
                atq(ii) = c_charge(cty)
            end if
          end do
          if (moltermid(imol,2).ne.1) then
            write(ilog,*) 'Fatal. Encountered unsupported terminus t&
 &ype for molecule ',imol,' in polar_groups() (3":',moltermid(imol,2),').'
            call fexit()
          end if
        else
          write(ilog,*) 'Fatal. Polar potential does not support C-t&
 &erminal residue ',resname,' (?) at the moment. Check back later.'
          call fexit()
        end if
!
!     and finally all middle-of-the-chain residues (these controls are purely to warn developers when adding unusual stuff)
!
      else if ((seqflag(rs).eq.2).OR.(seqflag(rs).eq.5).OR.(seqflag(rs).eq.8).OR.(seqflag(rs).eq.22)) then
        at(rs)%npol = 0
        do i=1,at(rs)%nbb+at(rs)%nsc
          if (i.le.at(rs)%nbb) then
            ii = at(rs)%bb(i)
          else
            ii = at(rs)%sc(i-at(rs)%nbb)
          end if
          cty = bio_ctyp(b_type(ii))
          if (cty.eq.0) cycle
            if ((abs(c_charge(cty)).gt.eps).or.(cty.eq.275))  then
                at(rs)%npol = at(rs)%npol + 1
                at(rs)%pol(at(rs)%npol) = ii
                atq(ii) = c_charge(cty)
            end if
        end do
!
      else
!
!       exit for any other residue
        write(ilog,*) 'Fatal. Polar potential does not support resid&
 &ue ',resname,' (?) at the moment. Check back later.'
        call fexit() 
!
      end if
!
    end do !over residues in a molecule
!
  end do !over molecules
! Martin : Here is the mod
! This goes with the implicit assumption that only one group per residue will be changed 
!
! apply primary amide corrections if wanted
  if (primamide_cis_chgshft.ne.0.0) then
!
    do imol=1,nmol
!
      do rs=rsmol(imol,1),rsmol(imol,2)
!
        resname = amino(seqtyp(rs))
        if (resname.eq.'ASN') then
          do j=1,at(rs)%npol
            if (at(rs)%pol(j).eq.at(rs)%sc(8-shf3)) &
 &    atq(at(rs)%sc(8-shf3)) = atq(at(rs)%sc(8-shf3)) - primamide_cis_chgshft
            if (at(rs)%pol(j).eq.at(rs)%sc(9-shf3)) &
 &    atq(at(rs)%sc(9-shf3)) = atq(at(rs)%sc(9-shf3)) + primamide_cis_chgshft
          end do
        end if
        if (resname.eq.'GLN') then
          do j=1,at(rs)%npol
            if (at(rs)%pol(j).eq.at(rs)%sc(11-shf5)) &
 &    atq(at(rs)%sc(11-shf5)) = atq(at(rs)%sc(11-shf5)) - primamide_cis_chgshft
            if (at(rs)%pol(j).eq.at(rs)%sc(12-shf5)) &
 &    atq(at(rs)%sc(12-shf5)) = atq(at(rs)%sc(12-shf5)) + primamide_cis_chgshft
          end do
        end if
        if (resname.eq.'NH2') then
          do j=1,at(rs)%npol
            if (at(rs)%pol(j).eq.at(rs)%bb(2)) atq(at(rs)%bb(2)) = atq(at(rs)%bb(2)) - primamide_cis_chgshft
            if (at(rs)%pol(j).eq.at(rs)%bb(3)) atq(at(rs)%bb(3)) = atq(at(rs)%bb(3)) + primamide_cis_chgshft
          end do
        end if
        if (resname.eq.'ACA') then
          do j=1,at(rs)%npol
            if (at(rs)%pol(j).eq.at(rs)%bb(4)) atq(at(rs)%bb(4)) = atq(at(rs)%bb(4)) - primamide_cis_chgshft
            if (at(rs)%pol(j).eq.at(rs)%bb(5)) atq(at(rs)%bb(5)) = atq(at(rs)%bb(5)) + primamide_cis_chgshft
          end do
        end if
        if (resname.eq.'PPA') then
          do j=1,at(rs)%npol
            if (at(rs)%pol(j).eq.at(rs)%bb(4)) atq(at(rs)%bb(4)) = atq(at(rs)%bb(4)) - primamide_cis_chgshft
            if (at(rs)%pol(j).eq.at(rs)%bb(5)) atq(at(rs)%bb(5)) = atq(at(rs)%bb(5)) + primamide_cis_chgshft
          end do
        end if
        if (resname.eq.'FOA') then
          do j=1,at(rs)%npol
            if (at(rs)%pol(j).eq.at(rs)%bb(5-shf)) &
 &     atq(at(rs)%bb(5-shf)) = atq(at(rs)%bb(5-shf)) - primamide_cis_chgshft
            if (at(rs)%pol(j).eq.at(rs)%bb(6-shf)) &
 &     atq(at(rs)%bb(6-shf)) = atq(at(rs)%bb(6-shf)) + primamide_cis_chgshft
          end do
        end if
!
      end do
    end do
  end if

! apply other patches if wanted
  allocate(tmpvec(n))
  call read_chargepatchfile(tmpvec,k)
! first add all patched charged atoms (if needed)
  do i=1,k
    rs = atmres(tmpvec(i))
    docycle = .false.
    do j=1,at(rs)%npol
      if (tmpvec(i).eq.at(rs)%pol(j)) then
        docycle = .true.
        exit ! found
      end if
    end do
    if (docycle.EQV..false.) then
      do j=1,at(rs)%npol
        if (at(rs)%pol(j).gt.tmpvec(i)) exit ! assume ordered
      end do
      if (j.le.at(rs)%npol) at(rs)%pol((j+1):(at(rs)%npol+1)) = at(rs)%pol(j:(at(rs)%npol))
      at(rs)%npol = at(rs)%npol + 1
      at(rs)%pol(j) = tmpvec(i)
    end if
  end do
! now check for missing parameters (patch is supposed to be able to override these types of errors)
  do i=1,n
    docycle = .false.
    cty = bio_ctyp(b_type(i))
    if (cty.eq.0) then
      do j=1,at(atmres(i))%npol
        if (at(atmres(i))%pol(j).eq.i) then
          docycle = .true.
          exit
        end if
      end do
      if (docycle.EQV..true.) cycle
      write(ilog,*) 'Fatal. Parameter file does not support polar potential for biotype&
 & ',b_type(i),' (',bio_code(b_type(i)),') yet (encountered for atom #',i,').'
      call fexit()
    end if
  end do
  
! lastly, remove (again) the ones with zero charge
  do i=1,k
    rs = atmres(tmpvec(i))
    ! Martin : change : 
    if ((abs(atq(tmpvec(i))).le.eps)) then !.and.(bio_ctyp(b_type(tmpvec(i))).ne.275)) then
      do j=1,at(rs)%npol
        if (tmpvec(i).eq.at(rs)%pol(j)) then
          if (j.lt.at(rs)%npol) at(rs)%pol(j:(at(rs)%npol-1)) = at(rs)%pol((j+1):at(rs)%npol)
          at(rs)%npol = at(rs)%npol - 1
          exit
        end if
      end do
    end if
  end do
  deallocate(tmpvec)
!
!
! populate netchg(1:nseq) (note that this is only used in setup routines)
!
 468  format('Warning. Residue #',i6,' has a net charge of ',f9.4,', which&
 & is generally illegal (fractional or out of range).')
 471  format('Atom #        Charge    Biotype #')
 472  format(i10,1x,f10.5,1x,i6)
  do rs=1,nseq
    foundit = .false.
    tc = 0.0
    do i=1,at(rs)%npol
      tc = tc + atq(at(rs)%pol(i))
    end do 
    do k=-20,20
      if (abs(tc-1.0*k).le.at(rs)%npol*eps) then
        foundit = .true.
        netchg(rs) = k
        exit
      end if
    end do
    if (foundit.EQV..false.) then
      write(ilog,468) rs,tc
      write(ilog,471)
      do i=1,at(rs)%npol
        write(ilog,472) at(rs)%pol(i),atq(at(rs)%pol(i)),b_type((at(rs)%pol(i)))
      end do
      if (tc.lt.0.0) netchg(rs) = floor(tc)
      if (tc.gt.0.0) netchg(rs) = ceiling(tc)
    end if
  end do
!
! get the monopole and dipole groups per residue
  call get_dipgrps()
! Martin : I have a a problem with the dipole groups not   
! now we can set up the interactions arrays accounting for connectivity
!
  allocate(skip(2*(MAXVLNC+MAXVLNC*(MAXVLNC-1)+&
 &              MAXVLNC*(MAXVLNC-1)*(MAXVLNC-1))))
!
  do i=1,n
    if (eqatm(i).lt.0) eqatm(i) = i
  end do
!
  
  allocate(at_temp(nseq))
  
  do rs=1,nseq
!   
    resname = amino(seqtyp(rs))
    nrpolintra(rs) = 0
    allocate(fudge(rs)%elin((at(rs)%npol*(at(rs)%npol-1))/2))
    allocate(iaa(rs)%polin((at(rs)%npol*(at(rs)%npol-1))/2,2))
    do k=1,(at(rs)%npol*(at(rs)%npol-1))/2
      fudge(rs)%elin(k) = 1.0
    end do
    if (rs.lt.nseq) then
      alcsz = at(rs)%npol*at(rs+1)%npol
      if (disulf(rs).gt.0) alcsz = max(alcsz,at(rs)%npol*at(disulf(rs))%npol)
      allocate(fudge(rs)%elnb(alcsz))
      allocate(iaa(rs)%polnb(alcsz,2))
      do k=1,alcsz
        fudge(rs)%elnb(k) = 1.0
      end do
    end if
!
!
!   now it gets a little complicated. on the outside we have to identify what assumptions
!   about rigidity are made (nbsr_model). this filters the available interactions
!   (see below).
!   within, we distinguish between the electrostatic models
!   for the sane model (only full dipoles/charges interacting), we use the automatically
!   determined charge groups to ensure only whole charge groups are interacting.
!   as soon as there is any pair of atoms 13-bonded the group interaction is
!   dismissed
!   if (quasi-)rigidity constraints are to be honored, we need to be able to
!   properly detect effective 13-interactions across rigid systems (like aromatic
!   rings), for which we use equivalence atoms which need to be set up in
!   proteus and sidechain (i.e., the builder) -> see there for details
!   for 14-fudges any group interaction including a 14-bound term will always be
!   14-fudged. for mode_14 = 2, this is extendend to all groups having any pair of
!   atoms being separated by only a single rotatable bond (happens in rigid systems
!   again, and is again handled through equivalence atoms) 
!   note that this function given the charge groups is fully automatic and makes no
!   local assumptions about rigid geometry or charge groups. it will only break
!   if charge groups cannot be described properly or if the use equivalence atoms is
!   impossible (such as for the omega bond in polypeptides).
!   if rigidity constraints are not present or to be ignored, every bond is assumed
!   to be rotatable and the parsing is much simpler (14-only)
!

     ! pKa : Here we have a problem with the dummy charge being assigned its own charged group
!            if ((myrank.eq.0).or.(myrank.eq.5)) then
!          write(ilog,*) "pre_at(rs)%ndpgrps",rs,at(rs)%ndpgrps
!          do i=1,at(rs)%ndpgrps
!              write(ilog,*) "pre_at(rs)%dpgrp(i)%nats",at(rs)%dpgrp(i)%nats
!              do j=1,at(rs)%dpgrp(i)%nats
!                  write(ilog,*) at(rs)%dpgrp(i)%ats(j),atq(at(rs)%dpgrp(i)%ats(j))
!              end do 
!          end do 
!        end if 
  foundit=.false.
  do i=1,at(rs)%ndpgrps
      if (at(rs)%dpgrp(i)%nats.eq.1) then 
          do j=1,at(rs)%ndpgrps
              do k=1,at(rs)%dpgrp(j)%nats
                  ii=at(rs)%dpgrp(i)%ats(1)
                  kk=at(rs)%dpgrp(j)%ats(k)
                  
                  do l=1,n12(kk)
                    if (i12(l,kk).eq.ii) then
                        foundit=.true.
                        exit
                    end if 
                  end do 
                  
                  if (foundit.eqv..false.) then 
                      do l=1,n12(ii)
                        if (i12(l,ii).eq.kk) then
                            foundit=.true.
                            exit
                        end if 
                      end do                      
                  end if    
                  if (foundit.eqv..true.) exit
              end do 
              if (foundit.eqv..true.) exit
          end do 
      end if 
      if (foundit.eqv..true.) exit
  end do 
           
  !At that point, j is the group where the atom is going to be , i is where it is  
  if (foundit.eqv..true.) then 
      ! I am going to create a array with a emtpy dpgrp, in case the dummy is not in the last
      !Then I can jsut skip it during reallocation
      at_temp(rs)%ndpgrps=at(rs)%ndpgrps
      allocate(at_temp(rs)%dpgrp(at_temp(rs)%ndpgrps))
      do l=1,at(rs)%ndpgrps
          if ((l.ne.i).and.(l.ne.j)) then
              at_temp(rs)%dpgrp(l)%nats=at(rs)%dpgrp(l)%nats
              allocate(at_temp(rs)%dpgrp(l)%ats(at(rs)%dpgrp(l)%nats))
              at_temp(rs)%dpgrp(l)%ats(:)=at(rs)%dpgrp(l)%ats(:)
              !at_temp(rs)%dpgrp(l)%cgn=at(rs)%dpgrp(l)%cgn
          end if
      end do 

      at_temp(rs)%dpgrp(j)%nats=at(rs)%dpgrp(j)%nats+1
      at_temp(rs)%dpgrp(i)%nats=at(rs)%dpgrp(i)%nats-1

      allocate(at_temp(rs)%dpgrp(j)%ats(at_temp(rs)%dpgrp(j)%nats))

      at_temp(rs)%dpgrp(j)%ats(1:at(rs)%dpgrp(j)%nats)=at(rs)%dpgrp(j)%ats(:)
      at_temp(rs)%dpgrp(j)%ats(at(rs)%dpgrp(j)%nats+1)=at(rs)%dpgrp(i)%ats(1)
      !at_temp(rs)%dpgrp(l)%cgn=at(rs)%dpgrp(l)%cgn


      deallocate(at(rs)%dpgrp)
      at(rs)%ndpgrps=at(rs)%ndpgrps-1
      
      
      allocate(at(rs)%dpgrp(at(rs)%ndpgrps))
      ! Can start reusing the i,j,k here
      offset=0
      do i=1,at_temp(rs)%ndpgrps
          if (at_temp(rs)%dpgrp(i)%nats.ne.0) then 
            at(rs)%dpgrp(i-offset)%nats=at_temp(rs)%dpgrp(i)%nats
            allocate(at(rs)%dpgrp(i-offset)%ats(at(rs)%dpgrp(i-offset)%nats))
            at(rs)%dpgrp(i-offset)%ats(:)=at_temp(rs)%dpgrp(i)%ats(:)

            at(rs)%dpgrp(i-offset)%nc=int(sum(atq(at(rs)%dpgrp(i-offset)%ats(:))))
            at(rs)%dpgrp(i-offset)%tc=sum(atq(at(rs)%dpgrp(i-offset)%ats(:)))
            at(rs)%dpgrp(i-offset)%qnm=sum(abs(atq(at(rs)%dpgrp(i-offset)%ats(:))))
            
            !at(rs)%dpgrp(i-offset)%cgn=at_temp(rs)%dpgrp(i-offset)%cgn
          else 
              offset=offset+1
          end if 
      end do 

    end if 

    if (nbsr_model.eq.1) then
      if (elec_model.eq.2) then !Martin : That's us 
!        if ((myrank.eq.0).or.(myrank.eq.5)) then
!          write(ilog,*) "at(rs)%ndpgrps",rs,at(rs)%ndpgrps
!          do i=1,at(rs)%ndpgrps
!              write(ilog,*) "at(rs)%dpgrp(i)%nats",at(rs)%dpgrp(i)%nats
!              do j=1,at(rs)%dpgrp(i)%nats
!                  write(ilog,*) at(rs)%dpgrp(i)%ats(j),atq(at(rs)%dpgrp(i)%ats(j))
!              end do 
!          end do 
!        end if 
        do dp1=1,at(rs)%ndpgrps
          do dp2=dp1+1,at(rs)%ndpgrps
            proceed = .true.
            its14 = .false.
            do i=1,at(rs)%dpgrp(dp1)%nats
              if (proceed.EQV..false.) exit
              ii = at(rs)%dpgrp(dp1)%ats(i)
!
              do k=1,at(rs)%dpgrp(dp2)%nats
                
                if (proceed.EQV..false.) exit
                kk = at(rs)%dpgrp(dp2)%ats(k)

                if (seqtyp(rs).eq.26) then
!
                  if ((abs(atq(ii)).le.eps).OR.(abs(atq(kk)).le.eps)) cycle ! Martin : Should I put one here as well ?
!
!                 standard topology checks and rotation list-based checks
                  call ia_rotlsts(ii,kk,proceed,its14)
!
                  if ((its14.EQV..true.).AND.(use_14.EQV..false.)) proceed = .false.
!
                  if (proceed.EQV..false.) exit
!
                else
!
                  if ((eqatm(ii).eq.eqatm(kk)).OR.&
                    &(eqatm(ii).eq.kk).OR.(eqatm(kk).eq.ii)) proceed = .false.
!
                  do j=1,n12(ii)
                    if (i12(j,ii).eq.kk) proceed = .false.
                    if (i12(j,ii).eq.eqatm(kk)) proceed = .false.
                  end do
                  if (proceed.EQV..false.) exit
                  if (eqatm(ii).ne.ii) then
                    do j=1,n12(eqatm(ii))
                      if (i12(j,eqatm(ii)).eq.kk) proceed = .false.
                    if (i12(j,eqatm(ii)).eq.eqatm(kk)) proceed=.false.
                    end do
                    if (proceed.EQV..false.) exit
                  end if
!
                  do j=1,n13(ii)
                    if (i13(j,ii).eq.kk) proceed = .false.
                    if (i13(j,ii).eq.eqatm(kk)) proceed = .false.
                  end do
                  if (proceed.EQV..false.) exit
                  if (eqatm(ii).ne.ii) then
                    do j=1,n13(eqatm(ii))
                      if (i13(j,eqatm(ii)).eq.kk) proceed = .false.
                    if (i13(j,eqatm(ii)).eq.eqatm(kk)) proceed=.false.
                    end do
                    if (proceed.EQV..false.) exit
                  end if
!
                  do j=1,n14(ii)
                    if (use_14.EQV..false.) then
                      if (i14(j,ii).eq.kk) proceed = .false.
                      if (i14(j,ii).eq.eqatm(kk)) proceed = .false.
                    else
                      if (i14(j,ii).eq.kk) its14 = .true.
                      if (mode_14.eq.2) then
                        if (i14(j,ii).eq.eqatm(kk)) its14 = .true.
                      end if
                    end if
                  end do
                  if (proceed.EQV..false.) exit
                  if (eqatm(ii).ne.ii) then
                    do j=1,n14(eqatm(ii))
                      if (use_14.EQV..false.) then
                        if (i14(j,eqatm(ii)).eq.kk) proceed = .false.
                      if(i14(j,eqatm(ii)).eq.eqatm(kk))proceed=.false.
                      else
                        if (mode_14.eq.2) then
                          if (i14(j,ii).eq.eqatm(kk)) its14 = .true.
                        if(i14(j,eqatm(ii)).eq.eqatm(kk)) its14=.true.
                        end if
                      end if
                    end do
                    if (proceed.EQV..false.) exit
                  end if
                end if
              end do
            end do
            ! Martin : pka: here things are modified only in the part that correspond to our situation. Use of different model
            !would request a rewritting
            if (proceed.EQV..true.) then
              do i=1,at(rs)%dpgrp(dp1)%nats
                ii = at(rs)%dpgrp(dp1)%ats(i)
                do k=1,at(rs)%dpgrp(dp2)%nats
                  kk = at(rs)%dpgrp(dp2)%ats(k)
                  test1=.false.
                  test2=.false.
                  if ((seqtyp(rs).eq.26).AND.((abs(atq(ii)).le.eps).OR.(abs(atq(kk)).le.eps)))  cycle 
                  nrpolintra(rs) = nrpolintra(rs) + 1 
                  iaa(rs)%polin(nrpolintra(rs),1) = ii
                  iaa(rs)%polin(nrpolintra(rs),2) = kk
                  if (its14.EQV..true.) then
                    fudge(rs)%elin(nrpolintra(rs)) = fudge_el_14
                  end if
                end do
              end do
            end if
          end do
        end do
        
!     in the blind model we handle rigid geometries analogously, but don't
!     respect charge groups whatsoever
!
      else if (elec_model.eq.1) then
!
        do i=1,at(rs)%npol
          ii = at(rs)%pol(i)
          do k=i+1,at(rs)%npol
            kk = at(rs)%pol(k)
            proceed = .true.
            its14 = .false.
!
            if (seqtyp(rs).eq.26) then
              if ((abs(atq(ii)).le.eps).OR.(abs(atq(kk)).le.eps)) cycle
!
!             standard topology checks and rotation list-based checks
              call ia_rotlsts(ii,kk,proceed,its14)
!
              if ((its14.EQV..true.).AND.(use_14.EQV..false.)) proceed = .false.
!
            else
!
              if ((eqatm(ii).eq.eqatm(kk)).OR.&
   &   (eqatm(ii).eq.kk).OR.(eqatm(kk).eq.ii)) proceed = .false.
!
              do j=1,n12(ii)
                if (i12(j,ii).eq.kk) proceed = .false.
                if (i12(j,ii).eq.eqatm(kk)) proceed = .false.
              end do
              if (eqatm(ii).ne.ii) then
                do j=1,n12(eqatm(ii))
                  if (i12(j,eqatm(ii)).eq.kk) proceed = .false.
                if (i12(j,eqatm(ii)).eq.eqatm(kk)) proceed=.false.
                end do
              end if
!
              do j=1,n13(ii)
                if (i13(j,ii).eq.kk) proceed = .false.
                if (i13(j,ii).eq.eqatm(kk)) proceed = .false.
              end do
              if (eqatm(ii).ne.ii) then
                do j=1,n13(eqatm(ii))
                  if (i13(j,eqatm(ii)).eq.kk) proceed = .false.
                if (i13(j,eqatm(ii)).eq.eqatm(kk)) proceed = .false.
                end do
              end if
!
              do j=1,n14(ii)
                if (use_14.EQV..false.) then
                  if (i14(j,ii).eq.kk) proceed = .false.
                  if (i14(j,ii).eq.eqatm(kk)) proceed = .false.
                else
                  if (i14(j,ii).eq.kk) its14 = .true.
                  if (mode_14.eq.2) then
                    if (i14(j,ii).eq.eqatm(kk)) its14 = .true.
                  end if
                end if
              end do
              if (eqatm(ii).ne.ii) then
                do j=1,n14(eqatm(ii))
                  if (use_14.EQV..false.) then
                    if (i14(j,eqatm(ii)).eq.kk) proceed = .false.
                    if (i14(j,eqatm(ii)).eq.eqatm(kk))proceed =.false.
                  else
                    if (mode_14.eq.2) then
                      if (i14(j,ii).eq.eqatm(kk)) its14 = .true.
                      if(i14(j,eqatm(ii)).eq.eqatm(kk)) its14 = .true.
                    end if
                  end if
                end do
              end if
            end if
!
            if (proceed.EQV..true.) then
              nrpolintra(rs) = nrpolintra(rs) + 1 
              iaa(rs)%polin(nrpolintra(rs),1) = ii
              iaa(rs)%polin(nrpolintra(rs),2) = kk
              if (its14.EQV..true.) then
                fudge(rs)%elin(nrpolintra(rs)) = fudge_el_14
              end if
            end if
          end do
        end do
!
      else
!
        write(ilog,*) 'Fatal. Encountered unknown electrostatic mode&
 &l in polar_groups(). This is most certainly a bug.'
        call fexit()
!
      end if
!
!   in this mode we assume every bond is rotatable
!
    else if ((nbsr_model.eq.2).OR.(nbsr_model.eq.3)) then
!
      if (elec_model.eq.2) then
        do dp1=1,at(rs)%ndpgrps
          do dp2=dp1+1,at(rs)%ndpgrps
            proceed = .true.
            its14 = .false.
            do i=1,at(rs)%dpgrp(dp1)%nats
              if (proceed.EQV..false.) exit
              ii = at(rs)%dpgrp(dp1)%ats(i)
!
              do k=1,at(rs)%dpgrp(dp2)%nats
                if (proceed.EQV..false.) exit
                kk = at(rs)%dpgrp(dp2)%ats(k)
!
                if ((seqtyp(rs).eq.26).AND.((abs(atq(ii)).le.eps).OR.(abs(atq(kk)).le.eps))) cycle
                do j=1,n12(ii)
                  if (i12(j,ii).eq.kk) proceed = .false.
                end do
                if (proceed.EQV..false.) exit
!
                do j=1,n13(ii)
                  if (i13(j,ii).eq.kk) proceed = .false.
                end do
                if (proceed.EQV..false.) exit
!
                do j=1,n14(ii)
                  if (use_14.EQV..false.) then
                    if (i14(j,ii).eq.kk) proceed = .false.
                  else
                    if (i14(j,ii).eq.kk) its14 = .true.
                  end if
                end do
                if (proceed.EQV..false.) exit
              end do
            end do
!
            if (proceed.EQV..true.) then
              do i=1,at(rs)%dpgrp(dp1)%nats
                ii = at(rs)%dpgrp(dp1)%ats(i)
                do k=1,at(rs)%dpgrp(dp2)%nats
                  kk = at(rs)%dpgrp(dp2)%ats(k)
                  if ((seqtyp(rs).eq.26).AND.((abs(atq(ii)).le.eps).OR.(abs(atq(kk)).le.eps))) cycle
                  nrpolintra(rs) = nrpolintra(rs) + 1 
                  iaa(rs)%polin(nrpolintra(rs),1) = ii
                  iaa(rs)%polin(nrpolintra(rs),2) = kk
                  if (its14.EQV..true.) then
                    fudge(rs)%elin(nrpolintra(rs)) = fudge_el_14
                  end if
                end do
              end do
            end if
          end do
        end do

      else if (elec_model.eq.1) then
!
        do i=1,at(rs)%npol
          ii = at(rs)%pol(i)
          do k=i+1,at(rs)%npol
            kk = at(rs)%pol(k)
            if ((seqtyp(rs).eq.26).AND.((abs(atq(ii)).le.eps).OR.(abs(atq(kk)).le.eps))) cycle
            proceed = .true.
            its14 = .false.
!
            do j=1,n12(ii)
              if (i12(j,ii).eq.kk) proceed = .false.
            end do
!
            do j=1,n13(ii)
              if (i13(j,ii).eq.kk) proceed = .false.
            end do
!
            do j=1,n14(ii)
              if (use_14.EQV..false.) then
                if (i14(j,ii).eq.kk) proceed = .false.
              else
                if (i14(j,ii).eq.kk) its14 = .true.
              end if
            end do
!
            if (proceed.EQV..true.) then
              nrpolintra(rs) = nrpolintra(rs) + 1 
              iaa(rs)%polin(nrpolintra(rs),1) = ii
              iaa(rs)%polin(nrpolintra(rs),2) = kk
              if (its14.EQV..true.) then
                fudge(rs)%elin(nrpolintra(rs)) = fudge_el_14
              end if
            end if
          end do
        end do
      else
!
        write(ilog,*) 'Fatal. Encountered unknown electrostatic mode&
 &l in polar_groups(). This is most certainly a bug.'
        call fexit()
!
      end if
!
    end if
    

!   cycle out if last residue in sequence
    if (rs.eq.nseq) cycle
!
!   if the residue is the last residue in a molecule, simply add
!   all possible interaction to the next one
    imol = molofrs(rs)
    

    
    if ((rs.eq.rsmol(imol,2)).AND.(rsmol(imol,2).lt.nseq)) then
      nrpolnb(rs) = 0
      resname = amino(seqtyp(rs))
      resname1 = amino(seqtyp(rs+1))
      do i=1,at(rs)%npol
        ii = at(rs)%pol(i)
        if (attyp(ii) .le. 0)  cycle
        do k=1,at(rs+1)%npol
          kk = at(rs+1)%pol(k)
          if (attyp(kk) .le. 0)  cycle
          nrpolnb(rs) = nrpolnb(rs) + 1 
          iaa(rs)%polnb(nrpolnb(rs),1) = ii
          iaa(rs)%polnb(nrpolnb(rs),2) = kk
        end do
      end do
      cycle
    end if
!
!   in the connected case, the situation is more complicated 
    nrpolnb(rs) = 0
    resname = amino(seqtyp(rs))
    resname1 = amino(seqtyp(rs+1))
!
!   similarly to the intra-residue case: if rigidity constraints are to be honored:
!   for the sane, group-based model we use an identical strategy with the
!   exception that equivalence atoms don't work across residues, which is to say
!   we cannot treat the omega bond as rigid through equivalence atoms
!   (remember that regardless of omega sampling, we'll always assume the
!   omega bond to be pseudo-rigid and let the torsional potential handle it)
!   equivalence atoms will still be used for the sidechain terms. interactions
!   involving CA(i)-CA(i+1) will be explicitly dismissed, and 15-interactions
!   spanning the omega bond as a potential rotatable bond will be assumed to be 14
!   instead (like N(i)-CA(i+1) or CA(i)-CB(i+1)) if mode_14 = 2
!   different minor fixes might be needed if the basic structure of the backbone
!   changes, i.e., if one goes to polyesters, polynucleotides, ...
!   if rigidity constraints are not present or to be ignored, every bond is assumed
!   to be rotatable and the parsing is much simpler (14-only)

    if (nbsr_model.eq.1) then
      if (elec_model.eq.2) then !MArtin that is our case
        do dp1=1,at(rs)%ndpgrps
          do dp2=1,at(rs+1)%ndpgrps
            proceed = .true.
            its14 = .false.
            !!!!!!!!!!!!!
            do i=1,at(rs)%dpgrp(dp1)%nats ! Martin : so for each polar atoms pairs combinations in the dipolar group, if only one has  
                                          ! a proceed equal false, nothing in the array will be populated
              if (proceed.EQV..false.) exit
              ii = at(rs)%dpgrp(dp1)%ats(i)
!
              do k=1,at(rs+1)%dpgrp(dp2)%nats
                if (proceed.EQV..false.) exit
                kk = at(rs+1)%dpgrp(dp2)%ats(k)
!
                if (((seqtyp(rs).eq.26).OR.(seqtyp(rs+1).eq.26)).AND.(molofrs(rs).eq.molofrs(rs+1))) then
!
                  if ((abs(atq(ii)).le.eps).OR.(abs(atq(kk)).le.eps)) cycle
!
!                 standard topology checks and rotation list-based checks
                  call ia_rotlsts(ii,kk,proceed,its14)
!
                  if ((its14.EQV..true.).AND.(use_14.EQV..false.)) proceed = .false.

                  if (proceed.EQV..false.) exit
!
                else
!
                  if ((eqatm(ii).eq.eqatm(kk)).OR.&
   &       (eqatm(ii).eq.kk).OR.(eqatm(kk).eq.ii)) proceed = .false.

                  do j=1,n12(ii)
                    if (i12(j,ii).eq.kk) proceed = .false.
                    if (i12(j,ii).eq.eqatm(kk)) proceed = .false.
                  end do
                  if (proceed.EQV..false.) exit
                  if (eqatm(ii).ne.ii) then
                    do j=1,n12(eqatm(ii))
                      if (i12(j,eqatm(ii)).eq.kk) proceed = .false.
                    if (i12(j,eqatm(ii)).eq.eqatm(kk)) proceed=.false.
                    end do
                    if (proceed.EQV..false.) exit
                  end if
                  
                  do j=1,n13(ii)
                    if (i13(j,ii).eq.kk) proceed = .false.
                    if (i13(j,ii).eq.eqatm(kk)) proceed = .false.
                  end do
                  if (proceed.EQV..false.) exit
                  if (eqatm(ii).ne.ii) then
                    do j=1,n13(eqatm(ii))
                      if (i13(j,eqatm(ii)).eq.kk) proceed = .false.
                    if (i13(j,eqatm(ii)).eq.eqatm(kk)) proceed=.false.
                    end do
                    if (proceed.EQV..false.) exit
                  end if
                  
                  do j=1,n14(ii)
                    if (use_14.EQV..false.) then
                      if (i14(j,ii).eq.kk) proceed = .false.
                      if (i14(j,ii).eq.eqatm(kk)) proceed = .false.
                    else
                      if (i14(j,ii).eq.kk) its14 = .true.
                      if (mode_14.eq.2) then
                        if (i14(j,ii).eq.eqatm(kk)) its14 = .true.
                      end if
                    end if
                  end do
                  if (proceed.EQV..false.) exit
                  
                  if (eqatm(ii).ne.ii) then
                    do j=1,n14(eqatm(ii))
                      if (use_14.EQV..false.) then
                        if (i14(j,eqatm(ii)).eq.kk) proceed = .false.
                      if(i14(j,eqatm(ii)).eq.eqatm(kk))proceed=.false.
                      else
                        if (mode_14.eq.2) then
                          if (i14(j,ii).eq.eqatm(kk)) its14 = .true.
                        if(i14(j,eqatm(ii)).eq.eqatm(kk)) its14=.true.
                        end if
                      end if
                    end do
                    if (proceed.EQV..false.) exit
                  end if
!
                end if
               
              end do
            end do
            !!!!!!!!!!!!!
!
            if (proceed.EQV..true.) then
              do i=1,at(rs)%dpgrp(dp1)%nats
                ii = at(rs)%dpgrp(dp1)%ats(i)
                do k=1,at(rs+1)%dpgrp(dp2)%nats
                  kk = at(rs+1)%dpgrp(dp2)%ats(k)
                  if (((seqtyp(rs).eq.26).AND.(abs(atq(ii)).le.eps)).OR.&
 &                  ((seqtyp(rs+1).eq.26).AND.(abs(atq(kk)).le.eps))) cycle
                  nrpolnb(rs) = nrpolnb(rs) + 1 
                  iaa(rs)%polnb(nrpolnb(rs),1) = ii
                  iaa(rs)%polnb(nrpolnb(rs),2) = kk
                  if (its14.EQV..true.) then
                    fudge(rs)%elnb(nrpolnb(rs)) = fudge_el_14
                  end if
                end do
              end do
            end if
          end do
        end do
      else if (elec_model.eq.1) then
        do i=1,at(rs)%npol
          ii = at(rs)%pol(i)
          do k=1,at(rs+1)%npol
            kk = at(rs+1)%pol(k)
            proceed = .true.
            its14 = .false.
!
            if (((seqtyp(rs).eq.26).OR.(seqtyp(rs+1).eq.26)).AND.(molofrs(rs).eq.molofrs(rs+1))) then
!
              if ((abs(atq(ii)).le.eps).OR.(abs(atq(kk)).le.eps)) cycle
!
!             standard topology checks and rotation list-based checks
              call ia_rotlsts(ii,kk,proceed,its14)
!
              if ((its14.EQV..true.).AND.(use_14.EQV..false.)) proceed = .false.
!
            else
!
              if ((eqatm(ii).eq.eqatm(kk)).OR.&
   &       (eqatm(ii).eq.kk).OR.(eqatm(kk).eq.ii)) proceed = .false.
!
              do j=1,n12(ii)
                if (i12(j,ii).eq.kk) proceed = .false.
                if (i12(j,ii).eq.eqatm(kk)) proceed = .false.
              end do
              if (eqatm(ii).ne.ii) then
                do j=1,n12(eqatm(ii))
                  if (i12(j,eqatm(ii)).eq.kk) proceed = .false.
                  if (i12(j,eqatm(ii)).eq.eqatm(kk)) proceed=.false.
                end do
              end if
!
              do j=1,n13(ii)
                if (i13(j,ii).eq.kk) proceed = .false.
                if (i13(j,ii).eq.eqatm(kk)) proceed = .false.
              end do
              if (eqatm(ii).ne.ii) then
                do j=1,n13(eqatm(ii))
                  if (i13(j,eqatm(ii)).eq.kk) proceed = .false.
                  if (i13(j,eqatm(ii)).eq.eqatm(kk)) proceed= .false.
                end do
              end if
!
              do j=1,n14(ii)
                if (use_14.EQV..false.) then
                  if (i14(j,ii).eq.kk) proceed = .false.
                  if (i14(j,ii).eq.eqatm(kk)) proceed = .false.
                else
                  if (i14(j,ii).eq.kk) its14 = .true.
                  if (mode_14.eq.2) then
                    if (i14(j,ii).eq.eqatm(kk)) its14 = .true.
                  end if
                end if
              end do
!
              if (eqatm(ii).ne.ii) then
                do j=1,n14(eqatm(ii))
                  if (use_14.EQV..false.) then
                    if (i14(j,eqatm(ii)).eq.kk) proceed = .false.
                  if (i14(j,eqatm(ii)).eq.eqatm(kk)) proceed = .false.
                  else
                    if (mode_14.eq.2) then
                      if (i14(j,ii).eq.eqatm(kk)) its14 = .true.
                    if (i14(j,eqatm(ii)).eq.eqatm(kk)) its14= .true.
                    end if
                  end if
                end do
              end if
            end if
!
            if (proceed.EQV..true.) then
              nrpolnb(rs) = nrpolnb(rs) + 1 
              iaa(rs)%polnb(nrpolnb(rs),1) = ii
              iaa(rs)%polnb(nrpolnb(rs),2) = kk
              if (its14.EQV..true.) then
                fudge(rs)%elnb(nrpolnb(rs)) = fudge_el_14
              end if
            end if
          end do
        end do
!
      end if
!
    else if ((nbsr_model.eq.2).OR.(nbsr_model.eq.3)) then
!
      if (elec_model.eq.2) then
        do dp1=1,at(rs)%ndpgrps
          do dp2=1,at(rs+1)%ndpgrps
            proceed = .true.
            its14 = .false.
            do i=1,at(rs)%dpgrp(dp1)%nats
              if (proceed.EQV..false.) exit
              ii = at(rs)%dpgrp(dp1)%ats(i)
!
              do k=1,at(rs+1)%dpgrp(dp2)%nats
                if (proceed.EQV..false.) exit
                kk = at(rs+1)%dpgrp(dp2)%ats(k)
                if (((seqtyp(rs).eq.26).AND.(abs(atq(ii)).le.eps)).OR.&
 &                  ((seqtyp(rs+1).eq.26).AND.(abs(atq(kk)).le.eps))) cycle
!
                do j=1,n12(ii)
                  if (i12(j,ii).eq.kk) proceed = .false.
                end do
                if (proceed.EQV..false.) exit
!
                do j=1,n13(ii)
                  if (i13(j,ii).eq.kk) proceed = .false.
                end do
                if (proceed.EQV..false.) exit
!
                do j=1,n14(ii)
                  if (use_14.EQV..false.) then
                    if (i14(j,ii).eq.kk) proceed = .false.
                  else
                    if (i14(j,ii).eq.kk) its14 = .true.
                  end if
                end do
                if (proceed.EQV..false.) exit
!
              end do
            end do
!
            if (proceed.EQV..true.) then
              do i=1,at(rs)%dpgrp(dp1)%nats
                ii = at(rs)%dpgrp(dp1)%ats(i)
                do k=1,at(rs+1)%dpgrp(dp2)%nats
                  kk = at(rs+1)%dpgrp(dp2)%ats(k)
                  if (((seqtyp(rs).eq.26).AND.(abs(atq(ii)).le.eps)).OR.&
 &                  ((seqtyp(rs+1).eq.26).AND.(abs(atq(kk)).le.eps))) cycle
                  nrpolnb(rs) = nrpolnb(rs) + 1 
                  iaa(rs)%polnb(nrpolnb(rs),1) = ii
                  iaa(rs)%polnb(nrpolnb(rs),2) = kk
                  if (its14.EQV..true.) then
                    fudge(rs)%elnb(nrpolnb(rs)) = fudge_el_14
                  end if
                end do
              end do
            end if
          end do
        end do
!
      else if (elec_model.eq.1) then
        do i=1,at(rs)%npol
          ii = at(rs)%pol(i)
          do k=1,at(rs+1)%npol
            kk = at(rs+1)%pol(k)
            if (((seqtyp(rs).eq.26).AND.(abs(atq(ii)).le.eps)).OR.&
 &              ((seqtyp(rs+1).eq.26).AND.(abs(atq(kk)).le.eps))) cycle
            proceed = .true.
            its14 = .false.
!
            do j=1,n12(ii)
              if (i12(j,ii).eq.kk) proceed = .false.
            end do
!
            do j=1,n13(ii)
              if (i13(j,ii).eq.kk) proceed = .false.
            end do
!
            do j=1,n14(ii)
              if (use_14.EQV..false.) then
                if (i14(j,ii).eq.kk) proceed = .false.
              else
                if (i14(j,ii).eq.kk) its14 = .true.
              end if
            end do
!
            if (proceed.EQV..true.) then
              nrpolnb(rs) = nrpolnb(rs) + 1 
              iaa(rs)%polnb(nrpolnb(rs),1) = ii
              iaa(rs)%polnb(nrpolnb(rs),2) = kk
              if (its14.EQV..true.) then
                fudge(rs)%elnb(nrpolnb(rs)) = fudge_el_14
              end if
            end if
          end do
        end do
!
      end if
!    
    end if
!
  end do

  deallocate(skip)
  deallocate(at_temp)
  if (nbsr_model.eq.3) then
    call GROMOS_excludes_polar()
  end if
  if (n_crosslinks.gt.0) then
    call crosslink_excludes_polar()
  end if
!
! finally populate the complementary list, i.e., the list of EXcluded
! polar interactions (used for Ewald sums, e.g.)

  do imol=1,nmol
! 
    do rs=rsmol(imol,1),rsmol(imol,2)
!
      nrexpolin(rs) = 0
      nrexpolnb(rs) = 0
!
      if (natres(rs).gt.1) then
        allocate(iaa(rs)%expolin((natres(rs)*(natres(rs)-1))/2,2))
      end if
      if (rs.lt.nseq) then
        alcsz = natres(rs)*natres(rs+1)
        if (disulf(rs).gt.0) alcsz = max(alcsz,natres(rs)*natres(disulf(rs)))
        allocate(iaa(rs)%expolnb(alcsz,2))
      end if
!
      do i=1,at(rs)%npol
        ii = at(rs)%pol(i)
        do k=i+1,at(rs)%npol
          kk = at(rs)%pol(k)
          proceed = .true. 
          do j=1,nrpolintra(rs)
            if (((iaa(rs)%polin(j,1).eq.ii).AND.&
            &(iaa(rs)%polin(j,2).eq.kk)).OR.&
            & ((iaa(rs)%polin(j,1).eq.kk).AND.(iaa(rs)%polin(j,2).eq.ii))) then
              proceed = .false.
              exit
            end if
          end do
          if (proceed.EQV..true.) then
            nrexpolin(rs) = nrexpolin(rs) + 1
            iaa(rs)%expolin(nrexpolin(rs),1) = ii
            iaa(rs)%expolin(nrexpolin(rs),2) = kk
          end if
        end do
      end do
!
      if (rs.eq.nseq) cycle
!
      do i=1,at(rs)%npol
        ii = at(rs)%pol(i)
        do k=1,at(rs+1)%npol
          kk = at(rs+1)%pol(k)
          proceed = .true. 
          do j=1,nrpolnb(rs)
            if (((iaa(rs)%polnb(j,1).eq.ii).AND.&
 &               (iaa(rs)%polnb(j,2).eq.kk)).OR.&
 & ((iaa(rs)%polnb(j,1).eq.kk).AND.(iaa(rs)%polnb(j,2).eq.ii))) then
              proceed = .false.
              exit
            end if
          end do
          if (proceed.EQV..true.) then
            nrexpolnb(rs) = nrexpolnb(rs) + 1
            iaa(rs)%expolnb(nrexpolnb(rs),1) = ii
            iaa(rs)%expolnb(nrexpolnb(rs),2) = kk
          end if
        end do
      end do
    end do
  end do

  do k=1,nmoltyp
    do rs=rsmol(moltyp(k,1),1),rsmol(moltyp(k,1),2)
      if (at(rs)%npol.gt.0) molischg(k) = .true.
    end do
  end do

  tc = 0.0
  npcat = 0
  npan = 0
  nccat = 0
  ncan = 0
  do imol=1,nmol
    if ((moltermid(imol,1).eq.1).AND.&
 &      (seqpolty(rsmol(imol,1)).eq.'P')) then
      npcat = npcat + 1
    end if
    if ((moltermid(imol,2).eq.1).AND.&
 &      (seqpolty(rsmol(imol,2)).eq.'P')) then
      npan = npan + 1
    end if
    do rs=rsmol(imol,1),rsmol(imol,2)
      resname = amino(seqtyp(rs))
      do i=1,at(rs)%npol
        tc = tc + atq(at(rs)%pol(i))
      end do
      if ((resname.eq.'ASP').OR.(resname.eq.'GLU').OR.&
 &        (resname.eq.'TYO').OR.(resname.eq.'CYX').OR.&
 &        (resname.eq.'S1P').OR.(resname.eq.'T1P').OR.&
 &        (resname.eq.'RPC').OR.(resname.eq.'DPC').OR.&
 &        (resname.eq.'RPU').OR.(resname.eq.'DPU').OR.&
 &        (resname.eq.'RPT').OR.(resname.eq.'DPT').OR.&
 &        (resname.eq.'RPA').OR.(resname.eq.'DPA').OR.&
 &        (resname.eq.'RPG').OR.(resname.eq.'DPG').OR.&
 &        (resname.eq.'R5P').OR.(resname.eq.'D5P').OR.&! MArtin : Added this in case or I did not...
 &        (resname.eq.'Y1P')) then
        npan = npan + 1
      else if ((resname.eq.'LYS').OR.(resname.eq.'ARG').OR.&
 &             (resname.eq.'ORN').OR.(resname.eq.'DAB').OR.&
 &             (resname.eq.'HIP')) then
        npcat = npcat + 1
      else if ((resname.eq.'NA+').OR.(resname.eq.'K+ ').OR.&
 &             (resname.eq.'CS+').OR.(resname.eq.'NH4').OR.&
 &             (resname.eq.'1MN').OR.(resname.eq.'2MN').OR.&
 &             (resname.eq.'GDN')) then
        nccat = nccat + 1
      else if ((resname.eq.'CL-').OR.(resname.eq.'BR-').OR.&
 &             (resname.eq.'I- ').OR.(resname.eq.'AC-').OR.&
 &             (resname.eq.'LCP').OR.(resname.eq.'NO3').OR.&
 &             (resname.eq.'CLX')) then
        ncan = ncan + 1
      
      end if
    end do
  end do

  if (abs(tc).gt.0.0001) then
      if (do_hsq.eqv..true.) then 
          if (abs(tc)-TOT_QX.gt.0.0001) then
              print *,"Doing HSQ with a net charge of ",TOT_QX
          end if 
      else
            write(ilog,*) 'Fatal. System is not net neutral (total charge is&
         &',tc,' units). Check sequence and add appropriate counterions.'
            if (be_unsafe.EQV..false.) call fexit()
      end if 
  end if
  if (npcat.gt.ncan) then
    write(ilog,*) 'Warning. Usually it is recommended to add at leas&
 &t one mobile counterion for every macromolecular charged group (he&
 &re: ',npcat-ncan,' units of surplus cation charge).'
    write(ilog,*)
  end if
  if (npan.gt.nccat) then
    write(ilog,*) 'Warning. Usually it is recommended to add at leas&
 &t one mobile counterion for every macromolecular charged group (he&
 &re: ',npan-nccat,' units of surplus anion charge).'
    write(ilog,*)
  end if
!
  if (elec_report .EQV..true.) then
 44     format('Residue ',i3,' (',a3,') has ',i3,' polar atoms:')
 43     format('Its total charge is ',f6.3,'.')
 45     format(i6,' (',a4,'): ',f9.6)
 46     format(i7,' ',i7,' ',f6.3)
 47     format('There are ',i4,' intra-residue interactions:')
 48     format('There are ',i4,' to-next-neighbor interactions:')
 49     format('Atom #1 Atom #2 Fudge-factor')
    write(ilog,*)
    write(ilog,*) '---Summary of polar interactions ---'
    write(ilog,*)
    do rs=1,nseq
      resname = amino(seqtyp(rs))
      write(ilog,44) rs,resname,at(rs)%npol
      tc = 0.0
      do i=1,at(rs)%npol
        write(ilog,45) at(rs)%pol(i),bio_code(b_type(at(rs)%pol(i)))&
 &                                      ,atq(at(rs)%pol(i))
        tc = tc + atq(at(rs)%pol(i))
      end do
      write(ilog,43) tc
      write(ilog,*)
      write(ilog,47) nrpolintra(rs)
      write(ilog,49)
      do i=1,nrpolintra(rs)
        write(ilog,46) iaa(rs)%polin(i,1),iaa(rs)%polin(i,2),&
 &fudge(rs)%elin(i)
      end do
!      do i=1,nrsintra(rs)
!        write(ilog,46) iaa(rs)%atin(i,1),iaa(rs)%atin(i,2),
! &fudge(rs)%rsin(i)
!      end do
      write(ilog,*)
      if (rs.lt.nseq) then
      write(ilog,48)nrpolnb(rs)
      write(ilog,49)
      do i=1,nrpolnb(rs)
        write(ilog,46) iaa(rs)%polnb(i,1),iaa(rs)%polnb(i,2),&
 &fudge(rs)%elnb(i)
      end do
!      do i=1,nrsnb(rs)
!        write(ilog,46) iaa(rs)%atnb(i,1),iaa(rs)%atnb(i,2),
! &fudge(rs)%rsnb(i)
!      end do
      write(ilog,*)
      end if
    end do
  end if
! finally set Coulomb screening factor if necessary
  if (use_IMPSOLV.EQV..true.) then
    if ((scrq_model.eq.5).OR.(scrq_model.eq.6).OR.&
 &      (scrq_model.eq.7).OR.(scrq_model.eq.8)) then
      coul_scr = 1.0 - 1.0/par_IMPSOLV(2)
    else
      coul_scr = 1.0 - 1.0/sqrt(par_IMPSOLV(2))
    end if
  end if
end
!
!-----------------------------------------------------------------------
!
subroutine get_dipgrps()
!
  use iounit
  use polypep
  use sequen
  use inter
  use aminos
  use atoms
  use energies
  use params
  use molecule
  use system
  use mpistuff
  use clusters
!
  implicit none
!
  integer i,j,k,ii,jj,rs,left,g,gg,hh,atmi,atpl,n1,n2,n3,n4
  integer ptv(maxval(at(1:nseq)%na)),dppdb,nccls,ccl,lefto
  integer, ALLOCATABLE:: donit(:)
  integer, ALLOCATABLE:: good_types1(:),good_types2(:),good_types3(:),good_types4(:)
  integer search_o,search_d,iters,pts,freeunit,swap_success,swap_dips
  logical proceed,exists
  character(100) fn
  character(3) resname
#ifdef ENABLE_MPI
  integer tl
  character(3) nod
#endif
  RTYPE tc,ceps,tar_c(maxval(at(1:nseq)%na))
  type(t_cnblst), ALLOCATABLE:: rsncp(:)
!
  ceps = 1.0*(10.0**(-precision(ceps)))
!
  allocate(rsncp(nseq))
  do rs=1,nseq
    allocate(rsncp(rs)%dis(at(rs)%na+1)) ! data structure is hijacked
    rsncp(rs)%dis(:) = 0.0
    rsncp(rs)%nbs = 0 
    rsncp(rs)%alsz = -1
  end do
!
  call read_ncpatchfile(rsncp)
!
  do rs=1,nseq
    resname = amino(seqtyp(rs))
    ! Martin : Added hard coded exceptions, this is in order to get around having differently defined dipolar groups for each 
    ! protonation forms, which would make the interpolation impossible
    
    at(rs)%ndpgrps = 0
        
    if (at(rs)%npol.eq.0) then
      cycle
    end if
    tc = sum(atq(at(rs)%pol(1:at(rs)%npol)))
    
    if ((dpgrp_corr.eqv..false.).and.((do_pka.eqv..true.).or.(do_hsq.eqv..true.))) then 
        print *,"You are going to need to use the Dipolar group correction (FMCSC_DPGRP_CORR) &
        &if you are to use q-cannonical features (that is if FMCSC_DO_PKA or FMCSC_DO_HSQ is equal to 1)"
        call fexit()
        
    else if (dpgrp_corr.eqv..true.) then        
        if ((resname.eq."HIP").or.(resname.eq."HIE").or.(resname.eq."HID").or.(resname.eq."HEX").or.&
            &(resname.eq."HDX").or.(resname.eq."TYX").or.(resname.eq."TYO").or.(resname.eq."TYR").or.&
            &(resname.eq."SXP").or.(resname.eq."S1P").or.(resname.eq."S2P").or.(resname.eq."HIX")) then 
            ! Doing a quick implementation
            if ((resname.eq."HIP").or.(resname.eq."HIE").or.(resname.eq."HID").or.(resname.eq."HEX").or.&
                &(resname.eq."HDX").or.(resname.eq."HIX")) then 
                at(rs)%ndpgrps=3
            else if ((resname.eq."TYX").or.(resname.eq."TYO").or.(resname.eq."TYR")) then 
                at(rs)%ndpgrps=4
            else if  ((resname.eq."SXP").or.(resname.eq."S1P").or.(resname.eq."S2P")) then 
                at(rs)%ndpgrps=3
            end if 
            
            allocate(at(rs)%dpgrp(at(rs)%ndpgrps))
            !first, assign the backbone dipolar group. O and CO is its own, NH-HN-CA-HA another one
            
            at(rs)%dpgrp(1)%nats=2
            at(rs)%dpgrp(2)%nats=4
            allocate(at(rs)%dpgrp(1)%ats(2))
            allocate(at(rs)%dpgrp(2)%ats(4))

            if ((resname.eq."HIE").or.(resname.eq."HID")) then 
                at(rs)%dpgrp(3)%nats=11
                allocate(at(rs)%dpgrp(3)%ats(at(rs)%dpgrp(3)%nats))

            else if ((resname.eq."HIP").or.(resname.eq."HEX").or.(resname.eq."HDX").or.(resname.eq."HIX")) then 
                at(rs)%dpgrp(3)%nats=12
                allocate(at(rs)%dpgrp(3)%ats(at(rs)%dpgrp(3)%nats))
            else if (resname.eq."TYO") then 
                at(rs)%dpgrp(3)%nats=3
                allocate(at(rs)%dpgrp(3)%ats(at(rs)%dpgrp(3)%nats))
                at(rs)%dpgrp(4)%nats=10
                allocate(at(rs)%dpgrp(4)%ats(at(rs)%dpgrp(4)%nats))
            else if ((resname.eq."TYR").or.(resname.eq."TYX")) then 
                at(rs)%dpgrp(3)%nats=3
                allocate(at(rs)%dpgrp(3)%ats(at(rs)%dpgrp(3)%nats))
                at(rs)%dpgrp(4)%nats=11
                allocate(at(rs)%dpgrp(4)%ats(at(rs)%dpgrp(4)%nats))
            else if ((resname.eq."SXP").or.(resname.eq."S1P")) then 
                at(rs)%dpgrp(3)%nats=9
                allocate(at(rs)%dpgrp(3)%ats(at(rs)%dpgrp(3)%nats))
            else if ((resname.eq."S2P")) then 
                at(rs)%dpgrp(3)%nats=8
                allocate(at(rs)%dpgrp(3)%ats(at(rs)%dpgrp(3)%nats))
            end if 
            
            allocate(good_types1(14))
            allocate(good_types2(28))
            allocate(good_types3(12+10+10+11+11+1+4+19))
            allocate(good_types4(8+7+5))

            ! Martin good_types has all the possible biotypes that can be use in part of the last dipolar group of his (to make
            ! coding simpler)            
            good_types1=(/159,161,176,178,192,194,1127,1129,122,120,1161,1163,1270,1272/)! CoHs
            good_types2=(/157,158,160,162,174,175,177,179,190,191,193,195,118,119,121,123,1125,1126,1128,1130,&
                         &1159,1160,1162,1164,1268,1269,1271,1273/)!SXP! Rest of backbone
            good_types3=(/163,164,165,165,166,167,168,169,170,171,172,173,& ! HIP
                         &180,181,182,183,184,185,186,187,188,189,& !HID
                         &196,197,198,199,200,201,202,203,204,205,& !HIE
                         &1229,1230,1231,1232,1233,1234,1235,1236,1237,1238,1239,& !HEX
                         &1240,1241,1242,1243,1244,1245,1246,1247,1248,1249,1250,& !HDX
                         &1208,& !Proton off
                         &1131,1132,124,125,&!,&!TYO's and TYR CH2
                         &1165,1166,1167,1168,1169,1170,1171,1274,1275,1276,1277,1278,1292,1293,1294,1295,1296,1297,1298/)! all of the side chain of SXP 
                         !(Yes is is a lot, but the partial charges in S1P won't allow any &
                        !different)/)
            good_types4=(/126,127,128,129,130,131,132,133,& !TYR
                         &1133,1134,1135,1136,1137,1138,1139,& !TYO
                         &1312,1313,1314,1315,1316,&!TYX
                         &1208/)! Proton off

            n1=0
            n2=0
            n3=0
            n4=0
            do i=1,at(rs)%npol 
                do j=1,size(good_types1)
                    if (good_types1(j).eq.b_type(at(rs)%pol(i))) then 
                        n1=n1+1
                        at(rs)%dpgrp(1)%ats(n1)=at(rs)%pol(i)
                        exit
                    end if 
                end do 
                do j=1,size(good_types2)
                    if (good_types2(j).eq.b_type(at(rs)%pol(i))) then 
                        n2=n2+1
                        at(rs)%dpgrp(2)%ats(n2)=at(rs)%pol(i)
                        exit 
                    end if 
                end do 
                do j=1,size(good_types3)
                    if (good_types3(j).eq.b_type(at(rs)%pol(i))) then 
                        n3=n3+1
                        at(rs)%dpgrp(3)%ats(n3)=at(rs)%pol(i)
                        exit
                    end if 
                end do 
                do j=1,size(good_types4)
                    if (good_types4(j).eq.b_type(at(rs)%pol(i))) then 
                        n4=n4+1
                        at(rs)%dpgrp(4)%ats(n4)=at(rs)%pol(i)
                        exit
                    end if 
                end do 
            end do     

            deallocate(good_types1)
            deallocate(good_types2)
            deallocate(good_types3)
            deallocate(good_types4)
            
            cycle 
        end if 
    end if 

    

    if (netchg(rs).ne.0) then
!     do nothing, this is currently a redundant test
    else if (abs(tc).gt.at(rs)%npol*ceps) then
      write(ilog,*) 'Fatal. Dipolar groups can only be set up with n&
 &et neutral residues (i.e., there can be no cross-residue polarizat&
 &ion (',tc,')'
      call fexit()
    end if

    atmi = 0
    atpl = 0
    do k=1,at(rs)%npol
      if (atq(at(rs)%pol(k)).lt.0.0) atmi = atmi + 1
      if (atq(at(rs)%pol(k)).gt.0.0) atpl = atpl + 1
    end do
!
!   allocate and initialize the temporary array 
    allocate(donit(at(rs)%npol))
    do i=1,at(rs)%npol
      donit(i) = 0
    end do
!
!   for dipole groups there is only one possibility if there is only a single
!   positive or negative charge (regardless of complexity otherwise)
!   for charge groups, there is only one possibility if there are no opp. charges
!   (- for + and + for -) left 
    if ((netchg(rs).eq.0).AND.((atmi.eq.1).OR.(atpl.eq.1)).AND.(rsncp(rs)%nbs.le.1)) then
      at(rs)%ndpgrps = at(rs)%ndpgrps + 1
      do k=1,at(rs)%npol
        donit(k) = at(rs)%ndpgrps
      end do
    else if ((netchg(rs).gt.0).AND.(atmi.eq.0).AND.(rsncp(rs)%nbs.le.1)) then
      at(rs)%ndpgrps = at(rs)%ndpgrps + 1
      do k=1,at(rs)%npol
        donit(k) = at(rs)%ndpgrps
      end do
    else if ((netchg(rs).lt.0).AND.(atpl.eq.0).AND.(rsncp(rs)%nbs.le.1)) then
      at(rs)%ndpgrps = at(rs)%ndpgrps + 1
      do k=1,at(rs)%npol
        donit(k) = at(rs)%ndpgrps
      end do
    else
      if (rsncp(rs)%nbs.gt.0) then
        if (rsncp(rs)%dis(rsncp(rs)%nbs).ne.0.0) then
          rsncp(rs)%nbs = rsncp(rs)%nbs + 1
          rsncp(rs)%dis(rsncp(rs)%nbs) = 0.0 ! may or may not come into play
        end if
      end if
      proceed = .false.
      left = at(rs)%npol
      search_o = 2
      search_d = 2
      iters = 0
!     the algorithm below really gets into trouble if there is multiple, pseudo-independent
!     charged groups in the residue -> define some exception rules if no patch provided
      if (rsncp(rs)%nbs.eq.0) then
        if (abs(netchg(rs)-tc).gt.at(rs)%npol*ceps) then
          write(ilog,*) 'Warning. Due to residue carrying a fractional net charge, charge group partitioning &
 &may suffer for residue ',rs,' (',resname,').'
          if (nint(tc).eq.0) then
            nccls = 2
            ccl = 1
            tar_c(1) = tc
            tar_c(2) = 0.0
          else
            nccls = 3
            ccl = 1
            tar_c(1) = tc - 1.0*nint(tc)
            tar_c(2) = nint(tc)*1.0
            tar_c(3) = 0.0
          end if
        else if (((resname.eq.'ASP').OR.(resname.eq.'CYX').OR.&
   &              (resname.eq.'S1P').OR.(resname.eq.'TYO').OR.&
   &              (resname.eq.'T1P').OR.(resname.eq.'Y1P'))&
   &               .AND.(moltermid(molofrs(rs),1).eq.1)&
   &               .AND.(rs.eq.rsmol(molofrs(rs),1))&
   &               .AND.(rs.ne.rsmol(molofrs(rs),2))&
   &               .AND.(netchg(rs).eq.0)) then
          nccls = 3
          ccl = 1
          tar_c(1) = -1.0
          tar_c(2) = 1.0
          tar_c(3) = 0.0
        else if (((resname.eq.'ARG').OR.(resname.eq.'LYS')&
   &          .OR.(resname.eq.'DAB').OR.(resname.eq.'ORN')&
   &          .OR.(resname.eq.'HIP'))&
   &               .AND.(moltermid(molofrs(rs),2).eq.1)&
   &               .AND.(rs.eq.rsmol(molofrs(rs),2))&
   &               .AND.(rs.ne.rsmol(molofrs(rs),1))&
   &               .AND.(netchg(rs).eq.0)) then
          nccls = 3
          ccl = 1
          tar_c(1) = -1.0
          tar_c(2) = 1.0
          tar_c(3) = 0.0
        else if (((resname.eq.'ASP').OR.(resname.eq.'GLU').OR.& 
   &              (resname.eq.'S1P').OR.(resname.eq.'TYO').OR.&
   &              (resname.eq.'T1P').OR.(resname.eq.'Y1P'))&
   &               .AND.(moltermid(molofrs(rs),2).eq.1)&
   &               .AND.(rs.eq.rsmol(molofrs(rs),2))&
   &               .AND.(rs.ne.rsmol(molofrs(rs),1))&
   &               .AND.(netchg(rs).eq.-2)) then
          nccls = 3
          ccl = 1
          tar_c(1) = -1.0
          tar_c(2) = -1.0
          tar_c(3) = 0.0
        else if (((resname.eq.'ARG').OR.(resname.eq.'LYS')&
   &          .OR.(resname.eq.'DAB').OR.(resname.eq.'ORN')&
   &          .OR.(resname.eq.'HIP'))&
   &               .AND.(moltermid(molofrs(rs),1).eq.1)&
   &               .AND.(rs.eq.rsmol(molofrs(rs),1))&
   &               .AND.(rs.ne.rsmol(molofrs(rs),2))&
   &               .AND.(netchg(rs).eq.2)) then
          nccls = 3
          ccl = 1
          tar_c(1) = 1.0
          tar_c(2) = 1.0
          tar_c(3) = 0.0
        else if (netchg(rs).ne.0) then
            
          nccls = 2
          ccl = 1
          tar_c(1) = 1.0*netchg(rs)
          tar_c(2) = 0.0
        else
          nccls = 1
          ccl = 1
          tar_c(1) = 0.0
        end if
      else
        ccl = 1
        nccls = rsncp(rs)%nbs
        tar_c(1:nccls) = rsncp(rs)%dis(1:nccls)
      end if
!
      do while (proceed.EQV..false.)
!       
        iters = iters + 1
        lefto = left
        tc = tar_c(ccl)
!
        if (search_d.eq.2) then
!
          do i=1,at(rs)%npol
            if (donit(i).gt.0) cycle
            ii = at(rs)%pol(i)
            pts = 0
            do j=1,n12(ii)
              jj = i12(j,ii)
              if (atmres(jj).ne.rs) cycle
              do k=1,at(rs)%npol
                if (at(rs)%pol(k).eq.jj) then
                  if (donit(k).eq.0) then
                    pts = pts + 1
                    ptv(pts) = k
                  end if
                end if
              end do
            end do
!
            call search_dips(pts,ptv,search_o,rs,i,donit,left,tc)
            if ((left.lt.lefto).AND.(tar_c(ccl).ne.0.0)) exit
!
          end do
!
          if (search_o.lt.(MAXVLNC+1)) then 
            search_o = search_o + 1
          else if (search_o.eq.(MAXVLNC+1)) then
            search_d = search_d + 1
            search_o = 2
          end if
!
        else if (search_d.eq.3) then
!
          do i=1,at(rs)%npol
            if (donit(i).gt.0) cycle
            ii = at(rs)%pol(i)
            do j=1,n12(ii)
              pts = 0
              jj = i12(j,ii)
              if (atmres(jj).ne.rs) cycle
              do k=1,at(rs)%npol
                if (at(rs)%pol(k).eq.jj) then
                  if (donit(k).eq.0) then
                    pts = pts + 1
                    ptv(pts) = k
                  end if
                end if
              end do
              do g=1,n12(ii)
                gg = i12(g,ii)
                if (atmres(gg).ne.rs) cycle
                if (gg.eq.jj) cycle
                do k=1,at(rs)%npol
                  if (at(rs)%pol(k).eq.gg) then
                    if (donit(k).eq.0) then
                      pts = pts + 1
                      ptv(pts) = k
                    end if
                  end if
                end do
              end do
              do g=1,n12(jj)
                gg = i12(g,jj)
                if (atmres(gg).ne.rs) cycle
                if (gg.eq.ii) cycle
                do k=1,at(rs)%npol
                  if (at(rs)%pol(k).eq.gg) then
                    if (donit(k).eq.0) then
                      pts = pts + 1
                      ptv(pts) = k
                    end if
                  end if
                end do
              end do
!
              call search_dips(pts,ptv,search_o,rs,i,donit,left,tc)
              if ((left.lt.lefto).AND.(tar_c(ccl).ne.0.0)) exit
!
            end do
!
            if ((left.lt.lefto).AND.(tar_c(ccl).ne.0.0)) exit
!
          end do
!
          if (search_o.lt.2*MAXVLNC) then 
            search_o = search_o + 1
          else if (search_o.eq.2*MAXVLNC) then
            search_d = search_d + 1
            search_o = 2
          end if
!
        else if (search_d.eq.4) then
!
          do i=1,at(rs)%npol
            if (donit(i).gt.0) cycle
            ii = at(rs)%pol(i)
            pts = 0
            do j=1,n12(ii)
              jj = i12(j,ii)
              if (atmres(jj).ne.rs) cycle
              do k=1,at(rs)%npol
                if (at(rs)%pol(k).eq.jj) then
                  if (donit(k).eq.0) then
                    pts = pts + 1
                    ptv(pts) = k
                  end if
                end if
              end do
            end do
            do j=1,n13(ii)
              jj = i13(j,ii)
              if (atmres(jj).ne.rs) cycle
              do k=1,at(rs)%npol
                if (at(rs)%pol(k).eq.jj) then
                  if (donit(k).eq.0) then
                    pts = pts + 1
                    ptv(pts) = k
                  end if
                end if
              end do
            end do
!
            call search_dips(pts,ptv,search_o,rs,i,donit,left,tc)
            if ((left.lt.lefto).AND.(tar_c(ccl).ne.0.0)) exit
!
          end do
!
          if (search_o.lt.(MAXVLNC*MAXVLNC+1)) then 
            search_o = search_o + 1
          else if (search_o.eq.(MAXVLNC*MAXVLNC+1)) then
            search_d = search_d + 1
            search_o = 2
          end if
!
        else if (search_d.eq.5) then
!         this depth is still missing (max. 1-6). placeholder like s.d. 4
!
          do i=1,at(rs)%npol
            if (donit(i).gt.0) cycle
            ii = at(rs)%pol(i)
            pts = 0
            do j=1,n12(ii)
              jj = i12(j,ii)
              if (atmres(jj).ne.rs) cycle
              do k=1,at(rs)%npol
                if (at(rs)%pol(k).eq.jj) then
                  if (donit(k).eq.0) then
                    pts = pts + 1
                    ptv(pts) = k
                  end if
                end if
              end do
            end do
            do j=1,n13(ii)
              jj = i13(j,ii)
              if (atmres(jj).ne.rs) cycle
              do k=1,at(rs)%npol
                if (at(rs)%pol(k).eq.jj) then
                  if (donit(k).eq.0) then
                    pts = pts + 1
                    ptv(pts) = k
                  end if
                end if
              end do
            end do
!
            call search_dips(pts,ptv,search_o,rs,i,donit,left,tc)
            if ((left.lt.lefto).AND.(tar_c(ccl).ne.0.0)) exit
!
          end do
!
          if (search_o.lt.(MAXVLNC*MAXVLNC+1)) then 
            search_o = search_o + 1
          else if (search_o.eq.(MAXVLNC*MAXVLNC+1)) then
            search_d = search_d + 1
            search_o = 2
          end if
!
        else if (search_d.eq.6) then
!
          do i=1,at(rs)%npol
            if (donit(i).gt.0) cycle
            ii = at(rs)%pol(i)
            pts = 0
            do j=1,n12(ii)
              jj = i12(j,ii)
              if (atmres(jj).ne.rs) cycle
              do k=1,at(rs)%npol
                if (at(rs)%pol(k).eq.jj) then
                  if (donit(k).eq.0) then
                    pts = pts + 1
                    ptv(pts) = k
                  end if
                end if
              end do
            end do
            do j=1,n13(ii)
              jj = i13(j,ii)
              if (atmres(jj).ne.rs) cycle
              do k=1,at(rs)%npol
                if (at(rs)%pol(k).eq.jj) then
                  if (donit(k).eq.0) then
                    pts = pts + 1
                    ptv(pts) = k
                  end if
                end if
              end do
            end do
            do j=1,n14(ii)
              jj = i14(j,ii)
              if (atmres(jj).ne.rs) cycle
              do k=1,at(rs)%npol
                if (at(rs)%pol(k).eq.jj) then
                  if (donit(k).eq.0) then
                    pts = pts + 1
                    ptv(pts) = k
                  end if
                end if
              end do
            end do
!
            call search_dips(pts,ptv,search_o,rs,i,donit,left,tc)
            if ((left.lt.lefto).AND.(tar_c(ccl).ne.0.0)) exit
!
          end do
!
          if (search_o.lt.(MAXVLNC*MAXVLNC*MAXVLNC+1)) then 
            search_o = search_o + 1
          else if (search_o.eq.(MAXVLNC*MAXVLNC*MAXVLNC+1)) then
            search_d = search_d + 1
            search_o = 2
          end if
!
        end if
!
        if (left.eq.0) then
          proceed = .true.
          cycle
        end if
!
        atmi = 0
        atpl = 0
        do k=1,at(rs)%npol
          if (donit(k).gt.0) cycle ! Martin : donit is litterally done it 
          if (atq(at(rs)%pol(k)).lt.0.0) atmi = atmi + 1 ! Martin : I think this is the neutral thing. 
          if (atq(at(rs)%pol(k)).gt.0.0) atpl = atpl + 1 !
  
        end do
        ! Martin : so in the following there is a case for net neutral and charge, but each of them assume an equal number of 
        ! oppositly hcraged residues
        ! First case is if the net charge is zero and the number of + and - is equivalent 
        if ((netchg(rs).eq.0).AND.((atmi.eq.1).OR.(atpl.eq.1)).AND.(ccl.eq.nccls)) then
          at(rs)%ndpgrps = at(rs)%ndpgrps + 1
          do k=1,at(rs)%npol
            if (donit(k).eq.0) then
              donit(k) = at(rs)%ndpgrps
            end if
          end do
          left = 0
        else if ((netchg(rs).gt.0).AND.(atmi.eq.0).AND.(ccl.eq.nccls)) then
          at(rs)%ndpgrps = at(rs)%ndpgrps + 1
          do k=1,at(rs)%npol
            if (donit(k).eq.0) then
              donit(k) = at(rs)%ndpgrps
            end if
          end do
          left = 0
        else if ((netchg(rs).lt.0).AND.(atpl.eq.0).AND.(ccl.eq.nccls)) then
          at(rs)%ndpgrps = at(rs)%ndpgrps + 1
          do k=1,at(rs)%npol
            if (donit(k).eq.0) then
              donit(k) = at(rs)%ndpgrps
            end if
          end do
          left = 0
        end if
        if (left.eq.0) then
          proceed = .true.
        else if (left.eq.1) then
          write(ilog,*) 'Automatic dipoles: inconsistent charges.'
        end if
!
 46     format(1x,a,i7,a,a3,a)
 48     format(1x,a,g10.4,a,i7,a,a3,a)
        if ((search_d.gt.6).AND.(ccl.eq.nccls)) then
          write(ilog,46) 'Automatic dipole group detection chokes on&
 & residue ',rs,' (',resname,').'
          write(ilog,*) 'Forming a single group out of the remaining&
 & ',left,' polar atoms.'
          at(rs)%ndpgrps = at(rs)%ndpgrps + 1
          pts = 0
          do k=1,at(rs)%npol
            if (donit(k).eq.0) then
              donit(k) = at(rs)%ndpgrps
              pts = pts + 1
              ptv(pts) = k
            end if
          end do
          swap_success = swap_dips(rs,at(rs)%ndpgrps,pts,ptv,at(rs)%npol,donit)
          if (swap_success.eq.1) then
            write(ilog,*) 'Sane dipole groups could be recovered by swaps of already&
 & assigned atoms.'
          end if
          left = 0
          proceed = .true.
!       for special targets, we have to manually advance the charge-counter in order to prevent pulling out further such units
        else if (search_d.gt.6) then
          write(ilog,48) 'Automatic dipole group detection chokes for group with target &
 & of ',tar_c(ccl),' in residue ',rs,' (',resname,'). Skipped (likely undesirable).'
          ccl = ccl + 1
          search_d = 2
          search_o = 2
        else if ((left.lt.lefto).AND.(ccl.lt.nccls)) then
          ccl = ccl + 1
          search_d = 2
          search_o = 2
        end if
        if (iters.gt.50000) then
          write(ilog,*) 'Automatic dipole group detection failed.'
          call fexit()
        end if
!
!       if we found a group with a net charge, we have to manually advance
!       the charge-counter in order to prevent pulling out non-sensical
!       net-charged units
!        if ((left.lt.lefto).AND.(ccl.lt.nccls).AND.&
! &          (abs(tar_c(ccl)).lt.at(rs)%npol*ceps)) then
!          ccl = ccl + 1
!        end if
      end do
!
    end if
!
    if (at(rs)%ndpgrps.eq.0) cycle
    allocate(at(rs)%dpgrp(at(rs)%ndpgrps))
    
    do i=1,at(rs)%ndpgrps
      hh = 0
      do k=1,at(rs)%npol
        if (donit(k).eq.i) then ! Martin : If the num of polar groups is equal to the current 
          hh = hh + 1
        end if
      end do
     
      allocate(at(rs)%dpgrp(i)%ats(hh))
      at(rs)%dpgrp(i)%nats = hh
      hh = 0
      at(rs)%dpgrp(i)%tc = 0.0
      do k=1,at(rs)%npol
        if (donit(k).eq.i) then
          hh = hh + 1
          at(rs)%dpgrp(i)%tc=at(rs)%dpgrp(i)%tc + atq(at(rs)%pol(k))
          at(rs)%dpgrp(i)%ats(hh) = at(rs)%pol(k)
        end if
      end do
!     init to be correct for the rare exceptions when tc is fractional
      if (at(rs)%dpgrp(i)%tc.lt.0.0) then
        at(rs)%dpgrp(i)%nc = floor(at(rs)%dpgrp(i)%tc)
      else if (at(rs)%dpgrp(i)%tc.gt.0.0) then
        at(rs)%dpgrp(i)%nc = ceiling(at(rs)%dpgrp(i)%tc)
      end if
      
!     override
      do k=-20,20
        if (abs(at(rs)%dpgrp(i)%tc-1.0*k).le.at(rs)%npol*ceps) then
          at(rs)%dpgrp(i)%nc = k
          exit
        end if
      end do
    end do
!
    deallocate(donit)
!
  end do
! Martin : This should be where the changes have to be made 
  do rs=1,nseq
!    write(*,*) 'here for ',rs,' and ',at(rs)%dpgrp(1)%nats
    do i=1,at(rs)%ndpgrps
!      at(rs)%dpgrp(i)%qnm = sum(abs(atq(at(rs)%dpgrp(i)%ats(1:at(rs)%dpgrp(i)%nats))))
!      which_dpg(at(rs)%dpgrp(i)%ats(1:at(rs)%dpgrp(i)%nats)) = i
      at(rs)%dpgrp(i)%qnm = 0.0
      do k=1,at(rs)%dpgrp(i)%nats
        at(rs)%dpgrp(i)%qnm = at(rs)%dpgrp(i)%qnm + abs(atq(at(rs)%dpgrp(i)%ats(k)))
        which_dpg(at(rs)%dpgrp(i)%ats(k)) = i
      end do
!     write(*,*) 'here for ',rs,' and ',i,' and ',at(rs)%dpgrp(i)%qnm
    end do
  end do
!
!
! build the global list of monopole groups and populate chgflag
  call assemble_cgplst(rsncp)
!
!
  do rs=1,nseq
    deallocate(rsncp(rs)%dis)
  end do
  deallocate(rsncp)
!
  if (dip_report.EQV..true.) then
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      if (myrank.ne.0) return
      tl = 3
      call int2str(myrank,nod,tl)
      fn = 'N_'//nod(1:tl)//'_DIPOLE_GROUPS.vmd'
    else if (use_MPIAVG.EQV..true.) then
      if (myrank.ne.0) return
      fn = 'DIPOLE_GROUPS.vmd'
    end if
#else
    fn = 'DIPOLE_GROUPS.vmd'
#endif
    call strlims(fn,ii,jj)
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      dppdb = freeunit()
      open(unit=dppdb,file=fn(ii:jj),status='old')
      close(unit=dppdb,status='delete')
    end if
    dppdb=freeunit()
    open(unit=dppdb,file=fn(ii:jj),status='new')
!
    write(dppdb,*) 'color Display Background black'
    write(dppdb,*) 'material change opacity Transparent 0.52'
    write(dppdb,*) 'mol load pdb ',basename(1:bleng),'_START.pdb'
    write(dppdb,*) 'mol modcolor Type'
    write(dppdb,*) 'mol modstyle 0 0 Licorice'
    gg = 1
    do rs=1,nseq
      if (at(rs)%ndpgrps.gt.gg) gg = at(rs)%ndpgrps
    end do
    do g=1,gg
      write(dppdb,*) ' mol selection "index \'
      do rs=1,nseq
        if (at(rs)%ndpgrps.ge.g) then 
          do k=1,at(rs)%dpgrp(g)%nats
            write(dppdb,*) at(rs)%dpgrp(g)%ats(k)-1,'\'
          end do
        end if
      end do
      write(dppdb,*) '"'
      write(dppdb,*) ' mol color colorID ',g
      write(dppdb,*) ' mol rep VdW'
      write(dppdb,*) ' mol material Transparent'
      write(dppdb,*) ' mol addrep 0'
    end do
    close(unit=dppdb)
  end if
!
end
!
!-----------------------------------------------------------------------
!
! there are prob. more elegant and clear ways to get all unique combinations for a set of
! varying size (# = binomial coefficient), but this version works and is reasonably compact
!
subroutine search_dips(pts,ptv,search_o,rs,i,donit,left,tar_c)
!
  use iounit
  use polypep
  use sequen
  use inter
  use aminos
  use atoms
  use energies
  use params
  use molecule
!
  implicit none
!
  integer rs,search_o,left,ii,i,donit(at(rs)%npol),ptv(maxval(at(1:nseq)%na))
  integer sup(search_o)
  integer pts,pt1,lk
  integer lk2
  logical keepit,chkdipset,itered
  RTYPE ceps,tc,tar_c
!
  ceps = 1.0*(10.0**(-precision(ceps)))
!
  ii = at(rs)%pol(i)
!  write(*,*) 'in with ',search_o,rs,ii,pts
!
  lk = 0
 66   format(100(i4))
  if ((pts+1).lt.search_o) return
  do lk=1,search_o-1
    sup(lk) = lk
  end do
!
  do while (sup(1).le.pts-(search_o-1)+1)
!
!   process our input vector for total charge
    tc = atq(ii)
    do pt1=1,search_o-1
      tc = tc + atq(at(rs)%pol(ptv(sup(pt1))))
    end do
!   if it matches the current target charge, check for inconsistency
    if (abs(tc-tar_c).le.at(rs)%npol*ceps) then
      keepit = chkdipset(search_o,pts,ptv,rs,ii,sup)
!     build in an extra check for full nucleotides (which cause problems frequently)
      if ((keepit.EQV..true.).AND.((tar_c+1.0).le.at(rs)%npol*ceps).AND.&
 &  (seqpolty(rs).eq.'N').AND.(seqflag(rs).ne.24)) then
        keepit = .false.
        do pt1=1,search_o-1
          if (at(rs)%pol(ptv(sup(pt1))).eq.nuci(rs,2)) then
            keepit = .true.
            exit
          end if
        end do
      end if
      if (keepit.EQV..true.) then
!       if accepted, add charge group and update arrays
        left = left - search_o
        at(rs)%ndpgrps = at(rs)%ndpgrps + 1
        donit(i) = at(rs)%ndpgrps
        do pt1=1,search_o-1
          donit(ptv(sup(pt1))) = at(rs)%ndpgrps
        end do
!        write(*,*) 'found gr with ',tar_c
!        write(*,66) ii,(at(rs)%pol(ptv(sup(pt1))),
! &                                     pt1=1,search_o-1)
        exit
      end if
    end if
!
!   if nothing was found, advance to the next possible member of the set
    itered = .false.
    do lk=search_o-1,1,-1
      if (sup(lk).lt.(pts-(search_o-1-lk))) then
        itered = .true.
        sup(lk) = sup(lk) + 1
        do lk2=lk+1,search_o-1
          sup(lk2) = lk2 - lk + sup(lk)
        end do
        exit
      end if
    end do
!
    if (itered.EQV..false.) exit
!
  end do
!
end
!
!-----------------------------------------------------------------------
!
! this function rejects a proposed charge group based on two criteria:
! 1) it uses only one of several atoms with the same biotype (chemically identical)
! 2) it "encapsulates" other polar atoms not part of the same group: while technically
!    ok, it makes little sense to have a charge group immersed in a larger charge group
! note that the error of the second type can still occur if the remaining atoms are grouped
! together as a single group (exception bail-out case)
!
function chkdipset(so,pts,ptv,rs,first,sup)
!
  use polypep
  use atoms
  use sequen
!
  implicit none
!
  integer so,pts,ptv(maxval(at(1:nseq)%na)),rs,pl,first,sup(so)
  integer k4
  integer i,btyr,pt(so),j,k,badone,i2,i3,i4,k2,ii,kk,k3
  logical chkdipset,it_13,it_14,it_15,it_12
!
  chkdipset = .true.
!
  pt(1) = first
  do i=1,so-1
    pt(i+1) = at(rs)%pol(ptv(sup(i)))
  end do
!
  do i=1,so
!   the reference biotype for the "chosen" one
    btyr = b_type(pt(i))
    do j=1,pts
!     now compare to all other biotypes in the set ...
      if (b_type(at(rs)%pol(ptv(j))).eq.btyr) then
!       and kill the group if it doesn't contain that second (third, ...) atom
!       of biotype btyr
        chkdipset = .false.
        do k=1,so
          if (at(rs)%pol(ptv(j)).eq.pt(k)) then
            chkdipset = .true.
            exit
          end if
        end do
        if (chkdipset.EQV..false.) return
      end if
    end do
  end do
!
! this part is a bit complicated:
! we search the whole set for pairs of atoms in the set which bracket a polar
! atom which does not belong to the set (e.g., X-Y-Z with Y polar, and X,Z part
! of the set)
! such an arrangement will always indicate a chemically non-sensical dipole group
! (even though formally it might work out fine).
  do i=1,so
    ii = pt(i)
    do k=i+1,so
      it_12 = .false.
      it_13 = .false.
      it_14 = .false.
      it_15 = .false.
      kk = pt(k)
      do i2=1,n12(ii)
        if (i12(i2,ii).eq.kk) then
          it_12 = .true.
          exit
        end if
      end do
      if (it_12.EQV..true.) cycle
      do i3=1,n13(ii)
        if (i13(i3,ii).eq.kk) then
          it_13 = .true.
          exit
        end if
      end do
      do i4=1,n14(ii)
        if (i14(i4,ii).eq.kk) then
          it_14 = .true.
          exit
        end if
      end do
      if ((it_13.EQV..false.).AND.(it_14.EQV..false.)) then
        do i3=1,n13(ii)
          if (kk.eq.i13(i3,ii)) exit
          do k3=1,n13(kk)
            if (ii.eq.i13(k3,kk)) exit
            if (i13(i3,ii).eq.i13(k3,kk)) then
              it_15 = .true.
              exit
            end if
          end do
        end do
      end if
!     we systematically scan all pairs in the set and now know the connectivity
!     (although 12, and >15 are obviously ignored) for atoms ii,kk
!     the requirement for 13 is that there is a new POLAR atom NOT part of the
!     set, which is directly bound to BOTH ii and kk. if such a thing is found
!     reject and exit immediately
      if (it_13.EQV..true.) then
        do i2=1,n12(ii)
          do pl=1,at(rs)%npol
            badone = -1
            if (i12(i2,ii).eq.at(rs)%pol(pl)) then
              badone = i12(i2,ii) !a polar atom directly bound to ii
              do j=1,so
                if (at(rs)%pol(pl).eq.pt(j)) then
                  badone = -1 !safe: it's part of the set
                  exit
                end if
              end do
            end if
            if (badone.gt.0) then
              do k2=1,n12(kk) !if not safe, scan atoms directly bound to kk
                if (i12(k2,kk).eq.badone) then !if found reject and exit
                  chkdipset = .false.
                  return
                end if
              end do
            end if
          end do
        end do
!     the requirement for 14 is that there is a new POLAR atom NOT part of the
!     set, which is directly bound to ii AND 13 to kk OR which is directly bound
!     to kk and 13 to ii.
      else if (it_14.EQV..true.) then
!       first 12 to ii, 13 to kk
        do i2=1,n12(ii)
          do pl=1,at(rs)%npol
            badone = -1
            if (i12(i2,ii).eq.at(rs)%pol(pl)) then
              badone = i12(i2,ii)
              do j=1,so
                if (at(rs)%pol(pl).eq.pt(j)) then
                  badone = -1
                  exit
                end if
              end do
            end if
            if (badone.gt.0) then
              do k3=1,n13(kk)
                if (i13(k3,kk).eq.badone) then
                  chkdipset = .false.
                  return
                end if
              end do
            end if
          end do
        end do
!       the other way around
        do i3=1,n13(ii)
          do pl=1,at(rs)%npol
            badone = -1
            if (i13(i3,ii).eq.at(rs)%pol(pl)) then
              badone = i13(i3,ii)
              do j=1,so
                if (at(rs)%pol(pl).eq.pt(j)) then
                  badone = -1
                  exit
                end if
              end do
            end if
            if (badone.gt.0) then
              do k2=1,n12(kk)
                if (i12(k2,kk).eq.badone) then
                  chkdipset = .false.
                  return
                end if
              end do
            end if
          end do
        end do
      else if (it_15.EQV..true.) then
!       first 12 to ii, 14 to kk
        do i2=1,n12(ii)
          do pl=1,at(rs)%npol
            badone = -1
            if (i12(i2,ii).eq.at(rs)%pol(pl)) then
              badone = i12(i2,ii)
              do j=1,so
                if (at(rs)%pol(pl).eq.pt(j)) then
                  badone = -1
                  exit
                end if
              end do
            end if
            if (badone.gt.0) then
              do k4=1,n14(kk)
                if (i14(k4,kk).eq.badone) then
                  chkdipset = .false.
                  return
                end if
              end do
            end if
          end do
        end do
!       now 13 to both
        do i3=1,n13(ii)
          do pl=1,at(rs)%npol
            badone = -1
            if (i13(i3,ii).eq.at(rs)%pol(pl)) then
              badone = i13(i3,ii)
              do j=1,so
                if (at(rs)%pol(pl).eq.pt(j)) then
                  badone = -1
                  exit
                end if
              end do
            end if
            if (badone.gt.0) then
              do k3=1,n13(kk)
                if (i13(k3,kk).eq.badone) then
                  chkdipset = .false.
                  return
                end if
              end do
            end if
          end do
        end do
!       finally 14 to ii, 12 to kk
        do i4=1,n14(ii)
          do pl=1,at(rs)%npol
            badone = -1
            if (i14(i4,ii).eq.at(rs)%pol(pl)) then
              badone = i14(i4,ii)
              do j=1,so
                if (at(rs)%pol(pl).eq.pt(j)) then
                  badone = -1
                  exit
                end if
              end do
            end if
            if (badone.gt.0) then
              do k2=1,n12(kk)
                if (i12(k2,kk).eq.badone) then
                  chkdipset = .false.
                  return
                end if
              end do
            end if
          end do
        end do
      end if
    end do
  end do
!
  return
end
!
!-----------------------------------------------------------------------
!
! a subroutine that - given a group assignment - tries to remedy inconsistencies
! by swapping select atoms
! returns 1 upon success, 0 else
!
function swap_dips(rs,dp1,pts,ptv,sob,donit)
!
  use polypep
  use iounit
  use atoms
  use sequen
!
  implicit none
!
  integer i,rs,dp1,so,sob,dp2,k,l,m,mm,buf2,first,pts,pts2,first2,so2,buf
  integer ptv(maxval(at(1:nseq)%na)),sup(sob),ptv2(maxval(at(1:nseq)%na)),donit(at(rs)%npol)
  integer level,maxlevel,swap_dips
  logical chkdipset
  RTYPE ceps
!
  swap_dips = 0
  ceps = 1.0*(10.0**(-precision(ceps)))
!
! nothing fixable here
  if (at(rs)%ndpgrps.le.1) return
!
  do i=1,sob-1
    sup(i) = i+1
  end do
  first = at(rs)%pol(ptv(1))
  so = pts
!
! shouldn't happen
  if (chkdipset(so,pts,ptv,rs,first,sup(1:so)).EQV..true.) then
    swap_dips = 1
    return
  end if
!
  maxlevel = 3
  do level=1,maxlevel
    do dp2=1,at(rs)%ndpgrps
      if (dp2.eq.dp1) cycle
      pts2 = 0
      first2 = 0
      do k=1,at(rs)%npol
        if (donit(k).eq.dp2) then
          if (first2.le.0) first2 = at(rs)%pol(k)
          pts2 = pts2 + 1
          ptv2(pts2) = k
        end if
      end do
      so2 = pts2
      if (level.eq.1) then
!       level 1: swap single atoms
        do k=1,pts
          do l=1,pts2
            if (abs(atq(at(rs)%pol(ptv(k)))-atq(at(rs)%pol(ptv2(l)))).le.at(rs)%npol*ceps) then
              buf = ptv(k)
              ptv(k) = ptv2(l)
              ptv2(l) = buf
              if ((chkdipset(so,pts,ptv,rs,first,sup(1:so)).EQV..true.).AND.&
     &            (chkdipset(so2,pts2,ptv2,rs,first2,sup(1:so2)).EQV..true.)) then
                buf = donit(ptv(k))
                donit(ptv(k)) = donit(ptv2(l))
                donit(ptv2(l)) = buf
                swap_dips = 1
                return
              else
                ptv2(l) = ptv(k)
                ptv(k) = buf
              end if
            end if
          end do
        end do
      else if (level.eq.2) then
!       level 2: swap single atom k vs. pairs in l
        do k=1,pts
          do l=1,pts2
            do m=l+1,pts2
              if (abs(atq(at(rs)%pol(ptv(k)))-atq(at(rs)%pol(ptv2(l)))-atq(at(rs)%pol(ptv2(m))))&
 &                     .le.at(rs)%npol*ceps) then
                buf = ptv(k)
                ptv(k) = ptv2(l)
                ptv2(l) = buf
                ptv(pts+1) = ptv2(m)
                buf2 = ptv2(m)
                do mm=m,pts2-1
                  ptv2(mm) = ptv2(mm+1)
                end do
                pts2 = pts2 - 1
                pts = pts + 1
                so = so + 1
                so2 = so2 - 1
                if ((chkdipset(so,pts,ptv,rs,first,sup(1:so)).EQV..true.).AND.&
 &                (chkdipset(so2,pts2,ptv2,rs,first2,sup(1:so2)).EQV..true.)) then
                  buf = donit(ptv(k))
                  donit(ptv(k)) = donit(ptv2(l))
                  donit(ptv2(l)) = buf
                  donit(ptv2(m)) = buf
                  swap_dips = 1
                  return
                else
                  ptv2(l) = ptv(k)
                  ptv(k) = buf
                  do mm=pts2,m,-1
                    ptv2(mm+1) = ptv2(mm)
                  end do
                  ptv2(m) = buf2
                  pts2 = pts2 + 1
                  so = so - 1
                  pts = pts - 1
                  so2 = so2 + 1
                end if
              end if
            end do
          end do
        end do
      else if (level.eq.3) then
!       level 2: swap single atom l vs. pairs in k
        do k=1,pts2
          do l=1,pts
            do m=l+1,pts
              if (abs(atq(at(rs)%pol(ptv2(k)))-atq(at(rs)%pol(ptv(l)))-atq(at(rs)%pol(ptv(m))))&
 &                     .le.at(rs)%npol*ceps) then
                buf = ptv2(k)
                ptv2(k) = ptv(l)
                ptv(l) = buf
                ptv2(pts2+1) = ptv(m)
                buf2 = ptv(m)
                do mm=m,pts-1
                  ptv(mm) = ptv(mm+1)
                end do
                pts = pts - 1
                pts2 = pts2 + 1
                so2 = so2 + 1
                so = so - 1
                if ((chkdipset(so,pts,ptv,rs,first,sup(1:so)).EQV..true.).AND.&
 &                (chkdipset(so2,pts2,ptv2,rs,first2,sup(1:so2)).EQV..true.)) then
                  buf = donit(ptv2(k))
                  donit(ptv2(k)) = donit(ptv(l))
                  donit(ptv(l)) = buf
                  donit(ptv(m)) = buf
                  swap_dips = 1
                  return
                else
                  ptv(l) = ptv2(k)
                  ptv2(k) = buf
                  do mm=pts,m,-1
                    ptv(mm+1) = ptv(mm)
                  end do
                  ptv(m) = buf2
                  pts = pts + 1
                  so2 = so2 - 1
                  pts2 = pts2 - 1
                  so = so + 1
                end if
              end if
            end do
          end do
        end do
      end if
    end do
  end do
!
  return
!
end
!
!------------------------------------------------------------------------------
!
! this routine builds that list of "monopole" groups and sets the residue level
! chgflag
!
subroutine assemble_cgplst(rsncp)
!
  use iounit
  use sequen
  use polypep
  use cutoffs
  use atoms
  use aminos
  use params
  use system
  use energies
  use mpistuff
  use clusters
!
  implicit none
!
  integer rs,i,j,kk,k,refk,dppdb,ncmi,ncma,freeunit,g,coli,ii,jj,nalcsz
  character(100) fn 
  RTYPE cvec(3),mind,dum,epsi,epsi2
  logical exists,heady,foundit
  RTYPE, ALLOCATABLE:: ncmlst(:)
#ifdef ENABLE_MPI
  character(3) nod
  integer tl
#endif
  type(t_cnblst) rsncp(nseq)
!
  heady = .false.
!
  epsi = 1.0e-3
!
! get the required allocation size
  nalcsz = 0
  do rs=1,nseq
    do i=1,at(rs)%ndpgrps
      if (at(rs)%dpgrp(i)%nc.ne.0) nalcsz = nalcsz + 1
    end do
  end do
!
  allocate(cglst%it(nalcsz))
  allocate(cglst%nc(nalcsz))
  allocate(cglst%tc(nalcsz))
  cglst%ncs = 0
!
  do rs=1,nseq
!
    do i=1,at(rs)%ndpgrps
      epsi2 = max(at(rs)%dpgrp(i)%nats*10.0*(10.0**(-precision(epsi2))),dpgrp_neut_tol)
      if (at(rs)%dpgrp(i)%nc.ne.0) then
!       here we set chgflag for all cases including those where the residue is net neutral, but has hidden charge moieties
!       for groups with fractional charge, there is a user-based tolerance
        if (abs(at(rs)%dpgrp(i)%tc).gt.epsi2) then
          chgflag(rs) = .true.
        else
          if (heady.EQV..false.) then
            write(ilog,*) 'Warning. At least one dipole group with a fractional net charge is treated &
 &as neutral due to setting for FMCSC_POLTOL. Results should be inspected with utmost care.'
            heady = .true.
          end if
          at(rs)%dpgrp(i)%nc = 0 ! pretends to be a dipole group
          cycle
        end if
!       find the "central" atom
        cglst%ncs = cglst%ncs + 1
        cglst%nc(cglst%ncs) = at(rs)%dpgrp(i)%nc
        cglst%tc(cglst%ncs) = at(rs)%dpgrp(i)%tc
        if (at(rs)%dpgrp(i)%nats.eq.1) then
          cglst%it(cglst%ncs) = at(rs)%dpgrp(i)%ats(1)
        else
          cvec(:) = 0.0
          do k=1,at(rs)%dpgrp(i)%nats
            kk = at(rs)%dpgrp(i)%ats(k)
            cvec(1) = cvec(1) + atq(kk)*x(kk)
            cvec(2) = cvec(2) + atq(kk)*y(kk)
            cvec(3) = cvec(3) + atq(kk)*z(kk)
          end do
          cvec(:) = cvec(:)/at(rs)%dpgrp(i)%tc
          mind = huge(mind)
          refk = 1
          do k=1,at(rs)%dpgrp(i)%nats
            kk = at(rs)%dpgrp(i)%ats(k)
            dum = sqrt((cvec(1)-x(kk))**2 + (cvec(2)-y(kk))**2 + &
 &  (cvec(3)-z(kk))**2)
            if (dum.lt.mind) then
              refk = k
              mind = dum
            end if
          end do
          cglst%it(cglst%ncs) = at(rs)%dpgrp(i)%ats(refk)
        end if
        at(rs)%dpgrp(i)%cgn = cglst%ncs
      end if
    end do
!
  end do
!
  if (nc_patched.EQV..true.) then
!
 24   format(' Residue #',i9,' was patched to no longer be treated as if carrying a net charge.')
 25   format(' Residue #',i9,' was patched to now be treated as if carrying a net charge.')
!
    heady = .false.
    do i=1,nseq
      if (rsncp(i)%alsz.lt.0) cycle
      if ((rsncp(i)%alsz.eq.0).AND.(chgflag(i).EQV..true.)) then
        if (heady.EQV..false.) then
          write(ilog,*)
          write(ilog,*) '-- Summary of Applied Patches to Residue-Level Charge Group Status --'
          write(ilog,*)
          heady = .true.
        end if
        write(ilog,24) i
        chgflag(i) = .false.
!      DISABLED
!      else if ((rsncp(i)%alsz.eq.1).AND.(chgflag(i).EQV..false.)) then
!        if (heady.EQV..false.) then
!          write(ilog,*)
!          write(ilog,*) '-- Summary of Applied Patches to Residue-Level Charge Group Status --'
!          write(ilog,*)
!          heady = .true.
!        end if
!        write(ilog,25) i
!        chgflag(i) = .true.
      end if
    end do
    write(ilog,*)
  end if
!
  if (cglst%ncs.le.0) return
!
  if (dip_report.EQV..true.) then
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      if (myrank.ne.0) return
      tl = 3
      call int2str(myrank,nod,tl)
      fn = 'N_'//nod(1:tl)//'_MONOPOLES.vmd'
    else if (use_MPIAVG.EQV..true.) then
      fn = 'MONOPOLES.vmd'
      if (myrank.ne.0) return
    end if
#else
    fn = 'MONOPOLES.vmd'
#endif
    call strlims(fn,ii,jj)
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      dppdb = freeunit()
      open(unit=dppdb,file=fn(ii:jj),status='old')
      close(unit=dppdb,status='delete')
    end if
    dppdb=freeunit()
  
    open(unit=dppdb,file=fn(ii:jj),status='new')
!
    write(dppdb,*) 'color Display Background black'
    write(dppdb,*) 'material change opacity Transparent 0.52'
    write(dppdb,*) 'mol load pdb ',basename(1:bleng),'_START.pdb'
    write(dppdb,*) 'mol modcolor Type'
    write(dppdb,*) 'mol modstyle 0 0 Licorice'
    ncmi = 0
    ncma = 0
    allocate(ncmlst(cglst%ncs))
    ncma = 0
    do i=1,cglst%ncs
      foundit = .false.
      do j=1,ncma
        if (abs(ncmlst(j)-cglst%tc(i)).le.epsi) then
          foundit = .true.
          exit
        end if
      end do
      if (foundit.EQV..false.) then
        ncma = ncma + 1
        ncmlst(ncma) = cglst%tc(i)
!      if (cglst%nc(i).lt.ncmi) ncmi = cglst%nc(i)
!      if (cglst%nc(i).gt.ncma) ncma = cglst%nc(i)
      end if
    end do
    coli = 0
    do g=1,ncma!ncmi,ncma
      heady = .false.
      do i=1,cglst%ncs
        if (abs(cglst%tc(i)-ncmlst(g)).le.epsi) then !(cglst%nc(i).eq.g) then
          if (heady.EQV..false.) then
            write(dppdb,*) ' mol selection "index \'
            heady = .true.
          end if
          write(dppdb,*) cglst%it(i)-1,'\'
        end if
      end do
      if (heady.EQV..true.) then
        coli = coli + 1
        write(dppdb,*) '"'
        write(dppdb,*) ' mol color colorID ',coli
        write(dppdb,*) ' mol rep VdW'
        write(dppdb,*) ' mol material Transparent'
        write(dppdb,*) ' mol addrep 0'
      end if
    end do
    close(unit=dppdb)
    deallocate(ncmlst)
  end if
!
end
!
!----------------------------------------------------------------------------
!
! GROMOS demands NB-excludes similar but not identical to nbsr_model 1, hence a separate
! and non-general fxn (see GROMOS_excludes in unbond.f90 as well)

  subroutine GROMOS_excludes_polar()
!
  use iounit
  use inter
  use sequen
  use aminos
  use atoms
  use polypep
  use molecule
  use system
!
  implicit none
! 
  integer i,j,rs,ii,kk,hb1,shf,shf2,shf3
  character(3) resname
  logical its14
!
  shf = 0
  shf2 = 0
  shf3 = 0
  if (ua_model.gt.0) then
    shf = 1
    shf2 = 2
    shf3 = 3
    if (ua_model.eq.2) then
      write(ilog,*) 'Fatal. The GROMOS exclusion rules are largely specific for hydrogens &
 &attached to aromatic rings and make no sense if those are not explicitly represented (choice &
 &of FMCSC_UAMODEL). Adjust setting for FMCSC_INTERMODEL.'
      call fexit()
    end if
  end if
!
  do rs=1,nseq
    do i=1,nrpolintra(rs)
!
      resname = amino(seqtyp(rs))
      ii = iaa(rs)%polin(i,1)
      kk = iaa(rs)%polin(i,2)
      if ((attyp(ii).eq.0).OR.(attyp(kk).eq.0)) cycle
!
      its14 = .false.
!
      do j=1,n14(ii)
        if (i14(j,ii).eq.kk) its14 = .true.
      end do
      if (its14.EQV..true.) then
        if ((resname.eq.'PHE').OR.(resname.eq.'TRP').OR.(resname.eq.'TYO').OR.&
 &          (resname.eq.'HIE').OR.(resname.eq.'HID').OR.(resname.eq.'HIP')) then
          if (ua_model.eq.0) then
            hb1 = 8
            if (resname.eq.'PHE') hb1 = 9
            if (resname.eq.'TYO') hb1 = 10
            if (resname.eq.'TRP') hb1 = 12
  !         both are sidechain atoms but not beta-hydrogens
            if (((ii.ge.at(rs)%sc(2)).AND.(ii.ne.at(rs)%sc(hb1)).AND.(ii.ne.at(rs)%sc(hb1+1))).AND.&
 &          ((kk.ge.at(rs)%sc(2)).AND.(kk.ne.at(rs)%sc(hb1)).AND.(kk.ne.at(rs)%sc(hb1+1)))) then
              fudge(rs)%elin(i) = 0.0
            end if
          else
!           both are sidechain atoms
            if ((ii.ge.at(rs)%sc(1)).AND.(kk.ge.at(rs)%sc(1))) then
              fudge(rs)%elin(i) = 0.0
            end if
          end if
        else if (resname.eq.'Y1P') then
          write(ilog,*) 'Fatal. Polar GROMOS exclusion rules are currently not supported for residue type Y1P. Please &
 &check back later.'
          call fexit()
        else if (resname.eq.'TYR') then
          hb1 = 10
          if (ua_model.eq.0) then
!           both are sidechain atoms but not beta-hydrogens and not phenolic hydrogen
            if (((ii.ge.at(rs)%sc(2)).AND.(ii.ne.at(rs)%sc(hb1)).AND.&
 &             (ii.ne.at(rs)%sc(hb1+1)).AND.(ii.ne.at(rs)%sc(16))).AND.&
 &            ((kk.ge.at(rs)%sc(2)).AND.(kk.ne.at(rs)%sc(hb1)).AND.&
 &             (kk.ne.at(rs)%sc(hb1+1)).AND.(kk.ne.at(rs)%sc(16)))) then
              fudge(rs)%elin(i) = 0.0
            end if
          else
!           both are sidechain atoms and not phenolic hydrogen
            if (((ii.ge.at(rs)%sc(1)).AND.(ii.ne.at(rs)%sc(at(rs)%nsc))).AND.&
 &              ((kk.ge.at(rs)%sc(1)).AND.(kk.ne.at(rs)%sc(at(rs)%nsc)))) then
              fudge(rs)%elin(i) = 0.0
            end if
          end if
        else if ((resname.eq.'RPA').OR.(resname.eq.'RIA')) then
          if ((ii.ge.at(rs)%sc(8-shf2)).AND.(ii.le.at(rs)%sc(19-shf2)).AND.&
 &            (kk.ge.at(rs)%sc(8-shf2)).AND.(kk.le.at(rs)%sc(19-shf2))) then
            fudge(rs)%elin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(8-shf2)).AND.(kk.le.at(rs)%sc(19-shf2))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(8-shf2)).AND.(ii.le.at(rs)%sc(19-shf2)))) then
            fudge(rs)%elin(i) = 0.0
          end if
        else if ((resname.eq.'RPU').OR.(resname.eq.'RIU')) then
          if ((ii.ge.at(rs)%sc(8-shf2)).AND.(kk.ge.at(rs)%sc(8-shf2))) then
            fudge(rs)%elin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(8-shf2))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(8-shf2)))) then
            fudge(rs)%elin(i) = 0.0
          end if
        else if ((resname.eq.'RPC').OR.(resname.eq.'RIC')) then
          if ((ii.ge.at(rs)%sc(8-shf2)).AND.(ii.le.at(rs)%sc(17-shf2)).AND.&
 &            (kk.ge.at(rs)%sc(8-shf2)).AND.(kk.le.at(rs)%sc(17-shf2))) then
            fudge(rs)%elin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(8-shf2)).AND.(kk.le.at(rs)%sc(17-shf2))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(8-shf2)).AND.(ii.le.at(rs)%sc(17-shf2)))) then
            fudge(rs)%elin(i) = 0.0
          end if
        else if ((resname.eq.'RPG').OR.(resname.eq.'RIG')) then
          if ((ii.ge.at(rs)%sc(8-shf2)).AND.(ii.le.at(rs)%sc(20-shf2)).AND.&
 &            (kk.ge.at(rs)%sc(8-shf2)).AND.(kk.le.at(rs)%sc(20-shf2))) then
            fudge(rs)%elin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(8-shf2)).AND.(kk.le.at(rs)%sc(20-shf2))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(8-shf2)).AND.(ii.le.at(rs)%sc(20-shf2)))) then
            fudge(rs)%elin(i) = 0.0
          end if
        else if ((resname.eq.'DPA').OR.(resname.eq.'DIA')) then
          if ((ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(18-shf3)).AND.&
 &            (kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(18-shf3))) then
            fudge(rs)%elin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(18-shf3))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(18-shf3)))) then
            fudge(rs)%elin(i) = 0.0
          end if
        else if ((resname.eq.'DPT').OR.(resname.eq.'DIT')) then
          if ((ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(17-shf3)).AND.&
 &            (kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(17-shf3))) then
            fudge(rs)%elin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(17-shf3))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(17-shf3)))) then
            fudge(rs)%elin(i) = 0.0
          end if
        else if ((resname.eq.'DPC').OR.(resname.eq.'DIC')) then
          if ((ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(16-shf3)).AND.&
 &            (kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(16-shf3))) then
            fudge(rs)%elin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(16-shf3))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(16-shf3)))) then
            fudge(rs)%elin(i) = 0.0
          end if
        else if ((resname.eq.'DPG').OR.(resname.eq.'DIG')) then
          if ((ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(19-shf3)).AND.&
 &            (kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(19-shf3))) then
            fudge(rs)%elin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(19-shf3))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(19-shf3)))) then
            fudge(rs)%elin(i) = 0.0
          end if
        end if
      end if
      if ((resname.eq.'RPA').OR.(resname.eq.'RIA')) then
        if (((ii.eq.at(rs)%sc(12-shf2)).AND.((kk.eq.at(rs)%sc(20-shf2)).OR.(kk.eq.at(rs)%sc(21-shf2)))).OR.&
            ((kk.eq.at(rs)%sc(12-shf2)).AND.((ii.eq.at(rs)%sc(20-shf2)).OR.(ii.eq.at(rs)%sc(21-shf2))))) then
          fudge(rs)%elin(i) = 0.0
        end if
      else if ((resname.eq.'DPA').OR.(resname.eq.'DIA')) then
        if (((ii.eq.at(rs)%sc(11-shf3)).AND.((kk.eq.at(rs)%sc(19-shf3)).OR.(kk.eq.at(rs)%sc(20-shf3)))).OR.&
            ((kk.eq.at(rs)%sc(11-shf3)).AND.((ii.eq.at(rs)%sc(19-shf3)).OR.(ii.eq.at(rs)%sc(20-shf3))))) then
          fudge(rs)%elin(i) = 0.0
        end if
      else if (resname.eq.'RPU') then
!       fortunately, no one besides the GROMOS developers needs to understand this
!       the two excludes are for the O3* with respect to the immediately preceding O2* and to the following H5 on the uracil ring
        if (((ii.eq.at(rs)%bb(1)).AND.(kk.eq.at(rs)%sc(18-shf2))).OR.&
            ((kk.eq.at(rs)%bb(1)).AND.(ii.eq.at(rs)%sc(18-shf2)))) then
          fudge(rs)%elin(i) = 0.0
        end if
        if (rs.eq.rsmol(molofrs(rs),2)) then ! O3* vs. O2* exclusion
          if (((ii.eq.at(rs)%sc(4)).AND.(kk.eq.at(rs)%bb(9))).OR.&
            ((kk.eq.at(rs)%sc(4)).AND.(ii.eq.at(rs)%bb(9)))) then
            fudge(rs)%elin(i) = 0.0
          end if
        end if
      end if
    end do
    do i=1,nrpolnb(rs)
!
      resname = amino(seqtyp(rs))
      ii = iaa(rs)%polnb(i,1)
      kk = iaa(rs)%polnb(i,2)
      if ((attyp(ii).eq.0).OR.(attyp(kk).eq.0)) cycle
!
      if (resname.eq.'RPU') then
        if (((ii.eq.at(rs)%sc(4)).AND.(kk.eq.at(rs+1)%bb(1))).OR.&
            ((kk.eq.at(rs)%sc(4)).AND.(ii.eq.at(rs+1)%bb(1)))) then
          fudge(rs)%elnb(i) = 0.0
        end if
      end if
    end do
  end do
!
end
!
!
!--------------------------------------------------------------------------------------------
!
! get exclusion lists for polar interactions across cross-links
!
subroutine crosslink_excludes_polar()
!
  use iounit
  use sequen
  use polypep
  use atoms
  use inter
  use params
  use energies
!
  implicit none
!
  integer i,j,k,l,rs1,rs2,ii,kk,dp1,dp2
  logical proceed,its14
!
  do j=1,n_crosslinks
    crosslink(j)%nrspol = 0
    if ((crosslink(j)%itstype.eq.1).OR.(crosslink(j)%itstype.eq.2)) then
      rs1 = crosslink(j)%rsnrs(1)
      rs2 = crosslink(j)%rsnrs(2)
      allocate(crosslink(j)%exclpol(at(rs1)%npol*at(rs2)%npol,2))
      allocate(crosslink(j)%cbpars(at(rs1)%npol*at(rs2)%npol))
      allocate(crosslink(j)%is14pol(at(rs1)%npol*at(rs2)%npol))
!
      if (elec_model.eq.2) then
        do dp1=1,at(rs1)%ndpgrps
          do dp2=1,at(rs2)%ndpgrps
            proceed = .false.
            its14 = .false.
            do i=1,at(rs1)%dpgrp(dp1)%nats
              ii = at(rs1)%dpgrp(dp1)%ats(i)
              do k=1,at(rs2)%dpgrp(dp2)%nats
                kk = at(rs2)%dpgrp(dp2)%ats(k)
!
                do l=1,n12(ii)
                  if (i12(l,ii).eq.kk) then
                    proceed = .true.
                    its14 = .false.
                  end if
                end do
!
                do l=1,n13(ii)
                  if (i13(l,ii).eq.kk) then
                    proceed = .true.
                    its14 = .false.
                  end if
                end do
!
                do l=1,n14(ii)
                  if (use_14.EQV..false.) then
                    if (i14(l,ii).eq.kk) proceed = .true.
                  else
                    if (i14(l,ii).eq.kk) then
                      if (proceed.EQV..false.) its14 = .true.
                      proceed = .true.
                    end if
                  end if
                end do
              end do
            end do
!
            if (proceed.EQV..true.) then
              do i=1,at(rs1)%dpgrp(dp1)%nats
                ii = at(rs1)%dpgrp(dp1)%ats(i)
                do k=1,at(rs2)%dpgrp(dp2)%nats
                  kk = at(rs2)%dpgrp(dp2)%ats(k)
                  crosslink(j)%nrspol = crosslink(j)%nrspol + 1
                  crosslink(j)%exclpol(crosslink(j)%nrspol,1) = ii
                  crosslink(j)%exclpol(crosslink(j)%nrspol,2) = kk
                  if (its14.EQV..true.) then 
                    crosslink(j)%cbpars(crosslink(j)%nrspol) = (1.0-fudge_el_14)*atq(ii)*atq(kk)
                    crosslink(j)%is14pol(crosslink(j)%nrspol) = .true.
                  else
                    crosslink(j)%cbpars(crosslink(j)%nrspol) = atq(ii)*atq(kk)
                    crosslink(j)%is14pol(crosslink(j)%nrspol) = .false.
                  end if
                end do
              end do
            end if
          end do
        end do
!
      else if (elec_model.eq.1) then
        do i=1,at(rs1)%npol
          ii = at(rs1)%pol(i)
          do k=1,at(rs2)%npol
            kk = at(rs2)%pol(k)
            proceed = .false.
            its14 = .false.
!
            do l=1,n12(ii)
              if (i12(l,ii).eq.kk) proceed = .true.
            end do
!
            if (proceed.EQV..false.) then
              do l=1,n13(ii)
                if (i13(l,ii).eq.kk) proceed = .true.
              end do
            end if
!
            if (proceed.EQV..false.) then
              do l=1,n14(ii)
                if (use_14.EQV..false.) then
                  if (i14(l,ii).eq.kk) proceed = .true.
                else if (fudge_el_14.ge.1.0) then
!                  do nothing
                else
                  if (i14(l,ii).eq.kk) then
                    proceed = .true.
                    its14 = .true.
                  end if
                end if
              end do
            end if
!
            if (proceed.EQV..true.) then
              crosslink(j)%nrspol = crosslink(j)%nrspol + 1
              crosslink(j)%exclpol(crosslink(j)%nrspol,1) = ii
              crosslink(j)%exclpol(crosslink(j)%nrspol,2) = kk
              if (its14.EQV..true.) then 
                crosslink(j)%cbpars(crosslink(j)%nrspol) = (1.0-fudge_el_14)*atq(ii)*atq(kk)
                crosslink(j)%is14pol(crosslink(j)%nrspol) = .true.
              else
                crosslink(j)%cbpars(crosslink(j)%nrspol) = atq(ii)*atq(kk)
                crosslink(j)%is14pol(crosslink(j)%nrspol) = .false.
              end if
            end if
          end do
        end do
!
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported crosslink type in crosslink_excludes_polar(...).&
 & This is most likely an omission bug.'
      call fexit()
    end if
  end do
!
end
!
!----------------------------------------------------------------------------
!
!  THE ROUTINE TO SETUP THE DIPOLE-GROUP BASED SCREENED ATOMIC CHARGES 
!
!----------------------------------------------------------------------------
!
function sq_ipol(ati)
!
  use atoms
  use iounit
  use energies

  implicit none
!
  integer ati
  RTYPE sq_ipol,rat
!
  
  rat = atsav(ati)/atbvol(ati)

  if (rat.ge.atsavmaxfr(ati)) then
    sq_ipol = 1.0
    
  else if (rat.gt.par_IMPSOLV(5)) then
    sq_ipol = 1.0/(1.0 + exp(-(rat-atsavprm(ati,4))/par_IMPSOLV(4)))
    sq_ipol = (sq_ipol-0.5)*atsavprm(ati,5) + atsavprm(ati,6)
!    sq_ipol = (1.0/(atsavmaxfr(ati)-atsavminfr2(ati)))*
! &                      (rat - atsavminfr2(ati))
  else
    sq_ipol = 0.0
  end if

!
  return
!
end
!
!
subroutine setup_scrqs2()
!
  use iounit
  use aminos
  use sequen
  use polypep
  use energies
  use molecule
  use atoms
!
  implicit none
!
  integer i,ii,imol,rs
  RTYPE fs_gr,sq_ipol
!
! in model #4 (pure distance-dependent dielectric), we don't need
! environmentally screened charges
  if (scrq_model.eq.4) return
!
! in models #2,5,8,9 we do this purely by atom, no group-considerations at all
!peinrt 
  if ((scrq_model.eq.2).OR.(scrq_model.eq.6).OR.(scrq_model.eq.8)&
 &     .OR.(scrq_model.eq.9)) then
    do i=1,n
      if (atq(i).ne.0.0) then
        scrq(i) = (1.0 - coul_scr*sq_ipol(i))
      else
        scrq(i) = 0.0
      end if
    end do
    return
  end if
!
!
! in (the clean) models #1,3,5,7 we do charge-weighted, group-consistent screening
! note that this effectively couples many atoms' solvation states together (just
! like the model compound-derived FOS calculation does, although there there are
! many more empirical rules about which atoms contribute what to the FOS).
!
  do imol=1,nmol
!
    do rs=rsmol(imol,1),rsmol(imol,2)
!
      do i=1,at(rs)%ndpgrps
        fs_gr = 0.0
        do ii=1,at(rs)%dpgrp(i)%nats
          fs_gr = fs_gr + abs(atq(at(rs)%dpgrp(i)%ats(ii)))*&
 &                              sq_ipol(at(rs)%dpgrp(i)%ats(ii))
        end do
        fs_gr = fs_gr/at(rs)%dpgrp(i)%qnm
        do ii=1,at(rs)%dpgrp(i)%nats
          scrq(at(rs)%dpgrp(i)%ats(ii)) = (1.0-coul_scr*fs_gr)
        end do
!
      end do
!
    end do ! over residues
!
  end do ! over molecules
!
end
!
!-----------------------------------------------------------------------
!
! a subroutine to setup (G)RF electrostatics treatment 
!
subroutine grf_setup()
!
  use energies
  use iounit
  use cutoffs
  use units
  use system
  use math
  use sequen
  use polypep
  use atoms
!
  implicit none
!
  RTYPE istr,epsr,t1,t2,kap,netQ
  integer i,rs,k
!
  if (use_POLAR.EQV..false.) then
    lrel_md = 1
    return
  end if
!
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    write(ilog,*) 'Fatal. Reaction-field method is currently not sup&
 &ported for ensembles with fluctuating volumes. Check back later.'
    call fexit()
  end if
!
  if (use_FEG.EQV..true.) then
    if (fegcbmode.ne.1) then
      write(ilog,*) 'Fatal. Reaction-field method is currently not s&
 &upported for free energy growth calculations with non-linear Coulo&
 &mb scaling. Check back later.'
      call fexit()
    end if
    if (scale_POLAR.ne.1.0) then
      write(ilog,*) 'Fatal. Reaction-field method is currently not s&
 &upported for background Hamiltonians with a scaled Coulomb term (i&
 &.e., FMCSC_SC_POLAR must be 1.0). Check back later.'
      call fexit()
    end if
    k = 0
    do i=1,nseq
      if (par_FEG(i).EQV..true.) then
        k = k + 1
      end if
      if (k.gt.1) then
!        write(ilog,*) 'Fatal. The reaction-field method in conjuncti&
! &n with simple linear Coulomb-scaling is only supported for a singl&
! &e residue to be grown in at the moment. Please check back later.'
!        call fexit()
      end if
    end do
  end if
!
  if (use_IMPSOLV.EQV..true.) then
    write(ilog,*) 'Fatal. Reaction-field method is currently not sup&
 &ported for the ABSINTH implicit solvent model. Check back later.'
    call fexit()
  end if
!
  if ((is_prflj.EQV..false.).AND.(is_fegprflj.EQV..false.)) then
    write(ilog,*) 'Fatal. Current Hamiltonian does not support react&
 &ion-field method. Please check back later.'
    call fexit()
  end if
!
  if (use_cutoffs.EQV..false.) then
    write(ilog,*) 'Fatal. Reaction-field method requires using and s&
 &etting a finite cutoff.'
    call fexit()
  end if
!
  if ((nbl_up.gt.1).AND.(mcel_cutoff.gt.mcnb_cutoff)) then
    write(ilog,*) 'Fatal. The use of reaction-field methods is not compatible with the twin-range &
 &cutoff setup. Either set neighbor list updates to 1 (FMCSC_NBL_UP) or set cutoffs to the same values &
 &(FMCSC_NBCUTOFF and FMCSC_ELCUTOFF).'
    call fexit()
  end if

  istr= 0.0
!
! we have at least three choices here:
! 1) ignore actual system and assign ionic strength
! 2) calculate ionic strength from all charge groups as independent ions
! 3) calculate ionic strength from all charged molecules as (macro)ions
! we go with 2) in general, since 1) is arbitrary and potentially inconsistent, and
! since 3) is typically a poor assumption as macroions are present in very low copy number (typically 1)
! and hence concentration and ionic strength are ill-defined
!
! other than that, we can only assume charges are more or less homogeneously distributed
! also note that fractional charges really break the idea of ionic strength
  do i=1,cglst%ncs
    istr = istr + cglst%tc(i)*cglst%tc(i)
  end do
  epsr = par_IMPSOLV(2)
  netQ = 0.0
  do rs=1,nseq
    do i=1,at(rs)%npol
      netQ = netQ + atq(at(rs)%pol(i))*atq(at(rs)%pol(i))
    end do
  end do
!
! standard RF (GRF limiting case with zero ionic strength)
  if (rf_mode.eq.2) then
    par_POLAR(1) = (epsr-1.0)/((2.0*epsr+1.0)*(mcel_cutoff**3.0))
    par_POLAR(2) = 1.0/mcel_cutoff + mcel_cutoff2*par_POLAR(1)
! GRF
  else
    kap = sqrt(4.0*PI*electric*invtemp*istr/ens%insV)
    t1 = 1.0 + kap*mcel_cutoff
    t2 = kap*kap*mcel_cutoff2
    par_POLAR(1) = (t1*(epsr-1.0) + 0.5*epsr*t2)/&
 &        (((2.0*epsr+1.0)*t1 + epsr*t2)*(mcel_cutoff**3.0))
    par_POLAR(2) = 1.0/mcel_cutoff + mcel_cutoff2*par_POLAR(1)
  end if
  par_POLAR(3) = istr/ens%insV
  rfcnst = -0.5*electric*scale_POLAR*netQ*par_POLAR(2)
!  write(*,*) istr,kap,par_POLAR(1),par_POLAR(2)
!
end
!
!----------------------------------------------------------------------------------
!
! a subroutine which computes and writes dipole moments to target vectors
!
subroutine get_dipoles(imol,rscn,rsdcnt,molcn,moldcnt)
!
  use molecule
  use atoms
  use polypep
  use sequen
!
  implicit none
!
  integer i,imol,rs,rsdcnt(rsmol(imol,2)-rsmol(imol,1)+1),moldcnt
  RTYPE rscn(rsmol(imol,2)-rsmol(imol,1)+1,3),molcn(3)
  RTYPE moltc,eps
  RTYPE, ALLOCATABLE:: rstc(:)
!
  allocate(rstc(rsmol(imol,2)-rsmol(imol,1)+1))
!
! initialize
  molcn(:) = 0.0
  moltc = 0.0
  moldcnt = 0
  eps = (atmol(imol,2)-atmol(imol,1)+1.0)*(10.0**(-precision(eps)))
  rstc(:) = 0.0
  rscn(:,:) = 0.0
  rsdcnt(:) = 0
!
! compute com
  do rs=rsmol(imol,1),rsmol(imol,2)
    do i=1,at(rs)%npol
      rstc(rs-rsmol(imol,1)+1) = rstc(rs-rsmol(imol,1)+1)&
 &                                     + atq(at(rs)%pol(i))
      moltc = moltc + atq(at(rs)%pol(i))
    end do
  end do
!
  if (abs(moltc).le.eps) then
    moldcnt = moldcnt + 1
  end if
  do rs=rsmol(imol,1),rsmol(imol,2)
    if (abs(rstc(rs-rsmol(imol,1)+1)).le.eps) then
      rsdcnt(rs-rsmol(imol,1)+1) = rsdcnt(rs-rsmol(imol,1)+1) + 1
      do i=1,at(rs)%npol
        rscn(rs-rsmol(imol,1)+1,1) = rscn(rs-rsmol(imol,1)+1,1)&
 & + atq(at(rs)%pol(i))*x(at(rs)%pol(i))
        rscn(rs-rsmol(imol,1)+1,2) = rscn(rs-rsmol(imol,1)+1,2)&
 & + atq(at(rs)%pol(i))*y(at(rs)%pol(i))
        rscn(rs-rsmol(imol,1)+1,3) = rscn(rs-rsmol(imol,1)+1,3)&
 & + atq(at(rs)%pol(i))*z(at(rs)%pol(i))
      end do
    end if
    if (abs(moltc).le.eps) then
      do i=1,at(rs)%npol
        molcn(1) = molcn(1) + atq(at(rs)%pol(i))*x(at(rs)%pol(i))
        molcn(2) = molcn(2) + atq(at(rs)%pol(i))*y(at(rs)%pol(i))
        molcn(3) = molcn(3) + atq(at(rs)%pol(i))*z(at(rs)%pol(i))
      end do
    end if
  end do
!
  deallocate(rstc)
!
end
!
!-----------------------------------------------------------------------
!
subroutine do_dipoles()
!
  use molecule
  use atoms
  use polypep
  use sequen
  use dipolavg
  use energies
  use grandensembles
  use system
  use pdb
!
  implicit none
!
  integer imol,rs,mt,k
  RTYPE, ALLOCATABLE:: rscn(:,:)
  RTYPE molcn(3),val,scal
  integer, ALLOCATABLE:: rsdcnt(:)
  integer moldcnt
  logical ismember
!
  if (use_POLAR.EQV..false.) return
!
  if (use_frame_weights.EQV..true.) then
    scal = 1.0*framewts(curframe)
  else
    scal = 1.0
  end if
!
  do imol=1,nmol
    mt = an_grp_mol(imol)
    if (molischg(moltypid(imol)).EQV..false.) cycle
!
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
!
    allocate(rscn(rsmol(imol,2)-rsmol(imol,1)+1,3))
    allocate(rsdcnt(rsmol(imol,2)-rsmol(imol,1)+1))
    call get_dipoles(imol,rscn,rsdcnt,molcn,moldcnt)
    if (moldcnt.gt.0) then
      nmtdavg(mt) = nmtdavg(mt) + moldcnt
      val = sqrt(molcn(1)**2+molcn(2)**2+molcn(3)**2)
      mtdavg(mt,1) = mtdavg(mt,1) + scal*val
      mtdavg(mt,2) = mtdavg(mt,2) + scal*val*val
      mtdavg(mt,3) = mtdavg(mt,3) + scal*molcn(1)
      mtdavg(mt,4) = mtdavg(mt,4) + scal*molcn(2)
      mtdavg(mt,5) = mtdavg(mt,5) + scal*molcn(3)
    end if
    do rs=rsmol(imol,1),rsmol(imol,2)
      k = rs-rsmol(imol,1)+1
      if (rsdcnt(k).gt.0) then
        nrsdavg(rs) = nrsdavg(rs) + rsdcnt(k)
        val = sqrt(rscn(k,1)**2+rscn(k,2)**2+rscn(k,3)**2)
        rsdavg(rs,1) = rsdavg(rs,1) + scal*val
        rsdavg(rs,2) = rsdavg(rs,2) + scal*val*val
        rsdavg(rs,3) = rsdavg(rs,3) + scal*rscn(k,1)
        rsdavg(rs,4) = rsdavg(rs,4) + scal*rscn(k,2)
        rsdavg(rs,5) = rsdavg(rs,5) + scal*rscn(k,3)
      end if
    end do 
!    write(*,*) debye*sqrt(molcn(1)**2+molcn(2)**2+molcn(3)**2)
    deallocate(rsdcnt)
    deallocate(rscn)
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine prt_dipoles()
!
  use molecule
  use atoms
  use polypep
  use sequen
  use dipolavg
  use energies
  use units
  use iounit
  use mpistuff
  use pdb
  use system
!
  implicit none
!
  integer rs,mt,effrs
  integer i,k,iu,freeunit,effmoltyp,ii,jj
  character(100) fn
  logical exists
#ifdef ENABLE_MPI
  integer tl
  character(3) nod
#endif
  RTYPE, ALLOCATABLE:: normer(:)
!
  allocate(normer(nangrps))
!
  if (use_frame_weights.EQV..true.) then
    normer(:) = 0.0
    k = 0
    do i=1,framecnt
      if ((framelst2(i).gt.nequil).AND.(mod(framelst2(i),dipcalc).eq.0)) then
        k = k + 1
        normer(1:nangrps) = normer(1:nangrps) + molangr(1:nangrps,2)*framewts2(i)
      end if
    end do
    if (k*molangr(1,2).ne.nmtdavg(1)) then
      write(ilog,*) 'Warning. Dipole moment averages have inconsistent normalization. This is &
 &a bug. Please report this error and ignore all output in MOLDIPOLES.dat or related files.'
    end if
  else
    normer(1:nangrps) = 1.0*nmtdavg(1:nangrps)
  end if

  effmoltyp = 0
  do mt=1,nangrps
    if (nmtdavg(mt).gt.0) then
      effmoltyp = effmoltyp + 1
      do i=1,DIPDIM
        mtdavg(mt,i) = mtdavg(mt,i)/normer(mt)
      end do
      mtdavg(mt,1) = debye*mtdavg(mt,1)
      mtdavg(mt,2) = debye*debye*mtdavg(mt,2)
    end if
  end do
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    tl = 3
    call int2str(myrank,nod,tl)
    fn = 'N_'//nod(1:tl)//'_MOLDIPOLES.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'MOLDIPOLES.dat'
  end if
#else
  fn = 'MOLDIPOLES.dat'
#endif
  call strlims(fn,ii,jj)
!
 24   format('# Analysis group ',i3,' (Ref.: Mol. ',i5,' / Res. ',i5,'-',&
 &i5,'):')
 89   format(100(g14.7,1x))
 88   format(i10,1x,100(g14.7,1x))
!
  if (effmoltyp.gt.0) then
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      iu = freeunit()
      open(unit=iu,file=fn(ii:jj),status='old')
      close(unit=iu,status='delete')
    end if
    iu=freeunit()
    open(unit=iu,file=fn(ii:jj),status='new')
!
    do mt=1,nangrps
      if (nmtdavg(mt).gt.0) then
        write(iu,24) mt,molangr(mt,1),rsmol(molangr(mt,1),1),rsmol(molangr(mt,1),2)
        write(iu,89) mtdavg(mt,1),mtdavg(mt,2),sqrt(max((mtdavg(mt,2)-mtdavg(mt,1)**2),0.0d0)),&
 &mtdavg(mt,3),mtdavg(mt,4),mtdavg(mt,5)
      end if
    end do
  end if
!
!
  effrs = 0
  do rs=1,nseq
    if (nrsdavg(rs).gt.0) then
      effrs = effrs + 1
      if (use_frame_weights.EQV..true.) then
       rsdavg(rs,1:DIPDIM) = 1.0*molangr(an_grp_mol(molofrs(rs)),2)*rsdavg(rs,1:DIPDIM)/normer(an_grp_mol(molofrs(rs)))
      else
       rsdavg(rs,1:DIPDIM) = rsdavg(rs,1:DIPDIM)/(1.0*nrsdavg(rs))
      end if
      rsdavg(rs,1) = debye*rsdavg(rs,1)
      rsdavg(rs,2) = debye*debye*rsdavg(rs,2)
    end if
  end do
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    tl = 3
    call int2str(myrank,nod,tl)
    fn = 'N_'//nod(1:tl)//'_RESDIPOLES.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'RESDIPOLES.dat'
  end if
#else
  fn = 'RESDIPOLES.dat'
#endif
!
  call strlims(fn,ii,jj)
  if (effrs.gt.0) then
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      iu = freeunit()
      open(unit=iu,file=fn(ii:jj),status='old')
      close(unit=iu,status='delete')
    end if
    iu=freeunit()
    open(unit=iu,file=fn(ii:jj),status='new')
!
    do rs=1,nseq
      if (nrsdavg(rs).gt.0) then
        write(iu,88) rs,rsdavg(rs,1),rsdavg(rs,2),&
 &sqrt(max((rsdavg(rs,2)-rsdavg(rs,1)**2),0.0d0)),&
 &rsdavg(rs,3),rsdavg(rs,4),rsdavg(rs,5)
      end if
    end do
  end if
!
  deallocate(normer)
!
end
!
!-----------------------------------------------------------------------