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
! obtain parameters from file
!
subroutine readprm()
!
  use iounit
  use params
  use keys
  use movesets
  use energies
  use sequen
!
  implicit none
!
  integer i,j,iprm,next,t1,t2,maxct,maxbt,maxljt
  integer maxbot,maxant,maxdit,maxcmt
  integer ia,ib,ic,id,ie,ig,aone
  integer kbo,kba,kdi,kdi2,kcm
  integer iomessage,freeunit
  integer, ALLOCATABLE:: tmplst(:)
  RTYPE rd,red2
  character(MAXKWLEN) keyword
  character(MAXKEYLEN) record,string,str2
  logical exists,newbo
  integer, allocatable:: b_types_used(:)
!
!
! some initialization
!
  aone = 1
!     
! set squared reduction factor (HSSCALE)
!     
  red2 = reduce**2
!
! process each line of the parameter file, first
! extract the keyword at the start of each line
!
  call strlims(paramfile,t1,t2)
  inquire(file=paramfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'Fatal. Cannot open specified parameter file (',&
 &paramfile(t1:t2),'). Wrong path?'
    call fexit()
  end if
  iprm = freeunit()
  open (unit=iprm,file=paramfile(t1:t2),status='old')
!
 50   format(a200)
!
! first a dry run to establish parameter dimensions
!
  maxbt = 0
  maxljt = 0
  maxct = 0
  maxbot = 0
  maxant = 0
  maxdit = 0
  maxcmt = 0
  maxcmpar = 10
  n_biotyp = 0
  n_ljtyp = 0
  n_ctyp = 0
  n_bondtyp = 0
  n_angltyp = 0
  n_torstyp = 0
  n_cmstyp = 0
  ba_lstsz =0 
  bo_lstsz = 0
  di_lstsz = 0
  cm_lstsz = 0
  impt_lstsz = 0
  do while (.true.)
!
    read (iprm,50,iostat=iomessage)  record
    if (iomessage.eq.-1) then
      exit
    else if (iomessage.eq.2) then
      write(ilog,*) 'Fatal. I/O error while reading parameter file (&
 &got: ',paramfile(t1:t2),').'
      call fexit()
    end if
    next = 1
    ia = 0
    call extract_str(record,keyword,next)
    call toupper(keyword)
    if (keyword(1:5) .eq. 'ATOM ') then
      call extract_int(record,ia,next)
      if (ia.le.0) then
        write(ilog,*) 'Bad parameter file. LJ classes are positive i&
 &ntegers (got ',ia,'). File is ',paramfile(t1:t2),'.'
        call fexit()
      end if 
      if (ia.gt.maxljt) maxljt = ia
      n_ljtyp = n_ljtyp + 1
    else if (keyword(1:7).eq.'CHARGE ') then
      call extract_int(record,ia,next)
      if (ia.le.0) then
        write(ilog,*) 'Bad parameter file. Charge types are positive&
 & integers (got ',ia,'). File is ',paramfile(t1:t2),'.'
        call fexit()
      end if
      if (ia.gt.maxct) maxct = ia
      n_ctyp = n_ctyp + 1
    else if (keyword(1:8).eq.'BIOTYPE ') then
      call extract_int(record,ia,next)
      if (ia.le.0) then
        write(ilog,*) 'Bad parameter file. Biotypes are positive int&
 &egers (got ',ia,'). File is ',paramfile(t1:t2),'.'
        call fexit()
      end if 
      if (ia.gt.maxbt) maxbt = ia
      n_biotyp = n_biotyp + 1
    else if (keyword(1:5).eq.'BOND ') then
      call extract_int(record,ia,next)
      if (ia.le.0) then
        write(ilog,*) 'Bad parameter file. Bond entries are positive&
 & integers (got ',ia,'). File is ',paramfile(t1:t2),'.'
        call fexit()
      end if 
      if (ia.gt.maxbot) maxbot = ia
      n_bondtyp = n_bondtyp + 1
    else if (keyword(1:6).eq.'ANGLE ') then
      call extract_int(record,ia,next)
      if (ia.le.0) then
        write(ilog,*) 'Bad parameter file. Angle entries are positiv&
 &e integers (got ',ia,'). File is ',paramfile(t1:t2),'.'
        call fexit()
      end if 
      if (ia.gt.maxant) maxant = ia
      n_angltyp = n_angltyp + 1
    else if (keyword(1:8).eq.'TORSION ') then
      call extract_int(record,ia,next)
      if (ia.le.0) then
        write(ilog,*) 'Bad parameter file. Torsion entries are posit&
 &ive integers (got ',ia,'). File is ',paramfile(t1:t2),'.'
        call fexit()
      end if 
      if (ia.gt.maxdit) maxdit = ia
      n_torstyp = n_torstyp + 1
    else if (keyword(1:5).eq.'CMAP ') then
      call extract_int(record,ia,next)
      if (ia.le.0) then
        write(ilog,*) 'Bad parameter file. CMAP entries are posit&
 &ive integers (got ',ia,'). File is ',paramfile(t1:t2),'.'
        call fexit()
      end if
      if (ia.gt.maxcmt) maxcmt = ia
      n_cmstyp = n_cmstyp + 1
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) ic,ig
      if (ig.gt.maxcmpar) maxcmpar = ig
    end if
  end do
!
! note that not all types are really needed for running the software, here we
! merely enforce the existence of at least one vanilla entry
! partially this is so we don't have to re-check for existence during deallocation
 47   format(' Number of biotypes defined         : ',i4/&
 &       ' Number of LJ-types defined         : ',i4/&
 &       ' Number of charge types defined     : ',i4/&
 &       ' Number of bond pot. types defined  : ',i4/&
 &       ' Number of angle pot. types defined : ',i4/&
 &       ' Number of tors. pot. types defined : ',i4/&
 &       ' Number of CMAP types defined       : ',i4/)
  if ((n_biotyp.le.0).OR.(n_ctyp.le.0).OR.&
 &    (n_ljtyp.le.0).OR.(n_bondtyp.le.0).OR.(n_angltyp.le.0).OR.&
 &    (n_torstyp.le.0)) then
    write(ilog,*)'Incomplete parameter file (',paramfile(t1:t2),').'
    write(ilog,47) n_biotyp,n_ljtyp,n_ctyp,n_bondtyp,n_angltyp,&
 &n_torstyp,n_cmstyp
    call fexit()
  end if
!
  if ((n_biotyp.ne.maxbt).OR.(n_ctyp.ne.maxct).OR.&
 &    (n_ljtyp.ne.maxljt).OR.(n_bondtyp.ne.maxbot).OR.&
 &    (n_angltyp.ne.maxant).OR.(n_torstyp.ne.maxdit).OR.(n_cmstyp.ne.maxcmt)) then
    write(ilog,*) 'Bad parameter file. All types (keywords "ATOM", "&
 &CHARGE", "BIOTYPE", "BOND", "ANGLE", "TORSION", and "CMAP") should be cont&
 &inuously numbered (',paramfile(t1:t2),').'
    write(ilog,*) 'BIOTYPE',(n_biotyp.ne.maxbt)
    write(ilog,*) 'CHARGE',(n_ctyp.ne.maxct)
    write(ilog,*) 'ATOM',(n_ljtyp.ne.maxljt)
    write(ilog,*) 'BOND',(n_bondtyp.ne.maxbot)
    write(ilog,*) 'ANGLE',(n_angltyp.ne.maxant)
    write(ilog,*) 'TORSION',(n_torstyp.ne.maxdit)
    write(ilog,*) 'CMAP',(n_cmstyp.ne.maxcmt)
    call fexit()
  end if
!
  if ((do_pka_2.eqv..true.).or.(do_hsq.eqv..true.)) then 
    allocate(transform_table(n_biotyp))
    do i=1, n_biotyp ! Martin : note that while adding elements to this list, erros will appear further downstream of the pregram 
                    ! if one of the element refers to an element that itself is tramformed 
        if ((i.lt.1208).OR.((i.ge.1251).AND.(i.le.1291))) then 
            transform_table(i)=i 
        else if (i.eq.1208) then 
            transform_table(i)=1208
        else if (i.eq.1209) then ! Glutamic acid titrable proton
            transform_table(i)=1087
        else if (i.eq.1210) then 
            transform_table(i)=1085
        else if (i.eq.1211) then 
            transform_table(i)=1086
        else if (i.eq.1212) then 
            transform_table(i)=1082 !End GLX
        else if (i.eq.1213) then 
            transform_table(i)=1084
        else if (i.eq.1214) then 
            transform_table(i)=1216 !CL
        else if (i.eq.1215) then 
            transform_table(i)=1223 !NA
        else if (i.eq.1216) then ! Self for salt off
            transform_table(i)=1216
        else if (i.eq.1217) then ! Lysine CE
            transform_table(i)=1112
        else if (i.eq.1218) then ! Lysine HE
            transform_table(i)=1113
        else if (i.eq.1219) then !Lysine NZ
            transform_table(i)=1114
        else if (i.eq.1220) then ! Lysine HZ
            transform_table(i)=1115
        else if (i.eq.1221) then ! Lysine HZ2 ! That is the one that should disappear
            transform_table(i)=1208!1208
        else if (i.eq.1222) then ! Lysine
            transform_table(i)=1083 ! Glutamic acid HG
        else if (i.eq.1223) then 
            transform_table(i)=1223    ! Salt off
        else if (i.eq.1224) then 
            transform_table(i)=1094
        else if (i.eq.1225) then 
            transform_table(i)=1096
        else if (i.eq.1226) then 
            transform_table(i)=1097 ! OD1 
        else if (i.eq.1227) then 
            transform_table(i)=1098 !OD2
        else if (i.eq.1228) then ! Aspartic acid titrable proton 
            transform_table(i)=1099 ! Martin  : changed again+
            ! From here on histidine tautomers
        else if (i.eq.1229) then 
            transform_table(i)=196 
        else if (i.eq.1230) then 
            transform_table(i)=197
        else if (i.eq.1231) then 
            transform_table(i)=198
        else if (i.eq.1232) then 
            transform_table(i)=199
        else if (i.eq.1233) then 
            transform_table(i)=1208 ! Martin : Can only reuse dummy atoms if the LJ param are interpolated from the same 
        else if (i.eq.1234) then 
            transform_table(i)=200
        else if (i.eq.1235) then 
            transform_table(i)=201
        else if (i.eq.1236) then 
            transform_table(i)=202
        else if (i.eq.1237) then 
            transform_table(i)=203
        else if (i.eq.1238) then 
            transform_table(i)=204
        else if (i.eq.1239) then 
            transform_table(i)=205 ! Martin : need to put a dummy here. Martin From now staring HDX
        else if (i.eq.1240) then 
            transform_table(i)=180 
        else if (i.eq.1241) then 
            transform_table(i)=181
        else if (i.eq.1242) then 
            transform_table(i)=182
        else if (i.eq.1243) then 
            transform_table(i)=183
        else if (i.eq.1244) then 
            transform_table(i)=184
        else if (i.eq.1245) then 
            transform_table(i)=185
        else if (i.eq.1246) then 
            transform_table(i)=186
        else if (i.eq.1247) then 
            transform_table(i)=187
        else if (i.eq.1248) then 
            transform_table(i)=188
        else if (i.eq.1249) then 
            transform_table(i)=189
        else if (i.eq.1250) then 
            transform_table(i)=1208
        ! Martin : start of the phosphoresidues    
        else if (i.eq.1292) then            ! Phosphoserine
            transform_table(i)=1165
        else if (i.eq.1293) then 
            transform_table(i)=1166
        else if (i.eq.1294) then 
            transform_table(i)=1167
        else if (i.eq.1295) then      !SXP P
            transform_table(i)=1168
        else if (i.eq.1296) then 
            transform_table(i)=1169
        else if (i.eq.1297) then 
            transform_table(i)=1170
        else if (i.eq.1298) then 
            transform_table(i)=1171
        else if (i.eq.1299) then ! Martin : phosphotyrosine
            transform_table(i)=1153 
        else if (i.eq.1300) then 
            transform_table(i)=1154
        else if (i.eq.1301) then 
            transform_table(i)=1155
        else if (i.eq.1302) then 
            transform_table(i)=1156
        else if (i.eq.1303) then 
            transform_table(i)=1157
        else if (i.eq.1304) then 
            transform_table(i)=1158  
        else if (i.eq.1305) then  ! phosphothreonine
            transform_table(i)=1178
        else if (i.eq.1306) then 
            transform_table(i)=1179
        else if (i.eq.1307) then 
            transform_table(i)=1180
        else if (i.eq.1308) then 
            transform_table(i)=1183
        else if (i.eq.1309) then 
            transform_table(i)=1184
        else if (i.eq.1310) then 
            transform_table(i)=1185
        else if (i.eq.1311) then 
            transform_table(i)=1186
        else if (i.eq.1312) then 
            transform_table(i)=129
        else if (i.eq.1313) then 
            transform_table(i)=130
        else if (i.eq.1314) then 
            transform_table(i)=131
        else if (i.eq.1315) then 
            transform_table(i)=132
        else if (i.eq.1316) then 
            transform_table(i)=133
        end if         
    end do
    
    allocate(lj_transform_table(n_ljtyp))
    ! Martin : The only purpose of this map is to get the equivalencies between type to overwrtie the 
    !special sigma and epsilon specified in the param file 
    ! This means that we need to interpolate the sigma and epsilon parameters, but not turn them to 0, since this is done by the 
    !replica parameters
    
    do i=1,n_ljtyp
        lj_transform_table(i)=i
    end do

    ! Martin : So the goal of this table is to get all the same special rules that apply to the LJ types are also applying exactly
    ! similarly to the titrable LJ types, without having to respecify them in the parameter files 
    ! That means this number should be set to whatever the charged states version of the non titrable LJ 
    
    lj_transform_table(31)=31 ! Martin : also equivalent to themselves, because can't set values to 0
    lj_transform_table(32)=32 !
    lj_transform_table(33)=33 !
    lj_transform_table(34)=4  ! Disapearing hydrogen
    lj_transform_table(35)=5  ! 29 to 6
    lj_transform_table(36)=13 ! Again 
    lj_transform_table(37)=12 !
    lj_transform_table(38)=5  !
    lj_transform_table(39)=29 !
    lj_transform_table(40)=8  !
    lj_transform_table(41)=1  ! Two different limits for the nitrogen atoms in order to have two different limits for histidne,                         
    lj_transform_table(42)=1  ! may not be necessary
    lj_transform_table(43)=11
    
    allocate(his_eqv_table(n_biotyp)) !This is define for all residues only to make code easier to read and understand
    
    ! Martin : changed the strategy : his_eqv is directly the third limit now 
    do i=1,n_biotyp
        his_eqv_table(i)=i
        if (i.eq.1229) then 
            his_eqv_table(i)=180
        else if (i.eq.1231) then 
            his_eqv_table(i)=182
        else if (i.eq.1232) then 
            his_eqv_table(i)=183
        else if (i.eq.1233) then 
            his_eqv_table(i)=184
        else if (i.eq.1234) then 
            his_eqv_table(i)=185
        else if (i.eq.1235) then 
            his_eqv_table(i)=186
        else if (i.eq.1236) then 
            his_eqv_table(i)=187            
        else if (i.eq.1237) then 
            his_eqv_table(i)=188
        else if (i.eq.1238) then 
            his_eqv_table(i)=189
        else if (i.eq.1239) then 
            his_eqv_table(i)=1208
        end if      
    end do 
  end if
  
  
  call allocate_params(aone)
  
  rewind(unit=iprm)
!
! now do the actual read-in
!
  do while (.true.)
!
    ia = 0
    read (iprm,50,iostat=iomessage) record
    if (iomessage.eq.-1) then
      exit
    else if (iomessage.eq.2) then
      write(ilog,*) 'Fatal. I/O error while reading parameter file (&
 &got: ',paramfile(t1:t2),').'
      call fexit()
    end if
    next = 1
    call extract_str(record,keyword,next)
    call toupper(keyword)
!
! LJ atom types (i.e., often real chemical types) 
!
    if (keyword(1:5) .eq. 'ATOM ') then
      call extract_int(record,ia,next)
      call extract_str(record,lj_symbol(ia),next)
      call extract_quo(record,lj_describe(ia),next)
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) lj_atnum(ia),lj_weight(ia),&
 &                                      lj_val(ia)
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
!
! group-based free energies of solvation
!
    else if (keyword(1:4).eq.'FOS ') then
      string = record(next:MAXKEYLEN)
      call read_fos(string)
!     
! charge type definitions
!     
    else if (keyword(1:7) .eq. 'CHARGE ') then
      call extract_int(record,ia,next)
      call extract_quo(record,c_note(ia),next)
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) c_charge(ia)
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
!     
! biotype definitions
!     
    else if (keyword(1:8) .eq. 'BIOTYPE ') then
      call extract_int(record,ia,next)
      call extract_abc(record,bio_code(ia),next)
      call extract_quo(record,bio_res(ia),next)
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) bio_ljtyp(ia),bio_ctyp(ia),&
 &bio_botyp(ia)
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
!
! bond definitions
!     
    else if (keyword(1:5) .eq. 'BOND ') then
      call extract_int(record,ia,next)
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) bo_typ(ia),&
 &                                (bo_par(ia,j),j=1,MAXBOPAR)
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
!
! bond definitions
!     
    else if (keyword(1:6) .eq. 'ANGLE ') then
      call extract_int(record,ia,next)
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) ba_typ(ia),&
 &                                (ba_par(ia,j),j=1,MAXBAPAR)
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
!
! torsion definitions
!     
    else if (keyword(1:8) .eq. 'TORSION ') then
      call extract_int(record,ia,next)
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) di_typ(ia),&
 &                                (di_par(ia,j),j=1,MAXDIPAR)
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
!
! CMAP definitions
!     
    else if (keyword(1:8) .eq. 'CMAP ') then
      call extract_int(record,ia,next)
      string = record(next:MAXKEYLEN)
      next = 1
      call extract_int(string,cm_typ(ia,1),next)
      str2 = string(next:MAXKEYLEN)
      next = 1
      call extract_int(str2,cm_typ(ia,2),next)
      cm_file(ia) = str2(next:MAXKEYLEN)
!      read (string,*,iostat=iomessage) cm_typ(ia,1),cm_typ(ia,2)
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
      if (cm_typ(ia,2).lt.4) then
        write(ilog,*) 'Fatal. CMAP potential specifications have to &
 &have at least four bins (in ',paramfile(t1:t2),').'
        call fexit()
      end if
      cm_par2(ia) = 360.0/(1.0*cm_typ(ia,2))
!
    end if
  end do
  
! read in CMAPs (separate files)
  call read_cmaps()
!
! some sanity checks
!
  do i=1,n_ljtyp
    if ((lj_weight(i).lt.0.0).OR.(lj_atnum(i).lt.0).OR.&
 &      (lj_val(i).lt.0)) then
      write(ilog,*) 'Fatal. Incomplete or disallowed parameters for &
 &LJ-type #',i,'. Multi-usage of a type in ',paramfile(t1:t2),'?'
      call fexit()
    end if
  end do


  allocate(tmplst(n_biotyp))
  n_botyp = 0
  do i=1,n_biotyp
!    write(*,*) i,bio_ljtyp(i),bio_ctyp(i)
    if ((bio_ljtyp(i).lt.0).OR.(bio_ctyp(i).lt.0)&
 &  .OR.(bio_botyp(i).lt.0)) then
      write(ilog,*) 'Fatal. Incomplete or disallowed parameters for &
 &biotype #',i,'. Multi-usage of a type in ',paramfile(t1:t2),'?'
      call fexit()
    end if
    
    
    
    newbo = .true.
    if (bio_botyp(i).eq.0) cycle
    do j=1,n_botyp
      if (bio_botyp(i).eq.tmplst(j)) then
        newbo = .false.
        exit
      end if
    end do
    if (newbo.EQV..true.) then
        
      n_botyp = n_botyp + 1
      tmplst(n_botyp) = bio_botyp(i)
    end if
  end do
  allocate(b_types_used(n_botyp+100))
  b_types_used(:)=0
  do i=1,n_biotyp
      if (bio_botyp(i).ne.0) then 
            b_types_used(bio_botyp(i))=1
      end if 
  end do 
   444   format(A8,I4,A,I3,A,I3,A,I3)
   do i=1,n_botyp
       
       if (b_types_used(i).eq.0) then 
           n_biotyp=n_biotyp+1
           write(ilog, 444 ) "biotype ",n_biotyp,'    XX " DUMMY     "',0,"  ",0,"  ",i
       end if 
   end do 
  
  do i=1,n_botyp
    if (tmplst(i).gt.n_botyp) then
      write(ilog,*) 'Fatal. Bonded types should be continuously numbe&
 &red, positive integers (in ',paramfile(t1:t2),').'
      call fexit()
    end if
  end do
!     
! if LJ-summary is requested, start providing that
!
 241  format('****   ',a8,': Arithmetic')
 242  format('****   ',a8,': Geometric')
 243  format('****   ',a8,': Harmonic')
 244  format('Distance for types ',i4,' and ',i4,': ',f5.2,' *',f8.4,&
 &' A')
 245  format('Strength for types ',i4,' and ',i4,': ',f8.4,' kcal/mol')
 246  format('Type ',i4,': Sigma   = ',f5.2,' *',f8.4,' A')
 247  format('         : Epsilon = ',f8.4,' kcal/mol')
 248  format('Type ',i4,': Radius  = ',f5.2,' A')
 249  format('Distance (14 only) for types ',i4,' and ',i4,': ',f5.2,' *',f8.4,&
 &' A')
 250  format('Strength (14 only) for types ',i4,' and ',i4,': ',f8.4,' kcal/mol')
!
  if (vdW_report.EQV..true.) then
    write(ilog,*)
    write(ilog,*) '--- Summary of vdW Parameters ---'
    write(ilog,*)
    write(ilog,*) 'Combination rules:'
    if (sigrule.eq.1) then
      write(ilog,241) 'Sigma   '
    else if (sigrule.eq.2) then
      write(ilog,242) 'Sigma   '
    else if (sigrule.eq.3) then
      write(ilog,243) 'Sigma   '
    end if
    if (epsrule.eq.1) then
      write(ilog,241) 'Epsilon '
    else if (epsrule.eq.2) then
      write(ilog,242) 'Epsilon '
    else if (epsrule.eq.3) then
      write(ilog,243) 'Epsilon '
    end if
    write(ilog,*)
    write(ilog,*) 'Exceptions to Combination Rules:'
  end if
!
! contact parameters for pairs of atom types
!
  rewind(unit=iprm)
  do while (.true.)
!
    read (iprm,50,iostat=iomessage) record
    if (iomessage.eq.-1) then
      exit
    else if (iomessage.eq.2) then
      write(ilog,*) 'Fatal. I/O error while reading parameter file (&
 &got: ',paramfile(t1:t2),').'
      call fexit()
    end if
    next = 1
    call extract_str(record,keyword,next)
    call toupper(keyword)
!
    if (keyword(1:8) .eq. 'CONTACT ') then
      ia = 0
      ib = 0
      rd = 0.0d0
      string = record(next:MAXKEYLEN)
      read (string,*,err=120,end=120)  ia,ib,rd
 120      continue
      if ((ia.gt.n_ljtyp).OR.(ib.gt.n_ljtyp).OR.&
 &         (ia.le.0).OR.(ib.le.0)) then
        write(ilog,*) 'Bad parameter file. Contact specification for&
 & illegal atom types (got ',ia,' and ',ib,' in ',paramfile,').'
        call fexit()
      else if (rd.lt.0.0) then
        write(ilog,*) 'Bad parameter file. Illegal contact specifica&
 &tion: ',rd,' for types ',ia,' and ',ib,' (in ',paramfile,').'
        call fexit()
      else if ((vdW_report.EQV..true.).AND.(ia.ne.ib)) then
        write(ilog,244) ia,ib,reduce,rd
      end if
      lj_sig(ia,ib) = red2*rd**2
      lj_sig(ib,ia) = red2*rd**2
!     
! interaction parameters for pairs of atom types
!     
    else if (keyword(1:9) .eq. 'INTERACT ') then
      ia = 0
      ib = 0
      rd = 0.0d0
      string = record(next:MAXKEYLEN)
      read (string,*,err=121,end=121)  ia,ib,rd
 121      continue
      if ((ia.gt.n_ljtyp).OR.(ib.gt.n_ljtyp).OR.&
 &         (ia.le.0).OR.(ib.le.0)) then
        write(ilog,*) 'Bad parameter file. Epsilon specification for&
 & illegal atom types (got ',ia,' and ',ib,' in ',paramfile,').'
        call fexit()
      else if (rd.lt.0.0) then
        write(ilog,*) 'Bad parameter file. Illegal epsilon specifica&
 &tion: ',rd,' for types ',ia,' and ',ib,' (in ',paramfile,').'
        call fexit()
      else if ((vdW_report.EQV..true.).AND.(ia.ne.ib)) then
        write(ilog,245) ia,ib,rd
      end if
      lj_eps(ia,ib) = rd
      lj_eps(ib,ia) = rd
!
!   the same for specific 14-contacts
    else if (keyword(1:11) .eq. 'CONTACT_14 ') then
      ia = 0
      ib = 0
      rd = 0.0d0
      string = record(next:MAXKEYLEN)
      read (string,*,err=123,end=123)  ia,ib,rd
 123      continue
      if ((ia.gt.n_ljtyp).OR.(ib.gt.n_ljtyp).OR.&
 &         (ia.le.0).OR.(ib.le.0)) then
        write(ilog,*) 'Bad parameter file. 14-contact specification for&
 & illegal atom types (got ',ia,' and ',ib,' in ',paramfile,').'
        call fexit()
      else if (rd.lt.0.0) then
        write(ilog,*) 'Bad parameter file. Illegal 14-contact specifica&
 &tion: ',rd,' for types ',ia,' and ',ib,' (in ',paramfile,').'
        call fexit()
      else if ((vdW_report.EQV..true.).AND.(ia.ne.ib)) then
        write(ilog,244) ia,ib,reduce,rd
      end if
      lj_sig_14(ia,ib) = red2*rd**2
      lj_sig_14(ib,ia) = red2*rd**2
!     
! interaction parameters for pairs of atom types
!     
    else if (keyword(1:12) .eq. 'INTERACT_14 ') then
      ia = 0
      ib = 0
      rd = 0.0d0
      string = record(next:MAXKEYLEN)
      read (string,*,err=124,end=124)  ia,ib,rd
 124      continue
      if ((ia.gt.n_ljtyp).OR.(ib.gt.n_ljtyp).OR.&
 &         (ia.le.0).OR.(ib.le.0)) then
        write(ilog,*) 'Bad parameter file. 14-epsilon specification for&
 & illegal atom types (got ',ia,' and ',ib,' in ',paramfile,').'
        call fexit()
      else if (rd.lt.0.0) then
        write(ilog,*) 'Bad parameter file. Illegal 14-epsilon specifica&
 &tion: ',rd,' for types ',ia,' and ',ib,' (in ',paramfile,').'
        call fexit()
      else if ((vdW_report.EQV..true.).AND.(ia.ne.ib)) then
        write(ilog,245) ia,ib,rd
      end if
      lj_eps_14(ia,ib) = rd
      lj_eps_14(ib,ia) = rd

!
! biotype pair bond assignments (checks only)
!     
    else if (keyword(1:17) .eq. 'BONDED_TYPE_BOND ') then
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) ib,ic,ia
      
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
      if ((ia.le.0).OR.(ia.gt.n_bondtyp)) then

        write(ilog,*) 'Bad parameter file. Illegal bond definition f&
 &or BONDED_TYPE_BOND request.'
        call fexit()
      end if
      if ((ib.le.0).OR.(ic.le.0).OR.&
 &      (ib.gt.n_botyp).OR.(ic.gt.n_botyp)) then
        write(ilog,*) 'Bad parameter file. Illegal bonded type reque&
 &st(s) for BONDED_TYPE_BOND request.'
        call fexit()
      end if
      bo_lstsz = bo_lstsz + 1
!
! biotype triple angle assignments (checks only)
!     
    else if (keyword(1:18) .eq. 'BONDED_TYPE_ANGLE ') then
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) ib,ic,id,ia
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
      if ((ia.le.0).OR.(ia.gt.n_angltyp)) then
        write(ilog,*) 'Bad parameter file. Illegal angle definition &
 &for BONDED_TYPE_ANGLE request.'
        call fexit()
      end if
      if ((ib.le.0).OR.(ic.le.0).OR.(id.le.0).OR.&
 & (ib.gt.n_botyp).OR.(ic.gt.n_botyp).OR.(id.gt.n_botyp)) then
        write(ilog,*) 'Bad parameter file. Illegal bonded type reque&
 &st(s) for BONDED_TYPE_ANGLE request.'
        call fexit()
      end if
      ba_lstsz = ba_lstsz + 1
!
! biotype quadruple torsion assignments (checks only)
!     
    else if (keyword(1:20) .eq. 'BONDED_TYPE_TORSION ') then
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) ib,ic,id,ie,ia
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
      if ((ia.le.0).OR.(ia.gt.n_torstyp)) then
        write(ilog,*) 'Bad parameter file. Illegal torsion definitio&
 &n for BONDED_TYPE_TORSION request.'
        call fexit()
      end if
      if ((ib.le.0).OR.(ic.le.0).OR.(id.le.0).OR.(ie.le.0).OR.&
 & (ib.gt.n_botyp).OR.(ic.gt.n_botyp).OR.(id.gt.n_botyp).OR.&
 & (ie.gt.n_botyp)) then
        write(ilog,*) 'Bad parameter file. Illegal bonded type reque&
 &st(s) for BONDED_TYPE_TORSION request.'
        call fexit()
      end if
      di_lstsz = di_lstsz + 1
!
! biotype quadruple improper torsion assignments (checks only)
!     
    else if (keyword(1:20) .eq. 'BONDED_TYPE_IMPTORS ') then
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) ib,ic,id,ie,ia
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
      if ((ia.le.0).OR.(ia.gt.n_torstyp)) then
        write(ilog,*) 'Bad parameter file. Illegal torsion definitio&
 &n for BONDED_TYPE_IMPTORS request.'
        call fexit()
      end if
      if ((ib.le.0).OR.(ic.le.0).OR.(id.le.0).OR.(ie.le.0).OR.&
 & (ib.gt.n_botyp).OR.(ic.gt.n_botyp).OR.(id.gt.n_botyp).OR.&
 & (ie.gt.n_botyp)) then
        write(ilog,*) 'Bad parameter file. Illegal bonded type reque&
 &st(s) for BONDED_TYPE_IMPTORS request.'
        call fexit()
      end if
      impt_lstsz = impt_lstsz + 1
!
! biotype quintuple improper torsion assignments (checks only)
!     
    else if (keyword(1:17) .eq. 'BONDED_TYPE_CMAP ') then
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) ib,ic,id,ie,ig,ia
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
      if ((ia.le.0).OR.(ia.gt.n_cmstyp)) then
        write(ilog,*) 'Bad parameter file. Illegal CMAP definitio&
 &n for BONDED_TYPE_CMAP request.'
        call fexit()
      end if
      if ((ib.le.0).OR.(ic.le.0).OR.(id.le.0).OR.(ie.le.0).OR.(ig.le.0).OR.&
 & (ib.gt.n_botyp).OR.(ic.gt.n_botyp).OR.(id.gt.n_botyp).OR.&
 & (ie.gt.n_botyp).OR.(ig.gt.n_botyp)) then
        write(ilog,*) 'Bad parameter file. Illegal bonded type reque&
 &st(s) for BONDED_TYPE_CMAP request.'
        call fexit()
      end if
      cm_lstsz = cm_lstsz + 1

!     
! hbond parameters for pairs of atom types
!     
!    else if (keyword(1:7) .eq. 'HBPAIR ') then
!       ia = 0
!       ib = 0
!       rd = 0.0d0
!       string = record(next:MAXKEYLEN)
!       read (string,*,err=130,end=130)  ia,ib,rd
! 130       continue
!       hbdist(ia,ib) = rd**2
!       hbdist(ib,ia) = rd**2
!
    end if
!
  end do


! #tyler pka: lets fill in the extra LJ tables (move to after these tables are defined
! n_ljtyp_temp+(1:3) is the used hydrogen,oxygen, oxygen on a GLX (from 4:6)


!  if (do_pka_2.eqv..true.) then ! Martin : Modified (simpler, easier to understand, only one loop)
!  do i=1,size(transform_table)
!    lj_eps_limits(bio_ljtyp(i),:,2)=lj_eps(bio_ljtyp(transform_table(i)),:)
!    lj_eps_limits(bio_ljtyp(i),:,1)=lj_eps(bio_ljtyp(i),:)
!    lj_eps_14_limits(bio_ljtyp(i),:,2)=lj_eps_14(bio_ljtyp(transform_table(i)),:)
!    lj_eps_14_limits(bio_ljtyp(i),:,1)=lj_eps_14(bio_ljtyp(i),:)
!    lj_sig_limits(bio_ljtyp(i),:,2)=lj_sig(bio_ljtyp(transform_table(i)),:)
!    lj_sig_limits(bio_ljtyp(i),:,1)=lj_sig(bio_ljtyp(i),:)
!    lj_sig_14_limits(bio_ljtyp(i),:,2)=lj_sig(bio_ljtyp(transform_table(i)),:)
!    lj_sig_14_limits(bio_ljtyp(i),:,1)=lj_sig(bio_ljtyp(i),:)
!  end do
!end if

! Martin : Old version for the sake of record keeping 

  
!if (do_pka_2.eqv..true.) then
!  do i=1,size(lj_limits_table,1)
!    lj_eps(n_ljtyp_temp+i,:)=lj_eps(lj_limits_table(i,2),:)
!    lj_eps(:,n_ljtyp_temp+i)=lj_eps(:,lj_limits_table(i,2))
!    lj_eps_14(n_ljtyp_temp+i,:)=lj_eps_14(lj_limits_table(i,2),:)
!    lj_eps_14(:,n_ljtyp_temp+i)=lj_eps_14(:,lj_limits_table(i,2))
!    lj_sig(n_ljtyp_temp+i,:)=lj_sig(lj_limits_table(i,2),:)
!    lj_sig(:,n_ljtyp_temp+i)=lj_sig(:,lj_limits_table(i,2))
!    lj_sig_14(n_ljtyp_temp+i,:)=lj_sig_14(lj_limits_table(i,2),:)
!    lj_sig_14(:,n_ljtyp_temp+i)=lj_sig_14(:,lj_limits_table(i,2))
!    do j=1,size(lj_limits_table,1)
!      lj_eps(n_ljtyp_temp+i,n_ljtyp_temp+j)=lj_eps(lj_limits_table(i,2),lj_limits_table(j,2))
!      lj_eps_14(n_ljtyp_temp+i,n_ljtyp_temp+j)=lj_eps_14(lj_limits_table(i,2),lj_limits_table(j,2))
!      lj_sig(n_ljtyp_temp+i,n_ljtyp_temp+j)=lj_sig(lj_limits_table(i,2),lj_limits_table(j,2))
!      lj_sig_14(n_ljtyp_temp+i,n_ljtyp_temp+j)=lj_sig_14(lj_limits_table(i,2),lj_limits_table(j,2))
!    end do
!  end do
!end if


!
! finally, mappings from bonded types to bond/angle/torsion types
!
  allocate(bo_lst(bo_lstsz,3))
  allocate(ba_lst(ba_lstsz,4))
  allocate(di_lst(di_lstsz,5))
  allocate(impt_lst(impt_lstsz,5))
  allocate(cm_lst(cm_lstsz,6))
  kbo = 0
  kba = 0
  kdi = 0
  kdi2 = 0
  kcm = 0
  rewind(unit=iprm)
  do while (.true.)
!
    read (iprm,50,iostat=iomessage) record
    if (iomessage.eq.-1) then
      exit
    else if (iomessage.eq.2) then
      write(ilog,*) 'Fatal. I/O error while reading parameter file (&
 &got: ',paramfile(t1:t2),').'
      call fexit()
    end if
    next = 1
    call extract_str(record,keyword,next)
    call toupper(keyword)
!
! biotype pair bond assignments (checks only)
!     
    if (keyword(1:17) .eq. 'BONDED_TYPE_BOND ') then
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) ib,ic,ia
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
      kbo = kbo + 1
      bo_lst(kbo,1) = ib
      bo_lst(kbo,2) = ic
      bo_lst(kbo,3) = ia
!
! biotype triple angle assignments
!     
    else if (keyword(1:18) .eq. 'BONDED_TYPE_ANGLE ') then
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) ib,ic,id,ia
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
      kba = kba + 1 ! Martin : here is what I want to use
      ba_lst(kba,1) = ib
      ba_lst(kba,2) = ic
      ba_lst(kba,3) = id
      ba_lst(kba,4) = ia

! biotype quadruple torsion assignments
     
    else if (keyword(1:20) .eq. 'BONDED_TYPE_TORSION ') then
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) ib,ic,id,ie,ia
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
      kdi = kdi + 1
      di_lst(kdi,1) = ib
      di_lst(kdi,2) = ic
      di_lst(kdi,3) = id
      di_lst(kdi,4) = ie
      di_lst(kdi,5) = ia
!
! biotype quadruple improper dihedral assignments
!     
    else if (keyword(1:20) .eq. 'BONDED_TYPE_IMPTORS ') then
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) ib,ic,id,ie,ia
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
      kdi2 = kdi2 + 1
      impt_lst(kdi2,1) = ib
      impt_lst(kdi2,2) = ic
      impt_lst(kdi2,3) = id
      impt_lst(kdi2,4) = ie
      impt_lst(kdi2,5) = ia
!
! biotype quintuple CMAP assignments
!     
    else if (keyword(1:17) .eq. 'BONDED_TYPE_CMAP ') then
      string = record(next:MAXKEYLEN)
      read (string,*,iostat=iomessage) ib,ic,id,ie,ig,ia
      if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading parameter file&
 & (got: ',paramfile(t1:t2),').'
        call fexit()
      end if
      kcm = kcm + 1
      cm_lst(kcm,1) = ib
      cm_lst(kcm,2) = ic
      cm_lst(kcm,3) = id
      cm_lst(kcm,4) = ie
      cm_lst(kcm,5) = ig
      cm_lst(kcm,6) = ia
!
    end if
  end do
!
! we must check those lists for double entries
  call check_bndprms()
!
  if (vdW_report.EQV..true.) then
    write(ilog,*)
    write(ilog,*) 'Self-Terms (From Which to Combine):'
  end if
  do ia=1,n_ljtyp
    if ((lj_sig(ia,ia).lt.0.0).OR.(lj_eps(ia,ia).lt.0.0)) then
      write(ilog,*) 'Fatal. Missing parameters for type ',ia,' (in ',paramfile,').'
      call fexit()
    end if
    if (vdW_report.EQV..true.) then
      write(ilog,246) ia,reduce,sqrt(lj_sig(ia,ia))/reduce
      write(ilog,247) lj_eps(ia,ia)
    end if
  end do
  if (vdW_report.EQV..true.) then
    write(ilog,*)
    write(ilog,*) 'Specific 14-Self-Terms (Override Only):'
  end if
  do ia=1,n_ljtyp
    if ((lj_sig_14(ia,ia).lt.0.0).AND.(lj_eps_14(ia,ia).lt.0.0)) then
!     do nothing (just yet)
    else if (((lj_sig_14(ia,ia).lt.0.0).AND.(lj_eps_14(ia,ia).ge.0.0)).OR.&
            &((lj_sig_14(ia,ia).ge.0.0).AND.(lj_eps_14(ia,ia).lt.0.0))) then
      write(ilog,*) 'Fatal. Please provide consistent sets of overrides for 14-LJ interaction parameters.'
      write(ilog,*) 'Either epsilon (INTERACT_14) or sigma (CONTACT_14) is missing for LJ type ',ia,'.'
      call fexit()
    else
      if (vdW_report.EQV..true.) then
        write(ilog,246) ia,reduce,sqrt(lj_sig_14(ia,ia))/reduce
        write(ilog,247) lj_eps_14(ia,ia)          
      end if
    end if
  end do
  do ia=1,n_ljtyp
    do ib=ia+1,n_ljtyp
      if (((lj_sig_14(ia,ib).lt.0.0).AND.(lj_eps_14(ia,ib).ge.0.0)).OR.&
         &((lj_sig_14(ia,ib).ge.0.0).AND.(lj_eps_14(ia,ib).lt.0.0))) then
        write(ilog,*) 'Fatal. Please provide consistent sets of overrides for 14-LJ interaction parameters.'
        write(ilog,*) 'Either epsilon (INTERACT_14) or sigma (CONTACT_14) is missing for LJ types ',ia,' and ',ib,'.'
        call fexit()
      end if
    end do
  end do
!
  if ((do_pka_2.eqv..false.).and.(do_hsq.eqv..false.)) then
    call clean_radik()
    call clean_epsik()
    ! This will be done later, in the case of pka calculatons should not interfer
  end if
 
! set the default parametric radii (never use 14-overrides here)
  do i=1,n_ljtyp
    lj_rad(i) = 0.5*sqrt(lj_sig(i,i)) ! Martin : I really do not get why the 0.5
  end do
 
!
  if (vdW_report.EQV..true.) then
    write(ilog,*)
    write(ilog,*) 'Radius Overrides (for Solvent-Accessibility):'
  end if
!
! and finally overwrites for atomic radii
!
  rewind(unit=iprm)
  do while (.true.)
!
    read (iprm,50,iostat=iomessage) record
    if (iomessage.eq.-1) then
      exit
    else if (iomessage.eq.2) then
      write(ilog,*) 'Fatal. I/O error while reading parameter file (&
 &got: ',paramfile(t1:t2),').'
      call fexit()
    end if
    next = 1
    call extract_str(record,keyword,next)
    call toupper(keyword)
!
    if (keyword(1:7) .eq. 'RADIUS ') then
       ia = 0
       rd = 0.0d0
       string = record(next:MAXKEYLEN)
       read (string,*,err=122,end=122)  ia,rd
 122       continue
       if ((ia.gt.n_ljtyp).OR.(ia.le.0)) then
         write(ilog,*) 'Bad parameter file. Radius specification for&
 & illegal atom type (got ',ia,' in ',paramfile,').'
         call fexit()
       else if (rd.le.0.0) then
         write(ilog,*) 'Bad parameter file. Illegal radius specifica&
 &tion: ',rd,' for type ',ia,' (in ',paramfile,').'
         call fexit()
       else if (vdW_report.EQV..true.) then
         write(ilog,248) ia,rd
       end if
       lj_rad(ia) = rd
     end if
!
  end do
!
  if (vdW_report.EQV..true.) then
    write(ilog,*)
  end if
!
  call fos_summary()


!
end
!
!
!
! ###################################################
! ##                                               ##
! ##  subroutine seteps -- sets pairwise epsilon's ##
! ##                                               ##
! ###################################################
!
!
! "seteps" calculates the pairwise epsilon parameters for
! each possible inter-atomic contact for use in the inverse
! power potential function
!
! Literature reference: "General Chemistry", by Linus Pauling, 
! Chapter 11, pg., 396. 
!
!  subroutine epsik(ia,ib)
!  implicit none
!  include 'sizes.i"
!  include 'inter.i"
!  include 'kcont.i"
!  integer ia,ib
!  RTYPE alpha(17),baseps(17),c2,const,radik,r6
!c      data alpha /0.90,1.03,1.34,0.4,0.85,0.69,3.0,0.97,0.4,0.9,0.9,
!c     &0.18,3.69/
!  data alpha /0.90,1.03,1.34,0.4,0.85,0.69,3.0,0.97,0.4,0.9,0.9,
! &0.18,3.69,0.69,0.4,0.69,0.4/ 
!  data baseps /0.1,0.1,0.1,0.025,0.1,0.1,0.3,0.1,0.025,0.1,0.1,0.1,
! &0.3,0.155406042,0.0,0.152072595,0.0/
!  parameter (const=334.60803d0)
!     
!  if (eps_type.eq.1) then
!     lj_eps(ia,ib) = 1.0d0
!     lj_eps(ib,ia) = 1.0d0
!  else if (eps_type.eq.2) then
!     lj_eps(ia,ib) = baseps(ia)*baseps(ib)
!     if (lj_eps(ia,ib).gt.0.0) then
!       lj_eps(ia,ib) = sqrt(lj_eps(ia,ib))
!     else
!       lj_eps(ia,ib) = 0.0
!     end if
!     lj_eps(ib,ia) = lj_eps(ia,ib)
!c         write(*,*)ia,' ',ib,':',lj_eps(ia,ib),' ',sqrt(lj_sig(ia,ib))
!  else
!     c2 = const*alpha(ia)*alpha(ib)
!     radik = lj_sig(ia,ib)
!     if (radik.le.0.0) then
!       lj_eps(ia,ib) = 0.0
!       lj_eps(ib,ia) = 0.0
!     else
!       r6 = radik**3
!       lj_eps(ia,ib) = 0.25d0*c2/r6
!       lj_eps(ib,ia) = lj_eps(ia,ib)
!     end if
!c         write(*,*)ia,' ',ib,':',lj_eps(ia,ib),' ',sqrt(lj_sig(ia,ib))
!  end if
!  return 
!  end
!c
!----------------------------------------------------------------------------
!
! this subroutine will fill in sigma-values according to the spec.d mean
! derived from the self-terms
! it will, however, not overwrite any explicitly read-in contact distances
! for specific 14 overrides (lj_sig_14) it will use the same combination rules and settings and
! assume the value for lj_sig if nothing is specified
!
subroutine clean_radik()
!
  use params
!
  implicit none
!
  integer ia,ib
!
  do ia=1,n_ljtyp
    do ib=ia+1,n_ljtyp
      if (lj_sig(ia,ib).lt.0.0) then
        if ((lj_sig(ia,ia).gt.0.0).AND.(lj_sig(ib,ib).gt.0.0))then
          if (sigrule.eq.1) then
            lj_sig(ia,ib) = 0.5*(sqrt(lj_sig(ia,ia))+sqrt(lj_sig(ib,ib)))
            lj_sig(ia,ib) = lj_sig(ia,ib)*lj_sig(ia,ib)
          else if (sigrule.eq.2) then
            lj_sig(ia,ib) = sqrt(lj_sig(ia,ia)*lj_sig(ib,ib))
          else if (sigrule.eq.3) then
            lj_sig(ia,ib) = 2.0/( (1.0/sqrt(lj_sig(ia,ia))) + (1.0/sqrt(lj_sig(ib,ib))) )
            lj_sig(ia,ib) = lj_sig(ia,ib)*lj_sig(ia,ib)
          end if
        else if ((lj_sig(ia,ia).gt.0.0).OR.(lj_sig(ib,ib).gt.0.0)) then
          if (sigrule.eq.1) then
            lj_sig(ia,ib) = 0.5*(sqrt(lj_sig(ia,ia))+sqrt(lj_sig(ib,ib)))
            lj_sig(ia,ib) = lj_sig(ia,ib)*lj_sig(ia,ib)
          else
            lj_sig(ia,ib) = 0.0
          end if
        else
          lj_sig(ia,ib) = 0.0
        end if
        lj_sig(ib,ia) = lj_sig(ia,ib)
      end if
    end do
  end do
! do the same for 14-specifics (but watch for preserving defaults)
  do ia=1,n_ljtyp
    do ib=ia+1,n_ljtyp
      if (lj_sig_14(ia,ib).lt.0.0) then
!       this case means both of them are set to non-zero values
        if ((lj_sig_14(ia,ia).gt.0.0).AND.(lj_sig_14(ib,ib).gt.0.0)) then
          if (sigrule.eq.1) then
            lj_sig_14(ia,ib) = 0.5*(sqrt(lj_sig_14(ia,ia))+sqrt(lj_sig_14(ib,ib)))
            lj_sig_14(ia,ib) = lj_sig_14(ia,ib)*lj_sig_14(ia,ib)
          else if (sigrule.eq.2) then
            lj_sig_14(ia,ib) = sqrt(lj_sig_14(ia,ia)*lj_sig_14(ib,ib))
          else if (sigrule.eq.3) then
            lj_sig_14(ia,ib) = 2.0/( (1.0/sqrt(lj_sig_14(ia,ia))) + (1.0/sqrt(lj_sig_14(ib,ib))) )
            lj_sig_14(ia,ib) = lj_sig_14(ia,ib)*lj_sig_14(ia,ib)
          end if
!       this case means only one ia is set to a non-zero value
        else if ((lj_sig_14(ia,ia).gt.0.0).AND.(lj_sig_14(ib,ib).lt.0.0)) then
          if (sigrule.eq.1) then
            lj_sig_14(ia,ib) = 0.5*(sqrt(lj_sig_14(ia,ia))+sqrt(lj_sig(ib,ib)))
            lj_sig_14(ia,ib) = lj_sig_14(ia,ib)*lj_sig_14(ia,ib)
          else if (sigrule.eq.2) then
            lj_sig_14(ia,ib) = sqrt(lj_sig_14(ia,ia)*lj_sig(ib,ib))
          else if (sigrule.eq.3) then
            lj_sig_14(ia,ib) = 2.0/( (1.0/sqrt(lj_sig_14(ia,ia))) + (1.0/sqrt(lj_sig(ib,ib))) )
            lj_sig_14(ia,ib) = lj_sig_14(ia,ib)*lj_sig_14(ia,ib)
          end if
!       this case means only one ib is set to a non-zero value
        else if ((lj_sig_14(ia,ia).lt.0.0).AND.(lj_sig_14(ib,ib).gt.0.0)) then
          if (sigrule.eq.1) then
            lj_sig_14(ia,ib) = 0.5*(sqrt(lj_sig(ia,ia))+sqrt(lj_sig_14(ib,ib)))
            lj_sig_14(ia,ib) = lj_sig_14(ia,ib)*lj_sig_14(ia,ib)
          else if (sigrule.eq.2) then
            lj_sig_14(ia,ib) = sqrt(lj_sig(ia,ia)*lj_sig_14(ib,ib))
          else if (sigrule.eq.3) then
            lj_sig_14(ia,ib) = 2.0/( (1.0/sqrt(lj_sig(ia,ia))) + (1.0/sqrt(lj_sig_14(ib,ib))) )
            lj_sig_14(ia,ib) = lj_sig_14(ia,ib)*lj_sig_14(ia,ib)
          end if
!       this case means both are explicitly set to zero
        else if ((lj_sig_14(ia,ia).eq.0.0).AND.(lj_sig_14(ib,ib).eq.0.0)) then
          lj_sig_14(ia,ib) = 0.0
!       this case means both are explicitly set but at least one of them is zero
        else if ((lj_sig_14(ia,ia).ge.0.0).AND.(lj_sig_14(ib,ib).ge.0.0)) then
          if (sigrule.eq.1) then
            lj_sig_14(ia,ib) = 0.5*(sqrt(lj_sig_14(ia,ia))+sqrt(lj_sig_14(ib,ib)))
            lj_sig_14(ia,ib) = lj_sig_14(ia,ib)*lj_sig_14(ia,ib)
          else
            lj_sig_14(ia,ib) = 0.0
          end if
!       this case means ia is set to zero (only arithmetic meaningful) with ib default
        else if ((lj_sig_14(ia,ia).eq.0.0).AND.(lj_sig_14(ib,ib).lt.0.0)) then
          if (sigrule.eq.1) then
            lj_sig_14(ia,ib) = 0.5*(sqrt(lj_sig_14(ia,ia))+sqrt(lj_sig(ib,ib)))
            lj_sig_14(ia,ib) = lj_sig_14(ia,ib)*lj_sig_14(ia,ib)
          else
            lj_sig_14(ia,ib) = 0.0 
          end if
!       this case means ib is set to zero (only arithmetic meaningful) with ia default
        else if ((lj_sig_14(ia,ia).lt.0.0).AND.(lj_sig_14(ib,ib).eq.0.0)) then
          if (sigrule.eq.1) then
            lj_sig_14(ia,ib) = 0.5*(sqrt(lj_sig(ia,ia))+sqrt(lj_sig_14(ib,ib)))
            lj_sig_14(ia,ib) = lj_sig_14(ia,ib)*lj_sig_14(ia,ib)
          else
            lj_sig_14(ia,ib) = 0.0 
          end if
        else
          lj_sig_14(ia,ib) = lj_sig(ia,ib)
        end if
        lj_sig_14(ib,ia) = lj_sig_14(ia,ib)
      end if
!      write(*,*) ia,ib,sqrt(lj_sig_14(ia,ib)),sqrt(lj_sig(ia,ib))
    end do
  end do
! now finalize by setting the missing self-terms to the default
  do ia=1,n_ljtyp
    if (lj_sig_14(ia,ia).lt.0.0) then
!     copy over value from lj_sig
      lj_sig_14(ia,ia) = lj_sig(ia,ia)
    end if
  end do 

!
end
!
!----------------------------------------------------------------------------
!
! this subroutine will fill in epsilon-values according to the spec.d mean
! derived from the self-terms
! it will, however, not overwrite any explicitly read-in interaction params
! in addition, net zero terms will be arranged such that the zero comes exclusively from *_eps*
!
subroutine clean_epsik()
!
  use params
  use iounit
  use energies
!
  implicit none
!
  integer ia,ib
!
 77 format('Warning. Pairwise LJ interaction parameters for atom types ',i6,' and ',i6,' are inconsistent with &
 &finite strength and zero contact distance. Buffering contact distance and setting strength to zero instead.')
 78 format('Warning. Pairwise LJ interaction parameters for atom types ',i6,' and ',i6,' are zero for &
 &both strength and contact distance. Buffering contact distance slightly.')
 79 format('Warning. Pairwise 14 overrides for LJ interaction parameters for atom types ',i6,' and ',i6,' are inconsistent with &
 &finite strength and zero contact distance. Buffering contact distance and setting strength to zero instead.')
 80 format('Warning. Pairwise 14 overrides for LJ interaction parameters for atom types ',i6,' and ',i6,' are zero for &
 &both strength and contact distance. Buffering contact distance slightly.')
!
  do ia=1,n_ljtyp
    do ib=ia+1,n_ljtyp
      if (lj_eps(ia,ib).lt.0.0) then ! Martin : why lt ? 
          
        if((lj_eps(ia,ia).gt.0.0).AND.(lj_eps(ib,ib).gt.0.0))then
          if (epsrule.eq.1) then
            lj_eps(ia,ib) = 0.5*(lj_eps(ia,ia) + lj_eps(ib,ib))
          else if (epsrule.eq.2) then
            lj_eps(ia,ib) = sqrt(lj_eps(ia,ia)*lj_eps(ib,ib))
          else if (epsrule.eq.3) then
            lj_eps(ia,ib) = 2.0/( (1.0/lj_eps(ia,ia)) +&
 &              (1.0/lj_eps(ib,ib)) )
          end if
        else if ((lj_eps(ia,ia).gt.0.0).OR.(lj_eps(ib,ib).gt.0.0))&
 &  then
          if (epsrule.eq.1) then
            lj_eps(ia,ib) = 0.5*(lj_eps(ia,ia) + lj_eps(ib,ib))
          else
            lj_eps(ia,ib) = 0.0
          end if
        else
          lj_eps(ia,ib) = 0.0
        end if
        lj_eps(ib,ia) = lj_eps(ia,ib)
!        write(*,*) ia,ib,lj_eps(ia,ib)
      end if
    end do
  end do
! do the same for 14-specifics
  do ia=1,n_ljtyp
    do ib=ia+1,n_ljtyp
      if (lj_eps_14(ia,ib).lt.0.0) then
        if((lj_eps_14(ia,ia).gt.0.0).AND.(lj_eps_14(ib,ib).gt.0.0))then
          if (epsrule.eq.1) then
            lj_eps_14(ia,ib) = 0.5*(lj_eps_14(ia,ia) + lj_eps_14(ib,ib))
          else if (epsrule.eq.2) then
            lj_eps_14(ia,ib) = sqrt(lj_eps_14(ia,ia)*lj_eps_14(ib,ib))
          else if (epsrule.eq.3) then
            lj_eps_14(ia,ib) = 2.0/( (1.0/lj_eps_14(ia,ia)) + (1.0/lj_eps_14(ib,ib)) )
          end if
        else if ((lj_eps_14(ia,ia).gt.0.0).AND.(lj_eps_14(ib,ib).lt.0.0)) then
          if (epsrule.eq.1) then
            lj_eps_14(ia,ib) = 0.5*(lj_eps_14(ia,ia) + lj_eps(ib,ib))
          else if (epsrule.eq.2) then
            lj_eps_14(ia,ib) = sqrt(lj_eps_14(ia,ia)*lj_eps(ib,ib))
          else if (epsrule.eq.3) then
            lj_eps_14(ia,ib) = 2.0/( (1.0/lj_eps_14(ia,ia)) + (1.0/lj_eps(ib,ib)) )
          end if
        else if ((lj_eps_14(ia,ia).lt.0.0).AND.(lj_eps_14(ib,ib).gt.0.0)) then
          if (epsrule.eq.1) then
            lj_eps_14(ia,ib) = 0.5*(lj_eps(ia,ia) + lj_eps_14(ib,ib))
          else if (epsrule.eq.2) then
            lj_eps_14(ia,ib) = sqrt(lj_eps(ia,ia)*lj_eps_14(ib,ib))
          else if (epsrule.eq.3) then
            lj_eps_14(ia,ib) = 2.0/( (1.0/lj_eps(ia,ia)) + (1.0/lj_eps_14(ib,ib)) )
          end if
        else if ((lj_eps_14(ia,ia).eq.0.0).AND.(lj_eps_14(ib,ib).eq.0.0)) then
          lj_eps_14(ia,ib) = 0.0
        else if ((lj_eps_14(ia,ia).ge.0.0).AND.(lj_eps_14(ib,ib).ge.0.0)) then
          if (epsrule.eq.1) then
            lj_eps_14(ia,ib) = 0.5*(lj_eps_14(ia,ia) + lj_eps_14(ib,ib))
          else
            lj_eps_14(ia,ib) = 0.0
          end if
        else if ((lj_eps_14(ia,ia).eq.0.0).AND.(lj_eps_14(ib,ib).lt.0.0)) then
          if (epsrule.eq.1) then
            lj_eps_14(ia,ib) = 0.5*(lj_eps_14(ia,ia) + lj_eps(ib,ib))
          else
            lj_eps_14(ia,ib) = 0.0
          end if
        else if ((lj_eps_14(ia,ia).lt.0.0).AND.(lj_eps_14(ib,ib).eq.0.0)) then
          if (epsrule.eq.1) then
            lj_eps_14(ia,ib) = 0.5*(lj_eps(ia,ia) + lj_eps_14(ib,ib))
          else
            lj_eps_14(ia,ib) = 0.0
          end if
        else
          lj_eps_14(ia,ib) = lj_eps(ia,ib)
        end if
        lj_eps_14(ib,ia) = lj_eps_14(ia,ib)
      end if
!      write(*,*) ia,ib,lj_eps_14(ia,ib),lj_eps(ia,ib)
    end do
  end do
! now finalize by setting the missing self-terms to the default
  do ia=1,n_ljtyp
     if (lj_eps_14(ia,ia).lt.0.0) then
!     copy over value from lj_eps
      lj_eps_14(ia,ia) = lj_eps(ia,ia)
    end if
  end do
! corrections for numerically problematic cases (when using the inverse of lj_sig)
!  if ((do_hsq.eqv..false.).and.(do_pka.eqv..false.)) then  
      do ia=1,n_ljtyp
        do ib=ia,n_ljtyp
          if (lj_sig(ia,ib).le.0.0) then
            if (lj_eps(ia,ib).le.0.0) then
              write(ilog,78) ia,ib
              lj_sig(ia,ib) = 0.01
              lj_sig(ib,ia) = 0.01 
            else
              write(ilog,77) ia,ib
              lj_sig(ia,ib) = 0.01
              lj_sig(ib,ia) = 0.01
            end if
          end if
          if (lj_sig_14(ia,ib).le.0.0) then
            if (lj_eps_14(ia,ib).le.0.0) then
              write(ilog,80) ia,ib
              lj_sig_14(ia,ib) = 0.01
              lj_sig_14(ib,ia) = 0.01
            else
              write(ilog,79) ia,ib
              lj_sig_14(ia,ib) = 0.01
              lj_sig_14(ib,ia) = 0.01
            end if
          end if
        end do
      end do
!  end if 
  


!
!
end
!
!--------------------------------------------------------------------------------
!
subroutine read_fos(string)
!
  use iounit
  use params
  use aminos
  use keys
  use fos
!
  implicit none
!
  integer i,j,seqt,length,next,t1,t2,kk,nvalsr
  character(MAXKWLEN) keyword
  character(MAXKEYLEN) string,string2
  character(3) resname
  RTYPE readv(3)
!
  call strlims(paramfile,t1,t2)
!
  next = 1
  call extract_str(string,keyword,next)
  call toupper(keyword)
  string2 = string(next:MAXKEYLEN)
!
  call strlims(keyword,t1,t2)
  length = t2-t1+1
!
  if ((length.eq.3).OR.(length.eq.2)) then
    seqt = 0
    if (length.eq.2) then
      keyword(3:3) = ' '
    end if
    i = 0
    do j = 1,MAXAMINO
      if(keyword(1:3) .eq. amino(j)) then
        i = i + 1
        seqt = j
      end if
    end do
!
    readv(:) = 0.0
    kk = 3
    next = 1
    call get_reals_from_str(string2,kk,readv,next,nvalsr)
    if (nvalsr.le.0) then
      write(ilog,*) 'Fatal. I/O error while reading FOS data in parameter file (got: ',paramfile(t1:t2),').'
      call fexit()
    end if
    if (nvalsr.eq.1) then
      readv(2) = readv(1)
      readv(3) = 0.0
    else if (nvalsr.eq.2) then
      readv(3) = 0.0
    end if
!
!   now: single-residue molecules and sidechain analogs first
    if (seqt.gt.0) then
      resname = amino(seqt)
    else
      resname = '   '
    end if
!   the peptide sidechain analogs (note no caps and no glycine) 
    if (resname.eq.'ALA') then
      FOS_ALA(1:3) = readv(1:3)
    else if (resname.eq.'VAL') then
      FOS_VAL(1:3) = readv(1:3)
    else if (resname.eq.'LEU') then
      FOS_LEU(1:3) = readv(1:3)
    else if (resname.eq.'ILE') then
      FOS_ILE(1:3) = readv(1:3)
    else if (resname.eq.'SER') then
      FOS_SER(1:3) = readv(1:3)
    else if (resname.eq.'THR') then
      FOS_THR(1:3) = readv(1:3)
    else if (resname.eq.'CYS') then
      FOS_CYS(1:3) = readv(1:3)
    else if (resname.eq.'PRO') then
      FOS_PRO(1:3) = readv(1:3)
    else if (resname.eq.'PHE') then
      FOS_PHE(1:3) = readv(1:3)
    else if (resname.eq.'TYR') then
      FOS_TYR(1:3) = readv(1:3)
    else if (resname.eq.'TRP') then
      FOS_TRP(1:3) = readv(1:3)
    else if (resname.eq.'HID') then
      FOS_HID(1:3) = readv(1:3)
    else if (resname.eq.'HIE') then
      FOS_HIE(1:3) = readv(1:3)
    else if (resname.eq.'HIP') then
      FOS_HIP(1:3) = readv(1:3)
    else if (resname.eq.'ASP') then
      FOS_ASP(1:3) = readv(1:3)
    else if (resname.eq.'ASN') then
      FOS_ASN(1:3) = readv(1:3)
    else if (resname.eq.'GLU') then
      FOS_GLU(1:3) = readv(1:3)
    else if (resname.eq.'GLN') then
      FOS_GLN(1:3) = readv(1:3)
    else if (resname.eq.'MET') then
      FOS_MET(1:3) = readv(1:3)
    else if (resname.eq.'LYS') then
      FOS_LYS(1:3) = readv(1:3)
    else if (resname.eq.'ARG') then
      FOS_ARG(1:3) = readv(1:3)
    else if (resname.eq.'ORN') then
      FOS_ORN(1:3) = readv(1:3)
    else if (resname.eq.'DAB') then
      FOS_DAB(1:3) = readv(1:3)
    else if (resname.eq.'AIB') then
      FOS_AIB(1:3) = readv(1:3)
    else if (resname.eq.'PCA') then
      FOS_PCA(1:3) = readv(1:3)
    else if (resname.eq.'GAM') then
      FOS_GAM(1:3) = readv(1:3)
    else if (resname.eq.'HYP') then
      FOS_HYP(1:3) = readv(1:3)
    else if (resname.eq.'ABA') then
      FOS_ABA(1:3) = readv(1:3)
    else if (resname.eq.'NVA') then
      FOS_NVA(1:3) = readv(1:3)
    else if (resname.eq.'NLE') then
      FOS_NLE(1:3) = readv(1:3)
    else if (resname.eq.'ASH') then
      FOS_ASH(1:3) = readv(1:3)
    else if (resname.eq.'GLH') then
      FOS_GLH(1:3) = readv(1:3)
    else if (resname.eq.'GLX') then
      FOS_GLX(1:3) = readv(1:3)
    else if (resname.eq.'TYO') then
      FOS_TYO(1:3) = readv(1:3)
    else if (resname.eq.'CYX') then
      FOS_CYX(1:3) = readv(1:3)
    else if (resname.eq.'LYD') then
      FOS_LYD(1:3) = readv(1:3)
    else if (resname.eq.'LYX') then
      FOS_LYX(1:3) = readv(1:3)
    else if (resname.eq.'S1P') then
      FOS_S1P(1:3) = readv(1:3)
    else if (resname.eq.'T1P') then
      FOS_T1P(1:3) = readv(1:3)
    else if (resname.eq.'Y1P') then
      FOS_Y1P(1:3) = readv(1:3)
    else if (resname.eq.'S2P') then
      FOS_S2P(1:3) = readv(1:3)
    else if (resname.eq.'T2P') then
      FOS_T2P(1:3) = readv(1:3)
    else if (resname.eq.'Y2P') then
      FOS_Y2P(1:3) = readv(1:3)
      !   nucleotide sidechain analogs
    else if (keyword(1:3).eq.'XXC') then
      FOS_XPC(1:3) = readv(1:3)
    else if (keyword(1:3).eq.'XXU') then
      FOS_XPU(1:3) = readv(1:3)
    else if (keyword(1:3).eq.'XXT') then
      FOS_XPT(1:3) = readv(1:3)
    else if (keyword(1:3).eq.'XXA') then
      FOS_XPA(1:3) = readv(1:3)
    else if (keyword(1:3).eq.'XXG') then
      FOS_XPG(1:3) = readv(1:3)
!   ions and small molecules (there should be some correspondence as
!   these are related!)
    else if (resname.eq.'NA+') then
      FOS_NA(1:3) = readv(1:3)
    else if (resname.eq.'CL-') then
      FOS_CL(1:3) = readv(1:3)
    else if (resname.eq.'URE') then
      FOS_URE(1:3) = readv(1:3)
    else if (resname.eq.'SPC') then
      FOS_SPC(1:3) = readv(1:3)
    else if (resname.eq.'T3P') then
      FOS_T3P(1:3) = readv(1:3)
    else if (resname.eq.'NMF') then
      FOS_NMF(1:3) = readv(1:3)
    else if (resname.eq.'NMA') then
      FOS_NMA(1:3) = readv(1:3)
    else if (resname.eq.'ACA') then
      FOS_ACA(1:3) = readv(1:3)
    else if (resname.eq.'PPA') then
      FOS_PPA(1:3) = readv(1:3)
    else if (resname.eq.'T4P') then
      FOS_T4P(1:3) = readv(1:3)
    else if (resname.eq.'T4E') then
      FOS_T4E(1:3) = readv(1:3)
    else if (resname.eq.'T5P') then
      FOS_T5P(1:3) = readv(1:3)
    else if (resname.eq.'CH4') then
      FOS_CH4(1:3) = readv(1:3)
    else if (resname.eq.'FOA') then
      FOS_FOA(1:3) = readv(1:3)
    else if (resname.eq.'DMA') then
      FOS_DMA(1:3) = readv(1:3)
    else if (resname.eq.'MOH') then
      FOS_MOH(1:3) = readv(1:3)
    else if (resname.eq.'PCR') then
      FOS_PCR(1:3) = readv(1:3)
    else if (resname.eq.'MSH') then
      FOS_MSH(1:3) = readv(1:3)
    else if (resname.eq.'EOH') then
      FOS_EOH(1:3) = readv(1:3)
    else if (resname.eq.'PRP') then
      FOS_PRP(1:3) = readv(1:3)
    else if (resname.eq.'IBU') then
      FOS_IBU(1:3) = readv(1:3)
    else if (resname.eq.'NBU') then
      FOS_NBU(1:3) = readv(1:3)
    else if (resname.eq.'EMT') then
      FOS_EMT(1:3) = readv(1:3)
    else if (resname.eq.'TOL') then
      FOS_TOL(1:3) = readv(1:3)
    else if (resname.eq.'IMD') then
      FOS_IMD(1:3) = readv(1:3)
    else if (resname.eq.'IME') then
      FOS_IME(1:3) = readv(1:3)
    else if (resname.eq.'MIN') then
      FOS_MIN(1:3) = readv(1:3)
    else if (resname.eq.'CYT') then
      FOS_CYT(1:3) = readv(1:3)
    else if (resname.eq.'THY') then
      FOS_THY(1:3) = readv(1:3)
    else if (resname.eq.'URA') then
      FOS_URA(1:3) = readv(1:3)
    else if (resname.eq.'GUA') then
      FOS_GUA(1:3) = readv(1:3)
    else if (resname.eq.'ADE') then
      FOS_ADE(1:3) = readv(1:3)
    else if (resname.eq.'K+ ') then
      FOS_K(1:3) = readv(1:3)
    else if (resname.eq.'BR-') then
      FOS_BR(1:3) = readv(1:3)
    else if (resname.eq.'CS+') then
      FOS_CS(1:3) = readv(1:3)
    else if (resname.eq.'I- ') then
      FOS_I(1:3) = readv(1:3)
    else if (resname.eq.'NH4') then
      FOS_NH4(1:3) = readv(1:3)
    else if (resname.eq.'AC-') then
      FOS_AC(1:3) = readv(1:3)
    else if (resname.eq.'GDN') then
      FOS_GDN(1:3) = readv(1:3)
    else if (resname.eq.'1MN') then
      FOS_1MN(1:3) = readv(1:3)
    else if (resname.eq.'2MN') then
      FOS_2MN(1:3) = readv(1:3)
    else if (resname.eq.'LCP') then
      FOS_LCP(1:3) = readv(1:3)
    else if (resname.eq.'NO3') then
      FOS_NO3(1:3) = readv(1:3)
    else if (resname.eq.'BEN') then
      FOS_BEN(1:3) = readv(1:3)
    else if (resname.eq.'NAP') then
      FOS_NAP(1:3) = readv(1:3)
    else if (resname.eq.'PUR') then
      FOS_PUR(1:3) = readv(1:3)
    else if (resname.eq.'O2 ') then
      FOS_O2(1:3) = readv(1:3)
   ! Martin : added the side chain parts   
    else if (resname.eq.'CH2') then
      FOS_CH2(1:3) = readv(1:3)
      
      
      
    end if
!
! and now the generic parameters: currently these are only for the
! polypeptides, i.e., for the backbone of polyamides
! future extensions are obviously nucleic acids for which we'd
! read the values for sugar and phosphate groups here
  else if (length.gt.4) then
!
    readv(:) = 0.0
    kk = 3
    next = 1
    call get_reals_from_str(string2,kk,readv,next,nvalsr)
    if (nvalsr.le.0) then
      write(ilog,*) 'Fatal. I/O error while reading FOS data in parameter file (got: ',paramfile(t1:t2),').'
      call fexit()
    end if
    if (nvalsr.eq.1) then
      readv(2) = readv(1)
      readv(3) = 0.0
    else if (nvalsr.eq.2) then
      readv(3) = 0.0
    end if
!
!   note that this covers all possible caps
    if (keyword(1:10).eq.'PEP_BB_NH2') then
      FOS_PEP_BB_NH2(1:3) = readv(1:3)
    else if (keyword(1:10).eq.'PEP_BB_FOR') then
      FOS_PEP_BB_FOR(1:3) = readv(1:3)
    else if (keyword(1:7).eq.'PEP_CNT') then
      FOS_PEP_CNT(1:3) = readv(1:3)
    else if (keyword(1:7).eq.'PEP_CCT') then
      FOS_PEP_CCT(1:3) = readv(1:3)
    else if (keyword(1:7).eq.'PEP_UNT') then
      FOS_PEP_UNT(1:3) = readv(1:3)
    else if (keyword(1:7).eq.'PEP_UCT') then
      FOS_PEP_UCT(1:3) = readv(1:3)
    else if (keyword(1:14).eq.'PEP_PRO_BB_FOR') then
      FOS_PEP_PROBB_FOR(1:3) = readv(1:3)
!   this is dirty, since PEP_PRO_BB Is a subset of PEP_PRO_BB_FOR
!   this option must be the second of the two
    else if (keyword(1:10).eq.'PEP_PRO_BB') then
      FOS_PEP_PROBB(1:3) = readv(1:3)
    else if (keyword(1:11).eq.'PEP_PRO_CNT') then
      FOS_PEP_PROCNT(1:3) = readv(1:3)
    else if (keyword(1:11).eq.'PEP_PRO_UNT') then
      FOS_PEP_PROUNT(1:3) = readv(1:3)
!   this is dirty, since PEP_BB Is a subset of PEP_BB_???: this option must be the last one
    else if (keyword(1:6).eq.'PEP_BB') then
      FOS_PEP_BB(1:3) = readv(1:3)
!   peptide sidechain crosslinks
    else if (keyword(1:7).eq.'LINK_CC') then
      FOS_LINK_CC(1:3) = readv(1:3)
!   finally, nucleic acid backbone values
    else if (keyword(1:7).eq.'NUC_PO4') then
      FOS_NUC_PO4D(1:3) = readv(1:3)
    else if (keyword(1:8).eq.'NUC_HPO4') then
      FOS_NUC_PO4M(1:3) = readv(1:3)
    else if (keyword(1:9).eq.'NUC_SUG_0') then
      FOS_NUC_MTHF(1:3) = readv(1:3)
    else if (keyword(1:9).eq.'NUC_SUG_1') then
      FOS_NUC_RIBO1(1:3) = readv(1:3)
    else if (keyword(1:9).eq.'NUC_SUG_2') then
      FOS_NUC_RIBO2(1:3) = readv(1:3)
    else if (keyword(1:9).eq.'NUC_SUG_3') then
      FOS_NUC_RIBO3(1:3) = readv(1:3)
    end if
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine fos_summary()
!
  use fos
  use iounit
  use params
!
  implicit none
!
  integer i
  RTYPE smv
!
! do some sanity checks and print out some warnings
!
 77   format('Sidechain vs. Model Compound Parameter Mismatch Detected:'&
 &    ,/,'Sidechain      (',a3,'): ',f9.4,' (kcal/mol or cal/mol/K)'&
 &    ,/,'Model compound (',a3,'): ',f9.4,' (kcal/mol or cal/mol/K)')
 76   format('Sidechain vs. Sidechain Parameter Mismatch Detected:'&
 &    ,/,'Sidechain #1   (',a3,'): ',f9.4,' (kcal/mol or cal/mol/K)'&
 &    ,/,'Sidechain #2   (',a3,'): ',f9.4,' (kcal/mol or cal/mol/K)')
 75   format('Backbone vs. Model Compound Parameter Mismatch Detected:'&
 &    ,/,a16,4x,': ',f9.4,' kcal/mol'&
 &    ,/,'Model compound (',a3,'): ',f9.4,' (kcal/mol or cal/mol/K)')
!
  smv = 3.0*(10.0**(-precision(smv)))
!
  do i=1,3
    if ((FOS_PEP_BB(i)-FOS_NMA(i)).gt.smv) then
      write(ilog,75) "GENERIC BACKBONE",FOS_PEP_BB(i),"NMA",FOS_NMA(i)
    end if
    if ((FOS_PEP_BB_NH2(i)-FOS_ACA(i)).gt.smv) then
      write(ilog,75) "1-AMINE BACKBONE",FOS_PEP_BB_NH2(i),"ACA",FOS_ACA(i)
    end if
    if ((FOS_PEP_BB_FOR(i)-FOS_NMF(i)).gt.smv) then
      write(ilog,75) "FORMYL-BACKBONE ",FOS_PEP_BB_FOR(i),"NMF",FOS_NMF(i)
    end if
    if ((FOS_PEP_PROBB(i)-FOS_DMA(i)).gt.smv) then
      write(ilog,75) "PROLINE BACKBONE",FOS_PEP_PROBB(i),"DMA",FOS_DMA(i)
    end if
    if ((FOS_PEP_CCT(i)-FOS_ASP(i)).gt.smv) then
      write(ilog,75) "CHARGED C-TERM. ",FOS_PEP_CCT(i),"AC-",FOS_AC(i)
    end if
    if ((FOS_NVA(i)-FOS_VAL(i)).gt.smv) then
      write(ilog,76) "NVA",FOS_NVA(i),"VAL",FOS_VAL(i)
    end if
    if ((FOS_NLE(i)-FOS_ILE(i)).gt.smv) then
      write(ilog,76) "NLE",FOS_NLE(i),"ILE",FOS_ILE(i)
    end if
    if ((FOS_HIE(i)-FOS_HID(i)).gt.smv) then
      write(ilog,76) "HIE",FOS_HIE(i),"HID",FOS_HID(i)
    end if
    if ((FOS_PRO(i)-FOS_VAL(i)).gt.smv) then
      write(ilog,76) "PRO",FOS_PRO(i),"VAL",FOS_VAL(i)
    end if
    if ((FOS_ALA(i)-FOS_CH4(i)).gt.smv) then
      write(ilog,77) "ALA",FOS_ALA(i),"CH4",FOS_CH4(i)
    end if
    if ((FOS_ASN(i)-FOS_ACA(i)).gt.smv) then
      write(ilog,77) "ASN",FOS_ASN(i),"ACA",FOS_ACA(i)
    end if
    if ((FOS_GLN(i)-FOS_PPA(i)).gt.smv) then
      write(ilog,77) "GLN",FOS_GLN(i),"PPA",FOS_PPA(i)
    end if
    if ((FOS_SER(i)-FOS_MOH(i)).gt.smv) then
      write(ilog,77) "SER",FOS_SER(i),"MOH",FOS_MOH(i)
    end if
    if ((FOS_TYR(i)-FOS_PCR(i)).gt.smv) then
      write(ilog,77) "TYR",FOS_TYR(i),"PCR",FOS_PCR(i)
    end if
!    if ((FOS_XPC(i)-FOS_CYT(i)).gt.smv) then
!      write(ilog,77) "  C",FOS_XPC(i),"CYT",FOS_CYT(i)
!    end if
!    if ((FOS_XPU(i)-FOS_URA(i)).gt.smv) then
!      write(ilog,77) "  U",FOS_XPU(i),"URA",FOS_URA(i)
!    end if
!    if ((FOS_XPT(i)-FOS_THY(i)).gt.smv) then
!      write(ilog,77) "  T",FOS_XPT(i),"THY",FOS_THY(i)
!    end if
!    if ((FOS_XPA(i)-FOS_ADE(i)).gt.smv) then
!      write(ilog,77) "  A",FOS_XPA(i),"ADE",FOS_ADE(i)
!    end if
!    if ((FOS_XPG(i)-FOS_GUA(i)).gt.smv) then
!      write(ilog,77) "  G",FOS_XPG(i),"GUA",FOS_GUA(i)
!    end if
   end do
!
 74   format(a3,': ',2(f9.4,' kcal/mol'),1x,f9.4,' cal/mol/K')
 73   format(a20,': ',2(f9.4,' kcal/mol'),1x,f9.4,' cal/mol/K')

  if (fos_report.EQV..true.) then
    write(ilog,*)
    write(ilog,*) '--- Summary of FOS Parameters (DG_298, DH, DCp) ---'
    write(ilog,*)
    write(ilog,*) '*** Small Model Compounds and Ions'
    write(ilog,*)
    write(ilog,74) "O2 ",FOS_O2(1:3)
    write(ilog,74) "SPC",FOS_SPC(1:3)
    write(ilog,74) "T3P",FOS_T3P(1:3)
    write(ilog,74) "T4P",FOS_T4P(1:3)
    write(ilog,74) "T5P",FOS_T5P(1:3)
    write(ilog,74) "T4E",FOS_T4E(1:3)
    write(ilog,74) "URE",FOS_URE(1:3)
    write(ilog,74) "NMF",FOS_NMF(1:3)
    write(ilog,74) "NMA",FOS_NMA(1:3)
    write(ilog,74) "ACA",FOS_ACA(1:3)
    write(ilog,74) "PPA",FOS_PPA(1:3)
    write(ilog,74) "FOA",FOS_FOA(1:3)
    write(ilog,74) "DMA",FOS_DMA(1:3)
    write(ilog,74) "CH4",FOS_CH4(1:3)
    write(ilog,74) "PRP",FOS_PRP(1:3)
    write(ilog,74) "IBU",FOS_IBU(1:3)
    write(ilog,74) "NBU",FOS_NBU(1:3)
    write(ilog,74) "EMT",FOS_EMT(1:3)
    write(ilog,74) "MOH",FOS_MOH(1:3)
    write(ilog,74) "EOH",FOS_EOH(1:3)
    write(ilog,74) "MSH",FOS_MSH(1:3)
    write(ilog,74) "BEN",FOS_BEN(1:3)
    write(ilog,74) "NAP",FOS_NAP(1:3)
    write(ilog,74) "TOL",FOS_TOL(1:3)
    write(ilog,74) "PCR",FOS_PCR(1:3)
    write(ilog,74) "IMD",FOS_IMD(1:3)
    write(ilog,74) "IME",FOS_IME(1:3)
    write(ilog,74) "MIN",FOS_MIN(1:3)
    write(ilog,74) "CYT",FOS_CYT(1:3)
    write(ilog,74) "URA",FOS_URA(1:3)
    write(ilog,74) "THY",FOS_THY(1:3)
    write(ilog,74) "PUR",FOS_PUR(1:3)
    write(ilog,74) "ADE",FOS_ADE(1:3)
    write(ilog,74) "GUA",FOS_GUA(1:3)
    write(ilog,74) "NA+",FOS_NA(1:3)
    write(ilog,74) "CL-",FOS_CL(1:3)
    write(ilog,74) "K+ ",FOS_K(1:3)
    write(ilog,74) "BR-",FOS_BR(1:3)
    write(ilog,74) "CS+",FOS_CS(1:3)
    write(ilog,74) "I- ",FOS_I(1:3)
    write(ilog,74) "NH4",FOS_NH4(1:3)
    write(ilog,74) "AC-",FOS_AC(1:3)
    write(ilog,74) "GDN",FOS_GDN(1:3)
    write(ilog,74) "1MN",FOS_1MN(1:3)
    write(ilog,74) "2MN",FOS_2MN(1:3)
    write(ilog,74) "LCP",FOS_LCP(1:3)
    write(ilog,74) "NO3",FOS_NO3(1:3)
    write(ilog,*) 
    write(ilog,*) '*** Polypeptide Sidechains'
    write(ilog,*)
    write(ilog,74) "ALA",FOS_ALA(1:3)
    write(ilog,74) "VAL",FOS_VAL(1:3)
    write(ilog,74) "LEU",FOS_LEU(1:3)
    write(ilog,74) "ILE",FOS_ILE(1:3)
    write(ilog,74) "SER",FOS_SER(1:3)
    write(ilog,74) "S1P",FOS_S1P(1:3)
    write(ilog,74) "THR",FOS_THR(1:3)
    write(ilog,74) "T1P",FOS_T1P(1:3)
    write(ilog,74) "CYS",FOS_CYS(1:3)
    write(ilog,74) "CYX",FOS_CYX(1:3)
    write(ilog,74) "PRO",FOS_PRO(1:3)
    write(ilog,74) "PHE",FOS_PHE(1:3)
    write(ilog,74) "TYR",FOS_TYR(1:3)
    write(ilog,74) "TYO",FOS_TYO(1:3)
    write(ilog,74) "Y1P",FOS_Y1P(1:3)
    write(ilog,74) "TRP",FOS_TRP(1:3)
    write(ilog,74) "HID",FOS_HID(1:3)
    write(ilog,74) "HIE",FOS_HIE(1:3)
    write(ilog,74) "HIP",FOS_HIP(1:3)
    write(ilog,74) "ASP",FOS_ASP(1:3)
    write(ilog,74) "ASH",FOS_ASH(1:3)
    write(ilog,74) "ASN",FOS_ASN(1:3)
    write(ilog,74) "GLU",FOS_GLU(1:3)
    write(ilog,74) "GLH",FOS_GLH(1:3)
    write(ilog,74) "GLX",FOS_GLX(1:3)
    write(ilog,74) "GLN",FOS_GLN(1:3)
    write(ilog,74) "MET",FOS_MET(1:3)
    write(ilog,74) "LYS",FOS_LYS(1:3)
    write(ilog,74) "LYD",FOS_LYD(1:3)
    write(ilog,74) "LYX",FOS_LYX(1:3)
    write(ilog,74) "ARG",FOS_ARG(1:3)
    write(ilog,74) "ORN",FOS_ORN(1:3)
    write(ilog,74) "DAB",FOS_DAB(1:3)
    write(ilog,74) "AIB",FOS_AIB(1:3)
    write(ilog,74) "PCA",FOS_PCA(1:3)
    write(ilog,74) "GAM",FOS_GAM(1:3)
    write(ilog,74) "HYP",FOS_HYP(1:3)
    write(ilog,74) "ABA",FOS_ABA(1:3)
    write(ilog,74) "NVA",FOS_NVA(1:3)
    write(ilog,74) "NLE",FOS_NLE(1:3)
    write(ilog,*)
    write(ilog,*) '*** Polypeptide Crosslinks'
    write(ilog,*)
    write(ilog,73) "Disulfide linkage   ",FOS_LINK_CC(1:3)
    write(ilog,*)
    write(ilog,*) '*** Polypeptide Backbone'
    write(ilog,*)
    write(ilog,73) "Generic Backbone    ",FOS_PEP_BB(1:3)
    write(ilog,73) "Proline Backbone    ",FOS_PEP_PROBB(1:3)
    write(ilog,73) "Formylated N-Term.  ",FOS_PEP_BB_FOR(1:3)
    write(ilog,73) "Proline Formyl-N-T. ",FOS_PEP_PROBB_FOR(1:3)
    write(ilog,73) "NH2-Amidated C-Term.",FOS_PEP_BB_NH2(1:3)
    write(ilog,73) "Charged N-Term.     ",FOS_PEP_CNT(1:3)
    write(ilog,73) "Proline charged N-T.",FOS_PEP_PROCNT(1:3)
    write(ilog,73) "Charged C-Term.     ",FOS_PEP_CCT(1:3)
    write(ilog,73) "Neutral N-Term.     ",FOS_PEP_UNT(1:3)
    write(ilog,73) "Proline neutral N-T.",FOS_PEP_PROUNT(1:3)
    write(ilog,73) "Neutral C-Term.     ",FOS_PEP_UCT(1:3)
    write(ilog,*)
    write(ilog,*) '*** Nucleic Acid Sidechains'
    write(ilog,*)
    write(ilog,74) "  C",FOS_XPC(1:3)
    write(ilog,74) "  U",FOS_XPU(1:3)
    write(ilog,74) "  T",FOS_XPT(1:3)
    write(ilog,74) "  A",FOS_XPA(1:3)
    write(ilog,74) "  G",FOS_XPG(1:3)
    write(ilog,*)
    write(ilog,*) '*** Nucleic Acid Backbone'
    write(ilog,*)
    write(ilog,73) "Phosphodiester      ",FOS_NUC_PO4D(1:3)
    write(ilog,73) "Phosphomonoester    ",FOS_NUC_PO4M(1:3)
    write(ilog,73) "2-Methyl-THF        ",FOS_NUC_MTHF(1:3)
    write(ilog,73) "2-Methyl-THF + 1xOH ",FOS_NUC_RIBO1(1:3)
    write(ilog,73) "2-Methyl-THF + 2xOH ",FOS_NUC_RIBO2(1:3)
    write(ilog,73) "2-Methyl-THF + 3xOH ",FOS_NUC_RIBO3(1:3)
    write(ilog,*)
    write(ilog,*) '---------------------------------'
  end if

!
end
!
!-----------------------------------------------------------------------
!
subroutine check_bndprms()
!
  use params
  use iounit
!
  implicit none
!
  integer i,j,i1,i2,i3,i4,i5,t1,t2,j1,j2,j3,j4,j5,k
!
  call strlims(paramfile,t1,t2)
!
  do i=1,n_bondtyp
    if ((bo_typ(i).gt.3).OR.(bo_typ(i).lt.1)) then
      write(ilog,*) 'Fatal. Unsupported bond type selection in list &
 &of unique bond length potentials.'
      call fexit()
    end if
  end do
  do i=1,bo_lstsz
    i1 = bo_lst(i,1)
    i2 = bo_lst(i,2)
    do j=i+1,bo_lstsz
      j1 = bo_lst(j,1)
      j2 = bo_lst(j,2)
      if (((i1.eq.j1).AND.(i2.eq.j2)).OR.&
 &        ((i1.eq.j2).AND.(i2.eq.j1))) then
        write(ilog,*) 'Fatal. Multiple definitions of bond length po&
 &tential for bonded types ',i1,' and ',i2,' (in: ',&
 &paramfile(t1:t2),').'
        call fexit()
      end if
    end do
  end do
!
  do i=1,n_angltyp
    if ((ba_typ(i).gt.3).OR.(ba_typ(i).lt.1)) then
      write(ilog,*) 'Fatal. Unsupported angle type selection in list&
 & of unique bond angle potentials.'
      call fexit()
    end if
  end do
  
  do i=1,ba_lstsz
    i1 = ba_lst(i,1)
    i2 = ba_lst(i,2)
    i3 = ba_lst(i,3)
    do j=i+1,ba_lstsz
      j1 = ba_lst(j,1)
      j2 = ba_lst(j,2)
      j3 = ba_lst(j,3)
      if (i2.eq.j2) then
        if (((i1.eq.j1).AND.(i3.eq.j3)).OR.&
 &          ((i1.eq.j3).AND.(i3.eq.j1))) then
          write(ilog,*) 'Fatal. Multiple definitions of bond angle p&
 &otential for bonded types ',i1,', ',i2,' and ',i3,' (in: ',&
 &paramfile(t1:t2),').'
          call fexit()
        end if
      end if
    end do
  end do
!
  do i=1,n_torstyp
    if ((di_typ(i).gt.3).OR.(di_typ(i).lt.1)) then
      write(ilog,*) 'Fatal. Unsupported torsion type selection in li&
 &st of unique dihedral angle potentials.'
      call fexit()
    end if
  end do
  do i=1,di_lstsz
    i1 = di_lst(i,1)
    i2 = di_lst(i,2)
    i3 = di_lst(i,3)
    i4 = di_lst(i,4)
    do j=i+1,di_lstsz
      j1 = di_lst(j,1)
      j2 = di_lst(j,2)
      j3 = di_lst(j,3)
      j4 = di_lst(j,4)
      if ((i2.eq.j2).AND.(i3.eq.j3)) then
        if ((i1.eq.j1).AND.(i4.eq.j4)) then
          write(ilog,*) 'Fatal. Multiple definitions of torsional po&
 &tential for bonded types ',i1,', ',i2,', ',i3,' and ',i4,' (in: ',&
 &paramfile(t1:t2),').'
          call fexit()
        end if
      end if
      if ((i2.eq.j3).AND.(i3.eq.j2)) then
        if ((i1.eq.j4).AND.(i4.eq.j1)) then
          write(ilog,*) 'Fatal. Multiple definitions of torsional po&
 &tential for bonded types ',i1,', ',i2,', ',i3,' and ',i4,' (in: ',&
 &paramfile(t1:t2),').'
          call fexit()
        end if
      end if
    end do
  end do
!
  do i=1,n_cmstyp
    if ((cm_typ(i,1).gt.4).OR.(cm_typ(i,1).lt.1)) then
      write(ilog,*) 'Fatal. Unsupported CMAP type selection in li&
 &st of unique CMAP potentials.'
      call fexit()
    end if
  end do
  do i=1,cm_lstsz
    i1 = cm_lst(i,1)
    i2 = cm_lst(i,2)
    i3 = cm_lst(i,3)
    i4 = cm_lst(i,4)
    i5 = cm_lst(i,5)
    do j=i+1,cm_lstsz
      j1 = cm_lst(j,1)
      j2 = cm_lst(j,2)
      j3 = cm_lst(j,3)
      j4 = cm_lst(j,4)
      j5 = cm_lst(j,5)
      if ((i2.eq.j2).AND.(i3.eq.j3).AND.(i4.eq.i4)) then
        if ((i1.eq.j1).AND.(i5.eq.j5)) then
          write(ilog,*) 'Fatal. Multiple definitions of CMAP po&
 &tential for bonded types ',i1,', ',i2,', ',i3,', ',i4,' and ',i5,' (in: ',&
 &paramfile(t1:t2),').'
          call fexit()
        end if
      end if
      if ((i2.eq.j4).AND.(i3.eq.j3).AND.(i4.eq.i2)) then
        if ((i1.eq.j5).AND.(i5.eq.j1)) then
          write(ilog,*) 'Fatal. Multiple definitions of CMAP po&
 &tential for bonded types ',i1,', ',i2,', ',i3,', ',i4,' and ',i5,' (in: ',&
 &paramfile(t1:t2),').'
          call fexit()
        end if
      end if
    end do
  end do

!
  do i=1,impt_lstsz
    i1 = impt_lst(i,1)
    i2 = impt_lst(i,2)
    i3 = impt_lst(i,3)
    i4 = impt_lst(i,4)
    if ((i2.ne.i3).AND.(i2.ne.i4).AND.(i3.ne.i4)) then
      do j=i+1,impt_lstsz
        j1 = impt_lst(j,1)
        j2 = impt_lst(j,2)
        j3 = impt_lst(j,3)
        j4 = impt_lst(j,4)
        if ((i1.eq.j1).AND.(i2.eq.j2)) then
          if ((i3.eq.j3).AND.(i4.eq.j4)) then
            write(ilog,*) 'Fatal. Multiple definitions of improper dih&
 &edral potential for bonded types ',i1,', ',i2,', ',i3,' and ',i4,'&
 & (in: ',paramfile(t1:t2),').'
            call fexit()
          end if
        end if
        if ((i1.eq.j1).AND.(i2.eq.j3)) then 
          if ((i3.eq.j2).AND.(i4.eq.j4)) then
            write(ilog,*) 'Fatal. Multiple definitions of improper dih&
 &edral potential for bonded types ',i1,', ',i2,', ',i3,' and ',i4,'&
 & (in: ',paramfile(t1:t2),').'
            call fexit()
          end if
        end if
      end do
    else if ((i2.eq.i3).AND.(i2.eq.i4)) then
      k = 0
      do j=1,impt_lstsz
        if (i.eq.j) cycle
        j1 = impt_lst(j,1)
        j2 = impt_lst(j,2)
        j3 = impt_lst(j,3)
        j4 = impt_lst(j,4)
        if ((j2.eq.j3).AND.(j2.eq.j4).AND.(j2.eq.i2)) k = k + 1
        if (k.gt.3) then
          write(ilog,*) 'Fatal. Too many redundant definitions of improper dih&
 &edral potential for bonded types ',i1,' and 3x ',i2,' (in: ',paramfile(t1:t2),').'
          call fexit()
        end if
      end do
    else
      if (i2.eq.i3) k = i4
      if (i2.eq.i4) k = i3
      if (i3.eq.i4) k = i2
      do j=1,impt_lstsz
        if (i.eq.j) cycle
        j1 = impt_lst(j,1)
        j2 = impt_lst(j,2)
        j3 = impt_lst(j,3)
        j4 = impt_lst(j,4)
!       this will allow x y z y and x z y y as separate entries in addition to x z z y
        if ((i1.eq.j1).AND.(i2.eq.j2)) then
          if ((i3.eq.j3).AND.(i4.eq.j4)) then
            write(ilog,*) 'Fatal. Multiple definitions of degenerate, improper dih&
 &edral potential for bonded types ',i1,', ',i2,', ',i3,' and ',i4,'&
 & (in: ',paramfile(t1:t2),').'
            call fexit()
          end if
        end if
      end do
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_cmaps()
!
  use iounit
  use params
!
  implicit none
!
  integer iomessage,iu,i,j,k,t11,t12,t21,t22,t31,t32,freeunit
  character(2*MAXSTRLEN) infi
  logical exists
!
  if (n_cmstyp.le.0) return
!
  call strlims(cm_dir,t11,t12)
  do i=1,n_cmstyp
    call strlims(cm_file(i),t21,t22)
    infi = cm_dir(t11:t12)//cm_file(i)(t21:t22)
    call strlims(infi,t31,t32)
    inquire(file=infi(t31:t32),exist=exists)
    if (exists.EQV..false.) then
      write(ilog,*) 'Fatal. File with CMAP corrections does not exist (&
 & for ',infi(t31:t32),'). Check directory and/or file names.'
      call fexit()
    end if
!
    iu = freeunit()
    open(unit=iu,file=infi(t31:t32),status='old')
!
    do j=1,cm_typ(i,2)
      read(iu,*,iostat=iomessage) (cm_par(i,j,k,1),k=1,cm_typ(i,2))
      if (iomessage.eq.-1) then
        write(ilog,*) 'Fatal. Got empty/incomplete file for CMAP corrections &
 &(',infi(t31:t32),'). Expected ',cm_typ(i,2),' rows and columns.'
        call fexit()
      else if (iomessage.eq.2) then
        write(ilog,*) 'Fatal. I/O error while processing input for CMAP corrections &
 &(',infi(t31:t32),').'
        call fexit()
      end if
    end do
!
    close(unit=iu)
!
!   setup basic (constant) CMAP parameters 
    call cmap_derivatives(i)
!
  end do
!
end
!
!----------------------------------------------------------------------------------
!
subroutine cubic_spline_PBC(lori,vv,spac,szi,fvals,fddvals,splval,splvald)
!
  use params
!
  implicit none
!
  integer flr,cel,szi
  RTYPE vv,lori,spac,dum,frac,frac2,splval,splvald
  RTYPE fddvals(szi),fvals(szi)
!
  dum = (vv - lori)/spac
  flr = floor(dum)
  frac = dum - flr
  frac2 = 1.0 - frac
!
  flr = flr + 1
  do while (flr.lt.0) 
    flr = flr + szi
  end do
  cel = flr + 1
  do while (cel.gt.szi) 
    cel = cel - szi
  end do
  if (cel.gt.1) flr = cel - 1
!
  splval = frac*fvals(cel) + frac2*fvals(flr) + (spac*spac/6.)*&
 & ((frac*frac*frac - frac)*fddvals(cel) + (frac2*frac2*frac2 - frac2)*fddvals(flr))
!
  splvald = (fvals(cel) - fvals(flr))/spac + (spac/6.)*&
 & ((3.0*frac*frac - 1.0)*fddvals(cel) - (3.0*frac2*frac2 - 1.0)*fddvals(flr))
!
end
!
!----------------------------------------------------------------------------------
!
subroutine cmap_derivatives(which)
!
! here, we wish to get (rough) approximations for each node point to dv/dx1, dv/dx2
! and dv/dx1dx2
!
  use params
!
  implicit none
!
  integer which,i,ii,jj,grhf,starti,iip,iim
  RTYPE, ALLOCATABLE:: arr1(:),arr2(:),arr3(:),gridder(:,:,:),helper2(:,:),helper1(:,:,:)
  RTYPE pvv,invres2,ori(2),val,grval
  logical atrue
!
  atrue = .true.
!
  allocate(arr1(2*cm_typ(which,2)))
  allocate(arr2(2*cm_typ(which,2)))
  allocate(arr3(2*cm_typ(which,2)))
  allocate(gridder(cm_typ(which,2),cm_typ(which,2),4))
!
! in general, the procedure does the following:
! - unwrap an individual axis (horizontal, vertical) by duplication to
!   avoid boundary artifacts
! - use an iterative procedure to get a consistent set of one-dimensional spline parameters
!   -> this resembles a forward finite difference equation (essentially Taylor expansion)
!      to get a "maximum information" estimate of the second derivative 
! - use that estimate to create a piecewise cubic spline fxn along that axis
! - analytically obtain the grid-point derivatives from the spline
! - repeat procedure usign grid-point derivatives as axis in orthogonal direction (-> cross-derivative)
! - use the complete set of grid-point derivatives (3*N*N) + original data
!   to construct the 16 parameters for the bicubic spline
! - store
!
!  if (cm_typ(which,1).eq.4) then
    ori(1) = -180.0
    ori(2) = -180.0
!  else
!    ori(1) = -180.0 + 0.5*cm_par2(which)
!    ori(2) = -180.0 + 0.5*cm_par2(which)
!  end if
!
  invres2 = 1.0/(cm_par2(which)*cm_par2(which))
  if (mod(cm_typ(which,2),2).eq.0) then
    grhf = cm_typ(which,2)/2
    starti = grhf + 1
  else
    grhf = (cm_typ(which,2)+1)/2
    starti = grhf
  end if
!
  do jj=1,cm_typ(which,2)
    arr1(1) = 0.0
    arr3(1) = 0.0
    do i=2,2*cm_typ(which,2)-1
      ii = i + grhf
      do while (ii.gt.cm_typ(which,2))
        ii = ii - cm_typ(which,2)
      end do
      iip = ii + 1
      if (iip.gt.cm_typ(which,2)) iip = iip - cm_typ(which,2)
      iim = ii - 1
      if (iim.le.0) iim = iim + cm_typ(which,2)
      pvv = 0.5*arr1(i-1) + 2.0
      arr1(i) = -0.5/pvv
      arr2(i) = invres2*(cm_par(which,iip,jj,1) - 2.0*cm_par(which,ii,jj,1) + cm_par(which,iim,jj,1))
      arr3(i) = (3.0*arr2(i) - 0.5*arr3(i-1))/pvv
    end do
    arr1(2*cm_typ(which,2)) = 0.0
    do i=2*cm_typ(which,2)-1,1,-1
      arr1(i) = arr1(i)*arr1(i+1) + arr3(i)
    end do 
    do i=1,cm_typ(which,2)
      arr2(i) = cm_par(which,i,jj,1)
    end do
    do i=1,cm_typ(which,2)
      grval = ori(1) + (i-1)*cm_par2(which)
      call cubic_spline_PBC(ori(1),grval,cm_par2(which),cm_typ(which,2),arr2(1:cm_typ(which,2)),&
 &arr1(starti:starti+cm_typ(which,2)-1),val,gridder(i,jj,1))
    end do
  end do
!
  do jj=1,cm_typ(which,2)
    arr1(1) = 0.0
    arr3(1) = 0.0
    do i=2,2*cm_typ(which,2)-1
      ii = i + grhf
      do while (ii.gt.cm_typ(which,2))
        ii = ii - cm_typ(which,2)
      end do
      iip = ii + 1
      if (iip.gt.cm_typ(which,2)) iip = iip - cm_typ(which,2)
      iim = ii - 1
      if (iim.le.0) iim = iim + cm_typ(which,2)
      pvv = 0.5*arr1(i-1) + 2.0
      arr1(i) = -0.5/pvv
      arr2(i) = invres2*(cm_par(which,jj,iip,1) - 2.0*cm_par(which,jj,ii,1) + cm_par(which,jj,iim,1))
      arr3(i) = (3.0*arr2(i) - 0.5*arr3(i-1))/pvv
    end do
    arr1(2*cm_typ(which,2)) = 0.0
    do i=2*cm_typ(which,2)-1,1,-1
      arr1(i) = arr1(i)*arr1(i+1) + arr3(i)
    end do
    do i=1,cm_typ(which,2)
      arr2(i) = cm_par(which,jj,i,1)
    end do
    do i=1,cm_typ(which,2)
      grval = ori(2) + (i-1)*cm_par2(which)
      call cubic_spline_PBC(ori(2),grval,cm_par2(which),cm_typ(which,2),arr2(1:cm_typ(which,2)),&
 &arr1(starti:starti+cm_typ(which,2)-1),val,gridder(jj,i,2))
    end do
  end do
!
  do jj=1,cm_typ(which,2)
    arr1(1) = 0.0 
    arr3(1) = 0.0
    do i=2,2*cm_typ(which,2)-1
      ii = i + grhf
      do while (ii.gt.cm_typ(which,2))
        ii = ii - cm_typ(which,2)
      end do
      iip = ii + 1
      if (iip.gt.cm_typ(which,2)) iip = iip - cm_typ(which,2)
      iim = ii - 1
      if (iim.le.0) iim = iim + cm_typ(which,2)
      pvv = 0.5*arr1(i-1) + 2.0
      arr1(i) = -0.5/pvv
      arr2(i) = invres2*(gridder(iip,jj,2) - 2.0*gridder(ii,jj,2) + gridder(iim,jj,2))
      arr3(i) = (3.0*arr2(i) - 0.5*arr3(i-1))/pvv
    end do
    arr1(2*cm_typ(which,2)) = 0.0
    do i=2*cm_typ(which,2)-1,1,-1
      arr1(i) = arr1(i)*arr1(i+1) + arr3(i)
    end do
    do i=1,cm_typ(which,2)
      arr2(i) = gridder(i,jj,2)
    end do
    do i=1,cm_typ(which,2)
      grval = ori(1) + (i-1)*cm_par2(which)
      call cubic_spline_PBC(ori(1),grval,cm_par2(which),cm_typ(which,2),arr2(1:cm_typ(which,2)),&
 &arr1(starti:starti+cm_typ(which,2)-1),val,gridder(i,jj,3))
    end do
  end do
!
  allocate(helper2(cm_typ(which,2),cm_typ(which,2)))
  helper2(:,:) = cm_par(which,1:cm_typ(which,2),1:cm_typ(which,2),1)
  allocate(helper1(cm_typ(which,2),cm_typ(which,2),16))
  call bicubic_spline_set_params(cm_typ(which,2),cm_typ(which,2),cm_par2(which),cm_par2(which),&
 &                               helper2(:,:),gridder(:,:,1),gridder(:,:,2),gridder(:,:,3),helper1(:,:,:),atrue)
  cm_par(which,1:cm_typ(which,2),1:cm_typ(which,2),2:17) = helper1(:,:,:)
  deallocate(helper1)
  deallocate(helper2)
!
  deallocate(arr1)
  deallocate(arr2)
  deallocate(arr3)
  deallocate(gridder)
!
 66 format(i6,1x,i6,2x,50(f12.6,1x))
 67 format(50(f12.6))
!
end
!
!----------------------------------------------------------------------------------
!

