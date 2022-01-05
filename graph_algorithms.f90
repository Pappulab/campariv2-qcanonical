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
! CONTRIBUTIONS: Nicolas Bloechliger                                       !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!--------------------------------------------------------------------------------------------
!
subroutine gen_graph_from_clusters(it,nbasins,prtclu)
!
  use clusters
  use mpistuff
  use iounit
  use system
!
  implicit none
!
  integer i,j,k,iu,freeunit,ii,jj,ll,nbasins,ik,jk,cpibu
  integer, ALLOCATABLE:: iv1(:)
  character(MAXSTRLEN) fn
  logical foundit,exists,prtclu
#ifdef ENABLE_MPI
  character(3) nod
  integer tl
#endif
  type(t_scluster) it(nbasins)
  RTYPE dis,dval
  RTYPE, ALLOCATABLE:: topss(:,:)
!
  if (prtclu.EQV..true.) then
!
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      tl = 3
      call int2str(myrank,nod,tl)
      fn = 'N_'//nod(1:tl)//'_STRUCT_CLUSTERING.clu'
    else if (use_MPIAVG.EQV..true.) then
      fn = 'STRUCT_CLUSTERING.clu'
    end if
#else
    fn = 'STRUCT_CLUSTERING.clu'
#endif
!
    call strlims(fn,ii,jj)
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      iu = freeunit()
      open(unit=iu,file=fn(ii:jj),status='old',position='append')
      close(unit=iu,status='delete')
    end if
    iu=freeunit()
    open(unit=iu,file=fn(ii:jj),status='new')
  end if
!
  allocate(iv1(cstored))
  iv1(:) = 0
  ll = 0
  do j=1,nbasins
    do k=1,it(j)%nmbrs
      if ((it(j)%snaps(k).le.0).OR.(it(j)%snaps(k).gt.cstored)) then
        write(ilog,*) 'Fatal. Encountered nonexistent snapshot in cluster member list. This is a bug.'
        call fexit()
      end if
      iv1(it(j)%snaps(k)) = j
    end do
    it(j)%nodewt = (1.0*it(j)%nmbrs)/(1.0*cstored)
  end do
 12 format(i12)
  if (prtclu.EQV..true.) then
    do i=1,cstored
      write(iu,12) iv1(i)
    end do
    close(unit=iu)
  end if
  cpibu = cprogindstart
  if ((cprogindstart.eq.-1).OR.(cprogindstart.eq.-3)) then
    cprogindstart = iv1(1) ! the cluster that the first snapshot belongs to
    write(ilog,*) 'Setting target node to be #',cprogindstart,' as it contains the first snapshot.'
  end if
!
  do i=1,nbasins
    if (allocated(it(i)%wghtsnb).EQV..true.) deallocate(it(i)%wghtsnb)
    if (allocated(it(i)%lensnb).EQV..true.) deallocate(it(i)%lensnb)
!    if (allocated(it(i)%fewtsnb).EQV..true.) deallocate(it(i)%fewtsnb)
    if (allocated(it(i)%map).EQV..true.) deallocate(it(i)%map)
    if (allocated(it(i)%lstnb).EQV..true.) deallocate(it(i)%lstnb)
    if (allocated(it(i)%flwnb).EQV..true.) deallocate(it(i)%flwnb)
    it(i)%nbalsz = 2
    allocate(it(i)%wghtsnb(it(i)%nbalsz,2))
    allocate(it(i)%lensnb(it(i)%nbalsz,2))
!    allocate(it(i)%fewtsnb(it(i)%nbalsz,2))
    it(i)%wghtsnb(:,:) = 0
    it(i)%lensnb(:,:) = 0.0
!    it(i)%fewtsnb(:,:) = 0.0
    allocate(it(i)%lstnb(it(i)%nbalsz))
    allocate(it(i)%map(it(i)%nbalsz))
    it(i)%nb = 0
  end do
  k = 1
  do while (iv1(k).eq.0)
    k = k + 1
  end do
  ll = iv1(k)
  itbrklst = 1
  do i=k+1,cstored
!
    if (ntbrks2.ge.itbrklst) then
      if (trbrkslst(itbrklst).eq.(i-1)) then
        itbrklst = itbrklst + 1
        ll = iv1(i)
        cycle
      end if
    end if
    call snap_to_snap_d(dis,i-1,i)
    dval = sqrt(dis*dis)
!
    if (iv1(i).eq.0) cycle
!
    if ((iv1(i).gt.0).AND.(ll.gt.0)) then
      foundit = .false.
      do j=1,it(iv1(i))%nb
        if (it(iv1(i))%lstnb(j).eq.ll) then
          it(iv1(i))%wghtsnb(j,2) = it(iv1(i))%wghtsnb(j,2) + 1
          it(iv1(i))%lensnb(j,2) = it(iv1(i))%lensnb(j,2) + dval
!          it(iv1(i))%fewtsnb(j,2) = it(iv1(i))%fewtsnb(j,2) + 1.0
          foundit = .true.
          ik = j
          exit
        end if
      end do
      if (foundit.EQV..false.) then
        it(iv1(i))%nb = it(iv1(i))%nb + 1
        ik = it(iv1(i))%nb
        if (it(iv1(i))%nb.gt.it(iv1(i))%nbalsz) then
          call scluster_resizenb(it(iv1(i)))
        end if
!        it(iv1(i))%map(ll) = it(iv1(i))%nb
        it(iv1(i))%lstnb(it(iv1(i))%nb) = ll
        it(iv1(i))%wghtsnb(it(iv1(i))%nb,2) = 1
        it(iv1(i))%wghtsnb(it(iv1(i))%nb,1) = 0
        it(iv1(i))%lensnb(it(iv1(i))%nb,2) = dval
        it(iv1(i))%lensnb(it(iv1(i))%nb,1) = 0.0
!        it(iv1(i))%fewtsnb(it(iv1(i))%nb,2) = 1.0
!        it(iv1(i))%fewtsnb(it(iv1(i))%nb,1) = 0.0
      end if
      foundit = .false.
      do j=1,it(ll)%nb
        if (it(ll)%lstnb(j).eq.iv1(i)) then
          it(ll)%wghtsnb(j,1) = it(ll)%wghtsnb(j,1) + 1
          it(ll)%lensnb(j,1) = it(ll)%lensnb(j,1) + dval
!          it(ll)%fewtsnb(j,1) = it(ll)%fewtsnb(j,1) + 1.0
          foundit = .true.
          jk = j
          exit
        end if
      end do
      if (foundit.EQV..false.) then
        it(ll)%nb = it(ll)%nb + 1
        jk = it(ll)%nb
        if (it(ll)%nb.gt.it(ll)%nbalsz) then
          call scluster_resizenb(it(ll))
        end if
!        it(ll)%map(iv1(i)) = it(ll)%nb
        it(ll)%lstnb(it(ll)%nb) = iv1(i)
        it(ll)%wghtsnb(it(ll)%nb,1) = 1
        it(ll)%wghtsnb(it(ll)%nb,2) = 0
        it(ll)%lensnb(it(ll)%nb,1) = dval
        it(ll)%lensnb(it(ll)%nb,2) = 0.0
!        it(ll)%fewtsnb(it(ll)%nb,1) = 1.0
!        it(ll)%fewtsnb(it(ll)%nb,2) = 0.0
      end if
!     we store the index for the list of node 2 (ll ; jk) in node 1 at the position of the neighbor (iv1(i) ; ik)
      it(iv1(i))%map(ik) = jk
!     we store the index for the list of node 1 (iv1(i) ; ik) in node 2 at the position of the neighbor (ll ; jk)
      it(ll)%map(jk) = ik
    end if
    ll = iv1(i)
  end do
! allocate flow double vector
  do i=1,nbasins
    if (it(i)%nb.gt.0) then
      allocate(it(i)%flwnb(it(i)%nbalsz,2))
      it(i)%flwnb(:,:) = 0
    end if
  end do
  call tarjan_scc(it,nbasins)
  deallocate(iv1)
!
! if requested, reequilibrate all SCCs separately
  if (cequil.gt.0) then  
    dval = 1.0e-7
    call equilibrate_MSM(it,nbasins,dval)
  end if
!
  if (ccfepmode.gt.0) then
    if (cprogindstart.gt.nbasins) then
      write(ilog,*) 'Warning. Requested node for computing cFEP (CPROGINDSTART) does not exist. Skipped ...'
    else if ((maxval(it(1:nbasins)%inscc).ne.minval(it(1:nbasins)%inscc)).AND.(cequil.le.0)) then
      write(ilog,*) 'Warning. Network does not represent a single strongly-connected component, which may &
 &alter meaning of cFEP.'
    else if ((maxval(it(1:nbasins)%inscc).ne.minval(it(1:nbasins)%inscc)).AND.(cequil.gt.0).AND.(cprogindstart.le.0)) then
      write(ilog,*) 'Remark. Computing separate cFEPs for individual and large enough SCCs using largest cluster within each SCC &
 &as reference (separate output files).'
    else if ((maxval(it(1:nbasins)%inscc).ne.minval(it(1:nbasins)%inscc)).AND.(cequil.gt.0).AND.(cprogindstart.gt.0)) then
      write(ilog,*) 'Warning. The computed cFEP will only include the SCC containing the chosen reference cluster.'
    end if 
    dval = 1.0e-7
    if ((maxval(it(1:nbasins)%inscc).ne.minval(it(1:nbasins)%inscc)).AND.(cequil.gt.0).AND.(cprogindstart.le.0)) then
      allocate(topss(maxval(it(1:nbasins)%inscc),3))
      topss(:,:) = 0.0
      do i=1,nbasins
        topss(it(i)%inscc,3) = topss(it(i)%inscc,3) + 1.0
        if (it(i)%nodewt.gt.topss(it(i)%inscc,2)) then
          topss(it(i)%inscc,2) = it(i)%nodewt
          topss(it(i)%inscc,1) = 1.0*i
        end if
      end do
      do i=1,maxval(it(1:nbasins)%inscc)
        if ((topss(i,1).gt.0.0).AND.(topss(i,3).gt.1.0)) then
          j = nint(topss(i,1))
          call cfep(it,nbasins,j,dval,ccfepmode)
        end if
      end do
      deallocate(topss)
    else
      i = max(1,cprogindstart)
      if (cprogindstart.le.nbasins) call cfep(it,nbasins,i,dval,ccfepmode)
    end if
  end if
  cprogindstart = cpibu
!
end
!
!--------------------------------------------------------------------------------
!
! get the number and list of unassigned snapshots that space transitions between basins
! WARNING: incomplete
!
subroutine get_graph_unassigned(it,nbasins,ivm)
!
  use clusters
!
  implicit none
!
  integer i,j,k,ll,nbasins,sep
  integer, ALLOCATABLE:: iv1(:)
  integer ivm(nbasins,nbasins)
  logical foundit
  type(t_scluster) it(nbasins)

  allocate(iv1(cstored))
  iv1(:) = 0
  ll = 0
  do i=1,cstored
    foundit = .false.
    do j=1,nbasins
      do k=1,it(j)%nmbrs
        if (i.eq.it(j)%snaps(k)) then
          foundit = .true.
          iv1(i) = j
          exit
        end if
      end do
      if (foundit.EQV..true.) exit
    end do
  end do
!
  k = 1
  do while (iv1(k).eq.0)
    k = k + 1
  end do
  ll = iv1(k)
  sep = 0
  do i=k+1,cstored
    if (iv1(i).eq.0) then
      sep = sep + 1
      cycle
    end if
!
    if ((iv1(i).gt.0).AND.(ll.gt.0)) then
      ivm(ll,iv1(i)) = ivm(ll,iv1(i)) + sep
      sep = 0
    end if
    ll = iv1(i)
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine graphml_helper_for_clustering(it,nclusters)
!
  use iounit
  use mpistuff
  use clusters
  use molecule
  use sequen
  use aminos
  use system
  use mcsums
  use pdb
!
  implicit none
!
  integer i,j,k,ii,jj,ipdbh,freeunit,nclusters,ilen,jlen,klen
  character(MAXSTRLEN) fn
#ifdef ENABLE_MPI
  character(3) nod
  integer tl
#endif
  character(MAXSTRLEN) istr,jstr,kstr,iistr
  logical exists
  type(t_scluster) it(nclusters)
!

#ifdef ENABLE_MPI
  tl = 3
  call int2str(myrank,nod,tl)
  if (use_REMC.EQV..true.) then
    fn = 'N_'//nod(1:tl)//'_STRUCT_CLUSTERING.graphml'
  else
    fn = 'STRUCT_CLUSTERING.graphml'
  end if
#else
  fn = 'STRUCT_CLUSTERING.graphml'
#endif
!
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    ipdbh = freeunit()
    open(unit=ipdbh,file=fn(ii:jj),status='old',position='append')
    close(unit=ipdbh,status='delete')
  end if
  ipdbh=freeunit()
  open(unit=ipdbh,file=fn(ii:jj),status='new')
!
 77  format(a)
 78 format(a,es11.5,a)
  write(ipdbh,77)  '<?xml version="1.0" encoding="UTF-8"?>'
  write(ipdbh,77) '<graphml xmlns="http://graphml.graphdrawing.org/xmlns"'
  write(ipdbh,*) '    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
  write(ipdbh,*) '    xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">'
  write(ipdbh,*) ' <key id="d0" for="node" attr.name="nweight" attr.type="double"/>'
  write(ipdbh,*) ' <key id="d1" for="node" attr.name="nradius" attr.type="double"/>'
  write(ipdbh,*) ' <key id="d2" for="node" attr.name="ndiamet" attr.type="double"/>'
  write(ipdbh,*) ' <key id="d3" for="node" attr.name="nqualit" attr.type="double"/>'
  write(ipdbh,*) ' <key id="i0" for="node" attr.name="ncenter" attr.type="int"/>'
  write(ipdbh,*) ' <key id="i2" for="node" attr.name="nsccidx" attr.type="int"/>'
  write(ipdbh,*) ' <key id="i1" for="edge" attr.name="eweight" attr.type="int"/>'
  write(ipdbh,*) ' <key id="d4" for="edge" attr.name="egeomsd" attr.type="double"/>'
!  write(ipdbh,*) ' <key id="d5" for="edge" attr.name="efloatw" attr.type="double"/>'
  write(ipdbh,*)
  write(ipdbh,*) ' <graph id="G" edgedefault="directed">'
!
  do i=1,nclusters
    ilen = log10(1.0*i) + 1
    write(istr,*) ilen
! 33 format(a,i<ilen>,a)
    write(ipdbh,'("    <node id=""n",i' // ADJUSTL(istr) // ',""">")') i
    write(ipdbh,78) '      <data key="d0">',it(i)%nodewt,'</data>'
    write(ipdbh,78) '      <data key="d1">',it(i)%radius,'</data>'
    write(ipdbh,78) '      <data key="d2">',it(i)%diam,'</data>'
    write(ipdbh,78) '      <data key="d3">',it(i)%quality,'</data>'
    write(iistr,*) int(log10(1.0*max(it(i)%center,1)) + 1)
    write(ipdbh,'("      <data key=""i0"">",i' // ADJUSTL(iistr) // ',"</data>")') it(i)%center
    write(iistr,*) int(log10(1.0*max(it(i)%inscc,1)) + 1)
    write(ipdbh,'("      <data key=""i2"">",i' // ADJUSTL(iistr) // ',"</data>")') it(i)%inscc
    write(ipdbh,*) '    </node>'
  end do
  k = 1
  do i=1,nclusters
    ilen = log10(1.0*i) + 1
    write(istr,*) ilen
    do j=1,it(i)%nb
      if (it(i)%wghtsnb(j,1).le.0) cycle
      if (k.le.0) then
        klen = 1
      else
        klen = log10(1.0*k) + 1
      end if
      jlen = log10(1.0*max(it(i)%lstnb(j),1)) + 1
      write(jstr,*) jlen
      write(kstr,*) klen
      write(iistr,*) int(log10(1.0*max(it(i)%wghtsnb(j,1),1)) + 1)
! 34   format(a,i<klen>,a,i<ilen>,a,i<jlen>,a)
!      write(ipdbh,34) '     <edge id="e',k,'" source="n',i,'" target="n',it(i)%lstnb(j),'">'
      write(ipdbh,'("     <edge id=""e",i' // ADJUSTL(kstr) // ',""" source=""n",i' // ADJUSTL(istr) // &
 &'""" target=""n",i' // ADJUSTL(jstr) // '""">")') k,i,it(i)%lstnb(j)
      write(ipdbh,'("      <data key=""i1"">",i' // ADJUSTL(iistr) // ',"</data>")') it(i)%wghtsnb(j,1)
      write(ipdbh,78) '      <data key="d4">',it(i)%lensnb(j,1)/(1.0*it(i)%wghtsnb(j,1)),'</data>'
!      write(ipdbh,78) '      <data key="d5">',it(i)%fewtsnb(j,1),'</data>'
!      write(ipdbh,*) '      <data key="i1">',it(i)%wghtsnb(j,1),'</data>'
      write(ipdbh,*) '    </edge>'
      k = k + 1
    end do 
  end do
  write(ipdbh,*) '  </graph>'
  write(ipdbh,77) '</graphml>'
  close(unit=ipdbh)
!
end
!
!-------------------------------------------------------------------------
!
subroutine cfep(it,nbasins,snode,thresher,nmode)
!
  use clusters
  use iounit
  use mpistuff
  use interfaces
  use system
!
  implicit none
!
  integer nbasins,snode,i,j,k,ipdbh,ii,jj,freeunit,tx,nmode,t1,t2
  RTYPE thresher,totew,progger,kr
  type(t_scluster) it(nbasins)
  RTYPE, ALLOCATABLE:: mfptv(:),incr(:),hlp1(:)!,mfpth(:),mfptmsd(:)
  integer, ALLOCATABLE:: cutter(:),iv3(:),iv4(:)
  logical upornot,exists
  character(8) nor
  character(MAXSTRLEN) fn,bitt
#ifdef ENABLE_MPI
  character(3) nod
  integer tl
#endif
!
  allocate(mfptv(nbasins))
!
  if (nmode.eq.1) then
    call iterative_mfpt(it,nbasins,snode,mfptv,thresher)
  else if (nmode.eq.2) then
    call reweight_edges(it,nbasins)
    call equilibrate_MSM(it,nbasins,thresher)
    call iterative_mfpt(it,nbasins,snode,mfptv,thresher)
  else if (nmode.le.7) then
    k = nmode - 2
    call dijkstra_sp(it,nbasins,snode,mfptv,k)
  end if
!
  totew = 0.0
  do j=1,nbasins
    if (it(j)%nmbrs.gt.0) then
      if (cequil.gt.0) then
        totew = totew + it(j)%nodewt
      else
        kr = 0.5*sum(it(j)%wghtsnb(1:it(j)%nb,1)+it(j)%wghtsnb(1:it(j)%nb,2))
        totew = totew + kr
      end if
    end if
  end do
  if (cequil.gt.0) totew = totew*cstored
!
  tx = 8
  call int2str(snode,nor,tx)
  bitt(:) = ' '
  bitt(1:5) = 'MFPT '
  if (nmode.gt.1) then
    if (nmode.eq.2) bitt(1:10) = 'DIFFUSION '
    if (nmode.gt.2) bitt(1:10) = 'SHORTPATH '
  end if
  call strlims(bitt,t1,t2)
#ifdef ENABLE_MPI
  tl = 3
  call int2str(myrank,nod,tl)
  if (use_REMC.EQV..true.) then
    fn = 'N_'//nod(1:tl)//'_'//bitt(t1:t2)//'_CFEP_'//nor(1:tx)//'.dat'
  else
    fn = bitt(t1:t2)//'_CFEP_'//nor(1:tx)//'.dat'
  end if
#else
  fn = bitt(t1:t2)//'_CFEP_'//nor(1:tx)//'.dat'
#endif
!
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    ipdbh = freeunit()
    open(unit=ipdbh,file=fn(ii:jj),status='old',position='append')
    close(unit=ipdbh,status='delete')
  end if
  ipdbh=freeunit()
  open(unit=ipdbh,file=fn(ii:jj),status='new')
!
  allocate(incr(nbasins))
  allocate(iv3(nbasins))
  allocate(iv4(nbasins))
  allocate(hlp1(nbasins))
  incr(:) = mfptv(:)
  do i=1,nbasins
    iv3(i) = i
  end do
  ii = 1
  jj = nbasins
  k = nbasins
  upornot = .true.
  call merge_sort(k,upornot,incr,hlp1,ii,jj,iv3,iv4)
!
  deallocate(incr)
  deallocate(iv3)
!
!  allocate(mfpth(100))
  allocate(cutter(nbasins))
!  allocate(mfptmsd(nbasins))
!
! for all edges 
  cutter(:) = 0
  hlp1(:) = 0.0
!  mfpth(:) = 0.0
!  mfptmsd(:) = 0.0
!  binsz = (maxval(mfptv(:))-minval(mfptv(:)))/100.
!  mfptmsd(1) = 1.0
  do i=1,nbasins
    if ((it(iv4(i))%inscc.ne.it(snode)%inscc).AND.(cequil.gt.0)) cycle
    do j=1,it(iv4(i))%nb
      jj = it(iv4(i))%lstnb(j) ! loop over all edges
      if ((it(jj)%inscc.ne.it(snode)%inscc).AND.(cequil.gt.0)) cycle
      do k=i+1,nbasins ! starting the summation at i insures that mfpt(iv4(i)) is always <= mfpt(iv4(k))
        if (mfptv(jj).ge.mfptv(iv4(k))) then
!         the total cut weight as the sum over all snap-to-snap transitions
          cutter(k) = cutter(k) + it(iv4(i))%wghtsnb(j,1)
!         the mean distance across the cut as a linear average over all snap-to-snap transitions 
          hlp1(k) = hlp1(k) + it(iv4(i))%lensnb(j,1)
!         the MSD in the variable used for ordering across the cut
!          mfptmsd(k) = mfptmsd(k) + it(iv4(i))%wghtsnb(j,1)*(mfptv(jj)-mfptv(iv4(k)))**2
        end if
      end do
    end do
!    k = max(1,min(100,ceiling((mfptv(iv4(i))-minval(mfptv(:)))/binsz)))
!    mfpth(k) = mfpth(k) + 1.0*it(iv4(i))%nmbrs
  end do
!
 14 format(i8,1x,i8,1x,g15.8,1x,i8,1x,g11.4,1x,g13.6,1x,g13.6,1x,g13.6,1x,g13.6)
  progger = it(snode)%nodewt
  jj = 1
  do i=2,nbasins
    if ((it(iv4(i))%inscc.ne.it(snode)%inscc).AND.(cequil.gt.0)) cycle
    jj = jj + 1
    cutter(i) = max(1,cutter(i)) ! we may have missed a single transition in opposite direction
    write(ipdbh,14) jj,it(iv4(i))%center,progger,cutter(i),-(1.0/invtemp)*log(cutter(i)/totew),mfptv(iv4(i)),&
 &                  hlp1(i)/(2.0*cstorecalc*cutter(i))!,&
! &                  -(1.0/invtemp)*log(mfptmsd(i)/normerr),mfptmsd(i)/(2.0*cstorecalc*cutter(i))
!  do i=1,100
!    write(ipdbh,14) i,101-i,minval(mfptv(:))+(i-0.5)*binsz,int(mfpth(i)),-(1.0/invtemp)*log((mfpth(i)+1.)/(1.0*cstored))
    progger = progger + it(iv4(i))%nodewt
  end do
!
  close(unit=ipdbh)
!
  deallocate(cutter)
  deallocate(iv4)
  deallocate(hlp1)
  deallocate(mfptv)
!  deallocate(mfptmsd)
!  deallocate(mfpth)
!
end
!
!--------------------------------------------------------------------------
!
! WARNING: incomplete
!
subroutine reweight_edges(it,nbasins)
!
  use clusters
  use iounit
!
  implicit none
!
  integer nbasins,i,j
  type(t_scluster) it(nbasins)
!
  do i=1,nbasins
    do j=1,it(i)%nb
      if (it(i)%lensnb(j,1).gt.0.0) it(i)%wghtsnb(j,1) = nint(1000.*it(i)%lensnb(j,1))
      if (it(i)%lensnb(j,2).gt.0.0) it(i)%wghtsnb(j,2) = nint(1000.*it(i)%lensnb(j,2))
    end do
  end do
!
end
!
!-------------------------------------------------------------------------------------
!
! assuming some changes have been made to the transition matrix, equilibrate a MSM to convergence (thresh)
!
subroutine equilibrate_MSM(it,nbasins,thresh)
!
  use clusters
  use iounit
!
  implicit none
!
  integer nbasins,iterations,i,j,nsccs
  type(t_scluster) it(nbasins)
  RTYPE thresh,normerr
  RTYPE, ALLOCATABLE:: pvec(:,:),pbase(:,:)
  logical, ALLOCATABLE:: hit(:)
  logical notdone
!
  write(ilog,*) 'Now (re)equilibrating network ...'
!
  nsccs = maxval(it(1:nbasins)%inscc)
  allocate(pvec(nbasins,2))
  allocate(hit(nbasins))
  allocate(pbase(nsccs,2))
  pbase(:,:) = 0.0
  do i=1,nbasins
    pbase(it(i)%inscc,1) = pbase(it(i)%inscc,1) + it(i)%nodewt
  end do
!
  pvec(:,1) = it(1:nbasins)%nodewt
  notdone = .true.
!
  iterations = 0
  do while (notdone.EQV..true.)
!
    iterations = iterations + 1
    pvec(:,2) = 0.0
    pbase(:,2) = 0.0
    hit(:) = .false.
    do i=1,nbasins
      do j=1,it(i)%nb
        if (it(it(i)%lstnb(j))%inscc.eq.it(i)%inscc) then
          pvec(it(i)%lstnb(j),2) = pvec(it(i)%lstnb(j),2) + it(i)%wghtsnb(j,1)/(1.0*cstored)
          hit(it(i)%lstnb(j)) = .true.
        end if
      end do
    end do
!   SCCs of size 1
    do i=1,nbasins
      if (hit(i).EQV..false.) then
        pvec(i,2) = it(i)%nodewt ! unreachable
      end if
      pbase(it(i)%inscc,2) = pbase(it(i)%inscc,2) + pvec(i,2)
    end do
    do i=1,nbasins
      pvec(i,2) = pbase(it(i)%inscc,1)*pvec(i,2)/pbase(it(i)%inscc,2)
    end do
    normerr = sum(abs(pvec(:,1)-pvec(:,2)))/(1.0*nbasins)
    if (normerr.le.thresh) notdone = .false.
    pvec(:,1) = pvec(:,2)
  end do
!
  it(:)%nodewt = pvec(:,1)
!
  deallocate(pvec)
  write(ilog,*) 'Done after ',iterations,' iterations.'
!
end
!
!-------------------------------------------------------------------------
!
! this subroutine computes different variants of shortest path measures for
! reaching arbitrary nodes from a source (snode)
! mode = 1: shortest path by number of edges
! mode = 2: shortest path by total geometric length of edges
! mode = 3: most probable path (minimal in negative logarithm) - this penalizes paths with low likelihood nodes on them
! mode = 4: most probable path assuming that all nodes are populated equally
! mode = 5: most probable path assuming equal capacity for all existing edges
! 
! name of the subroutine may be a misnomer
!
subroutine dijkstra_sp(it,nbasins,snode,plen,mode)
!
  use clusters
  use iounit
!
  implicit none
!
  integer nbasins,mode,snode,iterations,processed,i,j,hits,nohits
  RTYPE netlen,minlen,plen(nbasins),mtmp
  type(t_scluster) it(nbasins)
  integer, ALLOCATABLE:: touched(:)
  logical notdone
!
  write(ilog,*) 'Now computing shortest path (',mode,') from node ',snode,' ...'
!
  allocate(touched(nbasins))
!
  iterations = 0
  plen(:) = HUGE(plen)
  plen(snode) = 0.0
  if ((mode.eq.3).OR.(mode.eq.5)) plen(snode) = -log(it(snode)%nodewt)
  touched(:) = 0
  touched(snode) = 1
  notdone = .true.
  processed = nbasins - 1
  nohits = 0
  minlen = 0.0
!
  do while (notdone.EQV..true.)
!
    iterations = iterations + 1
    hits = 0
    mtmp = 0.001*HUGE(mtmp)
    do i=1,nbasins
      if (touched(i).ne.1) cycle
      do j=1,it(i)%nb
        if (it(i)%wghtsnb(j,1).le.0) cycle
        if (mode.eq.1) then
          netlen = 1.0
        else if (mode.eq.2) then
          netlen = it(i)%lensnb(j,1)/(1.0*it(i)%wghtsnb(j,1))
        else if (mode.eq.3) then
          netlen = -log(it(i)%nodewt*it(i)%wghtsnb(j,1)/(1.0*sum(it(i)%wghtsnb(1:it(i)%nb,1))))
        else if (mode.eq.4) then
          netlen = -log(1.0*it(i)%wghtsnb(j,1)/(1.0*sum(it(i)%wghtsnb(1:it(i)%nb,1))))
        else if (mode.eq.5) then
          netlen = -log(it(it(i)%lstnb(j))%nodewt/(1.0*it(i)%nb))
        end if
        mtmp = min(mtmp,netlen)
        if ((plen(i)+netlen).lt.plen(it(i)%lstnb(j))) then
          hits = hits + 1
          plen(it(i)%lstnb(j)) = plen(i) + netlen
          if (touched(it(i)%lstnb(j)).eq.0) processed = processed - 1
          touched(it(i)%lstnb(j)) = 1
        end if
      end do
      touched(i) = 2
    end do
    if (hits.eq.0) then
      nohits = nohits + 1
    else
      nohits = 0
    end if
    minlen = minlen + mtmp
    if ((maxval(plen(:)).lt.(minlen)).AND.(processed.eq.0)) exit
    if (nohits.eq.1000) exit ! emergency
  end do
!
! fix non-reachable
  netlen = 0.0
  do i=1,nbasins
    if (plen(i).eq.HUGE(plen(1))) cycle
    netlen = max(netlen,plen(i))
  end do
  j = 0
  do i=1,nbasins
    if (plen(i).eq.HUGE(plen(1))) then
      j = j + 1
      plen(i) = netlen + 1.0*j
      if ((j.eq.1).AND.(cequil.le.0)) then
        write(ilog,*) 'Warning. Path length values larger than ',netlen,' correspond to &
 &unreachable nodes (arbitrary values given).'
      end if
    end if
  end do
!
  deallocate(touched)
!
  write(ilog,*) 'Done after ',iterations,' iterations.'
!
end
!
!-------------------------------------------------------------------------
!
! iterative determination of the MFPT to snode
!
subroutine iterative_mfpt(it,nbasins,snode,mfptv,thresher)
!
  use clusters
  use iounit
!
  implicit none
!
  integer nbasins,snode,iterations,i,j,nsccs,dumsum
  RTYPE mfptv(nbasins),thresher,normerr,kr
  type(t_scluster) it(nbasins)
  logical notdone
  RTYPE, ALLOCATABLE:: invp(:),incr(:)
  logical, ALLOCATABLE:: hit(:)
  integer, ALLOCATABLE:: scclst(:)
!
  write(ilog,*) 'Now iteratively computing mean first passage time to node ',snode,' ...'
!
  nsccs = maxval(it(1:nbasins)%inscc)
  allocate(invp(nbasins))
  allocate(incr(nbasins))
  allocate(hit(nbasins))
  allocate(scclst(nsccs))
!
  do j=1,nbasins
    if (it(j)%nodewt.gt.0) then
      invp(j) = 1.0/(cstored*it(j)%nodewt)
    else
      invp(j) = 0.0
    end if
  end do
  notdone = .true.
  iterations = 0
  hit(:) = .false.
  scclst(:) = 0
! we need this to flag SCCs reachable from snode (must be excluded from MFPT iteration since they work as sinks -> infinite MFPT)
  if ((cequil.le.0).AND.(nsccs.gt.1)) then
    mfptv(:) = 0.0
    mfptv(snode) = 1.0
    dumsum = -1
    scclst(it(snode)%inscc) = 1
    do while (dumsum.ne.sum(scclst))
      dumsum = sum(scclst)
      do i=1,nbasins
        do j=1,it(i)%nb
          if (it(i)%wghtsnb(j,1).le.0) cycle
          if ((it(it(i)%lstnb(j))%inscc.ne.it(i)%inscc).AND.(scclst(it(i)%inscc).gt.0).AND.&
 &            (scclst(it(it(i)%lstnb(j))%inscc).eq.0)) then
            scclst(it(it(i)%lstnb(j))%inscc) = scclst(it(i)%inscc) + 1
          end if
        end do
      end do
    end do
  end if
  hit(:) = .false.
  mfptv(:) = 1.0
  mfptv(snode) = 0.0
  do while (notdone)
    iterations = iterations + 1
    incr(:) = 1.0
    do i=1,nbasins
      if (i.eq.snode) cycle
      if (cequil.gt.0) then
        if (it(i)%inscc.ne.it(snode)%inscc) cycle
      else
        if (scclst(it(i)%inscc).gt.1) cycle
      end if
      do j=1,it(i)%nb
        if (cequil.gt.0) then
          if (it(it(i)%lstnb(j))%inscc.ne.it(snode)%inscc) cycle
        else
          if (scclst(it(it(i)%lstnb(j))%inscc).gt.1) cycle
        end if
        if (it(i)%wghtsnb(j,1).gt.0) then
          incr(i) = incr(i) + it(i)%wghtsnb(j,1)*mfptv(it(i)%lstnb(j))*invp(i)
          hit(i) = .true.
        end if
      end do
    end do
    incr(snode) = 0.0
    normerr = sum(abs(incr(:)-mfptv(:)))/(1.0*nbasins)
    mfptv(:) = incr(:)
    if (normerr.lt.thresher) exit
    if (maxval(mfptv(:)).ge.cstored) exit ! emergency exit 
  end do
!
  kr = maxval(mfptv(:))
  j = 0
  do i=1,nbasins
    if ((hit(i).EQV..false.).AND.(i.ne.snode)) then
      if ((j.eq.0).AND.(cequil.le.0)) then
        write(ilog,*) 'Warning. MFPT values larger than ',kr,' correspond to &
 &ineligible nodes (arbitrary values given).'
      end if
      kr = kr + 1.0
      mfptv(i) = kr
      j = j + 1
    end if
  end do
!
  deallocate(scclst)
  deallocate(hit)
  deallocate(invp)
  deallocate(incr)
!
  write(ilog,*) 'Done after ',iterations,' iterations.'
!
end
!
!-------------------------------------------------------------------------------------
!
subroutine min_st_cut(it,nbasins,snode,tnode,cutmode,cutval,spart,spartsz)
!
  use clusters
  use iounit
!
  implicit none
!
  integer ii,i,j,k,tnodei
  integer nbasins,snode,tnode,cutval,spart(nbasins),spartsz
  integer startnode,compcut,setincr,oldsetsz,curnode,getexflow,setsz
  integer, ALLOCATABLE:: setmap(:)
  logical cutmode,notdone
  type(t_scluster) it(nbasins)
!
  if (snode.eq.tnode) then
    write(ilog,*) 'Fatal. Called min_st_cut(...) with identical source and target nodes. This is a bug.'
    call fexit()
  end if
!
  do i=1,nbasins
    it(i)%ldis = 0
    it(i)%active = 0
    it(i)%flwnb(:,:) = 0
    it(i)%rflw = 0
  end do
!
  
!
  startnode = -1
!  write(*,*) snode,tnode,it(snode)%nb,it(tnode)%nb
!  write(*,*) 'SL'
!  write(*,*) it(snode)%lstnb(1:it(snode)%nb)
!   write(*,*) 'SW'
!  write(*,*) it(snode)%wghtsnb(1:it(snode)%nb,:)
!  write(*,*) 'TL'
!  write(*,*) it(tnode)%lstnb(1:it(tnode)%nb)
!  write(*,*) 'TW'
!  write(*,*) it(tnode)%wghtsnb(1:it(tnode)%nb,:)


  do i=1,it(snode)%nb
    if (it(snode)%lstnb(i).eq.snode) cycle
    j = it(snode)%lstnb(i)
    if (j.eq.tnode) tnodei = i
    it(snode)%flwnb(i,1) = it(snode)%flwnb(i,1) + it(snode)%wghtsnb(i,1)
    k = it(snode)%map(i) ! k = it(j)%map(snode)
    it(j)%flwnb(k,2) = it(j)%flwnb(k,2) + it(snode)%wghtsnb(i,1)
    if ((it(snode)%wghtsnb(i,1).gt.0).AND.(j.ne.tnode)) then
      startnode = j
      it(j)%active = 1
    end if
  end do
  it(snode)%ldis = nbasins
!
  if (startnode.eq.-1) then
    cutval = it(snode)%wghtsnb(it(tnode)%map(tnodei),1) ! it(snode)%wghtsnb(it(snode)%map(tnode),1)
!   wait
  else
    notdone = .true.
    do while (notdone.EQV..true.) 
      call process_pushpreflow(it,startnode,snode,tnode,cutmode,curnode,nbasins)
      if (curnode.eq.-1) then
        notdone = .false.
        do i=1,nbasins
          if ((i.ne.snode).AND.(i.ne.tnode).AND.(it(i)%active.gt.0)) then
            startnode = i
            notdone = .true.
          end if
        end do
      else
        startnode = curnode
      end if
    end do
  end if
!
  setincr = 1
  it(tnode)%rflw = 1
  setsz = 1
  allocate(setmap(nbasins))
  setmap(:) = 0
  setmap(tnode) = 1
  do while (setincr.gt.0) 
    oldsetsz = setsz
    do i=1,nbasins
      if (it(i)%rflw.eq.0) cycle
      call process_pushresflow(it,i,nbasins,setmap,setsz)
    end do
    setincr = setsz - oldsetsz
  end do
  spartsz = 0
  do i=1,nbasins
    if (setmap(i).eq.0) then
      spartsz = spartsz + 1
      spart(spartsz) = i
    end if
  end do
!  write(*,*) 'SPART: ',spart(1:spartsz)
!
  cutval = 0
  do ii=1,nbasins
    do i=1,it(ii)%nb
      if (it(ii)%lstnb(i).eq.ii) cycle
      j = it(ii)%lstnb(i)
      if (j.lt.ii) then
        if ((setmap(ii).eq.0).AND.(setmap(j).eq.1)) then
          cutval = cutval + it(ii)%wghtsnb(i,1)
        else if ((setmap(j).eq.0).AND.(setmap(ii).eq.1)) then
          cutval = cutval + it(ii)%wghtsnb(i,2)
        end if
      end if
    end do
  end do
  compcut = getexflow(it,tnode)
  if (cutval.ne.compcut) then
    write(ilog,*) 'Fatal. Push-relabel algorithm failure in min_st_cut(...). &
 &Sum of cuts is ',cutval,'  but excess flow at target is ',compcut,'. This is a bug.'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------------------------
!
subroutine process_pushpreflow(it,inode,snode,tnode,cutmode,curnode,nbasins)
!  
  use clusters
  use iounit
!
  implicit none
!
  integer inode,snode,tnode,curnode,nbasins
  integer i,j,k,minldis,maxfl,avail,exflow
  logical cutmode
  integer block(nbasins)
  type(t_scluster) it(nbasins)
!
  curnode = -1
  exflow = 0
  avail = it(inode)%nb
!
  do i=1,it(inode)%nb
    if (it(inode)%lstnb(i).eq.inode) then
      block(i) = 0
      avail = avail - 1
      cycle
    end if
    j = it(inode)%lstnb(i) 
    block(i) = 0
    exflow = exflow + it(inode)%flwnb(i,2) - it(inode)%flwnb(i,1)
    if (it(j)%ldis.ne.(it(inode)%ldis-1)) then
      block(i) = 1
      avail = avail - 1
    end if
  end do
!  write(*,*) 'Excess flow for ',inode,' is ',exflow
!
  do while (exflow.ne.0)
    if (avail.gt.0) then
      do i=1,it(inode)%nb
        if (it(inode)%lstnb(i).eq.inode) cycle
        j = it(inode)%lstnb(i)
        k = it(inode)%map(i) ! it(j)%map(inode)
        if (block(i).eq.0) then
          maxfl = it(inode)%wghtsnb(i,1) - it(inode)%flwnb(i,1) + it(j)%flwnb(k,1)
!          write(*,*) 'Pushing ',maxfl,' / ',exflow,' down ',i,' ( ',j,').'
          if ((curnode.eq.-1).AND.(j.ne.tnode).AND.(j.ne.snode)) curnode = j
          if (j.ne.tnode) then
            if (cutmode.EQV..true.) then
              if (it(j)%ldis.lt.nbasins) it(j)%active = 1
            else
              it(j)%active = 1
            end if
          end if
          if (maxfl.ge.exflow) then
            it(j)%flwnb(k,2) = it(j)%flwnb(k,2) + exflow
            it(inode)%flwnb(i,1) = it(inode)%flwnb(i,1) + exflow
            exflow = 0
          else
            it(j)%flwnb(k,2) = it(j)%flwnb(k,2) + maxfl
            it(inode)%flwnb(i,1) = it(inode)%flwnb(i,1) + maxfl
            exflow = exflow - maxfl
            block(i) = 1
            avail = avail - 1
          end if
          if (exflow.eq.0) exit
        end if
      end do
    else
      if (exflow.eq.0) exit
      minldis = HUGE(minldis)
      do i=1,it(inode)%nb
        if (it(inode)%lstnb(i).eq.inode) cycle
        j = it(inode)%lstnb(i) 
        k = it(inode)%map(i) ! it(j)%map(inode)
        maxfl = it(inode)%wghtsnb(i,1) - it(inode)%flwnb(i,1) + it(j)%flwnb(k,1)
        if (maxfl.gt.0) then
!          write(*,*) 'eligible: ',j, '(',inode,').'
          if (it(j)%ldis.lt.minldis) minldis = it(j)%ldis
        end if
      end do
!       write(*,*) 'Relabeling ',inode,' to ',minldis + 1,'.'
      it(inode)%ldis = minldis + 1
      do i=1,it(inode)%nb
        if (it(inode)%lstnb(i).eq.inode) cycle
        j = it(inode)%lstnb(i)
        k = it(inode)%map(i) ! it(j)%map(inode)
        maxfl = it(inode)%wghtsnb(i,1) - it(inode)%flwnb(i,1) + it(j)%flwnb(k,1)
!        write(*,*) 'Max flow for ',j,' is ',maxfl,'.'
        if ((maxfl.gt.0).AND.(it(j)%ldis.eq.(it(inode)%ldis-1))) then
          block(i) = 0
          avail =  avail + 1
        end if
      end do
    end if
  end do
!
  it(inode)%active = 0
!
end
!
!-------------------------------------------------------------------------------------
!
subroutine process_pushresflow(it,inode,nbasins,setmap,setsz)
!  
  use clusters
  use iounit
!
  implicit none
!
  integer i,j,inode,nbasins,setmap(nbasins),setsz
  type(t_scluster) it(nbasins)
!
  do i=1,it(inode)%nb
    if (it(inode)%lstnb(i).eq.inode) cycle
    j = it(inode)%lstnb(i)
    if (it(inode)%wghtsnb(i,2).gt.(it(inode)%flwnb(i,2)-it(inode)%flwnb(i,1))) then
      it(j)%rflw = 1
      if (setmap(j).eq.0) then
        setsz = setsz + 1
        setmap(j) = 1
      end if
    end if
  end do
!
end
!
!--------------------------------------------------------------------------------------------
!
function getexflow(it,inode)
!
  use clusters
!
  implicit none
!
  integer getexflow,inode,i
  type(t_scluster) it(*)
!
  getexflow = 0
  do i=1,it(inode)%nb
    if (it(inode)%lstnb(i).eq.inode) cycle
    getexflow = getexflow + it(inode)%flwnb(i,2) - it(inode)%flwnb(i,1)
  end do
!
end
!
!------------------------------------------------------------------------------------------
!
subroutine tarjan_SCC(it,nbasins)
!
  use clusters
!
  implicit none
!
  integer nbasins,cursccid,glidx,i,idxcnt
  integer, ALLOCATABLE:: idxlst(:)
  type (t_scluster) it(nbasins)
!
  allocate(idxlst(nbasins+10))
!
  it(1:nbasins)%active = 0
  it(1:nbasins)%rflw = -1
  glidx = 0
  cursccid = 1
  it(1:nbasins)%inscc = 0
  idxcnt = 0
!
  do i=1,nbasins
    if (it(i)%inscc.gt.0) cycle
    call graph_dfs(it,nbasins,i,glidx,cursccid,idxlst,idxcnt)
  end do
!  
  deallocate(idxlst)
!
end
!
!-------------------------------------------------------------------------------------------
!
recursive subroutine graph_dfs(it,nbasins,which,glidx,cursccid,idxlst,idxcnt)
!
  use clusters
!
  implicit none
!
  integer glidx,i,nbasins,cursccid,which,idxcnt,idxlst(nbasins+10)
  type (t_scluster) it(nbasins)
!
  it(which)%ldis = glidx
  it(which)%active = glidx
  it(which)%rflw = 1
  glidx = glidx + 1
  idxcnt = idxcnt + 1
  idxlst(idxcnt) = which
  it(which)%inscc = cursccid
!
!  if (it(it(which)%lstnb(1))%nb.eq.1) then
!    write(*,*) it(which)%lstnb(1),' is a leaf'
!  else if ((it(it(which)%lstnb(1))%nb.eq.2).AND.&
! &    ((it(it(which)%lstnb(1))%wghtsnb(1,1).eq.0).AND.(it(it(which)%lstnb(1))%wghtsnb(2,2).eq.0)).OR.&
! &    ((it(it(which)%lstnb(1))%wghtsnb(1,2).eq.0).AND.(it(it(which)%lstnb(1))%wghtsnb(2,1).eq.0))) then
!    write(*,*) it(which)%lstnb(1),' is on a loop'
!  end if
  do i=1,it(which)%nb
    if (it(which)%lstnb(i).eq.which) cycle
    if (it(which)%wghtsnb(i,1).gt.0) then
      if (it(it(which)%lstnb(i))%active.le.0) then
        call graph_dfs(it,nbasins,it(which)%lstnb(i),glidx,cursccid,idxlst,idxcnt)
        it(which)%ldis = min(it(which)%ldis,it(it(which)%lstnb(i))%ldis)
      else if (it(it(which)%lstnb(i))%rflw.eq.1) then
        it(which)%ldis = min(it(which)%ldis,it(it(which)%lstnb(i))%active)
      end if
    end if
  end do
!
  if (it(which)%ldis.eq.it(which)%active) then
    do i=idxcnt,1,-1
      it(idxlst(i))%inscc = cursccid
      it(idxlst(i))%rflw = -1
      if (idxlst(i).eq.which) then
        idxcnt = i-1
        exit
      end if
    end do
    cursccid = cursccid + 1
  end if
!
end
!
!-------------------------------------------------------------------------------------------
!
subroutine mergetrees(ntrees,a,b,ib,snap2tree,lsnap2tree)
!
  use iounit
  use clusters
!
  implicit none
!
  integer ntrees
  type(t_progindextree) a,b ! component a is added to component b
  integer ib ! index of component b
  integer lsnap2tree
  integer snap2tree(lsnap2tree)
  integer, ALLOCATABLE:: snaps2(:)
  integer i ! local snapshot index, loop variable
  integer nsnaps2 ! new number of snapshots
!
  nsnaps2 = b%nsnaps + a%nsnaps
  allocate(snaps2(nsnaps2))
!
  ntrees = ntrees - 1
! update information for component b
  do i=1,b%nsnaps
    snaps2(i) = b%snaps(i)
  end do
  do i=(b%nsnaps+1),nsnaps2
    snaps2(i) = a%snaps(i)
  end do
  b%nsnaps = nsnaps2
  deallocate(b%snaps)
  allocate(b%snaps(nsnaps2))
  do i=1,nsnaps2
    b%snaps(i) = snaps2(i)
  end do
  if (a%connectedcomp.ne.b%connectedcomp) then
    write(ilog,*) 'NICO: Bug! a%connectedcomp = ',a%connectedcomp,' b%connectedcomp = ',b%connectedcomp
    call fexit()
  end if
! UPDATE SNAPS2COMP
  do i=1,a%nsnaps
    snap2tree(a%snaps(i)) = ib
  end do
! indicate component a as empty
  a%nsnaps = -1
  deallocate(a%snaps)
end
!
!-----------------------------------------------------------------------------------
!
subroutine heapify(heap,n,key,source,a)
!
  implicit none
!
  integer n,a,i,x,hleft,hright
  integer heap(n),source(n)
  real(KIND=4) key(n)
  logical isheap
!
  isheap = .false.
  i = a
  do while(isheap.EQV..false.)
    x = i
    if (hleft(i).le.n) then
      if (key(hleft(i)).lt.key(x)) then
        x = hleft(i)
      end if
    end if
    if (hright(i).le.n) then
      if (key(hright(i)).lt.key(x)) then
        x = hright(i)
      end if
    end if
    if (x.eq.i) then
      isheap = .true.
    else
      call hswap(heap,n,key,source,i,x)
      i = x
    end if
  end do
end
!
!-----------------------------------------------------------------------------------
!
subroutine hswap(heap,n,key,source,i,x)
!
  implicit none
!
  integer n,i,x
  integer heap(n),source(n)
  real(KIND=4) key(n)
  integer temp,tempsource
  real(KIND=4) tempkey
!
  temp = heap(i)
  tempkey = key(i)
  tempsource = source(i)
  heap(i) = heap(x)
  source(i) = source(x)
  key(i) = key(x)
  heap(x) = temp
  key(x) = tempkey
  source(x) = tempsource
end
!
!-----------------------------------------------------------------------------------
!
subroutine hbuild(heap,n,key,source)
!
  implicit none
!
  integer n,hparent,a
  integer heap(n),source(n)
  real(KIND=4) key(n)
!
  do a=hparent(n),1,-1
    call heapify(heap,n,key,source,a)
  end do
end
!
!-----------------------------------------------------------------------------------
!
function hparent(i)
!
  implicit none
!
  integer hparent,i
!
  if (i.eq.1) then
!    write(ilog,*) 'NICO: Bug! Parent of root of heap requested.'
!    call fexit()
  end if
  hparent = i/2
!
end
!
!-----------------------------------------------------------------------------------
!
function hleft(i)
!
  implicit none
!
  integer hleft,i
!
  hleft = 2*i
end
!
!-----------------------------------------------------------------------------------
!
function hright(i)
!
  implicit none
!
  integer hright,i
!
  hright = 2*i+1
end
!
!-----------------------------------------------------------------------------------
!
subroutine hprint(heap,n,key,source)
!
  use iounit
!
  implicit none
!
  integer n,i,hparent,hleft,hright
  integer heap(n),source(n)
  real(KIND=4) key(n)
!
  do i=1,n
    write(ilog,*) heap(i),key(i),source(i),hparent(i),hleft(i),hright(i)
  end do
end
!
!-----------------------------------------------------------------------------------
!
subroutine hinsert(heap,n,key,source,new,newkey,newsource)
!
  use clusters
  use iounit
!
  implicit none
!
  integer n,new,hparent,i,newsource
  integer heap(n+1),source(n+1)
  real(KIND=4) key(n+1)
  real(KIND=4) newkey
!
  if (n.ge.cstored) then
    write(ilog,*) 'NICO: Bug! Max heap size reached.'
    call fexit()
  end if
  n = n+1
  heap(n) = new
  key(n) = newkey
  source(n) = newsource
  i = n
  do while (i.gt.1)
    if (key(i).lt.key(hparent(i))) then
      call hswap(heap,n,key,source,i,hparent(i))
      i = hparent(i)
    else
      exit
    end if
  end do
end
!
!-----------------------------------------------------------------------------------
!
subroutine hremovemin(heap,n,key,source)
!
  implicit none
!
  integer n
  integer heap(n),source(n)
  real(KIND=4) key(n)
!
  heap(1) = heap(n)
  key(1) = key(n)
  source(1) = source(n)
  n = n-1
  call heapify(heap,n,key,source,1)
end
!
!-----------------------------------------------------------------------------------
!
subroutine pidxtree_growsiblings(t1)
!
  use clusters
!
  type(t_progindextree) t1
  integer tmp(max(t1%nsibalsz,1)),k
!
  if (t1%nsibalsz.eq.0) then
    t1%nsibalsz = 5
    allocate(t1%siblings(t1%nsibalsz))
  else
    tmp(:) = t1%siblings(1:t1%nsibalsz)
    deallocate(t1%siblings)
    k = t1%nsibalsz
    t1%nsibalsz = t1%nsibalsz*2
    allocate(t1%siblings(t1%nsibalsz))
    t1%siblings(1:k) = tmp(:)
  end if
!
end
!
!-----------------------------------------------------------------------------------------
!
  
