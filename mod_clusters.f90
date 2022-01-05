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
! cludata            : data extracted for clustering (all read into memory -> large)
! cstored            : number of snapshots stored for clustering (-> cludata)
! cmode              : choice of algorithm for clustering or similar tasks (1-5)
! cdis_crit          : the distance criterion underlying the clustering (1-8)
! cstorecalc         : frequency for collecting data for clustering
! cmaxsnaps          : allocation size for cludata in numbers of snapshots dimension
! calcsz,clstsz      : dimensionality size variables for stored data (-> cludata)
! cfilen             : input file for defining clustering subset
! cdofset/cl_imvec   : auxiliary variables for clustering data
! cradius            : general size threshold for clustering
! cprepmode          : options for signal (data) preprocessing
! cchangeweights     : for eligible cdis_crit, replace weights for distance evaluations
! cwcombination      : for locally adaptive weights, how to combine (L_? norm)
! cwwindowsz         : for locally adaptive weights, how big a window to use
! cdynwbuf           : for LAWs dependent on transition counts, buffer parameter for taking inverse
! cwacftau           : for ACF-based weights, fixed lag time
! csmoothie          : order for smoothing fxn (sets window size)
! refine_clustering  : general flag to trigger refinement in clustering
! birchtree,scluster : arrays to hold clustering result
! sumssz             : global size variable for (t_scluster)%sums(:,:)
! cleadermode        : if cmode = 1/2/5(4): processing direction flags
! clinkage           : if cmode = 3: linkage criterion (1-3)
! chardcut           : if cmode = 3(4): truncation cut for snapshot NB-list
! nblfilen,read_nbl_from_nc,cnblst,cnc_ids : if cmode = 3(4): snapshot NB-list variables
! cprogindex         : if cmode = 4: exact or approximate scheme (1-2)
! cprogindrmax       : if cmode = 4: number of nearest neighbor guesses for approximate scheme
! cprogindstart      : if cmode = 4: snapshot to start from
! cprogpwidth        : if cmode = 4: the (maximum) size of A/B partitions in localized cut
! cprogfold          : if cmode = 4: how often to collapse terminal vertices of the MST
! csivmax,csivmin    : if cmode = 4: used to identify starting points automatically
! c_nhier            : if cmode = 5(4): tree height
! cmaxrad            : if cmode = 5(4): root level size threshold
! ccfepmode          : if any clustering is performed, type of cFEP to produce
! align_for_clustering: if cdis_crit = 5/6 whether to do alignment
! cdofsetbnds        : if cdis_crit = 6, selectors for separating sets for alignment and distance compu.
! pcamode            : whether to do PCA and what output to produce (1-3)
! reduced_dim_clustering: whether to proceed with clustering in reduced dimensional space after PCA
! align              : structure used for trajectory alignment (ALIGNFILE)
!
module clusters
!
  type t_scluster
    integer nmbrs,center,geni,genidx,alsz,nb,ldis,active,rflw,nbalsz,chalsz,nchildren,parent,inscc ! misc.
    RTYPE radius,diam,nodewt,quality,sqsum ! properties
    RTYPE, ALLOCATABLE:: sums(:,:),lensnb(:,:) ! centroid data and mean geometric length associated with edges
    integer, ALLOCATABLE:: snaps(:),tmpsnaps(:) ! snapshot lists
    integer, ALLOCATABLE:: map(:),children(:),wghtsnb(:,:),lstnb(:),flwnb(:,:) ! for graph
  end type t_scluster
!
  type t_ctree
    type(t_scluster), ALLOCATABLE:: cls(:) ! nothing but an array of clusters
    integer ncls,nclsalsz ! real and allocation sizes
  end type t_ctree
!
  type t_cnblst
    integer nbs,alsz ! actual number of neighbors and current allocation size
    RTYPE, ALLOCATABLE:: dis(:) ! list of distances
    integer, ALLOCATABLE:: idx(:) ! and associated indices
    logical, ALLOCATABLE:: tagged(:) ! helper flag
  end type t_cnblst
!
  type t_progindextree
    integer nsnaps ! number of snapshots in tree
    integer nsiblings ! number of trees to merge with
    integer nsibalsz ! alloc size for that
    integer, ALLOCATABLE:: snaps(:) ! indices of snapshots in tree
    integer, ALLOCATABLE:: siblings(:) ! tree indices of tree to merge with
    integer mine(2) ! shortest edge leaving the tree
    RTYPE mind ! distance to nearest snapshot of tree, equals length(mine)
    integer connectedcomp ! giving the index of the connected component to which this tree will belong
  end type t_progindextree
!
  type t_progindexcomponent
    integer ntrees ! number of trees in connected component
    integer, ALLOCATABLE:: trees(:) ! indices of trees in connected component
  end type t_progindexcomponent
!
  type t_adjlist
    integer deg ! degree of vertex
    integer alsz ! allocation size
    integer, ALLOCATABLE:: adj(:) ! list of adjacent vertices
    real(KIND=4), ALLOCATABLE:: dist(:) ! distance to the adjacent vertices
  end type t_adjlist
!
  type(t_cnblst), ALLOCATABLE:: cnblst(:)
  type(t_scluster), ALLOCATABLE:: scluster(:)
  type(t_ctree), ALLOCATABLE:: birchtree(:)
  type(t_adjlist), ALLOCATABLE:: approxmst(:)
  RTYPE, ALLOCATABLE:: cludata(:,:),cl_imvec(:)
  integer, ALLOCATABLE:: trbrkslst(:) ! a user-populated list of transitions to delete in graph-related analyses
  integer cstored,cdis_crit,cstorecalc,calcsz,clstsz,sumssz,cmode,nstruccls,pcamode,cprogindex,cprogindrmax,cdofsbnds(4),cmaxsnaps
  integer csivmin,csivmax,cnc_ids(10),c_nhier,clinkage,cleadermode,reduced_dim_clustering,cprogindstart,cprogpwidth,ccfepmode
  integer ntbrks,ntbrks2,itbrklst,cequil,cchangeweights,cwacftau,cwcombination,cwwindowsz,csmoothie,cprepmode,cprogfold
  integer, ALLOCATABLE:: cdofset(:,:)
  RTYPE cradius,cmaxrad,chardcut,cdynwbuf
  character(MAXSTRLEN) cfilen,nblfilen,tbrkfilen
  logical read_nbl_from_nc,align_for_clustering,refine_clustering
!
! other
  type t_align
    integer, ALLOCATABLE:: set(:) ! atom set
    RTYPE, ALLOCATABLE:: refxyz(:),curxyz(:) ! coordinate helpers
    integer nr,calc,mmol ! size, frequency, ref. mol.
    logical yes,refset,haveset,instrmsd
    character(MAXSTRLEN) filen ! for providing custom alignment set
  end type t_align
!
  type(t_align) align
!
end module clusters
!
