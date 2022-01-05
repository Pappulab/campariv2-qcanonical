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
! use_cutoffs     : whether cutoffs are in use at all
! use_mcgrid      : whether grid-based cutoffs are in use (assignment through reference atom position)
! use_rescrit     : whether topology-based cutoffs are in use (usage of the known residue topology to pre-screen distances via
!                   reference atom distances)
! use_stericscreen: whether to skip out of a MC energy evaluation during an elementary step if steric clash occurs based on
!                   short-range energy for new conformation
! screenbarrier   : what quantifies a steric clash energetically
! use_waterloops  : whether to take advantage of special loops written for SPC/T3P and for T4P in dynamics
! rsw1            : the residue number for the first water molecule (waters must be at end of sequence file)
! mcnbcutoff      : short-range cutoff used in MC for the following potentials: WCA, IPP, attLJ, as well as for solvation states
!                   (via sav); in MD/LD short-range cutoff determines interactions re-calculated at every step regardless of what
!                   potentials are used
!                   global cutoff checks are always consistent but interior cutoff checks might not be (i.e., in vectorized inner
!                   loops distance is not re-rested, but in non-vectorized loops it is, and hence more interactions are EXcluded)
! mcnbcutoff2     : mcnbcutoff squared
! mcelcutoff      : long-range cutoff used in MC for the following potentials: POLAR, TABUL
!                   larger of the twin-range cutoffs LD/MD, i.e., interactions between mcnbcutoff and mcelcutoff are only
!                   computed every so often (nbl_up)
!                   there is never an interior cutoff check within the inner loops
! mcelcutoff2     : mcelcutoff squared
! imcel2          : inverse of mcelcutoff2
! nbl_up          : interval as to how often to re-compute the neighbor lists for MD/LD, also interval for re-computing the
!                   twin-range interactions between mcnbcutoff and mcelcutoff
! lr_up           : currently not in use
! lrel_md         : long-range electrostatics treatment flag for MD/LD (see documentation): 1) cut, 2) Ewald, 3) (G)RF,
!                   4) all monopole
! lrel_mc         : long-range electrostatics treatment flag for MC (see doc.): 1) all monopole, 2) just monopole-monopole,
!                   3) reduced
! cglst%ncs       : number of monopole groups in system
! cglst%it        : pointer array for all monopole groups in system to central atom (closest to monopole center)
! cglst%nc        : value of monopole charge for all monopole groups 
! rs_nbl          : neighbor list for every residue
! rs_nbl%nnbs     : number of fully coupled partner residues within mcnbcutoff
! rs_nbl%nnbats   : total number of atoms in fully coupled partner residues within mcnbcutoff
! rs_nbl%nb       : vector of fully coupled partner residues within mcnbcutoff
! rs_nbl%svec     : vector of shift vectors for fully coupled partner residues within mcnbcutoff (PBC)
! rs_nbl%nnbtrs   : number of fully coupled partner residues within mcelcutoff but outside of mcnbcutoff
! rs_nbl%nnbtrats : total number of atoms in fully coupled partner residues within mcelcutoff but outside of mcnbcutoff
! rs_nbl%nnblrs   : number of fully coupled partner residues outside of mcelcutoff
! rs_nbl%nnblrats : total number of atoms in fully coupled partner residues outside of mcelcutoff
! rs_nbl%nbtr     : vector of fully coupled partner residues within mcelcutoff but outside of mcnbcutoff
! rs_nbl%trsvec   : vector of shift vectors for fully coupled partner residues within mcelcutoff but outside of mcnbcutoff
! rs_nbl%nblr     : vector of fully coupled partner residues outside of mcelcutoff
! rs_nbl%lrsvec   : vector of shift vectors for fully coupled partner residues outside of mcelcutoff
! rs_nbl%ngnbs    : number of FEG-coupled partner residues within mcnbcutoff
! rs_nbl%ngnbats  : total number of atoms in FEG-coupled partner residues within mcnbcutoff
! rs_nbl%gnb      : vector of FEG-coupled partner residues within mcnbcutoff
! rs_nbl%gsvec    : vector of shift vectors for FEG-coupled partner residues within mcnbcutoff
! rs_nbl%ngnbtrs  : number of FEG-coupled partner residues within mcelcutoff but outside of mcnbcutoff
! rs_nbl%ngnbtrats: total number of atoms in FEG-coupled partner residues within mcelcutoff but outside of mcnbcutoff
! rs_nbl%ngnbtrs  : number of FEG-coupled partner residues outside of mcelcutoff
! rs_nbl%ngnbtrats: total number of atoms in FEG-coupled partner residues outside of mcelcutoff
! rs_nbl%gnbtr    : vector of FEG-coupled partner residues within mcelcutoff but outside of mcnbcutoff
! rs_nbl%gnblr    : vector of FEG-coupled partner residues outside of mcelcutoff
! rs_nbl%gtrsvec  : vector of shift vectors for FEG-coupled partner residues within mcelcutoff but outside of mcnbcutoff
! rs_nbl%glrsvec  : vector of shift vectors for FEG-coupled partner residues outside of mcelcutoff
! rs_nbl%ntabnbs  : total number of interacting partners through tabulated potential
! rs_nbl%ntabias  : total number of implied atom-atom interactions
! rs_nbl%tabnb    : vector of interacting partners through tabulated potential
!
module cutoffs
!
! neighbor lists
  type t_rs_nbl
    integer, ALLOCATABLE:: nb(:),nbtr(:),nblr(:),tabnb(:),tmpanb(:)
    RTYPE, ALLOCATABLE:: svec(:,:),trsvec(:,:),lrsvec(:,:)
    integer, ALLOCATABLE:: gnb(:),gnbtr(:),gnblr(:)
    RTYPE, ALLOCATABLE:: gsvec(:,:),gtrsvec(:,:),glrsvec(:,:)
    integer nnbs,nnbtrs,nnbats,nnbtrats,nnblrs,nnblrats,ntabnbs,ntabias
    integer ngnbs,ngnbtrs,ngnbats,ngnbtrats,ngnblrs,ngnblrats,ntmpanb,tmpalsz
    integer nbalsz,tralsz,lralsz,gnbalsz,gtralsz,glralsz
  end type t_rs_nbl
!
! temporary nb-list arrays 
  type t_tmp_nbl
    RTYPE, ALLOCATABLE:: dvec(:,:),svec(:,:),d1(:),svl(:,:)
    integer, ALLOCATABLE:: k1(:,:),k2(:,:),nbl(:),idx(:)
  end type t_tmp_nbl
!
! full charge list
  type t_cglst
    integer, ALLOCATABLE:: it(:),nc(:)
    integer ncs
    RTYPE, ALLOCATABLE:: tc(:),tc_limits(:,:)
  end type t_cglst
  
!
! cutoff variables
  RTYPE mcnb_cutoff,mcnb_cutoff2,mcel_cutoff,mcel_cutoff2,imcel2
  RTYPE screenbarrier
  logical use_cutoffs,use_mcgrid,use_rescrit,use_stericscreen
  logical use_waterloops
  integer, ALLOCATABLE:: rsp_mat(:,:),rsp_vec(:),glob_rslst(:)
  integer nsancheck,rsw1,nbl_up,lr_up,lrel_mc,lrel_md,nbl_maxnbs
  type(t_rs_nbl), ALLOCATABLE:: rs_nbl(:)
  type(t_cglst):: cglst
  type(t_tmp_nbl), ALLOCATABLE:: tmp_nbl(:)
!
end module cutoffs
!
