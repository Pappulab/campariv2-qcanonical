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
! esterms*     : energy terms: 1: IPP, 2: CORR, 3: attLJ, 4: IMPSOLV (DMFI), 5: WCA, 6: POLAR
!                7: TOR, 8: ZSEC, 9: TABUL, 10: DREST, 11: LCTOR (not in use), 12: boundary,
!                13: global SEV (not in use), 14: POLY, 15: BOND(1) = BL, 16: BOND(2) = BA,
!                17: BOND(3) = IM, 18: BOND(4) = DI, 19: DSSP, 20: CMAP
! nindx        : inverse power-potential exponent (repulsive wall)
! nhalf        : half of it
! is_*         : classifiers for usign dedicated branches in global force loops
! scale_*      : global outside scale factors for energy terms
! use_*        : global logicals for energy terms
!
module energies
!
! some remaining fix-size parameters (MAXTORLC(BINS)) should 
! be removed, but the feature is not supported currently anyway
!
  integer MAXTORLC,MAXTORLCBINS
  integer MAXENERGYTERMS,MAXENERGYPARAMS
  parameter (MAXTORLC=100)
  parameter (MAXTORLCBINS=20)
  parameter (MAXENERGYTERMS=25)
  parameter (MAXENERGYPARAMS=20)
!
! GENERAL VARIABLES
!
! global energy variables
  RTYPE esave,esavec,esaver,ekin
  RTYPE esterms(MAXENERGYTERMS)
  RTYPE esterms_tr(MAXENERGYTERMS),esterms_lr(MAXENERGYTERMS)
!
! simplifier flags
  logical is_plj,is_tablj,is_tab,is_lj,is_ev,ideal_run,is_impljp,is_implj
  logical is_fegev,is_fegplj,is_feglj,is_pewlj,is_prflj,is_fegprflj,do_n2loop
!
! flags and scaling factors for energy terms
  RTYPE scale_attLJ,scale_IPP,scale_CORR,scale_IMPSOLV,scale_WCA
  RTYPE scale_POLAR,scale_TOR,scale_ZSEC,scale_TABUL,scale_LCTOR
  RTYPE scale_POLY,scale_DREST,scale_FEGS(MAXENERGYTERMS),scale_DSSP
  RTYPE scale_BOND(5),scale_EMICRO
  logical use_IPP,use_attLJ,use_CORR,use_IMPSOLV,use_WCA,use_POLAR
  logical use_DREST,use_TOR,use_ZSEC,use_TABUL,use_LCTOR,use_POLY
  logical use_FEGS(MAXENERGYTERMS),use_FEG,use_BOND(5),use_DSSP,use_EMICRO
! says if I am using a required ZSEC length #tyler
  logical use_ZSEC_MinL(2)
! list of the required length for zsec #tyler
  integer par_ZSEC_L(2)

! #tyler pka
  logical do_pka
  RTYPE par_pka(10)
  
  ! Martin 
  integer par_top
  logical fos_corr,dpgrp_corr
  
  RTYPE save_E_HS(2,2)!,hsq_prog(3)

  RTYPE, allocatable :: par_hsq(:,:),par_hsq_save(:,:),hsq_sign(:)
  logical do_hsq,do_hs
  integer hsq_mode, TOT_QX
  RTYPE ener_bath
  character(MAXSTRLEN) charge_state_file
  logical read_charge_state
  integer n_charge_states
  integer n_res,n_titrable
  integer,allocatable :: q_states(:,:)
  
  integer,allocatable :: hs_residues(:)

  integer,allocatable :: titrable_res_index(:)
  
  logical E_rep_offset,FEG_FOS_OFFSET
  RTYPE E_rep_offset_val,FEG_FOS_OFFSET_VAL,FEG_FOS_OFFSET_FACT
  RTYPE BATH_ENER,DEBUG_ENER_SAVE,DEBUG_ENER_DRIFT
  RTYPE sc_alpha
  RTYPE sc_n
  RTYPE save_ener(2,2)
  integer hs_mode
  logical debug_pka
!  RTYPE, ALLOCATABLE:: ener_lambda_deriv(:)! This array will contain the memory of the derivative of the energy with regard to lambda
! Only problem is that I may need to calculate this derivative every step, even if I don't use it 
! parameters for energy terms
  logical use_hardsphere
  integer torlcmode,scrq_model,i_sqm,polbiasmode,nindx,nhalf
  integer fegljmode,fegcbmode,fegmode,rf_mode
  integer par_LCTOR2(MAXENERGYPARAMS)
  integer, ALLOCATABLE:: par_TOR2(:),par_POLY2(:)
  RTYPE primamide_cis_chgshft,dpgrp_neut_tol
  RTYPE par_IMPSOLV(MAXENERGYPARAMS),par_WCA(MAXENERGYPARAMS)
  RTYPE, ALLOCATABLE:: par_TOR(:,:),par_POLY(:,:)
  RTYPE par_ZSEC(MAXENERGYPARAMS),par_ZSEC2(MAXENERGYPARAMS)
  RTYPE coul_scr,lct_val(MAXTORLCBINS),par_LCTOR(MAXENERGYPARAMS)
  RTYPE lct_pot(MAXTORLC,MAXTORLCBINS),par_FEG2(MAXENERGYPARAMS)
  RTYPE lct_bnd(MAXTORLC,2),par_POLAR(MAXENERGYPARAMS),rfcnst
  character(MAXSTRLEN) polfile,torfile,torlcfile2,lctpotfile,fegfile
  logical tor_report,poly_report,feg_report,use_FEGS_FOS
  logical, ALLOCATABLE:: par_FEG(:),par_FEG3(:)
  RTYPE HS_attlj,HS_replj,HS_fos,HS_polar

! solvent-accessible volume and SA analysis (these are
! also used for energy terms)
  integer savcalc,nsavavg
  RTYPE savprobe,savavg,sexv,sysvol(3)
  RTYPE, ALLOCATABLE:: atsavavg(:,:)
  integer, ALLOCATABLE:: natsavavg(:)
!
end module energies