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
! -- GENERAL
! ewald_mode   : which Ewald treatment to use (1: PME, 2: Standard)
! ewpm         : the Ewald parameter determining weights of real vs. reciprocal space in 1/A
! ewpm_pre     : a temporary variable to allow the user to overwrite ewpm
! ewetol,ewftol: energy and force tolerances for Ewald, respectively
! ewfspac      : spacing in units of length (grid size in PME, prop. to inverse reciprocal space cutoff in Ewald)
! kdims        : lower and upper bounds on grid/image index (symm. around zero), directly derived from ewfspac
! ewcnst       : variable to hold the constant self-correction energy in the Ewald sum
! ewpm2        : ewpm*ewpm
! ewpite       : 2.0*ewpm/sqrt(PI)
! ewinvm       : inverse of total number of (partial) charges
! -- PME
! splor        : order of the B-spline interpolation in PME
! bspl,bspld   : B-splines and derivatives per partial charge
! bsplbu       : a temporary array used to hold B-splines for an individual charge and dimension (and thread)
! bsmx/y/z     : B-spline moduli along the first, second and third dimensions of the reciprocal lattice, respectively
! ewgfls       : grid coordinates per partial charge
! ewnetf       : 3D total force for non-conservation compensation in PME
! Qew1,Qew2    : arrays to hold splined charges and transforms / convolutions
! QewBC        : array to hold static grid variables in constant volume sim.s
! fftplanf     : FFTW plan (holding pointers Qew2->Qew1) which does the second transformation
! fftplanb     : FFTW plan (holding pointers Qew1->Qew2) which does the first transformation
!
module ewalds
!
  RTYPE ewpm,ewcnst,ewpm2,ewpite,ewetol,ewftol,ewfspac,ewpm_pre,ewinvm,ewnetf(3)
  RTYPE, ALLOCATABLE:: bsmx(:),bsmy(:),bsmz(:)
  RTYPE, ALLOCATABLE:: bspld(:,:,:),bsplbu(:,:),bspl(:,:,:)
  integer, ALLOCATABLE:: ewgfls(:,:)
  integer kdims(3,2),splor,ewald_mode
  integer(8) fftplanf,fftplanb
  complex(KIND=8), ALLOCATABLE:: Qew1(:,:,:),Qew2(:,:,:)
  RTYPE, ALLOCATABLE:: QewBC(:,:,:)
!
end module ewalds
!
