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
!
module zmatrix
!
  type t_rotlist
    integer alsz,alsz2 ! allocation size (first D of rotis)
    integer idx ! pointer to place in dc_di%recurs structure
    integer, ALLOCATABLE:: treevs(:) ! characterizers of local tree structure
    integer, ALLOCATABLE:: rotis(:,:) ! list of atom stretches
    integer, ALLOCATABLE:: diffis(:,:) ! list of stretches for difference set
  end type t_rotlist
!
  type(t_rotlist), ALLOCATABLE:: izrot(:)!,izrot_save(:)
!
  RTYPE, ALLOCATABLE:: blenpr(:),blen(:),bang(:),bangpr(:),ztor(:),ztorpr(:),ztor_save(:),blen_save(:),bang_save(:)
! #tyler pka
   RTYPE, ALLOCATABLE:: blen_limits(:,:), bang_limits(:,:)
   integer limits_cur
  integer nadd,dihed_wrncnt,dihed_wrnlmt,dblba_wrncnt,dblba_wrnlmt
  integer, ALLOCATABLE:: iz(:,:),iadd(:,:)
!
end module zmatrix