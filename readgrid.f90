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
! MAIN AUTHOR:   Hoang Tran                                                !
! CONTRIBUTIONS: Rohit Pappu, Andreas Vitalis                              !
!                                                                          !
!--------------------------------------------------------------------------!
!
!
#include "macros.i"
!
! ###################################################
! ##                                               ##
! ##  subroutine readgrid -- read in phi-psi grids ##
! ##                                               ##
! ###################################################
!
! "readgrid" reads in a predetermined grid of important phi,psi
! values into memory for later use in Markov chain sampling
!
subroutine readgrid()
!
  use grids
  use iounit
  use aminos
!
  implicit none
!
  integer funit,freeunit,next,i,t1,t2
  RTYPE ff,yy
  character(60) fyfile
  character(80) record
  character(3) aa3lc
  logical exists
!
  call strlims(griddir,t1,t2)
  write(ilog,*)
  write(ilog,*) 'Using grids in ',griddir(t1:t2)
!
  do i = 1,MAXAMINO
     aa3lc = amino(i)
     call tolower(aa3lc)
     fyfile=griddir(t1:t2)//aa3lc//'_grid.dat'
     inquire(file=fyfile,exist=exists)
     if (exists.EQV..false.) then
!
!       read in default grid
        fyfile=griddir(t1:t2)//'def_grid.dat'
     end if
     next = 1
!    
! read in grid centers from fyfile
     inquire(file=fyfile,exist=exists)
     if(exists.EQV..false.) then 
        call fexit()
     else
        funit=freeunit()
        open(unit=funit,file=fyfile,status='old')
        do while (.true.)
           read(funit,10,end=30) record
 10            format(a80)
           read(record,*,err=20,end=20)ff,yy
 20            continue
           stgr%it(i,next,1) = ff
           stgr%it(i,next,2) = yy
           next = next + 1
        end do
 30         continue
        next = next - 1 
        stgr%ngr(i) = next
        close(unit=funit)
     end if
      
     write(ilog,35) stgr%ngr(i),aa3lc//'_grid.dat',amino(i)
 35      format('     Read in [',i5,'] grid centers from: ',a20,&
 &          ' for ',a3)
  end do
!
  stgr%sfyc(1) = stgr%halfwindow
  write(ilog,*)
  do i = 1,MAXAMINO
     stgr%sfyc(i) = stgr%sfyc(1)
  end do
!
end
!
!-----------------------------------------------------------------------
!

