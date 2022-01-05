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
!-----------------------------------------------------------------------
!
!
!             ##################################
!             #                                #
!             #  A SET OF SUBROUTINES WHICH    #
!             #  DEAL WITH PROGRAM FLOW AND    #
!             #  FILE INTERFACES               #
!             #                                #
!             ##################################
!
!
!-----------------------------------------------------------------------------
!
! fatal exit: call when encountering an internal problem
!
subroutine fexit()
!
  use iounit
!
  implicit none
!
  integer k
!
#ifdef ENABLE_MPI
#include "mpif.h"
#endif
!
!
#ifdef ENABLE_MPI
  call makelogio(2)
  call MPI_Abort(MPI_COMM_WORLD,1)
#endif
!
  k = std_f_out
  write(k,*)
  write(k,*) '------------------------------------------------>'
  write(k,*)
  write(k,*) 'CAMPARI CRASHED UNEXPECTEDLY. PLEASE RECORD ANY '
  write(k,*) ' INFORMATION ABOUT THE PROBLEM PROVIDED ABOVE.'
  write(k,*)
  write(k,*) '<------------------------------------------------'
!
  stop 98
!
end
!
!-----------------------------------------------------------------------
!
! a trivial fxn to obtain command line arguments
!
subroutine grab_args()
!
  use commline
  use iounit
!
  implicit none
!
  integer argstat,i,k
!
! the FORTRAN intrinsics can handle this well in modern Fortran
! note the zero-index argument is the program call itself and is
! not counted in COMMAND_ARGUMENT_COUNT
!
  nargs = COMMAND_ARGUMENT_COUNT()
!
  if (nargs.eq.0) then
    k = std_f_out
    write(k,*) 'USAGE: CAMPARI requires a mandatory key-file with input options &
 &that specify the calculation to be attempted.'
    write(k,*)
    write(k,*) 'EXAMPLE: ${PATH_TO_CAMPARI}/bin/${ARCH}/campari -k ${PATH_TO_CAMPARI}/examples/tutorial1/TEMPLATE.key'
    call fexit()
  end if
!
  allocate(args(nargs))     
!
  do i=1,nargs
    call GET_COMMAND_ARGUMENT(i,args(i)%it,args(i)%leng,argstat)
  end do
     
!
end
!
!-----------------------------------------------------------------------
!
! a function which provides an empty file handle
! it just loops up until it finds an open one (tests explicitly)
!
function freeunit()
!
  use iounit
!
  implicit none
!
  integer freeunit
  logical inuse
!
! try each logical unit until an unopened one is found
!
  freeunit = 0
  inuse = .true.
  do while (inuse)
    freeunit = freeunit + 1
    if ((freeunit.ne.std_f_in).AND.(freeunit.ne.std_f_out)) then
      inquire (unit=freeunit,opened=inuse)
    end if
!
    if (freeunit.eq.huge(freeunit)) then
      write(ilog,*) 'Fatal. No free I/O unit could be found. This mo&
 &st likely reports on a more severe underlying problem with the fil&
 &e system.'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! this fxn saves about one line of code and checks whether a given
! unit is associated with an open file
!
function is_filehandle_open(fh)
!
  implicit none
!  
  integer fh
  logical is_filehandle_open
!
  inquire (unit=fh,opened=is_filehandle_open)
!
  return
!
end function is_filehandle_open
!
!-----------------------------------------------------------------------
!
! this routine is a shortcut for deleting a file, and reopening it 
! to the handle provided
! it is important that the passed string is properly stripped of blanks on the calling side
!
subroutine delete_then_openfile(fh,fpath)
!  
  implicit none
!  
  integer fh,freeunit
  character(len=*) fpath
  logical exists
!
  inquire(file=fpath,exist=exists)
  if (exists.EQV..true.) then
    fh = freeunit()
    open (unit=fh,file=fpath,status='old')
    close(unit=fh,status='delete')
  end if
!
  fh = freeunit()
  open (unit=fh,file=fpath,status='new') 
!
end subroutine delete_then_openfile
!
!-----------------------------------------------------------------------
!
! this routine is a shortcut for opening an existing file 
! it throws a warning upon failure to find such a file, but does not
! open a new file
! it is important that the passed string is properly stripped of blanks on the calling side
!
function openfile(fh,fpath)
!
  implicit none
!
  integer fh, freeunit
  character(len=*) fpath
  logical exists,openfile
!
  inquire(file=fpath,exist=exists)
  if (exists.EQV..true.) then
    fh = freeunit()
    open (unit=fh,file=fpath,status='old')
    openfile = .true.
  else
    write(*,*) 'WARNING: Cannot open input file (',fpath,') in openfile(...).'
    openfile = .false.
  end if
!
  return
!  
end function openfile
!
!-----------------------------------------------------------------------
!
! this is a redundant wrapper that saves zero lines of code
!
subroutine close_filehandle(fh)
!    
  implicit none
!
  integer fh
  logical is_filehandle_open
!
  if(is_filehandle_open(fh))  close(unit=fh)
!
end subroutine close_filehandle
!
!-----------------------------------------------------------------------
!
! further functions may be added
!
!   To be added. Common file manipulation support is lacking in fortran
!   subroutine copy_file(ffrom, fto)
!   
!     implicit none
! 
! !     character(len=*) :: fpath
! !     integer isrc, idst
! ! 
! !     OPEN(UNIT=ISRC, FILE='', ACCESS='DIRECT', STATUS='OLD', ACTION='READ', IOSTAT=IERR, RECL=1)
! ! OPEN(UNIT=IDST, FILE='', ACCESS='DIRECT', STATUS='REPLACE', ACTION='WRITE', IOSTATE=IERR, RE)
! ! IREC = 1
! ! DO
! !   READ(UNIT=ISRC, REC=IREC, IOSTAT=IERR) CHAR
! !   IF (IERR.NE.0) EXIT
! !   WRITE(UNIT=IDST, REC=I) CHAR
! !   IREC = IREC + 1
! ! END DO
! ! 
! ! 
!   end subroutine copy_file
!
