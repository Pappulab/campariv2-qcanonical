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
! MAIN AUTHOR:   Adam Steffen                                              !
! CONTRIBUTIONS: Andreas Vitalis, Nicholas Lyle                            !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! ##############################################################
! ##                                                          ##
! ## subroutine minimize -- Minimization over entire system   ##
! ##                                                          ##
! ##############################################################
!    
!
subroutine minimize(maxiter, grmsmin, istepsize, minimode)

  use iounit
  use molecule
  use forces
  use energies
  use movesets
  use system
  use mcsums
  use fyoc
  use cutoffs
  use mini
  use zmatrix
  use math
  use atoms
!      
  implicit none
!      
  integer maxiter,iter,nrdof,imol,j,i,j1,j2,t,pmode,iter_bfgs,wrnlmt(2),wrncnt(2),iter_avg
  integer ndump,minimode,mstore,incr,bound,evals,redux_count,redux_max,upaccept,badaccept,avgsz
  RTYPE grmsmin,grms,stepsize,eprevious,istepsize
  RTYPE force3,t1,t2,scalefactor,sigma,b_bfgs,curve,line_grad,grms_ahead,hnot_bfgs
  RTYPE, ALLOCATABLE:: cur_gradient_vec(:),prev_gradient_vec(:),d_pos_vec(:),peak_gradient_vec(:)
  RTYPE, ALLOCATABLE:: a_bfgs(:),q_bfgs(:),r_bfgs(:),s_bfgs(:,:),y_bfgs(:,:),normcg(:),avgstepsz(:)
  RTYPE, ALLOCATABLE:: mshfs(:,:),prev_normcg(:),peak_ahead_vec(:)
  logical postanalyze,atrue
!     
10   format(i8,1x,i8,1x,20000(f18.10,1x))
!
! possible globals
  mstore = min(mini_mem,10) ! threshold for accepting uphill steps in BFGS (max 10, min memory)
  redux_max = 50  ! number of stepsize divisions tried before exiting (usually indicates uphill attempt)
  avgsz = 100
!
! initializations
  evals = 0
  nrdof = 0
  ndump = 0
  grms = 0.0
  redux_count = 0
  iter = 0
  iter_bfgs = 0
  iter_avg = 0
  curve = 0.0
  upaccept = 0
  badaccept = 0
  wrncnt(:) = 0
  wrnlmt(:) = 1
  postanalyze = .false.
  atrue = .true.
  eprevious = esave
  stepsize = 1.0
!  
! "vectorize" all degrees of freedom into a single vector with dimension "nrdof"
  if (fycxyz.eq.1) then
    do imol=1,nmol
      do j=1,size(dc_di(imol)%frz)
        if (dc_di(imol)%frz(j).EQV..false.) then
          nrdof = nrdof + 1
        end if
      end do
      call update_rigidm(imol)
    end do 
  else if (fycxyz.eq.2) then
! SHAKE-type constraints and drift removal must be handled separately (not mappable)
    do i=1,n
      do j=1,3
        if (cart_frz(i,j).EQV..true.) cycle
        nrdof = nrdof + 1
      end do
    end do
    do imol=1,nmol
      call update_rigidm(imol)
    end do
  else
    write(ilog,*) 'Fatal. Called minimize(...) with unsupported system representation. This is an omission bug.'
    call fexit()
  end if
  if (nrdof.le.0) then
    write(ilog,*) 'Fatal. No remaining degrees of freedom to minimize. Relieve constraints.'
    call fexit()
  end if
!  
! allocations used by all methods
  allocate(cur_gradient_vec(nrdof))
  allocate(prev_gradient_vec(nrdof))
  allocate(d_pos_vec(nrdof))
  allocate(prev_normcg(nrdof))
  allocate(peak_gradient_vec(nrdof))
  allocate(peak_ahead_vec(nrdof))
  allocate(mshfs(3,nmol))
  allocate(normcg(nrdof))
  allocate(avgstepsz(avgsz))
  d_pos_vec(:) = 0.0
  normcg(:) = 0.0
  prev_normcg(:) = 0.0
  cur_gradient_vec(:) = 0.0
  avgstepsz(:) = 0.0
!  
! quasi-newton BFGS allocations/inits
  if(minimode.eq.3) then
    allocate(s_bfgs(nrdof,mini_mem))
    allocate(y_bfgs(nrdof,mini_mem))
    allocate(a_bfgs(mini_mem))
    allocate(q_bfgs(nrdof))
    allocate(r_bfgs(nrdof))
    s_bfgs(:,:) = 0.0
    y_bfgs(:,:) = 0.0
    a_bfgs(:) = 0.0
    q_bfgs(:) = 0.0
    r_bfgs(:) = 0.0
  end if  
!     
! ---------- MASTER LOOP -------------
!
  do while (evals.le.maxiter)
!
!   count the number of tried energies (peak and reduce method)
    evals = evals + 1
!
!   zero out internal coordinate forces
    if (fycxyz.eq.1) then
      do imol=1,nmol
        dc_di(imol)%f(:) = 0.0
      end do
    end if ! cart_f is zeroed elsewhere 
!
!   calculate energies and gradient at test point
    call CPU_time(t1)
    esave = force3(esterms,esterms_tr,esterms_lr,atrue)
    ens%insR(10) = ens%insU
    ens%insU = esave
    call CPU_time(t2)
    time_energy = time_energy + t2 - t1
    if (fycxyz.eq.1) call cart2int_f(skip_frz)
    peak_gradient_vec(:) = 0.0
    peak_ahead_vec(:) = 0.0
    call get_gradient(peak_ahead_vec,grms_ahead,nrdof)
! 
    line_grad = DOT_PRODUCT(peak_ahead_vec,cur_gradient_vec)
!
!   if the energy is lower (or nearly lower) then move on to a new point
    if ((evals.eq.1).OR.&
 &     ((minimode.eq.1).AND.(esave.lt.eprevious)).OR.&
 &     ((minimode.eq.2).AND.(esave.lt.eprevious).AND.(line_grad.gt.0.0)).OR.&
 &     ((minimode.eq.3).AND.(esave.lt.(eprevious + mini_uphill)).AND.((line_grad.gt.0.0).OR.(iter_bfgs.eq.1)))) then
      iter = iter + 1
      iter_bfgs = iter_bfgs + 1
      if (iter_avg.gt.0) then
        avgstepsz(mod(iter_avg,avgsz)+1) = stepsize
      end if
      iter_avg = iter_avg + 1
!
      nstep = nstep + 1
!
      pmode = 0
      peak_gradient_vec(:) = peak_ahead_vec(:)
      grms = grms_ahead
!      write(ilog,10) evals,iter,stepsize,esave,grms,curve
      call CPU_time(t1)
      call mcstat(iter,ndump)
      call CPU_time(t2)
      time_analysis = time_analysis + t2 - t1
!
!     store gradient information
      prev_gradient_vec(:) = cur_gradient_vec(:)
      prev_normcg(:) = normcg(:)
      cur_gradient_vec(:) = peak_gradient_vec(:)
      eprevious = esave
!
!     BFGS: update gradient stored value for l-bfgs
      if (minimode.eq.3) then
        j = MOD(iter_bfgs,mini_mem)
        if (j.eq.0) then
          j = mini_mem
        end if
        y_bfgs(:,j) = cur_gradient_vec(:)
      end if
!  
!     gradient RMS convergence check
      if (grms.lt.grmsmin) then
        write(ilog,*) 'Termination: Gradient RMS convergence reached.'
        postanalyze = .true.
        exit
      end if
      if(mod(iter,nsancheck).eq.0) then
 2  format('Step: ',i8,' Current Energy: ',g14.6,' G-RMS: ',g14.6,' Avg Step Scale: ',g14.6)
       write(ilog,2) iter,esave,grms,(1.0/(1.0*min(iter_avg,avgsz)))*sum(avgstepsz(1:min(iter_avg,avgsz)))
      end if
!
!     increase the stepsize if sequential minimization moves are accepted, reset if last step was rejected
!     note that due to step-size reductions it's rather likely for the energies in subsequent steps to be floating-point
!     identical (.ge.)
      if (esave.ge.eprevious) then
        badaccept = badaccept + 1
        if (esave.gt.eprevious) upaccept = upaccept + 1
!        
!       warn against real uphill problems
        if ((upaccept.gt.mstore)) then
          wrncnt(2) = wrncnt(2) + 1
          if (wrncnt(2).eq.wrnlmt(2)) then
            write(ilog,*) 'WARNING: Minimization has exceeded the permissible number of uphill steps.&
 & This means that the current direction(s) was (were) erroneous.'
            write(ilog,*) 'This was warning #',wrncnt(2),' of this type not all of which may be displayed.'
            wrnlmt(2) = wrnlmt(2)*10
          end if
        end if
!
!       reset
        if (badaccept.gt.mstore) then
          upaccept = 0
          badaccept = 0
          cur_gradient_vec(:) = peak_ahead_vec(:)
          if (mini_mode.eq.2) then
            prev_gradient_vec(:) = cur_gradient_vec(:)
            stepsize = (1.0/(1.0*min(iter_avg,avgsz)))*sum(avgstepsz(1:min(iter_avg,avgsz)))
          end if
          if (mini_mode.eq.3) then
            y_bfgs(:,:) = 0.0
            y_bfgs(:,1) = cur_gradient_vec(:)
            s_bfgs(:,:) = 0.0
            iter_bfgs = 1
            stepsize = (1.0/(1.0*min(iter_avg,avgsz)))*sum(avgstepsz(1:min(iter_avg,avgsz)))
          end if
        end if
      else
        upaccept = 0
        badaccept = 0
      end if
      if ((redux_count.eq.0).OR.(mini_mode.eq.2)) then
        stepsize = stepsize*1.618034
      else
        stepsize = (1.0/(1.0*min(iter_avg,avgsz)))*sum(avgstepsz(1:min(iter_avg,avgsz)))
      end if
      redux_count = 0
!
!   If peak energy is not lower, reset structure and try a smaller step
    else
!
!     reduce step size
      stepsize = stepsize*0.618034
      redux_count = redux_count + 1
      pmode = 1
!
!     explicitly undo last move (in internal space without overwriting ztorpr)
      call set_position(pmode,d_pos_vec,nrdof,mshfs)
      pmode = 0
!
!     give up if max reduction in stepsize is reached (usually indicates an uphill attempt for CG and BFGS)
      if (redux_count.gt.redux_max) then
        if ((minimode.eq.2).AND.(grms.gt.(1.0*grmsmin))) then
          wrncnt(1) = wrncnt(1) + 1
          if (wrncnt(1).eq.wrnlmt(1)) then
            write(ilog,*) 'WARNING: The CG method has exceed the permissible number of stepsize reductions.&
 & Resetting direction to steepest-descent and continuing minimization since GRMS is far away from requested value.'
            write(ilog,*) 'This was warning #',wrncnt(1),' of this type not all of which may be displayed.'
            wrnlmt(1) = wrnlmt(1)*10
          end if
          cur_gradient_vec(:) = peak_ahead_vec(:)
          prev_gradient_vec(:) = cur_gradient_vec(:)
          prev_normcg(:) = 0.0
          if (iter_avg.gt.0) then
            stepsize = (1.0/(1.0*min(iter_avg,avgsz)))*sum(avgstepsz(1:min(iter_avg,avgsz)))
          else
            stepsize = 1.0
          end if
          redux_count = 0
        else if ((minimode.eq.3).AND.(grms.gt.(1.0*grmsmin))) then
          wrncnt(1) = wrncnt(1) + 1
          if (wrncnt(1).eq.wrnlmt(1)) then
            write(ilog,*) 'WARNING: The BFGS method has exceed the permissible number of stepsize reductions.&
 & This means either that a local minimum was found or that the current direction was erroneous. Resetting&
 & Hessian and continuing minimization since GRMS is far away from requested value.'
            write(ilog,*) 'This was warning #',wrncnt(1),' of this type not all of which may be displayed.'
            wrnlmt(1) = wrnlmt(1)*10
          end if
          y_bfgs(:,:) = 0.0
          cur_gradient_vec(:) = peak_ahead_vec(:)
          y_bfgs(:,1) = cur_gradient_vec(:)
          s_bfgs(:,:) = 0.0
          iter_bfgs = 1
          if (iter_avg.gt.0) then
            stepsize = (1.0/(1.0*min(iter_avg,avgsz)))*sum(avgstepsz(1:min(iter_avg,avgsz)))
          else
            stepsize = 1.0
          end if
        else
          postanalyze = .true.
          esave = eprevious
          write(ilog,*) 'Termination: Maximum stepsize reductions reached during minimization. This indicates a&
 &n unreasonably small GRMS request, a nearly singular function in the energy landscape, alternative termination&
 & of the BFGS method (GRMS near target) or possibly a bug.'
          exit
        end if
      end if
!
    end if
!
!  
!   Make minimizing move, mode based on minimode:
!
!   1) steepest decent
    if (minimode.eq.1) then
      d_pos_vec(:) = stepsize*istepsize*cur_gradient_vec(:)
!
!    2) conj. grad. (Polak and Ribiere form)
    else if(minimode.eq.2) then
      if(iter.eq.1) then
        d_pos_vec(:) = stepsize*istepsize*cur_gradient_vec(:)
      else
        scalefactor = DOT_PRODUCT(cur_gradient_vec(:),cur_gradient_vec(:)) /&
    &    DOT_PRODUCT(prev_gradient_vec(:),prev_gradient_vec(:))
        normcg(:) = cur_gradient_vec(:) + scalefactor*prev_normcg(:)
        d_pos_vec(:) = istepsize*stepsize*normcg(:)
      end if       
!
!   3) quasi-newton L-BGFS method (Nocedal method)
    else if(minimode.eq.3) then
      if ((iter_bfgs.eq.1).OR.(iter.eq.1)) then 
        d_pos_vec(:) = stepsize*istepsize*cur_gradient_vec(:)  
        curve = 0.0
      else
      !
      ! Two Loop Recursion (Nocedal)
        if(iter_bfgs.le.(mini_mem+1)) then 
          incr = 0
          bound = iter_bfgs
        else
          incr = iter_bfgs - mini_mem
          bound = mini_mem
        end if
        q_bfgs(:) = cur_gradient_vec(:)
!
        j1 = MOD(iter_bfgs,mini_mem)
        if (j.eq.0) then
          j1 = mini_mem
        end if
        j2 = j1 - 1
        if (j2.eq.0) j2 = mini_mem
!
        do i=bound-1,1,-1
!
          j1 = MOD(i + incr,mini_mem)+1
          j2 = j1-1
          if (j1.eq.1) then
            j2 = mini_mem
          end if
!
          sigma = 0.0
          do t=1,nrdof
            sigma = sigma + (y_bfgs(t,j1)-y_bfgs(t,j2))*s_bfgs(t,j1)
          end do
          a_bfgs(i) = 0.0
          do t=1,nrdof
            a_bfgs(i) = a_bfgs(i) + (1.0/sigma)*s_bfgs(t,j1)*q_bfgs(t)
          end do
          q_bfgs(:) = q_bfgs(:) - a_bfgs(i)*(y_bfgs(:,j1)-y_bfgs(:,j2))
        end do
!
        hnot_bfgs = 0.0
        do t=1,nrdof
          hnot_bfgs = hnot_bfgs + (y_bfgs(t,j1)-y_bfgs(t,j2))*(y_bfgs(t,j1)-y_bfgs(t,j2))
        end do
        r_bfgs(:) = q_bfgs(:)*sigma/hnot_bfgs
!
!       second loop
        do i=1,bound-1
!
          j1 = MOD(i + incr,mini_mem)+1
          j2 = j1-1
          if (j1.eq.1) then
            j2 = mini_mem
          end if
!
          sigma = 0.0
          do t=1,nrdof
            sigma = sigma + (y_bfgs(t,j1)-y_bfgs(t,j2))*(s_bfgs(t,j1))
          end do
!          write(*,*) 'S2',sigma
          b_bfgs = 0.0
          do t=1,nrdof
            b_bfgs = b_bfgs + (1.0/sigma)*(y_bfgs(t,j1)-y_bfgs(t,j2))*r_bfgs(t)
          end do
          r_bfgs(:) = r_bfgs(:) + (a_bfgs(i)-b_bfgs)*(s_bfgs(:,j1))       
        end do
!
!     Normalize to the unit direction
      sigma = -stepsize
      d_pos_vec(:) = sigma*r_bfgs(:)/sqrt(DOT_PRODUCT(r_bfgs(:),r_bfgs(:)))
      curve = DOT_PRODUCT(r_bfgs(:),cur_gradient_vec(:))
      end if
!
    end if
!
!    BFGS: update stored position values
    if (minimode.eq.3) then
      j = MOD(iter_bfgs,mini_mem) + 1
      s_bfgs(:,j) = d_pos_vec(:)
    end if
!
!   set peak position, backing up current position
    call set_position(pmode,d_pos_vec,nrdof,mshfs)
!
  end do ! iterations of MASTER LOOP
!
! analysis if converged
  if (postanalyze.EQV..true.) then
    call CPU_time(t1)
    call mcstat(iter,ndump)
    call CPU_time(t2)
    time_analysis = time_analysis + t2 - t1
  end if
!    
! clean up
  if(minimode.eq.3) then
    deallocate(s_bfgs)
    deallocate(y_bfgs)
    deallocate(a_bfgs)
    deallocate(q_bfgs)
    deallocate(r_bfgs)
  end if
  deallocate(cur_gradient_vec)
  deallocate(prev_gradient_vec)
  deallocate(d_pos_vec)
  deallocate(prev_normcg)
  deallocate(peak_gradient_vec)
  deallocate(mshfs)
  deallocate(normcg)
  deallocate(avgstepsz)
!
! reached max iterations without convergence
  if (evals.gt.maxiter) then
    esave = eprevious
    do imol=1,nmol
      call getref_formol(imol)
    end do
    write(ilog,*) 'Termination: Max number of energy evaluations reached'
  end if
!
! Write summary
!      
  write(ilog,*) 'Minimization Summary:'
  write(ilog,*) 'Final Energy: ', esave
  write(ilog,*) 'Total iterations:   ', min(iter,maxiter)
  write(ilog,*) 'Energy evaluations: ', min(evals,maxiter)
!  
end
!
!----------------------------------------------------------------------
!
! The two following subroutines are helpers for the minimization functions
!
!----------------------------------------------------------------------
!
! populate the gradient vector (notice that the order of degrees of freedom in
! vec_tmp in internal coordinate space is nontrivial (follows recursive structure) 
!
subroutine get_gradient(vec_tmp,grms_tmp,size_vec)
!   
  use mini
  use molecule
  use forces
  use system
  use iounit
  use atoms
  use zmatrix
!
  implicit none
!
  integer size_vec,apos,imol,i,j,ttc
  RTYPE vec_tmp(size_vec),grms_tmp
!
  apos = 0
!
  if (fycxyz.eq.1) then
!
    do imol=1,nmol
      if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
         write(ilog,*) 'Fatal. Two-atom molecules are not yet support&
  &ed in minimization runs.'
         call fexit()
      end if
      do ttc=1,3
        if (dc_di(imol)%frz(ttc).EQV..true.) cycle
        apos = apos + 1
        vec_tmp(apos) = dc_di(imol)%f(ttc)*mini_xyzstep
      end do
!
      if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
        do i=-2,dc_di(imol)%maxntor
          if (i.le.0) then
            ttc = i+6 ! RB rotation
          else
            ttc =dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
          end if
          if (ttc.le.0) cycle
          if (dc_di(imol)%frz(ttc).EQV..true.) cycle
          apos = apos + 1
          if (i.le.0) then
            vec_tmp(apos) = dc_di(imol)%f(ttc)*mini_rotstep
          else
            vec_tmp(apos) = dc_di(imol)%f(ttc)*mini_intstep
          end if 
        end do
      end if
    end do
!
  else if (fycxyz.eq.2) then
    do i=1,n
      do j=1,3
        if (cart_frz(i,j).EQV..true.) cycle
        apos = apos + 1
        vec_tmp(apos) = cart_f(i,j)*mini_xyzstep
      end do
    end do
!
  end if
!
! bug check
  if (apos.ne.size_vec) then
    write(ilog,*) 'Fatal. Number of degrees of freedom does not ma&
 &tch up in get_gradient(...). This is most certainly a bug.'
    call fexit()
  end if
!
! calculate GRMS for convergence
  grms_tmp = sqrt(DOT_PRODUCT(vec_tmp,vec_tmp)/(1.0*size_vec))
!
end subroutine get_gradient
!
!----------------------------------------------------------------------
!
! propagate coordinates based on increments in d_pos_vec
!
subroutine set_position(mode,d_pos_vec,size_vec,mshfs)
!
  use iounit
  use molecule
  use forces
  use system
  use zmatrix
  use mini
  use cutoffs
  use movesets
  use math
  use atoms
!
  implicit none
! 
  integer mode,imol,i,ttc,size_vec,apos,aone,rs
  RTYPE d_pos_vec(size_vec),mshfs(3,nmol)
!
  aone = 1
!
  if (mode.eq.1) then
    d_pos_vec(:) = -1.0*d_pos_vec(:) ! for undoing a move
  else if (mode.eq.0) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Called set_position(...) with unsupported mode (related to minimize(...)). &
 &This is most certainly a bug.'
    call fexit()
  end if
!
  apos = 0
!
  if (fycxyz.eq.1) then
!
    do imol=1,nmol
      do ttc=1,3
        if (dc_di(imol)%frz(ttc).EQV..true.) then
          cur_trans(ttc) = 0.0
          cycle
        end if
        apos = apos + 1
        cur_trans(ttc) = d_pos_vec(apos)
      end do
!
      if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
        do i=-2,dc_di(imol)%maxntor
          if (i.le.0) then
            ttc = i+6 ! RB rotation
          else
            ttc =dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
          end if
          if (ttc.le.0) cycle
          if (dc_di(imol)%frz(ttc).EQV..true.) then
            dc_di(imol)%incr(ttc) = 0.0
          else
            apos = apos + 1
            dc_di(imol)%incr(ttc) = d_pos_vec(apos)
          end if
        end do
!       correct RB rotation increment
        cur_rot(1:3) = dc_di(imol)%incr(4:6)/RADIAN
!       note that this will edit the Z-matrix
        if (dc_di(imol)%maxntor.gt.0) call tmdmove(imol,mode)
      end if
!
      if (mode.eq.0) then
        call makeref_formol(imol)
!        call makeref_polym(imol)
        call IMD_prealign(imol,skip_frz)
        call makexyz_formol(imol)
!
!       finally: increment coordinates
!              note rotation is done first as it relies on comm
        if ((atmol(imol,2)-atmol(imol,1)).gt.0) then
          call rotxyzm(imol,aone)
        end if
        call transxyzm(imol,aone,mshfs(:,imol))
!       due to torsional moves we need to recompute rigid-body coordinates
        call update_rigidm(imol)
!       update grid association if necessary
        if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
          do rs=rsmol(imol,1),rsmol(imol,2)
            call updateresgp(rs)
          end do
        end if
      else
        call getref_formol(imol)
        call update_rigidm(imol)
      end if
!
    end do
!
  else if (fycxyz.eq.2) then
!
    if (mode.eq.0) then
      do imol=1,nmol
        call makeref_formol(imol)
      end do
      do i=1,n
        do ttc=1,3
          if (cart_frz(i,ttc).EQV..true.) cycle
          apos = apos + 1
          if (ttc.eq.1) x(i) = x(i) + d_pos_vec(apos)
          if (ttc.eq.2) y(i) = y(i) + d_pos_vec(apos)
          if (ttc.eq.3) z(i) = z(i) + d_pos_vec(apos)
        end do
      end do
!
!     loop over all molecules
      do imol=1,nmol 
        call update_rigidm(imol)
        call update_comv(imol)
        call genzmat(imol)
!       and shift molecule into central cell if necessary
        call update_image(imol)
!       update grid association if necessary
        if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
          do rs=rsmol(imol,1),rsmol(imol,2)
            call updateresgp(rs)
          end do
        end if
      end do
!     update pointer arrays such that analysis routines can work properly 
      call zmatfyc2()
!
    else if (mode.eq.1) then
      apos = size_vec
      do imol=1,nmol
        call getref_formol(imol)
        call update_rigidm(imol)
        call update_comv(imol)
        call genzmat(imol)
!       and shift molecule into central cell if necessary
        call update_image(imol)
!       update grid association if necessary
        if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
          do rs=rsmol(imol,1),rsmol(imol,2)
            call updateresgp(rs)
          end do
        end if
      end do
!     update pointer arrays such that analysis routines can work properly 
      call zmatfyc2()
    end if
!
  end if
!
!  bug check
  if (apos.ne.size_vec) then
    write(ilog,*) 'Fatal. Number of degrees of freedom does not ma&
 &tch up in set_position(...). This is most certainly a bug.'
    call fexit()
  end if
!    
end subroutine set_position
!
!----------------------------------------------------------------------
!
! this routine takes advantage of the built-in dynamics routines while externally modifying 
! the target bath temperature 
!
subroutine stochastic_min(maxiter)
!  
  use iounit
  use molecule
  use forces
  use system
  use mcsums
  use mini
  use energies
  use cutoffs
  use atoms
  use shakeetal
!    
  implicit none
!    
  integer i,istep,maxiter,collectend,ndump,nstepi,nrdof,schkbuf,imol,j,dumm
  RTYPE tmax,targetrate,instrate,t1,t2,grms,buke,butau
  RTYPE, ALLOCATABLE:: gradvec(:)
  logical postanalyze,afalse,atrue
!
  dumm = 100
  atrue = .true.
  afalse = .false.
  grms = 0.0
!    
! take a few steps of steepest descent to relax major energy barriers
  nstepi = maxiter
  if (mini_sc_sdsteps.gt.0) then
    if (no_shake.EQV..false.) then
      write(ilog,*) 'Warning. Skipping requested steepest-descent steps due to presence of holonomic constraints for &
 &stochastic minimizer.'
    else
      j = 1
      call minimize(mini_sc_sdsteps,grms,mini_stepsize,j)
    end if
  end if
!
  nrdof = 0
! "vectorize" all degrees of freedom into a single vector with dimension "nrdof"
  if (fycxyz.eq.1) then
    do imol=1,nmol
      do j=1,size(dc_di(imol)%frz)
        if (dc_di(imol)%frz(j).EQV..false.) then
          nrdof = nrdof + 1
        end if
      end do
    end do 
  else if (fycxyz.eq.2) then
! SHAKE-type constraints and drift removal must be handled separately (not mappable)
    do i=1,n
      do j=1,3
        if (cart_frz(i,j).EQV..true.) cycle
        nrdof = nrdof + 1
      end do
    end do
  end if
  if (nrdof.le.0) then
    write(ilog,*) 'Fatal. No remaining degrees of freedom to minimize. Relieve constraints.'
    call fexit()
  end if
  allocate(gradvec(nrdof))
!
  if (mini_sc_tbath.gt.kelvin) then
    write(ilog,*) 'Warning. For the stochastic minimizer, the supplied target temperature is larger than the &
 &initial one. This run will most likely terminate immediately after the number of steps defined by FMCSC_MINI_SC_HEAT.'
  end if
  buke = kelvin
  butau = tstat%params(1)
  kelvin = max(kelvin,mini_sc_tbath) ! target temperature (T of bath)
  collectend=floor(nstepi*mini_sc_heat)
  tmax=0.0
  ndump = 0
!    
  istep = 1
  if (fycxyz.eq.2) then
    call cart_mdmove(istep,dumm,atrue)
  else if (fycxyz.eq.1) then
    call int_mdmove(istep,dumm,atrue)
  else
    write(ilog,*) 'Fatal. Called stochastic_min(...) with unsupported system representation. This is an omission bug.'
    call fexit()
  end if
  do istep=2,nstepi
    schkbuf = nsancheck 
    nsancheck = nstepi+1 ! suppress ensemble output from MD
    if (fycxyz.eq.2) then
      call cart_mdmove(istep,dumm,afalse)
    else if (fycxyz.eq.1) then
      call int_mdmove(istep,dumm,afalse)
    end if
    nsancheck = schkbuf
    call get_gradient(gradvec,grms,nrdof)
    if(mod(istep,nsancheck).eq.0) then
 12  format('Step: ',i8,' Current Energy: ',g14.6,' Current Temperature: ',g14.6, ' GRMS: ',g14.6)
      write(ilog,12) istep,ens%insU,ens%insT,grms
    end if
!
    if(istep.lt.collectend) then
      if(tmax.lt.ens%insT) then  !T is initialized to FMCSC_TEMP, but it will rise as
        tmax = ens%insT          !PE is converted to KE.
      end if
    end if 
!      
! Monitor gradients here and terminate when criteria reached.
    if (istep.gt.collectend) then
      if (grms.lt.mini_econv) then
 14 format('Normalized RMS gradient is ',g10.4,' and has dropped below convergence criterion&
 & of ',g10.4,'.')
        write(ilog,14) grms,mini_econv
        postanalyze = .true.
        exit
      else if(ens%insT.lt.mini_sc_tbath) then
 13 format('System temperature is ',g10.4,' and has dropped below ',g10.4,'. Terminating.')
        write(ilog,13) ens%insT,mini_sc_tbath
        postanalyze = .true.
        exit
      else if(istep .eq. nstepi) then
        write(ilog,*) 'Maximum number of stochastic minimization steps exceeded. Forced termination.'
        postanalyze = .true.
        exit
      end if
    end if
!      
! Set target rate based on number of steps given
    if(istep.eq.collectend) then
      targetrate=1.25*(tmax-mini_sc_tbath)/(nstepi-collectend)
    end if
!      
! Dynamically update target temperature and rate for faster end convergence (should do after heating MD steps)
    if (istep.gt.collectend) then
      kelvin = mini_sc_tbath + buke*(exp(-((istep-collectend)/(0.2*(nstepi-collectend)))**2)) - buke*1.389e-11
      instrate=(tmax - ens%insT)/(istep-collectend)
      if ((instrate.gt.0.0).AND.(targetrate.gt.0.0)) then
        tstat%params(1) = tstat%params(1)*(instrate/targetrate)
        if (tstat%params(1).lt.(10.*dyn_dt)) then
          tstat%params(1) = 10.0*dyn_dt
        end if
        if (tstat%params(1).gt.butau) then
          tstat%params(1) = butau
        end if
      end if
    end if
!
    call CPU_time(t1)
    call mcstat(istep,ndump)
    call CPU_time(t2)
    time_analysis = time_analysis + t2 - t1
!     
  end do
!
  deallocate(gradvec)
!    
! last-step analysis if converged
  if (postanalyze.EQV..true.) then
    call CPU_time(t1)
    call mcstat(istep,ndump)
    call CPU_time(t2)
    time_analysis = time_analysis + t2 - t1
  end if
!   
! Write summary
!      
  write(ilog,*) 'Minimization Summary:'
  write(ilog,*) 'Final Energy: ', esave
  write(ilog,*) 'Total iterations: ', min(istep,nstepi)
!    
end
!
!---------------------------------------------------------------------------------------------------------
!
