! real roots dzeta_k of a comlex-valued function F(t)
! (usually F(t)=Delta(al), t=real(al) - denominator of 
! Green's matrix K(al) for a multilayered half-space (MultiKD))
	
		
module Halfc_
implicit none
	
	contains
    subroutine Halfc(F,tmin,tmax,ht,eps,Nmax,dz,Nx)
    implicit none
    integer Nmax,Nx,it,ir,ii

    real*8 tmin,tmax,t1,t2,ht,eps,dz(Nmax),rf1,rf2,if1,if2, &
	  epsf,signr,signi,u1,u2,u,sgr,sgi,rfu1,rfu2,ifu1,ifu2, &
	  rfu,ifu

	complex*16 F,ft

! initial values	
    Nx=0; it=1; epsf=eps/1d7
	t1=tmin;  ft=F(dcmplx(t1))
    rf1=real(ft);	if1=imag(ft)    
!	open(unit=2,file='RD2.dat',status='replace')
! ht steps

 1  t2=t1+ht 
    if(t2 > tmax)then
      t2=tmax; it=-1
	end if
    ft=F(dcmplx(t2)); rf2=real(ft);	if2=imag(ft)

!    write(2,33)t2,ft
! 33 format('  t,ft=',f11.5,2d12.5)

    if(abs(rf2) < epsf)then
	  signr=-1; ir=-1
	else
	  signr=rf1*rf2; ir=1
	end if

    if(abs(if2) < epsf)then
	  signi=-1; ii=-1
	else
	  signi=if1*if2; ii=1
	end if

! simultaneous change of the signs on the interval [t1,t2]

    if((signr < 0).and.(signi < 0))then

! special case ReF=ImF=0
	  if((ir < 0).and.(ii < 0))then
	    Nx=Nx+1; dz(Nx)=t2
		go to 2
	  end if

	  u1=t1; rfu1=rf1; ifu1=if1
	  u2=t2; rfu2=rf2; ifu2=if2	  	  

! root's localisation
   
  3   u=(u1+u2)/2;  ft=F(dcmplx(u))
	  rfu=real(ft); ifu=imag(ft)

      sgr=rfu1*rfu; sgi=ifu1*ifu

!    write(2,34)u,ft
! 34 format('      u,fu=',f11.5,2d12.5)

!   non-synchronous signs change
	  if((ir > 0).and.(ii > 0).and.(sgr*sgi < 0))go to 2
		
	  if(ir > 0)then
		if(sgr > 0)then
	      u1=u; rfu1=rfu; ifu1=ifu
		else
	      u2=u; rfu2=rfu; ifu2=ifu
		end if
	  else
	    if(sgi > 0)then
	      u1=u; rfu1=rfu; ifu1=ifu
		else
	      u2=u; rfu2=rfu; ifu2=ifu
		end if
	  end if
	  
	  if(abs(u1-u2) > eps)go to 3

! picking up the root selected
      Nx=Nx+1; dz(Nx)=u

!  write(2,35)Nx,u
! 35 format('  Nx,u=',i5,d12.5)

 	end if

! looking for the next roots
 
 2	t1=t2; rf1=rf2; if1=if2
    if((Nx < Nmax).and.(it > 0))go to 1

    end                       
        
endmodule Halfc_