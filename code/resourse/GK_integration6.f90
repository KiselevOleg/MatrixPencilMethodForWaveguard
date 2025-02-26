    module GK_integration6
        
        integer,parameter:: N_Int_Nodes=15
        integer Inf_Key/0/
        real*8 RV_eps_step_increment
        real*8,private:: GK_nodes(N_Int_Nodes),K_weights(N_Int_Nodes),G_weights(N_Int_Nodes)
        
    contains
        complex(8) function GK_integral_ab(functionName,accurate,a,b,upperPolesValue,depthOfAvoidingPoles) result(f)
        implicit none
            complex(8),external::functionName
            real(8),intent(in)::accurate
            real(8),intent(in)::a
            real(8),intent(in)::b
            real(8),intent(in)::upperPolesValue
            real(8),intent(in)::depthOfAvoidingPoles
            
            f=GK_integral(fun,accurate,b-a,upperPolesValue,depthOfAvoidingPoles)
        contains
            complex(8) function fun(alpha)
            implicit none
                complex(8),intent(in)::alpha
                
                fun=functionName(alpha+a)
            endfunction fun
        endfunction GK_integral_ab
        complex(8) function GK_integral(functionName,accurate,lengthOfIntegration,upperPolesValue,depthOfAvoidingPoles) result(f)
        implicit none
            complex(8),external::functionName
            real(8),intent(in)::accurate
            real(8),intent(in)::lengthOfIntegration
            real(8),intent(in)::upperPolesValue
            real(8),intent(in)::depthOfAvoidingPoles
        
            complex(8) Rd(1)
            
            call DINN5_GK(KQ_int,0d0,0d0,0d0,upperPolesValue,depthOfAvoidingPoles,0d0,accurate,1d-1,lengthOfIntegration,1,Rd)
            f=Rd(1)
        contains
            subroutine KQ_int(alpha1,mass,n_)
            implicit none
                complex(8),intent(in)::alpha1
                complex(8),intent(out)::mass(n_)
                integer,intent(in)::n_
                
                mass(1)=functionName(alpha1)
            endsubroutine KQ_int
        endfunction GK_integral
        subroutine GK_7_15_init
            GK_nodes(1)=-0.991455371120813d0; GK_nodes(2)=-0.949107912342759d0
            GK_nodes(3)=-0.864864423359769d0; GK_nodes(4)=-0.741531185599394d0;
            GK_nodes(5)=-0.586087235467691d0; GK_nodes(6)=-0.405845151377397d0;
            GK_nodes(7)=-0.207784955007898d0; 
            GK_nodes(8)= 0d0;
            GK_nodes(9)= 0.207784955007898d0; GK_nodes(10)=0.405845151377397d0;
            GK_nodes(11)=0.586087235467691d0; GK_nodes(12)=0.741531185599394d0;
            GK_nodes(13)=0.864864423359769d0; GK_nodes(14)=0.949107912342759d0;
            GK_nodes(15)=0.991455371120813d0;
            
            K_weights(1)=0.022935322010529d0; K_weights(2)=0.063092092629979d0
            K_weights(3)=0.104790010322250d0; K_weights(4)=0.140653259715525d0
            K_weights(5)=0.169004726639267d0; K_weights(6)=0.190350578064785d0
            K_weights(7)=0.204432940075298d0; 
            K_weights(8)=0.209482141084728d0
            K_weights(9)=0.204432940075298d0; K_weights(10)=0.190350578064785d0
            K_weights(11)=0.169004726639267d0; K_weights(12)=0.140653259715525d0
            K_weights(13)=0.104790010322250d0; K_weights(14)=0.063092092629979d0
            K_weights(15)=0.022935322010529d0;

            
            G_weights(1)=0d0; G_weights(2)=0.129484966168870d0
            G_weights(3)=0d0; G_weights(4)=0.279705391489277d0
            G_weights(5)=0d0; G_weights(6)=0.381830050505119d0
            G_weights(7)=0d0; 
            G_weights(8)=0.417959183673469d0
            G_weights(9)=0d0; G_weights(10)=0.381830050505119d0
            G_weights(11)=0d0; G_weights(12)=0.279705391489277d0
            G_weights(13)=0d0; G_weights(14)=0.129484966168870d0
            G_weights(15)=0d0;           
        end subroutine GK_7_15_init        
        !==========================!
        !==========================!
        subroutine GK_adaptive_int(int_func,a,b,int_h,eps,ret_arr,N)
        implicit none;
            integer N
            real*8 eps,eps_out,int_h
            complex*16 a,b,t_i_h,t_x
            complex*16 ret_arr(N),ret_arr_0(N)
            external:: int_func
    
                ret_arr=0d0
                
                t_x = a;
                t_i_h = dcmplx(int_h,0d0)
                if (abs(imag(b)-imag(a)).gt.eps) then 
                    if (imag(b).lt.imag(a)) then
                        t_i_h = dcmplx(0d0,-int_h)
                    else
                        t_i_h = dcmplx(0d0,int_h)
                    endif
                endif
                do while (abs(b-t_x).gt.eps)
                    if (abs(t_x + t_i_h).gt.abs(b)) t_i_h = b - t_x
                    call GK_int(int_func,t_x,t_x+t_i_h,ret_arr_0,eps_out,N)
                    if (eps_out.gt.eps) then
                        t_i_h=t_i_h*0.5d0
                    else
                        ret_arr=ret_arr+ret_arr_0
                        t_x=t_x+t_i_h
!                        print*,t_x
                        if (eps_out.lt.1d-3*eps) t_i_h=t_i_h*2d0
                    endif
                enddo
                int_h=abs(t_i_h)
        end subroutine GK_adaptive_int
        !==========================!
        !==========================!
        subroutine GK_adaptive_int_inf(int_func,a,b,int_h,eps,ret_arr,N)
        implicit none;
            integer N,i,ipr,it
            real*8 eps,int_h,temp_arr(N),t1,eps10,pm,pt,t,int_h_1,t_h
            complex*16 a,b,t_x_a,t_x_b
            complex*16 ret_arr(N),ret_arr_0(N)
            external:: int_func
    
                ret_arr=0d0
                
                t_x_b = a;  t=abs(b-a)
                eps10=eps*10d0; ipr=0; it=1

                ret_arr = 0d0
                
             	do while(it)
!                    print*,t_x_b
	                t_x_a=t_x_b;  t_x_b=t_x_b+int_h;  t1=abs(t_x_b-a)
	                if(t1 > t)then
	                  t_x_b=b; it=-1
                    end if
                    
                    int_h_1=int_h
                    call GK_adaptive_int(int_func,t_x_a,t_x_b,int_h_1,eps,ret_arr_0,N)
                    if (minval(abs(ret_arr_0)).lt.1d-10.and.abs(t_x_a).gt.50.0d0.and.maxval(abs(ret_arr_0)).lt.1d-7) return
            
            ! for cheking-up the convergence at infinity
            
                    t_h=abs(t_x_b-t_x_a)
            		temp_arr=abs(ret_arr_0)/t_h
            
                    ret_arr=ret_arr+ret_arr_0
                    
                    if(abs(int_h) < 10d0*abs(int_h_1)) then
	                    int_h=4*int_h;  else;   int_h=4*int_h_1
	                end if
            
            ! at infinity
            
                    pm=0
                    do i=1,N
                        if(abs(ret_arr(i)) > 1d-15)then
            		        pt=abs(temp_arr(i)*t1/ret_arr(i))
            		        if(pt > pm) pm=pt
            		    end if
                    end do
            
                    if(pm < eps10)then
                        ipr=ipr+1
            		    if(ipr > 4)return
            		else
            		    ipr=0
            	    end if ! pm < eps10
                enddo
        end subroutine GK_adaptive_int_inf
        !==========================!
        !==========================!
        subroutine GK_int(int_func,a,b,ret_arr,eps_out,N)
        implicit none;
            integer N,i
            real*8 eps_out
            complex*16 a,b,GK_nodes_arb(N_Int_Nodes),K_weights_arb(N_Int_Nodes),G_weights_arb(N_Int_Nodes)
            complex*16 ret_arr(N),temp_arr(N,2)
            external:: int_func
    
                GK_nodes_arb=0.5d0*(b-a)*GK_nodes+0.5d0*(b+a)
                K_weights_arb=K_weights*0.5d0*(b-a); G_weights_arb=G_weights*0.5d0*(b-a)
                
                temp_arr=0d0
                do i=1,N_Int_Nodes
                    call int_func(GK_nodes_arb(i),ret_arr,N)
                    temp_arr(:,1)=temp_arr(:,1)+K_weights_arb(i)*ret_arr
                    temp_arr(:,2)=temp_arr(:,2)+G_weights_arb(i)*ret_arr
                enddo
                eps_out=maxval(abs(temp_arr(:,1)-temp_arr(:,2)))
                ret_arr=temp_arr(:,1)
        end subroutine GK_int
        !==========================!
        !==========================!
        subroutine GK_int_old(int_func,a,b,ret_arr,eps_out,N)
        implicit none;
            integer N
            real*8 eps_out
            complex*16 a,b,GK_nodes_arb(N_Int_Nodes),K_weights_arb(N_Int_Nodes),G_weights_arb(N_Int_Nodes)
            complex*16 ret_arr(N),int_out(N,2)
            external:: int_func
    
                GK_nodes_arb=0.5d0*(b-a)*GK_nodes+0.5d0*(b+a)
                K_weights_arb=K_weights*0.5d0*(b-a); G_weights_arb=G_weights*0.5d0*(b-a)
                
                call int_func(GK_nodes_arb,K_weights_arb,G_weights_arb,N_Int_Nodes,int_out,N)
                eps_out=maxval(abs(int_out(:,1)-int_out(:,2)))
                ret_arr=int_out(:,1)
        end subroutine GK_int_old
        !==========================!
        !==========================!
        subroutine DINN5_GK(CF,t1,t2,t3,t4,tm,tp,eps,pr,gr,N,Rd)
        implicit none
        	integer N,ib,inf,i
        	real*8 t1,t2,t3,t4,tm,tp,eps,pr,gr, int_h
    	    complex*16 Rd(N),a,b,sb(N),h
            external CF
            common/simp5/h,ib,inf
            
            call GK_7_15_init
    	    Rd=0d0
            ! [0, t1]
            a=0; b=t1; int_h=0.25d0*abs(b-a)
            call GK_adaptive_int(CF,a,b,int_h,eps,sb,N)
            Rd=Rd+sb
        
    	    if(t3-t2 < eps)then   ! no inverse poles case
                ! [t1, t1-i*tm]    
                a=b; b=cmplx(t1,-tm); int_h=0.25d0*abs(b-a)
                call GK_adaptive_int(CF,a,b,int_h,eps,sb,N)
                Rd=Rd+sb
    
                !  [t1-i*tm, t4-i*tm]    
                a=b; b=cmplx(t4,-tm); int_h=0.1d0*abs(b-a)
                call GK_adaptive_int(CF,a,b,int_h,eps,sb,N)
                Rd=Rd+sb
    
                !  [t4-i*tm, t4]
                a=b; b=t4; int_h=0.25d0*abs(b-a)
                call GK_adaptive_int(CF,a,b,int_h,eps,sb,N)
                Rd=Rd+sb
            else ! t2 < t3 - an inverse pole case
                if (t2-t1 > eps) then
                    ! first deviation from below
                    !   [t1, t1-i*tm]    
                    a=b; b=cmplx(t1,-tm); int_h=0.5d0*abs(b-a)
                    call GK_adaptive_int(CF,a,b,int_h,eps,sb,N)
                    Rd=Rd+sb
    
                    !   [t1-i*tm,t2-i*tm]    
                    a=b; b=cmplx(t2,-tm); int_h=0.5d0*abs(b-a)
                    call GK_adaptive_int(CF,a,b,int_h,eps,sb,N)
                    Rd=Rd+sb
                end if ! t2 > t1
                
                ! diviation from above
                !   [b,t2+i*tp]    
                a=b; b=cmplx(t2,tp); int_h=0.5d0*abs(b-a)
                call GK_adaptive_int(CF,a,b,int_h,eps,sb,N)
                Rd=Rd+sb
    	   
                !  [t2+i*tp, t3+i*tp]
                a=b; b=cmplx(t3,tp); int_h=0.5d0*abs(b-a)
                call GK_adaptive_int(CF,a,b,int_h,eps,sb,N)
                Rd=Rd+sb
    
                !  [t3+i*tp, t3-i*tm]
    	        a=b; b=cmplx(t3,-tm); int_h=0.5d0*abs(b-a)
                call GK_adaptive_int(CF,a,b,int_h,eps,sb,N)
                Rd=Rd+sb
    
                ! second diviation from below
                !  [t3-i*tm, t4-i*tm] 
                a=b; b=cmplx(t4,-tm); int_h=0.25d0*abs(b-a)
                call GK_adaptive_int(CF,a,b,int_h,eps,sb,N)
                Rd=Rd+sb
    
                !  [t4-i*tm, t4]
                a=b; b=t4; int_h=0.5d0*abs(b-a)
                call GK_adaptive_int(CF,a,b,int_h,eps,sb,N)
                Rd=Rd+sb
            end if ! t3 > t2
    
            ! [t4, inf.] 
        	a=b; b=gr; inf=1; int_h=0.33d0 !0.25d0*abs(b-a)
            if (Inf_Key) then
                call GK_adaptive_int_inf(CF,a,b,int_h,eps*1d1,sb,N)
            else
                call GK_adaptive_int(CF,a,b,int_h,eps*1d1,sb,N)
            endif

            Rd=Rd+sb
        
!            h=0.5d0; ib=1;
!            a=b; b=gr; inf=1
!            call CDINN5(CF,a,b,eps*1d3,pr,N,sb,Rd)
            
        end subroutine DINN5_GK 
        !==========================!
        !==========================!        
        subroutine GK_adaptive_int_RealV(int_func,a,b,int_h,eps,ret_arr,N)
        implicit none;
            integer N
            real*8 eps,eps_out,int_h,a,b,t_i_h,t_x
            complex*16 ret_arr(N),ret_arr_0(N)
            external:: int_func
    
!                ret_arr=0d0
                
                t_x = a;
                t_i_h = int_h
                do while (b-t_x.gt.eps)
                    if (t_x + t_i_h.gt.b) t_i_h = b - t_x
                    call GK_int_RealV(int_func,t_x,t_x+t_i_h,ret_arr_0,eps_out,N)
                    if (eps_out.gt.eps) then
                        t_i_h=t_i_h*0.5d0
                    else
                        ret_arr=ret_arr+ret_arr_0
                        t_x=t_x+t_i_h
!                        print*,t_x
                        if (eps_out.lt.RV_eps_step_increment) t_i_h=t_i_h*1.5d0
                    endif
                enddo
!                int_h=t_i_h
        end subroutine GK_adaptive_int_RealV
        !==========================!
        !==========================!
        subroutine GK_int_RealV(int_func,a,b,ret_arr,eps_out,N)
        implicit none;
            integer N,i
            real*8 eps_out,a,b,GK_nodes_arb(N_Int_Nodes),K_weights_arb(N_Int_Nodes),G_weights_arb(N_Int_Nodes)
            complex*16 ret_arr(N),temp_arr(N,2)
            external:: int_func
    
                GK_nodes_arb=0.5d0*(b-a)*GK_nodes+0.5d0*(b+a)
                K_weights_arb=K_weights*0.5d0*(b-a); G_weights_arb=G_weights*0.5d0*(b-a)
                
                temp_arr=0d0
                do i=1,N_Int_Nodes
                    call int_func(GK_nodes_arb(i),ret_arr,N)
                    temp_arr(:,1)=temp_arr(:,1)+K_weights_arb(i)*ret_arr
                    temp_arr(:,2)=temp_arr(:,2)+G_weights_arb(i)*ret_arr
                enddo
                ret_arr=temp_arr(:,1)-temp_arr(:,2);
                eps_out=(200d0*maxval(abs(ret_arr)))**1.5d0
                ret_arr=temp_arr(:,1)
        end subroutine GK_int_RealV        
    end module  GK_integration6
    

!========================================================================!
!========================================================================!
