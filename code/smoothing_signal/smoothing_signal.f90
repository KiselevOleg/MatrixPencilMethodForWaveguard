module smoothing_signal
implicit none
    public::arithmetic_mean_smoothing,wavelet_transform_smoothing
    public::wavelet_transform,inverse_wavelet_transform
    public::mexican_hat,morlet_wavelet,meyer_wavelet
    
    interface wavelet_transform
        procedure wavelet_transform_delta
        procedure wavelet_transform_mass
    endinterface wavelet_transform
    interface inverse_wavelet_transform
        procedure inverse_wavelet_transform_delta
        procedure inverse_wavelet_transform_mass
    endinterface inverse_wavelet_transform
    
    private::wavelet_transform_delta,wavelet_transform_mass
    contains
    
    subroutine wavelet_transform_smoothing(signal_size,x,signal,a_min,da,a_max,b_min,db,b_max)
    implicit none
        integer(4),intent(in)::signal_size
        real(8),intent(in)::x(signal_size)
        real(8),intent(inout)::signal(signal_size)
        real(8),intent(in)::a_min
        real(8),intent(in)::da
        real(8),intent(in)::a_max
        real(8),intent(in)::b_min
        real(8),intent(in)::db
        real(8),intent(in)::b_max
        
        real(8),allocatable::wavelet(:,:)
        
        integer(4) a_size,b_size
        
        real(8) a,b
        integer(4) i,j
        
        a_size=(a_max-a_min)/da+1
        b_size=(b_max-b_min)/db+1
        allocate(wavelet(a_size,b_size))
        
        i=1
        do a=a_min,a_max,da
            j=1
            do b=b_min,b_max,db
                wavelet(i,j)=wavelet_transform(morlet_wavelet,a,b,signal_size,x,signal)
                
                j=j+1
            enddo
            
            i=i+1
            !print*,i,a_size
        enddo
        
        print*,"inverse_wavelet_transform"
        do i=1,signal_size
            signal(i)=inverse_wavelet_transform(morlet_wavelet,a_size,a_min,da,b_size,b_min,db,wavelet,x(i))
            !print*,x(i),i,signal_size
        enddo
        
        deallocate(wavelet)
    endsubroutine wavelet_transform_smoothing
    
    real(8) function inverse_wavelet_transform_mass(wavelet_function,a_size,a,b_size,b,wavelet_result,t) result(f)
    implicit none
        real(8),external::wavelet_function!real(8) function wavelet_function(real(8) x)
        integer(4),intent(in)::a_size
        real(8),intent(in)::a(a_size)
        integer(4),intent(in)::b_size
        real(8),intent(in)::b(b_size)
        real(8),intent(in)::wavelet_result(a_size,b_size)
        real(8),intent(in)::t
        
        integer(4) i,j
        
        f=0d0
        do i=1,a_size-1
            do j=1,b_size-1
                f=f+wavelet_result(i,j)*phi_ab(a(i),b(j),t)*(a(i+1)-a(i))*(b(j+1)-b(j))
            enddo
        enddo
        
    contains
        real(8) function phi_ab(a,b,x) result(f)
        implicit none
            real(8),intent(in)::a
            real(8),intent(in)::b
            real(8),intent(in)::x
            
            f=1d0/sqrt(a)*wavelet_function((x-b)/a)
        endfunction phi_ab
    endfunction inverse_wavelet_transform_mass
    real(8) function inverse_wavelet_transform_delta(wavelet_function,a_size,a_0,da,b_size,b_0,db,wavelet_result,t) result(f)
    implicit none
        real(8),external::wavelet_function!real(8) function wavelet_function(real(8) x)
        integer(4),intent(in)::a_size
        real(8),intent(in)::a_0
        real(8),intent(in)::da
        integer(4),intent(in)::b_size
        real(8),intent(in)::b_0
        real(8),intent(in)::db
        real(8),intent(in)::wavelet_result(a_size,b_size)
        real(8),intent(in)::t
        
        real(8) a,b
        integer(4) i,j
        
        f=0d0
        a=a_0
        do i=1,a_size
            b=b_0
            do j=1,b_size
                f=f+wavelet_result(i,j)*phi_ab(a,b,t)*da*db
                b=b+db
            enddo
            a=a+da
        enddo
    contains
        real(8) function phi_ab(a,b,x) result(f)
        implicit none
            real(8),intent(in)::a
            real(8),intent(in)::b
            real(8),intent(in)::x
            
            f=1d0/sqrt(a)*wavelet_function((x-b)/a)
        endfunction phi_ab
    endfunction inverse_wavelet_transform_delta
    
    real(8) function wavelet_transform_delta(wavelet_function,a,b,signal_size,t_0,dt,signal) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        real(8),external::wavelet_function!real(8) function wavelet_function(real(8) x)
        real(8),intent(in)::a
        real(8),intent(in)::b
        integer(4),intent(in)::signal_size
        real(8),intent(in)::t_0
        real(8),intent(in)::dt
        real(8),intent(in)::signal(signal_size)
        
        real(8) t
        integer(4) i
        
        if(a<=epsilon) call print_error("smoothing_signal.wavelet_transform_mass","a<=epsilon")
        
        t=t_0
        
        f=0d0
        do i=1,signal_size-1
            f=f+wavelet_function((t-b)/a)*signal(i)*dt
            
            t=t+dt
        enddo
        f=f/sqrt(abs(a))
    endfunction wavelet_transform_delta
    real(8) function wavelet_transform_mass(wavelet_function,a,b,signal_size,time,signal) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        real(8),external::wavelet_function!real(8) function wavelet_function(real(8) x)
        real(8),intent(in)::a
        real(8),intent(in)::b
        integer(4),intent(in)::signal_size
        real(8),intent(in)::time(signal_size)
        real(8),intent(in)::signal(signal_size)
        
        integer(4) i
        
        if(a<=epsilon) call print_error("smoothing_signal.wavelet_transform_mass","a<=epsilon")
        
        f=0d0
        do i=1,signal_size-1
            f=f+wavelet_function((time(i)-b)/a)*signal(i)*(time(i+1)-time(i))
        enddo
        f=f/sqrt(abs(a))
    endfunction wavelet_transform_mass
    
    real(8) function mexican_hat(x) result(f)
    implicit none
        real(8),intent(in)::x
        
        f=(1d0-x*x)*exp(-x*x*0.5d0)
    endfunction mexican_hat
    real(8) function morlet_wavelet(x) result(f)
    use math,only:pi
    implicit none
        real(8),intent(in)::x
        !real(8),optional,intent(in)::sigma_
        
        real(8) sigma
        real(8) kappa_sigma,c_sigma
        
        !if(.not.present(sigma_)) then
        !    sigma=7d0
        !else
        !    sigma=sigma_
        !endif
        sigma=7d0
        
        kappa_sigma=exp(-0.5d0*sigma*sigma)
        c_sigma=1d0/sqrt(1d0+exp(-sigma*sigma)-2d0*exp(-0.75d0*sigma*sigma))
        
        f=c_sigma*pi**(-0.25d0)*exp(-0.5d0*x*x)*(cos(sigma*x)-kappa_sigma)
    endfunction morlet_wavelet
    real(8) function meyer_wavelet(x) result(f)
    use math,only:pi
    implicit none
        real(8),intent(in)::x
        
        if(abs(x)*3d0<pi+pi) then
            f=1/sqrt(pi+pi)
            return
        endif
        if(pi+pi<3d0*abs(x).and.3d0*abs(x)<pi+pi+pi+pi) then
            f=1/sqrt(pi+pi)*cos(pi*0.5d0*nu(1.5d0*abs(x)/pi-1d0))
            return
        endif
        f=0d0
    contains
        real(8) function nu(x) result(f)
        implicit none
            real(8),intent(in)::x
            
            if(x<0d0.or.x>1d0) then
                f=0d0
                return 
            endif
            
            f=x**4*(35d0-84d0*x+70d0*x**2-20d0*x**3)
        endfunction nu
    endfunction meyer_wavelet
    
    subroutine arithmetic_mean_smoothing(strong, signal_size, signal, do_zeros_borders)
    use system,only:print_error,print_warning
    implicit none
        integer(4),intent(in)::strong
        integer(4),intent(in)::signal_size
        real(8),intent(inout)::signal(signal_size)
        logical(1),intent(in)::do_zeros_borders
        
        real(8) tail_sum_signal(strong)
        integer(4) tail_sum_signal_start
        real(8) sum,sum_arithmetic_mean,sum_del
        
        integer(4) i
        
        if(strong<1.or.strong>50) call print_error("smoothing_signal.arithmetic_mean_smoothing","strong<1.or.>strong>50")
        if(signal_size<=strong+strong) call print_error("smoothing_signal.arithmetic_mean_smoothing","signal_size<s=trong+strong")
        if(signal_size<50) call print_warning("smoothing_signal.arithmetic_mean_smoothing","signal_size<50")
        
        tail_sum_signal_start=1
        sum_del=strong+1+strong
        
        sum=0d0
        do i=1,strong
            sum=sum+signal(i)
            tail_sum_signal(i)=signal(i)
        enddo
        sum=sum+signal(strong+1)
        do i=strong+2,strong+1+strong
            sum=sum+signal(i)
        enddo
        
        do i=strong+1,signal_size-strong-1
            sum_arithmetic_mean=sum/sum_del
            
            sum=sum-tail_sum_signal(tail_sum_signal_start)+signal(i+strong+1)
            tail_sum_signal(tail_sum_signal_start)=signal(i)
            tail_sum_signal_start=tail_sum_signal_start+1
            if(tail_sum_signal_start==strong+1) tail_sum_signal_start=1
            
            signal(i)=sum_arithmetic_mean
        enddo
        
        sum_arithmetic_mean=sum/sum_del
        signal(signal_size)=sum_arithmetic_mean
        
        if(do_zeros_borders) then
            do i=1,strong
                signal(i)=0d0
            enddo
            do i=signal_size-strong,signal_size
                signal(i)=0d0
            enddo
        endif
    endsubroutine arithmetic_mean_smoothing
endmodule smoothing_signal
