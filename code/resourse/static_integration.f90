module static_integration
implicit none
    public::riemann_sum,simpson_sum
    
    interface riemann_sum
        procedure riemann_sum_complex_mass_x_complex_f
        procedure riemann_sum_complex_delta_x_complex_f
        procedure riemann_sum_real_mass_x_complex_f
        procedure riemann_sum_real_delta_x_complex_f
        procedure riemann_sum_real_mass_x_real_f
        procedure riemann_sum_real_delta_x_real_f
    endinterface riemann_sum
    interface simpson_sum
        procedure simpson_sum_complex_mass_x_complex_f
        procedure simpson_sum_complex_delta_x_complex_f
        procedure simpson_sum_real_mass_x_complex_f
        procedure simpson_sum_real_delta_x_complex_f
        procedure simpson_sum_real_mass_x_real_f
        procedure simpson_sum_real_delta_x_real_f
    endinterface simpson_sum
    
    private riemann_sum_complex_mass_x_complex_f,riemann_sum_complex_delta_x_complex_f
    private riemann_sum_real_mass_x_complex_f,riemann_sum_real_delta_x_complex_f
    private riemann_sum_real_mass_x_real_f,riemann_sum_real_delta_x_real_f
    
    private simpson_sum_complex_mass_x_complex_f,simpson_sum_complex_delta_x_complex_f
    private simpson_sum_real_mass_x_complex_f,simpson_sum_real_delta_x_complex_f
    private simpson_sum_real_mass_x_real_f,simpson_sum_real_delta_x_real_f
    contains
    
    complex(8) function riemann_sum_complex_mass_x_complex_f(N,x,fx) result(f)
    implicit none
        integer(4),intent(in)::N
        complex(8),intent(in)::x(N)
        complex(8),intent(in)::fx(N)
        
        integer(4) i
        
        f=0d0
        do i=1,N-1
            f=f+(x(i+1)-x(i))*fx(i)
        enddo
    endfunction riemann_sum_complex_mass_x_complex_f
    complex(8) function riemann_sum_complex_delta_x_complex_f(N,dx,fx) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        integer(4),intent(in)::N
        complex(8),intent(in)::dx
        complex(8),intent(in)::fx(N)
        
        integer(4) i
        
        if(abs(dx)<epsilon) call print_error("static_intergation.riemann_sum_complex_delta_x_complex_f","abs(dx)<epsilon")
        
        f=0d0
        do i=1,N
            f=f+dx*fx(i)
        enddo
    endfunction riemann_sum_complex_delta_x_complex_f
    
    complex(8) function riemann_sum_real_mass_x_complex_f(N,x,fx) result(f)
    implicit none
        integer(4),intent(in)::N
        real(8),intent(in)::x(N)
        complex(8),intent(in)::fx(N)
        
        integer(4) i
        
        f=0d0
        do i=1,N-1
            f=f+(x(i+1)-x(i))*fx(i)
        enddo
    endfunction riemann_sum_real_mass_x_complex_f
    complex(8) function riemann_sum_real_delta_x_complex_f(N,dx,fx) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        integer(4),intent(in)::N
        real(8),intent(in)::dx
        complex(8),intent(in)::fx(N)
        
        integer(4) i
        
        if(abs(dx)<epsilon) call print_error("static_intergation.riemann_sum_real_delta_x_complex_f","abs(dx)<epsilon")
        
        f=0d0
        do i=1,N
            f=f+dx*fx(i)
        enddo
    endfunction riemann_sum_real_delta_x_complex_f
    
    real(8) function riemann_sum_real_mass_x_real_f(N,x,fx) result(f)
    implicit none
        integer(4),intent(in)::N
        real(8),intent(in)::x(N)
        real(8),intent(in)::fx(N)
        
        integer(4) i
        
        f=0d0
        do i=1,N-1
            f=f+(x(i+1)-x(i))*fx(i)
        enddo
    endfunction riemann_sum_real_mass_x_real_f
    real(8) function riemann_sum_real_delta_x_real_f(N,dx,fx) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        integer(4),intent(in)::N
        real(8),intent(in)::dx
        real(8),intent(in)::fx(N)
        
        integer(4) i
        
        if(abs(dx)<epsilon) call print_error("static_intergation.riemann_sum_real_delta_x_real_f","abs(dx)<epsilon")
        
        f=0d0
        do i=1,N
            f=f+dx*fx(i)
        enddo
    endfunction riemann_sum_real_delta_x_real_f
    
    
    
    
    
    complex(8) function simpson_sum_complex_mass_x_complex_f(N,x,fx) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        integer(4),intent(in)::N
        complex(8),intent(in)::x(N)
        complex(8),intent(in)::fx(N)
        
        integer(4) i
        
        f=0d0
        i=1
        do while(i+2<=N)
            if(abs((x(i+1)-x(i))-(x(i+2)-x(i+1)))>=epsilon) &
                call print_error("static_integration.simpson_sum_complex_mass_x_complex_f",&
                    "abs((x(i+1)-x(i))-(x(i+2)-x(i+1)))>=epsilon   (require dx_1==dx_2)")
            
            f=f+(x(i+2)-x(i))/6d0*(fx(i)+fx(i+1)+fx(i+1)+fx(i+1)+fx(i+1)+fx(i+2))
            
            i=i+2
        enddo
    endfunction simpson_sum_complex_mass_x_complex_f
    complex(8) function simpson_sum_complex_delta_x_complex_f(N,dx,fx) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        integer(4),intent(in)::N
        complex(8),intent(in)::dx
        complex(8),intent(in)::fx(N)
        
        complex(8) dx26
        integer(4) i
        
        if(abs(dx)<epsilon) call print_error("static_intergation.simpson_sum_complex_delta_x_complex_f","abs(dx)<epsilon")
        
        dx26=dx/3d0
        
        f=0d0
        i=1
        do while(i+2<=N)
            f=f+dx26*(fx(i)+fx(i+1)+fx(i+1)+fx(i+1)+fx(i+1)+fx(i+2))
            
            i=i+2
        enddo
    endfunction simpson_sum_complex_delta_x_complex_f
    
    complex(8) function simpson_sum_real_mass_x_complex_f(N,x,fx) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        integer(4),intent(in)::N
        real(8),intent(in)::x(N)
        complex(8),intent(in)::fx(N)
        
        integer(4) i
        
        f=0d0
        i=1
        do while(i+2<=N)
            if(abs((x(i+1)-x(i))-(x(i+2)-x(i+1)))>=epsilon) &
                call print_error("static_integration.simpson_sum_real_mass_x_complex_f",&
                    "abs((x(i+1)-x(i))-(x(i+2)-x(i+1)))>=epsilon   (require dx_1==dx_2)")
            
            f=f+(x(i+2)-x(i))/6d0*(fx(i)+fx(i+1)+fx(i+1)+fx(i+1)+fx(i+1)+fx(i+2))
            
            i=i+2
        enddo
    endfunction simpson_sum_real_mass_x_complex_f
    complex(8) function simpson_sum_real_delta_x_complex_f(N,dx,fx) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        integer(4),intent(in)::N
        real(8),intent(in)::dx
        complex(8),intent(in)::fx(N)
        
        real(8) dx26
        integer(4) i
        
        if(abs(dx)<epsilon) call print_error("static_intergation.simpson_sum_real_delta_x_complex_f","abs(dx)<epsilon")
        
        dx26=dx/3d0
        
        f=0d0
        i=1
        do while(i+2<=N)
            f=f+dx26*(fx(i)+fx(i+1)+fx(i+1)+fx(i+1)+fx(i+1)+fx(i+2))
            
            i=i+2
        enddo
    endfunction simpson_sum_real_delta_x_complex_f
    
    real(8) function simpson_sum_real_mass_x_real_f(N,x,fx) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        integer(4),intent(in)::N
        real(8),intent(in)::x(N)
        real(8),intent(in)::fx(N)
        
        integer(4) i
        
        f=0d0
        i=1
        do while(i+2<=N)
            if(abs((x(i+1)-x(i))-(x(i+2)-x(i+1)))>=epsilon) &
                call print_error("static_integration.simpson_sum_real_mass_x_real_f",&
                    "abs((x(i+1)-x(i))-(x(i+2)-x(i+1)))>=epsilon   (require dx_1==dx_2)")
            
            f=f+(x(i+2)-x(i))/6d0*(fx(i)+fx(i+1)+fx(i+1)+fx(i+1)+fx(i+1)+fx(i+2))
            
            i=i+2
        enddo
    endfunction simpson_sum_real_mass_x_real_f
    real(8) function simpson_sum_real_delta_x_real_f(N,dx,fx) result(f)
    use math,only:epsilon
    use system,only:print_error
    integer(4),intent(in)::N
        real(8),intent(in)::dx
        real(8),intent(in)::fx(N)
        
        real(8) dx26
        integer(4) i
        
        if(abs(dx)<epsilon) call print_error("static_intergation.simpson_sum_real_delta_x_real_f","abs(dx)<epsilon")
        
        dx26=dx/3d0
        
        f=0d0
        i=1
        do while(i+2<=N)
            f=f+dx26*(fx(i)+fx(i+1)+fx(i+1)+fx(i+1)+fx(i+1)+fx(i+2))
            
            i=i+2
        enddo
    endfunction simpson_sum_real_delta_x_real_f
endmodule static_integration
