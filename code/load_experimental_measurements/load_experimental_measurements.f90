module load_experimental_measurements
implicit none
    public::init_load_experimental_measurements,destructor_load_experimental_measurements
    
    public::get_Nx,get_Nt,get_xi,get_tj,get_uij
    public::if_dx_const,get_dx,get_x0,if_dt_const,get_dt,get_t0
    public::get_sceptum_in_xi
    
    integer(4) Nx,Nt
    real(8),allocatable::x(:),t(:)
    real(8),allocatable::u(:,:)!xi,tj
    
    logical(1) is_dt_const,is_dx_const
    real(8) dx,dt
    real(8)::dt_max_difference_for_be_const=1d-7
    real(8)::dx_max_difference_for_be_const=1d-7
    
    private::Nx,Nt,x,t,u
    private::get_sceptum_in_xi_common,get_sceptum_in_xi_dt_const
    private::is_dt_const,is_dx_const,dt_max_difference_for_be_const,dx_max_difference_for_be_const
    contains
    
    subroutine init_load_experimental_measurements
    use smoothing_signal,only:arithmetic_mean_smoothing
    implicit none
        integer(4) file
        
        integer(4) i,j
        real(8),allocatable::signal(:)
        
        open(newunit=file,file="input/steel/_x.data")
        read(file,*),Nx
        allocate(x(Nx))
        do i=1,Nx
            read(file,*),x(i)
            x(i)=x(i)*1d2
        enddo
        close(file)
        
        open(newunit=file,file="input/steel/_t.data")
        read(file,*),Nt
        allocate(t(Nt))
        do j=1,Nt
            read(file,*),t(j)
            t(j)=t(j)*1d6
        enddo
        close(file)
        
        open(newunit=file,file="input/steel/_u.data")
        allocate(u(Nx,Nt))
        do i=1,Nx
            do j=1,Nt
                read(file,*),u(i,j)
                u(i,j)=u(i,j)*1d2
            enddo
        enddo
        close(file)
        
        is_dx_const=.true.
        dx=x(2)-x(1)
        do i=3,Nx
            if(abs(x(i)-x(i-1)-dx)>dx_max_difference_for_be_const) then
                is_dx_const=.false.
                dx=-1d0
                exit
            endif
        enddo
        
        is_dt_const=.true.
        dt=t(2)-t(1)
        do i=3,Nt
            if(abs(t(i)-t(i-1)-dt)>dt_max_difference_for_be_const) then
                is_dt_const=.false.
                dt=-1d0
                exit
            endif
        enddo
        
        allocate(signal(Nt))
        do i=1,Nx
            do j=1,Nt
                signal(j)=u(i,j)
            enddo
            call arithmetic_mean_smoothing(3,Nt,signal,.true.)
            do j=1,Nt
                u(i,j)=signal(j)
            enddo
        enddo
        deallocate(signal)
    endsubroutine init_load_experimental_measurements
    
    subroutine destructor_load_experimental_measurements
        deallocate(x)
        deallocate(t)
        deallocate(u)
    endsubroutine destructor_load_experimental_measurements
    
    integer(4) function get_Nx() result(f)
    implicit none
        f=Nx
    endfunction get_Nx
    integer(4) function get_Nt() result(f)
    implicit none
        f=Nt
    endfunction get_Nt
    real(4) function get_xi(i) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::i
        
        if(i<1.or.i>Nx) call print_error("load_experimental_measurements.get_xi","i<1.or.i>Nx")
        
        f=x(i)
    endfunction get_xi
    real(4) function get_tj(j) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::j
        
        if(j<1.or.j>Nt) call print_error("load_experimental_measurements.get_tj","j<1.or.j>Nt")
        
        f=t(j)
    endfunction get_tj
    real(4) function get_uij(i,j) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::i
        integer(4),intent(in)::j
        
        if(i<1.or.i>Nx) call print_error("load_experimental_measurements.get_uij","i<1.or.i>Nx")
        if(j<1.or.j>Nt) call print_error("load_experimental_measurements.get_uij","j<1.or.j>Nt")
        
        f=u(i,j)
    endfunction get_uij
    
    complex(8) function get_sceptum_in_xi(i,omega) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::i
        complex(8),intent(in)::omega
        
        if(i<1.or.i>Nx) call print_error("load_experimental_measurements.get_sceptum_in_xi","i<1.or.i>Nx")
        
        if(is_dt_const) then
            f=get_sceptum_in_xi_dt_const(i,omega)
        else
            f=get_sceptum_in_xi_common(i,omega)
        endif
    endfunction get_sceptum_in_xi
    complex(8) function get_sceptum_in_xi_common(i,omega) result(f)
    use math,only:ci,c0,pi
    implicit none
        integer(4),intent(in)::i
        complex(8),intent(in)::omega
        
        integer(4) tj
        
        f=u(i,1)*exp(ci*omega*t(1))*t(1)
        do tj=2,Nt
            f=f+u(i,tj)*exp(ci*omega*t(tj))*(t(tj)-t(tj-1))
        enddo
        f=f/sqrt(pi+pi)
    endfunction get_sceptum_in_xi_common
    complex(8) function get_sceptum_in_xi_dt_const(i,omega) result(f)
    use math,only:ci,c0,pi
    implicit none
        integer(4),intent(in)::i
        complex(8),intent(in)::omega
        
        integer(4) tj
        complex(8) exps,dexp
        
        exps=exp(ci*omega*t(1))
        dexp=exp(ci*omega*dt)
        
        f=u(i,1)*exps*t(1)
        do tj=2,Nt
            exps=exps*dexp
            f=f+u(i,tj)*exps*dt
        enddo
        f=f/sqrt(pi+pi)
    endfunction get_sceptum_in_xi_dt_const
    
    logical(1) function if_dt_const() result(f)
    implicit none
        f=is_dt_const
    endfunction if_dt_const
    real(8) function get_dt() result(f)
    use system,only:print_error
    implicit none
        if(.not.is_dt_const) call print_error("load_experimental_measurements.get_dt",".not.is_dt_const")
        
        f=dt
    endfunction get_dt
    real(8) function get_t0() result(f)
    use system,only:print_error
    implicit none
        if(.not.is_dt_const) call print_error("load_experimental_measurements.get_t0",".not.is_dt_const")
        
        f=t(1)
    endfunction get_t0
    logical(1) function if_dx_const() result(f)
    implicit none
        f=is_dx_const
    endfunction if_dx_const
    real(8) function get_dx() result(f)
    use system,only:print_error
    implicit none
        if(.not.is_dx_const) call print_error("load_experimental_measurements.get_dx",".not.is_dx_const")
        
        f=dx
    endfunction get_dx
    real(8) function get_x0() result(f)
    use system,only:print_error
    implicit none
        if(.not.is_dx_const) call print_error("load_experimental_measurements.get_x0",".not.is_dx_const")
        
        f=x(1)
    endfunction get_x0
endmodule load_experimental_measurements
