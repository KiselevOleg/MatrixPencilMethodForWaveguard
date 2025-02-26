module load_experimental_measurements
implicit none
    public::init_load_experimental_measurements,destructor_load_experimental_measurements
    
    public::get_Nx,get_Nt,get_xi,get_tj,get_uij
    
    integer(4)::Nx,Nt
    real(8),allocatable::x(:),t(:)
    real(8),allocatable::u(:,:)!xi,tj
    
    private::Nx,Nt,x,t,u
    contains
    subroutine init_load_experimental_measurements
        integer(4) file
        
        integer(4) i,j
        
        open(newunit=file,file="input/steel/_x.data")
        read(file,*),Nx
        allocate(x(Nx))
        do i=1,Nx
            read(file,*),x(i)
        enddo
        close(file)
        
        open(newunit=file,file="input/steel/_t.data")
        read(file,*),Nt
        allocate(t(Nt))
        do j=1,Nt
            read(file,*),t(j)
        enddo
        close(file)
        
        open(newunit=file,file="input/steel/_u.data")
        allocate(u(Nx,Nt))
        do i=1,Nx
            do j=1,Nt
                read(file,*),u(i,j)
            enddo
        enddo
        close(file)
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
endmodule load_experimental_measurements
