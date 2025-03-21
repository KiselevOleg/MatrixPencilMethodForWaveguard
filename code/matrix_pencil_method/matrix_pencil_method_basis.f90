module matrix_pencil_method_basis
implicit none
    public::init_matrix_pencil_method_basis,destructor_matrix_pencil_method_basis
    
    public::count_dispersion_numbers
    
    integer(4)::L
    integer(4)::Nx
    integer(4)::Nx_step
    complex(8)::omega
    
    complex(8),allocatable::H0(:,:),H1(:,:)
    complex(8),allocatable::H0_H0(:,:),H0Plus(:,:)
    complex(8),allocatable::H0Plus_H1(:,:)
    integer(4) res_size
    complex(8),allocatable::res(:),lambda(:),eigenvectors(:,:)
    
    private::generate_matrixes
    
    private::count_dispersion_numbers_,count_H,count_H0plus,count_H0Plus_H1,count_dispersion_numbers_for_counted_matrixes
    private::lambda,eigenvectors,res,res_size
    
    private::H0,H1,H0_H0,H0Plus,H0Plus_H1
    private::L,Nx,Nx_step,omega
    contains
    
    subroutine count_dispersion_numbers(omega,Nx_step,L,res,res_size)
    use load_experimental_measurements,only:get_Nx
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::omega
        integer(4),intent(in)::Nx_step
        integer(4),intent(in)::L
        complex(8),intent(out)::res(L)
        integer(4),intent(out)::res_size
        
        if(real(omega)<epsilon) call print_error("matrix_pencil_method_basis.count_dispersion_numbers","real(omega)<epsilon")
        if(Nx_step<1) call print_error("matrix_pencil_method_basis.count_dispersion_numbers","Nx_step<1")
        if(L<3) call print_error("matrix_pencil_method_basis.count_dispersion_numbers","L<3")
        if(L>=get_Nx()) call print_error("matrix_pencil_method_basis.count_dispersion_numbers","L>=get_Nx()")
        if(get_Nx()<4) call print_error("matrix_pencil_method_basis.count_dispersion_numbers","get_Nx()<4")
        
        call count_dispersion_numbers_(omega,Nx_step,L,res,res_size)
    endsubroutine count_dispersion_numbers
    subroutine count_dispersion_numbers_(omega_,Nx_step_,L_,res_,res_size_)
    use load_experimental_measurements,only:get_Nx
    implicit none
        complex(8),intent(in)::omega_
        integer(4),intent(in)::Nx_step_
        integer(4),intent(in)::L_
        complex(8),intent(out)::res_(L_)
        integer(4),intent(out)::res_size_
        
        integer(4) i
        
        if(.not.(L==L_.and.Nx==get_Nx().and.Nx_step==Nx_step_)) then
            L=L_
            omega=omega_
            Nx=get_Nx()
            Nx_step=Nx_step_
            call generate_matrixes()
        elseif(.not.omega==omega_) then
            omega=omega_
        else
            do i=1,L
                res_(i)=res(i)    
            enddo
            res_size_=res_size
            return
        endif
        
        call count_H()
        call count_H0plus()
        call count_H0Plus_H1()
        call count_dispersion_numbers_for_counted_matrixes()
        
        do i=1,res_size
            res_(i)=res(i)
        enddo
        res_size_=res_size
    endsubroutine count_dispersion_numbers_
    
    subroutine generate_matrixes()
    use system,only:print_error
    implicit none
        if(Nx_step<1) call print_error("matrix_pencil_method_basis.generate_matrixes","Nx_step<1")
        if(L<3) call print_error("matrix_pencil_method_basis.generate_matrixes","L<3")
        if(L>=Nx) call print_error("matrix_pencil_method_basis.generate_matrixes","L>=Nx")
        if(Nx<4) call print_error("matrix_pencil_method_basis.generate_matrixes","Nx<4")
        
        deallocate(H0)
        deallocate(H1)
        deallocate(H0_H0)
        deallocate(H0Plus)
        deallocate(H0Plus_H1)
        deallocate(lambda)
        deallocate(res)
        deallocate(eigenvectors)
        
        allocate(H0(Nx/Nx_step-L,L))
        allocate(H1(Nx/Nx_step-L,L))
        allocate(H0_H0(L,L))
        allocate(H0Plus(L,Nx/Nx_step-L))
        allocate(H0Plus_H1(L,L))
        allocate(lambda(L))
        allocate(res(L))
        allocate(eigenvectors(L,L))
    endsubroutine generate_matrixes
    
    
    
    
    
    subroutine count_H()
    use load_experimental_measurements,only:get_sceptum_in_xi
    implicit none
        integer(4) i,j
        
        do i=1,Nx/Nx_step-L
            do j=1,L
                H0(i,j)=get_sceptum_in_xi(Nx_step*(i+j-1)-Nx_step+1,omega)
                H1(i,j)=get_sceptum_in_xi(Nx_step*(i+j)-Nx_step+1,omega)
            enddo
        enddo
    endsubroutine count_H
    subroutine count_H0plus()
    use matrix_complex8,only:inverse
    use math,only:c0,ci
    implicit none
        integer(4) i,j,k
        
        do i=1,L
            do j=1,L
                H0_H0(i,j)=c0
                do k=1,Nx/Nx_step-L
                    H0_H0(i,j)=H0_H0(i,j)+conjg(H0(k,i))*H0(k,j)
                enddo
            enddo
            H0_H0(i,i)=H0_H0(i,i)+1d-9
        enddo
        
        call inverse(L,H0_H0)
        do i=1,L
            do j=1,Nx/Nx_step-L
                H0Plus(i,j)=c0
                do k=1,L
                    H0Plus(i,j)=H0Plus(i,j)+H0_H0(i,k)*conjg(H0(j,k))
                enddo
            enddo
        enddo
    endsubroutine count_H0plus
    subroutine count_H0Plus_H1()
    use matrix_complex8,only:multiply
    implicit none
        call multiply(L,Nx/Nx_step-L,L,H0Plus,H1,H0Plus_H1)
    endsubroutine count_H0Plus_H1
    subroutine count_dispersion_numbers_for_counted_matrixes()
    use load_experimental_measurements,only:get_dx
    use math,only:pi,ci,epsilon
    implicit none
        integer(4) i,j
        logical(1) filter_check
        complex(8) t
        
        call Ev_evaluate(H0Plus_H1,L,lambda,eigenvectors)
        
        res_size=0
        do i=1,L
            filter_check=abs(lambda(i))>epsilon
            if(.not.filter_check) cycle
            
            lambda(i)=log(lambda(i))*(-ci)/get_dx()/Nx_step
            
            filter_check=.false.
            filter_check=real(lambda(i))>epsilon
            if(.not.filter_check) cycle
            
            !filter_check=.false.
            !filter_check=aimag(lambda(i))<15d0
            !if(.not.filter_check) cycle
            
            filter_check=.false.
            filter_check=real(lambda(i))<pi/get_dx()/Nx_step
            if(.not.filter_check) cycle
            
            res_size=res_size+1
            res(res_size)=lambda(i)
        enddo
    endsubroutine count_dispersion_numbers_for_counted_matrixes
    
    
    
    
    
    subroutine init_matrix_pencil_method_basis()
    implicit none
        L=-1
        Nx=-1
        Nx_step=-1
        omega=-(1d0,0d0)
        res_size=-1
        allocate(H0(1,1))
        allocate(H1(1,1))
        allocate(H0_H0(1,1))
        allocate(H0Plus(1,1))
        allocate(H0Plus_H1(1,1))
        allocate(lambda(1))
        allocate(res(1))
        allocate(eigenvectors(1,1))
    endsubroutine init_matrix_pencil_method_basis
    subroutine destructor_matrix_pencil_method_basis
    implicit none
        deallocate(H0)
        deallocate(H1)
        deallocate(H0_H0)
        deallocate(H0Plus)
        deallocate(H0Plus_H1)
        deallocate(lambda)
        deallocate(res)
        deallocate(eigenvectors)
    endsubroutine destructor_matrix_pencil_method_basis
endmodule matrix_pencil_method_basis
