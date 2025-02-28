module matrix_pencil_method
implicit none
    public::init_matrix_pencil_method,destructor_matrix_pencil_method
    
    public::count_dispersion_numbers

    complex(8),allocatable::H0(:,:),H1(:,:)
    complex(8),allocatable::H0_H0(:,:),H0Plus(:,:)
    complex(8),allocatable::H0Plus_H1(:,:)
    integer(4)::L
    integer(4)::Nx
    complex(8)::omega
    integer(4) res_size
    complex(8),allocatable::res(:),mu(:),lambda(:),lambda_eigenvectors(:,:)
    
    private::count_dispersion_numbers_,generate_matrixes,count_H,count_H0plus,count_H0Plus_H1,count_dispersion_numbers_for_counted_matrixes
    private::count_H_for_mu,lambda,count_dispersion_numbers_for_counted_matrixes_for_mu,lambda_eigenvectors
    
    private::H0,H1,H0_H0,H0Plus,H0Plus_H1
    private::L,Nx,omega,res,res_size
    contains
    
    subroutine count_dispersion_numbers(omega,L,res,res_size)
    use load_experimental_measurements,only:get_Nx
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::omega
        integer(4),intent(in)::L
        complex(8),intent(out)::res(L)
        integer(4),intent(out)::res_size
        
        if(real(omega)<epsilon) call print_error("matrix_pencil_method.count_dispersion_numbers","real(omega)<epsilon")
        if(L<3) call print_error("matrix_pencil_method.count_dispersion_numbers","L<3")
        if(L>=get_Nx()) call print_error("matrix_pencil_method.count_dispersion_numbers","L>=get_Nx()")
        if(get_Nx()<4) call print_error("matrix_pencil_method.count_dispersion_numbers","get_Nx()<4")
        
        call count_dispersion_numbers_(omega,L,res,res_size)
    endsubroutine count_dispersion_numbers
    subroutine count_dispersion_numbers_(omega_,L_,res_,res_size_)
    use load_experimental_measurements,only:get_Nx
    implicit none
        complex(8),intent(in)::omega_
        integer(4),intent(in)::L_
        complex(8),intent(out)::res_(L)
        integer(4),intent(out)::res_size_
        
        integer(4) i
        
        if(.not.(L==L_.and.Nx==get_Nx())) then
            L=L_
            omega=omega_
            Nx=get_Nx()
            call generate_matrixes()
        elseif(.not.omega==omega_) then
            omega=omega_
        else
            do i=1,L
                res_(i)=res(i)    
            enddo
            return
        endif
        
        call count_H_for_mu()
        call count_H0plus()
        call count_H0Plus_H1()
        call count_dispersion_numbers_for_counted_matrixes_for_mu()
        
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
        if(L<3) call print_error("matrix_pencil_method.generate_matrixes","L<3")
        if(L>=Nx) call print_error("matrix_pencil_method.generate_matrixes","L>=Nx")
        if(Nx<4) call print_error("matrix_pencil_method.generate_matrixes","Nx<4")
        
        deallocate(H0)
        deallocate(H1)
        deallocate(H0_H0)
        deallocate(H0Plus)
        deallocate(H0Plus_H1)
        deallocate(lambda)
        deallocate(res)
        deallocate(lambda_eigenvectors)
        deallocate(mu)
        
        allocate(H0(Nx-L,L))
        allocate(H1(Nx-L,L))
        allocate(H0_H0(L,L))
        allocate(H0Plus(L,Nx-L))
        allocate(H0Plus_H1(L,L))
        allocate(lambda(L))
        allocate(res(L))
        allocate(lambda_eigenvectors(L,L))
        allocate(mu(L))
    endsubroutine generate_matrixes
    
    subroutine count_H()
    use load_experimental_measurements,only:get_sceptum_in_xi
    implicit none
        integer(4) i,j
        
        do i=1,Nx-L
            do j=1,L
                H0(i,j)=get_sceptum_in_xi(i+j-1,omega)
                H1(i,j)=get_sceptum_in_xi(i+j,omega)
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
                do k=1,Nx-L
                    H0_H0(i,j)=H0_H0(i,j)+conjg(H0(k,i))*H0(k,j)
                enddo
            enddo
            H0_H0(i,i)=H0_H0(i,i)+1d-9
        enddo
        
        call inverse(L,H0_H0)
        do i=1,L
            do j=1,Nx-L
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
        call multiply(L,Nx-L,L,H0Plus,H1,H0Plus_H1)
    endsubroutine count_H0Plus_H1
    subroutine count_dispersion_numbers_for_counted_matrixes()
    use load_experimental_measurements,only:get_dx
    use math,only:ci,epsilon
    implicit none
        integer(4) i,j
        logical(1) filter_check
        complex(8) t
        
        call Ev_evaluate(H0Plus_H1,L,lambda,lambda_eigenvectors)
        
        res_size=0
        do i=1,L
            filter_check=abs(lambda(i))>epsilon
            if(.not.filter_check) cycle
            
            filter_check=.false.
            do j=1,L
                if(abs((lambda(i)-mu(i))/lambda(i))<1d0) then
                    filter_check=.true.
                    exit
                endif
            enddo
            if(.not.filter_check) cycle
            
            lambda(i)=log(lambda(i))*(-ci)/get_dx()
            
            filter_check=.false.
            filter_check=real(lambda(i))>epsilon
            if(.not.filter_check) cycle
            
            res_size=res_size+1
            res(res_size)=lambda(i)
        enddo
    endsubroutine count_dispersion_numbers_for_counted_matrixes
    
    subroutine count_H_for_mu()
    use load_experimental_measurements,only:get_sceptum_in_xi
    implicit none
        integer(4) i,j
        
        do i=1,Nx-L
            do j=1,L
                H0(i,j)=get_sceptum_in_xi(i+j,omega)
                H1(i,j)=get_sceptum_in_xi(i+j-1,omega)
            enddo
        enddo
    endsubroutine count_H_for_mu
    subroutine count_dispersion_numbers_for_counted_matrixes_for_mu()
    use math,only:epsilon
    implicit none
        integer(4) i    
        
        call Ev_evaluate(H0Plus_H1,L,mu,lambda_eigenvectors)
        
        do i=1,L
            if(abs(mu(i))<epsilon) then
                mu(i)=1d100
                cycle
            endif
            
            mu(i)=1d0/mu(i)
        enddo
    endsubroutine count_dispersion_numbers_for_counted_matrixes_for_mu
    
    subroutine init_matrix_pencil_method()
    implicit none
        L=-1
        Nx=-1
        omega=-(1d0,0d0)
        res_size=-1
        allocate(H0(1,1))
        allocate(H1(1,1))
        allocate(H0_H0(1,1))
        allocate(H0Plus(1,1))
        allocate(H0Plus_H1(1,1))
        allocate(lambda(1))
        allocate(res(1))
        allocate(lambda_eigenvectors(1,1))
        allocate(mu(1))
    endsubroutine init_matrix_pencil_method
    subroutine destructor_matrix_pencil_method
    implicit none
        deallocate(H0)
        deallocate(H1)
        deallocate(H0_H0)
        deallocate(H0Plus)
        deallocate(H0Plus_H1)
        deallocate(lambda)
        deallocate(res)
        deallocate(lambda_eigenvectors)
        deallocate(mu)
    endsubroutine destructor_matrix_pencil_method
endmodule matrix_pencil_method
