module matrix_pencil_method
implicit none
    public::init_matrix_pencil_method,destructor_matrix_pencil_method

    complex(8),allocatable::H0(:,:),H1(:,:)
    complex(8),allocatable::H0_H0(:,:),H0Plus(:,:)
    complex(8),allocatable::H0Plus_H1(:,:)
    integer(4)::L
    integer(4)::Nx
    complex(8)::omega
    complex(8),allocatable::res(:),res_eigenvectors(:,:)
    
    private::count_dispersion_numbers_,generate_matrixes,count_H,count_H0plus,count_H0Plus_H1,count_dispersion_numbers_for_counted_matrixes
    
    private::H0,H1,H0_H0,H0Plus,H0Plus_H1
    private::L,Nx,omega,res,res_eigenvectors
    contains
    
    subroutine count_dispersion_numbers(omega,L,res)
    use load_experimental_measurements,only:get_Nx
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::omega
        integer(4),intent(in)::L
        complex(8),intent(out)::res(L)
        
        if(real(omega)<epsilon) call print_error("matrix_pencil_method.count_dispersion_numbers","real(omega)<epsilon")
        if(L<3) call print_error("matrix_pencil_method.count_dispersion_numbers","L<3")
        if(L>=get_Nx()) call print_error("matrix_pencil_method.count_dispersion_numbers","L>=get_Nx()")
        if(get_Nx()<4) call print_error("matrix_pencil_method.count_dispersion_numbers","get_Nx()<4")
        
        call count_dispersion_numbers_(omega,L,res)
    endsubroutine count_dispersion_numbers
    subroutine count_dispersion_numbers_(omega_,L_,res_)
    use load_experimental_measurements,only:get_Nx
    implicit none
        complex(8),intent(in)::omega_
        integer(4),intent(in)::L_
        complex(8),intent(out)::res_(L)
        
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
        
        call count_H()
        call count_H0plus()
        call count_H0Plus_H1()
        call count_dispersion_numbers_for_counted_matrixes()
        
        do i=1,L
            res_(i)=res(i)
        enddo
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
        deallocate(res)
        deallocate(res_eigenvectors)
        
        allocate(H0(Nx-L,L))
        allocate(H1(Nx-L,L))
        allocate(H0_H0(L,L))
        allocate(H0Plus(L,Nx-L))
        allocate(H0Plus_H1(L,L))
        allocate(res(L))
        allocate(res_eigenvectors(L,L))
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
    implicit none
        call Ev_evaluate(H0Plus_H1,L,res,res_eigenvectors)
    endsubroutine count_dispersion_numbers_for_counted_matrixes
    
    subroutine init_matrix_pencil_method()
    implicit none
        L=-1
        Nx=-1
        omega=-(1d0,0d0)
        allocate(H0(1,1))
        allocate(H1(1,1))
        allocate(H0_H0(1,1))
        allocate(H0Plus(1,1))
        allocate(H0Plus_H1(1,1))
        allocate(res_eigenvectors(1,1))
    endsubroutine init_matrix_pencil_method
    subroutine destructor_matrix_pencil_method
    implicit none
        deallocate(H0)
        deallocate(H1)
        deallocate(H0_H0)
        deallocate(H0Plus)
        deallocate(H0Plus_H1)
        deallocate(res_eigenvectors)
    endsubroutine destructor_matrix_pencil_method
endmodule matrix_pencil_method
