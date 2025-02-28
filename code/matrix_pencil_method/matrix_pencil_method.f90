module matrix_pencil_method
implicit none
    public::init_matrix_pencil_method,destructor_matrix_pencil_method
    
    public::count_dispersion_numbers
    
    integer(4)::L
    integer(4)::Nx
    complex(8)::omega
    
    complex(8),allocatable::H0(:,:),H1(:,:)
    complex(8),allocatable::H0_H0(:,:),H0Plus(:,:)
    complex(8),allocatable::H0Plus_H1(:,:)
    integer(4) res_size
    complex(8),allocatable::res(:),lambda(:),eigenvectors(:,:)
    
    complex(8),allocatable::H0_mu(:,:),H1_mu(:,:)
    complex(8),allocatable::H0_H0_mu(:,:),H0Plus_mu(:,:)
    complex(8),allocatable::H0Plus_H1_mu(:,:)
    integer(4) res_mu_size
    complex(8),allocatable::res_mu(:),mu(:)
    
    integer(4)::L_L
    complex(8),allocatable::H0_L(:,:),H1_L(:,:)
    complex(8),allocatable::H0_H0_L(:,:),H0Plus_L(:,:)
    complex(8),allocatable::H0Plus_H1_L(:,:)
    integer(4) res_L_size
    complex(8),allocatable::res_L(:),mu_L(:),eigenvectors_L(:,:)
    
    private::generate_matrixes
    
    private::count_dispersion_numbers_,count_H,count_H0plus,count_H0Plus_H1,count_dispersion_numbers_for_counted_matrixes
    private::lambda,eigenvectors,res,res_size
    
    private::count_H_for_mu,count_H0plus_for_mu,count_H0Plus_H1_for_mu,count_dispersion_numbers_for_counted_matrixes_for_mu
    private::mu
    
    private::L_L
    private::count_H_for_L,count_H0plus_for_L,count_H0Plus_H1_for_L,count_dispersion_numbers_for_counted_matrixes_for_L
    private::mu_L,eigenvectors_L
    
    private::H0,H1,H0_H0,H0Plus,H0Plus_H1
    private::H0_mu,H1_mu,H0_H0_mu,H0Plus_mu,H0Plus_H1_mu
    private::H0_L,H1_L,H0_H0_L,H0Plus_L,H0Plus_H1_L
    private::L,Nx,omega
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
        
        L_L=L+5
        
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
        call count_H0plus_for_mu()
        call count_H0Plus_H1_for_mu()
        call count_dispersion_numbers_for_counted_matrixes_for_mu()
        
        call count_H_for_L()
        call count_H0plus_for_L()
        call count_H0Plus_H1_for_L()
        call count_dispersion_numbers_for_counted_matrixes_for_L()
        
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
        if(L_L<3) call print_error("matrix_pencil_method.generate_matrixes","L_L<3")
        if(L>=Nx) call print_error("matrix_pencil_method.generate_matrixes","L>=Nx")
        if(L_L>=Nx) call print_error("matrix_pencil_method.generate_matrixes","L_L>=Nx")
        if(Nx<4) call print_error("matrix_pencil_method.generate_matrixes","Nx<4")
        
        deallocate(H0)
        deallocate(H1)
        deallocate(H0_H0)
        deallocate(H0Plus)
        deallocate(H0Plus_H1)
        deallocate(lambda)
        deallocate(res)
        deallocate(eigenvectors)
        
        deallocate(H0_mu)
        deallocate(H1_mu)
        deallocate(H0_H0_mu)
        deallocate(H0Plus_mu)
        deallocate(H0Plus_H1_mu)
        deallocate(mu)
        deallocate(res_mu)
        
        deallocate(H0_L)
        deallocate(H1_L)
        deallocate(H0_H0_L)
        deallocate(H0Plus_L)
        deallocate(H0Plus_H1_L)
        deallocate(mu_L)
        deallocate(res_L)
        deallocate(eigenvectors_L)
        
        allocate(H0(Nx-L,L))
        allocate(H1(Nx-L,L))
        allocate(H0_H0(L,L))
        allocate(H0Plus(L,Nx-L))
        allocate(H0Plus_H1(L,L))
        allocate(lambda(L))
        allocate(res(L))
        allocate(eigenvectors(L,L))
        
        allocate(H0_mu(Nx/2-L,L))
        allocate(H1_mu(Nx/2-L,L))
        allocate(H0_H0_mu(L,L))
        allocate(H0Plus_mu(L,Nx/2-L))
        allocate(H0Plus_H1_mu(L,L))
        allocate(mu(L))
        allocate(res_mu(L))
        
        allocate(H0_L(Nx-L_L,L_L))
        allocate(H1_L(Nx-L_L,L_L))
        allocate(H0_H0_L(L_L,L_L))
        allocate(H0Plus_L(L_L,Nx-L_L))
        allocate(H0Plus_H1_L(L_L,L_L))
        allocate(mu_L(L_L))
        allocate(res_L(L_L))
        allocate(eigenvectors_L(L_L,L_L))
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
            
            lambda(i)=log(lambda(i))*(-ci)/get_dx()
            
            filter_check=.false.
            do j=1,res_mu_size
                if(abs((lambda(i)-res_mu(j))/lambda(i))<0.01d0) then
                    filter_check=.true.
                    exit
                endif
            enddo
            if(.not.filter_check) cycle
            
            filter_check=.false.
            do j=1,res_L_size
                if(abs((lambda(i)-res_L(j))/lambda(i))<0.01d0) then
                    filter_check=.true.
                    exit
                endif
            enddo
            if(.not.filter_check) cycle
            
            filter_check=.false.
            filter_check=real(lambda(i))>epsilon
            if(.not.filter_check) cycle
            
            !filter_check=.false.
            !filter_check=aimag(lambda(i))<15d0
            !if(.not.filter_check) cycle
            
            filter_check=.false.
            filter_check=real(lambda(i))<pi/get_dx()
            if(.not.filter_check) cycle
            
            res_size=res_size+1
            res(res_size)=lambda(i)
        enddo
    endsubroutine count_dispersion_numbers_for_counted_matrixes
    
    
    
    
    
    subroutine count_H_for_mu()
    use load_experimental_measurements,only:get_sceptum_in_xi
    implicit none
        integer(4) i,j
        
        do i=1,Nx/2-L
            do j=1,L
                H0_mu(i,j)=get_sceptum_in_xi(2*(i+j-1)-1,omega)
                H1_mu(i,j)=get_sceptum_in_xi(2*(i+j)-1,omega)
            enddo
        enddo
    endsubroutine count_H_for_mu
    subroutine count_H0plus_for_mu()
    use matrix_complex8,only:inverse
    use math,only:c0,ci
    implicit none
        integer(4) i,j,k
        
        do i=1,L
            do j=1,L
                H0_H0_mu(i,j)=c0
                do k=1,Nx/2-L
                    H0_H0_mu(i,j)=H0_H0_mu(i,j)+conjg(H0_mu(k,i))*H0_mu(k,j)
                enddo
            enddo
            H0_H0_mu(i,i)=H0_H0_mu(i,i)+1d-9
        enddo
        
        call inverse(L,H0_H0_mu)
        do i=1,L
            do j=1,Nx/2-L
                H0Plus_mu(i,j)=c0
                do k=1,L
                    H0Plus_mu(i,j)=H0Plus_mu(i,j)+H0_H0_mu(i,k)*conjg(H0_mu(j,k))
                enddo
            enddo
        enddo
    endsubroutine count_H0plus_for_mu
    subroutine count_H0Plus_H1_for_mu()
    use matrix_complex8,only:multiply
    implicit none
        call multiply(L,Nx/2-L,L,H0Plus_mu,H1_mu,H0Plus_H1_mu)
    endsubroutine count_H0Plus_H1_for_mu
    subroutine count_dispersion_numbers_for_counted_matrixes_for_mu()
    use load_experimental_measurements,only:get_dx
    use math,only:pi,ci,epsilon
    implicit none
        integer(4) i,j
        logical(1) filter_check
        complex(8) t
        
        call Ev_evaluate(H0Plus_H1_mu,L,mu,eigenvectors)
        
        res_mu_size=0
        do i=1,L
            filter_check=abs(mu(i))>epsilon
            if(.not.filter_check) cycle
            
            mu(i)=log(mu(i))*(-ci)/get_dx()*0.5d0
            
            filter_check=.false.
            filter_check=real(mu(i))>epsilon
            if(.not.filter_check) cycle
            
            filter_check=.false.
            filter_check=real(mu(i))<pi/get_dx()*0.5d0
            if(.not.filter_check) cycle
            
            res_mu_size=res_mu_size+1
            res_mu(res_mu_size)=mu(i)
        enddo
    endsubroutine count_dispersion_numbers_for_counted_matrixes_for_mu
    
    
    
    
    
    subroutine count_H_for_L()
    use load_experimental_measurements,only:get_sceptum_in_xi
    implicit none
        integer(4) i,j
        
        do i=1,Nx-L_L
            do j=1,L_L
                H0_L(i,j)=get_sceptum_in_xi(i+j-1,omega)
                H1_L(i,j)=get_sceptum_in_xi(i+j,omega)
            enddo
        enddo
    endsubroutine count_H_for_L
    subroutine count_H0plus_for_L()
    use matrix_complex8,only:inverse
    use math,only:c0,ci
    implicit none
        integer(4) i,j,k
        
        do i=1,L_L
            do j=1,L_L
                H0_H0_L(i,j)=c0
                do k=1,Nx-L_L
                    H0_H0_L(i,j)=H0_H0_L(i,j)+conjg(H0_L(k,i))*H0_L(k,j)
                enddo
            enddo
            H0_H0_L(i,i)=H0_H0_L(i,i)+1d-9
        enddo
        
        call inverse(L_L,H0_H0_L)
        do i=1,L_L
            do j=1,Nx-L_L
                H0Plus_L(i,j)=c0
                do k=1,L_L
                    H0Plus_L(i,j)=H0Plus_L(i,j)+H0_H0_L(i,k)*conjg(H0_L(j,k))
                enddo
            enddo
        enddo
    endsubroutine count_H0plus_for_L
    subroutine count_H0Plus_H1_for_L()
    use matrix_complex8,only:multiply
    implicit none
        call multiply(L_L,Nx-L_L,L_L,H0Plus_L,H1_L,H0Plus_H1_L)
    endsubroutine count_H0Plus_H1_for_L
    subroutine count_dispersion_numbers_for_counted_matrixes_for_L()
    use load_experimental_measurements,only:get_dx
    use math,only:pi,ci,epsilon
    implicit none
        integer(4) i,j
        logical(1) filter_check
        complex(8) t
        
        call Ev_evaluate(H0Plus_H1_L,L_L,mu_L,eigenvectors_L)
        
        res_L_size=0
        do i=1,L_L
            filter_check=abs(mu_L(i))>epsilon
            if(.not.filter_check) cycle
            
            mu_L(i)=log(mu_L(i))*(-ci)/get_dx()
            
            filter_check=.false.
            filter_check=real(mu_L(i))>epsilon
            if(.not.filter_check) cycle
            
            filter_check=.false.
            filter_check=real(mu_L(i))<pi/get_dx()
            if(.not.filter_check) cycle
            
            res_L_size=res_L_size+1
            res_L(res_L_size)=mu_L(i)
        enddo
    endsubroutine count_dispersion_numbers_for_counted_matrixes_for_L
    
    
    
    
    
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
        allocate(eigenvectors(1,1))
        
        res_mu_size=-1
        allocate(H0_mu(1,1))
        allocate(H1_mu(1,1))
        allocate(H0_H0_mu(1,1))
        allocate(H0Plus_mu(1,1))
        allocate(H0Plus_H1_mu(1,1))
        allocate(mu(1))
        allocate(res_mu(1))
        
        L_L=-1
        res_L_size=-1
        allocate(H0_L(1,1))
        allocate(H1_L(1,1))
        allocate(H0_H0_L(1,1))
        allocate(H0Plus_L(1,1))
        allocate(H0Plus_H1_L(1,1))
        allocate(mu_L(1))
        allocate(res_L(1))
        allocate(eigenvectors_L(1,1))
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
        deallocate(eigenvectors)
        
        deallocate(H0_mu)
        deallocate(H1_mu)
        deallocate(H0_H0_mu)
        deallocate(H0Plus_mu)
        deallocate(H0Plus_H1_mu)
        deallocate(mu)
        deallocate(res_mu)
        
        deallocate(H0_L)
        deallocate(H1_L)
        deallocate(H0_H0_L)
        deallocate(H0Plus_L)
        deallocate(H0Plus_H1_L)
        deallocate(mu_L)
        deallocate(res_L)
        deallocate(eigenvectors_L)
    endsubroutine destructor_matrix_pencil_method
endmodule matrix_pencil_method
