module sigma_and_eigenvectors
implicit none
    public::get_sigma,get_eigenvector_component
    
    public::init_sigma_and_eigenvectors,destructor_sigma_and_eigenvectors
    
    complex(8),allocatable::sigma(:,:)!number_of_sigma,number_of_layer
    complex(8),allocatable::m(:,:,:)!number_of_sigma,number_of_index,number_of_layer
    
    complex(8),save::last_alpha=(1.182745372d99,1.91764672d101),last_beta=(1.183645274d98,1.8245282591d97)
    
    private::last_alpha,last_beta,sigma,m,count_sigma_and_eigenvectors
    contains
    
    complex(8) function get_sigma(alpha,beta,number_of_sigma,number_of_layer) result(f)
    use system,only:print_error
    use main_parameters,only:number_of_layers,get_anisotropic
    implicit none
        complex(8),intent(in)::alpha
        complex(8),intent(in)::beta
        integer(4),intent(in)::number_of_sigma
        integer(4),intent(in)::number_of_layer
        
        if(number_of_sigma<1.or.number_of_sigma>6) then
            call print_error("sigma_and_eigenvectors.get_sigma","number_of_sigma<1.or.number_of_sigma>6")
        endif
        if(number_of_layer<1.or.number_of_layer>number_of_layers) then
            call print_error("sigma_and_eigenvectors.get_sigma",&
                "number_of_layer<1.or.number_of_layer>number_of_layers")
        endif
        
        !if(get_anisotropic()==0.and..not.beta==0d0) call print_error("sigma_and_eigenvectors.get_sigma",&
        !    "(anisotropic==0.and..not.bera==0d0")
        
        call count_sigma_and_eigenvectors(alpha,beta)
        f=sigma(number_of_sigma,number_of_layer)
    endfunction get_sigma
    complex(8) function get_eigenvector_component(alpha,beta,number_of_sigma,index,number_of_layer) result(f)
    use system,only:print_error
    use main_parameters,only:number_of_layers,get_anisotropic
    implicit none
        complex(8),intent(in)::alpha
        complex(8),intent(in)::beta
        integer(4),intent(in)::number_of_sigma
        integer(4),intent(in)::index
        integer(4),intent(in)::number_of_layer
        
        if(number_of_sigma<1.or.number_of_sigma>6) then
            call print_error("sigma_and_eigenvectors.get_eigenvector_component","number_of_sigma<1.or.number_of_sigma>6")
        endif
        if(number_of_layer<1.or.number_of_layer>number_of_layers) then
            call print_error("sigma_and_eigenvectors.get_eigenvector_component",&
                "number_of_layer<1.or.number_of_layer>number_of_layers")
        endif
        if(index<1.or.index>3) call print_error("sigma_and_eigenvectors.get_eigenvector_component","index<1.or.index>3")
        
        !if(get_anisotropic()==0.and..not.beta==0d0) call print_error("sigma_and_eigenvectors.get_eigenvector_component",&
        !    "(anisotropic==0.and..not.bera==0d0")
        
        call count_sigma_and_eigenvectors(alpha,beta)
        f=m(number_of_sigma,index,number_of_layer)
    endfunction get_eigenvector_component
    
    subroutine count_sigma_and_eigenvectors(alpha,beta)
    use main_parameters,only:number_of_layers,Cijkl,rho,omega,get_anisotropic,lambda,mu
    use system,only:print_error
    use math,only:ci,c0,epsilon
    implicit none
        complex(8),intent(in)::alpha
        complex(8),intent(in)::beta
        
        integer(4) i,j,k,l,layer
        complex(8) A(6,6),B(3,3,0:2),sigma_(6),m_(6,6),q(6)
        complex(8) temp
        integer(4) is_more
        logical(1) all_sigma_is_real_aimag
        
        complex(8) alpha2
        
        !if(get_anisotropic()==0.and..not.beta==0d0) call print_error("sigma_and_eigenvectors.count_sigma_and_eigenvectors",&
        !    "(anisotropic==0.and..not.bera==0d0")
        
        if(alpha==last_alpha.and.beta==last_beta) then
            return
        else
            last_alpha=alpha
            last_beta=beta
        endif
        
        if(get_anisotropic()==0) then
            alpha2=alpha*alpha+beta*beta
            !complex(8),allocatable::sigma(:,:)!number_of_sigma,number_of_layer
            !complex(8),allocatable::m(:,:,:)!number_of_sigma,number_of_index,number_of_layer
            do layer=1,number_of_layers
                sigma(1,layer)=sqrt(alpha2-rho(layer)*omega*omega/(lambda(layer)+mu(layer)+mu(layer)))
                sigma(2,layer)=-sigma(1,layer)
                sigma(3,layer)=sqrt(alpha2-rho(layer)*omega*omega/mu(layer))
                sigma(4,layer)=-sigma(3,layer)
                sigma(5,layer)=sigma(3,layer)
                sigma(6,layer)=-sigma(3,layer)
                
                m(1,1,layer)=1d0
                !m(1,2,layer)=sigma(1,layer)
                m(1,2,layer)=0d0
                !m(1,4,layer)=0d0
                m(1,3,layer)=sigma(1,layer)
                !m(1,6,layer)=sigma(1,layer)*sigma(1,layer)
                
                m(2,1,layer)=1d0
                !m(2,2,layer)=-sigma(1,layer)
                m(2,2,layer)=0d0
                !m(2,4,layer)=0d0
                m(2,3,layer)=-sigma(1,layer)
                !m(2,6,layer)=sigma(1,layer)*sigma(1,layer)
                
                m(3,1,layer)=sigma(3,layer)
                !m(3,2,layer)=sigma(3,layer)*sigma(3,layer)
                m(3,2,layer)=0d0
                !m(3,4,layer)=0d0
                m(3,3,layer)=alpha2
                !m(3,6,layer)=alpha2*sigma(3,layer)
                
                m(4,1,layer)=-sigma(3,layer)
                !m(4,2,layer)=sigma(3,layer)*sigma(3,layer)
                m(4,2,layer)=0d0
                !m(4,4,layer)=0d0
                m(4,3,layer)=alpha2
                !m(4,6,layer)=-alpha2*sigma(3,layer)
                
                m(5,1,layer)=0d0
                !m(5,2,layer)=0d0
                m(5,2,layer)=1d0
                !m(5,4,layer)=sigma(3,layer)
                m(5,3,layer)=0d0
                !m(5,6,layer)=0d0
                
                m(6,1,layer)=0d0
                !m(6,2,layer)=0d0
                m(6,2,layer)=1d0
                !m(6,4,layer)=-sigma(3,layer)
                m(6,3,layer)=0d0
                !m(6,6,layer)=0d0
            enddo
            
            return
        endif
        
        do layer=1,number_of_layers
            do i=1,3
                do l=1,3
                    B(i,l,0)=&
                        Cijkl(i,1,1,l,layer)*alpha*alpha+&
                        Cijkl(i,1,2,l,layer)*alpha*beta+&
                        Cijkl(i,2,1,l,layer)*beta*alpha+&
                        Cijkl(i,2,2,l,layer)*beta*beta
                    B(i,l,1)=&
                        Cijkl(i,1,3,l,layer)*alpha+&
                        Cijkl(i,2,3,l,layer)*beta+&
                        Cijkl(i,3,1,l,layer)*alpha+&
                        Cijkl(i,3,2,l,layer)*beta
                    B(i,l,2)=&
                        Cijkl(i,3,3,l,layer)
                enddo
            enddo
            
            do i=1,3
                do j=1,3
                    if(.not.(i==j.and.abs(B(i,j,2))>epsilon.or.i/=j.and.abs(B(i,j,2))==c0)) then
                        call print_error("sigma_and_eigenvectors.count_sigma_with_eigenvectors",&
                            "doto common case for preparating B(i,l,2) before generating A")
                    endif
                enddo
            enddo
            
            do i=1,3
                B(i,i,0)=B(i,i,0)-rho(layer)*omega*omega
            enddo
            do i=1,3
                do j=1,3
                    B(i,j,0)=B(i,j,0)/B(i,i,2)
                    B(i,j,1)=B(i,j,1)/B(i,i,2)*ci
                enddo
                B(i,i,2)=1d0
            enddo
            
            !U_1 U_1'_z U_2 U_2'_z U_3 U_3'_z
            A(1,1)=0d0;      A(1,2)=1d0;      A(1,3)=0d0;      A(1,4)=0d0;      A(1,5)=0d0;      A(1,6)=0d0;
            A(2,1)=B(1,1,0); A(2,2)=B(1,1,1); A(2,3)=B(1,2,0); A(2,4)=B(1,2,1); A(2,5)=B(1,3,0); A(2,6)=B(1,3,1);
            A(3,1)=0d0;      A(3,2)=0d0;      A(3,3)=0d0;      A(3,4)=1d0;      A(3,5)=0d0;      A(3,6)=0d0;
            A(4,1)=B(2,1,0); A(4,2)=B(2,1,1); A(4,3)=B(2,2,0); A(4,4)=B(2,2,1); A(4,5)=B(2,3,0); A(4,6)=B(2,3,1);
            A(5,1)=0d0;      A(5,2)=0d0;      A(5,3)=0d0;      A(5,4)=0d0;      A(5,5)=0d0;      A(5,6)=1d0;
            A(6,1)=B(3,1,0); A(6,2)=B(3,1,1); A(6,3)=B(3,2,0); A(6,4)=B(3,2,1); A(6,5)=B(3,3,0); A(6,6)=B(3,3,1);
            
            call Ev_evaluate(A,6,sigma_,m_)
            
            all_sigma_is_real_aimag=.true.
            do i=1,6
                sigma(i,layer)=sigma_(i)
                if(real(sigma(i,layer))>=epsilon.and.aimag(sigma(i,layer))>=epsilon) all_sigma_is_real_aimag=.false.
                do j=1,3
                    m(i,j,layer)=m_(j+j-1,i)
                enddo
            enddo
            
            do i=1,6
                do j=i+1,6
                    if(all_sigma_is_real_aimag) then
                        is_more=0
                        
                        if(real(sigma(i,layer))>=epsilon) is_more=is_more+4
                        if(aimag(sigma(i,layer))<=-epsilon) is_more=is_more+3
                        if(real(sigma(i,layer))<=-epsilon) is_more=is_more+2
                        if(aimag(sigma(i,layer))>=epsilon) is_more=is_more+1
                        
                        if(real(sigma(j,layer))>=epsilon) is_more=is_more-4
                        if(aimag(sigma(j,layer))<=-epsilon) is_more=is_more-3
                        if(real(sigma(j,layer))<=-epsilon) is_more=is_more-2
                        if(aimag(sigma(j,layer))>=epsilon) is_more=is_more-1
                        
                        if(is_more>0.or.is_more==0.and..not.(&
                            real(sigma(i,layer))>=epsilon.and.real(sigma(i,layer))<real(sigma(j,layer)).or.&
                            aimag(sigma(i,layer))<=-epsilon.and.aimag(sigma(i,layer))<aimag(sigma(j,layer)).or.&
                            real(sigma(i,layer))<=-epsilon.and.real(sigma(i,layer))>real(sigma(j,layer)).or.&
                            aimag(sigma(i,layer))>=epsilon.and.aimag(sigma(i,layer))>aimag(sigma(j,layer))&
                            )) cycle
                    else
                        if(.not.(&
                            abs(real(sigma(i,layer)))>=epsilon.and.abs(real(sigma(j,layer)))>=epsilon.and.&
                                real(sigma(i,layer))<real(sigma(j,layer)).or.&
                            abs(real(sigma(i,layer)))<epsilon.and.abs(real(sigma(j,layer)))>=epsilon.or.&
                            abs(real(sigma(i,layer)))<epsilon.and.abs(real(sigma(j,layer)))<epsilon.and.&
                                aimag(sigma(i,layer))>aimag(sigma(j,layer))&
                            )) cycle
                    endif
                    
                    temp=sigma(i,layer)
                    sigma(i,layer)=sigma(j,layer)
                    sigma(j,layer)=temp
                    
                    do k=1,3
                        temp=m(i,k,layer)
                        m(i,k,layer)=m(j,k,layer)
                        m(j,k,layer)=temp
                    enddo
                enddo
            enddo
            
            if(all_sigma_is_real_aimag) then
                do i=1,3
                    if(.not.real(sigma(i,layer))<=-epsilon) cycle
                    do j=4,6
                        if(.not.real(sigma(j,layer))>=-epsilon) cycle
                        
                        temp=sigma(i,layer)
                        sigma(i,layer)=sigma(j,layer)
                        sigma(j,layer)=temp
                        
                        do k=1,3
                            temp=m(i,k,layer)
                            m(i,k,layer)=m(j,k,layer)
                            m(j,k,layer)=temp
                        enddo
                    enddo
                enddo
            endif
        enddo
    endsubroutine count_sigma_and_eigenvectors
    
    subroutine init_sigma_and_eigenvectors()
    use main_parameters,only:number_of_layers
    implicit none
        allocate(sigma(6,number_of_layers))
        allocate(m(6,3,number_of_layers))
    endsubroutine init_sigma_and_eigenvectors
    subroutine destructor_sigma_and_eigenvectors()
    implicit none
        deallocate(sigma)
        deallocate(m)
    endsubroutine destructor_sigma_and_eigenvectors
endmodule sigma_and_eigenvectors
