module count_K
implicit none
    public::init_count_K,destructor_count_K
    
    public::K
    public::get_det_A,get_det_AX,get_det_AY
    
    logical(1),save::already_counted(3)!number en for an elementary Q
    
    !anisotropic
    complex(8),allocatable::t(:,:,:)!number_of_t_for_layer,number_of_layers,number_of_en
    complex(8),save::S(3,6,2)!i,j,layer(1~up,2~down)
    complex(8),allocatable::C(:,:,:,:)!i,j,position(1~left,2~right),number_of_up_layer_for_border
    complex(8),allocatable::F(:)
    !isotropic
    complex(8),allocatable::Yt(:,:,:)!number_of_t_for_layer,number_of_layers,number_of_en
    complex(8),save::YS(2,4,2)!i,j,layer(1~up,2~down)
    complex(8),allocatable::YC(:,:,:,:)!i,j,position(1~left,2~right),number_of_up_layer_for_border
    complex(8),allocatable::YF(:)
    complex(8),allocatable::Xt(:,:,:)!number_of_t_for_layer,number_of_layers,number_of_en
    complex(8),save::XS(1,2,2)!i,j,layer(1~up,2~down)
    complex(8),allocatable::XC(:,:,:,:)!i,j,position(1~left,2~right),number_of_up_layer_for_border
    complex(8),allocatable::XF(:)
    
    
    complex(8),save::last_alpha=(1.182735375d98,1.91724672d99),last_beta=(1.183245574d100,1.8245262591d101)
    
    private::already_counted,t,S,C,F,last_alpha,last_beta,count_t,count_K_element,generate_matrics,solve_matrics_system,&
        generate_matrics_isotropic,solve_matrics_system_isotropic
    private::generate_summary_matrix_A,generate_summary_matrix_AY,generate_summary_matrix_AX
    
    !todo delete
    complex(8),save,allocatable::A(:,:),B(:,:),AY(:,:),BY(:,:),AX(:,:),BX(:,:)
    contains
    
    complex(8) function K(i,j,alpha,beta,z) result(f)
    use system,only:print_error
    use main_parameters,only:get_full_h,get_anisotropic
    implicit none
        integer(4),intent(in)::i
        integer(4),intent(in)::j
        complex(8),intent(in)::alpha
        complex(8),intent(in)::beta
        real(8),intent(in)::z
        
        if(i<1.or.i>3) call print_error("count_K.K","(i<1.or.i>3")
        if(j<1.or.j>3) call print_error("count_K.K","(j<1.or.j>3")
        if(z>0d0.or.z<-get_full_h()) call print_error("count_K.K","z>0d0.or.z<-get_full_h()")
        
        !if(get_anisotropic()==0.and..not.beta==0d0) call print_error("count_K.K","anisotropic==0.and..not.beta==0d0")
        
        if(alpha/=last_alpha.or.beta/=last_beta) then
            last_alpha=alpha
            last_beta=beta
            already_counted(1)=.false.
            already_counted(2)=.false.
            already_counted(3)=.false.
        endif
        
        if(.not.already_counted(j)) then
            call count_t(alpha,beta,j)
            already_counted(j)=.true.
        endif
        
        f=count_K_element(i,j,z)
    endfunction K
    
    complex(8) function count_K_element(i,j,z) result(f)
    use sigma_and_eigenvectors,only:get_sigma,get_eigenvector_component
    use main_parameters,only:number_of_layers,h,get_full_h,Calphabeta,rho,omega,Q,get_number_of_layer,Cijkl,get_anisotropic
    use system,only:print_error
    use math,only:ci
    implicit none
        integer(4),intent(in)::i
        integer(4),intent(in)::j
        real(8),intent(in)::z
        
        real(8) z_
        integer(4) number_of_layer,ind
        
        complex(8) phi,psi,w,sigma1,sigma2
        
        if(i<1.or.i>3) call print_error("count_K.count_K_element","(i<1.or.i>3")
        if(j<1.or.j>3) call print_error("count_K.count_K_element","(j<1.or.j>3")
        if(z>0d0.or.z<-get_full_h()) call print_error("count_K.count_K_element","z>0d0.or.z<-get_full_h()")
        
        number_of_layer=get_number_of_layer(z)
        z_=z
        do ind=1,number_of_layer-1
            z_=z_+h(ind)
        enddo
        
        f=0d0
        
        if(get_anisotropic()==0) then
            phi=0d0
            psi=0d0
            w=0d0
            
            sigma1=get_sigma(last_alpha,last_beta,1,number_of_layer)
            sigma2=get_sigma(last_alpha,last_beta,3,number_of_layer)
            
            ! get_eigenvector_component(alpha,beta,number_of_sigma,index,number_of_layer)
            !Yt(:,:,:)!number_of_t_for_layer,number_of_layers,number_of_en
            if(.not.i==3) then
                phi=phi+Yt(1,number_of_layer,j)*get_eigenvector_component(last_alpha,last_beta,1,1,number_of_layer)*exp(sigma1*z_)
                phi=phi+Yt(2,number_of_layer,j)*get_eigenvector_component(last_alpha,last_beta,2,1,number_of_layer)*exp(-sigma1*(h(number_of_layer)+z_))
                phi=phi+Yt(3,number_of_layer,j)*get_eigenvector_component(last_alpha,last_beta,3,1,number_of_layer)*exp(sigma2*z_)
                phi=phi+Yt(4,number_of_layer,j)*get_eigenvector_component(last_alpha,last_beta,4,1,number_of_layer)*exp(-sigma2*(h(number_of_layer)+z_))
                
                psi=psi+Xt(1,number_of_layer,j)*get_eigenvector_component(last_alpha,last_beta,5,2,number_of_layer)*exp(sigma2*z_)
                psi=psi+Xt(2,number_of_layer,j)*get_eigenvector_component(last_alpha,last_beta,6,2,number_of_layer)*exp(-sigma2*(h(number_of_layer)+z_))
            endif
            
            if(i==3) then
                w=w+Yt(1,number_of_layer,j)*get_eigenvector_component(last_alpha,last_beta,1,3,number_of_layer)*exp(sigma1*z_)
                w=w+Yt(2,number_of_layer,j)*get_eigenvector_component(last_alpha,last_beta,2,3,number_of_layer)*exp(-sigma1*(h(number_of_layer)+z_))
                w=w+Yt(3,number_of_layer,j)*get_eigenvector_component(last_alpha,last_beta,3,3,number_of_layer)*exp(sigma2*z_)
                w=w+Yt(4,number_of_layer,j)*get_eigenvector_component(last_alpha,last_beta,4,3,number_of_layer)*exp(-sigma2*(h(number_of_layer)+z_))
            endif
            
            
            if(i==1) then
                f=-ci*last_alpha*phi-ci*last_beta*psi
            elseif(i==2) then
                f=-ci*last_beta*phi+ci*last_alpha*psi
            else
                f=w
            endif
            
            return
        endif
        
        do ind=1,3
            f=f+&
                t(ind,number_of_layer,j)*&
                get_eigenvector_component(last_alpha,last_beta,ind,i,number_of_layer)*&
                exp(get_sigma(last_alpha,last_beta,ind,number_of_layer)*z_)
        enddo
        do ind=4,6
            f=f+&
                t(ind,number_of_layer,j)*&
                get_eigenvector_component(last_alpha,last_beta,ind,i,number_of_layer)*&
                exp(get_sigma(last_alpha,last_beta,ind,number_of_layer)*(h(number_of_layer)+z_))
        enddo
    endfunction count_K_element
    
    subroutine count_t(alpha,beta,number_of_en)
    use system,only:print_error
    implicit none
        complex(8),intent(in)::alpha
        complex(8),intent(in)::beta
        integer(4),intent(in)::number_of_en
        
        if(number_of_en<1.or.number_of_en>3) call print_error("count_K.count_t","number_of_en<1.or.number_of_en>3")
        
        call generate_matrics(alpha,beta,number_of_en)
        call solve_matrics_system(number_of_en)
    endsubroutine count_t
    
    subroutine generate_matrics(alpha,beta,number_of_en)
    use main_parameters,only:number_of_layers,Cijkl,h,get_anisotropic,&
        get_down_border_condition_type,&
        get_down_border_condition_type_fixed_border,get_down_border_condition_type_free_border,&
        get_down_border_condition_type_halfspace
    use sigma_and_eigenvectors,only:get_sigma,get_eigenvector_component
    use system,only:print_error
    use math,only:ci
    implicit none
        complex(8),intent(in)::alpha
        complex(8),intent(in)::beta
        integer(4),intent(in)::number_of_en
        
        integer(4) i,j,k,l,layer
        
        if(number_of_en<1.or.number_of_en>3) call print_error("count_K.generate_matrics","number_of_en<1.or.number_of_en>3")
        
        if(get_anisotropic()==0) then
            call generate_matrics_isotropic(alpha,beta,number_of_en)
            
            return
        endif
        
        do i=1,number_of_layers*6
            F(i)=0d0
        enddo
        F(number_of_en)=1d0
        
        if(already_counted(1).or.already_counted(2).or.already_counted(3)) return
        
        do i=1,3
            do j=1,3
                S(i,j,1)=(&
                    
                    Cijkl(i,3,1,1,1)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,1,1)+&
                    Cijkl(i,3,1,2,1)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,1,1)+&
                    Cijkl(i,3,1,3,1)*(get_sigma(alpha,beta,j,1))*get_eigenvector_component(alpha,beta,j,1,1)+&
                    
                    Cijkl(i,3,2,1,1)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,2,1)+&
                    Cijkl(i,3,2,2,1)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,2,1)+&
                    Cijkl(i,3,2,3,1)*(get_sigma(alpha,beta,j,1))*get_eigenvector_component(alpha,beta,j,2,1)+&
                    
                    Cijkl(i,3,3,1,1)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,3,1)+&
                    Cijkl(i,3,3,2,1)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,3,1)+&
                    Cijkl(i,3,3,3,1)*(get_sigma(alpha,beta,j,1))*get_eigenvector_component(alpha,beta,j,3,1)&
                    )!*exp(get_sigma(alpha,beta,j,1)*0)
            enddo
            do j=4,6
                S(i,j,1)=(&
                    
                    Cijkl(i,3,1,1,1)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,1,1)+&
                    Cijkl(i,3,1,2,1)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,1,1)+&
                    Cijkl(i,3,1,3,1)*(get_sigma(alpha,beta,j,1))*get_eigenvector_component(alpha,beta,j,1,1)+&
                    
                    Cijkl(i,3,2,1,1)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,2,1)+&
                    Cijkl(i,3,2,2,1)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,2,1)+&
                    Cijkl(i,3,2,3,1)*(get_sigma(alpha,beta,j,1))*get_eigenvector_component(alpha,beta,j,2,1)+&
                    
                    Cijkl(i,3,3,1,1)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,3,1)+&
                    Cijkl(i,3,3,2,1)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,3,1)+&
                    Cijkl(i,3,3,3,1)*(get_sigma(alpha,beta,j,1))*get_eigenvector_component(alpha,beta,j,3,1)&
                    )*exp(get_sigma(alpha,beta,j,1)*h(1))
            enddo
        enddo
        
        do layer=1,number_of_layers-1
            do i=1,3
                do j=1,3
                    C(i,j,1,layer)=get_eigenvector_component(alpha,beta,j,i,layer)*exp(-get_sigma(alpha,beta,j,layer)*h(layer))
                    C(i,j,2,layer)=-get_eigenvector_component(alpha,beta,j,i,layer+1)!*exp(get_sigma(alpha,beta,j,layer+1)*0d0)
                enddo
                do j=4,6
                    C(i,j,1,layer)=get_eigenvector_component(alpha,beta,j,i,layer)!*exp(get_sigma(alpha,beta,j,layer)*0d0)
                    C(i,j,2,layer)=-get_eigenvector_component(alpha,beta,j,i,layer+1)*exp(get_sigma(alpha,beta,j,layer+1)*h(layer+1))
                enddo
            enddo
            do i=1,3
                do j=1,3
                    C(i+3,j,1,layer)=(&
                        
                        Cijkl(i,3,1,1,layer)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,1,layer)+&
                        Cijkl(i,3,1,2,layer)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,1,layer)+&
                        Cijkl(i,3,1,3,layer)*(get_sigma(alpha,beta,j,layer))*get_eigenvector_component(alpha,beta,j,1,layer)+&
                        
                        Cijkl(i,3,2,1,layer)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,2,layer)+&
                        Cijkl(i,3,2,2,layer)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,2,layer)+&
                        Cijkl(i,3,2,3,layer)*(get_sigma(alpha,beta,j,layer))*get_eigenvector_component(alpha,beta,j,2,layer)+&
                        
                        Cijkl(i,3,3,1,layer)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,3,layer)+&
                        Cijkl(i,3,3,2,layer)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,3,layer)+&
                        Cijkl(i,3,3,3,layer)*(get_sigma(alpha,beta,j,layer))*get_eigenvector_component(alpha,beta,j,3,layer)&
                        )*exp(-get_sigma(alpha,beta,j,layer)*h(layer))
                    C(i+3,j,2,layer)=-(&
                        
                        Cijkl(i,3,1,1,layer+1)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,1,layer+1)+&
                        Cijkl(i,3,1,2,layer+1)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,1,layer+1)+&
                        Cijkl(i,3,1,3,layer+1)*(get_sigma(alpha,beta,j,layer+1))*get_eigenvector_component(alpha,beta,j,1,layer+1)+&
                        
                        Cijkl(i,3,2,1,layer+1)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,2,layer+1)+&
                        Cijkl(i,3,2,2,layer+1)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,2,layer+1)+&
                        Cijkl(i,3,2,3,layer+1)*(get_sigma(alpha,beta,j,layer+1))*get_eigenvector_component(alpha,beta,j,2,layer+1)+&
                        
                        Cijkl(i,3,3,1,layer+1)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,3,layer+1)+&
                        Cijkl(i,3,3,2,layer+1)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,3,layer+1)+&
                        Cijkl(i,3,3,3,layer+1)*(get_sigma(alpha,beta,j,layer+1))*get_eigenvector_component(alpha,beta,j,3,layer+1)&
                        )!*exp(get_sigma(alpha,beta,j,layer+1)*0d0)
                enddo
                do j=4,6
                    C(i+3,j,1,layer)=(&
                        
                        Cijkl(i,3,1,1,layer)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,1,layer)+&
                        Cijkl(i,3,1,2,layer)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,1,layer)+&
                        Cijkl(i,3,1,3,layer)*(get_sigma(alpha,beta,j,layer))*get_eigenvector_component(alpha,beta,j,1,layer)+&
                        
                        Cijkl(i,3,2,1,layer)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,2,layer)+&
                        Cijkl(i,3,2,2,layer)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,2,layer)+&
                        Cijkl(i,3,2,3,layer)*(get_sigma(alpha,beta,j,layer))*get_eigenvector_component(alpha,beta,j,2,layer)+&
                        
                        Cijkl(i,3,3,1,layer)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,3,layer)+&
                        Cijkl(i,3,3,2,layer)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,3,layer)+&
                        Cijkl(i,3,3,3,layer)*(get_sigma(alpha,beta,j,layer))*get_eigenvector_component(alpha,beta,j,3,layer)&
                        )!*exp(get_sigma(alpha,beta,j,layer)*0d0)
                    C(i+3,j,2,layer)=-(&
                        
                        Cijkl(i,3,1,1,layer+1)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,1,layer+1)+&
                        Cijkl(i,3,1,2,layer+1)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,1,layer+1)+&
                        Cijkl(i,3,1,3,layer+1)*(get_sigma(alpha,beta,j,layer+1))*get_eigenvector_component(alpha,beta,j,1,layer+1)+&
                        
                        Cijkl(i,3,2,1,layer+1)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,2,layer+1)+&
                        Cijkl(i,3,2,2,layer+1)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,2,layer+1)+&
                        Cijkl(i,3,2,3,layer+1)*(get_sigma(alpha,beta,j,layer+1))*get_eigenvector_component(alpha,beta,j,2,layer+1)+&
                        
                        Cijkl(i,3,3,1,layer+1)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,3,layer+1)+&
                        Cijkl(i,3,3,2,layer+1)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,3,layer+1)+&
                        Cijkl(i,3,3,3,layer+1)*(get_sigma(alpha,beta,j,layer+1))*get_eigenvector_component(alpha,beta,j,3,layer+1)&
                        )*exp(get_sigma(alpha,beta,j,layer+1)*h(layer+1))
                enddo
            enddo
        enddo
        
        if(get_down_border_condition_type()==get_down_border_condition_type_fixed_border()) then
            do i=1,3
                do j=1,3
                    S(i,j,2)=(&
                        get_eigenvector_component(alpha,beta,j,i,number_of_layers)&
                        )*exp(-get_sigma(alpha,beta,j,number_of_layers)*h(number_of_layers))
                enddo
                do j=4,6
                    S(i,j,2)=(&
                        get_eigenvector_component(alpha,beta,j,i,number_of_layers)&
                        )!*exp(get_sigma(alpha,beta,j,number_of_layers)*0)
                enddo
            enddo
        elseif(get_down_border_condition_type()==get_down_border_condition_type_free_border()) then
            do i=1,3
                do j=1,3
                    S(i,j,2)=(&
                        
                        Cijkl(i,3,1,1,number_of_layers)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,1,number_of_layers)+&
                        Cijkl(i,3,1,2,number_of_layers)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,1,number_of_layers)+&
                        Cijkl(i,3,1,3,number_of_layers)*(get_sigma(alpha,beta,j,number_of_layers))*get_eigenvector_component(alpha,beta,j,1,number_of_layers)+&
                        
                        Cijkl(i,3,2,1,number_of_layers)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,2,number_of_layers)+&
                        Cijkl(i,3,2,2,number_of_layers)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,2,number_of_layers)+&
                        Cijkl(i,3,2,3,number_of_layers)*(get_sigma(alpha,beta,j,number_of_layers))*get_eigenvector_component(alpha,beta,j,2,number_of_layers)+&
                        
                        Cijkl(i,3,3,1,number_of_layers)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,3,number_of_layers)+&
                        Cijkl(i,3,3,2,number_of_layers)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,3,number_of_layers)+&
                        Cijkl(i,3,3,3,number_of_layers)*(get_sigma(alpha,beta,j,number_of_layers))*get_eigenvector_component(alpha,beta,j,3,number_of_layers)&
                        )*exp(-get_sigma(alpha,beta,j,number_of_layers)*h(number_of_layers))
                enddo
                do j=4,6
                    S(i,j,2)=(&
                        
                        Cijkl(i,3,1,1,number_of_layers)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,1,number_of_layers)+&
                        Cijkl(i,3,1,2,number_of_layers)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,1,number_of_layers)+&
                        Cijkl(i,3,1,3,number_of_layers)*(get_sigma(alpha,beta,j,number_of_layers))*get_eigenvector_component(alpha,beta,j,1,number_of_layers)+&
                        
                        Cijkl(i,3,2,1,number_of_layers)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,2,number_of_layers)+&
                        Cijkl(i,3,2,2,number_of_layers)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,2,number_of_layers)+&
                        Cijkl(i,3,2,3,number_of_layers)*(get_sigma(alpha,beta,j,number_of_layers))*get_eigenvector_component(alpha,beta,j,2,number_of_layers)+&
                        
                        Cijkl(i,3,3,1,number_of_layers)*(-ci*alpha)*get_eigenvector_component(alpha,beta,j,3,number_of_layers)+&
                        Cijkl(i,3,3,2,number_of_layers)*(-ci*beta)*get_eigenvector_component(alpha,beta,j,3,number_of_layers)+&
                        Cijkl(i,3,3,3,number_of_layers)*(get_sigma(alpha,beta,j,number_of_layers))*get_eigenvector_component(alpha,beta,j,3,number_of_layers)&
                        )!*exp(get_sigma(alpha,beta,j,number_of_layers)*0)
                enddo
            enddo
        elseif(get_down_border_condition_type()==get_down_border_condition_type_halfspace()) then
            S(1,4,2)=1d0
            S(2,5,2)=1d0
            S(3,6,2)=1d0
        else
            call print_error("count_K.generate_matrics","a current get_down_border_condition_type for anisotropic is not supported yet")
        endif
    endsubroutine generate_matrics
    
    subroutine solve_matrics_system(number_of_en)
    use main_parameters,only:number_of_layers,get_anisotropic
    use linear_system,only:star5_
    use system,only:print_error
    implicit none
        integer(4),intent(in)::number_of_en
        
        integer(4) i,j,layer
        
        if(number_of_en<1.or.number_of_en>3) call print_error("count_K.solve_matrics_system","number_of_en<1.or.number_of_en>3")
        
        if(get_anisotropic()==0) then
            call solve_matrics_system_isotropic(number_of_en)
            
            return
        endif
        
        do i=1,number_of_layers*6
            B(i,1)=0d0
        enddo
        B(number_of_en,1)=1d0
        
        call generate_summary_matrix_A()
        
        call star5_(A,B)
        
        j=1
        do layer=1,number_of_layers
            do i=1,6
                t(i,layer,number_of_en)=B(j,1)
                j=j+1
            enddo
        enddo
        
        !integer(4) i,j,layer
        !real(8) norm_current,norm_last,min_diff
        !
        !do i=1,number_of_layers
        !    do j=1,6
        !        t(j,i,number_of_en)=0d0
        !    enddo
        !enddo
        !
        !
        !min_diff=1d-4;
        !norm_last=1d100
        !norm_current=1d90
        !
        !do while(abs(norm_current-norm_last)>=min_diff*norm_current)
        !    t(1,1,number_of_en)=F(1)&
        !        -t(2,1,number_of_en)*S(1,2,1)&
        !        -t(3,1,number_of_en)*S(1,3,1)&
        !        -t(4,1,number_of_en)*S(1,4,1)&
        !        -t(5,1,number_of_en)*S(1,5,1)&
        !        -t(6,1,number_of_en)*S(1,6,1)
        !    t(1,1,number_of_en)=t(1,1,number_of_en)/S(1,1,1)
        !    t(2,1,number_of_en)=F(2)&
        !        -t(1,1,number_of_en)*S(2,1,1)&
        !        -t(3,1,number_of_en)*S(2,3,1)&
        !        -t(4,1,number_of_en)*S(2,4,1)&
        !        -t(5,1,number_of_en)*S(2,5,1)&
        !        -t(6,1,number_of_en)*S(2,6,1)
        !    t(2,1,number_of_en)=t(2,1,number_of_en)/S(2,2,1)
        !    t(3,1,number_of_en)=F(3)&
        !        -t(1,1,number_of_en)*S(3,1,1)&
        !        -t(2,1,number_of_en)*S(3,2,1)&
        !        -t(4,1,number_of_en)*S(3,4,1)&
        !        -t(5,1,number_of_en)*S(3,5,1)&
        !        -t(6,1,number_of_en)*S(3,6,1)
        !    t(3,1,number_of_en)=t(3,1,number_of_en)/S(3,3,1)
        !    
        !    do layer=1,number_of_layers-1
        !        t(4,layer,number_of_en)=&
        !            -t(1,layer,number_of_en)*C(1,1,1,layer)&
        !            -t(2,layer,number_of_en)*C(1,2,1,layer)&
        !            -t(3,layer,number_of_en)*C(1,3,1,layer)&
        !            -t(5,layer,number_of_en)*C(1,5,1,layer)&
        !            -t(6,layer,number_of_en)*C(1,6,1,layer)&
        !            -t(1,layer+1,number_of_en)*C(1,1,2,layer)&
        !            -t(2,layer+1,number_of_en)*C(1,2,2,layer)&
        !            -t(3,layer+1,number_of_en)*C(1,3,2,layer)&
        !            -t(4,layer+1,number_of_en)*C(1,4,2,layer)&
        !            -t(5,layer+1,number_of_en)*C(1,5,2,layer)&
        !            -t(6,layer+1,number_of_en)*C(1,6,2,layer)
        !        t(4,layer,number_of_en)=t(4,layer,number_of_en)/C(1,4,1,layer)
        !        
        !        t(5,layer,number_of_en)=&
        !            -t(1,layer,number_of_en)*C(2,1,1,layer)&
        !            -t(2,layer,number_of_en)*C(2,2,1,layer)&
        !            -t(3,layer,number_of_en)*C(2,3,1,layer)&
        !            -t(4,layer,number_of_en)*C(2,4,1,layer)&
        !            -t(6,layer,number_of_en)*C(2,6,1,layer)&
        !            -t(1,layer+1,number_of_en)*C(2,1,2,layer)&
        !            -t(2,layer+1,number_of_en)*C(2,2,2,layer)&
        !            -t(3,layer+1,number_of_en)*C(2,3,2,layer)&
        !            -t(4,layer+1,number_of_en)*C(2,4,2,layer)&
        !            -t(5,layer+1,number_of_en)*C(2,5,2,layer)&
        !            -t(6,layer+1,number_of_en)*C(2,6,2,layer)
        !        t(5,layer,number_of_en)=t(5,layer,number_of_en)/C(2,5,1,layer)
        !            
        !        t(6,layer,number_of_en)=&
        !            -t(1,layer,number_of_en)*C(3,1,1,layer)&
        !            -t(2,layer,number_of_en)*C(3,2,1,layer)&
        !            -t(3,layer,number_of_en)*C(3,3,1,layer)&
        !            -t(4,layer,number_of_en)*C(3,4,1,layer)&
        !            -t(5,layer,number_of_en)*C(3,5,1,layer)&
        !            -t(1,layer+1,number_of_en)*C(3,1,2,layer)&
        !            -t(2,layer+1,number_of_en)*C(3,2,2,layer)&
        !            -t(3,layer+1,number_of_en)*C(3,3,2,layer)&
        !            -t(4,layer+1,number_of_en)*C(3,4,2,layer)&
        !            -t(5,layer+1,number_of_en)*C(3,5,2,layer)&
        !            -t(6,layer+1,number_of_en)*C(3,6,2,layer)
        !        t(6,layer,number_of_en)=t(6,layer,number_of_en)/C(3,6,1,layer)
        !        
        !        t(1,layer+1,number_of_en)=&
        !            -t(1,layer,number_of_en)*C(4,1,1,layer)&
        !            -t(2,layer,number_of_en)*C(4,2,1,layer)&
        !            -t(3,layer,number_of_en)*C(4,3,1,layer)&
        !            -t(4,layer,number_of_en)*C(4,4,1,layer)&
        !            -t(5,layer,number_of_en)*C(4,5,1,layer)&
        !            -t(6,layer,number_of_en)*C(4,6,1,layer)&
        !            -t(2,layer+1,number_of_en)*C(4,2,2,layer)&
        !            -t(3,layer+1,number_of_en)*C(4,3,2,layer)&
        !            -t(4,layer+1,number_of_en)*C(4,4,2,layer)&
        !            -t(5,layer+1,number_of_en)*C(4,5,2,layer)&
        !            -t(6,layer+1,number_of_en)*C(4,6,2,layer)
        !        t(1,layer+1,number_of_en)=t(1,layer+1,number_of_en)/C(4,1,2,layer)
        !        
        !        t(2,layer+1,number_of_en)=&
        !            -t(1,layer,number_of_en)*C(5,1,1,layer)&
        !            -t(2,layer,number_of_en)*C(5,2,1,layer)&
        !            -t(3,layer,number_of_en)*C(5,3,1,layer)&
        !            -t(4,layer,number_of_en)*C(5,4,1,layer)&
        !            -t(5,layer,number_of_en)*C(5,5,1,layer)&
        !            -t(6,layer,number_of_en)*C(5,6,1,layer)&
        !            -t(1,layer+1,number_of_en)*C(5,1,2,layer)&
        !            -t(3,layer+1,number_of_en)*C(5,3,2,layer)&
        !            -t(4,layer+1,number_of_en)*C(5,4,2,layer)&
        !            -t(5,layer+1,number_of_en)*C(5,5,2,layer)&
        !            -t(6,layer+1,number_of_en)*C(5,6,2,layer)
        !        t(2,layer+1,number_of_en)=t(2,layer+1,number_of_en)/C(5,2,2,layer)
        !        
        !        t(3,layer+1,number_of_en)=&
        !            -t(1,layer,number_of_en)*C(6,1,1,layer)&
        !            -t(2,layer,number_of_en)*C(6,2,1,layer)&
        !            -t(3,layer,number_of_en)*C(6,3,1,layer)&
        !            -t(4,layer,number_of_en)*C(6,4,1,layer)&
        !            -t(5,layer,number_of_en)*C(6,5,1,layer)&
        !            -t(6,layer,number_of_en)*C(6,6,1,layer)&
        !            -t(1,layer+1,number_of_en)*C(6,1,2,layer)&
        !            -t(2,layer+1,number_of_en)*C(6,2,2,layer)&
        !            -t(4,layer+1,number_of_en)*C(6,4,2,layer)&
        !            -t(5,layer+1,number_of_en)*C(6,5,2,layer)&
        !            -t(6,layer+1,number_of_en)*C(6,6,2,layer)
        !        t(3,layer+1,number_of_en)=t(3,layer+1,number_of_en)/C(6,3,2,layer)
        !    enddo
        !    
        !    t(4,number_of_layers,number_of_en)=&
        !        -t(1,number_of_layers,number_of_en)*S(1,1,2)&
        !        -t(2,number_of_layers,number_of_en)*S(1,2,2)&
        !        -t(3,number_of_layers,number_of_en)*S(1,3,2)&
        !        -t(5,number_of_layers,number_of_en)*S(1,5,2)&
        !        -t(6,number_of_layers,number_of_en)*S(1,6,2)
        !    t(4,number_of_layers,number_of_en)=t(4,number_of_layers,number_of_en)/S(1,4,2)
        !    t(5,number_of_layers,number_of_en)=&
        !        -t(1,number_of_layers,number_of_en)*S(2,1,2)&
        !        -t(2,number_of_layers,number_of_en)*S(2,2,2)&
        !        -t(3,number_of_layers,number_of_en)*S(2,3,2)&
        !        -t(4,number_of_layers,number_of_en)*S(2,4,2)&
        !        -t(6,number_of_layers,number_of_en)*S(2,6,2)
        !    t(5,number_of_layers,number_of_en)=t(5,number_of_layers,number_of_en)/S(2,5,2)
        !    t(6,number_of_layers,number_of_en)=&
        !        -t(1,number_of_layers,number_of_en)*S(3,1,2)&
        !        -t(2,number_of_layers,number_of_en)*S(3,2,2)&
        !        -t(3,number_of_layers,number_of_en)*S(3,3,2)&
        !        -t(4,number_of_layers,number_of_en)*S(3,4,2)&
        !        -t(5,number_of_layers,number_of_en)*S(3,5,2)
        !    t(6,number_of_layers,number_of_en)=t(6,number_of_layers,number_of_en)/S(3,6,2)
        !    
        !    norm_last=norm_current
        !    norm_current=0d0
        !    do i=1,number_of_layers
        !        do j=1,6
        !            if(norm_current<abs(t(j,i,number_of_en))) then
        !                norm_current=abs(t(j,i,number_of_en))
        !            endif
        !        enddo
        !    enddo
        !    !print*,"norm",norm_current
        !    !pause " "
        !    
        !    !do i=1,3
        !    !    do j=1,6
        !    !        print*,S(i,j,1)
        !    !    enddo
        !    !enddo
        !    !print*
        !    !do layer=1,number_of_layers-1
        !    !    do i=1,6
        !    !        do j=1,6
        !    !            print*,C(i,j,1,layer)
        !    !        enddo
        !    !    enddo
        !    !    print*
        !    !    do i=1,6
        !    !        do j=1,6
        !    !            print*,C(i,j,2,layer)
        !    !        enddo
        !    !    enddo
        !    !    print*
        !    !enddo
        !    !do i=1,3
        !    !    do j=1,6
        !    !        print*,S(i,j,2)
        !    !    enddo
        !    !enddo
        !    !print*
        !    !
        !    !do layer=1,number_of_layers
        !    !    do i=1,6
        !    !        print*,t(i,layer,number_of_en)
        !    !    enddo
        !    !    print*
        !    !enddo
        !    !
        !    !pause " "
        !enddo
    endsubroutine solve_matrics_system
    subroutine generate_summary_matrix_A()
    use main_parameters,only:number_of_layers
    implicit none
        integer(4) i,j,layer
        
        do i=1,number_of_layers*6
            do j=1,number_of_layers*6
                    A(i,j)=0d0
            enddo
        enddo
        
        do i=1,3
            do j=1,6
                A(i,j)=S(i,j,1)
            enddo
        enddo
        do layer=1,number_of_layers-1
            do i=1,6
                do j=1,6
                    A(3+(layer-1)*6+i,6*(layer-1)+j)=C(i,j,1,layer)
                    A(3+(layer-1)*6+i,6*layer+j)=C(i,j,2,layer)
                enddo
            enddo
        enddo
        do i=1,3
            do j=1,6
                A(number_of_layers*6-3+i,number_of_layers*6-6+j)=S(i,j,2)
            enddo
        enddo
    endsubroutine generate_summary_matrix_A
    
    subroutine generate_matrics_isotropic(alpha,beta,number_of_en)
    use main_parameters,only:number_of_layers,lambda,mu,h,get_anisotropic,&
        get_down_border_condition_type,&
        get_down_border_condition_type_fixed_border,get_down_border_condition_type_free_border,&
        get_down_border_condition_type_halfspace
    use sigma_and_eigenvectors,only:get_sigma,get_eigenvector_component
    use system,only:print_error
    use math,only:ci
    implicit none
        complex(8),intent(in)::alpha
        complex(8),intent(in)::beta
        integer(4),intent(in)::number_of_en
        
        integer(4) i,j,k,l,layer
        
        complex(8) sigma1,sigma2,alpha2,div_mul(4)
        
        if(number_of_en<1.or.number_of_en>3) call print_error("count_K.generate_matrics","number_of_en<1.or.number_of_en>3")
        
        alpha2=alpha*alpha+beta*beta
        
        !complex(8),allocatable::Yt(:,:,:)!number_of_t_for_layer,number_of_layers,number_of_en
        !complex(8),save::YS(2,4,2)!i,j,layer(1~up,2~down)
        !complex(8),allocatable::YC(:,:,:,:)!i,j,position(1~left,2~right),number_of_up_layer_for_border
        !complex(8),allocatable::YF(:)
        
        do i=1,4*number_of_layers
            YF(i)=0d0
        enddo
        if(number_of_en==1) then
            YF(2)=alpha
        elseif(number_of_en==2) then
            YF(2)=beta
        else
            YF(1)=1d0
        endif
        
        ! get_eigenvector_component(alpha,beta,number_of_sigma,index,number_of_layer)
        sigma1=get_sigma(alpha,beta,1,number_of_layer=1)
        sigma2=get_sigma(alpha,beta,3,number_of_layer=1)
        div_mul(1)=sigma1
        div_mul(2)=-sigma1
        div_mul(3)=sigma2
        div_mul(4)=-sigma2
        do j=1,4
            YS(1,j,1)=&
                -lambda(1)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=1)&
                +(lambda(1)+mu(1)+mu(1))*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=3,number_of_layer=1)*div_mul(j)
            YS(2,j,1)=&
                -ci*mu(1)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=1)*div_mul(j)&
                -ci*mu(1)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=3,number_of_layer=1)
            
            !YS(1,j,1)=&
            !    -lambda(1)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=1)&
            !    +(lambda(1)+mu(1)+mu(1))*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=6,number_of_layer=1)
            !YS(2,j,1)=&
            !    -ci*mu(1)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=2,number_of_layer=1)&
            !    -ci*mu(1)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=5,number_of_layer=1)
        enddo
        !YS(1,1,1)=YS(1,1,1)*exp(sigma1*0d0)
        YS(1,2,1)=YS(1,2,1)*exp(-sigma1*h(1))
        !YS(1,3,1)=YS(1,3,1)*exp(sigma2*0d0)
        YS(1,4,1)=YS(1,4,1)*exp(-sigma2*h(1))
        !YS(2,1,1)=YS(2,1,1)*exp(sigma1*0d0)
        YS(2,2,1)=YS(2,2,1)*exp(-sigma1*h(1))
        !YS(2,3,1)=YS(2,3,1)*exp(sigma2*0d0)
        YS(2,4,1)=YS(2,4,1)*exp(-sigma2*h(1))
        
        do layer=1,number_of_layers-1
            do j=1,4
                sigma1=get_sigma(alpha,beta,1,number_of_layer=layer)
                sigma2=get_sigma(alpha,beta,3,number_of_layer=layer)
                div_mul(1)=sigma1
                div_mul(2)=-sigma1
                div_mul(3)=sigma2
                div_mul(4)=-sigma2
                YC(1,j,1,layer)=&
                    -lambda(layer)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=layer)&
                    +(lambda(layer)+mu(layer)+mu(layer))*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=3,number_of_layer=layer)*div_mul(j)
                YC(2,j,1,layer)=&
                    -ci*mu(layer)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=layer)*div_mul(j)&
                    -ci*mu(layer)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=3,number_of_layer=layer)
                
                sigma1=get_sigma(alpha,beta,1,number_of_layer=layer+1)
                sigma2=get_sigma(alpha,beta,3,number_of_layer=layer+1)
                div_mul(1)=sigma1
                div_mul(2)=-sigma1
                div_mul(3)=sigma2
                div_mul(4)=-sigma2
                YC(1,j,2,layer)=&
                    +lambda(layer+1)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=layer+1)&
                    -(lambda(layer+1)+mu(layer+1)+mu(layer+1))*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=3,number_of_layer=layer+1)*div_mul(j)
                YC(2,j,2,layer)=&
                    +ci*mu(layer+1)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=layer+1)*div_mul(j)&
                    +ci*mu(layer+1)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=3,number_of_layer=layer+1)
                
                !sigma1=get_sigma(alpha,beta,1,number_of_layer=layer)
                !sigma2=get_sigma(alpha,beta,3,number_of_layer=layer)
                !YC(1,j,1,layer)=&
                !    -lambda(layer)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=layer)&
                !    +(lambda(layer)+mu(layer)+mu(layer))*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=6,number_of_layer=layer)
                !YC(2,j,1,layer)=&
                !    -ci*mu(layer)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=2,number_of_layer=layer)&
                !    -ci*mu(layer)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=5,number_of_layer=layer)
                
                !sigma1=get_sigma(alpha,beta,1,number_of_layer=layer+1)
                !sigma2=get_sigma(alpha,beta,3,number_of_layer=layer+1)
                !YC(1,j,2,layer)=&
                !    +lambda(layer+1)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=layer+1)&
                !    -(lambda(layer+1)+mu(layer+1)+mu(layer+1))*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=6,number_of_layer=layer+1)
                !YC(2,j,2,layer)=&
                !    +ci*alpha2*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=2,number_of_layer=layer+1)&
                !    +ci*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=5,number_of_layer=layer+1)
            enddo
            
            sigma1=get_sigma(alpha,beta,1,number_of_layer=layer)
            sigma2=get_sigma(alpha,beta,3,number_of_layer=layer)
            YC(1,1,1,layer)=YC(1,1,1,layer)*exp(-sigma1*h(layer))
            !YC(1,2,1,layer)=YC(1,2,1,layer)*exp(-sigma1*0d0)
            YC(1,3,1,layer)=YC(1,3,1,layer)*exp(-sigma2*h(layer))
            !YC(1,4,1,layer)=YC(1,4,1,layer)*exp(-sigma2*0d0)
            YC(2,1,1,layer)=YC(2,1,1,layer)*exp(-sigma1*h(layer))
            !YC(2,2,1,layer)=YC(2,2,1,layer)*exp(-sigma1*0d0)
            YC(2,3,1,layer)=YC(2,3,1,layer)*exp(-sigma2*h(layer))
            !YC((2,4,1,layer)=YC(2,4,1,layer)*exp(-sigma2*h(layer))
            
            sigma1=get_sigma(alpha,beta,1,number_of_layer=layer+1)
            sigma2=get_sigma(alpha,beta,3,number_of_layer=layer+1)
            !YC(1,1,2,layer)=YC(1,1,2,layer)*exp(sigma1*0d0)
            YC(1,2,2,layer)=YC(1,2,2,layer)*exp(-sigma1*h(layer+1))
            !YC(1,3,2,layer)=YC(1,3,2,layer)*exp(sigma2*0d0)
            YC(1,4,2,layer)=YC(1,4,2,layer)*exp(-sigma2*h(layer+1))
            !YC(2,1,2,layer)=YC(2,1,2,layer)*exp(sigma1*0d0)
            YC(2,2,2,layer)=YC(2,2,2,layer)*exp(-sigma1*h(layer+1))
            !YC(2,3,2,layer)=YC(2,3,2,layer)*exp(sigma2*0d0)
            YC(2,4,2,layer)=YC(2,4,2,layer)*exp(-sigma2*h(layer+1))
            
            do j=1,4
                sigma1=get_sigma(alpha,beta,1,number_of_layer=layer)
                sigma2=get_sigma(alpha,beta,3,number_of_layer=layer)
                YC(3,j,1,layer)=&
                    +get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=layer)
                YC(4,j,1,layer)=&
                    +get_eigenvector_component(alpha,beta,number_of_sigma=j,index=3,number_of_layer=layer)
                
                sigma1=get_sigma(alpha,beta,1,number_of_layer=layer+1)
                sigma2=get_sigma(alpha,beta,3,number_of_layer=layer+1)
                YC(3,j,2,layer)=&
                    -get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=layer+1)
                YC(4,j,2,layer)=&
                    -get_eigenvector_component(alpha,beta,number_of_sigma=j,index=3,number_of_layer=layer+1)
                
                !sigma1=get_sigma(alpha,beta,1,number_of_layer=layer)
                !sigma2=get_sigma(alpha,beta,3,number_of_layer=layer)
                !YC(3,j,1,layer)=&
                !    +get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=layer)
                !YC(4,j,1,layer)=&
                !    +get_eigenvector_component(alpha,beta,number_of_sigma=j,index=5,number_of_layer=layer)
                !
                !sigma1=get_sigma(alpha,beta,1,number_of_layer=layer+1)
                !sigma2=get_sigma(alpha,beta,3,number_of_layer=layer+1)
                !YC(3,j,2,layer)=&
                !    -get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=layer+1)
                !YC(4,j,2,layer)=&
                !    -get_eigenvector_component(alpha,beta,number_of_sigma=j,index=5,number_of_layer=layer+1)
            enddo
            
            sigma1=get_sigma(alpha,beta,1,number_of_layer=layer)
            sigma2=get_sigma(alpha,beta,3,number_of_layer=layer)
            YC(3,1,1,layer)=YC(3,1,1,layer)*exp(-sigma1*h(layer))
            !YC(3,2,1,layer)=YC(3,2,1,layer)*exp(-sigma1*0d0)
            YC(3,3,1,layer)=YC(3,3,1,layer)*exp(-sigma2*h(layer))
            !YC(3,4,1,layer)=YC(3,4,1,layer)*exp(-sigma2*0d0)
            YC(4,1,1,layer)=YC(4,1,1,layer)*exp(-sigma1*h(layer))
            !YC(4,2,1,layer)=YC(4,2,1,layer)*exp(-sigma1*0d0)
            YC(4,3,1,layer)=YC(4,3,1,layer)*exp(-sigma2*h(layer))
            !YC(4,4,1,layer)=YC(4,4,1,layer)*exp(-sigma2*h(layer))
            
            sigma1=get_sigma(alpha,beta,1,number_of_layer=layer+1)
            sigma2=get_sigma(alpha,beta,3,number_of_layer=layer+1)
            !YC(3,1,2,layer)=YC(3,1,2,layer)*exp(sigma1*0d0)
            YC(3,2,2,layer)=YC(3,2,2,layer)*exp(-sigma1*h(layer+1))
            !YC(3,3,2,layer)=YC(3,3,2,layer)*exp(sigma2*0d0)
            YC(3,4,2,layer)=YC(3,4,2,layer)*exp(-sigma2*h(layer+1))
            !YC(4,1,2,layer)=YC(4,1,2,layer)*exp(sigma1*0d0)
            YC(4,2,2,layer)=YC(4,2,2,layer)*exp(-sigma1*h(layer+1))
            !YC(4,3,2,layer)=YC(4,3,2,layer)*exp(sigma2*0d0)
            YC(4,4,2,layer)=YC(4,4,2,layer)*exp(-sigma2*h(layer+1))
        enddo
        
        sigma1=get_sigma(alpha,beta,1,number_of_layer=number_of_layers)
        sigma2=get_sigma(alpha,beta,3,number_of_layer=number_of_layers)
        if(get_down_border_condition_type()==get_down_border_condition_type_fixed_border()) then
            do j=1,4
                YS(1,j,2)=&
                    +get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=number_of_layers)
                YS(2,j,2)=&
                    +get_eigenvector_component(alpha,beta,number_of_sigma=j,index=3,number_of_layer=number_of_layers)
                
                !YS(1,j,2)=&
                !    +get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=number_of_layers)
                !YS(2,j,2)=&
                !    +get_eigenvector_component(alpha,beta,number_of_sigma=j,index=5,number_of_layer=number_of_layers)
            enddo
        elseif(get_down_border_condition_type()==get_down_border_condition_type_free_border()) then
            div_mul(1)=sigma1
            div_mul(2)=-sigma1
            div_mul(3)=sigma2
            div_mul(4)=-sigma2
            do j=1,4
                YS(1,j,2)=&
                    -lambda(number_of_layers)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=number_of_layers)&
                    +(lambda(number_of_layers)+mu(number_of_layers)+mu(number_of_layers))*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=3,number_of_layer=number_of_layers)*div_mul(j)
                YS(2,j,2)=&
                    -ci*mu(number_of_layers)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=number_of_layers)*div_mul(j)&
                    -ci*mu(number_of_layers)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=3,number_of_layer=number_of_layers)
                
                !YS(1,j,2)=&
                !    -lambda(number_of_layers)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=1,number_of_layer=number_of_layers)&
                !    +(lambda(number_of_layers)+mu(number_of_layers)+mu(number_of_layers))*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=6,number_of_layer=number_of_layers)
                !YS(2,j,2)=&
                !    -ci*mu(number_of_layers)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=2,number_of_layer=number_of_layers)&
                !    -ci*mu(number_of_layers)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=j,index=5,number_of_layer=number_of_layers)
            enddo
        elseif(get_down_border_condition_type()==get_down_border_condition_type_halfspace()) then
            YS(1,2,2)=1d0
            YS(2,4,2)=1d0
        else
            call print_error("count_K.generate_matrics_isotropic","a current get_down_border_condition_type for isotropic is not supported yet")
        endif
        YS(1,1,2)=YS(1,1,2)*exp(-sigma1*h(number_of_layers))
        !YS(1,2,2)=YS(1,2,2)*exp(-sigma1*0d0)
        YS(1,3,2)=YS(1,3,2)*exp(-sigma2*h(number_of_layers))
        !YS(1,4,2)=YS(1,4,2)*exp(-sigma2*0d0)
        YS(2,1,2)=YS(2,1,2)*exp(-sigma1*h(number_of_layers))
        !YS(2,2,2)=YS(2,2,2)*exp(-sigma1*0d0)
        YS(2,3,2)=YS(2,3,2)*exp(-sigma2*h(number_of_layers))
        !YS(2,4,2)=YS(2,4,2)*exp(-sigma2*0d0)
        
        !complex(8),allocatable::Xt(:,:,:)!number_of_t_for_layer,number_of_layers,number_of_en
        !complex(8),save::XS(1,2,2)!i,j,layer(1~up,2~down)
        !complex(8),allocatable::XC(:,:,:,:)!i,j,position(1~left,2~right),number_of_up_layer_for_border
        !complex(8),allocatable::XF(:)
        
        do i=1,2*number_of_layers
            XF(i)=0d0
        enddo
        if(number_of_en==1) then
            XF(1)=-beta
        elseif(number_of_en==2) then
            XF(1)=alpha
        else
            return
        endif
        
        ! get_eigenvector_component(alpha,beta,number_of_sigma,index,number_of_layer)
        sigma2=get_sigma(alpha,beta,3,number_of_layer=1)
        div_mul(1)=sigma2
        div_mul(2)=-sigma2
        do j=1,2
            XS(1,j,1)=&
                ci*mu(1)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=4+j,index=2,number_of_layer=1)*div_mul(j)
            !XS(1,j,1)=&
            !    ci*mu(1)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=4+j,index=4,number_of_layer=1)
        enddo
        !XS(1,1,1)=XS(1,1,1)*exp(sigma2*0d0)
        XS(1,2,1)=XS(1,2,1)*exp(-sigma2*h(1))
        
        do layer=1,number_of_layers-1
            do j=1,2
                sigma2=get_sigma(alpha,beta,3,number_of_layer=layer)
                div_mul(1)=sigma2
                div_mul(2)=-sigma2
                XC(1,j,1,layer)=&
                    ci*mu(layer)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=4+j,index=2,number_of_layer=layer)*div_mul(j)
                
                sigma2=get_sigma(alpha,beta,3,number_of_layer=layer+1)
                div_mul(1)=sigma2
                div_mul(2)=-sigma2
                XC(1,j,2,layer)=&
                    -ci*mu(layer+1)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=4+j,index=2,number_of_layer=layer+1)*div_mul(j)
                
                !sigma2=get_sigma(alpha,beta,3,number_of_layer=layer)
                !XC(1,j,1,layer)=&
                !    ci*mu(layer)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=4+j,index=4,number_of_layer=layer)
                !
                !sigma2=get_sigma(alpha,beta,3,number_of_layer=layer+1)
                !XC(1,j,2,layer)=&
                !    -ci*mu(layer+1)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=4+j,index=4,number_of_layer=layer+1)
            enddo
            
            sigma2=get_sigma(alpha,beta,3,number_of_layer=layer)
            XC(1,1,1,layer)=XC(1,1,1,layer)*exp(-sigma2*h(layer))
            !XC(1,2,1,layer)=XC(1,2,1,layer)*exp(-sigma2*0d0)
            
            sigma2=get_sigma(alpha,beta,3,number_of_layer=layer+1)
            !XC(1,1,2,layer)=XC(1,1,2,layer)*exp(sigma2*0d0)
            XC(1,2,2,layer)=XC(1,2,2,layer)*exp(-sigma2*h(layer+1))
            
            
            do j=1,2
                sigma2=get_sigma(alpha,beta,3,number_of_layer=layer)
                XC(2,j,1,layer)=&
                    get_eigenvector_component(alpha,beta,number_of_sigma=4+j,index=2,number_of_layer=layer)
                
                sigma2=get_sigma(alpha,beta,3,number_of_layer=layer+1)
                XC(2,j,2,layer)=&
                    -get_eigenvector_component(alpha,beta,number_of_sigma=4+j,index=2,number_of_layer=layer+1)
                
                !sigma2=get_sigma(alpha,beta,3,number_of_layer=layer)
                !XC(2,j,1,layer)=&
                !    get_eigenvector_component(alpha,beta,number_of_sigma=4+j,index=3,number_of_layer=layer)
                !
                !sigma2=get_sigma(alpha,beta,3,number_of_layer=layer+1)
                !XC(2,j,2,layer)=&
                !    -get_eigenvector_component(alpha,beta,number_of_sigma=4+j,index=3,number_of_layer=layer+1)
            enddo
            
            sigma2=get_sigma(alpha,beta,3,number_of_layer=layer)
            XC(2,1,1,layer)=XC(2,1,1,layer)*exp(-sigma2*h(layer))
            !XC(2,2,1,layer)=XC(2,2,1,layer)*exp(-sigma2*0d0)
            
            sigma2=get_sigma(alpha,beta,3,number_of_layer=layer+1)
            !XC(2,1,2,layer)=XC(2,1,2,layer)*exp(sigma2*0d0)
            XC(2,2,2,layer)=XC(2,2,2,layer)*exp(-sigma2*h(layer+1))
        enddo
        
        sigma2=get_sigma(alpha,beta,3,number_of_layer=number_of_layers)
        if(get_down_border_condition_type()==get_down_border_condition_type_fixed_border()) then
            do j=1,2
                XS(1,j,2)=&
                    +get_eigenvector_component(alpha,beta,number_of_sigma=4+j,index=2,number_of_layer=number_of_layers)
                
                !XS(1,j,2)=&
                !    +get_eigenvector_component(alpha,beta,number_of_sigma=4+j,index=3,number_of_layer=number_of_layers)
            enddo
        elseif(get_down_border_condition_type()==get_down_border_condition_type_free_border()) then
            div_mul(1)=sigma2
            div_mul(2)=-sigma2
            do j=1,2
                XS(1,j,2)=&
                    ci*mu(number_of_layers)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=4+j,index=2,number_of_layer=number_of_layers)*div_mul(j)
                !XS(1,j,2)=&
                !    ci*mu(number_of_layers)*alpha2*get_eigenvector_component(alpha,beta,number_of_sigma=4+j,index=4,number_of_layer=number_of_layers)
            enddo
        elseif(get_down_border_condition_type()==get_down_border_condition_type_halfspace()) then
            XS(1,2,2)=1d0
        else
            call print_error("count_K.generate_matrics_isotropic","a current get_down_border_condition_type for isotropic is not supported yet")
        endif
        XS(1,1,2)=XS(1,1,2)*exp(-sigma2*h(number_of_layers))
        !XS(1,2,2)=XS(1,2,2)*exp(-sigma20d0)
    endsubroutine generate_matrics_isotropic
    subroutine solve_matrics_system_isotropic(number_of_en)
    use main_parameters,only:number_of_layers,get_anisotropic
    use linear_system,only:star5_
    use system,only:print_error
    implicit none
        integer(4),intent(in)::number_of_en
        
        integer(4) i,j,layer
        
        if(number_of_en<1.or.number_of_en>3) call print_error("count_K.solve_matrics_system_isotropic","number_of_en<1.or.number_of_en>3")
        
        do i=1,number_of_layers*4
            BY(i,1)=YF(i)
        enddo
        
        call generate_summary_matrix_AY()
        
        call star5_(AY,BY)
        
        j=1
        do layer=1,number_of_layers
            do i=1,4
                Yt(i,layer,number_of_en)=BY(j,1)
                j=j+1
            enddo
        enddo
        
        do i=1,number_of_layers*2
            BX(i,1)=XF(i)
        enddo
        
        call generate_summary_matrix_AX()
        
        call star5_(AX,BX)
        
        j=1
        do layer=1,number_of_layers
            do i=1,2
                Xt(i,layer,number_of_en)=BX(j,1)
                j=j+1
            enddo
        enddo
    endsubroutine solve_matrics_system_isotropic
    subroutine generate_summary_matrix_AY()
    use main_parameters,only:number_of_layers
    implicit none
        integer(4) i,j,layer
        
        do i=1,number_of_layers*4
            do j=1,number_of_layers*4
                AY(i,j)=0d0
            enddo
        enddo
        
        do i=1,2
            do j=1,4
                AY(i,j)=YS(i,j,1)
            enddo
        enddo
        do layer=1,number_of_layers-1
            do i=1,4
                do j=1,4
                    AY(2+(layer-1)*4+i,4*(layer-1)+j)=YC(i,j,1,layer)
                    AY(2+(layer-1)*4+i,4*layer+j)=YC(i,j,2,layer)
                enddo
            enddo
        enddo
        do i=1,2
            do j=1,4
                AY(number_of_layers*4-2+i,number_of_layers*4-4+j)=YS(i,j,2)
            enddo
        enddo
    endsubroutine generate_summary_matrix_AY
    subroutine generate_summary_matrix_AX()
    use main_parameters,only:number_of_layers
    implicit none
        integer(4) i,j,layer
        
        do i=1,number_of_layers*2
            do j=1,number_of_layers*2
                AX(i,j)=0d0
            enddo
        enddo
        
        do i=1,1
            do j=1,2
                AX(i,j)=XS(i,j,1)
            enddo
        enddo
        do layer=1,number_of_layers-1
            do i=1,2
                do j=1,2
                    AX(1+(layer-1)*2+i,2*(layer-1)+j)=XC(i,j,1,layer)
                    AX(1+(layer-1)*2+i,2*layer+j)=XC(i,j,2,layer)
                enddo
            enddo
        enddo
        do i=1,1
            do j=1,2
                AX(number_of_layers*2-1+i,number_of_layers*2-2+j)=XS(i,j,2)
            enddo
        enddo
    endsubroutine generate_summary_matrix_AX
    
    
    
    complex(8) function get_det_A(alpha,beta,number_of_en) result(f)
    use main_parameters,only:number_of_layers,get_anisotropic
    use matrix_complex8,only:det
    use system,only:print_error
    implicit none
        complex(8),intent(in)::alpha
        complex(8),intent(in)::beta
        integer(4),intent(in)::number_of_en
        
        if(number_of_en<1.or.number_of_en>3) call print_error("count_K.get_det_A","number_of_en<1.or.number_of_en>3")
        if(get_anisotropic()==0) call print_error("count_K.get_det_A","get_anisotropic()==0")
        
        call generate_matrics(alpha,beta,number_of_en)
        call generate_summary_matrix_A()
        
        f=det(A,6*number_of_layers)
    endfunction get_det_A
    complex(8) function get_det_AY(alpha,beta,number_of_en) result(f)
    use main_parameters,only:number_of_layers,get_anisotropic
    use matrix_complex8,only:det
    use system,only:print_error
    implicit none
        complex(8),intent(in)::alpha
        complex(8),intent(in)::beta
        integer(4),intent(in)::number_of_en
        
        if(number_of_en<1.or.number_of_en>3) call print_error("count_K.get_det_AY","number_of_en<1.or.number_of_en>3")
        if(get_anisotropic()==1) call print_error("count_K.get_det_AY","get_anisotropic()==1")
        
        call generate_matrics(alpha,beta,number_of_en)
        call generate_summary_matrix_AY()
        
        f=det(AY,4*number_of_layers)
    endfunction get_det_AY
    complex(8) function get_det_AX(alpha,beta,number_of_en) result(f)
    use main_parameters,only:number_of_layers,get_anisotropic
    use matrix_complex8,only:det
    use system,only:print_error
    implicit none
        complex(8),intent(in)::alpha
        complex(8),intent(in)::beta
        integer(4),intent(in)::number_of_en
        
        if(number_of_en<1.or.number_of_en>3) call print_error("count_K.get_det_AX","number_of_en<1.or.number_of_en>3")
        if(get_anisotropic()==1) call print_error("count_K.get_det_AX","get_anisotropic()==1")
        
        call generate_matrics(alpha,beta,number_of_en)
        call generate_summary_matrix_AX()
        
        f=det(AX,2*number_of_layers)
    endfunction get_det_AX
    
    
    
    subroutine init_count_K()
    use main_parameters,only:number_of_layers,get_anisotropic
    implicit none
        if(get_anisotropic()==0) then
            allocate(Yt(4,number_of_layers,3))
            allocate(YC(4,4,2,number_of_layers-1))
            allocate(YF(4*number_of_layers))
            allocate(Xt(2,number_of_layers,3))
            allocate(XC(2,2,2,number_of_layers-1))
            allocate(XF(2*number_of_layers))
            
            !todo delete
            allocate(AY(number_of_layers*4,number_of_layers*4))
            allocate(BY(number_of_layers*4,1))
            allocate(AX(number_of_layers*2,number_of_layers*2))
            allocate(BX(number_of_layers*2,1))
        else
            allocate(t(6,number_of_layers,3))
            allocate(C(6,6,2,number_of_layers-1))
            allocate(F(number_of_layers*6))
            
            !todo delete
            allocate(A(number_of_layers*6,number_of_layers*6))
            allocate(B(number_of_layers*6,1))
        endif
    endsubroutine init_count_K
    subroutine destructor_count_K()
    use main_parameters,only:get_anisotropic
    implicit none
        if(get_anisotropic()==0) then
            deallocate(Yt)
            deallocate(YC)
            deallocate(YF)
            deallocate(Xt)
            deallocate(XC)
            deallocate(XF)
            
            !todo delte
            deallocate(AY)
            deallocate(BY)
            deallocate(AX)
            deallocate(BX)
        else
            deallocate(F)
            deallocate(C)
            deallocate(t)
            
            !todo delete
            deallocate(B)
            deallocate(A)
        endif
    endsubroutine destructor_count_K
endmodule count_K
