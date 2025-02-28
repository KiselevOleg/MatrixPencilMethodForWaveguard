module pre_and_post_runtime_actions
implicit none
    public::init,destructor
    public::h2,rho2,E2
    
    real(8),save::h2=1d0,rho2=3.168d0,E2=33d0
    
    private::init_main_parameters,init_main_parameters_isotropic,init_main_parameters_anisotropic,destructor_main_parameters
    contains
    
    subroutine init()
    use sigma_and_eigenvectors,only:init_sigma_and_eigenvectors
    use count_K,only:init_count_K
    use integral_solution,only:init_integral_solution
    use integral_solution_isotropic,only:init_integral_solution_isotropic
    use load_experimental_measurements,only:init_load_experimental_measurements
    implicit none
        call init_main_parameters
        
        call init_sigma_and_eigenvectors
        call init_count_K
        call init_integral_solution
        call init_integral_solution_isotropic
        call init_load_experimental_measurements
    endsubroutine init
    subroutine destructor
    use sigma_and_eigenvectors,only:destructor_sigma_and_eigenvectors
    use count_K,only:destructor_count_K
    use integral_solution,only:destructor_integral_solution
    use integral_solution_isotropic,only:destructor_integral_solution_isotropic
    use load_experimental_measurements,only:destructor_load_experimental_measurements
    implicit none
        call destructor_main_parameters
        
        call destructor_sigma_and_eigenvectors
        call destructor_count_K
        call destructor_integral_solution
        call destructor_integral_solution_isotropic
        call destructor_load_experimental_measurements
    endsubroutine destructor
    
    subroutine init_main_parameters()
    use main_parameters,only:number_of_layers,h,Calphabeta,rho,omega,Q,E,nu,lambda,mu,Cp,Cs,&
        set_anisotropic,get_anisotropic,set_down_border_condition_type,get_down_border_condition_type,&
        get_down_border_condition_type_fixed_border,get_down_border_condition_type_free_border,&
        get_down_border_condition_type_halfspace
    use system,only:print_error,print_warning,print_error
    use math,only:epsilon
    implicit none
        integer(4) i,j,k
        complex(8) sum
        
        !description of main_parameters.Q(res,alpha,beta)
        
        call set_anisotropic(0)
        !call set_down_border_condition_type(get_down_border_condition_type_fixed_border())
        call set_down_border_condition_type(get_down_border_condition_type_free_border())
        !call set_down_border_condition_type(get_down_border_condition_type_halfspace())
        
        number_of_layers=1
        allocate(h(number_of_layers))
        allocate(rho(number_of_layers))
        
        omega=3d0
        
        h(1)=2.8816d0
        rho(1)=7.743d0
        
        do k=1,number_of_layers
            if(rho(k)<=epsilon) call print_error("pre_and_post_runtime_actions.init_main_parameters","it exists k that rho(k)<=epsilon")
            if(rho(k)>=10d0) call print_warning("pre_and_post_runtime_actions.init_main_parameters","it looks like rho(k) is too large")
            if(h(k)<=epsilon) call print_error("pre_and_post_runtime_actions.init_main_parameters","it exists k that h(k)<=epsilon")
            if(h(k)<=0.001d0) call print_warning("pre_and_post_runtime_actions.init_main_parameters","it looks like h(k) is too small")
            if(h(k)>=10d0.and..not.(&
                get_down_border_condition_type()==get_down_border_condition_type_halfspace()&
                .and.k==number_of_layers&
                )) call print_warning("pre_and_post_runtime_actions.init_main_parameters","it looks like h(k) is too large")
        enddo
        if(omega<=epsilon) call print_error("pre_and_post_runtime_actions.init_main_parameters","omega<=epsilon")
        
        if(get_down_border_condition_type()==get_down_border_condition_type_halfspace()) then
            if(h(number_of_layers)<100d0) then
                call print_error("pre_and_post_runtime_actions.init_main_parameters",&
                    "get_down_border_condition_type()==get_down_border_condition_type_halfspace().and.h(number_of_layers)<100d0"&
                )
            endif
        endif
        
        if(get_anisotropic()==0) then
            call init_main_parameters_isotropic
        else if(get_anisotropic()==1) then
            call init_main_parameters_anisotropic
        else
            call print_warning("pre_and_post_runtime_actions.init_main_parameters","an incorreact anisotropic type")
        endif
    endsubroutine init_main_parameters
    subroutine destructor_main_parameters()
    use main_parameters,only:number_of_layers,h,Calphabeta,rho,omega,Q,E,nu,lambda,mu,Cp,Cs,get_anisotropic
    implicit none
        deallocate(rho)
        if(get_anisotropic()==1) then
            deallocate(Calphabeta)
        else
            deallocate(E)
            deallocate(nu)
            deallocate(lambda)
            deallocate(mu)
            deallocate(Cp)
            deallocate(Cs)
        endif
        deallocate(h)
    endsubroutine destructor_main_parameters
    
    subroutine init_main_parameters_isotropic()
    use main_parameters,only:number_of_layers,rho,E,nu,lambda,mu,Cp,Cs
    use system,only:print_error,print_warning
    use math,only:epsilon
    implicit none
        integer(4) k
        
        allocate(E(number_of_layers))
        allocate(nu(number_of_layers))
        allocate(lambda(number_of_layers))
        allocate(mu(number_of_layers))
        allocate(Cp(number_of_layers))
        allocate(Cs(number_of_layers))
        
        do k=1,number_of_layers
            E(k)=0d0
            nu(k)=0d0
            lambda(k)=0d0
            mu(k)=0d0
            Cp(k)=0d0
            Cs(k)=0d0
        enddo
        
        !lambda(1)=7.14285d0
        !mu(1)=1.78571d0
        !~
        !E(2)=5d0
        !nu(2)=0.4d0
        !~
        !Cp(3)=2.314550249431379d0
        !Cs(3)=0.944911182523068d0
        
        E(1)=2.06d0
        nu(1)=0.31d0
        
        do k=1,number_of_layers
            if(E(k)>epsilon) then
                if(lambda(k)>epsilon.or.Cp(k)>epsilon) call print_error("pre_and_post_runtime_actions.init_main_parameters_isotropic",&
                    "a multi initializated layer")
                
                lambda(k)=nu(k)*E(k)/(1d0+nu(k))/(1d0-2d0*nu(k))
                mu(k)=E(k)*0.5d0/(1d0+nu(k))
                
                Cp(k)=sqrt(E(k)/rho(k))*sqrt((1d0-nu(k))/(1d0+nu(k))/(1d0-2d0*nu(k)))
                Cs(k)=sqrt(E(k)/rho(k))*sqrt(0.5d0/(1d0+nu(k)))
            else if(lambda(k)>epsilon) then
                if(E(k)>epsilon.or.Cp(k)>epsilon) call print_error("pre_and_post_runtime_actions.init_main_parameters_isotropic",&
                    "a multi initializated layer")
                
                E(k)=mu(k)*(3d0*lambda(k)+2d0*mu(k))/(lambda(k)+mu(k))
                nu(k)=lambda(k)*0.5d0/(lambda(k)+mu(k))
                
                Cp(k)=sqrt((lambda(k)+mu(k)+mu(k))/rho(k))
                Cs(k)=sqrt(mu(k)/rho(k))
            else if(Cp(k)>epsilon) then
                if(lambda(k)>epsilon.or.E(k)>epsilon) call print_error("pre_and_post_runtime_actions.init_main_parameters_isotropic",&
                    "a multi initializated layer")
                
                lambda(k)=rho(k)*Cp(k)*Cp(k)-2d0*rho(k)*Cs(k)*Cs(k)
                mu(k)=rho(k)*Cs(k)*Cs(k)
                
                E(k)=Cs(k)*Cs(k)*rho(k)*2d0*(1d0+(2d0*Cs(k)*Cs(k)-Cp(k)*Cp(k))/(2d0*Cs(k)*Cs(k)-2d0*Cp(k)*Cp(k)))
                nu(k)=(2d0*Cs(k)*Cs(k)-Cp(k)*Cp(k))/(2d0*Cs(k)*Cs(k)-2d0*Cp(k)*Cp(k))
            else
                call print_error("pre_and_post_runtime_actions.init_main_parameters_isotropic","a not initializated layer")
            endif
        enddo
        
        do k=1,number_of_layers
            if(lambda(k)<epsilon.or.lambda(k)>1d3) call print_error("pre_and_post_runtime_actions.init_main_parameters_isotropic",&
                "lambda(k)<epsilon.or.lambda(k)>1d3")
            if(mu(k)<epsilon.or.mu(k)>1d3) call print_error("pre_and_post_runtime_actions.init_main_parameters_isotropic",&
                "mu(k)<epsilon.or.mu(k)>1d3")
            
            if(E(k)<epsilon.or.E(k)>1d3) call print_error("pre_and_post_runtime_actions.init_main_parameters_isotropic",&
                "E(k)<epsilon.or.E(k)>1d3")
            if(nu(k)<epsilon.or.nu(k)>0.5d0-epsilon) call print_error("pre_and_post_runtime_actions.init_main_parameters_isotropic",&
                "nu(k)<epsilon.or.nu(k)>0.5d0-epsilon")
            
            if(Cp(k)<epsilon.or.Cp(k)>1d3) call print_error("pre_and_post_runtime_actions.init_main_parameters_isotropic",&
                "Cp(k)<epsilon.or.Cp(k)>1d3")
            if(Cs(k)<epsilon.or.Cs(k)>1d3) call print_error("pre_and_post_runtime_actions.init_main_parameters_isotropic",&
                "Cs(k)<epsilon.or.Cs(k)>1d3")
            
            if(Cp(k)<Cs(k)) call print_error("pre_and_post_runtime_actions.init_main_parameters_isotropic","Cp(k)<Cs(k)")
        enddo
    endsubroutine init_main_parameters_isotropic
    
    subroutine init_main_parameters_anisotropic()
    use main_parameters,only:number_of_layers,Calphabeta
    use system,only:print_error,print_warning
    use math,only:epsilon
    implicit none
        complex(8) E(3),nu(3),lambda(3),mu(3)
        integer(4) i,j,k
        complex(8) sum
        
        allocate(Calphabeta(6,6,number_of_layers))
        
        do k=1,number_of_layers
            do i=1,6
                do j=1,6
                    Calphabeta(i,j,k)=0d0
                enddo
            enddo
        enddo
        
        !h(1)=1d0
        !rho(1)=8.905d0
        !Calphabeta(1,1,1)=25d0
        !Calphabeta(1,2,1)=16d0
        !Calphabeta(1,3,1)=16d0
        !Calphabeta(2,2,1)=25d0
        !Calphabeta(2,3,1)=16d0
        !Calphabeta(3,3,1)=25d0
        !Calphabeta(4,4,1)=11.85d0
        !Calphabeta(5,5,1)=11.85d0
        !Calphabeta(6,6,1)=11.85d0
        
        
        !h(1)=1d0
        !rho(1)=1.578d0
        !Calphabeta(1,1,1)=130.7d0
        !Calphabeta(1,2,1)=5.2d0
        !Calphabeta(1,3,1)=5.2d0
        !Calphabeta(2,2,1)=13d0
        !Calphabeta(2,3,1)=4.5d0
        !Calphabeta(3,3,1)=13d0
        !Calphabeta(4,4,1)=13.7d0
        !Calphabeta(5,5,1)=6d0
        !Calphabeta(6,6,1)=6d0
                
        !h(1)=0.1d0*10
        !rho(1)=2d0
        Calphabeta(1,1,1)=10d0
        Calphabeta(1,2,1)=6d0
        Calphabeta(1,3,1)=6d0
        Calphabeta(2,2,1)=10d0
        Calphabeta(2,3,1)=6d0
        Calphabeta(3,3,1)=10d0
        Calphabeta(4,4,1)=2d0
        Calphabeta(5,5,1)=2d0
        Calphabeta(6,6,1)=2d0
        
        Calphabeta(1,1,2)=11d0
        Calphabeta(1,2,2)=5d0
        Calphabeta(1,3,2)=5d0
        Calphabeta(2,2,2)=11d0
        Calphabeta(2,3,2)=5d0
        Calphabeta(3,3,2)=11d0
        Calphabeta(4,4,2)=3d0
        Calphabeta(5,5,2)=3d0
        Calphabeta(6,6,2)=3d0
        
        Calphabeta(1,1,3)=6d0
        Calphabeta(1,2,3)=4d0
        Calphabeta(1,3,3)=4d0
        Calphabeta(2,2,3)=6d0
        Calphabeta(2,3,3)=4d0
        Calphabeta(3,3,3)=6d0
        Calphabeta(4,4,3)=1d0
        Calphabeta(5,5,3)=1d0
        Calphabeta(6,6,3)=1d0
        
        E(1)=71d0
        nu(1)=0.33d0
        E(2)=4.76d0*7
        !E(2)=4.76d0*7
        nu(2)=0.25d0
        E(3)=131d0
        nu(3)=0.266d0
        
        do i=1,3
            lambda(i)=nu(i)*E(i)/(1d0+nu(i))/(1d0-2d0*nu(i))+(0d0,1d-2)
            mu(i)=E(i)*0.5d0/(1d0+nu(i))+(0d0,1d-2)
            
            Calphabeta(1,1,i)=lambda(i)+2*mu(i)
            Calphabeta(1,2,i)=lambda(i)
            Calphabeta(1,3,i)=lambda(i)
            Calphabeta(2,2,i)=lambda(i)+2*mu(i)
            Calphabeta(2,3,i)=lambda(i)
            Calphabeta(3,3,i)=lambda(i)+2*mu(i)
            Calphabeta(4,4,i)=mu(i)
            Calphabeta(5,5,i)=mu(i)
            Calphabeta(6,6,i)=mu(i)
        enddo
        
        !h(1)=1d0
        !rho(1)=1d0
        !Calphabeta(1,1,1)=119.14d0
        !Calphabeta(1,2,1)=3.57d0
        !Calphabeta(1,3,1)=3.57d0
        !Calphabeta(2,2,1)=10.56d0
        !Calphabeta(2,3,1)=4.96d0
        !Calphabeta(3,3,1)=10.56d0
        !Calphabeta(4,4,1)=2.8d0
        !Calphabeta(5,5,1)=4.45d0
        !Calphabeta(6,6,1)=4.45d0
        
        !h(1)=0.56d0
        !rho(1)=1.522d0
        !Calphabeta(1,1,1)=119.14d0
        !Calphabeta(1,2,1)=3.57d0
        !Calphabeta(1,3,1)=3.57d0
        !Calphabeta(2,2,1)=10.56d0
        !Calphabeta(2,3,1)=4.96d0
        !Calphabeta(3,3,1)=10.56d0
        !Calphabeta(4,4,1)=2.8d0
        !Calphabeta(5,5,1)=4.45d0
        !Calphabeta(6,6,1)=4.45d0
        !
        !h(2)=1.12d0
        !rho(2)=1.522d0
        !Calphabeta(1,1,2)=10.56d0
        !Calphabeta(1,2,2)=3.57d0
        !Calphabeta(1,3,2)=4.96d0
        !Calphabeta(2,2,2)=119.14d0
        !Calphabeta(2,3,2)=3.57d0
        !Calphabeta(3,3,2)=10.56d0
        !Calphabeta(4,4,2)=4.45d0
        !Calphabeta(5,5,2)=2.8d0
        !Calphabeta(6,6,2)=4.45d0
        !
        !h(3)=0.56d0
        !rho(3)=1.522d0
        !Calphabeta(1,1,3)=119.14d0
        !Calphabeta(1,2,3)=3.57d0
        !Calphabeta(1,3,3)=3.57d0
        !Calphabeta(2,2,3)=10.56d0
        !Calphabeta(2,3,3)=4.96d0
        !Calphabeta(3,3,3)=10.56d0
        !Calphabeta(4,4,3)=2.8d0
        !Calphabeta(5,5,3)=4.45d0
        !Calphabeta(6,6,3)=4.45d0
        
        !h(1)=1d0
        !rho(1)=8.932d0
        !Calphabeta(1,1,1)=222d0
        !Calphabeta(1,2,1)=71d0
        !Calphabeta(1,3,1)=123d0
        !Calphabeta(2,2,1)=222d0
        !Calphabeta(2,3,1)=123d0
        !Calphabeta(3,3,1)=170d0
        !Calphabeta(4,4,1)=75.5d0
        !Calphabeta(5,5,1)=75.5d0
        !Calphabeta(6,6,1)=23.5d0
        
        do k=1,number_of_layers
            do i=2,6
                do j=1,i
                    Calphabeta(i,j,k)=Calphabeta(j,i,k)
                enddo
            enddo
        enddo
        
        do k=1,number_of_layers
            sum=0d0
            do i=1,6
                do j=1,6
                    sum=sum+Calphabeta(i,j,k)
                    if(abs(Calphabeta(i,j,k))>300d0) then
                        call print_warning("pre_and_post_runtime_actions.init_main_parameters_anisotropic",&
                            "it looks like Calphabeta(i,j,k) is too large")
                    endif
                    if(abs(Calphabeta(i,j,k))>=epsilon.and.abs(Calphabeta(i,j,k))<1d0) then
                        call print_warning("pre_and_post_runtime_actions.init_main_parameters_anisotropic",&
                            "it looks like Calphabeta(i,j,k) is too small")
                    endif
                enddo
            enddo
            if(abs(sum)<epsilon) call print_error("pre_and_post_runtime_actions.init_main_parameters_anisotropic",&
                "Calphabeta not fully fulled")
        enddo
    endsubroutine init_main_parameters_anisotropic
endmodule pre_and_post_runtime_actions
