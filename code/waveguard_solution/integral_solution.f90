module integral_solution
implicit none
    public::u    
    
    public::init_integral_solution,destructor_integral_solution
    
    private::Q,K
    contains
    
    complex(8) function u(ind,x,y,z,accurate, lenght_of_integration,upper_poles_value,depthOfAvoidingPoles) result(f)
    use GK_integration,only:GK_integral
    use math,only:pi,epsilon
    use system,only:print_error
    implicit none
        integer(4),intent(in)::ind
        real(8),intent(in)::x
        real(8),intent(in)::y
        real(8),intent(in)::z
        real(8),intent(in)::accurate!1d-5
        real(8),intent(in)::lenght_of_integration!20d0
        real(8),intent(in)::upper_poles_value!4d0+4d0
        real(8),intent(in)::depthOfAvoidingPoles!1d-1
        
        complex(8) alpha_
        real(8) R
        
        if(ind<1.or.ind>3) call print_error("integral_solution.u","ind<1.or.ind>3")
        if(z>0d0) call print_error("integral_solution.u","z>0d0")
        
        R=sqrt(x*x+y*y+z*z)
        R=1d0
        
        f=0.25d0/R/(pi*pi)*GK_integral(functionName=external_integral_function,&
            accurate=accurate,lengthOfIntegration=lenght_of_integration,upperPolesValue=upper_poles_value,depthOfAvoidingPoles=depthOfAvoidingPoles)
    contains
        complex(8) function external_integral_function(alpha) result(f)
        use GK_integration,only:GK_integral
        implicit none
            complex(8),intent(in)::alpha
            
            alpha_=alpha
            f=alpha_*GK_integral(functionName=inner_integral_function,&
                accurate=accurate,lengthOfIntegration=pi+pi,upperPolesValue=pi,depthOfAvoidingPoles=epsilon)
        endfunction external_integral_function
        complex(8) function inner_integral_function(gamma) result(f)
        use math,only:ci
        implicit none
            complex(8),intent(in)::gamma
            
            complex(8) Qv
            integer(4) j
            
            f=0d0
            do j=1,3
                Qv=Q(j,alpha_,gamma)
                if(abs(Qv)<epsilon) cycle
                
                f=f+R*K(ind,j,alpha_,gamma,z)*Q(j,alpha_,gamma)*exp(-ci*alpha_*(x*sin(gamma)+y*cos(gamma)))
            enddo
        endfunction inner_integral_function
    endfunction u
    
    complex(8) function Q(ind,alpha,gamma) result(f)
    use main_parameters,only:Q_=>Q
    use math,only:pi,from_polar_to_cartesian2_coordinate_system
    use system,only:print_error
    implicit none
        integer(4),intent(in)::ind
        complex(8),intent(in)::alpha
        complex(8),intent(in)::gamma
        
        complex(8) alpha1,alpha2
        
        if(ind<1.or.ind>3) call print_error("integral_solution.Q","ind<1.or.ind>3")
        if(real(alpha)<0d0) call print_error("integral_solution.Q","real(alpha)<0d0")
        if(real(gamma)<0d0.or.real(gamma)>2d0*pi) call print_error("integral_solution.Q","real(gamma)<0d0.or.real(gamma)>2d0*pi")
        
        call from_polar_to_cartesian2_coordinate_system(alpha1,alpha2,alpha,gamma)
        
        f=Q_(ind,alpha1,alpha2)
    endfunction Q
    
    complex(8) function K(i,j,alpha,gamma,z) result(f)
    use count_K,only:K_=>K
    use math,only:pi,from_polar_to_cartesian2_coordinate_system
    use system,only:print_error
    implicit none
        integer(4),intent(in)::i
        integer(4),intent(in)::j
        complex(8),intent(in)::alpha
        complex(8),intent(in)::gamma
        real(8),intent(in)::z
        
        complex(8) alpha1,alpha2
        
        if(i<1.or.i>3) call print_error("integral_solution.K","i<1.or.i>3")
        if(j<1.or.j>3) call print_error("integral_solution.K","j<1.or.j>3")
        if(real(alpha)<0d0) call print_error("integral_solution.K","real(alpha)<0d0")
        if(real(gamma)<0d0.or.real(gamma)>2d0*pi) call print_error("integral_solution.K","real(gamma)<0d0.or.real(gamma)>2d0*pi")
        
        call from_polar_to_cartesian2_coordinate_system(alpha1,alpha2,alpha,gamma)
        
        f=K_(i,j,alpha1,alpha2,z)
    endfunction K
    
    subroutine init_integral_solution()
    use main_parameters,only:number_of_layers
    implicit none
    endsubroutine init_integral_solution
    subroutine destructor_integral_solution()
    implicit none
    endsubroutine destructor_integral_solution
endmodule integral_solution
