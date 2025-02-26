module integral_solution_isotropic
implicit none
    public::uz,uz_omega
    public::init_integral_solution_isotropic,destructor_integral_solution_isotropic
    
    real(8) r_,z_,t_
    
    private::r_,z_,t_,uz_inner_fun1,uz_inner_fun2
    contains
    
    complex(8) function uz(x,y,z,t) result(f)
    use main_parameters,only:get_anisotropic
    use GK_integration,only:GK_integral_ab
    use math,only:pi
    use system,only:print_error
    implicit none
        real(8),intent(in)::x
        real(8),intent(in)::y
        real(8),intent(in)::z
        real(8),intent(in)::t
        
        if(.not.get_anisotropic()==0) call print_error("integral_solution_isotropic.uz","material is not isotropic")
        if(z>0d0) call print_error("integral_solution_isotropic.uz","z>0d0")
        
        r_=sqrt(x*x+y*y)
        z_=z
        t_=t
        
        !complex(8) function GK_integral_ab(functionName,accurate,a,b,upperPolesValue,depthOfAvoidingPoles)
        !f=GK_integral_ab(uz_inner_fun1,1d-4,-5d0,5d0,0d0,0d0)
        f=GK_integral_ab(uz_inner_fun1,1d-4,-40d0,40d0,0d0,0d0)
        f=f/sqrt(pi+pi)
    endfunction uz
    
    complex(8) function uz_inner_fun1(omega_) result(f)
    use main_parameters,only:omega,Qomega
    use GK_integration,only:GK_integral_ab
    use math,only:ci
    implicit none
        complex(8),intent(in)::omega_
        
        omega=real(omega_)
        !complex(8) function GK_integral_ab(functionName,accurate,a,b,upperPolesValue,depthOfAvoidingPoles)
        !f=GK_integral_ab(uz_inner_fun2,1d-4,0d0,10d0,3d0,1d-2)
        f=GK_integral_ab(uz_inner_fun2,1d-4,0d0,300d0,250d0,1d-3)
        f=f*Qomega(omega_)*exp(-ci*omega*t_)
    endfunction uz_inner_fun1
    complex(8) function uz_inner_fun2(alpha) result(f)
    use main_parameters,only:Q
    use count_K,only:K
    use Jn,only:J0
    use math,only:c0
    implicit none
        complex(8),intent(in)::alpha
        
        f=K(3,3,alpha,c0,z_)*Q(3,alpha,(0d0,0d0))*alpha*J0(alpha*r_)
    endfunction uz_inner_fun2
    
    complex(8) function uz_omega(x,y,z,omega_) result(f)
    use main_parameters,only:omega,get_anisotropic
    use GK_integration,only:GK_integral_ab
    use math,only:pi,c0
    use system,only:print_error
    implicit none
        real(8),intent(in)::x
        real(8),intent(in)::y
        real(8),intent(in)::z
        real(8),intent(in)::omega_
        
        if(.not.get_anisotropic()==0) call print_error("integral_solution_isotropic.uz_omega","material is not isotropic")
        if(z>0d0) call print_error("integral_solution_isotropic.uz_omega","z>0d0")
        
        omega=omega_
        r_=sqrt(x*x+y*y)
        z_=z
        
        !complex(8) function GK_integral_ab(functionName,accurate,a,b,upperPolesValue,depthOfAvoidingPoles)
        f=GK_integral_ab(uz_inner_fun2,1d-5,0d0,300d0,250d0,1d-4)
    endfunction uz_omega
    
    subroutine init_integral_solution_isotropic()
    implicit none
    endsubroutine init_integral_solution_isotropic
    subroutine destructor_integral_solution_isotropic()
    implicit none
    endsubroutine destructor_integral_solution_isotropic
endmodule integral_solution_isotropic
