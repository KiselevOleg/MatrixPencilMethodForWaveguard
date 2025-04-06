module main_parameters
implicit none
    public::get_anisotropic_type,set_anisotropic_type,
        material_isotropic_type,material_anisotropic_type
    public::set_down_border_condition_type,get_down_border_condition_type,&
        down_border_condition_type_fixed_border,down_border_condition_type_free_border,&
        down_border_condition_type_halfspace
    
    public::set_omega,get_omega,set_f,get_f
    
    public::set_Q,get_Q,set_Qomega,get_Qomega
    
    public::get_number_of_layers,set_number_of_layers
    
    public::get_layer_parameter
    
    public::get_layer_E,get_layer_nu,get_layer_lambda,get_layer_mu,get_layer_Cp,get_layer_Cs
    public::get_layer_rho,get_layer_h,get_full_h
    public::set_layer_E_nu,set_layer_lambda_mu,set_layer_Cp_Cs
    public::set_layer_rho,set_layer_h
    
    public::get_layer_Cijkl_parameter,get_layer_Calphabeta_parameter
    public:set_layer_Calphabeta
    
    public::count_isotropic_material_parameters_form_lambda_mu,&
        count_isotropic_material_parameters_form_E_nu,&
        count_isotropic_material_parameters_form_Cp_Cs
    public::count_isotropic_material_parameter_E_form_lambda_mu,count_isotropic_material_parameter_nu_form_lambda_mu,
        count_isotropic_material_parameter_Cp_form_lambda_mu,count_isotropic_material_parameter_Cs_form_lambda_mu,
        count_isotropic_material_parameter_lambda_form_E_nu,count_isotropic_material_parameter_mu_form_E_nu,
        count_isotropic_material_parameter_Cp_form_E_nu,count_isotropic_material_parameter_Cs_form_E_nu,
        count_isotropic_material_parameter_lambda_form_Cp_Cs,count_isotropic_material_parameter_mu_form_Cp_Cs,
        count_isotropic_material_parameter_E_form_Cp_Cs,count_isotropic_material_parameter_nu_form_Cp_Cs
    
    
    
    integer(4),private::anisotropic_type=-1
    integer(4),private::down_border_condition_type=-1
    
    real(8),private::omega
    
    complex(8),private,external::Q,Qomega
    logical(1),private:Q_setted=.false.,Qomega_setted=.false.
    
    integer(4),private::number_of_layers=-1
    real(8),private,allocatable::h(:)!number_of_layer
    real(8),private,allocatable::rho(:)!number_of_layer
    
    complex(8),private,allocatable::Calphabeta(:,:,:)!number_of_layer,alpha,beta
    
    complex(8),private,allocatable::E(:)!number_of_layer
    real(8),private,allocatable::nu(:)!number_of_layer
    complex(8),private,allocatable::lambda(:)!number_of_layer
    complex(8),private,allocatable::mu(:)!number_of_layer
    complex(8),private,allocatable::Cp(:)!number_of_layer
    complex(8),private,allocatable::Cs(:)!number_of_layer
    
    logical(1),private::matrixes_created=.false.
    
    
    
    private::free_matrixes,create_matrixes
    contains
    
    pure integer(4) function get_anisotropic_type() result(f)
    use system,only::print_error
    implicit none
        if(anisotropic_type==-1) &
            call print_error("main_parameters.get_anisotropic_type","anisotripoc_type is not setted")
        
        f=anisotropic_type
    endfunction get_anisotropic_type
    subroutine set_anisotropic_type(anisotropic_type_) result(f)
    use system,only::print_error
    implicit none
        integer(4),intent(in)::anisotropic_type_
        
        if(.not.anisotropic_type==-1) &
            call print_error("main_parameters.set_anisotropic_type","anisotripoc_type has been already setted")
        if(.not.(anisotropic_type_==material_isotropic_type().or.anisotropic_type_==material_anisotropic_type())) &
            call print_error("main_parameters.set_anisotropic_type","incorrect anisotripoc_type")
        
        anisotropic_type=anisotropic_type_
    endsubroutine set_anisotropic_type
    pure integer(4) function material_isotropic_type() result(f)
    implicit none
        f=1
    endfunction material_isotropic_type
    pure integer(4) function material_anisotropic_type() result(f)
    implicit none
        f=2
    endfunction material_anisotropic_type
    
    
    
    subroutine set_down_border_condition_type(down_border_condition_type_)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::down_border_condition_type_
        
        if(.not.down_border_condition_type==-1) &
            call print_error("main_parameters.set_down_border_condition_type","down_border_condition_type has been already setted")
        if(.not.(&
            down_border_condition_type_==down_border_condition_type_fixed_border().or.&
            down_border_condition_type_==down_border_condition_type_free_border().or.&
            down_border_condition_type_==down_border_condition_type_halfspace().or.&
            )) call print_error("main_parameters.set_down_border_condition_type","incorrect down_border_condition_type_")
        
        down_border_condition_type=down_border_condition_type_
    endsubroutine set_down_border_condition_type
    pure integer(4) function get_down_border_condition_type() result(f)
    use system,only:print_error
    implicit none
        if(down_border_condition_type==-1) &
            call print_error("main_parameters.get_down_border_condition_type","down_border_condition_type is not setted")
        
        f=down_border_condition_type
    endfunction get_down_border_condition_type
    pure integer(4) function down_border_condition_type_fixed_border() result(f)
    implicit none
        f=1
    endfunction down_border_condition_type_fixed_border
    pure integer(4) function down_border_condition_type_free_border() result(f)
    implicit none
        f=2
    endfunction down_border_condition_type_free_border
    pure integer(4) function down_border_condition_type_halfspace() result(f)
    implicit none
        f=3
    endfunction down_border_condition_type_halfspace
    
    
    
    subroutine set_omega(omega_) result(f)
    implicit none
        real(8),intent(in)::omega_
        
        omega=omega_
    endsubroutine set_omega
    subroutine set_f(f_) result(f)
    use math,only:pi
    implicit none
        real(8),intent(in)::f_
        
        omega=f_*2d0*pi
    endsubroutine set_f
    pure real(8) function get_omega() result(f)
    implicit none
        f=omega
    endfunction get_omega
    pure real(8) function get_f() result(f)
    use math,only:pi
    implicit none
        f=omega*0.5d0/pi
    endfunction get_f
    
    
    
    subroutine set_Q(Q_)
    implicit none
        Q_setted=.true.
        Q=Q_
    endsubroutine set_Q
    pure complex(8) function get_Q(ind,alpha,beta) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::ind
        complex(8),intent(in)::alpha
        complex(8),intent(in)::beta
        
        if(.not.Q_setted) call print_error("main_parameters.set_Q","Q is not setted")
        if(.not.(1<=ind.and.ind<=3)) call print_error("main_parameters.set_Q","incorrect ind")
        
        f=Q(ind,alpha,beta)
    endfunction set_Q
    subroutine set_Qomega(Qomega_)
    implicit none
        Qomega_setted=.true.
        Qomega=Qomega_
    endsubroutine set_Q
    pure complex(8) function get_Qomega(omega) result(f)
    use system,only:print_error
    implicit none
        real(8),intent(in)::omega
        
        if(.not.Qomega_setted) call print_error("main_parameters.set_Qomega","Qomega is not setted")
        
        f=Qomega(omega)
    endfunction set_Q
    
    
    
    pure integer(4) function get_number_of_layers() result(f)
    use system,only:print_error
    implicit none
        if(number_of_layers==-1) call print_error("main_parameters.get_number_of_layers","number_of_layer is not setted")
        
        f=number_of_layers
    endfunction get_number_of_layers
    subroutine set_number_of_layers(number_of_layers_)
    use system,only:print_error,print_warning
    implicit none
        integer(4),intent(in)::set_number_of_layers_
        
        if(number_of_layers_<1) call print_error("mainr_parameters.set_number_of_layers","number_of_layers_<1")
        if(number_of_layers_>10) call print_warning("mainr_parameters.set_number_of_layers","number_of_layers_>10")
        
        call free_matrixes()
        call create_matrixes()
    endsubroutine set_number_of_layers
    
    subroutine free_matrixes()
    use system,only:print_error
    implicit none
        if(anisotropic_type==-1) call print_error("main_parameters.free_matrixes","anisotropic_type is not setted")
        if(.not.matrixes_created) call print_error("main_parameters.free_matrixes","matrixes is not created")
        
        matrixes_created=.false.
        
        deallocate(h)
        deallocate(rho)
        
        if(anisptripic_type==material_isotropic_type()) then
            deallocate(E)
            deallocate(nu)
            deallocate(lambda)
            deallocate(mu)
            deallocate(Cp)
            deallocate(Cs)
        elseif(anisotropic_type==material_anisotropic_type()) then
            deallocate(Calphabeta)
        endif
    endsubroutine free_matrixes
    subroutine create_matrixes()
    use system,only:print_error
    implicit none
        integer(4) i,j,k
        
        if(anisotropic_type==-1) call print_error("main_parameters.create_matrixes","anisotropic_type is not setted")
        if(matrixes_created) call print_error("main_parameters.create_matrixes","matrixes has been already created")
        if(number_Of_layers==-1) call print_error("main_parameters.create_matrixes","number_Of_layers is not setted")
        
        matrixes_created=.true.
        
        allocate(h(number_of_layers))
        allocate(rho(number_of_layers))
        do i=1,number_Of_layers
            h(i)=-1d0
            rho(i)=-1d0
        enddo
        
        if(anisptripic_type==material_isotropic_type()) then
            allocate(E(number_of_layers))
            allocate(nu(number_of_layers))
            allocate(lambda(number_of_layers))
            allocate(mu(number_of_layers))
            allocate(Cp(number_of_layers))
            allocate(Cs(number_of_layers))
            
            do i=1,number_of_layers
                E(i)=-1d0
                nu(i)=-1d0
                lambda(i)=-1d0
                mu(i)=-1d0
                Cp(i)=-1d0
                Cs(i)=-1d0
            enddo
        elseif(anisotropic_type==material_anisotropic_type()) then
            allocate(Calphabeta(number_of_layers,6,6))
            do i=1,number_Of_layers
                do j=1,6
                    do k=1,6
                        Calphabeta(i,j,k)=-1d0
                    enddo
                enddo
            enddo
        endif
    endsubroutine create_matrixes
    
    
    
    public::get_layer_parameter
    public::get_layer_E,get_layer_nu,get_layer_lambda,get_layer_mu,get_layer_Cp,get_layer_Cs
    public::get_layer_rho,get_layer_h,get_full_h
    public::set_layer_E_nu,set_layer_lambda_mu,set_layer_Cp_Cs
    public::set_layer_rho,set_layer_h
    public::get_layer_Cijkl_parameter,get_layer_Calphabeta_parameter
    public:set_layer_Calphabeta
endmodule main_parameters
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
module main_parameters_
implicit none
    public::number_of_layers,h,get_full_h,Calphabeta,rho,omega,Q,Qomega,get_number_of_layer,Cijkl
    public::E,nu,lambda,mu,Cp,Cs
    public::count_isotropic_material_parameters_form_lambda_mu,&
        count_isotropic_material_parameters_form_E_nu,&
        count_isotropic_material_parameters_form_Cp_Cs
    
    public::get_anisotropic,set_anisotropic
    public::set_down_border_condition_type,get_down_border_condition_type,&
        get_down_border_condition_type_fixed_border,get_down_border_condition_type_free_border,&
        get_down_border_condition_type_halfspace
    
    integer(4),save::anisotropic=-1
    integer(4),save::down_border_condition_type=-1
    
    integer(4),save::number_of_layers
    real(8),save,allocatable::h(:)!number_of_layer
    real(8),save,allocatable::rho(:)!number_of_layer
    
    complex(8),save,allocatable::Calphabeta(:,:,:)!alpha,beta,number_of_layer
    
    complex(8),save,allocatable::E(:)!number_of_layer
    real(8),save,allocatable::nu(:)!number_of_layer
    complex(8),save,allocatable::lambda(:)!number_of_layer
    complex(8),save,allocatable::mu(:)!number_of_layer
    complex(8),save,allocatable::Cp(:)!number_of_layer
    complex(8),save,allocatable::Cs(:)!number_of_layer
    
    real(8),save::omega
    
    real(8),save::full_h=-1d0
    private::full_h,anisotropic,down_border_condition_type
    contains
    
    complex(8) function Q(ind,alpha,beta) result(f)
    use math,only:pi
    use Jn,only:J2
    use system,only:print_error
    implicit none
        integer(4),intent(in)::ind
        complex(8),intent(in)::alpha
        complex(8),intent(in)::beta
        
        complex(8) alpha_
        alpha_=sqrt(alpha*alpha+beta*beta)
        
        if(ind<1.or.ind>3) call print_error("main_parameters.Q","ind<1.or.ind>3")
        
        if(ind==1) then
            f=0d0
        elseif(ind==2) then
            f=0d0
        else
            f=1d0
        endif
    endfunction Q
    complex(8) function Qomega(omega) result(f)
    use math,only:pi
    use system,only:print_error
    implicit none
        complex(8),intent(in)::omega
        
        f=1d0
    endfunction Qomega
    
    real(8) function get_full_h() result(f)
    implicit none
        integer(4) i
        
        if(full_h>0d0) then
            f=full_h
            return
        endif
        
        full_h=0d0
        do i=1,number_of_layers
            full_h=full_h+h(i)
        enddo
        
        f=full_h
        
        if(get_down_border_condition_type()==get_down_border_condition_type_halfspace()) then
            f=1d4
        endif
    endfunction get_full_h
    
    integer(4) function get_number_of_layer(z) result(f)
    use system,only:print_error
    implicit none
        real(8),intent(in)::z
        
        integer(4) n
        real(8) zr
        
        integer(4) i
        
        if(z>0d0) call print_error("main_parameters.get_number_of_layer","z>0d0")
        
        zr=z
        do i=1,number_of_layers
            zr=zr+h(i)
            
            if(zr>=0d0) then
                f=i
                return
            endif
        enddo
        
        if(get_down_border_condition_type()==get_down_border_condition_type_halfspace()) then
            f=number_of_layers
        else
            call print_error("main_parameters.get_number_of_layer","z<-\sum\limits_{n=1}^{number_of_layer}h(n)")
        endif
    endfunction get_number_of_layer
    
    complex(8) function Cijkl(i,j,k,l,number_of_layer) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::i
        integer(4),intent(in)::j
        integer(4),intent(in)::k
        integer(4),intent(in)::l
        integer(4),intent(in)::number_of_layer
        
        integer(4) alpha,beta
        
        if(.not.anisotropic==1) call print_error("main_parameters.Cijkl","not anisotropic materials")
        
        if(i<1.or.i>3) call print_error("main_parameters.Cijkl","i<1.or.i>3")
        if(j<1.or.j>3) call print_error("main_parameters.Cijkl","j<1.or.j>3")
        if(k<1.or.k>3) call print_error("main_parameters.Cijkl","k<1.or.k>3")
        if(l<1.or.l>3) call print_error("main_parameters.Cijkl","l<1.or.l>3")
        if(number_of_layer<1.and.number_of_layer>number_of_layers) then
            call print_error("main_parameters.Cijkl","number_of_layer<1.and.number_of_layer>number_of_layers")
        endif
        
        if(i==1.and.j==1) then
            alpha=1
        elseif(i==2.and.j==2) then
            alpha=2
        elseif(i==3.and.j==3) then
            alpha=3
        elseif(i==2.and.j==3.or.i==3.and.j==2) then
            alpha=4
        elseif(i==1.and.j==3.or.i==3.and.j==1) then
            alpha=5
        elseif(i==1.and.j==2.or.i==2.and.j==1) then
            alpha=6
        endif
        
        if(k==1.and.l==1) then
            beta=1
        elseif(k==2.and.l==2) then
            beta=2
        elseif(k==3.and.l==3) then
            beta=3
        elseif(k==2.and.l==3.or.k==3.and.l==2) then
            beta=4
        elseif(k==1.and.l==3.or.k==3.and.l==1) then
            beta=5
        elseif(k==1.and.l==2.or.k==2.and.l==1) then
            beta=6
        endif
        
        f=Calphabeta(alpha,beta,number_of_layer)
    endfunction Cijkl
    
    subroutine count_isotropic_material_parameters_form_lambda_mu(lambda,mu,E,nu,Cp,Cs,rho)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::lambda
        complex(8),intent(in)::mu
        complex(8),intent(out)::E
        complex(8),intent(out)::nu
        complex(8),intent(out)::Cp
        complex(8),intent(out)::Cs
        real(8),intent(in)::rho
        
        if(real(lambda)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_lambda_mu","real(lambda)<epsilon")
        if(real(mu)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_lambda_mu","real(mu)<epsilon")
        if(rho<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_lambda_mu","rho<epsilon")
        
        E=mu*(3d0*lambda+2d0*mu)/(lambda+mu)
        nu=lambda*0.5d0/(lambda+mu)
        
        Cp=sqrt((lambda+mu+mu)/rho)
        Cs=sqrt(mu/rho)
    endsubroutine count_isotropic_material_parameters_form_lambda_mu
    subroutine count_isotropic_material_parameters_form_E_nu(lambda,mu,E,nu,Cp,Cs,rho)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(out)::lambda
        complex(8),intent(out)::mu
        complex(8),intent(in)::E
        complex(8),intent(in)::nu
        complex(8),intent(out)::Cp
        complex(8),intent(out)::Cs
        real(8),intent(in)::rho
        
        if(real(E)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_E_nu","real(E)<epsilon")
        if(.not.(0d0<=real(nu).and.real(nu)<=0.5d0)) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_E_nu",".not.(0d0<=real(nu).and.real(nu)<=0.5d0)")
        if(rho<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_E_nu","rho<epsilon")
        
        lambda=nu*E/(1d0+nu)/(1d0-2d0*nu)
        mu=E*0.5d0/(1d0+nu)
        
        Cp=sqrt(E/rho)*sqrt((1d0-nu)/(1d0+nu)/(1d0-2d0*nu))
        Cs=sqrt(E/rho)*sqrt(0.5d0/(1d0+nu))
    endsubroutine count_isotropic_material_parameters_form_E_nu
    subroutine count_isotropic_material_parameters_form_Cp_Cs(lambda,mu,E,nu,Cp,Cs,rho)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(out)::lambda
        complex(8),intent(out)::mu
        complex(8),intent(out)::E
        complex(8),intent(out)::nu
        complex(8),intent(in)::Cp
        complex(8),intent(in)::Cs
        real(8),intent(in)::rho
        
        if(real(Cp)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_Cp_Cs","real(Cp)<epsilon")
        if(real(Cs)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_Cp_Cs","real(Cs)<epsilon")
        if(real(Cp)<real(Cs)) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_Cp_Cs","real(Cp)<real(Cs)")
        if(rho<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_Cp_Cs","rho<epsilon")
        
        lambda=rho*Cp*Cp-2d0*rho*Cs*Cs
        mu=rho*Cs*Cs
        
        E=Cs*Cs*rho*2d0*(1d0+(2d0*Cs*Cs-Cp*Cp)/(2d0*Cs*Cs-2d0*Cp*Cp))
        nu=(2d0*Cs*Cs-Cp*Cp)/(2d0*Cs*Cs-2d0*Cp*Cp)
    endsubroutine count_isotropic_material_parameters_form_Cp_Cs
    
    integer(4) function get_anisotropic() result(f)
    use system,only:print_error,print_warning
    implicit none
        if(anisotropic==-1) call print_error("main_parameters.get_anisotropic","anisotropic type is not defenited")
        f=anisotropic
    endfunction get_anisotropic
    subroutine set_anisotropic(anisotropic_v)
    use system,only:print_error,print_warning
    implicit none
        integer(4),intent(in)::anisotropic_v
        
        !if(.not.anisotropic==-1) call print_error("main_parameters.set_anisotropic","anisotropic type is already defenited")
        if(.not.anisotropic_v==0.and..not.anisotropic_v==1) call print_error("main_parameters.set_anisotropic",&
            ".not.anisotropic==0.and.not.anisotropic==1")
        
        anisotropic=anisotropic_v
    endsubroutine set_anisotropic
    
    integer(4) function get_down_border_condition_type_fixed_border() result(f)
    implicit none
        f=1
    endfunction get_down_border_condition_type_fixed_border
    integer(4) function get_down_border_condition_type_free_border() result(f)
    implicit none
        f=2
    endfunction get_down_border_condition_type_free_border
    integer(4) function get_down_border_condition_type_halfspace() result(f)
    implicit none
        f=3
    endfunction get_down_border_condition_type_halfspace
    integer(4) function get_down_border_condition_type() result(f)
    use system,only:print_error,print_warning
    implicit none
        if(down_border_condition_type==-1) &
            call print_error("main_parameters.get_down_border_condition_type","down_border_condition_type status is not defenited")
        f=down_border_condition_type
    endfunction get_down_border_condition_type
    subroutine set_down_border_condition_type(down_border_condition_type_v)
    use system,only:print_error,print_warning
    implicit none
        integer(4),intent(in)::down_border_condition_type_v
        
        !if(.not.down_border_condition_type==-1) &
        !    call print_error("main_parameters.set_down_border_condition_type","down_border_condition_type status is already defenited")
        if(.not.down_border_condition_type_v==get_down_border_condition_type_fixed_border()&
            .and..not.down_border_condition_type_v==get_down_border_condition_type_free_border()&
            .and..not.down_border_condition_type_v==get_down_border_condition_type_halfspace()) &
            call print_error("main_parameters.set_down_border_condition_type",&
            "down_border_condition_type is incorrect")
        
        down_border_condition_type=down_border_condition_type_v
    endsubroutine set_down_border_condition_type
endmodule main_parameters_
