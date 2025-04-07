module main_parameters
implicit none
    public::get_anisotropic_type,set_anisotropic_type,&
        material_isotropic_type,material_anisotropic_type
    public::set_down_border_condition_type,get_down_border_condition_type,&
        down_border_condition_type_fixed_border,down_border_condition_type_free_border,&
        down_border_condition_type_halfspace
    
    public::set_omega,get_omega,set_f,get_f
    
    public::set_Q,get_Q,set_Qomega,get_Qomega
    
    public::get_number_of_layers,set_number_of_layers
    
    public::get_layer_parameter
    
    public::get_layer_E,get_layer_nu,get_layer_lambda,get_layer_mu,get_layer_Cp,get_layer_Cs
    public::get_layer_rho,get_layer_h,get_full_h,get_number_of_layer_from_z
    public::set_layer_E_nu,set_layer_lambda_mu,set_layer_Cp_Cs
    public::set_layer_rho,set_layer_h
    
    public::get_layer_Cijkl_parameter,get_layer_Calphabeta_parameter
    public::set_layer_Calphabeta
    
    public::count_isotropic_material_parameters_form_lambda_mu,&
        count_isotropic_material_parameters_form_E_nu,&
        count_isotropic_material_parameters_form_Cp_Cs
    public::count_isotropic_material_parameter_E_form_lambda_mu,count_isotropic_material_parameter_nu_form_lambda_mu,&
        count_isotropic_material_parameter_Cp_form_lambda_mu,count_isotropic_material_parameter_Cs_form_lambda_mu,&
        count_isotropic_material_parameter_lambda_form_E_nu,count_isotropic_material_parameter_mu_form_E_nu,&
        count_isotropic_material_parameter_Cp_form_E_nu,count_isotropic_material_parameter_Cs_form_E_nu,&
        count_isotropic_material_parameter_lambda_form_Cp_Cs,count_isotropic_material_parameter_mu_form_Cp_Cs,&
        count_isotropic_material_parameter_E_form_Cp_Cs,count_isotropic_material_parameter_nu_form_Cp_Cs
    
    
    
    integer(4),private::anisotropic_type=-1
    integer(4),private::down_border_condition_type=-1
    
    real(8),private::omega
    
    abstract interface
        function Q_type(ind,alpha,beta) result(f)
            complex(8)::f
            integer(4),intent(in)::ind!1,2,3
            complex(8),intent(in)::alpha
            complex(8),intent(in)::beta
        endfunction Q_type
        function Qomega_type(omega) result(f)
            complex(8)::f
            complex(8),intent(in)::omega
        endfunction Qomega_type
    endinterface
    procedure(Q_type),pointer,private::Q
    procedure(Qomega_type),pointer,private::Qomega
    logical(1),private::Q_setted=.false.,Qomega_setted=.false.
    
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
    
    
    
    real(8),private::full_h=-1d0
    
    private::free_matrixes,create_matrixes
    contains
    
    integer(4) function get_anisotropic_type() result(f)
    use system,only:print_error
    implicit none
        if(anisotropic_type==-1) &
            call print_error("main_parameters.get_anisotropic_type","anisotripoc_type is not setted")
        
        f=anisotropic_type
    endfunction get_anisotropic_type
    subroutine set_anisotropic_type(anisotropic_type_)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::anisotropic_type_
        
        if(.not.anisotropic_type==-1) &
            call print_error("main_parameters.set_anisotropic_type","anisotripoc_type has been already setted")
        if(.not.(anisotropic_type_==material_isotropic_type().or.anisotropic_type_==material_anisotropic_type())) &
            call print_error("main_parameters.set_anisotropic_type","incorrect anisotripoc_type")
        
        anisotropic_type=anisotropic_type_
    endsubroutine set_anisotropic_type
    integer(4) function material_isotropic_type() result(f)
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
            down_border_condition_type_==down_border_condition_type_halfspace()&
            )) call print_error("main_parameters.set_down_border_condition_type","incorrect down_border_condition_type_")
        
        down_border_condition_type=down_border_condition_type_
    endsubroutine set_down_border_condition_type
    integer(4) function get_down_border_condition_type() result(f)
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
    
    
    
    subroutine set_omega(omega_)
    implicit none
        real(8),intent(in)::omega_
        
        omega=omega_
    endsubroutine set_omega
    subroutine set_f(f_)
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
        procedure(Q_type),pointer,intent(in)::Q_
        
        Q_setted=.true.
        Q=>Q_
    endsubroutine set_Q
    complex(8) function get_Q(ind,alpha,beta) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::ind
        complex(8),intent(in)::alpha
        complex(8),intent(in)::beta
        
        if(.not.Q_setted) call print_error("main_parameters.set_Q","Q is not setted")
        if(.not.(1<=ind.and.ind<=3)) call print_error("main_parameters.set_Q","incorrect ind")
        
        f=Q(ind,alpha,beta)
    endfunction get_Q
    subroutine set_Qomega(Qomega_)
    implicit none
        procedure(Qomega_type),pointer,intent(in)::Qomega_
        
        Qomega_setted=.true.
        Qomega=>Qomega_
    endsubroutine set_Qomega
    complex(8) function get_Qomega(omega) result(f)
    use system,only:print_error
    implicit none
        complex(8),intent(in)::omega
        
        if(.not.Qomega_setted) call print_error("main_parameters.set_Qomega","Qomega is not setted")
        
        f=Qomega(omega)
    endfunction get_Qomega
    
    
    
    integer(4) function get_number_of_layers() result(f)
    use system,only:print_error
    implicit none
        if(number_of_layers==-1) call print_error("main_parameters.get_number_of_layers","number_of_layer is not setted")
        
        f=number_of_layers
    endfunction get_number_of_layers
    subroutine set_number_of_layers(number_of_layers_)
    use system,only:print_error,print_warning
    implicit none
        integer(4),intent(in)::number_of_layers_
        
        if(number_of_layers_<1) call print_error("mainr_parameters.set_number_of_layers","number_of_layers_<1")
        if(number_of_layers_>10) call print_warning("mainr_parameters.set_number_of_layers","number_of_layers_>10")
        
        number_of_layers=number_of_layers_
        call free_matrixes()
        call create_matrixes()
        
        full_h=-1d0
    endsubroutine set_number_of_layers
    
    subroutine free_matrixes()
    use system,only:print_error
    implicit none
        if(anisotropic_type==-1) call print_error("main_parameters.free_matrixes","anisotropic_type is not setted")
        if(.not.matrixes_created) call print_error("main_parameters.free_matrixes","matrixes is not created")
        
        matrixes_created=.false.
        
        deallocate(h)
        deallocate(rho)
        
        if(anisotropic_type==material_isotropic_type()) then
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
        
        if(anisotropic_type==material_isotropic_type()) then
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
    
    
    
    
    complex(8) function get_layer_E(layer) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::layer
        
        if(.not.anisotropic_type==material_isotropic_type()) &
            call print_error("main_parameters.get_layer_E",".not.anisotropic_type==material_isotropic_type()")
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.get_layer_E",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        if(E(layer)==-1d0) call print_error("main_parameters.get_layer_E","E(layer) is not setted")
        
        f=E(layer)
    endfunction get_Layer_E
    real(8) function get_layer_nu(layer) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::layer
        
        if(.not.anisotropic_type==material_isotropic_type()) &
            call print_error("main_parameters.get_layer_nu",".not.anisotropic_type==material_isotropic_type()")
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.get_layer_nu",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        if(nu(layer)==-1d0) call print_error("main_parameters.get_layer_nu","nu(layer) is not setted")
        
        f=nu(layer)
    endfunction get_Layer_nu
    
    complex(8) function get_layer_lambda(layer) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::layer
        
        if(.not.anisotropic_type==material_isotropic_type()) &
            call print_error("main_parameters.get_layer_lambda",".not.anisotropic_type==material_isotropic_type()")
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.get_layer_lambda",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        if(lambda(layer)==-1d0) call print_error("main_parameters.get_layer_lambda","lambda(layer) is not setted")
        
        f=lambda(layer)
    endfunction get_Layer_lambda
    complex(8) function get_layer_mu(layer) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::layer
        
        if(.not.anisotropic_type==material_isotropic_type()) &
            call print_error("main_parameters.get_layer_mu",".not.anisotropic_type==material_isotropic_type()")
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.get_layer_mu",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        if(mu(layer)==-1d0) call print_error("main_parameters.get_layer_mu","mu(layer) is not setted")
        
        f=mu(layer)
    endfunction get_Layer_mu
    
    complex(8) function get_layer_Cp(layer) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::layer
        
        if(.not.anisotropic_type==material_isotropic_type()) &
            call print_error("main_parameters.get_layer_Cp",".not.anisotropic_type==material_isotropic_type()")
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.get_layer_Cp",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        if(Cp(layer)==-1d0) call print_error("main_parameters.get_layer_Cp","Cp(layer) is not setted")
        
        f=Cp(layer)
    endfunction get_Layer_Cp
    complex(8) function get_layer_Cs(layer) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::layer
        
        if(.not.anisotropic_type==material_isotropic_type()) &
            call print_error("main_parameters.get_layer_Cs",".not.anisotropic_type==material_isotropic_type()")
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.get_layer_Cs",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        if(Cs(layer)==-1d0) call print_error("main_parameters.get_layer_Cs","Cs(layer) is not setted")
        
        f=Cs(layer)
    endfunction get_Layer_Cs
    
    
    
    real(8) function get_layer_rho(layer) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::layer
        
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.get_layer_rho",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        if(rho(layer)==-1d0) call print_error("main_parameters.get_layer_rho","rho(layer) is not setted")
        
        f=rho(layer)
    endfunction get_layer_rho
    real(8) function get_layer_h(layer) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::layer
        
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.get_layer_h",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        if(h(layer)==-1d0) call print_error("main_parameters.get_layer_h","h(layer) is not setted")
        
        f=h(layer)
    endfunction get_layer_h
    real(8) function get_full_h() result(f)
    use system,only:print_error
    implicit none
        integer(4) i
        
        if(full_h==-1d0) then
            f=full_h
            return
        endif
        
        full_h=0d0
        do i=1,number_of_layers
            if(h(i)==-1d0) call print_error("main_parameters.get_full_h","h(i) is not setted")
            
            full_h=full_h+h(i)
        enddo
    endfunction get_full_h
    integer(4) function get_number_of_layer_from_z(z) result(f)
    use system,only:print_error
    implicit none
        real(8),intent(in)::z
        
        integer(4) n
        real(8) zr
        
        integer(4) i
        
        if(z>0d0) call print_error("main_parameters.get_number_of_layer_from_z","z>0d0")
        
        zr=z
        do i=1,number_of_layers
            zr=zr+h(i)
            
            if(zr>=0d0) then
                f=i
                return
            endif
        enddo
        
        if(get_down_border_condition_type()==down_border_condition_type_halfspace()) then
            f=number_of_layers
        else
            call print_error("main_parameters.get_number_of_layer_from_z","z<-\sum\limits_{n=1}^{number_of_layer}h(n)")
        endif
    endfunction get_number_of_layer_from_z
    
    
    
    subroutine set_layer_E_nu(layer,E_,nu_)
    use system,only:print_error
    use math,only:epsilon
    implicit none
        integer(4),intent(in)::layer
        complex(8),intent(in)::E_
        real(8),intent(in)::nu_
        
        if(.not.anisotropic_type==material_isotropic_type()) &
            call print_error("main_parameters.set_layer_E_nu",".not.anisotropic_type==material_isotropic_type()")
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.set_layer_E_nu",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        if(rho(layer)==-1d0) call print_error("main_parameters.set_layer_E_nu","rho(layer) is not setted")
        if(.not.(epsilon<real(E_).and.real(E_)<1000d0)) &
            call printerror("main_parameters.set_layer_E_nu",".not.(epsilon<real(E_).and.real(E_)<1000d0)")
        if(.not.(0d0<nu_.and.nu_<0.5d0)) call printerror("main_parameters.set_layer_E_nu",".not.(0d0<nu_.and.nu_<0.5d0)")
        
        E(layer)=E_
        nu(layer)=nu_
        call count_isotropic_material_parameters_form_E_nu(lambda(layer),mu(layer),E(layer),nu(layer),Cp(layer),Cs(layer),rho(layer))
    endsubroutine set_layer_E_nu
    
    subroutine set_layer_lambda_mu(layer,lambda_,mu_)
    use system,only:print_error
    use math,only:epsilon
    implicit none
        integer(4),intent(in)::layer
        complex(8),intent(in)::lambda_
        real(8),intent(in)::mu_
        
        if(.not.anisotropic_type==material_isotropic_type()) &
            call print_error("main_parameters.set_layer_lambda_mu",".not.anisotropic_type==material_isotropic_type()")
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.set_layer_lambda_mu",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        if(rho(layer)==-1d0) call print_error("main_parameters.set_layer_lambda_mu","rho(layer) is not setted")
        if(.not.(epsilon<real(lambda_).and.real(lambda_)<1000d0)) &
            call printerror("main_parameters.set_layer_lambda_mu",".not.(epsilon<real(lambda_).and.real(lambda_)<1000d0)")
        if(.not.(epsilon<real(mu_).and.real(mu_)<1000d0)) &
            call printerror("main_parameters.set_layer_lambda_mu",".not.(epsilon<real(mu_).and.real(mu_)<1000d0)")
        
        lambda(layer)=lambda_
        mu(layer)=mu_
        call count_isotropic_material_parameters_form_lambda_mu(lambda(layer),mu(layer),E(layer),nu(layer),Cp(layer),Cs(layer),rho(layer))
    endsubroutine set_layer_lambda_mu
    
    subroutine set_layer_Cp_Cs(layer,Cp_,Cs_)
    use system,only:print_error
    use math,only:epsilon
    implicit none
        integer(4),intent(in)::layer
        complex(8),intent(in)::Cp_
        real(8),intent(in)::Cs_
        
        if(.not.anisotropic_type==material_isotropic_type()) &
            call print_error("main_parameters.set_layer_Cp_Cs",".not.anisotropic_type==material_isotropic_type()")
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.set_layer_Cp_Cs",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        if(rho(layer)==-1d0) call print_error("main_parameters.set_layer_Cp_Cs","rho(layer) is not setted")
        if(.not.(epsilon<real(Cp_).and.real(Cp_)<1000d0)) &
            call printerror("main_parameters.set_layer_Cp_Cs",".not.(epsilon<real(Cp_).and.real(Cp_)<1000d0)")
        if(.not.(epsilon<real(Cs_).and.real(Cs_)<1000d0)) &
            call printerror("main_parameters.set_layer_Cp_Cs",".not.(epsilon<real(Cs_).and.real(Cs_)<1000d0)")
        if(real(cp_)<real(Cs_)) call print_error("main_parameters.set_layer_Cp_Cs","real(cp_)<real(Cs_)")
        
        Cp(layer)=Cp_
        Cs(layer)=Cs_
        call count_isotropic_material_parameters_form_lambda_mu(lambda(layer),mu(layer),E(layer),nu(layer),Cp(layer),Cs(layer),rho(layer))
    endsubroutine set_layer_Cp_Cs
    
    
    
    subroutine set_layer_h(layer,h_)
    use system,only:print_error,print_warning
    use math,only:epsilon
    implicit none
        integer(4),intent(in)::layer
        real(8),intent(in)::h_
        
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.set_layer_h",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        if(.not.(epsilon.le.h_.and.h_.le.100d0)) &
            call print_error("main_parameters.set_layer_h",".not.(epsilon.le.h_.and.h_.le.100d0)")
        
        full_h=-1d0
        h(layer)=h_
    endsubroutine set_layer_h
    subroutine set_layer_rho(layer,rho_)
    use system,only:print_error,print_warning
    use math,only:epsilon
    implicit none
        integer(4),intent(in)::layer
        real(8),intent(in)::rho_
        
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.set_layer_rho",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        if(.not.(epsilon.le.rho_.and.rho_.le.100d0)) &
            call print_error("main_parameters.set_layer_rho",".not.(epsilon.le.rho_.and.rho_.le.100d0)")
        
        rho(layer)=rho_
    endsubroutine set_layer_rho
    
    
    
    complex(8) function get_layer_Calphabeta_parameter(layer,alpha,beta) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::layer
        integer(4),intent(in)::alpha
        integer(4),intent(in)::beta
        
        if(.not.anisotropic_type==material_anisotropic_type()) &
            call print_error("main_parameters.get_layer_Calphabeta_parameter",".not.anisotropic_type==material_anisotropic_type()")
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.get_layer_Calphabeta_parameter",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        if(.not.(1.le.alpha.and.alpha.le.6)) &
            call print_error("main_parameters.get_layer_Calphabeta_parameter",".not.(1.le.alpha.and.alpha.le.6)")
        if(.not.(1.le.beta.and.beta.le.6)) &
            call print_error("main_parameters.get_layer_Calphabeta_parameter",".not.(1.le.beta.and.beta.le.6)")
        
        f=Calphabeta(layer,alpha,beta)
        if(f==-1d0) call print_error("main_parameters.get_layer_Calphabeta_parameter","Calphabeta(layer,alpha,beta)==-1d0")
    endfunction get_layer_Calphabeta_parameter
    complex(8) function get_layer_Cijkl_parameter(layer,i,j,k,l) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::layer
        integer(4),intent(in)::i
        integer(4),intent(in)::j
        integer(4),intent(in)::k
        integer(4),intent(in)::l
        
        integer(4) alpha,beta
        
        if(.not.anisotropic_type==material_anisotropic_type()) &
            call print_error("main_parameters.get_layer_Cijkl_parameter",".not.anisotropic_type==material_anisotropic_type()")
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.get_layer_Cijkl_parameter",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        if(.not.(1.le.i.and.i.le.3)) &
            call print_error("main_parameters.get_layer_Cijkl_parameter",".not.(1.le.i.and.i.le.3)")
        if(.not.(1.le.j.and.j.le.3)) &
            call print_error("main_parameters.get_layer_Cijkl_parameter",".not.(1.le.j.and.j.le.3)")
        if(.not.(1.le.k.and.k.le.3)) &
            call print_error("main_parameters.get_layer_Cijkl_parameter",".not.(1.le.k.and.k.le.3)")
        if(.not.(1.le.l.and.l.le.3)) &
            call print_error("main_parameters.get_layer_Cijkl_parameter",".not.(1.le.l.and.l.le.3)")
        
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
        
        f=get_layer_Calphabeta_parameter(layer,alpha,beta)
    endfunction get_layer_Cijkl_parameter
    
    
    
    
    subroutine set_layer_Calphabeta(layer,Calphabeta_)
    use system,only:print_error
    use math,only:epsilon
    implicit none
        integer(4),intent(in)::layer
        complex(8),intent(in)::Calphabeta_(6,6)
        
        complex(8) v
        integer(4) i,j
        
        if(.not.anisotropic_type==material_anisotropic_type()) &
            call print_error("main_parameters.get_layer_Cijkl_parameter",".not.anisotropic_type==material_anisotropic_type()")
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.get_layer_Cijkl_parameter",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        
        do i=1,6
            do j=1,6
                v=Calphabeta_(i,j)
                
                if(.not.(epsilon.le.real(v).and.real(v).le.1000d0)) &
                    call print_error("main_parameters.set_layer_Calphabeta",".not.(epsilon.le.real(v).and.real(v).le.1000d0)")
                
                Calphabeta(layer,i,j)=v
            enddo
        enddo
    endsubroutine set_layer_Calphabeta
    
    
    
    subroutine count_isotropic_material_parameters_form_lambda_mu(lambda,mu,E,nu,Cp,Cs,rho)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::lambda
        complex(8),intent(in)::mu
        complex(8),intent(out)::E
        real(8),intent(out)::nu
        complex(8),intent(out)::Cp
        complex(8),intent(out)::Cs
        real(8),intent(in)::rho
        
        if(real(lambda)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_lambda_mu","real(lambda)<epsilon")
        if(real(mu)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_lambda_mu","real(mu)<epsilon")
        if(rho<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_lambda_mu","rho<epsilon")
        
        E=count_isotropic_material_parameter_E_form_lambda_mu(lambda,mu)
        nu=count_isotropic_material_parameter_nu_form_lambda_mu(lambda,mu)
        
        Cp=count_isotropic_material_parameter_Cp_form_lambda_mu(lambda,mu,rho)
        Cs=count_isotropic_material_parameter_Cs_form_lambda_mu(lambda,mu,rho)
    endsubroutine count_isotropic_material_parameters_form_lambda_mu
    subroutine count_isotropic_material_parameters_form_E_nu(lambda,mu,E,nu,Cp,Cs,rho)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(out)::lambda
        complex(8),intent(out)::mu
        complex(8),intent(in)::E
        real(8),intent(in)::nu
        complex(8),intent(out)::Cp
        complex(8),intent(out)::Cs
        real(8),intent(in)::rho
        
        if(real(E)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_E_nu","real(E)<epsilon")
        if(.not.(0d0<=real(nu).and.real(nu)<=0.5d0)) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_E_nu",".not.(0d0<=real(nu).and.real(nu)<=0.5d0)")
        if(rho<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameters_form_E_nu","rho<epsilon")
        
        lambda=count_isotropic_material_parameter_lambda_form_E_nu(lambda,nu)
        mu=count_isotropic_material_parameter_mu_form_E_nu(lambda,nu)
        
        Cp=count_isotropic_material_parameter_Cp_form_E_nu(lambda,nu,rho)
        Cs=count_isotropic_material_parameter_Cs_form_E_nu(lambda,nu,rho)
    endsubroutine count_isotropic_material_parameters_form_E_nu
    subroutine count_isotropic_material_parameters_form_Cp_Cs(lambda,mu,E,nu,Cp,Cs,rho)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(out)::lambda
        complex(8),intent(out)::mu
        complex(8),intent(out)::E
        real(8),intent(out)::nu
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
        
        lambda=count_isotropic_material_parameter_lambda_form_Cp_Cs(Cp,Cs,rho)
        mu=count_isotropic_material_parameter_mu_form_Cp_Cs(Cp,Cs,rho)
        
        E=count_isotropic_material_parameter_E_form_Cp_Cs(Cp,Cs,rho)
        nu=count_isotropic_material_parameter_nu_form_Cp_Cs(Cp,Cs,rho)
    endsubroutine count_isotropic_material_parameters_form_Cp_Cs
    
    
    
    complex(8) function count_isotropic_material_parameter_E_form_lambda_mu(lambda,mu) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::lambda
        complex(8),intent(in)::mu
        
        if(real(lambda)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_E_form_lambda_mu","real(lambda)<epsilon")
        if(real(mu)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_E_form_lambda_mu","real(mu)<epsilon")
        
        f=mu*(3d0*lambda+2d0*mu)/(lambda+mu)
    endfunction count_isotropic_material_parameter_E_form_lambda_mu
    real(8) function count_isotropic_material_parameter_nu_form_lambda_mu(lambda,mu) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::lambda
        complex(8),intent(in)::mu
        
        if(real(lambda)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_nu_form_lambda_mu","real(lambda)<epsilon")
        if(real(mu)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_nu_form_lambda_mu","real(mu)<epsilon")
        
        f=lambda*0.5d0/(lambda+mu)
    endfunction count_isotropic_material_parameter_nu_form_lambda_mu
    complex(8) function count_isotropic_material_parameter_Cp_form_lambda_mu(lambda,mu,rho) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::lambda
        complex(8),intent(in)::mu
        real(8),intent(in)::rho
        
        if(real(lambda)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_Cp_form_lambda_mu","real(lambda)<epsilon")
        if(real(mu)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_Cp_form_lambda_mu","real(mu)<epsilon")
        if(rho<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_Cp_form_lambda_mu","rho<epsilon")
        
        f=sqrt((lambda+mu+mu)/rho)
    endfunction count_isotropic_material_parameter_Cp_form_lambda_mu
    complex(8) function count_isotropic_material_parameter_Cs_form_lambda_mu(lambda,mu,rho) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::lambda
        complex(8),intent(in)::mu
        real(8),intent(in)::rho
        
        if(real(lambda)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_Cs_form_lambda_mu","real(lambda)<epsilon")
        if(real(mu)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_Cs_form_lambda_mu","real(mu)<epsilon")
        if(rho<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_Cs_form_lambda_mu","rho<epsilon")
        
        f=sqrt(mu/rho)
    endfunction count_isotropic_material_parameter_Cs_form_lambda_mu
    
    complex(8) function count_isotropic_material_parameter_lambda_form_E_nu(E,nu) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::E
        real(8),intent(in)::nu
        
        if(real(E)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_lambda_form_E_nu","real(E)<epsilon")
        if(real(nu)<epsilon.or.real(nu)>0.5d0-epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_lambda_form_E_nu","real(nu)<epsilon.or.real(nu)>0.5d0-epsilon")
        
        f=nu*E/(1d0+nu)/(1d0-2d0*nu)
    endfunction count_isotropic_material_parameter_lambda_form_E_nu
    complex(8) function count_isotropic_material_parameter_mu_form_E_nu(E,nu) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::E
        real(8),intent(in)::nu
        
        if(real(E)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_mu_form_E_nu","real(E)<epsilon")
        if(real(nu)<epsilon.or.real(nu)>0.5d0-epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_mu_form_E_nu","real(nu)<epsilon.or.real(nu)>0.5d0-epsilon")
        
        f=E*0.5d0/(1d0+nu)
    endfunction count_isotropic_material_parameter_mu_form_E_nu
    complex(8) function count_isotropic_material_parameter_Cp_form_E_nu(E,nu,rho) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::E
        real(8),intent(in)::nu
        real(8),intent(in)::rho
        
        if(real(E)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_Cp_form_E_nu","real(E)<epsilon")
        if(real(nu)<epsilon.or.real(nu)>0.5d0-epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_Cp_form_E_nu","real(nu)<epsilon.or.real(nu)>0.5d0-epsilon")
        if(rho<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_Cp_form_E_nu","rho<epsilon")
        
        f=sqrt(E/rho)*sqrt((1d0-nu)/(1d0+nu)/(1d0-2d0*nu))
    endfunction count_isotropic_material_parameter_Cp_form_E_nu
    complex(8) function count_isotropic_material_parameter_Cs_form_E_nu(E,nu,rho) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::E
        real(8),intent(in)::nu
        real(8),intent(in)::rho
        
        if(real(E)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_Cs_form_E_nu","real(E)<epsilon")
        if(real(nu)<epsilon.or.real(nu)>0.5d0-epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_Cs_form_E_nu","real(nu)<epsilon.or.real(nu)>0.5d0-epsilon")
        if(rho<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_Cs_form_E_nu","rho<epsilon")
        
        f=sqrt(E/rho)*sqrt(0.5d0/(1d0+nu))
    endfunction count_isotropic_material_parameter_Cs_form_E_nu
    
    complex(8) function count_isotropic_material_parameter_lambda_form_Cp_Cs(Cp,Cs,rho) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::Cp
        complex(8),intent(in)::Cs
        real(8),intent(in)::rho
        
        if(real(Cp)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_lambda_form_Cp_Cs","real(Cp)<epsilon")
        if(real(Cs)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_lambda_form_Cp_Cs","real(Cs)<epsilon")
        if(real(Cp)<real(Cs)) &
            call print_error("main_parameters.count_isotropic_material_parameter_lambda_form_Cp_Cs","real(Cp)<real(Cs)")
        if(rho<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_lambda_form_Cp_Cs","rho<epsilon")
        
        f=rho*Cp*Cp-2d0*rho*Cs*Cs
    endfunction count_isotropic_material_parameter_lambda_form_Cp_Cs
    complex(8) function count_isotropic_material_parameter_mu_form_Cp_Cs(Cp,Cs,rho) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::Cp
        complex(8),intent(in)::Cs
        real(8),intent(in)::rho
        
        if(real(Cp)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_mu_form_Cp_Cs","real(Cp)<epsilon")
        if(real(Cs)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_mu_form_Cp_Cs","real(Cs)<epsilon")
        if(real(Cp)<real(Cs)) &
            call print_error("main_parameters.count_isotropic_material_parameter_mu_form_Cp_Cs","real(Cp)<real(Cs)")
        if(rho<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_mu_form_Cp_Cs","rho<epsilon")
        
        f=rho*Cs*Cs
    endfunction count_isotropic_material_parameter_mu_form_Cp_Cs
    complex(8) function count_isotropic_material_parameter_E_form_Cp_Cs(Cp,Cs,rho) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::Cp
        complex(8),intent(in)::Cs
        real(8),intent(in)::rho
        
        if(real(Cp)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_E_form_Cp_Cs","real(Cp)<epsilon")
        if(real(Cs)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_E_form_Cp_Cs","real(Cs)<epsilon")
        if(real(Cp)<real(Cs)) &
            call print_error("main_parameters.count_isotropic_material_parameter_E_form_Cp_Cs","real(Cp)<real(Cs)")
        if(rho<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_E_form_Cp_Cs","rho<epsilon")
        
        f=Cs*Cs*rho*2d0*(1d0+(2d0*Cs*Cs-Cp*Cp)/(2d0*Cs*Cs-2d0*Cp*Cp))
    endfunction count_isotropic_material_parameter_E_form_Cp_Cs
    real(8) function count_isotropic_material_parameter_nu_form_Cp_Cs(Cp,Cs,rho) result(f)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::Cp
        complex(8),intent(in)::Cs
        real(8),intent(in)::rho
        
        if(real(Cp)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_nu_form_Cp_Cs","real(Cp)<epsilon")
        if(real(Cs)<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_nu_form_Cp_Cs","real(Cs)<epsilon")
        if(real(Cp)<real(Cs)) &
            call print_error("main_parameters.count_isotropic_material_parameter_nu_form_Cp_Cs","real(Cp)<real(Cs)")
        if(rho<epsilon) &
            call print_error("main_parameters.count_isotropic_material_parameter_nu_form_Cp_Cs","rho<epsilon")
        
        f=(2d0*Cs*Cs-Cp*Cp)/(2d0*Cs*Cs-2d0*Cp*Cp)
    endfunction count_isotropic_material_parameter_nu_form_Cp_Cs
    
    
    
    
    complex(8) function get_layer_parameter(layer,parameter_name) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::layer
        character(len=10),intent(in)::parameter_name
        
        integer(4) alpha,beta,i,j,k,l
        
        if(.not.(1.le.layer.and.layer.le.number_of_Layers)) &
            call print_error("main_parameters.get_layer_parameter",".not.(1.le.layer.and.layer.le.number_of_Layers)")
        
        if(parameter_name=="E") then
            if(.not.anisotropic_type==material_isotropic_type()) &
                call print_error("main_parameters.get_Layer_parameter",".not.anisotropic_type==material_isotropic_type()")
            if(E(layer)==-1d0) call print_error("main_parameters.get_layer_parameter","E==-1d0")
            
            f=get_layer_E(layer)
            return
        endif
        if(parameter_name=="nu") then
            if(.not.anisotropic_type==material_isotropic_type()) &
                call print_error("main_parameters.get_Layer_parameter",".not.anisotropic_type==material_isotropic_type()")
            if(nu(layer)==-1d0) call print_error("main_parameters.get_layer_parameter","nu==-1d0")
            
            f=get_layer_nu(layer)
            return
        endif
        if(parameter_name=="lambda") then
            if(.not.anisotropic_type==material_isotropic_type()) &
                call print_error("main_parameters.get_Layer_parameter",".not.anisotropic_type==material_isotropic_type()")
            if(lambda(layer)==-1d0) call print_error("main_parameters.get_layer_parameter","lambda==-1d0")
            
            f=get_layer_lambda(layer)
            return
        endif
        if(parameter_name=="mu") then
            if(.not.anisotropic_type==material_isotropic_type()) &
                call print_error("main_parameters.get_Layer_parameter",".not.anisotropic_type==material_isotropic_type()")
            if(mu(layer)==-1d0) call print_error("main_parameters.get_layer_parameter","mu==-1d0")
            
            f=get_layer_mu(layer)
            return
        endif
        if(parameter_name=="Cp") then
            if(.not.anisotropic_type==material_isotropic_type()) &
                call print_error("main_parameters.get_Layer_parameter",".not.anisotropic_type==material_isotropic_type()")
            if(Cp(layer)==-1d0) call print_error("main_parameters.get_layer_parameter","Cp==-1d0")
            
            f=get_layer_Cp(layer)
            return
        endif
        if(parameter_name=="Cs") then
            if(.not.anisotropic_type==material_isotropic_type()) &
                call print_error("main_parameters.get_Layer_parameter",".not.anisotropic_type==material_isotropic_type()")
            if(Cs(layer)==-1d0) call print_error("main_parameters.get_layer_parameter","Cs==-1d0")
            
            f=get_layer_Cs(layer)
            return
        endif
        
        if(parameter_name=="h") then
            if(h(layer)==-1d0) call print_error("main_parameters.get_layer_parameter","h==-1d0")
            
            f=get_layer_h(layer)
            return
        endif
        if(parameter_name=="rho") then
            if(rho(layer)==-1d0) call print_error("main_parameters.get_layer_parameter","rho==-1d0")
            
            f=get_layer_rho(layer)
            return
        endif
        
        if(len_trim(parameter_name)==3.and.parameter_name(1:1)=='C'.and.&
            '1'.le.parameter_name(2:2).and.parameter_name(2:2).le.'6'.and.&
            '1'.le.parameter_name(3:3).and.parameter_name(3:3).le.'6') then
            if(.not.anisotropic_type==material_anisotropic_type()) &
                call print_error("main_parameters.get_Layer_parameter",".not.anisotropic_type==material_anisotropic_type()")
            
            read(parameter_name(2:2),'(i)'),alpha
            read(parameter_name(3:3),'(i)'),beta
            
            f=get_layer_Calphabeta_parameter(layer,alpha,beta)
            if(f==-1d0) call print_error("main_parameters.get_layer_parameter","Calphabeta(layer,i,j)==-1d0")
            return
        endif
        if(len_trim(parameter_name)==5.and.parameter_name(1:1)=='C'.and.&
            '1'.le.parameter_name(2:2).and.parameter_name(2:2).le.'3'.and.&
            '1'.le.parameter_name(3:3).and.parameter_name(3:3).le.'3'.and.&
            '1'.le.parameter_name(4:4).and.parameter_name(4:4).le.'3'.and.&
            '1'.le.parameter_name(5:5).and.parameter_name(5:5).le.'3') then
            if(.not.anisotropic_type==material_anisotropic_type()) &
                call print_error("main_parameters.get_Layer_parameter",".not.anisotropic_type==material_anisotropic_type()")
            
            read(parameter_name(2:2),'(i)'),i
            read(parameter_name(3:3),'(i)'),j
            read(parameter_name(4:4),'(i)'),k
            read(parameter_name(5:5),'(i)'),l
            
            f=get_layer_Cijkl_parameter(layer,i,j,k,l)
            if(f==-1d0) call print_error("main_parameters.get_layer_parameter","Cijkl(layer,i,j,k,l)==-1d0")
            return
        endif
        
        call print_error("main_parameters.get_layer_parameter","incorrect parameter name")
    endfunction get_Layer_parameter
endmodule main_parameters
