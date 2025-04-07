module pre_and_post_runtime_actions
implicit none
    public::init,destructor
    
    private::establish_main_parameters
    contains
    
    subroutine init()
    use main_parameters,only:init_main_parameters
    use sigma_and_eigenvectors,only:init_sigma_and_eigenvectors
    use count_K,only:init_count_K
    use integral_solution,only:init_integral_solution
    use integral_solution_isotropic,only:init_integral_solution_isotropic
    use load_experimental_measurements,only:init_load_experimental_measurements
    use matrix_pencil_method,only:init_matrix_pencil_method
    use matrix_pencil_method_basis,only:init_matrix_pencil_method_basis
    implicit none
        call init_main_parameters(establish_main_parameters)
        call init_sigma_and_eigenvectors
        call init_count_K
        call init_integral_solution
        call init_integral_solution_isotropic
        call init_load_experimental_measurements
        call init_matrix_pencil_method
        call init_matrix_pencil_method_basis
    endsubroutine init
    subroutine destructor
    use main_parameters,only:destructor_main_parameters
    use sigma_and_eigenvectors,only:destructor_sigma_and_eigenvectors
    use count_K,only:destructor_count_K
    use integral_solution,only:destructor_integral_solution
    use integral_solution_isotropic,only:destructor_integral_solution_isotropic
    use load_experimental_measurements,only:destructor_load_experimental_measurements
    use matrix_pencil_method,only:destructor_matrix_pencil_method
    use matrix_pencil_method_basis,only:destructor_matrix_pencil_method_basis
    implicit none
        call destructor_main_parameters
        call destructor_sigma_and_eigenvectors
        call destructor_count_K
        call destructor_integral_solution
        call destructor_integral_solution_isotropic
        call destructor_load_experimental_measurements
        call destructor_matrix_pencil_method
        call destructor_matrix_pencil_method_basis
    endsubroutine destructor
    
    subroutine establish_main_parameters()
    use main_parameters,only:&
        set_anisotropic_type,set_anisotropic_type,&
        material_isotropic_type,material_anisotropic_type,&
        
        set_down_border_condition_type,&
        down_border_condition_type_fixed_border,&
        down_border_condition_type_free_border,&
        down_border_condition_type_halfspace,&
        
        set_omega,set_f,&
        set_Q,set_qomega,&
        
        set_number_of_layers,&
        
        set_layer_rho,set_layer_h,&
        
        set_layer_E_nu,set_layer_lambda_mu,set_layer_Cp_Cs,&
        set_layer_Calphabeta,&
        
        check_correct_completing_parameters_establishment_throwable,&
        get_layer_parameter
    use math,only:c0
    use system,only:end_program_pause
    implicit none
        call set_anisotropic_type(material_isotropic_type())
        
        call set_down_border_condition_type(down_border_condition_type_free_border())
        
        call set_omega(3d0)
        call set_Q(Q)
        call set_Qomega(Qomega)
        
        call set_number_of_layers(1)
        
        call set_layer_h(1,1d0)
        call set_layer_rho(1,1d0)
        
        call set_layer_E_nu(1,1d0+c0,0.2d0)
        
        if(.not.check_correct_completing_parameters_establishment_throwable()) call end_program_pause()
        
        print*,"h",get_layer_parameter(layer=1,parameter_name_length=1,parameter_name="h")
        print*,"rho",get_layer_parameter(layer=1,parameter_name_length=3,parameter_name="rho")
        print*,"E",get_layer_parameter(layer=1,parameter_name_length=1,parameter_name="E")
        print*,"nu",get_layer_parameter(layer=1,parameter_name_length=2,parameter_name="nu")
        print*,"lambda",get_layer_parameter(layer=1,parameter_name_length=6,parameter_name="lambda")
        print*,"mu",get_layer_parameter(layer=1,parameter_name_length=2,parameter_name="mu")
        print*,"Cp",get_layer_parameter(layer=1,parameter_name_length=2,parameter_name="Cp")
        print*,"Cs",get_layer_parameter(layer=1,parameter_name_length=2,parameter_name="Cs")
        
    contains
        complex(8) function Q(ind,alpha,beta) result(f)
        use math,only:ci
        use system,only:print_error
        implicit none
            integer(4),intent(in)::ind
            complex(8),intent(in)::alpha
            complex(8),intent(in)::beta
            
            if(.not.(1.le.ind.and.ind.le.3)) &
                call print_error("pre_and_post_runtime_actions.establish_main_parameters",&
                ".not.(1.le.ind.and.ind.le.3)")
            
            if(ind==1) then
                f=0d0
            elseif(ind==2) then
                f=0d0
            else
                f=1d0
            endif
        endfunction Q
    endsubroutine establish_main_parameters
    pure complex(8) function Qomega(omega) result(f)
    use math,only:ci,c0
    implicit none
        complex(8),intent(in)::omega
        
        f=1d0
    endfunction Qomega
endmodule pre_and_post_runtime_actions
