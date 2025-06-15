module material_properties_by_dispersion_curves
implicit none
    public::find_material_properties
    
    integer(4)::parameter_type_count=4
    
    private::parameter_type_count,get_parameter_type_number,get_parameter_type_name
    contains
    
    subroutine find_material_properties(dispersion_curves_size,dispersion_curves,all_disoersion_points_is_difenetly_true,&
        parameters_for_detect_size,parameters_layer,parameters_type,parameters_min,parameters_max,dparameters,&
        res_max_size,res,res_value_of_right,res_size)
    use main_parameters,only:get_number_of_layers,get_anisotropic_type,material_isotropic_type
    use gradient_descent_method,only:find_strict_minimum_points_on_surf
    use math,only:epsilon
    use system,only:print_error
    implicit none
        integer(4),intent(in)::dispersion_curves_size
        real(8),intent(in)::dispersion_curves(dispersion_curves_size,2)![(omega_1,value_1),...]
        logical(1),intent(in)::all_disoersion_points_is_difenetly_true!if true then a right value = max diff between dispersion curves themseves else other
        
        integer(4),intent(in)::parameters_for_detect_size
        integer(4),intent(in)::parameters_layer(parameters_for_detect_size)
        character(len=3),intent(in)::parameters_type(parameters_for_detect_size)!types: "h","rho","Cp","Cs"
        real(8),intent(in)::parameters_min(parameters_for_detect_size)!start parameters value for gradient descent method
        real(8),intent(in)::dparameters(parameters_for_detect_size)!in [min, min+d, .., last value <= max]
        real(8),intent(in)::parameters_max(parameters_for_detect_size)
        
        integer(4),intent(in)::res_max_size
        real(8),intent(out)::res(res_max_size,parameters_for_detect_size)![(parameter_1_value,parameter_2_value,...),...]
        real(8),intent(out)::res_value_of_right(res_max_size)
        integer(4),intent(out)::res_size
        
        integer(4) current_layer,layer_type
        logical(1) layer_parameter_last(parameter_type_count)
        real(8) v,x(parameters_for_detect_size)
        integer(4) i,j,k,rho_detect_count
        
        if(.not.get_anisotropic_type()==material_isotropic_type()) &
            call print_error("material_properties_by_dispersion_curves.find_material_properties",&
            ".not.get_anisotropic_type()==material_isotropic_type()")
        
        if(dispersion_curves_size<3) call print_error("material_properties_by_dispersion_curves.find_material_properties",&
            "dispersion_curves_size<3")
        do i=1,dispersion_curves_size
            if(.not.(epsilon<dispersion_curves(i,1).and.dispersion_curves(i,1)<=100d0)) &
                call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                    ".not.(epsilon<dispersion_curves(i,1).and.dispersion_curves(i,1)<=100d0)")
            if(.not.(epsilon<dispersion_curves(i,2).and.dispersion_curves(i,2)<=100d0)) &
                call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                    ".not.(epsilon<dispersion_curves(i,2).and.dispersion_curves(i,2)<=100d0)")
        enddo
        do i=2,dispersion_curves_size
            if(dispersion_curves(i-1,1)>dispersion_curves(i,1)) &
                call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                    "dispersion_curves(i-1,1)>dispersion_curves(i,1)        dispersion curves omega points must be sorted")
        enddo
        
        do i=1,parameters_for_detect_size
            if(.not.(1.le.parameters_layer(i).and.parameters_layer(i).le.get_number_of_layers())) &
                call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                    ".not.(1.le.parameters_layer(i).and.parameters_layer(i).le.get_number_of_layers())")
            do j=1,i-1
                if(parameters_layer(i)<parameters_layer(j)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_layer(i)<parameters_layer(j),i>j        parameter_layer must be sorted from the smallest to the biggest")
            enddo
        enddo
        if(parameters_for_detect_size<1) call print_error("material_properties_by_dispersion_curves.find_material_properties",&
            "parameters_for_detect_size<1")
        rho_detect_count=0
        do i=1,parameters_for_detect_size
            if(.not.(parameters_type(i)=="h".or.parameters_type(i)=="rho".or.parameters_type(i)=="Cp".or.parameters_type(i)=="Cs")) &
                call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                    ".not.(parameters_type(i)=='h'.or.(parameters_type(i)=='rho'.or.(parameters_type(i)=='Cp'.or.(parameters_type(i)=='Cs')")
            
            if(parameters_type(i)=="h") then
                if(.not.(epsilon<parameters_min(i).and.parameters_min(i)<10d0)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='h'.and..not.(epsilon<parameters_min(i).and.parameters_min(i)<10d0)")
                if(.not.(epsilon<parameters_max(i).and.parameters_max(i)<10d0)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='h'.and..not.(epsilon<parameters_max(i).and.parameters_max(i)<10d0)")
                if(.not.(epsilon<dparameters(i).and.dparameters(i)<2d0)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='h'.and..not.(epsilon<dparameters(i).and.dparameters(i)<2d0)")
                if(parameters_min(i)>parameters_max(i)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='h'.and.parameters_min(i)>parameters_max(i)")
            endif
            
            if(parameters_type(i)=="rho") then
                rho_detect_count=rho_detect_count+1
                if(.not.(epsilon<parameters_min(i).and.parameters_min(i)<50d0)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='rho'.and..not.(epsilon<parameters_min(i).and.parameters_min(i)<50d0)")
                if(.not.(epsilon<parameters_max(i).and.parameters_max(i)<50d0)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='rho'.and..not.(epsilon<parameters_max(i).and.parameters_max(i)<50d0)")
                if(.not.(epsilon<dparameters(i).and.dparameters(i)<2d0)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='rho'.and..not.(epsilon<dparameters(i).and.dparameters(i)<2d0)")
                if(parameters_min(i)>parameters_max(i)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='rho'.and.parameters_min(i)>parameters_max(i)")
            endif
            
            !if(parameters_type(i)=="E") then
            !    if(.not.(epsilon<parameters_min(i).and.parameters_min(i)<1000d0)) &
            !        call print_error("material_properties_by_dispersion_curves.find_material_properties",&
            !            "parameters_type(i)=='E'.and..not.(epsilon<parameters_min(i).and.parameters_min(i)<1000d0)")
            !    if(.not.(epsilon<parameters_max(i).and.parameters_max(i)<1000d0)) &
            !        call print_error("material_properties_by_dispersion_curves.find_material_properties",&
            !            "parameters_type(i)=='E'.and..not.(epsilon<parameters_max(i).and.parameters_max(i)<1000d0)")
            !    if(.not.(epsilon<dparameters(i).and.dparameters(i)<10d0)) &
            !        call print_error("material_properties_by_dispersion_curves.find_material_properties",&
            !            "parameters_type(i)=='E'.and..not.(epsilon<dparameters(i).and.dparameters(i)<10d0)")
            !    if(parameters_min(i)>parameters_max(i)) &
            !        call print_error("material_properties_by_dispersion_curves.find_material_properties",&
            !            "parameters_type(i)=='E'.and.parameters_min(i)>parameters_max(i)")
            !endif
            !
            !if(parameters_type(i)=="nu") then
            !    if(.not.(0d0<parameters_min(i).and.parameters_min(i)<0.5d0)) &
            !        call print_error("material_properties_by_dispersion_curves.find_material_properties",&
            !            "parameters_type(i)=='nu'.and..not.(0d0<parameters_min(i).and.parameters_min(i)<0.5d0)")
            !    if(.not.(0d0<parameters_max(i).and.parameters_max(i)<0.5d0)) &
            !        call print_error("material_properties_by_dispersion_curves.find_material_properties",&
            !            "parameters_type(i)=='nu'.and..not.(0d0<parameters_max(i).and.parameters_max(i)<0.5d0)")
            !    if(.not.(epsilon<dparameters(i).and.dparameters(i)<0.2d0)) &
            !        call print_error("material_properties_by_dispersion_curves.find_material_properties",&
            !            "parameters_type(i)=='nu'.and..not.(epsilon<dparameters(i).and.dparameters(i)<0.2d0)")
            !    if(parameters_min(i)>parameters_max(i)) &
            !        call print_error("material_properties_by_dispersion_curves.find_material_properties",&
            !            "parameters_type(i)=='nu'.and.parameters_min(i)>parameters_max(i)")
            !endif
            
            if(parameters_type(i)=="Cp") then
                if(.not.(epsilon<parameters_min(i).and.parameters_min(i)<1000d0)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='Cp'.and..not.(epsilon<parameters_min(i).and.parameters_min(i)<1000d0)")
                if(.not.(epsilon<parameters_max(i).and.parameters_max(i)<1000d0)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='Cp'.and..not.(epsilon<parameters_max(i).and.parameters_max(i)<1000d0)")
                if(.not.(epsilon<dparameters(i).and.dparameters(i)<10d0)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='Cp'.and..not.(epsilon<dparameters(i).and.dparameters(i)<10d0)")
                if(parameters_min(i)>parameters_max(i)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='Cp'.and.parameters_min(i)>parameters_max(i)")
            endif
            
            if(parameters_type(i)=="Cs") then
                if(.not.(epsilon<parameters_min(i).and.parameters_min(i)<1000d0)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='Cs'.and..not.(epsilon<parameters_min(i).and.parameters_min(i)<1000d0)")
                if(.not.(epsilon<parameters_max(i).and.parameters_max(i)<1000d0)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='Cs'.and..not.(epsilon<parameters_max(i).and.parameters_max(i)<1000d0)")
                if(.not.(epsilon<dparameters(i).and.dparameters(i)<10d0)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='Cs'.and..not.(epsilon<dparameters(i).and.dparameters(i)<10d0)")
                if(parameters_min(i)>parameters_max(i)) &
                    call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                        "parameters_type(i)=='Cs'.and.parameters_min(i)>parameters_max(i)")
            endif
        enddo
        if(rho_detect_count==get_number_of_layers())  call print_error("material_properties_by_dispersion_curves.find_material_properties",&
             "rho_detect_count==get_number_of_layers() at least 1 layer must have known rho due to impossible to restore rho value (only rho_i/rho_j where i,j - layer numbers)")
        current_layer=-1
        do i=1,parameters_for_detect_size
            if(.not.current_layer==parameters_layer(i)) then
                current_layer=parameters_layer(i)
                do j=1,parameter_type_count
                    layer_parameter_last(j)=.false.
                enddo
            endif
            
            if(layer_parameter_last(get_parameter_type_number(parameters_type(i)))) &
                call print_error("material_properties_by_dispersion_curves.find_material_properties",&
                    "a parameter is set for one layer several times")
            
            layer_parameter_last(get_parameter_type_number(parameters_type(i)))=.true.
        enddo
        
        if(res_max_size<1) call print_error("material_properties_by_dispersion_curves.find_material_properties",&
            "res_max_size<1")
        
        !call find_strict_minimum_points_on_surf(parameters_for_detect_size,surf_by_dispersion_curves,&
        !    parameters_min,dparameters,parameters_max,&
        !    res_max_size,res,res_size,res_value_of_right,skip_point_function)
        call find_strict_minimum_points_on_surf(parameters_for_detect_size,surf_by_K,&
            parameters_min,dparameters,parameters_max,&
            res_max_size,res,res_size,res_value_of_right,skip_point_function)
        
        if(all_disoersion_points_is_difenetly_true) then
            do i=1,res_size
                do j=1,parameters_for_detect_size
                    x(j)=res(i,j)
                enddo
                res_value_of_right(i)=surf_by_dispersion_curves(x)
            enddo
        endif
        
        do i=1,res_size
            do j=i+1,res_size
                if(res_value_of_right(i)<=res_value_of_right(j)) cycle
                
                v=res_value_of_right(i)
                res_value_of_right(i)=res_value_of_right(j)
                res_value_of_right(j)=v
                do k=1,parameters_for_detect_size
                    v=res(i,k)
                    res(i,k)=res(j,k)
                    res(j,k)=v
                enddo
            enddo
        enddo
        
        contains
            real(8) function surf_by_K(x) result(f)
            use main_parameters,only:set_omega,get_omega,&
                lambda=>get_layer_lambda,mu=>get_layer_mu,E=>get_layer_E,nu=>get_layer_nu,&
                Cp=>get_layer_Cp,Cs=>get_layer_Cs,rho=>get_Layer_rho,h=>get_layer_h,&
                set_layer_Cp_Cs,set_layer_rho,set_layer_h
            use count_K,only:K
            use math,only:epsilon,c0
            implicit none
                real(8),intent(in)::x(parameters_for_detect_size)
                
                real(8) omega_current
                
                integer(4) current_layer
                complex(8) Cpv,Csv
                real(8) rhov
                real(8) v,fine_for_forbit
                integer(4) i,j
                
                f=0d0
                
                current_layer=-1
                Cpv=-1d0
                Csv=-1d0
                rhov=-1d0
                fine_for_forbit=0d0
                do i=1,parameters_for_detect_size
                    if(.not.current_layer==parameters_layer(i)) then
                        if(.not.current_layer==-1.and..not.(Cpv==-1d0.and.Csv==-1d0.and.rhov==-1d0)) then
                            if(Cpv==-1d0) Cpv=Cp(current_layer)
                            if(Csv==-1d0) Csv=Cs(current_layer)
                            if(rhov==-1d0) rhov=rho(current_layer)
                            if(real(Csv)*sqrt(2d0)>real(Cpv)) then
                                fine_for_forbit=fine_for_forbit+(Csv*sqrt(2d0)-Cpv)
                                Csv=Cpv/sqrt(2d0)-epsilon
                            endif
                            call set_layer_rho(current_layer,rhov)
                            call set_layer_Cp_Cs(current_layer,Cpv,Csv)
                            Cpv=-1d0
                            Csv=-1d0
                            rhov=-1d0
                        endif
                        
                        current_layer=parameters_layer(i)
                    endif
                    
                    if(parameters_type(i)=="h") then
                        call set_layer_h(current_layer,x(i))
                        
                        cycle
                    endif
                    if(parameters_type(i)=="rho") then
                        rhov=x(i)
                        
                        cycle
                    endif
                    if(parameters_type(i)=="Cp") then
                        Cpv=x(i)
                        
                        cycle
                    endif
                    if(parameters_type(i)=="Cs") then
                        Csv=x(i)
                        
                        cycle
                    endif
                enddo
                if(.not.current_layer==-1.and..not.(Cpv==-1d0.and.Csv==-1d0.and.rhov==-1d0)) then
                    if(Cpv==-1d0) Cpv=Cp(current_layer)
                    if(Csv==-1d0) Csv=Cs(current_layer)
                    if(rhov==-1d0) rhov=rho(current_layer)
                    if(real(Csv)*sqrt(2d0)>real(Cpv)) then
                        fine_for_forbit=fine_for_forbit+(Csv*sqrt(2d0)-Cpv)
                        Csv=Cpv/sqrt(2d0)-epsilon
                    endif
                    call set_layer_rho(current_layer,rhov)
                    call set_layer_Cp_Cs(current_layer,Cpv,Csv)
                endif
                
                omega_current=-1d0
                do i=1,dispersion_curves_size
                    if(abs(omega_current-dispersion_curves(i,1))>epsilon) then
                        call set_omega(dispersion_curves(i,1))
                        omega_current=get_omega()
                    endif
                    
                    f=f+1d0/max(1d-5,abs(K(i=3,j=3,alpha=(dispersion_curves(i,2)+c0),beta=c0,z=0d0)))
                enddo
                
                
                f=f/dispersion_curves_size
                !f=f*(1d0+fine_for_forbit)
                f=f+fine_for_forbit*100d0
            endfunction surf_by_K
            real(8) function surf_by_dispersion_curves(x) result(f)
            use main_parameters,only:set_omega,get_omega,&
                lambda=>get_layer_lambda,mu=>get_layer_mu,E=>get_layer_E,nu=>get_layer_nu,&
                Cp=>get_layer_Cp,Cs=>get_layer_Cs,rho=>get_Layer_rho,h=>get_layer_h,&
                set_layer_Cp_Cs,set_layer_rho,set_layer_h
            
            use dispersion_curves_for_K,only:count_poles
            use math,only:epsilon,c0
            implicit none
                real(8),intent(in)::x(parameters_for_detect_size)
                
                real(8) omega_current
                
                real(8) res(25)
                integer(4) res_size
                
                integer(4) current_layer
                complex(8) Cpv,Csv
                real(8) rhov
                real(8) v,fine_for_forbit
                integer(4) i,j
                
                f=0d0
                
                current_layer=-1
                Cpv=-1d0
                Csv=-1d0
                rhov=-1d0
                fine_for_forbit=0d0
                do i=1,parameters_for_detect_size
                    if(.not.current_layer==parameters_layer(i)) then
                        if(.not.current_layer==-1.and..not.(Cpv==-1d0.and.Csv==-1d0.and.rhov==-1d0)) then
                            if(Cpv==-1d0) Cpv=Cp(current_layer)
                            if(Csv==-1d0) Csv=Cs(current_layer)
                            if(rhov==-1d0) rhov=rho(current_layer)
                            if(real(Csv)*sqrt(2d0)>real(Cpv)) then
                                fine_for_forbit=fine_for_forbit+(Csv*sqrt(2d0)-Cpv)
                                Csv=Cpv/sqrt(2d0)-epsilon
                            endif
                            call set_layer_rho(current_layer,rhov)
                            call set_layer_Cp_Cs(current_layer,Cpv,Csv)
                            Cpv=-1d0
                            Csv=-1d0
                            rhov=-1d0
                        endif
                        
                        current_layer=parameters_layer(i)
                    endif
                    
                    if(parameters_type(i)=="h") then
                        call set_layer_h(current_layer,x(i))
                        
                        cycle
                    endif
                    if(parameters_type(i)=="rho") then
                        rhov=x(i)
                        
                        cycle
                    endif
                    if(parameters_type(i)=="Cp") then
                        Cpv=x(i)
                        
                        cycle
                    endif
                    if(parameters_type(i)=="Cs") then
                        Csv=x(i)
                        
                        cycle
                    endif
                enddo
                if(.not.current_layer==-1.and..not.(Cpv==-1d0.and.Csv==-1d0.and.rhov==-1d0)) then
                    if(Cpv==-1d0) Cpv=Cp(current_layer)
                    if(Csv==-1d0) Csv=Cs(current_layer)
                    if(rhov==-1d0) rhov=rho(current_layer)
                    if(real(Csv)*sqrt(2d0)>real(Cpv)) then
                        fine_for_forbit=fine_for_forbit+(Csv*sqrt(2d0)-Cpv)
                        Csv=Cpv/sqrt(2d0)-epsilon
                    endif
                    call set_layer_rho(current_layer,rhov)
                    call set_layer_Cp_Cs(current_layer,Cpv,Csv)
                endif
                
                omega_current=-1d0
                do i=1,dispersion_curves_size
                    if(abs(omega_current-dispersion_curves(i,1))>epsilon) then
                        call set_omega(dispersion_curves(i,1))
                        omega_current=get_omega()
                        
                        call count_poles(phi=0d0,number_of_en=3,res_max_size=25,res=res,res_size=res_size)
                    endif
                    
                    v=1d100
                    do j=1,res_size
                        !if(abs(res(j)-dispersion_curves(i,2))>3d0) cycle
                        
                        v=min(v,abs(real(res(j))-dispersion_curves(i,2)))
                    enddo
                    !f=f+v
                    f=max(f,v)
                enddo
                
                !f=f/dispersion_curves_size
                !f=f*(1d0+fine_for_forbit)
                f=f+fine_for_forbit*100d0
            endfunction surf_by_dispersion_curves
            
            logical(1) function skip_point_function(x) result(f)
            use main_parameters,only:get_layer_Cp,get_layer_Cs,get_layer_rho,&
                pure_check_correct_isotropic_parameters_for_Cp_Cs
            use math,only:c0
            implicit none
                real(8),intent(in)::x(parameters_for_detect_size)
                
                integer(4) layer
                real(8) Cp,Cs,rho
                
                integer(4) i
                
                !integer(4),intent(in)::parameters_for_detect_size
                !integer(4),intent(in)::parameters_layer(parameters_for_detect_size)
                !character(len=3),intent(in)::parameters_type(parameters_for_detect_size)!types: "h","rho","Cp","Cs"
                !real(8),intent(in)::parameters_min(parameters_for_detect_size)!start parameters value for gradient descent method
                !real(8),intent(in)::dparameters(parameters_for_detect_size)!in [min, min+d, .., last value <= max]
                !real(8),intent(in)::parameters_max(parameters_for_detect_size)
                
                layer=parameters_layer(1)
                Cp=-1d0
                Cs=-1d0
                rho=-1d0
                do i=1,parameters_for_detect_size
                    if(.not.layer==parameters_layer(i)) then
                        if(real(Cp)==-1d0) Cp=real(get_layer_Cp(layer))
                        if(real(Cs)==-1d0) Cs=real(get_layer_Cs(layer))
                        if(rho==-1d0) rho=get_layer_rho(layer)
                        
                        if(.not.pure_check_correct_isotropic_parameters_for_Cp_Cs(Cp+c0,Cs+c0,rho)) then
                            f=.false.
                            return
                        endif
                        
                        Cp=-1d0
                        Cs=-1d0
                        rho=-1d0
                    endif
                    
                    if(trim(parameters_type(i))=="Cp") Cp=x(i)
                    if(trim(parameters_type(i))=="Cs") Cs=x(i)
                    if(trim(parameters_type(i))=="rho") rho=x(i)
                enddo
                
                if(.not.layer==parameters_layer(i)) then
                    if(real(Cp)==-1d0) Cp=real(get_layer_Cp(layer))
                    if(real(Cs)==-1d0) Cs=real(get_layer_Cs(layer))
                    if(rho==-1d0) rho=get_layer_rho(layer)
                    
                    if(.not.pure_check_correct_isotropic_parameters_for_Cp_Cs(Cp+c0,Cs+c0,rho)) then
                        f=.false.
                        return
                    endif
                    
                    Cp=-1d0
                    Cs=-1d0
                    rho=-1d0
                endif
                
                f=.true.
            endfunction skip_point_function
        endsubroutine find_material_properties
        
        integer(4) function get_parameter_type_number(parameter_name) result(f)
        use system,only:print_error
        implicit none
            character(len=3),intent(in)::parameter_name
            
            if(parameter_name=="h") then
                f=1
                return
            endif
            if(parameter_name=="rho") then
                f=2
                return
            endif
            if(parameter_name=="Cp") then
                f=3
                return
            endif
            if(parameter_name=="Cs") then
                f=4
                return
            endif
            
            call print_error("material_properties_by_dispersion_curves.get_parameter_type_number","an incorrect parameter name")
        endfunction get_parameter_type_number
        character(len=3) function get_parameter_type_name(parameter_number) result(f)
        use system,only:print_error
        implicit none
            integer(4),intent(in)::parameter_number
            
            if(parameter_number==1) then
                f="h"
                return
            endif
            if(parameter_number==2) then
                f="rho"
                return
            endif
            if(parameter_number==3) then
                f="Cp"
                return
            endif
            if(parameter_number==4) then
                f="Cs"
                return
            endif
            
            call print_error("material_properties_by_dispersion_curves.get_parameter_type_name","an incorrect parameter number")
        endfunction get_parameter_type_name
endmodule material_properties_by_dispersion_curves
