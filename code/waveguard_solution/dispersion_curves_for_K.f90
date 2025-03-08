module dispersion_curves_for_K
implicit none
    public::count_poles,count_complex_poles

    real(8)::tmin=1d-2,tmax=35d0,ht=1d-3,eps=1d-7
    
    real(8)::Re_poles_min_value=1d-2,Re_poles_dvalue=5d-1,Re_poles_max_value=70d0
    real(8)::Im_poles_min_value=-3d0,Im_poles_dvalue=1d-1,Im_poles_max_value=3d0
    
    private::tmin,tmax,ht,eps
    private::count_poles_for_matrix_real_poles,count_poles_for_matrix_complex_poles
    
    contains
    
    subroutine count_poles(phi,number_of_en,res_max_size,res,res_size)
    use main_parameters,only:get_anisotropic
    use system,only:print_error
    use math,only:pi
    implicit none
        real(8),intent(in)::phi
        integer(4),intent(in)::number_of_en
        integer(4),intent(in)::res_max_size
        real(8),intent(out)::res(res_max_size)
        integer(4),intent(out)::res_size
        
        real(8) res_AX(res_max_size)
        integer(4) res_AX_size
        
        integer(4) i
        
        if(phi<0d0.or.phi>=pi+pi) call print_error("dispersion_curves_for_K.count_poles","phi<0d0.or.phi>=pi+pi")
        if(number_of_en<1.or.number_of_en>3) call print_error("dispersion_curves_for_K.count_poles","number_of_en<1.or.number_of_en>3")
        if(res_max_size<1) call print_error("dispersion_curves_for_K.count_poles","res_max_size<1")
        
        if(get_anisotropic()==0) then
            call count_poles_for_matrix_real_poles(in_f_with_AY,phi,number_of_en,res_max_size,res,res_size)
            call count_poles_for_matrix_real_poles(in_f_with_AX,phi,number_of_en,res_max_size,res_AX,res_AX_size)
            if(res_AX_size==res_max_size) then
                res_AX_size=0
            endif
            
            do i=res_size+1,min(res_max_size,res_size+res_AX_size)
                res(i)=res_AX(i-res_size)
            enddo
            res_size=min(res_max_size,res_size+res_AX_size)
        else
            call count_poles_for_matrix_real_poles(in_f_with_A,phi,number_of_en,res_max_size,res,res_size)
        endif
        
    contains
        complex(8) function in_f_with_A(alpha) result(f)
        use count_K,only:get_det_A
        implicit none
            complex(8),intent(in)::alpha
            
            f=get_det_A(alpha*cos(phi),alpha*sin(phi),number_of_en)
        endfunction in_f_with_A
        complex(8) function in_f_with_AY(alpha) result(f)
        use count_K,only:get_det_AY
        implicit none
            complex(8),intent(in)::alpha
            
            f=get_det_AY(alpha*cos(phi),alpha*sin(phi),number_of_en)
        endfunction in_f_with_AY
        complex(8) function in_f_with_AX(alpha) result(f)
        use count_K,only:get_det_AX
        implicit none
            complex(8),intent(in)::alpha
            
            f=get_det_AX(alpha*cos(phi),alpha*sin(phi),number_of_en)
        endfunction in_f_with_AX
    endsubroutine count_poles
    subroutine count_poles_for_matrix_real_poles(in_f,phi,number_of_en,res_max_size,res,res_size)
    use Halfc_,only:halfc
    implicit none
        complex(8),external::in_f!complex(8) function in_f(complex(8) alpha)
        real(8),intent(in)::phi
        integer(4),intent(in)::number_of_en
        integer(4),intent(in)::res_max_size
        real(8),intent(out)::res(res_max_size)
        integer(4),intent(out)::res_size
        
        call halfc(in_f,tmin,tmax,ht,eps,res_max_size,res,res_size)
    endsubroutine count_poles_for_matrix_real_poles
    
    subroutine count_complex_poles(phi,number_of_en,res_max_size,res,res_size)
    use main_parameters,only:get_anisotropic
    use count_K,only:K
    use system,only:print_error
    use math,only:pi
    implicit none
        real(8),intent(in)::phi
        integer(4),intent(in)::number_of_en
        integer(4),intent(in)::res_max_size
        complex(8),intent(out)::res(res_max_size)
        integer(4),intent(out)::res_size
        
        complex(8) res_AX(res_max_size)
        integer(4) res_AX_size
        
        complex(8) Kvalue
        integer(4) i,res_size_
        
        if(phi<0d0.or.phi>=pi+pi) call print_error("dispersion_curves_for_K.count_complex_poles","phi<0d0.or.phi>=pi+pi")
        if(number_of_en<1.or.number_of_en>3) call print_error("dispersion_curves_for_K.count_complex_poles","number_of_en<1.or.number_of_en>3")
        if(res_max_size<1) call print_error("dispersion_curves_for_K.count_complex_poles","res_max_size<1")
        
        if(get_anisotropic()==0) then
            call count_poles_for_matrix_complex_poles(in_f_with_AY,phi,number_of_en,res_max_size,res,res_size)
            call count_poles_for_matrix_complex_poles(in_f_with_AX,phi,number_of_en,res_max_size,res_AX,res_AX_size)
            if(res_AX_size==res_max_size) then
                res_AX_size=0
            endif
            
            do i=res_size+1,min(res_max_size,res_size+res_AX_size)
                res(i)=res_AX(i-res_size)
            enddo
            res_size=min(res_max_size,res_size+res_AX_size)
        else
            call count_poles_for_matrix_complex_poles(in_f_with_A,phi,number_of_en,res_max_size,res,res_size)
        endif
        
        res_size_=res_size
        res_size=0
        do i=1,res_size_
            Kvalue=K(i=3,j=3,alpha=res(i)*cos(phi),beta=res(i)*sin(phi),z=0d0)
            if(abs(Kvalue)<1d1) cycle
            
            res_size=res_size+1
            res(res_size)=res(i)
        enddo
        
    contains
        complex(8) function in_f_with_A(alpha) result(f)
        use count_K,only:get_det_A
        implicit none
            complex(8),intent(in)::alpha
            
            f=get_det_A(alpha*cos(phi),alpha*sin(phi),number_of_en)
        endfunction in_f_with_A
        complex(8) function in_f_with_AY(alpha) result(f)
        use count_K,only:get_det_AY
        implicit none
            complex(8),intent(in)::alpha
            
            f=get_det_AY(alpha*cos(phi),alpha*sin(phi),number_of_en)
        endfunction in_f_with_AY
        complex(8) function in_f_with_AX(alpha) result(f)
        use count_K,only:get_det_AX
        implicit none
            complex(8),intent(in)::alpha
            
            f=get_det_AX(alpha*cos(phi),alpha*sin(phi),number_of_en)
        endfunction in_f_with_AX
    endsubroutine count_complex_poles
    subroutine count_poles_for_matrix_complex_poles(in_f,phi,number_of_en,res_max_size,res,res_size)
    use gradient_descent_method,only:find_strict_minimum_points_on_surf
    use math,only:ci
    implicit none
        complex(8),external::in_f!complex(8) function in_f(complex(8) alpha)
        real(8),intent(in)::phi
        integer(4),intent(in)::number_of_en
        integer(4),intent(in)::res_max_size
        complex(8),intent(out)::res(res_max_size)
        integer(4),intent(out)::res_size
        
        real(8) x0_min(2),dx0(2),x0_max(2),res_(res_max_size,2)
        
        integer(4) i
        
        x0_min(1)=Re_poles_min_value; dx0(1)=Re_poles_dvalue; x0_max(1)=Re_poles_max_value
        x0_min(2)=Im_poles_min_value; dx0(2)=Im_poles_dvalue; x0_max(2)=Im_poles_max_value
        
        call find_strict_minimum_points_on_surf(2,surf,x0_min,dx0,x0_max,res_max_size,res_,res_size)
        
        do i=1,res_size
            res(i)=res_(i,1)+res_(i,2)*ci
        enddo
    contains
        real(8) function surf(x) result(f)
        use math,only:ci
        implicit none
            real(8) x(2)![Re alpha, Im alpha]
            
            f=abs(in_f(x(1)+x(2)*ci))
        endfunction surf
    endsubroutine count_poles_for_matrix_complex_poles
endmodule dispersion_curves_for_K
