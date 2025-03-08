module dispersion_curves_for_K
implicit none
    public::count_poles

    real(8)::tmin=1d-2,tmax=35d0,ht=1d-3,eps=1d-7
    
    private::tmin,tmax,ht,eps
    private::count_poles_for_matrix
    
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
        
        if(phi<0d0.or.phi>=pi+pi) call print_error("dispersion_curves_for_K.dispersion_curves_for_K","phi<0d0.or.phi>=pi+pi")
        if(number_of_en<1.or.number_of_en>3) call print_error("dispersion_curves_for_K.dispersion_curves_for_K","number_of_en<1.or.number_of_en>3")
        if(res_max_size<1) call print_error("dispersion_curves_for_K.dispersion_curves_for_K","res_max_size<1")
        
        if(get_anisotropic()==0) then
            call count_poles_for_matrix(in_f_with_AY,phi,number_of_en,res_max_size,res,res_size)
            call count_poles_for_matrix(in_f_with_AX,phi,number_of_en,res_max_size,res_AX,res_AX_size)
            if(res_AX_size==res_max_size) then
                res_AX_size=0
            endif
            
            do i=res_size+1,min(res_max_size,res_size+res_AX_size)
                res(i)=res_AX(i-res_size)
            enddo
            res_size=min(res_max_size,res_size+res_AX_size)
        else
            call count_poles_for_matrix(in_f_with_A,phi,number_of_en,res_max_size,res,res_size)
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
    
    subroutine count_poles_for_matrix(in_f,phi,number_of_en,res_max_size,res,res_size)
    use Halfc_,only:halfc
    implicit none
        complex(8),external::in_f!complex(8) function in_f(complex(8) alpha)
        real(8),intent(in)::phi
        integer(4),intent(in)::number_of_en
        integer(4),intent(in)::res_max_size
        real(8),intent(out)::res(res_max_size)
        integer(4),intent(out)::res_size
        
        call halfc(in_f,tmin,tmax,ht,eps,res_max_size,res,res_size)
    endsubroutine count_poles_for_matrix
endmodule dispersion_curves_for_K
