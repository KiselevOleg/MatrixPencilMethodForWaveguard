module gradient_descent_method
implicit none
    public::find_strict_minimum_points_on_surf,find_zeros_on_notnegative_surf
    
    interface find_strict_minimum_points_on_surf
        procedure find_strict_minimum_points_on_surf_with_values
        procedure find_strict_minimum_points_on_surf_without_values
    endinterface find_strict_minimum_points_on_surf
    
    real(8),private::res_min_diff_for_elements=1d-4
    
    real(8),private::derivative_delta=1d-4
    real(8),private::dx_min_for_step=1d-6,dx_max_for_step=1d-1
    
    real(8),private::grad_epsilon=1d-6
    
    real(8),private::max_value_for_surf_to_be_zero=1d-1
    
    private::find_strict_minimum_points_on_surf_with_values,find_strict_minimum_points_on_surf_without_values
    contains

    subroutine find_strict_minimum_points_on_surf_with_values(space_dimension,surf,x0_min,dx0,x0_max,res_size_max,res,res_size,surf_value_at_res)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        integer(4),intent(in)::space_dimension
        real(8),external::surf!real(8) function surf(x(space_dimension))
        real(8),intent(in)::x0_min(space_dimension)
        real(8),intent(in)::dx0(space_dimension)
        real(8),intent(in)::x0_max(space_dimension)
        integer(4),intent(in)::res_size_max
        real(8),intent(out)::res(res_size_max,space_dimension)
        integer(4),intent(out)::res_size
        real(8),intent(out)::surf_value_at_res(res_size_max)
        
        real(8) x0(space_dimension),x(space_dimension),x_(space_dimension)
        real(8) grad(space_dimension),abs_dx
        
        logical(1) enumeration_of_a_grid_of_values_started
        
        real(8) surf_value_at_x,surf_value_at_x_
        logical(1) correct_point
        integer(4) i,j,k
        
        if(space_dimension<1) call print_error("gradient_descent_method.find_strict_minimum_points_on_surf_with_values",&
            "space_dimension<1")
        do i=1,space_dimension
            if(x0_min(i)>=x0_max(i)) call print_error("gradient_descent_method.find_strict_minimum_points_on_surf_with_values",&
                "x0_min(i)>=x0_max(i)")
            if(dx0(i)<epsilon) call print_error("gradient_descent_method.find_strict_minimum_points_on_surf_with_values",&
                "dx0(i)<epsilon")
            if(dx0(i)>=x0_max(i)-x0_min(i)) call print_error("gradient_descent_method.find_strict_minimum_points_on_surf_with_values",&
                "dx0(i)>=x0_max(i)-x0_min(i)")
        enddo
        if(res_size_max<1) call print_error("gradient_descent_method.find_strict_minimum_points_on_surf_with_values",&
            "res_size_max<1")
        
        res_size=0
        enumeration_of_a_grid_of_values_started=.false.
        do while(next_grid_of_values_point())
            do i=1,space_dimension
                x(i)=x0(i)
            enddo
            
            if(.not.get_surf_value(x,surf_value_at_x)) cycle
            correct_point=.true.
            do while(.true.)
                abs_dx=0d0
                do i=1,space_dimension
                    x(i)=x(i)-derivative_delta
                    if(.not.get_surf_value(x,grad(i))) then
                        correct_point=.false.
                        exit
                    endif
                    x(i)=x(i)+derivative_delta+derivative_delta
                    if(.not.get_surf_value(x,surf_value_at_x_)) then
                        correct_point=.false.
                        exit
                    endif
                    grad(i)=(grad(i)-surf_value_at_x_)/derivative_delta*0.5d0
                    x(i)=x(i)-derivative_delta
                    abs_dx=abs_dx+grad(i)*grad(i)
                    
                    !x(i)=x(i)-derivative_delta
                    !grad(i)=f(x)
                    !x(i)=x(i)+derivative_delta+derivative_delta
                    !grad(i)=(grad(i)-f(x))/derivative_delta*0.5d0
                    !x(i)=x(i)-derivative_delta
                    !abs_dx=abs_dx+grad(i)*grad(i)
                enddo
                if(.not.correct_point) exit
                abs_dx=sqrt(abs_dx)
                if(abs_dx<grad_epsilon) then
                    if(.not.add_x_to_res()) return
                    
                    exit
                endif
                do i=1,space_dimension
                    grad(i)=grad(i)/abs_dx
                enddo
                
                surf_value_at_x_=surf_value_at_x+1d0
                abs_dx=dx_max_for_step+dx_max_for_step
                do while(surf_value_at_x<=surf_value_at_x_)
                    abs_dx=abs_dx*0.5d0
                    if(abs_dx<dx_min_for_step) then
                        abs_dx=0d0
                        surf_value_at_x_=surf_value_at_x
                        exit
                    endif
                    
                    do i=1,space_dimension
                        x_(i)=x(i)+grad(i)*abs_dx
                    enddo
                    if(.not.get_surf_value(x_,surf_value_at_x_)) then
                        correct_point=.false.
                        exit
                    endif
                enddo
                if(.not.correct_point) exit
                
                do i=1,space_dimension
                    x(i)=x_(i)
                enddo
                surf_value_at_x=surf_value_at_x_
                if(abs_dx<grad_epsilon) then
                    if(.not.add_x_to_res()) return
                    
                    exit
                endif
                
                if(.not.check_x_in_x0_area()) exit
            enddo
        enddo
        
    contains
        logical(1) function next_grid_of_values_point() result(f)
        implicit none
            integer(4) i
            
            if(.not.enumeration_of_a_grid_of_values_started) then
                enumeration_of_a_grid_of_values_started=.true.
                
                do i=1,space_dimension
                    x0(i)=x0_min(i)+derivative_delta+derivative_delta
                enddo
                
                f=.true.
                return
            endif
            
            do i=1,space_dimension
                if(x0(i)+dx0(i)<=x0_max(i)) then
                    x0(i)=x0(i)+dx0(i)
                    
                    f=.true.
                    return
                endif
                
                x0(i)=x0_min(i)+derivative_delta+derivative_delta
            enddo
            
            f=.false.
        endfunction next_grid_of_values_point
        
        logical(1) function add_x_to_res() result(f)
        implicit none
            real(8) diff
            integer(4) i,j
            
            do i=1,res_size
                diff=0d0
                do j=1,space_dimension
                    diff=max(diff,abs(res(i,j)-x(j)))
                enddo
                
                if(diff<res_min_diff_for_elements) then
                    f=.true.
                    return
                endif
            enddo
            
            res_size=res_size+1
            do i=1,space_dimension
                res(res_size,i)=x(i)
            enddo
            surf_value_at_res(res_size)=surf_value_at_x
            
            f=res_size/=res_size_max
        endfunction add_x_to_res
        
        logical(1) function get_surf_value(x,value) result(f)
        implicit none
            real(8),intent(in)::x(space_dimension)
            real(8),intent(out)::value
            
            f=check_correct(x)
            if(.not.f) return
            
            value=surf(x)
        endfunction get_surf_value
        
        logical(1) function check_correct(x) result(f)
        implicit none
            real(8),intent(in)::x(space_dimension)
            
            integer(4) i
            
            do i=1,space_dimension
                if(x0_min(i)<=x(i).and.x(i)<=x0_max(i)) cycle
                
                f=.false.
                return
            enddo
            
            f=.true.
        endfunction check_correct
        
        logical(1) function check_x_in_x0_area() result(f)
        implicit none
            integer(4) i
            
            do i=1,space_dimension
                if(x0(i)-dx0(i)*1.5d0<=x(i).and.x(i)<=x0(i)+dx0(i)*1.5d0) cycle
                
                f=.false.
                return
            enddo
            
            f=.true.
        endfunction check_x_in_x0_area
    endsubroutine find_strict_minimum_points_on_surf_with_values
    subroutine find_strict_minimum_points_on_surf_without_values(space_dimension,surf,x0_min,dx0,x0_max,res_size_max,res,res_size)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        integer(4),intent(in)::space_dimension
        real(8),external::surf!real(8) function surf(x(space_dimension))
        real(8),intent(in)::x0_min(space_dimension)
        real(8),intent(in)::dx0(space_dimension)
        real(8),intent(in)::x0_max(space_dimension)
        integer(4),intent(in)::res_size_max
        real(8),intent(out)::res(res_size_max,space_dimension)
        integer(4),intent(out)::res_size
        
        real(8) res_surf_value_at_resize(res_size_max)
        
        integer(4) i
        
        if(space_dimension<1) call print_error("gradient_descent_method.find_strict_minimum_points_on_surf_without_values",&
            "space_dimension<1")
        do i=1,space_dimension
            if(x0_min(i)>=x0_max(i)) call print_error("gradient_descent_method.find_strict_minimum_points_on_surf_without_values",&
                "x0_min(i)>=x0_max(i)")
            if(dx0(i)<epsilon) call print_error("gradient_descent_method.find_strict_minimum_points_on_surf_without_values",&
                "dx0(i)<epsilon")
            if(dx0(i)>=x0_max(i)-x0_min(i)) call print_error("gradient_descent_method.find_strict_minimum_points_on_surf_without_values",&
                "dx0(i)>=x0_max(i)-x0_min(i)")
        enddo
        if(res_size_max<1) call print_error("gradient_descent_method.find_strict_minimum_points_on_surf_without_values",&
            "res_size_max<1")
        
        call find_strict_minimum_points_on_surf_with_values(space_dimension,surf,x0_min,dx0,x0_max,res_size_max,res,res_size,res_surf_value_at_resize)
    endsubroutine find_strict_minimum_points_on_surf_without_values
    
    subroutine find_zeros_on_notnegative_surf(space_dimension,surf,x0_min,dx0,x0_max,res_size_max,res,res_size)
    use math,only:epsilon
    use system,only:print_error
    implicit none
        integer(4),intent(in)::space_dimension
        real(8),external::surf!real(8) function surf(x(space_dimension))
        real(8),intent(in)::x0_min(space_dimension)
        real(8),intent(in)::dx0(space_dimension)
        real(8),intent(in)::x0_max(space_dimension)
        integer(4),intent(in)::res_size_max
        real(8),intent(out)::res(res_size_max,space_dimension)
        integer(4),intent(out)::res_size
        
        real(8) res_(res_size_max,space_dimension)
        integer(4) res_size_
        
        real(8) surf_value(res_size_max)
        integer(4) i,j
        
        if(space_dimension<1) call print_error("gradient_descent_method.find_zeros_on_notnegative_surf",&
            "space_dimension<1")
        do i=1,space_dimension
            if(x0_min(i)>=x0_max(i)) call print_error("gradient_descent_method.find_zeros_on_notnegative_surf",&
                "x0_min(i)>=x0_max(i)")
            if(dx0(i)<epsilon) call print_error("gradient_descent_method.find_zeros_on_notnegative_surf",&
                "dx0(i)<epsilon")
            if(dx0(i)>=x0_max(i)-x0_min(i)) call print_error("gradient_descent_method.find_zeros_on_notnegative_surf",&
                "dx0(i)>=x0_max(i)-x0_min(i)")
        enddo
        if(res_size_max<1) call print_error("gradient_descent_method.find_zeros_on_notnegative_surf",&
            "res_size_max<1")
        
        call find_strict_minimum_points_on_surf(space_dimension,surf,x0_min,dx0,x0_max,res_size_max,res_,res_size_,surf_value)
        
        res_size=0
        do i=1,res_size_
            if(surf_value(i)>max_value_for_surf_to_be_zero) cycle
            
            res_size=res_size+1
            do j=1,space_dimension
                res(res_size,j)=res_(i,j)
            enddo
        enddo
    endsubroutine find_zeros_on_notnegative_surf
endmodule gradient_descent_method
