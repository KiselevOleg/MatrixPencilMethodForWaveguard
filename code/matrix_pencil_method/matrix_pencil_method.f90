module matrix_pencil_method
implicit none
    public::init_matrix_pencil_method,destructor_matrix_pencil_method
    
    public::count_dispersion_numbers
    public::get_dx_filter_strength,get_L_filter_strength,get_dL_filter_value
    public::set_dx_filter_strength,set_L_filter_strength,set_dL_filter_value
    
    public::get_result
    public::get_dx_filter_result,get_L_filter_result
    public::get_result_with_lower_filters
    
    integer(4),private::dx_filter_strength=1
    integer(4),private::L_filter_strength=1
    integer(4),private::dL_filter_value=5
    
    integer(4),private::L_
    complex(8),allocatable,private::res_(:)
    integer(4),private::res_size_
    
    complex(8),allocatable,private::dx_filter_res(:,:)!filter_number,res
    integer(4),allocatable,private::dx_filter_res_size(:)
    
    complex(8),allocatable,private::L_filter_res(:,:)!filter_number_res
    integer(4),allocatable,private::L_filter_res_size(:)
    
    logical(1),private::counted=.false.
    
    private::refresh_matrixes,count_result_with_filters
    contains
    
    subroutine count_dispersion_numbers(omega,L,res,res_size)
    use matrix_pencil_method_basis,only:count_dispersion_numbers_=>count_dispersion_numbers
    use load_experimental_measurements,only:get_Nx
    use math,only:epsilon
    use system,only:print_error
    implicit none
        complex(8),intent(in)::omega
        integer(4),intent(in)::L
        complex(8),intent(out)::res(L)
        integer(4),intent(out)::res_size
        
        integer(4) i,j
        
        if(real(omega)<epsilon) call print_error("matrix_pencil_method.count_dispersion_numbers","real(omega)<epsilon")
        if(L<3) call print_error("matrix_pencil_method.count_dispersion_numbers","L<3")
        if(L>=get_Nx()) call print_error("matrix_pencil_method.count_dispersion_numbers","L>=get_Nx()")
        if(get_Nx()<4) call print_error("matrix_pencil_method.count_dispersion_numbers","get_Nx()<4")

        counted=.false.
        L_=L
        call refresh_matrixes(L)
        
        do i=1,dx_filter_strength
            call count_dispersion_numbers_(omega,i+1,L,res_,res_size_)
            
            dx_filter_res_size(i)=res_size_
            do j=1,res_size_
                dx_filter_res(i,j)=res_(j)
            enddo
        enddo
        
        deallocate(res_)
        do i=1,L_filter_strength
            allocate(res_(L+dL_filter_value*i))
            
            call count_dispersion_numbers_(omega,1,L+dL_filter_value*i,res_,res_size_)
            
            L_filter_res_size(i)=res_size_
            do j=1,res_size_
                L_filter_res(i,j)=res_(j)
            enddo
            
            deallocate(res_)
        enddo
        allocate(res_(L))
        
        call count_dispersion_numbers_(omega,1,L,res_,res_size_)
        
        call count_result_with_filters(res,res_size,dx_filter_strength,L_filter_strength)
        counted=.true.
    endsubroutine count_dispersion_numbers
    subroutine count_result_with_filters(res,res_size,dx_filter_strength_,L_filter_strength_)
    use system,only:print_error
    implicit none
        complex(8),intent(inout)::res(L_)
        integer(4),intent(out)::res_size
        integer(4),intent(in)::dx_filter_strength_
        integer(4),intent(in)::L_filter_strength_
        
        integer(4) i,j,k
        logical(1) filter_check
        
        if(.not.(0.le.dx_filter_strength_.and.dx_filter_strength_.le.dx_filter_strength)) &
            call print_error("matrix_pencil_method.count_result_with_filters",&
            ".not.(0.le.dx_filter_strength_.and.dx_filter_strength_.le.dx_filter_strength)")
        if(.not.(0.le.L_filter_strength_.and.L_filter_strength_.le.L_filter_strength)) &
            call print_error("matrix_pencil_method.count_result_with_filters",&
            ".not.(0.le.L_filter_strength_.and.L_filter_strength_.le.L_filter_strength)")
        
        res_size=0
        do i=1,res_size_
            filter_check=dx_filter_strength_==0
            do j=1,dx_filter_strength_
                filter_check=.false.
                do k=1,dx_filter_res_size(j)
                    if(abs((res_(i)-dx_filter_res(j,k))/res_(i))<0.01d0) then
                        filter_check=.true.
                        exit
                    endif
                enddo
                if(.not.filter_check) exit
            enddo
            if(.not.filter_check) cycle
            
            filter_check=L_filter_strength_==0
            do j=1,L_filter_strength_
                filter_check=.false.
                do k=1,L_filter_res_size(j)
                    if(abs((res_(i)-L_filter_res(j,k))/res_(i))<0.01d0) then
                        filter_check=.true.
                        exit
                    endif
                enddo
                if(.not.filter_check) exit
            enddo
            if(.not.filter_check) cycle
            
            res_size=res_size+1
            res(res_size)=res_(i)
        enddo
    endsubroutine count_result_with_filters
    
    subroutine refresh_matrixes(L)
    implicit none
        integer(4),intent(in)::L
        
        deallocate(res_)
        
        deallocate(dx_filter_res)
        deallocate(dx_filter_res_size)
        
        deallocate(L_filter_res)
        deallocate(L_filter_res_size)
        
        allocate(res_(L))
        
        allocate(dx_filter_res(dx_filter_strength,L))
        allocate(dx_filter_res_size(dx_filter_strength))
        
        allocate(L_filter_res(L_filter_strength,L+5*L_filter_strength))
        allocate(L_filter_res_size(L_filter_strength))
        
        counted=.false.
    endsubroutine refresh_matrixes
    
    
    
    pure integer(4) function get_dx_filter_strength() result(f)
    implicit none
        f=dx_filter_strength
    endfunction get_dx_filter_strength
    pure integer(4) function get_L_filter_strength() result(f)
    implicit none
        f=L_filter_strength
    endfunction get_L_filter_strength
    pure integer(4) function get_dL_filter_value() result(f)
    implicit none
        f=dL_filter_value
    endfunction get_dL_filter_value
    
    subroutine set_dx_filter_strength(dx_filter_strength_)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::dx_filter_strength_
        
        if(dx_filter_strength_<0) call print_error("matrix_pencil_method.set_dx_filter_strength","dx_filter_strength_<0")
        
        dx_filter_strength=dx_filter_strength_
        counted=.false.
    endsubroutine set_dx_filter_strength
    subroutine set_L_filter_strength(L_filter_strength_)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::L_filter_strength_
        
        if(L_filter_strength_<0) call print_error("matrix_pencil_method.set_L_filter_strength","L_filter_strength_<0")
        
        L_filter_strength=L_filter_strength_
        counted=.false.
    endsubroutine set_L_filter_strength
    subroutine set_dL_filter_value(dL_filter_value_)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::dL_filter_value_
        
        if(dL_filter_value_<0) call print_error("matrix_pencil_method.set_dL_filter_value","dL_filter_value_<0")
        
        dL_filter_value=dL_filter_value_
        counted=.false.
    endsubroutine set_dL_filter_value
    
    
    
    subroutine get_result(res,res_size)
    use system,only:print_error
    implicit none
        complex(8),intent(out)::res(L_)
        integer(4),intent(out)::res_size
        
        if(.not.counted) &
            call print_error("matrix_pencil_method.get_result",&
            "result has not counted yet(use count_dispersion_numbers)")
        
        call get_result_with_lower_filters(res,res_size,dx_filter_strength,L_filter_strength)
    endsubroutine get_result
    subroutine get_dx_filter_result(res,res_size,filter_number)
    use system,only:print_error
    implicit none
        complex(8),intent(out)::res(L_)
        integer(4),intent(out)::res_size
        integer(4),intent(in)::filter_number
        
        integer(4) i
        
        if(.not.counted) &
            call print_error("matrix_pencil_method.get_dx_filter_result",&
            "result has not counted yet(use count_dispersion_numbers)")
        if(.not.(1.le.filter_number.and.filter_number.le.dx_filter_strength)) &
            call print_error("matrix_pencil_method.get_dx_filter_result",&
            ".not.(1.le.filter_number.and.filter_number.le.dx_filter_strength)")
        
        res_size=dx_filter_res_size(filter_number)
        do i=1,res_size
            res(i)=dx_filter_res(filter_number,i)
        enddo
    endsubroutine get_dx_filter_result
    subroutine get_L_filter_result(res,res_size,filter_number)
    use system,only:print_error
    implicit none
        complex(8),intent(out)::res(L_)
        integer(4),intent(out)::res_size
        integer(4),intent(in)::filter_number
        
        integer(4) i
        
        if(.not.counted) &
            call print_error("matrix_pencil_method.get_dL_filter_result",&
            "result has not counted yet(use count_dispersion_numbers)")
        if(.not.(1.le.filter_number.and.filter_number.le.L_filter_strength)) &
            call print_error("matrix_pencil_method.get_dL_filter_result",&
            ".not.(1.le.filter_number.and.filter_number.le.L_filter_strength)")
        
        res_size=L_filter_res_size(filter_number)
        do i=1,res_size
            res(i)=L_filter_res(filter_number,i)
        enddo
    endsubroutine get_L_filter_result
    subroutine get_result_with_lower_filters(res,res_size,dx_filter_strength_,dL_filter_strength_)
    use system,only:print_error
    implicit none
        complex(8),intent(inout)::res(L_)
        integer(4),intent(out)::res_size
        integer(4),intent(in)::dx_filter_strength_
        integer(4),intent(in)::dL_filter_strength_
        
        if(.not.counted) &
            call print_error("matrix_pencil_method.get_result_with_lower_filters",&
            "result has not counted yet(use count_dispersion_numbers)")
        if(.not.(0.le.dx_filter_strength_.and.dx_filter_strength_.le.dx_filter_strength)) &
            call print_error("matrix_pencil_method.get_result_with_lower_filters",&
            ".not.(0.le.dx_filter_strength_.and.dx_filter_strength_.le.dx_filter_strength)")
        if(.not.(0.le.dL_filter_strength_.and.dL_filter_strength_.le.L_filter_strength)) &
            call print_error("matrix_pencil_method.get_result_with_lower_filters",&
            ".not.(0.le.dL_filter_strength_.and.dL_filter_strength_.le.L_filter_strength)")
        
        call count_result_with_filters(res,res_size,dx_filter_strength_,dL_filter_strength_)
    endsubroutine get_result_with_lower_filters
    
    
    
    subroutine init_matrix_pencil_method()
    implicit none
        allocate(res_(1))
        res_size_=-1
        
        allocate(dx_filter_res(1,1))
        allocate(dx_filter_res_size(1))
        
        allocate(L_filter_res(1,1))
        allocate(L_filter_res_size(1))
        
        counted=.false.
    endsubroutine init_matrix_pencil_method
    subroutine destructor_matrix_pencil_method
    implicit none
       deallocate(res_)
       
       deallocate(dx_filter_res)
       deallocate(dx_filter_res_size)
       
       deallocate(L_filter_res)
       deallocate(L_filter_res_size)
    endsubroutine destructor_matrix_pencil_method
endmodule matrix_pencil_method
