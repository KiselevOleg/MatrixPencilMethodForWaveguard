module runtime_application_actions
implicit none
    public::test_K,load_experimental_measurements_view,count_matrix_pencil_method_distersion_curve_graphics,&
        count_distersion_curve_for_K_graphics,count_complex_distersion_curve_for_K_graphics,&
        count_complex_det_A_matrix_in_K_graphics
    
    contains
    
    subroutine count_complex_det_A_matrix_in_K_graphics()
    use math,only:ci
    implicit none
        real(8) phi
        integer(4) number_of_en
        
        real(8) Re,Re_min,dRe,Re_max
        real(8) Im,Im_min,dIm,Im_max
        
        complex(8) v
        
        integer(4) file1,file2,file3,file4
        integer(4) i
        
        phi=0d0
        number_of_en=3
        
        Re_min=0d0; dRe=1d-3; Re_max=3.90001d0
        Im_min=-0.01d0; dIm=1d-2; Im_max=0.0101d0
        
        open(newunit=file1,file="graphics/complex_det_A_matrix_in_K/sizes.data")
        open(newunit=file2,file="graphics/complex_det_A_matrix_in_K/A.data")
        open(newunit=file3,file="graphics/complex_det_A_matrix_in_K/AY.data")
        open(newunit=file4,file="graphics/complex_det_A_matrix_in_K/AX.data")
        
        i=(Re_max-Re_min)/dRe
        write(file1,*),i+1
        i=(Im_max-Im_min)/dIm
        write(file1,*),i+1
        
        do Re=Re_min,Re_max,dRe
            do Im=Im_min,Im_max,dIm
                !v=in_f_with_A(Re+Im*ci)
                !write(file2,*),Re,Im,abs(v)
                
                v=in_f_with_AY(Re+Im*ci)
                write(file3,*),Re,Im,abs(v)
                
                v=in_f_with_AX(Re+Im*ci)
                write(file4,*),Re,Im,abs(v)
            enddo
            
            print*,Re
        enddo
        
        close(file4)
        close(file3)
        close(file2)
        close(file1)
        
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
    endsubroutine count_complex_det_A_matrix_in_K_graphics
    
    subroutine count_complex_distersion_curve_for_K_graphics()
    use main_parameters,only:omega
    use dispersion_curves_for_K,only:count_complex_poles
    use math,only:pi
    implicit none
        real(8) omega_start,domega,omega_end
        
        real(8) phi
        
        complex(8) res(100)
        integer(4) res_size_max
        integer(4) res_size
        
        integer(4) i
        
        integer(4) file
        
        res_size_max=100;
        omega_start=0.3d0; domega=0.005d0; omega_end=12.499d0
        
        phi=0d0
        
        open(newunit=file,file="graphics/complex_distersion_curve_for_K/dispersion_curves.data")
        do omega=omega_start,omega_end,domega
            call count_complex_poles(phi,3,res_size_max,res,res_size)
            
            do i=1,res_size
                write(file,*),omega/pi*0.5d0,real(res(i)),aimag(res(i))
            enddo
            
            print*,omega,res_size
        enddo
        close(file)
    endsubroutine count_complex_distersion_curve_for_K_graphics
    subroutine count_distersion_curve_for_K_graphics()
    use main_parameters,only:omega
    use dispersion_curves_for_K,only:count_poles
    use math,only:pi
    implicit none
        real(8) omega_start,domega,omega_end
        
        real(8) phi
        
        real(8) res(100)
        integer(4) res_size_max
        integer(4) res_size
        
        integer(4) i
        
        integer(4) file
        
        res_size_max=100;
        omega_start=1d-1; domega=0.05d0; omega_end=12.499d0
        
        phi=0d0
        
        open(newunit=file,file="graphics/distersion_curve_for_K/dispersion_curves.data")
        do omega=omega_start,omega_end,domega
            call count_poles(phi,3,res_size_max,res,res_size)
            
            do i=1,res_size
                write(file,*),omega/pi*0.5d0,res(i)
            enddo
            
            print*,omega,res_size
        enddo
        close(file)
    endsubroutine count_distersion_curve_for_K_graphics
    
    subroutine count_matrix_pencil_method_distersion_curve_graphics()
    use matrix_pencil_method,only:count_dispersion_numbers
    use math,only:c0,pi
    implicit none
        integer(4) L,res_size
        complex(8),allocatable::res(:)
        
        real(8) omega,omega_start,domega,omega_end
        
        integer(4) i
        integer(4) file
        
        L=90
        allocate(res(L))
        
        omega_start=0.01; domega=0.01d0; omega_end=6.25d0*2
        
        open(newunit=file,file="graphics/matrix_pencil_method/dispersion_curve.data")
        do omega=omega_start,omega_end,domega
            call count_dispersion_numbers(omega+c0,L,res,res_size)
            
            do i=1,res_size
                write(file,*),omega*0.5d0/pi,real(res(i)),aimag(res(i))
            enddo
            
            print*,omega,res_size
        enddo
        close(file)
        
        deallocate(res)
    endsubroutine count_matrix_pencil_method_distersion_curve_graphics
    
    subroutine load_experimental_measurements_view
    use load_experimental_measurements,only:get_Nx,get_Nt,get_xi,get_tj,get_uij,get_sceptum_in_xi
    use math,only:pi,c0
    implicit none
        integer(4) file
        
        integer(4) xi,tj
        
        real(8) omega,omega_start,domega,omega_end
        
        open(newunit=file,file="graphics\load_experimental_measurements\xi=1_signal.data")
        do tj=1,get_Nt()
            write(file,*),get_tj(tj),get_uij(1,tj)
        enddo
        close(file)
        open(newunit=file,file="graphics\load_experimental_measurements\xi=Nx_div_2_signal.data")
        do tj=1,get_Nt()
            write(file,*),get_tj(tj),get_uij(get_Nx()/2,tj)
        enddo
        close(file)
        open(newunit=file,file="graphics\load_experimental_measurements\xi=Nx_signal.data")
        do tj=1,get_Nt()
            write(file,*),get_tj(tj),get_uij(get_Nx(),tj)
        enddo
        close(file)
        
        omega_start=0.01d0; domega=0.0025d0; omega_end=6.25d0*4
        
        open(newunit=file,file="graphics\load_experimental_measurements\xi=1_spectrum.data")
        do omega=omega_start,omega_end,domega
            write(file,*),omega/pi*0.5d0,abs(get_sceptum_in_xi(1,omega+c0))
        enddo
        close(file)
        open(newunit=file,file="graphics\load_experimental_measurements\xi=Nx_div_2_spectrum.data")
        do omega=omega_start,omega_end,domega
            write(file,*),omega/pi*0.5d0,abs(get_sceptum_in_xi(get_Nx()/2,omega+c0))
        enddo
        close(file)
        open(newunit=file,file="graphics\load_experimental_measurements\xi=Nx_spectrum.data")
        do omega=omega_start,omega_end,domega
            write(file,*),omega/pi*0.5d0,abs(get_sceptum_in_xi(get_Nx(),omega+c0))
        enddo
        close(file)
    endsubroutine load_experimental_measurements_view
    
    subroutine test_K()
    use count_K,only:K
    implicit none
        complex(8) alpha,beta
        real(8) z
        
        integer(4) i,j
        
        alpha=(1d0,0d0)
        beta=(0.5,0d0)
        z=0d0
        
        do i=1,3
            do j=1,3
                print*,K(i=i,j=j,alpha=alpha,beta=beta,z=z)
            enddo
            print*
        enddo
    endsubroutine test_K
endmodule runtime_application_actions
