module runtime_application_actions
implicit none
    public::test_K,load_experimental_measurements_view,count_matrix_pencil_method_distersion_curve_graphics
    
    contains
    
    subroutine count_matrix_pencil_method_distersion_curve_graphics()
    use matrix_pencil_method,only:count_dispersion_numbers
    use math,only:c0,pi
    implicit none
        integer(4) L,res_size
        complex(8),allocatable::res(:)
        
        real(8) omega,omega_start,domega,omega_end
        
        integer(4) i
        integer(4) file
        
        L=5
        allocate(res(L))
        
        omega_start=0.01; domega=0.025d0; omega_end=6.25d0*4
        
        open(newunit=file,file="graphics/matrix_pencil_method/dispersion_curve.data")
        do omega=omega_start,omega_end,domega
            call count_dispersion_numbers(omega+c0,L,res,res_size)
            
            do i=1,res_size
                write(file,*),omega*0.5d0/pi,real(res(i)),aimag(res(i))
            enddo
            
            !print*,omega,res_size
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
