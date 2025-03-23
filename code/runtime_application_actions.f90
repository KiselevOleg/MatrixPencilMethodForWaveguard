module runtime_application_actions
implicit none
    public::test_K,load_experimental_measurements_view,count_matrix_pencil_method_distersion_curve_graphics,&
        count_distersion_curve_for_K_graphics,count_complex_distersion_curve_for_K_graphics,&
        count_complex_det_A_matrix_in_K_graphics,count_wavelet_transform_smoothing,&
        find_material_properties_from_experimental_measurements
    
    contains
    
    subroutine find_material_properties_from_experimental_measurements()
    use material_properties_by_dispersion_curves,only:find_material_properties
    use math,only:pi
    implicit none
        integer(4) dispersion_curves_size
        real(8),allocatable::dispersion_curves(:,:)
        integer(4) parameters_for_detect_size
        integer(4),allocatable::parameters_layer(:)
        character(len=3),allocatable::parameters_type(:)
        real(8),allocatable::parameters_min(:),parameters_max(:),dparameters(:)
        integer(4) res_size,res_max_size
        real(8),allocatable::res(:,:)
        real(8),allocatable::res_value_of_right(:)
        
        integer(4) file
        integer(4) i,j
        
        !open(newunit=file,file="input/glass_experimental_dispersion_curves/_dispersion_curve_experimental.data")
        !open(newunit=file,file="input/glass_experimental_dispersion_curves/55000/with_first_points/_dispersion_curve_experimental.data")
        !open(newunit=file,file="input/glass_experimental_dispersion_curves/55000/with_first_points/_dispersion_curve_experimental_extended.data")
        open(newunit=file,file="input/glass_experimental_dispersion_curves/55000/with_first_points/accurate/_dispersion_curve_experimental.data")
        read(file,*),dispersion_curves_size
        allocate(dispersion_curves(dispersion_curves_size,2))
        do i=1,dispersion_curves_size
            read(file,*),dispersion_curves(i,1),dispersion_curves(i,2)
            dispersion_curves(i,1)=dispersion_curves(i,1)*2d0*pi
        enddo
        close(file)
        
        parameters_for_detect_size=2
        allocate(parameters_layer(parameters_for_detect_size))
        allocate(parameters_type(parameters_for_detect_size))
        allocate(parameters_min(parameters_for_detect_size))
        allocate(parameters_max(parameters_for_detect_size))
        allocate(dparameters(parameters_for_detect_size))
        
        parameters_layer(1)=1;  parameters_type(1)="E";   parameters_min(1)=0.5d0;    parameters_max(1)=1.00001d0;    dparameters(1)=0.25000d0/4
        parameters_layer(2)=1;  parameters_type(2)="nu";  parameters_min(2)=0.1d0;    parameters_max(2)=0.49999d0;    dparameters(2)=0.09999d0/4
        !parameters_layer(3)=1;  parameters_type(3)="h";   parameters_min(3)=0.25d0;   parameters_max(3)=0.30001d0;    dparameters(3)=0.05000d0/4
        !parameters_layer(4)=1;  parameters_type(4)="rho"; parameters_min(4)=2.30d0;   parameters_max(4)=2.50001d0;    dparameters(4)=0.10000d0/4
        !parameters_layer(1)=1;  parameters_type(1)="rho"; parameters_min(1)=2.30d0;   parameters_max(1)=2.50001d0;    dparameters(1)=0.10000d0/4
        
        res_max_size=1000
        allocate(res(res_max_size,parameters_for_detect_size))
        allocate(res_value_of_right(res_max_size))
        
        call find_material_properties(dispersion_curves_size,dispersion_curves,.true.,&
            parameters_for_detect_size,parameters_layer,parameters_type,parameters_min,parameters_max,dparameters,&
            res_max_size,res,res_value_of_right,res_size)
        
        print*,"res_size",res_size
        do i=1,min(20,res_size)
            print*
            print*,"res",i
            do j=1,parameters_for_detect_size
                print*,"layer",parameters_layer(j)
                print*,"type            ",parameters_type(j)
                print*,"value",res(i,j)
            enddo
            print*,"right",res_value_of_right(i)
            print*
        enddo
        
        deallocate(res_value_of_right)
        deallocate(res)
        deallocate(dparameters)
        deallocate(parameters_max)
        deallocate(parameters_min)
        deallocate(parameters_type)
        deallocate(parameters_layer)
    endsubroutine find_material_properties_from_experimental_measurements
    
    subroutine count_wavelet_transform_smoothing()
    use smoothing_signal,only:wavelet_transform_smoothing
    use load_experimental_measurements,only:get_Nx,get_Nt,get_xi,get_tj,get_uij
    implicit none
        real(8),allocatable::time(:),signal(:)
        
        real(8) a,a_min,da,a_max
        real(8) b,b_min,db,b_max
        
        real(8) t
        integer(4) xi,tj
        
        integer(4) xi_min,xi_max
        character(1024) filename
        integer(4) file
        integer(4) i,j
        
        xi_min=1; xi_max=get_Nx()
        
        read(*,*),xi_min,xi_max
        xi_max=min(xi_max,get_Nx())
        write(filename,*),"u_smoothing",xi_min,xi_max,".data"
        open(newunit=file,file=filename)
        !open(newunit=file,file="input/glass/u_smoothing.data")
        !open(newunit=file,file="graphics/load_experimental_measurements/u_smoothing.data")
        do xi=1,get_Nx()
            if(.not.(xi_min.le.xi.and.xi.le.xi_max)) cycle
            !if(xi/=get_Nx()/2) cycle
            print*,xi,get_Nx()
            
            allocate(time(get_Nt()))
            allocate(signal(get_Nt()))
            
            do tj=1,get_Nt()
                time(tj)=get_tj(tj)
                signal(tj)=get_uij(xi,tj)
            enddo
            
            a_min=1d-6; da=0.01d0; a_max=15.0001d0
            b_min=1d-6+30d0; db=0.5d0; b_max=100.0001d0
            
            a_min=1d-6; da=0.05d0; a_max=20.0001d0
            b_min=0d0; db=0.5d0; b_max=500.0001d0
            
            a_min=0.025d0; da=0.025d0; a_max=20.0001d0
            b_min=0d0; db=0.25d0; b_max=300.0001d0
            
            call wavelet_transform_smoothing(get_Nt(),time,signal,a_min,da,a_max,b_min,db,b_max)
            
            print*,"printing"
            do j=1,get_Nt()
                if(mod(j,2)==1) cycle
                t=get_tj(j)
                write(file,*),signal(j)
                !print*,t,j,get_Nt()
            enddo
            
            deallocate(time)
            deallocate(signal)
        enddo
        close(file)
        
        if(.not.xi_min==1) return
        
        !open(newunit=file,file="graphics/load_experimental_measurements/t_smoothing.data")
        open(newunit=file,file="t_smoothing.data")
        i=0
        do j=1,get_Nt()
            if(mod(j,2)==1) cycle
            i=i+1
        enddo
        write(file,*),i
        do j=1,get_Nt()
            if(mod(j,2)==1) cycle
            t=get_tj(j)
            write(file,*),t*1d-6
        enddo
        close(file)
    endsubroutine count_wavelet_transform_smoothing
    
    subroutine count_wavelet_transform_graphics()
    use smoothing_signal,only:wavelet_transform,mexican_hat,morlet_wavelet,meyer_wavelet
    use load_experimental_measurements,only:get_Nx,get_Nt,get_xi,get_tj,get_uij
    implicit none
        integer(4) xi,tj
        real(8),allocatable::time(:),signal(:)
        
        real(8) a,a_min,da,a_max
        real(8) b,b_min,db,b_max
        
        integer(4) file1,file2
        
        xi=get_Nx()/2
        allocate(time(get_Nt()))
        allocate(signal(get_Nt()))
        
        do tj=1,get_Nt()
            time(tj)=get_tj(tj)
            signal(tj)=get_uij(xi,tj)
        enddo
        
        open(newunit=file1,file="graphics/wavelet_transform/wavelet_transform_sizes.data")
        open(newunit=file2,file="graphics/wavelet_transform/wavelet_transform.data")
        
        a_min=1d-6; da=2d-1; a_max=40.0001d0
        b_min=1d-6; db=4d-1; b_max=300.0001d0
        
        !a_min=1d-6; da=0.0625d-1; a_max=0.25001d0
        !b_min=60d0; db=0.5d-1; b_max=100.0001d0
        
        write(file1,*),(a_max-a_min)/da+1
        write(file1,*),(b_max-b_min)/db+1
        
        do a=a_min,a_max,da
            do b=b_min,b_max,db
                write(file2,*),a,b,wavelet_transform(morlet_wavelet,a,b,get_Nt(),time,signal)
            enddo
            
            print*,a
        enddo
        
        close(file2)
        close(file1)
        
        deallocate(time)
        deallocate(signal)
    endsubroutine count_wavelet_transform_graphics
    
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
        omega_start=1d-1; domega=0.05d0; omega_end=12.499d0*1.5d0
        
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
        
        L=90/2
        L=40*2-20-20
        allocate(res(L))
        
        omega_start=0.01; domega=0.01d0*10d0; omega_end=6.25d0*2
        omega_start=0.1d0*20; domega=0.05d0*5*3/15; omega_end=6.25d0*3
        
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
