module runtime_application_actions
implicit none

    public::test_K,load_experimental_measurements_view
    contains
    
    subroutine load_experimental_measurements_view
    use load_experimental_measurements,only:get_Nx,get_Nt,get_xi,get_tj,get_uij
    implicit none
        integer(4) file
        
        integer(4) xi,tj
        
        real(8) omega,omega_start,domega,omega_end
        
        open(newunit=file,file="images\load_experimental_measurements\xi=1_signal.data")
        do tj=1,get_Nt()
            write(file,*),get_tj(tj),get_uij(1,tj)
        enddo
        close(file)
        open(newunit=file,file="images\load_experimental_measurements\xi=Nx_div_2_signal.data")
        do tj=1,get_Nt()
            write(file,*),get_tj(tj),get_uij(get_Nx()/2,tj)
        enddo
        close(file)
        open(newunit=file,file="images\load_experimental_measurements\xi=Nx_signal.data")
        do tj=1,get_Nt()
            write(file,*),get_tj(tj),get_uij(get_Nx(),tj)
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
