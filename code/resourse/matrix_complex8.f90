module matrix_complex8
implicit none
    public::det,inverse,transform,multiply,pseudoinverse
    
    private::inverse_notsafea
    contains
    
    complex(8) function det(A,n) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::n
        complex(8),intent(inout)::A(n,n)
        
        integer(4) i,j,k,sign
        complex(8) t
        
        if(n<1) call print_error("matrix_complex8.det","n<1")
        
        sign=0
        
        do i=1,n
            k=i
            do while(A(k,i)==0d0)
                k=k+1
                if(k==n+1) exit
            enddo
            if(k==n+1) then 
                f=0d0
                return
            endif
            if(k/=i) then
                sign=sign+1
                
                do j=1,n
                    t=A(j,i)
                    A(j,i)=A(k,i)
                    A(k,i)=t
                enddo
            endif
            
            do j=i+1,n
                t=-A(j,i)/A(i,i)
                do k=i,n
                    A(j,k)=A(j,k)+A(i,k)*t
                enddo
            enddo
        enddo
        
        f=1d0
        do i=1,n
            f=f*A(i,i)
        enddo
        if(mod(sign,2)/=0) f=-f
    endfunction det
    
    subroutine inverse(n,A)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::n
        complex(8),intent(inout)::A(n,n)
        
        complex(8) A1(n,n)
        
        integer(4) i,j
        
        if(n<1) call print_error("matrix_complex8.inverse","n<1")
        
        call inverse_notsafea(n,A,A1)
        
        do i=1,n
            do j=1,n
                A(i,j)=A1(i,j)
            enddo
        enddo
    endsubroutine inverse
    subroutine inverse_notsafea(n,a,a_1)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::n
        complex(8),intent(inout)::a(n,n)
        complex(8),intent(out)::a_1(n,n)
        
        integer(4) i,j,k,i_
        complex(8) t
        
        do i=1,n
            do j=1,n
                a_1(i,j)=0d0
            enddo
            a_1(i,i)=1d0
        enddo
        
        do i=1,n
            k=i
            do while(abs(a(k,i))<=1d-9)
                k=k+1
                if(k==n+1) call print_error("matrix_complex8.inverse","detA=0")
            enddo
            if(k/=i) then
                do j=1,n
                    t=a(i,j)
                    a(i,j)=a(k,j)
                    a(k,j)=t
                    
                    t=a_1(i,j)
                    a_1(i,j)=a_1(k,j)
                    a_1(k,j)=t
                enddo
            endif
            
            do j=i+1,n
                t=-a(j,i)/a(i,i)
                do k=1,n
                    a(j,k)=a(j,k)+a(i,k)*t
                    a_1(j,k)=a_1(j,k)+a_1(i,k)*t
                enddo
            enddo
        enddo
        
        do i=1,n
            t=a(i,i)
            do j=1,n
                a(i,j)=a(i,j)/t
                a_1(i,j)=a_1(i,j)/t
            enddo
        enddo
        
        do i_=-n,-1
            i=-i_
            do j=1,i-1
                t=a(j,i)
                do k=1,n
                    a_1(j,k)=a_1(j,k)-t*a_1(i,k)
                    a(j,k)=a(j,k)-t*a(i,k)
                enddo
            enddo
        enddo
    endsubroutine inverse_notsafea
    
    subroutine transform(n,m,a,a_t)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::n,m
        complex(8),intent(in)::a(n,m)
        complex(8),intent(out)::a_t(m,n)
        
        integer(4) i,j
        
        if(n<1) call print_error("matrix_complex8.transform","n<1")
        if(m<1) call print_error("matrix_complex8.transform","m<1")
        
        do i=1,m
            do j=1,n
                a_t(i,j)=a(j,i)
            enddo
        enddo
    endsubroutine transform
    
    subroutine multiply(n,m,k,A,B,AB)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::n
        integer(4),intent(in)::m
        integer(4),intent(in)::k
        complex(8),intent(in)::A(n,m)
        complex(8),intent(in)::B(m,k)
        complex(8),intent(out)::AB(n,k)
        
        integer(4) i,j,l
        
        if(n<1) call print_error("matrix_complex8.multiply","n<1")
        if(m<1) call print_error("matrix_complex8.multiply","m<1")
        if(k<1) call print_error("matrix_complex8.multiply","k<1")
        
        do i=1,n
            do j=1,k
                AB(i,j)=(0d0,0d0)
                do l=1,m
                    AB(i,j)=AB(i,j)+A(i,l)*B(l,j)
                enddo
            enddo
        enddo
    endsubroutine 
    
    subroutine pseudoinverse(n,m,a,a_plus)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::n
        integer(4),intent(in)::m
        complex(8),intent(in)::a(n,m)
        complex(8),intent(out)::a_plus(m,n)
        
        complex(8) a_star(m,n),a_stara(m,m),a_stara_plus_deltaI_inverse(m,m)
        real(8) delta
        integer(4) i,j,k
        
        if(n<1) call print_error("matrix_complex8.pseudoinverse","n<1")
        if(m<1) call print_error("matrix_complex8.pseudoinverse","m<1")
        
        delta=1d-9
        
        call transform(n,m,a,a_star)
        do i=1,m
            do j=1,n
                a_star(i,j)=real(a_star(i,j))-(0d0,1d0)*imag(a_star(i,j))
            enddo
        enddo
        
        do i=1,m
            do j=1,m
                a_stara(i,j)=0d0
                do k=1,n
                    a_stara(i,j)=a_stara(i,j)+a_star(i,k)*a(k,j)
                enddo
            enddo
        enddo
        
        do i=1,m
            a_stara(i,i)=a_stara(i,i)+delta
        enddo
        call inverse_notsafea(m,a_stara,a_stara_plus_deltaI_inverse)
        
        do i=1,m
            do j=1,n
                a_plus(i,j)=0d0
                do k=1,m
                    a_plus(i,j)=a_plus(i,j)+a_stara_plus_deltaI_inverse(i,k)*a_star(k,j)
                enddo
            enddo
        enddo
    endsubroutine pseudoinverse
endmodule matrix_complex8
