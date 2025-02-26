! Lapac procedure for matrix eigenproblem

subroutine EV_evaluate(A,N,e_values,e_vectors)
    implicit none
        integer N,ilo,ihi,info,lwork,ldz,m
        real*8 scale_arr(N),rwork(2*N)
        logical select(N)
        complex*16 A(N,N),e_values(N),e_vectors(N,N),A1(N,N),tau(N-1),z(N,N),Q(N,N),t_work(1)
        complex*16 vl(N,N),vr(N,N)
        complex*16,allocatable:: work(:)
            A1=A
            lwork=2*N
            call zgebal('B',N,A1,N,ilo,ihi,scale_arr,info)
                call zgehrd(N,ilo,ihi,A1,N,tau,t_work,-1,info)
                lwork=nint(abs(t_work(1))); allocate(work(lwork))
            call zgehrd(N,ilo,ihi,A1,N,tau,work,lwork,info)
                deallocate(work)
            Q=A1
                call zunghr(N,ilo,ihi,Q,N,tau,t_work,-1,info)
                lwork=nint(abs(t_work(1))); allocate(work(lwork))
            call zunghr(N,ilo,ihi,Q,N,tau,work,lwork,info)
                deallocate(work)
            ldz=N
                call zhseqr('S','V', N,ilo,ihi,A1,N,e_values,Q,ldz,t_work,-1,info)
                lwork=nint(abs(t_work(1))); allocate(work(lwork))
            call zhseqr('S','V', N,ilo,ihi,A1,N,e_values,Q,ldz,work,lwork,info)
                deallocate(work)
                allocate(work(2*N))
            call ztrevc('R','B',select,N,A1,N,vl,N,Q,N,N,m,work,rwork,info)
                deallocate(work)
            vr=Q
            call zgebak('B','R',N,ilo,ihi,scale_arr,N,vr,N,info)
            e_vectors=vr
end subroutine EV_evaluate