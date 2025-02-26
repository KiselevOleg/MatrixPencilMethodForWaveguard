!STAR5 is a solver for a complex linear algebraic system Ax=b 
!       by the double-orthogonalization method 
! N - system size; M - number of the right-hand parts;
! Nd - dimension (number of the lines) declared in fact; 
! A(Nd,N) - N*N complex matrix of the system; 
! B(Nd,M) - N*M complex matrix of the right-hand parts;
! C(Nd,N),R(Nd,M+3) - complex work arrays
! KORT - number of the ortogonalization cycles 
!(2 - 3 are usually quite enough)
 
module linear_system
implicit none
    public::star5_
    
    interface star5_
        module procedure star5_1
        module procedure star5_2
    endinterface star5_
    
    private::star5,star5_1,star5_2
    contains
    subroutine star5_saveA(A,B)
    implicit none
        complex(8),intent(in)::A(:,:)
        complex(8),intent(inout)::B(:,:)
        
        complex(8) A_(size(A,1),size(A,2))
        
        integer(4) i,j,n
        
        n=size(A,1)
        do i=1,n
            do j=1,n
                A_(i,j)=A(i,j)
            enddo
        enddo
        
        call star5_(A_,B)
    endsubroutine star5_saveA
    
    subroutine star5_2(A,B)
    implicit none
        !integer(4),intent(in)::N
        complex(8),intent(in)::A(:,:)
        complex(8),intent(inout)::B(:,:)
        
        if(.not.size(B(:,1))==size(A(1,:))) then
            pause "B has incorrect size"
            stop
        endif
        
        call star5_1(A,B,size(A(1,:)),size(B(1,:)),5)
    endsubroutine star5_2
    subroutine star5_1(A,B,N,M,kort)
    implicit none
        integer(4),intent(in)::kort,M,N
        complex(8),intent(in)::A(N,N)
        complex(8),intent(inout)::B(N,M)
        
        integer(4) ND
        complex(8),allocatable::C(:,:),R(:,:)
        
        ND=N
        allocate(C(ND,N))
        allocate(R(ND,M+3))
        call STAR5(A,B,C,R,ND,N,M,kort)
        deallocate(C)
        deallocate(R)
    endsubroutine star5_1
        
      SUBROUTINE STAR5(A,B,C,R,ND,N,M,kort)
      implicit none
!      implicit real*8(a-h,o-z)
      
	  integer Nd,N,M,kort,i,k,k1,kw,jr,im,Ms,Msm,Mt
	  
	  real*8 ra,rt

	  COMPLEX*16 A(Nd,N),B(Nd,M),C(Nd,N),R(Nd,M+3),TK,ST

	  Ms=M+1; Msm=M+2; Mt=M+3

! orthogonolization

      DO 1 K=1,N
        K1=K-1;  RA=0D0;  KW=0

        DO 2 I=1,N
          ST=A(I,K); 
		  R(I,Ms)=ST;  
		  R(I,Msm)=ST
  2     RA=RA+ST*CONJG(ST)
  !		print*,Ra
        RA=SQRT(RA)

        IF(K.NE.1)GOTO 12
        R(1,Mt)=RA

        DO 13 I=1,N
 13     C(I,1)=R(I,Ms)/RA
        GOTO 1

 12     DO 4 JR=1,K1
          ST=0D0
          DO 5 I=1,N
  5       ST=ST+R(I,Ms)*CONJG(C(I,JR))

          DO 7 I=1,N
  7       R(I,Msm)=R(I,Msm)-ST*C(I,JR)

  4     CONTINUE

        RT=0D0
        DO 14 I=1,N
 14     RT=RT+ABS(R(I,Msm))**2
        
!		print*, RT
        RT=SQRT(RT)
        IF(RT.GT.1D-13)GOTO 9
     
          DO 8 I=1,N
  8       R(I,Ms)=0D0

          ST=0D0
        GOTO 28
 
  9     KW=KW+1
        DO 10 I=1,N
          ST=R(I,Msm)/RT;  R(I,Ms)=ST
 10     R(I,Msm)=ST

        IF(KW.LT.KORT)GOTO 12
 
 34     ST=0D0
        DO 21 I=1,N
 21     ST=ST+R(I,Ms)*CONJG(A(I,K))
  
 28     R(K,Mt)=ST
 
        DO 23 I=1,N
 23     C(I,K)=R(I,Ms)

  1   CONTINUE
      
! backward steps
	
	K=N
      DO 44 IM=1,M
  3     DO 44 I=1,N
 44     R(I,IM)=B(I,IM)

 45     DO 49 I=1,N
 49     R(I,Ms)=C(I,K)

        TK=R(K,Mt)
        IF(ABS(TK).GT.1D-12)GOTO 50
          DO 55 IM=1,M
 55       B(K,IM)=0D0
        GOTO 56
 
 50     DO 51 IM=1,M
          ST=0D0
          DO 52 I=1,N
 52       ST=ST+R(I,Ms)*CONJG(R(I,IM))

 51     B(K,IM)=CONJG(ST/TK)
        
		IF(K.EQ.1) RETURN
      
	    DO 54 I=1,N
          TK=A(I,K)
          DO 54 IM=1,M
 54     R(I,IM)=R(I,IM)-B(K,IM)*TK

 56     K=K-1
        IF(K.GE.1)GOTO 45
      
	    RETURN
      END
      
endmodule linear_system
