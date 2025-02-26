module main_parameters
implicit none
    public::number_of_layers,h,get_full_h,Calphabeta,rho,omega,Q,Qomega,get_number_of_layer,Cijkl
    public::E,nu,lambda,mu,Cp,Cs
    
    public::get_anisotropic,set_anisotropic
    public::set_down_border_condition_type,get_down_border_condition_type,&
        get_down_border_condition_type_fixed_border,get_down_border_condition_type_free_border,&
        get_down_border_condition_type_halfspace
    
    integer(4),save::anisotropic=-1
    integer(4),save::down_border_condition_type=-1
    
    integer(4),save::number_of_layers
    real(8),save,allocatable::h(:)!number_of_layer
    real(8),save,allocatable::rho(:)!number_of_layer
    
    complex(8),save,allocatable::Calphabeta(:,:,:)!alpha,beta,number_of_layer
    
    real(8),save,allocatable::E(:)!number_of_layer
    real(8),save,allocatable::nu(:)!number_of_layer
    real(8),save,allocatable::lambda(:)!number_of_layer
    real(8),save,allocatable::mu(:)!number_of_layer
    real(8),save,allocatable::Cp(:)!number_of_layer
    real(8),save,allocatable::Cs(:)!number_of_layer
    
    real(8),save::omega
    
    real(8),save::full_h=-1d0
    private::full_h,anisotropic,down_border_condition_type
    contains
    
    complex(8) function Q(ind,alpha,beta) result(f)
    use math,only:pi
    use Jn,only:J2
    use system,only:print_error
    implicit none
        integer(4),intent(in)::ind
        complex(8),intent(in)::alpha
        complex(8),intent(in)::beta
        
        complex(8) alpha_
        alpha_=sqrt(alpha*alpha+beta*beta)
        
        if(ind<1.or.ind>3) call print_error("main_parameters.Q","ind<1.or.ind>3")
        
        if(ind==1) then
            f=0d0
        elseif(ind==2) then
            f=0d0
        else
            f=1d0
        endif
    endfunction Q
    complex(8) function Qomega(omega) result(f)
    use math,only:pi
    use system,only:print_error
    implicit none
        complex(8),intent(in)::omega
        
        f=1d0
    endfunction Qomega
    
    real(8) function get_full_h() result(f)
    implicit none
        integer(4) i
        
        if(full_h>0d0) then
            f=full_h
            return
        endif
        
        full_h=0d0
        do i=1,number_of_layers
            full_h=full_h+h(i)
        enddo
        
        f=full_h
        
        if(get_down_border_condition_type()==get_down_border_condition_type_halfspace()) then
            f=1d4
        endif
    endfunction get_full_h
    
    integer(4) function get_number_of_layer(z) result(f)
    use system,only:print_error
    implicit none
        real(8),intent(in)::z
        
        integer(4) n
        real(8) zr
        
        integer(4) i
        
        if(z>0d0) call print_error("main_parameters.get_number_of_layer","z>0d0")
        
        zr=z
        do i=1,number_of_layers
            zr=zr+h(i)
            
            if(zr>=0d0) then
                f=i
                return
            endif
        enddo
        
        if(get_down_border_condition_type()==get_down_border_condition_type_halfspace()) then
            f=number_of_layers
        else
            call print_error("main_parameters.get_number_of_layer","z<-\sum\limits_{n=1}^{number_of_layer}h(n)")
        endif
    endfunction get_number_of_layer
    
    complex(8) function Cijkl(i,j,k,l,number_of_layer) result(f)
    use system,only:print_error
    implicit none
        integer(4),intent(in)::i
        integer(4),intent(in)::j
        integer(4),intent(in)::k
        integer(4),intent(in)::l
        integer(4),intent(in)::number_of_layer
        
        integer(4) alpha,beta
        
        if(.not.anisotropic==1) call print_error("main_parameters.Cijkl","not anisotropic materials")
        
        if(i<1.or.i>3) call print_error("main_parameters.Cijkl","i<1.or.i>3")
        if(j<1.or.j>3) call print_error("main_parameters.Cijkl","j<1.or.j>3")
        if(k<1.or.k>3) call print_error("main_parameters.Cijkl","k<1.or.k>3")
        if(l<1.or.l>3) call print_error("main_parameters.Cijkl","l<1.or.l>3")
        if(number_of_layer<1.and.number_of_layer>number_of_layers) then
            call print_error("main_parameters.Cijkl","number_of_layer<1.and.number_of_layer>number_of_layers")
        endif
        
        if(i==1.and.j==1) then
            alpha=1
        elseif(i==2.and.j==2) then
            alpha=2
        elseif(i==3.and.j==3) then
            alpha=3
        elseif(i==2.and.j==3.or.i==3.and.j==2) then
            alpha=4
        elseif(i==1.and.j==3.or.i==3.and.j==1) then
            alpha=5
        elseif(i==1.and.j==2.or.i==2.and.j==1) then
            alpha=6
        endif
        
        if(k==1.and.l==1) then
            beta=1
        elseif(k==2.and.l==2) then
            beta=2
        elseif(k==3.and.l==3) then
            beta=3
        elseif(k==2.and.l==3.or.k==3.and.l==2) then
            beta=4
        elseif(k==1.and.l==3.or.k==3.and.l==1) then
            beta=5
        elseif(k==1.and.l==2.or.k==2.and.l==1) then
            beta=6
        endif
        
        f=Calphabeta(alpha,beta,number_of_layer)
    endfunction Cijkl
    
    integer(4) function get_anisotropic() result(f)
    use system,only:print_error,print_warning
    implicit none
        if(anisotropic==-1) call print_error("main_parameters.get_anisotropic","anisotropic type is not defenited")
        f=anisotropic
    endfunction get_anisotropic
    subroutine set_anisotropic(anisotropic_v)
    use system,only:print_error,print_warning
    implicit none
        integer(4),intent(in)::anisotropic_v
        
        !if(.not.anisotropic==-1) call print_error("main_parameters.set_anisotropic","anisotropic type is already defenited")
        if(.not.anisotropic_v==0.and..not.anisotropic_v==1) call print_error("main_parameters.set_anisotropic",&
            ".not.anisotropic==0.and.not.anisotropic==1")
        
        anisotropic=anisotropic_v
    endsubroutine set_anisotropic
    
    integer(4) function get_down_border_condition_type_fixed_border() result(f)
    implicit none
        f=1
    endfunction get_down_border_condition_type_fixed_border
    integer(4) function get_down_border_condition_type_free_border() result(f)
    implicit none
        f=2
    endfunction get_down_border_condition_type_free_border
    integer(4) function get_down_border_condition_type_halfspace() result(f)
    implicit none
        f=3
    endfunction get_down_border_condition_type_halfspace
    integer(4) function get_down_border_condition_type() result(f)
    use system,only:print_error,print_warning
    implicit none
        if(down_border_condition_type==-1) &
            call print_error("main_parameters.get_down_border_condition_type","down_border_condition_type status is not defenited")
        f=down_border_condition_type
    endfunction get_down_border_condition_type
    subroutine set_down_border_condition_type(down_border_condition_type_v)
    use system,only:print_error,print_warning
    implicit none
        integer(4),intent(in)::down_border_condition_type_v
        
        !if(.not.down_border_condition_type==-1) &
        !    call print_error("main_parameters.set_down_border_condition_type","down_border_condition_type status is already defenited")
        if(.not.down_border_condition_type_v==get_down_border_condition_type_fixed_border()&
            .and..not.down_border_condition_type_v==get_down_border_condition_type_free_border()&
            .and..not.down_border_condition_type_v==get_down_border_condition_type_halfspace()) &
            call print_error("main_parameters.set_down_border_condition_type",&
            "down_border_condition_type is incorrect")
        
        down_border_condition_type=down_border_condition_type_v
    endsubroutine set_down_border_condition_type
endmodule main_parameters
