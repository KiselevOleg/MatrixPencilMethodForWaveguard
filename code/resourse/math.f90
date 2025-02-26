module math
implicit none
    public::pi,e,ci,c0,epsilon,&
        from_cartesian2_to_polar_coordinate_system,from_polar_to_cartesian2_coordinate_system
    
    real(8),parameter::pi=3.14159265358979d0
    real(8),parameter::e=2.718281828459045d0
    complex(8),parameter::ci=(0d0,1d0)
    complex(8),parameter::c0=(0d0,0d0)
    real(8),parameter::epsilon=1d-7
    
    interface from_cartesian2_to_polar_coordinate_system
        procedure from_cartesian2_to_polar_coordinate_system_real
        procedure from_cartesian2_to_polar_coordinate_system_complex
    endinterface from_cartesian2_to_polar_coordinate_system
    interface from_polar_to_cartesian2_coordinate_system
        procedure from_polar_to_cartesian2_coordinate_system_real
        procedure from_polar_to_cartesian2_coordinate_system_complex
    endinterface from_polar_to_cartesian2_coordinate_system
    
    private::from_cartesian2_to_polar_coordinate_system_real,from_cartesian2_to_polar_coordinate_system_complex,&
        from_polar_to_cartesian2_coordinate_system_real,from_polar_to_cartesian2_coordinate_system_complex
    contains
    
    subroutine from_cartesian2_to_polar_coordinate_system_real(x,y,r,phi)
    implicit none
        real(8),intent(in)::x
        real(8),intent(in)::y
        real(8),intent(out)::r
        real(8),intent(out)::phi
        
        r=sqrt(x*x+y*y)
        
        if(abs(x)<epsilon) then
            if(y>=0d0) then
                phi=pi*0.5d0
            else
                phi=pi*1.5d0
            endif
        else if(abs(y)<epsilon) then
            if(x>=0) then
                phi=0d0
            else 
                phi=pi
            endif
        elseif(x>=0d0.and.y>=0d0) then
            phi=atan(y/x)
        elseif(x>=0d0.and.y<=0d0) then
            phi=pi+pi-atan(-y/x)
        elseif(x<=0d0.and.y>=0d0) then
            phi=pi-atan(-y/x)
        else
            phi=pi+atan(y/x)
        endif
    endsubroutine from_cartesian2_to_polar_coordinate_system_real
    subroutine from_cartesian2_to_polar_coordinate_system_complex(alpha1,alpha2,alpha,gamma)
    implicit none
        complex(8),intent(in)::alpha1
        complex(8),intent(in)::alpha2
        complex(8),intent(out)::alpha
        complex(8),intent(out)::gamma
        
        alpha=sqrt(alpha1*alpha1+alpha2*alpha2)
        
        if(abs(real(alpha1))<epsilon) then
            if(real(alpha2)>=0d0) then
                gamma=pi*0.5d0
            else
                gamma=pi*1.5d0
            endif
        else if(abs(real(alpha2))<epsilon) then
            if(real(alpha1)>=0) then
                gamma=0d0
            else 
                gamma=pi
            endif
        elseif(real(alpha1)>=0d0.and.real(alpha2)>=0d0) then
            gamma=atan(alpha2/alpha1)
        elseif(real(alpha1)>=0d0.and.real(alpha2)<=0d0) then
            gamma=pi+pi-atan(-alpha2/alpha1)
        elseif(real(alpha1)<=0d0.and.real(alpha2)>=0d0) then
            gamma=pi-atan(-alpha2/alpha1)
        else
            gamma=pi+atan(alpha2/alpha1)
        endif
    endsubroutine from_cartesian2_to_polar_coordinate_system_complex
    subroutine from_polar_to_cartesian2_coordinate_system_real(x,y,r,phi)
    use system,only:print_error
    implicit none
        real(8),intent(out)::x
        real(8),intent(out)::y
        real(8),intent(in)::phi
        real(8),intent(in)::r
        
        if(r<0d0) call print_error("math.from_polar_to_cartesian2_coordinate_system_real","r<=0d0")
        if(phi<0d0.or.phi>=2d0*pi) call print_error("math.from_polar_to_cartesian2_coordinate_system_real","phi<0d0.or.phi>=2d0*pi")
        
        x=r*cos(phi)
        y=r*sin(phi)
    endsubroutine from_polar_to_cartesian2_coordinate_system_real
    subroutine from_polar_to_cartesian2_coordinate_system_complex(alpha1,alpha2,alpha,gamma)
    use system,only:print_error
    implicit none
        complex(8),intent(out)::alpha1
        complex(8),intent(out)::alpha2
        complex(8),intent(in)::alpha
        complex(8),intent(in)::gamma
        
        if(real(alpha)<0d0) call print_error("math.from_polar_to_cartesian2_coordinate_system_complex","real(alpha)<0d0")
        if(real(gamma)<0d0.or.real(gamma)>=2d0*pi) then
            call print_error("math.from_polar_to_cartesian2_coordinate_system_complex","real(gamma)<0d0.or.real(gamma)>=2d0*pi")
        endif
        
        alpha1=alpha*cos(gamma)
        alpha2=alpha*sin(gamma)
    endsubroutine from_polar_to_cartesian2_coordinate_system_complex
endmodule math
