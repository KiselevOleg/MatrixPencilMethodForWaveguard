module Jn
implicit none
    public::J0
    
    contains
    
    complex(8) function J0(x) result(f)
    use GK_integration,only:GK_integral_ab
    use math,only:pi
    use system,only:print_error
    implicit none
        complex(8),intent(in)::x
        
        if(real(x)<0d0) call print_error("Jn.J0","(real(x)<0d0")
        
        !if(abs(x)>4d0) then
        if(abs(x)>10d0) then
            f=sqrt(2d0/pi/x)*cos(x-pi*0.25d0)
        else
            !complex(8) function GK_integral_ab(functionName,accurate,a,b,upperPolesValue,depthOfAvoidingPoles)
            f=GK_integral_ab(inner_fun,1d-6,-pi,pi,0d0,0d0)
            f=f*0.5d0/pi
        endif
    contains
        complex(8) function inner_fun(tau) result(f)
        use math,only:ci
        implicit none
            complex(8),intent(in)::tau
            
            f=exp(-ci*x*sin(tau))
        endfunction inner_fun
    endfunction J0
    
    complex(8) function J1(x) result(f)
    use GK_integration,only:GK_integral_ab
    use math,only:pi
    use system,only:print_error
    implicit none
        complex(8),intent(in)::x
        
        if(real(x)<0d0) call print_error("Jn.J2","(real(x)<0d0")
        
        if(abs(x)<0.05d0) then
        !if(abs(x)<0.01d0) then
            f=x*0.5d0
        elseif(abs(x)>10d0) then
        !elseif(abs(x)>20d0) then
            f=sqrt(2d0/pi/x)*cos(x-pi*0.75d0)
        else
            !complex(8) function GK_integral_ab(functionName,accurate,a,b,upperPolesValue,depthOfAvoidingPoles)
            f=GK_integral_ab(inner_fun,1d-6,-pi,pi,0d0,0d0)
            f=f*0.5d0/pi
        endif
    contains
        complex(8) function inner_fun(tau) result(f)
        use math,only:ci
        implicit none
            complex(8),intent(in)::tau
            
            f=exp(-ci*x*sin(tau)+ci*1d0*tau)
        endfunction inner_fun
    endfunction J1
    
    complex(8) function J2(x) result(f)
    use GK_integration,only:GK_integral_ab
    use math,only:pi
    use system,only:print_error
    implicit none
        complex(8),intent(in)::x
        
        if(real(x)<0d0) call print_error("Jn.J2","(real(x)<0d0")
        
        if(abs(x)<0.05d0) then
        !if(abs(x)<0.01d0) then
            f=x*x*0.125d0
        elseif(abs(x)>40d0) then
        !elseif(abs(x)>60d0) then
            f=sqrt(2d0/pi/x)*cos(x-pi*1.25d0)
        elseif(abs(x)<15d0) then
        !elseif(abs(x)<30d0) then
            f=GK_integral_ab(inner_fun,1d-6,-pi,pi,0d0,0d0)
            f=f*0.5d0/pi
        else
            f=2d0/x*J1(x)-J0(x)
        endif
        !if(abs(x)<0.5d0) then
        !    f=x*x*0.125d0
        !elseif(abs(x)>35d0) then
        !    f=sqrt(2d0/pi/x)*cos(x-pi*1.25d0)
        !else
        !    !complex(8) function GK_integral_ab(functionName,accurate,a,b,upperPolesValue,depthOfAvoidingPoles)
        !    f=GK_integral_ab(inner_fun,1d-6,-pi,pi,0d0,0d0)
        !    f=f*0.5d0/pi
        !endif
    contains
        complex(8) function inner_fun(tau) result(f)
        use math,only:ci
        implicit none
            complex(8),intent(in)::tau
            
            f=exp(-ci*x*sin(tau)+ci*2d0*tau)
        endfunction inner_fun
    endfunction J2
endmodule Jn
