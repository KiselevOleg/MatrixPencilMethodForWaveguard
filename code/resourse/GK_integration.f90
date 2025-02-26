module GK_integration
    public::GK_integral_ab
    
    integer(4),save::number_of_integral=1
    
    private::number_of_integral
contains
    complex(8) function GK_integral_ab(functionName,accurate,a,b,upperPolesValue,depthOfAvoidingPoles) result(f)
    use GK_integration1,only:GK_integral_ab1=>GK_integral_ab
    use GK_integration2,only:GK_integral_ab2=>GK_integral_ab
    use GK_integration3,only:GK_integral_ab3=>GK_integral_ab
    use GK_integration4,only:GK_integral_ab4=>GK_integral_ab
    use GK_integration5,only:GK_integral_ab5=>GK_integral_ab
    use GK_integration6,only:GK_integral_ab6=>GK_integral_ab
    use GK_integration7,only:GK_integral_ab7=>GK_integral_ab
    implicit none
        complex(8),external::functionName
        real(8),intent(in)::accurate
        real(8),intent(in)::a
        real(8),intent(in)::b
        real(8),intent(in)::upperPolesValue
        real(8),intent(in)::depthOfAvoidingPoles
        
        if(number_of_integral==1) then
            number_of_integral=number_of_integral+1
            f=GK_integral_ab1(functionName,accurate,a,b,upperPolesValue,depthOfAvoidingPoles)
            number_of_integral=number_of_integral-1
            return
        endif
        if(number_of_integral==2) then
            number_of_integral=number_of_integral+1
            f=GK_integral_ab2(functionName,accurate,a,b,upperPolesValue,depthOfAvoidingPoles)
            number_of_integral=number_of_integral-1
            return
        endif
        if(number_of_integral==3) then
            number_of_integral=number_of_integral+1
            f=GK_integral_ab3(functionName,accurate,a,b,upperPolesValue,depthOfAvoidingPoles)
            number_of_integral=number_of_integral-1
            return
        endif
        if(number_of_integral==4) then
            number_of_integral=number_of_integral+1
            f=GK_integral_ab4(functionName,accurate,a,b,upperPolesValue,depthOfAvoidingPoles)
            number_of_integral=number_of_integral-1
            return
        endif
        if(number_of_integral==5) then
            number_of_integral=number_of_integral+1
            f=GK_integral_ab5(functionName,accurate,a,b,upperPolesValue,depthOfAvoidingPoles)
            number_of_integral=number_of_integral-1
            return
        endif
        if(number_of_integral==6) then
            number_of_integral=number_of_integral+1
            f=GK_integral_ab6(functionName,accurate,a,b,upperPolesValue,depthOfAvoidingPoles)
            number_of_integral=number_of_integral-1
            return
        endif
        if(number_of_integral==7) then
            number_of_integral=number_of_integral+1
            f=GK_integral_ab7(functionName,accurate,a,b,upperPolesValue,depthOfAvoidingPoles)
            number_of_integral=number_of_integral-1
            return
        endif
        
        pause "error: number of counting integral in the same time is more then 7"
        stop
    endfunction GK_integral_ab
    
    complex(8) function GK_integral(functionName,accurate,lengthOfIntegration,upperPolesValue,depthOfAvoidingPoles) result(f)
    use GK_integration1,only:GK_integral1=>GK_integral
    use GK_integration2,only:GK_integral2=>GK_integral
    use GK_integration3,only:GK_integral3=>GK_integral
    use GK_integration4,only:GK_integral4=>GK_integral
    use GK_integration5,only:GK_integral5=>GK_integral
    use GK_integration6,only:GK_integral6=>GK_integral
    use GK_integration7,only:GK_integral7=>GK_integral
    implicit none
        complex(8),external::functionName
        real(8),intent(in)::accurate
        real(8),intent(in)::lengthOfIntegration
        real(8),intent(in)::upperPolesValue
        real(8),intent(in)::depthOfAvoidingPoles
        
        if(number_of_integral==1) then
            number_of_integral=number_of_integral+1
            f=GK_integral1(functionName,accurate,lengthOfIntegration,upperPolesValue,depthOfAvoidingPoles)
            number_of_integral=number_of_integral-1
            return
        endif
        if(number_of_integral==2) then
            number_of_integral=number_of_integral+1
            f=GK_integral2(functionName,accurate,lengthOfIntegration,upperPolesValue,depthOfAvoidingPoles)
            number_of_integral=number_of_integral-1
            return
        endif
        if(number_of_integral==3) then
            number_of_integral=number_of_integral+1
            f=GK_integral3(functionName,accurate,lengthOfIntegration,upperPolesValue,depthOfAvoidingPoles)
            number_of_integral=number_of_integral-1
            return
        endif
        if(number_of_integral==4) then
            number_of_integral=number_of_integral+1
            f=GK_integral4(functionName,accurate,lengthOfIntegration,upperPolesValue,depthOfAvoidingPoles)
            number_of_integral=number_of_integral-1
            return
        endif
        if(number_of_integral==5) then
            number_of_integral=number_of_integral+1
            f=GK_integral5(functionName,accurate,lengthOfIntegration,upperPolesValue,depthOfAvoidingPoles)
            number_of_integral=number_of_integral-1
            return
        endif
        if(number_of_integral==6) then
            number_of_integral=number_of_integral+1
            f=GK_integral6(functionName,accurate,lengthOfIntegration,upperPolesValue,depthOfAvoidingPoles)
            number_of_integral=number_of_integral-1
            return
        endif
        if(number_of_integral==7) then
            number_of_integral=number_of_integral+1
            f=GK_integral7(functionName,accurate,lengthOfIntegration,upperPolesValue,depthOfAvoidingPoles)
            number_of_integral=number_of_integral-1
            return
        endif
        
        pause "error: number of counting integral in the same time is more then 7"
        stop
    endfunction GK_integral
endmodule  GK_integration
