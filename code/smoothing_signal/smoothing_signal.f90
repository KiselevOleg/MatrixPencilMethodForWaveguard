module smoothing_signal
implicit none
    public::arithmetic_mean_smoothing
    
    contains
    
    subroutine arithmetic_mean_smoothing(strong, signal_size, signal)
    use system,only:print_error,print_warning
    implicit none
        integer(4),intent(in)::strong
        integer(4),intent(in)::signal_size
        real(8),intent(inout)::signal(signal_size)
        
        real(8) tail_sum_signal(strong)
        integer(4) tail_sum_signal_start
        real(8) sum,sum_arithmetic_mean,sum_del
        
        integer(4) i
        
        if(strong<1.or.strong>50) call print_error("smoothing_signal.arithmetic_mean_smoothing","strong<1.or.>strong>50")
        if(signal_size<=strong+strong) call print_error("smoothing_signal.arithmetic_mean_smoothing","signal_size<s=trong+strong")
        if(signal_size<50) call print_warning("smoothing_signal.arithmetic_mean_smoothing","signal_size<50")
        
        tail_sum_signal_start=1
        sum_del=strong+1+strong
        
        sum=0d0
        do i=1,strong
            sum=sum+signal(i)
            tail_sum_signal(i)=signal(i)
        enddo
        sum=sum+signal(strong+1)
        do i=strong+2,strong+1+strong
            sum=sum+signal(i)
        enddo
        
        do i=strong+1,signal_size-strong-1
            sum_arithmetic_mean=sum/sum_del
            
            sum=sum-tail_sum_signal(tail_sum_signal_start)+signal(i)
            tail_sum_signal(tail_sum_signal_start)=signal(i)
            tail_sum_signal_start=tail_sum_signal_start+1
            if(tail_sum_signal_start==strong+1) tail_sum_signal_start=1
            
            signal(i)=sum_arithmetic_mean
        enddo
        
        sum_arithmetic_mean=sum/sum_del
        signal(signal_size)=sum_arithmetic_mean
    endsubroutine arithmetic_mean_smoothing
endmodule smoothing_signal
