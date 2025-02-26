module system
implicit none
    public::end_program_pause,end_program,print_info,print_warning,print_error
    
    logical(1),parameter::with_logging=.true.
    
    interface end_program_pause
        procedure end_program_pause_default
        procedure end_program_pause_with_message
    endinterface end_program_pause
    
    interface end_program
        procedure end_program_default
        procedure end_program_with_message
    endinterface end_program
    
    private::end_program_pause_default,end_program_pause_with_message,&
        end_program_default,end_program_with_message,&
        logging,with_logging
    contains
    
    subroutine logging(type_log, function_log, message)
    implicit none
        character(len=*),intent(in)::type_log
        character(len=*),intent(in)::function_log
        character(len=*),intent(in)::message
        
        integer(4) logfile
        character(len=*),parameter::filename="log.txt"
        
        open(newunit=logfile, file=filename, status='unknown', action='write', position='append')
        write(logfile, '(A)', advance='no') type_log//" "//function_log//" "//message
        close(logfile)
    endsubroutine logging
    
    subroutine end_program_pause_default()
    implicit none
        call end_program_pause_with_message("print Enter to continue...")
    endsubroutine end_program_pause_default
    subroutine end_program_pause_with_message(message)
    implicit none
        character(len=*),intent(in)::message
        
        if(with_logging) call logging("info", "main", "end program: "//message)
        
        print*
        print*
        print*
        print*,message
        pause " "
    endsubroutine end_program_pause_with_message
    
    subroutine end_program_default()
    implicit none
        call end_program_with_message("print Enter to continue...")
    endsubroutine end_program_default
    subroutine end_program_with_message(message)
    implicit none
        character(len=*),intent(in)::message
        
        if(with_logging) call logging("info", "main", "end program: "//message)
    endsubroutine end_program_with_message
    
    subroutine print_info(function_name, message)
    implicit none
        character(len=*),intent(in)::function_name
        character(len=*),intent(in)::message
        
        if(with_logging) call logging("info", function_name, message)
        print*,"info in ",function_name,": ",message
    endsubroutine print_info
    subroutine print_warning(function_name, message)
    implicit none
        character(len=*),intent(in)::function_name
        character(len=*),intent(in)::message
        
        if(with_logging) call logging("warning", function_name, message)
        print*,"warning in ",function_name,": ",message
    endsubroutine print_warning
    subroutine print_error(function_name, message)
    implicit none
        character(len=*),intent(in)::function_name
        character(len=*),intent(in)::message
        
        if(with_logging) call logging("error", function_name, message)
        print*,"error in ",function_name,": ",message
        pause "print Enter to continue and stop program..."
        stop
    endsubroutine print_error
endmodule system
