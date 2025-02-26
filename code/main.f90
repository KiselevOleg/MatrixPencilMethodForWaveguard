program main
    use system,only:end_program_pause,end_program,print_info,print_warning,print_error
    use pre_and_post_runtime_actions,only:init,destructor
    use runtime_application_actions,only:test_K
implicit none
    call init()
    
    call test_K()
    
    call destructor()
    call end_program_pause()
    call end_program()
endprogram main
