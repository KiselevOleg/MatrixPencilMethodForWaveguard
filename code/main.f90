program main
    use system,only:end_program_pause,end_program,print_info,print_warning,print_error
    use pre_and_post_runtime_actions,only:init,destructor
    use runtime_application_actions,only:test_K,load_experimental_measurements_view,count_matrix_pencil_method_distersion_curve_graphics,&
        count_distersion_curve_for_K_graphics,count_complex_distersion_curve_for_K_graphics,count_complex_det_A_matrix_in_K_graphics,&
        count_wavelet_transform_graphics,count_wavelet_transform_smoothing,find_material_properties_from_experimental_measurements
implicit none
    call init()
    
    call find_material_properties_from_experimental_measurements()
    !call count_wavelet_transform_smoothing()
    !call count_wavelet_transform_graphics()
    !call count_complex_det_A_matrix_in_K_graphics()
    !call count_complex_distersion_curve_for_K_graphics()
    !!!call count_distersion_curve_for_K_graphics()
    !call count_matrix_pencil_method_distersion_curve_graphics()
    !call load_experimental_measurements_view()
    !call test_K()
    
    call destructor()
    call end_program_pause()
    call end_program()
endprogram main
