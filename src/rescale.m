function rescale(path_incoming, path_processed, base_str, mc,bc,m,b, r_mc,r_bc,r_m,r_b,k)
    approach_file = fullfile(path_incoming, [base_str, '_approach.mat']);
    approach= load(approach_file, 'height', 'deflection');
    retract_file = fullfile(path_incoming, [base_str, '_retract.mat']);
    retract=load(retract_file, 'height', 'deflection');

    shifted_approach_height =(approach.height + (-bc+ b)/(mc-m)); %0 should be where the fits of const compliance and baseline intersect
    scaled_approach_deflection =(approach.deflection-b)/-mc;%0 at baseline
    shifted_approach_height = shifted_approach_height-scaled_approach_deflection;%separation
    scaled_approach_deflection = k*scaled_approach_deflection;%force
            
    shifted_retract_height =(retract.height + (-r_bc+ r_b)/(r_mc-r_m)); %0 should be where the fits of const compliance and baseline intersect
    scaled_retract_deflection =(retract.deflection-r_b)/-r_mc;%0 at baseline
    shifted_retract_height = shifted_retract_height - scaled_retract_deflection;
    scaled_retract_deflection = k*scaled_retract_deflection;

      % Save rescaled approach data
    approach_scaled_file = fullfile(path_processed, [base_str, '_approach_scaled.mat']);
    shifted_height= shifted_approach_height;%+scaled_approach_deflection;
    scaled_deflection = scaled_approach_deflection;
    save(approach_scaled_file, 'shifted_height', 'scaled_deflection');
    
    % Save rescaled retract data
    retract_scaled_file = fullfile(path_processed, [base_str, '_retract_scaled.mat']);
    shifted_height= shifted_retract_height;%+scaled_retract_deflection;
    scaled_deflection = scaled_retract_deflection;
    save(retract_scaled_file, 'shifted_height', 'scaled_deflection');
end

