function extract(path_raw, path_incoming, base_str, channel)

    % Open the JPK force curve data
    force_curve = open_JPK(fullfile(path_raw, base_str));
    
    % Extract the relevant data for approach and deflection
    approach_height = force_curve(channel).extend.Data_nominal;
    approach_deflection = force_curve(5).extend.Data_distance;
    retract_height = force_curve(channel).retract.Data_nominal;
    retract_deflection = force_curve(5).retract.Data_distance;
    
    approach_height= approach_height-min(approach_height);
    retract_height=retract_height-min(retract_height);
    approach_file = fullfile(path_incoming, [base_str, '_approach.mat']);
    deflect_file = fullfile(path_incoming, [base_str, '_retract.mat']);
    
    % Save the approach and deflection data as .mat files
    height = approach_height;
    deflection = approach_deflection;
    save(approach_file, 'height', 'deflection');
    height = retract_height;
    deflection = retract_deflection;
    save(deflect_file, 'height', 'deflection');
end   
    
    

