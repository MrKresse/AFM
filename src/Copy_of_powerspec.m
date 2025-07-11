function Copy_of_powerspec(path_processed, base_str, velocity, baseline_min, baseline_max)
    load(fullfile(path_processed, [base_str, '_approach_scaled.mat']));
    
    % Time vector from height and velocity
    height_time = shifted_height ./ velocity;
    baseline_min_time = baseline_min / velocity;
    baseline_max_time = baseline_max / velocity;
    
    % Apply baseline mask
    mask = (height_time >= baseline_min_time) & (height_time <= baseline_max_time);
    force_curve = scaled_deflection(mask);
    
    % Compute the power spectrum (positive frequencies, excluding DC component)
    L = length(force_curve); % Number of points
    % Frequency vector 
    sampling_interval = abs(mean(diff(height_time(mask)))); % Approximate sampling interval
    Fs = 1 / sampling_interval; % Sampling frequency

    Y = fft(force_curve);
    P2 = abs(Y/L);
    P1 = P2(1:floor(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    frequencies = Fs/L*(0:(L/2));
    
    power_spectrum = P1;
    output_file = fullfile(path_processed, [base_str, '_power_spectrum.mat']);
    save(output_file, 'frequencies', 'power_spectrum');
end

