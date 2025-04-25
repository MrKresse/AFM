function powerspec(path_processed, base_str, velocity, baseline_min, baseline_max)
    load(fullfile(path_processed, [base_str, '_approach_scaled.mat']));
    
    % Time vector from height and velocity
    height_time = shifted_height ./ velocity;
    baseline_min_time = baseline_min / velocity;
    baseline_max_time = baseline_max / velocity;
    
    % Apply baseline mask
    mask = (height_time >= baseline_min_time) & (height_time <= baseline_max_time);
    force_curve = scaled_deflection(mask);
    
    % Compute the FFT of the force curve
    fft_result = fft(force_curve);
    
    % Compute the power spectrum (positive frequencies, excluding DC component)
    n = length(force_curve); % Number of points
    power_spectrum = (abs(fft_result(2:floor(n/2)+1))) / n; % Start from the second frequency
    
    % Frequency vector (positive frequencies only, excluding DC component)
    sampling_interval = mean(diff(height_time(mask))); % Approximate sampling interval
    fs = 1 / sampling_interval; % Sampling frequency
    frequencies = -(1:floor(n/2)) * (fs / n); % Start from the first positive frequency
    
    % Save the results to a .mat file
    output_file = fullfile(path_processed, [base_str, '_power_spectrum.mat']);
    save(output_file, 'frequencies', 'power_spectrum');
end


