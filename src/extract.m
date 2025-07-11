function extract(path_raw, path_incoming, base_str, heightName, deflectName)
    % Open the JPK force curve data
    force_curve = open_JPK(fullfile(path_raw, base_str));

    % Build cell-array of channel names and find indices
    chanNames = { force_curve.Channel_name };
    iH = find(strcmp(chanNames, heightName), 1);
    iD = find(strcmp(chanNames, deflectName), 1);

    % If name lookup fails, warn & exit early
    if isempty(iH)
        fprintf('Warning: height channel ''%s'' not found. Skipping extract.\n', heightName);
        return
    end
    if isempty(iD)
        fprintf('Warning: deflection channel ''%s'' not found. Skipping extract.\n', deflectName);
        return
    end

    % Helper to safely pull out Data_nominal or Data_Raw
    function vec = safeGet(idx, phase, field)
        vec = [];
        if isfield(force_curve(idx), phase) ...
           && isstruct(force_curve(idx).(phase)) ...
           && isfield(force_curve(idx).(phase), field)
            raw = force_curve(idx).(phase).(field);
            if ~isempty(raw)
                vec = raw;
            end
        end
    end

    % Extract or get empty if missing
    ah = safeGet(iH, 'extend', 'Data_nominal');
    ad = safeGet(iD, 'extend', 'Data_Raw');
    rh = safeGet(iH, 'retract', 'Data_nominal');
    rd = safeGet(iD, 'retract', 'Data_Raw');

    % Zero-offset
    if ~isempty(ah), ah = ah - min(ah); end
    if ~isempty(rh), rh = rh - min(rh); end

    % Save approach if both vectors exist
    if isempty(ah) || isempty(ad)
        fprintf('Warning: approach data for ''%s'' or ''%s'' is empty. Not saving approach.\n', heightName, deflectName);
    else
        height     = ah;
        deflection = ad;
        save(fullfile(path_incoming, [base_str '_approach.mat']), 'height', 'deflection');
    end

    % Save retract if both vectors exist
    if isempty(rh) || isempty(rd)
        fprintf('Warning: retract data for ''%s'' or ''%s'' is empty. Not saving retract.\n', heightName, deflectName);
    else
        height     = rh;
        deflection = rd;
        save(fullfile(path_incoming, [base_str '_retract.mat']), 'height', 'deflection');
    end
end


   
    
    

