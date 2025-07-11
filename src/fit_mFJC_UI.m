 function fit_mFJC_UI(path_processed)
    % Get the list of files matching the naming scheme *_rupture.mat
    fileList = dir(fullfile(path_processed, '*_rupture.mat'));
    fileNames = {fileList.name};
    if isempty(fileNames)
        errordlg('No rupture files found in the specified directory.');
        return;
    end

    %% Create Main Figure and Panels
    fig = figure('Name','mFJC Fitting UI','NumberTitle','off',...
        'Position',[100 100 1200 600]);
    zoom off; pan off; rotate3d off;
    
    % Left panel for controls
    controlPanel = uipanel(fig, 'Title','Controls','FontSize',12,...
        'Units','normalized','Position',[0.01 0.01 0.28 0.98]);
    
    % Right panel for plotting
    plotPanel = uipanel(fig, 'Title','Force Curve','FontSize',12,...
        'Units','normalized','Position',[0.3 0.01 0.69 0.98]);
    ax = axes(plotPanel, 'Units','normalized','Position',[0.1 0.1 0.85 0.85]);
    
    %% UI Controls in the Control Panel
    % File selection using a popup menu
    uicontrol(controlPanel, 'Style','text','String','Select File:',...
        'Units','normalized','Position',[0.05 0.92 0.4 0.05]);
    popupFile = uicontrol(controlPanel, 'Style','popupmenu',...
        'String',fileNames, 'Units','normalized',...
        'Position',[0.5 0.92 0.4 0.05], 'Callback', @popupCallback);
    
    % Next Curve button (loads the next file automatically)
    nextButton = uicontrol(controlPanel, 'Style','pushbutton',...
        'String','Next','Units','normalized',...
        'Position',[0.92 0.95 0.07 0.03],'Callback',@nextCurveCallback);
    % Parameter inputs
    uicontrol(controlPanel, 'Style','text','String','Velocity (m/s):',...
        'Units','normalized','Position',[0.05 0.84 0.4 0.05]);
    vField = uicontrol(controlPanel, 'Style','edit','String','5e-7',...
        'Units','normalized','Position',[0.5 0.84 0.4 0.05]);

    uicontrol(controlPanel, 'Style','text','String','Contuour Length (m):',...
        'Units','normalized','Position',[0.05 0.76 0.4 0.05]);
    l_contField = uicontrol(controlPanel, 'Style','edit','String','5e-8',...
        'Units','normalized','Position',[0.5 0.76 0.4 0.05]);
    
    uicontrol(controlPanel, 'Style','text','String','Temperature (K):',...
        'Units','normalized','Position',[0.05 0.68 0.4 0.05]);
    TField = uicontrol(controlPanel, 'Style','edit','String','298',...
        'Units','normalized','Position',[0.5 0.68 0.4 0.05]);
    
    uicontrol(controlPanel, 'Style','text','String','Kuhn length (m):',...
        'Units','normalized','Position',[0.05 0.60 0.4 0.05]);
    l_kuhnField = uicontrol(controlPanel, 'Style','edit','String','5e-8',...
        'Units','normalized','Position',[0.5 0.60 0.4 0.05]);
    
    uicontrol(controlPanel, 'Style','text','String','Segment elasticity (N/m):',...
        'Units','normalized','Position',[0.05 0.52 0.4 0.05]);
    k_segField = uicontrol(controlPanel, 'Style','edit','String','2',...
        'Units','normalized','Position',[0.5 0.52 0.4 0.05]);
    
    % Instructions for interval selection
    uicontrol(controlPanel, 'Style','text',...
        'String','Click on the force curve to mark interval boundaries (min and max).',...
        'Units','normalized','Position',[0.05 0.42 0.9 0.05]);
    
    % Load Curve and Delete Fits buttons
    loadButton = uicontrol(controlPanel, 'Style','pushbutton',...
        'String','Load Curve','Units','normalized',...
        'Position',[0.1 0.35 0.38 0.08],'Callback',@loadCurveCallback);
    
    deleteFitsButton = uicontrol(controlPanel, 'Style','pushbutton',...
        'String','Delete Fits','Units','normalized',...
        'Position',[0.52 0.35 0.38 0.08],'Callback',@deleteFitsCallback);
    
    % Start, Submit, and Dismiss buttons
    fitButton = uicontrol(controlPanel, 'Style','pushbutton',...
        'String','Start Fit','Units','normalized',...
        'Position',[0.1 0.26 0.8 0.08],'Callback',@startFitCallback);
    
    submitButton = uicontrol(controlPanel, 'Style','pushbutton',...
        'String','Submit Fit','Units','normalized',...
        'Position',[0.1 0.15 0.35 0.08],'Callback',@submitFitCallback, 'Enable','off');
    
    dismissButton = uicontrol(controlPanel, 'Style','pushbutton',...
        'String','Dismiss Fit','Units','normalized',...
        'Position',[0.55 0.15 0.35 0.08],'Callback',@dismissFitCallback, 'Enable','off');
    
    %% Variables to store state
    selectedIntervals = [];  % Each row: [min_sep, max_sep] in extension units
    fitActive = false;       % Flag to prevent overlapping fits
    currentFit = struct('base_str','','parameters',[],'loading_rate',[],...
        'unbinding_force',[],'intervals',[],'z_model',[],'F_model',[]);
    currentData = struct('F',[],'z',[]);
    selectedFile = '';       % Current file name
    prefitHandles = [];      % Handles for pre-existing fit markers on the plot
    
    %% Set up the figure's WindowButtonDownFcn to capture clicks for interval selection.
    set(fig, 'WindowButtonDownFcn', @figureClickCallback);
    
    %% Callback functions
    
    function popupCallback(~,~)
        % When a file is selected from the popup menu, load it automatically.
        loadCurveCallback();
    end

    function nextCurveCallback(~,~)
        % Advance the popup selection to the next file and load it.
        idx = get(popupFile, 'Value');
        idx = mod(idx, length(fileNames)) + 1;  % wrap-around after last file
        set(popupFile, 'Value', idx);
        loadCurveCallback();
    end
    % 
    function loadCurveCallback(~,~)
    % Get the selected file from the popup menu.
    idx = get(popupFile, 'Value');
    fileNamesCell = get(popupFile, 'String');
    selectedFile = fileNamesCell{idx};
    full_path = fullfile(path_processed, selectedFile);
    try
        data = load(full_path);
    catch ME
        errordlg(['Error loading file: ' ME.message]);
        return;
    end
    if ~isfield(data, 'scaled_deflection') || ~isfield(data, 'shifted_height')
        errordlg('Selected file does not contain rupture event data.');
        return;
    end
    currentData.F = data.scaled_deflection;
    currentData.z = data.shifted_height;

    % Clear plot & reset state
    selectedIntervals = [];
    cla(ax);
    hold(ax, 'on');
    plot(ax, currentData.z, currentData.F, 'b', 'HitTest','off');
    xlabel(ax, 'Separation (m)');
    ylabel(ax, 'Force (N)');
    grid(ax, 'on');
    submitButton.Enable = 'off';
    dismissButton.Enable = 'off';
    fitActive = false;
    prefitHandles = [];

    % Load existing fits
    fitFile = fullfile(path_processed, 'fitResults.mat');
    if exist(fitFile, 'file')
        loadedData = load(fitFile);
        if isfield(loadedData, 'fitResults')
            allFits  = loadedData.fitResults;
            fileFits = allFits(strcmp({allFits.base_str}, selectedFile));

            % --- NEW: reload last initials if present
            if ~isempty(fileFits) && isfield(fileFits(end), 'initialParams')
                lastInit = fileFits(end).initialParams;  % [v, Lc, T, lk, ks]
                set(vField,      'String', num2str(lastInit(1)));
                set(l_contField, 'String', num2str(lastInit(2)));
                set(TField,      'String', num2str(lastInit(3)));
                set(l_kuhnField, 'String', num2str(lastInit(4)));
                set(k_segField,  'String', num2str(lastInit(5)));
            end

            % plot any pre‑computed fits
            for i = 1:numel(fileFits)
                if isfield(fileFits(i), 'z_model') && isfield(fileFits(i), 'F_model')
                    hLine = plot(ax, fileFits(i).z_model, fileFits(i).F_model, ...
                                 'r--', 'LineWidth', 1.5, 'HitTest','off');
                    prefitHandles = [prefitHandles; hLine];
                end
            end
        end
    end

    hold(ax, 'off');
    legend(ax, {'Retract','mFJC model'});
end


    % 
    %  function loadCurveCallback(~,~)
    % % Get the containing figure and resize it
    % 
    % 
    % % Get the selected file
    % idx = get(popupFile, 'Value');
    % fileNamesCell = get(popupFile, 'String');
    % selectedFile = fileNamesCell{idx};
    % full_path = fullfile(path_processed, selectedFile);
    % 
    % % Load data
    % try
    %     data = load(full_path);
    % catch ME
    %     errordlg(['Error loading file: ' ME.message]);
    %     return;
    % end
    % if ~isfield(data, 'scaled_deflection') || ~isfield(data, 'shifted_height')
    %     errordlg('Selected file does not contain rupture event data.');
    %     return;
    % end
    % 
    % % Prepare data in nm / pN
    % currentData.z = data.shifted_height   * 1e9;   % m → nm
    % currentData.F = data.scaled_deflection * 1e12; % N → pN
    % 
    % % Clear and plot
    % cla(ax);
    % hold(ax, 'on');
    % hPlot = plot(ax, currentData.z, currentData.F, 'b', 'LineWidth', 1.5);
    % set(hPlot, 'HitTest', 'off');
    % 
    % % Axis formatting
    % xlim(ax, [0, 150]);
    % ylim(ax, [-250,500]);
    % xlabel(ax, 'Separation (nm)', 'FontSize', 22);
    % ylabel(ax, 'Force (pN)',     'FontSize', 22);
    % grid(ax, 'on');
    % set(ax, 'FontSize', 22);
    % box(ax, 'on');
    % 
    % % Disable buttons while plotting
    % submitButton.Enable = 'off';
    % dismissButton.Enable = 'off';
    % 
    % % Overlay any pre‑computed mFJC fits (also scaled)
    % fitFile = fullfile(path_processed, 'fitResults.mat');
    % if exist(fitFile, 'file')
    %     loadedData = load(fitFile);
    %     if isfield(loadedData, 'fitResults')
    %         allFits  = loadedData.fitResults;
    %         fileFits = allFits(strcmp({allFits.base_str}, selectedFile));
    %         for i = 1:numel(fileFits)
    %             if isfield(fileFits(i), 'z_model') && isfield(fileFits(i), 'F_model')
    %                 % scale fit to nm / pN
    %                 zfit = fileFits(i).z_model   * 1e9;
    %                 Ffit = fileFits(i).F_model   * 1e12;
    %                 plot(ax, zfit, Ffit, 'LineStyle', '--', ...
    %                 'Color',    [1, 0.5, 0],  'LineWidth', 1.5, 'HitTest', 'off');
    %             end
    %         end
    %     end
    % end
% 
%     hold(ax, 'off');
% 
%     % Legend & export
%     legend(ax, {'Retract', 'mFJC model'}, 'FontSize', 22, 'Location', 'best');
% 
%     drawnow;
% 
%     % Export just the axes content at 800×600px (150 DPI → 800×600 = 5.33×4 in)
%     exportgraphics(ax, 'mFJC_fit.png', ...
%                    'Resolution', 150, ...
%                    'BackgroundColor', 'white');
% end


    function deleteFitsCallback(~,~)
        if isempty(selectedFile)
            return;
        end
        if ~isempty(prefitHandles)
            delete(prefitHandles);
            prefitHandles = [];
        end
        fitFile = fullfile(path_processed, 'fitResults.mat');
        if exist(fitFile, 'file')
            loadedData = load(fitFile);
            if isfield(loadedData, 'fitResults')
                allFits = loadedData.fitResults;
                newFits = allFits(~strcmp({allFits.base_str}, selectedFile));
                fitResults = newFits; 
                save(fitFile, 'fitResults');
            end
        end
        msgbox('All fits for the current force curve have been deleted.');
    end

    function figureClickCallback(~,~)
        if ~isempty(selectedIntervals) && ~isnan(selectedIntervals(end,2))
            return;
        end
        cp = get(ax, 'CurrentPoint');
        xClick = cp(1,1);
        yClick = cp(1,2);
        xLimits = get(ax, 'XLim');
        yLimits = get(ax, 'YLim');
        if xClick < xLimits(1) || xClick > xLimits(2) || yClick < yLimits(1) || yClick > yLimits(2)
            return;
        end
        hold(ax, 'on');
        if isempty(selectedIntervals)
            selectedIntervals = [xClick, NaN];
            plot(ax, xClick, interp1(currentData.z, currentData.F, xClick, 'linear','extrap'),...
                'ro','MarkerSize',8,'LineWidth',2, 'HitTest', 'off');
        else
            selectedIntervals(end,2) = xClick;
            plot(ax, xClick, interp1(currentData.z, currentData.F, xClick, 'linear','extrap'),...
                'ro','MarkerSize',8,'LineWidth',2, 'HitTest', 'off');
            plot(ax, [selectedIntervals(end,1) selectedIntervals(end,2)],...
                [interp1(currentData.z, currentData.F, selectedIntervals(end,1), 'linear','extrap') ...
                 interp1(currentData.z, currentData.F, selectedIntervals(end,2), 'linear','extrap')],...
                'r--','LineWidth',1.5, 'HitTest', 'off');
        end
        hold(ax, 'off');
    end

  function startFitCallback(~,~)
    if fitActive, return; end
    if isempty(currentData.F)
        errordlg('No force curve loaded. Please load a curve first.');
        return;
    end
    if isempty(selectedIntervals)
        errordlg('Please select at least one interval by clicking on the force curve.');
        return;
    end

    % Read initial guesses from the UI
    v          = str2double(get(vField,      'String'));
    L_cont     = str2double(get(l_contField, 'String'));
    T_val      = str2double(get(TField,      'String'));
    l_kuhn_val = str2double(get(l_kuhnField, 'String'));
    k_seg_val  = str2double(get(k_segField,  'String'));

    % --- NEW: store *all five* initials for reproducibility
    currentFit.initialParams = [v, L_cont, T_val, l_kuhn_val, k_seg_val];

    % Interval bounds
    min_sep = min(selectedIntervals(1,:));
    max_sep = max(selectedIntervals(1,:));

    % Run the mFJC fit
    try
        [loading_rate, unbinding_force, z_model, F_model, p_fit] = ...
            fit_mFJC(path_processed, selectedFile, ...
                     min_sep, max_sep, L_cont, T_val, l_kuhn_val, k_seg_val, v, ax);
    catch ME
        errordlg(['Fit failed: ' ME.message]);
        return;
    end

    % Show message
    msgbox(sprintf(['Fit complete!\nLoading rate: %g N/s\nUnbinding force: %g N\n' ...
                    'Contour Length: %g\nKuhn Length: %g\nSegment elasticity: %g'], ...
                    loading_rate, unbinding_force, p_fit(1), p_fit(2), p_fit(3)));

    % Store fit results
    currentFit.base_str       = selectedFile;
    currentFit.parameters     = [p_fit(1), p_fit(2), p_fit(3), T_val, 1];
    currentFit.loading_rate   = loading_rate;
    currentFit.unbinding_force= unbinding_force;
    currentFit.z_model        = z_model;
    currentFit.F_model        = F_model;
    currentFit.intervals      = selectedIntervals;

    submitButton.Enable = 'on';
    dismissButton.Enable = 'on';
    fitActive = true;
end


    function submitFitCallback(~,~)
    if isempty(currentFit.base_str)
        errordlg('No fit to submit.');
        return;
    end

    savePath   = fullfile(path_processed, 'fitResults.mat');
    fitResults = [];

    % Load existing results if they exist
    if exist(savePath, 'file')
        loadedData = load(savePath, 'fitResults');
        if isfield(loadedData, 'fitResults')
            fitResults = loadedData.fitResults;
        end
    end

    % Ensure the old entries and currentFit have the same fields
    if ~isempty(fitResults)
        oldFields = fieldnames(fitResults);
        newFields = fieldnames(currentFit);

        % 1) Add any new fields (in currentFit but not in old) to all old entries
        missingInOld = setdiff(newFields, oldFields);
        for f = missingInOld'
            [fitResults.(f{1})] = deal([]);  %#ok<AGROW>
        end

        % 2) Add any old fields (in old entries but not in currentFit) to currentFit
        missingInNew = setdiff(oldFields, newFields);
        for f = missingInNew'
            currentFit.(f{1}) = [];
        end
    end

    % Now safe to concatenate
    fitResults = [fitResults; currentFit];  

    % Save back to disk
    save(savePath, 'fitResults');

    % Reset state
    currentFit        = struct('base_str','',  'initialParams',[],'parameters',[], ...
                               'loading_rate',[],'unbinding_force',[], ...
                               'intervals',[],'z_model',[],    'F_model',[]);
    selectedIntervals = [];
    fitActive         = false;
    submitButton.Enable  = 'off';
    dismissButton.Enable = 'off';
end



    function dismissFitCallback(~,~)
        currentFit = struct('base_str','','parameters',[],'loading_rate',[],...
            'unbinding_force',[],'intervals',[],'z_model',[],'F_model',[]);
        selectedIntervals = [];
        fitActive = false;
        submitButton.Enable = 'off';
        dismissButton.Enable = 'off';
        %msgbox('Fit dismissed.');
    end
end