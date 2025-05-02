clear;
% --- minimal change version ---
fname = './thermal-noise-data_vDeflection_2025.03.20-12.01.13.tnd';

% --- auto‐detect the first numeric line, then read data ---
if ~isfile(fname)
    error('Cannot find "%s" in %s.', fname, pwd);
end

% 1) Slurp all lines as text
fid = fopen(fname,'r');
allLines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
lines = allLines{1};

% 2) Find the first line that parses as 4 floats
isNumLine = cellfun(@(s) numel(sscanf(s, '%f')) == 4, lines);
startIdx  = find(isNumLine, 1);
if isempty(startIdx)
    error('No line with exactly 4 numeric values was found.');
end

% 3) Re‐open and skip up to that line, then read your 4 columns
fid = fopen(fname,'r');
C = textscan( ...
    fid, ...
    '%f %f %f %f', ...
    'HeaderLines', startIdx-1, ...   % skip all non-numeric lines
    'CollectOutput', true ...
);
fclose(fid);

% 4) Extract into matrix and save
data = C{1};  % Nx4
if isempty(data)
    error('textscan still returned no data—check file encoding or delimiters.');
end

disp(['Loaded ', num2str(size(data,1)), ' rows of data.']);
% --- after you’ve built `lines` and `isNumLine` ---

idx = find(isNumLine);          % indices of lines with 4 numbers
n   = numel(idx);
data = zeros(n,4);

for k = 1:n
    vals = sscanf(lines{idx(k)}, '%f');   % parse the k-th numeric line
    if numel(vals)==4
        data(k,:) = vals';                % store row
    else
        error('Line %d didn''t parse into 4 numbers.', idx(k));
    end
end

save('data.mat','data');
disp(['Loaded ',num2str(n),' rows of data.']);
