% Set up the environment
clear; clc;

% Get the current script's directory
currentDir = fileparts(mfilename('fullpath'));

% Add the src directory to the path
addpath(fullfile(currentDir, 'src'));

% Start the GUI
buildGUI();