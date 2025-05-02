clear;
load("data.mat")
% hold on;
% plot(data(:,1),log(data(:,2)))
% plot(data(:,1), log(data(:,4)), 'r-' )
% hold off
plot(data(3:1000,1), log(data(3:1000,3)))
addpath("./src/")
show_spectra(data(:,2))

