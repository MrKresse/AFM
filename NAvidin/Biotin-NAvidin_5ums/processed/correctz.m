load('fitResults_old.mat')

for i = 1:320
    fitResults(i).parameters(5) = 1
    fitResults(i).loading_rate=fitResults(i).loading_rate/5e-6
end

save("fitResults.mat", 'fitResults')


