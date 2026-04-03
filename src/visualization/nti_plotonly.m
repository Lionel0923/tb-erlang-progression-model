clear; clc; close all;

%% 1) Load saved posterior samples/results
load('postSamples_subclinical_uniform_adaptcov.mat', ...
    'chainThin','tData','k','yObs','x0_list','tPred','Ypred');

%% 2) Rebuild observational error bars from the original csv files
% This does NOT refit the model. It only reconstructs the plotted data.
dataFiles = {
    '~/data_nti1_preweight.csv'
    '~/data_nti2_preweight.csv'
    '~/data_nti3_preweight.csv'
};

nDatasets = numel(dataFiles);
ciLower = cell(nDatasets,1);
ciUpper = cell(nDatasets,1);

for d = 1:nDatasets
    T = readtable(dataFiles{d});
    atRisk   = T.atrisk;
    casesObs = T.cases;

    lambda   = casesObs;
    N0       = atRisk;

    lo_cases = 0.5 * chi2inv(0.025, 2*lambda);
    hi_cases = 0.5 * chi2inv(0.975, 2*(lambda+1));

    rate_lo  = lo_cases ./ N0;
    rate_hi  = hi_cases ./ N0;

    ciLower{d} = yObs{d} - rate_lo;
    ciUpper{d} = rate_hi - yObs{d};
end

%% 3) Posterior summary from saved Ypred
yMean = mean(Ypred, 1);
yLow  = quantile(Ypred, 0.025, 1);
yHigh = quantile(Ypred, 0.975, 1);

%% 4) Reproduce the same plot
figure('Position',[200 200 760 470]); 
hold on; grid on;

markers = {'^','o','s'};
colors  = {'b','m',[1 .5 0]};
labels  = {'NTI1','NTI2','NTI3'};

obsHandles = gobjects(nDatasets,1);
for d = 1:nDatasets
    obsHandles(d) = errorbar(tData{d}, yObs{d}, ciLower{d}, ciUpper{d}, ...
        markers{d}, ...
        'Color', colors{d}, ...
        'MarkerFaceColor', colors{d}, ...
        'LineWidth', 1.3);
end

bandHandle = fill([tPred; flipud(tPred)], [yLow'; flipud(yHigh')], ...
    'c', 'EdgeColor', 'none');
set(bandHandle, 'FaceAlpha', 0.2);

meanHandle = plot(tPred, yMean, '-k', 'LineWidth', 2);

xlabel('Years since infection');
ylabel('Annual net inflow to subclinical (rate)');
title('Fit to annual net inflow into subclinical (Gamma-distributed model)');

legend([obsHandles' bandHandle meanHandle], ...
    {'NTI1','NTI2','NTI3','95% posterior band','Posterior mean'}, ...
    'Location', 'northeast');
