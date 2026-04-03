% plot_subclinicalTB_posterior_from_saved_nonneg.m
% Load saved Erlang posterior + predictions and reproduce the same plot,
% with NON-NEGATIVE observation error bars using the same Poisson-chi2 logic
% as in bayes_calibration_subclinical_uniform_adaptcov.m

clear; clc; close all;

%% 0) (REQUIRED) Data files (same as in calibration script)
dataFiles = {
  '/Users/zhouzhichaozhichao/Downloads/Output/data_nti1_preweight.csv'
  '/Users/zhouzhichaozhichao/Downloads/Output/data_nti2_preweight.csv'
  '/Users/zhouzhichaozhichao/Downloads/Output/data_nti3_preweight.csv'
};

%% 1) Load saved posterior
S = load('postSamples_50_uniform_adaptcov.mat', ...
         'chainThin','tData','k','yObs','x0_list','tPred','Ypred');

tData   = S.tData;
yObs    = S.yObs;
x0_list = S.x0_list;
tPred   = S.tPred;
Ypred   = S.Ypred;
k       = S.k; %#ok<NASGU>

nDatasets = numel(tData);

%% 2) Recreate observation error bars using Poisson (chi-square) CI (NON-NEGATIVE)
% This matches the logic in your calibration script:
% lo_cases = 0.5*chi2inv(0.025, 2*lambda);
% hi_cases = 0.5*chi2inv(0.975, 2*(lambda+1));
% rate_lo  = lo_cases ./ N0;   rate_hi = hi_cases ./ N0;
% ciLower  = yObs - rate_lo;   ciUpper = rate_hi - yObs;
ciLower = cell(nDatasets,1);
ciUpper = cell(nDatasets,1);

for d = 1:nDatasets
    T        = readtable(dataFiles{d});
    atRisk   = T.atrisk;
    casesObs = T.cases;

    % (Optional safety) ensure consistency with saved yObs/tData
    % If you don't want checks, you can delete these lines.
    if ~isequal(T.year(:), tData{d}(:))
        warning('Dataset %d: CSV years do not match saved tData. Using CSV years for CI.', d);
    end

    yObs_csv = casesObs ./ atRisk;  % should match yObs{d}

    lambda   = casesObs;           % Poisson mean
    N0       = atRisk;

    lo_cases = 0.5 * chi2inv(0.025, 2*lambda);
    hi_cases = 0.5 * chi2inv(0.975, 2*(lambda+1));

    rate_lo  = lo_cases ./ N0;     % >= 0
    rate_hi  = hi_cases ./ N0;

    % Convert to errorbar "lengths" (must be non-negative)
    ciLower{d} = max(yObs_csv - rate_lo, 0);  % >=0
    ciUpper{d} = max(rate_hi - yObs_csv, 0);  % >=0

    % Also overwrite yObs{d} with CSV version to ensure exact match to CI
    yObs{d} = yObs_csv;
end

%% 3) Posterior summary for prediction band/mean (from saved Ypred)
yMean = mean(Ypred, 1);
yLow  = quantile(Ypred, 0.025, 1);
yHigh = quantile(Ypred, 0.975, 1);

%% 4) Plot (same look as your original Erlang figure) — non-negative axis
figure('Position',[200 200 760 470]); hold on; grid on;

markers = {'^','o','s'};
colors  = {'b','m',[1 .5 0]};
labels  = {'NTI1','NTI2','NTI3'};

for d = 1:nDatasets
    errorbar(tData{d}, yObs{d}, ciLower{d}, ciUpper{d}, markers{d}, ...
        'Color', colors{d}, 'MarkerFaceColor', colors{d}, 'LineWidth', 1.3);
end

h = fill([tPred; flipud(tPred)], [yLow'; flipud(yHigh')], ...
         'c', 'EdgeColor','none');
set(h, 'FaceAlpha', 0.2);

plot(tPred, yMean, '-k', 'LineWidth', 2);

% --- after plotting (before/after legend 都行) ---
set(gca,'YScale','log');

ylim([1e-4 5e-2]);

yticks([1e-4 3e-4 1e-3 3e-3 1e-2 3e-2]);

grid on; set(gca,'MinorGridLineStyle','-'); 



xlabel('Years since infection');
ylabel('Annual net inflow to subclinical (rate)');
legend([labels, {'95% posterior band','Posterior mean'}], 'Location', 'northeast');
title('Fit to annual net inflow into subclinical (Gamma-distributed model)');

% Enforce non-negative y-axis display (recommended for rates)
yl = ylim;
ylim([0, max(yl(2), 0.04)]);  % keep top similar; adjust 0.04 if you want
