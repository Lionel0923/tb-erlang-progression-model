% ====================== plot_simple_vs_erlang_overlay_clean_logsafe.m ======================
% Clean overlay of SIMPLE vs ERLANG posterior mean + 95% UI on ONE figure.
% Observations + CI are computed EXACTLY as in calibration (NO truncation).
% Log-safety is enforced ONLY at plotting time (display floor), without changing CI logic.

clear; clc; close all;

%% 0) Data files (same as calibration scripts)
dataFiles = {
  '~/data_nti1_preweight.csv'
  '~/data_nti2_preweight.csv'
  '~/data_nti3_preweight.csv'
};
nDatasets = numel(dataFiles);

%% 1) Load saved posteriors (Simple + Erlang)
S1 = load('postSamples_simpleTB_uniform_adaptcov.mat', 'tPred','Ypred');
S2 = load('postSamples_subclinical_uniform_adaptcov.mat', 'tPred','Ypred','k'); %#ok<NASGU>

tPred_simple = S1.tPred(:);
tPred_erlang = S2.tPred(:);

Ypred_simple = S1.Ypred;
Ypred_erlang = S2.Ypred;

%% 2) Observations + CI (EXACTLY match SIMPLE calibration: binofit CI; NO truncation)
tData   = cell(nDatasets,1);
yObs    = cell(nDatasets,1);
ciLower = cell(nDatasets,1);
ciUpper = cell(nDatasets,1);

for d = 1:nDatasets
    T          = readtable(dataFiles{d});
    atRisk     = T.atrisk;
    casesObs   = T.cases;
    tData{d}   = T.year;

    % EXACTLY as calibration:
    yObs{d}    = casesObs ./ atRisk;
    [~, pci]   = binofit(casesObs, atRisk, 0.05);
    ciLower{d} = yObs{d} - pci(:,1);
    ciUpper{d} = pci(:,2) - yObs{d};
end

%% 3) Posterior summaries (Simple)
yMean_simple = mean(Ypred_simple, 1);
yLow_simple  = quantile(Ypred_simple, 0.025, 1);
yHigh_simple = quantile(Ypred_simple, 0.975, 1);

%% 4) Posterior summaries (Erlang)
yMean_erlang = mean(Ypred_erlang, 1);
yLow_erlang  = quantile(Ypred_erlang, 0.025, 1);
yHigh_erlang = quantile(Ypred_erlang, 0.975, 1);

%% 5) Log-safe plotting (NO change to CI computation; only display-floor)
% If any yObs is 0, or implied CI lower hits 0, log-scale plotting can break.
% We DO NOT alter ciLower/ciUpper formulas; we only create "plot versions".

plotFloor = 1e-6;  % purely for display; choose << your y-limits lower bound

tData_plot   = tData;
yObs_plot    = yObs;
ciLower_plot = ciLower;
ciUpper_plot = ciUpper;

for d = 1:nDatasets
    y  = yObs{d};
    lo = y - ciLower{d};   % implied absolute lower bound (pci(:,1))
    hi = y + ciUpper{d};   % implied absolute upper bound (pci(:,2))

    % For log-scale display: replace any non-positive values with plotFloor
    y_plot  = y;
    lo_plot = lo;
    hi_plot = hi;

    y_plot(y_plot <= 0)   = plotFloor;
    lo_plot(lo_plot <= 0) = plotFloor;
    hi_plot(hi_plot <= 0) = plotFloor;

    % Convert back to errorbar lengths around the plotted y
    ciLower_plot{d} = y_plot - lo_plot;
    ciUpper_plot{d} = hi_plot - y_plot;
    yObs_plot{d}    = y_plot;
end

% Also make UI bands log-safe for display (rare but safe)
yLow_simple_plot  = max(yLow_simple(:),  plotFloor);
yHigh_simple_plot = max(yHigh_simple(:), plotFloor);
yMean_simple_plot = max(yMean_simple(:), plotFloor);

yLow_erlang_plot  = max(yLow_erlang(:),  plotFloor);
yHigh_erlang_plot = max(yHigh_erlang(:), plotFloor);
yMean_erlang_plot = max(yMean_erlang(:), plotFloor);

%% 6) Plot overlay (ONE figure)
figure('Position',[180 180 860 540]); hold on; box on;

% ---- Observations styling ----
markers   = {'^','o','s'};
obsColors = {
    [0.00 0.60 0.50];   
    [0.80 0.20 0.60];   % magent
    [0.95 0.55 0.05];   % vivid orange
};
labels    = {'NTI1','NTI2','NTI3'};

% >>> ADD: horizontal jitter (standard in papers to avoid overlap) <<<
xOffset = [-0.11, 0, 0.11];   % years; small enough to still look "same year"

hObs = gobjects(nDatasets,1);
for d = 1:nDatasets
    xShift = tData_plot{d} + xOffset(d);

    hObs(d) = errorbar(xShift, yObs_plot{d}, ciLower_plot{d}, ciUpper_plot{d}, markers{d}, ...
        'Color', obsColors{d}, ...
        'MarkerFaceColor', obsColors{d}, ...
        'MarkerEdgeColor', obsColors{d}, ...
        'LineWidth', 1.3, ...
        'CapSize', 6);
end

% ---- Model UI colors (stronger contrast, still clean) ----
meanColor_simple = [0.00 0.32 0.64];      % rich deep blue
meanColor_erlang = [0.75 0.10 0.10];      % strong medical red

uiColor_simple   = [0.60 0.78 0.95];      % light blue
uiColor_erlang   = [0.95 0.70 0.70];      % light red

alphaUI = 0.20;

% --- Simple UI band ---
hUI_simple = fill([tPred_simple; flipud(tPred_simple)], ...
                  [yLow_simple_plot; flipud(yHigh_simple_plot)], ...
                  uiColor_simple, 'EdgeColor','none');
set(hUI_simple, 'FaceAlpha', alphaUI);

% --- Erlang UI band ---
hUI_erlang = fill([tPred_erlang; flipud(tPred_erlang)], ...
                  [yLow_erlang_plot; flipud(yHigh_erlang_plot)], ...
                  uiColor_erlang, 'EdgeColor','none');
set(hUI_erlang, 'FaceAlpha', alphaUI);

% --- Mean lines ---
hMean_erlang = plot(tPred_erlang, yMean_erlang_plot, '-',  'Color', meanColor_erlang, 'LineWidth', 2.6);
hMean_simple = plot(tPred_simple, yMean_simple_plot, '-.', 'Color', meanColor_simple,  'LineWidth',2.6);

% ---- Axes formatting ----
set(gca,'YScale','log');
grid off;
set(gca,'MinorGridLineStyle','-');

xlabel('Years since infection');
ylabel('Annual net inflow to subclinical (rate)');

% Log-safe limits (NO zero!)
ylim([1e-4 5e-2]);
yticks([1e-4 3e-4 1e-3 3e-3 1e-2 3e-2]);

title('Fit to annual net inflow into subclinical: Simple vs Gamma-distributed');

legend([hObs(:); hUI_simple; hUI_erlang; hMean_simple; hMean_erlang], ...
       {labels{:}, ...
        '95% UI (Simple)', '95% UI (Gamma-distributed)', ...
        'Posterior mean (Simple)', 'Posterior mean (Gamma-distributed)'}, ...
       'Location','northeast');

% (Optional) white background
set(gcf,'Color','w'); set(gca,'Color','w');

