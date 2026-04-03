%% ================= Panel: % difference in cumulative TB cases averted =================
% Using MEANS only (from your tables). Definition:
%   Averted_y = BaselineMean_y - ACFMean_y
%   CumAverted_N = sum_{y=2024}^{2024+N-1} Averted_y
%   %Diff_N = (CumAverted_Erlang_N - CumAverted_Simple_N) / CumAverted_Simple_N * 100

% --- Inputs (means) ---
yrs_int = 2024:2026; % intervention years only (2024–26)

simple_base_mean = [167.650 161.387 155.502];
erlang_base_mean = [172.398 166.596 161.071];

simple_acf_mean  = [149.491 134.832 124.955];
erlang_acf_mean  = [135.161 119.336 110.846];

% --- Compute annual + cumulative averted (means) ---
simple_avert = simple_base_mean - simple_acf_mean;   % 2024,2025,2026
erlang_avert = erlang_base_mean - erlang_acf_mean;

simple_cum = cumsum(simple_avert);  % after 1y,2y,3y
erlang_cum = cumsum(erlang_avert);

pct_diff = (erlang_cum - simple_cum) ./ simple_cum * 100;  % % difference vs Simple

% --- Plot (make this your "additional panel") ---
% If you're embedding into an existing 2-panel Fig 3, call: subplot(1,2,2);
figure('Position',[120 120 700 420]); hold on; box off; grid on;

% Lancet-like (ggsci "lancet") deep blue
lancetBlue = [55 76 106]/255;   % #00468B

x = 1:3;
b = bar(x, pct_diff, 0.62, 'FaceColor', lancetBlue, 'EdgeColor', 'none');

yline(0,'k-','LineWidth',1);

set(gca, ...
    'XTick', x, ...
    'XTickLabel', {'After 1 year','After 2 years','After 3 years'}, ...
    'TickDir','out', ...
    'LineWidth',1, ...
    'FontName','Arial', ...
    'FontSize',11);

ylabel('% difference in cumulative cases averted');
title('ACF (2024–26): cumulative impact difference');

% Optional: annotate bar values
for i = 1:numel(pct_diff)
    text(x(i), pct_diff(i), sprintf('%.1f%%', pct_diff(i)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontName','Arial', 'FontSize',10);
end

% Tight y-limits with padding
yl = ylim;
ylim([min(0,yl(1)) - 5, yl(2) + 10]);

% If you want consistent panel label styling (e.g., "B"), uncomment:
% text(0.02, 0.98, 'B', 'Units','normalized', 'FontWeight','bold', 'FontSize',14);

% --- Quick sanity printout (optional) ---
disp(table((1:3)', simple_cum(:), erlang_cum(:), pct_diff(:), ...
    'VariableNames', {'YearsIntoIntervention','CumAverted_Simple','CumAverted_Erlang','PctDiff_ErlangVsSimple'}));
