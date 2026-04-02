%% ---- UI (95%) for forecasts (2024–2033) ----
% Helper: 95% UI using percentiles
pct = @(A) prctile(A, [2.5 97.5], 1);

% Baseline UIs
tmp = pct(Yf_simple);   fcast_simple_lo = tmp(1,:); fcast_simple_hi = tmp(2,:);
tmp = pct(Yf_erlang);   fcast_erlang_lo = tmp(1,:); fcast_erlang_hi = tmp(2,:);
tmp = pct(Yf_eFromS);   fcast_eFromS_lo = tmp(1,:); fcast_eFromS_hi = tmp(2,:);

% ACF UIs
tmp = pct(YfS_acf);     fcastS_acf_lo   = tmp(1,:); fcastS_acf_hi   = tmp(2,:);
tmp = pct(YfEr_acf);    fcastEr_acf_lo  = tmp(1,:); fcastEr_acf_hi  = tmp(2,:);
tmp = pct(YfE_acf);     fcastE_acf_lo   = tmp(1,:); fcastE_acf_hi   = tmp(2,:);

%% ---- Print key-year summary (mean + 95% UI) ----
keyYears = [2024 2025 2026 2027 2028 2029 2030 2031 2032 2033];
[~, keyIdx] = ismember(keyYears, yearsF);

fprintf('\n===== Forecast summary: mean (95%% UI) =====\n');
for ii = 1:numel(keyYears)
    k = keyIdx(ii);
    fprintf('\nYear %d\n', keyYears(ii));
    fprintf('  Simple baseline: %.3f (%.3f, %.3f)\n', fcast_simple_mean(k), fcast_simple_lo(k), fcast_simple_hi(k));
    fprintf('  Erlang baseline: %.3f (%.3f, %.3f)\n', fcast_erlang_mean(k), fcast_erlang_lo(k), fcast_erlang_hi(k));
    fprintf('  eFromS baseline: %.3f (%.3f, %.3f)\n', fcast_eFromS_mean(k), fcast_eFromS_lo(k), fcast_eFromS_hi(k));
    fprintf('  Simple ACF:      %.3f (%.3f, %.3f)\n', fcastS_acf_mean(k),  fcastS_acf_lo(k),  fcastS_acf_hi(k));
    fprintf('  Erlang ACF:      %.3f (%.3f, %.3f)\n', fcastEr_acf_mean(k), fcastEr_acf_lo(k), fcastEr_acf_hi(k));
    fprintf('  eFromS ACF:      %.3f (%.3f, %.3f)\n', fcastE_acf_mean(k),  fcastE_acf_lo(k),  fcastE_acf_hi(k));
end
fprintf('\n==========================================\n\n');
