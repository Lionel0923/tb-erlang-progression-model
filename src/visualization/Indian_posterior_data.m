clear; clc; close all; rng(20250609);

%% ---- Load calibrated outputs ----
S = load('postSamples_simple_clinical_AM.mat');   % Simple posterior: chainThin, inc1990_used, tData_save, yObs_save
E = load('postSamples_erlang_clinical_AM_50.mat');   % Erlang posterior: chainThin, inc1990_used, tData_save, yObs_save

yearsObs = S.tData_save(:)';      % 2010..2023
meanInc  = S.yObs_save(:)';       % WHO mean
lowInc  = [114 138 141 141 171 181 191 193 185 180 172 173 170 164];
uppInc  = [508 440 410 395 329 300 263 242 233 227 220 228 231 228];

inc1990 = S.inc1990_used;

%% ---- Demography g(t) (Option A) ----
yr_known  = [1990 1995 2000 2005 2010 2015 2017 2018 2019 2020 2024]';
pop_known = [864972000; 960301000; 1057920000; 1154680000; 1243480000; ...
             1328020000; 1359660000; 1374660000; 1389030000; 1402620000; ...
             1450940000];
years_all  = (1990:2023)';                       
pop_annual = round(interp1(yr_known, pop_known, years_all, 'pchip'));
g_year     = log(pop_annual(2:end) ./ pop_annual(1:end-1));           
midpts     = (years_all(1:end-1) + years_all(2:end))/2;               
g_interp   = griddedInterpolant(midpts, g_year, 'nearest','nearest'); 
gfun = @(t_calendar) g_interp( min(max(t_calendar, midpts(1)), midpts(end)) );

%% ---- Draws & ACF settings ----
n_draws_simple = min(400, size(S.chainThin,1));
n_draws_erlang = min(400, size(E.chainThin,1));
idxS = randperm(size(S.chainThin,1), n_draws_simple);
idxE = randperm(size(E.chainThin,1), n_draws_erlang);

acf.yStart = 2024;
acf.yEnd   = 2026;
acf.mult_delta_up = 1.40;   % M->R
acf.mult_zeta_up  = 1.40;   % S->M

%% ---- Indices (Simple/Erlang param order) ----
% Simple: [btrans rhoS rhoM alpha beta gamma delta epsilon zeta eta theta muC muM muS nu]
% Erlang: [btrans rhoS rhoM alpha beta gamma delta epsilon kappa zeta theta muC muM muS nu]
idx = struct();
idx.beta  = 5;  
idx.delta = 7; 
idx.eta   = 10; 
idx.theta = 11; 
idx.muC   = 12;

%% ---- Time-slope gates ----
b_beta  = 0.00;
b_delta = 0.00;

%% ---- Quadrature ----
nQuad = 24;

%% ---- Erlang structure ----
k_erlang = 50;

%% ---- 1) Reconstruct 2010–2023 (three tracks) ----
% 1a) Erlang-from-Simple Posterior (κ*)
Yfit_eFromS = zeros(n_draws_simple, numel(yearsObs));
for j=1:n_draws_simple
    pS = S.chainThin(idxS(j),:);
    pE = map_simple_to_erlang_kappaStar(pS, k_erlang, idx);  % κ*
    [x2010, ~] = spin_from_1990_erlang(pE, inc1990, 2010, [], gfun, nQuad, b_beta, b_delta);
    [y2010_23, ~] = simulate_annual_erlang(pE, 2010, 2023, x2010, gfun, nQuad, b_beta, b_delta, []);
    Yfit_eFromS(j,:) = y2010_23;
end
fit_eFromS_mean = mean(Yfit_eFromS,1);

% 1b) Simple (original posterior)
Yfit_simple = zeros(n_draws_simple, numel(yearsObs));
for j=1:n_draws_simple
    p = S.chainThin(idxS(j),:);
    [x2010, ~] = spin_from_1990_simple(p, inc1990, 2010, [], gfun, nQuad, b_beta, b_delta);
    [y2010_23, ~] = simulate_annual_simple(p, 2010, 2023, x2010, gfun, nQuad, b_beta, b_delta, []);
    Yfit_simple(j,:) = y2010_23;
end
fit_simple_mean = mean(Yfit_simple,1);

% 1c) Erlang (original posterior)
Yfit_erlang = zeros(n_draws_erlang, numel(yearsObs));
for j=1:n_draws_erlang
    p = E.chainThin(idxE(j),:);
    [x2010, ~] = spin_from_1990_erlang(p, inc1990, 2010, [], gfun, nQuad, b_beta, b_delta);
    [y2010_23, ~] = simulate_annual_erlang(p, 2010, 2023, x2010, gfun, nQuad, b_beta, b_delta, []);
    Yfit_erlang(j,:) = y2010_23;
end
fit_erlang_mean = mean(Yfit_erlang,1);

%% ---- 2) Baseline forecast (2024–2033) (three tracks) ----
yearsF = 2024:2033;

% 2a) eFromS (κ*) baseline
Yf_eFromS  = zeros(n_draws_simple, numel(yearsF));
for j=1:n_draws_simple
    pS = S.chainThin(idxS(j),:);
    pE = map_simple_to_erlang_kappaStar(pS, k_erlang, idx);
    [x2010, ~] = spin_from_1990_erlang(pE, inc1990, 2010, [], gfun, nQuad, b_beta, b_delta);
    [~, x2024] = simulate_annual_erlang(pE, 2010, 2023, x2010, gfun, nQuad, b_beta, b_delta, []);
    [yf, ~]    = simulate_annual_erlang(pE, yearsF(1), yearsF(end), x2024, gfun, nQuad, b_beta, b_delta, []);
    Yf_eFromS(j,:) = yf;
end
fcast_eFromS_mean = mean(Yf_eFromS,1);

% 2b) Simple baseline
Yf_simple  = zeros(n_draws_simple, numel(yearsF));
for j=1:n_draws_simple
    p = S.chainThin(idxS(j),:);
    [x2010, ~] = spin_from_1990_simple(p, inc1990, 2010, [], gfun, nQuad, b_beta, b_delta);
    [~, x2024] = simulate_annual_simple(p, 2010, 2023, x2010, gfun, nQuad, b_beta, b_delta, []);
    [yf, ~]    = simulate_annual_simple(p, yearsF(1), yearsF(end), x2024, gfun, nQuad, b_beta, b_delta, []);
    Yf_simple(j,:) = yf;
end
fcast_simple_mean = mean(Yf_simple,1);

% 2c) Erlang baseline
Yf_erlang  = zeros(n_draws_erlang, numel(yearsF));
for j=1:n_draws_erlang
    p = E.chainThin(idxE(j),:);
    [x2010, ~] = spin_from_1990_erlang(p, inc1990, 2010, [], gfun, nQuad, b_beta, b_delta);
    [~, x2024] = simulate_annual_erlang(p, 2010, 2023, x2010, gfun, nQuad, b_beta, b_delta, []);
    [yf, ~]    = simulate_annual_erlang(p, yearsF(1), yearsF(end), x2024, gfun, nQuad, b_beta, b_delta, []);
    Yf_erlang(j,:) = yf;
end
fcast_erlang_mean = mean(Yf_erlang,1);

%% ---- 3) ACF scenario (2024–2026: zeta↑, delta↑) (three tracks) ----
% 3a) eFromS (κ*) ACF
YfE_acf = zeros(n_draws_simple, numel(yearsF));
for j=1:n_draws_simple
    pS = S.chainThin(idxS(j),:);
    pE = map_simple_to_erlang_kappaStar(pS, k_erlang, idx);
    [x2010, ~] = spin_from_1990_erlang(pE, inc1990, 2010, [], gfun, nQuad, b_beta, b_delta);
    [~, x2024] = simulate_annual_erlang(pE, 2010, 2023, x2010, gfun, nQuad, b_beta, b_delta, []);
    [yf, ~]    = simulate_annual_erlang(pE, yearsF(1), yearsF(end), x2024, gfun, nQuad, b_beta, b_delta, acf);
    YfE_acf(j,:) = yf;
end
fcastE_acf_mean = mean(YfE_acf,1);

% 3b) Simple ACF
YfS_acf = zeros(n_draws_simple, numel(yearsF));
for j=1:n_draws_simple
    p = S.chainThin(idxS(j),:);
    [x2010, ~] = spin_from_1990_simple(p, inc1990, 2010, [], gfun, nQuad, b_beta, b_delta);
    [~, x2024] = simulate_annual_simple(p, 2010, 2023, x2010, gfun, nQuad, b_beta, b_delta, []);
    [yf, ~]    = simulate_annual_simple(p, yearsF(1), yearsF(end), x2024, gfun, nQuad, b_beta, b_delta, acf);
    YfS_acf(j,:) = yf;
end
fcastS_acf_mean = mean(YfS_acf,1);

% 3c) Erlang ACF
YfEr_acf = zeros(n_draws_erlang, numel(yearsF));
for j=1:n_draws_erlang
    p = E.chainThin(idxE(j),:);
    [x2010, ~] = spin_from_1990_erlang(p, inc1990, 2010, [], gfun, nQuad, b_beta, b_delta);
    [~, x2024] = simulate_annual_erlang(p, 2010, 2023, x2010, gfun, nQuad, b_beta, b_delta, []);
    [yf, ~]    = simulate_annual_erlang(p, yearsF(1), yearsF(end), x2024, gfun, nQuad, b_beta, b_delta, acf);
    YfEr_acf(j,:) = yf;
end
fcastEr_acf_mean = mean(YfEr_acf,1);

%% ---- 3.5) 95% UI (ONLY for intervention; ONLY for 2024+) ----
fcastS_acf_lo  = prctile(YfS_acf,  2.5, 1);
fcastS_acf_hi  = prctile(YfS_acf, 97.5, 1);
fcastEr_acf_lo = prctile(YfEr_acf, 2.5, 1);
fcastEr_acf_hi = prctile(YfEr_acf,97.5, 1);

fcastS_acf_lo  = fcastS_acf_lo(:)';   fcastS_acf_hi  = fcastS_acf_hi(:)';
fcastEr_acf_lo = fcastEr_acf_lo(:)';  fcastEr_acf_hi = fcastEr_acf_hi(:)';

%% ---- 3.6) Summaries for all trajectories ----

% ---------- Historical fits: 2010–2023 ----------
fit_simple_lo   = prctile(Yfit_simple,   2.5, 1);
fit_simple_hi   = prctile(Yfit_simple,  97.5, 1);

fit_erlang_lo   = prctile(Yfit_erlang,   2.5, 1);
fit_erlang_hi   = prctile(Yfit_erlang,  97.5, 1);

fit_eFromS_lo   = prctile(Yfit_eFromS,   2.5, 1);
fit_eFromS_hi   = prctile(Yfit_eFromS,  97.5, 1);

% ---------- Forecast baseline: 2024–2033 ----------
fcast_simple_lo   = prctile(Yf_simple,   2.5, 1);
fcast_simple_hi   = prctile(Yf_simple,  97.5, 1);

fcast_erlang_lo   = prctile(Yf_erlang,   2.5, 1);
fcast_erlang_hi   = prctile(Yf_erlang,  97.5, 1);

fcast_eFromS_lo   = prctile(Yf_eFromS,   2.5, 1);
fcast_eFromS_hi   = prctile(Yf_eFromS,  97.5, 1);

% ---------- Forecast ACF: 2024–2033 ----------
fcastS_acf_lo     = prctile(YfS_acf,    2.5, 1);
fcastS_acf_hi     = prctile(YfS_acf,   97.5, 1);

fcastEr_acf_lo    = prctile(YfEr_acf,   2.5, 1);
fcastEr_acf_hi    = prctile(YfEr_acf,  97.5, 1);

fcastE_acf_lo     = prctile(YfE_acf,    2.5, 1);
fcastE_acf_hi     = prctile(YfE_acf,   97.5, 1);

%% ---- 3.7) Print to command window ----

fprintf('\n================ HISTORICAL FITS (2010–2023) ================\n');

fprintf('\n---- Simple fit ----\n');
fprintf('Year\tMean\tLo95\tHi95\n');
for i = 1:numel(yearsObs)
    fprintf('%d\t%.3f\t%.3f\t%.3f\n', yearsObs(i), fit_simple_mean(i), fit_simple_lo(i), fit_simple_hi(i));
end

fprintf('\n---- Erlang fit ----\n');
fprintf('Year\tMean\tLo95\tHi95\n');
for i = 1:numel(yearsObs)
    fprintf('%d\t%.3f\t%.3f\t%.3f\n', yearsObs(i), fit_erlang_mean(i), fit_erlang_lo(i), fit_erlang_hi(i));
end

fprintf('\n---- Erlang-from-Simple fit ----\n');
fprintf('Year\tMean\tLo95\tHi95\n');
for i = 1:numel(yearsObs)
    fprintf('%d\t%.3f\t%.3f\t%.3f\n', yearsObs(i), fit_eFromS_mean(i), fit_eFromS_lo(i), fit_eFromS_hi(i));
end

fprintf('\n================ FORECAST BASELINE (2024–2033) ================\n');

fprintf('\n---- Simple baseline ----\n');
fprintf('Year\tMean\tLo95\tHi95\n');
for i = 1:numel(yearsF)
    fprintf('%d\t%.3f\t%.3f\t%.3f\n', yearsF(i), fcast_simple_mean(i), fcast_simple_lo(i), fcast_simple_hi(i));
end

fprintf('\n---- Erlang baseline ----\n');
fprintf('Year\tMean\tLo95\tHi95\n');
for i = 1:numel(yearsF)
    fprintf('%d\t%.3f\t%.3f\t%.3f\n', yearsF(i), fcast_erlang_mean(i), fcast_erlang_lo(i), fcast_erlang_hi(i));
end

fprintf('\n---- Erlang-from-Simple baseline ----\n');
fprintf('Year\tMean\tLo95\tHi95\n');
for i = 1:numel(yearsF)
    fprintf('%d\t%.3f\t%.3f\t%.3f\n', yearsF(i), fcast_eFromS_mean(i), fcast_eFromS_lo(i), fcast_eFromS_hi(i));
end

fprintf('\n================ FORECAST ACF (2024–2033) ================\n');

fprintf('\n---- Simple ACF (2024–26) ----\n');
fprintf('Year\tMean\tLo95\tHi95\n');
for i = 1:numel(yearsF)
    fprintf('%d\t%.3f\t%.3f\t%.3f\n', yearsF(i), fcastS_acf_mean(i), fcastS_acf_lo(i), fcastS_acf_hi(i));
end

fprintf('\n---- Erlang ACF (2024–26) ----\n');
fprintf('Year\tMean\tLo95\tHi95\n');
for i = 1:numel(yearsF)
    fprintf('%d\t%.3f\t%.3f\t%.3f\n', yearsF(i), fcastEr_acf_mean(i), fcastEr_acf_lo(i), fcastEr_acf_hi(i));
end

fprintf('\n---- Erlang-from-Simple ACF (2024–26) ----\n');
fprintf('Year\tMean\tLo95\tHi95\n');
for i = 1:numel(yearsF)
    fprintf('%d\t%.3f\t%.3f\t%.3f\n', yearsF(i), fcastE_acf_mean(i), fcastE_acf_lo(i), fcastE_acf_hi(i));
end

%% ---- 4) Plot (6 lines total) ----
figure('Position',[120 120 1100 620]); hold on; grid on; box on;

% Colors 
simpleColor = [0.00 0.35 0.90];   
erlangColor = [0.90 0.10 0.10];   
violet      = [0.55 0.00 0.85];   
obsBandColor = [0.78 0.86 0.62];
obsMeanColor = [0.20 0.50 0.10];

xF = yearsF;
bandAlpha = 0.16;

hUI_simpleACF = fill([xF fliplr(xF)], [fcastS_acf_lo fliplr(fcastS_acf_hi)], simpleColor, ...
    'EdgeColor','none','FaceAlpha',bandAlpha, 'DisplayName','Simple ACF 95% UI');
hUI_erlangACF = fill([xF fliplr(xF)], [fcastEr_acf_lo fliplr(fcastEr_acf_hi)], erlangColor, ...
    'EdgeColor','none','FaceAlpha',bandAlpha, 'DisplayName','Gamma Distributed ACF 95% UI');
uistack(hUI_simpleACF,'bottom');
uistack(hUI_erlangACF,'bottom');
% ===== END UI =====

% Observations
xHist = yearsObs;
fill([xHist fliplr(xHist)], [lowInc fliplr(uppInc)], obsBandColor, ...
    'EdgeColor','none','FaceAlpha',0.45, 'DisplayName','Observed 95% UI');
plot(xHist, meanInc, 'o-', 'Color',obsMeanColor, 'LineWidth',3, ...
     'MarkerFaceColor',obsMeanColor, 'DisplayName','Observed mean');

% Historical fits (hidden from legend)
plot(xHist, fit_simple_mean, '-', 'Color', simpleColor, 'LineWidth',3, 'HandleVisibility','off');
plot(xHist, fit_erlang_mean,  '-', 'Color', erlangColor, 'LineWidth',3, 'HandleVisibility','off');
plot(xHist, fit_eFromS_mean,  '-', 'Color', violet,      'LineWidth',3, 'HandleVisibility','off');

% Baselines (3)
plot(xF, fcast_simple_mean, '-', 'Color', simpleColor, 'LineWidth',3, 'DisplayName','Simple baseline');
plot(xF, fcast_erlang_mean, '-', 'Color', erlangColor, 'LineWidth',3, 'DisplayName','Gamma Distributed baseline');
plot(xF, fcast_eFromS_mean, '-', 'Color', violet,      'LineWidth',3, 'DisplayName','Gamma simple-posterior baseline');

% ACF (3)
plot(xF, fcastS_acf_mean,  '--', 'Color', simpleColor, 'LineWidth',3, 'DisplayName','Simple (ACF 2024–26)');
plot(xF, fcastEr_acf_mean, '--', 'Color', erlangColor, 'LineWidth',3, 'DisplayName','Gamm Distributed (ACF 2024–26)');
plot(xF, fcastE_acf_mean,  '--', 'Color', violet,      'LineWidth',3, 'DisplayName','Gamma simple-posterior (ACF 2024–26)');

% Connectors (hidden)
plot([2023 2024],[fit_simple_mean(end)   fcast_simple_mean(1)], '-',  'Color',simpleColor,'LineWidth',2,'HandleVisibility','off');
plot([2023 2024],[fit_erlang_mean(end)   fcast_erlang_mean(1)], '-',  'Color',erlangColor,'LineWidth',2,'HandleVisibility','off');
plot([2023 2024],[fit_eFromS_mean(end)   fcast_eFromS_mean(1)], '-',  'Color',violet,     'LineWidth',2,'HandleVisibility','off');
plot([2023 2024],[fit_simple_mean(end)   fcastS_acf_mean(1)],  '--',  'Color',simpleColor,'LineWidth',2,'HandleVisibility','off');
plot([2023 2024],[fit_erlang_mean(end)   fcastEr_acf_mean(1)], '--',  'Color',erlangColor,'LineWidth',2,'HandleVisibility','off');
plot([2023 2024],[fit_eFromS_mean(end)   fcastE_acf_mean(1)],  '--',  'Color',violet,     'LineWidth',2,'HandleVisibility','off');

xlabel('Year'); ylabel('Incidence per 100 000');
title('TB annual incidence in India');
xlim([2010, 2033]); 
vmax = @(x) max(x(:));
yMax = 1.10 * max([vmax(uppInc), vmax(meanInc), ...
                   vmax(fit_simple_mean), vmax(fit_erlang_mean), vmax(fit_eFromS_mean), ...
                   vmax(fcast_simple_mean), vmax(fcast_erlang_mean), vmax(fcast_eFromS_mean), ...
                   vmax(fcastS_acf_mean), vmax(fcastEr_acf_mean), vmax(fcastE_acf_mean)]);
ylim([0, yMax]);
legend('Location','northeast');
set(gca,'FontSize',12,'LineWidth',1);

%% ================== Helpers ==================
% ---------- year params ----------
function p_year = year_params_simple(p, year, b_beta, b_delta, acf)
    % Simple: [btr rhoS rhoM alpha beta gamma delta epsilon zeta eta theta muC muM muS nu]
    p_year = p;
    if year>=2010
        p_year(5) = p(5) * exp(b_beta  * (year-2010));   % beta
        p_year(7) = p(7) * exp(b_delta * (year-2010));   % delta
    end
    if nargin>=5 && ~isempty(acf)
        if year>=acf.yStart && year<=acf.yEnd
            p_year(7) = p_year(7) * acf.mult_delta_up;   % delta
            p_year(9) = p_year(9) * acf.mult_zeta_up;    % zeta
        end
    end
end

function p_year = year_params_erlang(p, year, b_beta, b_delta, acf)
    % Erlang: [btr rhoS rhoM alpha beta gamma delta epsilon kappa zeta theta muC muM muS nu]
    p_year = p;
    if year>=2010
        p_year(5) = p(5) * exp(b_beta  * (year-2010));   % beta
        p_year(7) = p(7) * exp(b_delta * (year-2010));   % delta
    end
    if nargin>=5 && ~isempty(acf)
        if year>=acf.yStart && year<=acf.yEnd
            p_year(7)  = p_year(7)  * acf.mult_delta_up; % delta
            p_year(10) = p_year(10) * acf.mult_zeta_up;  % zeta
        end
    end
end

% ---------- spin from 1990 ----------
function [x_at_target, yAnnual, years_list] = spin_from_1990_simple(p, inc1990, targetYear, recordYears, gfun, nQuad, b_beta, b_delta)
    theta0 = p(11); muC = p(12); muM = p(13); muS = p(14); nu = p(15);
    eta0 = p(10); zeta0 = p(9); delta0 = p(7); eps0 = p(8);
    C0 = (inc1990 / max(theta0 + nu + muC, 1e-9)) / 1e5;
    I0 = 1e-3;
  
    Ssub = max(((theta0 + nu + muC)/max(eta0,1e-9)) * C0, 0);
    M0   = max(zeta0*Ssub / max(delta0 + eps0 + nu + muM, 1e-9), 0);
    R0 = 0; D0 = 0; 
    S0pool = max(1-(I0+M0+Ssub+C0+R0), 1e-6);
    x_now = [S0pool; I0; M0; Ssub; C0; R0; D0];

    odef    = @(t,x,pp,aux) simpleTBmodel_demog(t,x,pp,aux);
    inc_fun = @(X,pp) pp(10)*X(4,:);  % eta * S

    years_list = recordYears(:)'; 
    collect = ~isempty(years_list);
    t_now = 0; 
    if collect, yAnnual = nan(1,numel(years_list)); else, yAnnual = []; end

    for yy = 1990:(targetYear-1)
        p_year = year_params_simple(p, yy, b_beta, b_delta, []);
        aux.gfun = gfun; aux.t0 = 1990;
        f = @(t,x) odef(t,x,p_year,aux);
        y_start = yy - 1990; 
        y_end   = y_start + 1;

        if t_now < y_start
            sol1 = ode15s(f, [t_now y_start], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
            x_now = deval(sol1, y_start); 
            t_now = y_start;
        end

        sol2 = ode15s(f, [y_start y_end], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
        tt = linspace(y_start, y_end, nQuad);
        Xs = deval(sol2, tt);

        if collect
            pos = find(years_list==yy, 1);
            if ~isempty(pos)
                yAnnual(pos) = trapz(tt, inc_fun(Xs, p_year)) * 1e5;
            end
        end

        x_now = deval(sol2, y_end); 
        t_now = y_end;
    end
    x_at_target = x_now;
end

function [x_at_target, yAnnual, years_list] = spin_from_1990_erlang(p, inc1990, targetYear, recordYears, gfun, nQuad, b_beta, b_delta)
    k = 50;
    theta = p(11); muC = p(12); muM = p(13); muS = p(14); nu = p(15);
    kappa = p(9);  zeta = p(10);

    C0 = (inc1990 / max(theta + nu + muC, 1e-9)) / 1e5;
    I0 = 1e-3;

    Sk = ((theta + nu + muC) / max(kappa,1e-9)) * C0;
    L  = (kappa + zeta + nu + muS);
    r  = kappa / max(L,1e-12);
    Svec = zeros(k,1);
    Svec(k) = Sk;
    for jj = k-1:-1:1
        Svec(jj) = r^(k-jj) * Sk;
    end

    delta0 = p(7); eps0 = p(8);
    M0 = max(zeta*sum(Svec) / max(delta0 + eps0 + nu + muM, 1e-9), 0);

    R0 = 0; D0 = 0;
    S0pool = max(1-(I0 + M0 + sum(Svec) + C0 + R0), 1e-6);

    x_now = [S0pool; I0; M0; Svec; C0; R0; D0];

    odef    = @(t,x,pp,aux) erlangODE_demog(t,x,pp,aux,k);
    inc_fun = @(X,pp) pp(9)*X(3+k,:);   % κ * S_k

    years_list = recordYears(:)'; 
    collect = ~isempty(years_list);
    t_now = 0; 
    if collect, yAnnual = nan(1,numel(years_list)); else, yAnnual = []; end

    for yy = 1990:(targetYear-1)
        p_year = year_params_erlang(p, yy, b_beta, b_delta, []);
        aux.gfun = gfun; aux.t0 = 1990;
        f = @(t,x) odef(t,x,p_year,aux);
        y_start = yy - 1990; 
        y_end   = y_start + 1;

        if t_now < y_start
            sol1 = ode15s(f, [t_now y_start], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
            x_now = deval(sol1, y_start); 
            t_now = y_start;
        end

        sol2 = ode15s(f, [y_start y_end], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
        tt = linspace(y_start, y_end, nQuad);
        Xs = deval(sol2, tt);

        if collect
            pos = find(years_list==yy, 1);
            if ~isempty(pos)
                yAnnual(pos) = trapz(tt, inc_fun(Xs, p_year)) * 1e5;
            end
        end

        x_now = deval(sol2, y_end); 
        t_now = y_end;
    end
    x_at_target = x_now;
end

% ---------- simulate annual ----------
function [yAnnual, x_after] = simulate_annual_simple(p, startYear, endYear, x_start, gfun, nQuad, b_beta, b_delta, acf)
    odef    = @(t,x,pp,aux) simpleTBmodel_demog(t,x,pp,aux);
    inc_fun = @(X,pp) pp(10)*X(4,:);  % eta * S

    yAnnual = nan(1, endYear-startYear+1);
    x_now = x_start; t_now = 0; baseYear = startYear;

    for yr = startYear:endYear
        p_year = year_params_simple(p, yr, b_beta, b_delta, acf);
        aux.gfun = gfun; aux.t0 = baseYear;
        f = @(t,x) odef(t,x,p_year,aux);

        y_start = (yr-startYear);
        y_end   = y_start + 1;

        if t_now < y_start
            sol1 = ode15s(f, [t_now y_start], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
            x_now = deval(sol1, y_start); 
            t_now = y_start;
        end

        sol2 = ode15s(f, [y_start y_end], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
        tt = linspace(y_start, y_end, nQuad);
        Xs = deval(sol2, tt);

        yAnnual(yr-startYear+1) = trapz(tt, inc_fun(Xs, p_year)) * 1e5;

        x_now = deval(sol2, y_end); 
        t_now = y_end;
    end
    x_after = x_now;
end

function [yAnnual, x_after] = simulate_annual_erlang(p, startYear, endYear, x_start, gfun, nQuad, b_beta, b_delta, acf)
    k = 50;
    odef    = @(t,x,pp,aux) erlangODE_demog(t,x,pp,aux,k);
    inc_fun = @(X,pp) pp(9)*X(3+k,:);   % κ * S_k

    yAnnual = nan(1, endYear-startYear+1);
    x_now = x_start; t_now = 0; baseYear = startYear;

    for yr = startYear:endYear
        p_year = year_params_erlang(p, yr, b_beta, b_delta, acf);
        aux.gfun = gfun; aux.t0 = baseYear;
        f = @(t,x) odef(t,x,p_year,aux);

        y_start = (yr-startYear); 
        y_end   = y_start + 1;

        if t_now < y_start
            sol1 = ode15s(f, [t_now y_start], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
            x_now = deval(sol1, y_start); 
            t_now = y_start;
        end

        sol2 = ode15s(f, [y_start y_end], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
        tt = linspace(y_start, y_end, nQuad);
        Xs = deval(sol2, tt);

        yAnnual(yr-startYear+1) = trapz(tt, inc_fun(Xs, p_year)) * 1e5;

        x_now = deval(sol2, y_end); 
        t_now = y_end;
    end
    x_after = x_now;
end

% ---------- ODEs ----------
function dx = simpleTBmodel_demog(t, x, p, aux)
% p = [btrans rhoS rhoM alpha beta gamma delta epsilon zeta eta theta muC muM muS nu]
btr = p(1); rhoS = p(2); rhoM = p(3);
alpha = p(4); beta = p(5); gamma = p(6); delta = p(7); epsilon = p(8);
zeta = p(9); eta = p(10); theta = p(11); muC = p(12); muM = p(13); muS = p(14); nu = p(15);
S0 = x(1); I = x(2); M = x(3); S = x(4); C = x(5); R = x(6); D = x(7); %#ok<NASGU>
Nlive = max(S0 + I + M + S + C + R, 1e-12);
g = aux.gfun(aux.t0 + t);
Pi = (g + nu) * Nlive;                     
lambda = btr * (C + rhoS*S + rhoM*M) / Nlive;
dS0 = Pi - lambda*S0 - nu*S0;
dI  = lambda*S0 - (alpha + beta + gamma + nu)*I;
dM  = beta*I - (delta + epsilon + nu + muM)*M + zeta*S;           
dS  = gamma*I + epsilon*M - (zeta + eta + nu + muS)*S + theta*C;
dC  = eta*S - (theta + nu + muC)*C;
dR  = alpha*I + delta*M - nu*R;                                    
dD  = nu*(S0 + I + M + S + C + R) + muC*C + muM*M + muS*S;
dx = [dS0; dI; dM; dS; dC; dR; dD];
end

function dx = erlangODE_demog(t, x, p, aux, k)
% p = [btrans rhoS rhoM alpha beta gamma delta epsilon kappa zeta theta muC muM muS nu]
btr = p(1); rhoS = p(2); rhoM = p(3);
alpha = p(4); beta = p(5); gamma = p(6); delta = p(7); epsilon = p(8);
kappa = p(9); zeta = p(10); theta = p(11); muC = p(12); muM = p(13); muS = p(14); nu = p(15);

S0 = x(1); I = x(2); M = x(3); S = x(4:3+k); C = x(3+k+1); R = x(3+k+2); D = x(3+k+3); %#ok<NASGU>
Nlive = max(S0 + I + M + sum(S) + C + R, 1e-12);
g = aux.gfun(aux.t0 + t);
Pi = (g + nu) * Nlive;
lambda = btr * (C + rhoS*sum(S) + rhoM*M) / Nlive;

dS0 = Pi - lambda*S0 - nu*S0;
dI  = lambda*S0 - (alpha + beta + gamma + nu)*I;
dM  = beta*I - (delta + epsilon + nu + muM)*M + zeta*sum(S);

dS  = zeros(k,1);
dS(1) = gamma*I + epsilon*M - (kappa + zeta + nu + muS)*S(1) + theta*C;
for j=2:k
    dS(j) = kappa*S(j-1) - (kappa + zeta + nu + muS)*S(j);
end

dC = kappa*S(k) - (theta + nu + muC)*C;
dR = alpha*I + delta*M - nu*R;
dD = nu*(S0 + I + M + sum(S) + C + R) + muC*C + muM*M + muS*sum(S);

dx = [dS0; dI; dM; dS; dC; dR; dD];
end

% ---------- Mapping κ* ----------
function pE = map_simple_to_erlang_kappaStar(pS, k, idx)
%   ((κ*/(κ* + L))^k) = η/(η + L), L = ζ + ν + μS
%   => r = (η/(η+L))^(1/k), κ* = L * r / (1 - r)
pE = zeros(1,15);

pE(1)  = pS(1);   % btrans
pE(2)  = pS(2);   % rhoS
pE(3)  = pS(3);   % rhoM
pE(4)  = pS(4);   % alpha
pE(5)  = pS(5);   % beta
pE(6)  = pS(6);   % gamma
pE(7)  = pS(7);   % delta
pE(8)  = pS(8);   % epsilon
pE(10) = pS(9);   % zeta
pE(11) = pS(idx.theta); % theta
pE(12) = pS(idx.muC);   % muC
pE(13) = pS(13);  % muM
pE(14) = pS(14);  % muS
pE(15) = pS(15);  % nu
% κ*
eta  = pS(idx.eta);
zeta = pS(9);
muS  = pS(14);
nu   = pS(15);
L    = zeta + nu + muS;
r    = (eta / max(eta + L, 1e-12))^(1/k);
kappa_star = L * r / max(1 - r, 1e-12);
pE(9) = kappa_star;
end
