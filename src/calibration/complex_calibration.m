% ====================== calibrate_who_india_fitonly_likeV2.m =======================
% Adaptive-Metropolis calibration of Simple and Erlang TB models
% to WHO India clinical incidence (per 100k) 2010–2023, NO forecast.
%
% Demography update (Option A births):
%   • Add S0 (susceptible stock) and demographic flows.
%   • Births Π_A(t) = (g(t)+ν) N(t), where g(t) is from WHO India population.
%   • Natural mortality ν acts on all living states; TB-excess deaths on C/M/S.
%   • Force of infection λ(t) = β_trans * ( C + ρS·ΣS_j + ρM·M ) / N(t).
%
% Changes per request (single 1990 incidence to try first):
%   • Fix the 1990 initial incidence at a chosen value (default 300 per 100k),
%     fit 2010–2023, and inspect the calibration outcome. No multi-candidate
%     search or selection — run exactly this initial value.
%   • Slope-gating kept ONLY for beta & delta as a simple *linear-on-log*
%     trend from 2010 onward. No COVID ramp; pre-2010 parameters are constant.
%   • Plots show 1990–2009 model history in a different color from 2010–2023.
%
% Unchanged from v2-like logic:
%   • Annual INCIDENCE likelihood by integrating (eta*S) or (kappa*S_k).
%   • Initial state mapping (quasi-steady via C0 = inc/(theta+muC), etc.).
%   • WHO-style calibration look/feel and saved filenames/fields.
% -----------------------------------------------------------------------------------
clear; clc; close all; rng(20250609);

%% Interval-width knob (likelihood only)
sigma_add_like = 10;   % per 100k; larger => wider intervals via the likelihood

%% Choose ONE 1990 incidence to run (per 100k)
inc1990_fixed = 2000;  % default; may be auto-overridden below

%% A) WHO India clinical data (hard-coded)
yearsObs = (2010:2023)';               % observed calendar years
yearsPre = (1990:2009)';               % model prehistory (no likelihood)

% WHO mean and 95% UI (per 100k)
meanInc = [276 268 258 249 243 237 225 217 208 202 195 200 199 195]'; % mean per 100k
lowInc  = [114 138 141 141 171 181 191 193 185 180 172 173 170 164]'; % lower
uppInc  = [508 440 410 395 329 300 263 242 233 227 220 228 231 228]'; % upper

% Observation SD from 95% CI (normal approximation)
sdObs = max((uppInc - lowInc) / 3.92, 1e-3);

% Time axis in "years since 2010" for likelihood integration
tData = yearsObs - yearsObs(1);        % 0..13 for 2010..2023

%% B) WHO India population (manual anchors -> 1990–2023 annual g(t))
% Anchors from your screenshot (persons, no commas)
yr_known  = [1990 1995 2000 2005 2010 2015 2017 2018 2019 2020 2024]';
pop_known = [864972000; 960301000; 1057920000; 1154680000; 1243480000; ...
             1328020000; 1359660000; 1374660000; 1389030000; 1402620000; ...
             1450940000];

% Build annual 1990–2023 population by pchip; then continuous annual growth g_y
years_all = (1990:2023)';
pop_annual = round(interp1(yr_known, pop_known, years_all, 'pchip'));
g_year = log(pop_annual(2:end) ./ pop_annual(1:end-1));    % 1990..2022
midpts = (years_all(1:end-1) + years_all(2:end))/2;        % year midpoints
g_interpolant = griddedInterpolant(midpts, g_year, 'nearest','nearest'); % piecewise-const
gfun = @(t_calendar) g_interpolant(t_calendar);            % calendar year -> per-year growth

%% C) Priors / bounds (updated with new parameters)
%（Simple）:
% pS = [btrans, rhoS, rhoM, alpha, beta, gamma, delta, epsilon, zeta, eta, theta, muC, muM, muS, nu]
%
%（Erlang）:
% pE = [btrans, rhoS, rhoM, alpha, beta, gamma, delta, epsilon, kappa, zeta, theta, muC, muM, muS, nu]
LB_simple = @() [0.05 0.00 0.00 0.15 0.08 0.08 0.08 0.08 0.08 0.08 0.08 0.10 0.00 0.00 0.005];
UB_simple = @() [4.00 1.00 1.00 5.85 2.93 2.93 2.93 2.93 2.93 2.93 2.93 1.00 0.20 0.20 0.030];

k_erlang = 5;  % Erlang shape
LB_erlang_local = @(k) [0.05  0.00 0.00   0.15 0.08 0.08 0.08 0.08   3*k/2 0.08 0.08   0.10 0.00 0.00 0.005];
UB_erlang_local = @(k) [4.00  1.00 1.00   5.85 2.93 2.93 2.93 2.93   4*k   2.93 2.93   1.00 0.20 0.20 0.030];
LB_erlang = LB_erlang_local(k_erlang);
UB_erlang = UB_erlang_local(k_erlang);

% Logistic transform helpers
logistic = @(u) 1./(1+exp(-u));
u2p     = @(u,lb,ub) lb + (ub-lb).*logistic(u);
logJac  = @(u,lb,ub) sum(log(ub-lb) + log(logistic(u)) + log(1-logistic(u)));

%% === Plan A ===
use_auto_inc1990 = true;      
if use_auto_inc1990
    fprintf('\n[Auto-pick inc1990] Finding a reasonable 1990 incidence by anchoring 2010 level...\n');
    [pMAP_erlang, okMAP] = AUTO_get_MAP_for_erlang( ...
        LB_erlang, UB_erlang, k_erlang, tData, meanInc, sdObs, ...
        sigma_add_like, yearsPre, gfun);
    if ~okMAP
        warning('Auto-pick inc1990 (Erlang): MAP warm start failed; using mid-point params.');
        pMAP_erlang = (LB_erlang + UB_erlang)/2;
    end

    y2010_obs = meanInc(1);
    obj = @(inc90) ( AUTO_simulate_2010_inc_from_inc1990( ...
                        pMAP_erlang, k_erlang, inc90, yearsPre, gfun, ...
                        @(t,x,p,aux) erlangODE_demog(t,x,p,aux,k_erlang)) ...
                 - y2010_obs ).^2;

    lo = 100; hi = 4000;
    opts = optimset('Display','off');
    [inc_best, ~] = fminbnd(obj, lo, hi, opts);
    inc1990_fixed = max(50, min(5000, inc_best));
    fprintf('[Auto-pick inc1990] Using inc1990_fixed = %.0f per 100k\n', inc1990_fixed);
else
    fprintf('[Auto-pick inc1990] Disabled. Using manual inc1990_fixed = %.0f\n', inc1990_fixed);
end

%% D) Adaptive Metropolis settings (same style)
nIter = 2e5; burnIn = 5e4; thin = 10;
% s2 tuned by parameter dimension inside runner

%% E) Fit SIMPLE model (k=1) at the fixed 1990 incidence
fprintf('\n=== Calibrating SIMPLE+Demography (Option A) @ inc1990=%g ===\n', inc1990_fixed);
[chainThin_simple, Yhist_simple, Ypre_simple, ~] = run_AM_hist_likeV2( ...
    @simpleTBmodel_demog, LB_simple(), UB_simple(), 1, tData, meanInc, sdObs, ...
    nIter, burnIn, thin, u2p, logJac, sigma_add_like, inc1990_fixed, yearsPre, gfun);

%% F) Fit ERLANG model (k=2) at the fixed 1990 incidence
fprintf('\n=== Calibrating ERLANG(k=%d)+Demography (Option A) @ inc1990=%g ===\n', k_erlang, inc1990_fixed);
erlang_handle = @(t,x,p,aux) erlangODE_demog(t,x,p,aux,k_erlang);
[chainThin_erlang, Yhist_erlang, Ypre_erlang, ~] = run_AM_hist_likeV2( ...
    erlang_handle, LB_erlang, UB_erlang, k_erlang, tData, meanInc, sdObs, ...
    nIter, burnIn, thin, u2p, logJac, sigma_add_like, inc1990_fixed, yearsPre, gfun);

%% G) WHO-style calibration plots (v2 colors/labels/font) + prehistory
make_calibration_plot_with_prehistory(yearsPre, Ypre_simple, ...
    yearsObs, meanInc, lowInc, uppInc, Yhist_simple, ...
    sprintf('Calibration (1990–2023) — Simple+Demog (incidence); inc1990=%g', inc1990_fixed));

make_calibration_plot_with_prehistory(yearsPre, Ypre_erlang, ...
    yearsObs, meanInc, lowInc, uppInc, Yhist_erlang, ...
    sprintf('Calibration (1990–2023) — Erlang(k=%d)+Demog (incidence); inc1990=%g', k_erlang, inc1990_fixed));

%% H) Save — keep v2 filenames/fields; put empty forecast placeholders
tData_save = yearsObs;           % calendar years (v2 style)
yObs_save  = meanInc;            % incidence per 100k (likelihood scale)
x0_list    = {[]};
tPred      = [];                 % no forecast
Ypred      = [];                 % no forecast

% SIMPLE
chainThin = chainThin_simple; k = 1; inc1990_used = inc1990_fixed; %#ok<NASGU>
save('postSamples_simple_clinical_AM.mat', ...
     'chainThin','tData_save','k','yObs_save','x0_list','tPred','Ypred','inc1990_used');

% ERLANG
chainThin = chainThin_erlang; k = k_erlang; inc1990_used = inc1990_fixed; %#ok<NASGU>
save('postSamples_erlang_clinical_AM.mat', ...
     'chainThin','tData_save','k','yObs_save','x0_list','tPred','Ypred','inc1990_used');

% Legacy filename (Erlang content) — kept for compatibility
save('postSamples_subclinical_uniform_adaptcov.mat', ...
     'chainThin','tData_save','k','yObs_save','x0_list','tPred','Ypred','inc1990_used');

fprintf('\nSaved:\n postSamples_simple_clinical_AM.mat\n postSamples_erlang_clinical_AM.mat\n postSamples_subclinical_uniform_adaptcov.mat\n');



%% ====================== core: v2-like fit, history only ======================
function [chainThin, Yhist, YpreHist, mapLogPost] = run_AM_hist_likeV2( ...
    odefun, LB, UB, k, tData, yObs, sdObs, ...
    nIter, burnIn, thin, u2p, logJac, sigma_add_like, inc1990_per100k, yearsPre, gfun)

nParam = numel(LB);
s2 = 0.7*(2.38^2) / nParam;   % dimension-safe AM scale
epsilon = 1e-5;                % nugget for adaptive covariance
adaptIntv = 1000;              % more frequent selfadapt, improve burn-in efficient
nQuad = 24;                    % annual integral evaluation points
% ---- Linear-from-2010 slope gating (ONLY beta & delta) ----
% log theta_y = log theta0 + b * (year - 2010) * 1_{year >= 2010}
b_beta  = -0.00;   % slow decline for beta
b_delta = +0.00;   % slow improvement for delta

if k==1
    idx_beta = 5; idx_delta = 7; idx_eta = 10; % simple indices
else
    idx_beta = 5; idx_delta = 7; idx_kappa = 9; % erlang indices
end

beta_piecewise   = @(beta0,year)  beta0  .* exp(b_beta  * max(year-2010,0));
delta_piecewise  = @(delta0,year) delta0 .* exp(b_delta * max(year-2010,0));

% ---- v2-style initial state mapping at an arbitrary calendar year ----
    function x0_here = make_x0_from_params_v2(p, k_local, inc_per100k, year0)
        % Compute a quasi-steady C0 from incidence = (theta+muC)*C
        theta0 = p(11); muC = p(12);
        C0 = (inc_per100k / max(theta0 + muC, 1e-9)) / 1e5;  % fraction
        I0 = 1e-3;  % tiny seed

        if k_local==1
            eta0   = p(idx_eta);
            delta0 = p(7);
            zeta0  = p(9);
            eps0   = p(8);
            % S and M from steady-like balances
            S0sub  = max( (eta0>0) * ((theta0+muC)/eta0) * C0, 0); % inverse of incidence balance
            M0     = max( zeta0*S0sub / max(delta0 + eps0,1e-9), 0);
            R0     = 0;
            D0     = 0;
            % Allocate S0 (susceptible pool) as the rest of living mass ~ 1
            S0pool = max(1 - (I0 + M0 + S0sub + C0 + R0), 1e-6);
            x0_here = [S0pool; I0; M0; S0sub; C0; R0; D0];
        else
            kappa0 = p(9); zeta  = p(10); delta0 = p(7);
            Svec = zeros(k_local,1);
            % set Sk from incidence balance: (theta+muC)*C = kappa*Sk  => Sk = (theta+muC)/kappa * C
            Sk = ((theta0+muC)/max(kappa0,1e-9)) * C0;
            Svec(k_local) = Sk;
            for jj=k_local-1:-1:1
                Svec(jj) = (kappa0/max(kappa0+zeta,1e-9))^(k_local-jj) * Sk;
            end
            M0  = max( zeta*sum(Svec) / max(delta0 + p(8),1e-9), 0);
            R0  = 0; D0 = 0;
            I0  = 1e-3;
            S0pool = max(1 - (I0 + M0 + sum(Svec) + C0 + R0), 1e-6);
            x0_here = [S0pool; I0; M0; Svec; C0; R0; D0];
        end
    end

% ---- spin-up from 1990 to 2010; also return 1990–2009 annual incidence ----
    function [x2010, yPre_row] = spinup_to_2010(p, k_local, inc1990)
        % initial state at 1990 using mapping with year0=1990
        x_now = make_x0_from_params_v2(p, k_local, inc1990, 1990);
        t_now = 0;   % time in "years since 1990"
        yPre = nan(1, numel(yearsPre));
        for ii = 1:numel(yearsPre)     % 1990..2009
            yCal = yearsPre(ii);
            if k_local==1
                beta_y  = beta_piecewise( p(idx_beta),  yCal);
                delta_y = delta_piecewise(p(idx_delta), yCal);
                eta_y   = p(idx_eta);
                p_year = p; p_year(idx_beta)=beta_y; p_year(idx_delta)=delta_y; % eta const
                aux.gfun = gfun; aux.t0 = 1990;
                ode_year = @(t,x) odefun(t,x,p_year,aux);
                y_start = (yCal - 1990); y_end = y_start + 1;
                if t_now < y_start
                    sol1 = ode15s(ode_year, [t_now y_start], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
                    x_now = deval(sol1, y_start); t_now = y_start;
                end
                sol2 = ode15s(ode_year, [y_start y_end], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
                tt   = linspace(y_start, y_end, nQuad);
                Xs   = deval(sol2, tt);
                Ssub = Xs(4,:); inc_rt = eta_y.*Ssub;
                yPre(ii) = trapz(tt, inc_rt) * 1e5;
                x_now = deval(sol2, y_end); t_now = y_end;
            else
                beta_y  = beta_piecewise( p(idx_beta),  yCal);
                delta_y = delta_piecewise(p(idx_delta), yCal);
                kappa_y = p(idx_kappa);
                p_year = p; p_year(idx_beta)=beta_y; p_year(idx_delta)=delta_y; p_year(idx_kappa)=kappa_y;
                aux.gfun = gfun; aux.t0 = 1990;
                ode_year = @(t,x) odefun(t,x,p_year,aux);
                y_start = (yCal - 1990); y_end = y_start + 1;
                if t_now < y_start
                    sol1 = ode15s(ode_year, [t_now y_start], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
                    x_now = deval(sol1, y_start); t_now = y_start;
                end
                sol2 = ode15s(ode_year, [y_start y_end], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
                tt   = linspace(y_start, y_end, nQuad);
                Xs   = deval(sol2, tt);
                kLoc = k_local; Sk = Xs(3+kLoc,:); inc_rt = kappa_y.*Sk;
                yPre(ii) = trapz(tt, inc_rt) * 1e5;
                x_now = deval(sol2, y_end); t_now = y_end;
            end
        end
        x2010 = x_now;                         % state at start of 2010
        yPre_row = reshape(yPre,1,[]);
    end

% ---- annual-incidence integrator for 2010–2023 (returns row vector 1×T) ----
    function yInc_row = simulate_incidence_annual_from2010(times_year_idx, p, x0_2010)
        if isempty(times_year_idx), yInc_row = []; return; end
        ty = times_year_idx(:)'; x_now = x0_2010; t_now = 0; 
        yInc = nan(1, numel(ty));
        for ii = 1:numel(ty)                     % years since 2010
            y_start = ty(ii); y_end = y_start + 1;
            yCal = 2010 + y_start;               % 2010..2023
            if k==1
                beta_y  = beta_piecewise( p(idx_beta),  yCal);
                delta_y = delta_piecewise(p(idx_delta), yCal);
                eta_y   = p(idx_eta);
                p_year = p; p_year(idx_beta)=beta_y; p_year(idx_delta)=delta_y;
                aux.gfun = gfun; aux.t0 = 2010;
                ode_year = @(t,x) odefun(t,x,p_year,aux);
                if t_now < y_start
                    sol1 = ode15s(ode_year, [t_now y_start], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
                    x_now = deval(sol1, y_start); t_now = y_start;
                end
                sol2 = ode15s(ode_year, [y_start y_end], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
                tt   = linspace(y_start, y_end, nQuad);
                Xs   = deval(sol2, tt);
                Ssub = Xs(4,:); inc_rt = eta_y.*Ssub;
                yInc(ii) = trapz(tt, inc_rt) * 1e5;
                x_now = deval(sol2, y_end); t_now = y_end;
            else
                beta_y  = beta_piecewise( p(idx_beta),  yCal);
                delta_y = delta_piecewise(p(idx_delta), yCal);
                kappa_y = p(idx_kappa);
                p_year = p; p_year(idx_beta)=beta_y; p_year(idx_delta)=delta_y; p_year(idx_kappa)=kappa_y;
                aux.gfun = gfun; aux.t0 = 2010;
                ode_year = @(t,x) odefun(t,x,p_year,aux);
                if t_now < y_start
                    sol1 = ode15s(ode_year, [t_now y_start], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
                    x_now = deval(sol1, y_start); t_now = y_start;
                end
                sol2 = ode15s(ode_year, [y_start y_end], x_now, odeset('RelTol',1e-6,'AbsTol',1e-8));
                tt   = linspace(y_start, y_end, nQuad);
                Xs   = deval(sol2, tt);
                kLoc = k; Sk = Xs(3+kLoc,:); inc_rt = kappa_y.*Sk;
                yInc(ii) = trapz(tt, inc_rt) * 1e5;
                x_now = deval(sol2, y_end); t_now = y_end;
            end
        end
        yInc_row = reshape(yInc, 1, []);
    end

% ---- log-posterior (v2-like) ----
    function lp = logPosterior_u(u)
        p = u2p(u, LB, UB);
        [x2010, ~] = spinup_to_2010(p, k, inc1990_per100k);      % get 2010 state
        yM = simulate_incidence_annual_from2010(tData, p, x2010);% 1×T
        sdEff = sqrt(sdObs.^2 + sigma_add_like^2);
        r     = (yObs - yM.')./ sdEff;                           % column minus column
        ll    = -0.5*sum(r.^2 + log(2*pi*sdEff.^2));
        lp    = ll + logJac(u, LB, UB);
    end

% ---- MAP warm start ----
uCurr = zeros(1,nParam);
logPostCurr = logPosterior_u(uCurr);
mapLogPost  = logPostCurr;      % keep for model comparison
try
    neglogpost = @(uu) -logPosterior_u(uu);
    opts = optimoptions('fmincon','Algorithm','interior-point', ...
        'Display','off','MaxFunctionEvaluations',20000);
    [uMAP, fval, exitflag] = fmincon(neglogpost, zeros(1,nParam), [],[],[],[],[],[],[],opts);
    if isfinite(fval) && exitflag>0
        uCurr = uMAP; logPostCurr = -fval; mapLogPost = -fval;
        fprintf('MAP warm start used (exitflag=%d).\n', exitflag);
    else
        fprintf('MAP warm start failed (exitflag=%d); using zero start.\n', exitflag);
    end
catch ME
    fprintf('MAP warm start not available (%s); using zero start.\n', ME.identifier);
end

% ---- Adaptive Metropolis ----
chainU = zeros(nIter, nParam); Sigma = eye(nParam);
accTot = 0; accWin = 0;
for i=1:nIter
    uProp = uCurr + mvnrnd(zeros(1,nParam), s2*Sigma);
    logProp = logPosterior_u(uProp);
    if rand < exp(logProp - logPostCurr)
        uCurr = uProp; logPostCurr = logProp;
        accTot = accTot+1; accWin = accWin+1;
    end
    chainU(i,:) = uCurr;
    if i<=burnIn && mod(i,adaptIntv)==0
        Sigma = cov(chainU(1:i,:)) + epsilon*eye(nParam);
        fprintf('Iter %6d accWin=%4d (%.2f%%) max|u|=%.2f\n', ...
            i, accWin, 100*accWin/adaptIntv, max(abs(uCurr)));
        accWin = 0;
    end
end
fprintf('Overall acceptance: %.1f%%%%\n', 100*accTot/nIter);

% ---- robust thinning ----
thin = max(1, round(thin));
if nIter <= burnIn
    warning('Burn-in >= nIter; adjusting burnIn to nIter-1.'); 
    burnIn = max(0, nIter-1);
end
idx = (burnIn+1):thin:nIter;
if isempty(idx)
    warning('Thinning produced no samples; falling back to the final draw.');
    idx = nIter;
end
chainP    = u2p(chainU, LB, UB);
chainThin = chainP(idx, :);

% ---- posterior sims on history (for calibration plot) ----
nPost    = size(chainThin,1);
Yhist    = nan(nPost, numel(tData));          % 2010–2023
YpreHist = nan(nPost, numel(yearsPre));       % 1990–2009
for j=1:nPost
    p         = chainThin(j,:);
    [x2010, yPreRow] = spinup_to_2010(p, k, inc1990_per100k);
    yrow      = simulate_incidence_annual_from2010(tData, p, x2010);   % 1×T
    Yhist(j,:)= yrow;
    YpreHist(j,:) = yPreRow;
end
end  % run_AM_hist_likeV2



%% ====================== plotting (v2 style, with prehistory) ======================
function make_calibration_plot_with_prehistory(yearsPre, Ypre, yearsObs, yObsMean, yObsLow, yObsHigh, Yhist, ttl)
    greenLine = [0.38 0.55 0.18]; greenBand = [0.78 0.86 0.62];
    blueLine  = [0.12 0.46 0.70]; blueBand  = [0.68 0.85 0.95];
    grayLine  = [0.40 0.40 0.40]; grayBand  = [0.85 0.85 0.85];

    xPre   = yearsPre(:)';                       % 1990..2009
    xHist  = yearsObs(:)';                       % 2010..2023
    yMeanH = yObsMean(:)'; yLowH = yObsLow(:)'; yHighH = yObsHigh(:)';

    yMeanModel_hist = mean(Yhist,1);
    yLoModel_hist   = quantile(Yhist,0.025,1);
    yHiModel_hist   = quantile(Yhist,0.975,1);

    yMeanModel_pre = mean(Ypre,1);
    yLoModel_pre   = quantile(Ypre,0.025,1);
    yHiModel_pre   = quantile(Ypre,0.975,1);

    figure('Position',[140 140 980 560]); hold on; grid on; box on;

    % Observed band/line (2010–2023)
    fill([xHist, fliplr(xHist)], [yLowH, fliplr(yHighH)], greenBand, ...
        'EdgeColor','none','FaceAlpha',0.45, 'DisplayName','Observed 95% UI');
    plot(xHist, yMeanH, '-', 'Color', greenLine, 'LineWidth', 4, ...
        'DisplayName','Observed mean');

    % Model prehistory (1990–2009)
    fill([xPre, fliplr(xPre)], [yLoModel_pre, fliplr(yHiModel_pre)], grayBand, ...
        'EdgeColor','none','FaceAlpha',0.35, 'DisplayName','Model 95% (1990–2009)');
    plot(xPre, yMeanModel_pre, '-', 'Color', grayLine, 'LineWidth', 3, ...
        'DisplayName','Model mean (1990–2009)');

    % Model 2010–2023
    fill([xHist, fliplr(xHist)], [yLoModel_hist, fliplr(yHiModel_hist)], blueBand, ...
        'EdgeColor','none','FaceAlpha',0.35, 'DisplayName','Model 95% (2010–2023)');
    plot(xHist, yMeanModel_hist, '-', 'Color', blueLine, 'LineWidth', 4, ...
        'DisplayName','Model mean (2010–2023)');

    set(gca,'FontSize',12,'LineWidth',1);
    xlabel('Year'); ylabel('Rate per 100 000 population');
    title(ttl, 'FontSize',16, 'FontWeight','bold');
    yMax = max([yHighH(:); yHiModel_hist(:); yHiModel_pre(:); yMeanH(:); yMeanModel_hist(:); yMeanModel_pre(:)])*1.08;
    ylim([0, yMax]);
    xlim([xPre(1)-0.5, xHist(end)+0.5]);
    legend('Location','northeast');
end



%% ======================= ODE models with demography (Option A) =======================
function dx = simpleTBmodel_demog(t, x, p, aux)
% p = [btrans rhoS rhoM alpha beta gamma delta epsilon zeta eta theta muC muM muS nu
btr = p(1); rhoS = p(2); rhoM = p(3);
alpha = p(4); beta = p(5); gamma = p(6); delta = p(7); epsilon = p(8);
zeta = p(9); eta = p(10); theta = p(11); muC = p(12); muM = p(13); muS = p(14); nu = p(15);

S0 = x(1); I = x(2); M = x(3); S = x(4); C = x(5); R = x(6); D = x(7); %#ok<NASGU>

Nlive = max(S0 + I + M + S + C + R, 1e-12);
g = aux.gfun(aux.t0 + t);
Pi = (g + nu) * Nlive;          % births (Option A)
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
for j=2:k-1
    dS(j) = kappa*S(j-1) - (kappa + zeta + nu + muS)*S(j);
end

dS(k) = kappa*S(k-1) - (kappa + zeta + nu + muS)*S(k) ;
dC = kappa*S(k) - (theta + nu + muC)*C;

dR = alpha*I + delta*M - nu*R;
dD = nu*(S0 + I + M + sum(S) + C + R) + muC*C + muM*M + muS*sum(S);

dx = [dS0; dI; dM; dS; dC; dR; dD];
end


%% =================== [AUTO helpers for Scheme A] ===================
function [pMAP, ok] = AUTO_get_MAP_for_simple(LB, UB, tData, yObs, sdObs, sigma_add_like, yearsPre, gfun)
    ok = true;
    logistic = @(u) 1./(1+exp(-u));
    u2p_loc = @(u) LB + (UB-LB).*logistic(u);
    logJac_loc = @(u) sum(log(UB-LB) + log(logistic(u)) + log(1-logistic(u)));

    neglogpost = @(uu) -AUTO_logPosterior_simple_u(uu, LB, UB, tData, yObs, sdObs, ...
        sigma_add_like, yearsPre, gfun, u2p_loc, logJac_loc);

    u0 = zeros(1, numel(LB));
    try
        opts = optimoptions('fmincon','Algorithm','interior-point', ...
                'Display','off','MaxFunctionEvaluations',20000);
        [uMAP, fval, exitflag] = fmincon(neglogpost, u0, [],[],[],[],[],[],[],opts); %#ok<ASGLU>
        if exitflag>0 && isfinite(fval)
            pMAP = u2p_loc(uMAP);
        else
            ok = false; pMAP = (LB+UB)/2;
        end
    catch
        ok = false; pMAP = (LB+UB)/2;
    end
end

function lp = AUTO_logPosterior_simple_u(u, LB, UB, tData, yObs, sdObs, sigma_add_like, yearsPre, gfun, u2p, logJac)
    nQuad = 24;
    p = u2p(u);
    k = 1;
    [x2010, ~] = AUTO_spinup_to_2010(p, k, 1200, yearsPre, gfun, @simpleTBmodel_demog, nQuad); 
    yM = AUTO_simulate_from2010(tData, p, x2010, k, gfun, @simpleTBmodel_demog, nQuad);
    sdEff = sqrt(sdObs.^2 + sigma_add_like^2);
    r = (yObs - yM.')./sdEff;
    ll = -0.5*sum(r.^2 + log(2*pi*sdEff.^2));
    lp = ll + logJac(u);
end

function y2010 = AUTO_simulate_2010_inc_from_inc1990(p, k, inc1990, yearsPre, gfun, odefun)
    nQuad = 24;
    [x2010, ~] = AUTO_spinup_to_2010(p, k, inc1990, yearsPre, gfun, odefun, nQuad);
    y = AUTO_simulate_from2010(0, p, x2010, k, gfun, odefun, nQuad); % times_year_idx=0
    y2010 = y(1);
end

function [x2010, yPre_row] = AUTO_spinup_to_2010(p, k_local, inc1990, yearsPre, gfun, odefun, nQuad)
    b_beta = -0.00; b_delta = +0.00;
    if k_local==1
        idx_beta=5; idx_delta=7; idx_eta=10;
    else
        idx_beta=5; idx_delta=7; idx_kappa=9;
    end
    beta_piece = @(b0,yr)  b0.*exp(b_beta * max(yr-2010,0));
    delt_piece = @(d0,yr)  d0.*exp(b_delta* max(yr-2010,0));

    theta0=p(11); muC=p(12); C0=(inc1990/max(theta0+muC,1e-9))/1e5; I0=1e-3;
    if k_local==1
        eta0=p(idx_eta); delta0=p(7); zeta0=p(9); eps0=p(8);
        Ssub=max( ((theta0+muC)/max(eta0,1e-9))*C0 ,0);
        M0 =max( zeta0*Ssub / max(delta0+eps0,1e-9), 0);
        R0=0; D0=0; S0pool=max(1-(I0+M0+Ssub+C0+R0),1e-6);
        x_now=[S0pool; I0; M0; Ssub; C0; R0; D0];
    else
        kappa0=p(9); zeta=p(10); delta0=p(7);
        Svec=zeros(k_local,1);
        Sk=((theta0+muC)/max(kappa0,1e-9))*C0; Svec(k_local)=Sk;
        for jj=k_local-1:-1:1
            Svec(jj)=(kappa0/max(kappa0+zeta,1e-9))^(k_local-jj)*Sk;
        end
        M0=max( zeta*sum(Svec) / max(delta0+p(8),1e-9), 0);
        R0=0; D0=0; I0=1e-3; S0pool=max(1-(I0+M0+sum(Svec)+C0+R0),1e-6);
        x_now=[S0pool; I0; M0; Svec; C0; R0; D0];
    end

    t_now=0; yPre_row=nan(1,numel(yearsPre));
    for ii=1:numel(yearsPre)
        yCal=yearsPre(ii); y_start=(yCal-1990); y_end=y_start+1;
        if k_local==1
            p_year=p; p_year(idx_beta)=beta_piece(p(idx_beta),yCal);
            p_year(idx_delta)=delt_piece(p(idx_delta),yCal);
            aux.gfun=gfun; aux.t0=1990;
            f=@(t,x) odefun(t,x,p_year,aux);
            if t_now<y_start
                sol1=ode15s(f,[t_now y_start],x_now,odeset('RelTol',1e-6,'AbsTol',1e-8));
                x_now=deval(sol1,y_start); t_now=y_start;
            end
            sol2=ode15s(f,[y_start y_end],x_now,odeset('RelTol',1e-6,'AbsTol',1e-8));
            tt=linspace(y_start,y_end,nQuad); X=deval(sol2,tt);
            inc_rt=p(idx_eta)*X(4,:); yPre_row(ii)=trapz(tt,inc_rt)*1e5;
            x_now=deval(sol2,y_end); t_now=y_end;
        else
            p_year=p; p_year(idx_beta)=beta_piece(p(idx_beta),yCal);
            p_year(idx_delta)=delt_piece(p(idx_delta),yCal);
            p_year(idx_kappa)=p(idx_kappa);
            aux.gfun=gfun; aux.t0=1990;
            f=@(t,x) odefun(t,x,p_year,aux);
            if t_now<y_start
                sol1=ode15s(f,[t_now y_start],x_now,odeset('RelTol',1e-6,'AbsTol',1e-8));
                x_now=deval(sol1,y_start); t_now=y_start;
            end
            sol2=ode15s(f,[y_start y_end],x_now,odeset('RelTol',1e-6,'AbsTol',1e-8));
            tt=linspace(y_start,y_end,nQuad); X=deval(sol2,tt);
            kLoc=k_local; yPre_row(ii)=trapz(tt, p(idx_kappa)*X(3+kLoc,:))*1e5;
            x_now=deval(sol2,y_end); t_now=y_end;
        end
    end
    x2010 = x_now;
end

function yInc_row = AUTO_simulate_from2010(times_year_idx, p, x0_2010, k, gfun, odefun, nQuad)
    b_beta = -0.00; b_delta = +0.00;
    if k==1
        idx_beta=5; idx_delta=7; idx_eta=10;
    else
        idx_beta=5; idx_delta=7; idx_kappa=9;
    end
    beta_piece = @(b0,yr)  b0.*exp(b_beta * max(yr-2010,0));
    delt_piece = @(d0,yr)  d0.*exp(b_delta* max(yr-2010,0));

    if isempty(times_year_idx), yInc_row=[]; return; end
    ty=times_year_idx(:)'; x_now=x0_2010; t_now=0; yInc=nan(1,numel(ty));
    for ii=1:numel(ty)
        y_start=ty(ii); y_end=y_start+1; yCal=2010+y_start;
        if k==1
            p_year=p; p_year(idx_beta)=beta_piece(p(idx_beta),yCal);
            p_year(idx_delta)=delt_piece(p(idx_delta),yCal);
            aux.gfun=gfun; aux.t0=2010;
            f=@(t,x) odefun(t,x,p_year,aux);
            if t_now<y_start
                sol1=ode15s(f,[t_now y_start],x_now,odeset('RelTol',1e-6,'AbsTol',1e-8));
                x_now=deval(sol1,y_start); t_now=y_start;
            end
            sol2=ode15s(f,[y_start y_end],x_now,odeset('RelTol',1e-6,'AbsTol',1e-8));
            tt=linspace(y_start,y_end,nQuad); X=deval(sol2,tt);
            yInc(ii)=trapz(tt, p(idx_eta)*X(4,:))*1e5;
            x_now=deval(sol2,y_end); t_now=y_end;
        else
            p_year=p; p_year(idx_beta)=beta_piece(p(idx_beta),yCal);
            p_year(idx_delta)=delt_piece(p(idx_delta),yCal);
            p_year(idx_kappa)=p(idx_kappa);
            aux.gfun=gfun; aux.t0=2010;
            f=@(t,x) odefun(t,x,p_year,aux);
            if t_now<y_start
                sol1=ode15s(f,[t_now y_start],x_now,odeset('RelTol',1e-6,'AbsTol',1e-8));
                x_now=deval(sol1,y_start); t_now=y_start;
            end
            sol2=ode15s(f,[y_start y_end],x_now,odeset('RelTol',1e-6,'AbsTol',1e-8));
            tt=linspace(y_start,y_end,nQuad); X=deval(sol2,tt);
            kLoc=k; yInc(ii)=trapz(tt, p(idx_kappa)*X(3+kLoc,:))*1e5;
            x_now=deval(sol2,y_end); t_now=y_end;
        end
    end
    yInc_row = reshape(yInc,1,[]);
end

% ---- function ----
function [pMAP, ok] = AUTO_get_MAP_for_erlang(LB, UB, k, tData, yObs, sdObs, sigma_add_like, yearsPre, gfun)
    ok = true;
    logistic = @(u) 1./(1+exp(-u));
    u2p_loc = @(u) LB + (UB-LB).*logistic(u);
    logJac_loc = @(u) sum(log(UB-LB) + log(logistic(u)) + log(1-logistic(u)));

    neglogpost = @(uu) -AUTO_logPosterior_erlang_u(uu, LB, UB, k, tData, yObs, sdObs, ...
        sigma_add_like, yearsPre, gfun, u2p_loc, logJac_loc);

    u0 = zeros(1, numel(LB));
    try
        opts = optimoptions('fmincon','Algorithm','interior-point', ...
                'Display','off','MaxFunctionEvaluations',20000);
        [uMAP, fval, exitflag] = fmincon(neglogpost, u0, [],[],[],[],[],[],[],opts); %#ok<ASGLU>
        if exitflag>0 && isfinite(fval)
            pMAP = u2p_loc(uMAP);
        else
            ok = false; pMAP = (LB+UB)/2;
        end
    catch
        ok = false; pMAP = (LB+UB)/2;
    end
end

function lp = AUTO_logPosterior_erlang_u(u, LB, UB, k, tData, yObs, sdObs, sigma_add_like, yearsPre, gfun, u2p, logJac)
    nQuad = 24;                     
    p = u2p(u);
    odefun = @(t,x,pp,aux) erlangODE_demog(t,x,pp,aux,k);

    
    inc1990_neutral = 1200;
    [x2010, ~] = AUTO_spinup_to_2010(p, k, inc1990_neutral, yearsPre, gfun, odefun, nQuad);

    yM = AUTO_simulate_from2010(tData, p, x2010, k, gfun, odefun, nQuad);
    sdEff = sqrt(sdObs.^2 + sigma_add_like^2);
    r = (yObs - yM.')./sdEff;
    ll = -0.5*sum(r.^2 + log(2*pi*sdEff.^2));
    lp = ll + logJac(u);
end
