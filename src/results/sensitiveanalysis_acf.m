% ====================== forecast_ACF_10yr_SA_printOnly.m =======================
% goal：
% - 
% - keep print：Simple ACF、Erlang ACF (2024–26) 的 2024–2033 prediction (Mean + 95% UI)
% - acf.mult_delta_up / acf.mult_zeta_up ：
%       multipliers = [1.2, 1.5, 1.7, 2.0]
%
% - 1990 spin-up -> 2010-2023 -> x2024 -> 2024-2033 ACF
% - baseline，no history UI
%
clear; clc; close all; rng(20250609);

%% ---- Load calibrated outputs ----
S = load('postSamples_simple_clinical_AM.mat');        % Simple posterior: chainThin, inc1990_used, tData_save, yObs_save
E = load('postSamples_erlang_clinical_AM_50.mat');     % Erlang posterior: chainThin, inc1990_used, tData_save, yObs_save

yearsObs = S.tData_save(:)';      %#ok<NASGU>  % 2010..2023 
meanInc  = S.yObs_save(:)';       %#ok<NASGU>  % WHO mean 
lowInc  = [114 138 141 141 171 181 191 193 185 180 172 173 170 164]; %#ok<NASGU>
uppInc  = [508 440 410 395 329 300 263 242 233 227 220 228 231 228]; %#ok<NASGU>

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

%% ---- Draws ----
n_draws_simple = min(400, size(S.chainThin,1));
n_draws_erlang = min(400, size(E.chainThin,1));
idxS = randperm(size(S.chainThin,1), n_draws_simple);
idxE = randperm(size(E.chainThin,1), n_draws_erlang);

%% ---- ACF settings (sensitivity multipliers) ----
acf_base.yStart = 2024;
acf_base.yEnd   = 2026;

mult_list = [1.2, 1.5, 1.7, 2.0];   % sensitivity analysis grid

%% ---- Time-slope gates ----
b_beta  = 0.00;
b_delta = 0.00;

%% ---- Quadrature ----
nQuad = 24;

%% ---- Erlang structure ----
k_erlang = 50; %#ok<NASGU>  % simulate_annual_erlang 内部固定 k=50


fprintf('Precomputing x2024 for Simple (%d draws) and Erlang (%d draws)...\n', ...
    n_draws_simple, n_draws_erlang);

x2024_simple = cell(n_draws_simple, 1);
for j=1:n_draws_simple
    p = S.chainThin(idxS(j),:);
    [x2010, ~] = spin_from_1990_simple(p, inc1990, 2010, [], gfun, nQuad, b_beta, b_delta);
    [~, x2024] = simulate_annual_simple(p, 2010, 2023, x2010, gfun, nQuad, b_beta, b_delta, []);
    x2024_simple{j} = x2024;
end

x2024_erlang = cell(n_draws_erlang, 1);
for j=1:n_draws_erlang
    p = E.chainThin(idxE(j),:);
    [x2010, ~] = spin_from_1990_erlang(p, inc1990, 2010, [], gfun, nQuad, b_beta, b_delta);
    [~, x2024] = simulate_annual_erlang(p, 2010, 2023, x2010, gfun, nQuad, b_beta, b_delta, []);
    x2024_erlang{j} = x2024;
end

fprintf('Done.\n');


yearsF = 2024:2033;

for mm = 1:numel(mult_list)
    mult = mult_list(mm);

    acf = acf_base;
    acf.mult_delta_up = mult;   % M->R
    acf.mult_zeta_up  = mult;   % S->M

    fprintf('\n============================================================\n');
    fprintf('ACF sensitivity: mult_delta_up = %.2f, mult_zeta_up = %.2f (years %d-%d)\n', ...
        acf.mult_delta_up, acf.mult_zeta_up, acf.yStart, acf.yEnd);
    fprintf('Forecast years: %d-%d\n', yearsF(1), yearsF(end));
    fprintf('============================================================\n');

    % ---- Simple ACF ----
    YfS_acf = zeros(n_draws_simple, numel(yearsF));
    for j=1:n_draws_simple
        p = S.chainThin(idxS(j),:);
        x0 = x2024_simple{j};
        [yf, ~] = simulate_annual_simple(p, yearsF(1), yearsF(end), x0, gfun, nQuad, b_beta, b_delta, acf);
        YfS_acf(j,:) = yf;
    end

    fcastS_acf_mean = mean(YfS_acf,1);
    fcastS_acf_lo   = prctile(YfS_acf,  2.5, 1);
    fcastS_acf_hi   = prctile(YfS_acf, 97.5, 1);

    print_UI_table(sprintf('Simple ACF (mult=%.2f, 2024–26)', mult), yearsF, ...
        fcastS_acf_mean, fcastS_acf_lo, fcastS_acf_hi);

    % ---- Erlang ACF ----
    YfEr_acf = zeros(n_draws_erlang, numel(yearsF));
    for j=1:n_draws_erlang
        p = E.chainThin(idxE(j),:);
        x0 = x2024_erlang{j};
        [yf, ~] = simulate_annual_erlang(p, yearsF(1), yearsF(end), x0, gfun, nQuad, b_beta, b_delta, acf);
        YfEr_acf(j,:) = yf;
    end

    fcastEr_acf_mean = mean(YfEr_acf,1);
    fcastEr_acf_lo   = prctile(YfEr_acf,  2.5, 1);
    fcastEr_acf_hi   = prctile(YfEr_acf, 97.5, 1);

    print_UI_table(sprintf('Erlang ACF (mult=%.2f, 2024–26)', mult), yearsF, ...
        fcastEr_acf_mean, fcastEr_acf_lo, fcastEr_acf_hi);

end

fprintf('\nAll sensitivity scenarios finished.\n');

%% ================== Helpers ==================
function print_UI_table(label, years, meanVec, loVec, hiVec)
    fprintf('\n==== %s ====\n', label);
    fprintf('%6s %12s %12s %12s\n','Year','Mean','Lo95','Hi95');
    for i = 1:numel(years)
        fprintf('%6d %12.3f %12.3f %12.3f\n', ...
            years(i), meanVec(i), loVec(i), hiVec(i));
    end
end

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
