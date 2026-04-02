% evaluate_who_india_metrics.m
% -------------------------------------------------------------------------
% Evaluate DIC / WAIC / (optional) PSIS-LOOIC / Avg Log Score for the
% Simple & Erlang models fitted by calibrate_who_india.m
% IMPORTANT: Recompute per-year log-likelihood using ANNUAL INCIDENCE
% (eta*S for Simple; kappa*S_k for Erlang), with parameter-consistent x0.
% Do NOT use Ypred (C stock) to form likelihood.
% -------------------------------------------------------------------------
clear; clc; close all;

%% --- MUST MATCH your calibration knob -----------------------------------
sigma_add_like = 10;    % per 100k, same as in calibrate_who_india.m

%% --- WHO India clinical incidence (2010–2023) ---------------------------
years = (2010:2023)';
meanInc = [276 268 258 249 243 237 225 217 208 202 195 200 199 195]';
lowInc  = [114 138 141 141 171 181 191 193 185 180 172 173 170 164]';
uppInc  = [508 440 410 395 329 300 263 242 233 227 220 228 231 228]';
sdObs   = max((uppInc - lowInc)/3.92, 1e-3);   % normal approx to 95% CI
sdEff   = sqrt(sdObs.^2 + sigma_add_like^2);  % effective obs SD (per 100k)
nObs    = numel(years);
tData_idx = (0:(nObs-1))';   % 0..13 for 2010..2023
firstYear = years(1);        % 2010

%% --- Files produced by calibrate_who_india.m ----------------------------
files = { ...
  'postSamples_simple_clinical_AM.mat', ...
  'postSamples_erlang_clinical_AM.mat' ...
};

modelNames = {'Simple (k=1)','Erlang (k=2)'};

outRows = [];

for m = 1:numel(files)
    if ~isfile(files{m})
        warning('Cannot find %s. Skipping.', files{m});
        continue;
    end

    S = load(files{m}, 'chainThin','tData','k','yObs','x0_list','tPred','Ypred');
    chainThin = S.chainThin;   % posterior draws [S x nParam]
    k         = S.k;           % 1 or 2 (your calib used k=2 for Erlang)
    nPost     = size(chainThin,1);

    % Build ODE handle (match calib)
    if k==1
        odefun = @simpleTBmodel;
    else
        odefun = @(t,x,p) erlangODE(t,x,p,k);
    end

    % ===================== Recompute annual incidence =====================
    IncMat = nan(nPost, nObs);   % [S x nObs]
    for s = 1:nPost
        p = chainThin(s,:);
        IncMat(s,:) = simulate_incidence_annual(odefun, p, k, tData_idx, meanInc(1));
    end

    % ---- Log-likelihood matrix on incidence scale ------------------------
    Y  = repmat(meanInc', nPost, 1);   % [S x nObs]
    SD = repmat(sdEff',  nPost, 1);
    logLik = -0.5 * ( ((Y - IncMat)./SD).^2 + log(2*pi*(SD.^2)) );  % [S x nObs]

    % ---- WAIC (pointwise) ------------------------------------------------
    lppd_i   = logmeanexp(logLik, 1);      % 1 x nObs
    pWAIC_i  = var(logLik, 0, 1);          % 1 x nObs  (pWAIC2)
    lppd     = sum(lppd_i);
    pWAIC    = sum(pWAIC_i);
    WAIC     = -2 * (lppd - pWAIC);

    % ---- DIC (DIC2) ------------------------------------------------------
    mu_post_mean = mean(IncMat, 1)';                 % [nObs x 1]
    Dbar  = mean(-2 * sum(logLik, 2));
    logLik_hat_i = -0.5 * ( ((meanInc - mu_post_mean)./sdEff).^2 + log(2*pi*(sdEff.^2)) );
    Dhat = -2 * sum(logLik_hat_i);
    pD   = Dbar - Dhat;
    DIC  = 2*Dbar - Dhat;

    % ---- Avg Log Score (per observation) ---------------------------------
    avgLogScore = mean(lppd_i);

    % ---- Optional: PSIS-LOO ---------------------------------------------
    doPSIS = false;   
    if doPSIS
        try
            [elpd_loo, se_elpd_loo, looic, pareto_k] = psisloo(logLik);
        catch ME
            warning('PSIS-LOO failed for %s: %s\nTurning off doPSIS.', modelNames{m}, ME.message);
            doPSIS = false; looic = NaN; se_elpd_loo = NaN; pareto_k = NaN(size(lppd_i));
        end
    else
        looic = NaN; se_elpd_loo = NaN; pareto_k = NaN(size(lppd_i));
    end

    % ---- Pack results ----------------------------------------------------
    res.model         = modelNames{m};
    res.nPost         = nPost;
    res.DIC           = DIC;
    res.pD            = pD;
    res.WAIC          = WAIC;
    res.lppd          = lppd;
    res.pWAIC         = pWAIC;
    res.AvgLogScore   = avgLogScore;
    res.LOOIC         = looic;
    res.se_elpd_loo   = se_elpd_loo;
    res.max_pareto_k  = max(pareto_k);
    outRows = [outRows; res]; %#ok<AGROW>

    % Quick console view
    fprintf('\n== %s ==\n', modelNames{m});
    fprintf('nPost = %d\n', nPost);
    fprintf('DIC   = %.3f   (pD = %.3f)\n', DIC, pD);
    fprintf('WAIC  = %.3f   (pWAIC = %.3f)\n', WAIC, pWAIC);
    fprintf('Avg Log Score = %.4f (per obs)\n', avgLogScore);
    if doPSIS
        fprintf('LOOIC = %.3f   (se elpd_loo = %.3f)   max Pareto k = %.3f\n', ...
            looic, se_elpd_loo, max(pareto_k));
    end
end

%% --- Nicely formatted comparison table ----------------------------------
if ~isempty(outRows)
    T = struct2table(outRows);
    disp(' ');
    disp('==== Model comparison (smaller is better for DIC/WAIC/LOOIC) ====');
    disp(T(:, {'model','DIC','WAIC','LOOIC','AvgLogScore','pD','pWAIC','nPost','max_pareto_k'}));
end

%% ======================= helpers ========================================
function y = logmeanexp(A, dim)
    % log(mean(exp(A), dim)) computed stably
    if nargin < 2, dim = 1; end
    M = max(A, [], dim);
    y = M + log(mean(exp(A - M), dim));
end

% =========== Annual-incidence simulator (with parameter-consistent x0) ===
function yInc = simulate_incidence_annual(odefun, p, k, times_year_idx, inc0_per100k)
    if isempty(times_year_idx)
        yInc = [];
        return;
    end
    ty = times_year_idx(:)';
    x_now = make_x0_from_params(p, k, inc0_per100k);
    t_now = 0;
    yInc  = nan(1, numel(ty));
    for ii = 1:numel(ty)
        y_start = ty(ii);
        if t_now < y_start
            opts1 = odeset('RelTol',1e-8,'AbsTol',1e-10);
            sol1  = ode45(@(t,x) odefun(t,x,p), [t_now y_start], x_now, opts1);
            x_now = deval(sol1, y_start);
            t_now = y_start;
        end
        y_end  = y_start + 1;
        opts2  = odeset('RelTol',1e-8,'AbsTol',1e-10);
        sol2   = ode45(@(t,x) odefun(t,x,p), [y_start y_end], x_now, opts2);
        tt     = linspace(y_start, y_end, 200);
        Xs     = deval(sol2, tt);
        if k==1
            S      = Xs(3,:); 
            eta    = p(7);
            inc_rt = eta .* S;
        else
            Sk     = Xs(2+k,:);
            kappa  = p(6);
            inc_rt = kappa .* Sk;
        end
        yInc(ii) = trapz(tt, inc_rt) * 1e5;   % annual integral per 100k
        x_now = deval(sol2, y_end);
        t_now = y_end;
    end
    yInc = yInc(:)';   % row vector
end

% ---------------- make_x0_from_params (shared with calib logic) ----------
function x0_here = make_x0_from_params(p, k, inc0_per100k)
    C0 = (inc0_per100k / max(p(end-1)+p(end), 1e-6)) / 1e5;
    if k==1
        eta = p(7);
        S0  = ((p(end-1)+p(end))/max(eta,1e-9)) * C0;
        M0  = (p(6)*S0) / max(p(4)+p(5),1e-9);
        x0_here = [0.001; M0; S0; C0; 0; 0];
    else
        kappa = p(6); zeta = p(7);
        Sk = ((p(end-1)+p(end))/max(kappa,1e-9)) * C0;
        Svec = zeros(k,1);
        for j=1:k-1
            Svec(j) = (kappa/max(kappa+zeta,1e-9))^(k-j) * Sk;
        end
        Svec(k) = Sk;
        M0  = (zeta*sum(Svec)) / max(p(4)+p(5),1e-9);
        x0_here = zeros(k+5,1);
        x0_here(1)=0.001; x0_here(2)=M0; x0_here(3:2+k)=Svec; x0_here(3+k)=C0;
    end
end

%% ======================= ODE models (copied) ============================
function dx = simpleTBmodel(~, x, p)
    % p = [alpha beta gamma delta epsilon zeta eta theta mu]
    alpha   = p(1); beta = p(2); gamma = p(3); delta = p(4); epsilon = p(5);
    zeta    = p(6); eta   = p(7); theta = p(8); mu      = p(9);
    I = x(1); M = x(2); S = x(3); C = x(4); R = x(5); D = x(6); %#ok<NASGU>
    dI = -(alpha + beta + gamma) * I;
    dM =  beta * I - (delta + epsilon) * M + zeta * S;
    dS =  gamma * I + epsilon * M - (zeta + eta) * S;
    dC =  eta   * S - (theta + mu) * C;
    dR =  alpha * I + delta * M;
    dD =  mu    * C;
    dx = [dI; dM; dS; dC; dR; dD];
end

function dx = erlangODE(~, x, p, k)
    % p = [alpha beta gamma delta epsilon kappa zeta theta mu]
    alpha = p(1); beta = p(2); gamma = p(3); delta = p(4); epsilon = p(5);
    kappa = p(6); zeta = p(7); theta = p(8); mu = p(9);
    I = x(1); M = x(2); S = x(3:2+k); C = x(3+k); R = x(4+k); D = x(5+k); %#ok<NASGU>
    dI = -(alpha + beta + gamma)*I;
    dM = beta*I - (delta + epsilon)*M + zeta*sum(S);
    dS = zeros(k,1);
    dS(1) = gamma*I + epsilon*M - (kappa + zeta)*S(1);
    for j = 2:k
        if j < k
            dS(j) = kappa*S(j-1) - (kappa + zeta)*S(j);
        else
            dS(j) = kappa*S(j-1) - (kappa + zeta)*S(j) + theta*C;
        end
    end
    dC = kappa*S(k) - (theta + mu)*C;
    dR = alpha*I + delta*M;
    dD = mu*C;
    dx = [dI; dM; dS; dC; dR; dD];
end

% --------------- Optional PSIS-LOO (unchanged) ---------------------------
function [elpd_loo, se_elpd_loo, looic, k_hat] = psisloo(log_lik)
    [S, n] = size(log_lik);
    lw = bsxfun(@minus, log_lik, max(log_lik,[],1));
    r = -lw; w = exp(-r); w = bsxfun(@rdivide, w, mean(w,1));
    k_hat = zeros(1,n); elpd_i = zeros(1,n);
    tail_frac = 0.2; m = max(10, round(tail_frac * S));
    for i = 1:n
        wi = w(:,i);
        [wi_sorted, ~] = sort(wi, 'descend');
        tail = wi_sorted(1:m); u = tail(end); excess = tail - u;
        try
            [parmHat, ~] = gpfit(excess); k = parmHat(1);
        catch
            k = NaN;
        end
        k_hat(i) = k;
        if ~isnan(k) && k < 1
            q = ((1:m)' - 0.5) / m;
            sm_tail = u + gpinv(q, k, parmHat(2));
            wi_sm   = wi_sorted; wi_sm(1:m) = sm_tail;
        else
            wi_sm = wi_sorted;
        end
        wi2 = wi; wi2(1:m) = wi_sm(1:m);
        elpd_i(i) = -log(mean(wi2));
    end
    elpd_loo = sum(elpd_i);
    se_elpd_loo = sqrt(n * var(elpd_i));
    looic = -2 * elpd_loo;
end
