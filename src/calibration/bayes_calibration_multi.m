clear; clc; close all;
rng(20250609);

%% A) Load data ----------------------------------------------------------
dataFiles = {
  'data/data_nti1_preweight.csv'
  'data/data_nti2_preweight.csv'
  'data/data_nti3_preweight.csv'
};
nDatasets = numel(dataFiles);
tData = cell(nDatasets,1);
yObs = cell(nDatasets,1);
ciLower = cell(nDatasets,1);
ciUpper = cell(nDatasets,1);
x0_list = cell(nDatasets,1);

k = 50;  % Erlang subclinical stages (number of compartments)
         % Controls the shape of the subclinical duration distribution:
         % k = 1 → exponential (memoryless)
         % larger k → more concentrated (less variance, more symmetric)
         % k can be adjusted depending on the desired level of biological realism
         % while keeping the mean duration fixed via kappa = k * eta
for d = 1:nDatasets
    T        = readtable(dataFiles{d});
    atRisk   = T.atrisk;   
    casesObs = T.cases;   
    tData{d} = T.year;     
    yObs{d}  = casesObs ./ atRisk;  

    
    lambda   = casesObs; 
    N0       = atRisk;   
    lo_cases = 0.5*chi2inv(0.025, 2*lambda);        % Poisson
    hi_cases = 0.5*chi2inv(0.975, 2*(lambda+1));    % Poisson
    rate_lo  = lo_cases ./ N0;
    rate_hi  = hi_cases ./ N0;
    ciLower{d} = yObs{d} - rate_lo;
    ciUpper{d} = rate_hi - yObs{d};

 
    x0_list{d} = [atRisk(1); zeros(k+4,1)];
end

%% B) Prior bounds (NOW AS ROW-VECTORS) ---------------------------------
nParam  = 9;
priorLB = [0.93 0.04 0.01 0.14 0.08 1.5*k 1.24 0.46 0.28];  % 1×9
priorUB = [3.3 0.23 0.1   0.23 1.50  4*k  2.03 0.72 0.38]; % 1×9
% p = [alpha beta gamma delta epsilon kappa zeta theta mu]
% zeta: Avoid being forced to reflux early on(different from simple)
logistic  = @(u) 1 ./ (1+exp(-u));                  % element-wise
p2u = @(p) log( (p-priorLB) ./ (priorUB-priorLB) ); % works for rows or matrices
u2p = @(u) priorLB + (priorUB-priorLB) .* logistic(u); % ditto

%% C) MCMC settings ------------------------------------------------------
nIter  = 2e5;  burnIn = 5e4;  thin = 10;
s2     = 0.8*(2.38^2) / nParam;  % AM scaling
epsilon= 1e-5;                   % nugget
adaptIntv = 1000;

%% D) Initialise ---------------------------------------------------------
uCurr = zeros(1,nParam);  % row vector in u-space
pCurr = u2p(uCurr);

LPsum = 0;
for d = 1:nDatasets
    LPsum = LPsum + logPosterior(pCurr, tData{d}, yObs{d}, x0_list{d}, k);
end
logPostCurr = LPsum/nDatasets + logJacobian(uCurr, priorLB, priorUB);

chainU = zeros(nIter, nParam);
Sigma  = eye(nParam);
accTot = 0;  accWin = 0;

%% E) Adaptive-Metropolis ------------------------------------------------
for i = 1:nIter
    uProp = uCurr + mvnrnd(zeros(1,nParam), s2*Sigma);
    pProp = u2p(uProp);

    LPsum = 0;
    for d = 1:nDatasets
        LPsum = LPsum + logPosterior(pProp, tData{d}, yObs{d}, x0_list{d}, k);
    end
    logProp = LPsum/nDatasets + logJacobian(uProp, priorLB, priorUB);

    if rand < exp(logProp - logPostCurr)
        uCurr = uProp;  pCurr = pProp;  logPostCurr = logProp;
        accTot = accTot + 1;  accWin = accWin + 1;
    end
    chainU(i,:) = uCurr;

    if i <= burnIn && mod(i,adaptIntv) == 0
        Sigma = cov(chainU(1:i,:)) + epsilon*eye(nParam);
        fprintf('Iter %6d  accWin=%4d (%.2f)  max|u|=%.2f\n', ...
                i, accWin, accWin/adaptIntv, max(abs(uCurr)));
        accWin = 0;
    end
end
fprintf('\nOverall acceptance: %.1f%%\n', 100*accTot/nIter);

%% F) Back-transform & thin ---------------------------------------------
chainP    = u2p(chainU);                     % now dimension-safe
chainThin = chainP(burnIn+1:thin:end, :);

postMean = mean(chainThin,1);
postCI   = quantile(chainThin, [0.025 0.975]);

paramNames = {'alpha','beta','gamma','delta','epsilon', ...
              'kappa','zeta','theta','mu'};
Tpost = table(paramNames(:), postMean(:), postCI(1,:)', postCI(2,:)', ...
              'VariableNames', {'Param','Mean','CI_2.5%','CI_97.5%'} );
disp(Tpost);

%% G) Forecast & plotting -----------------------------------------------
tPred = (1:10)';                
nPost = size(chainThin,1);
Ypred = nan(nPost, numel(tPred));
x0_f  = x0_list{1};

for j = 1:nPost
    Ypred(j,:) = annual_inflow_rate(chainThin(j,:), k, x0_f, tPred(end));
end
yMean = mean(Ypred,1);
yLow  = quantile(Ypred, 0.025, 1);
yHigh = quantile(Ypred, 0.975, 1);

figure('Position',[200 200 760 470]); hold on; grid on;
markers = {'^','o','s'};  colors = {'b','m',[1 .5 0]};  labels = {'NTI1','NTI2','NTI3'};
for d = 1:nDatasets
    errorbar(tData{d}, yObs{d}, ciLower{d}, ciUpper{d}, markers{d}, ...
             'Color', colors{d}, 'MarkerFaceColor', colors{d}, 'LineWidth', 1.3);
end
h = fill([tPred; flipud(tPred)], [yLow'; flipud(yHigh')], 'c', 'EdgeColor', 'none');
set(h, 'FaceAlpha', 0.2);
plot(tPred, yMean, '-k', 'LineWidth', 2);
xlabel('Years since infection');
ylabel('Annual net inflow to subclinical (rate)');
legend([labels, {'95% posterior band','Posterior mean'}], 'Location', 'northeast');
title('Fit to annual net inflow into subclinical (Erlang model)');

%% H)  Save --------------------------------------------------------------
save('postSamples_subclinical_uniform_adaptcov.mat', ...
     'chainThin','tData','k','yObs','x0_list','tPred','Ypred');

%% ---- helpers ----------------------------------------------------------
function lj = logJacobian(u, lb, ub)
    sig = 1./(1+exp(-u));
    lj  = sum( log(ub-lb) + log(sig) + log(1-sig) );
end

% --------- LOG-POSTERIOR (incidence target) -----------------------------
function lp = logPosterior(p, tGridYears, yObs, x0, k)
% p: [alpha beta gamma delta epsilon kappa zeta theta mu]
    N0      = x0(1);
    r_model = annual_inflow_rate(p, k, x0, max(tGridYears));
    r_model = r_model(tGridYears);                            

    cases_hat = yObs(:) * N0;
    pi_hat    = max(min(r_model(:), 1-1e-9), 1e-9);           

    ll = cases_hat .* log(pi_hat) + (N0 - cases_hat) .* log(1 - pi_hat);

    lp = sum(ll);
end

% --------- Compute annual inflow rate vector ----------------------------
function r = annual_inflow_rate(p, k, x0, T_max)
% r(y) = (1/N0) * ∫_{y-1}^{y} [gamma*I(t) + epsilon*M(t)] dt
    N0 = x0(1);
    tspan   = [0 T_max];
    opts    = odeset('RelTol',1e-7,'AbsTol',1e-9);
    sol     = ode45(@(t,x) erlangODE(t,x,p,k), tspan, x0, opts);

    npt_per_year = 200;
    r = nan(T_max,1);
    for y = 1:T_max
        tt  = linspace(y-1, y, npt_per_year);
        XY  = deval(sol, tt);
        I   = XY(1,:);              % I(t)
        M   = XY(2,:);              % M(t)
        gamma  = p(3);  epsilon = p(5);
        Fin = gamma.*I + epsilon.*M;                   
        inflow = trapz(tt, Fin);                        
        r(y)   = inflow / N0;                          
    end
end

% --------- ODE system ---------------------------------------------------
function dx = erlangODE(~, x, p, k)
% States order: [I, M, S1..Sk, C, R, D]  (length = k+5)
% Params: p = [alpha beta gamma delta epsilon kappa zeta theta mu]
    alpha  = p(1);
    beta   = p(2);
    gamma  = p(3);
    delta  = p(4);
    epsl   = p(5);   % epsilon
    kappa  = p(6);
    zeta   = p(7);
    theta  = p(8);
    mu     = p(9);

    I = x(1);
    M = x(2);
    S = x(3:2+k);
    C = x(3+k);
    % R = x(4+k); 
    % D = x(5+k);

    dx = zeros(size(x));

    % dI/dt = -(alpha + beta + gamma) I
    dx(1) = -(alpha + beta + gamma) * I;

    % dM/dt = beta I - (delta + epsl) M + zeta * sum(Sj)
    dx(2) = beta*I - (delta + epsl)*M + zeta*sum(S);

    % dS1/dt + theta C
    dx(3) = gamma*I + epsl*M - (kappa + zeta)*S(1) + theta * x(k+3);

    % dSj/dt, j=2..k-1
    if k >= 3
        for j = 2:(k-1)
            dx(2+j) = kappa*S(j-1) - (kappa + zeta)*S(j);
        end
    end

    % dSk/dt = kappa S_{k-1} - (kappa + zeta) S_k
    dx(2+k) = kappa*S(k-1) - (kappa + zeta)*S(k);

    % dC/dt = kappa S_k - (theta + mu) C
    dx(3+k) = kappa*S(k) - (theta + mu)*x(3+k);

    % dR/dt = alpha I + delta M
    dx(4+k) = alpha*I + delta*M;

    % dD/dt = mu C
    dx(5+k) = mu*x(3+k);
end
