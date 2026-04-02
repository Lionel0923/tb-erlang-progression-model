clear; clc; close all;
rng(20250609);

%% A) Load data -------------------------------------------------------------
dataFiles = {
  'data/data_nti1_preweight.csv'
  'data/data_nti2_preweight.csv'
  'data/data_nti3_preweight.csv'
};
nDatasets = numel(dataFiles);
tData     = cell(nDatasets,1);
yObs      = cell(nDatasets,1);
ciLower   = cell(nDatasets,1);
ciUpper   = cell(nDatasets,1);
x0_list   = cell(nDatasets,1);

for d = 1:nDatasets
    T          = readtable(dataFiles{d});
    atRisk     = T.atrisk;
    casesObs   = T.cases;
    tData{d}   = T.year;
    yObs{d}    = casesObs ./ atRisk;  
    
    [~,pci]    = binofit(casesObs, atRisk, 0.05);
    ciLower{d} = yObs{d} - pci(:,1);
    ciUpper{d} = pci(:,2) - yObs{d};
    % [I; M; S; C; R; D]
    x0_list{d} = [atRisk(1); zeros(5,1)];
end

%% B) Uniform prior bounds (row vectors) ------------------------------------
nParam  = 9;
priorLB = [0.93 0.04 0.01 0.14 0.08 0.08 1.24 0.46 0.28];  % 1×9
priorUB = [3.3 0.23 0.1   0.23 1.50 2.93 2.03   0.72 0.38]; % 1×9
% p = [alpha beta gamma delta epsilon zeta eta theta mu]

% logit-transform functions
logistic = @(u) 1./(1+exp(-u));
p2u      = @(p) log( (p - priorLB) ./ (priorUB - p) );  
u2p      = @(u) priorLB + (priorUB - priorLB) .* logistic(u);

%% C) MCMC settings ---------------------------------------------------------
nIter     = 2e5;
burnIn    = 5e4;
thin      = 10;
s2        = (2.38^2)/nParam;
epsilon   = 1e-5;
adaptIntv = 1000;

%% D) Initialize ------------------------------------------------------------
uCurr = zeros(1,nParam);         % start at midpoint -> p = (LB+UB)/2
pCurr = u2p(uCurr);

% compute initial log-posterior
LPsum = 0;
for d = 1:nDatasets
    LPsum = LPsum + logPosterior_simple(pCurr, tData{d}, yObs{d}, x0_list{d});
end
logPostCurr = LPsum/nDatasets + logJacobian(uCurr, priorLB, priorUB);

chainU = zeros(nIter, nParam);
Sigma  = eye(nParam);
accTot = 0;
accWin = 0;

%% E) Adaptive Metropolis over u-space -------------------------------------
for i = 1:nIter
    % propose in u-space
    uProp = uCurr + mvnrnd(zeros(1,nParam), s2*Sigma);
    pProp = u2p(uProp);

    % compute log-posterior for pProp
    LPsum = 0;
    for d = 1:nDatasets
        LPsum = LPsum + logPosterior_simple(pProp, tData{d}, yObs{d}, x0_list{d});
    end
    logProp = LPsum/nDatasets + logJacobian(uProp, priorLB, priorUB);

    % accept/reject
    if rand < exp(logProp - logPostCurr)
        uCurr       = uProp;
        pCurr       = pProp;
        logPostCurr = logProp;
        accTot      = accTot + 1;
        accWin      = accWin + 1;
    end
    chainU(i,:) = uCurr;

    % adapt covariance during burn-in
    if i <= burnIn && mod(i, adaptIntv) == 0
        Sigma = cov(chainU(1:i,:)) + epsilon*eye(nParam);
        fprintf('Iter %6d  accWin=%4d (%.2f)\n', i, accWin, accWin/adaptIntv);
        accWin = 0;
    end
end
fprintf('\nOverall acceptance: %.1f%%\n', 100*accTot/nIter);

%% F) Back-transform, thin, and summarize ---------------------------------
chainP    = u2p(chainU);                     
chainThin = chainP(burnIn+1:thin:end, :);

postMean = mean(chainThin,1);
postCI   = quantile(chainThin, [0.025 0.975]);

paramNames = {'alpha','beta','gamma','delta','epsilon',...
              'zeta','eta','theta','mu'};
Tpost = table(paramNames(:), postMean(:), postCI(1,:)', postCI(2,:)', ...
    'VariableNames', {'Param','Mean','CI_2.5%','CI_97.5%'}); %#ok<NASGU>
disp(Tpost);

%% G) Forecast & plot (simple TB model) ------------------------------------
tPred = (1:10)';     
nPost = size(chainThin,1);
Ypred = nan(nPost, numel(tPred));
x0_f  = x0_list{1};

for j = 1:nPost
    % r(y) = (1/N0) ∫_{y-1}^{y} [gamma*I + epsilon*M] dt
    Ypred(j,:) = annual_inflow_rate_simple(chainThin(j,:), x0_f, tPred(end));
end

yMean = mean(Ypred,1);  
yLow  = quantile(Ypred,0.025,1);  
yHigh = quantile(Ypred,0.975,1);

figure('Position',[200 200 700 450]); hold on; grid on;
markers = {'^','o','s'}; colors = {'b','m',[1 .5 0]}; labels = {'NTI1','NTI2','NTI3'};
for d = 1:nDatasets
    errorbar(tData{d}, yObs{d}, ciLower{d}, ciUpper{d}, markers{d}, ...
        'Color',colors{d},'MarkerFaceColor',colors{d},'LineWidth',1.3);
end
h = fill([tPred; flipud(tPred)], [yLow'; flipud(yHigh')], 'c','EdgeColor','none');
set(h,'FaceAlpha',0.2);
plot(tPred, yMean,'-k','LineWidth',2);
xlabel('Years since infection');
ylabel('Annual net inflow to subclinical (rate)');
legend([labels, {'95% posterior band','Posterior mean'}],'Location','northeast');
title('Fit to annual net inflow into subclinical (simple TB model)');

%% H) Save samples ---------------------------------------------------------
save('postSamples_simpleTB_uniform_adaptcov.mat', ...
     'chainThin','tData','yObs','x0_list','tPred','Ypred');

%% --- helper functions ---------------------------------------------------
function lj = logJacobian(u, lb, ub)
    sig = 1./(1+exp(-u));
    lj  = sum( log(ub-lb) + log(sig) + log(1-sig) );
end

% ---- log-posterior (annual net inflow target) ----------------------------
function lp = logPosterior_simple(p, tGridYears, yObs, x0)
% p: [alpha beta gamma delta epsilon zeta eta theta mu]
    N0       = x0(1);
    r_model  = annual_inflow_rate_simple(p, x0, max(tGridYears));  
    r_select = r_model(tGridYears);                                 

    cases_hat = yObs(:) * N0;
    pi_hat    = max(min(r_select(:), 1-1e-9), 1e-9);                

    % Binomial log-likelihood with fractional counts (quasi-likelihood)
    ll = cases_hat .* log(pi_hat) + (N0 - cases_hat) .* log(1 - pi_hat);

    lp = sum(ll);  
end

% ---- compute annual inflow rate vector (simple model) --------------------
function r = annual_inflow_rate_simple(p, x0, T_max)
% r(y) = (1/N0) * ∫_{y-1}^{y} [gamma*I(t) + epsilon*M(t)] dt
    N0   = x0(1);
    opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
    sol  = ode45(@(t,x) simpleTBmodel(t, x, p), [0 T_max], x0, opts);

    npt_per_year = 200;
    r = nan(T_max,1);
    gamma  = p(3);
    epsilon= p(5);
    for y = 1:T_max
        tt = linspace(y-1, y, npt_per_year);
        XY = deval(sol, tt);
        I  = XY(1,:); 
        M  = XY(2,:);
        Fin     = gamma.*I + epsilon.*M;  
        inflow  = trapz(tt, Fin);        
        r(y)    = inflow / N0;            
    end
end

% ---- simple TB ODE (states: [I; M; S; C; R; D]) --------------------------
function dx = simpleTBmodel(~, x, p)
% Params: p = [alpha beta gamma delta epsilon zeta eta theta mu]
    alpha = p(1);  beta  = p(2);  gamma = p(3);
    delta = p(4);  epsl  = p(5);  zeta  = p(6);
    eta   = p(7);  theta = p(8);  mu    = p(9);

    I = x(1); M = x(2); S = x(3); C = x(4);

    dx = zeros(6,1);
    % dI/dt = -(alpha + beta + gamma) I
    dx(1) = -(alpha + beta + gamma) * I;

    % dM/dt = beta I - (delta + epsilon) M + zeta S
    dx(2) = beta*I - (delta + epsl)*M + zeta*S;

    % dS/dt = gamma I + epsilon M - (eta + zeta) S
    dx(3) = gamma*I + epsl*M - (eta + zeta)*S + theta*C;

    % dC/dt = eta S - (theta + mu) C
    dx(4) = eta*S - (theta + mu)*C;

    % dR/dt = alpha I + delta M
    dx(5) = alpha*I + delta*M;

    % dD/dt = mu C
    dx(6) = mu*C;
end
