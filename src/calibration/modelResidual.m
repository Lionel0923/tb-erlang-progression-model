function res = modelResidual(p, tData, yObs, x0, k)
    % integrate from t=0
    sol = ode45(@(t,x) erlangODE(t,x,p,k), [0 tData(end)], x0);
    % extract modeled C(t) at each observation year
    Cmodel = deval(sol, tData, 3 + k);
    % residual = model – data
    res = Cmodel(:) - yObs(:);
end
