function lp = logPosterior_simple(p, tData, yObs, x0)
    % Solve the ODE from t=0 to tData(end) using simpleTBmodel
    sol = ode45(@(t,x) simpleTBmodel(t, x, p), [0 tData(end)], x0);

    % Evaluate all states at each observation time
    % Xmod is 6×numel(tData): rows = [I; M; S; C; R; D]
    Xmod = deval(sol, tData);

    % Extract subclinical counts (single S compartment)
    Smod = Xmod(3, :)';    

    % Convert to proportion of the initial at-risk population
    pred = Smod ./ x0(1);

    % Ensure observed proportions is a column
    yObs = yObs(:);

    % Gaussian log‐likelihood on proportions
    sigma2 = 5e-6;
    ll = -0.5 * sum( (yObs - pred).^2 ) / sigma2;

    lp = ll;
end
