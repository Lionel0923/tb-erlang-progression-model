function dx = simpleTBmodel(t, x, p)
    % unpack parameters
    alpha   = p(1);
    beta    = p(2);
    gamma   = p(3);
    delta   = p(4);
    epsilon = p(5);
    zeta    = p(6);
    eta     = p(7);
    theta   = p(8);
    mu      = p(9);

    % unpack states
    I = x(1);
    M = x(2);
    S = x(3);    % single subclinical compartment
    C = x(4);
    R = x(5);
    D = x(6);

    % derivatives
    dI = -(alpha + beta + gamma) * I;
    dM =  beta * I - (delta + epsilon) * M + zeta * S;
    dS =  gamma * I + epsilon * M - (zeta + eta) * S;
    dC =  eta   * S - (theta + mu) * C;
    dR =  alpha * I + delta * M;
    dD =  mu    * C;

    % pack
    dx = [dI; dM; dS; dC; dR; dD];
end
