function dx = erlangODE(t, x, p, k)
    % unpack parameters
    alpha   = p(1);
    beta    = p(2);
    gamma   = p(3);
    delta   = p(4);
    epsilon = p(5);
    kappa   = p(6);
    zeta    = p(7);
    theta   = p(8);
    mu      = p(9);

    % unpack states
    I = x(1);
    M = x(2);
    S = x(3 : 2 + k);      % S1…Sk
    C = x(3 + k);
    R = x(4 + k);
    D = x(5 + k);

    % derivatives
    dI = -(alpha + beta + gamma)*I;
    dM = beta*I - (delta + epsilon)*M + zeta*sum(S);

    dS = zeros(k,1);
    dS(1) = gamma*I + epsilon*M - (kappa + zeta)*S(1);
    for j = 2:k
        if j < k
            dS(j) = kappa*S(j-1) - (kappa + zeta)*S(j);
        else
            % j == k
            dS(j) = kappa*S(j-1) - (kappa + zeta)*S(j) + theta*C;
        end
    end

    dC = kappa*S(k) - (theta + mu)*C;
    dR = alpha*I + delta*M;
    dD = mu*C;

    dx = [dI; dM; dS; dC; dR; dD];
end
