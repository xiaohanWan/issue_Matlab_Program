% Tamed-Euler Test stong convergence of tamed Euler : vectorized
%
% Solves dX = 2X - lambda*X^5 dt + beta*X dW, X(0) = Xzero,
%
% Discretized Brownian path over [0, 1] has dt = 2^(-11)
% uses timesteps 128*dt, 64*dt, 32*dt, 16*dt (also dt for reference)
%
% Examines strong convergence at T = 1: E|X_L - X_T|.
% Code is vectorized: all paths computed simultaneously.

randn('state', 101)
lambda = 1; beta = 1; Xzero = 1;         % problem parameters
index = 5;
T = 1; N = 2^(15); dt = T/N;                    %
M = 5000;                                        % number of paths sampled
R = [1; 16; 32; 64; 128; 256];                       % Milstein stepsizes are R*dt

dW = sqrt(dt)*randn(M, N);                      % Brownian increments
Xtamed_eul = zeros(M, 6);                       % preallocate array
for p = 1:6
    Dt = R(p)*dt; L = N/R(p);                   % L timesteps of size Dt = R dt
    Xtemp = Xzero*ones(M, 1);
    for j = 1:L
        Winc = sum( dW(:, R(p)*(j-1)+1:R(p)*j), 2);
        semi_tamed_term = 1 + abs(lambda*Xtemp.^index)*Dt;
        Xtemp = Xtemp + 2*Xtemp*Dt - (lambda*Xtemp.^index*Dt)./semi_tamed_term + beta*Xtemp.*Winc...
            + 0.5*beta^2*Xtemp.*(Winc.^2 - Dt);
    end
    Xtamed_eul(:, p) = Xtemp;
end

Xref = Xtamed_eul(:, 1);
Xerr = abs(Xtamed_eul(:, 2:6) - repmat(Xref, 1, 5));
mean(Xerr);
Dtvals = dt*R(2:6);

loglog(Dtvals, mean(Xerr), 'b*-'), hold on
loglog(Dtvals, Dtvals, 'r--'), hold off
axis([1e-4 1 1e-4 1])
xlabel('\Delta t')
ylabel('Sample average of |X(T) - X_L|')
title('Tamed Euler strong convergence', 'Fontsize', 10)
