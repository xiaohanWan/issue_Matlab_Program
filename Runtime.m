randn('state', 1)
beta = 1; Xzero = 1;         % problem parameters
T = 1; N = 2^(19); dt = T/N;                    %
M = 1000;                     % number of paths sampled
lambda = 1;
index = 5;
R = [1;2^1;2^2;2^3;2^4; 2^5; 2^6; 2^7;2^8;2^9;2^10;2^11;2^12;2^13;2^14];                    % Milstein stepsizes are R*dt
runtime = zeros(13,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dW = sqrt(dt)*randn(M, N);                      % Brownian increments
Xmil = zeros(M, 15);                             % preallocate array
for p = 1:15
    Dt = R(p)*dt; L = N/R(p);                   % L timesteps of size Dt = R dt
    Xtemp = Xzero*ones(M, 1);
    mean_runtime = zeros(10,1);
    for iterator = 1:10
        tic;
    for j = 1:L
        Winc = sum( dW(:, R(p)*(j-1)+1:R(p)*j), 2);
        semi_tamed_term = 1 + abs(lambda*Xtemp.^index)*Dt;
        Xtemp = Xtemp + 2*Xtemp*Dt - (lambda*Xtemp.^index*Dt)./semi_tamed_term + beta*Xtemp.*Winc...
            + 0.5*beta^2*Xtemp.*(Winc.^2 - Dt);
    end
    mean_runtime(iterator) = toc;
    end
    runtime(p) = mean(mean_runtime);
    Xmil(:, p) = Xtemp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
Xref = Xmil(:, 1);
Xerr = abs(Xmil(:, 2:15) - repmat(Xref, 1, 14));
mean_error = sqrt(mean(Xerr.^2));
runtime = runtime(2:15);
loglog(runtime(3:8), mean_error(3:8),'b*-');
