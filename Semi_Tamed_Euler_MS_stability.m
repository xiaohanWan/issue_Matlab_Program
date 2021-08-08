% STAB Mean-Square stability test for semi tamed E-M
%
% SDE is dX = -lambda*X-sigma*X^5 dt + mu*X dW, X(0) = Xzero,
%       where lambda and mu are constant and Xzero = 1.
%

randn('state', 103)
T = 20; M = 5000; Xzero = 1;
ltype = {'b-','r--', 'm-'};

lambda = 2; mu = sqrt(2); sigma = 1;

% dt = 2^(-3);
% Num = T/dt;
% dW = sqrt(2^(-3))*randn(1 ,Num);
% W = cumsum(dW);
% Xtrue = (Xzero*exp((lambda - 0.5*mu^2)*([dt:dt:T]) + mu*W)).^2;
% plot([0:dt:T], [Xzero, Xtrue], 'm*-'), hold on

for k = 1:3
    Dt = 2^(-1-k);
    N = T/Dt;
    Xms = zeros(1,N); Xtemp = Xzero*ones(M,1);
    for j = 1:N
        Winc = sqrt(Dt)*randn(M,1);
        tamed_term = 1 + abs( - sigma*Xtemp.^5)*Dt;
        Xtemp = Xtemp - Dt*lambda*Xtemp -  (Dt*sigma*Xtemp.^5)./tamed_term + mu*Xtemp.*Winc;
        Xms(j) = mean(Xtemp.^2);
    end
    semilogy([0:Dt:T], [Xzero, Xms], ltype{k}, 'Linewidth',1), hold on
end
legend('\Delta t = 1/4', '\Delta t = 1/8', '\Delta t = 1/16')
title('Mean-Square stability')
ylabel('E[X^2]', 'FontSize', 12), axis([0, T, 1e-20, 1e+20]), hold off