function [x, tau2, mse] = amp(A, y, P, u_g, v_g, x_true, inter_max)

% amp - This is an implementation of Approximate Message Passing algorithm
% described in Arian Maleki's thesis
% 
% Inputs:
% A         - measurement matrix
% y         - measurement vector
% lambda    - denoiser parameter in the noisy case
% x_true    - the uknown vector (here used to get mse accross iterations)
% amp0      - boolean indicator (true for noiseless case)
% inter_max - maximum number of amp iterations
%
% Outputs:
% x         - the last estimate
% tau2      - vector of effective variances accross interations
% mse       - vector of mean squared errors accross interations

% Example:
% amp(A, y, 0.5, x_true, true, 20)

% Author:   Osman Musa
% email:    osman.musa@nt.tuwien.ac.at
% Website:  https://www.nt.tuwien.ac.at/about-us/staff/osman-musa/
% Last revision: 01-Aug-2017

[M, N] = size(A);
delta = M/N;
tau2 = zeros(inter_max, 1);
mse = zeros(inter_max, 1);

z(:,1) = y;
x(:,1) = zeros(N,1);
tau2(1) = 1/M*(norm(z(:,1),2)^2);
mse(1) = norm(x_true - x,2).^2/N;


for t=2:inter_max
    
%     x(:,t) = wthresh(A'*z(:,t-1) + x(:,t-1),'s', deniser_parameter);
    r = A'*z(:,t-1) + x(:,t-1);
    [u_post, v_post] = BG_MMSE_denoiser(r, tau2(t-1), P, u_g, v_g, N);
    x(:,t) = u_post;
    
%     eta_prime = 1/N * nnz(wthresh(A'*z(:,t-1) + x(:,t-1),'s', deniser_parameter));
    eta_prime = v_post / tau2(t-1);
    z(:,t) = y - A*x(:,t) + 1/delta * z(:,t-1) * eta_prime;

 
    tau2(t) = 1/M*(norm(z(:,t),2)^2);
    mse(t) = norm(x_true - x(:,t),2).^2/N; 

    
%     fprintf('AMP MSE = %f \n', 10*log10(ampsim_MSE(i+1))); % optional printing
end

% fprintf('MSE = %f \n', 10*log10(1/N*norm(x_true - x,2).^2)); % optional printing
