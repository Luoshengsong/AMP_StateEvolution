function [ se_tau2, se_mse ] = state_evolution(tau2, delta, rho, u_g, v_g, sigmaw2)

% state_evolution - this function outputs state evolution predictions
%
% Inputs:
% tau2      - the effective noise variance from the previous interation
% delta     - the measurement ratio (M/N)
% rho       - the probability of a nonzero component, i.e., sparsity ratio (K/N)
% sigmaw2   - the variance of the noise
% lambda    - denoiser parameter in the noisy case
%
% Outputs:
% tau2      - the effective noise variance
% mse       - the mean squared error
%
% Author:   Osman Musa
% email:    osman.musa@nt.tuwien.ac.at
% Website:  https://www.nt.tuwien.ac.at/about-us/staff/osman-musa/
% Last revision: 01-Aug-2017


% % Gauss part
% f1 = @(x,z) (rho.*reshape(normpdf(x(:),0,1),size(x)) ) ...
%             .* reshape(normpdf(z(:),0,sqrt(tau2)),size(z)) ...
%             .* (wthresh(x + z,'s',deniser_parameter) - x).^2;


% Gauss part
f1 = @(x,z) (rho.*reshape(normpdf(x(:),0,1),size(x)) ) ...
            .* reshape(normpdf(z(:),0,sqrt(tau2)),size(z)) ...
            .* (BG_MMSE_denoiser_se(x + z, tau2, rho, u_g, v_g)- x).^2;

% % Dirac part        
% f2 = @(z) (1-rho) ...
%             .* reshape(normpdf(z(:),0,sqrt(tau2)),size(z)) ...
%             .* (wthresh( z,'s',deniser_parameter)).^2;        

% Dirac part        
f2 = @(z) (1-rho) ...
            .* reshape(normpdf(z(:),0,sqrt(tau2)),size(z)) ...
            .* (BG_MMSE_denoiser_se(0 + z, tau2, rho, u_g, v_g)- 0).^2;        
         
lim = inf;

se_mse = (integral2(f1,-lim,lim,-lim,lim) + integral(f2,-lim,lim));


se_tau2 = sigmaw2 + 1/delta * se_mse;   


end

function u_post = BG_MMSE_denoiser_se(u, v, P, u_g, v_g)

    % prior: (1-P) * delta_0 + P * N(u_g, v_g)

    EXP_MAX = 80;
    EXP_MIN = -80;
%     ug = u_g * ones(N, 1);
    ug = u_g;
    vg = v_g;
    % p1
    a = sqrt((v + vg) / v);
    b = 0.5 * ((u - ug).^2 / (v + vg) - (u.^2) / v);
    % set threshold
    b(b > EXP_MAX) = EXP_MAX;
    b(b < EXP_MIN) = EXP_MIN;
    c = (1 - P) / P;
    p1 = 1 ./ (1 + a * exp(b) * c);
    % Gaussian addition
    v1 = (vg^(-1) + v^(-1))^(-1);
    u1 = v1 * (vg^(-1) * ug + v^(-1) * u);
    % post u and v
    u_post = p1 .* u1;
    % v_post = mean(((p1 - p1.^2) .* (u1.^2) + p1 * v1));
end
