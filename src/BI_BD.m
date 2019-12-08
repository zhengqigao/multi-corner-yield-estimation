function   beta_star = BI_BD(data)
% BI-BD algorithim, edited by Zhengqi Gao, Sep 8th.
% data[matrix R^{K*N}]:corner number * total samples
% esti[vector R^{K*1}]: yield estimation result
%##################################################
[K, N] = size(data);
Mk = sum(data,2);
MAX_ITER1 = 3000;
MAX_ITER2 = 3000;
%EPS1 = 1/N * 0.1;
%EPS2 = 1/N * 0.1;
EPS1 = 1/N * 0.1;
EPS2 = 1/N * 0.1;

EPS1 = 1/N * 0.1;
EPS2 = 1/N * 0.1;

% if N <= 1000
%     EPS1 = 1e-4;
%     EPS2 = 1e-3;
% else
%     EPS1 = 1e-5;
%     EPS2 = 1e-4;
% end

% MAX_ITER1 = 5000;
% MAX_ITER2 = 5000;
% EPS1 = 1e-5;
% EPS2 = 1e-5;


info = [];
% intialize mu, Sigma
mu = mean(Mk/N)*ones(K,1);
% mu = 0.1 * ones(K,1);
ini_beta = mu(1);
ratio = (N*ini_beta^2-2*Mk(1)*ini_beta+Mk(1))/(ini_beta*(1-ini_beta))^2;
Sigma = 1/ratio* (0.5*ones(K,K)+diag(0.5*ones(K,1)));
beta_star = zeros(K,1) * nan;
out_count = 0;
% begin algorithim
while out_count < MAX_ITER1
    beta = zeros(K,1);
    % always initialize near Mk, but a little larger
    for kk=1:1:K
        if N == Mk(kk)
            %disp('[Warning]: no suitable inital guess...');
            difff = -Mk(kk) + mu(kk)*N;
            %pause();
        else
            if N-1==Mk(kk)
                difff = 0.5;
            else  
                difff = randperm(N-Mk(kk)-1,1);
            end
        end
        beta(kk)= (Mk(kk)+difff)/N;
    end
    ins_count = 0; beta_opt = [];H_opt = [];
    while ins_count < MAX_ITER2
        b = (N*beta-Mk)./beta./(1-beta);
        B = diag((N.*beta.^2-2.*Mk.*beta+Mk)./(beta.*(1-beta)).^2);
        grad = inv(Sigma)*(beta-mu)+b;
        H = inv(Sigma)+B;
        beta_new = beta - inv(H)*grad;
        if abs(norm(beta_new - beta,2))/norm(beta,2) < EPS1
            beta_opt = beta_new;H_opt = H; 
            info = [info ' inside breaks successfully'];
            break;
        else
            ins_count = ins_count + 1;
            beta = beta_new;
        end
    end
    mu_new = beta_opt; Sigma_new = inv(H_opt);
    if abs(norm(mu_new - mu,2))/norm(mu,2) < EPS2
        beta_star = mu_new;
        info = [info ' outside breaks successfully'];
        break;
    else
        out_count = out_count + 1;
        mu = mu_new; Sigma = Sigma_new;
    end
end

