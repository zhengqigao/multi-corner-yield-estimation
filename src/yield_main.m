% Last modifed by zqgao on Dec. 8, 2019.
% The code for TCAD paper 'Efficient Parametric Yield Estimation Over 
% Multiple Process Corners via Bayesian Inference Based on Bernoulli Distribution'

close all;clear all;
data_cell = importdata('data_cell.mat');
res_cell = importdata('res_cell.mat');
experi = 1:1:30;% expri repeated times, should less than size of res_cell
N_list = 50:50:250;

K=size(data_cell{1},1);
ERR_cell = cell(length(experi),length(N_list));
err_bibd = zeros(K,length(N_list));
MCERR_cell = cell(length(experi),length(N_list));
err_mc = zeros(K,length(N_list));
for j =1:length(experi)
    cur_experi = experi(j);
    for i = 1:length(N_list)
        cur_N = N_list(i);
        data = data_cell{cur_experi};data = data(:,1:cur_N);
        res = res_cell{cur_experi};
        esti = BI_BD(data);
        if isnan(esti)
            sprintf('error occurs in algorithim');
        end
        if size(esti,1)==1; esti = esti';end
        errk = abs(esti-res)./res;
        ERR_cell{j,i} = errk;
        err_bibd(:,i) = err_bibd(:,i) + errk/length(experi);
        MCERR_cell{j,i} = abs(sum(data,2)/cur_N - res)./res;
        err_mc(:,i) = err_mc(:,i) + abs(sum(data,2)/cur_N - res)./res/length(experi);
    end
end
figure;
plot(N_list, sum(err_bibd),'r*-','linewidth',1.5);hold on ;
plot(N_list, sum(err_mc),'b*-','linewidth',1.5);
legend('BI-BD', 'MC');
disp(err_bibd); % BI-BD error
disp(err_mc); % MC error