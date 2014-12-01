import_data_GEN; %takes care of formatting, an autogenned file

%% Simple Calculations

y = ics_all;
x = [ones(length(y), 1), pago_r_all, pexp_r_all, news_r_all, ratex_r_all, px1_mean_all ...
    px1_var_all, px1_std_all, govt_r_all, dur_r_all, veh_r_all, vehrn_np_all];

n = length(y);
p = size(x);
p = p(2);

b = inv( transp(x) * x ) * transp(x) * y; %#ok<MINV>

u = x * b;

e = y - u;
scatter(1:length(e), e)
fprintf('STD DEV of e: %f \n', std(e));
fprintf('MEAN of e: %f \n', mean(e));

s2 = transp(e) * e / (n - p - 1);

var = inv(transp(x) * x);

%% ANOVA
%creating the anova table

SSR = transp(b) * transp(x) * y - n * mean(y)^2;
SSE = transp(y - x * b) * (y - x * b);
SST = sum((y - mean(y)).^2);
fprintf('SSE: %f\n', SSE);
fprintf('SSR: %f\n', SSR);
if (abs(SST - (SSR + SSE)) > 1)
    fprintf('SUMS DO NOT MATCH, ERROR\n');
    fprintf('%f != %f\n', SST, SSR + SSE);
else
    fprintf('SST = SSR + SSE, proceeding\n');
    fprintf('%f = %f\n', SST, SSR + SSE);
end

MSR = SSR/p;
MSE = SSE/(n-p-1);
fprintf('MSR: %f, MSE: %f\n', MSR, MSE);

F = MSR/MSE;
fprintf('F: %f\n', F);
pval = fpdf( F, n-p-1, p);
fprintf('P value: %f\n', pval);
if (pval > 0.05)
    fprintf('Fail to reject null hypothesis, as our p value is so large\n');
else
    fprintf('Reject the null hypothesis, as our p value is so small\n');
end

%% Testing Important Regressors
% testing if financial situation compared to a year ago is an important regressor

T = b(2) / sqrt(s2 * var(2,2));

pval = 2*tcdf( T, n - p -1, 'upper');


fprintf('P-value when testing financial situation compared to a year ago %f\n', pval);


%% Prediction 
% make a prediction given a relatively reasonable data set

test_val = [1; 112; 117; 78; 55; 3.9; 40; 5; 100; 127; 130; 8];
t_val = tinv(.95, n-p-1)* sqrt(s2) * sqrt(transp(test_val) * inv(1+transp(x) * x) * test_val);

pred_int = [transp(test_val) * b - t_val, transp(test_val) * b + t_val];
fprintf('Prediction interval: %f - %f\n', pred_int(1), pred_int(2));
