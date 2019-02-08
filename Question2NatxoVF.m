%%% DEVELOPMENT ECONOMICS. HOMEWORK 2 %%%
%%% QUESTION 2 %%%
%% 2.0. PARAMETERS AND CONSUMPTION LEVELS FOR PART 1 %%
% Clear and set the seed
clc; clear all; close all;
rng(7777777)
% Determine global Beta and Kappa
global BETA KAPPA
% Set the parameters
N = 1000; 
M=12;
T=40;
N = 1000; 
TM = T*M;
var_u = 0.2;
var_e = 0.2;
sd_e = sqrt(var_e);
sd_u = sqrt(var_u);
beta = 0.99^(1/M);
KAPPA = 0.6*3.5*(28.5*30/7)^(-2); % Choose kappa from FOCS considering nu is 1
% Determine Beta matrix
BETA = beta.^(12:TM+11);
% Permanent level of consumption and labor
z_c = repmat(exp(-var_u/2+sd_u.*randn(N,1)),1,TM);
z_l = repmat(exp(-var_u/2+sd_u.*randn(N,1)),1,TM);
% Ideosincratic nonseasonal components
inc_c = kron(exp(-var_e/2+sd_e.*randn(N,T)),ones(1,M));
inc_l = kron(exp(-var_e/2+sd_e.*randn(N,T)),ones(1,M));
% Table 1 Deterministic seasonal component
gm_med = [-0.147, -0.37, 0.141, 0.131, 0.09, 0.058, 0.036, 0.036, 0.036, 0.002, -0.033, -0.082];
gm_low = [-0.073, -0.185, 0.071, 0.066, 0.045, 0.029, 0.018, 0.018, 0.018, 0.001, -0.017, -0.041];
gm_high = [-0.293, -0.739, 0.282, 0.262, 0.18, 0.116, 0.072, 0.072, 0.072, 0.004, -0.066, -0.164];
% Table 2 Stochastic seasonal component
var_m_med = [0.085, 0.068, 0.290, 0.283, 0.273, 0.273, 0.239, 0.205, 0.188, 0.188, 0.171, 0.137];
var_m_high = [0.171, 0.137, 0.580, 0.567, 0.546, 0.546, 0.478, 0.410, 0.376, 0.376, 0.341, 0.273];
var_m_low = [0.043, 0.034, 0.145, 0.142, 0.137, 0.137, 0.119, 0.102, 0.094, 0.094, 0.085, 0.068];
% Deterministic seasonal component of consumption
e_gm_low = exp(repmat(gm_low,N,T));
e_gm_med = exp(repmat(gm_med,N,T));
e_gm_high = exp(repmat(gm_high,N,T));
% Inverse of the deterministic seasonal component
inv_e_gm_low_neg = 1./e_gm_low;
inv_e_gm_mid_neg = 1./e_gm_med;
inv_e_gm_hig_neg = 1./e_gm_high;
% Generate positive correlation of the stochastic seasonal components of consumption and labor
for i = 1:N
    for k = 1:12
        mu_low = [0 0];
        mu_med = [0 0];
        mu_high = [0 0];
        VAR_low = [var_m_low(k),0.03;0.03,var_m_low(k)]; %0.03 is the biggest coefficient that allows for matrix inversion
        VAR_med = [var_m_med(k),0.03;0.03,var_m_med(k)];
        VAR_high = [var_m_high(k),0.03;0.03,var_m_high(k)];
        LN_low = mvnrnd(mu_low,VAR_low,T);
        LN_med = mvnrnd(mu_med,VAR_med,T);
        LN_high = mvnrnd(mu_high,VAR_high,T);
        
        for j = 1:(T)
            ssc_low_c_pos(i,k+(j-1)*M)=exp(-var_m_low(k)/2+LN_low(j,1));
            ssc_med_c_pos(i,k+(j-1)*M)=exp(-var_m_med(k)/2+LN_med(j,1));
            ssc_high_c_pos(i,k+(j-1)*M)=exp(-var_m_high(k)/2+LN_high(j,1));
            ssc_low_l_pos(i,k+(j-1)*M)=exp(-var_m_low(k)/2+LN_low(j,2));
            ssc_med_l_pos(i,k+(j-1)*M)=exp(-var_m_med(k)/2+LN_med(j,2));
            ssc_high_l_pos(i,k+(j-1)*M)=exp(-var_m_high(k)/2+LN_high(j,2));
        end
    end
end
% Generate negative correlated stochastic seasonal components of consumption and labor
for i = 1:N
    for k = 1:12
        mu_low = [0 0];
        mu_med = [0 0];
        mu_high = [0 0];
        VAR_low = [var_m_low(k),-.03;-.03,var_m_low(k)];
        VAR_med = [var_m_med(k),-.03;-.03,var_m_med(k)];
        VAR_high = [var_m_high(k),-.03;-.03,var_m_high(k)];
        LN_low = mvnrnd(mu_low,VAR_low,T);
        LN_med = mvnrnd(mu_med,VAR_med,T);
        LN_high = mvnrnd(mu_high,VAR_high,T);
        
        for j = 1:(T)
            ssc_low_c_neg(i,k+(j-1)*M)=exp(-var_m_low(k)/2+LN_low(j,1));
            ssc_med_c_neg(i,k+(j-1)*M)=exp(-var_m_med(k)/2+LN_med(j,1));
            ssc_high_c_neg(i,k+(j-1)*M)=exp(-var_m_high(k)/2+LN_high(j,1));
            ssc_low_l_neg(i,k+(j-1)*M)=exp(-var_m_low(k)/2+LN_low(j,2));
            ssc_med_l_neg(i,k+(j-1)*M)=exp(-var_m_med(k)/2+LN_med(j,2));
            ssc_hig_l_neg(i,k+(j-1)*M)=exp(-var_m_high(k)/2+LN_high(j,2));
        end
    end
end
% Household consumptions
C4_low = z_c.*e_gm_low.*ssc_low_c_pos.*inc_c;
C4_med = z_c.*e_gm_med.*ssc_med_c_pos.*inc_c;
C4_hig = z_c.*e_gm_high.*ssc_high_c_pos.*inc_c;
C4_low_neg_inv = z_c.*e_gm_low.*ssc_low_c_neg.*inc_c;
C4_med_neg_inv = z_c.*e_gm_med.*ssc_med_c_neg.*inc_c;
C4_high_neg_inv = z_c.*e_gm_high.*ssc_high_c_neg.*inc_c;
C4_low_neg = z_c.*inv_e_gm_low_neg.*ssc_low_c_neg.*inc_c;
C4_med_neg = z_c.*inv_e_gm_mid_neg.*ssc_med_c_neg.*inc_c;
C4_high_neg = z_c.*inv_e_gm_hig_neg.*ssc_high_c_neg.*inc_c;
C1 = z_c.*inc_c;

% Household labor supplies
L4_low = 115*z_l.*e_gm_low.*ssc_low_l_pos.*inc_l;
L4_med = 115*z_l.*e_gm_med.*ssc_med_l_pos.*inc_l;
L4_high = 115*z_l.*e_gm_high.*ssc_high_l_pos.*inc_l;
L4_low_neg_inv = 115*z_l.*inv_e_gm_low_neg.*ssc_low_l_neg.*inc_l;
L4_med_neg_inv = 115*z_l.*inv_e_gm_mid_neg.*ssc_med_l_neg.*inc_l;
L4_high_neg_inv = 115*z_l.*inv_e_gm_hig_neg.*ssc_hig_l_neg.*inc_l;
L4_low_neg = 115*z_l.*e_gm_low.*ssc_low_l_neg.*inc_l;
L4_med_neg = 115*z_l.*e_gm_med.*ssc_med_l_neg.*inc_l;
l4_high_neg = 115*z_l.*e_gm_high.*ssc_hig_l_neg.*inc_l;
L1 = 115*z_l.*inc_l;

%% 2.A) POSITIVE CORRELATION
% Total effects
for i = 1:N
    E_total_low(i) = wg(C4_low,L4_low,C1,L1,i);
    E_total_med(i) = wg(C4_med,L4_med,C1,L1,i);
    E_total_high(i) = wg(C4_hig,L4_high,C1,L1,i);
end
% Consumption effects
for i = 1:N
    E_c_low(i) = wg(C4_low,L4_low,C1,L4_low,i);
    E_c_med(i) = wg(C4_med,L4_med,C1,L4_med,i);
    E_c_high(i) = wg(C4_hig,L4_high,C1,L4_high,i);
end
% Labour effects
for i = 1:N
    E_l_low(i) = wg(C1,L4_low,C1,L1,i);
    E_l_med(i) = wg(C1,L4_med,C1,L1,i);
    E_l_hig(i) = wg(C1,L4_high,C1,L1,i);
end

disp('Part (a) Result')
MeanWelfareGains = {'Total effects';'Effects of consumption';'Effects of labor'};
Low = [mean(E_total_low);mean(E_c_low);mean(E_l_low)];
Middle = [mean(E_total_med);mean(E_c_med);mean(E_l_med)];
High = [mean(E_total_high);mean(E_c_high);mean(E_l_hig)];
Q2A = table(MeanWelfareGains,Low,Middle,High)

%% 2.B) NEGATIVE CORRELATION
% Consumption behaves as usual but labor is adjusted to find negative
% correlation
% Effects of both
for i = 1:N
    E2_total_low(i) = wg(C4_low_neg,L4_low_neg,C1,L1,i);
    E2_total_med(i) = wg(C4_med_neg,L4_med_neg,C1,L1,i);
    E2_total_high(i) = wg(C4_high_neg,l4_high_neg,C1,L1,i);
end

% Effects of consumption
for i = 1:N
    E2_c_low(i) = wg(C4_low_neg,L4_low_neg,C1,L4_low_neg,i);
    E2_c_med(i) = wg(C4_med_neg,L4_med_neg,C1,L4_med_neg,i);
    E2_c_high(i) = wg(C4_high_neg,l4_high_neg,C1,l4_high_neg,i);
end

% Effects of labor
for i = 1:N
    E2_l_low(i) = wg(C1,L4_low_neg,C1,L1,i);
    E2_l_med(i) = wg(C1,L4_med_neg,C1,L1,i);
    E2_l_high(i) = wg(C1,l4_high_neg,C1,L1,i);
end

disp('Part (b) Result')
MeanWelfareGains = {'Effects of both';'Effects of consumption';'Effects of labor'};
Low2 = [mean(E2_total_low);mean(E2_c_low);mean(E2_l_low)];
Middle2 = [mean(E2_total_med);mean(E2_c_med);mean(E2_l_med)];
High2 = [mean(E2_total_high);mean(E2_c_high);mean(E2_l_high)];
Q2B = table(MeanWelfareGains,Low2,Middle2,High2)

%% The function to compute welfare gains
function g = wg(a,b,c,d,i)
global BETA KAPPA
tmp = @(g) abs(sum(BETA.*(log(a(i,:)*(1+g))-KAPPA*b(i,:).^(1+1)/(1+1))-BETA.*(log(c(i,:))-KAPPA*d(i,:).^(1+1)/(1+1))));
g = fminbnd(tmp,-30,30);
end