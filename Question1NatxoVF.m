%%% DEVELOPMENT ECONOMICS. HOMEWORK 2 %%%
%%% QUESTION 1 %%%
%% 1.1.0. PARAMETERS AND CONSUMPTION LEVELS FOR PART 1 %%
% Clear and set the seed
clc; clear all; close all;
rng(7777777)
% Determine global Beta and Eta variables 
global BETA ETA
% Set the main parameters
M=12;
T=40;
N = 1000; 
TM = T*M; 
var_u = .2; 
var_e = .2;
sd_e = sqrt(var_e);
sd_u = sqrt(var_u);
beta = .99^(1/12);
eta1 = 1;
eta2 = 2;
eta3 = 4;
% Determine Beta matrix
BETA = beta.^(12:TM+11);
% Permanent level of consumption
Z = repmat(exp(-var_u/2+sd_u.*randn(N,1)),1,TM);
% Table 1 Deterministic seasonal component
gm_med = [-0.147, -0.37, 0.141, 0.131, 0.09, 0.058, 0.036, 0.036, 0.036, 0.002, -0.033, -0.082];
gm_low = [-0.073, -0.185, 0.071, 0.066, 0.045, 0.029, 0.018, 0.018, 0.018, 0.001, -0.017, -0.041];
gm_high = [-0.293, -0.739, 0.282, 0.262, 0.18, 0.116, 0.072, 0.072, 0.072, 0.004, -0.066, -0.164];
% Table 2 Stochastic seasonal component
var_m_med = [0.085, 0.068, 0.290, 0.283, 0.273, 0.273, 0.239, 0.205, 0.188, 0.188, 0.171, 0.137];
var_m_high = [0.171, 0.137, 0.580, 0.567, 0.546, 0.546, 0.478, 0.410, 0.376, 0.376, 0.341, 0.273];
var_m_low = [0.043, 0.034, 0.145, 0.142, 0.137, 0.137, 0.119, 0.102, 0.094, 0.094, 0.085, 0.068];
% Deterministic seasonal component of consumption
e_gm_med = exp(repmat(gm_med,N,T));
e_gm_low = exp(repmat(gm_low,N,T));
e_gm_high = exp(repmat(gm_high,N,T));
% Ideosincratic nonseasonal component
inc = kron(exp(-var_e/2+sd_e.*randn(N,T)),ones(1,12));
% Stochastic seasonal component
SSC_low = exp(repmat(-var_m_low/2,N,T)+repmat(sqrt(var_m_low),N,T).*randn(N,TM));
SSC_med = exp(repmat(-var_m_med/2,N,T)+repmat(sqrt(var_m_med),N,T).*randn(N,TM));
SSC_high = exp(repmat(-var_m_high/2,N,T)+repmat(sqrt(var_m_high),N,T).*randn(N,TM));
% Total consumption matrices with all risk components
C3_low = Z.*e_gm_low.*inc;
C3_med = Z.*e_gm_med.*inc;
C3_high = Z.*e_gm_high.*inc;
% Total consumption matrix removing nonseasonal consumption risk
C2_low = Z.*e_gm_low;
C2_med = Z.*e_gm_med;
C2_high = Z.*e_gm_high;
% Total consumption removing seasonal consumption risk
C_2 = Z.*inc;
% Total consumption (with stochastic seasonal component)
C3 = {};
C3{1} = C3_low;
C3{2} = C3_med;
C3{3} = C3_high;
C2 = {};
C2{1} = C2_low;
C2{2} = C2_med;
C2{3} = C2_high;
ssc = {};
ssc{1} = SSC_low;
ssc{2} = SSC_med;
ssc{3} = SSC_high;
con_ssc = {};
con_norisk_ssc = {};
for i = 1:3
    for j = 1:3
        con_ssc{i,j} = C3{i}.*ssc{j};
        con_norisk_ssc{i,j} = C2{i}.*ssc{j};
    end
end
% Auxiliary for display
Eta = []; 
Low = []; 
Med = []; 
High = [];
Low_Low = [];
Low_Med = [];
Low_High = [];
Med_Low = [];
Med_Med = [];
Med_High = [];
High_Low = [];
High_Med = [];
High_High = [];

%% QUESTION 1.1.
%% QUESTION 1.a) and 1.b)
ETA = eta1
    for i = 1:N
        gl1_1(i) = wg(C3_low,C_2,i);
        gm1_1(i) = wg(C3_med,C_2,i);
        gh1_1(i) = wg(C3_high,C_2,i);
        gl2_1(i) = wg(C3_low,C2_low,i);
        gm2_1(i) = wg(C3_med,C2_med,i);
        gh2_1(i) = wg(C3_high,C2_high,i);
    end
    Low_a = [Low;mean(gl1_1);mean(gl2_1)];
    Med_a = [Med;mean(gm1_1);mean(gm2_1)];
    High_a = [High;mean(gh1_1);mean(gh2_1)];
    Low_sd = [Low;std(gl1_1);std(gl2_1)];
    Med_sd = [Med;std(gm1_1);std(gm2_1)];
    High_sd = [High;std(gh1_1);std(gh2_1)];
% Results
disp('Part 1. Mean welfare gains for Eta 1')
Break = {'Removing seasonal component';'Removing nonseasonal component';};
A_Q1AB = table(Break,Low_a,Med_a,High_a)
SD_Q1AB = table(Break,Low_sd,Med_sd,High_sd)

%% QUESTION 1.c)
% Comments on page

%% QUESTION 1.d)
% Eta 2
ETA = eta2
    for i = 1:N
        gl1_2(i) = wg(C3_low,C_2,i);
        gm1_2(i) = wg(C3_med,C_2,i);
        gh1_2(i) = wg(C3_high,C_2,i);
        gl2_2(i) = wg(C3_low,C2_low,i);
        gm2_2(i) = wg(C3_med,C2_med,i);
        gh2_2(i) = wg(C3_high,C2_high,i);
    end
    Low_a = [Low;mean(gl1_2);mean(gl2_2)];
    Med_a = [Med;mean(gm1_2);mean(gm2_2)];
    High_a = [High;mean(gh1_2);mean(gh2_2)];
    Low_sd = [Low;std(gl1_2);std(gl2_2)];
    Med_sd = [Med;std(gm1_2);std(gm2_2)];
    High_sd = [High;std(gh1_2);std(gh2_2)];

% Results
disp('Part 1. Mean welfare gains for Eta 2')
Break = {'Removing seasonal component';'Removing nonseasonal component';};
A_Q1D_2 = table(Break,Low_a,Med_a,High_a)
SD_Q1D_2 = table(Break,Low_sd,Med_sd,High_sd)
% Eta 4
ETA = eta3
    for i = 1:N
        gl1_3(i) = wg(C3_low,C_2,i);
        gm1_3(i) = wg(C3_med,C_2,i);
        gh1_3(i) = wg(C3_high,C_2,i);
        gl2_3(i) = wg(C3_low,C2_low,i);
        gm2_3(i) = wg(C3_med,C2_med,i);
        gh2_3(i) = wg(C3_high,C2_high,i);
    end
    Low_a = [Low;mean(gl1_3);mean(gl2_3)];
    Med_a = [Med;mean(gm1_3);mean(gm2_3)];
    High_a = [High;mean(gh1_3);mean(gh2_3)];
    Low_sd = [Low;std(gl1_3);std(gl2_3)];
    Med_sd = [Med;std(gm1_3);std(gm2_3)];
    High_sd = [High;std(gh1_3);std(gh2_3)];
% Results
disp('1. Mean welfare gains for Eta 4')
Break = {'Removing seasonal component';'Removing nonseasonal component';};
A_Q1D_4 = table(Break,Low_a,Med_a,High_a)
SD_Q1D_4 = table(Break,Low_sd,Med_sd,High_sd)
% Histograms comparing gains from different ETAs
hold on
histogram(gl2_1,25);
hold on
histogram(gl2_2,25);
hold on
histogram(gl2_3,25);
xlabel('Welfare gains distribution')
legend('?=1','?=2','?=4')
print('Question 1.d)','-dpng')

%% QUESTION 1.2.
for ETA = [eta1,eta2,eta3]
    for i = 1:N
        gll1(i) = wg(con_ssc{1,1},C_2,i);
        glm1(i) = wg(con_ssc{1,2},C_2,i);
        glh1(i) = wg(con_ssc{1,3},C_2,i);
        gml1(i) = wg(con_ssc{2,1},C_2,i);
        gmm1(i) = wg(con_ssc{2,2},C_2,i);
        gmh1(i) = wg(con_ssc{2,3},C_2,i);
        ghl1(i) = wg(con_ssc{3,1},C_2,i);
        ghm1(i) = wg(con_ssc{3,2},C_2,i);
        ghh1(i) = wg(con_ssc{3,3},C_2,i);
        gll2(i) = wg(con_ssc{1,1},con_norisk_ssc{1,1},i);
        glm2(i) = wg(con_ssc{1,2},con_norisk_ssc{1,2},i);
        glh2(i) = wg(con_ssc{1,3},con_norisk_ssc{1,3},i);
        gml2(i) = wg(con_ssc{2,1},con_norisk_ssc{2,1},i);
        gmm2(i) = wg(con_ssc{2,2},con_norisk_ssc{2,2},i);
        gmh2(i) = wg(con_ssc{2,3},con_norisk_ssc{2,3},i);
        ghl2(i) = wg(con_ssc{3,1},con_norisk_ssc{3,1},i);
        ghm2(i) = wg(con_ssc{3,2},con_norisk_ssc{3,2},i);
        ghh2(i) = wg(con_ssc{3,3},con_norisk_ssc{3,3},i);
    end
    
    Eta = [Eta;ETA;ETA];
    Low_Low = [Low_Low;mean(gll1);mean(gll2)]
    Low_Med = [Low_Med;mean(glm1);mean(glm2)];
    Low_High = [Low_High;mean(glh1);mean(glh2)];
    Med_Low = [Med_Low;mean(gml1);mean(gml2)];
    Med_Med = [Med_Med;mean(gmm1);mean(gmm2)];
    Med_High = [Med_High;mean(gmh1);mean(gmh2)];
    High_Low = [High_Low;mean(ghl1);mean(ghl2)];
    High_Med = [High_Med;mean(ghm1);mean(ghm2)];
    High_High = [High_High;mean(ghh1);mean(ghh2)];
end

disp('2. Welfare Gains')
Break = {'Removing seasonal component';'Removing nonseasonal component';};
Q2 = table(Low/Low,Low/Med,Low/High,Med/Low,Med/Med,Med/High,High/Low,High/Med,High/High)

%% Welfare gain function. Allows to compute the welfare gains (g) considering the consumption levels compared
function g = wg(a,b,c)
global BETA ETA
if ETA == 1
    tmp = @(g) abs(sum((BETA.*log(a(c,:)*(1+g)))-(BETA.*log(b(c,:)))));
else
    tmp = @(g) abs(sum((BETA.*((a(c,:)*(1+g)).^(1-ETA)/(1-ETA)))-(BETA.*(b(c,:).^(1-ETA)/(1-ETA)))));
end
g = fminbnd(tmp,-5,5);
end
