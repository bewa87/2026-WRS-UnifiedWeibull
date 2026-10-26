%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of Random Data Sets %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clearing and Closing

clear all;
close all;
clc;

% Reproducibility

rng(20261024, "twister");

% Parameters

N       = 1000;  % Sample Size
N_data  = 1000; % Number of Data Sets
k_shape = 1.5;  % Form Parameter
lambda  = 2.5;  % Scale Parameter

% Uniformly Distribution Random Numbers

U_weibull = rand(N, N_data);
U_inverse = rand(N, N_data);

% Weibull-Distribution

X_weibull = lambda.*((-log(U_weibull)).^(1/k_shape));

% Inverse-Weibull-Distribution

X_inverse = (lambda./((-log(U_inverse)).^(1/k_shape)));

% Control

fprintf('Weibull: Size of the Matrix = %d x %d\n', ...
        size(X_weibull,1), size(X_weibull,2));

fprintf('Inverse Weibull: Size of the Matrix = %d x %d\n', ...
        size(X_inverse,1), size(X_inverse,2));

fprintf('\nFirst Weibull-Data Set:\n');
disp(X_weibull(:,1)');

fprintf('\nFirst inverse Weibull-Data Set:\n');
disp(X_inverse(:,1)');

% Save Generated Data Sets

save('example02_weibull_datasets_1000.mat', 'X_weibull');
save('example02_inverse_weibull_datasets_1000.mat', 'X_inverse');
