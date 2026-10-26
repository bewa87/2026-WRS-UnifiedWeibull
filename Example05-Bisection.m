%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step I: Initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clearing and Closing

clear all;
close all;
clc;

% Data Vector X and Choice of Probability Density Function

X = [3.071 3.140 3.147 3.256 3.759...
     3.955 4.025 4.158 4.629 4.912...
     5.170 5.205 6.321 6.735 8.250...
     8.777 9.773 10.218 11.909 12.446]';
s = - 1; % s = 1: Weibull, s = -1: inverse Weibull

% Logarithmic Data Vector Y, mean and the range of Y

Y = log(X);
n = numel(X);
Ymean = mean(Y);
Yrange = max(Y) - min(Y);

% Score Functions: G_s and dG_s are evaluated

evalGH = @(z) local_evalGH(z, Y, s, n, Ymean);

% Separate Function, if only G or dG_s are needed

G_s  = @(z) evalGH(z)(:,1);
dG_s = @(z) evalGH(z)(:,2);

% Initial Bounds for the Initial Interval

if(s == -1)
  a_init = 1/(Ymean - min(Y));
endif

if(s == 1)
  a_init = 1/(max(Y) - Ymean);
endif

w = exp(s*a_init*Y);
H_a = sum((w/sum(w)).*Y);

if(s == -1)
  b_init = min([ 1/(s*(H_a - Ymean)), ...
                 (n^2 + sqrt(n^4 + 4*n^2)) / (2*Yrange)]);
endif

if(s == 1)
  b_init = min([ 1/(s*(H_a - Ymean)), ...
                 (n^2 + sqrt(n^4 + 4*n^2)) / (2*Yrange) ]);
endif

if(~(isfinite(a_init) && isfinite(b_init) && ...
     a_init > 0 && b_init > a_init))
  error('Unvalid initial interval.');
endif

% Maximum Number of Iterations and Tolerance

N_max = 50;
tol   = 1e-12;

% Mao/Li-Schranken
%MaoM  = @(a) (1./a.^2) + Yrange^2;
%Mao_m = @(b) (1./b.^2);

% Wacker-Schranken
%WackerM  = @(a) (1./a.^2) + (Yrange^2/4);
%Wacker_m = @(b) (1./b.^2);

% Check if the Initial Interval Includes the Unique Zero

if(G_s(a_init) < 0 || G_s(b_init) > 0)
    error('No change of sign in this initial interval.');
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step II: Bisection Method %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_vec    = zeros(N_max+1,1);
B_vec    = zeros(N_max+1,1);

a = a_init;
b = b_init;

A_vec(1) = a_init;
B_vec(1) = b_init;

ga = G_s(a);
gb = G_s(b);

tic;
for k = 1:1:N_max
  x = a + (b-a)/2;
  gx = G_s(x);

  if(gx == 0)
    a = x;
    b = x;

    A_vec(k+1) = a;
    B_vec(k+1) = b;

    break;
  elseif(gx > 0)
    a = x;
    ga = gx;

    A_vec(k+1) = x;
    B_vec(k+1) = b;

  else
    b = x;
    gb = gx;

    A_vec(k+1) = a;
    B_vec(k+1) = b;

  endif

  if(b-a <= tol)
    break;
  endif
endfor
time_bis = toc;

z_bis_star = a + (b-a)/2;
width_bis  = b-a;

fprintf('Bisection: z = %.16g, Width = %.3e, Time = %.6g s\n', ...
        z_bis_star, width_bis, time_bis);

% Plot of G_{s,Y}

figure(1)
plot(1:0.01:4,G_s(1:0.01:4),'linewidth',1.0)
hold on
plot(1:0.01:4,0*ones(length(1:0.01:4),1),'linewidth',1.0,'color','red','linestyle','--')
hold off
title('Plot of G_{s,Y}','fontsize',18,'interpreter','tex')
xlabel('Inverse Weibull parameter \beta','fontsize',14,'interpreter','tex')
ylabel('G_{s,Y}(\beta)','fontsize',14,'interpreter','tex')

% Semi-Logarithmic Error Plot

figure(2)
semilogy(0:1:41,(B_vec(1:1:42)-A_vec(1:1:42)),'linewidth',1.0)
title('Semi-Logarithmic Error Plot (Bisection)','fontsize',18)
yticks([10^0 10^(-4) 10^(-8) 10^(-12) 10^(-16)])
xlabel('Iteration l','fontsize',14,'interpreter','tex')
ylabel('Error |b_{l} - a_{l}|','fontsize',14,'interpreter','tex')
