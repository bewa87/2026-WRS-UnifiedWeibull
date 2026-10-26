%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step I: Initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clearing and Closing

clear all;
close all;
clc;

% Data Vector X and Choice of Probability Density Function

X = [0.35 0.59 0.96 0.99 1.69...
     1.97 2.07 2.58 2.71 2.90...
     3.67 3.99 5.35 13.77 25.50]';
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step II: Classical Interval-Newton with Mao/Li Derivative Bounds %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_vec_mao    = zeros(N_max+1,1);
B_vec_mao    = zeros(N_max+1,1);

A_vec_wac    = zeros(N_max+1,1);
B_vec_wac    = zeros(N_max+1,1);

a = a_init;
b = b_init;

A_vec_mao(1) = a_init;
B_vec_mao(1) = b_init;

A_vec_wac(1) = a_init;
B_vec_wac(1) = b_init;

% Loop for Mao/Li-Bounds in Interval-Newton

for k = 1:1:N_max
  x  = a + (b-a)/2;
  gx = G_s(x);

  if(gx == 0)
    a = x;
    b = x;

    A_vec_mao(k+1) = a;
    B_vec_mao(k+1) = b;

    break;
  endif

  % Mao/Li-Bounds for -G_s'(z)
  M = (1/a^2) + Yrange^2;
  m = (1/b^2);

  if(~(isfinite(M) && isfinite(m) && M >= m && m > 0))
    error('Unvalid Mao/Li-Bounds.');
  endif

  % Interval-Newton
  t1 = x + gx/M;
  t2 = x + gx/m;

  Nlow  = min(t1, t2);
  Nhigh = max(t1, t2);

  % Intersection with current [a, b]
  anew = max(a, Nlow);
  bnew = min(b, Nhigh);

  A_vec_mao(k+1) = anew;
  B_vec_mao(k+1) = bnew;

  % If the Intersection is Empty, Apply Bisection
  if(anew > bnew)
    if(gx > 0)
      anew  = x;
      bnew  = b;

      A_vec_mao(k+1) = anew;
      B_vec_mao(k+1) = bnew;

    else
      anew  = a;
      bnew  = x;

      A_vec_mao(k+1) = anew;
      B_vec_mao(k+1) = bnew;

    endif
  endif

  a = anew;
  b = bnew;

  if(b-a <= tol)
    break;
  endif
endfor
time_newton_mao = toc;

z_newton_mao_star = a + (b-a)/2;
width_newton_mao  = b-a;

fprintf(['Interval-Newton Mao/Li: z = %.16g, ', ...
         'Width = %.3e, Time = %.6g s\n'], ...
        z_newton_mao_star, width_newton_mao, time_newton_mao);

% Loop for Wacker-Bounds in Interval-Newton

a = a_init;
b = b_init;

for k = 1:1:N_max
  x  = a + (b-a)/2;
  gx = G_s(x);

  if(gx == 0)
    a = x;
    b = x;

    A_vec_wac(k+1) = a;
    B_vec_wac(k+1) = b;

    break;
  endif

  % Wacker-Bounds for -G_s'(z)
  M = (1/a^2) + (Yrange^2/4);
  m = (1/b^2);

  if(~(isfinite(M) && isfinite(m) && M >= m && m > 0))
    error('Unvalid Wacker-Bounds.');
  endif

  % Interval-Newton
  t1 = x + gx/M;
  t2 = x + gx/m;

  Nlow  = min(t1, t2);
  Nhigh = max(t1, t2);

  % Intersection with current [a, b]
  anew = max(a, Nlow);
  bnew = min(b, Nhigh);

  A_vec_wac(k+1) = anew;
  B_vec_wac(k+1) = bnew;

  % If the Intersection is Empty, Apply Bisection
  if(anew > bnew)
    if(gx > 0)
      anew  = x;
      bnew  = b;

      A_vec_wac(k+1) = anew;
      B_vec_wac(k+1) = bnew;

    else
      anew  = a;
      bnew  = x;

      A_vec_wac(k+1) = anew;
      B_vec_wac(k+1) = bnew;

    endif
  endif

  a = anew;
  b = bnew;

  if(b-a <= tol)
    break;
  endif
endfor
time_newton_wac = toc;

z_newton_wac_star = a + (b-a)/2;
width_newton_wac  = b-a;

fprintf(['Interval-Newton Wacker: z = %.16g, ', ...
         'Width = %.3e, Time = %.6g s\n'], ...
        z_newton_wac_star, width_newton_wac, time_newton_wac);

% Semi-Logarithmic Error Plot

figure(3)
semilogy(0:1:12,(B_vec_mao(1:1:13)-A_vec_mao(1:1:13)),'linewidth',1.0,'color','blue','marker','+')
hold on
semilogy(0:1:12,(B_vec_wac(1:1:13)-A_vec_wac(1:1:13)),'linewidth',1.0,'color','red','marker','*')
hold on
title('Semi-Logarithmic Error Plot (Interval Newton)','fontsize',18)
yticks([10^0 10^(-4) 10^(-8) 10^(-12) 10^(-16)])
xlabel('Iteration l','fontsize',14,'interpreter','tex')
ylabel('Error |b_{l} - a_{l}|','fontsize',14,'interpreter','tex')
leg = legend({'Mao/Li','Wacker'})
legend(leg,'location','eastoutside')
set(leg,'fontsize',14)
hold off
