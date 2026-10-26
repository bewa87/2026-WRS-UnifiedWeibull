%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 02: Bisection          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

% Loading of Data Sets For Comparison

%load('example02_weibull_datasets_1000.mat');
load('example02_inverse_weibull_datasets_1000.mat');

N_sets         = 125;

if(N_sets > 1000)
  error('You choose more data sets than available.');
endif

z_mean_bis     = zeros(N_sets,1); %zeros(length(X_weibull),1);
time_mean_bis  = zeros(N_sets,1); %zeros(length(X_weibull),1);
time_med_bis   = zeros(N_sets,1); %zeros(length(X_weibull),1);
it_mean_bis    = zeros(N_sets,1); %zeros(length(X_weibull),1);
it_med_bis     = zeros(N_sets,1); %zeros(length(X_weibull),1);
width_mean_bis = zeros(N_sets,1); %zeros(length(X_weibull),1);

% Choice of Probability Density Function

s= - 1; %s = 1;

% Loop over All Randomly Distributed Data Sets

for p = 1:1:125 % p = 1:1:length(X_weibull)

  X       = X_inverse(:,p); %X_weibull(:,p);
  Y       = log(X);
  n       = numel(X);
  Ymean   = mean(Y);
  Yrange  = max(Y) - min(Y);

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

  if (~(isfinite(a_init) && isfinite(b_init) && ...
      a_init > 0 && b_init > a_init))
    error('Unvalid initial interval.');
  endif

  % Maximum Number of Iterations and Tolerance

  N_max = 50;
  tol   = 1e-12;

  % Check if the Initial Interval Includes the Unique Zero

  if(G_s(a_init) < 0 || G_s(b_init) > 0)
    error('No change of sign in this initial interval.');
  endif

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Application of Bisection Method %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  N_reps    = 1000;
  z_vec     = zeros(N_reps,1);
  width_vec = zeros(N_reps,1);
  time_vec  = zeros(N_reps,1);
  it_vec    = zeros(N_reps,1);

  for j = 1:1:N_reps

    a = a_init;
    b = b_init;

    ga = G_s(a);
    gb = G_s(b);

    tic;
    for k = 1:1:N_max
      x = a + (b-a)/2;
      gx = G_s(x);

      if(gx == 0)
        a = x;
        b = x;
        break;
      elseif(gx > 0)
        a = x;
        ga = gx;
      else
        b = x;
        gb = gx;
      endif

      if(b-a <= tol)
        break;
      endif
    endfor
    time_bis = toc;

    z_vec(j)     = a + (b-a)/2;
    width_vec(j) = b-a;
    it_vec(j)    = k;
    time_vec(j)  = time_bis;

  endfor

z_mean_bis(p)     = mean(z_vec);
time_mean_bis(p)  = mean(time_vec);
time_med_bis(p)   = median(time_vec);
it_mean_bis(p)    = mean(it_vec);
it_med_bis(p)     = median(it_vec);
width_mean_bis(p) = mean(width_vec);

endfor

z_mean_final     = mean(z_mean_bis);
time_mean_final  = mean(time_mean_bis);
time_med_final   = median(time_med_bis);
it_mean_final    = mean(it_mean_bis);
it_med_final     = median(it_med_bis);
width_mean_final = mean(width_mean_bis);
