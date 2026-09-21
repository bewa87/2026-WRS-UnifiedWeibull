% Define the corresponding data set (Jana 2022, Drug 3 from Table 6)

X    =     [3.071 3.140 3.147 3.256 3.759...
            3.955 4.025 4.158 4.629 4.912...
            5.170 5.205 6.321 6.735 8.250...
            8.777 9.773 10.218 11.909 12.446]';

% s  = 1 (Weibull), s = -1 (inverse Weibull)

s    = -1;

% Define functions G_s, dG_s and H_s

G_s  = @(z) (1./z) + (s./length(X)).*sum(log(X)) - s.*sum(((X.^(s.*z))./(sum(X.^(s.*z)))).*log(X));
dG_s = @(z) - (1./(z.^2)) - (1./((sum(X.^(s.*z))).^2)).*((sum(log(X).*log(X).*(X.^(s.*z)))).*(sum(X.^(s.*z)))-((sum(log(X).*(X.^(s.*z)))).^2));
H_s  = @(z) sum(((X.^(s.*z))./(sum(X.^(s.*z)))).*log(X));

% Define initial Interval

a_init = 1/(log(max(X)) - (1./length(X)).*sum(log(X)));
b_init = min([1/(s*(H_s(a_init) - (1./length(X)).*sum(log(X)))) (length(X).*length(X)+sqrt((length(X)).^4+4.*(length(X)).^2))./(2*(max(log(X))-min(log(X))))]);

% Start the loop

N_max     = 20;
tolerance = 10^(-16);
N_points  = 100;

A      = zeros(1,(N_max+1));
B      = zeros(1,(N_max+1));

A(1)   = a_init;
B(1)   = b_init;

if(s == 1)
tic
  for j = 1:1:N_max
    m = (A(j) + B(j))/2;
    if(G_s(m) == 0)
      A(j + 1) = m;
      B(j + 1) = m;
      break
    elseif(G_s(m) < 0)
      A(j + 1) = max([A(j) (m - G_s(m)./(max(dG_s(A(j):(B(j)-A(j))/N_points:B(j)))) )]);
      B(j + 1) = min([B(j) (m - G_s(m)./(min(dG_s(A(j):(B(j)-A(j))/N_points:B(j)))) )]);
    elseif(G_s(m) > 0)
      A(j + 1) = max([A(j) (m - G_s(m)./(min(dG_s(A(j):(B(j)-A(j))/N_points:B(j)))) )]);
      B(j + 1) = min([B(j) (m - G_s(m)./(max(dG_s(A(j):(B(j)-A(j))/N_points:B(j)))) )]);
    endif

    if(B(j+1) - A(j+1) < tolerance)
      break
    endif
  endfor
toc
endif

if(s == -1)
tic
  for j = 1:1:N_max
    m = (A(j) + B(j))/2;
    if(G_s(m) == 0)
      A(j + 1) = m;
      B(j + 1) = m;
      break
    elseif(G_s(m) < 0)
      A(j + 1) = max([A(j) (m - G_s(m)./(max(dG_s(A(j):(B(j)-A(j))/N_points:B(j)))) )]);
      B(j + 1) = min([B(j) (m - G_s(m)./(min(dG_s(A(j):(B(j)-A(j))/N_points:B(j)))) )]);
    elseif(G_s(m) > 0)
      A(j + 1) = max([A(j) (m - G_s(m)./(min(dG_s(A(j):(B(j)-A(j))/N_points:B(j)))) )]);
      B(j + 1) = min([B(j) (m - G_s(m)./(max(dG_s(A(j):(B(j)-A(j))/N_points:B(j)))) )]);
    endif

    if(B(j+1) - A(j+1) < tolerance)
      break
    endif

  endfor
toc
endif

z_sought = (B(j+1) + A(j+1))/2;

figure(3)
semilogy(0:1:j,(B(1:1:j+1)-A(1:1:j+1)),'linewidth',1.0)
title('Semi-Logarithmic Error Plot (Interval Newton)','fontsize',18)
yticks([10^0 10^(-4) 10^(-8) 10^(-12) 10^(-16)])
xlabel('Iteration l','fontsize',14,'interpreter','tex')
ylabel('Error |b_{l} - a_{l}|','fontsize',14,'interpreter','tex')
