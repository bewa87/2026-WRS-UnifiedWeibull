% Define the corresponding data set

X    =     [0.35 0.59 0.96 0.99 1.69...
            1.97 2.07 2.58 2.71 2.90...
            3.67 3.99 5.35 13.77 25.50]';

% s  = 1 (Weibull), s = -1 (inverse Weibull)

s    = -1;

% Define functions G_s and H_s

G_s  = @(z) (1./z) + (s./length(X)).*sum(log(X)) - s.*sum(((X.^(s.*z))./(sum(X.^(s.*z)))).*log(X));
H_s  = @(z) sum(((X.^(s.*z))./(sum(X.^(s.*z)))).*log(X));

% Define initial Interval

a_init = 1/(log(max(X)) - (1./length(X)).*sum(log(X)));
b_init = min([1/(s*(H_s(a_init) - (1./length(X)).*sum(log(X)))) (length(X).*length(X)+sqrt((length(X)).^4+4.*(length(X)).^2))./(2*(max(log(X))-min(log(X))))]);

% Start the loop

N_max     = 50;
tolerance = 10^(-16);

A      = zeros(1,(N_max+1));
B      = zeros(1,(N_max+1));

A(1)   = a_init;
B(1)   = b_init;

tic

for j = 1:1:N_max
  m = (A(j) + B(j))/2;
  if(G_s(m)==0)
    A(j + 1) = m;
    B(j + 1) = m;
    break
  elseif(G_s(A(j))*G_s(m)<0)
    A(j + 1) = A(j);
    B(j + 1) = m;
  elseif(G_s(m)*G_s(B(j))<0)
    A(j + 1) = m;
    B(j + 1) = B(j);
  endif

  if(B(j+1) - A(j+1) < tolerance)
    break
  endif

endfor

toc

z_sought = (B(j+1) + A(j+1))/2;

figure(1)
plot(a_init:0.01:b_init,G_s(a_init:0.01:b_init),'linewidth',1.0)
hold on
plot(a_init:0.01:b_init,0*ones(length(a_init:0.01:b_init),1),'linewidth',1.0,'color','red','linestyle','--')
hold off
title('Plot of G_s','fontsize',18)
xlabel('Inverse Weibull parameter \beta','fontsize',14,'interpreter','tex')
ylabel('G_s(\beta)','fontsize',14,'interpreter','tex')

figure(2)
semilogy(0:1:j,(B(1:1:j+1)-A(1:1:j+1)),'linewidth',1.0)
title('Semi-Logarithmic Error Plot (Bisection)','fontsize',18)
yticks([10^0 10^(-4) 10^(-8) 10^(-12) 10^(-16)])
xlabel('Iteration l','fontsize',14,'interpreter','tex')
ylabel('Error |b_{l} - a_{l}|','fontsize',14,'interpreter','tex')
