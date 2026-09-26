% Step 1: Load Interval-Package For GNU Octave 8.4.0

pkg load interval

% Step 2: Define Corresponding Data Set

X    =     [0.35 0.59 0.96 0.99 1.69...
            1.97 2.07 2.58 2.71 2.90...
            3.67 3.99 5.35 13.77 25.50]';

% Step 3: Choice Of Probability Distribution [s  = 1 (Weibull), s = -1 (inverse Weibull)]

s    = -1;

% Step 4: Define Functions G_s, dG_s And H_s

G_s  = @(z) (1./z) + (s./length(X)).*sum(log(X)) - s.*sum(((X.^(s.*z))./(sum(X.^(s.*z)))).*log(X));
dG_s = @(z) - (1./(z.^2)) - (1./((sum(X.^(s.*z))).^2)).*((sum(log(X).*log(X).*(X.^(s.*z)))).*(sum(X.^(s.*z)))-((sum(log(X).*(X.^(s.*z)))).^2));
H_s  = @(z) sum(((X.^(s.*z))./(sum(X.^(s.*z)))).*log(X));

% Step 5: Define Initial Interval (Guaranteed Convergence Region)

a_init = 1/(log(max(X)) - (1./length(X)).*sum(log(X)));
b_init = min([1/(s*(H_s(a_init) - (1./length(X)).*sum(log(X)))) (length(X).*length(X)+sqrt((length(X)).^4+4.*(length(X)).^2))./(2*(max(log(X))-min(log(X))))]);

% Step 6: Interval Newton Method (GNU Octave Interval-Package)

sol = fzero(G_s, infsup("[0.4280028657631536, 2.203152719856162]"), dG_s)
