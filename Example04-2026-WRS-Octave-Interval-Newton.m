% Step 1: Load Interval-Package For GNU Octave 8.4.0

pkg load interval

% Step 2: Define Corresponding Data Set (Jana 2022, Drug 3 from Table 6)

X    =     [3.071 3.140 3.147 3.256 3.759...
            3.955 4.025 4.158 4.629 4.912...
            5.170 5.205 6.321 6.735 8.250...
            8.777 9.773 10.218 11.909 12.446]';

% Step 3: Choice Of Probability Distribution [s  = 1 (Weibull), s = -1 (inverse Weibull)]

s    = -1;

% Step 4: Define Functions G_s, dG_s And H_s

G_s  = @(z) (1./z) + (s./length(X)).*sum(log(X)) - s.*sum(((X.^(s.*z))./(sum(X.^(s.*z)))).*log(X));
dG_s = @(z) - (1./(z.^2)) - (1./((sum(X.^(s.*z))).^2)).*((sum(log(X).*log(X).*(X.^(s.*z)))).*(sum(X.^(s.*z)))-((sum(log(X).*(X.^(s.*z)))).^2));
H_s  = @(z) sum(((X.^(s.*z))./(sum(X.^(s.*z)))).*log(X));

% Step 5: Define Initial Interval

a_init = 1/(log(max(X)) - (1./length(X)).*sum(log(X)));
b_init = min([1/(s*(H_s(a_init) - (1./length(X)).*sum(log(X)))) (length(X).*length(X)+sqrt((length(X)).^4+4.*(length(X)).^2))./(2*(max(log(X))-min(log(X))))]);

% Step 6: Interval Newton Method (GNU Octave Interval-Package)

sol = fzero(G_s, infsup("[1.229589403549322, 4.622760482352577]"), dG_s)
