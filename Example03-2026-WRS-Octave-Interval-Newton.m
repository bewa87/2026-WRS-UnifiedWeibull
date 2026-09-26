% Step 1: Load Interval-Package For GNU Octave 8.4.0

pkg load interval

% Step 2: Preparation of Velocity Vector

pkg load statistics

fid = fopen('ReadData-01639.txt','r');
C   = textscan(fid,'%f32 %f32 %f32 %f32 %f32 %s','delimiter',',');
fclose(fid);

vel_proc       = cell2mat(C(4));
final_vel_proc = sort(vel_proc);
ind_proc       = find(final_vel_proc>=0.09);
final_vel_proc = final_vel_proc(ind_proc);

% Step 3: Define Corresponding Data Set

X    =     final_vel_proc;

% Step 4: Choice Of Probability Distribution [s  = 1 (Weibull), s = -1 (inverse Weibull)]

s    = 1;

% Step 5: Define Functions G_s, dG_s And H_s

G_s  = @(z) (1./z) + (s./length(X)).*sum(log(X)) - s.*sum(((X.^(s.*z))./(sum(X.^(s.*z)))).*log(X));
dG_s = @(z) - (1./(z.^2)) - (1./(








































































































(sum(X.^(s.*z))).^2)).*((sum(log(X).*log(X).*(X.^(s.*z)))).*(sum(X.^(s.*z)))-((sum(log(X).*(X.^(s.*z)))).^2));
H_s  = @(z) sum(((X.^(s.*z))./(sum(X.^(s.*z)))).*log(X));

% Step 6: Define Initial Interval

a_init = 1/(log(max(X)) - (1./length(X)).*sum(log(X)));
b_init = min([1/(s*(H_s(a_init) - (1./length(X)).*sum(log(X)))) (length(X).*length(X)+sqrt((length(X)).^4+4.*(length(X)).^2))./(2*(max(log(X))-min(log(X))))]);

% Step 7: Interval Newton Method (GNU Octave Interval-Package)

sol = fzero(G_s, infsup("[0.46680114, 3.6508296]"), dG_s)
