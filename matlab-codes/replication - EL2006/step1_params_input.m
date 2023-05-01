% Define relevant parameters

clear all

global beta rho n g dep alpha permval permprob R W AlphaVec

beta = 0.96;
rho = 2;
n = 20;
g = 1.01;
dep = 0.90;
alpha = 0.36;

permval = [0.9, 1, 1.1];
permprob = [0.25, 0.5, 0.25];

R = @(k) 1 + ( alpha .*k .^(alpha-1) );
W = @(k) (1-alpha).*(k.^alpha);

kSS = ( ( ((g^rho)/(beta*dep) )-1) / alpha )^( 1/(alpha-1) );
aSS = kSS*g/dep;
AlphaVec = exp( linspace(0.00, log(3*aSS), n) ) - 1; % set up grid points










