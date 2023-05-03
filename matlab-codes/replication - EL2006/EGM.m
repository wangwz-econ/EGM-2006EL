clear all


beta = 0.96;   % discount factor
rho = 2;       % RRA of the utility function
tau = 0.90;    % 1 - depreciation rate
alpha = 0.36;  % Cobb-Douglous production function

g = 1.01;                  % derterministic part of labor growth
phi = [0.9, 1, 1.1];       % stochastic part of labor growth
P_phi = [0.25, 0.5, 0.25]; % corresponding probability of different states of labor growth
Delta = g * phi;           % gross labor growth (row vector)
P_Delta = P_phi;           % corresponding prob of each gross labor growth state (row vector)


NGrid = 20;    % number of grid points
kSS = ( ( ((g^rho)/(beta*tau) )-1) / alpha )^( 1/(alpha-1) );
aSS = kSS*g/tau;
savings = transpose(exp(linspace(0, log(3*aSS), NGrid))-1); % set up grid points (column vector)

T = 99;       % number of periods
consump = NaN(NGrid, T);
wealth  = NaN(NGrid, T);

u_D_inv = @(u) u.^(-1/rho); % the inverse of u'
u_D     = @(c) c.^(-rho);   % u' (first derivative of u)

capital_tplus1 = tau * savings * (1 ./ Delta);
wealth_tplus1 = capital_tplus1 + (capital_tplus1 .^ alpha);
return_capital = 1 + alpha * (capital_tplus1.^(alpha-1)); 

for t = T:-1:1
  if t == T
    consump(:, t) = savings;
    wealth(:, t) = savings;
  else
    consump_tplus1 = interp1(wealth(:, t+1), consump(:, t+1), wealth_tplus1, 'linear','extrap');
    disp(consump_tplus1)
    RHS = tau*beta* P_Delta * ...
      transpose( (repmat(Delta.^(-rho),NGrid,1) .* u_D(consump_tplus1) .* return_capital ) );
    consump(:, t) = transpose(u_D_inv(RHS));
    wealth(:, t)  = consump(:, t) + savings;
  end
end

for t = (T):-1:1
	plot(wealth(:, t), consump(:, t))
	hold all
end
xlabel('Wealth')
ylabel('Optimal consumption')
title('My EGM Algorithm')

saveas(gcf,'Policy Function (c wrt w) in Different Periods.jpg')

