%
global beta rho n AlphaVec u_D_inv u_D C M


u_D_inv = @(u) u.^(-1/rho); % the inverse of u'
u_D     = @(c) c.^(-rho);   % u' (first derivative of u)

T = 99;
C = (0:(n-1))';
M = (0:(n-1))';

for t = 1:T
    ChiVec = u_D_inv(fun_GothVP(AlphaVec));
    MuVec = AlphaVec + ChiVec;
    M = [M, MuVec'];
    C = [C, ChiVec'];
end







