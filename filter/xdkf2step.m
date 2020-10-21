function [x,Px,K,s,Ps,Psx,Kth] = xdkf2step(x0,th0,P0x,s0,P0s,P0sx,t,u,y,R,Q,gamma)

% XDKF 
% - with assumption of independent gain sensitivity

% System dimensions
nth = size(th0,1);
nx = size(x0,1);
ny = size(y,1);
nu = size(u,1);
if nu<1
    u = 0;
end

% Estimate of zero noise mean 
w = zeros(nx,1);
e = zeros(ny,1);

% Reshape
x0 = reshape(x0,nx,1);
th0 = reshape(th0,nth,1);
y = reshape(y,ny,1);

% Evaluate system
A  = full(genmod('dfddx',t,x0,u,th0,w));
Ath_col = full(genmod('ddfddxdth',t,x0,u,th0,w));
C = full(genmod('dgdx',t,x0,u,th0,e));

% weight
alpha = 1-sum(gamma);

% Optimal gain
Knum = alpha*A*P0x*C';
Kden = alpha*(C*P0x*C'+R);
for p = 1:nth
    P0sp  = P0s(:,:,p);
    P0sxp = P0sx(:,:,p);
    P0xsp = P0sxp';
    Knum = Knum + gamma(p)*( A*P0sp*C' - (A*P0sxp*C')/(C*P0x*C' + R)*(C*P0xsp*C' ) );
    Kden = Kden + gamma(p)*( C*P0sp*C' - (C*P0sxp*C' )/(C*P0x*C' + R)*(C*P0xsp*C' ) );
end
K = Knum/Kden;

% Update state and error covariance
x = full(genmod('fd',t,x0,u,th0,w)) + K*(y-C*x0);
Px = (A-K*C)*P0x*(A-K*C)' + Q + K*R*K';

% Update error sensitivity covariance
Kth = zeros(nx,ny,nth);
Ps = P0s;
Psx = P0sx;
for p = 1:nth
    Athp = reshape(Ath_col(:,p),nx,nx);
    P0sp  = P0s(:,:,p);
    P0sxp = P0sx(:,:,p);
    P0xsp = P0sxp';
    Kth(:,:,p) = ( (A-K*C)*P0sxp*C')/(C*P0x*C' + R);
    Psx(:,:,p) = (-Kth(:,:,p)*C)*P0x*(A-K*C)' + (A-K*C)*P0sxp*(A-K*C)' + Kth(:,:,p)*R*K';
    Ps(:,:,p) = (-Kth(:,:,p)*C)*P0x*(-Kth(:,:,p)*C)' + (A-K*C)*P0sp*(A-K*C)' ...
        + (-Kth(:,:,p)*C)*P0xsp*(A-K*C)' + (A-K*C)*P0sxp*(-Kth(:,:,p)*C)' ...
        + Kth(:,:,p)*R*Kth(:,:,p)'...
        + (Athp*x0)*(Athp*x0)';
end

% Sensitivity update -- for analysis purpose only
s = s0;
for p = 1:nth
    Athp = reshape(Ath_col(:,p),nx,nx);
    s(:,p) = (A-K*C)*s0(:,p) - Kth(:,:,p)*(y-C*x0) - Athp*x0;
end

end