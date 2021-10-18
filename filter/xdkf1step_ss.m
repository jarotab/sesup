function [x,Px,K,s,Ps,Psx,Kth] = xdkf1step_ss(x0,th0,P0x,s0,P0s,P0sx,t,u,y,R,Q,gamma)

% XDKF 
% - with assumption of zero gain sensitivity 
% - steady-state gain

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
% DARE config
drR = alpha*R;
drQ = alpha*Q;
drS = zeros(nx,ny);
for p = 1:nth
    Athp = reshape(Ath_col(:,p),nx,nx);
    drS = drS - gamma(p)*(Athp*x0)*s0(:,p)'*C';
    drQ = drQ + gamma(p)*((Athp*x0)*(Athp*x0)'...
        - (Athp*x0)*s0(:,p)'*A' - A*s0(:,p)*(Athp*x0)');
end
drB = C';
drA = A';
[~,drK,~] = idare(drA,drB,drQ,drR,drS,[]);
K = drK';
            
% Update state and error covariance
x = full(genmod('fd',t,x0,u,th0,w)) + K*(y-C*x0);
Px = (A-K*C)*P0x*(A-K*C)' + Q + K*R*K';

% Update error sensitivity covariance
Kth = zeros(nx,ny,nth);
Ps = P0s;
Psx = P0sx;
for p = 1:nth
    P0sp  = P0s(:,:,p);
    P0sxp = P0sx(:,:,p);
    Psx(:,:,p) = (A-K*C)*P0sxp*(A-K*C)';
    Ps(:,:,p) = (A-K*C)*P0sp*(A-K*C)' + (Athp*x0)*(Athp*x0)';
end

% Sensitivity update -- for analysis purpose only
s = s0;
for p = 1:nth
    Athp = reshape(Ath_col(:,p),nx,nx);
    s(:,p) = (A-K*C)*s0(:,p) - Athp*x0;
end

end