function [x,Px,K,s,Ps,Psx,Kth] = xdkf1step(x0,th0,P0x,s0,P0s,P0sx,t,u,y,R,Q,gamma)

% XDKF 
% - with assumption of zero gain sensitivity 

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
f = full(genmod('fd',t,x0,u,th0,w));
g = full(genmod('g',t,x0,u,th0,e));

% weight
alpha = 1-sum(gamma);

% Optimal gain
Psum = alpha*P0x;
Ssum = zeros(nx,ny);
for p = 1:nth
    Athp = reshape(Ath_col(:,p),nx,nx);
    Psum = Psum + gamma(p)*P0s(:,:,p);
    Ssum = Ssum + gamma(p)*(Athp*x0)*s0(:,p)'*C';
end
K = (A*Psum*C' - Ssum )/(C*Psum*C' + alpha*R);

% Update state and error covariance
x = f + K*(y-g);
Px = (A-K*C)*P0x*(A-K*C)' + Q + K*R*K';

% Update error sensitivity covariance
Kth = zeros(nx,ny,nth);
Ps = P0s;
Psx = P0sx;
for p = 1:nth
    Athp = reshape(Ath_col(:,p),nx,nx);
    P0sp  = P0s(:,:,p);    
    P0sxp = P0sx(:,:,p);
    Psx(:,:,p) = (A-K*C)*P0sxp*(A-K*C)';
    Ps(:,:,p) = (A-K*C)*P0sp*(A-K*C)' + (Athp*x0)*(Athp*x0)'...
        - (A-K*C)*s0(:,p)*(Athp*x0)' - (Athp*x0)*s0(:,p)'*(A-K*C)';
end

% Sensitivity update
s = s0;
for p = 1:nth
    Athp = reshape(Ath_col(:,p),nx,nx);
    s(:,p) = (A-K*C)*s0(:,p) - Athp*x0;
end

end