function [x,Px,K,s,Ps,Psx,Kth] = xdkf1nstep(x0,th0,P0x,s0,P0s,P0sx,t,u,y,R,Q,gamma)

% XDKF-Z N
% - with assumption of zero gain sensitivity 
% - normalized objectives

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
if trace(P0s) > 0
% weight normalization
g_max = 0.99;
K_g0 = (A*P0x*C')/(C*P0x*C' + R);
[~,~,K_gmax] = xdkf1step(x0,th0,P0x,s0,P0s,P0sx,t,u,y,R,Q,g_max);
JxU = optimizationCriteria(K_g0,0,A,Ath_col,C,x0,th0,P0x,s0,P0s,P0sx,R,Q);
JxN = optimizationCriteria(K_gmax,0,A,Ath_col,C,x0,th0,P0x,s0,P0s,P0sx,R,Q);
JsU = optimizationCriteria(K_gmax,g_max,A,Ath_col,C,x0,th0,P0x,s0,P0s,P0sx,R,Q);
JsN = optimizationCriteria(K_g0,g_max,A,Ath_col,C,x0,th0,P0x,s0,P0s,P0sx,R,Q);
NormAlpha = 1/(JxN-JxU);
NormGamma = 1/(JsN-JsU);
alpha = alpha*NormAlpha+(1-g_max)*NormGamma;
gamma = gamma*g_max*NormGamma;
end

% Optimal gain
Psum = alpha*P0x;
Knum = zeros(nx,ny);
for p = 1:nth
    Athp = reshape(Ath_col(:,p),nx,nx);
    Psum = Psum + gamma(p)*P0s(:,:,p);
    Knum = Knum + gamma(p)*(Athp*x0)*s0(:,p)'*C';
end
K = (A*Psum*C' - Knum )/(C*Psum*C' + alpha*R);

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

function J = optimizationCriteria(K,gamma,A,Ath_col,C,x0,th0,P0x,s0,P0s,P0sx,R,Q)

% System dimensions
nth = size(th0,1);
nx = size(x0,1);

Px = (A-K*C)*P0x*(A-K*C)' + Q + K*R*K';

Popt = (1-gamma)*Px;

Ps = P0s;
Psx = P0sx;
for p = 1:nth
    Athp = reshape(Ath_col(:,p),nx,nx);
    P0sp  = P0s(:,:,p);    
    P0sxp = P0sx(:,:,p);
    Psx(:,:,p) = (A-K*C)*P0sxp*(A-K*C)';
    Ps(:,:,p) = (A-K*C)*P0sp*(A-K*C)' + (Athp*x0)*(Athp*x0)'...
        - (A-K*C)*s0(:,p)*(Athp*x0)' - (Athp*x0)*s0(:,p)'*(A-K*C)';
    Popt = Popt + gamma(p)*Ps(:,:,p);
end

J = trace(Popt);

end
