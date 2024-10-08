function [x,Px,K,s,Ps,Psx,Kth] = xdkfstep(x0,th0,P0x,s0,P0s,P0sx,t,u,y,R,Q,gamma)

% XDKF
% - optimal

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
% if trace(P0s) > 0
% % weight normalization
% K_g0 = (A*P0x*C' )/(C*P0x*C' + R);
% g_max = 0.99;
% a_max = 1-g_max;
% Psum = a_max*P0x;
% Ssum = zeros(nx,ny);
% Kth_gmax = zeros(nx,ny,nth);
% for p = 1:nth
%     Psum = Psum + g_max(p)*P0s(:,:,p);
% end
% for p = 1:nth
%     Athp = reshape(Ath_col(:,p),nx,nx);
%     Kth_gmax(:,:,p) = Athp*Psum*C'/(C*Psum*C'+alpha*R);
%     Ssum = Ssum - g_max(p)*(Kth_gmax(:,:,p)*C*P0sx(:,:,p)'*C' + (Athp*x0)*s0(:,p)'*C');
% end
% K_gmax = (A*Psum*C' + Ssum )/(C*Psum*C' + a_max*R);
% JxU = optimizationCriteria(K_g0,0,0,A,Ath_col,C,x0,th0,P0x,s0,P0s,P0sx,R,Q);
% JxN = optimizationCriteria(K_gmax,Kth_gmax,0,A,Ath_col,C,x0,th0,P0x,s0,P0s,P0sx,R,Q);
% JsU = optimizationCriteria(K_gmax,Kth_gmax,g_max,A,Ath_col,C,x0,th0,P0x,s0,P0s,P0sx,R,Q);
% JsN = optimizationCriteria(K_g0,0,g_max,A,Ath_col,C,x0,th0,P0x,s0,P0s,P0sx,R,Q);
% NormAlpha = 1/(JxN-JxU);
% NormGamma = 1/(JsN-JsU);
% alpha = (NormAlpha+a_max*NormGamma)*alpha;
% gamma = (g_max*NormGamma)*gamma;
% end

% Optimal gain
Kth = zeros(nx,ny,nth);
Ssum = zeros(nx,ny);
Psum = alpha*P0x;
for p = 1:nth
    Psum = Psum + gamma(p)*P0s(:,:,p);
end
for p = 1:nth
    Athp = reshape(Ath_col(:,p),nx,nx);
    Kth(:,:,p) = Athp*Psum*C'/(C*Psum*C'+alpha*R);
    Ssum = Ssum - gamma(p)*(Kth(:,:,p)*C*P0sx(:,:,p)'*C' + (Athp*x0)*s0(:,p)'*C');
end
K = (A*Psum*C' + Ssum)/(C*Psum*C' + alpha*R);

% Update state and error covariance
x = full(genmod('fd',t,x0,u,th0,w)) + K*(y-C*x0);
Px = (A-K*C)*P0x*(A-K*C)' + Q + K*R*K';

% Update error sensitivity covariance
Ps = P0s;
Psx = P0sx;
for p = 1:nth
    Athp = reshape(Ath_col(:,p),nx,nx);
    P0sp  = P0s(:,:,p);
    P0sxp = P0sx(:,:,p);
    P0xsp = P0sxp';
    Psx(:,:,1) = (-Kth(:,:,p)*C)*P0x*(A-K*C)' + (A-K*C)*P0sxp*(A-K*C)' + Kth*R*K';
    Ps(:,:,1) = (-Kth(:,:,p)*C)*P0x*(-Kth(:,:,p)*C)' + (-Kth*C)*P0xsp*(A-K*C)' ...
        + (A-K*C)*P0sxp*(-Kth(:,:,p)*C)' + (A-K*C)*P0sp*(A-K*C)' ...
        + Kth(:,:,p)*R*Kth(:,:,p)' ...
        + Athp*x0*(Athp*x0)' ...
        - (A-K*C)*s0(:,p)*(Athp*x0)' - (Athp*x0)*s0(:,p)'*(A-K*C)';
end

% Sensitivity update
s = s0;
for p = 1:nth
    Athp = reshape(Ath_col(:,p),nx,nx);
    s(:,p) = (A-K*C)*s0(:,p) - Kth(:,:,p)*(y-C*x0) - Athp*x0;
end

end
