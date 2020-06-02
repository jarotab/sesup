function [x,Px,K,Ps,Psx,Kr] = xdkf1step(x0,th0,P0x,P0s,P0sx,t,u,y,R,Q,Rss,Rys,gamma)

% RES-KF-A1 algorithm

nth = size(th0,1);
nx = size(x0,1);
ny = size(y,1);
nu = size(u,1);
if nu<1
    u = 0;
end

w = zeros(nx,1);
e = zeros(ny,1);

x0 = reshape(x0,nx,1);
th0 = reshape(th0,nth,1);
y = reshape(y,ny,1);

A  = full(genmod('dfddx',t,x0,u,th0,w));
Ath_col = full(genmod('ddfddxdth',t,x0,u,th0,w));
C = full(genmod('dgdx',t,x0,u,th0,e));

% Additional noise 
V = Rss;
S = Rys;

alpha = 1-sum(gamma);
% Optimal gain
Knum = alpha*A*P0x*C';
Kden = alpha*(C*P0x*C'+R);
for p = 1:nth
    Athp = reshape(Ath_col(:,p),nx,nx);
    P0sp  = P0s(:,:,p);
    P0sxp = P0sx(:,:,p);
    P0xsp = P0sxp';
    Knum = Knum + gamma(p)*(Athp*P0xsp*C' + A*P0sp*C');
    Kden = Kden + gamma(p)*(C*P0sp*C'+ V);
end
K = Knum/Kden;
% Update state and error covariance
x = full(genmod('fd',t,x0,u,th0,w)) + K*(y-C*x0);
Px = (A-K*C)*P0x*(A-K*C)' + Q + K*R*K';
% Update error sensitivity covariance
Kr = zeros(nx,ny,nth);
Ps = P0s;
Psx = P0sx;
for p = 1:nth
    Athp = reshape(Ath_col(:,p),nx,nx);
    P0sp  = P0s(:,:,p);
    P0sxp = P0sx(:,:,p);
    P0xsp = P0sxp';
    
    Psx(:,:,p) = (Athp-Kr(:,:,p)*C)*P0x*(A-K*C)' + (A-K*C)*P0sxp*(A-K*C)' ...
               + Kr(:,:,p)*R*K' + K*S*K';
    Ps(:,:,p) = (Athp-Kr(:,:,p)*C)*P0x*(Athp-Kr(:,:,p)*C)' + (A-K*C)*P0sp*(A-K*C)' ...
        + (Athp-Kr(:,:,p)*C)*P0xsp*(A-K*C)' + (A-K*C)*P0sxp*(Athp-Kr(:,:,p)*C)' ...
        + Kr(:,:,p)*R*Kr(:,:,p)' + Kr(:,:,p)*S*K' + K*S'*Kr(:,:,p)' + K*V*K';
end

end