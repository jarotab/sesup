function [x,Px,K,s,Ps,Psx,Kth] = exdkf1step(x0,th0,P0x,s0,P0s,P0sx,t,u,y,R,Q,gamma)

% Extended XDKF 
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
Fx  = full(genmod('dfddx',t,x0,u,th0,w));
fth = full(genmod('dfddth',t,x0,u,th0,w));
Fxth_col = full(genmod('ddfddxdth',t,x0,u,th0,w));
Gx = full(genmod('dgdx',t,x0,u,th0,e));
gth = full(genmod('dgdth',t,x0,u,th0,e));
Gxth_col = full(genmod('ddgdxdth',t,x0,u,th0,e));

f = full(genmod('fd',t,x0,u,th0,w));
g = full(genmod('g',t,x0,u,th0,e));

% weight
alpha = 1-sum(gamma);

% Gain
Kn = alpha*(Fx*P0x*Gx');
Kd = alpha*(Gx*P0x*Gx' + R);
for p = 1:nth
    Fxththp = reshape(Fxth_col(:,p),nx,nx);
    Gxththp = reshape(Gxth_col(:,p),ny,nx);
    Kn = Kn + gamma(p)*(...
        Fxththp*P0x*Gxththp' + Fx*P0s(:,:,p)*Gx'...
        + Fxththp*P0sx(:,:,p)'*Gx' + Fx*P0sx(:,:,p)*Gxththp'...
        - Fx*s0(:,p)*gth(:,p)' - fth(:,p)*s0(:,p)'*Gx' + fth(:,p)*gth(:,p)'...
        );
    Kd = Kd + gamma(p)*(...
        Gxththp*P0x*Gxththp' + Gx*P0s(:,:,p)*Gx'...
        + Gxththp*P0sx(:,:,p)'*Gx' + Gx*P0sx(:,:,p)*Gxththp'...
        - Gx*s0(:,p)*gth(:,p)' - gth(:,p)*s0(:,p)'*Gx' + gth(:,p)*gth(:,p)'...
        );    
end
K = Kn/Kd;
Kth = K*0;

% Update state and error covariance
x = f + K*(y-g);
Fx_KGx = Fx-K*Gx;
Px = Fx_KGx*P0x*Fx_KGx' + Q + K*R*K';

% Update error sensitivity covariance and Sensitivity
Ps = P0s;
Psx = P0sx;
s = s0;
for p = 1:nth
    Fxthp = reshape(Fxth_col(:,p),nx,nx);
    Gxthp = reshape(Gxth_col(:,p),ny,nx);
    P0sp  = P0s(:,:,p);    
    P0sxp = P0sx(:,:,p);    
    Fxthp_KGxthp = Fxthp - K*Gxthp;
    fth_Kgth = fth(:,p) - K*gth(:,p);
    Psx(:,:,p) = Fx_KGx*P0sxp*Fx_KGx' ...
                + Fxthp_KGxthp*P0x*Fx_KGx';    
    Ps(:,:,p) = Fx_KGx*P0sp*Fx_KGx' + Fxthp_KGxthp*P0x*Fxthp_KGxthp'...
        + Fxthp_KGxthp*P0sxp'*Fx_KGx' + Fx_KGx*P0sxp*Fxthp_KGxthp'...
        - Fx_KGx*s0(:,p)*fth_Kgth' - fth_Kgth*s0(:,p)'*Fx_KGx'...
        + fth_Kgth*fth_Kgth';
    s(:,p) = Fx_KGx*s0(:,p) - fth_Kgth;
end

end
