function [mse_est,rmse_est,x_mu,x_est_mu,Kx_mu,Px_mu,Kth_mu,Ps_mu,Psx_mu] = sesupsimulate(app)

mse_est  = app.mse_est;
rmse_est = app.rmse_est;
x_mu     = app.x_mu;
x_est_mu = app.x_est_mu;
Kx_mu    = app.Kx_mu;
Kth_mu   = app.Kth_mu;
Px_mu    = app.Px_mu;
Ps_mu    = app.Ps_mu;
Psx_mu   = app.Psx_mu;

nt = numel(app.t);
nthr = numel(app.ParReal);
if any(strcmp(app.KfTab.FilterName,'DKF'))
    Ks = sym('k',[app.nx,app.ny]);
else
    Ks = zeros(app.nx,app.ny);
end

% Time to end estimation
MeanSimTime = 0;
ProgressDialog = uiprogressdlg(app.UIFigure,'Title','Simulation in progress',...
    'Message','Starting simulation');
ProgressDialog.Value = 0;

% Local variables definitions due to parallel toolbox
t = app.t;
x0 = app.x0;
th0 = app.th0;
P0 = app.P0;
y0 = app.y0;
R = app.R;
Q = app.Q;
nx = app.nx;
nth = app.nth;
ny = app.ny;
nu = app.nu;
nf = app.nf;
FilterName = app.KfTab.FilterName;
Par = app.KfTab.Par;

% Additional nosie parameters for XDKF
Rss = R*10;
Rys = R;

for thi = 1:numel(app.ParReal)
    
    tmp_mse_est = mse_est(:,:,:,thi);
    tmp_rmse_est = rmse_est(:,:,:,thi);
    tmp_x_mu = x_mu(:,:,thi);
    tmp_x_est_mu = x_est_mu(:,:,:,thi);
    tmp_Kx_mu = Kx_mu(:,:,:,:,thi);
    tmp_Kth_mu = Kth_mu(:,:,:,:,thi);
    tmp_Px_mu = Px_mu(:,:,:,:,thi);
    tmp_Ps_mu = Ps_mu(:,:,:,:,thi);
    tmp_Psx_mu = Psx_mu(:,:,:,:,thi);
    
    
    th_sim = app.ParReal(:,thi);
    
    ShowText = sprintf('Running %d experiments with parameter id %d/%d', app.ExpNum, thi, numel(app.ParReal));
    ProgressDialog.Message = ShowText;    
    if thi>1
        ShowText = sprintf('%s\n Estimated time to end:  %s',ShowText,datestr(TimeToFinish/86400,'HH:MM:SS'));
        ProgressDialog.Message = ShowText;
    end
    
    % Intial variables   
    x_est_init = zeros(nx,1,nt,nf);
    P_est_init = zeros(nx,nx,nt,nf);
    Kx_init = zeros(nx,ny,nt,nf);
    y_init = zeros(ny,1,nt);
    x_init = zeros(nx,1,nt);
    u_init = app.u;
    x_init(:,:,1) = x0;
    for kfi = 1:nf
        x_est_init(:,:,1,kfi) = x0;
        P_est_init(:,:,1,kfi) = P0;
    end
    Pz_skf_init =  repmat(blkdiag(P_est_init(:,:,1,1),eye(nth)),1,1,nf);
    for d = 1:nf
        if strcmp(FilterName{d},'SKF')
            Pz_skf_init(:,:,d) =  blkdiag(P_est_init(:,:,1,1),diag(Par(d)));
        end
    end
    Ps_xdkf_init  = repmat(P0,1,1,nth,nt,nf);
    Psx_xdkf_init = repmat(P0,1,1,nth,nt,nf);
    s_dkf_init = zeros(nx,nth,nf);
    S_sdkf_init = zeros(nx,nth,nf);    
    Kr_xdkf_init = zeros(nx,ny,nth,nt,nf);
    
    tic
    parfor n = 1:app.ExpNum
%     for n = 1:app.ExpNum
        % Intialize data matrices/vectors
        x_est    = x_est_init;
        P_est    = P_est_init;
        Kx       = Kx_init;
        y        = y_init;
        x        = x_init;
        u        = u_init;        
        % Additional variables
        Pz_skf   = Pz_skf_init;
        Ps_xdkf  = Ps_xdkf_init;
        Psx_xdkf = Psx_xdkf_init;
        Kr_xdkf  = Kr_xdkf_init; 
        s_dkf    = s_dkf_init;
        S_sdkf   = S_sdkf_init;
        
        % Simulate experiment
        for k = 1:nt
            
            % Data generator
            y(:,:,k) = full(genmod('g',k,x(:,:,k),u(:,k),th_sim,chol(R)'*randn(ny,1)));
            if k<nt
            x(:,:,k+1) = full(genmod('fd',k,x(:,:,k),u(:,k),th_sim,chol(Q)'*randn(nx,1)));
            end
            
            Tmpx = x_est(:,:,k,:);
            TmpPx = P_est(:,:,k,:);
            TmpPs = Ps_xdkf(:,:,:,k,:);
            TmpPsx = Psx_xdkf(:,:,:,k,:);
            for kfi = 1:nf
                KfPar = Par(kfi);
                switch FilterName(kfi)
                    case "Optimal KF"
                        % Optimal KF
                        [Tmpx(:,:,1,kfi),TmpPx(:,:,1,kfi),Kx(:,:,k,kfi)] = kfstep(Tmpx(:,:,1,kfi),th_sim,TmpPx(:,:,1,kfi),k,u(:,k),y(:,:,k),R,Q);
                    case "KF"
                        % KF
                        [Tmpx(:,:,1,kfi),TmpPx(:,:,1,kfi),Kx(:,:,k,kfi)] = kfstep(Tmpx(:,:,1,kfi),th0,TmpPx(:,:,1,kfi),k,u(:,k),y(:,:,k),R,Q);
                    case "SKF"
                        % SKF
                        [Tmpx(:,:,1,kfi),TmpPx(:,:,1,kfi),Kx(:,:,k,kfi),Pz_skf(:,:,kfi)] ...
                            = skfstep(Tmpx(:,:,1,kfi),th0,Pz_skf(:,:,kfi),k,u(:,k),y(:,:,k),R,Q);
                    case "DKF"
                        % DKF
                        [Tmpx(:,:,1,kfi),TmpPx(:,:,1,kfi),Kx(:,:,k,kfi),s_dkf(:,:,kfi)] ...
                            = dkfstep(Tmpx(:,:,1,kfi),th0,TmpPx(:,:,1,kfi),s_dkf(:,:,kfi),k,u(:,k),y(:,:,k),R,Q,KfPar,Ks);
                    case "SDKF"
                        % SDKF
                        [Tmpx(:,:,1,kfi),TmpPx(:,:,1,kfi),Kx(:,:,k,kfi),S_sdkf(:,:,kfi)] ...
                            = sdkfstep(Tmpx(:,:,1,kfi),th0,TmpPx(:,:,1,kfi),S_sdkf(:,:,kfi),k,u(:,k),y(:,:,k),R,Q,KfPar);
                    case "XDKF"
                        % XDKF
                        [Tmpx(:,:,1,kfi),TmpPx(:,:,1,kfi),Kx(:,:,k,kfi),TmpPs(:,:,:,1,kfi),TmpPsx(:,:,:,1,kfi),Kr_xdkf(:,:,:,k,kfi)] ...
                            = xdkfstep(Tmpx(:,:,1,kfi),th0,TmpPx(:,:,1,kfi),TmpPs(:,:,:,1,kfi),TmpPsx(:,:,:,1,kfi),k,u(:,k),y(:,:,k),R,Q,Rss,Rys,KfPar);
                    case "XDKF-S1"
                        % XDKF-S1
                        [Tmpx(:,:,1,kfi),TmpPx(:,:,1,kfi),Kx(:,:,k,kfi),TmpPs(:,:,:,1,kfi),TmpPsx(:,:,:,1,kfi),Kr_xdkf(:,:,:,k,kfi)] ...
                            = xdkf1step(Tmpx(:,:,1,kfi),th0,TmpPx(:,:,1,kfi),TmpPs(:,:,:,1,kfi),TmpPsx(:,:,:,1,kfi),k,u(:,k),y(:,:,k),R,Q,Rss,Rys,KfPar);
                    case "XDKF-S2"
                        % XDKF-S2
                        [Tmpx(:,:,1,kfi),TmpPx(:,:,1,kfi),Kx(:,:,k,kfi),TmpPs(:,:,:,1,kfi),TmpPsx(:,:,:,1,kfi),Kr_xdkf(:,:,:,k,kfi)] ...
                            = xdkf2step(Tmpx(:,:,1,kfi),th0,TmpPx(:,:,1,kfi),TmpPs(:,:,:,1,kfi),TmpPsx(:,:,:,1,kfi),k,u(:,k),y(:,:,k),R,Q,Rss,Rys,KfPar);
                    case "Hinf"
                        % Hinf (TmpPx := M_est)
                        if k < 2
                            % M0 definition
                            TmpPx(:,:,1,kfi) = inv(inv(TmpPx(:,:,1,kfi)) + KfPar^(-2)*eye(nx));
                        end
                        [Tmpx(:,:,1,kfi),TmpPx(:,:,1,kfi),Kx(:,:,k,kfi)] = hinfstep(Tmpx(:,:,1,kfi),th0,TmpPx(:,:,1,kfi),k,u(:,k),y(:,:,k),KfPar);
                end
            end
            
            if k<nt
            x_est(:,:,k+1,:)  = Tmpx;
            P_est(:,:,k+1,:)  = TmpPx;
            Ps_xdkf(:,:,:,k+1,:)  = TmpPs;
            Psx_xdkf(:,:,:,k+1,:)  = TmpPsx;
            end
            
        end
        
        x = reshape(x,nx,nt);
        x_est = reshape(x_est,nx,nt,nf);
        
        % RMSE Sum
        tmp_mse_est   = tmp_mse_est + cumsum((x-x_est).^2,2)./(t+1);        
        % RMSE Sum
        tmp_rmse_est  = tmp_rmse_est + sqrt(cumsum((x-x_est).^2,2)./(t+1));        
        % est Sum
        tmp_x_mu      = tmp_x_mu + x;
        tmp_x_est_mu  = tmp_x_est_mu + x_est;
        tmp_Kx_mu     = tmp_Kx_mu + Kx;
        tmp_Px_mu     = tmp_Px_mu + P_est;
        tmp_Kth_mu    = tmp_Kth_mu + reshape(Kr_xdkf(:,:,1,:,:),nx,ny,nt,nf);
        tmp_Ps_mu     = tmp_Ps_mu + reshape(Ps_xdkf(:,:,1,:,:),nx,nx,nt,nf);
        tmp_Psx_mu    = tmp_Psx_mu + reshape(Psx_xdkf(:,:,1,:,:),nx,nx,nt,nf);
        
    end
    
    % Mean in experiments
    mse_est(:,:,:,thi)  = tmp_mse_est(:,:,:)/app.ExpNum;
    rmse_est(:,:,:,thi) = tmp_rmse_est(:,:,:)/app.ExpNum;
    x_mu(:,:,thi)       = tmp_x_mu(:,:)/app.ExpNum;
    x_est_mu(:,:,:,thi) = tmp_x_est_mu(:,:,:)/app.ExpNum;
    Kx_mu(:,:,:,:,thi)  = tmp_Kx_mu(:,:,:,:)/app.ExpNum;
    Px_mu(:,:,:,:,thi)  = tmp_Px_mu(:,:,:,:)/app.ExpNum;
    Kth_mu(:,:,:,:,thi)  = tmp_Kth_mu(:,:,:,:)/app.ExpNum;
    Ps_mu(:,:,:,:,thi)  = tmp_Ps_mu(:,:,:,:)/app.ExpNum;
    Psx_mu(:,:,:,:,thi) = tmp_Psx_mu(:,:,:,:)/app.ExpNum;
    
    SimTime = toc;
    MeanSimTime = (MeanSimTime*(thi-1) + SimTime)/thi;
    TimeToFinish = (nthr-thi)*MeanSimTime;      
    ProgressDialog.Value = 1 - (nthr-thi)/nthr;
    
end
end

