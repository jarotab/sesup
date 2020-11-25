function sesupPlotResults(app,ax,validKfId,xn,thn,Type,YScale,XScale,plot_opt)

if numel(validKfId)<1
    validKfId = 1;
end

leg_opt = {'Location','northeast','NumColumns',min([5,ceil(numel(validKfId)/2)])};

if strcmp(XScale,'Log')
    t_s = log10(app.t);
    tLabel = 'Samples in log. scale (-)';
else
    t_s = app.t;
    tLabel = 'Samples (-)';
end

if strcmp(YScale,'Log')
    x_mu_s     = 10*log10(app.x_mu);
    x_est_mu_s = 10*log10(app.x_est_mu);
    s_est_mu_s = 10*log10(app.s_est_mu);
    rmse_s = 10*log10(app.rmse_est);
    mad_s = 10*log10(app.mad_est);
    mse_s = 10*log10(app.mse_est);
    Kx_mu_s     = 10*log10(app.Kx_mu);
    yUnit = "(dB)";
else
    x_mu_s = app.x_mu;
    x_est_mu_s = app.x_est_mu;
    s_est_mu_s = app.s_est_mu;
    rmse_s = app.rmse_est;
    mad_s = app.mad_est;
    mse_s = app.mse_est;
    Kx_mu_s = app.Kx_mu;
    yUnit = "(-)";
end
Px_mu_s = app.Px_mu;
Ps_mu_s = app.Ps_mu;
Psx_mu_s = app.Psx_mu;
Kth_mu_s = app.Kth_mu;

switch Type
    case 'State'
        plot(ax,t_s,x_mu_s(xn(1),:,thn),plot_opt{:})
        for i = 1:numel(validKfId)
            plot(ax,t_s,x_est_mu_s(xn(1),:,validKfId(i),thn),plot_opt{:})
        end
        legend(ax,'Data',app.KfTab.FullName{validKfId},leg_opt{:})
        xlabel(ax,tLabel)
        ylabel(ax,strcat("Amplitude ",yUnit));
    case 'Sensitivity'
        for i = 1:numel(validKfId)
            plot(ax,t_s,reshape(s_est_mu_s(xn(1),1,:,validKfId(i),thn),1,numel(t_s)),plot_opt{:})
        end
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
        xlabel(ax,tLabel)
        ylabel(ax,strcat("Amplitude ",yUnit));
    case 'Norm (sensitivity)'
        for i = 1:numel(validKfId)
            s_tmp = reshape(s_est_mu_s(:,1,:,validKfId(i),thn),app.nx,numel(t_s));
            norm_s = sqrt(diag(s_tmp'*s_tmp));
            plot(ax,t_s,norm_s,plot_opt{:})
        end
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
        xlabel(ax,tLabel)
        ylabel(ax,strcat("Amplitude ",yUnit));
    case '2nd moment (state)'
        for i = 1:numel(validKfId)
            plot(ax,t_s,reshape(Px_mu_s(xn(1),xn(2),:,validKfId(i),thn),1,numel(t_s)),plot_opt{:})
        end
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
        xlabel(ax,tLabel)
        ylabel(ax,strcat("Amplitude ",yUnit));
    case 'Trace 2nd moment (state)'
        for i = 1:numel(validKfId)
            trP = t_s*0;
            for di = 1:app.nx
                trP = trP + reshape(Px_mu_s(di,di,:,validKfId(i),thn),1,numel(t_s));
            end
            plot(ax,t_s,trP,plot_opt{:})
        end
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
        xlabel(ax,tLabel)
        ylabel(ax,strcat("Amplitude ",yUnit));
    case 'Trace 2nd moment (sensitivity)'
        for i = 1:numel(validKfId)
            trP = t_s*0;
            for di = 1:app.nx
                trP = trP + reshape(Ps_mu_s(di,di,:,validKfId(i),thn),1,numel(t_s));
            end
            plot(ax,t_s,trP,plot_opt{:})
        end
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
        xlabel(ax,tLabel)
        ylabel(ax,strcat("Amplitude ",yUnit));
    case '2nd moment (sensitivity)'
        for i = 1:numel(validKfId)
            plot(ax,t_s,reshape(Ps_mu_s(xn(1),xn(2),:,validKfId(i),thn),1,numel(t_s)),plot_opt{:})
        end
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
        xlabel(ax,tLabel)
        ylabel(ax,strcat("Amplitude ",yUnit));
    case '2nd moment (cross)'
        for i = 1:numel(validKfId)
            plot(ax,t_s,reshape(Psx_mu_s(xn(1),xn(2),:,validKfId(i),thn),1,numel(t_s)),plot_opt{:})
        end
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
        xlabel(ax,tLabel)
        ylabel(ax,strcat("Amplitude ",yUnit));
    case 'Gain (sensitivity)'
        for i = 1:numel(validKfId)
            plot(ax,t_s,reshape(Kth_mu_s(xn(1),xn(2),:,validKfId(i),thn),1,numel(t_s)),plot_opt{:})
        end
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
        xlabel(ax,tLabel)
        ylabel(ax,strcat("Gain sensitivity ",yUnit));
    case 'Gain'
        for i = 1:numel(validKfId)
            plot(ax,t_s,reshape(Kx_mu_s(xn(1),xn(2),:,validKfId(i),thn),1,numel(t_s)),plot_opt{:})
        end
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
        xlabel(ax,tLabel)
        ylabel(ax,strcat("Gain ",yUnit));
    case 'MSE'
        for i = 1:numel(validKfId)
            plot(ax,t_s,mse_s(xn(1),:,validKfId(i),thn),plot_opt{:})
        end
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
        xlabel(ax,tLabel)
        ylabel(ax,strcat("MSE ",yUnit));
    case 'MAD'
        for i = 1:numel(validKfId)
            plot(ax,t_s,mad_s(xn(1),:,validKfId(i),thn),plot_opt{:})
        end
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
        xlabel(ax,tLabel)
        ylabel(ax,strcat("MAD ",yUnit));
    case 'RMSE'
        for i = 1:numel(validKfId)
            plot(ax,t_s,rmse_s(xn(1),:,validKfId(i),thn),plot_opt{:})
        end
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
        xlabel(ax,tLabel)
        ylabel(ax,strcat("RMSE ",yUnit));
    case 'Final RMSE'
        for i = 1:numel(validKfId)
            plot(ax,app.ParReal,reshape(rmse_s(xn(1),end,validKfId(i),:),1,numel(app.ParReal)),'-o',plot_opt{:})
        end
        xlabel(ax,'Parameter Value')
        ylabel(ax,strcat("Total RMSE ",yUnit));
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
    case 'Final MAD'
        for i = 1:numel(validKfId)
            plot(ax,app.ParReal,reshape(mad_s(xn(1),end,validKfId(i),:),1,numel(app.ParReal)),'-o',plot_opt{:})
        end
        xlabel(ax,'Parameter Value')
        ylabel(ax,strcat("MAD ",yUnit));
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
    case 'Final Norm (sensitivity)'
        for i = 1:numel(validKfId)
            norm_i = zeros(1,numel(app.ParReal));
            for j = 1:numel(app.ParReal)
                s_tmp = reshape(s_est_mu_s(:,1,end,validKfId(i),j),app.nx,1);
                norm_i(j) = sqrt(diag(s_tmp'*s_tmp));
            end
            plot(ax,app.ParReal,norm_i,'-o',plot_opt{:})
        end
        xlabel(ax,'Parameter Value')
        ylabel(ax,strcat("Norm ",yUnit))
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
    case 'Final MAD Zero (sensitivity)'
        for i = 1:numel(validKfId)
            madz_i = zeros(1,numel(app.ParReal));
            for j = 1:numel(app.ParReal)
                madz_i(j) = median(abs(reshape(s_est_mu_s(xn(1),1,:,validKfId(i),j),1,numel(t_s))));
            end
            plot(ax,app.ParReal,madz_i,'-o',plot_opt{:})
        end
        xlabel(ax,'Parameter Value')
        ylabel(ax,'median(|X|)')
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
    case 'Final Gain'
        for i = 1:numel(validKfId)
            plot(ax,app.ParReal,reshape(Kx_mu_s(xn(1),xn(2),end,validKfId(i),:),1,numel(app.ParReal)),'-o',plot_opt{:})
        end
        xlabel(ax,'Parameter Value')
        ylabel(ax,strcat("Gain ",yUnit))
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
    case 'Final Gain (sensitivity)'
        for i = 1:numel(validKfId)
            plot(ax,app.ParReal,reshape(Kth_mu_s(xn(1),xn(2),end,validKfId(i),:),1,numel(app.ParReal)),'-o',plot_opt{:})
        end
        xlabel(ax,'Parameter Value')
        ylabel(ax,strcat("Gain sensitivity ",yUnit))
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
    case 'Final 2nd moment (state)'
        for i = 1:numel(validKfId)
            plot(ax,app.ParReal,reshape(Px_mu_s(xn(1),xn(2),end,validKfId(i),:),1,numel(app.ParReal)),'-o',plot_opt{:})
        end
        xlabel(ax,'Parameter Value')
        ylabel(ax,strcat("Amplitude ",yUnit))
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
    case 'Final Trace (state)'
        for i = 1:numel(validKfId)
            trP = app.ParReal*0;
            for thi = 1:numel(app.ParReal)
                trP(thi)=trace(Px_mu_s(:,:,end,validKfId(i),thi));
            end
            plot(ax,app.ParReal,trP,'-o',plot_opt{:})
        end
        xlabel(ax,'Parameter Value')
        ylabel(ax,strcat("Amplitude ",yUnit))
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
    case 'Final Trace (sensitivity)'
        for i = 1:numel(validKfId)
            trP = app.ParReal*0;
            for thi = 1:numel(app.ParReal)
                trP(thi)=trace(Ps_mu_s(:,:,end,validKfId(i),thi));
            end
            plot(ax,app.ParReal,trP,'-o',plot_opt{:})
        end
        xlabel(ax,'Parameter Value')
        ylabel(ax,strcat("Amplitude ",yUnit))
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
    case 'Final 2nd moment (sensitivity)'
        for i = 1:numel(validKfId)
            plot(ax,app.ParReal,reshape(Ps_mu_s(xn(1),xn(2),end,validKfId(i),:),1,numel(app.ParReal)),'-o',plot_opt{:})
        end
        xlabel(ax,'Parameter Value')
        ylabel(ax,strcat("Amplitude ",yUnit))
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
    case 'Final 2nd moment (cross)'
        for i = 1:numel(validKfId)
            plot(ax,app.ParReal,reshape(Psx_mu_s(xn(1),xn(2),end,validKfId(i),:),1,numel(app.ParReal)),'-o',plot_opt{:})
        end
        xlabel(ax,'Parameter Value')
        ylabel(ax,strcat("Amplitude ",yUnit))
        legend(ax,app.KfTab.FullName{validKfId},leg_opt{:})
end
end