function fd = discretizeStochastic(fc,inputs,input_names,ix,opt)
    % !!!!! Assumes that the last input is a stocatic term that is not in the definition of fc
    
        import casadi.*
        
        Td = opt.Ts / opt.SubSamples;
        
        fc = Function('fc',inputs,{fc(inputs{1:end-1})},input_names,{'fc'}); 
        f = Function('f',inputs,{inputs{ix}},input_names,{'f'}); % auxliary function for sub-sampling 
        for isubsample = 1 : opt.SubSamples
            switch lower(opt.DiscretizationMethod)
                case 'euler'
                    % fd = Function('fd',{x,u},...
                    %      {f(x,u) + Td*fc(f(x,u),u)},...
                    %      {'x','u'},{'fd'});
                    inputs1 = inputs;
                    inputs1{ix} = f(inputs{:});
                    fd = Function('fd',inputs,...
                        { f(inputs{:}) + Td*fc(inputs1{:}) + sqrt(Td)*inputs{end}},...
                        input_names,{'fd'});                                      
            end
            f = Function('f',inputs,{fd(inputs{:})},input_names,{'f'});
        end
        
end