function scriptV2closed2DRWHistBatch(num, L, W, wTop, wBot)
        
tic
        batch = [100 250 500 750 1000 2500 5000 7500 10000 25000 50000 75000 100000 250000 500000 750000 1000000];
            
        fprintf('Initializing...');
        mkdir(strcat('/home/amg188/MATLAB/results/Group', num2str(num)))
        for a = 1:length(batch)
        
        num = num;
        L = L;
        W = W;
        wTop = wTop;
        wBot = wBot;
        
        histArray = [];
        
        dimensions = 2;         % two dimensional simulation
        tau = .001;               % time interval in seconds
        
        R    = .145e-6;              % radius in meters 
        eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K 
        kB   = 1.38e-23;            % Boltzmann constant T    = 293;                 % Temperature in degrees Kelvin
        eta  = 1.827e-5;              % viscosity of air in SI units (Pascal-seconds) at 293 K 
        kB   = 1.38e-23;            % Boltzmann constant
        T    = 293;                 % Temperature in degrees Kelvin
        
        D    = kB * T / (6 * pi * eta * R);
        k = sqrt(D * dimensions * tau);
        
        
        for j = 1:batch(a)
            i = 1;
            x = -W/2 + (W/2-(-W/2))*rand(1,1);
            y = L/2;
            
            while ~(x(i) > -W/2 && x(i) < W/2 && y(i) >= L + L/2) && ~(x(i) > -W/2 && x(i) < W/2 && y(i) <= -L/2)
                
                i = i + 1;
                
                dx = k * randn(1,1);
                dy = k * randn(1,1);
                
                x = [x, x(i - 1) + dx];
                y = [y, y(i - 1) + dy];
                
                if ~(x(i) > -wTop/2 && x(i) < wTop/2)
                    if x(i) < -W/2
                        
                        
                        x(i) = x(i) + 2 * (abs(x(i)) - W/2);
                       
                    elseif x(i) > W/2
                        
                        
                        x(i) = x(i) - 2 * (abs(x(i)) - W/2);
                        
                       
                    end
                    
                    if (y(i) < 0 && y(i-1) > 0) || (y(i) > 0 && y(i-1) < 0)
                        
                        
                        if y(i) < 0
                            y(i) = y(i) + 2 * abs(y(i));
                        else
                            y(i) = y(i) - 2 * abs(y(i));
                        end
                        
                    elseif (y(i) > L && y(i-1) < L) || (y(i) < L && y(i-1) > L)
                        
                        if y(i) > L
                            y(i) = y(i) - 2 * (abs(y(i) - L));
                        else
                            y(i) = y(i) + 2 * (abs(y(i) - L));
                        end
                        
                       
                    end
                end
            end
            histArray = [histArray, i];
            if (int8(j/batch(a) * 100) > 0 && int8(j/batch(a) * 100) < 2) || (int8(j/batch(a) * 100) > 24 && int8(j/batch(a) * 100) < 26) || (int8(j/batch(a) * 100) > 49 && int8(j/batch(a) * 100) < 51) ||  (int8(j/batch(a) * 100) > 74 && int8(j/batch(a) * 100) < 76) || (int8(j/batch(a) * 100) > 98 && int8(j/batch(a) * 100) < 101) 
            	fprintf('\rCurrent: %i%% ; Overall: %i%%', int8(j/batch(a) * 100), int8(a/length(batch) * 100));
       	    elseif int8(a/length(batch) * 100) == 25 || int8(a/length(batch) * 100) == 50 || int8(a/length(batch) * 100) == 75 || int8(a/length(batch) * 100) == 100
                fprintf('\rCurrent: %i%% ; Overall: %i%%', int8(j/batch(a) * 100), int8(a/length(batch) * 100));
            end
	end
        
        histArray = histArray * tau;
        bins = 0.00005:size(histArray)/20:max(histArray);
        filename = strcat('/home/amg188/MATLAB/results/Group', ' ', num2str(num), '/hist', num2str(batch(a)), '_1.txt');
        save(filename, 'histArray', 'bins', 'D', 'tau', '-double', '-ascii'); 
        
        end
        toc
    end
