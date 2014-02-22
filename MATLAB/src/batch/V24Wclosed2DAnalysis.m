% Adam Gleichsner (amg188@case.edu)
% Case Western Reserve University

function V24Wclosed2DAnalysis()

x = 0;
y = 0;

seed = '1';

f3 = figure('Position', [1 1 900 800], 'Toolbar', 'figure', 'Visible', 'off', 'MenuBar', 'none', 'Name', 'Closed 2D Analysis');
movegui(f3, 'center');
axhan3 = axes('Units', 'Pixels', 'Position', [50, 100 , 600, 600]);
staticN = uicontrol('Style', 'text', 'Position', [675 700 100 30], 'String', 'Steps(N) =');
editN = uicontrol('Style', 'edit', 'String', '1000', 'Position', [775 700 100 31]);
staticD = uicontrol('Style', 'text', 'Position', [675 650 100 30], 'String', 'Diffusion Coef(D) =');
editD = uicontrol('Style', 'edit', 'String', '.000045', 'Position', [775 650 100 31]);
staticH = uicontrol('Style', 'text', 'Position', [675 550 100 30], 'String', 'Batch Size = ');
editH = uicontrol('Style', 'edit', 'String', '100', 'Position', [775 550 100 31]);
staticW = uicontrol('Style', 'text', 'Position', [675 500 50 30], 'String', 'W = ');
editW = uicontrol('Style', 'edit', 'String', '.00001', 'Position', [725 500 50 31]);
staticL = uicontrol('Style', 'text', 'Position', [775 500 50 30], 'String', 'L = ');
editL = uicontrol('Style', 'edit', 'String', '.00001', 'Position', [825 500 50 31]);
staticUpWin = uicontrol('Style', 'text', 'Position', [675 450 50 30], 'String', 'UpWin Width = ');
editUpWin = uicontrol('Style', 'edit', 'String', '.000002', 'Position', [725 450 50 31]);
staticBotWin = uicontrol('Style', 'text', 'Position', [775 450 50 30], 'String', 'BotWin Width = ');
editBotWin = uicontrol('Style', 'edit', 'String', '.000002', 'Position', [825 450 50 31]);
staticeW = uicontrol('Style', 'text', 'Position', [675 300 50 30], 'String', 'eW = ');
editeW = uicontrol('Style', 'edit', 'String', '.01', 'Position', [725 300 50 31]);
staticeL = uicontrol('Style', 'text', 'Position', [775 300 50 30], 'String', 'eL = ');
editeL = uicontrol('Style', 'edit', 'String', '.01', 'Position', [825 300 50 31]);
buttonHist = uicontrol('Style', 'pushbutton', 'String', 'Histogram', 'Position',[675 350 200 40], 'Callback', @histogram);
buttonCOT = uicontrol('Style', 'pushbutton', 'String', 'Change Over Time', 'Position',[675 250 200 40], 'Callback', @changeOverTime);
buttonWCOT = uicontrol('Style', 'pushbutton', 'String', 'Window Change Over Time', 'Position',[675 150 200 40], 'Callback', @windowChangeOverTime);
staticetw = uicontrol('Style', 'text', 'Position', [675 200 50 30], 'String', 'etw = ');
editetw = uicontrol('Style', 'edit', 'String', '.000005', 'Position', [725 200 50 31]);
staticebw = uicontrol('Style', 'text', 'Position', [775 200 50 30], 'String', 'ebw = ');
editebw = uicontrol('Style', 'edit', 'String', '.000005', 'Position', [825 200 50 31]);
buttontau = uicontrol('Style', 'pushbutton', 'String', 'MET as function of Tau', 'Position',[675 50 200 40], 'Callback', @tauChange);
staticinitTau = uicontrol('Style', 'text', 'Position', [675 100 50 30], 'String', 'init Tau = ');
editinittau = uicontrol('Style', 'edit', 'String', '100', 'Position', [725 100 50 31]);
staticfintau = uicontrol('Style', 'text', 'Position', [775 100 50 30], 'String', 'fin Tau = ');
editfintau = uicontrol('Style', 'edit', 'String', '.001', 'Position', [825 100 50 31]);

set(f3, 'Visible', 'on');

title('Particle Track of a Single Simulated Particle');
xlabel('Average Time ');
ylabel('Y Position');


    function histogram(source, eventdata)
        tic
        
        batch = [100 250 500 750 1000 2500 5000 7500 10000 25000 50000 75000 100000 250000 500000 750000 1000000];
            
        for a = 1:length(batch)
        N = str2num(get(editN, 'String'));
        D = str2num(get(editD, 'String'));
        L = str2num(get(editL, 'String'));
        W = str2num(get(editW, 'String'));
        wTop = str2num(get(editUpWin, 'String'));
        wBot = str2num(get(editBotWin, 'String'));
        randArray = randn(N, 2);
        
        histArray = [];
        
        dimensions = 2;         % two dimensional simulation
        tau = .0001;               % time interval in seconds
        time = tau * 1:N;       % create a time vector for plotting
        
        R    = .145e-6;              % radius in meters 
        eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K 
        kB   = 1.38e-23;            % Boltzmann constant T    = 293;                 % Temperature in degrees Kelvin
        eta  = 1.827e-5;              % viscosity of air in SI units (Pascal-seconds) at 293 K 
        kB   = 1.38e-23;            % Boltzmann constant
        T    = 293;                 % Temperature in degrees Kelvin
        
        D    = kB * T / (6 * pi * eta * R);
        k = sqrt(D * dimensions * tau);
        progressbar;
        
        for j = 1:batch(a)
            
            i = 1;
            randomChoose = 0 + (2-(0))*rand(1,1);
            if randomChoose >= 1
                x = -W/2 + (W/2-(-W/2))*rand(1,1);
                y = L/2;
            else
                x = 0;
                y = 0 + (L-(0))*rand(1,1);
            end
            
            while x(i) > -W && x(i) < W && y(i) < L + L/2 && y(i) > -L/2
            
                i = i + 1;

                dx = k * randn(1,1);
                dy = k * randn(1,1);

                x = [x, x(i - 1) + dx];
                y = [y, y(i - 1) + dy];

                if x(i) > -W && x(i) < W && y(i) < L + L/2 && y(i) > -L/2
                    if ~(y(i) < L/2 + wTop/2 && y(i) > L/2 - wTop/2) && ~(y(i) < L + L/2 && y(i) > L + L/2 - wTop/2) && ~(y(i) > -L/2 && y(i) < -L/2 + wTop/2)
                        if (x(i) < -W/2 && x(i-1) > -W/2) || (x(i) > -W/2 && x(i-1) < -W/2)
                            if x(i) < -W/2
                                x(i) = x(i) + 2 * (abs(x(i)) - W/2);
                            else
                                x(i) = x(i) - 2 * (abs(x(i)) - W/2);
                            end
                        elseif (x(i) > W/2 && x(i-1) < W/2) || (x(i) < W/2 && x(i-1) > W/2)
                            if x(i) > W/2
                                x(i) = x(i) - 2 * (abs(x(i)) - W/2);
                            else
                                x(i) = x(i) + 2 * (abs(x(i)) - W/2);
                            end
                        end
                    end

                    if ~(x(i) > -wTop/2 && x(i) < wTop/2) && ~(x(i) > -W && x(i) < -W + wBot/2) && ~(x(i) < W && x(i) >  W - wBot/2)
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
            end
            histArray = [histArray, i];
            progressbar(a/length(batch), j/batch(a));
        end
        
        histArray = histArray * tau;
        bins = 0.00005:size(histArray)/20:max(histArray);
        figure(f3);
        hist(histArray, bins);
        h = findobj(gca,'Type','patch');
        set(h, 'EdgeColor','w')
        
        title(sprintf('Occurances versus Time for \nLarge Batch Diffusion (AVG = %f)', mean(histArray)));
        xlabel('Time');
        ylabel('Number of Occurances');
        
        filename = strcat('./2D Analysis/iMET/Square/4 Window/Eta = 1.000e-3/Group 1/hist', num2str(batch(a)), '_1.fig');
        saveas(f3, filename, 'fig');
        
        end
        
        toc
    end

    function changeOverTime(source, eventdata)
        tic
        N = str2num(get(editN, 'String'));
        D = str2num(get(editD, 'String'));
        sL = str2num(get(editL, 'String'));
        sW = str2num(get(editW, 'String'));
        eL = str2num(get(editeL, 'String'));
        eW = str2num(get(editeW, 'String'));
        wTop = str2num(get(editUpWin, 'String'));
        wBot = str2num(get(editBotWin, 'String'));
        randArray = randn(N, 2);
        
        step = .0005;
        histArray = [];
        avgArray = [];
        
        progressbar('Overall', 'Step Test');

        dimensions = 2;         % two dimensional simulation
        tau = .001;               % time interval in seconds
        time = tau * 1:N;       % create a time vector for plotting
        
        R    = .145e-6;              % radius in meters eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K kB   = 1.38e-23;            % Boltzmann constant T    = 293;                 % Temperature in degrees Kelvin
        eta  = 1.827e-5;              % viscosity of air in SI units (Pascal-seconds) at 293 K 
        kB   = 1.38e-23;            % Boltzmann constant
        T    = 293;                 % Temperature in degrees Kelvin
        
        D    = kB * T / (6 * pi * eta * R);
        k = sqrt(D * dimensions * tau);
        
        L = sL;
        W = sW;
        
        for yy = sL:step:eL
            
            L = yy;
            W = yy;
            
            for j = 1:str2num(get(editH, 'String'))
                i = 1;
                x = 0;
                y = 0;
                while ~(x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= L)
                    
                    i = i + 1;
                    
                    dx = k * randn(1,1);
                    dy = k * randn(1,1);
                    
                    x = [x, x(i - 1) + dx];
                    y = [y, y(i - 1) + dy];
                    
                    if ~(x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= L)
                        if x(i) < -W/2
                            x(i) = x(i) + 2 * (abs(x(i)) - W/2);
                        elseif x(i) > W/2
                            x(i) = x(i) - 2 * (abs(x(i)) - W/2);
                        end
                        
                        if y(i) < 0
                            y(i) = y(i) + 2 * abs(y(i));
                        elseif y(i) > L
                            y(i) = y(i) - 2 * (abs(y(i) - L));
                        end
                    end
                    progressbar(yy/eL, j/str2num(get(editH, 'String')));                
                end
                histArray = [histArray, i];
            end
            
            avgArray = [avgArray, mean(histArray)];
            histArray = [];
        end
        progressbar(1, 1);
        avgArray = avgArray * tau
        plot([sW:step:eW], avgArray);
        
        title('Average Time of Diffusion versus Change in Void Dimension');
        xlabel('Current Dimension Length (W = , L =)');
        ylabel('Average Time per Batch');
        
        toc
    end

    function tauChange(source, eventdata)
        tic
        N = str2num(get(editN, 'String'));
        D = str2num(get(editD, 'String'));
        L = str2num(get(editL, 'String'));
        W = str2num(get(editW, 'String'));
        wTop = str2num(get(editUpWin, 'String'));
        wBot = str2num(get(editBotWin, 'String'));
        randArray = randn(N, 2);
                endTau = str2num(get(editfintau, 'String'));

        step = -.0001;
        histArray = [];
        avgArray = [];
        
        progressbar('Overall', 'Step Test');

        dimensions = 2;         % two dimensional simulation
        initTau = str2num(get(editinittau, 'String'));               % time interval in seconds
        endTau = str2num(get(editfintau, 'String'));
        %time = tau * 1:N;       % create a time vector for plotting
        
        R    = .145e-6;              % radius in meters eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K kB   = 1.38e-23;            % Boltzmann constant T    = 293;                 % Temperature in degrees Kelvin
        eta  = 1.827e-5;              % viscosity of air in SI units (Pascal-seconds) at 293 K 
        kB   = 1.38e-23;            % Boltzmann constant
        T    = 293;                 % Temperature in degrees Kelvin
        
        D    = kB * T / (6 * pi * eta * R);
        
        tauC = initTau;
        x1 = [];
        while(tauC > (1.01* endTau))
            
            k = sqrt(D * dimensions * tauC);
            
            for j = 1:str2num(get(editH, 'String'))
                i = 1;
                x = 0;
                y = 0;
                while ~(x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= L)
                    
                    i = i + 1;
                    
                    dx = k * randn(1,1);
                    dy = k * randn(1,1);
                    
                    x = [x, x(i - 1) + dx];
                    y = [y, y(i - 1) + dy];
                    
                    if ~(x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= L)
                        if x(i) < -W/2
                            x(i) = x(i) + 2 * (abs(x(i)) - W/2);
                        elseif x(i) > W/2
                            x(i) = x(i) - 2 * (abs(x(i)) - W/2);
                        end
                        
                        if y(i) < 0
                            y(i) = y(i) + 2 * abs(y(i));
                        elseif y(i) > L
                            y(i) = y(i) - 2 * (abs(y(i) - L));
                        end
                    end
                    progressbar((tauC - initTau)/(endTau - initTau) - .01, j/str2num(get(editH, 'String')));                
                end
                histArray = [histArray, i];
            end
            x1 = [x1, tauC]
            avgArray = [avgArray, mean(histArray) * tauC];
            histArray = [];
            tauC = .5 * tauC;
        end
        progressbar(1, 1);
        plot(x1, avgArray);
        
        title('Average Time of Diffusion versus Change in Void Dimension');
        xlabel('Value of dT');
        ylabel('Average Time per Batch');
        
        toc
    end

    function windowChangeOverTime(source, eventdata)
        
        tic
        N = str2num(get(editN, 'String'));
        D = str2num(get(editD, 'String'));
        %seed = str2num(get(editS, 'String'));
        stW = str2num(get(editUpWin, 'String'));
        sbW = str2num(get(editBotWin, 'String'));
        etW = str2num(get(editetw, 'String'));
        ebW = str2num(get(editebw, 'String'));
        wTop = stW;
        wBot = sbW;
        W = str2num(get(editW, 'String'));
        L = str2num(get(editL, 'String'));        
        
        step = .1e-6;
        histArray = [];
        avgArray = [];
        
        %randn('seed', seed);
        progressbar('Overall', 'Step Test');
        
        dimensions = 2;         % two dimensional simulation
        tau = .001;               % time interval in seconds
        time = tau * 1:N;       % create a time vector for plotting
        
        R    = .145e-6;              % radius in meters eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K kB   = 1.38e-23;            % Boltzmann constant T    = 293;                 % Temperature in degrees Kelvin
        eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K
        kB   = 1.38e-23;            % Boltzmann constant
        T    = 293;                 % Temperature in degrees Kelvin
        
        D    = kB * T / (6 * pi * eta * R);
        k = sqrt(D * dimensions * tau);
        
        for yy = sbW:step:ebW
            
            wTop = yy;
            wBot = yy;
            
            for j = 1:str2num(get(editH, 'String'))
                i = 1;
                x = 0;
                y = 0;
                while ~(x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= L)
                    
                    i = i + 1;
                    
                    dx = k * randn(1,1);
                    dy = k * randn(1,1);
                    
                    x = [x, x(i - 1) + dx];
                    y = [y, y(i - 1) + dy];
                    
                    if ~(x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= L)
                        if x(i) < -W/2
                            x(i) = x(i) + 2 * (abs(x(i)) - W/2);
                        elseif x(i) > W/2
                            x(i) = x(i) - 2 * (abs(x(i)) - W/2);
                        end
                        
                        if y(i) < 0
                            y(i) = y(i) + 2 * abs(y(i));
                        elseif y(i) > L
                            y(i) = y(i) - 2 * (abs(y(i) - L));
                        end
                    end
                end
                histArray = [histArray, i];
                progressbar((yy - sbW)/(ebW - sbW), j/str2num(get(editH, 'String')));
            end
            
            avgArray = [avgArray, mean(histArray)];
            histArray = [];
        end
        progressbar(1, 1);
        avgArray = avgArray * tau
        plot([sbW:step:ebW], avgArray);
        
        title('Average Time of Diffusion versus Change in Void Dimension');
        xlabel('Current Dimension Length (W = , L =)');
        ylabel('Average Time per Batch');
        toc
    end
    
end