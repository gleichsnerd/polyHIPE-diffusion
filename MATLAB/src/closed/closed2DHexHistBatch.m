% Adam Gleichsner (amg188@case.edu)
% Case Western Reserve University

function closed2DHexHistBatch()

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
staticLs = uicontrol('Style', 'text', 'Position', [675 500 50 30], 'String', 'Ls = ');
editLs = uicontrol('Style', 'edit', 'String', '.000005', 'Position', [725 500 50 31]);
staticUpWin = uicontrol('Style', 'text', 'Position', [675 450 50 30], 'String', 'UpWin Width = ');
editUpWin = uicontrol('Style', 'edit', 'String', '.000001', 'Position', [725 450 50 31]);
staticBotWin = uicontrol('Style', 'text', 'Position', [775 450 50 30], 'String', 'BotWin Width = ');
editBotWin = uicontrol('Style', 'edit', 'String', '.000001', 'Position', [825 450 50 31]);
% staticeW = uicontrol('Style', 'text', 'Position', [675 300 50 30], 'String', 'eW = ');
% editeW = uicontrol('Style', 'edit', 'String', '.01', 'Position', [725 300 50 31]);
% staticeL = uicontrol('Style', 'text', 'Position', [775 300 50 30], 'String', 'eL = ');
% editeL = uicontrol('Style', 'edit', 'String', '.01', 'Position', [825 300 50 31]);
buttonHist = uicontrol('Style', 'pushbutton', 'String', 'Histogram', 'Position',[675 350 200 40], 'Callback', @histogram);
% buttonCOT = uicontrol('Style', 'pushbutton', 'String', 'Change Over Time', 'Position',[675 250 200 40], 'Callback', @changeOverTime);
% buttonWCOT = uicontrol('Style', 'pushbutton', 'String', 'Window Change Over Time', 'Position',[675 150 200 40], 'Callback', @windowChangeOverTime);
% staticetw = uicontrol('Style', 'text', 'Position', [675 200 50 30], 'String', 'etw = ');
% editetw = uicontrol('Style', 'edit', 'String', '.000005', 'Position', [725 200 50 31]);
% staticebw = uicontrol('Style', 'text', 'Position', [775 200 50 30], 'String', 'ebw = ');
% editebw = uicontrol('Style', 'edit', 'String', '.000005', 'Position', [825 200 50 31]);
% buttontau = uicontrol('Style', 'pushbutton', 'String', 'MET as function of Tau', 'Position',[675 50 200 40], 'Callback', @tauChange);
% staticinitTau = uicontrol('Style', 'text', 'Position', [675 100 50 30], 'String', 'init Tau = ');
% editinittau = uicontrol('Style', 'edit', 'String', '100', 'Position', [725 100 50 31]);
% staticfintau = uicontrol('Style', 'text', 'Position', [775 100 50 30], 'String', 'fin Tau = ');
% editfintau = uicontrol('Style', 'edit', 'String', '.001', 'Position', [825 100 50 31]);

set(f3, 'Visible', 'on');

title('Particle Track of a Single Simulated Particle');
xlabel('Average Time ');
ylabel('Y Position');


    function histogram(source, eventdata)
        tic
        
        batch = [100 250 500 750 1000 2500 5000 7500 10000 25000 50000 75000 100000 250000 500000 750000 1000000];
            
        for a = 1:length(batch) 
        
        hold off;
        
        N = str2num(get(editN, 'String'));
        Ls = str2num(get(editLs, 'String'));
        wTop = str2num(get(editUpWin, 'String'));
        wBot = str2num(get(editBotWin, 'String'));
        
        placeHolder = 1;
        
        dimensions = 2;         % two dimensional simulation
        tau = .001 %str2num(get(editTau, 'String'));               % time interval in seconds
        
        R    = .145e-6;              % radius in meters
        eta  = 1.0e-3;              % viscosity of air in SI units (Pascal-seconds) at 293 K
        kB   = 1.38e-23;            % Boltzmann constant
        dT = .001;
        T    = 293;                 % Temperature in degrees Kelvin
        V = 4/3 * pi * R^3; %m^3
        p = 1.05 * 1000; % kg/cm
        
        m1 = p * V; %kg
        gamma = 6 * pi * eta * R;
        %tau = m/gamma
        
        R    = .145e-6;              % radius in meters eta  = 1.0e-3;   
        eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K
        kB   = 1.38e-23;            % Boltzmann constant
        T    = 293;                 % Temperature in degrees Kelvin
        
        D    = kB * T / (6 * pi * eta * R);
        k = sqrt(D * dimensions * tau);
        
        histArray = [];
        progressbar;
        for j = 1:batch(a)
    %    if get(rUntilExit, 'Value') == 1
           
            i = 1;
            x = 0;
            y = 0;
            
            while ~(x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= 3 ^ .5 * Ls)
                
                i = i + 1;
                
                dx = k * randn(1,1);
                dy = k * randn(1,1);
                
                x = [x, x(i - 1) + dx];
                y = [y, y(i - 1) + dy];
                
                if ~(x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= 3 ^ .5 * Ls)
                    if y(i) >= 3 ^ .5 * Ls / 2
                        if y(i) > 3 ^ .5 * Ls
                            y(i) = y(i) - 2 * (abs(y(i) - 3 ^ .5 * Ls));
                        end
                        if x(i) >= .5 * Ls
                            yc = (-3 ^ .5) * x(i) + 3 * (3 ^ .5) * Ls / 2;
                            
                            if y(i) > yc
                                m1 = (y(i) - y(i - 1))/(x(i) - x(i - 1));
                                xp = (3 * 3 ^ .5 * Ls / 2 + m1 * x(i - 1) - y(i - 1)) / ( m1 + 3 ^ .5);
                                yp = - 3 ^ .5 * xp + 3 * 3 ^ .5 * Ls / 2;
                                d = sqrt((x(i) - x(i - 1)) ^ 2 + (y(i) - y(i - 1)) ^ 2) - sqrt((xp - x(i - 1)) ^ 2 + (yp - y(i - 1)) ^ 2);
                                kappa =  2 * pi / 3 + pi / 6 + pi / 2 - atan(abs(m1));
                                
                                x(i) = xp + d * cos(kappa);
                                y(i) = yp + d * sin(kappa);
                                
                            end
                            
                        elseif x(i) <= -.5 * Ls
                            yc = 3 ^ .5 * x(i) + 3 * 3 ^ .5 * Ls / 2;

                            if y(i) > yc
                                m1 = (y(i) - y(i - 1))/(x(i) - x(i - 1));
                                xp = (3 * 3 ^ .5 * Ls / 2 + m1 * x(i - 1) - y(i - 1)) / ( m1 - 3 ^ .5);
                                yp = 3 ^ .5 * xp + 3 * 3 ^ .5 * Ls / 2;
                                d = sqrt((x(i) - x(i - 1)) ^ 2 + (y(i) - y(i - 1)) ^ 2) - sqrt((xp - x(i - 1)) ^ 2 + (yp - y(i - 1)) ^ 2);
                                kappa =  atan(abs(m1)) - pi / 3;

                                x(i) = xp + d * cos(kappa);
                                y(i) = yp + d * sin(kappa);
                                
                            end

                        end
                    elseif y(i) < 3 ^ .5 * Ls / 2
                        if y(i) < 0
                            y(i) = y(i) + 2 * abs(y(i));
                        end
                        if x(i) >= .5 * Ls
                            yc = 3 ^ .5 * x(i) - 3 ^ .5 * Ls / 2;
                            
                            if y(i) < yc
                                m1 = (y(i) - y(i - 1))/(x(i) - x(i - 1));
                                xp = (m1 * x(i - 1) - 3 ^ .5 * Ls / 2 - y(i - 1)) / ( m1 - 3 ^ .5);
                                yp = 3 ^ .5 * xp - 3 ^ .5 * Ls / 2;
                                d = sqrt((x(i) - x(i - 1)) ^ 2 + (y(i) - y(i - 1)) ^ 2) - sqrt((xp - x(i - 1)) ^ 2 + (yp - y(i - 1)) ^ 2);
                                kappa =  3 * pi / 2 - 2 * pi / 6 - pi / 2 + atan(abs(m1));
                                
                                x(i) = xp + d * cos(kappa);
                                y(i) = yp + d * sin(kappa);
                                
                            end

                        elseif x(i) <= -.5 * Ls
                            yc = -3 ^ .5 * x(i) -  3 ^ .5 * Ls / 2;

                            if y(i) < yc
                                m1 = (y(i) - y(i - 1))/(x(i) - x(i - 1));
                                xp = (m1 * x(i - 1) - 3 ^ .5 * Ls / 2 - y(i - 1)) / ( m1 + 3 ^ .5);
                                yp = -3 ^ .5 * xp - 3 ^ .5 * Ls / 2;
                                d = sqrt((x(i) - x(i - 1)) ^ 2 + (y(i) - y(i - 1)) ^ 2) - sqrt((xp - x(i - 1)) ^ 2 + (yp - y(i - 1)) ^ 2);
                                kappa =  3 * pi / 2 + 2 * pi / 6 + pi / 2 - atan(abs(m1));
                                
                                x(i) = xp + d * cos(kappa);
                                y(i) = yp + d * sin(kappa);
                            end
                        end
                    end
                end
            end
            histArray = [histArray, i];
            progressbar(a/length(batch), j/batch(a));
        end
        histArray = histArray * tau;
        bins = min(histArray):size(histArray)/10:max(histArray);
        figure(f3);
        hist(histArray, bins);
%         h = findobj(gca,'Type','patch');
%         set(h, 'EdgeColor','w')
        
        title(sprintf('Occurances versus Time for \nLarge Batch Diffusion (AVG = %f)', mean(histArray)));
        xlabel('Time');
        ylabel('Number of Occurances');
        
        filename = strcat('./2D Analysis/Hexagon/Eta = 1.000e-3/SA 1/Group 1 (25 square microns R1)/hist', num2str(batch(a)), '_1.fig');
        saveas(f3, filename, 'fig');
        
        end 
        
        toc
            
        
    end

    function changeOverTime(source, eventdata)
        tic
        N = str2num(get(editN, 'String'));
        D = str2num(get(editD, 'String'));
        sL = str2num(get(editL, 'String'));
        sW = str2num(get(editLs, 'String'));
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
        
        R    = .145e-6;              % radius in meters eta  = 1.562e-5;              % viscosity of water in SI units (Pascal-seconds) at 293 K kB   = 1.38e-23;            % Boltzmann constant T    = 293;                 % Temperature in degrees Kelvin
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
        W = str2num(get(editLs, 'String'));
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
        
        R    = .145e-6;              % radius in meters eta  = 1.562e-5;              % viscosity of water in SI units (Pascal-seconds) at 293 K kB   = 1.38e-23;            % Boltzmann constant T    = 293;                 % Temperature in degrees Kelvin
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
        W = str2num(get(editLs, 'String'));
        L = str2num(get(editL, 'String'));        
        
        step = .1e-6;
        histArray = [];
        avgArray = [];
        
        %randn('seed', seed);
        progressbar('Overall', 'Step Test');
        
        dimensions = 2;         % two dimensional simulation
        tau = .001;               % time interval in seconds
        time = tau * 1:N;       % create a time vector for plotting
        
        R    = .145e-6;              % radius in meters eta  = 1.562e-5;              % viscosity of water in SI units (Pascal-seconds) at 293 K kB   = 1.38e-23;            % Boltzmann constant T    = 293;                 % Temperature in degrees Kelvin
        eta  = 1.827e-5;              % viscosity of air in SI units (Pascal-seconds) at 293 K 
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