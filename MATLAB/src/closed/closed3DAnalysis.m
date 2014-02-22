% Adam Gleichsner (amg188@case.edu)
% Case Western Reserve University

function closed3DAnalysis()

x = 0;
y = 0;

seed = '1';

f3 = figure('Position', [1 1 900 800], 'Visible', 'off', 'MenuBar', 'none', 'Name', 'Closed 3D Analysis');
movegui(f3, 'center');
axhan3 = axes('Units', 'Pixels', 'Position', [150, 150 , 400, 400]);
staticN = uicontrol('Style', 'text', 'Position', [675 700 100 30], 'String', 'Steps(N) =');
editN = uicontrol('Style', 'edit', 'String', '1000', 'Position', [775 700 100 31]);
staticD = uicontrol('Style', 'text', 'Position', [675 650 100 30], 'String', 'Diffusion Coef(D) =');
editD = uicontrol('Style', 'edit', 'String', '.000045', 'Position', [775 650 100 31]);
%staticS = uicontrol('Style', 'text', 'Position', [675 600 100 30], 'String', 'Random Seed =');
%editS = uicontrol('Style', 'edit', 'String', seed, 'Position', [775 600 100 31]);
%buttonS = uicontrol('Style', 'pushbutton', 'String', 'Random Seed', 'Position',[675 400 200 40], 'Callback', @randomSeed);
staticB = uicontrol('Style', 'text', 'Position', [675 550 100 30], 'String', 'Batch Size = ');
editB = uicontrol('Style', 'edit', 'String', '100', 'Position', [775 550 100 31]);
staticW = uicontrol('Style', 'text', 'Position', [675 500 33 30], 'String', 'W = ');
editW = uicontrol('Style', 'edit', 'String', '.1', 'Position', [708 500 33 31]);
staticL = uicontrol('Style', 'text', 'Position', [742 500 33 30], 'String', 'L = ');
editL = uicontrol('Style', 'edit', 'String', '.1', 'Position', [770 500 33 31]);
staticH = uicontrol('Style', 'text', 'Position', [803 500 33 30], 'String', 'H = ');
editH = uicontrol('Style', 'edit', 'String', '.1', 'Position', [836 500 33 31]);
staticUpWin = uicontrol('Style', 'text', 'Position', [675 450 50 30], 'String', 'Up w = ');
editUpWin = uicontrol('Style', 'edit', 'String', '.01', 'Position', [725 450 50 31]);
staticBotWin = uicontrol('Style', 'text', 'Position', [775 450 50 30], 'String', 'Bot w = ');
editBotWin = uicontrol('Style', 'edit', 'String', '.01', 'Position', [825 450 50 31]);
staticUpL = uicontrol('Style', 'text', 'Position', [675 400 50 30], 'String', 'Up l = ');
editUpL = uicontrol('Style', 'edit', 'String', '.01', 'Position', [725 400 50 31]);
staticBotL = uicontrol('Style', 'text', 'Position', [775 400 50 30], 'String', 'Bot l = ');
editBotL = uicontrol('Style', 'edit', 'String', '.01', 'Position', [825 400 50 31]);
staticeW = uicontrol('Style', 'text', 'Position', [675 300 50 30], 'String', 'eW = ');
editeW = uicontrol('Style', 'edit', 'String', '.01', 'Position', [725 300 50 31]);
staticeL = uicontrol('Style', 'text', 'Position', [775 300 50 30], 'String', 'eL = ');
editeL = uicontrol('Style', 'edit', 'String', '.01', 'Position', [825 300 50 31]);
buttonHist = uicontrol('Style', 'pushbutton', 'String', 'Histogram', 'Position',[675 350 200 40], 'Callback', @histogram);
buttonCOT = uicontrol('Style', 'pushbutton', 'String', 'Change Over Time', 'Position',[675 250 200 40], 'Callback', @changeOverTime);
buttonWCOT = uicontrol('Style', 'pushbutton', 'String', 'Window Change Over Time', 'Position',[675 150 200 40], 'Callback', @windowChangeOverTime);
staticetw = uicontrol('Style', 'text', 'Position', [675 200 50 30], 'String', 'etw = ');
editetw = uicontrol('Style', 'edit', 'String', '.01', 'Position', [725 200 50 31]);
staticebw = uicontrol('Style', 'text', 'Position', [775 200 50 30], 'String', 'ebw = ');
editebw = uicontrol('Style', 'edit', 'String', '.01', 'Position', [825 200 50 31]);
%h = uibuttongroup('visible','off','Position',[0 0 0.0001 0.1], 'SelectionChangeFcn', @switchButton);
%rUponExit = uicontrol('Style','Radio','String','Stop plot upon exit', 'pos',[675 200 200 30], 'parent',h);
%rUntilExit = uicontrol('Style','Radio','String','Plot until exit', 'pos',[675 175 200 30], 'parent',h);


set(f3, 'Visible', 'on');
%set(h,'SelectedObject',[]);
%set(h,'Visible','on');

title('Particle Track of a Single Simulated Particle');
xlabel('Average Time ');
ylabel('Y Position');


    function histogram(source, eventdata)
        
        N = str2num(get(editN, 'String'));
        D = str2num(get(editD, 'String'));
        L = str2num(get(editL, 'String'));
        W = str2num(get(editW, 'String'));
        H = str2num(get(editH, 'String'));
        wTop = str2num(get(editUpWin, 'String'));
        wBot = str2num(get(editBotWin, 'String'));
        lTop = str2num(get(editUpL, 'String'));
        lBot = str2num(get(editBotL, 'String'));
        randArray = randn(N, 2);
        
        dimensions = 3;         % two dimensional simulation
        tau = .1;               % time interval in seconds
        time = tau * 1:N;       % create a time vector for plotting
        
        histArray = [];
        
        N = N * 100;
        k = sqrt(D * dimensions * tau);
        
        for i = 1:str2num(get(editB, 'String'))
            
            dx = [L/2];
            dy = [0];
            dz = [L/2];
            dx2 = k * randn(N - 1,1);
            dy2 = k * randn(N - 1,1);
            dz2 = k * randn(N - 1,1);
            
            for i = 2:N
                dx(i) = dx2(i - 1);
                dy(i) = dy2(i - 1);
                dz(i) = dz2(i - 1);
            end
            
            x = cumsum(dx);
            y = cumsum(dy);
            z = cumsum(dz);
            
            j = 1;
            
            while(j < N)
                if x(j) > W/2 - wTop/2 && x(j) < W/2 + wTop/2 && y(j) >= L && z(j) > H/2 - lTop/2 && z(j) < H/2 + lTop/2
                    histArray = [histArray, j]
                    j
                    j = N - 1;
                else
                    if x(j) < -W/2
                        xDif = 2 * (abs(x(j)) - W/2);
                        for k = j:N
                            x(k) = x(k) + xDif;
                        end
                    end
                    
                    if x(j) > W/2
                        xDif = 2 * (abs(x(j)) - W/2);
                        for k = j:N
                            x(k) = x(k) - xDif;
                        end
                    end
                    
                    if y(j) > L
                        yDif = 2 * (abs(y(j)) - L);
                        for k = j:N
                            y(k) = y(k) - yDif;
                        end
                    end
                    
                    if y(j) < 0
                        yDif = 2 * abs(y(j));
                        for k = j:N
                            y(k) = y(k) + yDif;
                        end
                    end
                    
                    if z(j) > H
                        zDif = 2 * (abs(z(j)) - H);
                        for k = j:N
                            z(k) = z(k) - zDif;
                        end
                    end
                    
                    if z(j) < 0
                        zDif = 2 * (abs(z(j)));
                        for k = j:N
                            z(k) = z(k) + zDif;
                        end
                    end
                    fprintf('ping %d\n', j);
                end
                j = j + 1;
            end
        end
        
        
        
        histArray = histArray * k
        bins = [0:1:max(histArray)];
        hist(histArray, bins);
        
        title(sprintf('Occurances versus Time for \nLarge Batch Diffusion (AVG = %f)', mean(histArray)));
        xlabel('Time');
        ylabel('Number of Occurances');
        
    end

    function changeOverTime(source, eventdata)
        
        N = str2num(get(editN, 'String'));
        D = str2num(get(editD, 'String'));
        %seed = str2num(get(editS, 'String'));
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
        
        %randn('seed', seed);
        
        dimensions = 2;         % two dimensional simulation
        tau = .1;               % time interval in seconds
        time = tau * 1:N;       % create a time vector for plotting
        
        
        N = N * 10;
        k = sqrt(D * dimensions * tau);
        
        L = sL;
        W = sW;
        
        for yy = sL:step:eL
            
            L = yy;
            W = yy;
            
            for z = 1:str2num(get(editB, 'String'))
                
                N = str2num(get(editN, 'String'));
                N = N * 10;
                dx = [0];
                dy = [0];
                dx2 = k * randn(N - 1,1);
                dy2 = k * randn(N - 1,1);
                
                for i = 2:N
                    dx(i) = dx2(i - 1);
                    dy(i) = dy2(i - 1);
                end
                
                x = cumsum(dx);
                y = cumsum(dy);
                
                i = 1;
                
                while(i < N)
                    if x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= L
                        N = i;
                        histArray = [histArray, N];
                    else
                        if x(i) < -W/2
                            xDif = 2 * (abs(x(i)) - W/2);
                            for j = i:N
                                x(j) = x(j) + xDif;
                            end
                        end
                        
                        if x(i) > W/2
                            xDif = 2 * (abs(x(i)) - W/2);
                            for j = i:N
                                x(j) = x(j) - xDif;
                            end
                        end
                        
                        if y(i) > L
                            yDif = 2 * (abs(y(i)) - L);
                            for j = i:N
                                y(j) = y(j) - yDif;
                            end
                        end
                        
                        if y(i) < 0
                            yDif = 2 * abs(y(i));
                            for j = i:N
                                y(j) = y(j) + yDif;
                            end
                        end
                    end
                    i = i + 1;
                end
            end
            
            avgArray = [avgArray, mean(histArray)];
            histArray = [];
            yy
        end
        avgArray = avgArray * k
        plot([sW:step:eW], avgArray);
        
        title('Average Time of Diffusion versus Change in Void Dimension');
        xlabel('Current Dimension Length (W = , L =)');
        ylabel('Average Time per Batch');
    end

    function windowChangeOverTime(source, eventdata)
        
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
        randArray = randn(N, 2);
        
        
        step = .0001;
        histArray = [];
        avgArray = [];
        
        %randn('seed', seed);
        
        dimensions = 2;         % two dimensional simulation
        tau = .1;               % time interval in seconds
        time = tau * 1:N;       % create a time vector for plotting
        
        
        N = N * 10;
        k = sqrt(D * dimensions * tau);
        
        for yy = sbW:step:ebW
            
            wTop = yy;
            wBot = yy;
            
            for z = 1:str2num(get(editB, 'String'))
                
                N = str2num(get(editN, 'String'));
                N = N * 10;
                dx = [0];
                dy = [0];
                dx2 = k * randn(N - 1,1);
                dy2 = k * randn(N - 1,1);
                
                for i = 2:N
                    dx(i) = dx2(i - 1);
                    dy(i) = dy2(i - 1);
                end
                
                x = cumsum(dx);
                y = cumsum(dy);
                
                i = 1;
                
                while(i < N)
                    if x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= L
                        N = i;
                        histArray = [histArray, N];
                    else
                        if x(i) < -W/2
                            xDif = 2 * (abs(x(i)) - W/2);
                            for j = i:N
                                x(j) = x(j) + xDif;
                            end
                        end
                        
                        if x(i) > W/2
                            xDif = 2 * (abs(x(i)) - W/2);
                            for j = i:N
                                x(j) = x(j) - xDif;
                            end
                        end
                        
                        if y(i) > L
                            yDif = 2 * (abs(y(i)) - L);
                            for j = i:N
                                y(j) = y(j) - yDif;
                            end
                        end
                        
                        if y(i) < 0
                            yDif = 2 * abs(y(i));
                            for j = i:N
                                y(j) = y(j) + yDif;
                            end
                        end
                    end
                    i = i + 1;
                end
            end
            
            avgArray = [avgArray, mean(histArray)];
            histArray = [];
            yy
        end
        avgArray = avgArray * k
        plot([sbW:step:ebW], avgArray);
        
        title('Average Time of Diffusion versus Change in Void Dimension');
        xlabel('Current Dimension Length (W = , L =)');
        ylabel('Average Time per Batch');
    end


end