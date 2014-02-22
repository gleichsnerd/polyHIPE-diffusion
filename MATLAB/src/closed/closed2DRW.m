% Adam Gleichsner (amg188@case.edu)
% Case Western Reserve University

function closed2DRW()

x = 0;
y = 0;

seed = '1.000';
N = 1000;


f2 = figure('Position', [1 1 900 800], 'Toolbar', 'figure', 'Visible', 'off', 'MenuBar', 'none', 'Name', 'Random Walk in Closed 2D');
movegui(f2, 'center')
axhan2 = axes('Units', 'Pixels', 'Position', [50, 150 , 600, 600]);
staticN = uicontrol('Style', 'text', 'Position', [675 700 100 30], 'String', 'Steps(N) =');
editN = uicontrol('Style', 'edit', 'String', '1000', 'Position', [775 700 100 31]);
staticD = uicontrol('Style', 'text', 'Position', [675 650 100 30], 'String', 'Diffusion Coef(D) =');
editD = uicontrol('Style', 'edit', 'String', '.000045', 'Position', [775 650 100 31]);
staticS = uicontrol('Style', 'text', 'Position', [675 600 100 30], 'String', 'Random Seed =');
editS = uicontrol('Style', 'edit', 'String', seed, 'Position', [775 600 100 31]);
buttonS = uicontrol('Style', 'pushbutton', 'String', 'Random Seed', 'Position',[675 400 200 40], 'Callback', @randomSeed);
staticW = uicontrol('Style', 'text', 'Position', [675 500 50 30], 'String', 'W = ');
editW = uicontrol('Style', 'edit', 'String', '.00001', 'Position', [725 500 50 31]);
staticL = uicontrol('Style', 'text', 'Position', [775 500 50 30], 'String', 'L = ');
editL = uicontrol('Style', 'edit', 'String', '.00001', 'Position', [825 500 50 31]);
staticUpWin = uicontrol('Style', 'text', 'Position', [675 450 50 30], 'String', 'UpWin Width = ');
editUpWin = uicontrol('Style', 'edit', 'String', '.000005', 'Position', [725 450 50 31]);
staticBotWin = uicontrol('Style', 'text', 'Position', [775 450 50 30], 'String', 'BotWin Width = ');
editBotWin = uicontrol('Style', 'edit', 'String', '.000005', 'Position', [825 450 50 31]);
buttonUpdate = uicontrol('Style', 'pushbutton', 'String', 'Update', 'Position',[675 350 200 40], 'Callback', @updateGraph);
buttonDisp = uicontrol('Style', 'pushbutton', 'String', 'Plot Displacement', 'Position',[675 250 200 40], 'Callback', @plotDisplacement);
h = uibuttongroup('visible','off','Position',[0 0 0.0001 0.1], 'SelectionChangeFcn', @switchButton);
rUponExit = uicontrol('Style','Radio','String','Stop plot upon exit', 'pos',[675 200 200 30], 'parent',h);
rUntilExit = uicontrol('Style','Radio','String','Plot until exit', 'Selected', 'on','pos',[675 175 200 30], 'parent',h);
staticTau = uicontrol('Style', 'text', 'Position', [675 100 50 30], 'String', 'dT = ');
editTau = uicontrol('Style', 'edit', 'String', '.1', 'Position', [725 100 50 31]);

set(f2, 'Visible', 'on');
set(h,'SelectedObject',[]);
set(h,'Visible','on');

title('Particle Track of a Single Simulated Particle');
xlabel('X Position');
ylabel('Y Position');

    function plotDisplacement(source, eventdata)
        
        N = str2num(get(editN, 'String'));
        D = str2num(get(editD, 'String'));
        seed = str2num(get(editS, 'String'));
        L = str2num(get(editL, 'String'));
        W = str2num(get(editW, 'String'));
        wTop = str2num(get(editUpWin, 'String'));
        wBot = str2num(get(editBotWin, 'String'));
        randArray = randn(N, 2);
 
        
        randn('seed', seed);
        
        dimensions = 2;         % two dimensional simulation
        tau = str2num(get(editTau, 'String'));               % time interval in seconds
        time = tau * 1:N;       % create a time vector for plotting
        
        R    = .145e-6;              % radius in meters
        eta  = 1.0e-3;              % viscosity of air in SI units (Pascal-seconds) at 293 K
        kB   = 1.38e-23;            % Boltzmann constant
        T    = 293;                 % Temperature in degrees Kelvin
        V = 4/3 * pi * R^3; %m^3
        p = 1.05 * 1000; % kg/cm
        
        m = p * V; %kg
        gamma = 6 * pi * eta * R;
        %tau = m/gamma;
        
        R    = .145e-6;              % radius in meters eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K kB   = 1.38e-23;            % Boltzmann constant T    = 293;                 % Temperature in degrees Kelvin
        eta  = 1.827e-5;              % viscosity of air in SI units (Pascal-seconds) at 293 K 
        kB   = 1.38e-23;            % Boltzmann constant
        T    = 293;                 % Temperature in degrees Kelvin
        
        D    = kB * T / (3 * pi * eta * d);
        
        k = sqrt(D * dimensions * tau);
        
        if get(rUntilExit, 'Value') == 1
            
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
            
            set(editN, 'String', num2str(i));
        
        elseif get(rUponExit, 'Value') == 1
            
            i = 1;
            x = 0;
            y = 0;
            
            while ~(x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= L) && i < N
                
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
            
            set(editN, 'String', num2str(i));
            
        else
            i = 1;
            x = 0;
            y = 0;
            
            while i < N
                
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
        end
        
        displacement = sqrt(x .^ 2 + y .^ 2);
        hold off;
        plot(displacement);
        axis([0 length(displacement) 0 max(displacement) + .05 * max(displacement)]);
        
        ylabel('Displacement');
        xlabel('Step');
        title('Displacement versus Steps of Simulated Particle Diffusion');
        
    end

    function switchButton(source, eventdata)
        
        event = get(eventdata.NewValue,'String');
        if strcmp(event, 'Stop plot upon exit')
            set(rUntilExit, 'Selected', 'off');
        elseif strcmp(event, 'Plot until exit')
            set(rUponExit, 'Selected', 'off');
        end
    end


    function updateGraph(source, eventdata)
        
        N = str2num(get(editN, 'String'));
        D = str2num(get(editD, 'String'));
        seed = str2num(get(editS, 'String'));
        L = str2num(get(editL, 'String'));
        W = str2num(get(editW, 'String'));
        wTop = str2num(get(editUpWin, 'String'));
        wBot = str2num(get(editBotWin, 'String'));
        randArray = randn(N, 2);
 
        hold off;
        
        placeHolder = 1;
        
        randn('seed', seed);
        
        dimensions = 2;         % two dimensional simulation
        tau = .001;%str2num(get(editTau, 'String'));               % time interval in seconds
        time = tau * 1:N;       % create a time vector for plotting
        
        R    = .145e-6;              % radius in meters
        eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K
        kB   = 1.38e-23;            % Boltzmann constant
        dT = .001;
        T    = 293;                 % Temperature in degrees Kelvin
        V = 4/3 * pi * R^3; %m^3
        p = 1.05 * 1000; % kg/cm
        
        m = p * V; %kg
        gamma = 6 * pi * eta * R;
        %tau = m/gamma
        
        R    = .145e-6;              % radius in meters eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K kB   = 1.38e-23;            % Boltzmann constant T    = 293;                 % Temperature in degrees Kelvin
        eta  = 1.00e-3;              % viscosity of air in SI units (Pascal-seconds) at 293 K
        kB   = 1.38e-23;            % Boltzmann constant
        T    = 293;                 % Temperature in degrees Kelvin
        
        D    = kB * T / (6 * pi * eta * R);
        
        k = sqrt(D * dimensions * tau);
        
        if get(rUntilExit, 'Value') == 1
            
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
                        slope = (y(i) - y(i-1))/(x(i) - x(i - 1));
                        b = y(i) - slope * x(i);
                        newY = slope * -W/2 + b;
                        
                        plot(x(placeHolder:i-1), y(placeHolder: i - 1));
                        hold on;
                        plot([x(i - 1), -W/2], [y(i - 1), newY]);
                        
                        x(i) = x(i) + 2 * (abs(x(i)) - W/2);
                        
                        plot([-W/2, x(i)], [newY, y(i)])
                        placeHolder = i;
                    elseif x(i) > W/2
                        slope = (y(i) - y(i-1))/(x(i) - x(i - 1));
                        b = y(i) - slope * x(i);
                        newY = slope * W/2 + b;
                        
                        plot(x(placeHolder:i-1), y(placeHolder: i - 1));
                        hold on;
                        plot([x(i - 1), W/2], [y(i - 1), newY]);
                        
                        x(i) = x(i) - 2 * (abs(x(i)) - W/2);
                        
                        plot([W/2, x(i)], [newY, y(i)])
                        placeHolder = i;
                    end
                    
                    if y(i) < 0
                        slope = (y(i) - y(i-1))/(x(i) - x(i - 1));
                        b = y(i) - slope * x(i);
                        newX = (0 - b) / slope;
                        
                        plot(x(placeHolder:i-1), y(placeHolder: i - 1));
                        hold on;
                        plot([x(i - 1), newX], [y(i - 1), 0]);
                        
                        y(i) = y(i) + 2 * abs(y(i));
                        
                        plot([newX, x(i)], [0, y(i)])
                        placeHolder = i;
                    elseif y(i) > L
                        slope = (y(i) - y(i-1))/(x(i) - x(i - 1));
                        b = y(i) - slope * x(i);
                        newX = (L - b) / slope;
                        
                        plot(x(placeHolder:i-1), y(placeHolder: i - 1));
                        hold on;
                        plot([x(i - 1), newX], [y(i - 1), L]);
                        
                        y(i) = y(i) - 2 * (abs(y(i) - L));
                        
                        plot([newX, x(i)], [L, y(i)])
                        placeHolder = i;
                    end
                else
                    slope = (y(i) - y(i-1))/(x(i) - x(i - 1));
                    b = y(i) - slope * x(i);
                    newX = (L - b) / slope;
                    
                    if~(newX > -wTop/2 && newX < wTop/2)
                        plot(x(placeHolder:i-1), y(placeHolder: i - 1));
                        hold on;
                        plot([x(i - 1), newX], [y(i - 1), L]);
                        
                        y(i) = y(i) - 2 * (abs(y(i) - L));
                        
                        plot([newX, x(i)], [L, y(i)])
                        placeHolder = i;
                    end
                end
            end
            
            set(editN, 'String', num2str(i));
        
        elseif get(rUponExit, 'Value') == 1
            
            i = 1;
            x = 0;
            y = 0;
            
            while ~(x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= L) && i < N
                
                i = i + 1;
                
                dx = k * randn(1,1);
                dy = k * randn(1,1);
                
                x = [x, x(i - 1) + dx];
                y = [y, y(i - 1) + dy];
                
                if ~(x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= L)
                    if x(i) < -W/2
                        slope = (y(i) - y(i-1))/(x(i) - x(i - 1));
                        b = y(i) - slope * x(i);
                        newY = slope * -W/2 + b;
                        
                        plot(x(placeHolder:i-1), y(placeHolder: i - 1));
                        hold on;
                        plot([x(i - 1), -W/2], [y(i - 1), newY]);
                        
                        x(i) = x(i) + 2 * (abs(x(i)) - W/2);
                        
                        plot([-W/2, x(i)], [newY, y(i)])
                        placeHolder = i;
                    elseif x(i) > W/2
                        slope = (y(i) - y(i-1))/(x(i) - x(i - 1));
                        b = y(i) - slope * x(i);
                        newY = slope * W/2 + b;
                        
                        plot(x(placeHolder:i-1), y(placeHolder: i - 1));
                        hold on;
                        plot([x(i - 1), W/2], [y(i - 1), newY]);
                        
                        x(i) = x(i) - 2 * (abs(x(i)) - W/2);
                        
                        plot([W/2, x(i)], [newY, y(i)])
                        placeHolder = i;
                    end
                    
                    if y(i) < 0
                        slope = (y(i) - y(i-1))/(x(i) - x(i - 1));
                        b = y(i) - slope * x(i);
                        newX = (0 - b) / slope;
                        
                        plot(x(placeHolder:i-1), y(placeHolder: i - 1));
                        hold on;
                        plot([x(i - 1), newX], [y(i - 1), 0]);
                        
                        y(i) = y(i) + 2 * abs(y(i));
                        
                        plot([newX, x(i)], [0, y(i)])
                        placeHolder = i;
                    elseif y(i) > L
                        slope = (y(i) - y(i-1))/(x(i) - x(i - 1));
                        b = y(i) - slope * x(i);
                        newX = (L - b) / slope;
                        
                        plot(x(placeHolder:i-1), y(placeHolder: i - 1));
                        hold on;
                        plot([x(i - 1), newX], [y(i - 1), L]);
                        
                        y(i) = y(i) - 2 * (abs(y(i) - L));
                        
                        plot([newX, x(i)], [L, y(i)])
                        placeHolder = i;
                    end
                else
                    slope = (y(i) - y(i-1))/(x(i) - x(i - 1));
                    b = y(i) - slope * x(i);
                    newX = (L - b) / slope;
                    
                    if~(newX > -wTop/2 && newX < wTop/2)
                        plot(x(placeHolder:i-1), y(placeHolder: i - 1));
                        hold on;
                        plot([x(i - 1), newX], [y(i - 1), L]);
                        
                        y(i) = y(i) - 2 * (abs(y(i) - L));
                        
                        plot([newX, x(i)], [L, y(i)])
                        placeHolder = i;
                    end
                end
            end
            
            set(editN, 'String', num2str(i));
            
        end
       
        
        plot(-wTop/2, L, 'ro', wTop/2, L, 'ro', -wBot/2, 0, 'ro', wBot/2, 0, 'ro');
        hold on;
        plot([-W/2,-W/2,W/2,W/2,-W/2], [0, L, L, 0, 0], 'r');
        hold on;
        plot(x(placeHolder:length(x)), y(placeHolder:length(y)));
        hold on;
        axis([-W/2 W/2 0 L]);
        
        title('Particle Track of a Single Simulated Particle');
        xlabel('X Position');
        ylabel('Y Position');
        
        
    end

    function randomSeed(source, eventdata)
        
        seed = num2str(rand(1,1)*100);
        set(editS, 'String', seed);
        
    end

end