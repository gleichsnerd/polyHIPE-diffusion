% Adam Gleichsner (amg188@case.edu)
% Case Western Reserve University

function closed3DRW()

x = 0;
y = 0;
z = 0;

seed = '1.000'
N = 1000;
dimensions = 3;         % two dimensional simulation
tau = .001;               % time interval in seconds
time = tau * 1:N;       % create a time vector for plotting

R    = .145e-6;              % radius in meters
eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K
kB   = 1.38e-23;            % Boltzmann constant
T    = 293;                 % Temperature in degrees Kelvin

D    = kB * T / (6 * pi * eta * R);


seedFile = strcat(num2str(seed), '.txt');

f2 = figure('Position', [1 1 900 800], 'Toolbar', 'figure', 'Visible', 'off', 'MenuBar', 'none', 'Name', 'Random Walk in Closed 3D');
movegui(f2, 'center')
axhan2 = axes('Units', 'Pixels', 'Position', [50, 150 , 600, 600]);
staticN = uicontrol('Style', 'text', 'Position', [675 700 100 30], 'String', 'Steps(N) =');
editN = uicontrol('Style', 'edit', 'String', '1000', 'Position', [775 700 100 31]);
staticD = uicontrol('Style', 'text', 'Position', [675 650 100 30], 'String', 'Diffusion Coef(D) =');
editD = uicontrol('Style', 'edit', 'String', num2str(D), 'Position', [775 650 100 31]);
staticS = uicontrol('Style', 'text', 'Position', [675 600 100 30], 'String', 'Random Seed =');
editS = uicontrol('Style', 'edit', 'String', seed, 'Position', [775 600 100 31]);
buttonS = uicontrol('Style', 'pushbutton', 'String', 'Random Seed', 'Position',[675 350 200 40], 'Callback', @randomSeed);
staticW = uicontrol('Style', 'text', 'Position', [675 500 33 30], 'String', 'W = ');
editW = uicontrol('Style', 'edit', 'String', '.00001', 'Position', [708 500 33 31]);
staticL = uicontrol('Style', 'text', 'Position', [742 500 33 30], 'String', 'L = ');
editL = uicontrol('Style', 'edit', 'String', '.00001', 'Position', [770 500 33 31]);
staticH = uicontrol('Style', 'text', 'Position', [803 500 33 30], 'String', 'H = ');
editH = uicontrol('Style', 'edit', 'String', '.00001', 'Position', [836 500 33 31]);
staticUpWin = uicontrol('Style', 'text', 'Position', [675 450 50 30], 'String', 'Up w = ');
editUpWin = uicontrol('Style', 'edit', 'String', '.000005', 'Position', [725 450 50 31]);
staticBotWin = uicontrol('Style', 'text', 'Position', [775 450 50 30], 'String', 'Bot w = ');
editBotWin = uicontrol('Style', 'edit', 'String', '.000005', 'Position', [825 450 50 31]);
staticUpL = uicontrol('Style', 'text', 'Position', [675 400 50 30], 'String', 'Up l = ');
editUpL = uicontrol('Style', 'edit', 'String', '.000005', 'Position', [725 400 50 31]);
staticBotL = uicontrol('Style', 'text', 'Position', [775 400 50 30], 'String', 'Bot l = ');
editBotL = uicontrol('Style', 'edit', 'String', '.000005', 'Position', [825 400 50 31]);
buttonUpdate = uicontrol('Style', 'pushbutton', 'String', 'Update', 'Position',[675 300 200 40], 'Callback', @updateGraph);
buttonDisp = uicontrol('Style', 'pushbutton', 'String', 'Plot Displacement', 'Position',[675 200 200 40], 'Callback', @plotDisplacement);
h = uibuttongroup('visible','off','Position',[0 0 0.0001 0.1], 'SelectionChangeFcn', @switchButton);
rUponExit = uicontrol('Style','Radio','String','Stop plot upon exit', 'pos',[675 150 200 30], 'parent',h);
rUntilExit = uicontrol('Style','Radio','String','Plot until exit', 'pos',[675 125 200 30], 'parent',h);


set(f2, 'Visible', 'on');
set(h,'SelectedObject',[]);
set(h,'Visible','on');

title('Particle Track of a Single Simulated Particle');
xlabel('X Position');
ylabel('Y Position');

    function plotDisplacement(source, eventdata)
        
        N = str2num(get(editN, 'String'));
        seed = str2num(get(editS, 'String'));
        L = str2num(get(editL, 'String'));
        W = str2num(get(editW, 'String'));
        H = str2num(get(editH, 'String'));
        wTop = str2num(get(editUpWin, 'String'));
        wBot = str2num(get(editBotWin, 'String'));
        lTop = str2num(get(editUpL, 'String'));
        lBot = str2num(get(editBotL, 'String'));
        
        randn('seed', seed);
        
        dimensions = 3;         % two dimensional simulation
        tau = .1;               % time interval in seconds
        time = tau * 1:N;       % create a time vector for plotting
        
        R    = .145e-6;              % radius in meters eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K kB   = 1.38e-23;            % Boltzmann constant T    = 293;                 % Temperature in degrees Kelvin
        eta  = 1.827e-5;              % viscosity of air in SI units (Pascal-seconds) at 293 K 
        kB   = 1.38e-23;            % Boltzmann constant
        T    = 293;                 % Temperature in degrees Kelvin
        
        D    = kB * T / (3 * pi * eta * d);
        
        k = sqrt(D * dimensions * tau);
        
        x = W/2;
        y = 0;
        z = L/2;
        
        if get(rUntilExit, 'Value') == 1
            
            i = 1;
            
            while ~(x(i) > W/2 - wTop/2 && x(i) < W/2 + wTop/2 && y(i) >= L && z(i) > H/2 - lTop/2 && z(i) < H/2 + lTop/2)
                
                i = i + 1;
                
                dx = k * randn(1,1);
                dy = k * randn(1,1);
                dz = k * randn(1,1);
                
                x = [x, x(i - 1) + dx];
                y = [y, y(i - 1) + dy];
                z = [z, z(i - 1) + dz];
                
                if ~(x(i) > W/2 - wTop/2 && x(i) < W/2 + wTop/2 && y(i) >= L && z(i) > H/2 - lTop/2 && z(i) < H/2 + lTop/2)
                    if x(i) < 0
                        x(i) = x(i) + 2 * abs(x(i));
                    elseif x(i) > W
                        x(i) = x(i) - 2 * (abs(x(i)) - W);
                    end
                    
                    if y(i) < 0
                        y(i) = y(i) + 2 * abs(y(i));
                    elseif y(i) > L
                        y(i) = y(i) - 2 * (abs(y(i) - L));
                    end
                    
                    if z(i) < 0
                        z(i) = z(i) + 2 * abs(z(i));
                    elseif z(i) > H
                        z(i) = z(i) - 2 * (abs(z(i)) - H);
                    end
                end
            end
            
            set(editN, 'String', num2str(i));
            
        elseif get(rUponExit, 'Value') == 1
            
            i = 1;
            
            while ~(x(i) > W/2 - wTop/2 && x(i) < W/2 + wTop/2 && y(i) >= L && z(i) > H/2 - lTop/2 && z(i) < H/2 + lTop/2) && i < N
                
                i = i + 1;
                
                dx = k * randn(1,1);
                dy = k * randn(1,1);
                dz = k * randn(1,1);
                
                x = [x, x(i - 1) + dx];
                y = [y, y(i - 1) + dy];
                z = [z, z(i - 1) + dz];
                
                if ~(x(i) > W/2 - wTop/2 && x(i) < W/2 + wTop/2 && y(i) >= L && z(i) > H/2 - lTop/2 && z(i) < H/2 + lTop/2)
                    if x(i) < 0
                        x(i) = x(i) + 2 * abs(x(i));
                    elseif x(i) > W
                        x(i) = x(i) - 2 * (abs(x(i)) - W);
                    end
                    
                    if y(i) < 0
                        y(i) = y(i) + 2 * abs(y(i));
                    elseif y(i) > L
                        y(i) = y(i) - 2 * (abs(y(i) - L));
                    end
                    
                    if z(i) < 0
                        z(i) = z(i) + 2 * abs(z(i));
                    elseif z(i) > H
                        z(i) = z(i) - 2 * (abs(z(i)) - H);
                    end
                end
            end
            
            set(editN, 'String', num2str(i));
            
        else
            i = 1;
            
            while i < N
                
                i = i + 1;
                
                dx = k * randn(1,1);
                dy = k * randn(1,1);
                dz = k * randn(1,1);
                
                x = [x, x(i - 1) + dx];
                y = [y, y(i - 1) + dy];
                z = [z, z(i - 1) + dz];
                
                if ~(x(i) > W/2 - wTop/2 && x(i) < W/2 + wTop/2 && y(i) >= L && z(i) > H/2 - lTop/2 && z(i) < H/2 + lTop/2)
                    if x(i) < 0
                        x(i) = x(i) + 2 * abs(x(i));
                    elseif x(i) > W
                        x(i) = x(i) - 2 * (abs(x(i)) - W);
                    end
                    
                    if y(i) < 0
                        y(i) = y(i) + 2 * abs(y(i));
                    elseif y(i) > L
                        y(i) = y(i) - 2 * (abs(y(i) - L));
                    end
                    
                    if z(i) < 0
                        z(i) = z(i) + 2 * abs(z(i));
                    elseif z(i) > H
                        z(i) = z(i) - 2 * (abs(z(i)) - H);
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
        seed = str2num(get(editS, 'String'));
        L = str2num(get(editL, 'String'));
        W = str2num(get(editW, 'String'));
        H = str2num(get(editH, 'String'));
        wTop = str2num(get(editUpWin, 'String'));
        wBot = str2num(get(editBotWin, 'String'));
        lTop = str2num(get(editUpL, 'String'));
        lBot = str2num(get(editBotL, 'String'));
        
        randn('seed', seed);
        
        dimensions = 3;         % two dimensional simulation
        tau = .001;               % time interval in seconds
        time = tau * 1:N;       % create a time vector for plotting
        
        R    = .145e-6;              % radius in meters eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K kB   = 1.38e-23;            % Boltzmann constant T    = 293;                 % Temperature in degrees Kelvin
        eta  = 1.827e-5;              % viscosity of air in SI units (Pascal-seconds) at 293 K 
        kB   = 1.38e-23;            % Boltzmann constant
        T    = 293;                 % Temperature in degrees Kelvin
        
        D    = kB * T / (6 * pi * eta * R);
        
        k = sqrt(D * dimensions * tau);
        
        x = W/2;
        y = 0;
        z = L/2;
        
        if get(rUntilExit, 'Value') == 1
            
            i = 1;
            
            while ~(x(i) > W/2 - wTop/2 && x(i) < W/2 + wTop/2 && y(i) >= L && z(i) > H/2 - lTop/2 && z(i) < H/2 + lTop/2)
                
                i = i + 1;
                
                dx = k * randn(1,1);
                dy = k * randn(1,1);
                dz = k * randn(1,1);
                
                x = [x, x(i - 1) + dx];
                y = [y, y(i - 1) + dy];
                z = [z, z(i - 1) + dz];
                
                if ~(x(i) > W/2 - wTop/2 && x(i) < W/2 + wTop/2 && y(i) >= L && z(i) > H/2 - lTop/2 && z(i) < H/2 + lTop/2)
                    if x(i) < 0
                        x(i) = x(i) + 2 * abs(x(i));
                    elseif x(i) > W
                        x(i) = x(i) - 2 * (abs(x(i)) - W);
                    end
                    
                    if y(i) < 0
                        y(i) = y(i) + 2 * abs(y(i));
                    elseif y(i) > L
                        y(i) = y(i) - 2 * (abs(y(i) - L));
                    end
                    
                    if z(i) < 0
                        z(i) = z(i) + 2 * abs(z(i));
                    elseif z(i) > H
                        z(i) = z(i) - 2 * (abs(z(i)) - H);
                    end
                end
            end
            
            set(editN, 'String', num2str(i));
            
        elseif get(rUponExit, 'Value') == 1
            
            i = 1;
            
            while ~(x(i) > W/2 - wTop/2 && x(i) < W/2 + wTop/2 && y(i) >= L && z(i) > H/2 - lTop/2 && z(i) < H/2 + lTop/2) && i < N
                
                i = i + 1;
                
                dx = k * randn(1,1);
                dy = k * randn(1,1);
                dz = k * randn(1,1);
                
                x = [x, x(i - 1) + dx];
                y = [y, y(i - 1) + dy];
                z = [z, z(i - 1) + dz];
                
                if ~(x(i) > W/2 - wTop/2 && x(i) < W/2 + wTop/2 && y(i) >= L && z(i) > H/2 - lTop/2 && z(i) < H/2 + lTop/2)
                    if x(i) < 0
                        x(i) = x(i) + 2 * abs(x(i));
                    elseif x(i) > W
                        x(i) = x(i) - 2 * (abs(x(i)) - W);
                    end
                    
                    if y(i) < 0
                        y(i) = y(i) + 2 * abs(y(i));
                    elseif y(i) > L
                        y(i) = y(i) - 2 * (abs(y(i) - L));
                    end
                    
                    if z(i) < 0
                        z(i) = z(i) + 2 * abs(z(i));
                    elseif z(i) > H
                        z(i) = z(i) - 2 * (abs(z(i)) - H);
                    end
                end
            end
            
            set(editN, 'String', num2str(i));
            
        else
            i = 1;
            
            while i < N
                
                i = i + 1;
                
                dx = k * randn(1,1);
                dy = k * randn(1,1);
                dz = k * randn(1,1);
                
                x = [x, x(i - 1) + dx];
                y = [y, y(i - 1) + dy];
                z = [z, z(i - 1) + dz];
                
                if ~(x(i) > W/2 - wTop/2 && x(i) < W/2 + wTop/2 && y(i) >= L && z(i) > H/2 - lTop/2 && z(i) < H/2 + lTop/2)
                    if x(i) < 0
                        x(i) = x(i) + 2 * abs(x(i));
                    elseif x(i) > W
                        x(i) = x(i) - 2 * (abs(x(i)) - W);
                    end
                    
                    if y(i) < 0
                        y(i) = y(i) + 2 * abs(y(i));
                    elseif y(i) > L
                        y(i) = y(i) - 2 * (abs(y(i) - L));
                    end
                    
                    if z(i) < 0
                        z(i) = z(i) + 2 * abs(z(i));
                    elseif z(i) > H
                        z(i) = z(i) - 2 * (abs(z(i)) - H);
                    end
                end
            end
        end
        
        
        w1x = [W/2 - wTop/2, W/2 + wTop/2, W/2 + wTop/2, W/2 - wTop/2, W/2 - wTop/2];
        w1y = [0, 0, 0, 0, 0];
        w1z = [H/2 + lTop/2, H/2 + lTop/2, H/2 - lTop/2, H/2 - lTop/2, H/2 + lTop/2];
        
        w2x = w1x;
        w2y = [L, L, L, L, L];
        w2z = w1z;
        
        hold off;
        plot3(w1x, w1y, w1z, 'r');
        hold on;
        plot3(w2x, w2y, w2z, 'r');
        hold on;
        plot3(x, y, z);
        grid on;
        axis([0 W 0 L 0 H]);
        grid on;
        
        title('Particle Track of a Single Simulated Particle');
        xlabel('X Position');
        ylabel('Y Position');
        zlabel('Z Position');
        
        
    end


    function randomSeed(source, eventdata)
        seed = num2str(rand(1,1)*100);
        set(editS, 'String', seed);
    end

end