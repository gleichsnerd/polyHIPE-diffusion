function einstein()

    f2 = figure('Position', [1 1 900 800], 'Toolbar', 'figure', 'Visible', 'off', 'MenuBar', 'none', 'Name', 'Einsteinian Random Walk in Closed 2D');
    movegui(f2, 'center')
    axhan2 = axes('Units', 'Pixels', 'Position', [50, 150 , 600, 600]);
    staticN = uicontrol('Style', 'text', 'Position', [675 700 100 30], 'String', 'Steps(N) =');
    editN = uicontrol('Style', 'edit', 'String', '1000', 'Position', [775 700 100 31]);
    staticT = uicontrol('Style', 'text', 'Position', [675 650 50 30], 'String', 'dT = ');
    editT = uicontrol('Style', 'edit', 'String', '100', 'Position', [725 650 50 31]);
    staticT2 = uicontrol('Style', 'text', 'Position', [775 650 50 30], 'String', '* tau');
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
    rUntilExit = uicontrol('Style','Radio','String','Plot until exit', 'pos',[675 175 200 30], 'parent',h);


    set(f2, 'Visible', 'on');
    set(h,'SelectedObject',[]);
    set(h,'Visible','on');

    title('Particle Track of a Single Simulated Particle');
    xlabel('X Position');
    ylabel('Y Position');
    
    function updateGraph(source, eventdata)
        
        R    = .145e-6;              % radius in meters
        eta  = 1.0e-3;              % viscosity of air in SI units (Pascal-seconds) at 293 K
        kB   = 1.38e-23;            % Boltzmann constant
        dT = .001
        T    = 293;                 % Temperature in degrees Kelvin
        V = 4/3 * pi * R^3; %m^3
        p = 1.05 * 1000; % kg/cm
        
        m = p * V; %kg
        gamma = 6 * pi * eta * R;
        tau = m/gamma;
        D = kB * T / gamma;

        N = str2num(get(editN, 'String'));
        L = str2num(get(editL, 'String'));
        W = str2num(get(editW, 'String'));
        wTop = str2num(get(editUpWin, 'String'));
        wBot = str2num(get(editBotWin, 'String'));
        placeHolder = 1;

        x = 0; %m
        y = 0; %m
        

        N = 1000;
        hold off;
        
        if get(rUntilExit, 'Value') == 1
            
            i = 1;
            x = 0;
            y = 0;
            
            while ~(x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= L)
                
                i = i + 1;

                x = [x, x(i-1) + sqrt(2 * D * dT) * randn(1)];
                y = [y, y(i-1) + sqrt(2 * D * dT) * randn(1)];
                
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
                
                x = [x, x(i-1) + sqrt(2 * D * dT) * randn(1)];
                y = [y, y(i-1) + sqrt(2 * D * dT) * randn(1)];
                
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

    function switchButton(source, eventdata)
        
        event = get(eventdata.NewValue,'String');
        if strcmp(event, 'Stop plot upon exit')
            set(rUntilExit, 'Selected', 'off');
        elseif strcmp(event, 'Plot until exit')
            set(rUponExit, 'Selected', 'off');
        end
    end

end