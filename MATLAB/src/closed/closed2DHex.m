% Adam Gleichsner (amg188@case.edu)
% Case Western Reserve University

function closed2DHex()

f = figure('Position', [1 1 900 800], 'Toolbar', 'figure', 'Visible', 'off', 'MenuBar', 'none', 'Name', 'Random Walk in Closed 2D');
movegui(f, 'center')
axhan2 = axes('Units', 'Pixels', 'Position', [50, 150 , 600, 600]);
staticN = uicontrol('Style', 'text', 'Position', [675 700 100 30], 'String', 'Steps(N) =');
editN = uicontrol('Style', 'edit', 'String', '1000', 'Position', [775 700 100 31]);
staticLs = uicontrol('Style', 'text', 'Position', [675 500 50 30], 'String', 'Ls = ');
editLs = uicontrol('Style', 'edit', 'String', '.00001', 'Position', [725 500 50 31]);
staticUpWin = uicontrol('Style', 'text', 'Position', [675 450 50 30], 'String', 'UpWin Width = ');
editUpWin = uicontrol('Style', 'edit', 'String', '.000005', 'Position', [725 450 50 31]);
staticBotWin = uicontrol('Style', 'text', 'Position', [775 450 50 30], 'String', 'BotWin Width = ');
editBotWin = uicontrol('Style', 'edit', 'String', '.000005', 'Position', [825 450 50 31]);
buttonUpdate = uicontrol('Style', 'pushbutton', 'String', 'Update', 'Position',[675 350 200 40], 'Callback', @updateGraph);
%buttonDisp = uicontrol('Style', 'pushbutton', 'String', 'Plot Displacement', 'Position',[675 250 200 40], 'Callback', @plotDisplacement);
h = uibuttongroup('visible','off','Position',[0 0 0.0001 0.1], 'SelectionChangeFcn', @switchButton);
rUponExit = uicontrol('Style','Radio','String','Stop plot upon exit', 'pos',[675 200 200 30], 'parent',h);
rUntilExit = uicontrol('Style','Radio','String','Plot until exit', 'Selected', 'on','pos',[675 175 200 30], 'parent',h);
staticTau = uicontrol('Style', 'text', 'Position', [675 100 50 30], 'String', 'dT = ');
editTau = uicontrol('Style', 'edit', 'String', '.001', 'Position', [725 100 50 31]);

set(f, 'Visible', 'on');
set(h,'SelectedObject',[]);
set(h,'Visible','on');

title('Particle Track of a Single Simulated Particle');
xlabel('X Position');
ylabel('Y Position');

    function plotDisplacement(source, eventdata)
        
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
        
        hold off;
        
        N = str2num(get(editN, 'String'));
        Ls = str2num(get(editLs, 'String'));
        wTop = str2num(get(editUpWin, 'String'));
        wBot = str2num(get(editBotWin, 'String'));
        
        placeHolder = 1;
        
        dimensions = 2;         % two dimensional simulation
        tau = .001;%str2num(get(editTau, 'String'));               % time interval in seconds
        
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
                                
                                plot(x(placeHolder:i - 1), y(placeHolder:i - 1));
                                hold on;
                                plot([x(i-1), xp], [y(i - 1), yp]);
                                hold on;
                                plot([xp, x(i)], [yp, y(i)]);
                                hold on;
                                placeHolder = i;
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
                                
                                plot(x(placeHolder:i - 1), y(placeHolder:i - 1));
                                hold on;
                                plot([x(i-1), xp], [y(i - 1), yp]);
                                hold on;
                                plot([xp, x(i)], [yp, y(i)]);
                                hold on;
                                placeHolder = i;
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
                                
                                plot(x(placeHolder:i - 1), y(placeHolder:i - 1));
                                hold on;
                                plot([x(i-1), xp], [y(i - 1), yp]);
                                hold on;
                                plot([xp, x(i)], [yp, y(i)]);
                                hold on;
                                placeHolder = i;
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
                                
                                plot(x(placeHolder:i - 1), y(placeHolder:i - 1));
                                hold on;
                                plot([x(i-1), xp], [y(i - 1), yp]);
                                hold on;
                                plot([xp, x(i)], [yp, y(i)]);
                                hold on;
                                placeHolder = i;
                            end

                        end
                    end
                end
            end
            
            set(editN, 'String', num2str(i));
%         
%         elseif get(rUponExit, 'Value') == 1
%             
%             i = 1;
%             x = 0;
%             y = 0;
%             
%             while ~(x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= 3 ^ .5 * Ls) && i < N
%                 
%                 i = i + 1;
%                 
%                 dx = k * randn(1,1);
%                 dy = k * randn(1,1);
%                 
%                 x = [x, x(i - 1) + dx];
%                 y = [y, y(i - 1) + dy];
%                 
%                 if ~(x(i) > -wTop/2 && x(i) < wTop/2 && y(i) >= L)
%                     if x(i) < -W/2
%                         slope = (y(i) - y(i-1))/(x(i) - x(i - 1));
%                         b = y(i) - slope * x(i);
%                         newY = slope * -W/2 + b;
%                         
%                         plot(x(placeHolder:i-1), y(placeHolder: i - 1));
%                         hold on;
%                         plot([x(i - 1), -W/2], [y(i - 1), newY]);
%                         
%                         x(i) = x(i) + 2 * (abs(x(i)) - W/2);
%                         
%                         plot([-W/2, x(i)], [newY, y(i)])
%                         placeHolder = i;
%                     elseif x(i) > W/2
%                         slope = (y(i) - y(i-1))/(x(i) - x(i - 1));
%                         b = y(i) - slope * x(i);
%                         newY = slope * W/2 + b;
%                         
%                         plot(x(placeHolder:i-1), y(placeHolder: i - 1));
%                         hold on;
%                         plot([x(i - 1), W/2], [y(i - 1), newY]);
%                         
%                         x(i) = x(i) - 2 * (abs(x(i)) - W/2);
%                         
%                         plot([W/2, x(i)], [newY, y(i)])
%                         placeHolder = i;
%                     end
%                     
%                     if y(i) < 0
%                         slope = (y(i) - y(i-1))/(x(i) - x(i - 1));
%                         b = y(i) - slope * x(i);
%                         newX = (0 - b) / slope;
%                         
%                         plot(x(placeHolder:i-1), y(placeHolder: i - 1));
%                         hold on;
%                         plot([x(i - 1), newX], [y(i - 1), 0]);
%                         
%                         y(i) = y(i) + 2 * abs(y(i));
%                         
%                         plot([newX, x(i)], [0, y(i)])
%                         placeHolder = i;
%                     elseif y(i) > L
%                         slope = (y(i) - y(i-1))/(x(i) - x(i - 1));
%                         b = y(i) - slope * x(i);
%                         newX = (L - b) / slope;
%                         
%                         plot(x(placeHolder:i-1), y(placeHolder: i - 1));
%                         hold on;
%                         plot([x(i - 1), newX], [y(i - 1), L]);
%                         
%                         y(i) = y(i) - 2 * (abs(y(i) - L));
%                         
%                         plot([newX, x(i)], [L, y(i)])
%                         placeHolder = i;
%                     end
%                 else
%                     slope = (y(i) - y(i-1))/(x(i) - x(i - 1));
%                     b = y(i) - slope * x(i);
%                     newX = (L - b) / slope;
%                     
%                     if~(newX > -wTop/2 && newX < wTop/2)
%                         plot(x(placeHolder:i-1), y(placeHolder: i - 1));
%                         hold on;
%                         plot([x(i - 1), newX], [y(i - 1), L]);
%                         
%                         y(i) = y(i) - 2 * (abs(y(i) - L));
%                         
%                         plot([newX, x(i)], [L, y(i)])
%                         placeHolder = i;
%                     end
%                 end
%             end
%             
%             set(editN, 'String', num2str(i));
%             
%         end
        
        plot(-wTop/2, 3 ^ .5 * Ls, 'ro', wTop/2, 3 ^ .5 * Ls, 'ro', -wBot/2, 0, 'ro', wBot/2, 0, 'ro');
        hold on;
        plot([-Ls/2, -Ls, -Ls/2, Ls/2, Ls, Ls/2, -Ls/2], [0, 3 ^ .5 * Ls / 2, 3 ^ .5 * Ls, 3 ^ .5 * Ls, 3 ^ .5 * Ls / 2, 0, 0], 'r');
        hold on;
        plot(x(placeHolder:length(x)), y(placeHolder:length(y)));
        hold on;
        axis([-Ls Ls -.15*Ls 1.85*Ls]);
        
        title('Particle Track of a Single Simulated Particle');
        xlabel('X Position');
        ylabel('Y Position');
        
        
    end

    function randomSeed(source, eventdata)
        
        seed = num2str(rand(1,1)*100);
        set(editS, 'String', seed);
        
    end

end