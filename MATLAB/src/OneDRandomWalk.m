% Adam Gleichsner (amg188@case.edu)
% Case Western Reserve University

function OneDRandomWalk()

f = figure('Position', [1 1 900 800], 'Toolbar', 'figure', 'Visible', 'off', 'MenuBar', 'none', 'Name', 'Random Walk in Closed 2D');
movegui(f, 'center')
axhan = axes('Units', 'Pixels', 'Position', [50, 100 , 600, 600]);
staticNx = uicontrol('Style', 'text', 'Position', [675 650 100 30], 'String', 'Nx =');
editNx = uicontrol('Style', 'edit', 'String', '1000', 'Position', [775 650 100 31]);
staticn = uicontrol('Style', 'text', 'Position', [675 700 100 30], 'String', 'n =');
editn = uicontrol('Style', 'edit', 'String', '1000', 'Position', [775 700 100 31]);
staticT = uicontrol('Style', 'text', 'Position', [675 600 100 30], 'String', 'Time Step =');
editT = uicontrol('Style', 'edit', 'String', '.35', 'Position', [775 600 100 31]);
buttonUpdate = uicontrol('Style', 'pushbutton', 'String', 'Update', 'Position',[675 350 200 40], 'Callback', @updateGraph);

set(f, 'Visible', 'on');

title('One Dimensional Random Walk');
xlabel('Step');
ylabel('Position');

    function updateGraph (source, eventdata)
        
        n = str2num(get(editn, 'String'));
        Nx = str2num(get(editNx, 'String'));
        p = 0;
        
        for i = 1:n
        
            direction = -1 + (1 - (-1)) * rand(1,1);
            if direction > 0
                if p(i) + 1 > Nx/2
                    p = [p, p(i)];
                else
                    p = [p, p(i) + 1];
                end
            elseif direction < 0
                if p(i) - 1 < -Nx / 2
                    p = [p, p(i)];
                else
                    p = [p, p(i) - 1];
                end
            else
                p = [p, p(i)];
            end
            progressbar(i/n);
        end
            
            hold off;
            bar(p);
            hold on;
            axis([0 n (min(p)-5) (max(p)+5)]);
        
    end

end