% Adam Gleichsner (amg188@case.edu)
% Case Western Reserve University

cTheta = 0:.1:2 * pi;
r = .0000025;
wW = .000001;
cx = r * cos(cTheta);
cy = r * sin(cTheta);
bX = [wW/2, -wW/2, -wW/2, wW/2];
bTheta = acos(wW/(2*r));
bDif = pi/2 - bTheta;
bY = [r * sin(bTheta), r * sin(bTheta + 2 * bDif), r * sin(bTheta + pi), r * sin(bTheta + 2 * bDif + pi)];

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
        
placeHolder = 1;

hold off;
plot(cx,cy);
axis([-r r -r r]);
hold on;
plot(bX, bY, 'o');

i = 1;
x = 0;
y = -r;
d = r;

while ~(x(i) > (d * cos(bTheta + bDif))/2 && x(i) < (d * cos(bTheta - bDif))/2 && y(i) >0 && d > r)

    i = i + 1;

    dx = k * randn(1,1);
    dy = k * randn(1,1);

    x = [x, x(i - 1) + dx];
    y = [y, y(i - 1) + dy];
    
    d = sqrt((x(i) - 0)^2 + (y(i) - 0)^2);

    if ~(x(i) > (d * cos(bTheta + bDif))/2 && x(i) < (d * cos(bTheta - bDif))/2 && y(i) >0) && d > r
        if y(i) > 0
            m = (y(i) - y(i-1))/(x(i) - x(i - 1));
            b = y(i) - m * x(i);
            
            yp = (m * sqrt(m^2 * r^2 - b^2 + r^2) + b)/(m^2 + 1);
            xp = sqrt(r^2 + yp^2);
            
            m2 = (-xp)/sqrt(r^2 + xp^2);
            b2 = yp - m2 * xp;
            
            uV = [x(i) - x(i - 1), y(i) - y(i - 1)];
            vV = [xp - x(i), yp - y(i)];

            vTheta = acos((uV(1) * vV(1) + uV(2) * vV(2))/(sqrt(uV(1)^2 + uV(2)^2) * sqrt(vV(1)^2 + vV(2)^2)));
            
            
            
            plot(x(placeHolder:i-1), y(placeHolder: i - 1));
            hold on;
            plot([x(i - 1), -W/2], [y(i - 1), newY]);

            x(i) = x(i) + 2 * (abs(x(i)) - W/2);

            plot([-W/2, x(i)], [newY, y(i)])
            placeHolder = i;
        elseif t(i) < 0
            m = (y(i) - y(i-1))/(x(i) - x(i - 1));
            b = y(i) - m * x(i);
            newY = m * W/2 + b;

            plot(x(placeHolder:i-1), y(placeHolder: i - 1));
            hold on;
            plot([x(i - 1), W/2], [y(i - 1), newY]);

            x(i) = x(i) - 2 * (abs(x(i)) - W/2);

            plot([W/2, x(i)], [newY, y(i)])
            placeHolder = i;
        end

        if y(i) < 0
            m = (y(i) - y(i-1))/(x(i) - x(i - 1));
            b = y(i) - m * x(i);
            newX = (0 - b) / m;

            plot(x(placeHolder:i-1), y(placeHolder: i - 1));
            hold on;
            plot([x(i - 1), newX], [y(i - 1), 0]);

            y(i) = y(i) + 2 * abs(y(i));

            plot([newX, x(i)], [0, y(i)])
            placeHolder = i;
        elseif y(i) > L
            m = (y(i) - y(i-1))/(x(i) - x(i - 1));
            b = y(i) - m * x(i);
            newX = (L - b) / m;

            plot(x(placeHolder:i-1), y(placeHolder: i - 1));
            hold on;
            plot([x(i - 1), newX], [y(i - 1), L]);

            y(i) = y(i) - 2 * (abs(y(i) - L));

            plot([newX, x(i)], [L, y(i)])
            placeHolder = i;
        end
    else
        m = (y(i) - y(i-1))/(x(i) - x(i - 1));
        b = y(i) - m * x(i);
        newX = (L - b) / m;

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


% Line of reflection is found by the derivative of the curve
% Angles between before/after trajectory and wall is equal
% Breach is determined by distance of particle from center > radius
% Exit when distance > radius and x is within window
% Determine scaling of window dependent on distance from center



