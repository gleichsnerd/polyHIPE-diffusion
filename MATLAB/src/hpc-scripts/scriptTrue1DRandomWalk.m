% Case Western Reserve University
% Adam Gleichsner (amg188)

function scriptTrue1DRandomWalk(columns, particles)

tic
    TotalEscapeTime = [];
    dimensions = 2;         % two dimensional simulation
    tau = .001;               % time interval in seconds

    R    = .145e-6;              % radius in meters 
    %eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K 
    kB   = 1.38e-23;            % Boltzmann constant T    = 293;                 % Temperature in degrees Kelvin
    eta  = 1.827e-5;              % viscosity of air in SI units (Pascal-seconds) at 293 K 
    kB   = 1.38e-23;            % Boltzmann constant
    T    = 293;                 % Temperature in degrees Kelvin

    D    = kB * T / (6 * pi * eta * R);
    k = sqrt(D * dimensions * tau);

    i = 1;
    
    x(i) = 1;
    
    for a = 1:particles
        while(x(i) <= rows)
            direction = round(0 + (2-(0))*rand(1,1));
            if direction == 1 || x(i) == 1
                x(i) = x(i) + 1;
            else
                x(i) = x(i) - 1;
            end
        end
        TotalEscapeTime = [TotalEscapeTime, length(x)];
    end
    
    filename = strcat('/home/amg188/MATLAB/results/Group', ' ', num2str(num), '/1DWalk', '/walk', num2str(batch(a)), '_1.txt');
    save(filename, 'histArray', 'bins', 'D', 'tau', '-double', '-ascii');
toc

end