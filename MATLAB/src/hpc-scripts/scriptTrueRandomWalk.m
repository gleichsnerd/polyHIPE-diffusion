% Case Western Reserve University
% Adam Gleichsner (amg188)

function scriptTrue1DRandomWalk(columns)

tic

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
    
    x(i) = sx;
    y(i) = sy;
    
    while(x(i) <= rows)
        bwall = true;
        twall = true;
        lwall = true;
        rwall = true;
        wallFactor = 0;
        if x(i) = 0
            lwall = false;
            wallFactor = wallFactor + 1;
        end
        if x(i) = columns
            rwall = false;
            wallFactor = wallFactor + 1;
        end
        if y(i) = 0
            bwall = false;
            wallFactor = wallFactor + 1;
        end
        if y(i) = rows
            twall = false;
            wallFactor = wallFactor + 1;
        end
        
        numberWindows = round(1 + (2-(1))*rand(1,1)) * 2;
    end
    
    
toc

end