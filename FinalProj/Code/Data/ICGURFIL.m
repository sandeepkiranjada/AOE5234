        %
        % Start date and time of simulation
        %
        Mjd_UTC_Epoch = Mjday(2015, 01, 01, 00, 00, 00);
        %
        % Spacecraft properties
        %
        Cd = 2.2;                                   % Drag coefficient of the spacecraft
        AMR = 0.02;                                 % S/m: Area to mass ratio of the spacecraft [m^2/kg]
        delta = 0.5*AMR*Cd;                         % Ballistic coefficient;
        %
        % Orbit properties
        %
        hp0 = 250*1e3;                              % Initial perigee height [m]
        ha0 = 35943*1e3;                            % Initial apogee hegiht [m]
        r_p0 = Re+hp0;                              % Initial perigee radius [m]
        r_a0 = Re+ha0;                              % Initial apogee radius [m]
        %
        a0 = (r_a0+r_p0)/2;                         % Initial semi-major axis [m]
        e0 = 1-r_p0/a0;                             % Initial eccentricity
        i0 = deg2rad(6);                            % Initial inclination [rad]
        argp0 = deg2rad(178);                       % Initial argument of perigee [rad]
        raan0 = deg2rad(60);                        % Initial RAAN [rad]
        M0 = 0;                                     % Initial mean anomaly [rad]
        %
        H0 = sqrt(a0*mu_earth*(1-e0^2));            % Initial angular momentum
        %
        % Atmospheric properties
        %
        we = 7.2921159e-5;                          % Angular velocity of the earth [rad/s]
        wa = we;                                    % Angular velocity of the atmosphere in z-direction [rad/s]
        %