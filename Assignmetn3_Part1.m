%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        ELEC 4700 A - Winter 2022                        %
%                              Assignment 3                               %
%        Author: Julie-Anne Chaine             Student #: 101104568       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% Constants:
mo = 9.1093837015e-31;              % Electron rest mass in kg
mn = 0.26*mo;                       % Effective mass of electrons
k = 1.38e-23;                       % Boltzmann's constant 
W = 200e-9;                         % Nominal region width in m
L = 100e-9;                         % Nominal region length in m
tau_mn = 0.2e-12;                   % Mean time between collisions
q = 1.602e-19;                      % Electron charge

%Variables
T = 300;                            % Room temperature in K
elecpop = 1000;                    % Number of particles to simulate
dt = 1e-14;                         % Time step in s
x = zeros(elecpop,1);               % Electron y positions
y = zeros(elecpop,1);               % Electron x positions
oldx = zeros(elecpop,1);            % Previous electron x position
oldy = zeros(elecpop,1);            % Previous electron y position
Vth = sqrt(2*k*T/mn);               % Thermal velocity in m/s
mfp = Vth*tau_mn;                   % Mean free path in m
samplepop = 10;                     % # of particles to plot
samp = randi(elecpop,samplepop,1);  % Random particles to observe
iter = 100;                         % # of iterations 
Pscatter = 1 - exp(-dt/tau_mn);     % Scattering probability

E = (0.1 - 0)/W;    % (Vleft - Vright)/d
F = E*q;
a = F/mn;

concentration = 10e15;

% Initializing positions and velocities 
x(:,1) = rand(elecpop,1)*W;
y(:,1) = rand(elecpop,1)*L;
vx = Vth*randn(elecpop,1)/sqrt(2);
vy = Vth*randn(elecpop,1)/sqrt(2); 
    
    for j = 1:iter 
        oldx = x;
        oldy = y;
        PartScatter = Pscatter > rand(elecpop,1);
        vx(PartScatter) = Vth*randn(sum(PartScatter),1)/sqrt(2);
        vy(PartScatter) = Vth*randn(sum(PartScatter),1)/sqrt(2);
        vx = vx + a*dt;
        x = oldx + vx*dt + 1/2*a*dt^2;       % Previous x position + delta L
        y = oldy + vy*dt;       % Previous y position + delta L    
        oldx(x<0) = W;          % Making the particles on the left boundary, appear on the right
        oldx(x>W) = 0;          % Making the particles on the right boundary, appear on the left
        x(x<0) = x(x<0) + W;    % All points passed the left get pushed to the other side
        x(x>W) = x(x>W) - W;    % All points passed the right get pushed to the other side
        vy(y>L) = -vy(y>L);     % Particle Y direction gets flipped if it hits the top
        vy(y<0) = -vy(y<0);     % Particle Y direction gets flipped
        
        for m = 1:samplepop         % For plotting
            figure(2)
            plot([oldx(samp(m)), x(samp(m))],[oldy(samp(m)),y(samp(m))],'SeriesIndex',m)
            title('Sample Particles Trajectories')
            xlabel('W (m) or x position')
            ylabel('L (m) or y position')
            axis([0 200e-9 0 100e-9])
            hold on
        end
        Vmean = mean( sqrt( vx.^2 + vy.^2 ) );
        V = sqrt( vx.^2 + vy.^2 );
        J(j) = q*concentration*Vmean*L;
    end 
        plot(J)
        title('Current in x direction')
        xlabel('Time')
        ylabel('Current')

% Electron Density map:
nx = 20; ny = 30;
BinX = ceil( x*nx/W );
BinY = ceil( y*ny/L );
for i = 1:nx
    for j = 1:ny
        match = BinX == i & BinY == j;
        Ebox(i,j) = sum((match));
    end
end
figure
surf(Ebox)
xlabel('Width (m)')
ylabel('Length (m)')
zlabel('Density')
title('Electron density map')       

% Temperature map
for i = 1:nx
    for j = 1:ny
        match = BinX == i & BinY == j;
        Tbox(i,j) = sum( V(match).^2/(2*k) );
    end
end
figure
surf(Tbox)
title('Temperature Density')
xlabel('Width (m)')
ylabel('Length (m)')
zlabel('Temperature (K)')