%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        ELEC 4700 A - Winter 2022                        %
%                              Assignment 3                               %
%        Author: Julie-Anne Chaine             Student #: 101104568       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
ny = 30; nx = 20;
W = 20e-9; L = 30e-9;
G = sparse(nx*ny, nx*ny);
V = zeros(nx, ny);
F = zeros(nx*ny, 1);
vMap = zeros(nx,ny);
Vo = 1;

cMap = ones(nx, ny);
for i = 8:12
    for j = 1:7
        cMap(i,j) = 0.5;
    end
    for j = 22:30
        cMap(i,j) = 0.5;
    end
end

for i = 1:nx 
    for j = 1:ny
        nxm = j + (i-2)*ny;
        nxp = j + i*ny;
        nyp = j + 1 + (i-1)*ny;
        nym = j - 1 + (i-1)*ny;
        n = j + (i-1)*ny;
        if i == 1           % Left
            G(n, n) = cMap(i,j);
            F(n) = 0.1;
        elseif i == nx      % Right
            G(n, n) = cMap(i,j);
            F(n) = 0;
        elseif j == 1       % Bottom     
            G(n, nxm) = cMap(i-1,j);         
            G(n, nxp) = cMap(i+1,j);          
            G(n, nyp) = cMap(i,j+1); 
            G(n, n) = -(G(n, nxp)+G(n, nxm)+G(n, nyp));
        elseif j == ny      % Top
            G(n, nxm) = cMap(i-1,j);         
            G(n, nxp) = cMap(i+1,j);                   
            G(n, nym) = cMap(i,j-1); 
            G(n, n) = -(G(n, nxp)+G(n, nxm)+G(n, nym));
        else
            G(n, nxm) = cMap(i-1,j);         
            G(n, nxp) = cMap(i+1,j);          
            G(n, nyp) = cMap(i,j+1);         
            G(n, nym) = cMap(i,j-1);  
            G(n, n) = -(G(n, nxp)+G(n, nxm)+G(n, nym)+G(n, nyp));
        end
    end
end
V = G\F;

for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        vMap(i,j) = V(n);
    end
end
[Ey,Ex] = gradient(vMap);
Ex = -Ex*nx/W;
Ey = -Ey*ny/L;

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
elecpop = 1000;                     % Number of particles to simulate
dt = 1e-14;                         % Time step in s
x = zeros(elecpop,1);               % Electron y positions
y = zeros(elecpop,1);               % Electron x positions
oldx = zeros(elecpop,1);            % Previous electron x position
oldy = zeros(elecpop,1);            % Previous electron y position
Vth = sqrt(2*k*T/mn);               % Thermal velocity in m/s 
mfp = Vth*tau_mn;                   % Mean free path in m
samplepop = 10;                     % # of particles to plot
samp = randi(elecpop,samplepop,1);  % Random particles to observe
iter = 1000;                         % # of iterations 
Pscatter = 1 - exp(-dt/tau_mn);     % Scattering probability
ax = zeros(elecpop,1);
ay = zeros(elecpop,1);
concentration = 10e15;
%Comment out to test bottleneck size effect on current!!!!
% top_bottleneck = 0.6e-7;
% bottom_bottleneck = 0.4e-7;
% top_bottleneck = 0.55e-7;
% bottom_bottleneck = 0.45e-7;
top_bottleneck = 0.2e-7;
bottom_bottleneck = 0.8e-7;


% Initializing positions and velocities 
x(:,1) = rand(elecpop,1)*W;
y(:,1) = rand(elecpop,1)*L;
vx = Vth*randn(elecpop,1)/sqrt(2);
vy = Vth*randn(elecpop,1)/sqrt(2); 
    for z = 1:elecpop   % Generating points on the left side
       x(z,1) = 0.2e-7 + (0.8e-7 - 0.2e-7)*rand();
    end
    for z = elecpop/2+1:elecpop     % Generating points on the right
       x(z,1) = 1.2e-7 + (2e-7 - 1.2e-7)*rand();        
    end


    
%Boxes
Boxes = {};
Boxes{1}.X = [0.8e-7, 0.4e-7];
Boxes{1}.Y = [0e-7, bottom_bottleneck];
Boxes{2}.X = [0.8e-7, 0.4e-7];
Boxes{2}.Y = [top_bottleneck, 1e-7];    

    for j = 1:iter 
        
        Ybin = discretize(y,ny);
        Xbin = discretize(x,nx);
        for b = 1:elecpop
            Ex_particle = Ex(Xbin(b), Ybin(b));
            Ey_particle = Ey(Xbin(b), Ybin(b));
            Fx = Ex_particle*q;
            ax(b) = Fx/mn;
            Fy = Ey_particle*q;
            ay(b) = Fy/mn;
        end
        %Creating box reflection
        vx(oldx>=1.2e-7 & x<=1.2e-7 & (y>top_bottleneck | y<bottom_bottleneck)) = -vx(oldx>=1.2e-7 & x<=1.2e-7 & (y>top_bottleneck | y<bottom_bottleneck));
        vx(oldx<=0.8e-7 & x>=0.8e-7 & (y>top_bottleneck | y<bottom_bottleneck)) = -vx(oldx<=0.8e-7 & x>=0.8e-7 & (y>top_bottleneck | y<bottom_bottleneck));        
        vy(x>0.8e-7 & x<1.2e-7 & oldy>=bottom_bottleneck & y<=bottom_bottleneck) = -vy(x>0.8e-7 & x<1.2e-7 & oldy>=bottom_bottleneck & y<=bottom_bottleneck);
        vy(x>0.8e-7 & x<1.2e-7 & oldy<=top_bottleneck & y>=top_bottleneck) = -vy(x>0.8e-7 & x<1.2e-7 & oldy<=top_bottleneck & y>=top_bottleneck);
        
        oldx = x;
        oldy = y;
        PartScatter = Pscatter > rand(elecpop,1);
        vx(PartScatter) = Vth*randn(sum(PartScatter),1)/sqrt(2);
        vy(PartScatter) = Vth*randn(sum(PartScatter),1)/sqrt(2); 
        x = oldx + vx*dt + 1/2*ax*dt^2;       % Previous x position + delta L + acc
        y = oldy + vy*dt + 1/2*ay*dt^2;       % Previous y position + delta L + acc       
        vx = vx + ax*dt;
        vy = vy + ay*dt;
        oldx(x<0) = W;          % Making the particles on the left boundary, appear on the right
        oldx(x>W) = 0;          % Making the particles on the right boundary, appear on the left
        x(x<0) = x(x<0) + W;    % All points passed the left get pushed to the other side
        x(x>W) = x(x>W) - W;    % All points passed the right get pushed to the other side
        vy(y>L) = -vy(y>L);     % Particle Y direction gets flipped if it hits the top
        vy(y<0) = -vy(y<0);     % Particle Y direction gets flipped
        % If the particle enters the box:
        x(oldx>=1.2e-7 & x<=1.2e-7 & (y>0.6e-7 | y<0.4e-7)) = 1.2e-7;
        x(oldx<=0.8e-7 & x>=0.8e-7 & (y>0.6e-7 | y<0.4e-7)) = 0.8e-7;
        y(oldy>=0.4e-7 & y<=0.4e-7 & x<=1.2e-7 & x>0.8e-7) = 0.4e-7;
        y(oldy<=0.6e-7 & y>=0.6e-7 & x<=1.2e-7 & x>0.8e-7) = 0.6e-7;
        
        V = mean( sqrt( vx.^2 + vy.^2 ) );
        J(j) = q*concentration*V*L;
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
