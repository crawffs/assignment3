%Assignment 3-part 3
%Mike Crawford 100952432
clearvars
clearvars -GLOBAL
close all
format shortE
global C
global Ecount
global Vx Vy Vtotal x y
Ecount =1000;   %setting amount of electrons used
C.mo = 9.10938215e-31; %Electron mass
C.k = 1.3806504e-23; %Boltzmann constant
q = -1.60217662e-19;
T =300; %Temperature (K)
mn = 0.26*C.mo;
L = 200e-9;  %Length of Border
W = 100e-9; %Width of Border
Vth = sqrt((2*C.k*T)/mn); %Calculation of Thermal Velocity
dt = 10e-15;  %timestep
Stop = 100*dt; %timeframe
x = zeros(Ecount, 2);  %array of x positions
y = zeros(Ecount, 2);  %array of y positions
Temperature = zeros(1,2); %array of temperatures
Time = 0;
VisibleEcount = 50; %visibble electrons
tmn = 0.2e-12; %mean time between collisions
PScat = 1 - exp(-dt/tmn);  %scattering probability equation
VTotalHistogram = zeros(Ecount, 1); %array of thermal velocities
DrawBoxX = [80e-9 80e-9 120e-9 120e-9 80e-9]; %X border limits of boxes
%(Same X values for both)
DrawBoxY1 = [100e-9 60e-9 60e-9 100e-9 100e-9];%Y border limits of top box
DrawBoxY2 = [40e-9 0 0 40e-9 40e-9]; %Y border limits of bottom box
Specular = true; %boolean indicating whether specular or diffusive
InsideBox = true; %boolean showing if electron spawned inside boxes
MappingStep = 10e-9; %size of each suare used on mapping grid
%setting arrays of density and temperature based off border demesions and
%mapping step
DensityMap = zeros(W/MappingStep, L/MappingStep);
TemperatureMap = zeros(W/MappingStep, L/MappingStep);

nx = 30;
ny = 20;
dx = L/nx;
dy = W/ny;
conductOut = 1;
conductBox = 01e-2;
conduct = zeros(nx,ny);
G = sparse (nx*ny, nx*ny);
V = zeros(1, nx*ny);
Voltagex = 0.1;

for i = 1:nx
    for j = 1:ny
        %sets boundaries of boxes
        if (i > (0.3*nx) || i < (0.6*nx)) && (j > (0.6*ny) || j < (0.3*ny))
            conduct(i,j) = conductBox;
        else
            conduct(i,j) = conductOut;
        end
        
    end
end
%same procedure as the previous parts,
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1)*ny;
        nx1 = j + ((i-1) - 1)*ny;
        nx2 = j + ((i+1) - 1)*ny;
        ny1 = (j-1) + (i - 1)*ny;
        ny2 = (j+1) + (i - 1)*ny;
       
        if i == 1
            V(n) = Voltagex;
            G(n,n) = 1;
        elseif i == nx
            V(n) = 0;
            G(n,n) =1;
        elseif j == 1 
            Cx1 = (conduct(i,j) + conduct(i-1,j))/2;
            Cx2 = (conduct(i,j) + conduct(i+1,j))/2;
        
            Cy2 = (conduct(i,j) + conduct(i,j+1))/2;
            G(n,n) = -Cx1 - Cx2 - Cy2;
            G(n, nx1) = Cx1;
            G(n,nx2) = Cx2;
            G(n, ny2) = Cy2;
        elseif j == ny
             Cx1 = (conduct(i,j) + conduct(i-1,j))/2;
             Cx2 = (conduct(i,j) + conduct(i+1,j))/2;
             Cy1 = (conduct(i,j) + conduct(i,j-1))/2;
            
            G(n,n) = -Cx1 - Cx2 - Cy1;
            G(n,nx1) = Cx1;
            G(n,nx2) = Cx2;
            G(n,ny1) = Cy1;
        else
             Cx1 = (conduct(i,j) + conduct(i-1,j))/2;
            Cx2 = (conduct(i,j) + conduct(i+1,j))/2;
            Cy1 = (conduct(i,j) + conduct(i,j-1))/2;
            Cy2 = (conduct(i,j) + conduct(i,j+1))/2;
            G(n,n) = -Cx1 - Cx2 - Cy1 - Cy2;
            G(n,nx1) = Cx1;
            G(n,nx2) = Cx2;
            G(n,ny1) = Cy1;
            G(n,ny2) = Cy2;
        end
    end
end

S = G\V';
surface = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1)*ny;
        surface(i,j) = S(n);
    end
end
[Ex, Ey] = gradient(-surface);

Fx = q*Ex;
Fy = q*Ey;

Accx = Fx /mn;
Accy = Fy /mn;

for i = 1:Ecount
    x(i,1) = rand()*200e-9; %setting random electron positions
    y(i,1) = rand()*100e-9;
    InsideBox = true;
    %checks if set electron is inside box, if so reinitialize its position
    while InsideBox == true
        if (x(i) >= 40e-9 && x(i) <= 120e-9) && (y(i) >= 60e-9 ||...
                y(i) <= 40e-9)
            x(i,1) = rand * 200e-9;
            y(i,1) = rand * 100e-9;
        else
            InsideBox = false;
        end
    end
    
end

for i = 1:Ecount
    
Vx(1:Ecount) = Vth * randn; %setting velocities
Vy(1:Ecount) = Vth * randn;
end

%plot of electrons with boxes
figure(1)
subplot(2,1,1);
plot(DrawBoxX, DrawBoxY1, DrawBoxX, DrawBoxY2)
axis([0 L 0 W]);
title('Part 3');
xlabel('X');
ylabel('Y');
hold on;

while Time < Stop
    subplot(2,1,1)
    %looping electrons
    for j = 1:Ecount
        %set to assume leaking to start condition for diffusive bounces
        leaking = true;
        %scattering probability condition
        if PScat> rand
                Vx(j) = Vth * randn;
                Vy(j) = Vth * randn;
        end
        
        %updating old/new x and y positions
        Axindex = round((x(j,2)/L) * 30);
        Ayindex = round((y(j,2)/W)*20);
        if Axindex < 1
            Axindex = 1;
        elseif Axindex > 30 
                Axindex = 30;
        end
        if Ayindex < 1
            Ayindex = 1;
        elseif Ayindex > 20
            Ayindex = 20;
        end
        
        Vx(j) =  Vx(j) + Accx(Axindex,Ayindex)*dt;
        Vy(j) =  Vy(j) + Accy(Axindex,Ayindex)*dt;
        x(j,2) = x(j,1);
        y(j,2) = y(j,1);
        x(j,1) = x(j,1) + (dt * Vx(j));
        y(j,1) = y(j,1) + (dt * Vy(j));
        %checking if collision with top box
        if (x(j,1) >= 80e-9 && x(j,1) <= 120e-9) && y(j,1) >= 60e-9
            %checking which side of box collision came from
            %if left side, reflect X velocity and reset positions outside 
            %of box
                if y(j,2) < 60e-9
                    Vy(j) = -Vy(j);
                    y(j,1) = 60e-9;
                    y(j,2) = 60e-9;
                    %if right side, reflect X velocity and reset position
                    %outside of box
                elseif x(j,2) < 80e-9
                    Vx(j) = -Vx(j);
                    x(j,1) = 80e-9;
                    x(j,2) = 80e-9;
                    %if bottom side, reflect Y velocity and reset position
                    %outside of box
                elseif x(j,2) > 120e-9
                    Vx(j) = -Vx(j);
                    x(j,1) = 120e-9;
                    x(j,2) = 120e-9;
                end
                %check if specular or diffusive collisions
            if Specular == true
                %if specular, simply update positions
                x(j,1) = x(j,2) + Vx(j)*dt;
                y(j,1) = y(j,2) + Vy(j)*dt;
            else
                %if diffusive, update velocities with random value
             Vx(j) = Vth * randn;
             Vy(j) = Vth * randn;
             %assume random velocity will direct it inside the box
             while leaking == true
                 %check if new velocity is directed towards the box by 
                 %examining position point and new velocity
                 if(x(j,2) < 80e-9 && Vx(j) >= 0) || ...
                         (x(j,2) > 120e-9 && Vx(j) <= 0) || ...
                         (y(j,2) < 60e-9 && Vy(j) >= 0)
                     %if new velocity conflicts with boxes, reinitialize 
                     %the velocity
                     Vx(j) = Vth * randn;
                     Vy(j) = Vth * randn;
                 else
                     leaking = false;
                 end
             end
             %update positions with random velocity
             x(j,1) = x(j,2) + Vx(j)*dt;
             y(j,1) = y(j,2) + Vy(j)*dt;
            end
        end
         %checking if collision with bottom box
        if (x(j,1) >= 80e-9 && x(j,1) <= 120e-9) && y(j,1) <= 40e-9
             %checking which side of box collision came from
            %if top side, reflect Y velocity and reset positions outside 
            %of box
                if y(j,2) > 40e-9
                    Vy(j) = -Vy(j);
                    y(j,1) = 40e-9;
                    y(j,2) = 40e-9;
                   %if left side, reflect X velocity and reset positions outside 
            %of box  
                elseif x(j,2) < 80e-9
                    Vx(j) = -Vx(j);
                    x(j,1) = 80e-9;
                    x(j,2) = 80e-9;
                     %if right side, reflect X velocity and reset position
                    %outside of box
                elseif x(j,2) > 120e-9
                    Vx(j) = -Vx(j);
                    x(j,1) = 120e-9;
                    x(j,2) = 120e-9;
                end
                 %check if specular or diffusive collisions
            if Specular == true
                %if specular, simply update positions
                x(j,1) = x(j,2) + Vx(j)*dt;
                y(j,1) = y(j,2) + Vy(j)*dt;
            else
                 %if diffusive, update velocities with random value
             Vx(j) = Vth * randn;
             Vy(j) = Vth * randn;
              %assume random velocity will direct it inside the box
             while leaking == true
                 %check if new velocity is directed towards the box by 
                 %examining position point and new velocity
                 if(x(j,2) < 80e-9 && Vx(j) >= 0) || ...
                         (x(j,2) > 120e-9 && Vx(j) <= 0) || ...
                         (y(j,2) > 40e-9 && Vy(j) <= 0)
                      %if new velocity conflicts with boxes, reinitialize 
                     %the velocity
                     Vx(j) = Vth * randn;
                     Vy(j) = Vth * randn;
                 else
                     leaking = false;
                 end
             end
              %update positions with random velocity
             x(j,1) = x(j,2) + Vx(j)*dt;
             y(j,1) = y(j,2) + Vy(j)*dt;
            end
        end
          %check right wall border
        if x(j,1) > L
            x(j,2) = 0;
            x(j,1) = dt * Vx(j);
        end
        %check left wall border
        if x(j,1) < 0
            x(j,2) = L;
            x(j,1) = x(j,2) + (dt * Vx(j));
        end
         %check roof/floor border
        if y(j,1) > W || y(j,1) < 0
            Vy(j) = -Vy(j);
        end
         %set line vectores for x and y positions
        XPlot = [x(j,2) x(j,1)];
        YPlot = [y(j,2) y(j,1)];
         %plot visible x and y line vectors
        if j < VisibleEcount
        plot(XPlot,YPlot);
        end
        
        %sum of themal velocities and temperatures
       VTotal = sqrt(Vx(j)^2 + Vy(j)^2);
       %Temperature(1,2) = Temperature(1,2) + (mn*Vtotal(j)^2)/(2*C.k);
         
    end
     %plotting average temperature and ressetting values for next timestep
    AvgTemperature = Temperature(1,2)/Ecount;
    TemperaturePlot = [Temperature(1,1) AvgTemperature];
    TimePlot = [(Time - dt) Time];
    subplot(2,1,2);
    plot(TimePlot, TemperaturePlot);
    Temperature(1,1) = AvgTemperature;
    AvgTemperature = 0;
    Temperature(1,2) = 0;
    pause(1e-19)
    Time = Time + dt;
end 
for i = 1:(L/MappingStep)
    %looping through Y indexes of the Mapping grid
    for j = 1:(W/MappingStep)
        %looping through electrons
        for m = 1:Ecount
            %Checks to see if the electons position is with the map grip
            %being indexed. does this by checking all 4 sides of the 
            %square grid index
            if(x(m,1) > MappingStep*(i -1)) && ...
                    (x(m,1) < MappingStep*(i)) && ...
                    (y(m,1) > MappingStep*(j - 1)) && ...
                    (y(m,1) < MappingStep*(j))
                %sum of x and y velocities
                Vtotal(m) = sqrt(Vx(m)^2 + Vy(m)^2);
                %counts increased on the indexed location of both Maps
                DensityMap(j, i) = DensityMap(j, i) + 1;
                TemperatureMap(j, i) = TemperatureMap(j,i) + ...
                    (mn*Vtotal(m)^2)/(2*C.k);
            end
            %Averages the temperature values given in each grid index
            TemperatureMap(j,i) = TemperatureMap(j,i)/DensityMap(j,i);
        end
    end
end

figure(3)
imagesc(DensityMap)
title('Density map of Electrons')
xlabel('X (5nm)');
ylabel('Y (5nm)');
set(gca, 'Ydir', 'Normal')
c = colorbar;
title(c, 'Electron Count')