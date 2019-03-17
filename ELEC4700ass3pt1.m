
%Assignment 3-part 1
%Mike Crawford
%reset
clearvars
clearvars -GLOBAL
close all
format shortE
global C
global Ecount
global Vx Vy Vtotal x y
Ecount =100;   %setting amount of electrons used
C.mo = 9.10938215e-31;  %Electron mass
C.k = 1.3806504e-23;    %Boltzmann constant
q = -1.60217662e-19;
T =300;     %Temperature (K)
mn = 0.26*C.mo;
Voltagex = 0.1; %initial voltage condition

L = 200e-9; %Length of Border
W = 100e-9; %Width of Border
Edensity = 1e15*100^2;

Vth = sqrt(2*(C.k*T)/mn); %Calculation of Thermal Velocity
dt = 10e-15; %timestep
Stop =500*dt; %timeframe
x = zeros(Ecount, 2);   %array of x positions
y = zeros(Ecount, 2);   %array of y positions

Temperature = zeros(1,2);   %array of temperatures
Time = 0;
VisibleEcount = 10; %number of visible electrons
Ex = Voltagex/L;    %E field in x and y direction
Ey =  0;

Fx = q*Ex;          %force in x and y direction
Fy = q*Ey;

Accx = Fx /mn;      %Acceleration in x and y direction
Accy = Fy /mn;

Pscatter = 1 - exp(dt/0.2e-12);
Ix = zeros(1,Ecount); %initialize matrix for current values
MappingStep = 10e-9; %size of each suare used on mapping grid
%setting arrays of density and temperature based off border demesions and
%mapping step
DensityMap = zeros(W/MappingStep, L/MappingStep);
TemperatureMap = zeros(W/MappingStep, L/MappingStep);
for i = 1:Ecount %setting random locations for eact electron
    x(i,1) = rand()*200e-9; 
    y(i,1) = rand()*100e-9;
end
for i = 1:Ecount
    dir = 2*pi*rand;
Vx(1:Ecount) = Vth * cos(dir); %setting velocities
Vy(1:Ecount) = Vth * sin(dir);
end
for i = 1:Ecount
    Vtotal(i) = sqrt(Vx(i)^2 + Vy(i)^2); %summing x and y velocities
end


figure(1)   %plotting electrons
subplot(2,1,1);
axis([0 L 0 W]);
title('Trajectory of electron in Silicon');
xlabel('X');
ylabel('Y');
hold on;

subplot(2,1,2); %plotting avg temperature
axis([0 Stop 0 500]);
title('Semiconductor Temperature');
xlabel('Time (s)');
ylabel('Temperature (K)');
hold on;

for i = 1:Ecount %Sum of all initial temperatures
   Temperature(1,2) = Temperature(1,2) + (mn*Vtotal(i)^2)/(2*C.k);
end
Vdrift = mean(Vx); %Calculate aveage drift velocity
Ix(1) = L * W * q * Edensity * Vdrift;
%intitial average temperature
AvgTemperature = Temperature(1,2)/Ecount;
TemperaturePlot = [300 AvgTemperature]; %setting temperature line to plot
TimePlot = [0 Time];    %setting time line to plot
plot(Time, AvgTemperature);

%reseting values before loop
VTotal = 0;
Temperature(1,2) = 0;
AvgTemperature = 0;

%looping through timesteps

for i = 2:Ecount
    subplot(2,1,1)
    for j = 1:Ecount %looping through electrons
        if Pscatter > rand
            Vx(j) = Vt * randn;
            Vy(j) = Vt * randn;
        end
        %updating velocity with new acceleration
        Vx(j) =  Vx(j) + Accx*dt;
        Vy(j) =  Vy(j) + Accy*dt;
        x(j,2) = x(j,1); %updating previous positions
        y(j,2) = y(j,1);
        x(j,1) = x(j,1) + (dt * Vx(j));%updating new positions
        y(j,1) = y(j,1) + (dt * Vy(j));
        if x(j,1) > L %setting right side boundary condition
            x(j,2) = 0;
            x(j,1) = dt * Vx(j);
        end
        if x(j,1) < 0 %setting left side boundary condition
            x(j,2) = L;
            x(j,1) = x(j,2) + (dt * Vx(j));
        end 
        %setting roof/floor boundary condition
        if y(j,1) > W || y(j,1) < 0 
            Vy(j) = -Vy(j);
        end
        %setting up position line vectors to plot
        XPlot = [x(j,2) x(j,1)];
        YPlot = [y(j,2) y(j,1)];
        if j < VisibleEcount
            %plot visible position line vectors
        plot(XPlot,YPlot);
        end
        
        %sum of x and y  velocities
       VTotal = sqrt(Vx(j)^2 + Vy(j)^2);
       
        %sum of all temperatures
       Temperature(1,2) = Temperature(1,2) + (mn*Vtotal(j)^2)/(2*C.k);
         
    end
    %averaging temperature sum and plotting against each timestep
    Vdrift = mean(Vx(:));
    Ix(i) = L * W * q * Edensity * Vdrift;
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
%setup time axis for current plot
TimePlot = (0:1:Ecount-1)*dt;
figure(2)
plot(TimePlot, Ix)
xlabel('Time (s)')
ylabel('current')
title('Drift Current')

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
%uses histogram to model density with a 3rd demension of number of
%electrons. Would have ideally reindexed the position matrixes to account
%for all the different time locations, but didnt have time, would be the
%next thing to change
figure(3)
hist3([x(:,1) y(:,1)], [200 100])
title('3D Density map of Electrons')
xlabel('X (5nm)');
ylabel('Y (5nm)');
set(gca, 'Ydir', 'Normal')
c = colorbar;
title(c, 'Electron Count')

%Couldnt figure out how to make temperature 3d. thought i could plot it as
%the Z component but could not get the dimensions to match
% [Tx,Ty] = meshgrid(1:10,1:20);
% Tm = DensityMap(Tx, Ty)
figure(4)
% surf(Tx,Ty,Tm)
imagesc(TemperatureMap)
title('temperature map of Electrons');
xlabel('X (5nm)')
ylabel('Y (5nm)')
set(gca, 'Ydir', 'Normal')
c = colorbar;
title(c, 'Temperature (K)')
