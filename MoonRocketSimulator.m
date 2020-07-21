%% Moon Sim (Doesn't Work But Makes Cool Graphs)
clear;
clc; 

%Details on the Rocket and Sim
G = 6.67e-11;
M = 7.34767309e22;
Rm = 1.737e6;
mnaught = 4700;
percentfuel = 0.543;
consumptionrate = 2500; 
ejectvel = 3000;
step = 0.01;
timelength = 1000;

%Details on the Environment
R = [Rm+10];
Vr = [0];
ct = [0];
r = Rm+10;
newvr = 0;
vr = 0;
phidot = 0.5; 
phi = 0;
Phi = [0];
thetadot = 0;
theta = pi/2;
Theta = [pi/2]; 



%Calculate the Burn Time
burntime = (percentfuel * mnaught) / consumptionrate;


for t = step:step:1000
    %Burn
    if t <= burntime
       %Radial calculations
       newvdotr = (ejectvel / ((mnaught/consumptionrate) - t)) - ((G*M)/ r^2);
       newvr = vr + (newvdotr * step);
       newr = r + newvr * step;
       
       %Phi Calculations
       newphi = phi + (phidot * step);
       
       %Theta Calculations 
       newtheta = theta + (thetadot * step);
       
    
    %Non-Burn
    else
       newvdotr = -(G*M)/ r^2;
       newvr = vr + (newvdotr * step);
       newr = r + (newvr * step); 
       
       %Phi Calculations
       newphi = phi + (phidot * step);
       
       %Theta Calculations 
       newtheta = theta + (thetadot * step);
     
    end
    
    %Data Collection
    R = [R, newr];
    Phi = [Phi, newphi];
    Theta = [Theta, newtheta]; 
    
    vr = newvr;
    r = newr;
    phi = newphi; 
    theta= newtheta;
    ct = [ct, t];
    


end
%Coordinate Transformations


X = R .* sin(Theta) .* cos(Phi);
Y = R .* sin(Theta) .* sin(Phi);
Z = R .* cos(pi/2);

%Plotting
figure(1);
plot3(X, Y, Z);



figure(2);
plot(ct, R);
ylim([0, inf])

%% Moon Sim


clear;
clc; 

%Details on the Rocket and Sim
G = 6.67e-11;
M = 7.34767309e22;
Rm = 1.737e6;
mnaught = 4700;
percentfuel = 0.543;
consumptionrate = 8; 
ejectvel = 2000;
step = 0.01;
timelength = 10000;
m = mnaught; 
newm = 0;
Mass = [m, zeros(1, (timelength/step)-1)];
Gravity = [(G*M)/Rm^2, zeros(1, (timelength/step)-1)];
Ca = [zeros(1, timelength/step)];
totalke = 0;
KE = [zeros(1, timelength/step)];
Velocity = [zeros(1, timelength/step)];

%Details on the Environment

Rdot = [0];
ct = [zeros(1, timelength/step)];
cT = 1;
m = mnaught;
r = Rm;


rdot =0;
newrdot = 0;
R = [Rm, zeros(1, (timelength/step)-1)];
Rdoubledot= [0];

phidot = 2.4e-6;
newphidot = 0;
Phidot = [0];

phi = 0;
newphi = 0;
Phi = [0];

thetadot = 0;
newthetadot = 0;
theta = pi/2;
newtheta = pi/2;
Theta = [pi/2]; 

X = [Rm, zeros(1, timelength/step - 1)];
Y = zeros(1, timelength/step);
Z = zeros(1, timelength/step);
thrusttracker = zeros(1, timelength/step); 



%Calculate the Burn Time
burntime = (percentfuel * mnaught) / consumptionrate;

%Begin Simulation
for t = step:step:timelength 
    
    gravity = (M*G)/(r^2);
    centripedala = r * phidot^2;
    thrust = (consumptionrate * ejectvel)/m;
    velocity = sqrt(phidot^2 + rdot^2);
    kinetice = 0.5 * m * velocity^2;
    totalke = totalke + kinetice;
    vhatrad = rdot/sqrt(rdot^2 + (r^2 * phidot^2));
    vhatphi = (r*phidot)/sqrt(rdot^2 + (r^2 * phidot^2));
    
 
    if t < 438 && t <= burntime 

        rdoubledot = -gravity + centripedala + thrust;
        phidoubledot = (-2 * rdot * phidot)/r; 
        newphidot = phidot + (phidoubledot * step);
        newrdot = rdot + (rdoubledot * step);
        newphi = phi + newphidot * step;
        newr = r + newrdot * step; 
        newm = m - (consumptionrate * step);
    end
    
    if t < 438 && t >= burntime 

        rdoubledot = -gravity + centripedala;
        phidoubledot = (-2 * rdot * phidot)/r; 
        newphidot = phidot + (phidoubledot * step);
        newrdot = rdot + (rdoubledot * step);
        newphi = phi + newphidot * step;
        newr = r + newrdot * step; 
        newm = m - (consumptionrate * step);
    end
  
    if t >= 438 && t < 4013
        
        rdoubledot = -gravity + centripedala;
        phidoubledot = (-2 * rdot * phidot)/r; 
        newphidot = phidot + (phidoubledot * step);
        newrdot = rdot + (rdoubledot * step);
        newphi = phi + newphidot * step;
        newr = r + newrdot * step;
    
    end
        
    if  t >= 4013 && t < 4058 

        rdoubledot = -gravity + centripedala;
        phidoubledot = (-2 * rdot * phidot)/r ; 
        newphidot = phidot + (phidoubledot * step);
        newrdot = rdot + (rdoubledot * step);
        newphi = phi + newphidot * step;
        newr = r + newrdot * step; 
        
   
    end 
    
   
    
    if  t >= 4058 
        
        rdoubledot = -gravity + centripedala;
        phidoubledot = (-2 * rdot * phidot)/r; 
        newphidot = phidot + (phidoubledot * step);
        newrdot = rdot + (rdoubledot * step);
        newphi = phi + newphidot * step;
        newr = r + newrdot * step; 
        
        
    end
    
    if newr < Rm 
        break
    end
    
    newx = newr * sin(newtheta) * cos(newphi);
    newy = newr * sin(newtheta) * sin(newphi);
    newz = 0;
    
    X(cT) = newx;
    Y(cT) = newy;
    Z(cT) = newz;
    
    R(cT) = newr - Rm;
    Phi(cT) = newphi;
    Theta(cT) = newtheta; 
    thrusttracker(cT) = thrust;
    Mass(cT) = m;
    Rdoubledot(cT) = rdoubledot;
    Gravity(cT) = gravity;
    Ca(cT) = centripedala;
    Phidot(cT) = phidot;
    KE(cT) = totalke;
    Velocity(cT) = velocity;
    
    rdot = newrdot;
    r = newr;
    phidot = newphidot; 
    phi = newphi; 
    theta= newtheta;
    m = newm;
    ct(cT) = t;
    cT =cT+1;
end

%Some Data Collection
[MaxV, tforMaxV] = max(Velocity);
[MaxR, tforMaxR] = max(R);


%Matrix Preallocation Correction
X = X(1:cT-1);
Y = Y(1:cT-1);
Z = Z(1:cT-1);
R = R(1:cT-1);
ct = ct(1:cT-1);

%Plotting
[X1, Y1, Z1] = sphere;
X2 = X1*Rm;
Y2 = Y1*Rm;
Z2 = Z1*Rm;

figure(1);
Plot = plot3(X, Y, Z);
hold on
Moon = mesh(X2, Y2, Z2);
hold off

figure(2);
plot(ct, R);
ylim([0, inf]);


