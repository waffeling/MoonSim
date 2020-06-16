%% Moon Sim (Doesn't Work But Makes Cool Graphs)
clear;
clc; 

%Details on the Rocket and Sim
G = 6.67e-11;
M = 7.34767309e22;
Rm = 1.737e6;
mnaught = 1500;
percentfuel = 0.1;
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
mnaught = 15000;
percentfuel = 0.3;
consumptionrate = 5; 
ejectvel = 3000;
step = 0.1;
timelength = 1000000;

%Details on the Environment

Rdot = [0];
ct = [zeros(1, timelength/step)];
cT = 1;
m = mnaught;
r = Rm+100000;


rdot = 0;
newrdot = 0;
R = [Rm+100000, zeros(1, (timelength/step)-1)];

phidot = 9.6e-4;
newphidot = 0;

phi = 0;
newphi = 0;
Phi = [0];

thetadot = 0;
newthetadot = 0;
theta = pi/2;
newtheta = pi/2;
Theta = [pi/2]; 

X = [Rm+100000, zeros(1, timelength/step - 1)];
Y = zeros(1, timelength/step);
Z = zeros(1, timelength/step);

%Calculate the Burn Time
burntime = (percentfuel * mnaught) / consumptionrate;


for t = step:step:timelength 
    
    gravity = (M*G)/(r^2);
    centripedala = r * phidot^2;
    thrust = ejectvel/ ((mnaught/consumptionrate) - t); 
    
    if t<= burntime
        rdoubledot = -gravity + centripedala + thrust;
        phidoubledot = (-2 * rdot * phidot)/r; 
        newphidot = phidot + (phidoubledot * step);
        newrdot = rdot + (rdoubledot *step);
        newphi = phi + newphidot * step;
        newr = r + newrdot * step; 
    
    else 
        rdoubledot = -gravity + centripedala; 
        phidoubledot = (-2 * rdot * phidot)/r; 
        newphidot = phidot + (phidoubledot * step);
        newrdot = rdot + (rdoubledot *step);
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
    
    R(cT) = newr;
    Phi(cT) = newphi;
    Theta(cT) = newtheta; 
    
    rdot = newrdot;
    r = newr;
    phidot = newphidot; 
    phi = newphi; 
    theta= newtheta;
    ct(cT) = t;
    cT =cT+1;
end

X = X(1:cT-1);
Y = Y(1:cT-1);
Z = Z(1:cT-1);
R = R(1:cT-1);
ct = ct(1:cT-1);
%Matrix Preallocation Correction


%Plotting
figure(1);
plot3(X, Y, Z);



figure(2);
plot(ct, R);
ylim([0, inf])

