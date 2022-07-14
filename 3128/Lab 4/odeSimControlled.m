function x = odeSimControlled(t,s,C)
g = C{1}; m = C{2}; R = C{3}; km = C{4}; Ix = C{5}; Iy = C{6}; Iz = C{7}; vCoef = C{8}; mu = C{9}; k1_lat = C{10}; k2_lat = C{11}; k1_long = C{12}; k2_long = C{13};
   
   
   
   %S Vector x y z phi theta psi u v w p q r
   
   [X,Y,Z] = aeroDragForce(vCoef,s(7),s(8),s(9));
   [L,M,N] = aeroMoments(s(10),s(11),s(12),mu);
   
   MomentControl = [-k1_lat*s(10)-k2_lat*s(4);-k1_long*s(11)-k2_long*s(5); -0.004*s(12)]; % control moment
 
   
   Zc = controlZc(m,g); 
   
   eq_1 = firstEq(s(4),s(5),s(6),s(7),s(8),s(9));
   eq_2 = secondEq(s(4),s(5),s(6),s(10),s(11),s(12));
    
   eq_3 = thirdEq(s(4),s(5),s(6),s(10),s(11),s(12),s(7),s(8),s(9),X,Y,Z,Zc,g,m);   
   eq_4 = fourthEq(Ix,Iy,Iz,s(10),s(11),s(12),L,M,N,MomentControl(1),MomentControl(2),MomentControl(3)); 
   
   
   deltaF = [0 ; -0.004 * eq_4];
    
   
   

    
   x = [eq_1;eq_2;eq_3;eq_4];
end

function Zc = controlZc(m,g)
    Zc = m*g;
end

function [X,Y,Z] = aeroDragForce(vCoef,u,v,w)
    X = -vCoef*sqrt(u^2+v^2+w^2)*u;
    Y = -vCoef*sqrt(u^2+v^2+w^2)*v;
    Z = -vCoef*sqrt(u^2+v^2+w^2)*w;
    %out = vCoef*(u^2+v^2+w^2)*[u;v;w];
end

function [L,M,N] = aeroMoments(p,q,r,mu)
    L = -mu*(sqrt(p^2+q^2+r^2))*p;
    M = -mu*(sqrt(p^2+q^2+r^2))*q;
    N = -mu*(sqrt(p^2+q^2+r^2))*r;
    %out = -mu*(sqrt(p^2+q^2+r^2))*[p;q;r];    
end

function out = firstEq(phi,theta,psi,uE,vE,wE)
%Linear Velocities

    velocities = [uE;vE;wE];
    
    A = [cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);...
        cos(theta)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);...
        -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta);];
   
    out = A*[uE;vE;wE];
end

function out = secondEq(phi,theta,psi,p,q,r)
%Angular Rates
    A = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);...
        0 cos(phi) -sin(phi);...
        0 sin(phi)*sec(theta) cos(phi)*sec(theta)];
    out = A*[p;q;r];
end

function out = thirdEq (phi,theta,psi,p,q,r,uE,vE,wE,X,Y,Z,Zc,g,m)
%Acceleration
    A = [r*vE - q*wE; p*wE - r*uE; q*uE - p*vE]; %Angular Part
    B = g*[-sin(theta); cos(theta)*sin(phi); -cos(theta)*cos(phi)]; %Gravity Part
    C = [X; Y; Z]/m;                     %Aero Force Part
    D = [0; 0; Zc]/m;   %Control Force Part
    out = A + B + C + D;
%     out(3,1) = 0;
end

function out = fourthEq(Ix,Iy,Iz,p,q,r,L,M,N,Lc,Mc,Nc)
%Angular Accel. w/ moments
    A = [((Iy-Iz)/Ix)*q*r; ((Iz-Ix)/Iy)*p*r; ((Ix-Iy)/Iz)*p*q];
    B = [L/Ix; M/Iy; N/Iz];
    C = [Lc/Ix; Mc/Iy; Nc/Iz];
    out = A + B + C;
    
end