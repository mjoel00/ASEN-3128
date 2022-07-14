% Matthew Pabin
% ASEN 3128
% Simulate Equations of Motion


tspan = [0 1];
y0 = [1 1 1];

[tout,yout] = ode45(@(tout,yout) odefun(yout),tspan,y0);

subplot(3,1,1)
plot(yout(:,1),tout)

subplot(3,1,2)
plot(yout(:,2),tout)

subplot(3,1,3)
plot(yout(:,3),tout)

function yout = odefun(yin)
x = yin(1);
y = yin(2);
z = yin(3);
yout(1) = x + (2*y) + z;
yout(2) = x - (5*z);
yout(3) = (x*y) - y^2 + 3*(z^2);

yout = yout';
end

