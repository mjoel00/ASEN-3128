function PlotAircraftSim(t,s,control, col,fig,Zc)


    figure( fig(1))    %Position
    subplot(3,1,1), plot(t,s(:,1),col); hold on;
    title('X Position (North +)');
    xlabel('Time [seconds]');
    ylabel('X Position [m]');
    subplot(3,1,2), plot(t,s(:,2),col); hold on;
    title('Y Position (East +)');
    xlabel('Time [seconds]');
    ylabel('Y Position [m]');
    subplot(3,1,3), plot(t,s(:,3),col); hold on;
    title('Z Position (Down +)');
    xlabel('Time [seconds]');
    ylabel('Z Position [m]');
    sgtitle(sprintf('Position of %s',control))
    
    figure(fig(2))    %Angular Position
    subplot(3,1,1), plot(t,s(:,4),col); hold on;
    title('Roll(\phi)');
    xlabel('Time [seconds]');
    ylabel('Roll(\phi) [rad]');
    subplot(3,1,2), plot(t,s(:,5),col); hold on;
    title('Pitch(\theta)');
    xlabel('Time [seconds]');
    ylabel('Pitch(\theta) [rad]');
    subplot(3,1,3), plot(t,s(:,6),col); hold on;
    title('Yaw(\psi)');
    xlabel('Time [seconds]');
    ylabel('Yaw(\psi)Position [rad]');
    sgtitle(sprintf('Angular Position of %s',control))
    
    figure(fig(3))    %Linear Velocity
    subplot(3,1,1), plot(t,s(:,7),col); hold on;
    title('u Inertial Velocity (North +)');
    xlabel('Time [seconds]');
    ylabel('u Inertial Velocity [m/s]');
    subplot(3,1,2), plot(t,s(:,8),col); hold on;
    title('v Velocity (East +)');
    xlabel('Time [seconds]');
    ylabel('v Inertial Velocity [m/s]');
    subplot(3,1,3), plot(t,s(:,9),col); hold on;
    title('w Velocity (Down +)');
    xlabel('Time [seconds]');
    ylabel('w Inertial Velocity [m/s]');
    sgtitle(sprintf('Linear Velocity of %s',control))
    
    figure(fig(4))    %Angular Velocity
    subplot(3,1,1), plot(t,s(:,10),col); hold on;
    title('p Angular Velocity');
    xlabel('Time [seconds]');
    ylabel('p Angular Velocity [rad/s]');
    subplot(3,1,2), plot(t,s(:,11),col); hold on;
    title('q Angular Velocity'); 
    xlabel('Time [seconds]');
    ylabel('q Angular Velocity [rad/s]');
    subplot(3,1,3), plot(t,s(:,12),col); hold on;
    title('r Angular Velocity');
    xlabel('Time [seconds]');
    ylabel('r Angular Velocity [rad/s]');
 sgtitle(sprintf('Angular Velocity of %s',control))
   
 
 %% Plot control forces and moments vs time
%     figure(fig(5))
%     subplot(2,2,1)
%     plot(t,Zc(:,1),col); hold on;
%     grid on
%     title([control ' Control Force vs. time']);
%     ylabel('Control Force');
%     xlim([t(1) t(end)]);
%     subplot(2,2,2)
%     plot(t,Zc(:,2),col,'linewidth',2); hold on;
%     grid on
%     title([control ' x Control Moment vs. time']);
%     ylabel('x Control Moment');
%     xlim([t(1) t(end)]);
%     subplot(2,2,3)
%     plot(t,Zc(:,3),col,'linewidth',2); hold on;
%     grid on
%     title([control ' y Control Moment vs. time']);
%     ylabel('y Control Moment');
%     xlabel('time [s]');
%     xlim([t(1) t(end)]);
%     subplot(2,2,4)
%     plot(t,Zc(:,4),col,'linewidth',2); hold on;
%     grid on
%     title([control ' z Control Moment vs. time']);
%     xlabel('time [s]');
%     ylabel('z Control Moment');
%     xlim([t(1) t(end)]);
    
    %% Plot flight path of the simulated quadrotor
    figure(fig(6))
    plot3(s(:,1),s(:,2),s(:,3),col,'linewidth',2); hold on;
    plot3(s(1,1),s(1,2),s(1,3),'g.','markersize',20); hold on;
    plot3(s(end,1),s(end,2),s(end,3),'r.','markersize',20); hold on;
    set(gca, 'YDir','reverse')
    set(gca, 'ZDir','reverse')
    grid on
    title([control ' Flight Path']);
    xlabel('x');
    ylabel('y');
    zlabel('z');

 
end