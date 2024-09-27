%% Midterm P1-1)
% DO NOT MODIFY THIS CELL
clear
clc
close all

load P1_b_j2k.mat

DTR = pi / 180;
RTD = 180 / pi;

J = [100 30 -10     % Moment of Inertia, kg m^2
    30 150 20
    -10 20 170];

m_max = 20; % Am^2, The maximum dipole moment of the torque rod

q0 = [0;0;0;1];
w0 = [5;5;5] * DTR;

X = zeros(7,N);
X(1:4,1) = q0;
X(5:7,1) = w0;
B_body = zeros(3, N);
m_mt = [0;0;0];
T_mag = [0;0;0];
for i = 2:1:N
    % Attitude Propagation
    X(:,i) = RK4(@eq_xdot, t(i), X(:,i-1), J, T_mag, dt);
    X(1:4,i) = X(1:4,i) / norm(X(1:4,i));

    % Geomagnetic Field Transformation to Body Frame
    q = X(1:4,i);
    w = X(5:7,i);
    T_j2k_to_body = quaternion_to_DCM(q);
    B_body(:,i) = T_j2k_to_body * B_j2k(:,i);
    
    % MT Command
    if(i > 2)
        m_mt = cal_mt_cmd_bdot(B_body(:,i), B_body(:,i-1), m_max, dt);
        
        % Checking the MT command validity
        if(norm(m_mt) > m_max*1.0001)
            error('ERROR! norm(m_mt) > m_max')
        end
    end
    
    T_mag = cross(m_mt, B_body(:,i)) * 1e-9;
end


w = X(5:7,:);
w_mag = zeros(1,N);
H = J * w;
H_mag = zeros(1,N);
for i = 1:1:N
    w_mag(1,i) = norm(w(:,i));
    H_mag(i) = norm(H(:,i));
end

figure
subplot(2,1,1)
plot(t, w_mag * RTD)
xlabel('t(sec)')
ylabel('w mag (deg/sec)')
grid on
title('P1-1 Bdot, Body Rate')

subplot(2,1,2)
for i = 1:1:3
    plot(t, w(i,:) * RTD)
    hold on
end
xlabel('t(sec)')
ylabel('w body (deg/sec)')
legend('x','y','z')
grid on


figure
subplot(2,1,1)
plot(t, H_mag)
xlabel('t(sec)')
ylabel('H mag (Nms)')
grid on
title('P1-1 Bdot, Angular Momentum')

subplot(2,1,2)
for i = 1:1:3
    plot(t, H(i,:))
    hold on
end
xlabel('t(sec)')
ylabel('H (Nms)')
legend('x','y','z')
grid on

%% TO DO: Complete the function
function m_mt = cal_mt_cmd_bdot(b_i, b_i_1, m_max, dt)
% INPUT
% b_i   : the geo-magnetic field vector in body frame at now (nT) [3 x 1]
% b_i_1 : the geo-magnetic field vector in body frame at the previous step
% (nT) [3 x 1]
% m_max : the maximum dipole moment of your magnetorquer (Am^2) (3D vector norm)
% [3 x 1]
% dt    : the time step between not and the previous step (sec)
%
% OUTPUT
% m_mt  : the magnetorquer dipole moment command (Am^2) [3 x 1]
md= (b_i_1-b_i)/dt;
m_mt=m_max/norm(md)*md;
end
%% DO NOT MODIFY THIS CELL
function A = quaternion_to_DCM(q)

A = [   q(1)^2-q(2)^2-q(3)^2+q(4)^2,    2 * (q(1) * q(2) + q(3) * q(4)),    2 * (q(1) * q(3) - q(2) * q(4));
        2 * (q(1) * q(2) - q(3) * q(4)),-q(1)^2 + q(2)^2 - q(3)^2 + q(4)^2, 2 * (q(2) * q(3) + q(1) * q(4));
        2 * (q(1) * q(3) + q(2) * q(4)),2 * (q(2) * q(3) - q(1) * q(4)),    -q(1)^2 - q(2)^2 + q(3)^2 + q(4)^2];
    
end

function [qdot] = quaternion_kinematics(q,w)
Omega = [0 w(3,1) -w(2,1) w(1,1);-w(3,1) 0 w(1,1) w(2,1); w(2,1) -w(1,1) 0 w(3,1); -w(1,1) -w(2,1) -w(3,1) 0];
qdot = 1/2*Omega * q;
end


function xdot = eq_xdot(t, x,J, T)
q = x(1:4,1);
w = x(5:7,1);
qdot = quaternion_kinematics(q, w);
wdot = inv(J) * (T - cross(w, J*w));
xdot = [qdot;wdot];
end

function [x] = RK4(Func, t, x,k,c, delt)

k1  = Func(t, x,k,c)*delt;
k2  = Func(t + 0.5*delt, x + k1*0.5,k,c)*delt;
k3  = Func(t + 0.5*delt, x + k2*0.5,k,c)*delt;
k4  = Func(t + delt, x + k3,k,c)*delt;

x = x + (k1 + 2*(k2 + k3) + k4)/6.0;
end