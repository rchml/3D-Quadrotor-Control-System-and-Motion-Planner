function [u1, u2] = controller(t, state, des_state, params)

%   state: The current state of the robot with the following fields:
%   state.pos = [x; y; z], 
%   state.vel = [x_dot; y_dot; z_dot],
%   state.rot = [phi; theta; psi], 
%   state.omega = [p; q; r]
%
%   des_state: The desired states are:
%   des_state.pos = [x; y; z], 
%   des_state.vel = [x_dot; y_dot; z_dot],
%   des_state.acc = [x_ddot; y_ddot; z_ddot], 
%   des_state.yaw,
%   des_state.yawdot
%
%   params: robot parameters


% -------------------- |         global params       | --------------------

m = params.mass;
g = params.gravity;

% State vars
x = state.pos(1);
y = state.pos(2);
z = state.pos(3);

x_vel = state.vel(1);
y_vel = state.vel(2);
z_vel = state.vel(3);

phi = state.rot(1);
theta = state.rot(2);
psi = state.rot(3);

p = state.omega(1);
q = state.omega(2);
r = state.omega(3);

% Desired vars
x_des = des_state.pos(1);
y_des = des_state.pos(2);
z_des = des_state.pos(3);

x_vel_des = des_state.vel(1);
y_vel_des = des_state.vel(2);
z_vel_des = des_state.vel(3);

x_acc_des = des_state.acc(1);
y_acc_des = des_state.acc(2);
z_acc_des = des_state.acc(3);

yaw_des = des_state.yaw;
yaw_dot_des = des_state.yawdot;

% p,q,r, psi desired
p_des = 0;
q_des = 0;
psi_des = yaw_des;
r_des = yaw_dot_des;


% ------------------------- |         U1       | --------------------------


kpz = 800;
kdz = 30;

kpx = 32;
kdx = 2.4;

kpy = 32;
kdy = 3.2;

u1 = m * (g + z_acc_des + ( kdz * ( z_vel - z_vel )) +  kpz * ( z_des - z ) );


% ----------------------- |         U2       | -------------------------

% U2 Phi
Kpphi = 160;
Kdphi = 2;

% U2 theta
Kptheta = 160;
Kdtheta = 2;

% U2 psi
Kppsi = 160;
Kdpsi = 2;


m1 = x_acc_des + kdx * (x_vel_des - x_vel) + kpx * (x_des - x);
m2 = y_acc_des + kdy * (y_vel_des - y_vel) + kpy * (y_des - y);

% psi desired
phi_des = (1/g) * ( m1 * yaw_des - m2  );

% theta desired
theta_des = (1/g) * ( m1 +  m2 * yaw_des );

u2_1 = Kpphi * (phi_des - phi)          + Kdphi * (p_des - p);
u2_2 = Kptheta * (theta_des - theta)    + Kdtheta * (q_des - q );
u2_3 = Kppsi * (yaw_des - psi)          + Kdpsi * (yaw_dot_des - r);


u2 = [u2_1;u2_2;u2_3];



end
