k1 = 20;
k2 = 30;

kt = 1/(1/k1 + 1/k2);

sigma = k2/k1;

r1 = 15;
r2 = 10;

kappa = r1/r2;

d = 60;
x_zero = d - r1 - r2;

x1_zero = x_zero * sigma / (1 + sigma);
x2_zero = x_zero - x1_zero;

vec_2 = [0 -d];
vec_1 = [0 0];

ang_vel1 = pi/4;
ang_vel2 = ang_vel1*kappa;

syms theta
%eqn = (acos((d-r2*cos(theta*kappa) - d)^2 - r2^2 - ((r2*sin(kappa*theta) + r1*sin(theta))^2 + (r2*cos(kappa*theta) - d + r1*cos(theta))^2)))/(-2*r2*((r2*sin(kappa*theta) + r1*sin(theta))^2 + (r2*cos(kappa*theta) - d + r1*cos(theta))^2)^(1/2)) == 0;
eqn = atan((-r2*cos(kappa*theta) + d - r1*cos(theta)) / (r2*sin(kappa*theta) + r1*sin(theta))) + pi/2 == theta*kappa;
endtheta = solve(eqn, theta);

theta1 = 0:0.01:endtheta;
theta2 = theta1*kappa;
time = theta1/ang_vel1;

x1 = -r1 * sin(theta1);
y1 = -r1 * cos(theta1);

x2 = r2 * sin(theta2);
y2 = r2 * cos(theta2) - d;

xm = x1 + (x2-x1)*sigma/(1 + sigma);
ym = y1 + (y2-y1)*sigma/(1 + sigma);

newtheta = endtheta*kappa + 0.01;
newtime = newtheta/ang_vel2;
time_diff = 0.02;
ang_vel = 0;
vel_damping = 5;
moment_of_inertia = 25;

x1_n = -r1 * sin(endtheta + 0.01/kappa);
y1_n = -r1 * cos(endtheta + 0.01/kappa);

p1 = [x1_n y1_n];

disp(double(endtheta*180/pi));
disp(p1);
disp(p2);

x2_n = r2 * sin(newtheta);
y2_n = r2 * cos(newtheta) - d;
    
p2 = [x2_n y2_n];

disp_thetas = [];

mass_vel = 0;
mass = 0.5;
mass_damping = 3;

spring_vec = p1 - p2;

spring1_x = norm(spring_vec) * sigma / (sigma + 1);
spring2_x = norm(spring_vec) - spring1_x;

dispx = [];
dispt = [];
dispx_len = [];

dispxm = [];
dispym = [];

dispp2x = [];
dispp2y = [];

simlen = 3;

for t = newtime:time_diff:newtime+simlen;
    newtheta = newtheta + ang_vel*time_diff;

    x2_n = r2 * sin(newtheta);
    y2_n = r2 * cos(newtheta) - d;
    
    p2 = [x2_n y2_n];
    
    spring_vec = p1 - p2;
    
    force = (sqrt(sum(spring_vec .* spring_vec)) - x_zero) * kt;
    force_vec = (spring_vec/norm(spring_vec)) * force;
    
    r2_vec = p2 - vec_2;
    
    r2_unit = r2_vec/norm(r2_vec);
    moment_arm = force_vec - dot(force_vec, r2_unit) * r2_unit;
    
    moment_arm_3d = [moment_arm 0];
    r2_vec_3d = [r2_vec 0];
    
    moment_vec = cross(moment_arm_3d, r2_vec_3d);
    
    moment = moment_vec(3);
    acceleration = moment / moment_of_inertia;
    delta_v = acceleration * time_diff;
    ang_vel = ang_vel + delta_v;
    
    ang_vel = ang_vel - ang_vel * vel_damping * time_diff;
    newtime = newtime + time_diff;
    
    disp_thetas = [disp_thetas double(newtheta * 180 / pi)];
    
    % Mass spring system
    spring1_x = spring1_x - mass_vel * time_diff;
    spring2_x = norm(spring_vec) - spring1_x;
    
    force1 = (spring1_x - x1_zero) * k1;
    force2 = (spring2_x - x2_zero) * k2;
    
    net_force = force1 - force2;
    mass_accel = net_force / mass;
    mass_delta_v = mass_accel * time_diff;
    mass_vel = mass_vel + mass_delta_v;
    
    mass_vel = mass_vel - mass_vel * mass_damping * time_diff;
    
    dispx = [dispx spring1_x];
    dispt = [dispt t];
    
    length = spring1_x + spring2_x;
    dispx_len = [dispx_len length];
    
    ratio = spring1_x / (spring1_x + spring2_x);
    
    m_pos = p1 - spring_vec * ratio;
    
    %disp(p1)
    %disp(spring_vec);
    %disp(ratio);
    %disp(m_pos);
    
    dispxm = [dispxm m_pos(1)];
    dispym = [dispym m_pos(2)];
    
    dispp2x = [dispp2x p2(1)];
    dispp2y = [dispp2y p2(2)];
end

subplot(2,2,1);
plot(dispt, dispx);

subplot(2,2,2);
plot(dispt, dispx_len);

subplot(2,2,3);
dispp2x = [x2 dispp2x];
dispp2y = [y2 dispp2y];
%plot(dispp2x, dispp2y);

% catching up
theta1_c = endtheta:0.01:newtheta/kappa;
time_c = theta1_c/ang_vel1 + simlen;

dispt = [time dispt time_c];

x1_c = -r1 * sin(theta1_c);
y1_c = -r1 * cos(theta1_c);

numelem = size(x1_c);
numelem = numelem(2);

x2_c = x2_n * ones(1, numelem);
y2_c = y2_n * ones(1, numelem);

xm_c = x1_c + (x2_c-x1_c)*sigma/(1 + sigma);
ym_c = y1_c + (y2_c-y1_c)*sigma/(1 + sigma);

dispxm = [xm dispxm xm_c];
dispym = [ym dispym ym_c];

plot(dispt, dispxm);

subplot(2,2,4);
plot(dispt, dispym);

%{
N = NaN(size(dispxm));

h = plot(N, N);
xdata = get(h, 'XData');
ydata = get(h, 'YData');

axis([-13 13 -45 0])
numelem = size(dispxm);
numelem = numelem(2);

for k = 1:numelem
    xdata(k) = dispxm(k);
    ydata(k) = dispym(k);
    % Update the plot
    set(h, 'XData', xdata, 'YData', ydata);
    drawnow
end
%}


