%% IDM_simulation.m
% Intelligent Driver Model (IDM) demo: lead + ego vehicle longitudinal simulation
% Saves plots in ./results/
% Author: <Your Name> -- ready for GitHub

clear; close all; clc;

% Create results folder
if ~exist('results','dir'), mkdir('results'); end

%% Simulation parameters
dt = 0.05;            % time step (s)
Tsim = 60;            % total simulation time (s)
time = 0:dt:Tsim;

% Vehicle geometry
car_length = 4.5;     % effective vehicle length (m) (lead and ego same for gap calc)

% Lead vehicle trajectory (you can modify for scenarios)
% Example: lead drives at constant speed, then brakes at t=20 for 5s, then recovers
v_lead0 = 20;         % initial lead speed (m/s)
v_lead = v_lead0 * ones(size(time));

% make a braking maneuver from t=20 to t=25:
t_brake_start = 20;
t_brake_end   = 25;
for k=1:length(time)
    t = time(k);
    if t >= t_brake_start && t <= t_brake_end
        % linear decel to slow down by 8 m/s over 5s (~-1.6 m/s^2)
        frac = (t - t_brake_start)/(t_brake_end - t_brake_start);
        v_lead(k) = v_lead0 - 8*frac;
    elseif t > t_brake_end
        v_lead(k) = v_lead0 - 8; % slowed down permanently
    end
end

% Integrate lead positions
x_lead = zeros(size(time));
for k=2:length(time)
    x_lead(k) = x_lead(k-1) + v_lead(k-1)*dt;
end

%% IDM parameters (default / baseline)
params.v0    = 30;      % desired speed (m/s)
params.T     = 1.2;     % desired time headway (s)
params.s0    = 2;       % minimum gap (m)
params.a_max = 1.0;     % max acceleration (m/s^2)
params.b     = 2.0;     % comfortable deceleration (m/s^2)
params.delta = 4;       % exponent

% Initial conditions for ego
x0_ego = -20;           % start 20 m behind lead
v0_ego = 20;            % initial ego speed (m/s)

% Pre-allocate
N = length(time);
x_ego = zeros(1,N);
v_ego = zeros(1,N);
a_ego = zeros(1,N);
gap   = zeros(1,N);
TTC   = inf(1,N);

x_ego(1) = x0_ego;
v_ego(1) = v0_ego;

%% Helper: IDM acceleration as function handle
idm_acc = @(v, s, dv, p) p.a_max*( 1 - (v./p.v0).^p.delta - ( ( p.s0 + v.*p.T + (v.*dv)./(2*sqrt(p.a_max*p.b)) )./s ).^2 );

%% Simulation loop (baseline)
for k = 1:N-1
    s_net = (x_lead(k) - x_ego(k) - car_length); % net gap
    if s_net < 0.01
        s_net = 0.01; % prevent division by zero / collisions numerically
    end
    dv = v_ego(k) - v_lead(k); % relative speed (positive if ego faster)
    a_now = idm_acc(v_ego(k), s_net, dv, params);
    
    % clip unrealistic accelerations
    a_now = max(a_now, -9.0); % don't allow >9 m/s2 braking for stability
    a_now = min(a_now, 5.0);  % limit accel
    
    % Euler integration
    v_ego(k+1) = max(0, v_ego(k) + a_now*dt);
    x_ego(k+1) = x_ego(k) + v_ego(k)*dt + 0.5*a_now*dt^2;
    a_ego(k)    = a_now;
    gap(k)      = s_net;
    
    % TTC estimate
    if dv > 0
        TTC(k) = s_net / dv;
    else
        TTC(k) = inf;
    end
end
a_ego(end) = a_ego(end-1);
gap(end)   = (x_lead(end) - x_ego(end) - car_length);
TTC(end)   = TTC(end-1);

%% Plotting: baseline case
figure('Units','normalized','Position',[0.1 0.1 0.6 0.7]);
subplot(4,1,1); plot(time, x_lead, 'b', 'LineWidth', 1.5); hold on; plot(time, x_ego, 'r', 'LineWidth', 1.2);
ylabel('Position (m)'); legend('Lead','Ego'); title('Positions');
subplot(4,1,2); plot(time, v_lead, 'b', 'LineWidth', 1.5); hold on; plot(time, v_ego, 'r', 'LineWidth', 1.2);
ylabel('Velocity (m/s)'); legend('Lead','Ego'); title('Velocities');
subplot(4,1,3); plot(time, a_ego, 'r', 'LineWidth', 1.2); ylabel('Acceleration (m/s^2)'); title('Ego Acceleration');
subplot(4,1,4); plot(time, gap, 'k', 'LineWidth', 1.2); hold on; plot(time, TTC, '--'); ylabel('Gap (m) / TTC (s)'); xlabel('Time (s)');
legend('Gap','TTC'); title('Gap and TTC');

% Save figure
saveas(gcf, fullfile('results','baseline_plots.png'));

%% Parameter sensitivity study (varying T and b)
T_values = [0.8, 1.2, 1.8];        % time headways
b_values = [1.0, 2.0, 3.5];        % comfortable decels

% record min gaps and min TTC
minGap = zeros(length(T_values), length(b_values));
minTTC = zeros(length(T_values), length(b_values));
collision = false(length(T_values), length(b_values));

for i=1:length(T_values)
    for j=1:length(b_values)
        p = params;
        p.T = T_values(i);
        p.b = b_values(j);
        
        % reset ego initial
        xe = x0_ego;
        ve = v0_ego;
        
        min_gap_local = inf;
        min_ttc_local = inf;
        collided = false;
        
        for k=1:N-1
            s_net = (x_lead(k) - xe - car_length);
            if s_net <= 0
                collided = true;
                s_net = 0.01;
            end
            dv = ve - v_lead(k);
            a_now = idm_acc(ve, s_net, dv, p);
            a_now = max(a_now,-9.0); a_now = min(a_now,5.0);
            % integrate
            ve = max(0, ve + a_now*dt);
            xe = xe + ve*dt + 0.5*a_now*dt^2;
            % statistics
            min_gap_local = min(min_gap_local, s_net);
            if dv > 0
                ttc_k = s_net/dv;
                min_ttc_local = min(min_ttc_local, ttc_k);
            end
        end
        
        minGap(i,j) = min_gap_local;
        minTTC(i,j) = min_ttc_local;
        collision(i,j) = collided;
    end
end

% Display sensitivity table (text)
fprintf('\nParameter sweep results (min gap in meters):\n');
for i=1:length(T_values)
    fprintf('T = %.2f: ', T_values(i));
    fprintf('%.2f  ', minGap(i,:));
    fprintf('\n');
end
fprintf('\nCollision occurred (1=yes):\n');
disp(collision);

% Plot heatmap of minGap
figure;
imagesc(b_values, T_values, minGap);
colorbar; xlabel('b (m/s^2)'); ylabel('T (s)');
title('Min gap during simulation (m)');
set(gca,'YDir','normal'); % keep y low->high upward
saveas(gcf, fullfile('results','sensitivity_mingap.png'));

% Plot collision flags
figure;
imagesc(b_values, T_values, double(collision));
colormap(parula);
colorbar; xlabel('b (m/s^2)'); ylabel('T (s)');
title('Collision flag (1 = collision observed)');
set(gca,'YDir','normal');
saveas(gcf, fullfile('results','sensitivity_collision.png'));

%% Save important data
save(fullfile('results','simulation_data.mat'), 'time', 'x_lead','v_lead','x_ego','v_ego','a_ego','gap','T_values','b_values','minGap','minTTC','collision','params');

disp('Simulation complete. Results saved to ./results/');
