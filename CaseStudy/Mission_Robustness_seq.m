function [negative_rob,xx,yy,zz] = Mission_Robustness_seq(var,optParams,d,w_opt_rob)
%%
import casadi.*
type_of = isfloat(var); %0 for casadi
optParams.type_of = type_of;
H = optParams.H_formula; %better be even
w = var(1:numel(var)/2);
v = var(numel(var)/2+1:end);

% assign intervals for goals
N_per_T = optParams.N_per_T;
% 50 works for 2 drones, 2 works for rural mission
%C = 30.0; %const for smooth min/max operation %for 2 drones, use 10 with a period of 5s, casadi is unstable numerically
C1 =30.0; %const for smooth max %20 works for 10 drones, det init points
C2 =30.0; %5 makes no numerical instblty in 12 drones
%optParams.C = C;
C = optParams.C;
optParams.C1 = C1;
optParams.C2 = C2;

Nobs = size(optParams.obs,1); %number of obstacles
dT = optParams.T/optParams.N_per_T:optParams.T/optParams.N_per_T:optParams.T;
dv = optParams.dv;
da = optParams.da;
M1 = optParams.M1;
T = optParams.T;
Clen = optParams.Clen;

if(type_of) %if double input
    temp_x = zeros(1, numel(dT));
    temp_y = zeros(1, numel(dT));
    temp_z = zeros(1, numel(dT));
    xx = zeros(numel(dT)+1,1);
    yy = zeros(numel(dT)+1,1);
    zz = zeros(numel(dT)+1,1);
    rho_unsafe = zeros(1,1);
    rho_goal = zeros(1,1);
    mutual_distances = zeros(numel(dT)*optParams.H_formula+1,1);
    if (d > 1)
        dists = zeros(d-1,1);
    end
else
    temp_x = MX.zeros(1, numel(dT));
    temp_y = MX.zeros(1, numel(dT));
    temp_z = MX.zeros(1, numel(dT));
    xx = MX.zeros(numel(dT)*optParams.H_formula+1,1);
    yy = MX.zeros(numel(dT)*optParams.H_formula+1,1);
    zz = MX.zeros(numel(dT)*optParams.H_formula+1,1);
    rho_unsafe = MX.zeros(1,1);
    rho_goal = MX.zeros(1,1);
    mutual_distances = MX.zeros(numel(dT)*optParams.H_formula+1,1);
    if (d > 1)
        dists = MX.zeros(d-1,1);
    end
end

% For each Drone
%for d = 1:optParams.N_drones
    %init posns
    xx(1,1) = w(1);
    yy(1,1) = w(2);
    zz(1,1) = w(3);
    
    %get all sampled splines
    for k = 1:optParams.H_formula
        
        % dp for all axes
        dp_x = w(k*3+1) - w((k-1)*3+1) - T*v((k-1)*3+1);
        dp_y = w(k*3+2) - w((k-1)*3+2) - T*v((k-1)*3+2);
        dp_z = w(k*3+3) - w((k-1)*3+3) - T*v((k-1)*3+3);
        
        % constants for all 3 axes
        al_x = M1(1,:)*[dp_x;dv;da];
        be_x = M1(2,:)*[dp_x;dv;da];
        gam_x = M1(3,:)*[dp_x;dv;da];
        al_y = M1(1,:)*[dp_y;dv;da];
        be_y = M1(2,:)*[dp_y;dv;da];
        gam_y = M1(3,:)*[dp_y;dv;da];
        al_z = M1(1,:)*[dp_z;dv;da];
        be_z = M1(2,:)*[dp_z;dv;da];
        gam_z = M1(3,:)*[dp_z;dv;da];
        
        temp_x(1,:) = (al_x/120)*dT.^5 + (be_x/24)*dT.^4 + ...
            (gam_x/6)*dT.^3   + w((k-1)*3+1) + ...
            dT*v((k-1)*3+1); %fix w points
        
        xx(2+(k-1)*N_per_T:k*N_per_T+1,1) = temp_x(1,:)';
        
        temp_y(1,:) = (al_y/120)*dT.^5 + (be_y/24)*dT.^4 + ...
            (gam_y/6)*dT.^3   + w((k-1)*3+2) + ...
            dT*v((k-1)*3+2); % fix w points
        
        yy(2+(k-1)*N_per_T:k*N_per_T+1,1) = temp_y(1,:)';
        
        temp_z(1,:) = (al_z/120)*dT.^5 + (be_z/24)*dT.^4 + ...
            (gam_z/6)*dT.^3   + w((k-1)*3+3) + ...
            dT*v((k-1)*3+3); %fix w points
        
        zz(2+(k-1)*N_per_T:k*N_per_T+1,1) = temp_z(1,:)'; 
        
    end
    
    % Robusteness unsafe set
    %rho_unsafe(d) = robustness_unsafe(xx,yy,zz,d,optParams);
    rho_unsafe(1) = robustness_unsafe_seq(xx,yy,zz,d,optParams);
    
    
    % Robustness Get to goal in interval
    i = 1;
    rho = [];
    drone_goals = optParams.drone_goals{d};
    if size(drone_goals)
        for g = drone_goals(:,1)'
            I = 1+drone_goals(i,2)*N_per_T/T:1+drone_goals(i,3)*N_per_T/T;
            rho = [rho; robustness_goal_seq(xx,yy,zz,d,g,I,optParams)];
            i = i + 1;
        end
    end
    rho_goal(1) = SmoothMin(rho,C);
    
%     negative_rob = -SmoothMin([rho_unsafe;rho_goal],C);

%% pairwise distances
if (d > 1)
    for p = 1:(d-1)
        w_opt=w_opt_rob(:,p);
        [w_opt_xx,w_opt_yy,w_opt_zz] = Mission_Robustness_exact_seq(w_opt,optParams);
        for k=1:size(xx,1) %for all time steps
            % this drone
            pa = [xx(k,1);yy(k,1);zz(k,1)];
            % other drone
            pb = [w_opt_xx(k,1);w_opt_yy(k,1);w_opt_zz(k,1)];
            mutual_distances(k) = norm(pa-pb,2)-optParams.d_min;
        end
        %minimum distance
        dists(p) = SmoothMin(mutual_distances,C2);
    end
end

if (d > 1)
    negative_rob = -SmoothMin([rho_unsafe;rho_goal;dists],C);
else
    negative_rob = -SmoothMin([rho_unsafe;rho_goal],C);
end
