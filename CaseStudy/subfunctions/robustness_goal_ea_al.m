function [ rho_goal ] = robustness_goal_ea_al(xx,yy,zz,d,g,I,optParams)
% with respect to goal g
% 
import casadi.*
I1=I{1,2};
I2=I{2,2};

type_of = optParams.type_of;
if(type_of)
        %temp = zeros(numel(I),1);
        temp_e = zeros(numel(I1),1);
        temp_a = zeros(numel(I2),1);
        temp_unsafe = zeros(size(optParams.obs,1),1);
    else
%         temp = MX.sym('temp',numel(I),1);
%         temp_unsafe = MX.sym('temp_unsafe',size(optParams.obs,1),1);
          %temp = MX.zeros(numel(I),1);
          temp_e = MX.zeros(numel(I1),1);
          temp_a = MX.zeros(numel(I2),1);
          temp_unsafe = MX.zeros(size(optParams.obs,1),1);

end


C1 = optParams.C1;
C = optParams.C;

% eventually goal in x y z
    rho_lb_xx = xx(I1,d)-optParams.goal{g}.goal_lb_N(I1,1);
    rho_ub_xx = optParams.goal{g}.goal_ub_N(I1,1)-xx(I1,d);
    rho_lb_yy = yy(I1,d)-optParams.goal{g}.goal_lb_N(I1,2);
    rho_ub_yy = optParams.goal{g}.goal_ub_N(I1,2)-yy(I1,d);
    rho_lb_zz = zz(I1,d)-optParams.goal{g}.goal_lb_N(I1,3);
    rho_ub_zz = optParams.goal{g}.goal_ub_N(I1,3)-zz(I1,d);
    
    
    % make this more efficient
    for i = 1:numel(rho_lb_xx)
        temp_vec = [rho_lb_xx(i) rho_ub_xx(i) rho_lb_yy(i) rho_ub_yy(i) ...
            rho_lb_zz(i) rho_ub_zz(i)];
        temp_e(i) = SmoothMin(temp_vec,C);
        
    end
    rho_goal = SmoothMax(temp_e,C1);
    
% always goal in x y z
    rho_lb_xx = xx(I2,d)-optParams.goal{g}.goal_lb_N(I2,1);
    rho_ub_xx = optParams.goal{g}.goal_ub_N(I2,1)-xx(I2,d);
    rho_lb_yy = yy(I2,d)-optParams.goal{g}.goal_lb_N(I2,2);
    rho_ub_yy = optParams.goal{g}.goal_ub_N(I2,2)-yy(I2,d);
    rho_lb_zz = zz(I2,d)-optParams.goal{g}.goal_lb_N(I2,3);
    rho_ub_zz = optParams.goal{g}.goal_ub_N(I2,3)-zz(I2,d);
    
    
    % make this more efficient
    for i = 1:numel(rho_lb_xx)
        temp_vec = [rho_lb_xx(i) rho_ub_xx(i) rho_lb_yy(i) rho_ub_yy(i) ...
            rho_lb_zz(i) rho_ub_zz(i)];
        temp_a(i) = SmoothMin(temp_vec,C);
        
    end
    rho_goal = [rho_goal;SmoothMin(temp_a,C1)];



end

