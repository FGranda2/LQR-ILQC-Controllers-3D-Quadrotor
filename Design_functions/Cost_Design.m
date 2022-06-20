% Cost_Design: Definition of cost functions for reaching goal state and/or
% passing through via-point.
%
% Control for Robotics
% AER1517 Spring 2022
% Assignment 2
%
% --
% University of Toronto Institute for Aerospace Studies
% Dynamic Systems Lab
%
% Course Instructor:
% Angela Schoellig
% schoellig@utias.utoronto.ca
%
% Teaching Assistant: 
% SiQi Zhou
% siqi.zhou@robotics.utias.utoronto.ca
% Lukas Brunke
% lukas.brunke@robotics.utias.utoronto.ca
% Adam Hall
% adam.hall@robotics.utias.utoronto.ca
%
% This script is adapted from the course on Optimal & Learning Control for
% Autonomous Robots at the Swiss Federal Institute of Technology in Zurich
% (ETH Zurich). Course Instructor: Jonas Buchli. Course Webpage:
% http://www.adrlab.org/doku.php/adrl:education:lecture:fs2015
%
% --
% Revision history
% [20.01.31]    first version

function Cost = Cost_Design( m_quad, Task)
%COST_DESIGN Creates a cost function for LQR and for ILQC
%   .Q_lqr 
%   .R_lqr
%
%   A ILQC cost J = h(x) + sum(l(x,u)) is defined by:
%   .h   - continuous-time terminal cost
%   .l - continous-time intermediate cost

% Quadcopter system state x
syms qxQ qyQ qzQ qph qth qps dqxQ dqyQ dqzQ dqph dqth dqps real
x_sym = [ qxQ qyQ qzQ ...    % position x,y,z
          qph qth qps ...    % roll, pitch, yaw
          dqxQ dqyQ dqzQ ... % velocity x,y,z
          dqph dqth dqps ]'; % angular velocity roll, pitch, yaw   
      
% Quadcopter input u (Forces / Torques)
syms Fz Mx My Mz real; u_sym = [ Fz Mx My Mz ]';

% Time variable for time-varying cost
syms t_sym real;
Cost = struct;


%% LQR cost function
% LQR controller
% Q needs to be symmetric and positive semi-definite
Cost.Q_lqr = diag([  1   1    1   ...   % penalize positions
                     3   3    3   ...   % penalize orientations
                     0.1  0.1    2   ...   % penalize linear velocities
                     1   1    1 ]);     % penalize angular velocities         
% R needs to be positive definite
Cost.R_lqr = 10*diag([1 1 1 1]);        % penalize control inputs

%% ILQC cost function
Cost.Qm  = Cost.Q_lqr;                  % it makes sense to redefine these
Cost.Rm  = Cost.R_lqr;
Cost.Qmf = Cost.Q_lqr;

% Reference states and input the controller is supposed to track
gravity = 9.81;
f_hover = m_quad*gravity; % keep input close to the one necessary for hovering
Cost.u_eq = [ f_hover ; zeros(3,1) ];
Cost.x_eq = Task.goal_x;

% alias for easier referencing
x = x_sym; 
u = u_sym;
x_goal = Task.goal_x;

%ilqc_type = 'goal_state'; % Choose 'goal_state or 'via_point'
ilqc_type = 'via_point'; % Choose 'goal_state or 'via_point'
fprintf('ILQC cost function type: %s \n', ilqc_type);
switch ilqc_type
    case 'goal_state'     
        Cost.h = simplify((x-x_goal)'*Cost.Qmf*(x-x_goal));
        Cost.l = simplify( (x-Cost.x_eq)'*Cost.Qm*(x-Cost.x_eq) ...
             + (u-Cost.u_eq)'*Cost.Rm*(u-Cost.u_eq));
               
    case 'via_point'
%% [Problem 1.2 (e)] Include via_point p1 in the ILQC cost function formulation
        %p1 = Task.vp1;
        p1 = Task.vp2;% p1 = Task.vp2 also try this one
        t1 = Task.vp_time;

        % =================================================================
        % [Todo] Define an appropriate weighting for way points (see
        % handout Eqn.(19)) Hint: Which weightings must be zero for the
        % algorithm to determine optimal values?
        %
        % Q_vp = ...;
        % =================================================================

        % don't penalize position deviations, drive system with final cost
        Cost.Qm(1:3,1:3) = zeros(3);
        Q_vp = Cost.Qm;
        Q_vp(1:6,1:6) = zeros(6);
        Q_vp = eye(12);
        % =================================================================
        % [Todo] Define symbolic cost function. 
        % Note: Use function "viapoint(.)" below
        %
        % viapoint_cost = viapoint(t1,p1,x,t_sym,Q_vp);
        % Cost.h = ...;
        % Cost.l = ... + viapoint_cost;
        % =================================================================       
        viapoint_cost = viapoint(t1,p1,x,t_sym,Q_vp);
        Cost.h = simplify((x-x_goal)'*Cost.Qmf*(x-x_goal));
        Cost.l = simplify( (x-Cost.x_eq)'*Cost.Qm*(x-Cost.x_eq) ...
             + (u-Cost.u_eq)'*Cost.Rm*(u-Cost.u_eq)) + viapoint_cost;
        
      otherwise
        error('Unknown ilqc cost function mode');
end 

Cost.x  = x;
Cost.u  = u;
Cost.t  = t_sym; 
end


function viapoint_cost = viapoint(vp_t,vp,x,t,Qm_vp)
% WAYPOINT Creates cost depending on deviation of the system state from a
% desired state vp at time t. Doesn't need to be modified. 
% For more details see:
% http://www.adrlab.org/archive/p_14_mlpc_quadrotor_cameraready.pdf
%    vp_t  :   time at which quadrotor should be at viapoint wp
%    vp    :   position where quad should be at given time instance
%    x,t   :   symbolic references to state and time
%    Qm_vp :   The weighting matrix for deviation from the via-point

prec = 3; % how 'punctual' does the quad have to be? The bigger the
          % number, the harder the time constraint
          
viapoint_cost = (x-vp)'*Qm_vp*(x-vp) ...
                *exp(-0.5*prec*(t-vp_t)^2) /sqrt(2*pi/prec);

end

