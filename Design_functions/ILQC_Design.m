% ILQC_Design: Implementation of the ILQC controller.
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
% Modified and completed by Francisco Granda

function [Controller,cost] = ILQC_Design(Model,Task,Controller,Simulator)
% ILQC_DESIGN Implements the Iterative Linear Quadratic Controller (ILQC)
%    (see Ch. 4 notes for a formal description of the algorithm)

% Define functions that return the quadratic approximations of the cost 
% function at specific states and inputs (Eqns. (6)-(7) in handout)
% Example usage:
%     xn = [ x y z ... ]'; 
%     un = [ Fx Mx My Mz]';
%     t  = t;
%     Qm(xn,un) = Qm_fun(t,xn,un); 

% stage cost (l) quadratizations
l_ = Task.cost.l*Task.dt;
q_fun   = matlabFunction(  l_,'vars',{Task.cost.t,Task.cost.x,Task.cost.u});
% dl/dx
l_x = jacobian(Task.cost.l,Task.cost.x)'*Task.dt; % cont -> discr. time
Qv_fun  = matlabFunction( l_x,'vars',{Task.cost.t,Task.cost.x,Task.cost.u});
% ddl/dxdx
l_xx = jacobian(l_x,Task.cost.x);
Qm_fun  = matlabFunction(l_xx,'vars',{Task.cost.t,Task.cost.x,Task.cost.u});
% dl/du
l_u = jacobian(Task.cost.l,Task.cost.u)'*Task.dt; % cont -> discr. time
Rv_fun  = matlabFunction( l_u,'vars',{Task.cost.t,Task.cost.x,Task.cost.u});
% ddl/dudu
l_uu = jacobian(l_u,Task.cost.u);
Rm_fun  = matlabFunction(l_uu,'vars',{Task.cost.t,Task.cost.x,Task.cost.u});
% ddl/dudx
l_xu = jacobian(l_x,Task.cost.u)';
Pm_fun  = matlabFunction(l_xu,'vars',{Task.cost.t,Task.cost.x,Task.cost.u});

% terminal cost (h) quadratizations
h_ = Task.cost.h;
qf_fun  = matlabFunction(  h_,'vars',{Task.cost.x});
% dh/dx
h_x = jacobian(Task.cost.h,Task.cost.x)';
Qvf_fun = matlabFunction( h_x,'vars',{Task.cost.x});
% ddh/dxdx
h_xx = jacobian(h_x,Task.cost.x);
Qmf_fun = matlabFunction(h_xx,'vars',{Task.cost.x});

% dimensions
n = length(Task.cost.x); % dimension of state space
m = length(Task.cost.u); % dimension of control input
N  = (Task.goal_time-Task.start_time)/Task.dt + 1; % number of time steps

% desired value function V* is of the form Eqn.(9) in handout
% V*(dx,n) = s + dx'*Sv + 1/2*dx'*Sm*dx
s    = zeros(1,N);
Sv   = zeros(n,N);
Sm   = zeros(n,n,N);

% Initializations
theta_temp = init_theta();
duff = zeros(m,1,N-1);
K    = repmat(theta_temp(2:end,:)', 1, 1, N-1);
sim_out.t = zeros(1, N);
sim_out.x = zeros(n, N);
sim_out.u = zeros(m, N-1);
X0 = zeros(n, N);
U0 = zeros(m, N-1);

% Shortcuts for function pointers to linearize systems dynamics:
% e.g. Model_Alin(x,u,Model_Param)
Model_Param = Model.param.syspar_vec;
Model_Alin  = Model.Alin{1}; 
Model_Blin  = Model.Blin{1}; 

%% Implementation of ILQC controller
% Each ILQC iteration approximates the cost function as quadratic around the
% current states and inputs and solves the problem using DP.
% Linearization around equilibrium and goal state
x_lin1 = Task.cost.x_eq;
u_lin1 = Task.cost.u_eq;
A_lin1 = Model_Alin(x_lin1,u_lin1,Model_Param);
B_lin1 = Model_Blin(x_lin1,u_lin1,Model_Param);

% Discretization of system
dt = Task.dt;
A_k1 = (eye(size(A_lin1)) + A_lin1*dt);
B_k1 = B_lin1*dt;
i  = 1;
%while ( i <= Task.max_iteration && ( norm(squeeze(duff)) > 0.01 || i == 1 ))
while ( i <= Task.max_iteration)    
    %% Forward pass / "rollout" of the current policy
    % Forward Pass
    % simulation
    sim_out = Simulator(Model,Task,Controller);
    % pause if cost diverges
    cost(i) = Calculate_Cost(sim_out, q_fun, qf_fun);
    fprintf('Cost of Iteration %2d (metric: ILQC cost function!): %6.4f \n', i-1, cost(i));
    
    if ( i > 1 && cost(i) > 2*cost(i-1) )
        fprintf('It looks like the solution may be unstable. \n')
        fprintf('Press ctrl+c to interrupt iLQG, or any other key to continue. \n')
        pause
    end
    
    %% Solve Riccati-like equations backwards in time
	% =====================================================================
    % result of forward pass
    X0 = sim_out.x;
    U0 = sim_out.u;
    T0 = sim_out.t;
    % Quadratization Final Cost
    xf = X0(:,end);
    s_N(N) = qf_fun(xf);
    s_NN(:,N) = Qvf_fun(xf);
    S_N(:,:,N) = Qmf_fun(xf);

    for k = N-1:-1:1
        % Actual State
        x0 = X0(:,k);
        u0 = U0(:,k);
        t0 = T0(:,k);
        % Linearization and discretization
        A_lin = Model_Alin(x0,u0,Model_Param);
        B_lin = Model_Blin(x0,u0,Model_Param);
        A_k = (eye(size(A_lin)) + A_lin*dt);
        B_k = B_lin*dt;
        % Stage cost quadratization
        q_k = q_fun(t0,x0,u0);
        q_kk = Qv_fun(t0,x0,u0);
        Q_k = Qm_fun(t0,x0,u0);
        r_k = Rv_fun(t0,x0,u0);
        R_k = Rm_fun(t0,x0,u0);
        P_k = Pm_fun(t0,x0,u0);
        % Control dependant terms
        s_N_in = s_N(k+1);
        s_NN_in = s_NN(:,k+1);
        S_N_in = S_N(:,:,k+1);
        g_k = r_k + B_k.' * s_NN_in;
        G_k = P_k + B_k.' * S_N_in*A_k;
        H_k = R_k + B_k.' * S_N_in*B_k;
        H_k = (H_k+H_k')/2;
        % S compute
        du_ff = -inv(H_k)*g_k;
        K_k = -inv(H_k)*G_k;
        s_N(k) = q_k + s_N_in + 1/2*du_ff.'*H_k*du_ff + du_ff.'*g_k;
        s_NN(:,k) = q_kk + A_k.'*s_NN_in+K_k.'*H_k*du_ff+...
            K_k.'*g_k+G_k.'*du_ff;
        S_N(:,:,k) = Q_k + A_k.'*S_N_in*A_k + K_k.'*H_k*K_k + ...
            K_k.'*G_k + G_k.'*K_k;
        theta_ff = du_ff + u0 -(K_k*x0);
        theta_k = [theta_ff,K_k].';
        Controller.theta(:,:,k) = theta_k;
    end
    i = i+1;
end

% simulating for the last update just to calculate the final cost
sim_out    = Simulator(Model,Task,Controller);
cost(i) = Calculate_Cost(sim_out, q_fun, qf_fun);
fprintf('Cost of Iteration %2d: %6.4f \n', i-1, cost(i));
end



function theta = Update_Controller(X0,U0,dUff,K)
% UPDATE_CONTROLLER Updates the controller after every ILQC iteration
%
%  X0  - state trajectory generated with controller from previous 
%        ILQC iteration.
%  UO  - control input generated from previous ILQC iteration.
%  dUff- optimal update of feedforward input found in current iteration
%  K   - optimal state feedback gain found in current iteration
%
%  The updated control policy has the following form:
%  U1 = U0 + dUff + K(X - X0)
%     = U0 + dUff - K*X0 + K*X
%     =      Uff         + K*x
%  
%  This must be brought into the form 
%  U1 = theta' * [1,x']   

%% Update Controller
% input and state dimensions
n  = size(X0,1); % dimension of state
m = size(U0,1); % dimension of control input
N  = size(X0,2); % number of time steps

% initialization
theta = init_theta();
theta_fb = zeros(n, m, N-1);
theta_ff = repmat(theta(1,:), 1, 1, N-1);

% =========================================================================
% feedforward control input

% feedback gain of control input (size: n * m * N-1)
theta_fb = permute(K,[2 1 3]);      

% puts below (adds matrices along first(=row) dimension). 
% (size: (n+1) * m * N-1)
theta = [theta_ff; theta_fb];  

end


function cost = Calculate_Cost(sim_out, q_fun, qf_fun)
% CALCULATE_COST: calcules the cost of current state and input trajectory
% for the current ILQC cost function. Not neccessarily the same as the LQR
% cost function.

X0 = sim_out.x(:,1:end-1);
xf = sim_out.x(:,end);
U0 = sim_out.u;
T0 = sim_out.t(1:end-1);

cost = sum(q_fun(T0,X0,U0)) + qf_fun(xf);
end