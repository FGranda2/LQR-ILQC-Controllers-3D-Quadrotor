% Compute_Input: Compute control inputs.
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
% [14.05.25]    first version

function u = Compute_Input(t,x,Controller,param)

% Controller: Controller structure
%             u(x) = theta([t]) ' * BaseFnc(x)
%             .theta:   Parameter matrix (each column is associated to one
%                       control input)
%             .BaseFnc
%             .time
%
% tau_sat = [tau_min(:) tau_max(:)]
% NOTE: t here is the time instances for which we require the input and not
% the simulation time. Hence the input length is length(t), not length(t)-1

if isempty(Controller.time)     % Steady-State Controller
    theta = repmat(Controller.theta(:,:,1),[1 1 length(t)]);
else                            % Controller changes over time
    % select correct controller: note that ODE sometimes jumps back in time
    % both time vectors need to be row vectors
    index = sum(repmat(Controller.time,length(t),1) <= repmat(t',1,length(Controller.time)),2);
    index = min( index , size(Controller.theta,3) );
    theta = Controller.theta(:,:,index);
end

% Adding noise to input parameters
rand_var = Controller.rand_var;
index = sum(repmat(rand_var.t,length(t),1) <= repmat(t',1,length(rand_var.t)),2);
rand_theta = rand_var.v_theta(:,:,index);
theta = theta + rand_theta;

% Computing control inputs (n_inputs = size(theta,2))
% u = "tau" = (Fz, Mx, My, Mz)' = theta' * x_aug = uff + K'*x
temp = theta .* permute( repmat(Controller.BaseFnc(t,x),[1 1 size(theta,2)]) ,[1 3 2] );
u = permute(sum(temp,1),[2 3 1]);

% Compute Control Input
% Thrust (=Force) of every rotor individually
thrust = Compute_Thrust(u,param);

% Adding noise to individual thrusts
index = sum(repmat(rand_var.t,length(t),1) <= repmat(t',1,length(rand_var.t)),2);
rand_thrust = rand_var.v_input(:,index);
thrust = thrust + rand_thrust;

% Control Saturation
thrust(thrust<param.Fsat(1)) = param.Fsat(1); % lower bound
thrust(thrust>param.Fsat(2)) = param.Fsat(2); % upper bound

% transform thrust back to (Fz, Mx, My, Mz), since that is what the model
% works with
a1 = param.kM/param.kF;
u = [ 1 1 1 1 ; 
    0 param.La 0 -param.La ; 
    -param.La 0 param.La 0 ;
    a1 -a1 a1 -a1 ]*thrust;

end