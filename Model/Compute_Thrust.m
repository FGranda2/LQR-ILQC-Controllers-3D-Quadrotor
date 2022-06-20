% Compute_Thrust: Computes rotor thrust from force and moments
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

function Ft = Compute_Thrust(tav,param)
% COMPUTE_THRUST: returns thrust of each rotor when given the force and
% moments generated, when assuming Ft(i) = kF*w_i^2, where w_i is the
% rotational speed of rotor i
%
% tav = [ Fz ; Mx; My; Mz ] in body-fixed frame
% param: struct with
%           kF: thrust coefficient, Ft(i) = kF*wi^2
%           kM: moment coefficient, i.e. Mz = kM*(w1^2-w2^2+w3^2-w4^2);
%           La: length of quadcopter arm
%

a1 = 1/(2*param.La);
a2 = param.kF/(4*param.kM);

Ft = [ 0.25 0 -a1 a2 ; 
    0.25 a1 0 -a2 ; 
    0.25 0 a1 a2 ; 
    0.25 -a1 0 -a2 ]*tav;
end