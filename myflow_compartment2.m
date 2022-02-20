function [cv_compartments2] = myflow_compartment2(t)
% t;

global flow_compartments2

% to apply noise to flow rates. optional. default: no noise

cv_compartments2 = flow_compartments2+0*randn(size(t));
end

