function [cv_compartments1] = myflow_compartment1(t)
% t;

global flow_compartments1

% to apply noise to flow rates (optional, default: no noise)

cv_compartments1 = flow_compartments1+0*randn(size(t));
end

