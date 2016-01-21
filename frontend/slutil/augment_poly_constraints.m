function [C_total,d_total] = augment_poly_constraints(scsbHandle,pthb)

% Augment constraints to the overall system order
%
% Syntax:
%   "[C_total,d_total] = augment_poly_constraints(scsbHandle,pthb)"
%
% Description:
%   "augment_poly_constraints(scsbHandle,pthb)" returns the constraints
%   from the polyhedral threshold block (PTHB) "pthb" with respect to the
%   overall system order.  "scsbHandle" is a list of handles to all
%   switched continuous system blocks (SCSBs) in the system, and is
%   used to determine the overall order of the system.
%
% See Also:
%   piha

% Get the (C,d) matrix-vector pair for the specified PTHB
poly_pthb = evalin('base',get_param(pthb,'polyhedron'));
[CE,dE,C_pthb,d_pthb] = linearcon_data(poly_pthb);
if ~isempty(CE) | ~isempty(dE)
  blockname = get_param(pthb,'name');
  error(['Invalid constraints for polyhedral threshold block ''' ...
	blockname '''.'])
end

% Find the SCSBs feeding to the input of the specified PTHB
input_scsb = trace_pthb_input(pthb);

% Extract the columns of C_pthb corresponding to each input SCSB
C_input_scsb = {};
nx_total = 0;
for k = 1:length(input_scsb)
  nx = eval(get_param(input_scsb(k),'nx'));
  if strcmp(get_param(input_scsb(k),'use_sd'),'on')
      nx = nx+eval(get_param(input_scsb(k),'nz'));
  end
  C_input_scsb{k} = C_pthb(:,nx_total+1:nx_total+nx);
  nx_total = nx_total + nx;
end

% Rearrange the columns according to the overall SCSB in the system.
C_total = [];
for k = 1:length(scsbHandle)
  % Check if the kth SCSB is part of the input to pthb
  found = 0;
  for l = 1:length(input_scsb)
    if (scsbHandle(k) == input_scsb(l))
      found = 1;
      break;
    end
  end
  if found
    % if found, append the columns corresponding to the kth SCSB
    C_total = [C_total C_input_scsb{l}];
  else
    % if not, append the zero columns
    nx = eval(get_param(scsbHandle(k),'nx'));
    C_total = [C_total zeros(size(C_pthb,1),nx)];
  end
end

d_total = d_pthb;
return
