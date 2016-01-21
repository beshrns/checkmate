function st = check_model_syntax(sys)

% Check syntax of a CheckMate model
%
% Syntax:
%   "st = check_model_syntax(sys)"
%
% Description:
%   "check_model_syntax(sys)" returns "1" if "sys" conforms to proper
%   CheckMate syntax, otherwise, "check_model_syntax(sys)" returns "0".
%
% Note:
%   Currently, proper CheckMate syntax requires that the names of
%   polyhedral threshold and finite state machine blocks (PTHBs) &
%   (FSMBs) begin with a lower case letter and contain only letters,
%   numbers, or underscores. 

% Set default return status to 1
st = 1;

% ----------------------
% Check PTHB block names
% ----------------------

pthb_list = ...
    find_system(sys,'SearchDepth',2,'MaskType','PolyhedralThreshold');
for k = 1:length(pthb_list)
  name = get_param(pthb_list{k},'Name');
  if ~valid_block_name(name)
    fprintf(1,'\007Invalid polyhedral threshold block name ''%s''.\n',name)
    fprintf(1,['PTHB name must be a lower case letter followed by' ...
	  ' letters, numbers, or _''s.'])
    st = 0;
    return
  end
end

% ----------------------
% Check FSMB block names
% ----------------------

fsmb_list = ...
    find_system(sys,'SearchDepth',2,'MaskType','Stateflow');
for k = 1:length(fsmb_list)
  name = get_param(fsmb_list{k},'Name');
  if ~valid_block_name(name)
    fprintf(1,'\007Invalid finite state machine block name ''%s''.\n',name)
    fprintf(1,['FSMB name must be a lower case letter followed by' ...
	  ' letters, numbers, or _''s.'])
    st = 0;
    return
  end
end




if 0 % begin DISABLED code
  
  % -----------------------------
  % Check Composition Order Block
  % -----------------------------

  compose_block = check_param_block(sys,'CompositionOrder');
  if isempty(compose_block)
    st = 0;
    return
  end

  % Check validity of SCS Block order
  scsb_order = str2cell(get_param(compose_block,'scsb_order'));
  scsbH = ...
      find_system(sys,'SearchDepth',2,'MaskType','SwitchedContinuousSystem');
  if length(scsb_order) ~= length(scsbH)
    fprintf(1,'\007Incorrect number of SCS blocks specified in ''%s''.\n', ...
	compose_block)
    st = 0;
    return
  end
  for k = 1:length(scsb_order)
    if ~ismember([sys '/' scsb_order{k}],scsbH)
      fprintf(1,'\007Invalid SCS Block ''%s'' specified in ''%s''.\n', ...
	  [sys '/' scsb_order{k}],compose_block)
      st = 0;
      return
    end
  end

  % Check validity of FSM Block order
  fsmb_order = str2cell(get_param(compose_block,'fsmb_order'));
  fsmbH = find_system(sys,'SearchDepth',2,'MaskType','Stateflow');
  if length(fsmb_order) ~= length(fsmbH)
    fprintf(1,'\007Incorrect number of FSM blocks specified in ''%s''.\n', ...
	compose_block)
    st = 0;
    return
  end
  for k = 1:length(fsmb_order)
    if ~ismember([sys '/' fsmb_order{k}],fsmbH)
      fprintf(1,'\007Invalid FSM Block ''%s'' specified in ''%s''.\n', ...
	  [sys '/' fsmb_order{k}],compose_block)
      st = 0;
      return
    end
  end

  % ---------------------------
  % Check Analysis Region Block
  % ---------------------------

  AR_block = check_param_block(sys,'AnalysisRegion');
  if isempty(AR_block)
    st = 0;
    return
  end

  % ----------------------------------
  % Check Initial Continuous Set Block
  % ----------------------------------

  ICS_block = check_param_block(sys,'InitialContinuousSet');
  if isempty(ICS_block)
    st = 0;
    return
  end

end % end DISABLED code




% ---------------------------
% Check Start Event Generator
% ---------------------------

%The following has been commented out because new version of Checkmate does not have 'start event' generators
%******************************************************************************************
%start_block = find_system(sys,'SearchDepth',1,'BlockType','Step', ...
%    'Name','start');
%if isempty(start_block)
%  fprintf(1,'\007Start event generator is missing.\n')
%  st = 0;
%  return
%end
%start_block = start_block{1};

% check that start event generator is hooked up to each FSM properly
% (to be done later)
%*******************************************************************************************

return



% -----------------------------------------------------------------------------

function st = valid_block_name(name)

st = 1;
for l = 1:length(name)
  ch = name(l);
  if (l == 1) 
    % first character
    if ~islowercase(ch)
      st = 0;
      return
    end
  else
    % second character on
    if ~((ch == '_') | isletter(ch) | isnumber(ch))
      st = 0;
      return
    end
  end
end
return

% -----------------------------------------------------------------------------

function st = islowercase(ch)

st = ('a' <= ch) & (ch <= 'z');
return

% -----------------------------------------------------------------------------

function st = isnumber(ch)

st = ('0' <= ch) & (ch <= '9');
return

% -----------------------------------------------------------------------------

function param_block = check_param_block(sys,param_type)

param_block = find_system(sys,'SearchDepth',2,'MaskType',param_type);

% Must have exactly one parameter block for each type.
if isempty(param_block)
  fprintf(1,'\007%s block missing.\n',param_type)
  param_block = '';
  return
end

if length(param_block) > 1
  fprintf(1,'\007Found more than one %s block. ',param_type)
  fprintf(1,'Please delete all excess blocks.\n')
  param_block = '';
  return
end

param_block = param_block{1};
return
