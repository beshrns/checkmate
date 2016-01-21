function [destination,null_event,time_limit,out_of_bound,terminal] = ...
    compute_mapping_no_SD(X0,srcloc,srccell)

% Perform reachability analysis to compute the destination locations and
% cells from the given polytope in the specified location and cell.
%
% Syntax:
%   "[destination,null_event,time_limit,out_of_bound,terminal] =
%    compute_mapping_no_SD(X0,srcloc,srccell)"
%
% Description:
%   The inputs to this function are
%
%   * "X0": a "linearcon" object representing the initial continuous set
%     from which the reachability is to be computed. It is assumed that "X0"
%     is full-dimensional and resides in the given source cell which in
%     turns is part of the interior region for the source location.
%   * "srcloc": current location.
%   * "srccell": current polyhedral cell where "X0" resides
%
%   The outputs of this function are
%
%   * "destination": a cell array of structures containing the
%     information about next destination locations and cells that are
%     reachable from the given initial condition (see note below for more
%     detail on the structure).
%   * "null_event": a boolean flag indicating that the flow pipe computation
%     was terminated because it can be concluded that the subsequent flow
%     pipe segments will remain inside "INV" forever.
%   * "time_limit": a boolean flag indicating that the flow pipe computation
%     was terminated because the time limit "max_time" was exceeded.
%   * "out_of_bound": a boolean flag indicating that the flow pipe computation
%     resulted in some mapppings polytopes that went out of the analysis
%     region.
%   * "terminal": a cell array of FSMB state (rows) vectors indicating
%     the terminal FSM states that can be reached from the given initial
%     condition.
%
% Note:
%   The output argument "destination" is a cell array of structures.  Each
%   element in the cell array correspond to a destination with the following
%   fields.
%
%   * ".location" - destination location index
%   * ".cell" - destination cell index
%   * ".mapping" - a single polytope representing the reachable set from
%                  the given initial condition in the destination location
%                  and cell.
%   * ".Tstamp" - not used here since no sampled data is involved.
%
% Implementation:
%   Call the mapping computation routine corresponding to the type of the
%   composite continuous dynamics of all SCSBs for the given location
%   ("fs_lin_map" for `linear` (affine) dynamics, "fs_nonlin_map" for
%   `nonlinear` dynamics, and "clk_map" for `clock` dynamics) to map the
%   continuous trajectory from "X0" to the boundaries of the current
%   cell. For each cell boundary that has been reached determine the
%   following.
%
%   1. If the mapping is entering another interior cell, then add a
%   destination corresponding to that interior cell. The destination
%   location and mapping remain unchanged.
%
%   2. If the mapping is entering a guard cell and at least one transition
%   in one FSM is enabled, compute all possible destination locations. For
%   each possible destination location, apply the reset transformation to
%   the mapping as needed and find the interior cells in the destination
%   location that can be reached by the transformed mapping. Add all
%   combinations of destination interior cells and locations as well as the
%   correponding mapping to the destination cell array.
%
%   3. If the mapping does not enable any transition, find all guard cells
%   that can be reached from the cell partitions of all transitions. Then
%   eliminate all the guard cells that are already covered by some other
%   guard cells. Add the remaining cells to the destination list. As with
%   case 1, the destination location and the mapping remain unchanged.
%
%   To speed up the subsequent computation, the mapping polytopes are over
%   approximate by a single polytope before the result is stored in the
%   ".mapping" field above.
%
% See Also:
%   fs_nonlin_map,fs_lin_map,clk_map

global GLOBAL_PIHA
global GLOBAL_APPROX_PARAM


% Default return values.
null_event = 0;
time_limit = 0;
out_of_bound = 0;
destination = {};
terminal = {};

% Do nothing if the initial condition is empty.
if isempty(X0)
  return
end

% Get the polytope for the current cell.
INV = return_invariant(srccell);

% Get the constraints on the parameters for the given location.
Pcon = return_parameter_cons();

% Compute the type of the composite dynamics of all SCSBs for the FSM
% state in the source location.
overall_dynamics = check_overall_dynamics(GLOBAL_PIHA.SCSBlocks, ...
                                          GLOBAL_PIHA.Locations{srcloc}.q);

% Apply different methods for computing the mapping polytopes for each
% type of dynamics.
switch overall_dynamics,
 
 case 'linear',
  [A,b] = overall_system_matrix(GLOBAL_PIHA.SCSBlocks, ...
                                GLOBAL_PIHA.Locations{srcloc}.q);
  [mapping,null_event,time_limit] = fs_lin_map(A,b,X0,INV);
 
 case 'clock',
  v = overall_system_clock(GLOBAL_PIHA.SCSBlocks, ...
                           GLOBAL_PIHA.Locations{srcloc}.q);
  % Check to make sure that the clock vector is not zero before we
  % proceed to compute the mapping by quantifier elimination.
  if all(v == 0)
    error('CheckMate:ComputeMapping:ZeroClock', ...
        ['Found zero clock vector for location ' num2str(srcloc) '.']);
    return
  end
  % Compute the mapping by projecting X0 along the clock direction and
  % intersect it with the current cell polytope.
  mapping = clk_map(X0,v,INV);
  % As long as the clock vector is not zero, the system cannot remain in
  % the same location forever, hence the null event flag is always 0 in
  % this case.
  null_event = 0;
  % Since the above mapping computation by quantifier elimination always
  % terminates without having to check for the time-limit, the time-limit
  % flag is always 0 in this case.
  time_limit = 0;
 
 case 'nonlinear',
  % Why is q an ode_param here? Need to confirm this with Ansgar.
  ode_param = GLOBAL_PIHA.Locations{srcloc}.q;
  [mapping,null_event,time_limit] = ...
      fs_nonlin_map('overall_system_ode',ode_param,X0,INV,Pcon, ...
                    GLOBAL_APPROX_PARAM.T,GLOBAL_APPROX_PARAM.max_time);
 
 otherwise,
  error('CheckMate:ComputeMapping:WrongDynamics', ...
      ['Unknown dynamics type ''' overall_dynamics '''.'])

end

% For each cell face with non-empty mapping, find where the mapping leads
% (i.e. what location and cell) as described in the comments at the start of
% the file.
destination = [];
terminal = [];
NAR = GLOBAL_PIHA.NAR;
transitions = GLOBAL_PIHA.Locations{srcloc}.transitions;
interior_cells = GLOBAL_PIHA.Locations{srcloc}.interior_cells;
for k = 1:length(mapping)
  % Do nothing and skip to the next cell face if the mapping for the current
  % cell face is empty.
  if isempty(mapping{k})
    continue;
  end
  
  % Get the global hyperplane index for the current cell face.
  cell_face_hp_idx = GLOBAL_PIHA.Cells{srccell}.boundary(k);
  
  % If the current cell face is a boundary of the analysis region, set the
  % out-of-bound flag to 1.
  if cell_face_hp_idx <= NAR
    out_of_bound = 1;
    % Since there can't possibly be another neighboring cell when the
    % current face is a boundary of the analysis region, there is no
    % destination cell. Therefore, we can skip the rest of the loop and
    % proceed to the next cell face.
    continue;
  end
  
  % Find destinations within the interior cells for the current location for
  % the mapping on the current face of the current cell.
  interior_dst = find_interior_destination( ...
      cell_face_hp_idx,interior_cells,srcloc,srccell,mapping{k});
  destination = [destination interior_dst];

  % Find destinations outside the interior cells.
  [guard_dst,terminal_k] = find_guard_destination( ...
      cell_face_hp_idx,transitions,srcloc,srccell,mapping{k});
  destination = [destination guard_dst];
  terminal = [terminal; terminal_k];
end

% Convert the return variable "terminal" to a cell array of row vectors as
% required by the caller function.
temp = cell(size(terminal,1), 1);
for k = 1:size(terminal,1)
  temp{k} = terminal(k,:);
end
terminal = temp;

% Convert the return variable "destination" to a cell array of structures,
% which is inefficient, but required by the caller function. Need to remove
% all cell arrays of structures from the code and replace them with
% structure array.
temp = cell(size(destination,1), 1);
for k = 1:length(destination)
  temp{k} = destination(k);
end
destination = temp;

return

% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

function destination = find_interior_destination( ...
    cell_face_hp_idx,interior_cells,srcloc,srccell,mapping)

% Find destinations within the interior region for the current location for
% the mapping on the current face of the current cell.

global GLOBAL_PIHA

% Over-approximate the mapping by a single polytope.
box = bounding_box(mapping);

destination = [];
% Find the neigboring cells (other than the current cell) in the interior
% region for the current location that have been hit by the mapping
% polytopes.
for i = 1:length(interior_cells)
  if interior_cells(i) ~= srccell
    % Check if the current cell face for the source cell is also a boundary
    % face of the current interior cell.
    [dum1,boundary_idx] = ...
        intersect(GLOBAL_PIHA.Cells{interior_cells(i)}.boundary, ...
                  cell_face_hp_idx);
    if ~isempty(boundary_idx)
      % If so, check further if any mapping polytope actually lands on the
      % neighboring cell face.
      dst_cell_inv = return_invariant(interior_cells(i));
      for j = 1:length(mapping)
        % If a mapping polytope overlaps with the invariant, then we can
        % conclude that the current interior cell can be reached. Add the
        % current interior cell to the destination list. 
        if isfeasible(box, dst_cell_inv)
          % Since reaching another interior cell cannot trigger a transition
          % the location must remain the same. To reduce the computational
          % complexity later on, we also over approximate the mapping
          % polytopes by a single polytope the covers all the mapping
          % polytopes. Since we are interested in finding the mapping
          % polytopes that reach each destination cell, the resulting
          % approximation of the mapping is also intersected with the
          % destination cell.
          temp.cell = interior_cells(i);
          temp.location = srcloc; 
          temp.mapping = box & dst_cell_inv;
          temp.Tstamp = [];
          destination = [destination temp];
          break;
        end % if ~isempty(intersection)
      end % for j
    end % if ~isempty(boundary_idx)
  end % if interior(i) ~= srccell
end % for i

return
 
% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

function [destination,terminal] = ...
    find_guard_destination(cell_face_hp_idx,transitions,srcloc,srccell,mapping)

% Find destination cells in cases where a location transition has been
% taken.

global GLOBAL_PIHA

% Find the list of transitions being enabled by the mapping polytopes for
% the current cell face. Keep a separate list for each FSMB. Also, find the
% list of (possibly overlapping) guard cells in the cell partition for each
% transition that have been reached by the mapping polytopes.
q = GLOBAL_PIHA.Locations{srcloc}.q;
enabled = cell(1,length(q));
guard_cells_reached = [];
enabled_transition_found = 0;
for j = 1:length(transitions)
  % Get all guard cells associated with the current transition.
  guard_cells = transitions{j}.guard;
  guard_cell_event_flags = transitions{j}.guard_cell_event_flags;
  for i = 1:length(guard_cells)
    % Check if the current cell face is also a boundary face of the current
    % guard cell.
    [dum1,boundary_idx] = ...
        intersect(GLOBAL_PIHA.Cells{guard_cells(i)}.boundary, ...
                  cell_face_hp_idx);
    if ~isempty(boundary_idx)
      % If so, check further if any mapping polytope actually lands on the
      % guard cell face.
      dst_cell_inv = return_invariant(guard_cells(i));
      for k = 1:length(mapping)
        if isfeasible(mapping{k}, dst_cell_inv) 
          % If the event flag set to 1 for the current guard cell face being
          % hit, then we know that the parent transition for this guard cell
          % must be enabled. Add the enabled transition to the list for
          % its parent FSMB.
          if guard_cell_event_flags{i}(boundary_idx)
            fsm_idx = transitions{j}.idx;
            enabled{fsm_idx} = [enabled{fsm_idx} j];
            enabled_transition_found = 1;
          end
          % Keep the list of guard cells reached by the mapping.
          guard_cells_reached = union(guard_cells_reached,guard_cells(i));
          break;
        end
      end
    end
  end
end

destination = [];
terminal=[];

if enabled_transition_found
  
  % If there are enabled transitions, find all next possible
  % locations. 
  [dest_loc,dest_scs_reset,terminal] = ...
      find_destination_location(srcloc,transitions,enabled);

  % For each non-terminal destination location, find the destination
  % interior cells and a new entry into the destination list for each
  % destination cell found.
  for i = 1:length(dest_loc)
    % Apply the reset transformation for the destination location to
    % the mapping polytopes.
    q = GLOBAL_PIHA.Locations{dest_loc(i)}.q;
    scs_reset_indices = dest_scs_reset{i};
    [T,v] = overall_system_reset(GLOBAL_PIHA.SCSBlocks,q,scs_reset_indices);
    reset_mapping = cell(length(mapping),1);
    for k = 1:length(mapping)
      % The new transform routine should give full dimensional polytope even
      % if T is singular.
      reset_mapping{k} = transform(mapping{k},T,v);
    end
    
    % Intersect each interior cell with the reset mapping polytopes. Add
    % each non-empty result to the destination list.
    found = 0;
    interior_cells = GLOBAL_PIHA.Locations{dest_loc(i)}.interior_cells;
    for j = 1:length(interior_cells)
      dst_cell_inv = return_invariant(interior_cells(j));
      intersection = {};
      for k = 1:length(reset_mapping)
        if isfeasible(reset_mapping{k}, dst_cell_inv)
          intersection{end+1} = reset_mapping{k};
        end
      end
      if ~isempty(intersection)
        % Over approximate the intersection between the reset mapping
        % polytopes and the destination cell by a single polytope and
        % intersect it with the destination cell again before adding
        % the mapping to the destination list.
        temp.cell = interior_cells(j);
        temp.location = dest_loc(i); 
        temp.mapping = bounding_box(intersection) & dst_cell_inv;
        temp.Tstamp = [];
        destination = [destination temp];
        % Turn on the flag indicating that we have found at least one
        % destination cell in the interior region.
        found = 1;
      end
    end
      
    % If no destination interior cell is found, issue an error message to
    % the user about this.
    if ~found
      src_name = location_name(srcloc);
      dst_name = location_name(dest_loc(i));
      msg = [sprintf('\n\n') 'Transition from location ' src_name ...
             ' to location ' dst_name ' lands the system completely ' ...
             'outside of ' sprintf('\n') 'the interior region of the ' ...
             'destination location. Reachability not reliable!'];
      error('CheckMate:ComputeMapping:TransitionNotFound', msg);
    end
  end %for i

else

  % If no transition is enabled but some guard cells have been reached,
  % warn the user about this.
  if ~isempty(guard_cells_reached)
    msg = ['Warning: Reachability analysis indicates that the system may ' ...
           'exit the interior region' sprintf('\n') 'without taking any ' ...
           'transition from cell ' num2str(srccell) ' of location ' ...
           location_name(srcloc) '.'];
    fprintf(1,'%s\n',msg);
    return
  end
  
end
  
return

% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

function [dest_loc,dest_reset,terminal] = ...
    find_destination_location(srcloc,transitions,enabled)
                 
% Given the list of enabled transitions for all FSMBs, find the list of
% next possible destination locations and the indices of SCSBs to be
% reset for each destination location.
%
% Inputs:
%   * "srcloc" - Source location index.
%   * "transitions" - A structure representing the list of FSM transitions
%      for the current location.
%   * "enabled" - A cell array of the same length as the number of
%     FSMBs. Each cell "enabled{i}" is a vector containing the indices for
%     the transitions that have been enabled for the i-th FSMB.
%
% Outputs:
%   * "dest_loc" - A vector of location indices representing the list of
%     next possible locations.
%   * "dest_reset" - A cell array of vectors. Each element "dest_reset{i}" 
%     of the cell array contains the indices of the SCSBs to be reset if
%     the transition to the location "dest_loc(i)" is taken.
%   * "terminal" - A matrix of row vectors. Each row contains the FSMB
%     state vector corresponding to a terminal FSM state that can be
%     reached.

global GLOBAL_PIHA 

q = GLOBAL_PIHA.Locations{srcloc}.q;

% From the list of enabled transitions, find the unique combinations of the
% next state and reset flag for each FSMB.
dest = cell(1,length(q));
for k = 1:length(q)
  if isempty(enabled{k})
    % If no transition is enabled for the current FSMB, then the only
    % possible combination is that the next state is the same as the current
    % state and the reset flag is zero.
    dest{k} = [q(k) 0];
  else
    dest{k} = [];
    for i = 1:length(enabled{k})
      % Find the next q and reset flag for the current transition for the
      % current FSMB.
      next_q = transitions{enabled{k}(i)}.destination;
      next_reset = transitions{enabled{k}(i)}.reset_flag;
      next_combo = [next_q next_reset];
      % Verify that the next state and reset flag combo is unique before
      % including it in the destination list for the current FSMB.
      found = 0;
      for j = 1:size(dest{k},1)
        if all(next_combo == dest{k}(j,:))
          found = 1;
          break;
        end
      end
      if ~found
        dest{k} = [dest{k}; next_combo];
      end
    end
  end
end

% Compute all possible combinations (cross products of the next 'q' vectors
% and 'reset' vectors across all FSMBs. When the loop below finishes,
% "dest_q" and "dest_reset_fsm" are matrices of the same size. "dest_q(i,:)"
% and "dest_reset_fsm(i,:)" go together to represent a combination being
% generated below. For the i-th combination, "dest_q(i,j)" is the next state
% and "dest_reset_fsm(i,j)" is the reset flag for the j-th FSMB.
dest_q = dest{1}(:,1);
dest_reset_fsm = dest{1}(:,2);
for k = 2:length(q)
  dest_q_new = [];
  dest_reset_fsm_new = [];
  m = size(dest_q,1);
  for i = 1:size(dest{k},1)
    dest_q_new = [dest_q_new; [dest_q dest{k}(i,1)*ones(m,1)]];
    dest_reset_fsm_new = [dest_reset_fsm_new;
                    [dest_reset_fsm dest{k}(i,2)*ones(m,1)]];
  end
  dest_q = dest_q_new;
  dest_reset_fsm = dest_reset_fsm_new;
end
clear dest_q_new
clear dest_reset_fsm_new

% Get the list of SCSBs to be reset by each FSMB. This list is identical for
% all transitions under the same FSM, so we can just grab the list from the
% first transition found under each FSM.

for k = 1:length(q)
  if isempty(enabled{k})
    reset_scs_idx{k} = [];
  else
    reset_scs_idx{k} = transitions{enabled{k}(1)}.reset_scs_index;
  end
end

% Convert the reset flag vector for each FSMB in each row of
% "dest_reset_fsm" into an index vector for SCSB to be reset.
dest_reset_scs_idx = cell(size(dest_reset_fsm,1),1);
for i = 1:size(dest_reset_fsm,1)
  for j = 1:size(dest_reset_fsm,2)
    if dest_reset_fsm(i,j)
      dest_reset_scs_idx{i} = ...
          union(dest_reset_scs_idx{i},reset_scs_idx{j});
    end
  end
end

% Find the destination location in GLOBAL_PIHA that corresponds to each
% row of "dest_q".
dest_loc = [];
dest_reset = {};
terminal = [];
for i = 1:size(dest_q,1)
  found = 0;
  for j = 1:length(GLOBAL_PIHA.Locations)
    if all(GLOBAL_PIHA.Locations{j}.q == dest_q(i,:))
      found = 1;
      break;
    end
  end
  if found
    % If a matching location is found, add the location index and the
    % corresponding list of SCSB to be reset to output list.
    dest_loc = [dest_loc; j];
    dest_reset = [dest_reset; {dest_reset_scs_idx{i}}];
  else
    % If a matching location is not found, add the 'q' vector to the
    % terminal state list.
    terminal = [terminal; dest_q(i,:)];
  end
end

return

% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

function str = location_name(loc)

global GLOBAL_PIHA

% The state names stored in the .state field of the locations are
% character matrices, i.e. each row of .state contains the name
% of each component state.
name = GLOBAL_PIHA.Locations{loc}.state;
str = '(';
for k = 1:size(name,1)
  if (k > 1)
    str = [str ','];
  end
  str = [str deblank(name(k,:))];
end
str = [str ')'];

return

% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

function box = bounding_box(mapping)

% New version of this routine: uses the chosen option for hull computation to
% determine the bounding box.
%
% --OS--06/21/02

global GLOBAL_APPROX_PARAM

V = vertices;
if length(mapping) > 1
  for j = 1:length(mapping)
    vj = vertices(mapping{j});
    V = V | vj;
  end   
  box = linearcon(polyhedron(V));
else
%   Changed from
%  box = mapping{1};
%   to the following codes by Dong Jia, Mar. 27.
     [CE,dE,CI,dI]=linearcon_data(mapping{1});
%     if length(CE)==0
%         box=mapping{1};
%     else
         CI=[CI;CE;-CE];
         dI=[dI;[dE;-dE]+ones(2*length(dE),1)*GLOBAL_APPROX_PARAM.poly_bloat_tol];
         box=linearcon([],[],CI,dI);
         box=clean_up(box);
%     end
end
return
