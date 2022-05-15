function init_pos = initialize_agents(num_agents, bounds)
%INITIALIZE_AGENTS inializes agents randomly inside axis-aliged, convex 
% boundaries.

init_pos = bounds(:,1) + (bounds(:,2)-bounds(:,1)).*rand(2,num_agents);

end