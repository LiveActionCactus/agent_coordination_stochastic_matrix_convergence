function [dist_mat, adj_mat, deg_mat] = calc_distance_matrix(agent_posns, comm_dist)
%CALC_DISTANCE_MATRIX generates an nxn distance matrix and mask based on
%communication capabilities of each agent, given 2xn position matrix.

diag_num_agents = eye( length(agent_posns) );       % repeated operation on matrix diagonal

dist_mat = hypot( agent_posns(1,:)-agent_posns(1,:).', agent_posns(2,:)-agent_posns(2,:).');

adj_mat = ( dist_mat <= comm_dist ) - diag_num_agents;
deg_mat = diag_num_agents .* sum(adj_mat, 2);

end