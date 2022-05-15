function [matrices, idx] = run_sim(sim_steps, bounds, num_agents, comm_dist, trees)
%RUN_SIM run a single iteration of simulation, breaks when spanning tree
%first found.

    matrices = cell(2,sim_steps);                % row 1: adjacency matrices; row 2: degree matrices
    G = eye(num_agents);                       % initialize the G matrix
    idx = zeros(1,trees);                          % indicies for each tree
    t = 0;                                       % tree counter
    
    % TODO: I THINK I AM CALCULATING THE UNION OF GRAPHS INCORRECTLY!
    % SHOULD BE SUBTRACTING OUT THE INTERSECTION... DO THE 1s ON THE
    % DIAGONAL SUM? 
    % It doesn't seem like the union operation should be creating weighted
    % edges... something to explore though. Let's do max(0,1) for now
    for i = 1:sim_steps
        init_pos = initialize_agents(num_agents, bounds);
        [dist_mat, adj_mat, deg_mat] = calc_distance_matrix(init_pos, comm_dist);
        matrices{1,i} = adj_mat;
        matrices{2,i} = deg_mat;
    
        G = max(G,adj_mat);
        G_r = rank(G.*-eye(num_agents));
        matrices{3,i} = G_r;
        %matrices{4,i} = deg_mat-adj_mat

        % TODO: WE'RE GETTING STUCK WITH RANK 4 MATRICES; DOES THIS MEAN
        % SPANNING TREE AS WELL? OR AM I DOING SOMETHING WRONG WITH THE
        % CHECK? OR WRONG WITH WHAT TO CHECK?
        %
        % I THINK I SHOULD BE CHECKING THE STOCHASTIC DIRECTED MATRIX D FOR
        % A SPANNING TREE; SEE THEOREM 3.2; SPANNING TREE FOR G IS
        % NECESSARY CONDITION FOR SPANNING TREE FOR D
        %if isequal(rank(L), num_agents-1)    % check for spanning tree
        if rank(G) >= num_agents-1
            t = t + 1;
            idx(1,t) = i;
            test = G;
            G = zeros(num_agents);

            if isequal(t, trees)
                matrices{3,i} = G_r;
                break
            end
        end % end if rank denotes spanning tree
    
    end % end for sim_steps 

%     if ~isequal(sum(any(idx,1)), 3)
%         disp(test)
%         disp(L_r)
%         disp(num_agents)
%         disp(t)
%         disp(idx)
%         disp(i)
%         assignin('base', 'test', test)
%         assignin('base', 'L', L)
%         assignin('base', 'adj', adj_mat)
%         pause()
%     end

end % end run_sim()