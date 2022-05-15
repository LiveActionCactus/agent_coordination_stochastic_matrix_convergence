%% Consensus for agents on time-varying (undirected) graph -- naive-case 
% 
% By: Patrick Ledzian
% Date: 12 May 2022
%

% TODO:
% [DONE] 1) initialize agents within axis-aligned boundaries
% [DONE] 2) generate distance matrix between agents
% [DONE] 3) mask the distance matrix according to comms distances
% [DONE] 4) iterate the union of graphs
% [DONE] 5) generate Laplacian matrix and check for spanning tree
% [DONE] 6) multi-start and find average time-to-spanning tree 
% 7) check discrete simulation of average-consensus converges for a given
% sequence of graphs. (Not relationship between init pos and convergence... Fiedler eigenvalue?)
% 8) store and check all graphs making up 4) in every multi-start and see
% if they converge, analyze the ones that don't
%

clear all
close all
clc

%%
start_itrs = 1000;
sim_steps = 1000;               % upper bound on steps in single iteration

bounds = [0 30; 0 30];          % axis-aligned bounds
num_agents = 4;
comm_dist = 2.0;
trees = 5;                    % how many spanning trees need to occur before stopping

multi_start_results = cell(2,start_itrs);

for i = 1:start_itrs
    [multi_start_results{1,i}, multi_start_results{2,i}] = run_sim(sim_steps, bounds, num_agents, comm_dist, trees);
end

%% Sim spanning tree analysis

%
% Analysis:
% For axis aligned square boundaries the "time-to-spanning-tree" appears as
% a chi-squared distribution.
%

idxs = cat(1,multi_start_results{2,:});

ave = sum(idxs) / start_itrs;
std_dev = std(idxs);

figure()
subplot(1,2,1)
hist(idxs(:,1))
xlabel("Sim steps")
ylabel("Occurances")
title("Sim steps before first spanning tree; " + "worst case=" + max(idxs(:,1)) + " steps")

subplot(1,2,2)
hist(idxs(:,3))
xlabel("Sim steps")
ylabel("Occurances")
title("Sim steps before third spanning tree; " + "worst case=" + max(idxs(:,3)) + " steps")

sgtitle(num_agents + " agents; " + "boundary: " + bounds(1,2) + "x" + bounds(2,2))

%% Sim consensus analysis -- row stochastic (SIA) matrices approach
%
% ref: Consensus of Information Under Dynamically Changing Interaction
% Topologies; W. Ren
%
% Analysis:
% convergence of SIA matrix seems to be very sensitive to the uniform
% interval requirement. I should do more investigating.

%
% TODO
% 1) sum rows of final matrices see if all 1.0 (one necessary condition for SIA)
% 2) find eignevalues of final matrices, should only have 1 lambda=1
% (another necessary condition for SIA)
% 3) plot dynamics over time for some of them
%
% The SIA matrix converges asymptotically if there are successive time
% intervals where the union of graphs in each interval has a spanning tree

% WHAT IS THE CONNECTION BETWEEN THE SIA AND LAPLACIAN MATRICES??
% right now I find a spanning tree via the rank of the laplacian, but I am
% simulating the dynamics of an SIA matrix. If the laplacian is rank 3,
% then is the SIA matrix rank 3? EVIDENCE POINTS TO NO...
%
% Am I properly creating the SIA matrices from stored adjacency matrices?
%

S_tot = cell(1,start_itrs);                     % S{1,:} SIA matrices, S{2,:} sum of rows, S{3,:} eigenvalues
for itr = 1:start_itrs                          % simulation iteration to examine

    steps = multi_start_results{2, itr}(1,end);        % sim steps req'd for convergence
    S = cell(1,steps);                          

    for k = 1:steps
        G = multi_start_results{1, itr}{1,k}+eye(num_agents);
        G = G ./ sum(G, 2);
        if isequal(k,1)
            S{1,k} = G;
        else
            S{1,k} = G*S{1,(k-1)};
        end

    end % for steps in single sim iteration

    S_tot{1,itr} = S{1,steps};
    S_tot{2,itr} = sum(S{1,steps},2);
    S_tot{3,itr} = eig(S{1,steps});
    S_tot{4,itr} = S;

end % for each sim iteration

%% Plotting consensus results
%
% Analysis
% I was expecting to see only 1 zero-eigenvalue, as it seems to indicate in
% the paper. I'm not sure if there is a bug in my code or what.
%
% via the paper it seems the SIA condition is that must be only 1
% 1-eigenvalue, all others must be less.
%

figure()
bar(cat(1,S_tot{2,:}))
xlabel("Row sum (x" + num_agents + " agents)")
ylabel("Frequency")
title({"Row sum of stochastic matrices over all iterations", " (" ...
    + num_agents + " agents; " + bounds(1,2) + "x" + bounds(2,2) + " bounds; sim itrs=" + start_itrs + ")" })

check_eig = sum(abs(1-cat(2, S_tot{3,:})) < 1e-12, 1);

figure()
hist(check_eig)
xlabel("Number of 1-eigenvalues")
ylabel("Occurances")
title({"1-eigenvalue occurances in final SIA matrix", num_agents + " agents; " + "boundary: " + bounds(1,2) + "x" + bounds(2,2)})

%% Check sim convergence for all sim runs
%
% I think I could generate a distribution with # of sims in convergence vs
% comm step 
%

convergence = zeros(1,itr);
for r = 1:itr
    S_t = S_tot{4,r};
    len_S = length(S_t);
    init_pos = linspace(0,1,num_agents)';
    result = zeros(num_agents, len_S+1);
    result(:,1) = init_pos;

    for i = 1:len_S
        result(:,i+1) = S_t{1,i}*result(:,i);
    end
    
    convergence(1,r) = norm(result(:,len_S+1) - result(:,len_S));
    conv_f = sum(convergence > 0.05);
end

%% Plot sim convergence

figure()
for q = 1:num_agents
    k = randi(itr,1);
    S_t = S_tot{4,k};
    len_S = length(S_t);
    
    init_pos = linspace(0,1,num_agents)';
    result = zeros(num_agents, len_S+1);
    result(:,1) = init_pos;
    
    for i = 1:len_S
        result(:,i+1) = S_t{1,i}*result(:,i);
    end
    
    x = linspace(0, len_S+1, len_S+1);
    subplot(2,2,q)
    hold on
    for i=1:num_agents
        plot(x, result(i,:), "--o")
    end
    ylim([0.0,1.0])
    grid on
    
    spans = multi_start_results{2,k};
    for k=1:trees
        xline(spans(1,k), "-.", "stree")
    end
    
    title(bounds(1,2) + "x" + bounds(2,2) + " bounds; " + num_agents + " agents; " + trees + " spanning trees")
    xlabel("Communication step / re-initialization of agent")
    ylabel("Information value")
end

sgtitle("DT information consensus vs communication step")