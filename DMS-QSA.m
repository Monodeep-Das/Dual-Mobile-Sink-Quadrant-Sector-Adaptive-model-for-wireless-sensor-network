% =========================================================================
% IEEE-Quality WSN Energy Simulation Framework
% WITH MATHEMATICAL FORMULATION (MTE, MRE, PLBTC, EEDRB, DMS-QSA, SPIRAL METHOD)
% =========================================================================
clear; clc; close all;

%% 1. Network & Radio Model Parameters
num_nodes = 100;
field_x = 100; field_y = 100;
initial_energy = 0.5;
packet_size = 4000;
max_rounds = 3000;
N_runs = 100; % Monte Carlo runs

E_elec = 50e-9;
E_fs = 10e-12;
E_mp = 0.0013e-12;
d0 = sqrt(E_fs/E_mp);
E_rx = E_elec;
E_da = 5e-9;
max_hop_range = 40;

%% 2. Legacy Sink Trajectory
% A=(0,100), B=(100,100), C=(100,0), D=(0,0)
waypoints = [0,100; 100,100; 100,0; 0,0; 0,100];
num_waypoints = size(waypoints,1);
segment_dists = sqrt(sum(diff(waypoints).^2,2));
total_dist = sum(segment_dists);
cum_dist = [0;cumsum(segment_dists)];

%% Data Storage (Updated for 6 protocols to include SPIRAL METHOD)
alive_runs = zeros(6,max_rounds,N_runs);
energy_runs = zeros(6,max_rounds,N_runs);
throughput_runs = zeros(6,max_rounds,N_runs);

protocol_names = {'MTE','MRE','PLBTC','Proposed EEDRB', 'Proposed DMS-QSA', 'SPIRAL METHOD'};

%% ================= MONTE CARLO =================
for run = 1:N_runs
    fprintf('Run %d/%d\n',run,N_runs);
    rng(42+run); % Ensure repeatable but distinct layouts per run
   
    nodes_x = rand(1,num_nodes)*field_x;
    nodes_y = rand(1,num_nodes)*field_y;
    [X,Y] = meshgrid(nodes_x,nodes_y);
    dist_matrix = sqrt((X-X').^2+(Y-Y').^2);
   
    % Give nodes a wide variance in initial energy to prevent simultaneous deaths
    real_initial_energy = initial_energy .* (0.5 + (0.5 * rand(1, num_nodes)));
    protocol_energies = repmat({real_initial_energy},1,6);
   
    %% INITIALIZE DUAL SINKS (NOW FOR DMS-QSA)
    A = [0,100]; B = [100,100]; C = [100,0]; D = [0,0];
    P = [50,100]; Q = [100,50]; R = [50,0]; S = [0,50];
    O = [50,50];
   
    % Sink 1: START Q → B
    path1 = [Q; B; P; O; S; A; P; O; Q];
    % Sink 2: START D → S
    path2 = [D; S; O; R; C; Q; O; R; D];
   
    idx1 = 1; idx2 = 1;
    pos1 = path1(1,:);
    pos2 = path2(1,:);
    speed = 2;
   
    pause_duration = 5; 
    pause_points = [P; Q; R; S; O];
    pause_counter1 = pause_duration; % Sink 1 starts at Q
    pause_counter2 = 0; % Sink 2 starts at D
   
    %% ROUNDS
    for r = 1:max_rounds
        %% Mobile Sink Position (Legacy - used for MTE, MRE, PLBTC, and now EEDRB)
        cycle_length = 200;
        dist_cycle = mod(r,cycle_length)/cycle_length * total_dist;
        seg_idx = find(cum_dist <= dist_cycle, 1, 'last');
        if seg_idx == num_waypoints, seg_idx = seg_idx - 1; end
       
        t_seg = (dist_cycle - cum_dist(seg_idx)) / segment_dists(seg_idx);
        bs_x = waypoints(seg_idx,1) + t_seg*(waypoints(seg_idx+1,1)-waypoints(seg_idx,1));
        bs_y = waypoints(seg_idx,2) + t_seg*(waypoints(seg_idx+1,2)-waypoints(seg_idx,2));
       
        % Base distance array for legacy sink
        dist_to_bs = sqrt((nodes_x - bs_x).^2 + (nodes_y - bs_y).^2);
       
        %% RUN PROTOCOLS
        for p_idx = 1:6
            energy = protocol_energies{p_idx};
            alive_idx = find(energy > 0);
            if isempty(alive_idx), continue; end
           
            bs_idx = num_nodes + 1;
            parent = zeros(1,num_nodes);
           
            %% ==========================================================
            %% ROUTING (MTE, MRE, PLBTC)
            %% ==========================================================
            if p_idx <= 3
                for i = alive_idx
                    best_cost = Inf;
                    best_j = -1;
                   
                    for j = alive_idx
                        if j == i, continue; end
                        d = dist_matrix(i,j);
                        if d > max_hop_range, continue; end
                        if dist_to_bs(j) >= dist_to_bs(i), continue; end
                       
                        Etx_ij = calculate_tx_energy(packet_size, d, d0, E_elec, E_fs, E_mp);
                       
                        switch p_idx
                            case 1 % MTE
                                cost = Etx_ij;
                            case 2 % MRE
                                cost = Etx_ij / sqrt(energy(j));
                            case 3 % PLBTC
                                Etx_j_bs = calculate_tx_energy(packet_size, dist_to_bs(j), d0, E_elec, E_fs, E_mp);
                                cost = (E_rx + Etx_j_bs + E_da*packet_size) / sqrt(energy(j));
                        end
                       
                        if cost < best_cost
                            best_cost = cost;
                            best_j = j;
                        end
                    end
                   
                    %% Direct transmission to sink
                    Etx_bs = calculate_tx_energy(packet_size, dist_to_bs(i), d0, E_elec, E_fs, E_mp);
                    if p_idx == 1
                        cost_bs = Etx_bs;
                    elseif p_idx == 2
                        cost_bs = Etx_bs / energy(i);
                    elseif p_idx == 3
                        cost_bs = (Etx_bs + E_da*packet_size) / energy(i);
                    end
                   
                    if cost_bs < best_cost || best_j == -1
                        best_j = bs_idx;
                    end
                    parent(i) = best_j;
                end
               
            elseif p_idx == 4
                %% ======================================================
                %% PROPOSED EEDRB (NOW USING SINGLE PERIMETER SINK)
                %% ======================================================
                % Uses the default dist_to_bs calculated at start of round
                
                parent = zeros(1,num_nodes);
                for i = alive_idx
                    best_cost = Inf;
                    best_j = -1;
                   
                    current_max_hop = max_hop_range;
                    if energy(i) > (initial_energy * 0.3) % If battery is above 30%
                        current_max_hop = max_hop_range * 1.5; % Expand range to 60m
                    end
                   
                    for j = alive_idx
                        if i == j, continue; end
                        d = dist_matrix(i,j);
                       
                        if d > current_max_hop, continue; end
                        if dist_to_bs(j) >= dist_to_bs(i), continue; end
                       
                        Etx_ij = calculate_tx_energy(packet_size, d, d0, E_elec, E_fs, E_mp);
                        Etx_j_bs = calculate_tx_energy(packet_size, dist_to_bs(j), d0, E_elec, E_fs, E_mp);
                        E_j_process = E_rx + (E_da * packet_size);
                       
                        cost = (Etx_ij + E_j_process + Etx_j_bs) / sqrt(energy(j));
                       
                        if cost < best_cost
                            best_cost = cost;
                            best_j = j;
                        end
                    end
                   
                    Etx_bs = calculate_tx_energy(packet_size, dist_to_bs(i), d0, E_elec, E_fs, E_mp);
                    cost_bs = (Etx_bs + E_da*packet_size) / energy(i);
                   
                    if cost_bs < best_cost || best_j == -1
                        best_j = bs_idx;
                    end
                    parent(i) = best_j;
                end
               
            elseif p_idx == 5
                %% ======================================================
                %% PROPOSED DMS-QSA (NOW USING DUAL GRID-BASED SINKS)
                %% ======================================================
                % Move Sink 1 with Pause Logic
                if pause_counter1 > 0
                    pause_counter1 = pause_counter1 - 1;
                else
                    target1 = path1(idx1+1,:);
                    dir1 = target1 - pos1;
                    if norm(dir1) <= speed
                        pos1 = target1;
                        idx1 = idx1 + 1;
                        if idx1 >= size(path1,1), idx1 = 1; end
                       
                        for k = 1:size(pause_points,1)
                            if norm(pos1 - pause_points(k,:)) < 0.1
                                pause_counter1 = pause_duration;
                                break;
                            end
                        end
                    else
                        dir1 = dir1 / norm(dir1);
                        pos1 = pos1 + dir1 * speed;
                    end
                end
                ms1_x = pos1(1); ms1_y = pos1(2);
               
                % Move Sink 2 with Pause Logic
                if pause_counter2 > 0
                    pause_counter2 = pause_counter2 - 1;
                else
                    target2 = path2(idx2+1,:);
                    dir2 = target2 - pos2;
                    if norm(dir2) <= speed
                        pos2 = target2;
                        idx2 = idx2 + 1;
                        if idx2 >= size(path2,1), idx2 = 1; end
                       
                        for k = 1:size(pause_points,1)
                            if norm(pos2 - pause_points(k,:)) < 0.1
                                pause_counter2 = pause_duration;
                                break;
                            end
                        end
                    else
                        dir2 = dir2 / norm(dir2);
                        pos2 = pos2 + dir2 * speed;
                    end
                end
                ms2_x = pos2(1); ms2_y = pos2(2);
               
                % Update distance array for DMS-QSA dual sinks
                for i = alive_idx
                    d1 = sqrt((nodes_x(i)-ms1_x)^2 + (nodes_y(i)-ms1_y)^2);
                    d2 = sqrt((nodes_x(i)-ms2_x)^2 + (nodes_y(i)-ms2_y)^2);
                    dist_to_bs(i) = min(d1, d2);
                end
               
                parent = zeros(1,num_nodes);
                for i = alive_idx
                    best_cost = Inf;
                    best_j = -1;
                   
                    for j = alive_idx
                        if i == j, continue; end
                        d = dist_matrix(i,j);
                       
                        if d > max_hop_range, continue; end
                        if dist_to_bs(j) >= dist_to_bs(i), continue; end % Guarantee DAG for energy flow calc
                       
                        % Orthogonal Relay Rule: Move closer to axis crosshairs
                        dist_j_axis = min(abs(nodes_x(j)-50), abs(nodes_y(j)-50));
                        dist_i_axis = min(abs(nodes_x(i)-50), abs(nodes_y(i)-50));
                       
                        if dist_j_axis >= dist_i_axis, continue; end
                       
                        cost = dist_j_axis / (energy(j) + 1e-4);
                       
                        if cost < best_cost
                            best_cost = cost;
                            best_j = j;
                        end
                    end
                   
                    Etx_bs = calculate_tx_energy(packet_size, dist_to_bs(i), d0, E_elec, E_fs, E_mp);
                    cost_bs = (Etx_bs + E_da*packet_size) / energy(i);
                   
                    if cost_bs < best_cost || best_j == -1
                        best_j = bs_idx;
                    end
                    parent(i) = best_j;
                end
               
            elseif p_idx == 6
                %% ======================================================
                %% NEW: SPIRAL METHOD (INWARD SPIRAL SINK)
                %% ======================================================
                % Sink starts at edge and spirals inwards to the center
                rotations = 5; % 5 full rotations over the simulation
                theta_spiral = (r / max_rounds) * 2 * pi * rotations;
                rad_spiral = 50 * (1 - r / max_rounds); % Decreasing radius
               
                spiral_x = 50 + rad_spiral * cos(theta_spiral);
                spiral_y = 50 + rad_spiral * sin(theta_spiral);
               
                % Update distance array for SPIRAL METHOD
                for i = alive_idx
                    dist_to_bs(i) = sqrt((nodes_x(i)-spiral_x)^2 + (nodes_y(i)-spiral_y)^2);
                end
               
                parent = zeros(1,num_nodes);
                for i = alive_idx
                    best_cost = Inf;
                    best_j = -1;
                   
                    for j = alive_idx
                        if i == j, continue; end
                        d = dist_matrix(i,j);
                       
                        if d > max_hop_range, continue; end
                        if dist_to_bs(j) >= dist_to_bs(i), continue; end
                       
                        Etx_ij = calculate_tx_energy(packet_size, d, d0, E_elec, E_fs, E_mp);
                        Etx_j_bs = calculate_tx_energy(packet_size, dist_to_bs(j), d0, E_elec, E_fs, E_mp);
                        E_j_process = E_rx + (E_da * packet_size);
                       
                        % Cost dynamically weighting distance to spiral sink and residual energy
                        cost = (Etx_ij + E_j_process + Etx_j_bs) / sqrt(energy(j));
                       
                        if cost < best_cost
                            best_cost = cost;
                            best_j = j;
                        end
                    end
                   
                    Etx_bs = calculate_tx_energy(packet_size, dist_to_bs(i), d0, E_elec, E_fs, E_mp);
                    cost_bs = (Etx_bs + E_da*packet_size) / energy(i);
                   
                    if cost_bs < best_cost || best_j == -1
                        best_j = bs_idx;
                    end
                    parent(i) = best_j;
                end
            end
           
            %% ==========================================================
            %% ENERGY CALCULATION
            %% ==========================================================
            % Nodes have an 85% chance to transmit, 15% chance to be idle/drop packet
            packets_tx = double(rand(1,num_nodes) <= 0.85);
            [~, order] = sort(dist_to_bs(alive_idx), 'descend');
            sorted_alive = alive_idx(order);
           
            for i = sorted_alive
                par = parent(i);
                if par <= num_nodes && par > 0
                    packets_tx(par) = packets_tx(par) + packets_tx(i);
                end
            end
           
            packets_to_sink = 0;
            for i = alive_idx
                par = parent(i);
                if par == bs_idx
                    d = dist_to_bs(i);
                    packets_to_sink = packets_to_sink + packets_tx(i);
                else
                    d = dist_matrix(i, par);
                end
               
                Etx = calculate_tx_energy(packet_size, d, d0, E_elec, E_fs, E_mp);
                drain = Etx + (packets_tx(i)-1)*E_rx + packets_tx(i)*packet_size*E_da;
                energy(i) = max(0, energy(i) - drain);
            end
           
            protocol_energies{p_idx} = energy;
            alive_runs(p_idx,r,run) = sum(energy > 0);
            energy_runs(p_idx,r,run) = sum(energy);
            throughput_runs(p_idx,r,run) = packets_to_sink;
        end
    end
end

%% ================= 3. MAIN SIMULATION RESULTS (PLOTS) =================
alive_mean = mean(alive_runs,3);
energy_mean = mean(energy_runs,3);
throughput_mean = mean(throughput_runs,3);
dead_mean = num_nodes - alive_mean;

% Set up a strong smoothing window to create gradual curves
smooth_window = 150;
styles = {':k','-.r','--g','-c','-b', '-m'}; % Added solid magenta for SPIRAL METHOD
widths = [2 2 2 2 3 3]; % Ensure both proposed methods and Spiral stand out

figure; hold on; grid on;
for p = 1:6
    smooth_dead = smoothdata(dead_mean(p,:), 'gaussian', smooth_window);
    plot(smooth_dead, styles{p}, 'LineWidth', widths(p));
end
title('Dead Nodes vs Rounds');
legend(protocol_names);
xlabel('Rounds'); ylabel('Number of Dead Nodes');

figure; hold on; grid on;
for p = 1:6
    smooth_energy = smoothdata(energy_mean(p,:), 'gaussian', smooth_window);
    plot(smooth_energy, styles{p}, 'LineWidth', widths(p));
end
title('Residual Energy vs Rounds');
legend(protocol_names);
xlabel('Rounds'); ylabel('Total Residual Energy (J)');

figure; hold on; grid on;
for p = 1:6
    smooth_throughput = smoothdata(throughput_mean(p,:), 'gaussian', smooth_window);
    plot(smooth_throughput, styles{p}, 'LineWidth', widths(p));
end
title('Throughput');
legend(protocol_names, 'Location', 'northwest');
xlabel('Rounds'); ylabel('Packets Delivered to Sink');

%% ================= 4. PDR (Packet Delivery Ratio) GRAPH =================
% Calculate expected packets generated based on 85% transmission probability
expected_generated = alive_mean * 0.85;

% Calculate PDR (%) safely avoiding division by zero
PDR_mean = zeros(size(throughput_mean));
valid_idx = expected_generated > 0;
PDR_mean(valid_idx) = (throughput_mean(valid_idx) ./ expected_generated(valid_idx)) * 100;

figure('Name', 'PDR Performance'); hold on; grid on;
for p = 1:6
    % Smooth the PDR curve for visual clarity (matches previous graphs)
    smooth_PDR = smoothdata(PDR_mean(p,:), 'gaussian', smooth_window);
    
    % Prevent PDR from exceeding 100% due to smoothing artifacts
    smooth_PDR(smooth_PDR > 100) = 100;
    % Prevent dropping below 0
    smooth_PDR(smooth_PDR < 0) = 0; 
    
    plot(smooth_PDR, styles{p}, 'LineWidth', widths(p));
end
title('Packet Delivery Ratio (PDR) vs Rounds');
legend(protocol_names, 'Location', 'southwest');
xlabel('Rounds'); 
ylabel('Packet Delivery Ratio (%)');
ylim([0 105]); % Cap Y-axis slightly above 100% for readability

%% ================= 5. PERFORMANCE SUMMARY TABLE =================
% Initialize matrices to hold metrics for every single run
num_protocols = 6;
fnd_runs = zeros(num_protocols, N_runs);
hnd_runs = zeros(num_protocols, N_runs);
cum_pkts_runs = zeros(num_protocols, N_runs);
peak_pdr_runs = zeros(num_protocols, N_runs); 

for p = 1:num_protocols
    for run = 1:N_runs
        % 1. Extract alive nodes and throughput for this specific run
        alive_in_run = alive_runs(p, :, run);
        throughput_in_run = throughput_runs(p, :, run);
        
        % 2. Calculate FND (First Node Dead)
        dead_idx = find(alive_in_run < num_nodes, 1);
        if isempty(dead_idx)
            fnd_runs(p, run) = max_rounds; % Node never died
        else
            fnd_runs(p, run) = dead_idx;
        end
        
        % 3. Calculate HND (Half Nodes Dead)
        half_dead_idx = find(alive_in_run <= num_nodes/2, 1);
        if isempty(half_dead_idx)
            hnd_runs(p, run) = max_rounds;
        else
            hnd_runs(p, run) = half_dead_idx;
        end
        
        % 4. Calculate Cumulative Packets
        cum_pkts_runs(p, run) = sum(throughput_in_run);
        
        % 5. Calculate Peak PDR (%) for this run
        % Expected packets generated (85% probability per alive node)
        exp_pkts = alive_in_run * 0.85; 
        pdr_array = zeros(1, max_rounds);
        valid = exp_pkts > 0;
        pdr_array(valid) = (throughput_in_run(valid) ./ exp_pkts(valid)) * 100;
        
        % Cap at 100% to handle minor statistical overshoots in probability
        pdr_array(pdr_array > 100) = 100; 
        peak_pdr_runs(p, run) = max(pdr_array);
    end
end

% Compute Means and Standard Deviations across all N_runs
fnd_mean = mean(fnd_runs, 2);
fnd_std = std(fnd_runs, 0, 2);
hnd_mean = mean(hnd_runs, 2);
hnd_std = std(hnd_runs, 0, 2);
peak_pdr_mean = mean(peak_pdr_runs, 2);
cum_pkts_mean = mean(cum_pkts_runs, 2) / 100000; % Scale to 10^5

% Clean up protocol names for the table output to match Table II exactly
display_names = {'MTE', 'MRE', 'PLBTC', 'EEDRB', 'DMS-QSA', 'Spiral'};

%% ---------------- PRINT TABLE COMMAND WINDOW ----------------
fprintf('\nTABLE II\n');
fprintf('PERFORMANCE SUMMARY (MEAN ± STD DEV)\n');
fprintf('---------------------------------------------------------------------\n');
fprintf('%-12s %-12s %-12s %-15s %-15s\n', 'Protocol', 'FND', 'HND', 'Peak PDR (%)', 'Cum. Pkts ×10^5');
fprintf('---------------------------------------------------------------------\n');

for p = 1:num_protocols
    % Handle the custom display names array (re-mapping indices to match Table order)
    if p == 5
        disp_name = display_names{6}; % Print Spiral 5th
        actual_p_idx = 6;
    elseif p == 6
        disp_name = display_names{5}; % Print DMS-QSA 6th
        actual_p_idx = 5;
    else
        disp_name = display_names{p};
        actual_p_idx = p;
    end
    
    % Format standard deviation strings (e.g., 295±38)
    fnd_str = sprintf('%.0f±%.0f', fnd_mean(actual_p_idx), fnd_std(actual_p_idx));
    hnd_str = sprintf('%.0f±%.0f', hnd_mean(actual_p_idx), hnd_std(actual_p_idx));
    
    % Format standard decimals
    pdr_str = sprintf('%.1f', peak_pdr_mean(actual_p_idx));
    pkts_str = sprintf('%.2f', cum_pkts_mean(actual_p_idx));
    
    % Print row
    fprintf('%-12s %-12s %-12s %-15s %-15s\n', disp_name, fnd_str, hnd_str, pdr_str, pkts_str);
end
fprintf('---------------------------------------------------------------------\n');

%% =========================================================================
%% 6. MATHEMATICAL FORMULATION: SPIRAL METHOD EVALUATION
%% =========================================================================
%% 1. Simulation Parameters (From Table 2)
M = 200; % Network area side length in meters (200x200m^2)
E_initial = 5000; % Initial energy (5 kJ)
e_elec_math = 50e-9; % Energy for electronics circuit (50 nJ/bit)
e_amp_math = 10e-12; % Energy for amplification (10 pJ/bit/m^2)
Rt = 50; % Transmission range (50 m)
R_sense = 25; % Sensing range (25 m)
K = 0.1; % Information generation rate (0.1 bits/sec)
pkt_size_math = 500; % Packet size in bits
L = 3; % Optimal Sub-grid side length chosen from empirical results

% Test scenarios for varying sensor nodes
N_total_array = 100:100:500;
network_lifetimes = zeros(size(N_total_array));
total_energy_consumptions = zeros(size(N_total_array));

%% 2. Grid Construction (Equations 2 & 3)
% Calculate grid division factor (beta)
beta = (sqrt(5) * M) / Rt;

% Determine Grid Size (G x G) based on beta
if beta <= 6
    G = 6;
elseif beta > 6 && beta <= 9
    G = 9;
elseif beta > 9 && beta <= 12
    G = 12;
elseif beta > 12 && beta <= 15
    G = 15;
else
    G = 18;
end

%% 3. Sub-grid and Corona Architecture Mapping (Equations 4, 6, 7)
% Total number of sub-grid coronas
N_sgc = (G * G) / (L * L);
% Radius of innermost corona
R1 = Rt / N_sgc;
% Calculate number of coronas (n) in each L x L sub-grid
n = floor((3 * M * N_sgc) / (2 * G * Rt));
l_n = n;

% Calculate constant k for PDF (Equation 10)
sum_pdf_term = 0;
for i = 1:n
    sum_pdf_term = sum_pdf_term + ((2*i - 1)^2) / (i^4);
end
k_pdf = (l_n^2) / (pi * (R1^2) * sum_pdf_term);

%% 4. Node Distribution and Energy Evaluation Loop
for idx = 1:length(N_total_array)
    N_total = N_total_array(idx);
   
    % Total number of sensor nodes in one sub-grid (Equation 8)
    N_subgrid = N_total / N_sgc;
   
    % Arrays to store properties for each layer i
    Pi = zeros(1, n);
    Ni = zeros(1, n);
    Di = zeros(1, n);
    Ri = zeros(1, n);
    E_TX = zeros(1, n);
    E_RX = zeros(1, n);
    E_total_i = zeros(1, n);
    LT_i = zeros(1, n);
   
    % Energy consumed to transmit/receive a single bit over transmission range
    e_T = (pkt_size_math * e_elec_math) + (pkt_size_math * e_amp_math * (Rt^2));
    e_R = pkt_size_math * e_elec_math;
   
    % Step A: Calculate Deployment Densities (Equations 9, 11, 12)
    for i = 1:n
        Ri(i) = i * R1;
        % Probability of deployment in i-th layer
        Pi(i) = (k_pdf * ((2*i - 1)^2) * pi * (R1^2)) / ((l_n^2) * (i^4));
       
        % Number of nodes in i-th layer
        Ni(i) = Pi(i) * N_subgrid;
       
        % Node density in i-th layer
        Di(i) = Ni(i) / (pi * (Ri(i)^2));
    end
   
    % Step B: Calculate Energy Consumption Rate (ECR) per node (Equations 18, 19, 24, 25)
    for i = 1:n
        if i == n
            % Outermost corona only transmits its own sensed data
            E_total_i(i) = (K * e_T) / Di(i);
        else
            % Inner coronas transmit own data + forward outer corona data
            sum_Dx = 0;
            for x = i+1:n
                sum_Dx = sum_Dx + Di(x);
            end
           
            % E_TX + E_RX
            E_total_i(i) = (K * (e_T + (e_R * sum_Dx))) / Di(i);
        end
       
        % Lifetime of individual sensor node in i-th corona (Equation 29)
        LT_i(i) = E_initial / E_total_i(i);
    end
   
    % Step C: Calculate Network Lifetime & Total ECR (Equations 26, 30)
    % Network lifetime is determined by the corona that depletes energy fastest
    network_lifetimes(idx) = min(LT_i);
   
    % Total Energy Consumption Rate for the entire network
    total_ecr = 0;
    for j = 1:N_sgc
        for i = 1:n
            total_ecr = total_ecr + (Ni(i) * E_total_i(i));
        end
    end
    total_energy_consumptions(idx) = total_ecr;
end

%% 5. Plotting Results for SPIRAL METHOD Math Evaluation
figure('Name', 'SPIRAL METHOD Performance Metrics');

% Plot 1: Network Lifetime vs Number of Sensor Nodes
subplot(1,2,1);
semilogy(N_total_array, network_lifetimes, '-ms', 'LineWidth', 2, 'MarkerSize', 8); % Magenta to match main graph
grid on;
xlabel('Number of Sensor Nodes');
ylabel('Network Lifetime (in Rounds)');
title('Effect of Number of Sensor Nodes on Network Lifetime');
legend('SPIRAL METHOD', 'Location', 'best');

% Plot 2: Total Energy Consumption vs Number of Sensor Nodes
subplot(1,2,2);
semilogy(N_total_array, total_energy_consumptions, '-ms', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Number of Sensor Nodes');
ylabel('Energy Consumption in Joules');
title('Effect of Node Count on Energy Consumption');
legend('SPIRAL METHOD', 'Location', 'best');

%% ================= MAIN TX FUNCTION =================
function E_tx = calculate_tx_energy(k, d, d0, E_elec, E_fs, E_mp)
    if d < d0
        E_tx = k*E_elec + k*E_fs*d^2;
    else
        E_tx = k*E_elec + k*E_mp*d^4;
    end
end

% =========================================================================
% APPENDIX: FULL NEW DMS-QSA APPROACH (STANDALONE SIMULATION)
% =========================================================================
%% ========================= 1. PARAMETERS =========================
x = 400;
density = 0.0025;
N_per_quad = round(density * x^2);
total_nodes = 4 * N_per_quad; % 1600
node_range = 60; % member-to-CH radio range [m]
ch_range = 200; % CH-to-CH / CH-to-sink range [m]
sink_speed = 5; % [m/s]
sojourn_time = 30; % [s]
num_cycles = 50;
dt = 5; % simulation time step [s]
T_data = 15; % member data generation interval [s]
buf_max = 20; % CH buffer capacity [packets]
packet_size_bits = 4000;

E_elec = 50e-9; % [J/bit]
E_fs = 10e-12; % [J/bit/m^2]
E_mp = 0.0013e-12; % [J/bit/m^4]
E_DA = 5e-9; % [J/bit] aggregation energy
E_rx = E_elec; % [J/bit] receive energy
E_th = sqrt(E_fs/E_mp);% crossover distance ~87.7 m
E_init = 2; % [J] base initial energy (normal nodes)
E_sleep_th = 0.05*E_init; % sleep threshold = 5% of E_init

% SINR model constants
P_tx_ref = 1; % normalised transmit power
N0 = 1e-12; % noise floor

% FIX 4: Sink-side reception energy cost
E_sink_rx = E_rx; % [J/bit] sink receive energy (same as node Rx)

% -------------------------------------------------------------------------
% UPGRADE 1: Probabilistic Link Model parameters
SINR_th_dB = 10; % [dB] threshold (same value as before, now in dB)
k_logistic = 0.3; % FIX 3: reduced from 0.7 -> smoother sigmoid

%% =================== 2. RANDOM NODE DEPLOYMENT ===================
rng(42);
pos_all = rand(total_nodes,2) .* [2*x, 2*x];

% -------------------------------------------------------------------------
% UPGRADE 2: Heterogeneous Node Energy Model
type_rand = rand(total_nodes, 1); % uniform random for type draw
nodes = struct();
for i = 1:total_nodes
    nodes(i).id = i;
    nodes(i).pos = pos_all(i,:);
    nodes(i).cluster_id = 0;
    nodes(i).is_CH = false;
    nodes(i).is_backup = false;
    nodes(i).primary_CH_id = 0;
    nodes(i).backup_CH_id = 0;
    nodes(i).is_sleeping = false;
    nodes(i).buffer = 0;
    nodes(i).buf_time_sum = 0;
    
    % Assign node type and its specific initial energy
    if type_rand(i) <= 0.70
            nodes(i).type = 'normal';
            nodes(i).E_node_init = E_init; % 2 J
    elseif type_rand(i) <= 0.90
            nodes(i).type = 'advanced';
            nodes(i).E_node_init = 1.5 * E_init; % 3 J
    else
            nodes(i).type = 'super';
            nodes(i).E_node_init = 2.0 * E_init; % 4 J
    end
    nodes(i).energy = nodes(i).E_node_init; % start fully charged
end

% Print heterogeneous energy distribution summary
n_normal = sum(strcmp({nodes.type}, 'normal'));
n_advanced = sum(strcmp({nodes.type}, 'advanced'));
n_super = sum(strcmp({nodes.type}, 'super'));
fprintf('=== Node Type Distribution ===\n');
fprintf(' Normal (%5.1f%%): %4d nodes @ %.1f J\n', 100*n_normal/total_nodes, n_normal, E_init);
fprintf(' Advanced (%5.1f%%): %4d nodes @ %.1f J\n', 100*n_advanced/total_nodes, n_advanced, 1.5*E_init);
fprintf(' Super (%5.1f%%): %4d nodes @ %.1f J\n', 100*n_super/total_nodes, n_super, 2.0*E_init);
fprintf(' Total initial energy: %.1f J\n\n', sum([nodes.E_node_init]));

%% ========================= 3. RP DEFINITIONS =========================
rp_labels = {'P','M','S','T','U','W','X','Z','R','Y','Q','V','F','H','D','G'};
rp_pos = [
    x/2, 0; 3*x/2, 0;
    0, x/2; 0, 3*x/2;
    x/2, 3*x/2; 3*x/2, 3*x/2;
    2*x, 3*x/2; 2*x, x/2;
    x/2, x; 3*x/2, x;
    x, x/2; x, 3*x/2;
    x, 3*x/2; 2*x, x;
    0, x; 2*x, 2*x
];
num_rp = size(rp_pos,1);

%% ====================== 4. SINK PATHS ======================
s1_labels = {'S','R','Q','M','Z','Y','Q','P'};
s2_labels = {'X','W','V','R','T','U','V','Y'};
s1_idx = cellfun(@(l) find(strcmp(rp_labels,l)), s1_labels);
s2_idx = cellfun(@(l) find(strcmp(rp_labels,l)), s2_labels);

s1_path = rp_pos(s1_idx,:);
s2_path = rp_pos(s2_idx,:);

T_cycle = max(cycle_dur(s1_path, sink_speed, sojourn_time), ...
              cycle_dur(s2_path, sink_speed, sojourn_time));
t_vec = 0:dt:T_cycle;
n_steps = length(t_vec);

%% ================== 5. INITIAL CLUSTER ASSIGNMENT ==================
for i = 1:total_nodes
    d = sqrt(sum((rp_pos - nodes(i).pos).^2, 2));
    [~, nodes(i).cluster_id] = min(d);
end

%% ========================= 6. MAIN LOOP =========================
alive_log = zeros(1, num_cycles);
energy_log = zeros(total_nodes, num_cycles);
pkt_log = zeros(1, num_cycles);
delay_log = zeros(1, num_cycles);
pkt_loss_log = zeros(1, num_cycles); % probabilistic link losses
sink_rx_e_log = zeros(1, num_cycles); % FIX 4: sink Rx energy per cycle [J]

for cycle = 1:num_cycles
    % -------------------------------------------------------------------
    % UPGRADE 3: Multi-Iteration K-Means Clustering
    if mod(cycle, 2) == 1
            nodes = adaptive_cluster_kmeans(nodes, total_nodes, num_rp, rp_pos, 5);
    end
    
    mean_e_ratio = mean([nodes.energy]) / E_init;
    nodes = elect_CH_adaptive(nodes, num_rp, E_init, node_range, mean_e_ratio);
    nodes = promote_backup(nodes);
    
    CH_alive = [nodes.is_CH] & ([nodes.energy] > 0) & ~[nodes.is_sleeping];
    CH_ids = find(CH_alive);
    n_ch = length(CH_ids);
    
    if n_ch > 0
            CH_pos = reshape([nodes(CH_ids).pos], 2, [])';
    else
            CH_pos = zeros(0,2);
    end
    
    max_hops = max(8, ceil(2*sqrt(2)*x / ch_range) + 2);
    
    for i = 1:total_nodes
            nodes(i).buffer = 0;
            nodes(i).buf_time_sum = 0;
    end
    
    c_pkts = 0;
    c_loss = 0;
    c_d_sum = 0;
    c_d_cnt = 0;
    c_sink_rx = 0; % FIX 4: accumulated sink Rx energy this cycle [J]
    
    %==================================================================%
    % TIME-STEP LOOP %
    %==================================================================%
    for ti = 1:n_steps
            t = t_vec(ti);
            
            sn1 = sink_pos_at(t, s1_path, sink_speed, sojourn_time);
            sn2 = sink_pos_at(t, s2_path, sink_speed, sojourn_time);
            
            % Energy-based sleep / wake transitions
            for i = 1:total_nodes
                if nodes(i).energy <= 0
                                nodes(i).is_sleeping = true;
                elseif ~nodes(i).is_CH && nodes(i).energy < E_sleep_th
                                nodes(i).is_sleeping = true;
                elseif nodes(i).energy >= E_sleep_th * 1.5
                                nodes(i).is_sleeping = false;
                end
            end
            
            % Buffered data model – members generate at interval T_data
            if mod(round(t), T_data) == 0
                for i = 1:total_nodes
                    if nodes(i).energy <= 0 || nodes(i).is_sleeping || nodes(i).is_CH
                        continue;
                    end
                    ch = nodes(i).primary_CH_id;
                    if ch <= 0 || nodes(ch).energy <= 0, continue; end
                    d_mc = norm(nodes(i).pos - nodes(ch).pos);
                    if d_mc > node_range, continue; end
                    
                    E_m = E_elec*packet_size_bits + E_fs*packet_size_bits*d_mc^2;
                    nodes(i).energy = max(0, nodes(i).energy - E_m);
                    nodes(ch).energy = max(0, nodes(ch).energy - E_rx*packet_size_bits);
                    
                    if nodes(ch).buffer < buf_max
                                        nodes(ch).buffer = nodes(ch).buffer + 1;
                                        nodes(ch).buf_time_sum = nodes(ch).buf_time_sum + t;
                    end
                end
            end
            
            % CH transmits whenever sink is within range
            tx_list = [];
            for ki = 1:n_ch
                        ch = CH_ids(ki);
                if nodes(ch).energy <= 0 || nodes(ch).is_sleeping || nodes(ch).buffer == 0
                    continue;
                end
                if norm(nodes(ch).pos - sn1) <= ch_range || ...
                               norm(nodes(ch).pos - sn2) <= ch_range
                                tx_list(end+1) = ki; %#ok<AGROW>
                end
            end
            
            if isempty(tx_list), continue; end
            
            % -----------------------------------------------------------------
            % UPGRADE 1: Probabilistic Link Model
            for ti_ch = 1:length(tx_list)
                        ki = tx_list(ti_ch);
                        ch_id = CH_ids(ki);
                if nodes(ch_id).energy <= 0, continue; end
                
                % Choose closer sink
                d1 = norm(nodes(ch_id).pos - sn1);
                d2 = norm(nodes(ch_id).pos - sn2);
                if d1 <= d2; snk_pos = sn1; d_snk = d1;
                else; snk_pos = sn2; d_snk = d2;
                end
                
                % Signal power from this CH at the chosen sink
                d_sig = max(1, d_snk);
                alpha = 2 + 2*(d_sig >= E_th); % 2 (LOS) or 4 (NLOS)
                P_sig = P_tx_ref / d_sig^alpha;
                
                % Aggregate interference from all other concurrent transmitters
                P_int = 0;
                for tj_ch = 1:length(tx_list)
                    if tj_ch == ti_ch, continue; end
                                    cj = CH_ids(tx_list(tj_ch));
                                    dij = max(1, norm(nodes(cj).pos - snk_pos));
                                    a = 2 + 2*(dij >= E_th);
                                    P_int = P_int + P_tx_ref / dij^a;
                end
                
                % Compute SINR
                SINR = P_sig / (N0 + P_int);
                SINR = max(SINR, 1e-12); % FIX 6: realistic noise floor
                SINR_dB = 10 * log10(SINR);
                
                % Logistic success probability
                P_success = 1 / (1 + exp(-k_logistic * (SINR_dB - SINR_th_dB)));
                
                % Dijkstra multi-hop routing cost
                [E_tx, nodes] = dijkstra_tx(ch_id, snk_pos, nodes, ...
                    CH_ids, CH_pos, n_ch, ch_range, max_hops, ...
                    E_elec, E_rx, E_fs, E_mp, E_th, packet_size_bits, E_init);
                E_agg = E_DA * max(1, nodes(ch_id).buffer) * packet_size_bits;
                
                % ---- Probabilistic outcome decision ----
                if rand() < P_success
                    % SUCCESS
                    nodes(ch_id).energy = max(0, nodes(ch_id).energy - E_tx - E_agg);
                    n_buf = nodes(ch_id).buffer;
                    c_sink_rx = c_sink_rx + E_sink_rx * packet_size_bits * max(1, n_buf);
                    
                    if n_buf > 0
                                        avg_gen = nodes(ch_id).buf_time_sum / n_buf;
                                        c_d_sum = c_d_sum + (t - avg_gen) * n_buf;
                                        c_d_cnt = c_d_cnt + n_buf;
                    end
                    c_pkts = c_pkts + n_buf;
                    nodes(ch_id).buffer = 0;
                    nodes(ch_id).buf_time_sum = 0;
                else
                    % FAILURE
                    nodes(ch_id).energy = max(0, nodes(ch_id).energy - 0.05*E_tx);
                    drop_ratio = 1 - P_success;
                    n_buf_old = nodes(ch_id).buffer;
                    n_buf_new = round(n_buf_old * (1 - drop_ratio));
                    n_dropped = n_buf_old - n_buf_new;
                    c_loss = c_loss + n_dropped;
                    nodes(ch_id).buffer = n_buf_new;
                    
                    if n_buf_new == 0
                                        nodes(ch_id).buf_time_sum = 0;
                    elseif n_buf_old > 0
                                        nodes(ch_id).buf_time_sum = nodes(ch_id).buf_time_sum ...
                                            * (n_buf_new / n_buf_old);
                    end
                end
            end % transmitting CH loop
    end % time-step loop
    
    % ---- Record metrics ----
    alive_log(cycle) = sum([nodes.energy] > 0);
    energy_log(:,cycle) = [nodes.energy]';
    pkt_log(cycle) = c_pkts;
    pkt_loss_log(cycle) = c_loss;
    sink_rx_e_log(cycle) = c_sink_rx;
    
    if c_d_cnt > 0
            delay_log(cycle) = c_d_sum / c_d_cnt;
    end
    
    fprintf('Cycle %2d | Alive: %4d | Mean E: %.4f J | Pkts: %4d | Lost: %4d | SinkRx: %.4f J | Delay: %5.1fs\n', ...
        cycle, alive_log(cycle), mean([nodes.energy]), c_pkts, c_loss, c_sink_rx, delay_log(cycle));
end

%% ============================= PLOTS ==============================
figure('Name','WSN Upgraded Model – v4 Refined','Position',[80 80 960 1350]);

subplot(6,1,1);
plot(1:num_cycles, alive_log,'b-o','LineWidth',2,'MarkerSize',4);
xlabel('Cycle'); ylabel('Alive Nodes');
title('Network Lifetime – Alive Nodes per Cycle'); grid on; xlim([1 num_cycles]);

subplot(6,1,2);
plot(1:num_cycles, mean(energy_log,1),'r-s','LineWidth',2,'MarkerSize',4);
xlabel('Cycle'); ylabel('Mean Energy [J]');
title('Average Node Energy per Cycle'); grid on; xlim([1 num_cycles]);

subplot(6,1,3);
bar(1:num_cycles, pkt_log,'FaceColor',[0.18 0.63 0.33]);
xlabel('Cycle'); ylabel('Packets Delivered');
title('Successful Packet Delivery per Cycle (Probabilistic Link, k=0.3)'); grid on;

subplot(6,1,4);
bar(1:num_cycles, pkt_loss_log,'FaceColor',[0.85 0.20 0.20]);
xlabel('Cycle'); ylabel('Packets Lost');
title('Packet Loss per Cycle (P_{success}-Proportional Drop)'); grid on;

subplot(6,1,5);
bar(1:num_cycles, sink_rx_e_log * 1e6, 'FaceColor', [0.20 0.45 0.80]);
xlabel('Cycle'); ylabel('Sink Rx Energy [\muJ]');
title('Sink Reception Energy per Cycle [FIX 4]'); grid on;

subplot(6,1,6);
plot(1:num_cycles, delay_log,'m-d','LineWidth',2,'MarkerSize',4);
xlabel('Cycle'); ylabel('Avg Delay [s]');
title('End-to-End Packet Delay per Cycle (Buffer Model)'); grid on; xlim([1 num_cycles]);

fnd = find(alive_log < total_nodes, 1);
fprintf('\n=== SIMULATION COMPLETE ===\n');
fprintf('FND (first node death) : cycle %d\n', fnd);
fprintf('Alive at end : %d / %d\n', alive_log(end), total_nodes);
fprintf('Total packets delivered: %d\n', sum(pkt_log));
fprintf('Total packets lost : %d\n', sum(pkt_loss_log));
fprintf('Packet delivery ratio : %.2f%%\n', ...
    100 * sum(pkt_log) / max(1, sum(pkt_log) + sum(pkt_loss_log)));
fprintf('Total sink Rx energy : %.4f J\n', sum(sink_rx_e_log));
fprintf('Average packet delay : %.1f s\n', mean(delay_log(delay_log>0)));

%% ========================= LOCAL FUNCTIONS ==========================
function T = cycle_dur(path, spd, soj)
    n = size(path,1);
    T = n * soj;
    for s = 2:n
            T = T + norm(path(s,:) - path(s-1,:)) / spd;
    end
end

function pos = sink_pos_at(t, path, spd, soj)
    n = size(path,1);
    t_cur = 0;
    for s = 1:n
            t_soj = t_cur + soj;
            if t <= t_soj
                        pos = path(s,:); return;
            end
            t_cur = t_soj;
            
            if s < n
                        d_seg = norm(path(s+1,:) - path(s,:));
                        t_travel = d_seg / spd;
                        t_end = t_cur + t_travel;
                if t <= t_end
                                frac = (t - t_cur) / t_travel;
                                pos = path(s,:) + frac * (path(s+1,:) - path(s,:));
                                return;
                end
                        t_cur = t_end;
            end
    end
    pos = path(end,:);
end

function nodes = adaptive_cluster_kmeans(nodes, total_nodes, K, rp_pos, num_iter)
    centers = rp_pos;
    for iter = 1:num_iter
        for i = 1:total_nodes
            if nodes(i).energy <= 0, continue; end
                        d = sqrt(sum((centers - nodes(i).pos).^2, 2));
                        [~, nodes(i).cluster_id] = min(d);
        end
        new_centers = centers;
        for c = 1:K
                    mem = find(([nodes.cluster_id] == c) & ([nodes.energy] > 0));
            if ~isempty(mem)
                                pm = reshape([nodes(mem).pos], 2, [])';
                                new_centers(c,:) = mean(pm, 1);
            end
        end
        centers = new_centers;
    end
    
    for i = 1:total_nodes
        if nodes(i).energy <= 0, continue; end
            d = sqrt(sum((centers - nodes(i).pos).^2, 2));
            [~, nodes(i).cluster_id] = min(d);
    end
end

function nodes = elect_CH_adaptive(nodes, num_rp, E_init, node_range, mean_e_ratio)
    if mean_e_ratio > 0.7; w=[0.40, 0.40, 0.20];
    elseif mean_e_ratio > 0.3; w=[0.60, 0.30, 0.10];
    else; w=[0.80, 0.15, 0.05];
    end
    
    N = length(nodes);
    for i = 1:N
            nodes(i).is_CH = false;
            nodes(i).is_backup = false;
            nodes(i).primary_CH_id = 0;
            nodes(i).backup_CH_id = 0;
    end
    
    cids = unique([nodes.cluster_id]);
    cids = cids(cids > 0);
    
    for c = cids
            alive = find(([nodes.cluster_id]==c) & ([nodes.energy]>0) & ~[nodes.is_sleeping]);
            if length(alive) < 2, continue; end
            
            pos_c = reshape([nodes(alive).pos], 2, [])';
            center = mean(pos_c, 1);
            scores = zeros(length(alive),1);
            
            for j = 1:length(alive)
                        nd = nodes(alive(j));
                        d_c = norm(nd.pos - center);
                        d_nb = sqrt(sum((pos_c - nd.pos).^2, 2));
                        n_nb = (sum(d_nb <= node_range) - 1) / max(1, length(alive)-1);
                        
                switch nd.type
                    case 'super'; type_mult = 1.20;
                    case 'advanced'; type_mult = 1.10;
                    otherwise; type_mult = 1.00;
                end
                        e_score = nd.energy / nd.E_node_init;
                        scores(j) = type_mult * ...
                            (w(1)*e_score + w(2)/(d_c+1) + w(3)*n_nb);
            end
            
            [~,i1] = max(scores); prim = alive(i1);
            scores(i1) = -inf;
            [~,i2] = max(scores); back = alive(i2);
            
            nodes(prim).is_CH = true;
            nodes(back).is_backup = true;
            
            for j = 1:length(alive)
                        nid = alive(j);
                if nid ~= prim && nid ~= back
                                    nodes(nid).primary_CH_id = prim;
                                    nodes(nid).backup_CH_id = back;
                end
            end
    end
end

function nodes = promote_backup(nodes)
    for i = 1:length(nodes)
            pch = nodes(i).primary_CH_id;
        if pch > 0 && nodes(pch).energy <= 0
                    bch = nodes(i).backup_CH_id;
            if bch > 0 && nodes(bch).energy > 0
                                nodes(i).primary_CH_id = bch;
                                nodes(bch).is_CH = true;
                                nodes(bch).is_backup = false;
            else
                                nodes(i).primary_CH_id = 0;
            end
        end
    end
end

function [E_total, nodes] = dijkstra_tx(src_id, dst_pos, nodes, CH_ids, CH_pos, ...
        n_ch, ch_range, max_hops, E_elec, E_rx, E_fs, E_mp, E_th, pkt_bits, E_init)
        
    src_k = find(CH_ids == src_id, 1);
    if isempty(src_k) || n_ch == 0
            E_total = tx_e(nodes(src_id).pos, dst_pos, E_elec, E_fs, E_mp, E_th, pkt_bits);
            return;
    end
    
    SINK = n_ch + 1;
    INF = 1e18;
    dist = INF * ones(1, SINK);
    prev = zeros(1, SINK);
    vis = false(1, SINK);
    dist(src_k) = 0;
    
    for iter = 1:SINK
            tmp = dist; tmp(vis) = INF;
            [d_u, u] = min(tmp);
            if d_u >= INF, break; end
            vis(u) = true;
            if u == SINK, break; end
            
            u_pos = CH_pos(u,:);
            for v = 1:n_ch
                if vis(v), continue; end
                            cj = CH_ids(v);
                if nodes(cj).energy <= 0, continue; end
                            d_uv = norm(u_pos - CH_pos(v,:));
                if d_uv > ch_range, continue; end
                
                            e_fac = E_init / max(1e-9, nodes(cj).energy);
                            w_uv = d_uv * e_fac;
                if dist(u) + w_uv < dist(v)
                                    dist(v) = dist(u) + w_uv;
                                    prev(v) = u;
                end
            end
            
            d_snk = norm(u_pos - dst_pos);
            if d_snk <= ch_range
                if dist(u) + d_snk < dist(SINK)
                                    dist(SINK) = dist(u) + d_snk;
                                    prev(SINK) = u;
                end
            end
    end
    
    path_k = [];
    cur = SINK;
    hops = 0;
    while cur ~= 0 && cur ~= src_k && hops <= max_hops
            path_k(end+1) = cur;
            cur = prev(cur);
            hops = hops + 1;
    end
    path_k(end+1) = src_k;
    path_k = fliplr(path_k);
    
    if length(path_k) < 2
            E_total = tx_e(nodes(src_id).pos, dst_pos, E_elec, E_fs, E_mp, E_th, pkt_bits);
            return;
    end
    
    E_total = 0;
    for p = 1:length(path_k)-1
            uk = path_k(p);
            vk = path_k(p+1);
            if uk > n_ch, break; end
            
            u_pos = CH_pos(uk,:);
            if vk == SINK
                        v_pos = dst_pos;
            elseif vk <= n_ch
                        v_pos = CH_pos(vk,:);
                        nodes(CH_ids(vk)).energy = max(0, nodes(CH_ids(vk)).energy - E_rx*pkt_bits);
            else
                break;
            end
            E_total = E_total + tx_e(u_pos, v_pos, E_elec, E_fs, E_mp, E_th, pkt_bits);
    end
    
    if E_total == 0
            E_total = tx_e(nodes(src_id).pos, dst_pos, E_elec, E_fs, E_mp, E_th, pkt_bits);
    end
end

function E = tx_e(src, dst, E_elec, E_fs, E_mp, E_th, bits)
    d = norm(src - dst);
    if d < E_th
            E = E_elec*bits + E_fs*bits*d^2;
    else
            E = E_elec*bits + E_mp*bits*d^4;
    end
end