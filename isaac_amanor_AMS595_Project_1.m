%% PROJECT_MONTE_CARLO
%% Task 1 (minimal): Monte Carlo estimate of pi using a FOR loop over fixed N.    % Brief description of the task
% Outputs: (1) pi_hat vs N, (2) |pi_hat - pi| vs N, (3) time vs N, (4) |error| vs time. % Required plots for the report

%% Housekeeping
clear; clc; close all;% Clear variables, clear command window, and close all open figures
rng(0); % Fix the random-number seed so results are reproducible
% LaTeX text rendering (for nicer math in labels/titles)
set(groot,'defaultTextInterpreter','latex');            % Use LaTeX for text (titles, axis labels)
set(groot,'defaultAxesTickLabelInterpreter','latex');   % Use LaTeX for tick labels
set(groot,'defaultLegendInterpreter','latex');          % Use LaTeX for legend entries

%% Configuration
N_values = round(logspace(3,9,10));% Vector of sample sizes N: 10^3, ..., 10^9 (10 log-spaced points)
true_pi  = pi;  % MATLAB's built-in value of pi (used only for computing error)
tiny     = 1e-12; % Small positive number to avoid log(0) on log-scale plots
maxBatch = 5e6; % Upper bound on one-shot random draws (keeps memory usage in check)

%% Preallocate results (preallocation avoids reallocations in the loop and speeds things up)
numN      = numel(N_values); % Number of different N values we will test
pi_hat    = zeros(numN,1);% Placeholder for Monte Carlo estimates of pi at each N
abs_error = zeros(numN,1);% Placeholder for absolute error |pi_hat - true_pi| at each N
run_time  = zeros(numN,1);% Placeholder for elapsed time (seconds) at each N

%% FOR-loop over fixed N (project requirement)
for k = 1:numN % Loop over the index of N values
    N = N_values(k); % Pick the current sample size N

    % Compute estimate using a single vectorized draw if small enough, else stream in batches
    tStart = tic; % Start a stopwatch timer for this iteration
    if N <= maxBatch % If N is small enough, draw all points at once
        x = rand(N,1);  % Draw N x-coordinates uniformly in [0,1]
        y = rand(N,1);% Draw N y-coordinates uniformly in [0,1]
        count_inside = sum((x.^2 + y.^2) <= 1); % Count points that land in the unit quarter-circle (x^2+y^2 <= 1)
    else % Otherwise, stream to limit memory
        remaining = N;% How many points are left to draw
        count_inside = 0; % Initialize inside-counter
        while remaining > 0 % Keep drawing until we reach N samples
            b = min(maxBatch, remaining);  % Batch size: the smaller of maxBatch or what remains
            xb = rand(b,1); % Draw a batch of x-coordinates
            yb = rand(b,1); % Draw a batch of y-coordinates
            count_inside = count_inside + sum((xb.^2 + yb.^2) <= 1); % Accumulate inside-count for this batch
            remaining = remaining - b; % Decrease the remaining sample count
        end % End streaming loop
    end % End conditional (vectorized vs streaming)

    pi_hat(k) = 4 * (count_inside / N); % Monte Carlo estimate of pi = 4 * (fraction of points inside)
    run_time(k) = toc(tStart);% Record elapsed time for this N
    abs_error(k) = abs(pi_hat(k) - true_pi);% Compute absolute error against the reference pi (for evaluation only)
end % End FOR-loop over N

%% Log-scale safety (avoid warnings or -Inf due to zeros on log axes)
abs_error(abs_error <= 0) = tiny; % Replace zeros (if any) by tiny positive number for log plotting
run_time(run_time <= 0)   = tiny; % Replace zeros (if any) by tiny positive number for log plotting

%% Plot 1: pi_hat vs N (with true pi reference)
figure('Color','w'); % Create a new figure with white background (better for export)
plot(N_values, pi_hat, 'o-', 'LineWidth', 1.2); hold on; % Plot estimated pi against N on linear y / log x
h = yline(true_pi, '--', 'LineWidth', 1.2); % Add a horizontal reference line at true pi
h.Label = 'True $\pi$'; % Label the reference line using LaTeX
h.Interpreter = 'latex';% Ensure the label uses LaTeX interpreter
hold off; % Release the plot hold
set(gca,'XScale','log'); grid on; % Use log scale on x-axis; turn on the grid
xlabel('Number of points $N$ (log)'); % x-axis label
ylabel('Estimated $\pi$'); % y-axis label
title('Estimated $\pi$ vs $N$'); % Title for the figure
saveas(gcf,'fig_task1_min_pi_vs_N.png');    % Save a PNG (handy for quick viewing)
print(gcf,'-dpdf','fig_task1_min_pi_vs_N.pdf'); % Save a vector PDF (best for the report)

%% Plot 2: Absolute Error vs N (log-log)
figure('Color','w'); % New figure with white background
loglog(N_values, abs_error, 's-', 'LineWidth', 1.2); grid on; % Log-log plot of absolute error vs N
xlabel('Number of points $N$ (log)');       % x-axis label
ylabel('Absolute error $|\hat{\pi}-\pi|$ (log)'); % y-axis label with LaTeX math
title('Absolute error vs $N$');             % Title for the figure
saveas(gcf,'fig_task1_min_error_vs_N.png'); % Save a PNG copy
print(gcf,'-dpdf','fig_task1_min_error_vs_N.pdf'); % Save a vector PDF copy

%% Plot 3: Run time vs N (log-log)
figure('Color','w'); % New figure with white background
loglog(N_values, run_time, 'd-', 'LineWidth', 1.2); grid on; % Log-log plot of runtime vs N
xlabel('Execution time (s, log)');% x-axis label
ylabel('Execution time (s, log)');% y-axis label
title('Execution time vs $N$'); % Title for the figure
saveas(gcf,'fig_task1_min_time_vs_N.png');  % Save a PNG copy
print(gcf,'-dpdf','fig_task1_min_time_vs_N.pdf'); % Save a vector PDF copy

%% Plot 4: Absolute Error vs Run Time (precision–cost)
figure('Color','w'); % New figure with white background
loglog(run_time, abs_error, 'o-', 'LineWidth', 1.2); grid on; % Log-log plot of error vs time (precision vs cost)
xlabel('Execution time (s, log)'); % x-axis label
ylabel('Absolute error $|\hat{\pi}-\pi|$ (log)'); % y-axis label with LaTeX math
title('Precision vs computational cost');   % Title for the figure
saveas(gcf,'fig_task1_absolute_error_vs_time.png'); % Save a PNG copy
print(gcf,'-dpdf','fig_task1_absolute_error_vs_time.pdf'); % Save a vector PDF copy

%% Save results table for report (valid MATLAB variable names)
results_tbl = table( ...   % Construct a table to store results
    N_values(:), ...   % Column 1: the N values (as a column vector)
    pi_hat(:), ...  % Column 2: Monte Carlo estimates of pi
    abs_error(:), ...  % Column 3: absolute errors
    run_time(:), ...  % Column 4: execution times (seconds)
    'VariableNames', {'N','pi_hat','abs_error','time_sec'}); % Use MATLAB-safe variable names for CSV

% Optional human-readable descriptions (metadata; helpful when inspecting the table)
results_tbl.Properties.VariableDescriptions = { ... % Add descriptions for each column
    'Sample size N', ... % Description for column N
    'Monte Carlo estimate \hat{\pi}', ...   % Description for column pi_hat
    'Absolute error |\hat{\pi}-\pi|', ...   % Description for column abs_error
    'Execution time (seconds)'}; % Description for column time_sec

writetable(results_tbl, 'isaac_amanor_task1_Forloop.csv'); % Write the table to CSV for your report / repository

%% Console summary
disp(' Task 1 (minimal) Summary'); % Header line for console
disp(results_tbl);  % Display the results table in the command window
disp('Saved figures (PNG + PDF): fig_task1_min_pi_vs_N, fig_task1_min_error_vs_N, fig_task1_min_time_vs_N, fig_task1_absolute_error_vs_time'); % List of saved figures
disp('Saved table: isaac_amanor_task1_Forloop.csv'); % Confirm the saved CSV file


%% project1_task2_while_precision_basic_eol.m%Section header & filename hint
% Task 2 — Monte Carlo estimate of pi with a WHILE loop and auto-stop (no true pi).%High-level description
% Stop rule: successive running estimates must agree to s significant figures for STABILITY_OK checks.%Stop condition

clear; clc; close all;%Reset workspace, console, and figures
rng(1,'twister');%Reproducible RNG seed

%% Configuration (no plotting)%Config section
s_list       = [2 3 4 5 6 7];%Target significant figures to test
BATCH        = 5e6;%Samples per while-iteration
CHECK_EVERY  = 1e6;%Check stability every this many new samples
STABILITY_OK = 5;%Consecutive matches required to stop
MAX_POINTS   = 1e9;%Safety cap on total samples
save_out     = true;%Enable CSV save

out_dir = "task2_outputs";%Output directory
if save_out && ~exist(out_dir,'dir'), mkdir(out_dir); end%Create directory if missing
out_csv = fullfile(out_dir, "task2_results_" + datestr(now,'yyyymmdd_HHMMSS') + ".csv");%Timestamped CSV filename

results = table('Size',[numel(s_list) 6], ...%Preallocate table
    'VariableTypes',{'double','double','double','double','double','double'}, ...%Column types
    'VariableNames',{'s_sigfigs','N_used','pi_hat','time_s','iters_batches','checks'});%Column names

fprintf('Task 2: while-loop with stability-based stopping (no true pi used)\n');%Banner

%% Main loop over requested precisions%Loop header
for idx = 1:numel(s_list)%For each target s
    s = s_list(idx);%Current target significant figures
    fprintf('\nTarget: %d significant figures\n', s);%Announce target

    n_inside   = 0;%Running count of hits inside quarter-circle
    N          = 0;%Total samples seen
    t0         = tic;%Start timer
    last_pi    = NaN;%Last estimate at a check
    consec_ok  = 0;%Consecutive agreements counter
    next_check = CHECK_EVERY;%Next absolute N to check
    iters      = 0;%Batches processed
    chk_count  = 0;%Checks performed

    while true%Stream until stable or capped
        b = min(BATCH, MAX_POINTS - N);%Batch size obeying cap
        if b <= 0%Cap reached
            warning('Reached MAX_POINTS without stability for s=%d.', s);%Warn user
            break;%Exit while
        end%End cap guard
        iters = iters + 1;%Count this batch

        x = rand(b,1);%Uniform x in [0,1]
        y = rand(b,1);%Uniform y in [0,1]
        n_inside = n_inside + sum(x.*x + y.*y <= 1.0);%Accumulate hits
        N = N + b;%Update total samples

        while N >= next_check%Perform due checks (may be multiple)
            chk_count = chk_count + 1;%Count this check
            p_hat  = n_inside / N;%Hit probability estimate
            pi_now = 4 * p_hat;%Running pi estimate

            if ~isnan(last_pi)%If there is a previous estimate
                if agree_to_s_sigfigs(pi_now, last_pi, s)%Check s-sig-fig agreement
                    consec_ok = consec_ok + 1;%Extend streak
                else%Mismatch
                    consec_ok = 0;%Reset streak
                end%End agreement test
            end%End previous-estimate guard

            last_pi   = pi_now;%Cache for next comparison
            next_check = next_check + CHECK_EVERY;%Schedule next check

            if consec_ok >= STABILITY_OK%Stable enough
                break;%Exit inner check loop
            end%End stability condition
        end%End while over due checks

        if consec_ok >= STABILITY_OK%If stable after this batch
            break;%Exit main while
        end%End post-batch stability check
    end%End streaming while

    elapsed    = toc(t0);%Elapsed time for this s
    pi_hat_val = 4 * (n_inside / N);%Final estimate from counts

    results.s_sigfigs(idx)     = s;%Log s
    results.N_used(idx)        = N;%Log N used
    results.pi_hat(idx)        = pi_hat_val;%Log pi_hat
    results.time_s(idx)        = elapsed;%Log time
    results.iters_batches(idx) = iters;%Log batches
    results.checks(idx)        = chk_count;%Log checks

    fprintf('  Done: N=%s, pi_hat=%.10f, time=%.2fs, batches=%d, checks=%d\n', ...%One-line summary
        format_int(N), pi_hat_val, elapsed, iters, chk_count);%Pretty numbers
end%End for over s_list

%% Save CSV (new file per run)%CSV section
if save_out%If saving enabled
    writetable(results, out_csv);%Write table to CSV
    fprintf('\nSaved: %s\n', out_csv);%Report path
end%End save guard

%% Summary%Console summary
disp('Summary (Task 2):');%Label
disp(results);%Show results


%% Helpers (single, shared) 
function tf = agree_to_s_sigfigs(a,b,s)%Return true if a and b agree to s significant figures
    if ~isfinite(a) || ~isfinite(b), tf = false; return; end%Reject NaN/Inf
    if a == 0 || b == 0, tf = (a == b); return; end%Zero-edge case
    ea = floor(log10(abs(a)));%Order of magnitude of a
    eb = floor(log10(abs(b)));%Order of magnitude of b
    if ea ~= eb, tf = false; return; end%Different magnitudes cannot match to s sig figs
    tol = 0.5 * 10^(ea - (s-1));%Half-unit at the s-th significant digit (absolute tolerance)
    tf  = abs(a - b) <= tol;%Within tolerance => agree
end%End agree_to_s_sigfigs

function t = format_int(N)%Add thousands separators
    t = regexprep(sprintf('%.0f',N), '(\d)(?=(\d{3})+(?!\d))', '$1,');%Comma grouping
end%End format_int

function s = sigfig_string(x, n) % Format numeric x to n significant figures as a string (fixed point where possible)
    if x==0, s = sprintf('%.*f', n-1, 0); return; end % Special case: zero formatted with n-1 decimals
    e = floor(log10(abs(x))); % Order of magnitude of x
    f = max(n - 1 - e, 0); % Number of decimals to show to reach n significant figures
    s = sprintf(['%.' num2str(f) 'f'], x); % Build formatted fixed-point string
end % End sigfig_string


%% Task 3 function 
function Pi_estimate_n = mcpi_task3_live(precision_sf) % Function entry point; returns Monte Carlo estimate of pi
% mcpi_task3_live  Live Monte Carlo estimator of pi with sig-fig stability stop (CPU-only). % One-line description
% USAGE: Pi_estimate_n = mcpi_task3_live(6);  % target 6 significant figures % Example usage for users
% Stop rule (no true pi): running estimates must agree to s sig figs for K consecutive checks (default K=4). % Stopping criterion summary

%  Inputs & tunables % Section: input parsing and knobs
if nargin==0 || isempty(precision_sf), precision_sf = 7; end % Default to 7 significant figures if input omitted/empty
validateattributes(precision_sf, {'numeric'}, {'scalar','integer','>=',1}); % Validate input: positive integer scalar

BATCH        = 2e5; % Samples per streamed update (tune for machine)
CHECK_EVERY  = 1e5; % Perform stability checks every this many new samples
STABILITY_OK = 4; % Require this many consecutive matches at s significant figures
MAX_POINTS   = 1e9;% Hard cap to avoid runaway runs
SHOW_LAST    = 3e4;% Render only the most recent SHOW_LAST points (keeps UI responsive)
rng(0,'twister');% Reproducible random stream

%  Figure & visuals % Section: plotting setup
fig = figure('Color','w','Name',sprintf('Monte Carlo pi (%d-sf)',precision_sf)); % Create figure window with plain-text name
axis square; hold on; box on; axis([0 1 0 1]); % Square axes, keep [0,1]x[0,1] viewport, enable box
xlabel('x'); ylabel('y'); % Axis labels
title(sprintf('Estimating $\\pi$ to %d significant figures (stability stop)', precision_sf), 'Interpreter','latex'); % LaTeX title with \pi in math mode

tt  = linspace(0,pi/2,400); % Parameter for quarter-circle arc (0 to pi/2)
hB  = plot(cos(tt), sin(tt), 'k-','LineWidth',1.2); % Draw quarter-circle boundary

hIn  = scatter(nan,nan,6,'filled','MarkerFaceColor',[0.20 0.60 1.00],'MarkerEdgeColor','none'); % Scatter for points inside circle
hOut = scatter(nan,nan,6,'filled','MarkerFaceColor',[0.85 0.20 0.20],'MarkerEdgeColor','none'); % Scatter for points outside circle
legend([hB hIn hOut], {'Quarter circle','Inside $(x^2+y^2\le 1)$','Outside'}, 'Interpreter','latex','Location','southoutside'); % Legend with LaTeX

statusTxt = text(0.02,0.98,'','Units','normalized','VerticalAlignment','top','FontSize',11,'Interpreter','latex'); % Status text in top-left

%  Circular ring buffer (preallocated; no growth) % Section: rolling buffer for fast UI
Xbuf = nan(SHOW_LAST,1); % X buffer (NaNs initialize as empty)
Ybuf = nan(SHOW_LAST,1); % Y buffer (NaNs initialize as empty)
Cbuf = false(SHOW_LAST,1); % Class buffer (true=inside, false=outside)
head = 0; % Write index (1..SHOW_LAST rolling)

%  Streaming state % Section: counters and stopping state
n_inside       = 0; % Total hits inside quarter circle
N              = 0; % Total samples processed
next_check     = CHECK_EVERY; % Absolute sample index at which to perform next check
consec_matches = 0; % Consecutive sig-fig agreement counter
last_pi        = NaN; % Previous checked estimate (NaN means none yet)
Pi_estimate_n  = NaN; % Return value placeholder

% - Main loop % Section: generate, update, check, render
while true % Loop until stability achieved, cap reached, or user closes figure
    b = min(BATCH, MAX_POINTS - N); % Choose batch size without exceeding MAX_POINTS
    if b <= 0 % If no room left to sample more points
        warning('Reached MAX_POINTS without meeting the stability rule.'); % Notify that cap was hit
        break; % Exit loop due to cap
    end % End cap check

    x = rand(b,1);  y = rand(b,1); % Draw b uniform points in [0,1]^2
    inside = (x.*x + y.*y) <= 1.0; % Logical mask of points inside quarter circle
    n_inside = n_inside + sum(inside); % Accumulate count of hits
    N = N + b; % Update total samples

    kvis = min(b, SHOW_LAST); % Determine how many of the latest points to show
    if kvis > 0 % If there are points to write into the ring buffer
        xs = x(end-kvis+1:end); % Latest kvis x-values
        ys = y(end-kvis+1:end); % Latest kvis y-values
        cs = inside(end-kvis+1:end); % Latest kvis class flags
        idxWrite = mod((head + (1:kvis)) - 1, SHOW_LAST) + 1; % Compute write positions (wrap around)
        Xbuf(idxWrite) = xs; % Write x slice into buffer
        Ybuf(idxWrite) = ys; % Write y slice into buffer
        Cbuf(idxWrite) = cs; % Write class slice into buffer
        head = idxWrite(end); % Advance head to last written index
    end % End ring-buffer write
    didCheck = false; % Flag to indicate at least one scheduled check ran
    while N >= next_check % Perform all checks whose thresholds were crossed by this batch
        didCheck = true; % Mark that we performed a check
        Pi_estimate_n = 4 * (n_inside / N); % Current running estimate of pi from hits ratio

        if head < SHOW_LAST % Build chronological order indices for plotting
            order = [head+1:SHOW_LAST, 1:head]; % If buffer not yet wrapped fully
        else
            order = 1:SHOW_LAST; % If buffer is exactly at end index
        end % End order construction
        maskIn  = Cbuf(order); % Inside mask in chronological order
        Xin  = Xbuf(order(maskIn)); % Chronological inside X
        Yin  = Ybuf(order(maskIn)); % Chronological inside Y
        Xout = Xbuf(order(~maskIn)); % Chronological outside X
        Yout = Ybuf(order(~maskIn)); % Chronological outside Y

        if isvalid(hIn) && isvalid(hOut) && isvalid(statusTxt) % Ensure graphics handles are still valid
            set(hIn,  'XData', Xin,  'YData', Yin); % Update inside scatter data
            set(hOut, 'XData', Xout, 'YData', Yout); % Update outside scatter data
            set(statusTxt,'String', sprintf('N = %s\n$\\hat{\\pi}$ = %.10f\nrounded (%d sf): %s', format_int(N), Pi_estimate_n, precision_sf, sigfig_string(Pi_estimate_n, precision_sf))); % Update status text with N, raw estimate, rounded
            drawnow limitrate; % Throttled UI refresh
        else % If figure was closed
            break; % Exit check loop early due to closed figure
        end % End handle-validity guard

        if ~isnan(last_pi) && agree_to_s_sigfigs(Pi_estimate_n, last_pi, precision_sf) % Test agreement with previous estimate to s sig figs
            consec_matches = consec_matches + 1; % Increase consecutive-agreement streak
        else % Else no agreement or first check
            consec_matches = 0; % Reset streak
        end % End agreement update
        last_pi = Pi_estimate_n; % Store current estimate for next comparison
        if consec_matches >= STABILITY_OK % If required consecutive agreements reached
            break; % Exit check loop (stability achieved)
        end % End stability check
        next_check = next_check + CHECK_EVERY; % Schedule next check at a later absolute N
    end % End while over scheduled checks
    if (didCheck && consec_matches >= STABILITY_OK) || ~ishandle(fig) % If stability achieved or figure closed
        break; % Exit main loop
    end % End exit condition
end % End main while loop
if ~isnan(Pi_estimate_n) % If we produced a valid estimate
    final_str = sigfig_string(Pi_estimate_n, precision_sf); % Build rounded string at requested sig figs
    text(0.02,0.06, sprintf('Final $\\hat{\\pi}$ = %s  (N = %s)', final_str, format_int(N)), 'Units','normalized','FontSize',11,'Interpreter','latex','BackgroundColor','w'); % Annotate final result on plot
    fprintf('Final pi estimate (to %d sig figs): %s  (N = %s)\n', precision_sf, final_str, format_int(N)); % Print final result to Command Window
else % If no estimate was produced (e.g., cap reached immediately)
    warning('No estimate produced.'); % Warn user that result is missing
end % End finalization block
end % End function mcpi_task3_live
