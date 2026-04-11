% DRIVER_SET2  Program 3 — Set 2 (General Trends) 

% This driver was heavily modeled after Dr. Gallivan's Program 1 driver. 
% I used ChatGPT to modify and understand what each line does.
%                                        
% uses my_lu, applies permutations, solves once, records metrics
%
% plots mean/max curves vs n and gamma histograms                               
                                                    

% ========================= USER PARAMETERS (EDIT THESE) ========================= 

rngseed    = 1;                              % random seed so results are repeatable
nmin       = 30;                             % smallest matrix size to test
nmax       = 250;                            % largest matrix size to test
ndelta     = 10;                             % step between sizes
ntests     = 5;                              % how many random problems per (n, class)

tol        = 1e-12;                          % pivot tolerance for your LU
modes_to_test = [0 1 2];                     % pivot modes to sweep: 0=none, 1=partial, 2=complete
classes_to_test = [1 2 3 4];                 % matrix classes: 1=LU, 2=Sym Pos Def, 3=Row Diag Dom, 4=Permuted Row Diag Dom

% Output style (match Program 1):
%   1 = print one line per test problem
%   2 = print one line per n (mean & max for that n)
%   3 = make plots (log10 of mean and max vs n) + gamma histograms

outputtype = 3;                              % choose 1, 2, or 3

% ------------------- choose a norm for all metrics -------------------

% Set normtype to 1, inf, or 'fro' for the 1-/∞-/Frobenius-norm everywhere.
% (2 is allowed but more expensive; assignment suggests finite norms in drivers.)
normtype   = 1;                              % use 1-norm across the board (fast and stable for these summaries)



% ========================= NO USER EDITS BELOW THIS LINE ========================  

rng(rngseed);                                % set RNG for reproducibility

nlist = (nmin:ndelta:nmax)';                 % column vector of sizes to test
nexperiments = numel(nlist);                 % number of different sizes
nModes = numel(modes_to_test);               % number of pivot modes we will test
nClasses = numel(classes_to_test);           % number of matrix classes we will test

% ---- Work array dimensions  ---------------------------------- 

Mrdim = nmax;                                % rows in work array M (≥ any n)
Mcdim = nmax;                                % cols in work array M
Acomprdim = nmax;                            % rows in scratch array Acomp
Acompcdim = nmax;                            % cols in scratch array Acomp
xydim = nmax;                                % length for x, b, xhat vectors

% ---- Overall maxima across all classes/modes/sizes  -------- 

overall_max_relF = 0.0;                      % max relative factorization error seen
overall_max_relX = 0.0;                      % max relative solution error seen
overall_max_relR = 0.0;                      % max relative residual seen
overall_max_gam  = 0.0;                      % max gamma seen

% ---- Header print  -------------------------------------------- 

fprintf(1,'\n**************************************** \n');                                          % visual break
fprintf(1,'Set 2 Driver — sizes n=%g:%g:%g, ntests=%g, tol=%g \n', nmin, ndelta, nmax, ntests, tol); % run setup line
fprintf(1,'Pivot modes tested: [%s]\n', num2str(modes_to_test));                                     % show modes
fprintf(1,'Matrix classes tested: [%s]\n', num2str(classes_to_test));                                % show classes
fprintf(1,'**************************************** \n');                                            % visual break

% ========================= LOOP OVER MATRIX CLASSES ============================ 

for ci = 1:nClasses                           % loop over classes
  ptype = classes_to_test(ci);                % current class code (1..4)

  % --- Per-class per-mode summary arrays (mean/max vs n) ----------------------- 

  mean_relF_modes = zeros(nexperiments, nModes);  % mean rel factor error per n, per mode
  max_relF_modes  = zeros(nexperiments, nModes);  % max  rel factor error per n, per mode
  mean_relX_modes = zeros(nexperiments, nModes);  % mean rel solution error per n, per mode
  max_relX_modes  = zeros(nexperiments, nModes);  % max  rel solution error per n, per mode
  mean_relR_modes = zeros(nexperiments, nModes);  % mean rel residual per n, per mode
  max_relR_modes  = zeros(nexperiments, nModes);  % max  rel residual per n, per mode
  mean_gam_modes  = zeros(nexperiments, nModes);  % mean gamma per n, per mode
  max_gam_modes   = zeros(nexperiments, nModes);  % max  gamma per n, per mode

  %collect all gamma by mode (per class) to make ONE overlaid histogram later

  GammaAll_modes  = cell(1, nModes);          % will hold NaN-filtered gamma across all n/tests for each mode

  fprintf(1,'\n===== Class %d =====\n', ptype);        % print class section header

  % ===================== LOOP OVER PIVOT MODES ================================= 

  for mi = 1:nModes                           % loop over pivot modes
    mode = modes_to_test(mi);                 % current pivot mode (0,1,2)

    fprintf(1,'-- Mode %d --\n', mode);       % label for this mode

    % --- Collect all gamma values across sizes/tests (for histogram) ---------- 

    GammaAll = [];                            % will append per-test gamma here

    experimentindex = 0;                      % index for position in mean/max arrays

    % ================= LOOP OVER SIZES n ====================================== 

    for n = nlist'                            % loop over n values (transpose to iterate values)
      experimentindex = experimentindex + 1;  % advance index for summaries

      % --- Allocate per-test arrays for this n -------------- 

      RelFactorErr = zeros(ntests,1);         % relative factorization error for each test at this n
      RelSolErr    = zeros(ntests,1);         % relative solution error for each test at this n
      RelResid     = zeros(ntests,1);         % relative residual for each test at this n
      GammaVal     = zeros(ntests,1);         % gamma for each test at this n

      if (outputtype == 1) || (outputtype == 2)         
        fprintf(1,'\nTests for n = %g (class %d, mode %d)\n', n, ptype, mode); % section header
      end

      % ================= LOOP OVER RANDOM TESTS =============================== 

      for itest = 1:ntests                    % repeat ntests times for randomness

        % --- Reinitialize work arrays  ----------------------- 

        M     = zeros(Mrdim, Mcdim);          % packed LU will be stored in top-left n×n block
        Acomp = zeros(Acomprdim, Acompcdim);  % scratch 
        xhat  = zeros(xydim,1);               % computed solution vector (top n entries used)
        b     = zeros(xydim,1);               % right-hand side vector (top n entries used)
        xtrue = zeros(xydim,1);               % true solution vector (top n entries used)

        % --- Build one test matrix A0 of size n×n for this class --------------  

        if ptype == 1                          % Class 1: A = L*U (random, full rank)
          L = tril(randn(n), -1);              % random strict-lower part
          L = 0.5 * L;                         % scale to keep off-diagonals modest
          L(1:n+1:end) = 1;                    % make L unit lower (diagonal ones)
          U = triu(randn(n));                  % random upper part
          d = 1 + rand(n,1);                   % diagonal entries in [1,2]
          U(1:n+1:end) = d;                    % set U diagonal
          A0 = L * U;                          % form A0 = L*U

        elseif ptype == 2                      % Class 2: SPD (A0 = Lt*Lt')
          Lt = tril(randn(n));                 % random lower-triangular
          Lt(1:n+1:end) = 1 + rand(n,1);       % positive diagonal
          A0 = Lt * Lt.';                      % symmetric positive definite

        elseif ptype == 3                      % Class 3: strictly row-diagonally dominant
          R = randn(n);                        % random dense matrix
          s = sum(abs(R), 2);                  % row-wise sums of absolute values
          A0 = R + diag(s + 1.0);              % boost diagonal to enforce dominance

        elseif ptype == 4                      % Class 4: permuted row-diagonally dominant
          R  = randn(n);                       % random dense matrix
          s  = sum(abs(R), 2);                 % row-wise sums
          A1 = R + diag(s + 1.0);              % make a diag-dominant matrix
          rp0 = randperm(n); cp0 = randperm(n);% random permutations
          A0 = A1(rp0, cp0);                   % apply permutations to obscure dominance

        else                                               % safety check for class code
          fprintf(1,'Unknown class ptype = %g\n', ptype);  % error message
          return;                                          % stop if bad input
        end

        % --- Pick a true x and compute b = A0*x (so we know the exact answer) -- 

        xtrue(1:n) = randn(n,1);               % random ground-truth solution
        b(1:n)     = A0 * xtrue(1:n);          % exact right-hand side

        

        % ----- store ||A|| separately and then floor with "at least 1" (NA) for denominators

        A0norm = norm(A0,      normtype);           % compute ||A|| in the chosen norm
        NA     = max(A0norm,   1);                  % denominator for factorization error (avoid divide by small)
        NB     = max(norm(b(1:n),  normtype), 1);   % denominator for residual (avoid divide by small)
        NX     = max(norm(xtrue(1:n), normtype), 1);% denominator for solution error (avoid divide by small)

        % --- FACTOR with your LU (packed) ------------------------------------- 

        % Use try/catch so a tiny/zero pivot doesn't stop the sweep;
        %      if an error is produced here, we mark this trial as failed (NaN metrics) and continue.

        lu_ok = true;                          % assume factorization succeeds
        lu_errmsg = '';                        % capture error string for printing outside catch
        try
          [M(1:n,1:n), rp, cp] = my_lu(A0, mode, tol);      % call: [M,rp,cp] = my_lu(A,mode,tol)
        catch ME
          lu_ok = false;                       % LU failed (e.g., near-zero pivot); mark as failed
          lu_errmsg = ME.message;              % save the error message for printing
        end
        if ~lu_ok
          RelFactorErr(itest) = NaN;           % mark all metrics NaN for this trial
          RelSolErr(itest)    = NaN;           % "
          RelResid(itest)     = NaN;           % "
          GammaVal(itest)     = NaN;           % "
          if outputtype == 1                   % optional notice for per-test output
            fprintf(1,'n=%4d  (LU failed: %s)\n', n, lu_errmsg); % use captured message
          end
          continue                             % go to next trial (do not attempt solve)
        end
        
        % --- Form Pr*A*Pc using your permutation routine ---------------------- % 

        Aperm = apply_perms(A0, rp, cp);       % Aperm = P_r * A0 * P_c

        % --- Rebuild L and U from packed M (Program 1 storage) ---------------- 

        Lhat = tril(M(1:n,1:n), -1) + eye(n);  % unit-lower L from strict-lower of M + I on diagonal
        Uhat = triu(M(1:n,1:n));               % upper-triangular U from upper of M

        % --- Factorization error: ||Pr*A*Pc - L*U|| / max(||A||,1) ------------ 

        Ahat = Lhat * Uhat;                                        % product of factors
        RelFactorErr(itest) = norm(Aperm - Ahat,  normtype) / NA;  % relative factorization error

        % --- Growth factor gamma: || |L||U| || / ||A|| ------------------------
        %         use growth_factor(M,A0,normtype); it returns || |L||U| || / ||A||

        gamma0 = growth_factor(M(1:n,1:n), A0, normtype);  % compute growth factor for this trial
        GammaVal(itest) = gamma0;                          % store gamma for summaries and histogram

        % --- Solve once using the factors (Ly = Pr*b, Uz = y, undo Pc) -------- 

        % wrap triangular solves in try/catch so a tiny diag(U) doesn't stop the sweep;
        %      on failure we mark this trial as failed (NaN metrics) and continue.
        
        solve_ok = true;                        % assume solves succeed
        solve_errmsg = '';                      % capture error string for printing outside catch
        try
          bhat = b(rp);                               % apply row permutation to b
          y    = forward_unit(M(1:n,1:n), bhat(1:n)); % solve L*y = bhat with unit-lower solver
          z    = back_sub(M(1:n,1:n), y, tol);        % solve U*z = y with back-sub routine
        catch ME
          solve_ok = false;                    % solve failed (e.g., "Zero/near-zero pivot at row k")
          solve_errmsg = ME.message;           % save error text
        end
        if ~solve_ok
          RelFactorErr(itest) = NaN;           % keep all metrics for this trial as NaN
          RelSolErr(itest)    = NaN;           % "
          RelResid(itest)     = NaN;           % "
          GammaVal(itest)     = NaN;           % "
          if outputtype == 1                   % optional notice for per-test output
            fprintf(1,'n=%4d  (Solve failed: %s)\n', n, solve_errmsg); % print failure reason
          end
          continue                             % go to next trial
        end

        if mode == 2                           % if complete pivoting (columns permuted)
          xhat(1:n) = 0;                       % initialize solution vector
          for k = 1:n                          % invert the column permutation
            xhat(cp(k)) = z(k);                % place z(k) into position cp(k) so xhat = P_c * z
          end
        else                                   % if no column permutation (modes 0 or 1)
          xhat(1:n) = z;                       % z is already in the correct order
        end

        % --- Solution metrics: rel solution error and rel residual ------------ 

        RelSolErr(itest) = norm(xtrue(1:n) - xhat(1:n),  normtype) / NX;   % relative solution error
        RelResid(itest)  = norm(b(1:n)     - A0*xhat(1:n), normtype) / NB; % relative residual

        % --- Optional per-test printing --------------------- 

        if outputtype == 1                     % print one line per test if requested
          fprintf(1,'n=%4d  relF=%9.2e  relX=%9.2e  relR=%9.2e  gamma=%9.2e\n', ...
                    n, RelFactorErr(itest), RelSolErr(itest), RelResid(itest), GammaVal(itest));
        end

        % --- Append gamma to the all-tests collector (for histogram) ---------- 

        GammaAll = [GammaAll; GammaVal(itest)]; % append (simple concatenation, fine at this scale)

        % --- Update overall maxima across everything ---------- 

        if ~isnan(RelFactorErr(itest))         % guard against NaNs from failed trials
          overall_max_relF = max(overall_max_relF, RelFactorErr(itest)); % track worst case so far
        end
        if ~isnan(RelSolErr(itest))
          overall_max_relX = max(overall_max_relX, RelSolErr(itest));    % "
        end
        if ~isnan(RelResid(itest))
          overall_max_relR = max(overall_max_relR, RelResid(itest));      % "
        end
        if ~isnan(GammaVal(itest))
          overall_max_gam  = max(overall_max_gam,  GammaVal(itest));      % "
        end

      end  % end loop over itest

      % --- After ntests, compute per-n mean and max for this mode ------------- 
      % NOTE: use 'omitnan' so failed trials (NaN) are ignored in the summary.

      mean_relF_modes(experimentindex, mi) = mean(RelFactorErr(1:ntests), 'omitnan'); % mean relF for this n/mode
      max_relF_modes(experimentindex,  mi) = max(RelFactorErr(1:ntests), [], 'omitnan'); % max relF for this n/mode
      mean_relX_modes(experimentindex, mi) = mean(RelSolErr(1:ntests),    'omitnan'); % mean relX "
      max_relX_modes(experimentindex,  mi) = max(RelSolErr(1:ntests),    [], 'omitnan'); % max relX "
      mean_relR_modes(experimentindex, mi) = mean(RelResid(1:ntests),     'omitnan'); % mean relR "
      max_relR_modes(experimentindex,  mi) = max(RelResid(1:ntests),     [], 'omitnan'); % max relR "
      mean_gam_modes(experimentindex,  mi) = mean(GammaVal(1:ntests),     'omitnan'); % mean gamma  "
      max_gam_modes(experimentindex,   mi) = max(GammaVal(1:ntests),     [], 'omitnan'); % max gamma  "

      % --- Optional per-n printing -------------------------- 

      if outputtype == 2
        fprintf(1,'n=%4d  mean[relF,relX,relR,gam]=[%9.2e %9.2e %9.2e %9.2e]  ', ...
                  n, mean_relF_modes(experimentindex,mi), mean_relX_modes(experimentindex,mi), ...
                     mean_relR_modes(experimentindex,mi), mean_gam_modes(experimentindex,mi)); % print means
        fprintf(1,'max[relF,relX,relR,gam]=[%9.2e %9.2e %9.2e %9.2e]\n', ...
                  max_relF_modes(experimentindex,mi), max_relX_modes(experimentindex,mi), ...
                  max_relR_modes(experimentindex,mi), max_gam_modes(experimentindex,mi)); % print maxima
      end

    end % end loop over n

    % --- If plotting requested, (REPLACED) per-mode gamma histogram ---------------
    
    if outputtype == 3
      GammaAll_plot = GammaAll(~isnan(GammaAll));    % keep only real numbers (drop failed trials)
      GammaAll_modes{mi} = GammaAll_plot;            % save this mode's gamma for the class-level plot
    end

  end % end loop over modes

  % ===================== COMPACT SUMMARY PLOTS PER CLASS ========================
  if outputtype == 3

    % --- Utility to floor at 1e-17 so log10 is safe (avoid -Inf)  ----------

    floorlog = @(y) max(-17, log10(max(y, 1e-17)));  % safe base-10 log for tiny values

    % --- class name for titles -----------------------------------

    switch ptype
      case 1, cname = 'Class 1 — A = L·U (full-rank product)';      % title text for class 1
      case 2, cname = 'Class 2 — Symmetric positive definite (SPD)';% title text for class 2
      case 3, cname = 'Class 3 — Strictly row-diagonally dominant'; % title text for class 3
      case 4, cname = 'Class 4 — Permuted row-diagonally dominant'; % title text for class 4
      otherwise, cname = sprintf('Class %d', ptype);                 % fallback label
    end

    % --- labels for pivoting modes (legend) -----------------------
    mode_labels = cell(1, nModes);                                   % preallocate legend labels
    for mi2 = 1:nModes                                               % fill in one label per tested mode
      switch modes_to_test(mi2)
        case 0, mode_labels{mi2} = 'Mode 0: no pivoting';            % name for mode 0
        case 1, mode_labels{mi2} = 'Mode 1: partial (row) pivoting'; % name for mode 1
        case 2, mode_labels{mi2} = 'Mode 2: complete (row+col) pivoting'; % name for mode 2
        otherwise, mode_labels{mi2} = sprintf('Mode %d', modes_to_test(mi2)); % catch-all text
      end
    end

    % === Figure 1: DASHBOARD — mean/max curves for relF, relX, relR ==========
    % Panels:
    %   Row 1: log10(mean relative factorization / solution / residual error)
    %   Row 2: log10(max  relative factorization / solution / residual error)
    fig1 = figure;                                                    % open a new figure window
    t1 = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');% 2x3 grid with tight spacing
    title(t1, sprintf('%s — Error metrics vs matrix dimension n (modes overlaid; base-10 log scale)', cname)); % figure title

    labelsY = { ...                                                   % y-axis text for each tile
      'log_{10} (mean relative factorization error)', ...
      'log_{10} (mean relative solution error)', ...
      'log_{10} (mean relative residual)'; ...
      'log_{10} (max relative factorization error)', ...
      'log_{10} (max relative solution error)', ...
      'log_{10} (max relative residual)'};

    % NOTE on definitions used throughout:
    %   relative factorization error  : ||Pr·A·Pc − L·U|| / max(||A||, 1)
    %   relative solution error       : ||x_true − x_hat|| / max(||x_true||, 1)
    %   relative residual             : ||b − A·x_hat|| / max(||b||, 1)

    dataMean = {mean_relF_modes, mean_relX_modes, mean_relR_modes};  % handles to mean arrays
    dataMax  = {max_relF_modes,  max_relX_modes,  max_relR_modes};   % handles to max arrays

    for row = 1:2                                                     % loop over the two rows (mean / max)
      for col = 1:3                                                   % loop over the three metrics (relF, relX, relR)
        nexttile; hold on; grid on;                                   % go to next panel; keep plots; show grid
        if row == 1
          A = dataMean{col};                                          % pick the mean data for this metric
        else
          A = dataMax{col};                                           % pick the max data for this metric
        end
        for mi = 1:nModes                                             % draw one line per pivoting mode
          plot(nlist, floorlog(A(:,mi)), 'o-', 'LineWidth', 1);       % plot log10 data vs n
        end
        xlabel('matrix dimension n');                                 % clear x-axis label
        ylabel(labelsY{row,col});                                     % clear y-axis label for this panel
        if (row == 1) && (col == 3)                                   % only add legend once (top-right panel)
          legend(mode_labels, 'Location','best', 'FontSize',9);       % legend with mode names
        end
        hold off;                                                     % release the panel
      end
    end

    % === Figure 2: GROWTH-FACTOR (gamma) SUMMARY =================================
    % Left:  log10(mean gamma) and log10(max gamma) vs n for each pivoting mode
    % Right: Overlaid histograms of gamma (probability density) by pivoting mode

    fig2 = figure;                                                    % open a second figure window
    t2 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');% 1x2 grid, tight layout
    title(t2, sprintf('%s — Growth factor \\gamma summary (\\gamma = || |L||U| || / ||A||)', cname)); % figure title

    % Left tile: mean (solid) and max (dashed) gamma vs n
    nexttile; hold on; grid on;                                       % go to left panel
    for mi = 1:nModes
      plot(nlist, floorlog(mean_gam_modes(:,mi)), 'o-',  'LineWidth',1); % draw mean gamma curve for this mode
    end
    for mi = 1:nModes
      plot(nlist, floorlog(max_gam_modes(:,mi)),  '--',  'LineWidth',1); % draw max gamma curve (dashed) for this mode
    end
    xlabel('matrix dimension n');                                     % x-axis text
    ylabel('log_{10} growth factor \gamma (solid = mean across tests; dashed = max)'); % y-axis text
    % Build explicit legend entries
    leg_entries = [ ...                                               % combine labels for mean and max per mode
      cellfun(@(s) ['mean \gamma — ' s], mode_labels, 'UniformOutput', false), ...
      cellfun(@(s) ['max \gamma — '  s], mode_labels, 'UniformOutput', false)];
    legend(leg_entries, 'Location','bestoutside','FontSize',9);       % show legend outside to keep panel clear
    hold off;                                                         % done with left panel

    % Right tile: Overlaid gamma histograms (shared bin edges)
    nexttile; hold on; grid on;                                       % go to right panel
    % Collect all gamma to set common edges robustly
    allg = [];                                                        % start with empty list
    for mi = 1:nModes
      if ~isempty(GammaAll_modes{mi})                                 % skip modes with no data
        allg = [allg; GammaAll_modes{mi}(:)];                         % append all gamma values for this mode
      end
    end
    edges = [];                                                       % default: let histogram pick bins
    if ~isempty(allg)                                                 % if we have any gamma at all
      gmin = min(allg); gmax = max(allg);                             % find range for gamma
      if isfinite(gmin) && isfinite(gmax) && (gmax > gmin)            % make sure the range is valid
        edges = linspace(gmin, gmax, 30);                             % build shared bin edges (30 bins)
      end
    end
    for mi = 1:nModes                                                 % draw one histogram per mode
      g = GammaAll_modes{mi};                                         % gamma values for this mode
      if ~isempty(g)
        if isempty(edges)
          histogram(g, 30, 'Normalization','pdf', 'DisplayStyle','stairs', 'LineWidth',1.2); % auto bins; area = 1
        else
          histogram(g, edges, 'Normalization','pdf', 'DisplayStyle','stairs', 'LineWidth',1.2); % shared bins for fair comparison
        end
      end
    end
    xlabel('\gamma (growth factor)');                                  % x-axis 
    ylabel('Probability density (normalized histogram)');              % y-axis 
    legend(mode_labels,'Location','best');                             % legend with mode names
    hold off;                                                          % done with right panel

  end % end plotting per class

end % end loop over classes

% ---- Print overall maxima across all classes/modes/sizes/tests --------------- 

fprintf(1,'\nOverall maxima across ALL tests:\n');                            % header for final summary
fprintf(1,'  max rel factorization error = %9.2e\n', overall_max_relF);       % report worst relF
fprintf(1,'  max rel solution error      = %9.2e\n', overall_max_relX);       % report worst relX
fprintf(1,'  max rel residual            = %9.2e\n', overall_max_relR);       % report worst relR
fprintf(1,'  max gamma                   = %9.2e\n', overall_max_gam);        % report worst gamma
