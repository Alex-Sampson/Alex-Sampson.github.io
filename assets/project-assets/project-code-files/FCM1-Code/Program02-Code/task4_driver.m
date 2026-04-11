function task4_driver()
% TASK 4: DIA with three bands [ -k, 0, +k ] — sweep n with fixed k list
% Uses:
%   [DIAG,IOFF] = DIA_triple(n,k,d_main,d_sub,d_super)
%   y = DIA_spmv(DIAG,IOFF,x)

rng(4);
n_list  = 30:15:900;         % <<< RANGE HERE
k_list  = [1 2 4 8 16 32];   % change if desired
ntrials = 3;

% We'll make one figure per metric (curves for each k over n)
meanRel = zeros(numel(n_list), numel(k_list));
meanTsp = zeros(numel(n_list), numel(k_list));
meanTdn = zeros(numel(n_list), numel(k_list));
effNnz  = zeros(numel(n_list), numel(k_list));

for ti = 1:numel(n_list)
    n = n_list(ti);
    for kj = 1:numel(k_list)
        k = k_list(kj);
        if k >= n, meanRel(ti,kj)=NaN; meanTsp(ti,kj)=NaN; meanTdn(ti,kj)=NaN; effNnz(ti,kj)=NaN; continue; end

        rels=zeros(ntrials,1); tsp=zeros(ntrials,1); tdn=zeros(ntrials,1);

        for it = 1:ntrials
            d_main  = 1 + randn(n,1);
            d_sub   = randn(n-k,1);
            d_super = randn(n-k,1);
            [DIAG, IOFF] = DIA_triple(n, k, d_main, d_sub, d_super);
            x = randn(n,1);

            tic; ysp = DIA_spmv(DIAG, IOFF, x); tsp(it) = toc;

            C = diag(d_main,0) + diag(d_super,+k) + diag(d_sub,-k);
            tic; ydn = C*x; tdn(it) = toc;

            rels(it) = relerr(ysp, ydn);
        end

        meanRel(ti,kj)=mean(rels);
        meanTsp(ti,kj)=mean(tsp);
        meanTdn(ti,kj)=mean(tdn);
        effNnz(ti,kj) = nnz(DIAG(:,1)) + nnz(DIAG(:,2)) + nnz(DIAG(:,3));
    end
end

% ---------- plots ----------
figure('Name','Task 4 — Relative error vs n');
hold on;
for kj=1:numel(k_list)
    semilogy(n_list, meanRel(:,kj),'-o');
end
hold off; grid on; xlabel('n'); ylabel('mean relative error');
legend(compose('k=%d',k_list),'Location','best');
title('Task 4: DIA SpMV — error vs n');
saveas(gcf,'task4_relerr_vs_n.png');

figure('Name','Task 4 — SpMV timing vs n');
hold on;
for kj=1:numel(k_list)
    plot(n_list, meanTsp(:,kj),'-o');
end
hold off; grid on; xlabel('n'); ylabel('mean time (s)');
legend(compose('k=%d',k_list),'Location','northwest');
title('Task 4: DIA SpMV — timing vs n');
saveas(gcf,'task4_spmv_timing_vs_n.png');

figure('Name','Task 4 — Dense timing vs n');
hold on;
for kj=1:numel(k_list)
    plot(n_list, meanTdn(:,kj),'-s');
end
hold off; grid on; xlabel('n'); ylabel('mean time (s)');
legend(compose('k=%d',k_list),'Location','northwest');
title('Task 4: Dense C*x — timing vs n');
saveas(gcf,'task4_dense_timing_vs_n.png');

figure('Name','Task 4 — Effective nnz vs n');
hold on;
for kj=1:numel(k_list)
    plot(n_list, effNnz(:,kj),'-d');
end
hold off; grid on; xlabel('n'); ylabel('effective nnz in stored bands');
legend(compose('k=%d',k_list),'Location','best');
title('Task 4: DIA nnz vs n');
saveas(gcf,'task4_nnz_vs_n.png');

fprintf('Task 4: saved task4_relerr_vs_n.png, task4_spmv_timing_vs_n.png, task4_dense_timing_vs_n.png, task4_nnz_vs_n.png\n');

function e = relerr(y,yref)
    nr=norm(yref); if nr==0, e=norm(y-yref); else, e=norm(y-yref)/nr; end
end
end
