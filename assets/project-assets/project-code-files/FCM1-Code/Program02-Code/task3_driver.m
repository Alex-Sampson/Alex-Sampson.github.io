function task3_driver()
% TASK 3: symmetric CSR(lower+diag) — error & timing over n with plots
% Uses:
%   [IA,JA,AA] = Sym_to_CSR(B)
%   y = Sym_CSR_spmv(IA,JA,AA,x)

rng(3);
n_list  = 30:15:900;   % <<< RANGE HERE
ntrials = 3;
kmax    = 6;

meanRel = zeros(numel(n_list),1);
meanTsp = zeros(numel(n_list),1);
meanTdn = zeros(numel(n_list),1);

for ti = 1:numel(n_list)
    n = n_list(ti);
    rels = zeros(ntrials,1); tsp=zeros(ntrials,1); tdn=zeros(ntrials,1);

    for it = 1:ntrials
        % build sparse symmetric B
        B = zeros(n,n);
        for i = 1:n
            rnnz = randi([1, min(kmax,i)]);
            cols = randperm(i, rnnz);
            vals = 2*rand(rnnz,1)-1;
            for s=1:rnnz
                j = cols(s); v = vals(s);
                B(i,j) = B(i,j) + v; B(j,i) = B(i,j);
            end
        end

        x = randn(n,1);
        [IA, JA, AA] = Sym_to_CSR(B);

        tic; ysp = Sym_CSR_spmv(IA, JA, AA, x); tsp(it) = toc;
        tic; ydn = B*x;                         tdn(it) = toc;
        rels(it) = relerr(ysp, ydn);
    end

    meanRel(ti)=mean(rels); meanTsp(ti)=mean(tsp); meanTdn(ti)=mean(tdn);
end

% plots
figure('Name','Task 3 — Relative error');
semilogy(n_list, meanRel,'-o'); grid on;
xlabel('n'); ylabel('mean relative error');
title('Task 3: Symmetric CSR SpMV — error vs n');
saveas(gcf,'task3_relerr.png');

figure('Name','Task 3 — Timing');
plot(n_list, meanTsp,'-o'); hold on;
plot(n_list, meanTdn,'-s'); hold off; grid on;
xlabel('n'); ylabel('mean time (s)');
legend('Sym-CSR SpMV','Dense B*x','Location','northwest');
title('Task 3: Timing vs n');
saveas(gcf,'task3_timing.png');

fprintf('Task 3: saved task3_relerr.png and task3_timing.png\n');

function e = relerr(y,yref)
    nr=norm(yref); if nr==0, e=norm(y-yref); else, e=norm(y-yref)/nr; end
end
end
