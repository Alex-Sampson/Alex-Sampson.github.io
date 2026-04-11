function task2_driver()
% TASK 2: unordered COO -> CSR+elbow; insertions; SpMV vs dense; plots over n
% Uses your signatures:
%   [AA,JA,IA,NEXT] = coo_to_csr_elbow(JR,JC,AAc,n,Kpad)
%   y = spmv_csr_elbow(AA,JA,IA,NEXT,x)
%   [AA,JA,NEXT,ok] = csr_insert(i,j,v,AA,JA,IA,NEXT)

rng(2);
n_list  = 30:15:900;   % <<< RANGE HERE
ntrials = 3;
kmax    = 6;
Kpad    = 4;

meanRel_before = zeros(numel(n_list),1);
meanRel_after  = zeros(numel(n_list),1);
meanT_sp_bef   = zeros(numel(n_list),1);
meanT_sp_aft   = zeros(numel(n_list),1);
meanT_dn_bef   = zeros(numel(n_list),1);
meanT_dn_aft   = zeros(numel(n_list),1);
meanInsOK      = zeros(numel(n_list),1);
meanCapFull    = zeros(numel(n_list),1);
meanDupSkip    = zeros(numel(n_list),1);

for ti = 1:numel(n_list)
    n = n_list(ti);

    rel_bef = zeros(ntrials,1);
    rel_aft = zeros(ntrials,1);
    tsp_bef = zeros(ntrials,1);
    tsp_aft = zeros(ntrials,1);
    tdn_bef = zeros(ntrials,1);
    tdn_aft = zeros(ntrials,1);
    insOK   = zeros(ntrials,1);
    capFull = zeros(ntrials,1);
    dupSkip = zeros(ntrials,1);

    for it = 1:ntrials
        % Unordered COO
        [JR, JC, AAc] = gen_coo_inrow(n, kmax);
        p = randperm(numel(AAc)); JR = JR(p); JC = JC(p); AAc = AAc(p);
        A = coo_to_dense(JR, JC, AAc, n);
        x = randn(n,1);

        % CSR-with-elbow
        [AA, JA, IA, NEXT] = coo_to_csr_elbow(JR, JC, AAc, n, Kpad);

        % BEFORE
        tic; ysp = spmv_csr_elbow(AA, JA, IA, NEXT, x); tsp_bef(it) = toc;
        tic; ydn = A*x;                                 tdn_bef(it) = toc;
        rel_bef(it) = relerr(ysp, ydn);

        % Insert random entries
        attempts = max(2*n, 60);  % simple scale
        ok=0; cap=0; dup=0;
        for s = 1:attempts
            i = randi(n); j = randi(n);
            if A(i,j) ~= 0, dup=dup+1; continue; end
            v = 2*rand()-1;
            [AA, JA, NEXT, good] = csr_insert(i, j, v, AA, JA, IA, NEXT);
            if good, A(i,j)=v; ok=ok+1; else, cap=cap+1; end
        end
        insOK(it)=ok; capFull(it)=cap; dupSkip(it)=dup;

        % AFTER
        tic; ysp2 = spmv_csr_elbow(AA, JA, IA, NEXT, x); tsp_aft(it) = toc;
        tic; ydn2 = A*x;                                 tdn_aft(it) = toc;
        rel_aft(it) = relerr(ysp2, ydn2);
    end

    meanRel_before(ti)=mean(rel_bef); meanRel_after(ti)=mean(rel_aft);
    meanT_sp_bef(ti)=mean(tsp_bef);   meanT_sp_aft(ti)=mean(tsp_aft);
    meanT_dn_bef(ti)=mean(tdn_bef);   meanT_dn_aft(ti)=mean(tdn_aft);
    meanInsOK(ti)=mean(insOK);        meanCapFull(ti)=mean(capFull); meanDupSkip(ti)=mean(dupSkip);
end

% ---------- plots ----------
figure('Name','Task 2 — Relative error');
semilogy(n_list, meanRel_before,'-o'); hold on;
semilogy(n_list, meanRel_after,'-s'); hold off; grid on;
xlabel('n'); ylabel('mean relative error');
legend('before insert','after insert','Location','best');
title('Task 2: CSR-with-elbow SpMV — error vs n');
saveas(gcf,'task2_relerr.png');

figure('Name','Task 2 — Timing');
plot(n_list, meanT_sp_bef,'-o'); hold on;
plot(n_list, meanT_sp_aft,'-s');
plot(n_list, meanT_dn_bef,'-^');
plot(n_list, meanT_dn_aft,'-d'); hold off; grid on;
xlabel('n'); ylabel('mean time (s)');
legend('SpMV before','SpMV after','Dense before','Dense after','Location','northwest');
title('Task 2: Timing vs n');
saveas(gcf,'task2_timing.png');

figure('Name','Task 2 — Insertion stats');
bar(n_list, [meanInsOK meanCapFull meanDupSkip], 'grouped');
xlabel('n'); ylabel('count (mean)'); grid on;
legend('ok inserts','capacity full','duplicates skipped','Location','northwest');
title('Task 2: insert outcomes per n');
saveas(gcf,'task2_inserts.png');

fprintf('Task 2: saved task2_relerr.png, task2_timing.png, task2_inserts.png\n');

% ---------- helpers ----------
function e = relerr(y, yref)
    nr = norm(yref); if nr==0, e=norm(y-yref); else, e=norm(y-yref)/nr; end
end
function [JR, JC, AAc] = gen_coo_inrow(n, kmax)
    JR=[]; JC=[]; AAc=[];
    for i=1:n
        rnnz = randi([0, min(kmax,n)]);
        if rnnz>0
            cols = randperm(n,rnnz);
            vals = 2*rand(rnnz,1)-1;
            JR=[JR; i*ones(rnnz,1)]; JC=[JC; cols(:)]; AAc=[AAc; vals];
        end
    end
    if isempty(AAc), JR=1; JC=1; AAc=1; end
end
function A = coo_to_dense(JR,JC,AAc,n)
    A = zeros(n,n); for k=1:numel(AAc), A(JR(k),JC(k))=AAc(k); end
end
end
