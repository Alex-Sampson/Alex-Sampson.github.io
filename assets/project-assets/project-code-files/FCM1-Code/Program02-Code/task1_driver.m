function task1_driver()
% TASK 1 DRIVER (data-driven, with plots)
% Uses your provided COO triplets (JR, JC, AA) and vectors v, v2, v3.
% Compares CSR / modified CSR / ELLPACK SpMV against dense A*x.
% Saves two figures per matrix: *_relerr.png and *_timing.png

rng(1);

% =======================
% --- Test Problem 1 ---
AA = [10,-1,3,11,5,12,7,13,9,2,3];
JR = [1,1,2,2,3,3,4,4,5,5,5];
JC = [1,5,1,2,2,3,3,4,4,1,5]; % column out-of-order at the end
v  = [1; 0; -1; 0; 2];

run_case('Matrix1', JR, JC, AA, v);

% =======================
% --- Test Problem 2 ---
AA2 = [1, 5, -2, 4, 1, 2, 1, 2, 3];
JR2 = [1, 1, 2, 2, 3, 3, 4, 4, 4];
JC2 = [1, 4, 1, 2, 1, 4, 2, 3, 4];
v2  = [3; 5; -1; -1];

run_case('Matrix2', JR2, JC2, AA2, v2);

% =======================
% --- Test Problem 3 ---
AA3 = [1, 2, 3, 2, -3, 1, 2, 3, -1, 1, 2, 1, -1, 4, 1, 1, -4, 5];
JR3 = [1, 1, 1, 2, 2, 3, 3, 4, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8];
JC3 = [1, 3, 5, 2, 6, 3, 5, 4, 2, 5, 3, 4, 5, 6, 2, 7, 4, 8];
v3  = [-2; 2; -2; 2; -2; 2; -2; 2];

run_case('Matrix3', JR3, JC3, AA3, v3);

fprintf('\nTask 1: finished. Plots saved next to this script.\n');

% ======================================================================
function run_case(tag, JR, JC, AA, x)
    n = max([JR(:); JC(:)]);              % infer dimension
    A = coo_to_dense(JR, JC, AA, n);      % dense reference
    x = x(:);
    if length(x) ~= n
        error('%s: length(x)=%d but n=%d inferred from indices.', tag, length(x), n);
    end

    % ---------- CSR ----------
    [AA_csr, JA_csr, IA_csr] = toCSR(JR, JC, AA, n);
    t_csr = tic; y_csr = spmv_csr(IA_csr, JA_csr, AA_csr, x); t_csr = toc(t_csr);

    % ---------- modified CSR ----------
    [AA_m, JA_m, IA_m, D_m] = toModCSR(JR, JC, AA, n);
    t_m = tic; y_m = spmv_modcsr(IA_m, JA_m, AA_m, D_m, x); t_m = toc(t_m);

    % ---------- ELLPACK ----------
    [VAL, COL, K] = toELLPACK(JR, JC, AA, n);
    t_ell = tic; y_ell = spmv_ellpack(VAL, COL, K, x); t_ell = toc(t_ell);

    % ---------- Dense reference ----------
    t_dn = tic; y_ref = A*x; t_dn = toc(t_dn);

    % ---------- Errors ----------
    rel = [relerr(y_csr,y_ref), relerr(y_m,y_ref), relerr(y_ell,y_ref)];
    fprintf('%s  |  rel.err  CSR=%.2e  mCSR=%.2e  ELL=%.2e\n', tag, rel(1), rel(2), rel(3));

    % ---------- Plots ----------
    % Errors
    f1 = figure('Name',[tag ' — Relative error'],'Visible','off');
    bar(categorical({'CSR','mCSR','ELLPACK'}), rel);
    ylabel('relative error'); grid on;
    title([tag ': Relative error (vs dense)']);
    saveas(f1, [lower(tag) '_relerr.png']); close(f1);

    % Timing
    f2 = figure('Name',[tag ' — Timing'],'Visible','off');
    bar(categorical({'CSR SpMV','mCSR SpMV','ELLPACK SpMV','Dense A*x'}), [t_csr t_m t_ell t_dn]);
    ylabel('time (s)'); grid on;
    title([tag ': Timing']);
    saveas(f2, [lower(tag) '_timing.png']); close(f2);
end

% ======================================================================
% ---------- helper: dense from COO ----------
function A = coo_to_dense(JR, JC, AA, n)
    A = zeros(n,n);
    for k = 1:numel(AA)
        A(JR(k), JC(k)) = A(JR(k), JC(k)) + AA(k);
    end
end

% ---------- helper: robust relative error ----------
function e = relerr(y, yref)
    nr = norm(yref);
    if nr == 0, e = norm(y - yref);
    else,        e = norm(y - yref)/nr;
    end
end

% ---------- wrappers that tolerate output-order differences ----------
function [AA_csr, JA_csr, IA_csr] = toCSR(I,J,V,n)
    % Try common orders: [AA,JA,IA] then [IA,JA,AA]
    try
        [AA_csr, JA_csr, IA_csr] = coo_to_csr(I, J, V, n);
        % sanity check sizes
        assert(numel(IA_csr)==max(I) + 1 || numel(IA_csr)==n+1);
    catch
        [IA_csr, JA_csr, AA_csr] = coo_to_csr(I, J, V, n);
    end
end

function [AA_m, JA_m, IA_m, D_m] = toModCSR(I,J,V,n)
    % Try [AA,JA,IA,Diag] then [IA,JA,AA,Diag]
    try
        [AA_m, JA_m, IA_m, D_m] = coo_to_modcsr(I, J, V, n);
        assert(numel(IA_m)==n+1);
    catch
        [IA_m, JA_m, AA_m, D_m] = coo_to_modcsr(I, J, V, n);
    end
end

function [VAL, COL, K] = toELLPACK(I,J,V,n)
    [VAL, COL, K] = coo_to_ellpack(I, J, V, n);
end
end
