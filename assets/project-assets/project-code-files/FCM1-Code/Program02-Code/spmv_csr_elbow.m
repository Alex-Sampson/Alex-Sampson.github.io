function y = spmv_csr_elbow(AA, JA, IA, NEXT, x)
% SPMV_CSR_ELBOW  y = A*x using CSR-with-elbow (only used part of each row).
    n = length(NEXT);     % tells us how many rows are in AA
    y = zeros(n,1);       % initialize y as a zero storage array
    for i = 1:n           % loop over the number of rows.
        s = 0;            % initialize sum at zero
        for k = IA(i):NEXT(i)-1     % stop before unallocated storage
            s = s + AA(k)*x(JA(k)); % same sum as before.
        end
        y(i) = s;         % return i-th entry of y.
    end
end