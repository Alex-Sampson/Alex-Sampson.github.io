function run_hw2()
% RUN_HW2  Simple menu to run HW2 drivers.
% Make sure task1_driver.m, task2_driver.m, task3_driver.m, task4_driver.m
% are on the MATLAB path (ideally in the same folder).

clc;
disp('========== FCM HW2 – Driver Menu ==========');
disp(' 1) Task 1  (data-driven COO tests)');
disp(' 2) Task 2  (CSR-with-elbow + inserts; plots over n)');
disp(' 3) Task 3  (Symmetric CSR; plots over n)');
disp(' 4) Task 4  (DIA three bands; plots over n)');
disp(' 5) Run ALL (1 → 4)');
disp('===========================================');

choice = input('Select 0–5: ');

switch choice
    case 1
        task1_driver();
    case 2
        task2_driver();
    case 3
        task3_driver();
    case 4
        task4_driver();
    case 5
        task1_driver();
        task2_driver();
        task3_driver();
        task4_driver();
    otherwise
        disp('Invalid selection.');
end
end
