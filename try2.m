n = 200; % Matrix size
A = rand(n); % Generate a large random matrix
B = rand(n); % Generate another matrix

tic;
A1=A';
elapsed_time = toc;

num_operations = 2 * n^3; % Number of floating-point operations in matrix multiplication
flops = num_operations / elapsed_time;

disp(['Estimated FLOPS: ', num2str(flops), ' FLOPS']);
