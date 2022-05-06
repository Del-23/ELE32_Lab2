N = 1e6;
L = floor(N / 4);
x = randi([0 1], 4, L);
probabilities = [0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001 5e-4 2e-4 1e-4];
prob_biterr_hamm = zeros(size(probabilities));

Gt = [1 1 0 1;
     1 0 1 1;
     1 0 0 0;
     0 1 1 1;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];

H = [1 0 1 0 1 0 1;
     0 1 1 0 0 1 1;
     0 0 0 1 1 1 1];

R = [0 0 1 0 0 0 0;
     0 0 0 0 1 0 0;
     0 0 0 0 0 1 0;
     0 0 0 0 0 0 1];


for i = 1:length(probabilities)
    z = Gt * x;
    zmod = mod(z, 2);
    nz = binarychannel(zmod, probabilities(i));
    y = H * nz;
    ymod = mod(y, 2);
    rows_ymod = size(ymod, 1);
    for col=1:L
        num = calcnum(ymod, rows_ymod, col);
        if num ~= 0
            nz(num, col) = ~nz(num, col);
        end
    end
    y = R * nz;
    prob_biterr_hamm(i) = double(1 - sum(y == x, 'all') / (L*4));
end

function x = binarychannel(data, prob)
    data_size = size(data);
    rand_vals = rand(data_size);
    err = double(rand_vals < prob);
    x = cast(xor(data, err), 'like', data);
end

function num = calcnum(ymod, rows, col)
    num = 0;
    for i = 1:rows
        num = num + ymod(i, col) * 2^(i-1);
    end
end