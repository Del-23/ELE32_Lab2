ns = [12 14 15 18 20];
ks = [7 8 9 10 12];
distances = zeros(1, length(ns));

for i=1:length(ns)
    fprintf("n = %d\n", ns(i));
    fprintf("k = %d\n\n", ks(i));
    gs = cyclpoly(ns(i), ks(i), 'all');
    for j=1:size(gs, 1)
        g = gs(j, :);
        c = [g(1) zeros(1, ks(i) - 1)];
        r = [g zeros(1, ks(i) - 1)];
        G = toeplitz(c, r);
        u = ff2n(ks(i));
        u = u(2:end, :);
        v = mod(u * G, 2);
        distance = min(sum(v, 2));
        if distance > distances(i)
            distances(i) = distance;
        end
    end
end


% split 