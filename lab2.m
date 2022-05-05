ns = [9 11 12 14 16 18 20];
ks = [5 6 7 8 9 10 11];
polynomials = containers.Map('KeyType', 'int32', 'ValueType', 'any');

for i=1:length(ns)
    fprintf("n = %d\n", ns(i));
    fprintf("k = %d\n", ks(i));
    polynomials(i) = cyclpoly(ns(i), ks(i), 'all')
end

% split 