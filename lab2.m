ns = [12 14 15 18 20];
ks = [7 8 9 10 12];
distances = zeros(1, length(ns));
polygs = {};

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
            polygs(i) = {g};
        end
    end
end

probabilities = [0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001 5e-4 2e-4 1e-4];
probabilities_biterr = zeros(size(ns, 2), size(probabilities, 2));
N = 1e5;
for i=1:length(polygs)
    L = floor(N/ks(i));
    u = randi([0 1], L, ks(i));
    v = coder(polygs(i), u);
    for l=1:length(probabilities)
        err = rand(size(u, 1), ns(i)) < probabilities(l);
        base_synd = decoder(error_first(ns(i), 1, -1), polygs(i));
        r = mod(v + err, 2);
        s = decoder(r, polygs(i));
        for m=1:size(r, 1)
            if any(s(m, :))
                rotations = 0;
                while ~equal_base_synd(s(m, :), base_synd) && rotations < ns(i)
                    r(m, :) = circshift(r(m, :), 1);
                    s(m, :) = decoder(r(m, :), polygs(i));
                    rotations = rotations + 1;
                end
                r(m, 1) = ~r(m, 1);
                r(m, :) = circshift(r(m, :), -rotations);
            end
        end
        probabilities_biterr(i, l) = double(1 - sum(r == v, 'all')/numel(r));
    end
end

function v = coder(g_cell, u)
    g = cell2mat(g_cell);
    v = zeros(size(u, 1), size(u, 2) + size(g, 2) - 1);
    for i=1:size(u, 1)
        v(i, :) = mod(conv(g, u(i, :)), 2);
    end
end

function s = decoder(r, g_cell)
    g = cell2mat(g_cell);
    s = [];
    for i=1:size(r, 1)
        [~, synd] = deconv(r(i, :), g);
        s(i, :) = mod(synd, 2);
    end
end

function e = error_first(n, k2, p)
    e = zeros(k2, n);
    e(:, 1) = rand(k2, 1) > p;
end

function val = equal_base_synd(synd, base_synd)
    val = all(synd == base_synd);
end

% item 3:
% fun????o que recebe como par??metros:u(D), que ?? a informa????o,
% e g(D)'s gerados no item 2 (tantos quantos forem os pares (n,k))

%item 4:
% ap??s v(D) passar pelo canal, sai do canal r(D) = v(D) + e(D)

% o resto de r(D)/g(D) ?? o mesmo resto de e(D)/g(D) e se chama s??ndrome
% (s(D))

% compara-se essa s??ndrome com a s??ndrome de quando h?? erro na
% primeira posi????o de r(D)

% se a compara????o der verdadeira, significa que h?? erro na primeira posi????o
% de r(D). Troca-se o primeiro bit de r(D) (para corrigir o erro) e tem-se
% o vetor desejado

% se a compara????o der falsa, significa que n??o h?? erro na primeira posi????o
% de r(D). Nesse caso, rotaciona-se r(D), obtendo-se r_1(D). 
% Novamente compara-se s_1(D) com o conjunto de s??ndromes de quando
% h?? erro na primeira posi????o. Caso encontre erro, corrige-se o erro e esse ?? o
% vetor desejado. Caso n??o encontre erro, o processo de rota????o continua at??
% o erro ser encontrado

% ?? preciso armazenar o n??mero de rota????es feitas em uma vari??vel para
% poder "desrotacionar" o vetor e recolocar as posi????es dele (corrigidas)
% nos devidos lugares

% precisa-se de fun????es para rotacionar r(D) em rela????o a (1+D^N) e s(D)
% em rela????o a g(D) 