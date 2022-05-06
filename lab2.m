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
N = 1e6;
for i=1:length(gs)
    L = floor(N/ks(i));
    u = randi([0 1], ks(i), L);
    v = coder(gs(i), u);
    for l=length(probabilities)
        err = error(ns(i), size(u, 1), probabilities(l));
        base_synd = error(ns(i), 1, -1);
        r = mod(v + err, 2);
        s = decoder(r, gs(i));
        for m=1:size(r, 1)
            if s(m, :) ~= 0
                rotations = 0;
                while ~equal_base_synd(s(m, :), base_synd)
                    r(m, :) = circshift(r(m, :), 1);
                    s(m, :) = decoder(r(m, :), gs(i));
                    rotations = rotations+1;
                end
                r(m, 1) = ~r(m, :);
                r(m, :) = circshift(r(m, :), -rotations);
            end
        end
    end
end

function v = coder(g_cell, u)
    g = cell2mat(g_cell);
    v = zeros(size(u, 1), size(u, 2) + size(g, 2) - 1);
    for i=1:size(u, 1)
        v(i, :) = conv(g, u(i, :));
    end
end

function s = decoder(r, g_cell)
    g = cell2mat(g_cell);
    s = [];
    for i=1:size(r, 1)
        [~, synd] = deconv(r(i), g);
        s(i, :) = mod(synd, 2);
    end
end

function e = error(n, k2, p)
    e = zeros(k2, n);
    e(:, 1) = rand(k2, 1) > p;
end

function val = equal_base_synd(synd, base_synd)
    val = (synd == base_synd);
end

% item 3:
% função que recebe como parâmetros:u(D), que é a informação,
% e g(D)'s gerados no item 2 (tantos quantos forem os pares (n,k))

%item 4:
% após v(D) passar pelo canal, sai do canal r(D) = v(D) + e(D)

% o resto de r(D)/g(D) é o mesmo resto de e(D)/g(D) e se chama síndrome
% (s(D))

% compara-se essa síndrome com a síndrome de quando há erro na
% primeira posição de r(D)

% se a comparação der verdadeira, significa que há erro na primeira posição
% de r(D). Troca-se o primeiro bit de r(D) (para corrigir o erro) e tem-se
% o vetor desejado

% se a comparação der falsa, significa que não há erro na primeira posição
% de r(D). Nesse caso, rotaciona-se r(D), obtendo-se r_1(D). 
% Novamente compara-se s_1(D) com o conjunto de síndromes de quando
% há erro na primeira posição. Caso encontre erro, corrige-se o erro e esse é o
% vetor desejado. Caso não encontre erro, o processo de rotação continua até
% o erro ser encontrado

% é preciso armazenar o número de rotações feitas em uma variável para
% poder "desrotacionar" o vetor e recolocar as posições dele (corrigidas)
% nos devidos lugares

% precisa-se de funções para rotacionar r(D) em relação a (1+D^N) e s(D)
% em relação a g(D) 