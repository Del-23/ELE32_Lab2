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


%item 3:
% função que recebe como parâmetros:u(D), que é a informação,
% e g(D)'s gerados no item 2 (tantos quantos forem os pares (n,k))

%item 4:
% após v(D) passar pelo canal, sai do canal r(D) = v(D) + e(D)

% o resto de r(D)/g(D) é o mesmo resto de e(D)/g(D) e se chama síndrome
% (s(D))

% compara-se essa síndrome com o conjunto de síndromes estabelecidas como
% as síndromes de quando há erro na primeira posição de r(D)

% se a comparação der verdadeira, significa que há erro na primeira posição
% de r(D). Troca-se o primeiro bit de r(D) (para corrigir o erro) e tem-se
% o vetor desejado

% se a comparação der falsa, significa que não há erro na primeira posição
%de r(D). Nesse caso, rotaciona-se r(D) e s(D), obtendo-se r_1(D) e
%s_1(D). Novamente compara-se s_1(D) com o conjunto de síndromes de quando
%há erro na primeira posição. Caso encontre erro, corrige-se o erro e esse é o
%vetor desejado. Caso não encontre erro, o processo de rotação continua até
%o erro ser encontrado

% é preciso armazenar o número de rotações feitas em uma variável para
% poder "desrotacionar" o vetor e recolocar as posições dele (corrigidas)
% nos devidos lugares

%precisa-se de funções para rotacionar r(D) em relação a (1+D^N) e s(D)
%em relação a g(D) 