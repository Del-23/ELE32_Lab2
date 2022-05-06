lgd = cell(size(probabilities_biterr, 1)+2, 1);
lgd{1} = "Sem codificação";
loglog(probabilities, probabilities)
hold on
lgd{2} = "Hamming(7,4)";
loglog(probabilities, prob_biterr_hamm)
for i=1:size(probabilities_biterr, 1)
    loglog(probabilities, probabilities_biterr(i, :))
    lgd{i+2} = strcat("n=", num2str(ns(i)), ", k=", num2str(ks(i)));
end
legend(lgd)
xlabel("Probabilidade de erro de bit no canal")
ylabel("Probabilidade de erro de bit na saída")
set(gca, 'xdir', 'reverse')
