function x = dd(xs, i)
% prema diferenca
% xs ... seznam
n = length(xs);
x = 0;
if 1 <= i && i < n
    x = xs(i+1) - xs(i);
end
end

