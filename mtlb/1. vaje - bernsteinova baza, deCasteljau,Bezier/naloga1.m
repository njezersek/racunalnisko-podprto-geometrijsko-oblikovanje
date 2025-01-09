% 1. naloga:

% 1.1)
% polinom v potenčni bazi
p = [4 7 2 8 9 1 4];

% prevedba v Bernsteinovo bazo
b = power2bernstein(p)

% prevedba nazaj v potenčno bazo
p = bernstein2power(b)

% 1.2)
power2bernstein(1)
power2bernstein([0 1])
power2bernstein([0 0 1])
power2bernstein([1 0])
power2bernstein([0 1 0])
power2bernstein([0 0 1 0])

% 1.3)
n = 5;
x = linspace(0,1);
clf
hold on
for i = 0:n
    b = zeros(1,n+1);
    b(i+1) = 1;
    plot(x,polyval(bernstein2power(b),x));
end
hold off