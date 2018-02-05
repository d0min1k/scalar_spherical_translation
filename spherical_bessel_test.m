x = linspace(0, 10, 200);

plot(x, sphbes(0, x));
hold on
for n = 1: 5
   plot(x, sphbes(n, x));
end    
    
grid

% [Passed]