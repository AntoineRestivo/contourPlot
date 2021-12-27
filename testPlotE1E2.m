%O = compileE1E2(sdpsettings('solver', 'linprog'));
O = compileE1E2(sdpsettings('solver', 'mosek', 'verbose', 0, 'debug', 0));
pt = [1 1];
C = ContourPlot(@(x,y) O(x,y), [-1.2 1.2 -1.2 1.2], pt);
figure;
hold on;
axis([-1.2 1.2 -1.2 1.2]);
i = 1;
while ~isempty(C.todo)
    C.step;
    if mod(i, 20) == 1
        drawnow
    end
    i = i + 1;
end

savefig('Plot.fig')