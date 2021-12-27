pt = [0 1];
C = ContourPlot(@(x,y) 1-x^2-y^2, [-1.1 1.1 -1.1 1.1], pt);
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
