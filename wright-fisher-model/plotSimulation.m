function plotSimulation(gen, dataSimulation)

    n = size(dataSimulation, 1);
    hold off
    for i = 1:n
        plot(gen, dataSimulation(i, :))
        hold on
    end

end