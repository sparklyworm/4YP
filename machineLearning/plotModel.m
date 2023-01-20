function plotModel(x_data, y_data, true_y, model_y, model_f)
    N = length(model_f);
    max_x = max(x_data);
    x_plot = linspace(1,max_x,N)';

    n = length(model_y);
    sd = sum((y_data - model_y).^2)/n;
    top_sd = model_y + sd;
    bottom_sd = model_y - sd;

    figure
    patch('Faces', linspace(1,2*n,2*n),'Vertices',[[x_data; flip(x_data)], [top_sd; flip(bottom_sd)]], 'FaceColor', 'magenta','EdgeColor', 'none','FaceAlpha', '0.4')
    hold on
    plot(x_plot, model_f, 'magenta', 'LineWidth',2);
    hold on
    plot(x_data, y_data, 'bo')
    hold on
    plot(x_plot, true_y, 'g')

end