function u = FDMWithSource(xmin, xmax, nx, dx, ymin, ymax, ny, dy, nt, dt, x, y, c, BC)
    % Initialize u
    u = zeros(nt, nx, ny);
    u(1, :, :) = 0; % Initial condition is set to zero

   function s = source_function(x, y, t)
       x_s = 2.5;
       y_s = 2.5;
       xx_s = x; %x_d(2:end-1, 2:end-1);
       yy_s = y; %y_d(2:end-1, 2:end-1);
       sigma_s = 0.001;
       x_part_s = (xx_s-x_s).^2/(2*sigma_s.^2);
       y_part_s = (yy_s-y_s).^2/(2*sigma_s.^2);
       source_location = exp(-(x_part_s + y_part_s));

       s = 20*sin(30*pi*t/20)*source_location;
       mesh(s )
   end

    % Apply source function to initial time step
%     u(1, 51, 51) = source_function(x(51, 51), y(51, 51),0);

    tic;
    for n = 2:nt-1
        if mod(n, 100) == 0
            fprintf('>>>>> FDM computing... n: %d, nt: %d\n', n, nt);
        end
        for i = 2:nx-1
            for j = 2:ny-1
                u(n + 1, i, j) = 2 * u(n, i, j) - u(n - 1, i, j) ...
                                + (c * dt / dx)^2 * (u(n, i + 1, j) - 2 * u(n, i, j) + u(n, i - 1, j)) ...
                                + (c * dt / dy)^2 * (u(n, i, j + 1) - 2 * u(n, i, j) + u(n, i, j - 1)) ;
          
            end
        end

       % u(n + 1, 51, 51) =  dt.^2 * (20*sin(30*pi*(n*dt)/20)); %* exp(-(x - 2).^2) .* exp(-(y - 2).^2);
        u(n + 1, :, :)  =  dt.^2 * source_function(x, y , n*dt);

        if strcmp(BC, 'Dir')
            for i = 2:nx-1
                u(n + 1, i,  1) = 0;
                u(n + 1, i, end) = 0;
            end
            for j = 2:ny-1
                u(n + 1,  1, j) = 0;
                u(n + 1, end, j) = 0;
            end
            u(n,  1,  1) = 0;
            u(n, end,  1) = 0;
            u(n,  1, end) = 0;
            u(n, end, end) = 0;
        elseif strcmp(BC, 'Neu')
            for i = 2:nx-1
                u(n + 1, i,  1) = u(n + 1, i,  2);
                u(n + 1, i, end) = u(n + 1, i, end - 1);
            end
            for j = 2:ny-1
                u(n + 1,  1, j) = u(n + 1,  2, j);
                u(n + 1, end, j) = u(n + 1, end - 1, j);
            end
            u(n,  1,  1) = (u(n,  2,  1) + u(n,  1,  2)) / 2;
            u(n, end,  1) = (u(n, end - 1,  1) + u(n, end,  2)) / 2;
            u(n,  1, end) = (u(n,  2, end) + u(n,  1, end - 1)) / 2;
            u(n, end, end) = (u(n, end - 1, end) + u(n, end, end - 1)) / 2;
        end
    end
    t1 = toc;
    fprintf('>>>>> elapse time for FDM (sec): %f\n', t1);
    fprintf('>>>>> elapse time for FDM (min): %f\n', t1 / 60);

    % Plot FDM solutions
    figure('Position', [100, 100, 1600, 400]);
    subplot(1, 3, 1);
    surf(x, y, squeeze(u(65, :, :)), 'EdgeColor', 'none');
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);
%     zlim([-1, 1]);
    xlabel('x');
    ylabel('y');
    zlabel('u (x, y)');
    title('Snapshot 1');

    subplot(1, 3, 2);
    surf(x, y, squeeze(u(95, :, :)), 'EdgeColor', 'none');
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);
%     zlim([-1, 1]);
    xlabel('x');
    ylabel('y');
    zlabel('u (x, y)');
    title('Snapshot 2');

    subplot(1, 3, 3);
    surf(x, y, squeeze(u(170, :, :)), 'EdgeColor', 'none');
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);
%     zlim([-1, 1]);
    xlabel('x');
    ylabel('y');
    zlabel('u (x, y)');
    title('Snapshot 3');

    sgtitle('FDM Solutions');
end

