function [] = coherent_sets(alpha, v, epsilon, dom_ev, dom_vec, G_mode_dic, method, threshold, name, cumulative)
    % Computes, animates and evaluates coherent sets for a given
    % vectorfield v, based on the eigenfunction dom_vec of the augmented
    % generator with eigenvalue dom_ev. Visualization only works for p_dim=2.
    % After the simulation is done, there is the option to compute the
    % cumulative survival probabiliy. 
    %
    % There are a number of parameters that can be set for the simulation
    % below.
    %
    % The input is:
    % alpha: d_dim x 1 array indicating the quasiperiodic rotation direction.
    % v: (d_dim+p_dim) x n -> p_dim x n function handle of a divergence-free 
    % vectorfield. 
    % epsilon: The strength of the noise.
    % dom_ev: A dominant eigenvalue of the discrete generator G. 
    % dom_ev: A corresponding dominant right eigenvector of the discrete
    % generator G. 
    % G_mode_dic: A mode dictionary of the Fourier mode used in the galerkin
    % approximation of the generator.
    % method: a string that selects the method to compute coherent sets:
    %   "positive": Identifies the positive support as a coherent set.
    %   "real": Identifies points with absolute real part above threshold as 
    %         a coherent set.
    %   "abs": Identifies points with absolute value above threshold as a 
    %        coherent set.
    %   If 'positive' is chosen, threshold has no relevance.
    % name: The name of the example as a string.
    % cumulative: A boolean variable whether or not to compute the
    % cumulative survival probability. WARNING: setting this to true may be
    % very time consuming.

    d_dim = length(alpha);
    p_dim = size(G_mode_dic{1}, 2) - d_dim;

    % Setting the parameters:
    theta_0 = zeros(d_dim, 1); % The initial point in the driving system
    particle_res = 150; % Resolution of the particles in each dimension
    particle_size = 30; % Size of the particles in the scatter plots.
    T = 10; % Timespan of the simulation.
    dt = 0.001; % Time resolution of the simulation.
    render_time = 0.02; % Time resolution for the resulting animation.
    snapshots = 3; % Number of snapshots taken from the simulation and saved as a pdf.
    % The colors of points sith status 1, 2, 3
    colors = [0.54 0.54 0.54; 0, 0.5, 1; 0.54, 0.27, 0.27]; 

    % Initialization
    frames = ceil(T / render_time);
    Gid2m = G_mode_dic{1};
    Gm2id = G_mode_dic{2};

    % Initialize positions
    particles = particle_res^2;
    I = linspace(0,1, particle_res + 1);
    I = I(1:end-1);
    [X, Y] = ndgrid(I, I);

    ppos = [X(:)'; Y(:)'];  % p_dim x particles array of particle positions.
    theta_0 = mod(theta_0, 1);
    theta = theta_0;        % d_dim x 1 array, the positition of the driving.

    % Basepoints for contours
    contour_res = 40;
    I = linspace(0,1, contour_res);
    [cX, cY] = ndgrid(I, I);
    cpos = [cX(:)'; cY(:)'];

    % Boundary times for the driving plot
    if d_dim == 2
        BT1 = ((1: ceil(T * alpha(1))) - theta(1)) / alpha(1);
        BT2 = ((1: ceil(T * alpha(2))) - theta(2)) / alpha(2);
        BT = sort([BT1, BT2]);
    end
    
    % Initializing the coherent set. During the simulation each particle
    % has one of 3 status:
    % Status 1: Started outside the coherent set.
    % Status 2: Started inside the coherent set and is still inside.
    % Status 3: Started inside the coherent set but has left at some point.
    sparsity = abs(dom_vec(:, 1)) > 10^-12; 
    sF = transpose(dom_vec(sparsity, 1));
    sGid2m = Gid2m(sparsity, :);

    cf = sF * exp(2*pi*1i*(sGid2m(:, 1:d_dim)*theta_0 + sGid2m(:, d_dim+1:end)*cpos));
    z = reshape(cf, contour_res, contour_res);
    interf = griddedInterpolant(cX,cY,z);
    f = interf(ppos');

    status = ones(particles, 1);    
    switch method
        case "positive"
            f = real(f);
            status(f > 0) = 2;
        case "real"
            f = real(f) * length(cf) / sum(abs(real(cf)));
            status(abs(f) > threshold) = 2;
        case "abs"
            f = abs(f) * length(cf) / sum(abs(cf));
            status(abs(f) > threshold) = 2;
        otherwise
            disp('Invalid method chosen. Use "positive", "real" or "abs".');
            return
    end

    % Keep track of survival probability
    survival_prob = zeros(frames, 1);
    total = sum(status == 2);
    cumulative_sp = 0;
    
    % Create folder and prepare snapshots
    mkdir(name + "/" + method);
    snapshot_frames = floor(linspace(1, frames, snapshots));
    snapshotcount = 0;

    % Setting up movie
    fig=figure(2);
    set(gcf,'color','w');
    winsize = [0 0 900 600];
    fig.Position = winsize;
    movmat = moviein(frames,fig,winsize);
    set(fig,'NextPlot','replacechildren');

    function plot_particles()
        % Plot particles
        scatter(ppos(1,status == 1), ppos(2,status == 1), particle_size, status(status == 1), 'filled');
        hold on
        scatter(ppos(1,status == 3), ppos(2,status == 3), particle_size, status(status == 3), 'filled');
        scatter(ppos(1,status == 2), ppos(2,status == 2), particle_size, status(status == 2), 'filled');
        colormap(colors);
        clim([1 3]);

        % Plot contour   
        switch method
            case "positive"
                f = real(cf);
                Z = reshape(f, contour_res, contour_res);
                contour(cX,cY,Z, [0, 0], '-', 'Linewidth', 3, 'Color', 'black');
            case "real"
                f = abs(real(cf)) * length(cf) / sum(abs(real(cf)));
                Z = reshape(f, contour_res, contour_res);
                contour(cX,cY,Z, [threshold, threshold], '-', 'Linewidth', 3, 'Color', 'black');
            case "abs"
                f = abs(cf) * length(cf) / sum(abs(cf));
                Z = reshape(f, contour_res, contour_res);
                contour(cX,cY,Z, [threshold, threshold], '-', 'Linewidth', 3, 'Color', 'black');
        end

        hold off
        xlim([0,1] - 1/(2*particle_res));
        ylim([0,1] - 1/(2*particle_res));
        box on
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
    end
    
    % Simulation loop
    t = 0; % Time
    count = 1; % The number of the next frame
    elapsed = render_time; % Elapsed time since the last rendering
    video_done = false;
    while count <= frames || (cumulative && sum(status == 2) > 0)
        if elapsed >= render_time
            % Update status of particles. To optimize runtime, only an
            % interpolation is used
            cf = sF * exp(2*pi*1i*(sGid2m(:, 1:d_dim)*theta + sGid2m(:, d_dim+1:end)*cpos));
            cf = cf * exp(t*imag(dom_ev)*1i);
            z = reshape(cf, contour_res, contour_res);
            interf = griddedInterpolant(cX,cY,z);
            f = interf(ppos');
    
            switch method
                case "positive"
                    f = real(f);
                    killed = (f < 0);
                case "real"
                    f = real(f) * length(cf) / sum(abs(real(cf)));
                    killed = (abs(f) < threshold);
                case "abs"
                    f = abs(f) * length(cf) / sum(abs(cf));
                    killed = (abs(f) < threshold);
            end

            status(killed(:) & (status(:)==2)) = 3;
            cumulative_sp = cumulative_sp + render_time*sum(status == 2) / total;
        
            if ~video_done
                % The simulation plot
                set(0,'CurrentFigure',fig)
                subplot('Position', [0.05 * 2/3, 0.05, 0.9 * 2/3, 0.9])
                cla();
                plot_particles();
                
                % Survival rate plot
                set(0,'CurrentFigure',fig)
                subplot('Position', [1.05*2/3, 0.55, 0.4*2/3, 0.4]);
                survival_prob(count) = sum(status == 2) / total;
                if survival_prob(count) > 0
                    cla();
                    semilogy([0,T], exp(2*real(dom_ev) * [0, T]), 'black');
                    hold on
                    semilogy([0,T], exp(real(dom_ev) * [0, T]), 'black--');
                    semilogy(linspace(0, t, count), survival_prob(1:count), 'color', colors(2, :));
                    
                    xlim([0, ceil(T)]);
                    ylim([min(exp(2*real(dom_ev)*T), survival_prob(count)), 1]);
                    set(gca, 'YScale', 'log')
                    box on;
                    title('Survival probability');
                end
    
                % Driving plot
                set(0,'CurrentFigure',fig)
                subplot('Position', [1.05*2/3, 0.05, 0.4*2/3, 0.4]);
                cla();
                if d_dim == 2
                    nB = sum(BT < t);
                    bt = [0, BT(1:nB), t];
                    for k=1:nB+1
                        if bt(k+1) - bt(k) > 10^-4
                            p1 = mod(theta_0 + bt(k)*alpha + 10^-5, 1);
                            p2 = mod(theta_0 + bt(k+1)*alpha - 10^-5, 1);
                            line([p1(1), p2(1)], [p1(2), p2(2)]);
                            hold on
                        end
                    end
                elseif d_dim == 1
                    angle = 0:pi/50:2*pi;
                    x = 0.5*cos(angle) + 0.5;
                    y = 0.5*sin(angle) + 0.5;
                    plot(x, y);
                    hold on
                    x = 0.5*cos(2*pi*theta) + 0.5;
                    y = 0.5*sin(2*pi*theta) + 0.5;
                    scatter(x,y)
                end
                hold off
                xlim([0,1]);
                ylim([0,1]);
                box on
                title("Driving position");
                
                movmat(:,count)=getframe(gcf);
    
                % Take snapshot
                if any(snapshot_frames(:) == count)
                    snapfig=figure(3);
                    cla();
                    snapfig.Position = [0 0 900 900];
                    plot_particles();
    
                    filename = name + "/" + method + "/Snapshot" + int2str(snapshotcount) + ".pdf";
                    exportgraphics(snapfig, filename, 'ContentType', 'vector');
                    close;
                    snapshotcount = snapshotcount + 1;
                end
    
                count = count + 1;
    
                if count > frames
                    % create movie file
                    writerObj = VideoWriter(name + "/" + method + "/Simulation",'MPEG-4');
                    writerObj.FrameRate = 24;
                    open(writerObj);
                    writeVideo(writerObj,movmat); 
                    close(writerObj);
                    
                    % Generate final survival probability plot
                    spfig = figure(3);
                    cla();
                    spfig.Position = [0 0 400 300];
                    
                    semilogy([0,T], exp(2*real(dom_ev) * [0, T]), 'black', LineWidth=1);
                    hold on
                    semilogy([0,T], exp(real(dom_ev) * [0, T]), 'black--', LineWidth=1);
                
                    nz = nnz(survival_prob);
                    xspace = linspace(0, T, length(survival_prob));
                    semilogy(xspace(1:nz), survival_prob(1:nz), 'color', colors(2,:), LineWidth=2);
                    hold off
                
                    box on;
                    xlim([0, ceil(T)]);
                    lowLim = min(exp(2*real(dom_ev)*T), survival_prob(nz));
                    ylim([lowLim, 1]);
                    set(gca, 'YScale', 'log')
                    xlabel('Time');
                    ylabel('Survival probability');
                    title('Survival probability');
                    filename = name + "/" + method + "/Survival_Probability.pdf";
                    exportgraphics(spfig, filename, 'ContentType', 'vector');
                    filename = name + "/" + method + "/Survival_Probability.fig";
                    savefig(filename);
                    close
                    video_done = true;
                end
            end

            elapsed = 0;
        end
        
        if video_done
            idx = (status == 2);
            ppos = ppos(:, idx);
            status = status(idx);
            particles = sum(idx);
        end

        % Update particles using Simpson integration
        dat1 = [repmat(theta, 1, particles); ppos];
        k1 = v(dat1);
        dat2 = [repmat(theta+dt/2.*alpha, 1, particles); ppos + dt/2*k1];
        k2 = v(dat2);
        dat3 = [repmat(theta+dt.*alpha, 1, particles); ppos - dt*k1 + 2*dt*k2];
        k3 = v(dat3);
        ppos = ppos + dt*(k1/6 + 4*k2/6 + k3/6);
    
        diff = sqrt(dt) * randn(2, particles);
        ppos = ppos + epsilon * diff;
        ppos = mod(ppos, 1);

        % Update driving 
        theta = mod(theta + alpha * dt, 1);

        % Update time
        t = t + dt;
        elapsed = elapsed + dt;
    end
    
    if cumulative
        fprintf("Cumulative survival probability: %f5.\n", cumulative_sp);
    end
end