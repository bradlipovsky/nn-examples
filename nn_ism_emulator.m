clear;

RunTheseSteps = [3];

p = LoadParameters;
%
% Run the ice sheet model once and plot the solution.
%
if any(RunTheseSteps == 1)
    
    u0 = ssaflowline(p,p.H0);

    figure(1);
    subplot(211); plot(p.x/1e3,p.H0/1e3); 
    hold on; ylabel('Thickness (km)'); xlabel('Distance (km)');

    subplot(212); plot(p.x/1e3,u0); 
    hold on; ylabel('Velocity (m/y)'); xlabel('Distance (km)');
    save('u0','u0');
end

%
% Now run the ice sheet model many times with many different perturbations.
%
if any(RunTheseSteps == 2)
    tic;
    load u0
    
    NumberOfSimulations = 100;
    u = zeros(p.J+1,NumberOfSimulations); H = u;
    tic; figure(2);
    for i = 1:NumberOfSimulations
        location = rand*p.L;
        width = rand*(p.L-p.L/p.J* 2) + p.L/p.J* 2;
        height = rand*100;

        H(:,i) = p.H0 + height*exp(-(p.x-location).^2 / width^2 );
        u(:,i) = ssaflowline(p,H(:,i))';
    end
    toc_sim = toc/NumberOfSimulations;

    dU = (mean(u) - mean(u0))/mean(u0) ;
    Hscaled = (H-min(H(:)))/(max(H(:))-min(H(:)));

    save('ISM_Output.mat','dU','Hscaled','toc_sim');
    
    figure;
    subplot(211); plot(p.x/1e3,H/1e3); hold on;
    ylabel('Thickness (km)'); xlabel('Distance (km)');
    subplot(212); plot(p.x/1e3,u); 
    hold on; ylabel('Velocity (m/y)'); xlabel('Distance (km)');
    
    figure;
    histogram(dU*100,10); xlabel('Percent Glacier Vel. Change, dU');xlim(100*[-1 1]);
end

%
% Train neural network on the output of the ice sheet model.
%
if any(RunTheseSteps == 3)
    load ISM_Output
    load u0
    
    net = feedforwardnet([5 5 5]);
    [net,tr] = train(net,Hscaled,dU);
    
    tic;
    err = net(Hscaled) - dU;
    toc_nn = toc / numel(dU);
    Residual = mean( abs(err) );
    
    figure;
    histogram(log10(abs(err))); xlabel('Log10 |Error|');
    title(['Mean Scaled Error = ' num2str(Residual*mean(u0)) ' (Glacier Vel., m/a)']);
    hold on; line(log10([3e-3 3e-3]),ylim,'color','r');
    legend('Error Histogram','Simulation Accuracy');
    
    disp('=================================================');
    disp(['Total time to run simulations:  ' num2str(toc_sim)]);
    disp(['Total time to run NN:  ' num2str(toc_nn)]);
    disp(['Speed up by a factor of ' num2str(round(toc_sim/toc_nn)) ]);
end