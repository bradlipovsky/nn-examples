clear;

RunTheseSteps = [1 2 3];

%
% Run the ice sheet model once and plot the solution.
%
if any(RunTheseSteps == 1)
    p = LoadParameters;
    Hinit = p.C*p.uleft / (p.rho*p.g*p.alpha);
    H0 = Hinit * (ones(p.J+1,1) -  (1:p.J+1)' * p.alpha);
    u0 = ssaflowline(p,H0);

    figure(1);
    subplot(211); plot(p.x/1e3,H0/1e3); 
    hold on; ylabel('Thickness (km)'); xlabel('Distance (km)');

    subplot(212); plot(p.x/1e3,u0); 
    hold on; ylabel('Velocity (m/y)'); xlabel('Distance (km)');
end

%
% Now run the ice sheet model many times with many different perturbations.
%
if any(RunTheseSteps == 2)
    nit = 1000;
    u = zeros(p.J+1,nit); H = u;
    tic; figure(2);
    for i = 1:nit
        location = rand*p.L;
        width = rand*(p.L-p.L/p.J* 2) + p.L/p.J* 2;
        height = rand*100;

        H(:,i) = H0 + height*exp(-(p.x-location).^2 / width^2 );
        u(:,i) = ssaflowline(p,H(:,i))';

        disp(i); toc
    end


    subplot(211); plot(p.x/1e3,H/1e3); hold on;
    ylabel('Thickness (km)'); xlabel('Distance (km)');

    subplot(212); plot(p.x/1e3,u); 
    hold on; ylabel('Velocity (m/y)'); xlabel('Distance (km)');

    drawnow;

    dU = (mean(u) - mean(u0))/mean(u0) ;
    Hscaled = (H-min(H(:)))/(max(H(:))-min(H(:)));

    save('ISM_Output.mat','dU','Hscaled');
end

%
% Train neural network on the output of the ice sheet model.
%
if any(RunTheseSteps == 3)
    load ISM_Output
    tic;
    net = feedforwardnet([10 10]);

    [net,tr] = train(net,Hscaled,dU);

    err = net(Hscaled) - dU;
    Residual = mean( abs(err) );
    figure;
    subplot(211);histogram(log10(abs(err))); xlabel('Log10 |Error|');
    title(['Mean Scaled Error = ' num2str(Residual*100) ' (Percent Speedup)']);
    subplot(212);histogram(dU*100,10); xlabel('Percent Speedup, dU');xlim(100*[-1 1]);
    toc
end