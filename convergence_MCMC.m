
for setting = 1:3
    switch setting
        case 1
            dataName = '_Per';
            figString = 'UGL';
        case 2
            dataName = '_simple_model';
            figString = '4DIR';
        case 3
            dataName = '_better_simple_model';
            figString = 'ANYDIR';     
    end
    
    dataDir = fullfile('Results-onlytop','sub-1', ['MCMC2D' dataName], 'everything.mat');
    load(dataDir, 'alphaVec', 'betaVec')
    
    figure(1)
    subplot(2,3,setting)
    
    plot(alphaVec', 'LineWidth', 2)
    title(['Convergence of \alpha_k with ', figString, ' method'])
    
    subplot(2,3,setting+3)
  
    plot(betaVec')
    title(['Convergence of \beta_p with ', figString, ' method'])
    
end

%%
set(gcf, 'Position', [0,0,2000,1200])
