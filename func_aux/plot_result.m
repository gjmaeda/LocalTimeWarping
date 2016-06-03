function [] = plot_result( dataRef, dataBefore, dataAfter, k )

    hd = figurew(['aligned_data' num2str(k)]);
    set_fig_position([0.331 0.106 0.195 0.792]);
    subplot(3,1,1); grid on; hold on; ylabel 'x';
    subplot(3,1,2); grid on; hold on; ylabel 'y';
    subplot(3,1,3); grid on; hold on; ylabel 'z'; xlabel 'time steps';

    for j=1:3
        subplot(3,1,j);
        plot(linspace( 0,1,numel( dataRef(:,j) )),   dataRef(:,j), sty([0 0 1], [], 3)  );
        plot(linspace( 0,1,numel( dataBefore(:,j) )),   dataBefore(:,j), sty([0.65 0.65 0.65], [], 3)  );
        plot(linspace( 0,1,numel( dataAfter(:,j) )),   dataAfter(:,j), sty([1 0 0 ], [], 3)  );
    end
    subplot(3,1,1);
    legend({'Reference', 'Before', 'After'});
    
    
    
    
end