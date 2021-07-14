function draw_raster(rgc, start_trial_n, end_trial_n) 


% Modify expression to add input arguments.
% Example:
%   draw_raster(RCG ID, Initial Trial, last trial);

load spikes_onset
figure

spike_train_n =size(spikes_onset,1);
spike_n = zeros(spike_train_n,1);


for i=1:total_stim_n-1
    spike_n(i) = size(spikes_onset{rgc}{i},1);
end


spike_train_n = size(spikes_onset,1);  %including all cells 

f_size=15; %font size



for i=start_trial_n:end_trial_n
    for j=1:spike_n(i)
        
        line([spikes_onset{rgc}{i}(j) spikes_onset{rgc}{i}(j)],[i-.5 i+.5],'Color','k')
    
    end
end


ylim([start_trial_n-1 end_trial_n+1])
xlim([0 7])
xlabel('time(s)')
set(gca,'FontSize',f_size)
ylabel('trials')
set(gcf,'color','w')
title(['RGC #' num2str(rgc)]);

    
end



