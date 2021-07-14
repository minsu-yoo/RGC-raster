clear

load('Data.mat')
%it load the data from sorted spikes

channel_labels = who;


nChannels = size(channel_labels,1);  %total number of channels


spikeCh_n = nChannels-1;

last_ch = channel_labels(nChannels);

Spikes = cell(nChannels,3);
str_ch = cell(nChannels,1);



chID = last_ch{1};
chID(4:end-1);
nCh = str2num(chID(4:end-1));


RGCid = 1;

%%

RGCid = 0;

chID = channel_labels{2};
firstCh = str2num(chID(4:end-1));

current_ch = firstCh;

dibuk = cell(11,4); %to debug
chan_neuron_id = zeros(nChannels-1,2);


for i= 2:nChannels
    
    chID = channel_labels{i};
    nCh = str2num(chID(4:end-1));
    
    %next_nCh = str2num(next_chID(4:end-1));
    %next_chID = channel_labels{i+1};
    %current_ch = firstCh;
    
    %to_next_ch = next_nCh - nCh;
    
    if (current_ch==nCh)
        RGCid = RGCid+1;
        Spikes{nCh, RGCid} = eval(chID);
    else
        
        to_next_ch = nCh - current_ch;
        %channel should go up 
        
        current_ch = current_ch + to_next_ch;
        RGCid = 1;
        Spikes{nCh, RGCid} = eval(chID);
        
    end
    dibuk(i,:) = {chID,current_ch, nCh,RGCid};
    chan_neuron_id(i-1,:) =  [nCh,RGCid];

      
 
end

%%

chN = 12;
neuron_id = 1;
f_size = 15;

spike_n = size(Spikes{chN,neuron_id},1);
firings = Spikes{chN,neuron_id};



for j=1:100
    
line([firings(j) firings(j)],[1-.5 1+.5],'Color','k')
hold on

end

ylim([0 2])
%xlim([0 0.6])
xlabel('time(s)')
 set(gca,'FontSize',f_size)
ylabel('trials')
 set(gcf,'color','w')

 
%%
  close all
spikeCh_n = nChannels-1;

total_Ch_N = size(Spikes,1);
total_cell_N = size(Spikes,2);

stim_time = A1a;
spikes_in_order = Spikes(~cellfun('isempty',Spikes));
sorted_Spikes = cell(spikeCh_n,1);
count = 1;


for i = 1:total_Ch_N
    for j=1:total_cell_N
    
    if ~isempty(Spikes{i,j})
        sorted_Spikes{count} = Spikes{i,j};
        count = count+1 ;
    else
    end
    
    end
    
end


%%  segrate RGC responses by onsets of visual stimulus


spike_train_n = size(sorted_Spikes,1);  %including all cells 
total_stim_n = size(stim_time,1);  % total number of onsets
fire_in_cell = cell(total_stim_n,1); % spikes are saved as the data strucutre of 
%cell 

count = 1; 
spikes_by_stim = cell(spike_train_n,1);

for i=1:spike_train_n
    spike = sorted_Spikes{i};
    first_trial_idx = find(spike > stim_time(1));
    
    
    for j = 2:total_stim_n
        
        fire_in_cell{1} = spike(first_trial_idx);
        idx = find (spike>stim_time(j-1)& spike<stim_time(j));
        fire_in_cell{j} = spike(idx);
    end
    
     spikes_by_stim{i} = fire_in_cell;
    
end


%%

spikes_onset = cell(spike_train_n,1); % spikes are saved as the data strucutre of 
   
for i=1:spike_train_n
    
    
    for j = 2:total_stim_n
        
        %spikes_onset{i}{1} = spikes_by_stim{i}{1};
        spikes_onset{i}{j-1} = spikes_by_stim{i}{j} - stim_time(j-1);
        
        
    end
    
end


save spikes_onset

        
%% draw a neural firing based on cell number


f_size=15; %font size



spike_n = zeros(spike_train_n,1);


for i=1:total_stim_n-1
spike_n(i) = size(spikes_onset{1}{i},1);
end


trial_n = 20;

for i=1:trial_n
for j=1:spike_n(i)
line([spikes_onset{1}{i}(j) spikes_onset{1}{i}(j)],[i-.5 i+.5],'Color','k')
hold on
end
end

ylim([0 trial_n+1])
xlim([0 4])
xlabel('time(s)')
set(gca,'FontSize',f_size)
ylabel('trials')
set(gcf,'color','w')

    



