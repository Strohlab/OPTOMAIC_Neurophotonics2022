function [DC_VisStim_signal,DC_VisStim_onsets,DC_Carec_signal,DC_Carec_OnOffsets,DC_sync_delay_s,err]=read_Spike2(g_MRI_cutoff_stst_ms,g_DC_Spike2_SamplingRate,datum,prot_fid,analyse_animalfolder,animal_ident,DC_ProtImages);
 
s=['- \n']; fprintf(prot_fid,s);disp(s)
analyse_function_name_version='FUNCTION DC_read_Spike2 V20220519';
analyse_function_author='Dirk Cleppien';
s=['### (' analyse_function_name_version ' - ' analyse_function_author ')  \n']; fprintf(prot_fid,s);disp(s)

%% Function to preprocess SPIKE2 data

%% Load SPIKE2 data
[fname,pname]=uigetfile('*.mat','Select the matfile of SPIKE2 to be transformed',['../data/' analyse_animalfolder] );
old_folder=cd(pname);
file_fields=fields(load(fname));
myfile=load(fname);
%%
%% Load Spike2 data
%% Data format of Spike2 file:
%   MR gradient signal:             gradient    --> gradient
%   Calcium transient:              Carec1      --> carec1
%   stimulation pulse time-course:  trigger     --> trigger
%   breathing                       breathin    --> breath
%   blood saturation O2             spo2        --> spo2
spo2=[]

disp( file_fields)

for rr=1:length(file_fields) %loop on the different channels

    eval(['myname=myfile.',file_fields{rr},';']);
    myname.title

    if (strcmp(myname.title,'CaRec1')),
        carec1=myname.values;
    elseif (strcmp(myname.title,'gradient')),
        gradient=myname.values;
    elseif (strcmp(myname.title,'visual')),
        visual=myname.values;
    elseif (strcmp(myname.title,'fiber2')),
        trigger=myname.values;
    elseif (strcmp(myname.title,'spo2')),
        spo2=myname.values;
    elseif (strcmp(myname.title(1:4),'brea')),
        breath=myname.values;
    end
end
size(carec1)
size(gradient)
size(breath)
size(spo2)

%% Depiction of signal-time courses
%% Initialization of analysis figure
fig=figure('Name',s,'units','normalized','outerposition',[0 0 1 1]);
fig.ToolBar='none';
fig.NumberTitle='off';
fig.Name=[animal_ident ' - BOLD response analysis'];
fig.FileName=[animal_ident '-BOLDResponseAnalysis'];

fig_row=4;
fig_col=5;
fig_counter=1;
%% Depiction of measured Spike2 time courses
subplot(fig_row,fig_col,fig_counter)
    plot(gradient)
    title(['Fig.' num2str(fig_counter) ': Measured MR gradient signal']);
    ylabel('arbitrary units')
    xlabel('time')
    fig_counter=fig_counter+1;
gradient_max=max(gradient)
gradient_std1000=std(gradient(1:1000))
gradient_std100000=std(gradient(50000:100000))
subplot(fig_row,fig_col,fig_counter)
    plot(gradient(1:100000)./(max(gradient)),'b')
    hold on 
    plot(trigger(1:100000)./max(trigger),'r')
    hold off
    title(['Fig.' num2str(fig_counter) ': -"- (zoomed  to start)']);
    ylabel('arbitrary units')
    xlabel('time')
    fig_counter=fig_counter+1;
subplot(fig_row,fig_col,fig_counter)
    plot(carec1);
    title(['Fig.' num2str(fig_counter) ': Measured Calcium signal'])
    ylabel('arbitrary units')
    xlabel('time')
    fig_counter=fig_counter+1;
subplot(fig_row,fig_col,fig_counter)
    plot(carec1(1:100000)./max(carec1))
    title(['Fig.' num2str(fig_counter) ': -"- (zoomed to start)']);
    hold on 
    plot(trigger(1:100000)./max(trigger),'r')
    hold off
    ylabel('arbitrary units')
    xlabel('time')
    fig_counter=fig_counter+1;
subplot(fig_row,fig_col,fig_counter)
    plot(trigger)
    title(['Fig.' num2str(fig_counter) ': Measured Stimulation pulses'])
    ylabel('arbitrary units')
    xlabel('time')
    fig_counter=fig_counter+1;

%% Extraction of start point of MR pulse sequence
step=100;
max_step=100000
intervall_std_before=max(gradient)
counter_mri=0;
mri_start=0;
intervall_min=0;
%%
for zi=1:step:max_step,
    zi
    intervall_std=std(gradient(zi:zi+step))
    if (intervall_std>2*intervall_std_before),
        counter_mri=counter_mri+1
        if (counter_mri==2),
            mri_start=zi
            dy=gradient(zi:zi+step-1)<intervall_min*2
            dx=1:step;
            dx=dx(dy);
            min(dx)
            mri_start=min(dx)+zi;
            %keyboard
        end
    end
    intervall_std_before=intervall_std
    intervall_min=min(gradient(zi:zi+step))
end
mri_start   % timestamp of start of mri sequence

%% Extract all measured trace interval simultaneous to MRI measurement
% Cut out measurement interval before starting Ca-measurement 
% from trigger and carec1
trigger1=trigger;
trigger(1:mri_start)=0;

%% 
dx=1:length(trigger);
trigger_start_1=min(dx(trigger>3)) % first stimulation pulse (signal>3[V])
%  !!!! Shift by g_cutoff_stst - this is the time the MR signal needs to
%  get into steady state, actually one stimulation period
trigger_start=trigger_start_1+g_MRI_cutoff_stst_ms;
trigger_end=max(dx(trigger>3)) % last stimulation pulse

%%
carec1_start=min(dx(carec1<0))+10
messung_start=trigger_start;
messung_end=trigger_end;
carec_messung=carec1(messung_start:messung_end);
DC_sync_delay_s=(messung_start-mri_start-g_MRI_cutoff_stst_ms)/g_DC_Spike2_SamplingRate 
breath_messung=breath(messung_start:messung_end);

%% Depiction of extraction
subplot(fig_row,fig_col,fig_counter)
    plot(carec1(trigger_start:trigger_end))
    title(['Fig.' num2str(fig_counter) ': Cutted Calcium signal'])
    ylabel('arbitrary units')
    xlabel('time')
    fig_counter=fig_counter+1;
subplot(fig_row,fig_col,fig_counter)
    plot(trigger(messung_start:messung_end)./max(trigger))
    hold on
    plot(carec1(messung_start:messung_end)./max(carec1(messung_start:messung_end)),'g');
    hold off
    title(['Fig.' num2str(fig_counter) ': Combined Stim- and Calcium time course'])
    ylabel('arbitrary units')
    xlabel('time')
    fig_counter=fig_counter+1;
subplot(fig_row,fig_col,fig_counter)
    plot(gradient(mri_start:messung_end )./(max(gradient)))
    hold on 
    plot(trigger(mri_start:messung_end )./max(trigger),'r')
    plot(abs(carec1(mri_start:messung_end)./max(abs(carec1(messung_start:messung_end)))/200),'g')
    hold off
    title({['Fig.' num2str(fig_counter) ': Combined depiction of sync. signals'];['Offset(MR-Stim) =' num2str(DC_sync_delay_s) 's']})
    ylabel('arbitrary units')
    xlabel('time')
    fig_counter=fig_counter+1;
subplot(fig_row,fig_col,fig_counter)
    plot(trigger(mri_start:messung_end )./max(trigger),'r')
    plot(abs(carec1(mri_start:messung_end)./max(abs(carec1(messung_start:messung_end)))/200),'g')
    hold off
    title({['Fig.' num2str(fig_counter) ': Combined depiction of sync. signals'];['Offset(MR-Stim) =' num2str(DC_sync_delay_s) 's']})
    ylabel('arbitrary units')
    xlabel('time')
    fig_counter=fig_counter+1;% DC 1.9.2020 Depiction of breathing
subplot(fig_row,fig_col,fig_counter)
    xb=(mri_start:messung_end)-mri_start+1;
    b_data=breath(mri_start:messung_end )./ceil(max(breath));
    plot(xb(1:min(length(b_data),40000)),b_data(1:min(length(b_data),40000)));
    title({['Fig.' num2str(fig_counter) ': Breathing Curve'];['Zoomed to start']})
    fix_fig=fig_counter;
    fig_counter=fig_counter+1;% DC 1.9.2020 Depiction of breathing
%% Calculation of Slow wave onsets
cd (old_folder)

%% 
loop_opt=1;
while (loop_opt),
    fig_counter_loop=fig_counter;
    %% Parameter 
    prompt2 = {'min length (ms)','min refrac (ms)','min highA (%std)','min highB (%std)',...
        'addon (ms)','sampl freq (kHz)','Window Base (ms)','Window ExpAvg (ms)','downsampling'};
    dlg_title2 ='Variables for identifying onsets:(if u dont know DONT TOUCH)';
    num_lines2 = 1;
    defaultans2 = {'300','300','30','10','0','2','3000','25','0'};
    answer2 = inputdlg(prompt2,dlg_title2,[1, length(dlg_title2)+30],defaultans2);
    min_length=str2double(answer2{1,1}); %(event in ms)
    min_refrac=str2double(answer2{2,1}); %(event in ms)
    min_highA=str2double(answer2{3,1});
    min_highB=str2double(answer2{4,1});
    addon=str2double(answer2{5,1});
    sampl_freq=str2double(answer2{6,1}); %(khz)
    Window_Base=str2double(answer2{7,1});
    Window_ExpAvg=str2double(answer2{8,1});
    downsampling=str2double(answer2{9,1});
    %trsnform above from ms to samplingrate
    min_length=min_length*sampl_freq;
    min_refrac=min_refrac*sampl_freq;
    Window_Base=Window_Base*sampl_freq;
    Window_ExpAvg=Window_ExpAvg*sampl_freq;
    addon=addon*sampl_freq;

    %% written in analysis figure
    subplot(fig_row,fig_col,fig_counter_loop)
    title('Variables for identifying onsets');
    text(0.01,5,prompt2);
    text(0.7,5,answer2);
    ylim([0 10])
    xlim([0 1])
    fig_counter_loop=fig_counter_loop+1;

     %% flipping of calcium signal for depiction
     data=-carec_messung*2000;

    %% Deleting of one small spike due to some weird effects (no influence on detection)
    Signal1=data;
    std_Signal1=std(Signal1);
    mean_Signal1=mean(Signal1);
    thres_spike=mean_Signal1+5*std_Signal1;
    d=(Signal1<thres_spike);
    [d_y,d_min]=min(d);
    d_Signal1=Signal1;
    d_Signal1(d_min:d_min+10)=mean_Signal1;

    %% depiction in analysis figure
    subplot(fig_row,fig_col,fig_counter_loop)
    plot(Signal1)
    hold on
    plot(d_Signal1)
    hold off
    Signal1=d_Signal1;
    title(['Fig.' num2str(fig_counter_loop) ': Flipped Calcium signal'])
    fig_counter_loop=fig_counter_loop+1;

    %%
    ant=input('Remove a Spike? [j/n]','s');
    if (ant=='j'),
        data=Signal1;
    end

    %% Downsampling if needed
    if downsampling ~= 0
        for i=1:downsampling
            if mod(length(data),2) ~= 0
                data(1,:)=[];
            end
            data=(DC_dyaddown(data,0)+DC_dyaddown(data,1))/2;
            min_length=min_length/2;
            min_refrac=min_refrac/2;
            Window_Base=Window_Base/2;
            Window_ExpAvg=Window_ExpAvg/2;
        end
    end
    Signal_raw=data;
    % Round the "downsampled" values to whole numbers
    min_length=round(min_length);
    min_refrac=round(min_refrac);
    Window_Base=round(Window_Base);
    Window_ExpAvg=round(Window_ExpAvg);
    addon=round(addon);

    %% Flatten the baseline
    % Signal1 is processed calcium signal
    Signal1=DC_Baseline(Signal_raw,Window_Base,Window_Base/2);
    
    %% Depiction of downsampling and flattening
    hold on
    plot(Signal_raw,'b')
    plot(Signal1+min(Signal_raw(:)),'g')
    title(['Fig.' num2str(fig_counter_loop) ': Flipped+flattened calcium signal'])
    ylabel('arbitrary units')
    xlabel('time')
    hold off

    %% Onset calculation
    % now the onsets of carec
    % Moving average with exponential weighting (s. Seamari, PlosOne, 2007)
    thresh_start1=Signal1;
    thresh_start1(:)=0; % will contain the threshold to define the onset of an up-state
    box1=Signal1; %binary vector, 0= down, 1 = up, preliminary
    box1(:)=0;
    onset1=box1; % binary vector containing only the onsets of up-states
    box2=box1;
    onset2=box2;
    
    %%
    [n,xout] = hist(Signal1(2500:numel(Signal1)-2500),ceil((numel(Signal1)-5000)/1000));
    sum_n=cumsum(n);

    off_hist=find(sum_n >= max(sum_n)*(1-min_highA/100),1,'first');
    min_high1=xout(off_hist); % new definition of min_high1 based on the calculated results
    off_hist=find(sum_n >= max(sum_n)*(1-min_highB/100),1,'first');
    min_high2=xout(off_hist); % new definition of min_high2 based on the calculated results
    
    %Moving average with exponential weighting (s. Seamari, PlosOne, 2007)
    thresh_start1(:)=min_high1; %min_high1
    Signal_filt = DC_tsmovavg(Signal1, 'e', Window_ExpAvg, 1);
    
    %%
    subplot(fig_row,fig_col,fig_counter_loop)
    plot(Signal1,'b')
    hold on
    plot(Signal_filt,'g')
    hold off
    title(['Fig.' num2str(fig_counter) ': Filtered Calcium Signal'])
    ylabel('arbitrary units')
    xlabel('time')
    fig_counter=fig_counter+1;
    
    %% box finding
    % calculate a preliminary on-off-vactor
    box1(Signal_filt >= thresh_start1)=0.001;
    %%
    subplot(fig_row,fig_col,fig_counter_loop)
    plot(box1/0.001*500)
    hold on
    plot(Signal1*50,'g')
    %%
    % remove too short gaps between upstates SIGNAL1
    % here the variable min_refrac is used
    count=0;
    for i=2:numel(Signal1)
        %start to count when upstate ends
        if box1(i)==0
            count=count+1;
        end
        %if the gap is shorter than min_refrac
        %set the gap to 1
        if (count < min_refrac) && (box1(i)==0.001) && (box1(i-1)==0)
            box1(i-count:i-1)=0.001;
            count=0;
        end
        %if the gap is longer than min_refrac
        %reset the counter
        if (count >= min_refrac) && (box1(i)==0.001) && (box1(i-1)==0)
            count=0;
        end
    end
    %
    plot(box1/0.001*600,'g')
    j=1;
    % correct offsets SIGNAL1
    for z=2:numel(Signal1)
        %search for an offset
        if box1(z)==0 && box1(z-1)==0.001
            j=z;
            %keep the vector on 1 as long as the signal is high enough
            if j+z <= numel(Signal1) && Signal1(j+1) >= min_high1*0.15
                j=j+1;
                box1(z)=0.001;
            end
        end
        z=j+1;
    end
    %
    plot(box1/0.001*700,'r')
    % remove suspicious upstates Signal1
    count=0;
    for u=2:numel(Signal1)
        %start counting when an upstate starts
        if box1(u)==0.001
            count=count+1;
        end
        %if the upstate is too short (shorter than min_length)
        %set it to 0
        if (box1(u)==0) && (box1(u-1)>0) && (count < min_length)
            box1(u-count:u-1)=0;
            count=0;
        end
        %if the upstate is long enough, test if it is high enough
        %(higher than min_high2)
        %otherwise set it to 0
        if (box1(u)==0) && (box1(u-1)>0) && (count > min_length)
            if max(Signal1(u-count:u-1)) < min_high2
                box1(u-count:u-1)=0;
            end
            count=0;
        end
    end
    %
    plot(box1/0.001*800,'k')
    xlim([200000 400000])
    s_opt=input('Slow wave detection ok [y/n]?','s')
    if (s_opt=='y'),loop_opt=0;end
    
end
fig_counter=fig_counter_loop+1;

%% Collection of ONSETS AND OFFSETS
onsets=[];
offsets=[];

%%
for i=2:length(box1)
    % find an onset when "box" switches from 0 to 1   
    if box1(i) == 0.001 && box1(i-1) == 0
           onsets=vertcat(onsets,i);
    end
end
%
for i=1:(length(box1)-1)
    % find an onset when "box" switches from 0 to 1   
    if box1(i) == 0.001 && box1(i+1) == 0
           offsets=vertcat(offsets,i); 
    end
end
%
if (size(onsets)~=size(offsets)), disp('Do not fit!'); return;end
size(onsets)

%% Analysing VisStim signal
% finding of stimulation intervals
x=1:length(Signal1);
s_version=version;
    if (s_version(1:19)=='7.12.0.635 (R2011a)'),
    [trig_pk,trig_lk]=findpeaks(trigger(messung_start:messung_end),'MinPeakDistance',2000,'MinPeakHeight',0);
else
    [trig_pk,trig_lk]=findpeaks(trigger(messung_start:messung_end),x,'MinPeakDistance',2000,'MinPeakHeight',0);
end%
disp(['Number of stimulation intervals: ' num2str(size(trig_lk,1))]);

%% ########
%% Slow Waves within response window 30-300ms
% Attention! Sampling rate is 2kHz, 2points per ms!
response=[];
zj=1;
for zi=1:size(onsets),
    %onsets(zi)
    %trig_lk(zj+1)
    % not all vis stim has a detected response!
    while (onsets(zi)>trig_lk(zj+1)), zj=zj+1; disp(['While ' num2str(zj)]); end
    
    %zj
    %keyboard
    if (zj<length(trig_lk)),
        dummy=onsets(zi)-trig_lk(zj);
        response=[response dummy];
    end
    disp([num2str(zi) ' ' num2str(zj) ' ' num2str(onsets(zi)) ' ' num2str(offsets(zi)) ' ' num2str(trig_lk(zj)) ' ' num2str(trig_lk(zj+1)) ' ' num2str(dummy)]);
    
end
response_interval_min=30*2; % 30ms
response_interval_max=300*2; % 300ms

%%
[response_y,response_x]=hist(response,max(response)/g_DC_Spike2_SamplingRate*10);
%%
subplot(fig_row,fig_col,fig_counter)
    bar((response_x/g_DC_Spike2_SamplingRate),response_y/sum(response_y))
    title(['Fig.' num2str(fig_counter) ': Histogram of SlWaves response time'])
    ylabel('arbitrary units')
    xlabel('time [s]')
    ylim([0 0.5])
    fig_counter=fig_counter+1;


%%
d_response=response.*(response<response_interval_max)
d1_response=d_response.*(d_response>response_interval_min)

resp_events=sum(d1_response>0)

%% Function tail
% On- Offsets in sec
DC_VisStim_signal=trigger(messung_start:messung_end);
    size(DC_VisStim_signal)
DC_VisStim_onsets=trig_lk;
    size(DC_VisStim_onsets)
DC_Carec_signal=Signal1;
    size(DC_Carec_signal)
DC_Carec_OnOffsets=vertcat(onsets',offsets');
    size(DC_Carec_OnOffsets)
DC_sync_delay_s
err=0;

%% for visualization filtering of Carec-signal
size(Signal1)
s1=1;
s2=length(Signal1);
Medfilter_length=60;
DC_data_Medfilter=medfilt1(Signal1,Medfilter_length);
%%
subplot(fig_row,fig_col,[16:20])
hold on
offset=2.2
m_Signal1=mean(Signal1)
m_Filter=mean(DC_data_Medfilter)
sDC_data_Medfilter=DC_data_Medfilter;
plot(sDC_data_Medfilter(s1:s2)-m_Filter+offset,'k')
x=[0 s2-s1]
y=zeros(size(x,2),1)+offset
stim_pulse=(DC_VisStim_onsets>s1 & DC_VisStim_onsets<s2).*DC_VisStim_onsets
ystim=zeros(sum(stim_pulse>0),1)+0.8
stim_pulse(stim_pulse>0)

bar(stim_pulse(stim_pulse>0)-s1,ystim,0.02,'k')
ylim([-3 10])
xlim([1 s2+100])
sig=zeros(length(onsets),1)+1;
    plot(trig_lk,trig_pk-2.9+1.6*offset,'g+')
    plot(onsets,sig+1.5*offset,'rx')
   plot(xb,6*(b_data-1.8*mean(b_data)))
title(['Fig.' num2str(fig_counter) ': Overview Slow Wave detection'])
hold off
%%

image_suffix='Subroutine_Analysis_SPIKE2-file';
old=cd(['../data/' analyse_animalfolder '/']);
    DC_hg_name=[ animal_ident '_' datum '_' image_suffix '.emf' ];
saveas(fig,DC_hg_name,'emf')
cd(old)

%% - End of function - 
s=['### (' analyse_function_name_version ') - end \n']; fprintf(prot_fid,s);disp(s)

