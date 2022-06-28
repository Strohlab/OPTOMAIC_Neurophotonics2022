function [DC_data_filter,DC_ProtImages, err_log]=DC_filter_data(DC_data_long,Rep_Stim,MR_SamplingRate,datum,DC_ProtImages,err_log,prot_fid,analyse_animalfolder,animal_ident);
%% 
s='- \n'; fprintf(prot_fid,s);disp(s)
analyse_function_name_version='FUNCTION DC_filter_data; V20200402_breathingV20201104';
analyse_function_author='Dirk Cleppien';
s=['### (' analyse_function_name_version ' - ' analyse_function_author ')  \n']; fprintf(prot_fid,s);disp(s)

%%
%% Update 9.12.2020
%% filters
%% 1. raw
%% 2. baseline correction
%% 3. Filter 1
%% 4. Median filter
%% Depiction: S-t-course per voxel / spectral analysis per filter step
err_log=1;

%% Filter 
s=['### Avaiable filters:  \n']; fprintf(prot_fid,s);disp(s)
s=['### Baseline detrending (get rid of heating effects)  \n']; fprintf(prot_fid,s);disp(s)
s=['### LowPass-Filter  \n']; fprintf(prot_fid,s);disp(s)
s=['### Median-Filter in time:  \n']; fprintf(prot_fid,s);disp(s)
%
num_filter=4;s=['### Var: Number available filtered data num_filter =  ' num2str(num_filter) ' \n']; fprintf(prot_fid,s);disp(s)
s=['### Raw data size = ' num2str(size(DC_data_long)) '  \n']; fprintf(prot_fid,s);disp(s)
size_DC_data=size(DC_data_long);s=['### Var: size_DC_data =  ' num2str(size_DC_data) ' \n']; fprintf(prot_fid,s);disp(s)

% data arrays for calculation
DC_data_filter_baseline=DC_data_long;s=['### Baseline correction, size = ' num2str(size(DC_data_filter_baseline)) '  \n']; fprintf(prot_fid,s);disp(s)
DC_data_filter_FIR=DC_data_long;s=['### Low Pass Filter size = ' num2str(size(DC_data_long)) '  \n']; fprintf(prot_fid,s);disp(s)
DC_data_filter_median=DC_data_long;s=['### Median Filter size = ' num2str(size(DC_data_long)) '  \n']; fprintf(prot_fid,s);disp(s)
DC_data_filter_HPass=DC_data_long;s=['### Median Filter size = ' num2str(size(DC_data_long)) '  \n']; fprintf(prot_fid,s);disp(s)

% returned data arr: DC_data_filter
DC_data_filter=zeros(size(DC_data_long,1),size(DC_data_long,2),num_filter);s=['### Var: size_DC_data_filter =  ' num2str(size(DC_data_filter)) ' \n']; fprintf(prot_fid,s);disp(s)

% the voxel for depiction: pixel with max signal
Signal_profile=mean(mean(DC_data_long(:,:),3),2);
[py,px]=max(Signal_profile);
depicted_pixel=px;s=['### local Var: depicted_pixel = ' num2str(depicted_pixel) ' \n']; fprintf(prot_fid,s);disp(s)
depicted_timepoint=600;s=['### local Var: depicted_timepoint = ' num2str(depicted_timepoint) ' \n']; fprintf(prot_fid,s);disp(s)
Medfilter_length=20;s=['### local Var Medfilter_length: ' num2str(Medfilter_length) ' \n']; fprintf(prot_fid,s);disp(s)
xlim_spectrum=[0 2];

%% perform filtering and store it in DC_data_filter (:,:,n); n depending on kind of filter
% DC_data_filter: n=1 -> raw data
DC_data_filter(:,:,1)=DC_data_long;s=['### DC_data_filter (:,:,1)= raw data [DC_data_long]  \n']; fprintf(prot_fid,s);disp(s)

%% Spectral Analysis
Fs=MR_SamplingRate %sampling per second
NFFT=length(DC_data_filter(depicted_pixel,:,1));
F2fft=(0:1/NFFT:1/2-1/NFFT)*Fs;

%% Analysis figure
fig=figure('units','normalized','outerposition',[0 0 1 1]);
fig_row=num_filter+1;

fig_col=5;
fig.ToolBar='none';
fig.NumberTitle='off';
fig.Name=['Subroutine FilterData'];
fig.FileName=['Subroutine FilterData'];
fig_counter=1;
figrow_counter=1;

%% subplot shows unfiltered and filtered data of depicted_pixel for comparison
subplot(fig_row,fig_col,fig_counter); %#1
    plot(DC_data_filter(depicted_pixel,:,1))
    hold on
    title(['Complete S-t courses of chosen voxel (' num2str(depicted_pixel) ') by Smax']);
    DC_legend='Raw signal    ';
    ylim([0 max(DC_data_filter(depicted_pixel,:,1))])
    fig_counter=fig_counter+1;

%% subplot shows beginning of unfiltered and filtered data of depicted_pixel for comparison
subplot(fig_row,fig_col,fig_counter); %#2
    plot(DC_data_filter(depicted_pixel,1:depicted_timepoint,1))
    hold on
    title(['Beginning of S-t courses of voxel (' num2str(depicted_pixel) ')']);
    DC_legend='Raw signal    ';
    ylim([0 max(DC_data_filter(depicted_pixel,:,1))])
 
%% subplot shows unfiltered 2D data
% Filter 1 = raw
z_img=1;
subplot(fig_row,fig_col,[fig_col*z_img+1:fig_col*z_img+2])
    imagesc(DC_data_filter(:,:,1));
    title('Raw data');

%% temporary dummy figure
subplot(fig_row,fig_col,3)
%% Baseline correction (n=2) and LowPass filter (n=3)
%% GUI    
prompt2 = {'Filter typ','low frequency pass','low freqency stop',...
        'Amplitude pass','Amplitude stop',...
        'high frequency stop','high frequency stop',...
        };
    dlg_title2 ='Parameter for filtering (if You dont know, DONT TOUCH!)';
    num_lines2 = 1;
    defaultans2 = {'1','0.3','0.32','0.1','50','0.1','0.075'};%0.475,0.5
    answer2 = inputdlg(prompt2,dlg_title2,[1, length(dlg_title2)+30],defaultans2);

%% Definition of filter response figure
str2double(answer2{1,1})
if (str2double(answer2{1,1})==1)
hd=fdesign.lowpass('Fp,Fst,Ap,Ast',...
    str2double(answer2{2,1}),str2double(answer2{3,1}),...
    str2double(answer2{4,1}),str2double(answer2{5,1}),...
    MR_SamplingRate);
    s=['Filter ' hd.Response ':' num2str(hd.Fpass) ';' num2str(hd.Fstop) ...
        ';' num2str(hd.Apass) ';' num2str(hd.Astop) ...
        ';' num2str(hd.Fs)]
else
hd=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
    str2double(answer2{7,1}),str2double(answer2{6,1}),...
    str2double(answer2{2,1}),str2double(answer2{3,1}),...
    str2double(answer2{5,1}),str2double(answer2{4,1}),str2double(answer2{5,1}),...
    MR_SamplingRate);%bandpass filter F:frequency; A:amplitude
    s=['Filter ' hd.Response ':' ...
          num2str(hd.Fstop1) ';' num2str(hd.Fpass1) ...
        ';' num2str(hd.Fpass2) ';' num2str(hd.Fstop2) ...
        ';' num2str(hd.Apass) ';' num2str(hd.Astop1)  ';' num2str(hd.Astop2) ...
        ';' num2str(hd.Fs)]
end

disp('calculating ...'); % next steps take some time

dfilt=design(hd,'equiripple');
fvtool(dfilt)
xlim([0 2])
ylim([-50 5])
[h,w]=freqz(dfilt,4096)
%%

subplot(fig_row,fig_col,4)
    plot(w*MR_SamplingRate/(2*pi),db(h)); %x-axis: normalized frequency * Nyquist/pi
    xlim([0 2])
    ylim([-50 5])
    title(s)
    xlabel('Frequency (Hz)');
    ylabel('dB');
    grid on

    %%
for zi=1:size_DC_data(1)
    hel=DC_data_long(zi,:);
    hel=transpose(hel);
    helb=DC_Baseline(hel,Rep_Stim,Rep_Stim/2);
    helf=filtfilt(dfilt.Numerator,1,helb); %zero-phase filtering

    DC_data_filter_baseline(zi,:)=helb;
    DC_data_filter_FIR(zi,:)=helf;
end

%%
% DC_data_filter: n=2 -> baseline corrected
DC_data_filter(:,:,2)=DC_data_filter_baseline;
s1=['### DC_data_filter (:,:,2)= Baseline Corrected  \n']; fprintf(prot_fid,s1);disp(s1)

% DC_data_filter: n=3 -> low pass filtered
DC_data_filter(:,:,3)=DC_data_filter_FIR;
s1=['### DC_data_filter (:,:,3)= ' s ' \n']; fprintf(prot_fid,s1);disp(s1)

%% output of calculation 
% subplot shows unfiltered and filtered data of depicted_pixel for comparison
figure(fig)
subplot(fig_row,fig_col,1)  % S-t-curve
    plot(DC_data_filter(depicted_pixel,:,2),'g')
    plot(DC_data_filter(depicted_pixel,:,3),'r')
% subplot shows beginning of filtered data of depicted_pixel for comparison
subplot(fig_row,fig_col,2)  % S-t-curve at beginning
    plot(DC_data_filter(depicted_pixel,1:depicted_timepoint,2),'g')
    plot(DC_data_filter(depicted_pixel,1:depicted_timepoint,3),'r')
    ylim([0 max(DC_data_filter(depicted_pixel,:,2))])

% subplot shows baseline filtered 2D data
z_img=z_img+1;
% #2: Baseline Correction
subplot(fig_row,fig_col,[fig_col*z_img+1:fig_col*z_img+2])
    d=round(mean(DC_data_filter(:,:,2),'all'));
    ds=round(mean(std(DC_data_filter(:,:,2))))
    clims=[d d+ds]
    imagesc( DC_data_filter(:,:,2),clims)
    title('Baseline Correction');
    colorbar

% % subplot shows filtered 2D data
z_img=z_img+1;
% #3: filtered data
subplot(fig_row,fig_col,[fig_col*z_img+1:fig_col*z_img+2])
    clims=[d d+ds]

    imagesc(DC_data_filter(:,:,3),clims)
    title('Low Pass Filter');
    colorbar

%% Median filter in time
s=['### - Median filter in time:  \n']; fprintf(prot_fid,s);disp(s)
disp('calculating ...');
for zx=1:size_DC_data(1)
    DC_data_filter_median(zx,:)=medfilt1(DC_data_filter(zx,:,3),Medfilter_length);
end
DC_data_filter(:,:,4)=DC_data_filter_median;
s=['### DC_data_filter (:,:,4)= Median filtered in time  \n']; fprintf(prot_fid,s);disp(s)

%% subplot shows Median filtered 2D data
z_img=z_img+1;
%Filter 4 = Median Filter
subplot(fig_row,fig_col,[fig_col*z_img+1:fig_col*z_img+2])
    clims=[d d+ds]
imagesc(DC_data_filter(:,:,4),clims)
    title('Medianfilter in time');
    colorbar

%% subplot shows unfiltered and filtered data of depicted_pixel for comparison
subplot(fig_row,fig_col,1)
plot(DC_data_filter(depicted_pixel,:,4),'-k')

%% subplot shows beginning of unfiltered and filtered data of depicted_pixel for comparison
subplot(fig_row,fig_col,2)
plot(DC_data_filter(depicted_pixel,1:depicted_timepoint,4),'-k')
ylim([0 2*mean(DC_data_filter(depicted_pixel,1:depicted_timepoint,4))])

%% subplot shows signal profile over space of unfiltered and filtered data of depicted_pixel
Signal_profile_DC_data_long=mean(DC_data_long(:,:),3);
size(Signal_profile_DC_data_long)
subplot(fig_row,fig_col,3)
    plot(Signal_profile_DC_data_long/max(Signal_profile_DC_data_long),'b')
    title('Signal profile of un-/ filtered data');
hold on


Signal_profile_DC_data_filter_Median=mean(DC_data_filter(:,:,4),2);
size( Signal_profile_DC_data_filter_Median)
plot( Signal_profile_DC_data_filter_Median/max( Signal_profile_DC_data_filter_Median),'r')
DC_legend=vertcat(DC_legend,'Median time F ');
hold off

%% 2d spectral analysis depiction in db
for zi=1:size(DC_data_filter,3),
    subplot(fig_row,fig_col,[fig_col*zi+3:fig_col*zi+4])
    dummy_filter=DC_data_filter(:,:,zi);size(dummy_filter)
    d2_fft=fft(dummy_filter,NFFT,2);size(d2_fft)
    logd2_fft=20*log10(abs(d2_fft(:,2:end/2)));
    if (zi==1),
        c2=round(max(max(logd2_fft(:,200:800))))
        c1=round(mean (mean(logd2_fft(:,200:800)))+1.25*mean(std(logd2_fft(:,200:800))))
        ydb_max=logd2_fft(depicted_pixel,2);
    end
        clims=[ c1 c2]
        imagesc((logd2_fft(:,1:NFFT*0.1)),clims)
        colorbar
        ylabel('voxel')
        title(['Spectrum in db [0-2Hz]'])
    subplot(fig_row,fig_col,fig_col*zi+5)
    plot(F2fft(2:end),logd2_fft(depicted_pixel,:))
        ylim([ydb_max/2 ydb_max])
        xlabel('Hz');
        title(['Spectrum: Data of voxel (' num2str(depicted_pixel) ')']);
        xlim(xlim_spectrum)
    if (zi==1), spect_raw=logd2_fft(depicted_pixel,:); end
    if (zi==3), spect_filtered=logd2_fft(depicted_pixel,:); end
end

%% save analysis figure
image_suffix='Subroutine_DC_filter_data';
old=cd(['../data/' analyse_animalfolder '/']);
    DC_hg_name=[ animal_ident '_' datum '_' image_suffix '.emf' ];
saveas(fig,DC_hg_name,'emf')
cd(old)
%%
% returned data arr: DC_data_filter (voxel,s-t-course,filter)
err_log=0;
s=['### (' analyse_function_name_version ') - end \n']; fprintf(prot_fid,s);disp(s)

