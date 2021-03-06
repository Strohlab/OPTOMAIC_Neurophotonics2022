%% OPTOMAIC 
%% - Analysis Program - 
%% Version 19.5.2022 
%% Author: Dirk Cleppien

%% data structure:
%% ../[folder study_name]
%% ../[folder study_name]/matlab                    matlab scripts
%% ../[folder study_name]/data                      all measured animal data
%% ../[folder study_name]/data/[animal name]        single measured animal data + Spike2 file
%% ../[folder study_name]/data/[animal name]/dicom  dicom files of animal

%% Initialisation
clear all
close all
fclose all

%% --- Initialization of Protocol file ---
[datum,prot_fid,prot_filename,analyse_animalfolder]=DC_init_prot;

%% --- Create a backup analysis data file ---
animal_ident=input('Short identifier of the animal:','s');
m_file=['../data/' analyse_animalfolder '/' animal_ident '-OPTOMAIC-' datum '.mat']
save (m_file);

%%
%file_results=['../data/' analyse_animalfolder '/Results.txt']

%% - global parameter
s=['### Global Variables: \n']; fprintf(prot_fid,s);disp(s)
g_Stim_Interval=10;s=['### Global Var: g_Stim_Interval = ' num2str(g_Stim_Interval) ' (Optical stimulation interval [s] ) \n']; fprintf(prot_fid,s);disp(s)
g_DC_MR_SamplingRate=0.05;s=['### Global Var: g_DC_MR_SamplingRate = ' num2str(g_DC_MR_SamplingRate) ' (Sampling Rate [s] MR signal ) \n']; fprintf(prot_fid,s);disp(s)
g_DC_Spike2_SamplingRate=2000;s=['### Global Var: g_DC_Spike2_SamplingRate = ' num2str(g_DC_Spike2_SamplingRate) ' (Sampling Rate [ms] in Spike - should be the same for all measurements) \n']; fprintf(prot_fid,s);disp(s)
g_DC_sync_StimMRI=100;s=['### Global Var: g_DC_sync_StimMRI = ' num2str(g_DC_sync_StimMRI) ' (Delay of MRI to stimulation pulse - individual for every measurement) \n']; fprintf(prot_fid,s);disp(s)
g_Rep_stim=g_Stim_Interval/g_DC_MR_SamplingRate;s=['### Global Var: g_Rep_stim = ' num2str(g_Rep_stim) ' (one stimulation interval) \n']; fprintf(prot_fid,s);disp(s)
g_cutoff_stst=g_Rep_stim;s=['### Global Var: g_cutoff_stst = ' num2str(g_cutoff_stst) ' (one stimulation interval - check, if it is ok!) \n']; fprintf(prot_fid,s);disp(s)
g_MRI_cutoff_stst_ms=g_Stim_Interval*g_DC_Spike2_SamplingRate;s=['### Global Var: g_MRI_cutoff_stst_ms = ' num2str(g_MRI_cutoff_stst_ms) ' (Cutoff used for MRI data adapting for Spike2 data [s] ) \n']; fprintf(prot_fid,s);disp(s)
g_DC_stim_offset=20;s=['### Global Var: g_DC_stim_offset = ' num2str(g_DC_stim_offset) ' (Baseline before stimulation in time course - for every measurement the same) \n']; fprintf(prot_fid,s);disp(s)
g_DC_baseline_interval=15;s=['### Global Var: g_DC_baseline_interval = ' num2str(g_DC_baseline_interval) ' (Baseline before stimulation in time course - for every measurement the same) \n']; fprintf(prot_fid,s);disp(s)
g_DC_signal_offset=100;s=['### Global Var: g_DC_signal_offset = ' num2str(g_DC_signal_offset) ' (Baseline before stimulation in time course - for every measurement the same) \n']; fprintf(prot_fid,s);disp(s)
g_DC_signal_interval=10;s=['### Global Var: g_DC_signal_interval = ' num2str(g_DC_signal_interval) ' (Baseline before stimulation in time course - for every measurement the same) \n']; fprintf(prot_fid,s);disp(s)
g_DC_ProtImages={};   %List of generated protocol images of the analysis

% Input necessary due to individual value
% ant_Sync_StimMRI_val=input('Time delay of MRI to stimulation pulse [in s]: '); 

%%----------------------------
%% ########################### 
%% OPTOMAIC: Step 1 + 2 - Stimulus and Calcium signal traces
%% ###########################
%% --- Preparation of data ---
%% --- Read of Spike2 data ---
[DC_VisStim_signal,DC_VisStim_onsets,DC_Carec_signal,DC_Carec_OnOffsets,DC_sync_delay_s,err]=DC_read_Spike2(g_MRI_cutoff_stst_ms,g_DC_Spike2_SamplingRate,datum,prot_fid,analyse_animalfolder,animal_ident,g_DC_ProtImages);
% DC_Breathfilter_BS: in kHz
% DC_Breathdata in 2kHz sampling
g_DC_sync_StimMRI=round(DC_sync_delay_s/g_DC_MR_SamplingRate)
save(m_file); %Backup
whos('-file',m_file)
disp(['Workspace saved in ' m_file]);

%% downsampling DC_CaSignal to MRI time resolution
DC_downsampling=g_DC_Spike2_SamplingRate*g_DC_MR_SamplingRate
DC_CaSignal_dsampled=zeros(round(size(DC_Carec_signal,1)/DC_downsampling),1);
size(DC_CaSignal_dsampled)
DC_Breathdata_dsampled=zeros(round(size(DC_Carec_signal,1)/DC_downsampling),1);
size(DC_Breathdata_dsampled)

for zi=2:size(DC_CaSignal_dsampled,1),
    zi;
   DC_CaSignal_dsampled(zi)=mean(DC_Carec_signal((zi-1)*100:zi*100)) ;
end

%% --- Analysis of slow waves per stimulation interval
err_log=0;
g_SW_DetectWindow=[800 17000] ;
if (err_log==0),[Slwaves_arr, Slwaves_arr_sel, err_log]=DC_analysis_Slwaves(DC_CaSignal_dsampled,DC_Carec_OnOffsets,DC_VisStim_onsets,g_DC_Spike2_SamplingRate,g_SW_DetectWindow,datum,prot_fid,analyse_animalfolder,animal_ident,g_DC_ProtImages);end
size(Slwaves_arr)
% 1: corresponding stim interval
% 2: response time after stim 
% 3: interval between current and previous slwave onset
% 4: interval between current and previous slwave offset
% 5: number of slow waves per stim interval
%% ########################### 
%% OPTOMAIC: Step 3 - Preproc fMRI signal trace
%% ###########################
%% FMRI Preprocessing
%% --- read of images
[DC_Bilder, err_log]=DC_read_Bilder(prot_fid,analyse_animalfolder);
save(m_file); %Backup

%% - reshape of images to data_long and sync with stimulation pulse
if (err_log==0),[DC_data_long,mean_DC_spatial_profile, g_DC_ProtImages, err_log]=DC_reshape_Bilder(DC_Bilder,g_DC_sync_StimMRI,g_DC_stim_offset,g_cutoff_stst,g_Rep_stim,datum,g_DC_ProtImages,err_log,prot_fid,analyse_animalfolder,animal_ident);end;
save (m_file); %Backup

%% - filter of images
if (err_log==0),[DC_data_filter_arr,g_DC_ProtImages, err_log]=DC_filter_data(DC_data_long,g_Rep_stim,1/g_DC_MR_SamplingRate,datum,g_DC_ProtImages,err_log,prot_fid,analyse_animalfolder,animal_ident);end;
size(DC_data_filter_arr)    % ehemals DC_data_filter!!!
save (m_file); %Backup

%% global standardization of every filter data
glob_SI=zeros(size(DC_data_filter_arr,1),size(DC_data_filter_arr,3));
size(glob_SI)
%%
for zfilt=1:size(DC_data_filter_arr,3),
   for zv=1: size(DC_data_filter_arr,1),
        glob_SI(zv)=mean(DC_data_filter_arr(zv,:,zfilt));
   end
end
%%
DC_glob_SI=mean(DC_data_filter_arr(:,:,1),2);
size(DC_glob_SI)
%% standardization of every filter data
sDC_data_filter_arr=DC_data_filter_arr;
sDC_data_filter_arr(:)=0;
%%
for zfilt=1:size(DC_data_filter_arr,3),
   for zv=1: size(DC_data_filter_arr,1),
        sDC_data_filter_arr(zv,:,zfilt)=DC_data_filter_arr(zv,:,zfilt)./glob_SI(zv);
   end
end

%% - collection of all stimuli into 3D stack - MRI
%% cut into stimulation intervals
 
if (err_log==0),[DC_data_collect_arr,g_DC_ProtImages, err_log]=DC_collect_data_stim_MRI(sDC_data_filter_arr,g_Rep_stim,datum,g_DC_ProtImages,err_log,prot_fid,analyse_animalfolder,animal_ident);end;
save (m_file); %Backup
%% ########################### 
%% OPTOMAIC: Step 4 - Neural event x BOLD response
%% ###########################

%% GUI Onset Analysis
BOLD_condition_end=1;
prompt2 = {'1. Printout Figures[j/n]:','2. Threshold(voxels per event) [% of FOV]: ','3. Filter Nr: ',...
        '4. Frontal Offset (voxel number):','5. FOV interval (velocity fit):','6. Voxel distance [mm]:'...
        '7. Baseline window end:','8. Cond.#2:BOLD-Window - start','9. Cond.#2: BOLD-Window - end',...
        '10. Cond.#2: S(BOLD)-ratio start','11. Cond.#2: S(BOLD)-ratio increment','12. Cond.#2: Cond.#2: S(BOLD)-ratio end'...
        '13. FOV -lower boundary','14. FOV -upper boundary','15. Threshold(events per voxel)',...
        '16. R^2 (HRF-Fit):'};
    dlg_title2 ='Parameter for onset calculation (if You dont know, DONT TOUCH!)';
    num_lines2 = 1;
    defaultans2 = {'j','30','3','9','15','0.3125',...
        '20','80','150','1.03','0.1','1.03'...
        '11','46','20','0.8'};%0.475,0.5
    answer2 = inputdlg(prompt2,dlg_title2,[1, length(dlg_title2)+30],defaultans2);
%%
show_fig=1;
ant_BOLD_condition=1

%% Initialisation of Analysis
FOV_voxel=[str2double(answer2{13,1}) str2double(answer2{14,1})]
z_filter=str2double(answer2{3,1})
DC_Events=DC_data_collect_arr(FOV_voxel(1):FOV_voxel(2),:,:,z_filter);
size(DC_Events)
anz_DC_Events=size(DC_Events,3)
anz_DC_voxels=size(DC_Events,1)
size(DC_Events)

%%         Pure events (no spontaneous slwaves in interval before
%%         [depends on menue]
d_timeInterval=Slwaves_arr_sel(:,1).*(Slwaves_arr_sel(:,6)~=0);% Choice of selected ONset-ONset interval
DC_slwaves_sel=zeros(length(d_timeInterval),1)+1;
size(DC_slwaves_sel)
d_sel=DC_VisStim_onsets(d_timeInterval(d_timeInterval>0));
DC_slwaves_sel=d_timeInterval(d_timeInterval>0);       % all selected stim intervals!

[i,j]=max(DC_slwaves_sel)
if (i>size(DC_data_collect_arr,3)),
    DC_slwaves_sel=DC_slwaves_sel(1:j-1);
end
if (ant_BOLD_condition==1),
    sel_DC_Events=DC_data_collect_arr(FOV_voxel(1):FOV_voxel(2),:,DC_slwaves_sel,z_filter);
else
    sel_DC_Events=DC_data_collect_arr(FOV_voxel(1):FOV_voxel(2),:,:,z_filter);
end
%
size(sel_DC_Events)
anz_DC_Events=[anz_DC_Events size(sel_DC_Events,3)]
anz_DC_voxels=[anz_DC_voxels size(sel_DC_Events,1)]

%%         per event: BOLD signal higher factor*Baseline signal
%%         [factor z_ratio from menue]
z_ratio=str2double(answer2{10,1})
size(sel_DC_Events)
matrix_sel_Events=zeros(size(sel_DC_Events,1),size(sel_DC_Events,3));
size(matrix_sel_Events)
% prepare of ratio matrix
for z_voxel=1:size(sel_DC_Events,1)
    for z_event=1:size(sel_DC_Events,3),
        s_Baseline=mean(sel_DC_Events(z_voxel,...
            1:str2double(answer2{7,1}),z_event))
        s_BOLD=mean(sel_DC_Events(z_voxel,str2double(answer2{8,1}):...
            str2double(answer2{9,1}),z_event))
        matrix_sel_Events(z_voxel,z_event)=s_BOLD/s_Baseline;
    end
end
% matrix of collected events over threshold
DC_matrix_thres_sel_Events=(matrix_sel_Events>z_ratio);
size(DC_matrix_thres_sel_Events)
% collected events per voxel
Events_per_voxel=sum(DC_matrix_thres_sel_Events,2)
size(Events_per_voxel)
mEvents_per_voxel=round(mean(Events_per_voxel))
anz_DC_Events=[anz_DC_Events mEvents_per_voxel]

%%          mean image of collected events per voxel in FOV
msel_DC_Events_BOLD=zeros(size(sel_DC_Events,1),size(sel_DC_Events,2));
size(msel_DC_Events_BOLD)
for z_voxel=1:size(sel_DC_Events,1),
    msel_DC_Events_BOLD(z_voxel,:)=mean(sel_DC_Events...
    (z_voxel,:,DC_matrix_thres_sel_Events(z_voxel,:)),3);
end
    
%% Result depiction
s=[animal_ident ': Mean image (selected events) per filter (' datestr(now,'yymmdd-HHMMSS') ')']
fig=figure('Name',s,'units','normalized','outerposition',[0 0 1 1]);
    fig.ToolBar='none';
fig.NumberTitle='off';
fig.Name=[animal_ident ' - Mean image (selected events) per Filter']; 
fig.FileName=[animal_ident '-MeanImage_SelEvents_PerFilter']; 
fig_row=4;
fig_col=3;

%%
fig_counter=1;
subplot(fig_row,fig_col,fig_counter)
    clims=[0.8 1.2]
    imagesc(matrix_sel_Events,clims)
    colorbar
    title([' Ratio(S-BOLD/S-Baseline); Events (' num2str(anz_DC_Events(2)) ...
        '/' num2str(anz_DC_Events(1)) ')'])
    xlabel('time [arbitrary]')
    ylabel('FOV[voxel]')
    fig_counter=fig_counter+1;
subplot(fig_row,fig_col,fig_counter)
    clims=[0 1]
    imagesc((DC_matrix_thres_sel_Events),clims)
    colorbar
    title([' Ratio(S-BOLD/S-Baseline) - selected events(yellow)'])
    xlabel('time [arbitrary]')
    ylabel('FOV[voxel]')
    fig_counter=fig_counter+1;
subplot(fig_row,fig_col,fig_counter);
    plot(sum(DC_matrix_thres_sel_Events,1))
    hold on
    ylim([1 (str2double(answer2{14,1})-str2double(answer2{13,1}))])
    title([' Number of voxels per selected event'])
    xlabel('event')
    ylabel('Number of voxel')
    fig_counter=fig_counter+1;
subplot(fig_row,fig_col,fig_counter)
    tx=1:size(Events_per_voxel,1);
    plot(Events_per_voxel)
    hold on
    ty=zeros(size(Events_per_voxel,1))+str2double(answer2{15,1})
    plot(tx,ty)
    ty=zeros(size(Events_per_voxel,1))+mEvents_per_voxel
    plot(tx,ty,'--')
    title(['Number of events per voxel ; Events (Mean(' num2str(anz_DC_Events(3)) ...
        ')/' num2str(anz_DC_Events(1)) ')'])
    xlabel('FOV (voxel)')
    ylabel('number of events')
    ylim([0 max(Events_per_voxel)+5])
    xlim([1 size(Events_per_voxel,1)])
    fig_counter=fig_counter+1;
     text(-10,max(Events_per_voxel)+20,'A) selected events per voxel')
subplot(fig_row,fig_col,fig_counter)
    clims=[0.9 1.05]
    imagesc(msel_DC_Events_BOLD,clims)
    colorbar
    title([' Ratio(S-BOLD/S-Baseline)'])
    xlabel('time [arbitrary]')
    ylabel('FOV[voxel]')
    fig_counter=fig_counter+1;
    fig_counter=fig_counter+1;

size(msel_DC_Events_BOLD)
anz_DC_Events

%% %% ########################### 
%% OPTOMAIC: Step 5 - BOLD response matrix
%% ###########################

Eventthreshold_voxel=str2double(answer2{15,1})
thres_Events_BOLD=zeros(size(Events_per_voxel,1),1);
size(thres_Events_BOLD)
thres_Events_BOLD= Events_per_voxel> Eventthreshold_voxel

%%
anz_DC_voxels=[anz_DC_voxels sum( thres_Events_BOLD)]

%%
thres_msel_DC_Events_BOLD=msel_DC_Events_BOLD;
% set lines under threshold to zero
for zvoxel=1:size(Events_per_voxel,1),
    if ( thres_Events_BOLD(zvoxel)==0),
        thres_msel_DC_Events_BOLD(zvoxel,:)=0; 
    end
end

%% Depiction of results
subplot(fig_row,fig_col,fig_counter)
    tx=1:size(Events_per_voxel,1);
    plot(Events_per_voxel)
    hold on
    plot(thres_Events_BOLD.*Events_per_voxel,'*-r');
    ty=zeros(size(Events_per_voxel,1))+str2double(answer2{15,1})
    plot(tx,ty)
    ty=zeros(size(Events_per_voxel,1))+mEvents_per_voxel
    plot(tx,ty,'--')
    title(['Number of events per voxel; voxels = (' num2str(anz_DC_voxels(3))...
        '/' num2str(anz_DC_voxels(2)) ')'])
    xlabel('FOV (voxel)')
    ylabel('number of events')
    ylim([0 max(Events_per_voxel)+5])
    xlim([1 size(Events_per_voxel,1)])
    fig_counter=fig_counter+1;
    text(-10,max(Events_per_voxel)+20,['B) selected voxels after event threshold (=' ...
        answer2{15,1} ')'])
subplot(fig_row,fig_col,fig_counter)
    clims=[0.9 1.05]
    imagesc(thres_msel_DC_Events_BOLD,clims)
    colorbar
    title([' Ratio(S-BOLD/S-Baseline) - Event threshold'])
    xlabel('time [arbitrary]')
    ylabel('FOV[voxel]')
    fig_counter=fig_counter+1;
    fig_counter=fig_counter+1;
    
    %% ########################### 
    %% OPTOMAIC: Step 6 - Mean BOLD response
    %% ###########################
    % Depiction of Signal time curve per voxel
    %%
    s1='Onset detection per voxel'
    if (show_fig==1),
        figfit=figure('units','normalized','outerposition',[0 0 1 1]);
        figfit.ToolBar='none';
        figfit.NumberTitle='off';

        figfit.Name=[animal_ident s ': ' s1];
        figfit.FileName=[animal_ident s '_' s1];
        figfit_row=6;
        figfit_col=7;
        figfit_counter=1;

        subplot(figfit_row,figfit_col,figfit_counter)
    end

    %%    
    HRF_thres_Events_BOLD=thres_Events_BOLD;
    HRF_fitted_arr=msel_DC_Events_BOLD;
    HRF_fitted_arr(:)=0;
    size(HRF_fitted_arr)
    R2_arr=[];
    ylim_max_arr=[];
    HRF_param=[];
    dummy=msel_DC_Events_BOLD;
    
    %%
    for zvoxel=1:size(Events_per_voxel,1),
        %%
        zvoxel
        subplot(figfit_row,figfit_col,zvoxel)
        hold on

        %%  HRF fitting
        [y_min,x_min]=min( dummy(zvoxel,:))
        [y_max,x_max]=max(dummy(zvoxel,:))
        x_HRF=1:size(dummy,2)
        x_HRF=x_HRF'
        ydata_HRF=dummy(zvoxel,:)-y_min;
        hold on
        plot(x_HRF,ydata_HRF)
        A=0.5;T0=70;b=1;a=3;
        f_HRF_1=@(p,x_HRF) ((x_HRF-p(2)).^(a-1));
        f_HRF_2=@(p,x_HRF) f_HRF_1(p,x_HRF).*p(3).^a;
        f_HRF_3=@(p,x_HRF) f_HRF_2(p,x_HRF)./gamma(a);
        f_HRF_4=@(p,x_HRF) f_HRF_3(p,x_HRF).*exp(-(x_HRF-p(2)).*p(3));
        f_HRF_5=@(p,x_HRF) p(1)*f_HRF_4(p,x_HRF);
        f_HRF_6=@(p,x_HRF) f_HRF_5(p,x_HRF).*(x_HRF > (p(2)-1));
        % plot(x_HRF,f_HRF_6([A,T0,b,a],x_HRF))
        mdl=fitnlm(x_HRF,ydata_HRF,f_HRF_6,[A,T0,b,a])
        HRF_fitted_arr(zvoxel,:)=predict(mdl,x_HRF);
            plot(x_HRF,HRF_fitted_arr(zvoxel,:),'-r');
        
        if (mdl.Rsquared.Adjusted<str2double(answer2{16,1})),
            HRF_thres_Events_BOLD(zvoxel)=0;
        end

        R2_arr=[R2_arr mdl.Rsquared.Adjusted]
        ylim_max_arr=[ylim_max_arr max(ydata_HRF)];
        HRF_param=vertcat(HRF_param, mdl.Coefficients.Estimate')
        title({[num2str(zvoxel) ': R^2(adj)=' num2str(R2_arr(zvoxel))];...
            ['!' num2str(HRF_thres_Events_BOLD(zvoxel)) '; ' ...
            num2str(HRF_param(zvoxel,1)) ';'...
            num2str(HRF_param(zvoxel,2)) ';'...
            num2str(HRF_param(zvoxel,3)) ';'...
            num2str(HRF_param(zvoxel,4))]})
        d_arr=x_HRF;d_arr(1:20)=0;d_arr(21:end)=max(HRF_fitted_arr(zvoxel,:));plot(d_arr,'--k')
        hold off

    end
    
    %%
    HRF_thres_msel_DC_Events_BOLD=thres_msel_DC_Events_BOLD;
    for zvoxel=1:size(HRF_thres_msel_DC_Events_BOLD,1),
        if (HRF_thres_Events_BOLD(zvoxel)==0),
            HRF_thres_msel_DC_Events_BOLD(zvoxel,:)=0;
        end
    end
    anz_DC_voxels=[anz_DC_voxels sum(HRF_thres_Events_BOLD)]
   
    %%
    figure(fig)
    subplot(fig_row,fig_col,fig_counter)
        tx=1:size(Events_per_voxel,1);
        plot(Events_per_voxel)
        hold on
        plot(HRF_thres_Events_BOLD.*Events_per_voxel,'*-k');
        ty=zeros(size(Events_per_voxel,1))+str2double(answer2{15,1})
        plot(tx,ty)
        ty=zeros(size(Events_per_voxel,1))+mEvents_per_voxel
        plot(tx,ty,'--')
                title(['Number of events per voxel; voxels = (' num2str(anz_DC_voxels(4))...
            '/' num2str(anz_DC_voxels(2)) ')'])

        xlabel('FOV (voxel)')
        ylabel('number of events')
        ylim([0 max(Events_per_voxel)+5])
        xlim([1 size(Events_per_voxel,1)])
        fig_counter=fig_counter+1;
             text(-10,max(Events_per_voxel)+20,['C) selected voxels after HRF analysis'])        
   subplot(fig_row,fig_col,fig_counter)
        clims=[0.9 1.05]
        imagesc(HRF_thres_msel_DC_Events_BOLD,clims)
        colorbar
        title([' Ratio(S-BOLD/S-Baseline) - Event threshold'])
        xlabel('time [arbitrary]')
        ylabel('FOV[voxel]')
        fig_counter=fig_counter+1;
        fig_counter=fig_counter+1;
%% ########################### 
%% OPTOMAIC: Step 7 - BOLD response Onset detection
%% ###########################
T50perc_arr=zeros(size(HRF_thres_Events_BOLD,1),1);
T50perc_arr_HRFfit=T50perc_arr;
size(T50perc_arr)

%%
for zvoxel=1:size(HRF_thres_Events_BOLD,1),
    %%
    zvoxel
    figure(figfit)
    subplot(figfit_row,figfit_col,zvoxel)
    hold on

    %%  HRF fitting
    d_plot=zeros(size(dummy,2),1);
    %%
    dummy=msel_DC_Events_BOLD(zvoxel,:);
    dummy_HRF=HRF_fitted_arr(zvoxel,:);
    [y_min,x_min]=min(dummy)
    [y_max,x_max]=max(dummy)
    S50perc=y_min+0.5*(y_max-y_min)
    S50perc_HRFfit=min(dummy_HRF)+0.5*(max(dummy_HRF)-min(dummy_HRF))
    tx=1:size(dummy,2);
    T50perc_arr(zvoxel)=min(tx(dummy>S50perc))
    T50perc_arr_HRFfit(zvoxel)=min(tx(dummy_HRF>S50perc_HRFfit))
    %plot(dummy)
    hold on
    d_plot(:)=0;d_plot(T50perc_arr(zvoxel):end)=y_max-y_min;plot(d_plot,'g-');
    d_plot(:)=0;d_plot(T50perc_arr_HRFfit(zvoxel):end)=y_max-y_min;plot(d_plot,'b--');

    if (T50perc_arr(zvoxel)<str2double(answer2{7,1})),HRF_thres_Events_BOLD(zvoxel)=0;end% baseline end

    %  title({[num2str(zvoxel) ': '];...
    %     ['  ']})
    title({[num2str(zvoxel) ':' num2str(HRF_thres_Events_BOLD(zvoxel)) ';R2=' num2str(R2_arr(zvoxel))...
        '; T=' num2str(T50perc_arr(zvoxel))];...
        [ ...
        num2str(HRF_param(zvoxel,1)) ';'...
        num2str(HRF_param(zvoxel,2)) ';'...
        num2str(HRF_param(zvoxel,3)) ';'...
        num2str(HRF_param(zvoxel,4)) ...
        ]})

    hold off
    %%
    %end
    if (HRF_thres_Events_BOLD(zvoxel)==0),
        hold on
        cross_graph=0:ylim_max_arr(zvoxel)/size(msel_DC_Events_BOLD,2):ylim_max_arr(zvoxel);plot(cross_graph,'.k')
        ylim([0 ylim_max_arr(zvoxel)])
    end

end

%% Figures until here
figonset=fig;
figfit;
%%
%% Result depiction
s=[animal_ident ': Onset_depiction (' datestr(now,'yymmdd-HHMMSS') ')']
fig=figure('Name',s,'units','normalized');
fig.ToolBar='none';
fig.NumberTitle='off';
fig.Name=[animal_ident ' - Onset depiction'];
fig.FileName=[animal_ident '-Onset_depiction'];

%%
x_lin=1:size(T50perc_arr,1);
x_lin=x_lin(HRF_thres_Events_BOLD==1)
Tdata=T50perc_arr(HRF_thres_Events_BOLD==1)
Tdata_HRFfit=T50perc_arr_HRFfit(HRF_thres_Events_BOLD==1)
x_start=str2double(answer2{4,1});
x_end=size(Tdata,1);
plot(T50perc_arr,'--b','LineWidth',2)
hold on
plot(x_lin(x_start:x_end),Tdata(x_start:x_end),'-*b','LineWidth',3)
ylim([0 100])
hold off
title(['Measured data'])

%% save of analysis
save(m_file)

%% - end of analysis
%% ########################### 
%% OPTOMAIC: end
%% ###########################

%% Figures for NeuroPhotonics 2022
%% Simultaneous depiction of Calcium and Stimulus trace

s=[animal_ident ': Calcium + Stimulus trace (' datestr(now,'yymmdd-HHMMSS') ')']
    fig=figure('Name',s,'units','normalized');
        fig.ToolBar='none';
    fig.NumberTitle='off';
    fig.Name=[animal_ident ' - CalciumTrace+Stim']; 
    fig.FileName=[animal_ident '-CalciumStim']; 
    fig_row=1;
    fig_col=1;
    fig_counter=1;
hold off
stim=zeros(length(DC_VisStim_onsets),1)+0.2;
% shift bars to onset start
onset_interval=(DC_VisStim_onsets(4)-DC_VisStim_onsets(2))/4/g_DC_Spike2_SamplingRate;
zoom_factor=0.05;
onset_shift=onset_interval*(zoom_factor);
hold on
xlim([430 640])
x_CA=(1:length(DC_CaSignal_dsampled))*g_DC_MR_SamplingRate; % for Rat73 it was 20min
% filter for depiction moving-average filter
windowSize=20;
b=(1/windowSize)*ones(1,windowSize);
a=1;
y_CA=filter(b,a,DC_CaSignal_dsampled);
my_CA=mean(y_CA)
plot(x_CA,y_CA/(my_CA)-0.75,'g','Linewidth',2)
ylim([-0.3 1.5])
bar((DC_VisStim_onsets(:)/g_DC_Spike2_SamplingRate)+onset_shift-0.3,stim,zoom_factor,'k');
% part: resulting BOLD intervals
% --


d_timeInterval=Slwaves_arr_sel(:,1).*(Slwaves_arr_sel(:,6)~=0);% Choice of selected ONset-ONset interval
DC_slwaves_sel=zeros(length(d_timeInterval),1)+1;
size(DC_slwaves_sel)
d_sel=DC_VisStim_onsets(d_timeInterval(d_timeInterval>0));
DC_slwaves_sel=d_timeInterval(d_timeInterval>0);       % all selected stim intervals!
stim=zeros(length(d_sel),1)+0.8;

plot(d_sel/g_DC_Spike2_SamplingRate,stim,'b*','LineWidth',2);
ylim([-0.2 2])
xlim([-20 x_CA(end)+20])
title(['Combined Depiction of MRI Signal, Calcium Recordings and Stimulation pulses ']);
hold off
     image_suffix=['DepictionCombinedSignals_Voxel'];
%     comment=['Depiction of synchronised signal. Unnecessary signal is cut off.'];
old=cd(['../data/' analyse_animalfolder '/']);
    DC_hg_name=[ animal_ident '_' datum '_' image_suffix '.emf' ];
saveas(fig,DC_hg_name,'emf')
cd(old)

s=[animal_ident ': BOLD response intervals (' datestr(now,'yymmdd-HHMMSS') ')']
    fig=figure('Name',s,'units','normalized');
    fig.ToolBar='none';
    fig.NumberTitle='off';
    fig.Name=[animal_ident ' - BOLD response intervals']; 
    fig.FileName=[animal_ident '-BOLDintervals']; 
    fig_row=1;
    fig_col=1;
    fig_counter=1;
d_sel=DC_VisStim_onsets(d_timeInterval(d_timeInterval>0));
y_fMRI=zeros(round(d_sel(end)/g_DC_Spike2_SamplingRate),3);
size(y_fMRI)
for zi=1:size(d_sel,2)-1,
    y_val=round(d_sel(zi)/g_DC_Spike2_SamplingRate)
    y_fMRI(y_val:y_val+8,2)=-2;
end
for zi=10:10:size(y_fMRI,1)-1,
    y_fMRI(zi:zi+8,2)=y_fMRI(zi:zi+8,2)+1;
end
size(y_fMRI)
colormap(jet)
imagesc(y_fMRI')
xlim([430 640])
%% Standardized Mean BOLD all vs. Mean BOLD sel.

size(DC_Events)
Image_Mean_all=mean(DC_Events,3);
Image_Mean_all(:,20)=max(Image_Mean_all(:))
cmap=[1 1.03]
s=[animal_ident ': MeanBOLD_all (' datestr(now,'yymmdd-HHMMSS') ')']
    fig=figure('Name',s,'units','normalized');
    fig.ToolBar='none';
    fig.NumberTitle='off';
    fig.Name=[animal_ident ' - MeanBOLD_all']; 
    fig.FileName=[animal_ident '-MeanBOLD_all']; 
    fig_row=1;
    fig_col=1;
    fig_counter=1;
    imagesc(Image_Mean_all,cmap)
    colorbar
    title([' Mean Image (all)'])
    xlabel('time [arbitrary]')
    ylabel('FOV[voxel]')
s=[animal_ident ': MeanBOLD_sel (' datestr(now,'yymmdd-HHMMSS') ')']
    fig=figure('Name',s,'units','normalized');
    fig.ToolBar='none';
    fig.NumberTitle='off';
    fig.Name=[animal_ident ' - MeanBOLD_sel']; 
    fig.FileName=[animal_ident '-MeanBOLD_sel']; 
    fig_row=1;
    fig_col=1;
    fig_counter=1;
    dummy_arr=HRF_thres_msel_DC_Events_BOLD;
    dummy_arr(:,20)=max(dummy_arr(:));
    imagesc(dummy_arr,cmap)
    colorbar
    title([' Mean Image (selected'])
    xlabel('time [arbitrary]')
    ylabel('FOV[voxel]')

%% BOLD response fit to Mean BOLD (voxel) all vs. selected
s=[animal_ident ': BOLD response fit (' datestr(now,'yymmdd-HHMMSS') ')']
    fig=figure('Name',s,'units','normalized');
    fig.ToolBar='none';
    fig.NumberTitle='off';
    fig.Name=[animal_ident ' - BOLD response fit']; 
    fig.FileName=[animal_ident '-BOLDfit']; 
    fig_row=1;
    fig_col=1;
    fig_counter=1;
    for zNP=1:2,
        if (zNP==1),s_Farbe='k-';
        else s_Farbe='b-';end
        R2_arr_NP=[];
        if (zNP==1),
            %dummy=Image_Mean_all(3:36,:);
            dummy=Image_Mean_all(:,:);
        else
            dummy=HRF_thres_msel_DC_Events_BOLD;
        end
        size(dummy)
        %keyboard
        for zvoxel=1:size(dummy,1),
            x_HRF=1:size(dummy,2);
            [y_min,x_min]=min( dummy(zvoxel,:));
            [y_max,x_max]=max(dummy(zvoxel,:));
            x_HRF=x_HRF'
            ydata_HRF=dummy(zvoxel,:)-y_min;
            plot(x_HRF,ydata_HRF,s_Farbe,'LineWidth',5)
            hold on
            A=0.5;T0=70;b=1;a=3;
            f_HRF_1=@(p,x_HRF) ((x_HRF-p(2)).^(a-1));
            f_HRF_2=@(p,x_HRF) f_HRF_1(p,x_HRF).*p(3).^a;
            f_HRF_3=@(p,x_HRF) f_HRF_2(p,x_HRF)./gamma(a);
            f_HRF_4=@(p,x_HRF) f_HRF_3(p,x_HRF).*exp(-(x_HRF-p(2)).*p(3));
            f_HRF_5=@(p,x_HRF) p(1)*f_HRF_4(p,x_HRF);
            f_HRF_6=@(p,x_HRF) f_HRF_5(p,x_HRF).*(x_HRF > (p(2)-1));
            % plot(x_HRF,f_HRF_6([A,T0,b,a],x_HRF))
            mdl=fitnlm(x_HRF,ydata_HRF,f_HRF_6,[A,T0,b,a])
            HRF_fitted_arr(zvoxel,:)=predict(mdl,x_HRF);
            plot(x_HRF,HRF_fitted_arr(zvoxel,:),':r','LineWidth',6);
            xlim([-10 210])
            if (zNP==1),ylim([-0.02 0.13]);else ylim([-0.05 0.25]);end
            R2_arr_NP=[R2_arr_NP mdl.Rsquared.Adjusted]
            hold off
        end
        if (zNP==1),
            R2_all=R2_arr_NP;
        else
            R2_sel=R2_arr_NP;
        end
    end
    R2_test=vertcat(R2_all,R2_sel)
    figure(50)
    b=boxchart(R2_test','Notch','on')
    b.LineWidth=3
    ylim([0 1.1])
    [p,tbl,stats] = kruskalwallis(R2_test')
%% R2 depiction all vs. sel events
s=[animal_ident ': BOLD response fit - R^2(' datestr(now,'yymmdd-HHMMSS') ')']
    fig=figure('Name',s,'units','normalized');
    fig.ToolBar='none';
    fig.NumberTitle='off';
    fig.Name=[animal_ident ' - BOLD response fit - R^2']; 
    fig.FileName=[animal_ident '-BOLDfit_R2']; 
    fig_row=1;
    fig_col=1;
    fig_counter=1;
plot(R2_all,'k:','LineWidth',2)
hold on
plot(R2_all,'ko','LineWidth',2,'MarkerSize',6)
hold on
plot(R2_sel,'b-*','LineWidth',2)
plot(R2_sel,'b*','LineWidth',2,'MarkerSize',6)

ylim([0.5 1])
[y1,x1]=max(R2_all)
[y2,x2]=max(R2_sel)
ylim([-0.1 1.1])
%% depiction of BOLD reponse matrix
% detected events blue, other red
s=[animal_ident ': BOLD response matrix (' datestr(now,'yymmdd-HHMMSS') ')']
    fig=figure('Name',s,'units','normalized');
    fig.ToolBar='none';
    fig.NumberTitle='off';
    fig.Name=[animal_ident ' - BOLD response matrix']; 
    fig.FileName=[animal_ident '-BOLDmatrix']; 
    fig_row=1;
    fig_col=1;
    fig_counter=1;

clims=[0 1]
colormap(jet)
        imagesc(abs(abs((DC_matrix_thres_sel_Events)-1)-0.2),clims)

ylim([-2 40])
xlim([-3 110])

%% end of file