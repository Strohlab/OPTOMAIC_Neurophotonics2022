function [Slwaves_arr, Slwaves_arr_sel,err]=DC_analysis_SlWaves(DC_CaSignal_dsampled,DC_Carec_OnOffsets,DC_VisStim_onsets,g_DC_Spike2_SamplingRate,g_SW_DetectWindow,datum,prot_fid,analyse_animalfolder,animal_ident,g_DC_ProtImages);
%% 
s=['- \n']; fprintf(prot_fid,s);disp(s)
analyse_function_name_version='FUNCTION DC_analysis_Slwaves V20220519';
analyse_function_author='Dirk Cleppien';
s=['### (' analyse_function_name_version ' - ' analyse_function_author ')  \n']; fprintf(prot_fid,s);disp(s)

%% Global variables
g_SW_DetectWindow
g_DC_Spike2_SamplingRate

%%
DC_h_max=5; % max. number of depicted number of slow waves per stim interval

%% Initialization of analysis figure
s='Slow wave analysis'
    fig=figure('Name',s,'units','normalized','outerposition',[0 0 1 1]);
        fig.ToolBar='none';
    fig.NumberTitle='off';
    fig.Name=[animal_ident ' - Slow wave analysis']; 
    fig.FileName=[animal_ident '-SlowWaveAnalysis']; 
    fig_row=2;
    fig_col=3;
    fig_counter=1;

%% slow wave array
Slwaves_arr=zeros(size(DC_Carec_OnOffsets,2),5);
size(Slwaves_arr)
% 1: corresponding stim interval
% 2: response time after stim 
% 3: interval between current and previous slwave onset
% 4: interval between current and previous slwave offset
% 5: number of slow waves per stim interval

%% Identification and fill in slow wave array
zj=1;
zn=1;
for zi=2:size(Slwaves_arr,1),
    % stim interval counter zj
        d_sw=DC_Carec_OnOffsets(1,zi);
       % not all vis stim has a detected response!
    while (d_sw>DC_VisStim_onsets(zj+1)), zj=zj+1; disp(['While ' num2str(zj)]);
     zn=1;   
    end
    % #1: Slwave is in zj stim interval
        Slwaves_arr(zi,1)=zj;    
    % #2: response time on stimulus    
        d_vis=DC_VisStim_onsets(zj);
        diff_vis_sw=d_sw-d_vis;
        Slwaves_arr(zi,2)=diff_vis_sw;
    % #3: time interval before slwave (ONset(zi-1) to onset(zi))
        d_sw_prev=DC_Carec_OnOffsets(1,zi-1);
        diff_sw=d_sw-d_sw_prev;
        Slwaves_arr(zi,3)=diff_sw;
    % #4: time interval before slwave (OFFset(zi-1) to onset(zi))
        d_sw_prev_on=DC_Carec_OnOffsets(2,zi-1);
        diff_sw_on=d_sw-d_sw_prev_on;
        Slwaves_arr(zi,4)=diff_sw_on;
    % #5: number of slow wave per interval
        Slwaves_arr(zi,5)=zn;
        zn=zn+1;

     disp([num2str(zi) ' ' num2str(zj) ' ' num2str(d_sw) ' - ' num2str(d_vis) ... 
        ' = ' num2str(diff_vis_sw) ' ; ' ...
        num2str(d_sw) ' - ' num2str(d_sw_prev) ' = ' num2str(diff_sw) ' ; ' ...
        num2str(d_sw) ' - ' num2str(d_sw_prev_on) ' = ' num2str(diff_sw_on) ' ; ' ...
        'Slwaves_arr: ' num2str(Slwaves_arr(zi,:))]);
end

%% analysis figure
subplot(fig_row,fig_col,fig_counter)
d=sort(Slwaves_arr(:,3))
h=plot(1:40,d(1:40),'r-','LineWidth',5)
hold on
h=plot(41:size(d,1),d(41:end),'g-','LineWidth',5)
hold off
ylim([0 2.5*10^4])
    title(['Fig.' num2str(fig_counter) ':Interval between adjacent slow waves']);
    ylabel('time')
    xlabel('number slow wave')
    xlim([-10 150])
    fig_counter=fig_counter+1;

%% Number of slow waves per stim interval
subplot(fig_row,fig_col,fig_counter)
    h=histogram(Slwaves_arr(:,5),DC_h_max);
    title({['Fig.' num2str(fig_counter) ': Number of slow waves per stimulation interval'];['Number of stimulation intervals = ' num2str(max(Slwaves_arr(:,1)))]});
    ylabel('Number slow waves')
    xlabel('Number per stimulation interval')
    xlim([0 DC_h_max]);
    fig_counter=fig_counter+1;

    %% Choice: in time window (30-300ms) after stim + long time interval to previous slow wave
d1_arr=squeeze(Slwaves_arr(:,2)<g_SW_DetectWindow(1)); 
% 2: response time after stim; in stim time window
d2_arr=squeeze(Slwaves_arr(:,3)>g_SW_DetectWindow(2));
% 3: interval between current and previous slwave ONSET
d3_arr=squeeze(Slwaves_arr(:,4)>g_SW_DetectWindow(2));
% 4: interval between current and previous slwave OFFSET
%
d_res_on=d1_arr.*d2_arr;
sum(d_res_on)
d_res_off=d1_arr.*d3_arr;
disp(['Number of locked slow waves longer than ' num2str(g_SW_DetectWindow(2)/2) ...
    ' ms = ' num2str(sum(d_res_on)) '(onset to onset) --- = ' num2str(sum(d_res_off)) ' (offset to onset)']);

%% Selection of slow waves
%Slwaves_arr_sel
d=horzcat(Slwaves_arr,d_res_on,d_res_off);
Slwaves_arr_sel=d;  % assign selected slow waves to main program variable
size(Slwaves_arr_sel)
% 1: corresponding stim interval
% 2: response time after stim 
% 3: interval between current and previous slwave onset
% 4: interval between current and previous slwave offset
% 5: number of slow waves per stim interval
% 6: onset previous slow wave longer than g_SW_DetectWindow(2)
% 7: offset previous slow wave longer than g_SW_DetectWindow(2)
%% Depiction in analysis figure
subplot(fig_row,fig_col,fig_counter)
dummy=Slwaves_arr_sel(:,6);
dummy=dummy(:)+Slwaves_arr_sel(:,7)*2;
imagesc(dummy');
colorbar
    title({['Fig.' num2str(fig_counter) ': Number of stim. slow waves with inter slow wave interval >' num2str(g_SW_DetectWindow(2)/g_DC_Spike2_SamplingRate) 's'];['by Onsets=' num2str(sum(Slwaves_arr_sel(:,6))) '; by Offsets= ' num2str(sum(Slwaves_arr_sel(:,7))) ' of all ' num2str(size(Slwaves_arr_sel(:,7),1))]});
    ylabel('[arbitrary units]')
    xlabel('Number slow waves')
    fig_counter=fig_counter+1;

%% Slow waves sorted by inter slow wave intervals and depict in analysis figure
DC_InterSlWavesInterval=[0 3 5 7 9 10 20]
anz_DC_InterSlWavesInterval=size(DC_InterSlWavesInterval,2)
dummy=zeros(anz_DC_InterSlWavesInterval-1,size(Slwaves_arr_sel,1));
size(dummy)
for zi=1:anz_DC_InterSlWavesInterval-1,
    d1=squeeze(Slwaves_arr(:,3)>(g_DC_Spike2_SamplingRate*DC_InterSlWavesInterval(zi)))*2;
    d2=squeeze(Slwaves_arr(:,3)<(g_DC_Spike2_SamplingRate*DC_InterSlWavesInterval(zi+1)));
    dummy(zi,:)=((d1-1)==d2)
end
%%
subplot(fig_row,fig_col,fig_counter)
    bar(sum(dummy'))
    ylim([0  max(Slwaves_arr(:,1))])
    title({['Fig.' num2str(fig_counter) ': Number of slow waves per inter slow wave interval'];['Intervals within: (' num2str(DC_InterSlWavesInterval) ')s']});
    ylabel('Number slow waves')
    xlabel('Inter slow wave interval')
    fig_counter=fig_counter+1;
subplot(fig_row,fig_col,fig_counter)
    imagesc((dummy))
    title({['Fig.' num2str(fig_counter) ': Distribution of slow waves sorted by inter slow wave interval'];['Intervals within: (' num2str(DC_InterSlWavesInterval) ')s']});
    ylabel('Inter slow wave interval')
    xlabel('Number slow wave')
    fig_counter=fig_counter+1;

%% function tail
err=0;
image_suffix='Subroutine_Analysis_Slwaves'
old=cd(['../data/' analyse_animalfolder '/']);
    DC_hg_name=[ animal_ident '_' datum '_' image_suffix '.emf' ];
saveas(fig,DC_hg_name,'emf')
cd(old)

%% end of function
s=['### (' analyse_function_name_version ') - end \n']; fprintf(prot_fid,s);disp(s)
