function [all_reps_geo,rep_struct] = FI_size_normalization(Data, styler, err, norm_plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Written by Daniel Brink and Celina Tufvegren 2015-2016                              %
%                                                                                       %
%   Requirements: Knijnenburg morphology correction model:                              %
%                   Knijnenburg, T.A., Roda, O., Wan, Y., Nolan, G. P.,                 %                
%                   Aitchison, J.D., & Shmulevich, I. (2011).                           %
%                   A regression model approach to enable cell morphology               %
%                   correction in high-throughput flow cytometry.                       %
%                   Molecular systems biology, 7(1).                                    %
%                                                                                       %
%   Input:  Data = The filepath to your experiment data. Example:                       % 
%               '\Data\microplate' to get to the microplate data                        %
%               The data is required to be in Folders according to                      %
%               the following structure:                                                %
%               \Strains\Replicates\Conditions\.fcs files                               %
%                                                                                       %
%           styler = Choose the plotting style; barplots: 'bar',                        %
%               or a scatter plot FI vs time: 'time'.                                   %
%                                                                                       %
%            err = Whether or not you would like to work with                           %               
%               the replicates. If err = 1 you will get two outputs:                    %
%               all_reps_geo, which is the same as with err = 0, and                    %
%               reps_struct that has all the the means for the replicates.              %         
%               If err = 1 the file will also produce one                               %      
%               figure per strain, in the mode you've chosen ('bar' or 'time'),         %
%               showing the means of the replicates as well as error bars               % 
%               (=+/- standard deviation of the replicates)                             %
%                                                                                       %
%            norm_plot= Gives plots of the raw data histograms versus                   %
%                       the normalized histograms when norm_plot=1                      %
%                                                                                       %
%   Output: all_reps_data = a structure with each strain, replicate, and                %
%               condition. all_reps_data(n).mean will give you a matrix with            % 
%               [times means] for each experiment. Unless similarity == 1;              %
%               rep_struct is only applicable if err = 1. If so, the struct             % 
%               will show the means for the replicates at each condition and            % 
%               for each strain.                                                        %
%                                                                                       %
%                                                                                       %
%   Last altered: 2016-05-04 by Daniel Brink                                            %
%                                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION
warning off;
close all;
format short g  

%In-arguments
if nargin < 4
    norm_plot=0;
    if nargin < 3
        err = 0; 
        rep_struct = [];
        if nargin < 2
            styler = 'bar';
        end
    end
end


% Init. variables
home=pwd;% Your current working directory
cd (home); addpath(genpath(pwd));
cd ([home Data])            %CHECK WHAT THIS DOES

fig_count=1;
norm_plot_fig_count=0; norm_plot_fig_count_start=10;
simi = 0;
all_reps_geo = struct;
dynamic_y_axis=[];
folders = dir(pwd);

for k = length(folders):-1:1
    % remove non-folders from count
    if ~folders(k).isdir
        folders(k) = [ ];
        continue
    end

    % remove folders starting with . from count
    fname = folders(k).name;
    if fname(1) == '.'
        folders(k) = [ ];
    end
end

DirInfoS = folders;
s=length(folders);
g = 0;


%% BEGIN NORMALIZATION AND FIGURE PLOTTING PIPELINE
for is = 1:s
    
    strains{is,1} = DirInfoS(is).name;
    strain_string = strains{is};
    cd ([home Data '/' strain_string])
    folders = dir(pwd);

    for k = length(folders):-1:1
        % remove non-folders
        if ~folders(k).isdir
            folders(k) = [ ];
            continue
        end

        % remove folders starting with .
        fname = folders(k).name;
        if fname(1) == '.'
            folders(k) = [ ];
        end
    end
    DirInfoR = folders;
    rep = length(folders);
    clearvars replicates
    for ir = 1:rep
        geo_all_experiments = struct;
        times = {};
        replicates{ir,1} = DirInfoR(ir).name;
        rep_string = replicates{ir};
        cd ([home Data '/' strains{is,1} '/' rep_string])
        folders = dir(pwd);

        for k = length(folders):-1:1
            % remove non-folders
            if ~folders(k).isdir
                folders(k) = [ ];
                continue
            end

            % remove folders starting with .
            fname = folders(k).name;
            if fname(1) == '.'
                folders(k) = [ ];
            end
        end

        DirInfoC = folders;
        c = length(folders);
        clearvars conditions
%% NORMALIZE DATA        
        for ic = 1:c                
            conditions{ic,1} = DirInfoC(ic).name;
            condition_string = conditions{ic};
            cd([home Data '/' strain_string '/' rep_string '/' condition_string '/']);
            folder=[];
            skip_count=0; skip_count2=0;

            [DATA,N,Z,FName,time]=normalizer_caller(strain_string,condition_string,rep,home); %call the subscripts that do the actual normalization
            %%REMEMBER THAT THE SUBSCRIPT HERE MUST ALSO BE POLISHED
            
            
%% EXTRACT DATA FROM THE NORMALIZED SET                               
             for j = 1:N;
              	geo(j,1)=time(j);   
               	geo(j,2)=geomean(abs(Z{1,j}));  %Calculate the geometrical mean of the normalized histograms 
             end
             geo_sorted=sortrows(geo,1);

             
             geo_all_experiments(ic).strain = strain_string;    % Struct that keeps the data and background information of all normalization runs
             geo_all_experiments(ic).replicate = rep_string;
             geo_all_experiments(ic).condition = condition_string;
             geo_all_experiments(ic).mean = geo_sorted;
             
             simi = simi+1;                                     % Code block for the similarity subroutine
             all_reps_geo(simi).strain = strain_string;
             all_reps_geo(simi).replicate = rep_string;
             all_reps_geo(simi).condition = condition_string;
             all_reps_geo(simi).mean = geo_sorted; 
             all_reps_geo(simi).rep = rep;

             cd(home)
                    
%% HISTOGRAM NORMALIZATION PLOTS (code block from Knijnenburg et al.(2011))
                    % Plots histogram normalization plots if they were 
                    % chosen when calling the script (norm_plot=1)
                    
                    %NB! THIS HAS BEEN WRITTEN TO ADD THE 24H SAMPLE LAST.
                    %ONLY APPLIES TO DATASETS WITH 1,2,3,4,5,6,7,8,24H
             if norm_plot==1        
                XL = linspace(0,20000,length(DATA{1})/250);
                if norm_plot_fig_count>0
                    norm_plot_fig_count=norm_plot_fig_count;
                else
                    norm_plot_fig_count=fig_count+10;
                end
                                            
                h = figure(norm_plot_fig_count);
                set(h,'Color',[1 1 1]);
                for n = 1:N;
                    if n<=6
                        sub_pos=n;
                        mygca_sub_pos=n;
                        if any(regexp(FName{n},'24h.fcs$'))
                            skip_pos=n;
                            skip_count=1;
                            skip_count2=1;                                    
                            continue
                        end
                        if skip_count==1
                            sub_pos=n-1;
                            mygca_sub_pos=n-1;
                        end
                        figure(norm_plot_fig_count)
                        subplotfill(2,6,sub_pos,0.02,0.05);
                        bar(XL,histc(DATA{n}(:,3),XL));
                        mygcaDATA(mygca_sub_pos) = gca;
                        th = title([FName{n} ' Raw']);grid
                        set(th,'Interpreter','none','FontSize',12)
                        ylim auto
                        xlabel('FI'); ylabel('Count');
                        subplotfill(2,6,sub_pos+6,0.02,0.05);
                        bar(XL,histc(Z{n},XL));
                        mygcaZ(mygca_sub_pos) = gca;
                        th = title([FName{n} ' Norm.']);grid
                        ylim auto
                        xlabel('FI'); ylabel('Count');
                        norm_plot_max_y_DATA(n)=max(histc(DATA{n}(:,3),XL));
                        norm_plot_max_y_Z(n)=max(histc(Z{n},XL));
                        norm_plot_max_y_total=max([norm_plot_max_y_DATA norm_plot_max_y_Z])+0.1*max([norm_plot_max_y_DATA norm_plot_max_y_Z]);
                    end

                    if n>6 %use loop count not N-size, thus the loop enters this if-statement after it has done 6 in the above statement
                        sub_pos2=n;
                        mygca_sub_pos2=n;
                        if skip_count2==1
                            sub_pos2=n-1;
                            mygca_sub_pos2=n-1;
                        end
                        figure(norm_plot_fig_count+1)
                        subplotfill(2,6,sub_pos2-6+1,0.02,0.05); %have to be to n-6
                        bar(XL,histc(DATA{n}(:,3),XL));  %have to be to n
                        mygcaDATA(mygca_sub_pos2) = gca;
                        th = title([FName{n} ' Raw']);grid  
                        set(th,'Interpreter','none','FontSize',12)
                        ylim auto
                        xlabel('FI'); ylabel('Count');
                        subplotfill(2,6,sub_pos2+1,0.02,0.05);    %have to be to n
                        bar(XL,histc(Z{n-6},XL));        %have to be to n-6
                        mygcaZ(mygca_sub_pos2) = gca;
                        th = title([FName{n} ' Norm.']);grid
                        ylim auto
                        xlabel('FI'); ylabel('Count');
                        norm_plot_max_y_DATA(n)=max(histc(DATA{n}(:,3),XL));
                        norm_plot_max_y_Z(n)=max(histc(Z{n},XL));
                        norm_plot_max_y_total=max([norm_plot_max_y_DATA norm_plot_max_y_Z])+0.1*max([norm_plot_max_y_DATA norm_plot_max_y_Z]);
                    end
                            
                    if n==10
                        figure(norm_plot_fig_count+1)
                        subplotfill(2,6,sub_pos2-6+2,0.02,0.05); %have to change this to n-6
                        bar(XL,histc(DATA{skip_pos}(:,3),XL));  %have to change this to n-6
                        mygcaDATA(n) = gca;
                        th = title([FName{skip_pos} ' Raw']);grid                            
                        set(th,'Interpreter','none','FontSize',12)
                        ylim auto
                        xlabel('FI')
                        ylabel('Count')
                        subplotfill(2,6,sub_pos2+2,0.02,0.05);    %have to change this to n
                        bar(XL,histc(Z{skip_pos},XL));        %have to change this to n-6
                        mygcaZ(n) = gca;
                        th = title([FName{skip_pos} ' Norm.']);grid
                        ylim auto
                        xlabel('FI'); ylabel('Count');
                        norm_plot_max_y_DATA(n)=max(histc(DATA{n}(:,3),XL));
                        norm_plot_max_y_Z(n)=max(histc(Z{n},XL));
                        norm_plot_max_y_total=max([norm_plot_max_y_DATA norm_plot_max_y_Z])+0.1*max([norm_plot_max_y_DATA norm_plot_max_y_Z]);
                    end
                 end
                 norm_plot_fig_count=norm_plot_fig_count+2;
                 for n2=1:length(mygcaZ)
                    set(mygcaZ(n2), 'Ylim', [0, norm_plot_max_y_total])
                    set(mygcaDATA(n2), 'Ylim', [0, norm_plot_max_y_total])
                 end
             end

        end 
            
%% PLOTS WITHOUT ERRORBARS
            % Plots histogram normalization plots if they were 
            % chosen when calling the script (err=0),
            % plot the normalized data with one replicate per plot
            if err == 0
                for barcount=1:c   
                    geo_vec = geo_all_experiments(barcount).mean;                            
                    ybar(barcount,:) = geo_vec(:,2);
                    time_legend(1) = {num2str(geo_vec(1,1))};
                    for tcount = 2:length(geo_vec(:,2))
                        time_legend = [time_legend num2str(geo_vec(tcount,1))];
                    end
                end
                [n,m] = size(ybar);
                if n == 1
                    ybar = [ybar; zeros(1,m)];
                end
                
                figure(fig_count)

                if strcmp(styler,'bar') == 1
                    dynamic_y_axis=[dynamic_y_axis,max(max(ybar))];
                    bar(ybar)
                    for labelcount=1:c
                        Labels{labelcount}=geo_all_experiments(labelcount).condition;
                    end
                    set(gca, 'XTick', 1:c, 'XTickLabel', Labels);
                    title(['Strain TMB' all_reps_geo(simi).strain ' ' all_reps_geo(simi).replicate])
                    ylim([0,50000])
                    legend(time_legend)

                    elseif strcmp(styler,'time') == 1                           
                        for timecount = 1:c
                            xtime(:,timecount) = all_reps_geo(simi+timecount-c).mean(:,1);
                            ytime(:,timecount) = all_reps_geo(simi+timecount-c).mean(:,2); 
                            dynamic_y_axis=[dynamic_y_axis,max(max(geo_vec))];
                            if timecount ==1
                                cond_legend(1) = {all_reps_geo(1).condition}; 
                            else                                    
                                cond_legend = [cond_legend {all_reps_geo(timecount).condition}];
                            end
                        end
                    plot(xtime,ytime,'-v')
                    title(['Strain TMB' all_reps_geo(simi).strain ' ' all_reps_geo(simi).replicate])
                    ylim([0,50000]); xlabel('Time (h)')
                    legend(cond_legend)
                end
            end
                        
            clearvars xtime; clearvars ytime;
%% PLOTS WITH ERRORBARS
            % Plots histogram normalization plots if they were 
            % chosen when calling the script (err=1),
            % Gives the mean of the replicates with standard deviation 
            % in one plot per strain
            
            if err == 1
                if ir == rep
                    rep_struct(is).strain = strain_string;
                    rep_struct(is).time = geo_sorted(:,1);
                                
                    for icrep = 1:c 
                        for irep = 1:ir
                            gg = g+(irep-1)*c+icrep;
                                value(:,irep) = all_reps_geo(gg).mean(:,2);
                                cond_name = strrep(all_reps_geo(icrep).condition,' ','');
                                rep_struct(1,is).condition(1,icrep).condition = struct(cond_name,value);
                        end
                    end
                    g = gg;

                    for barcount=1:c                          
                        condition_str(barcount) = fieldnames(rep_struct(is).condition(barcount).condition);
                        geo_vec = cell2mat(struct2cell(rep_struct(1,is).condition(1,barcount).condition));
                        ybar(barcount,:) = mean(geo_vec,2);
                        y_error(barcount,:) = std(geo_vec,0,2);
                    end
                    time_legend(1) = {num2str(rep_struct(is).time(1))};
                    for tcount = 2:length(rep_struct(is).time)
                        time_legend = [time_legend num2str(rep_struct(is).time(tcount))];
                    end

                    [n,m] = size(ybar);
                    if n == 1
                        ybar = [ybar; zeros(1,m)];
                        y_error = [y_error; zeros(1,m)];
                    end

                    figure(fig_count)

                    if strcmp(styler,'bar') == 1
                        for labelcount=1:c
                            Labels{labelcount}=geo_all_experiments(labelcount).condition;
                        end
                        dynamic_y_axis=[dynamic_y_axis,max(max(ybar))];
                        h = bar(ybar);
                        set(h,'BarWidth',1);    % The side of the bars will now touch each other
                        set(gca,'YGrid','on')
                        set(gca,'GridLineStyle','-')
                        set(get(gca,'YLabel'),'String','Fluorescence','FontWeight','bold')

                        hold on;
                        numgroups = size(ybar, 1); 
                        numbars = size(ybar, 2); 
                        groupwidth = min(0.8, numbars/(numbars+1.5));
                        for i = 1:numbars
                            % Code block based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
                            x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
                            errorbar(x, ybar(:,i), y_error(:,i), 'k', 'linestyle', 'none'); %PLOTS THE ERRORBAR GRAPH
                        end
                        legend(time_legend)%set(time_legend,'Location','BestOutside','Orientation','horizontal')

                        set(gca, 'XTick', 1:c, 'XTickLabel', Labels);
                        title(['Strain TMB' all_reps_geo(simi).strain ' n = ' num2str(rep)])
                        ylim([0,50000])
                        legend(time_legend)
                        clearvars time_legend
                        
                        elseif strcmp(styler,'time') == 1                           
                            for timecount = 1:c
                                geo_vec = cell2mat(struct2cell(rep_struct(1,is).condition(1,timecount).condition));
                                xtime(:,timecount) = rep_struct(is).time;
                                ytime(:,timecount) = mean(geo_vec,2); 
                                ytime_error(:,timecount) = std(geo_vec,0,2);
                                dynamic_y_axis=[dynamic_y_axis,max(max(geo_vec))];
                                if timecount == 1
                                    cond_legend(1) = condition_str(1); 
                                else                                    
                                    cond_legend = [cond_legend condition_str(timecount)];
                                end
                                errorbar(xtime,ytime,ytime_error,'-v');  %PLOTS THE ERRORBAR GRAPH
                                %hold on    
                                title(['Strain TMB' all_reps_geo(simi).strain ' n = ' num2str(rep)])
                                xlim([min(rep_struct(1,is).time),max(rep_struct(1,is).time)])
                                xlabel('Time (h)')
                                ylabel('Mean Normalized GFP FI')
                                legend(cond_legend)  
                                        
                            end 
                            clearvars cond_legend
                    end
                end
            end
            
            fig_count=fig_count+1;
            clearvars geo geo_all_experiments; clearvars ybar y_error cond_legend;
            clearvars ytime xtime ytime_error time_legend; clearvars value geo_vec;
    end

%% DYNAMIC Y-AXIS
        % Sets the y-axes of all plots (except norm plot histograms) to
        % display the same max value, based on the max FI-value in the
        % normalized dataset
        fig_handles = get(0,'Children');
        for dynamic_y_count=1:length(fig_handles)
            if fig_handles(dynamic_y_count).Number < norm_plot_fig_count_start
                figure(fig_handles(dynamic_y_count).Number)
                ylim([0,max(dynamic_y_axis)+2000]); 
            end
        end
        
%% NORMALIZATION LOOP END
end


function [DATA,N,Z,FName,time]=normalizer_caller(strain_string,condition_string,replicates,home)

% This file is modified by D. Brink from files provided by:
% Knijnenburg, T. A., Roda, O., Wan, Y., Nolan, G. P., Aitchison, J. D., & Shmulevich, I. (2011). 
% A regression model approach to enable cell morphology correction in high?throughput flow cytometry. 
% Molecular systems biology, 7(1), 531.

%% Load data - a bunch of fcs files
folder = [pwd '/']; %for UNIX/Mac OS
D = dir(folder);
p = 0;
FData = cell(0);
FName = cell(0);
for n = 1:length(D); 
    if length(D(n).name)>4;
        if strcmpi(D(n).name(end-3:end),'.fcs')
            
            time(n-2,:) = str2num(D(n).name(end-6:end-5));  %extracts the time from the filename
            
            p = p + 1;
            [Data, Header] = fca_readfcs([folder D(n).name]);
            Pars = cell(Header.NumOfPar,1);
            for j = 1:Header.NumOfPar;
                Pars{j} = Header.par(j).name;
            end
            FSC = strmatch('FSC',Pars);
            SSC = strmatch('SSC',Pars);
            FL = strmatch('FL1',Pars);
            FData{p} = Data(:,[FSC(1) SSC(1) FL(1)]);
            FName{p} = Header.filename;
        end
    end
end
DATA = FData;  
%Each cell in DATA contains a matrix with three columns, column 1 = FSC,
%column 2 = SSC, column 3 = FL, the number of rows is the number of cells (or events) that are measured 

%% Number of experiments
N = length(DATA);

%% Step 1 Preprocessing according to Newman, 2006, Nature
[DATA,s,val1,val2,val3] = FC_preprocess(DATA);

%% Step 2 Determining dense scatter area CONTINUOUS
[y,c] = selectdensescatterarea_cont(DATA,val1,val2);

%% Step 3 Regression computing fluorescence compensated for scatters
Z = doregress_constrained(DATA,y,c,val1,val2,'lin');
