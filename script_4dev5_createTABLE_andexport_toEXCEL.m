%% 
% last modified: 20200606

clear all
close all

%% 
do_plot = 1
do_save = 0


%%
%cd('F:\=MEDEA_backup20191122\=MEDEA_BOING_eyetracker')
current_path = pwd;
pathsep = '\';
regexp(current_path, pathsep,'split');   
current_path = ans;

drive_dir = current_path{1,1};
exp_folder = fullfile(drive_dir, current_path{1,2}, current_path{1,3});
data_dir = fullfile(exp_folder, 'DATA_PILOT');
%video_dir = fullfile(exp_folder, 'VIDEO_def_trigger_code(20190703)');
disp(data_dir)


%% LIST  DATA pilot (only ADult data) in matlab conversion
cd(data_dir)
file_list = dir('*.xlsx')
file_todel = [];
for jj = 1:length(file_list)
    if file_list(jj).name(1:2) == '~$'
        file_todel = [ file_todel; jj ];
    elseif contains(file_list(jj).name, 'bam')
        file_todel = [ file_todel; jj ];
    end
end
disp(file_todel)
file_list(file_todel) = [];


%% VIDEO NAME
video_name_list = { 
                      '110.avi'; '111.avi'; '112.avi'; '113.avi'; '114.avi'; ...
                      '120.avi'; '121.avi'; '122.avi'; '123.avi'; '124.avi'; ...
                      '210.avi'; '211.avi'; '212.avi'; '213.avi'; '214.avi'; ...
                      '220.avi';
                      %'221.avi'
                      '221(to check).avi'
                      '222.avi'
                      '223.avi'
                      '224.avi'
                     %'215_control(60sec).avi'
                     %'225_control(60sec).avi'
                     %'115_control(60sec).avi'
                     %'125_control(60sec).avi'
                     };
%video_name_list = {'221(to check).avi'}
%video_name_list = {'221.avi'}
%video_name_list = {'111.avi'}


%% CYCLE FOR ALL THE VIDEO (??)
%for ii = 1:length(video_name_list)
    ii = 2;
    disp(video_name_list{ii})
    
    i_video = str2num(video_name_list{ii}(1:3))
    
    % LOAD DOT position in one video (as table_video)
    load([ 'table_video_' num2str(i_video) ])
    disp(table_video(1:3,:))
    
    video_duration_frame = size(table_video,1)
    video_sample_rate = 60 %Hz
    video_msec_xframe = 1/video_sample_rate*1000
    
    video_duration_msec = video_msec_xframe * video_duration_frame 
    
    dot1_ydown = table_video.dot1_ydown;
    dot1_yup = table_video.dot1_yup;
    dot2_ydown = table_video.dot2_ydown;
    dot2_yup = table_video.dot2_yup;
   
    figure; plot(dot1_ydown)
            hold on; plot(dot1_yup)
            hold on; plot(dot2_ydown, '--');  plot(dot2_yup, '--') 
    bounce_down_frame  = find(dot1_yup == max(dot1_yup))
    bounce_up_frame = find(dot1_ydown == min(dot1_ydown))
   
    n_bounce_down = length(bounce_down_frame)
    n_bounce_up = length(bounce_up_frame)
    
    bounce_dot1_frame = sort([bounce_up_frame; bounce_down_frame])
    n_bounce = length(bounce_dot1_frame)
    
    bounce_dot2_frame = sort([ find(dot2_yup == max(dot2_yup)); ...
                               find(dot2_ydown == min(dot2_ydown)) ]);
    set(gca, 'Ydir', 'reverse')
    legend('dot1 Ydown', 'dot1 Yup', 'dot2 Ydown', 'dot2 Yup')
   
    
    %% LOAD all subjects data for 1 video
    %table_eyetrack = [];
    
    %initialize - - - - -
    media_width_height = []; media_width = []; media_height = [];
    screen_x_pixel = []; screen_y_pixel = [];
    time_vector_msec = []; fixation_order = [];
    gaze_event_type = []; gaze_event_msec = [];
    fixation_x = []; fixation_y = [];
    pupil_size_left = []; pupil_size_right = [];
    gaze_point_x = []; gaze_point_y = [];
    
    % !!! do not consider INFANT subj: 'Boing_Test1_01bam01.xlsx'
    % LOAD 
    group_struct = [];
    i_subj = 1; subj_excluded = 0
    %for i_subj = 2:length(file_list)
    for ii = 2:length(file_list)
        subj_name = file_list(ii).name
        subj_name = subj_name(1:end-5)
        
        %subj_tag = [ 's' num2str(i_subj) ]
        
        % LOAD - - - - -
        load(subj_name)
        video_numb = [trial_struct(:).video_numb];
        
        row_idx = find(video_numb == i_video)
        %ismember(trial_struct.video_numb, 110)
        
        if trial_struct(row_idx).monitor_XYpixel(1) == 1920 &&...
            trial_struct(row_idx).monitor_XYpixel(2) == 1080
            
            group_struct(i_subj).subj_name              = subj_name
            group_struct(i_subj).video_numb             = trial_struct(row_idx).video_numb;
            group_struct(i_subj).monitor_XYpixel        = trial_struct(row_idx).monitor_XYpixel;  
            group_struct(i_subj).media_width_height     = trial_struct(row_idx).media_width_height; 
            group_struct(i_subj).media_xy_origin        = trial_struct(row_idx).media_xy_origin;
            group_struct(i_subj).n_sample       = length(trial_struct(row_idx).time_vector_msec);
            group_struct(i_subj).time_vector_msec       = trial_struct(row_idx).time_vector_msec;
            group_struct(i_subj).trial_duration_msec   = trial_struct(row_idx).time_vector_microsec(end);
            group_struct(i_subj).gaze_event_type        = trial_struct(row_idx).gaze_event_type;
            group_struct(i_subj).gaze_event_msec        = trial_struct(row_idx).gaze_event_msec;
            group_struct(i_subj).saccade_idx            = trial_struct(row_idx).saccade_idx;         
            group_struct(i_subj).fixation_order         = trial_struct(row_idx).fixation_order;
            group_struct(i_subj).fixationXY             = trial_struct(row_idx).fixationXY;   
            group_struct(i_subj).pupil_size             = trial_struct(row_idx).pupil_size;
            %group_struct(i_subj).    gaze_leftXY: [2190×2 double]
            % group_struct(i_subj).   gaze_rightXY: [2190×2 double]
            group_struct(i_subj).gazepointXY_MCS        = trial_struct(row_idx).gazepointXY_MCS
            %group_struct(i_subj).gazepointXY_ADCS: [2190×2 double]
        
            i_subj = i_subj +1;
        else
            subj_excluded = subj_excluded +1;
            disp(['!!! subject NOT INCLUDED'])
        end
    end
    disp(['!!! subject NOT INCLUDED = ' num2str(subj_excluded)])

    
    %% TIME VECTOR of different length 
    video_duration_msec 
    video_duration_frame*2
    
    n_sample_vector = []; 
    trial_duration_msec = [];
    for ii = 1:length(group_struct)
        n_sample_vector = [ n_sample_vector; group_struct(ii).n_sample ];
        trial_duration_msec = [ trial_duration_msec ; group_struct(ii).time_vector_msec(end) ];
    end
    n_sample_min = min(n_sample_vector)
    n_sample_max = max(n_sample_vector)
    
    if n_sample_min > video_duration_frame*2
        n_resample = video_duration_frame*2
    else
        n_resample = n_sample_min;
    end
    resample_vector = [];
    for ii = 1:length(group_struct)
        n_sample_diff = n_sample_vector(ii)-n_resample;
        sample_diff = round(linspace(1,n_sample_vector(ii),n_sample_diff));
        resample_vector_tmp = 1:n_sample_vector(ii);
        resample_vector_tmp(sample_diff) = [];
        
        resample_vector(ii,:) = resample_vector_tmp;
    end
    
%     %n_sample_mean = round(n_sample_max - (n_sample_max - n_sample_min)/2)
%     sample_vector = []
%     for ii = 1:length(group_struct)
%         %ii = 4
%        %tmp = n_sample_max(ii) - n_sample_vector(ii)
%        tmp = n_sample_vector(ii) - n_sample_min
%        if tmp == 0; 
%            sample_vector(ii,:) = 1:n_sample_vector(ii);
%            time_vector_msec = group_struct(ii).time_vector_msec;
%        else
%            sample_vector(ii,:) = ceil(tmp/2) : (n_sample_vector(ii) -floor(tmp/2)-1);
%        end
%     end

% DUPLICATE duration of video to match the sampling rate of eyetracker
video_table_ext = zeros(video_duration_frame*2,8);
video_table_ext(1:2:end) = [ table_video.dot1_xleft, table_video.dot1_xright, ...
                    table_video.dot1_yup, table_video.dot1_ydown, ...
                     table_video.dot2_xleft, table_video.dot2_xright, ...
                    table_video.dot2_yup, table_video.dot2_ydown ];
video_table_ext(2:2:end) = [ table_video.dot1_xleft, table_video.dot1_xright, ...
                    table_video.dot1_yup, table_video.dot1_ydown, ...
                     table_video.dot2_xleft, table_video.dot2_xright, ...
                    table_video.dot2_yup, table_video.dot2_ydown ];

                
    %% create VIDEO_table
    fixation_XY_allsubj = [];
    gaze_XY_allsubj = [];
    for iii = 1:length(group_struct)
        sample_vector_tmp = resample_vector(iii,:);
        
        % FIXATION mean
        fixation_XY_allsubj = [ fixation_XY_allsubj,  group_struct(iii).fixationXY(sample_vector_tmp,1), ...
                                                group_struct(iii).fixationXY(sample_vector_tmp,2) ];
        % GAZE POINT mean
        gaze_Xtmp = []; gaze_Ytmp = []; %group_struct(iii).gazepointXY_MCS;
        if isstring(group_struct(iii).gazepointXY_MCS(1))
           for jj=1:size(group_struct(iii).gazepointXY_MCS,1)
               if strcmp(group_struct(iii).gazepointXY_MCS(jj,1), "")
                    gaze_Xtmp(jj,1) = NaN; gaze_Ytmp(jj,1) = NaN;
               else
                    gaze_Xtmp(jj,1) = str2num(group_struct(iii).gazepointXY_MCS(jj,1));
                    gaze_Ytmp(jj,1) = str2num(group_struct(iii).gazepointXY_MCS(jj,2));
               end
           end
        else
            gaze_Xtmp = group_struct(iii).gazepointXY_MCS(:,1);
            gaze_Ytmp = group_struct(iii).gazepointXY_MCS(:,2);
        end
        
        gaze_XY_allsubj = [ gaze_XY_allsubj  gaze_Xtmp(sample_vector_tmp)  gaze_Ytmp(sample_vector_tmp) ];
    
        %PUPIL SIZE
    end 
    
    %fixation_XY_table = [ video_table_ext, fixation_XY_allsubj ];  
    trial_duration_tmp = 1:video_duration_frame*2';
    fixation_XY_table = [ trial_duration_tmp' ,video_table_ext, fixation_XY_allsubj ];

    table_header = {'frame', 'dot1'}
    T = table(fixation_XY_table(:,2:3));
    T.Properties.VariableNames = {'frame', 'dot1'};
    
    %%    
    
    
    LastName = {'Smith';'Johnson';'Williams';'Jones';'Brown'};
Age = [38;43;38;40;49];
Height = [71;69;64;67;64];
Weight = [176;163;131;133;119];
BloodPressure = [124 93; 109 77; 125 83; 117 75; 122 80];
T = table(Age,Height,Weight,BloodPressure,...
    'RowNames',LastName)
writetable(T,'cond111_fixation.xls') %,'FileType', 'spreadsheet',1)
    
    % SAVE in EXCEL
    xls_name = [ 'cond' num2str(i_video) '_fixation' ]
    xlswrite(xls_name,fixation_XY_table);
    
    T = 
    % PLOT
    %figure; plot(fixation_XY_allsubj(:,1:2:end))
    %hold on; plot(nanmean(fixation_XY_allsubj(:,1:2:end),2))
    figure; plot(gaze_XY_allsubj(:,1:2:end))
    hold on; plot(nanmean(gaze_XY_allsubj(:,1:2:end),2))
    
        %monitor_XYpixel
        media_width_height(i_subj,:) = trial_struct(row_idx).media_width_height;
        media_width(i_subj,:) = trial_struct(row_idx).media_width_height(1);
        media_height(i_subj,:) = trial_struct(row_idx).media_width_height(2);
        %media_xy_origin(i_subj,:) = trial_struct(row_idx).media_xy_origin;
   
        screen_x_pixel(i_subj,:) = trial_struct(row_idx).monitor_XYpixel(1);
        screen_y_pixel(i_subj,:) = trial_struct(row_idx).monitor_XYpixel(2);
        
        
        time_vector_msec(i_subj,:) = trial_struct(row_idx).time_vector_msec;
        eyetrack_msec_xframe = (time_vector_msec(end) / length(time_vector_msec))
        eyetrack_sample_rate = round(1/eyetrack_msec_xframe *1000)
        trial_duration_msec = time_vector_msec(i_subj,end)
        
        gaze_event_type(i_subj,:) = trial_struct(row_idx).gaze_event_type;
        gaze_event_msec(i_subj,:) = trial_struct(row_idx).gaze_event_msec;
   
        %saccade_idx(i_subj,:) =  trial_struct(row_idx).saccade_idx;
        fixation_order(i_subj,:) =  trial_struct(row_idx).fixation_order;
        
        fixation_x(i_subj,:) = trial_struct(row_idx).fixationXY(:,1);
        fixation_y(i_subj,:) = trial_struct(row_idx).fixationXY(:,2);
    
        pupil_size = trial_struct(row_idx).pupil_size;
        pupil_size_left(i_subj,:) = pupil_size(:,1);
        pupil_size_right(i_subj,:) = pupil_size(:,2);
        
        %replace -1 with NAN
        pupil_size_left(i_subj,find(pupil_size_left(i_subj,:) == -1)) = NaN;
        pupil_size_right(i_subj,find(pupil_size_right(i_subj,:) == -1)) = NaN;
        
        gaze_left_tmp = trial_struct(row_idx).gaze_leftXY;
        gaze_right_tmp = trial_struct(row_idx).gaze_rightXY;
        
        gaze_point_tmp = trial_struct(row_idx).gazepointXY_MCS;
        gaze_point_x(i_subj,:) = gaze_point_tmp(:,1);
        gaze_point_y(i_subj,:) = gaze_point_tmp(:,2);
%         for cc = 1:length(gaze_point_tmp)
%             if strcmp(gaze_point_tmp(cc,1), ""); 
%                 gaze_point_x(cc) = NaN; gaze_point_y(cc) = NaN;
%             elseif strcmp(gaze_point_tmp(cc,1), 'NaN');
%                 gaze_point_x(cc) = str2num(gaze_point_tmp(cc,1));
%                  gaze_point_y(cc) = str2num(gaze_point_tmp(cc,2));
%             else
%                 gaze_point_x(cc) = gaze_point_tmp(cc,1);
%                  gaze_point_y(cc) = gaze_point_tmp(cc,2);
%             end
%         end
    end
    
    
    %% 1-PLOT fixation_position over the frame
    % with a sliding bar
    
    %figure; plot(~isnan(fixation_x)', fixation_y')
    figure; scatter(fixation_x', fixation_y')
    figure; subplot 311; plot(fixation_x)
            subplot 312; plot(fixation_y)
            subplot 313; hold on; 
            %plot(pupil_size_left); plot(pupil_size_right);
            scatter(1:length(pupil_size_left), pupil_size_left, 1);
            scatter(1:length(pupil_size_right), pupil_size_right, 1);
            %plot(pupil_size_right); 
    
    % - - - - - - - - -  - - - -
    % LOAD video frame_struct
    figure(10); 
    subplot_order = [ 9, 2, 11, 4, 13, 6, 15, 8 ]
    %if n_bounce_down = n_bounce_up +1
    for jj = 1:n_bounce -1
        idx_frame = bounce_dot1_frame(jj)
        idx_msec = idx_frame * frame_duration_msec
        
        [ tmp, eyetrack_idx_tmp ] = min(abs(time_vector_msec - idx_msec))
        %[ tmp, eyetrack_idx(jj) ] = min(abs(time_vector_msec - idx_msec))
        
        gaze_duration_tmp = gaze_event_msec(eyetrack_idx_tmp)
        gaze_type_tmp = gaze_event_type(eyetrack_idx_tmp)
        % 2 = fixation
        
        %fixation_idx = fixation_order(eyetrack_idx_tmp)
        gaze_duration_max = round(max(gaze_event_msec) / eyetrack_msec_xframe);
         
        for hh = 1:gaze_duration_max  %hh=2
            if gaze_event_msec(eyetrack_idx_tmp - hh) - gaze_event_msec(eyetrack_idx_tmp) ~= 0
                gaze_onset_idx = eyetrack_idx_tmp - hh
                break
                %hh = gaze_duration_max;
            end
        end
        for hh = 1:gaze_duration_max  %hh=2
            if gaze_event_msec(eyetrack_idx_tmp + hh) - gaze_event_msec(eyetrack_idx_tmp) ~= 0
                gaze_offset_idx = eyetrack_idx_tmp + hh
                break
            end
        end
        
        %gaze_offset_idx - gaze_onset_idx 
        % duration of fixation before the bounce
        fixation_prebounce_msec = (eyetrack_idx_tmp - gaze_onset_idx) * eyetrack_msec_xframe
        
        subplot(2,8,subplot_order(jj)); hold on
        %scatter(fixation_x(eyetrack_idx(jj)), fixation_y(eyetrack_idx(jj)))
        scatter(fixation_x(eyetrack_idx_tmp), fixation_y(eyetrack_idx_tmp), gaze_duration_tmp)
        scatter(fixation_x(eyetrack_idx_tmp), fixation_y(eyetrack_idx_tmp), fixation_prebounce_msec, 'filled')
        xlim([0 media_width]); ylim([0 media_height]);
        set(gca, 'YDir','reverse')
    end
    
    
    %% PLOT  HEAT-map  for each AREA of INTEREST (aoi)
    bounce_dot1_frame
    %bounce_dot2_frame
    
    % TIME before the bounce (of the leading dot)
    % = time after the lagging bounce
    interval_pre_msec = 100
    interval_pre_frame = round(interval_pre_msec / eyetrack_msec_xframe)
    interval_post_frame = interval_pre_frame;
    
    % TIME window to select before dot1 bounce
    % (and after dot2 bounce) 
    bounce_interval_videoframe = mean(bounce_dot2_frame - bounce_dot1_frame)
    bounce_interval_msec = bounce_interval_videoframe * video_msec_xframe
    bounce_interval_eyetrackframe = round(bounce_interval_msec/eyetrack_msec_xframe) 
    
    %area_ofinterest_x = 1:media_width;
    %area_ofinterest_ydown = media_height/3*2+1 : media_height;
    %area_ofinterest_yup = 1: media_height/3*1;
    
    area_ofinterest_ydown = screen_y_pixel/3*2;
    area_ofinterest_yup = screen_y_pixel/3;
    
    dot1_x_left = table_video.dot1_xleft(1) + trial_struct(1).media_xy_origin(1)
    dot1_x_right = table_video.dot1_xright(1) + trial_struct(1).media_xy_origin(1)
    dot2_x_left = table_video.dot2_xleft(1) + trial_struct(1).media_xy_origin(1)
    dot2_x_right = table_video.dot2_xright(1) + + trial_struct(1).media_xy_origin(1)
    
    for jj = 1:n_bounce -1
        if rem(jj,2) > 0; 
            title_tmp = 'bounce down '; bounce_down = 1
        else; 
            title_tmp = 'bounce up '; bounce_down = 0; 
        end
        
        idx_frame = bounce_dot1_frame(jj)
        idx_msec = idx_frame * video_msec_xframe
          
        [ tmp, eyetrack_idx_tmp ] = min(abs(time_vector_msec - idx_msec))
       
        window_vector_frame = [ (eyetrack_idx_tmp -  interval_pre_frame) : ...
                                (eyetrack_idx_tmp + bounce_interval_eyetrackframe + interval_post_frame) ];
        window_duration_msec = length(window_vector_frame) * eyetrack_msec_xframe
        
%         figure; subplot 312; plot(fixation_y(window_vector_frame))
%             hold on; plot(gaze_left_tmp(window_vector_frame,2))
%             hold on; plot(gaze_right_tmp(window_vector_frame,2))
% 
%             subplot 313;  plot(fixation_x(window_vector_frame)) 
%             hold on; plot(gaze_left_tmp(window_vector_frame,1))
%             hold on; plot(gaze_right_tmp(window_vector_frame,1))
% 
%             subplot 311;  plot(fixation_x(window_vector_frame)) 
%             hold on; plot(gaze_point_y(window_vector_frame))
            %
        heat_map = zeros(media_width, media_height);
        heat_map_allscreen = zeros(screen_x_pixel, screen_y_pixel);
        
        fixation_left_nan = []; fixation_left_out = [];
        fixation_right_nan = []; fixation_right_out = [];
        
        subplot_order = [ 9, 2, 11, 4, 13, 6, 15, 8 ]
    
        for ww = 1:length(window_vector_frame)
            %ww = 12
            i_frame = window_vector_frame(ww)  %for eyetracker
            
            % heat map ALL SCREEN size with gaze_point measure
            if ~isnan(gaze_left_tmp(i_frame,2)) 
                if bounce_down && gaze_left_tmp(i_frame,2) > area_ofinterest_ydown(1)
                        gaze_left_tmp(i_frame,2)
                    %if ismember(fixation_y(i_frame), area_ofinterest_ydown)  
                        heat_map_allscreen(gaze_left_tmp(i_frame,1), gaze_left_tmp(i_frame,2)) = ...
                                heat_map_allscreen(gaze_left_tmp(i_frame,1), gaze_left_tmp(i_frame,2)) +10;
                    %end
                elseif bounce_down == 0 && gaze_left_tmp(i_frame,2) < area_ofinterest_yup(end) 
                        gaze_left_tmp(i_frame,2)
                     %if ismember(fixation_y(i_subj,i_frame), area_ofinterest_yup)  
                        heat_map_allscreen(gaze_left_tmp(i_frame,1), gaze_left_tmp(i_frame,2)) = ...
                                heat_map_allscreen(gaze_left_tmp(i_frame,1), gaze_left_tmp(i_frame,2)) +1;
                     %end
                else
                    fixation_left_out(end+1) = i_frame
                end
                
                figure(4); subplot(2,8,subplot_order(jj))
                hold on; scatter(gaze_left_tmp(i_frame,1), gaze_left_tmp(i_frame,2))
                %heat_map_allscreen(fixation_x(i_subj,i_frame), fixation_y(i_subj,i_frame))
            else
                fixation_left_nan(end+1) = i_frame;
            end
            
            if ~isnan(gaze_right_tmp(i_frame,2)) 
                if bounce_down && gaze_right_tmp(i_frame,2) > area_ofinterest_ydown(1)
                        gaze_right_tmp(i_frame,2)
                        heat_map_allscreen(gaze_right_tmp(i_frame,1), gaze_right_tmp(i_frame,2)) = ...
                                heat_map_allscreen(gaze_right_tmp(i_frame,1), gaze_right_tmp(i_frame,2)) +10;
                elseif bounce_down == 0 && gaze_right_tmp(i_frame,2) < area_ofinterest_yup(end) 
                     %gaze_left_tmp(i_frame,2)
                        heat_map_allscreen(gaze_right_tmp(i_frame,1), gaze_right_tmp(i_frame,2)) = ...
                                heat_map_allscreen(gaze_right_tmp(i_frame,1), gaze_right_tmp(i_frame,2)) +1;
                else
                    fixation_right_out(end+1) = i_frame
                end
                
                figure(4); hold on; 
                    subplot(2,8,subplot_order(jj)); scatter(gaze_right_tmp(i_frame,1), gaze_right_tmp(i_frame,2), 'filled')
                %heat_map_allscreen(fixation_x(i_subj,i_frame), fixation_y(i_subj,i_frame))
            else
                fixation_right_nan(end+1) = i_frame;
            end
        end
        
        xline(dot1_x_left); xline(dot1_x_right); 
        xline(dot2_x_left); xline(dot2_x_right); 
        set(gca, 'YDir','reverse')
        
        %title([ 'bounce ' num2str(jj) ' (fixation OUT: ' num2str(fixation_left_out + fixation_right_out) ]);
        title([ title_tmp  num2str(jj) newline ...
            '(fixation OUT: ' num2str(length(fixation_left_out) + length(fixation_right_out)) newline ...
            '(fixation NaN: ' num2str(length(fixation_left_nan) + length(fixation_right_nan)) ]);

    end
       
    %figure; imagesc(heat_map_allscreen); colorbar
    
    
    %%
    
%             %figure; imagesc(heat_map_allscreen); colorbar
%             figure; fig = pcolor(heat_map_allscreen)
%             fig.FaceColor = 'interp'; colorbar
%             find(heat_map_allscreen > 0)
%             % heat map MCS  with fixation y - - - - - 
%             if ~isnan(fixation_y(i_frame))
%                 if bounce_down
%                     %if ismember(fixation_y(i_frame), area_ofinterest_ydown)  
%                         heat_map(fixation_x(i_frame), fixation_y(i_frame)) = ...
%                                 heat_map(fixation_x(i_frame), fixation_y(i_frame)) +1;
%                     end
%                 else
%                      if ismember(fixation_y(i_subj,i_frame), area_ofinterest_yup)  
%                         heat_map(fixation_x(i_subj,i_frame), fixation_y(i_subj,i_frame)) = ...
%                                                     heat_map(fixation_x(i_subj,i_frame), fixation_y(i_subj,i_frame)) +1;
%                      end
%                 end
%                 heat_map(fixation_x(i_subj,i_frame), fixation_y(i_subj,i_frame))
%             end
%         end
%         
%     end 
%     find(heat_map > 0)
%     
% figure; colormap gray;
% imagesc(heat_map'); colorbar
        
        
%         for ww = 1:length(window_vector_frame)
%             %ww = 12
%             i_frame = window_vector_frame(ww)  %for eyetracker
%             
%             if ~isnan(fixation_y(i_subj,i_frame))
%                 if bounce_down
%                     if ismember(fixation_y(i_subj,i_frame), area_ofinterest_ydown)  
%                         heat_map(fixation_x(i_subj,i_frame), fixation_y(i_subj,i_frame)) = ...
%                                                     heat_map(fixation_x(i_subj,i_frame), fixation_y(i_subj,i_frame)) +1;
%                     end
%                 else
%                      if ismember(fixation_y(i_subj,i_frame), area_ofinterest_yup)  
%                         heat_map(fixation_x(i_subj,i_frame), fixation_y(i_subj,i_frame)) = ...
%                                                     heat_map(fixation_x(i_subj,i_frame), fixation_y(i_subj,i_frame)) +1;
%                      end
%                 end
%                 heat_map(fixation_x(i_subj,i_frame), fixation_y(i_subj,i_frame))
%             end
%         end
     
  
