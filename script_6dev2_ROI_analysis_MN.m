% Tabella per analisi Eye Tracker - Adults (TN)
% to do: fitting retta scatter line 238
% done: figure di ogni soggetto

load('cond_111_fixation.mat')
X_subj = (T.XY(:,1:2:end));
Y_subj = (T.XY(:,2:2:end));
sbj=17;
dots=2;
bounces=8;
screen_center_x=640;

% calculate pre-stimulus - change msec
bounce_msec=50; %50 o 100 msec
ifi=(1/120)*1000;
pre_bounce_frame=bounce_msec/ifi;
post_bounce_frame=bounce_msec/ifi;

% bounce vector of every dot
% max = down, min = up
[idx_max, tmp_max]=max(T.dot1_Yup);
max_dot1=find(T.dot1_Yup==idx_max);
max_dot1(2:2:end)=[];
[idx_min, tmp_min]=min(T.dot1_Yup);
min_dot1=find(T.dot1_Yup==idx_min);
min_dot1(2:2:end)=[];
idx_bounces_dot1=sort([min_dot1; max_dot1])

[idx_max, tmp_max]=max(T.dot2_Yup);
max_dot2=find(T.dot2_Yup==idx_max);
max_dot2(2:2:end)=[];
[idx_min, tmp_min]=min(T.dot2_Yup);
min_dot2=find(T.dot2_Yup==idx_min);
min_dot2(2:2:end)=[];
idx_bounces_dot2=sort([min_dot2; max_dot2]);

% dove guardava il soggetto - asse x (sx-dx)?
pre_bounce_dot1=[]; pre_bounce_dot2=[];
mean_pre_bounce_dot1=[]; mean_pre_bounce_dot2=[];
post_bounce_dot1=[]; post_bounce_dot2=[];
mean_post_bounce_dot1=[]; mean_post_bounce_dot2=[];
table_predot1=[]; table_predot2=[];
table_postdot1=[]; table_postdot2=[];
table_predot1_left=[]; table_predot1_right=[]; table_predot1_center=[];
table_predot2_left=[]; table_predot2_right=[]; table_predot2_center=[];
for i_sbj=1:sbj
    x=X_subj(:,i_sbj);
    
    %calculate pre-stimulus fixation - every bounce
    % dot1: left
    % dot2: right
    idx_prebounce_dot1=idx_bounces_dot1-pre_bounce_frame;
    idx_prebounce_dot2=idx_bounces_dot2-pre_bounce_frame;
    idx_postbounce_dot1=idx_bounces_dot1+1+pre_bounce_frame;
    idx_postbounce_dot2=idx_bounces_dot2+1+pre_bounce_frame;
    
    %calculate pre-stimulus fixation - all bounces
    pre_bounce_dot1=x(idx_prebounce_dot1); %figure;plot(pre_bounce_dot1)
    pre_bounce_dot2=x(idx_prebounce_dot2);
    mean_pre_bounce_dot1(end+1)=nanmean(pre_bounce_dot1);
    mean_pre_bounce_dot2(end+1)=nanmean(pre_bounce_dot2);
    
    for i_bounce=1:bounces
        %pre - mean fixation
        vector_prebounce_idx_dot1=idx_prebounce_dot1(i_bounce):idx_bounces_dot1(i_bounce)
        mean_prebounce_fixation_dot1=nanmean(x(vector_prebounce_idx_dot1));
        table_predot1(i_sbj, i_bounce)=mean_prebounce_fixation_dot1;
        
        vector_prebounce_idx_dot2=idx_prebounce_dot2(i_bounce):idx_bounces_dot2(i_bounce);
        mean_prebounce_fixation_dot2=nanmean(x(vector_prebounce_idx_dot2));
        table_predot2(i_sbj, i_bounce)=mean_prebounce_fixation_dot2;

        %pre - laterality    
        %dot1
        left_frames=0; center_frames=0; right_frames=0; nan_frames=0;
        for i_frame=1:length(vector_prebounce_idx_dot1)
            if ~isnan(x(vector_prebounce_idx_dot1(i_frame)))
                if x(vector_prebounce_idx_dot1(i_frame))<screen_center_x
                    left_frames=left_frames+1;
                elseif x(vector_prebounce_idx_dot1(i_frame))==screen_center_x
                    center_frames=center_frames+1;
                elseif x(vector_prebounce_idx_dot1(i_frame))>screen_center_x
                    right_frames=right_frames+1;      
                end
            elseif isnan(x(vector_prebounce_idx_dot1(i_frame)))
                nan_frames=nan_frames+1;
            end
        end
        if nan_frames>0 && left_frames>0
            table_predot1_left(i_sbj, i_bounce)=left_frames;
            table_predot1_right(i_sbj, i_bounce)=NaN;
            table_predot1_center(i_sbj, i_bounce)=NaN;
        elseif nan_frames>0 && right_frames>0
            table_predot1_left(i_sbj, i_bounce)=NaN;
            table_predot1_right(i_sbj, i_bounce)=right_frames;
            table_predot1_center(i_sbj, i_bounce)=NaN;
        elseif nan_frames==0
            table_predot1_left(i_sbj, i_bounce)=left_frames;
            table_predot1_right(i_sbj, i_bounce)=right_frames;
            table_predot1_center(i_sbj, i_bounce)=center_frames;
        end
        
        %dot 2
        left_frames=0; center_frames=0; right_frames=0; nan_frames=0;
        for i_frame=1:length(vector_prebounce_idx_dot2)
            if ~isnan(x(vector_prebounce_idx_dot2(i_frame)))
                if x(vector_prebounce_idx_dot2(i_frame))<screen_center_x
                    left_frames=left_frames+1;
                elseif x(vector_prebounce_idx_dot2(i_frame))==screen_center_x
                    center_frames=center_frames+1;
                elseif x(vector_prebounce_idx_dot2(i_frame))>screen_center_x
                    right_frames=right_frames+1;      
                end
            elseif isnan(x(vector_prebounce_idx_dot2(i_frame)))
                nan_frames=nan_frames+1;
            end
        end
        if nan_frames>0 && left_frames>0
            table_predot2_left(i_sbj, i_bounce)=left_frames;
            table_predot2_right(i_sbj, i_bounce)=NaN;
            table_predot2_center(i_sbj, i_bounce)=NaN;
        elseif nan_frames>0 && right_frames>0
            table_predot2_left(i_sbj, i_bounce)=NaN;
            table_predot2_right(i_sbj, i_bounce)=right_frames;
            table_predot2_center(i_sbj, i_bounce)=NaN;
        elseif nan_frames==0
            table_predot2_left(i_sbj, i_bounce)=left_frames;
            table_predot2_right(i_sbj, i_bounce)=right_frames;
            table_predot2_center(i_sbj, i_bounce)=center_frames;
        end
        
        %post - mean fixation
        vector_postbounce_idx_dot1=idx_bounces_dot1(i_bounce)+1:idx_postbounce_dot1(i_bounce);
        mean_postbounce_fixation_dot1=nanmean(x(vector_postbounce_idx_dot1));
        table_postdot1(i_sbj, i_bounce)=mean_postbounce_fixation_dot1;
        
        vector_postbounce_idx_dot2=idx_bounces_dot2(i_bounce)+1:idx_postbounce_dot2(i_bounce);
        mean_postbounce_fixation_dot2=nanmean(x(vector_postbounce_idx_dot2));
        table_postdot2(i_sbj, i_bounce)=mean_postbounce_fixation_dot2;

        %post - laterality        
        %dot1
        left_frames=0; center_frames=0; right_frames=0; nan_frames=0;
        for i_frame=1:length(vector_postbounce_idx_dot1)
            if ~isnan(x(vector_postbounce_idx_dot1(i_frame)))
                if x(vector_postbounce_idx_dot1(i_frame))<screen_center_x
                    left_frames=left_frames+1;
                elseif x(vector_postbounce_idx_dot1(i_frame))==screen_center_x
                    center_frames=center_frames+1;
                elseif x(vector_postbounce_idx_dot1(i_frame))>screen_center_x
                    right_frames=right_frames+1;      
                end
            elseif isnan(x(vector_postbounce_idx_dot1(i_frame)))
                nan_frames=nan_frames+1;
            end
        end
        if nan_frames>0 && left_frames>0
            table_postdot1_left(i_sbj, i_bounce)=left_frames;
            table_postdot1_right(i_sbj, i_bounce)=NaN;
            table_postdot1_center(i_sbj, i_bounce)=NaN;
        elseif nan_frames>0 && right_frames>0
            table_postdot1_left(i_sbj, i_bounce)=NaN;
            table_postdot1_right(i_sbj, i_bounce)=right_frames;
            table_postdot1_center(i_sbj, i_bounce)=NaN;
        elseif nan_frames==0
            table_postdot1_left(i_sbj, i_bounce)=left_frames;
            table_postdot1_right(i_sbj, i_bounce)=right_frames;
            table_postdot1_center(i_sbj, i_bounce)=center_frames;
        end
        
        %dot 2
        left_frames=0; center_frames=0; right_frames=0; nan_frames=0;
        for i_frame=1:length(vector_postbounce_idx_dot2)
            if ~isnan(x(vector_postbounce_idx_dot2(i_frame)))
                if x(vector_postbounce_idx_dot2(i_frame))<screen_center_x
                    left_frames=left_frames+1;
                elseif x(vector_postbounce_idx_dot2(i_frame))==screen_center_x
                    center_frames=center_frames+1;
                elseif x(vector_postbounce_idx_dot2(i_frame))>screen_center_x
                    right_frames=right_frames+1;      
                end
            elseif isnan(x(vector_postbounce_idx_dot2(i_frame)))
                nan_frames=nan_frames+1;
            end
        end
        if nan_frames>0 && left_frames>0
            table_postdot2_left(i_sbj, i_bounce)=left_frames;
            table_postdot2_right(i_sbj, i_bounce)=NaN;
            table_postdot2_center(i_sbj, i_bounce)=NaN;
        elseif nan_frames>0 && right_frames>0
            table_postdot2_left(i_sbj, i_bounce)=NaN;
            table_postdot2_right(i_sbj, i_bounce)=right_frames;
            table_postdot2_center(i_sbj, i_bounce)=NaN;
        elseif nan_frames==0
            table_postdot2_left(i_sbj, i_bounce)=left_frames;
            table_postdot2_right(i_sbj, i_bounce)=right_frames;
            table_postdot2_center(i_sbj, i_bounce)=center_frames;
        end
    end
   
%figure bounces dot 1 - single subject
% X_prebounces=nanmean(table_predot1(i_sbj,:),1);
% figure;scatter(1:8, X_prebounces); 
% figure_name=strcat('Sbj', num2str(i_sbj), '_predot1_fixation_bounces.jpg')
% saveas(gcf,figure_name)
% 
% X_postbounces=nanmean(table_postdot1(i_sbj,:),1);
% figure;scatter(1:8, X_postbounces); 
% figure_name=strcat('Sbj', num2str(i_sbj), '_postdot1_fixation_bounces.jpg')
% saveas(gcf,figure_name)
end

%dove sta guardando prima del rimbalzo di ogni pallina?
%table all, rows=sbj, columns=2dot(pre-dot1, pre-dot2)*8bounces
table_pre_dot12=[table_predot1 table_predot2];
table_post_dot12=[table_postdot1 table_postdot2];

X_down=table_pre_dot12(:,1);
Y_down=table_pre_dot12(:,7);
[H,P,CI,STATS]=ttest(X_down,Y_down)

X_up=table_pre_dot12(:,2);
Y_up=table_pre_dot12(:,8);
[H2,P2,CI2,STATS2]=ttest(X_up,Y_up)

X_mean=nanmean([table_pre_dot12(:,1),table_pre_dot12(:,2)], 2);
Y_mean=nanmean([table_pre_dot12(:,7),table_pre_dot12(:,8)], 2);
[H3,P3,CI3,STATS3]=ttest(X_mean,Y_mean)

%scatter
X_mean_down=nanmean([table_pre_dot12(:,1:2:8)], 1);
X_mean_up=nanmean([table_pre_dot12(:,2:2:8)], 1);

% figure;scatter(1:4, X_mean_down); hold on;
% scatter(1:4, X_mean_up); legend('mean down-bounces', 'mean up-bounces');
% title('Mean fixation across bounces')

% fitting della retta by AV
X_prebounces=nanmean(table_pre_dot12(:,1:8),1);
figure;scatter(1:8, X_prebounces); 
% hold on;
% scatter(1:4, X_mean_up);
X_postbounces=nanmean(table_post_dot12(:,1:8),1);
figure;scatter(1:8, X_postbounces); 

X_prebounces=nanmean(table_pre_dot12(:,9:16),1);
figure;scatter(1:8, X_prebounces); 
% hold on;
% scatter(1:4, X_mean_up);
X_postbounces=nanmean(table_post_dot12(:,9:16),1);
figure;scatter(1:8, X_postbounces); 


%lateralità dx-sx rispetto alle palline
%table all, rows=sbj, columns=2dot*2side(left-right)*8bounces
table_pre_dot12_side=[table_predot1_left table_predot1_right table_predot2_left table_predot2_right];
table_post_dot12_side=[table_postdot1_left table_postdot1_right table_postdot2_left table_postdot2_right];

X=table_pre_dot12_side(:,1);
Y=table_pre_dot12_side(:,9);
[H,P,CI,STATS]=ttest(X,Y)

X_down=table_pre_dot12_side(:,1+7);
Y_down=table_pre_dot12_side(:,9+7);
[H,P,CI,STATS]=ttest(X_down,Y_down)

X=table_pre_dot12_side(:,17+7);
Y=table_pre_dot12_side(:,25+7);
[H,P,CI,STATS]=ttest(X,Y)