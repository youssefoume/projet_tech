clear;clc;
close all;
base_path='U:\base_donnes';
% videos={'ball';'basketball';'board';'book';'bus';'bus2';'campus';'car';'car2';'car3';'card';'coin';'coke';'drive';'excavator';'face';'face2';'forest';'forest2';'fruit';'hand';'kangaroo';'paper';'pedestrain';'pedestrian2';'player';'playground';'rider1';'rider2';'rubik';'student';'toy1';'toy2';'trucker';'worker'};OPs = zeros(numel(videos),1);
videos={'bus'};
FPSs = zeros(numel(videos),1);
videoNum=numel(videos);
distance_rec=zeros(videoNum,50);
PASCAL_rec=zeros(videoNum,50);
average_cle_rec=zeros(videoNum,50);
learning_rate =0.0023;
for vid = 1:videoNum
    close all;
    video_path = [base_path '/' videos{vid}];
    [seq, ground_truth] = load_video_info(video_path);
    seq.VidName = videos{vid};
    st_frame = 1;
    en_frame = seq.len;
    seq.st_frame = st_frame;
    seq.en_frame = en_frame;
    gt_boxes = [ground_truth(:,1:2), ground_truth(:,1:2) + ground_truth(:,3:4) - ones(size(ground_truth,1), 2)];
    
    % Run MHT- main function
    
    results = run_MHT(seq, video_path, learning_rate);
    results.gt = gt_boxes;
    pd_boxes = results.res;
    %perf mesures
    [distance_rec(vid,:),PASCAL_rec(vid,:),average_cle_rec(vid,:),~,~]= computeMetric2(pd_boxes,ground_truth);
    FPSs = results.fps;
    result=matfile(sprintf(strcat(videos{vid},'trackingMHT.mat')),'writable',true);
    result.results=results;
end
save(strcat('trackingMHT','.mat'),'FPSs','PASCAL_rec','average_cle_rec','distance_rec');




