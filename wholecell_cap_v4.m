% wholecell_cap description, Joshua B Kelley, Elston Lab, UNC Chapel Hill
%
% Inputs for wholecell_cap are 3 tifs, mask, gfp, cherry
% •	maskdir, directory location of *mask.tif
% •	gfpdir, directory location of *gfp.tif
% •	cherrydir, directory location of *cherry.tif
% 
% 
% The threshold percentile 'threshper' and the lower threshhold percentile,
% 'threshperl' are used to calculate the polar cap angle of orientation.
% If the angle of orientation does not seem to be working, attempt
% different threshold values.
%
% Useful Outputs from wholecell_cap
% •	TLmask, a (time x 1) cell of sorted and labeled cell masks [x,y]
% •	gfpin, gfp images
% •	cherryin, cherry images
% •	gfpcell, a (cell# x 1) cell of matrices [x,y,t] which is just masked gfp data for each image
% •	celldata, a cell of (cell# x 1) containing matrices of [cellx,celly,peakx,peaky,angle,isreal] columns 3,4,5 are not currently used
% •	centroids, a cell of (cell# x 1) containing matrices of [highx,highy,lowx,lowy] which are the x and y values for the centroids of the high and low thresholded bem1 to determine cap angle
% •	cellangles, a cell (cell# x 1) of matrices (time x 1) containing angle data in degrees

%% Pick images to be input for analysis

maskdir = 'Pos1_Mask of MAX_Merged.tif'; %'Result of Erode3_TRACKING_Mask_Pos1_8.7.18.tif'
gfpdir = 'AVG_4.8.19_sst2_s539a_fusgfp_bem1ruby_300nM_12hr_20min_pos1_rfp.tif'; % bem1
cherrydir = 'GFP_Sub.tif'; % fus3 to be aligned, has had nucleus subtracted

preimin = cell(145,1);

tmin = 1;
tmax = [];
tmax = size(imfinfo(maskdir),1);
times = [tmin:1:tmax];

threshper = 99;
threshperl = 88;

maskin01 = [];
trackvar = [];

%% load images in
for c = 1:3;
    if c == 1;
        gadir = maskdir;
    elseif c==2;
        gadir = gfpdir;
    elseif c == 3;
        gadir = cherrydir;
    end
    
        
    for i = 1:tmax;
        preimin{i,1} = imread(gadir,i);
    end


    for i = 1:tmax
        if c ==1;
            maskin{i,1} = preimin{i,1};
        elseif c == 2;
            gfpin{i,1}= preimin{i,1};
        elseif c == 3;
            cherryin{i,1}=preimin{i,1};
        end
    end
end

    


 figure();
 subplot(1,3,1), imagesc(maskin{tmax,1});
 subplot(1,2,1), imagesc(maskin{tmax,1});
 subplot(1,3,3), imagesc(cherryin{tmax,1});

tnum = tmax - tmin + 1;
max_field = zeros(tnum,1);
peak_thresh = cell(tmax,1);
currmax = [];
low_thresh_clean = cell(tnum,1);
low_thresh_count = zeros(tnum,1);
peak_thresh_count = zeros(tnum,1);
labeledmask_peak = cell(tnum,1);
labeledmask = cell(tnum,1);


% find the maximum value present in each time frame
% for i = 1:tnum;
%     max_field(i,1) = max(max(imin{i,1}));
% end

% convert the mask to a logical
for i = 1:tnum;
    maskin01{i,1} = ~maskin{i,1} > 0; %mask out of imagej is inverted
    cellmask{i,1} = uint16(maskin01{i,1});
end



figure();
subplot(1,2,1), 
imagesc(gfpin{tmax,1}.*cellmask{tmax,1});
subplot(1,2,2), imagesc(cherryin{tmax,1}.*cellmask{tmax,1});


%% label the mask and keep track of number of cells in each time point
for i = 1:tmax;
    labeledmask{i,1} = bwlabel(maskin01{i,1});
    cellcount(i,1) = max(max(labeledmask{i,1}));
end



%% use jtrack to sort the cells.  Use x,y,time,label
total = 1;
for i = 1:tnum;
    for j = 1:cellcount(i,1);
        currmask = [];
        currmask = labeledmask{i,1}==j;
        [ys,xs] = find(currmask); 
        trackvar(total,1) = mean(xs);
        trackvar(total,2) = mean(ys);
        trackvar(total,3) = i;
        trackvar(total,4) = j;
        total = total+1; %keeps track through the larger loop
    end
end

tracked = jtrackv3(trackvar,0,100); %had to up the max distance to 30 pixels for the tracking to work

%% relabel the mask based on cell identity. T(rack)L(abel)mask
% start with tracked cell 1, and assign all of its indices to 1 in the new
% mask.
TLmask = cell(tmax,1);
for i = 1:tmax;
    TLmask{i,1} = zeros(size(labeledmask{i,1},1),size(labeledmask{i,1},2));
end

for i = 1:size(tracked,1);
    for t = 1:tnum;
        if tracked{i,1}(t,4)==1;
        currentID = tracked{i,1}(t,3);
        currmask = labeledmask{t,1} == currentID;
        TLmask{t,1}(currmask) = i;
        else
        end
    end
end

clear labeledmask

%% Output a TIFF of the tracked mask for trouble shooting/validation
imwrite(uint8(TLmask{1,1}),'tracked_mask.tif', 'WriteMode', 'OverWrite');
for i = 2:tmax;
    imwrite(uint8(TLmask{i,1}),'tracked_mask.tif','WriteMode', 'append');
end

%% The following sections need cleaning up

% need to find the high and low thresholds in gfp and store the centroids
% of them in cell data


hicents = [];
lowcents = [];
for i = 1:size(tracked, 1);
    for t = 1:tnum;
        gfpcell = [];
        gfpcell= gfpin{t,1}.*uint16(TLmask{t,1}==i);
        if size(nonzeros(gfpcell),1) > 2;
            temprop = [];
            temprop =regionprops(gfpcell > jthresh(gfpcell, threshper));
            hicents{t,i} = temprop.Centroid;
            temprop = [];
            temprop = regionprops(gfpcell > jthresh(gfpcell, threshperl));
            lowcents{t,i} = temprop.Centroid;
%             gfpthresh_hi{i,1}(:,:,t) = gfpcell{i,1}(:,:,t) > jthresh(gfpcell{i,1}(:,:,t), threshper); 
%             gfpthresh_low{i,1}(:,:,t) = gfpcell{i,1}(:,:,t) > jthresh(gfpcell{i,1}(:,:,t), threshperl);
            per_threshs{i,1}(t,1) = jthresh(gfpcell, threshper); %store the values of the jthresh output
            per_threshs{i,1}(t,2) = jthresh(gfpcell, threshperl);
        else
            hicents{t,i} = [0,0];
            lowcents{t,i} = [0,0];
%             gfpthresh_hi{i,1}(:,:,t) = zeros(size(gfpcell{i,1}(:,:,t),1),size(gfpcell{i,1}(:,:,t),2)); 
%             gfpthresh_low{i,1}(:,:,t) = zeros(size(gfpcell{i,1}(:,:,t),1),size(gfpcell{i,1}(:,:,t),2));
            per_threshs{i,1}(t,1) = 0; %store the values of the jthresh output
            per_threshs{i,1}(t,2) = 0;   
        end
        
    end
end

%% Threshold the GFP within each cell

% separate each cell's gfp image out
% gfpcell = [];
% 
% for i = 1:size(tracked,1);
%     for t = 1:tnum;
%         gfpcell{i,1}(:,:,t)= gfpin{t,1}.*uint16(TLmask{t,1}==i);
%     end
% end
% %% Calculate thresholded images for a high value and low value(per cell per timepoint)and
% %  a single value which will work for all timepoints
% 
% gfpthresh_hi = [];
% gfpthresh_low = [];
% per_threshs =[];
% for i = 1:size(tracked, 1);
%     for t = 1:tnum;
%         if size(nonzeros(gfpcell{i,1}(:,:,t)),1) > 2;
%             gfpthresh_hi{i,1}(:,:,t) = gfpcell{i,1}(:,:,t) > jthresh(gfpcell{i,1}(:,:,t), threshper); 
%             gfpthresh_low{i,1}(:,:,t) = gfpcell{i,1}(:,:,t) > jthresh(gfpcell{i,1}(:,:,t), threshperl);
%             per_threshs{i,1}(t,1) = jthresh(gfpcell{i,1}(:,:,t), threshper); %store the values of the jthresh output
%             per_threshs{i,1}(t,2) = jthresh(gfpcell{i,1}(:,:,t), threshperl);
%         else
%             gfpthresh_hi{i,1}(:,:,t) = zeros(size(gfpcell{i,1}(:,:,t),1),size(gfpcell{i,1}(:,:,t),2)); 
%             gfpthresh_low{i,1}(:,:,t) = zeros(size(gfpcell{i,1}(:,:,t),1),size(gfpcell{i,1}(:,:,t),2));
%             per_threshs{i,1}(t,1) = 0; %store the values of the jthresh output
%             per_threshs{i,1}(t,2) = 0;   
%         end
%         
%     end
% end


% gfpconstantthresh = [];
% 
% 
% %% Find Area of thresholded Bem1
% 
% 
% for i = 1:size(tracked,1);
%     for t=1:tnum;
%         if sum(sum(gfpcell{i,1}(:,:,t)))>0;
%         gfpconstantthresh{i,1}(:,:,t) = gfpcell{i,1}(:,:,t) > min(min(per_threshs{i,1}));
%         else
%         end
%     end
% end

        


%% find the centroid of each cell
celldata = cell(size(tracked,1),1); % cellx, celly, peakx, peaky, angle, isreal

for i = 1:size(tracked,1);
    for t=1:tnum;
        % find the centroid of the cell mask        
        if sum(find(TLmask{t,1}==i))>0;
            [ys,xs] = find(TLmask{t,1} == i);
            celldata{i,1}(t,1) = mean(xs);
            celldata{i,1}(t,2) = mean(ys);
            celldata{i,1}(t,6) = 1;
        else
           celldata{i,1}(t,1) = 0;
           celldata{i,1}(t,2) = 0; 
           celldata{i,1}(t,6) = 0;
        end
    end
end

figure();
imagesc(TLmask{end,1});
hold on
for i = 1:size(tracked,1);
text(celldata{i,1}(end,1),celldata{i,1}(end,2), ['Cell',num2str(i)], 'FontSize', 14);
end


%% calculate the angle of the polar cap in each cell based on the high and 
% low thresholds

centroids = cell(size(tracked,1),1); % centroids cell of [highx, highy, lowx, lowy];

for i = 1:size(tracked, 1);
    for t = 1:tnum;
        centroids{i,1}(t,1:2) = hicents{t,i};
        centroids{i,1}(t, 3:4) = lowcents{t,i};
    end
end
% 
% 
% for i = 1:size(tracked,1);
%     for t = 1:tnum;
%         if sum(sum(gfpthresh_hi{i,1}(:,:,t))) > 0;
%         [ys,xs] = find(gfpthresh_hi{i,1}(:,:,t));
%         centroids{i,1}(t,1)=mean(xs);
%         centroids{i,1}(t,2)=mean(ys);
%         else
%         centroids{i,1}(t,1) = 0;
%         centroids{i,1}(t,2) = 0;
%         end
%         if sum(sum(gfpthresh_low{i,1}(:,:,t))) > 0;
%         [ys,xs] = find(gfpthresh_low{i,1}(:,:,t));
%         centroids{i,1}(t,3) = mean(xs);
%         centroids{i,1}(t,4) = mean(ys);
%         else
%         centroids{i,1}(t,3) = 0;
%         centroids{i,1}(t,4) = 0;    
%         end
%         
%     end
% end

%% calculate the angles from the xy positions

cellangles = cell(size(tracked,1),1);

for i = 1:size(tracked,1);
    for t = 1:tnum;
        if sum(sum(centroids{i,1}))>0;
        cellangles{i,1}(t,1) = atan2d(centroids{i,1}(t,2)-centroids{i,1}(t,4),centroids{i,1}(t,1)-centroids{i,1}(t,3));
        else
            cellangles{i,1}(t,1) = NaN;
        end
    end
end

        

% This script will take the angle output from wholecell_cap and adjust the
% angles when the angle crosses the 180/-180 degree line.  If an angle is
% positive and it crosses 180, it will be converted to an angle larger than
% 180, rather than an angle between 0 and -180.  The opposite is true for
% an angle that starts out negative.

angle_diff = cell(size(tracked,1),1);
cellsfixed = cellangles;

for i= 2:size(cellangles,1);
    for t = 2:size(cellangles{i,1},1);
        if abs(cellsfixed{i,1}(t,1)-cellsfixed{i,1}(t-1,1)) > 180;
           anglenew = (180-abs(cellangles{i,1}(t,1))) + (180-abs(cellangles{i,1}(t-1,1))); % calculate the absolute difference of both angles from 180 
           
           if cellangles{i,1}(t-1,1) >= 0;
               cellsfixed{i,1}(t,1) = cellsfixed{i,1}(t-1,1) + anglenew; % if the angle is positive, keep it going positive
           else
               cellsfixed{i,1}(t,1) = cellsfixed{i,1}(t-1,1) - anglenew; % if the angle is negative, keep it going negative
           end
       else
        end
    end
end



function [ output_args ] = jthresh( image, percentile )
%jthresh takes an image and a precentile value, and returns the value that
%falls at that percentile.  This value can then be used to threhold the
%image to leave the pixels above the indicated percentile.
    immax = max(max(nonzeros(image)));
    immin = min(min(nonzeros(image)));
           
    bins = [immin:(immax-immin)/98:immax];
    imagehist = histc(nonzeros(image),bins);
    if size(imagehist,1) > 0;
        target = sum(imagehist)*(percentile/100);
        currsum = 0;
        for i = 1:100;
            currsum = currsum + imagehist(i,1);
            if currsum >= target;
                targetbin = i;
                break;
            else
                targetbin = 100;
            end
        end
        output_args = bins(1,targetbin);
    else
        output_args = 0;
    

end

end