% combokeeper combines the selected cells, collating angle and gfp/cherry kymographs
% creates "combined" which is a cell of {[angle], [GFP kymo], [cherry
% kymo], centroids[highx, highy, lowx, lowy]]}

% Pick "keepers", a variable containing the cell number of the cells you wish to analyze.
keepers = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51];

% pick linethickness for the width of the perimeter averaging
linethickness = 5;

kept=[];
keptfixed = [];
keptmat = [];
areas = [];

% Combined cell contains {cellangles,gfpkymo,cherrykymo}
combined = cell(size(keepers,2),4);
% profiles is an intermediate to hold raw kymographs
profiles = cell(size(keepers,2),2);

%load with excess zeros so data can fit, don't know maximum size of the
%kymograph until its been generated
for i = 1:size(keepers,2);
    profiles{i,1} = zeros(400,tmax);
    profiles{i,2} = zeros(400,tmax);
end

%% Generate Line profile from the masks and GFP and Cherry data
for i = 1:size(keepers,2);
    combined{i,1}= cellsfixed{keepers(i),1};
    combined{i,4} = centroids{keepers(i),1};
    for t = 1:size(TLmask,1);
        gfptemp = [];
        cherrytemp =[];
        %t
        %c
        gfptemp = lineprof(gfpin{t,1},(TLmask{t,1}==keepers(i)),linethickness,0);
        cherrytemp = lineprof(cherryin{t,1},(TLmask{t,1}==keepers(i)), linethickness,0); 
        profiles{i,1}(1:size(gfptemp,1),t) = gfptemp;
        profiles{i,2}(1:size(cherrytemp,1),t) = cherrytemp;
    end
    
end

%% find maximum profile length for each cell and center
maxsize = [];
tempsize = [];
for i = 1:size(profiles,1);
    for t = 1:size(profiles{i,1},2);
        if sum(sum(profiles{i,1}(:,t))) > 1;
        tempsize(i,t) = find(profiles{i,1}(:,t),1,'last');
        else
            tempsize(i,t) = 0;
        end
        
    end
    maxsize(i,1) = max(tempsize(i,:));
    profiles{i,1} = profiles{i,1}(1:maxsize(i,1),:);
    profiles{i,2} = profiles{i,2}(1:maxsize(i,1),:);
end

% Fill empy matrices with correct number of zeros to load centered data
centeredprofiles = cell(size(profiles,1),2);
for i = 1:size(centeredprofiles,1);
    centeredprofiles{i,1}=zeros(maxsize(i,1),size(profiles{i,1},2));
    centeredprofiles{i,2}=zeros(maxsize(i,1),size(profiles{i,1},2));
end

% Center the data

for i = 1:size(profiles,1);
    for t = 1 : size(profiles{i,1},2);
    currentzeros = sum(profiles{i,1}(:,t)==0);
        if currentzeros > 1;
            prespace = floor(currentzeros/2);
            postspace = currentzeros - prespace;
            postindex = size(profiles{i,1},1)-postspace;
            centeredprofiles{i,1}(prespace:postindex,t)= profiles{i,1}(1:postindex-prespace+1,t);
            centeredprofiles{i,2}(prespace:postindex,t)= profiles{i,2}(1:postindex-prespace+1,t);
        else
            centeredprofiles{i,1}(:,t) = profiles{i,1}(:,t);
            centeredprofiles{i,2}(:,t) = profiles{i,2}(:,t);
        end
    end
    
   
end

%% Don't Normalize the Data!
% normprofiles = [];
for i = 1:size(centeredprofiles,1);
%     normprofiles{i,1} = (centeredprofiles{i,1} - min(nonzeros(centeredprofiles{i,1}))) ./ max(max(centeredprofiles{i,1})).*double(centeredprofiles{i,1}>0);
%     normprofiles{i,2} = (centeredprofiles{i,2} - min(nonzeros(centeredprofiles{i,2}))) ./ max(max(centeredprofiles{i,2})).*double(centeredprofiles{i,2}>0);
    combined{i,2} = centeredprofiles{i,1};
    combined{i,3} = centeredprofiles{i,2};
end


%% Example of output

doplot = 0;

if doplot == 1;
    for i = 1:size(combined,1);
        figure();
        subplot(1,3,1), plot(combined{i,1});
        xlim([0,size(combined{i,1},1)]);
        title('Angles')
        subplot(1,3,2), imagesc(combined{i,2});
        title('Bem1')
        subplot(1,3,3), imagesc(combined{i,3});
        title('Cdc3')
    end
end





function [ lineout ] = lineprof(image,mask,width,graph)
%lineprof uses a mask "mask" to do a linescan on the periphery of that mask in
%image "image".  The linescan will go in to the mask by "width" pixels


% create line


if sum(sum(mask)) > 0;
    %make the mask binary
    mask = mask > 0;

    masksmall = [];

    % structuring element for the erosion
    struct = strel('diamond', 1);

    masksmall{1,1} = imerode(mask,struct);

    % find the edge by subtracting the eroded mask from the mask, leaving a 1
    % pixel width line
    objedge= [];
    objedge{1,1} = double(mask)- double(masksmall{1,1});


    % determine the start point for the line, find the minimum x value in the
    % edge, then find the yvalues at that minumum x. 
    [y,x]=find(objedge{1,1});
    minx = min(x);
    [ymin,xmin] = find(objedge{1,1}(:,minx));

    %make boundary ROI (broi), a series of xy values defining lines for
    %improfile to use as the selection to measure, using minimum x and first y
    %associated with that x as a starting point
    broi = [];

    broi{1,1} = bwtraceboundary(objedge{1,1},[ymin(1),minx],'N');

    % define size of the line to use as number of points in the first (largest) line
    linesize = size(x,1);

    c1(:,1) = improfile(image,broi{1,1}(:,2),broi{1,1}(:,1),linesize);

    % figure();
    % plot(c1(:,1));

    %% to do a line of more than one pixel, repeat above steps
    if width > 1;
        for i = 2:width;
           masksmall{i,1}=imerode(masksmall{i-1,1},struct);
           objedge{i,1} = double(masksmall{i-1}) - double(masksmall{i,1});
           [yi,xi]=find(objedge{i,1});
           minxi = min(xi);
           [yimin,ximin]=find(objedge{i,1}(:,minxi));
           try
           broi{i,1}=bwtraceboundary(objedge{i,1},[yimin(1),minxi],'N');
           %makes the profile use the length of the longest line so that the
           %vectors can be averaged
           % use try statement, because if the perimeter collapses from
           % erosion, it breaks the boundary trace, which will still return a
           % boundary, but it won't work in improfile.  In situations where
           % this occurs, the script will just average the number of lines that
           % it can
           catch
               display(['Line Broken at ',num2str(i)]);
               break;
           end
           
           try
            c1(:,i)=improfile(image,broi{i,1}(:,2),broi{i,1}(:,1),linesize);
           catch
               display(['line broken at ',num2str(i)]); % will let the user know if this error is happening
               break;
           end
        end


    else
    end



    % output the average profile
    lineout = mean(c1,2);


    if graph >0;
        figure();
        plot(lineout);
        title('Line Out');
    else
    end
else
    lineout = 0;
end




end

