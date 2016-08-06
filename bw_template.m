clear all
run('~/Documents/MATLAB/CompVis/project/constants.m')

%%% CODE TO FIND BINARY TEMPLATE INDIVIDUAL PIECES

global TEMPLATE1            TEMPLATE1_MOMENTS...
       TEMPLATE1_GRADIENT   TEMPLATE1_GRAD_MOMENTS...
       TEMPLATE2            TEMPLATE2_GRADIENT...
       TEMPLATE_RADII       TEMPLATE1_GRAD_MOMENTS_BLURRED...
       TEMPLATE_CENTERS

TEMPLATE1_BLURRED = imgaussfilt(TEMPLATE1,2.5);
TEMPLATE_BW = im2bw(TEMPLATE1_BLURRED,graythresh(TEMPLATE1));
figure(100), imshow(TEMPLATE_BW,[])

[LW,nw] = bwlabel(TEMPLATE_BW);
whiteStats = regionprops(TEMPLATE_BW,...
    'Area','Centroid','BoundingBox','Perimeter','FilledArea','Solidity');
[Lb,nb] = bwlabel(~TEMPLATE_BW);
blackStats = regionprops(~TEMPLATE_BW,...
    'Area','Centroid','BoundingBox','Perimeter','FilledArea','Solidity');

% for k = 1:nb
%     rectangle('Position', blackStats(k).BoundingBox, ...
%                     'EdgeColor', 'r', 'Linewidth', 2);
%     pause
% end

% 2, 5, 10, 13, 15, 17   ---   White blobs
% 2, 3, 5,  6,  9,  11   ---   Black blobs
windex = [2 5 10 13 15 17];
whiteStats = whiteStats(windex);
whiteRatios = zeros(1,6);
for i = 1:6
    whiteRatios(i) = (whiteStats(i).Solidity + 2) ...
        * whiteStats(i).FilledArea / whiteStats(i).Perimeter;
end


% [iy,ix] = size(TEMPLATE_BW);
% polarity = [];
% for iy = 1:size(TEMPLATE_BW,1)
%     for ix = 1:size(TEMPLATE_BW,2)
%         if (ix - TEMPLATE_CENTERS(i,1))^2 ...
%                 + (iy - TEMPLATE_CENTERS(i,2))^2 <= TEMPLATE_RADII(1)^2
%             polarity = [polarity; TEMPLATE_BW(iy,ix)];
%         end
%     end
% end
% 
% polarity = mode(mode(polarity));

% Finds blobs that are of the same color as the piece icon
[Lw,nw] = bwlabel(TEMPLATE_BW);
[Lb,nb] = bwlabel(~TEMPLATE_BW);
whiteStats = regionprops(Lw,'Area','Centroid', ...
    'BoundingBox','Perimeter','FilledArea','Solidity');
blackStats = regionprops(Lb,'Area','Centroid', ...
    'BoundingBox','Perimeter','FilledArea','Solidity');
blobStats = [whiteStats; blackStats];
nBlobs = nw + nb;

wrongBlobs = [];
% Removes blobs that are a large portion of the image
for l = 1:nBlobs   
    bbArea = blobStats(l).BoundingBox(3)*blobStats(l).BoundingBox(4);
    if bbArea > 0.9*numel(TEMPLATE_BW)            % CHANGE THIS TOLERANCE
        wrongBlobs = [wrongBlobs, l];
    end
end    

% Determines which of those blobs is nearest to the center
% of the circle surrounding the piece (or, if no circle was
% found, to the center of the image)
diffs = zeros(nBlobs,1);
for l = 1:nBlobs
    diffs(l) = norm(blobStats(l).Centroid - TEMPLATE_CENTERS(i));
end                 
diffs(wrongBlobs) = NaN;

mindex = find(diffs == min(diffs));


figure(99), imshow(TEMPLATE1,[])
binaryMoments = zeros(12,7);
correspCenters = zeros(12,2);
q = 1;
for p = 1:nBlobs
    if any(p == [2, 5, 10, 13, 15, 17, 22, 23, 25, 26, 29, 31])
    
    bb = blobStats(p).BoundingBox;
    
    a = 0;
    if q > 6
        a = 1;
    end
    
    miniTemplate = a + zeros(125,125);
    miniTemplate( 62-bb(4)/2:62+bb(4)/2-1, 62-bb(3)/2:62+bb(3)/2-1) = TEMPLATE_BW( bb(2):bb(2)+bb(4)-1, bb(1):bb(1)+bb(3)-1 );
    figure(45), imshow(miniTemplate,[])
    eta = SI_Moment(miniTemplate,ones(size(miniTemplate)));
    inv_moments = Hu_Moments(eta);
    binaryMoments(q,:) = real(log(inv_moments));
    correspCenters(q,:) = blobStats(p).Centroid;
    
    q = q+1;
    end
end



    
 



