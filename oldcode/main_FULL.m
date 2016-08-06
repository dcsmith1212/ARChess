%%% Fourier descriptors

clear variables
close all
imtool close all

warning off MATLAB:colon:nonIntegerIndex
warning off images:initSize:adjustingMag

run('~/Documents/MATLAB/CompVis/project/constants.m')

global TEMPLATE1            TEMPLATE1_MOMENTS...
       TEMPLATE1_GRADIENT   TEMPLATE1_GRAD_MOMENTS...
       TEMPLATE2            TEMPLATE2_GRADIENT...
       TEMPLATE_RADII       TEMPLATE_CENTERS...
       TEMPLATE1_GRAD_MOMENTS_BLURRED...
       TEMPLATE1_BINARY_MOMENTS_WHITE    CORRESP_CENTERS...
       TEMPLATE_BW          TEMPLATE_RATIOS...
       TEMPLATE1_BINARY_MOMENTS_BLACK
   
% Open movie file.
% movieObj = VideoReader('board.mp4');
% nFrames = movieObj.NumberOfFrames;
% fprintf('Opening movie file with %d images\n', nFrames);
% Go through movie. We don't need to process every frame.
nFrames = 1;

% The size of orthophoto of the board
N = 1000;

all_moments = [];

squareCenters = zeros(64,2);

pieceNames = {'white pawn','white bishop','white knight','white rook','white queen','white king';
              'black pawn','black bishop','black knight','black rook','black queen','black king'};
                          
% Will store locations of pieces amd their ID at each frame
boardLayout = cell(8,8,2);

for iFrame=1:10:nFrames
    I = imread('test5.jpg');        % Or open single frame
    
    
%     I = read(movieObj,iFrame);
%     fprintf('Frame %d\n', iFrame);
    % Reduce image size; is faster and we don't need full size to find
    % board.
    if size(I,2)>1280
        I = imresize(I, 1280/size(I,2));
    end

    % Find the checkerboard. Return the four outer corners as a 4x2 array,
    % in the form [ [x1,y1]; [x2,y2]; ... ].
    [corners, nMatches, avgErr, intersections] = findCheckerBoard(I);
    
    figure(87), imshow(I);
    line([corners(1,1) corners(2,1)], [corners(1,2) corners(2,2)], 'Color', 'g', 'Linewidth', 2);
    line([corners(1,1) corners(4,1)], [corners(1,2) corners(4,2)], 'Color', 'g', 'Linewidth', 2);
    line([corners(3,1) corners(2,1)], [corners(3,2) corners(2,2)], 'Color', 'g', 'Linewidth', 2);
    line([corners(3,1) corners(4,1)], [corners(3,2) corners(4,2)], 'Color', 'g', 'Linewidth', 2);
    title('Identified chess board with square intersections')

    %%% Displays red dots at centers of each square
    hold on
    k = 1;
    % plot(intersections(:,1), intersections(:,2), 'co');
%     for i = 1:length(intersections)-9
%         if mod(i,9) ~= 0
%             avg = (intersections(i,:) + intersections(i+1,:) + intersections(i+9,:) + intersections(i+10,:))/4;
%             plot(avg(1), avg(2), '.r', 'MarkerSize', 8);
%             squareCenters(k,:) = avg;
%             k = k + 1;
%         end
%     end
    plot(intersections(:,1),intersections(:,2),'c.','MarkerSize',20)
    hold off
    %%%
    
    orthoPts = [
        1, 1;
        N, 1;
        N, N;
        1, N
        ];
    
blah = 1;
    % Compute transform, from corresponding control points
    Tform1 = fitgeotrans(corners,orthoPts,'projective');
    % Transform input image to output image
    ref2Doutput = imref2d( ...
        [N, N], [0,N], [0,N]);
    Iout1 = imwarp(I,Tform1,'OutputView',ref2Doutput);

    % Transforms the intersection points found from findCheckBoard to the
    % orthophoto's coordinates
    [ortho_intersect_x, ortho_intersect_y] = transformPointsForward(...
        Tform1, intersections(:,1), intersections(:,2));
    [ortho_centers_x, ortho_centers_y] = transformPointsForward(Tform1, squareCenters(:,1), squareCenters(:,2));
    
    % Display the orthophoto of the board, with intersections and square
    % centers
    figure(87), imshow(Iout1,[]);
    hold on 
    plot(ortho_intersect_x,ortho_intersect_y,'.b','MarkerSize',20);
    %plot(sC2x,sC2y,'.r','MarkerSize',8);
    hold off

    for i = 35:length(ortho_intersect_x)-9
        if mod(i,9) ~= 0
            x0 = ortho_intersect_x(i);  xf = ortho_intersect_x(i+10);
            y0 = ortho_intersect_y(i);  yf = ortho_intersect_y(i+10);
            %rectangle('Position', [x0 y0 xf-x0 yf-y0], 'EdgeColor', 'g', 'LineWidth', 2);          
            
            square = Iout1(y0:yf,x0:xf);
            sub = square(end/4:3*end/4,end/4:3*end/4);
            %square = square(end/4:3*end/4,end/4:3*end/4);
            

            if max(max(sub)) - min(min(sub)) > 100  % For those squares that have a piece on them

%                 % ---------------------------------------------------------
%                 % Method 0: Use SIFT to match test square to template. This
%                 % method doesn't work very, except for kings and queens
%                 % which have lots of features. It fails every time for
%                 % matching pawns, as they are relatively featureless
%                 % ---------------------------------------------------------
%                 correctSpots1 = SiftComparer(square,TEMPLATE1);             
%                 determinePieceName(TEMPLATE1,correctSpots1,TEMPLATE1,correctSpots1);


                % Hoff's Code
%                 E = edge(square, 'Canny');
%                 sensitivity = 0.98;  % default is 0.85
%                 [centers_d,radii_d, metrics_d] = imfindcircles(E,[40 48], ...
%                 'Sensitivity', sensitivity, ...
%                 'ObjectPolarity','dark');
% 
%                 [~,ind] = max(metrics_d);   % Pick the strongest circle
%                 figure(21), imshow(E, []);
%                 if ~isempty(ind) 
%                     % Draw on top of the subimage.
%                     viscircles(centers_d(ind,:),radii_d(ind),'Color','r');
%                 end
% 
%  
%                 
%                 if ~isempty(ind)
%                 pieceCenter = centers_d(ind,:);
%                 pieceRadius = radii_d(ind);
%                 % Creates a mask of the test square, with ones where the
%                 % piece is (using the circle just found)
%                 cx = pieceCenter(1);       cy = pieceCenter(2);
%                 ix = size(square,1);       iy = size(square,2);
%                 r = pieceRadius-5;
%                 [x,y] = meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
%                 mask = ((x.^2+y.^2)<=r^2);              
% 
%                 % Find a threshold for the pixels inside the mask area.
%                 colors = square(:) .* uint8(mask(:));
%                 t = 255*graythresh(colors);              
% 
%                 figure(699), imshow(square,[])
%                 % Detect whether this is a white piece or black piece.
%                 
%                 color = -99;
%                 if mean(colors) > 30
%                     % Black piece
%                     color = 0;
%                     sqMask = (square < t) & mask;
%                     figure(23), imshow(sqMask, []);
%                 else
%                     % White piece
%                     color = 1;
%                     sqMask = (square > t) & mask;
%                     figure(23), imshow(sqMask, []);
%                 end
%                 
%                 
%                
%                 eta = SI_Moment(sqMask);
%                 inv_moments = Hu_Moments(eta);  
%                 all_moments = [all_moments; inv_moments];
%                 
%                 
%                 inv_moments = real(log(repmat(inv_moments,6,1)));
% %                 inv_moments = inv_moments(1:2);
% %                 inv_moments(2) = inv_moments(2) / inv_moments(1);
% %                 test_moment = inv_moments(2);
% 
%                 difference = [];
%                 if color == 1
%                     difference = abs(TEMPLATE1_BINARY_MOMENTS_WHITE(:,1:2) - inv_moments(:,1:2));
%                 else
%                     difference = abs(TEMPLATE1_BINARY_MOMENTS_BLACK(:,1:2) - inv_moments(:,1:2));
%                 end
%                 
%                 difference = sum(difference,2);
%                 match_id = find(difference == min(difference));
%                 
%                 perimSum = 0;
%                 areaSum = 0;
%                 extentSum = 0;
%                 majL = 0;
%                 minL = 0;
%                 eccSum = 0;
%                 cvArea = 0;
%                 %extrema = [];
%                 bbArea = 0;
% 
%                 
%                 S = strel('disk',3,0);
%                 
%                 if any(match_id == [1,3,4,5])
%                     closedSq = imclose(sqMask,S);
%                     figure(124), imshow(closedSq,[])
%                     
%                     [L,n] = bwlabel(closedSq);
%                     stats = regionprops(L,'Area','Perimeter','Centroid','BoundingBox', 'Eccentricity', ...
%                         'Extent','MajorAxisLength','MinorAxisLength','ConvexArea','ConvexImage');
%                     for l = 1:n
%                         perimSum = perimSum + stats(l).Perimeter;
%                         areaSum  = areaSum + stats(l).Area;
%                         extentSum = extentSum + stats(l).BoundingBox(4)*stats(l).BoundingBox(3);
%                         majL = majL + stats(l).MajorAxisLength;
%                         minL = minL + stats(l).MinorAxisLength;
%                         eccSum = eccSum + stats(l).Eccentricity;
%                         %extrema = [extrema; stats(l).Extrema];
%                         cvArea = cvArea + stats(l).ConvexArea;
%                         bbArea = bbArea + stats(l).BoundingBox(3)*stats(l).BoundingBox(4);
%                     end
%                     
%                     ratios = cvArea/areaSum
%                     %ratios = repmat(ratios,2,1);
%                     
%                     match_id
%                     
%                     testRats = [];
%                     if match_id == 4
%                         if color == 1
%                             testRats = [TEMPLATE_RATIOS(3); TEMPLATE_RATIOS(7)];
%                         else
%                             testRats = [TEMPLATE_RATIOS(4); TEMPLATE_RATIOS(8)];
%                         end
%                         
%                         ratDiffs = abs(ratios - testRats);
%                         %ratDiffs = sum(ratDiffs,2);
%                         new_id = find(ratDiffs == min(ratDiffs));
%                         if new_id == 1
%                             match_id = 2;
%                         end
%                     elseif match_id == 3
%                         testRats = [TEMPLATE_RATIOS(6); TEMPLATE_RATIOS(12)];
%                         
%                         ratDiffs = abs(ratios - testRats);
%                         %ratDiffs = sum(ratDiffs,2);
%                         new_id = find(ratDiffs == min(ratDiffs));
%                         if new_id == 2
%                             match_id = 6;
%                         end
%                     elseif match_id == 5
%                         testRats = [TEMPLATE_RATIOS(9); TEMPLATE_RATIOS(11)];
%                         
%                         ratDiffs = abs(ratios - testRats);
%                         %ratDiffs = sum(ratDiffs,2);
%                         new_id = find(ratDiffs == min(ratDiffs));
%                         if new_id == 2
%                             match_id = 6;
%                         end
%                     elseif match_id == 1
%                         testRats = [TEMPLATE_RATIOS(2); TEMPLATE_RATIOS(6)];
%                         
%                         ratDiffs = abs(ratios - testRats);
%                         %ratDiffs = sum(ratDiffs,2);
%                         new_id = find(ratDiffs == min(ratDiffs));
%                         if new_id == 2
%                             match_id = 3;
%                         end
%                     end   
% 
%                 end
%                 
%                 
%                 ix = 0; iy = 0;
%                 if color == 1
%                     iy = TEMPLATE1_BINARY_MOMENTS_WHITE(match_id,4);
%                     ix = TEMPLATE1_BINARY_MOMENTS_WHITE(match_id,5);
%                 else
%                     iy = TEMPLATE1_BINARY_MOMENTS_BLACK(match_id,4);
%                     ix = TEMPLATE1_BINARY_MOMENTS_BLACK(match_id,5);
%                 end
%                 
%                 pieceNames(iy,ix)
% 
%                 
%                 end
%                 
                
                


                % -----------------------------------------------------
                % This code finds the mask of the square that's used to
                % find image moments
                % -----------------------------------------------------
                
                % Gradient of the test image
                [squareGradient,~] = imgradient(square,'prewitt');
                
                % Finds the circle around the piece
                [pieceCenter,pieceRadius] = imfindcircles(squareGradient,[43 50],...
                    'Sensitivity',0.99, 'ObjectPolarity','bright', 'EdgeThreshold',0);
                pieceCenter = mean(pieceCenter,1);
                pieceRadius = mean(pieceRadius);

                % Creates a mask of the test square, with ones where the
                % piece is (using the circle just found)
                cx = pieceCenter(1);       cy = pieceCenter(2);
                ix = size(square,1);       iy = size(square,2);
                r = pieceRadius-5;
                [x,y] = meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
                mask = ((x.^2+y.^2)<=r^2);

                % For plotting masks
                figure(69)
                subplot(3,1,2)
                imshow(mask,[])
                title('Image mask')
                
                subplot(3,1,1)
                imshow(square,[])
                viscircles(pieceCenter,pieceRadius,'Color','r')
                title('Test  square')
                
                subplot(3,1,3)
                imshow(uint8(mask).*square,[])
                
                
                
                
                
                             
                % ---------------------------------------------------------
                % Method 1: Using NCC on gradient of test square and
                % template images (I did the same process with the original
                % images and it didn't work even as well as this
                %
                % Current success rate:   12/16 white pieces
                % (This doesn't into account guessing the incorrect color;
                % I've written code that can accurately determine the color
                % of a piece, so this isn't a problem.
                %
                % Also, regardless of the piece being tested, the NCC
                % always matches quite well with the queen, which causes
                % most black pieces to be matched incorrectly to the queen.
                % I've thought about doing this algorithm with a template
                % that excludes the queen and king pieces, so that this
                % isn't an issue. But I started working on Method 3
                % instead.
                %----------------------------------------------------------
                
                % Algorithm runs NCC over increasingly blurrier template
                % images
                sig = 0.5:0.5:3;
                
                % squareNew is the test square image scaled so that the
                % piece's size matches those in the template
                scaleUpRatio = mean(TEMPLATE_RADII) / pieceRadius;
                squareNew = uint8(squareGradient).*uint8(mask);
                
                subplot(3,1,3)
                imshow(squareNew,[])
                title('Masked gradient of test square')
                
                squareNew = imresize(squareNew,scaleUpRatio);
                
                
                pieceMaxes = zeros(length(sig),3);
                tempIt = 1;
                blurredTemp = imgaussfilt(TEMPLATE1_GRADIENT,0.01);
                for l = 1%:length(sig)       % Runcs NCC for templates of various blur levels
                    currentMax = 0;
                    for k = 0:10:360        % Runs NCC for all orientations of the test square image
                        rotated = imrotate(squareNew,k,'crop');
                        figure(55), imshowpair(rotated,blurredTemp,'montage');
                        
                        figure(66), imshow(rotated,[]);
                        title(sprintf('Test image against template,  at %d degrees',k))


                        tempIt = tempIt + 1;
                        c = normxcorr2(rotated,blurredTemp);
%                         figure(49), surf(c), shading flat
%                         xlabel('x')
%                         ylabel('y')
%                         zlabel('NCC')
%                         title('Normalized Cross-Correlation of test image with template')
                        
                        % Finds best correlation over all orientations of
                        % test square
                        if max(c(:)) > currentMax
                            currentMax = max(c(:));
                            [ypeak, xpeak] = find(c == currentMax);                          
                        end
                    end
                    % pieceMaxes stores best correlation and its location
                    % for each level of the gaussian filter used
                    pieceMaxes(l,1)   = currentMax;
                    pieceMaxes(l,2:3) = [ypeak,xpeak];
                
                    % Applies more blurring to the template
                    blurredTemp = imgaussfilt(TEMPLATE1_GRADIENT,sig(l));
                end 
                
                % Accounts for the shift in NCC
                pieceMaxes(:,3) = pieceMaxes(:,3) - size(squareNew,2)/2;
                pieceMaxes(:,2) = pieceMaxes(:,2) - size(squareNew,1)/2;

                % Shows each of the best correlations found on the template
                figure(2), imshowpair(TEMPLATE1,square,'montage')
                hold on
                plot(pieceMaxes(:,3),pieceMaxes(:,2),'r.','MarkerSize',12)
                hold off

                % Rough way of determining which piece the best
                % correlations fall onto (this assumes pieces in template
                % image are evenly spaced over image
                pieceWidth = size(TEMPLATE1,2) / 6;
                vertMid = size(TEMPLATE1,1) / 2;
                
                % Finds the corresponding index of best matches
                indx = ceil(pieceMaxes(:,3) / pieceWidth);
                indy = pieceMaxes(:,2) / vertMid > 1;   indy = indy + 1;
                
                % Determines which piece was correlated to most often
                indx_m = mode(indx);
                indy_m = mode(indy);
                
%                 % Used to display the matched piece name
%                 pieceNames = {'white pawn','white bishop','white knight','white rook','white queen','white king';
%                               'black pawn','black bishop','black knight','black rook','black queen','black king'};
%                 pieceNames(indy_m,indx_m)
%                 



                blah = blah + 1;


                
                
                
%                 % ------------------------------------------------------
%                 % Method 2: Matching piece with Hu's moments on original
%                 % test square and template images
%                 % Currently this method has a horrible success rate
%                 % ------------------------------------------------------
%                 
%                 figure(2), imshow(square,[])
%                 % Finds moments of the gradient of the test square
%                 % Uses the mask made earlier to specify which portion of
%                 % image to find moments for
%                 eta = SI_Moment(squareGradient,mask);
%                 inv_moments = Hu_Moments(eta);
%                 % Converts the moments to a useful form
%                 inv_moments = real(log(repmat(inv_moments,12,1)));
%                 
%                 % Compares the moments to the moments of each piece in the
%                 % template images gradient (with gaussian filter)
%                 diff = abs(TEMPLATE1_GRAD_MOMENTS_BLURRED - inv_moments);
%                 diff = sum(diff(:,1:4),2);      %Only compares the first 4 moments
%                 
%                 % Finds the best match
%                 ind = find( diff == min(diff) );
%                 figure(3), imshow(TEMPLATE1,[]);
%                 hold on
%                 plot(CORRESP_CENTERS(ind,1),CORRESP_CENTERS(ind,2),'r.','MarkerSize',30)
%                 hold off

                
                


                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
%                 % -------------------------------------------------------
%                 % Method 3: Using binary images of template and pieces to
%                 % match based on Hu's moments or on ratios of connected
%                 % components' properties
%                 % Works fairly well for the white piece, except pawns. I'm
%                 % thinking of running this for a few blurriness levels of
%                 % the binary template, but if you have any other ideas let
%                 % me know.
%                 % -------------------------------------------------------
% 
%                 % Converts test square to binary
%                 bwSquare = im2bw(square, 0.95*graythresh(square));
%                 
%                 % Finds piece outlines for either color
%                 [centers_d,radii_d] = imfindcircles(bwSquare,[30 60],'ObjectPolarity','dark');
%                 [centers_b,radii_b] = imfindcircles(bwSquare,[30 60],'ObjectPolarity','bright');
%                 
%                 % Stores all circles in one variable (should only be one
%                 % circle per piece)
%                 cents = [centers_d; centers_b];
%                 rads  = [radii_d; radii_b];
%                 
%                 % Displays circle on binary image
%                 figure(21), imshow(bwSquare,[])
%                 viscircles(cents,rads,'Color','g');
%                 
% 
%                 while true              
%                     % Finds the color of the piece by taking mode of binary
%                     % image (or of just the piece for cases in which a circle
%                     % was detected)
%                     [iy,ix] = size(bwSquare);
%                     polarity = [];
%                     if ~isempty(rads)
%                         for iy = 1:size(bwSquare,1)
%                             for ix = 1:size(bwSquare,2)
%                                 if (ix - cents(1,1))^2 + (iy - cents(1,2))^2 <= rads(1)^2
%                                     polarity = [polarity; bwSquare(iy,ix)];
%                                 end
%                             end
%                         end
%                         polarity = mode(mode(polarity));
%                     else
%                         polarity = mode(mode(bwSquare));
%                     end
% 
%                     % Finds blobs that are of the same color as the piece icon
%                     if polarity == 0
%                         [LBlobs,nBlobs] = bwlabel(bwSquare);
%                     else
%                         [LBlobs,nBlobs] = bwlabel(~bwSquare);
%                     end
%                     blobStats = regionprops(LBlobs,'Area','Centroid', ...
%                         'BoundingBox','Perimeter','FilledArea','Solidity');
% 
%                     wrongBlobs = [];
%                     % Removes blobs that are a large portion of the image
%                     for l = 1:nBlobs   
%                         bbArea = blobStats(l).BoundingBox(3)*blobStats(l).BoundingBox(4);
%                         if bbArea > 0.5*numel(bwSquare)            % CHANGE THIS TOLERANCE
%                             wrongBlobs = [wrongBlobs, l];
%                         end
%                     end
% 
%                     % Determines which of those blobs is nearest to the center
%                     % of the circle surrounding the piece (or, if no circle was
%                     % found, to the center of the image)
%                     diffs = zeros(nBlobs,1);
%                     if ~isempty(rads)
%                         for l = 1:nBlobs
%                             diffs(l) = norm(blobStats(l).Centroid - cents);
%                         end
%                     else
%                         for l = 1:nBlobs
%                             diffs(l) = norm(blobStats(l).Centroid - size(bwSquare)/2);
%                         end 
%                     end                         
%                     diffs(wrongBlobs) = NaN;
% 
%                     mindex = find(diffs == min(diffs));
% 
%                     % If the best blob found is tiny, run the process again
%                     % with a lower graythresh value to fill out the piece
%                     % icon
%                     % (Need to define these tolerance values elsewhere)
%                     if blobStats(mindex).Area < 0.01*numel(bwSquare)
%                         bwSquare = im2bw(square,0.75*graythresh(square));
%                     else 
%                         break;
%                     end
%                     
%                 end
%                 
%                              
%                 % Displays the bounding box of that blob
%                 rectangle('Position', blobStats(mindex).BoundingBox, ...
%                     'EdgeColor', 'r', 'Linewidth', 2);
% 
%                 
%                 
%                 
%                 % FOR HU'S MOMENTS
%                 bb = blobStats(mindex).BoundingBox;
%                 
%                 % Transfer the test blob to a square whose background is of
%                 % the same color as the piece's background. This is to get
%                 % rid of background noise
%                 % NOTE: Because I do this with the bounding box,
%                 % occasionaly additional blobs will be transferred over,
%                 % which could be negatively affecting the moment
%                 % calculations. There must be a way to only take the info
%                 % from one blob and transfer that, rather than taking the
%                 % whole rectangle contained in the bounding box. Either
%                 % this, or I could run bwlabel again and filter out the
%                 % blobs that got transferred over
%                 filteredSquare = zeros(size(bwSquare)) + polarity;
%                 filteredSquare( bb(2):bb(2)+bb(4)-1, bb(1):bb(1)+bb(3)-1 ) = ...
%                     bwSquare(bb(2):bb(2)+bb(4)-1, bb(1):bb(1)+bb(3)-1);
%                 
%                 figure(89), imshow(filteredSquare,[])    
%                 
%                 % Finds the moments of the filteredSquare
%                 eta = SI_Moment(filteredSquare,ones(size(filteredSquare)));
%                 inv_moments = Hu_Moments(eta);
%                 % Stores the moments in a useful way
%                 inv_moments = real(log(repmat(inv_moments,12,1)));
%                 
%                 % Compares the moments to those of each piece in the binary
%                 % template image 
%                 diff = abs(TEMPLATE1_BINARY_MOMENTS - inv_moments);
%                 diff = sum(diff,2);
%                 
%                 % Finds the best match and displays it on the template
%                 ind = find( diff == min(diff) );
%                 figure(54), imshow(TEMPLATE_BW,[]);
%                 hold on
%                 plot(CORRESP_CENTERS(ind,1),CORRESP_CENTERS(ind,2),'r.','MarkerSize',30)
%                 hold off
                                

            end
            pause
        end
    end
    
    % pause;
end