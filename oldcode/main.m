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
       TEMPLATE_BW          TEMPLATE_RATIOS_WHITE...
       TEMPLATE1_BINARY_MOMENTS_BLACK    TEMPLATE_RATIOS_BLACK
   
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
    I = imread('test6.jpg');        % Or open single frame
    
    
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
    
    figure, imshow(I);
    line([corners(1,1) corners(2,1)], [corners(1,2) corners(2,2)], 'Color', 'g', 'Linewidth', 2);
    line([corners(1,1) corners(4,1)], [corners(1,2) corners(4,2)], 'Color', 'g', 'Linewidth', 2);
    line([corners(3,1) corners(2,1)], [corners(3,2) corners(2,2)], 'Color', 'g', 'Linewidth', 2);
    line([corners(3,1) corners(4,1)], [corners(3,2) corners(4,2)], 'Color', 'g', 'Linewidth', 2);


    
    orthoPts = [
        1, 1;
        N, 1;
        N, N;
        1, N
        ];
    

    % Compute transform, from corresponding control points
    Tform1 = fitgeotrans(corners,orthoPts,'projective');
    % Transform input image to output image
    ref2Doutput = imref2d( ...
        [N, N], [0,N], [0,N]);
    Iout1 = imwarp(I,Tform1,'OutputView',ref2Doutput);

    % Transforms the intersection points found from findCheckBoard to the
    % orthophoto's coordinates
    [ortho_intersect_x, ortho_intersects_y] = transformPointsForward(...
        Tform1, intersections(:,1), intersections(:,2));
    [ortho_centers_x, ortho_centers_y] = transformPointsForward(Tform1, squareCenters(:,1), squareCenters(:,2));
    


    for i = 1:length(ortho_intersect_x)-9
        if mod(i,9) ~= 0
            x0 = ortho_intersect_x(i);  xf = ortho_intersect_x(i+10);
            y0 = ortho_intersects_y(i);  yf = ortho_intersects_y(i+10);
            %rectangle('Position', [x0 y0 xf-x0 yf-y0], 'EdgeColor', 'g', 'LineWidth', 2);          
            
            square = Iout1(y0:yf,x0:xf);
            sub = square(end/4:3*end/4,end/4:3*end/4);
            %square = square(end/4:3*end/4,end/4:3*end/4);

            if max(max(sub)) - min(min(sub)) > 100  % For those squares that have a piece on them


                % Hoff's Code
                E = edge(square, 'Canny');
                sensitivity = 0.98;  % default is 0.85
                [centers_d,radii_d, metrics_d] = imfindcircles(E,[40 48], ...
                'Sensitivity', sensitivity, ...
                'ObjectPolarity','dark');

                [~,ind] = max(metrics_d);   % Pick the strongest circle
                figure(21), imshow(E, []);
                if ~isempty(ind) 
                    % Draw on top of the subimage.
                    viscircles(centers_d(ind,:),radii_d(ind),'Color','r');
                end

 
                
                if ~isempty(ind)
                pieceCenter = centers_d(ind,:);
                pieceRadius = radii_d(ind);
                % Creates a mask of the test square, with ones where the
                % piece is (using the circle just found)
                cx = pieceCenter(1);       cy = pieceCenter(2);
                ix = size(square,1);       iy = size(square,2);
                r = pieceRadius-5;
                [x,y] = meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
                mask = ((x.^2+y.^2)<=r^2);              

                % Find a threshold for the pixels inside the mask area.
                colors = square(:) .* uint8(mask(:));
                t = 255*graythresh(colors);              

                figure(699), imshow(square,[])
                % Detect whether this is a white piece or black piece.
                
                color = -99;
                if mean(colors) > 30
                    % Black piece
                    color = 0;
                    sqMask = (square < t) & mask;
                    figure(23), imshow(sqMask, []);
                else
                    % White piece
                    color = 1;
                    sqMask = (square > t) & mask;
                    figure(23), imshow(sqMask, []);
                end
                
                
               
                eta = SI_Moment(sqMask);
                inv_moments = Hu_Moments(eta);      
                inv_moments = real(log(repmat(inv_moments,6,1)));
%                 inv_moments = inv_moments(1:2);
%                 inv_moments(2) = inv_moments(2) / inv_moments(1);
%                 test_moment = inv_moments(2);

            
                
                
                % For ratio matching
                perimSum = 0;
                areaSum = 0;
                extentSum = 0;
                majL = 0;
                minL = 0;
                eccSum = 0;
                cvArea = 0;
                %extrema = [];
                bbArea = 0;

                
                S = strel('disk',3,0);
                closedSq = imclose(sqMask,S);
                figure(124), imshow(closedSq,[])
                title('closed test square')

                [L,n] = bwlabel(closedSq);
                stats = regionprops(L,'Area','Perimeter','Centroid','BoundingBox', 'Eccentricity', ...
                    'Extent','MajorAxisLength','MinorAxisLength','ConvexArea','ConvexImage');
                for l = 1:n
                    perimSum = perimSum + stats(l).Perimeter;
                    areaSum  = areaSum + stats(l).Area;
                    extentSum = extentSum + stats(l).BoundingBox(4)*stats(l).BoundingBox(3);
                    majL = majL + stats(l).MajorAxisLength;
                    minL = minL + stats(l).MinorAxisLength;
                    eccSum = eccSum + stats(l).Eccentricity;
                    %extrema = [extrema; stats(l).Extrema];
                    cvArea = cvArea + stats(l).ConvexArea;
                    bbArea = bbArea + stats(l).BoundingBox(3)*stats(l).BoundingBox(4);
                end

                ratios = cvArea/areaSum;
                ratios = repmat(ratios,6,1);
                
                
                
                A = 0.25;
                B = 1 - A;
                
                momentDiffs = [];
                ratioDiffs = [];
                if color == 1
                    momentDiffs = abs(TEMPLATE1_BINARY_MOMENTS_WHITE(:,1:2) - inv_moments(:,1:2));
                    ratioDiffs = abs(TEMPLATE_RATIOS_WHITE - ratios);
                else
                    momentDiffs = abs(TEMPLATE1_BINARY_MOMENTS_BLACK(:,1:2) - inv_moments(:,1:2));
                    ratioDiffs = abs(TEMPLATE_RATIOS_BLACK - ratios);
                end
                
                momentDiffs = sum(momentDiffs,2);
                finalDiffs  = A*momentDiffs + B*ratioDiffs;              
                matchID = find(finalDiffs == min(finalDiffs));
                
          
                ix = 0; iy = 0;
                if color == 1
                    iy = TEMPLATE1_BINARY_MOMENTS_WHITE(matchID,4);
                    ix = TEMPLATE1_BINARY_MOMENTS_WHITE(matchID,5);
                else
                    iy = TEMPLATE1_BINARY_MOMENTS_BLACK(matchID,4);
                    ix = TEMPLATE1_BINARY_MOMENTS_BLACK(matchID,5);
                end
                
                pieceNames(iy,ix)

                end
                
                pause
            end                
        end
    end
end