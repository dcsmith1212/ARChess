clear variables
close all
imtool close all

warning off MATLAB:colon:nonIntegerIndex
warning off images:initSize:adjustingMag

% Runs file to define global constants
run('~/Documents/MATLAB/CompVis/project/constants.m')
   
% Just dealing with a single picture here
nFrames = 1;

% The size of orthophoto of the board
N = 1000;

% Stores the moments for each piece on the board
all_moments = [];

% Stores the locations of the centers of the board in the image plane (from
% the input image)
imgCenters = zeros(64,2);

% Stores piece names and IDs for later identification
pieceNames(:,:,1) = {'white pawn','white bishop','white knight', ...
                            'white rook','white queen','white king';
                     'black pawn','black bishop','black knight', ...
                            'black rook','black queen','black king'};
pieceNames(:,:,2) = {[1, 1], [1, 2], [1, 3], [1, 4], [1, 5], [1, 6];
                     [2, 1], [2, 2], [2, 3], [2, 4], [2, 5], [2, 6]};
                          
% Will store locations of pieces amd their ID at each frame
boardLayout = zeros(8,8,5);
  
tic
for iFrame=1:10:nFrames
    I = imread('test12.jpg');
    
    % Reduce image size; is faster and we don't need full size to find
    % board.
    if size(I,2)>1280
        I = imresize(I, 1280/size(I,2));
    end

    % Find the checkerboard. Return the four outer corners as a 4x2 array,
    % in the form [ [x1,y1]; [x2,y2]; ... ].
    [corners, nMatches, avgErr, intersections] = findCheckerBoard(I);
    
    % Displays green lines around board border
    figure, imshow(I);
    line([corners(1,1) corners(2,1)], [corners(1,2) corners(2,2)], ...
        'Color', 'g', 'Linewidth', 2);
    line([corners(1,1) corners(4,1)], [corners(1,2) corners(4,2)], ...
        'Color', 'g', 'Linewidth', 2);
    line([corners(3,1) corners(2,1)], [corners(3,2) corners(2,2)], ...
        'Color', 'g', 'Linewidth', 2);
    line([corners(3,1) corners(4,1)], [corners(3,2) corners(4,2)], ...
        'Color', 'g', 'Linewidth', 2);

    % Displays red dots at centers of each square
    hold on
    k = 1;
    for i = 1:length(intersections)-9
        if mod(i,9) ~= 0
            avg = (intersections(i,:) + intersections(i+1,:) ...
                + intersections(i+9,:) + intersections(i+10,:))/4;
            plot(avg(1), avg(2), '.r', 'MarkerSize', 6);
            imgCenters(k,:) = avg;
            k = k + 1;
        end
    end
    hold off
    
    % Stores the centers of the square in boardLayout
    boardLayout(:,:,1) = reshape(imgCenters(:,1),[8 8]);
    boardLayout(:,:,2) = reshape(imgCenters(:,2),[8 8]);

    % Defines boundary of orthophoto
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
    

    for i = 1:length(ortho_intersect_x)-9
        if mod(i,9) ~= 0
            x0 = ortho_intersect_x(i);  xf = ortho_intersect_x(i+10);
            y0 = ortho_intersects_y(i);  yf = ortho_intersects_y(i+10);       
            
            square = Iout1(y0:yf,x0:xf);
            sub = square(end/4:3*end/4,end/4:3*end/4);

            if max(max(sub)) - min(min(sub)) > 100  % For those squares that 
                                                    % have a piece on them

                figure(99), imshow(square,[])
                
                % Hoff's Code
                E = edge(square, 'Canny');
                sensitivity = 0.98;  % default is 0.85
                [centers_d,radii_d, metrics_d] = imfindcircles(E,[40 48], ...
                'Sensitivity', sensitivity, ...
                'ObjectPolarity','dark');

                [~,ind] = max(metrics_d);   % Pick the strongest circle
 
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

                % Detect whether this is a white piece or black piece.
                color = -99;
                if mean(colors) > 30
                    % Black piece
                    color = 0;
                    sqMask = (square < t) & mask;
                else
                    % White piece
                    color = 1;
                    sqMask = (square > t) & mask;
                end
               
                % Finds Hu's moments
                eta = SI_Moment(sqMask);
                inv_moments = Hu_Moments(eta);      
                
                % Stores moments
                all_moments = [all_moments; inv_moments(1:3)];
                
                % Converts moments to a more useful format
                inv_moments = real(log(repmat(inv_moments,6,1)));
        
                % For ratio matching
                areaSum = 0;
                cvArea = 0;

                % Structuring element for image closing
                S = strel('disk',3,0);
                closedSq = imclose(sqMask,S);

                % Finds solidity value for blob of piece
                [L,n] = bwlabel(closedSq);
                stats = regionprops(L,'Area','Perimeter','Centroid', ...
                    'ConvexArea','ConvexImage','BoundingBox');
                for l = 1:n
                    areaSum  = areaSum + stats(l).Area;
                    cvArea = cvArea + stats(l).ConvexArea;
                end
                
                % Stores solidity value
                cvRatios = cvArea/areaSum;
                cvRatios = repmat(cvRatios,6,1);              
                
                % Weights for Hu's moments and solidity, respectively
                A = 1/4;
                B = 3/4;
                
                % Stores Hu's moments and solidity differences,
                % respectively 
                momentDiffs = [];
                cvDiffs = [];
                
                % Only searches for matches between pieces of the same
                % color as the unknown
                if color == 1
                    momentDiffs = abs(TEMPLATE1_BINARY_MOMENTS_WHITE(:,1:3) ...
                        - inv_moments(:,1:3));
                    cvDiffs = abs(TEMPLATE_RATIOS_WHITE - cvRatios);
                else
                    momentDiffs = abs(TEMPLATE1_BINARY_MOMENTS_BLACK(:,1:3) ...
                        - inv_moments(:,1:3));
                    cvDiffs = abs(TEMPLATE_RATIOS_BLACK - cvRatios);
                end
                
                % Sums the differences for the first three moments
                momentDiffs = sum(momentDiffs,2);
                
                % Weighted average of the moment and solidity differences
                finalDiffs  = A*momentDiffs + B*cvDiffs;         
                matchID = find(finalDiffs == min(finalDiffs));
                
                % Finds the location of the matched piece in the cell array
                % named pieceNames
                ix = 0; iy = 0;
                if color == 1
                    iy = TEMPLATE1_BINARY_MOMENTS_WHITE(matchID,4);
                    ix = TEMPLATE1_BINARY_MOMENTS_WHITE(matchID,5);
                else
                    iy = TEMPLATE1_BINARY_MOMENTS_BLACK(matchID,4);
                    ix = TEMPLATE1_BINARY_MOMENTS_BLACK(matchID,5);
                end
                
                % Prints the matched piece and stores the match
                strName = pieceNames{iy,ix};
                sprintf('Piece detected: %s', strName)

                end
                
                % Converts the location of the test square back to 8x8 grid
                % location
                xl = (x0 + xf)/2;
                yl = (y0 + yf)/2;
                piece_x = ceil(8 * xl / size(Iout1,2));
                piece_y = ceil(8 * yl / size(Iout1,1));
                pieceID = pieceNames{iy,ix,2};
                
                % Stores the identified piece in boardLayout for further
                % processing
                boardLayout(piece_y,piece_x,3:4) = pieceID;
                
            end                
        end
    end
    toc
    
    pause
    
    % Once all the pieces have been identified, open up the function that
    % registers user clicks and displays piece moves
    mouseRegister(I,boardLayout,imgCenters);
    
end