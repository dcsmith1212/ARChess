clear all
run('~/Documents/MATLAB/CompVis/project/constants.m')

%%% CODE TO FIND BINARY TEMPLATE INDIVIDUAL PIECES

global TEMPLATE1            TEMPLATE1_MOMENTS...
       TEMPLATE1_GRADIENT   TEMPLATE1_GRAD_MOMENTS...
       TEMPLATE2            TEMPLATE2_GRADIENT...
       TEMPLATE_RADII       TEMPLATE1_GRAD_MOMENTS_BLURRED...
       TEMPLATE_CENTERS
   
w = size(TEMPLATE1,2);
h = size(TEMPLATE1,1);

template_moments = [];
closed_moments = [];
regionProps = [];

for i = 0:5
    for j = 0:1
        square = TEMPLATE1( round(j*h/2 + (1:h/2)) , round(i*w/6 + (1:w/6)) );
        figure(47), imshow(square,[])
        
        E = edge(square, 'Canny');
        sensitivity = 0.99;  % default is 0.85
        [centers_d,radii_d, metrics_d] = imfindcircles(E,[45 60], ...
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
        r = pieceRadius-1;
        [x,y] = meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
        mask = ((x.^2+y.^2)<=r^2);              

        % Find a threshold for the pixels inside the mask area.
        colors = square(:) .* uint8(mask(:));
        t = 255*graythresh(colors);              

        % Detect whether this is a white piece or black piece.
        if sum(colors) / sum(mask(:)) > t
            % White piece
            sqMask = (square < t) & mask;
            figure(23), imshow(sqMask, []);
        else
            % Black piece
            sqMask = (square > t) & mask;
            figure(23), imshow(sqMask, []);
        end
        
        if i == 4 && j == 0
            sqMask = (square > t) & mask;
            figure(23), imshow(sqMask, []);
        end
        
        eta = SI_Moment(sqMask);
        inv_moments = Hu_Moments(eta);
        template_moments = [template_moments; inv_moments, j+1, i+1];
        
        
        % For differentiating between rooks and bishops
        
        S = strel('disk',3,0);
        closedSq = imclose(sqMask,S);
        eta2 = SI_Moment(closedSq);
        inv_moments2 = Hu_Moments(eta);
        closed_moments = [closed_moments; inv_moments2, j+1, i+1];
        
        figure(89), imshow(closedSq,[])
        

        perimSum = 0;
        areaSum = 0;
        extentSum = 0;
        majL = 0;
        minL = 0;
        eccSum = 0;
        cvArea = 0;
        bbArea = 0;
        
        [L,n] = bwlabel(closedSq);
        
        stats = regionprops(L,'Area','Perimeter','Centroid','BoundingBox', ...
            'Eccentricity','Extent','MajorAxisLength','MinorAxisLength','ConvexArea');
        for l = 1:n
            perimSum = perimSum + stats(l).Perimeter;
            areaSum  = areaSum + stats(l).Area;
            majL = majL + stats(l).MajorAxisLength;
            minL = minL + stats(l).MinorAxisLength;
            eccSum = eccSum + stats(l).Eccentricity;
            cvArea = cvArea + stats(l).ConvexArea;
            bbArea = bbArea + stats(l).BoundingBox(3)*stats(l).BoundingBox(4);
        end
        
        regionProps = [regionProps; perimSum/areaSum, majL/minL];
        
        pause;
        

        end
      
    end
end

template_moments = [template_moments(:,1:3) template_moments(:,end-1:end)];