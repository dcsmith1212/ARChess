function mouseRegister(I,boardLayout,imgCenters)
% Takes the current board layout and displays peice moves if a user clicks
% on a specific square

% Opens image for clicking
figure(432), imshow(I);
title('Click a piece to display its moves.')

% Constantly registers user clicks
while true
    dists = zeros(length(imgCenters),1);
    
    % Stores the image points for the user's click
    [Ximg,Yimg] = ginput(1);  clickLocation = [Ximg, Yimg];
    
    % Determines which square the user clikced on by finding the closest
    % square center
    for i = 1:length(imgCenters)
        dists(i) = norm(imgCenters(i,:) - clickLocation);
    end
    dists = reshape(dists,[8 8]);
    [yMin, xMin] = find(dists == min(dists(:)));
    
    % Delets any existing figures
    if exist('h1','var'); delete(h1); end
    if exist('h2','var'); delete(h2); end
    if exist('h3','var'); delete(h3); end
    
    % Displays a cyan circle on the square that the user clicked
    hold on
    h1 = plot(boardLayout(yMin,xMin,1),boardLayout(yMin,xMin,2),...
        'co','MarkerSize',22);
    hold off
    
    % Determines the potential moves based on what the corresponding piece
    % ID is in boardLayout
    pieceId = boardLayout(yMin,xMin,4);
    switch pieceId
        case 1
            possibleMoves = [
                xMin + 1, yMin
                ];
                if xMin == 2
                    possibleMoves = [
                        possibleMoves;
                        xMin + 2, yMin
                        ];
                end
        case 2
            possibleMoves = [
                [xMin+(1:7)', yMin+(1:7)'];
                [xMin+(1:7)', yMin-(1:7)'];
                [xMin-(1:7)', yMin+(1:7)'];
                [xMin-(1:7)', yMin-(1:7)'];
                ];
            
        case 3
            possibleMoves = [
                xMin + 1, yMin - 2;
                xMin + 1, yMin + 2;
                xMin - 1, yMin - 2;
                xMin - 1, yMin + 2;
                xMin + 2, yMin - 1;
                xMin + 2, yMin + 1;
                xMin - 2, yMin - 1;
                xMin - 2, yMin + 1;
                ];
                
        case 4
            possibleMoves = [
                [xMin*ones(7,1), yMin+(1:7)'];
                [xMin*ones(7,1), yMin-(1:7)'];
                [xMin+(1:7)', yMin*ones(7,1)];
                [xMin-(1:7)', yMin*ones(7,1)];
                ];
            
        case 5
            possibleMoves = [
                [xMin*ones(7,1), yMin+(1:7)'];
                [xMin*ones(7,1), yMin-(1:7)'];
                [xMin+(1:7)', yMin*ones(7,1)];
                [xMin-(1:7)', yMin*ones(7,1)];
                [xMin+(1:7)', yMin+(1:7)'];
                [xMin+(1:7)', yMin-(1:7)'];
                [xMin-(1:7)', yMin+(1:7)'];
                [xMin-(1:7)', yMin-(1:7)'];
                ];
            
        case 6
            possibleMoves = [
                xMin,     yMin + 1;
                xMin,     yMin - 1;
                xMin + 1, yMin + 1;
                xMin + 1, yMin;
                xMin + 1, yMin - 1;
                xMin - 1, yMin + 1;
                xMin - 1, yMin;
                xMin - 1, yMin - 1;
                ];
            
        otherwise

    end
    
    if exist('possibleMoves','var')
        valid = find( possibleMoves(:,1)>=1 &...
                      possibleMoves(:,1)<=8 &...
                      possibleMoves(:,2)>=1 &...
                      possibleMoves(:,2)<=8);
        possibleMoves = possibleMoves(valid,:);
        
        enemies = [];
        validMoves = sub2ind([8 8],possibleMoves(:,2),possibleMoves(:,1));
        otherPieces = find(boardLayout(:,:,4));
        collisions = intersect(otherPieces, validMoves);
        
        % Removes pieces that current peice collides with
        invalidLocs = ismember(validMoves,collisions);
        c = find(invalidLocs == 1);
        validMoves(c) = [];
        
        % Removes piece behind collision pieces
        [yV,xV] = ind2sub([8 8],validMoves);
        [yC,xC] = ind2sub([8 8],collisions);
        
        % Removes square from validMoves if they collide with potential
        % moves
        for i = 1:length(collisions)
            for k = 1:length(validMoves)
                p1 = [yMin   xMin];
                p2 = [yC(i) xC(i)];
                p3 = [yV(k) xV(k)];
                if collinear(p1,p2,p3) && p2InMiddle(p1,p2,p3)
                    validMoves(k) = -99;
                end
            end

            if boardLayout(yMin,xMin,3) ~= boardLayout(yC(i),xC(i),3)
                enemies = [
                    enemies;
                    collisions(i);
                    ];
            end
        end
        
        % Marks colliding enemy pieces separately
        [yE,xE] = ind2sub([8 8],enemies);
        
        for i = 1:length(collisions)
            for k = 1:length(enemies)
                p1 = [yMin   xMin];
                p2 = [yC(i) xC(i)];
                p3 = [yE(k) xE(k)];
                if collinear(p1,p2,p3) && p2InMiddle(p1,p2,p3)
                    enemies(k) = -99;
                end
            end
        end
        
        validMoves(find(validMoves == -99)) = [];
        
        % Plots green circles for valid moves and red if the move is to
        % take an enemy
        hold on
        h2 = plot( imgCenters(validMoves(:),1), imgCenters(validMoves(:),2),...
            'g.','MarkerSize',27);
        h3 = plot( imgCenters(enemies(:),1), imgCenters(enemies(:),2),...
            'r.','MarkerSize',29);
        hold off
        
        clearvars 'possibleMoves'
    end
end

% Determines if three points are collinear
function tf = collinear(p1,p2,p3)
mat = [p1(1)-p3(1) p1(2)-p3(2); ...           
        p2(1)-p3(1) p2(2)-p3(2)];
tf = det(mat) == 0;
end

% Determines if p2 is in between p1 and p8
function tf = p2InMiddle(p1,p2,p3)
    n1 = norm(p1-p2) + norm(p1-p3);
    n2 = norm(p2-p1) + norm(p2-p3);
    n3 = norm(p3-p1) + norm(p3-p2);
    tf = n2 < n1 && n2 < n3;
end

end