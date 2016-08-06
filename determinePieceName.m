function determinePieceName(template1,finalMatches1,template2,finalMatches2)

    pieceCount = zeros(2,6);
    
    % I = imread('template1.jpg');
    % I = imresize(I,0.25);
    for k = 1
        if k == 1
            template = template1;
            finalMatches = finalMatches1;
        else
            template = template2;
            finalMatches = finalMatches2;
        end
    
        pieceWidth = size(template,2) / 6;
        vertMid = size(template,1) / 2;


%     finalMatches = [
%       408.3378  443.7082  446.6647  439.5700  181.9366  181.9366;
%        43.7286   44.1203   97.4600   94.2068  173.3499  173.3499
%     ];

%     figure, imshow(I,[]);
%     hold on
%     plot(finalMatches(1,:),finalMatches(2,:),'r.','MarkerSize',12);
%     hold off
    
        indx = ceil(finalMatches1(1,:) / pieceWidth);
        indy = finalMatches1(2,:) / vertMid > 1;  indy = indy + 1;

        for i = 1:size(finalMatches1,2)
            pieceCount(indy(i),indx(i)) = pieceCount(indy(i),indx(i)) + 1;
        end
    end

[maxy,maxx] = find(pieceCount == max(pieceCount(:)));

pieceNames = {'white pawn','white bishop','white knight','white rook','white queen','white king';
              'black pawn','black bishop','black knight','black rook','black queen','black king'};

pieceNames(maxy,maxx)
end