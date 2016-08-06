function [f1match, f2match] = findSiftMatches(I1,I2)
    if size(I1,3)>1; I1 = rgb2gray(I1); end
    I1 = single(I1);                      % Convert to single precision floating point
    if size(I2,3)>1; I2 = rgb2gray(I2); end
    I2 = single(I2);                      % Convert to single precision floating point

    % Find SIFT features for both images
    [f1,d1] = findSiftFeatures(I1);
    [f2,d2] = findSiftFeatures(I2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Threshold for matching
    % Descriptor D1 is matched to a descriptor D2 only if the distance d(D1,D2)
    % multiplied by THRESH is not greater than the distance of D1 to all other
    % descriptors
    thresh = 1.0; % default = 1.5; increase to limit matches
    [matches, ~] = vl_ubcmatch(d1, d2, thresh);
    fprintf('Number of matching frames (features): %d\n', size(matches,2));
    indices1 = matches(1,:); % Get matching features
    f1match = f1(:,indices1);
    d1match = d1(:,indices1);
    indices2 = matches(2,:);
    f2match = f2(:,indices2);
    d2match = d2(:,indices2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Show matches
    figure(799), imshow([I1,I2],[]);
    o = size(I1,2) ;
    line([f1match(1,:);f2match(1,:)+o], ...
        [f1match(2,:);f2match(2,:)]) ;
    for i=1:size(f1match,2)
        x = f1match(1,i);
        y = f1match(2,i);
        text(x,y,sprintf('%d',i), 'Color', 'r');
    end
    for i=1:size(f2match,2)
        x = f2match(1,i);
        y = f2match(2,i);
        text(x+o,y,sprintf('%d',i), 'Color', 'r');
    end
    title('Initial SIFT Matches, peak\_thresh = 4.5, edge\_thresh = 8')
end