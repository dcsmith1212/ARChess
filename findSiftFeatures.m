function [f,d] = findSiftFeatures(I)
    if size(I,3)>1; I = rgb2gray(I); end
    I = single(I);                      % Convert to single precision floating point

    % These parameters limit the number of features detected
    peak_thresh = 4.5;                    % increase to limit; default is 0
    edge_thresh = 8;                   % decrease to limit; default is 10
    [f,d] = vl_sift(I, ...
        'PeakThresh', peak_thresh, ...
        'EdgeThresh', edge_thresh );
    fprintf('Number of frames (features) detected: %d\n', size(f,2));

%     % Show all SIFT features detected
%     h = vl_plotframe(f);
%     set(h,'color','y','linewidth',2);
end