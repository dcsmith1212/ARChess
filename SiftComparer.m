function [correctSpots] = SiftComparer(square,ref)
    square_ext = uint8(zeros(size(ref,1),size(ref,2),3));
    square = imresize(square, 'outputSize', [size(ref,1) NaN]);
    square_ext(:,1:size(square,2),:) = square;

    [f1match,f2match] = findSiftMatches(square_ext,ref);

    sigma = 1;                    % uncertainty of image points at lowest scale
    maxIterations = 1000;           % don't do more than this many iterations
    Ps = 0.95;                      % desired confidence level (prob of success); value from 0..1
                
    [tform_1_2, indices, rmsErr] = fitHomographyRansac(...
        f1match, f2match, ...
        size(ref,1), size(ref,2), ...
        sigma, maxIterations, Ps);


    correctSpots = f2match(1:2,find(indices));
    figure(69), imshow(ref,[]);
    hold on
    plot(correctSpots(1,:),correctSpots(2,:),'r.','MarkerSize',12)
    hold off
    title('Correct SIFT matches (using RANSAC)')
end