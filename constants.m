global TEMPLATE1...
       TEMPLATE1_GRADIENT...
       TEMPLATE_RADII...
       TEMPLATE_CENTERS     TEMPLATE1_BINARY_MOMENTS_WHITE...
       TEMPLATE_BW...
       TEMPLATE1_BINARY_MOMENTS_BLACK   TEMPLATE_RATIOS_WHITE...
       TEMPLATE_RATIOS_BLACK

TEMPLATE1_orig = imread('template1.jpg');
TEMPLATE1_orig = imresize(TEMPLATE1_orig,0.25);
TEMPLATE1 = rgb2gray(TEMPLATE1_orig);
[TEMPLATE1_GRADIENT,~] = imgradient(TEMPLATE1,'prewitt');

TEMPLATE2_orig = imread('template2.jpg');
TEMPLATE2_orig = imresize(TEMPLATE2_orig,0.25);
TEMPLATE2 = rgb2gray(TEMPLATE2_orig);
[TEMPLATE2_GRADIENT,~] = imgradient(TEMPLATE2,'prewitt');

TEMPLATE1_BLURRED = imgaussfilt(TEMPLATE1,2.5);
TEMPLATE_BW = im2bw(TEMPLATE1_BLURRED,graythresh(TEMPLATE1));

% Centers and radii of 12 pieces in the template image
TEMPLATE_CENTERS = [
    182.8081  180.6276;
    560.1821  187.4813;
    306.8586  187.8617;
    431.5000  184.5404;
     60.3374  188.0670;
    685.0000  190.5480;
    554.7895   57.9027;
    687.5255   61.9312;
     60.2908   60.7179;
    305.0062   57.9827;
    184.3244   64.4859;
    432.0982   63.1064
    ];

TEMPLATE_RADII = [
    58.2265;
    59.7275;
    58.4230;
    58.3462;
    57.6622;
    57.9711;
    57.0657;
    62.9484;
    57.5575;
    57.1582;
    58.7638;
    57.0462
    ];

% % Code for creating template moments
% TEMPLATE1_MOMENTS = zeros(12,7);
% for k = 1:length(TEMPLATE_RADII)
%     cx = TEMPLATE_CENTERS(k,1);   cy = TEMPLATE_CENTERS(k,2);
%     [iy,ix] = size(TEMPLATE1);
%     r = TEMPLATE_RADII(k);
%     [x,y] = meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
%     mask = ((x.^2+y.^2)<=r^2);
% 
%     eta = SI_Moment(TEMPLATE1,mask);
%     TEMPLATE1_MOMENTS(k,:) = Hu_Moments(eta);
% end

TEMPLATE1_BINARY_MOMENTS_WHITE = [
   0.204003740806440   0.006491636155188   0.002867988685693 ...
   1.000000000000000   1.000000000000000;
   0.336660818521282   0.021256347270486   0.004310362118452 ...
   1.000000000000000   2.000000000000000;
   0.205095478843948   0.003245716401621   0.000336170868333 ...
   1.000000000000000   3.000000000000000;
   0.292356522555595   0.021554796108091   0.000308535985715 ...
   1.000000000000000   4.000000000000000;
   0.231846875696146   0.000273648598620   0.000809686179691 ...
   1.000000000000000   5.000000000000000;
   0.342782019278974   0.000157987559688   0.000090908453714 ...
   1.000000000000000   6.000000000000000
   ];

TEMPLATE1_BINARY_MOMENTS_BLACK = [
   0.207642535774075   0.006901635018905   0.003689808887304 ...
   2.000000000000000   1.000000000000000;
   0.306137175249630   0.012826824462872   0.010106251653722 ...
   2.000000000000000   2.000000000000000;
   0.250557677080537   0.005387910213969   0.000580275115315 ...
   2.000000000000000   3.000000000000000;
   0.297434134166280   0.013108201325016   0.000415417588049 ...
   2.000000000000000   4.000000000000000;
   0.271069298327684   0.000136212249348   0.001041296589194 ...
   2.000000000000000   5.000000000000000;
   0.285713294859744   0.000233634679513   0.000214476253981 ...
   2.000000000000000   6.000000000000000
   ];

TEMPLATE1_BINARY_MOMENTS_WHITE(:,1:3) = ...
    real(log(TEMPLATE1_BINARY_MOMENTS_WHITE(:,1:3)));
TEMPLATE1_BINARY_MOMENTS_BLACK(:,1:3) = ...
    real(log(TEMPLATE1_BINARY_MOMENTS_BLACK(:,1:3)));

TEMPLATE_RATIOS_WHITE = [
   1.103286384976526;
   1.727065959059894;
   1.212885154061625;
   1.397714907508161;
   1.231809932354397;
   1.712087912087912
   ];

TEMPLATE_RATIOS_BLACK = [
   1.120000000000000;
   1.446002076843198; 
   1.168270944741533;
   1.314007782101167;
   1.461387691101529;
   1.213091922005571
   ];