clc;
clear all;

% compare images acquired with each quadrant of the BSD detector

display ('select images for quadrants');
[fnameQ1, pathname1] = uigetfile('*.tif','Select quadrant 1 image');
[fnameQ2, pathname2] = uigetfile('*.tif','Select quadrant 2 image');
[fnameQ3, pathname3] = uigetfile('*.tif','Select quadrant 3 image');
[fnameQ4, pathname4] = uigetfile('*.tif','Select quadrant 4 image');

% read images

Q1 = imread(fullfile(pathname1,fnameQ1));
Q2 = imread(fullfile(pathname2,fnameQ2));
Q3 = imread(fullfile(pathname3,fnameQ3));
Q4 = imread(fullfile(pathname4,fnameQ4));

% Normalize images

nQ1 = double(Q1);nQ1=nQ1-mean(nQ1(:));nQ1=nQ1/std(nQ1(:));
nQ2 = double(Q2);nQ2=nQ2-mean(nQ2(:));nQ2=nQ2/std(nQ2(:));
nQ3 = double(Q3);nQ3=nQ3-mean(nQ3(:));nQ3=nQ3/std(nQ3(:));
nQ4 = double(Q4);nQ4=nQ4-mean(nQ4(:));nQ4=nQ4/std(nQ4(:));

%plot Normalized images
figure;
subplot(221),hist(nQ1(:),255), Legend('Q1');
subplot(222),hist(nQ2(:),255), Legend('Q2');
subplot(223),hist(nQ3(:),255), Legend('Q3');
subplot(224),hist(nQ4(:),255), Legend('Q4');

figure; % histograms overlayed
[n1,xout1] = hist(nQ1(:),255);
[n2,xout2] = hist(nQ2(:),255);
[n3,xout3] = hist(nQ3(:),255);
[n4,xout4] = hist(nQ4(:),255);

plot(xout1,n1,'.r');
hold on
plot(xout2,n1,'.g');
plot(xout3,n1,'.b');
plot(xout4,n1,'.c');

legend ('Q1', 'Q2', 'Q3', 'Q4', 'Location','NorthEast');