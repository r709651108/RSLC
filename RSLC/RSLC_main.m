
 close all; clc; clear;

I=imread('img/test.bmp');
GT=imread('img/gt.bmp');

tic
I=double(I);
I(I==0)=0.01;
[row,col]= size(I);

%% edge region smoothing
dW=7; %dection template width
sW=5; %smoothing template width
N1=4;  % iter
[dGauMeanI,NdirMap] = edgeRegionSM(I,dW,sW,N1);
toc
%% homogeneous regions smoothing
difSW=5;      %gaussian filter W
medW=5;       %median filter W
N2=2;         %iter
[Im,sigmaMap] = homoRegionSM(I,NdirMap,medW,difSW,N2);
toc
%% class label correction
C=5;        %clusters
W=11;      %search width
[finalI] = labelCorrection(dGauMeanI,Im,sigmaMap,C,W);
toc
tim=toc;
%  figure();imagesc(finalI);

%% calACC
[grayMap,rgbMap,Acc] = calAcc(finalI,GT);
sprintf('acc:%.2f, time:%.2f s ',Acc,tim)

figure();
subplot(2,3,1);imshow(uint8(I)); title('test image');
subplot(2,3,2);imshow(uint8(GT));title('gray groundtruth');
subplot(2,3,3);imagesc(GT); title('rgb groundtruth');
subplot(2,3,5);imshow(uint8(grayMap)); title('gray result');
subplot(2,3,6);imagesc(finalI); title('rgb result');




