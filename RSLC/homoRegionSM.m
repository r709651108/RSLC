function [Im,sigmaMap] = homoRegionSM(I,sigmaMap,medW,difSW,N2)  
%Gaussian + Median filter

    sigmaMap2=sigmaMap.^2;    
    Im=I;
    for k=1:N2
        [difSigI]=difSigmaFilter(sigmaMap2,Im,difSW);
        [Im]=medfilt2(difSigI,[medW,medW]);
    end     
    
end


function [difSigI]=difSigmaFilter(sigmaMap,I, W)
% use diff sigma gaussian filter    
    r=(W-1)/2;
    [xx,yy] = meshgrid(-r:r);
    [row,col] = size(I);

    img_pad = padarray(I,[r r],'symmetric');
    difSigI = zeros(row,col);
    for i = 1:row
        for j = 1:col
            x=r+i;
            y=r+j;
            cut = img_pad(x-r:x+r,y-r:y+r);
            sigma=sigmaMap(i,j)+0.1;
            w = exp(-(xx.^2+yy.^2)/(2*sigma^2));
            s = cut.*w;
            difSigI(i,j) = sum(s(:))/sum(w(:));
        end
    end
end