function [dGauMeanI,NdirMap] = edgeRegionSM(I,dW,sW,N)
%edge region smoothing  

    dirNum=8; %direction number
    [dKerList]=getDetectTemplate(dW,dirNum);
    [sKerList]=getSmoothTemplate(sW,dirNum);

    [row,col]=size(I);
    NdirMap=zeros(row,col,N);    
    NimgMap=zeros(row,col,N);
    for m=1:N
        [dirMap]=detecDirection(I,dKerList);
        [I]= edgeSmoothing(I,dirMap,sKerList,sW);
        NdirMap(:,:,m)=dirMap; 
        NimgMap(:,:,m)=I;
    end
    dGauMeanI=I;
    
%     [showEdgMap]= showTemplate(dKerList,dW,dWNumber);
%     [showEdgMap]= showTemplate(sKerList,sW,sWNumber);
    
end

%% smoothing
function [DirGauMap]= edgeSmoothing(I,dirMap,sKerList,sW)
    [row,col]=size(I);   
    [~,~,len]=size(sKerList);
    
    DirGauMap=zeros(row,col);
    ConvMap=zeros(row,col,len);
    gausFilter=fspecial('gaussian',[sW, sW],2);
    
    for n=1:len
        kernel=sKerList(:,:,n).*gausFilter;  
        kernel=kernel/sum(kernel(:));
        ConvMap(:,:,n)=imfilter(I,kernel, 'replicate'); 
    end    
    for i=1:row
        for j=1:col
            sel=dirMap(i,j);
            DirGauMap(i,j)=ConvMap(i,j,sel);
        end
    end

end


%% detection
function [dirMap]=detecDirection(I,decKerList)
    [row,col]=size(I);    
    [~,~,len]=size(decKerList);
    
    ConvMap=zeros(row,col,len);
    for n=1:len
        kernel=decKerList(:,:,n);
        w=sum(abs(kernel(:)))/2;
        kernel=kernel/w;
        ConvMap(:,:,n)=imfilter(I,kernel, 'replicate'); 
    end        

    dirMap =zeros(row,col);   
    for x=1:row
        for y=1:col
            Dir=-1;
            Val=0;
            for k=1:len
                P=abs(ConvMap(x,y,k));
                if P>=Val
                   Val=P;
                   Dir=k;
                end
            end
            dirMap(x,y)=Dir;           
        end        
    end
    
end

 %% construct smoothing templates 
function [sKerList]=getSmoothTemplate(sW,dirNum)

    sR=(sW-1)/2;
    sKerList=zeros(sW,sW,dirNum); 
    edgMapTemp=zeros(sW*2+1,sW*2+1);
    edgMapTemp(:,sW+1)=2;  
    
    for k=0:dirNum
        theta=k*(180/dirNum);
        edgMap=imrotate(edgMapTemp, theta-45,'nearest','crop'); 
        sKerList(:,:,k+1)=edgMap(sW+1-sR:sW+1+sR,sW+1-sR:sW+1+sR); 
    end
    
end

%% construct detection templates
function [decKerList]=getDetectTemplate(dW,dirNum)

    dR=(dW-1)/2;
    decKerList=zeros(dW,dW,dirNum*3); 
    
    %first template
    edgMapTemp=zeros(dW*2+1,dW*2+1);
    edgMapTemp(1:dW+1,1:dW+1)=1; 
    edgMapTemp(dW+1:end,dW+1:end)=-1;
    edgMapTemp(dW+1,dW+1)=0;
    
    for k=0:dirNum-1
        theta=k*(180/dirNum);
        edgMap=imrotate(edgMapTemp, theta,'nearest','crop'); 
        decKerList(:,:,k+1)=edgMap(dW+1-dR:dW+1+dR,dW+1-dR:dW+1+dR); 
    end    
    
end


%% show templates
function [showEdgMap]= showTemplate(kernelList,W,dirNumber)
    
    numW=dirNumber;          %direction number
    colN=8;                  %set show col number
    rowN=ceil(numW/colN);    %set row number
    
    addW=2;
    maxW=W+addW;
    showEdgMap=zeros(maxW*rowN,maxW*colN);  %show map

    EdgTempList=zeros(maxW,maxW,numW);
    
    pw=(maxW-W)/2;
    for k=1:dirNumber
        edgCut=kernelList(:,:,k);
        EdgTempList(:,:,k)=padarray(edgCut,[pw,pw],200);
    end
    
    cnt=0;
    for i=0:rowN-1
        for j=0:colN-1
            cnt=cnt+1;
            showEdgMap(i*maxW+1:i*maxW+maxW,j*maxW+1:j*maxW+maxW)=(EdgTempList(:,:,cnt)+1)*100;
        end
    end
    
    edgMap=kron(showEdgMap,ones(5,5));       
    figure();imshow(uint8(edgMap));

end