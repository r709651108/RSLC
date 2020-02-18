function [CLabelMap,edgCLabelMap2,finalI] = labelCorrection(dirGauMeanI,medI,sigmaMap,C,W)
     
    difSigMedI=(medI.*sigmaMap+dirGauMeanI)./(sigmaMap+1); % combine
    [CLabelMap] = mKMeans(single(difSigMedI),C); 
    
    edgI=edge(difSigMedI,'canny');
    edgCLabelMap=CLabelMap;
    edgCLabelMap(edgI~=0)=Inf;
    
    fr=(W-1)/2;
    [row,col]=size(medI);
    edgCLabelMap2=edgCLabelMap;
    [corrMap]= getCorrMap(edgCLabelMap,CLabelMap,W);

    for i=1:row
        for j=1:col
            if edgCLabelMap(i,j)~=Inf && corrMap(i,j)
                X0= max(1,i-fr);
                Xn= min(row,i+fr);
                Y0= max(1,j-fr);
                Yn= min(col,j+fr);
                [mostLabel]=findMostSimilarLabel(edgCLabelMap(X0:Xn,Y0:Yn),C,i-X0+1,j-Y0+1);
                edgCLabelMap2(i,j)=mostLabel;
            end
        end
    end
    [finalI]=transEdgToLabel(edgCLabelMap2,dirGauMeanI);

end

function [resMap]=transEdgToLabel(edgLabelMap,Img)
    adj=[-1,0;0,-1;1,0;0,1];
    [row,col]=size(edgLabelMap);
    resMap=edgLabelMap;
    
    [outerEdgMap,hasEdg]=findOuterEdg(edgLabelMap);
    while(hasEdg)
        for i=1:row
            for j=1:col
                if outerEdgMap(i,j)==1
                   mindis=Inf;
                   mink=1;
                   for k=1: length(adj)
                        xx=i+adj(k,1);
                        yy=j+adj(k,2);
                        if xx>=1 && xx<=row && yy>=1 && yy <=col && edgLabelMap(xx,yy)~=Inf
                            d=abs(Img(xx,yy)-Img(i,j));
                            if d<mindis
                               mindis=d;
                               mink=k;
                            end
                        end
                   end
                   resMap(i,j)=edgLabelMap(i+adj(mink,1),j+adj(mink,2));
                end
            end
        end  
        edgLabelMap=resMap;
        [outerEdgMap,hasEdg]=findOuterEdg(edgLabelMap);  
    end 
end

function [outerEdgMap,hasEdg]=findOuterEdg(edgLabelMap)
    adj=[-1,0;0,-1;1,0;0,1;1,1;1,-1;-1,1;-1,-1];
    [row,col]=size(edgLabelMap);
    outerEdgMap=zeros(row,col);
    
    hasEdg=0;
    for i=1:row
        for j=1:col
            if edgLabelMap(i,j)==Inf
                for k=1:length(adj)
                   xx=i+adj(k,1);
                   yy=j+adj(k,2);
                   if xx>=1 && xx<=row && yy>=1 && yy <=col && edgLabelMap(xx,yy)~=Inf
                      outerEdgMap(i,j)=1; 
                      if hasEdg==0
                         hasEdg=1; 
                         break;
                      end
                   end
                end
            end
        end
    end

end
    
function [mostLabel]=findMostSimilarLabel(cut,C,i,j)
    adj=[-1,0;0,-1;1,0;0,1];
    [row,col]=size(cut);
    labelVec=zeros(C,1);     
    roadMap=zeros(row,col);
    
    Xarray=zeros(row*col,1); 
    Yarray=zeros(row*col,1);
    AftPoint=1;                  
    BefPoint=1;                  
    Xarray(AftPoint)=i;
    Yarray(AftPoint)=j;

    while BefPoint <= AftPoint
        x=Xarray(BefPoint);
        y=Yarray(BefPoint);
        BefPoint=BefPoint+1;
        roadMap(x,y)=1;
        labelVec(cut(x,y),1)=labelVec(cut(x,y),1)+1;        
        for k=1:length(adj)
           xx=x+adj(k,1);
           yy=y+adj(k,2);
           if xx>=1 && xx<=row && yy>=1 && yy <=col && cut(xx,yy)~=Inf
               if (roadMap(xx,yy)==0)  
                   AftPoint=AftPoint+1;
                   Xarray(AftPoint)=xx;
                   Yarray(AftPoint)=yy;
                   roadMap(xx,yy)=2;                                   
               end
           end
        end
    end
    [maxnum,mostLabel]=max(labelVec);
end


function [corrMap]= getCorrMap(edgCLabelMap,CLabelMap, W)

    [rLabelMap,LabelNum]= getReLabelMap(edgCLabelMap);
    [row,col]=size(rLabelMap);
    LabelNumVec=zeros(LabelNum,1);
    for i=1:row
        for j=1:col
            label=rLabelMap(i,j);
            if label~=0 && label~=Inf
                LabelNumVec(label)=LabelNumVec(label)+1;
            end
        end
    end
    LabelNumMap=zeros(row,col);
    for i=1:row 
        for j=1:col
            tempLabel=rLabelMap(i,j);
            if tempLabel~=0 && tempLabel~=Inf
                LabelNumMap(i,j)=LabelNumVec(tempLabel);
            end
        end
    end
    
    kW = round(W*2/3);
    if mod(kW,2)==0
       kW=kW+1; 
    end
    kernel = ones(kW,kW);
    kernel((kW+1)/2,(kW+1)/2)= -(kW*kW-1);
    resMap = imfilter(CLabelMap,kernel, 'replicate');
    corrMap=(resMap~=0)& LabelNumMap<(col*row/10);

end

    
function [reEdgLabelMap,labelCnt]= getReLabelMap(edgCLabelMap)

    adj=[-1,0;0,-1;1,0;0,1];    
    [row,col]=size(edgCLabelMap);
    globalRoadMap=zeros(row,col);
    reEdgLabelMap=zeros(row,col);
    labelCnt=0;    
    Xarray=zeros(row*col,1); 
    Yarray=zeros(row*col,1);
    %relabel 
    for i=1:row 
        for j=1:col
            if  edgCLabelMap(i,j)~=Inf && globalRoadMap(i,j)==0 
                globalRoadMap(i,j)=1;
                labelCnt=labelCnt+1;                
                AftPoint=1;                   
                BefPoint=1;                   
                Xarray(AftPoint)=i;
                Yarray(AftPoint)=j;                
                while BefPoint <= AftPoint
                    x=Xarray(BefPoint);
                    y=Yarray(BefPoint);
                    BefPoint=BefPoint+1;
                    
                    currentLabel=edgCLabelMap(x,y);
                    globalRoadMap(x,y)=1;
                    reEdgLabelMap(x,y)=labelCnt;  %relabel                    
                    for k=1:length(adj)
                       xx=x+adj(k,1);
                       yy=y+adj(k,2);
                       if xx>=1 && xx<=row && yy>=1 &&yy <=col && edgCLabelMap(xx,yy)~=Inf
                           if globalRoadMap(xx,yy)==0   && edgCLabelMap(xx,yy)==currentLabel 
                               AftPoint=AftPoint+1;
                               Xarray(AftPoint)=xx;
                               Yarray(AftPoint)=yy;
                               globalRoadMap(xx,yy)=2;                                   
                           end
                       end
                    end
                end
            end
            %
        end
    end
    
end    

function [resMap] = mKMeans(I,C)
    
    [Idx,cent]=kmeans(I(:),C);
    CMap=reshape(Idx,size(I));
    [~,ind]=sort(cent);
    [~,ind3]=sort(ind);
    [row,col]=size(I);
    
    resMap=zeros(row,col);
    for i=1:row
        for j=1:col
            nexLabel=ind3(CMap(i,j));
            resMap(i,j)=nexLabel;
        end
    end

end