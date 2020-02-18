function [grayMap,rgbMap,Acc] = calAcc(resMap,groundT)
%acc
    [row,col]=size(resMap);
    [valList,~]=sort(unique(groundT(:)));

    grayMap=zeros(row,col,'uint8');    
    for i=1:row
        for j=1:col
            grayMap(i,j)=valList(resMap(i,j));
        end
    end
    
    color=[ 61, 38,168; %deep blue
            39,150,235; %light blue
           128,203, 88; %green            
           249,250, 20; %yellow
           208,191, 39; %brown
           205,  0,  0; %red
        ];
    rgbMap=zeros(row,col,3,'uint8');
    for i=1:row
        for j=1:col
            cur=resMap(i,j);
            rgbMap(i,j,1)=color(cur,1);
            rgbMap(i,j,2)=color(cur,2);
            rgbMap(i,j,3)=color(cur,3);
        end
    end
    Acc=100*sum(sum(grayMap==groundT))/(row*col); 
     
end

