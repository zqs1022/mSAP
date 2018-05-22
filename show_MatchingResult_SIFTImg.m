function show_MatchingResult_SIFTImg(CldName,x,Para,TheLink)
ID_frame=sscanf(CldName(1:3),'%d');
ImageSetName=CldName(4:end);
load(sprintf('./mSAP_produce_SIFTImg/mat/feature/%s_feature.mat',ImageSetName),'FeatureSet_positive');
rawdata=FeatureSet_positive(ID_frame).rawdata;
n=size(FeatureSet_positive(ID_frame).rawdata.locs,1);
TheX=x(x<=n);
rawdata.locs=zeros(size(x,2),4);
rawdata.locs(x<=n,:)=FeatureSet_positive(ID_frame).rawdata.locs(TheX,:);
figure
showkeys(rawdata.image,FeatureSet_positive(ID_frame).rawdata.locs);
hold on;
for i=1:size(FeatureSet_positive(ID_frame).rawdata.locs,1)
    len=FeatureSet_positive(ID_frame).rawdata.locs(i,3);
    PlotEllipse(FeatureSet_positive(ID_frame).rawdata.locs(i,1),FeatureSet_positive(ID_frame).rawdata.locs(i,2),len,len/2,FeatureSet_positive(ID_frame).rawdata.locs(i,4),'c');
end
for i=1:size(TheX,2)
    len=FeatureSet_positive(ID_frame).rawdata.locs(TheX(i),3);
    PlotEllipse(FeatureSet_positive(ID_frame).rawdata.locs(TheX(i),1),FeatureSet_positive(ID_frame).rawdata.locs(TheX(i),2),len,len/2,FeatureSet_positive(ID_frame).rawdata.locs(TheX(i),4),'m');
end
if(true)
    num=size(TheLink,2);
    for i=1:num
        if(x(i)>n)
            continue;
        end
        for j=find(TheLink(i).to==1)'
            if(x(j)>n)
                continue;
            end
            if(TheLink(j).to(i)==1)
                line([rawdata.locs(i,2),rawdata.locs(j,2)],[rawdata.locs(i,1),rawdata.locs(j,1)],'Color','m','LineWidth',1);
            else
                line([rawdata.locs(i,2),rawdata.locs(j,2)],[rawdata.locs(i,1),rawdata.locs(j,1)],'Color','w','LineWidth',1);
            end
        end
    end
end
axis off
hold off
%saveas(gcf,sprintf('%s.fig',CldName),'fig'); 


function PlotEllipse(OY,OX,Long,Short,Angle,col)
tmp=linspace(0,2*pi,100);
tmpY=sin(tmp).*Short.*2;
tmpX=cos(tmp).*Long.*2;
template=[cos(-Angle),-sin(-Angle);sin(-Angle),cos(-Angle)];
A=template*[tmpX;tmpY];
plot(OX+A(1,:),OY+A(2,:),col,'LineWidth',2);
