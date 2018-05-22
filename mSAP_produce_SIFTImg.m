function [model,ARGSet]=mSAP_produce_SIFTImg(ImageSetName)
addpath('./mSAP_produce_SIFTImg/siftDemoV4');
load(sprintf('./mSAP_produce_SIFTImg/mat/feature/%s_feature.mat',ImageSetName),'FeatureSet_positive');
ARG_num=size(FeatureSet_positive,2);
ARGSet(ARG_num).arg=[];
for i=1:ARG_num
    ARGSet(i).arg=FeatureSet_positive(i).arg;
end
%model=model_produce(FeatureSet_positive,ImageSetName);
load(sprintf('./mSAP_produce_SIFTImg/mat/feature/%s_model.mat',ImageSetName),'model');


function model=model_produce(FeatureSet,ImageSetName)
ARG_num=size(FeatureSet,2);
model=[];
figure;
for i=1:ARG_num
    %FeatureSet(i).rawdata.locs=FeatureSet(i).rawdata.locs(1,:);
    showkeys(rgb2gray(FeatureSet(i).rawdata.image),FeatureSet(i).rawdata.locs);
    locs=FeatureSet(i).rawdata.locs;
    fNum=size(locs,1);
    fprintf('Frame %d fNum %d\n',i,fNum);
    A=input('Input ''y'' to choose the model','s');
    if((strcmp(A,'y')==1)||(strcmp(A,'Y')==1))
        model=ModelInteraction(FeatureSet(i));
        save(sprintf('./mSAP_produce_SIFTImg/mat/feature/%s_model.mat',ImageSetName),'model');
        break;
    end
end


function model=ModelInteraction(FeatureSet)
list=[];
[h,w,~]=size(FeatureSet.rawdata.image);
OriNorm=2*pi*2;
hold on
while(1)
    waitforbuttonpress;
    point1=get(gca,'CurrentPoint');
    rbbox;
    point2=get(gca,'CurrentPoint');
    point1=point1(1,1:2);
    point2=point2(1,1:2);
    MinH=min(point1(2),point2(2));
    MaxH=max(point1(2),point2(2));
    MinW=min(point1(1),point2(1));
    MaxW=max(point1(1),point2(1));
    if((MinH>h)||(MaxH<0)||(MinW>w)||(MaxW<0))
        model.nodeNum=size(list,2);
        model.name=FeatureSet.arg.name;
        model.list=list;
        for atri=1:size(FeatureSet.arg.unaryF,2)
            model.unaryF(atri).f=FeatureSet.arg.unaryF(atri).f(:,list);
            model.unaryF(atri).w=1;
        end
        for atri=1:size(FeatureSet.arg.pairwiseF,2)
            model.pairwiseF(atri).f=FeatureSet.arg.pairwiseF(atri).f(:,list,list);
            model.pairwiseF(atri).w=1;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        model.penalty.unmatch=100;
        model.penalty.unmatchPair=100;
        model.penalty.large=100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:model.nodeNum
            model.link(i).to=ones(model.nodeNum,1);
            model.link(i).to(i)=0;
        end
        hold off
        break;
    end
    tmp=point2-point1;
    if(tmp(2)<0)
        ori=acos(tmp(1)/sqrt(sum(tmp.^2)));
    else
        ori=-acos(tmp(1)/sqrt(sum(tmp.^2)));
    end
    tmp_f=[point1(2)/h,point1(1)/w,ori/OriNorm];
    locs=FeatureSet.rawdata.locs(:,[1,2,4]);
    locs_norm=locs;
    locs_norm(:,1)=locs_norm(:,1)./h;
    locs_norm(:,2)=locs_norm(:,2)./w;
    locs_norm(:,3)=locs_norm(:,3)./OriNorm;
    [~,index]=sort(sum((locs_norm-repmat(tmp_f,[size(locs_norm,1),1])).^2,2));
    index=index(1);
    list=[list,index];
    len=FeatureSet.rawdata.locs(index,3)*2;
    plot([locs(index,2),locs(index,2)+len*cos(locs(index,3))],[locs(index,1),locs(index,1)-len*sin(locs(index,3))],'r-','LineWidth',2);
    fprintf('nodeNum: %d\n',size(list,2));
end
hold off
