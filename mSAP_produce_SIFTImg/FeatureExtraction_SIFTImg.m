function FeatureExtraction_SIFTImg(ImageSetName)
addpath('./siftDemoV4');
FeatureSet_positive=getFeatures(ImageSetName,'');
filename=sprintf('./mat/feature/%s_feature.mat',ImageSetName);
save(filename,'FeatureSet_positive');


function FeatureSet=getFeatures(ImageSetName,PreName)
Small=0.00000001;
root=sprintf('./mat/ImgSet/%s/%s',ImageSetName,PreName);
TotalImgNo=0;
filter=fspecial('gaussian',[15,15],3);
while(true)
    filename=sprintf('%s%03d.jpg',root,TotalImgNo+1);
    if(exist(filename,'file')~=2)
        break;
    end
    TotalImgNo=TotalImgNo+1;
end
FeatureSet(TotalImgNo).rawdata=[];
FeatureSet(TotalImgNo).arg=[];
for ImgNo=1:TotalImgNo
    filename=sprintf('%s%03d.jpg',root,ImgNo);
    rawdata.image=imread(filename);
    smoothed_img=imfilter(rawdata.image,filter);
    [~,rawdata.descriptors,rawdata.locs]=sift(smoothed_img);
    FeatureSet(ImgNo).rawdata=rawdata;
    %showkeys(rawdata.image,rawdata.locs);
    showkeys(smoothed_img,rawdata.locs);
    arg.name=sprintf('%03d%s%s',ImgNo,ImageSetName,PreName);
    arg.nodeNum=size(rawdata.descriptors,1);
    arg.unaryF(1).f=rawdata.descriptors';
    arg.unaryF(2).f=rawdata.locs(:,4)';
    tmp_1=repmat(rawdata.locs(:,1),[1,arg.nodeNum])-repmat(rawdata.locs(:,1)',[arg.nodeNum,1]);
    tmp_2=repmat(rawdata.locs(:,2),[1,arg.nodeNum])-repmat(rawdata.locs(:,2)',[arg.nodeNum,1]);
    dist=sqrt(tmp_1.^2+tmp_2.^2);
    centerline_ori=acos(tmp_2./max(dist,Small));
    centerline_ori(tmp_1>0)=-centerline_ori(tmp_1>0);
    scale=rawdata.locs(:,3);
    ori=rawdata.locs(:,4);
    ori_matrix=pi-abs(mod(repmat(ori,[1,arg.nodeNum])-repmat(ori',[arg.nodeNum,1]),2*pi)-pi);
    arg.pairwiseF(1).f=reshape(ori_matrix,[1,arg.nodeNum,arg.nodeNum]);
    ori_matrix=pi-abs(mod(repmat(ori,[1,arg.nodeNum])-centerline_ori,2*pi)-pi);
    arg.pairwiseF(2).f=reshape(ori_matrix,[1,arg.nodeNum,arg.nodeNum]);
    ori_matrix=pi-abs(mod(repmat(ori',[arg.nodeNum,1])-centerline_ori,2*pi)-pi);
    arg.pairwiseF(3).f=reshape(ori_matrix,[1,arg.nodeNum,arg.nodeNum]);
    arg.pairwiseF(4).f=log(reshape(sqrt(repmat(scale,[1,arg.nodeNum]).^2+repmat(scale',[arg.nodeNum,1]).^2)./max(dist,Small),[1,arg.nodeNum,arg.nodeNum]));
    arg.pairwiseF(5).f=log(reshape(repmat(scale,[1,arg.nodeNum])./max(repmat(scale',[arg.nodeNum,1]),Small),[1,arg.nodeNum,arg.nodeNum]));
    arg.pairwiseF(6).f=reshape(centerline_ori,[1,arg.nodeNum,arg.nodeNum]);
    FeatureSet(ImgNo).arg=arg;
    fprintf('%d / %d\n',ImgNo,TotalImgNo);
end
