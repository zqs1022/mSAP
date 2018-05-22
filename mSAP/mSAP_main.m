function mSAP_main(model,ARGSet,ModelName,Para,TargetDir)
Para.MaxMRFIterNum=1000;
tic;
record.model=model;
record.Log.Insert=0;
record.Log.Delete=0;
filename=sprintf('%s%s_model.mat',TargetDir,ModelName);
save(filename,'record');
load(sprintf('%s%s_model.mat',TargetDir,model.name),'record');
model=record(end).model;
clear record
for Iter=1:Para.MaxIterationNum
    ARG_num=size(ARGSet,2);
    if(ARG_num<=1)
        delete(filename);
        error('bad initial labeling.');
    end
    model=ModelAttributeEstimation(model,ARGSet,Para,TargetDir);
    [model,IsChange1,argFeature,MaxUnmatchNum]=EliminateNode(model,ARGSet,Para,TargetDir);
    [model,IsChange2]=AddNode(argFeature,MaxUnmatchNum,model,ARGSet,Para,TargetDir);
    % if((IsChange1==0)&&(IsChange2==0))
    %     break;
    % end
end
load(sprintf('%s%s_model.mat',TargetDir,model.name),'record');
time=toc;
save(sprintf('%s%s_model.mat',TargetDir,model.name),'record','ARGSet','time');


function model=ModelAttributeEstimation(model,ARGSet,Para,TargetDir)
TheModel=model;
nodeNum=TheModel.nodeNum;
model.name=TheModel.name;
model.nodeNum=nodeNum;
unaryNum=size(model.unaryF,2);
pairwiseNum=size(model.pairwiseF,2);
for atri=1:unaryNum
    model.unaryF(atri).f=zeros(size(TheModel.unaryF(atri).f));
end
for atri=1:pairwiseNum
    model.pairwiseF(atri).f=zeros(size(TheModel.pairwiseF(atri).f));
end
UnaryCount=zeros(1,nodeNum);
PairwiseCount=zeros(1,nodeNum,nodeNum);
argNum=size(ARGSet,2);
IsShown=0;
match=getMatch(TheModel,ARGSet,IsShown,Para);
for argNo=1:argNum
    x=match(:,argNo)';
    argNodeNum=ARGSet(argNo).arg.nodeNum;
    list=find(x<=argNodeNum);
    TheX=x(list);
    if(size(list,2)>0)
        Attribute=getAttribute(TheX,TheX,ARGSet(argNo).arg);
        for atri=1:unaryNum
            model.unaryF(atri).f(:,list)=model.unaryF(atri).f(:,list)+Attribute.unaryF(atri).f;
        end
        for atri=1:pairwiseNum
            model.pairwiseF(atri).f(:,list,list)=model.pairwiseF(atri).f(:,list,list)+Attribute.pairwiseF(atri).f;
        end
        UnaryCount(:,list)=UnaryCount(:,list)+1;
        PairwiseCount(:,list,list)=PairwiseCount(:,list,list)+1;
    end
end
UnaryCount=max(UnaryCount,1);
PairwiseCount=max(PairwiseCount,1);
for atri=1:unaryNum
    model.unaryF(atri).f=model.unaryF(atri).f./repmat(UnaryCount,[size(model.unaryF(atri).f,1),1]);
end
for atri=1:pairwiseNum
    model.pairwiseF(atri).f=model.pairwiseF(atri).f./repmat(PairwiseCount,[size(model.pairwiseF(atri).f,1),1,1]);
end
load(sprintf('%s%s_model.mat',TargetDir,TheModel.name),'record');
Iteration=size(record,2)+1;
record(Iteration).model=model;
record(Iteration).Log=record(Iteration-1).Log;
save(sprintf('%s%s_model.mat',TargetDir,TheModel.name),'record');


function [TheModel,IsChange,argFeature,MaxUnmatchNum]=EliminateNode(TheModel,ARGSet,Para,TargetDir)
N=size(ARGSet,2);
UnaryPenalty=zeros(TheModel.nodeNum,1);
PairwisePenalty=zeros(TheModel.nodeNum,TheModel.nodeNum);
argFeature.TheX=zeros(TheModel.nodeNum,N);
argFeature.argNodeNum=zeros(N,1);
MaxUnmatchNum=0;
IsShown=0;
[match,TheF1,TheF2]=getMatch(TheModel,ARGSet,IsShown,Para);
for argNo=1:N
    x=match(:,argNo)';
    f1=TheF1(argNo).matrix;
    f2=TheF2(argNo).matrix;
    argNodeNum=ARGSet(argNo).arg.nodeNum;
    list=find(x<=argNodeNum);
    TheX=x(list);
    RawX=x;
    x(x>argNodeNum)=-1;
    argFeature.TheX(:,argNo)=x';
    argFeature.argNodeNum(argNo)=argNodeNum;
    MaxUnmatchNum=max(MaxUnmatchNum,argNodeNum-size(TheX,2));
    tmp=(0:TheModel.nodeNum-1).*(argNodeNum+1)+RawX;
    UnaryPenalty=UnaryPenalty+f1(tmp)';
    t=0;
    for i=1:TheModel.nodeNum-1
        for j=i+1:TheModel.nodeNum
            t=t+1;
            tmp=f2(RawX(i),RawX(j),t);
            PairwisePenalty(i,j)=PairwisePenalty(i,j)+tmp;
            PairwisePenalty(j,i)=PairwisePenalty(j,i)+tmp;
        end
    end
end
UnaryPenalty=UnaryPenalty./N;
PairwisePenalty=PairwisePenalty./N;
TotalPenalty=UnaryPenalty+sum(PairwisePenalty,2)./(TheModel.nodeNum-1);
[Value,NonMatchNode]=max(TotalPenalty);
disp([TotalPenalty,UnaryPenalty,TotalPenalty-UnaryPenalty]);
%pause
load(sprintf('%s%s_model.mat',TargetDir,TheModel.name),'record');
if((Value>Para.PenaltyThreshold)&&(TheModel.nodeNum>2))
    IsChange=1;
    Remain=setdiff((1:TheModel.nodeNum)',NonMatchNode);
    argFeature.TheX=argFeature.TheX(Remain,:);
    MaxUnmatchNum=MaxUnmatchNum+1;
    TheModel=NodeEliminate(NonMatchNode,TheModel);
    record(end).model=TheModel;
    record(end).Log.Delete=record(end).Log.Delete+1;
else
    IsChange=0;
    record(end).model=TheModel;
end
save(sprintf('%s%s_model.mat',TargetDir,TheModel.name),'record');


function [TheModel,IsChange]=AddNode(argFeature,MaxUnmatchNum,TheModel,ARGSet,Para,TargetDir)
unaryNum=size(TheModel.unaryF,2);
pairwiseNum=size(TheModel.pairwiseF,2);
N=size(ARGSet,2);
argFeature.ModelNodeNum=TheModel.nodeNum;
argFeature.LabelNum=MaxUnmatchNum;
for atri=1:unaryNum
    argFeature.unaryF(atri).f=zeros(size(TheModel.unaryF(atri).f,1),MaxUnmatchNum,N);
end
for atri=1:pairwiseNum
    argFeature.pairwiseF(atri).f=zeros(size(TheModel.pairwiseF(atri).f,1),TheModel.nodeNum,MaxUnmatchNum,N);
end
argFeature.TargetNodeNum=zeros(N,1);
for argNo=1:N
    list=find(argFeature.TheX(:,argNo)~=-1);
    if(size(list,1)>0)
        list2=argFeature.TheX(list,argNo)';
        list1=setdiff(1:argFeature.argNodeNum(argNo),list2);
        Attribute=getAttribute(list1,list2,ARGSet(argNo).arg);
        targetNum=size(list1,2);
        argFeature.TargetNodeNum(argNo)=targetNum;
        for atri=1:unaryNum
            argFeature.unaryF(atri).f(:,1:targetNum,argNo)=Attribute.unaryF(atri).f;
        end
        for atri=1:pairwiseNum
            argFeature.pairwiseF(atri).f(:,list,1:targetNum,argNo)=Attribute.pairwiseF(atri).f;
        end
    end
end
%get the attributes and penalty of the new node (the average attribute of the predicted matched nodes in N images)
[~,UnaryPenalty,NodeRelatedPairPenalty,Add]=MRF_ForNewNode(argFeature,TheModel,Para,MaxUnmatchNum,[]);
PairwisePenalty=mean(NodeRelatedPairPenalty);
Penalty=UnaryPenalty+PairwisePenalty;
disp([Penalty,UnaryPenalty,PairwisePenalty]);
%pause;
%add the new node
IsChange=0;
if(Penalty<=Para.PenaltyThreshold)
    IsChange=1;
    load(sprintf('%s%s_model.mat',TargetDir,TheModel.name),'record');
    TheModel.nodeNum=TheModel.nodeNum+1;
    for atri=1:unaryNum
        TheModel.unaryF(atri).f=[TheModel.unaryF(atri).f,Add.unaryF(atri).f];
    end
    for atri=1:pairwiseNum
        tmp=size(TheModel.pairwiseF(atri).f,1);
        TheModel.pairwiseF(atri).f(:,TheModel.nodeNum,1:TheModel.nodeNum-1)=reshape(Add.pairwiseF(atri).f,[tmp,TheModel.nodeNum-1]);
        TheModel.pairwiseF(atri).f(:,1:TheModel.nodeNum-1,TheModel.nodeNum)=reshape(Add.pairwiseF(atri).f,[tmp,TheModel.nodeNum-1]);
        TheModel.pairwiseF(atri).f(:,TheModel.nodeNum,TheModel.nodeNum)=0;
    end
    for i=1:TheModel.nodeNum-1
        TheModel.link(i).to=[TheModel.link(i).to;1];
    end
    TheModel.link(TheModel.nodeNum).to=ones(TheModel.nodeNum,1);
    TheModel.link(TheModel.nodeNum).to(TheModel.nodeNum)=0;
    record(end).model=TheModel;
    record(end).Log.Insert=record(end).Log.Insert+1;
    save(sprintf('%s%s_model.mat',TargetDir,TheModel.name),'record');
    %pause;
end


function [x,UnaryPenalty,NodeRelatedPairPenalty,Add]=MRF_ForNewNode(argFeature,TheModel,Para,MaxUnmatchNum,TheLink)
unaryNum=size(argFeature.unaryF,2);
pairwiseNum=size(argFeature.pairwiseF,2);
N=size(argFeature.unaryF(1).f,3);
V=0:N-1;
Enum=N*(N-1)/2;
E=zeros(2,Enum);
t=0;
for i=0:N-2
    len=N-1-i;
    E(:,t+1:t+len)=[ones(1,len).*i;i+1:i+len];
    t=t+len;
end
LabelNum=argFeature.LabelNum;
f1=ones(LabelNum,N);
f2=zeros(LabelNum,LabelNum,Enum);
t=0;
TheDiag=zeros(LabelNum,LabelNum);
TheDiag(linspace(1,LabelNum^2,LabelNum))=TheModel.penalty.large;
for i=1:N-1
    for j=i+1:N
        UnaryPenalty=zeros(LabelNum,LabelNum);
        for atri=1:unaryNum
            tmpSize=size(argFeature.unaryF(atri).f,1);
            tmpDiff=repmat(reshape(argFeature.unaryF(atri).f(:,:,i),[tmpSize,LabelNum,1]),[1,1,LabelNum])-repmat(reshape(argFeature.unaryF(atri).f(:,:,j),[tmpSize,1,LabelNum]),[1,LabelNum,1]);
            UnaryPenalty=UnaryPenalty+reshape(sum(tmpDiff.^2,1),size(UnaryPenalty)).*(TheModel.unaryF(atri).w/(2*N^2));
        end
        ModelNodeSet1=find(argFeature.TheX(:,i)'~=-1);
        ModelNodeSet2=find(argFeature.TheX(:,j)'~=-1);
        ValidModelNodeSet=intersect(ModelNodeSet1,ModelNodeSet2);
        PairwisePenalty=zeros(LabelNum,LabelNum,TheModel.nodeNum);
        for k=ValidModelNodeSet
            CountValidARG=sum(argFeature.TheX(k,:)'~=-1);
            term_NonNone=zeros(LabelNum,LabelNum);
            for atri=1:pairwiseNum
                tmpSize=size(argFeature.pairwiseF(atri).f,1);
                tmpDiff=repmat(reshape(argFeature.pairwiseF(atri).f(:,k,:,i),[tmpSize,LabelNum,1]),[1,1,LabelNum])-repmat(reshape(argFeature.pairwiseF(atri).f(:,k,:,j),[tmpSize,1,LabelNum]),[1,LabelNum,1]);
                term_NonNone=term_NonNone+reshape(sum(tmpDiff.^2,1),[LabelNum,LabelNum]).*(TheModel.pairwiseF(atri).w/(2*N*CountValidARG));
            end
            PairwisePenalty(:,:,k)=term_NonNone;
        end
        for k=setdiff(1:TheModel.nodeNum,ValidModelNodeSet)
            CountValidARG=sum(argFeature.TheX(k,:)'~=-1);
            term_None=TheModel.penalty.unmatchPair/(N*(N+CountValidARG));
            PairwisePenalty(:,:,k)=term_None;
        end
        Transfer=UnaryPenalty+sum(PairwisePenalty,3)./TheModel.nodeNum;
        Transfer(argFeature.TargetNodeNum(i)+1:LabelNum,:)=TheModel.penalty.large;
        Transfer(:,argFeature.TargetNodeNum(j)+1:LabelNum)=TheModel.penalty.large;
        t=t+1;
        f2(:,:,t)=max(Transfer,TheDiag);
    end
end
x=zeros(N,1);
tmp=TRWSProcess(V,E,f1,f2,Para);
x(:,1)=tmp';
x=double(x');
[UnaryPenalty,NodeRelatedPairPenalty,Add]=getInfo_NewNode(x,MaxUnmatchNum,argFeature,TheModel);


function [UnaryPenalty,NodeRelatedPairPenalty,Add]=getInfo_NewNode(x,MaxUnmatchNum,argFeature,TheModel)
unaryNum=size(argFeature.unaryF,2);
pairwiseNum=size(argFeature.pairwiseF,2);
N=size(argFeature.unaryF(1).f,3);
UnaryPenalty=0;
tmp=x+(0:N-1).*MaxUnmatchNum;
for atri=1:unaryNum
    theSize=size(argFeature.unaryF(atri).f);
    tmpF=reshape(argFeature.unaryF(atri).f,[theSize(1),theSize(2)*theSize(3)]);
    Add.unaryF(atri).f=mean(tmpF(:,tmp),2);
    tmpDiff=tmpF(:,tmp)-repmat(Add.unaryF(atri).f,[1,N]);
    UnaryPenalty=UnaryPenalty+sum(sum(tmpDiff.^2,1),2)*TheModel.unaryF(atri).w/N;
end
for i=1:pairwiseNum
    Add.pairwiseF(atri).f=zeros(size(TheModel.pairwiseF(atri).f,1),TheModel.nodeNum);
end
NodeRelatedPairPenalty=zeros(TheModel.nodeNum,1);
for i=1:TheModel.nodeNum
    ValidModelNodeSet=find(argFeature.TheX(i,:)~=-1);
    ValidSize=size(ValidModelNodeSet,2);
    tmp=i+TheModel.nodeNum.*(x(ValidModelNodeSet)-1)+(TheModel.nodeNum*MaxUnmatchNum).*(ValidModelNodeSet-1);
    for atri=1:pairwiseNum
        theSize=size(argFeature.pairwiseF(atri).f);
        tmpF=reshape(argFeature.pairwiseF(atri).f,[theSize(1),theSize(2)*theSize(3)*theSize(4)]);
        Add.pairwiseF(atri).f(:,i)=mean(tmpF(:,tmp),2);
        tmpDiff=tmpF(:,tmp)-repmat(Add.pairwiseF(atri).f(:,i),[1,ValidSize]);
        NodeRelatedPairPenalty(i)=NodeRelatedPairPenalty(i)+sum(sum(tmpDiff.^2,1),2)*TheModel.pairwiseF(atri).w/N;
    end
end


function TheModel=NodeEliminate(Nodes,TheModel)
unaryNum=size(TheModel.unaryF,2);
pairwiseNum=size(TheModel.pairwiseF,2);
n=size(Nodes,1);
Remain=setdiff((1:TheModel.nodeNum)',Nodes);
TheModel.nodeNum=TheModel.nodeNum-n;
for atri=1:unaryNum
    TheModel.unaryF(atri).f=TheModel.unaryF(atri).f(:,Remain);
end
for atri=1:pairwiseNum
    TheModel.pairwiseF(atri).f=TheModel.pairwiseF(atri).f(:,Remain,Remain);
end
TheModel.link=TheModel.link(Remain);
for i=1:TheModel.nodeNum
    TheModel.link(i).to=TheModel.link(i).to(Remain);
end


function Attribute=getAttribute(xList1,xList2,TheARG)
unaryNum=size(TheARG.unaryF,2);
pairwiseNum=size(TheARG.pairwiseF,2);
for atri=1:unaryNum
    Attribute.unaryF(atri).f=TheARG.unaryF(atri).f(:,xList1);
end
for atri=1:pairwiseNum
    Attribute.pairwiseF(atri).f=TheARG.pairwiseF(atri).f(:,xList2,xList1);
end


function GoodARGSet=selectGoodARGs(ARGSet,model,Para)
IsShown=0;
[match,TheF1,TheF2]=getMatch(model,ARGSet,IsShown,Para);
list=[];
energy=[];
for argNo=1:size(match,2)
    argNodeNum=ARGSet(argNo).arg.nodeNum;
    if(sum(match(:,argNo)<=argNodeNum)>=max(Para.GoodMatchRate*model.nodeNum,2))
        list=[list,argNo];
        x=match(:,argNo)';
        tmp=(0:model.nodeNum-1).*(argNodeNum+1)+x;
        UnaryEnergy=sum(TheF1(argNo).matrix(tmp));
        PairwisePenalty=zeros(model.nodeNum,model.nodeNum);
        t=0;
        for i=1:model.nodeNum-1
            for j=i+1:model.nodeNum
                t=t+1;
                tmp=TheF2(argNo).matrix(x(i),x(j),t);
                PairwisePenalty(i,j)=tmp;
                PairwisePenalty(j,i)=tmp;
            end
        end
        TotalPairwisePenalty=sum(sum(PairwisePenalty,2),1)./((model.nodeNum-1)*2);
        energy=[energy,UnaryEnergy+TotalPairwisePenalty];
    end
end
[~,tmp]=sort(energy);
GoodARGSet=ARGSet(list(tmp));
MaxARGNum=min(size(energy,2),Para.MaxARGNumPerIteration);
GoodARGSet=GoodARGSet(1:MaxARGNum);
