function [match,TheF1,TheF2]=getMatch(TheModel,ARGSet,IsShown,Para)
Para.MaxMRFIterNum=1000;
num=size(ARGSet,2);
Para.ThreadNum=min(Para.ThreadNum,num);
nodeNum=TheModel.nodeNum;
match=zeros(nodeNum,num);
TheF1(num).matrix=[];
TheF2(num).matrix=[];
if(Para.ThreadNum>1)
    tmp=round(linspace(1,num+1,Para.ThreadNum+1));
    list1=tmp(1:end-1);
    list2=tmp(2:end)-1;
    Division(Para.ThreadNum).subARGSet=[];
    Division(Para.ThreadNum).match=[];
    Division(Para.ThreadNum).TheF1=[];
    Division(Para.ThreadNum).TheF2=[];
    for i=1:Para.ThreadNum
        len=list2(i)-list1(i)+1;
        Division(i).subARGSet=ARGSet(list1(i):list2(i));
        Division(i).match=zeros(nodeNum,len);
        Division(i).TheF1(len).matrix=[];
        Division(i).TheF2(len).matrix=[];
    end
    clear ARGSet
    parfor par=1:Para.ThreadNum
        len=size(Division(par).subARGSet,2);
        for i=1:len
            if(Division(par).subARGSet(i).arg.nodeNum>100)
                error('An ARG is too large, i.e. its node number is larger than 100.');
            end
            [x,f1,f2]=MRF(TheModel,Division(par).subARGSet(i).arg,IsShown,Para);
            Division(par).match(:,i)=x';
            Division(par).TheF1(i).matrix=f1;
            Division(par).TheF2(i).matrix=f2;
            Division(par).subARGSet(i).arg=[];
        end
    end
    for i=1:Para.ThreadNum
        match(:,list1(i):list2(i))=Division(i).match;
        c=0;
        for j=list1(i):list2(i)
            c=c+1;
            TheF1(j).matrix=Division(i).TheF1(c).matrix;
            TheF2(j).matrix=Division(i).TheF2(c).matrix;
        end
    end 
else
    for argNo=1:num
        [x,f1,f2]=MRF(TheModel,ARGSet(argNo).arg,IsShown,Para);
        match(:,argNo)=x';
        TheF1(argNo).matrix=f1;
        TheF2(argNo).matrix=f2;
    end
end


function [x,f1,f2]=MRF(TheModel,TheARG,IsShown,Para)
unaryNum=size(TheModel.unaryF,2);
pairwiseNum=size(TheModel.pairwiseF,2);
UnaryPenalty=zeros(TheModel.nodeNum,TheARG.nodeNum);
for atri=1:unaryNum
    for j=1:TheModel.nodeNum
        if(TheARG.nodeNum>0)
            UnaryPenalty(j,:)=UnaryPenalty(j,:)+sum((repmat(TheModel.unaryF(atri).f(:,j),[1,TheARG.nodeNum])-TheARG.unaryF(atri).f).^2,1).*TheModel.unaryF(atri).w;
        end
    end
end
for i=1:TheModel.nodeNum
    for j=1:TheModel.nodeNum
        if(j==i)
            continue;
        end
        PairwisePenalty=zeros(TheARG.nodeNum,TheARG.nodeNum);
        if(TheARG.nodeNum>1)
            for atri=1:pairwiseNum
                PairwisePenalty=PairwisePenalty+reshape(sum((repmat(TheModel.pairwiseF(atri).f(:,i,j),[1,TheARG.nodeNum,TheARG.nodeNum])-TheARG.pairwiseF(atri).f).^2,1),size(PairwisePenalty)).*TheModel.pairwiseF(atri).w;
            end
        end
        MRF(i,j).Transfer=PairwisePenalty;
    end
end
V=0:TheModel.nodeNum-1;
Enum=TheModel.nodeNum*(TheModel.nodeNum-1)/2;
E=zeros(2,Enum);
t=0;
for i=0:TheModel.nodeNum-2
    len=TheModel.nodeNum-1-i;
    E(:,t+1:t+len)=[ones(1,len).*i;i+1:i+len];
    t=t+len;
end
f1=[UnaryPenalty';TheModel.penalty.unmatch.*ones(1,TheModel.nodeNum)];
f2=ones(TheARG.nodeNum+1,TheARG.nodeNum+1,Enum).*TheModel.penalty.unmatchPair;
f2_normalized=ones(TheARG.nodeNum+1,TheARG.nodeNum+1,Enum).*TheModel.penalty.unmatchPair;
t=0;
TheDiag=zeros(TheARG.nodeNum,TheARG.nodeNum);
TheDiag(linspace(1,TheARG.nodeNum^2,TheARG.nodeNum))=TheModel.penalty.large;
for i=1:TheModel.nodeNum-1
    for j=i+1:TheModel.nodeNum
        t=t+1;
        tmp=(TheModel.link(i).to(j)/sum(TheModel.link(i).to)+TheModel.link(j).to(i)/sum(TheModel.link(j).to))/2;
        f2(1:TheARG.nodeNum,1:TheARG.nodeNum,t)=MRF(i,j).Transfer;
        f2_normalized(:,:,t)=f2(:,:,t).*tmp;
        ManyToOneIndex=linspace(1,(TheARG.nodeNum+1)^2,TheARG.nodeNum+1);
        ManyToOneIndex=ManyToOneIndex(1:TheARG.nodeNum)+((TheARG.nodeNum+1)^2)*(t-1);
        f2(ManyToOneIndex)=TheModel.penalty.large;
        f2_normalized(ManyToOneIndex)=TheModel.penalty.large;
    end
end
f1_normalized=f1-repmat(min(f1,[],1),[TheARG.nodeNum+1,1]);
x=TRWSProcess(V,E,f1_normalized,f2_normalized,Para)';
TheX=x;TheX(TheX>TheARG.nodeNum)=-1;disp(TheX);
