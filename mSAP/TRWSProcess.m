function x=TRWSProcess(V,E,f1,f2,Para)
%eps=0.00000000000000001;
eps=0.000000000001;
if(size(E,2)==1)
    len=size(f1,1);
    tmp=repmat(f1(:,1),[1,len])+repmat(f1(:,2),[1,len])'+f2;
    [v,r1]=min(tmp,[],1);
    [~,r2]=min(v,[],2);
    x=[r1(r2),r2]-1;
else
    x=TRWS_mex(E,f1,f2,[eps,Para.MaxMRFIterNum],0)';
end
x=double(x')+1;
