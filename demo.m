% Developed by Quanshi Zhang in June, 2014.

% License & disclaimer

% Copyright 2014 Quanshi Zhang (zqs1022@gmail.com)
% This software can be used for research purposes only. 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

% Related papers

% 1. Quanshi Zhang, Xuan Song, Xiaowei Shao, Huijing Zhao, and Ryosuke Shibasaki,
% "Attributed Graph Mining and Matching: An Attempt to Define and Extract Soft
% Attributed Patterns", in Proc. of IEEE International Conference on Computer
% Vision and Pattern Recognition (CVPR) 2014.

% 2. Quanshi Zhang, Xuan Song, Xiaowei Shao, Huijing Zhao, and Ryosuke Shibasaki,
% "Object Discovery: Soft Attributed Graph Mining", to appear in IEEE Transactions
% on Pattern Analysis and Machine Intelligence (TPAMI).


function demo

ImageSetName='panda';

compile; % compile the code

startup;

% load the ARGs and the initial graph template
[model,ARGSet]=mSAP_produce_SIFTImg(ImageSetName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Para.PenaltyThreshold=2.3; % threshold tau mentioned in the paper
Para.MaxIterationNum=10; % the maximum iteration number M mentioned in the paper
try
    p=gcp;
    Para.ThreadNum=p.NumWorkers; % the number of workers for parallel computing
catch
    Para.ThreadNum=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings of attribute weights, P_none, and Q_none
model.unaryF(1).w=3;
model.unaryF(2).w=1;
model.pairwiseF(1).w=0.5;
model.pairwiseF(2).w=0.5;
model.pairwiseF(3).w=1.0;
model.pairwiseF(4).w=1.0;
model.pairwiseF(5).w=0.5;
model.pairwiseF(6).w=0.5;
model.penalty.unmatch=3;
model.penalty.unmatchPair=5;
model.penalty.large=100; % infty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TargetDir='./mat/Models/';
mSAP_main(model,ARGSet,model.name,Para,TargetDir); % attributed graph mining

% show results
load(sprintf('%s%s_model.mat',TargetDir,model.name),'record','ARGSet','time');
model=record(end).model;
match=getMatch(model,ARGSet,1,Para);
for i=1:size(match,2)
    show_MatchingResult_SIFTImg(ARGSet(i).arg.name,match(:,i)',Para,model.link);
end
