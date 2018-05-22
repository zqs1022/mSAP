# Developed by Quanshi Zhang

# License & disclaimer

Copyright 2014 Quanshi Zhang (zqs1022@gmail.com, zhangqs@g.ucla.edu)
This software can be used for research purposes only. 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

This code applies "Tree-reweighted max-product message passing algorithm (TRW-S), verstion 1.3", which is developed by Microsoft Research. You need to agree the Microsoft Research Shared Source license agreement ("MSR-SSLA") before using it. Please read "./TRW_S/LICENSE.TXT" for details.


# Compile the code.

run compile;

#Use of the software.
Please directly run the demo. The initial model only has three nodes (object parts). The mined maximal-size SAP will be saved in the directory of "./mat/Models/".

run demo;


# Parallel computing

This mining process (mSAP_main.p) supports parallel computing. You can simply open some matlab workers before running the program.
For example,

run matlabpool(2);
run demo;


# Design of attributes

You can define your own unary and pairwise attributes (features) for your own applications. These features can be used as the input of the graph-mining algorithm. Please see mSAP_produce_SIFTImg.m to understand the data structure of the ARG.


# Citation

You should cite the following paper, if you use this software.

1. Quanshi Zhang, Xuan Song, Xiaowei Shao, Huijing Zhao, and Ryosuke Shibasaki, "Attributed Graph Mining and Matching: An Attempt to Define and Extract Soft Attributed Patterns", in Proc. of IEEE International Conference on Computer Vision and Pattern Recognition (CVPR) 2014.

2. Quanshi Zhang, Xuan Song, Xiaowei Shao, Huijing Zhao, and Ryosuke Shibasaki, "Object Discovery: Soft Attributed Graph Mining", to appear in IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI).


# Others

In fact, you can design your own unary and pairwise features, and use these features as the input of the graph-mining algorithm. In our demo, we apply the SIFT keypoint detector developed by David Lowe. If you use "FeatureExtraction_SIFTImg.m" to extract features, please cite David G. Lowe, "Distinctive image features from scale-invariant keypoints," International Journal of Computer Vision, 60, 2 (2004), pp. 91-110. Please see the website (http://www.cs.ubc.ca/~lowe/keypoints/) for details


In order to extract the SIFT features (the preprocessing of constructing the ARGs in this demo), please
1. go to the "mSAP_produce_SIFTImg" directory
2. run FeatureExtraction_SIFTImg('panda'); % only based on the windows system

In order to define your own unary and pairwise features, please see mSAP_produce_SIFTImg.m to understand the data structure of the ARG.
