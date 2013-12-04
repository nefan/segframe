#include "cljmex.hpp"
 
#define cljmex_start() \
cljmex_mex_fun();\
cljmex_check_inargs(8,"gaussianTKsC");\
cljmex_check_varoutargs(6,"gaussianTKsC");\
cljmex_setup();\
\
mxAssert(!mxIsComplex(pargin[0]),"q_a_i should be of type :real");\
mxAssert(mxGetM(pargin[0]) == dim,"q_a_i should have dim rows");\
mxAssert(mxGetN(pargin[0]) > 0 && mxGetM(pargin[0]) > 0,"q_a_i should be non-empty");\
cljmex_real_matrix q_a_i;\
q_a_i.rows = mxGetM(pargin[0]);\
q_a_i.cols = mxGetN(pargin[0]);\
q_a_i.x = (cljmexEntry*)mxGetData(pargin[0]);\
\
mxAssert(!mxIsComplex(pargin[1]),"p_a_i should be of type :real");\
mxAssert(mxGetM(pargin[1]) == dim,"p_a_i should have dim rows");\
mxAssert(mxGetN(pargin[1]) > 0 && mxGetM(pargin[1]) > 0,"p_a_i should be non-empty");\
cljmex_real_matrix p_a_i;\
p_a_i.rows = mxGetM(pargin[1]);\
p_a_i.cols = mxGetN(pargin[1]);\
p_a_i.x = (cljmexEntry*)mxGetData(pargin[1]);\
\
mxAssert(!mxIsComplex(pargin[2]),"scales2 should be of type :real");\
mxAssert(mxGetM(pargin[2]) == R,"scales2 should have R rows");\
mxAssert(mxGetM(pargin[2]) == 1,"scales2 should be row vector");\
cljmex_real_matrix scales2;\
scales2.rows = mxGetM(pargin[2]);\
scales2.cols = mxGetN(pargin[2]);\
scales2.x = (cljmexEntry*)mxGetData(pargin[2]);\
\
mxAssert(!mxIsComplex(pargin[3]),"scaleweight2 should be of type :real");\
mxAssert(mxGetM(pargin[3]) == R,"scaleweight2 should have R rows");\
mxAssert(mxGetM(pargin[3]) == 1,"scaleweight2 should be row vector");\
cljmex_real_matrix scaleweight2;\
scaleweight2.rows = mxGetM(pargin[3]);\
scaleweight2.cols = mxGetN(pargin[3]);\
scaleweight2.x = (cljmexEntry*)mxGetData(pargin[3]);\
\
mxAssert(mxIsInt64(pargin[4]),"dim should be of type :int");\
mxAssert(mxGetN(pargin[4]) == 1 && mxGetM(pargin[4]) == 1,"dim should be single element");\
cljmexInt dim = (cljmexInt)(*((cljmexInt*)mxGetData(pargin[4])));\
\
mxAssert(mxIsInt64(pargin[5]),"Lq should be of type :int");\
mxAssert(mxGetN(pargin[5]) == 1 && mxGetM(pargin[5]) == 1,"Lq should be single element");\
cljmexInt Lq = (cljmexInt)(*((cljmexInt*)mxGetData(pargin[5])));\
\
mxAssert(mxIsInt64(pargin[6]),"Lp should be of type :int");\
mxAssert(mxGetN(pargin[6]) == 1 && mxGetM(pargin[6]) == 1,"Lp should be single element");\
cljmexInt Lp = (cljmexInt)(*((cljmexInt*)mxGetData(pargin[6])));\
\
mxAssert(mxIsInt64(pargin[7]),"R should be of type :int");\
mxAssert(mxGetN(pargin[7]) == 1 && mxGetM(pargin[7]) == 1,"R should be single element");\
cljmexInt R = (cljmexInt)(*((cljmexInt*)mxGetData(pargin[7])));\
\
pargout[0] = mxCreateDoubleMatrix(Lq,Lp,mxREAL);\
cljmex_real_matrix K__ij;\
K__ij.rows = Lq;\
K__ij.cols = Lp;\
K__ij.x = (cljmexEntry*)mxGetPr(pargout[0]);\
\
int D1Ks__ijbdimsa[] = {Lq,Lp,dim};\
cljmexInt D1Ks__ijbdimsIa[] = {Lq,Lq*Lp};\
pargout[1] = mxCreateNumericArray(3,D1Ks__ijbdimsa,mxDOUBLE_CLASS,mxREAL);\
cljmex_real_multidimarray D1Ks__ijb;\
D1Ks__ijb.dims = D1Ks__ijbdimsa;\
D1Ks__ijb.dimsI = D1Ks__ijbdimsIa;\
D1Ks__ijb.x = (cljmexEntry*)mxGetPr(pargout[1]);\
\
int D2Ks__ijbgdimsa[] = {Lq,Lp,dim,dim};\
cljmexInt D2Ks__ijbgdimsIa[] = {Lq,Lq*Lp,Lq*Lp*dim};\
pargout[2] = mxCreateNumericArray(4,D2Ks__ijbgdimsa,mxDOUBLE_CLASS,mxREAL);\
cljmex_real_multidimarray D2Ks__ijbg;\
D2Ks__ijbg.dims = D2Ks__ijbgdimsa;\
D2Ks__ijbg.dimsI = D2Ks__ijbgdimsIa;\
D2Ks__ijbg.x = (cljmexEntry*)mxGetPr(pargout[2]);\
\
int D3Ks__ijbgddimsa[] = {Lq,Lp,dim,dim,dim};\
cljmexInt D3Ks__ijbgddimsIa[] = {Lq,Lq*Lp,Lq*Lp*dim,Lq*Lp*dim*dim};\
pargout[3] = mxCreateNumericArray(5,D3Ks__ijbgddimsa,mxDOUBLE_CLASS,mxREAL);\
cljmex_real_multidimarray D3Ks__ijbgd;\
D3Ks__ijbgd.dims = D3Ks__ijbgddimsa;\
D3Ks__ijbgd.dimsI = D3Ks__ijbgddimsIa;\
D3Ks__ijbgd.x = (cljmexEntry*)mxGetPr(pargout[3]);\
\
int D4Ks__ijbgdedimsa[] = {Lq,Lp,dim,dim,dim,dim};\
cljmexInt D4Ks__ijbgdedimsIa[] = {Lq,Lq*Lp,Lq*Lp*dim,Lq*Lp*dim*dim,Lq*Lp*dim*dim*dim};\
pargout[4] = mxCreateNumericArray(6,D4Ks__ijbgdedimsa,mxDOUBLE_CLASS,mxREAL);\
cljmex_real_multidimarray D4Ks__ijbgde;\
D4Ks__ijbgde.dims = D4Ks__ijbgdedimsa;\
D4Ks__ijbgde.dimsI = D4Ks__ijbgdedimsIa;\
D4Ks__ijbgde.x = (cljmexEntry*)mxGetPr(pargout[4]);\
\
int D5Ks__ijbgdepdimsa[] = {Lq,Lp,dim,dim,dim,dim,dim};\
cljmexInt D5Ks__ijbgdepdimsIa[] = {Lq,Lq*Lp,Lq*Lp*dim,Lq*Lp*dim*dim,Lq*Lp*dim*dim*dim,Lq*Lp*dim*dim*dim*dim};\
pargout[5] = mxCreateNumericArray(7,D5Ks__ijbgdepdimsa,mxDOUBLE_CLASS,mxREAL);\
cljmex_real_multidimarray D5Ks__ijbgdep;\
D5Ks__ijbgdep.dims = D5Ks__ijbgdepdimsa;\
D5Ks__ijbgdep.dimsI = D5Ks__ijbgdepdimsIa;\
D5Ks__ijbgdep.x = (cljmexEntry*)mxGetPr(pargout[5]);\
\

 
#define cljmex_end() \
\
\
\
\
\
\
\
}\

 
