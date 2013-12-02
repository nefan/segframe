#include "cljmex.hpp"
 
#define cljmex_start() \
cljmex_mex_fun();\
cljmex_check_inargs(6,"gaussianTKsC");\
cljmex_check_varoutargs(6,"gaussianTKsC");\
cljmex_setup();\
\
mxAssert(!mxIsComplex(pargin[0]),"q0_a_i should be of type :real");\
mxAssert(mxGetM(pargin[0]) == dim,"q0_a_i should have dim rows");\
mxAssert(mxGetN(pargin[0]) > 0 && mxGetM(pargin[0]) > 0,"q0_a_i should be non-empty");\
cljmex_real_matrix q0_a_i;\
q0_a_i.rows = mxGetM(pargin[0]);\
q0_a_i.cols = mxGetN(pargin[0]);\
q0_a_i.x = (cljmexEntry*)mxGetData(pargin[0]);\
\
mxAssert(!mxIsComplex(pargin[1]),"scales2 should be of type :real");\
mxAssert(mxGetM(pargin[1]) == R,"scales2 should have R rows");\
mxAssert(mxGetM(pargin[1]) == 1,"scales2 should be row vector");\
cljmex_real_matrix scales2;\
scales2.rows = mxGetM(pargin[1]);\
scales2.cols = mxGetN(pargin[1]);\
scales2.x = (cljmexEntry*)mxGetData(pargin[1]);\
\
mxAssert(!mxIsComplex(pargin[2]),"scaleweight2 should be of type :real");\
mxAssert(mxGetM(pargin[2]) == R,"scaleweight2 should have R rows");\
mxAssert(mxGetM(pargin[2]) == 1,"scaleweight2 should be row vector");\
cljmex_real_matrix scaleweight2;\
scaleweight2.rows = mxGetM(pargin[2]);\
scaleweight2.cols = mxGetN(pargin[2]);\
scaleweight2.x = (cljmexEntry*)mxGetData(pargin[2]);\
\
mxAssert(mxIsInt64(pargin[3]),"dim should be of type :int");\
mxAssert(mxGetN(pargin[3]) == 1 && mxGetM(pargin[3]) == 1,"dim should be single element");\
cljmexInt dim = (cljmexInt)(*((cljmexInt*)mxGetData(pargin[3])));\
\
mxAssert(mxIsInt64(pargin[4]),"L should be of type :int");\
mxAssert(mxGetN(pargin[4]) == 1 && mxGetM(pargin[4]) == 1,"L should be single element");\
cljmexInt L = (cljmexInt)(*((cljmexInt*)mxGetData(pargin[4])));\
\
mxAssert(mxIsInt64(pargin[5]),"R should be of type :int");\
mxAssert(mxGetN(pargin[5]) == 1 && mxGetM(pargin[5]) == 1,"R should be single element");\
cljmexInt R = (cljmexInt)(*((cljmexInt*)mxGetData(pargin[5])));\
\
pargout[0] = mxCreateDoubleMatrix(L,L,mxREAL);\
cljmex_real_matrix K__ij;\
K__ij.rows = L;\
K__ij.cols = L;\
K__ij.x = (cljmexEntry*)mxGetPr(pargout[0]);\
\
int D1Ks__ijbdimsa[] = {L,L,dim};\
cljmexInt D1Ks__ijbdimsIa[] = {L,L*L};\
pargout[1] = mxCreateNumericArray(3,D1Ks__ijbdimsa,mxDOUBLE_CLASS,mxREAL);\
cljmex_real_multidimarray D1Ks__ijb;\
D1Ks__ijb.dims = D1Ks__ijbdimsa;\
D1Ks__ijb.dimsI = D1Ks__ijbdimsIa;\
D1Ks__ijb.x = (cljmexEntry*)mxGetPr(pargout[1]);\
\
int D2Ks__ijbgdimsa[] = {L,L,dim,dim};\
cljmexInt D2Ks__ijbgdimsIa[] = {L,L*L,L*L*dim};\
pargout[2] = mxCreateNumericArray(4,D2Ks__ijbgdimsa,mxDOUBLE_CLASS,mxREAL);\
cljmex_real_multidimarray D2Ks__ijbg;\
D2Ks__ijbg.dims = D2Ks__ijbgdimsa;\
D2Ks__ijbg.dimsI = D2Ks__ijbgdimsIa;\
D2Ks__ijbg.x = (cljmexEntry*)mxGetPr(pargout[2]);\
\
int D3Ks__ijbgddimsa[] = {L,L,dim,dim,dim};\
cljmexInt D3Ks__ijbgddimsIa[] = {L,L*L,L*L*dim,L*L*dim*dim};\
pargout[3] = mxCreateNumericArray(5,D3Ks__ijbgddimsa,mxDOUBLE_CLASS,mxREAL);\
cljmex_real_multidimarray D3Ks__ijbgd;\
D3Ks__ijbgd.dims = D3Ks__ijbgddimsa;\
D3Ks__ijbgd.dimsI = D3Ks__ijbgddimsIa;\
D3Ks__ijbgd.x = (cljmexEntry*)mxGetPr(pargout[3]);\
\
int D4Ks__ijbgdedimsa[] = {L,L,dim,dim,dim,dim};\
cljmexInt D4Ks__ijbgdedimsIa[] = {L,L*L,L*L*dim,L*L*dim*dim,L*L*dim*dim*dim};\
pargout[4] = mxCreateNumericArray(6,D4Ks__ijbgdedimsa,mxDOUBLE_CLASS,mxREAL);\
cljmex_real_multidimarray D4Ks__ijbgde;\
D4Ks__ijbgde.dims = D4Ks__ijbgdedimsa;\
D4Ks__ijbgde.dimsI = D4Ks__ijbgdedimsIa;\
D4Ks__ijbgde.x = (cljmexEntry*)mxGetPr(pargout[4]);\
\
int D5Ks__ijbgdepdimsa[] = {L,L,dim,dim,dim,dim,dim};\
cljmexInt D5Ks__ijbgdepdimsIa[] = {L,L*L,L*L*dim,L*L*dim*dim,L*L*dim*dim*dim,L*L*dim*dim*dim*dim};\
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

 
