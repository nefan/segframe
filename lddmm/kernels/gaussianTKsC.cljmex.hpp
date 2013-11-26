#include "cljmex.hpp"
 
#define cljmex_start() \
cljmex_mex_fun();\
cljmex_check_argsLE(6,4,"gaussianTKsC");\
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

 
#define cljmex_end() \
\
\
}\

 
