%%% compile mex files

% -D__QUIET__ 
% no openmp: -D_NO_OPENMP
% OPTIMFLAGS = /O2 /Oy- /DNDEBUG /openmp
% for more speed compile externally with other compiler

% download eigen library first:

% with openmp:
mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -I./ -I./Source/eigen/ -I./Source/helper/ -output ./ofmex ./Source/IRLS_Mex2.cpp

% no openmp:
% mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG" -O -I./ -D_NO_OPENMP -I./Source/eigen/ -I./Source/helper/ -output ./mex/ofmex ./Source/IRLS_Mex2.cpp 

%%%% no eigen lib: %%%% slower less memory

% with openmp:
%mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -I./ -D_noeigen_ -I./Source/helper/ -output ./mex/ofmex ./Source/IRLS_Mex2.cpp

% no openmp
%mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG" -O -I./ -D_noeigen_ -D_NO_OPENMP -I./Source/helper/ -output ./mex/ofmex ./Source/IRLS_Mex2.cpp 
