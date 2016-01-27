disp('COMPILING MEX FILES...')

cd mex

disp('==> compiling ''Poses2Trans_mex.cpp'' (1 out of 4)');
mex Poses2Trans_mex.cpp COMPFLAGS="/openmp $COMPFLAGS"

disp('==> compiling ''CreateSet_mex.cpp'' (2 out of 4)');
mex CreateSet_mex.cpp

disp('==> compiling ''EvaluateEa_color_mex.cpp'' (3 out of 4)');
mex EvaluateEa_color_mex.cpp COMPFLAGS="/openmp $COMPFLAGS"

disp('==> compiling ''EvaluateEa_photo_mex.cpp'' (4 out of 4)');
mex EvaluateEa_photo_mex.cpp COMPFLAGS="/openmp $COMPFLAGS"

disp('==> DONE!');
fprintf('\n');
cd ..