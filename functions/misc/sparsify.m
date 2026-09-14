function [newA,doSparse,sTime,fTime] = sparsify(A,opt)
% SPARSIFY Returns the matrix A in sparse/full form depending on which form is more efficient.
% 
% If sparse(A) requires more storage/memory than full(A), returns full(A).
% Otherwise if sparse(A)*b is faster than full(A)*b, returns sparse(A).
% 
% The optional arguments "testFunc" and "testTime" can be used to change the timing comparison 
% function and the allotted comparison time respectively.
arguments
	A double {mustBeMatrix}
	opt.testTime (1,1) double {mustBePositive} = 0.5;
	opt.testFunc function_handle = @mtimes;
end

% check the test function
testFunc = opt.testFunc;
numArgs = nargin(testFunc);
if (numArgs < 1) || (numArgs > 2)
	error("ERROR: the test function must take either 1 or 2 arguments");
elseif (numArgs == 1)
	% check valid evaluation for matrix functions
	for nn = [2,3,4,5,10,100]
		try
			tmp = sprand(nn,nn,0.25,"double");
			[~] = testFunc(full(tmp));
			[~] = testFunc(tmp);
		catch ME
			error("ERROR: the test function fails for matrix inputs");
		end
	end
elseif (numArgs == 2)
	% check valid evaluation for matrix-vector functions
	for nn = [2,3,4,5,10,100]
		try
			tmp = sprand(nn,nn,0.25,"double");
			[~] = testFunc(full(tmp),randn(nn,1));
			[~] = testFunc(tmp,randn(nn,1));
		catch ME
			error("ERROR: the test function fails for matrix inputs");
		end
	end
end

% handle "out of memory" errors when converting between full/sparse
try
	fA = full(A);
	sA = sparse(A);
catch ME
	% default to current format
	newA = A;
	doSparse = issparse(A);
	sTime = NaN;
	fTime = NaN;
	return
end

% check the memory requirements
sMem = whos("sA").bytes;
fMem = whos("fA").bytes;

if (sMem >= fMem)
	% defualt to the full-form
	doSparse = false;
	sTime = NaN;
	fTime = NaN;
else
	% compare the efficiency of the test-function for sparse/full inputs
	sTime = 0;
	fTime = 0;
	
	iters = 0;
	isOverTime = false;
	
	tic;
	while ~isOverTime
		% compare compute time for a ...
		if (numArgs == 1)
			% matrix function
			t0 = toc;
			[~] = testFunc(sA);
			t1 = toc;
			[~] = testFunc(fA);
			t2 = toc;
		elseif (numArgs == 2)
			% matrix-vector function
			vec = randn(size(A,2),1);
			t0 = toc;
			[~] = testFunc(sA,vec);
			t1 = toc;
			[~] = testFunc(fA,vec);
			t2 = toc;
		end
		
		% add to total time
		sTime = sTime + (t1 - t0);
		fTime = fTime + (t2 - t1);
		
		% exit loop if next test would exceed allotted time
		iters = iters + 1;
		isOverTime = (toc * (iters+1)/iters > opt.testTime);
	end

	% return the average sparse/full timings
	sTime = sTime/iters;
	fTime = fTime/iters;

	% return the optimal format
	doSparse = (sTime <= fTime);
end

% return the optimal matrix format
if doSparse; newA = sA; else; newA = fA; end

% DONE
end