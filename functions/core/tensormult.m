function [C] = tensormult(A,indexA,B,indexB,indexC)
% Computes the indexed tensor product
% Labeling an index with the characters '~','-', or '_' will denote the corresponding
% dimension is singular, i.e. the length of the tensor along that dimension is one.
%
% If an index occurs in A and B, and does NOT occur in C, the result is an inner-product.
% If an index occurs in A xor B, and also occurs in C, the result is an outer-product.
% If an index occurs in A and B, and also occurs in C, the result is an element-wise product.
%
% for example: 
%   C = tensormult(A,'abjn',B,'abkn','njk'); 
% yields the tensor:
%    C(n,j,k) = sum_{a,b} { A(a,b,j,n) * B(a,b,k,n) }
%
% the matrix-vector product x=A*b is given by:
%   x = tensormult(A,'ij',b,'j','i');
%
% whereas x = b*A is given by:
%   x = tensormult(b,'~i',A,'ij','~j');
%
% SEE ALSO: TENSORPROD
arguments
	A double {mustBeNonempty}
	indexA (1,:) {mustBeA(indexA,{'string','char','numeric'})}
	B double {mustBeNonempty}
	indexB (1,:) {mustBeA(indexB,{'string','char','numeric'})}
	indexC (1,:) {mustBeA(indexC,{'string','char','numeric'})}
end

% define faster method for unique(x,"stable")
isunique = @(x) ~any(triu(x(:).'==x(:),1),1);

%% convert index labels to string arrays
%===========================================================
if ~isstring(indexA)
	tmp = strings(size(indexA));
	for ii = 1:numel(tmp)
		tmp(ii) = string(indexA(ii));
	end
	indexA = tmp;
end
if ~isstring(indexB)
	tmp = strings(size(indexB));
	for ii = 1:numel(tmp)
		tmp(ii) = string(indexB(ii));
	end
	indexB = tmp;
end
if ~isstring(indexC)
	tmp = strings(size(indexC));
	for ii = 1:numel(tmp)
		tmp(ii) = string(indexC(ii));
	end
	indexC = tmp;
end

%% remove singular dimensions
%===========================================================
% check indices for singular-dimension labels
singularA = matches(indexA,["~","-","_"]);
singularB = matches(indexB,["~","-","_"]);
singularC = matches(indexC,["~","-","_"]);

% check that singular labels agree with the tensor sizes
if any(size(A,find(singularA))~=1)
	error("ERROR: 'A' contains non-singular dimensions in locations given by 'indexA'.")
end
if any(size(B,find(singularB))~=1)
	error("ERROR: 'B' contains non-singular dimensions in locations given by 'indexB'.")
end

% remove the singular dimensions of A and indexA
if any(singularA)
	indexA = indexA(~singularA);
	A = permute(A,[find(~singularA),find(singularA)]);
end

% remove the singular dimensions of B and indexB
if any(singularB)
	indexB = indexB(~singularB);
	B = permute(B,[find(~singularB),find(singularB)]);
end

% remove the singular dimensions of indexC
outIndex = indexC;
indexC = indexC(~singularC);


%% determine the product-type for each index label
%===========================================================
% check that the index labels are unique
if any(~isunique(indexA)); error("ERROR: non-unique labels given in 'indexA'."); end
if any(~isunique(indexB)); error("ERROR: non-unique labels given in 'indexB'."); end
if any(~isunique(indexC)); error("ERROR: non-unique labels given in 'indexC'."); end

% define index ordering
allIndx = [indexC,indexA,indexB];
allIndx = allIndx(isunique(allIndx));

% define index locations
in_A = any(allIndx == indexA(:),1);
in_B = any(allIndx == indexB(:),1);
in_C = any(allIndx == indexC(:),1);

% check for unpaired index labels in the input/output
if any(~or(in_A,in_B) & in_C)
	% unpaired labels in the output
	error("ERROR: index labels in 'indexC' are missing from 'indexA' and/or 'indexB'.")
end
if any(xor(in_A,in_B) & ~in_C)
	% unpaired labels in the inputs
	error("ERROR: index labels in either 'indexA' or 'indexB' are missing from 'indexC'.")
end

% define the product-type for each index label 
isJoint = and(in_A, in_B) &  in_C;
isInner = and(in_A, in_B) & ~in_C;
isOuter = xor(in_A, in_B) &  in_C;

%% get the tensor dims/sizes for each product type
%===========================================================
% inner products
innerIndx = allIndx(isInner);
[iDimsA,~] = find(innerIndx == indexA(:));
[iDimsB,~] = find(innerIndx == indexB(:));
iSizeA = size(A,iDimsA);
iSizeB = size(B,iDimsB);

% outer products
outerIndx = allIndx(isOuter);
[oDimsA,~] = find(outerIndx == indexA(:));
[oDimsB,~] = find(outerIndx == indexB(:));
oSizeA = size(A,oDimsA);
oSizeB = size(B,oDimsB);

% joint products
jointIndx = allIndx(isJoint);
[jDimsA,~] = find(jointIndx == indexA(:));
[jDimsB,~] = find(jointIndx == indexB(:));
jSizeA = size(A,jDimsA);
jSizeB = size(B,jDimsB);

% check that the size of the tensors agree along the inner dimensions
if ~isempty(innerIndx) & any((iSizeA ~= iSizeB) & (iSizeA~=1) & (iSizeB~=1))
	error("ERROR: the sizes of the contracted dimensions in 'A' and 'B' do not agree")
else
	iSizes = max(iSizeA,iSizeB);
end

% check that the size of the tensors agree along the joint dimensions
if ~isempty(jointIndx) && any((jSizeA ~= jSizeB) & (jSizeA~=1) & (jSizeB~=1))
	error("ERROR: the sizes of the non-contracted dimensions in 'A' and 'B' do not agree.")
else
	jSizes = max(jSizeA,jSizeB);
end

%% expand inner-product dimensions
%===========================================================
iRepsA = iSizes./iSizeA;
if prod([iRepsA,1]) > 1
	% expand A's inner-dims
	repeats = ones(1,max([ndims(A);iDimsA]));
	repeats(iDimsA) = iRepsA;
	A = repmat(A,repeats);
end

iRepsB = iSizes./iSizeB;
if prod([iRepsB,1]) > 1
	% expand B's inner-dims
	repeats = ones(1,max([ndims(B);iDimsB]));
	repeats(iDimsB) = iRepsB;
	B = repmat(B,repeats);
end

%% expand joint-product dimensions
%===========================================================
jRepsA = jSizes./jSizeA;
if prod([jRepsA,1]) > 1
	% expand A's joint-dims
	repeats = ones(1,max([ndims(A);jDimsA]));
	repeats(jDimsA) = jRepsA;
	A = repmat(A,repeats);
end

jRepsB = jSizes./jSizeB;
if prod([jRepsB,1]) > 1
	% expand B's joint-dims
	repeats = ones(1,max([ndims(B);jDimsB]));
	repeats(jDimsB) = jRepsB;
	B = repmat(B,repeats);
end

%% re-order tensors for multiplication
%===========================================================
% permute tensors so inner-product dims are contracted by pagemtimes, 
% while joint-product dims form the pages of each tensor
A = safe_permute(A,[oDimsA; iDimsA; jDimsA]);
A = safe_reshape(A,[prod([oSizeA,1]),prod([iSizes,1]),prod([jSizes,1])]);

B = safe_permute(B,[iDimsB; oDimsB; jDimsB]);
B = safe_reshape(B,[prod([iSizes,1]),prod([oSizeB,1]),prod([jSizes,1])]);

%% multiply and re-expand the outer/joint product dims
%===========================================================
if isempty(jDimsA) && isempty(jDimsB)
	% use basic matrix-multiplication (no element-wise products)
	C = A*B;
else
	% use page-wise multiplication
	if issparse(A) || issparse(B)
		warning("The inputs cannot be sparse when computing element-wise " + ...
			"products along some dimension. Converting to a full matrix.")
		A = full(A);
		B = full(B);
	end
	C = pagemtimes(A,B);
end
C = safe_reshape(C,[oSizeA,oSizeB,jSizes]);

%% re-order dims to match the output indices
%===========================================================
% find the correct dimension order
currIndexOrder = [indexA(oDimsA),indexB(oDimsB),jointIndx];
[oldDimLoc,newDimLoc] = find(currIndexOrder(:) == outIndex);
dimOrderC = accumarray(newDimLoc,oldDimLoc,[numel(outIndex),1]).';

% add in any singular dimensions
dimOrderC(~dimOrderC) = (max(dimOrderC)+1:numel(dimOrderC));

% apply new ordering;
C = safe_permute(C,dimOrderC);

%% DONE
end

% external functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A] = safe_permute(A,dimOrder)
% permutes safely if dimOrder is missing trailing dims
% really only necessary because permute requires dimOrder 
% to have at least two elements, even if A is a vector
trailingDims = ((numel(dimOrder)+1):ndims(A));
dimOrder = [dimOrder(:)', trailingDims];
A = permute(A,dimOrder);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A] = safe_reshape(A,dimSizes)
% reshapes safely if dimSizes is missing trailing dims
% really only necessary because reshape requires dimSizes 
% to have at least two elements, even if A is a vector 
trailingSizes = ones(1,max(0,ndims(A)-numel(dimSizes)));
dimSizes = [dimSizes(:).',trailingSizes];
A = reshape(A,dimSizes);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%