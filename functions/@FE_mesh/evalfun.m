function [data] = evalfun(obj,expr)
% FE_mesh.func_eval: Evaluates the given function over the mesh
arguments
	obj FE_mesh
	expr {mustBeA(expr,{'function_handle','char','string','numeric','logical'})}
end

% skip evaluation if already numeric
if isnumeric(expr) || islogical(expr)
	data = double(expr);
else
	% evaluate the expression at the simplex vertices
	pts = obj.ptsPoly;
	tri = obj.triPoly;
	
	% load all pre-defined variables for functional expressions
	def_var = {{"pts",pts},{"tri",tri},{"reg",obj.reg}};
	numDim = obj.num("dim");
	
	% using actual coordinates
	def_var{end+1} = {"r",vecnorm(pts,2,1)};
	if numDim >= 1; def_var{end+1} = {"x", pts(1,:)}; end
	if numDim >= 2; def_var{end+1} = {"y", pts(2,:)}; end
	if numDim >= 3; def_var{end+1} = {"z", pts(3,:)}; end
	
	% using normalized/centered coordinates
	ptsMax = max(obj.pts,[],2);
	ptsMin = min(obj.pts,[],2);
	normpts = 2/max(ptsMax - ptsMin) * (pts - (ptsMax+ptsMin)/2);
	def_var{end+1} = {"rr",vecnorm(normpts,2,1)};
	if numDim >= 1; def_var{end+1} = {"xx", normpts(1,:)}; end
	if numDim >= 2; def_var{end+1} = {"yy", normpts(2,:)}; end
	if numDim >= 3; def_var{end+1} = {"zz", normpts(3,:)}; end
	
	% attempt evaluation of the inputted coefficient
	[expr,args,isValid] = fparse(expr,def_var);
	if isValid
		data = feval(expr,args{:});
	else
		error("ERROR: the specified function is not evaluatable over the mesh")
	end
end
% DONE
end