function [funcHandle, funcArgs, isValid] = fparse(func, argList)
% Checks if the function can be evaluated from the set of inputs 'argList', and
% returns the consolidated function arguments
% 
% Passes the defined variables in argList into the inputted function and 
% attempts to evaluate the function body. If unable to evaluate, returns the 
% reduced expression for the function.
arguments
	func {mustBeA(func,{'function_handle','char','string'})}
end
arguments (Repeating)
	argList (1,:) cell
end

% check validity of var_list
if isscalar(argList) && iscell(argList{1})
	argList = argList{1};
end

if ~all(cellfun(@numel,argList)==2)
	error("ERROR: the argument list has an invalid size")
elseif ~all(cellfun(@(x) isStringScalar(x{1}), argList))
	error("ERROR: the first arguments of argList must be a string scalars")
end

% extract the argument names/values from var_list
argNames = strings(1,numel(argList));
argValues = cell(1,numel(argList));
for jj = 1:numel(argList)
	argNames(jj) = argList{jj}{1};
	argValues{jj} = argList{jj}{2};
end
% check for repeats
[~,indx] = unique(argNames);
argNames = argNames(indx);
argValues = argValues(indx);


% append any internalized workspace variables to the input args
if isa(func,'function_handle')
	f_work = functions(func).workspace{1};
	workNames = string(fieldnames(f_work)).';
	workValues = struct2cell(f_work).';

	% add or replace input arguments, with priority given to workspace vars
	for ii = 1:numel(workNames)
		isTaken = (workNames(ii)==argNames);
		if any(isTaken)
			argValues(isTaken) = workValues{ii};
		else
			Narg = numel(argNames);
			argNames(Narg+1) = workNames(ii);
			argValues{Narg+1} = workValues{ii};
		end
	end
end

% convert function to string
if isa(func,'function_handle')
	f_string = string(functions(func).function);
else
	f_string = string(func);
end
f_string = regexprep(f_string,"\s*","");

% extract possible valid variable names from the function body
f_body = regexp(f_string,"(?<=@\([^)]+\)).*","match");
[f_args,f_notarg] = regexp(f_body,"([a-z_A-Z]+(\d+|))+","match","split");

% find variables in the function body that are defined by argList
def_args = intersect(f_args,argNames);

% find variables in the function body that are not defined or built-in/existing functions
undef_args = setdiff(f_args,def_args);
isafunc = arrayfun(@(var) any(exist(var)==[2,3,5]), undef_args);
undef_args = undef_args(~isafunc);

% replace any defined variables with an index into the cell array f_inputs
if ~isempty(def_args)
	f_inputs = cell(1,numel(def_args));
	for ii = 1:numel(def_args)
		f_args(f_args==def_args(ii)) = "varargin{"+ii+"}";
		f_inputs{ii} = argValues{argNames==def_args(ii)};
	end
	% rejoin function body
	f_body = f_notarg(1) + join(reshape([f_args; f_notarg(2:end)],1,[]),"");
end

% return the outputs
if isempty(f_body)
	% the function does not have a valid form
	funcHandle = [];
	funcArgs = {};
	isValid = false;
elseif isempty(def_args) & isempty(undef_args)
	% the function requires no arguments for evaluation
	funcHandle = str2func("@(~) " + f_body);
	funcArgs = {};
	isValid = true;
elseif ~isempty(def_args) & isempty(undef_args)
	% all of the required arguments have been defined
	funcHandle = str2func("@(varargin) " + f_body);
	funcArgs = f_inputs;
	isValid = true;
elseif ~isempty(def_args) & ~isempty(undef_args)
	% some of the required arguments have been defined
	funcHandle = str2func("@(" + join([undef_args,"varargin"],",") + ") " + f_body);
	funcArgs = f_inputs;
	isValid = false;	
elseif isempty(def_args) & ~isempty(undef_args)
	% none of the required arguments have been defined
	funcHandle = str2func("@(" + join(undef_args,",") + ") " + f_body);
	funcArgs = {};
	isValid = false;
end
% DONE
end