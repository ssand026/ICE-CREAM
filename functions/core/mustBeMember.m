function mustBeMember(A,B)
% MUSTBEMEMBER(overloaded) Validates values are members of the given set
%	MUSTBEMEMBER(A,B) throws an error if any A are not members of the set B.
%   MATLAB calls ismember(A,B) to determine if A is a member of B.
%
% Overloads the default MUSTBEMEMBER, using a faster method for validating
% string/character inputs while maintaining autocompletion support
%
% SEE ALSO: ISMEMBER, MATCHES

try
	% check if sets have compatible text formats
	isStrA = isstring(A);
	isStrB = isstring(B);
	isCharA = ischar(A) | iscellstr(A); %#ok
	isCharB = ischar(B) | iscellstr(B); %#ok

	if (isStrA && isStrB) || (isCharA && isCharB)
		% sets A and B have compatible text formats, validate membership using "matches"
		isValid = all(matches(A,B),"all");
	% elseif (~isStrA && isStrB) || (~isCharA && isCharB)
	% 	% sets A and B have incompatible text formats
	% 	isValid = false;
	else
		% use default validation method
		isValid = all(ismember(A,B),"all");
	end
catch ME
	isValid = false;
end

% throw error if elements of A are not valid members of B
if ~isValid
	throwAsCaller(...
		matlab.internal.validation.util.createValidatorExceptionWithValue( ...
		matlab.internal.validation.util.createPrintableList(B), ...
		'MATLAB:validators:mustBeMemberGenericText',...
		'MATLAB:validators:mustBeMember')...
		);
end
% DONE
end