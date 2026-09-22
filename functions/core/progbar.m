function [dispBar] = progbar(label,style,barWidth)
% Returns a function-handle that displays a progress bar in the command window. 
% 
% Calling dispBar(0) will initialize the progress bar, and subsequent calls
% dispBar(r) will update the displayed text so that the displayed progress bar 
% corresponds to the fractional progress ratio "r". 
%
% Unfortunately, any print statements between subsequent calls of "dispBar"
% can ruin the progress bar formatting. However this method is useful when 
% logging code outputs or when the overhead of waitbar is too large. 
%
% SEE ALSO: WAITBAR
arguments
	label char {mustBeTextScalar}
	style (1,1) {mustBeMember(style,["simple","fancy","explicit"])}
	barWidth (1,1) double {mustBeInteger,mustBePositive}
end

switch style
	case "simple"
		% simplest progress bar
		stepNum = 0;
		dispBar = @update_simple;
		
	case "fancy"
		% fanciest progress bar
		len = 0;
		stepNum  = 0;
		dispBar = @update_fancy;
		
	case "explicit"
		% explicit progress bar
		len = 0;
		dispBar = @update_explicit;
end

% internal functions
%***********************************************************
	function update_simple(r)
	% simple progress bar with no backspace characters
	if (r * barWidth) > (stepNum+1)
		% intermediate display
		fprintf('#')
		stepNum = stepNum + 1;
	elseif (r == 0)
		% initial display
		fprintf([repmat(' ',1,length(sprintf(label))),' ',repmat('_',1,barWidth),'\n']);
		fprintf([label,'['])
	elseif (r == 1)
		% final display
		fprintf('#')
		fprintf('] [100%%] \n')
	end
	% DONE
	end
%***********************************************************
	function update_fancy(r)
	% fancy progress bar with inline updating
	seq = [' ','▏','▎','▍','▌','▋','▊','▉','█']; % display sequence
	n = length(seq); % number of sequence steps
	
	% check if bar requires update
	step = floor(r * max((n-1)*barWidth,1000));
	doUpdate = (step ~= stepNum);
	stepNum = step; % update step number

	if (r == 0)
		% initial display
		pcnt = '  0';
		bar = repmat(seq(1),1,barWidth);
		bar = ['|',bar,'|  [',pcnt,'%%',']\n'];
		
		fprintf(repmat('\b',1,len))
		len = fprintf([label,bar]);
		
	elseif (r == 1)
		% final display
		pcnt = '100';
		bar = repmat(seq(n),1,barWidth);
		bar = ['|',bar,'|  [',pcnt,'%%',']\n'];
		
		fprintf(repmat('\b',1,len))
		len = fprintf([label,bar]);
		
	elseif doUpdate
		% intermediate display
		pcnt = sprintf('%4.3g',100*r);
		
		numL = floor(r * barWidth);
		numR = (barWidth-numL-1);
		seqIndex = 1 + floor((n-1) * mod(r*barWidth,1));
		
		bar = [repmat(seq(n),1,numL),seq(seqIndex),repmat(seq(1),1,numR)];
		bar = ['|',bar,'|  [',pcnt,'%%',']\n'];
		
		fprintf(repmat('\b',1,len))
		len = fprintf([label,bar]);
	end
	% DONE
	end
%***********************************************************
	function update_explicit(r)
	% explicit display with inline updating
	% when using this method, barLen should equal total iteration number
	fprintf(repmat('\b',1,len));
	len = fprintf(label+"[ "+round(r*barWidth)+" / "+barWidth+" ]\n");
	% DONE
	end
%***********************************************************
% DONE
end