function [] = fullwidth(fig, full_height)
% Make a figure full width or full screen. 
% Inputs:
%	fig,	figure number or handle to figure (optional)
%	full_height,	logical, default: false

	if ~exist('fig', 'var') || isempty(fig)
		fig = gcf;
	end
	if ~isobject(fig)
		try
			fig = figobj(fig);
		catch me
			full_height = fig;
			fig = gcf;
		end
	end
	
	if ~exist('full_height', 'var') || isempty(full_height)
		full_height = false;
	end
	
    set(fig, 'units', 'normalized');
    pos = get(fig, 'position');
	if full_height
		set(fig, 'position', [0 0 1 1]);
	else
		set(fig, 'position', [0 pos(2) 1 pos(4)]);
	end
end
