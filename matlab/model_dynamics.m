function F = model_dynamics(~, stateMat, paramStruct, model, p)
% state should be a matrix


if ~exist('model', 'var') || isempty(model)
	model = 'HH';
end

if ~exist('p', 'var') || isempty(p)
	p = default_params(model);
end


switch model
	case 'Izh'  % Izhikevich model

		%% Parse input
		v = stateMat(1, :);
		u = stateMat(2, :);

		if exist('paramStruct', 'var')
			for f = fieldnames(paramStruct)'
				p.(f{:}) = paramStruct.(f{:});
			end
		end

		%% Calculate changes
		vdot = 0.04 * v.^2 + 5 * v + 140 - u + p.I;
		udot = p.a .* (p.b .* v - u);

		F = [vdot; udot]; % state change
		
		
	case 'HH'  % Hodgkin-Huxley model
		
		%% Parse input
		V = stateMat(1, :);
		n = stateMat(2, :);
		h = stateMat(3, :);
		B = stateMat(4, :);

		for f = fieldnames(paramStruct)'
			p.(f{:}) = paramStruct.(f{:});
		end

		%% Set parameters

		mInf = 1./(1 + exp((-V-34.5)/10));	% HH options
		nInf = 1./(1 + exp((-V-29.5)/10));
		hInf = 1./(1 + exp((V+59.4)/10.7));
		tauN = 0.25 + 4.35*exp(-abs(V+10)/10);
		tauH = 0.15 + 1.15./(1 + exp((V+33.5)/15));
		BInf = 1./(1 + exp(-(V-p.VBth)./p.SB));
		tauB = p.tauB;

		%% Calculate changes
		Vdot = (... % F1(V, n, h, B)
			p.I - p.gK .* n.^4 .* (V - p.EK) ... % drive current minus Potassium current
			- p.gNa .* mInf.^3 .* h .* (V - p.ENa) ... % Sodium current
			- p.gB .* B .* (V - p.EB) ... % mystery current
			- p.gL .* (V - p.EL)) / p.C; ... % Leak

		ndot = (nInf - n) ./ tauN; % F2(n)

		hdot = (hInf - h) ./ tauH; % F3(h)

		Bdot = (BInf - B) ./ tauB; % F4(B)

		try 
			F = [Vdot; ndot; hdot; Bdot]; % state change
		catch ME
			disp('bug')
		end

end
