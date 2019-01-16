function window = updateWindow(prediction, W, transitionFcn)

	window = zeros([size(prediction) W]);

	window(:, :, 1) = prediction;
	for i = 2:W
		window(:, :, i) = transitionFcn(window(:, :, i - 1));
	end
		
end