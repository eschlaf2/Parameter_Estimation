
xTrue = [abs(EEGanalytic).*cos(angle(EEGanalytic)) abs(EEGanalytic) unwrap(angle(EEGanalytic))];
fs = 3e3;
tSpan = (0:(length(xTrue)-1)) / fs;
sqrtR = 0.04;
yMeas = xTrue(:,1) + sqrtR*randn(numel(tSpan),1);

inds = 10001:20000;
observations = xTrue(inds, [2 3])';
diffs = diag([10;1.1]) * diff(observations, [], 2);

myPF = particleFilter(...
	@(particles) transitionFcn(particles, 'distribution', diffs), ...
	@(particles, measurement) likelihoodFcn(particles, measurement, ...
										'observations', observations));
									

initialize(myPF, 1000, mean(observations, 2), diag([1;10]));
myPF.ResamplingMethod = 'systematic';
myPF.ResamplingPolicy.MinEffectiveParticleRatio = 0.2;

xEst = zeros(size(xTrue, 1), 2);
for k=1:1e4
    xEst(k,:) = correct(myPF,yMeas(k));
    predict(myPF);
end

figure(3);
plot(tSpan(1:k), xTrue(1:k, 1)','b',tSpan(1:k), xEst(1:k,1)'.*cos(xEst(1:k,2))','r');
legend('True','Estimated');