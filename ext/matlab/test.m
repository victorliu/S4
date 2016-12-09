clear all;

S = S4Simulation([1 0; 0 1], 25);
mVac = S.setMaterial(1);
mSi = S.setMaterial(12);

L0 = S.addLayer(0, mVac);
Lslab = S.addLayer(0.5, mSi);
Lslab.setRegion(mVac, 'circle', 0.2);
LN = S.addLayer(0, mVac);

S.setPlanewave([0 0 1], [0 1 0], 1.0, 0.0);

freqs = 0.2:0.001:0.8;
tra = zeros(size(freqs));
ref = zeros(size(freqs));
for ifreq = 1:length(freqs)
	S.setFrequency(freqs(ifreq));
	[f,b] = LN.getPowerFlux();
	tra(ifreq) = f;
	ref(ifreq) = b;
end

delete(S);
plot(freqs, tra);
