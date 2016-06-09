function [G] = computeIntraMatrix(intra);

G = zeros(6);
n_intra = size(intra, 1);

for k = 1 : n_intra,
  i1 = intra(k, 1);           % Source cell index.
  j1 = intra(k, 2);           % Target cell index.
  strength = intra(k, 3);
  i2 = mod(i1 + 2, 6) + 1;    % Source cell index for symmetric connection.
  j2 = mod(j1 + 2, 6) + 1;    % Target cell index for symmetric connection.
  
  G(j1, i1) = strength;
  G(j2, i2) = strength;
end
end