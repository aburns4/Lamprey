function [l, w, Hs, psi] = Htune3_01(phi, H1, H2, H3, psi_bar, dw)
% [l w Hs psi] = Htune3_01(phi, H1, H2, H3, psi_bar, dw)
%
% Tuning with three coupling functions.
% Weighted average of coupling functions is adjusted with length so that
% stable zero is proportional to length.
%
% phi       Relative phases at which coupling functions are computed.
% H1,H2,H3  np-by-1 coupling functions.  Assumptions:
%             H2 has stable zero psi2.
%             r*H1 + (1 - r*H2) has stable zero psi(r) < psi2 for small r > 0.
%             r*H3 + (1 - r*H2) has stable zero psi(r) > psi2 for small r > 0.
% psi_bar   Desired phase lag per segment.
% dw        Weight step size used during tuning.
%
% l         nl-by-1 vector of lenghts at which tuning was achieved.
% w         nl-by-3 matrix of weights.
% Hs        np-by-nl matrix of tuned coupling functions.
% psi       nl-by-1 vector of stable zeros of tuned coupling functions.

% 10-Jun-2007  Tim Kiemel


r     = (0 : dw : 1)';
q     = length(r);

psi12 = zeros(q, 1);
psi23 = zeros(q, 1);

for i = 1 : q,
  psi12(i) = Hz_01(phi, (1 - r(i))*H1 + r(i)*H2);
  psi23(i) = Hz_01(phi, (1 - r(i))*H2 + r(i)*H3);
end

psi_min = min(psi12);
psi2    = Hz_01(phi, H2);
psi_max = max(psi23);

lmin = ceil(psi_min/psi_bar);
lmax = floor(psi_max/psi_bar);

l  = (lmin : lmax)';
nl = length(l);
w  = zeros(nl, 3);

for i = 1 : nl,
 psii = l(i)*psi_bar;

 if psii < psi2
   ri = Hz_01(r, psi12 - psii);
   w(i,:) = [1-ri  ri  0];
 else
   ri = Hz_01(r, psi23 - psii);
   w(i,:) = [0  1-ri  ri];
 end
end

Hs = [H1 H2 H3]*w';

psi = zeros(nl, 1);
for i = 1 : nl,
  psi(i) = Hz_01(phi, Hs(:,i));
end
