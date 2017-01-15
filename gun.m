function gun
% The "gun" problem is to find the eigenvalues of
% F(z) = K - z*M + 1i*sqrt(z-sigma1^2)*W1  (sigma1 = 0) 
%                + 1i*sqrt(z-sigma2^2)*W2  (sigma2 = 108.8774).
% F(z) is a nonlinear matrix-valued function, so the
% gun problem is a nonlinear eigenvalue problem.
% The eigenvalues can be computed by conversion to a polynomial
% eigenvalue problem as shown in Appendix B of thesis.pdf.
% The eigenvalues nearest the shift 146.71 are computed in 
% this program.

close all

%% Problem definition

[coeffs,fun] = nlevp('gun');
K = coeffs{1};
M = coeffs{2};
W1 = coeffs{3};
W2 = coeffs{4};
sigma2 = 108.8774;
    function f1z = f1(z)
        x = fun(z);
        f1z = x(3);
    end
    function f2z = f2(z)
        x = fun(z);
        f2z = x(4);
    end
F = @(L) K - L*M + f1(L)*W1 + f2(L)*W2; % original variable L
T = @(z) F(z^2);

%% Change of variable
% Noticing that 
% sqrt(z^2 - sigma2^2) = z*sqrt(1 - (sigma2/z)^2), 
% we turn
% this into a polynomial eigenvalue problem by a trigonometric change of
% variables. 
% 
% First, let sigma2/z = cos(theta). Then z = sigma2/cos(theta)
% and sqrt(1 - (sigma2/z)^2) = sin(theta). Then
% cos(theta)^2*T(sigma2/cos(theta)) = cos(theta)^2*K - sigma2^2*M +
% 1i*sigma2*cos(theta)*W1 + 1i*sigma2*cos(theta)*sin(theta)*W2.
% Rewriting the sines and cosines in terms of exp( +/- i theta) gives a
% linear combination of exp(-2i*theta) through exp(2i*theta), so
% multiplying by exp(2i*theta) gives a fourth-degree polynomial in
% exp(1i*theta). With w = exp(1i*theta), we get P, a fourth-degree
% polynomial in w. Altogether,
% z = sigma2/cos(theta) = 2*sigma2/(w + 1/w).

A0 = K/4 - 1i*sigma2*W2/4i;
A1 = 1i*sigma2*W1/2;
A2 = K/2 - sigma2^2*M;
A3 = 1i*sigma2*W1/2;
A4 = K/4 + 1i*sigma2*W2/4i;
P = @(w) A0 + A1*w + A2*w^2 + A3*w^3 + A4*w^4;

%%
% Check change of variables

z0 = rand+1i*rand;
L0 = z0^2;
theta0 = acos(sigma2/z0);
omega0 = exp(1i*theta0);
T0 = T(z0); 
F0 = F(L0); 
err = abs(F0 - T0);
fprintf('Max error between F and T: %4.2e\n', full(max(err(:))) );
P0 = P(omega0);
err = abs( P0 - exp(2i*theta0)*cos(theta0)^2*T0 );
fprintf('Max error between T and P: %4.2e\n', full(max(err(:))) );

%%
% Companion matrix formulation of the polynomial eigenvalue problem.
% Then we can use eigs().

I = speye(size(A0,1));
A_companion = [-A3, -A2, -A1, -A0; ...
                 I, 0*I, 0*I, 0*I; ...
               0*I,   I, 0*I, 0*I; ...
               0*I, 0*I,   I, 0*I];
B_companion = blkdiag(A4,I,I,I);

%% Finding sqrt(lambda) near 146.71
% We want z near c = 146.71. 
% Recall z = 2*sigma2/(w + 1/w), so z*w^2  - 2*sigma2*w + z0 = 0. Therefore
% w = (sigma2 +/- sqrt(sigma2^2 - z^2) )/z. So, if z near c, then w is near
% one of
% w = (sigma2 +/- sqrt(sigma2^2 - c^2) )/c.

zc = 146.71; % point of interest in z variable
wc1 = (sigma2 + sqrt(sigma2^2 - zc^2) )/zc; % " " " " w "
wc2 = (sigma2 - sqrt(sigma2^2 - zc^2) )/zc;

%%
% Compute eigs of P near wc1. The eigs of P near wc2 have large residual
% with respect to the original problem, so we don't want those.
tic
w_eigs1 = eigs(A_companion, B_companion, 10, wc1);
fprintf('Took %.2f s to find 10 eigs of P near %f\n', toc, wc1);
fprintf('That is, 10 eigs of T near %f\n', zc);
z_eigs1 = 2*sigma2./(w_eigs1 + 1./w_eigs1);

% fprintf('     sqrt(lambda)        residual   Q (want > 10)\n');
% fprintf('--------------------------------------------------\n');
% for z = z_eigs1.'
%     L = z^2;
%     resid = svds(F(L),1,0);
%     Q = real(z)/imag(z)/2;
%     fprintf('%+6.2e + %+6.2ei   %4.2e   %4.2e\n', real(z),imag(z),resid,Q);
% end
% 
% % These aren't the ones we want--large residual.
% tic
% w_eigs2 = eigs(A_companion, B_companion, 10, wc2);
% fprintf('Took %.2f s to find 10 eigs of P near wc2\n', toc);
% z_eigs2 = 2*sigma2./(w_eigs2 + 1./w_eigs2);
% for z = z_eigs2.'
%     L = z^2;
%     resid = svds(F(L),1,0);
%     fprintf('%+6.2e + %+6.2ei, %4.2e\n', real(z),imag(z),resid);
% end

%%
% Compute many more eigs, plot, and save.
if exist('gun_eigs.mat','file') == 2
    load('gun_eigs','z_eigs');
else
    N_eigs = 100;
    tic
    w_eigs = eigs(A_companion, B_companion, N_eigs, wc1);
    fprintf('Took %.2f s to find %d eigs of T near %f\n', toc, N_eigs, zc);
    z_eigs = 2*sigma2./(w_eigs + 1./w_eigs);
    save('gun_eigs','z_eigs');
end

figure, hold all
plot(real(z_eigs ),imag(z_eigs ),'*c');
plot(real(z_eigs1),imag(z_eigs1),'*k');
plot(real(zc     ),imag(zc     ),'*r');

end
