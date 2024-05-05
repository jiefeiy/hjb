clear, clc
% numerical tests:
% [1] Forsyth, Labahn. (2007) Numerical methods for controlled
% Hamilton-Jacobi-Bellman PDEs in finance. 
%
% Example 2.3: American option 
% \min( -V_t - {\sigma^2 S^2/2 V_{SS} + rS V_S - rV}, V - V^* ) = 0.
% 
% This is equivalent to a penalized control problem 
% V_t + 
% \sup_{\mu\in\{0, 1\}} {\sigma^2 S^2/2 V_{SS} + rS V_S - rV + \mu(V^* - V)/\varepsilon} = 0.
%
% Let L^Q V := a(S, \tau, Q) V_SS + b(S, \tau, Q)V_s - c(S, \tau, Q)V.
% The HJB has the form of 
% V_\tau = sup_Q {L^Q V + d(S, \tau, Q)}. 
% 
% - a(S, \tau, Q) = \sigma^2 S^2/2,
% - b(S, \tau, Q) = rS, 
% - c(S, \tau, Q) = r + \mu/\varepsilon
% - d(S, \tau, Q) = \mu/varepsilon
%
% Aim to solve the problem on [0,T]\times[0, S_{max}]










