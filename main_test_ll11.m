clear all; close all; clc;

% size_tens = [30 31 22 23];
size_tens = [30 30 30 5];

N = [2 2]; % Order of the tensors
L = [5 5]; % Rank type 1
P = [1 1]; % Rank type 2
U = lp_rnd(size_tens,N,L,P);

T = lpgen(U);

SNR = 10;
Tnoisy = noisy(T,SNR);

%% Let us compare the different versions

options = struct;
options.Display = 1;

% Random version
disp('Random version')
U_random_init = lp_rnd(size_tens,N,L,P);
U_sol_random_init = lp_nls(Tnoisy,U_random_init,N,options);


% GEVD version
disp('GEVD version')
U_gevd_init = ll11_gevd_rank1(Tnoisy,L);
U_sol_gevd_init = lp_nls(Tnoisy,U_gevd_init,N,options);
