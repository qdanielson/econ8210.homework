%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Housekeeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% add the Miranda-Fackler toolbox
addpath('../compecon/CEtools');

% add the Elenev et al. toolbox
addpath('../');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rrho = .95;
ssigma = .007;

% Rouwenhorst method
nz = 3;
mu_uncond = -.5 * ssigma^2/(1-rrho^2);
[transmat, zgrid] = Helpers.rouwen(rrho, mu_uncond, ssigma, 0, nz);
transmat = transmat';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Steady state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [1, 0.1]';

xsol = fsolve(@ssres, x0);
[~, css] = ssres(xsol);
kss = xsol(1);
lss = xsol(2);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Projection with Chebyshev
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define grids
nk = 6;
kmin = 0.01*kss;
kmax = 2*kss;

knodes = chebnode(nk,kmin,kmax);
fspace = fundef({'cheb',nk,kmin,kmax});
T = chebbas(nk,kmin,kmax,knodes,0);

% initial guess
coefs0 = funfitxy(fspace,knodes,lss*ones(nk,1));
coefs0 = repmat(coefs0,1,nz);
coefs0 = reshape(coefs0,nk*nz,1);

f2solve = @(xx) projres(xx, T, fspace, knodes, zgrid, transmat);
opts = optimoptions('fsolve','Display','iter');
coefsopt = fsolve(f2solve, coefs0, opts);

coefsopt = reshape(coefsopt,nk,nz);


% unpack the solved policy function l(k,z)
coefz1 = coefsopt(:,1);
coefz2 = coefsopt(:,2);
coefz3 = coefsopt(:,3);


% plot
kfinegrid = linspace(kmin,kmax,100)';
lz1 = funeval(coefz1,fspace,kfinegrid);
lz2 = funeval(coefz2,fspace,kfinegrid);
lz3 = funeval(coefz3,fspace,kfinegrid);

% get consumption policy
cz1 = 1./(lz1.^1.33) *.67*exp(zgrid(1)).*kfinegrid.^.33;
cz2 = 1./(lz2.^1.33) *.67*exp(zgrid(2)).*kfinegrid.^.33;
cz3 = 1./(lz3.^1.33) *.67*exp(zgrid(3)).*kfinegrid.^.33;


% report mean EE error
error = EE_err(kfinegrid,coefsopt,zgrid,transmat,fspace);
disp(['Mean error:', num2str(mean(error))]);


% plot labor policy
figure;
subplot(2,1,1);
plot(kfinegrid, lz1, 'DisplayName', 'z1');
hold on;
plot(kfinegrid, lz2, 'DisplayName', 'z2');
plot(kfinegrid, lz3, 'DisplayName', 'z3');
title('Labor Policy');
xlabel('Capital (k)');
ylabel('Labor (l)');
legend;
hold off;

% plot consumption policy
subplot(2,1,2);
plot(kfinegrid, cz1, 'DisplayName', 'z1');
hold on;
plot(kfinegrid, cz2, 'DisplayName', 'z2');
plot(kfinegrid, cz3, 'DisplayName', 'z3');
title('Consumption Policy');
xlabel('Capital (k)');
ylabel('Consumption (c)');
legend;
hold off;



% show euler equation error for different number of Chebyshev polynomials
n_values = [6, 10, 20, 50, 100];
mean_errors = zeros(length(n_values), 1);

opts = optimoptions('fsolve','Display','off');

for i = 1:length(n_values)
    nk = n_values(i);
    knodes = chebnode(nk, kmin, kmax);
    fspace = fundef({'cheb', nk, kmin, kmax});
    T = chebbas(nk, kmin, kmax, knodes, 0);

    coefs0 = funfitxy(fspace, knodes, lss * ones(nk, 1));
    coefs0 = repmat(coefs0, 1, nz);
    coefs0 = reshape(coefs0, nk * nz, 1);

    f2solve = @(xx) projres(xx, T, fspace, knodes, zgrid, transmat);
    coefsopt = fsolve(f2solve, coefs0, opts);

    coefsopt = reshape(coefsopt, nk, nz);
    error = EE_err(kfinegrid, coefsopt, zgrid, transmat, fspace);
    mean_errors(i) = mean(error);
end

% create and print table
table_data = table(n_values', mean_errors, 'VariableNames', {'n', 'Mean_EE_Error'});
disp(table_data);