function error = EE_err(kgrid, coefs, zgrid, transmat, fspace)

    nk = length(kgrid);
    nz = length(zgrid);

    error = zeros(nk*nz,1);

    for iz =1:nz
        z = exp(zgrid(iz));
        znext = exp(zgrid)';

        % C(k,z), l(k,z)
        l = funeval(coefs(:,iz),fspace,kgrid);
        C = 1./(l.^1.33) *.67*z.*kgrid.^.33;

        % kp(k,z)
        kp = z*kgrid.^.33.*l.^.67 + .9*kgrid - C;

        lnext = zeros(nk,nz);
        for jz = 1:nz
            lnext(:,jz) = funeval(coefs(:,jz),fspace,kp);
        end
        
        Cnext = 1./(lnext.^1.33) *.67.*znext.*kp.^.33;
        MRnext = .97 * C./Cnext .* (.33*znext.*kp.^(-.67).*lnext.^.67 + .9);
        error((iz-1)*nk+1:iz*nk) = MRnext * transmat(iz, :)' - 1;
    end
end