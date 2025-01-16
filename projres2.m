function res = projres2(coefs, T, fspace, knodes, zgrid, transmat)

    nz = length(zgrid);
    nk = length(knodes);
    coefs = reshape(coefs,nk,nz);

    res = zeros(nk*nz, 1);

    for iz =1:nz
        z = exp(zgrid(iz));
        znext = exp(zgrid)';

        % C(k,z), l(k,z)
        l = T*coefs(:,iz);
        C = 1./(l.^1.33) *.67*z.*knodes.^.33;

        % kp(k,z)
        kp = z*knodes.^.33.*l.^.67 + .9*knodes - C;

        lnext = zeros(nk,nz);
        for jz = 1:nz
            lnext(:,jz) = funeval(coefs(:,jz),fspace,kp);
        end
        
        Cnext = 1./(lnext.^1.33) *.67.*znext.*kp.^.33;
        MRnext = .97 * C./Cnext .* (.33*znext.*kp.^(-.67).*lnext.^.67 + .9);
        res((iz-1)*nk+1:iz*nk) = MRnext * transmat(iz, :)' - 1;
    end

end