;+
;
; NAME: rescale_ivar
;
; PURPOSE: Rescale inverse variance based on input sky
;  residual RMS and flux - synflux residuals
;
; USAGE:
;   newivar = rescale_ivar(resflux=resflux, invvar=invvar, $
;    loglam=loglam, skyloglam=skyloglam, skyrenorm=skyrenorm)
;
; ARGUMENTS: (all keywords, all mandatory)
;   resflux: npix X nspec array of flux - synflux residuals
;   invvar: initial inverse-variance array (also npix X nspec)
;   loglam: npix vector of log-Angstrom wavelengths for resflux
;   skyrenorm: factor by which to rescale error estimates based
;     sky-subtraction residual analysis.
;   skyloglam: log-Angstrom wavelength baseline for skyrenorm
;
; RETURNS: a new inverse-variance array.
;
; Copyright 2008 Adam S. Bolton
;
;-

function rescale_ivar, resflux=resflux, invvar=invvar, loglam=loglam, $
 skyloglam=skyloglam, skyrenorm=skyrenorm

npix = (size(resflux))[1]
nspec = (size(resflux, /n_dimen) eq 2) ? (size(resflux))[2] : 1L

; First rescale on the basis of sky residuals:

thisnorm = interpol(skyrenorm, skyloglam, loglam)

outivar = float(invvar * ((1./thisnorm^2) # replicate(1., nspec)))

; Now rescale on the basis of spectrum residuals:

specscale = replicate(1., nspec)
minpix = 100.
twosiglo = 0.025
twosighi = 0.975
for i = 0L, nspec-1 do begin
    resthis = resflux[*,i]
    ivarthis = outivar[*,i]
    whgood = where(ivarthis gt 0., ngood)
    if (ngood ge minpix) then begin
        resthis = resthis[whgood]
        ivarthis = ivarthis[whgood]
        resthis = resthis * sqrt(ivarthis)
        resthis = resthis[sort(resthis)]
        resrank = findgen(ngood) / float(ngood-1)
        res_2siglo = interpol(resthis, resrank, twosiglo)
        res_2sighi = interpol(resthis, resrank, twosighi)
        specscale[i] = ((res_2sighi - res_2siglo) / 4.) > 1.
    endif
endfor

outivar = outivar * (replicate(1., npix) # (1./specscale^2))

return, outivar
end

