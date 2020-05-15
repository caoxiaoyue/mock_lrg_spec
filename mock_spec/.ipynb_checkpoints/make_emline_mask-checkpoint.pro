;+
;
; NAME: make_emline_mask
;
; PURPOSE: make byte mask for rest-frame emission lines
;
; USAGE:
;   emask = make_emline_mask(wave=wave, z=z, linelist=linelist, hw=hw)
;
; ARGUMENTS (all keywords, all mandatory):
;   loglam: log10-wavelength baseline vector (MUST BE MONOTONIC!)
;   z: vector of redshifts
;   linelist: list of emission log10-wavelengths to mask
;   hw: half-width of masking region, in pixels
;
; RETURNS:
;   an nwave X nz bytarr, =1 where not emission-line masked
;
; Copyright 2008 Adam S. Bolton
;
;-

function make_emline_mask, loglam=loglam, z=z, linelist=linelist, hw=hw

npix = n_elements(loglam)
nspec = n_elements(z)
nline = n_elements(linelist)

; Make a buffered log-wavelength vector to catch overflow:
wstep_lo = loglam[1] - loglam[0]
wstep_hi = loglam[npix-1] - loglam[npix-2]
wbuff_lo = loglam[0] - wstep_lo * reverse(findgen(hw+1L) + 1.)
wbuff_hi = loglam[npix-1] + wstep_hi * (findgen(hw+1L) + 1.)

floglam = [wbuff_lo, loglam, wbuff_hi]
nfpix = n_elements(floglam)

emask = replicate(1B, nfpix, nspec)

; Loop over spectra and lines (I know, it's inelegant...)
for i = 0L, nline-1 do begin & $
    for j = 0L, nspec-1 do begin & $
        loglamobs = linelist[i] + alog10(1. + z[j]) & $
        junk = min(abs(floglam - loglamobs), whmin) & $
        pix_lo = (whmin - hw) > 0 & $
        pix_hi = (whmin + hw) < (nfpix - 1) & $
        emask[pix_lo:pix_hi,j] = 0B & $
    endfor & $
endfor

emask = emask[hw+1:hw+npix,*]

return, emask
end




