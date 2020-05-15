;+
;
; NAME: galtempfit
;
; PURPOSE: perform fit of rest-frame template set to observed galaxy
;   spectrum
;
; USAGE:
;   yfit = galtempfit(objflux=objflux, objivar=objivar, objloglam=objloglam, z=z, $
;                     tempset=tempset, temploglam=temploglam, npoly=npoly)
;
; ARGUMENTS (all keywords, all mandatory):
;   objflux: object flux vector
;   objivar: object inverse-variance vector
;   objloglam: object log10-wavelength vector
;   z: galaxy redshift
;   tempset: ntpix X ntemp array of eigen-galaxy spectra
;   temploglam: template log10-wavelength (rest-frame) baseline vector
;   npoly: number of polynomial terms to incorporate in addition
;
; RETURNS:
;   Fitted model to the spectrum (that's all for now).
;
; Copyright 2008 Adam S. Bolton
;
; Added zeroing-out of non-covered red-end rest wavelengths
; via tempmask (also output), ASBfeb2010
;
;-

function galtempfit, objflux=objflux, objivar=objivar, objloglam=objloglam, z=z, $
 tempset=tempset, temploglam=temploglam, npoly=npoly, tempmask=tempmask, coeff=coeff

ntemp = (size(tempset))[2]
npix = n_elements(objloglam)

; Interpolate the templates to the frame of the galaxy:

thistemp = fltarr(npix, ntemp)
restloglam = objloglam - alog10(1. + z)
maxrest = max(temploglam)
for i = 0L, ntemp-1 do thistemp[*,i] = $
 interpol(tempset[*,i], temploglam, restloglam)

; Normalize the templates:
for i = 0L, ntemp-1 do thistemp[*,i] = thistemp[*,i] / sqrt(mean(thistemp[*,i]^2))

; Add on some polynomial terms:
pbase = findgen(npix) / float(npix-1)  ;[0.0,......,1.0] array
polyset = fpoly(pbase, npoly)
;pbase = [0.1,0.2,0.3,0.4]
;IDL> print,fpoly(pbase, 3)
;      1.00000      1.00000      1.00000      1.00000  --x^0
;     0.100000     0.200000     0.300000     0.400000  --x^1
;    0.0100000    0.0400000    0.0900000     0.160000  --x^2

thistemp = [[thistemp], [polyset]] ;stack along the axis-0 (in python axis order convention)
;shape [npix,ntemp+npoly]

; Do the fitting:
tempmask = (restloglam le maxrest)  ;[1,1,1,1,....,0,0] or [1,1,1,1,...,1]
ithistemp = thistemp * ((objivar * tempmask) # replicate(1., ntemp+npoly))  
;[npix,ntemp+npoly]*[npix,ntemp+npoly] = [npix,ntemp+npoly]
;IDL> print,x
;       1       2       3
;IDL> print,y
;      1.00000      1.00000
;IDL> print,x#y
;      1.00000      2.00000      3.00000
;      1.00000      2.00000      3.00000

alpha = transpose(ithistemp) # thistemp ; [ntemp+npoly,npix]#[npix,ntemp+npoly] --> [ntemp+npoly,ntemp+npoly]
beta = transpose(ithistemp) # objflux  ;[ntemp+npoly,npix]#[npix] -> [ntemp+npoly,npix]#[npix，1] -> [ntemp+npoly,1] -> [ntemp+npoly]
;IDL> print,z
;      1.00000      2.00000      3.00000
;      1.00000      2.00000      3.00000
;IDL> print,z#[1,2]
;      3.00000      6.00000      9.00000

;https://www.harrisgeospatial.com/docs/SVDC.html
;https://www.harrisgeospatial.com/docs/SVSOL.html
svdc, alpha, w, u, v ;do svd for alpha matrix, w,u,v are decomposed matrix
coeff = svsol(u, w, v, beta) ;use results of svd, to solve the matrix equation alpha*coeff = beta
;coeff [ntemp+npoly] -> 7 in this case

;ialpha = invert(alpha)
;coeff = ialpha # beta
yfit = thistemp # coeff  ;[npix,ntemp+npoly]#[ntemp+npoly] -> [npix,1] -> [npix]
;stop
return, yfit
end
