; Get sky-residual rescaling info:
in_file = FILEPATH('boss_skyrenorm_' + getenv('RUN2D') + '.fits', ROOT_DIR='../data/')
skystruc = mrdfits(in_file,1)

; Get the plates from the skystruc:
plate = skystruc.platelist
mjd = skystruc.mjdlist

nmain = n_elements(plate)

splog, 'Found ', nmain, ' good plates'

; Get the emission lines:
readcol, '../emlines.txt', emwave, emname, format='F,A'

; Get the empirical absorption-line info:
readcol, '../abslines.txt', abswave, absname, format='F,A', comment='#'

; Get empirical skyline info:
readcol, '../skylines.txt', skywave, skyname, format='F,A', comment='#'

; Use Yanmei's templates:
gtfile = getenv('IDLSPEC2D_DIR') + '/templates/spEigenGal-55740.fits'
;gtfile = 'spMLpcaGal-55331-7temp.fits'

;-------------------------------------------
galtemp = mrdfits(gtfile,0,galhdr)
ntpix = (size(galtemp))[1]
ntemp = (size(galtemp))[2]
gtloglam = sxpar(galhdr, 'COEFF0') + sxpar(galhdr, 'COEFF1') * findgen(ntpix)

num_lrgs = 0
; Loop over plates:
for i = 0L, nmain-1 do begin
    filename_h5 = STRCOMPRESS(string(plate[i])+'_'+string(mjd[i])+'.h5',/remove_all)
    out_file = FILEPATH(filename_h5, ROOT_DIR='./sky_data/')
    
    ;if (plate[i] ne 6759) or (mjd[i] ne 56416) then begin
    ;    continue
    ;endif
    
    splog, 'plate number ===', plate[i], 'mjd number =',mjd[i]
    splog, 'Scanning plate ', i+1, ' of ', nmain

    readspec, plate[i], mjd=mjd[i], flux=flux, $
      synflux=synflux, $
      invvar=invvar, plug=plug, wave=wave, loglam=loglam, $
      zans=zans, zhdr=zhdr, $
      andmask=andmask, ormask=ormask, /align


      newivar = skymask(invvar, andmask)

; Use sky fibers instead, if requested:
    whgal = where(strtrim(zans.objtype,2) eq 'SKY', ngal)

    ;print,'----ngal-----',ngal
    if (ngal gt 0) then begin

       ; Structure to hold info about the plates under consideration:
        fibermask = keyword_set(sdss1) ? replicate(0B, 640) : replicate(0B, 1000)
        fibermask[whgal] = 1B
        plate_struc = {plate: plate[i], mjd: mjd[i], fibermask: fibermask, $
                       wave: 10.^skystruc.skyloglam, med_lnoise: 0. * skystruc.skyloglam}

        flux = flux[*,whgal]
        synflux = synflux[*,whgal]
        invvar = invvar[*,whgal]
        newivar = newivar[*,whgal]
        andmask = andmask[*,whgal]
        ormask = ormask[*,whgal]
        zans = zans[whgal]
        plug = plug[whgal]

        zans.z_noqso *= 0.

; Generate emission/absorption-line mask:
        emask = make_emline_mask(loglam=loglam, z=zans.z_noqso, linelist=alog10([emwave, abswave]), hw=10)

; Kluge observed-frame sky feature masks onto this (zeroing out redshifts):
        emask = emask * make_emline_mask(loglam=loglam, z=0.*zans.z_noqso, linelist=alog10(skywave), hw=10)
        
        
        ;write hdf5
        mg_h5_putdata, out_file, 'flux', flux
        mg_h5_putdata, out_file, 'loglam', loglam
        mg_h5_putdata, out_file, 'zans', zans

        print,out_file+' hdf5 file written finished!'
        
        num_lrgs += (size(flux))[2]
        if (num_lrgs gt 10000) then break

    endif
endfor

end

