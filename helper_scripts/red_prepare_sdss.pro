@red_photometry.pro

PRO red_prepare_sdss, path, ts_file, ra, dec, im, imerr, post, psfpos, $
                      ctrpos, tmpdir, box, sky, pre=pre, sdss_psf=sdss_psf
;Prepare the SDSS images for all 5 filters.
;For a given PATH and TS-FILE, and source position (RA,DEC), calculate
;the extinction corrected image (IM) and error-image (IMERR).
;The output variables IM and IMERR are in cts. IM has softbias and sky
;(from header) removed

  red_sdss_params, exptimelo, scllo, dum0, path, ts_file, dum1, dum2, dum3, $
                   gainlo, darkvarlo, dum6, zplo

  red_string2sdss, ts_file, run, camcol, rerun, field

  fwhm = red_shift_match_psf(path, ts_file, post, psfpos, ctrpos, tmpdir, $
                             box, sdss_psf=sdss_psf)

;store images in one 3d array: im[*,*,0:4]: u,g,r,i,z-band
  filter = ['u', 'g', 'r', 'i', 'z']
  IF keyword_set(pre) THEN $
     print, strrep(strrep(red_sdss2string(run, camcol, rerun, field, $
                                          fpc = filter[0]), $
                          '.fit', post+'.fit'), 'fpC', pre) $
  ELSE pre = 'fpC'
  fits_read, path+strrep(strrep(red_sdss2string(run, camcol, rerun, field, $
                                                fpc = filter[0]), $
                                '.fit', post+'.fit'), 'fpC', pre), im
  im = float(im)
  FOR i = 1, 4 DO BEGIN
    fits_read, path+strrep(strrep(red_sdss2string(run, camcol, rerun, field, $
                                                  fpc = filter[i]), $
                                  '.fit', post+'.fit'), 'fpC', pre), im0
    im = [[[im]], [[im0]]]
  ENDFOR

;the value for the sky is in the image header
  tsf = readfits(path+ts_file, tsfhd, /exten, /silent)
;sky and skyerr
  sky = maggies2cts(tbget(tsfhd, tsf, 'sky'), exptimelo, zplo)*scllo^2
  skyerr = maggies2cts(tbget(tsfhd, tsf, 'skyerr'), exptimelo, zplo)*scllo^2

;calculate the error of the counts-image
;error(counts) = sqrt([counts+sky]/gain + Npix*(dark_variance+skyErr))
  imerr = im
  FOR i = 0, 4 DO $
     imerr[*, *, i] = sqrt((im[*, *, i]-1000.)/gainlo[i]+ $
                           (darkvarlo[i]+skyerr[i]))

;apply exctinction correction
;  tso = readfits(path+red_sdss2string(run, camcol, rerun, field, /tso), $
;                 tsohd, /exten, /silent)
;  extinction = tbget(tsohd, tso, 'reddening')
;  IF ra LT 0 THEN BEGIN
;    s = where(tbget(tsohd, tso, 'id') EQ dec)
;    ra = tbget(tsohd, tso, 'ra') & ra = ra[s]
;    dec = tbget(tsohd, tso, 'dec') & dec = dec[s]
;  ENDIF ELSE BEGIN
;    tso_ra = tbget(tsohd, tso, 'ra')
;    tso_de = tbget(tsohd, tso, 'dec')
;    srccor, tso_ra/15., tso_de, ra/15., dec, 50, s, d, $
;            /spherical, option = 1, /silent
;    IF n_elements(d) GT 1 OR d[0] EQ -1 THEN $
;      message, 'Something went wrong in the cross-correlation!'
;  ENDELSE
;  extinction = reform(extinction[*, s], 5)

;==============================================================================
;applying an extinction correction results in a systematic bias as a
;function of redshift
;==============================================================================
;extinction is in luptitudes, image is in cts!!!
;convert image to luptitudes, apply extinction correction, then go
;back to counts
;  filt = ['u', 'g', 'r', 'i', 'z']
;  FOR i = 0, 4 DO BEGIN
;    im[*, *, i] -= (sky[i]+1000.)
;    im[*, *, i] = cts2maggies(im[*, *, i], exptimelo, zplo[i])
;    im[*, *, i] = maggies2lup(im[*, *, i], filt[i])
;    im[*, *, i] -= extinction[i]
;    im[*, *, i] = lup2maggies(im[*, *, i], filt[i])
;    im[*, *, i] = maggies2cts(im[*, *, i], exptimelo, zplo[i])
;  ENDFOR
;==============================================================================

  delta_m = [-0.036, 0.012, 0.010, 0.028, 0.040]
;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
;get SDSS counts onto AB system and subtract sky
;  sky = [1038.7558, 1085.2558, 1134.156, 1172.496, 1120.596];save
  sky = fltarr(5)
  FOR i = 0, 4 DO BEGIN
;    sky[i] = red_find_sky(im[*, *, i], 50, /twosided, sigma = [-3, 3])
     sky[i] = ring_sky(im[*, *, i], 20, 5)
    im[*, *, i] -= sky[i]
    im[*, *, i] = cts2maggies(im[*, *, i], exptimelo, zplo[i])
    im[*, *, i] *= 10^(-0.4*(delta_m[i]));4 conversion done!!!!!!!!!!!!!!!
    im[*, *, i] = maggies2cts(im[*, *, i], exptimelo, zplo[i])
  ENDFOR
;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

;images are now in cts and sky subtracted
END
