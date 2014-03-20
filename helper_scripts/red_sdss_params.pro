PRO red_sdss_params, exptime, scllo, skylo, path, ts_file, aa, kk, airmass, $
                     gain, dark_variance, sky_err, zp, ra, dec, extinction
;PATH: path to the tsField-??????-?-??-????.fit file
;TS_FILE: e.g. 'tsField-001412-2-40-0303.fit'
;ra/dec: position of current object
;the remaining parameters are OUTPUT
;  path = '/home/barden/redshifting_galaxies/best20/'
;  ts_file = 'tsField-001412-2-40-0303.fit'
  exptime = 53.907456
  scllo = 0.396
  skylo = 134.156

  IF n_params() GT 3 THEN BEGIN
    tab = readfits(path+ts_file, hd, exten_no = 1, /silent)
    aa = tbget(hd, tab, 'aa')
    kk = tbget(hd, tab, 'kk')
    airmass = tbget(hd, tab, 'airmass')
    gain = tbget(hd, tab, 'gain')
    dark_variance = tbget(hd, tab, 'dark_variance')
    sky_err = tbget(hd, tab, 'skyerr')/tbget(hd, tab, 'sky')
    zp = -reform(aa+kk*airmass, 5)
  ENDIF

  IF n_params() EQ 15 THEN BEGIN
    obj2 = 'tsObj'+strmid(ts_file, 7, 21)
    tab = readfits(path+obj2, hd, exten_no = 1, /silent)
    extinction = tbget(hd, tab, 'reddening')
    IF ra LT 0 THEN BEGIN
      s = where(tbget(hd, tab, 'id') EQ dec)
      ra = tbget(hd, tab, 'ra') & ra = ra[s]
      dec = tbget(hd, tab, 'dec') & dec = dec[s]
    ENDIF ELSE BEGIN
      obj_ra = tbget(hd, tab, 'ra')
      obj_de = tbget(hd, tab, 'dec')
      srccor, obj_ra/15., obj_de, ra/15., dec, 50, s, d, $
              /spherical, option = 1, /silent
      IF n_elements(d) GT 1 OR d[0] EQ -1 THEN BEGIN
        print, 'Something went wrong in the cross-correlation!'
        stop
      ENDIF
    ENDELSE
    extinction = reform(extinction[*, s], 5)
  ENDIF
END
