FUNCTION red_read_sdss_psf, psf_file, tmpdir, psfpos, lambda, z, out = out
;psf_file: SDSS psf file, e.g. psField-002243-3-0284.fit (including path)
;tmpdir: working directory for temporary psf storage of "psf.fits"
;(not deleted upon exit)
;psfpos: pixel position where the PSF should be calculated ([x,y]-array)
;==============================================================================
;lambda: wavelength in angstroem of the output filter
;z: redshift of the current source
;==============================================================================
;if z is not given lambda is assumed to be an integer specifying the
;filter: 0,1,2,3,4 for u,g,r,i,z
;==============================================================================
;out: location of output psf (including path)

  print, 'Reading SDSS PSF...'

  IF n_params() EQ 5 THEN BEGIN
    delta = abs(red_lambda_SDSS()-lambda/(1.+z))
    psf_filter = where(delta EQ min(delta))+1
  ENDIF ELSE psf_filter = lambda+1

  read_PSF_path = '/home/barden/Apps/photoop-v1_0/'

;  psf_filter = 5
;  print, psf_filter
;  print, round2(xpsf, 1, /string)+' '+round2(ypsf, 1, /string)

  spawn, read_PSF_path+'read_PSF ' +psf_file+' '+ $
         round2(psf_filter, 0, /string)+' '+ $
         round2(psfpos[0], 1, /string)+' '+ $
         round2(psfpos[1], 1, /string)+' '+tmpdir+'psf.fits'
  fits_read, tmpdir+'psf.fits', psf

  psf = sdss_psf_recon(mrdfits(psf_file, psf_filter), psfpos[0], psfpos[1])

  psf = psf[10:40, 10:40]
;  psf -= 1000.
  psf /= total(psf)
  fits_write, tmpdir+'psf_org.fits', psf

  str = ''
  openr, 1, '/home/barden/IDL/red_read_sdss_psf.obj'
  openw, 2, tmpdir+'psf.obj'
  ct = 0
  WHILE NOT eof(1) DO BEGIN
    readf, 1, str
    IF ct EQ 0 THEN printf, 2, 'A) '+tmpdir+'psf.fits'
    IF ct EQ 1 THEN printf, 2, 'B) '+tmpdir+'psf.gf.fits'
    IF ct GT 1 THEN printf, 2, str
    ct++
  ENDWHILE
  close, 1
  close, 2

  fits_write, tmpdir+'psf.fits', psf

  spawn, '/home/barden/galfit/galfit3 '+tmpdir+'psf.obj', dum
  res = red_read_sersic_results(tmpdir+'psf.gf.fits')

  openr, 1, '/home/barden/IDL/red_read_sdss_psf.obj'
  openw, 2, tmpdir+'psf.obj'
  ct = 0
  WHILE NOT eof(1) DO BEGIN
    readf, 1, str
    IF ct EQ 22 THEN BREAK
    printf, 2, str
    ct++
  ENDWHILE
  printf, 2, '0) sersic'
  printf, 2, '1) '+round2(res[10], 3, /str)+' '+round2(res[12], 3, /str)+' 1 1'
  printf, 2, '3) '+round2(res[0], 2, /str)+' 1'
  printf, 2, '4) '+round2(res[2], 2, /str)+' 1'
  printf, 2, '5) '+round2(res[4], 2, /str)+' 1'
  printf, 2, '9) '+round2(res[6], 2, /str)+' 1'
  printf, 2, '10) '+round2(res[8], 2, /str)+' 1'
  printf, 2, 'Z) 0'
  printf, 2
  printf, 2, '0) sersic'
  printf, 2, '1) '+round2(res[10], 3, /str)+' '+round2(res[12], 3, /str)+' 1 1'
  printf, 2, '3) '+round2(res[0]+2, 2, /str)+' 1'
  printf, 2, '4) '+round2(res[2], 2, /str)+' 1'
  printf, 2, '5) '+round2(res[4], 2, /str)+' 0'
  printf, 2, '9) '+round2(res[6], 2, /str)+' 1'
  printf, 2, '10) '+round2(res[8], 2, /str)+' 1'
  printf, 2, 'Z) 0'
  close, 1
  close, 2

  spawn, '/home/barden/galfit/galfit3 '+tmpdir+'psf.obj', dum

;  fits_read, tmpdir+'psf.gf.fits', psf, hd, exten = 3
;  ds9_array, psf
;  stop
  fits_read, tmpdir+'psf.gf.fits', psf, hd, exten = 2

  nx1 = sxpar(hd, 'NAXIS1')
  nx2 = sxpar(hd, 'NAXIS2')
  IF NOT (nx1 MOD 2) THEN BEGIN
    nx1++
    new = fltarr(nx1, nx2)
    new[1:*, *] = psf
    psf = new
  ENDIF
  IF NOT (nx2 MOD 2) THEN BEGIN
    nx2++
    new = fltarr(nx1, nx2)
    new[*, 1:*] = psf
    psf = new
  ENDIF

  dum = gauss2dfit(psf, par)
  psf = sshift2d(psf, 15-par[4:5])

  psf /= total(psf)

  IF keyword_set(out) THEN fits_write, out, psf

;  CASE psf_filter OF
;    1: f = 'u'
;    2: f = 'g'
;    3: f = 'r'
;    4: f = 'i'
;    5: f = 'z'
;  ENDCASE
;  fits_write, '/home/barden/redshifting_galaxies/sample/sim/psf_000752_'+ $
;              f+'3_0511.fits', psf
;stop

  print, 'done!'

  return, psf
END
