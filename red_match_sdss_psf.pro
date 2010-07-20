FUNCTION symmetrise_psf, psf
  psfo = psf
  psf_s1 = rotate(psf, 1)
  psf_s2 = rotate(psf, 2)
  psf_s3 = rotate(psf, 3)
  FOR i = 0, n_elements(psf)-1 DO BEGIN
    resistant_mean, [psf[i], psf_s1[i], psf_s2[i], psf_s3[i]], 3, mn
    psfo[i] = mn
  ENDFOR
  return, psfo
END

FUNCTION amoeba_func, sig
  COMMON amoeba_func_par, amo_x, amo_y, amo_narrow, amo_wide

  z = sig[1]+sig[2]*exp(-0.5*(((amo_x-15)/sig[0])^2+((amo_y-15)/sig[0])^2))
  c = convolve(amo_narrow, z/total(z))
  return, total(abs(amo_wide-c))
END

PRO red_match_sdss_psf, path, ts_file, post, psfpos, tmpdir
  COMMON amoeba_func_par, amo_x, amo_y, amo_narrow, amo_wide
;red_match_sdss_psf, '/home/barden/redshifting_galaxies/sample/data/', 'tsField-001896-2-41-0065.fit', '.sm', [1551.2723, 1031.1382], '/home/barden/redshifting_galaxies/sample/data/tmp_batch3/'

  file_mkdir, tmpdir

  red_string2sdss, ts_file, run, camcol, rerun, field
  psfield = red_sdss2string(run, camcol, rerun, field, /psf)

  filters = ['u', 'g', 'r', 'i', 'z']

  psf_tab = readfits(path+ts_file, psf_tab_hd, /exten, /silent)
  psf_width = tbget(psf_tab_hd, psf_tab, 'psf_width')

  loidx = (where(psf_width EQ max(psf_width)))[0]
  psflo = red_read_sdss_psf(path+psfield, tmpdir, psfpos+1, loidx)

  fpc = red_sdss2string(run, camcol, rerun, field, fpc = filters[loidx])
  file_copy, path+fpc, path+strrep(fpc, '.fit', post+'.fit'), /overwrite

  FOR i = 0, 4 DO BEGIN
    IF i EQ loidx THEN CONTINUE
    psfhi = red_read_sdss_psf(path+psfield, tmpdir, psfpos+1, i)

    nx = 31 & ny = 31
    amo_x = findgen(nx)#replicate(1, ny)
    amo_y = replicate(1, nx)#findgen(ny)
    amo_narrow = symmetrise_psf(psfhi)
    amo_wide = symmetrise_psf(psflo)
    sig = amoeba(1.0e-5, scale = [1e-3, 1e-6, 1e-3], p0 = [1, 0, 1], $
                 FUNCTION_name = 'amoeba_func')

    fpc = red_sdss2string(run, camcol, rerun, field, fpc = filters[i])
    fits_read, path+fpc, in, hd

    z = exp(-0.5*(((amo_x-15)/sig[0])^2+((amo_y-15)/sig[0])^2))
    c = convolve(amo_narrow, z/total(z))
    print, total(abs(amo_wide-c)), sig

    out = convolve(in, z/total(z))
    fits_write, path+strrep(fpc, '.fit', post+'.fit'), out, hd

    fits_write, path+strrep(strrep(fpc, '.fit', post+'.fit'), 'fpC', 'psf'), $
                convolve(psfhi, z/total(z))
  ENDFOR
END

PRO red_shift_match_psf, path, ts_file, post, psfpos, tmpdir, pos, box, $
                         symmetrise = symmetrise
  COMMON amoeba_func_par, amo_x, amo_y, amo_narrow, amo_wide
;red_shift_match_psf, '/home/barden/redshifting_galaxies/sample/data/', 'tsField-001896-2-41-0065.fit', '.sm', [1551.2723, 1031.1382], '/home/barden/redshifting_galaxies/sample/data/tmp_batch3/', [533, 444], 30

  file_mkdir, tmpdir

  red_string2sdss, ts_file, run, camcol, rerun, field
  psfield = red_sdss2string(run, camcol, rerun, field, /psf)

  filters = ['u', 'g', 'r', 'i', 'z']

;find the band with the worst seeing
  psf_tab = readfits(path+ts_file, psf_tab_hd, /exten, /silent)
  psf_width = tbget(psf_tab_hd, psf_tab, 'psf_width')
  loidx = (where(psf_width EQ max(psf_width)))[0]

;read the worst PSF and corresponding image
  psflo = red_read_sdss_psf(path+psfield, tmpdir, psfpos+1, loidx)
  fpc = red_sdss2string(run, camcol, rerun, field, fpc = filters[loidx])
  fits_read, path+fpc, im, hd

;make slightly smoothed copies
  smo = 1.
  nx = 31 & ny = 31
  amo_x = findgen(nx)#replicate(1, ny)
  amo_y = replicate(1, nx)#findgen(ny)
  z = exp(-0.5*(((amo_x-15)/smo)^2+((amo_y-15)/smo)^2))
;symmetrise the input PSF?
  IF keyword_set(symmetrise) THEN psflo = symmetrise(psflo)
;convolve and save
  psflo = convolve(psflo, z/total(z))
  im = convolve(im, z/total(z))
  fits_write, path+strrep(strrep(fpc, '.fit', post+'.fit'), 'fpC', 'psf'), $
              psflo
  fits_write, path+strrep(fpc, '.fit', post+'.fit'), im, hd

;fit the reference position with a 2d Gaussian
  sz = [sxpar(hd, 'NAXIS1'), sxpar(hd, 'NAXIS2')]
  x = [pos[0]-1-box/2 > 0 < (sz[0]-1), pos[0]-1+box/2 > 0 < (sz[0]-1)]
  y = [pos[1]-1-box/2 > 0 < (sz[1]-1), pos[1]-1+box/2 > 0 < (sz[1]-1)]
  flo = gauss2dfit(im[x[0]:x[1], y[0]:y[1]], plo)
  plo = plo[4:5]

  print, plo

;loop over all other bands
  FOR i = 0, 4 DO BEGIN
    IF i EQ loidx THEN CONTINUE
;read PSF
    psfhi = red_read_sdss_psf(path+psfield, tmpdir, psfpos+1, i)

;symmetrise the input PSF?
    IF keyword_set(symmetrise) THEN BEGIN
      amo_narrow = symmetrise_psf(psfhi)
      amo_wide = symmetrise_psf(psflo)
    ENDIF ELSE BEGIN
      amo_narrow = psfhi
      amo_wide = psflo
    ENDELSE

;fit the sigma of the narrow PSF
    sig = amoeba(1.0e-5, scale = [1e-3, 1e-6, 1e-3], p0 = [1, 0, 1], $
                 FUNCTION_name = 'amoeba_func')

;read the image
    fpc = red_sdss2string(run, camcol, rerun, field, fpc = filters[i])
    fits_read, path+fpc, in, hd

;fit the position
    fhi = gauss2dfit(in[x[0]:x[1], y[0]:y[1]], phi)
    phi = plo-phi[4:5]

    z = exp(-0.5*(((amo_x-15-phi[0])/sig[0])^2+((amo_y-15-phi[1])/sig[0])^2))
    c = convolve(amo_narrow, z/total(z))
;    print, total(abs(amo_wide-c)), sig

    out = convolve(in, z/total(z))
    fits_write, path+strrep(fpc, '.fit', post+'.fit'), out, hd

    fhi = gauss2dfit(out[x[0]:x[1], y[0]:y[1]], phi)
    phi = plo-phi[4:5]
    print, sig[0], sqrt(total(phi^2)), phi

    z = exp(-0.5*(((amo_x-15)/sig[0])^2+((amo_y-15)/sig[0])^2))

    fits_write, path+strrep(strrep(fpc, '.fit', post+'.fit'), 'fpC', 'psf'), $
                convolve(psfhi, z/total(z))
  ENDFOR
END
