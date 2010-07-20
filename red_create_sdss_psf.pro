FUNCTION red_create_sdss_psf, path, ts_file, pos, box
;red_create_sdss_psf, '/home/barden/redshifting_galaxies/sample/data/', 'tsField-001896-2-41-0065.fit', [459, 269], 20
  percent = 0.2

  red_string2sdss, ts_file, run, camcol, rerun, field

  filters = ['u', 'g', 'r', 'i', 'z']
  
  ctr = box/2

  FOR i = 0, 4 DO BEGIN
    fpc = red_sdss2string(run, camcol, rerun, field, fpc = filters[i])
    fits_read, path+fpc, im, hd

    sz = [sxpar(hd, 'NAXIS1'), sxpar(hd, 'NAXIS2')]
    x = [pos[0]-1-box/2 > 0 < (sz[0]-1), pos[0]-1+box/2 > 0 < (sz[0]-1)]
    y = [pos[1]-1-box/2 > 0 < (sz[1]-1), pos[1]-1+box/2 > 0 < (sz[1]-1)]

    dum = gauss2dfit(im[x[0]:x[1], y[0]:y[1]], par)
;    print, par[2:5]

    out = sshift2d(im[x[0]:x[1], y[0]:y[1]], ctr-par[4:5])
    dum = gauss2dfit(out, par)
;    print, par[2:5]

    sz = [1, 1]*box
    IF NOT (box MOD 2) THEN sz += 1
    dist_circle, arr, sz, ctr, ctr

    ord = sort(arr)
    n = n_elements(arr)
    nout = fix(n*percent)
    idx = ord[n-1-nout:n-1]
    resistant_mean, out[idx], 3, m
    out -= m
    out /= total(out)

;    plot, arr, out, psym = 3
;    hor, 0, col = 150
;    ver, min(arr[idx])
;    ds9_array, out

    IF i EQ 0 THEN psf = out ELSE psf = [[[psf]], [[out]]]
  ENDFOR

  return, psf
END
