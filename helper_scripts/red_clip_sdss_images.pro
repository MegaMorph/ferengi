FUNCTION shift_images, im0, hd0, im1, hd1, out_hd = out_hd
;im0: reference image
;im1: image that is to be shifted
;hd0,hd1: image headers

;  path = '/home/barden/redshifting_galaxies/best20/good/'
;  static_image = 'fpC-003919-r5-0047.fit'
;  movable_image = 'fpC-003919-g5-0047.fit'
;  fits_read, path+static_image, im0, hd0
;  fits_read, path+movable_image, im1, hd1

  sz = (size(im0))[1:2]
;  x0 = sxpar(hd0, 'CRPIX1')
;  y0 = sxpar(hd0, 'CRPIX2')
  x0 = sz[0]/2
  y0 = sz[1]/2
  xyad, hd0, x0, y0, ra0, de0
  adxy, hd1, ra0, de0, x1, y1

  struc = alignstr()
  struc.obj_origin = [x1-x0, y1-y0]

  sz = size(im1)
  new = alignmap(im1, struc, [0, sz[1]-1], [0, sz[2]-1])

  out_hd = hd1
  extast, hd0, astr
  putast, out_hd, astr
;  fits_write, path+'new_'+movable_image, new, hd1
;  hlp, im0, im1, new

  return, new
END

PRO string2sdss, ts_file, run, camcol, rerun, field
  run = fix(strmid(ts_file, 8, 6))
  camcol = fix(strmid(ts_file, 15, 1))
  rerun = fix(strmid(ts_file, 17, 2))
  field = fix(strmid(ts_file, 20, 4))
END

FUNCTION sdss2string, run, camcol, rerun, field, fpc = fpc, psf = psf, $
  extra = extra, tso = tso
;fpc: flag that actually contains the required filter, e.g.: 'r'
  out = ''
  IF NOT keyword_set(extra) THEN extra = ''
  FOR i = 0, n_elements(run)-1 DO BEGIN
    s1 = strtrim(run[i], 2)
    WHILE strlen(s1) LT 6 DO s1 = '0'+s1
    s2 = strtrim(rerun[i], 2)
    WHILE strlen(s2) LT 2 DO s2 = '0'+s2
    s3 = strtrim(field[i], 2)
    WHILE strlen(s3) LT 4 DO s3 = '0'+s3
    IF keyword_set(psf) THEN BEGIN
      out = [out, 'psField-'+s1+'-'+strtrim(camcol[i], 2)+'-'+s3+'.fit']
      CONTINUE
    ENDIF
    IF keyword_set(tso) THEN BEGIN
      out = [out, 'tsObj-'+s1+'-'+strtrim(camcol[i], 2)+'-'+s2+'-'+ $
             s3+extra+'.fit']
      CONTINUE
    ENDIF
    IF keyword_set(fpc) THEN $
      out = [out, 'fpC-'+s1+'-'+fpc+strtrim(camcol[i], 2)+ $
             '-'+s3+extra+'.fit'] $
    ELSE out = [out, 'tsField-'+s1+'-'+strtrim(camcol[i], 2)+'-'+s2+'-'+ $
                s3+extra+'.fit']
  ENDFOR

  return, out[1:*]
END

PRO red_clip_sdss_images, inpath, ts_file, outpath, coords = coords, fpc = fpc
;  inpath = '/home/barden/redshifting_galaxies/best_hi/'
;  outpath = '/home/barden/redshifting_galaxies/best_hi/test/'
;  ts_file = 'tsField-003972-3-40-0035.fit'
;if fpc is set then ts_file may be any fpc file name
;
;this routine interactively clips all 5 SDSS images corresponding to
;one object and aligns them updating the astrometry headers properly.
;
;if coords is set (should contain a size array: left-edge, right-edge,
;bottom-edge, top-edge) non-interactive mode is chosen

;red_clip_sdss_images, '/home/barden/redshifting_galaxies/red_gems/SDSS/', 'fpC-001462-g6-0305.fit', '/home/barden/redshifting_galaxies/red_gems/SDSS/cut/', /fpc


  IF keyword_set(fpc) THEN BEGIN
    left = strmid(ts_file, 0, 11)
    right = strmid(ts_file, 12, 10)
    u_file = left+'u'+right
    g_file = left+'g'+right
    r_file = left+'r'+right
    i_file = left+'i'+right
    z_file = left+'z'+right
  ENDIF ELSE BEGIN
    string2sdss, ts_file, run, camcol, rerun, field
    u_file = sdss2string(run, camcol, rerun, field, fpc = 'u')
    g_file = sdss2string(run, camcol, rerun, field, fpc = 'g')
    r_file = sdss2string(run, camcol, rerun, field, fpc = 'r')
    i_file = sdss2string(run, camcol, rerun, field, fpc = 'i')
    z_file = sdss2string(run, camcol, rerun, field, fpc = 'z')
  ENDELSE

  fits_read, inpath+u_file, u, uhd, /noscale
  fits_read, inpath+g_file, g, ghd, /noscale
  fits_read, inpath+r_file, r, rhd, /noscale
  fits_read, inpath+i_file, i, ihd, /noscale
  fits_read, inpath+z_file, z, zhd, /noscale

  us = shift_images(r, rhd, u, uhd, out_hd = us_hd)
  gs = shift_images(r, rhd, g, ghd, out_hd = gs_hd)
  is = shift_images(r, rhd, i, ihd, out_hd = is_hd)
  zs = shift_images(r, rhd, z, zhd, out_hd = zs_hd)

  IF keyword_set(coords) THEN BEGIN
    minx = coords[0]
    maxx = coords[1]
    miny = coords[2]
    maxy = coords[3]
  ENDIF ELSE BEGIN
    stack = us+gs+r+is+zs
    temp_file = inpath+strcompress(systime(), /remove_all)+'.fits'
    fits_write, temp_file, stack, rhd
    clip_image, temp_file, temp_file, minx, maxx, miny, maxy, /interactive
    spawn, 'rm '+temp_file
  ENDELSE

  fits_write, outpath+u_file, us, us_hd
  fits_write, outpath+g_file, gs, gs_hd
  fits_write, outpath+r_file, r, rhd
  fits_write, outpath+i_file, is, is_hd
  fits_write, outpath+z_file, zs, zs_hd

  clip_image, outpath+u_file, outpath+u_file, minx, maxx, miny, maxy
  clip_image, outpath+g_file, outpath+g_file, minx, maxx, miny, maxy
  clip_image, outpath+r_file, outpath+r_file, minx, maxx, miny, maxy
  clip_image, outpath+i_file, outpath+i_file, minx, maxx, miny, maxy
  clip_image, outpath+z_file, outpath+z_file, minx, maxx, miny, maxy
END

