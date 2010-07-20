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

FUNCTION red_shift_match_psf, path, ts_file, post, psfpos0, ctrpos, tmpdir, $
                              box0, sdss_psf=sdss_psf, fwhmonly=fwhmonly
;red_shift_match_psf, '/home/barden/redshifting_galaxies/sample/data/', 'tsField-001896-2-41-0065.fit', '.sm', [1551.2723, 1031.1382], '/home/barden/redshifting_galaxies/sample/data/tmp_batch3/', 20

;psfpos: array with positions of stars in 5 bands: 
;[[xu,yu],[xg,yg],[xr,yr],[xi,yi],[xz,yz]]
;ctrpos: position used for centring

   IF box0 MOD 2 THEN box = box0 ELSE box = box0+1

   IF n_elements(psfpos0) LT 10 THEN $
    psfpos = [[psfpos0[0], psfpos0[1]], $
              [psfpos0[0], psfpos0[1]], $
              [psfpos0[0], psfpos0[1]], $
              [psfpos0[0], psfpos0[1]], $
              [psfpos0[0], psfpos0[1]]] $
   ELSE psfpos = psfpos0

   file_mkdir, tmpdir

   red_string2sdss, ts_file, run, camcol, rerun, field
   psfield = red_sdss2string(run, camcol, rerun, field, /psf)

   filters = ['u', 'g', 'r', 'i', 'z']

;read a PSF at position PSFPOS in each band
   IF keyword_set(sdss_psf) THEN BEGIN
      psf = red_read_sdss_psf(path+psfield, tmpdir, psfpos[*, 0], 0)
      FOR i=1, 4 DO BEGIN
         psf = [[[psf]], $
                [[red_read_sdss_psf(path+psfield, tmpdir, psfpos[*, i], i)]]]
      ENDFOR
   ENDIF ELSE BEGIN
      fpc = red_sdss2string(run, camcol, rerun, field, fpc = filters[0])
      fits_read, path+fpc, psf
      pos = fix(round(psfpos))
      pbox = 15
      psf = psf[pos[0, 0]-pbox:pos[0, 0]+pbox, pos[1, 0]-pbox:pos[1, 0]+pbox]
      FOR i = 1, 4 DO BEGIN
         fpc = red_sdss2string(run, camcol, rerun, field, fpc = filters[i])
         fits_read, path+fpc, im
         psf = [[[psf]], [[im[pos[0, i]-pbox:pos[0, i]+pbox, $
                              pos[1, i]-pbox:pos[1, i]+pbox]]]]
      ENDFOR
      psf = float(psf)
      FOR i = 0, 4 DO BEGIN
;centre the psf
         dum = gauss2dfit(psf[*, *, i], par)
         psf[*, *, i] = sshift2d(psf[*, *, i], 15-par[4:5])
;subtract sky and normalise
         sky = ring_sky(psf[*, *, i], box*2, box*0.3, /nw)
         psf[*, *, i] -= sky
         psf[*, *, i] /= total(psf[*, *, i])
      ENDFOR
   ENDELSE

;find the band with the worst seeing and define PSFlo
   worst = 0
   FOR i = 0, 4 DO BEGIN
      dum = gauss2dfit(psf[*, *, i], par)
      sig = mean(par[2:3])
;print, i, sig
      IF sig GT worst THEN BEGIN
         loidx = i
         worst = sig
         psflo = psf[*, *, i]
      ENDIF
   ENDFOR

;read the image corresponding to worst PSF
   fpc = red_sdss2string(run, camcol, rerun, field, fpc = filters[loidx])
   fits_read, path+fpc, im, hd

;make PSF round and a little smoother
   nx = 31 & ny = 31
   amo_x = findgen(nx)#replicate(1, ny)
   amo_y = replicate(1, nx)#findgen(ny)
   dum = gauss2dfit(psflo, par, /tilt)
   smo = max(par[2:3])*1.1
   IF keyword_set(fwhmonly) THEN GOTO, finish
   psfnom = exp(-0.5*(((amo_x-15)/smo)^2+((amo_y-15)/smo)^2))
   trans = red_deconvolve(psfnom, psflo)

   openw, 1, path+'psfFWHMpix'+strmid(ts_file, 8, 16)
   printf, 1, smo*2.3548
   close, 1

;smooth image with new PSF
   psflo = convolve(psflo, trans)
   psflo /= total(psflo)
   im = convolve(im, trans)

   fits_write, path+strrep(strrep(fpc, '.fit', post+'.fit'), 'fpC', 'psf'), $
               psflo
   fits_write, path+strrep(fpc, '.fit', post+'.fit'), im, hd

;fit the reference position with a 2d Gaussian
   sz = [sxpar(hd, 'NAXIS1'), sxpar(hd, 'NAXIS2')]
   x = [ctrpos[0]-1-box/2 > 0 < (sz[0]-1), ctrpos[0]-1+box/2 > 0 < (sz[0]-1)]
   y = [ctrpos[1]-1-box/2 > 0 < (sz[1]-1), ctrpos[1]-1+box/2 > 0 < (sz[1]-1)]
   dum = gauss2dfit(im[x[0]:x[1], y[0]:y[1]], plo)
   plo = plo[2:5]

;  print, plo

;loop over all other bands
   FOR i = 0, 4 DO BEGIN
      IF i EQ loidx THEN CONTINUE
;read PSF and define image
      psfhi = psf[*, *, i]
      psfhi /= total(psfhi)
      fpc = red_sdss2string(run, camcol, rerun, field, fpc = filters[i])
      fits_read, path+fpc, im, hd
      trans = red_deconvolve(psfnom, psfhi)
      psfhi = convolve(psfhi, trans)
      psfhi /= total(psfhi)

;smooth image
      im = convolve(im, trans)
;find shift
      dum = gauss2dfit(im[x[0]:x[1], y[0]:y[1]], phi)
      phi = phi[2:5]

;shift image
      im = sshift2d(im, plo[2:3]-phi[2:3])
      dum = gauss2dfit(im[x[0]:x[1], y[0]:y[1]], phi)
      phi = phi[2:5]

;    print, phi

      fits_write, path+strrep(strrep(fpc, '.fit', post+'.fit'), 'fpC', $
                              'psf'), psfhi
      fits_write, path+strrep(fpc, '.fit', post+'.fit'), im, hd
   ENDFOR

finish:
   return, 2.3548*smo
END
