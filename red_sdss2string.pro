;convert run, camcol, rerun, field numbers into other file names
FUNCTION red_sdss2string, run, camcol, rerun, field, fpc = fpc, psf = psf, $
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

  out = out[1:*]
  IF n_elements(out) EQ 1 THEN out = out[0]

  return, out
END
