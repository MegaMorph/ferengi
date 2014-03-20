;routine to calculate a sky value from a postage stamp
;routine should be replaced with the proper sky finder
FUNCTION red_find_sky, im, width, twosided = twosided, sigma = sigma
  !except = 0
  IF NOT keyword_set(sigma) THEN sigma = [-1.5, 0.5]
  sz = size(im)
  IF 2*width GE sz[1] THEN width = (sz[1]-1)/2
  IF 2*width GE sz[2] THEN width = (sz[2]-1)/2
  
  e1 = im[0:*, 0:width-1]
  e2 = im[0:*, sz[2]-width:sz[2]-1]
  e3 = im[0:width-1, width:sz[2]-width-1]
  e4 = im[sz[1]-width:sz[1]-1, width:sz[2]-width-1]
  e1 = reform(e1, n_elements(e1), /overwrite)
  e2 = reform(e2, n_elements(e2), /overwrite)
  e3 = reform(e3, n_elements(e3), /overwrite)
  e4 = reform(e4, n_elements(e4), /overwrite)
  sky = [e1, e2, e3, e4]

;  plothist, sky, xr = [1110, 1160], /peak

;to calculate the sky fit a one-sided Gaussian to the good pixels
  resistant_mean, sky, 3, mn, s, n
  sig = s*sqrt(n_elements(sky)-1-n)
  xrange = sigma*sig+min(mn)
  bin = (xrange[1]-xrange[0])/50.
  plothist, sky, x, y, bin = bin, xrange = xrange, xstyle = 1, /noplot
  idx = where(x GE xrange[0] AND x LE xrange[1], ct)
  IF ct GT 0 THEN BEGIN
    x = x[idx]
    y = y[idx]
  ENDIF
  yfit = gaussfit(x, y, gf, nterms = 3)
  sky = gf[1]

; the following section is for a two-sided Gaussian
  IF keyword_set(twosided) THEN BEGIN
    par = [gf, gf[2]]
    yfit = curvefit(x, y, noweight, par, sigma, $
                    FUNCTION_name = 'curvefit_gauss2side', /noderivative, $
                    status = status, iter = iter)
    
    sky = par[1]
;  skysigma = par[2:3]
  ENDIF

;  xx = findgen(1000)/1000.*(xrange[1]-xrange[0])+xrange[0]
;  oplot, xx, gauss(xx, gf[2], x0 = gf[1], norm = 1, /sigma), col = 50
;  ver, sky
;  par[0] = 1
;  IF keyword_set(twosided) THEN oplot, xx, gauss2side(xx, par), col = 250

  !except = 1
  return, sky
END
