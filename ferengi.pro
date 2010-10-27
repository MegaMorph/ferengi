;+
; NAME: 
;   ferengi
;
; PURPOSE: 
;   apply the effects of redshift to a local galaxy image
;
; EXPLANATION: 
;   "Full and Efficient Redshifting of Ensembles of Nearby Galaxy Images"
;   see the Webpage (http://www.mpia-hd.mpg.de/FERENGI/), paper
;   (http://www.mpia-hd.mpg.de/FERENGI/paper/ferengi_V20070725.pdf),
;   and README for details.
;
; CALLING SEQUENCE: 
;   ferengi, sky, im, imerr, psflo, err0_mag, psfhi, $
;            lambda_lo, filter_lo, zlo, scllo, zplo, tlo, $
;            lambda_hi, filter_hi, zhi, sclhi, zphi, thi, $
;            im_out_file, psf_out_file, $
;            noflux=noflux, evo=evo, noconv=noconv, countsout=countsout
;
; INPUTS:
;   sky = high redshift sky background image
;         2d float array in units [counts/sec]
;         layout is: [x,y] (pixels x, pixels y)
;   im = local input image cube
;        3d float array in units [cts]
;        layout is: [x,y,bands] (pixels x, pixels y, filters)
;   imerr = local input error image cube (Poisson noise image)
;           3d float array in units [cts]
;           layout is: [x,y,bands] (pixels x, pixels y, filters)
;   psflo = low redshift PSF
;           3d float array
;           layout is: [x,y,bands] (pixels x, pixels y, filters)
;           [total flux normalised to 1!]
;   err0_mag = minimum errors to apply to K-correction
;              1d float array in magnitudes
;              layout is: [bands] (filters)
;              for SDSS: err0_mag=[0.05, 0.02, 0.02, 0.02, 0.03]
;   psfhi = high redshift PSF
;           2d float array
;           layout is: [x,y] (pixels x, pixels y)
;           [total flux normalised to 1!]
;   lambda_lo = central/effective wavelength of local filters
;               1d float array in units [angstroems]
;               layout is: [bands] (filters)
;   filter_lo = local filters
;               1d string array
;               layout is: [bands] (filters)
;               for details see description of parameter FILTERLIST to KCORRECT
;   zlo = local redshift
;         float
;   scllo = local pixel scale
;           float in units [arcseconds per pixel]
;   zplo = local zeropoints
;          1d float array in units [magnitudes]
;          layout is: [bands] (filters)
;   tlo = local exposure time (array: [nbands])
;         1d float array in units [seconds]
;         layout is: [bands] (filters)
;   lambda_hi = high redshift wavelength
;               float in units [angstroems]
;   filter_hi = high redshift filters
;               string
;               for details see description of parameter FILTERLIST to KCORRECT
;   zhi = high redshift
;         float
;   sclhi = high redshift pixel scale
;           float in units [arcseconds per pixel]
;   zphi = high redshift zeropoint
;          float in units [magnitudes]
;   thi = high redshift exposure time
;         float in units [seconds]
;
; OPTIONAL INPUTS: 
;   None
;
; OPTIONAL KEYWORD PARAMETERS: 
;   /noflux - do not dim the total flux of the input
;             image. Effectively disable surface brightness
;             dimming. Note, sizes are still correctly scaled down
;   evo - appply an evolution correction. Fluxes are scaled according
;         to: M' = evo × z + M (M': corrected magnitude, evo:
;         evolutionary correction in magnitudes, z: redshift, M:
;         uncrorrected magnitude), i.e. evo=-1 essentially brightens
;         the galaxy by 1 magnitude at redshift z=1
;   /noconv - the output image is not convolved with the output PSF,
;             and no noise is added. This option is used for testing
;             only.
;   /countsout - output image in counts (i.e. same as the input image)
;                rather than counts/sec
;
; OUTPUTS: 
;   im_out_file = filename (including path) for the output image
;                 string
;   psf_out_file = filename (including path) for the output PSF image
;                  string
;
; SIDE EFFECTS: 
;   If the output images exist in advance, they are overwritten.
;
; NOTES: 
;   If im,imerr,psflo,err0_mag,lambda_lo,zplo contain only 1 band, no
;   K-correction is applied.
;   All local input images must be registered to the same pixel frame
;   (source positions must match; pixel scales must be the same).
;
; EXAMPLE:           
;   Redshift a galaxy taken with SDSS out to redshift 0.1 and
;   reobserve it with the SDSS telescope, but with ten times the
;   integration time.
;
;   Read in the input images and compose them into the array form
;   required by FERENGI
;   IDL> path = '/home/test/'
;   IDL> fits_read, path+'fpC-003530-u6-0314.sm.fit', u
;   IDL> fits_read, path+'fpC-003530-g6-0314.sm.fit', g
;   IDL> fits_read, path+'fpC-003530-r6-0314.sm.fit', r
;   IDL> fits_read, path+'fpC-003530-i6-0314.sm.fit', i
;   IDL> fits_read, path+'fpC-003530-z6-0314.sm.fit', z
;   
;   IDL> im = [[[u]], [[g]], [[r]], [[i]], [[z]]]
;
;   For the sake of simplicity we simply take the sqrt() of the input
;   image as the error image
;   IDL> imerr = sqrt(abs(im))
;
;   Read in the corresponding PSF images and compose the 3d array
;   IDL> fits_read, path+'psf-003530-u6-0314.sm.fit', u
;   IDL> fits_read, path+'psf-003530-g6-0314.sm.fit', g
;   IDL> fits_read, path+'psf-003530-r6-0314.sm.fit', r
;   IDL> fits_read, path+'psf-003530-i6-0314.sm.fit', i
;   IDL> fits_read, path+'psf-003530-z6-0314.sm.fit', z
;   
;   IDL> psflo = [[[u]], [[g]], [[r]], [[i]], [[z]]]
;
;   Run the redshifter. As background sky image take an empty array
;   (1000x1000 pixel). The zeropoints are arbitrarily set to 25.
;   IDL> ferengi, fltarr(1000, 1000), im, imerr, psflo, $
;                 [0.05, 0.02, 0.02, 0.02, 0.03], r, $
;                 [3561., 4718, 6185, 7501, 8962], $
;                 ['sdss_u0.par', 'sdss_g0.par', 'sdss_r0.par', $
;                  'sdss_i0.par', 'sdss_z0.par'], 2433/300000., 0.396,
;                 fltarr(5)+25., fltarr(5)+60., 6185., 'sdss_r0.par',
;                 0.025, 0.396, 25., 600., $
;                 path+'imout.fits', path+'psfout.fits'
;
; PROCEDURES CALLED:
;   RSI: FILEPATH, INT_TABULATED, PATH_SEP, STRSPLIT, UNIQ,
;        DOC_LIBRARY, CONGRID, GAUSS2DFIT, CURVEFIT, GAUSSFIT,
;        POLY_FIT, RSTRPOS, REVERSE
;   ASTROLIB: DELVARX, FITS_WRITE, RESISTANT_MEAN, ASINH, AVG,
;             CONVOLVE, FREBIN, LUMDIST, MRDFITS, ROBUST_LINEFIT,
;             SXPAR, FITS_READ, FITS_OPEN, FITS_CLOSE, DIST_ELLIPSE,
;             FXADDPAR, HOST_TO_IEEE, IEEE_TO_HOST, COSMO_PARAM,
;             MATCH, MRD_HREAD, MRD_SKIP, QSIMP, FXMOVE, FXPAR,
;             FXPARPOS, FXPOSIT, SXADDPAR, SXDELPAR, DETABIFY, GETTOK,
;             IS_IEEE_BIG, MRD_STRUCT, ROBUST_SIGMA, ROB_CHECKFIT,
;             VALID_NUM, HPRINT, TRAPZD, ZPARCHECK, EXPAND_TILDE,
;             NINT, NUMLINES, REPSTR, FDECOMP, SPEC_DIR
;   KCORRECT: KCORRECT, K_RECONSTRUCT_MAGGIES, KCORRECT_SO_EXT,
;             K_FIT_NONNEG, K_SOLAR_MAGNITUDES, LF_DISTMOD, LF_T2Z,
;             LF_Z2T, KLOG, K_ABFIX, K_LOAD_FILTERS, K_LOAD_VMATRIX,
;             K_MINERROR, K_PROJECTION_TABLE, K_READ_BASEL, K_SDSSFIX,
;             K_LAMBDA_TO_EDGES, K_LUPS2MAGGIES, K_PROJECT_FILTERS,
;             K_SDSS_ERR2IVAR, K_STR_SEP, K_READ_ASCII_TABLE
;   IDLUTILS V5.0.0+: SETLOG, SPLOG, SSHIFT2D, STR_SEP, YANNY_READONE,
;                     YANNY_FREE, YANNY_READ, FILEANDPATH,
;                     HOGG_UNQUOTED_REGEX, HOGG_STRSPLIT
;
; MODIFICATION HISTORY: 
;       Written by M. Barden, Innsbruck, 2007 <marco.barden@uibk.ac.at>
;
;-

FUNCTION nu2lam, nu             ;[in Hz]
   return, 299792458./nu*1e-10  ;[in Angstroem]
END

FUNCTION lam2nu, lam            ;[in Angstroem]
   return, 299792458./lam*1e10  ;[in Hz]
END

FUNCTION maggies2mags, maggies
   return, -2.5*alog10(maggies)
END

FUNCTION mags2maggies, mags
   return, 10^(-0.4*mags)
END

FUNCTION maggies2fnu, maggies
   return, 3631e-23*maggies     ;[erg s-1 Hz-1 cm-2]
END

FUNCTION fnu2maggies, fnu
   return, 3631e23*fnu
END

FUNCTION fnu2flam, fnu, lam     ;[erg s-1 Hz-1 cm-2], [angstroem]
   return, 299792458.*1e10/lam^2.*fnu ;[erg s-1 cm-2 A-1]
END

FUNCTION flam2fnu, flam, lam    ;[erg s-1 cm-2 A-1], [angstroem]
   return, flam/299792458./1e10*lam^2. ;[erg s-1 Hz-1 cm-2]
END

FUNCTION lambda_eff, lam, trans
;calculate effective wavelength for a filter
   idx = where(lam NE 0, ct)
   IF ct EQ 0 THEN message, 'ERROR: no non-zero wavelengths'
   return, int_tabulated(reform(lam[idx]), reform(lam[idx]*trans[idx])) $
           /int_tabulated(reform(lam[idx]), reform(trans[idx]))
END

FUNCTION cts2mags, cts, expt, zp
;zp: zeropoint is a positive number
   return, maggies2mags(cts2maggies(cts, expt, zp))
END

FUNCTION cts2maggies, cts, expt, zp
;zp: zeropoint is a positive number
   return, cts / expt * 10^(-0.4*zp)
END

FUNCTION mags2cts, mags, expt, zp
;zp: zeropoint is a positive number
   return, maggies2cts(mags2maggies(mags), expt, zp)
END

FUNCTION maggies2cts, maggies, expt, zp
;zp: zeropoint is a positive number
   return, maggies * expt / 10^(-0.4*zp)
END

FUNCTION maggies2lup, maggies, filter
;filter is a single element string (not an array)
   CASE filter OF 
      'u': b = 1.4e-10
      'g': b = 0.9e-10
      'r': b = 1.2e-10
      'i': b = 1.8e-10
      'z': b = 7.4e-10
   ENDCASE

   return, -2.5/alog(10)*(asinh(maggies/b*0.5)+alog(b))
END

FUNCTION lup2maggies, lup, filter
;filter is a single element string (not an array)
   maggies = lup

   CASE filter OF 
      'u': b = 1.4e-10
      'g': b = 0.9e-10
      'r': b = 1.2e-10
      'i': b = 1.8e-10
      'z': b = 7.4e-10
   ENDCASE

   return, 2*b*sinh(-0.4*alog(10)*lup-alog(b))
END

FUNCTION random_indices, len, n_in
;produces for an array with LEN elements a set of N_IN indices. The
;returend indices do not contain duplicates.
;Taken from: http://www.dfanning.com/code_tips/randomindex.html

   swap = n_in gt len/2
   IF swap THEN n = len-n_in ELSE n = n_in
   inds = LonArr(n, /NOZERO)
   M = n
   WHILE n GT 0 DO BEGIN
      inds[M-n] = Long( RandomU(seed, n)*len )
      inds = inds[Sort(inds)]
      u = Uniq(inds)
      n = M-n_elements(u)
      inds[0] = inds[u]
   ENDWHILE
   
   IF swap THEN inds = Where(Histogram(inds,MIN=0,MAX=len-1) EQ 0)
   RETURN, inds
END

FUNCTION edge_index, a, rx, ry
;the routine creates an index of a ring with width 1 around the
;centre at radius rx and ry. E.g. rx=1 & ry=0:
;
;       0       1       2       3       4
;       5       6       7       8       9
;      10      11      12      13      14
;      15      16      17      18      19
;      20      21      22      23      24
;      25      26      27      28      29
;
;       0       1       2       3       4
;       5       6       7       8       9
;      10      -1      -1      -1      14
;      15      -1      -1      -1      19
;      20      21      22      23      24
;      25      26      27      28      29
  on_error, 2

  sz = size(a)
  IF sz[1] MOD 2 THEN px = 0 ELSE px = 0.5
  IF sz[2] MOD 2 THEN py = 0 ELSE py = 0.5
  b = fix(abs((lindgen(sz[1], sz[2]) MOD sz[1])-sz[1]/2+px))
  c = fix(abs(transpose((lindgen(sz[2], sz[1]) MOD sz[2])-sz[2]/2+py)))
  i = where(b EQ rx AND c LE ry OR c EQ ry AND b LE rx)

  return, i
END

FUNCTION ring_sky, image, width0, nap, x, y, q, pa, rstart=rstart, nw=nw
;for an image (if 1 element assumed filename, else array) measure
;flux in aperture around position x,y in rings with axis ratio q and
;position angle pa (measured from up counterclockwise)
;nap: number of apertures used for sky slope measurement
;rstart: starting radius
;if ellipse parameters are not known: use circular aperture and image centre
;nw: if set, width is not the width of an annulus in pixels, but the
;    number of annuli calculated

;number of apertures must be greater than 3, otherwise resistant_mean
;complains
   nap >= 4

;read image, prepare array
   IF n_elements(image) EQ 1 THEN fits_read, image, im ELSE im = image

   sz = (size(im))[1:2]

;if not specified assume circular parameters for ellipse
   IF n_params() LT 7 THEN BEGIN
      x = sz[0]*0.5
      y = sz[1]*0.5
      q = 1.
      pa = 0.
      rstart = min(sz)*0.05
   ENDIF

;create radius array
   dist_ellipse, rad, sz, x, y, 1./q, pa

   max_rad = max(rad)*0.95

   IF keyword_set(nw) THEN width = max_rad/float(width0) $
   ELSE width = width0

;get global mean and sigma
   resistant_mean, im, 3, mean, sig, nrej
   sig *= sqrt(n_elements(im)-1-nrej)

;maximum radius of current aperture
   IF keyword_set(rstart) THEN rhi = rstart ELSE rhi = width
;setup array for flux,radius curve
   flux = (r = 0.)
;counter
   i = 0
;slope of last nap apertures is negativ
;(also defines minimum number of positive slope measurements:
;-1=2measurements, 0=1measurement)
   sign = -1
   WHILE rhi LE max_rad DO BEGIN
      extra = 0
      ct = 0
;select aperture and reject 3sigma outliers
      WHILE ct LT 10 DO BEGIN
         idx = where(rad LE rhi+extra AND rad GE rhi-width-extra AND $
                     im LE mean+3*sig AND im GE mean-3*sig, ct)
         extra++
         IF extra GT max(sz)*2 THEN BREAK
      ENDWHILE
;if enough pixels left get robust mean flux in current aperture
      IF ct LT 5 THEN sky = flux[n_elements(flux)-1] $
      ELSE resistant_mean, im[idx], 3, sky
;enlarge array with radius and flux
      r = [r, rhi-0.5*width]
      flux = [flux, sky]
;increase counter
      i++
      IF n_elements(flux) GT nap THEN BEGIN
;if enough measurements available, measure slope
         co = robust_linefit(r[i-nap+1:i], flux[i-nap+1:i])
;stop if slope is positive for the first time stop
         IF sign GT 0 AND co[1] GT 0 THEN BREAK
;if slope is positive remember
         IF co[1] GT 0 THEN sign++
      ENDIF
;increase radius
      rhi += width
   ENDWHILE

;calculate robust sky over last nap apertures
   resistant_mean, flux[(i-nap+1 > 0):i], 3, sky

   return, sky
END

PRO ferengi_make_psf_same, psf1, psf2
;enlarge smaller psf to the size of the larger one (zero-padding)
  flag_c = 1
  IF n_elements(psf1) GT n_elements(psf2) THEN BEGIN
    big = psf1
    small = psf2
  ENDIF ELSE BEGIN
    big = psf2
    small = psf1
    flag_c = 0
  ENDELSE
  szbig = size(big)
  szsmall = size(small)
  mid = szbig[1]/2
  len = szsmall[1]
  lo = mid-len/2
  hi = mid+len/2
  small0 = big*0
  small0[lo:hi, lo:hi] = small
  IF flag_c THEN psf2 = small0 ELSE psf1 = small0
END

FUNCTION ferengi_psf_centre, psf0, xshift, yshift
  psf = psf0
  !except = 0

;resize to odd number of pixels
  sz = size(psf)
  IF (sz[1]+1) MOD 2 THEN psf = [psf, fltarr(1, sz[2])]
  sz = size(psf)
  IF (sz[2]+1) MOD 2 THEN psf = [[psf], [fltarr(sz[1])]]

;centre PSF
  sz = (size(psf))[1:2]/2
  dum = gauss2dfit(psf, par)
  psf = sshift2d(psf, sz-par[4:5])

  xshift = (sz-par[4:5])[0] & yshift = (sz-par[4:5])[1]

  !except = 1
  return, psf
END

FUNCTION ferengi_deconvolve, wide, narrow
;make sure WIDE and NARROW have same size, odd number of pixels, are
;centred and normalised

  sz = max([(size(wide))[1:2], (size(narrow))[1:2]])
  bigsz = 2
  WHILE bigsz LT sz DO bigsz *= 2
  IF bigsz GT 2048 THEN message, 'Requested PSF array is larger than 2x2k!'
;bigsz = 2048

  psf_n_2k = dblarr(bigsz, bigsz)
  psf_w_2k = dblarr(bigsz, bigsz)
  sz = size(narrow)
  psf_n_2k[0:sz[1]-1, 0:sz[2]-1] = narrow
  psf_w_2k[0:sz[1]-1, 0:sz[2]-1] = wide

;get FFT FOR both
  psf_n_2k = dcomplex(psf_n_2k, 0)
  psf_w_2k = dcomplex(psf_w_2k, 0)
  fft_n = FFT(psf_n_2k, /double)
  fft_w = FFT(psf_w_2k, /double)

  delvarx, psf_n_2k, psf_w_2k
  fft_n = (ABS(fft_n)/(ABS(fft_n)+0.000000001))*fft_n
  fft_w = (ABS(fft_w)/(ABS(fft_w)+0.000000001))*fft_w

;calculate ratio
  psfrat = fft_w/fft_n
  delvarx, fft_w, fft_n

;create transformation psf
  psfhlp = double(FFT(psfrat, /double))
  psfcorr = fltarr(sz[1], sz[2])
  lo = bigsz-sz[1]/2
  hi = sz[1]/2
  psfcorr[0:hi-1,0:hi-1] = psfhlp[lo:bigsz-1, lo:bigsz-1]
  psfcorr[hi:sz[1]-1,0:hi-1] = psfhlp[0:hi, lo:bigsz-1]
  psfcorr[hi:sz[1]-1, hi:sz[1]-1] = psfhlp[0:hi, 0:hi]
  psfcorr[0:hi-1, hi:sz[1]-1] = psfhlp[lo:bigsz-1, 0:hi]
  delvarx, psfhlp
  psfcorr = rotate(psfcorr, 2)

  return, psfcorr/total(psfcorr)
END

PRO ferengi_clip_edge, npix, im, auto_frac = auto_frac, $
                       clip_also = clip_also, norm = norm
;npix: number of outer pixels to be clipped
;im: image array that is to be clipped
;auto_frac: fraction of the radius of the image at which automatic
;clipping is started (default 50%, i.e. auto_frac=2)
;clip_also: clip an additional array as well (provide array variable)
;norm: flag. If set im (and clip_also) are normalised

  IF NOT keyword_set(auto_frac) THEN auto_frac = 2
  sz = size(im)
  rx = fix(sz[1]/2/auto_frac)
  ry = fix(sz[2]/2/auto_frac)
  sig = 0.
  r = 0
  WHILE 1 DO BEGIN
    i = edge_index(im, rx, ry)
    IF i[0] EQ -1 THEN BREAK
    resistant_mean, im[i], 3, mn, sg, nr
    sg *= sqrt(n_elements(i)-1-nr)
    sig = [sig, sg]
    r = [r, rx]
    ++rx
    ++ry
  ENDWHILE
  r = r[1:*]
  sig = sig[1:*]
  resistant_mean, sig, 3, mn, sg, nr
  sg = sg*sqrt(n_elements(sig)-1-nr)
  i = where(sig GT mn+10*sg, ct)
  IF ct GT 0 THEN BEGIN
    lim = min(r[i])-1
    IF ct GT nr*3 THEN message, 'Large gap?'
    npix = round(sz[1]/2.-lim)

    IF keyword_set(clip_also) THEN $
      clip_also = clip_also[npix:sz[1]-1-npix, npix:sz[2]-1-npix]
    im = im[npix:sz[1]-1-npix, npix:sz[2]-1-npix]
  ENDIF

  IF keyword_set(norm) THEN BEGIN
    IF keyword_set(clip_also) THEN clip_also /= total(clip_also)
    im /= total(im)
  ENDIF
END

FUNCTION ferengi_downscale, im_lo, z_lo, z_hi, p_lo, p_hi, $
                            upscl = upscl, nofluxscl = nofluxscl, evo = evo
;calculate the scaling in flux and size for a given set of redshifts
;without applying any K-corrections, and without PSF convolution
;im_lo: the input image
;z_lo,z_hi: the redshift of the local and the distant object
;p_lo,p_hi: sizes of the low and high redshift pixels
;upscl: enlarge instead of scale down
;nofluxscl: do not scale the flux
;evo: evolutionary correction: scale the flux evo*zhi magnitudes up
;(negative values of evo increase the flux at high redshift)
   IF NOT keyword_set(evo) THEN evo_fact = 1 $
   ELSE evo_fact = 10^(-0.4*evo*z_hi)

   d_lo = lumdist(z_lo, /silent)
   d_hi = lumdist(z_hi, /silent)

;the magnification (size correction)
   magnification = (d_lo/d_hi*(1.+z_hi)^2/(1.+z_lo)^2*p_lo/p_hi)[0]
   IF keyword_set(upscl) THEN magnification = 1./magnification
;the flux scaling (surface brightness dimming)
   IF keyword_set(nofluxscl) THEN flux_ratio = 1. $
   ELSE flux_ratio = (d_lo/d_hi)[0]^2

;calculate the number of pixels in the input image
   sz_lo = (size(im_lo))[1:2]
;calculate the number of pixels of the scaled image
   nx_hi = round(sz_lo[0]*magnification)
   ny_hi = round(sz_lo[1]*magnification)
;create a new scaled version of the input image roughly in the centre
;of an array of the size of the input image and apply the surface
;brightness dimming correction

   return, frebin(im_lo, nx_hi, ny_hi, /total)*flux_ratio*evo_fact
END

FUNCTION centroid, array
  s = Size(array, /Dimensions)
  totalMass = Total(array)
  xcm = Total( Total(array, 2) * Indgen(s[0]) ) / totalMass
  ycm = Total( Total(array, 1) * Indgen(s[1]) ) / totalMass
  RETURN, [xcm, ycm]
END

FUNCTION ferengi_odd_n_square, psf0, centre=centre
;make the input PSF array square, the number of pixels along each axis
;odd and centre the result
;centre: central position (do not centre using Gauss-Fit)
   psf = psf0

;resize to odd number of pixels
   sz = size(psf)
   IF (sz[1]+1) MOD 2 THEN psf = [psf, fltarr(1, sz[2])]
   sz = size(psf)
   IF (sz[2]+1) MOD 2 THEN psf = [[psf], [fltarr(sz[1])]]

;make array square
   sz = size(psf)
   IF sz[1] GT sz[2] THEN psf = [[psf], [fltarr(sz[1], sz[1]-sz[2])]]
   IF sz[2] GT sz[1] THEN psf = [psf, fltarr(sz[2]-sz[1], sz[2])]

;centre array
   IF n_elements(centre) EQ 2 THEN BEGIN
      psf = sshift2d(psf, centre)
   ENDIF ELSE BEGIN
      sz = (size(psf))[1:2]

      CASE 1 OF
         sz[0] LE 1: psf = 1.
         sz[0] EQ 3 OR sz[0] EQ 5: BEGIN
            psf = rebin(psf, sz[0]*3, sz[0]*3, /sample)
            dum = gauss2dfit(transpose([transpose([psf, psf, psf]*0), $
                                        transpose([psf*0, psf, psf*0]), $
                                        transpose([psf, psf, psf]*0)]), par)
            psf = sshift2d(psf, 3*sz/2-par[4:5]+sz*3)
            psf = rebin(psf, sz[0], sz[0], /sample)
         END
         ELSE: BEGIN
            dum = gauss2dfit(transpose([transpose([psf, psf, psf]*0), $
                                        transpose([psf*0, psf, psf*0]), $
                                        transpose([psf, psf, psf]*0)]), par)
            sigma = par[2:3]
            shft = sz/2-par[4:5]+sz
            c = centroid(psf)
            cshft = (sz/2)-c
            IF abs(shft[0]) GE sz[0]/2 OR abs(shft[1]) GE sz[1]/2 OR $
               sigma[0] LE 0 OR sigma[1] LE 0 OR $
               sigma[0] GE sz[0]/8 OR sigma[1] GE sz[1]/8 THEN BEGIN
               print, 'Warning: psf centering fit failed, using centroid!'
               shft = cshft
            ENDIF 
            IF abs(shft[0]) LT 0.01 THEN shft[0] = 0
            IF abs(shft[1]) LT 0.01 THEN shft[1] = 0
            IF total(abs(shft)) GT 0 THEN psf = sshift2d(psf, shft)
         END
      ENDCASE
   ENDELSE

   return, psf
END

FUNCTION ferengi_transformation_psf, psf_s0, psf_c0, z_lo, z_hi, p_lo, p_hi, $
                                     same_size = same_size
;calculate transformation psf
;psf_s0: local psf (has to be narrower than psf_c0 after down-scaling)
;psf_c0: highz psf
;z_lo,z_hi: local,highz redshift
;p_lo,p_hi: local,highz pixel scales
  psf_s = psf_s0
  psf_c = psf_c0

;make size odd & make square & centre for both psfs
  psf_s = ferengi_odd_n_square(psf_s)
  psf_c = ferengi_odd_n_square(psf_c)

  d_lo = lumdist(z_lo, /silent)
  d_hi = lumdist(z_hi, /silent)
  insz = (size(psf_s))[1]
  add = 0
  outsz = round((d_lo/d_hi*(1.+z_hi)^2/(1.+z_lo)^2*p_lo/p_hi)[0]*(insz+add))
  WHILE (outsz MOD 2 EQ 0 or outsz le 2) DO BEGIN
    add += 2
    psf_s = [[psf_s], [fltarr((size(psf_s))[1], 2)]]
    outsz = round((d_lo/d_hi*(1.+z_hi)^2/(1.+z_lo)^2*p_lo/p_hi)[0]*(insz+add))
    IF add GT insz*3 THEN message, 'enlarging PSF failed!'
  ENDWHILE
  psf_s = ferengi_odd_n_square(psf_s)

;downscale the local PSF
  psf_s = ferengi_downscale(psf_s, z_lo, z_hi, p_lo, p_hi, /nofluxscl)

;make size the same
  psf_s = ferengi_odd_n_square(psf_s)
  ferengi_make_psf_same, psf_c, psf_s

;make size odd & make square & centre for both psfs
  psf_s = ferengi_odd_n_square(psf_s)
  psf_c = ferengi_odd_n_square(psf_c)

;normalise
  psf_s /= total(psf_s)
  psf_c /= total(psf_c)

  IF keyword_set(same_size) THEN BEGIN
    psf_s0 = psf_s
    psf_c0 = psf_c
  ENDIF

  return, ferengi_deconvolve(psf_c, psf_s)
END

FUNCTION ferengi_convolve_plus_noise, im, psf, sky, exptime, $
                                      nonoise = nonoise, $
                                      border_clip = border_clip, $
                                      extend = extend
;clip the borders of the PSF? If yes, the input psf is changed!
   IF NOT keyword_set(border_clip) THEN border_clip = 0
   sz_psf = (size(psf))[1:2]
   psf = psf[border_clip:sz_psf[0]-1-border_clip, $
             border_clip:sz_psf[1]-1-border_clip]

;in order to convolve every pixel properly the image has to be enlarged
   sz_psf = (size(psf))[1:2]
   sz_im = (size(im))[1:2]
   out = [fltarr(sz_psf[0], sz_im[1]), im, fltarr(sz_psf[0], sz_im[1])]
   sz_out = (size(out))[1:2]
   out = transpose([transpose(fltarr(sz_out[0], sz_psf[1])), transpose(out), $
                    transpose(fltarr(sz_out[0], sz_psf[1]))])

;convolve with the PSF
   out = convolve(out, psf/total(psf))

;remove the excess border
   sz_out = (size(out))[1:2]
   IF NOT keyword_set(extend) THEN $
    out = out[sz_psf[0]:sz_out[0]-1-sz_psf[0], sz_psf[1]:sz_out[1]-1-sz_psf[1]]

;add poisson noise to the image and add the sky
   sz_out = (size(out))[1:2]

;the output image is in counts/sec, and so is the sky image
   IF NOT keyword_set(nonoise) THEN $
    out += sky[0:sz_out[0]-1, 0:sz_out[1]-1]+ $
           sqrt(abs(out*exptime))*randomn(1, sz_out[0], sz_out[1])/exptime

   return, out
END

PRO ferengi, sky, im, imerr, psflo, err0_mag, psfhi, $
             lambda_lo, filter_lo, zlo, scllo, zplo, tlo, $
             lambda_hi, filter_hi, zhi, sclhi, zphi, thi, $
             im_out_file, psf_out_file, $
             noflux=noflux, evo=evo, noconv=noconv, countsout=countsout

;the size of the output sky image
   sz_sky = (size(sky))[1:2]

;number of input filters
   nbands = n_elements(im[0, 0, *])
   IF nbands EQ 1 THEN nok = 1 ELSE nok = 0

;convert from cts (input frame) to maggies and back to cts (output frame)
   IF nok THEN $
    im_nok = maggies2cts(cts2maggies(im, tlo, zplo), thi, zphi)

;select best matching PSF for output redshift
   IF nok THEN psf_lo = psflo ELSE BEGIN
      dz = abs(lambda_hi/lambda_lo-1)
      idx_bestfilt = where(dz EQ min(dz), ct)
      IF ct GT 1 THEN idx_bestfilt = idx_bestfilt[0]
      psf_lo = psflo[*, *, idx_bestfilt]
   ENDELSE

;scale the image down
   IF nok THEN BEGIN
      im_ds = ferengi_downscale(im_nok, zlo, zhi, scllo, sclhi, $
                                nofluxscl=noflux, evo=evo)
   ENDIF ELSE BEGIN
;weight the closest filters in rest-frame more
      dz1 = abs(lambda_hi-lambda_lo)
      ord = sort(dz1)
      weight = fltarr(nbands)+1
      IF dz1[ord[0]] EQ 0 THEN BEGIN
         IF nbands EQ 2 THEN weight[ord] = [10, 4]
         IF nbands EQ 3 THEN weight[ord] = [10, 4, 4]
         IF nbands GE 4 THEN weight[ord] = [10, 4, 4, fltarr(nbands-3)+1]
      ENDIF ELSE BEGIN
         IF nbands EQ 2 THEN weight[ord] = [10, 8]
         IF nbands EQ 3 OR nbands EQ 4 THEN $
          weight[ord] = [10, 8, fltarr(nbands-2)+4]
         IF nbands GT 4 THEN weight[ord] = [10, 8, 4, 4, fltarr(nbands-4)+1]
      ENDELSE

      im_ds = ferengi_downscale(im[*, *, 0], zlo, zhi, $
                                scllo, sclhi, nofluxscl=noflux, evo=evo)
      FOR j=1, nbands-1 DO $
       im_ds = [[[im_ds]], $
                [[ferengi_downscale(im[*, *, j], zlo, zhi, scllo, $
                                    sclhi, nofluxscl=noflux, evo=evo)]]]

;*******************************************************************************
;subtracting sky here: the user may replace this with his/her own routine
      FOR j=0, nbands-1 DO $
       im_ds[*, *, j] -= ring_sky(im_ds[*, *, j], 50, 15, /nw)
;*******************************************************************************

;and the counts-error-image
      imerr_ds = ferengi_downscale(imerr[*, *, 0], zlo, zhi, $
                                   scllo, sclhi, nofluxscl=noflux, evo=evo)
      FOR j=1, nbands-1 DO $
       imerr_ds = [[[imerr_ds]], $
                   [[ferengi_downscale(imerr[*, *, j], zlo, zhi, scllo, $
                                       sclhi, nofluxscl=noflux, evo=evo)]]]
      
;convert the error from cts to mags
      FOR j=0, nbands-1 DO $
       imerr_ds[*, *, j] = 2.5/alog(10)*imerr_ds[*, *, j]/ $
       im_ds[*, *, j]

;calculate the flux in each pixel (convert image from cts to maggies)
      FOR j=0, nbands-1 DO $
       im_ds[*, *, j] = cts2maggies(im_ds[*, *, j], tlo[j], zplo[j])

;+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
;siglim defines the sigma above which K-corrections are calculated
;+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
      siglim = 2
      sig = fltarr(nbands)
      npix = n_elements(im_ds[*, *, 0])

;find the index of the input filter that matches the output wavelength
;best
      zmin = abs(lambda_hi/lambda_lo-1-zhi)
      filt_i = where(zmin EQ min(zmin))

;select the pixels above nsig with resistant_mean
;create a sigma map
      nsig = im_ds*0
      nhi = im_ds[*, *, 0]*0
      FOR j=0, nbands-1 DO BEGIN
         resistant_mean, im_ds[*, *, j], 3, m, s, n
         sig[j] = s*sqrt(npix-1-n)
         nsig[*, *, j] = median(im_ds[*, *, j], 3)/sig[j]
         hi = where(abs(nsig[*, *, j]) GT siglim, ct)
         IF ct GT 0 THEN nhi[hi] += 1
      ENDFOR

;from the "closest" filter select good pixels
      nhi = transpose(reform(nhi, n_elements(nhi), 1))
      good1 = where(abs(nsig[*, *, filt_i]) GT 0.25 AND $
                    abs(nsig[*, *, filt_i]) LE siglim, ct1)
;select only 50% of all pixels with 0.25 < nsig < siglim
      IF ct1 GT 0 THEN good1 = good1[random_indices(ct1, round(ct1*0.5))]
      good = where(nhi GE 3 AND abs(nsig[*, *, filt_i]) GT siglim, ct)
      IF ct GE 0 THEN print, '3+ filters have high sigma pixels'
      IF ct EQ 0 THEN BEGIN
         print, 'less than 3 filters have high sigma pixels'
         good = where(nhi GE 2 AND abs(nsig[*, *, filt_i]) GT siglim, ct)
      ENDIF
      IF ct EQ 0 THEN BEGIN
         print, 'less than 2 filters have high sigma pixels'
         good = where(nhi GE 1 AND abs(nsig[*, *, filt_i]) GT siglim, ct)
      ENDIF
      IF ct EQ 0 THEN BEGIN
         print, 'NO filter has high sigma pixels'
         good = where(nhi GE 0 AND abs(nsig[*, *, filt_i]) GT siglim, ct)
      ENDIF
      IF ct1 GT 0 THEN good = [good, good1]
      ngood = n_elements(good)
      IF ngood EQ 1 AND good[0] EQ -1 THEN GOTO, nopixels

;setup the arrays for the pixels that are to be K-corrected
      maggies = (err = fltarr(nbands, ngood))
      FOR j=0, nbands-1 DO maggies[j, *] = (im_ds[*, *, j])[good]
      FOR j=0, nbands-1 DO err[j, *] = (imerr_ds[*, *, j])[good]
      nsig_2d = fltarr(nbands, ngood)
      FOR j=0, nbands-1 DO nsig_2d[j, *] = (nsig[*, *, j])[good]

;remove infinite values in the error image
      inf = where(finite(err) NE 1, ct)
      IF ct GT 0 THEN err[inf] = 99999
      err = abs(err) < 99999

;reset K-correct
      delvarx, rmatrix, zvals, wavel, vmatrix, coeffs

;setup array with minimum errors for SDSS
      err0 = err0_mag#(fltarr(n_elements(maggies[0, *]))+1)
      wei = weight#(fltarr(n_elements(maggies[0, *]))+1)

;add image errors and minimum errors in quadrature
      err = sqrt(err0^2+err^2)/wei

;convert errors from magnitudes 
      ivar = (2.5/alog(10)/err/double(maggies))^2
      inf = where(finite(ivar) NE 1, ct)
      IF ct GT 0 THEN ivar[inf] = max(ivar[where(finite(ivar) EQ 1)])

      z_tmp = fltarr(n_elements(maggies[0, *]))+zlo

      print, systime()
;calculate K-corrections for redshift z=0
      kcorrect, maggies, ivar, z_tmp, k, $
                filterlist=filter_lo, rmatrix=rmatrix, zvals=zvals, $
                lambda=wavel, vmatrix=vmatrix, coeffs=coeffs
;  plot, wavel, vmatrix#coeffs, xrange = [2000., 12000.], /ylog, ystyle = 1
      
;calculate a new array of output redshifts
      z_tmp = fltarr(n_elements(maggies[0, *]))+zhi

;reconstruct magnitudes in a certain filter at a certain redshift
      k_reconstruct_maggies, coeffs, z_tmp, maggies, vmatrix=vmatrix, $
                             lambda=wavel, filterlist=filter_hi
      print, systime()

nopixels:
;as background choose closest in redshift-space
      bg = (im_ds = im_ds[*, *, filt_i]/(1.+zhi))

;put in K-corrections
      IF ngood GT 0 AND good[0] NE -1 THEN im_ds[good] = maggies/(1.+zhi)

;convert image back to cts
      im_ds = maggies2cts(im_ds, thi, zphi)
      bg = maggies2cts(bg, thi, zphi)
   ENDELSE

;remove infinite pixels: replace with median (3x3)
   med = median(im_ds, 3)
   idx = where(finite(im_ds) NE 1, ct)
   IF ct GT 0 THEN im_ds[idx] = med[idx]

;replace 0-value pixels with median (3x3)
   idx = where(im_ds EQ 0, ct)
   IF ct GT 0 THEN im_ds[idx] = med[idx]

   IF nok EQ 0 THEN BEGIN
      resistant_mean, im_ds, 3, m, sig, nrej
      sig *= sqrt(n_elements(im_ds)-1-nrej)
      idx = where(abs(im_ds) GT 10*sig, ct)
      IF ct GT 0 THEN BEGIN
         fit = robust_linefit(abs(bg[idx]), abs(im_ds[idx]), /bisect)
         if (fit EQ 0) THEN message, 'Fit failed: not enough sky in image?'
         delta = abs(im_ds[idx])-(fit[0]+fit[1]*abs(bg[idx]))
         resistant_mean, delta, 3, m, sig, nrej
         sig *= sqrt(n_elements(im_ds)-1-nrej)
         idx1 = where(delta/sig GT 50, ct)
         IF ct GT 0 THEN im_ds[idx[idx1]] = med[idx[idx1]]
      ENDIF
   ENDIF
   
;*******************************************************************************
;subtracting sky here: the user may replace this with his/her own routine
   im_ds -= ring_sky(im_ds, 50, 15, /nw)
;*******************************************************************************

   IF keyword_set(noconv) THEN BEGIN
      im_ds /= thi
      recon = psf_lo/total(psf_lo)
      GOTO, write_out
   ENDIF

;the output sky image might be too small to put in the galaxy
   sz_im_ds = (size(im_ds))[1:2]
   IF sz_im_ds[0] GT sz_sky[0] OR sz_im_ds[1] GT sz_sky[1] THEN $
    message, 'sky image not big enough'

;calculate the transformation PSF
   psf_hi = psfhi
   psf_t = ferengi_transformation_psf(psf_lo, psf_hi, zlo, zhi, scllo, sclhi, $
                                      /same)

   recon = ferengi_odd_n_square(convolve(psf_lo, psf_t))

;get rid of potential bad pixels around the edges of the PSF
   rem = 3                      ;remove 3 pixels
   ferengi_clip_edge, rem, recon, clip_also = psf_hi, /norm

;normalise reconstructed PSF
   recon /= total(recon)

;convolve the high redshift image with the transformation PSF
   im_ds = ferengi_convolve_plus_noise(im_ds/thi, psf_t, sky, thi, $
                                       border_clip = 3, extend = extend)

; optionally output as counts, rather than counts/sec
   IF keyword_set(countsout) THEN im_ds *= thi

write_out:
   fits_write, im_out_file, im_ds
   fits_write, psf_out_file, recon

alldone:
END
