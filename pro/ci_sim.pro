function ci_par_struc

   par = {width: 3072, height: 2048, $
          nominal_zeropoint: 26.56, $
          ci_extnames: ['CIE', 'CIN', 'CIC', 'CIS', 'CIW']}

   return, par

end

pro _cache_nominal_astrometry

  COMMON _CI_ASTROM, astrom
  if n_elements(astrom) EQ 0 then begin
      fname = '../etc/viewer_astrom_index-as_designed.bigtan.fits'
      astrom = mrdfits(fname, 1)
  end

end

function get_nominal_astrometry, extname

  check_valid_extname, extname

  _cache_nominal_astrometry
  COMMON _CI_ASTROM, astrom

  w = where(astrom.extname EQ extname, nw)

  if (nw NE 1) then stop

  astr = astrom[w[0]]

  return, astr
end

function get_astrometry_radec, telra, teldec, extname

  check_valid_extname, extname

  if n_elements(telra) NE 1 then stop
  if n_elements(teldec) NE 1 then stop

  if (teldec LT -90) OR (teldec GT 90) then stop
  if (telra LT 0) OR (telra GE 360) then stop

  if size(telra, /type) NE 5 then stop
  if size(teldec, /type) NE 5 then stop

; this is centered at (ra, dec) = (0, 0)
  astr = get_nominal_astrometry(extname)

  astr.crval = [telra, teldec]

  return, astr
end

function dark_current_rate, t_celsius
    
  ; t_celsius - temperature in deg celsius
  ; I = I(0 Celsius)*2^(T_celsius/dT), dT = doubling rate in deg C
  ; I(0 Celsius) and hence output will have units of e-/sec/pix
  ; parameters I(0 Celsius) and dT determined from fit_dark_doubling_rate()

  ; should work for both scalar and array t_celsius inputs

    I0 = 0.0957
    dT = 6.774

    I = I0*2^(t_celsius/dT)

    return, I
end

function nominal_pixel_area, extname

; return output in units of square asec !!

  check_valid_extname, extname

  astr = get_nominal_astrometry(extname)

  cd = astr.cd

  w = where(cd NE 0, nw)
  if nw NE 2 then stop

  area_sq_asec = abs((cd[w[0]]*3600.0d)*(cd[w[1]]*3600.0d))

  return, area_sq_asec

end

function sky_mag_to_e_per_s, mag_per_sq_asec, extname

; input should be sky brightness in mag per sq asec
; returns electrons per second **per pixel**

  if n_elements(mag_per_sq_asec) NE 1 then stop

  par = ci_par_struc()
  e_per_s_per_sq_asec = 10^((par.nominal_zeropoint-mag_per_sq_asec)/2.5)

  e_per_s_per_pixel = e_per_s_per_sq_asec*nominal_pixel_area(extname)

  return, e_per_s_per_pixel
end

function get_gain, extname

; e-/ADU

  check_valid_extname, extname

  case extname of
      'CIW': gain = 1.64
      'CIS': gain = 1.65
      'CIC': gain = 1.64
      'CIN': gain = 1.61
      'CIE': gain = 1.67
  endcase

  return, gain

end

function sky_mag_to_adu_per_s, mag_per_sq_asec, extname

  check_valid_extname, extname

  e_per_s_per_pixel = sky_mag_to_e_per_s(mag_per_sq_asec, extname)

  gain = get_gain(extname)

  adu_per_s_per_pixel = e_per_s_per_pixel/gain

  return, adu_per_s_per_pixel
end

function get_readnoise_electrons, extname

; readnoise in electrons per pixel

  check_valid_extname, extname

  case extname of
      'CIW': readnoise = 12.8
      'CIS': readnoise = 13.3
      'CIC': readnoise = 13.7
      'CIN': readnoise = 13.2
      'CIE': readnoise = 14.2
  endcase

  return, readnoise

end

pro check_valid_extname, extname

  if n_elements(extname) NE 1 then stop
  if size(extname, /type) NE 7 then stop

  par = ci_par_struc()
  extnames = par.ci_extnames

  if total(extnames EQ extname) NE 1 then stop

end

pro _cache_ci_bias

  COMMON _CI_BIAS, bias_cie, bias_cin, bias_cic, bias_cis, bias_ciw 
  if n_elements(bias_cie) EQ 0 then begin
      fname_bias = concat_dir(getenv('CI_REDUCE_ETC'), $
          'CI_master_bias.fits.gz')
      bias_cin = readfits(fname_bias, h_cin)
      if strtrim(sxpar(h_cin, 'EXTNAME'), 2) NE 'CIN' then stop
      bias_ciw = readfits(fname_bias, h_ciw, ex=1)
      if strtrim(sxpar(h_ciw, 'EXTNAME'), 2) NE 'CIW' then stop
      bias_cic = readfits(fname_bias, h_cic, ex=2)
      if strtrim(sxpar(h_cic, 'EXTNAME'), 2) NE 'CIC' then stop
      bias_cie = readfits(fname_bias, h_cie, ex=3)
      if strtrim(sxpar(h_cie, 'EXTNAME'), 2) NE 'CIE' then stop
      bias_cis = readfits(fname_bias, h_cis, ex=4)
      if strtrim(sxpar(h_cis, 'EXTNAME'), 2) NE 'CIS' then stop
  endif

end

function _get_ci_bias, extname

  check_valid_extname, extname

  _cache_ci_bias
  COMMON _CI_BIAS, bias_cie, bias_cin, bias_cic, bias_cis, bias_ciw
  case extname of
      'CIE' : bias = bias_cie
      'CIN' : bias = bias_cin
      'CIC' : bias = bias_cic
      'CIS' : bias = bias_cis
      'CIW' : bias = bias_ciw
  endcase

  return, bias
end

pro _cache_ci_flat

  COMMON _CI_FLAT, flat_cie, flat_cin, flat_cic, flat_cis, flat_ciw 
  if n_elements(flat_cie) EQ 0 then begin
      fname_flat = concat_dir(getenv('CI_REDUCE_ETC'), $
          'CI_master_flat.fits.gz')
      flat_cin = readfits(fname_flat, h_cin)
      if strtrim(sxpar(h_cin, 'EXTNAME'), 2) NE 'CIN' then stop
      flat_ciw = readfits(fname_flat, h_ciw, ex=1)
      if strtrim(sxpar(h_ciw, 'EXTNAME'), 2) NE 'CIW' then stop
      flat_cic = readfits(fname_flat, h_cic, ex=2)
      if strtrim(sxpar(h_cic, 'EXTNAME'), 2) NE 'CIC' then stop
      flat_cie = readfits(fname_flat, h_cie, ex=3)
      if strtrim(sxpar(h_cie, 'EXTNAME'), 2) NE 'CIE' then stop
      flat_cis = readfits(fname_flat, h_cis, ex=4)
      if strtrim(sxpar(h_cis, 'EXTNAME'), 2) NE 'CIS' then stop
  endif

end

function _get_ci_flat, extname

  check_valid_extname, extname

  _cache_ci_flat
  COMMON _CI_FLAT, flat_cie, flat_cin, flat_cic, flat_cis, flat_ciw
  case extname of
      'CIE' : flat = flat_cie
      'CIN' : flat = flat_cin
      'CIC' : flat = flat_cic
      'CIS' : flat = flat_cis
      'CIW' : flat = flat_ciw
  endcase

  return, flat

end

function ci_header_1extname, extname, im, acttime, t_celsius, $
                             primary=primary, sky_mag=sky_mag, seed=seed

  extend = keyword_set(primary)

  mkhdr, header, im, /IMAGE, extend=extend

  sxaddpar, header, 'EXTNAME', extname, 'CI camera name'
  sxaddpar, header, 'CAMTEMP', t_celsius, 'T (deg Celsius)'
  sxaddpar, header, 'ACTTIME', acttime, 'actual exposure time'

; this would be problematic if sky_mag were set to 0.0, although that
; seems like an input unlikely to ever be specified

; obviously, the real raw data won't come with a sky magnitude
; specified in its header, but this may be useful for debugging
  if keyword_set(sky_mag) then $
      sxaddpar, header, 'SKYMAG', sky_mag, 'sky mag per sq asec AB'

  if keyword_set(seed) then $
      sxaddpar, header, 'SEED', seed, 'input seed for random generation'

  return, header
end

function ci_sim_1extname, extname, sky_mag=sky_mag, acttime=acttime, $
                          t_celsius=t_celsius, seed=seed

; sky_mag should be **mags per sq asec**

  if ~keyword_set(sky_mag) then sky_mag = 20.6
; DESI-2549, IN.CI-91010 "Assuming nominal 5 sec exposures"
  if ~keyword_set(acttime) then acttime = 5.0
  if ~keyword_set(t_celsius) then t_celsius = 10.0

  check_valid_extname, extname

  par = ci_par_struc()

  im_electrons = fltarr(par.width, par.height)
; simulate in electrons, then convert to ADU at the very end
; add bias*gain
  bias = _get_ci_bias(extname)
  gain = get_gain(extname)

  im_electrons += bias*gain

; add readnoise
  readnoise_electrons = get_readnoise_electrons(extname)

  im_readnoise_electrons = $
      readnoise_electrons*randomn(seed, size(im_electrons, /dim))

  im_electrons += im_readnoise_electrons

; add constant sky multiplied by flat field
;     for now don't add in poisson noise associated with sky
  sky_e_per_pix = $
      sky_mag_to_e_per_s(sky_mag, extname)*_get_ci_flat(extname)*acttime

  im_electrons += sky_e_per_pix

; add in dark current at some point
;    for now don't add in poisson noise associated with dark current
  total_dark_current_e = dark_current_rate(t_celsius)*acttime ; for now a scalar

  im_electrons += total_dark_current_e

; this isn't formally the precisely right thing to do but w/e
  poisson_noise_e = sqrt(sky_e_per_pix + total_dark_current_e)*randomn(seed, $
      size(im_electrons, /dim))

  im_electrons += poisson_noise_e

; convert to ADU by dividing by the gain
  gain = get_gain(extname)
  im_adu = im_electrons/gain

; convert the image to the right type of integer
  im_adu = uint(round(im_adu)) ; round to avoid biasing low

  return, im_adu
end

; presumably i'll want to have a way of inputting the (ra, dec) once
; i start adding in actual sources ...
pro ci_sim, outname, sky_mag=sky_mag, acttime=acttime, t_celsius=t_celsius, $
            seed=seed

; sky_mag should be **mags per sq asec**

  if ~keyword_set(outname) then stop
  if file_test(outname) then stop ; don't overwrite anything

  if ~keyword_set(sky_mag) then sky_mag = 20.6
; DESI-2549, IN.CI-91010 "Assuming nominal 5 sec exposures"
  if ~keyword_set(acttime) then acttime = 5.0
  if ~keyword_set(t_celsius) then t_celsius = 10.0  

; store seed's original value in order to store it in the output image
; headers for debugging purposes
  if keyword_set(seed) then _seed = seed

  par = ci_par_struc()

  for i=0L, n_elements(par.ci_extnames)-1 do begin
      extname = (par.ci_extnames)[i]
      print, 'Working on ' + extname
      print, sky_mag, acttime, t_celsius
      im = ci_sim_1extname(extname, sky_mag=sky_mag, acttime=acttime, $
                           t_celsius=t_celsius, seed=seed)

      primary = (i EQ 0)
      h = ci_header_1extname(extname, im, acttime, t_celsius, $
                             primary=primary, sky_mag=sky_mag, seed=_seed)
      print, transpose(h)
      writefits, outname, im, h, append=(~primary)
  endfor

end
