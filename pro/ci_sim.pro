function ci_par_struc

   par = {width: 3072, height: 2048, $
          nominal_zeropoint: 26.56}

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

pro sky_mag_to_e_per_s, mag_per_sq_asec

; input should be sky brightness in mag per sq asec

  if n_elements(mag_per_sq_asec) NE 1 then stop

  par = ci_par_struc()
  e_per_sq_asec = 10^((par.nominal_zeropoint-mag_per_sq_asec)/2.5)

  print, e_per_sq_asec

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

  extnames = ['CIE', 'CIN', 'CIC', 'CIS', 'CIW']

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

pro ci_header_1extname, extname, acttime=acttime, t_celsius=t_celsius

end

pro ci_sim_1extname, extname, sky_mag=sky_mag, acttime=acttime, $
                     t_celsius=t_celsius

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

  im_readnoise_electrons = readnoise*randomn(seed, size(im_electrons, /dim))

  im_electrons += im_readnoise_electrons

; add constant sky multiplied by flat field
;     for now don't add in poisson noise associated with sky


; add in dark current at some point


; convert to ADU by dividing by the gain !!!
; convert the image to the right type of integer !!!

end

; presumably i'll want to have a way of inputting the (ra, dec) once
; i start adding in actual sources ...
pro ci_sim, outname, sky_mag=sky_mag, acttime=acttime, t_celsius=t_celsius

  if ~keyword_set(sky_mag) then sky_mag = 20.6
; DESI-2549, IN.CI-91010 "Assuming nominal 5 sec exposures"
  if ~keyword_set(acttime) then acttime = 5.0
  if ~keyword_set(t_celsius) then t_celsius = 10.0

  

end
