pro _cache_ci_bias

  COMMON _CI_BIAS, bias_cie, bias_cin, bias_cic, bias_cis, bias_ciw 
  if n_elements(bias_cie) EQ 0 then begin
      fname_bias = concat_dir(getenv('CI_REDUCE_ETC'), 'CI_master_bias.fits')
      bias_cin = readfits(fname_bias)
      bias_ciw = readfits(fname_bias, ex=1)
      bias_cic = readfits(fname_bias, ex=2)
      bias_cie = readfits(fname_bias, ex=3)
      bias_cis = readfits(fname_bias, ex=4)
  endif

end

function _get_ci_bias, extname

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
      fname_flat = concat_dir(getenv('CI_REDUCE_ETC'), 'CI_master_flat.fits')
      flat_cin = readfits(fname_flat)
      flat_ciw = readfits(fname_flat, ex=1)
      flat_cic = readfits(fname_flat, ex=2)
      flat_cie = readfits(fname_flat, ex=3)
      flat_cis = readfits(fname_flat, ex=4)
  endif

end

function _get_ci_flat

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

pro ci_sim_1extname

; simulate in electrons, then convert to ADU at the very end

; add bias*gain
; add readnoise
; add constant sky multiplied by flat field
;     for now don't add in poisson noise associated with sky


end

; presumably i'll want to have a way of inputting the (ra, dec) once
; i start adding in actual sources ...
pro ci_sim, outname, sky_mag=sky_mag, acttime=acttime, t_celsius=t_celsius

  if ~keyword_set(sky_mag) then sky_mag = 20.3
; DESI-2549, IN.CI-91010 "Assuming nominal 5 sec exposures"
  if ~keyword_set(acttime) then acttime = 5.0
  if ~keyword_set(t_celsius) then t_celsius = 10.0

  

end
