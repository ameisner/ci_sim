pro _cache_nominal_astrometry

  COMMON _CI_ASTROM, astrom
  if n_elements(astrom) EQ 0 then begin
      fname = concat_dir(getenv('CI_REDUCE_ETC'), $
          'viewer_astrom_index-as_designed.bigtan.fits')
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

function ci_bdy_coords

  par = ci_par_struc()

  x_bottom = dindgen(par.width+1) - 0.5d
  y_bottom = dblarr(par.width+1) - 0.5d

  x_top = dindgen(par.width+1) - 0.5d
  y_top = dblarr(par.width+1) + (par.height - 1.0d) + 0.5d

  x_left = dblarr(par.height+1) - 0.5d
  y_left = dindgen(par.height+1) - 0.5d

  x_right = dblarr(par.height+1) + (par.width - 1.0d) + 0.5d
  y_right = dindgen(par.height+1) - 0.5

  x_bdy = [x_bottom, x_right, reverse(x_top), reverse(x_left)]
  y_bdy = [y_bottom, y_right, reverse(y_top), reverse(y_left)]

  outstr = replicate({x_bdy: 0.0d, y_bdy: 0.0d}, n_elements(x_bdy))
  outstr.x_bdy = x_bdy
  outstr.y_bdy = y_bdy

  return, outstr
end

function _wrap_ra, ra
  return, ra + (360.0d)*(ra LT 0)
end

pro distance_to_neighbor_source, cat, astr

; want to have a metric related to isolation of the injected sources
; so that downstream tests can conveniently ignore blends when 
; that is appropriate/useful

; calculate distance to nearest neighbor in asec and append that
; as a column in input cat (which thus gets modified)

  if size(cat, /type) NE 8 then stop
  if size(astr, /type) NE 8 then stop

  xy2ad, cat.x, cat.y, astr, ra, dec

  ang_max = 1.0d/60.0d ; deg (this is 1 arcmin)
  spherematch, _wrap_ra(ra), dec, _wrap_ra(ra), dec, ang_max, m, _m, dist, $
      maxmatch=0

  w = where(m NE _m, nw) ; don't care about self-matches

  if nw EQ 0 then begin
      addstr = replicate({nearest_neighbor_asec: -1.0}, n_elements(cat))
      cat = struct_addtags(cat, addstr)
      return
  endif

  m = m[w]
  _m = _m[w]
  dist = dist[w]

  sind = sort(m)

  m = m[sind]
  _m = _m[sind]
  dist = dist[sind]

  ind_bdy = uniq(m)

  n_bdy = n_elements(ind_bdy)

; -1.0 will be dummy value indicating no neighbors
  d_nearest_asec = fltarr(n_elements(cat)) - 1.0
  for i=0L, n_bdy-1 do begin
      ind_l = (i EQ 0) ? 0 : (ind_bdy[i-1] + 1)
      ind_u = ind_bdy[i]
      d_nearest_asec[m[ind_bdy[i]]] = min(dist[ind_l:ind_u])*3600.0
  endfor

  addstr = replicate({nearest_neighbor_asec: 0.0}, n_elements(cat))
  addstr.nearest_neighbor_asec = d_nearest_asec
  cat = struct_addtags(cat, addstr)

end

function mjd_lst_lookup_table

; 07:00 UTC on 1 Jan 2019
  mjd0 = 58484.29166667d
  mjds = mjd0 + dindgen(365)

  jds = mjds + 2400000.5d

  longitude = -111.59989d ; KPNO deg E copied from DESI-4000
  timezone = -7 ; KPNO

  lsts = dblarr(n_elements(jds))
  for i=0L, n_elements(jds)-1 do begin
      ct2lst, lst, longitude, timezone, jds[i]
      lsts[i] = lst ; what are the units of lst returned ???
  endfor

  outstr = replicate({mjd: 0.0d, jd: 0.0d, lst: 0.0d}, n_elements(jds))
  outstr.mjd = mjds
  outstr.jd = jds
  outstr.lst = lsts*15.0d ; put LST in degrees rather than hours

  return, outstr
end

pro _cache_lst

  COMMON LST_LOOKUP, lst_str
  if n_elements(lst_str) EQ 0 then begin
      lst_str = mjd_lst_lookup_table()
  endif

end

function ra_to_mjd, ra, lst=lst

  if n_elements(ra) NE 1 then stop

  _cache_lst
  COMMON LST_LOOKUP, lst_str
  _ra = _wrap_ra(ra)

  _ = min(abs(lst_str.lst - _ra), indmin)

  lst = lst_str[indmin].lst ; optional output

  return, lst_str[indmin].mjd
end

function psf_stamp_size

; as a function of how bright a source is, return the 
; necessary PSF stamp size in pixels

  return, -1
end

function scale_psf, stamp, flux

  if total(stamp LT -1.0d-8) NE 0 then stop

  return, (stamp/total(stamp))*flux

end

function centered_psf_stamp, fwhm_pix, sidelen=sidelen

; sidelen keyword argument intended to be optional input

  if n_elements(fwhm_pix) NE 2 then stop

; demand that sidelength always be *odd* integer

  if ~keyword_set(sidelen) then $
      sidelen = long(round(max(fwhm_pix)*11)) + $
              ((long(round(max(fwhm_pix)*11)) MOD 2) EQ 0)

  if n_elements(sidelen) NE 1 then stop
  if round(sidelen) NE sidelen then stop
  if (sidelen MOD 2) EQ 0 then stop

  psf = psf_gaussian(npixel=[sidelen, sidelen], fwhm=fwhm_pix, /normalize)

  return, psf

end

function shifted_psf_stamp, fwhm_pix, frac_shift, sidelen=sidelen

  if n_elements(frac_shift) NE 2 then stop
  if total(abs(frac_shift) GT 0.5) NE 0 then stop

  centered_psf = centered_psf_stamp(fwhm_pix, sidelen=sidelen)

  shifted_psf = sshift2d(centered_psf, frac_shift)
; renormalize when done shifting
  shifted_psf = scale_psf(shifted_psf, 1.0)

  return, shifted_psf

end

pro _mwr_dummy, fname, exptime=exptime, ra=ra, dec=dec

    openw, 2, fname, /swap_if_little

    compile_opt idl2,hidden

    fxaddpar, header, 'SIMPLE', 'T','Dummy Created by MWRFITS'
    fxaddpar, header, 'BITPIX', 16, 'Dummy primary header created by MWRFITS'
    fxaddpar, header, 'NAXIS', 0, 'No data is associated with this header'
    fxaddpar, header, 'EXTEND', 'T', 'Extensions may be present'

    fxaddpar, header, 'EXTNAME', 'CI'
    fxaddpar, header, 'OBSERVAT', 'KPNO', 'Observatory name'
    fxaddpar, header, 'OBS-LAT', 31.96403, '[deg] Observatory latitude'
    fxaddpar, header, 'OBS-LONG', -111.59989, '[deg] Observatory east longitude'

    fxaddpar, header, 'OBS-ELEV', 2097, '[m] Observatory elevation'
    fxaddpar, header, 'TELESCOP', 'KPNO 4.0m telescope', 'Telescope name'
    fxaddpar, header, 'INSTRUME', 'CI', 'Instrument name'

    fxaddpar, header, 'PROGRAM', 'commiss', 'Program name'
    fxaddpar, header, 'FLAVOR', 'science', 'Observation type'

    if arg_present(exptime) then $
    fxaddpar, header, 'EXPTIME', exptime, '[s] Actual exposure time'

    if arg_present(ra) then $
    fxaddpar, header, 'SKYRA', ra, $
        '[deg] Telescope right ascension (pointing on sky)'
    fxaddpar, header, 'SKYDEC', dec, $
        '[deg] Telescope declination (pointing on sky)'

    COMMON _ci_frame_info, _tileid, _pass, _expid, _time_specific
    
    if n_elements(_tileid) EQ 1 then $
        fxaddpar, header, 'TILEID', _tileid

    if n_elements(_pass) EQ 1 then $
        fxaddpar, header, 'PASS', _pass

    if n_elements(_expid) EQ 1 then $
        fxaddpar, header, 'EXPID', _expid, 'Exposure number'

    if n_elements(_time_specific) NE 0 then begin
        fxaddpar, header, 'MJD-OBS', _time_specific.mjd, $
            'Modified Julian Date of observation'
        fxaddpar, header, 'AIRMASS', _time_specific.airmass, $
            'Airmass'
        fxaddpar, header, 'NIGHT', _time_specific.night, $
            'Observing night'
        fxaddpar, header, 'DATE-OBS', _time_specific.date_obs, $
            '[UTC] Observation date start time'
        fxaddpar, header, 'MOONRA', _time_specific.moonra, $
            '[deg] Moon RA at start of exposure'
        fxaddpar, header, 'MOONDEC', _time_specific.moondec, $
            '[deg] Moon declination at start of exposure'
        fxaddpar, header, 'MOONANGL', _time_specific.moon_dist_deg, $
            '[deg] Distance to moon at start of exposure'
    endif

    mwr_header, 2, header

    close, 2

end

function fwhm_asec_to_pix, fwhm_asec, extname, force_symmetric=force_symmetric

; force_symmetric is optional boolean argument meant to allow
; forcing of PSF to be symmetric in terms of fwhm_x, fwhm_y - this
; may be useful for downstream testing purposes

; note that there will be an x direction fwhm in pixels
; and also a y direction fwhm in pixels

  if n_elements(fwhm_asec) NE 1 then stop

  check_valid_extname, extname

  astr = get_nominal_astrometry(extname)
  
  w = where(astr.cd NE 0)


  platescale_x = abs((astr.cd)[w[0]])*3600.0d ; asec/pix
  platescale_y = abs((astr.cd)[w[1]])*3600.0d ; asec/pix

  fwhm_x = fwhm_asec/platescale_x
  fwhm_y = fwhm_asec/platescale_y

  if keyword_set(force_symmetric) then begin
      fwhm_mean = (fwhm_x + fwhm_y)/2
      fwhm_x = fwhm_mean
      fwhm_y = fwhm_mean
  endif

  return, [fwhm_x, fwhm_y]
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
      'CIE': gain = 1.64
      'CIN': gain = 1.65
      'CIC': gain = 1.64
      'CIS': gain = 1.61
      'CIW': gain = 1.67
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
      'CIE': readnoise = 12.8
      'CIN': readnoise = 13.3
      'CIC': readnoise = 13.7
      'CIS': readnoise = 13.2
      'CIW': readnoise = 14.2
  endcase

  return, readnoise

end

function ci_min_edge_dist, cat

  par = ci_par_struc()

  min_edge_dist = 10000L

  min_edge_dist <= (cat.x+0.5)
  min_edge_dist <= (cat.y+0.5)
  min_edge_dist <= (par.width - 0.5 - cat.x)
  min_edge_dist <= (par.height - 0.5 - cat.y)

  return, min_edge_dist

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

pro _cache_desi_tiles

  COMMON _DESI_TILES, desi_tiles
  if n_elements(desi_tiles) EQ 0 then begin
      fname_desi_tiles = concat_dir(getenv('CI_REDUCE_ETC'), 'desi-tiles.fits')
      desi_tiles = mrdfits(fname_desi_tiles, 1)
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

function sources_only_image, fwhm_pix, acttime, astr, $
                             do_gaia_sources=do_gaia_sources, $
                             source_catalog=cat

; source_catalog keyword argument meant to be an optional output

  par = ci_par_struc()

  if ~keyword_set(do_gaia_sources) then begin
      ; for initial testing purposes, add just two well-separated
      ; moderately bright sources near the middle of the image
      cat = replicate({x: 0.0d, y: 0.0d, mag_ab: 0.0, extname: astr.extname}, 2)
      cat[0].x = 1500.0d & cat[0].y = 1000.0d & cat[0].mag_ab = 17.0
      cat[1].x = 1800.0d & cat[1].y = 1200.0d & cat[1].mag_ab = 18.0
  endif else begin
      cat = gaia_sources_1cam(astr)
      print, 'Adding ', n_elements(cat), ' Gaia sources'
      addstr = replicate({mag_ab: 0.0, extname: astr.extname}, n_elements(cat))
      addstr.mag_ab = cat.phot_g_mean_mag
      cat = struct_addtags(cat, addstr)
  endelse
; construct this image in units of electrons !!!

  im_electrons = fltarr(par.width, par.height)

  ; this is the total flux in (detected)
  ; electrons of the entire source if its entire profile were to
  ; "land" within the CI image's boundaries, and assuming
  ; the flat field is perfectly uniform across the CCD

  ; this can be different from the amount of this source's flux that
  ; actually lands inside of the CI image's boundaries if
  ; the source is near an edge

  total_flux_electrons = $
      (10^((par.nominal_zeropoint - cat.mag_ab)/2.5))*acttime

  for i=0L, n_elements(cat)-1 do begin
      ix = long(round(cat[i].x))
      iy = long(round(cat[i].y))

      frac_shift = [cat[i].x-ix, cat[i].y-iy]

      psf = shifted_psf_stamp(fwhm_pix, frac_shift, sidelen=sidelen)

      psf = scale_psf(psf, total_flux_electrons[i])

      half = long(sidelen)/2

      xmin_int = ix - half
      xmax_int = ix + half

      ymin_int = iy - half
      ymax_int = iy + half

      xmin_add = (xmin_int > 0)
      xmax_add = (xmax_int < (par.width-1))

      ymin_add = (ymin_int > 0)
      ymax_add = (ymax_int < (par.height-1))

      im_electrons[xmin_add:xmax_add, ymin_add:ymax_add] += $
          psf[(xmin_add-xmin_int):(sidelen-(xmax_int-xmax_add)-1), $
              (ymin_add-ymin_int):(sidelen-(ymax_int-ymax_add)-1)] 
  endfor

  delvarx, sidelen

  total_flux_adu = total_flux_electrons/get_gain(astr.extname)

  addstr = replicate({total_flux_adu_flatfielded: 0.0, $
                      min_edge_dist: 0.0}, n_elements(cat))
  addstr.total_flux_adu_flatfielded = total_flux_adu
  addstr.min_edge_dist = ci_min_edge_dist(cat)

  cat = struct_addtags(cat, addstr)

; add extra column to catalog with distance to its nearest neighbor
  distance_to_neighbor_source, cat, astr

  return, im_electrons
end

function ci_header_1extname, extname, im, acttime, t_celsius, $
                             primary=primary, sky_mag=sky_mag, seed=seed, $
                             telra=telra, teldec=teldec, fwhm_pix=fwhm_pix, $
                             astr=astr, fwhm_asec=fwhm_asec

; make astr a required argument rather than an optional keyword argument?

  check_valid_extname, extname

  if ~keyword_set(telra) then telra = 0.0d
  if ~keyword_set(teldec) then teldec = 0.0d

  extend = keyword_set(primary)

  mkhdr, header, im, /IMAGE, extend=extend
  putast, header, astr

  sxaddpar, header, 'EXTNAME', extname, 'Extension name' ; CI camera name
  sxaddpar, header, 'CAMTEMP', t_celsius, '[deg] Camera temperature' ; Celsius
  sxaddpar, header, 'EXPTIME', acttime, '[s] Actual exposure time'
  sxaddpar, header, 'REQRA', telra, $
      '[deg] Requested right ascension (observer input)'
  sxaddpar, header, 'REQDEC', teldec, $
      '[deg] Requested declination (observer input)'
  sxaddpar, header, 'TARGTRA', telra, $
      '[deg] Target right ascension (to TCS)'
  sxaddpar, header, 'TARGTDEC', teldec, $
      '[deg] Target declination (to TCS)'
  sxaddpar, header, 'SKYRA', telra, $
      '[deg] on sky right ascension (from TCS) '
  sxaddpar, header, 'SKYDEC', teldec, $
      '[deg] on sky declination (from TCS)'

; this would be problematic if sky_mag were set to 0.0, although that
; seems like an input unlikely to ever be specified

; obviously, the real raw data won't come with a sky magnitude
; specified in its header, but this may be useful for debugging
  if keyword_set(sky_mag) then $
      sxaddpar, header, 'SKYMAG', sky_mag, 'sky mag per sq asec AB'

  if keyword_set(seed) then $
      sxaddpar, header, 'SEED', seed, 'input seed for random generation'

  if keyword_set(fwhm_pix) then begin
      sxaddpar, header, 'FWHMPX', fwhm_pix[0], 'x FWHM in pixels'
      sxaddpar, header, 'FWHMPY', fwhm_pix[1], 'y FWHM in pixels'
  endif

  if keyword_set(fwhm_asec) then begin
      sxaddpar, header, 'FWHMASEC', fwhm_asec, 'FWHM in asec'
  endif

  sxdelpar, header, 'HISTORY'

  return, header
end

function ci_sim_1extname, extname, sky_mag=sky_mag, acttime=acttime, $
                          t_celsius=t_celsius, seed=seed, $
                          fwhm_asec=fwhm_asec, fwhm_pix=fwhm_pix, $
                          astr=astr, do_gaia_sources=do_gaia_sources, $
                          source_catalog=source_catalog, $
                          force_symmetric=force_symmetric, $
                          galaxy_catalog=galaxy_catalog

; fwhm_pix meant to be used as optional **output**
; sky_mag should be **mags per sq asec**
; source_catalog meant to be used as optional output
; galaxy_catalog meant to be used as optional output

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
  flat = _get_ci_flat(extname)
  sky_e_per_pix = $
      sky_mag_to_e_per_s(sky_mag, extname)*flat*acttime

  im_electrons += sky_e_per_pix

; add in dark current at some point
;    for now don't add in poisson noise associated with dark current
  total_dark_current_e = dark_current_rate(t_celsius)*acttime ; for now a scalar

  im_electrons += total_dark_current_e

  fwhm_pix = fwhm_asec_to_pix(fwhm_asec, extname, $
      force_symmetric=force_symmetric)
  im_sources_e = sources_only_image(fwhm_pix, acttime, $
      astr, do_gaia_sources=do_gaia_sources, $ 
      source_catalog=source_catalog)*flat

  galaxy_catalog = galaxies_1cam(astr)
help,galaxy_catalog

  im_electrons += im_sources_e
; this isn't formally the precisely right thing to do but w/e
  poisson_noise_e = sqrt(sky_e_per_pix + total_dark_current_e + $
      im_sources_e)*randomn(seed, $
      size(im_electrons, /dim))

  im_electrons += poisson_noise_e

; convert to ADU by dividing by the gain
  gain = get_gain(extname)
  im_adu = im_electrons/gain

; convert the image to the right type of integer
  im_adu = uint(round(im_adu) < 65535) ; round to avoid biasing low

  return, im_adu
end

pro ci_sim, outname, telra=telra, teldec=teldec, sky_mag=sky_mag, $
            acttime=acttime, t_celsius=t_celsius, seed=seed, $
            fwhm_asec=fwhm_asec, do_gaia_sources=do_gaia_sources, $
            dummy_ext=dummy_ext, force_symmetric=force_symmetric

; sky_mag should be **mags per sq asec**

  if ~keyword_set(outname) then stop
  if file_test(outname) then stop ; don't overwrite anything

  if ~keyword_set(sky_mag) then sky_mag = 20.6
; DESI-2549, IN.CI-91010 "Assuming nominal 5 sec exposures"
  if ~keyword_set(acttime) then acttime = 5.0
  if ~keyword_set(t_celsius) then t_celsius = 10.0
  if ~keyword_set(telra) then telra = 0.0d
  if ~keyword_set(teldec) then teldec = 0.0d  
; my guess at median r band seeing with Mayall ...
  if ~keyword_set(fwhm_asec) then fwhm_asec = 1.25

; store seed's original value in order to store it in the output image
; headers for debugging purposes
  if keyword_set(seed) then _seed = seed

  par = ci_par_struc()

  if (n_elements(sky_mag) NE 1) AND $
     (n_elements(sky_mag) NE n_elements(par.ci_extnames)) then stop

  sky_mag_array = (n_elements(sky_mag) EQ n_elements(par.ci_extnames))

  if keyword_set(dummy_ext) then _mwr_dummy, outname, exptime=acttime, $
      ra=telra, dec=teldec

  for i=0L, n_elements(par.ci_extnames)-1 do begin
      extname = (par.ci_extnames)[i]
      print, 'Working on ' + extname
      if sky_mag_array then _sky_mag = sky_mag[i] else _sky_mag = sky_mag
      print, _sky_mag, acttime, t_celsius

      astr = get_astrometry_radec(telra, teldec, extname)

      im = ci_sim_1extname(extname, sky_mag=_sky_mag, acttime=acttime, $
                           t_celsius=t_celsius, seed=seed, $
                           fwhm_asec=fwhm_asec, fwhm_pix=fwhm_pix,astr=astr, $
                           do_gaia_sources=do_gaia_sources, $
                           source_catalog=source_catalog, $
                           force_symmetric=force_symmetric, $
                           galaxy_catalog=galaxy_catalog)

      primary = (i EQ 0)
      h = ci_header_1extname(extname, im, acttime, t_celsius, $
                             primary=primary, sky_mag=_sky_mag, seed=_seed, $
                             telra=telra, teldec=teldec, fwhm_pix=fwhm_pix, $
                             astr=astr, fwhm_asec=fwhm_asec)
      print, transpose(h)
      writefits, outname, im, h, append=(~primary OR keyword_set(dummy_ext))

; deal with the source catalog if there is one
      if size(source_catalog, /type) EQ 8 then begin
          if ~keyword_set(joint_catalog) then $
              joint_catalog = source_catalog else $
              joint_catalog = struct_append(joint_catalog, source_catalog)
      endif
      if size(galaxy_catalog, /type) EQ 8 then begin
          if ~keyword_set(joint_g_catalog) then $
              joint_g_catalog = galaxy_catalog else $
              joint_g_catalog = struct_append(joint_g_catalog, galaxy_catalog)
      endif
      delvarx, source_catalog
  endfor

; append one extra extension for debugging purposes only, containing
; the full table of injected sources including all CI cameras (rather
; than have one such extension per CI camera)
  if size(joint_catalog, /type) EQ 8 then begin
      outname_cat = repstr(outname, '.fits', '.cat.fits')
      if file_test(outname_cat) then stop
      mwrfits, joint_catalog, outname_cat
      if size(joint_g_catalog, /type) EQ 8 then begin
          mwrfits, joint_g_catalog, outname_cat
      endif
  endif

end

pro sim_desi_pointing, desi_tiles_row, outdir=outdir, dummy_ext=dummy_ext, $
                       force_symmetric=force_symmetric, expnum=expnum

; desi_tiles_row should be one row of the desi-tiles.fits table

  if ~keyword_set(outdir) then $
   outdir = '/project/projectdirs/desi/users/ameisner/CI/ci_data_challenge/sims'

  if size(desi_tiles_row, /type) NE 8 then stop
  if n_elements(desi_tiles_row) NE 1 then stop

  if (desi_tiles_row.pass NE 0) OR ~desi_tiles_row.in_desi then stop

  telra = desi_tiles_row.ra
  teldec = desi_tiles_row.dec

  if ~keyword_set(expnum) then expnum = desi_tiles_row.tileid

  outname = 'desi-' + string(expnum, format='(I06)') + '.fits'

  outname = concat_dir(outdir, outname)

  if file_test(outname) then stop

  print, 'Working on: TILEID = ', desi_tiles_row.tileid, $
      ' , pass = ', desi_tiles_row.pass, $
      ' , ra = ', desi_tiles_row.ra, ' , ', desi_tiles_row.dec, $
      ' , output name = ', outname

  seed = long(desi_tiles_row.tileid)
  ci_sim, outname, telra=desi_tiles_row.ra, teldec=desi_tiles_row.dec, $
      sky_mag=sky_mag, acttime=acttime, t_celsius=t_celsius, seed=seed, $
      fwhm_asec=fwhm_asec, /do_gaia_sources, dummy_ext=dummy_ext, $
      force_symmetric=force_symmetric

end

pro sim_desi_pointings, indstart=indstart, nproc=nproc, outdir=outdir, $
                        dummy_ext=dummy_ext, force_symmetric=force_symmetric

; wrapper for sim_desi_pointing

  _cache_desi_tiles
  COMMON _DESI_TILES, all_tiles

  desi_tiles_pass0 = all_tiles[where((all_tiles.pass EQ 0) AND $
                                      all_tiles.in_desi, n_sim)]

; randomize ordering of pass=0 tiles
  seed = 99L
  sind = sort(randomu(seed, n_sim))

  desi_tiles_pass0 = desi_tiles_pass0[sind]

  if ~keyword_set(indstart) then indstart = 0L
  if ~keyword_set(nproc) then nproc = n_sim

  indstart = long(indstart)
  nproc = long(nproc)

  indend = (indstart + nproc - 1) < (n_sim - 1)

  for i=indstart, indend do begin
      sim_desi_pointing, desi_tiles_pass0[i], outdir=outdir, $
          dummy_ext=dummy_ext, force_symmetric=force_symmetric
  endfor

end

pro desi_1pointing, tileid, outdir=outdir, dummy_ext=dummy_ext, $
                    force_symmetric=force_symmetric

; alternative wrapper for sim_desi_pointing

  if n_elements(tileid) NE 1 then stop

  _cache_desi_tiles
  COMMON _DESI_TILES, all_tiles

  w = where(all_tiles.tileid EQ tileid, nw)
  if (nw NE 1) then stop

  sim_desi_pointing, all_tiles[w[0]], outdir=outdir, dummy_ext=dummy_ext, $
      force_symmetric=force_symmetric

end

pro _cache_ci_tiles

  COMMON _CI_TILES, ci_tiles
  if n_elements(ci_tiles) EQ 0 then begin
      ; when this becomes more finalized put a copy in CI_REDUCE_ETC
      fname = '/project/projectdirs/desi/cmx/ci/tiles/v1/ci-tiles-v1.fits'
      ci_tiles = mrdfits(fname, 1)
  endif

end

pro _cache_ci_frame_info, tileid, pass, time_specific

  COMMON _CI_FRAME_INFO, _tileid, _pass, _expid, _time_specific
  if n_elements(tileid) NE 0 then _tileid = tileid
  if n_elements(pass) NE 0 then _pass = pass
  if n_elements(tileid) NE 0 then _expid = tileid + 50000L ; HACK
  if n_elements(time_specific) NE 0 then _time_specific = time_specific

end

function time_specific_metadata, ra, dec

  mjd = ra_to_mjd(ra, lst=lst)

  kpno_lat = 31.96403d
  z = djs_diff_angle(ra, dec, lst, kpno_lat)

  airmass = 1/cos(z/(180.0d/!dpi))

  if (airmass LT 1) then stop

  date_string = strmid(date_conv(mjd +  2400000.5d, 'FITS'), 0, 10)

; subtract half a day to figure out what "night" this MJD belongs to
  night_string = strmid(date_conv(mjd + 2400000.5d - 0.5d, 'FITS'), 0, 10)

; ra_moon, dec_moon are in degrees
  moonpos, mjd + 2400000.5d, ra_moon, dec_moon

  moon_dist = djs_diff_angle(ra_moon, dec_moon, ra, dec)

  outstr = {mjd: mjd, airmass: airmass, date_obs: date_string, $
            night: night_string, moonra: ra_moon, moondec: dec_moon, $
            moon_dist_deg: moon_dist}

  return, outstr
end

; use 47002 as a test case
; .COM mwrfits
pro ci_1pointing, tileid, outdir=outdir, force_symmetric=force_symmetric

  if ~keyword_set(outdir) then outdir = '$SCRATCH/sims'

  dummy_ext = 1 ; HARDCODED !

  if n_elements(tileid) NE 1 then stop

  _cache_ci_tiles
  COMMON _CI_TILES, ci_tiles

  w = where(ci_tiles.tileid EQ tileid, nw)
  if (nw NE 1) then stop

  _cache_ci_frame_info, ci_tiles[w[0]].tileid, ci_tiles[w[0]].pass, $
      time_specific_metadata(ci_tiles[w[0]].ra, ci_tiles[w[0]].dec)

  COMMON _CI_FRAME_INFO, _tileid, _pass, _expid, _time_specific

  sim_desi_pointing, ci_tiles[w[0]], outdir=outdir, dummy_ext=dummy_ext, $
      force_symmetric=force_symmetric, expnum=_expid

end

; .COM mwrfits
pro ci_pointings, indstart, nproc

  _cache_ci_tiles
  COMMON _CI_TILES, ci_tiles

  ci_tiles = ci_tiles[where(ci_tiles.in_desi)]
  ntiles = n_elements(ci_tiles)

  seed = 99L
  sind = sort(randomu(seed, ntiles))
  ci_tiles = ci_tiles[sind]

  indend = (indstart + nproc - 1) < (ntiles - 1)

  for i=indstart, indend do begin
      ci_1pointing, ci_tiles[i].tileid
  endfor

end
