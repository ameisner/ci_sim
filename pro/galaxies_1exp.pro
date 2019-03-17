pro _cache_leda

  COMMON _LEDA, leda
  if n_elements(leda) EQ 0 then begin
      fname = $
          '/project/projectdirs/cosmo/work/wise/etc/leda-logd25-0.05.fits.gz'
      leda = mrdfits(fname, 1)
      leda[where(leda.ba EQ -999)].ba = 1.0 ; circular
      leda[where(leda.pa EQ -999)].pa = 0.0 ; ??
      leda = leda[where(leda.d25 LE 60)] ; 1 arcmin
      leda = leda[where((leda.bmag NE -999) OR (leda.imag NE -999))]
      mag = (((leda.bmag NE -999)*leda.bmag + $
              (leda.imag NE -999)*leda.imag))/((leda.bmag NE -999) + $
                                               (leda.imag NE -999))
      addstr = replicate({mag_ab: 0.0}, n_elements(mag))
      addstr.mag_ab = mag
      leda = struct_addtags(leda, addstr)
  end

end

function galaxies_1exp, racen, deccen

  ang_max = 1.7 ; degrees

;  racen = astr.crval[0]
;  deccen = astr.crval[1]

  _cache_leda

  COMMON _LEDA, leda

  w = where((leda.dec GT (deccen-ang_max)) AND $
            (leda.dec LT (deccen+ang_max)), nw)

  if nw EQ 0 then return, -1

  leda_ = leda[w]

  dangle = djs_diff_angle(leda_.ra, leda_.dec, racen, deccen)

  w = where(dangle LT ang_max, nw)

  if nw EQ 0 then return, -1

  leda_ = leda_[w]

  return, leda_
end
