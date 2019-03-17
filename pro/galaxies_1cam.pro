pro _cache_galaxies_1exp, racen, deccen

  COMMON _GALAXIES_EXP, ra_cache, dec_cache, galaxies
  if (n_elements(galaxies) EQ 0) || $
     ((ra_cache NE racen) OR (dec_cache NE deccen)) then begin
      galaxies = galaxies_1exp(racen, deccen)
      ra_cache = racen
      dec_cache = deccen
  endif

end

function galaxies_1cam, astr

  racen = astr.crval[0]
  deccen = astr.crval[1]

  _cache_galaxies_1exp, racen, deccen
  COMMON _GALAXIES_EXP, _, __, galaxies

  if size(galaxies, /type) NE 8 then return, -1

  ad2xy, galaxies.ra, galaxies.dec, astr, x, y

  good = ((x GT -0.5) AND (x LT 3071.5) AND (y GT -0.5) AND $
          (y LT 2047.5))

  wgood = where(good, nw)
  if nw EQ 0 then return, -1

  g = galaxies[wgood]
  extname = (total(tag_names(astr) EQ 'EXTNAME') GT 0) ? astr.extname : ''
  addstr = replicate({extname: extname, x: 0.0d, y: 0.0d}, n_elements(g))

  addstr.x = x[wgood]
  addstr.y = y[wgood]

  outstr = struct_addtags(g, addstr)
  return, outstr

end
