function gaia_sources_1cam, astr

  bdy_pixel_coords = ci_bdy_coords()

  xy2ad, bdy_pixel_coords.x_bdy, bdy_pixel_coords.y_bdy, $
      astr, ra_bdy, dec_bdy

  gaia = read_gaia_cat(ra_bdy, dec_bdy)
 
; for now ignore the issue of stars off edge contributing light into
; the chip, and downselect to the sources with centroids that 
; fall inside the chip

; ignore gaia motions for now

  ad2xy, gaia.ra, gaia.dec, astr, x, y

  par = ci_par_struc()
  good = (x GT -0.5) AND (x LT (par.width-0.5d)) AND $ 
         (y GT -0.5) AND (y LT (par.height-0.5d))

; this should probably never happen ??
  if total(good) EQ 0 then return, -1

  result = gaia[where(good, ngood)]

  x = x[where(good)]
  y = y[where(good)]

  addstr = replicate({x: 0.0d, y: 0.0d}, $
      ngood)

  addstr.x = x
  addstr.y = y

  result = struct_addtags(result, addstr)

  return, result

end
