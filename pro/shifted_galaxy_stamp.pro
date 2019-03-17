function shifted_galaxy_stamp, leda_row, sidelen=sidelen, fwhm_pix=fwhm_pix

; sidelen meant to be optional output

  sidelen = long(round((leda_row.d25/0.1)/0.4))

  sidelen += ((sidelen MOD 2) EQ 0)

; put major axis along +Y to start out with
  galaxy = psf_gaussian(npixel=sidelen, $
      fwhm=[leda_row.d25*leda_row.ba, leda_row.d25])

; the -1.0 is necessary because rot() apparently rotates clockwise
  galaxy = rot(galaxy, -1.0*leda_row.pa, /interp)

  case leda_row.extname of
      'CIC' : begin
                  ; flip both vertical and horizontal but no rotation
                  galaxy = reverse(galaxy, 1)
                  galaxy = reverse(galaxy, 2)
              end
      'CIS' : 
      'CIN' : begin
                  ; flip both vertical and horizontal but no rotation
                  galaxy = reverse(galaxy, 1)
                  galaxy = reverse(galaxy, 2)
              end
      'CIE' : galaxy = rotate(galaxy, 1)
      'CIW' : galaxy = rotate(galaxy, 3)
  endcase

; then do the shifting

  ix = long(round(leda_row.x))
  iy = long(round(leda_row.y))
  galaxy = sshift2d(galaxy, [leda_row.x-ix, leda_row.y-iy])

  return, galaxy

end
