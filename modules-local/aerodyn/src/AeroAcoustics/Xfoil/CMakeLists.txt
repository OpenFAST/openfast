
add_library(RunXfoil  dplot_noise.f get_airfoil_coords.f naca.f  plutil_noise.f  profil.f  sort.f  spline.f  userio.f xbl.f  xblsys.f xfoil_noise.f xgdes_noise.f  xgeom_noise.f xoper_noise.f  xpanel.f xsolve.f xtcam_noise.f xutils.f xfoil_noise_mods.f90 XFOIL.INC XBL.INC CIRCLE.INC PINDEX.INC)

target_include_directories (RunXfoil PUBLIC)
