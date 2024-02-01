SELECT 
	-- Basic information
	f1.object_id, f1.skymap_id, f1.ra, f1.dec,
	f1.a_g, f1.a_r, f1.a_i, f1.a_z, f1.a_y,
	f1.g_extendedness_value, f1.r_extendedness_value, f1.i_extendedness_value, 
	f1.z_extendedness_value, f1.y_extendedness_value,
	m2.g_blendedness_abs, m2.r_blendedness_abs, 
	m2.i_blendedness_abs, m2.z_blendedness_abs, m2.y_blendedness_abs, 
	m2.g_blendedness_flag, m2.r_blendedness_flag, m2.i_blendedness_flag, 
	m2.z_blendedness_flag, m2.y_blendedness_flag, 
	f3.g_apertureflux_10_flux, f3.g_apertureflux_10_flag,
	f3.r_apertureflux_10_flux, f3.r_apertureflux_10_flag, f3.i_apertureflux_10_flux, f3.i_apertureflux_10_flag,
	f3.z_apertureflux_10_flux, f3.z_apertureflux_10_flag, f3.y_apertureflux_10_flux, f3.y_apertureflux_10_flag,

	-- Flags for later selection
	-- 1. Source is outside usable exposure region
	f1.g_pixelflags_edge,
	f1.r_pixelflags_edge,
	f1.i_pixelflags_edge,
	f1.z_pixelflags_edge,
	f1.y_pixelflags_edge,
	
	-- 2. Saturated or interpolated pixels on the footprint (not center)
	f1.g_pixelflags_saturated,
	f1.r_pixelflags_saturated,
	f1.i_pixelflags_saturated,
	f1.z_pixelflags_saturated,
	f1.y_pixelflags_saturated,
	f1.g_pixelflags_interpolated,
	f1.r_pixelflags_interpolated,
	f1.i_pixelflags_interpolated,
	f1.z_pixelflags_interpolated,
	f1.y_pixelflags_interpolated,
	
	-- 3. Center is close to a clipped or suspicious pixel
	f1.g_pixelflags_saturatedcenter,
	f1.r_pixelflags_saturatedcenter,
	f1.i_pixelflags_saturatedcenter,
	f1.z_pixelflags_saturatedcenter,
	f1.y_pixelflags_saturatedcenter,
	f1.g_pixelflags_interpolatedcenter,
	f1.r_pixelflags_interpolatedcenter,
	f1.i_pixelflags_interpolatedcenter,
	f1.z_pixelflags_interpolatedcenter,
	f1.y_pixelflags_interpolatedcenter,
	f1.g_pixelflags_suspectcenter,
	f1.r_pixelflags_suspectcenter,
	f1.i_pixelflags_suspectcenter,
	f1.z_pixelflags_suspectcenter,
	f1.y_pixelflags_suspectcenter,
	f1.g_pixelflags_clippedcenter,
	f1.r_pixelflags_clippedcenter,
	f1.i_pixelflags_clippedcenter,
	f1.z_pixelflags_clippedcenter,
	f1.y_pixelflags_clippedcenter,
	f1.g_pixelflags_bright_object,
	f1.r_pixelflags_bright_object,
	f1.i_pixelflags_bright_object,
	f1.z_pixelflags_bright_object,
	f1.y_pixelflags_bright_object,
	
	-- 4. Failed SDSS centroid algorithm
	-- It means the center is not well defined.
	f2.g_sdsscentroid_flag,
	f2.r_sdsscentroid_flag,
	f2.i_sdsscentroid_flag,
	f2.z_sdsscentroid_flag,
	f2.y_sdsscentroid_flag,
	
	-- 5. Number of images contributing at the center
	f1.g_inputcount_value as g_input_count,
	f1.r_inputcount_value as r_input_count,
	f1.i_inputcount_value as i_input_count,
	f1.z_inputcount_value as z_input_count,
	f1.y_inputcount_value as y_input_count,

	-- 6. Flag for measurements from the HSC-R2/I2 filter
	f1.merge_measurement_r2,
	f1.merge_measurement_i2,
	c.corr_rmag,
	c.corr_imag

FROM s20a_wide.forced as f1
	LEFT JOIN s20a_wide.forced2 as f2 USING (object_id)
	LEFT JOIN s20a_wide.forced3 as f3 USING (object_id)
	LEFT JOIN s20a_wide.meas2 as m2 USING (object_id)
	LEFT JOIN s20a_wide.mag_corr as c USING (object_id)
	
WHERE
  tractSearch(object_id, 8524)
  AND f1.isprimary
;
