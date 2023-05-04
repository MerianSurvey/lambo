SELECT object_id, skymap_id, ra, dec, 
	m_g_blendedness_abs, m_r_blendedness_abs, 
	m_i_blendedness_abs, m_z_blendedness_abs, m_y_blendedness_abs, 
	m_g_blendedness_flag, m_r_blendedness_flag, m_i_blendedness_flag, 
	m_z_blendedness_flag, m_y_blendedness_flag, g_apertureflux_10_flux, g_apertureflux_10_flag,
	r_apertureflux_10_flux, r_apertureflux_10_flag, i_apertureflux_10_flux, i_apertureflux_10_flag,
	z_apertureflux_10_flux, z_apertureflux_10_flag, y_apertureflux_10_flux, y_apertureflux_10_flag,
    r_extendedness_value,i_extendedness_value
FROM s20a_wide.summary
WHERE skymap_id = 10000000
