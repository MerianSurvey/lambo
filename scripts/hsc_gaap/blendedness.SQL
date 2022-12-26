SELECT object_id, skymap_id, ra, dec, 
	m_g_blendedness_abs, m_r_blendedness_abs, 
	m_i_blendedness_abs, m_z_blendedness_abs, m_y_blendedness_abs, 
	m_g_blendedness_flag, m_r_blendedness_flag, m_i_blendedness_flag, 
	m_z_blendedness_flag, m_y_blendedness_flag
FROM s20a_wide.summary
WHERE skymap_id = 10000000