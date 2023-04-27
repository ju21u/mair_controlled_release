 &time_control
 run_days                            = 0,
 run_hours                           = {{.RH}},
 run_minutes                         = 0,
 run_seconds                         = 0,
 start_year                          = {{.YYYY}}, {{.YYYY}}, {{.YYYY}}, {{.YYYY}}, {{.YYYY}}, {{.YYYY}},
 start_month                         = {{.MM}}, {{.MM}},  {{.MM}}, {{.MM}}, {{.MM}}, {{.MM}},
 start_day                           = {{.SD}},  {{.SD}}, {{.SD}}, {{.SD}}, {{.SD}}, {{.SD}},
 start_hour                          = {{.SH}},  {{.SH}},   {{.SH}}, {{.SH}},  {{.SH}}, {{.SH}},
 start_minute                        = 00,   00,   00,	 00,   00,   00,
 start_second                        = 00,   00,   00,	 00,   00,   00,
 end_year                            = {{.YYYY}}, {{.YYYY}}, {{.YYYY}}, {{.YYYY}}, {{.YYYY}}, {{.YYYY}},
 end_month                           = {{.MM}}, {{.MM}},  {{.MM}}, {{.MM}}, {{.MM}}, {{.MM}},
 end_day                             = {{.ED}},  {{.ED}}, {{.ED}}, {{.ED}}, {{.ED}}, {{.ED}},
 end_hour                            = {{.EH}},  {{.EH}}, {{.EH}}, {{.EH}}, {{.EH}}, {{.EH}},
 end_minute                          = 00,   00,   00,	 00,   00,   00,
 end_second                          = 00,   00,   00,	 00,   00,   00,
 auxinput1_inname                    = 'met_em.d<domain>.<date>',
 interval_seconds                    = 3600,
 input_from_file                     = .true.,.true.,.true.,.true.,.true.,.true.,
 history_interval                    = {{.history_interval}}
 frames_per_outfile                  = 100,    100,    100,    100,    100,   1000,
 restart                             = .false.,
 restart_interval                    = 60,
 io_form_history                     = 2
 io_form_restart                     = 2
 io_form_input                       = 2
 io_form_boundary                    = 2
 auxinput5_inname                    = 'wrfchemi_d<domain>_<date>'
 auxinput5_interval_m                = 720, 720, 720, 720,  720, 720,     ! Anthropogenic emissions input
 io_form_auxinput5                   = 2,
 frames_per_auxinput5                = 1,
 force_use_old_data		     = T,
 /

 &domains
 time_step                           = 6,
 time_step_fract_num                 = 0,
 time_step_fract_den                 = 1,
 max_dom                             = 4,
 e_we                                = 121,   121,   121,   121,   121,   121,
 e_sn                                = 121,   121,   121,   121,   121,   121,
 e_vert                              = 51,    51,    51,    51,    51,    51,
 p_top_requested                     = 5000,
 num_metgrid_levels                  = {{.num_metgrid_levels}},
 num_metgrid_soil_levels             = {{.num_metgrid_soil_levels}},
 dx                                  = 3000,
 dy                                  = 3000,
 grid_id                             = 1,     2,     3,     4,     5,     6,
 parent_id                           = 0,     1,     2,     3,     4,     5, 
 i_parent_start                      = 1,     41,    41,    41,    41,    41,
 j_parent_start                      = 1,     41,    41,    41,    41,    41,
 parent_grid_ratio                   = 1,     3,     3,     3,     3,     3, 
 parent_time_step_ratio              = 1,     3,     3,     3,     3,     3,
 feedback                            = 1,
 smooth_option                       = 0
 /

 &physics
 mp_physics                          = 8,     8,     8,    8,     8,     8,
 ra_lw_physics                       = 4,     4,     4,    4,     4,     4,
 ra_sw_physics                       = 4,     4,     4,    4,     4,     4,
 radt                                = 1,     1,     1,    1,     1,     1,
 sf_sfclay_physics                   = 1,     1,     1,    1,     1,     1,
 sf_surface_physics                  = 2,     2,     2,    2,     2,     2,
 bl_pbl_physics                      = 1,     1,     1,    0,     0,     0, 
 bldt                                = 0,     0,     0,    0,     0,     0,
 ysu_topdown_pblmix                  = 1,
 cu_physics                          = 3,     3,     3,    3,     3,     3,
 cugd_avedx                          = 1,
 cudt                                = 0,     0,     0,    0,     0,     0,
 cu_diag                             = 1,      1,     1,    1,     1,    1,
 cu_rad_feedback                     = .true.,  .true., .true.,  .true., .true.,  .false.,
 icloud                              = 1,
 num_land_cat                        = 21,
 sf_urban_physics                    = 0,     0,     0,    0,     0,     0,
 surface_input_source                = 1,
 num_soil_cat                        = 16,
 usemonalb                           = .true.
 rdmaxalb                            = .true.
 rdlai2d                             = .false.
 sst_update                          = 0,
 isfflx                              = 1,
 ifsnow                              = 0,
 icloud                              = 1,
 prec_acc_dt                         = 180.,
 ishallow                            = 0,

 maxiens                             = 1,
 maxens                              = 3,
 maxens2                             = 3,
 maxens3                             = 16,
 ensdim                              = 144, 
/

 &dynamics
 hybrid_opt                          = 2, 
 w_damping                           = 0,
 diff_opt                            = 2,      2,      2,      2,      2,      2,
 km_opt                              = 4,      4,      4,      3,      3,      3,
 diff_6th_opt                        = 0,      0,      0,      0,      0,      0,
 diff_6th_factor                     = 0.12,   0.12,   0.12,   0.12,   0.12,   0.12,
 base_temp                           = 290.
 damp_opt                            = 3,
 zdamp                               = 5000.,  5000.,  5000.,  5000.,  5000.,  5000.,
 dampcoef                            = 0.2,    0.2,    0.2,    0.2,    0.2,    0.2, 
 khdif                               = 0,      0,      0,      0,      0,      0,
 kvdif                               = 0,      0,      0,      0,      0,      0,
 non_hydrostatic                     = .true., .true., .true., .true., .true., .true.,
 moist_adv_opt                       = 2,      2,      2,      2,      2,      2,
 scalar_adv_opt                      = 2,      2,      2,      2,      2,      2,
 /

&chem
 chem_opt                            = 14,       14,      14, 14,  14, 14,
 chem_in_opt                         = 0,        0,       0, 0, 0, 0,

 have_bcs_chem                       = .true., .false., .false., .false., .false., .false.,
 have_bcs_tracer                     = .false., .false., .false., .false., .false., .false.,

 emi_inname                          = "wrfchemi_d<domain>_<date>",
 input_chem_inname                   = "wrfchemi_d<domain>_<date>",

 gaschem_onoff                       = 0,        1,       1, 1, 1, 1,
 aerchem_onoff                       = 0,        1,       1, 1, 1, 1,

 vertmix_onoff                       = 1,        1,       1, 1, 1, 1,
 chem_conv_tr                        = 1,        1,       1, 1, 0, 0,
 gas_drydep_opt                      = 0,        0,       0, 0, 0, 0,
 aer_drydep_opt                      = 0,        0,       0, 0, 0, 0,
 diagnostic_chem                     = 0,        0,       0, 0, 0, 0,

 chemdt                              = 0.033,    0,       0, 0, 0, 0,
 bioemdt                             = 30,       30,      30, 0, 0, 0,

 emiss_inpt_opt                      = 1,        1,       1, 1, 1, 1,
 emiss_opt                           = 12,       12,      12, 12, 12, 12,

 kemit                               = 1,
 io_style_emissions                  = 2,
 aircraft_emiss_opt                  = 0,

 bio_emiss_opt                       = 0,        0,     0, 0, 0, 0,

 phot_opt                            = 0,         0,      0, 0, 0, 0,
 photdt                              = 30,        30,     30, 30, 30, 30,

 wetscav_onoff                       = 0,         0,      0, 0,  0, 0,
 cldchem_onoff                       = 0,         0,      0, 0,  0, 0,

conv_tr_aqchem                      = 0,         0,      0, 0,  0, 0,
 conv_tr_wetscav                     = 0,         0,      0, 0,  0, 0,

 seas_opt                            = 0,
 dust_opt                            = 0,
 dmsemis_opt                         = 0,

 biomass_burn_opt                    = 0,         0,       0, 0,  0, 0,
 plumerisefire_frq                   = 30,        30,      30, 30,  30, 30,

 gas_bc_opt                          = 16,        16,      16, 16,  16, 16,
 gas_ic_opt                          = 16,        16,      16, 16,  16, 16,

 aer_bc_opt                          = 1,         1,       1, 1,   1, 1,
 aer_ic_opt                          = 1,         1,       1, 1,   1, 1,

 aer_ra_feedback                     = 0,         0,       0, 0, 0, 0,
 opt_pars_out                        = 0,
/
 &bdy_control
 spec_bdy_width                      = 5,
 specified                           = .true.,
 /

 &namelist_quilt
 nio_tasks_per_group = 0,
 nio_groups = 1,
 /
