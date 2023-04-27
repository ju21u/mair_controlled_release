 &time_control
 run_days                            = 0,
 run_hours                           = {{.RH}},
 run_minutes                         = 0,
 run_seconds                         = 0,
 start_year                          = {{.YYYY}}, {{.YYYY}}, {{.YYYY}}, {{.YYYY}}, {{.YYYY}}, {{.YYYY}},
 start_month                         = {{.MM}}, {{.MM}},  {{.MM}}, {{.MM}}, {{.MM}}, {{.MM}},
 start_day                           = {{.SD}},  {{.SD}}, {{.SD}}, {{.SD}}, {{.SD}}, {{.SD}},
 start_hour                          = {{.SH}},  {{.SH}},   {{.SH}}, {{.SH}},  {{.SH}}, {{.SH}},
 start_minute                        = 00,   00,   00, 00,  00, 00,
 start_second                        = 00,   00,   00, 00,  00, 00,
 end_year                            = {{.YYYY}}, {{.YYYY}}, {{.YYYY}}, {{.YYYY}}, {{.YYYY}}, {{.YYYY}},
 end_month                           = {{.MM}}, {{.MM}},  {{.MM}}, {{.MM}}, {{.MM}}, {{.MM}},
 end_day                             = {{.ED}},  {{.ED}}, {{.ED}}, {{.ED}}, {{.ED}}, {{.ED}},
 end_hour                            = {{.EH}},  {{.EH}}, {{.EH}}, {{.EH}}, {{.EH}}, {{.EH}},
 end_minute                          = 00,   00,   00, 00, 00, 00,
 end_second                          = 00,   00,   00, 00, 00, 00,
 input_from_file                     = .true.,.true.,.true., .true.,.true., .true.,
 fine_input_stream                   =  0,   0,     0, 0,   0, 0,
 
 auxinput1_inname                    = 'met_em.d<domain>.<date>',
 interval_seconds                    = 21600,
 history_interval                    = 60,  60,   60,  60,   60, 1,  
 frames_per_outfile                  = 1000, 1000, 1000,  1000, 1000, 1000,
 restart                             = .false.,
 restart_interval                    = 60,
 io_form_history                     = 2
 io_form_restart                     = 2
 io_form_input                       = 2
 io_form_boundary                    = 2
 debug_level                         = 100
 auxinput5_inname                    = 'wrfchemi_d<domain>_<date>'
 auxinput5_interval_m                = 720, 720, 720, 720,  720, 720,     ! Anthropogenic emissions input
 io_form_auxinput5                   = 2,
 frames_per_auxinput5                = 1,
 /

 &domains
 time_step                           = 180,! 6, 1, 0, 
 time_step_fract_num                 = 0,! 0, 1, 1,
 time_step_fract_den                 = 1,! 1, 2, 2,
 max_dom                             = 6,
 e_we                                = 103, 103, 103, 103, 103, 103,
 e_sn                                = 103, 103, 103, 103, 103, 103,
 e_vert                              = 51,    51,    51, 51, 51, 51,
 p_top_requested                     = 5000,
 num_metgrid_levels                  = 34,
 num_metgrid_soil_levels             = 4,
 dx                                  = 27000, 9000, 3000, 1000,  333.3333, 111.1111,
 dy                                  = 27000, 9000, 3000, 1000,  333.3333, 111.1111,
 grid_id                             = 1,     2,     3, 4, 5, 6,
 parent_id                           = 0,     1,     2, 3, 4, 5, 
 i_parent_start                      = 1,    35,    35, 35, 35, 35,
 j_parent_start                      = 1,    35,    35, 35, 35, 35,
 parent_grid_ratio                   = 1,     3,     3, 3,  3, 3,
 parent_time_step_ratio              = 1,     3,     3, 3,  3, 3,
 feedback                            = 0,
 smooth_option                       = 0,
 eta_levels                          =  1.0000,
                                        0.9975,0.9937,0.9884,0.9821,0.9758,0.9695,0.9632,
                                        0.9568,0.9505,0.9442,0.9368,0.9295,0.9221,0.9147,
                                        0.9074,0.9000,0.8926,0.8853,0.8779,0.8705,0.8632,
                                        0.8558,0.8474,0.8389,0.8263,0.8084,0.7895,0.7684,
                                        0.7474,0.7211,0.6947,0.6579,0.6211,0.5842,0.5474,
                                        0.5105,0.4737,0.4337,0.3895,0.3463,0.3053,0.2674,
                                        0.2326,0.2000,0.1684,0.1389,0.1105,0.0758,0.0368
                                        0.0000, 
 sfcp_to_sfcp 			     = T,
 num_traj			     = 0,
 /

 &physics
 mp_physics                          = 8,      4,     4,     4,     4, 4,
 progn                               = 0,     
 
 ra_lw_physics                       = 4,      4,     4,     4,     4, 4,
 ra_sw_physics                       = 4,      4,     4,     4,     4, 4,
 radt                                = 5,   12,    12, 12, 12, 12,
 
 swrad_scat                          = 1
 shadlen                             = 25000.
 slope_rad                           = 0
 topo_shading                        = 0
 
 sf_sfclay_physics                   = 5,     5,     5, 5,  5, 5,
 bl_pbl_physics                      = 5,     5,      5,     5,   0, 0,
 bldt                                = 0,     0,      0,     0,     0, 0,
 iz0tlnd                             = 1,
 
 sf_urban_physics                    = 0,      0,     0,     0,      0, 0,
 pxlsm_smois_init                    = 0,      0,     0,     0,      0, 0,
 fractional_seaice                   = 0
 
 topo_wind                           = 0,     0,     0,      0,      0, 0,
 grav_settling                       = 0
 
 sf_surface_physics                  = 2,     2,     2, 2,   2, 2,
 num_soil_layers                     = 4,
 
 surface_input_source                = 1,
 num_land_cat                        = 21,
 num_soil_cat                        = 16,      
 usemonalb                           = .true.
 rdmaxalb                            = .true.
 rdlai2d                             = .false.
 sst_update                          = 0,
 
 isfflx                              = 1,
 ifsnow                              = 0,
 icloud                              = 1,
 
 cu_physics                          = 3,      3,     3,     3,    3,     3,
 cugd_avedx                          = 1,
 cu_diag                             = 1,      1,     1,    1,     1,    1,
 
 prec_acc_dt                         = 180.,
 cudt                                = 0,      0,     0,    0,      0,     0,
 cu_rad_feedback                     = .true.,  .true., .true.,  .true., .true.,  .false.,
 ishallow                            = 0,
 
 maxiens                             = 1,
 maxens                              = 3,
 maxens2                             = 3,
 maxens3                             = 16,
 ensdim                              = 144,
 ! traj_opt			      = 0, 
 /

 &fdda
 /

 &dynamics
 rk_ord                              = 3
 c_s                                 = 0.25
 
 w_damping                           = 1
 diff_opt                            = 2
 
 km_opt                              = 4, 4,  4, 4, 2, 2,
 
 damp_opt                            = 3,
 zdamp                               = 5000.,  5000., 5000.,  5000.,  5000., 5000.,
 dampcoef                            = 0.2,    0.2,   0.2,    0.2,   0.2, 0.2,
 gwd_opt                             = 0,
 
 khdif                               = 0,      0,      0, 0,    0, 0,
 kvdif                               = 0,      0,      0, 0,    0, 0,
 smdiv                               = 0.1
 emdiv                               = 0.01
 epssm                               = 0.1
 
 do_avgflx_em                        = 1,      1,       1, 1,   1, 1,
 do_avgflx_cugd                      = 1,
 
 mix_full_fields                     = .false.,
 
 mix_isotropic                       = 0,      0,       0, 0,   0, 0,
 mix_upper_bound                     = 0.1,    0.1,     0.1, 0.1, 0.1, 0.1,
 
 diff_6th_opt                        = 2,      2,       2, 2,    2, 2,
 diff_6th_factor                     = 0.12,     0.12,     0.12, 0.12,  0.12, 0.12,
 
 base_temp                           = 290.
 non_hydrostatic                     = .true., .true., .true., .true.,.true., .true.,
 
 moist_adv_opt                       = 2,      2,       2, 2, 2, 2,
 chem_adv_opt                        = 2,      2,       2, 2, 2, 2,
 scalar_adv_opt                      = 2,      2,       2, 2, 2, 2,
 tke_adv_opt                         = 2,      2,       2, 2, 2, 2,
 
 time_step_sound                     = 4,      4,       4, 4, 4, 4,
 h_mom_adv_order                     = 5,      5,       5, 5, 5, 5,
 v_mom_adv_order                     = 3,      3,       3, 3, 3, 3,
 h_sca_adv_order                     = 5,      5,       5, 5, 5, 5,
 v_sca_adv_order                     = 3,      3,       3, 3, 3, 3,
 
 tke_drag_coefficient                = 0.,     0.,      0., 0.,  0., 0.,
 tke_heat_flux                       = 0.,     0.,      0., 0.,  0., 0.,
 
 top_lid                             = .false., .false., .false., .false.,  .false., .false.,
 &bdy_control
 spec_bdy_width                      = 5,
 spec_zone                           = 1,
 relax_zone                          = 4,
 specified                           = .true., .false.,.false., .false.,  .false., .false.,
 nested                              = .false., .true., .true., .true., .true., .true.,

 /
 &grib2
 /

 &namelist_quilt
 nio_tasks_per_group = 0,
 nio_groups = 1,
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
