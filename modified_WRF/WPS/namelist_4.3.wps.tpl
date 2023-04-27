&share
 wrf_core = 'ARW',
 max_dom = 4,
 start_date = '{{.start_YMD_HMS}}', '{{.start_YMD_HMS}}', '{{.start_YMD_HMS}}', '{{.start_YMD_HMS}}', '{{.start_YMD_HMS}}', '{{.start_YMD_HMS}}',
 end_date   = '{{.end_YMD_HMS}}', '{{.end_YMD_HMS}}', '{{.end_YMD_HMS}}', '{{.end_YMD_HMS}}', '{{.end_YMD_HMS}}', '{{.end_YMD_HMS}}',
 interval_seconds = 3600,
 io_form_geogrid = 2,
/

&geogrid
 parent_id         =   0,   1,  2, 3, 4,
 parent_grid_ratio =   1,   3,  3, 3, 3, 
 i_parent_start    =   1,  41, 41, 41, 41,
 j_parent_start    =   1,  41, 41, 41, 41,
 e_we              =  121, 121, 121, 121,
 e_sn              =  121, 121, 121, 121, 
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT NOTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! The default datasets used to produce the HGT_M, GREENFRAC, 
 ! and LU_INDEX/LANDUSEF fields have changed in WPS v3.8. The HGT_M field
 ! is now interpolated from 30-arc-second USGS GMTED2010, the GREENFRAC 
 ! field is interpolated from MODIS FPAR, and the LU_INDEX/LANDUSEF fields 
 ! are interpolated from 21-class MODIS.
 !
 ! To match the output given by the default namelist.wps in WPS v3.7.1, 
 ! the following setting for geog_data_res may be used:
 !
 ! geog_data_res = 'gtopo_10m+usgs_10m+nesdis_greenfrac+10m','gtopo_2m+usgs_2m+nesdis_greenfrac+2m',
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT NOTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 geog_data_res = 'default','default',
 dx = 3000, !1000, 333.33333333, 111.11111111,
 dy = 3000, !1000, 333.33333333, 111.11111111,
 map_proj = 'lambert',
 ref_lat   = {{.target_lat}},
 ref_lon   = {{.target_lon}},
 truelat1  =  32,
 stand_lon = -102,

 geog_data_path = '{{.GEOG_PATH}}'
/

&ungrib
 out_format = 'WPS',
 prefix = 'FILE_ERA5',
/

&metgrid
 fg_name = 'FILE_ERA5'
 io_form_metgrid = 2, 
/
