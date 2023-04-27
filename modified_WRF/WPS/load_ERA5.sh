#! /bin/csh -f
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-02:10          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p test   # Partition to submit to
#SBATCH --mem=36000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o myoutput_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-type=ALL
#SBATCH --mail-user=achulakadabba@g.harvard.edu


#
#hell script to download selected files from rda.ucar.edu using Wget
# NOTE: if you want to run under a different shell, make sure you change
#       the 'set' commands according to your shell's syntax
# after you save the file, don't forget to make it executable
#   i.e. - "chmod 755 <name_of_script>"
#
# Experienced Wget Users: add additional command-line flags here
#   Use the -r (--recursive) option with care
#   Do NOT use the -b (--background) option - simultaneous file downloads
#       can cause your data access to be blocked
set opts = "-N"
#
# Replace xxxxxx with your rda.ucar.edu password on the next uncommented line
# IMPORTANT NOTE:  If your password uses a special character that has special
#                  meaning to csh, you should escape it with a backslash
#                  Example:  set passwd = "my\!password"
set passwd = 'YOUR_PASSWORD'
set num_chars = `echo "$passwd" |awk '{print length($0)}'`
if ($num_chars == 0) then
  echo "You need to set your password before you can continue"
  echo "  see the documentation in the script"
  exit
endif
@ num = 1
set newpass = ""
while ($num <= $num_chars)
  set c = `echo "$passwd" |cut -b{$num}-{$num}`
  if ("$c" == "&") then
    set c = "%26";
  else
    if ("$c" == "?") then
      set c = "%3F"
    else
      if ("$c" == "=") then
        set c = "%3D"
      endif
    endif
  endif
  set newpass = "$newpass$c"
  @ num ++
end
set passwd = "$newpass"
#
set cert_opt = ""
# If you get a certificate verification error (version 1.10 or higher),
# uncomment the following line:
#set cert_opt = "--no-check-certificate"
#
if ("$passwd" == "xxxxxx") then
  echo "You need to set your password before you can continue - see the documentation in the script"
  exit
endif
#
# authenticate - NOTE: You should only execute this command ONE TIME.
# Executing this command for every data file you download may cause
# your download privileges to be suspended.
wget $cert_opt -O auth_status.rda.ucar.edu --save-cookies auth.rda.ucar.edu.$$ --post-data="email={YOUR_EMAIL}&passwd=$passwd&action=login" https://rda.ucar.edu/cgi-bin/login
#
# download the file(s)
# NOTE:  if you get 403 Forbidden errors when downloading the data files, check
#        the contents of the file 'auth_status.rda.ucar.edu'


wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_015_aluvp.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_016_aluvd.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_017_alnip.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_018_alnid.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_031_ci.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_032_asn.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_033_rsn.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_034_sstk.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_035_istl1.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_036_istl2.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_037_istl3.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_038_istl4.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_039_swvl1.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_040_swvl2.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_041_swvl3.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_042_swvl4.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_059_cape.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_066_lailv.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_067_laihv.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_078_tclw.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_079_tciw.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_134_sp.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_136_tcw.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_137_tcwv.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_139_stl1.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_141_sd.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_148_chnk.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_151_msl.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_159_blh.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_164_tcc.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_165_10u.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_166_10v.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_167_2t.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_168_2d.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_170_stl2.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_183_stl3.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_186_lcc.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_187_mcc.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_188_hcc.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_198_src.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_206_tco3.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_229_iews.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_230_inss.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_231_ishf.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_232_ie.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_235_skt.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_236_stl4.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_238_tsn.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_243_fal.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_244_fsr.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.128_245_flsr.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.228_010_lblt.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.228_011_ltlt.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.228_012_lshf.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.228_013_lict.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.228_014_licd.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.228_089_tcrw.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.228_090_tcsw.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.228_131_u10n.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.228_132_v10n.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.228_246_100u.ll025sc.2020030100_2020033123.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.sfc/202003/e5.oper.an.sfc.228_247_100v.ll025sc.2020030100_2020033123.grb

wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.pl/202003/e5.oper.an.pl.128_060_pv.ll025sc.2020030200_2020030223.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.pl/202003/e5.oper.an.pl.128_075_crwc.ll025sc.2020030200_2020030223.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.pl/202003/e5.oper.an.pl.128_076_cswc.ll025sc.2020030200_2020030223.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.pl/202003/e5.oper.an.pl.128_129_z.ll025sc.2020030200_2020030223.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.pl/202003/e5.oper.an.pl.128_130_t.ll025sc.2020030200_2020030223.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.pl/202003/e5.oper.an.pl.128_131_u.ll025uv.2020030200_2020030223.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.pl/202003/e5.oper.an.pl.128_132_v.ll025uv.2020030200_2020030223.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.pl/202003/e5.oper.an.pl.128_133_q.ll025sc.2020030200_2020030223.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.pl/202003/e5.oper.an.pl.128_135_w.ll025sc.2020030200_2020030223.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.pl/202003/e5.oper.an.pl.128_138_vo.ll025sc.2020030200_2020030223.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.pl/202003/e5.oper.an.pl.128_155_d.ll025sc.2020030200_2020030223.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.pl/202003/e5.oper.an.pl.128_157_r.ll025sc.2020030200_2020030223.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.pl/202003/e5.oper.an.pl.128_203_o3.ll025sc.2020030200_2020030223.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.pl/202003/e5.oper.an.pl.128_246_clwc.ll025sc.2020030200_2020030223.grb
wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds633.0/e5.oper.an.pl/202003/e5.oper.an.pl.128_247_ciwc.ll025sc.2020030200_2020030223.grb
#
#
#
#
#
# clean up
rm auth.rda.ucar.edu.$$ auth_status.rda.ucar.edu
