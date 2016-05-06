4c4
< **              get_chunk_and_distribute
---
> **              get_chunk_and_distribute_subdaily
7,10c7,10
< **              Extract a chunk (num_rows x total_cols) from the AWAP
< **              meteorological data and create a GDAY spinup and forcing
< **              simulation file for each row and col within the chunk of
< **              Australia. The code then divides the pixels between processors
---
> **              Extract a chunk (Nrows x Ncols) from the AWAP meteorological
> **              data and create GDAY spinup and forcing simulation files
> **              for each row and col within the chunk of Australia. The code
> **              divides each of the row,col pairs between different processors
17c17
< **              This version creates daily files.
---
> **              This version creates 30 minute files.
21c21
< **              mpirun -np 8 get_chunk_and_distribute -rs 300 -re 301 \
---
> **              mpirun -np 8 get_chunk_and_distribute_subdaily -rs 300 -re 301 \
29c29
< ** DATE:        6th May, 2016
---
> ** DATE:        4th April, 2016
33c33
< #include "get_chunk_and_distribute.h"
---
> #include "get_chunk_and_distribute_subdaily.h"
215,222c215,223
<         Each set of npairs is the total met data array divided by the number of
<         cores, so we are working on a chunk here. For each chunk we will loop
<         over the i,j pair unpack each met var and write the i,j met driving file
<         for GDAY.
< 
<         The unpacking logic, is reversed from the packing logic, but works fine
<         because the num_day_offset has been appropriately moved to match the
<         k loop
---
>     ** Each set of npairs is the total met data array divided by the number of
>     ** cores, so we are working on a chunk here. For each chunk we will loop
>     ** over the i,j pair unpack each met var and write the i,j met driving file
>     ** for GDAY.
>     **
>     ** The unpacking logic, is reversed from the packing logic, but works fine
>     ** because the num_day_offset has been appropriately moved to match the
>     **
>     k loop
363a365
> 
680d681
< 
1044,1053c1045,1050
<     int   k=0, kk, yr_to_get, st_idx, en_idx, ndays, year;
<     float co2=0.0, ndep=0.0, wind_sp=0.0, atpress=0.0, wind_am=0.0;
<     float wind_pm=0.0, vpd_avg=0.0, par_day=0.0, sw_am=0.0;
<     float Tmean=0.0, Tsoil=0.0, vpd_am=0.0, vpd_pm=0.0;
<     float sw_pm=0.0, sw=0.0, rainfall=0.0, day_length;
<     float tmin_tomorrow;
<     float Tam, Tpm, SEC_TO_DAY, Tavg, sw_w_m2;
<     float MJ_TO_J = 1.0 / 1.0E-6;
<     float J_TO_UMOL = 4.6;
<     float SW_TO_PAR = 0.48;
---
>     int   k=0, kk, yr_to_get, st_idx, en_idx, ndays, year, hod;
>     float co2=0.0, ndep=0.0, wind=0.0, press=0.0;
>     float vpd=0.0, tsoil=0.0;
>     float sw=0.0, day_length;
>     float vph09_tomorrow, vph15_yesterday;
>     float vph[NHRS], rain[NHRS], tair[NHRS], par[NHRS];
1065a1063,1077
>     long odays = 0;
>     for (k = 0; k < len_shuffled_yrs; k++) {
>         yr_to_get = shuffled_yrs[k];
>         if (is_leap_year(yr_to_get)) {
>             odays += 366;
>         } else {
>             odays += 365;
>         }
>     }*/
>     long odays = 10958;
>     int   ovars = 12;
>     long  ocnt;
>     float odata[ovars * odays * NHRS];
> 
>     /*
1070c1082
<     sprintf(ofname, "met_data/spinup/met_spinup_%d_%d.csv", i, j);
---
>     sprintf(ofname, "met_data/spinup/met_spinup_%d_%d.bin", i, j);
1072c1084,1087
< 
---
>     if (ofp == NULL) {
>         fprintf(stderr, "Error opening file for write\n");
>         exit(EXIT_FAILURE);
>     }
1074a1090,1092
> 
>     /*sprintf(ofname, "met_data/spinup/met_spinup_%d_%d.csv", i, j);
>     ofp = fopen(ofname, "wb");
1078c1096
<     fprintf(ofp, "# Daily met: Row:%d x Col:%d ;  Lat:%f x Lon:%f\n", i, j,
---
>     fprintf(ofp, "# Sub-Daily met: Row:%d x Col:%d ;  Lat:%f x Lon:%f\n", i, j,
1082,1086c1100,1103
<     fprintf(ofp, "#--,--,mj/m2/day,c,mm,c,c,c,kPa,kPa,kPa,ppm,t/ha/year,");
<     fprintf(ofp, "m/s,kPa,umol/m2/d,m/s,m/s,mj/m2/am,mj/m2/pm\n");
<     fprintf(ofp, "#year,doy,sw_rad,tair,rain,tsoil,tam,tpm,vpd_am,vpd_pm,");
<     fprintf(ofp, "vpd_avg,co2,ndep,wind,atmos_press,par,wind_am,wind_pm,");
<     fprintf(ofp, "sw_rad_am,sw_rad_pm\n");
---
>     fprintf(ofp, "#--,--,--,mm/30min,umol/m2/s,degC,degC,kPa,ppm,t/ha/30min,");
>     fprintf(ofp, "m/s,kPa,\n");
>     fprintf(ofp, "#year,doy,hod,rain,par,tair,tsoil,vpd,co2,ndep,wind,press\n");
>     */
1088c1105
<     co2 = 285.0;
---
>     co2 = 350.0;        /* spin up using pre 1990 value (1989 = 351.69 */
1090,1093c1107,1108
<     wind_sp = 3.0; /* Haverd et al. 2012 */
<     atpress = 100.0; /* 1000 mb -> kPa, Haverd et al. 2012 */
<     wind_am = wind_sp;
<     wind_pm = wind_sp;
---
>     wind = 3.0; /* Haverd et al. 2012 */
>     press = 100.0; /* 1000 mb -> kPa, Haverd et al. 2012 */
1094a1110
>     ocnt = 0;
1120,1123d1135
<             if (kk+1 > en_idx)
<                 tmin_tomorrow = tmin_ij[kk];
<             else
<                 tmin_tomorrow = tmin_ij[kk+1];
1125,1126c1137,1141
<             calc_tam_tpm(&Tam, &Tpm, tmin_ij[kk], tmin_tomorrow,
<                          tmax_ij[kk], day_length);
---
>             if (kk+1 > en_idx) {
>                 vph09_tomorrow = vph09_ij[kk];
>             } else {
>                 vph09_tomorrow = vph09_ij[kk+1];
>             }
1128,1130c1143,1147
<             Tavg = (tmin_ij[kk] + tmax_ij[kk]) / 2.0;
<             Tsoil = Tavg;
<             Tmean = Tavg;
---
>             if (kk == st_idx) {
>                 vph15_yesterday = vph15_ij[kk];
>             } else {
>                 vph15_yesterday = vph15_ij[kk-1];
>             }
1132,1136c1149,1150
<             /*
<             1 MJ m-2 d-1 = 1000000 J m-2 d-1 / 86400 s d-1
<                            = 11.574 J m-2 s-1
<                            = 11.574 W m-2
<             */
---
> 
>             /* dissagregate drivers */
1141,1148d1154
<             sw_am = sw / 2.0;
<             sw_pm = sw / 2.0;
<             sw_w_m2 = sw * 11.574;
< 
<             SEC_TO_DAY = 3600. * day_length;
<             /*
<                 Convert radiation from W/m2 -> umol/m2/s (PAR).
<                 2.3 umol/J for conversion from sw -> PAR (Monteith & Unsworth).
1150,1156c1156,1167
<             par_day = sw_w_m2 * 2.3 * SEC_TO_DAY;
<             */
<             par_day = sw * MJ_TO_J * J_TO_UMOL * SW_TO_PAR;
< 
<             vpd_am = calc_vpd(Tam, vph09_ij[kk]);
<             vpd_pm = calc_vpd(Tpm, vph09_ij[kk]);
<             vpd_avg = (vpd_am + vpd_pm) / 2.0;
---
>             estimate_dirunal_par(latitude, longitude, doy_cnt+1, sw, &(par[0]));
>             estimate_diurnal_vph(vph09_ij[kk], vph15_ij[kk], vph09_tomorrow,
>                                  vph15_yesterday, &(vph[0]));
>             disaggregate_rainfall(rain_ij[kk], &(rain[0]));
>             estimate_diurnal_temp(tmin_ij[kk], tmax_ij[kk], day_length,
>                                   &(tair[0]));
> 
>             tsoil = 0.0;
>             for (hod = 0; hod < NHRS; hod++) {
>                 tsoil += tair[hod];
>             }
>             tsoil /= (float)NHRS;
1158c1169
<             rainfall = rain_ij[kk];
---
>             for (hod = 0; hod < NHRS; hod++) {
1160,1164c1171
<             fprintf(ofp,
<             "%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
<              year, doy_cnt+1, sw, Tmean, rainfall, Tsoil, Tam, Tpm, vpd_am,
<              vpd_pm, vpd_avg, co2, ndep, wind_sp, atpress, par_day, wind_am,
<              wind_pm, sw_am, sw_pm);
---
>                 vpd = calc_vpd(tair[hod], vph[hod]);
1166c1173,1193
<              doy_cnt++;
---
>                 /* save everything and do a single big dump at the end */
>                 odata[ocnt] = (float)year;
>                 odata[ocnt+1] = (float)doy_cnt+1;
>                 odata[ocnt+2] = (float)hod;
>                 odata[ocnt+3] = rain[hod];
>                 odata[ocnt+4] = par[hod];
>                 odata[ocnt+5] = tair[hod];
>                 odata[ocnt+6] = tsoil;
>                 odata[ocnt+7] = vpd;
>                 odata[ocnt+8] = co2;
>                 odata[ocnt+9] = ndep;
>                 odata[ocnt+10] = wind;
>                 odata[ocnt+11] = press;
>                 /*
>                 fprintf(ofp, "%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
>                         year, doy_cnt+1, hod, rain[hod], par[hod], tair[hod],
>                         tsoil, vpd, co2, ndep, wind, press);
>                 */
>                 ocnt += ovars;
>             }
>             doy_cnt++;
1169a1197,1203
> 
>     if (fwrite(odata, sizeof(float), ovars * odays * NHRS, ofp) !=\
>                                      ovars * odays * NHRS) {
>         fprintf(stderr, "Error writing spinup file\n");
> 	    exit(EXIT_FAILURE);
>     }
> 
1170a1205
>     free(odata);
1188,1199c1223,1235
<     int k=0, kk, jj, yr_to_get, st_idx, en_idx, ndays, doy_cnt, year,st_idx_rad;
<     float co2=0.0, ndep=0.0, wind_sp=0.0, atpress=0.0, wind_am=0.0;
<     float wind_pm=0.0, vpd_avg=0.0, par_day=0.0, sw_am=0.0;
<     float Tmean=0.0, Tsoil=0.0, vpd_am=0.0, vpd_pm=0.0;
<     float sw_pm=0.0, sw=0.0, rainfall=0.0, day_length;
<     float tmin_tomorrow;
<     float Tam, Tpm, SEC_TO_DAY, Tavg, sw_w_m2;
< 
<     float MJ_TO_J = 1.0 / 1.0E-6;
<     float J_TO_UMOL = 4.6;
<     float SW_TO_PAR = 0.48;
<     sprintf(ofname, "met_data/forcing/met_forcing_preindustco2_%d_%d.csv", i, j);
---
>     int k=0, kk, jj, hod, yr_to_get, st_idx, en_idx, ndays, doy_cnt, year;
>     int st_idx_rad, co2_index;
>     float ndep=0.0, wind=0.0, press=0.0;
>     float tsoil=0.0, vpd=0.0;
>     float sw=0.0, day_length;
>     float vph09_tomorrow, vph15_yesterday;
>     float vph[NHRS], rain[NHRS], tair[NHRS], par[NHRS];
> 
>     /* 1990-2011 */
>     float co2[] = {352.97, 354.37, 355.33, 356.0, 357.68, 359.837, 361.462,
>                    363.155, 365.322, 367.348, 368.865, 370.467, 372.522,
>                    374.76, 376.812, 378.812, 380.827, 382.777, 384.8,
>                    387.001, 389.285, 391.563};
1201c1237,1249
<     ofp = fopen(ofname, "wb");
---
>     /*
>     long odays = 0;
>     for (k = c->start_yr_forcing; k <= c->end_yr_forcing; k++) {
>         if (is_leap_year(k)) {
>             odays += 366;
>         } else {
>             odays += 365;
>         }
>     } */
>     long odays = 8035;
>     int   ovars = 12;
>     long  ocnt;
>     float odata[ovars * odays * NHRS];
1202a1251,1256
>     sprintf(ofname, "met_data/forcing/met_forcing_%d_%d.bin", i, j);
>     ofp = fopen(ofname, "wb");
>     if (ofp == NULL) {
>         fprintf(stderr, "Error opening file for write\n");
>         exit(EXIT_FAILURE);
>     }
1205a1260,1264
>     /*
>     sprintf(ofname, "met_data/forcing/met_forcing_%d_%d.csv", i, j);
>     ofp = fopen(ofname, "wb");
> 
> 
1210c1269
<     fprintf(ofp, "# Daily met: Row:%d x Col:%d ;  Lat:%f x Lon:%f\n", i, j,
---
>     fprintf(ofp, "# Sub-Daily met: Row:%d x Col:%d ;  Lat:%f x Lon:%f\n", i, j,
1214,1219c1273,1276
<     fprintf(ofp, "#--,--,mj/m2/day,c,mm,c,c,c,kPa,kPa,kPa,ppm,t/ha/year,");
<     fprintf(ofp, "m/s,kPa,umol/m2/d,m/s,m/s,mj/m2/am,mj/m2/pm\n");
<     fprintf(ofp, "#year,doy,sw_rad,tair,rain,tsoil,tam,tpm,vpd_am,vpd_pm,");
<     fprintf(ofp, "vpd_avg,co2,ndep,wind,atmos_press,par,wind_am,wind_pm,");
<     fprintf(ofp, "sw_rad_am,sw_rad_pm\n");
< 
---
>     fprintf(ofp, "#--,--,--,mm/30min,umol/m2/s,degC,degC,kPa,ppm,t/ha/30min,");
>     fprintf(ofp, "m/s,kPa,\n");
>     fprintf(ofp, "#year,doy,hod,rain,par,tair,tsoil,vpd,co2,ndep,wind,press\n");
>     */
1221d1277
<     co2 = 285.0;
1223,1227c1279,1280
<     wind_sp = 3.0; /* Haverd et al. 2012 */
<     atpress = 100.0; /* 1000 mb -> kPa, Haverd et al. 2012 */
<     wind_am = wind_sp;
<     wind_pm = wind_sp;
< 
---
>     wind = 3.0; /* Haverd et al. 2012 */
>     press = 100.0; /* 1000 mb -> kPa, Haverd et al. 2012 */
1228a1282,1283
>     co2_index = 0;
>     ocnt = 0;
1236,1237c1291,1292
<             year = m->tmax_dates[date_offset];
<             if (year == yr_to_get) {
---
>             year = (int)m->tmax_dates[date_offset];
>             if (year == (int)yr_to_get) {
1256,1257c1311,1312
<                 year = m->rad_dates[date_offset];
<                 if (year == yr_to_get) {
---
>                 year = (int)m->rad_dates[date_offset];
>                 if (year == (int)yr_to_get) {
1270,1271d1324
< 
< 
1273a1327,1331
>             if (kk+1 > en_idx) {
>                 vph09_tomorrow = vph09_ij[kk];
>             } else {
>                 vph09_tomorrow = vph09_ij[kk+1];
>             }
1275,1295c1333,1337
<             if (kk+1 > en_idx)
<                 tmin_tomorrow = tmin_ij[kk];
<             else
<                 tmin_tomorrow = tmin_ij[kk+1];
< 
<             calc_tam_tpm(&Tam, &Tpm, tmin_ij[kk], tmin_tomorrow,
<                          tmax_ij[kk], day_length);
< 
<             Tavg = (tmin_ij[kk] + tmax_ij[kk]) / 2.0;
<             Tsoil = Tavg;
<             Tmean = Tavg;
< 
<             /*1 MJ m-2 d-1 = 1000000 J m-2 d-1 / 86400 s d-1
<                            = 11.574 J m-2 s-1
<                            = 11.574 W m-2 */
<            if (year < 1990 && ndays == 365)
<                sw = rad_clim_nonleap_ij[doy_cnt];
<            else if (year < 1990 && ndays == 366)
<                sw = rad_clim_leap_ij[doy_cnt];
<            else
<                sw = rad_ij[jj];
---
>             if (kk == st_idx) {
>                 vph15_yesterday = vph15_ij[kk];
>             } else {
>                 vph15_yesterday = vph15_ij[kk-1];
>             }
1298,1305c1340,1341
<             /*
<                 There are a sequence (as much as 12 days, perhaps more) of bad
<                 PAR data in the AWAP data for certain pixels. If we hit one of
<                 these instances we are going to infill based on the climatology.
<                 Because it looks like long sequences are missing it makes no
<                 sense to attempt to fill with days around the bad day I think
<             */
<             if (sw < 0.0 && ndays == 365) {
---
>             /* dissagregate drivers */
>             if (year < 1990 && ndays == 365)
1307c1343
<             } else if (sw < 0.0 && ndays == 366) {
---
>             else if (year < 1990 && ndays == 366)
1309c1345,1346
<             }
---
>             else
>                 sw = rad_ij[jj];
1310a1348,1353
>             estimate_dirunal_par(latitude, longitude, doy_cnt+1, sw, &(par[0]));
>             estimate_diurnal_vph(vph09_ij[kk], vph15_ij[kk], vph09_tomorrow,
>                                  vph15_yesterday, &(vph[0]));
>             disaggregate_rainfall(rain_ij[kk], &(rain[0]));
>             estimate_diurnal_temp(tmin_ij[kk], tmax_ij[kk], day_length,
>                                   &(tair[0]));
1312,1314d1354
<             sw_am = sw / 2.0;
<             sw_pm = sw / 2.0;
<             sw_w_m2 = sw * 11.574;
1316,1319c1356,1360
<             SEC_TO_DAY = 3600. * day_length;
<             /*
<                 Convert radiation from W/m2 -> umol/m2/s (PAR).
<                 2.3 umol/J for conversion from sw -> PAR (Monteith & Unsworth).
---
>             tsoil = 0.0;
>             for (hod = 0; hod < NHRS; hod++) {
>                 tsoil += tair[hod];
>             }
>             tsoil /= (float)NHRS;
1321,1322c1362,1376
<             par_day = sw_w_m2 * 2.3 * SEC_TO_DAY;
<             */
---
>             for (hod = 0; hod < NHRS; hod++) {
> 
>                 /*
>                 ** There are a sequence (as much as 12 days, perhaps more) of
>                 ** bad PAR data in the AWAP data for certain pixels. If we hit
>                 ** one of these instances we are going to infill based on the
>                 ** climatology. Because it looks like long sequences are
>                 ** missing it makes nosense to attempt to fill with days
>                 ** around the bad day I think
>                 */
>                 if (par[hod] < 0.0 && ndays == 365) {
>                     par[hod] = rad_clim_nonleap_ij[doy_cnt] * SW_2_PAR;
>                 } else if (sw < 0.0 && ndays == 366) {
>                     par[hod] = rad_clim_leap_ij[doy_cnt] * SW_2_PAR;
>                 }
1323a1378
>                 vpd = calc_vpd(tair[hod], vph[hod]);
1325c1380,1399
<             par_day = sw * MJ_TO_J * J_TO_UMOL * SW_TO_PAR;
---
>                 /* save everything and do a single big dump at the end */
>                 odata[ocnt] = (float)year;
>                 odata[ocnt+1] = (float)doy_cnt+1;
>                 odata[ocnt+2] = (float)hod;
>                 odata[ocnt+3] = rain[hod];
>                 odata[ocnt+4] = par[hod];
>                 odata[ocnt+5] = tair[hod];
>                 odata[ocnt+6] = tsoil;
>                 odata[ocnt+7] = vpd;
>                 odata[ocnt+8] = co2[co2_index];
>                 odata[ocnt+9] = ndep;
>                 odata[ocnt+10] = wind;
>                 odata[ocnt+11] = press;
>                 /*
>                 fprintf(ofp, "%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
>                         year, doy_cnt+1, hod, rain[hod], par[hod], tair[hod],
>                         tsoil, vpd, co2[co2_index], ndep, wind, press);
>                 */
>                 ocnt += ovars;
>             }
1327,1337d1400
<             vpd_am = calc_vpd(Tam, vph09_ij[kk]);
<             vpd_pm = calc_vpd(Tpm, vph09_ij[kk]);
<             vpd_avg = (vpd_am + vpd_pm) / 2.0;
< 
<             rainfall = rain_ij[kk];
< 
<             fprintf(ofp,
<             "%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
<              year, doy_cnt+1, sw, Tmean, rainfall, Tsoil, Tam, Tpm, vpd_am,
<              vpd_pm, vpd_avg, co2, ndep, wind_sp, atpress, par_day, wind_am,
<              wind_pm, sw_am, sw_pm);
1341a1405,1406
>         co2_index++;
>     }
1342a1408,1411
>     if (fwrite(odata, sizeof(float), ovars * odays * NHRS, ofp) !=\
>                                      ovars * odays * NHRS) {
> 	   fprintf(stderr, "Error writing forcing file\n");
> 	   exit(EXIT_FAILURE);
1388,1425d1456
< void calc_tam_tpm(float *Tam, float *Tpm, float Tmin, float Tmin_tomorrow,
<                   float Tmax, float daylength) {
<     /*
<     The diurnal pattern of air temperature T(t) is calculated from Tmax
<     and Tmin on the assumption of a sinusoidal pattern with T = Tmin at
<     sunrise and T = (Tmin + Tmax) / 2 at sunset.
< 
<     Ross assumed that there was a 3/4 sinusoid in temperature from dawn to
<     dusk, ie started at Tmin, went to Tmax 2/3 of the way through the day,
<     then decayed to Tav at dusk.
< 
<     If Tav = (Tmin + Tmax) / 2 and Tampl = (Tmax - Tmin) / 2 then the
<     time course of daytime temperature is described by:
< 
<     Tav - Tampl * cos(t)
< 
<     where t goes from 0 (dawn) to 3 * pi / 2 (dusk).
< 
<     To get morning and afternoon averages you integrate this formula to get
<     Tam and Tpm
< 
<     Reference
<     ----------
<     * McMurtie et al (1990) Modelling the Yield of Pinus radiata on a Site
<       Limited by Water and Nitrogen. Foremainder Ecology and Management, 30,
<       381-413.
<     */
<     float Tav, Tampl;
< 
<     Tav = (Tmin + Tmax) / 2.0;
<     Tampl = (Tmax - Tmin) / 2.0;
< 
<     *Tam = Tav - Tampl * (1.0 / sqrt(2.0)) / (3.0 * M_PI / 4.0);
<     *Tpm = Tav + Tampl * (1.0 + 1.0 / sqrt(2.0)) / (3.0 * M_PI / 4.0);
< 
<     return;
< }
< 
1441,1442c1472
<     float hPa_2_kPa = 0.1;
<     float DEG_TO_KELVIN = 273.15;
---
> 
1471a1502,2028
> 
> void estimate_diurnal_vph(float vph09, float vph15, float vph09_next,
>                           float vph15_prev, float *vph) {
>     /*
>     Interpolate VPH between 9am and 3pm values to generate diurnal VPD
>     following the method of Haverd et al. This seems reasonable, vapour pressure
>     plotted aginst time of day often does not reveal consistent patterns, with
>     small fluctuations (see Kimball and Bellamy, 1986).
>     Reference:
>     ---------
>     * Haverd et al. (2013) Multiple observation types reduce uncertainty in
>       Australia's terrestrial carbon and water cycles. Biogeosciences, 10,
>       2011-2040.
>     */
>     /* number of hours gap, i.e. 3pm to 9am the next day */
>     float gap = 18.0;
>     float hour;
>     int    i;
> 
>     for (i = 1; i < NHRS+1; i++) {
>         /* first zero values */
>         *(vph+i) = 0.0;
> 
>         hour = (float)i / 2.0;
> 
>         if (hour <= 9.0) {
>            *(vph+(i-1)) = vph15_prev + (vph09 - vph15_prev) * (9.0 + hour) / gap;
>        } else if (hour > 9.0 && hour <= 15.0) {
>            *(vph+(i-1)) = vph09 + (vph15 - vph09) * (hour - 9.0) / (15.0 - 9.0);
>         } else if (hour > 15.0) {
>             *(vph+(i-1)) =  vph15 + (vph09_next - vph15) * (hour - 15.0) / gap;
>         }
>     }
> 
>     return;
> }
> 
> int rand_int(unsigned int min, unsigned int max) {
> 
>     int   value;
>     float scaled;
> 
>     scaled = (float)rand() / RAND_MAX;
>     value = (int)(max - min + 1) * scaled + min;
> 
>     return (value);
> }
> 
> 
> void disaggregate_rainfall(float rain_day, float *rain) {
>     /*
>     Assign daily PPT total to hours of the day, following MAESTRA, which follows
>     algorithm from GRAECO (model of D. Loustau).
>     Reference:
>     * Loustau, D., F. Pluviaud, A. Bosc, A. Porte, P. Berbigier, M. Deque
>       and V. Perarnaud. 2001. Impact of a regional 2 x CO2 climate scenario
>       on the water balance, carbon balance and primary production
>       of maritime pine in southwestern France. In Models for the Sustainable
>       Management of Plantation Forests. Ed. M. Tome. European
>       Cultivated Forest Inst., EFI Proc. No. 41D, Bordeaux, pp 45-58.
>     */
>     int   i, j, hour_index, num_hrs_with_rain;
>     float rate;
> 
>     /* zero everything before we start */
>     for (i = 0; i < NHRS; i++) {
>         *(rain+i) = 0.0;
>     }
> 
>     if (rain_day <= 2.0) {
>         /* All rain falls in one hour for light storms (<2 mm) */
>         hour_index = rand_int(0, 47);
>         *(rain+hour_index) = rain_day;
> 
>     } else if (rain_day > 46.0) {
>         /* All rain falls in 24 hours for storms >46 mm */
>         for (i = 0; i < NHRS; i++) {
>             *(rain+i) = rain_day / (float)NHRS;
>         }
> 
>     } else {
>         /*
>         ** Aim if for all rain to fall at 2mm/hour at a random time of the day.
>         ** If we generate the same random number, then we increase rainfall
>         ** for this hour
>         */
>         num_hrs_with_rain = (int)(rain_day / 2.0);
>         rate = rain_day / (float)num_hrs_with_rain;
> 
>         for (j = 0; j < num_hrs_with_rain; j++) {
>             hour_index = rand_int(0, 47);
>             *(rain+hour_index) += rate;
>         }
>     }
> 
>     return;
> }
> 
> 
> void estimate_diurnal_temp(float tmin, float tmax, float day_length,
>                            float *tair) {
>     /*
>     Calculate diurnal temperature following Parton and Logan
>     the day is divided into two segments and using a truncated sine wave
>     in the daylight and an exponential decrease in temperature
>     at night.
>     TO DO:
>     - Hours between 00:00 and sunrise should be modelled using the previous
>       days information.
>     References:
>     ----------
>     * Parton and Logan (1981) A model for dirunal variation in soil and
>        air temperature. Agricultural Meteorology, 23, 205--216.
>     * Kimball and Bellamy (1986) Energy in Agriculture, 5, 185-197.
> 
>     */
>     /* 1.5 m air temperature values from Parton and Logan, table 1 */
>     float a = 1.86;
>     float b = 2.2;     /* nighttime coeffcient */
>     float c = -0.17;   /* lag of the min temp from the time of runrise */
> 
>     float night_length = 24.0 - day_length;
>     float sunrise = 12.0 - day_length / 2.0 + c;
>     float sunset = 12.0 + day_length / 2.0;
>     float m, n, d, tset, hour;
>     int   i;
> 
>     /* temperature at sunset */
>     m = sunset - sunrise + c;
>     tset = (tmax - tmin) * sin(M_PI * m / (day_length + 2.0 * a)) + tmin;
> 
> 
>     for (i = 1; i < NHRS+1; i++) {
> 
>         hour = (float)i / 2.0;
> 
>         /* hour - time of the minimum temperature (accounting for lag time) */
>         m = hour - sunrise + c;
>         if (hour >= sunrise && hour <= sunset) {
>             *(tair+(i-1)) = tmin + (tmax - tmin) * \
>                         sin((M_PI * m) / (day_length + 2.0 * a));
>         } else {
>             if (hour > sunset) {
>                 n = hour - sunset;
>             } else if (hour < sunrise) {
>                 n = (24.0 + hour) - sunset;
>             }
> 
>             d = (tset - tmin) / (exp(b) - 1.0);
> 
>             /* includes missing displacement to allow T to reach Tmin, this
>             ** removes a discontinuity in the original Parton and Logan eqn.
>             ** See Kimball and Bellamy (1986) Energy in Agriculture, 5, 185-197
>             **/
>             *(tair+(i-1)) = (tmin -d) + (tset - tmin - d) * \
>                         exp(-b * n / (night_length + c));
>         }
>     }
> 
>     return;
> }
> 
> 
> void estimate_dirunal_par(float lat, float lon, int doy, float sw_rad_day,
>                           float *par) {
>     /*
>         Calculate daily course of incident PAR from daily totals using routine
>         from MAESTRA
>     */
>     int   i;
>     float cos_zenith[NHRS];
>     float tau = 0.76;            /* Transmissivity of atmosphere */
>     float direct_frac, diffuse_frac;
>     float cos_bm[NHRS], cos_df[NHRS], sum_bm, sum_df;
>     float zenith, rddf, rdbm, par_day, beam_rad, diffuse_rad;
> 
>     /* MJ m-2 d-1 -> J m-2 s-1 = W m-2 -> umol m-2 s-1 -> MJ m-2 d-1 */
>     par_day = sw_rad_day * MJ_TO_J * DAY_2_SEC * SW_2_PAR * \
>               UMOL_TO_J * J_TO_MJ * SEC_2_DAY;
> 
>     calculate_solar_geometry(doy, lat, lon, &(cos_zenith[0]));
>     diffuse_frac = spitters(doy, par_day, cos_zenith);
>     direct_frac = 1.0 - diffuse_frac;
> 
>     /* daily total beam PAR (MJ m-2 d-1) */
>     beam_rad = par_day * direct_frac;
> 
>     /* daily total diffuse PAR (MJ m-2 d-1) */
>     diffuse_rad = par_day * diffuse_frac;
> 
>     sum_bm = 0.0;
>     sum_df = 0.0;
>     for (i = 0; i < NHRS; i++) {
>         cos_bm[i] = 0.0;
>         cos_df[i] = 0.0;
> 
>         if (cos_zenith[i] > 0.0) {
>             zenith = acos(cos_zenith[i]);
> 
>             /* set FBM = 0.0 for ZEN > 80 degrees */
>             if (zenith < (80.0 * M_PI / 180.0)) {
>                 cos_bm[i] = cos_zenith[i] * pow(tau, (1.0 / cos_zenith[i]));
>             } else {
>                 cos_bm[i] = 0.0;
>             }
>             cos_df[i] = cos_zenith[i];
>             sum_bm += cos_bm[i];
>             sum_df += cos_df[i];
>         }
>     }
> 
>     for (i = 0; i < NHRS; i++) {
> 
>         if (sum_bm > 0.0) {
>             rdbm = beam_rad * cos_bm[i] / sum_bm;
>         } else {
>             rdbm = 0.0;
>         }
> 
>         if (sum_df > 0.0) {
>             rddf = diffuse_rad * cos_df[i] / sum_df;
>         } else {
>             rddf = 0.0;
>         }
> 
>         /* MJ m-2 30min-1 -> J m-2 s-1 -> umol m-2 s-1 */
>         *(par+i) = (rddf + rdbm) * MJ_TO_J * J_TO_UMOL * HLFHR_2_SEC;
>     }
> 
>     return;
> }
> 
> void calculate_solar_geometry(int doy, float latitude, float longitude,
>                               float *cos_zenith) {
>     /*
>     The solar zenith angle is the angle between the zenith and the centre
>     of the sun's disc. The solar elevation angle is the altitude of the
>     sun, the angle between the horizon and the centre of the sun's disc.
>     Since these two angles are complementary, the cosine of either one of
>     them equals the sine of the other, i.e. cos theta = sin beta. I will
>     use cos_zen throughout code for simplicity.
> 
>     Arguments:
>     ----------
>     doy : float
>         day of year
>     latitude : float
>         latitude (degrees)
>     longitude : float
>         longitude (degrees)
> 
>     References:
>     -----------
>     * De Pury & Farquhar (1997) PCE, 20, 537-557.
>     */
>     int   i;
>     float rdec, et, t0, h, gamma, rlat, sin_beta;
>     float hod;
> 
>     for (i = 1; i < NHRS+1; i++) {
> 
>         /* need to convert 30 min data, 0-47 to 0-23.5 */
>         hod = i / 2.0;
> 
>         gamma = day_angle(doy);
>         rdec = calculate_solar_declination(doy, gamma);
>         et = calculate_eqn_of_time(gamma);
>         t0 = calculate_solar_noon(et, longitude);
>         h = calculate_hour_angle(hod, t0);
>         rlat = latitude * M_PI / 180.0;
> 
>         /* A13 - De Pury & Farquhar */
>         sin_beta = sin(rlat) * sin(rdec) + cos(rlat) * cos(rdec) * cos(h);
>         /* The same thing, going to use throughout */
>         *(cos_zenith+(i-1)) = sin_beta;
>         if (*(cos_zenith+(i-1)) > 1.0) {
>             *(cos_zenith+(i-1)) = 1.0;
>         } else if (cos_zenith[i-1] < 0.0) {
>             *(cos_zenith+(i-1)) = 0.0;
>         }
>         /*zenith = 180.0 / M_PI * acos(cos_zenith[i-1]);
>         elevation = 90.0 - zenith;*/
>     }
>     return;
> }
> 
> float day_angle(int doy) {
>     /* Calculation of day angle - De Pury & Farquhar, '97: eqn A18
> 
>     Reference:
>     ----------
>     * De Pury & Farquhar (1997) PCE, 20, 537-557.
>     * J. W. Spencer (1971). Fourier series representation of the position of
>       the sun.
> 
>     Returns:
>     ---------
>     gamma - day angle in radians.
>     */
> 
>     return (2.0 * M_PI * ((float)doy - 1.0) / 365.0);
> }
> 
> float calculate_solar_declination(int doy, float gamma) {
>     /*
>     Solar Declination Angle is a function of day of year and is indepenent
>     of location, varying between 23deg45' to -23deg45'
> 
>     Arguments:
>     ----------
>     doy : int
>         day of year, 1=jan 1
>     gamma : float
>         fractional year (radians)
> 
>     Returns:
>     --------
>     dec: float
>         Solar Declination Angle [radians]
> 
>     Reference:
>     ----------
>     * De Pury & Farquhar (1997) PCE, 20, 537-557.
>     * Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.
>     * J. W. Spencer (1971). Fourier series representation of the position of
>       the sun.
>     */
>     float decl;
> 
>     /* declination (radians) */
>     /*decl = 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) - \
>            0.006758 * cos(2.0 * gamma) + 0.000907 * sin(2.0 * gamma) -\
>            0.002697 * cos(3.0 * gamma) + 0.00148 * sin(3.0 * gamma);*/
> 
> 
>     /* (radians) A14 - De Pury & Farquhar  */
>     decl = -23.4 * (M_PI / 180.) * cos(2.0 * M_PI * ((float)doy + 10.) / 365.);
> 
>     return (decl);
> 
> }
> 
> float calculate_eqn_of_time(float gamma) {
>     /* Equation of time - correction for the difference btw solar time
>     and the clock time.
> 
>     Arguments:
>     ----------
>     doy : int
>         day of year
>     gamma : float
>         fractional year (radians)
> 
>     References:
>     -----------
>     * De Pury & Farquhar (1997) PCE, 20, 537-557.
>     * Campbell, G. S. and Norman, J. M. (1998) Introduction to environmental
>       biophysics. Pg 169.
>     * J. W. Spencer (1971). Fourier series representation of the position of
>       the sun.
>     * Hughes, David W.; Yallop, B. D.; Hohenkerk, C. Y. (1989),
>       "The Equation of Time", Monthly Notices of the Royal Astronomical
>       Society 238: 1529–1535
>     */
>     float et;
> 
>     /* radians */
>     et = 0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) -\
>          0.014615 * cos(2.0 * gamma) - 0.04089 * sin(2.0 * gamma);
> 
>     /* radians to minutes */
>     et *= 229.18;
> 
>     /* radians to hours */
>     /*et *= 24.0 / (2.0 * M_PI);*/
> 
>     /* minutes - de Pury and Farquhar, 1997 - A17 */
>     /*et = (0.017 + 0.4281 * cos(gamma) - 7.351 * sin(gamma) - 3.349 *
>           cos(2.0 * gamma) - 9.731  * sin(gamma));*/
> 
>     return (et);
> }
> 
> float calculate_solar_noon(float et, float longitude) {
>     /* Calculation solar noon - De Pury & Farquhar, '97: eqn A16
> 
>     Reference:
>     ----------
>     * De Pury & Farquhar (1997) PCE, 20, 537-557.
> 
>     Returns:
>     ---------
>     t0 - solar noon (hours).
>     */
>     float t0, Ls;
> 
>     /* all international standard meridians are multiples of 15deg east/west of
>        greenwich */
>     Ls = round_to_value(longitude, 15.);
>     t0 = 12.0 + (4.0 * (Ls - longitude) - et) / 60.0;
> 
>     return (t0);
> }
> 
> float calculate_hour_angle(float t, float t0) {
>     /* Calculation solar noon - De Pury & Farquhar, '97: eqn A15
> 
>     Reference:
>     ----------
>     * De Pury & Farquhar (1997) PCE, 20, 537-557.
> 
>     Returns:
>     ---------
>     h - hour angle (radians).
>     */
>     return (M_PI * (t - t0) / 12.0);
> 
> }
> 
> 
> float spitters(int doy, float par, float *cos_zenith) {
>     /*
>     Spitters algorithm to estimate the diffuse component from the total daily
>     incident radiation.
> 
>     NB. Eqns. 2a-d, not 20a-d
> 
>     Parameters:
>     ----------
>     doy : int
>         day of year
>     par : float
>         daily total photosynthetically active radiation (MJ m-2 d-1)
>     cos_zenith : float
>         cosine of zenith angle (radians)
> 
>     Returns:
>     -------
>     diffuse : float
>         diffuse component of incoming radiation
> 
>     References:
>     ----------
>     * Spitters, C. J. T., Toussaint, H. A. J. M. and Goudriaan, J. (1986)
>       Separating the diffuse and direct component of global radiation and its
>       implications for modeling canopy photosynthesis. Part I. Components of
>       incoming radiation. Agricultural Forest Meteorol., 38:217-229.
>     */
> 
>     /* Fraction of global radiation that is PAR */
>     float fpar = 0.5;
>     float conv = SEC_2_HFHR * J_TO_MJ;
>     float S0, tau, diffuse_frac;
>     int   i;
> 
> 
>     /* Calculate extra-terrestrial radiation */
>     S0 = 0.0;
>     for (i = 1; i < NHRS+1; i++) {
>         S0 += calc_extra_terrestrial_rad(doy, *(cos_zenith+(i-1))) * conv;
>     }
> 
>     /* atmospheric transmisivity */
>     tau = (par / fpar) / S0;
> 
>     /* Spitter's formula (Eqns. 2a-d) */
>     if (tau < 0.07) {
>         diffuse_frac = 1.0;
>     } else if (tau < 0.35) {
>         diffuse_frac = 1.0 - 2.3 * (tau - 0.07) * (tau - 0.07);
>     } else if (tau < 0.75) {
>         diffuse_frac = 1.33 - 1.46 * tau;
>     } else {
>         diffuse_frac = 0.23;
>     }
> 
>     return (diffuse_frac);
> }
> 
> float calc_extra_terrestrial_rad(int doy, float cos_zenith) {
>     /* Solar radiation incident outside the earth's atmosphere, e.g.
>     extra-terrestrial radiation. The value varies a little with the earths
>     orbit.
> 
>     Using formula from Spitters not Leuning!
> 
>     Arguments:
>     ----------
>     doy : double
>         day of year
>     cos_zenith : double
>         cosine of zenith angle (radians)
> 
>     Returns:
>     --------
>     So : float
>         solar radiation normal to the sun's bean outside the Earth's atmosphere
>         (J m-2 s-1)
> 
>     Reference:
>     ----------
>     * Spitters et al. (1986) AFM, 38, 217-229, equation 1.
>     */
> 
>     float So, Sc;
> 
>     /* Solar constant (J m-2 s-1) */
>     Sc = 1370.0;
> 
>     if (cos_zenith > 0.0) {
>         /*
>         ** remember sin_beta = cos_zenith; trig funcs are cofuncs of each other
>         ** sin(x) = cos(90-x) and cos(x) = sin(90-x).
>         */
>         So = Sc * (1.0 + 0.033 * cos((float)doy / 365.0 * 2.0 * M_PI)) *\
>                 cos_zenith;
>     } else {
>         So = 0.0;
>     }
> 
>     return (So);
> 
> }
> 
> float round_to_value(float number, float roundto) {
>     return (round(number / roundto) * roundto);
> }
