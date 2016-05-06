0a1,3
> #ifndef GET_CHUNK_H
> #define GET_CHUNK_H
> 
15d17
< 
64a67,68
>     int    start_yr_all;
>     int    end_yr_all;
104c108
< void write_spinup_file(int, int, control *, met *, float *, float *,
---
> void   write_spinup_file(int, int, control *, met *, float *, float *,
106,110c110,116
< void write_forcing_file(int, int, control *, met *, float *, float *,
<                         float *, float *, float *, float *, float *, float *);
< 
< 
< 
---
> void   write_forcing_file(int, int, control *, met *, float *, float *,
>                           float *, float *, float *, float *, float *, float *);
> void   estimate_dirunal_par(float, float, int, float, float *);
> void   disaggregate_rainfall(float, float *rain);
> void   estimate_diurnal_temp(float, float, float, float *);
> void   estimate_diurnal_vph(float, float, float, float, float *);
> int    rand_int(unsigned int, unsigned int);
112d117
< void   calc_tam_tpm(float *, float *, float, float, float, float);
114a120,131
> float  spitters(int, float, float *);
> float  day_angle(int);
> float  calculate_solar_declination(int, float);
> float  calculate_eqn_of_time(float);
> float  calculate_solar_noon(float, float);
> float  calculate_hour_angle(float, float);
> float  calc_extra_terrestrial_rad(int, float);
> float  round_to_value(float, float);
> void   calculate_solar_geometry(int, float, float, float *);
> 
> 
> #endif
