#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>


#define STRING_LENGTH 2000
#define TRUE 1
#define FALSE 0

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#define SW_2_PAR 2.3
#define MJ_TO_J 1E6
#define SEC_2_DAY 86400.0
#define DAY_2_SEC 1.0 / SEC_2_DAY
#define J_TO_UMOL 4.57
#define UMOL_TO_J 1.0 / J_TO_UMOL
#define J_TO_MJ 1E-6
#define hPa_2_kPa 0.1
#define DEG_TO_KELVIN 273.15
#define SEC_2_HFHR 1800.0
#define NHRS 48
#define HLFHR_2_SEC 1.0 / 1800.0


typedef struct  {
    int    row_start;
    int    row_end;
    int    col_start;
    int    col_end;
    int    nrows_in_slice;
    int    ncols_in_slice;
    int    ncols;
    int    nrows;
    int    num_land_pixels;
    char   land_mask_fn[STRING_LENGTH];
    int    root_processor;
    int    rank;
    int    size;
    int    nsize;
    int    remainder;
    float  land_id;
    float  cellsize;
    float  xllcorner;
    float  yllcorner;
    float  xurcorner;
    float  yurcorner;
    int    start_yr;
    int    end_yr;
    int    start_yr_forcing;
    int    end_yr_forcing;
    int    start_yr_rad;
    int    end_yr_rad;
    int    start_yr_all;
    int    end_yr_all;
    char   fdir[STRING_LENGTH];
} control;

typedef struct  {
    int   *tmax_dates;
    int   *tmin_dates;
    int   *rain_dates;
    int   *rad_dates;
    int   *vph09_dates;
    int   *vph15_dates;
    int   tmax_size;
    int   tmin_size;
    int   rain_size;
    int   rad_size;
    int   vph09_size;
    int   vph15_size;
    int   tmax_ndays;
    int   tmin_ndays;
    int   rain_ndays;
    int   rad_ndays;
    int   vph09_ndays;
    int   vph15_ndays;
    float *tmax_slice;
    float *tmin_slice;
    float *rain_slice;
    float *rad_slice;
    float *vph09_slice;
    float *vph15_slice;
 } met;


void   clparser(int, char **, control *);
void   initialise_stuff(control *);
void   mask_ij(control *, float *, int *);
void   read_met_data_slice(control *, met *, int *);
void   get_data(control *, char *, int, float **, int **, int *);
int    distribute_ij(control *, int *, int **);
int    distribute(control *, int *, float *, float **, int, int);
void   build_radiation_clim(control *, int *, float *, float **, float **);
void   write_spinup_file(int, int, control *, met *, float *, float *,
                        float *, float *, float *, float *, float *);
void   write_forcing_file(int, int, control *, met *, float *, float *,
                          float *, float *, float *, float *, float *, float *);
void   estimate_dirunal_par(float, float, int, float, float *);
void   disaggregate_rainfall(float, float *rain);
void   estimate_diurnal_temp(float, float, float, float *);
void   estimate_diurnal_vph(float, float, float, float, float *);
int    rand_int(unsigned int, unsigned int);
float  calc_day_length(int, int, float);
float  calc_vpd(float, float);
int    is_leap_year(int);
float  spitters(int, float, float *);
float  day_angle(int);
float  calculate_solar_declination(int, float);
float  calculate_eqn_of_time(float);
float  calculate_solar_noon(float, float);
float  calculate_hour_angle(float, float);
float  calc_extra_terrestrial_rad(int, float);
float  round_to_value(float, float);
void   calculate_solar_geometry(int, float, float, float *);
void   calc_tam_tpm(float *, float *, float, float, float, float);
