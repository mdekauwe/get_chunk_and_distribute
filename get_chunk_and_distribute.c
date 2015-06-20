/*
** 
** PROGRAM:
**              get_chunk_and_distribute
**
** DESCRIPTION:
**              Extract a chunk (num_rows x total_cols) from the AWAP 
**              meteorological data and create a GDAY spinup and forcing 
**              simulation file for each row and col within the chunk of 
**              Australia. The code then divides the pixels between processors
**              using the MPI libraries.
**
**              We don't have overlapping radiation data, i.e. nothing before 
**              1990, so we need to build a climatology and in fill with that 
**              for our spinup files.
**
** NOTES:
**              run cmd:
**              mpirun -np 8 get_chunk_and_distribute -s 300 -e 301
**                     
**
** AUTHOR:      Martin De Kauwe with MPI assistance from Andrey Bliznyuk @NCI
**
** EMAIL:       mdekauwe@gmail.com
**
** DATE:        15th June, 2015
**
*/

#include "get_chunk_and_distribute.h"

int main(int argc, char **argv)
{
    int    i, j, k, doy, idx;
    long   offset = 0, date_offset, offset_1, offset_2;
    float *land_mask = NULL;
    float *tmax_ij = NULL;
    float *tmin_ij = NULL;
    float *rain_ij = NULL;
    float *vph09_ij = NULL;
    float *vph15_ij = NULL;
    float *rad_ij = NULL;
    float *rad_clim_leap_ij = NULL;
    float *rad_clim_nonleap_ij = NULL;
    int   *years_ij_rad = NULL;
    int   *land_ij = NULL;
    int   *pairs = NULL;
    int    mpi_err =0;
    int    npairs;
    long   pixel_count;
    
    /*MPI_Status status;*/
    FILE *land_mask_fp = NULL;

    control *c;
    c = (control *)malloc(sizeof (control));
	if (c == NULL) {
		fprintf(stderr, "control structure: Not allocated enough memory!\n");
		exit(1);
	}

	met *m;
	m = (met *)malloc(sizeof (met));
	if (m == NULL) {
		fprintf(stderr, "met structure: Not allocated enough memory!\n");
		exit(1);
	}

    initialise_stuff(c);
    clparser(argc, argv, c);


    if (c->debug) {
        c->start_yr = 1990;
        c->end_yr = 1990;
        c->start_yr_rad = 1990;
        c->end_yr_rad = 1990;
        strcpy(c->fdir, "test_files");
    }


    if ((c->row_start == -999) || (c->row_end == -999)) {
        fprintf(stderr,"You need to set row start/end on command line\n");
		exit(EXIT_FAILURE);
    }
    c->nrows_in_slice = c->row_end - c->row_start;
    c->root_processor = 0;
    mpi_err = MPI_Init(&argc, &argv);
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &(c->rank));
    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &(c->size));


    /* do some work as process 0 */
    if (c->rank == c->root_processor) {

        /* read land mask file */
        if ((land_mask_fp = fopen(c->land_mask_fn, "r")) == NULL) {
		    fprintf(stderr, "%s: couldn't open land mask file %s for read\n",
		            argv[0], c->land_mask_fn);
		    MPI_Abort(MPI_COMM_WORLD, -1);
	    }

	    if ((land_mask = (float *)calloc(c->nrows * c->ncols,
	                      sizeof(float))) == NULL ) {
		    fprintf(stderr,"%s: error allocating space for land mask array\n",
		            argv[0]);
		    MPI_Abort(MPI_COMM_WORLD, -1);
	    }

        if (fread(land_mask, sizeof(float), c->nrows * c->ncols,
                 land_mask_fp) != c->nrows * c->ncols) {
            fprintf(stderr,"%s: error in reading land mask file\n", *argv);
		    MPI_Abort(MPI_COMM_WORLD, -1);
	    }
	    fclose(land_mask_fp);

	    c->num_land_pixels = 0;
	    for (i = 0; i < c->nrows_in_slice; i++) {
	        for (j = 0; j < c->ncols; j++) {
	            offset = (i + c->row_start) * c->ncols + j;
	            if (land_mask[offset] > c->land_id) {
	                c->num_land_pixels++;
                }
            }
        }

        if ((land_ij = (int *)calloc(c->num_land_pixels * 2,
                       sizeof(int))) == NULL ) {
		    fprintf(stderr,"%s: error allocating space for land_ij array\n",
		            argv[0]);
		    MPI_Abort(MPI_COMM_WORLD, -1);
	    }
        mask_ij(c, land_mask, land_ij);
    }

    if (MPI_Bcast(&(c->num_land_pixels), 1, MPI_INT, c->root_processor,
                  MPI_COMM_WORLD) != MPI_SUCCESS) {
        fprintf(stderr, "Error broadcasting num_land_pixels\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }

    c->nsize = c->num_land_pixels / c->size;
    c->rest = c->num_land_pixels - (c->nsize * c->size);
    read_met_data_slice(c, m, land_ij);


    npairs = distribute_ij(c, land_ij, &pairs);
    if (npairs < 0) {
        fprintf(stderr,"Error in distribute_ij\n");
		MPI_Abort(MPI_COMM_WORLD, mpi_err);
    }

    /* allocate space for each ij pair */
    if ((tmax_ij = (float *)calloc(m->tmax_ndays, sizeof(float))) == NULL) {
        fprintf(stderr,"Error allocating space for tmax_ij array\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if ((tmin_ij = (float *)calloc(m->tmin_ndays, sizeof(float))) == NULL) {
        fprintf(stderr,"Error allocating space for tmin_ij array\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if ((rain_ij = (float *)calloc(m->rain_ndays, sizeof(float))) == NULL) {
        fprintf(stderr,"Error allocating space for rain_ij array\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if ((vph09_ij = (float *)calloc(m->vph09_ndays, sizeof(float))) == NULL) {
        fprintf(stderr,"Error allocating space for vph09_ij array\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if ((vph15_ij = (float *)calloc(m->vph15_ndays, sizeof(float))) == NULL) {
        fprintf(stderr,"Error allocating space for vph15_ij array\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if ((rad_ij = (float *)calloc(m->rad_ndays, sizeof(float))) == NULL) {
        fprintf(stderr,"Error allocating space for rad_ij array\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if ((years_ij_rad = (int *)calloc(m->rad_ndays,
                        sizeof(int))) == NULL) {
        fprintf(stderr,"Error allocating space for years_ij_rad array\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if ((rad_clim_nonleap_ij = (float *)calloc(365,
                                sizeof(float))) == NULL) {
        fprintf(stderr,"Error allocating space for rad_clim_nonleap_ij array\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if ((rad_clim_leap_ij = (float *)calloc(366,
                             sizeof(float))) == NULL) {
        fprintf(stderr,"Error allocating space for rad_clim_leap_ij array\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    /* 
        Each set of npairs is the total met data array / number of cores,
        so we are working on a chunk here. For each chunk we will loop over
        the i,j pair unpack each met var and write the i,j met driving file
        for GDAY.
    */
    pixel_count = 0;
    for (k = 0; k < npairs*2; k+=2) {

        i = pairs[k];
        j = pairs[k+1];
        for (doy = 0; doy < m->tmax_ndays; doy++) {
            if ((i == 299) && (j == 321)) {
                
                offset = doy * npairs + pixel_count;
                
                printf("%f\n", m->tmax_slice[offset]);
            }
        
        }
        pixel_count++;
    }
    
    
    /* Stop this process */
    mpi_err = MPI_Finalize();

    
    free(c);
    free(m);
    free(land_mask);
    free(land_ij);
    free(years_ij_rad);
    free(tmax_ij);
    free(tmin_ij);
    free(rain_ij);
    free(vph09_ij);
    free(vph15_ij);
    free(rad_ij);
    free(rad_clim_nonleap_ij);
    free(rad_clim_leap_ij);
    
    return (EXIT_SUCCESS);
}

void clparser(int argc, char **argv, control *c) {
	int i;

	for (i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			if (!strncasecmp(argv[i], "-lm", 3)) {
			    strcpy(c->land_mask_fn, argv[++i]);
			} else if (!strncasecmp(argv[i], "-s", 2)) {
			    c->row_start = atoi(argv[++i]);
			    c->row_start--; /* correct for 0 index */
			} else if (!strncasecmp(argv[i], "-e", 2)) {
			    c->row_end = atoi(argv[++i]);
			    c->row_end--; /* correct for 0 index */
			} else if (!strncasecmp(argv[i], "-d", 2)) {
			    c->debug = TRUE;
			} else {
                fprintf(stderr,"%s: unknown argument on command line: %s\n",
                        argv[0], argv[i]);
				exit(EXIT_FAILURE);
     		}
        }
    }
	return;
}

void initialise_stuff(control *c) {

    strcpy(c->land_mask_fn, "DATA/land_mask/AWAP_land_mask.flt");

    c->nrows = 681;
    c->ncols = 841;
	c->row_end = -999;
	c->row_start = -999;
    c->land_id = 0.5; /* Anything > 0.000001 */
    c->nsize = -999;
    c->rest = -999;
    c->debug = FALSE;
    c->cellsize = 0.05;
    c->xllcorner = 111.975;
    c->yllcorner = -44.025;
    
    /*
    c->start_yr = 1950;
    c->end_yr = 2011;
    c->start_yr_rad = 1990;
    c->end_yr_rad = 2011;
    */
    
    c->start_yr = 1950;
    c->end_yr = 1952;
    c->start_yr_rad = 1990;
    c->end_yr_rad = 2011;
    
    strcpy(c->fdir, "AWAP_data");


    return;
}

void mask_ij(control *c, float *land_mask, int *land_ij) {

    int  i, j;
    long out_offset = 0, in_offset = 0;

    for (i = 0; i < c->nrows_in_slice; i++) {
	    for (j = 0; j < c->ncols; j++) {
            in_offset = (i + c->row_start) * c->ncols + j;
	        if (land_mask[in_offset] > c->land_id) {
	            land_ij[out_offset] = c->row_start + i;
	            land_ij[out_offset+1] = j ;
                out_offset += 2;
            }
        }
    }
}

void read_met_data_slice(control *c, met *m, int *land_ij) {

    /*int    i, j, k;
    uint64_t   offset;*/
    float *met_data = NULL;
    int    total_days, total_days_rad;
    int    yr;
    int tag1 = 101, tag2 = 102, tag3 = 103, tag4 = 104, tag5 = 105, tag6 = 106;

    /*
        Count the number of days to size arrays
    */
    total_days = 0;
    for (yr = c->start_yr; yr <= c->end_yr; yr++) {
        if (is_leap_year(yr)) {
            total_days += 366;
        } else {
            total_days += 365;
        }
    }
    
    /*
        Count the number of days to size arrays
        For the Rad data the timeseries is shorter, 1990-2011, as opposed to 
        1950-2011 so we need to build a shorter array.
    */
    total_days_rad = 0;
    for (yr = c->start_yr_rad; yr <= c->end_yr_rad; yr++) {
        if (is_leap_year(yr)) {
            total_days_rad += 366;
        } else {
            total_days_rad += 365;
        }
    }

    /*
    ** Tmax
    */
    m->tmax_ndays = total_days;
    if ((m->tmax_dates = (int *)calloc(total_days * 3, sizeof(int))) == NULL) {
        fprintf(stderr, "Error allocating space for tmax dates array\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if (c->rank == c->root_processor) {
        get_data(c, "tmax", m->tmax_ndays, &met_data, &m->tmax_dates, land_ij);
    }

    if (MPI_Bcast(&(*m->tmax_dates), m->tmax_ndays*3, MPI_INT, c->root_processor,
                  MPI_COMM_WORLD) != MPI_SUCCESS) {
        fprintf(stderr, "Error broadcasting tmax_dates\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }
    m->tmax_size = distribute(c, land_ij, met_data, &m->tmax_slice, tag1,
                              m->tmax_ndays);
    free(met_data);
    met_data = NULL;


    /*
    ** Tmin
    */
    m->tmin_ndays = total_days;
    if ((m->tmin_dates = (int *)calloc(total_days * 3, sizeof(int))) == NULL) {
        fprintf(stderr, "Error allocating space for tmin dates array\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if (c->rank == c->root_processor) {
        get_data(c, "tmin", m->tmin_ndays, &met_data, &m->tmin_dates, land_ij);
    
    }

    if (MPI_Bcast(&(*m->tmin_dates), m->tmin_ndays*3, MPI_INT, c->root_processor,
                  MPI_COMM_WORLD) != MPI_SUCCESS) {
        fprintf(stderr, "Error broadcasting tmin_dates\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }
    m->tmin_size = distribute(c, land_ij, met_data, &m->tmin_slice, tag2,
                              m->tmin_ndays);
    free(met_data);
    met_data = NULL;


    /*
    ** Rain
    */
    m->rain_ndays = total_days;
    if ((m->rain_dates = (int *)calloc(total_days * 3, sizeof(int))) == NULL) {
        fprintf(stderr, "Error allocating space for rain dates array\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if (c->rank == c->root_processor) {
        get_data(c, "rain", m->rain_ndays, &met_data, &m->rain_dates, land_ij);
    }

    if (MPI_Bcast(&(*m->rain_dates), m->rain_ndays*3, MPI_INT, c->root_processor,
                  MPI_COMM_WORLD) != MPI_SUCCESS) {
        fprintf(stderr, "Error broadcasting rain_dates\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }
    m->rain_size = distribute(c, land_ij, met_data, &m->rain_slice, tag3,
                              m->rain_ndays);
    free(met_data);
    met_data = NULL;


    /*
    ** Vph09
    */
    m->vph09_ndays = total_days;
    if ((m->vph09_dates = (int *)calloc(total_days * 3, sizeof(int))) == NULL) {
        fprintf(stderr, "Error allocating space for vph09 dates array\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if (c->rank == c->root_processor) {
        get_data(c, "vph09", m->vph09_ndays, &met_data, &m->vph09_dates, land_ij);
    }

    if (MPI_Bcast(&(*m->vph09_dates), m->vph09_ndays*3, MPI_INT, c->root_processor,
                  MPI_COMM_WORLD) != MPI_SUCCESS) {
        fprintf(stderr, "Error broadcasting vph09_dates\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }
    m->vph09_size = distribute(c, land_ij, met_data, &m->vph09_slice, tag4,
                               m->vph09_ndays);
    free(met_data);
    met_data = NULL;


    /*
    ** Vph15
    */
    m->vph15_ndays = total_days;
    if ((m->vph15_dates = (int *)calloc(total_days * 3, sizeof(int))) == NULL) {
        fprintf(stderr, "Error allocating space for vph15 dates array\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if (c->rank == c->root_processor) {
        get_data(c, "vph15", m->vph15_ndays, &met_data, &m->vph15_dates, land_ij);
    }

    if (MPI_Bcast(&(*m->vph15_dates), m->vph15_ndays*3, MPI_INT, c->root_processor,
                  MPI_COMM_WORLD) != MPI_SUCCESS) {
        fprintf(stderr, "Error broadcasting vph15_dates\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }
    m->vph15_size = distribute(c, land_ij, met_data, &m->vph15_slice, tag5,
                              m->vph15_ndays);
    free(met_data);
    met_data = NULL;



    /*
    ** Rad - timeseries is shorter, 1990-2011, as opposed to 1950-2011
    **       so we need to build a shorter array. We will use a climatology
    **       to infill data for spinup purposes.
    */

    m->rad_ndays = total_days_rad;
    if ((m->rad_dates = (int *)calloc(total_days_rad * 3,
                         sizeof(int))) == NULL) {
        fprintf(stderr, "Error allocating space for rad dates array\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if (c->rank == c->root_processor) {
        get_data(c, "rad", m->rad_ndays, &met_data, &m->rad_dates, land_ij);
    }
    
    if (MPI_Bcast(&(*m->rad_dates), m->rad_ndays*3, MPI_INT, c->root_processor,
                  MPI_COMM_WORLD) != MPI_SUCCESS) {
        fprintf(stderr, "Error broadcasting rad_dates\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
    }
    m->rad_size = distribute(c, land_ij, met_data, &m->rad_slice, tag6,
                              m->rad_ndays);
    
    free(met_data);
    met_data = NULL;
    
    return;
}


void get_data(control *c, char *met_var, int total_days, float **met_data,
              int **dates, int *land_ij) {

    int    day_cnt, i, j, k, dt_cnt;
    long   in_offset, out_offset;
    float *met_data_day = NULL;

    FILE *fp = NULL;
    int   start_yr, end_yr;
    int   days_in_month[] = {0,31,28,31,30,31,30,31,31,30,31,30,31};
    int   yr, mth, day, ndays;
    char  imth[3];
    char  iday[3];
    long pixel_count;
    char  infname[STRING_LENGTH];
    
    
    if ((*met_data = (float *)calloc(total_days * c->num_land_pixels,
                      sizeof(float))) == NULL ) {
        fprintf(stderr, "Error allocating space for met_data array\n");
    }
    
    if ((met_data_day = (float *)calloc(c->nrows * c->ncols,
                         sizeof(float))) == NULL) {
         fprintf(stderr, "Error allocating space for met array\n");
         MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if ((strncmp(met_var, "rad", 3) == 0)) {
        start_yr = c->start_yr_rad;
        end_yr = c->end_yr_rad;
    } else {
        start_yr = c->start_yr;
        end_yr = c->end_yr;
    }

    dt_cnt = 0;
    day_cnt = 0;
    for (yr = start_yr; yr <= end_yr; yr++) {
        if (is_leap_year(yr)) {
            ndays = 366;
            days_in_month[2] = 29;
        } else {
            ndays = 365;
            days_in_month[2] = 28;
        }
	    for (mth = 1; mth <= 12; mth++) {
	        for (day = 1; day <= days_in_month[mth]; day++) {
	            if (day < 10)
	                sprintf(iday, "0%d", day);
	            else
	                sprintf(iday, "%d", day);

	            if (mth < 10)
	                sprintf(imth, "0%d", mth);
	            else
	                sprintf(imth, "%d", mth);

	            if ((strncmp(met_var, "vph09", 5) == 0) ||
                    (strncmp(met_var, "vph15", 5) == 0)) {

	                sprintf(infname,
	                        "%s/%s/bom-%s-day-%d%s%s-%d%s%s.flt",
	                         c->fdir, met_var, met_var, yr, imth, iday, yr,
	                         imth, iday);
	            } else {
	                sprintf(infname, "%s/%s/%d%s%s_%s.flt",
	                        c->fdir, met_var, yr, imth, iday, met_var);
	            }

                /* read the file file */
                if ((fp = fopen(infname, "r")) == NULL) {
                    fprintf(stderr, "Couldn't open met file %s for read\n",
                            infname);
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }

                if (fread(met_data_day, sizeof(float), c->nrows * c->ncols, fp)
                          != c->nrows * c->ncols) {
                    fprintf(stderr, "Error in reading met file\n");
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
                fclose(fp);

                pixel_count = 0;
                for (k = 0; k < c->num_land_pixels * 2; k+=2) {
                    i = land_ij[k],
                    j = land_ij[k+1];
                    
                    in_offset = i * c->ncols + j;
                    out_offset = day_cnt * c->num_land_pixels + pixel_count;
                    
                    (*met_data)[out_offset] = met_data_day[in_offset];
                    
                    /*
                    if ((strncmp(met_var, "tmax", 3) == 0)) {
                        if ((i == 299) && (j == 321)) {
                            printf("%f %ld\n", met_data_day[in_offset], pixel_count);  
                        }
                    }
                    */
                    pixel_count++;
                }
                
                (*dates)[dt_cnt] = yr;
                (*dates)[dt_cnt+1] = mth;
                (*dates)[dt_cnt+2] = day;
                dt_cnt +=3;
                day_cnt++;
	        }
	    }
    }

    free(met_data_day);
    
    
    
    return ;
}



int  distribute(control *c, int *land_ij, float *met_data, float **my_data,
               int tag, int ndays) {

    int     ii, k, mpi_err;
    int     rest = 0, index;
    long     size_to_send = 0, size_my_data;
    MPI_Status status;

    /* Enough space to fit all data on a process*/
    size_my_data = ndays * (c->nsize+1);

    if ((*my_data = (float *)calloc(size_my_data, sizeof(float))) == NULL ) {
        fprintf(stderr, "Error allocating space for pairs array\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if (c->rank == c->root_processor) {

        /* Divide up the chunks between our available processors */
        if (c->rest > 0) {
            for (k = 0; k < (c->nsize+1) * ndays; k++) {
                /*(*my_data)[k] = slice[k];*/
                (*my_data)[k] = met_data[k];
            }
            index = (c->nsize + 1) * ndays;
            rest = c->rest - 1;
        } else {
            for (k = 0; k < c->nsize * ndays; k++) {
                (*my_data)[k] = met_data[k];
            }
            index = c->nsize * ndays;
            rest = c->rest - 1;
        }

        size_to_send = ndays * (c->nsize+1);
        for (ii = 0; ii < rest; ii++) {
            mpi_err = MPI_Send(&(met_data[index]), size_to_send, MPI_FLOAT, ii+1,
                              tag, MPI_COMM_WORLD);
            if (mpi_err != MPI_SUCCESS) {
                fprintf(stderr, "Error sending MPI slice");
                MPI_Abort(MPI_COMM_WORLD, mpi_err);
            }
            index += size_to_send;
        }

        size_to_send = ndays * c->nsize;
        for (ii = rest; ii < c->size-1; ii++) {
           mpi_err = MPI_Send(&(met_data[index]), size_to_send, MPI_FLOAT, ii+1,
                              tag, MPI_COMM_WORLD);

            if (mpi_err != MPI_SUCCESS) {
                fprintf(stderr, "Error sending MPI slice");
                MPI_Abort(MPI_COMM_WORLD, mpi_err);
            }
            index += size_to_send;
        }
    } else {
        mpi_err = MPI_Recv(*my_data, size_my_data, MPI_FLOAT, MPI_ANY_SOURCE,
                           tag, MPI_COMM_WORLD, &status);
        if (mpi_err != MPI_SUCCESS) {
            fprintf(stderr, "Error in recieving in distribute\n");
            MPI_Abort(MPI_COMM_WORLD, mpi_err);
        }
    }

    if (c->rank < c->rest )
       return size_my_data;
    else
       return size_my_data - ndays;
}





int distribute_ij(control *c, int *land_ij, int **pairs) {

    int  i, mpi_err, index;
    int *position = NULL; /* displacement index for various processors */
    int *send_count = NULL; /* number of elements sent to each processor */
    int  pairs_size = 2 * (c->nsize + 1);

    if ((*pairs = (int *)calloc(pairs_size, sizeof(int))) == NULL) {
        fprintf(stderr,"Error allocating space for pairs array\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if (c->rank == c->root_processor) {
        if ((send_count = (int *)calloc(c->size, sizeof(int))) == NULL) {
	        fprintf(stderr,"Error allocating space for send_count array\n");
		    MPI_Abort(MPI_COMM_WORLD, -1);
        }

        if ((position = (int *)calloc(c->size, sizeof(int))) == NULL) {
	        fprintf(stderr,"Error allocating space for position array\n");
		    MPI_Abort(MPI_COMM_WORLD, -1);
        }

        index = 0;
        for(i = 0; i < c->rest; i++) {
            send_count[i] = 2 * (c->nsize + 1);
            position[i] = index;
            index += 2 * (c->nsize+1);
        }

        for (i = c->rest; i < c->size; i++) {
            send_count[i] = 2 * c->nsize;
            position[i] = index;
            index += 2 * c->nsize;
        }

    }

    /* divide the data among processes */
    mpi_err = MPI_Scatterv(land_ij, send_count, position, MPI_INT, (*pairs),
                           pairs_size, MPI_INT, c->root_processor,
                           MPI_COMM_WORLD);
    if (mpi_err != MPI_SUCCESS) {
        fprintf(stderr, "Error in Scatterv array\n");
        MPI_Abort(MPI_COMM_WORLD, mpi_err);
    }

    free(position);
    free(send_count);

    if (mpi_err) {
        return -1;
    } else {
        if(c->rank < c->rest) {
            return 2*(c->nsize+1);
        } else {
            return 2*c->nsize;
        }
    }
}

void build_radiation_clim(control *c, int *rad_dates, float *rad,
                          float **rad_clim_nonleap, float **rad_clim_leap) {

    /* The radiation data does not overlap (1990-2011) the other met
       forcing, so we need to build a climatology to effectively gap fill for
       the spin up files
    */
    long  date_offset, date_offset2;
    int   doy, yr, year, month, day, ndays;
    int   jan_ndays = 0, feb_ndays = 0, mar_ndays = 0, apr_ndays = 0;
    int   may_ndays = 0, jun_ndays = 0, jul_ndays = 0, aug_ndays = 0;
    int   sep_ndays = 0, oct_ndays = 0, nov_ndays = 0, dec_ndays = 0;
    float jan = 0.0, feb = 0.0, mar = 0.0, apr = 0.0, may = 0.0, jun = 0.0;
    float jul = 0.0, aug = 0.0, sep = 0.0, oct = 0.0, nov = 0.0, dec = 0.0;

    date_offset = 0;
    date_offset2 = 0;
    for (yr = c->start_yr_rad; yr <= c->end_yr_rad; yr++) {
        if (is_leap_year(yr)) {
            ndays = 366;
        } else {
            ndays = 365;
        }

        for (doy = 0; doy < ndays; doy++) {
            year = rad_dates[date_offset];
            month = rad_dates[date_offset+1];
            day = rad_dates[date_offset+2];

            date_offset += 3;

            if (month == 1) {
                jan += rad[date_offset2];
                jan_ndays++;
                /*printf("%f %d %ld\n", rad[date_offset2], jan_ndays, date_offset2);*/
            } else if (month == 2) {
                feb += rad[date_offset2];
                feb_ndays++;
            } else if (month == 3) {
                mar += rad[date_offset2];
                mar_ndays++;
            } else if (month == 4) {
                apr += rad[date_offset2];
                apr_ndays++;
            } else if (month == 5) {
                may += rad[date_offset2];
                may_ndays++;
            } else if (month == 6) {
                jun += rad[date_offset2];
                jun_ndays++;
            } else if (month == 7) {
                jul += rad[date_offset2];
                jul_ndays++;
            } else if (month == 8) {
                aug += rad[date_offset2];
                aug_ndays++;
            } else if (month == 9) {
                sep += rad[date_offset2];
                sep_ndays++;
            } else if (month == 10) {
                oct += rad[date_offset2];
                oct_ndays++;
            } else if (month == 11) {
                nov += rad[date_offset2];
                nov_ndays++;
            } else if (month == 12) {
                dec += rad[date_offset2];
                dec_ndays++;
            }
            date_offset2++;
        }
    }

    for (doy = 0; doy < 365; doy++) {
        if (doy >= 0 && doy < 31) {
            (*rad_clim_nonleap)[doy] = jan / (float)jan_ndays;
        } else if (doy >= 31 && doy < 59) {
            (*rad_clim_nonleap)[doy] = feb / (float)feb_ndays;
        } else if (doy >= 59 && doy < 90) {
            (*rad_clim_nonleap)[doy] = mar / (float)mar_ndays;
        } else if (doy >= 90 && doy < 120) {
            (*rad_clim_nonleap)[doy] = apr / (float)apr_ndays;
        } else if (doy >= 120 && doy < 151) {
            (*rad_clim_nonleap)[doy] = may / (float)may_ndays;
        } else if (doy >= 151 && doy < 181) {
            (*rad_clim_nonleap)[doy] = jun / (float)jun_ndays;
        } else if (doy >= 181 && doy < 212) {
            (*rad_clim_nonleap)[doy] = jul / (float)jul_ndays;
        } else if (doy >= 211 && doy < 243) {
            (*rad_clim_nonleap)[doy] = aug / (float)aug_ndays;
        } else if (doy >= 243 && doy < 273) {
            (*rad_clim_nonleap)[doy] = sep / (float)sep_ndays;
        } else if (doy >= 273 && doy < 304) {
            (*rad_clim_nonleap)[doy] = oct / (float)oct_ndays;
        } else if (doy >= 304 && doy < 334) {
            (*rad_clim_nonleap)[doy] = nov / (float)nov_ndays;
        } else if (doy >= 334 && doy < 365) {
            (*rad_clim_nonleap)[doy] = dec / (float)dec_ndays;
        }
    }
    
    

    for (doy = 0; doy < 366; doy++) {
        if (doy >= 0 && doy < 31) {
            (*rad_clim_leap)[doy] = jan / (float)jan_ndays;
        } else if (doy >= 31 && doy < 60) {
            (*rad_clim_leap)[doy] = feb / (float)feb_ndays;
        } else if (doy >= 60 && doy < 91) {
            (*rad_clim_leap)[doy] = mar / (float)mar_ndays;
        } else if (doy >= 91 && doy < 121) {
            (*rad_clim_leap)[doy] = apr / (float)apr_ndays;
        } else if (doy >= 121 && doy < 152) {
            (*rad_clim_leap)[doy] = may / (float)may_ndays;
        } else if (doy >= 152 && doy < 182) {
            (*rad_clim_leap)[doy] = jun / (float)jun_ndays;
        } else if (doy >= 182 && doy < 213) {
            (*rad_clim_leap)[doy] = jul / (float)jul_ndays;
        } else if (doy >= 212 && doy < 244) {
            (*rad_clim_leap)[doy] = aug / (float)aug_ndays;
        } else if (doy >= 244 && doy < 274) {
            (*rad_clim_leap)[doy] = sep / (float)sep_ndays;
        } else if (doy >= 274 && doy < 305) {
            (*rad_clim_leap)[doy] = oct / (float)oct_ndays;
        } else if (doy >= 305 && doy < 335) {
            (*rad_clim_leap)[doy] = nov / (float)nov_ndays;
        } else if (doy >= 335 && doy < 366) {
            (*rad_clim_leap)[doy] = dec / (float)dec_ndays;
        }
    }

    return;
}


void write_spinup_file(int i, int j, control *c, met *m, float *tmax_ij,
                       float *tmin_ij, float *rain_ij, float *vph09_ij,
                       float *vph15_ij, float *rad_clim_nonleap_ij,
                       float *rad_clim_leap_ij) {

    float  latitude, longitude;
    char   ofname[STRING_LENGTH];
    time_t current_time;
    char*  c_time_string;
    FILE  *ofp;
    long  date_offset;
    int   doy_cnt;
    int   k=0, kk, yr_to_get, st_idx, en_idx, ndays, year;
    float co2=0.0, ndep=0.0, wind_sp=0.0, atpress=0.0, wind_am=0.0;
    float wind_pm=0.0, vpd_avg=0.0, par_day=0.0, sw_am=0.0;
    float Tmean=0.0, Tsoil=0.0, vpd_am=0.0, vpd_pm=0.0;
    float sw_pm=0.0, sw=0.0, rainfall=0.0, day_length;
    float tmin_tomorrow;
    float Tam, Tpm, SEC_TO_DAY, Tavg, sw_w_m2;
    float MJ_TO_J = 1.0 / 1.0E-6;
    float J_TO_UMOL = 4.6;
    float SW_TO_PAR = 0.48;

    /*
    int len_shuffled_yrs = 1;
    int shuffled_yrs[] = {1990};
    */
    
    /*
    int shuffled_yrs[] = {1968,1958,1981,1983,1999,1965,1994,2004,1995,
                          1974,1959,2009,1964,1985,1957,1980,2010,2007,
                          1990,1997,1973,1993,1989,1984,1955,1987,1966,
                          1953,2003,1971,1976,1952,2008,1979,1975,1970,
                          1963,2005,1960,1992,2000,1998,1954,1977,1962,
                          1969,2006,1986,1950,1988};
    int len_shuffled_yrs = 50;
    */
    
    
    int len_shuffled_yrs = 3;
    int shuffled_yrs[] = {1951,1950,1952};
    
    
    sprintf(ofname, "met_data/spinup/met_spinup_%d_%d.csv", i, j);
    ofp = fopen(ofname, "wb");

    latitude = c->yllcorner + (i * c->cellsize);
    longitude = c->xllcorner + (j  * c->cellsize);
    current_time = time(NULL);
    c_time_string = ctime(&current_time);

    fprintf(ofp, "# Daily met: Row:%d x Col:%d\n", i, j);
    fprintf(ofp, "#            Lat:%f x Lon:%f\n", latitude, longitude);
    fprintf(ofp, "# Created by Martin De Kauwe: %s", c_time_string);
    fprintf(ofp, "#--',--,mj/m2/day,c,mm,c,c,c,kPa,kPa,kPa,ppm,t/ha/year,");
    fprintf(ofp, "m/s,kPa,umol/m2/d,m/s,m/s,mj/m2/am,mj/m2/pm\n");
    fprintf(ofp, "year,doy,sw_rad,tair,rain,tsoil,tam,tpm,vpd_am,vpd_pm,");
    fprintf(ofp, "vpd_avg,co2,ndep,wind,atmos_press,par,wind_am,wind_pm,");
    fprintf(ofp, "sw_rad_am,sw_rad_pm\n");

    co2 = 285.0;
    ndep = -9999.9;
    wind_sp = 3.0; /* Haverd et al. 2012 */
    atpress = 100.0; /* 1000 mb -> kPa, Haverd et al. 2012 */
    wind_am = wind_sp;
    wind_pm = wind_sp;

    for (k = 0; k < len_shuffled_yrs; k++) {
        yr_to_get = shuffled_yrs[k];
        st_idx = -999;
        en_idx = -999;
        date_offset = 0;
        for (kk = 0; kk < m->tmax_ndays; kk++) {
            year = m->tmax_dates[date_offset];
            /*printf("*%d %d %ld\n", yr_to_get, year, date_offset);*/
            if (year == yr_to_get) {
                st_idx = kk;
                if (is_leap_year(yr_to_get)) {
                    ndays = 366;
                } else {
                    ndays = 365;
                }

                en_idx = kk + ndays;
                break;
            }
            date_offset += 3;
        }
        doy_cnt = 0;
        for (kk = st_idx; kk < en_idx; kk++) {
            /*printf("**%d %d\n", st_idx, en_idx);*/
            day_length = calc_day_length(kk, ndays, latitude);
            if (kk+1 > en_idx)
                tmin_tomorrow = tmin_ij[kk];
            else
                tmin_tomorrow = tmin_ij[kk+1];

            calc_tam_tpm(&Tam, &Tpm, tmin_ij[kk], tmin_tomorrow,
                         tmax_ij[kk], day_length);

            Tavg = (tmin_ij[kk] + tmax_ij[kk]) / 2.0;
            Tsoil = Tavg;
            Tmean = Tavg;

            /*
            1 MJ m-2 d-1 = 1000000 J m-2 d-1 / 86400 s d-1
                           = 11.574 J m-2 s-1
                           = 11.574 W m-2
            */
            if (ndays == 365)
                sw = rad_clim_nonleap_ij[doy_cnt];
            else
                sw = rad_clim_leap_ij[doy_cnt];
            sw_am = sw / 2.0;
            sw_pm = sw / 2.0;
            sw_w_m2 = sw * 11.574;

            SEC_TO_DAY = 3600. * day_length;
            /*
                Convert radiation from W/m2 -> umol/m2/s (PAR).
                2.3 umol/J for conversion from sw -> PAR (Monteith & Unsworth).

            par_day = sw_w_m2 * 2.3 * SEC_TO_DAY;
            */
            par_day = sw * MJ_TO_J * J_TO_UMOL * SW_TO_PAR;

            vpd_am = calc_vpd(Tam, vph09_ij[kk]);
            vpd_pm = calc_vpd(Tpm, vph09_ij[kk]);
            vpd_avg = (vpd_am + vpd_pm) / 2.0;

            rainfall = rain_ij[kk];

            fprintf(ofp,
            "%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
             year, doy_cnt+1, sw, Tmean, rainfall, Tsoil, Tam, Tpm, vpd_am,
             vpd_pm, vpd_avg, co2, ndep, wind_sp, atpress, par_day, wind_am,
             wind_pm, sw_am, sw_pm);

             doy_cnt++;
        }
        
        
        
    }
    fclose(ofp);

    return;
}


void write_forcing_file(int i, int j, control *c, met *m, float *tmax_ij,
                        float *tmin_ij, float *rain_ij, float *vph09_ij,
                        float *vph15_ij, float *rad_ij,
                        float *rad_clim_nonleap_ij, float *rad_clim_leap_ij) {

    float latitude, longitude;
    char  ofname[STRING_LENGTH];
    time_t current_time;
    char*  c_time_string;
    FILE *ofp;
    
    long date_offset;
    int k=0, kk, jj, yr_to_get, st_idx, en_idx, ndays, doy_cnt, year,st_idx_rad;
    float co2=0.0, ndep=0.0, wind_sp=0.0, atpress=0.0, wind_am=0.0;
    float wind_pm=0.0, vpd_avg=0.0, par_day=0.0, sw_am=0.0;
    float Tmean=0.0, Tsoil=0.0, vpd_am=0.0, vpd_pm=0.0;
    float sw_pm=0.0, sw=0.0, rainfall=0.0, day_length;
    float tmin_tomorrow;
    float Tam, Tpm, SEC_TO_DAY, Tavg, sw_w_m2;

    float MJ_TO_J = 1.0 / 1.0E-6;
    float J_TO_UMOL = 4.6;
    float SW_TO_PAR = 0.48;
    sprintf(ofname, "met_data/forcing/met_forcing_%d_%d.csv", i, j);

    ofp = fopen(ofname, "wb");

    latitude = c->yllcorner + (i * c->cellsize);
    longitude = c->xllcorner + (j * c->cellsize);


    current_time = time(NULL);
    c_time_string = ctime(&current_time);

    fprintf(ofp, "# Daily met: Row:%d x Col:%d\n", i, j);
    fprintf(ofp, "#            Lat:%f x Lon:%f\n", latitude, longitude);
    fprintf(ofp, "# Created by Martin De Kauwe: %s", c_time_string);
    fprintf(ofp, "#--,--,mj/m2/day,c,mm,c,c,c,kPa,kPa,kPa,ppm,t/ha/year,");
    fprintf(ofp, "m/s,kPa,umol/m2/d,m/s,m/s,mj/m2/am,mj/m2/pm\n");
    fprintf(ofp, "year,doy,sw_rad,tair,rain,tsoil,tam,tpm,vpd_am,vpd_pm,");
    fprintf(ofp, "vpd_avg,co2,ndep,wind,atmos_press,par,wind_am,wind_pm,");
    fprintf(ofp, "sw_rad_am,sw_rad_pm\n");


    co2 = 285.0;
    ndep = -9999.9;
    wind_sp = 3.0; /* Haverd et al. 2012 */
    atpress = 100.0; /* 1000 mb -> kPa, Haverd et al. 2012 */
    wind_am = wind_sp;
    wind_pm = wind_sp;


    for (k = c->start_yr; k <= c->end_yr; k++) {
        yr_to_get = k;

        st_idx = -999;
        en_idx = -999;
        date_offset = 0;
        for (kk = 0; kk < m->tmax_ndays; kk++) {
            year = m->tmax_dates[date_offset];
            if (year == yr_to_get) {
                st_idx = kk;
                if (is_leap_year(yr_to_get)) {
                    ndays = 366;
                } else {
                    ndays = 365;
                }

                en_idx = kk + ndays;
                break;
            }
            date_offset+=3;
        }

        /* Need to do the same thing for the rad data - remember array
           length is shorter */
        if (yr_to_get >= 1990) {
            date_offset = 0;
            for (jj = 0; jj < m->rad_ndays; jj++) {
                year = m->rad_dates[date_offset];
                if (year == yr_to_get) {
                    st_idx_rad = jj;
                    break;
                }
                date_offset+=3;
            }
        }


        doy_cnt = 0;
        jj = st_idx_rad;
        for (kk = st_idx; kk < en_idx; kk++) {
            /*printf("%d/%d/%d\n", years[kk], months[kk], days[kk]);*/


            day_length = calc_day_length(kk, ndays, latitude);


            if (kk+1 > en_idx)
                tmin_tomorrow = tmin_ij[kk];
            else
                tmin_tomorrow = tmin_ij[kk+1];

            calc_tam_tpm(&Tam, &Tpm, tmin_ij[kk], tmin_tomorrow,
                         tmax_ij[kk], day_length);

            Tavg = (tmin_ij[kk] + tmax_ij[kk]) / 2.0;
            Tsoil = Tavg;
            Tmean = Tavg;

            /*1 MJ m-2 d-1 = 1000000 J m-2 d-1 / 86400 s d-1
                           = 11.574 J m-2 s-1
                           = 11.574 W m-2 */
           if (year < 1990 && ndays == 365)
               sw = rad_clim_nonleap_ij[doy_cnt];
           else if (year < 1990 && ndays == 366)
               sw = rad_clim_leap_ij[doy_cnt];
           else
               sw = rad_ij[jj];

            sw_am = sw / 2.0;
            sw_pm = sw / 2.0;
            sw_w_m2 = sw * 11.574;

            SEC_TO_DAY = 3600. * day_length;
            /*
                Convert radiation from W/m2 -> umol/m2/s (PAR).
                2.3 umol/J for conversion from sw -> PAR (Monteith & Unsworth).

            par_day = sw_w_m2 * 2.3 * SEC_TO_DAY;
            */


            par_day = sw * MJ_TO_J * J_TO_UMOL * SW_TO_PAR;

            vpd_am = calc_vpd(Tam, vph09_ij[kk]);
            vpd_pm = calc_vpd(Tpm, vph09_ij[kk]);
            vpd_avg = (vpd_am + vpd_pm) / 2.0;

            rainfall = rain_ij[kk];

            fprintf(ofp,
            "%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
             year, doy_cnt+1, sw, Tmean, rainfall, Tsoil, Tam, Tpm, vpd_am,
             vpd_pm, vpd_avg, co2, ndep, wind_sp, atpress, par_day, wind_am,
             wind_pm, sw_am, sw_pm);

            doy_cnt++;
            jj++;
        }

    }

    fclose(ofp);

    return;
}

float calc_day_length(int doy, int yr_days, float latitude) {

    /*

    Daylength in hours

    Eqns come from Leuning A4, A5 and A6, pg. 1196

    Reference:
    ----------
    Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.

    Parameters:
    -----------
    doy : int
        day of year, 1=jan 1
    yr_days : int
        number of days in a year, 365 or 366
    latitude : float
        latitude [degrees]

    Returns:
    --------
    dayl : float
        daylength [hrs]

    */
    float deg2rad, latr, sindec, a, b;

    deg2rad = M_PI / 180.0;
    latr = latitude * deg2rad;
    sindec = -sin(23.5 * deg2rad) * cos(2.0 * M_PI * (doy + 10.0) / yr_days);
    a = sin(latr) * sindec;
    b = cos(latr) * cos(asin(sindec));

    return 12.0 * (1.0 + (2.0 / M_PI) * asin(a / b));
}

void calc_tam_tpm(float *Tam, float *Tpm, float Tmin, float Tmin_tomorrow,
                  float Tmax, float daylength) {
    /*
    The diurnal pattern of air temperature T(t) is calculated from Tmax
    and Tmin on the assumption of a sinusoidal pattern with T = Tmin at
    sunrise and T = (Tmin + Tmax) / 2 at sunset.

    Ross assumed that there was a 3/4 sinusoid in temperature from dawn to
    dusk, ie started at Tmin, went to Tmax 2/3 of the way through the day,
    then decayed to Tav at dusk.

    If Tav = (Tmin + Tmax) / 2 and Tampl = (Tmax - Tmin) / 2 then the
    time course of daytime temperature is described by:

    Tav - Tampl * cos(t)

    where t goes from 0 (dawn) to 3 * pi / 2 (dusk).

    To get morning and afternoon averages you integrate this formula to get
    Tam and Tpm

    Reference
    ----------
    * McMurtie et al (1990) Modelling the Yield of Pinus radiata on a Site
      Limited by Water and Nitrogen. Forest Ecology and Management, 30,
      381-413.
    */
    float Tav, Tampl;

    Tav = (Tmin + Tmax) / 2.0;
    Tampl = (Tmax - Tmin) / 2.0;

    *Tam = Tav - Tampl * (1.0 / sqrt(2.0)) / (3.0 * M_PI / 4.0);
    *Tpm = Tav + Tampl * (1.0 + 1.0 / sqrt(2.0)) / (3.0 * M_PI / 4.0);

    return;
}

float calc_vpd(float temp, float ea) {

    /*
    Empirical equation following Tetens (1930), using the form of Murray
    (1967) as given in Montieth and Unsworth (1990), pg. 10.

    Parameters:
    -----------
    temp : float
        air temperature (deg C)

    References:
    ----------
    * Monteith JL & Unsworth MH (1990) Principles of environmental physics.
    */
    float hPa_2_kPa = 0.1;
    float DEG_TO_KELVIN = 273.15;
    float Tk, A, T_star, T_dash, es_T_star, esat, vpd;

    Tk = temp + DEG_TO_KELVIN;

    A = 17.27;
    T_star = 273.0; /* K */
    T_dash = 36.0;  /* K */
    es_T_star = 0.611; /* kPa */

    /* saturation vapor pressure (kPa)
       Values of saturation vapour pressure from the Tetens formula are
       within 1 Pa of the exact values.
    */
    esat = es_T_star * exp(A * (Tk - T_star) / (Tk - T_dash));

    /* VPD is the saturated vapour pressure - the actual vapour pressure */
    vpd = MAX(0.05, esat - (ea * hPa_2_kPa));

    return vpd;
}

int is_leap_year(int yr) {

    if (yr % 400 == 0 || (yr % 100 != 0 && yr % 4 == 0)) {
        return TRUE;
    } else {
        return FALSE;
    }
}
