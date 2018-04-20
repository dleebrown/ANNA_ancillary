/*
 *  
 * This program will interpolate between lattice points in the
 * Kurucz 1992  model files.  A user provides (1) a temperature, (2) a
 * surface gravity, and (3) a metallicity.  If the input values fall
 * between lattice points then an interpolated model is returned
 * otherwise the lattice point model is returned.  The Kurucz file
 * database is incomplete in that Kurucz did not run all of the available
 * gravity and temperture combinations possible.  Although the official
 * T, g, [Fe/H] bounds for this routine are [3500,10000], [0.0,5.0], &
 * [-3.5,0.5] several gaps exist in the lattice.  I have ID'd most of the
 * holes (I hope) and their existence (or non-existence) is noted in the
 * "key" array below.  The interpolation is a simple (linear) Numerical
 * Recipes routine for 3-D data cubes.  
 *
 *
 * compile this using: gcc -O3 -g mint.c -lm -o mint
 *
 *
 * Created:  Mar. 5, 1996 Modified: Oct. 10, 1996 Alex Stephens
 *
 * BUG -- The low metallicity edge of the lattice is riddled with holes
 * and upon closer inspection some duplicate models.  The new kurucz data
 * files "na???k2.dat" have been edited to remove duplicate models, and
 * holes in the data grid will be filled (soon). 
 *
 * I have pasted together a model for the T/g/Fe = 5750/4.0/-3.5, but
 * this is the only "interpolation" or spackled model.
 *
 * November 16, 1998:  Included the model grid of [Fe/H] = -1.5 which
 * we didn't have on file before.  It was checked for holes.  I then
 * added it to the list of "na???k2.dat" files. Also updated the help.
 *
 */

/*********************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Effective temperature bounds, stepsize, and number of available temps */
#define max_temp 7750
#define min_temp 3500
#define del_temp 250
#define num_temp 1 + (max_temp - min_temp)/del_temp

/* Surface gravity bounds, stepsize, and number of available gravities */
#define max_grav 5
#define min_grav 0
#define del_grav 0.5 
#define num_grav 1 + (max_grav - min_grav)*2

/* Metallicity bounds and total number of available metallicities */
#define max_met +0.5
#define min_met -4.0
#define num_met 16 /* CHANGED +1 11/16/98, +1 12/04/98 */

/* Microturbulence Parameter Bounds */
#define max_vturb 5.00
#define min_vturb 0.00

/* Define the irregular metallicity grid here -- ADD -1.5 11/16/98 */
/* Define the irregular metallicity grid here -- ADD -4.0 12/04/98 */
double 	met_arr[num_met] = {-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, \
	 	            -0.3, -0.2, -0.1,  0.0, +0.1, +0.2, +0.3, +0.5}; 

#define kl	72    /* number of data lines in a Kurucz 91 file */
#define kc	7     /* number of data cols. in a Kurucz 91 file */	

/* database of files that contain the Kurucz 91 model data -- ADD -1.5 11/16/98 */
static char *filelist[] = {"./atm/AM40K2.DAT",  
			"./atm/AM35K2.DAT",  
			"./atm/AM30K2.DAT",  
			"./atm/AM25K2.DAT", 
			"./atm/AM20K2.DAT",  
			"./atm/AM15K2.DAT", 
			"./atm/AM10K2.DAT",  
			"./atm/AM05K2.DAT", 
			"./atm/AM03K2.DAT", 
			"./atm/AM02K2.DAT",  
			"./atm/AM01K2.DAT",  
			"./atm/AP00K2.DAT",  
			"./atm/AP01K2.DAT", 
			"./atm/AP02K2.DAT",
			"./atm/AP03K2.DAT",  
			"./atm/AP05K2.DAT",   
			NULL};

/* This will hold the "key" to the Kurucz file availability */
int	key[num_grav][num_met][num_temp];

/* These hold will hold the 8 models used in the interpolation routine. */
double	mod1[kc][kl], mod2[kc][kl], mod3[kc][kl], mod4[kc][kl];
double	mod5[kc][kl], mod6[kc][kl], mod7[kc][kl], mod8[kc][kl];

/* This is an array that will hold the final, interpolated model */
double  out_matrix[kc][kl];

/* These hold the available temperatures and gravities */
double  temp_arr[num_temp];
double  grav_arr[num_met]; 

double 	temp, grav, met, vturb;  /*  will hold the user input parameters */

char	*program_name;

/****************************************************************************/

/* subroutine definitions */

void interp3D();    /* Interpolate between the 8 models */
void load_grav();   /* Load up the available gravities */
void load_temp();   /* Load up the available temperatures */
void load_model();  /* load an array with kurucz file data */
void init_matrix(); /* initialize an array to all zeros */

void load_key();       /* Load up the array that tell us if a model exists */ 	
int get_depth();      /* Get depth into a model file for a set of t,g,[Fe/H] */
int get_p_filecode(); /* Get the position in the met_arr for a [Fe/H]>0 */
int get_m_filecode(); /* Get the position in the met_arr for a [Fe/H]<=0 */

void write_model();   /* Write the interpolated model to disk */
void print_err();     /* An error message */
void print_err2();    /* Another error message */
void print_instructions(); /* Write instructions if user makes a boo-boo */


/****************************************************************************/
/*		Begin Main Routine: Command Line Input Required             */
/****************************************************************************/

main(int argc, char *argv[])

{
	int 	i, list_num;
	int	dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
	int	tot;
	int	ltpos, htpos, lgpos, hgpos, lmpos, hmpos;
	char	metcode;
	char	modname[80];
	char 	*vturbs, *temps, *gravs, *mets;
	
	// if program doesn't like this size, just increase it
	char 	vturbs2[8];
	char	*out_file = "Default_output";
	double 	low_t = 0 , high_t = 0, low_g = 0, high_g = 0, low_m = 0, high_m = 0;
	double 	t, u, v;

	program_name = argv[0];  /* Save for future use */

	 if (argc <= 4)
       	 print_instructions(program_name); 

	/*** Read off the flags set by the user ***/
	while ((argc > 1) && (argv[1][0] == '-')) {
	 switch (argv[1][1]) {
	   case 't': temps = &argv[1][2];
	   	     temp  = atof(temps);      /* Read the temperature */
	     break;
	   case 'g': gravs = &argv[1][2]; 
		     grav  = atof(gravs);      /* Read the surface gravity */
	     break;
           case 'v': vturbs = &argv[1][2];         /* Read the microturbulence */
	   /* Convert the vturbs to a number, and to a string + E+05  */	
	             vturb  = atof(vturbs);
	 	     strcpy(vturbs2, vturbs); 
		     strcat(vturbs2,"E+05");
             break;
	   case 'w': out_file = &argv[1][2];  /* Get the output filename */
	     break;	
	   case 'p': mets = &argv[1][2]; 
		     met  = 1.0*atof(mets);  /* Read [Fe/H] if >= 0 */
	     metcode = 'p';	
	     break;
	   case 'm': mets = &argv[1][2];
		     met  = -1.0*atof(mets); /* Read [Fe/H] if < 0  */
	     metcode = 'm';
	       if (met == 0.0) {
		printf("You must enter [Fe/H] = 0.0 with the -p flag. \n");
		exit(8);
               }	
	     break;
	   case 'h': print_instructions(program_name);
	     break;
	   default: fprintf(stderr, "Bad flag %s \n", argv[1]);
	     print_instructions(program_name);
         } /* end switch */
        argv++;
        argc--; 
        } /*** end of command line loop ***/

	/*
	 * Echo  Input
	 */
	 printf("\nRequested Model = \n");
	 printf(" - Temperature     = %2.2f \n", temp);
	 printf(" - Surface Gravity = %2.2f \n", grav);
	 printf(" - Metallicity     = %2.2f \n", met);
	 printf(" - Microturbulence = %2.2f \n\n", vturb); 


	/*
	 * Create a default output filename 
	 */
	strcpy(modname, temps); 
	strcat(modname,"_"); 
	strcat(modname,gravs);
	strcat(modname,"_");
	 if (metcode == 'm') {
	  strcat(modname, "m"); 
	 } else {
	  strcat(modname, "p");
	 } 
	strcat(modname, mets); 
	strcat(modname,"_");
	strcat(modname,vturbs);
	strcat(modname,".nm");

      
printf("NAME OF THE DEAULT OUTPUT FILE = %s \n", modname);

 
	/* Begin the routine only if the input values are within the */
	/* temperature, surface gravity, and metallicity bounds      */	
	
    	if ( (temp <= max_temp)   &&  (temp >= min_temp) && \
             (grav <= max_grav)   &&  (grav >= min_grav) && \
             (vturb <= max_vturb) &&  (vturb >= min_vturb) && \
	     (met <= max_met)     &&  (met >= min_met)) {


	load_key(key); /* load up available model martix */


	/* Load up the available temperatures into an array. Also */
	/* note the position of the high and low temps that       */
	/* bracket the user input temperature.			  */

	load_temp(temp_arr);
	ltpos = -1; htpos = -1;

		for (i=0; i<num_temp; i++) {
		  if (temp_arr[i] <= temp)  {
		    low_t = temp_arr[i];
		    ltpos = i;
		       if (temp_arr[i] == temp) {
                         high_t = temp; 
                         htpos = i;
                       } else {
		    high_t = temp_arr[i+1];
		    htpos = i+1;
		    }
		  }
		} 

printf("DEBUG:: Low Temp = %2.2f, High Temp =  %2.2f \n", low_t, high_t);	

	/* Load up the available gravities into an array. Also	*/
	/* note the position in the array of the high and low  	*/
	/* gravity that bracket the user input gravity.		*/

	load_grav(grav_arr);
	lgpos = -1; hgpos = -1;

		for (i=0; i<num_grav; i++) {
                  if (grav_arr[i] <= grav)  {
                    low_g = grav_arr[i];
		    lgpos = i;
		       if (grav_arr[i] == grav) {
		         high_g = grav;
		         hgpos = i; 
		       } else {
                    high_g = grav_arr[i+1];
		    hgpos = i+1;
		    }
                  }
                }
printf("DEBUG:: Low Grav = %2.2f, High Grav = %2.2f \n", low_g, high_g);	

	/* The metallicity array already exist so just scan through */
	/* the array and find the high and low metallicities that   */
	/* brackect the user input metallicity.			    */

	lmpos = -1; hmpos = -1;

		for (i=0; i<num_met; i++) {
		  if (met_arr[i] <= met)  {
                    low_m = met_arr[i];
                    lmpos = i;
                       if (met_arr[i] == met) {
                         high_m = met;
                         hmpos = i;
                       } else {
                    high_m = met_arr[i+1];
                    hmpos = i+1;
                    }
                  }
                }
printf("DEBUG:: Low Met = %2.2f, High Met = %2.2f \n", low_m, high_m);

printf("\n");

/* debugging */

 printf("Low Temp code = %d,  High Temp code = %d \n", ltpos, htpos);
 printf("Low Grav code = %d,  High Grav code = %d \n", lgpos, hgpos);
 printf("Low Met  code = %d,  High Met  code = %d \n", lmpos, hmpos); 
 printf("\n");


	/* Check to see if the chosen values have high & low models */
	/* Need the braketing models to do the interpolation        */

	tot = key[lgpos][lmpos][ltpos] + key[lgpos][lmpos][htpos] + \
	      key[lgpos][hmpos][ltpos] + key[lgpos][hmpos][htpos] + \
	      key[hgpos][lmpos][ltpos] + key[hgpos][lmpos][htpos] + \
	      key[hgpos][hmpos][ltpos] + key[hgpos][hmpos][htpos];	

	if (tot == 8)  {  /* we have all of the models so start the interp. */
	
	printf("%s -- All %i models available:\n\n", program_name, tot);	

	/* The next 8 block do the following:  Since the bracketing T, g,   */
	/* and [Fe/H] values are known, we need to extract these bracketing */
	/* models from the Kurucz files.  So these blocks locate the place  */
	/* inside the Kurucz file where the model resides, it then load a   */
	/* 2-D array with the rhox,T,Pg,Ne,accrad,abross,&vturb model data  */ 


	/** (1) -- grab the info for the low_t, low_g, low_m model **/	
	 if (metcode == 'm') list_num = get_m_filecode(low_m);
	  else if (metcode == 'p') list_num = get_p_filecode(low_m); 

	printf("DEBUG::%d -- get_p_filecode \n", list_num);	

		dn1 = get_depth(lgpos, lmpos, ltpos, key);
	 	init_matrix(mod1); /* set all values to zero */
	 	load_model(list_num, dn1, mod1);

	/** (2) -- grab the info for the high_t, low_g, low_m model **/
         if (metcode == 'm') list_num = get_m_filecode(low_m);
          else if (metcode == 'p') list_num = get_p_filecode(low_m);
	  	dn2 = get_depth(lgpos, lmpos, htpos, key);
         	init_matrix(mod2); /* set all values to zero */
         	load_model(list_num, dn2, mod2);

        /** (3) -- grab the info for the high_t, high_g, low_m model **/
         if (metcode == 'm') list_num = get_m_filecode(low_m);
          else if (metcode == 'p') list_num = get_p_filecode(low_m);
         	dn3 = get_depth(hgpos, lmpos, htpos, key);
         	init_matrix(mod3); /* set all values to zero */
         	load_model(list_num, dn3, mod3);

        /** (4) -- grab the info for the low_t, high_g, low_m model **/
         if (metcode == 'm') list_num = get_m_filecode(low_m);
          else if (metcode == 'p') list_num = get_p_filecode(low_m);
         	dn4 = get_depth(hgpos, lmpos, ltpos, key);
         	init_matrix(mod4); /* set all values to zero */
         	load_model(list_num, dn4, mod4);

        /** (5) -- grab the info for the low_t, low_g, high_m model **/
         if (metcode == 'm') list_num = get_m_filecode(high_m);
          else if (metcode == 'p') list_num = get_p_filecode(high_m);
         	dn5 = get_depth(lgpos, hmpos, ltpos, key);
         	init_matrix(mod5); /* set all values to zero */
         	load_model(list_num, dn5, mod5);

        /** (6) -- grab the info for the high_t, low_g, high_m model **/
         if (metcode == 'm') list_num = get_m_filecode(high_m);
          else if (metcode == 'p') list_num = get_p_filecode(high_m);
         	dn6 = get_depth(lgpos, hmpos, htpos, key);
         	init_matrix(mod6); /* set all values to zero */
         	load_model(list_num, dn6, mod6);

        /** (7) -- grab the info for the high_t, high_g, high_m model **/
         if (metcode == 'm') list_num = get_m_filecode(high_m);
          else if (metcode == 'p') list_num = get_p_filecode(high_m);
         	dn7 = get_depth(hgpos, hmpos, htpos, key);
         	init_matrix(mod7); /* set all values to zero */
         	load_model(list_num, dn7, mod7);

        /** (8) -- grab the info for the low_t, high_g, high_m model **/
         if (metcode == 'm') list_num = get_m_filecode(high_m);
          else if (metcode == 'p') list_num = get_p_filecode(high_m);
         	dn8 = get_depth(hgpos, hmpos, ltpos, key);
         	init_matrix(mod8); /* set all values to zero */
         	load_model(list_num, dn8, mod8);

	
	/** Have to make sure we dont divide by zero when interpolating **/

	/* t is the fractional distance between the low and high temp */	
		if (low_t == high_t)  t = 0; 
	 	else  t = (temp - low_t)/(high_t - low_t);
	
	/* u is the fractional distance between the low and high grav */
		if (low_g == high_g)  u = 0; 
	 	else  u = (grav - low_g)/(high_g - low_g); 

	/* v is the fractional distance betwee the low and high [Fe/H] */	
        	if (low_m == high_m)  v = 0;
         	else  v = (met - low_m)/(high_m - low_m);

/* debugging */
/*
printf("User met = %f, Low Met = %f, High Met = %f \n", met, low_m, high_m);
printf("fractions t=%f u=%f v=%f\n", t, u, v);
*/

	/* We have the models and the fractional distances: So interpolate. */  

 	interp3D(mod1, mod2, mod3, mod4, \
		 	mod5, mod6, mod7, mod8, \
		 		t, u, v, out_matrix); 
	
	/* Write the interpolated file to disk */

 	write_model(out_matrix, temp, grav, met, vturbs2, out_file); 	
	printf("\nKurucz/MOOG model atm. written to --> %s\n\n", out_file);

	}  /** end of if statement that began with "if (tot == 8)" **/ 
	else {         
	 print_err2(); /* if dont have all models write an err. msg. */
	}
	
	} /*** end of if that began with parameter boundary test  ***/  
	else {
	 print_err(); /* the input was out-of-bounds so write error msg. */
	}
exit(1); 

}  /****  End of Main ****/


/*****************************************************************************/
/*---------------------------------------------------------------------------*/
/*****************************************************************************/
/*---------------------------------------------------------------------------*/
/*****************************************************************************/


 		 /**** All Subroutines are defined below. ****/


/*****************************************************************************/
/*				init_matix				     */
/*									     */
/*  Initialize a model matrix by setting all values equal to zero            */
/*									     */
/*  Input (inti_matrix) = the array to be initialized.	                     */
/*  Output = the zeroed array 					             */
/* 									     */	
/*****************************************************************************/

void init_matrix(double data[][kl])
{
	int	i,j;

	for (i=0; i<kc; i++) {
	  for (j=0; j<kl; j++)
	  data[i][j] = 0.0; 
	}

} /* end of subroutine init_matix */


/*****************************************************************************/
/*                              load_temp                                    */
/*									     */
/*  Subroutine that loads up the available model temperatures into an array. */
/*									     */
/*  Input (temp_in) = the temperature array.			             */
/*  Output = the temperature array with the model values in it. (duh).	     */
/*									     */ 
/*****************************************************************************/

void load_temp(double temp_in[num_temp])
{
	int	i;

	for (i=0; i<num_temp; i++) {
	  temp_in[i] = min_temp + (i*del_temp);
	}

} /* end of load_temp */


/*****************************************************************************/
/*                              load_grav                                    */
/*								             */	
/*  Subroutine that loads up the available model gravities into an array.    */
/*                                                                           */
/*  Input (grav_in) = the gravity array.                                     */
/*  Output = the gravity array with the model values in it. (duh).           */
/*                                                                           */
/*****************************************************************************/

void load_grav(double grav_in[num_grav])
{
	int	i;

	for (i=0; i<num_grav; i++) {
 	  grav_in[i] = (min_grav/1.0) + (i*del_grav);
	}

}  /* end of load_grav */ 


/*****************************************************************************/
/*				load_model				     */
/*									     */
/*  Subroutine to load up a model array. The model array (kc by kl elements) */
/*  is passed to the subroutine.  The subroutine will fill the array with    */
/*  the relevant parameters rhox, Temp, Pgas   , e- density, rosseland       */
/*  absorption coefficient, radiative acceleration, and a microturbulent     */
/*  velocity from a data  file that matches the Teff, g, and [Fe/H] for      */
/*  that particular model.						     */
/*									     */
/*  Input (file_code) =  the name of the "na??k2.dat" data file to be used in */
/*  	 loading up the model array.  (mod_num) = the number of files into   */
/*	 "a??k2.dat" where I can find the model that corresponds to a bound  */
/*       of the user input data (i.e. a model with a particular [Fe/H], g,   */
/*       Teff that bounds the user input g, Teff, [Fe/H]).  (model) =  the   */
/*       model array to be filled.					     */
/*  Output = the loaded model array.					     */
/*								             */
/*****************************************************************************/

void load_model(int file_code, int mod_num, double model[][kl])
{

        int     k = 0;
        int     thresh, correct, counter; 
	char	line[255];
	char    *str_ptr;
	char	FILE_NAME[255];
	FILE    *in_file;
	double  rhox[kl], T[kl], P[kl], Xne[kl];
	double  abross[kl], accrad[kl], vturb[kl];
	double  dummy1[kl], dummy2[kl];

	/* Search for the Kurucz model file of interest */
 	strcpy(FILE_NAME, filelist[file_code]);	

	/* Check to see if file exists */
        in_file = fopen(FILE_NAME, "r");
            if (in_file == NULL) {
              printf("Cant open the file -> %s.\n", FILE_NAME);
              exit(8);
            }

	counter = 1;

	/* If the file is there, read the data into arrays */
        while(1) {

	 /* This "thresh" represents the number of lines down into a model */
	 /* file -- where the model of interest actually starts. There are */
	 /* 97 lines per model, and the model is mod_num files down.       */

	 thresh = ((mod_num-1)*97); /* subtract 1 b/c numbering starts @ 1*/

	   /* skip over all of the data before the relevant part */
	   /* There are 24 lines of useless information.  Also   */
	   /* stop after the end of the file is read.            */

	   if ( (counter <= (thresh+24)) || (counter >= (thresh+97)) ) {
             str_ptr = fgets(line, sizeof(line), in_file);
	           if (counter >= (thresh+97)) {
       		   break;
		   }
	     counter++;
	   } else {  

	   /* If have reached the relevant data then extract */
	   str_ptr = fgets(line, sizeof(line), in_file);
           correct = sscanf(line, "%lE %lE %lE %lE %lE %lE %lE %lE %lE", \
                           &rhox[k], &T[k], &P[k], &Xne[k], \
                            &abross[k], &accrad[k], &vturb[k], &dummy1[k], &dummy2[k]);
	   
           /* If hit the EOF exit */
           if (str_ptr == NULL) 
            break;
      
	  /* For clarity assign the atmospheric parameters to the */
	  /* model here instead of burying it in the above lines  */

	    model[0][k] = rhox[k]; 
	    model[1][k] = T[k];
	    model[2][k] = P[k];
	    model[3][k] = Xne[k];
	    model[4][k] = abross[k];
	    model[5][k] = accrad[k];
	    model[6][k] = vturb[k]; 

	    counter++;
	   k++;
	} /* end of the else */
      
      } /* end of the while loop */
        
  fclose(in_file);
  printf("File %s closed. Read %d lines \n", FILE_NAME, k);

} /* end of subroutine load_model */


/*****************************************************************************/
/*				write_model				     */
/*									     */
/*  Write out a MOOG readable atmosphere file to disk.			     */
/*									     */
/* MOOG/KURUCZ files have the following format:				     */
/*	line 1:  KURUCZ							     */
/*	line 2:  Kurucz 1991: Teff, lg[g], [A/H]			     */
/*      line 3:  NTAU           kl				             */
/*	line 4-68: data from the interpolation routines			     */
/*      line 69: microturbulence parameter (assumed = 1.0e+05 km/s)	     */
/*      line 70: NATOMS           1   (metallicty)			     */
/* 	line 71:       4.00     1.11 (the atom of interest and an abund)     */
/*	line 72:						             */	
/*****************************************************************************/

void write_model(double model[][kl], double temp_par, double grav_par, \
                  double met_par, char *vturb_str, char *filename)
{

	int	j;
	char	line1[255], line2[255], line3[255], line4[255], line5[255];
	char	line6[255], line7[255], line8[255], line9[255], line10[255];
	char	line11[255];
	char    line12[255], line13[255];
	FILE	*out_file;

	printf("OUTPUT FILE = %s \n", filename);

	/* Check to see that we can write to the file. */	
        
	out_file = fopen(filename, "w");
	if (out_file == NULL) {
                printf("Cant open the file -> %s.\n", filename);
                exit(8);
        }

	/*********************************************************/
	/* Header Information for a MOOG/Kurucz Atmosphere file. */
	/*********************************************************/
	
	strcpy(line1,"KURUCZ");
        strcpy(line2,"#Kss72: T=");
	strcpy(line3,"[g]=");
	strcpy(line4,"[Fe/H]=");
	strcpy(line11,"vt=");
        strcpy(line5,"NTAU            72");
	strcpy(line6,"     ");
	strcat(line6, vturb_str);
	strcpy(line7,"NATOMS           1 "); 
  	strcpy(line8,"      3.00     3.30"); /* Meteoritic Li Abundance */
	strcpy(line9,"NMOL             19");
	strcpy(line10,"     606.0     106.0     607.0     608.0      107.0     108.0     112.0 707.0");
	strcpy(line12,"     708.0     808.0      12.1   60808.0    10108.0     101.0       6.1 7.1");
	strcpy(line13,"       8.1     822.0      22.1");
	
	/* Write the headers to a file */
	fprintf(out_file, "%s\n", line1);
	fprintf(out_file, "%s %4.0f,%s%1.2f,%s%1.2f,%s%s\n", \
	                   line2, temp_par, line3, grav_par, line4, \
			   met_par, line11, vturb_str);
	fprintf(out_file, "%s\n", line5);
		
	/* Write out the model information to a file:  Again for a kurucz */
	/* model index 0=rhox, 1=temp, 2=Pgas, 3=Xne, 4=abross, 5=accrad  */
	/* 6=vturb.  This info will eventually be the result of an interp */

	for (j=0; j<kl; j++)  {
	  fprintf(out_file, " %1.8E   %4.1f %1.3E %1.3E %1.3E %1.3E %1.3E \n", \
	     model[0][j], model[1][j], model[2][j], model[3][j], model[4][j], \
             model[5][j], model[6][j]);
	}

	/* write out MOOG trailing information to the file */	
	fprintf(out_file, "%s\n", line6);
	fprintf(out_file, "%s %1.2f\n", line7, met_par);
	fprintf(out_file, "%s\n", line8);
	fprintf(out_file, "%s\n", line9);
	fprintf(out_file, "%s\n", line10);
	fprintf(out_file, "%s\n",line12);
	fprintf(out_file, "%s\n",line13);
	
	fclose(out_file);
 
} /* end of subroutine write_model */


/*****************************************************************************/
/*				interp3D				     */
/*									     */
/*  Interpolate between points using the numerical recipes formulation.  In  */
/*  other words do a simple linear interpolation between stellar parameter   */
/*  lattice points using the model files and the user input parameters.      */
/*									     */
/*  Input (m*[][]) = the model arrays that contain the data that are the     */
/*  	 	stellar parameter lattice points (in T, g, [Fe/H] space).    */
/*  	 	(T, U, V) = the fractional difference between one lattice    */
/*  		point in one dimension (T for temp, U for grav, V for met)   */
/*  Output (output)  = an array that contains the interpolated information.  */
/*									     */
/*****************************************************************************/

void interp3D(double m1[][kl], double m2[][kl], double m3[][kl], \
		       double m4[][kl], double m5[][kl], double m6[][kl], \
			double m7[][kl], double m8[][kl], double T, double U, \
		         double V, double output[][kl]) 
{
	int	i, j;

	/******************************************************/
	/* Found that the linear interpolation way is the     */
	/* preferable way to determine the interpolation.     */
	/* First transforming from the data value to the      */
	/* logarithm of the data value and then back produces */
	/* spurious results.	                              */
	/******************************************************/
			
	for (j=0; j<kc; j++) {
	 for (i=0; i<kl; i++) {

	/* Interpolation formula sim to NR pg. 104 : C manual */
	 output[j][i] =   ((1-T)*(1-U)*(1-V)*m1[j][i]) + \
	 	              (T*(1-U)*(1-V)*m2[j][i]) + \
			          (T*U*(1-V)*m3[j][i]) + \
			      ((1-T)*U*(1-V)*m4[j][i]) + \
		              ((1-T)*(1-U)*V*m5[j][i]) + \
			          (T*(1-U)*V*m6[j][i]) + \
			  	      (T*U*V*m7[j][i]) + \
			          ((1-T)*U*V*m8[j][i]); 
	 }
	}

} /* end of subroutine interp3D */


/*****************************************************************************/
/*				load_key				     */
/*									     */
/*  Using information regarding the completeness of the kurucz stellar atm.  */
/*  data files, this routine will take a 3D array with dimensions equal to   */
/*  the number of gravs, mets, and temps available.  So each array point     */
/*  represents a 'model' as it were with a specific set of stellar params.   */
/*  this routine goes to each point in the 3D array and sets the value at    */
/*  that point = 1 if the model for that point exists, otherwise the value   */
/*  at that points is set to 0.  This key array come in handy when trying    */
/*  to figure out the position of a particular model within a Kurucz         */
/*  ap???k2.dat file.						             */
/*                                                                           */
/*  Input (cube_in) = the 3D array.  Output = the loaded array.              */
/*                                                                           */
/*****************************************************************************/

void load_key(int cube_in[num_grav][num_met][num_temp])
{
	int i,j,k;

/* 
*

Need to load up the key array which is basically a map of the holes in
the Kurucz model atmosphere grid.  When one searched through the low
metallicity end of the grid one finds that there are several (1) typos
(that's troubling) (2) duplicate models where there shouldnt be any
and (3) just basic omissions.  So the key reflects my search through
the grid and a '0' in the key array stands for an omission, etc. while
a '1' represents an available model.  Since I have had to go throuhg
this more than once I include below a key which shows how the array
subscripts correspond to  T, g, and [Fe/H] values in the Kurucz grid.

*
*/

/*****************************************************************************/
/*

                    Array nomenclature:
		    (ADD -1.5=Fe on 11/16/98, reordering j array)

     j = num_met subscript: 	      i = num_grav subscript:
    0 = -4.0
    1 = -3.5	; 9  = -0.2		0 = 0.0	;  6 = 3.0
    2 = -3.0	; 10 = -0.1		1 = 0.5	;  7 = 3.5
    3 = -2.5 	; 11 =  0.0		2 = 1.0 ;  8 = 4.0
    4 = -2.0 	; 12 = +0.1		3 = 1.5 ;  9 = 4.5
    5 = -1.5 <----- NEW ADDITION
    6 = -1.0	; 13 = +0.2 		4 = 2.0	; 10 = 5.0
    7 = -0.5	; 14 = +0.3		5 = 2.5
    8 = -0.3	; 15 = +0.5

             k = num_temp subscript
    0 = 3500	;  9 = 5750	; 18 = 8000
    1 = 3750 	; 10 = 6000	; 19 = 8250
    2 = 4000	; 11 = 6250	; 20 = 8500
    3 = 4250 	; 12 = 6500	; 21 = 8750
    4 = 4500	; 13 = 6750	; 22 = 9000
    5 = 4750    ; 14 = 7000	; 23 = 9250
    6 = 5000	; 15 = 7250	; 24 = 9500
    7 = 5250 	; 16 = 7500	; 25 = 9750
    8 = 5500	; 17 = 7750	; 26 = 10,000

 NOTE NOTE NOTE: FOR MSPAWN72, ONLY GO UP TO 8000 K or 18.
*/
/*****************************************************************************/


	/* step through each point in the array */

	for (i=0; i<num_grav; i++) 
	 for (j=0; j<num_met; j++)
	  for (k=0; k<num_temp; k++) {


           /* if the model for a particular i,j,k exists then set the */
	   /* value at i,j,k == 1:  otherwise set to zero (model DNE) */

	   /* default is to set the key = 1, as if a model exists */ 
	   cube_in[i][j][k] = 1;

	    /*** general trends in the model grid ***/

	    /* No models when T>=6250 and g=0.0 */
	    if (k>=11 && i==0) cube_in[i][j][k] = 0;

	    /* No models when T>=7750 and g=0.5 */
	    if (k>=17 && i==1) cube_in[i][j][k] = 0;

	    /* No models when T>=8750 and g=1.0 */
	    if (k>=21 && i==2) cube_in[i][j][k] = 0;

	    /* No models when T>=9250 and g=1.5 */
	    if (k>=23 && i==3) cube_in[i][j][k] = 0;

            /********************************************/
	    /*** Special cases that happen to crop up ***/
	    /********************************************/

	    /*  Special [Fe/H] = +0.5 cases */  
	    if (k>=15 && i==1  && j==15)  cube_in[i][j][k] = 0;
	   /* if (k>=20 && i==2  && j==15)  cube_in[i][j][k] = 0; */
	    
	    /*  Special [Fe/H] = +0.3 cases */  
	    if (k>=16 && i==1  && j==14)  cube_in[i][j][k] = 0;
	   
	    /*  Special [Fe/H] = -1.0 cases  !!! ONLY CHECKED UP TO k=18 */  
	    /* if (k>=20 && i==2  && j==6)   cube_in[i][j][k] = 0; */
     	   
            /*  Special [Fe/H] = -1.5 cases */
           /* if (k>=19 && i==2  && j==5)   cube_in[i][j][k] = 0; */
           /* if (k>=21 && i==3  && j==5)   cube_in[i][j][k] = 0; */
 
	    /*  Special [Fe/H] = -2.0 cases */  
	    /* if (k>=19 && i==2  && j==4)   cube_in[i][j][k] = 0; */
	    /* if (k>=21 && i==3  && j==4)   cube_in[i][j][k] = 0; */

	    /*  Special [Fe/H] = -2.5 cases */	    
	    if (k==0  && i==9  && j==3)   cube_in[i][j][k] = 0; 
	    if (k==0  && i==10 && j==3)   cube_in[i][j][k] = 0;
	    if (k==1  && i==10 && j==3)   cube_in[i][j][k] = 0;
	    if (k==2  && i==10 && j==3)   cube_in[i][j][k] = 0;
	    if (k>=18 && i==2  && j==3)   cube_in[i][j][k] = 0;
	    /* if (k>=21 && i==3  && j==3)   cube_in[i][j][k] = 0; */

	    /*  Special [Fe/H] = -3.0 cases */
	    if (k==0  && i==9  && j==2)   cube_in[i][j][k] = 0;
	    if (k==0  && i==10 && j==2)   cube_in[i][j][k] = 0;
	 /* if (k==1  && i==9  && j==2)   cube_in[i][j][k] = 0; */
	 /* if (k==1  && i==10 && j==2)   cube_in[i][j][k] = 0; */
	    if (k==6  && i==2  && j==2)   cube_in[i][j][k] = 0;
         /* if (k==7  && i==1  && j==2)   cube_in[i][j][k] = 0; */
/* HOLE FILLED WITH INTERPOLATED MODEL -- PRINT WARNING */
/* if (k==9  && i==8  && j==1)   cube_in[i][j][k] = 0;  */
/* if (k==9  && i==8  && j==1)  printf("WARN: 5750,4.0,-3.0 model is a filled Kurucz hole.\n");*/
	    if (k>=16 && i==1  && j==2)   cube_in[i][j][k] = 0;
	    if (k>=17 && i==2  && j==2)   cube_in[i][j][k] = 0;
         /* if (k>=21 && i==3  && j==2)   cube_in[i][j][k] = 0; */


	    /*  Special [Fe/H] = -3.5 cases */	
	    if (k==0  && i==8  && j==1)   cube_in[i][j][k] = 0;
	    if (k==0  && i==9  && j==1)   cube_in[i][j][k] = 0;
	    if (k==0  && i==10 && j==1)   cube_in[i][j][k] = 0;
            if (k==1  && i==9  && j==1)   cube_in[i][j][k] = 0;
	    if (k==1  && i==10 && j==1)   cube_in[i][j][k] = 0;
            if (k==2  && i==10 && j==1)   cube_in[i][j][k] = 0;
	    if (k==3  && i==10 && j==1)   cube_in[i][j][k] = 0;
	    if (k>=16 && i==1  && j==1)   cube_in[i][j][k] = 0;
	    if (k>=17 && i==2  && j==1)   cube_in[i][j][k] = 0;
	   /*  if (k==3  && i==9  && j==1)   cube_in[i][j][k] = 0; */
	   /*  if (k==7  && i==7  && j==1)   cube_in[i][j][k] = 0; */
	   /*  if (k>=21 && i==3  && j==1)   cube_in[i][j][k] = 0; */


             /*  Special [Fe/H] = -4.0 cases */
            if (k==0  && i==7  && j==0)   cube_in[i][j][k] = 0;
            if (k==0  && i==8  && j==0)   cube_in[i][j][k] = 0;
            if (k==0  && i==9  && j==0)   cube_in[i][j][k] = 0;
            if (k==0  && i==10 && j==0)   cube_in[i][j][k] = 0;
            if (k==1  && i==8  && j==0)   cube_in[i][j][k] = 0;
            if (k==1  && i==9  && j==0)   cube_in[i][j][k] = 0;
            if (k==1  && i==10 && j==0)   cube_in[i][j][k] = 0;
            if (k==2  && i==10 && j==0)   cube_in[i][j][k] = 0;
	    if (k>=15 && i==1  && j==0)   cube_in[i][j][k] = 0;
	    if (k>=17 && i==2  && j==0)   cube_in[i][j][k] = 0;

	   /* Two special insertions */
	    if (k==22 && i==3  && j==2)   cube_in[i][j][k] = 1;
	    if (k==22 && i==3  && j==1)   cube_in[i][j][k] = 1;  

	 }
	
} /* end of subroutione load_key */


/****************************************************************************/
/* 				check_bounds			            */
/*                                                                          */
/*  This subroutine finds the position in the stellar parameter array that  */
/*  matches the position of the passed low and high stellar parameter.      */
/*									    */
/*  Input (low_val, high_val, num_val) the low and high bounds and the # of */
/*          elements in a particular value array. (val_arra) = the array    */
/*          that contains the stellar data.                                 */
/*  Output (lv_pos, hv_pos) = the low and high postions of the values that  */
/*	    bound the user input. 					    */
/*                                                                          */
/****************************************************************************/

int check_bounds(double low_val, double high_val, int num_val, \
		     double val_arr[], int lv_pos, int hv_pos) 
{
        int	count = 0;
        int	ll = 0;
        int	hh = 0;
      
        /* Check to see if the chosen values have high & low models */
	
	while (((ll == 0) || (hh == 0)) && (count < num_val)) {
         
	    if (low_val == val_arr[count]) {
               lv_pos = count;
	        printf("%d %d \n", count, lv_pos);
                ll = 1; /* set a switch */
               printf("dingl %f \n", val_arr[count]);
            }
	    
            if (high_val == val_arr[count]) {
               hv_pos = count;
	        printf("%d %d \n", count, hv_pos);
                hh = 1; /* set a switch */
               printf("dingh %f \n", val_arr[count]);
             }
        count++;
        }

 return(lv_pos, hv_pos);

} /* end of subroutine check_bounds */

/****************************************************************************/
/*				get_depth				    */
/*                                                                          */
/*  Subroutine that will take the 3D matrix 'key' that tells us whether or  */
/*  not a model of a certain temp, grav, and metallicity exists (i.e. a     */
/*  key = 1 if a model with a [grav][met][temp] exisits, otherwise key =0   */
/*  at that position) and then tell us the position of the model of interest*/
/*  in a file that holds the models for a given metallicity.  Here is an    */
/*  Say we want the model with [g][m][t] = [0.5][-0.1][3500].  Each of the  */
/*  Input parameters represents a position in the grav_arr, met_arr, and    */
/*  temp_arr.  Knowing the metallicity means we know which model file to    */
/*  search (i.e. apm01k2.dat) but we need to know where the model with g=0.5*/
/*  and t=3500 is within this file.  This is why key exists.  What I do is  */
/*  assign a temperature, then step through all the gravities for that temp */
/*  adding up the ones and zeros that correspond to a file (or no file)     */
/*  until I reach the temperature and gravity in question.  The sum of 1s   */
/*  is the number of models into the file where the model with my g and temp*/
/*  rests.								    */
/*								            */
/*  Input (gpos, mpos, tpos) = The position in the grav_arr, met_arr, and   */
/*		      temp_arr that represents the g, T, [Fe/H]. Also the   */
/*		      key matrix (key).					    */
/*  Output (depth) =  The number of models into a file where the model      */
/*	              of interest resides.				    */ 
/*                                                                          */
/****************************************************************************/

int get_depth(int gpos, int mpos, int tpos, \
	       int key_arr[num_grav][num_met][num_temp]) 
{
	int 	i,j;

	int	depth = 0;

	for (i=0; i<=tpos; i++) {
	 for (j=0; j<num_grav; j++) {
 		
		/* do nothing once exceed the grav and temp position */
		/* otherwise add up the values as you go	     */ 

 		 if ((i == tpos) && (j > gpos)) { 
		 } else {
		 depth = (key_arr[j][mpos][i] + depth);
		 } 
	 }
	}	
	
	/* echo the depth into the model file */
	printf("Found model # %i\n", depth);
	
  return(depth);

} /* end of get_depth subroutine */


/****************************************************************************/
/*				get_p_filecode 			            */
/*                                                                          */
/*  Subroutine that will take an input metallicity bound (low or high) and  */
/*  will return the code (i.e. the number of the value in globally defined  */
/*  filelist) which represents the names of the 2MB+ files that contain all */
/*  of the models available for a certain metallicity. Probably a lame way  */
/*  accessing a filename, but simple. For positive metallicities only.      */
/*                                                                          */
/*  Input (param) = the value of the metallicity bound.                     */
/*  Output (code) = the position in the filelist array that equals a file   */
/*                  which has the (param) metallicity value in its filename */
/*                                                                          */
/****************************************************************************/

int get_p_filecode(double param) 
{
        int     code;

        /* for positive metallicities */
             if (param == 0.5)  code = 15;
        else if (param == 0.3)  code = 14;
        else if (param == 0.2)  code = 13;
        else if (param == 0.1)  code = 12;
        else if (param == 0.0)  code = 11;
        
	/* print an error message if the file doesnt match an    */
        /* available model file.  Shouldn't need resort to this  */
        /* trap door because bad values should be filtered out   */
        /* long before this function.                            */

	else {
             printf("Error: the requested file does not exist. \n");
             exit(8);
	     }
 
  return code;

} /* end of get_p_filecode subroutine */


/****************************************************************************/
/*				get_m_filecode 			            */
/*									    */
/*  Subroutine that will take an input metallicity bound (low or high) and  */
/*  will return the code (i.e. the number of the value in globally defined  */
/*  filelist) which represents the names of the 2MB+ files that contain all */
/*  of the models available for a certain metallicity. Probably a lame way  */
/*  accessing a filename, but simple. For negative metallicities only.	    */
/*									    */
/*  Input (param) = the value of the metallicity bound.			    */
/*  Output (code) = the position in the filelist array that equals a file   */
/*		    which has the (param) metallicity value in its filename */ 
/*									    */
/****************************************************************************/
 
int get_m_filecode(double param) 
{
        int     code;

        /* for negative metallicities */
	     if (param ==  0.0) code = 10;
	
	/* the above line is a cheap fix (i hope) for when a user inputs  */
	/* a negative metallicity -0.1 < [Fe/H] < 0.0.   The upper limit  */
	/* for this star is _actually_ +0.0 or p00 in the filename scheme */
	/* this way I call to the negative metallicty routine but can     */
	/* still assign a positive value when the upper bound happens > 0 */
       
	/* ADDED -1.5 FILE and reordered list on 11/16/98 */ 
	else if (param == -0.1) code = 10;
        else if (param == -0.2) code = 9;
        else if (param == -0.3) code = 8;
        else if (param == -0.5) code = 7;
        else if (param == -1.0) code = 6;
        else if (param == -1.5) code = 5;
        else if (param == -2.0) code = 4;
        else if (param == -2.5) code = 3;
        else if (param == -3.0) code = 2;
        else if (param == -3.5) code = 1;
        else if (param == -4.0) code = 0;


	/* print an error message if the file doesnt match an    */
        /* available model file.  Shouldn't need resort to this  */
        /* trap door because bad values should be filtered out   */
        /* long before this function.                            */

        else {
             printf("Error: the requested file does not exist. \n");
	     exit(8);
	     }
 
 return code;

} /* end of get_m_filecode subroutine */


/****************************************************************************/
/*				print_err 			            */
/*									    */
/*   Subroutine that prints an error message when the user inputs a temp,   */
/*   gravity, or metallicity that is not within the acceptable range        */
/*									    */
/*   No input or output required.					    */
/*									    */ 
/****************************************************************************/

void print_err(void)
{
	
	printf("\n");
	printf("*********************************************************\n");
	printf("You have entered a value that exceeds the program bounds.\n");
	printf("*********************************************************\n");
	printf("\nBounds: \n");
	printf("	Temp = [3500,7750] \n");
	printf("	Grav = [0.00,5.00] \n");
	printf("	Met  = [-4.0,+0.5] \n");
	printf("\n\n");
	exit(8);

} /* end of print_err subroutine */


/****************************************************************************/
/*                        	print_err2			            */
/*									    */
/*  Subroutine prints an error message if the user inputs a valid set of    */
/*  parameters, but the model grid doesnt have all 8 of the models that     */
/*  bound or bracket the input values.  The model grid is just incomplete   */
/*  for some low gravity cases.					            */
/*   No input or output required.                                           */
/*                                                                          */
/****************************************************************************/

void print_err2(void)
{

	printf("One (or more) of the gridpoints needed for the interpolation \n");
	printf("do not exist.  Retry the input, in case you made a typo. \n");
	exit(8);

} /* end of print_err2 subroutine */


/****************************************************************************/
/*                           print_instructions                             */
/*                                                                          */
/*  Subroutine writes out instructions to the screen if the user doesnt set */
/*  enough flags, the correct flags, or if just makes a boo-boo.            */
/*                                                                          */
/*   No input or output required.                                           */
/*                                                                          */
/****************************************************************************/

void print_instructions(char *name)
{

        printf("\n\t **** Usage::  %s [-w<filename>] [-tgmpoh] ****  \n", name);
	printf("\n");
        printf("*NOTE*: Must set -t -g and -p (or -m) for %s to run. \n\n", name);
	printf("*NOTE*: Must use -w flag FIRST!!!!  \n\n");
	printf(" -w<filename> = Filename of the output model file. \n");
        printf(" -t<number>   = Effective temperature <TTTT> \n");
        printf(" -g<number>   = Surface gravity <g.gg>  \n");
        printf(" -p<number>   = Metallicity <p.pp>. \n");
        printf("                 For positive metallicities or 0.0. \n");
        printf(" -m<number>   = Metallicity <m.mm>. \n");
        printf("                 For negative metallicites. _not_ 0.0. \n");
        printf(" -v<number>   = Microturbulence Parameter <v.vv> \n");
	printf(" -h           = This help menu. \n\n");
	printf(" Bounds : Temp = [3500,7750] (FOR MSPAWN72, sorry) \n");
	printf("          Grav = [0.0, 5.0] \n");
	printf("	  Met  = [-3.5, 0.5] \n\n");        
	exit(1);

} /* End of subroutine print_instructions */
