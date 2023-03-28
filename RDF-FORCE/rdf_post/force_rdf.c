
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define GAS_CONST 8.3144626181
#define KCAL2JOULE 4184.0
#define MAXLINE 1024
#define BinNum 4000
#define Delta 0.01

int parse_dump_timestep(FILE *stream, long *Timestep, long *NumberOfParticles, 
                        double *Box, double *Rxx[], double *Ryy[], 
                        double *Rzz[], double *Fxx[], double *Fyy[], double *Fzz[]) 
{
    /**
     * @brief Parse a single timestep from a LAMMPS dump file created with 
                    "custom ${Nrun} dump id x y z fx fy fz"
     * @param stream File pointer of the dump file with its offset at the 
                        beginning of the timestep, offset is always moved to 
                        the beginning of the next timestep
     * @param NumberOfParticles Return addres for the number of particles in the timestep
     * @param Box Return address for the box size, box is assumed to be cubic 
                    otherwise the function exits and returns 1
     * @param n Size of the arrays Rxx, Ryy, Rzz, Fxx, Fyy, Fzz, 
                    if the number of particles in the timestep is larger than 
                    n, the function exits and returns 2
     * @param Rxx Return array for the x coordinates of the particles
     * @param Ryy Return array for the y coordinates of the particles
     * @param Rzz Return array for the z coordinates of the particles
     * @param Fxx Return array for the x components of the forces on the particles
     * @param Fyy Return array for the y components of the forces on the particles
     * @param Fzz Return array for the z components of the forces on the particles
     * @return 0 if the timestep was parsed successfully, 1 if an error occured
     */
    char buffer[MAXLINE];
    double boxLeft[3], boxRight[3];
    long i, ID;

    fgets(buffer, MAXLINE, stream);
    if (feof(stream)) {
        return 1;
    } else if (strstr(buffer, "ITEM: TIMESTEP") == NULL) {
        fprintf(stderr, "Expected next line to be 'ITEM: TIMESTEP', but instead \
                    is:\n%s.\n", buffer);
        exit(1);
    } else { 
        fgets(buffer, MAXLINE, stream);
        *Timestep = atoi(buffer);
    }

    fgets(buffer, MAXLINE, stream);
    if (strstr(buffer, "ITEM: NUMBER OF ATOMS") == NULL) {
        fprintf(stderr, "Expected next line to be 'ITEM: NUMBER OF ATOMS', but \
                    instead is:\n%s.\n", buffer);
        exit(1);
    } else {
        fgets(buffer, MAXLINE, stream);
        *NumberOfParticles = atoi(buffer);
    }

    fgets(buffer, MAXLINE, stream);
    if (strstr(buffer, "ITEM: BOX BOUNDS") == NULL) {
        fprintf(stderr, "Expected next line to be 'ITEM: BOX BOUNDS', \
                        but instead is:\n%s.\n", buffer);
        exit(1);
    } else {
        fgets(buffer, MAXLINE, stream);
        sscanf(buffer, "%lf %lf", &boxLeft[0], &boxRight[0]);
        fgets(buffer, MAXLINE, stream);
        sscanf(buffer, "%lf %lf", &boxLeft[1], &boxRight[1]);
        fgets(buffer, MAXLINE, stream);
        sscanf(buffer, "%lf %lf", &boxLeft[2], &boxRight[2]);
        // TODO: this is comparing floats, maybe change to checking "if close"
        if (boxLeft[0] != boxLeft[1] || boxLeft[0] != boxLeft[2] || 
            boxRight[0] != boxRight[1] || boxRight[0] != boxRight[2]) {
            fprintf(stderr, "Box is not a cube.\n");
            exit(1);
        }
    }

    *Box = boxRight[0] - boxLeft[0];

    fgets(buffer, MAXLINE, stream);
    if (strstr(buffer, "ITEM: ATOMS") == NULL) {
        fprintf(stderr, "Expected next line to be 'ITEM: ATOMS', but instead \
                    is:\n%s.\n", buffer);
        exit(1);
    } else{
        // free memory if it was allocated before
        if (*Rxx != NULL || *Ryy != NULL || *Rzz != NULL || *Fxx != NULL || 
            *Fyy != NULL || *Fzz != NULL) {
            free(*Rxx), free(*Ryy), free(*Rzz);
            free(*Fxx), free(*Fyy), free(*Fzz);
        } 
        
        *Rxx = (double *) malloc(*NumberOfParticles * sizeof(double));
        *Ryy = (double *) malloc(*NumberOfParticles * sizeof(double));
        *Rzz = (double *) malloc(*NumberOfParticles * sizeof(double));
        *Fxx = (double *) malloc(*NumberOfParticles * sizeof(double));
        *Fyy = (double *) malloc(*NumberOfParticles * sizeof(double));
        *Fzz = (double *) malloc(*NumberOfParticles * sizeof(double));

        if (*Rxx == NULL || *Ryy == NULL || *Rzz == NULL || *Fxx == NULL || 
            *Fyy == NULL || *Fzz == NULL) {
            free(*Rxx), free(*Ryy), free(*Rzz);
            free(*Fxx), free(*Fyy), free(*Fzz);
            fprintf(stderr, "Error allocating memory for the arrays to hold atom \
                posistion and forces.\n");
            exit(1);
            }

        for (i = 0; i < *NumberOfParticles; i++) {
            fgets(buffer, MAXLINE, stream);
            sscanf(buffer, "%ld %lf %lf %lf %lf %lf %lf", &ID, &((*Rxx)[i]), 
                    &((*Ryy)[i]), &((*Rzz)[i]), &((*Fxx)[i]), &((*Fyy)[i]), &((*Fzz)[i]));
        }
    }
    return 0;
}


double get_temp(long timestep, char *fname) 
{
    /**
     * @brief Get the temperature from a LAMMPS fix ave/time file
     * @param timestep Timestep to get the temperature for
     * @param fname Name of the log file
     * @return Temperature at the given timestep
     */
    FILE *fp;
    char buffer[MAXLINE];
    double Temp;
    long Time;

    fp = fopen(fname, "r");
    if (fp == NULL) {
      fprintf(stderr, "Could not open logfile %s.\n", fname);
      exit(1);
    }

    while (fgets(buffer, MAXLINE, fp) != NULL) {
        if (buffer[0] != '#') {
            sscanf(buffer, "%ld %lf", &Time, &Temp);
            if (Time == timestep) {
                fclose(fp);
                return Temp;
            }
        }
    }   
    fclose(fp);
    fprintf(stderr, "Could not find temperature for timestep %ld.\n", timestep);
    exit(1);
}


void sample_rdf(const long NumberOfParticles, const double Box, double g[], 
                const double Rxx[], const double Ryy[], const double Rzz[], 
                const double Fxx[], const double Fyy[], const double Fzz[]) 
{
    /**
     * @brief Sample the radial distribution function of the system in array g 

     * @param NumberOfParticles Number of particles in the system
     * @param Box Size of the box
     * @param g Array to store the sampled rdf, array values are zeroed before sampling.
     * @param Rxx Array of x coordinates of the particles
     * @param Ryy Array of y coordinates of the particles
     * @param Rzz Array of z coordinates of the particles
     * @param Fxx Array of x components of the forces on the particles
     * @param Fyy Array of y components of the forces on the particles
     * @param Fzz Array of z components of the forces on the particles
     */
    long i, j, k;
    double Dx, Dy, Dz, Rij, Ff, sum=0.0;
    const double HalfBox=Box/2.0;

    for (i = 0; i < BinNum; i++) {
        g[i] = 0.0;
    }

    for(i=0;i<NumberOfParticles;i++) {
        for(j=i+1;j<NumberOfParticles;j++) { 
            // Calculate the distance between the two particles, R, accounting
            // for periodic boundary conditions
            Dx = Rxx[i] - Rxx[j];
            Dy = Ryy[i] - Ryy[j];
            Dz = Rzz[i] - Rzz[j];
 
            if (Dx > HalfBox) {
               Dx = Dx - Box;
            } else if (Dx < -HalfBox) {
               Dx = Dx + Box;
            }
 
            if (Dy > HalfBox) {
               Dy = Dy - Box;
            } else if (Dy< -HalfBox) {
               Dy = Dy + Box;
            }
 
            if (Dz > HalfBox) {
               Dz = Dz - Box;
            } else if (Dz < -HalfBox) {
               Dz = Dz + Box;
            }

            Rij = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);

            // RDF is only valid till half the box size
            if (Rij < HalfBox) {
                Ff = ((Fxx[i]-Fxx[j])*Dx + (Fyy[i]-Fyy[j])*Dy + (Fzz[i]-Fzz[j])*Dz) / (Rij*Rij*Rij);
                
                // R = (k+0.5)*Delta
                // smallest k for R > Rij
                k = ceil(Rij / Delta - 0.5);
                
                if (k < BinNum) { 
                    g[k] += Ff; 
                    // only add Ff to the smallest k for now, should be added to all k's where R>Rij but this is done later all at once
                }
            }
        }
    }

    // Cumulative sum because of the heaviside function
    for (i = 0; i < BinNum; i++) {
        sum += g[i];
        g[i] = sum;
    }
}


int main(int argc, char *argv[])
{
    FILE *fp;
    int eof;
    long NumberOfParticles, i, k, timestep, nTimesteps, minTimestep=0;
    double Temp, c, Box, R, *Rxx, *Ryy, *Rzz, *Fxx, *Fyy, *Fzz, *g, *gavg;
    Rxx = (double *) NULL, Ryy = (double *) NULL, Rzz = (double *) NULL;
    Fxx = (double *) NULL, Fyy = (double *) NULL, Fzz = (double *) NULL;

    if (argc < 4) {
        fprintf(stderr, "Usage: %s <dump file> <Temp log file> <output file> \
                        (optional)<minimum timestep>\n", argv[0]);
        exit(1);
    }

    if (argc == 5) {
        minTimestep = atol(argv[4]);
    }

    fp = fopen(argv[1], "r");
    if (fp == NULL) {
        fprintf(stderr, "Could not open dump file %s'\n", argv[1]);
        exit(1);
    }

    g = (double *) malloc(BinNum * sizeof(double));
    gavg = (double *) calloc(BinNum, sizeof(double));    

    nTimesteps = 0;
    while (feof(fp) == 0) { 
        eof = parse_dump_timestep(fp, &timestep, &NumberOfParticles, &Box, 
                                    &Rxx, &Ryy, &Rzz, &Fxx, &Fyy, &Fzz);

        if ((eof == 0) && (timestep >= minTimestep)) {
            sample_rdf(NumberOfParticles, Box, g, Rxx, Ryy, Rzz, Fxx, Fyy, Fzz);
            Temp = get_temp(timestep, argv[2]);

            c = KCAL2JOULE/GAS_CONST * Box*Box*Box / 
            (4 * M_PI * Temp * NumberOfParticles*NumberOfParticles); 

            for (k = 0; k < BinNum; k++) {
                gavg[k] += g[k] * c;
            }
            nTimesteps++;
        }
    }

    if (nTimesteps == 0) {
        fprintf(stderr, "No timesteps larger than the minimum timestep found in dump file.\n");
        exit(1);
    } else {
        for (k = 0; k < BinNum; k++) {
            gavg[k] /= nTimesteps; // Average the rdf
        }
    } 
    printf("nTimesteps = %ld\n", nTimesteps);
    fclose(fp);

    // Write the rdf to a file
    fp = fopen(argv[3], "w");
    if (fp == NULL) {
        fprintf(stderr, "Could not open output file %s'\n", argv[3]);
        exit(1);
    }
    for (i = 0; i < BinNum; i++) {
        R = (i + 0.5) * Delta;
        if (R < Box/2.0){
            fprintf(fp, "%lf %lf\n", R, gavg[i]);
        }
    }
    fclose(fp);

    free(g), free(gavg);
    free(Rxx), free(Ryy), free(Rzz);
    free(Fxx), free(Fyy), free(Fzz);
    return 0;
}