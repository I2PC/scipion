/*
        rwMRC.h
        Header file for reading and writing MRC files
        Format: 3D crystallographic image file format for the MRC package
        Author: Bernard Heymann
        Created: 19990321       Modified: 20030723
*/

#ifndef RWMRC_H
#define RWMRC_H

#define MRCSIZE    1024 // Minimum size of the MRC header (when nsymbt = 0)

struct MRCheadold {          // file header for MRC data
        int nx;              //  0   0       image size
        int ny;              //  1   4
        int nz;              //  2   8
        int mode;            //  3           0=uchar,1=short,2=float
        int nxStart;         //  4           unit cell offset
        int nyStart;         //  5
        int nzStart;         //  6
        int mx;              //  7           unit cell size in voxels
        int my;              //  8
        int mz;              //  9
        float a;             // 10   40      cell dimensions in A
        float b;             // 11
        float c;             // 12
        float alpha;         // 13           cell angles in degrees
        float beta;          // 14
        float gamma;         // 15
        int mapc;            // 16           column axis
        int mapr;            // 17           row axis
        int maps;            // 18           section axis
        float amin;          // 19           minimum density value
        float amax;          // 20   80      maximum density value
        float amean;         // 21           average density value
        int ispg;            // 22           space group number
        int nsymbt;          // 23           bytes used for sym. ops. table
        float extra[29];     // 24           user-defined info
        float xOrigin;       // 53           phase origin in pixels
        float yOrigin;       // 54
        int nlabl;           // 55           number of labels used
        char labels[10][80]; // 56-255       10 80-character labels
} ;

struct MRChead {             // file header for MRC data
        int nx;              //  0   0       image size
        int ny;              //  1   4
        int nz;              //  2   8
        int mode;            //  3           0=char,1=short,2=float
        int nxStart;         //  4           unit cell offset
        int nyStart;         //  5
        int nzStart;         //  6
        int mx;              //  7           unit cell size in voxels
        int my;              //  8
        int mz;              //  9
        float a;             // 10   40      cell dimensions in A
        float b;             // 11
        float c;             // 12
        float alpha;         // 13           cell angles in degrees
        float beta;          // 14
        float gamma;         // 15
        int mapc;            // 16           column axis
        int mapr;            // 17           row axis
        int maps;            // 18           section axis
        float amin;          // 19           minimum density value
        float amax;          // 20   80      maximum density value
        float amean;         // 21           average density value
        int ispg;            // 22           space group number
        int nsymbt;          // 23           bytes used for sym. ops. table
        float extra[25];     // 24           user-defined info
        float xOrigin;       // 49           phase origin in pixels
        float yOrigin;       // 50
        float zOrigin;       // 51
        char map[4];         // 52       identifier for map file ("MAP ")
        char machst[4];      // 53           machine stamp
        float arms;          // 54       RMS deviation
        int nlabl;           // 55           number of labels used
        char labels[800];    // 56-255       10 80-character labels
} ;

// I/O prototypes
int readMRC()
{
    FILE        *fimg;
    if ( ( fimg = fopen(filename.c_str(), "r") ) == NULL ) return(-1);

        MRChead*        header = (MRChead *) balloc(sizeof(MRChead));
        if ( fread( header, MRCSIZE, 1, fimg ) < 1 ) return(-2);

    // Determine byte order and swap bytes if from little-endian machine
    char*       b = (char *) header;
    int         i, swap = 0;
    if ( ( abs( header->mode ) > SWAPTRIG ) || ( abs(header->nz) > SWAPTRIG ) ) {
#ifdef DEBUG        
        fprintf(stderr, "Warning: Swapping header byte order for 4-byte types\n");
#endif
        swap = 1;
        int     extent = MRCSIZE - 800; // exclude labels from swapping
        for ( i=0; i<extent; i+=4 ) swapbytes(b+i, 4);
    }

    // Convert VAX floating point types if necessary
    if ( header->amin > header->amax )
        REPORT_ERROR(1,"readMRC: VAX floating point conversion unsupported");

    // Map the parameters
    XSIZE(data) = header->nx;
    YSIZE(data) = header->ny;
    ZSIZE(data) = header->nz;
    switch ( header->mode%5 ) {
        case 0: datatype = SChar; break;
        case 1: datatype = Short; break;
        case 2: datatype = Float; break;
        case 3: datatype = ComplexShort; break;
        case 4: datatype = ComplexFloat; break;
        default: datatype = UChar; break;
    }
    offset = MRCSIZE + header->nsymbt;
    if ( header->mode%5 > 2 && header->mode%5 < 5 ) {
        transform = CentHerm;
        fseek(fimg, 0, SEEK_END);
        if ( ftell(fimg) > offset + 0.8*XSIZE(data)*YSIZE(data)*ZSIZE(data)*gettypesize(datatype) )
            XSIZE(data) = 2*(XSIZE(data) - 1);
        if ( header->mx%2 == 1 ) XSIZE(data) += 1;     // Quick fix for odd x-size maps
    }

    min = header->amin;
    max = header->amax;
    avg = header->amean;
    std = header->arms;
    ua = header->a;
    ub = header->b;
    uc = header->c;
    alf = header->alpha;
    bet = header->beta;
    gam = header->gamma;
    if ( header->mx ) ux = header->a/header->mx;
    if ( header->my ) uy = header->b/header->my;
    if ( header->mz ) uz = header->c/header->mz;
    spacegroup = header->ispg;
    
    // Allocating the single sub-image and setting its origin
    image = new SubImage [1];
    if (image == NULL)
            REPORT_ERROR(1001, "Allocate: No space left for subimage structure");
    //image->xoff = -header->nxStart;	// Old header
    //image->yoff = -header->nyStart;
    //image->zoff = -header->nzStart;
    image->xoff = -header->xOrigin/ux;	// New header
    image->yoff = -header->yOrigin/uy;
    image->zoff = -header->zOrigin/uz;
    image->rot = image->tilt = image->psi = image->flip = 0.;
    image->weight = 1.;

    bfree(header, sizeof(MRChead));

    readData(fimg, -1, swap, 0);

    fclose(fimg);

    return(0);
}
/*
int     writeMRC(Bimage* p)
{
        img_RGB2gray(p);

    switch ( p->datatype ) {
        case UChar: img_to_signed_byte(p); break;
//      case SChar: img_to_byte(p); break;
        case UShort: img_to_short(p); break;
        case Int: img_to_float(p); break;
        case ComplexInt: img_to_complex_float(p); break;
        default: break;
    }

        if ( p->transform != NoTransform )
                img_convert_fourier(p, CentHerm);

        MRChead*        header = (MRChead *) balloc(sizeof(MRChead));

        // Map the parameters
        strncpy(header->map, "MAP ", 4);
        set_CCP4_machine_stamp(header->machst);
        header->nx = p->x;
        header->ny = p->y;
        header->nz = p->z;
        if ( p->transform == CentHerm ) header->nx = p->x/2 + 1;        // If a transform, physical storag
e is nx/2 + 1
        switch ( p->datatype ) {
                case UChar:
                case SChar: header->mode = 0; break;
                case Short: header->mode = 1; break;
                case Float: header->mode = 2; break;
                case ComplexShort: header->mode = 3; break;
                case ComplexFloat: header->mode = 4; break;
                default: header->mode = 0; break;
        }
        header->nxStart = (int) (-p->image->ox - 0.5);
        header->nyStart = (int) (-p->image->oy - 0.5);
        header->nzStart = (int) (-p->image->oz - 0.5);
        header->mx = (int) (p->ua/p->ux + 0.5);
        header->my = (int) (p->ub/p->uy + 0.5);
        header->mz = (int) (p->uc/p->uz + 0.5);
        header->mapc = 1;
        header->mapr = 2;
        header->maps = 3;
        header->amin = p->min;
        header->amax = p->max;
        header->amean = p->avg;
        header->arms = p->std;
        header->a = p->ua;
        header->b = p->ub;
        header->c = p->uc;
//      header->xOrigin = p->image->ox;
//      header->yOrigin = p->image->oy;
//      header->zOrigin = p->image->oz;
        header->xOrigin = -p->image->ox*p->ux;
        header->yOrigin = -p->image->oy*p->uy;
        header->zOrigin = -p->image->oz*p->uz;

        // This is a band-aid to overcome the limitations of the image format
        if ( fabs(p->ua - p->ux*header->mx) > 0.001 || fabs(p->ub - p->uy*header->my) > 0.001 ||
                        fabs(p->uc - p->uz*header->mz) > 0.001 ) {
                header->a = p->ux*header->mx;
                header->b = p->uy*header->my;
                header->c = p->uz*header->mz;
                if ( verbose )
                        fprintf(stderr, "Warning: Resetting the unit cell to: %g %g %g A\n", header->a, he
ader->b, header->c );
        }

       header->alpha = p->alf*180/M_PI;
        header->beta = p->bet*180/M_PI;
        header->gamma = p->gam*180/M_PI;
        header->ispg = p->spacegroup;

        int                     nsym = 0;
        Bstring         temp;
        char*           symop = NULL;

#ifndef NOSYMOP
        if ( p->spacegroup > 0 ) symop = read_symop(temp, p->spacegroup, &nsym);
#endif

        header->nsymbt = nsym*80;

        if ( verbose & VERB_DEBUG )
                printf("DEBUG rwMRC: Nsymbt = %d\n", header->nsymbt);

        long                    i;

        header->nlabl = 10;
        strncpy(header->labels, p->label.c_str(), 799);

        if ( verbose & VERB_DEBUG ) {
                printf("DEBUG rwMRC: Nlabels = %d\n", header->nlabl);
                for ( i=0; i<header->nlabl; i++ )
                        if ( header->labels[i*80] ) printf("%-80s\n", &header->labels[i*80]);
        }

        p->offset = MRCSIZE + header->nsymbt;

       long                    datatypesize = gettypesize(p->datatype);
        unsigned long   datasize = (unsigned long) header->nx*header->ny*header->nz*datatypesize;

        if ( verbose & VERB_DEBUG )
                printf("DEBUG rwMRC: Offset = %ld,  Typesize = %ld,  Datasize = %ld\n",
                                p->offset, datatypesize, datasize);

    FILE        *fimg;
    if ( ( fimg = fopen(p->filename.c_str(), "w") ) == NULL ) return(-1);

        fwrite( header, MRCSIZE, 1, fimg );
        if ( header->nsymbt ) fwrite( symop, header->nsymbt, 1, fimg );
        if ( p->dataflag ) fwrite( p->data, datasize, 1, fimg );

        fclose(fimg);

        if ( symop ) bfree(symop, header->nsymbt*sizeof(char));
        bfree(header, sizeof(MRChead));

        return(0);
}
*/
#endif
