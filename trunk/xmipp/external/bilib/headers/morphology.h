/**@name Morphology */
//@{
/*--------------------------------------------------------------------------*/
/** Perform the grey-scale morphological operation known as bright top hat.
    VolumeSource is the grey-scale volume to process.
    VolumeDestination is the resulting grey-scale volume.
    Both VolumeSource and VolumeDestination have size (Nx, Ny, Nz).
    Kernel is the grey-scale structural element by which VolumeSource
    is processed. (Kx, Ky, Kz) is the size (bounding box) of the
    structural element. (Ox, Oy, Oz) is the origin of the structural element.

    Convention is the boundary convention applied to VolumeSource before
    processing.

    success: return(!ERROR); failure: return(ERROR) */
extern int  BrightTopHatFloat
    (
        float *VolumeSource,  /* data to process */
        float *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        float *Kernel,   /* structural element */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        long Ox,     /* kernel X origin */
        long Oy,     /* kernel Y origin */
        long Oz,     /* kernel Z origin */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform a morphological grey-scale bright top hat filter using a cuboid
    as structural element.
    Volume is the grey-scale volume to process in-place.
    Volume has size (Nx, Ny, Nz).
    (Kx, Ky, Kz) is the size of the structural element that is filled with
    the value 1.0F.
    The origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L,
    (Kz - 1L) / 2L).
    Convention is the boundary convention applied to Volume before processing.
    success: return(!ERROR); failure: return(ERROR) */
extern int  BrightTopHatFloatCuboid
    (
        float *Volume,   /* data to process in-place */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform the grey-scale morphological operation known as bright top hat.
    VolumeSource is the grey-scale volume to process.
    VolumeDestination is the resulting grey-scale volume.
    Both VolumeSource and VolumeDestination have size (Nx, Ny, Nz).
    Kernel is the grey-scale structural element by which VolumeSource is
    processed. (Kx, Ky, Kz) is the size (bounding box) of the structural
    element. (Ox, Oy, Oz) is the origin of the structural element.
    Convention is the boundary convention applied to VolumeSource before
    processing.
    success: return(!ERROR); failure: return(ERROR) */
extern int  BrightTopHatShort
    (
        short *VolumeSource,  /* data to process */
        short *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        short *Kernel,   /* structural element */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        long Ox,     /* kernel X origin */
        long Oy,     /* kernel Y origin */
        long Oz,     /* kernel Z origin */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform a morphological grey-scale bright top hat filter using a cuboid
    as structural element.
    Volume is the grey-scale volume to process in-place.
    Volume has size (Nx, Ny, Nz).
    (Kx, Ky, Kz) is the size of the structural element that is filled with
    the value 1.0F. The origin of the structural element is ((Kx - 1L) / 2L,
    (Ky - 1L) / 2L, (Kz - 1L) / 2L).
    Convention is the boundary convention applied to Volume before processing.

    success: return(!ERROR); failure: return(ERROR) */
extern int  BrightTopHatShortCuboid
    (
        short *Volume,   /* data to process in-place */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform the grey-scale morphological operation known as closing.
    VolumeSource is the grey-scale volume to close.
    VolumeDestination is the resulting grey-scale volume after closing.
    Both VolumeSource and VolumeDestination have size (Nx, Ny, Nz).
    Kernel is the grey-scale structural element by which VolumeSource is
    closed. (Kx, Ky, Kz) is the size (bounding box) of the structural element.
    (Ox, Oy, Oz) is the origin of the structural element.
    Convention is the boundary convention applied to VolumeSource before
    closing.

    success: return(!ERROR); failure: return(ERROR) */
extern int  ClosingFloat
    (
        float *VolumeSource,  /* data to process */
        float *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        float *Kernel,   /* structural element */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        long Ox,     /* kernel X origin */
        long Oy,     /* kernel Y origin */
        long Oz,     /* kernel Z origin */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform a morphological grey-scale closing filter using a cuboid as
    structural element.
    Volume is the grey-scale volume to close in-place.
    Volume has size (Nx, Ny, Nz).
    (Kx, Ky, Kz) is the size of the structural element that is filled with
    the value 1.0F.
    The origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L,
    (Kz - 1L) / 2L).
    Convention is the boundary convention applied to Volume before closing.

    success: return(!ERROR); failure: return(ERROR) */
extern int  ClosingFloatCuboid
    (
        float *Volume,   /* data to process in-place */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform the grey-scale morphological operation known as closing.
    VolumeSource is the grey-scale volume to close.
    VolumeDestination is the resulting grey-scale volume after closing.
    Both VolumeSource and VolumeDestination have size (Nx, Ny, Nz).
    Kernel is the grey-scale structural element by which VolumeSource is
    closed. (Kx, Ky, Kz) is the size (bounding box) of the structural element.
    (Ox, Oy, Oz) is the origin of the structural element.
    Convention is the boundary convention applied to VolumeSource before
    closing.

    success: return(!ERROR); failure: return(ERROR) */
extern int  ClosingShort
    (
        short *VolumeSource,  /* data to process */
        short *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        short *Kernel,   /* structural element */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        long Ox,     /* kernel X origin */
        long Oy,     /* kernel Y origin */
        long Oz,     /* kernel Z origin */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform a morphological grey-scale closing filter using a cuboid as
    structural element.
    Volume is the grey-scale volume to close in-place.
    Volume has size (Nx, Ny, Nz).
    (Kx, Ky, Kz) is the size of the structural element that is filled with
    the value 1.0F.
    The origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L,
    (Kz - 1L) / 2L).
    Convention is the boundary convention applied to Volume before closing.

    success: return(!ERROR); failure: return(ERROR) */
extern int  ClosingShortCuboid
    (
        short *Volume,   /* data to process in-place */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform the grey-scale morphological operation known as dark top hat.
    VolumeSource is the grey-scale volume to process.
    VolumeDestination is the resulting grey-scale volume.
    Both VolumeSource and VolumeDestination have size (Nx, Ny, Nz).
    Kernel is the grey-scale structural element by which VolumeSource
    is processed. (Kx, Ky, Kz) is the size (bounding box) of the structural
    element. (Ox, Oy, Oz) is the origin of the structural element.
    Convention is the boundary convention applied to VolumeSource before
    processing.

    success: return(!ERROR); failure: return(ERROR) */
extern int  DarkTopHatFloat
    (
        float *VolumeSource,  /* data to process */
        float *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        float *Kernel,   /* structural element */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        long Ox,     /* kernel X origin */
        long Oy,     /* kernel Y origin */
        long Oz,     /* kernel Z origin */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform a morphological grey-scale dark top hat filter using a cuboid
    as structural element.
    Volume is the grey-scale volume to process in-place.
    Volume has size (Nx, Ny, Nz).
    (Kx, Ky, Kz) is the size of the structural element that is filled with
    the value 1.0F. The origin of the structural element is ((Kx - 1L) / 2L,
    (Ky - 1L) / 2L, (Kz - 1L) / 2L).
    Convention is the boundary convention applied to Volume before processing
    success: return(!ERROR); failure: return(ERROR) */
extern int  DarkTopHatFloatCuboid
    (
        float *Volume,   /* data to process in-place */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform the grey-scale morphological operation known as dark top hat.
    VolumeSource is the grey-scale volume to process.
    VolumeDestination is the resulting grey-scale volume.
    Both VolumeSource and VolumeDestination have size (Nx, Ny, Nz).
    Kernel is the grey-scale structural element by which VolumeSource is
    processed. (Kx, Ky, Kz) is the size (bounding box) of the structural
    element. (Ox, Oy, Oz) is the origin of the structural element.
    Convention is the boundary convention applied to VolumeSource before
    processing.

    success: return(!ERROR); failure: return(ERROR) */
extern int  DarkTopHatShort
    (
        short *VolumeSource,  /* data to process */
        short *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        short *Kernel,   /* structural element */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        long Ox,     /* kernel X origin */
        long Oy,     /* kernel Y origin */
        long Oz,     /* kernel Z origin */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform a morphological grey-scale dark top hat filter using a cuboid
    as structural element.
    Volume is the grey-scale volume to process in-place.
    Volume has size (Nx, Ny, Nz).
    (Kx, Ky, Kz) is the size of the structural element that is filled
    with the value 1.0F.
    The origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L,
    (Kz - 1L) / 2L).
    Convention is the boundary convention applied to Volume before processing.

    success: return(!ERROR); failure: return(ERROR) */
extern int  DarkTopHatShortCuboid
    (
        short *Volume,   /* data to process in-place */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform the grey-scale morphological operation known as dilatation.
    VolumeSource is the grey-scale volume to dilate.
    VolumeDestination is the resulting grey-scale volume after dilatation.
    Both VolumeSource and VolumeDestination have size (Nx, Ny, Nz).
    Kernel is the grey-scale structural element by which VolumeSource is
    dilated. (Kx, Ky, Kz) is the size (bounding box) of the structural element.
    (Ox, Oy, Oz) is the origin of the structural element.
    Convention is the boundary convention applied to VolumeSource before
    dilatation.

    success: return(!ERROR); failure: return(ERROR) */
extern int  DilateFloat
    (
        float *VolumeSource,  /* data to process */
        float *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        float *Kernel,   /* structural element */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        long Ox,     /* kernel X origin */
        long Oy,     /* kernel Y origin */
        long Oz,     /* kernel Z origin */
        enum TBoundaryConvention
        Convention   /* boundary convention */
    );

/*--------------------------------------------------------------------------*/
/** Perform the grey-scale morphological operation known as dilatation.
    VolumeSource is the grey-scale volume to dilate.
    VolumeDestination is the resulting grey-scale volume after dilatation.
    Both VolumeSource and VolumeDestination have size (Nx, Ny, Nz).
    Kernel is the grey-scale structural element by which VolumeSource is
    dilated. (Kx, Ky, Kz) is the size (bounding box) of the structural element.
    (Ox, Oy, Oz) is the origin of the structural element.
    Convention is the boundary convention applied to VolumeSource before
    dilatation.

    success: return(!ERROR); failure: return(ERROR) */
extern int  DilateShort
    (
        short *VolumeSource,  /* data to process */
        short *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        short *Kernel,   /* structural element */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        long Ox,     /* kernel X origin */
        long Oy,     /* kernel Y origin */
        long Oz,     /* kernel Z origin */
        enum TBoundaryConvention
        Convention   /* boundary convention */
    );

/*--------------------------------------------------------------------------*/
/** Perform the grey-scale morphological operation known as erosion.
    VolumeSource is the grey-scale volume to erode.
    VolumeDestination is the resulting grey-scale volume after erosion.
    Both VolumeSource and VolumeDestination have size (Nx, Ny, Nz).
    Kernel is the grey-scale structural element by which VolumeSource is
    eroded. (Kx, Ky, Kz) is the size (bounding box) of the structural element.
    (Ox, Oy, Oz) is the origin of the structural element.
    Convention is the boundary convention applied to VolumeSource before
    erosion.
    success: return(!ERROR); failure: return(ERROR) */
extern int  ErodeFloat
    (
        float *VolumeSource,  /* data to process */
        float *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        float *Kernel,   /* structural element */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        long Ox,     /* kernel X origin */
        long Oy,     /* kernel Y origin */
        long Oz,     /* kernel Z origin */
        enum TBoundaryConvention
        Convention   /* boundary convention */
    );

/*--------------------------------------------------------------------------*/
/** Perform the grey-scale morphological operation known as erosion.
    VolumeSource is the grey-scale volume to erode.
    VolumeDestination is the resulting grey-scale volume after erosion.
    Both VolumeSource and VolumeDestination have size (Nx, Ny, Nz).
    Kernel is the grey-scale structural element by which VolumeSource is eroded.
    (Kx, Ky, Kz) is the size (bounding box) of the structural element.
    (Ox, Oy, Oz) is the origin of the structural element.
    Convention is the boundary convention applied to VolumeSource before erosion.

    success: return(!ERROR); failure: return(ERROR) */
extern int  ErodeShort
    (
        short *VolumeSource,  /* data to process */
        short *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        short *Kernel,   /* structural element */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        long Ox,     /* kernel X origin */
        long Oy,     /* kernel Y origin */
        long Oz,     /* kernel Z origin */
        enum TBoundaryConvention
        Convention   /* boundary convention */
    );

/*--------------------------------------------------------------------------*/
/** Perform a morphological grey-scale max filter using a cuboid as
    structural element.
    Volume is the grey-scale volume to dilate in-place.
    Volume has size (Nx, Ny, Nz).
    (Kx, Ky, Kz) is the size of the structural element that is filled with
    the value 1.0F. The origin of the structural element is ((Kx - 1L) / 2L,
    (Ky - 1L) / 2L, (Kz - 1L) / 2L).
    Convention is the boundary convention applied to Volume before dilatation.

    success: return(!ERROR); failure: return(ERROR) */
extern int  MaxFilterFloatCuboid
    (
        float *Volume,   /* data to process in-place */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform a morphological grey-scale max filter using a cuboid as
    structural element.
    Volume is the grey-scale volume to dilate in-place.
    Volume has size (Nx, Ny, Nz).
    (Kx, Ky, Kz) is the size of the structural element that is filled
    with the value 1.0F.
    The origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L,
    (Kz - 1L) / 2L).
    Convention is the boundary convention applied to Volume before dilatation.

    success: return(!ERROR); failure: return(ERROR) */
extern int  MaxFilterShortCuboid
    (
        short *Volume,   /* data to process in-place */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform a morphological grey-scale min filter using a cuboid as
    structural element.
    Volume is the grey-scale volume to erode in-place.
    Volume has size (Nx, Ny, Nz).
    (Kx, Ky, Kz) is the size of the structural element that is filled with
    the value 1.0F.
    The origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L,
    (Kz - 1L) / 2L).
    Convention is the boundary convention applied to Volume before erosion.

    success: return(!ERROR); failure: return(ERROR) */
extern int  MinFilterFloatCuboid
    (
        float *Volume,   /* data to process in-place */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform a morphological grey-scale min filter using a cuboid as
    structural element.
    Volume is the grey-scale volume to erode in-place.
    Volume has size (Nx, Ny, Nz).
    (Kx, Ky, Kz) is the size of the structural element that is filled with
    the value 1.0F.
    The origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L,
    (Kz - 1L) / 2L).
    Convention is the boundary convention applied to Volume before erosion.

    success: return(!ERROR); failure: return(ERROR) */
extern int  MinFilterShortCuboid
    (
        short *Volume,   /* data to process in-place */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform the grey-scale morphological operation known as morphological
    gradient.
    VolumeSource is the grey-scale volume to process.
    VolumeDestination is the resulting grey-scale volume after the gradient
    operation.
    Both VolumeSource and VolumeDestination have size (Nx, Ny, Nz).
    Kernel is the grey-scale structural element by which VolumeSource is
    processed. (Kx, Ky, Kz) is the size (bounding box) of the structural
    element. (Ox, Oy, Oz) is the origin of the structural element.
    Convention is the boundary convention applied to VolumeSource before
    the gradient operation.

    success: return(!ERROR); failure: return(ERROR) */
extern int  MorphologicalGradientFloat
    (
        float *VolumeSource,  /* data to process */
        float *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        float *Kernel,   /* structural element */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        long Ox,     /* kernel X origin */
        long Oy,     /* kernel Y origin */
        long Oz,     /* kernel Z origin */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform a morphological grey-scale gradient filter using a cuboid as
    structural element.
    Volume is the grey-scale volume to process in-place.
    Volume has size (Nx, Ny, Nz).
    (Kx, Ky, Kz) is the size of the structural element that is filled with
    the value 1.0F.
    The origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L,
    (Kz - 1L) / 2L).
    Convention is the boundary convention applied to Volume before processing.

    success: return(!ERROR); failure: return(ERROR) */
extern int  MorphologicalGradientFloatCuboid
    (
        float *Volume,   /* data to process in-place */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform the grey-scale morphological operation known as morphological
    gradient.
    VolumeSource is the grey-scale volume to process.
    VolumeDestination is the resulting grey-scale volume after the gradient
    operation.
    Both VolumeSource and VolumeDestination have size (Nx, Ny, Nz).
    Kernel is the grey-scale structural element by which VolumeSource is
    processed. (Kx, Ky, Kz) is the size (bounding box) of the structural
    element. (Ox, Oy, Oz) is the origin of the structural element.
    Convention is the boundary convention applied to VolumeSource before
    the gradient operation.

    success: return(!ERROR); failure: return(ERROR) */
extern int  MorphologicalGradientShort
    (
        short *VolumeSource,  /* data to process */
        short *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        short *Kernel,   /* structural element */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        long Ox,     /* kernel X origin */
        long Oy,     /* kernel Y origin */
        long Oz,     /* kernel Z origin */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform a morphological grey-scale gradient filter using a cuboid as
    structural element.
    Volume is the grey-scale volume to process in-place.
    Volume has size (Nx, Ny, Nz).
    (Kx, Ky, Kz) is the size of the structural element that is filled with
    the value 1.0F.
    The origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L,
    (Kz - 1L) / 2L).
    Convention is the boundary convention applied to Volume before processing.

    success: return(!ERROR); failure: return(ERROR) */
extern int  MorphologicalGradientShortCuboid
    (
        short *Volume,   /* data to process in-place */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform the grey-scale morphological operation known as opening.
    VolumeSource is the grey-scale volume to open.
    VolumeDestination is the resulting grey-scale volume after opening.
    Both VolumeSource and VolumeDestination have size (Nx, Ny, Nz).
    Kernel is the grey-scale structural element by which VolumeSource is
    opened. (Kx, Ky, Kz) is the size (bounding box) of the structural element.
    (Ox, Oy, Oz) is the origin of the structural element.
    Convention is the boundary convention applied to VolumeSource before opening.

    success: return(!ERROR); failure: return(ERROR) */
extern int  OpeningFloat
    (
        float *VolumeSource,  /* data to process */
        float *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        float *Kernel,   /* structural element */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        long Ox,     /* kernel X origin */
        long Oy,     /* kernel Y origin */
        long Oz,     /* kernel Z origin */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform a morphological grey-scale opening filter using a cuboid as
    structural element.
    Volume is the grey-scale volume to open in-place.
    Volume has size (Nx, Ny, Nz).
    (Kx, Ky, Kz) is the size of the structural element that is filled with
    the value 1.0F.
    The origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L,
    (Kz - 1L) / 2L).
    Convention is the boundary convention applied to Volume before opening.

    success: return(!ERROR); failure: return(ERROR) */
extern int  OpeningFloatCuboid
    (
        float *Volume,   /* data to process in-place */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform the grey-scale morphological operation known as opening.
    VolumeSource is the grey-scale volume to open.
    VolumeDestination is the resulting grey-scale volume after opening.
    Both VolumeSource and VolumeDestination have size (Nx, Ny, Nz).
    Kernel is the grey-scale structural element by which VolumeSource is
    opened. (Kx, Ky, Kz) is the size (bounding box) of the structural element.
    (Ox, Oy, Oz) is the origin of the structural element.
    Convention is the boundary convention applied to VolumeSource before
    opening.

    success: return(!ERROR); failure: return(ERROR) */
extern int  OpeningShort
    (
        short *VolumeSource,  /* data to process */
        short *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        short *Kernel,   /* structural element */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        long Ox,     /* kernel X origin */
        long Oy,     /* kernel Y origin */
        long Oz,     /* kernel Z origin */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Perform a morphological grey-scale opening filter using a cuboid as
    structural element.
    Volume is the grey-scale volume to open in-place.
    Volume has size (Nx, Ny, Nz).
    (Kx, Ky, Kz) is the size of the structural element that is filled with
    the value 1.0F
    The origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L,
    (Kz - 1L) / 2L)
    Convention is the boundary convention applied to Volume before opening.

    success: return(!ERROR); failure: return(ERROR) */
extern int  OpeningShortCuboid
    (
        short *Volume,   /* data to process in-place */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long Kx,     /* width of the kernel */
        long Ky,     /* height of the kernel */
        long Kz,     /* depth of the kernel */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        int  *Status    /* error management */
    );
//@}
