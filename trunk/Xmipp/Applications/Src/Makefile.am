OPTIONAL_PROGS =

if BUILD_MPI
OPTIONAL_PROGS += MPIArt MPI_WBP 
endif

if BUILD_TIFF
OPTIONAL_PROGS += Tiff2Raw 
endif

if BUILD_QT
if BUILD_OGL
OPTIONAL_PROGS += Phan3D
endif
OPTIONAL_PROGS += Mark Show XMask 
endif

SUBDIRS = 2DInterpolation \
	3DInterpolation \
	Adapt_for_Spider \
	Add_noise \
	Adjust_CTF \
	Adjust_surface \
	Align3D \
	Align2D \
	Ang_distribution \
	Angular_distance \
	Angular_predict \
	Angular_predict_continuous \
	Angular_refinement \
	ApplyGeo \
	Are_different \
	Art \
	Assign_angles \
	Assign_CTF \
	Assign_multiple_CTFs \
	BatchSOM \
	CCP4Spider \
	Centilt \
	CodeStatis \
	Combine_stats \
	Compare_Selfiles \
	CorrectPhase \
	Correlation \
	CpSel \
	Create_phantom \
	CTF_profile \
	Data2Img \
	Data2Vol \
	deconvoluteCTF \
	Deform \
	Denoising \
	Denoising_micrograph \
	DocFile_Histogram \
	DownSample \
	Draw_surface \
	DWT \
	Euler \
	Evaluate \
	Evaluate_FSCs \
	FCMeans \
	FFT \
	FindCenter3D \
	FKCN \
	Flip \
	FourierFilter \
	FSOM \
	Headerinfo \
	Histogram \
	IDR_ART \
	IDWT \
	IFFT \
	Img2Data \
	Infogeo \
	KCMeans \
	KerDenSOM \
	Lattice_vectors \
	Mask \
	Microscope \
	MLalign2D \
	MLalign2D_combine \
	Morphology \
	MvSel \
	Normalize \
	Operate \
	PCA \
	Phantom_transform \
	Project \
	Projection_Neighbourhood \
	Projection_matching \
	ProjectPCA \
	Pyramid \
	Radon_Transform \
	RandomMapping \
	Random_phantom \
	Range_adjust \
	Raw22spi \
	Recons_test \
	Reverse_endian \
	Resolution \
	RmSel \
	Rotate \
	Sammon \
	Sample_size \
	Scale \
	Scissor \
	Segment \
	Select_Images \
	SelStats \
	Separate_Objects \
	Shift \
	Skew \
	SOM \
	Split_selfile \
	SpAR \
	SpARMA \
	Spots2RealSpace2D \
	SSNR \
	Statis \
	Surface \
	Symmetrize \
	Taxa_tilt \
	Threshold \
	ToPDB \
	ToPhantom \
	Visualize_FFT \
	Vol2Data \
	VolSampler \
	Voxels22Blobs \
	Wbp \
	Window \
	Window_micrograph \
	$(OPTIONAL_PROGS)

