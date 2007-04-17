OPTIONAL_PROGS =

if BUILD_MPI
OPTIONAL_PROGS += \
	MPIArt \
	MPI_Angular_Predict \
	MPI_Angular_Predict_Continuous \
	MPI_MLalign2D \
	MPI_MLrefine3D \
	MPI_Projection_matching \
	MPI_WBP \
	MPI_Run
endif

if BUILD_TIFF
OPTIONAL_PROGS += Tiff2Raw
endif

if BUILD_QT
OPTIONAL_PROGS += \
	CTFViewer \
	Edit \
	Mark \
	Show \
	XMask
endif

SUBDIRS = \
	add_noise \
	adjust_volume_grey_levels \
	Align3D \
	align2D \
	Ang_distribution \
	Angular_distance \
	Angular_predict \
	Angular_predict_continuous \
	Angular_predict_tomography \
	ApplyGeo \
	AppendDocFile \
	Art \
	Assign_CTF \
	BatchSOM \
	Break_sym \
	Centilt \
	Compare_Selfiles \
	CorrectPhase \
	Correlation \
	CpSel \
	Create_phantom \
	CTF_profile \
	Data2Img \
	Data2Vol \
	Deform \
	Denoising \
	Denoising_micrograph \
	DocFile_Histogram \
	DownSample \
	EM2Angles \
	Enhance_PSD \
	Evaluate \
	FCMeans \
	FindCenter3D \
	FKCN \
	Flip \
	FourierFilter \
	FSOM \
	Headerinfo \
	Histogram \
	IDR_ART \
	Img2Data \
	Infogeo \
	KCMeans \
	KerDenSOM \
	Lattice_vectors \
	MakeSpectra \
	Mask \
	Markhan \
	Microscope \
	MLalign2D \
	MLalign2D_combine \
	MLrefine3D \
	Morphology \
	MvSel \
	Normalize \
	Operate \
	PCA \
	PDBPhantom \
	Phantom_transform \
	Project \
	Projection_Neighbourhood \
	Projection_matching \
	ProjectPCA \
	Pyramid \
	Range_adjust \
	Raw22spi \
	Recons_test \
	Reverse_endian \
	Resolution \
	RmSel \
	Rotate \
	Sammon \
	Scale \
	Scissor \
	Segment \
	Select_Images \
	SelStats \
	Separate_Objects \
	Shift \
	Skew \
	SOM \
	Sort_junk \
	Split_selfile \
	Spi22em \
	Spider22CCP4 \
	Spots2RealSpace2D \
	SSNR \
	Statis \
	Symmetrize \
	Threshold \
	ToPDB \
	ToPhantom \
	Unbend \
	Visualize_FFT \
	Vol2Data \
	VolSampler \
	Voxels22Blobs \
	Wbp \
	Window \
	Window_micrograph \
	$(OPTIONAL_PROGS)

