################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../libraries/parallel/mpi_angular_class_average.cpp \
../libraries/parallel/mpi_angular_gcar_commonlines.cpp \
../libraries/parallel/mpi_angular_projection_matching.cpp \
../libraries/parallel/mpi_classify_CL2D_core_analysis.cpp \
../libraries/parallel/mpi_classify_FTTRI.cpp \
../libraries/parallel/mpi_image_rotational_pca.cpp \
../libraries/parallel/mpi_ml_align2d.cpp \
../libraries/parallel/mpi_performance_test.cpp \
../libraries/parallel/mpi_project_XR.cpp \
../libraries/parallel/mpi_reconstruct_art.cpp \
../libraries/parallel/mpi_reconstruct_fourier.cpp \
../libraries/parallel/mpi_reconstruct_wbp.cpp \
../libraries/parallel/xmipp_mpi.cpp 

OBJS += \
./libraries/parallel/mpi_angular_class_average.o \
./libraries/parallel/mpi_angular_gcar_commonlines.o \
./libraries/parallel/mpi_angular_projection_matching.o \
./libraries/parallel/mpi_classify_CL2D_core_analysis.o \
./libraries/parallel/mpi_classify_FTTRI.o \
./libraries/parallel/mpi_image_rotational_pca.o \
./libraries/parallel/mpi_ml_align2d.o \
./libraries/parallel/mpi_performance_test.o \
./libraries/parallel/mpi_project_XR.o \
./libraries/parallel/mpi_reconstruct_art.o \
./libraries/parallel/mpi_reconstruct_fourier.o \
./libraries/parallel/mpi_reconstruct_wbp.o \
./libraries/parallel/xmipp_mpi.o 

CPP_DEPS += \
./libraries/parallel/mpi_angular_class_average.d \
./libraries/parallel/mpi_angular_gcar_commonlines.d \
./libraries/parallel/mpi_angular_projection_matching.d \
./libraries/parallel/mpi_classify_CL2D_core_analysis.d \
./libraries/parallel/mpi_classify_FTTRI.d \
./libraries/parallel/mpi_image_rotational_pca.d \
./libraries/parallel/mpi_ml_align2d.d \
./libraries/parallel/mpi_performance_test.d \
./libraries/parallel/mpi_project_XR.d \
./libraries/parallel/mpi_reconstruct_art.d \
./libraries/parallel/mpi_reconstruct_fourier.d \
./libraries/parallel/mpi_reconstruct_wbp.d \
./libraries/parallel/xmipp_mpi.d 


# Each subdirectory must supply rules for building sources it contributes
libraries/parallel/%.o: ../libraries/parallel/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


