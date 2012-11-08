################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../libraries/bindings/matlab/tom_xmipp_adjust_ctf_wrapper.cpp \
../libraries/bindings/matlab/tom_xmipp_align2d_wrapper.cpp \
../libraries/bindings/matlab/tom_xmipp_ctf_correct_phase_wrapper.cpp \
../libraries/bindings/matlab/tom_xmipp_mask_wrapper.cpp \
../libraries/bindings/matlab/tom_xmipp_mirror_wrapper.cpp \
../libraries/bindings/matlab/tom_xmipp_morphology_wrapper.cpp \
../libraries/bindings/matlab/tom_xmipp_normalize_wrapper.cpp \
../libraries/bindings/matlab/tom_xmipp_psd_enhance_wrapper.cpp \
../libraries/bindings/matlab/tom_xmipp_resolution_wrapper.cpp \
../libraries/bindings/matlab/tom_xmipp_rotate_wrapper.cpp \
../libraries/bindings/matlab/tom_xmipp_scale_pyramid_wrapper.cpp \
../libraries/bindings/matlab/tom_xmipp_scale_wrapper.cpp \
../libraries/bindings/matlab/tom_xmipp_volume_segment_wrapper.cpp 

OBJS += \
./libraries/bindings/matlab/tom_xmipp_adjust_ctf_wrapper.o \
./libraries/bindings/matlab/tom_xmipp_align2d_wrapper.o \
./libraries/bindings/matlab/tom_xmipp_ctf_correct_phase_wrapper.o \
./libraries/bindings/matlab/tom_xmipp_mask_wrapper.o \
./libraries/bindings/matlab/tom_xmipp_mirror_wrapper.o \
./libraries/bindings/matlab/tom_xmipp_morphology_wrapper.o \
./libraries/bindings/matlab/tom_xmipp_normalize_wrapper.o \
./libraries/bindings/matlab/tom_xmipp_psd_enhance_wrapper.o \
./libraries/bindings/matlab/tom_xmipp_resolution_wrapper.o \
./libraries/bindings/matlab/tom_xmipp_rotate_wrapper.o \
./libraries/bindings/matlab/tom_xmipp_scale_pyramid_wrapper.o \
./libraries/bindings/matlab/tom_xmipp_scale_wrapper.o \
./libraries/bindings/matlab/tom_xmipp_volume_segment_wrapper.o 

CPP_DEPS += \
./libraries/bindings/matlab/tom_xmipp_adjust_ctf_wrapper.d \
./libraries/bindings/matlab/tom_xmipp_align2d_wrapper.d \
./libraries/bindings/matlab/tom_xmipp_ctf_correct_phase_wrapper.d \
./libraries/bindings/matlab/tom_xmipp_mask_wrapper.d \
./libraries/bindings/matlab/tom_xmipp_mirror_wrapper.d \
./libraries/bindings/matlab/tom_xmipp_morphology_wrapper.d \
./libraries/bindings/matlab/tom_xmipp_normalize_wrapper.d \
./libraries/bindings/matlab/tom_xmipp_psd_enhance_wrapper.d \
./libraries/bindings/matlab/tom_xmipp_resolution_wrapper.d \
./libraries/bindings/matlab/tom_xmipp_rotate_wrapper.d \
./libraries/bindings/matlab/tom_xmipp_scale_pyramid_wrapper.d \
./libraries/bindings/matlab/tom_xmipp_scale_wrapper.d \
./libraries/bindings/matlab/tom_xmipp_volume_segment_wrapper.d 


# Each subdirectory must supply rules for building sources it contributes
libraries/bindings/matlab/%.o: ../libraries/bindings/matlab/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


