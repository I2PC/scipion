################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/tiff-3.9.4/contrib/pds/tif_imageiter.c \
../external/tiff-3.9.4/contrib/pds/tif_pdsdirread.c \
../external/tiff-3.9.4/contrib/pds/tif_pdsdirwrite.c 

OBJS += \
./external/tiff-3.9.4/contrib/pds/tif_imageiter.o \
./external/tiff-3.9.4/contrib/pds/tif_pdsdirread.o \
./external/tiff-3.9.4/contrib/pds/tif_pdsdirwrite.o 

C_DEPS += \
./external/tiff-3.9.4/contrib/pds/tif_imageiter.d \
./external/tiff-3.9.4/contrib/pds/tif_pdsdirread.d \
./external/tiff-3.9.4/contrib/pds/tif_pdsdirwrite.d 


# Each subdirectory must supply rules for building sources it contributes
external/tiff-3.9.4/contrib/pds/%.o: ../external/tiff-3.9.4/contrib/pds/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


