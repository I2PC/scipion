################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../external/tiff-3.9.4/contrib/addtiffo/addtiffo.o \
../external/tiff-3.9.4/contrib/addtiffo/tif_overview.o \
../external/tiff-3.9.4/contrib/addtiffo/tif_ovrcache.o 

C_SRCS += \
../external/tiff-3.9.4/contrib/addtiffo/addtiffo.c \
../external/tiff-3.9.4/contrib/addtiffo/tif_overview.c \
../external/tiff-3.9.4/contrib/addtiffo/tif_ovrcache.c 

OBJS += \
./external/tiff-3.9.4/contrib/addtiffo/addtiffo.o \
./external/tiff-3.9.4/contrib/addtiffo/tif_overview.o \
./external/tiff-3.9.4/contrib/addtiffo/tif_ovrcache.o 

C_DEPS += \
./external/tiff-3.9.4/contrib/addtiffo/addtiffo.d \
./external/tiff-3.9.4/contrib/addtiffo/tif_overview.d \
./external/tiff-3.9.4/contrib/addtiffo/tif_ovrcache.d 


# Each subdirectory must supply rules for building sources it contributes
external/tiff-3.9.4/contrib/addtiffo/%.o: ../external/tiff-3.9.4/contrib/addtiffo/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


