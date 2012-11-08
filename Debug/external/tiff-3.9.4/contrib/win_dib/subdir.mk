################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../external/tiff-3.9.4/contrib/win_dib/Tiffile.cpp 

C_SRCS += \
../external/tiff-3.9.4/contrib/win_dib/tiff2dib.c 

OBJS += \
./external/tiff-3.9.4/contrib/win_dib/Tiffile.o \
./external/tiff-3.9.4/contrib/win_dib/tiff2dib.o 

C_DEPS += \
./external/tiff-3.9.4/contrib/win_dib/tiff2dib.d 

CPP_DEPS += \
./external/tiff-3.9.4/contrib/win_dib/Tiffile.d 


# Each subdirectory must supply rules for building sources it contributes
external/tiff-3.9.4/contrib/win_dib/%.o: ../external/tiff-3.9.4/contrib/win_dib/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

external/tiff-3.9.4/contrib/win_dib/%.o: ../external/tiff-3.9.4/contrib/win_dib/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


