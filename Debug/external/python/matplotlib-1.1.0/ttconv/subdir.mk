################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../external/python/matplotlib-1.1.0/ttconv/pprdrv_tt.cpp \
../external/python/matplotlib-1.1.0/ttconv/pprdrv_tt2.cpp \
../external/python/matplotlib-1.1.0/ttconv/ttutil.cpp 

OBJS += \
./external/python/matplotlib-1.1.0/ttconv/pprdrv_tt.o \
./external/python/matplotlib-1.1.0/ttconv/pprdrv_tt2.o \
./external/python/matplotlib-1.1.0/ttconv/ttutil.o 

CPP_DEPS += \
./external/python/matplotlib-1.1.0/ttconv/pprdrv_tt.d \
./external/python/matplotlib-1.1.0/ttconv/pprdrv_tt2.d \
./external/python/matplotlib-1.1.0/ttconv/ttutil.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/matplotlib-1.1.0/ttconv/%.o: ../external/python/matplotlib-1.1.0/ttconv/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


