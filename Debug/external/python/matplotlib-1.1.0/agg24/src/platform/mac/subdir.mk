################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../external/python/matplotlib-1.1.0/agg24/src/platform/mac/agg_mac_pmap.cpp \
../external/python/matplotlib-1.1.0/agg24/src/platform/mac/agg_platform_support.cpp 

OBJS += \
./external/python/matplotlib-1.1.0/agg24/src/platform/mac/agg_mac_pmap.o \
./external/python/matplotlib-1.1.0/agg24/src/platform/mac/agg_platform_support.o 

CPP_DEPS += \
./external/python/matplotlib-1.1.0/agg24/src/platform/mac/agg_mac_pmap.d \
./external/python/matplotlib-1.1.0/agg24/src/platform/mac/agg_platform_support.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/matplotlib-1.1.0/agg24/src/platform/mac/%.o: ../external/python/matplotlib-1.1.0/agg24/src/platform/mac/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


