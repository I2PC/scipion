################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../external/python/matplotlib-1.1.0/agg24/src/platform/win32/agg_platform_support.cpp \
../external/python/matplotlib-1.1.0/agg24/src/platform/win32/agg_win32_bmp.cpp 

OBJS += \
./external/python/matplotlib-1.1.0/agg24/src/platform/win32/agg_platform_support.o \
./external/python/matplotlib-1.1.0/agg24/src/platform/win32/agg_win32_bmp.o 

CPP_DEPS += \
./external/python/matplotlib-1.1.0/agg24/src/platform/win32/agg_platform_support.d \
./external/python/matplotlib-1.1.0/agg24/src/platform/win32/agg_win32_bmp.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/matplotlib-1.1.0/agg24/src/platform/win32/%.o: ../external/python/matplotlib-1.1.0/agg24/src/platform/win32/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


