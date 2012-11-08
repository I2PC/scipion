################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../external/inria/convert.cc \
../external/inria/curvature.cc \
../external/inria/extrema.cc \
../external/inria/recbuffer.cc \
../external/inria/recline.cc 

OBJS += \
./external/inria/convert.o \
./external/inria/curvature.o \
./external/inria/extrema.o \
./external/inria/recbuffer.o \
./external/inria/recline.o 

CC_DEPS += \
./external/inria/convert.d \
./external/inria/curvature.d \
./external/inria/extrema.d \
./external/inria/recbuffer.d \
./external/inria/recline.d 


# Each subdirectory must supply rules for building sources it contributes
external/inria/%.o: ../external/inria/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


