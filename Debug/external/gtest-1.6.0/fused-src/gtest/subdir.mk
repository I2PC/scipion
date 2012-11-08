################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../external/gtest-1.6.0/fused-src/gtest/gtest-all.cc \
../external/gtest-1.6.0/fused-src/gtest/gtest_main.cc 

OBJS += \
./external/gtest-1.6.0/fused-src/gtest/gtest-all.o \
./external/gtest-1.6.0/fused-src/gtest/gtest_main.o 

CC_DEPS += \
./external/gtest-1.6.0/fused-src/gtest/gtest-all.d \
./external/gtest-1.6.0/fused-src/gtest/gtest_main.d 


# Each subdirectory must supply rules for building sources it contributes
external/gtest-1.6.0/fused-src/gtest/%.o: ../external/gtest-1.6.0/fused-src/gtest/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


