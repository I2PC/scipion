################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/numpy-1.6.1/numpy/distutils/tests/swig_ext/src/example.c 

CC_SRCS += \
../external/python/numpy-1.6.1/numpy/distutils/tests/swig_ext/src/zoo.cc 

OBJS += \
./external/python/numpy-1.6.1/numpy/distutils/tests/swig_ext/src/example.o \
./external/python/numpy-1.6.1/numpy/distutils/tests/swig_ext/src/zoo.o 

C_DEPS += \
./external/python/numpy-1.6.1/numpy/distutils/tests/swig_ext/src/example.d 

CC_DEPS += \
./external/python/numpy-1.6.1/numpy/distutils/tests/swig_ext/src/zoo.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/numpy-1.6.1/numpy/distutils/tests/swig_ext/src/%.o: ../external/python/numpy-1.6.1/numpy/distutils/tests/swig_ext/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

external/python/numpy-1.6.1/numpy/distutils/tests/swig_ext/src/%.o: ../external/python/numpy-1.6.1/numpy/distutils/tests/swig_ext/src/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


