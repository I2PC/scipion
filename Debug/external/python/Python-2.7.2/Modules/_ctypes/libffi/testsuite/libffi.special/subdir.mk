################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/testsuite/libffi.special/unwindtest.cc \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/testsuite/libffi.special/unwindtest_ffi_call.cc 

OBJS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/testsuite/libffi.special/unwindtest.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/testsuite/libffi.special/unwindtest_ffi_call.o 

CC_DEPS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/testsuite/libffi.special/unwindtest.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/testsuite/libffi.special/unwindtest_ffi_call.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/_ctypes/libffi/testsuite/libffi.special/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi/testsuite/libffi.special/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


