################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/matplotlib-1.1.0/CXX/Python2/cxxextensions.c 

CXX_SRCS += \
../external/python/matplotlib-1.1.0/CXX/Python2/IndirectPythonInterface.cxx \
../external/python/matplotlib-1.1.0/CXX/Python2/cxx_extensions.cxx \
../external/python/matplotlib-1.1.0/CXX/Python2/cxxsupport.cxx 

OBJS += \
./external/python/matplotlib-1.1.0/CXX/Python2/IndirectPythonInterface.o \
./external/python/matplotlib-1.1.0/CXX/Python2/cxx_extensions.o \
./external/python/matplotlib-1.1.0/CXX/Python2/cxxextensions.o \
./external/python/matplotlib-1.1.0/CXX/Python2/cxxsupport.o 

C_DEPS += \
./external/python/matplotlib-1.1.0/CXX/Python2/cxxextensions.d 

CXX_DEPS += \
./external/python/matplotlib-1.1.0/CXX/Python2/IndirectPythonInterface.d \
./external/python/matplotlib-1.1.0/CXX/Python2/cxx_extensions.d \
./external/python/matplotlib-1.1.0/CXX/Python2/cxxsupport.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/matplotlib-1.1.0/CXX/Python2/%.o: ../external/python/matplotlib-1.1.0/CXX/Python2/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

external/python/matplotlib-1.1.0/CXX/Python2/%.o: ../external/python/matplotlib-1.1.0/CXX/Python2/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


