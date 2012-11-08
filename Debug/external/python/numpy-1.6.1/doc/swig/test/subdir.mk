################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CXX_SRCS += \
../external/python/numpy-1.6.1/doc/swig/test/Array1.cxx \
../external/python/numpy-1.6.1/doc/swig/test/Array2.cxx \
../external/python/numpy-1.6.1/doc/swig/test/Farray.cxx \
../external/python/numpy-1.6.1/doc/swig/test/Fortran.cxx \
../external/python/numpy-1.6.1/doc/swig/test/Matrix.cxx \
../external/python/numpy-1.6.1/doc/swig/test/Tensor.cxx \
../external/python/numpy-1.6.1/doc/swig/test/Vector.cxx 

OBJS += \
./external/python/numpy-1.6.1/doc/swig/test/Array1.o \
./external/python/numpy-1.6.1/doc/swig/test/Array2.o \
./external/python/numpy-1.6.1/doc/swig/test/Farray.o \
./external/python/numpy-1.6.1/doc/swig/test/Fortran.o \
./external/python/numpy-1.6.1/doc/swig/test/Matrix.o \
./external/python/numpy-1.6.1/doc/swig/test/Tensor.o \
./external/python/numpy-1.6.1/doc/swig/test/Vector.o 

CXX_DEPS += \
./external/python/numpy-1.6.1/doc/swig/test/Array1.d \
./external/python/numpy-1.6.1/doc/swig/test/Array2.d \
./external/python/numpy-1.6.1/doc/swig/test/Farray.d \
./external/python/numpy-1.6.1/doc/swig/test/Fortran.d \
./external/python/numpy-1.6.1/doc/swig/test/Matrix.d \
./external/python/numpy-1.6.1/doc/swig/test/Tensor.d \
./external/python/numpy-1.6.1/doc/swig/test/Vector.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/numpy-1.6.1/doc/swig/test/%.o: ../external/python/numpy-1.6.1/doc/swig/test/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


