################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/tcl8.5.10/unix/dltest/pkga.c \
../external/python/tcl8.5.10/unix/dltest/pkgb.c \
../external/python/tcl8.5.10/unix/dltest/pkgc.c \
../external/python/tcl8.5.10/unix/dltest/pkgd.c \
../external/python/tcl8.5.10/unix/dltest/pkge.c \
../external/python/tcl8.5.10/unix/dltest/pkgua.c 

OBJS += \
./external/python/tcl8.5.10/unix/dltest/pkga.o \
./external/python/tcl8.5.10/unix/dltest/pkgb.o \
./external/python/tcl8.5.10/unix/dltest/pkgc.o \
./external/python/tcl8.5.10/unix/dltest/pkgd.o \
./external/python/tcl8.5.10/unix/dltest/pkge.o \
./external/python/tcl8.5.10/unix/dltest/pkgua.o 

C_DEPS += \
./external/python/tcl8.5.10/unix/dltest/pkga.d \
./external/python/tcl8.5.10/unix/dltest/pkgb.d \
./external/python/tcl8.5.10/unix/dltest/pkgc.d \
./external/python/tcl8.5.10/unix/dltest/pkgd.d \
./external/python/tcl8.5.10/unix/dltest/pkge.d \
./external/python/tcl8.5.10/unix/dltest/pkgua.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/tcl8.5.10/unix/dltest/%.o: ../external/python/tcl8.5.10/unix/dltest/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


