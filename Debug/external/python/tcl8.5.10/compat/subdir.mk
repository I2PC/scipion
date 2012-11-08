################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/tcl8.5.10/compat/fixstrtod.c \
../external/python/tcl8.5.10/compat/gettod.c \
../external/python/tcl8.5.10/compat/memcmp.c \
../external/python/tcl8.5.10/compat/opendir.c \
../external/python/tcl8.5.10/compat/strncasecmp.c \
../external/python/tcl8.5.10/compat/strstr.c \
../external/python/tcl8.5.10/compat/strtod.c \
../external/python/tcl8.5.10/compat/strtol.c \
../external/python/tcl8.5.10/compat/strtoul.c \
../external/python/tcl8.5.10/compat/waitpid.c 

OBJS += \
./external/python/tcl8.5.10/compat/fixstrtod.o \
./external/python/tcl8.5.10/compat/gettod.o \
./external/python/tcl8.5.10/compat/memcmp.o \
./external/python/tcl8.5.10/compat/opendir.o \
./external/python/tcl8.5.10/compat/strncasecmp.o \
./external/python/tcl8.5.10/compat/strstr.o \
./external/python/tcl8.5.10/compat/strtod.o \
./external/python/tcl8.5.10/compat/strtol.o \
./external/python/tcl8.5.10/compat/strtoul.o \
./external/python/tcl8.5.10/compat/waitpid.o 

C_DEPS += \
./external/python/tcl8.5.10/compat/fixstrtod.d \
./external/python/tcl8.5.10/compat/gettod.d \
./external/python/tcl8.5.10/compat/memcmp.d \
./external/python/tcl8.5.10/compat/opendir.d \
./external/python/tcl8.5.10/compat/strncasecmp.d \
./external/python/tcl8.5.10/compat/strstr.d \
./external/python/tcl8.5.10/compat/strtod.d \
./external/python/tcl8.5.10/compat/strtol.d \
./external/python/tcl8.5.10/compat/strtoul.d \
./external/python/tcl8.5.10/compat/waitpid.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/tcl8.5.10/compat/%.o: ../external/python/tcl8.5.10/compat/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


