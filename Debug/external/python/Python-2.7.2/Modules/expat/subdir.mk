################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/expat/xmlparse.c \
../external/python/Python-2.7.2/Modules/expat/xmlrole.c \
../external/python/Python-2.7.2/Modules/expat/xmltok.c \
../external/python/Python-2.7.2/Modules/expat/xmltok_impl.c \
../external/python/Python-2.7.2/Modules/expat/xmltok_ns.c 

OBJS += \
./external/python/Python-2.7.2/Modules/expat/xmlparse.o \
./external/python/Python-2.7.2/Modules/expat/xmlrole.o \
./external/python/Python-2.7.2/Modules/expat/xmltok.o \
./external/python/Python-2.7.2/Modules/expat/xmltok_impl.o \
./external/python/Python-2.7.2/Modules/expat/xmltok_ns.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/expat/xmlparse.d \
./external/python/Python-2.7.2/Modules/expat/xmlrole.d \
./external/python/Python-2.7.2/Modules/expat/xmltok.d \
./external/python/Python-2.7.2/Modules/expat/xmltok_impl.d \
./external/python/Python-2.7.2/Modules/expat/xmltok_ns.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/expat/%.o: ../external/python/Python-2.7.2/Modules/expat/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


