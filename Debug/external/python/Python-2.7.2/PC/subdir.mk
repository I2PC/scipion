################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/PC/WinMain.c \
../external/python/Python-2.7.2/PC/_msi.c \
../external/python/Python-2.7.2/PC/_subprocess.c \
../external/python/Python-2.7.2/PC/_winreg.c \
../external/python/Python-2.7.2/PC/config.c \
../external/python/Python-2.7.2/PC/dl_nt.c \
../external/python/Python-2.7.2/PC/empty.c \
../external/python/Python-2.7.2/PC/frozen_dllmain.c \
../external/python/Python-2.7.2/PC/generrmap.c \
../external/python/Python-2.7.2/PC/getpathp.c \
../external/python/Python-2.7.2/PC/import_nt.c \
../external/python/Python-2.7.2/PC/make_versioninfo.c \
../external/python/Python-2.7.2/PC/msvcrtmodule.c \
../external/python/Python-2.7.2/PC/w9xpopen.c \
../external/python/Python-2.7.2/PC/winsound.c 

OBJS += \
./external/python/Python-2.7.2/PC/WinMain.o \
./external/python/Python-2.7.2/PC/_msi.o \
./external/python/Python-2.7.2/PC/_subprocess.o \
./external/python/Python-2.7.2/PC/_winreg.o \
./external/python/Python-2.7.2/PC/config.o \
./external/python/Python-2.7.2/PC/dl_nt.o \
./external/python/Python-2.7.2/PC/empty.o \
./external/python/Python-2.7.2/PC/frozen_dllmain.o \
./external/python/Python-2.7.2/PC/generrmap.o \
./external/python/Python-2.7.2/PC/getpathp.o \
./external/python/Python-2.7.2/PC/import_nt.o \
./external/python/Python-2.7.2/PC/make_versioninfo.o \
./external/python/Python-2.7.2/PC/msvcrtmodule.o \
./external/python/Python-2.7.2/PC/w9xpopen.o \
./external/python/Python-2.7.2/PC/winsound.o 

C_DEPS += \
./external/python/Python-2.7.2/PC/WinMain.d \
./external/python/Python-2.7.2/PC/_msi.d \
./external/python/Python-2.7.2/PC/_subprocess.d \
./external/python/Python-2.7.2/PC/_winreg.d \
./external/python/Python-2.7.2/PC/config.d \
./external/python/Python-2.7.2/PC/dl_nt.d \
./external/python/Python-2.7.2/PC/empty.d \
./external/python/Python-2.7.2/PC/frozen_dllmain.d \
./external/python/Python-2.7.2/PC/generrmap.d \
./external/python/Python-2.7.2/PC/getpathp.d \
./external/python/Python-2.7.2/PC/import_nt.d \
./external/python/Python-2.7.2/PC/make_versioninfo.d \
./external/python/Python-2.7.2/PC/msvcrtmodule.d \
./external/python/Python-2.7.2/PC/w9xpopen.d \
./external/python/Python-2.7.2/PC/winsound.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/PC/%.o: ../external/python/Python-2.7.2/PC/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


