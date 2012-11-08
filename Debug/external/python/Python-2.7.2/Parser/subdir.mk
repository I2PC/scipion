################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../external/python/Python-2.7.2/Parser/acceler.o \
../external/python/Python-2.7.2/Parser/bitset.o \
../external/python/Python-2.7.2/Parser/firstsets.o \
../external/python/Python-2.7.2/Parser/grammar.o \
../external/python/Python-2.7.2/Parser/grammar1.o \
../external/python/Python-2.7.2/Parser/listnode.o \
../external/python/Python-2.7.2/Parser/metagrammar.o \
../external/python/Python-2.7.2/Parser/myreadline.o \
../external/python/Python-2.7.2/Parser/node.o \
../external/python/Python-2.7.2/Parser/parser.o \
../external/python/Python-2.7.2/Parser/parsetok.o \
../external/python/Python-2.7.2/Parser/pgen.o \
../external/python/Python-2.7.2/Parser/pgenmain.o \
../external/python/Python-2.7.2/Parser/printgrammar.o \
../external/python/Python-2.7.2/Parser/tokenizer.o \
../external/python/Python-2.7.2/Parser/tokenizer_pgen.o 

C_SRCS += \
../external/python/Python-2.7.2/Parser/acceler.c \
../external/python/Python-2.7.2/Parser/bitset.c \
../external/python/Python-2.7.2/Parser/firstsets.c \
../external/python/Python-2.7.2/Parser/grammar.c \
../external/python/Python-2.7.2/Parser/grammar1.c \
../external/python/Python-2.7.2/Parser/intrcheck.c \
../external/python/Python-2.7.2/Parser/listnode.c \
../external/python/Python-2.7.2/Parser/metagrammar.c \
../external/python/Python-2.7.2/Parser/myreadline.c \
../external/python/Python-2.7.2/Parser/node.c \
../external/python/Python-2.7.2/Parser/parser.c \
../external/python/Python-2.7.2/Parser/parsetok.c \
../external/python/Python-2.7.2/Parser/pgen.c \
../external/python/Python-2.7.2/Parser/pgenmain.c \
../external/python/Python-2.7.2/Parser/printgrammar.c \
../external/python/Python-2.7.2/Parser/tokenizer.c \
../external/python/Python-2.7.2/Parser/tokenizer_pgen.c 

OBJS += \
./external/python/Python-2.7.2/Parser/acceler.o \
./external/python/Python-2.7.2/Parser/bitset.o \
./external/python/Python-2.7.2/Parser/firstsets.o \
./external/python/Python-2.7.2/Parser/grammar.o \
./external/python/Python-2.7.2/Parser/grammar1.o \
./external/python/Python-2.7.2/Parser/intrcheck.o \
./external/python/Python-2.7.2/Parser/listnode.o \
./external/python/Python-2.7.2/Parser/metagrammar.o \
./external/python/Python-2.7.2/Parser/myreadline.o \
./external/python/Python-2.7.2/Parser/node.o \
./external/python/Python-2.7.2/Parser/parser.o \
./external/python/Python-2.7.2/Parser/parsetok.o \
./external/python/Python-2.7.2/Parser/pgen.o \
./external/python/Python-2.7.2/Parser/pgenmain.o \
./external/python/Python-2.7.2/Parser/printgrammar.o \
./external/python/Python-2.7.2/Parser/tokenizer.o \
./external/python/Python-2.7.2/Parser/tokenizer_pgen.o 

C_DEPS += \
./external/python/Python-2.7.2/Parser/acceler.d \
./external/python/Python-2.7.2/Parser/bitset.d \
./external/python/Python-2.7.2/Parser/firstsets.d \
./external/python/Python-2.7.2/Parser/grammar.d \
./external/python/Python-2.7.2/Parser/grammar1.d \
./external/python/Python-2.7.2/Parser/intrcheck.d \
./external/python/Python-2.7.2/Parser/listnode.d \
./external/python/Python-2.7.2/Parser/metagrammar.d \
./external/python/Python-2.7.2/Parser/myreadline.d \
./external/python/Python-2.7.2/Parser/node.d \
./external/python/Python-2.7.2/Parser/parser.d \
./external/python/Python-2.7.2/Parser/parsetok.d \
./external/python/Python-2.7.2/Parser/pgen.d \
./external/python/Python-2.7.2/Parser/pgenmain.d \
./external/python/Python-2.7.2/Parser/printgrammar.d \
./external/python/Python-2.7.2/Parser/tokenizer.d \
./external/python/Python-2.7.2/Parser/tokenizer_pgen.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Parser/%.o: ../external/python/Python-2.7.2/Parser/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


