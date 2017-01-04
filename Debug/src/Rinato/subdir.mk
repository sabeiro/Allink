################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Rinato/Cubo.cpp \
../src/Rinato/Disegna.cpp \
../src/Rinato/Elabora.cpp \
../src/Rinato/Palline.cpp \
../src/Rinato/Particella.cpp \
../src/Rinato/Rinato.cpp \
../src/Rinato/Variabili.cpp 

C_SRCS += \
../src/Rinato/ran1.c 

OBJS += \
./src/Rinato/Cubo.o \
./src/Rinato/Disegna.o \
./src/Rinato/Elabora.o \
./src/Rinato/Palline.o \
./src/Rinato/Particella.o \
./src/Rinato/Rinato.o \
./src/Rinato/Variabili.o \
./src/Rinato/ran1.o 

C_DEPS += \
./src/Rinato/ran1.d 

CPP_DEPS += \
./src/Rinato/Cubo.d \
./src/Rinato/Disegna.d \
./src/Rinato/Elabora.d \
./src/Rinato/Palline.d \
./src/Rinato/Particella.d \
./src/Rinato/Rinato.d \
./src/Rinato/Variabili.d 


# Each subdirectory must supply rules for building sources it contributes
src/Rinato/%.o: ../src/Rinato/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/Rinato/%.o: ../src/Rinato/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


