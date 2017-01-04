################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../src/Dinamica/Dinamica.o \
../src/Dinamica/Forces.o \
../src/Dinamica/ForcesBoundary.o \
../src/Dinamica/ForcesCalcEnergies.o \
../src/Dinamica/ForcesCreate.o \
../src/Dinamica/ForcesDraw.o \
../src/Dinamica/ForcesForceField.o \
../src/Dinamica/ForcesIntegration.o \
../src/Dinamica/ForcesLoop.o \
../src/Dinamica/ForcesTens.o 

CPP_SRCS += \
../src/Dinamica/Dinamica.cpp \
../src/Dinamica/Forces.cpp \
../src/Dinamica/ForcesBoundary.cpp \
../src/Dinamica/ForcesCalcEnergies.cpp \
../src/Dinamica/ForcesCreate.cpp \
../src/Dinamica/ForcesDraw.cpp \
../src/Dinamica/ForcesForceField.cpp \
../src/Dinamica/ForcesIntegration.cpp \
../src/Dinamica/ForcesLoop.cpp \
../src/Dinamica/ForcesTens.cpp 

C_SRCS += \
../src/Dinamica/Prova.c 

OBJS += \
./src/Dinamica/Dinamica.o \
./src/Dinamica/Forces.o \
./src/Dinamica/ForcesBoundary.o \
./src/Dinamica/ForcesCalcEnergies.o \
./src/Dinamica/ForcesCreate.o \
./src/Dinamica/ForcesDraw.o \
./src/Dinamica/ForcesForceField.o \
./src/Dinamica/ForcesIntegration.o \
./src/Dinamica/ForcesLoop.o \
./src/Dinamica/ForcesTens.o \
./src/Dinamica/Prova.o 

C_DEPS += \
./src/Dinamica/Prova.d 

CPP_DEPS += \
./src/Dinamica/Dinamica.d \
./src/Dinamica/Forces.d \
./src/Dinamica/ForcesBoundary.d \
./src/Dinamica/ForcesCalcEnergies.d \
./src/Dinamica/ForcesCreate.d \
./src/Dinamica/ForcesDraw.d \
./src/Dinamica/ForcesForceField.d \
./src/Dinamica/ForcesIntegration.d \
./src/Dinamica/ForcesLoop.d \
./src/Dinamica/ForcesTens.d 


# Each subdirectory must supply rules for building sources it contributes
src/Dinamica/%.o: ../src/Dinamica/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/Dinamica/%.o: ../src/Dinamica/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


