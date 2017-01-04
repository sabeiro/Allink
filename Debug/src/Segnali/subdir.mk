################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Segnali/Segnali.cpp \
../src/Segnali/VarSegnali.cpp 

OBJS += \
./src/Segnali/Segnali.o \
./src/Segnali/VarSegnali.o 

CPP_DEPS += \
./src/Segnali/Segnali.d \
./src/Segnali/VarSegnali.d 


# Each subdirectory must supply rules for building sources it contributes
src/Segnali/%.o: ../src/Segnali/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


