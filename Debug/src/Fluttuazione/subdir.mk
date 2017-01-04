################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Fluttuazione/Fluttuazione.cpp \
../src/Fluttuazione/VarFluttuazione.cpp 

OBJS += \
./src/Fluttuazione/Fluttuazione.o \
./src/Fluttuazione/VarFluttuazione.o 

CPP_DEPS += \
./src/Fluttuazione/Fluttuazione.d \
./src/Fluttuazione/VarFluttuazione.d 


# Each subdirectory must supply rules for building sources it contributes
src/Fluttuazione/%.o: ../src/Fluttuazione/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


