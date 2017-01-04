################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_UPPER_SRCS += \
../Macro/root/AutunniteRoot.C \
../Macro/root/CampioniRoot.C \
../Macro/root/Funzioni.C \
../Macro/root/RispostaRoot.C \
../Macro/root/Sostituto.C \
../Macro/root/asc2root.C \
../Macro/root/asctoroot.C 

OBJS += \
./Macro/root/AutunniteRoot.o \
./Macro/root/CampioniRoot.o \
./Macro/root/Funzioni.o \
./Macro/root/RispostaRoot.o \
./Macro/root/Sostituto.o \
./Macro/root/asc2root.o \
./Macro/root/asctoroot.o 

C_UPPER_DEPS += \
./Macro/root/AutunniteRoot.d \
./Macro/root/CampioniRoot.d \
./Macro/root/Funzioni.d \
./Macro/root/RispostaRoot.d \
./Macro/root/Sostituto.d \
./Macro/root/asc2root.d \
./Macro/root/asctoroot.d 


# Each subdirectory must supply rules for building sources it contributes
Macro/root/%.o: ../Macro/root/%.C
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


