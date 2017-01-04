################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../src/ElPoly/ElPoly.o \
../src/ElPoly/ElPolyDraw.o \
../src/ElPoly/ElPolyDrawControl.o \
../src/ElPoly/ElPolyDrawSkinCGAL.o \
../src/ElPoly/ElPolyDrawSurf.o \
../src/ElPoly/ElPolyDrawTriangCGAL.o \
../src/ElPoly/ElPolyEl.o \
../src/ElPoly/ElPolyMeasure.o \
../src/ElPoly/ElPolyOutput.o \
../src/ElPoly/ElPolyProfDens.o \
../src/ElPoly/ElPolyProfPre.o \
../src/ElPoly/ElPolyProfile.o \
../src/ElPoly/ElPolyRepr.o 

CPP_SRCS += \
../src/ElPoly/ElPoly.cpp \
../src/ElPoly/ElPolyDraw.cpp \
../src/ElPoly/ElPolyDrawCGAL.cpp \
../src/ElPoly/ElPolyDrawControl.cpp \
../src/ElPoly/ElPolyDrawDef.cpp \
../src/ElPoly/ElPolyDrawSkinCGAL.cpp \
../src/ElPoly/ElPolyDrawSurf.cpp \
../src/ElPoly/ElPolyDrawSurfCGAL.cpp \
../src/ElPoly/ElPolyDrawSurfMarchCubes.cpp \
../src/ElPoly/ElPolyDrawTriangCGAL.cpp \
../src/ElPoly/ElPolyEl.cpp \
../src/ElPoly/ElPolyMeasure.cpp \
../src/ElPoly/ElPolyOutput.cpp \
../src/ElPoly/ElPolyProfDens.cpp \
../src/ElPoly/ElPolyProfPre.cpp \
../src/ElPoly/ElPolyProfile.cpp \
../src/ElPoly/ElPolyRepr.cpp 

OBJS += \
./src/ElPoly/ElPoly.o \
./src/ElPoly/ElPolyDraw.o \
./src/ElPoly/ElPolyDrawCGAL.o \
./src/ElPoly/ElPolyDrawControl.o \
./src/ElPoly/ElPolyDrawDef.o \
./src/ElPoly/ElPolyDrawSkinCGAL.o \
./src/ElPoly/ElPolyDrawSurf.o \
./src/ElPoly/ElPolyDrawSurfCGAL.o \
./src/ElPoly/ElPolyDrawSurfMarchCubes.o \
./src/ElPoly/ElPolyDrawTriangCGAL.o \
./src/ElPoly/ElPolyEl.o \
./src/ElPoly/ElPolyMeasure.o \
./src/ElPoly/ElPolyOutput.o \
./src/ElPoly/ElPolyProfDens.o \
./src/ElPoly/ElPolyProfPre.o \
./src/ElPoly/ElPolyProfile.o \
./src/ElPoly/ElPolyRepr.o 

CPP_DEPS += \
./src/ElPoly/ElPoly.d \
./src/ElPoly/ElPolyDraw.d \
./src/ElPoly/ElPolyDrawCGAL.d \
./src/ElPoly/ElPolyDrawControl.d \
./src/ElPoly/ElPolyDrawDef.d \
./src/ElPoly/ElPolyDrawSkinCGAL.d \
./src/ElPoly/ElPolyDrawSurf.d \
./src/ElPoly/ElPolyDrawSurfCGAL.d \
./src/ElPoly/ElPolyDrawSurfMarchCubes.d \
./src/ElPoly/ElPolyDrawTriangCGAL.d \
./src/ElPoly/ElPolyEl.d \
./src/ElPoly/ElPolyMeasure.d \
./src/ElPoly/ElPolyOutput.d \
./src/ElPoly/ElPolyProfDens.d \
./src/ElPoly/ElPolyProfPre.d \
./src/ElPoly/ElPolyProfile.d \
./src/ElPoly/ElPolyRepr.d 


# Each subdirectory must supply rules for building sources it contributes
src/ElPoly/%.o: ../src/ElPoly/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


