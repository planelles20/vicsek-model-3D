LINKER_FLAGS = -lGLEW -lglfw3 -lGL -lX11 -lXi -lXrandr -lXxf86vm -lXinerama -lXcursor -ldl -Wno-deprecated-gpu-targets #-pthread
COMPILER_FLAGS = -std=c++11

# The compiler
NC = nvcc

all:
	$(NC) $(COMPILER_FLAGS) -o main.exe main.cpp ./common/display.cpp ./common/shaderGeometry.cpp ./common/shader.cpp systemSPH.cu $(LINKER_FLAGS)

clean:
	rm *.exe
