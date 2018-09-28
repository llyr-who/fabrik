all: GLnixDEMO2

#For debugging
OPT=-g -Wall
LDFLAGS=-lm -lX11 -lGL -lGLU -lXext -lXrender
OBJECTS=GLnixAPP.o myextloader.h GeometryGen.o MathHelper.o cloth.o antmath.o Timer.o
		
#For optimistaion
#OPT=-O

#All objects (except main) come from cpp and hpp 
%.o:	%.cpp %.hpp
	g++ ${OPT} -c -o $@ $<
#use_vectors relies on objects which rely on headers
GLnixDEMO2:	GLnixDEMO2.cpp ${OBJECTS}
		g++ ${OPT} ${OBJECTS} -o GLnixDEMO2 GLnixDEMO2.cpp ${LDFLAGS}
clean:
	rm -f *.o *~ GLnixDEMO2
