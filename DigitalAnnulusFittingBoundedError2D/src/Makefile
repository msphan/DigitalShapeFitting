CXX := g++
OUT := DigitalAnnulusFitting
LIBS := -lglut -lGL -lGLU
$(OUT) : CircleFitting.o OctCube.o DPlane.o Plane.o Point3D.o
	$(CXX) $^ -o $(OUT) $(LIBS)
Point3D.o : Point3D.cpp
	$(CXX) -c $^ -o $@ $(LIBS)
Plane.o : Plane.cpp
	$(CXX) -c $^ -o $@ $(LIBS) 
DPlane.o : DPlane.cpp
	$(CXX) -c $^ -o $@ $(LIBS) 
OctCube.o : OctCube.cpp
	$(CXX) -c $^ -o $@ $(LIBS)
CircleFitting.o : CircleFitting.cpp
	$(CXX) -c $^ -o $@ $(LIBS) 

.PHONY : clean
clean ::
	 rm -f *.o
