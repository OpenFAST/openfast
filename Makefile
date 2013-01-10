all: gnb

gnb: quadrature.o dlegen.o dlagint.o NodeLoc.o nodeMat.o Amatrix.o Dmatrix.o Chimatrix.o RHS_element.o AssembleRHS.o Norm.o Ximatrix.o KT_element.o AssembleKT.o CGSolver.o LineSearch.o Newton.o driver.o
	gfortran -O3 -o gnb quadrature.o dlegen.o dlagint.o NodeLoc.o nodeMat.o Amatrix.o Dmatrix.o Chimatrix.o RHS_element.o AssembleRHS.o Norm.o Ximatrix.o KT_element.o AssembleKT.o CGSolver.o LineSearch.o Newton.o driver.o
	
quadrature.o: quadrature.f90
	gfortran -O3 -c quadrature.f90

dlegen.o: dlegen.f90
	gfortran -O3 -c dlegen.f90
		
dlagint.o: dlagint.f90
	gfortran -O3 -c dlagint.f90
	
NodeLoc.o: NodeLoc.f90
	gfortran -O3 -c NodeLoc.f90
	
NodeMat.o: NodeMat.f90
	gfortran -O3 -c NodeMat.f90
	
Amatrix.o: Amatrix.f90
	gfortran -O3 -c Amatrix.f90
	
Dmatrix.o: Dmatrix.f90
	gfortran -O3 -c Dmatrix.f90
	
Chimatrix.o: Chimatrix.f90
	gfortran -O3 -c Chimatrix.f90
	
RHS_element.o: RHS_element.f90
	gfortran -O3 -c RHS_element.f90
	
AssembleRHS.o: AssembleRHS.f90
	gfortran -O3 -c AssembleRHS.f90
	
Norm.o: Norm.f90
	gfortran -O3 -c Norm.f90
	
Ximatrix.o: Ximatrix.f90
	gfortran -O3 -c Ximatrix.f90
	
KT_element.o: KT_element.f90
	gfortran -O3 -c KT_element.f90
	
AssembleKT.o: AssembleKT.f90
	gfortran -O3 -c AssembleKT.f90
	
CGSolver.o: CGSolver.f90
	gfortran -O3 -c CGSolver.f90
	
LineSearch.o: LineSearch.f90
	gfortran -O3 -c LineSearch.f90
	
Newton.o: Newton.f90
	gfortran -O3 -c Newton.f90
	
driver.o: driver.f90
	gfortran -O3 -c driver.f90