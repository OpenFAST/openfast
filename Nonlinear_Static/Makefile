all: gnb

#gnb: gen_gll.o gen_deriv.o NodeLoc.o NodeMat.o Amatrix.o Dmatrix.o Chimatrix.o RHS_element.o AssembleRHS.o Norm.o Ximatrix.o KT_element.o AssembleKT.o CGSolver.o LineSearch.o NewtonRaphson.o driver.o
#	gfortran -fdefault-real-8 -fdefault-double-8 -o gnb gen_gll.o gen_deriv.o NodeLoc.o nodeMat.o Amatrix.o Dmatrix.o Chimatrix.o RHS_element.o AssembleRHS.o Norm.o Ximatrix.o KT_element.o AssembleKT.o CGSolver.o LineSearch.o NewtonRaphson.o driver.o

gnb: gen_gll.o gen_deriv.o NodeLoc.o NodeMat.o Amatrix.o Dmatrix.o Chimatrix.o RHS_element.o AssembleRHS.o Norm.o Ximatrix.o KT_element.o AssembleKT.o CGSolver.o Newton_New.o filename.o outputvtk.o driver.o
	gfortran -fdefault-real-8 -fdefault-double-8 -o gnb gen_gll.o gen_deriv.o NodeLoc.o nodeMat.o Amatrix.o Dmatrix.o Chimatrix.o RHS_element.o AssembleRHS.o Norm.o Ximatrix.o KT_element.o AssembleKT.o CGSolver.o Newton_New.o filename.o outputvtk.o driver.o
	
gen_gll.o: gen_gll.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c gen_gll.f90

gen_deriv.o: gen_deriv.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c gen_deriv.f90
		
dlagint.o: dlagint.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c dlagint.f90
	
NodeLoc.o: NodeLoc.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c NodeLoc.f90
	
NodeMat.o: NodeMat.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c NodeMat.f90
	
Amatrix.o: Amatrix.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c Amatrix.f90
	
Dmatrix.o: Dmatrix.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c Dmatrix.f90
	
Chimatrix.o: Chimatrix.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c Chimatrix.f90
	
RHS_element.o: RHS_element.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c RHS_element.f90
	
AssembleRHS.o: AssembleRHS.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c AssembleRHS.f90
	
Norm.o: Norm.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c Norm.f90
	
Ximatrix.o: Ximatrix.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c Ximatrix.f90
	
KT_element.o: KT_element.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c KT_element.f90
	
AssembleKT.o: AssembleKT.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c AssembleKT.f90
	
CGSolver.o: CGSolver.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c CGSolver.f90
	
#LineSearch.o: LineSearch.f90
#	gfortran -fdefault-real-8 -fdefault-double-8 -c LineSearch.f90
	
#NewtonRaphson.o: NewtonRaphson.f90
#	gfortran -fdefault-real-8 -fdefault-double-8 -c NewtonRaphson.f90

filename.o: filename.f
	gfortran -fdefault-real-8 -fdefault-double-8 -c filename.f

outputvtk.o: outputvtk.f
	gfortran -fdefault-real-8 -fdefault-double-8 -c outputvtk.f

Newton_New.o: Newton_New.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c Newton_New.f90
	
driver.o: driver.f90
	gfortran -fdefault-real-8 -fdefault-double-8 -c driver.f90
