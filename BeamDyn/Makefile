all: BeamStaticTest

BeamStaticTest: Tilde.o OuterProduct.o gen_gll.o gen_deriv.o CrvMatrixR.o CrvMatrixH.o CrvCompose.o NodeLoc.o ElemNodalDisp.o NodalRelRot.o NodalDataAt0.o NodalData.o ElasticForce.o ElementMatrix.o BeamStatic.o StaticSolution.o DriverStatic.o
    gfortran -O3 -o BeamStaticTest Tilde.o OuterProduct.o gen_gll.o gen_deriv.o CrvMatrixR.o CrvMatrixH.o CrvCompose.o NodeLoc.o ElemNodalDisp.o NodalRelRot.o NodalDataAt0.o NodalData.o ElasticForce.o ElementMatrix.o BeamStatic.o StaticSolution.o DriverStatic.o

Tilde.o: Tilde.f90
    gfortran -O3 -c Tilde.f90
    
OuterProduct.o: OuterProduct.f90
    gfortran -O3 -c OuterProduct.f90
    
gen_gll.o: gen_gll.f90
    gfortran -O3 -c gen_gll.f90
    
gen_deriv.o: gen_deriv.f90
    gfortran -O3 -c gen_deriv.f90

CrvMatrixR.o: CrvMatrixR.f90
    gfortran -O3 -c CrvMatrixR.f90
    
CrvMatrixH.o: CrvMatrixH.f90
    gfortran -O3 -c CrvMatrixH.f90
    
CrvCompose.o: CrvCompose.f90
    gfortran -O3 -c CrvCompose.f90
    
NodeLoc.o: NodeLoc.f90
    gfortran -O3 -c NodeLoc.f90
    
ElemNodalDisp.o: ElemNodalDisp.f90
    gfortran -O3 -c ElemNodalDisp.f90
    
NodalRelRot.o: NodalRelRot.f90
    gfortran -O3 -c NodalRelRot.f90
    
NodalDataAt0.o: NodalDataAt0.f90
    gfortran -O3 -c NodalDataAt0.f90
    
NodalData.o: NodalData.f90
    gfortran -O3 -c NodalData.f90
    
ElasticForce.o: ElasticForce.f90
    gfortran -O3 -c ElasticForce.f90
    
ElementMatrix.o: ElementMatrix.f90
    gfortran -O3 -c ElementMatrix.f90
    
BeamStatic.o: BeamStatic.f90
    gfortran -O3 -c BeamStatic.f90
    
StaticSolution.o: StaticSolution.f90
    gfortran -O3 -c StaticSolution.f90
    
DriverStatic.o: DriverStatic.f90
    gfortran -O3 -c DriverStatic.f90
