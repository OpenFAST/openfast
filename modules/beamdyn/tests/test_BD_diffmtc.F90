module test_BD_diffmtc

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer                    :: n, i
    type(BD_ParameterType)     :: parametertype
    real(BDKi), allocatable    :: test_shape(:,:), test_shapederivative(:,:)
    real(BDKi), allocatable    :: baseline_shape(:,:), baseline_shapederivative(:,:)
    real(BDKi), allocatable    :: gll_nodes(:)
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    character(1024)            :: testname
    real(BDKi)                 :: tolerance

    ! mathematica code for calculating shape function values and derivative values
    !ClearAll[p, \[Xi]j, xj, h, j, \[Xi], x, sol];
    !(* find p+1 GLL points *)
    !p = 4;
    !test[x_] = (1 - x^2) D[LegendreP[p, x], x];
    !\[Xi]j = Array[f, p + 1]; (* initialize list *)
    !For[j = 0, j <= p, j++,
    !sol = x /. FindRoot[test[x], {x, -Cos[(Pi (j))/p]}];
    ! \[Xi]j[[j + 1]] = sol
    !]
    !h[\[Xi]_, xj_] = -(((1. - \[Xi]^2) D[ LegendreP[p, \[Xi]], \[Xi]])/(
    !p (p + 1) LegendreP[p, xj] (\[Xi] - xj)));
    !hd[\[Xi]_, xj_] = D[h[\[Xi], xj], \[Xi]]
    !Plot[h[x, \[Xi]j[[1]]], {x, -1, 1}, PlotRange -> All]
    !h[0.5, \[Xi]j[[2]]] 

contains

    @test
    subroutine test_BD_diffmtc_2node()
        ! branches to test
        ! - 2 node, 1 element
        
        ! initialize NWTC_Num constants
        call SetConstants()
        
        tolerance = 1e-14
        
        ! --------------------------------------------------------------------------
        testname = "2-node element: evaluate shape/shapederivative at nodes"
        
        ! the shape functions should be:
        ! h1(-1) = 1,  h1(+1) = 0
        ! h2(-1) = 0,  h2(+1) = 1
        !
        ! this is satisfied by these linear equations
        ! h1(s) = 0.5*(1-s)
        ! h2(s) = 0.5*(1+s)
        ! therefore,
        ! h1` = -0.5
        ! h2` =  0.5
        !
        ! the expected result of BD_diffmtc is
        ! Shp - the shape function evaluated at the GLL nodes
        ! ShpDer - the shape function derivative evaluated at the GLL nodes
        !
        ! shp(1,:) = 1.0,  0.0
        ! shp(2,:) = 0.0,  1.0
        ! shpder(1,:) = -0.5, -0.5
        ! shpder(2,:) =  0.5,  0.5
        
        parametertype = simpleParameterType(1,2,2,0,1)
        n = parametertype%nodes_per_elem
        
        call AllocAry(test_shape, parametertype%nodes_per_elem, parametertype%nqp, "test_shape", ErrStat, ErrMsg)
        call AllocAry(test_shapederivative, parametertype%nodes_per_elem, parametertype%nqp, "test_shapederivative", ErrStat, ErrMsg)
        
        call AllocAry(parametertype%QPtN, parametertype%nodes_per_elem, 'QPtN', ErrStat, ErrMsg)
        parametertype%QPtN = (/ -1.0, 1.0 /)
        
        call AllocAry(gll_nodes, n, "GLL points array", ErrStat, ErrMsg)
        gll_nodes = (/ -1.0, 1.0 /)
        
        call BD_diffmtc(parametertype%nodes_per_elem, gll_nodes, parametertype%QPtN, parametertype%nqp, test_shape, test_shapederivative)
        
        call AllocAry(baseline_shape, parametertype%nodes_per_elem,parametertype%nqp, "baseline_shape", ErrStat, ErrMsg)
        call AllocAry(baseline_shapederivative, parametertype%nodes_per_elem,parametertype%nqp, "baseline_shapederivative", ErrStat, ErrMsg)
        baseline_shape(1,:) = (/ 1.0, 0.0 /)
        baseline_shape(2,:) = (/ 0.0, 1.0 /)
        baseline_shapederivative(1,:) = (/ -0.5, -0.5 /)
        baseline_shapederivative(2,:) = (/  0.5,  0.5 /)
        
        @assertEqual(baseline_shape, test_shape, tolerance, testname)
        @assertEqual(baseline_shapederivative, test_shapederivative, tolerance, testname)
       
        deallocate(test_shape) 
        deallocate(test_shapederivative) 
        deallocate(gll_nodes)
        deallocate(baseline_shape) 
        deallocate(baseline_shapederivative) 

        call BD_DestroyParam(parametertype, ErrStat, ErrMsg)
 
    end subroutine

    @test
    subroutine test_BD_diffmtc_5node()
        ! branches to test
        ! - 5 node, 1 element
        
        ! initialize NWTC_Num constants
        call SetConstants()
        
        tolerance = 1e-14
        
        ! --------------------------------------------------------------------------
        testname = "5-node element: evaluate shape/shapederivative at nodes and at three non-node locations"
        
        parametertype = simpleParameterType(1,5,8,0,1)
        !parametertype%nodes_per_elem = 5
        !parametertype%nqp = 8  ! testing the GLL nodes and three non-nodal locations (-0.8, 0.1, 0.4)
        n = parametertype%nodes_per_elem
        
        call AllocAry(test_shape, parametertype%nodes_per_elem, parametertype%nqp, "test_shape", ErrStat, ErrMsg)
        call AllocAry(test_shapederivative, parametertype%nodes_per_elem, parametertype%nqp, "test_shapederivative", ErrStat, ErrMsg)
        
        call AllocAry(parametertype%QPtN, parametertype%nodes_per_elem, 'QPtN', ErrStat, ErrMsg)
        ! in following, first five points are the 4th-order GLL locations; last three points are randomly chosen off-node points
        parametertype%QPtN = (/  -1.0, -0.6546536707079771, 0.0, 0.6546536707079771, 1.0, -0.8, 0.1, 0.4 /)
        
        call AllocAry(gll_nodes, n, "GLL points array", ErrStat, ErrMsg)
        gll_nodes = (/ -1.0, -0.6546536707079771, 0.0, 0.6546536707079771, 1.0 /)
        
        call BD_diffmtc(parametertype%nodes_per_elem, gll_nodes, parametertype%QPtN, parametertype%nqp, test_shape, test_shapederivative)
        
        call AllocAry(baseline_shape,           parametertype%nodes_per_elem, parametertype%nqp, "baseline_shape", ErrStat, ErrMsg)
        call AllocAry(baseline_shapederivative, parametertype%nodes_per_elem, parametertype%nqp, "baseline_shapederivative", ErrStat, ErrMsg)

        ! baseline_shape and baseline_shapederivative were calculated in Mathematica for fourth order Legendre spectral FE

        baseline_shape(1,:) = (/   1., 0., 0., 0., 0., 0.266400000000000,          0.0329625, 0.0564/)
        baseline_shape(2,:) = (/   0., 1., 0., 0., 0., 0.855336358376291,  -0.1121093731918499, -0.1746924181056723/)
        baseline_shape(3,:) = (/   0., 0., 1., 0., 0., -0.1776000000000001,              0.9669, 0.5263999999999998/)
        baseline_shape(4,:) = (/   0., 0., 0., 1., 0., 0.0854636416237095, 0.1525343731918499, 0.7234924181056726/)
        baseline_shape(5,:) = (/   0., 0., 0., 0., 1., -0.0296000000000000,           -0.0402875, -0.1316/)
        baseline_shapederivative(1,:) = (/-5., -1.240990253030983,0.375,-0.2590097469690172,0.5,-2.497,   0.27725, -0.1210000000000001/)
        baseline_shapederivative(2,:) = (/ 6.756502488724241,0.,-1.336584577695454,0.7637626158259736,-1.410164177942427,2.144324478146486,-0.896320373697923, 0.4156426862650312/)
        baseline_shapederivative(3,:) = (/-2.666666666666667,   1.74574312188794,0.,-1.74574312188794,2.666666666666667,0.5546666666666656,-0.6573333333333338, -2.069333333333333/)
        baseline_shapederivative(4,:) = (/ 1.410164177942427,-0.7637626158259736, 1.336584577695454,0.,-6.756502488724241,-0.31499114481315,1.696653707031257, 1.805690647068303/)
        baseline_shapederivative(5,:) = (/-0.5, 0.2590097469690172,-0.375,1.240990253030983,5.,0.1129999999999999,-0.42025, -0.03099999999999978/)
        
        @assertEqual(baseline_shape, test_shape, tolerance, testname)
        @assertEqual(baseline_shapederivative, test_shapederivative, tolerance, testname)

        deallocate(test_shape) 
        deallocate(test_shapederivative) 
        deallocate(gll_nodes)
        deallocate(baseline_shape) 
        deallocate(baseline_shapederivative) 

        call BD_DestroyParam(parametertype, ErrStat, ErrMsg)
        
    end subroutine
end module
