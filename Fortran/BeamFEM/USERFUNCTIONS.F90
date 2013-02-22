MODULE USERFUNCTIONS

CONTAINS


!==========================================================
! SOLVE THE LINEAR EQUATIONS Ax=b
! where A is a symmetric matrix
! Using conjugate gradient method
! Output: x
!----------------------------------------------------------
SUBROUTINE CG(A, b, N, X, iA, jA, Nz)
	IMPLICIT NONE
	
	INTEGER:: N, k, Nz
	INTEGER:: kmax
	INTEGER:: jA(Nz), iA(N+1)
	REAL(8):: A(Nz), b(N), x(N)
	REAL(8):: p(N), r(N), w(N)
	REAL(8):: norm0, norm1, beta, alpha, ErrLimit
!	REAL(8):: SparseMul(N)
	
!	write(*, *) " Begin CG "

	kmax = 100000
	ErrLimit = 1.0d-9

	! initial guess
	x = 0;
	x(1) = 1
	! r: resisual
	r = b - SparseMul(A, x, iA, jA, Nz, N)
	norm0 = DOT_PRODUCT(r, r)
	norm1 = norm0
	k = 0

	DO WHILE( (SQRT(norm1)/SQRT(norm0) > ErrLimit) .AND. (k<kmax) )
		IF (k==0) THEN
			p = r
		ELSE
			beta = norm1/norm0
			p = r + beta*p
			
		END IF
		w = SparseMul(A, p, iA, jA, Nz, N)
		alpha = norm1/DOT_PRODUCT(p, w)
		x = x + alpha*p
		r = r - alpha*w
		norm0 = norm1
		norm1 = DOT_PRODUCT(r, r)
		k = k+1
	END DO

	write(*, *) "k=",k
!	write(*, *) " end CG "

	If (k==kmax) Then
		write(*, *) " " 
		write(*, *) " The linear system doesn't converge!  "
		write(*, *) " "
	END IF

END SUBROUTINE CG
!==========================================================

!==========================================================
! Subroutine IndexArray(Arr, ArrIndex)
!
! Indexes an array arr, i.e., outputs the array index
! of length N such that arr(index(j)) is in ascending 
! order for j = 1, 2, . . . ,N. The input quantity arr
! is not changed.
!-----------------------------------------------------------
Subroutine IndexArray(Arr, ArrIndex, N)
IMPLICIT NONE

integer:: i, n, j, k
real(8):: arr(n), newarr(n)
real(8):: ta, temp
integer:: ArrIndex(n), ti

newarr = arr

DO i = 1,n
	ArrIndex(i) = i
END DO

DO i = 1, n
	temp = newarr(i)
	k = i
	DO j = i+1, n
		IF (temp > newarr(j)) THEN
			k = j
			temp = newarr(j)
		END IF
	END DO
	
	ta = newarr(k)
	newarr(k) = newarr(i)
	newarr(i) = ta

	ti = ArrIndex(k)
	ArrIndex(k) = ArrIndex(i)
	ArrIndex(i) = ti

END DO

!write(*, *) "finish subroutine index array!"

End Subroutine IndexArray
!==========================================================


!==========================================================
! Sparse matrix multiply vector A*x
!
! Output: SparseMul
!-----------------------------------------------------------
FUNCTION SparseMUL(A, x, iA, jA, Nz, N)
	IMPLICIT NONE

	INTEGER:: N, Nz
	REAL(8):: SparseMUL(N)
	REAL(8):: x(N), A(Nz)
	INTEGER:: jA(Nz), iA(N+1)
	INTEGER:: i, j, k1, k2

	SparseMUL = 0
	DO i = 1, N
		k1 = iA(i)
		k2 = iA(i+1)-1
		DO j = k1, k2
			SparseMUL(i) = SparseMUL(i) + A(j)*x(jA(j))
			IF (i /= jA(j)) THEN ! count for the symmetric part
				SparseMUL(jA(j)) = SparseMUL(jA(j)) + A(j)*x(i)
			END IF
		END DO

	END DO



END FUNCTION SparseMUL
!==========================================================

!==========================================================
! SOLVE THE LINEAR EQUATIONS Ax=b
! where A is a symmetric matrix
! Using conjugate gradient method
! Output: x
!----------------------------------------------------------
SUBROUTINE CG0(A, b, N, X)
	IMPLICIT NONE
	
	INTEGER:: N, k, Nz
	INTEGER:: kmax
	REAL(8):: A(N, N), b(N), x(N)
	REAL(8):: p(N), r(N), w(N)
	REAL(8):: norm0, norm1, beta, alpha, ErrLimit
!	REAL(8):: SparseMul(N)

	kmax = 100000
	ErrLimit = 0.0000000001

	! initial guess
	x = 0.00001;
	! r: resisual
	r = b - MATMUL(A, X)
	norm0 = DOT_PRODUCT(r, r)
	norm1 = norm0
	k = 0

	DO WHILE( (SQRT(norm1)/SQRT(norm0) > ErrLimit) .AND. (k<kmax) )
		IF (k==0) THEN
			p = r
		ELSE
			beta = norm1/norm0
			p = r + beta*p
			
		END IF
		w = MATMUL(A, p)
		alpha = norm1/DOT_PRODUCT(p, w)
		x = x + alpha*p
		r = r - alpha*w
		norm0 = norm1
		norm1 = DOT_PRODUCT(r, r)
		k = k+1
	END DO

	write(*, *) "k=",k
	If (k==kmax) Then
		write(*, *) " " 
		write(*, *) " The linear system doesn't converge!  "
		write(*, *) " "
	END IF

END SUBROUTINE CG0
!==========================================================


!==========================================================
! Sort node by ascending Y coordnate and then
! by Z coordinate
!
! - output: OrderNodes (order of nodes)
!----------------------------------------------------------
SUBROUTINE SortNodesYZ(OrderNodes, A, N)

	IMPLICIT NONE

	INTEGER:: i, j, N, k
	INTEGER:: out
	INTEGER:: OrderNodes(N)
	REAL(8):: A(N, 3), As(N, 3), temp(1,3)

   out = 65;
!	write(out, *) " Sort nodes according to Y then Z"
!	write(out, '(1x, I6, f15.7, f15.7)') (INT(A(i, 1)), A(i, 2), A(i, 3), i=1,N)
	write(out, *) " "
	
	As = A
	DO i = N, 2, -1
		k = i
		DO j = i-1, 1, -1
			IF ( As(j, 2) == As(k, 2) ) THEN
				IF ( As(j, 3) > As(k, 3)) THEN
					k = j
				END IF
			ELSE IF ( As(j, 2) > As(k, 2) ) THEN
					k = j
			END IF
		END DO

!		write(out, *) "i = ", i,  " k = ", k
		temp(1, :) = As(i, :)
		As(i, :) = As(k, :)
		As(k, :) = temp(1, :)
	END DO
	OrderNodes = INT(As(:, 1))

!	write(out, *) " After sorting "
!	write(out, '(1x, I8, f15.7, f15.7)') (OrderNodes(i), As(i,2), As(i,3), i=1,N)

END SUBROUTINE SortNodesYZ
!==========================================================


!==================================================================
! Calculate Eigenvalues and Eigenvectors
! 
!------------------------------------------------------------------
SUBROUTINE CalEigenValueVector(KK, MM, jKK, jMM, iKK, iMM, Nkz, Nmz, N, r, w, v)

   USE NWTC_Library
	
	IMPLICIT NONE

	integer, parameter:: oute = 65
	real(REKI), parameter:: err = 0.01, hmax = 10
	integer:: N, Nkz, Nmz, r
	real(REKI):: KK(Nkz), MM(Nmz)
	integer:: jKK(Nkz), jMM(Nmz)
	integer:: ikk(N+1), imm(N+1)

	real(8):: x(N, r), t(N), y(N, r), y1(N, r)
	integer:: i, j, tIndex(N), h, ti(N), FLAG
	real(REKI):: Kww(r, r), Mww(r, r)
	real(REKI):: w(r), v(N, r), phi(r, r), w_err(r), wold(r), w_err_max


	OPEN(oute, FILE = 'EigenOut.txt', STATUS='unknown', ACTION = 'write')

	!------ initialize matrix X -------
	x = 0
	x(:, 1) = 1

	! find maximum Mjj/Kjj
	DO i = 1, N
		t(i) = MM(imm(i))/KK(ikk(i)) 
	END DO
	
	CALL IndexArray(t, tIndex, N)
	CALL IndexArray(DBLE(tIndex), ti, N)

	DO i = 2, r
		x(ti(N-i+2), i) = 1
	END DO
	
	! calculate y=Mx
	Y = 0
	DO i = 1, r
		Y(:, i) = SparseMUL(DBLE(MM), x(:, i), iMM, jMM, NMz, N)
	END DO
	
	wold = 0
	DO h = 1, hmax

		! (1) Solve Kx = y
		Do i = 1, r
			CALL CG(DBLE(KK), Y(:, i), N, X(:, i), iKK, jKK, NKz)
		END DO

		! (2) calculate Kww
		Kww = MATMUL(Transpose(X), Y)

		! (3) calculate y1, y1 = M*x
		DO i = 1, r
			Y1(:, i) = SparseMUL(DBLE(MM), x(:, i), iMM, jMM, NMz, N)
		END DO

		! (4) calculate Mww=x^T*y1
		Mww = MATMUL(Transpose(X), Y1)

		! (5) find eigenvalues and eigenvectors for general eigen value
		!     problem Kww*phi = Mww*phi*w
		CALL JacobiEigen(DBLE(Kww), DBLE(Mww), r, DBLE(w), DBLE(phi), FLAG)
		IF (FLAG == 1) THEN
			v = MATMUL(x, phi)
			return
		END IF

		! (6) check iteration criteria
		DO i = 1, r
			w_err(i) = (w(i)-wold(i))/w(i)
		END DO
		
		! Find maximum w_err
		w_err_max = abs(w_err(1))
		Do i = 2, r-1
			w_err_max = Max(w_err_max, abs(w_err(i)) )
		END DO

	write(oute, *) " h = ", h
	write(oute, *) "w_err_max = ", w_err_max
	write(oute, *) " w " 
	write(oute, '(1x, e15.5)') (w(i), i = 1, r)


		IF (w_err_max <= err) THEN
			v = MATMUL(x, phi)
			return
		ELSE
			Y = MATMUL(Y1, phi)
			wold = w
		END IF

	END DO

	close(oute)

END SUBROUTINE CalEigenValueVector
!==================================================================

!==================================================================
! Jacobi Method for general eigenvalue problem 
!
!------------------------------------------------------------------
SUBROUTINE JacobiEigen(Kww, Mww, r, w, phi, FLAG)
	
	IMPLICIT NONE
	
	integer:: r
	REAL(8):: Kww(r, r), Mww(r, r), w(r), phi(r, r), phi0(r, r)
	REAL(8):: A(r, r), A1(r, r), err
	integer:: i, j, p, q, k, kmax, FLAG
	real(8):: QQ(r, r), maxapq, t, z, signz
	real(8):: cosphi, sinphi, cos2phi, sin2phi, ctgphi

	
	kmax = 100000
	err = 0.0000001
	! calculate Cholesky Factorization of Kww=Q^T*Q
	! and get the inverse of L

	CALL CholeskyAndInverse(Mww, r, QQ, FLAG)
	IF (FLAG == 1) THEN
		RETURN
	END IF

	A = MATMUL(MATMUL(QQ, Kww), TRANSPOSE(QQ))

!	A(1, 1:3) = (/3, 1, 1/)	
!	A(2, 1:3) = (/1, 2, 1/)
!	A(3, 1:3) = (/1, 1, 5/)


	phi0 = 0
	DO i = 1, r
		phi0(i, i) = 1
	END DO
	phi = phi0

	DO k = 1, kmax
		! find the maximum 
		maxapq = 0
		DO i = 1, r-1
			Do j = i+1, r
				IF (maxapq < ABS(A(i, j))) THEN
					maxapq = ABS(A(i, j))
					p = i
					q = j
				END IF
			END DO
		END DO
		
		z = 0.5*( A(p,p)-A(q,q) )/A(p, q)
		t = 1/z

		if (abs(t) < 1) THEN
			cos2phi = 1.0/sqrt(1+t*t)
			sin2phi = t/sqrt(1+t*t)
		ELSE
			cos2phi = abs(z)/sqrt(1+z*z)
			IF (z>=0) THEN
				signz = 1
			else
				signz = -1
			end if
			sin2phi = signz/sqrt(1+z*z)
		END IF

		cosphi = sqrt(0.5*(1+cos2phi))
		sinphi = 0.5*sin2phi/cosphi

		A1(p, p) = A(p, p)*cosphi*cosphi + A(q, q)*sinphi*sinphi + 2*A(p, q)*cosphi*sinphi
		A1(q, q) = A(p, p)*sinphi*sinphi + A(q, q)*cosphi*cosphi - 2*A(p, q)*cosphi*sinphi
		
		DO i = 1, r
			IF ((i/=p).and.(i/=q)) THEN
				A1(p, i) =  A(p, i)*cosphi + A(q, i)*sinphi
				A1(q, i) = -A(p, i)*sinphi + A(q, i)*cosphi
				A1(i, p) = A1(p, i)
				A1(i, q) = A1(q, i)
			END IF
		END DO 

		DO i = 1, r
			IF ((i/=p).and.(i/=q)) THEN
				DO j = 1, r
					IF ((j/=p).and.(j/=q)) THEN
						A1(i, j) = A(i, j)
						A1(j, i) = A1(i, j)
					END IF				
				END DO
			END IF
		END DO

		A1(p, q) = 0.5*(A(q, q) - A(p,p))*sin2phi+A(p, q)*cos2phi
		A1(q, p) = A1(p, q)

		DO i = 1, r
			phi(i, p) =  phi0(i, p)*cosphi + phi0(i, q)*sinphi
			phi(i, q) = -phi0(i, p)*sinphi + phi0(i, q)*cosphi
		END DO

		maxapq = 0
		DO i = 1, r-1
			Do j = i+1, r
				IF (maxapq < ABS(A1(i, j))) THEN
					maxapq = ABS(A1(i, j))
				END IF
			END DO
		END DO

		IF (maxapq < err) THEN
			DO i = 1, r
				w(i) = A1(i, i)
			END DO
			return
		ELSE
			A = A1
			phi0 = phi
		END IF

	 END DO

	
END SUBROUTINE JacobiEigen
!------------------------------------------------------------------


!==================================================================
! Cholesky factorization A = L*L^T
! and inverse of Q=inverse(L)
!
!------------------------------------------------------------------
SUBROUTINE CholeskyAndInverse(A, N, Q, FLAG)
	
	IMPLICIT NONE
	INTEGER:: N
	REAL(8):: A(N, N), Q(N, N), L(N, N), B(N), y(N)
!	REal(8):: MatrixTempA(N, N), MatrixTempL(N, N)

	INTEGER:: i, j, k, FLAG
	real(8):: temp

!	A = 0
!	A(1, 1) = 4; A(1,2)=-1;A(1,3)=1;
!	A(2, 1) = -1; A(2,2) = 4.25; A(2,3) = 2.75;
!	A(3, 1) =1; A(3, 2) = 2.75; A(3,3) = 3.5;
!	A(4,4)

	L = 0
	L(1, 1) = sqrt(A(1,1))

	DO i = 2, N
		L(i, 1) = A(i, 1)/L(1, 1)	
	END DO

	DO j = 2, N-1
		temp = 0
		DO k = 1, j-1
			temp = temp + L(j, k)*L(j, k) 
		END DO

		IF ( (A(j, j) - temp) < 0 ) THEN
			write(*, *) " sqrt < 0 !!! ", A(j, j), temp
			FLAG = 1
			RETURN
		ELSE
			L(j, j) = sqrt(A(j, j) - temp)
			FLAG = 0
		END IF
		
		DO i = j+1, N
			temp = 0
			Do k = 1, j-1
				temp = temp + L(i, k)*L(j, K)
			END DO
			L(i, j) = (A(i, j)-temp)/L(j,j)
		END DO
	END DO
	
	temp = 0
	DO k = 1, N-1
		temp = temp + L(N, k)*L(N, k)
	END DO
	IF ( (A(N, N) - temp) < 0 ) THEN
		write(*, *) " sqrt < 0 !!! ", A(N, N), temp
		FLAG = 1
		RETURN
	ELSE
		L(N, N) = sqrt(A(N, N) - temp)
		FLAG = 0
	END IF

	DO j = 1, N
		B = 0
		B(j) = 1
		DO i = 1, N
			temp = 0
			DO k = 1, i-1
				temp = temp + L(i, k)*y(k)
			END DO
			y(i) = (B(i)-temp)/L(i, i)
		END DO
		Q(:, j) = Y
	ENd DO

!	MatrixTempA = A - MATMUL(L, TRanspose(L))
!	MatrixTempL = MATMUL(Q, L)
!	write(*, *) MatrixTempA

END SUBROUTINE CholeskyAndInverse
!------------------------------------------------------------------



!==================================================================
! Test teh sparse matrix storage and CG method 
! Inverse Power method
!------------------------------------------------------------------
SUBROUTINE TESTSPARCECG

	
	IMPLICIT NONE

	integer, parameter:: oute = 65
	integer, parameter:: r = 3

	integer, parameter:: Nkz=9, Nmz=5, N=5

	real(8):: KK(Nkz), MM(Nmz)
	integer:: jKK(Nkz), jMM(Nmz)
	integer:: ikk(N+1), imm(N+1)
	integer:: i

	real(8):: beta(r), phi(N, r)

	KK =  (/2,-1,2,-1,2,-1,2,-1,2/)
	jKK = (/1,2,2,3,3,4,4,5,5/)
	iKK = (/1,3,5,7,9,10/)

	MM =  (/1.,1.,1.,1.,0.5/)
	jMM = (/1,2,3,4,5/)
	iMM = (/1,2,3,4,5,6/)
	
	OPEN(oute, FILE = 'Testsparcecg.txt', STATUS='unknown', ACTION = 'write')

	CALL LanczosEigen(KK, MM, jKK, jMM, iKK, iMM, Nkz, Nmz, N, r, beta, phi)

	write(oute, *) " "
	write(oute, '(e16.7)') beta
	write(oute, *) " "
	write(oute, '(3(e16.7))') phi
	
	close(oute)

END SUBROUTINE TESTSPARCECG
!==================================================================

!==================================================================
! get A(i, j)
!
!------------------------------------------------------------------
Function GETAIJ(A, iA, jA, N, nnz, i, j)

	IMPLICIT NONE
	integer:: N, nnz
	real(8):: A(nnz), GetAij
	integer:: jA(nnz), iA(N+1)
	integer:: i, j, k1, k2, jj

	k1 = iA(i)
	k2 = iA(i+1)-1

	DO jj = k1, k2
		IF (j == jA(jj)) THEN
			GetAij = A(jj)
			return;
		END IF
	END DO
		
	GetAij = 0


END function

!------------------------------------------------------------------

!==================================================================
! FUNCTION tqli (From Numerical recipe)
! finds eigenvalues and eigenvectors for tridiagonal matrix
!------------------------------------------------------------------
SUBROUTINE tqli(d,e,n,np,z)

	IMPLICIT NONE

	REAL(8), PARAMETER:: ONE = 1
	INTEGER n,np
	REAL(8):: d(np),e(np),z(np,np)

	INTEGER i,iter,k,l,m
	REAL(8):: b,c,dd,f,g,p,r,s

	do i=2,n 
		e(i-1)=e(i)
	end do

	e(n)=0.
	do l=1,n
		iter=0

1		do m=l,n-1
			dd=abs(d(m))+abs(d(m+1)) 
			if (abs(e(m))+dd.eq.dd) goto 2
		end do 

		m=n
2		if(m.ne.l)then
			if(iter.eq.30) pause 

			iter=iter+1
			g=(d(l+1)-d(l))/(2.*e(l))
			r=pythag(g,ONE)
			g=d(m)-d(l)+e(l)/(g+sign(r,g))
			s=1.
			c=1.
			p=0.

			do i=m-1,l,-1
				f=s*e(i)
				b=c*e(i)
				r=pythag(f,g)
				e(i+1)=r

				if(r.eq.0.) then
					d(i+1)=d(i+1)-p
					e(m)=0.
					goto 1
				end if
				s=f/r
				c=g/r
				g=d(i+1)-p
				r=(d(i)-g)*s+2.*c*b
				p=s*r
				d(i+1)=g+p
				g=c*r-b

				do k=1,n
					f=z(k,i+1)
					z(k,i+1)=s*z(k,i)+c*f
					z(k,i)=c*z(k,i)-s*f
				end do 
			end do 

			d(l)=d(l)-p
			e(l)=g
			e(m)=0.
			goto 1
		end if
	end do 

END SUBROUTINE TQLI

!==================================================================
! FUNCTION pythag (From Numerical recipe)
! finds sqrt(a**2 + b**2) without overflow or destructive underflow
!------------------------------------------------------------------
FUNCTION pythag(a, b)

	IMPLICIT NONE

	REAL(8):: a, b
	REAL(8):: fn_val
	REAL(8):: r, s, t, u, pythag

	fn_val = MAX(ABS(a), ABS(b))
	IF (fn_val /= 0) THEN
	  r = (MIN(ABS(a), ABS(b))/fn_val) ** 2
	  DO
		t = 4.0 + r
		IF (t /= 4.0) THEN
		  s = r / t
		  u = 1 + 2.0 * s
		  fn_val = u * fn_val
		  r = (s/u) ** 2 * r
		  CYCLE
		END IF
		EXIT
	  END DO
	END IF
	
	pythag = fn_val

END FUNCTION pythag
!------------------------------------------------------------------


!==================================================================
! Lanczos method to calculate eigenvector and eigenvalues
! 
!------------------------------------------------------------------
SUBROUTINE LanczosEigen(KK, MM, jKK, jMM, iKK, iMM, Nkz, Nmz, N, r, omega, phi)
	
!	USE ALLPARAMETERS, ONLY: tol_abs, tol_rel, itr_max, rr
	USE MGMRES

	IMPLICIT NONE

	integer, parameter:: oute = 77
	integer:: itr_max
	integer:: Nkz, Nmz, N, r, rr

	real(8):: KK(Nkz), MM(Nmz)
	integer:: jKK(Nkz), jMM(Nmz)
	integer:: ikk(N+1), imm(N+1), IndexOmega(r)
	real(8):: omega(r), phi(N, r)

	real(8):: x(N, r+1), y(N), xw(N), x_hat(N), x_hat_new(N)
	real(8):: alpha(r), beta(r+1), z(r, r), epsilon, err, temp, tol_rel, tol_abs
	integer:: i, k, kmax, j

!	OPEN(oute, FILE = 'TestEigenOut.txt', STATUS='unknown', ACTION = 'write')

   itr_max = 10000
   tol_abs = 1.0E-6
   tol_rel = 1.0E-6
     
	rr = N - 1
	kmax = 100
	err = 0.00001
	! initialize
	x = 0; z = 0;
	alpha = 0; beta = 0;
	x(:, 1) = 1
	
	! generate x1
	y = SparseMUL(MM, x(:, 1), iMM, jMM, Nmz, N)
	beta(1) = sqrt(Dot_product(x(:, 1), y))
	x(:, 1) = x(:, 1)/beta(1)

	! generate xi, i = 2, r
	DO i = 2, r+1
		y = SparseMUL(MM, x(:, (i-1)), iMM, jMM, Nmz, N)

		xw = 0
		xw(i-1) = 1
	    call mgmres_st (N, Nkz, ikk, jkk, KK, &
		& xw, y, itr_max, rr, tol_abs, tol_rel )

!		CALL CG(KK, y, N, xw, iKK, jKK, Nkz)
		
		alpha(i-1) = dot_product(xw, y)
		IF (i>2) THEN
			x_hat = xw - alpha(i-1)*x(:, i-1) - beta(i-1)*x(:, i-2)
		else
			x_hat = xw - alpha(i-1)*x(:, i-1)
		end if


		! re-orthogonalization
		x_hat_new = x_hat
		loop1: Do k = 1, kmax
			DO j = 1, i-1
				y = SparseMUL(MM, x(:, j), iMM, jMM, Nmz, N)
				epsilon = dot_product(x_hat_new, y)
				x_hat_new = x_hat_new - epsilon*x(:, j)
				if (abs(epsilon) <= err) then
					write(*, *) " epsilon = ", epsilon, " k = ", k
					exit loop1
				end if
			END DO
		end do loop1

		if (k == kmax) Then
			write(*, *) " "
			write(*, *) " re-orthognal fails !"
			write(*, *) " "
		end if

		y = SparseMUL(MM, x_hat_new, iMM, jMM, Nmz, N)

		temp = Dot_product(x_hat_new, y)
		if (temp >= 0) then
			beta(i) = sqrt(temp)
		else
			write(*, *)
			write(*, *) " error: sqrt < 0 ! ", temp
			write(*, *)
			return

		end if
		x(:, i) = x_hat/beta(i)

	END DO

	DO i = 1, r
		z(i, i) = 1
	end do
	
	CALL tqli(alpha, beta(1:r), r, r, z)

	DO i = 1, r
		omega(i) = 1/alpha(i)
	END DO
	
	phi(:, 1:r) = MATMUL(x(:, 1:r), z)


	CALL IndexArray(Omega, IndexOmega, R)
	Omega = Omega(IndexOmega)
	phi = phi(:, IndexOmega)

!	write(oute, *) " omega "
!	write(oute, '(e16.7)') omega
!	write(oute, *) " "
!	write(oute, *) " phi "
!	write(oute, '(3(e16.7))') phi

!	close(oute)


END SUBROUTINE LanczosEigen
!------------------------------------------------------------------

END MODULE USERFUNCTIONS