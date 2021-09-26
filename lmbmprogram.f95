MODULE r_precision ! Precision for reals.

IMPLICIT NONE

INTEGER, PARAMETER, PUBLIC :: prec = SELECTED_REAL_KIND(12)
! INTEGER, PARAMETER, PUBLIC :: prec = SELECTED_REAL_KIND(6)   ! do not work very well

END MODULE r_precision


MODULE param ! Parameters

USE r_precision, ONLY : prec     ! Precision for reals.

IMPLICIT NONE

!Parameters
INTEGER, PARAMETER, PUBLIC :: maxeps = 20, maxnrs = 2000
REAL(KIND=prec), PARAMETER, PUBLIC :: &
zero    = 0.0_prec,    & ! 
half    = 0.5_prec,    & ! 
one     = 1.0_prec,    & ! 
large   = 3.40282347*10.**38,  & ! HUGE(zero) ja seuraavalla rivilla TINY(zero)
small   = 1.17549435*10.**(-38)  ! Maaritetty erikseen

END MODULE param


MODULE initializat   ! Initialization of parameters.

USE r_precision, ONLY : prec   ! Precision for reals.
USE param, ONLY : large,small

IMPLICIT NONE

! Parameters
INTEGER, PARAMETER :: &
na     =    2, &        ! Size of the bundle na >= 2.
mcu    =    15, &       ! Upper limit for maximum number of stored corrections, mcu >= 3.
mcinit =    7           ! Initial maximum number of stored corrections, mcu >= mcinit >= 3. If mcinit <= 0, the default value mcinit = 3 will be used. However, the value mcinit = 7 is recommented.

! Real parameters (if parameter value <= 0.0 the default value of the parameter will be used).
REAL(KIND=prec), SAVE :: &
tolf  = 1.0E-4_prec, &  ! Tolerance for change of function values (default = 1.0E-8), pie 1.0E-7 reg, OLLUT 1.0e-4
tolf2 = -1.0_prec, & ! Second tolerance for change of function values. If tolf2 < 0 the the parameter and the corresponding termination criterion will be ignored. If tolf2 = 0 the default value 1.0E+4 will be used.
tolb  = -large + small, & ! Tolerance for the function value (default = -large).
tolg  = 1.0E-4_prec, &  ! Tolerance for the first termination criterion (default = 1.0E-6). kannattaa kokeilla viela pienempaa, pie 1.0E-7 reg OLLUT 1.0e-4
tolg2 = 1.0E-4_prec, &  ! Tolerance for the second termination criterion (default = tolg). tai sitten pienentaa tata, pie 1.0e-6 reg OLLUT 1.0e-4
eta   = 0.5_prec, &  ! Distance measure parameter, eta >= 0. If eta < 0  the default value 0.5 will be used. pie 0.5 reg
epsl  = 0.24_prec, &  ! Line search parameter, 0 < epsl < 0.25 (default = 1.0E-4). kannattaa ehka pienentaa, pie 0.24 reg
xmax  = 1000_prec     ! Maximum stepsize, 1 < XMAX (default = 1.5). voi vaikuttaa 2-1000, pie 1000 reg RATKAISEVA 2 svm pienille ongelmille

! Integer parameters (if parameter value <= 0.0 the default value of the parameter will be used).
INTEGER, SAVE :: &
mit    = 10000, &      ! Maximun number of iterations (default = 10000). 500 ihan hyva, OLLUT 10 000!
mfe    = 500000000, &  ! Maximun number of function evaluations (default = n*mit).
mtesf  = 10, &         ! Maximum number of iterations with changes of function values smaller than tolf (default = 10). OLLUT 50, NOPEUTTAMISEKSI 20
iscale = 0             ! Selection of the scaling: 0  - Scaling at every iteration with STU/UTU (default). 1  - Scaling at every iteration with STS/STU. 2  - Interval scaling with STU/UTU. 3  - Interval scaling with STS/STU. 4  - Preliminary scaling with STU/UTU. 5  - Preliminary scaling with STS/STU. 6  - No scaling.

END MODULE initializat


MODULE lmbm_sub  ! Subprograms for LMBM

USE r_precision, ONLY : prec   ! Precision for reals.

IMPLICIT NONE

CONTAINS

FUNCTION vdot(n,x,y) RESULT(xty)  ! Dot product of two vectors.

USE param, ONLY : zero

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
x,y         ! Input vectors.

! Scalar Arguments
INTEGER, INTENT(IN) :: &
n           ! Vectors dimension.

! Local scalars
REAL(KIND=prec) xty
INTEGER :: i

xty = zero

DO i = 1,n
xty = xty + x(i-1)*y(i-1)
END DO

END FUNCTION vdot

SUBROUTINE vneg(n,x,y)  ! Change the signs of vector elements.

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
x           ! Input vector.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(OUT) :: &
y           ! Output vector y:= -x.

! Scalar arguments
INTEGER, INTENT(IN) :: &
n           ! Vectors dimension.

! Local scalars
INTEGER :: i

!PM O(n)=n
DO i = 1,n
y(i-1) = -x(i-1)
END DO

END SUBROUTINE vneg

SUBROUTINE scalex(n,a,x,y)   ! Scaling a vector y:= a*x.

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
x           ! Input vector.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(OUT) :: &
y           ! Output vector y:= a*x.

! Scalar arguments
REAL(KIND=prec), INTENT(IN) :: &
a           ! Scaling parameter.
INTEGER, INTENT(IN) :: &
n           ! Vectors dimension.

! Local scalars
INTEGER :: i

DO i = 1,n
y(i-1) = a*x(i-1)
END DO

END SUBROUTINE scalex

SUBROUTINE xdiffy(n,x,y,z)   ! Difference of two vectors z:= x - y.

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
x,y         ! Input vector.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(OUT) :: &
z           ! Output vector z:= x - y.

! Scalar arguments
INTEGER, INTENT(IN) :: &
n           ! Vectors dimension.

! Local scalars
INTEGER :: i

DO  i = 1,n
z(i-1) = x(i-1) - y(i-1)
END DO

END SUBROUTINE xdiffy

SUBROUTINE xsumy(n,x,y,z)   ! Sum of two vectors z:= x + y.

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
x,y         ! Input vectors.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(OUT) :: &
z           ! Output vector z:= x + y.

! Scalar arguments
INTEGER, INTENT(IN) :: &
n           ! Vectors dimension.

! Local scalars
INTEGER :: i

DO  i = 1,n
z(i-1) = x(i-1) + y(i-1)
END DO

END SUBROUTINE xsumy

SUBROUTINE scdiff(n,a,x,y,z)   ! Difference of the scaled vector and a vector z:= a*x - y.

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
x,y           ! Input vector.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(OUT) :: &
z           ! Output vector z:= a*x - y.

! Scalar arguments
REAL(KIND=prec), INTENT(IN) :: &
a           ! Scaling factor.
INTEGER, INTENT(IN) :: &
n           ! Vectors dimension.

! Local scalars
INTEGER :: i

DO  i = 1,n
z(i-1) = a*x(i-1) - y(i-1)
END DO

END SUBROUTINE scdiff

SUBROUTINE scsum(n,a,x,y,z)   ! Sum of a vector and the scaled vector z:= y + a*x.

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
x,y         ! Input vectors.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(OUT) :: &
z           ! Output vector z:= a*x + y.

! Scalar arguments
REAL(KIND=prec), INTENT(IN) :: &
a           ! Scaling factor.
INTEGER, INTENT(IN) :: &
n           ! Vectors dimension.

! Local scalars
INTEGER :: i

DO  i = 1,n
z(i-1) = a*x(i-1) + y(i-1)
END DO

END SUBROUTINE scsum

SUBROUTINE copy(n,x,y)  ! Copying a vector y:= x.

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
x           ! Input vector.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(OUT) :: &
y           ! Output vector.

! Scalar arguments
INTEGER, INTENT(IN) :: &
n           ! Vectors dimension.

! Local scalars
INTEGER :: i

DO i = 1,n
y(i-1) = x(i-1)
END DO

END SUBROUTINE copy

SUBROUTINE copy2(n,x,y,z,v)  ! Copying of two vectors: y:=x, v:=z.

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
x,z         ! Input vector.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(OUT) :: &
y,v         ! Output vectors.

! Scalar arguments
INTEGER, INTENT(IN) :: &
n           ! Vectors dimension.

! Local scalars
INTEGER :: i

DO i = 1,n
y(i-1) = x(i-1)
v(i-1) = z(i-1)
END DO

END SUBROUTINE copy2

SUBROUTINE vxdiag(n,d,x,y)   ! Vector is multiplied by a diagonal matrix y:=d*x.

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
x, &        ! Input vector.
d           ! Diagonal matrix stored as a vector with n elements.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(OUT) :: &
y           ! Output vector y:= d*x.

! Scalar arguments
INTEGER, INTENT(IN) :: &
n           ! Vectors dimension.

! Local scalars
INTEGER :: i

DO  i = 1,n
y(i-1) = x(i-1)*d(i-1)
END DO

END SUBROUTINE vxdiag

SUBROUTINE symax(n,m,iold,a,x,y)  ! Multiplication of a dense symmetric matrix A by a vector x.

USE param, ONLY : zero

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:m-1), INTENT(IN) :: &
x           ! Input vector stored in a circular order.
REAL(KIND=prec), DIMENSION(0:(n*(n+1)/2)-1), INTENT(IN) :: &
a           ! Dense symmetric matrix stored in the packed form: a(n*(n+1)/2).
REAL(KIND=prec), DIMENSION(0:m-1), INTENT(OUT) :: &
y           ! Output vector y:= a*x. Vector y has the same circular order than x.

! Scalar arguments
INTEGER, INTENT(IN) :: &
n, &         ! Order of matrix A.
m, &         ! Length of vector x, m >= n, note that only n components from vector x are used.
iold         ! Index, which controlls the circular order of the vector x.

! Local scalars
INTEGER :: i,j,k,l

DO j=1,n
l=j+iold-1
IF (l > m) l=l-m
y(l-1) = zero
k=l
DO i=j,n
y(l-1) = a((i-1)*i/2+j-1)*x(k-1)+y(l-1)
k=k+1
IF (k > m) k=k-m
END DO
END DO

DO j=2,n
l=j+iold-1
IF (l > m) l=l-m
k=iold
DO i=1,j-1
IF (k > m) k=k-m
y(l-1) = a((j-1)*j/2+i-1)*x(k-1)+y(l-1)
k=k+1
END DO
END DO

END SUBROUTINE symax

SUBROUTINE cwmaxv(n,m,a,x,y)  ! Multiplication of a columnwise stored dense rectangular matrix A by a vector x.

USE param, ONLY : zero

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n*m-1), INTENT(IN) :: &
a           ! Rectangular matrix stored columnwise in the one-dimensional array (dimension n*m).
REAL(KIND=prec), DIMENSION(0:m-1), INTENT(IN) :: &
x           ! Input vector (dimension m).
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(OUT) :: &
y           ! Output vector equal to s*a*x. If m = 0 y is a zero vector. 

! Scalar arguments
INTEGER, INTENT(IN) :: &
n, &        ! Number of rows of the matrix A.
m           ! Number of columns of the matrix A.

! Local scalars
INTEGER :: i,j,k

DO i = 1,n
y(i-1) = zero
END DO

k = 1

DO j = 1,m
CALL scsum(n,x(j-1),a(k-1:),y,y)
k = k + n
END DO

END SUBROUTINE cwmaxv

SUBROUTINE rwaxv2(n,m,a,b,x,y,v,w)  ! Multiplication of two rowwise stored dense rectangular matrices A and B by vectors x and y.

USE param, ONLY : zero

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
x,y         ! Input vectors (dimension n).
REAL(KIND=prec), DIMENSION(0:n*m-1), INTENT(IN) :: &
a,b         ! Rectangular matrices stored rowwise in the one-dimensional array (dimension n*m).
REAL(KIND=prec), DIMENSION(0:m-1), INTENT(OUT) :: &
v,w         ! Output vectors v=a*x and w=b*y. 

! Scalar arguments
INTEGER, INTENT(IN) :: &
n, &        ! Number of columns of the matrices A and B.
m           ! Number of rows of the matrices A and B.

! Local scalars
REAL(KIND=prec) :: tmp1,tmp2
INTEGER :: i,j,k

k = 0

DO i = 1,m
tmp1 = zero
tmp2 = zero
DO j = 1,n
tmp1 = tmp1 + a(k+j-1)*x(j-1)
tmp2 = tmp2 + b(k+j-1)*y(j-1)
END DO
v(i-1) = tmp1
w(i-1) = tmp2
k = k + n
END DO

END SUBROUTINE rwaxv2

SUBROUTINE trlieq(n,m,iold,u,x,y,job,ierr)  ! Solving x from linear equation u*x=y or u'*x=y, where u is an upper triangular matrix.

USE param, ONLY : small

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:m-1), INTENT(IN) :: &
y           ! Input vector stored in a circular order.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
u           ! Triangular matrix.
REAL(KIND=prec), DIMENSION(0:m-1), INTENT(OUT) :: &
x           ! Output vector y:= a*x. Vector y has the same circular order than x. Note that x may be equal to y in calling sequence.

! Scalar arguments
INTEGER, INTENT(IN) :: &
n, &         ! Order of matrix U.
m, &         ! Length of vectors x and y, m >= n, note that only n components from vectors are used.
iold, &      ! Index, which controlls the circular order of the vectors x and y.
job          ! Option: 0  - x:=(u')**(-1)*y, u upper triangular. 1  - x:=u**(-1)*y, u upper triangular.
INTEGER, INTENT(OUT) :: &
ierr         ! Error indicador: 0   - Everything is ok. -3   - Error; 0 at diagonal.

! Local scalars
INTEGER :: i,ii,ij,j,k,l,ji

! Intrinsic functions
INTRINSIC ABS

ierr = -3

DO i=1,m
x(i-1)=y(i-1)
END DO

IF (job == 0) THEN

! x=u'**(-1)*y, u' = [u1         ] is lower triangular.
!                    [u2 u3      ]
!                    [u4 u5 u6   ]
!                    [.  .  .  . ]

ii = 0

DO  i = 1,n
ii=ii+i
l=i+iold-1
IF (l > m) l=l-m
IF (ABS(u(ii-1)) <= small) RETURN
x(l-1) = x(l-1)/u(ii-1)
DO j = i+1,n
ji = (j-1)*j/2+i
k=j+iold-1
IF (k > m) k=k-m
x(k-1) = x(k-1) - u(ji-1)*x(l-1)
END DO
END DO

ELSE IF (job == 1) THEN

! x=u**(-1)*y, u = [u1 u2 u4 . ] is upper triangular.
!                  [   u3 u5 . ]
!                  [      u6 . ]
!                  [         . ]

ii = n* (n+1)/2

DO i = n,1,-1
l=i+iold-1
IF (l > m) l=l-m
IF (ABS(u(ii-1)) <= small) RETURN
ij = ii
DO j = i + 1,n
k=j+iold-1
IF (k > m) k=k-m
ij = ij + j - 1
x(l-1) = x(l-1) - u(ij-1)*x(k-1)
END DO
x(l-1)=x(l-1)/u(ii-1)
ii = ii - i
END DO
ELSE
RETURN
END IF

ierr = 0

END SUBROUTINE trlieq

SUBROUTINE lineq(n,m,iold,a,x,y,ierr)  ! Solving x from linear equation A*x=y. Positive definite matrix A+E is given using the factorization A+E=L*D*L' obtained by the subroutine mxdpgf.

USE param, ONLY : small

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:m-1), INTENT(IN) :: &
y           ! Input vector stored in a circular order (dimension m).
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
a           ! Factorization a+e=l*d*l' obtained by the subroutine mxdpgf.
REAL(KIND=prec), DIMENSION(0:m-1), INTENT(OUT) :: &
x           ! Output vector y:= a*x. Vector x has the same circular order than y. Note that x may be equal to y in calling sequence.

! Scalar arguments
INTEGER, INTENT(IN) :: &
n, &        ! Order of matrix a.
m, &        ! Length of vectors x and y, m >= n, note that only n components from vectors are used.
iold        ! Index, which controlls the circular order of the vectors x and y.
INTEGER, INTENT(OUT) :: &
ierr        ! Error indicador: 0   - Everything is ok. -2   - Error; indefinite matrix.

! Local scalars
INTEGER :: i,ii,ij,j,k,l

ierr = -2

! Phase 1: x=l**(-1)*x

ij = 0
DO i = 1,n
l=i+iold-1
IF (l > m) l=l-m
x(l-1) = y(l-1)
DO j = 1,i - 1
ij = ij + 1
k=j+iold-1
IF (k > m) k=k-m
x(l-1) = x(l-1) - a(ij-1)*x(k-1)
END DO
ij = ij + 1
END DO

! Phase 2: x:=d**(-1)*x

ii = 0
DO i = 1,n
ii = ii + i
IF (a(ii-1) <= small) RETURN
l=i+iold-1
IF (l > m) l=l-m
x(l-1) = x(l-1)/a(ii-1)
END DO

! Phase 3: x:=trans(l)**(-1)*x

ii = n* (n-1)/2
DO i = n - 1,1,-1
ij = ii
l=i+iold-1
IF (l > m) l=l-m
DO j = i + 1,n
k=j+iold-1
IF (k > m) k=k-m
ij = ij + j - 1
x(l-1) = x(l-1) - a(ij-1)*x(k-1)
END DO
ii = ii - i
END DO

ierr = 0

END SUBROUTINE lineq

SUBROUTINE mxdpgf(n,a,inf,alf,tau)  ! Factorization A+E=L*D*trans(L) of a dense symmetric positive definite matrix A+E, where D and E are diagonal positive definite matrices and L is a lower triangular matrix. If A is sufficiently positive definite then E=0.


USE param, ONLY : zero,one

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(INOUT) :: &
a         ! On input: Dense symmetric matrix stored in the packed form. On output: factorization A+E=L*D*trans(L).

! Scalar arguments
REAL(KIND=prec), INTENT(INOUT) :: &
alf       ! On input a desired tolerance for positive definiteness. On output the most negative diagonal element used in the factorization process (if inf>0).
REAL(KIND=prec), INTENT(OUT) :: &
tau       ! Maximum diagonal element of matrix E.
INTEGER, INTENT(IN) :: &
n         ! Order of matrix a.
INTEGER, INTENT(OUT) :: &
inf       ! An information obtained in the factorization process: inf=0  - A is sufficiently positive definite and E=0. inf<0  - A is not sufficiently positive definite and E>0. inf>0  - A is indefinite and inf is an index of the most negative diagonal element used in the factorization process.

! Local scalars
REAL(KIND=prec) :: bet,del,gam,rho,sig,tol
INTEGER :: i,ij,ik,j,k,kj,kk,l

! Intrinsic functions
INTRINSIC ABS,MAX

l = 0
inf = 0
tol = alf

! Estimation of the matrix norm

alf = zero
bet = zero
gam = zero
tau = zero
kk = 0

DO k = 1,n
kk = kk + k
bet = MAX(bet,ABS(a(kk-1)))
kj = kk
DO j = k + 1,n
kj = kj + j - 1
gam = MAX(gam,ABS(a(kj-1)))
END DO
END DO

bet = MAX(tol,bet,gam/n)
del = tol*MAX(bet,one)
kk = 0

DO k = 1,n
kk = kk + k
! Determination of a diagonal correction
sig = a(kk-1)
IF (alf > sig) THEN
alf = sig
l = k
END IF
gam = zero
kj = kk
DO j = k + 1,n
kj = kj + j - 1
gam = MAX(gam,ABS(a(kj-1)))
END DO
gam = gam*gam
rho = MAX(ABS(sig),gam/bet,del)
IF (tau < rho-sig) THEN
tau = rho - sig
inf = -1
END IF
! Gaussian elimination
a(kk-1) = rho
kj = kk
DO j = k + 1,n
kj = kj + j - 1
gam = a(kj-1)
a(kj-1) = gam/rho
ik = kk
ij = kj
DO i = k + 1,j
ik = ik + i - 1
ij = ij + 1
a(ij-1) = a(ij-1) - a(ik-1)*gam
END DO
END DO
END DO

IF (l > 0 .AND. ABS(alf) > del) inf = l

END SUBROUTINE mxdpgf

SUBROUTINE calq(n,m,iold,a,x,y)  ! Solving x from linear equation A*x=y.

USE param, ONLY : small

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:m-1), INTENT(IN) :: &
y           ! Input vector stored in a circular order (dimension m).
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(INOUT) :: &
a           ! On input: Dense symmetric matrix stored in the packed form. On output: factorization A+E=L*D*trans(L).
REAL(KIND=prec), DIMENSION(0:m-1), INTENT(OUT) :: &
x           ! Output vector y:= a*x. Vector x has the same circular order than y. Note that x may be equal to y in calling sequence.

! Scalar arguments
INTEGER, INTENT(IN) :: &
n, &        ! Order of matrix a.
m, &        ! Length of vectors x and y, m >= n, note that only n components from vectors are used.
iold        ! Index, which controlls the circular order of the vectors x and y.
! INTEGER, INTENT(OUT) :: &
! ierr        ! Error indicador: 0   - Everything is ok. -2   - Error; indefinite matrix.

! Local scalars
REAL(KIND=prec) :: eta,bet
INTEGER :: inf,ierr

eta = small+small

CALL mxdpgf(n,a,inf,eta,bet)
CALL lineq(n,m,iold,a,x,y,ierr)

END SUBROUTINE calq

END MODULE lmbm_sub


MODULE lmbm_mod  ! Limited memory bundle method

USE r_precision, ONLY : prec  ! Precision for reals.

IMPLICIT NONE

CONTAINS

!***********************************************************************
!*                                                                     *
!*     * SUBROUTINE init_lmbm *                                        *
!*                                                                     *
!*     Initialization for limited memory bundle subroutine for         *
!*     large-scale unconstrained nonsmooth optimization.               *
!*                                                                     *
!***********************************************************************

SUBROUTINE init_lmbm(n,mc,iterm)

USE param, ONLY : zero,half,small,large
USE initializat, ONLY: na,mcu,iscale,tolf,tolf2,tolb,tolg,tolg2,xmax,eta,epsl,mtesf,mit,mfe

IMPLICIT NONE

! Scalar arguments
INTEGER, INTENT(IN) :: n         ! Number of the variables.
INTEGER, INTENT(INOUT) :: mc     ! Initial maximum number of stored corrections.
INTEGER, INTENT(OUT) :: iterm    ! Cause of termination: 0  - Everything is ok. -5  - Invalid input parameters.

! Initialization and error checking

! PM muokkaus
iterm = 0

IF (n <= 0) THEN
iterm = -5
RETURN 
END IF

IF (na < 2) THEN
iterm = -5
RETURN
END IF

IF (epsl >= 0.25_prec) THEN
iterm = -5
RETURN
END IF

IF (mcu <= 3) THEN
iterm = -5
RETURN
END IF

IF (mc > mcu) THEN
mc = mcu
END IF

! Default values

IF (mc    <= 0) mc       = 3                   ! Initial maximum number of corrections.
IF (mit   <= 0) mit      = 10000               ! Maximum number of iterations.
IF (mfe   <= 0) mfe      = n*mit               ! Maximum number of function evalutions.
IF (tolf  <= zero) tolf  = 1.0E-08_prec        ! Tolerance for change of function values.
IF (tolf2 == zero) tolf2 = 1.0E+04_prec        ! Second tolerance for change of function values.
IF (tolb  == zero) tolb  = -large + small      ! Tolerance for the function value.
IF (tolg  <= zero) tolg  = 1.0E-06_prec        ! Tolerance for the first termination criterion.
IF (tolg2 <= zero) tolg2 = tolg                ! Tolerance for the second termination criterion.
IF (xmax  <= zero) xmax  = 1.5_prec            ! Maximum stepsize.
IF (eta   <  zero) eta   = half                ! Distance measure parameter.
IF (epsl  <= zero) epsl  = 1.0E-04_prec        ! Line search parameter.
IF (mtesf <= 0) mtesf    = 10                  ! Maximum number of iterations with changes of function values smaller than tolf.
IF (iscale > 6 .OR. iscale < 0) iscale = 0     ! Selection of the scaling.

END SUBROUTINE init_lmbm

SUBROUTINE restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk,alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)   ! Initialization and reinitialization.

USE param, ONLY : zero,one
USE lmbm_sub, ONLY : copy,vneg

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
gp        ! Basic subgradient of the objective function.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(INOUT) :: &
g         ! Current (auxiliary) subgradient of the objective function.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(OUT) :: &
d         ! Search direction.

! Scalar arguments
INTEGER, INTENT(IN) :: &
n, &      ! Number of variables.
mcinit    ! Initial maximum number of stored corrections.
INTEGER, INTENT(OUT) :: &
mc, &     ! Current maximum number of stored corrections.
mcc, &    ! Current number of stored corrections.
inew, &   ! Index for the circular arrays.
ibun, &   ! Index for the circular arrays in bundle updating.
ibfgs, &  ! Index of the type of BFGS update.
nnk, &    ! Consecutive null steps counter.
ic, &     ! Correction indicator.
icn, &    ! Correction indicator for null steps.
mal, &    ! Current size of the bundle.
iflag     ! Index for adaptive version.
INTEGER, INTENT(INOUT) :: &
iters, &  ! Null step indicator. 0  - Null step. 1  - Serious step.
ncres     ! Number of restarts.
REAL(KIND=prec), INTENT(OUT) :: & 
alfn, &   ! Locality measure.
alfv, &   ! Aggregate locality measure.
gamma     ! Scaling parameter.

! Restart
mc    = mcinit
mcc   = 0
inew  = 1
ibun  = 1
ibfgs = 0
ic    = 0
icn   = 0
mal   = 0
ncres = ncres + 1
iflag = 0

IF (iters == 0) THEN
CALL copy(n,gp,g)
iters = 1
nnk = 0
alfv=zero
alfn=zero
END IF

gamma = one
CALL vneg(n,g,d)

END SUBROUTINE restar

SUBROUTINE dobun(n,ma,mal,x,g,fu,ay,ag,af,iters,ibun)   ! Bundle construction.

USE lmbm_sub, ONLY : copy2

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
g, &      ! Subgradient of the objective function.
x         ! Vector of variables
REAL(KIND=prec), DIMENSION(0:(n*ma)-1), INTENT(INOUT) :: &
ay, &     ! Matrix whose columns are bundle points (stored in one-dimensional n*ma -array).
ag        ! Matrix whose columns are bundle subgradients (stored in one-dimensional n*ma -array).
REAL(KIND=prec), DIMENSION(0:ma-1), INTENT(INOUT) :: &
af        ! Vector of values of bundle functions (stored in one-dimensional ma -array).

! Scalar arguments
INTEGER, INTENT(IN) :: &
n, &      ! Number of variables.
iters, &  ! Null step indicator. 0  - Null step. 1  - Serious step.
ma        ! Maximum size of the bundle.
INTEGER, INTENT(INOUT) :: &
ibun, &   ! Index for the circular arrays in bundle updating.
mal       ! Current size of the bundle.
REAL(KIND=prec), INTENT(IN) :: & 
fu        ! Value of the objective function.

! Local scalars
INTEGER :: i,j

IF (iters == 1) THEN

! Serious step

af(ibun-1) = fu
i = (ibun-1)*n+1
CALL copy2(n,g,ag(i-1:),x,ay(i-1:))

ELSE

! Null step

IF (mal < ma) THEN
af(ibun-1) = af(mal-1)
af(mal-1) = fu
i = mal*n + 1
CALL copy2(n,ag(i-n-1:),ag(i-1:),ay(i-n-1:),ay(i-1:))
CALL copy2(n,g,ag(i-n-1:),x,ay(i-n-1:))
ELSE
i = ibun-1
IF (i < 1) i = mal
af(ibun-1) = af(i-1)
af(i-1) = fu
i = (ibun-2)*n + 1
IF (i < 1) i = (mal-1)*n + 1
j = (ibun-1)*n + 1
CALL copy2(n,ag(i-1:),ag(j-1:),ay(i-1:),ay(j-1:))
CALL copy2(n,g,ag(i-1:),x,ay(i-1:))
END IF
END IF

mal = mal + 1
IF (mal > ma) mal = ma
ibun = ibun + 1
IF (ibun > ma) ibun = 1

END SUBROUTINE dobun

SUBROUTINE destep(n,ma,mal,x,af,ag,ay,ibun,d,fu,df,t,eta,iterm)   ! Stepsize selection.

USE param, ONLY : zero,half,one,large

IMPLICIT NONE

! Scalar arguments
REAL(KIND=prec), INTENT(INOUT) :: & 
t         ! Initial stepsize
REAL(KIND=prec), INTENT(IN) :: & 
df, &     ! Directional derivative.
fu, &     ! Value of the objective function.
eta       ! Distance measure parameter.
INTEGER, INTENT(IN) :: & 
n, &      ! Number of variables
ma, &     ! Maximum size of the bundle.
mal, &    ! Current size of the bundle.
ibun      ! Index for the circular arrays in bundle updating.
INTEGER, INTENT(OUT) :: &
iterm     ! Cause of termination: 0  - Everything is ok. -6  - Error.

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
x, &      ! Vector of variables (n array).
d         ! Direction vector (n array).
REAL(KIND=prec), DIMENSION(0:(n*ma)-1), INTENT(IN) :: &
ay, &     ! Matrix whose columns are bundle points (stored in one-dimensional n*ma -array).
ag        ! Matrix whose columns are bundle subgradients (stored in one-dimensional n*ma -array).
REAL(KIND=prec), DIMENSION(0:ma-1), INTENT(IN) :: &
af        ! Vector of values of bundle functions (stored in one-dimensional ma -array).

! Local arrays
REAL(KIND=prec), DIMENSION(0:(2*ma)-1) :: &
tmparray  ! Auxiliary array.

! Local scalars
REAL(KIND=prec) :: alf,alfl,alfr,bet,betl,betr,dx,q,r,w
INTEGER :: i,j,jn,k,l,lq,ib

! Intrinsic functions
INTRINSIC ABS,REAL,MAX,MIN,SQRT

iterm = 0
alfl = zero
betl = zero
w = df*t* (one - t*half)

! Initial choice of possibly active lines

k = 0
l = -1
jn = (ibun-1)*n
betr = - large

DO j=1,mal-1
ib = ibun - 1 + j
IF (ib > mal) ib = ib - mal
IF (jn >= mal*n) jn = jn - mal*n
r = zero
bet = zero
alfl = af(ib-1) - fu
DO i=1,n
dx = x(i-1) - ay(jn+i-1)
q = ag(jn+i-1)
r = r + dx*dx
alfl = alfl + dx*q
bet = bet + d(i-1)*q
END DO
alf = MAX(ABS(alfl),eta*r)
r = one - bet/df
IF (r*r + (alf+alf)/df > 1.0E-6_prec) THEN
k = k + 1
tmparray(k-1) = alf
tmparray(ma+k-1) = bet
r = t*bet - alf
IF (r > w) THEN
w = r
l = k
END IF
END IF
betr = MAX(betr,bet-alf)
jn = jn + n
END DO

lq = -1
IF (betr <= df*half) RETURN
lq = 1
IF (l <= 0) RETURN
betr = tmparray(ma+l-1)

IF (betr <= zero) THEN
IF (t < one .OR. betr == zero) RETURN
lq = 2
END IF

alfr = tmparray(l-1)

! Iteration loop

ds_iteration: DO

IF (lq >= 1) THEN
q = one - betr/df
r = q + SQRT(q*q + (alfr+alfr)/df)
IF (betr >= zero) r = - (alfr+alfr)/ (df*r)
r = MIN(1.95_prec,MAX(zero,r))
ELSE
IF (ABS(betr-betl)+ABS(alfr-alfl)< -1.0E-4_prec*df) RETURN
IF (betr-betl  == zero) THEN
iterm = -6
RETURN
END IF
r = (alfr-alfl)/ (betr-betl)
END IF

IF (ABS(t-r) < 1.0E-4_prec) RETURN
t = r
tmparray(l-1) = - one
w = t*betr - alfr
l = -1

DO j = 1,k
alf = tmparray(j-1)
IF (alf < zero) EXIT
bet = tmparray(ma+j-1)
r = t*bet - alf
IF (r > w) THEN
w = r
l = j
END IF
END DO

IF (l < 0) RETURN
bet = tmparray(ma+l-1)
IF (bet == zero) RETURN

! New interval selection

alf = tmparray(l-1)

IF (bet < zero) THEN
IF (lq == 2) THEN
alfr = alf
betr = bet
ELSE
alfl = alf
betl = bet
lq = 0
END IF
ELSE
IF (lq == 2) THEN
alfl = alfr
betl = betr
lq = 0
END IF
alfr = alf
betr = bet
END IF

END DO ds_iteration

END SUBROUTINE destep

SUBROUTINE nulstep(n,ma,mal,x,af,ag,ay,ibun,d,fu,df,t,eta,iterm)   ! Stepsize selection.

USE param, ONLY : zero,one,large

IMPLICIT NONE

! Scalar arguments
REAL(KIND=prec), INTENT(INOUT) :: & 
t         ! Initial stepsize.
REAL(KIND=prec), INTENT(IN) :: & 
df, &     ! Directional derivative.
fu, &     ! Value of the objective function.
eta       ! Distance measure parameter.
INTEGER, INTENT(IN) :: &
n, &      ! Number of variables.
ma, &     ! Maximum size of the bundle.
mal, &    ! Current size of the bundle.
ibun      ! Index for the circular arrays in bundle updating.
INTEGER, INTENT(OUT) :: &
iterm     ! Cause of termination: 0  - Everything is ok. -6  - Error.

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
x, &      ! Vector of variables (n array).
d         ! Direction vector (n array).
REAL(KIND=prec), DIMENSION(0:(n*ma)-1), INTENT(IN) :: &
ay, &     ! Matrix whose columns are bundle points (stored in one-dimensional n*ma -array).
ag        ! Matrix whose columns are bundle subgradients (stored in one-dimensional n*ma -array).
REAL(KIND=prec), DIMENSION(0:(4*ma)-1), INTENT(IN) :: &
af        ! Vector of values of bundle functions (stored in one-dimensional 4*ma -array).

! Local arrays
REAL(KIND=prec), DIMENSION(0:(2*ma)-1) :: &
tmparray  ! Auxiliary array.

! Local scalars
REAL(KIND=prec) :: alf,alfl,alfr,bet,betl,betr,dx,q,r,w
INTEGER :: i,j,jn,k,l,ib

! Intrinsic functions
INTRINSIC ABS,REAL,MAX,MIN,SQRT

! Initial choice of possibly active parabolas

iterm = 0
w = df*t
k = 0
l = -1
jn = (ibun-1)*n
betr = - large

DO j = 1,mal - 1
ib = ibun - 1 + j
IF (ib > mal) ib = ib - mal
IF (jn >= mal*n) jn = jn - mal*n
bet = zero
r = zero
alfl = af(ib-1) - fu
DO i = 1,n
dx = x(i-1) - ay(jn+i-1)
r = r + dx*dx
q = ag(jn+i-1)
alfl = alfl + dx*q
bet = bet + d(i-1)*q
END DO
alf = MAX(ABS(alfl),eta*r)
betr = MAX(betr,bet-alf)
IF (alf < bet-df) THEN
k = k + 1
r = t*bet - alf
tmparray(k-1) = alf
tmparray(ma+k-1) = bet
IF (r > w) THEN
w = r
l = k
END IF
END IF
jn = jn + n
END DO

IF (l <= 0) RETURN
betr = tmparray(ma+l-1)
alfr = tmparray(l-1)
alf = alfr
bet = betr
alfl = zero
betl = df

! Iteration loop

ns_iteration: DO

w = bet/df
IF (ABS(betr-betl)+ABS(alfr-alfl)< -1.0E-4_prec*df) RETURN
IF (betr-betl  == zero) THEN
iterm = -6
RETURN
END IF
r = (alfr-alfl)/ (betr-betl)
IF (ABS(t-w) < ABS(t-r)) r = w
q = t
t = r
IF (ABS(t-q) < 1.0E-3_prec) RETURN
tmparray(l-1) = - one
w = t*bet - alf
l = -1

DO j=1,k
alf = tmparray(j-1)
IF (alf < zero) EXIT
bet = tmparray(ma+j-1)
r = t*bet - alf
IF (r > w) THEN
w = r
l = j
END IF
END DO

IF (l <= 0) RETURN
bet = tmparray(ma+l-1)
q = bet - t*df
IF (Q == zero) RETURN

! New interval selection

alf = tmparray(l-1)

IF (q < zero) THEN
alfl = alf
betl = bet
ELSE
alfr = alf
betr = bet
END IF

END DO ns_iteration

END SUBROUTINE nulstep

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE tinit *                                             *
!*                                                                      *
!*     Initial stepsize selection for limited memory bundle method      *
!*                                                                      *
!************************************************************************

SUBROUTINE tinit(n,na,mal,x,af,ag,ay,ibun,d,fu,p,t,tmax,tmin,eta,iters,iterm)

USE param, ONLY : zero,one

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
x, &      ! Vector of variables (n array).
d         ! Direction vector (n array).
REAL(KIND=prec), DIMENSION(0:(n*na)-1), INTENT(IN) :: &
ay, &     ! Matrix whose columns are bundle points (stored in one-dimensional n*ma -array).
ag        ! Matrix whose columns are bundle subgradients (stored in one-dimensional n*ma -array).
REAL(KIND=prec), DIMENSION(0:na-1), INTENT(INOUT) :: &
af        ! Vector of values of bundle functions (stored in one-dimensional ma -array).

! Scalar arguments
REAL(KIND=prec), INTENT(OUT) :: & 
t         ! Initial stepsize
REAL(KIND=prec), INTENT(IN) :: & 
p, &      ! Directional derivative.
eta, &    ! Distance measure parameter.
fu, &     ! Value of the objective function.
tmax, &   ! Upper limit for stepsize parameter.
tmin      ! Lower limit for stepsize parameter.
INTEGER, INTENT(IN) :: &
n, &      ! Number of variables.
na, &     ! Maximum size of the bundle.
mal, &    ! Current size of the bundle.
ibun, &   ! Index for the circular arrays in bundle updating.
iters     ! Null step indicator. 0  - Null step. 1  - Serious step.
INTEGER, INTENT(OUT) :: &
iterm     ! Cause of termination: 0  - Everything is ok. -6  - Error.

! Intrinsic functions
INTRINSIC MAX,MIN

t = MIN(one,tmax)

IF (p == zero) RETURN

IF (iters == 1) THEN
CALL destep(n,na,mal,x,af,ag,ay,ibun,d,fu,p,t,eta,iterm)
ELSE
CALL nulstep(n,na,mal,x,af,ag,ay,ibun,d,fu,p,t,eta,iterm)
END IF

t = MIN(MAX(t,tmin),tmax)

END SUBROUTINE tinit

FUNCTION qint(tu,fl,fuv,xnorm,kappa) RESULT(t)  ! Quadratic interpolation.

USE param, ONLY : half,one

IMPLICIT NONE

! Scalar arguments
REAL(KIND=prec), INTENT(IN) :: & 
fl, &     ! Value of the objective function.
fuv, &    ! Value of the objective function for t=tu.
xnorm, &  ! Directional derivative.
tu, &     ! Upper value of the stepsize parameter.
kappa     ! Interpolation parameter.
REAL(KIND=prec) :: &
t         ! Stepsize.

! Local scalars
REAL(KIND=prec) :: tmp1,tmp2

! Intrinsic functions
INTRINSIC MAX

tmp1 = (fuv-fl)/ (-xnorm*tu)

! Quadratic interpolation with one directional derivative

tmp2 = 2.0_prec * (one - tmp1)

IF (tmp2 > one) THEN

! Interpolation accepted

t = MAX(kappa*tu,tu/tmp2)
RETURN
END IF

! Bisection

t = half*tu

END FUNCTION qint

FUNCTION sclpar(mcc,iscale,sts,stu,utu) RESULT(spar)  ! Calculation of the scaling parameter.

USE param, ONLY : small,one,half

IMPLICIT NONE

! Scalar arguments
REAL(KIND=prec), INTENT(IN) :: &
sts, &    ! sts = trans(s)*s. 
stu, &    ! stu = trans(s)*u. 
utu       ! utu = trans(u)*u. 
REAL(KIND=prec) :: &
spar      ! Scaling parameter.
INTEGER, INTENT(IN) :: & 
mcc, &    ! Current number of stored corrections. SEURAAVALLA RIVILLA OLI ALUNPERIN method, & ! Selection of the method: 0  - Limited memory bundle method. 1  - L-BFGS bundle method.
iscale    ! Selection of the scaling: 0  - Scaling at every iteration with STU/UTU. 1  - Scaling at every iteration with STS/STU. 2  - Interval scaling with STU/UTU. 3  - Interval scaling with STS/STU. 4  - Preliminary scaling with STU/UTU. 5  - Preliminary scaling with STS/STU. 6  - No scaling.      

! Intrinsic functions
INTRINSIC SQRT

! Computation of scaling parameter.

SELECT CASE(iscale)

! Scaling parameter = STU/UTU

CASE(0,2,4)
IF (utu < SQRT(small)) THEN
spar = one
RETURN
ELSE
spar = stu/utu
END IF

! Scaling parameter = STS/UTU

CASE(1,3,5)
IF (stu < SQRT(small)) THEN
spar = one
RETURN
ELSE
spar = sts/stu
END IF

! No scaling

CASE DEFAULT
spar = one
RETURN
END SELECT

! Scaling

IF (mcc == 0) THEN
IF (spar < 0.01_prec) spar=0.01_prec
IF (spar > 100.0_prec) spar=100.0_prec
ELSE
SELECT CASE(iscale)

! Interval scaling

CASE(2)
IF (spar < 0.6_prec .OR. spar > 6.0_prec) spar = one

CASE(3)
IF (spar < half .OR. spar > 5.0_prec) spar = one

! Preliminary scaling

CASE(4,5)
spar = one

! Scaling at every iteration

CASE DEFAULT
CONTINUE
END SELECT
END IF

IF (spar < 1.0E+03_prec*small) spar = 1.0E+03_prec*small

END FUNCTION sclpar

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE dlbfgs *                                            *
!*                                                                      *
!*     Matrix update and computation of the search direction d = -dm*g  *
!*     by the limited memory BFGS update.                               *
!*                                                                      *
!************************************************************************

SUBROUTINE dlbfgs(n,mc,mcc,inew,ibfgs,iflag,d,g,gp,s,u,sm,um,rm,umtum,cm,smtgp,umtgp,gamma,tmpn1,iscale)

USE param, ONLY : zero,small,one
USE lmbm_sub, ONLY: copy2,rwaxv2,vneg,xdiffy,vxdiag,symax,trlieq,scsum,cwmaxv,scdiff,xsumy,vdot

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
g, &      ! Current subgradient of the objective function.
gp        ! Previous subgradient of the objective function.
REAL(KIND=prec), DIMENSION(0:(n*(mcc+1))-1), INTENT(INOUT) :: &
sm, &     ! Matrix whose columns are stored corrections.
um        ! Matrix whose columns are stored subgradient differences.
REAL(KIND=prec), DIMENSION(0:((mcc+2)*(mcc+1)/2)-1), INTENT(INOUT) :: &
rm, &     ! Upper triangular matrix.
umtum     ! Matrix umtum = trans(um) * um.
REAL(KIND=prec), DIMENSION(0:mcc), INTENT(INOUT) :: &
cm, &     ! Diagonal matrix.
smtgp, &  ! Vector smtgp = trans(sm)*gp.
umtgp     ! Vector umtgp = trans(um)*gp.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(INOUT) :: &
s, &      ! Difference of current and previous variables.
u, &      ! Difference of current and previous subgradients.
tmpn1     ! Auxiliary array. On input: previous aggregate subgradient.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(OUT) :: &
d         ! Direction vector.

! Scalar arguments
REAL(KIND=prec), INTENT(INOUT) :: &
gamma     ! Scaling parameter.
INTEGER, INTENT(IN) :: &
n, &      ! Number of variables.
mc, &     ! Declared number of stored corrections. SEURAAVALLA RIVILLA OLI ALUNPERIN method, & ! Selection of the method: 0  - Limited memory bundle method. 1  - L-BFGS bundle method.
iscale    ! Selection of the scaling: 0  - Scaling at every iteration with STU/UTU. 1  - Scaling at every iteration with STS/STU. 2  - Interval scaling with STU/UTU. 3  - Interval scaling with STS/STU. 4  - Preliminary scaling with STU/UTU. 5  - Preliminary scaling with STS/STU. 6  - No scaling.      
INTEGER, INTENT(INOUT) :: &
mcc, &    ! Current number of stored corrections.
inew, &   ! Index for circular arrays.
iflag     ! Index for adaptive version: 0  - Maximum number of stored corrections has not been changed at this iteration. 1  - Maximum number of stored corrections has been changed at this iteration.
INTEGER, INTENT(OUT) :: & 
ibfgs     ! Index of the type of BFGS update: 1  - BFGS update: the corrections are stored. 2  - BFGS update: the corrections are not stored. 3  - BFGS update is skipped.

! Local arrays
REAL(KIND=prec), DIMENSION(0:mcc) :: &
tmpmc1,tmpmc2,tmpmc3,tmpmc4

! Local scalars
REAL(KIND=prec) :: &
stu, &    ! stu = trans(s)*u. 
sts       ! sts = trans(s)*s. 
INTEGER :: i,j,k, &
mcnew, &  ! Current size of vectors.
iold, &   ! Index of the oldest corrections.
iflag2, & ! Index for adaptive version.
ierr      ! Error indicator

! Intrinsic functions
INTRINSIC SQRT,MIN,MAX

ierr = 0
ibfgs = 0
iflag2 = 0
stu = vdot(n,s,u)
sts = vdot(n,s,s)

! Positive definiteness

IF (stu > zero) THEN
IF (-vdot(n,d,u)-vdot(n,tmpn1,s) < -small) THEN

! Update matrices

ibfgs = 1

! Initialization of indices

CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,ibfgs)

! Update sm and um

CALL copy2(n,s,sm((inew-1)*n:),u,um((inew-1)*n:))

! Computation of trans(sm)*g and trans(um)*g

IF (inew >= mcnew) THEN
CALL rwaxv2(n,mcnew,sm((inew-mcnew)*n:),um((inew-mcnew)*n:),g,g,tmpmc1(iold-1:),tmpmc2(iold-1:))
ELSE
CALL rwaxv2(n,inew,sm,um,g,g,tmpmc1,tmpmc2)
CALL rwaxv2(n,mcnew-inew,sm((iold-1)*n:),um((iold-1)*n:),g,g,tmpmc1(iold-1:),tmpmc2(iold-1:))
END IF

! Computation of trans(sm)*u and trans(um)*u

IF (inew >= mcnew) THEN
DO i=iold,inew-1
tmpmc3(i-1) = tmpmc1(i-1) - smtgp(i-1)
smtgp(i-1)  = tmpmc1(i-1)
tmpmc4(i-1) = tmpmc2(i-1) - umtgp(i-1)
umtgp(i-1)  = tmpmc2(i-1)
END DO
ELSE
DO i=1,inew-1
tmpmc3(i-1) = tmpmc1(i-1) - smtgp(i-1)
smtgp(i-1)  = tmpmc1(i-1)
tmpmc4(i-1) = tmpmc2(i-1) - umtgp(i-1)
umtgp(i-1)  = tmpmc2(i-1)
END DO
DO i=iold,mcnew+1
tmpmc3(i-1) = tmpmc1(i-1) - smtgp(i-1)
smtgp(i-1)  = tmpmc1(i-1)
tmpmc4(i-1) = tmpmc2(i-1) - umtgp(i-1)
umtgp(i-1)  = tmpmc2(i-1)
END DO
END IF
tmpmc3(inew-1) = tmpmc1(inew-1) - vdot(n,s,gp)
smtgp(inew-1)  = tmpmc1(inew-1)
tmpmc4(inew-1) = tmpmc2(inew-1) - vdot(n,u,gp)
umtgp(inew-1)  = tmpmc2(inew-1)

! Update rm and umtum

IF (mcc >= mc .AND. iflag2 /= 1) THEN
DO i=1,mcnew-1
j=(i-1)*i/2+1
k=i*(i+1)/2+2
CALL copy2(i,rm(k-1:),rm(j-1:),umtum(k-1:),umtum(j-1:))
END DO
END IF
IF (inew >= mcnew) THEN
CALL copy2(mcnew,tmpmc3(iold-1:),rm((mcnew-1)*mcnew/2:),tmpmc4(iold-1:),umtum((mcnew-1)*mcnew/2:))
ELSE
CALL copy2(mcnew-inew,tmpmc3(iold-1:),rm((mcnew-1)*mcnew/2:),tmpmc4(iold-1:),umtum((mcnew-1)*mcnew/2:))
CALL copy2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew:),tmpmc4,umtum((mcnew-1)*mcnew/2+mcnew-inew:))
END IF

! Update cm

cm(inew-1) = stu

! Computation of gamma

gamma = sclpar(mcc,iscale,sts,stu,tmpmc4(inew-1))
inew = inew + 1
IF (inew > mc + 1) inew = 1
IF (iflag == 0 .AND. mcc < mc + 1) mcc = mcc + 1
ELSE

! BFGS update, corrections are not saved.

ibfgs = 2

! Initialization of indices.

CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,ibfgs)

! Update sm and um

CALL copy2(n,s,sm((inew-1)*n:),u,um((inew-1)*n:))

! Computation of trans(sm)*g and trans(um)*g

CALL rwaxv2(n,mcnew,sm,um,g,g,tmpmc1,tmpmc2)

! Computation of trans(sm)*u and trans(um)*u

IF (iold /= 1) THEN
DO i=1,inew-1
tmpmc3(i-1) = tmpmc1(i-1) - smtgp(i-1)
smtgp(i-1)  = tmpmc1(i-1)
tmpmc4(i-1) = tmpmc2(i-1) - umtgp(i-1)
umtgp(i-1)  = tmpmc2(i-1)
END DO
DO i=iold,mcnew
tmpmc3(i-1) = tmpmc1(i-1) - smtgp(i-1)
smtgp(i-1)  = tmpmc1(i-1)
tmpmc4(i-1) = tmpmc2(i-1) - umtgp(i-1)
umtgp(i-1)  = tmpmc2(i-1)
END DO
ELSE
DO i=1,mcnew-1
tmpmc3(i-1) = tmpmc1(i-1) - smtgp(i-1)
smtgp(i-1)  = tmpmc1(i-1)
tmpmc4(i-1) = tmpmc2(i-1) - umtgp(i-1)
umtgp(i-1)  = tmpmc2(i-1)
END DO
END IF
tmpmc3(inew-1) = tmpmc1(inew-1) - vdot(n,s,gp)
smtgp(inew-1)  = tmpmc1(inew-1)
tmpmc4(inew-1) = tmpmc2(inew-1) - vdot(n,u,gp)
umtgp(inew-1)  = tmpmc2(inew-1)

! Update rm and umtum

IF (iold /= 1) THEN
CALL copy2(mcnew-inew,tmpmc3(iold-1:),rm((mcnew-1)*mcnew/2:),tmpmc4(iold-1:),umtum((mcnew-1)*mcnew/2:))
CALL copy2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew:),tmpmc4,umtum((mcnew-1)*mcnew/2+mcnew-inew:))
ELSE
CALL copy2(mcnew,tmpmc3,rm((mcnew-1)*mcnew/2:),tmpmc4,umtum((mcnew-1)*mcnew/2:))
END IF

! Update cm

cm(inew-1) = stu

! Computation of gamma

gamma = sclpar(mcc,iscale,sts,stu,tmpmc4(inew-1))
END IF
ELSE

! BFGS update is skipped

ibfgs = 3
IF (mcc == 0) THEN
iflag = 0
CALL vneg(n,g,d)
RETURN
END IF

! Initialization of indices.

CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,ibfgs)

! Computation of gamma

IF (iscale >= 4) gamma = one

! Computation of trans(sm)*g and trans(um)*g and the two intermediate values

IF (iold <= 2) THEN
CALL rwaxv2(n,mcnew,sm((iold-1)*n:),um((iold-1)*n:),g,g,smtgp(iold-1:),umtgp(iold-1:))
ELSE
CALL rwaxv2(n,inew-1,sm,um,g,g,smtgp,umtgp)
CALL rwaxv2(n,mcnew-inew+1,sm((iold-1)*n:),um((iold-1)*n:),g,g,smtgp(iold-1:),umtgp(iold-1:))
END IF
END IF

! Computation of two intermediate values tmpmc1 and tmpmc2

IF (iold == 1 .OR. ibfgs == 2) THEN
CALL trlieq(mcnew,mcnew,iold,rm,tmpmc1,smtgp,1,ierr)
CALL symax(mcnew,mcnew,iold,umtum,tmpmc1,tmpmc3)
CALL vxdiag(mcnew,cm,tmpmc1,tmpmc2)
CALL scsum(mcnew,gamma,tmpmc3,tmpmc2,tmpmc2)
CALL scsum(mcnew,-gamma,umtgp,tmpmc2,tmpmc3)
CALL trlieq(mcnew,mcnew,iold,rm,tmpmc2,tmpmc3,0,ierr)
ELSE IF (iflag == 0) THEN
CALL trlieq(mcnew,mc+1,iold,rm,tmpmc1,smtgp,1,ierr)
CALL symax(mcnew,mc+1,iold,umtum,tmpmc1,tmpmc3)
CALL vxdiag(mc+1,cm,tmpmc1,tmpmc2)
CALL scsum(mc+1,gamma,tmpmc3,tmpmc2,tmpmc2)
CALL scsum(mc+1,-gamma,umtgp,tmpmc2,tmpmc3)
CALL trlieq(mcnew,mc+1,iold,rm,tmpmc2,tmpmc3,0,ierr)
ELSE
CALL trlieq(mcnew,mc,iold,rm,tmpmc1,smtgp,1,ierr)
CALL symax(mcnew,mc,iold,umtum,tmpmc1,tmpmc3)
CALL vxdiag(mc,cm,tmpmc1,tmpmc2)
CALL scsum(mc,gamma,tmpmc3,tmpmc2,tmpmc2)
CALL scsum(mc,-gamma,umtgp,tmpmc2,tmpmc3)
CALL trlieq(mcnew,mc,iold,rm,tmpmc2,tmpmc3,0,ierr)
END IF

! Computation of the search direction d

IF (iold == 1 .OR. ibfgs == 2) THEN
CALL cwmaxv(n,mcnew,um,tmpmc1,d)
CALL xdiffy(n,d,g,d)
CALL cwmaxv(n,mcnew,sm,tmpmc2,tmpn1)
CALL scdiff(n,gamma,d,tmpn1,d)
ELSE
CALL cwmaxv(n,inew-1,um,tmpmc1,d)
CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n:),tmpmc1(iold-1:),tmpn1)
CALL xsumy(n,d,tmpn1,d)
CALL xdiffy(n,d,g,d)
CALL cwmaxv(n,inew-1,sm,tmpmc2,tmpn1)
CALL scdiff(n,gamma,d,tmpn1,d)
CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n:),tmpmc2(iold-1:),tmpn1)
CALL xdiffy(n,d,tmpn1,d)
END IF

END SUBROUTINE dlbfgs

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE dlsr1 *                                             *
!*                                                                      *
!*     Matrix update and computation of the search direction d = -dm*ga *
!*     by the limited memory SR1 update.                                *
!*                                                                      *
!************************************************************************

SUBROUTINE dlsr1(n,mc,mcc,inew,isr1,iflag,d,gp,ga,s,u,sm,um,rm,umtum,cm,smtgp,umtgp,gamma,tmpmc1,tmpmc2,tmpn1,nnk)

USE param, ONLY : zero,small,one
USE lmbm_sub, ONLY : vdot,vneg,scalex,xdiffy,scdiff,xsumy,cwmaxv,rwaxv2,calq,copy,copy2

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
ga, &     ! Current aggregate subgradient of the objective function.
gp        ! Basic subgradient of the objective function.
REAL(KIND=prec), DIMENSION(0:(n*(mcc+1))-1), INTENT(INOUT) :: &
sm, &     ! Matrix whose columns are stored corrections.
um        ! Matrix whose columns are stored subgradient differences.
REAL(KIND=prec), DIMENSION(0:((mcc+2)*(mcc+1)/2)-1), INTENT(INOUT) :: &
rm, &     ! Upper triangular matrix.
umtum     ! Matrix umtum = trans(um) * um.
REAL(KIND=prec), DIMENSION(0:mcc), INTENT(INOUT) :: &
cm, &     ! Diagonal matrix.
smtgp, &  ! Vector smtgp = trans(sm)*gp.
umtgp     ! Vector umtgp = trans(um)*gp.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(INOUT) :: &
s, &      ! Difference of current and previous variables.
u, &      ! Difference of current and previous subgradients.
tmpn1     ! Auxiliary array. On input: previous aggregate subgradient.
REAL(KIND=prec), DIMENSION(0:mcc), INTENT(OUT) :: &
tmpmc1, & ! Auxiliary array. On output: trans(sm)*ga.
tmpmc2    ! Auxiliary array. On output: trans(um)*ga.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(OUT) :: &
d         ! Direction vector.

! Scalar arguments
REAL(KIND=prec), INTENT(INOUT) :: & 
gamma     ! Scaling parameter.
INTEGER, INTENT(IN) :: &
mc, &     ! Declared number of stored corrections.
nnk       ! Consecutive null steps counter.
INTEGER, INTENT(INOUT) :: &
n, &      ! Number of variables.
mcc, &    ! Current number of stored corrections.
inew, &   ! Index for circular arrays.
iflag     ! Index for adaptive version: 0  - Maximum number of stored corrections has not been changed at this iteration. 1  - Maximum number of stored corrections has been changed at this iteration.
INTEGER, INTENT(OUT) :: & 
isr1      ! Index of the type of L-SR1 update: 1  - SR1 update: the corrections are stored. 3  - SR1 update is skipped.

! Local arrays
REAL(KIND=prec), DIMENSION(0:n-1) :: tmpn2
REAL(KIND=prec), DIMENSION(0:((mcc+1)*(mcc+2)/2)-1) :: tmpmat
REAL(KIND=prec), DIMENSION(0:mcc) :: tmpmc3,tmpmc4,tmpmc5,tmpmc6

! Local scalars
REAL(KIND=prec) :: &
stu, &    ! stu = trans(s)*u. 
a, &      ! a = trans(ga) dm_(k-1) ga.
b         ! b = trans(ga) dm_k ga.
INTEGER :: i,j,k, &
mcnew, &  ! Current size of vectors.
iold, &   ! Index of the oldest corrections.
iflag2    ! Index for adaptive version.

iflag2 = 0
isr1 = 0 

! Initialization of indices

CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,3)

! Computation of gamma

gamma = one

! Computation of trans(sm)*ga and trans(um)*ga

IF (iold <= 2) THEN
CALL rwaxv2(n,mcnew,sm((iold-1)*n:),um((iold-1)*n:),ga,ga,tmpmc1(iold-1:),tmpmc2(iold-1:))
ELSE
CALL rwaxv2(n,inew-1,sm,um,ga,ga,tmpmc1,tmpmc2)
CALL rwaxv2(n,mcnew-inew+1,sm((iold-1)*n:),um((iold-1)*n:),ga,ga,tmpmc1(iold-1:),tmpmc2(iold-1:))
END IF

! Positive definiteness

IF (-vdot(n,d,u) - vdot(n,tmpn1,s) >= -small) THEN

! SR1 update is skipped

isr1 = 3
IF (mcc == 0) THEN
iflag = 0
CALL vneg(n,ga,d)
RETURN
END IF
ELSE
stu = vdot(n,s,u)
tmpmc1(inew-1) = vdot(n,s,ga)
tmpmc2(inew-1) = vdot(n,u,ga)

! Convergence conditions

IF ((nnk == 1 .OR. mcc < mc) .OR. (iflag == 1 .AND. (inew == 1 .OR. inew == mc))) THEN

! SR1 update

isr1 = 1

! Initialization of indices

CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,1)        
IF (iflag2 == 1 .AND. iold == 2) THEN
tmpmc1(inew-1) = tmpmc1(0)
tmpmc2(inew-1) = tmpmc2(0)
END IF

! Update sm and um

CALL copy2(n,s,sm((inew-1)*n:),u,um((inew-1)*n:))

! Update trans(sm)*gp and trans(um)*gp

smtgp(inew-1) = vdot(n,s,gp)
umtgp(inew-1) = vdot(n,u,gp)

! Computation of trans(sm)*u and trans(um)*u

IF (iold <= 2) THEN
CALL rwaxv2(n,mcnew-1,sm((iold-1)*n:),um((iold-1)*n:),u,u,tmpmc3(iold-1:),tmpmc4(iold-1:))
ELSE
CALL rwaxv2(n,inew-1,sm,um,u,u,tmpmc3,tmpmc4)
CALL rwaxv2(n,mcnew-inew,sm((iold-1)*n:),um((iold-1)*n:),u,u,tmpmc3(iold-1:),tmpmc4(iold-1:))
END IF

tmpmc3(inew-1) = stu
tmpmc4(inew-1) = vdot(n,u,u)

! Update rm and umtum

IF (mcc >= mc .AND. iflag2 /= 1) THEN
DO i=1,mcnew-1
j=(i-1)*i/2+1
k=i*(i+1)/2+2
CALL copy2(i,rm(k-1:),rm(j-1:),umtum(k-1:),umtum(j-1:))
END DO
END IF

IF (inew >= mcnew) THEN
CALL copy2(mcnew,tmpmc3(iold-1:),rm((mcnew-1)*mcnew/2:),tmpmc4(iold-1:),umtum((mcnew-1)*mcnew/2:))
ELSE
CALL copy2(mcnew-inew,tmpmc3(iold-1:),rm((mcnew-1)*mcnew/2:),tmpmc4(iold-1:),umtum((mcnew-1)*mcnew/2:))
CALL copy2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew:),tmpmc4,umtum((mcnew-1)*mcnew/2+mcnew-inew:))
END IF

! Update cm

cm(inew-1) = stu
inew = inew + 1
IF (inew > mc + 1) inew = 1
IF (iflag == 0 .AND. mcc < mc + 1) mcc = mcc + 1
ELSE

! Calculation of matrix (umtum-rm-trans(rm)+cm) from previous iteration

DO  i=1,mcnew*(mcnew+1)/2
tmpmat(i-1)= gamma * umtum(i-1) - rm(i-1)
END DO

! Computation of tmpmat*tmpmc4 = gamma*trans(um)*ga-trans(sm)*ga

IF (iold == 1) THEN
CALL scdiff(mcnew,gamma,tmpmc2,tmpmc1,tmpmc5)
CALL calq(mcnew,mcnew,iold,tmpmat,tmpmc4,tmpmc5)
ELSE IF (iflag == 0) THEN
CALL scdiff(mc+1,gamma,tmpmc2,tmpmc1,tmpmc5)
CALL calq(mcnew,mc+1,iold,tmpmat,tmpmc4,tmpmc5)
ELSE
CALL scdiff(mc,gamma,tmpmc2,tmpmc1,tmpmc5)
CALL calq(mcnew,mc,iold,tmpmat,tmpmc4,tmpmc5)
END IF

! Computation of a = -trans(ga)*dm_(k-1)*ga

IF (iold <= 2) THEN
CALL scalex(mcnew,gamma,tmpmc4(iold-1:),tmpmc3(iold-1:))
CALL cwmaxv(n,mcnew,sm((iold-1)*n:),tmpmc4(iold-1:),tmpn1)
CALL scdiff(n,-gamma,ga,tmpn1,tmpn2)
CALL cwmaxv(n,mcnew,um((iold-1)*n:),tmpmc3(iold-1:),tmpn1)
CALL xsumy(n,tmpn2,tmpn1,tmpn2)
ELSE
CALL scalex(mcc,gamma,tmpmc4,tmpmc3)
CALL cwmaxv(n,inew-1,sm,tmpmc4,tmpn1)
CALL scdiff(n,-gamma,ga,tmpn1,tmpn2)
CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n:),tmpmc4(iold-1:),tmpn1)
CALL xdiffy(n,tmpn2,tmpn1,tmpn2)
CALL cwmaxv(n,inew-1,um,tmpmc3,tmpn1)
CALL xsumy(n,tmpn2,tmpn1,tmpn2)
CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n:),tmpmc3(iold-1:),tmpn1)
CALL xsumy(n,tmpn2,tmpn1,tmpn2)
END IF

a = vdot(n,ga,tmpn2)

IF (iflag == 0) THEN
mcnew = mc
iold = inew + 2
IF (iold > mc+1) iold = iold - mc - 1
ELSE
mcnew = mc - 1
iold = inew + 2
IF (iold > mc) iold = iold - mc
END IF

! Calculation of the new canditate for search direction
! Updates are not necessarily saved

! Update sm and um

CALL copy2(n,s,sm((inew-1)*n:),u,um((inew-1)*n:))

! Computation of trans(sm)*u and trans(um)*u

IF (iold == 1 .OR. iold == 2) THEN
CALL rwaxv2(n,mcnew-1,sm((iold-1)*n:),um((iold-1)*n:),u,u,tmpmc3(iold-1:),tmpmc4(iold-1:))
ELSE
CALL rwaxv2(n,inew-1,sm,um,u,u,tmpmc3,tmpmc4)
CALL rwaxv2(n,mcnew-inew,sm((iold-1)*n:),um((iold-1)*n:),u,u,tmpmc3(iold-1:),tmpmc4(iold-1:))
END IF

tmpmc3(inew-1) = stu
tmpmc4(inew-1) = vdot(n,u,u)

! Calculation of matrix (umtum-rm-trans(rm)+cm) without updating matrices rm, umtum and cm

DO i=1,mcnew*(mcnew+1)/2
tmpmat(i-1)= gamma * umtum(i-1) - rm(i-1)
END DO

DO i=1,mcnew-1
j=(i-1)*i/2+1
k=i*(i+1)/2+2
CALL copy(i,tmpmat(k-1:),tmpmat(j-1:))
END DO

CALL scdiff(mcnew+1,gamma,tmpmc4,tmpmc3,tmpmc5)

IF (inew >= mcnew) THEN
CALL copy(mcnew,tmpmc5(iold-1:),tmpmat((mcnew-1)*mcnew/2:))
ELSE
CALL copy(mcnew-inew,tmpmc5(iold-1:),tmpmat((mcnew-1)*mcnew/2:))
CALL copy(inew,tmpmc5,tmpmat((mcnew-1)*mcnew/2+mcnew-inew:))
END IF

IF (iflag == 0) THEN
CALL scdiff(mc+1,gamma,tmpmc2,tmpmc1,tmpmc5)
CALL calq(mcnew,mc+1,iold,tmpmat,tmpmc5,tmpmc5)
ELSE
CALL scdiff(mc,gamma,tmpmc2,tmpmc1,tmpmc5)
CALL calq(mcnew,mc,iold,tmpmat,tmpmc5,tmpmc5)
END IF

! Calculation of the new canditate for search direction d = -dm_k*ga and computation of b = -trans(ga)*dm_k*ga

IF (iold <= 2) THEN
CALL scalex(mcnew,gamma,tmpmc5(iold-1:),tmpmc6(iold-1:))
CALL cwmaxv(n,mcnew,sm((iold-1)*n:),tmpmc5(iold-1:),tmpn1)
CALL scdiff(n,-gamma,ga,tmpn1,d)
CALL cwmaxv(n,mcnew,um((iold-1)*n:),tmpmc6(iold-1:),tmpn1)
CALL xsumy(n,d,tmpn1,d)
ELSE
CALL scalex(mcnew+1,gamma,tmpmc5,tmpmc6)
CALL cwmaxv(n,inew,sm,tmpmc5,tmpn1)
CALL scdiff(n,-gamma,ga,tmpn1,d)
CALL cwmaxv(n,mcnew-inew,sm((iold-1)*n:),tmpmc5(iold-1:),tmpn1)
CALL xdiffy(n,d,tmpn1,d)
CALL cwmaxv(n,inew,um,tmpmc6,tmpn1)
CALL xsumy(n,d,tmpn1,d)
CALL cwmaxv(n,mcnew-inew,um((iold-1)*n:),tmpmc6(iold-1:),tmpn1)
CALL xsumy(n,d,tmpn1,d)
END IF

b = vdot(n,ga,d)

! Checking the convergence conditions

IF (b - a < zero) THEN
isr1 = 3
CALL copy(n,tmpn2,d)
ELSE
isr1 = 1

! Update trans(sm)*gp and trans(um)*gp

smtgp(inew-1) = vdot(n,s,gp)
umtgp(inew-1) = vdot(n,u,gp)

! Update rm and umtum

DO i=1,mcnew-1
j=(i-1)*i/2+1
k=i*(i+1)/2+2
CALL copy2(i,rm(k-1:),rm(j-1:),umtum(k-1:),umtum(j-1:))
END DO

IF (inew >= mcnew) THEN
CALL copy2(mcnew,tmpmc3(iold-1:),rm((mcnew-1)*mcnew/2:),tmpmc4(iold-1:),umtum((mcnew-1)*mcnew/2:))
ELSE
CALL copy2(mcnew-inew,tmpmc3(iold-1:),rm((mcnew-1)*mcnew/2:),tmpmc4(iold-1:),umtum((mcnew-1)*mcnew/2:))
CALL copy2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew:),tmpmc4,umtum((mcnew-1)*mcnew/2+mcnew-inew:))
END IF

! Update cm

cm(inew-1) = stu
inew = inew + 1
IF (inew > mc + 1) inew = 1
IF (iflag == 0 .AND. mcc < mc + 1) mcc = mcc + 1
END IF
RETURN

END IF

END IF

DO i=1,mcnew*(mcnew+1)/2
tmpmat(i-1)= gamma * umtum(i-1) - rm(i-1)
END DO

! Computation of tmpmat*tmpmc4 = gamma*trans(um)*ga-trans(sm)*ga

IF (iold == 1) THEN
CALL scdiff(mcnew,gamma,tmpmc2,tmpmc1,tmpmc4)
CALL calq(mcnew,mcnew,iold,tmpmat,tmpmc4,tmpmc4)
ELSE IF (iflag == 0) THEN
CALL scdiff(mc+1,gamma,tmpmc2,tmpmc1,tmpmc4)
CALL calq(mcnew,mc+1,iold,tmpmat,tmpmc4,tmpmc4)
ELSE
CALL scdiff(mc,gamma,tmpmc2,tmpmc1,tmpmc4)
CALL calq(mcnew,mc,iold,tmpmat,tmpmc4,tmpmc4)
END IF

! Computation of the search direction d

IF (iold <= 2) THEN
CALL scalex(mcnew,gamma,tmpmc4(iold-1:),tmpmc3(iold-1:))
CALL cwmaxv(n,mcnew,sm((iold-1)*n:),tmpmc4(iold-1:),tmpn1)
CALL scdiff(n,-gamma,ga,tmpn1,d)
CALL cwmaxv(n,mcnew,um((iold-1)*n:),tmpmc3(iold-1:),tmpn1)
CALL xsumy(n,d,tmpn1,d)
ELSE
CALL scalex(mcc,gamma,tmpmc4,tmpmc3)
CALL cwmaxv(n,inew-1,sm,tmpmc4,tmpn1)
CALL scdiff(n,-gamma,ga,tmpn1,d)
CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n:),tmpmc4(iold-1:),tmpn1)
CALL xdiffy(n,d,tmpn1,d)
CALL cwmaxv(n,inew-1,um,tmpmc3,tmpn1)
CALL xsumy(n,d,tmpn1,d)
CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n:),tmpmc3(iold-1:),tmpn1)
CALL xsumy(n,d,tmpn1,d)
END IF

END SUBROUTINE dlsr1

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE indic1 *                                            *
!*                                                                      *
!*     Initialization of indices.                                       *
!*                                                                      *
!************************************************************************

SUBROUTINE indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,itype)

IMPLICIT NONE

! Scalar arguments
INTEGER, INTENT(IN) :: & 
mc, &     ! Declared number of stored corrections.
mcc, &    ! Current number of stored corrections.
itype     ! Type of Initialization: 1  - corrections are stored, 2  - corrections are not stored, 3  - update is skipped.
INTEGER, INTENT(INOUT) :: & 
inew, &   ! Index for circular arrays.
iflag, &  ! Index for adaptive version: 0  - Maximum number of stored corrections has not been changed at this iteration. 1  - Maximum number of stored corrections has been changed at this iteration.
iflag2    ! Index for adaptive version. 0  - iflag has not been changed. 1  - iflag has been changed.
INTEGER, INTENT(OUT) :: & 
mcnew, &  ! Current size of vectors.
iold      ! Index of the oldest corrections.

IF (itype == 1) THEN
IF (mcc < mc) THEN
mcnew = mcc + 1
iold = 1
iflag = 0
ELSE
IF (iflag == 0) THEN
mcnew = mc
iold = inew + 2
IF (iold > mc+1) THEN
iold = iold - mc - 1
END IF
ELSE
IF (inew == 1) THEN
inew = mc + 1
mcnew = mc
iold = 2
iflag = 0
iflag2 = 1
ELSE IF (inew == mc) THEN
mcnew = mc
iold = 1
iflag = 0
iflag2 = 1
ELSE
mcnew = mc - 1
iold = inew + 2
IF (iold > mc) THEN
iold = iold - mc
END IF
END IF
END IF
END IF
ELSE IF (itype == 2) THEN
IF (mcc < mc) THEN
mcnew = mcc + 1
iold = 1
iflag = 0
ELSE
IF (iflag == 0) THEN
mcnew = mc + 1
iold = inew + 1
IF (iold > mc + 1) THEN
iold = 1
END IF
ELSE
mcnew = mc
iold = inew + 1
IF (iold > mc) THEN
iold = 1
END IF
END IF
END IF
ELSE
IF (mcc < mc) THEN
mcnew = mcc
iold = 1
iflag = 0
ELSE
IF (iflag == 0) THEN
mcnew = mc
iold = inew + 1
IF (iold > mc + 1) THEN
iold = 1
END IF
ELSE
mcnew = mc - 1
iold = inew + 1
IF (iold > mc) THEN
iold = 1
END IF
END IF
END IF
END IF

END SUBROUTINE indic1

!************************************************************************
!*
!*     * SUBROUTINE agbfgs *
!*
!*     Computation of aggregate values by the limited memory BFGS update.
!*
!************************************************************************

SUBROUTINE agbfgs(n,mc,mcc,inew,ibfgs,iflag,g,gp,ga,u,d,sm,um,rm,cm,umtum,alfn,alfv,gamma,ic,rho)

USE param, ONLY : zero,half,one
USE lmbm_sub, ONLY : symax,rwaxv2,trlieq,vdot

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
d, &      ! Direction vector.
g, &      ! Current (auxiliary) subgradient of the objective function.
gp, &     ! Previous subgradient of the objective function.
u         ! Difference of trial and aggregate gradients.
REAL(KIND=prec), DIMENSION(0:(n*(mcc+1))-1), INTENT(IN) :: &
sm, &     ! Matrix whose columns are stored corrections.
um        ! Matrix whose columns are stored subgradient differences.
REAL(KIND=prec), DIMENSION(0:((mcc+2)*(mcc+1)/2)-1), INTENT(IN) :: &
rm, &     ! Upper triangular matrix.
umtum     ! Matrix umtum = trans(um) * um.
REAL(KIND=prec), DIMENSION(0:mcc), INTENT(IN) :: &
cm        ! Diagonal matrix.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(OUT) :: &
ga        ! Next aggregate subgradient of the objective function.

! Scalar arguments
REAL(KIND=prec), INTENT(OUT) :: & 
alfv      ! Aggregate locality measure.
REAL(KIND=prec), INTENT(IN) :: & 
gamma, &  ! Scaling parameter.
alfn, &   ! Locality measure.
rho       ! Correction parameter.
INTEGER, INTENT(IN) :: &
n, &      ! Number of variables.
mc, &     ! Declared number of stored corrections.
mcc, &    ! Current number of stored corrections.
inew, &   ! Index for circular arrays.
ibfgs, &  ! Index of the type of BFGS update.
ic        ! Correction indicator.
INTEGER, INTENT(INOUT) :: &
iflag     ! Index for adaptive version: 0  - Maximum number of stored corrections has not been changed at this iteration. 1  - Maximum number of stored corrections has been changed at this iteration.

! Local arrays
REAL(KIND=prec), DIMENSION(0:mcc) :: tmpmc1, tmpmc2

! Local scalars
REAL(KIND=prec) :: &
p, &      ! p = trans(d)*u - alfn.
q, &      ! q = trans(u)*dm*u, where dm is the inverse approximation of the Hessian calculated by using the L-BFGS formula.
lam, &    ! Multiplier used to calculate aggregate values.
w         ! Correction.
INTEGER :: i, &
mcnew, &  ! Current size of vectors.
iold, &   ! Index of the oldest corrections.
ierr      ! Error indicador.

! Intrinsic functions
INTRINSIC MAX,MIN,SIGN

ierr = 0

IF (mcc < mc) THEN
IF (ibfgs == 2) THEN
mcnew = mcc + 1
ELSE
mcnew = mcc
END IF
iold = 1
ELSE
IF (iflag == 0) THEN
IF (ibfgs == 2) THEN
mcnew = mc + 1
ELSE
mcnew = mc
END IF
iold = inew + 1
IF (iold > mc+1) iold = 1
ELSE
IF (ibfgs == 2) THEN
mcnew = mc
ELSE
mcnew = mc - 1
END IF
iold = inew + 1
IF (iold > mc) iold = 1
END IF
END IF

! Computation of trans(d)*u-alfn

p = vdot(n,d,u) - alfn
q = vdot(n,u,u)

IF (ic == 1) THEN
w = rho * q
ELSE
w = zero
END IF

! Computation of the product trans(u)*dm*u

IF (mcc > 0 .OR. ibfgs == 2) THEN
IF (iold == 1 .OR. ibfgs == 2) THEN
CALL rwaxv2(n,mcnew,sm,um,u,u,tmpmc1,tmpmc2)
CALL trlieq(mcnew,mcnew,iold,rm,tmpmc1,tmpmc1,1,ierr)
q = q - 2.0_prec*vdot(mcnew,tmpmc2,tmpmc1)
q = gamma*q
DO i=1,mcnew
tmpmc2(i-1) = cm(i-1)*tmpmc1(i-1)
END DO
q = q + vdot(mcnew,tmpmc1,tmpmc2)
CALL symax(mcnew,mcnew,iold,umtum,tmpmc1,tmpmc2)
q = q + gamma*vdot(mcnew,tmpmc1,tmpmc2)
ELSE
CALL rwaxv2(n,inew-1,sm,um,u,u,tmpmc1,tmpmc2)
CALL rwaxv2(n,mcc-inew,sm((iold-1)*n:),um((iold-1)*n:),u,u,tmpmc1(iold-1:),tmpmc2(iold-1:))
CALL trlieq(mcnew,mcc,iold,rm,tmpmc1,tmpmc1,1,ierr)
q = q - 2.0_prec*(vdot(mcc-inew,tmpmc2(iold-1:),tmpmc1(iold-1:)) + vdot(inew-1,tmpmc2,tmpmc1))
q = gamma*q
DO i=1,mcc
tmpmc2(i-1) = cm(i-1)*tmpmc1(i-1)
END DO
q = q + vdot(mcc-inew,tmpmc1(iold-1:),tmpmc2(iold-1:)) + vdot(inew-1,tmpmc1,tmpmc2)
CALL symax(mcnew,mcc,iold,umtum,tmpmc1,tmpmc2)
q = q + gamma*(vdot(mcc-inew,tmpmc1(iold-1:),tmpmc2(iold-1:)) + vdot(inew-1,tmpmc1,tmpmc2))
END IF
END IF

q = q + w
lam = half + SIGN(half,p)
IF (q > zero) lam = MIN(one,MAX(zero,p/q))

! Computation of the aggregate values

p = one - lam

DO i=1,n
ga(i-1)=lam*g(i-1) + p*gp(i-1)
END DO

alfv = lam*alfn

END SUBROUTINE agbfgs

!************************************************************************
!*
!*     * SUBROUTINE aggsr1 *
!*
!*     Computation of aggregate values by the limited memory SR1 update.
!*
!************************************************************************

SUBROUTINE aggsr1(n,mc,mcc,inew,iflag,g,gp,ga,d,alfn,alfv,umtum,rm,gamma,smtgp,umtgp,smtga,umtga,sm,um,icn,rho)

USE param, ONLY : zero,one,small
USE lmbm_sub, ONLY : vdot,scalex,xsumy,xdiffy,scsum,scdiff,rwaxv2,cwmaxv,lineq,calq

IMPLICIT NONE

! Array arguments
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(IN) :: &
d, &      ! Direction vector.
g, &      ! Current (auxiliary) subgradient of the objective function.
gp        ! Previous subgradient of the objective function.
REAL(KIND=prec), DIMENSION(0:(n*(mcc+1))-1), INTENT(IN) :: &
sm, &     ! Matrix whose columns are stored corrections.
um        ! Matrix whose columns are stored subgradient differences.
REAL(KIND=prec), DIMENSION(0:((mcc+2)*(mcc+1)/2)-1), INTENT(IN) :: &
rm, &     ! Upper triangular matrix.
umtum     ! Matrix umtum = trans(um) * um.
REAL(KIND=prec), DIMENSION(0:mcc), INTENT(IN) :: &
smtgp, &  ! Vector smtgp = trans(sm)*gp.
umtgp, &  ! vector umtgp = trans(um)*gp.
smtga, &  ! vector smtga = trans(sm)*ga.
umtga     ! vector umtga = trans(um)*ga.
REAL(KIND=prec), DIMENSION(0:n-1), INTENT(INOUT) :: &
ga        ! Aggregate subgradient of the objective function.

! Scalar arguments
REAL(KIND=prec), INTENT(INOUT) :: & 
alfv      ! Aggregate locality measure.
REAL(KIND=prec), INTENT(IN) :: & 
gamma, &  ! Scaling parameter.
alfn, &   ! Locality measure.
rho       ! Correction parameter.
INTEGER, INTENT(IN) :: &
n, &      ! Number of variables
mc, &     ! Declared number of stored corrections.
mcc, &    ! Current number of stored corrections.
inew, &   ! Index for circular arrays.
icn       ! Correction indicator.
INTEGER, INTENT(INOUT) :: &
iflag     ! Index for adaptive version: 0  - Maximum number of stored corrections has not been changed at this iteration. 1  - Maximum number of stored corrections has been changed at this iteration.

! Local arrays
REAL(KIND=prec), DIMENSION(0:n-1) :: tmpn2,tmpn3,tmpn4
REAL(KIND=prec), DIMENSION(0:((mcc+1)*(mcc)/2)-1) :: tmpmat
REAL(KIND=prec), DIMENSION(0:mcc) :: tmpmc3, tmpmc4

! Local scalars
REAL(KIND=prec) :: &
pr, &     ! pr = trans(gp-ga) dm (gp-ga), where dm presents the L-SR1- approximation of Hessian.
rrp, &    ! rrp = trans(gp-ga) dm ga - alfv.
prqr, &   ! prqr = trans(gp-ga) dm (g-ga).
rrq, &    ! rrq = trans(g-ga) dm ga - alfv + alfn.
qr, &     ! qr = trans(g-ga) dm (g-ga).
pq, &     ! pq = trans(g-gp) dm (g-gp).
qqp, &    ! qqp = trans(g-gp) dm g + alfn.
lam1, &   ! Multiplier used to calculate aggregate values.
lam2, &   ! Multiplier used to calculate aggregate values.
w, &      ! Correction.
tmp1, &   ! Auxiliary scalar.
tmp2      ! Auxiliary scalar.
INTEGER :: i, &
mcnew, &  ! Current size of vectors.
iold, &   ! Index of the oldest corrections.
ierr      ! Error indicador.

! Intrinsic functions
INTRINSIC MIN,MAX

ierr = 0

IF (mcc < mc) THEN
iold = 1
mcnew = mcc
ELSE IF (iflag == 0) THEN
mcnew = mc
iold = inew + 1
IF (iold > mc+1) iold = 1
ELSE
mcnew = mc - 1
iold = inew + 1
IF (iold > mc) iold = 1
END IF

CALL xdiffy(n,gp,ga,tmpn2)

! Calculation of tmpn3 = trans(gp - ga)dm

IF (mcc > 0) THEN
DO i=1,mcnew*(mcnew+1)/2
tmpmat(i-1)= gamma * umtum(i-1) - rm(i-1)
END DO
IF (iold == 1) THEN
CALL xdiffy(mcnew,umtgp,umtga,tmpmc4)
CALL scdiff(mcnew,gamma,tmpmc4,smtgp,tmpmc4)
CALL xsumy(mcnew,tmpmc4,smtga,tmpmc4)
CALL calq(mcnew,mcnew,iold,tmpmat,tmpmc3,tmpmc4)
CALL scalex(mcnew,gamma,tmpmc3,tmpmc4)
CALL cwmaxv(n,mcnew,sm,tmpmc3,tmpn4)
CALL scsum(n,gamma,tmpn2,tmpn4,tmpn3)
CALL cwmaxv(n,mcnew,um,tmpmc4,tmpn4)
CALL xdiffy(n,tmpn3,tmpn4,tmpn3)
ELSE
CALL xdiffy(mcc,umtgp,umtga,tmpmc4)
CALL scdiff(mcc,gamma,tmpmc4,smtgp,tmpmc4)
CALL xsumy(mcc,tmpmc4,smtga,tmpmc4)
CALL calq(mcnew,mcc,iold,tmpmat,tmpmc3,tmpmc4)
CALL scalex(mcc,gamma,tmpmc3,tmpmc4)
CALL cwmaxv(n,inew-1,sm,tmpmc3,tmpn4)
CALL scsum(n,gamma,tmpn2,tmpn4,tmpn3)
CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n:),tmpmc3(iold-1:),tmpn4)
CALL xsumy(n,tmpn3,tmpn4,tmpn3)
CALL cwmaxv(n,inew-1,um,tmpmc4,tmpn4)
CALL xdiffy(n,tmpn3,tmpn4,tmpn3)
CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n:),tmpmc4(iold-1:),tmpn4)
CALL xdiffy(n,tmpn3,tmpn4,tmpn3)
END IF
IF (icn == 1) THEN
CALL scsum(n,rho,tmpn2,tmpn3,tmpn3)
END IF
pr = vdot(n,tmpn3,tmpn2)
rrp = vdot(n,tmpn3,ga) 
CALL xdiffy(n,g,ga,tmpn4)
prqr = vdot(n,tmpn3,tmpn4)
rrq = -vdot(n,tmpn4,d)
ELSE
pr = vdot(n,tmpn2,tmpn2)
rrp = vdot(n,tmpn2,ga) 
CALL xdiffy(n,g,ga,tmpn4)
prqr = vdot(n,tmpn2,tmpn4)
rrq = -vdot(n,tmpn4,d)
END IF

! Calculation of qr = trans(g - ga) dm (g - ga)

qr = vdot(n,tmpn4,tmpn4)

IF (icn == 1) THEN
w = rho*qr
ELSE
w = zero
END IF

IF (mcc > 0) THEN
qr = gamma*qr
IF (iold == 1) THEN
CALL rwaxv2(n,mcnew,sm,um,tmpn4,tmpn4,tmpmc4,tmpmc3)
CALL scsum(mcnew,-gamma,tmpmc3,tmpmc4,tmpmc4)
CALL lineq(mcnew,mcnew,iold,tmpmat,tmpmc3,tmpmc4,ierr)  
qr = qr - vdot(mcnew,tmpmc4,tmpmc3) + w
ELSE
CALL rwaxv2(n,inew-1,sm,um,tmpn4,tmpn4,tmpmc4,tmpmc3)
CALL rwaxv2(n,mcnew-inew+1,sm((iold-1)*n:),um((iold-1)*n:),tmpn4,tmpn4,tmpmc4(iold-1:),tmpmc3(iold-1:))
CALL scsum(mcc,-gamma,tmpmc3,tmpmc4,tmpmc4)
CALL lineq(mcnew,mcc,iold,tmpmat,tmpmc3,tmpmc4,ierr)
qr = qr - vdot(mcc-inew,tmpmc4(iold-1:),tmpmc3(iold-1:)) - vdot(inew-1,tmpmc4,tmpmc3) + w
END IF
END IF

pq = qr - prqr - prqr + pr
qqp = pq + prqr + rrq - pr - rrp + alfn
rrp = rrp - alfv
rrq = rrq + alfn - alfv

! Computation of multipliers lam1 and lam2

IF (pr > zero .AND. qr > zero) THEN
tmp1 = rrq/qr
tmp2 = prqr/qr
w = pr - prqr*tmp2
IF (w /= zero) THEN
lam1 = (tmp1*prqr - rrp)/w
lam2 = -tmp1 - lam1*tmp2
IF (lam1*(lam1 - one) < zero .AND. lam2*(lam1 + lam2 - one) < zero) GO TO 200
END IF
END IF

! Minimum on the boundary

100 CONTINUE

lam1 = zero
lam2 = zero
IF (alfn <= alfv) lam2 = one
IF (qr > zero) lam2 = MIN(one,MAX(zero,-rrq/qr))
w = (lam2*qr + rrq+rrq)*lam2
! w = (lam2*qr + 2.0_prec*rrq)*lam2
tmp1 = zero
IF (alfv >= zero) tmp1 = one
IF (pr > zero) tmp1 = MIN(one,MAX(zero,-rrp/pr))
! tmp2 = (tmp1*pr + 2.0_prec*rrp)*tmp1
tmp2 = (tmp1*pr + rrp+rrp)*tmp1

IF (tmp2 < w) THEN
w = tmp2
lam1 = tmp1
lam2 = zero
END IF

IF (qqp*(qqp - pq) < zero) THEN
IF (qr + rrq + rrq - qqp*qqp/pq < W) THEN
lam1 = qqp/pq
lam2 = one - lam1
END IF
END IF
    
200 CONTINUE

IF (lam1 == zero .AND. lam2*(lam2 - one) < zero .AND. -rrp - lam2*prqr > zero .AND. pr > zero) THEN
lam1 = MIN(one - lam2, (-rrp-lam2*prqr)/pr)
END IF

! Computation of the aggregate values

tmp1 = one - lam1 - lam2
DO i=1,n
ga(i-1)=lam1*gp(i-1)+lam2*g(i-1)+tmp1*ga(i-1)
END DO

alfv = lam2*alfn + tmp1*alfv

END SUBROUTINE aggsr1

!***********************************************************************
!*                                                                     *
!*     * SUBROUTINE lmbm *                                             *
!*                                                                     *
!*     Limited memory bundle subroutine for nonsmooth optimization.    *
!*                                                                     *
!***********************************************************************

SUBROUTINE lmbm(n,x,fu,func,lopa,vali,g,subgra,mc,nit,nfe,nge,iterm)

USE param, ONLY : small,large,zero,half,one
USE initializat, ONLY : na,mcu,iscale,tolf,tolf2,tolb,tolg,tolg2,xmax,eta,epsl,mtesf,mit,mfe
USE lmbm_sub, ONLY : vdot,xdiffy,copy,copy2,scsum

IMPLICIT NONE

! Arrays/scalar arguments/external functions coming from Python
INTEGER n
REAL(KIND=prec) fu,func,lopa,vali
REAL(KIND=prec) x(n)
REAL(KIND=prec) g(n)
EXTERNAL func
EXTERNAL vali
EXTERNAL subgra

! Scalar arguments
INTEGER, INTENT(INOUT) :: &
mc           ! Maximum number of stored corrections.
INTEGER, INTENT(OUT) :: & 
nit, &       ! Number of iterations.
nfe, &       ! Number of function evaluations.
nge, &       ! Number of subgradient evaluations.
iterm        ! Cause of termination: 1  - The problem has been solved with desired accuracy. 2  - Changes in function values < tolf in mtesf subsequent iterations. 3  - Changes in function value < tolf*small*MAX(|fu_k|,|fu_(k-1)|,1), where small is the smallest positive number such that 1.0 + small > 1.0. 4  - Number of function calls > mfe. 5  - Number of iterations > mit. 6  - Time limit exceeded. 7  - fu < tolb. -1  - Two consecutive restarts. -2  - Number of restarts > maximum number of restarts. -3  - Failure in function or subgradient calculations (assigned by the user). -4  - Failure in attaining the demanded accuracy. -5  - Invalid input parameters. -6  - Unspecified error.

! Local arrays
REAL(KIND=prec), DIMENSION(0:n-1) :: &
xo, &        ! Previous vector of variables.
gp, &        ! Previous subgradient of the ohjective function.
ga, &        ! Aggregate subgradient.
s, &         ! Difference of current and previous variables.
u, &         ! Difference of current and previous subgradients.
d, &         ! Direction vector.
tmpn1        ! Auxiliary array.
REAL(KIND=prec), DIMENSION(0:(n*(mcu+1))-1) :: &
sm,um        ! Matrises whose columns are stored differences of variables (sm) and subgradients (um).
REAL(KIND=prec), DIMENSION(0:((mcu+2)*(mcu+1)/2)-1) :: &
rm, &        ! Upper triangular matrix stored columnwise in the one-dimensional array.
umtum        ! Matrix whose columns are stored subgradient differences.
REAL(KIND=prec), DIMENSION(0:mcu) :: &
cm, &        ! Diagonal matrix.
smtgp, &     ! smtgp = trans(sm)*gp.
umtgp, &     ! umtgp = trans(um)*gp.
tmpmc1, &    ! Auxiliary array.
tmpmc2       ! Auxiliary array.
REAL(KIND=prec), DIMENSION(0:(n*na)-1) :: &
ax,ag        ! Matrix whose columns are bundle points and subgradients.
REAL(KIND=prec), DIMENSION(0:na-1) :: &
af           ! Vector of bundle values.

! Local scalars
REAL(KIND=prec) :: &
alfn, &      ! Locality measure.
alfv, &      ! Aggregate locality measure.
epsr, &      ! Line search parameter.
dnorm, &     ! Euclidean norm of the direction vector.
gnorm, &     ! Euclidean norm of the aggregate subgradient.
xnorm, &     ! Stopping criterion.
pxnorm, &    ! Previous stopping criterion.
p, &         ! Directional derivative.
tmax, &      ! Maximum stepsize.
t, &         ! Stepsize.
theta, &     ! Correction parameter for stepsize.
fo, &        ! Previous value of the objective.
gamma, &     ! Scaling parameter.
tl, &        ! Lower limit for t used in interpolation.
tu, &        ! Upper limit for t used in interpolation.
fl, &        ! Value of the objective function for t=tl.
fuv, &       ! Value of the objective function for t=tu.
epsa, &      ! Line search parameter.
epst, &      ! Line search parameter.
epslk, &     ! Line search parameter.
epsrk, &     ! Line search parameter.
thdnorm, &   ! Auxiliary scalar.
epsawk, &    ! Auxiliary scalar.
epstwk, &    ! Auxiliary scalar.
epslwk, &    ! Auxiliary scalar.
epsrwk       ! Auxiliary scalar.
INTEGER :: i, &
mcinit, &    ! Initial maximum number of stored corrections.
mcc, &       ! Current number of stored corrections.
inew, &      ! Index for the circular arrays.
ibfgs, &     ! Index of the type of BFGS update.
isr1, &      ! Index of the type of SR1 update.
iters, &     ! Null step indicator. 0  - Null step. 1  - Serious step.
nnk, &       ! Consecutive null steps counter.
ibun, &      ! Index for the circular arrays in bundle updating.
mal, &       ! Current size of the bundle.
ic, &        ! Correction indicator.
icn, &       ! Correction indicator for null steps.
iflag, &     ! Index for adaptive version.
neps, &      ! Number of consecutive equal stopping criterions.
ntesf, &     ! Number of tests on function decrease.
ncres, &     ! Number of restarts.
nres, &      ! Number of consecutive restarts.
nress, &     ! Number of consecutive restarts in case of tmax < tmin.
nout, &      ! Auxilary printout specification.
nin          ! Number of interpolations.

! Intrinsic functions
INTRINSIC ABS,MAX,SQRT

! Parameters
REAL(KIND=prec), PARAMETER :: &
fmin    = -large, &        ! Smallest acceptable value of the function.
tmin    = 1.0E-12_prec, &  ! Minimum stepsize.
lengthd = 1.0E+20_prec, &  ! Direction vector length.
rho     = 1.0E-12_prec     ! Correction parameter.
INTEGER, PARAMETER :: &
maxeps = 20, &             ! Maximum number of consecutive equal stopping criterions.
maxnrs = 2000, &           ! Maximum number of restarts.
maxint = 200               ! Maximum number of interpolations.

! Initialization

nout   = 0
nit    = 0
nfe    = 0
nge    = 0
ntesf  = 0
nres   = 1
ncres  = -1
nress  = 0
neps   = 0
iterm  = 0
iters  = 1
nnk    = 0
isr1   = 0
alfn   = zero
alfv   = zero
mcinit = mc
tmax   = xmax
xnorm  = large
epsr   = 0.25_prec+small

IF (epsl+epsl >= epsr) THEN
epsr = epsl+epsl + small
END IF

! Computation of the value and the subgradient of the objective function and the search direction for the first iteration

iterm = 0
fu = func(n,x)
iterm = 0
CALL subgra(n,x,g)
nfe = nfe + 1
nge = nge + 1

CALL restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk,alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)
CALL dobun(n,na,mal,x,g,fu,ax,ag,af,iters,ibun)

! Start of the iteration

iteration: DO

!PRINT*,nit

lopa = vali(n,x)

IF (lopa == 1.0) THEN
	iterm = -6
	EXIT iteration
END IF

! Computation of norms

IF (iters > 0) THEN
gnorm = vdot(n,g,g)
dnorm = SQRT(vdot(n,d,d))
p = vdot(n,g,d)
ELSE
gnorm = vdot(n,ga,ga)
dnorm = SQRT(vdot(n,d,d))
p = vdot(n,ga,d)
END IF

! Test on descent direction

IF (p+small*SQRT(gnorm)*dnorm <= zero) THEN
nres = 0
ELSE
nres = nres + 1
IF (nres == 1) THEN
CALL restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk,alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)
IF (ncres > maxnrs) THEN
nout = maxnrs
iterm = -2
EXIT iteration
END IF
CALL dobun(n,na,mal,x,g,fu,ax,ag,af,iters,ibun)
CYCLE iteration
END IF
nout = -1
iterm = -1
EXIT iteration
END IF

! Stopping criterion

nit = nit + 1
pxnorm = xnorm
xnorm = -p + 2.0_prec*alfv

! Tests for termination

IF (xnorm <= 1.0E+03_prec*tolg .AND. (mcc > 0 .OR. ibfgs == 2)) THEN
IF(half*gnorm + alfv <= tolg2 .AND. xnorm <= tolg) THEN
iterm = 1
EXIT iteration
END IF
IF (mc < mcu .AND. iflag == 0) THEN
mc = mc+1
iflag = 1
END IF
END IF

IF (nfe >= mfe) THEN
nout = mfe
iterm = 4
EXIT iteration
END IF

IF (nit >= mit) THEN
nout = mit
iterm = 5
EXIT iteration
END IF

IF (fu <= tolb) THEN
iterm = 7
EXIT iteration
END IF

IF (iters == 0) THEN
IF (ABS(xnorm - pxnorm) <= small) THEN
neps = neps + 1
IF (neps > maxeps) THEN
iterm = -4
EXIT iteration
END IF
ELSE
neps = 0
END IF
ELSE
neps = 0
END IF

! Correction

IF (-p < rho*gnorm .OR. icn == 1) THEN
xnorm = xnorm + rho*gnorm
dnorm = SQRT(dnorm*dnorm-2.0_prec*rho*p+rho*rho*gnorm)
IF (iters > 0) THEN
DO i=1,n
d(i-1) = d(i-1)-rho*g(i-1)
END DO
ELSE
DO i=1,n
d(i-1) = d(i-1)-rho*ga(i-1)
END DO
icn = 1
END IF
ic = 1
ELSE
ic = 0
END IF

! Preparation of line search

fo = fu

IF (iters > 0) THEN
CALL copy2(n,x,xo,g,gp)
END IF

IF (dnorm > zero) tmax = xmax/dnorm

IF (tmax > tmin) THEN
nress = 0
ELSE
nress = nress + 1
IF (nress == 1) THEN
CALL restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk,alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)
IF (ncres > maxnrs) THEN
nout = maxnrs
iterm = -2
EXIT iteration
END IF
CALL dobun(n,na,mal,x,g,fu,ax,ag,af,iters,ibun)
CYCLE iteration
END IF
iterm = -1
EXIT iteration
END IF

! Initial step size

CALL tinit(n,na,mal,x,af,ag,ax,ibun,d,fu,p,t,tmax,tmin,eta,iters,iterm)
IF (iterm /= 0) EXIT iteration

! Line search with directional derivatives which allows null steps

theta = one
IF (dnorm > lengthd) THEN
theta=lengthd/dnorm
END IF

! Old lls subroutine begins
!************************************************************************
!*                                                                      *
!*     * SUBROUTINE lls *                                               *
!*                                                                      *
!*     Special line search for limited memory bundle method             *
!*                                                                      *
!************************************************************************

! Initialization of old lls subroutine

nin = 0
epst = epsl+epsl
epsa = half*(epsr - epst)
thdnorm = theta*dnorm

tl = zero
tu = t
fl = fo

IF (theta < one) THEN
epst  = theta*epst
epsa  = theta*epsa
epslk = epsl
epsl  = theta*epsl
epsrk = epsr
epsr  = theta*epsr
END IF

epsawk   = epsa*xnorm
epslwk   = epsl*xnorm
epsrwk   = epsr*xnorm
epstwk   = epst*xnorm

! Function evalution at a new point

lmbm_iteration: DO

CALL scsum(n,theta*t,d,xo,x)
iterm = 0
fu = func(n,x,fu)
nfe = nfe + 1

IF (iterm /= 0) RETURN

! Null/descent step test (iters=0/1)

iters = 1
IF (fu <= fo - t*epstwk) THEN
tl = t
fl = fu
ELSE
tu = t
fuv = fu
END IF

! Additional interpolation

IF (fu > fo .AND. tu-tl >= tmin*0.1_prec .AND. nnk >= 1 .AND. nin < maxint) THEN
nin=nin+1
IF (tl == zero .AND. xnorm > zero) THEN
t = qint(tu,fl,fuv,xnorm,one-half/(one-epst))
ELSE
t = half*(tu+tl)
END IF
CYCLE lmbm_iteration
END IF

iterm = 0
CALL subgra(n,x,g)
nge = nge + 1

IF (iterm /= 0) RETURN

p = theta*vdot(n,g,d)
alfn = MAX(ABS(fo-fu+p*t),eta*(t*thdnorm)**2)

! Serious step

IF (fu <= fo - t*epslwk .AND. (t >= tmin .OR. alfn > epsawk)) EXIT lmbm_iteration

! Null step

IF (p-alfn >= -epsrwk .OR. tu-tl < tmin*0.1_prec .OR. nin >= maxint) THEN
iters = 0
EXIT lmbm_iteration
END IF

! Interpolation

nin=nin+1
IF (tl == zero .AND. xnorm > zero) THEN
t = qint(tu,fl,fuv,xnorm,one-half/(one-epst))
ELSE
t = half*(tu+tl)
END IF

END DO lmbm_iteration

IF (theta /= one) THEN
epsl = epslk
epsr = epsrk
END IF

! Old lls subroutine ends

IF (iterm /= 0) EXIT iteration

IF (tolf2 >= 0) THEN
IF (ABS(fo-fu) <= tolf2*small*MAX(ABS(fu),ABS(fo),one) .AND. iters == 1) THEN
iterm = 3
EXIT iteration
END IF
END IF

IF (ABS(fo-fu) <= tolf) THEN
ntesf = ntesf + 1
IF (ntesf >= mtesf .AND. iters == 1) THEN
iterm = 2
EXIT iteration
END IF
ELSE
ntesf = 0
END IF

! Bundle updating

CALL dobun(n,na,mal,x,g,fu,ax,ag,af,iters,ibun)

! Computation of variables difference

CALL xdiffy(n,x,xo,s)

! Computation of aggregate values and gradients difference

IF (iters == 0) THEN
nnk = nnk + 1
IF (nnk == 1) THEN
CALL copy(n,gp,tmpn1)
CALL xdiffy(n,g,gp,u)
CALL agbfgs(n,mc,mcc,inew,ibfgs,iflag,g,gp,ga,u,d,sm,um,rm,cm,umtum,alfn,alfv,gamma,ic,rho)
ELSE
CALL copy(n,ga,tmpn1)
CALL aggsr1(n,mc,mcc,inew,iflag,g,gp,ga,d,alfn,alfv,umtum,rm,gamma,smtgp,umtgp,tmpmc1,tmpmc2,sm,um,icn,rho)
CALL xdiffy(n,g,gp,u)
END IF
CALL copy(n,xo,x)
fu = fo
ELSE
IF (nnk /= 0) THEN
CALL copy(n,ga,tmpn1)
ELSE
CALL copy(n,gp,tmpn1)
END IF
nnk = 0
CALL xdiffy(n,g,gp,u)
END IF

! Serious step initialization

IF (iters > 0) THEN
icn = 0
alfn = zero
alfv = zero
END IF

! Direction finding

IF (iters > 0) THEN

! BFGS update and direction determination

CALL dlbfgs(n,mc,mcc,inew,ibfgs,iflag,d,g,gp,s,u,sm,um,rm,umtum,cm,smtgp,umtgp,gamma,tmpn1,iscale)
ELSE

! SR1 update and direction determination

CALL dlsr1(n,mc,mcc,inew,isr1,iflag,d,gp,ga,s,u,sm,um,rm,umtum,cm,smtgp,umtgp,gamma,tmpmc1,tmpmc2,tmpn1,nnk)
ibfgs=0
END IF

END DO iteration

!PRINT*,x

END SUBROUTINE lmbm

END MODULE lmbm_mod