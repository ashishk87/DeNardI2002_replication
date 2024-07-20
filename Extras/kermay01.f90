PROGRAM forlife

! US3y.f90 
! US calibration, based on swe3y.f90 (april 4th) code, with 
! mistake on agini corrected. 
! KERNEL CONVERGENCE ON BEQDIST
! lifespan uncertainty from AGE 60. estate taxation.
! 20 years old inherits parent's prod. when parent is 40.
! observes parent's prod. when the parent is 40.
! adult-child consumption equivalents (efc)

USE numerical_libraries

IMPLICIT NONE

INTEGER (KIND=1), PARAMETER :: tlength=5 ! period length
INTEGER (KIND=1), PARAMETER :: dage=5 ! parent-child age difference
INTEGER (KIND=1), PARAMETER :: dy=3 ! number of income states 
INTEGER (KIND=1), PARAMETER :: def=45 ! size of eflife
INTEGER (KIND=1), PARAMETER :: TL=9 ! TL=def/tlength # of working periods
INTEGER (KIND=1), PARAMETER :: dsur=65 ! alive from 20 to 65
INTEGER (KIND=1), PARAMETER :: T=14 !=dsur/tlength+1, max. # periods alive
INTEGER (KIND=1), PARAMETER :: TR=T-TL ! retirement length

! preferences
REAL, PARAMETER :: sig=1.5 ! risk aversion
DOUBLE PRECISION, PARAMETER :: bet=.9455**tlength ! disc. factor
! joy of bequeathing
REAL, PARAMETER :: phi1=-9.5
REAL, PARAMETER :: phi2=8.0

! income and inheritance processes
REAL, PARAMETER :: gam=0.83 ! income autoregr. coef.
REAL, PARAMETER :: hbet=0.67 ! inheritance autoregr. coeff.

DOUBLE PRECISION, PARAMETER :: gpop=1.012**tlength ! pop. growth rate
DOUBLE PRECISION, PARAMETER :: pi=3.1415926

!government & tech
REAL, PARAMETER :: r=1.06**tlength-1 ! exogenous interest rate
REAL, PARAMETER :: taua=.20          ! tax rate on capital income
REAL, PARAMETER :: taub=.10           ! tax rate on estate
REAL, PARAMETER :: exb=8.5   ! # years avg earn for estate taxes exempt
REAL, PARAMETER :: ggovt=.18         ! public expenditure
REAL, PARAMETER :: pens=.40          ! pensions: repl rate*average income
REAL, PARAMETER :: delt=1-.94**tlength ! depr rate, 6% yearly

! age-cons. requirement profile, correcting for kids by G.Weber
! efc=1+numch*[0; 0; 0; 0; 0; zeros(TL-1-5,1)]
!REAL, PARAMETER, DIMENSION(TL-1) :: efc=1+numch*(/0.0,0.06,0.07,0.08,0.0,0.0,0.0/)

! grids 
REAL, PARAMETER :: maxa=50.0   ! maximum level of capital
REAL, PARAMETER :: mina=0.0   ! minimum capital level. do not change
INTEGER, PARAMETER :: da=120 ! # of grid points for capital
! policy fns will be defined on a finer grid, then will interpolate
INTEGER, PARAMETER :: dar=120 ! # of grid points for retired agents
! consumption grid for working agents (we are substituting a_{t+1})
INTEGER, PARAMETER :: dc=dar ! # of grid points for consumption 
REAL, PARAMETER :: minc=0.0   ! min cons proportion. do not change
REAL, PARAMETER :: maxc=1.0   ! max cons proportion. do not change

! polynomial extrapolation (for val fn)
INTEGER (KIND=2), PARAMETER :: dagg=30 ! # of points added on the grid 
INTEGER (KIND=1), PARAMETER :: dpol=2  ! dpol-1 is the degree of the polynomial

! relaxation parameter for the iteration on the govt b.c.   
REAL, PARAMETER :: relax=.7   

! transition matrices M and MT
! # number of non-zero elements in the transition matrix M
! 20-25 yrs old: possible values today: (da*dy*dy)
! possible values tomorrow: (dy*2), (because of the interpolation and
! how we attribute the mass to 2 points)
INTEGER (KIND=4), PARAMETER :: nonzeroM=(da*dy*dy*3)*(dy*2)+(da*dy*dy)*(dy*2+dy*da)+ &
&	(da*dy*dy*4)*(dy*2+dy*da)+(da*dy*4)*(dy*2)+ &
&	(da*dy*dy)*da+(da*dy)*2+(da*4)*2
! for 20 and 25 yrs old
! 30 yrs old: next period might inherit
! 35 to 55 yrs old
! 60 yrs old
! 65 to 85 yrs old
! # number of non-zero elements in the transition matrix MT
! the 90 yrs old guys
INTEGER (KIND=2), PARAMETER :: nonzeroMT=da*2
! size of the invariant distribution invm
INTEGER (KIND=4), PARAMETER :: ninvm=4*dy*dy*da+5*dy*(dy+1)*da+5*da
! starting row index -1 for the retired
INTEGER (KIND=4), PARAMETER :: indret=4*dy*dy*da+5*dy*(dy+1)*da
! starting row index -1 for the people too young to have inherited
INTEGER (KIND=2), PARAMETER :: indnoh=dy*dy*da*4 
! define dy*dy*da 
INTEGER (KIND=2), PARAMETER :: dydyda=dy*dy*da
! define dy*(dy+1)*da 
INTEGER (KIND=2), PARAMETER :: dydy1da=dy*(dy+1)*da

! kkkkk  KERNEL STUFF !kkkk**** modified apr23 2000
INTEGER, PARAMETER :: numkernel=31 ! number of kernel functions
INTEGER, PARAMETER :: groupoint=da/(numkernel-1) ! number of points between
! the mean of 1 kernel and the mean of the next
! the first kernel is special (match probability at zero)
INTEGER, PARAMETER :: middlepoint=(groupoint+1)/2 
INTEGER, DIMENSION(numkernel-1) :: indmeank,upint,lowint
REAL, PARAMETER :: numstdev=0.5 !********************* check old stuff REAL
REAL, DIMENSION(numkernel-1) ::stdevk
REAL, DIMENSION(numkernel-1) :: meankernel
REAL, DIMENSION(da,numkernel) :: kernelfun ! value of the kernel
REAL, DIMENSION(numkernel*dy*(TL-3)) :: expkernf,expkernf2
REAL, DIMENSION(numkernel,dy,TL-3) :: kercoeff,kercoeff2 
REAL, DIMENSION(numkernel,numkernel) :: kertemp
REAL, DIMENSION(numkernel,da) :: kernelreg

REAL :: mny,vary,mnyh,varyh,eigval,epsi,epsilon,epsilon2
REAL :: area,aarea,gini,agini,gini1,gini2
REAL :: numch,aveinc,numret
REAL :: taul,taul2,gbal
REAL :: netearn,weal
REAL :: tempreal
REAL :: avgk,beqtax,wgini,wginic,fpoor
REAL :: trwealth,gnp,k2gnp,kshare,kotsum

! age-cons. requirement profile, correcting for kids by G.Weber
! efc=1+numch*[0; 0; 0; 0; 0; zeros(TL-1-5,1)]
REAL, DIMENSION(TL-1) :: efc

REAL, DIMENSION(def) :: eflife
REAL, DIMENSION(TL) :: eff,endeff
REAL, DIMENSION(dsur) :: surlife
REAL, DIMENSION(T) :: sur,sur1,pfrac
REAL, DIMENSION(T+4) :: csur
REAL, DIMENSION(da) :: a,anet
REAL, DIMENSION(dar) :: ar,phi
REAL, DIMENSION(dy,dy) :: p,Qy,Qyh,Qher,jt4020 
REAL, DIMENSION(dy) :: invdy,invdyh,cvary,cvaryh,grid,wt
REAL, DIMENSION(dy) :: y,cmny,cmnyh,eigvec
REAL, DIMENSION(dy) :: invher,inv40
REAL, DIMENSION(dy+1) :: yy,dist,cumdist,cumearn
REAL, DIMENSION(dy+1) :: dist1,dist2,cumdist1,cumdist2,cumearn1,cumearn2
REAL, DIMENSION(TL*dy) :: eq,eminv,eminv1,aearn,aearn1
REAL, DIMENSION(TL*dy+1) :: aearn0,adist,cumadist,cumaearn
REAL, DIMENSION(da) :: vti
REAL, DIMENSION(da,dy) :: vtai
REAL, DIMENSION(TL,dy) :: ageyn
REAL, DIMENSION(TL*dy,TL*dy) :: EM
REAL, DIMENSION(dpol,dpol) :: xinv
REAL, DIMENSION(da,TR+1) :: vr
REAL, DIMENSION(da,dy,dy+1,TL) :: vt,sopt !,copt
REAL, DIMENSION(dar,da) :: cr,ur
REAL, DIMENSION(dar) :: vrr,valr
REAL, DIMENSION(dc+1) :: cc
REAL, DIMENSION(dc) :: c,u,st,at,vst,phit,vat,evat
REAL, DIMENSION(nonzeroM) :: valM
REAL, DIMENSION(nonzeroMT) :: valMT
REAL, DIMENSION(ninvm) :: invm,invm1
REAL, DIMENSION(dar) :: rest
REAL, DIMENSION(da,TR) :: restr
REAL, DIMENSION(da,dy,TL-4) :: restt
REAL, DIMENSION(da,dy,dy,TL) :: restt2
REAL, DIMENSION(da) :: restt3
REAL, DIMENSION(dy*dy*da) :: qinit,q
REAL, DIMENSION(da,dy,TL-3) :: beqdist,beqdist2
REAL, DIMENSION(da*dy*(dy+1)) :: temprealv,temprealv2
REAL, DIMENSION(dy) :: temprealv4
REAL, DIMENSION(da) :: temprealda,temprealda2,temprealda4
REAL, DIMENSION(da*dy*(dy+1),dy) :: stbeq,stbeq2
REAL, DIMENSION(da,dy) :: stbeqret,stbeqret2
REAL, DIMENSION(da) :: dista,cumdista,cumwealth
REAL, DIMENSION(da) :: distac,cumdistac,cumwealthc
REAL, DIMENSION(da,2) :: topwel,topwelc
REAL, DIMENSION(da,T) :: agew
REAL, DIMENSION(da) :: agewT
REAL, DIMENSION(7) :: sp, squant

DOUBLE PRECISION :: sigz,sigh,scal

INTEGER (KIND=2) :: j,j0,j1,j2,j4,j5,j6,j7,k,jjj,jj,dyy,OpenStatus
INTEGER (KIND=4) :: i,tempint1,tempint2
INTEGER (KIND=4) :: numsort,counter
INTEGER (KIND=2) :: counter2,counterbeq
INTEGER, DIMENSION(TL*dy) :: iperm
INTEGER (KIND=2), DIMENSION(1) :: imax
INTEGER (KIND=2), DIMENSION(da,TR) :: polr,indexr
INTEGER (KIND=2), DIMENSION(da,dy,dy+1,TL) :: polt
INTEGER (KIND=2), DIMENSION(da,dy,TL-4) :: indext
INTEGER (KIND=2), DIMENSION(da,dy,dy,TL) :: indext2
INTEGER (KIND=2), DIMENSION(da) :: indext3

INTEGER (KIND=4), DIMENSION(nonzeroM) :: colM,rowM
INTEGER (KIND=2), DIMENSION(nonzeroMT) :: colMT,rowMT
INTEGER (KIND=2), DIMENSION(dar) :: inde

INTEGER :: IHOUR, IMIN, ISEC

COMPLEX, DIMENSION(dy) :: EVAL
COMPLEX, DIMENSION(dy,dy) :: EVEC

INTERFACE

SUBROUTINE cumsumv(AA,BB)
! matlab couterpart of cumsum.f90
! NOW THE SPECIAL CASE WHERE AA IS A VECTOR
IMPLICIT NONE
REAL,INTENT(IN),DIMENSION(:) :: AA
REAL,INTENT(OUT),DIMENSION(:) :: BB
INTEGER i
END SUBROUTINE cumsumv

SUBROUTINE  interplin(l,x,y,n,z,dpol,xinv,v)
! linear interpolation. extrapolates if z
! falls outside the grid of x
IMPLICIT NONE
INTEGER, INTENT(IN) :: l
INTEGER, INTENT(IN) :: n
INTEGER (KIND=1), INTENT(IN) :: dpol
REAL, DIMENSION(l), INTENT(IN) :: x
REAL, DIMENSION(l), INTENT(IN) :: y
REAL, DIMENSION(n), INTENT(IN) :: z
REAL, DIMENSION(dpol,dpol), INTENT(in) :: xinv
REAL, DIMENSION(n), INTENT(OUT) :: v
INTEGER :: i,ind
INTEGER, DIMENSION(1) :: k
REAL :: diff
REAL, DIMENSION(dpol) :: zpoly
END SUBROUTINE interplin

SUBROUTINE linspace (xmin,xmax,npoints,lspace)
REAL, INTENT(IN) :: xmin,xmax
INTEGER, INTENT(IN) :: npoints
REAL, DIMENSION(:), INTENT(OUT) :: lspace
END SUBROUTINE linspace

END INTERFACE

! initialize tax rate on labor income
taul=0.2823769

sigz=SQRT(.412) ! income std dev
sigh=SQRT(.42) ! inher. std dev
scal=(sigz+sigh)/2 !scaling for hermite quadr.

! discretize income and inherit. process using HERMITE QUADRATURE
! want them on the same grid y, take some combi of the condit. distr.
! (normal with mean 0 and variance a linear combination of the two
! variances) define a scaling variable for the variance of the weigthing fn
!!INCOME PROCESS y, Qy
! load grid points and weights from a code that used the nag routines
! $markov.dat contains already sorted grid and wt data. write sort routine later$.
OPEN(1,file='markov.dat',status='old',form='formatted',IOSTAT=OpenStatus)
READ (1,*) grid, wt
IF(OpenStatus>0) STOP 'cannot open markov.dat'
CLOSE(1)

! renormalize since the routine does not take into account pi
wt=wt/SQRT(2*pi)

! construct transition matrix for the markov process
p=0.0 ! initialization
DO i=1,dy
    DO j=1,dy
        p(i,j)=(scal/sigz)*EXP(-(grid(j)-gam*grid(i))**2*&
        scal**2/(2*sigz**2))*wt(j)*EXP(grid(j)**2/2)
    END DO ! for j
    p(i,:)=p(i,:)/SUM(p(i,:),dim=1)
END DO ! for i

! map z_t**star to z_t
grid=grid*scal
! conditional mean and variance
cmny=MATMUL(p,grid)
cvary=MATMUL(p,grid**2)-cmny**2

! find eigenvalues and vectors of P'
CALL EVCRG (dy, TRANSPOSE(p), dy, EVAL, EVEC, dy)
IF((REAL(EVAL(1))-1>1e-5) .OR. (REAL(EVAL(1)-1)>1e-5)) STOP 'problems with Qy'
! EVCRG provides the eigenvalues in decreasing order of magnitude
! If the first one is not close enough to one, or the second one
! is also close to one, stop: missing unit root or multiple ones
eigval=REAL(EVAL(1))
eigvec=REAL(EVEC(:,1))

invdy=eigvec/SUM(eigvec)

mny=DOT_PRODUCT(grid,invdy)
vary=DOT_PRODUCT(grid**2,invdy)-mny**2
y=exp(grid)
Qy=p

! inheritance process y Qyh
! undo the tranfo on grid to return to the original one
grid=grid/scal
p=0 !initialization
DO i=1,dy
    DO j=1,dy
        p(i,j)=(scal/sigh)*EXP(-(grid(j)-hbet*grid(i))**2*&
        scal**2/(2*sigh**2))*wt(j)*EXP(grid(j)**2/2)
    END DO ! for j
    p(i,:)=p(i,:)/sum(p(i,:),dim=1)
END DO ! for i
! map z_t**star to z_t
grid=grid*scal
! conditional mean and variance
cmnyh=MATMUL(p,grid)
cvaryh=MATMUL(p,grid**2)-cmnyh**2
! find eigenvalues and vectors of P'
CALL EVCRG (dy, TRANSPOSE(p), dy, EVAL, EVEC, dy)
IF((REAL(EVAL(1))-1>1e-5) .OR. (REAL(EVAL(1)-1)>1e-5)) STOP 'problems with Qy'
! EVCRG provides the eigenvalues in decreasing order of magnitude
! If the first one is not close enough to one, or the second one
! is also close to one, stop: missing unit root or multiple ones
eigval=REAL(EVAL(1))
eigvec=REAL(EVEC(:,1))
invdyh=eigvec/SUM(eigvec)

mnyh=DOT_PRODUCT(grid,invdyh)
varyh=DOT_PRODUCT(grid**2,invdyh)-mnyh**2
Qyh=p

! matrix Qyh: inheritance of y: 40 years old parent-20 years old child
! want trans.  matrix Qher:
! 20 years old parent- her future 20 years old child
! (because want the invariant distr. at birth of these char.in pop.)
! have: c_20'=p_40'*Qyh, rewrite as
! c_20'=p_35'*Qy*Qyh, as
! and p_40'=p_35'*Qy.
Qher=MATMUL(MATMUL(Qy,MATMUL(Qy,MATMUL(Qy,Qy))),Qyh)
!compute invariant distribution
CALL EVCRG (dy, TRANSPOSE(Qher), dy, EVAL, EVEC, dy)
IF((REAL(EVAL(1))-1>1e-5) .OR. (REAL(EVAL(1)-1)>1e-5)) STOP 'problems with Qy'
! EVCRG provides the eigenvalues in decreasing order of magnitude
! If the first one is not close enough to one, or the second one
! is also close to one, stop: missing unit root or multiple ones
eigval=REAL(EVAL(1))
eigvec=REAL(EVEC(:,1))
invher=eigvec/SUM(eigvec)


! distribution of the characteristics at 40 years of age
inv40=MATMUL(invher,MATMUL(Qy,MATMUL(Qy,MATMUL(Qy,Qy))))
! joint distr. on the parent's charact. at 40 and the kid's ones at 20
! prob(y40,y)=prob(y|y40)*prob(y40), i.e.
! prob(y40=yj,y=yi)=prob(y=yi|y40=yj)*prob(y40=yj)
DO i=1,dy
    DO j=1,dy
        jt4020(i,j)=Qyh(i,j)*inv40(i)
    END DO
END DO

! construct age-efficiency profile for workers
! from 20 to 65
OPEN(2,file='eflife.dat',status='old',form='formatted')
READ (2,*) eflife
CLOSE(2)
! transform yearly age-efficiency profile according to our period length
DO i=1,def,tlength
    eff((i-1)/tlength+1)=sum(eflife(i:i+tlength-1))
END DO
eff=eff*TL/SUM(eff)

! construct conditional survival probabilities
OPEN(3,file='surlife.dat',status='old',form='formatted')
READ (3,*) surlife
CLOSE(3)
! surlife.dat contains conditional surv. prob.from 20 to 90
DO i=1,dsur,tlength
    sur((i-1)/tlength+1)=PRODUCT(surlife(i:i+tlength-1))
END DO
! people do not die before age 60
sur(1:8)=1.0
sur(T)=0.0
! define the cumulative survival probabilities (up to age 20 it is 1)
csur(1:4)=1.0
csur(5:T+4)=sur
DO i=5,T+4
    csur(i)=PRODUCT(sur(1:i-4))
END DO

! statistics for the productivity markov process
! compute the gini index for productivity.
! the number of people for each income level is given by the invariant
! distribution of the markov process we are considering.
dist(1)=0
dist(2:dy+1)=invdy
yy(1)=0
yy(2:dy+1)=y
dyy=SIZE(yy)
! lorenz curve plot (cumulative productivity plotted against cumulative
! number of agents both normalized to one)
CALL cumsumv(dist,cumdist)
CALL cumsumv(yy*dist,cumearn)
cumearn=cumearn/SUM(yy*dist)
! compute the gini index
! the area under the lorenz curve is composed by trapeziums:
!          |c
! b |      |
!   |      |
! a |      |d
! whose area is (cd+ba)*ad/2.
! the formula for the gini index given by the ratio of the
! area between the 45 degree line and the the lorenz curve divided by
! the total area under the 45 degree line:
! (1/2 - \sum_{i=1}**{dyy-1}(e(i+1)+e(i))*(cumdist(i+1)-cumdist(i))/2)/(1/2)=
! 1-\sum_{i=1}**{dyy-1}(e(i+1)+e(i))*(cumdist(i+1)-cumdist(i))
area=0
DO i=1,dyy-1
   area=area+(cumearn(i+1)+cumearn(i))*(cumdist(i+1)-cumdist(i))
END DO
gini=1-area
! we also want to compute the gini index for earnings taking into account
! the age-productivity profile and how the inheritance process interacts
! with the earnings distribution. call it "adjusted gini". we thus need to
! perform some computation. form the transition matrix:
!                                    age 1        age 2    .... age TL
!                              prod1....prod dy  .....         ........
! age 1 productivity level 1
!        productivity level 2
!        ....................
!        productivity level dy
! age 2 productivity level 1
!        productivity level 2
!        ....................
!        productivity level dy
! ...........................
! age TL productivity level 1
!        ....................
!        productivity level dy
EM=0
DO i=1,TL-1 ! today's age
    DO j=1,dy ! today's productivity
        EM((i-1)*dy+j,dy+(i-1)*dy+1:dy+i*dy)=sur(i)*Qy(j,:)
    END DO
END DO
eq=0
eq(1:dy)=invher
! now compute the invariant distribution
epsi=1
eminv=eq
DO WHILE (epsi>1e-08)
   eminv1=eminv
   eminv=MATMUL(TRANSPOSE(EM),eminv)/gpop+eq
   epsi=MAXVAL(ABS(eminv1-eminv))
END DO !for DO WHILE (epsi>1e-08)
! renormalize the invariant distribution
eminv=eminv/SUM(eminv)

! now the we have the invariant distribution over age and earnings
! we can compute the "adjusted" gini index that takes into account
! also the age efficiency profile.
! form a matrix with age on the rows and productivity on the columns
! containing the number of people in each cell in order to compute
! the age-efficiency profile IMPLIED by our income process (if we normalize
! the mean of the log income process wo do not have the right to other
! normalizations)
ageyn=TRANSPOSE(RESHAPE(SOURCE=eminv,SHAPE=(/dy,TL/)))
! transform the matrix in terms of fractions of people of a given age
DO i=1,TL
    ageyn(i,:)=ageyn(i,:)/sum(ageyn(i,:))
END DO
! compute the average age-productivity profile implied by our income process
endeff=MATMUL(ageyn,y)
! readjust the income process so that the implied age-efficiency profile
! is the one given in the eflife matrix
eff=eff/endeff
! adjusted gini index:construct y*eff for each possible
! y and eff and sort them in ascending order, keeping track of the
! corresponding number of people in that cell
DO i=1,TL
    DO j=1,dy
        aearn((i-1)*dy+j)=eff(i)*y(j)
    END DO
END DO
numsort=TL*dy
!create vector of indexes to be permuted for the sorting fn
iperm(1)=1
DO i=2,numsort
    iperm(i)=iperm(i-1)+1
END DO

! sort vector
CALL SVRGP(numsort,aearn,aearn1,iperm)

aearn0(1)=0.0
aearn0(2:dy*TL+1)=aearn1(1:dy*TL)
adist(1)=0.0
adist(2:dy*TL+1)=eminv(iperm(1:dy*TL)) 
! ********* corretto errore in indice di iperm

CALL cumsumv(adist,cumadist) !cumulative # of agents
CALL cumsumv(aearn0*adist/SUM(aearn0*adist),cumaearn)! cum. adj. productivity

! compute the gini index
aarea=0.0
DO i=1,TL*dy
    aarea=aarea+(cumaearn(i+1)+cumaearn(i))*(cumadist(i+1)-cumadist(i))
END DO
agini=1-aarea

! compute now the gini index for productivity for just born people and
! for people who are TL years old
dist1(1)=0.0
dist1(2:dy+1)=eminv(1:dy)

dist2(1)=0
dist2(2:dy+1)=eminv((TL-1)*dy+1:TL*dy)

yy(1)=0
yy(2:dy+1)=y(1:dy)

! lorenz curve plot (cumulative productivity plotted against cumulative
! number of agents both normalized to one)
CALL cumsumv(dist1/SUM(dist1),cumdist1) !cum. # of agents
CALL cumsumv(dist2/SUM(dist2),cumdist2)
CALL cumsumv(yy*dist1/SUM(yy*dist1),cumearn1) ! cum. productivity
CALL cumsumv(yy*dist2/SUM(yy*dist2),cumearn2)
! compute the gini index for new agents
area=0
DO i=1,dyy-1
    area=area+(cumearn1(i+1)+cumearn1(i))*(cumdist1(i+1)-cumdist1(i))
END DO
gini1=1-area
! compute the gini index for Tl years old agents
area=0
DO i=1,dyy-1
   area=area+(cumearn2(i+1)+cumearn2(i))*(cumdist2(i+1)-cumdist2(i))
END DO
gini2=1-area

! compute per capita number of children
! we have: newborns_{t+1}=gpop*newborns_t (1); pop. law of motion
! newborns_t=numch*parents_t (2); numch is the per capita number of children
! (2): total number of children
! parents_t=newborns_{t-dage}*sur(1) (3)
! (3): the number of parents (at 25) equals the number of newborns 5 periods
! before times their probability of survival (one until they are 20 and sur(1)
! between 20 and 25).
! from (2) and (3): newborns_t=numch*newborns_{t-dage}*sur(1)  (4)
! from (1): newborns_t=gpop^(dage)*newborns_{t-dage}           (5)
! equate (4) and (5), simplify and get:
numch=gpop**dage/sur(1)
!efc=1+numch*(/0.5,0.06,0.07,0.08,0.0,0.0,0.0,0.0/)
efc=1.0

sur1(1)=1.0
sur1(2:T)=sur(1:T-1)
DO i=1,T
    pfrac(i)=PRODUCT((sur1(1:i)/gpop))
END DO
pfrac=pfrac/sum(pfrac)

! aggregate labor income (which is exogenous).
aveinc=SUM(pfrac(1:TL)*eff*endeff)
numret=SUM(pfrac(TL+1:T))

! form assets grid for which the value fn will be defined
! how the grid is formed is crucial for later transformations
! to map one grid to another
!a=real( (/ (i,i=0,da-1) /) )/real(da-1)*(maxa-mina)+mina
CALL linspace(mina,SQRT(maxa),da,a)
a=a**2
! define RECEIVED ASSETS NET OF ESTATE TAXES PER PERSON
! later on we will also divide it by csur(j0+4), j0 being age
WHERE (a .LE. exb)
    anet=a/numch
ELSEWHERE
    anet=(a-taub*(a-exb))/numch
END WHERE

! kkkkkk KERNEL STUFF !kkkk**** modified apr23 2000
! find location of mean and each upper and lower bound
! ofreference intervals
DO i=1,numkernel-1
	indmeank(i)=middlepoint+(i-1)*groupoint
	upint(i)=groupoint+(i-1)*groupoint
	lowint(i)=1+(i-1)*groupoint
END DO ! DO i
meankernel=a(indmeank)
stdevk=(a(upint)-a(lowint))/(2*numstdev)
kernelfun=0.0 !initialize 
kernelfun(1,1)=1.0 ! first kernel is degenerate on zero
! approximate the normal distributions that constitute the 
! kernel on our grid
DO j=1,numkernel-1
	kernelfun(1,j+1)=((a(2)-meankernel(j))/(a(2)-a(1)))* &
		&	(ANORDF((a(2)-meankernel(j))/stdevk(j))- &
		&	ANORDF((a(1)-meankernel(j))/stdevk(j)))+ &
		&	(stdevk(j)/(SQRT(2.0*pi)*(a(2)-a(1))))* &
		&	(EXP(-(a(2)-meankernel(j))**2/(2.0*stdevk(j)**2))- &
		&	EXP(-(a(1)-meankernel(j))**2/(2.0*stdevk(j)**2)))
	DO i=2,da-1
		kernelfun(i,j+1)=((meankernel(j)-a(i-1))/(a(i)-a(i-1)))* &
		&	(ANORDF((a(i)-meankernel(j))/stdevk(j))- &
		&	ANORDF((a(i-1)-meankernel(j))/stdevk(j)))- &
		&	(stdevk(j)/(SQRT(2.0*pi)*(a(i)-a(i-1)))) * &
		&	(EXP(-(a(i)-meankernel(j))**2/(2.0*stdevk(j)**2))- &
		&	EXP(-(a(i-1)-meankernel(j))**2/(2.0*stdevk(j)**2)))+ &
		&	((a(i+1)-meankernel(j))/(a(i+1)-a(i)))* &
		&	(ANORDF((a(i+1)-meankernel(j))/stdevk(j))- &
		&	ANORDF((a(i)-meankernel(j))/stdevk(j)))+ &
		&	(stdevk(j)/(SQRT(2.0*pi)*(a(i+1)-a(i))))* &
		&	(EXP(-(a(i+1)-meankernel(j))**2/(2.0*stdevk(j)**2))- &
		&	EXP(-(a(i)-meankernel(j))**2/(2.0*stdevk(j)**2)))		
	END DO ! DO i
	kernelfun(da,j+1)=((meankernel(j)-a(da-1))/(a(da)-a(da-1)))* &
		&	(ANORDF((a(da)-meankernel(j))/stdevk(j))- &
		&	ANORDF((a(da-1)-meankernel(j))/stdevk(j)))- &
		&	(stdevk(j)/(SQRT(2.0*pi)*(a(da)-a(da-1)))) * &
		&	(EXP(-(a(da)-meankernel(j))**2/(2.0*stdevk(j)**2))- &
		&	EXP(-(a(da-1)-meankernel(j))**2/(2.0*stdevk(j)**2)))
END DO ! DO j
!kkkkkkkkkkkkkkk
! save stuff in matrices to check kernel functions
OPEN(UNIT=20, FILE= "USkerfun19.dat", STATUS="REPLACE", ACTION= "WRITE", &
&   POSITION= "REWIND", IOSTAT= OpenStatus)
IF (OpenStatus>0) STOP "***** Cannot open file *****"
DO k=1,da
	DO j0=1,numkernel
		WRITE(20,*) kernelfun(k,j0)					
	END DO ! DO j
END DO ! DO k
CLOSE(20)
! kkkkkk modified May 6th
kertemp=MATMUL(TRANSPOSE(kernelfun),kernelfun)           
CALL LINRG (numkernel, kertemp, numkernel, kertemp, numkernel)
kernelreg=MATMUL(kertemp,TRANSPOSE(kernelfun))
		
! will need to take care of the possibility that when interpolating
! the value fn later on we fall out of the range of values we
! have for the value fn itself.
! polynomial extrapolation, how do we do it?
! consider the quadratic case as an example
! the equation is y=a+bx+cx^2, where y=[y_1 y_2 y_3]' and analogously
! for x. we can rewrite it in matrix form as
! |y_1|   | 1 x_1  x_1^2 |   | a |
! |y_2| = | 1 x_2  x_2^2 |   | b |
! |y_3|   | 1 x_3  x_3^2 |   | c |
! invert it to get the vector of coefficients:
!|a |   | 1 x_1  x_1^2 |-1   | y_1 |
!|b | = | 1 x_2  x_2^2 |     | y_2 |
!|c |   | 1 x_3  x_3^2 |     | y_3 |
! now use:
! |y_4|   | 1 x_4  x_5^2 |   | a |
! |y_5| = | 1 x_4  x_5^2 |   | b |
! |y_6|   | 1 x_4  x_5^2 |   | c |
! substituiting the expression for the coefficients derived above, obtain:
! |y_4|   | 1 x_4  x_5^2 |   | 1 x_1  x_1^2 |-1   | y_1 |
! |y_5| = | 1 x_4  x_5^2 |   | 1 x_2  x_2^2 |     | y_2 |
! |y_6|   | 1 x_4  x_5^2 |   | 1 x_3  x_4^2 |     | y_3 |
!         call this temp2    call this xinv
! temp2*xinv is what we will use when we will have to extrapolate the
! value fn later on (in the interpolation/extrapolation routine)
! when we will have the relevant y's.
! CONSTRUCT xinv ONLY HERE
DO i=1,dpol ! row index
    DO j=1,dpol ! column index
        xinv(i,j)=(a(da-dpol+i))**(j-1)
    END DO ! for j
END DO ! for i
! invert xinv
CALL LINRG (dpol, xinv, dpol, xinv, dpol)

! grid on tomorrow's assets the agent can choose.
! It is bigger than the grid for a, for which the value fn
! is defined to save on computation time while not constraining too much
! the agent's choice on a narrow grid. Will interpolate the value fn
! on ar later on (using the value fn defined on a) to compute
! the agent's maximization problem for each period of time
! retired agents
CALL linspace(mina,SQRT(maxa),dar,ar)
ar=ar**2

! working agents: need a grid for consumption
!(we are subtituting a_{t+1}) which will be different for different
! asset level, productivity level and age
! c = PROPORTION of this period's income consumed
CALL linspace(minc,maxc,dc+1,cc) !! this creates a linear grid between minc and maxc with dc+1 points

! grid for bequests: same one as for assets (a)
! initialize a distribution for the bequests which depends
! on agent's age (the agent can inherit up to the first year of
! retirement) and on her parent's productivity (at 40 years old).
! SINCE PEOPLE DO NOT DIE BEFORE 60 IN THIS VERSION OF THE CODE,
! PEOPLE YOUNGER THAN 40 DO NOT INHERIT, (AT THE EARLIEST THEY INHERIT
! AT 40; IN THE OLD VERSION OF THE CODE THEY COULD FIRST INHERIT AT 20,
! THEREFORE THE VECTOR FOR BEQDIST WILL BE FOUR PERIOD SHORTER: IT WAS
! TL+1 LONG, IT IS NOW TL-3).
! BEQDIST CONTAINS THE BEQUEST DISTRIBUTION GROSS OF ESTATE TAXES.
! initializition: they do not expect to receive any bequest
beqdist=0
beqdist(1,:,:)=1
!load beqdist

! joy of bequeathing
! assume the child's preferences are the parent's ones.
! Assume that the child will consume 1 each period of her life
! plus the bequest she receives in 4 equal amounts, starting
! the moment she gets them. Her utility:
! beta^2*(1+a/4)^(1-sigma)/(1-sigma)+beta^3*(1+a/4)^(1-sigma)/(1-sigma)+...
! beta^4*(1+a/4)^(1-sigma)/(1-sigma)+beta^5*(1+a/4)^(1-sigma)/(1-sigma)-...
! beta^2/(1-sigma)-beta^3/(1-sigma)-beta^4/(1-sigma)-beta^5/(1-sigma)
! which we can rewrite as:
! 1/(1-sigma)*((1+a/4)^(1-sigma)-1)*beta^2*(1+beta+beta^2+beta^3)=
! 1/(1-sigma)*((1+a/4)^(1-sigma)-1)*beta^2*(1-beta^4)/(1-beta)=
! k1*((1+a/k2)^(1-sigma)-1)
! phi1=3*(1/(1-sig))*bet.^2.*(1-bet.^4)./(1-bet)
! NOW ESTATES ABOVE EXB ARE TAXED AT TAUB. UTILITY FROM WARM GLOW IS FROM
! THE NET AMOUNT
phi=((1+(ar-taub*MAX(ar-exb,0*ar))/phi2)**(1-sig)-1)*phi1

! value and policy function initializations
! retired agents
vr=0
! set the last period's value fn (the agent will be dead for
! sure)to the value of the "joy of bequeathing" for that period.
! TAKE INTO ACCOUNT AFTER TAX BEQUEST VALUE
vr(:,TR+1)=((1+(a-taub*MAX(a-exb,0*a))/phi2)**(1-sig)-1)*phi1
! initialization of the retired agent's policy function
polr=0
! working people
vt=0
polt=0 ! policy function index
!copt=0 ! optimal consumption policy function
sopt=0 ! optimal savings policy function. When the guy can inherit,
       ! there is a probability distribution on next period's assets,
       ! linked to the relevant probability distribution on bequests.

! initialize q: the vector for the newborns
! NOTE THAT SINCE NO PARENT DIES BEFORE AGE 60, NO CHILD INHERITS BEFORE
! AGE 40. HENCE ALL NEWBORNS START OFF FROM ZERO WEALTH
! define its size
DO j=1,dy ! parent's productivity
    DO k=1,dy ! child's productivity
        qinit((j-1)*da*dy+(k-1)*da+1)=jt4020(j,k)
    END DO ! DO k
END DO ! DO j

! initialize invariant distribution
invm=0

! RETIREES: their VALUE FN does not depend on taul and from
! received bequests=does not need to stay in the loops!
! the same is true for filling up MT and their section of M
! SOLVE THE VALUE FUNTION FOR THE RETIREES
! solving backward: retired agents
! take care of negative consumption and define utility from
! different choices of tomorrow's k
DO i=1,da ! today's assets
    cr(:,i)=(1+r*(1-taua))*a(i)+pens-ar(:)
END DO ! for i
!take care of negative consumption
WHERE (cr>0) !take care of negative consumption
    ur=(cr**(1-sig))/(1-sig)
ELSEWHERE
    ur=-1e+10
END WHERE
! RETIREES iterate on the value fn, starting from the last period
DO j=TR,1,-1
    ! interpolate value fn
    CALL interplin(da,a,vr(:,j+1),dar,ar,dpol,xinv,vrr)
    DO i=1,da ! evaluate value fn given today's assets
        valr=ur(:,i)+sur(TL+j)*bet*vrr+(1-sur(TL+j))*phi
        imax=MAXLOC(valr)
        polr(i,j)=imax(1) ! optimal asset choice: index tells the savings level at which value fn is maximized
        vr(i,j)=valr(imax(1))
    END DO
END DO ! for j

! COMPUTE MT AND FILL THE SECTION OF M FOR THE RETIREES
! STRUCTURE OF M: age-parents' productivity-agent's productivity-assets
! the first variable is the agent's age.
! for ages from TL+1 to T the agent is retired, therefore for each age
! we only distinguish by assets (parent's productivity is not relevant because
! they do not inherit any more and current productivity is not relevant
! because they do not work any more but just receive pensions).
!initialize
rowM=0
colM=0
valM=0

! retired people
! to compute M we need first to map the optimal asset choice
! given by the policy fns on the same grid for the assets on which
! keep track of the value function (which is coarser to save on memory
! space).
! map the ar indexes to the a indexes: note that this
! transformation and the analogous ones that will follow crucially
! depend on how the grids
! were formed. How do we do it?
! With a linear grid we have: index=a*x+b, hence x=c*index+d.
! We do not have a linear grid, and more in general:
! index=f(x), x=g(index)
! look how we constructed our grid: linspace(0,sqrt(maxa),da) and then
! square it. Lets' invert the mapping first using linearity and then the
! nonlinear transformation.
! The linear mapping is y=h(index)=h_0+h_1*index
! have y(0)=0, y(da)=sqrt(maxa), hence h_0+h_1=0, h_0+h_1*da=sqrt(maxa)
! solving we get h_0=-sqrt(maxa)/(da-1), h_1=sqrt(maxa)/(da-1).
! The nonlinear mapping is x=y^2, where y is defined above, hence
! x=y^2=(-sqrt(maxa)/(da-1) + sqrt(maxa)/(da-1) index)^2
! invert the mapping and take the floor to get index in terms of x:
inde=FLOOR((da-1)*SQRT(ar)/SQRT(maxa)+1)
!take care of the endpoint
IF (ar(dar)==a(da)) THEN
    inde(dar)=da-1
END IF
! find the distance from the smallest closer element to the element
! in sopt that we are considering.
! This is because we will approximate the transition matrix and
! the corresponding invariant distribution using the grid a and
! giving weight to the two level of assets between which each optimal
! choice is included. These weight will be the distances form those
! two point (the total distance is obviously normalized to one since
! we are talking about weights).
rest=(ar-a(inde))/(a(inde+1)-a(inde))
! take the optimal policy function for retired people (which lies on the
! ar grid) and convert each of its points in the corresponding smaller,
! closer value of assets on the grid a (indexr) and its distance from it
! (restr)
DO k=1,TR
    indexr(:,k)=inde(polr(:,k))
    restr(:,k)=rest(polr(:,k))
END DO
! fill in M for the retired people
! and store the assets left as bequests by people that die for sure next
! period (we do not need them for M, but for the bequest distribution we
! will compute later on
! Note that the
! agent's age always increases by one period. The rows indicate the
! agent's state
! during the current period, the columns indicate her possible state
! next period). Therefore we shift by one group of columns for jj with
! respect to the rows for the agent's aging. Also, given the agent's
! age and assets,  the agent (if she survives), will move into next
! period one period older and with optimal assets
! that depend on the agent policy function. As described before,
! when the agent's policy fn does not coincide with the asset grid
! we have for the value
! fn, we attribute weight to the immediately smaller and bigger assets
! amount in a, proport to their difference to the actual agent's assets
! choice.
counter=1
DO jj=1,TR-1 !age
    DO jjj=1,da !assets
        !attributing mass to smaller and larger asset level tomorrow
        rowM(counter)=indret+(jj-1)*da+jjj
        colM(counter)=indret+jj*da+indexr(jjj,jj)
        valM(counter)=sur(jj+TL)*(1-restr(jjj,jj))
        counter=counter+1
        rowM(counter)=indret+(jj-1)*da+jjj
        colM(counter)=indret+jj*da+indexr(jjj,jj)+1
        valM(counter)=sur(jj+TL)*restr(jjj,jj)
        counter=counter+1
    END DO ! for jjj=1:da,
END DO ! for jj=1:TR-1,

! now consider "oldies" that die for sure to record their assets
counter2=1
DO jjj=1,da
    rowMT(counter2)=jjj
    colMT(counter2)=indexr(jjj,T-TL)
    valMT(counter2)=1-restr(jjj,T-TL)
    counter2=counter2+1
    rowMT(counter2)=jjj
    colMT(counter2)=indexr(jjj,T-TL)+1
    valMT(counter2)=restr(jjj,T-TL)
    counter2=counter2+1
END DO ! DO jjj=1,da,

! loop iterating on the govt budget constraint
epsilon2=1
DO WHILE (epsilon2>1e-04)
	counterbeq=1
    !loop iterating on the bequest distribution
    epsilon=1 ! initialize conv par to compute bequests distribution
	counter2=1
    DO WHILE (epsilon>1e-04)
        CALL TIMDY (IHOUR, IMIN, ISEC)
        WRITE (*,*) IHOUR,IMIN,ISEC
        ! WORKERS, AGE TL: next period the guy is retired and her
        ! value fn has only 3 states
        DO jjj=1,dy   ! jjj, index for y_t
            netearn=(1-taul)*eff(TL)*y(jjj)
            DO j2=1,da ! j2,index for today's assets
                ! compute today's wealth
                weal=(1+r*(1-taua))*a(j2)+netearn
                ! compute utility matrix for each value of consumption
                ! actual consumption values (c is PROPORTION of wealth)
                c=cc(2:dc+1)*weal     !! cc is a vector of dimension dc + 1, where dc = 120. It was constructed as a linear grid between minc and maxc with dc+1 points. The interpretation is the proprtion of this period's income consumed
                u=c**(1-sig)/(1-sig)
                ! WHEN THE GUY IS NOT INHERITING (INHERITED ALREADY)
                ! compute tomorrow's assets (=savings)
                ! for each possible value of cons today
                st=weal-c
                ! interpolate value fn
                CALL interplin(da,a,vr(:,1),dar,st,dpol,xinv,vst)
                ! compute the joy of bequeathing
                ! warm glow is on bequests net of estate taxes
                phit=((1+(st-taub*MAX(st-exb,0*st))/phi2)**(1-sig)-1)*phi1
                !compute value fn (NOT INHERITING)
                valr=u+sur(TL)*bet*vst+(1-sur(TL))*phit
                imax=MAXLOC(valr)
                polt(j2,jjj,dy+1,TL)=imax(1)
                vt(j2,jjj,dy+1,TL)=valr(imax(1))
                sopt(j2,jjj,dy+1,TL)=(1-(REAL(imax(1))/REAL(dc)))*weal
                !copt(j2,jjj,dy+1,TL)=(REAL(imax(1))/REAL(dc))*weal             
                ! WHEN THE GUY IS INHERITING (DID NOT IN THE PAST, AGE TL)
                ! compute tomorrow's assets for each possible
                ! level of bequests received tomorrow (same grid as assets)
                ! and each possible level of consumption today
                DO jj=1,dy ! jj, index for yb_t, jj=1:dy, not inherited
                    evat=0 ! initialize evat
					DO j1=1,da ! each possible bequest level, net of tax
                      at=st+anet(j1)/csur(TL+4)
                      ! compute value fn WHEN INHERITING
                      ! compute the value fn for all the possible levels
                      ! of assets tomorrow (at)
                      ! interpolate value fn
                      CALL interplin(da,a,vr(:,1),dar,at,dpol,xinv,vat)
                      evat=evat+beqdist(j1,jj,TL-3)*vat ! expected value
                   END DO ! DO j1=1,da, j4 is for each beq level
                   valr=u+sur(TL)*bet*evat+(1-sur(TL))*phit
                   imax=MAXLOC(valr)
                   polt(j2,jjj,jj,TL)=imax(1)
                   vt(j2,jjj,jj,TL)=valr(imax(1))
                   sopt(j2,jjj,jj,TL)=(1-(REAL(imax(1))/REAL(dc)))*weal
                   !copt(j2,jjj,jj,TL)=(REAL(imax(1))/REAL(dc))*weal
                END DO ! DO jj=1,dy, jj is for hy_t
            END DO ! DO j2=1,da, j2 is for today's assets
        END DO ! DO jjj=1,dy, jjj is for y_t
     
        ! workers AGE 55 TO AGE 35: CAN INHERIT (NEXT
        ! PERIOD); BUT CANNOT DIE (SO NO JOY OF GIVING NOW)
        DO j0=TL-1,4,-1 ! j0, AGE index IT GOES UNTIL 35 NOW
            write (*,*) j0
            DO jjj=1,dy ! jjj, index for yt
                ! GUYS THAT HAVE INHERITED ALREADY
                ! tomorrow's value fn, depends on tomorrow's
                ! productivity shock, take expectation wrt prod
                vti= MATMUL(vt(:,:,dy+1,j0+1),Qy(jjj,:))				
                DO jj=1,dy ! parent's productivity
					vtai(:,jj)= MATMUL(vt(:,:,jj,j0+1),Qy(jjj,:))
				END DO ! DO jj
				netearn=(1-taul)*eff(j0)*y(jjj)
                DO j2=1,da ! j2,index for today's assets
                    ! compute today's wealth
                    weal=(1+r*(1-taua))*a(j2)+netearn
                    ! utility for each value of consumption
                    c=cc(2:dc+1)*weal
                    u=c**(1-sig)/(1-sig)
                    ! compute tomorrow's assets (=savings)
                    ! for each possible value of cons today
                    st=weal-c
                    ! interpolate value fn
                    CALL interplin(da,a,vti,dar,st,dpol,xinv,vst)
                    !compute value fn (jj=dy+1, inherited in the past)
                    valr=u+bet*vst
                    imax=MAXLOC(valr)
                    polt(j2,jjj,dy+1,j0)=imax(1)
                    vt(j2,jjj,dy+1,j0)=valr(imax(1))
                    sopt(j2,jjj,dy+1,j0)=(1-(REAL(imax(1))/REAL(dc)))*weal
                    !copt(j2,jjj,dy+1,j0)=(REAL(imax(1))/REAL(dc))*weal
                    ! GUYS THAT HAVE NOT YET INHERITED
                    ! compute tomorrow's assets for each possible
                    ! level of bequests received tomorrow
                    ! and each possible level of consumption today
	                DO jj=1,dy ! jj, index for yb_t, jj=1:dy, not inherited
						evat=0 ! initialize evat
						DO j1=1,da ! each possible bequest level, net of tax
							at=st+anet(j1)/csur(j0+4)
							! compute value fn WHEN INHERITING
		                    ! compute the value fn for all the possible levels
		                    ! of assets tomorrow (at)
		                    ! interpolate value fn
		                    CALL interplin(da,a,vti,dar,at,dpol,xinv,vat)
		                    evat=evat+beqdist(j1,jj,j0-3)*vat ! expected value
	                    END DO ! DO j1=1,da, j4 is for each beq level
						CALL interplin(da,a,vtai(:,jj),dar,st,dpol,xinv,vst)                    							
						valr=u+(1-sur(j0+dage))*bet*evat+bet*sur(j0+dage)*vst
	                    imax=MAXLOC(valr)
	                    polt(j2,jjj,jj,j0)=imax(1)
	                    vt(j2,jjj,jj,j0)=valr(imax(1))
	                    sopt(j2,jjj,jj,j0)=(1-(REAL(imax(1))/REAL(dc)))*weal
	                    !copt(j2,jjj,jj,j0)=(REAL(imax(1))/REAL(dc))*weal
	                END DO ! DO jj=1,dy, jj is for hy_t
				END DO ! DO j2=1,da
            END DO ! DO jjj=1,dy, jjj is for h_t
        END DO ! DO j0=TL-1,4,-1 j0 is for age


        ! WORKERS 20 to 30: DO NOT DIE AND DO NOT EXPECT TO INHERIT NOW
        DO j0=3,1,-1 ! j0 age index, start from 30 yrs old
          write (*,*) j0
          DO jjj=1,dy ! jjj index for yt
                netearn=(1-taul)*eff(j0)*y(jjj)
                ! extrapolate value function
                DO j2=1,da ! j2,index for today's assets
                    weal=(1+r*(1-taua))*a(j2)+netearn
                    ! utility for each value of consumption
                    c=cc(2:dc+1)*weal
                    u=efc(j0)*(c/efc(j0))**(1-sig)/(1-sig)
                    ! compute tomorrow's assets (=savings)
                    ! for each possible value of cons today
                    st=weal-c
                    DO jj=1,dy ! jj index for hy_t, not inherited yet
                        ! compute expectation
                        vti=MATMUL(vt(:,:,jj,j0+1),Qy(jjj,:))
                        ! interpolate value fn (do not inherit)
                        CALL interplin(da,a,vti,dar,st,dpol,xinv,vst)
                        ! YOUNG VALUE FN
                        valr=u+bet*vst
                        imax=MAXLOC(valr)
                        polt(j2,jjj,jj,j0)=imax(1)
                        vt(j2,jjj,jj,j0)=valr(imax(1))
                        sopt(j2,jjj,jj,j0)=(1-(REAL(imax(1))/REAL(dc)))*weal
                        !copt(j2,jjj,jj,j0)=(REAL(imax(1))/REAL(dc))*weal
                    END DO !jj=1,dy index for hy_t, not inherited yet
                END DO ! j2=1,da today's assets
            END DO ! jjj=1,dy index for yt
        END DO ! j0=3,1,-1 age index, start from 30 yrs old

    	WRITE (*,*) "val fun computed"
        CALL TIMDY (IHOUR, IMIN, ISEC)
        WRITE (*,*) IHOUR,IMIN,ISEC
	
		!!!! compute the invariant distribution
        ! since we have agents that die and are born at each period, the rows
        ! of the transition matrix M do not sum up to one (people that die
        ! disappear from the system and at each period a new cohort is born):
        ! the people's distribution over the state variables evolves according
        ! to: n_{t+1}' = n_t' M + q_{t+1}'
        ! where q_{t+1} is the distribution of the people that are born at t+1

        ! how is M structured?
        ! the first variable is the agent's age.
        ! For the first TL ages we distinguish by parent's productivity at 40,
        ! current productivity and assets.
        ! For ages from TL+1 to T the agent is retired, therefore for each age
        ! we only distinguish by assets (parent's productivity is not relevant because
        ! they do not inherit any more and current productivity is not relevant
        ! because they do not work any more but just receive pensions).
        !initialize
        ! NOTE THAT WE FILLED IN THE BLOCK FOR RETIRED *OUTSIDE* THE LOOPS
        ! FOR TAUL AND BEQDIST, so keep thos elements and clean up the others
        ! from the previous iteration
        rowM(4*da*2+1:)=0
        colM(4*da*2+1:)=0
        valM(4*da*2+1:)=0
        ! reinitialize counter, keeping the retired guys
        counter=4*da*2+1

        ! working people who have already inherited
        ! AT THE EARLIEST ONE INHERITS AT 40, THEREFORE THESE
        ! GUYS MUST BE AT LEAST 40, SINCE THEY ALREADY INHERITED
        ! like for the retired people, the assets grid on which the working
        ! agents are allowed to choose and the one on which we keep track of
        ! their value function are different. For working people it is easier
        ! to solve the value function substituting for the assets instead of
        ! consumption (as I did for the retired people). Therefore there is
        ! a grid of dc points on consumption from which, in the previous value
        ! fn computation, given the optimal value for consumption, I
        ! compute the value of next's period optimal savings(sopt).
        ! The grid for consumption is different for different ages and
        ! productivity, therefore the same is true for next period's optimal
        ! savings.
        ! find the index corresponding to the optimal savings using the same
        ! mapping described above for the retired agents.
        ! NOTE AGE: FROM 5 (40) TO TL
        indext=FLOOR((da-1)*SQRT(sopt(:,:,dy+1,5:TL))/SQRT(maxa)+1)
        ! set upper bound to da-1
        indext=-(indext-da+1)*(indext<da)+da-1

        ! now find the distance form the smallest closer element to the
        ! element in sopt that we are considering. Note the we put on it
        ! an upper bound of 1
        ! in the case we would go over the maximum allowed grid for a. This is
        ! because we will approximate the transition matrix and the corresponding
        ! invariant distribution using the grid a and giving weight to the two
        ! level of assets between which each optimal choice is included. These
        ! weight will be the distances form those two point (the total distance
        ! is obviously normalized to one since we are talking about weights).
        ! NOTE AGE: FROM 5 (40) TO TL

        DO jj=1,dy
            DO jjj=5,TL ! note different size of restt,sopt,indext (-4)
                restt(:,jj,jjj-4)=(sopt(:,jj,dy+1,jjj)-a(indext(:,jj,jjj-4)))/ &
                &   (a(indext(:,jj,jjj-4)+1)-a(indext(:,jj,jjj-4)))
            END DO !DO jjj=5,TL
        END DO !DO jj=1,dy
        ! take min(restt,1)
        restt=(1-restt)*(restt<1)+1

        ! fill in the blocks for the working people that have already
        ! inherited (yb_t=0), but are not going to retire next period
        ! (those require special treatment).
        DO jj=5,TL-1 ! NOTE AGE: FROM 5 (40) TO TL-1
            DO jjj=1,dy ! today's income shock
                DO j4=1,da ! today's assets
                    DO j5=1,dy ! tomorrow's income shock
                        ! attribute the relevant mass
                        ! to the closest smaller value for tomorrow's assets
                        rowM(counter)=indnoh+(jj-5)*dydy1da+ &  ! Ashish's comment:  parent's productivity does not matter, still it appears in the state space so dy^2*da is added to add all possible combinations of parent's productivity with agents productivity and assets
                        &   dy**2*da+(jjj-1)*da+j4
                        colM(counter)=indnoh+(jj-4)*dydy1da+ &
                        &   dy**2*da+(j5-1)*da+indext(j4,jjj,jj-4)
                        valM(counter)=Qy(jjj,j5)*(1-restt(j4,jjj,jj-4))
                        counter=counter+1
                        !attribute the relevant mass
                        ! to the closest bigger value for tomorrow assets
                        rowM(counter)=indnoh+(jj-5)*dydy1da+ &
                        &   dy**2*da+(jjj-1)*da+j4
                        colM(counter)=indnoh+(jj-4)*dydy1da+ &
                        &   dy**2*da+(j5-1)*da+ &
                        & indext(j4,jjj,jj-4)+1
                        valM(counter)=Qy(jjj,j5)*restt(j4,jjj,jj-4)
                        counter=counter+1
                    END DO ! DO j5=1,dy
                END DO ! DO j4=1,da
            END DO ! DO jjj=1,da
        END DO !DO jj=1,TR-1

        ! fill in the blocks for the working people that have already
        ! inherited (yb_t=0), and are are going to retire next period
        ! (if they will still be alive)
        DO jjj=1,dy
            DO j4=1,da
                ! attribute the relevant mass
                !to the closest smaller value for tomorrow's assets
                rowM(counter)=indnoh+(TL-5)*dydy1da+dy**2*da+(jjj-1)*da+j4
                colM(counter)=indret+indext(j4,jjj,TL-4)
                valM(counter)=sur(TL)*(1-restt(j4,jjj,TL-4))
                counter=counter+1
                rowM(counter)=indnoh+(TL-5)*dydy1da+dy**2*da+(jjj-1)*da+j4
                colM(counter)=indret+indext(j4,jjj,TL-4)+1
                valM(counter)=sur(TL)*restt(j4,jjj,TL-4)
                counter=counter+1
            END DO ! DO j4=1,da
        END DO ! DO jjj=1,da

        ! working people that have not yet inherited
        ! consider two cases: those who do not inherit next period, and thus
        ! remain in the same category and those who do inherit and thus move
        ! to yb_t=0

        ! people who keep not inheriting for another period
        ! take care of their value function first:
        ! with the same logic used above for people who have already inherited,
        ! find the smaller, closest element in the grid for a for each element
        ! in sopt.
        indext2=FLOOR((da-1)*SQRT(sopt(:,:,1:dy,:))/SQRT(maxa)+1)
        indext2=-(indext2-da+1)*(indext2<da)+da-1
        DO jj=1,dy
            DO jjj=1,dy
                DO j4=1,TL
                    restt2(:,jj,jjj,j4)=(sopt(:,jj,jjj,j4)- &
                    &   a(indext2(:,jj,jjj,j4)))/ &
                    &   (a(indext2(:,jj,jjj,j4)+1)-a(indext2(:,jj,jjj,j4)))
                END DO ! DO j4
            END DO ! DO jjj
        END DO ! DO jj
        restt2=(1-restt2)*(restt2<1)+1
        ! age 20 to 35 included
        !fill in the corresponding blocks in the transition matrix M
        DO jj=1,4 ! until 4 because if they have not inherited yet,
            DO j0=1,dy ! for productivity at birth until inheritance
                DO jjj=1,dy
                    DO j4=1,da
                        DO j5=1,dy ! tomorrow's y shock
                            rowM(counter)=(jj-1)*dydyda+(j0-1)*da*dy+ &
                            &   (jjj-1)*da+j4
                            colM(counter)=jj*dydyda+(j0-1)*da*dy+ &
                            &   (j5-1)*da+indext2(j4,jjj,j0,jj)
                            valM(counter)=Qy(jjj,j5)*sur(jj+dage)* &
                            &   (1-restt2(j4,jjj,j0,jj))
                            counter=counter+1
                            rowM(counter)=(jj-1)*dydyda+(j0-1)*da*dy+ &
                            &   (jjj-1)*da+j4
                            colM(counter)=jj*dydyda+(j0-1)*da*dy+(j5-1)*da+ &
                            &   indext2(j4,jjj,j0,jj)+1
                            valM(counter)=Qy(jjj,j5)*sur(jj+dage)* &
                            &   restt2(j4,jjj,j0,jj)
                            counter=counter+1
                        END DO ! DO j5
                    END DO ! DO j4
                END DO ! DO jjj
            END DO ! DO j0
        END DO ! jj

        ! age 40 to 55, included
        !   fill in the corresponding blocks in the transition matrix M
        DO jj=5,TL-1
            DO j0=1,dy ! for productivity at birth until inheritance
                DO jjj=1,dy
                    DO j4=1,da
                        DO j5=1,dy ! tomorrow's y shock
                            rowM(counter)=indnoh+(jj-5)*dydy1da+ &
                            &   (j0-1)*da*dy+(jjj-1)*da+j4
                            colM(counter)=indnoh+(jj-4)*dydy1da+ &
                            &   (j0-1)*da*dy+(j5-1)*da+indext2(j4,jjj,j0,jj)
                            valM(counter)=Qy(jjj,j5)*sur(jj+dage)* &
                            &   (1-restt2(j4,jjj,j0,jj))
                            counter=counter+1
                            rowM(counter)=indnoh+(jj-5)*dydy1da+ &
                            &   (j0-1)*da*dy+(jjj-1)*da+j4
                            colM(counter)=indnoh+(jj-4)*dydy1da+ &
                            &   (j0-1)*da*dy+(j5-1)*da+ &
                            &   indext2(j4,jjj,j0,jj)+1
                            valM(counter)=Qy(jjj,j5)*sur(jj+dage)* &
                            &   restt2(j4,jjj,j0,jj)
                            counter=counter+1
                        END DO ! DO j5
                    END DO ! DO j4
                END DO ! DO jjj
            END DO ! DO j0
        END DO ! jj

        ! working people that inherit
        ! take AGE 4=35 set jj=4
        DO j0=1,dy ! for productivity at birth until inheritance
            DO jjj=1,dy ! today's income
                DO j4=1,da ! today's assets
                    !take care of the grid
                    indext3=FLOOR((da-1)*SQRT(sopt(j4,jjj,j0,4)+ &
                    &   (a-taub*MAX(a-exb,a*0)))/SQRT(maxa)+1)
                    indext3=-(indext3-da+1)*(indext3<da)+da-1
                    restt3=(sopt(j4,jjj,j0,4)+(a-taub*MAX(a-exb,a*0)) &
                    &   -a(indext3))/(a(indext3+1)-a(indext3))
                    restt3=(1-restt3)*(restt3<1)+1 !take min(restt3,1)
                    DO j5=1,dy ! tomorrow's income
                        DO j6=1,da ! inherited assets
                            valM(counter+indext3(j6)-1)= &
                            &   valM(counter+indext3(j6)-1)+ &
                            &   sur(4)*(1-sur(4+dage))*Qy(jjj,j5)* &
                            &   (1-restt3(j6))*beqdist(j6,j0,1)
                            valM(counter+indext3(j6))= &
                            &   valM(counter+indext3(j6))+ &
                            &   sur(4)*(1-sur(4+dage))*Qy(jjj,j5)* &
                            &   restt3(j6)*beqdist(j6,j0,1)
                        END DO ! j6
                        rowM(counter:counter+da-1)=3*dydyda+(j0-1)*da*dy+ &
                        &   (jjj-1)*da+j4
                        colM(counter:counter+da-1)=indnoh+dydyda+(j5-1)*da+ &
                        &   (/ (j7,j7=1,da) /)
                        counter=counter+da
                    END DO ! j5
                END DO ! j4n
            END DO ! jjj
        END DO ! j0

        ! working people that inherit but are not retiring next period
        ! AGE 5=40 YEARS AND OLDER
        DO jj=5,TL-1  ! note TL-1. 4 and TL are special
            tempint1=indnoh+(jj-5)*dydy1da ! index used in loops
            tempint2=indnoh+(jj-4)*dydy1da ! index used in loops
            DO j0=1,dy ! for productivity at birth until inheritance
                DO jjj=1,dy ! today's income
                    DO j4=1,da ! today's assets
                        !take care of the grid
                        indext3=FLOOR((da-1)*SQRT(sopt(j4,jjj,j0,jj)+ &
                        &   (a-taub*MAX(a-exb,a*0)))/SQRT(maxa)+1)
                        indext3=-(indext3-da+1)*(indext3<da)+da-1
                        restt3=(sopt(j4,jjj,j0,jj)+(a-taub*MAX(a-exb,a*0))-a(indext3))/ &
                        &   (a(indext3+1)-a(indext3))
                        restt3=(1-restt3)*(restt3<1)+1 !take min(restt3,1)
                        DO j5=1,dy ! tomorrow's income
                            DO j6=1,da ! inherited assets
                                valM(counter+indext3(j6)-1)= &
                                &   valM(counter+indext3(j6)-1)+ &
                                &   sur(jj)*(1-sur(jj+dage))*Qy(jjj,j5)* &
                                &   (1-restt3(j6))*beqdist(j6,j0,jj-3)
                                valM(counter+indext3(j6))= &
                                &   valM(counter+indext3(j6))+ &
                                &   sur(jj)*(1-sur(jj+dage))*Qy(jjj,j5)* &
                                &   restt3(j6)*beqdist(j6,j0,jj-3)
                            END DO ! j6
                            rowM(counter:counter+da-1)=tempint1+ &
                            &   (j0-1)*da*dy+(jjj-1)*da+j4
                            colM(counter:counter+da-1)=tempint2+ &
                            &   dydyda+(j5-1)*da+(/ (j7,j7=1,da) /)
                            counter=counter+da
                        END DO ! j5
                    END DO ! j4
                END DO ! jjj
            END DO ! j0
        END DO !jjx
        ! working people that inherit and retire next period
        tempint1=indnoh+(TL-5)*dydy1da
        DO j0=1,dy ! for productivity at birth until inheritance
            DO jjj=1,dy
                DO j4=1,da
                    !take care of the grid
                    indext3=FLOOR((da-1)*SQRT(sopt(j4,jjj,j0,TL)+ &
                    &   (a-taub*MAX(a-exb,a*0)))/SQRT(maxa)+1)
                    indext3=-(indext3-da+1)*(indext3<da)+da-1
                    restt3=(sopt(j4,jjj,j0,TL)+(a-taub*MAX(a-exb,a*0)) &
                    &   -a(indext3))/(a(indext3+1)-a(indext3))
                    restt3=(1-restt3)*(restt3<1)+1
                    DO j6=1,da ! inherited assets
                        valM(counter+indext3(j6)-1)=valM(counter+indext3(j6)-1)+ &
                        &   sur(TL)*(1-restt3(j6))*beqdist(j6,j0,TL-3)
                        valM(counter+indext3(j6))=valM(counter+indext3(j6))+ &
                        &   sur(TL)*restt3(j6)*beqdist(j6,j0,TL-3)
                    END DO ! j6
                    rowM(counter:counter+da-1)=tempint1+(j0-1)*da*dy+ &
                    &   (jjj-1)*da+j4
                    colM(counter:counter+da-1)=indret+(/ (j7,j7=1,da) /)
                    counter=counter+da
                END DO ! j4
            END DO ! jjj
        END DO !j0
		
        q=qinit
		epsi=1
        DO WHILE (epsi>1e-08)
            invm1=invm
            invm=0.0
            ! this do loop is the product M'*invm, M sparse
            ! to transpose M, we simply use colM instead of rowM
            ! invm= new inv distr
            ! invm1= old inv distr
            DO i=1,nonzeroM
                invm(colM(i))=invm(colM(i))+(invm1(rowM(i))*valM(i))/gpop
            END DO
            invm(1:dydyda)=invm(1:dydyda)+q ! add newborns
            epsi=MAXVAL(ABS(invm-invm1))            
        END DO
        invm=invm/SUM(invm)
        q=q/SUM(invm)
       
        ! compute bequests left
        beqdist=0
        ! the 20 years old child observes the parent's productivity
        ! when the parent is 40. therefore the bequest
        ! distribution will be conditional on the parent's productivity
        ! at 40 years of age: beqdis(a|yp40,t).
        ! to contruct this object we need to go through several steps:
        ! 1) find the asset distribution of the parent at 40 years of age:
        ! m(a,y40,yp40)   ! Ashish's comment: yp40 here refers to parents productivity of parents, i.e. grandparents of the child, coz parents asset distribution when they are 40 also depend on their parents productivity when they are 40
        ! to do so, pick the right age in invariant distribution
        temprealv2=invm(4*dy*dy*da+1:4*dy*dy*da+dy*(dy+1)*da)
        ! 2) we want this distribution to be conditional on y40, the parents'
        ! income at 40, because this is what the child USES when making
        ! inference m(a, yp40|y40). to compute, use the definition:
        !
        ! pr(yp=yp_j,a=a_k|y=y_i)=
        ! pr(yp=yp_j,a=a_k,y=y_i)/(sum_(j,k) pr(yp=yp_j,a=a_k,y=y_i)
        !
        ! stack these elements into vector temprealv ordered as
        ! y_p, y_40, a (convenient later, same order as in invm)
        temprealv=0
		temprealv4=0 ! marginal of y40
        DO i=1,dy ! y40
            DO j=1,dy+1 ! productivity of the grandparent at 40 (might be dead)
                DO k=1,da ! parent's asset at 40
                    temprealv4(i)=temprealv4(i)+temprealv2((j-1)*dy*da+(i-1)*da+k)
                END DO ! k
            END DO ! DO j
        END DO ! DO i
        DO j=1,dy+1 ! productivity of the grandparent at 40 (might be dead)
            DO k=1,da ! parent's asset at 40
                DO i=1,dy ! productivity of the parent at 40
                     temprealv((j-1)*dy*da+(i-1)*da+k)= &
                     &  temprealv2((j-1)*dy*da+(i-1)*da+k)/temprealv4(i)
                END DO ! DO i
            END DO ! DO k
        END DO ! DO j
			 
        ! 2) now we need m(a,y40,yp40|y40) (*):
        ! the matrix of characteristics of the parent at 40
        ! conditional on what the child observes, i.e. the parent's y.
        ! this is because we need also m(a',y45,yp40|y40),m(a'',y50,yp40|y40)
        ! ecc... for all ages and to obtain them we moltiply (*) by
        ! the appropriate part of the transition matrix.
        ! so now: construct m(a,y40,yp40|y40) from m(a,yp40|y40)
        ! Pr(yb,y,a|parent's y)=Pr(parent's yb,a|parent's y)*indicator fn of y
        ! stbeq is of size(da*dy*(dy+1),dy): stack a vector of probabilities,
        ! conditional on y_40 in each column
        stbeq=0
        DO i=1,dy ! productivity of the parent at 40
            DO jj=1,dy+1 ! productivity of the grandparent (might be dead)
                stbeq((jj-1)*dy*da+(i-1)*da+(/ (j7,j7=1,da) /),i)= &
                &   temprealv((jj-1)*dy*da+(i-1)*da+(/ (j7,j7=1,da) /))
            END DO ! jj
        END DO ! DO i

		! now we have the parent's characteristics at 40, want to have them
        ! as she ages. start aging the parent to 60, when her child does
        ! not inherit, in matlab this was:
        ! for jj=dage:TL-1
        ! stbeq=full(stbeq*M((jj-1)*dy*(dy+1)*da+(1:dy*(dy+1)*da),...
        !      jj*dy*(dy+1)*da+(1:dy*(dy+1)*da)))/sur(dage)
        ! end
        ! here it is more complicated because of the sparse storage
        ! remember how invm, M were constructed and find the right indexes
        ! dage = 5 (parent-child age difference)
        DO jj=dage,TL-1  ! for loop from dage = 5 to TL-1 = 8
            stbeq2=0 ! initialization
            ! first, working that have already inherited

            ! i here is a proxy for counter. Understanding that will help a lot.
            ! My comment: Transition matrix started for retirees, so 4*da*2 is added to index everywhere to get across the retirees part of the transition matrix, i.e. by coming across the 4*da*2 elemeents in the rows of transition matrix, the states of working people start
            DO i=4*da*2+(jj-5)*dydyda*2+1,4*da*2+(jj-4)*dydyda*2 ! this i is in reference to "counter" in the previous code. Check the line where "counter is defined"
                tempint1=colM(i)-indnoh-(jj-4)*dydy1da
                tempint2=rowM(i)-indnoh-(jj-5)*dydy1da
                stbeq2(tempint1,:)=stbeq2(tempint1,:)+ &
                &   stbeq(tempint2,:)*valM(i)
            END DO ! i
            ! second, people that have not yet inherited and keep not inheriting
            DO i=4*da*2+4*dydyda*2+dy*da*2+4*dydyda*dy*2+ & 
                &   (jj-5)*dydyda*dy*2+1, &  !For case-1, check how many 1's were added to counter relative to 4*da*2. That should be the first element of i for this case
                &    4*da*2+4*dydyda*2+dy*da*2+4*dydyda*dy*2+ &
                &   (jj-4)*dydyda*dy*2
                tempint1=colM(i)-indnoh-(jj-4)*dydy1da
                tempint2=rowM(i)-indnoh-(jj-5)*dydy1da
                stbeq2(tempint1,:)=stbeq2(tempint1,:)+ &
                &   stbeq(tempint2,:)*valM(i)
            END DO ! i
            ! lastly, people that have not inherited but inherit next period
            DO i=4*da*2+4*dydyda*2+dy*da*2+8*dydyda*dy*2+dydyda*dy*da+ &
                &   (jj-5)*dydyda*dy*da+1, &
                  4*da*2+4*dydyda*2+dy*da*2+8*dydyda*dy*2+dydyda*dy*da+ &
                &   (jj-4)*dydyda*dy*da
                tempint1=colM(i)-indnoh-(jj-4)*dydy1da
                tempint2=rowM(i)-indnoh-(jj-5)*dydy1da
                stbeq2(tempint1,:)=stbeq2(tempint1,:)+ &
                &   stbeq(tempint2,:)*valM(i)
            END DO ! i
            stbeq=stbeq2/sur(jj)
        END DO !jj
       	
		! age to 60 (the parent that is going to retire next period)
        stbeqret=0
        stbeqret2=0
        ! first people that have already inherited
        DO i=4*da*2+4*dydyda*2+1,4*da*2+4*dydyda*2+dy*da*2
            tempint1=colM(i)-indret
            tempint2=rowM(i)-indnoh-(TL-5)*dydy1da
            stbeqret2(tempint1,:)=stbeqret2(tempint1,:)+ &
            &   stbeq(tempint2,:)*valM(i)
        END DO ! i
        ! now people that have NOT YET inherited
        DO i=4*da*2+4*dydyda*2+dy*da*2+8*dydyda*dy*2+5*dydyda*dy*da+1, &
            &   nonzeroM
            tempint1=colM(i)-indret
            tempint2=rowM(i)-indnoh-(TL-5)*dydy1da
            stbeqret2(tempint1,:)=stbeqret2(tempint1,:)+ &
            &   stbeq(tempint2,:)*valM(i)
        END DO ! i
        stbeqret=stbeqret2/sur(TL)
      	beqdist(:,:,TL-dage-3)=stbeqret       
        ! parent that is retired but does not die for sure next period
        DO jj=TL+1,T-1
            stbeqret2=0
            DO i=(jj-TL-1)*da*2+1,(jj-TL)*da*2
                tempint1=colM(i)-indret-(jj-TL)*da
                tempint2=rowM(i)-indret-(jj-TL-1)*da
                stbeqret2(tempint1,:)=stbeqret2(tempint1,:)+ &
                &   stbeqret(tempint2,:)*valM(i)
            END DO ! i
            stbeqret=stbeqret2/sur(jj)
            beqdist(:,:,jj-dage-3)=stbeqret
        END DO ! jj
        ! parent that is T years old and is going to die for sure next period
        stbeqret2=0
        DO i=1,nonzeroMT
            stbeqret2(colMT(i),:)=stbeqret2(colMT(i),:)+ &
            &   stbeqret(rowMT(i),:)*valMT(i)
        END DO
        beqdist(:,:,T-dage-3)=stbeqret2
		
		!kkkkkkkkkkk ***** modified May 6th
		! smooth beqdist with kernels
		kercoeff2=kercoeff
		beqdist2=beqdist
		kercoeff=0.0
		DO i=1,TL-3
			kercoeff(:,:,i)=MATMUL(kernelreg,beqdist(:,:,i))
			beqdist(:,:,i)=MATMUL(kernelfun,kercoeff(:,:,i))
		END DO !DO i	
		!set to zero negative values
		WHERE (beqdist<0) beqdist=0.0 
		! kkkkkk ****** renormalize smoothed distribution
		DO i=1,TL-3
			DO jj=1,dy
				beqdist(:,jj,i)=beqdist(:,jj,i)/SUM(beqdist(:,jj,i))
			END DO !DO jj
		END DO !DO i	
		epsilon=MAXVAL(ABS(kercoeff-kercoeff2))
        WRITE (*,*) "counterbeq=",counterbeq
        WRITE (*,*) "epsilon=", epsilon		
		counterbeq=counterbeq+1

		! save stuff in matrices to make graphs
		OPEN(UNIT=14, FILE= "UStemp19.dat", STATUS="REPLACE", ACTION= "WRITE", &
		&   POSITION= "REWIND", IOSTAT= OpenStatus)
		IF (OpenStatus>0) STOP "***** Cannot open file *****"
		WRITE (14,*) da
		WRITE (14,*) dar
		WRITE (14,*) dy
		WRITE (14,*) T
		WRITE (14,*) TL
		WRITE (14,*) dage
		WRITE (14,*) ninvm
		WRITE (14,*) taua
		WRITE (14,*) taub
		WRITE (14,*) exb
		WRITE (14,*) indret
		WRITE (14,*) indnoh
		DO i=1,dy
			WRITE (14,*) y(i)
		END DO ! DO i
		DO i=1,da
			WRITE (14,*) a(i)
		END DO ! DO i
		DO i=1,T+4
			WRITE (14,*) csur(i)
		END DO ! DO i
		DO i=1,T
			WRITE (14,*) sur(i)
		END DO ! DO i
		DO i=1,ninvm
			WRITE (14,*) invm(i)
		END DO ! DO i
		DO k=1,TL
			DO j=1,dy+1
				DO j0=1,dy
					DO i=1,da
						WRITE(14,*) sopt(i,j0,j,k)
					END DO ! DO i
				END DO ! DO j0
			END DO ! DO j
		END DO ! DO k
		DO k=1,TL-3
			DO j0=1,dy
				DO i=1,da
					WRITE(14,*) beqdist(i,j0,k)
				END DO ! DO i
			END DO ! DO j0
		END DO ! DO k
		DO k=1,TL-3
			DO j0=1,dy
				DO i=1,da
					WRITE(14,*) beqdist2(i,j0,k)
				END DO ! DO i
			END DO ! DO j0
		END DO ! DO k
		DO i=1,dar
			WRITE (14,*) ar(i)
		END DO ! DO i
		DO i=1,TL*dy
			WRITE (14,*) eminv(i)
		END DO
		CLOSE(14)
	END DO ! WHILE (BEQUESTS)

    WRITE (*,*) "beqdist computed"
    CALL TIMDY (IHOUR, IMIN, ISEC)
    WRITE (*,*) IHOUR,IMIN,ISEC
	
    ! compute the government budget constraint
    ! compute distribution of people over each asset level
    dista=0
    DO i=1,4*dy*dy+5*(dy+1)*dy+5 ! sum people up
        dista=dista+invm((i-1)*da+1:i*da)
    END DO !i

    ! compute average capital (multiply each asset level for
    ! the number of people owning it, working or retired)
    avgk=SUM(dista*a)

    ! construct a matrix with assets on the rows and age on the columns
    !        age 1  age 2 ......age TL.....ageT
    !  a1
    !  a2
    !  ...
    !  ada
    agew=0
    DO j=1,4 ! age 20-35
        DO i=1,dy*dy
            agew(:,j)=agew(:,j)+invm((i-1)*da+(j-1)*dydyda+1: &
            &   i*da+(j-1)*dydyda)
        END DO ! DO i
    END DO ! DO j
    DO  j=5,TL !age 40-60
        DO i=1,dy*(dy+1)
            agew(:,j)=agew(:,j)+invm(indnoh+(i-1)*da+(j-5)*dydy1da+1: &
            &   indnoh+i*da+(j-5)*dydy1da)
        END DO ! DO i
    END DO ! DO j
    ! age 65-90
    DO j=TL+1,T
        agew(:,j)=invm(indret+(j-TL-1)*da+1:indret+(j-TL)*da)
    END DO ! j
    beqtax=0
    tempint1=COUNT(a<exb)
    ! compute the asset distribution of the 90 years old
    agewT=0
    DO i=1,nonzeroMT
        agewT(colMT(i))=agewT(colMT(i))+agew(rowMT(i),T)*valMT(i)
    END DO ! i
    DO i=tempint1+1,da
        tempreal=taub*(a(i)-exb)
        DO j=TL+1,T
            beqtax=beqtax+tempreal*agew(i,j)*(1-sur(j-1))/sur(j-1)
        END DO ! j
        ! 90 years old
        beqtax=beqtax+tempreal*agewT(i)
    END DO ! i

    ! government balance
	write (*,*) "taua=",taua,"r=",r,"avgk=",avgk,"taul=",taul,"aveinc=",aveinc
	write (*,*) "pens=",pens,"numret=",numret,"ggovt=",ggovt,"beqtax=",beqtax	
    gbal=taua*r*avgk+taul*aveinc-pens*numret-ggovt*(aveinc+r*avgk)+beqtax   
	taul2=taul
    taul=taul-relax*(gbal/aveinc)
    epsilon2=ABS(taul-taul2)
    WRITE (*,*) "new taul computed",taul

    CALL TIMDY (IHOUR, IMIN, ISEC)
    WRITE (*,*) IHOUR,IMIN,ISEC
END DO ! WHILE TAUL

! compute statistics on inequality for the whole population
! compute vectors for lorenz curve
CALL cumsumv(dista,cumdista) ! cumulative number of agents
CALL cumsumv(a*dista,cumwealth)
cumwealth=cumwealth/SUM(a*dista) ! cumulative wealth
! compute gini for whole population
area=0
DO i=1,da-1
    area=area+(cumwealth(i+1)+cumwealth(i))*(cumdista(i+1)-cumdista(i))
END DO ! i
wgini=1-area

! compute statistics on inequality for the whole population
! compute distribution of people over assets, EXCLUDING 20 YEARS OLD
distac=0
DO i=1,3*dy*dy+5*(dy+1)*dy+5 ! sum people up
    distac=distac+invm(dydyda+(i-1)*da+1:dydyda+i*da)/(1-pfrac(1))
END DO !i
! compute vectors for lorenz curve
CALL cumsumv(distac,cumdistac) ! cumulative number of agents
CALL cumsumv(a*distac,cumwealthc)
cumwealthc=cumwealthc/SUM(a*distac) ! cumulative wealth
! compute gini, population without the youngsters
area=0
DO  i=1,da-1
    area=area+(cumwealthc(i+1)+cumwealthc(i))*(cumdistac(i+1)-cumdistac(i))
END DO ! i
wginic=1-area

! top percentiles for wealth for the population
! topwel=[cumsum(flipud(a*dist/sum(a*dist))) cumsum(flipud(dist))]
temprealda=a*dista/SUM(a*dista)
! flipud
DO i=1,da
    temprealda2(i)=temprealda(da-i+1)
	temprealda4(i)=dista(da-i+1)
END DO !i
! cumsumv
CALL cumsumv(temprealda2,temprealda)
topwel(:,1)=temprealda
CALL cumsumv(temprealda4,temprealda)
topwel(:,2)=temprealda

! top percentiles for wealth for the population excluding 20 years old
!topwelc=[cumsum(flipud(a.*distc/sum(a.*distc))) cumsum(flipud(distc))]
temprealda=a*distac/SUM(a*distac)
! flipud
DO i=1,da
    temprealda2(i)=temprealda(da-i+1)
	temprealda4(i)=distac(da-i+1)
END DO !i
! cumsumv
CALL cumsumv(temprealda2,temprealda)
topwelc(:,1)=temprealda
CALL cumsumv(temprealda4,temprealda)
topwelc(:,2)=temprealda

! correct the number of people that are at zero
! excluding those that are there because of a coarse grid for capital
fpoor=0
! do workers first
DO jj=1,4 ! age 20-35
    DO j=1,dy !parent's productivity
        DO k=1,dy ! agent's productivity
            DO i=1,da !assets
                IF (sopt(i,k,j,jj)>a(1) .AND. sopt(i,k,j,jj)<a(2)) THEN
                    fpoor=fpoor+(1-sopt(i,k,j,jj)/a(2))* &
					&	invm(i+(k-1)*da+(j-1)*da*dy+(jj-1)*dydyda)*sur(jj)
                END IF
            END DO !i
        END DO !k
    END DO ! j
END DO ! jj
DO jj=5,9 ! age 40-60
    DO j=1,dy+1 !parent's productivity
        DO k=1,dy ! agent's productivity
            DO i=1,da !assets
                IF (sopt(i,k,j,jj)>a(1) .AND. sopt(i,k,j,jj)<a(2)) THEN
                    fpoor=fpoor+(1-sopt(i,k,j,jj)/a(2))* &
					&	invm(4*dydyda+i+(k-1)*da+(j-1)*da*dy+ &
					&	(jj-5)*dydy1da)*sur(jj)
                END IF
            END DO !i
        END DO !k
    END DO ! j
END DO ! jj
! now retired
DO jj=1,5 !age 65-85
    DO i=1,da
        IF (polr(i,jj)>a(1) .AND. polr(i,jj)<a(2)) THEN
			fpoor=fpoor+(1-polr(i,jj)/a(2))*invm(4*dydyda+5*dydy1da+(jj-1)*da+i)*sur(jj)
        END IF
    END DO ! i
END DO !jj

! compute transfer wealth to perform Kotlikoff-Summers' computation
! general idea of how we are going to do it:
! we start from the parents that leave the bequest. this
! is much easier: since we are not interested in tracing who
! gets the transfer wealth, but only how much it is in the aggregate,
! we avoid having to keep track of all the state variables of the
! children.
! we are not computing the average bequest received, but the average
! bequest LEFT. furthermore, it is NOT true that the fraction
! of people who die and the fraction of people who are born are
! the same: population grows!
! for the parents that die at each age, we compute the bequest
! that they leave and we weigh it by the appropriate measure
! to get it in per capita terms. now this wealth is in the hands of
! the children: Kotlikoff and Summers assume that each period they are
! alive this transfer wealth increases with the net-of-taxes interest
! rate. Of course, this only works for those who are still alive:
! the children who die leave a bequest to other people and do not
! count anymore. Here we use the stationarity of our problem:
! we follow the wealth left by the parents AS IF it was given
! to the children who gradually age; instead, what we are
! interested is a cross-sectional computation. we use the fact that
! the people who inherited last period are similar to how the
! people that inherit this period will be in the next period.
! OK, let us start.
trwealth=0
! compute average inheritance left in the population
! use the fact that, given the age, the probability
! of dying does not depend on the asset level.
! hence the average conditional bequest is related to the assets
! distribution in the population, weighted with the relevant prob. of dying
! and conditional on probability of receiving a bequest (related to the
! fraction of people dying).
! bequests are left by retired + retired that are T years old
! (use the MT matrix keeps track of the evolution of their assets)
! the assets left as bequest are the asset the
! parent would have if she were still alive next period
! further note: the bequests are expressed in per capita
! terms, counting only the adults; this is OK, because
! the same is true for capital
DO j=TL,T-1 ! age of the parent just before dying
   ! these are people that would be retired in
   ! the following period, if they survived;
   tempreal=(1-sur(j))* &
   &	DOT_PRODUCT(a,invm(indret+(j-TL)*da+1:indret+(j-TL+1)*da))
   ! temp is the amount of bequests these guys leave. Now we
   ! need to keep track of how it evolves when the children have it
   ! the children inherit this bequest at age j+1-dage; add 3 to
   ! get the probability that they are still alive at the various
   ! following ages (divide by the probability that they are still
   ! alive at the age of the bequest, because the ones who get
   ! the bequest are indeed alive),
   ! and divide by gpop to take into account
   ! that the population grows and hence that the fraction of
   ! this people (and their wealth) becomes gradually smaller.
   trwealth=trwealth+tempreal*SUM((csur(j+1-dage+3:T-1+3)/ &
   &   csur(j+1-dage+3))*((1+r*(1-taua))/gpop)** &
   &   (/ (i,i=0,T-1+3-j-1+dage-3) /))
   ! I left a lot of -1-3+1 etc. to make things clearer; one should
   ! just sum them out; also, the indexes are quite complicated,
   ! but it is easy to think of the extremes (say, start with j=T-1, i.e.
   ! the parent that dies between 85 and 90; the child inherits at 65
   ! the second term should correspond to the probability that it
   ! survives to 70 years conditional on having arrived at 65, and so on...).
END DO
! add now  retired that are T years old (use the MT matrix that
! keeps track of the evolution of their assets and the invariant distr invm
temprealda2=invm(indret+(T-TL-1)*da+1:indret+(T-TL)*da)
temprealda=0
DO i=1,nonzeroMT
	temprealda(colMT(i))=temprealda(colMT(i))+temprealda2(rowMT(i))*valMT(i)
END DO ! DO i
tempreal=DOT_PRODUCT(temprealda,a)
trwealth=trwealth+tempreal*SUM((csur(T+1-dage+3:T-1+3)/ &
&	csur(T+1-dage+3))*((1+r*(1-taua))/gpop)** &
&	(/ (i,i=0,T-1+3-T-1+dage-3) /))
! take the ratio of transfer wealth to average capital
kotsum=trwealth/avgk

! note that we have net national product so far, we want to compute GNP,
! given the depreciation rate for capital.
gnp=aveinc+(r+delt)*avgk
k2gnp=(avgk*tlength)/gnp
kshare=((r+delt)*avgk)/gnp

DATA sp/.01, .05, .10, .20, .40, .60, .80/
CALL interplin(da,topwelc(:,2),topwelc(:,1),7,sp,dpol,xinv,squant)

! save stuff in matrices to make graphs
OPEN(UNIT=14, FILE= "USdat19.dat", STATUS="REPLACE", ACTION= "WRITE", &
&   POSITION= "REWIND", IOSTAT= OpenStatus)
IF (OpenStatus>0) STOP "***** Cannot open file *****"
WRITE (14,*) da
WRITE (14,*) dar
WRITE (14,*) dy
WRITE (14,*) T
WRITE (14,*) TL
WRITE (14,*) dage
WRITE (14,*) ninvm
WRITE (14,*) taua
WRITE (14,*) taub
WRITE (14,*) exb
WRITE (14,*) indret
WRITE (14,*) indnoh
DO i=1,dy
	WRITE (14,*) y(i)
END DO ! DO i
DO i=1,da
	WRITE (14,*) a(i)
END DO ! DO i
DO j=1,T
	DO i=1,da
		WRITE (14,*) agew(i,j)
	END DO ! DO i
END DO ! DO j	 
DO i=1,da
	WRITE (14,*) agewT(i)
END DO ! DO i
DO i=1,T+4
	WRITE (14,*) csur(i)
END DO ! DO i
DO i=1,T
	WRITE (14,*) sur(i)
END DO ! DO i
DO i=1,ninvm
	WRITE (14,*) invm(i)
END DO ! DO i
DO k=1,TL
	DO j=1,dy+1
		DO j0=1,dy
			DO i=1,da
				WRITE(14,*) sopt(i,j0,j,k)
			END DO ! DO i
		END DO ! DO j0
	END DO ! DO j
END DO ! DO k
DO k=1,TL-3
	DO j0=1,dy
		DO i=1,da
			WRITE(14,*) beqdist(i,j0,k)
		END DO ! DO i
	END DO ! DO j0
END DO ! DO k
DO i=1,dar
	WRITE (14,*) ar(i)
END DO ! DO i
CLOSE(14)

! save stuff in matrices to make graphs
OPEN(UNIT=16, FILE= "USvf19.dat", STATUS="REPLACE", ACTION= "WRITE", &
&   POSITION= "REWIND", IOSTAT= OpenStatus)
IF (OpenStatus>0) STOP "***** Cannot open file *****"
DO k=1,TR
	DO i=1,da
		WRITE(16,*) vr(i,k)
	END DO ! DO i	
END DO ! DO k
DO k=1,TL
	DO j=1,dy+1
		DO j0=1,dy
			DO i=1,da
				WRITE(16,*) vt(i,j0,j,k)
			END DO ! DO i
		END DO ! DO j0
	END DO ! DO j
END DO ! DO k
CLOSE(16)


!kkkkkkkkkkkkkkk
! save stuff in matrices to make graphs
OPEN(UNIT=18, FILE= "USkerco19.dat", STATUS="REPLACE", ACTION= "WRITE", &
&   POSITION= "REWIND", IOSTAT= OpenStatus)
IF (OpenStatus>0) STOP "***** Cannot open file *****"
DO k=1,TL-3
	DO j=1,dy
		DO j0=1,numkernel
			WRITE(18,*) kercoeff(j0,j,k)			
		END DO ! DO j0
	END DO ! DO j
END DO ! DO k
DO k=1,TL-3
	DO j=1,dy
		DO j0=1,numkernel
			WRITE(18,*) kercoeff2(j0,j,k)			
		END DO ! DO j0
	END DO ! DO j
END DO ! DO k
CLOSE(18)


OPEN(UNIT=12, FILE= "USdker19.dat", STATUS="REPLACE", ACTION= "WRITE", &
&   POSITION= "REWIND", IOSTAT= OpenStatus)
IF (OpenStatus>0) STOP "***** Cannot open file *****"
WRITE (12,*) "------------------------ run USD18 --------------------------- "
WRITE (12,*) "-------------------------------------------------------------- "
WRITE (12,*) " parameter values"
WRITE (12,*) "-------------------------------------------------------------- "
WRITE (12,*) "sig=",sig,"bet=",bet**(1/REAL(tlength)),"dy=",dy
WRITE (12,*) "gam=",gam,"varz=",sigz**2
WRITE (12,*) "hbet=",hbet,"varh=",sigh**2
WRITE (12,*) "gpop=",gpop**(1.0/tlength),"r=",(r+1)**(1/REAL(tlength))-1
WRITE (12,*) "ggovt=",ggovt,"pens=",pens
WRITE (12,*) "taua=",taua,"taub=",taub,"exb=",exb
WRITE (12,*) "numret=",numret,"phi1=",phi1,"phi2=",phi2
WRITE (12,*) "maxa=",maxa,"da=",da,"dc=",dc,"dar=",dar
WRITE (12,*) "y=",y
WRITE (12,*) "cvary=",cvary
WRITE (12,*) "cvaryh=",cvaryh
WRITE (12,*) "sigy=",SQRT(vary),"sigyh=",SQRT(varyh)
WRITE (12,*) "-------------------------------------------------------------- "
WRITE (12,*) "gini indexes for productivity"
WRITE (12,*) "gini=",gini,"agini=",agini
WRITE (12,*) "gini1=",gini1,"gini2=",gini2
WRITE (12,*) "-------------------------------------------------------------- "
WRITE (12,*) "output values"
WRITE (12,*) "taul=",taul,"k2gnp=",k2gnp,"kshare=",kshare
WRITE (12,*) "kotsum=",kotsum,"rnet=",(1-taua)*((r+1)**(1.0/tlength)-1)
WRITE (12,*) "wgini=",wgini,"wginic=",wginic
WRITE (12,*) "poor=",dista(1),"poorc=", distac(1)-fpoor/(1-pfrac(1))
WRITE (12,*) "-------------------------------------------------------------- "
WRITE (12,*) "wealth quantiles"
WRITE (12,*) "quantiles=",sp
WRITE (12,*) "wealth=",squant*100
WRITE (12,*) "-------------------------------------------------------------- "
WRITE (12,*) "topwel"
WRITE (12,*) topwel(1:7,2)*100
WRITE (12,*) "-------------------------------------------------------------- "
CLOSE (UNIT=12)
        
END PROGRAM forlife

SUBROUTINE cumsumv(AA,BB)
! cumulated sum
! SPECIAL CASE WHERE AA IS A VECTOR
IMPLICIT NONE
REAL,INTENT(IN),DIMENSION(:) :: AA
REAL,INTENT(OUT),DIMENSION(:) :: BB
INTEGER i
DO i=1,SIZE(AA)
    BB(i)=SUM(AA(1:i))
END DO
END SUBROUTINE cumsumv

SUBROUTINE  interplin(l,x,y,n,z,dpol,xinv,v)
! interpolate, from values of the fn given as ordinates of
! input data points in an x-y plane and for a given set of
! x values (abscissas), the values of a single-valued fn y=f(x)
! input parameters are
!     l  = # of data points (2 or greater)
!     x  = array of dimension l storing the x values (abscissas)
!          of input data points (in ascending order)
!     y  = array of dimension l storing the y values (ordinates)
!          of input data points
!     n  = # of points at which interpolation of the y value (ordinate)
!          (ordinate) is needed (1 or greater)
!     z  = array of dimension n storing the x values (abscissas) of
!          desired points
!     dpol-1 = degree of polynomial extrapolation
!     xinv = matrix for extrapolation computed in main code to save
!            time
! output parameter
!     v  = array of dimension n where the interpolated y values
!          (ordinates) are to be displayed
IMPLICIT NONE

INTEGER, INTENT(IN) :: l
INTEGER, INTENT(IN) :: n
INTEGER (KIND=1), INTENT(IN) :: dpol
REAL, DIMENSION(l), INTENT(IN) :: x
REAL, DIMENSION(l), INTENT(IN) :: y
REAL, DIMENSION(n), INTENT(IN) :: z
REAL, DIMENSION(dpol,dpol), INTENT(in) :: xinv
REAL, DIMENSION(n), INTENT(OUT) :: v

INTEGER :: i,ind
INTEGER, DIMENSION(1) :: k
REAL :: diff
REAL, DIMENSION(dpol) :: zpoly

DO i=1,n
    k=MAXLOC(-ABS(z(i)-x))
    ind=k(1)
    diff=z(i)-x(ind)
    IF (abs(diff)<1e-06) THEN
        v(i)=y(ind)
    ELSE
        IF ((ind.EQ.l).AND.(diff.GE.0)) THEN
            zpoly=( (/ (z**i,i=0,dpol-1) /) )
            v(i)=DOT_PRODUCT(zpoly,MATMUL(xinv,y(l-dpol+1:l)))
        ELSE
            IF (diff<0) THEN
                v(i)=y(ind-1)+(z(i)-x(ind-1))/(x(ind)-x(ind-1))&
                     & *(y(ind)-y(ind-1))
            ELSE
                v(i)=y(ind)+(z(i)-x(ind))/(x(ind+1)-x(ind))*(y(ind+1)-y(ind))
            END IF
        END IF
    END IF
END DO ! DO i=1,n
END SUBROUTINE interplin

SUBROUTINE linspace(xmin,xmax,npoints,lspace)
IMPLICIT NONE
REAL, INTENT(IN) :: xmin,xmax
INTEGER, INTENT(IN) :: npoints
REAL, DIMENSION(:), INTENT(OUT) :: lspace
INTEGER :: i
lspace=real( (/ (i,i=0,npoints) /) )/real(npoints-1)*(xmax-xmin)+xmin
END SUBROUTINE linspace

