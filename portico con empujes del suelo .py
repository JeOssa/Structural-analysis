import sympy as sy 
import numpy as np 
import matplotlib.pyplot as plt

xE,xA,xB,xD,xC,EI,Q,x,xi,xiA,xiB,xE,xiE,xiD,Lambda,E,I,A,L,EI,AE,C1,C2,C1e,C2e =sy.symbols ('xE,xA,xB,xD,xC,EI,Q,x,xi,xiA,xiB,xE,xiE,xiD,Lambda,E,I,A,L,EI,AE,C1,C2,C1e,C2e')

#  1.0 Propiedades del modelo

E=1.95e7           #Modulo de elasticidad
R=0.9              #Radio de la seccion transversal de los elementos tipo pilas AB
k=42000
#   2.0 Definición de longitud de elementos 
LA=5
LB=5
LC=9
LD=2*(5)**0.5
LE=(29)**0.5

# Dfinicion de cargas en los elementos D y E 
QE=30
QD=60

cos=4/LD
sen=2/LD

Rd= C1+C2*xD

E1=Rd.subs({xD:LD})
E2=Rd.subs({xD:LD/2})+QD


sol= sy.linsolve([E1,E2],[C1,C2])


C1sol=sol.args[0][0]
C2sol=sol.args[0][1]

 
Rd=Rd.subs({C1:C1sol,C2:C2sol})


qD= Rd*cos*cos 
pD=Rd*cos*sen

LE=(5**2+2**2)**0.5

thetaE=sy.atan2(-2,5)



Qe= C1e+C2e*xE

E3= Qe.subs({xE:0})+QE
E4=Qe.subs({xE:LE})

sol1= sy.linsolve([E3,E4],[C1e,C2e])

C1esol=sol1.args[0][0]
C2esol=sol1.args[0][1]



Qe=Qe.subs({C1e:C1esol,C2e:C2esol})
pE=Qe*sy.Abs(sy.sin(thetaE))*sy.cos(thetaE)

qE=-Qe*sy.Abs(sy.sin(thetaE))*sy.sin(thetaE)


# fUNCIONES DE FORMA 

pi1=(1-x/L)
pi2=(1-3*(x/L)**2+2*(x/L)**3)
pi3=((x/L-2*(x/L)**2+(x/L)**3)*L)
pi4=x/L
pi5=3*(x/L)**2-2*(x/L)**3
pi6= (-(x/L)**2+(x/L)**3)*L 

# Funcionde forma del Elemento D

pi1D=pi1.subs({x:xD,L:LD})
pi2D=pi2.subs({x:xD,L:LD})
pi3D=pi3.subs({x:xD,L:LD})
pi4D=pi4.subs({x:xD,L:LD})
pi5D=pi5.subs({x:xD,L:LD})
pi6D=pi6.subs({x:xD,L:LD})

# Funciondes de forma del elemento E

pi1E=pi1.subs({x:xE,L:LE})
pi2E=pi2.subs({x:xE,L:LE})
pi3E=pi3.subs({x:xE,L:LE})
pi4E=pi4.subs({x:xE,L:LE})
pi5E=pi5.subs({x:xE,L:LE})
pi6E=pi6.subs({x:xE,L:LE})

# Funciones de forma del elemento axial A 

pi1A=pi1.subs({x:xA,L:LA})
pi4A=pi4.subs({x:xA,L:LA})

# Funciones de forma del elemento axial B 


pi1B=pi1.subs({x:xB,L:LB})
pi4B=pi4.subs({x:xB,L:LB})

# Funciondes de forma del elemento C

pi1C=pi1.subs({x:xE,L:LC})
pi2C=pi2.subs({x:xE,L:LC})
pi3C=pi3.subs({x:xE,L:LC})
pi4C=pi4.subs({x:xE,L:LC})
pi5C=pi5.subs({x:xE,L:LC})
pi6C=pi6.subs({x:xE,L:LC})


# funciones de Green 

Gxx1=L/(AE)*pi4*pi1.subs({x:xi})
Gxx2=L/(AE)*pi1*pi4.subs({x:xi})

Gyy1=((L**3)/(6*EI))*(-(x/L)**3*pi2.subs({x:xi})+3*(x/L)**2*pi3.subs({x:xi})/L)
Gyy2=((L**3)/(6*EI))*(-(1-x/L)**3*pi5.subs({x:xi})-3*(1-x/L)**2*pi6.subs({x:xi})/L)

# Elemento E 

Gxx1E=Gxx1.subs({x:xE,xi:xiE,L:LE})
Gxx2E=Gxx2.subs({x:xE,xi:xiE,L:LE})

Gyy1E=Gyy1.subs({x:xE,xi:xiE,L:LE})
Gyy2E=Gyy2.subs({x:xE,xi:xiE,L:LE})

# Elementos D 


Gxx1D=Gxx1.subs({x:xD,xi:xiD,L:LD})
Gxx2D=Gxx2.subs({x:xD,xi:xiD,L:LD})

Gyy1D=Gyy1.subs({x:xD,xi:xiD,L:LD})
Gyy2D=Gyy2.subs({x:xD,xi:xiD,L:LD})

# Seciones tranversales de elementos 

BaseC=0.55   # base del elemento C 
AlturaC=0.25  # Altura del eleemtno C

BaseD=0.45    # base del elemento D 
AlturaD=0.35   # Altura del eleemtno D

BaseE=0.45     # base del elemento E 
AlturaE=0.35    # Altura del eleemtno E


#  2.1 Elemento A (Pila izquierda)
    #  2.1.1 Formulacion en coordenadas locales
A=np.pi*R**2     # AREA de la seccion tranversal del elemento
I=np.pi/4*R**4
L=LA
AE=A*E
EI=E*I
ThetaA=np.pi/2
L=LA
kA=2*R*k
Lambda=(0.25*kA/EI)**0.25  #[1/L]

s=np.sin(Lambda*L)
c=np.cos(Lambda*L)
sh=np.sinh(Lambda*L)
ch=np.cosh(Lambda*L)

k22=4*EI*Lambda**3*(s*c+sh*ch)/(sh**2-s**2)
k23=2*EI*Lambda**2*(s**2+sh**2)/(sh**2-s**2)
k25=-4*EI*Lambda**3*(s*ch+c*sh)/(sh**2-s**2)
k26=4*EI*Lambda**2*s*sh/(sh**2-s**2)
k33=2*EI*Lambda*(sh*ch-s*c)/(sh**2-s**2)
k36=2*EI*Lambda*(s*ch-c*sh)/(sh**2-s**2)

kALoc=np.zeros([6,6])

kALoc[0,0]=AE/L
kALoc[1,0]=0
kALoc[2,0]=0
kALoc[3,0]=-AE/L
kALoc[4,0]=0
kALoc[5,0]=0

kALoc[0,1]=0
kALoc[1,1]=k22
kALoc[2,1]=k23
kALoc[3,1]=0
kALoc[4,1]=k25
kALoc[5,1]=k26

kALoc[0,2]=0
kALoc[1,2]=k23
kALoc[2,2]=k33
kALoc[3,2]=0
kALoc[4,2]=-k26
kALoc[5,2]=k36

kALoc[0,3]=-AE/L
kALoc[1,3]=0
kALoc[2,3]=0
kALoc[3,3]=AE/L
kALoc[4,3]=0
kALoc[5,3]=0

kALoc[0,4]=0
kALoc[1,4]=k25
kALoc[2,4]=-k26
kALoc[3,4]=0
kALoc[4,4]=k22
kALoc[5,4]=-k23

kALoc[0,5]=0
kALoc[1,5]=k26
kALoc[2,5]=k36
kALoc[3,5]=0
kALoc[4,5]=-k23
kALoc[5,5]=k33 

# Matriz de tranformacion  para el elemento A 

TA=np.zeros([6,6])

TA[0,0]= np.cos(ThetaA)
TA[1,0]=-np.sin(ThetaA)

TA[0,1]=np.sin(ThetaA)
TA[1,1]=np.cos(ThetaA)

TA[2,2]=1

TA[3,3]= np.cos(ThetaA)
TA[4,3]=-np.sin(ThetaA)

TA[3,4]=np.sin(ThetaA)
TA[4,4]=np.cos(ThetaA)

TA[5,5]=1

kAGlo=np.matrix.transpose(TA)@kALoc@TA

##  2.2 Elemento B (Pila izquierda) 

L=LB
A=np.pi*R**2
I=np.pi/4*R**4
AE=A*E
EI=E*I
ThetaB=np.pi/2
kB=2*R*k

Lambda=(0.25*kB/EI)**0.25  #[1/L]

s=np.sin(Lambda*L)
c=np.cos(Lambda*L)
sh=np.sinh(Lambda*L)
ch=np.cosh(Lambda*L)

k22=4*EI*Lambda**3*(s*c+sh*ch)/(sh**2-s**2)
k23=2*EI*Lambda**2*(s**2+sh**2)/(sh**2-s**2)
k25=-4*EI*Lambda**3*(s*ch+c*sh)/(sh**2-s**2)
k26=4*EI*Lambda**2*s*sh/(sh**2-s**2)
k33=2*EI*Lambda*(sh*ch-s*c)/(sh**2-s**2)
k36=2*EI*Lambda*(s*ch-c*sh)/(sh**2-s**2)

kBLoc=np.zeros([6,6])

kBLoc[0,0]=AE/L
kBLoc[1,0]=0
kBLoc[2,0]=0
kBLoc[3,0]=-AE/L
kBLoc[4,0]=0
kBLoc[5,0]=0

kBLoc[0,1]=0
kBLoc[1,1]=k22
kBLoc[2,1]=k23
kBLoc[3,1]=0
kBLoc[4,1]=k25
kBLoc[5,1]=k26

kBLoc[0,2]=0
kBLoc[1,2]=k23
kBLoc[2,2]=k33
kBLoc[3,2]=0
kBLoc[4,2]=-k26
kBLoc[5,2]=k36

kBLoc[0,3]=-AE/L
kBLoc[1,3]=0
kBLoc[2,3]=0
kBLoc[3,3]=AE/L
kBLoc[4,3]=0
kBLoc[5,3]=0

kBLoc[0,4]=0
kBLoc[1,4]=k25
kBLoc[2,4]=-k26
kBLoc[3,4]=0
kBLoc[4,4]=k22
kBLoc[5,4]=-k23

kBLoc[0,5]=0
kBLoc[1,5]=k26
kBLoc[2,5]=k36
kBLoc[3,5]=0
kBLoc[4,5]=-k23
kBLoc[5,5]=k33

#  2.4.2 Formulacion en coordenadas globales

TB=np.zeros([6,6])

TB[0,0]= np.cos(ThetaB)
TB[1,0]=-np.sin(ThetaB)

TB[0,1]=np.sin(ThetaB)
TB[1,1]=np.cos(ThetaB)

TB[2,2]=1

TB[3,3]= np.cos(ThetaB)
TB[4,3]=-np.sin(ThetaB)

TB[3,4]=np.sin(ThetaB)
TB[4,4]=np.cos(ThetaB)

TB[5,5]=1

kBGlo=np.matrix.transpose(TB)@kBLoc@TB

# Elemtnos C en coordenadas ya globales 

L=LC
A=AlturaC*BaseC
I=BaseC*AlturaC**3/12
AE=A*E
EI=E*I
kC=BaseC*k

Lambda=(0.25*kC/EI)**0.25  #[1/L]

s=np.sin(Lambda*L)
c=np.cos(Lambda*L)
sh=np.sinh(Lambda*L)
ch=np.cosh(Lambda*L)

k22=4*EI*Lambda**3*(s*c+sh*ch)/(sh**2-s**2)
k23=2*EI*Lambda**2*(s**2+sh**2)/(sh**2-s**2)
k25=-4*EI*Lambda**3*(s*ch+c*sh)/(sh**2-s**2)
k26=4*EI*Lambda**2*s*sh/(sh**2-s**2)
k33=2*EI*Lambda*(sh*ch-s*c)/(sh**2-s**2)
k36=2*EI*Lambda*(s*ch-c*sh)/(sh**2-s**2)

kCLoc=np.zeros([6,6])

kCLoc[0,0]=AE/L
kCLoc[1,0]=0
kCLoc[2,0]=0
kCLoc[3,0]=-AE/L
kCLoc[4,0]=0
kCLoc[5,0]=0

kCLoc[0,1]=0
kCLoc[1,1]=k22
kCLoc[2,1]=k23
kCLoc[3,1]=0
kCLoc[4,1]=k25
kCLoc[5,1]=k26

kCLoc[0,2]=0
kCLoc[1,2]=k23
kCLoc[2,2]=k33
kCLoc[3,2]=0
kCLoc[4,2]=-k26
kCLoc[5,2]=k36

kCLoc[0,3]=-AE/L
kCLoc[1,3]=0
kCLoc[2,3]=0
kCLoc[3,3]=AE/L
kCLoc[4,3]=0
kCLoc[5,3]=0

kCLoc[0,4]=0
kCLoc[1,4]=k25
kCLoc[2,4]=-k26
kCLoc[3,4]=0
kCLoc[4,4]=k22
kCLoc[5,4]=-k23

kCLoc[0,5]=0
kCLoc[1,5]=k26
kCLoc[2,5]=k36
kCLoc[3,5]=0
kCLoc[4,5]=-k23
kCLoc[5,5]=k33 

# Elemtnos D en coordenadas ya globales 

L=LD
A=BaseD*AlturaD
I=BaseD*AlturaD**3/12
AE=A*E
EI=E*I

ThetaD=np.arctan2(2,4)

kDLoc=np.zeros([6,6])

kDLoc[0,0]=AE/L
kDLoc[1,0]=0
kDLoc[2,0]=0
kDLoc[3,0]=-AE/L
kDLoc[4,0]=0
kDLoc[5,0]=0          

kDLoc[0,1]=0
kDLoc[1,1]=12*EI/L**3
kDLoc[2,1]=6*EI/L**2
kDLoc[3,1]=0
kDLoc[4,1]=-12*EI/L**3
kDLoc[5,1]=6*EI/L**2

kDLoc[0,2]=0
kDLoc[1,2]=6*EI/L**2
kDLoc[2,2]=4*EI/L
kDLoc[3,2]=0    
kDLoc[4,2]=-6*EI/L**2
kDLoc[5,2]=2*EI/L

kDLoc[0,3]=-AE/L
kDLoc[1,3]=0
kDLoc[2,3]=0
kDLoc[3,3]=AE/L
kDLoc[4,3]=0
kDLoc[5,3]=0  

kDLoc[0,4]=0  
kDLoc[1,4]=-12*EI/L**3
kDLoc[2,4]=-6*EI/L**2
kDLoc[3,4]=0  
kDLoc[4,4]=12*EI/L**3
kDLoc[5,4]=-6*EI/L**2

kDLoc[0,5]=0
kDLoc[1,5]=6*EI/L**2
kDLoc[2,5]=2*EI/L
kDLoc[3,5]=0
kDLoc[4,5]=-6*EI/L**2
kDLoc[5,5]=4*EI/L

FEmpDLoc= sy.zeros(6,1)


FEmpDLoc[0,0]=-sy.integrate(pD*pi1D,(xD,LD/2,LD))
FEmpDLoc[1,0]=-sy.integrate(qD*pi2D,(xD,LD/2,LD))
FEmpDLoc[2,0]=-sy.integrate(qD*pi3D,(xD,LD/2,LD))
FEmpDLoc[3,0]=-sy.integrate(pD*pi4D,(xD,LD/2,LD))
FEmpDLoc[4,0]=-sy.integrate(qD*pi5D,(xD,LD/2,LD))
FEmpDLoc[5,0]=-sy.integrate(qD*pi6D,(xD,LD/2,LD)) 


#  2.5.2 Formulacion en coordenadas globales

TD=np.zeros([6,6])

TD[0,0]= np.cos(ThetaD)
TD[1,0]=-np.sin(ThetaD)

TD[0,1]=np.sin(ThetaD)
TD[1,1]=np.cos(ThetaD)

TD[2,2]=1

TD[3,3]= np.cos(ThetaD)
TD[4,3]=-np.sin(ThetaD)

TD[3,4]=np.sin(ThetaD)
TD[4,4]=np.cos(ThetaD)

TD[5,5]=1

kDGlo=np.matrix.transpose(TD)@kDLoc@TD

FEmDGlo=np.matrix.transpose(TD)@FEmpDLoc


# Elemtnos E en coordenadas globales 

L=LE
A=AlturaE*BaseE
I=BaseE*AlturaE**3/12
AE=A*E
EI=E*I

ThetaE=np.arctan2(-2,5)

kEGLoc=np.zeros([6,6])

kEGLoc[0,0]=AE/L
kEGLoc[1,0]=0
kEGLoc[2,0]=0
kEGLoc[3,0]=-AE/L
kEGLoc[4,0]=0
kEGLoc[5,0]=0          

kEGLoc[0,1]=0
kEGLoc[1,1]=12*EI/L**3
kEGLoc[2,1]=6*EI/L**2
kEGLoc[3,1]=0
kEGLoc[4,1]=-12*EI/L**3
kEGLoc[5,1]=6*EI/L**2

kEGLoc[0,2]=0
kEGLoc[1,2]=6*EI/L**2
kEGLoc[2,2]=4*EI/L
kEGLoc[3,2]=0    
kEGLoc[4,2]=-6*EI/L**2
kEGLoc[5,2]=2*EI/L

kEGLoc[0,3]=-AE/L
kEGLoc[1,3]=0
kEGLoc[2,3]=0
kEGLoc[3,3]=AE/L
kEGLoc[4,3]=0
kEGLoc[5,3]=0  

kEGLoc[0,4]=0  
kEGLoc[1,4]=-12*EI/L**3
kEGLoc[2,4]=-6*EI/L**2
kEGLoc[3,4]=0  
kEGLoc[4,4]=12*EI/L**3
kEGLoc[5,4]=-6*EI/L**2

kEGLoc[0,5]=0
kEGLoc[1,5]=6*EI/L**2
kEGLoc[2,5]=2*EI/L
kEGLoc[3,5]=0
kEGLoc[4,5]=-6*EI/L**2
kEGLoc[5,5]=4*EI/L

FEmpELoc= sy.zeros(6,1)


FEmpELoc[0,0]=-sy.integrate(pE*pi1E,(xE,0,LE))
FEmpELoc[1,0]=-sy.integrate(qE*pi2E,(xE,0,LE))
FEmpELoc[2,0]=-sy.integrate(qE*pi3E,(xE,0,LE))
FEmpELoc[3,0]=-sy.integrate(pE*pi4E,(xE,0,LE))
FEmpELoc[4,0]=-sy.integrate(qE*pi5E,(xE,0,LE))
FEmpELoc[5,0]=-sy.integrate(qE*pi6E,(xE,0,LE)) 

#  2.5.2 Formulacion en coordenadas globales

TE=np.zeros([6,6])

TE[0,0]= np.cos(ThetaE)
TE[1,0]=-np.sin(ThetaE)

TE[0,1]=np.sin(ThetaE)
TE[1,1]=np.cos(ThetaE)

TE[2,2]=1

TE[3,3]= np.cos(ThetaE)
TE[4,3]=-np.sin(ThetaE)

TE[3,4]=np.sin(ThetaE)
TE[4,4]=np.cos(ThetaE)

TE[5,5]=1

kEGlo=np.matrix.transpose(TE)@kEGLoc@TE

FEmEGlo=np.matrix.transpose(TE)@FEmpELoc


# Calculo de los desplazamiento nodales 

Kdes=np.zeros([10,10])

Kdes[0,0]=kAGlo[3,3]+kCLoc[0,0]+kDGlo[0,0]
Kdes[1,0]=kAGlo[4,3]+kCLoc[1,0]+kDGlo[1,0]
Kdes[2,0]=kAGlo[5,3]+kCLoc[2,0]+kDGlo[2,0]
Kdes[3,0]=kDGlo[3,0]
Kdes[4,0]=kDGlo[4,0]
Kdes[5,0]=kDGlo[5,0]
Kdes[7,0]=kCLoc[3,0]
Kdes[8,0]=kCLoc[4,0]
Kdes[9,0]=kCLoc[5,0]


Kdes[0,1]=kAGlo[3,4]+kCLoc[0,1]+kDGlo[0,1]
Kdes[1,1]=kAGlo[4,4]+kCLoc[1,1]+kDGlo[1,1]
Kdes[2,1]=kAGlo[5,4]+kCLoc[2,1]+kDGlo[2,1]
Kdes[3,1]=kDGlo[3,1]
Kdes[4,1]=kDGlo[4,1]
Kdes[5,1]=kDGlo[5,1]
Kdes[7,1]=kCLoc[3,1]
Kdes[8,1]=kCLoc[4,1]
Kdes[9,1]=kCLoc[5,1]


Kdes[0,2]=kAGlo[3,5]+kCLoc[0,2]+kDGlo[0,2]
Kdes[1,2]=kAGlo[4,5]+kCLoc[1,2]+kDGlo[1,2]
Kdes[2,2]=kAGlo[5,5]+kCLoc[2,2]+kDGlo[2,2]
Kdes[3,2]=kDGlo[3,2]
Kdes[4,2]=kDGlo[4,2]
Kdes[5,2]=kDGlo[5,2]
Kdes[7,2]=kCLoc[3,2]
Kdes[8,2]=kCLoc[4,2]
Kdes[9,2]=kCLoc[5,2]


Kdes[0,3]=kDGlo[0,3]
Kdes[1,3]=kDGlo[1,3]
Kdes[2,3]=kDGlo[2,3]
Kdes[3,3]=kDGlo[3,3]+kEGlo[0,0]
Kdes[4,3]=kDGlo[4,3]+kEGlo[1,0]
Kdes[5,3]=kDGlo[5,3]
Kdes[6,3]=kEGlo[2,0]
Kdes[7,3]=kEGlo[3,0]
Kdes[8,3]=kEGlo[4,0]
Kdes[9,3]=kEGlo[5,0]


Kdes[0,4]=kDGlo[0,4]
Kdes[1,4]=kDGlo[1,4]
Kdes[2,4]=kDGlo[2,4]
Kdes[3,4]=kDGlo[3,4]+kEGlo[0,1]
Kdes[4,4]=kDGlo[4,4]+kEGlo[1,1]
Kdes[5,4]=kDGlo[5,4]
Kdes[6,4]=kEGlo[2,1]
Kdes[7,4]=kEGlo[3,1]
Kdes[8,4]=kEGlo[4,1]
Kdes[9,4]=kEGlo[5,1]



Kdes[0,5]=kDGlo[0,5]
Kdes[1,5]=kDGlo[1,5]
Kdes[2,5]=kDGlo[2,5]
Kdes[3,5]=kDGlo[3,5]
Kdes[4,5]=kDGlo[4,5]
Kdes[5,5]=kDGlo[5,5]


Kdes[3,6]=kEGlo[0,2]
Kdes[4,6]=kEGlo[1,2]
Kdes[6,6]=kEGlo[2,2]
Kdes[7,6]=kEGlo[3,2]
Kdes[8,6]=kEGlo[4,2]
Kdes[9,6]=kEGlo[5,2]


Kdes[0,7]=kCLoc[0,3]
Kdes[1,7]=kCLoc[1,3]
Kdes[2,7]=kCLoc[2,3]
Kdes[3,7]=kEGlo[0,3]
Kdes[4,7]=kEGlo[1,3]
Kdes[6,7]=kEGlo[2,3]
Kdes[7,7]=kBGlo[3,3]+kCLoc[3,3]+kEGlo[3,3]
Kdes[8,7]=kBGlo[4,3]+kCLoc[4,3]+kEGlo[4,3]
Kdes[9,7]=kBGlo[5,3]+kCLoc[5,3]+kEGlo[5,3]


Kdes[0,8]=kCLoc[0,4]
Kdes[1,8]=kCLoc[1,4]
Kdes[2,8]=kCLoc[2,4]
Kdes[3,8]=kEGlo[0,4]
Kdes[4,8]=kEGlo[1,4]
Kdes[6,8]=kEGlo[2,4]
Kdes[7,8]=kBGlo[3,4]+kCLoc[3,4]+kEGlo[3,4]
Kdes[8,8]=kBGlo[4,4]+kCLoc[4,4]+kEGlo[4,4]
Kdes[9,8]=kBGlo[5,4]+kCLoc[5,4]+kEGlo[5,4]


Kdes[0,9]=kCLoc[0,5]
Kdes[1,9]=kCLoc[1,5]
Kdes[2,9]=kCLoc[2,5]
Kdes[3,9]=kEGlo[0,5]
Kdes[4,9]=kEGlo[1,5]
Kdes[6,9]=kEGlo[2,5]
Kdes[7,9]=kBGlo[3,5]+kCLoc[3,5]+kEGlo[3,5]
Kdes[8,9]=kBGlo[4,5]+kCLoc[4,5]+kEGlo[4,5]
Kdes[9,9]=kBGlo[5,5]+kCLoc[5,5]+kEGlo[5,5]

kdes=Kdes-np.transpose(Kdes)

# VECTOR DE FUERZAS DE EMPOTRAMIENTO 

FEmpDes=np.zeros([10,1])


FEmpDes[0,0]= FEmDGlo[0]
FEmpDes[1,0]= FEmDGlo[1]
FEmpDes[2,0]= FEmDGlo[2]
FEmpDes[3,0]= FEmDGlo[3]+FEmEGlo[0]
FEmpDes[4,0]= FEmDGlo[4]+FEmEGlo[1]
FEmpDes[5,0]= FEmDGlo[5]
FEmpDes[6,0]= FEmEGlo[2]
FEmpDes[7,0]= FEmEGlo[3]
FEmpDes[8,0]= FEmEGlo[4]
FEmpDes[9,0]= FEmEGlo[5]


# Calculo de desplazamiento nodales en coordenadas Globlales 

DesNod=np.linalg.solve(Kdes,-FEmpDes)


u2=DesNod[ 0,0]
v2=DesNod[ 1,0]
theta2=DesNod[ 2,0]
u3=DesNod[ 3,0]
v3=DesNod[ 4,0]
theta3D=DesNod[ 5,0]
theta3E=DesNod[ 6,0]
u4=DesNod[ 7,0]
v4=DesNod[ 8,0]
theta4=DesNod[ 9,0]

u1=0
v1=0
theta1=0
u5=0
v5=0
theta5=0

# Calculo de reacciones 

KRea= np.zeros([6,10]) 


KRea[0,0]= kAGlo[0,3]
KRea[0,1]= kAGlo[0,4]
KRea[0,2]= kAGlo[0,5] 

KRea[1,0]= kAGlo[1,3]
KRea[1,1]= kAGlo[1,4]
KRea[1,2]= kAGlo[1,5] 


KRea[2,0]= kAGlo[2,3]
KRea[2,1]= kAGlo[2,4]
KRea[2,2]= kAGlo[2,5] 


KRea[3,7]= kBGlo[0,3]
KRea[3,8]= kBGlo[0,4]
KRea[3,9]= kBGlo[0,5] 


KRea[4,7]= kBGlo[1,3]
KRea[4,8]= kBGlo[1,4]
KRea[4,9]= kBGlo[1,5] 

KRea[5,7]= kBGlo[2,3]
KRea[5,8]= kBGlo[2,4]
KRea[5,9]= kBGlo[2,5] 


Rea=KRea@DesNod  

FX1=Rea[0,0]
FY1=Rea[1,0]
M1=Rea[2,0]
FX5=Rea[3,0]
FY5=Rea[4,0]
M5=Rea[5,0]


#  5.1 Elemento A 

[ui,vi,qi,uj,vj,qj]=TA@[u1,v1,theta1,u2,v2,theta2]


A=np.pi*R**2     # AREA de la seccion tranversal del elemento
I=(np.pi*R**4)/4
L=LA
AEA=A*E
EIA=E*I
ThetaA=np.pi/2
L=LA
kA=2*R*k


LambdaA=(kA/(4*EIA))**0.25

s =np.sin(LambdaA*LA)
c =np.cos(LambdaA*LA)
sh=np.sinh(LambdaA*LA)
ch=np.cosh(LambdaA*LA)

# 2. Preliminares

# 2.1 Funciones de forma 


pi1A=(1-xA/LA)
pi4A=xA/LA

Psi2A=(-(s**2*ch**2+c**2*sh**2)*sy.sin(LambdaA*xA)*sy.sinh(LambdaA*xA)+(s*c+sh*ch)*sy.sin(LambdaA*xA)*sy.cosh(LambdaA*xA)-(s*c+sh*ch)*sy.cos(LambdaA*xA)*sy.sinh(LambdaA*xA)+(sh**2-s**2)*sy.cos(LambdaA*xA)*sy.cosh(LambdaA*xA))/(sh**2-s**2)
Psi3A=(1/LambdaA)*((s*c-sh*ch)*sy.sin(LambdaA*xA)*sy.sinh(LambdaA*xA)+sh**2*sy.sin(LambdaA*xA)*sy.cosh(LambdaA*xA)-s**2*sy.cos(LambdaA*xA)*sy.sinh(LambdaA*xA))/(sh**2-s**2)
Psi5A=(2*s*sh*sy.sin(LambdaA*xA)*sy.sinh(LambdaA*xA)-(s*ch+c*sh)*sy.sin(LambdaA*xA)*sy.cosh(LambdaA*xA)+(s*ch+c*sh)*sy.cos(LambdaA*xA)*sy.sinh(LambdaA*xA))/(sh**2-s**2)
Psi6A=(1/LambdaA)*((c*sh-s*ch)*sy.sin(LambdaA*xA)*sy.sinh(LambdaA*xA)+s*sh*sy.sin(LambdaA*xA)*sy.cosh(LambdaA*xA)-s*sh*sy.cos(LambdaA*xA)*sy.sinh(LambdaA*xA))/(sh**2-s**2)

# Desplazamiento axiales 

uA= pi1A*ui+pi4A*uj

# Desplazamiento  perpendiculares al elemento 

vA= Psi2A*vi+Psi3A*qi+Psi5A*vj+Psi6A*qj


fA=-kA*vA

MA=EIA*sy.diff(vA,xA,2)

VA=-EIA*sy.diff(vA,xA,3)

PA=AEA*sy.diff(uA,xA,1)


#  5.1 Elemento B


[ui,vi,qi,uj,vj,qj]=TB@[u5,v5,theta5,u4,v4,theta4]


A=np.pi*R**2     # AREA de la seccion tranversal del elemento
I=(np.pi*R**4)/4
L=LB
AEB=A*E
EIB=E*I
ThetaB=np.pi/2
L=LB
kB=2*R*k


LambdaB=(kB/(4*EIB))**0.25

s =np.sin(LambdaB*LB)
c =np.cos(LambdaB*LB)
sh=np.sinh(LambdaB*LB)
ch=np.cosh(LambdaB*LB)

# 2. Preliminares

# 2.1 Funciones de forma


pi1B=(1-xB/LB)
pi4B=xB/LB

Psi2B=(-(s**2*ch**2+c**2*sh**2)*sy.sin(LambdaB*xB)*sy.sinh(LambdaB*xB)+(s*c+sh*ch)*sy.sin(LambdaB*xB)*sy.cosh(LambdaB*xB)-(s*c+sh*ch)*sy.cos(LambdaB*xB)*sy.sinh(LambdaB*xB)+(sh**2-s**2)*sy.cos(LambdaB*xB)*sy.cosh(LambdaB*xB))/(sh**2-s**2)
Psi3B=(1/LambdaB)*((s*c-sh*ch)*sy.sin(LambdaB*xB)*sy.sinh(LambdaB*xB)+sh**2*sy.sin(LambdaB*xB)*sy.cosh(LambdaB*xB)-s**2*sy.cos(LambdaB*xB)*sy.sinh(LambdaB*xB))/(sh**2-s**2)
Psi5B=(2*s*sh*sy.sin(LambdaB*xB)*sy.sinh(LambdaB*xB)-(s*ch+c*sh)*sy.sin(LambdaB*xB)*sy.cosh(LambdaB*xB)+(s*ch+c*sh)*sy.cos(LambdaB*xB)*sy.sinh(LambdaB*xB))/(sh**2-s**2)
Psi6B=(1/LambdaB)*((c*sh-s*ch)*sy.sin(LambdaB*xB)*sy.sinh(LambdaB*xB)+s*sh*sy.sin(LambdaB*xB)*sy.cosh(LambdaB*xB)-s*sh*sy.cos(LambdaA*xB)*sy.sinh(LambdaA*xB))/(sh**2-s**2)


# Des plazamiento axiales 

uB= pi1B*ui+pi4B*uj

#Desplazamiento perpendiculares al elemento 

vB= Psi2B*vi+Psi3B*qi+Psi5B*vj+Psi6B*qj


# Cortante momento fuerzas axial 


VB=-EIB*sy.diff(vB,xB,3)
MB= EIB*sy.diff(vB,xB,2)
PB= AEB*sy.diff(uB,xB,1)

fB=-kB*vB


#  5.1 Elemento C


[ui,vi,qi,uj,vj,qj]=[u2,v2,theta2,u4,v4,theta4]


L=LC
AC=AlturaC*BaseC
I=BaseC*AlturaC**3/12
AEC=AC*E
EIC=E*I
kC=BaseC*k


LambdaC=(kC/(4*EIC))**0.25

s =np.sin(LambdaC*LC)
c =np.cos(LambdaC*LC)
sh=np.sinh(LambdaC*LC)
ch=np.cosh(LambdaC*LC)

# 2. Preliminares

# 2.1 Funciones de forma


pi1C=(1-xC/LC)
pi4C=xC/LC


Psi2C=(-(s**2*ch**2+c**2*sh**2)*sy.sin(LambdaC*xC)*sy.sinh(LambdaC*xC)+(s*c+sh*ch)*sy.sin(LambdaC*xC)*sy.cosh(LambdaC*xC)-(s*c+sh*ch)*sy.cos(LambdaC*xC)*sy.sinh(LambdaC*xC)+(sh**2-s**2)*sy.cos(LambdaC*xC)*sy.cosh(LambdaC*xC))/(sh**2-s**2)
Psi3C=(1/LambdaC)*((s*c-sh*ch)*sy.sin(LambdaC*xC)*sy.sinh(LambdaC*xC)+sh**2*sy.sin(LambdaC*xC)*sy.cosh(LambdaC*xC)-s**2*sy.cos(LambdaC*xC)*sy.sinh(LambdaC*xC))/(sh**2-s**2)
Psi5C=(2*s*sh*sy.sin(LambdaC*xC)*sy.sinh(LambdaC*xC)-(s*ch+c*sh)*sy.sin(LambdaC*xC)*sy.cosh(LambdaC*xC)+(s*ch+c*sh)*sy.cos(LambdaC*xC)*sy.sinh(LambdaC*xC))/(sh**2-s**2)
Psi6C=(1/LambdaC)*((c*sh-s*ch)*sy.sin(LambdaC*xC)*sy.sinh(LambdaC*xC)+s*sh*sy.sin(LambdaC*xC)*sy.cosh(LambdaC*xC)-s*sh*sy.cos(LambdaC*xC)*sy.sinh(LambdaC*xC))/(sh**2-s**2)


# Des plazamiento axiales 

uC= pi1C*ui+pi4C*uj

# Desplazamiento perpendiculaes al elemento 

vC= Psi2C*vi+Psi3C*qi+Psi5C*vj+Psi6C*qj

fC=-kC*vC

# Momentos , cortante y axial 


VC=-EIC*sy.diff(vC,xC,3)
MC=EIC*sy.diff(vC,xC,2)
PC=AEC*sy.diff(uC,xC,1)


Sx= FX1+FX5-sy.integrate(fA,(xA,0,LA))-sy.integrate(fB,(xB,0,LB))-30

Sy=FY1+FY5+sy.integrate(fC,(xC,0,LC))-60


#Elemento D 

L=LD
A=AlturaD*BaseD
I=BaseD*AlturaD**3/12
AED=A*E
EID=E*I

# funciones de Green 



pi1D=(1-xD/LD)
pi4D=xD/LD

pi2=(1-3*(x/L)**2+2*(x/L)**3)
pi3=((x/L-2*(x/L)**2+(x/L)**3)*L)
pi5=3*(x/L)**2-2*(x/L)**3
pi6= (-(x/L)**2+(x/L)**3)*L 


pi2D=pi2.subs({x:xD,L:LD})
pi3D=pi3.subs({x:xD,L:LD})
pi5D=pi5.subs({x:xD,L:LD})
pi6D=pi6.subs({x:xD,L:LD})


# funciones de Green 


Gxx1=L/(AED)*pi4*pi1.subs({x:xi})
Gxx2=L/(AED)*pi1*pi4.subs({x:xi})


Gyy1=(sy.Rational((L**3),(6*EID)))*(-(x/L)**3*pi2.subs({x:xi})+3*(x/L)**2*pi3.subs({x:xi})/L)
Gyy2=(sy.Rational((L**3),(6*EID)))*(-(1-x/L)**3*pi5.subs({x:xi})-3*(1-x/L)**2*pi6.subs({x:xi})/L)

L=sy.symbols('L')

Gxx1D=Gxx1.subs({pi4:pi4D,pi1:pi1D,L:LD,x:xD,xi:xiD})
Gxx2D=Gxx1.subs({pi4:pi4D,pi1:pi1D,L:LD,x:xD,xi:xiD})

Gyy1D=Gyy1.subs({pi2:pi2D,pi3:pi3D,L:LD,x:xD,xi:xiD})
Gyy2D=Gyy2.subs({pi5:pi5D,pi6:pi6D,L:LD,x:xD,xi:xiD})

[ui,vi,qi,uj,vj,qj]=TD@[u2,v2,theta2,u3,v3,theta3D]

# Desplazamiento axiales 

uDh= ui*pi1D+uj*pi4D

uDf1= sy.integrate(Gxx1D*pD.subs({xD:xiD}),(xiD,LD/2,LD))

uDf2= sy.integrate(Gxx2D*pD.subs({xD:xiD}),(xiD,LD/2,xD))+sy.integrate(Gxx1D*pD.subs({xD:xiD}),(xiD,xD,LD))

uD1= uDh+uDf1
uD2= uDh+uDf2 

# Desplazameitno perpendiculares al elemento 

vDh= vi*pi2D+qi*pi3D+vj*pi5D+qj*pi6D

vDf1= sy.integrate(Gyy1D*qD.subs({xD:xiD}),(xiD,LD/2,LD))

vDf2= sy.integrate(Gyy2D*qD.subs({xD:xiD}),(xiD,LD/2,xD))+sy.integrate(Gyy1D*qD.subs({xD:xiD}),(xiD,xD,LD))

vD1= vDh+vDf1
vD2= vDh +vDf2

# Carga axial 

PD1= AED*sy.diff(uD1,xD,1)
PD2= AED*sy.diff(uD2,xD,1)


#Momentos del elemento D 
MD1= EID*sy.diff(vD1,xD,2 )
MD2= EID*sy.diff(vD2,xD,2 )


# Cortante del elemento D 

VD1=-EID*sy.diff(vD1,xD,3 )
VD2=-EID*sy.diff(vD2,xD,3 )


#Elemento E

L=LE
A=AlturaE*BaseE
I=BaseE*AlturaE**3/12
AEE=A*E
EIE=E*I


# funciones de Green 


pi1E=(1-xE/LE)
pi4E=xE/LE

pi2=(1-3*(x/L)**2+2*(x/L)**3)
pi3=((x/L-2*(x/L)**2+(x/L)**3)*L)
pi5=3*(x/L)**2-2*(x/L)**3
pi6= (-(x/L)**2+(x/L)**3)*L 


pi2E=pi2.subs({x:xE,L:LE})
pi3E=pi3.subs({x:xE,L:LE})
pi5E=pi5.subs({x:xE,L:LE})
pi6E=pi6.subs({x:xE,L:LE})


# funciones de Green 

L=sy.symbols('L')

Gxx1=L/(AEE)*pi4*pi1.subs({x:xi})
Gxx2=L/(AEE)*pi1*pi4.subs({x:xi})


Gyy1=(sy.Rational((LE**3),(6*EIE)))*(-(x/L)**3*pi2.subs({x:xi})+3*(x/L)**2*pi3.subs({x:xi})/LE)
Gyy2=(sy.Rational((LE**3),(6*EIE)))*(-(1-x/L)**3*pi5.subs({x:xi})-3*(1-x/L)**2*pi6.subs({x:xi})/LE)



Gxx1E=Gxx1.subs({L:LE,x:xE,xi:xiE})
Gxx2E=Gxx2.subs({pi4:pi4E,pi1:pi1E,L:LE,x:xE,xi:xiE})



Gyy1E=Gyy1.subs({pi2:pi2E,pi3:pi3E,L:LE,x:xE,xi:xiE})
Gyy2E=Gyy2.subs({pi5:pi5E,pi6:pi6E,L:LE,x:xE,xi:xiE})



[ui,vi,qi,uj,vj,qj]=TE@[u3,v3,theta3E,u4,v4,theta4]

# Desplazamiento axiales 

uEh= ui*pi1E+uj*pi4E

uEf= sy.integrate(Gxx2E*pE.subs({xE:xiE}),(xiE,0,xE))+sy.integrate(Gxx1E*pE.subs({xE:xiE}),(xiE,xE,LE))


uE= uEh+uEf 

# Desplazameitno perpendiculares al elemento 

vEh= vi*pi2E+qi*pi3E+vj*pi5E+qj*pi6E
vEf= sy.integrate(Gyy2E*qE.subs({xE:xiE}),(xiE,0,xE))+sy.integrate(Gyy1E*qE.subs({xE:xiE}),(xiE,xE,LE))
vE= vEh+vEf

# Campo Axial 

PE= AEE*sy.diff(uE,xE,1)

#Momentos del elemento E

ME= EIC*sy.diff(vE,xE,2)

# Cortante del elemento E

VE=-EIC*sy.diff(vE,xE,3)



# Sumatoria de fuerzas 


Sx= FX1+FX5-sy.integrate(fA,(xA,0,LA))-sy.integrate(fB,(xB,0,LB))-30

Sy=FY1+FY5+sy.integrate(fC,(xC,0,LC))-60

M0= sy.integrate(fA*xA,(xA,0,LA))+sy.integrate(fB*xB,(xB,0,LB))-60*(LD*sy.cos(ThetaD)-2/3*2)+(2/3*2+LB)*30+M1+M5+FY5*LC + sy.integrate(fC*xC,(xC,0,LC))



# DIAGRAMAS DE DESPLAZAMIENTOS 

EscDiaDes=1000
EscDiaAxi=0.01
EscDiaCor=0.02
EscDiaMom=0.02

Nx=100

xAN=np.linspace(0,LA,Nx)
xBN=np.linspace(0,LB,Nx)
xCN=np.linspace(0,LC,Nx)
xDN=np.linspace(0,LD,Nx)
xEN=np.linspace(0,LE,Nx)

uANum=np.zeros([Nx])
uBNum=np.zeros([Nx])
uCNum=np.zeros([Nx])
uDNum=np.zeros([Nx])
uENum=np.zeros([Nx])

vANum=np.zeros([Nx])
vBNum=np.zeros([Nx])
vCNum=np.zeros([Nx])
vDNum=np.zeros([Nx])
vENum=np.zeros([Nx])

PANum=np.zeros([Nx])
PBNum=np.zeros([Nx])
PCNum=np.zeros([Nx])
PDNum=np.zeros([Nx])
PENum=np.zeros([Nx])

MANum=np.zeros([Nx])
MBNum=np.zeros([Nx])
MCNum=np.zeros([Nx])
MDNum=np.zeros([Nx])
MENum=np.zeros([Nx])

VANum=np.zeros([Nx])
VBNum=np.zeros([Nx])
VCNum=np.zeros([Nx])
VDNum=np.zeros([Nx])
VENum=np.zeros([Nx])

for ix in range(Nx):
    
    uANum[ix]=uA.subs(xA,xAN[ix])
    uBNum[ix]=uB.subs(xB,xBN[ix])
    uCNum[ix]=uC.subs(xC,xCN[ix])
    uENum[ix]=uE.subs(xE,xEN[ix])
    
    vANum[ix]=vA.subs(xA,xAN[ix])
    vBNum[ix]=vB.subs(xB,xBN[ix])
    vCNum[ix]=vC.subs(xC,xCN[ix])
    vENum[ix]=vE.subs(xE,xEN[ix])

    PANum[ix]=PA.subs(xA,xAN[ix])
    PBNum[ix]=PB.subs(xB,xBN[ix])
    PCNum[ix]=PC.subs(xC,xCN[ix])
    PENum[ix]=PE.subs(xE,xEN[ix])

    MANum[ix]=MA.subs(xA,xAN[ix])
    MBNum[ix]=MB.subs(xB,xBN[ix])
    MCNum[ix]=MC.subs(xC,xCN[ix])
    MENum[ix]=ME.subs(xE,xEN[ix])

    VANum[ix]=VA.subs(xA,xAN[ix])
    VBNum[ix]=VB.subs(xB,xBN[ix])
    VCNum[ix]=VC.subs(xC,xCN[ix])
    VENum[ix]=VE.subs(xE,xEN[ix])
    

    if xDN[ix]<LD/2:
        vDNum[ix]=vD1.subs(xD,xDN[ix])
        uDNum[ix]=uD1.subs(xD,xDN[ix])
        MDNum[ix]=MD1.subs(xD,xDN[ix])
        PDNum[ix]=PD1.subs(xD,xDN[ix])
        VDNum[ix]=VD1.subs(xD,xDN[ix])
        
    else:
        vDNum[ix]=vD2.subs(xD,xDN[ix])
        uDNum[ix]=uD2.subs(xD,xDN[ix])
        MDNum[ix]=MD2.subs(xD,xDN[ix])
        PDNum[ix]=PD2.subs(xD,xDN[ix])
        VDNum[ix]=VD2.subs(xD,xDN[ix])


plt.figure(1,figsize=(8,4))

plt.plot([-2,LC+2],[0,0],color='g',linewidth=1)
plt.plot([-2,LC+2],[LA,LA],color='g',linewidth=1)

plt.plot([0,0],[0,LA],color='k',linewidth=1)
plt.plot([0,4],[LA,7],color='k',linewidth=1)
plt.plot([4,LC],[7,LB],color='k',linewidth=1)
plt.plot([LC,LC],[0,LB],color='k',linewidth=1)
plt.plot([0,LC],[LB,LB],color='k',linewidth=1)

#   Elemento A
xGloA=0+(xAN+uANum*EscDiaDes)*np.cos(ThetaA)-vANum*np.sin(ThetaA)*EscDiaDes
yGloA=0+(xAN+uANum*EscDiaDes)*np.sin(ThetaA)+vANum*np.cos(ThetaA)*EscDiaDes    

#   Elemento B

xGloB=LC+(xBN+uBNum*EscDiaDes)*np.cos(ThetaB)-vBNum*np.sin(ThetaB)*EscDiaDes
yGloB=0+(xBN+uBNum*EscDiaDes)*np.sin(ThetaB)+vBNum*np.cos(ThetaB)*EscDiaDes 

#   Elemento C

xGloC=(xCN+uCNum*EscDiaDes)
yGloC=LB+vCNum*EscDiaDes

#   Elemento D
# 
xGloD=(xDN+uDNum*EscDiaDes)*np.cos(ThetaD)-vDNum*np.sin(ThetaD)*EscDiaDes
yGloD=LB+(xDN+uDNum*EscDiaDes)*np.sin(ThetaD)+vDNum*np.cos(ThetaD)*EscDiaDes  
# 
#   Elemento E

xGloE=4+(xEN+uENum*EscDiaDes)*np.cos(ThetaE)-vENum*np.sin(ThetaE)*EscDiaDes
yGloE=LB+2+(xEN+uENum*EscDiaDes)*np.sin(ThetaE)+vENum*np.cos(ThetaE)*EscDiaDes


plt.plot(xGloA,yGloA,color='k',linestyle='--')
plt.plot(xGloB,yGloB,color='k',linestyle='--')
plt.plot(xGloC,yGloC,color='k',linestyle='--')
plt.plot(xGloD,yGloD,color='k',linestyle='--')
plt.plot(xGloE,yGloE,color='k',linestyle='--')

plt.text(0.1,0.1,'(0,0)')
plt.text(LC+0.1,0.1,'(0,0)')
plt.text(-2,LA+0.1,'('+str('%.3e' % u2)+','+str('%.3e' % v2)+')')
plt.text(3,LA+2.2,'('+str('%.3e' % u3)+','+str('%.3e' % v3)+')')
plt.text(LC+0.1,LA+0.1,'('+str('%.3e' % u4)+','+str('%.3e' % v4)+')')

#plt.text(LB+0.2,LB+2+2,'('+str('%.3e' % u6)+','+str('%.3e' % v6)+')')

plt.grid('on')
plt.xlabel(r'$x$ [m]',fontsize=16)
plt.ylabel(r'$y$ [m]',fontsize=16)
plt.tick_params(labelsize=16)
plt.axis('equal')

#%%   8.2 Figura del campo de fuerza axial

plt.figure(2,figsize=(8,4))

plt.plot([-2,LC+2],[0,0],color='g',linewidth=1)
plt.plot([-2,LC+2],[LA,LA],color='g',linewidth=1)

plt.plot([0,0],[0,LA],color='k',linewidth=1)
plt.plot([0,4],[LA,7],color='k',linewidth=1)
plt.plot([4,LC],[7,LB],color='k',linewidth=1)
plt.plot([LC,LC],[0,LB],color='k',linewidth=1)
plt.plot([0,LC],[LB,LB],color='k',linewidth=1)

#   Elemento A
xGloA=0+xAN*np.cos(ThetaA)-PANum*np.sin(ThetaA)*EscDiaAxi
yGloA=0+xAN*np.sin(ThetaA)+PANum*np.cos(ThetaA)*EscDiaAxi    

#   Elemento B

xGloB=LC+xBN*np.cos(ThetaB)-PBNum*np.sin(ThetaB)*EscDiaAxi
yGloB=0+xBN*np.sin(ThetaB)+PBNum*np.cos(ThetaB)*EscDiaAxi

#   Elemento C

xGloC=xCN
yGloC=LB+PCNum*EscDiaAxi

#   Elemento D

xGloD=xDN*np.cos(ThetaD)-PDNum*np.sin(ThetaD)*EscDiaAxi
yGloD=LB+xDN*np.sin(ThetaD)+PDNum*np.cos(ThetaD)*EscDiaAxi

#   Elemento E

xGloE=4+xEN*np.cos(ThetaE)-PENum*np.sin(ThetaE)*EscDiaAxi
yGloE=LB+2+xEN*np.sin(ThetaE)+PENum*np.cos(ThetaE)*EscDiaAxi


plt.plot(xGloA,yGloA,color='k',linestyle='--')
plt.plot(xGloB,yGloB,color='k',linestyle='--')
plt.plot(xGloC,yGloC,color='k',linestyle='--')
plt.plot(xGloD,yGloD,color='k',linestyle='--')
plt.plot(xGloE,yGloE,color='k',linestyle='--')

plt.text(0.1,0.1,PA.subs({xA:0}))

plt.text(LC+0.1,0.1,PB.subs({xB:0}))

plt.text(0.6,LA- 0.2,''+str ('%.3e' % PD1.subs({xD:0})))

plt.text(3,LA+2.2,'('+str('%.3e' % PC.subs({xC:LC}))+','+str('%.3e' % PE.subs({xE:LE}))+')')

plt.text(LC+0.1,LA+0.1,'('+str('%.3e' % PD2.subs({xD:LD}))+','+str('%.3e' % PE.subs({xE:0}))+')')

#plt.text(LB+0.2,LB+2+2,'('+str('%.3e' % u6)+','+str('%.3e' % v6)+')')

plt.grid('on')
plt.xlabel(r'$x$ [KN]',fontsize=16)
plt.ylabel(r'$y$ [KN]',fontsize=16)
plt.tick_params(labelsize=16)
plt.axis('equal')
#%%   8.2 Figura del campo de fuerza CORTANTES


plt.figure(3,figsize=(8,4))
plt.plot([-2,LC+2],[0,0],color='g',linewidth=1)
plt.plot([-2,LC+2],[LA,LA],color='g',linewidth=1)
plt.plot([0,0],[0,LA],color='k',linewidth=1)
plt.plot([0,4],[LA,7],color='k',linewidth=1)
plt.plot([4,LC],[7,LB],color='k',linewidth=1)
plt.plot([LC,LC],[0,LB],color='k',linewidth=1)
plt.plot([0,LC],[LB,LB],color='k',linewidth=1)

#   Elemento A
xGloA=0+xAN*np.cos(ThetaA)-VANum*np.sin(ThetaA)*EscDiaCor
yGloA=0+xAN*np.sin(ThetaA)+VANum*np.cos(ThetaA)*EscDiaCor 

#   Elemento B

xGloB=LC+xBN*np.cos(ThetaB)-VBNum*np.sin(ThetaB)*EscDiaCor
yGloB=0+xBN*np.sin(ThetaB)+VBNum*np.cos(ThetaB)*EscDiaCor

#   Elemento C

xGloC=xCN
yGloC=LB+VCNum*EscDiaCor

#   Elemento D

xGloD=xDN*np.cos(ThetaD)-VDNum*np.sin(ThetaD)*EscDiaCor
yGloD=LB+xDN*np.sin(ThetaD)+VDNum*np.cos(ThetaD)*EscDiaCor

#   Elemento E

xGloE=4+xEN*np.cos(ThetaE)-VENum*np.sin(ThetaE)*EscDiaCor
yGloE=LB+2+xEN*np.sin(ThetaE)+VENum*np.cos(ThetaE)*EscDiaCor


plt.plot(xGloA,yGloA,color='k',linestyle='--')
plt.plot(xGloB,yGloB,color='k',linestyle='--')
plt.plot(xGloC,yGloC,color='k',linestyle='--')
plt.plot(xGloD,yGloD,color='k',linestyle='--')
plt.plot(xGloE,yGloE,color='k',linestyle='--')



plt.text(0.1,0.1,VA.subs({xA:0}))

plt.text(LC+0.1,0.1,VB.subs({xB:0}))

plt.text(0.6,LA- 0.2,''+str ('%.3e' % VD1.subs({xD:0})))

plt.text(3,LA+2.2,'('+str('%.3e' % VC.subs({xC:LC}))+','+str('%.3e' % VE.subs({xE:LE}))+')')

plt.text(LC+0.1,LA+0.1,'('+str('%.3e' % VD2.subs({xD:LD}))+','+str('%.3e' % VE.subs({xE:0}))+')')
#plt.text(LB+0.2,LB+2+2,'('+str('%.3e' % u6)+','+str('%.3e' % v6)+')')

plt.grid('on')
plt.xlabel(r'$x$ [KN]',fontsize=16)
plt.ylabel(r'$y$ [KN]',fontsize=16)
plt.tick_params(labelsize=16)
plt.axis('equal')

#%%   8.2 Figura del campo de fuerza MOMENTOS 

plt.figure(4,figsize=(8,4))
plt.plot([-2,LC+2],[0,0],color='g',linewidth=1)
plt.plot([-2,LC+2],[LA,LA],color='g',linewidth=1)
plt.plot([0,0],[0,LA],color='k',linewidth=1)
plt.plot([0,4],[LA,7],color='k',linewidth=1)
plt.plot([4,LC],[7,LB],color='k',linewidth=1)
plt.plot([LC,LC],[0,LB],color='k',linewidth=1)
plt.plot([0,LC],[LB,LB],color='k',linewidth=1)


#Elemento A
xGloA=0+xAN*np.cos(ThetaA)-MANum*np.sin(ThetaA)*EscDiaMom
yGloA=0+xAN*np.sin(ThetaA)+MANum*np.cos(ThetaA)*EscDiaMom   

#Elemento B

xGloB=LC+xBN*np.cos(ThetaB)-MBNum*np.sin(ThetaB)*EscDiaMom
yGloB=xBN*np.sin(ThetaB)+MBNum*np.cos(ThetaB)*EscDiaMom

#Elemento C

xGloC=xCN
yGloC=LB+MCNum*EscDiaMom

#Elemento D

xGloD=xDN*np.cos(ThetaD)-MDNum*np.sin(ThetaD)*EscDiaMom
yGloD=LB+xDN*np.sin(ThetaD)+MDNum*np.cos(ThetaD)*EscDiaMom

#Elemento E


xGloE=4+xEN*np.cos(ThetaE)-MENum*np.sin(ThetaE)*EscDiaMom
yGloE=LB+2+xEN*np.sin(ThetaE)+MENum*np.cos(ThetaE)*EscDiaMom


plt.plot(xGloA,yGloA,color='k',linestyle='--')
plt.plot(xGloB,yGloB,color='k',linestyle='--')
plt.plot(xGloC,yGloC,color='k',linestyle='--')
plt.plot(xGloD,yGloD,color='k',linestyle='--')
plt.plot(xGloE,yGloE,color='k',linestyle='--')


plt.text(0.1,0.1,MA.subs({xA:0}))

plt.text(LC+0.1,0.1,MB.subs({xB:0}))

plt.text(0.6,LA- 0.2,''+str ('%.3e' % MD1.subs({xD:0})))

plt.text(3,LA+2.2,'('+str('%.3e' % MC.subs({xC:LC}))+','+str('%.3e' % ME.subs({xE:LE}))+')')

plt.text(LC+0.1,LA+0.1,'('+str('%.3e' % MD2.subs({xD:LD}))+','+str('%.3e' % ME.subs({xE:0}))+')')

#plt.text(LB+0.2,LB+2+2,'('+str('%.3e' % u6)+','+str('%.3e' % v6)+')')

plt.grid('on')
plt.axis('equal')
plt.xlabel(r'$x$ [KN-m]',fontsize=16)
plt.ylabel(r'$y$ [KN-m]',fontsize=16)
plt.tick_params(labelsize=16)

#%%


# Fact=20000
# num= 11 

# xN= np.linspace(0,LE,num)


# vDes=np.zeros((num,2))

# N=np.array([[0,0],
#             [0,5],
#             [4,7],
#             [9,5],
#             [9,0]]) 



# LN=np.array([[0,1],[1,2],[]])



#for i in range (num):
    
#     xGlo= uE.subs({xE:xN[i]})*Fact+xN[i]
#     yGlo= vE.subs({xE:xN[i]})*Fact
    
#     vDes[i,0],vDes[i,1]= np.transpose(TE[:2,:2])@[[xGlo],[yGlo]]+[[4],[7]]



# plt.plot(vDes[:,0],vDes[:,1])





# #   Elemento D

# xGloD=(xDN+uDNum*EscDiaDes)*np.cos(ThetaD)-vDNum*np.sin(ThetaD)*EscDiaDes

# yGloD=LA+(xDN+uDNum*EscDiaDes)*np.sin(ThetaD)+vDNum*np.cos(ThetaD)*EscDiaDes


# plt.plot(xGloD,yGloD,color='k',linestyle='--')

# # Elemento A


# xGloA=0+(xAN+uANum*EscDiaDes)*np.cos(ThetaA)-vANum*np.sin(ThetaA)*EscDiaDes
# yGloA=0+(xAN+uANum*EscDiaDes)*np.sin(ThetaA)+vANum*np.cos(ThetaA)*EscDiaDes 


# # Elemento C 

# xGloC= LB+(xCN+uCNum*EscDiaDes)-vCNum*EscDiaDes
# yGloC= (xCN+uCNum*EscDiaDes)+vCNum*EscDiaDes


# plt.plot(xGloD,yGloD,color='k',linestyle='--')
# plt.plot(xGloA,yGloA,color='k',linestyle='--')
# plt.plot(xGloC,yGloC,color='k',linestyle='--')


# plt.text(0,LA+0.1*LA-1,'%4.6g' % VD1.subs(xD,0))
# plt.text(4,7,'%4.6g' % VD1.subs(xD,LD))

# plt.text(0.2,0.1,'(0,0)')
# plt.text(5.2,0.1,'(0,0)')
# plt.text(-3.4,4+0.1,'('+str('%.3e' % u3)+','+str('%.3e' % v3)+')')
# plt.text(5-3.6,4+0.1,'('+str('%.3e' % u4)+','+str('%.3e' % v4)+')')
# plt.text(-3.4,4+2+1,'('+str('%.3e' % u5)+','+str('%.3e' % v5)+')')


# plt.grid('on')
# plt.axis('equal')
# plt.xlabel(r'$x$ [m]',fontsize=16)
# plt.ylabel(r'$y$ [m]',fontsize=16)
# plt.tick_params(labelsize=16)




