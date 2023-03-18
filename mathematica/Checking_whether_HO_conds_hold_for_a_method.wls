#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
(*PacletInstall["http://github.com/jkuczm/MathematicaCellsToTeX/releases/download/v0.2.2/CellsToTeX-0.2.2"]*)

(*This package turns Mathematica input into Tex code*)(* Needs["CellsToTeX`"]*)
(*NCAlgebra allows for non commutative expr*)
PacletInstall["https://github.com/NCAlgebra/NC/blob/v6.0.0/NCAlgebra-6.0.0.paclet?raw=true"];
<<NCAlgebra`;

(* ::Input::Initialization:: *)
(* To clear all variables *)
Clear["Global`*"]

(*Modify the number of stages*)
stages = 3
s=stages;
(*Building the Butcher Tableau and vector of c_i*)
a=Table[Symbol["a"<>ToString[i]<>ToString[j]],{i,0,stages-1},{j,0,i}]
c=Table[Symbol["c"<>ToString[i]],{i,0,stages-1}]
(*Evaluating each phi function at each c_i*)
phi1c=Table[phi1[-c[[i]]*L*h],{i,1,stages}];
phi2c=Table[phi2[-c[[i]]*L*h],{i,1,stages}];
phi3c=Table[phi3[-c[[i]]*L*h],{i,1,stages}];
phi4c=Table[phi4[-c[[i]]*L*h],{i,1,stages}];

(*Stiff conditions as in 05 Hochbruck paper*)SetNonCommutative[J,K,L,W,a]
order1={Sum[a[[s]][[i]],{i,stages}]==phi1[-L*h]};
order2a={Sum[a[[s]][[i]]*c[[i]],{i,stages}]==phi2[-L*h]};
order2b=Table[Sum[a[[i]][[k]],{k,1,i}]==phi1c[[i+1]]*c[[i+1]],{i,1,stages-1}];
order3a={Sum[a[[s]][[i]]*c[[i]]^2/(2!),{i,1,stages}]==phi3[-L*h]};
order3b={Sum[a[[s]][[i]]**J**(phi2c[[i]]*c[[i]]^2-Sum[a[[i-1]][[j]]*c[[j]],{j,2,i-1}]),{i,2,stages}]==0};
order4a={Sum[a[[s]][[i]]*c[[i]]^3/(3!),{i,stages}]==phi4[-L*h]};
order4b={Sum[a[[s]][[i]]**J**(phi3c[[i]]*c[[i]]^3-Sum[a[[i-1]][[j]]*c[[j]]^2/(2!),{j,2,i-1}]),{i,2,stages}]==0};
order4c={Sum[a[[s]][[i]]**J**Sum[a[[i-1]][[j]]**J**(phi2c[[j]]*c[[j]]^2-Sum[a[[j-1]][[k]]*c[[k]],{k,2,j-1}]),{j,2,i-1}],{i,2,stages}]==0};
order4d={Sum[a[[s]][[i]]*c[[i]]**K**(phi2c[[i]]*c[[i]]^2-Sum[a[[i-1]][[j]]*c[[j]],{j,2,i-1}]),{i,2,stages}]==0};

(* Stiff order conditions from Luan14 paper*)
order5a = {Sum[a[[s]][[i]]*c[[i]]^4/(4!),{i,stages}]==phi5[-L*h]};
order5b={Sum[a[[s]][[i]]**J**(phi4c[[i]]*c[[i]]^4-Sum[a[[i-1]][[j]]*c[[j]]^3/(3!),{j,2,i-1}]),{i,2,stages}] == 0};
order5c = {Sum[a[[s]][[i]]**J**Sum[a[[i-1]][[j]]**J**(phi3c[[j]]*c[[j]]^3-Sum[a[[j-1]][[k]]*c[[k]]^2/(2!),{k,2,j-1}]),{j,2,i-1}],{i,2,stages}] == 0};
order5d = {Sum[a[[s]][[i]]**J**Sum[a[[i-1]][[j]]**J**Sum[a[[j-1]][[k]]**J**(phi2c[[k]]*c[[k]]^2-Sum[a[[k-1]][[l]]*c[[l]],{l,2,k-1}]),{k,2,j-1}],{j,2,i-1}],{i,2,stages}] == 0};
order5e = {Sum[a[[s]][[i]]**J**Sum[(a[[i-1]][[j]]*c[[j]])**K**((phi2c[[j]]*c[[j]]^2 )-Sum[a[[j-1]][[k]]*c[[k]],{k,2,j-1}]),{j,2,i-1}],{i,2,stages}]== 0};
order5f={Sum[a[[s]][[i]]*c[[i]]**K**(phi3c[[i]]*c[[i]]^3-Sum[a[[i-1]][[j]]*c[[j]]^2/(2!),{j,2,i-1}]),{i,2,stages}] == 0};
order5g = {Sum[a[[s]][[i]]*c[[i]]**K**Sum[a[[i-1]][[j]]**J**(phi2c[[j]]*c[[j]]^2-Sum[a[[j-1]][[k]]*c[[k]],{k,2,j-1}]),{j,2,i-1}],{i,2,stages}] == 0};

(* Condition 17 is satisfied like it is satisfied in the paper by Lu&Os14 *)
order5h=Table[a[[s]][[i]]==0 ||phi2c[[i]]*c[[i]]^2-Sum[a[[i-1]][[j]]*c[[j]],{j,2,i-1}]==0,{i,2,stages}] ;

order5i={Sum[a[[s]][[i]]*c[[i]]^2**W**(phi2c[[i]]*c[[i]]^2-Sum[a[[i-1]][[j]]*c[[j]],{j,2,i-1}]),{i,2,stages}] == 0};


(* ::Input:: *)
(*phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);phi5[0]=1/(5!);*)
(**)
(*allorder3=Join[order3a,order3b];*)
(*allorder3/.{L->0}*)
(*ExpandNonCommutativeMultiply[allorder3/.{L->0}]*)


(* ::Input:: *)
(*c1 a32*a11+c1 a33*a21+c2 a33*a22 -1/6==0*)


(* ::Input:: *)
(**)
(*ExpandNonCommutativeMultiply[allorder3]*)


(* ::Input::Initialization:: *)
(* Testing the second order method from ERKBS32 *)
c0=0;c1=1/2;c2=3/4;c3=1;
a00=c1*phi1[-c1*L*h];
a11=9/8 *phi2[-c2*L*h] + 3/8 *phi2[-c1*L*h];
a10=c2*phi1[-c2*L*h]-a11;
a21=1/3 * phi1[-L*h];
a22=4/3 * phi2[-L*h]-2/9 * phi1[-L*h];
a20=phi1[-L*h]-a21-a22;

a30=phi1[-L*h]-17/12*phi2[-L*h];
a31=1/2 * phi2[-L*h];
a32=2/3 * phi2[-L*h];
a33=1/4 * phi2[-L*h];

(*  Checking whether conditions hold exactly for orders 1 and 2 *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]

phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);phi5[0]=1/(5!);

(* Checking that this method is never 3rd order *)
allorder3=Join[order3a,order3b];
Print["21b"]
order3b
FullSimplify[allorder3]
allorder3/.{L->0}
Print["###################"]
a30=a30/.L->0;
a31=a31/.L->0;
a32=a32/.L->0;
a33=a33/.L->0;
FullSimplify[order3b]


(* ::Input::Initialization:: *)
(* Testing the third order method from ERKBS32 *)
c0=0;c1=1/2;c2=3/4;c3=1;

a00=c1*phi1[-c1*L*h];
a11=9/8 *phi2[-c2*L*h] + 3/8 *phi2[-c1*L*h];
a10=c2*phi1[-c2*L*h]-a11;
a21=1/3 * phi1[-L*h];
a22=4/3 * phi2[-L*h]-2/9 * phi1[-L*h];
a20=phi1[-L*h]-a21-a22;

(*  Checking whether conditions hold exactly for order 1 *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]

(*   Checking whether condition 3 holds weakly i.e. b_i(0) instead       of b_i(-hL)  *)
(*  Always reread Theorem 4.7 in the 05 paper when applying this step, "holding weakly" has a very particular meaning *) 

a20=a20/.L->0;
a21=a21/.L->0;
a22=a22/.L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);phi5[0]=1/(5!);

(* The following line takes care of psi_p(0)=0 in "holding weakly" *)
order3a = order3a/.{L->0};
allorder3 = Join[order3a,order3b];

FullSimplify[allorder3]


(* ::Input::Initialization:: *)
(* Testing the second order method from ERK32BP *)
c0=0;c1=1/2;c2=3/4;c3=1;
a00=c1*phi1[-c1*L*h];
a11=3/2*phi2[-c2*L*h];
a10=c2*phi1[-c2*L*h]-a11;
a21=2 * phi2[-L*h]-4*phi3[-L*h];
a22=8/3 * phi3[-L*h];
a20=phi1[-L*h]-2 * phi2[-L*h]+4/3 * phi3[-L*h];

a30=phi1[-L*h]-17/12*phi2[-L*h];
a31=1/2 * phi2[-L*h];
a32=2/3 * phi2[-L*h];
a33=1/4 * phi2[-L*h];

(*  Checking whether conditions hold exactly for orders 1 and 2 *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]

phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);phi5[0]=1/(5!);

(* Checking that this method is never 3rd order *)
allorder3=Join[order3a,order3b];
allorder3/.{L->0}


(* ::Input::Initialization:: *)
(* Testing the third order method from ERK32BP *)
c0=0;c1=1/2;c2=3/4;c3=1;

a00=c1*phi1[-c1*L*h];
a11=3/2*phi2[-c2*L*h];
a10=c2*phi1[-c2*L*h]-a11;
a21=2 * phi2[-L*h]-4*phi3[-L*h];
a22=8/3 * phi3[-L*h];
a20=phi1[-L*h]-2 * phi2[-L*h]+4/3 * phi3[-L*h];

(*  Checking whether conditions hold exactly for order 1 *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]

(*   Checking whether condition 3 holds weakly i.e. b_i(0) instead       of b_i(-hL)  *)
(*  Always reread Theorem 4.7 in the 05 paper when applying this step, "holding weakly" has a very particular meaning *) 

a20=a20/.L->0;
a21=a21/.L->0;
a22=a22/.L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);phi5[0]=1/(5!);

(* The following line takes care of psi_p(0)=0 in "holding weakly" *)
order3a = order3a/.{L->0};
allorder3 = Join[order3a,order3b];

FullSimplify[allorder3]


(* ::Input::Initialization:: *)
(* Testing the second order method from ERK32robust *)
c0=0;c1=1/2;c2=3/4;c3=1;

a00=1/2 phi1[-((h L)/2)];a10=3/4 phi1[-((3 h L)/4)]-9/8 phi2[-((3 h L)/4)]-3/8 phi2[-((h L)/2)];
a11=9/8 phi2[-((3 h L)/4)]+3/8 phi2[-((h L)/2)];
a20=phi1[-h L]-19/12 phi2[-h L]+1/12 phi3[-h L];
a21=3/4 phi2[-h L]-1/4 phi3[-h L];
a22=5/6 phi2[-h L]+1/6 phi3[-h L];
a30=29/18 phi1[-h L]+7/6 phi1[-((3 h L)/4)]+9/14 phi1[-((h L)/2)]+3/4 phi2[-h L]+2/7 phi2[-((3 h L)/4)]+1/12 phi2[-((h L)/2)]-8083/420 phi3[-h L]+11/30 phi3[-((h L)/2)];
a31=-(1/9) phi1[-h L]-1/6 phi1[-((3 h L)/4)]-1/2 phi2[-h L]-1/7 phi2[-((3 h L)/4)]-1/3 phi2[-((h L)/2)]+1/6 phi3[-h L]+1/6 phi3[-((h L)/2)];
a32=2/3 phi1[-h L]-1/2 phi1[-((3 h L)/4)]-1/7 phi1[-((h L)/2)]+1/3 phi2[-h L]-1/7 phi2[-((3 h L)/4)]-1/5 phi3[-((h L)/2)];
a33=-(7/6) phi1[-h L]-1/2 phi1[-((3 h L)/4)]-1/2 phi1[-((h L)/2)]-7/12 phi2[-h L]+1/4 phi2[-((h L)/2)]+2671/140 phi3[-h L]-1/3 phi3[-((h L)/2)];

(*  Checking whether conditions hold exactly for orders 1 *)
allorder1 = order1;
FullSimplify[allorder1]

(*   Checking whether condition 3 holds weakly i.e. b_i(0) instead       of b_i(-hL)  *)
(*  Always reread Theorem 4.7 in the 05 paper when applying this step, "holding weakly" has a very particular meaning *) 

a30=a30/.L->0;
a31=a31/.L->0;
a32=a32/.L->0;
a33=a33/.L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);phi5[0]=1/(5!);

(* The following line takes care of psi_p(0)=0 in "holding weakly" *)
order2a = order2a/.{L->0};
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]

(* Checking that this method is never 3rd order *)
allorder3=Join[order3a,order3b];
allorder3/.{L->0}


(* ::Input::Initialization:: *)
(* Testing the third order method from ERK32robust*)
(* norm E4 for this estimate is 1.07026 *)
c0=0;c1=1/2;c2=3/4;c3=1;

a00=1/2 phi1[-((h L)/2)];a10=3/4 phi1[-((3 h L)/4)]-9/8 phi2[-((3 h L)/4)]-3/8 phi2[-((h L)/2)];
a11=9/8 phi2[-((3 h L)/4)]+3/8 phi2[-((h L)/2)];
a20=phi1[-h L]-19/12 phi2[-h L]+1/12 phi3[-h L];
a21=3/4 phi2[-h L]-1/4 phi3[-h L];
a22=5/6 phi2[-h L]+1/6 phi3[-h L];

(*  Checking whether conditions hold exactly for order 1 *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]

(*   Checking whether condition 3 holds weakly i.e. b_i(0) instead       of b_i(-hL)  *)
(*  Always reread Theorem 4.7 in the 05 paper when applying this step, "holding weakly" has a very particular meaning *) 

a20=a20/.L->0;
a21=a21/.L->0;
a22=a22/.L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);phi5[0]=1/(5!);

(* The following line takes care of psi_p(0)=0 in "holding weakly" *)
order3a = order3a/.{L->0};
allorder3 = Join[order3a,order3b];

FullSimplify[allorder3]


(* ::Input::Initialization:: *)
(*Testing the method ERK4HO5*)
c0=0;c1=1/2;c2=1/2;c3=1;c4=1/2;
a00=1/2*phi1[-c1*L*h];
a10=1/2*phi1[-c2*L*h] -phi2[-c2*L*h];
a11=phi2[-c2*L*h];
a20=phi1[-c3*L*h]-2*phi2[-c3*L*h];
a21=phi2[-c3*L*h];
a22=phi2[-c3*L*h];
a31=1/2*phi2[-c4*L*h]-phi3[-c3*L*h]+1/4*phi2[-c3*L*h]-1/2phi3[-c4*L*h];
a32=a31;
a33=1/4*phi2[-c4*L*h]-a31;
a30=1/2*phi1[-c4*L*h]-2*a31-a33;

a40=phi1[-L*h]-3*phi2[-L*h]+4*phi3[-L*h];
a41=0;
a42=0;
a43=-phi2[-L*h]+4*phi3[-L*h];
a44=4*phi2[-L*h]-8*phi3[-L*h];

(*  Checking whether conditions hold exactly for orders 1,2,3  *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]
allorder3 = Join[order3a,order3b];
FullSimplify[ExpandNonCommutativeMultiply[allorder3]]

(*   Checking whether condition 4 holds weakly i.e. b_i(0) instead       of b_i(-hL)  *)
(*  Always reread Theorem 4.7 in the 05 paper when applying this step, "holding weakly" has a very particular meaning *) 
a40=a40/.L->0;
a41=a41/.L->0;
a42=a42/.L->0;
a43=a43/.L->0;
a44=a44/.L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);phi5[0]=1/(5!);
(* The following line takes care of psi_p(0)=0 in "holding weakly" *)
order4a = order4a/.{L->0};
allorder4=Join[order4a,order4b,order4c,order4d];
FullSimplify[ExpandNonCommutativeMultiply[allorder4]]


(* ::Input::Initialization:: *)
(*Testing the method ERK4K*)
c0=0;c1=1/2;c2=1/2;c3=1;
a00=1/2*phi1[-c1*L*h];

a10=1/2*phi1[-c2*L*h]-phi2[-c2*L*h];
a11=phi2[-c2*L*h];

a20=phi1[-c3*L*h]-2*phi2[-c3*L*h];
a21=0;
a22=2*phi2[-c3*L*h];

a30=phi1[-L*h]-3*phi2[-L*h]+4*phi3[-L*h];
a31=2*phi2[-L*h]-4*phi3[-L*h];
a32=2*phi2[-L*h]-4*phi3[-L*h];
a33=4*phi3[-L*h]-phi2[-L*h];

(*Checking whether conditions hold exactly for orders 1,2,3*)

allorder1=order1;
FullSimplify[allorder1]
allorder2=Join[order2a,order2b];
FullSimplify[allorder2]
allorder3=Join[order3a,order3b];
FullSimplify[ExpandNonCommutativeMultiply[allorder3]]

(*Checking whether condition 4 holds weakly i.e.b_i(0) instead of b_i(-hL)*)(*Always reread Theorem 4.7 in the 05 paper when applying this step,"holding weakly" has a very particular meaning*)
a40=a40/. L->0;
a41=a41/. L->0;
a42=a42/. L->0;
a43=a43/. L->0;
a44=a44/. L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);
phi5[0]=1/(5!);
(*The following line takes care of psi_p(0)=0 in "holding weakly"*)

order4a=order4a/. {L->0};
allorder4=Join[order4a,order4b,order4c,order4d];
FullSimplify[ExpandNonCommutativeMultiply[allorder4]]


(* ::Input::Initialization:: *)
(* The third order method of ERK43ZB *)
c0=0;c1=1/6;c2=1/2;c3=1/2;c4=1;
a00=1/6*phi1[-L*h/6];

a11=3/2*phi2[-L*h/2]+1/2*phi2[-L*h/6];
a10=1/2*phi1[-L*h/2] -a11;

a22 = -19/180*phi1[-L*h]-1/6*phi1[-L*h/2]-1/6*phi1[-L*h/6]-1/6*phi2[-L*h/2]+1/9*phi2[-L*h/6]-1/5*phi3[-L*h/2];
a21 =  19/60*phi1[-L*h]+1/2*phi1[-L*h/2]+1/2*phi1[-L*h/6]+2*phi2[-L*h/2]+13/6*phi2[-L*h/6]+3/5*phi3[-L*h/2];
a20 = 1/2*phi1[-L*h/2] - a21-a22;

a33=phi2[-L*h]+phi2[-L*h/2]-6*phi3[-L*h]-3*phi3[-L*h/2];
a31=3*phi2[-L*h]-9/2*phi2[-L*h/2]-5/2*phi2[-L*h/6]+6*a33+a21;
a32=6*phi3[-L*h]+3*phi3[-L*h/2]-2*a33+a22;
a30=phi1[-L*h]-a31-a32-a33;

(*  Checking whether conditions hold exactly for orders 1,2  *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]

(*Checking whether condition 3 holds weakly i.e.b_i(0) instead of b_i(-hL)*)(*Always reread Theorem 4.7 in the 05 paper when applying this step,"holding weakly" has a very particular meaning*)a30=a30/. L->0;
a31=a31/. L->0;
a32=a32/. L->0;
a33=a33/. L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);
phi5[0]=1/(5!);

(*The following line takes care of psi_p(0)=0 in "holding weakly"*)

order3a=order3a/. {L->0};
allorder3=Join[order3a,order3b];

FullSimplify[allorder3]

(* Checking that this method is never 4th order *)
allorder4=Join[order4a,order4b,order4c,order4d];
allorder4/.{L->0}


(* ::Input:: *)
(**)
(*a/.{L->0}*)


(* ::Input::Initialization:: *)
(* The fourth order method of ERK43ZB *)
c0=0;c1=1/6;c2=1/2;c3=1/2;c4=1;c5=1;
a00=1/6*phi1[-L*h/6];

a11=3/2*phi2[-L*h/2]+1/2*phi2[-L*h/6];
a10=1/2*phi1[-L*h/2] -a11;

a22 = -19/180*phi1[-L*h]-1/6*phi1[-L*h/2]-1/6*phi1[-L*h/6]-1/6*phi2[-L*h/2]+1/9*phi2[-L*h/6]-1/5*phi3[-L*h/2];
a21 =  19/60*phi1[-L*h]+1/2*phi1[-L*h/2]+1/2*phi1[-L*h/6]+2*phi2[-L*h/2]+13/6*phi2[-L*h/6]+3/5*phi3[-L*h/2];
a20 = 1/2*phi1[-L*h/2] - a21-a22;

a33=phi2[-L*h]+phi2[-L*h/2]-6*phi3[-L*h]-3*phi3[-L*h/2];
a31=3*phi2[-L*h]-9/2*phi2[-L*h/2]-5/2*phi2[-L*h/6]+6*a33+a21;
a32=6*phi3[-L*h]+3*phi3[-L*h/2]-2*a33+a22;
a30=phi1[-L*h]-a31-a32-a33;

a44=4/3*phi3[-L*h]-1/9*phi2[-L*h];
a43=7/9*phi2[-L*h]-10/3*phi3[-L*h];
a42=26/3*phi3[-L*h]-11/9*phi2[-L*h];
a41=8*phi2[-L*h]-24*phi3[-L*h];
a40=phi1[-L*h]-67/9*phi2[-L*h]+52/3*phi3[-L*h];

(*  Checking whether conditions hold exactly for orders 1,2,3  *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]
allorder3 = Join[order3a,order3b];
FullSimplify[ExpandNonCommutativeMultiply[allorder3]]

(*   Checking whether condition 4 holds weakly i.e. b_i(0) instead       of b_i(-hL)  *)
(*  Always reread Theorem 4.7 in the 05 paper when applying this step, "holding weakly" has a very particular meaning *) 
a40=a40/.L->0;
a41=a41/.L->0;
a42=a42/.L->0;
a43=a43/.L->0;
a44=a44/.L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);phi5[0]=1/(5!);
(* The following line takes care of psi_p(0)=0 in "holding weakly" *)
order4a = order4a/.{L->0};
allorder4=Join[order4a,order4b,order4c,order4d];
FullSimplify[ExpandNonCommutativeMultiply[allorder4]]


(* ::Input:: *)
(**)
(*a/.{L->0}*)


(* ::Input::Initialization:: *)
(* Checking the third order method of ERK43ZBrobust *)
c0=0;c1=1/6;c2=1/2;c3=1/2;c4=1;c5=1;

a00=1/6 phi1[-((h L)/6)];

a10=1/2 phi1[-((h L)/2)]-3/2 phi2[-((h L)/2)]-1/2 phi2[-((h L)/6)];
a11=3/2 phi2[-((h L)/2)]+1/2 phi2[-((h L)/6)];

a20=-(19/90) phi1[-h L]+1/6 phi1[-((h L)/2)]-1/3 phi1[-((h L)/6)]-11/6 phi2[-((h L)/2)]-41/18 phi2[-((h L)/6)]-2/5 phi3[-((h L)/2)];
a21=19/60 phi1[-h L]+1/2 phi1[-((h L)/2)]+1/2 phi1[-((h L)/6)]+2 phi2[-((h L)/2)]+13/6 phi2[-((h L)/6)]+3/5 phi3[-((h L)/2)];
a22=-(19/180) phi1[-h L]-1/6 phi1[-((h L)/2)]-1/6 phi1[-((h L)/6)]-1/6 phi2[-((h L)/2)]+1/9 phi2[-((h L)/6)]-1/5 phi3[-((h L)/2)];

a30=71/90 phi1[-h L]-1/3 phi1[-((h L)/2)]-1/3 phi1[-((h L)/6)]-8 phi2[-h L]-7/3 phi2[-((h L)/2)]+2/9 phi2[-((h L)/6)]+24 phi3[-h L]+58/5 phi3[-((h L)/2)];
a31=19/60 phi1[-h L]+1/2 phi1[-((h L)/2)]+1/2 phi1[-((h L)/6)]+9 phi2[-h L]+7/2 phi2[-((h L)/2)]-1/3 phi2[-((h L)/6)]-36 phi3[-h L]-87/5 phi3[-((h L)/2)];
a32=-(19/180) phi1[-h L]-1/6 phi1[-((h L)/2)]-1/6 phi1[-((h L)/6)]-2 phi2[-h L]-13/6 phi2[-((h L)/2)]+1/9 phi2[-((h L)/6)]+18 phi3[-h L]+44/5 phi3[-((h L)/2)];
a33=phi2[-h L]+phi2[-((h L)/2)]-6 phi3[-h L]-3 phi3[-((h L)/2)];


(*  Checking whether conditions hold exactly for orders 1,2  *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]

(*Checking whether condition 3 holds weakly i.e.b_i(0) instead of b_i(-hL)*)(*Always reread Theorem 4.7 in the 05 paper when applying this step,"holding weakly" has a very particular meaning*)a30=a30/. L->0;
a31=a31/. L->0;
a32=a32/. L->0;
a33=a33/. L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);
phi5[0]=1/(5!);

(*The following line takes care of psi_p(0)=0 in "holding weakly"*)

order3a=order3a/. {L->0};
allorder3=Join[order3a,order3b];

FullSimplify[allorder3]

(* Checking that this method is never 4th order *)
order4a=order4a/. {L->0};
allorder4=Join[order4a,order4b,order4c,order4d];

allorder4/.{L->0}


(* ::Input::Initialization:: *)
(* Checking the fourth order method of ERK43ZBrobust *)
c0=0;c1=1/6;c2=1/2;c3=1/2;c4=1;c5=1;

a00=1/6 phi1[-((h L)/6)];

a10=1/2 phi1[-((h L)/2)]-3/2 phi2[-((h L)/2)]-1/2 phi2[-((h L)/6)];
a11=3/2 phi2[-((h L)/2)]+1/2 phi2[-((h L)/6)];

a20=-(19/90) phi1[-h L]+1/6 phi1[-((h L)/2)]-1/3 phi1[-((h L)/6)]-11/6 phi2[-((h L)/2)]-41/18 phi2[-((h L)/6)]-2/5 phi3[-((h L)/2)];
a21=19/60 phi1[-h L]+1/2 phi1[-((h L)/2)]+1/2 phi1[-((h L)/6)]+2 phi2[-((h L)/2)]+13/6 phi2[-((h L)/6)]+3/5 phi3[-((h L)/2)];
a22=-(19/180) phi1[-h L]-1/6 phi1[-((h L)/2)]-1/6 phi1[-((h L)/6)]-1/6 phi2[-((h L)/2)]+1/9 phi2[-((h L)/6)]-1/5 phi3[-((h L)/2)];

a30=71/90 phi1[-h L]-1/3 phi1[-((h L)/2)]-1/3 phi1[-((h L)/6)]-8 phi2[-h L]-7/3 phi2[-((h L)/2)]+2/9 phi2[-((h L)/6)]+24 phi3[-h L]+58/5 phi3[-((h L)/2)];
a31=19/60 phi1[-h L]+1/2 phi1[-((h L)/2)]+1/2 phi1[-((h L)/6)]+9 phi2[-h L]+7/2 phi2[-((h L)/2)]-1/3 phi2[-((h L)/6)]-36 phi3[-h L]-87/5 phi3[-((h L)/2)];
a32=-(19/180) phi1[-h L]-1/6 phi1[-((h L)/2)]-1/6 phi1[-((h L)/6)]-2 phi2[-h L]-13/6 phi2[-((h L)/2)]+1/9 phi2[-((h L)/6)]+18 phi3[-h L]+44/5 phi3[-((h L)/2)];
a33=phi2[-h L]+phi2[-((h L)/2)]-6 phi3[-h L]-3 phi3[-((h L)/2)];

a40=phi1[-h L]-67/9 phi2[-h L]+52/3 phi3[-h L];
a41=8 phi2[-h L]-24 phi3[-h L];
a42=-(11/9) phi2[-h L]+26/3 phi3[-h L];
a43=7/9 phi2[-h L]-10/3 phi3[-h L];
a44=-(1/9) phi2[-h L]+4/3 phi3[-h L];

(*  Checking whether conditions hold exactly for orders 1,2,3  *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]
allorder3 = Join[order3a,order3b];
FullSimplify[ExpandNonCommutativeMultiply[allorder3]]

(*   Checking whether condition 4 holds weakly i.e. b_i(0) instead       of b_i(-hL)  *)
(*  Always reread Theorem 4.7 in the 05 paper when applying this step, "holding weakly" has a very particular meaning *) 
a40=a40/.L->0;
a41=a41/.L->0;
a42=a42/.L->0;
a43=a43/.L->0;
a44=a44/.L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);phi5[0]=1/(5!);
(* The following line takes care of psi_p(0)=0 in "holding weakly" *)
order4a = order4a/.{L->0};
allorder4=Join[order4a,order4b,order4c,order4d];
FullSimplify[ExpandNonCommutativeMultiply[allorder4]]


(* ::Input::Initialization:: *)
(* The third order method of ERK43ZBa *)
c0=0;c1=1/2;c2=1/2;c3=1;c4=1/2;c5=1;

a00=1/2 phi1[-((h L)/2)];

a10=1/2 phi1[-((h L)/2)]-2/3 phi2[-((h L)/2)];
a11=2/3 phi2[-((h L)/2)];

a20=phi1[-h L]-2 phi2[-h L];
a21=phi2[-((h L)/2)];
a22=2 phi2[-h L]-phi2[-((h L)/2)];

a30=1/2 phi1[-((h L)/2)]-1/84 phi2[-h L]-43/84 phi2[-((h L)/2)]+1/21 phi3[-h L]+1/42 phi3[-((h L)/2)];
a31=1/84 phi2[-h L]+10/21 phi2[-((h L)/2)]-1/42 phi3[-h L]-1/84 phi3[-((h L)/2)];
a32=1/84 phi2[-h L]+1/21 phi2[-((h L)/2)]-1/14 phi3[-h L]-1/28 phi3[-((h L)/2)];
a33=-(1/84) phi2[-h L]-1/84 phi2[-((h L)/2)]+1/21 phi3[-h L]+1/42 phi3[-((h L)/2)];

a40=phi1[-h L]-3 phi2[-h L]+4 phi3[-h L];
a41=-2 phi1[-((h L)/2)]-2 phi2[-((h L)/2)]-2 phi3[-((h L)/2)];
a42=-6 phi1[-((h L)/2)]-6 phi2[-((h L)/2)]-6 phi3[-((h L)/2)];
a43=-phi2[-h L]+4 phi3[-h L];
a44=8 phi1[-((h L)/2)]+4 phi2[-h L]+8 phi2[-((h L)/2)]-8 phi3[-h L]+8 phi3[-((h L)/2)];

a50=phi1[-h L]-2 phi2[-h L]+83/42 phi3[-h L];
a51=1/42 *phi3[-hL];
a52=2*phi2[-hL]-167/42 *phi3[-hL];
a53=83/42*phi3[-hL];
a54=0;
a55=0;

(*
a50=5/4 phi1[-h L]+11/12 phi1[-((h L)/2)]-9/4 phi2[-h L]+1/4 phi2[-((h L)/2)]-1929/560 phi3[-h L]-1/2 phi3[-((h L)/2)];
a51=1/2 phi1[-((h L)/2)]-2/3 phi2[-h L]-4/5 phi2[-((h L)/2)]+202/105 phi3[-h L]-1/2 phi3[-((h L)/2)];
a52=-(1/2) phi1[-h L]-25/12 phi1[-((h L)/2)]+19/6 phi2[-h L]+11/10 phi2[-((h L)/2)]+3961/840 phi3[-h L]+13/6 phi3[-((h L)/2)];
a53=-(1/4) phi1[-h L]+1/4 phi1[-((h L)/2)]-1/4 phi2[-h L]+1/2 phi2[-((h L)/2)]-1649/560 phi3[-h L]-1/2 phi3[-((h L)/2)];
a54=-(1/4) phi1[-((h L)/2)]-4/5 phi2[-((h L)/2)]+1/4 phi3[-h L]-2/3 phi3[-((h L)/2)];
a55=1/2 phi1[-h L]+2/3 phi1[-((h L)/2)]-1/4 phi2[-((h L)/2)]-1/2 phi3[-h L];
*)
(*  Checking whether conditions hold exactly for orders 1,2  *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]


(*Checking whether condition 3 holds weakly i.e.b_i(0) instead of b_i(-hL)*)(*Always reread Theorem 4.7 in the 05 paper when applying this step,"holding weakly" has a very particular meaning*)a50=a50/. L->0;
a51=a51/. L->0;
a52=a52/. L->0;
a53=a53/. L->0;
a54=a54/. L->0;
a55=a55/. L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);
phi5[0]=1/(5!);

(*The following line takes care of psi_p(0)=0 in "holding weakly"*)
order3a=order3a/. {L->0};
allorder3=Join[order3a,order3b];

FullSimplify[allorder3]

(* Checking that this method is never 4th order *)
allorder4=Join[order4a,order4b,order4c,order4d];
allorder4/.{L->0}


(* ::Input::Initialization:: *)
(* The fourth order method of ERK43ZBa *)
c0=0;c1=1/2;c2=1/2;c3=1;c4=1/2;c5=1;

a00=1/2 phi1[-((h L)/2)];

a10=1/2 phi1[-((h L)/2)]-2/3 phi2[-((h L)/2)];
a11=2/3 phi2[-((h L)/2)];

a20=phi1[-h L]-2 phi2[-h L];
a21=phi2[-((h L)/2)];
a22=2 phi2[-h L]-phi2[-((h L)/2)];

a30=1/2 phi1[-((h L)/2)]-1/84 phi2[-h L]-43/84 phi2[-((h L)/2)]+1/21 phi3[-h L]+1/42 phi3[-((h L)/2)];
a31=1/84 phi2[-h L]+10/21 phi2[-((h L)/2)]-1/42 phi3[-h L]-1/84 phi3[-((h L)/2)];
a32=1/84 phi2[-h L]+1/21 phi2[-((h L)/2)]-1/14 phi3[-h L]-1/28 phi3[-((h L)/2)];
a33=-(1/84) phi2[-h L]-1/84 phi2[-((h L)/2)]+1/21 phi3[-h L]+1/42 phi3[-((h L)/2)];a40=phi1[-h L]-3 phi2[-h L]+4 phi3[-h L];
a41=-2 phi1[-((h L)/2)]-2 phi2[-((h L)/2)]-2 phi3[-((h L)/2)];
a42=-6 phi1[-((h L)/2)]-6 phi2[-((h L)/2)]-6 phi3[-((h L)/2)];
a43=-phi2[-h L]+4 phi3[-h L];
a44=8 phi1[-((h L)/2)]+4 phi2[-h L]+8 phi2[-((h L)/2)]-8 phi3[-h L]+8 phi3[-((h L)/2)];

(*  Checking whether conditions hold exactly for orders 1,2,3  *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]
allorder3 = Join[order3a,order3b];
FullSimplify[ExpandNonCommutativeMultiply[allorder3]]

(*   Checking whether condition 4 holds weakly i.e. b_i(0) instead       of b_i(-hL)  *)
(*  Always reread Theorem 4.7 in the 05 paper when applying this step, "holding weakly" has a very particular meaning *) 
a40=a40/.L->0;
a41=a41/.L->0;
a42=a42/.L->0;
a43=a43/.L->0;
a44=a44/.L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);phi5[0]=1/(5!);
(* The following line takes care of psi_p(0)=0 in "holding weakly" *)
order4a = order4a/.{L->0};
allorder4=Join[order4a,order4b,order4c,order4d];
FullSimplify[ExpandNonCommutativeMultiply[allorder4]]


(* ::Input::Initialization:: *)
(* The fourth order method of ERK54ZBa *)
c0=0;c1=1/2;c2=1/2;c3=1/4;c4=1/2;
c5=1/5;
c6=2/3;
c7=1;
c8=1;

a00=1/2 phi1[-((h L)/2)];

a10=1/2 phi1[-((h L)/2)]-1/2 phi2[-((h L)/2)];
a11=1/2 phi2[-((h L)/2)];

a20=1/4 phi1[-((h L)/4)]-1/8 phi2[-((h L)/4)];
a21=0;
a22=1/8 phi2[-((h L)/4)];

a30=1/2 phi1[-((h L)/2)]-3/2 phi2[-((h L)/2)]+8/15 phi2[-((h L)/4)]+2 phi3[-((h L)/2)]-8/15 phi3[-((h L)/4)];
a31=0;
a32=-(1/2) phi2[-((h L)/2)]+8/15 phi2[-((h L)/4)]+2 phi3[-((h L)/2)]-8/15 phi3[-((h L)/4)];
a33=2 phi2[-((h L)/2)]-16/15 phi2[-((h L)/4)]-4 phi3[-((h L)/2)]+16/15 phi3[-((h L)/4)];

a40=1/5 phi1[-((h L)/5)]-6/25 phi2[-((h L)/5)]+16/125 phi3[-((h L)/5)];
a41=0;
a42=0;
a43=8/25 phi2[-((h L)/5)]-32/125 phi3[-((h L)/5)];
a44=-(2/25) phi2[-((h L)/5)]+16/125 phi3[-((h L)/5)];

a50=2/3 phi1[-((2 h L)/3)]-28/9 phi2[-((2 h L)/3)]-41/162 phi2[-((h L)/4)]+160/27 phi3[-((2 h L)/3)];
a51=0;
a52=0;
a53=-(164/81) phi2[-((h L)/4)];
a54=-(16/27) phi2[-((2 h L)/3)]+41/243 phi2[-((h L)/4)]+320/81 phi3[-((2 h L)/3)];a55=100/27 phi2[-((2 h L)/3)]+1025/486 phi2[-((h L)/4)]-800/81 phi3[-((2 h L)/3)];

a60=phi1[-h L]-69/10 phi2[-h L]+36/175 phi2[-((2 h L)/3)]+41/35 phi2[-((h L)/4)]-9/70 phi2[-((h L)/5)]+101/5 phi3[-h L]+48/175 phi3[-((2 h L)/3)]+3/25 phi3[-((h L)/5)]-18 phi4[-h L]-576/35 phi4[-((2 h L)/3)]-18/175 phi4[-((h L)/5)];
a61=0;
a62=0;
a63=512/25 phi2[-h L]+36864/875 phi2[-((2 h L)/3)]+328/35 phi2[-((h L)/4)]+312/175 phi2[-((h L)/5)]-6656/25 phi3[-h L]-319488/875 phi3[-((2 h L)/3)]-5856/875 phi3[-((h L)/5)]+4608/5 phi4[-h L]+147456/175 phi4[-((2 h L)/3)]+4608/875 phi4[-((h L)/5)];
a64=-(48/5) phi2[-h L]-2976/175 phi2[-((2 h L)/3)]-82/105 phi2[-((h L)/4)]-38/35 phi2[-((h L)/5)]+624/5 phi3[-h L]+26752/175 phi3[-((2 h L)/3)]+544/175 phi3[-((h L)/5)]-432 phi4[-h L]-13824/35 phi4[-((2 h L)/3)]-432/175 phi4[-((h L)/5)];
a65=-(50/7) phi2[-h L]-1560/49 phi2[-((2 h L)/3)]-205/21 phi2[-((h L)/4)]-50/49 phi2[-((h L)/5)]+1150/7 phi3[-h L]+13280/49 phi3[-((2 h L)/3)]+230/49 phi3[-((h L)/5)]-4500/7 phi4[-h L]-28800/49 phi4[-((2 h L)/3)]-180/49 phi4[-((h L)/5)];
a66=1107/350 phi2[-h L]+(39852 phi2[-((2 h L)/3)])/6125+(1107 phi2[-((h L)/5)])/2450-7533/175 phi3[-h L]-(361584 phi3[-((2 h L)/3)])/6125-(7533 phi3[-((h L)/5)])/6125+5994/35 phi4[-h L]+(191808 phi4[-((2 h L)/3)])/1225+(5994 phi4[-((h L)/5)])/6125;


(*  Checking whether conditions hold exactly for orders 1,2,3  *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]
allorder3 = Join[order3a,order3b];
FullSimplify[ExpandNonCommutativeMultiply[allorder3]]

(*   Checking whether condition 4 holds weakly i.e. b_i(0) instead       of b_i(-hL)  *)
(*  Always reread Theorem 4.7 in the 05 paper when applying this step, "holding weakly" has a very particular meaning *) 
a60=a60/.L->0;
a61=a61/.L->0;
a62=a62/.L->0;
a63=a63/.L->0;
a64=a64/.L->0;
a65=a65/.L->0;
a66=a66/.L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);phi5[0]=1/(5!);
(* The following line takes care of psi_p(0)=0 in "holding weakly" *)
order4a = order4a/.{L->0};
allorder4=Join[order4a,order4b,order4c,order4d];
FullSimplify[ExpandNonCommutativeMultiply[allorder4]]

(*Checking if this method is ever 5th order*)
allorder5 = Join[order5a,ExpandNonCommutativeMultiply[order5b],order5c,order5d,order5e,order5f,order5g,order5h,order5i];

allorder5/.{L->0}


(* ::Input::Initialization:: *)
(* Checking 5th order method*)
c0=0;c1=1/2;c2=1/2;c3=1/4;c4=1/2;
c5=1/5;
c6=2/3;
c7=1;
c8=1;


a00=1/2 *phi1[-((h L)/2)];

a10=1/2* phi1[-((h L)/2)]-1/2 *phi2[-((h L)/2)];
a11=1/2 *phi2[-((h L)/2)];

a20=1/4* phi1[-((h L)/4)]-1/8* phi2[-((h L)/4)];
a21=0;
a22=1/8 *phi2[-((h L)/4)];

a30=1/2* phi1[-((h L)/2)]-3/2* phi2[-((h L)/2)]+2* phi3[-((h L)/2)];
a31=0;
a32=-(1/2) *phi2[-((h L)/2)]+2 *phi3[-((h L)/2)];
a33=2 *phi2[-((h L)/2)]-4* phi3[-((h L)/2)];

a40=1/5 *phi1[-((h L)/5)]-6/25* phi2[-((h L)/5)]+16/125* phi3[-((h L)/5)];
a41=0;
a42=0;
a43=8/25* phi2[-((h L)/5)]-32/125* phi3[-((h L)/5)];
a44=-(2/25)* phi2[-((h L)/5)]+16/125* phi3[-((h L)/5)];

a50=2/3* phi1[-((2 h L)/3)]-28/9* phi2[-((2 h L)/3)]-5/162 *phi2[-((h L)/5)]+160/27 *phi3[-((2 h L)/3)]+2/81* phi3[-((h L)/5)];
a51=0;
a52=0;
a53=-(20/81)* phi2[-((h L)/5)]+16/81* phi3[-((h L)/5)];
a54=-(16/27) *phi2[-((2 h L)/3)]+5/243* phi2[-((h L)/5)]+320/81 *phi3[-((2 h L)/3)]-4/243* phi3[-((h L)/5)];
a55=100/27* phi2[-((2 h L)/3)]+125/486* phi2[-((h L)/5)]-800/81* phi3[-((2 h L)/3)]-50/243* phi3[-((h L)/5)];

a60=phi1[-h L]-17/2* phi2[-h L]-108/35 *phi2[-((2 h L)/3)]-3/14 *phi2[-((h L)/5)]+41* phi3[-h L]+144/5* phi3[-((2 h L)/3)]+3/5* phi3[-((h L)/5)]-90* phi4[-h L]-576/7 *phi4[-((2 h L)/3)]-18/35* phi4[-((h L)/5)];
a61=0;
a62=0;
a63=0;
a64=-(16/3)* phi2[-h L]-288/35* phi2[-((2 h L)/3)]-4/7* phi2[-((h L)/5)]+208/3* phi3[-h L]+384/5 *phi3[-((2 h L)/3)]+8/5 *phi3[-((h L)/5)]-240* phi4[-h L]-1536/7* phi4[-((2 h L)/3)]-48/35 *phi4[-((h L)/5)];
a65=250/21 *phi2[-h L]+360/49 *phi2[-((2 h L)/3)]+25/49* phi2[-((h L)/5)]-250/3* phi3[-h L]-480/7 *phi3[-((2 h L)/3)]-10/7* phi3[-((h L)/5)]+1500/7 *phi4[-h L]+9600/49 *phi4[-((2 h L)/3)]+60/49* phi4[-((h L)/5)];
a66=27/14* phi2[-h L]+972/245* phi2[-((2 h L)/3)]+27/98 *phi2[-((h L)/5)]-27* phi3[-h L]-1296/35* phi3[-((2 h L)/3)]-27/35 *phi3[-((h L)/5)]+810/7 *phi4[-h L]+5184/49 *phi4[-((2 h L)/3)]+162/245* phi4[-((h L)/5)];

a70=phi1[-h L]-15/2 *phi2[-h L]+28 *phi3[-h L]-45* phi4[-h L];
a71=0;a72=0;a73=0;
a74=0;
a75=125/14* phi2[-h L]-625/14* phi3[-h L]+1125/14* phi4[-h L];
a76=-(27/14)*phi2[-h L]+162/7 *phi3[-h L]-405/7 *phi4[-h L];
a77=1/2 *phi2[-h L]-13/2 *phi3[-h L]+
45/2 *phi4[-h L];

(*  Checking whether conditions hold exactly for orders 1,2,3  *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]
allorder3 = Join[order3a,order3b];
FullSimplify[ExpandNonCommutativeMultiply[allorder3]]
allorder4 = Join[order4a,order4b,order4c,order4d];
FullSimplify[ExpandNonCommutativeMultiply[allorder4]]

(*   Checking whether condition for order 5 holds weakly i.e. b_i(0) instead       of b_i(-hL)  *)
(*  Always reread Theorem 4.7 in the 05 paper when applying this step, "holding weakly" has a very particular meaning *) 

a70=a70/.L->0;
a71=a71/.L->0;
a72=a72/.L->0;
a73=a73/.L->0;
a74=a74/.L->0;
a75=a75/.L->0;
a76=a76/.L->0;
a77=a77/.L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);phi5[0]=1/(5!);

(* The following line takes care of psi_p(0)=0 in "holding weakly" *)
order5a = order5a/.{L->0};
allorder5 = Join[order5a,ExpandNonCommutativeMultiply[order5b],order5c,order5d,order5e,order5f,order5g,order5h,order5i];

FullSimplify[allorder5]


(* ::Input::Initialization:: *)
(* Checking 5th order method from ERK54ZBa*)
c0=0;c1=1/2;c2=1/2;c3=1/4;c4=1/2;
c5=1/5;
c6=2/3;
c7=1;
c8=1;

a00=1/2 phi1[-((h L)/2)];

a10=1/2 phi1[-((h L)/2)]-1/2 phi2[-((h L)/2)];
a11=1/2 phi2[-((h L)/2)];

a20=1/4 phi1[-((h L)/4)]-1/8 phi2[-((h L)/4)];
a21=0;
a22=1/8 phi2[-((h L)/4)];

a30=1/2 phi1[-((h L)/2)]-3/2 phi2[-((h L)/2)]+8/15 phi2[-((h L)/4)]+2 phi3[-((h L)/2)]-8/15 phi3[-((h L)/4)];
a31=0;
a32=-(1/2) phi2[-((h L)/2)]+8/15 phi2[-((h L)/4)]+2 phi3[-((h L)/2)]-8/15 phi3[-((h L)/4)];
a33=2 phi2[-((h L)/2)]-16/15 phi2[-((h L)/4)]-4 phi3[-((h L)/2)]+16/15 phi3[-((h L)/4)];

a40=1/5 phi1[-((h L)/5)]-6/25 phi2[-((h L)/5)]+16/125 phi3[-((h L)/5)];
a41=0;
a42=0;
a43=8/25 phi2[-((h L)/5)]-32/125 phi3[-((h L)/5)];
a44=-(2/25) phi2[-((h L)/5)]+16/125 phi3[-((h L)/5)];

a50=2/3 phi1[-((2 h L)/3)]-28/9 phi2[-((2 h L)/3)]-41/162 phi2[-((h L)/4)]+160/27 phi3[-((2 h L)/3)];
a51=0;
a52=0;
a53=-(164/81) phi2[-((h L)/4)];
a54=-(16/27) phi2[-((2 h L)/3)]+41/243 phi2[-((h L)/4)]+320/81 phi3[-((2 h L)/3)];a55=100/27 phi2[-((2 h L)/3)]+1025/486 phi2[-((h L)/4)]-800/81 phi3[-((2 h L)/3)];

a60=phi1[-h L]-69/10 phi2[-h L]+36/175 phi2[-((2 h L)/3)]+41/35 phi2[-((h L)/4)]-9/70 phi2[-((h L)/5)]+101/5 phi3[-h L]+48/175 phi3[-((2 h L)/3)]+3/25 phi3[-((h L)/5)]-18 phi4[-h L]-576/35 phi4[-((2 h L)/3)]-18/175 phi4[-((h L)/5)];
a61=0;
a62=0;
a63=512/25 phi2[-h L]+36864/875 phi2[-((2 h L)/3)]+328/35 phi2[-((h L)/4)]+312/175 phi2[-((h L)/5)]-6656/25 phi3[-h L]-319488/875 phi3[-((2 h L)/3)]-5856/875 phi3[-((h L)/5)]+4608/5 phi4[-h L]+147456/175 phi4[-((2 h L)/3)]+4608/875 phi4[-((h L)/5)];
a64=-(48/5) phi2[-h L]-2976/175 phi2[-((2 h L)/3)]-82/105 phi2[-((h L)/4)]-38/35 phi2[-((h L)/5)]+624/5 phi3[-h L]+26752/175 phi3[-((2 h L)/3)]+544/175 phi3[-((h L)/5)]-432 phi4[-h L]-13824/35 phi4[-((2 h L)/3)]-432/175 phi4[-((h L)/5)];
a65=-(50/7) phi2[-h L]-1560/49 phi2[-((2 h L)/3)]-205/21 phi2[-((h L)/4)]-50/49 phi2[-((h L)/5)]+1150/7 phi3[-h L]+13280/49 phi3[-((2 h L)/3)]+230/49 phi3[-((h L)/5)]-4500/7 phi4[-h L]-28800/49 phi4[-((2 h L)/3)]-180/49 phi4[-((h L)/5)];
a66=1107/350 phi2[-h L]+(39852 phi2[-((2 h L)/3)])/6125+(1107 phi2[-((h L)/5)])/2450-7533/175 phi3[-h L]-(361584 phi3[-((2 h L)/3)])/6125-(7533 phi3[-((h L)/5)])/6125+5994/35 phi4[-h L]+(191808 phi4[-((2 h L)/3)])/1225+(5994 phi4[-((h L)/5)])/6125;

a70=phi1[-h L]-15/2 phi2[-h L]+28 phi3[-h L]-45 phi4[-h L];
a71=0;
a72=0;
a73=0;
a74=0;
a75=125/14 phi2[-h L]-625/14 phi3[-h L]+1125/14 phi4[-h L];
a76=-(27/14) phi2[-h L]+162/7 phi3[-h L]-405/7 phi4[-h L];
a77=1/2 phi2[-h L]-13/2 phi3[-h L]+45/2 phi4[-h L];
(*  Checking whether conditions hold exactly for orders 1,2,3  *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]
allorder3 = Join[order3a,order3b];
FullSimplify[ExpandNonCommutativeMultiply[allorder3]]
allorder4 = Join[order4a,order4b,order4c,order4d];
FullSimplify[ExpandNonCommutativeMultiply[allorder4]]

(*   Checking whether condition for order 5 holds weakly i.e. b_i(0) instead       of b_i(-hL)  *)
(*  Always reread Theorem 4.7 in the 05 paper when applying this step, "holding weakly" has a very particular meaning *) 

a70=a70/.L->0;
a71=a71/.L->0;
a72=a72/.L->0;
a73=a73/.L->0;
a74=a74/.L->0;
a75=a75/.L->0;
a76=a76/.L->0;
a77=a77/.L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);phi5[0]=1/(5!);

(* The following line takes care of psi_p(0)=0 in "holding weakly" *)
order5a = order5a/.{L->0};
allorder5 = Join[order5a,ExpandNonCommutativeMultiply[order5b],ExpandNonCommutativeMultiply[order5c],order5d,order5e,order5f,order5g,order5h,order5i];

FullSimplify[allorder5]
